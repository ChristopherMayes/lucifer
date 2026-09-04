!+
! Module fel_device_mod
!
! The Fortran side of the device seam: the only unit that calls the C interface of
! lucifer_device.h, behind which one backend file sits per build (the Metal backend
! where the toolchain can carry it, a refusing stub in every other build, including a
! macOS one whose Objective-C++ compiler is a GNU one). Everything device-shaped
! crosses here: the chart conversions, the residency bookkeeping, the lockstep
! instrument's device role, and the refusals. No Metal type or call appears outside
! device/lucifer_metal.mm, and no luc_dev call appears outside this module. The FP64
! comparison quantities the instrument's rows are measured against (the roll-off, the
! source deposit, the median) are fel_fp32_mod's own module procedures, called from
! here rather than repeated: they are transcriptions of fel_field_step's conventions,
! and a second copy is a second thing to find the day a convention changes.
!
! The device chart. Transverse coordinates and the energy offset are FP32; the
! longitudinal state is a 64-bit fixed-point phase in ticks of 2 pi / 2^32, held off a
! static FP64 per-slice reference z_ref. The conversions are exact FP64 both ways:
!
!   up:    goff = gamma - gamma0,  uphase = nint((ks z / beta - ks z_ref) * ticks)
!   down:  gamma = gamma0 + goff, pz = (sqrt(gamma^2-1) - p0_mc)/p0_mc,
!          z = beta (z_ref + delta / ks) with delta = uphase / ticks
!
! (the down form is fel_advance's own exit chart z = -beta (phi0 - theta)/ks with
! theta = phi0 + ks z_ref + delta, so the boundary is the CPU chart evaluated in
! FP64). Bucket wraps in the accumulator are exact integer arithmetic -- the migration
! operation a later landing needs -- and fel_device_setup asserts that exactly on the
! device itself, not to a tolerance. The kernel-side reformulations (energy offset,
! detuning difference, base rotator times small angle) mirror fel_fp32_mod, whose
! lockstep levels price them; the fixed-point phase replaces the FP32 residual dzr,
! whose ulp floor forced that module's moving reference (FINDINGS 7.39): with a
! uniform 1.5e-9 rad quantum the reference can stay static for a whole element.
!
! Two roles. With fp32_check off the device is the run inside averaged FEL elements:
! beam and field upload at element entry, stay resident through the element, and come
! back at the comb's stats positions and the element end (fel_track_line_mod owns that
! schedule). With fp32_check = 'lockstep' or 'freerun' the device takes the CPU twin's
! role in the instrument: the FP64 path runs untouched, the device advances the same
! step from the shared state, and the rows, guard, stream and ceilings are
! fel_fp32_mod's own machinery with the device's arithmetic under test. The guard
! column is then the median per-step phase increment in ticks of the fixed-point
! quantum rather than in FP32 ulps; the same floor applies.
!
! The device deposit accumulates with atomic adds whose ordering is not fixed, so two
! runs of the same step differ in the source's last bit or two. Both reference
! backends behave the same way (manual/GPU.md, gpu/metal-engine 4919b01), so no device
! output is asserted byte-identical; the ceilings absorb it.
!
! The field set. The resident field is every member of the run's set (the fundamental
! and its harmonics, one grid) with one or two planes each (Ex, or the (Ex, Ey) pair
! when two polarizations are live). The routines here take the set as ff(:), the
! fel_field_struct of each member in set order, and the step's constants carry each member's
! harmonic, coupling and deposit scale plus the element's polarization pair
! (fel_device_par_struct). Every member rides the one fixed-point phase, at h times
! it, and every plane slips together, as on the CPU. One member and one plane is the
! single-field path it was before the set.
!-

module fel_device_mod

use fel_beam_mod
use fel_field_mod
use fel_fp32_mod
use wavefront_mod

use, intrinsic :: iso_c_binding
use, intrinsic :: ieee_exceptions, only: ieee_usual, ieee_get_flag, ieee_set_flag

implicit none

! Ticks of phase: 2^32 per radiation period, the fixed-point quantum. The FP64 forms
! here are exact; the kernels carry their own FP32 roundings of the same constants,
! which is part of the arithmetic the instrument prices.

real(rp), parameter :: fel_dev_ticks_per_rad$ = 4294967296.0_rp / twopi
real(rp), parameter :: fel_dev_rad_per_tick$ = twopi / 4294967296.0_rp

! The field set's size bound, LUC_DEV_MAX_FIELD in lucifer_device.h and the
! harmonics(9) list's length. The two mirrors move together.

integer, parameter :: fel_dev_max_field$ = 9

!+
! Struct fel_device_par_struct
!
! One step's constants for luc_dev_step, mirroring luc_dev_step_par of
! lucifer_device.h field for field (all doubles, then the 32-bit ints). Editing
! either mirror alone skews the layout silently; keep them together.
!-

type, bind(c) :: fel_device_par_struct
  real(c_double) :: dz, ks, ku, aw, qres, e0, gam0, p0_mc
  real(c_double) :: kx, ky, ax, ay, cos_t, sin_t
  real(c_double) :: k1x, k1y, gridmax, dgrid
  real(c_double) :: harm(fel_dev_max_field$)    ! Per member: the harmonic number as a factor.
  real(c_double) :: rtmp(fel_dev_max_field$)    ! Per member: fc(h) / (sqrt(2) m_e).
  real(c_double) :: scl_w(fel_dev_max_field$)   ! Per member: fel_field_step's scl_w at h.
  real(c_double) :: pol_re(2), pol_im(2)        ! The element's polarization pair.
  integer(c_int) :: first, helical, mutate, nfield, npol, pad
end type

!+
! Struct fel_device_struct
!
! The device run state: the knob's answer, the residency flag the walk consults, the
! static per-slice references, the kernel-upload key (fp32_kernel_cache's key), and
! the staging arrays the transfers reuse. su0 holds the previous step's phase
! accumulators for the instrument's guard statistic.
!-

type fel_device_struct
  logical :: on = .false.               ! A backend is armed for this run.
  ! resident means the device state IS the run: only the production role sets it, and
  ! the walk's readbacks key on it. The twin's carried freerun state gets its own
  ! flag, because a readback of twin state into the FP64 arrays would be the
  ! instrument steering the run it observes -- the exact failure the read-only proof
  ! exists to catch, and the first thing it caught.
  logical :: resident = .false.         ! Production role: beam and field live on the device.
  logical :: twin_live = .false.        ! Instrument role: the freerun twin state carries.
  integer :: nslice = 0, npart = 0, ngrid = 0
  integer :: nfield = 1, npol = 1       ! The set: members, and planes per member.
  integer, allocatable :: harm(:)       ! Each member's harmonic number.
  real(rp), allocatable :: z_ref(:)     ! Static FP64 longitudinal reference per slice [m].
  real(rp), allocatable :: k_key(:,:)   ! Propagator upload key, (4, nfield).
  character(64) :: name = ''            ! The device, for the log line.
  logical :: wrap_exact = .false.       ! The exact-wrap assertion's verdict.
  integer :: nstep = 0                  ! Device steps encoded this run.
  ! The instrument's per-member record: run-level worst phasor, source and field
  ! divergence of each member, the footer's attribution when the set has more than
  ! one member (the ten-column stream itself carries the worst over members).
  real(rp), allocatable :: worst_h(:,:) ! (3, nfield).
  ! Staging, sized once: per-slice particle arrays padded to npart with zero-weight
  ! entries (the device wants a rectangular array; padding radiates nothing and is
  ! never read back), one field plane, and the per-slice phase rotators per member.
  real(c_float), allocatable :: sx(:), spx(:), sy(:), spy(:), sg(:), sw(:)
  integer(c_int64_t), allocatable :: su(:)
  integer(c_int64_t), allocatable :: su0(:,:)
  complex(c_float_complex), allocatable :: se(:,:), sbase(:), sbdep(:)
end type

! The C interface, lucifer_device.h. One implementation per build sits behind it.

interface

  function luc_dev_available (name, name_len, reason, reason_len) &
                              bind(c, name = 'luc_dev_available') result (ok)
    import c_char, c_int
    character(kind=c_char) :: name(*), reason(*)
    integer(c_int), value :: name_len, reason_len
    integer(c_int) ok
  end function

  function luc_dev_init (nslice, npart, ngrid, nfield, npol, reason, reason_len) &
                         bind(c, name = 'luc_dev_init') result (ierr)
    import c_char, c_int
    integer(c_int), value :: nslice, npart, ngrid, nfield, npol, reason_len
    character(kind=c_char) :: reason(*)
    integer(c_int) ierr
  end function

  subroutine luc_dev_close () bind(c, name = 'luc_dev_close')
  end subroutine

  subroutine luc_dev_upload_slice (is, n, x, px, y, py, goff, uphase, w) &
                                   bind(c, name = 'luc_dev_upload_slice')
    import c_int, c_float, c_int64_t
    integer(c_int), value :: is, n
    real(c_float) :: x(*), px(*), y(*), py(*), goff(*), w(*)
    integer(c_int64_t) :: uphase(*)
  end subroutine

  subroutine luc_dev_download_slice (is, n, x, px, y, py, goff, uphase) &
                                     bind(c, name = 'luc_dev_download_slice')
    import c_int, c_float, c_int64_t
    integer(c_int), value :: is, n
    real(c_float) :: x(*), px(*), y(*), py(*), goff(*)
    integer(c_int64_t) :: uphase(*)
  end subroutine

  subroutine luc_dev_upload_field_slice (im, ip, is, e) bind(c, name = 'luc_dev_upload_field_slice')
    import c_int, c_float_complex
    integer(c_int), value :: im, ip, is
    complex(c_float_complex) :: e(*)
  end subroutine

  subroutine luc_dev_download_field_slice (im, ip, is, e) bind(c, name = 'luc_dev_download_field_slice')
    import c_int, c_float_complex
    integer(c_int), value :: im, ip, is
    complex(c_float_complex) :: e(*)
  end subroutine

  subroutine luc_dev_zero_field_slice (im, ip, is) bind(c, name = 'luc_dev_zero_field_slice')
    import c_int
    integer(c_int), value :: im, ip, is
  end subroutine

  subroutine luc_dev_download_source_slice (im, is, s) bind(c, name = 'luc_dev_download_source_slice')
    import c_int, c_float_complex
    integer(c_int), value :: im, is
    complex(c_float_complex) :: s(*)
  end subroutine

  subroutine luc_dev_set_kernel (im, expk) bind(c, name = 'luc_dev_set_kernel')
    import c_int, c_float_complex
    integer(c_int), value :: im
    complex(c_float_complex) :: expk(*)
  end subroutine

  subroutine luc_dev_set_slice_phases (base, base_dep) bind(c, name = 'luc_dev_set_slice_phases')
    import c_float_complex
    complex(c_float_complex) :: base(*), base_dep(*)
  end subroutine

  function luc_dev_step (par, cret_ticks, reason, reason_len) &
                         bind(c, name = 'luc_dev_step') result (ierr)
    import c_char, c_int, c_int64_t
    import fel_device_par_struct
    type (fel_device_par_struct) :: par
    integer(c_int64_t), value :: cret_ticks
    character(kind=c_char) :: reason(*)
    integer(c_int), value :: reason_len
    integer(c_int) ierr
  end function

  subroutine luc_dev_sync () bind(c, name = 'luc_dev_sync')
  end subroutine

  function luc_dev_wrap_check (bucket_ticks) bind(c, name = 'luc_dev_wrap_check') result (ierr)
    import c_int, c_int64_t
    integer(c_int64_t), value :: bucket_ticks
    integer(c_int) ierr
  end function

  function luc_dev_seconds () bind(c, name = 'luc_dev_seconds') result (s)
    import c_double
    real(c_double) s
  end function

  function luc_dev_bytes () bind(c, name = 'luc_dev_bytes') result (b)
    import c_int64_t
    integer(c_int64_t) b
  end function

end interface

! A generic name, kept inside the module: every caller of it is here, and the seam's
! strings are the only reason it exists.

private from_c

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_query (usable, device, reason)
!
! Routine to report what device backend this build carries, without arming anything.
! luc_dev_available is the one call a stub build answers, and it allocates nothing, so
! this is safe and free in every build. It is what the usage banner reports and what a
! caller reads to decide whether device = "metal" is worth asking for at all.
!
! Two things can make a backend unusable, and they are different: a build compiled
! against the stub carries none at all, and a build carrying the Metal backend can
! still run on a machine with no Metal device or no unified memory. The reason
! distinguishes them. Neither is a platform test: a caller that keys on the platform
! rather than on this answer will skip the device on a machine that has one.
!
! The floating-point exception flags are saved and restored across the call. Metal's
! own device enumeration raises overflow on this machine, and a program whose last act
! is this query would then report a signalling overflow at exit, attributing the
! framework's arithmetic to the run's own. A query answers a question and leaves no
! trace.
!
! Output:
!   usable -- logical: True if device = "metal" would be accepted here.
!   device -- character(*): The device's own name, or 'none' when there is no backend.
!   reason -- character(*): Why it is unusable. Blank when usable is True.
!-

subroutine fel_device_query (usable, device, reason)

logical usable
character(*) device, reason
character(kind=c_char) c_name(64), c_reason(256)
logical flag_was(size(ieee_usual))

!

call ieee_get_flag (ieee_usual, flag_was)
usable = (luc_dev_available(c_name, 64, c_reason, 256) /= 0)
call ieee_set_flag (ieee_usual, flag_was)

call from_c (c_name, device)
call from_c (c_reason, reason)
if (usable) reason = ''

end subroutine fel_device_query

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine from_c (c_str, f_str)
!
! Routine to copy a null-terminated C string into a Fortran one, truncating at the
! Fortran length. Every reason and name the seam returns arrives this way.
!
! Input:
!   c_str -- character(kind=c_char)(*): The C string.
!
! Output:
!   f_str -- character(*): The Fortran string, blank padded.
!-

subroutine from_c (c_str, f_str)

character(kind=c_char) c_str(*)
character(*) f_str
integer i

!

f_str = ''
do i = 1, len(f_str)
  if (c_str(i) == c_null_char) exit
  f_str(i:i) = c_str(i)
enddo

end subroutine from_c

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_setup (dev, device_req, beam, ngrid, sample, harm, npol, fp32_iu, err_flag)
!
! Routine to arm the device from the input knob. '' or 'off' leaves it dark. 'metal'
! asks for the one backend this tree knows; a build without it (the stub) refuses with
! the stub's own reason, and any other value is refused, with the error listing what
! this tree knows. Allocates the resident buffers for this run shape (the grid
! refusal, naming the nearest supported size, comes back from the backend), runs the
! exact-wrap assertion on the device, and writes the device header lines into the
! instrument's stream when one is open.
!
! Input:
!   device_req -- character(*): global%device.
!   beam       -- fel_beam_struct: The beam (slice fills size the rectangular array).
!   ngrid      -- integer: Transverse grid points per side, one grid for the set.
!   sample     -- real(rp): Slice spacing in radiation wavelengths (integer by check).
!   harm(:)    -- integer: The set's harmonic numbers, in set order. size(harm) is the member count.
!   npol       -- integer: Planes per member, 1 or 2.
!   fp32_iu    -- integer: The instrument's stream unit, 0 when dark.
!
! Output:
!   dev        -- fel_device_struct: Armed (or left dark).
!   err_flag   -- logical: Set True on any refusal. False otherwise.
!-

subroutine fel_device_setup (dev, device_req, beam, ngrid, sample, harm, npol, fp32_iu, err_flag)

type (fel_device_struct) dev
type (fel_beam_struct) beam
character(*) device_req
integer ngrid, fp32_iu, npol
integer harm(:)
real(rp) sample
logical err_flag

character(kind=c_char) c_reason(256)
character(256) reason
character(64) name
integer is, np, ierr
integer(c_int64_t) bucket
logical usable
character(*), parameter :: r_name = 'fel_device_setup'

!

err_flag = .false.

select case (device_req)
case ('', 'off')
  return
case ('metal')
case default
  call out_io (s_error$, r_name, 'UNRECOGNIZED DEVICE: "' // trim(device_req) // '".', &
                                 'THIS TREE KNOWS off AND metal.')
  err_flag = .true.
  return
end select

call fel_device_query (usable, name, reason)
if (.not. usable) then
  call out_io (s_error$, r_name, 'DEVICE = "metal" REFUSED: ' // trim(reason) // '.')
  err_flag = .true.
  return
endif
dev%name = name

! The rectangular particle array: padded to the largest fill, so uneven slice fills
! cost padding rather than a refusal. The padding has zero weight and is never read
! back.

np = 0
do is = 1, size(beam%slice)
  np = max(np, beam%slice(is)%n)
enddo
if (np == 0) then
  call out_io (s_error$, r_name, 'DEVICE = "metal" WITH NO PARTICLES.')
  err_flag = .true.
  return
endif

dev%nfield = size(harm)
dev%npol = npol
dev%harm = harm

ierr = luc_dev_init(size(beam%slice), np, ngrid, dev%nfield, npol, c_reason, 256)
if (ierr /= 0) then
  call from_c (c_reason, reason)
  call out_io (s_error$, r_name, 'DEVICE = "metal" REFUSED: ' // trim(reason) // '.')
  err_flag = .true.
  return
endif

dev%nslice = size(beam%slice)
dev%npart = np
dev%ngrid = ngrid

! The exact-wrap assertion, on the device's own arithmetic: a bucket shift and its
! return must be bit-exact, and the extracted phase must never see the shift. Wraps
! are modular arithmetic, so this asserts exactly rather than to a tolerance. The
! bucket is sample whole periods, which requires an integer sample.

if (abs(sample - nint(sample)) > 0) then
  call out_io (s_error$, r_name, 'DEVICE = "metal" NEEDS AN INTEGER WINDOW_SAMPLE FOR EXACT', &
                                 'BUCKET ARITHMETIC; THIS RUN HAS \es10.2\ .', r_array = [sample])
  err_flag = .true.
  return
endif
bucket = int(nint(sample), c_int64_t) * 4294967296_c_int64_t
if (luc_dev_wrap_check(bucket) /= 0) then
  call out_io (s_error$, r_name, 'THE DEVICE EXACT-WRAP ASSERTION FAILED: A BUCKET SHIFT DID NOT', &
                                 'RETURN BIT-EXACTLY. THIS IS A KERNEL BUG BY DEFINITION.')
  err_flag = .true.
  return
endif
dev%wrap_exact = .true.

allocate (dev%z_ref(dev%nslice))
dev%z_ref = 0
allocate (dev%k_key(4, dev%nfield))
dev%k_key = -1
allocate (dev%worst_h(3, dev%nfield))
dev%worst_h = 0
allocate (dev%sx(np), dev%spx(np), dev%sy(np), dev%spy(np), dev%sg(np), dev%sw(np))
allocate (dev%su(np), dev%su0(np, dev%nslice))
allocate (dev%se(ngrid, ngrid))
allocate (dev%sbase(dev%nfield * dev%nslice), dev%sbdep(dev%nfield * dev%nslice))
dev%su0 = 0

dev%on = .true.
call out_io (s_info$, r_name, 'Device: ' // trim(dev%name) // ', \i0\ MB resident, ' // &
             '\i0\ field member(s) of \i0\ plane(s).', &
             i_array = [int(luc_dev_bytes() / 2**20), dev%nfield, npol])

if (fp32_iu /= 0) then
  write (fp32_iu, '(a)') '# device = ' // trim(dev%name) // ' (the device holds the twin''s role;'
  write (fp32_iu, '(a)') '#   guard_ulp is the median per-step phase increment in ticks of 2 pi / 2^32)'
  write (fp32_iu, '(a, l1)') '# device_wrap_exact ', dev%wrap_exact
endif

end subroutine fel_device_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_element_begin (dev, beam, wf)
!
! Routine to make the beam and the field set resident at an FEL element's entry: the
! static per-slice reference is set to the slice's FP64 mean z, every slice converts
! through the module header's chart, and every plane of every member rounds to FP32
! in ring order. From here to fel_device_element_end the device state is the run
! inside the element; the host arrays are stale until a readback refreshes them.
!
! Input:
!   beam  -- fel_beam_struct: The beam.
!   ff(:) -- fel_field_struct: The field set, in set order.
!
! Output:
!   dev   -- fel_device_struct: Resident.
!-

subroutine fel_device_element_begin (dev, beam, ff)

type (fel_device_struct) dev
type (fel_beam_struct) beam
type (fel_field_struct) ff(:)
integer is

!

do is = 1, size(beam%slice)
  dev%z_ref(is) = 0
  if (beam%slice(is)%n > 0) dev%z_ref(is) = &
                            sum(beam%slice(is)%z(1:beam%slice(is)%n)) / beam%slice(is)%n
  call fel_device_stage_slice (dev, beam, is, dev%z_ref(is))
  call luc_dev_upload_slice (is-1, dev%npart, dev%sx, dev%spx, dev%sy, dev%spy, &
                             dev%sg, dev%su, dev%sw)
enddo
call fel_device_upload_fields (dev, ff)
dev%resident = .true.

end subroutine fel_device_element_begin

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_upload_fields (dev, ff)
!
! Routine to round every plane of every member of the set to FP32 and upload it in
! ring order: plane 1 is Ex and plane 2 is Ey, present only when the run carries two
! polarizations, in which case every member carries both.
!-

subroutine fel_device_upload_fields (dev, ff)

type (fel_device_struct) dev
type (fel_field_struct) ff(:)
integer im, ip, is

!

do im = 1, dev%nfield
  do ip = 1, dev%npol
    do is = 1, size(ff(im)%wf%Ex, 3)
      if (ip == 1) then
        dev%se = cmplx(ff(im)%wf%Ex(:,:,is), kind = c_float_complex)
      else
        dev%se = cmplx(ff(im)%wf%Ey(:,:,is), kind = c_float_complex)
      endif
      call luc_dev_upload_field_slice (im-1, ip-1, is-1, dev%se)
    enddo
  enddo
enddo

end subroutine fel_device_upload_fields

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_set_kernel (dev, im, exp_k2)
!
! Routine to upload member im's propagator exp(K2 dz), rounded from the FP64 kernel
! and keyed on the values themselves (fp32_kernel_cache's own key, per member), so
! an unchanged kernel costs no transfer. Each member diffracts at its own wavelength.
!
! Input:
!   im     -- integer: The member, 1-based in set order.
!   exp_k2 -- complex(rp): The member's FP64 propagator this step.
!-

subroutine fel_device_set_kernel (dev, im, exp_k2)

type (fel_device_struct) dev
integer im
complex(rp) exp_k2(:,:)
real(rp) key(4)
integer ng

!

ng = size(exp_k2, 1)
key = [real(ng, rp), real(exp_k2(1,1), rp), aimag(exp_k2(1,1)), aimag(exp_k2(ng/2+1, ng/2+1))]
if (any(dev%k_key(:, im) /= key)) then
  dev%se = cmplx(exp_k2, kind = c_float_complex)
  call luc_dev_set_kernel (im-1, dev%se)
  dev%k_key(:, im) = key
endif

end subroutine fel_device_set_kernel

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_stage_slice (dev, beam, is, z_ref)
!
! Routine to fill the staging arrays with one slice's device image: the exact FP64
! chart conversions of the module header, rounded once into the device
! representation, padded to the rectangular width with zero-weight entries.
!-

subroutine fel_device_stage_slice (dev, beam, is, z_ref)

type (fel_device_struct) dev
type (fel_beam_struct) beam
integer is, ip, n
real(rp) z_ref, p0_mc, gamma0, ks, p_mc, gam, beta, delta

!

p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)
ks = twopi / beam%wavelength
n = beam%slice(is)%n

do ip = 1, n
  dev%sx(ip) = real(beam%slice(is)%x(ip), c_float)
  dev%spx(ip) = real(beam%slice(is)%px(ip), c_float)
  dev%sy(ip) = real(beam%slice(is)%y(ip), c_float)
  dev%spy(ip) = real(beam%slice(is)%py(ip), c_float)
  dev%sw(ip) = real(beam%slice(is)%weight(ip), c_float)
  p_mc = p0_mc * (1 + beam%slice(is)%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam
  dev%sg(ip) = real(gam - gamma0, c_float)
  delta = ks * beam%slice(is)%z(ip) / beta - ks * z_ref
  dev%su(ip) = nint(delta * fel_dev_ticks_per_rad$, c_int64_t)
enddo
do ip = n+1, dev%npart
  dev%sx(ip) = 0;  dev%spx(ip) = 0;  dev%sy(ip) = 0;  dev%spy(ip) = 0
  dev%sg(ip) = 0;  dev%sw(ip) = 0;  dev%su(ip) = 0
enddo
dev%su0(:, is) = dev%su

end subroutine fel_device_stage_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_step_run (dev, par, phi0, phi0_new, ks, err_flag)
!
! Routine to encode one resident integration step: the per-slice phase rotators of
! every member (the push reads e^{-i h (phi0 + ks z_ref)}, the deposit its phi0_new
! counterpart), the exact tick count of the common phase advance, and one command
! buffer holding the whole step. The propagators are the caller's, uploaded through
! fel_device_set_kernel before this. Nothing is waited on here; the next host touch
! of a buffer drains it.
!
! Input:
!   par       -- fel_device_par_struct: The step's constants, filled by the caller.
!   phi0      -- real(rp): The common phase at step entry.
!   phi0_new  -- real(rp): The common phase at step exit.
!   ks        -- real(rp): The fundamental radiation wavenumber [1/m].
!
! Output:
!   dev       -- fel_device_struct: One more step encoded.
!   err_flag  -- logical: Set True if the backend refuses. False otherwise.
!-

subroutine fel_device_step_run (dev, par, phi0, phi0_new, ks, err_flag)

type (fel_device_struct) dev
type (fel_device_par_struct) par
real(rp) phi0, phi0_new, ks
logical err_flag

real(rp) phi_ref, h
integer is, im, ierr, k
integer(c_int64_t) cret
character(kind=c_char) c_reason(256)
character(256) reason
character(*), parameter :: r_name = 'fel_device_step_run'

!

err_flag = .false.

! The per-slice rotators per member, FP64 once per slice per step: the push works
! against the entry phase, the deposit against the exit phase, exactly the two epochs
! fel_fp32_twin_slice and fel_fp32_field_twin use, each at h times the fundamental's
! phase (fel_field_step's harm*theta). Member-major, the seam's order.

do im = 1, dev%nfield
  h = par%harm(im)
  do is = 1, dev%nslice
    k = (im - 1) * dev%nslice + is
    phi_ref = h * (phi0 + ks * dev%z_ref(is))
    dev%sbase(k) = cmplx(cmplx(cos(phi_ref), -sin(phi_ref), rp), kind = c_float_complex)
    phi_ref = h * (phi0_new + ks * dev%z_ref(is))
    dev%sbdep(k) = cmplx(cmplx(cos(phi_ref), -sin(phi_ref), rp), kind = c_float_complex)
  enddo
enddo
call luc_dev_set_slice_phases (dev%sbase, dev%sbdep)

cret = nint((phi0_new - phi0) * fel_dev_ticks_per_rad$, c_int64_t)
ierr = luc_dev_step (par, cret, c_reason, 256)
if (ierr /= 0) then
  call from_c (c_reason, reason)
  call out_io (s_error$, r_name, 'DEVICE STEP REFUSED: ' // trim(reason) // '.')
  err_flag = .true.
  return
endif
dev%nstep = dev%nstep + 1

end subroutine fel_device_step_run

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_readback (dev, beam, wf)
!
! Routine to refresh the host FP64 state from the resident device state: the beam
! through the module header's exact down-chart, the field record slice for slice in
! ring order. The comb's stats positions and the element end call this; between them
! the host arrays are stale by design ("stats read back at the comb's positions
! only"). A readback observes and never steers: the device state is untouched.
!
! Parallel over slices, and the measurement that made it so is in
! doc/performance.md's device section: at the finest comb the serial conversion loops
! were half the cost of every stats row, ~4.3 ms a record on the 96 x 8192 case. The
! conversions are elementwise with no accumulation and each iteration writes only its
! own slice's arrays or its own field plane, so any thread order computes the same
! bits. Two structural points hold that. The device drains once, serially, before the
! parallel regions, after which every transfer is a pure disjoint copy from the
! shared-storage buffers (lucifer_device.h's post-drain contract). And the staging is
! block-local per slice rather than the module-level scratch, which is sized for one
! slice at a time and is not thread-safe by design.
!-

subroutine fel_device_readback (dev, beam, ff)

type (fel_device_struct) dev
type (fel_beam_struct) beam
type (fel_field_struct) ff(:)
real(rp) p0_mc, gamma0, ks
integer is, im, ip

!

if (.not. dev%resident) return

p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)
ks = twopi / beam%wavelength

! The one serial drain: everything encoded completes here, and the downloads inside
! the parallel loops below never touch the backend's encoder state.

call luc_dev_sync ()

!$OMP parallel do
do is = 1, size(beam%slice)
  block
    real(c_float), allocatable :: bx(:), bpx(:), by(:), bpy(:), bg(:)
    integer(c_int64_t), allocatable :: bu(:)
    real(rp) gam, p_mc, beta, delta
    integer ip, n
    allocate (bx(dev%npart), bpx(dev%npart), by(dev%npart), bpy(dev%npart), bg(dev%npart))
    allocate (bu(dev%npart))
    n = beam%slice(is)%n
    call luc_dev_download_slice (is-1, dev%npart, bx, bpx, by, bpy, bg, bu)
    do ip = 1, n
      beam%slice(is)%x(ip) = real(bx(ip), rp)
      beam%slice(is)%px(ip) = real(bpx(ip), rp)
      beam%slice(is)%y(ip) = real(by(ip), rp)
      beam%slice(is)%py(ip) = real(bpy(ip), rp)
      gam = gamma0 + real(bg(ip), rp)
      p_mc = sqrt(gam**2 - 1)
      beam%slice(is)%pz(ip) = (p_mc - p0_mc) / p0_mc
      beta = p_mc / gam
      delta = real(bu(ip), rp) * fel_dev_rad_per_tick$
      beam%slice(is)%z(ip) = beta * (dev%z_ref(is) + delta / ks)
    enddo
  end block
enddo
!$OMP end parallel do

do im = 1, dev%nfield
  do ip = 1, dev%npol
    !$OMP parallel do
    do is = 1, size(ff(im)%wf%Ex, 3)
      block
        complex(c_float_complex), allocatable :: be(:,:)
        allocate (be(size(ff(im)%wf%Ex, 1), size(ff(im)%wf%Ex, 2)))
        call luc_dev_download_field_slice (im-1, ip-1, is-1, be)
        if (ip == 1) then
          ff(im)%wf%Ex(:,:,is) = cmplx(be, kind = wf_rp)
        else
          ff(im)%wf%Ey(:,:,is) = cmplx(be, kind = wf_rp)
        endif
      end block
    enddo
    !$OMP end parallel do
  enddo
enddo

end subroutine fel_device_readback

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_element_end (dev, beam, wf)
!
! Routine to end residency at the element's exit: the final readback, after which the
! host FP64 state is the run again (rounded through the device representation, the
! residency boundary the design brief names) and the walk's interludes proceed on the
! CPU. The next FEL element re-uploads.
!-

subroutine fel_device_element_end (dev, beam, ff)

type (fel_device_struct) dev
type (fel_beam_struct) beam
type (fel_field_struct) ff(:)

!

if (.not. dev%resident) return
call fel_device_readback (dev, beam, ff)
dev%resident = .false.

end subroutine fel_device_element_end

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_device_slice_energy (dev, im, ifld, dx) result (energy)
!
! Routine to read one resident field slice's energy-like sum, sum|E|^2 * dx^2/(2 Z0)
! in FP64 over the FP32 record, both planes when two are live: the escape accounting
! the device-resident slippage needs (fel_device_apply_slippage in fel_track_mod owns
! the bookkeeping; this module owns the seam crossing). im is the member and ifld the
! 1-based record index.
!-

function fel_device_slice_energy (dev, im, ifld, dx) result (energy)

type (fel_device_struct) dev
integer im, ifld, ip
real(rp) dx, energy

!

energy = 0
do ip = 1, dev%npol
  call luc_dev_download_field_slice (im-1, ip-1, ifld-1, dev%se)
  energy = energy + sum(real(real(dev%se, sp), rp)**2 + real(aimag(dev%se), rp)**2)
enddo
energy = energy * dx**2 / (2 * (mu_0_vac * c_light))

end function fel_device_slice_energy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_zero_slice (dev, im, ifld)
!
! Routine to zero one resident field slice of member im, every plane, the
! non-periodic fill of the slippage rotation. ifld is the 1-based record index.
!-

subroutine fel_device_zero_slice (dev, im, ifld)

type (fel_device_struct) dev
integer im, ifld, ip

!

do ip = 1, dev%npol
  call luc_dev_zero_field_slice (im-1, ip-1, ifld-1)
enddo

end subroutine fel_device_zero_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_twin_begin (dev, fp32, beam, wf, s0, par, phi0, phi0_new, ks, err_flag)
!
! Routine to put the device in the twin's role for one step, part one: uploads and
! the step encode, called after the FP64 particle advance and BEFORE the production
! field solve (the lockstep field image must round from the pre-solve record). In
! lockstep the shared state is the caller's pre-step snapshot s0 and the reference
! moves to the slice mean; in freerun the resident state carries and only the first
! step uploads. The encoded step then runs concurrently with the production solve.
! The propagators are the caller's, set through fel_device_set_kernel before this.
!
! Input:
!   ff(:)     -- fel_field_struct: The field set, in set order.
!   s0(:,:,:) -- real(rp): Pre-step FP64 state, (6, npart, nslice) as x, px, y, py,
!                  z, pz, snapshotted before the leading transverse half step.
!-

subroutine fel_device_twin_begin (dev, fp32, beam, ff, s0, par, phi0, phi0_new, ks, err_flag)

type (fel_device_struct) dev
type (fel_fp32_struct) fp32
type (fel_beam_struct) beam
type (fel_field_struct) ff(:)
type (fel_device_par_struct) par
real(rp) s0(:,:,:), phi0, phi0_new, ks
logical err_flag

real(rp) p0_mc, gamma0, p_mc, gam, beta, delta
integer is, ip, n, nslice
logical upload

!

err_flag = .false.
p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)
nslice = size(beam%slice)

upload = .not. fp32%freerun .or. .not. dev%twin_live

if (upload) then
  do is = 1, nslice
    n = beam%slice(is)%n
    dev%z_ref(is) = 0
    if (n > 0) dev%z_ref(is) = sum(s0(5, 1:n, is)) / n
    do ip = 1, n
      dev%sx(ip) = real(s0(1, ip, is), c_float)
      dev%spx(ip) = real(s0(2, ip, is), c_float)
      dev%sy(ip) = real(s0(3, ip, is), c_float)
      dev%spy(ip) = real(s0(4, ip, is), c_float)
      dev%sw(ip) = real(beam%slice(is)%weight(ip), c_float)
      p_mc = p0_mc * (1 + s0(6, ip, is))
      gam = sqrt(p_mc**2 + 1)
      beta = p_mc / gam
      dev%sg(ip) = real(gam - gamma0, c_float)
      delta = ks * s0(5, ip, is) / beta - ks * dev%z_ref(is)
      dev%su(ip) = nint(delta * fel_dev_ticks_per_rad$, c_int64_t)

      ! The check's mutation reaches the shared state here in lockstep, the residual
      ! coarsening fel_fp32_twin_slice applies to dzr, in ticks: without it the qres
      ! perturbation alone would carry falsifiability, and the goal is two hooks
      ! (state and kernel), each of which must move a recorded level.

      if (fp32%mutate) dev%su(ip) = 65536_c_int64_t * (dev%su(ip) / 65536_c_int64_t)
    enddo
    do ip = n+1, dev%npart
      dev%sx(ip) = 0;  dev%spx(ip) = 0;  dev%sy(ip) = 0;  dev%spy(ip) = 0
      dev%sg(ip) = 0;  dev%sw(ip) = 0;  dev%su(ip) = 0
    enddo
    dev%su0(:, is) = dev%su
    call luc_dev_upload_slice (is-1, dev%npart, dev%sx, dev%spx, dev%sy, dev%spy, &
                               dev%sg, dev%su, dev%sw)
  enddo

  ! The field images, every plane of every member, rounded from the pre-solve FP64
  ! records in ring order.

  call fel_device_upload_fields (dev, ff)
  dev%twin_live = .true.
endif

call fel_device_step_run (dev, par, phi0, phi0_new, ks, err_flag)

end subroutine fel_device_twin_begin

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_twin_rows (dev, fp32, beam, ff, par, phi0, phi0_new, ks)
!
! Routine to put the device in the twin's role, part two: after the production solve,
! read the device state back and fill fel_fp32_mod's rows against the FP64 path.
! Rows 1 to 4 carry the device's FP32 transverse maps (the CPU twin refreshed those
! from FP64, so its recorded levels there were pure rounding; the device's are
! arithmetic and get their own recorded levels). Rows 5 to 7 follow
! fel_fp32_twin_slice's conventions with the phase compared in the tick chart's FP64
! image; rows 8 and 9 follow fel_fp32_field_twin's, against a fresh FP64 deposit of
! the shared post-step state. The guard statistic is the median per-step phase
! increment in ticks (uniform quantum, so no spacing division). The device's own
! phase accumulators are kept as the next step's guard baseline, which is what makes
! the freerun guard cost no extra sync.
!
! The field set: the phasor, source and field rows are each the worst over the set's
! members, so the stream keeps its ten columns and its recorded ceilings, and every
! member's own worst is kept in dev%worst_h for the footer's attribution. A member's
! phasor and deposit are compared at h times the fundamental phase, fel_field_step's
! harm*theta, through fel_fp32_deposit64 unchanged (h*theta is its theta, the
! member's scale its scl_w). Every member's source and field rows are normalized by
! the fundamental's post-solve field norm over its planes, the run's power scale. A
! harmonic's own norm cannot serve: on a quiet start the CPU cancels the harmonic
! bunching to FP64 roundoff while the device's FP32 phase leaves it at ~1e-7 of the
! charge, and a dark member's row divided by its own roundoff-level field measured
! that noise against nothing (the first attempt read 6 to 8 on a deposit whose
! end-to-end power agreed at 1e-3). The harmonic's fidelity lives in its phasor row,
! which needs no field normalization, and in the end-to-end harmonic power against
! the CPU that check_device.py holds. With one member this is the normalization the
! single-field rows always had.
!
! Input:
!   ff(:) -- fel_field_struct: The field set, in set order, post-solve. The ring
!              offset is ff(1)%slip%first, one for the set.
!-

subroutine fel_device_twin_rows (dev, fp32, beam, ff, par, phi0, phi0_new, ks)

type (fel_device_struct) dev
type (fel_fp32_struct), target :: fp32
type (fel_beam_struct), target :: beam
type (fel_field_struct) ff(:)
type (fel_device_par_struct) par
real(rp) phi0, phi0_new, ks
integer first

type (fel_slice_struct), pointer :: sl
complex(rp), allocatable :: s64(:,:), e64(:,:)
complex(sp), allocatable :: sdev(:,:), edev(:,:)
complex(rp) p32sum(dev%nfield), p64sum(dev%nfield), e_ip
real(rp), allocatable :: incr(:)
real(rp) p0_mc, gamma0, gam, p_mc, beta, pzd, delta, theta_rel, phi_ref, h, th_h
real(rp) x, y, awloc_d, awloc_c, wsum, enorm, w_s, w_f, tick_med, pw64, pw32, v
real(rp) dstat(fel_fp32_nq$), sc(fel_fp32_nq$), d7, d8, d9
integer(c_int64_t) cret
integer is, ip, n, ng, ifld, ix, iy, im, ipl

!

p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)
ng = size(ff(1)%wf%Ex, 1)
first = ff(1)%slip%first
cret = nint((phi0_new - phi0) * fel_dev_ticks_per_rad$, c_int64_t)
allocate (s64(ng, ng), e64(ng, ng), sdev(ng, ng), edev(ng, ng), incr(dev%npart))

do is = 1, size(beam%slice)
  sl => beam%slice(is)
  n = sl%n
  ifld = 1 + mod(is - 1 + first, size(ff(1)%wf%Ex, 3))   ! fel_field_index's mapping.

  call luc_dev_download_slice (is-1, dev%npart, dev%sx, dev%spx, dev%sy, dev%spy, &
                               dev%sg, dev%su)

  dstat = 0
  p32sum = 0
  p64sum = 0
  wsum = 0
  do ip = 1, n
    x = real(dev%sx(ip), rp)
    y = real(dev%sy(ip), rp)
    gam = gamma0 + real(dev%sg(ip), rp)
    p_mc = sqrt(gam**2 - 1)
    pzd = (p_mc - p0_mc) / p0_mc
    delta = real(dev%su(ip), rp) * fel_dev_rad_per_tick$

    ! The step's own phase advance in ticks: the accumulator moved by the RK
    ! increment minus the common cret, so adding cret back isolates the physics.

    incr(ip) = abs(real(dev%su(ip) - dev%su0(ip, is) + cret, rp))

    dstat(1) = max(dstat(1), abs(x - sl%x(ip)))
    dstat(2) = max(dstat(2), abs(real(dev%spx(ip), rp) - sl%px(ip)))
    dstat(3) = max(dstat(3), abs(y - sl%y(ip)))
    dstat(4) = max(dstat(4), abs(real(dev%spy(ip), rp) - sl%py(ip)))
    dstat(5) = max(dstat(5), abs(pzd - sl%pz(ip)))

    ! The CPU phase in the same chart: theta relative to the new reference.

    p_mc = p0_mc * (1 + sl%pz(ip))
    beta = p_mc / sqrt(p_mc**2 + 1)
    theta_rel = ks * sl%z(ip) / beta - ks * dev%z_ref(is)
    dstat(6) = max(dstat(6), abs(delta - theta_rel))

    ! The source phasor of every member, both sides in FP64 from their own states at
    ! that member's phase h*theta, the twin's full-charge normalization.

    awloc_d = fel_fp32_faw (x, y, par%kx, par%ky, par%ax, par%ay, par%cos_t, par%sin_t)
    awloc_c = fel_fp32_faw (sl%x(ip), sl%y(ip), par%kx, par%ky, par%ax, par%ay, &
                            par%cos_t, par%sin_t)
    do im = 1, dev%nfield
      h = par%harm(im)
      phi_ref = h * (phi0_new + ks * dev%z_ref(is))
      e_ip = cmplx(cos(phi_ref + h * delta), -sin(phi_ref + h * delta), rp)
      p32sum(im) = p32sum(im) + sl%weight(ip) * awloc_d * e_ip
      e_ip = cmplx(cos(phi_ref + h * theta_rel), -sin(phi_ref + h * theta_rel), rp)
      p64sum(im) = p64sum(im) + sl%weight(ip) * awloc_c * e_ip
    enddo
    wsum = wsum + sl%weight(ip) * awloc_c
  enddo
  dev%su0(1:dev%npart, is) = dev%su(1:dev%npart)

  sc(1) = maxval(abs(sl%x(1:n))) + 1e-30_rp
  sc(2) = maxval(abs(sl%px(1:n))) + 1e-30_rp
  sc(3) = maxval(abs(sl%y(1:n))) + 1e-30_rp
  sc(4) = maxval(abs(sl%py(1:n))) + 1e-30_rp
  sc(5) = maxval(abs(sl%pz(1:n))) + 1e-30_rp
  sc(6) = 1
  sc(7) = wsum + 1e-30_rp

  fp32%div_slice(1:6, is) = dstat(1:6) / sc(1:6)
  d7 = 0
  do im = 1, dev%nfield
    v = abs(p32sum(im) - p64sum(im)) / sc(7)
    d7 = max(d7, v)
    dev%worst_h(1, im) = max(dev%worst_h(1, im), v)
  enddo
  fp32%div_slice(7, is) = d7
  fp32%bmag64(is) = abs(p64sum(1)) / sc(7)
  fp32%bmag32(is) = abs(p32sum(1)) / sc(7)

  call fel_fp32_median (incr, n, tick_med)
  fp32%ulp_slice(is) = tick_med

  ! The field rows per member: the device's own source grid and post-solve planes
  ! against a fresh FP64 deposit of the shared post-step state at h*theta and the
  ! post-solve records, fel_fp32_field_twin's conventions (the source counted times 2
  ! as the field sees it), every member normalized by the fundamental's post-solve
  ! field norm over its planes. The end-to-end power sums are the fundamental's.

  enorm = 0
  do ipl = 1, dev%npol
    if (ipl == 1) then
      e64 = ff(1)%wf%Ex(:,:,ifld)
    else
      e64 = ff(1)%wf%Ey(:,:,ifld)
    endif
    enorm = enorm + sum(real(e64, rp)**2 + aimag(e64)**2)
  enddo
  enorm = sqrt(enorm) + 1e-30_rp

  d8 = 0
  d9 = 0
  pw64 = 0
  pw32 = 0
  do im = 1, dev%nfield
    h = par%harm(im)
    call luc_dev_download_source_slice (im-1, ifld-1, sdev)

    s64 = 0
    do ip = 1, n
      p_mc = p0_mc * (1 + sl%pz(ip))
      gam = sqrt(p_mc**2 + 1)
      beta = p_mc / gam
      th_h = h * (beam%phi0 + ks * sl%z(ip) / beta)
      call fel_fp32_deposit64 (s64, sl%x(ip), sl%y(ip), th_h, sl%weight(ip) / gam, &
                               par%scl_w(im), par%gridmax, par%dgrid, &
                               par%kx, par%ky, par%ax, par%ay, par%cos_t, par%sin_t)
    enddo

    w_s = 0
    do iy = 1, ng
      do ix = 1, ng
        w_s = w_s + abs(cmplx(sdev(ix,iy), kind=rp) - s64(ix,iy))**2
      enddo
    enddo

    w_f = 0
    do ipl = 1, dev%npol
      call luc_dev_download_field_slice (im-1, ipl-1, ifld-1, edev)
      if (ipl == 1) then
        e64 = ff(im)%wf%Ex(:,:,ifld)
      else
        e64 = ff(im)%wf%Ey(:,:,ifld)
      endif
      do iy = 1, ng
        do ix = 1, ng
          w_f = w_f + abs(cmplx(edev(ix,iy), kind=rp) - e64(ix,iy))**2
        enddo
      enddo
      if (im == 1) then
        pw64 = pw64 + sum(real(e64, rp)**2 + aimag(e64)**2)
        pw32 = pw32 + sum(real(real(edev, sp), rp)**2 + real(aimag(edev), rp)**2)
      endif
    enddo

    v = sqrt(w_s) * 2 / enorm
    d8 = max(d8, v)
    dev%worst_h(2, im) = max(dev%worst_h(2, im), v)
    v = sqrt(w_f) / enorm
    d9 = max(d9, v)
    dev%worst_h(3, im) = max(dev%worst_h(3, im), v)
  enddo
  fp32%div_slice(8, is) = d8
  fp32%div_slice(9, is) = d9
  fp32%pow64(is) = pw64
  fp32%pow32(is) = pw32
enddo

end subroutine fel_device_twin_rows

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_close_run (dev, fp32)
!
! Routine to report what the device did and release it: the step count and the
! device-busy seconds from the backend's own command-buffer timestamps, the honest
! pair a wall clock is judged against (the dispatch floor is their difference). When
! the device held the twin's role over a set of more than one member, its per-member
! worst levels go into the instrument's footer first, since the ten-column stream
! carries only the worst over members. The caller closes that stream afterwards.
!-

subroutine fel_device_close_run (dev, fp32)

type (fel_device_struct) dev
type (fel_fp32_struct) fp32
integer im
character(*), parameter :: r_name = 'fel_device_close_run'

!

if (.not. dev%on) return
call luc_dev_sync ()

if (fp32%on .and. fp32%iu /= 0 .and. dev%nfield > 1) then
  write (fp32%iu, '(a)') '# per member, the device twin''s worst over the run:'
  do im = 1, dev%nfield
    write (fp32%iu, '(a, i0, a, es13.4)') 'worst_phasor_h', dev%harm(im), ' ', dev%worst_h(1, im)
    write (fp32%iu, '(a, i0, a, es13.4)') 'worst_source_h', dev%harm(im), ' ', dev%worst_h(2, im)
    write (fp32%iu, '(a, i0, a, es13.4)') 'worst_field_h', dev%harm(im), ' ', dev%worst_h(3, im)
  enddo
endif

call out_io (s_info$, r_name, 'Device: \i0\ steps encoded, device busy \f10.3\ s.', &
             i_array = [dev%nstep], r_array = [real(luc_dev_seconds(), rp)])
call luc_dev_close ()

end subroutine fel_device_close_run

end module fel_device_mod
