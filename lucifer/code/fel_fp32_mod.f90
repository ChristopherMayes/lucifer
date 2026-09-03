!+
! Module fel_fp32_mod
!
! The single-precision particle path and the lockstep instrument that prices it
! (doc/performance.md, doc/validation.md). Apple GPUs have no FP64, so before any
! kernel is written a single-precision form of the averaged FEL advance must exist
! with its divergence from the FP64 path a measured number. Nothing outside the
! instrument reads the FP32 state, and the FP64 path is untouched: this module is a
! measuring device, not a production mode.
!
! The FP32 state is the packed six-vector (x, px, y, py, dzr, pz), 24 bytes per
! particle. dzr is the longitudinal residual of the split z = z_ref + dzr, with z_ref
! one FP64 number per slice. Storing Bmad's z absolutely in FP32 fails silently (the
! per-step increment rounds to zero and bunching never forms), so the residual form is
! load-bearing and the runtime guard below watches it. The reference must also move:
! the beam's common z drift is secular, a static reference lets the residual grow until
! its ulp swallows the physics, and the guard caught exactly that on this instrument's
! first run (step 514 of the steady-state example, 28 ulps against the floor of 32).
! z_ref therefore refreshes to the slice's FP64 mean each step, one host number per
! slice, and in freerun mode the persistent residuals rebase across the reference move
! by fel_fp32_renorm, which is the migration operation doing double duty.
!
! Reformulations, each forced by an FP32 quantum and each below FP64 notice:
!
!   energy   The working variable is goff = gamma - gamma0. The conversion is
!            goff = pz * p0_mc^2/gamma0, exact to O(1/gamma^2) ~ 4e-9, below FP32
!            resolution. Absolute FP32 gamma has a quantum of 1.35e-3 at these
!            energies against a per-step change of 3.0e-4, and sqrt(gamma^2 - 1)
!            loses the 1 entirely (the quantum of gamma^2 is about 15).
!   phase    The working variable is delta = theta - phi0, formed as ks * dzr
!            (the 1/beta corrections are ~4e-9 of theta, below FP32 resolution and
!            part of what the instrument measures). The phase factor is
!            e^{-i theta} = e^{-i phi0} * e^{-i delta}: the base rotator is one FP64
!            evaluation per slice, the per-particle sincos is FP32 on the small
!            angle. This is the angle-addition form a device kernel would use.
!   detuning The FP64 ODE forms sqrt(1 - btper0/gamma^2), which is 1 minus ~1.3e-8:
!            in FP32 that is exactly 1 and the phase equation loses the energy
!            dependence completely. The FP32 ODE forms the detuning as the
!            difference 0.5*ks*(qres - q) with q = btper0/gamma^2 and
!            qres = 2*ku/ks, so the resonant cancellation happens between two
!            small like quantities. The remaining floor is btper0's own FP32
!            representation, and the recorded levels carry it.
!
! The instrument runs both precisions from a shared state at fel_advance's own
! sequence point: the FP64 slice is copied before its
! advance, the twin advances the FP32 image of that copy, and the comparison reads
! both results before anything else touches the slice. In lockstep mode the FP32
! state is rebuilt from FP64 every step, so a wrong formula shows as a jump rather
! than as growth. In freerun mode dzr, pz and the FP32 field record persist across
! steps, the deposit feeds the twin's own field and the gather reads it, so the twin is
! a complete single-precision run beside the FP64 one; the transverse coordinates are
! refreshed from FP64 in both modes, since the transverse maps stay FP64 and enter
! rounded. A freerun window is one slice by refusal: the twin keeps no slippage
! rotation of its own, and slippage moves light between slices. The source phasor sum(w * awloc * e^{-i theta}) is compared at the
! same point; the production deposit runs one FP64 transverse half-step later, which
! is precision-neutral between the sides.
!
! Everything here is serial-safe inside the caller's parallel slice loop: the twin
! touches only its slice's state and its slice's row of the divergence arrays, so the
! instrument's own output is identical at any thread count.
!-

module fel_fp32_mod

use fel_beam_mod
use wavefront_mod

implicit none

integer, parameter :: fel_fp32_nq$ = 9   ! x, px, y, py, pz, theta, phasor, source, field.

! The field twin transforms with FFTW's single-precision interface, and not every
! toolchain carries it: the conda environment does, and the off-site distribution
! builds FFTW from source in double precision only. The build detects the library and
! defines LUCIFER_HAVE_FFTW3F when it is there (lucifer/CMakeLists.txt). Without it the
! instrument refuses at setup rather than transforming in another precision, which
! would measure something other than what it reports.

#ifdef LUCIFER_HAVE_FFTW3F
logical, parameter :: fel_fp32_have_fftw3f$ = .true.
#else
logical, parameter :: fel_fp32_have_fftw3f$ = .false.
#endif

character(8), parameter :: fel_fp32_qname(fel_fp32_nq$) = &
      [character(8) :: 'x', 'px', 'y', 'py', 'pz', 'theta', 'phasor', 'source', 'field']

!+
! Struct fel_fp32_slice_struct
!
! One slice's packed FP32 particle set: the six-vector plus the working longitudinal
! chart the twin evolves. Allocated to the slice's fill count at first use.
!-

type fel_fp32_slice_struct
  real(sp), allocatable :: x(:), px(:), y(:), py(:)   ! Transverse, Bmad conventions.
  real(sp), allocatable :: dzr(:)                     ! z residual off z_ref [m].
  real(sp), allocatable :: pz(:)                      ! (p - p0)/p0, already an offset.
  real(rp) :: z_ref = 0                               ! The slice's FP64 z reference [m].
end type

!+
! Struct fel_fp32_struct
!
! The instrument's run state: mode, the per-slice FP32 sets, the per-step divergence
! table the parallel loop fills and the serial epilogue reduces, and the guard and
! worst-case accumulators the footer and the check read.
!-

type fel_fp32_struct
  logical :: on = .false.
  logical :: freerun = .false.       ! dzr/pz/field persist across steps (compounding measured).
  logical :: mutate = .false.        ! Check hook: truncate dzr harder so the check can fail.
  integer :: iu = 0                  ! The .fp32.txt stream.
  integer :: istep = 0               ! Steps instrumented.
  type (fel_fp32_slice_struct), allocatable :: sl32(:)
  real(rp), allocatable :: div_slice(:,:)    ! (nq, nslice): this step, filled in parallel.
  real(rp), allocatable :: ulp_slice(:)      ! (nslice): median per-step |d delta| in ulps.
  real(rp) :: worst(fel_fp32_nq$) = 0        ! Run-level worst per quantity.
  real(rp) :: ulp_min = huge(1.0_rp)         ! Run-level worst (smallest) guard statistic.
  logical :: checked = .false.               ! A configuration check ran (first step).
  ! The field twin: one FP32 field record per beam slice. In lockstep it is the FP64
  ! record rounded each step before the solve; in freerun it carries, fed by its own
  ! deposit, and the particle twin gathers from it.
  complex(sp), allocatable :: e32(:,:,:)     ! (ngrid, ngrid, nslice).
  complex(sp), allocatable :: k32(:,:)       ! The propagator, rounded from the FP64 kernel.
  real(rp) :: k32_key(4) = -1                ! (ngrid, dgrid, ks, dz) the rounding matches.
  real(rp), allocatable :: pow32(:), pow64(:)   ! Post-solve field power sums, this step.
  real(rp), allocatable :: bmag32(:), bmag64(:) ! |phasor|/charge, this step (bunching).
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_setup (fp32, mode, mutate, nslice, out_root, err_flag)
!
! Routine to arm the instrument from the input knob: '' or 'off' leaves it dark,
! 'lockstep' rebuilds the FP32 state from FP64 every step, 'freerun' lets the
! longitudinal state carry across steps. Opens the .fp32.txt stream.
!
! Input:
!   mode     -- character(*): global%fp32_check.
!   mutate   -- logical: global%fp32_mutate, the check's own failure hook.
!   nslice   -- integer: Slice count.
!   out_root -- character(*): Output root for the stream file.
!
! Output:
!   fp32     -- fel_fp32_struct: Armed (or left dark).
!   err_flag -- logical: Set True on an unrecognized mode. False otherwise.
!-

subroutine fel_fp32_setup (fp32, mode, mutate, nslice, bucket_shift, out_root, err_flag)

type (fel_fp32_struct) fp32
character(*) mode, out_root
logical mutate, err_flag
integer nslice
real(rp) bucket_shift
real(sp), allocatable :: rres(:), rorig(:)
real(rp) renorm_worst
integer i, nprobe
character(*), parameter :: r_name = 'fel_fp32_setup'

!

err_flag = .false.

select case (mode)
case ('', 'off')
  return
case ('lockstep')
  fp32%freerun = .false.
case ('freerun')
  fp32%freerun = .true.
case default
  call out_io (s_error$, r_name, 'UNRECOGNIZED FP32_CHECK MODE: "' // trim(mode) // '".', &
                                 'THE MODES ARE off, lockstep AND freerun.')
  err_flag = .true.
  return
end select

! The field twin's transform is FFTW's own single-precision interface, so a build
! without that library cannot run the instrument in either role and says so here.

if (.not. fel_fp32_have_fftw3f$) then
  call out_io (s_error$, r_name, 'FP32_CHECK NEEDS THE SINGLE-PRECISION FFTW (libfftw3f), AND THIS BUILD', &
        'DOES NOT CARRY IT. THE FIELD TWIN TRANSFORMS WITH FFTW''S OWN fftwf INTERFACE,', &
        'AND ANOTHER PRECISION WOULD MEASURE SOMETHING OTHER THAN WHAT IT REPORTS.')
  err_flag = .true.
  return
endif

! Freerun now carries the FP32 field, and the field twin keeps one record per beam
! slice with no slippage bookkeeping of its own, so a multi-slice freerun would hold
! each slice's light fixed where the real dynamics rotate it across slices.

if (fp32%freerun .and. nslice > 1) then
  call out_io (s_error$, r_name, 'FP32_CHECK = "freerun" CARRIES THE FP32 FIELD AND COVERS A', &
        'SINGLE-SLICE WINDOW ONLY: SLIPPAGE MOVES LIGHT BETWEEN SLICES AND THE TWIN', &
        'KEEPS NO ROTATION OF ITS OWN. USE "lockstep" FOR A TIME-DEPENDENT WINDOW.')
  err_flag = .true.
  return
endif

fp32%on = .true.
fp32%mutate = mutate
allocate (fp32%sl32(nslice))
allocate (fp32%div_slice(fel_fp32_nq$, nslice))
allocate (fp32%ulp_slice(nslice))
allocate (fp32%pow32(nslice), fp32%pow64(nslice), fp32%bmag32(nslice), fp32%bmag64(nslice))
fp32%pow32 = 0;  fp32%pow64 = 0;  fp32%bmag32 = 0;  fp32%bmag64 = 0

open (newunit = fp32%iu, file = trim(out_root) // '.fp32.txt', action = 'write')
write (fp32%iu, '(a)') '# FP32 lockstep instrument: per-step worst relative divergence per quantity,'
write (fp32%iu, '(a)') '# FP32 twin against the FP64 path from a shared state at fel_advance. theta is'
write (fp32%iu, '(a)') '# absolute [rad], guard_ulp is the median per-step residual phase increment in'
write (fp32%iu, '(a)') '# ulps of the residual (the silent-z guard: small means FP32 cannot resolve the'
write (fp32%iu, '(a)') '# step and the run refuses).'
write (fp32%iu, '(a, l1)') '# freerun = ', fp32%freerun

! Migration in the residual representation is a bucket renormalization: a mover's
! residual re-references to the new slice origin by the FP64 bucket shift. The
! round trip (shift out and back) is measured here on residuals spanning the slice,
! in ulps of the residual, and recorded beside the run's own levels.

nprobe = 4097
allocate (rres(nprobe), rorig(nprobe))
do i = 1, nprobe
  rres(i) = real((real(i-1, rp) / (nprobe-1) - 0.5_rp) * bucket_shift, sp)
enddo
rorig = rres
call fel_fp32_renorm (rres, bucket_shift)
call fel_fp32_renorm (rres, -bucket_shift)
renorm_worst = 0
do i = 1, nprobe
  renorm_worst = max(renorm_worst, real(abs(rres(i) - rorig(i)), rp) / &
                     max(real(spacing(max(abs(rorig(i)), real(bucket_shift, sp))), rp), tiny(1.0_rp)))
enddo
write (fp32%iu, '(a, es13.4, a, es13.4, a)') '# renorm_roundtrip_ulp ', renorm_worst, &
      '   (bucket shift ', bucket_shift, ' m)'
write (fp32%iu, '(a)') '#   step          x            px           y            py           pz         theta        phasor       source        field      guard_ulp'

end subroutine fel_fp32_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_renorm (dzr, shift)
!
! Routine to re-reference FP32 residuals across a bucket boundary: the migration
! operation in the split representation. The shift is FP64 (it is the slice bucket,
! exactly beta * slice_spacing) and rounds once into the residual's own precision.
!
! Input:
!   dzr(:) -- real(sp): Residuals off the old reference.
!   shift  -- real(rp): The bucket shift [m], signed.
!
! Output:
!   dzr(:) -- real(sp): Residuals off the new reference.
!-

subroutine fel_fp32_renorm (dzr, shift)

real(sp) dzr(:)
real(rp) shift

dzr = dzr - real(shift, sp)

end subroutine fel_fp32_renorm

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_fp32_faw (x, y, kx, ky, ax, ay, cos_t, sin_t) result (value)
!
! The FP64 comparison side's transverse roll-off, faw of fel_track_mod on loose
! scalars rather than on an fel_und_struct, which the instrument and the device seam
! do not carry. Shared by both of the instrument's roles: the CPU twin compares
! against it and the device twin does the same, and one transcription of the
! expression is the point.
!
! Input:
!   x, y            -- real(rp): Transverse position [m].
!   kx, ky          -- real(rp): Natural-focusing roll-off [1/m^2].
!   ax, ay          -- real(rp): Undulator field offset [m].
!   cos_t, sin_t    -- real(rp): Wiggle-plane tilt.
!
! Output:
!   value           -- real(rp): aw(x,y)/aw, the rolled-off factor.
!-

function fel_fp32_faw (x, y, kx, ky, ax, ay, cos_t, sin_t) result (value)

real(rp) x, y, kx, ky, ax, ay, cos_t, sin_t
real(rp) value, ddx, ddy, ddt

!

ddx = x - ax
ddy = y - ay
if (sin_t /= 0) then
  ddt = cos_t * ddx + sin_t * ddy
  ddy = -sin_t * ddx + cos_t * ddy
  ddx = ddt
endif
value = 1 + 0.5_rp * (kx * ddx*ddx + ky * ddy*ddy)

end function fel_fp32_faw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_deposit64 (src, x, y, theta, wgt, scl_w, gridmax, dgrid, &
!                                   kx2, ky2, ax, ay, cos_t, sin_t)
!
! One particle's FP64 source deposit, fel_field_step's own expressions:
! part = sqrt(faw2) * scl_w * wgt with faw2 = 1 + kx*dx^2 + ky*dy^2 (no half, which
! is Genesis's own roll-off transcribed), cpart = (sin theta + i cos theta) * part,
! bilinear scatter onto the lower-left cell and its three neighbours.
!
! This is the reference the FP32 source row is measured against, and both roles of the
! instrument need it: the CPU field twin and the device twin. It lives here, once,
! because it is a transcription of a convention -- a second copy would have to be
! found and changed the day fel_field_step's deposit changes.
!
! Input:
!   src(:,:)        -- complex(rp): The accumulating source grid.
!   x, y            -- real(rp): Transverse position [m].
!   theta           -- real(rp): Ponderomotive phase [rad].
!   wgt             -- real(rp): Charge weight over gamma [C].
!   scl_w           -- real(rp): The deposit scale of fel_field_step.
!   gridmax, dgrid  -- real(rp): Grid half width and spacing [m].
!   kx2, ky2        -- real(rp): faw2's roll-off coefficients [1/m^2].
!   ax, ay          -- real(rp): Undulator field offset [m].
!   cos_t, sin_t    -- real(rp): Wiggle-plane tilt.
!
! Output:
!   src(:,:)        -- complex(rp): This particle's contribution added.
!-

subroutine fel_fp32_deposit64 (src, x, y, theta, wgt, scl_w, gridmax, dgrid, &
                               kx2, ky2, ax, ay, cos_t, sin_t)

complex(rp) src(:,:)
real(rp) x, y, theta, wgt, scl_w, gridmax, dgrid, kx2, ky2, ax, ay, cos_t, sin_t
real(rp) f2, ppart, wwx, wwy, sth, cth, ddx, ddy, ddt
complex(rp) cpart
integer jx, jy

!

if (.not. (x > -gridmax .and. x < gridmax .and. y > -gridmax .and. y < gridmax)) return
wwx = (x + gridmax) / dgrid
wwy = (y + gridmax) / dgrid
jx = int(floor(wwx));  jy = int(floor(wwy))
wwx = 1 + real(jx, rp) - wwx
wwy = 1 + real(jy, rp) - wwy
jx = jx + 1;  jy = jy + 1
ddx = x - ax;  ddy = y - ay
if (sin_t /= 0) then
  ddt = cos_t * ddx + sin_t * ddy
  ddy = -sin_t * ddx + cos_t * ddy
  ddx = ddt
endif
f2 = 1 + kx2 * ddx*ddx + ky2 * ddy*ddy
ppart = sqrt(f2) * scl_w * wgt
sth = sin(theta);  cth = cos(theta)
cpart = cmplx(sth, cth, rp) * ppart
src(jx,   jy)   = src(jx,   jy)   + (wwx * wwy) * cpart
src(jx+1, jy)   = src(jx+1, jy)   + ((1-wwx) * wwy) * cpart
src(jx,   jy+1) = src(jx,   jy+1) + (wwx * (1-wwy)) * cpart
src(jx+1, jy+1) = src(jx+1, jy+1) + ((1-wwx) * (1-wwy)) * cpart

end subroutine fel_fp32_deposit64

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_median (v, nn, med)
!
! The median of the first nn entries, scratch reordered in place. Heapsort:
! O(n log n), no recursion, no degenerate case. This runs per slice per step in both
! of the instrument's roles, so it must be cheap and it must terminate on any input
! (a quickselect written here first did not, and hung a run). The sift is inlined
! twice because a contained procedure cannot hold one of its own.
!
! Input:
!   v(:)  -- real(rp): The values. Reordered in place.
!   nn    -- integer: How many of them count.
!
! Output:
!   v(:)  -- real(rp): Reordered.
!   med   -- real(rp): The median.
!-

subroutine fel_fp32_median (v, nn, med)

real(rp) v(:), med, tmp
integer nn, i, iend, root, child

!

if (nn < 1) then
  med = 0
  return
endif

do i = nn/2, 1, -1
  root = i
  do while (2*root <= nn)
    child = 2*root
    if (child < nn) then
      if (v(child+1) > v(child)) child = child + 1
    endif
    if (v(root) >= v(child)) exit
    tmp = v(root);  v(root) = v(child);  v(child) = tmp
    root = child
  enddo
enddo

do iend = nn, 2, -1
  tmp = v(1);  v(1) = v(iend);  v(iend) = tmp
  root = 1
  do while (2*root <= iend - 1)
    child = 2*root
    if (child < iend - 1) then
      if (v(child+1) > v(child)) child = child + 1
    endif
    if (v(root) >= v(child)) exit
    tmp = v(root);  v(root) = v(child);  v(child) = tmp
    root = child
  enddo
enddo

if (mod(nn, 2) == 1) then
  med = v((nn+1)/2)
else
  med = 0.5_rp * (v(nn/2) + v(nn/2+1))
endif

end subroutine fel_fp32_median

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_twin_slice (fp32, is, x0, y0, px0, py0, z0, pz0, sl, beam, &
!                                    exfld, dx, dy, gridmax, aw, ku, kx, ky, ax, ay, &
!                                    cos_t, sin_t, rtmp, delz, xks, xku, phi0, phi0_new)
!
! Routine to advance one slice's FP32 twin from the shared pre-advance state and fill
! this slice's divergence row. The x0.. arrays are the FP64 state exactly as
! fel_advance received it, sl is the slice after the FP64 advance, and the twin
! mirrors fel_advance's arithmetic under the module header's reformulations:
! fundamental field only, no collective terms (the setup refusals hold that).
!
! Runs inside the caller's parallel loop: everything written is indexed by is.
!-

subroutine fel_fp32_twin_slice (fp32, is, x0, y0, px0, py0, z0, pz0, sl, beam, &
                                exfld, dx, dy, gridmax, aw, ku, kx, ky, ax, ay, &
                                cos_t, sin_t, rtmp, delz, xks, xku, phi0, phi0_new)

type (fel_fp32_struct), target :: fp32
type (fel_fp32_slice_struct), pointer :: t
type (fel_slice_struct) sl
type (fel_beam_struct) beam
integer is
real(rp) x0(:), y0(:), px0(:), py0(:), z0(:), pz0(:)
complex(rp) exfld(:,:)     ! The FP64 field record. The gather reads the FP32 record
                           !   fp32%e32(:,:,is), which lockstep rounded from this one.
real(rp) dx, dy, gridmax, aw, ku, kx, ky, ax, ay, cos_t, sin_t, rtmp, delz, xks, xku
real(rp) phi0, phi0_new

! FP32 constants of the step, converted once.

real(sp) aw32, ku32, kx32, ky32, ax32, ay32, ct32, st32, rt32, dz32, ks32, gam032
real(sp) e032, p032, qres32, gmax32, dxg32, dyg32, cret32
complex(sp) base32

! FP32 per-particle working set.

real(sp) x, y, px, py, pz, dzr, goff, delta, gam, awloc, btpar, px_g, py_g
real(sp) wx, wy, q, s_d, c_d, d_new
complex(sp) cpart, rpart, p32sum
complex(rp) p64sum, e_ip
real(rp) wsum
real(rp) gamma0, p0_mc, theta64, beta64, p_mc64, aw64loc, ulp_med, phi_ref, zref_new
real(rp) dstat(fel_fp32_nq$), sc(fel_fp32_nq$), dinc
real(rp), allocatable :: incr(:)
integer ip, n, ix, iy
logical on_grid

!

t => fp32%sl32(is)
n = sl%n
p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)

if (.not. allocated(fp32%e32)) then
  allocate (fp32%e32(size(exfld,1), size(exfld,2), size(fp32%sl32)))
  fp32%e32 = 0
endif
if (.not. fp32%freerun .or. all(fp32%e32(:,:,is) == 0)) then
  fp32%e32(:,:,is) = cmplx(exfld, kind=sp)
endif

if (.not. allocated(t%x)) then
  allocate (t%x(n), t%px(n), t%y(n), t%py(n), t%dzr(n), t%pz(n))
  t%z_ref = sum(z0(1:n)) / n
  t%pz = real(pz0(1:n), sp)
  t%dzr = real(z0(1:n) - t%z_ref, sp)
endif

! The shared state: transverse always from FP64 (the transverse maps stay FP64), the
! longitudinal from FP64 in lockstep mode and carried in freerun mode. The FP64
! reference moves to the slice mean each step in both modes: one host number per
! slice, which is the split's working principle. In freerun the persistent residuals
! rebase across the move, the renormalization the migration path uses.

t%x = real(x0(1:n), sp);   t%px = real(px0(1:n), sp)
t%y = real(y0(1:n), sp);   t%py = real(py0(1:n), sp)
if (fp32%freerun) then
  zref_new = sum(z0(1:n)) / n
  call fel_fp32_renorm (t%dzr, zref_new - t%z_ref)
  t%z_ref = zref_new
else
  t%z_ref = sum(z0(1:n)) / n
  t%pz = real(pz0(1:n), sp)
  t%dzr = real(z0(1:n) - t%z_ref, sp)
endif

! The check's own mutation: a deliberately coarser residual, so the recorded level
! must move and the check is falsifiable.

if (fp32%mutate) t%dzr = anint(t%dzr / (256 * spacing(t%dzr))) * (256 * spacing(t%dzr))

! Step constants. e0 = p0^2/gamma0 converts pz to goff = gamma - gamma0, exact to
! O(1/gamma^2) ~ 4e-9 (module header). qres = 2*ku/ks is the resonance in the
! detuning difference. cret = phi0_new - phi0 is the chart-return base, small.

aw32 = real(aw, sp);  ku32 = real(ku, sp);  kx32 = real(kx, sp);  ky32 = real(ky, sp)
ax32 = real(ax, sp);  ay32 = real(ay, sp);  ct32 = real(cos_t, sp);  st32 = real(sin_t, sp)
rt32 = real(rtmp, sp);  dz32 = real(delz, sp);  ks32 = real(xks, sp)
gam032 = real(gamma0, sp);  p032 = real(p0_mc, sp)
e032 = real(p0_mc**2 / gamma0, sp)
qres32 = real(2 * xku / xks, sp)
gmax32 = real(gridmax, sp);  dxg32 = real(dx, sp);  dyg32 = real(dy, sp)
! The chart return is z_ref-free: theta and dzr are both referenced to phi_ref, so the
! reference cancels and cret stays the small per-step phase advance.

cret32 = real(phi0_new - phi0, sp)
phi_ref = phi0 + xks * t%z_ref
base32 = cmplx(cmplx(cos(phi_ref), -sin(phi_ref), rp), kind=sp)   ! e^{-i(phi0 + ks z_ref)}, FP64 once.

allocate (incr(n))
p32sum = 0
p64sum = 0
wsum = 0
dstat = 0

do ip = 1, n

  x = t%x(ip);  y = t%y(ip);  px = t%px(ip);  py = t%py(ip)
  pz = t%pz(ip);  dzr = t%dzr(ip)

  goff = pz * e032
  delta = ks32 * dzr
  gam = gam032 + goff

  awloc = faw32(x, y)
  px_g = px * p032
  py_g = py * p032
  btpar = 1 + px_g*px_g + py_g*py_g + aw32*aw32*awloc*awloc

  call gather32 (x, y, cpart, on_grid)
  if (on_grid) then
    rpart = rt32 * awloc * conjg(cpart)
  else
    rpart = 0
  endif

  call rk32 (goff, delta)

  ! Chart return, the reformulated forms: pz from goff through e0, the residual from
  ! delta with the beta corrections below FP32 resolution (module header).

  d_new = delta
  pz = goff / e032
  dzr = -(cret32 - d_new) / ks32

  incr(ip) = abs(real(d_new - ks32 * t%dzr(ip), rp))
  t%pz(ip) = pz
  t%dzr(ip) = dzr

  ! The FP32 side of the source phasor, at this same sequence point.

  call sincos32 (d_new, s_d, c_d)
  p32sum = p32sum + real(sl%weight(ip), sp) * awloc * base32 * cmplx(c_d, -s_d, sp)

  ! The FP64 comparison values, from the slice fel_advance just advanced.

  p_mc64 = p0_mc * (1 + sl%pz(ip))
  beta64 = p_mc64 / sqrt(p_mc64**2 + 1)
  theta64 = (phi0_new - phi0) + xks * sl%z(ip) / beta64 - xks * t%z_ref   ! delta convention.
  aw64loc = fel_fp32_faw (sl%x(ip), sl%y(ip), kx, ky, ax, ay, cos_t, sin_t)
  e_ip = cmplx(cos(phi_ref + theta64), -sin(phi_ref + theta64), rp)
  p64sum = p64sum + sl%weight(ip) * aw64loc * e_ip
  wsum = wsum + sl%weight(ip) * aw64loc

  ! Divergences: representation for the transverse set, arithmetic for the rest.

  dstat(1) = max(dstat(1), abs(real(x, rp) - x0(ip)))
  dstat(2) = max(dstat(2), abs(real(px, rp) - px0(ip)))
  dstat(3) = max(dstat(3), abs(real(y, rp) - y0(ip)))
  dstat(4) = max(dstat(4), abs(real(py, rp) - py0(ip)))
  dstat(5) = max(dstat(5), abs(real(pz, rp) - sl%pz(ip)))
  dstat(6) = max(dstat(6), abs(real(d_new, rp) - theta64))
enddo

! Scales: relative where a scale exists, absolute radians for theta.

sc(1) = maxval(abs(x0(1:n))) + 1e-30_rp
sc(2) = maxval(abs(px0(1:n))) + 1e-30_rp
sc(3) = maxval(abs(y0(1:n))) + 1e-30_rp
sc(4) = maxval(abs(py0(1:n))) + 1e-30_rp
sc(5) = maxval(abs(sl%pz(1:n))) + 1e-30_rp
sc(6) = 1

! The phasor scale is the slice's full charge weight, not |P64|: before bunching forms
! the phasor is a noise-level sum and a relative error against it means nothing. On
! this scale the row is the error in the bunching factor itself.

sc(7) = wsum + 1e-30_rp

fp32%div_slice(1:6, is) = dstat(1:6) / sc(1:6)
fp32%div_slice(7, is) = abs(cmplx(p32sum, kind=rp) - p64sum) / sc(7)

! The bunching magnitudes, refreshed every step so run end holds the exit values.

fp32%bmag64(is) = abs(p64sum) / sc(7)
fp32%bmag32(is) = abs(cmplx(p32sum, kind=rp)) / sc(7)

! The guard statistic: the median per-step phase-residual increment in ulps of the
! residual's own magnitude. The silent failure is this number reaching zero.

call fel_fp32_median (incr, n, ulp_med)
fp32%ulp_slice(is) = ulp_med / max(real(spacing(maxval(abs(ks32 * t%dzr(1:n)))), rp), tiny(1.0_rp))

!------------------------------------------------------------------------------
contains

! The FP32 image of faw (fel_track_mod), the same expressions in sp.

function faw32 (xx, yy) result (value)
real(sp) xx, yy, value, ddx, ddy, ddt
ddx = xx - ax32
ddy = yy - ay32
if (st32 /= 0) then
  ddt = ct32 * ddx + st32 * ddy
  ddy = -st32 * ddx + ct32 * ddy
  ddx = ddt
endif
value = 1 + 0.5_sp * (kx32 * ddx*ddx + ky32 * ddy*ddy)
end function faw32

! The FP32 gather: fel_grid_weights_pre's expressions in sp, the four corner values
! converted to FP32 before the blend, as a resident FP32 field would hold them.

subroutine gather32 (xx, yy, cp, ong)
real(sp) xx, yy, wwx, wwy
complex(sp) cp
logical ong
integer jx, jy
if (xx > -gmax32 .and. xx < gmax32 .and. yy > -gmax32 .and. yy < gmax32) then
  wwx = (xx + gmax32) / dxg32
  wwy = (yy + gmax32) / dyg32
  jx = int(floor(wwx));  jy = int(floor(wwy))
  wwx = 1 + real(jx, sp) - wwx
  wwy = 1 + real(jy, sp) - wwy
  jx = jx + 1;  jy = jy + 1

  ! FP32 rounding at the outer boundary can land one cell past where the FP64 form
  ! stops. That is an off-grid particle at this precision, and it counts as one.

  if (jx < 1 .or. jy < 1 .or. jx+1 > size(exfld,1) .or. jy+1 > size(exfld,2)) then
    cp = 0
    ong = .false.
    return
  endif
  cp =      fp32%e32(jx,   jy,   is) * wwx * wwy
  cp = cp + fp32%e32(jx+1, jy,   is) * (1-wwx) * wwy
  cp = cp + fp32%e32(jx,   jy+1, is) * wwx * (1-wwy)
  cp = cp + fp32%e32(jx+1, jy+1, is) * (1-wwx) * (1-wwy)
  ong = .true.
else
  cp = 0
  ong = .false.
endif
end subroutine gather32

! The FP32 image of fel_runge_kutta on (goff, delta), the same stage structure.

subroutine rk32 (g, d)
real(sp) g, d
real(sp) k2gg, k2pp, k3gg, k3pp, stpz
k2gg = 0;  k2pp = 0
call ode32 (g, d, k2gg, k2pp)
stpz = 0.5_sp * dz32
g = g + stpz * k2gg;  d = d + stpz * k2pp
k3gg = k2gg;  k3pp = k2pp
k2gg = 0;  k2pp = 0
call ode32 (g, d, k2gg, k2pp)
g = g + stpz * (k2gg - k3gg);  d = d + stpz * (k2pp - k3pp)
k3gg = k3gg / 6;  k3pp = k3pp / 6
k2gg = k2gg * (-0.5_sp);  k2pp = k2pp * (-0.5_sp)
call ode32 (g, d, k2gg, k2pp)
stpz = dz32
g = g + stpz * k2gg;  d = d + stpz * k2pp
k3gg = k3gg - k2gg;  k3pp = k3pp - k2pp
k2gg = k2gg * 2;  k2pp = k2pp * 2
call ode32 (g, d, k2gg, k2pp)
g = g + stpz * (k3gg + k2gg / 6.0_sp)
d = d + stpz * (k3pp + k2pp / 6.0_sp)
end subroutine rk32

! The FP32 image of fel_ode under the detuning reformulation (module header):
! k2pp gains 0.5*ks*(qres - q) in place of ks*(1 - 1/btpar0) + ku, whose sqrt
! argument is exactly 1 in FP32. b0 is 1 to FP32 resolution and is not formed.

subroutine ode32 (g, d, kg, kp)
real(sp) g, d, kg, kp
real(sp) s_t, c_t, gam_l, q_l, btper
complex(sp) ctmp
call sincos32 (d, s_t, c_t)
ctmp = rpart * (base32 * cmplx(c_t, -s_t, sp))
gam_l = gam032 + g
btper = btpar + (-2.0_sp / ks32) * real(ctmp, sp)
q_l = btper / (gam_l * gam_l)
kp = kp + 0.5_sp * ks32 * (qres32 - q_l)
kg = kg + aimag(ctmp) / gam_l
end subroutine ode32

! The FP32 pair, plain intrinsics on the small angle.

subroutine sincos32 (a, s, c)
real(sp) a, s, c
s = sin(a)
c = cos(a)
end subroutine sincos32

end subroutine fel_fp32_twin_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_field_twin (fp32, is, sl, beam, exfld_post, exp_k2, scl_w, &
!                                    dx, gridmax, kx2, ky2, ax, ay, cos_t, sin_t, xks)
!
! Routine to run one slice's FP32 field step and fill the source and field rows. Serial,
! from the caller's epilogue after the production solve loop: the FP64 record already
! carries this step's solve, and the twin recomputes the same step in FP32 from its own
! field record (rounded pre-solve in lockstep, carried in freerun).
!
! The source comparison uses the shared post-step FP64 particle state on both sides
! (the FP32 side through rounded values), so the rows price the field arithmetic alone:
! the deposit's weights, phases and FP32 accumulation, the single-precision transform
! pair, and the rounded propagator. The freerun deposit instead uses the twin's own
! longitudinal state, since there the FP32 field must be fed by the FP32 run. The
! deposit expressions mirror fel_field_step exactly: part = sqrt(faw2)*scl_w*w/gamma
! with faw2 = 1 + kx*dx^2 + ky*dy^2 (no half: Genesis's own roll-off, transcribed),
! cpart = (sin theta + i cos theta) * part = i e^{-i theta} * part, bilinear scatter,
! transform, kernel multiply, inverse transform, 1/N, plus 2*source in real space.
!
! Both rows are normalized by the post-solve FP64 field's norm: the source alone can
! sit at noise level on a dark start, and the post-solve field bounds it from below.
!-

subroutine fel_fp32_field_twin (fp32, is, sl, beam, exfld_post, exp_k2, scl_w, &
                                dx, gridmax, kx2, ky2, ax, ay, cos_t, sin_t, xks)

type (fel_fp32_struct), target :: fp32
type (fel_fp32_slice_struct), pointer :: t
type (fel_slice_struct) sl
type (fel_beam_struct) beam
integer is
complex(rp) exfld_post(:,:), exp_k2(:,:)
real(rp) scl_w, dx, gridmax, kx2, ky2, ax, ay, cos_t, sin_t, xks

complex(rp), allocatable :: s64(:,:)
complex(sp), allocatable :: s32(:,:)
complex(sp) cbase
real(rp) p0_mc, p_mc, gam, beta, theta, phi_dep, enorm
real(rp) w_s, w_f
real(sp) x_s, y_s, pz_s, gam_s, del_s
real(sp) scl32, gmax32, dx32, kx32, ky32, ax32, ay32, ct32, st32, gam032, p032, e032, ks32
integer ip, ix, iy, ng, n

!

t => fp32%sl32(is)
n = sl%n
ng = size(exfld_post, 1)
p0_mc = fel_p0_mc(beam)

allocate (s64(ng, ng), s32(ng, ng))
s64 = 0
s32 = 0

scl32 = real(scl_w, sp);  gmax32 = real(gridmax, sp);  dx32 = real(dx, sp)
kx32 = real(kx2, sp);  ky32 = real(ky2, sp);  ax32 = real(ax, sp);  ay32 = real(ay, sp)
ct32 = real(cos_t, sp);  st32 = real(sin_t, sp)
gam032 = real(fel_gamma0(beam), sp);  p032 = real(p0_mc, sp)
e032 = real(p0_mc**2 / fel_gamma0(beam), sp)
ks32 = real(xks, sp)

! The deposit phase base, FP64 once per slice: i e^{-i(phi0 + ks z_ref)}, so the
! per-particle FP32 angle is the small residual phase, the particle twin's convention.

phi_dep = beam%phi0 + xks * t%z_ref
cbase = cmplx(cmplx(sin(phi_dep), cos(phi_dep), rp), kind=sp)   ! i e^{-i phi} = (sin phi + i cos phi).

do ip = 1, n

  ! The FP64 side, fel_field_step's own expressions on the post-step state.

  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam
  theta = beam%phi0 + xks * sl%z(ip) / beta
  call fel_fp32_deposit64 (s64, sl%x(ip), sl%y(ip), theta, sl%weight(ip) / gam, &
                           scl_w, gridmax, dx, kx2, ky2, ax, ay, cos_t, sin_t)

  ! The FP32 side. Lockstep prices the field arithmetic from the same shared state, so
  ! it rounds the FP64 particle; freerun feeds the field from the twin's own run, so it
  ! uses the twin's longitudinal state with the transverse rounded at this sequence
  ! point (the transverse maps stay FP64, the standing model).

  x_s = real(sl%x(ip), sp)
  y_s = real(sl%y(ip), sp)
  if (fp32%freerun) then
    del_s = ks32 * t%dzr(ip)
    pz_s = t%pz(ip)
  else
    del_s = real(theta - phi_dep, sp)
    pz_s = real(sl%pz(ip), sp)
    ! The check's mutation hook reaches the field side here: in lockstep the deposit
    ! reads the rounded FP64 state, so the residual coarsening alone would not move
    ! these rows, and a falsifiable check needs it to.
    if (fp32%mutate) del_s = anint(del_s / (256 * spacing(del_s))) * (256 * spacing(del_s))
  endif
  gam_s = gam032 + pz_s * e032
  call dep32 (x_s, y_s, del_s, real(sl%weight(ip), sp) / gam_s)
enddo

! Transform, propagate, add: the FP64 record already holds the production result, and
! the twin applies the same step to its FP32 record with the single-precision transform
! and the rounded propagator.

call fp32_kernel_cache (fp32, exp_k2, ng)
call fft32_solve (fp32%e32(:,:,is), fp32%k32, ng)
fp32%e32(:,:,is) = fp32%e32(:,:,is) + 2 * s32

! The rows, both against the post-solve FP64 field's norm.

enorm = sqrt(sum(real(exfld_post, rp)**2 + aimag(exfld_post)**2)) + 1e-30_rp
w_s = 0
w_f = 0
do iy = 1, ng
  do ix = 1, ng
    w_s = w_s + abs(cmplx(s32(ix,iy), kind=rp) - s64(ix,iy))**2
    w_f = w_f + abs(cmplx(fp32%e32(ix,iy,is), kind=rp) - exfld_post(ix,iy))**2
  enddo
enddo
fp32%div_slice(8, is) = sqrt(w_s) * 2 / enorm     ! As the field sees it: the source adds times 2.
fp32%div_slice(9, is) = sqrt(w_f) / enorm

! The exit observables, kept fresh every step so the last step's values are the run's.

fp32%pow64(is) = sum(real(exfld_post, rp)**2 + aimag(exfld_post)**2)
fp32%pow32(is) = sum(real(real(fp32%e32(:,:,is), rp))**2 + real(aimag(fp32%e32(:,:,is)), rp)**2)

!------------------------------------------------------------------------------
contains

subroutine dep32 (xx, yy, del, wg)
real(sp) xx, yy, del, wg, ppart, wwx, wwy, ddx, ddy, ddt, ss_l, cc_l
complex(sp) cp_l
integer jx, jy
if (.not. (xx > -gmax32 .and. xx < gmax32 .and. yy > -gmax32 .and. yy < gmax32)) return
wwx = (xx + gmax32) / dx32
wwy = (yy + gmax32) / dx32
jx = int(floor(wwx));  jy = int(floor(wwy))
wwx = 1 + real(jx, sp) - wwx
wwy = 1 + real(jy, sp) - wwy
jx = jx + 1;  jy = jy + 1
if (jx < 1 .or. jy < 1 .or. jx+1 > ng .or. jy+1 > ng) return
ddx = xx - ax32;  ddy = yy - ay32
if (st32 /= 0) then
  ddt = ct32 * ddx + st32 * ddy
  ddy = -st32 * ddx + ct32 * ddy
  ddx = ddt
endif
ppart = sqrt(1 + kx32 * ddx*ddx + ky32 * ddy*ddy) * scl32 * wg
ss_l = sin(del);  cc_l = cos(del)

! (sin(phi+del) + i cos(phi+del)) = (sin phi + i cos phi) * (cos del - i sin del).

cp_l = cbase * cmplx(cc_l, -ss_l, sp) * ppart
s32(jx,   jy)   = s32(jx,   jy)   + (wwx * wwy) * cp_l
s32(jx+1, jy)   = s32(jx+1, jy)   + ((1-wwx) * wwy) * cp_l
s32(jx,   jy+1) = s32(jx,   jy+1) + (wwx * (1-wwy)) * cp_l
s32(jx+1, jy+1) = s32(jx+1, jy+1) + ((1-wwx) * (1-wwy)) * cp_l
end subroutine dep32

end subroutine fel_fp32_field_twin

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fp32_kernel_cache (fp32, exp_k2, ng)
!
! Routine to hold the FP32 image of the FP64 propagator, rebuilt only when the FP64
! kernel it rounds from changes (keyed by the values themselves: the caller's kernel
! cache already resolves grid, wavelength and step).
!-

subroutine fp32_kernel_cache (fp32, exp_k2, ng)

type (fel_fp32_struct) fp32
complex(rp) exp_k2(:,:)
integer ng
real(rp) key(4)

key = [real(ng, rp), real(exp_k2(1,1), rp), aimag(exp_k2(1,1)), aimag(exp_k2(ng/2+1, ng/2+1))]
if (allocated(fp32%k32)) then
  if (all(fp32%k32_key == key)) return
  deallocate (fp32%k32)
endif
allocate (fp32%k32(ng, ng))
fp32%k32 = cmplx(exp_k2, kind=sp)
fp32%k32_key = key

end subroutine fp32_kernel_cache

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fft32_solve (e, k32, ng)
!
! Routine to apply one field-solve step in single precision: forward transform,
! propagator multiply, inverse transform, 1/N. FFTW's single-precision interface,
! plans cached per grid size with FFTW_ESTIMATE (FFTW_MEASURE picks its algorithm by
! timing and is nondeterministic at the ulp level, the same decision wavefront_mod
! records). Serial by design: the caller's epilogue runs slices one at a time, so one
! plan pair and one aligned buffer suffice.
!
! A build without the single-precision library compiles the body out, since a call to
! fftwf_execute_dft would fail to link at all. Setup refuses such a build by name
! (fel_fp32_have_fftw3f$), so this routine is unreachable there, and the message says
! so rather than pretending to a fallback.
!-

subroutine fft32_solve (e, k32, ng)

use, intrinsic :: iso_c_binding

complex(sp) e(:,:), k32(:,:)
integer ng

#ifdef LUCIFER_HAVE_FFTW3F

include 'fftw3.f03'

type (c_ptr), save :: plan_f = c_null_ptr, plan_b = c_null_ptr, pbuf = c_null_ptr
integer, save :: ng_plan = 0
complex(c_float_complex), pointer, save :: buf(:,:) => null()
character(*), parameter :: r_name = 'fft32_solve'

!

if (ng /= ng_plan) then
  if (c_associated(plan_f)) then
    call fftwf_destroy_plan (plan_f)
    call fftwf_destroy_plan (plan_b)
    call fftwf_free (pbuf)
  endif
  pbuf = fftwf_alloc_complex (int(ng * ng, c_size_t))
  call c_f_pointer (pbuf, buf, [ng, ng])
  plan_f = fftwf_plan_dft_2d (ng, ng, buf, buf, FFTW_FORWARD,  FFTW_ESTIMATE)
  plan_b = fftwf_plan_dft_2d (ng, ng, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE)
  ng_plan = ng
endif

buf = e
call fftwf_execute_dft (plan_f, buf, buf)
buf = buf * k32
call fftwf_execute_dft (plan_b, buf, buf)
e = buf / real(ng * ng, sp)

#else

character(*), parameter :: r_name = 'fft32_solve'

! The stop stays, for the reason fel_assert_wiggler_sane records: err_exit's traceback
! bomb does not trap on arm64 and its bare stop exits zero, so a refusal that must be
! seen has to stop with a nonzero status itself.

call out_io (s_fatal$, r_name, 'THE FP32 FIELD TWIN NEEDS THE SINGLE-PRECISION FFTW (libfftw3f),', &
      'WHICH THIS BUILD DOES NOT CARRY. SETUP REFUSES SUCH A RUN, SO REACHING HERE', &
      'IS A BUG. PLEASE REPORT THIS!')
stop 1

#endif

end subroutine fft32_solve

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_step_close (fp32, err_flag)
!
! Routine to reduce the step's per-slice divergences, write the row, and run the
! guard. Serial, after the caller's parallel loop. The guard refuses when the median
! per-step residual increment falls below 32 ulps of the residual: below that the
! FP32 representation is absorbing the physics, which is the silent failure the
! design brief records.
!-

subroutine fel_fp32_step_close (fp32, err_flag)

type (fel_fp32_struct) fp32
logical err_flag
real(rp) row(fel_fp32_nq$), gulp
integer iq
character(*), parameter :: r_name = 'fel_fp32_step_close'

!

err_flag = .false.
fp32%istep = fp32%istep + 1

do iq = 1, fel_fp32_nq$
  row(iq) = maxval(fp32%div_slice(iq, :))
  fp32%worst(iq) = max(fp32%worst(iq), row(iq))
enddo
gulp = minval(fp32%ulp_slice)
fp32%ulp_min = min(fp32%ulp_min, gulp)

write (fp32%iu, '(i8, 10es13.4)') fp32%istep, row, gulp

if (gulp < 32) then
  call out_io (s_error$, r_name, 'FP32 RESIDUAL GUARD: THE MEDIAN PER-STEP PHASE INCREMENT IS \es10.2\ ULPS', &
        'OF THE RESIDUAL, BELOW THE FLOOR OF 32. AT THIS GRANULARITY FP32 ABSORBS THE PHYSICS', &
        'SILENTLY AND A RUN LOOKS CLEAN WITH NO BUNCHING (doc/performance.md), SO IT REFUSES INSTEAD.', &
        r_array = [gulp])
  err_flag = .true.
endif

end subroutine fel_fp32_step_close

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_close (fp32)
!
! Routine to write the run summary block and close the stream.
!-

subroutine fel_fp32_close (fp32)

type (fel_fp32_struct) fp32
real(rp) pw
integer iq

!

if (.not. fp32%on) return
if (fp32%iu == 0) return

write (fp32%iu, '(a)') '#'
do iq = 1, fel_fp32_nq$
  write (fp32%iu, '(2a, es13.4)') 'worst_', trim(fel_fp32_qname(iq)), fp32%worst(iq)
enddo
write (fp32%iu, '(a, es13.4)') 'guard_ulp_min ', fp32%ulp_min

! The end-to-end exit observables: the last step's field power and bunching, FP32
! against FP64. In freerun this is the whole-run compounding a device port will be
! judged against; in lockstep it restates the last per-step rows and says so.

if (allocated(fp32%pow64)) then
  pw = 0
  do iq = 1, size(fp32%pow64)
    if (fp32%pow64(iq) > 0) pw = max(pw, abs(fp32%pow32(iq) - fp32%pow64(iq)) / fp32%pow64(iq))
  enddo
  write (fp32%iu, '(a, es13.4)') 'endtoend_power_rel ', pw
  write (fp32%iu, '(a, es13.4)') 'endtoend_power_total_rel ', &
        abs(sum(fp32%pow32) - sum(fp32%pow64)) / max(sum(fp32%pow64), tiny(1.0_rp))
  write (fp32%iu, '(a, es13.4)') 'endtoend_bunching_abs ', maxval(abs(fp32%bmag32 - fp32%bmag64))
endif
write (fp32%iu, '(a, i0)') 'steps ', fp32%istep
close (fp32%iu)
fp32%iu = 0

end subroutine fel_fp32_close

end module fel_fp32_mod
