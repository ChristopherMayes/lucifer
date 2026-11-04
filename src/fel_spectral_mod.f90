!+
! Module fel_spectral_mod
!
! The spectral field propagator, field.propagator = "spectral" (docs/spectral-field.md).
! The slice propagator (fel_field_step, fel_field_diffract) advances every slice of the
! window by the paraxial kernel exp(-i dz k_perp^2 / 2 k0) at the one reference
! wavenumber k0, and the slippage moves the record by whole slices once 0.8 of a slice
! has accumulated (fel_apply_slippage). This module transforms a member's whole window
! in (x, y, tau) instead. Two things change and the rest does not. The kernel becomes
! exp(-i dz k_perp^2 / 2 k(kappa)) with k(kappa) = k0 + kappa, where kappa is the
! wavenumber conjugate to the window coordinate tau and c kappa is the frequency offset
! from the carrier. That is the term a chirped or short pulse gets wrong at relative
! order dOmega/omega under one k0. And the slippage becomes the phase ramp
! exp(-i kappa delta), delta the step's slippage as a length, exact and continuous, so
! a step may slip a fraction of a slice and the field is never a fraction of a slice
! behind the beam. The push, the deposit and the averaged equations are untouched.
!
! The window coordinate. Window slice i (fel-physics.md sec-window) sits at
! tau = (i - 1) * slice_spacing, so tau grows toward the head, the end light slips
! toward. The envelope rides the carrier exp(i (k0 z - omega0 t)), the convention the
! kernel's sign fixes, so a mode exp(i kappa tau) of the envelope has the carrier's
! frequency plus c kappa and the wavenumber k0 + kappa. The discrete transform of the
! window carries kappa_m = 2 pi m / (n_t slice_spacing) for m up to half of n_t and
! negative beyond, as the transverse table does (fel_field_kernel_init).
!
! The head and the bank. The transform along tau is periodic, so light the ramp moves
! past the head would reappear at the tail. The window is followed along tau by a pad of
! fel_spectral_n_free$ free slots and a cos^2 sponge of fel_spectral_n_sponge$ slots
! that falls to zero at the pad's end, and every transform runs over the window and the
! pad together. Light the ramp carries past the head enters the pad, diffracts there
! under the kernel as escaped light does in vacuum, and is attenuated by the sponge
! after every slip pass, so the field is smooth at the periodic seam and nothing reaches
! the tail. The escape bank takes its planes from the pad's first slot. The slippage
! accumulates in slip%ramp_slip, and each time it crosses a whole slice the pass is
! split there in the transform domain, the field at that slot is read by a direct sum
! over the modes, and the plane goes to the bank as fel_apply_slippage's transmitted
! slice does, its energy added to u_escaped. Between two crossings the ramp moves the
! field by exactly one slot, which the transform does exactly, so the bank receives the
! window's samples one slice apart, each once, and the pulse file is assembled as
! before. The sponge removes light the bank already holds and counts nothing. One pass
! slips at most one slice, so light walks through every slot of the sponge and reaches
! its end at the product of the weights, below 1e-5 of its amplitude two slots before
! the seam, and a larger slippage, an interlude's, is taken in such passes. A pass of
! several slices was measured to jump light into the sponge's last slots, where the
! next pass's interpolation carried 2e-3 of the exit peak power round to the tail slice
! of the SASE tier's window.
!
! Note: A dump carries the window and not the pad. A continuation starts with an empty
! pad, and every plane that had entered it was banked when it crossed, so the ledger
! loses nothing.
!
! The kernel factors. exp(-i dz (m_x^2 + m_y^2) dk^2 / 2 k(kappa)) is the product of a
! table in (m_x, kappa) and one in (m_y, kappa), each n_grid by n_t, so the whole kernel
! costs one complex multiply a point and the tables are rebuilt at every call.
!
! The transverse transforms are wavefront_fft2's, per slice and parallel over slices.
! The transform along tau runs on each line of the transverse spectrum, parallel over
! the second index, through a plan cache of this module's own in wavefront_fft2's form:
! one plan pair and one buffer per thread, planned serially before the parallel loop
! and executed inside it.
!-

module fel_spectral_mod

use fel_track_mod
use fel_device_mod
use, intrinsic :: iso_c_binding

implicit none

private spectral_fft1, spectral_fft1_plan, spectral_fft1_plan_threads, bank_append

! The pad beyond the window's head (module header), in slices.

integer, parameter :: fel_spectral_n_free$ = 4      ! Free slots, the first where the bank reads its planes.
integer, parameter :: fel_spectral_n_sponge$ = 12   ! The cos^2 sponge after them.

! The field's sampling along tau under the spectral propagator, field.n_wavelength
! (fel_field_sample, docs/spectral-field.md): the fraction of the cooperation length
! lambda/(4 pi rho) a derived plane spans, the fraction above which a stated one is
! refused, and the planes a window keeps at least. The two fractions are the sweep's.

real(rp), parameter :: fel_field_coop_fraction$ = 0.1_rp
real(rp), parameter :: fel_field_coop_cliff$ = 0.25_rp
integer, parameter :: fel_field_planes_min$ = 8

!+
! Struct fel_spectral_struct
!
! The spectral propagator's switch and the pad it sized for the window.
!-

type fel_spectral_struct
  logical :: on = .false.      ! field.propagator = "spectral".
  integer :: n_free = 0        ! Free slots of the pad beyond the head.
  integer :: n_sponge = 0      ! Sponge slots after them.
  integer :: n_pad = 0         ! Their sum. 0 on a window of one slice, which never slips.
end type

! The plan cache for the transform along tau, one per thread (see the module header).

type(C_PTR), private, save :: sp_plan_fwd = C_NULL_PTR
type(C_PTR), private, save :: sp_plan_bwd = C_NULL_PTR
type(C_PTR), private, save :: sp_buf_cptr = C_NULL_PTR
complex(wf_rp), private, pointer, save :: sp_buf(:) => null()
integer, private, save :: sp_cache_n = 0
!$OMP threadprivate(sp_plan_fwd, sp_plan_bwd, sp_buf_cptr, sp_buf, sp_cache_n)

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_setup (sp, propagator)
!
! Routine to resolve field.propagator into the switch. The key table has already
! refused a name outside the key's set. fel_spectral_size then sizes the pad.
!
! Input:
!   propagator -- character(*): field.propagator, "slices" or "spectral".
!
! Output:
!   sp         -- fel_spectral_struct: The switch set.
!-

subroutine fel_spectral_setup (sp, propagator)

type (fel_spectral_struct) sp
character(*) propagator

!

sp%on = (propagator == 'spectral')
sp%n_free = 0
sp%n_sponge = 0
sp%n_pad = 0

end subroutine fel_spectral_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_size (sp, nslice)
!
! Routine to size the pad (module header). A window of one slice, the steady state,
! takes none, since nothing slips.
!
! Input:
!   nslice -- integer: The window's slices.
!
! Output:
!   sp     -- fel_spectral_struct: n_free, n_sponge and n_pad set.
!-

subroutine fel_spectral_size (sp, nslice)

type (fel_spectral_struct) sp
integer nslice

!

if (nslice == 1) then
  sp%n_free = 0
  sp%n_sponge = 0
else
  sp%n_free = fel_spectral_n_free$
  sp%n_sponge = fel_spectral_n_sponge$
endif
sp%n_pad = sp%n_free + sp%n_sponge

end subroutine fel_spectral_size

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_step (sp, und, beam, ff, dz, slippage, keep, err_flag)
!
! Routine to do the field half of one undulator step under the spectral propagator,
! in the slice propagator's order (fel_field_step, then fel_apply_slippage): each
! member's window is propagated over the step by the kernel, the source is deposited
! slice by slice from the just-advanced particles (fel_spectral_deposit), and the
! slippage then moves the field and the fresh source together as the ramp, with the
! planes that cross a slice banked (fel_spectral_slip). fel_track_und_step has pushed
! the particles and skipped its own field solve. The caller drains each member's bank
! afterwards (fel_drain_bank). The two halves are fel_spectral_step_field and
! fel_spectral_step_slip, which the walk calls apart when the instrument's rows are
! read between them (fel_track_und_twin_close).
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   und      -- fel_und_struct: The segment.
!   beam     -- fel_beam_struct: The beam, pushed over the step.
!   ff(:)    -- fel_field_struct: The field set.
!   dz       -- real(rp): Step length [m].
!   slippage -- real(rp): The step's slippage [fundamental wavelengths].
!   keep     -- logical: Bank the planes that cross? False adds their energy alone.
!
! Output:
!   ff(:)    -- fel_field_struct: Every member advanced, its bank holding the planes.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_step (sp, und, beam, ff, dz, slippage, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
real(rp) dz, slippage
logical keep, err_flag

!

call fel_spectral_step_field (sp, und, beam, ff, dz, keep, err_flag)
if (err_flag) return
call fel_spectral_step_slip (sp, ff, slippage, keep, err_flag)

end subroutine fel_spectral_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_step_field (sp, und, beam, ff, dz, keep, err_flag)
!
! Routine to do the first half of fel_spectral_step, each member's kernel pass over
! the step and its deposit from the just-advanced particles, without the slippage.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   und      -- fel_und_struct: The segment.
!   beam     -- fel_beam_struct: The beam, pushed over the step.
!   ff(:)    -- fel_field_struct: The field set.
!   dz       -- real(rp): Step length [m].
!   keep     -- logical: Bank the planes that cross? Passed to the kernel pass, which
!                 crosses nothing.
!
! Output:
!   ff(:)    -- fel_field_struct: Every member propagated with its source added.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_step_field (sp, und, beam, ff, dz, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
real(rp) dz
logical keep, err_flag

real(rp) ks
integer io, is, nslice_f
logical any_err, err

!

err_flag = .true.
ks = twopi / ff(1)%wf%wavelength

do io = 1, size(ff)
  call fel_spectral_apply (sp, ff(io), dz, 0.0_rp, .true., keep, err)
  if (err) return
  nslice_f = size(ff(io)%wf%Ex, 3)
  any_err = .false.
  if (ff(io)%n_per > 1) then
    call fel_spectral_deposit_coarse (und, beam, ff(io), dz, ks, err)
    if (err) return
  else
    !$OMP parallel do private(err) reduction(.or.: any_err)
    do is = 1, size(beam%slice)
      call fel_spectral_deposit (und, beam, beam%slice(is), ff(io)%wf, &
                                 fel_field_index(ff(io)%slip, is, nslice_f), dz, ff(io)%harm, ks, err)
      any_err = any_err .or. err
    enddo
    !$OMP end parallel do
    if (any_err) return
  endif
enddo

err_flag = .false.

end subroutine fel_spectral_step_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_step_slip (sp, ff, slippage, keep, err_flag)
!
! Routine to do the second half of fel_spectral_step, each member's slippage as the
! ramp with the planes that cross a slice banked (fel_spectral_slip).
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   ff(:)    -- fel_field_struct: The field set.
!   slippage -- real(rp): The step's slippage [fundamental wavelengths].
!   keep     -- logical: Bank the planes that cross? False adds their energy alone.
!
! Output:
!   ff(:)    -- fel_field_struct: Every member slipped, its bank holding the planes.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_step_slip (sp, ff, slippage, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: ff(:)
real(rp) slippage
logical keep, err_flag
integer io

!

err_flag = .true.
do io = 1, size(ff)
  call fel_spectral_slip (sp, ff(io), slippage, keep, err_flag)
  if (err_flag) return
enddo
err_flag = .false.

end subroutine fel_spectral_step_slip

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_drift (sp, ff, dz, slippage, keep, err_flag)
!
! Routine to propagate every member of the field set through an interlude under the
! spectral propagator: the kernel over the interlude's length with no source, then
! the slippage the schedule charges it, as the ramp with the bank (fel_spectral_slip).
! This is what wavefront_drift and fel_apply_slippage do for the slice propagator, with
! the kernel at k(kappa) and the slippage continuous. The caller drains each member's
! bank afterwards (fel_drain_bank).
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   ff(:)    -- fel_field_struct: The field set.
!   dz       -- real(rp): The light's path [m].
!   slippage -- real(rp): The interlude's slippage [fundamental wavelengths].
!   keep     -- logical: Bank the planes that cross? False adds their energy alone.
!
! Output:
!   ff(:)    -- fel_field_struct: Every member advanced, its bank holding the planes.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_drift (sp, ff, dz, slippage, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: ff(:)
real(rp) dz, slippage
logical keep, err_flag
integer io

!

err_flag = .true.
do io = 1, size(ff)
  call fel_spectral_apply (sp, ff(io), dz, 0.0_rp, .true., keep, err_flag)
  if (err_flag) return
  call fel_spectral_slip (sp, ff(io), slippage, keep, err_flag)
  if (err_flag) return
enddo
err_flag = .false.

end subroutine fel_spectral_drift

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_slip (sp, fs, slippage, keep, err_flag)
!
! Routine to slip one member by the ramp, in passes along tau alone of at most one
! slice each (module header), the plane that crosses a slice banked by the pass that
! crosses it. The member's bank is emptied first and holds every plane of the call on
! return, in transmission order, for the caller to drain. A steady-state record never
! slips.
!
! Note: Backward slippage is refused, as fel_apply_slippage refuses it with the bank:
! a plane that leaves at the tail has no place ahead of the window in the pulse file,
! and the pad is at the head.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   fs       -- fel_field_struct: The member.
!   slippage -- real(rp): The slippage to apply [fundamental wavelengths].
!   keep     -- logical: Bank the planes that cross? False adds their energy alone.
!
! Output:
!   fs       -- fel_field_struct: The window slipped, the bank filled, u_escaped and
!                 slip%ramp_slip advanced.
!   err_flag -- logical: Set True on a transform failure or backward slippage.
!-

subroutine fel_spectral_slip (sp, fs, slippage, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
real(rp) slippage, chunk, left
logical keep, err_flag
integer n_chunk, ic
character(*), parameter :: r_name = 'fel_spectral_slip'

!

err_flag = .false.
fs%bank%n = 0
if (.not. fs%slip%timerun) return
if (slippage == 0) return

if (slippage < 0) then
  call out_io (s_error$, r_name, 'BACKWARD SLIPPAGE UNDER FIELD.PROPAGATOR = "spectral".', &
                                 'THE PAD THAT FEEDS THE ESCAPED-FIELD BANK IS AT THE WINDOW''S HEAD,', &
                                 'AND A PLANE THAT LEAVES AT THE TAIL HAS NO PLACE IN THE PULSE FILE.', &
                                 'POSSIBLE SOLUTION: FIELD.PROPAGATOR = "slices" WITH KEEP_ESCAPED_FIELD = F.')
  err_flag = .true.
  return
endif

n_chunk = max(1, ceiling(slippage / fs%slip%sample - 1e-9_rp))
left = slippage
do ic = 1, n_chunk
  chunk = slippage / n_chunk
  if (ic == n_chunk) chunk = left
  left = left - chunk
  call fel_spectral_apply (sp, fs, 0.0_rp, chunk, .false., keep, err_flag)
  if (err_flag) return
enddo

end subroutine fel_spectral_slip

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_apply (sp, fs, dz, slippage, transverse, keep, err_flag)
!
! Routine to apply the kernel over dz to one member's whole window and pad, when
! transverse is True, or the ramp for the slippage alone when it is False, which is the
! slip pass (fel_spectral_slip). The record is read and written in window order through
! fel_field_index, so a rotated record is handled as the slice propagator handles it,
! and the pad follows it.
!
! The passes, per polarization: the forward transverse transform of every plane; on
! every line of the transverse spectrum the forward transform along tau, the kernel
! factors, and the inverse; the inverse transverse transform of every plane. The
! transforms are unnormalized, so the result is scaled by n_t and by n_grid^2. The
! ramp alone is diagonal in tau's conjugate and needs no transverse pass, so it runs
! on the lines of the field itself. A slip pass is split at each whole slice the
! accumulated slippage crosses: the ramp to the crossing, the field at the pad's first
! slot read by a direct sum over the modes into a plane for the bank, then the ramp on,
! and the sponge on the pad afterwards (module header).
!
! Input:
!   sp         -- fel_spectral_struct: The propagator's switch and pad.
!   fs         -- fel_field_struct: The member: its wavefront, pad, slippage state and bank.
!   dz         -- real(rp): Step length [m].
!   slippage   -- real(rp): The slippage to apply [fundamental wavelengths], 0 with
!                   transverse True.
!   transverse -- logical: Apply the kernel, with the transverse passes. False applies
!                   the ramp alone.
!   keep       -- logical: Bank the planes that cross? False adds their energy alone.
!
! Output:
!   fs         -- fel_field_struct: The window and pad advanced, the bank grown by the
!                   planes that crossed, u_escaped and slip%ramp_slip advanced.
!   err_flag   -- logical: Set True on a transform failure. False otherwise.
!-

subroutine fel_spectral_apply (sp, fs, dz, slippage, transverse, keep, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
type (wavefront_struct), pointer :: wf
real(rp) dz, slippage
logical transverse, keep, err_flag

complex(rp), allocatable :: ax(:,:), ay(:,:), ramp(:,:), samp(:)
complex(wf_rp), allocatable :: planes(:,:,:), planes_y(:,:,:)
real(rp), allocatable :: shift(:)
real(rp) k0, dk, dtau, delta, kappa, kt, e_scale, t_slice, w
integer, allocatable :: idx(:)
integer ngrid_arr(3), ngrid, nz, nt, ix, m, mx, mt, j, n_cross, k
logical err
character(*), parameter :: r_name = 'fel_spectral_apply'

!

err_flag = .true.
wf => fs%wf
ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
nz = ngrid_arr(3)
nt = nz + sp%n_pad
k0 = twopi / wf%wavelength
dk = twopi / (ngrid * wf%dx)
dtau = wf%dz

if (transverse .and. slippage /= 0) then
  call out_io (s_error$, r_name, 'A KERNEL PASS WAS ASKED TO SLIP. PLEASE REPORT THIS!')
  return
endif
if (sp%n_pad > 0 .and. .not. allocated(fs%pad)) then
  call out_io (s_error$, r_name, 'THE PAD WAS NOT ALLOCATED. PLEASE REPORT THIS!')
  return
endif

! The crossings of the slip pass (module header): the shift to each, then the rest.
! The slippage and the slice are in fundamental wavelengths, and a shift becomes a
! length through the member's own wavelength times its harmonic, the fundamental's.
! A steady-state record never slips.

n_cross = 0
if (fs%slip%timerun .and. slippage /= 0) then
  n_cross = floor((fs%slip%ramp_slip + slippage) / fs%slip%sample + 1e-9_rp)
endif
allocate (shift(n_cross + 1))
if (n_cross > 0) then
  shift(1) = fs%slip%sample - fs%slip%ramp_slip
  shift(2:n_cross) = fs%slip%sample
  shift(n_cross + 1) = slippage - (n_cross * fs%slip%sample - fs%slip%ramp_slip)
else
  shift(1) = slippage
endif
if (fs%slip%timerun .and. slippage /= 0) then
  fs%slip%ramp_slip = max(0.0_rp, fs%slip%ramp_slip + slippage - n_cross * fs%slip%sample)
endif

! The energy of one plane, as fel_apply_slippage takes it: the sum of |E|^2 over the grid
! times dx^2 / 2 Z0 is the slice's power, and times the slice's light time its energy. A
! harmonic member's wavelength is the fundamental's over harm, and the window is the
! fundamental's, so the light time is scaled back up.

e_scale = wf%dx**2 / (2 * (mu_0_vac * c_light))
t_slice = fs%slip%sample * wf%wavelength * fs%harm / c_light

! The two kernel tables (module header), each column one kappa and each row one
! transverse mode, the ramp of each segment of the slip pass, and the phasor that reads
! the pad's first slot, window index nz, from the unnormalized spectrum.

allocate (ax(ngrid, nt), ay(ngrid, nt), ramp(nt, n_cross + 1), samp(nt), idx(nz))
do m = 1, nt
  mt = m - 1;  if (2*mt > nt-1) mt = mt - nt
  kappa = twopi * mt / (nt * dtau)
  kt = k0 + kappa
  do ix = 1, ngrid
    mx = ix - 1;  if (2*mx > ngrid-1) mx = mx - ngrid
    ax(ix, m) = exp(cmplx(0.0_rp, -dz * (mx * dk)**2 / (2 * kt), rp))
    ay(ix, m) = ax(ix, m)
  enddo
  do k = 1, n_cross + 1
    delta = shift(k) * wf%wavelength * fs%harm
    ramp(m, k) = exp(cmplx(0.0_rp, -kappa * delta, rp))
  enddo
  samp(m) = exp(cmplx(0.0_rp, twopi * real(mt, rp) * nz / nt, rp)) / nt
enddo
do j = 1, nz
  idx(j) = fel_field_index(fs%slip, j, nz)
enddo

! The plans every thread will execute, made serially (wavefront_fft2_plan_threads).

if (transverse) then
  call wavefront_fft2_plan_threads (ngrid, ngrid, err);  if (err) return
endif
call spectral_fft1_plan_threads (nt, err);  if (err) return

allocate (planes(ngrid, ngrid, n_cross))
call apply_component (wf%Ex, fs%pad, planes, err);  if (err) return
if (allocated(wf%Ey)) then
  allocate (planes_y(ngrid, ngrid, n_cross))
  call apply_component (wf%Ey, fs%pad_y, planes_y, err);  if (err) return
endif

! The sponge, after a slip pass: the weight falls from just under 1 on the sponge's
! first slot to 0 on the pad's last, so what the next pass carries round to the tail is
! zero. The light it removes was banked when it crossed, so nothing is counted.

if (.not. transverse .and. sp%n_sponge > 0) then
  do j = 1, sp%n_sponge
    w = cos(0.5_rp * pi * real(j, rp) / sp%n_sponge)**2
    fs%pad(:,:,sp%n_free + j) = fs%pad(:,:,sp%n_free + j) * w
    if (allocated(fs%pad_y)) fs%pad_y(:,:,sp%n_free + j) = fs%pad_y(:,:,sp%n_free + j) * w
  enddo
endif

! The planes that crossed, in order: their energy to the ledger and, when kept, the
! planes themselves to the bank.

do k = 1, n_cross
  fs%slip%u_escaped = fs%slip%u_escaped + &
                sum(real(planes(:,:,k), rp)**2 + aimag(planes(:,:,k))**2) * e_scale * t_slice
  if (allocated(planes_y)) then
    fs%slip%u_escaped = fs%slip%u_escaped + &
                sum(real(planes_y(:,:,k), rp)**2 + aimag(planes_y(:,:,k))**2) * e_scale * t_slice
  endif
  if (keep) then
    if (allocated(planes_y)) then
      call bank_append (fs%bank, planes(:,:,k), planes_y(:,:,k))
    else
      call bank_append (fs%bank, planes(:,:,k))
    endif
  endif
enddo

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine apply_component (e, pad, pl, err)
!
! Routine to put one polarization's planes, the window's and the pad's, through the
! passes of the header, reading the plane at each crossing into pl.
!-

subroutine apply_component (e, pad, pl, err)

complex(wf_rp) e(:,:,:), pl(:,:,:)
complex(wf_rp), allocatable :: pad(:,:,:)
logical err
integer is_l, ix_l, iy_l, m_l, j_l, k_l
logical err_l, any_err

!

err = .true.
any_err = .false.
if (transverse) then
  !$OMP parallel do private(err_l) reduction(.or.: any_err)
  do is_l = 1, nt
    if (is_l <= nz) then
      call wavefront_fft2 (e(:,:,is_l), wf_fft_forward$, err_l)
    else
      call wavefront_fft2 (pad(:,:,is_l-nz), wf_fft_forward$, err_l)
    endif
    any_err = any_err .or. err_l
  enddo
  !$OMP end parallel do
  if (any_err) return
endif

!$OMP parallel do private(ix_l, j_l, m_l, k_l, err_l) reduction(.or.: any_err)
do iy_l = 1, ngrid
  do ix_l = 1, ngrid
    do j_l = 1, nz
      sp_buf(j_l) = e(ix_l, iy_l, idx(j_l))
    enddo
    do j_l = nz + 1, nt
      sp_buf(j_l) = pad(ix_l, iy_l, j_l - nz)
    enddo
    call spectral_fft1 (wf_fft_forward$, err_l)
    any_err = any_err .or. err_l
    if (transverse) then
      do m_l = 1, nt
        sp_buf(m_l) = sp_buf(m_l) * (ax(ix_l, m_l) * ay(iy_l, m_l))
      enddo
    endif
    do k_l = 1, n_cross
      do m_l = 1, nt
        sp_buf(m_l) = sp_buf(m_l) * ramp(m_l, k_l)
      enddo
      pl(ix_l, iy_l, k_l) = sum(sp_buf(1:nt) * samp(1:nt))
    enddo
    if (shift(n_cross + 1) /= 0) then
      do m_l = 1, nt
        sp_buf(m_l) = sp_buf(m_l) * ramp(m_l, n_cross + 1)
      enddo
    endif
    call spectral_fft1 (wf_fft_backward$, err_l)
    any_err = any_err .or. err_l
    do j_l = 1, nz
      e(ix_l, iy_l, idx(j_l)) = sp_buf(j_l) / real(nt, rp)
    enddo
    do j_l = nz + 1, nt
      pad(ix_l, iy_l, j_l - nz) = sp_buf(j_l) / real(nt, rp)
    enddo
  enddo
enddo
!$OMP end parallel do
if (any_err) return

if (transverse) then
  !$OMP parallel do private(err_l) reduction(.or.: any_err)
  do is_l = 1, nt
    if (is_l <= nz) then
      call wavefront_fft2 (e(:,:,is_l), wf_fft_backward$, err_l)
      e(:,:,is_l) = e(:,:,is_l) / real(ngrid*ngrid, rp)
    else
      call wavefront_fft2 (pad(:,:,is_l-nz), wf_fft_backward$, err_l)
      pad(:,:,is_l-nz) = pad(:,:,is_l-nz) / real(ngrid*ngrid, rp)
    endif
    any_err = any_err .or. err_l
  enddo
  !$OMP end parallel do
  if (any_err) return
endif

err = .false.

end subroutine apply_component

end subroutine fel_spectral_apply

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_device_slip (sp, dev, ff, slippage, keep, twin, err_flag)
!
! Routine to slip every member of the resident field set by the ramp on the device, as
! fel_spectral_slip slips it on the CPU: passes along tau of at most one slice each, the
! plane that crosses a slice read at the pad's first slot into a bank slot by the pass
! that crosses it (lucifer_device.h, luc_dev_tau_slip). The crossings are
! fel_spectral_apply's arithmetic on the fundamental's slippage state, which the host
! keeps in FP64 for the device as for the CPU, and every member slips together, as the
! members share the window. The kernel pass and the deposit ran in the device's own
! step (fel_track_und_step).
!
! When the device is the run (twin False) the crossed planes are downloaded after the
! passes, in transmission order, their energy added to u_escaped as fel_spectral_apply
! adds it, and, when kept, the planes appended to each member's bank for the caller to
! drain (fel_drain_bank). In the twin's role (twin True) the passes are encoded alone,
! so the twin's record slips as the CPU's does, and the state is left for the CPU's own
! step to advance.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switch and pad.
!   dev      -- fel_device_struct: The device, the field set resident on it.
!   ff(:)    -- fel_field_struct: The field set.
!   slippage -- real(rp): The slippage to apply [fundamental wavelengths].
!   keep     -- logical: Bank the planes that cross? False adds their energy alone.
!   twin     -- logical: The device is the instrument's twin, and the state stays.
!
! Output:
!   ff(:)    -- fel_field_struct: Each member's bank holding its planes, u_escaped and
!                 slip%ramp_slip advanced, unless twin.
!   err_flag -- logical: Set True on a refusal or backward slippage. False otherwise.
!-

subroutine fel_spectral_device_slip (sp, dev, ff, slippage, keep, twin, err_flag)

type (fel_spectral_struct) sp
type (fel_device_struct) dev
type (fel_field_struct), target :: ff(:)
real(rp) slippage
logical keep, twin, err_flag

type (fel_device_slip_par_struct) spar
complex(wf_rp), allocatable :: plane(:,:), plane_y(:,:)
real(rp) chunk, left, ramp_slip, sample, shift(2)
integer n_chunk, ic, n_cross, io, ngrid_arr(3), slot
character(*), parameter :: r_name = 'fel_spectral_device_slip'

!

err_flag = .false.
do io = 1, size(ff)
  ff(io)%bank%n = 0
enddo
if (.not. ff(1)%slip%timerun) return
if (slippage == 0) return

if (slippage < 0) then
  call out_io (s_error$, r_name, 'BACKWARD SLIPPAGE UNDER FIELD.PROPAGATOR = "spectral".', &
                                 'THE PAD THAT FEEDS THE ESCAPED-FIELD BANK IS AT THE WINDOW''S HEAD,', &
                                 'AND A PLANE THAT LEAVES AT THE TAIL HAS NO PLACE IN THE PULSE FILE.', &
                                 'POSSIBLE SOLUTION: FIELD.PROPAGATOR = "slices" WITH KEEP_ESCAPED_FIELD = F.')
  err_flag = .true.
  return
endif

ngrid_arr = wavefront_shape(ff(1)%wf)
allocate (plane(ngrid_arr(1), ngrid_arr(2)))
if (allocated(ff(1)%wf%Ey)) allocate (plane_y(ngrid_arr(1), ngrid_arr(2)))

! The passes, fel_spectral_slip's chunks of at most one slice, and in each the
! crossing and the shifts of fel_spectral_apply, at most one crossing a pass. A shift in
! fundamental wavelengths becomes a shift in planes over the plane's wavelengths.

sample = ff(1)%slip%sample
ramp_slip = ff(1)%slip%ramp_slip
n_chunk = max(1, ceiling(slippage / sample - 1e-9_rp))
left = slippage
slot = 0
do ic = 1, n_chunk
  chunk = slippage / n_chunk
  if (ic == n_chunk) chunk = left
  left = left - chunk

  n_cross = floor((ramp_slip + chunk) / sample + 1e-9_rp)
  if (n_cross > 0) then
    shift(1) = sample - ramp_slip
    shift(2) = chunk - (n_cross * sample - ramp_slip)
  else
    shift(1) = chunk
    shift(2) = 0
  endif
  ramp_slip = max(0.0_rp, ramp_slip + chunk - n_cross * sample)

  spar%s_a = shift(1) / ff(1)%n_wavelength
  spar%s_b = shift(2) / ff(1)%n_wavelength
  spar%n_cross = n_cross
  spar%slot = slot
  call fel_device_slip_pass (dev, spar, err_flag)
  if (err_flag) return

  if (n_cross > 0) slot = slot + 1
  if (slot > 0 .and. (slot == fel_dev_bank_slots$ .or. ic == n_chunk)) then
    if (.not. twin) call drain_slots (slot)
    slot = 0
  endif
enddo

if (.not. twin) then
  do io = 1, size(ff)
    ff(io)%slip%ramp_slip = ramp_slip
  enddo
endif

!------------------------------------------------------------------------------
contains

!+
! Subroutine drain_slots (n_slot)
!
! Routine to download the crossed planes in the first n_slot bank slots of every
! member, in transmission order, add each one's energy to the member's ledger as
! fel_spectral_apply adds it, and append it to the member's bank when kept.
!-

subroutine drain_slots (n_slot)

integer n_slot, k, io_l
real(rp) e_scale, t_slice

!

do io_l = 1, size(ff)
  e_scale = ff(io_l)%wf%dx**2 / (2 * (mu_0_vac * c_light))
  t_slice = ff(io_l)%slip%sample * ff(io_l)%wf%wavelength * ff(io_l)%harm / c_light
  do k = 0, n_slot - 1
    call fel_device_bank_plane (dev, io_l, 1, k, plane)
    ff(io_l)%slip%u_escaped = ff(io_l)%slip%u_escaped + &
                  sum(real(plane, rp)**2 + aimag(plane)**2) * e_scale * t_slice
    if (allocated(plane_y)) then
      call fel_device_bank_plane (dev, io_l, 2, k, plane_y)
      ff(io_l)%slip%u_escaped = ff(io_l)%slip%u_escaped + &
                  sum(real(plane_y, rp)**2 + aimag(plane_y)**2) * e_scale * t_slice
    endif
    if (keep) then
      if (allocated(plane_y)) then
        call bank_append (ff(io_l)%bank, plane, plane_y)
      else
        call bank_append (ff(io_l)%bank, plane)
      endif
    endif
  enddo
enddo

end subroutine drain_slots

end subroutine fel_spectral_device_slip

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine bank_append (bank, plane, plane_y)
!
! Routine to append one transmitted plane to a member's escape bank, growing it as
! fel_apply_slippage grows it: four planes at first, then doubling.
!-

subroutine bank_append (bank, plane, plane_y)

type (fel_bank_struct) bank
complex(wf_rp) plane(:,:)
complex(wf_rp), optional :: plane_y(:,:)
complex(wf_rp), allocatable :: grow(:,:,:)

!

if (.not. allocated(bank%plane)) then
  allocate (bank%plane(size(plane,1), size(plane,2), 4))
elseif (bank%n == size(bank%plane, 3)) then
  call move_alloc (bank%plane, grow)
  allocate (bank%plane(size(grow,1), size(grow,2), 2*size(grow,3)))
  bank%plane(:,:,1:size(grow,3)) = grow
  deallocate (grow)
endif
if (present(plane_y) .and. .not. allocated(bank%plane_y)) then
  allocate (bank%plane_y(size(plane,1), size(plane,2), size(bank%plane,3)))
elseif (allocated(bank%plane_y) .and. size(bank%plane_y,3) < size(bank%plane,3)) then
  call move_alloc (bank%plane_y, grow)
  allocate (bank%plane_y(size(grow,1), size(grow,2), size(bank%plane,3)))
  bank%plane_y(:,:,1:size(grow,3)) = grow
  deallocate (grow)
endif
bank%n = bank%n + 1
bank%plane(:,:,bank%n) = plane
if (present(plane_y)) bank%plane_y(:,:,bank%n) = plane_y

end subroutine bank_append

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_deposit (und, beam, sl, wf, ifld, delz, harm, xks1, err_flag)
!
! Routine to deposit one slice's source into its field record with no propagation.
! This is fel_field_step's deposit, the per-particle scatter with the source filter when
! it is on, kept apart from that routine so the slice propagator stays as it was. The
! source scale, the phase, the bilinear weights and the filter's transform pair are
! fel_field_step's (fel-physics.md sec-field, sec-source-filter), and so is the add of
! 2*crsource, through the element's polarization pair when both planes are live. The
! coherent source model is not carried here, and setup refuses it under the spectral
! propagator.
!
! Input:
!   und      -- fel_und_struct: Undulator parameters.
!   beam     -- fel_beam_struct: The beam.
!   sl       -- fel_slice_struct: One slice's packed particles.
!   wf       -- wavefront_struct: One field of the set.
!   ifld     -- integer: Field-record index of the slice.
!   delz     -- real(rp): Step length [m].
!   harm     -- integer: This field's harmonic number.
!   xks1     -- real(rp): Fundamental radiation wavenumber [1/m].
!
! Output:
!   wf       -- wavefront_struct: Source added to the slice's field record.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_deposit (und, beam, sl, wf, ifld, delz, harm, xks1, err_flag)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (wavefront_struct), target :: wf
integer ifld, harm
real(rp) delz, xks1
logical err_flag

real(rp) dgrid, scl_w, part, theta, beta, wx, wy, gam, p_mc, p0_mc, s_t, c_t
complex(rp) cpart
complex(rp), allocatable :: crsource(:,:)
real(rp), pointer :: sigmoid(:,:)
integer ip, ix, iy, ngrid_arr(3), ngrid
logical on_grid, err
character(*), parameter :: r_name = 'fel_spectral_deposit'

!

err_flag = .true.

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
dgrid = wf%dx
p0_mc = fel_p0_mc(beam)

allocate (crsource(ngrid, ngrid))
crsource = 0

scl_w = fel_und_coupling(und, harm) * (mu_0_vac * c_light) * sqrt(2.0_rp) * c_light * delz
scl_w = scl_w / (4 * dgrid * dgrid * beam%slice_spacing)

if (scl_w /= 0) then
  do ip = 1, sl%n
    call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
    if (.not. on_grid) cycle

    p_mc = p0_mc * (1 + sl%pz(ip))
    gam = sqrt(p_mc**2 + 1)
    beta = p_mc / gam

    theta = beam%phi0 + xks1 * sl%z(ip) / beta
    if (harm /= 1) theta = harm * theta

    part = sqrt(faw2(und, sl%x(ip), sl%y(ip))) * scl_w * sl%weight(ip) / gam
    call fel_sincos (theta, s_t, c_t)
    cpart = cmplx(s_t, c_t, rp) * part

    crsource(ix,   iy)   = crsource(ix,   iy)   + (wx * wy) * cpart
    crsource(ix+1, iy)   = crsource(ix+1, iy)   + ((1-wx) * wy) * cpart
    crsource(ix,   iy+1) = crsource(ix,   iy+1) + (wx * (1-wy)) * cpart
    crsource(ix+1, iy+1) = crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cpart
  enddo
endif

! The source filter multiplies the source's transverse spectrum by the sigmoid the
! kernel cache holds for this wavelength (fel_field_kernel_init). The transform pair is
! the source's own here, since nothing else is transformed.

if (und%filter%on) then
  sigmoid => fel_field_kernel_sigmoid(twopi / wf%wavelength)
  if (.not. associated(sigmoid)) then
    call out_io (s_error$, r_name, 'SOURCE FILTER IS ON BUT ITS SIGMOID WAS NOT BUILT FOR THIS WAVELENGTH.', &
                                   'PLEASE REPORT THIS!')
    return
  endif
  call wavefront_fft2 (crsource, wf_fft_forward$, err);  if (err) return
  crsource = crsource * sigmoid
  call wavefront_fft2 (crsource, wf_fft_backward$, err);  if (err) return
  crsource = crsource / real(ngrid*ngrid, rp)
endif

if (allocated(wf%Ey)) then
  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) + 2 * und%pol(1) * crsource
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) + 2 * und%pol(2) * crsource
else
  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) + 2 * crsource
endif

err_flag = .false.

end subroutine fel_spectral_deposit

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_deposit_coarse (und, beam, fs, delz, xks1, err_flag)
!
! Routine to deposit every slice's source into a field sampled coarser than the slices
! (fel_field_sample): plane by plane, in parallel over the planes, each plane summing
! the sources of the slices on it and of the slices between it and its neighbors, each
! slice's source weighted by the slice's distance from the plane, 1 - w on the plane
! below it and w on the plane above, and the sum divided by the slices a plane spans.
! The weights on a plane add to that count, and the source the envelope equation takes
! at a plane is the source per slice there, so the division makes the deposit the
! readback's adjoint under each grid's own measure: a uniform source over the slices
! gives every plane that one slice's source, as one plane a slice does. The source
! filter then runs once on the plane's sum, since the filter is linear, so a step costs
! one transform pair a plane and not one a slice. The per-particle scatter is
! fel_spectral_deposit's, and the scale, the phase and the add of 2*crsource are
! fel_field_step's.
!
! Input:
!   und      -- fel_und_struct: Undulator parameters.
!   beam     -- fel_beam_struct: The beam.
!   fs       -- fel_field_struct: The member, the fundamental.
!   delz     -- real(rp): Step length [m].
!   xks1     -- real(rp): Fundamental radiation wavenumber [1/m].
!
! Output:
!   fs       -- fel_field_struct: The sources added to the record.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_deposit_coarse (und, beam, fs, delz, xks1, err_flag)

type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: fs
type (wavefront_struct), pointer :: wf
real(rp) delz, xks1
logical err_flag

real(rp) dgrid, scl_w, p0_mc
real(rp), pointer :: sigmoid(:,:)
integer ngrid_arr(3), ngrid, nz, nf, m, j, is, is1, is2, ifld
logical any_err
character(*), parameter :: r_name = 'fel_spectral_deposit_coarse'

!

err_flag = .true.
wf => fs%wf
ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
nf = ngrid_arr(3)
nz = size(beam%slice)
m = fs%n_per
dgrid = wf%dx
p0_mc = fel_p0_mc(beam)

scl_w = fel_und_coupling(und, fs%harm) * (mu_0_vac * c_light) * sqrt(2.0_rp) * c_light * delz
scl_w = scl_w / (4 * dgrid * dgrid * beam%slice_spacing)

sigmoid => null()
if (und%filter%on) then
  sigmoid => fel_field_kernel_sigmoid(twopi / wf%wavelength)
  if (.not. associated(sigmoid)) then
    call out_io (s_error$, r_name, 'SOURCE FILTER IS ON BUT ITS SIGMOID WAS NOT BUILT FOR THIS WAVELENGTH.', &
                                   'PLEASE REPORT THIS!')
    return
  endif
endif

any_err = .false.
!$OMP parallel do private(is, is1, is2, ifld) reduction(.or.: any_err)
do j = 1, nf
  call plane_source (j)
enddo
!$OMP end parallel do
err_flag = any_err

!------------------------------------------------------------------------------
contains

!+
! Subroutine plane_source (j)
!
! Routine to sum, filter and add the source of window plane j: the slices from
! (j-1) m + 1 to j m sit on it or above it, with the weight 1 - w, and the slices from
! (j-2) m + 1 to (j-1) m below it, with the weight w of fel_field_sample.
!-

subroutine plane_source (j)

complex(rp), allocatable :: crsource(:,:)
integer j, is_l, is1_l, is2_l, ifld_l
real(rp) w_l
logical err_l

!

allocate (crsource(ngrid, ngrid))
crsource = 0

if (scl_w /= 0) then
  is1_l = (j - 1) * m + 1
  is2_l = min(nz, j * m)
  do is_l = is1_l, is2_l
    w_l = real(mod(is_l - 1, m), rp) / m
    call scatter (beam%slice(is_l), (1 - w_l) / m, crsource)
  enddo
  is1_l = (j - 2) * m + 1
  is2_l = min(nz, (j - 1) * m)
  do is_l = max(1, is1_l), is2_l
    w_l = real(mod(is_l - 1, m), rp) / m
    if (w_l /= 0) call scatter (beam%slice(is_l), w_l / m, crsource)
  enddo
endif

if (associated(sigmoid)) then
  call wavefront_fft2 (crsource, wf_fft_forward$, err_l)
  if (err_l) then
    any_err = .true.;  return
  endif
  crsource = crsource * sigmoid
  call wavefront_fft2 (crsource, wf_fft_backward$, err_l)
  if (err_l) then
    any_err = .true.;  return
  endif
  crsource = crsource / real(ngrid*ngrid, rp)
endif

ifld_l = fel_field_index(fs%slip, j, nf)
if (allocated(wf%Ey)) then
  wf%Ex(:,:,ifld_l) = wf%Ex(:,:,ifld_l) + 2 * und%pol(1) * crsource
  wf%Ey(:,:,ifld_l) = wf%Ey(:,:,ifld_l) + 2 * und%pol(2) * crsource
else
  wf%Ex(:,:,ifld_l) = wf%Ex(:,:,ifld_l) + 2 * crsource
endif

end subroutine plane_source

!+
! Subroutine scatter (sl, weight, crsource)
!
! Routine to add one slice's source times weight to a plane's sum, the per-particle
! scatter of fel_spectral_deposit.
!-

subroutine scatter (sl, weight, crsource)

type (fel_slice_struct) sl
complex(rp) crsource(:,:)
real(rp) weight, part, theta, beta, wx, wy, gam, p_mc, s_t, c_t
complex(rp) cpart
integer ip, ix, iy
logical on_grid

!

do ip = 1, sl%n
  call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
  if (.not. on_grid) cycle

  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam

  theta = beam%phi0 + xks1 * sl%z(ip) / beta
  if (fs%harm /= 1) theta = fs%harm * theta

  part = sqrt(faw2(und, sl%x(ip), sl%y(ip))) * scl_w * sl%weight(ip) / gam * weight
  call fel_sincos (theta, s_t, c_t)
  cpart = cmplx(s_t, c_t, rp) * part

  crsource(ix,   iy)   = crsource(ix,   iy)   + (wx * wy) * cpart
  crsource(ix+1, iy)   = crsource(ix+1, iy)   + ((1-wx) * wy) * cpart
  crsource(ix,   iy+1) = crsource(ix,   iy+1) + (wx * (1-wy)) * cpart
  crsource(ix+1, iy+1) = crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cpart
enddo

end subroutine scatter

end subroutine fel_spectral_deposit_coarse

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine spectral_fft1 (direction, err_flag)
!
! Routine to transform the calling thread's line buffer sp_buf in place, unnormalized,
! by its cached plan. The buffer was filled by the caller and the plans made by
! spectral_fft1_plan_threads before the parallel loop (see wavefront_fft2 for the
! alignment and planner rules this follows).
!
! Input:
!   direction -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   err_flag  -- logical: Set True when the thread holds no plan. False otherwise.
!-

subroutine spectral_fft1 (direction, err_flag)

include 'fftw3.f03'

integer direction
logical err_flag

!

err_flag = .true.
if (sp_cache_n == 0) return

if (direction == wf_fft_forward$) then
  call fftw_execute_dft (sp_plan_fwd, sp_buf, sp_buf)
else
  call fftw_execute_dft (sp_plan_bwd, sp_buf, sp_buf)
endif

err_flag = .false.

end subroutine spectral_fft1

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine spectral_fft1_plan (n, cache_ok)
!
! Routine to hold the calling thread's plan pair and buffer for a transform of length
! n, remaking them when the length changes. FFTW_ESTIMATE, so the plan is a pure
! function of the length and a run reproduces (wavefront_fft2_plan).
!-

subroutine spectral_fft1_plan (n, cache_ok)

include 'fftw3.f03'

integer n
logical cache_ok

!

!$OMP CRITICAL (wavefront_fft_plan_lock)

if (n /= sp_cache_n) then
  if (C_ASSOCIATED(sp_plan_fwd)) call fftw_destroy_plan(sp_plan_fwd)
  if (C_ASSOCIATED(sp_plan_bwd)) call fftw_destroy_plan(sp_plan_bwd)
  if (C_ASSOCIATED(sp_buf_cptr)) call fftw_free(sp_buf_cptr)
  sp_plan_fwd = C_NULL_PTR
  sp_plan_bwd = C_NULL_PTR
  sp_buf_cptr = C_NULL_PTR
  sp_buf => null()
  sp_cache_n = 0

  sp_buf_cptr = fftw_alloc_complex(int(n, C_SIZE_T))
  if (C_ASSOCIATED(sp_buf_cptr)) then
    call C_F_POINTER (sp_buf_cptr, sp_buf, [n])
    sp_plan_fwd = fftw_plan_dft_1d(n, sp_buf, sp_buf, FFTW_FORWARD,  FFTW_ESTIMATE)
    sp_plan_bwd = fftw_plan_dft_1d(n, sp_buf, sp_buf, FFTW_BACKWARD, FFTW_ESTIMATE)
    if (C_ASSOCIATED(sp_plan_fwd) .and. C_ASSOCIATED(sp_plan_bwd)) sp_cache_n = n
  endif
endif

cache_ok = (sp_cache_n == n)

!$OMP END CRITICAL (wavefront_fft_plan_lock)

end subroutine spectral_fft1_plan

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine spectral_fft1_plan_threads (n, err_flag)
!
! Routine to fill every thread's plan cache for a transform of length n, serially,
! before a parallel loop executes them (wavefront_fft2_plan_threads).
!-

subroutine spectral_fft1_plan_threads (n, err_flag)

integer n
logical err_flag
logical any_bad, cache_ok
character(*), parameter :: r_name = 'spectral_fft1_plan_threads'

!

any_bad = .false.

!$OMP PARALLEL private(cache_ok) reduction(.or.: any_bad)
call spectral_fft1_plan (n, cache_ok)
any_bad = any_bad .or. .not. cache_ok
!$OMP END PARALLEL

if (any_bad) then
  call out_io (s_error$, r_name, 'FFTW COULD NOT ALLOCATE A WORK BUFFER OR CREATE A PLAN FOR A TRANSFORM OF \i0\ .', &
               i_array = [n])
endif

err_flag = any_bad

end subroutine spectral_fft1_plan_threads

end module
