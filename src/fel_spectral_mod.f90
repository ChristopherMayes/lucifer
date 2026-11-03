!+
! Module fel_spectral_mod
!
! The spectral field propagator, a research prototype (docs/spectral-field.md). The
! slice propagator (fel_field_step, fel_field_diffract) advances every slice of the
! window by the paraxial kernel exp(-i dz k_perp^2 / 2 k0) at the one reference
! wavenumber k0, and the slippage moves the record by whole slices once 0.8 of a slice
! has accumulated (fel_apply_slippage). This module transforms a member's whole window
! in (x, y, tau) instead. Two things change and the rest does not. The kernel becomes
! exp(-i dz k_perp^2 / 2 k(kappa)) with k(kappa) = k0 + kappa, where kappa is the
! wavenumber conjugate to the window coordinate tau and c kappa is the frequency offset
! from the carrier. That is the term a chirped or short pulse gets wrong at relative
! order dOmega/omega under one k0. And the slippage becomes the phase ramp
! exp(-i kappa delta), delta the step's slippage as a length, exact and continuous, or
! stays the whole-slice rotation under field.spectral_slippage = "slices". The push, the
! deposit and the averaged equations are untouched.
!
! The window coordinate. Window slice i (fel-physics.md sec-window) sits at
! tau = (i - 1) * slice_spacing, so tau grows toward the head, the end light slips
! toward. The envelope rides the carrier exp(i (k0 z - omega0 t)), the convention the
! kernel's sign fixes, so a mode exp(i kappa tau) of the envelope has the carrier's
! frequency plus c kappa and the wavenumber k0 + kappa. The discrete transform of the
! window carries kappa_m = 2 pi m / (n_t slice_spacing) for m up to half of n_t and
! negative beyond, as the transverse table does (fel_field_kernel_init).
!
! The transform along tau is periodic, so light the ramp moves past the head would
! reappear at the tail. A derived window holds the line's slippage ahead of the bunch
! and nothing reaches the tail. Three forms meet the head on a shorter window
! (field.spectral_boundary): the headroom alone, with the run refused once the head
! slice holds more than field.spectral_head_tolerance of the window's energy; a cos^2
! sponge over the last field.spectral_sponge_slices slices, applied before each shift
! of at most a slice, the slippage taken in such chunks, whose removed energy joins the
! member's escaped energy; and zero padding of the transform along tau, eight slices
! beyond the largest shift one pass applies, with what crosses into the padding removed
! and its energy joined the same way. Under "slices"
! slippage the rotation zero-fills the slice it transmits and the slice propagator's
! bank (fel_apply_slippage) keeps the light itself, so nothing wraps and the headroom
! form checks nothing, while the padding then meets the kernel's own spread alone,
! which is below one wavelength a step.
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
use, intrinsic :: iso_c_binding

implicit none

private spectral_fft1, spectral_fft1_plan, spectral_fft1_plan_threads, plane_energy

! The boundary forms (field.spectral_boundary) and the padding of the third.

integer, parameter :: fel_spectral_headroom$ = 1, fel_spectral_sponge$ = 2, fel_spectral_pad$ = 3
integer, parameter :: fel_spectral_n_pad$ = 8    ! Slices of zero padding along tau beyond the
                                                 !   schedule's largest shift (fel_setup_lattice).

!+
! Struct fel_spectral_struct
!
! The spectral propagator's switches, resolved at setup from the [field] keys, and what
! it measured on the way: the largest share of the window's energy the head slice held,
! and where.
!-

type fel_spectral_struct
  logical :: on = .false.               ! field.propagator = "spectral".
  logical :: ramp = .true.              ! The slippage as the phase ramp. False: whole-slice rotation.
  integer :: boundary = fel_spectral_headroom$
  integer :: n_sponge = 0               ! Sponge slices at the head.
  integer :: n_pad = 0                  ! Zero-padding slices along tau.
  real(rp) :: head_tolerance = 1e-3_rp  ! The headroom form's refusal threshold.
  real(rp) :: head_max = 0              ! Largest head-slice share seen.
  real(rp) :: z_head_max = 0            ! Where [m].
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
! Subroutine fel_spectral_setup (sp, propagator, slippage, boundary, n_sponge, head_tolerance)
!
! Routine to resolve the [field] keys of the spectral propagator into its switches. The
! key table has already refused a name outside each key's set. fel_spectral_size then
! sizes the sponge and the padding from the window and the schedule.
!
! Input:
!   propagator     -- character(*): field.propagator, "slices" or "spectral".
!   slippage       -- character(*): field.spectral_slippage, "ramp" or "slices".
!   boundary       -- character(*): field.spectral_boundary, "headroom", "sponge" or "pad".
!   n_sponge       -- integer: field.spectral_sponge_slices, 0 derived.
!   head_tolerance -- real(rp): field.spectral_head_tolerance.
!
! Output:
!   sp             -- fel_spectral_struct: The resolved switches.
!-

subroutine fel_spectral_setup (sp, propagator, slippage, boundary, n_sponge, head_tolerance)

type (fel_spectral_struct) sp
character(*) propagator, slippage, boundary
integer n_sponge
real(rp) head_tolerance

!

sp%on = (propagator == 'spectral')
sp%ramp = (slippage /= 'slices')
select case (boundary)
case ('sponge');   sp%boundary = fel_spectral_sponge$
case ('pad');      sp%boundary = fel_spectral_pad$
case default;      sp%boundary = fel_spectral_headroom$
end select
sp%head_tolerance = head_tolerance
sp%n_sponge = n_sponge
sp%n_pad = 0

end subroutine fel_spectral_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_size (sp, shift_max, nslice, err_flag)
!
! Routine to size the sponge and the padding. The padding is eight slices beyond the
! largest shift one pass applies, an element's end slippage with its last step's, in
! slices, since the ramp takes it in one pass. The sponge takes the slippage in chunks
! of a slice (fel_spectral_shift), so its width is the window's alone: a count of 0
! derives an eighth of the window, at least 4 slices, and a count not below the window
! is refused. A window of one slice, the steady state, takes no boundary, since
! nothing slips.
!
! Input:
!   shift_max -- real(rp): The largest shift one pass applies [slices].
!   nslice    -- integer: The window's slices.
!
! Output:
!   sp        -- fel_spectral_struct: n_sponge and n_pad set.
!   err_flag  -- logical: Set True on a sponge not below the window.
!-

subroutine fel_spectral_size (sp, shift_max, nslice, err_flag)

type (fel_spectral_struct) sp
real(rp) shift_max
integer nslice, n_shift
logical err_flag
character(*), parameter :: r_name = 'fel_spectral_size'

!

err_flag = .true.
n_shift = ceiling(shift_max)

if (nslice == 1) then
  sp%n_sponge = 0
  sp%n_pad = 0
  err_flag = .false.
  return
endif

if (sp%boundary == fel_spectral_sponge$) then
  if (sp%n_sponge == 0) sp%n_sponge = max(4, nslice / 8)
  if (sp%n_sponge >= nslice) then
    call out_io (s_error$, r_name, 'THE SPONGE OF \i0\ SLICES IS NOT BELOW THE WINDOW''S \i0\ SLICES.', &
                 i_array = [sp%n_sponge, nslice])
    return
  endif
else
  sp%n_sponge = 0
endif

if (sp%boundary == fel_spectral_pad$) sp%n_pad = fel_spectral_n_pad$ + n_shift

err_flag = .false.

end subroutine fel_spectral_size

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_step (sp, und, beam, ff, dz, slippage, err_flag)
!
! Routine to do the field half of one undulator step under the spectral propagator,
! in the slice propagator's order (fel_field_step, then fel_apply_slippage): each
! member's window is propagated over the step by the kernel, the source is deposited
! slice by slice from the just-advanced particles (fel_spectral_deposit), and the
! slippage then moves the field and the fresh source together, as the phase ramp with
! the boundary (fel_spectral_shift), or by the caller's rotation under "slices"
! slippage. fel_track_und_step has pushed the particles and skipped its own field solve.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switches.
!   und      -- fel_und_struct: The segment.
!   beam     -- fel_beam_struct: The beam, pushed over the step.
!   ff(:)    -- fel_field_struct: The field set.
!   dz       -- real(rp): Step length [m].
!   slippage -- real(rp): The step's slippage [fundamental wavelengths].
!
! Output:
!   ff(:)    -- fel_field_struct: Every member advanced.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_step (sp, und, beam, ff, dz, slippage, err_flag)

type (fel_spectral_struct) sp
type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
real(rp) dz, slippage
logical err_flag

real(rp) ks
integer io, is, nslice_f
logical any_err, err

!

err_flag = .true.
ks = twopi / ff(1)%wf%wavelength

do io = 1, size(ff)
  call fel_spectral_apply (sp, ff(io), dz, 0.0_rp, .true., err)
  if (err) return
  nslice_f = size(ff(io)%wf%Ex, 3)
  any_err = .false.
  !$OMP parallel do private(err) reduction(.or.: any_err)
  do is = 1, size(beam%slice)
    call fel_spectral_deposit (und, beam, beam%slice(is), ff(io)%wf, &
                               fel_field_index(ff(io)%slip, is, nslice_f), dz, ff(io)%harm, ks, err)
    any_err = any_err .or. err
  enddo
  !$OMP end parallel do
  if (any_err) return
  call fel_spectral_shift (sp, ff(io), slippage, err)
  if (err) return
enddo

err_flag = .false.

end subroutine fel_spectral_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_drift (sp, ff, dz, slippage, err_flag)
!
! Routine to propagate every member of the field set through an interlude under the
! spectral propagator: the kernel over the interlude's length with no source, then
! the slippage the schedule charges it with the boundary (fel_spectral_shift). This is
! what wavefront_drift does for the slice propagator, with the kernel at k(kappa) and
! the slippage as the ramp.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switches.
!   ff(:)    -- fel_field_struct: The field set.
!   dz       -- real(rp): The light's path [m].
!   slippage -- real(rp): The interlude's slippage [fundamental wavelengths].
!
! Output:
!   ff(:)    -- fel_field_struct: Every member advanced.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_spectral_drift (sp, ff, dz, slippage, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: ff(:)
real(rp) dz, slippage
logical err_flag
integer io

!

err_flag = .true.
do io = 1, size(ff)
  call fel_spectral_apply (sp, ff(io), dz, 0.0_rp, .true., err_flag)
  if (err_flag) return
  call fel_spectral_shift (sp, ff(io), slippage, err_flag)
  if (err_flag) return
enddo
err_flag = .false.

end subroutine fel_spectral_drift

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_shift (sp, fs, slippage, err_flag)
!
! Routine to slip one member by the ramp with the boundary form in force, in passes
! along tau alone. Under the sponge the slippage is taken in chunks of at most one
! slice with the sponge applied before each, so what a chunk carries past the head is
! the sponge's zero on the head slice and nothing else; an interlude's whole slippage,
! tens of slices on a line like Aramis, would otherwise carry its unsponged slices
! round to the tail. The other forms take the slippage in one pass, the padding sized
! to hold it, and the headroom check follows. Under "slices" slippage the sponge alone
! runs, before the caller's rotation, and nothing else.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switches.
!   fs       -- fel_field_struct: The member.
!   slippage -- real(rp): The slippage to apply [fundamental wavelengths].
!
! Output:
!   fs       -- fel_field_struct: The window slipped, u_escaped grown by what the
!                 boundary removed.
!   err_flag -- logical: Set True on a transform failure or the headroom refusal.
!-

subroutine fel_spectral_shift (sp, fs, slippage, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
real(rp) slippage, chunk, left
logical err_flag
integer n_chunk, ic

!

err_flag = .false.
if (.not. sp%ramp) then
  call fel_spectral_sponge (sp, fs)
  return
endif
if (.not. fs%slip%timerun) return

if (sp%boundary == fel_spectral_sponge$) then
  n_chunk = max(1, ceiling(abs(slippage) / fs%slip%sample - 1e-9_rp))
  left = slippage
  do ic = 1, n_chunk
    chunk = slippage / n_chunk
    if (ic == n_chunk) chunk = left
    left = left - chunk
    call fel_spectral_sponge (sp, fs)
    call fel_spectral_apply (sp, fs, 0.0_rp, chunk, .false., err_flag)
    if (err_flag) return
  enddo
else
  call fel_spectral_apply (sp, fs, 0.0_rp, slippage, .false., err_flag)
  if (err_flag) return
  call fel_spectral_head_check (sp, fs, err_flag)
endif

end subroutine fel_spectral_shift

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_apply (sp, fs, dz, slippage, transverse, err_flag)
!
! Routine to apply the kernel over dz and the ramp for the slippage to one member's
! whole window, under sp%ramp, or the ramp alone when transverse is False, which is the
! step's pass after the deposit (fel_spectral_step). The record is read and written
! in window order through fel_field_index, so a rotated record is handled as the slice
! propagator handles it, and under "slices" slippage the caller rotates it afterwards
! (fel_apply_slippage) as before.
!
! The passes, per polarization: the forward transverse transform of every slice; on
! every line of the transverse spectrum the forward transform along tau, the kernel
! and ramp factors, and the inverse; the inverse transverse transform of every slice.
! The transforms are unnormalized, so the result is scaled by n_t and by n_grid^2. The
! ramp alone is diagonal in tau's conjugate and needs no transverse pass, so it runs
! on the lines of the field itself.
!
! Input:
!   sp         -- fel_spectral_struct: The propagator's switches.
!   fs         -- fel_field_struct: The member: its wavefront, slippage state and harmonic.
!   dz         -- real(rp): Step length [m].
!   slippage   -- real(rp): The slippage to apply [fundamental wavelengths].
!   transverse -- logical: Apply the kernel, with the transverse passes. False applies
!                   the ramp alone.
!
! Output:
!   fs         -- fel_field_struct: The window advanced, u_escaped grown by what the
!                   padding removed.
!   err_flag   -- logical: Set True on a transform failure. False otherwise.
!-

subroutine fel_spectral_apply (sp, fs, dz, slippage, transverse, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
type (wavefront_struct), pointer :: wf
real(rp) dz, slippage
logical transverse, err_flag

complex(rp), allocatable :: ax(:,:), ay(:,:)
real(rp) k0, dk, dtau, delta, kappa, kt, e_scale, t_slice, e_pad, grid_norm
integer, allocatable :: idx(:)
integer ngrid_arr(3), ngrid, nz, nt, ix, m, mx, mt, j
logical err, any_err

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

! The slippage as a length, in the window. A steady-state record never slips.

delta = 0
if (sp%ramp .and. fs%slip%timerun) delta = slippage * wf%wavelength * fs%harm

! The energy of one plane, as fel_field_diag and fel_apply_slippage take it: the sum of
! |E|^2 over the grid times dx^2 / 2 Z0 is the slice's power, and times the slice's
! light time its energy. A harmonic member's wavelength is the fundamental's over harm,
! and the window is the fundamental's, so the light time is scaled back up. In the
! transverse spectrum the sum is n_grid^2 times the field's (Parseval).

e_scale = wf%dx**2 / (2 * (mu_0_vac * c_light))
t_slice = fs%slip%sample * wf%wavelength * fs%harm / c_light
grid_norm = 1
if (transverse) grid_norm = real(ngrid*ngrid, rp)

! The two kernel tables (module header). Each column m is one kappa, each row one
! transverse mode, and the ramp rides the x table. The ramp alone is a column vector.

allocate (ax(ngrid, nt), ay(ngrid, nt), idx(nz))
do m = 1, nt
  mt = m - 1;  if (2*mt > nt-1) mt = mt - nt
  kappa = twopi * mt / (nt * dtau)
  kt = k0 + kappa
  do ix = 1, ngrid
    mx = ix - 1;  if (2*mx > ngrid-1) mx = mx - ngrid
    if (transverse) then
      ax(ix, m) = exp(cmplx(0.0_rp, -dz * (mx * dk)**2 / (2 * kt) - kappa * delta, rp))
      ay(ix, m) = exp(cmplx(0.0_rp, -dz * (mx * dk)**2 / (2 * kt), rp))
    else
      ax(ix, m) = exp(cmplx(0.0_rp, -kappa * delta, rp))
      ay(ix, m) = 1
    endif
  enddo
enddo
do j = 1, nz
  idx(j) = fel_field_index(fs%slip, j, nz)
enddo

! The plans every thread will execute, made serially (wavefront_fft2_plan_threads).

if (transverse) then
  call wavefront_fft2_plan_threads (ngrid, ngrid, err);  if (err) return
endif
call spectral_fft1_plan_threads (nt, err);  if (err) return

e_pad = 0
call apply_component (wf%Ex, err);  if (err) return
if (allocated(wf%Ey)) then
  call apply_component (wf%Ey, err);  if (err) return
endif
fs%slip%u_escaped = fs%slip%u_escaped + e_pad * e_scale * t_slice / grid_norm

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine apply_component (e, err)
!
! Routine to put one polarization's planes through the passes of the header. The
! padding's energy is summed into e_pad over the lines, in the units of the pass.
!-

subroutine apply_component (e, err)

complex(wf_rp) e(:,:,:)
logical err
integer is_l, ix_l, iy_l, m_l, j_l
logical err_l

!

err = .true.
any_err = .false.
if (transverse) then
  !$OMP parallel do private(err_l) reduction(.or.: any_err)
  do is_l = 1, nz
    call wavefront_fft2 (e(:,:,is_l), wf_fft_forward$, err_l)
    any_err = any_err .or. err_l
  enddo
  !$OMP end parallel do
  if (any_err) return
endif

!$OMP parallel do private(ix_l, j_l, m_l, err_l) reduction(.or.: any_err) reduction(+: e_pad)
do iy_l = 1, ngrid
  do ix_l = 1, ngrid
    do j_l = 1, nz
      sp_buf(j_l) = e(ix_l, iy_l, idx(j_l))
    enddo
    if (nt > nz) sp_buf(nz+1:nt) = 0
    call spectral_fft1 (wf_fft_forward$, err_l)
    any_err = any_err .or. err_l
    do m_l = 1, nt
      sp_buf(m_l) = sp_buf(m_l) * (ax(ix_l, m_l) * ay(iy_l, m_l))
    enddo
    call spectral_fft1 (wf_fft_backward$, err_l)
    any_err = any_err .or. err_l
    do j_l = 1, nz
      e(ix_l, iy_l, idx(j_l)) = sp_buf(j_l) / real(nt, rp)
    enddo
    do j_l = nz + 1, nt
      e_pad = e_pad + (real(sp_buf(j_l), rp)**2 + aimag(sp_buf(j_l))**2) / real(nt, rp)**2
    enddo
  enddo
enddo
!$OMP end parallel do
if (any_err) return

if (transverse) then
  !$OMP parallel do private(err_l) reduction(.or.: any_err)
  do is_l = 1, nz
    call wavefront_fft2 (e(:,:,is_l), wf_fft_backward$, err_l)
    e(:,:,is_l) = e(:,:,is_l) / real(ngrid*ngrid, rp)
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
! Subroutine fel_spectral_sponge (sp, fs)
!
! Routine to apply the cos^2 sponge at the window's head when it is the boundary
! (module header), before each shift of at most a slice, so what the shift carries
! past the head is the sponge's zero. The weight falls from just under 1 on the sponge's
! first slice to 0 on the head slice, in window order, and the field is multiplied by
! it, so a slice keeps w^2 of its energy and the rest joins u_escaped. A window of one
! slice takes no sponge.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switches.
!   fs       -- fel_field_struct: The member.
!
! Output:
!   fs       -- fel_field_struct: The sponge applied, u_escaped grown by what it removed.
!-

subroutine fel_spectral_sponge (sp, fs)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
type (wavefront_struct), pointer :: wf

real(rp) e_scale, t_slice, w, e_before
integer ngrid_arr(3), nz, j, ifld

!

if (sp%boundary /= fel_spectral_sponge$) return
wf => fs%wf
ngrid_arr = wavefront_shape(wf)
nz = ngrid_arr(3)
if (nz == 1) return

e_scale = wf%dx**2 / (2 * (mu_0_vac * c_light))
t_slice = fs%slip%sample * wf%wavelength * fs%harm / c_light

do j = nz - sp%n_sponge + 1, nz
  ifld = fel_field_index(fs%slip, j, nz)
  w = cos(0.5_rp * pi * real(j - (nz - sp%n_sponge), rp) / sp%n_sponge)**2
  e_before = plane_energy(wf, ifld)
  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * w
  if (allocated(wf%Ey)) wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * w
  fs%slip%u_escaped = fs%slip%u_escaped + (e_before - plane_energy(wf, ifld)) * e_scale * t_slice
enddo

end subroutine fel_spectral_sponge

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_spectral_head_check (sp, fs, err_flag)
!
! Routine to hold the window's head after a pass under the headroom form (module
! header): the run is refused once the head slice holds more than sp%head_tolerance of
! the window's energy, and the largest share seen is recorded with its position. The
! check is the ramp's alone, since the whole-slice rotation banks and zero-fills what
! reaches the head, and a window of one slice takes none.
!
! Input:
!   sp       -- fel_spectral_struct: The propagator's switches.
!   fs       -- fel_field_struct: The member.
!
! Output:
!   sp       -- fel_spectral_struct: head_max and z_head_max updated.
!   err_flag -- logical: Set True on the refusal. False otherwise.
!-

subroutine fel_spectral_head_check (sp, fs, err_flag)

type (fel_spectral_struct) sp
type (fel_field_struct), target :: fs
type (wavefront_struct), pointer :: wf
logical err_flag

real(rp) e_head, e_total
integer ngrid_arr(3), nz, j
character(*), parameter :: r_name = 'fel_spectral_head_check'

!

err_flag = .false.
if (sp%boundary /= fel_spectral_headroom$ .or. .not. sp%ramp) return
wf => fs%wf
ngrid_arr = wavefront_shape(wf)
nz = ngrid_arr(3)
if (nz == 1) return

e_head = plane_energy(wf, fel_field_index(fs%slip, nz, nz))
e_total = 0
do j = 1, nz
  e_total = e_total + plane_energy(wf, j)
enddo
if (e_total <= 0) return

if (e_head / e_total > sp%head_max) then
  sp%head_max = e_head / e_total
  sp%z_head_max = wf%s_position
endif
if (e_head > sp%head_tolerance * e_total) then
  call out_io (s_error$, r_name, &
        'THE WINDOW''S HEAD SLICE HOLDS \es10.2\ OF THE WINDOW''S ENERGY AT s = \f0.3\ m, ABOVE', &
        'FIELD.SPECTRAL_HEAD_TOLERANCE = \es10.2\ . THE TRANSFORM ALONG TAU IS PERIODIC, SO LIGHT', &
        'PAST THE HEAD WOULD REAPPEAR AT THE TAIL. POSSIBLE SOLUTION: A LONGER WINDOW, OR', &
        'FIELD.SPECTRAL_BOUNDARY = "sponge" OR "pad".', &
        r_array = [e_head / e_total, wf%s_position, sp%head_tolerance])
  err_flag = .true.
endif

end subroutine fel_spectral_head_check

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function plane_energy (wf, ifld) result (e_sum)
!
! Routine to sum |E|^2 over one record plane, both polarizations, in grid units.
!-

function plane_energy (wf, ifld) result (e_sum)

type (wavefront_struct) wf
integer ifld
real(rp) e_sum

!

e_sum = sum(real(wf%Ex(:,:,ifld), rp)**2 + aimag(wf%Ex(:,:,ifld))**2)
if (allocated(wf%Ey)) e_sum = e_sum + sum(real(wf%Ey(:,:,ifld), rp)**2 + aimag(wf%Ey(:,:,ifld))**2)

end function plane_energy

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
