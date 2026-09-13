!+
! Module fel_unaveraged_mod
!
! The unaveraged mode (fel-physics.md sec-unaveraged):
! particles integrated through the undulator's real field (the full Newton-Lorentz
! quiver, no period averaging) with the radiation field as a co-evolving kick. There
! is no resonance approximation. What the cost per step buys is physics the averaged
! map cannot reach: the full quiver dynamics, the energy accounting the beam actually
! pays, polarization-agnostic coupling and arbitrary harmonic content. The mode is
! also an independent check on the averaged path, since the coupling factor fc, the harmonic
! content and the entry/exit behavior of the averaged mode become measured outputs
! here instead of assumed inputs. MINERVA (Freund and van der Slot) is the published
! existence proof for this physics, and only its published work was used here. This mode differs by keeping the grid field and
! the Lorentz force where MINERVA evaluates modal fields. The cost is priced in
! fel-physics.md and measured in examples/saturation_demo.
!
! DO NOT introduce fc, faw, or any other period-averaged coupling quantity into this
! path. Those are what this mode measures, and their appearance here would make the
! check circular. The harness greps this file for them.
!
! The step, per substep delta (Strang split, second order):
!
!   1. half magnetic push: RK4 through the analytic normalized field b = curl(a),
!      a = e*A/(m_e c), with gamma exact (B does no work). The vector potential carries
!      a sin^2 amplitude envelope g(s) over n_ramp periods at each end, with g' terms
!      retained in b, so the quiver builds adiabatically and vanishes at the segment
!      ends -- the averaged<->unaveraged handoff (the K/gamma chart hazard) happens
!      where the two momentum conventions coincide. Each handoff also applies the
!      ramp's slippage compensation as a discrete phase jump (unavg_ramp_phase_jump:
!      the ramped ends slip ~3 rad of optical phase less than the contracted
!      hard-edge element. Uncompensated, that scrambles the bunching-to-field phase
!      of every pre-bunched segment entry).
!   2. radiation kick + source deposit at the substep midpoint:
!        dgamma/ds = -Re[W conj(j)]/(u_s m_e),   W = -i*Ehat*e^{i Psi}
!      with Psi = (phi0 - ku*s_local) - ks*tau the optical phase (phi0 is the beam's
!      ponderomotive reference, advanced at the averaged rate, so diagnostics stay
!      comparable), and the polarization-basis current
!        j = u_x (planar),  (u_x - i u_y)/sqrt(2) (helical).
!      The source is the same SVEA deposit as the averaged solver with the coupling
!      Removed and the actual quiver current in its place:
!        src += i e^{-i Psi} * j * (Z0 c dz /(2 dgrid^2 Ds)) * w/u_s
!      added to the record as +2*src and the pair then carried through the shared pure
!      diffraction (fel_field_diffract): E' = D (E + 2 src). The /u_s (where the
!      averaged solver has Genesis's /gamma) makes the kick/deposit pair exact energy
!      duals per substep (same operands, same bilinear weights, and the source landing
!      on the record the kick read before the unitary diffraction), so the ledger closes
!      to the physical spontaneous-emission term and rounding, by construction. With the
!      source added after the diffraction the two sides carried different cross terms,
!      first order in the substep's diffraction phase (FINDINGS 7.86). Period-averaging the
!      pair reproduces the averaged mode's fc to O(1-beta_par) ~ 5e-9 (the JJ factor
!      emerges from the figure-8).
!   3. half magnetic push.
!
! Units: E in V/m (wavefront convention), m_electron in eV, b in 1/m. The physical
! field of the scalar envelope: planar E_x = Re[-i Ehat e^{i Psi}]; helical
! (E_x, E_y) = (Re[-i Ehat e^{i Psi}], Re[Ehat e^{i Psi}])/sqrt(2), the pair (1, +i)
! that rotates with the electron (jhat = (u_x - i u_y)/sqrt(2) below is its dual). Both
! give intensity |Ehat|^2/(2 Z0), so the power diagnostic is mode-independent.
!
! The magnetic push is classical RK4 on the exact z-ODEs in kinetic variables,
! chosen on merit: the currency here is short-probe
! accuracy, and 4th order is what makes fc measurable at 6e-4 with 20
! steps/period (a 2nd-order symplectic scheme needs ~100 steps/period to match, and
! no explicit symplectic method exists for this non-separable Hamiltonian without
! paying implicit iterations). The structural cost is measured, not argued: gamma is
! conserved exactly by construction (B does no work: gamma changes only in the kick),
! and over the longest benchmark segment (266 periods, 5320 steps) the dark-run
! emittance drifts by <= 3.3e-6, orders below every check. If production-length
! unaveraged runs ever appear (oscillator passes), revisit with a symplectic
! composition. The ballistic check is the instrument that will say when.
!
! Parallel over slices with the averaged step's own guarantees (the OpenMP design's
! design): disjoint particle arrays and field slices per iteration, serial kernel
! init, threadprivate FFT plans, per-slice energy summed in fixed order. Results
! are bit-identical across thread counts, and the harness checks it.
!
! The device path (fel_device_mod, doc/validation.md val-device-unaveraged). With a
! backend armed the beam and the field stay resident for the whole segment and one
! record step is one command buffer of nsub substeps. Three kernels are this mode's own
! and the rest of a substep is the kernels the averaged path already had. Residency
! starts inside the first record step rather than at the walk's element entry, because
! the entry handoff moves z and the chart conversion has to see the moved value, and it
! ends inside the last one for the same reason. The quantities that do not depend on the
! particle are built here in FP64 and uploaded once a substep, which is the hoist
! FINDINGS 7.37 records. With fp32_check on, the device takes the instrument's twin role
! instead and the FP64 path below runs untouched.
!-

module fel_unaveraged_mod

use fel_track_mod

implicit none

!+
! Struct fel_unavg_struct
!
! Per-segment state of the unaveraged tracker, carried by the driver across the
! element's record steps: position into the segment (the carrier phase ku*s and the
! ramp envelope are functions of it), the substep count, and the running energy
! ledger.
!-

type fel_unavg_struct
  real(rp) :: s = 0                ! Distance into the segment [m].
  real(rp) :: l = 0                ! Segment length [m].
  real(rp) :: l_ramp = 0           ! sin^2 ramp length at each end [m].
  integer :: nsub = 1              ! Substeps per record step.
  real(rp) :: dsub = 0             ! Substep length [m].
  logical :: active = .false.      ! Between entry and exit handoff.
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_unavg_setup (und, ustate, l, dz_record, steps_per_period, ramp_periods, err_flag)
!
! Routine to size one segment's unaveraged state: substeps per record step from the
! requested steps per period (10 is the floor our own fc convergence supports, 20 the
! default), ramps
! from the requested periods. Refuses a ramp pair longer than the segment and a record
! step that does not hold an integer substep count.
!
! Input:
!   und              -- fel_und_struct: Undulator parameters.
!   l                -- real(rp): Segment length [m].
!   dz_record        -- real(rp): Record step (the element's ds_step) [m].
!   steps_per_period -- integer: Substeps per undulator period (floor 10).
!   ramp_periods     -- real(rp): sin^2 end-ramp length in periods (0 = hard edge).
!
! Output:
!   ustate           -- fel_unavg_struct: Substep grid, ramp geometry, work arrays.
!   err_flag         -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_unavg_setup (und, ustate, l, dz_record, steps_per_period, ramp_periods, err_flag)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
real(rp) l, dz_record, ramp_periods, lambda_w
integer steps_per_period
logical err_flag
character(*), parameter :: r_name = 'fel_unavg_setup'

!

err_flag = .true.
lambda_w = twopi / und%ku

ustate%s = 0
ustate%l = l
ustate%l_ramp = ramp_periods * lambda_w
ustate%active = .false.

if (2 * ustate%l_ramp > l) then
  call out_io (s_error$, r_name, 'UNAVERAGED RAMPS (\es10.2\ m EACH) DO NOT FIT THE SEGMENT (\es10.2\ m).', &
                                 r_array = [ustate%l_ramp, l])
  return
endif

if (steps_per_period < 10) then
  call out_io (s_error$, r_name, 'UNAVERAGED_STEPS_PER_PERIOD BELOW THE CONVERGENCE FLOOR OF 10: \i0\ ', &
                                 i_array = [steps_per_period])
  return
endif

ustate%nsub = max(1, nint(dz_record / (lambda_w / steps_per_period)))
ustate%dsub = dz_record / ustate%nsub

err_flag = .false.

end subroutine fel_unavg_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_unavg_envelope (ustate, s, gp) result (g)
!
! The undulator amplitude envelope at s into the segment: sin^2 up over l_ramp,
! flat 1, sin^2 down over the last l_ramp. Amplitude and slope are continuous (the
! slope gp feeds the ramp-induced field terms in fel_unavg_bfield). l_ramp = 0 is
! the hard-edge mutation configuration. The handoff check exists to catch it.
!
! Input:
!   ustate -- fel_unavg_struct: Ramp geometry.
!   s      -- real(rp): Position inside the segment [m].
!
! Output:
!   gp     -- real(rp): The envelope derivative dg/ds [1/m].
!   g      -- real(rp): The field envelope g(s) (sin^2 ramps, 1 in the body).
!-

function fel_unavg_envelope (ustate, s, gp) result (g)

type (fel_unavg_struct) ustate
real(rp) s, g, gp, arg

!

g = 1;  gp = 0
if (ustate%l_ramp <= 0) return

if (s < ustate%l_ramp) then
  arg = pi * s / (2 * ustate%l_ramp)
  g = sin(arg)**2
  gp = (pi / ustate%l_ramp) * sin(arg) * cos(arg)
elseif (s > ustate%l - ustate%l_ramp) then
  arg = pi * (ustate%l - s) / (2 * ustate%l_ramp)
  g = sin(arg)**2
  gp = -(pi / ustate%l_ramp) * sin(arg) * cos(arg)
endif

end function fel_unavg_envelope

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_unavg_bfield (und, ustate, x, y, s, bx, by, bz)
!
! The normalized magnetostatic field b = e*B/(m_e c) = curl(a) [1/m] at s into the
! segment, from the analytic vector potential with the ramp envelope g(s):
!
!   planar:  a_x = a0 g(s) cos(ku s) cosh(ku y),          a0 = sqrt(2) aw
!   helical: a_x = a0 g(s) cos(ku s) (1 + ku^2 r^2/4),    a0 = aw
!            a_y = a0 g(s) sin(ku s) (1 + ku^2 r^2/4)
!
! (aw rms; peak = rms*sqrt(2) planar, = rms helical). The curl retains the g' ramp
! terms. The transverse profiles are the near-axis models whose ponderomotive-average
! focusing reproduces the averaged mode's natural-focusing split exactly (planar
! kx = 0, ky = ku^2; helical kx = ky = ku^2/2) -- checked in sec-unaveraged.
!
! The four quantities that depend on s and not on the particle arrive as arguments:
! the envelope g and its slope gp, and cos(ku s), sin(ku s). Every particle at one RK
! stage shares them, so the caller evaluates them once per stage rather than once per
! particle per stage (FINDINGS 7.37). Taking them as arguments rather than computing
! them here is what makes that structural: this routine can no longer be the place a
! per-particle transcendental hides.
!
! Input:
!   und        -- fel_und_struct: Undulator parameters (helicity, tilt frame).
!   x, y       -- real(rp): Transverse position [m].
!   g, gp      -- real(rp): The ramp envelope and its slope at this s (fel_unavg_envelope).
!   c_u, s_u   -- real(rp): cos(und%ku * s) and sin(und%ku * s) at this s.
!
! Output:
!   bx, by, bz -- real(rp): The analytic undulator field B = curl(a), with the
!                   envelope's g' terms so the ramps stay divergence-free [T].
!-

subroutine fel_unavg_bfield (und, x, y, g, gp, c_u, s_u, bx, by, bz)

type (fel_und_struct) und
real(rp) x, y, g, gp, c_u, s_u, bx, by, bz
real(rp) a0, fperp, xl, yl, bt

!

! A tilted planar element: evaluate the untilted potential in the wiggle frame
! (coordinates rotated in), rotate b back out. sin_t = 0 skips both rotations.

xl = x;  yl = y
if (und%sin_t /= 0) then
  xl =  und%cos_t * x + und%sin_t * y
  yl = -und%sin_t * x + und%cos_t * y
endif

if (und%helical) then
  a0 = und%aw
  fperp = 1 + (und%ku**2 / 4) * (xl*xl + yl*yl)
  ! b = curl(a) with a_x = a0 g c_u fperp, a_y = a0 g s_u fperp:
  bx = -a0 * (gp * s_u + g * und%ku * c_u) * fperp                        ! -d(a_y)/ds
  by =  a0 * (gp * c_u - g * und%ku * s_u) * fperp                        !  d(a_x)/ds
  bz =  a0 * g * (und%ku**2 / 2) * (s_u * xl - c_u * yl)                  !  d(a_y)/dx - d(a_x)/dy
else
  a0 = sqrt(2.0_rp) * und%aw
  bx = 0
  by = a0 * (gp * c_u - g * und%ku * s_u) * cosh(und%ku * yl)           ! d(a_x)/ds
  bz = -a0 * g * c_u * und%ku * sinh(und%ku * yl)                       ! -d(a_x)/dy
endif

if (und%sin_t /= 0) then        ! Rotate the transverse field components back out.
  bt = und%cos_t * bx - und%sin_t * by
  by = und%sin_t * bx + und%cos_t * by
  bx = bt
endif

end subroutine fel_unavg_bfield

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_unavg_step (und, ustate, beam, wf, slip, dz_record, first, last, dE_beam, err_flag)
!
! Routine to advance the beam and field one record step of the unaveraged mode:
! nsub Strang substeps of (half B push, radiation kick + deposit, half B push, field
! diffract + add source), per slice, serially. first/last mark the segment's record
! boundaries: entry asserts and takes the averaged chart (the ramp guarantees a = 0
! there, so kinetic momentum equals the stored canonical-convention px), exit hands
! it back and restores the flag. beam%phi0 advances at the averaged ponderomotive
! rate so every diagnostic row stays comparable. The optical carrier used internally
! is Psi = (phi0 - ku*s) - ks*tau. dE_beam returns the weighted particle energy
! change of this step [J] for the energy-ledger check.
!
! Input:
!   und       -- fel_und_struct: Undulator parameters.
!   ustate    -- fel_unavg_struct: Substep grid and work arrays.
!   beam      -- fel_beam_struct: The beam in the quiver chart.
!   ff(:)     -- fel_field_struct: The field set. This mode carries one member of one
!                  plane, so ff(1) holds the field and the rotating record.
!   dz_record -- real(rp): The record step to advance by [m].
!   first     -- logical: True on the segment's first record step (entry handoff).
!   last      -- logical: True on the last (exit handoff and ramp phase jump).
!   fp32      -- fel_fp32_struct: The lockstep instrument, dark when it is off.
!   dev       -- fel_device_struct: The device backend, dark when it is off. With the
!                  instrument off it runs the step; with it on it is the twin.
!
! Output:
!   beam      -- fel_beam_struct: Advanced by Newton-Lorentz RK4 through the field.
!   ff(:)     -- fel_field_struct: Sources deposited, records diffracted.
!   dE_beam   -- real(rp): The step's kick-side beam energy change [J] (ledger).
!   dU_spont  -- real(rp): The step's spontaneous source energy [J] (ledger).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_unavg_step (und, ustate, beam, ff, dz_record, first, last, dE_beam, dU_spont, fp32, dev, err_flag)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_slip_struct), pointer :: slip
type (fel_fp32_struct), target :: fp32
type (fel_device_struct) dev

real(rp) dz_record, dE_beam, dU_spont
logical first, last, err_flag

real(rp), allocatable :: ux(:), uy(:), xx(:), yy(:), tau(:), gam(:), dE_slice(:), dU_sp_slice(:)
real(rp), allocatable :: dEturn_slice(:), dev_s0(:,:,:)
real(rp), allocatable :: kst(:,:,:)
real(rp), allocatable :: cx(:), cy(:), cpx(:), cpy(:), cz(:), cpz(:)
complex(rp), allocatable :: crsource(:,:), crsource_y(:,:)
real(rp) p0_mc, gamma0b, inv_beta0, ks, dsub, s_sub, phi0_rate_avg, scl_u, dgrid
real(rp) u_s, wx, wy, psi_mid, dgam, p_mc, beta
real(rp) fq(4,4)     ! The stage field factors, rebuilt per push. Private per thread.
complex(rp), pointer :: exp_k2_p(:,:)
complex(rp) ehat, jhat, wphasor, cdep, ehat_y, cph
logical two_pol
integer is, ip, isub, nslice, ifld, ix, iy, ngrid_arr(3), ngrid
logical on_grid, err, any_err, dev_twin
integer n_s0
character(*), parameter :: r_name = 'fel_unavg_step'

!

err_flag = .true.
dE_beam = 0
dU_spont = 0

wf => ff(1)%wf
slip => ff(1)%slip
nslice = size(wf%Ex, 3)
two_pol = allocated(wf%Ey)
if (size(beam%slice) /= nslice) then
  call out_io (s_error$, r_name, 'BEAM HAS \i0\ SLICES BUT THE FIELD RECORD HAS \i0\ .', &
                                 i_array = [size(beam%slice), nslice])
  return
endif

p0_mc = fel_p0_mc(beam)
gamma0b = sqrt(p0_mc**2 + 1)
inv_beta0 = gamma0b / p0_mc
ks = twopi / wf%wavelength
dsub = ustate%dsub
if (abs(ustate%nsub * dsub - dz_record) > 1e-12_rp * dz_record) then
  call out_io (s_error$, r_name, 'RECORD STEP DOES NOT MATCH THE SEGMENT SETUP.')
  return
endif

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
dgrid = wf%dx
call fel_field_kernel_init (ngrid, dgrid, ks, dsub)

! The twin's own field. Its propagator is this mode's, one substep and not one record
! step, which is the whole reason the averaged instrument's field row does not carry over
! here: sixty of these run where that one runs once, so the roundings are sixty deep.
! Both caches are filled from this serial context, the single-precision planner not being
! thread safe where its executor is.

if (fp32%on) then
  call fel_fp32_field_prep (fp32, ngrid, ngrid)
  call fel_fp32_fft_plan_threads (ngrid, err)
  if (err) then
    err_flag = .true.
    return
  endif
  exp_k2_p => fel_field_kernel_exp_k2 (ks)
  if (associated(exp_k2_p)) call fel_fp32_kernel_cache (fp32, exp_k2_p, ngrid)
endif

! The averaged ponderomotive rate: phi0 advances exactly as the averaged mode's, so
! theta-derived diagnostics and the downstream bookkeeping see one convention. The
! optical carrier is Psi = phi0 - ku*s_local (the ku part of phi0 undone: it belongs
! to the undulator, whose cos(ku s) the particles now actually ride).

phi0_rate_avg = fel_phi0_rate(ks, und%ku, p0_mc)

! The unaveraged source scale: the averaged scl_w with the coupling factor removed
! and the 1/2 period-average projection undone (fel-physics.md sec-unaveraged).

scl_u = (mu_0_vac * c_light) * c_light * dsub / (2 * dgrid * dgrid * beam%slice_spacing)

! Entry handoff: the ramp makes a = 0 here, so the stored (averaged-convention) px is
! the kinetic momentum. The flag records that px now carries the quiver.

if (first) then
  call fel_assert_averaged_chart (beam, r_name // ' (segment entry)', err)
  if (err) return
  beam%quiver_in_px = .true.
  ustate%active = .true.
  call unavg_ramp_phase_jump ()
endif

! The device path (fel_device_mod). The beam and the field stay resident for the whole
! segment and one record step is one command buffer of nsub substeps. Residency starts
! here rather than at the walk's element entry because the entry handoff above moves z,
! and the chart conversion has to see the moved value.

if (dev%on .and. .not. fp32%on) then
  call unavg_device_step (err)
  if (err) return
  err_flag = .false.
  return
endif

allocate (dE_slice(nslice), dU_sp_slice(nslice), dEturn_slice(nslice))
dE_slice = 0
dU_sp_slice = 0
dEturn_slice = 0
any_err = .false.

! The device in the instrument's twin role. It advances the same record step from the
! state this one received, encoded before the FP64 loop below so the two overlap and so
! its field image rounds from the pre-step record. The FP64 path is untouched, which the
! harness asserts the way it does for the CPU twin.

dev_twin = dev%on .and. fp32%on
if (dev_twin) then
  allocate (dev_s0(6, maxval(beam%slice(:)%n), nslice))
  do is = 1, nslice
    n_s0 = beam%slice(is)%n
    dev_s0(1, 1:n_s0, is) = beam%slice(is)%x(1:n_s0)
    dev_s0(2, 1:n_s0, is) = beam%slice(is)%px(1:n_s0)
    dev_s0(3, 1:n_s0, is) = beam%slice(is)%y(1:n_s0)
    dev_s0(4, 1:n_s0, is) = beam%slice(is)%py(1:n_s0)
    dev_s0(5, 1:n_s0, is) = beam%slice(is)%z(1:n_s0)
    dev_s0(6, 1:n_s0, is) = beam%slice(is)%pz(1:n_s0)
  enddo
  call fel_device_unavg_twin_begin (dev, beam, ff, dev_s0, ustate%nsub, fp32%mutate, err)
  if (err) return
  call unavg_dev_encode (err)
  if (err) return
endif

! Parallel over slices, the averaged step's own design: each slice
! touches only its own particle arrays and its own field slice (the beam-to-field
! mapping is a bijection), the kernel cache is read-only here (initialized serially
! above), the FFT plan cache is threadprivate, and the per-slice energy lands in
! dE_slice so the final sum is a fixed-order serial reduction. Results are
! bit-identical across thread counts, and the harness checks that.

!$OMP parallel do private(sl, ifld, ux, uy, xx, yy, tau, gam, crsource, crsource_y, s_sub, isub, ip, &
!$OMP&   cx, cy, cpx, cpy, cz, cpz, &
!$OMP&   p_mc, beta, psi_mid, u_s, jhat, wx, wy, ix, iy, on_grid, ehat, ehat_y, wphasor, dgam, cdep, cph, err, fq, kst) &
!$OMP&   reduction(.or.: any_err)
do is = 1, nslice
  sl => beam%slice(is)
  ifld = fel_field_index(slip, is, nslice)

  allocate (crsource(ngrid, ngrid))
  if (two_pol) allocate (crsource_y(ngrid, ngrid))
  allocate (ux(sl%n), uy(sl%n), xx(sl%n), yy(sl%n), tau(sl%n), gam(sl%n))
  allocate (kst(5, sl%n, 4))     ! The push's RK stage arrays, one allocation per slice.

  ! The instrument's copy of the state this step received, taken before the FP64 advance
  ! touches it. The twin runs from this afterwards and the comparison reads both.

  if (fp32%on) then
    fp32%e32(:,:,is) = cmplx(wf%Ex(:,:,ifld), kind = sp)
    allocate (cx(sl%n), cy(sl%n), cpx(sl%n), cpy(sl%n), cz(sl%n), cpz(sl%n))
    cx = sl%x(1:sl%n);    cy = sl%y(1:sl%n);    cpx = sl%px(1:sl%n)
    cpy = sl%py(1:sl%n);  cz = sl%z(1:sl%n);    cpz = sl%pz(1:sl%n)
  endif
  do ip = 1, sl%n
    p_mc = p0_mc * (1 + sl%pz(ip))
    gam(ip) = sqrt(p_mc**2 + 1)
    beta = p_mc / gam(ip)
    ux(ip) = sl%px(ip) * p0_mc
    uy(ip) = sl%py(ip) * p0_mc
    xx(ip) = sl%x(ip)
    yy(ip) = sl%y(ip)
    tau(ip) = -sl%z(ip) / beta
  enddo

  s_sub = ustate%s

  do isub = 1, ustate%nsub

    crsource = 0
    if (two_pol) crsource_y = 0

    call unavg_field_quartet (s_sub, dsub/2, fq)
    call unavg_push_all (dsub/2, sl%n, xx, yy, ux, uy, tau, gam, fq, kst)

    ! Radiation kick + deposit at the substep midpoint, phi0 advanced to it.

    psi_mid = beam%phi0 + (s_sub - ustate%s + dsub/2) * phi0_rate_avg &
              - und%ku * (s_sub + dsub/2)

    do ip = 1, sl%n
      u_s = sqrt(gam(ip)**2 - 1 - ux(ip)**2 - uy(ip)**2)

      if (und%helical) then
        jhat = cmplx(ux(ip), -uy(ip), rp) / sqrt(2.0_rp)
      else
        jhat = cmplx(ux(ip), 0.0_rp, rp)
      endif

      call fel_grid_weights (wf, xx(ip), yy(ip), ix, iy, wx, wy, on_grid)
      if (on_grid .and. two_pol) then

        ! Two live polarizations: Component-wise duals. The instantaneous kinetic
        ! momenta u_x, u_y are real numbers, each working against and depositing into
        ! its own field component. No polarization convention enters at all (manual
        ! sec-field vector convention). The scalar branch below keeps the folded
        ! ĵ-convention verbatim for single-polarization lines.

        ehat =        wf%Ex(ix,   iy,   ifld) * wx * wy
        ehat = ehat + wf%Ex(ix+1, iy,   ifld) * (1-wx) * wy
        ehat = ehat + wf%Ex(ix,   iy+1, ifld) * wx * (1-wy)
        ehat = ehat + wf%Ex(ix+1, iy+1, ifld) * (1-wx) * (1-wy)
        ehat_y =          wf%Ey(ix,   iy,   ifld) * wx * wy
        ehat_y = ehat_y + wf%Ey(ix+1, iy,   ifld) * (1-wx) * wy
        ehat_y = ehat_y + wf%Ey(ix,   iy+1, ifld) * wx * (1-wy)
        ehat_y = ehat_y + wf%Ey(ix+1, iy+1, ifld) * (1-wx) * (1-wy)

        cph = exp(cmplx(0.0_rp, psi_mid - ks*tau(ip), rp))
        dgam = -dsub * real(cmplx(0.0_rp, -1.0_rp, rp) * (ehat * ux(ip) + ehat_y * uy(ip)) * cph, rp) &
                     / (u_s * m_electron)
        dE_slice(is) = dE_slice(is) + sl%weight(ip) * dgam * m_electron
        if (fp32%on) dEturn_slice(is) = dEturn_slice(is) + sl%weight(ip) * abs(dgam) * m_electron
        gam(ip) = gam(ip) + dgam

        cdep = cmplx(0.0_rp, 1.0_rp, rp) * conjg(cph) * scl_u * sl%weight(ip) / u_s
        crsource(ix,   iy)   = crsource(ix,   iy)   + (wx * wy) * cdep * ux(ip)
        crsource(ix+1, iy)   = crsource(ix+1, iy)   + ((1-wx) * wy) * cdep * ux(ip)
        crsource(ix,   iy+1) = crsource(ix,   iy+1) + (wx * (1-wy)) * cdep * ux(ip)
        crsource(ix+1, iy+1) = crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cdep * ux(ip)
        crsource_y(ix,   iy)   = crsource_y(ix,   iy)   + (wx * wy) * cdep * uy(ip)
        crsource_y(ix+1, iy)   = crsource_y(ix+1, iy)   + ((1-wx) * wy) * cdep * uy(ip)
        crsource_y(ix,   iy+1) = crsource_y(ix,   iy+1) + (wx * (1-wy)) * cdep * uy(ip)
        crsource_y(ix+1, iy+1) = crsource_y(ix+1, iy+1) + ((1-wx) * (1-wy)) * cdep * uy(ip)

      elseif (on_grid) then
        ehat =        wf%Ex(ix,   iy,   ifld) * wx * wy
        ehat = ehat + wf%Ex(ix+1, iy,   ifld) * (1-wx) * wy
        ehat = ehat + wf%Ex(ix,   iy+1, ifld) * wx * (1-wy)
        ehat = ehat + wf%Ex(ix+1, iy+1, ifld) * (1-wx) * (1-wy)

        wphasor = cmplx(0.0_rp, -1.0_rp, rp) * ehat * exp(cmplx(0.0_rp, psi_mid - ks*tau(ip), rp))
        dgam = -dsub * real(wphasor * conjg(jhat), rp) / (u_s * m_electron)
        dE_slice(is) = dE_slice(is) + sl%weight(ip) * dgam * m_electron
        if (fp32%on) dEturn_slice(is) = dEturn_slice(is) + sl%weight(ip) * abs(dgam) * m_electron
        gam(ip) = gam(ip) + dgam

        ! /u_s, not Genesis's averaged /gamma: the source and the E.v force are exact
        ! duals of one wave equation, and using the same u_s the kick used makes the
        ! per-substep energy exchange cancel identically (the diffraction between
        ! substeps is unitary), leaving only physical spontaneous emission in the
        ! ledger. The period-averaged limit shifts by beta_par ~ 5e-9 -- five orders
        ! below the fc checks. A merit choice, not a transcription (sec-unaveraged).
        cdep = cmplx(0.0_rp, 1.0_rp, rp) * exp(cmplx(0.0_rp, -(psi_mid - ks*tau(ip)), rp)) &
               * jhat * scl_u * sl%weight(ip) / u_s
        crsource(ix,   iy)   = crsource(ix,   iy)   + (wx * wy) * cdep
        crsource(ix+1, iy)   = crsource(ix+1, iy)   + ((1-wx) * wy) * cdep
        crsource(ix,   iy+1) = crsource(ix,   iy+1) + (wx * (1-wy)) * cdep
        crsource(ix+1, iy+1) = crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cdep
      endif
    enddo

    call unavg_field_quartet (s_sub + dsub/2, dsub/2, fq)
    call unavg_push_all (dsub/2, sl%n, xx, yy, ux, uy, tau, gam, fq, kst)

    ! The source lands on the record the kick read, and the pair diffracts together:
    ! E' = D (E + 2 src). The kick charged the beam the cross term 2 Re<E, dE> against
    ! this E, and with D unitary the record's energy moves by exactly that plus the
    ! deposit's own 4|src|^2, so the ledger closes to rounding. Diffracting first and
    ! adding after, E' = D E + 2 src, charges the beam against E while the record gains
    ! the cross term against D E, and the difference, 4 Re<(D - I) E, dE>, is first
    ! order in the substep's diffraction phase: 1.2 percent of the turnover at 20
    ! substeps a period on a seed of 100 um waist, halving with the substep (FINDINGS
    ! 7.86). The stored record is then the substep's midpoint field, which is the same
    ! Strang split read at the other half step, and the device's solve does the same.

    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) + 2 * crsource
    if (two_pol) wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) + 2 * crsource_y
    call fel_field_diffract (wf, ifld, dsub, err)
    any_err = any_err .or. err

    ! The same step on the twin's record, in single precision, from the source the FP64
    ! side just built: the row then prices the transform pair, the rounded propagator and
    ! the accumulation over the substeps, and not the deposit that fed them.

    if (fp32%on .and. allocated(fp32%k32)) then
      fp32%e32(:,:,is) = fp32%e32(:,:,is) + 2 * cmplx(crsource, kind = sp)
      call fft32_solve (fp32%e32(:,:,is), fp32%k32, ngrid)
    endif

    ! The deposit's own energy |dE|^2 = 4|src|^2: the one term of the field-energy
    ! increment the kick/deposit duality does not charge to the beam (the beam pays the
    ! cross term 2 Re<E, dE> exactly, manual eq:ledger). Physically this is the
    ! spontaneous emission of the substep. Numerically it is banked here so the
    ! time-dependent ledger closes exactly: E_beam + U_window + U_escaped - U_spont.

    dU_sp_slice(is) = dU_sp_slice(is) + 4 * sum(real(crsource, rp)**2 + aimag(crsource)**2) &
                        * dgrid**2 / (2 * (mu_0_vac * c_light)) * (beam%slice_spacing / c_light)
    if (two_pol) then
      dU_sp_slice(is) = dU_sp_slice(is) + 4 * sum(real(crsource_y, rp)**2 + aimag(crsource_y)**2) &
                          * dgrid**2 / (2 * (mu_0_vac * c_light)) * (beam%slice_spacing / c_light)
    endif

    s_sub = s_sub + dsub
  enddo

  ! Back to the stored chart. px carries the quiver mid-segment (the flag says so).
  ! z = -beta*tau with the full-momentum beta, the chart's own convention.

  do ip = 1, sl%n
    p_mc = sqrt(gam(ip)**2 - 1)
    sl%pz(ip) = (p_mc - p0_mc) / p0_mc
    sl%px(ip) = ux(ip) / p0_mc
    sl%py(ip) = uy(ip) / p0_mc
    sl%x(ip) = xx(ip)
    sl%y(ip) = yy(ip)
    beta = p_mc / gam(ip)
    sl%z(ip) = -beta * tau(ip)
  enddo

  deallocate (ux, uy, xx, yy, tau, gam, crsource, kst)
  if (two_pol) deallocate (crsource_y)
  if (fp32%on .and. .not. dev_twin) then
    call unavg_twin_slice (is, sl%n, ifld, cx, cy, cpx, cpy, cz, cpz)
    deallocate (cx, cy, cpx, cpy, cz, cpz)
  elseif (fp32%on) then
    deallocate (cx, cy, cpx, cpy, cz, cpz)
  endif

enddo
!$OMP end parallel do
if (any_err) return

dE_beam = sum(dE_slice)     ! Fixed-order serial sums: thread-count independent.
dU_spont = sum(dU_sp_slice)

beam%phi0 = beam%phi0 + dz_record * phi0_rate_avg
ustate%s = ustate%s + dz_record

! Exit handoff: the ramp has closed (a = 0), kinetic equals canonical again.

! The device twin's rows, read before the exit handoff: the ramp's phase jump moves z
! on the host and the twin advanced the record step, not the handoff, so a comparison
! taken after it would price the jump instead of the arithmetic. The CPU twin runs
! inside the loop above and sits before the handoff for the same reason.

if (dev_twin) call fel_device_unavg_twin_rows (dev, fp32, beam, ff, ks, dE_slice, dEturn_slice)

if (last) then
  beam%quiver_in_px = .false.
  ustate%active = .false.
  call unavg_ramp_phase_jump ()
endif

! The instrument's serial epilogue, which reduces this step's per-slice rows and writes
! the stream. Outside the parallel region by contract, as the averaged path's is.

if (fp32%on) then
  call fel_fp32_step_close (fp32, err)
  if (err) then
    err_flag = .true.
    return
  endif
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

! One record step on the device. The quantities that do not depend on the particle are
! built here in FP64 and uploaded once a substep rather than once a particle, which is
! the hoist FINDINGS 7.37 records, and everything else is the kernels'. What the device
! does not carry the caller has refused at setup, so nothing here falls back.

subroutine unavg_device_step (uerr)

real(rp) du_now
logical uerr

uerr = .true.

if (first) then
  dev%unavg = .true.
  dev%spont_prev = 0
  call fel_device_element_begin (dev, beam, ff)
  call fel_device_unavg_begin (dev, ustate%nsub, err)
  if (err) return
endif

call unavg_dev_encode (err)
if (err) return

! The ledger's two device-side terms. The command buffer has to drain before the next
! record step's stage factors can be written anyway, so reading them here costs the
! step nothing it was not already paying.

call fel_device_unavg_ledger (dev, beam, dE_beam, du_now)
dU_spont = du_now - dev%spont_prev
dev%spont_prev = du_now

beam%phi0 = beam%phi0 + dz_record * phi0_rate_avg
ustate%s = ustate%s + dz_record

! Exit handoff. The state comes back to the host first, since the ramp phase jump and
! the chart flag are the host arrays' business and the walk's slippage and stats read
! them from here on.

if (last) then
  call fel_device_readback (dev, beam, ff)
  call fel_device_release (dev)
  dev%unavg = .false.
  if (dev%dep_breach) return
  beam%quiver_in_px = .false.
  ustate%active = .false.
  call unavg_ramp_phase_jump ()
endif

uerr = .false.

end subroutine unavg_device_step

! One record step's encode: the propagator, the per-substep stage factors and carrier
! rotators, the deposit's bound, and the command buffer itself. The production path and
! the instrument's twin role both go through here, so the arithmetic under test is one
! arithmetic.

subroutine unavg_dev_encode (uerr)

type (fel_device_unavg_par_struct) upar
real(rp), allocatable :: fqs(:,:,:,:)
complex(rp), allocatable :: cbase(:,:)
real(rp) s0, psi_m, psi_s, s_bound, q_tot
integer j, is_d
logical uerr

uerr = .true.

! A readback found a particle whose |j|/u_s broke the deposit's bound, so the bound no
! longer holds and nothing further is deposited. The walk stops the run at the comb
! position that read it; this is the second reading, as fill_device_par's is.

if (dev%dep_breach) return

! The propagator at the substep, the entry fel_field_diffract reads. The cache was
! filled for this ngrid, spacing and substep above, so a miss is a bug.

exp_k2_p => fel_field_kernel_exp_k2 (ks)
if (.not. associated(exp_k2_p)) then
  call out_io (s_error$, r_name, 'THE UNAVERAGED SUBSTEP PROPAGATOR IS NOT IN THE KERNEL CACHE.', &
               'PLEASE REPORT THIS!')
  return
endif
call fel_device_set_kernel (dev, 1, exp_k2_p)

! The stage field factors of both half pushes, and the optical carrier's base rotator
! per substep per slice. psi_mid reaches some 1700 radians over a segment, which no
! float carries, so it is reduced modulo 2 pi here and crosses the seam as a rotator.
! The lag reference is the slice's own mean, held in FP64 by the device state.

allocate (fqs(4, 4, 2, ustate%nsub), cbase(nslice, ustate%nsub))
do j = 1, ustate%nsub
  s0 = ustate%s + (j - 1) * dsub
  call unavg_field_quartet (s0, dsub/2, fq)
  fqs(:,:,1,j) = fq
  call unavg_field_quartet (s0 + dsub/2, dsub/2, fq)
  fqs(:,:,2,j) = fq
  psi_m = beam%phi0 + (s0 - ustate%s + dsub/2) * phi0_rate_avg - und%ku * (s0 + dsub/2)
  do is_d = 1, nslice
    psi_s = modulo(psi_m - ks * dev%z_ref(is_d), twopi)
    cbase(is_d, j) = cmplx(cos(psi_s), sin(psi_s), rp)
  enddo
enddo

upar%dsub = dsub
upar%ks = ks
upar%ku = und%ku
upar%aw = und%aw
upar%gam0 = gamma0b
upar%beta0 = p0_mc / gamma0b
upar%g0inv2 = 1 / gamma0b**2
upar%cos_t = und%cos_t
upar%sin_t = und%sin_t
upar%gridmax = (ngrid - 1) * dgrid / 2
upar%dgrid = dgrid
upar%scl_u = scl_u
upar%m_electron = m_electron
upar%spont_fac = 4 * dgrid**2 / (2 * (mu_0_vac * c_light)) * (beam%slice_spacing / c_light)
upar%nsub = ustate%nsub
upar%first = slip%first
upar%helical = merge(1, 0, und%helical)
upar%pad = 0

! The deposit's bound. A contribution here carries |j|/u_s where the averaged deposit
! carries a roll-off over gamma, so the bound is the run's whole charge times the
! largest ratio the kernel is allowed to meet. The quiver sets that ratio: |j| reaches
! sqrt(2) aw at the peak of a planar wiggle and less on a helical one, the betatron
! momentum adds a hundredth of that on the checked decks, and u_s falls no lower than
! the gamma floor the device measured at setup. One rest momentum of headroom on the
! numerator covers a transverse momentum no FEL beam has. The kernel checks the ratio
! on every particle and refuses the run through the same fault the averaged deposit's
! gamma floor uses, so the bound is checked rather than assumed.

q_tot = dev%dep_qbound * dev%dep_gam_floor
upar%dep_u_bound = (sqrt(2.0_rp) * und%aw + 1) / dev%dep_gam_floor
s_bound = scl_u * q_tot * upar%dep_u_bound
call fel_device_dep_scale (dev, s_bound, 'the run''s charge, the substep''s source ' // &
        'scale and the largest quiver-to-longitudinal momentum ratio the kernel allows', &
        upar%dep_scale, err)
if (err) return

upar%mutate = merge(1, 0, fp32%mutate)
call fel_device_unavg_step (dev, upar, fqs, cbase, err)
if (err) return

uerr = .false.

end subroutine unavg_dev_encode

! The ramp's slippage compensation, applied as one discrete jump per segment end,
! at the handoffs where the envelope is exactly zero and nothing couples. The
! element's contract is L meters of full-strength undulator. The entry/exit ramps
! are a numerical device (an adiabatic switch-on), and an electron under a ramped
! quiver <u_perp^2> = g^2 aw^2 lags the wave less than the contracted hard-edge
! element would have it, by dtau = (gamma aw^2 / 2 u_s^3) INT (1-g^2) ds, which
! is ks*dtau ~ 2.6 rad of optical phase per end at the benchmark parameters
! (sin^2 envelope: INT (1-g^2) = (5/8) l_ramp per end). Uncompensated, every
! segment after the first receives a pre-bunched beam with its bunching-to-field
! phase rotated by that much. The first segment is immune (nothing is bunched
! yet), and that is exactly how it was caught: per-segment ln-power deviations vs
! the averaged mode of {0.0000, +0.08, +0.005, +0.02, +0.01, +0.13}. The jump must
! be discrete and at the ends: compensating continuously inside the ramp detunes
! the live interaction where the coupling is already substantial (measured -0.9%
! gain on the first segment, doubling with ramp length). In hardware terms this is
! the phase shifter that makes a tapered-end segment equivalent to its ideal
! hard-edged length. In the stored chart z = -beta*tau, so tau += dtau is
! z -= aw^2 INT(1-g^2) / (2 p^2), per particle with its own momentum.

subroutine unavg_ramp_phase_jump ()

real(rp) ramp_int
integer is_j, ip_j

if (ustate%l_ramp <= 0) return
ramp_int = 0.625_rp * ustate%l_ramp * und%aw**2 / 2

!$OMP parallel do private(ip_j)
do is_j = 1, size(beam%slice)
  do ip_j = 1, beam%slice(is_j)%n
    beam%slice(is_j)%z(ip_j) = beam%slice(is_j)%z(ip_j) - &
                     ramp_int / (p0_mc * (1 + beam%slice(is_j)%pz(ip_j)))**2
  enddo
enddo
!$OMP end parallel do

end subroutine unavg_ramp_phase_jump

! The s-dependent field factors at the four RK stage positions of a push over h from
! s0: the envelope and its slope, and the undulator phase's cosine and sine. None of
! them depends on the particle, so one call here replaces four per particle. The stage
! positions are the expressions the push evaluated inline before the hoist, in the
! same order, so every value is the same bits it was. Stages two and three share a
! position and the third is a copy of the second, which says so.

subroutine unavg_field_quartet (s0, h, fq)

real(rp) s0, h, fq(4,4)
real(rp) s_st(4), g, gp
integer j

s_st = [s0, s0 + h/2, s0 + h/2, s0 + h]

do j = 1, 4
  if (j == 3) then
    fq(:,3) = fq(:,2)
    cycle
  endif
  g = fel_unavg_envelope(ustate, s_st(j), gp)
  fq(1,j) = g
  fq(2,j) = gp
  fq(3,j) = cos(und%ku * s_st(j))
  fq(4,j) = sin(und%ku * s_st(j))
enddo

end subroutine unavg_field_quartet

! One RK4 magnetic push of a whole slice over step h, stages outermost: each RK stage
! is one loop over the particles into its stage array, and the combination is a last
! loop. Per particle this is the same operations on the same values in the same order
! as the particle-outermost form it replaces (nothing couples particles inside a push,
! and the stage field factors are per stage already), so the results are byte-identical
! and the identity checks hold it there. Stages outermost is also the shape a GPU
! kernel transcribes: array work per stage, no call tree per particle.
!
! All per-particle state passes by argument: host-associated variables privatized by
! the caller's OMP region are not redirected inside called procedures, so nothing
! mutable may be host-associated here (und/ustate/inv_beta0 are read-only shared).
! gamma is untouched: B does no work, exactly. kst is caller scratch, sized (5, n, 4),
! so the allocation is paid once per slice rather than once per push.

subroutine unavg_push_all (h, n, xx, yy, ux, uy, tau, gam, fq, kst)

integer n, ip
real(rp) h, xx(n), yy(n), ux(n), uy(n), tau(n), gam(n)
real(rp) fq(4,4), kst(5,n,4)
real(rp) yt(5)

do ip = 1, n
  call unavg_ode ([xx(ip), yy(ip), ux(ip), uy(ip), tau(ip)], fq(:,1), gam(ip), kst(:,ip,1))
enddo

do ip = 1, n
  yt = [xx(ip), yy(ip), ux(ip), uy(ip), tau(ip)] + (h/2) * kst(:,ip,1)
  call unavg_ode (yt, fq(:,2), gam(ip), kst(:,ip,2))
enddo

do ip = 1, n
  yt = [xx(ip), yy(ip), ux(ip), uy(ip), tau(ip)] + (h/2) * kst(:,ip,2)
  call unavg_ode (yt, fq(:,3), gam(ip), kst(:,ip,3))
enddo

do ip = 1, n
  yt = [xx(ip), yy(ip), ux(ip), uy(ip), tau(ip)] + h * kst(:,ip,3)
  call unavg_ode (yt, fq(:,4), gam(ip), kst(:,ip,4))
enddo

! The combination, elementwise exactly as the per-particle form wrote it.

do ip = 1, n
  xx(ip)  = xx(ip)  + (h/6) * (kst(1,ip,1) + 2*kst(1,ip,2) + 2*kst(1,ip,3) + kst(1,ip,4))
  yy(ip)  = yy(ip)  + (h/6) * (kst(2,ip,1) + 2*kst(2,ip,2) + 2*kst(2,ip,3) + kst(2,ip,4))
  ux(ip)  = ux(ip)  + (h/6) * (kst(3,ip,1) + 2*kst(3,ip,2) + 2*kst(3,ip,3) + kst(3,ip,4))
  uy(ip)  = uy(ip)  + (h/6) * (kst(4,ip,1) + 2*kst(4,ip,2) + 2*kst(4,ip,3) + kst(4,ip,4))
  tau(ip) = tau(ip) + (h/6) * (kst(5,ip,1) + 2*kst(5,ip,2) + 2*kst(5,ip,3) + kst(5,ip,4))
enddo

end subroutine unavg_push_all

! The exact z-ODEs of ballistic motion in the magnetostatic field:
!   dx/ds = u_x/u_s, du_x/ds = b_y - u_y b_z/u_s, du_y/ds = -b_x + u_x b_z/u_s,
!   dtau/ds = gamma/u_s - 1/beta0.

subroutine unavg_ode (y, fq, gamma, dyds)

real(rp) y(5), fq(4), gamma, dyds(5)
real(rp) bx, by, bz, us_l

call fel_unavg_bfield (und, y(1), y(2), fq(1), fq(2), fq(3), fq(4), bx, by, bz)
us_l = sqrt(gamma**2 - 1 - y(3)**2 - y(4)**2)

dyds(1) = y(3) / us_l
dyds(2) = y(4) / us_l
dyds(3) = by - y(4) * bz / us_l
dyds(4) = -bx + y(3) * bz / us_l

dyds(5) = gamma / us_l - inv_beta0

end subroutine unavg_ode

!------------------------------------------------------------------------------
! contains
!+
! Subroutine unavg_twin_slice (is, n, x0, y0, px0, py0, z0, pz0)
!
! Routine to advance one slice's single-precision twin over this record step's
! substeps and fill the slice's divergence row. The FP64 arrays are the state as this
! step received it, and beam%slice(is) already carries the FP64 result, so the two are
! compared where the record step ends. That granularity is the choice: a substep is
! internal to the integrator, where a record step is the point the FP64 path and any
! device port synchronize, so it measures what a port would have to match.
!
! Both sides gather the field from the FP64 record. What the twin therefore prices is
! the particle path, which is where every reformulation below lives. The field's own
! single-precision accumulation is a separate quantity and is not measured here: this
! mode diffracts nsub times a record step where the averaged mode diffracts once, so it
! pays that many roundings, and the averaged instrument's field row does not carry over.
! doc/validation.md says so where the levels are recorded.
!
! What single precision destroys in this advance, each measured on a real mid-segment
! state of the one-segment deck rather than argued:
!
!   slippage  dtau/ds = gamma/u_s - 1/beta0. Both terms are one to within 1.3e-8 and
!             the difference is 1.9e-9, which is 0.016 of the quantum of one in single
!             precision, so the naive difference is exactly zero and the slippage is
!             gone. Formed through the identity in fel_fp32_unavg_ode it comes back to
!             1.2e-7. This is the reformulation the mode lives on.
!   energy    gamma is 11358 with a quantum of 9.8e-4 against a per-substep change of
!             1.0e-6, which is a thousandth of a quantum. The working variable is
!             goff = gamma - gamma0, whose quantum is 6.0e-8 and which carries the same
!             change at 18 quanta.
!   lag       tau is the particle's lag and its per-substep change is 3.1e-12 m, which
!             the guard below measures rather than infers. Most of that is the quiver's
!             own longitudinal oscillation and not the slippage drift, so a figure taken
!             from two states many substeps apart understates it by nearly three orders.
!             Carrying the whole lag, one slice resolves the change at 4.4e5 quanta and
!             32 slices at 1.4e4, but a 351-slice window at 164 nm reaches 5.8e-5 m and
!             resolves it at 0.84 of a quantum, which is lost. The residual off a
!             per-slice FP64 reference is one slice's spread whatever the window, so it
!             does not move with the slice count. The margin is the reason the guard
!             watches this and not something else.
!   phase     Psi carried whole per particle costs 2.7e-5 rad at a segment's end and
!             grows with s. The base is one FP64 number a slice a substep, which is
!             what a device kernel would upload, and only the residual ks*dtaur is
!             single precision, at 2.4e-7 rad.
!
! u_s itself needs no reformulation, which measurement rather than expectation settled:
! sqrt(gamma^2 - 1 - ux^2 - uy^2) does lose the 1 and the ux^2 entirely in single
! precision, but they are 6.7e-9 of the result, so the naive and the reformed values
! agree to 5.0e-8. The same value serves the kick and the deposit, which is what keeps
! the pair exact energy duals.
!
! Runs inside the caller's parallel slice loop: everything written is indexed by is.
!-

subroutine unavg_twin_slice (is, n, ifl, x0, y0, px0, py0, z0, pz0)

! ifl is passed and not host-associated, and the slice is reached through beam, which
! is shared: this module's push header states the rule, that a variable the caller's
! OMP region privatizes is not redirected inside a called procedure.

integer is, n, ifl
real(rp) x0(:), y0(:), px0(:), py0(:), z0(:), pz0(:)

type (fel_fp32_unavg_slice_struct), pointer :: t
type (fel_fp32_unavg_struct) uc
real(sp) fq32(4,4), h32, psi_b32, ks32, ehat_r, ehat_i, cs, sn, ang
real(sp) us32, dg32, jr, ji, wr, wi
real(rp) fqd(4,4), s_t, psi_b, tau_r, gam_l, beta_l, p_mc_l
real(rp) dstat(fel_fp32_nq$), sc(fel_fp32_nq$), wsum, tau64
real(rp), allocatable :: dtau_ulp(:)
real(sp), allocatable :: tsave(:)
real(rp) gmed, dE32, dEturn
real(rp) p32r, p32i, p64r, p64i, th64, th32
integer ip, jsub, ixl, iyl
logical og
real(rp) wxd, wyd
real(sp) wxl, wyl

!

t => fp32%un32(is)
if (.not. allocated(t%xx)) then
  allocate (t%xx(n), t%yy(n), t%ux(n), t%uy(n), t%dtaur(n), t%goff(n))
elseif (size(t%xx) /= n) then
  deallocate (t%xx, t%yy, t%ux, t%uy, t%dtaur, t%goff)
  allocate (t%xx(n), t%yy(n), t%ux(n), t%uy(n), t%dtaur(n), t%goff(n))
endif

! The constants the twin's arithmetic reads, rounded once. beta0 and 1/gamma0^2 are the
! slippage identity's two halves, and they are formed in FP64 here because they are
! per-run scalars: a device would upload them the same way.

uc%aw = real(und%aw, sp);        uc%ku = real(und%ku, sp)
uc%cos_t = real(und%cos_t, sp);  uc%sin_t = real(und%sin_t, sp)
uc%gam0 = real(gamma0b, sp)
uc%beta0 = real(p0_mc / gamma0b, sp)
uc%g0inv2 = real(1 / gamma0b**2, sp)
uc%helical = und%helical

! The lag reference is this slice's own FP64 mean, so the residual the twin carries is a
! slice's spread and not a window's offset.

tau_r = 0
do ip = 1, n
  p_mc_l = p0_mc * (1 + pz0(ip))
  gam_l = sqrt(p_mc_l**2 + 1)
  tau_r = tau_r - z0(ip) / (p_mc_l / gam_l)
enddo
tau_r = tau_r / max(1, n)
t%tau_ref = tau_r

do ip = 1, n
  p_mc_l = p0_mc * (1 + pz0(ip))
  gam_l = sqrt(p_mc_l**2 + 1)
  beta_l = p_mc_l / gam_l
  t%xx(ip) = real(x0(ip), sp)
  t%yy(ip) = real(y0(ip), sp)
  t%ux(ip) = real(px0(ip) * p0_mc, sp)
  t%uy(ip) = real(py0(ip) * p0_mc, sp)
  t%goff(ip) = real(gam_l - gamma0b, sp)
  t%dtaur(ip) = real(-z0(ip) / beta_l - tau_r, sp)
enddo

ks32 = real(ks, sp)
h32 = real(dsub, sp)
s_t = ustate%s

! The guard. The lag residual is the one single-precision quantity here whose per-substep
! change can fall under its own quantum, and whether it does depends on the window and
! not on the physics: the residual is a slice's own spread, where the whole lag grows with
! the slice count. The statistic is that change in units of the residual's spacing, summed
! over the substeps here and reduced to a median at the end, which is the shape and the
! floor the averaged instrument's guard already uses.

allocate (dtau_ulp(n), tsave(n))
dtau_ulp = 0
dE32 = 0
dEturn = 0

do jsub = 1, ustate%nsub

  tsave = t%dtaur
  call unavg_field_quartet (s_t, dsub/2, fqd)
  fq32 = real(fqd, sp)
  call fel_fp32_unavg_push (h32/2, n, t, uc, fq32)

  ! The kick, at the substep midpoint. The base phase is one FP64 number for the whole
  ! slice, so the per-particle single-precision angle is the residual alone.

  psi_b = beam%phi0 + (s_t - ustate%s + dsub/2) * phi0_rate_avg - und%ku * (s_t + dsub/2) &
          - ks * tau_r
  psi_b32 = real(modulo(psi_b, twopi), sp)

  do ip = 1, n
    call fel_grid_weights (wf, real(t%xx(ip), rp), real(t%yy(ip), rp), ixl, iyl, wxd, wyd, og)
    if (.not. og) cycle
    wxl = real(wxd, sp);  wyl = real(wyd, sp)
    ehat_r = real(real(wf%Ex(ixl, iyl, ifl), rp) * wxd * wyd &
                + real(wf%Ex(ixl+1, iyl, ifl), rp) * (1-wxd) * wyd &
                + real(wf%Ex(ixl, iyl+1, ifl), rp) * wxd * (1-wyd) &
                + real(wf%Ex(ixl+1, iyl+1, ifl), rp) * (1-wxd) * (1-wyd), sp)
    ehat_i = real(aimag(wf%Ex(ixl, iyl, ifl)) * wxd * wyd &
                + aimag(wf%Ex(ixl+1, iyl, ifl)) * (1-wxd) * wyd &
                + aimag(wf%Ex(ixl, iyl+1, ifl)) * wxd * (1-wyd) &
                + aimag(wf%Ex(ixl+1, iyl+1, ifl)) * (1-wxd) * (1-wyd), sp)
    ang = psi_b32 - ks32 * t%dtaur(ip)
    cs = cos(ang);  sn = sin(ang)
    us32 = sqrt((uc%gam0 + t%goff(ip))**2 - 1.0_sp - t%ux(ip)**2 - t%uy(ip)**2)
    if (uc%helical) then
      jr = t%ux(ip) / sqrt(2.0_sp);  ji = -t%uy(ip) / sqrt(2.0_sp)
    else
      jr = t%ux(ip);  ji = 0.0_sp
    endif
    ! W = -i Ehat e^{i ang}; dgamma = -dsub Re[W conj(j)] / (u_s m_e).
    wr =  ehat_r * sn + ehat_i * cs
    wi = -ehat_r * cs + ehat_i * sn
    dg32 = -h32 * (wr * jr + wi * ji) / (us32 * real(m_electron, sp))
    dE32 = dE32 + beam%slice(is)%weight(ip) * real(dg32, rp) * m_electron
    dEturn = dEturn + beam%slice(is)%weight(ip) * abs(real(dg32, rp)) * m_electron
    t%goff(ip) = t%goff(ip) + dg32
  enddo

  call unavg_field_quartet (s_t + dsub/2, dsub/2, fqd)
  fq32 = real(fqd, sp)
  call fel_fp32_unavg_push (h32/2, n, t, uc, fq32)
  do ip = 1, n
    if (spacing(t%dtaur(ip)) > 0) dtau_ulp(ip) = dtau_ulp(ip) + &
              abs(real(t%dtaur(ip) - tsave(ip), rp)) / real(spacing(t%dtaur(ip)), rp)
  enddo
  s_t = s_t + dsub
enddo

! The rows. Every one is a worst per-particle difference, which sees a common offset:
! ux and uy hold the quiver, which is an offset a slice shares and not a spread, so a
! statistic blind to an offset would report health on a wrong chart.

dstat = 0
wsum = 0
p32r = 0;  p32i = 0;  p64r = 0;  p64i = 0
do ip = 1, n
  p_mc_l = p0_mc * (1 + beam%slice(is)%pz(ip))
  gam_l = sqrt(p_mc_l**2 + 1)
  beta_l = p_mc_l / gam_l
  tau64 = -beam%slice(is)%z(ip) / beta_l
  dstat(1) = max(dstat(1), abs(real(t%xx(ip), rp) - beam%slice(is)%x(ip)))
  dstat(2) = max(dstat(2), abs(real(t%ux(ip), rp) / p0_mc - beam%slice(is)%px(ip)))
  dstat(3) = max(dstat(3), abs(real(t%yy(ip), rp) - beam%slice(is)%y(ip)))
  dstat(4) = max(dstat(4), abs(real(t%uy(ip), rp) / p0_mc - beam%slice(is)%py(ip)))
  dstat(5) = max(dstat(5), abs((real(t%goff(ip), rp) + gamma0b - gam_l) / p0_mc))
  th32 = -ks * (real(t%dtaur(ip), rp) + tau_r)
  th64 = -ks * tau64
  dstat(6) = max(dstat(6), abs(modulo(th32 - th64 + pi, twopi) - pi))
  wsum = wsum + beam%slice(is)%weight(ip)
  p32r = p32r + beam%slice(is)%weight(ip) * cos(th32)
  p32i = p32i - beam%slice(is)%weight(ip) * sin(th32)
  p64r = p64r + beam%slice(is)%weight(ip) * cos(th64)
  p64i = p64i - beam%slice(is)%weight(ip) * sin(th64)
enddo

sc(1) = maxval(abs(beam%slice(is)%x(1:n))) + 1e-30_rp
sc(2) = maxval(abs(beam%slice(is)%px(1:n))) + 1e-30_rp
sc(3) = maxval(abs(beam%slice(is)%y(1:n))) + 1e-30_rp
sc(4) = maxval(abs(beam%slice(is)%py(1:n))) + 1e-30_rp
sc(5) = maxval(abs(beam%slice(is)%pz(1:n))) + 1e-30_rp
sc(6) = 1
sc(7) = wsum + 1e-30_rp

fp32%div_slice(1:6, is) = dstat(1:6) / sc(1:6)
fp32%div_slice(7, is) = sqrt((p32r - p64r)**2 + (p32i - p64i)**2) / sc(7)
! The energy the twin's kicks took from the beam, against the FP64 step's own. This is
! the ledger's own term and not a coordinate, and it is the row a single-precision path
! that quietly stopped conserving would move first, the kick and the deposit being exact
! duals only while they share operands.
!
! The scale is the energy the step actually moved and not the energy it netted. Over a
! record step the gains and the losses very nearly cancel, so the net is a small
! difference of large exchanges and a relative error against it means nothing: on the
! benchmark segment that ratio reaches 3, which is a statement about the cancellation
! rather than about single precision. The turnover is the sum of the absolute exchanges,
! which is what the ledger's own check normalizes by for the same reason, and what the
! averaged instrument's phasor row does when it scales to charge rather than to a
! noise-level sum.

fp32%div_slice(8, is) = abs(dE32 - dE_slice(is)) / (dEturn + 1e-30_rp)
! The field row, the twin's record against the FP64 one, normalized by the FP64 field's
! own norm as the averaged instrument's is.

if (allocated(fp32%k32)) then
  fp32%div_slice(9, is) = sqrt(sum(abs(cmplx(fp32%e32(:,:,is), kind = rp) - wf%Ex(:,:,ifl))**2)) / &
                          (sqrt(sum(abs(wf%Ex(:,:,ifl))**2)) + 1e-30_rp)
else
  fp32%div_slice(9, is) = 0
endif
fp32%bmag64(is) = sqrt(p64r**2 + p64i**2) / sc(7)
fp32%bmag32(is) = sqrt(p32r**2 + p32i**2) / sc(7)
call fel_fp32_median (dtau_ulp, n, gmed)
fp32%ulp_slice(is) = gmed / max(1, ustate%nsub)
deallocate (dtau_ulp, tsave)

end subroutine unavg_twin_slice

end subroutine fel_unavg_step

end module fel_unaveraged_mod
