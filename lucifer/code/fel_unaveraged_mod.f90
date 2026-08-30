!+
! Module fel_unaveraged_mod
!
! The unaveraged mode (fel-physics.md sec-unaveraged):
! particles integrated through the undulator's real field (the full Newton-Lorentz
! quiver, no period averaging) with the radiation field as a co-evolving kick. There
! is no resonance approximation. What the cost per step buys is physics the averaged
! map cannot reach: the full quiver dynamics, the energy accounting the beam actually
! pays, polarization-agnostic coupling and arbitrary harmonic content. The mode is
! also the referee for the averaged path, since the coupling factor fc, the harmonic
! content and the entry/exit behavior of the averaged mode become measured outputs
! here instead of assumed inputs. MINERVA (Freund and van der Slot) is the production
! existence proof for this physics. This mode differs by keeping the grid field and
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
!      REMOVED and the actual quiver current in its place:
!        src += i e^{-i Psi} * j * (Z0 c dz /(2 dgrid^2 Ds)) * w/u_s
!      followed by the shared pure diffraction (fel_field_diffract) and the +2*src
!      convention. The /u_s (where the averaged solver has Genesis's /gamma) makes the
!      kick/deposit pair EXACT energy duals per substep (same operands, same bilinear
!      weights, unitary diffraction between), so the ledger closes to the physical
!      spontaneous-emission term and rounding, by construction. Period-averaging the
!      pair reproduces the averaged mode's fc to O(1-beta_par) ~ 5e-9 (the JJ factor
!      emerges from the figure-8).
!   3. half magnetic push.
!
! Units: E in V/m (wavefront convention), m_electron in eV, b in 1/m. The physical
! field of the scalar envelope: planar E_x = Re[-i Ehat e^{i Psi}]; helical
! (E_x, E_y) = (Re[-i Ehat e^{i Psi}], -Re[Ehat e^{i Psi}])/sqrt(2). Both give
! intensity |Ehat|^2/(2 Z0), so the power diagnostic is mode-independent.
!
! The magnetic push is classical RK4 on the exact z-ODEs in kinetic variables,
! chosen on merit, not because MINERVA uses it: the currency here is short-probe
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
  real(rp) :: l_ramp = 0           ! sin^2 ramp length at EACH end [m].
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
! requested steps per period (MINERVA's envelope: 10 floor, 20-30 recommended), ramps
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
  call out_io (s_error$, r_name, 'fel_steps_per_period BELOW MINERVA''S FLOOR OF 10: \i0\ ', &
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
! flat 1, sin^2 down over the last l_ramp. Amplitude AND slope are continuous (the
! slope gp feeds the ramp-induced field terms in fel_unavg_bfield). l_ramp = 0 is
! the hard-edge MUTATION configuration. The handoff check exists to catch it.
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
! kx = 0, ky = ku^2; helical kx = ky = ku^2/2) -- checked in sec:unaveraged.
!
! Input:
!   und        -- fel_und_struct: Undulator parameters (helicity, tilt frame).
!   ustate     -- fel_unavg_struct: Ramp geometry.
!   x, y, s    -- real(rp): Position [m].
!
! Output:
!   bx, by, bz -- real(rp): The analytic undulator field B = curl(a), with the
!                   envelope's g' terms so the ramps stay divergence-free [T].
!-

subroutine fel_unavg_bfield (und, ustate, x, y, s, bx, by, bz)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
real(rp) x, y, s, bx, by, bz
real(rp) a0, g, gp, c_u, s_u, fperp, xl, yl, bt

!

g = fel_unavg_envelope(ustate, s, gp)
c_u = cos(und%ku * s)
s_u = sin(und%ku * s)

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
!   wf        -- wavefront_struct: The fundamental field.
!   slip      -- fel_slip_struct: The rotating field record.
!   dz_record -- real(rp): The record step to advance by [m].
!   first     -- logical: True on the segment's first record step (entry handoff).
!   last      -- logical: True on the last (exit handoff and ramp phase jump).
!
! Output:
!   beam      -- fel_beam_struct: Advanced by Newton-Lorentz RK4 through the field.
!   wf        -- wavefront_struct: Sources deposited, records diffracted.
!   dE_beam   -- real(rp): The step's kick-side beam energy change [J] (ledger).
!   dU_spont  -- real(rp): The step's spontaneous source energy [J] (ledger).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_unavg_step (und, ustate, beam, wf, slip, dz_record, first, last, dE_beam, dU_spont, err_flag)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (wavefront_struct), target :: wf
type (fel_slip_struct) slip

real(rp) dz_record, dE_beam, dU_spont
logical first, last, err_flag

real(rp), allocatable :: ux(:), uy(:), xx(:), yy(:), tau(:), gam(:), dE_slice(:), dU_sp_slice(:)
complex(rp), allocatable :: crsource(:,:), crsource_y(:,:)
real(rp) p0_mc, gamma0b, inv_beta0, ks, dsub, s_sub, phi0_rate_avg, scl_u, dgrid
real(rp) u_s, wx, wy, psi_mid, dgam, p_mc, beta
complex(rp) ehat, jhat, wphasor, cdep, ehat_y, cph
logical two_pol
integer is, ip, isub, nslice, ifld, ix, iy, ngrid_arr(3), ngrid
logical on_grid, err, any_err
character(*), parameter :: r_name = 'fel_unavg_step'

!

err_flag = .true.
dE_beam = 0
dU_spont = 0

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

! The averaged ponderomotive rate: phi0 advances exactly as the averaged mode's, so
! theta-derived diagnostics and the downstream bookkeeping see one convention. The
! optical carrier is Psi = phi0 - ku*s_local (the ku part of phi0 undone: it belongs
! to the undulator, whose cos(ku s) the particles now actually ride).

phi0_rate_avg = fel_phi0_rate(ks, und%ku, p0_mc)

! The unaveraged source scale: the averaged scl_w with the coupling factor REMOVED
! and the 1/2 period-average projection undone (fel-physics.md sec-unaveraged).

scl_u = (mu_0_vac * c_light) * c_light * dsub / (2 * dgrid * dgrid * beam%slice_spacing)

! Entry handoff: the ramp makes a = 0 here, so the stored (averaged-convention) px IS
! the kinetic momentum. The flag records that px now carries the quiver.

if (first) then
  call fel_assert_averaged_chart (beam, r_name // ' (segment entry)', err)
  if (err) return
  beam%quiver_in_px = .true.
  ustate%active = .true.
  call unavg_ramp_phase_jump ()
endif

allocate (dE_slice(nslice), dU_sp_slice(nslice))
dE_slice = 0
dU_sp_slice = 0
any_err = .false.

! Parallel over slices, the averaged step's own design: each slice
! touches only its own particle arrays and its own field slice (the beam-to-field
! mapping is a bijection), the kernel cache is read-only here (initialized serially
! above), the FFT plan cache is threadprivate, and the per-slice energy lands in
! dE_slice so the final sum is a fixed-order serial reduction. Results are
! bit-identical across thread counts, and the harness checks that.

!$OMP parallel do private(sl, ifld, ux, uy, xx, yy, tau, gam, crsource, crsource_y, s_sub, isub, ip, &
!$OMP&   p_mc, beta, psi_mid, u_s, jhat, wx, wy, ix, iy, on_grid, ehat, ehat_y, wphasor, dgam, cdep, cph, err) &
!$OMP&   reduction(.or.: any_err)
do is = 1, nslice
  sl => beam%slice(is)
  ifld = fel_field_index(slip, is, nslice)

  allocate (crsource(ngrid, ngrid))
  if (two_pol) allocate (crsource_y(ngrid, ngrid))
  allocate (ux(sl%n), uy(sl%n), xx(sl%n), yy(sl%n), tau(sl%n), gam(sl%n))
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

    do ip = 1, sl%n
      call unavg_push (s_sub, dsub/2, xx(ip), yy(ip), ux(ip), uy(ip), tau(ip), gam(ip))
    enddo

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

        ! Two live polarizations: COMPONENT-WISE duals. The instantaneous kinetic
        ! momenta u_x, u_y are real numbers, each working against and depositing into
        ! its own field component. No polarization convention enters at all (manual
        ! sec:field vector convention). The scalar branch below keeps the folded
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
        gam(ip) = gam(ip) + dgam

        ! /u_s, not Genesis's averaged /gamma: the source and the E.v force are exact
        ! duals of one wave equation, and using the SAME u_s the kick used makes the
        ! per-substep energy exchange cancel identically (the diffraction between
        ! substeps is unitary), leaving only physical spontaneous emission in the
        ! ledger. The period-averaged limit shifts by beta_par ~ 5e-9 -- five orders
        ! below the fc checks. A merit choice, not a transcription (sec:unaveraged).
        cdep = cmplx(0.0_rp, 1.0_rp, rp) * exp(cmplx(0.0_rp, -(psi_mid - ks*tau(ip)), rp)) &
               * jhat * scl_u * sl%weight(ip) / u_s
        crsource(ix,   iy)   = crsource(ix,   iy)   + (wx * wy) * cdep
        crsource(ix+1, iy)   = crsource(ix+1, iy)   + ((1-wx) * wy) * cdep
        crsource(ix,   iy+1) = crsource(ix,   iy+1) + (wx * (1-wy)) * cdep
        crsource(ix+1, iy+1) = crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cdep
      endif
    enddo

    do ip = 1, sl%n
      call unavg_push (s_sub + dsub/2, dsub/2, xx(ip), yy(ip), ux(ip), uy(ip), tau(ip), gam(ip))
    enddo

    call fel_field_diffract (wf, ifld, dsub, err)
    any_err = any_err .or. err
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) + 2 * crsource
    if (two_pol) wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) + 2 * crsource_y

    ! The deposit's own energy |dE|^2 = 4|src|^2: the ONE term of the field-energy
    ! increment the kick/deposit duality does not charge to the beam (the beam pays the
    ! cross term 2 Re<E, dE> exactly, manual eq:ledger). Physically this is the
    ! spontaneous emission of the substep. Numerically it is banked here so the
    ! time-dependent ledger closes EXACTLY: E_beam + U_window + U_escaped - U_spont.

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

  deallocate (ux, uy, xx, yy, tau, gam, crsource)
  if (two_pol) deallocate (crsource_y)
enddo
!$OMP end parallel do
if (any_err) return

dE_beam = sum(dE_slice)     ! Fixed-order serial sums: thread-count independent.
dU_spont = sum(dU_sp_slice)

beam%phi0 = beam%phi0 + dz_record * phi0_rate_avg
ustate%s = ustate%s + dz_record

! Exit handoff: the ramp has closed (a = 0), kinetic equals canonical again.

if (last) then
  beam%quiver_in_px = .false.
  ustate%active = .false.
  call unavg_ramp_phase_jump ()
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

! The ramp's slippage compensation, applied as one discrete jump per segment end,
! at the handoffs where the envelope is exactly zero and nothing couples. The
! element's contract is L meters of full-strength undulator. The entry/exit ramps
! are a numerical device (an adiabatic switch-on), and an electron under a ramped
! quiver <u_perp^2> = g^2 aw^2 lags the wave LESS than the contracted hard-edge
! element would have it, by dtau = (gamma aw^2 / 2 u_s^3) INT (1-g^2) ds, which
! is ks*dtau ~ 2.6 rad of optical phase per end at the benchmark parameters
! (sin^2 envelope: INT (1-g^2) = (5/8) l_ramp per end). Uncompensated, every
! segment after the first receives a pre-bunched beam with its bunching-to-field
! phase rotated by that much. The first segment is immune (nothing is bunched
! yet), and that is exactly how it was caught: per-segment ln-power deviations vs
! the averaged mode of {0.0000, +0.08, +0.005, +0.02, +0.01, +0.13}. The jump must
! be discrete AND at the ends: compensating continuously inside the ramp detunes
! the live interaction where the coupling is already substantial (measured -0.9%
! gain on the FIRST segment, doubling with ramp length). In hardware terms this is
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

! One RK4 magnetic push of one particle over step h from segment position s0. All
! per-particle state passes by argument: host-associated variables privatized by the
! caller's OMP region are not redirected inside called procedures, so nothing mutable
! may be host-associated here (und/ustate/inv_beta0 are read-only shared). gamma is
! untouched: B does no work, exactly.

subroutine unavg_push (s0, h, x1, y1, u1, v1, t1, gamma)

real(rp) s0, h, x1, y1, u1, v1, t1, gamma
real(rp) y0(5), k1(5), k2(5), k3(5), k4(5)

y0 = [x1, y1, u1, v1, t1]
call unavg_ode (y0,                s0,         gamma, k1)
call unavg_ode (y0 + (h/2) * k1,   s0 + h/2,   gamma, k2)
call unavg_ode (y0 + (h/2) * k2,   s0 + h/2,   gamma, k3)
call unavg_ode (y0 + h * k3,       s0 + h,     gamma, k4)
y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

x1 = y0(1);  y1 = y0(2)
u1 = y0(3);  v1 = y0(4)
t1 = y0(5)

end subroutine unavg_push

! The exact z-ODEs of ballistic motion in the magnetostatic field:
!   dx/ds = u_x/u_s, du_x/ds = b_y - u_y b_z/u_s, du_y/ds = -b_x + u_x b_z/u_s,
!   dtau/ds = gamma/u_s - 1/beta0.

subroutine unavg_ode (y, s0, gamma, dyds)

real(rp) y(5), s0, gamma, dyds(5)
real(rp) bx, by, bz, us_l

call fel_unavg_bfield (und, ustate, y(1), y(2), s0, bx, by, bz)
us_l = sqrt(gamma**2 - 1 - y(3)**2 - y(4)**2)

dyds(1) = y(3) / us_l
dyds(2) = y(4) / us_l
dyds(3) = by - y(4) * bz / us_l
dyds(4) = -bx + y(3) * bz / us_l

dyds(5) = gamma / us_l - inv_beta0

end subroutine unavg_ode

end subroutine fel_unavg_step

end module fel_unaveraged_mod
