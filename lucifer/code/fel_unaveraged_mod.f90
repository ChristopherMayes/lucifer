!+
! Module fel_unaveraged_mod
!
! The unaveraged verification mode (design brief 6.6; fel-physics.tex sec:unaveraged):
! particles integrated through the undulator's REAL field -- the full Newton-Lorentz
! quiver, no period averaging -- with the radiation field as a co-evolving kick, so the
! coupling factor fc, the harmonic content, and the entry/exit behavior of the averaged
! mode become measured OUTPUTS instead of assumed inputs. MINERVA (Freund & van der
! Slot; minerva-code-analysis.md) is the production existence proof for this physics;
! this mode differs by keeping the GRID field and the Lorentz force where MINERVA
! evaluates modal fields, and it is a verification mode: run on a handful of slices
! where the ~100x cost against the averaged path is irrelevant.
!
! DELIBERATELY ABSENT from this path: fc, faw, and every other period-averaged
! coupling quantity -- those are what this mode measures, and their appearance here
! would make the check circular. The harness greps this file for them.
!
! The step, per substep delta (Strang split, second order):
!
!   1. half magnetic push: RK4 through the analytic normalized field b = curl(a),
!      a = e*A/(m_e c), with gamma exact (B does no work). The vector potential carries
!      a sin^2 amplitude envelope g(s) over n_ramp periods at each end, with g' terms
!      retained in b, so the quiver builds adiabatically and vanishes at the segment
!      ends -- the averaged<->unaveraged handoff (brief 6.6's K/gamma hazard) happens
!      where the two momentum conventions coincide.
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
!      kick/deposit pair EXACT energy duals per substep -- same operands, same bilinear
!      weights, unitary diffraction between -- so the ledger closes to the physical
!      spontaneous-emission term and rounding, by construction. Period-averaging the
!      pair reproduces the averaged mode's fc to O(1-beta_par) ~ 5e-9 (the JJ factor
!      emerges from the figure-8).
!   3. half magnetic push.
!
! Units: E in V/m (wavefront convention), m_electron in eV, b in 1/m. The physical
! field of the scalar envelope: planar E_x = Re[-i Ehat e^{i Psi}]; helical
! (E_x, E_y) = (Re[-i Ehat e^{i Psi}], -Re[Ehat e^{i Psi}])/sqrt(2) -- both give
! intensity |Ehat|^2/(2 Z0), so the power diagnostic is mode-independent.
!
! The magnetic push is classical RK4 on the exact z-ODEs in kinetic variables --
! chosen ON MERIT, not because MINERVA uses it: for a verification mode the currency
! is short-probe ACCURACY, and 4th order is what makes fc measurable at 6e-4 with 20
! steps/period (a 2nd-order symplectic scheme needs ~100 steps/period to match, and
! no explicit symplectic method exists for this non-separable Hamiltonian without
! paying implicit iterations). The structural cost is measured, not argued: gamma is
! conserved EXACTLY by construction (B does no work; gamma changes only in the kick),
! and over the longest benchmark segment -- 266 periods, 5320 steps -- the dark-run
! emittance drifts by <= 3.3e-6, orders below every check. If production-length
! unaveraged runs ever appear (oscillator passes), revisit with a symplectic
! composition; the ballistic check is the instrument that will say when.
!
! Serial over slices BY DECISION: this is a verification mode run on few slices; the
! thread-independence property of the production path is untouched because this path
! never runs inside it.
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
! flat 1, sin^2 down over the last l_ramp -- continuous amplitude AND slope (the
! slope gp feeds the ramp-induced field terms in fel_unavg_bfield). l_ramp = 0 is
! the hard-edge MUTATION configuration; the handoff check exists to catch it.
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
!-

subroutine fel_unavg_bfield (und, ustate, x, y, s, bx, by, bz)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
real(rp) x, y, s, bx, by, bz
real(rp) a0, g, gp, c_u, s_u, fperp

!

g = fel_unavg_envelope(ustate, s, gp)
c_u = cos(und%ku * s)
s_u = sin(und%ku * s)

if (und%helical) then
  a0 = und%aw
  fperp = 1 + (und%ku**2 / 4) * (x*x + y*y)
  ! b = curl(a) with a_x = a0 g c_u fperp, a_y = a0 g s_u fperp:
  bx = -a0 * (gp * s_u + g * und%ku * c_u) * fperp                        ! -d(a_y)/ds
  by =  a0 * (gp * c_u - g * und%ku * s_u) * fperp                        !  d(a_x)/ds
  bz =  a0 * g * (und%ku**2 / 2) * (s_u * x - c_u * y)                    !  d(a_y)/dx - d(a_x)/dy
else
  a0 = sqrt(2.0_rp) * und%aw
  bx = 0
  by = a0 * (gp * c_u - g * und%ku * s_u) * cosh(und%ku * y)            ! d(a_x)/ds
  bz = -a0 * g * c_u * und%ku * sinh(und%ku * y)                        ! -d(a_x)/dy
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
! rate so every diagnostic row stays comparable; the optical carrier used internally
! is Psi = (phi0 - ku*s) - ks*tau. dE_beam returns the weighted particle energy
! change of this step [J] for the energy-ledger check.
!-

subroutine fel_unavg_step (und, ustate, beam, wf, slip, dz_record, first, last, dE_beam, err_flag)

type (fel_und_struct) und
type (fel_unavg_struct) ustate
type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (wavefront_struct), target :: wf
type (fel_slip_struct) slip

real(rp) dz_record, dE_beam
logical first, last, err_flag

real(rp), allocatable :: ux(:), uy(:), xx(:), yy(:), tau(:), gam(:)
complex(rp), allocatable :: crsource(:,:)
real(rp) p0_mc, gamma0b, inv_beta0, ks, dsub, s_sub, phi0_rate_avg, scl_u, dgrid
real(rp) u_s, wx, wy, psi_mid, dgam, p_mc, beta
complex(rp) ehat, jhat, wphasor, cdep
integer is, ip, isub, nslice, ifld, ix, iy, ngrid_arr(3), ngrid
logical on_grid, err
character(*), parameter :: r_name = 'fel_unavg_step'

!

err_flag = .true.
dE_beam = 0

nslice = size(wf%Ex, 3)
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
! and the 1/2 period-average projection undone (fel-physics.tex sec:unaveraged).

scl_u = (mu_0_vac * c_light) * c_light * dsub / (2 * dgrid * dgrid * beam%slice_spacing)

! Entry handoff: the ramp makes a = 0 here, so the stored (averaged-convention) px IS
! the kinetic momentum; the flag records that px now carries the quiver.

if (first) then
  call fel_assert_averaged_chart (beam, r_name // ' (segment entry)', err)
  if (err) return
  beam%quiver_in_px = .true.
  ustate%active = .true.
endif

allocate (crsource(ngrid, ngrid))

do is = 1, nslice
  sl => beam%slice(is)
  ifld = fel_field_index(slip, is, nslice)

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

    do ip = 1, sl%n
      call unavg_push (ip, s_sub, dsub/2)
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
      if (on_grid) then
        ehat =        wf%Ex(ix,   iy,   ifld) * wx * wy
        ehat = ehat + wf%Ex(ix+1, iy,   ifld) * (1-wx) * wy
        ehat = ehat + wf%Ex(ix,   iy+1, ifld) * wx * (1-wy)
        ehat = ehat + wf%Ex(ix+1, iy+1, ifld) * (1-wx) * (1-wy)

        wphasor = cmplx(0.0_rp, -1.0_rp, rp) * ehat * exp(cmplx(0.0_rp, psi_mid - ks*tau(ip), rp))
        dgam = -dsub * real(wphasor * conjg(jhat), rp) / (u_s * m_electron)
        dE_beam = dE_beam + sl%weight(ip) * dgam * m_electron
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
      call unavg_push (ip, s_sub + dsub/2, dsub/2)
    enddo

    call fel_field_diffract (wf, ifld, dsub, err);  if (err) return
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) + 2 * crsource

    s_sub = s_sub + dsub
  enddo

  ! Back to the stored chart. px carries the quiver mid-segment (the flag says so);
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

  deallocate (ux, uy, xx, yy, tau, gam)
enddo

beam%phi0 = beam%phi0 + dz_record * phi0_rate_avg
ustate%s = ustate%s + dz_record

! Exit handoff: the ramp has closed (a = 0), kinetic equals canonical again.

if (last) then
  beam%quiver_in_px = .false.
  ustate%active = .false.
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

! One RK4 magnetic push of particle jp over step h starting at segment position s0.
! gamma is untouched: B does no work, exactly.

subroutine unavg_push (jp, s0, h)

integer jp
real(rp) s0, h
real(rp) y0(5), k1(5), k2(5), k3(5), k4(5)

y0 = [xx(jp), yy(jp), ux(jp), uy(jp), tau(jp)]
call unavg_ode (y0,                s0,         gam(jp), k1)
call unavg_ode (y0 + (h/2) * k1,   s0 + h/2,   gam(jp), k2)
call unavg_ode (y0 + (h/2) * k2,   s0 + h/2,   gam(jp), k3)
call unavg_ode (y0 + h * k3,       s0 + h,     gam(jp), k4)
y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

xx(jp) = y0(1);  yy(jp) = y0(2)
ux(jp) = y0(3);  uy(jp) = y0(4)
tau(jp) = y0(5)

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
