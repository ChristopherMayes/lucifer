!+
! Module fel_collective_mod
!
! In-undulator collective effects for the FEL tracker: wakefields applied as a
! per-slice energy-loss rate, and longitudinal space charge entering the pendulum
! equation as a per-particle ez. Physics, Genesis provenance, and validation:
! lucifer/doc/fel-physics.tex (sec:wakes, sec:spacecharge, sec:seamwake).
!
! Placement in the step: the wake's gamma decrement lands BETWEEN the longitudinal
! advance and the second transverse half step. ez is computed per slice before the RK
! loop and held fixed through the stages, entering dgamma/dz as -ez. Both act in every
! element, interludes included -- the chamber does not end where the undulator does.
!
! fel_resistive_wall_wake is kept a clean, separable routine BY DECISION: it is a
! future port target into Bmad proper as a wake source, and nothing in it knows about
! the FEL. Space charge is transcribed FOR CONSISTENCY with Genesis (directly
! testable), behind this module's interface BY DECISION: Bmad's slice space-charge
! method is suspected the better model long-term, and a Bmad-slice implementation of
! fel_shortrange_ez / fel_longrange_esc is an explicit future task.
!
! Weights: every particle enters the sources with its own charge. The short-range
! source term scales per particle as c*w_j/slice_spacing where Genesis has
! current/npart (identical for uniform weights), and the long-range and wake current
! profiles are the weighted slice currents. Thread safety: fel_shortrange_ez uses
! per-call locals only (callable from the parallel slice loop). The wake update is
! serial at the caller's barrier.
!
! Constants are Bmad's. The Genesis-comparison floors this creates are tabulated in
! the manual (sec:numerics) and doc/validation.md.
!
! Deliberately absent: Genesis's transient wake option (&wake transient/ztrans), the
! incoherent-synchrotron module, harmonics beyond what the solver provides, MPI.
!-

module fel_collective_mod

use fel_beam_mod

implicit none

!+
! Struct fel_wake_struct
!
! Wake configuration and state: the single-particle kernels at wavelength resolution
! over the full window, the external loss, and the per-slice eloss they produce. The
! kernels are built once. The convolution is hoisted when currents cannot change and
! recomputed at the caller's migration stride when they can (the hoist's premise
! predates migration).
!-

type fel_wake_struct
  logical :: on = .false.
  ! Namelist parameters, Genesis &wake names:
  real(rp) :: loss = 0             ! External loss [eV/m], uniform (Genesis allows profiles).
  real(rp) :: radius = 2.5e-3_rp   ! Chamber radius, or half gap if flat [m].
  real(rp) :: conductivity = 0     ! DC conductivity [1/(Ohm m)]. 0: no resistive wake.
  real(rp) :: relaxation = 0       ! AC relaxation distance c*tau [m].
  logical :: roundpipe = .true.    ! Round chamber; false: flat (parallel plates).
  character(8) :: material = ''    ! 'CU' or 'AL' shortcut for conductivity+relaxation.
  real(rp) :: gap = 0              ! Undulator gap [m]. 0: no geometric wake.
  real(rp) :: lgap = 1             ! Period of the gaps [m].
  real(rp) :: hrough = 0           ! Roughness amplitude [m]. 0: no roughness wake.
  real(rp) :: lrough = 1           ! Roughness period [m].
  ! State:
  integer :: ns = 0                ! Kernel length: window slices * sample.
  real(rp) :: ds = 0               ! Kernel resolution: the radiation wavelength [m].
  real(rp), allocatable :: wakeres(:), wakegeo(:), wakerou(:)  ! Kernels [eV/(m electron)].
  real(rp), allocatable :: eloss(:)                            ! Per-slice loss [eV/m].
  logical :: needs_update = .true.
end type

!+
! Struct fel_efield_struct
!
! Space-charge configuration, Genesis &efield names. The solver itself is stateless per
! call (thread safe). This carries only the run facts.
!-

type fel_efield_struct
  logical :: on = .false.          ! Any space charge at all (shortrange if nz*nphi*ngrid set).
  real(rp) :: rmax = 0             ! Radial grid extent scale [m]. Grows adaptively as Genesis's.
  integer :: ngrid = 100           ! Radial grid points.
  integer :: nz = 0                ! Longitudinal harmonics. 0 disables the short-range solve.
  integer :: nphi = 0              ! Azimuthal modes m = -nphi..nphi.
  logical :: longrange = .false.   ! The whole-window longESC term.
end type

!+
! Struct fel_collective_struct
!
! Everything collective the step needs, in one bag: wake state and space-charge
! configuration, plus the per-slice longESC scratch the step refreshes. Default
! constructed = everything off = every existing check bit-identical.
!-

type fel_collective_struct
  type (fel_wake_struct) :: wake
  type (fel_efield_struct) :: efield
  real(rp), allocatable :: long_esc(:)   ! Per-slice longESC [eV/m], refreshed per step.
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_resistive_wall_wake (radius, conductivity, relaxation, roundpipe, ns, ds, wake)
!
! The single-particle resistive-wall wake w(i*ds), i = 0..ns-1, in eV per meter per
! electron (negative = loss), from the NUMERICAL impedance of Bane & Stupakov
! SLAC-PUB-10707 (fel-physics.tex sec:wakes), transcribed from
! Wake::singleWakeResistive with Genesis's exact numerics: k in [0, 100/s0] on 1000
! intervals, s0 = (2 a^2/(Z0 sigma0))^(1/3), flat x-integral on [0,15] with 20000
! points, trapezoid half weights at every endpoint. Kept separable: this routine is
! the future port target into Bmad proper as a wake source, and nothing in it knows
! about the FEL.
!
! Input:
!   radius       -- real(rp): Chamber radius (round) or half-gap (flat) [m].
!   conductivity -- real(rp): Wall conductivity [1/(Ohm m)].
!   relaxation   -- real(rp): AC relaxation time [s] (0 = DC).
!   roundpipe    -- logical: Round chamber (True) or parallel plates (False).
!   ns           -- integer: Number of kernel samples.
!   ds           -- real(rp): Kernel sample spacing [m].
!
! Output:
!   wake(:)      -- real(rp): Single-particle resistive-wall kernel [eV/(m electron)].
!                     The s = 0 entry carries the Bane self-slice half factor.
!-

subroutine fel_resistive_wall_wake (radius, conductivity, relaxation, roundpipe, ns, ds, wake)

real(rp) radius, conductivity, relaxation, ds
logical roundpipe
integer ns
real(rp) wake(:)

integer, parameter :: nk = 1000, nq = 20000
real(rp), parameter :: kappa_max = 100.0_rp, xmax = 15.0_rp
real(rp) z0, a, s0, kmax, dk, k, s_i, acc, dx, x, cos_pref
real(rp), allocatable :: re_z(:), ks_grid(:)
complex(rp) zk, zeta, den, accz
integer i, j, jq

!

wake(1:ns) = 0
if (conductivity <= 0) return

z0 = mu_0_vac * c_light
a = radius
s0 = (2 * a**2 / (z0 * conductivity))**(1.0_rp/3.0_rp)
kmax = kappa_max / s0
dk = kmax / nk

allocate (re_z(0:nk), ks_grid(0:nk))

do j = 0, nk
  k = dk * j
  ks_grid(j) = k
  if (k == 0) then
    re_z(j) = 0
    cycle
  endif
  zeta = cmplx(1.0_rp, -1.0_rp, rp) * sqrt(cmplx(k, 0.0_rp, rp) / &
              (2 * (conductivity / cmplx(1.0_rp, -k * relaxation, rp)) * z0))
  if (roundpipe) then
    den = 1 / zeta - cmplx(0.0_rp, 0.5_rp * k * a, rp)
    zk = (z0 / (twopi * a)) / den
  else
    dx = xmax / (nq - 1)
    accz = 0
    do jq = 0, nq - 1
      x = dx * jq
      den = cmplx(cosh(x), 0.0_rp, rp) * (cosh(x) / zeta - cmplx(0.0_rp, k * a * sinhc(x), rp))
      if (den /= den) cycle                         ! NaN guard, as the original.
      if (jq == 0 .or. jq == nq - 1) then
        accz = accz + 0.5_rp * (z0 / (twopi * a)) / den
      else
        accz = accz + (z0 / (twopi * a)) / den
      endif
    enddo
    zk = accz * dx
  endif
  re_z(j) = real(zk, rp)
enddo

re_z(0) = re_z(0) * 0.5_rp
re_z(nk) = re_z(nk) * 0.5_rp
cos_pref = (2 * c_light / pi) * dk

do i = 1, ns
  s_i = (i - 1) * ds
  acc = 0
  do j = 0, nk
    acc = acc + re_z(j) * cos(ks_grid(j) * s_i)
  enddo
  wake(i) = -e_charge * acc * cos_pref
enddo

!------------------------------------------------------------------------------
contains

function sinhc (xx) result (v)
real(rp) xx, v, x2
if (abs(xx) < 1e-5_rp) then
  x2 = xx * xx
  v = 1 + x2 * (1.0_rp/6.0_rp + x2 / 120.0_rp)
else
  v = sinh(xx) / xx
endif
end function sinhc

end subroutine fel_resistive_wall_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_wake_init (wake, nslice, sample, wavelength, err_flag)
!
! Build the single-particle kernels over the window at wavelength resolution
! (Wake::init): material shortcuts, the three kernels, the self-loading half weight on
! the s = 0 bin (fel-physics.tex sec:wakes).
!
! Input:
!   wake       -- fel_wake_struct: The wake description (radius, materials, roughness).
!   nslice     -- integer: Number of beam slices.
!   sample     -- integer: Slice spacing in wavelengths (Genesis's sample).
!   wavelength -- real(rp): Radiation wavelength [m].
!
! Output:
!   wake       -- fel_wake_struct: Kernels built (%wakeres/%wakegeo/%wakerou, %ns, %ds).
!   err_flag   -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_wake_init (wake, nslice, sample, wavelength, err_flag)

type (fel_wake_struct) wake
integer nslice, sample
real(rp) wavelength
logical err_flag

real(rp) coef, rrough, tau
integer i
character(*), parameter :: r_name = 'fel_wake_init'

!

err_flag = .true.

select case (wake%material)
case ('CU', 'Cu', 'cu')
  wake%conductivity = 5.813e7_rp
  wake%relaxation = 8.1e-6_rp
case ('AL', 'Al', 'al')
  wake%conductivity = 3.571e7_rp
  wake%relaxation = 2.4e-6_rp
case ('')
case default
  call out_io (s_error$, r_name, 'UNKNOWN WAKE MATERIAL: ' // trim(wake%material))
  return
end select

wake%ns = nslice * sample
wake%ds = wavelength
if (allocated(wake%wakeres)) deallocate (wake%wakeres, wake%wakegeo, wake%wakerou, wake%eloss)
allocate (wake%wakeres(wake%ns), wake%wakegeo(wake%ns), wake%wakerou(wake%ns), wake%eloss(nslice))
wake%eloss = 0
wake%needs_update = .true.

call fel_resistive_wall_wake (wake%radius, wake%conductivity, wake%relaxation, &
                              wake%roundpipe, wake%ns, wake%ds, wake%wakeres)

! Geometric (gap) wake, convolved with dI/ds downstream (Wake::singleWakeGeometric).

wake%wakegeo = 0
if (wake%gap > 0) then
  coef = -(mu_0_vac * c_light) * c_light * e_charge / (pi**2 * wake%radius * wake%lgap) &
         * 2 * sqrt(0.5_rp * wake%gap)
  if (.not. wake%roundpipe) coef = coef * 0.956_rp
  do i = 1, wake%ns
    wake%wakegeo(i) = coef * sqrt(wake%ds * (i-1))
  enddo
endif

! Roughness wake (Wake::singleWakeRoughness): complex-q contour in four trapezoid
! segments, Genesis's exact panelization.

wake%wakerou = 0
if (wake%hrough > 0) then
  rrough = pi**3 / wake%lrough**3 * wake%hrough**2 * wake%radius
  coef = rrough / pi * 4 / wake%radius**2 * e_charge / (4 * pi * eps_0_vac)
  do i = 1, wake%ns
    tau = twopi * wake%ds * (i-1) / wake%lrough
    wake%wakerou(i) = coef * (rough_seg(cmplx(0.0_rp, 0.0_rp, rp), cmplx(0.0_rp, 2e-3_rp, rp), 128, tau, rrough) &
                            + rough_seg(cmplx(0.0_rp, 2e-3_rp, rp), cmplx(1.0_rp, 2e-3_rp, rp), 1024, tau, rrough) &
                            + rough_seg(cmplx(1.0_rp, 2e-3_rp, rp), cmplx(1.0_rp, 0.0_rp, rp), 128, tau, rrough) &
                            + rough_seg(cmplx(1.0_rp, 0.0_rp, rp), cmplx(100.0_rp, 0.0_rp, rp), 1024, tau, rrough))
  enddo
endif

! Self-loading theorem: the s = 0 bin of every kernel carries half weight. Genesis
! halves all three (the geometric one is zero there anyway). fel-physics.tex sec:wakes.

wake%wakeres(1) = wake%wakeres(1) * 0.5_rp
wake%wakegeo(1) = wake%wakegeo(1) * 0.5_rp
wake%wakerou(1) = wake%wakerou(1) * 0.5_rp

err_flag = .false.

!------------------------------------------------------------------------------
contains

! One trapezoid segment of the roughness kernel integral: KernelRoughness (the
! integrand with endpoint half weights) followed by TrapIntegrateRoughness.

function rough_seg (q1, q2, n, tau_a, rrough_a) result (val)

complex(rp) q1, q2, dq, q, s, kj, iunit, acc
real(rp) tau_a, rrough_a, val
integer n, j

iunit = cmplx(0.0_rp, 1.0_rp, rp)
dq = (q2 - q1) / (n - 1)
acc = 0

do j = 0, n - 1
  q = q1 + j * dq
  s = (sqrt(2*q + 1) - iunit * sqrt(2*q - 1)) * q / sqrt(4*q*q - 1)
  kj = (s + 1) / (1 - iunit * rrough_a * q * s) / (1 + iunit * rrough_a * q)
  if (j == 0 .or. j == n - 1) kj = kj * 0.5_rp
  acc = acc + exp(-iunit * q * tau_a) * kj
enddo

val = real(acc * dq, rp)

end function rough_seg

end subroutine fel_wake_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_wake_update (wake, beam)
!
! The convolution of Collective::update, one shared-memory node: interpolate the
! per-slice currents to wavelength resolution (zero-padded past the head), convert to
! electrons per bin and the derivative for the geometric term, then for each slice sum
! causally from the evaluation point TOWARD THE HEAD (a trailing slice collects the
! wakes of the charge ahead of it), averaging over the sample steps within the slice.
! Serial by design: the caller holds the barrier. Call once when currents cannot change
! (Genesis's hoisting), and at the migration stride when they can.
!
! Input:
!   wake -- fel_wake_struct: Built kernels.
!   beam -- fel_beam_struct: The beam (per-slice charge -> current profile).
!
! Output:
!   wake -- fel_wake_struct: %eloss per slice, the causal convolution of the kernels
!             with the current profile.
!-

subroutine fel_wake_update (wake, beam)

type (fel_wake_struct) wake
type (fel_beam_struct) beam

real(rp), allocatable :: cur(:), current(:), dcurrent(:)
real(rp) s, wei, dscur
integer nslice, sample, is, ic, i, j, is0, idx

!

if (.not. wake%on) return

nslice = size(beam%slice)
sample = nint(beam%slice_spacing / beam%wavelength)
dscur = beam%slice_spacing

allocate (cur(0:nslice), current(0:wake%ns-1), dcurrent(0:wake%ns-1))

do ic = 1, nslice
  cur(ic-1) = c_light * sum(beam%slice(ic)%weight(1:beam%slice(ic)%n)) / beam%slice_spacing
enddo
cur(nslice) = 0            ! Zero pad past the head (fel-physics.tex sec:wakes).

do is = 0, wake%ns - 1
  s = wake%ds * is
  idx = int(floor(s / dscur))
  wei = 1 - (s - idx * dscur) / dscur
  current(is) = (wei * cur(idx) + (1 - wei) * cur(idx+1)) * wake%ds / (e_charge * c_light)
  dcurrent(is) = -(cur(idx+1) - cur(idx)) * wake%ds / (e_charge * c_light) / dscur
enddo

do ic = 1, nslice
  wake%eloss(ic) = 0
  is0 = (ic - 1) * sample
  do j = 0, sample - 1
    is = is0 + j
    do i = 0, wake%ns - 1 - is
      wake%eloss(ic) = wake%eloss(ic) + current(is+i) * (wake%wakeres(i+1) + wake%wakerou(i+1))
      wake%eloss(ic) = wake%eloss(ic) + dcurrent(is+i) * wake%wakegeo(i+1)
    enddo
  enddo
  wake%eloss(ic) = wake%eloss(ic) / sample + wake%loss
enddo

wake%needs_update = .false.

end subroutine fel_wake_update

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_wake_apply (wake, beam, delz)
!
! Collective::apply in this port's chart: every particle of slice ic loses
! eloss(ic)*delz/m_electron of gamma. Genesis changes gamma and leaves theta untouched.
! Here theta is derived, theta = phi0 + ks*z/beta, so z rescales by beta_new/beta_old
! to hold the phase fixed through the kick (the same bookkeeping fel_advance does at
! its exit). Parallel-safe per slice (pure per-slice arithmetic). The caller may place
! it inside the slice loop or outside.
!
! Input:
!   wake -- fel_wake_struct: %eloss per slice.
!   beam -- fel_beam_struct: The beam.
!   delz -- real(rp): Step length [m].
!
! Output:
!   beam -- fel_beam_struct: Every slice's particles decelerated by eloss*delz.
!-

subroutine fel_wake_apply (wake, beam, delz)

type (fel_wake_struct) wake
type (fel_beam_struct), target :: beam
real(rp) delz
integer ic

!

if (.not. wake%on) return
do ic = 1, size(beam%slice)
  call fel_wake_apply_slice (wake, beam, ic, delz)
enddo

end subroutine fel_wake_apply

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_wake_apply_slice (wake, beam, ic, delz)
!
! One slice's share of fel_wake_apply, per-slice pure so it can sit inside the parallel
! slice loop at Genesis's position in the step (between the longitudinal advance and the
! second transverse half).
!
! Input:
!   wake -- fel_wake_struct: %eloss per slice.
!   beam -- fel_beam_struct: The beam.
!   ic   -- integer: Slice index.
!   delz -- real(rp): Step length [m].
!
! Output:
!   beam -- fel_beam_struct: Slice ic's particles decelerated (chart z rescaled to
!             hold theta, Genesis's convention).
!-

subroutine fel_wake_apply_slice (wake, beam, ic, delz)

type (fel_wake_struct) wake
type (fel_beam_struct), target :: beam
integer ic
real(rp) delz

type (fel_slice_struct), pointer :: sl
real(rp) p0_mc, dg, gam, beta_old, beta_new, p_mc
integer ip

!

if (.not. wake%on) return

sl => beam%slice(ic)
dg = wake%eloss(ic) * delz / m_electron
if (dg == 0) return

p0_mc = fel_p0_mc(beam)
do ip = 1, sl%n
  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta_old = p_mc / gam
  gam = gam + dg
  p_mc = sqrt(gam**2 - 1)
  sl%pz(ip) = (p_mc - p0_mc) / p0_mc
  beta_new = p_mc / gam
  sl%z(ip) = sl%z(ip) * (beta_new / beta_old)     ! Hold theta fixed through the kick.
enddo

end subroutine fel_wake_apply_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_longrange_esc (ef, beam, gamma0, aw, long_esc)
!
! EFieldSolver::longRange, one node: the per-slice long-range space-charge field from
! the whole-window weighted current and rms-size profiles, in Genesis's longESC units
! (eV/m). The CALLER converts at use exactly as Genesis does (fel-physics.tex
! sec:spacecharge): the per-particle ODE ez is fel_shortrange_ez(ip) - long_esc(is)/m_electron.
!
! Input:
!   ef       -- fel_efield_struct: Space-charge configuration.
!   beam     -- fel_beam_struct: The beam.
!   gamma0   -- real(rp): Reference Lorentz factor.
!   aw       -- real(rp): Undulator parameter (parallel-velocity correction).
!
! Output:
!   long_esc(:) -- real(rp): Long-range space-charge Ez per slice [eV/m scale].
!-

subroutine fel_longrange_esc (ef, beam, gamma0, aw, long_esc)

type (fel_efield_struct) ef
type (fel_beam_struct), target :: beam
real(rp) gamma0, aw, long_esc(:)

type (fel_slice_struct), pointer :: sl
real(rp), allocatable :: fcur(:), fsize(:)
real(rp) gamma, scl, efld, dsl, coef, sgn, wsum, x1, y1, x2, y2, w
integer nslice, i, j, ip

!

nslice = size(beam%slice)
long_esc(1:nslice) = 0
if (.not. (ef%on .and. ef%longrange)) return

gamma = gamma0 / sqrt(1 + aw**2)
allocate (fcur(nslice), fsize(nslice))

! Weighted current and transverse size per slice. Genesis's getSize is the PRODUCT of
! the rms sizes, sigma_x*sigma_y: an effective area scale, not a variance sum
! (transcribed wrong once, caught by the SC tier at 1.7e-1, sec:spacecharge). Weighted
! moments where Genesis counts particles: identical for uniform weights, correct
! otherwise. Zero-size guard as the original.

do i = 1, nslice
  sl => beam%slice(i)
  wsum = sum(sl%weight(1:sl%n))
  fcur(i) = c_light * wsum / beam%slice_spacing
  x1 = 0; y1 = 0; x2 = 0; y2 = 0
  do ip = 1, sl%n
    w = sl%weight(ip)
    x1 = x1 + w * sl%x(ip);   x2 = x2 + w * sl%x(ip)**2
    y1 = y1 + w * sl%y(ip);   y2 = y2 + w * sl%y(ip)**2
  enddo
  if (wsum > 0) then
    x1 = x1/wsum; x2 = x2/wsum
    y1 = y1/wsum; y2 = y2/wsum
    fsize(i) = sqrt(abs(x2 - x1*x1)) * sqrt(abs(y2 - y1*y1))
  else
    fsize(i) = 0
  endif
  if (fsize(i) <= 0) fsize(i) = 1
enddo

scl = beam%slice_spacing / pi / c_light / 2 / eps_0_vac

do i = 1, nslice
  efld = 0
  do j = 1, nslice
    dsl = (j - i) * beam%slice_spacing * gamma
    coef = 1 - sqrt(dsl*dsl / (dsl*dsl + fsize(j)))
    if (j > i) then
      sgn = -1
    elseif (j < i) then
      sgn = 1
    else
      sgn = 0
    endif
    efld = efld + sgn * coef * scl * fcur(j) / fsize(j)
  enddo
  long_esc(i) = efld
enddo

end subroutine fel_longrange_esc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_shortrange_ez (ef, beam, sl, gz2, ks, ez)
!
! EFieldSolver::shortRange for one slice: center and radially bin the particles
! (analyseBeam), build the radial Laplace operator (constructLaplaceOperator), and for
! every azimuthal mode m = -nphi..nphi and longitudinal harmonic l = 1..nz solve the
! tridiagonal system for the harmonic potential, accumulating per particle
! ez(ip) += 2*Re(e^{im phi} e^{il theta} u(r_ip)). Units of m_e c^2 per meter, entering
! dgamma/dz as -ez.
!
! Weighted: the source term carries c*w_j/slice_spacing where Genesis has
! current/npart. Thread safe: locals only, callable from the parallel slice loop. The
! adaptive rmax growth is per call (a local), matching the physics without mutating
! shared state. Genesis grows a member and keeps it grown, a difference that only
! affects which slices trigger the growth message.
!
! This is the interface a Bmad-slice space-charge implementation can later fill
! (module header). Callers depend on (beam, slice, gz2, ks) -> ez only.
!
! Input:
!   ef    -- fel_efield_struct: Space-charge configuration (grid, harmonics).
!   beam  -- fel_beam_struct: The beam.
!   sl    -- fel_slice_struct: One slice's packed particles.
!   gz2   -- real(rp): Longitudinal gamma^2 (with the undulator correction).
!   ks    -- real(rp): Radiation wavenumber [1/m].
!
! Output:
!   ez(:) -- real(rp): Short-range space-charge Ez at each particle [eV/m scale].
!-

subroutine fel_shortrange_ez (ef, beam, sl, gz2, ks, ez)

type (fel_efield_struct) ef
type (fel_beam_struct) beam
type (fel_slice_struct) sl
real(rp) gz2, ks, ez(:)

real(rp) xcen, ycen, rbound, rmax_l, dr, tx, ty, radi, coef
real(rp), allocatable :: vol(:), ldig(:), rlog(:), lmid(:), theta_p(:), econst_p(:)
complex(rp), allocatable :: cwork(:), csrc(:), clow(:), cmid(:), cupp(:), celm(:), gam_w(:), cph(:,:)
integer, allocatable :: idxr(:)
integer np, ngrid, m, l, i, ip

!

np = sl%n
ez(1:np) = 0
if (.not. ef%on .or. ef%nz < 1 .or. np < 1) return

ngrid = ef%ngrid

allocate (cwork(np), idxr(np), theta_p(np), econst_p(np), cph(np, ef%nz))
allocate (vol(ngrid), ldig(ngrid+1), rlog(ngrid), lmid(ngrid))
allocate (csrc(ngrid), clow(ngrid), cmid(ngrid), cupp(ngrid), celm(ngrid), gam_w(ngrid))

! analyseBeam: slice centroid, radial extent, bins and azimuthal phases.

xcen = sum(sl%x(1:np)) / np
ycen = sum(sl%y(1:np)) / np

rbound = 0
do ip = 1, np
  tx = sl%x(ip) - xcen
  ty = sl%y(ip) - ycen
  rbound = max(rbound, tx*tx + ty*ty)
enddo
rbound = sqrt(rbound)

rmax_l = ef%rmax
if (rbound > rmax_l) rmax_l = rbound * 1.5_rp
dr = rmax_l / (ngrid - 1)

do ip = 1, np
  tx = sl%x(ip) - xcen
  ty = sl%y(ip) - ycen
  radi = sqrt(tx*tx + ty*ty)
  if (radi > 0) then
    cwork(ip) = cmplx(tx/radi, ty/radi, rp)
  else
    cwork(ip) = cmplx(1.0_rp, 0.0_rp, rp)
  endif
  idxr(ip) = int(floor(radi / dr)) + 1
  if (idxr(ip) > ngrid) idxr(ip) = ngrid
enddo

! The ponderomotive phase, the source weight, and the harmonic phasors e^{il theta}
! are all (m, l)-invariant: computed once per particle here, not 2*(2*nphi+1)*nz
! times inside the mode loops. (Genesis reads a stored theta coordinate there for
! free. This port's theta is derived, so it is cached.) The source side uses the
! phasors' exact conjugates.

do ip = 1, np
  theta_p(ip) = fel_theta(beam, sl, ip, ks)
  ! Weighted source: c*w/slice_spacing where Genesis has current/npart.
  econst_p(ip) = (mu_0_vac * c_light) / m_electron * (c_light * sl%weight(ip) / beam%slice_spacing) / ks
enddo

do l = 1, ef%nz
  do ip = 1, np
    cph(ip, l) = cmplx(cos(l * theta_p(ip)), sin(l * theta_p(ip)), rp)
  enddo
enddo

! constructLaplaceOperator.

vol(1) = pi * dr * dr
rlog(1) = 0.5_rp
ldig(1) = 0
do i = 2, ngrid
  vol(i) = pi * dr * dr * (2*(i-1) + 1)
  ldig(i) = twopi * (i-1)
  rlog(i) = log(real(i, rp) / real(i-1, rp))
enddo
ldig(ngrid+1) = 0

coef = -gz2 / ks**2

do m = -ef%nphi, ef%nphi
  do i = 1, ngrid
    lmid(i) = -ldig(i) - ldig(i+1) - twopi * m * m * rlog(i)
  enddo
  lmid(ngrid) = lmid(ngrid) - twopi * ngrid

  do l = 1, ef%nz

    csrc = 0
    do ip = 1, np
      csrc(idxr(ip)) = csrc(idxr(ip)) + econst_p(ip) * cwork(ip)**(-m) * conjg(cph(ip, l))
    enddo

    do i = 1, ngrid
      csrc(i) = csrc(i) * cmplx(0.0_rp, 1.0_rp / l / vol(i), rp)
      clow(i) = cmplx(coef * ldig(i)   / l / l / vol(i), 0.0_rp, rp)
      cmid(i) = cmplx(1 + coef * lmid(i) / l / l / vol(i), 0.0_rp, rp)
      cupp(i) = cmplx(coef * ldig(i+1) / l / l / vol(i), 0.0_rp, rp)
    enddo

    ! tridiag (Thomas algorithm, Genesis's exact recurrence).

    celm(1) = csrc(1) / cmid(1)
    gam_w(1) = 0
    block
      complex(rp) bet
      bet = cmid(1)
      do i = 2, ngrid
        gam_w(i) = cupp(i-1) / bet
        bet = cmid(i) - clow(i) * gam_w(i)
        celm(i) = (csrc(i) - clow(i) * celm(i-1)) / bet
      enddo
    end block
    do i = ngrid - 1, 1, -1
      celm(i) = celm(i) - gam_w(i+1) * celm(i+1)
    enddo

    do ip = 1, np
      ez(ip) = ez(ip) + 2 * real(cwork(ip)**m * cph(ip, l) * celm(idxr(ip)), rp)
    enddo
  enddo
enddo

end subroutine fel_shortrange_ez

end module fel_collective_mod
