!+
! Module fel_track_mod
!
! The FEL step: period-averaged particle push, source deposition and FFT field solve
! inside an undulator segment. The step works on the Bmad-coordinate packed beam of
! fel_beam_mod. It runs steady state (one slice) or time dependent (many slices with
! slippage).
!
! The physics is transcribed from Genesis 1.3 Version 4, which GPL permits.
! docs/fel-physics.md gives the physics, the Genesis source of each routine and
! the validation (sec-core, sec-transverse, sec-field, sec-slippage).
!
! The step composition (transverse half step, longitudinal advance, transverse half
! step, field step, slippage) is Beam::track plus Gencore's steps 4 and 5.
!
! Time dependence. Beam slice is couples to field slice
!   fel_field_index(slip, is, nslice) = 1 + mod(is-1 + first, nslice)
! where first is Genesis's Field::first. The field record is a circular buffer over the
! wavefront's third index. Slippage rotates the buffer and does not move the data
! (fel_apply_slippage). The caller owns the schedule: what slips when, including the
! drift autophasing of Lattice::calcSlippage. This module owns the mechanics.
!
! Note: Anything reading the field record in time order (dumps, per-slice field
! diagnostics) must unrotate through fel_field_index, as Genesis's writers do
! (sec-slippage).
!
! Coordinates. The stored state is Bmad's (x, px/p0, y, py/p0, z, pz) plus weight (see
! fel_beam_mod). The longitudinal RK4 runs in Genesis's chart (theta, gamma). These are
! per-particle working variables, derived at step entry and converted back at step exit:
!
!   entry:  gamma = sqrt((p0_mc*(1+pz))^2 + 1),  theta = phi0 + ks*z/beta
!   exit:   pz = (sqrt(gamma^2-1) - p0_mc)/p0_mc,  z = -beta*(phi0_new - theta)/ks
!
! theta <-> (phi0, z) is an affine change of variables with a state-independent shift,
! and RK4 is exactly invariant under affine maps. So integrating in theta is the same
! discretization as integrating in z, and the verbatim RK4 transcription is unchanged.
! The per-step gamma <-> pz round trip costs about 1 ulp of gamma per step (the
! subtraction in pz is exact once P is formed). That is the order of the step's own
! rounding. The common reference phase phi0 advances once per step per beam
! (fel_phi0_rate in fel_beam_mod). The split keeps z small and bunch-scale.
!
! Weights. The source deposition scales per particle as c*w_j/slice_spacing where
! Genesis has current/N. The two agree for uniform weights. Nonuniform weights are
! handled since each macroparticle radiates its own charge.
!
! Units. The field is in V/m throughout, the wavefront_struct convention. There is no
! internal unit system. The formulas follow from Genesis's internal-unit forms via the
! exact relation
!   u = E*ks/(sqrt(2)*m_electron)
! which composes the dump scale dfl = u*dgrid*eev/(ks*sqrt(Z0)) with
! E = dfl*sqrt(2*Z0)/dgrid. With this relation:
!   energy exchange coupling = fc*conj(E)/(sqrt(2)*m_electron)
!   source scale             = fc*Z0*sqrt(2)*c*delz*w_j/(4*dgrid^2*slice_spacing)
!   power                    = sum|E|^2*dgrid^2/(2*Z0)
! Constants are Bmad's. m_electron equals Genesis's eev exactly. Z0 = mu_0_vac*c_light
! differs from Genesis's truncated 376.73 by 8e-7 relative. That is the floor of the
! Genesis comparison.
!
! Not in this module: space charge, wakes, incoherent synchrotron radiation, harmonics
! beyond the coupling formula, orbit or field errors, one4one sorting, slice migration,
! and Genesis's periodic time-window boundary. The non-periodic zero fill is what is
! transcribed.
!-

module fel_track_mod

use fel_beam_mod
use fel_collective_mod
use fel_field_mod
use fel_fp32_mod
use fel_device_mod
use wavefront_mod
use fel_timer_mod

use, intrinsic :: iso_c_binding, only: c_double

implicit none

!+
! Struct fel_source_filter_struct
!
! The angular filter on the source term, transcribed from Genesis4's source_filter
! (FieldSolverFFT.cpp:129-146 and 158-170, release 4.6.12). The source is transformed,
! multiplied by a sigmoid in normalized transverse spatial frequency, and only then
! added to the field. The field's own propagation is untouched. So the filter removes
! what a point-like beamlet radiates into the wide angles the grid carries, and does
! not change how the existing field diffracts (docs/startup-noise.md).
!-

type fel_source_filter_struct
  logical :: on = .false.
  real(rp) :: xcut = 1, ycut = 1   ! Sigmoid edge per plane, in half-Nyquist units.
  real(rp) :: width = 1            ! Sigmoid width, same units. Genesis4's sigmoid.
  logical :: mutate = .false.      ! The check's self-test: filter the field, not the source.
  real(rp) :: angle = 0            ! The edge as an angle [rad], 0 when set by xcut instead.
end type

!+
! Struct fel_und_struct
!
! One undulator segment, with constant parameters along it. That is one contiguous run
! of aw > 0 steps in Genesis's unrolled lattice. The segment has nstep equal steps of
! dz = l/nstep, with nstep = round(l/delz_target). kx, ky carry Genesis's unroll
! scaling by ku^2 (fel-physics.md sec-element).
!-

type fel_und_struct
  real(rp) :: aw = 0          ! rms undulator parameter.
  real(rp) :: ku = 0          ! Undulator wavenumber twopi/lambdau [1/m].
  real(rp) :: kx = 0          ! Natural focusing, Genesis4's kx * ku^2 [1/m^2].
  real(rp) :: ky = 0          ! Natural focusing, Genesis4's ky * ku^2 [1/m^2].
  real(rp) :: ax = 0, ay = 0  ! Transverse offset of the undulator field [m].
  real(rp) :: tilt = 0        ! Wiggle-plane tilt [rad]: planar wiggles along
  real(rp) :: cos_t = 1, sin_t = 0   !   e = (cos_t, sin_t). Cached for the maps.
  complex(rp) :: pol(2) = [(1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp)]
                              ! Polarization 2-vector on the (Ex, Ey) pair
                              ! (fel-physics.md sec-field): planar (cos t, sin t),
                              ! helical (1, -i)/sqrt(2). The averaged kick reads
                              ! E_eff = conj(pol).E and the deposit writes pol*src.
                              ! With one polarization this is the scalar path exactly.
  logical :: helical = .false.
  logical :: bmad_transport = .false.  ! Transverse maps: True for the flattened Bmad
                                       !   periodic-wiggler kernel, False for the
                                       !   transcribed TrackBeam. Set by the walk from
                                       !   run.transport_model (default 'bmad').
  integer :: nstep = 0        ! Number of integration steps over the segment.
  real(rp) :: dz = 0          ! Step length [m].
  integer :: source_model = 0 ! fel_source_deposit$ (default) or fel_source_coherent$.
  type (fel_source_filter_struct) :: filter   ! The source-term angular filter, off by default.
end type

!+
! Struct fel_coherent_struct
!
! One slice's coherent-source summary (fel-physics.md sec-coherent-source), computed
! each step from the just-advanced particles. It holds:
!   The source phasor S. This is exactly the sum the deposit would scatter, and
!     SUM crsource = S is the normalization.
!   The charge-weighted centroid and 2x2 second moments (widths and tilt). The
!     centering and tilt extend the paper, which assumes an on-axis untilted beam.
!   The Laguerre-Gauss order-0/1 sums for the global kappa fit.
!   The Gaussianity guard metric (excess kurtosis, worst plane).
!-

type fel_coherent_struct
  complex(rp) :: S = 0                   ! Source phasor: SUM part_j (sin+icos theta_j).
  real(rp) :: x0 = 0, y0 = 0             ! Charge-weighted centroid [m].
  real(rp) :: sxx = 0, sxy = 0, syy = 0  ! Charge-weighted central second moments [m^2].
  real(rp) :: g2sum = 0                  ! |g0|^2 + |g1|^2 (Eq 20 sums, charge-weighted).
  real(rp) :: b2 = 0                     ! |B|^2, same weighting (kappa's Eq 26 ratio).
  real(rp) :: kurt = 0                   ! max plane |excess kurtosis| (guard metric).
  real(rp) :: m_ind = 1                  ! Independent transverse samples at load time
                                         !   (beamlet copies share coordinates). The
                                         !   guard tests significance and not raw kurtosis.
                                         !   Sample kurtosis spreads as sqrt(24/m_ind).
  real(rp) :: wsum = 0                   ! Slice charge [C], which the guard tests. Edge
                                         !   slices resample from few source points and
                                         !   look degenerate while their source is
                                         !   negligible.
  logical :: ok = .false.                ! Enough charge and spread to fit.
end type

! The field set's types (fel_slip_struct, fel_bank_struct, fel_field_struct) live in
! fel_field_mod. This module uses and re-exports them. The device seam takes the set
! whole, so it needs to see the types.


! Cached kernels for fel_field_step, as in FieldSolverFFT::init. There is one entry per
! wavelength (per harmonic), rebuilt when its grid or step changes. This is module
! state. fel_field_kernel_init builds it serially before any parallel slice loop, and
! it is read-only after that.
! fel_field_step never rebuilds. It looks its entry up (fel_kernel_index) and errors on
! a miss, since a rebuild inside a parallel loop would race. The source accumulator is
! a local of fel_field_step, so concurrent slices deposit into their own.

type fel_kernel_struct
  complex(rp), allocatable :: k2(:,:)      ! -i (kx^2+ky^2)/(2 ks), FFT order.
  complex(rp), allocatable :: exp_k2(:,:)  ! exp(K2 * dz), the step propagator.
  real(rp), allocatable :: sigmoid(:,:)    ! The source filter's sigmoid, FFT order. Built
                                           !   only when a caller asks for it.
  integer :: ngrid = 0
  real(rp) :: dgrid = 0, ks = 0, dz = 0
end type

! The source model (fel-physics.md sec-coherent-source): the standard per-particle
! deposit, or the SIMPLEX-hybrid coherent-Gaussian source (Tanaka, PRAB 27, 030703
! (2024), arXiv:2310.20197). The coherent source drops the spatially incoherent
! artifact. The slice bunch factor keeps the physical shot noise. The transverse shape
! is a guarded Gaussian from phasor sums and charge moments.

integer, parameter :: fel_source_deposit$ = 0
integer, parameter :: fel_source_coherent$ = 1

! The FEL mode of an element is its fel_method attribute. Bmad names the two values,
! averaged$ and unaveraged$ (bmad_struct.f90), so there is no mode enum here. The
! element's tracking_method still sets how one particle crosses it outside the FEL walk.

type (fel_kernel_struct), allocatable, target, private, save :: fel_kernels(:)

!+
! Subroutine fel_sincos (theta, s, c)
!
! Routine to compute the (sin, cos) pair with one libm sincos call (the C routine in
! fel_sincos.c).
!
! Note: The result is *not* bit-identical to gfortran's sin and cos intrinsics. The two
! differ by one ulp on about 2e-6 of arguments of the theta domain (docs/validation.md,
! "The particle-path cost").
!-

interface
  subroutine fel_sincos (theta, s, c) bind(c, name = 'fel_sincos_c')
    import c_double
    real(c_double), value :: theta
    real(c_double) :: s, c
  end subroutine
end interface

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_ele_as_wiggler (orbit, ele, param, err_flag, finished, track)
!
! Routine for the track1_custom hook. Outside a driver's own FEL walk, an FEL element (a
! wiggler tracked by an FEL method) is a periodic wiggler, and Bmad's standard kernel
! tracks it. The reference time and energy pass inside bmad_parser and any seam-side
! track1 come through here. So the element carries the resonant undulation delay from
! Bmad's own code.
!
! Note: This routine must stay at module scope. gfortran implements pointers to
! internal procedures with stack trampolines. The non-executable stack of Apple Silicon
! turns those into a segfault at the first call.
!
! Input:
!   orbit     -- coord_struct: Starting orbit.
!   ele       -- ele_struct: The FEL wiggler/undulator element.
!   param     -- lat_param_struct: Lattice parameters.
!   track     -- track_struct, optional: Ignored (no step-by-step recording here).
!
! Output:
!   orbit     -- coord_struct: Orbit at the element exit.
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!   finished  -- logical: Set True (the hook fully handles the element).
!-

subroutine fel_ele_as_wiggler (orbit, ele, param, err_flag, finished, track)

type (coord_struct) orbit
type (ele_struct) ele
type (lat_param_struct) param
logical err_flag, finished
type (track_struct), optional :: track

!

err_flag = .false.
finished = .true.
call track_a_wiggler (orbit, ele, param)

end subroutine fel_ele_as_wiggler

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_assert_wiggler_sane (ele)
!
! Routine to check that an FEL element is one this tracker can track. It is called from
! the setup pass that reads fel_method, and from nowhere else. Each message names the
! attribute so a lattice author knows what to fix.
! An FEL element tracks as a periodic wiggler outside the walk, so Bmad handles every
! field_calc on its own and the setup pass is early enough.
!
! Input:
!   ele -- ele_struct: The FEL wiggler/undulator element to validate.
!
! Output:
!   None. Every violation is fatal: the message names the attribute and the program
!   stops with a nonzero exit.
!-

subroutine fel_assert_wiggler_sane (ele)

type (ele_struct) ele

character(*), parameter :: r_name = 'fel_assert_wiggler_sane'

! These refusals must end the program with a nonzero exit, and the validation suite's refusal
! checks test for it. err_exit cannot be used here. Its integer-divide traceback does
! not trap on arm64, and its final bare stop exits 0. So the stops below stay.

! A superimposed element cuts an undulator into super_slaves whose field_calc refers to
! the lord. Nothing here reads the lord: the undulator's parameters, its method and its
! ramps are the element's own. Such a lattice is refused before any frame is written.
! So no dump from a sliced endpoint inside an undulator exists, and
! fel_checkpoint_eligible guards that boundary should one ever exist.

if (ele%slave_status == super_slave$ .or. ele%field_calc == refer_to_lords$) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT IS A SUPER_SLAVE, ITS FIELD REFERRED TO ITS LORD: ' // trim(ele%name), &
                                 'A SUPERIMPOSED ELEMENT CANNOT CUT AN UNDULATOR THIS TRACKER TRACKS.', &
                                 'POSSIBLE SOLUTION: PLACE THE SUPERIMPOSED ELEMENT BETWEEN UNDULATORS.')
  stop 1
endif
if (ele%field_calc /= planar_model$ .and. ele%field_calc /= helical_model$) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT FIELD_CALC MUST BE PLANAR_MODEL OR HELICAL_MODEL', &
                                 '(A FIELDMAP GETS NO FOCUSING HERE): ' // trim(ele%name))
  stop 1
endif
if (ele%value(b_max$) <= 0) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT HAS ZERO B_MAX (NO FIELD, NO RESONANCE,', &
                                 'AND BMAD ITSELF WOULD NOT WARN): ' // trim(ele%name))
  stop 1
endif
if (ele%value(l_period$) <= 0) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT HAS ZERO L_PERIOD (OSC_AMPLITUDE WOULD BE', &
                                 'SILENTLY ZERO): ' // trim(ele%name))
  stop 1
endif
if (ele%value(kx$) /= 0) then
  call out_io (s_fatal$, r_name, 'THE BMAD KX ROLL-OFF ATTRIBUTE IS NOT YET MAPPED TO THE FEL FOCUSING SPLIT.', &
                                 'POSSIBLE SOLUTION: SET KX = 0 ON: ' // trim(ele%name))
  stop 1
endif
if (ele%value(tilt$) /= 0 .and. ele%field_calc == helical_model$) then
  call out_io (s_fatal$, r_name, 'TILT ON A HELICAL FEL ELEMENT IS A ROTATION OF A CIRCULARLY', &
                                 'SYMMETRIC FIELD -- A NO-OP THAT READS AS A MISTAKE: ' // trim(ele%name))
  stop 1
endif

! The integration step. With neither ds_step nor num_steps set, Bmad's bookkeeper falls
! back to bmad.default_ds_step (0.2 m, about thirteen periods of a 15 mm undulator).
! That step gives wrong FEL physics with no complaint, and Bmad only warns. So an FEL
! element must set its step.

if (ele%value(ds_step$) == bmad_com%default_ds_step) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT HAS NO INTEGRATION STEP; SET DS_STEP (GENESIS''S DELZ)', &
                                 'OR NUM_STEPS ON: ' // trim(ele%name), &
                                 '(ITS DS_STEP EQUALS BMAD.DEFAULT_DS_STEP, THE UNSET FALLBACK. IF YOU', &
                                 'REALLY WANT THAT EXACT VALUE, SET IT EXPLICITLY VIA NUM_STEPS.)')
  stop 1
endif

end subroutine fel_assert_wiggler_sane

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_mat6_as_wiggler (ele, param, start_orb, end_orb, err_flag)
!
! Routine for the make_mat6_custom hook, the transfer-matrix companion of
! fel_ele_as_wiggler. The standard periodic-wiggler kernel tracks with matrix
! propagation and fills ele%mat6, ele%vec0 and end_orb as make_mat6_bmad does.
!
! Note: An FEL element's mat6_calc_method resolves to custom (auto follows an FEL
! tracking method), and Bmad's bookkeeping calls make_mat6_custom_ptr unconditionally.
! A program that leaves the pointer null segfaults at a jump to address zero.
!
! Input:
!   ele        -- ele_struct: The FEL wiggler/undulator element.
!   param      -- lat_param_struct: Lattice parameters.
!   start_orb  -- coord_struct: Orbit at the element entrance.
!
! Output:
!   ele        -- ele_struct: %mat6, %vec0 filled per the make_mat6_bmad convention.
!   end_orb    -- coord_struct: Orbit at the element exit.
!   err_flag   -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_mat6_as_wiggler (ele, param, start_orb, end_orb, err_flag)

type (ele_struct), target :: ele
type (coord_struct) start_orb, end_orb
type (lat_param_struct) param
logical err_flag

!

err_flag = .false.
end_orb = start_orb
call mat_make_unit (ele%mat6)
call track_a_wiggler (end_orb, ele, param, ele%mat6, make_matrix = .true.)
ele%vec0 = end_orb%vec - matmul(ele%mat6, start_orb%vec)

end subroutine fel_mat6_as_wiggler

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_field_index (slip, is, nslice) result (ifld)
!
! Routine to map beam slice is (1-based) to its field slice in the rotated record:
!   ifld = 1 + mod(is-1 + first, nslice)
! This is the Fortran form of Genesis's (is + field->first) % field.size()
! (fel-physics.md sec-slippage). The same mapping restores time order when reading the
! record out. Field slice fel_field_index(slip, is, nslice) is the field at time window
! position is.
!-

elemental function fel_field_index (slip, is, nslice) result (ifld)

type (fel_slip_struct), intent(in) :: slip
integer, intent(in) :: is, nslice
integer ifld

!

ifld = 1 + mod(is - 1 + slip%first, nslice)

end function fel_field_index

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_apply_slippage (slip, wf, slippage, bank, harm, err_flag)
!
! Routine to account the slippage of one step and rotate the field record when it
! crosses a slice boundary.
! Transcribed from Control::applySlippage (fel-physics.md sec-slippage), reduced to one
! shared-memory node. The MPI ring exchange between adjacent ranks is the identity when
! the ring has one member. This node is both rank 0 and rank size-1, so the
! non-periodic zero fill always applies to the transmitted slice. What remains is
! Genesis's exact bookkeeping. Slippage accumulates in radiation wavelengths. While
! |accuslip| exceeds 0.8*sample, the record rotates one slice: the slice coupled to the
! head of the time window is zeroed and becomes the new tail. That is, radiation leaves
! the window at the head and fresh vacuum enters behind the bunch. Backward slippage is
! transcribed the same way with the direction reversed.
!
! The caller owns the schedule (per-step undulator slippage and the drift autophasing
! of Lattice::calcSlippage). The caller calls this after the field step, Gencore's
! step 5. A steady-state record (timerun false) is left untouched, as in Genesis.
!
! Note: A slice transmitted backward leaves at the tail. The escaped file and the pulse
! file place banked light only ahead of the window. So backward slippage with a bank
! present is refused.
!
! Input:
!   slip        -- fel_slip_struct: Slippage state of this field record.
!   wf          -- wavefront_struct: The field record.
!   slippage    -- real(rp): Slippage of this step [radiation wavelengths].
!   bank        -- fel_bank_struct, optional: Receives the transmitted planes when present.
!   harm        -- integer, optional: Harmonic number of this field. Default is 1.
!
! Output:
!   slip, wf    -- Updated state and rotated record.
!   err_flag    -- logical, optional: Set True on the refusal above. False otherwise.
!-

subroutine fel_apply_slippage (slip, wf, slippage, bank, harm, err_flag)

type (fel_slip_struct) slip
type (wavefront_struct) wf
type (fel_bank_struct), optional :: bank
integer, optional :: harm
logical, optional :: err_flag
real(rp) slippage, t_slice
character(*), parameter :: r_name = 'fel_apply_slippage'
complex(wf_rp), allocatable :: grow(:,:,:)
integer nslice, last, direction

!

! The transmitted slice's light time for the escape bank: slice_spacing / c. For a
! harmonic field wf%wavelength is the fundamental's / harm, while sample stays in
! fundamental units (every field shares the window and rotates in lockstep). So the
! product is scaled back up by harm.

t_slice = (slip%sample * wf%wavelength / c_light)
if (present(harm)) then
  if (harm /= 1) t_slice = t_slice * harm
endif

if (present(err_flag)) err_flag = .true.
if (present(bank)) bank%n = 0
if (.not. slip%timerun) then
  if (present(err_flag)) err_flag = .false.
  return
endif

slip%accuslip = slip%accuslip + slippage
nslice = size(wf%Ex, 3)

do while (abs(slip%accuslip) > slip%sample * 0.8_rp)
  direction = 1
  if (slip%accuslip < 0) direction = -1
  if (direction < 0 .and. present(bank)) then
    call out_io (s_error$, r_name, 'BACKWARD SLIPPAGE WITH THE ESCAPED-FIELD BANK ON, WHICH IT IS BY DEFAULT.', &
                                   'A SLICE TRANSMITTED FROM THE TAIL HAS NO PLACE AHEAD OF THE WINDOW', &
                                   'IN THE ESCAPED AND PULSE FILES.', &
                                   'POSSIBLE SOLUTION: SET KEEP_ESCAPED_FIELD = F.')
    return
  endif
  slip%accuslip = slip%accuslip - slip%sample * direction

  ! The transmitted slice: the last of the record in time order for forward slippage,
  ! and the first for backward. 0-based here as in Genesis (sec-slippage).

  last = mod(slip%first + nslice - 1, nslice)
  if (direction < 0) last = mod(last + 1, nslice)

  ! One node: send and receive to self, then the non-periodic zero fill of the
  ! transmitted slice. Its energy leaves the simulation here. It is banked first so the
  ! time-dependent energy ledger can close. The banked energy is the slice's power (as
  ! in fel_field_diag) times its light time slice_spacing/c. Wakes would be a second,
  ! unbanked exit channel. The unaveraged mode, whose ledger this feeds, refuses them.

  slip%u_escaped = slip%u_escaped + sum(real(wf%Ex(:,:,last+1), rp)**2 + aimag(wf%Ex(:,:,last+1))**2) &
                     * wf%dx**2 / (2 * (mu_0_vac * c_light)) * t_slice
  if (allocated(wf%Ey)) then
    slip%u_escaped = slip%u_escaped + sum(real(wf%Ey(:,:,last+1), rp)**2 + aimag(wf%Ey(:,:,last+1))**2) &
                       * wf%dx**2 / (2 * (mu_0_vac * c_light)) * t_slice
  endif

  ! Bank the light itself when asked: the transmitted slice, copied before the zero.

  if (present(bank)) then
    if (.not. allocated(bank%plane)) then
      allocate (bank%plane(size(wf%Ex,1), size(wf%Ex,2), 4))
    elseif (bank%n == size(bank%plane, 3)) then
      call move_alloc (bank%plane, grow)
      allocate (bank%plane(size(grow,1), size(grow,2), 2*size(grow,3)))
      bank%plane(:,:,1:size(grow,3)) = grow
      deallocate (grow)
    endif
    if (allocated(wf%Ey) .and. .not. allocated(bank%plane_y)) then
      allocate (bank%plane_y(size(wf%Ex,1), size(wf%Ex,2), size(bank%plane,3)))
    elseif (allocated(bank%plane_y) .and. size(bank%plane_y,3) < size(bank%plane,3)) then
      call move_alloc (bank%plane_y, grow)
      allocate (bank%plane_y(size(grow,1), size(grow,2), size(bank%plane,3)))
      bank%plane_y(:,:,1:size(grow,3)) = grow
      deallocate (grow)
    endif
    bank%n = bank%n + 1
    bank%plane(:,:,bank%n) = wf%Ex(:,:,last+1)
    if (allocated(wf%Ey)) bank%plane_y(:,:,bank%n) = wf%Ey(:,:,last+1)
  endif

  wf%Ex(:, :, last+1) = 0
  if (allocated(wf%Ey)) wf%Ey(:, :, last+1) = 0

  ! The transmitted slice becomes the start of the record.

  slip%first = last
  if (direction < 0) slip%first = mod(last + 1, nslice)
enddo

if (present(err_flag)) err_flag = .false.

end subroutine fel_apply_slippage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_device_apply_slippage (dev, im, slip, wf, slippage, harm)
!
! Routine to account one step's slippage while the field is resident on the device.
! The bookkeeping is fel_apply_slippage's. The rotation's two data motions act on the
! resident record through fel_device_mod. The transmitted slice's energy is read back
! (one slice) and its device slice is zeroed. The host FP64 record stays stale, as
! between any two readbacks. There is no escape bank, since a device run keeps none
! (fel_setup). Each member of the field set takes one call.
!
! Input:
!   dev         -- fel_device_struct: The resident device state.
!   im          -- integer: The field-set member.
!   slip        -- fel_slip_struct: Slippage state of this member's record.
!   wf          -- wavefront_struct: The field record. Only its geometry is read.
!   slippage    -- real(rp): Slippage of this step [radiation wavelengths].
!   harm        -- integer, optional: The member's harmonic number, as fel_apply_slippage
!                    takes it. Default 1.
!
! Output:
!   dev, slip   -- Updated state. The device record is rotated.
!-

subroutine fel_device_apply_slippage (dev, im, slip, wf, slippage, harm)

type (fel_device_struct) dev
type (fel_slip_struct) slip
type (wavefront_struct) wf
real(rp) slippage, t_slice
integer im, nslice, last, direction
integer, optional :: harm

!

! The transmitted slice's light time, as in fel_apply_slippage. sample counts
! fundamental wavelengths and a harmonic's wf%wavelength is the fundamental's over harm.
! So the product is scaled back up by harm. Without the factor a third harmonic's
! escaped energy would be a third of what left the record.

t_slice = (slip%sample * wf%wavelength / c_light)
if (present(harm)) then
  if (harm /= 1) t_slice = t_slice * harm
endif

if (.not. slip%timerun) return

slip%accuslip = slip%accuslip + slippage
nslice = size(wf%Ex, 3)

do while (abs(slip%accuslip) > slip%sample * 0.8_rp)
  direction = 1
  if (slip%accuslip < 0) direction = -1
  slip%accuslip = slip%accuslip - slip%sample * direction

  last = mod(slip%first + nslice - 1, nslice)
  if (direction < 0) last = mod(last + 1, nslice)

  slip%u_escaped = slip%u_escaped + fel_device_slice_energy(dev, im, last+1, wf%dx) * t_slice
  call fel_device_zero_slice (dev, im, last+1)

  slip%first = last
  if (direction < 0) slip%first = mod(last + 1, nslice)
enddo

end subroutine fel_device_apply_slippage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_und_coupling (und, h) result (fc)
!
! Routine to return the coupling factor fc(h), transcribed from Undulator::fc.
! Helical couples the fundamental with strength aw and nothing else. Planar couples odd
! harmonics through
!   JJ = J_{(h-1)/2}(xi) - J_{(h+1)/2}(xi),  xi = h/2 * aw^2/(1+aw^2)
! with sign (-1)^((h-1)/2). Planar does not couple even harmonics (a polarization
! limitation of the planar coupling).
!
! Input:
!   und -- fel_und_struct: Undulator parameters (aw, helicity).
!   h   -- integer: Harmonic number.
!
! Output:
!   fc  -- real(rp): The coupling factor aw*JJ(h) (planar) or aw (helical, h = 1).
!            Zero for helical harmonics.
!-

function fel_und_coupling (und, h) result (fc)

type (fel_und_struct) und
integer h, h0, h1
real(rp) fc, xi, coup

!

coup = und%aw

if (und%helical) then
  if (h == 1) then
    fc = coup
  else
    fc = 0
  endif
  return
endif

if (mod(h, 2) == 1) then
  xi = und%aw**2
  xi = 0.5_rp * xi / (1 + xi) * h
  h0 = (h - 1) / 2
  h1 = h0 + 1
  fc = coup * (bessel_jn(h0, xi) - bessel_jn(h1, xi)) * (-1)**h0
else
  fc = 0
endif

end function fel_und_coupling

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function faw (und, x, y) result (value)
!
! Routine to return the transverse dependence of the undulator field, the first order
! form used in the particle gather. Transcribed from Undulator::faw.
!
! Input:
!   und  -- fel_und_struct: Undulator parameters (roll-off kx/ky, offsets, tilt).
!   x, y -- real(rp): Transverse position [m].
!
! Output:
!   value -- real(rp): aw(x,y), the transversely rolled-off undulator parameter.
!-

function faw (und, x, y) result (value)

type (fel_und_struct) und
real(rp) x, y, value, dx, dy, dt

!

dx = x - und%ax
dy = y - und%ay
if (und%sin_t /= 0) then     ! The roll-off lives in the wiggle frame (tilted planar).
  dt = und%cos_t * dx + und%sin_t * dy
  dy = -und%sin_t * dx + und%cos_t * dy
  dx = dt
endif
value = 1 + 0.5_rp * (und%kx * dx*dx + und%ky * dy*dy)

end function faw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function faw2 (und, x, y) result (value)
!
! Routine to return the square of the transverse dependence, used in the source
! deposition. Transcribed from Undulator::faw2.
!
! Input:
!   und  -- fel_und_struct: Undulator parameters (roll-off kx/ky, offsets, tilt).
!   x, y -- real(rp): Transverse position [m].
!
! Output:
!   value -- real(rp): aw^2(x,y), Genesis's faw2 (the squared roll-off used by the
!            ponderomotive phase).
!-

function faw2 (und, x, y) result (value)

type (fel_und_struct) und
real(rp) x, y, value, dx, dy, dt

!

dx = x - und%ax
dy = y - und%ay
if (und%sin_t /= 0) then     ! Wiggle-frame roll-off, as in faw.
  dt = und%cos_t * dx + und%sin_t * dy
  dy = -und%sin_t * dx + und%cos_t * dy
  dx = dt
endif
value = 1 + und%kx * dx*dx + und%ky * dy*dy

end function faw2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_und_envelope (l, l_ramp, s, gp) result (g)
!
! Routine to return the undulator amplitude envelope at s into the segment: sin^2 up
! over l_ramp, flat 1, sin^2 down over the last l_ramp. Amplitude and slope are
! continuous. The slope gp feeds the ramp-induced field terms in fel_unavg_bfield.
! The unaveraged mode integrates through this envelope.
!
! Note: l_ramp = 0 is a hard edge, which the handoff check flags.
!
! Input:
!   l      -- real(rp): Segment length [m].
!   l_ramp -- real(rp): Ramp length at each end [m]. Zero is a hard edge.
!   s      -- real(rp): Position inside the segment [m].
!
! Output:
!   gp     -- real(rp): The envelope derivative dg/ds [1/m].
!   g      -- real(rp): The field envelope g(s) (sin^2 ramps, 1 in the body).
!-

function fel_und_envelope (l, l_ramp, s, gp) result (g)

real(rp) l, l_ramp, s, g, gp, arg

!

g = 1;  gp = 0
if (l_ramp <= 0) return

if (s < l_ramp) then
  arg = pi * s / (2 * l_ramp)
  g = sin(arg)**2
  gp = (pi / l_ramp) * sin(arg) * cos(arg)
elseif (s > l - l_ramp) then
  arg = pi * (l - s) / (2 * l_ramp)
  g = sin(arg)**2
  gp = -(pi / l_ramp) * sin(arg) * cos(arg)
endif

end function fel_und_envelope

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_und_step (und, beam, ff, coll, err_flag, fp32, dev)
!
! Routine to advance the whole beam and the field set one integration step of length
! und%dz inside an undulator segment. The order is Genesis's:
!   1. Transverse half step.
!   2. Longitudinal RK4. Every field is gathered once and the couplings are held
!      through the stages.
!   3. Transverse half step.
!   4. Source deposition and field solve per field, harmonic loop outermost.
! The common phase phi0 advances once per step, with the undulator's ku as the
! reference wavenumber. Each beam slice couples to its field slice through
! fel_field_index per field.
!
! Note: Slippage is *not* applied here. The caller owns the schedule and calls
! fel_apply_slippage after this step (Gencore's step 5).
!
! Input:
!   und         -- fel_und_struct: Segment parameters.
!   beam        -- fel_beam_struct: The beam. Steady state: one slice.
!   ff          -- fel_field_struct(:): The field set [V/m]. Entry 1 is the fundamental,
!                    one field slice per beam slice each.
!   coll        -- fel_collective_struct: Collective-effects state.
!   fp32        -- fel_fp32_struct, optional: The FP32 lockstep check state.
!   dev         -- fel_device_struct, optional: The device backend state.
!
! Output:
!   beam, ff    -- Advanced one step.
!   err_flag    -- logical: Set True on error.
!   fp32        -- fel_fp32_struct, optional: This step's rows recorded.
!   dev         -- fel_device_struct, optional: Device state advanced.
!-

subroutine fel_track_und_step (und, beam, ff, coll, err_flag, fp32, dev)

type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
type (fel_collective_struct) coll
type (fel_fp32_struct), optional :: fp32
type (fel_device_struct), optional :: dev
logical err_flag

type (fel_device_par_struct) dpar
real(rp) ks, phi0_new, phi0_old, kappa, b2_tot, g2_tot
real(rp) fp_rtmp, fp_gridmax
real(rp), allocatable :: dev_s0(:,:,:)
type (fel_coherent_struct), allocatable :: coh(:)
integer is, io, nslice, nslice_f, ngrid_arr(3)
logical any_err, err, ok, fp_on, dev_on, dev_twin, dev_prod
character(*), parameter :: r_name = 'fel_track_und_step'

!

err_flag = .true.

call fel_assert_averaged_chart (beam, r_name, err)
if (err) return

nslice = size(ff(1)%wf%Ex, 3)
if (size(beam%slice) /= nslice) then
  call out_io (s_error$, r_name, 'BEAM HAS \i0\ SLICES BUT THE FIELD RECORD HAS \i0\ .', &
                                 i_array = [size(beam%slice), nslice])
  return
endif

! ks, phi0 and the slippage schedule are all fundamental quantities (field 1). A
! harmonic field enters only through its coupling fc(h), its phase h*theta and its
! own diffraction kernel.

ks = twopi / ff(1)%wf%wavelength
phi0_new = beam%phi0 + und%dz * fel_phi0_rate(ks, und%ku, fel_p0_mc(beam))

! Everything cross-slice happens serially, before and between the parallel loops: the
! phi0 advance above, the kernel builds and the long-range space-charge profile below,
! and slippage in the caller. Inside the loops each slice touches only its own
! particle arrays and its own field slice (the beam-to-field mapping is a bijection).
! Each slice's arithmetic does not depend on the thread that runs it. So results are
! bit-identical across thread counts, and the validation suite checks this.

call fel_tic (fel_t_und_prep$)
do io = 1, size(ff)
  ngrid_arr = wavefront_shape(ff(io)%wf)
  call fel_field_kernel_init (ngrid_arr(1), ff(io)%wf%dx, twopi / ff(io)%wf%wavelength, und%dz, und%filter)
enddo
call fel_toc (fel_t_und_prep$)

call fel_tic (fel_t_sc_profile$)
if (.not. allocated(coll%long_esc)) allocate (coll%long_esc(nslice))
call fel_longrange_esc (coll%efield, beam, fel_gamma0(beam), und%aw, coll%long_esc)
call fel_toc (fel_t_sc_profile$)

! The FP32 lockstep check (fel_fp32_mod). Its twin mirrors the fundamental-only,
! collective-free advance. So any configuration outside that is refused at first use.
! With the device on, the twin role passes to the device (fel_device_mod) under the
! same refusals. Without the FP32 check the device is the run itself, and this routine
! reduces to the serial prep plus one resident step.

fp_on = .false.
if (present(fp32)) fp_on = fp32%on
dev_on = .false.
if (present(dev)) dev_on = dev%on
dev_twin = dev_on .and. fp_on
dev_prod = dev_on .and. .not. fp_on

if (dev_on) then
  call fill_device_par (dpar, ok)
  if (.not. ok) return
  if (fp_on) dpar%mutate = merge(1, 0, fp32%mutate)
endif

! The production device step. The beam and the field set are resident (the walk
! uploaded them at element entry). So the whole step is one command buffer: transverse
! maps, the push in the offset-energy tick-phase chart gathering every member, the
! deposit into every member and the FP32 field solve of every plane. Nothing returns to
! the host here. Setup refuses everything this step does not cover (collectives,
! radiation, migration, the coherent source).

if (dev_prod) then
  call set_device_kernels ()
  call fel_tic (fel_t_device$)
  call fel_device_step_run (dev, dpar, beam%phi0, phi0_new, ks, err)
  call fel_toc (fel_t_device$)
  if (err) return
  beam%phi0 = phi0_new
  err_flag = .false.
  return
endif

! The CPU twin covers the fundamental of one polarization. The device twin covers the
! whole set. So the first two refusals below apply only to the CPU twin.

if (fp_on .and. .not. fp32%checked) then
  if (size(ff) > 1 .and. .not. dev_twin) then
    call out_io (s_error$, r_name, 'FP32_CHECK DOES NOT COVER HARMONIC FIELD SETS.')
    return
  endif
  if (allocated(ff(1)%wf%Ey) .and. .not. dev_twin) then
    call out_io (s_error$, r_name, 'FP32_CHECK DOES NOT COVER TWO-POLARIZATION FIELDS.')
    return
  endif
  if (und%source_model == fel_source_coherent$) then
    call out_io (s_error$, r_name, 'FP32_CHECK DOES NOT COVER THE COHERENT SOURCE MODEL.')
    return
  endif
  if (coll%wake%on) then
    call out_io (s_error$, r_name, 'FP32_CHECK DOES NOT COVER WAKES.')
    return
  endif
  if (coll%efield%active) then
    call out_io (s_error$, r_name, 'FP32_CHECK DOES NOT COVER SPACE CHARGE.')
    return
  endif
  fp32%checked = .true.
endif
if (fp_on) then
  fp_rtmp = fel_und_coupling(und, 1) / (sqrt(2.0_rp) * m_electron)
  ngrid_arr = wavefront_shape(ff(1)%wf)
  fp_gridmax = (ngrid_arr(1) - 1) * ff(1)%wf%dx / 2
endif
if (dev_twin) then

  ! Sized to the slices as they are now, and not to the device's rectangle. Migration at
  ! the previous element's end can leave a slice fuller than the rectangle.
  ! fel_device_twin_begin grows the rectangle only when it stages this snapshot.

  allocate (dev_s0(6, maxval(beam%slice(:)%n), size(beam%slice)))
endif

! Per-slice sequence in Genesis's order (Beam::track): transverse half, longitudinal
! advance (ez inside the RK), the wake's gamma decrement, transverse half. All four act
! on one slice only, so one loop is arithmetic-identical to Genesis's four sweeps.

! The twin's field record spans the window. It is allocated here because two threads
! could reach a first-use allocation inside the loop below at once.

if (fp_on .and. .not. dev_twin) call fel_fp32_field_prep (fp32, size(ff(1)%wf%Ex, 1), &
                                                                size(ff(1)%wf%Ex, 2))

call fel_tic (fel_t_particles$)
any_err = .false.
!$OMP parallel do private(err) reduction(.or.: any_err)
do is = 1, size(beam%slice)
  if (dev_twin) then

    ! The shared state the device twin advances from, taken before the leading half
    ! step. The device runs the whole step, transverse maps included. The CPU twin below
    ! copies later and takes the transverse state from FP64.

    block
      integer nn
      nn = beam%slice(is)%n
      dev_s0(1, 1:nn, is) = beam%slice(is)%x(1:nn)
      dev_s0(2, 1:nn, is) = beam%slice(is)%px(1:nn)
      dev_s0(3, 1:nn, is) = beam%slice(is)%y(1:nn)
      dev_s0(4, 1:nn, is) = beam%slice(is)%py(1:nn)
      dev_s0(5, 1:nn, is) = beam%slice(is)%z(1:nn)
      dev_s0(6, 1:nn, is) = beam%slice(is)%pz(1:nn)
    end block
  endif
  if (und%bmad_transport) then
    call fel_transverse_track_bmad (und, beam, beam%slice(is), und%dz/2, .true.)
  else
    call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
  endif
  if (fp_on .and. .not. dev_twin) then

    ! The shared state the twin advances from, copied at fel_advance's own sequence
    ! point. The copy is block-local, like the seam's scratch bunch, one per slice per step.

    block
      real(rp), allocatable :: cx(:), cy(:), cpx(:), cpy(:), cz(:), cpz(:)
      integer nn
      nn = beam%slice(is)%n
      cx = beam%slice(is)%x(1:nn);    cy = beam%slice(is)%y(1:nn)
      cpx = beam%slice(is)%px(1:nn);  cpy = beam%slice(is)%py(1:nn)
      cz = beam%slice(is)%z(1:nn);    cpz = beam%slice(is)%pz(1:nn)
      call fel_advance (und, beam, beam%slice(is), ff, und%dz, phi0_new, coll, is, err)
      any_err = any_err .or. err
      call fel_fp32_twin_slice (fp32, is, cx, cy, cpx, cpy, cz, cpz, beam%slice(is), beam, &
            ff(1)%wf%Ex(:, :, fel_field_index(ff(1)%slip, is, size(ff(1)%wf%Ex, 3))), &
            ff(1)%wf%dx, ff(1)%wf%dy, fp_gridmax, und%aw, und%ku, und%kx, und%ky, &
            und%ax, und%ay, und%cos_t, und%sin_t, fp_rtmp, und%dz, ks, und%ku, &
            beam%phi0, phi0_new)
    end block
  else
    call fel_advance (und, beam, beam%slice(is), ff, und%dz, phi0_new, coll, is, err)
    any_err = any_err .or. err
  endif
  call fel_wake_apply_slice (coll%wake, beam, is, und%dz)
  if (und%bmad_transport) then
    call fel_transverse_track_bmad (und, beam, beam%slice(is), und%dz/2, .false.)
  else
    call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
  endif
enddo
!$OMP end parallel do
call fel_toc (fel_t_particles$)

! A refused collective solve stops the step before the field is touched.

if (any_err) then
  err_flag = .true.;  return
endif

! The FP32 check's step epilogue (rows, stream, guard) runs after the field twin at the
! end of this routine. So the field rows of this step join the same line.

phi0_old = beam%phi0
beam%phi0 = phi0_new

! The device twin's uploads and step encode come before the production solve. The
! lockstep field image must round from the pre-solve FP64 record. The encoded step
! then overlaps the CPU's own solve. The rows read back after both finish.

if (dev_twin) then
  call set_device_kernels ()
  call fel_tic (fel_t_device$)
  call fel_device_twin_begin (dev, fp32, beam, ff, dev_s0, dpar, phi0_old, phi0_new, ks, err)
  call fel_toc (fel_t_device$)
  if (err) return
endif

! The coherent source (fel-physics.md sec-coherent-source). First the per-slice
! phasor, moments and LG sums from the just-advanced particles, slice-parallel. Then
! the one global kappa (Tanaka Eq 26). Its integrals run over the whole window, so this
! is the cross-slice serial point. Then the guard. A slice whose transverse charge
! profile is significantly non-Gaussian is refused, and the message gives its number. The
! Gaussian model would bias the gain there. The distorted-beam checks set the
! thresholds.

if (und%source_model == fel_source_coherent$) then
  call fel_tic (fel_t_source$)
  allocate (coh(size(beam%slice)))
  !$OMP parallel do
  do is = 1, size(beam%slice)
    call fel_coherent_prep (und, beam, beam%slice(is), ks, coh(is))
  enddo
  !$OMP end parallel do

  ! The kappa fit's window sums (Eq 26), charge-weighted, since g_q and B are built from
  ! normalized weights. An unweighted sum would let near-empty edge slices spoil the
  ! global width fit, since their normalized phasors are O(1) noise. wsum^2 keeps both
  ! sums quadratic in the same weighting.

  b2_tot = 0;  g2_tot = 0
  do is = 1, size(beam%slice)
    b2_tot = b2_tot + coh(is)%b2 * coh(is)%wsum**2
    g2_tot = g2_tot + coh(is)%g2sum * coh(is)%wsum**2
  enddo
  do is = 1, size(beam%slice)

    ! An edge slice (Gaussian tail, resampled from few source points) looks degenerate
    ! while its source is negligible. Only slices with at least 0.05 of the largest
    ! slice charge are held to the Gaussian profile.

    if (coh(is)%wsum < 0.05_rp * maxval(coh%wsum)) cycle
    if (coh(is)%ok .and. coh(is)%kurt > max(0.5_rp, 5 * sqrt(24.0_rp / coh(is)%m_ind))) then

      ! The guard tests significance and not raw kurtosis. The sample statistic spreads
      ! as sqrt(24/m_ind) for m_ind independent transverse samples. So a small-M Gaussian
      ! beam does not trip the guard on noise. A structured profile (double horn, about
      ! -1.5) trips it at any M where it is resolved.

      call out_io (s_error$, r_name, 'COHERENT SOURCE: slice \i0\ transverse charge profile is NOT', &
        'GAUSSIAN ENOUGH (excess kurtosis \es10.2\ , threshold \es10.2\ at this sample size):', &
        'the coherent-Gaussian model would bias the gain here. Use source_model = "deposit".', &
        i_array = [is], r_array = [coh(is)%kurt, max(0.5_rp, 5 * sqrt(24.0_rp / coh(is)%m_ind))])
      return
    endif
  enddo
  kappa = 1
  if (g2_tot > 0) kappa = sqrt(b2_tot / (4 * pi * g2_tot))
  call fel_toc (fel_t_source$)
endif

! Field solve, harmonic loop outermost. Each pass is self-contained (fel-physics.md
! sec-field-set). The deposit reads the just-advanced particles, so every field sees
! the same beam state.

call fel_tic (fel_t_solve$)
do io = 1, size(ff)
  nslice_f = size(ff(io)%wf%Ex, 3)
  any_err = .false.
  if (und%source_model == fel_source_coherent$) then
    !$OMP parallel do private(err) reduction(.or.: any_err)
    do is = 1, size(beam%slice)
      call fel_field_step (und, beam, beam%slice(is), ff(io)%wf, &
                           fel_field_index(ff(io)%slip, is, nslice_f), und%dz, ff(io)%harm, ks, err, &
                           coh(is), kappa)
      any_err = any_err .or. err
    enddo
    !$OMP end parallel do
  else
    !$OMP parallel do private(err) reduction(.or.: any_err)
    do is = 1, size(beam%slice)
      call fel_field_step (und, beam, beam%slice(is), ff(io)%wf, &
                           fel_field_index(ff(io)%slip, is, nslice_f), und%dz, ff(io)%harm, ks, err)
      any_err = any_err .or. err
    enddo
    !$OMP end parallel do
  endif
  if (any_err) return
enddo
call fel_toc (fel_t_solve$)

! The FP32 field twin, serial after the production solve (fel_fp32_mod). It makes one
! FP32 image of each slice's deposit-transform-propagate step against the FP64 record
! that now carries it. scl_w and the kernel come from the same sources the production
! step used. With the device as the twin, its step, encoded before the solve, is read
! back here instead and fills every row of the same stream.

if (fp_on) then
  if (dev_twin) then
    call fel_tic (fel_t_device$)
    call fel_device_twin_rows (dev, fp32, beam, ff, dpar, phi0_old, phi0_new, ks)
    call fel_toc (fel_t_device$)
  else
    block
      real(rp) fp_scl
      integer ikk2, ifl
      fp_scl = fel_und_coupling(und, 1) * (mu_0_vac * c_light) * sqrt(2.0_rp) * c_light * und%dz
      fp_scl = fp_scl / (4 * ff(1)%wf%dx * ff(1)%wf%dx * beam%slice_spacing)
      ikk2 = fel_kernel_index(ngrid_arr(1), ff(1)%wf%dx, twopi / ff(1)%wf%wavelength, und%dz)
      do is = 1, size(beam%slice)
        ifl = fel_field_index(ff(1)%slip, is, size(ff(1)%wf%Ex, 3))
        call fel_fp32_field_twin (fp32, is, beam%slice(is), beam, ff(1)%wf%Ex(:,:,ifl), &
              fel_kernels(ikk2)%exp_k2, fp_scl, ff(1)%wf%dx, fp_gridmax, und%kx, und%ky, &
              und%ax, und%ay, und%cos_t, und%sin_t, ks)
      enddo
    end block
  endif
  call fel_fp32_step_close (fp32, err)
  if (err) return
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine fill_device_par (par, ok)
!
! Routine to fill one step's device constants from the same sources the FP64 step
! reads. These are the undulator segment, the beam's references, the set's one grid
! and the element's polarization pair. Per member they are fel_advance's coupling
! fc(h)/(sqrt(2) m_e) and fel_field_step's source scale at that harmonic. The Bmad
! transverse map's k1 locals are fel_transverse_track_bmad's forms. ok is False when
! the deposit bound cannot be set.
!-

subroutine fill_device_par (par, ok)

type (fel_device_par_struct) par
real(rp) g_max, p0_mc_l, s_bound, r_max
integer ngrid_l(3), io
logical ok

!

p0_mc_l = fel_p0_mc(beam)
par%dz = und%dz
par%ks = ks
par%ku = und%ku
par%aw = und%aw
par%qres = 2 * und%ku / ks
par%e0 = p0_mc_l**2 / fel_gamma0(beam)
par%gam0 = fel_gamma0(beam)
par%p0_mc = p0_mc_l
par%kx = und%kx;  par%ky = und%ky
par%ax = und%ax;  par%ay = und%ay
par%cos_t = und%cos_t;  par%sin_t = und%sin_t

if (und%helical) then
  g_max = und%aw * und%ku / p0_mc_l
  par%k1x = -0.5_rp * g_max**2
  par%k1y = par%k1x
else
  g_max = und%aw * sqrt(2.0_rp) * und%ku / p0_mc_l
  par%k1x = 0
  par%k1y = -0.5_rp * g_max**2
endif

ngrid_l = wavefront_shape(ff(1)%wf)
par%gridmax = (ngrid_l(1) - 1) * ff(1)%wf%dx / 2
par%dgrid = ff(1)%wf%dx

par%harm = 0
par%rtmp = 0
par%scl_w = 0
do io = 1, size(ff)
  par%harm(io) = ff(io)%harm
  par%rtmp(io) = fel_und_coupling(und, ff(io)%harm) / (sqrt(2.0_rp) * m_electron)
  par%scl_w(io) = fel_und_coupling(und, ff(io)%harm) * (mu_0_vac * c_light) * sqrt(2.0_rp) * c_light * und%dz
  par%scl_w(io) = par%scl_w(io) / (4 * ff(1)%wf%dx * ff(1)%wf%dx * beam%slice_spacing)
enddo
par%pol_re = real(und%pol, rp)
par%pol_im = aimag(und%pol)

! The deposit's fixed-point scale. The bound is the charge the run loaded over the
! gamma floor (fel_device_setup), times this element's deposit scale, times the
! roll-off at the grid's corner. That is every macroparticle in one cell at the worst
! coupling the mesh carries. fel_device_dep_scale turns the bound into a scale and
! refuses a bound it cannot carry.

ok = .false.

! A readback found a particle below the gamma the bound was built on. The bound no
! longer holds, so nothing further is deposited. fel_device_readback printed the gamma.

if (dev%dep_breach) return

! A bound on the roll-off the kernel applies. The kernel takes the particle's offset
! from the undulator axis, rotates it when the segment is tilted, and forms
! 1 + kx*ddx^2 + ky*ddy^2. A particle can sit anywhere in the box, so the offset reaches
! gridmax plus the axis offset. A rotation can move that whole radius onto either axis,
! against the larger of kx and ky. A bound at gridmax per axis, with no offset and no
! rotation, is short by up to a factor of 4.5 in the roll-off term.

r_max = par%gridmax + max(abs(par%ax), abs(par%ay))
s_bound = sqrt(1 + 2 * max(abs(par%kx), abs(par%ky)) * r_max**2) * &
          maxval(par%scl_w(1:size(ff))) * dev%dep_qbound
par%dep_gam_floor = dev%dep_gam_floor
call fel_device_dep_scale (dev, s_bound, 'the run''s charge at an eighth of its lowest ' // &
        'loaded gamma, the element''s source scale and the roll-off at the box corner', &
        par%dep_scale, err)
if (err) return

ok = .true.

par%first = ff(1)%slip%first
par%helical = merge(1, 0, und%helical)
par%mutate = 0
par%source_filter = merge(1, 0, und%filter%on)
par%nfield = size(ff)
par%npol = merge(2, 1, allocated(ff(1)%wf%Ey))
par%pad = 0

end subroutine fill_device_par

!+
! Subroutine set_device_kernels ()
!
! Routine to hand the device every member's propagator for this step, each at its own
! wavelength, from the cache fel_field_step reads. fel_device_set_kernel keys the
! upload, so an unchanged kernel is not uploaded again.
!-

subroutine set_device_kernels ()

integer io, ik

!

do io = 1, size(ff)
  ik = fel_kernel_index(size(ff(io)%wf%Ex, 1), ff(io)%wf%dx, twopi / ff(io)%wf%wavelength, und%dz)
  call fel_device_set_kernel (dev, io, fel_kernels(ik)%exp_k2)
  if (und%filter%on) call fel_device_set_filter (dev, io, fel_kernels(ik)%sigmoid)
enddo

end subroutine set_device_kernels

end subroutine fel_track_und_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_interlude_genesis (qf, length, beam, ff, coll, err_flag)
!
! Routine to advance one field-free interlude element (a drift or a quadrupole) as
! Genesis does, in one integration step:
!   1. Transverse half step.
!   2. Longitudinal advance, with the path-length term sampled at mid element and
!      Genesis's drift reference xku = ks*0.5/gamma0/gamma0. The division order is
!      kept for bit identity (fel-physics.md sec-interlude).
!   3. The wake's gamma decrement.
!   4. Transverse half step.
!   5. Field diffraction with zero source.
!
! The longitudinal advance has two paths. With space charge off, the RK4 collapses
! exactly (slope independent of theta, gamma constant). The collapsed step is kept for
! bit identity with the code before collective effects. With space charge on, gamma
! changes inside the step. So the full RK4 runs with rpart = 0 and the per-particle
! ez. That is Genesis's own path in drifts (BeamSolver::advance with aw = 0).
!
! This is the Genesis interlude model, transcribed. A full-lattice run with it
! isolates what the Bmad seam changes. The seam integrates the path-length term
! exactly through the quad map, where this samples it once (docs/validation.md). The
! production configuration is the seam.
!
! Note: Slippage is *not* applied here. The caller schedules it after the step, as
! with fel_track_und_step.
!
! Input:
!   qf       -- real(rp): Quadrupole focusing strength of the interlude [1/m^2].
!   length   -- real(rp): Interlude length [m].
!   beam     -- fel_beam_struct: The beam.
!   ff(:)    -- fel_field_struct: The field set (phi0 advances with the drift).
!   coll     -- fel_collective_struct: Collective terms (space charge acts here).
!
! Output:
!   beam     -- fel_beam_struct: Advanced through the interlude.
!   ff(:)    -- fel_field_struct: phi0 advanced.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_track_interlude_genesis (qf, length, beam, ff, coll, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_collective_struct) coll
real(rp) qf, length
logical err_flag

type (fel_und_struct) und0
real(rp) xks, xku, qquad, phi0_new, gamma0, p0_mc
real(rp) px_g, py_g, gam, beta, theta, btpar, btpar0, slope, q_hat, p_mc
integer is, io, ip, nslice, nslice_f, ngrid_arr(3)
logical err, any_err, sc_active

!

err_flag = .true.

wf => ff(1)%wf                                    ! The fundamental sets the beam side.
gamma0 = fel_gamma0(beam)                         ! Genesis's gammaref.
p0_mc = fel_p0_mc(beam)
xks = twopi / wf%wavelength
xku = xks * 0.5_rp / gamma0 / gamma0
qquad = qf * gamma0                               ! fel-physics.md sec-interlude.
q_hat = qquad / p0_mc
phi0_new = beam%phi0 + length * fel_phi0_rate(xks, xku, p0_mc)
nslice = size(wf%Ex, 3)

sc_active = coll%efield%active .and. (coll%efield%nz >= 1 .or. coll%efield%longrange)
if (.not. allocated(coll%long_esc)) allocate (coll%long_esc(size(beam%slice)))
call fel_longrange_esc (coll%efield, beam, gamma0, 0.0_rp, coll%long_esc)

! The per-particle temporaries live at routine scope, so the parallel loop must make
! them private explicitly. A missed one is a race, which the validation suite's
! thread-count-independence check catches.

any_err = .false.
!$OMP parallel do private(sl, ip, gam, beta, theta, px_g, py_g, btpar, btpar0, slope, p_mc, err) &
!$OMP     reduction(.or.: any_err)
do is = 1, size(beam%slice)
  sl => beam%slice(is)

  call interlude_transverse_half (sl, q_hat, length/2)

  if (sc_active) then

    call interlude_advance_full_rk (sl, is, err)
    any_err = any_err .or. err

  else

    ! theta advance, the collapsed RK4. With no field and no space charge the slope does
    ! not depend on theta and is constant through the stages. So RK4 reduces to one
    ! exact step. btpar = 1 + px_g^2 + py_g^2 (aw = 0), px_g = gamma*beta_x = px * p0_mc.

    do ip = 1, sl%n
      p_mc = p0_mc * (1 + sl%pz(ip))
      gam = sqrt(p_mc**2 + 1)
      beta = p_mc / gam
      theta = beam%phi0 + xks * sl%z(ip) / beta

      px_g = sl%px(ip) * p0_mc
      py_g = sl%py(ip) * p0_mc
      btpar = 1 + px_g*px_g + py_g*py_g
      btpar0 = sqrt(1 - btpar / (gam * gam))
      slope = xks * (1 - 1/btpar0) + xku
      theta = theta + length * slope

      ! Back to z: tau = (phi0_new - theta)/ks. pz does not change (no energy change), so
      ! beta does not either.
      sl%z(ip) = -beta * (phi0_new - theta) / xks
    enddo
  endif

  call fel_wake_apply_slice (coll%wake, beam, is, length)

  call interlude_transverse_half (sl, q_hat, length/2)
enddo
!$OMP end parallel do
if (any_err) return

beam%phi0 = phi0_new

! Field diffraction: Genesis's kernel with zero source. aw = 0 makes the coupling
! exactly zero, as Genesis skips the deposition outside an undulator. The beam slice to
! field slice mapping is a bijection, so looping over beam slices covers every field
! slice once. The kernels are built serially first. A lattice can open with an
! interlude, so no undulator step may have built them yet. Harmonic loop outermost,
! each field diffracting at its own wavelength.

do io = 1, size(ff)
  ngrid_arr = wavefront_shape(ff(io)%wf)
  call fel_field_kernel_init (ngrid_arr(1), ff(io)%wf%dx, twopi / ff(io)%wf%wavelength, length)
enddo

und0%aw = 0
do io = 1, size(ff)
  nslice_f = size(ff(io)%wf%Ex, 3)
  any_err = .false.
  !$OMP parallel do private(err) reduction(.or.: any_err)
  do is = 1, size(beam%slice)
    call fel_field_step (und0, beam, beam%slice(is), ff(io)%wf, &
                         fel_field_index(ff(io)%slip, is, nslice_f), length, ff(io)%harm, xks, err)
    any_err = any_err .or. err
  enddo
  !$OMP end parallel do
  if (any_err) return
enddo

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine interlude_transverse_half (sl, q_hat, dz)
!
! Routine to apply half an interlude's transverse map to one slice. That is the
! focusing map per plane at the interlude's quad strength (Genesis's drift/quad
! composite).
!-

subroutine interlude_transverse_half (sl, q_hat, dz)

type (fel_slice_struct) sl
real(rp) q_hat, dz, gz_hat
integer ip

do ip = 1, sl%n
  gz_hat = sqrt((1 + sl%pz(ip))**2 - sl%px(ip)**2 - sl%py(ip)**2)
  call fel_apply_focus (dz,  q_hat, sl%x(ip), sl%px(ip), gz_hat)
  call fel_apply_focus (dz, -q_hat, sl%y(ip), sl%py(ip), gz_hat)
enddo

end subroutine interlude_transverse_half

!------------------------------------------------------------------------------
!+
! Subroutine interlude_advance_full_rk (sl, is, serr)
!
! Routine to advance one slice's longitudinal coordinates through the interlude.
! This is Genesis's drift path (BeamSolver::advance with aw = 0): the full RK4 with
! rpart = 0 and the per-particle space-charge ez held through the stages. gamma changes
! here, so the chart round trip at exit is fel_advance's. The routine uses locals only
! and is safe inside the parallel slice loop. serr is set when the slice's short-range
! solve is refused.
!-

subroutine interlude_advance_full_rk (sl, is, serr)

type (fel_slice_struct) sl
integer is
logical serr

real(rp), allocatable :: ez(:)
real(rp) gam, beta, theta, px_g, py_g, btpar, esc_loss, p_mc
integer ip

allocate (ez(max(1, sl%n)))
call fel_shortrange_ez (coll%efield, beam, sl, gamma0**2, xks, ez, serr)   ! gz2 at aw = 0.
if (serr) return          ! No force for this slice, so nothing to integrate.
esc_loss = -coll%long_esc(is) / m_electron

do ip = 1, sl%n
  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam
  theta = beam%phi0 + xks * sl%z(ip) / beta

  px_g = sl%px(ip) * p0_mc
  py_g = sl%py(ip) * p0_mc
  btpar = 1 + px_g*px_g + py_g*py_g

  call fel_runge_kutta (length, xks, xku, btpar, cmplx(0.0_rp, 0.0_rp, rp), &
                        ez(ip) + esc_loss, gam, theta)

  p_mc = sqrt(gam**2 - 1)
  sl%pz(ip) = (p_mc - p0_mc) / p0_mc
  sl%z(ip) = -(p_mc / gam) * (phi0_new - theta) / xks
enddo

end subroutine interlude_advance_full_rk

end subroutine fel_track_interlude_genesis

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_transverse_track (und, beam, sl, delz)
!
! Routine to advance the transverse coordinates over delz inside an undulator: drift
! plus natural focusing. Transcribed from TrackBeam::track. The quad, corrector and
! chicane branches are dropped since they never overlap undulator fields here.
!
! In the stored normalization the per particle longitudinal factor is
!
!   gz_hat = gammaz/p0 = sqrt((1+pz)^2 - px^2 - py^2 - (aw/p0_mc)^2)
!
! That is Genesis's gammaz = sqrt(gamma^2 - 1 - aw^2 - px_g^2 - py_g^2) divided by p0.
! The maps read the same with px in Bmad units:
!
!   drift:  x += px*delz/gz_hat
!   quad:   foc^2 = q_hat/gz_hat with q_hat = q/p0_mc
!           x'  = a1 x + a2 px/gz_hat
!           px' = a3 x gz_hat + a1 px
!
! The effective strengths (fel-physics.md sec-natfocus), with Genesis's reference gamma:
!   qnat_{x,y} = k_{x,y}*aw^2/(gamma0*betpar0),  betpar0 = sqrt(1 - (1+aw^2)/gamma0^2)
!
! Input:
!   und  -- fel_und_struct: Undulator parameters.
!   beam -- fel_beam_struct: The beam (reference energy).
!   sl   -- fel_slice_struct: One slice's packed particles.
!   delz -- real(rp): Step length [m].
!
! Output:
!   sl   -- fel_slice_struct: Transverse coordinates advanced by delz.
!-

subroutine fel_transverse_track (und, beam, sl, delz)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
real(rp) delz
real(rp) betpar0, qx_hat, qy_hat, aw2, aw_p0_sq, gz_hat, gamma0, p0_mc
integer ip

!

gamma0 = fel_gamma0(beam)
p0_mc = fel_p0_mc(beam)
aw2 = und%aw**2
betpar0 = sqrt(1 - (1 + aw2)/gamma0**2)

qx_hat = und%kx * aw2 / gamma0 / betpar0 / p0_mc
qy_hat = und%ky * aw2 / gamma0 / betpar0 / p0_mc
aw_p0_sq = (und%aw / p0_mc)**2

do ip = 1, sl%n
  gz_hat = sqrt((1 + sl%pz(ip))**2 - sl%px(ip)**2 - sl%py(ip)**2 - aw_p0_sq)
  call fel_apply_focus (delz, qx_hat, sl%x(ip), sl%px(ip), gz_hat)
  call fel_apply_focus (delz, qy_hat, sl%y(ip), sl%py(ip), gz_hat)
enddo

end subroutine fel_transverse_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_transverse_track_bmad (und, beam, sl, delz, leading)
!
! Routine to advance the transverse coordinates with the transverse maps of Bmad's own
! periodic-wiggler kernel (track_a_wiggler.f90:90-186), flattened (fel-physics.md
! sec-element). This is the default transport. Quadrupole bodies go through
! quad_mat2_calc with the per-particle 1/rel_p^2 chromatic scaling and the
! half-octupole edge kicks. They use the tracking-local k1 values, and never the stored
! k1x/k1y attributes, whose helical sign disagrees.
!
! The kernel's z bookkeeping (path-length dz terms, low-energy correction, the
! end-of-element undulation factor) is left out. The ponderomotive phase evolution,
! including the aw^2 and px^2+py^2 path terms, lives in fel_advance's RK. Applying it
! here too would double-count.
!
! Every leading half-step applies half the octupole kick, then the quad body. Every
! trailing half-step applies the body, then half the kick. Adjacent half-kicks between
! steps sum to Bmad's full inter-step kick. The segment faces get the half kick, as in
! Bmad's n_step loop. Two model differences from the transcribed TrackBeam: Bmad's
! 1/rel_p^2 chromaticity against Genesis's 1/gz_hat, and the octupole edge kicks that
! Genesis does not have.
!
! g_max follows from aw: c*b_max = K*ku*m_e with K = aw (helical) or aw*sqrt(2)
! (planar), so g_max = K*ku/p0_mc.
!
! Input:
!   und     -- fel_und_struct: Undulator parameters.
!   beam    -- fel_beam_struct: The beam (reference energy).
!   sl      -- fel_slice_struct: One slice's packed particles.
!   delz    -- real(rp): Step length [m].
!   leading -- logical: True for the leading half-step (entrance fringe applies).
!
! Output:
!   sl      -- fel_slice_struct: Transverse coordinates advanced by delz.
!-

subroutine fel_transverse_track_bmad (und, beam, sl, delz, leading)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
real(rp) delz
logical leading

real(rp) p0_mc, g_max, k1x_loc, k1y_loc, rel_p, k1xx, k1yy, k3l, kz
real(rp) m2(2,2), dz_c(3), ddz_c(3)
integer ip

!

p0_mc = fel_p0_mc(beam)
kz = und%ku

if (und%helical) then
  g_max = und%aw * und%ku / p0_mc
  k1x_loc = -0.5_rp * g_max**2
  k1y_loc = k1x_loc
else
  g_max = und%aw * sqrt(2.0_rp) * und%ku / p0_mc
  k1x_loc = 0                              ! kx attribute is asserted zero at setup.
  k1y_loc = -0.5_rp * g_max**2
endif

do ip = 1, sl%n
  rel_p = 1 + sl%pz(ip)
  k1yy = k1y_loc / rel_p**2
  k1xx = k1x_loc / rel_p**2
  k3l = 2 * delz * k1yy                    ! Half of Bmad's per-step kick. See header.

  ! A tilted planar element focuses in its own wiggle frame: rotate in, apply the
  ! untilted map, rotate out. This is Bmad's tilt_coords composition. sin_t = 0 (every
  ! untilted element) skips both rotations.

  if (und%sin_t /= 0) call rot_xy (sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), und%cos_t, und%sin_t)

  if (leading) call octupole_kick ()

  call quad_mat2_calc (k1xx, delz, rel_p, m2, dz_c, ddz_c)
  call apply_mat2 (sl%x(ip), sl%px(ip))
  call quad_mat2_calc (k1yy, delz, rel_p, m2, dz_c, ddz_c)
  call apply_mat2 (sl%y(ip), sl%py(ip))

  if (.not. leading) call octupole_kick ()

  if (und%sin_t /= 0) call rot_xy (sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), und%cos_t, -und%sin_t)
enddo

!------------------------------------------------------------------------------
contains

!+
! Subroutine octupole_kick ()
!
! Routine to apply the wiggle-plane octupole-like kick of Bmad's wiggler body map
! (the k3l term of track_a_wiggler, transcribed to the packed arrays).
!-

subroutine octupole_kick ()
sl%py(ip) = sl%py(ip) + k3l * rel_p * kz**2 * sl%y(ip)**3 / 3
if (und%helical) then
  sl%px(ip) = sl%px(ip) + k3l * rel_p * kz**2 * sl%x(ip)**3 / 3
endif
end subroutine octupole_kick

!+
! Subroutine rot_xy (x, px, y, py, c, s)
!
! Routine to rotate one particle's transverse coordinates by the tilt angle whose
! cosine and sine are (c, s). Passing s > 0 rotates into the wiggle frame and s < 0
! rotates out of it.
!-

subroutine rot_xy (x, px, y, py, c, s)

real(rp) x, px, y, py, c, s, t1

t1 = c * x + s * y;   y = -s * x + c * y;   x = t1
t1 = c * px + s * py; py = -s * px + c * py; px = t1

end subroutine rot_xy

!+
! Subroutine apply_mat2 (v, vp)
!
! Routine to apply the cached 2x2 plane map m2 to one (position, momentum) pair.
!-

subroutine apply_mat2 (v, vp)
real(rp) v, vp, v1, v2
v1 = v; v2 = vp
v  = m2(1,1) * v1 + m2(1,2) * v2
vp = m2(2,1) * v1 + m2(2,2) * v2
end subroutine apply_mat2

end subroutine fel_transverse_track_bmad

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_apply_focus (delz, q_hat, x, px, gz_hat)
!
! Routine to apply one transverse plane's map over delz. Transcribed from
! TrackBeam::applyDrift for q_hat = 0, applyFQuad for q_hat > 0 and applyDQuad for
! q_hat < 0. The variables are those of fel_transverse_track (px = P_x/p0,
! gz_hat = gammaz/p0, q_hat = q/p0).
!
! Input:
!   delz   -- real(rp): Step length [m].
!   q_hat  -- real(rp): Scaled focusing strength (sign selects focus/defocus/drift).
!   x, px  -- real(rp): One transverse plane's position and momentum.
!   gz_hat -- real(rp): Longitudinal Lorentz factor scale.
!
! Output:
!   x, px  -- real(rp): Advanced through the focusing map.
!-

subroutine fel_apply_focus (delz, q_hat, x, px, gz_hat)

real(rp) delz, q_hat, x, px, gz_hat
real(rp) foc, omg, a1, a2, a3, xtmp

!

if (q_hat == 0) then
  x = x + px * delz / gz_hat
elseif (q_hat > 0) then
  foc = sqrt(q_hat/gz_hat)
  omg = foc * delz
  a1 = cos(omg)
  a2 = sin(omg)/foc
  a3 = -a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gz_hat
  px = a3 * xtmp * gz_hat + a1 * px
else
  foc = sqrt(-q_hat/gz_hat)
  omg = foc * delz
  a1 = cosh(omg)
  a2 = sinh(omg)/foc
  a3 = a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gz_hat
  px = a3 * xtmp * gz_hat + a1 * px
endif

end subroutine fel_apply_focus

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_advance (und, beam, sl, ff, delz, phi0_new, coll, is, err_flag)
!
! Routine to advance the longitudinal plane of every particle over delz. Transcribed
! from BeamSolver::advance.
!   1. Gather the field at (x, y) by bilinear interpolation from field slice ifld (the
!      rotated-record index, fel_field_index).
!   2. Form rpart = (fc(h)/ks)*faw*conj(E_h) per field of the set ff. The fundamental is
!      entry 1, and a single-entry set takes the scalar path.
!   3. Integrate (theta, gamma) by the verbatim RK4. The rpart set, px, py, faw and the
!      space-charge ez are held fixed through the stages. The harmonic phases h*theta
!      change per stage (fel_ode_multi).
! ez per particle is fel_shortrange_ez(ip) - long_esc(is)/m_electron. That is Genesis's
! ez = getEField(ip) + eloss (fel-physics.md sec-eom, sec-spacecharge).
!
! (theta, gamma) are derived at entry from the stored (z, pz) and written back at exit
! using phi0_new, the common phase at the end of this step. The module header says why
! this chart change is exact for RK4 and costs about 1 ulp for the energy.
!
! Input:
!   und      -- fel_und_struct: Undulator parameters (incl. source model).
!   beam     -- fel_beam_struct: The beam (reference energy, wavelength).
!   sl       -- fel_slice_struct: One slice's packed particles.
!   ff(:)    -- fel_field_struct: The field set (gathered once per particle per step).
!   delz     -- real(rp): Step length [m].
!   phi0_new -- real(rp): The field's phi0 after this step.
!   coll     -- fel_collective_struct: Collective terms (wake/space-charge Ez).
!   is       -- integer: This slice's beam index (for the field and collective lookups).
!
! Output:
!   err_flag -- logical: Set True when the collective solve is refused. False otherwise.
!   sl       -- fel_slice_struct: gamma/theta (and chart z) advanced by delz.
!-

subroutine fel_advance (und, beam, sl, ff, delz, phi0_new, coll, is, err_flag)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_collective_struct) coll
integer is
real(rp) delz, phi0_new
logical err_flag

real(rp) xks, xku, aw, rtmp, awloc, btpar, gamma, theta, beta, wx, wy, px_g, py_g, p_mc, p0_mc
real(rp) gz2, ez_ip, esc_loss
real(rp), allocatable :: ez(:)
real(rp) rtmp_h(size(ff)), rharm(size(ff))
real(rp) gridmax
integer ngrid_arr(3)
complex(rp) cpart, rpart
complex(rp) rpart_h(size(ff))
integer ip, ix, iy, ifld, io, nf
integer ifld_h(size(ff))
logical on_grid

!

nf = size(ff)
wf => ff(1)%wf
ifld = fel_field_index(ff(1)%slip, is, size(wf%Ex, 3))

p0_mc = fel_p0_mc(beam)
xks = twopi / wf%wavelength
xku = und%ku
aw = und%aw

! Space charge of this slice, computed once per step and held through the RK stages
! (BeamSolver::advance's order). The short-range solve uses call-local storage, so it
! is safe in the parallel slice loop. The caller refreshed long_esc serially.

err_flag = .false.
allocate (ez(max(1, sl%n)))
gz2 = fel_gamma0(beam)**2 / (1 + aw**2)
call fel_shortrange_ez (coll%efield, beam, sl, gz2, xks, ez, err_flag)
if (err_flag) return       ! No force for this slice, so nothing to integrate.
esc_loss = 0
if (allocated(coll%long_esc)) esc_loss = -coll%long_esc(is) / m_electron

! Coupling coefficient for the energy exchange: fc/(sqrt(2)*m_electron), the V/m form of
! Genesis's fc/ks (see the module header for the unit relation). rpart then has units of
! 1/m and dgamma/dz is per meter directly.

rtmp = fel_und_coupling(und, 1) / (sqrt(2.0_rp) * m_electron)

! Grid bounds are loop-invariant, so they are computed once per call. The expression is
! the one fel_grid_weights uses, so the result is the same bit for bit.

ngrid_arr = wavefront_shape(wf)
gridmax = (ngrid_arr(1) - 1) * wf%dx / 2

if (nf == 1) then

  ! Single field: the scalar path.

  do ip = 1, sl%n
    p_mc = p0_mc * (1 + sl%pz(ip))         ! gamma and beta = P/gamma share one sqrt,
    gamma = sqrt(p_mc**2 + 1)              ! bit-identical to fel_gamma_of/fel_beta_of.
    beta = p_mc / gamma
    theta = beam%phi0 + xks * sl%z(ip) / beta

    awloc = faw(und, sl%x(ip), sl%y(ip))
    px_g = sl%px(ip) * p0_mc
    py_g = sl%py(ip) * p0_mc
    btpar = 1 + px_g*px_g + py_g*py_g + aw*aw*awloc*awloc

    call fel_grid_weights_pre (gridmax, wf%dx, wf%dy, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
    if (on_grid) then
      cpart =         wf%Ex(ix,   iy,   ifld) * wx * wy
      cpart = cpart + wf%Ex(ix+1, iy,   ifld) * (1-wx) * wy
      cpart = cpart + wf%Ex(ix,   iy+1, ifld) * wx * (1-wy)
      cpart = cpart + wf%Ex(ix+1, iy+1, ifld) * (1-wx) * (1-wy)

      ! Two live polarizations: the element couples to E_eff = conj(pol).E
      ! (fel-physics.md sec-field). With one polarization this branch does not run.

      if (allocated(wf%Ey)) then
        cpart = conjg(und%pol(1)) * cpart
        cpart = cpart + conjg(und%pol(2)) * (wf%Ey(ix,   iy,   ifld) * wx * wy &
                                           + wf%Ey(ix+1, iy,   ifld) * (1-wx) * wy &
                                           + wf%Ey(ix,   iy+1, ifld) * wx * (1-wy) &
                                           + wf%Ey(ix+1, iy+1, ifld) * (1-wx) * (1-wy))
      endif
      rpart = rtmp * awloc * conjg(cpart)
    else
      rpart = 0
    endif

    ez_ip = ez(ip) + esc_loss     ! Short range plus the long-range loss (sec-spacecharge).
    call fel_runge_kutta (delz, xks, xku, btpar, rpart, ez_ip, gamma, theta)

    ! Back to the stored chart. pz comes from gamma (the subtraction is exact once p_mc
    ! is formed). z comes from tau = (phi0_new - theta)/ks with the updated beta.

    p_mc = sqrt(gamma**2 - 1)
    sl%pz(ip) = (p_mc - p0_mc) / p0_mc
    beta = p_mc / gamma
    sl%z(ip) = -beta * (phi0_new - theta) / xks
  enddo

else

  ! The field set: every field acts on the beam at once, as in BeamSolver::advance
  ! (BeamSolver.cpp:24-30,150). The per-field coupling fc(h) and record index are set
  ! outside the particle loop. The per-field rpart (coupling times that field's
  ! interpolated amplitude, no phase) is held through the RK stages. fel_ode_multi
  ! applies the harmonic phase h*theta per stage.

  do io = 1, nf
    rtmp_h(io) = fel_und_coupling(und, ff(io)%harm) / (sqrt(2.0_rp) * m_electron)
    rharm(io) = ff(io)%harm
    ifld_h(io) = fel_field_index(ff(io)%slip, is, size(ff(io)%wf%Ex, 3))
  enddo

  do ip = 1, sl%n
    p_mc = p0_mc * (1 + sl%pz(ip))
    gamma = sqrt(p_mc**2 + 1)
    beta = p_mc / gamma
    theta = beam%phi0 + xks * sl%z(ip) / beta

    awloc = faw(und, sl%x(ip), sl%y(ip))
    px_g = sl%px(ip) * p0_mc
    py_g = sl%py(ip) * p0_mc
    btpar = 1 + px_g*px_g + py_g*py_g + aw*aw*awloc*awloc

    do io = 1, nf
      call fel_grid_weights (ff(io)%wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
      if (on_grid) then
        cpart =         ff(io)%wf%Ex(ix,   iy,   ifld_h(io)) * wx * wy
        cpart = cpart + ff(io)%wf%Ex(ix+1, iy,   ifld_h(io)) * (1-wx) * wy
        cpart = cpart + ff(io)%wf%Ex(ix,   iy+1, ifld_h(io)) * wx * (1-wy)
        cpart = cpart + ff(io)%wf%Ex(ix+1, iy+1, ifld_h(io)) * (1-wx) * (1-wy)
        rpart_h(io) = rtmp_h(io) * awloc * conjg(cpart)
      else
        rpart_h(io) = 0
      endif
    enddo

    ez_ip = ez(ip) + esc_loss
    call fel_runge_kutta_multi (delz, xks, xku, btpar, rpart_h, rharm, ez_ip, gamma, theta)

    p_mc = sqrt(gamma**2 - 1)
    sl%pz(ip) = (p_mc - p0_mc) / p0_mc
    beta = p_mc / gamma
    sl%z(ip) = -beta * (phi0_new - theta) / xks
  enddo

endif

end subroutine fel_advance

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_runge_kutta (delz, xks, xku, btpar, rpart, ez, gamma, theta)
!
! Routine to take one RK4 step. The stage bookkeeping is BeamSolver::RungeKutta,
! verbatim (fel-physics.md sec-eom). ez is held fixed through the stages, as Genesis
! holds it.
!
! Note: The in-place stage algebra is kept exactly for bit identity. Do not rewrite it.
!
! Input:
!   delz   -- real(rp): Step length [m].
!   xks    -- real(rp): Radiation wavenumber [1/m].
!   xku    -- real(rp): Undulator wavenumber [1/m].
!   btpar  -- real(rp): Parallel velocity beta_par.
!   rpart  -- complex(rp): The particle's gathered field phasor times coupling.
!   ez     -- real(rp): Collective longitudinal field at the particle [eV/m scale].
!   gamma  -- real(rp): Particle Lorentz factor.
!   theta  -- real(rp): Ponderomotive phase.
!
! Output:
!   gamma  -- real(rp): Advanced by the RK4 step.
!   theta  -- real(rp): Advanced by the RK4 step.
!-

subroutine fel_runge_kutta (delz, xks, xku, btpar, rpart, ez, gamma, theta)

real(rp) delz, xks, xku, btpar, ez, gamma, theta
complex(rp) rpart
real(rp) k2gg, k2pp, k3gg, k3pp, stpz

! first step

k2gg = 0
k2pp = 0

call fel_ode (gamma, theta, xks, xku, btpar, rpart, ez, k2gg, k2pp)

! second step

stpz = 0.5_rp * delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k2gg
k3pp = k2pp

k2gg = 0
k2pp = 0

call fel_ode (gamma, theta, xks, xku, btpar, rpart, ez, k2gg, k2pp)

! third step

gamma = gamma + stpz * (k2gg - k3gg)
theta = theta + stpz * (k2pp - k3pp)

k3gg = k3gg / 6
k3pp = k3pp / 6

k2gg = k2gg * (-0.5_rp)
k2pp = k2pp * (-0.5_rp)

call fel_ode (gamma, theta, xks, xku, btpar, rpart, ez, k2gg, k2pp)

! fourth step

stpz = delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k3gg - k2gg
k3pp = k3pp - k2pp

k2gg = k2gg * 2
k2pp = k2pp * 2

call fel_ode (gamma, theta, xks, xku, btpar, rpart, ez, k2gg, k2pp)

gamma = gamma + stpz * (k3gg + k2gg / 6.0_rp)
theta = theta + stpz * (k3pp + k2pp / 6.0_rp)

end subroutine fel_runge_kutta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_ode (tgam, tthet, xks, xku, btpar, rpart, ez, k2gg, k2pp)
!
! Routine to evaluate the longitudinal equations of motion, fundamental only.
! Transcribed from BeamSolver::ODE (fel-physics.md sec-eom). ez is the per-particle
! space-charge term in m_e c^2 per meter (short-range harmonics plus the long-range
! loss). ez is zero when space charge is off, and subtracting zero is exact.
!
! Input:
!   tgam, tthet -- real(rp): Stage values of gamma and theta.
!   xks, xku    -- real(rp): Radiation and undulator wavenumbers [1/m].
!   btpar       -- real(rp): Parallel velocity beta_par.
!   rpart       -- complex(rp): Field phasor times coupling.
!   ez          -- real(rp): Collective longitudinal field.
!
! Output:
!   k2gg, k2pp  -- real(rp): The stage derivatives of gamma and theta.
!-

subroutine fel_ode (tgam, tthet, xks, xku, btpar, rpart, ez, k2gg, k2pp)

real(rp) tgam, tthet, xks, xku, btpar, ez, k2gg, k2pp
complex(rp) rpart, ctmp
real(rp) ztemp1, btper0, btpar0, s_t, c_t

!

ztemp1 = -2.0_rp / xks
call fel_sincos (tthet, s_t, c_t)
ctmp = rpart * cmplx(c_t, -s_t, rp)

btper0 = btpar + ztemp1 * real(ctmp, rp)
btpar0 = sqrt(1 - btper0 / (tgam * tgam))

k2pp = k2pp + xks * (1 - 1/btpar0) + xku
k2gg = k2gg + aimag(ctmp) / btpar0 / tgam - ez

end subroutine fel_ode

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_runge_kutta_multi (delz, xks, xku, btpar, rpart, rharm, ez, gamma, theta)
!
! Routine to take one RK4 step over a field set. The stage bookkeeping is that of
! fel_runge_kutta (BeamSolver::RungeKutta). The slope sums every field's contribution
! at that field's harmonic phase (fel_ode_multi). A single-field walk does not call
! this, since fel_advance keeps the scalar path.
!
! Input:
!   delz     -- real(rp): Step length [m].
!   xks, xku -- real(rp): Fundamental radiation and undulator wavenumbers [1/m].
!   btpar    -- real(rp): Parallel velocity beta_par.
!   rpart(:) -- complex(rp): Per-field gathered phasors times couplings.
!   rharm(:) -- real(rp): The harmonic number of each field, as a real factor.
!   ez       -- real(rp): Collective longitudinal field.
!   gamma    -- real(rp): Particle Lorentz factor.
!   theta    -- real(rp): Fundamental ponderomotive phase.
!
! Output:
!   gamma    -- real(rp): Advanced by the RK4 step.
!   theta    -- real(rp): Advanced by the RK4 step.
!-

subroutine fel_runge_kutta_multi (delz, xks, xku, btpar, rpart, rharm, ez, gamma, theta)

real(rp) delz, xks, xku, btpar, ez, gamma, theta
real(rp) rharm(:)
complex(rp) rpart(:)
real(rp) k2gg, k2pp, k3gg, k3pp, stpz

! first step

k2gg = 0
k2pp = 0

call fel_ode_multi (gamma, theta, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)

! second step

stpz = 0.5_rp * delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k2gg
k3pp = k2pp

k2gg = 0
k2pp = 0

call fel_ode_multi (gamma, theta, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)

! third step

gamma = gamma + stpz * (k2gg - k3gg)
theta = theta + stpz * (k2pp - k3pp)

k3gg = k3gg / 6
k3pp = k3pp / 6

k2gg = k2gg * (-0.5_rp)
k2pp = k2pp * (-0.5_rp)

call fel_ode_multi (gamma, theta, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)

! fourth step

stpz = delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k3gg - k2gg
k3pp = k3pp - k2pp

k2gg = k2gg * 2
k2pp = k2pp * 2

call fel_ode_multi (gamma, theta, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)

gamma = gamma + stpz * (k3gg + k2gg / 6.0_rp)
theta = theta + stpz * (k3pp + k2pp / 6.0_rp)

end subroutine fel_runge_kutta_multi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_ode_multi (tgam, tthet, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)
!
! Routine to evaluate the longitudinal equations of motion over a field set.
! Transcribed from BeamSolver::ODE (BeamSolver.cpp:150). ctmp sums
!   rpart(i) * exp(-i * rharm(i) * theta)
! over the fields at the current stage theta. The couplings are held fixed through the
! stages and the harmonic phases change. xks is the fundamental wavenumber. The theta
! equation is the fundamental's, and harmonics enter the slope only through ctmp.
!
! Input:
!   tgam, tthet -- real(rp): Stage values of gamma and the fundamental theta.
!   xks, xku    -- real(rp): Fundamental radiation and undulator wavenumbers [1/m].
!   btpar       -- real(rp): Parallel velocity beta_par.
!   rpart(:)    -- complex(rp): Per-field phasors times couplings.
!   rharm(:)    -- real(rp): Harmonic numbers as real factors (phase h*theta).
!   ez          -- real(rp): Collective longitudinal field.
!
! Output:
!   k2gg, k2pp  -- real(rp): The stage derivatives of gamma and theta.
!-

subroutine fel_ode_multi (tgam, tthet, xks, xku, btpar, rpart, rharm, ez, k2gg, k2pp)

real(rp) tgam, tthet, xks, xku, btpar, ez, k2gg, k2pp
real(rp) rharm(:)
complex(rp) rpart(:), ctmp
real(rp) ztemp1, btper0, btpar0, s_t, c_t
integer i

!

ztemp1 = -2.0_rp / xks
ctmp = 0
do i = 1, size(rpart)
  call fel_sincos (rharm(i) * tthet, s_t, c_t)
  ctmp = ctmp + rpart(i) * cmplx(c_t, -s_t, rp)
enddo

btper0 = btpar + ztemp1 * real(ctmp, rp)
btpar0 = sqrt(1 - btper0 / (tgam * tgam))

k2pp = k2pp + xks * (1 - 1/btpar0) + xku
k2gg = k2gg + aimag(ctmp) / btpar0 / tgam - ez

end subroutine fel_ode_multi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_grid_weights (wf, x, y, ix, iy, wx, wy, on_grid)
!
! Routine to map a particle position to its lower-left grid cell corner and bilinear
! weights. Transcribed from Field::getLLGridpoint. wx is the weight of the lower x point.
! A particle outside |x|,|y| < gridmax neither feels the field nor radiates into it.
!
! Input:
!   wf      -- wavefront_struct: The field (grid geometry).
!   x, y    -- real(rp): Transverse position [m].
!
! Output:
!   ix, iy  -- integer: Lower-left grid cell of the bilinear stencil.
!   wx, wy  -- real(rp): Fractional weights within the cell.
!   on_grid -- logical: False when the particle is outside the grid.
!-

subroutine fel_grid_weights (wf, x, y, ix, iy, wx, wy, on_grid)

type (wavefront_struct) wf
real(rp) x, y, wx, wy, gridmax
integer ix, iy, ngrid_arr(3), ngrid
logical on_grid

!

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
gridmax = (ngrid - 1) * wf%dx / 2       ! Genesis's gridmax. dgrid = 2*gridmax/(ngrid-1).

if (x > -gridmax .and. x < gridmax .and. y > -gridmax .and. y < gridmax) then
  wx = (x + gridmax) / wf%dx
  wy = (y + gridmax) / wf%dy
  ix = int(floor(wx))
  iy = int(floor(wy))
  wx = 1 + floor(wx) - wx
  wy = 1 + floor(wy) - wy
  ix = ix + 1                           ! To Fortran 1-based.
  iy = iy + 1
  on_grid = .true.
else
  on_grid = .false.
endif

end subroutine fel_grid_weights

!------------------------------------------------------------------------------
!+
! Subroutine fel_grid_weights_pre (gridmax, dx, dy, x, y, ix, iy, wx, wy, on_grid)
!
! Routine to do what fel_grid_weights does, with the grid bounds precomputed by the
! caller. The particle loops call this per particle, and wavefront_shape and gridmax
! are loop-invariant. The arithmetic is that of fel_grid_weights, so the results agree
! bit for bit when the caller passes gridmax = (ngrid - 1) * wf%dx / 2.
!
! Input:
!   gridmax -- real(rp): Grid half-width [m].
!   dx, dy  -- real(rp): Grid spacings [m].
!   x, y    -- real(rp): Transverse position [m].
!
! Output:
!   ix, iy  -- integer: Lower-left grid cell of the bilinear stencil.
!   wx, wy  -- real(rp): Fractional weights within the cell.
!   on_grid -- logical: False when the particle is outside the grid.
!-

subroutine fel_grid_weights_pre (gridmax, dx, dy, x, y, ix, iy, wx, wy, on_grid)

real(rp) gridmax, dx, dy, x, y, wx, wy
integer ix, iy
logical on_grid

!

if (x > -gridmax .and. x < gridmax .and. y > -gridmax .and. y < gridmax) then
  wx = (x + gridmax) / dx
  wy = (y + gridmax) / dy
  ix = int(floor(wx))
  iy = int(floor(wy))
  wx = 1 + floor(wx) - wx
  wy = 1 + floor(wy) - wy
  ix = ix + 1                           ! To Fortran 1-based.
  iy = iy + 1
  on_grid = .true.
else
  on_grid = .false.
endif

end subroutine fel_grid_weights_pre

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_coherent_prep (und, beam, sl, xks1, coh)
!
! Routine to compute one slice's coherent-source summary (fel-physics.md
! sec-coherent-source). Implemented from Tanaka, PRAB 27, 030703 (2024). Computed:
!   The source phasor S, as the deposit would form it (same theta, same part factors,
!     post-advance particles). The normalization is SUM crsource = S.
!   The charge-weighted centroid and central second moments. These extend the paper to
!     offset, mismatched and tilted beams.
!   The Laguerre-Gauss order-0/1 sums (Eqs 19-20, charge-normalized weights):
!     g_q = SUM w_j Psi_q((x-x0)/sx, (y-y0)/sy) e^{-i theta_j}
!     Psi_q(u,v) = (1/sqrt(pi)) L_q(u^2+v^2) exp(-(u^2+v^2)/2)
!   The guard metric (worst-plane excess kurtosis).
!
! Note: Scale factors constant across the slice (scl_w) are *not* folded into S here.
! The caller applies them, so S matches the deposit's sum with part = sqrt(faw2)*w/gam.
!
! Input:
!   und  -- fel_und_struct: Undulator parameters (source model, tilt frame).
!   beam -- fel_beam_struct: The beam (reference energy, wavelength).
!   sl   -- fel_slice_struct: One slice's packed particles.
!   xks1 -- real(rp): Fundamental radiation wavenumber [1/m].
!
! Output:
!   coh  -- fel_coherent_struct: The slice's phasor-weighted moments, LG reduction
!             sums, Gaussianity statistics and validity flag.
!-

subroutine fel_coherent_prep (und, beam, sl, xks1, coh)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (fel_coherent_struct) coh
real(rp) xks1

real(rp) p0_mc, p_mc, gam, beta, theta, part, w, wsum, w2sum, s_t, c_t
real(rp) sx, sy, u, v, r2, e2, psi0, psi1, x4, y4
complex(rp) g0, g1, bsum, ephase
integer ip

!

coh = fel_coherent_struct()
if (sl%n < 8) return
p0_mc = fel_p0_mc(beam)

! Pass 1: the source phasor (as the deposit forms it) and the charge moments.

wsum = 0
w2sum = 0
do ip = 1, sl%n
  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam
  theta = beam%phi0 + xks1 * sl%z(ip) / beta
  part = sqrt(faw2(und, sl%x(ip), sl%y(ip))) * sl%weight(ip) / gam
  call fel_sincos (theta, s_t, c_t)
  coh%S = coh%S + cmplx(s_t, c_t, rp) * part
  w = sl%weight(ip)
  wsum = wsum + w
  w2sum = w2sum + w*w
  coh%x0 = coh%x0 + w * sl%x(ip)
  coh%y0 = coh%y0 + w * sl%y(ip)
enddo
if (wsum <= 0) return
coh%wsum = wsum
coh%x0 = coh%x0 / wsum
coh%y0 = coh%y0 / wsum

do ip = 1, sl%n
  w = sl%weight(ip)
  coh%sxx = coh%sxx + w * (sl%x(ip) - coh%x0)**2
  coh%sxy = coh%sxy + w * (sl%x(ip) - coh%x0) * (sl%y(ip) - coh%y0)
  coh%syy = coh%syy + w * (sl%y(ip) - coh%y0)**2
enddo
coh%sxx = coh%sxx / wsum;  coh%sxy = coh%sxy / wsum;  coh%syy = coh%syy / wsum
if (coh%sxx <= 0 .or. coh%syy <= 0 .or. coh%sxx * coh%syy - coh%sxy**2 <= 0) return
sx = sqrt(coh%sxx)
sy = sqrt(coh%syy)

! Pass 2: the LG sums (Eqs 19-20, principal widths, centered) and the guard's fourth
! moments. B has the same charge weighting, so the Eq 26 ratio is consistent.

g0 = 0;  g1 = 0;  bsum = 0;  x4 = 0;  y4 = 0
do ip = 1, sl%n
  p_mc = p0_mc * (1 + sl%pz(ip))
  gam = sqrt(p_mc**2 + 1)
  beta = p_mc / gam
  theta = beam%phi0 + xks1 * sl%z(ip) / beta
  call fel_sincos (theta, s_t, c_t)
  ephase = cmplx(c_t, -s_t, rp)                 ! e^{-i theta}
  w = sl%weight(ip) / wsum
  u = (sl%x(ip) - coh%x0) / sx
  v = (sl%y(ip) - coh%y0) / sy
  r2 = u*u + v*v
  e2 = exp(-r2 / 2) / sqrt(pi)
  psi0 = e2
  psi1 = (1 - r2) * e2
  g0 = g0 + w * psi0 * ephase
  g1 = g1 + w * psi1 * ephase
  bsum = bsum + w * ephase
  x4 = x4 + w * u**4
  y4 = y4 + w * v**4
enddo
coh%g2sum = abs(g0)**2 + abs(g1)**2
coh%b2 = abs(bsum)**2
coh%kurt = max(abs(x4 - 3), abs(y4 - 3))        ! Excess kurtosis, worst plane.
coh%m_ind = max(1.0_rp, sl%m_ind)      ! Independent transverse samples, counted at load.
coh%ok = .true.

end subroutine fel_coherent_prep

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_step (und, beam, sl, wf, ifld, delz, harm, xks1, err_flag, coh, kappa)
!
! Routine to advance field slice ifld (the rotated-record index, fel_field_index) one
! step. The source is built from the already pushed particles of beam slice sl. The
! field is propagated by exp(K2 delz) in transverse Fourier space and the source is
! added. Transcribed from FieldSolverFFT::advance and FFT.
!
! The source scale, weighted, in V/m (fel-physics.md sec-field derives it from
! Genesis's internal-unit form):
!
!   scl_w = fc * Z0 * sqrt(2) * c * delz / (4 * dgrid^2 * slice_spacing)
!
! The rest energy and the wavenumber cancel. Per particle
! part = sqrt(faw2(x,y))*scl_w*w_j/gamma, deposited as (sin theta + i cos theta)*part
! with the bilinear weights. Without the source filter the source is added times 2 in
! real space after the transform pair.
!
! Input:
!   und      -- fel_und_struct: Undulator parameters (incl. source model).
!   beam     -- fel_beam_struct: The beam (reference energy, wavelength).
!   sl       -- fel_slice_struct: One slice's packed particles.
!   wf       -- wavefront_struct: One field of the set.
!   ifld     -- integer: Field-record index of the slice.
!   delz     -- real(rp): Step length [m].
!   harm     -- integer: This field's harmonic number.
!   xks1     -- real(rp): Fundamental radiation wavenumber [1/m].
!   coh      -- fel_coherent_struct, optional: The slice's coherent-source moments
!                 (source_model = 'coherent' only).
!   kappa    -- real(rp), optional: The fitted coherent width ratio.
!
! Output:
!   wf       -- wavefront_struct: Source deposited into the slice's field record.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_field_step (und, beam, sl, wf, ifld, delz, harm, xks1, err_flag, coh, kappa)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (wavefront_struct), target :: wf
integer ifld, harm
real(rp) delz, xks1
logical err_flag
type (fel_coherent_struct), optional :: coh    ! Coherent source (sec-coherent-source).
real(rp), optional :: kappa                    !   Both present when und%source_model says so.

real(rp) xks, dgrid, scl_w, part, theta, beta, wx, wy, gam, p_mc, p0_mc, s_t, c_t
complex(rp) cpart
complex(rp), allocatable :: crsource(:,:)    ! Per-call source accumulator, thread safe.
integer ip, ix, iy, ngrid_arr(3), ngrid, ik
logical on_grid, err
character(*), parameter :: r_name = 'fel_field_step'

!

err_flag = .true.

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
xks = twopi / wf%wavelength
dgrid = wf%dx
p0_mc = fel_p0_mc(beam)

! The kernel is read-only here, since a rebuild would race with concurrent slices. The
! caller builds it serially (fel_field_kernel_init). A miss is a caller bug.

ik = fel_kernel_index(ngrid, dgrid, xks, delz)
if (ik == 0) then
  call out_io (s_error$, r_name, 'FIELD KERNEL NOT INITIALIZED FOR THIS GRID AND STEP. ' // &
                                 'CALL fel_field_kernel_init FIRST (SERIALLY).')
  return
endif

allocate (crsource(ngrid, ngrid))
crsource = 0

scl_w = fel_und_coupling(und, harm) * (mu_0_vac * c_light) * sqrt(2.0_rp) * c_light * delz
scl_w = scl_w / (4 * dgrid * dgrid * beam%slice_spacing)

if (scl_w /= 0 .and. und%source_model == fel_source_coherent$) then

  ! The coherent-Gaussian source (fel-physics.md sec-coherent-source). The slice's
  ! whole phasor S is deposited as one analytic Gaussian at the charge centroid, with
  ! covariance kappa^2 * (second-moment matrix). S carries the physical shot noise
  ! through B(s), which the Fawley loading sets to <|B|^2> N_lambda = 1. The discrete
  ! sum is normalized to exactly S, as the deposit has SUM crsource = SUM cpart. So
  ! edge truncation loses nothing.

  call coherent_gaussian_source ()

elseif (scl_w /= 0) then
  do ip = 1, sl%n
    call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
    if (.not. on_grid) cycle

    p_mc = p0_mc * (1 + sl%pz(ip))
    gam = sqrt(p_mc**2 + 1)
    beta = p_mc / gam

    ! The fundamental ponderomotive phase, scaled to harm*theta for a harmonic source
    ! (FieldSolver.cpp's theta = harm * particle.theta). The caller passes xks1, the
    ! fundamental wavenumber, so a harmonic field never rebuilds it from the rounding
    ! of its own wavelength.

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

! The source filter (fel-physics.md sec-source-filter). The source gets a forward
! transform of its own and the sigmoid multiplies it there. It is added to the field in
! Fourier space before the one inverse transform. FieldSolverFFT::FFT carries both
! branches for the same reason. Without the filter the transform pair is linear, and
! IFFT(FFT(crsource))/ngrid^2 is crsource again. So the unfiltered branch adds the
! source in real space, one FFT less per field per step. That branch is kept whole, so
! results with the filter off do not change.

if (und%filter%on) then

  if (.not. allocated(fel_kernels(ik)%sigmoid)) then
    call out_io (s_error$, r_name, 'SOURCE FILTER IS ON BUT ITS SIGMOID WAS NOT BUILT FOR ' // &
                                   'THIS GRID AND STEP.', 'PLEASE REPORT THIS!')
    return
  endif

  ! The filter's self-test switch puts the sigmoid on the propagated field instead of the
  ! source, the error check_source_filter.py must catch. That also suppresses wide
  ! angles, so the run completes with a plausible power curve. Only the Genesis4
  ! comparison tells the two apart.

  if (.not. und%filter%mutate) then
    call wavefront_fft2 (crsource, wf_fft_forward$, err);  if (err) return
    crsource = crsource * fel_kernels(ik)%sigmoid
  else
    call wavefront_fft2 (crsource, wf_fft_forward$, err);  if (err) return
  endif

  call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_forward$, err);  if (err) return
  if (und%filter%mutate) wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%sigmoid

  if (allocated(wf%Ey)) then
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%exp_k2 + 2 * und%pol(1) * crsource
    call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp)
    call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_forward$, err);  if (err) return
    if (und%filter%mutate) wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * fel_kernels(ik)%sigmoid
    wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * fel_kernels(ik)%exp_k2 + 2 * und%pol(2) * crsource
    call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_backward$, err);  if (err) return
    wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) / real(ngrid*ngrid, rp)
  else
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%exp_k2 + 2 * crsource
    call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp)
  endif

else

  ! Propagate and add: FFT, multiply by the cached exp(K2 delz), inverse FFT, normalize,
  ! and add 2*crsource in real space.

  call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_forward$, err);  if (err) return
  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%exp_k2
  call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return

  if (allocated(wf%Ey)) then

    ! Two live polarizations: the element's source lands as pol * src on the pair, the
    ! exact dual of the kick's conj(pol).E. Ey diffracts the same way as Ex.

    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * und%pol(1) * crsource
    call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_forward$, err);  if (err) return
    wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * fel_kernels(ik)%exp_k2
    call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_backward$, err);  if (err) return
    wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * und%pol(2) * crsource
  else
    wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * crsource
  endif

endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine coherent_gaussian_source ()
!
! Routine to deposit the coherent-Gaussian source (fel-physics.md sec-coherent-source).
! The slice's bunching phasor is deposited as a Gaussian of the charge second moments
! scaled by kappa. It is normalized so the deposited sum equals the particle deposit's.
!-

subroutine coherent_gaussian_source ()

real(rp) det, a11, a12, a22, xg, yg, dx_c, dy_c, q, gsum, gridmax
real(rp), allocatable :: gk(:,:)
complex(rp) s_scaled
integer i, j

! Covariance kappa^2 * [[sxx, sxy],[sxy, syy]], and its inverse for the quad form.

if (.not. (present(coh) .and. present(kappa))) return
s_scaled = scl_w * coh%S
if (.not. coh%ok .or. abs(s_scaled) == 0) return

det = (kappa**2)**2 * (coh%sxx * coh%syy - coh%sxy**2)
a11 =  kappa**2 * coh%syy / det
a22 =  kappa**2 * coh%sxx / det
a12 = -kappa**2 * coh%sxy / det

allocate (gk(ngrid, ngrid))
gridmax = (ngrid - 1) * dgrid / 2
gsum = 0
do j = 1, ngrid
  yg = (j - 1) * dgrid - gridmax
  dy_c = yg - coh%y0
  do i = 1, ngrid
    xg = (i - 1) * dgrid - gridmax
    dx_c = xg - coh%x0
    q = 0.5_rp * (a11 * dx_c*dx_c + 2 * a12 * dx_c*dy_c + a22 * dy_c*dy_c)
    if (q < 60.0_rp) then
      gk(i, j) = exp(-q)
      gsum = gsum + gk(i, j)
    else
      gk(i, j) = 0
    endif
  enddo
enddo
if (gsum <= 0) return
crsource = crsource + (s_scaled / gsum) * gk

end subroutine coherent_gaussian_source

end subroutine fel_field_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_field_kernel_exp_k2 (ks) result (p)
!
! Routine to point at the cached step propagator for this wavelength. A caller outside
! this module can then round it into another precision. The result is null when no
! entry has been built, that is, when the caller did not run fel_field_kernel_init first.
!-

function fel_field_kernel_exp_k2 (ks) result (p)

complex(rp), pointer :: p(:,:)
real(rp) ks
integer ik

!

p => null()
if (.not. allocated(fel_kernels)) return
do ik = 1, size(fel_kernels)
  if (fel_kernels(ik)%ks == ks .and. allocated(fel_kernels(ik)%exp_k2)) then
    p => fel_kernels(ik)%exp_k2
    return
  endif
enddo

end function fel_field_kernel_exp_k2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_kernel_init (ngrid, dgrid, ks, dz, filter)
!
! Routine to build the K2 kernel when the grid or the step changes, as
! FieldSolverFFT::init builds it:
!   dk = twopi/(ngrid*dgrid), integer offsets from the grid center,
!   fftshift index mapping, K2 = -i*(dx^2+dy^2)*dk^2/(2*ks)
! The step propagator exp(K2*dz) is built too. K2 and dz are constant across each
! caller's parallel loop. So the ngrid^2 complex exponentials are computed once per
! (grid, ks, dz) and not once per slice per step.
!
! Note: Call this serially. fel_track_und_step and fel_track_interlude_genesis call it
! before their parallel slice loops. fel_field_step only reads the kernel and errors on
! a miss, so nothing writes module state concurrently.
!
! Input:
!   ngrid  -- integer: Grid points per side.
!   dgrid  -- real(rp): Grid spacing [m].
!   ks     -- real(rp): Radiation wavenumber [1/m].
!   dz     -- real(rp): Step length [m].
!   filter -- fel_source_filter_struct, optional: The source-term angular filter. The
!               entry gains its sigmoid table when this is present and on. Absent or
!               off, no table is built and fel_field_step adds the source unfiltered.
!
! Output:
!   None directly. The module kernel cache (fel_kernels) gains an entry, and the FFTW
!   plans are made serially, so the parallel loops only execute them.
!-

subroutine fel_field_kernel_init (ngrid, dgrid, ks, dz, filter)

integer ngrid
real(rp) dgrid, ks, dz
type (fel_source_filter_struct), optional :: filter
real(rp) dk, shift, dx, dy, xf, yf, r
integer mx, my
type (fel_kernel_struct), allocatable :: grow(:)
type (fel_kernel_struct), pointer :: kn
integer ix, iy, iix, iiy, ik
logical err, want_sigmoid

!

! Make every thread's FFT plans from this serial context. No FFTW planner then runs
! concurrently with transforms in the parallel slice loops that follow (see
! wavefront_fft2_plan_threads). This runs even on a kernel-cache hit, since the thread
! team may have grown since the plans were made. A planning failure is reported there,
! and the first transform then returns an error.

call wavefront_fft2_plan_threads (ngrid, ngrid, err)

! One cache entry per wavelength. Find this ks's entry and rebuild it when the grid or
! the step changed, or append a new one. A single-field run uses entry 1 only.

if (.not. allocated(fel_kernels)) allocate (fel_kernels(0))

want_sigmoid = .false.
if (present(filter)) want_sigmoid = filter%on

do ik = 1, size(fel_kernels)
  if (fel_kernels(ik)%ks /= ks) cycle

  ! A hit must also carry the table this caller needs. An interlude reaches the same
  ! four keys with no filter. An entry it built first would leave the undulator step's
  ! filter with no sigmoid. Rebuilding here is safe since this routine runs serially.

  if (fel_kernels(ik)%ngrid == ngrid .and. fel_kernels(ik)%dgrid == dgrid .and. &
      fel_kernels(ik)%dz == dz .and. &
      (allocated(fel_kernels(ik)%sigmoid) .or. .not. want_sigmoid)) return
  exit
enddo

if (ik > size(fel_kernels)) then
  call move_alloc (fel_kernels, grow)
  allocate (fel_kernels(size(grow) + 1))
  fel_kernels(1:size(grow)) = grow
  deallocate (grow)
endif

kn => fel_kernels(ik)
if (allocated(kn%k2)) deallocate(kn%k2, kn%exp_k2)
allocate (kn%k2(ngrid, ngrid), kn%exp_k2(ngrid, ngrid))

dk = twopi / (ngrid * dgrid)

! The transverse wavenumbers a discrete transform of ngrid points carries: 0, 1, ... up
! to just under half of ngrid, then the negatives, each a whole multiple of dk. Genesis
! builds the table as (ix - (ngrid-1)/2) * dk and rotates it into transform order
! (FieldSolverFFT.cpp, FieldSolverFFT::init). For an odd ngrid that offset is a whole
! number and the two tables are the same. For an even ngrid the offset is a half.
! Genesis then carries every mode at (m + 1/2) * dk, and the constant mode at dk/2
! instead of zero. The error in k^2 is m + 1/4, which grows with the mode. The
! propagated field then picks up a phase that no beam produced. So the table here
! differs from Genesis's on even grids.

do iy = 0, ngrid-1
  my = iy;  if (2*my > ngrid-1) my = my - ngrid
  dy = my
  do ix = 0, ngrid-1
    mx = ix;  if (2*mx > ngrid-1) mx = mx - ngrid
    dx = mx
    kn%k2(ix+1, iy+1) = cmplx(0.0_rp, -(dx*dx + dy*dy) * dk * dk / 2 / ks, rp)
  enddo
enddo

kn%exp_k2 = exp(kn%k2 * dz)

! The source filter's sigmoid, in the same FFT order as K2 and on Genesis's own
! normalized frequency axis: FieldSolverFFT.cpp:133-144, where x runs over
! (ix + shift)/ngrid/xcut.
!
! That axis is the shifted one K2 does not use. So on an even grid the sigmoid's edge
! sits half a wavenumber out from the modes it multiplies, and it is left where Genesis
! puts it. The filter is a transcribed shape, and check_source_filter tests it against
! Genesis. The propagator above carries physics and is corrected instead. Half a
! wavenumber is far inside the sigmoid's own width, which defaults to a twentieth of an
! edge tens of wavenumbers out.
!
! The shifted index spans half of ngrid either side of the axis. So x = 1 is twice the
! Nyquist index, and xcut = 1 puts the sigmoid's edge at the angle lambda/dx, off the
! grid. Genesis leaves the exponent unguarded. That is safe at its own cuts but
! overflows for a cut small enough to put the grid's corner far outside the edge. The
! clamp at 700 keeps the value unchanged, since 1/(1 + exp(700)) is zero in double
! precision.

if (allocated(kn%sigmoid)) deallocate (kn%sigmoid)

if (want_sigmoid) then
  allocate (kn%sigmoid(ngrid, ngrid))
  shift = -0.5_rp * (ngrid - 1)
  do iy = 0, ngrid-1
    yf = (iy + shift) / ngrid / filter%ycut
    do ix = 0, ngrid-1
      xf = (ix + shift) / ngrid / filter%xcut
      iiy = mod(iy + (ngrid+1)/2, ngrid)
      iix = mod(ix + (ngrid+1)/2, ngrid)
      r = (sqrt(xf*xf + yf*yf) - 1) / filter%width
      kn%sigmoid(iix+1, iiy+1) = 1 / (1 + exp(min(r, 700.0_rp)))
    enddo
  enddo
endif

kn%ngrid = ngrid
kn%dgrid = dgrid
kn%ks = ks
kn%dz = dz

end subroutine fel_field_kernel_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_filter_angles (beam, und, lambda, th_mode, th_rho)
!
! Routine to compute the two angles a source filter's edge can be set from
! (fel-physics.md sec-source-filter). Both come from the loaded beam and the undulator
! it enters, so an input file states neither.
!
! th_mode is four diffraction angles of the fundamental mode, 4 lambda / (2 pi sigma),
! with sigma the rms transverse size of the beam. The radiation that can couple to the
! mode lies inside it.
!
! th_rho is sqrt(2 rho lambda / lambda_u). At this angle the resonant wavelength
! red-shifts by rho, and beyond it emission is outside the bandwidth of the
! interaction. rho is the one-dimensional Pierce parameter in Genesis's form, with the
! undulator's coupling and the peak current of the loaded beam.
!
! The ratio of the two goes as sqrt(z_R / L_g). That is 15 on the Benchmark1-SASE line and
! near 1 on a diffraction-dominated machine. So neither angle alone serves as a default.
!
! Input:
!   beam    -- fel_beam_struct: The loaded beam, for the transverse size and the current.
!   und     -- fel_und_struct: The undulator, for aw, ku and the coupling.
!   lambda  -- real(rp): Fundamental radiation wavelength [m].
!
! Output:
!   th_mode -- real(rp): Four mode diffraction angles [rad]. Zero if the beam has no size.
!   th_rho  -- real(rp): The rho angle [rad]. Zero if the beam carries no charge.
!-

subroutine fel_filter_angles (beam, und, lambda, th_mode, th_rho)

type (fel_beam_struct), target :: beam
type (fel_und_struct) und
type (fel_slice_struct), pointer :: sl
real(rp) lambda, th_mode, th_rho
real(rp) q, qx, qy, qxx, qyy, sigma, cur, cur_max, gam, rho3, fc, i_alfven
integer is

!

th_mode = 0
th_rho = 0

! The weighted rms transverse size about the beam's centroid, over both planes, and the
! peak current over slices. A steady-state run is one slice, and its current is that
! slice's.

q = 0;  qx = 0;  qy = 0;  qxx = 0;  qyy = 0;  cur_max = 0
do is = 1, size(beam%slice)
  sl => beam%slice(is)
  if (sl%n == 0) cycle
  q = q + sum(sl%weight(1:sl%n))
  qx = qx + sum(sl%weight(1:sl%n) * sl%x(1:sl%n))
  qy = qy + sum(sl%weight(1:sl%n) * sl%y(1:sl%n))
  qxx = qxx + sum(sl%weight(1:sl%n) * sl%x(1:sl%n)**2)
  qyy = qyy + sum(sl%weight(1:sl%n) * sl%y(1:sl%n)**2)
  cur = sum(sl%weight(1:sl%n)) * c_light / beam%slice_spacing
  cur_max = max(cur_max, cur)
enddo
if (q <= 0) return

sigma = sqrt(max(0.5_rp * ((qxx - qx**2 / q) + (qyy - qy**2 / q)) / q, 0.0_rp))
if (sigma > 0) th_mode = 4 * lambda / (twopi * sigma)

! rho^3 = (I / I_A) fc^2 / (8 gamma^3 sigma^2 ku^2). This is Genesis's form, with the
! undulator's own coupling, which already carries aw. The Alfven current is
! 4 pi eps0 m c^3 / e. With the rest energy in eV it is 4 pi m_electron / (mu_0 c),
! 1.70e4 A.

gam = fel_gamma0(beam)
fc = fel_und_coupling(und, 1)
if (sigma <= 0 .or. cur_max <= 0 .or. und%ku <= 0 .or. fc == 0 .or. gam <= 0) return

i_alfven = 4 * pi * m_electron / (mu_0_vac * c_light)
rho3 = (cur_max / i_alfven) * fc**2 / (8 * gam**3 * sigma**2 * und%ku**2)
if (rho3 <= 0) return

! lambda_u = twopi / ku, so 2 rho lambda / lambda_u is rho lambda ku / pi.

th_rho = sqrt(rho3**(1.0_rp/3.0_rp) * lambda * und%ku / pi)

end subroutine fel_filter_angles

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_filter_margin (tolerance, width) result (margin)
!
! Routine to return the factor that places a source filter's edge outside the angle it
! protects (fel-physics.md sec-source-filter). The edge is the protected angle over this
! margin, so a margin of one leaves the edge on that angle.
!
! The sigmoid of Eq. eq-sigmoid passes amplitude t at the protected angle when the edge
! sits at that angle over 1 - width log[t/(1-t)]. The intensity there is t^2. So a
! stated largest intensity loss fixes t = sqrt(1 - tolerance).
!
! The logarithm is taken as log t + log(1+t) - log(tolerance). That is the same number,
! since tolerance = (1-t)(1+t), and it never forms 1 - t. Formed directly, 1 - t is
! exactly zero for any tolerance below about 2.2e-16, where 1 - tolerance rounds to one.
! Yet at a width of 0.01 such a tolerance still has a margin of 0.60 and gives a valid
! filter. This form also gives a margin of exactly one for a tolerance of 0.75, since
! the three logarithms cancel to the bit. So the edge derived there is the edge derived
! with no tolerance set.
!
! Note: The caller checks that the tolerance lies strictly between zero and one. Both
! ends are singular here. At one the margin is positive infinity, and not negative, so
! a caller testing only the sign of the result would accept it.
!
! Input:
!   tolerance -- real(rp): Largest fraction of source intensity the filter may take
!                 anywhere inside the protected angle. Strictly between 0 and 1.
!   width     -- real(rp): The sigmoid's width as a fraction of its edge.
!
! Output:
!   margin    -- real(rp): The protected angle divided by the edge. At or below zero the
!                 width is too large for the tolerance and no edge meets it.
!-

function fel_filter_margin (tolerance, width) result (margin)

real(rp) tolerance, width, margin, t

!

t = sqrt(1 - tolerance)
margin = 1 - width * (log(t) + log(1 + t) - log(tolerance))

end function fel_filter_margin

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_kernel_index (ngrid, dgrid, ks, dz) result (ik)
!
! Routine to look up the cached kernel entry matching all four keys exactly. A miss
! returns 0. The lookup is read-only, so it can be called inside the parallel slice loops.
!
! Input:
!   ngrid -- integer: Grid points per side.
!   dgrid -- real(rp): Grid spacing [m].
!   ks    -- real(rp): Radiation wavenumber [1/m].
!   dz    -- real(rp): Step length [m].
!
! Output:
!   ik    -- integer: Index into the module kernel cache (0 = not cached).
!-

function fel_kernel_index (ngrid, dgrid, ks, dz) result (ik)

integer ngrid, ik
real(rp) dgrid, ks, dz

!

if (allocated(fel_kernels)) then
  do ik = 1, size(fel_kernels)
    if (fel_kernels(ik)%ngrid == ngrid .and. fel_kernels(ik)%dgrid == dgrid .and. &
        fel_kernels(ik)%ks == ks .and. fel_kernels(ik)%dz == dz) return
  enddo
endif
ik = 0

end function fel_kernel_index

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_diffract (wf, ifld, dz, err_flag)
!
! Routine to diffract field slice ifld in free space over the step dz: FFT, multiply by
! the cached exp(K2 dz), inverse FFT, normalize. There is no source, coupling factor or
! undulator here. This is the field half that the unaveraged mode (fel-physics.md
! sec-unaveraged) shares with the averaged solver. The caller deposits its own source.
!
! Note: fel_field_kernel_init must have built the kernel serially for this grid,
! wavelength and step. A miss is an error, as in fel_field_step.
!
! Input:
!   wf       -- wavefront_struct: The field.
!   ifld     -- integer: Field-record index to diffract.
!   dz       -- real(rp): Step length [m].
!
! Output:
!   wf       -- wavefront_struct: The record advanced by the cached exp(K2 dz) kernel.
!   err_flag -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_field_diffract (wf, ifld, dz, err_flag)

type (wavefront_struct), target :: wf
integer ifld
real(rp) dz
logical err_flag

real(rp) xks, dgrid
integer ngrid_arr(3), ngrid, ik
logical err
character(*), parameter :: r_name = 'fel_field_diffract'

!

err_flag = .true.

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
xks = twopi / wf%wavelength
dgrid = wf%dx

ik = fel_kernel_index(ngrid, dgrid, xks, dz)
if (ik == 0) then
  call out_io (s_error$, r_name, 'FIELD KERNEL NOT INITIALIZED FOR THIS GRID. ' // &
                                 'CALL fel_field_kernel_init FIRST (SERIALLY).')
  return
endif

call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_forward$, err);  if (err) return
wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%exp_k2
call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return
wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp)

if (allocated(wf%Ey)) then      ! Diffraction is polarization-diagonal.
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_forward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * fel_kernels(ik)%exp_k2
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_backward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) / real(ngrid*ngrid, rp)
endif

err_flag = .false.

end subroutine fel_field_diffract

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_diffract_by (wf, ifld, dz, err_flag)
!
! Routine to diffract one field record over a step the kernel cache does not hold. The
! phase table exp(K2 dz) is formed here from the cached wavenumbers of this grid and
! wavelength. The cache holds one step per wavelength, the unaveraged mode's substep.
! The two half steps that put the stored record on its named plane (fel_unaveraged_mod)
! would otherwise evict that entry at every element boundary. dz may be negative, which
! undoes a forward step to rounding.
!
! Input:
!   wf       -- wavefront_struct: The field.
!   ifld     -- integer: Field-record index to diffract.
!   dz       -- real(rp): Step length [m], either sign.
!
! Output:
!   wf       -- wavefront_struct: The record advanced by exp(K2 dz).
!   err_flag -- logical: Set True if no kernel of this grid and wavelength is cached.
!-

subroutine fel_field_diffract_by (wf, ifld, dz, err_flag)

type (wavefront_struct), target :: wf
integer ifld
real(rp) dz
logical err_flag

complex(rp), allocatable :: ph(:,:)
real(rp) xks
integer ngrid_arr(3), ngrid, ik
logical err
character(*), parameter :: r_name = 'fel_field_diffract_by'

!

err_flag = .true.
ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
xks = twopi / wf%wavelength

ik = 0
if (allocated(fel_kernels)) then
  do ik = 1, size(fel_kernels)
    if (fel_kernels(ik)%ngrid == ngrid .and. fel_kernels(ik)%dgrid == wf%dx .and. &
        fel_kernels(ik)%ks == xks .and. allocated(fel_kernels(ik)%k2)) exit
  enddo
  if (ik > size(fel_kernels)) ik = 0
endif
if (ik == 0) then
  call out_io (s_error$, r_name, 'FIELD KERNEL NOT INITIALIZED FOR THIS GRID AND WAVELENGTH. ' // &
                                 'CALL fel_field_kernel_init FIRST (SERIALLY).')
  return
endif

allocate (ph(ngrid, ngrid))
ph = exp(fel_kernels(ik)%k2 * dz)

call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_forward$, err);  if (err) return
wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * ph
call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return
wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp)
if (allocated(wf%Ey)) then
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_forward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * ph
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_backward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) / real(ngrid*ngrid, rp)
endif

err_flag = .false.

end subroutine fel_field_diffract_by

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_diag (wf, ifld, power, on_axis_intensity)
!
! Routine to compute the field diagnostics of field slice ifld from a field in V/m. The
! definitions are Genesis's DiagField::calc in physical units:
!   power             = sum|E|^2 * dgrid^2/(2*Z0) [W]  (cell intensity times cell area)
!   on-axis intensity = |E(center)|^2/(2*Z0) [W/m^2]
! The accumulation order is y outer, x inner.
!
! ifld is the raw record index. For the field at time window position is, pass
! fel_field_index(slip, is, nslice), the unrotation of fel-physics.md sec-slippage.
!
! Input:
!   wf    -- wavefront_struct: The field.
!   ifld  -- integer: Field-record index to evaluate.
!
! Output:
!   power             -- real(rp): Radiation power of the record [W].
!   on_axis_intensity -- real(rp): Genesis's on-axis intensity diagnostic.
!-

subroutine fel_field_diag (wf, ifld, power, on_axis_intensity)

type (wavefront_struct) wf
integer ifld
real(rp) power, on_axis_intensity
real(rp) scl, wei
integer ix, iy, ngrid_arr(3), ngrid, ic

!

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
scl = wf%dx**2 / (2 * (mu_0_vac * c_light))

power = 0
do iy = 1, ngrid
  do ix = 1, ngrid
    wei = real(wf%Ex(ix,iy,ifld), rp)**2 + aimag(wf%Ex(ix,iy,ifld))**2
    power = power + wei
  enddo
enddo

! With a live second polarization, power and intensity are totals (|Ex|^2 + |Ey|^2).

if (allocated(wf%Ey)) then
  do iy = 1, ngrid
    do ix = 1, ngrid
      wei = real(wf%Ey(ix,iy,ifld), rp)**2 + aimag(wf%Ey(ix,iy,ifld))**2
      power = power + wei
    enddo
  enddo
endif
power = power * scl

ic = ngrid/2 + 1
on_axis_intensity = (real(wf%Ex(ic,ic,ifld), rp)**2 + aimag(wf%Ex(ic,ic,ifld))**2) &
                    / (2 * (mu_0_vac * c_light))
if (allocated(wf%Ey)) then
  on_axis_intensity = on_axis_intensity + &
                    (real(wf%Ey(ic,ic,ifld), rp)**2 + aimag(wf%Ey(ic,ic,ifld))**2) &
                    / (2 * (mu_0_vac * c_light))
endif

end subroutine fel_field_diag

end module fel_track_mod
