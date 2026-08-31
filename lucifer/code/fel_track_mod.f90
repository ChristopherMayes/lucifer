!+
! Module fel_track_mod
!
! The FEL step: period-averaged particle push, source deposition and FFT field solve
! inside an undulator segment, operating on the Bmad-coordinate packed beam of
! fel_beam_mod. Runs steady state (one slice) or time dependent (many slices with
! slippage).
! The physics is transcribed from Genesis 1.3 Version 4, which GPL permits. Physics,
! Genesis provenance (routine by routine), and validation: lucifer/doc/fel-physics.md
! (sec:core, sec:transverse, sec:field, sec:slippage).
!
! The step composition (transverse half step, longitudinal advance, transverse half
! step, field step, slippage) is Beam::track plus Gencore's steps 4 and 5.
!
! Time dependence. Beam slice is couples to field slice fel_field_index(slip, is, nslice)
! = 1 + mod(is-1 + first, nslice), where first is Genesis's Field::first: the field
! record is a circular buffer over the wavefront's third index, rotated by slippage
! rather than moved. Slippage accumulates in units of the radiation wavelength and
! rotates the record one slice whenever it exceeds 0.8 of the slice spacing. The field
! slice that was coupled to the head of the time window is discarded and re-enters at
! the tail zeroed. Radiation slips out of the window at the head and fresh vacuum enters
! behind the bunch. The caller owns the schedule (what slips when, including the drift
! autophasing of Lattice::calcSlippage). This module owns the mechanics. Anything
! reading the field record in time order (dumps, per-slice field diagnostics) must
! unrotate through fel_field_index, exactly as Genesis's writers do (sec:slippage).
!
! Coordinates. The stored state is Bmad's (x, px/p0, y, py/p0, z, pz) plus weight (see
! fel_beam_mod). The longitudinal RK4 runs in Genesis's chart (theta, gamma) as
! per-particle working variables, derived at step entry and converted back at step exit:
!
!   entry:  gamma = sqrt((p0_mc*(1+pz))^2 + 1),  theta = phi0 + ks*z/beta
!   exit:   pz = (sqrt(gamma^2-1) - p0_mc)/p0_mc,  z = -beta*(phi0_new - theta)/ks
!
! Two facts make this clean rather than a compromise. First, theta <-> (phi0, z) is an
! affine change of variables with a state-independent shift, and RK4 is exactly invariant
! under affine maps. Integrating in theta is therefore the same discretization as
! integrating in z, and the audited, verbatim RK4 transcription survives unchanged.
! Second, the per-step gamma <-> pz round trip costs ~1 ulp of gamma per step (the
! subtraction in pz is exact once P is formed). That is the same order as the step's own
! rounding. The common reference phase phi0 advances once per step per beam
! (fel_phi0_rate in fel_beam_mod). The split keeps z small and bunch-scale.
!
! Weights: the source deposition scales per particle as c*w_j/slice_spacing where Genesis
! has current/N. The algebra is identical for uniform weights, and correct for nonuniform
! ones (each macroparticle radiates its own charge).
!
! Units. The field is in V/m throughout, the wavefront_struct convention. There is no
! internal unit system. The formulas here follow from Genesis's internal-unit forms via
! the exact relation u = E*ks/(sqrt(2)*m_electron) (composing the dump scale
! dfl = u*dgrid*eev/(ks*sqrt(Z0)) with E = dfl*sqrt(2*Z0)/dgrid). Under this relation the
! coupling, source and power expressions come out simpler than the originals: the energy
! exchange coupling is fc*conj(E)/(sqrt(2)*m_electron), the source scale is
! fc*Z0*sqrt(2)*c*delz*w_j/(4*dgrid^2*slice_spacing), and power is sum|E|^2*dgrid^2/(2*Z0).
! Constants are Bmad's. m_electron equals Genesis's eev exactly. Z0 = mu_0_vac*c_light
! differs from Genesis's truncated 376.73 by 8e-7 relative, the accepted floor of the
! Genesis comparison.
!
! Deliberately absent, per the deliverable: space charge, wakes, incoherent synchrotron
! radiation, harmonics beyond the coupling formula, the source filter, orbit or field
! errors, one4one sorting and slice migration, and the periodic time-window boundary
! (Genesis's periodic option: the non-periodic zero fill is what is transcribed).
!-

module fel_track_mod

use fel_beam_mod
use fel_collective_mod
use wavefront_mod

use, intrinsic :: iso_c_binding, only: c_double

implicit none

!+
! Struct fel_und_struct
!
! One undulator segment, constant parameters along it: one contiguous run of aw > 0 steps
! in Genesis's unrolled lattice, divided into nstep equal steps of dz = l/nstep with
! nstep = round(l/delz_target). kx, ky carry Genesis's unroll scaling by ku^2
! (fel-physics.md sec-element).
!-

type fel_und_struct
  real(rp) :: aw = 0          ! rms undulator parameter.
  real(rp) :: ku = 0          ! Undulator wavenumber twopi/lambdau [1/m].
  real(rp) :: kx = 0          ! Natural focusing, deck kx * ku^2 [1/m^2].
  real(rp) :: ky = 0          ! Natural focusing, deck ky * ku^2 [1/m^2].
  real(rp) :: ax = 0, ay = 0  ! Transverse offset of the undulator field [m].
  real(rp) :: tilt = 0        ! Wiggle-plane tilt [rad]: planar wiggles along
  real(rp) :: cos_t = 1, sin_t = 0   !   e = (cos_t, sin_t). Cached for the maps.
  complex(rp) :: pol(2) = [(1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp)]
                              ! Polarization 2-vector on the (Ex, Ey) pair (manual
                              ! sec:field vector convention): planar (cos t, sin t),
                              ! helical (1, -i)/sqrt(2). The averaged kick reads
                              ! E_eff = conj(pol).E and the deposit writes pol*src.
                              ! This reproduces the scalar path exactly when only
                              ! one polarization is live.
  logical :: helical = .false.
  logical :: bmad_transport = .false.  ! Transverse maps: transcribed TrackBeam (default)
                                       !   or the flattened Bmad periodic-wiggler kernel.
  integer :: nstep = 0        ! Number of integration steps over the segment.
  real(rp) :: dz = 0          ! Step length [m].
  integer :: source_model = 0 ! fel_source_deposit$ (default) or fel_source_coherent$.
end type

!+
! Struct fel_coherent_struct
!
! One slice's coherent-source summary (fel-physics.md sec-coherent-source), computed from the
! just-advanced particles each step: the source phasor S (exactly the sum the deposit
! would scatter: SUM crsource = S is the normalization contract), the charge-weighted
! centroid and 2x2 second moments (widths + tilt: the centering and tilt are this
! port's extension of the paper, which assumes an on-axis untilted beam), the
! Laguerre-Gauss order-0/1 sums feeding the global kappa fit, and the Gaussianity
! guard metric (excess kurtosis, worst plane).
!-

type fel_coherent_struct
  complex(rp) :: S = 0                   ! Source phasor: SUM part_j (sin+icos theta_j).
  real(rp) :: x0 = 0, y0 = 0             ! Charge-weighted centroid [m].
  real(rp) :: sxx = 0, sxy = 0, syy = 0  ! Charge-weighted central second moments [m^2].
  real(rp) :: g2sum = 0                  ! |g0|^2 + |g1|^2 (Eq 20 sums, charge-weighted).
  real(rp) :: b2 = 0                     ! |B|^2, same weighting (kappa's Eq 26 ratio).
  real(rp) :: kurt = 0                   ! max plane |excess kurtosis| (guard metric).
  real(rp) :: m_ind = 1                  ! Independent transverse samples: N_eff/nbins
                                         !   (beamlet copies share coordinates). The
                                         !   guard tests SIGNIFICANCE, not raw kurtosis
                                         !   (sample kurtosis spreads as sqrt(24/m_ind)).
  real(rp) :: wsum = 0                   ! Slice charge [C]: the guard tests on it
                                         !   (edge slices resample from few source
                                         !   points and look degenerate while their
                                         !   source is negligible).
  logical :: ok = .false.                ! Enough charge and spread to fit.
end type

!+
! Struct fel_slip_struct
!
! The slippage state of one field record: Genesis's Field::first and Field::accuslip,
! plus the two run facts they are meaningless without. One per wavefront. The default is
! the steady state: timerun false makes fel_apply_slippage a no-op (as Genesis's)
! and first = 0 makes fel_field_index the identity.
!-

type fel_slip_struct
  logical :: timerun = .false.  ! Time-dependent run? False: slippage is a no-op.
  integer :: first = 0          ! Rotation offset of the field record, 0-based (Field::first).
  real(rp) :: accuslip = 0      ! Accumulated slippage [radiation wavelengths] (Field::accuslip).
  real(rp) :: sample = 1        ! Slice spacing / radiation wavelength (Control::sample).
  real(rp) :: u_escaped = 0     ! Energy [J] transmitted out of the window by slippage
                                ! (summed over every zero-filled slice). The TD energy
                                ! ledger's escape column, so E_beam + U_window + U_escaped
                                ! closes in a wake-free run.
end type

!+
! Struct fel_bank_struct
!
! Scratch carrier for the field slices one fel_apply_slippage call transmits out of the
! window: the caller passes it when it wants the light itself, not just its banked
! energy (keep_escaped_field). Reset and refilled per call. The caller drains it
! immediately (streams to file), so peak memory is a handful of grid planes -- one per
! rotation of the call, ~1 inside undulators, ~10 over an interlude.
!-

type fel_bank_struct
  complex(wf_rp), allocatable :: plane(:,:,:)  ! Transmitted planes, in transmission order.
  complex(wf_rp), allocatable :: plane_y(:,:,:)  ! Ey planes, allocated when Ey is live.
  integer :: n = 0                             ! How many this call transmitted.
end type

!+
! Struct fel_field_struct
!
! One radiation field of the walk's field set: the harmonic number, the wavefront
! record, and that field's own slippage state and escape bank. The walk carries an
! ordered set of these (Genesis's vector<Field*>), with the FUNDAMENTAL always
! entry 1. The ponderomotive phase, the phi0 advance and the slippage schedule are
! all defined against the fundamental. A harmonic field couples through fc(h) and the
! phase h*theta and diffracts at its own wavelength. A single-entry set is the
! pre-harmonic walk, bit for bit. wf%wavelength = fundamental wavelength / harm. Every
! field shares the time window, so the slippage state advances in fundamental-wavelength
! units for all of them (Genesis's one Control::sample) and the records rotate in
! lockstep.
!-

type fel_field_struct
  integer :: harm = 1               ! Harmonic number h (Field::harm).
  type (wavefront_struct) :: wf
  type (fel_slip_struct) :: slip
  type (fel_bank_struct) :: bank
end type

! Cached kernels for fel_field_step, mirroring FieldSolverFFT::init: one entry per
! WAVELENGTH (per harmonic), each rebuilt when its grid or step changes. All module
! state, built SERIALLY by fel_field_kernel_init before any parallel slice loop and
! read-only ever after. fel_field_step never rebuilds: it looks its entry up
! (fel_kernel_index) and errors on a miss, because a rebuild from inside a parallel
! loop would race. The per-call source accumulator is a local of fel_field_step, one
! per invocation, so concurrent slices deposit into their own.

type fel_kernel_struct
  complex(rp), allocatable :: k2(:,:)      ! -i (kx^2+ky^2)/(2 ks), FFT order.
  complex(rp), allocatable :: exp_k2(:,:)  ! exp(K2 * dz), the step propagator.
  integer :: ngrid = 0
  real(rp) :: dgrid = 0, ks = 0, dz = 0
end type

! Named values of the fel_tracking lattice attribute (fel-physics.md sec-element), in Bmad's
! named-integer convention. Lattices use the same names via one-line variable
! definitions (e.g. "fel_unaveraged = 1" before the element that sets it).

! The source model (fel-physics.md sec-coherent-source): the standard per-particle deposit
! (the referee) or the SIMPLEX-hybrid coherent-Gaussian source (Tanaka, PRAB 27,
! 030703 (2024); arXiv:2310.20197). In the coherent source the spatially incoherent
! artifact is dropped. The slice bunch factor keeps the physical shot noise. The
! transverse shape is a guarded Gaussian from phasor sums and charge moments.

integer, parameter :: fel_source_deposit$ = 0
integer, parameter :: fel_source_coherent$ = 1

integer, parameter :: fel_transcribed$ = -1   ! Transcribed-Genesis transverse maps
                                              !   (validation-internal, Genesis tiers).
integer, parameter :: fel_averaged$ = 0       ! Averaged, bmad_standard kernel maps
                                              !   (the unset default).
integer, parameter :: fel_unaveraged$ = 1     ! The unaveraged mode.

type (fel_kernel_struct), allocatable, target, private, save :: fel_kernels(:)

! One libm call for the (sin, cos) pair (lucifer/code/fel_sincos.c). NOT bit-identical
! to gfortran's own sin/cos intrinsics: they differ from libm by one ulp on ~2e-6 of
! arguments (73 mismatches in a 44M-point sweep of the theta domain). Adopting
! this was therefore a named value change -- every benchmark tier re-measured and
! re-recorded (doc/validation.md, "The particle-path cost").

!+
! Subroutine fel_sincos (theta, s, c)
!
! Routine to compute the (sin, cos) pair with one libm sincos call (C shim,
! fel_sincos.c). See the named-value-change note above.
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
! track1_custom hook: outside a driver's own FEL walk, an FEL element (a wiggler with
! tracking_method = custom) is just a periodic wiggler. Delegate to Bmad's standard
! kernel. The reference time/energy pass inside bmad_parser and any seam-side track1
! resolve through this, so the element carries the resonant undulation delay
! from Bmad's own code. Kept at module scope deliberately: gfortran implements pointers
! to internal procedures with stack trampolines. Apple Silicon's non-executable
! stack turns those into a segfault at the first call.
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

call fel_assert_wiggler_sane (ele)

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
! The brief's 7.5 assertions, enforced at the first touch of the element: the
! reference time/energy pass INSIDE bmad_parser, through the hooks above. Enforcing
! them any later is too late. A missing b_max parses cleanly and only fails downstream
! with an unrelated message. A fieldmap field_calc segfaults track_a_wiggler during
! the parse itself. Refusal is by name so a lattice author knows which attribute to fix.
! (The reference pass runs before lat_sanity_check, so these fire first. Bmad's own
! sanity check would also refuse a missing l_period by name if this were removed.)
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

! These refusals fire at parse time from the tracking hooks and must end the program
! with a NONZERO exit (the harness's refusal checks assert it). err_exit cannot be
! used here: its deliberate integer-divide traceback bomb does not trap on arm64 and
! its final bare stop exits 0. The stops below stay.

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

! The integration step: with neither ds_step nor num_steps set, Bmad's bookkeeper
! falls back to bmad_com%default_ds_step (0.2 m -- thirteen periods of this device,
! garbage FEL physics that would run without complaint). Bmad itself only warns.
! An FEL element must say its step.

if (ele%value(ds_step$) == bmad_com%default_ds_step) then
  call out_io (s_fatal$, r_name, 'FEL ELEMENT HAS NO INTEGRATION STEP; SET DS_STEP (GENESIS''S DELZ)', &
                                 'OR NUM_STEPS ON: ' // trim(ele%name), &
                                 '(ITS DS_STEP EQUALS BMAD_COM%DEFAULT_DS_STEP, THE UNSET FALLBACK. IF YOU', &
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
! make_mat6_custom hook, the transfer-matrix companion of fel_ele_as_wiggler: an FEL
! element's mat6_calc_method resolves to custom (auto follows tracking_method = custom),
! and Bmad's bookkeeping calls through make_mat6_custom_ptr unconditionally. A
! program that leaves the pointer null segfaults at a jump to address zero. Delegate to the
! standard periodic-wiggler kernel with matrix propagation, filling ele%mat6, ele%vec0
! and end_orb per the make_mat6_bmad convention.
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

call fel_assert_wiggler_sane (ele)

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
! ifld = 1 + mod(is-1 + first, nslice), the Fortran form of Genesis's
! (is + field->first) % field.size() (fel-physics.md sec-slippage). The same mapping
! restores time order when reading the record out: field slice
! fel_field_index(slip, is, nslice) is the field at time window position is.
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
! Subroutine fel_apply_slippage (slip, wf, slippage, bank, harm)
!
! Routine to account the slippage of one step and rotate the field record when it
! crosses a slice boundary. Transcribed from Control::applySlippage (fel-physics.md
! sec:slippage) reduced to one shared-memory node. The MPI ring exchange between
! adjacent ranks is the identity when the ring has one member. This node is
! simultaneously rank 0 and rank size-1, so the non-periodic zero fill always applies
! to the transmitted slice. What remains is Genesis's exact bookkeeping: accumulate in
! units of the radiation wavelength, and while |accuslip| exceeds 0.8*sample, rotate
! the record one slice -- the slice coupled to the head of the time window is zeroed
! and becomes the new tail (radiation leaves the window at the head and fresh vacuum
! enters behind the bunch). Backward slippage is transcribed with the same fidelity,
! direction reversed.
!
! The caller owns the schedule (per-step undulator slippage and the drift autophasing
! of Lattice::calcSlippage) and calls this after the field step, Gencore's step 5.
! A steady-state record (timerun false) is left untouched, exactly as Genesis's.
!
! Input:
!   slip        -- fel_slip_struct: Slippage state of this field record.
!   wf          -- wavefront_struct: The field record.
!   slippage    -- real(rp): Slippage of this step [radiation wavelengths].
!
! Output:
!   slip, wf    -- Updated state and rotated record.
!-

subroutine fel_apply_slippage (slip, wf, slippage, bank, harm)

type (fel_slip_struct) slip
type (wavefront_struct) wf
type (fel_bank_struct), optional :: bank
integer, optional :: harm
real(rp) slippage, t_slice
complex(wf_rp), allocatable :: grow(:,:,:)
integer nslice, last, direction

!

! The transmitted slice's light time, for the escape bank: slice_spacing / c. For a
! harmonic field wf%wavelength is the fundamental's / harm while sample stays in
! fundamental units (every field shares the window and rotates in lockstep). The
! product must therefore be scaled back up by harm. harm absent or 1 leaves the
! fundamental's expression untouched.

t_slice = (slip%sample * wf%wavelength / c_light)
if (present(harm)) then
  if (harm /= 1) t_slice = t_slice * harm
endif

if (present(bank)) bank%n = 0
if (.not. slip%timerun) return

slip%accuslip = slip%accuslip + slippage
nslice = size(wf%Ex, 3)

do while (abs(slip%accuslip) > slip%sample * 0.8_rp)
  direction = 1
  if (slip%accuslip < 0) direction = -1
  slip%accuslip = slip%accuslip - slip%sample * direction

  ! The transmitted slice: the last of the record in time order for forward slippage,
  ! the first for backward, 0-based here as in Genesis (sec:slippage).

  last = mod(slip%first + nslice - 1, nslice)
  if (direction < 0) last = mod(last + 1, nslice)

  ! One node: send/receive to self, then the non-periodic zero fill of the transmitted
  ! slice. Its energy leaves the simulation here. Bank it first, so the time-dependent
  ! energy ledger can close: the slice's power (same convention as fel_field_diag) times
  ! its light-time slice_spacing/c. Wakes would be a second, unbanked exit channel from
  ! the beam. The unaveraged mode, whose ledger this feeds, refuses them by name.

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

end subroutine fel_apply_slippage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_und_coupling (und, h) result (fc)
!
! The coupling factor fc(h), transcribed from Undulator::fc. Helical couples the
! fundamental with strength aw and nothing else. Planar couples odd harmonics through
! JJ = J_{(h-1)/2}(xi) - J_{(h+1)/2}(xi), xi = h/2 * aw^2/(1+aw^2), sign (-1)^((h-1)/2),
! and even harmonics not at all (a polarisation limitation of the planar coupling).
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
! Transverse dependence of the undulator field, first order form used in the particle
! gather. Undulator::faw.
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
! Square of the transverse dependence, used in the source deposition. Undulator::faw2.
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
! Subroutine fel_track_und_step (und, beam, ff, coll, err_flag)
!
! Routine to advance the whole beam and the field set one integration step of length
! und%dz inside an undulator segment, in Genesis's order: transverse half step,
! longitudinal RK4 (gathering every field of this moment and holding the couplings
! through the stages), transverse half step, then source deposition and field solve per
! field, harmonic loop outermost. The common phase phi0 advances once per step, with
! the undulator's ku as the reference wavenumber. Each beam slice couples to its field
! slice through fel_field_index per field. Slippage is NOT applied here -- the caller
! owns the schedule and calls fel_apply_slippage after this step (Gencore's step 5).
!
! Input:
!   und         -- fel_und_struct: Segment parameters.
!   beam        -- fel_beam_struct: The beam. Steady state: one slice.
!   ff          -- fel_field_struct(:): The field set [V/m]. Entry 1 is the fundamental,
!                    one field slice per beam slice each.
!   coll        -- fel_collective_struct: Collective-effects state.
!
! Output:
!   beam, ff    -- Advanced one step.
!   err_flag    -- logical: Set True on error.
!-

subroutine fel_track_und_step (und, beam, ff, coll, err_flag)

type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
type (fel_collective_struct) coll
logical err_flag

real(rp) ks, phi0_new, kappa, b2_tot, g2_tot
type (fel_coherent_struct), allocatable :: coh(:)
integer is, io, nslice, nslice_f, ngrid_arr(3)
logical any_err, err
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

! ks, phi0 and the slippage schedule are all FUNDAMENTAL quantities (field 1). A
! harmonic field enters only through its coupling fc(h), its phase h*theta and its own
! diffraction kernel.

ks = twopi / ff(1)%wf%wavelength
phi0_new = beam%phi0 + und%dz * fel_phi0_rate(ks, und%ku, fel_p0_mc(beam))

! Everything cross-slice happens serially here, before and between the parallel loops:
! the phi0 advance above, the kernel builds and the long-range space-charge profile
! below, and slippage in the caller. Inside the loops each slice touches only its own
! particle arrays and its own field slice (the beam-to-field mapping is a bijection).
! Each slice's arithmetic is independent of which thread runs it, so results are
! bit-identical across thread counts -- the check the benchmark harness holds.

do io = 1, size(ff)
  ngrid_arr = wavefront_shape(ff(io)%wf)
  call fel_field_kernel_init (ngrid_arr(1), ff(io)%wf%dx, twopi / ff(io)%wf%wavelength, und%dz)
enddo

if (.not. allocated(coll%long_esc)) allocate (coll%long_esc(nslice))
call fel_longrange_esc (coll%efield, beam, fel_gamma0(beam), und%aw, coll%long_esc)

! Per-slice sequence in Genesis's order (Beam::track): transverse half, longitudinal
! advance (ez inside the RK), the wake's gamma decrement, transverse half. All four are
! per-slice pure, so folding them into one loop is arithmetic-identical to Genesis's
! four sweeps.

!$OMP parallel do
do is = 1, size(beam%slice)
  if (und%bmad_transport) then
    call fel_transverse_track_bmad (und, beam, beam%slice(is), und%dz/2, .true.)
  else
    call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
  endif
  call fel_advance (und, beam, beam%slice(is), ff, und%dz, phi0_new, coll, is)
  call fel_wake_apply_slice (coll%wake, beam, is, und%dz)
  if (und%bmad_transport) then
    call fel_transverse_track_bmad (und, beam, beam%slice(is), und%dz/2, .false.)
  else
    call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
  endif
enddo
!$OMP end parallel do

beam%phi0 = phi0_new

! The coherent source (fel-physics.md sec-coherent-source): per-slice phasor, moments and LG
! sums from the just-advanced particles (slice-parallel, each slice its own data),
! then the ONE global kappa (Tanaka Eq 26: the integrals run over the whole window,
! so this is the cross-slice serial point), then the guard. A slice whose transverse
! charge profile is measurably non-Gaussian is refused by name with its number: the
! Gaussian model would bias the gain there. Thresholds are set by the distorted-beam
! checks.

if (und%source_model == fel_source_coherent$) then
  allocate (coh(size(beam%slice)))
  !$OMP parallel do
  do is = 1, size(beam%slice)
    call fel_coherent_prep (und, beam, beam%slice(is), ks, coh(is))
  enddo
  !$OMP end parallel do

  ! The kappa fit's window sums (Eq 26), CHARGE-WEIGHTED: g_q and B are built from
  ! normalized weights. An unweighted sum would let near-empty edge slices (their
  ! normalized phasors are O(1) noise) poison the global width fit for any
  ! real bunch profile. wsum^2 keeps both sums quadratic in the same measure.

  b2_tot = 0;  g2_tot = 0
  do is = 1, size(beam%slice)
    b2_tot = b2_tot + coh(is)%b2 * coh(is)%wsum**2
    g2_tot = g2_tot + coh(is)%g2sum * coh(is)%wsum**2
  enddo
  do is = 1, size(beam%slice)

    ! Charge-weighted: an edge slice (Gaussian tail, resampled from few source points)
    ! looks degenerate while contributing negligible source. Only slices carrying
    ! real charge are held to the Gaussian-profile requirement.

    if (coh(is)%wsum < 0.05_rp * maxval(coh%wsum)) cycle
    if (coh(is)%ok .and. coh(is)%kurt > max(0.5_rp, 5 * sqrt(24.0_rp / coh(is)%m_ind))) then

      ! Significance, not raw kurtosis: the sample statistic spreads as sqrt(24/m_ind)
      ! (m_ind independent transverse samples). A small-M Gaussian beam thus does not
      ! trip the guard on noise, while a genuinely structured profile (double horn
      ! ~ -1.5) does at any M where it is resolvable.

      call out_io (s_error$, r_name, 'COHERENT SOURCE: slice \i0\ transverse charge profile is NOT', &
        'GAUSSIAN ENOUGH (excess kurtosis \es10.2\ , threshold \es10.2\ at this sample size):', &
        'the coherent-Gaussian model would bias the gain here. Use source_model = "deposit".', &
        i_array = [is], r_array = [coh(is)%kurt, max(0.5_rp, 5 * sqrt(24.0_rp / coh(is)%m_ind))])
      return
    endif
  enddo
  kappa = 1
  if (g2_tot > 0) kappa = sqrt(b2_tot / (4 * pi * g2_tot))
endif

! Field solve, harmonic loop outermost (each pass is self-contained, fel-physics.md sec-field-set). The
! deposit reads the just-advanced particles, so every field sees the same beam state.

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

err_flag = .false.

end subroutine fel_track_und_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_interlude_genesis (qf, length, beam, ff, coll, err_flag)
!
! Routine to advance one field-free interlude element (a drift or a quadrupole) the
! way Genesis does it, as one integration step: transverse half step, the longitudinal
! advance with the path-length term sampled at mid element and Genesis's drift
! reference xku = ks*0.5/gamma0/gamma0 (division order kept for bit identity,
! fel-physics.md sec-interlude), the wake's gamma decrement, transverse half step,
! field diffraction with zero source.
! Slippage is NOT applied here. The caller schedules it after the step, as with
! fel_track_und_step.
!
! The longitudinal advance has two paths. With space charge off, the RK4 collapses
! exactly (slope theta-independent, gamma constant) and the collapsed step is kept for
! bit-identity with the pre-collective code. With space charge on, gamma changes inside
! the step, so the full RK4 runs with rpart = 0 and the per-particle ez. That is
! Genesis's actual code path in drifts (BeamSolver::advance with aw = 0).
!
! This is the Genesis interlude model, transcribed. Running the full lattice with
! it isolates what the Bmad seam changes: the seam integrates the path-length term
! exactly through the quad map where this samples it once. See doc/validation.md.
! The production configuration is the seam.
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

sc_active = coll%efield%on .and. (coll%efield%nz >= 1 .or. coll%efield%longrange)
if (.not. allocated(coll%long_esc)) allocate (coll%long_esc(size(beam%slice)))
call fel_longrange_esc (coll%efield, beam, gamma0, 0.0_rp, coll%long_esc)

! The per-particle temporaries live at routine scope, so the parallel loop must make
! them private explicitly. A missed one here is a race, which is what the harness's
! thread-count-independence check exists to catch.

!$OMP parallel do private(sl, ip, gam, beta, theta, px_g, py_g, btpar, btpar0, slope, p_mc)
do is = 1, size(beam%slice)
  sl => beam%slice(is)

  call interlude_transverse_half (sl, q_hat, length/2)

  if (sc_active) then

    call interlude_advance_full_rk (sl, is)

  else

    ! theta advance, the collapsed RK4: with no field and no space charge the slope is
    ! theta independent and constant through the stages, so RK4 reduces to one exact
    ! step. btpar = 1 + px_g^2 + py_g^2 (aw = 0), px_g = gamma*beta_x = px * p0_mc.

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

      ! Back to z: tau = (phi0_new - theta)/ks. pz is unchanged (no energy change), so
      ! beta is unchanged too.
      sl%z(ip) = -beta * (phi0_new - theta) / xks
    enddo
  endif

  call fel_wake_apply_slice (coll%wake, beam, is, length)

  call interlude_transverse_half (sl, q_hat, length/2)
enddo
!$OMP end parallel do

beam%phi0 = phi0_new

! Field diffraction: Genesis's kernel, zero source (aw = 0 makes the coupling exactly
! zero, matching Genesis skipping the deposition when not in an undulator). The beam
! slice to field slice mapping is a bijection, so looping beam slices covers every field
! slice exactly once. Kernels built serially first: a lattice can open with an interlude,
! so this routine cannot rely on an undulator step having initialized them. Harmonic
! loop outermost, each field diffracting at its own wavelength.

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
! Routine to apply half an interlude's transverse map to one slice: the focusing
! map per plane at the interlude's quad strength (Genesis's drift/quad composite).
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
! Subroutine interlude_advance_full_rk (sl, is)
!
! Routine to advance one slice's longitudinal coordinates through the interlude:
! the field-free RK4 step (rpart = 0) plus the collective Ez, Genesis's method.
!-

subroutine interlude_advance_full_rk (sl, is)

! Genesis's actual drift path (BeamSolver::advance with aw = 0): the full RK4 with
! rpart = 0 and the per-particle space-charge ez held through the stages. gamma changes
! here, so the chart round trip at exit is fel_advance's. Locals only, safe inside the
! parallel slice loop.

type (fel_slice_struct) sl
integer is

real(rp), allocatable :: ez(:)
real(rp) gam, beta, theta, px_g, py_g, btpar, esc_loss, p_mc
integer ip

allocate (ez(max(1, sl%n)))
call fel_shortrange_ez (coll%efield, beam, sl, gamma0**2, xks, ez)   ! gz2 at aw = 0.
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
! Routine to advance the transverse coordinates over delz inside an undulator: drift plus
! natural focusing, transcribed from TrackBeam::track (quad, corrector and chicane
! branches dropped: they never overlap undulator fields here).
!
! In the stored normalization the per particle longitudinal factor is
!
!   gz_hat = gammaz/p0 = sqrt((1+pz)^2 - px^2 - py^2 - (aw/p0_mc)^2)
!
! (Genesis's gammaz = sqrt(gamma^2 - 1 - aw^2 - px_g^2 - py_g^2), divided by p0), and the
! maps read identically with px in Bmad units:
!
!   drift:  x += px*delz/gz_hat
!   quad:   foc^2 = q_hat/gz_hat with q_hat = q/p0_mc;  x' = a1 x + a2 px/gz_hat;
!           px' = a3 x gz_hat + a1 px
!
! The effective strengths (fel-physics.md sec-natfocus): qnat_{x,y} = k_{x,y}*aw^2/
! (gamma0*betpar0), betpar0 = sqrt(1 - (1+aw^2)/gamma0^2), with Genesis's reference gamma.
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
! The priced transport alternative (fel-physics.md sec-element): the transverse maps of Bmad's own
! periodic-wiggler kernel (track_a_wiggler.f90:90-186), flattened. Quadrupole bodies go
! via quad_mat2_calc with the per-particle 1/rel_p^2 chromatic scaling and the
! half-octupole edge kicks, using the TRACKING-LOCAL k1 values (7.5: never the stored
! k1x/k1y attributes, whose helical sign disagrees). The z bookkeeping of the kernel
! (path-length dz terms, low-energy correction, the end-of-element undulation factor)
! is deliberately ABSENT: the ponderomotive phase evolution, including the aw^2 and
! px^2+py^2 path terms, lives in fel_advance's RK. Applying it here too would
! double-count.
!
! Structure: every leading half-step applies half the octupole kick then the quad body.
! Every trailing half-step applies the body then half the kick. Adjacent half-kicks
! between steps sum to Bmad's full inter-step kick. The segment faces get the half kick,
! exactly Bmad's n_step loop. Model differences from the transcribed TrackBeam (measured,
! not argued): Bmad's 1/rel_p^2 chromaticity against Genesis's 1/gz_hat, and the octupole
! edge kicks Genesis does not have.
!
! g_max reconstructs from aw: c*b_max = K*ku*m_e with K = aw (helical) or aw*sqrt(2)
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
  ! untilted map, rotate out -- exactly Bmad's tilt_coords composition. sin_t = 0
  ! (every untilted element) skips both rotations, arithmetic untouched.

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
! cosine/sine are (c, s).
!-

subroutine rot_xy (x, px, y, py, c, s)

! Rotate the transverse pair into (s > 0 passed) or out of (s < 0) the wiggle frame.

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
! One transverse plane's map over delz: TrackBeam::applyDrift for q_hat = 0, applyFQuad
! for q_hat > 0, applyDQuad for q_hat < 0, in the normalized variables of
! fel_transverse_track (px = P_x/p0, gz_hat = gammaz/p0, q_hat = q/p0).
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
! Subroutine fel_advance (und, beam, sl, wf, ifld, delz, phi0_new, coll, is)
!
! Routine to advance the longitudinal plane of every particle over delz. Transcribed from
! BeamSolver::advance. Gather the field at (x, y) by bilinear interpolation from field
! slice ifld (the rotated-record index, fel_field_index). Form
! rpart = (fc(h)/ks)*faw*conj(E_h) per field of the set ff (the fundamental is entry 1,
! and a single-entry set takes the pre-harmonic scalar path verbatim). Integrate
! (theta, gamma) by the verbatim RK4 with the rpart set, px, py, faw AND the
! space-charge ez held fixed through the stages. The harmonic phases h*theta are live
! per stage (fel_ode_multi). ez per particle is
! fel_shortrange_ez(ip) - long_esc(is)/m_electron, exactly Genesis's
! ez = getEField(ip) + eloss (fel-physics.md sec-eom, sec:spacecharge). is is this
! slice's beam index.
!
! (theta, gamma) are derived at entry from the stored (z, pz) and written back at exit
! using phi0_new, the common phase at the end of this step. See the module header for why
! this chart change is exact for RK4 and ~1 ulp for the energy.
!
! Input:
!   und      -- fel_und_struct: Undulator parameters (incl. source model).
!   beam     -- fel_beam_struct: The beam (reference energy, wavelength).
!   sl       -- fel_slice_struct: One slice's packed particles.
!   ff(:)    -- fel_field_struct: The field set (gathered per particle each stage).
!   delz     -- real(rp): Step length [m].
!   phi0_new -- real(rp): The field's phi0 after this step.
!   coll     -- fel_collective_struct: Collective terms (wake/space-charge Ez).
!   is       -- integer: Slice index (for the collective lookups).
!
! Output:
!   sl       -- fel_slice_struct: gamma/theta (and chart z) advanced by delz.
!-

subroutine fel_advance (und, beam, sl, ff, delz, phi0_new, coll, is)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_collective_struct) coll
integer is
real(rp) delz, phi0_new

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
! (BeamSolver::advance's order). The short-range solve is per-call-local, so this is
! parallel-slice safe. long_esc was refreshed serially by the caller.

allocate (ez(max(1, sl%n)))
gz2 = fel_gamma0(beam)**2 / (1 + aw**2)
call fel_shortrange_ez (coll%efield, beam, sl, gz2, xks, ez)
esc_loss = 0
if (allocated(coll%long_esc)) esc_loss = -coll%long_esc(is) / m_electron

! Coupling coefficient for the energy exchange: fc/(sqrt(2)*m_electron), the V/m form of
! Genesis's fc/ks (see the module header for the unit relation). rpart then has units of
! 1/m and dgamma/dz is per meter directly.

rtmp = fel_und_coupling(und, 1) / (sqrt(2.0_rp) * m_electron)

! Grid bounds are loop-invariant: hoisted once per call (bit-for-bit: the same
! expression fel_grid_weights computed per particle).

ngrid_arr = wavefront_shape(wf)
gridmax = (ngrid_arr(1) - 1) * wf%dx / 2

if (nf == 1) then

  ! Single field (all pre-harmonic decks): the pre-field-set body, verbatim.

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

      ! Two live polarizations: the element couples to E_eff = conj(pol).E (manual
      ! sec:field). With one polarization pol = (1,0) and this branch never runs, so
      ! single-polarization arithmetic is untouched.

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

    ez_ip = ez(ip) + esc_loss     ! Short range plus the long-range loss (sec:spacecharge).
    call fel_runge_kutta (delz, xks, xku, btpar, rpart, ez_ip, gamma, theta)

    ! Back to the stored chart: pz from gamma (the subtraction is exact once p_mc is
    ! formed), z from tau = (phi0_new - theta)/ks with the updated beta.

    p_mc = sqrt(gamma**2 - 1)
    sl%pz(ip) = (p_mc - p0_mc) / p0_mc
    beta = p_mc / gamma
    sl%z(ip) = -beta * (phi0_new - theta) / xks
  enddo

else

  ! The field set: every field works the beam at once, BeamSolver::advance's structure.
  ! Per-field coupling fc(h) and record index held outside the particle loop, the
  ! per-field rpart (coupling times that field's interpolated amplitude, no phase) held
  ! through the RK stages, and the harmonic phase h*theta applied per stage inside
  ! fel_ode_multi (BeamSolver.cpp:24-30,150).

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
! The RK4 stage bookkeeping of BeamSolver::RungeKutta, VERBATIM -- the in-place stage
! algebra is kept exactly for bit identity, do not "clean up" (fel-physics.md sec-eom).
! ez is held fixed through the stages, as Genesis holds it.
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
! The longitudinal equations of motion, BeamSolver::ODE (fel-physics.md sec-eom),
! fundamental only. ez is the per-particle space-charge term in m_e c^2 per meter
! (short-range harmonics plus the long-range loss), zero when space charge is off.
! In that case the arithmetic is bit-identical to the pre-collective code because
! subtracting a literal zero is exact.
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
! The field-set form of fel_runge_kutta: identical RK4 stage bookkeeping
! (BeamSolver::RungeKutta), with the slope summing every field's contribution at that
! field's harmonic phase (fel_ode_multi). The single-field walk never calls this
! (fel_advance keeps the scalar path verbatim), so the pre-harmonic arithmetic is
! untouched by construction.
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
! The longitudinal equations of motion over a field set, BeamSolver::ODE verbatim:
! ctmp sums rpart(i) * exp(-i * rharm(i) * theta) over the fields at the CURRENT stage
! theta: the couplings held fixed through the stages, the harmonic phases live
! (BeamSolver.cpp:150). xks is the FUNDAMENTAL wavenumber (the theta equation is the
! fundamental's, and harmonics enter the slope only through ctmp).
!
! Input:
!   tgam, tthet -- real(rp): Stage values of gamma and the fundamental theta.
!   xks, xku    -- real(rp): Fundamental radiation and undulator wavenumbers [1/m].
!   btpar       -- real(rp): Parallel velocity beta_par.
!   rpart(:)    -- complex(rp): Per-field phasors times couplings.
!   rharm(:)    -- real(rp): Harmonic numbers as real factors (phase h*theta live).
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
! weights. Transcribed from Field::getLLGridpoint. wx is the weight of the LOWER x point.
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
! fel_grid_weights with the grid bounds PRECOMPUTED by the caller (the particle
! loops call this per particle: wavefront_shape and gridmax are loop-invariant).
! The arithmetic is character-identical to fel_grid_weights -- bit-for-bit --
! and the caller passes gridmax = (ngrid - 1) * wf%dx / 2 exactly as computed there.
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
! One slice's coherent-source summary (fel-physics.md sec-coherent-source; Tanaka PRAB 27,
! 030703 (2024), implemented from the paper). Computed: the source phasor S, EXACTLY as
! the deposit would (same theta, same part factors, post-advance particles: the
! normalization contract is SUM crsource = S), the charge-weighted centroid and
! central second moments (this port's extension for offset/mismatched/tilted beams),
! the Laguerre-Gauss order-0/1 sums g_q = SUM w_j Psi_q((x-x0)/sx, (y-y0)/sy) e^{-i
! theta_j} with Psi_q(u,v) = (1/sqrt(pi)) L_q(u^2+v^2) exp(-(u^2+v^2)/2) (Eqs 19-20,
! charge-normalized weights), and the guard metric (worst-plane excess kurtosis).
! Scale factors constant across the slice (scl_w) are NOT folded into S here. The
! caller applies them, so S matches the deposit's sum with part = sqrt(faw2)*w/gam.
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

! Pass 1: the source phasor (deposit-identical) and the charge moments.

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

! Pass 2: the LG sums (Eqs 19-20: principal widths, centered) and the guard's fourth
! moments. B in the same charge weighting, so the Eq 26 ratio is weight-consistent.

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
coh%m_ind = max(1.0_rp, wsum**2 / w2sum / max(1, beam%nbins))
coh%ok = .true.

end subroutine fel_coherent_prep

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_step (und, beam, sl, wf, ifld, delz, err_flag)
!
! Routine to advance field slice ifld (the rotated-record index, fel_field_index) one
! step: build the source from the (already pushed) particles of beam slice sl, propagate
! exp(K2 delz) in transverse Fourier space, add the source. Transcribed from
! FieldSolverFFT::advance and FFT, unfiltered path.
!
! The source scale, weighted, in V/m (derivation from Genesis's internal-unit form:
! fel-physics.md sec-field):
!
!   scl_w = fc * Z0 * sqrt(2) * c * delz / (4 * dgrid^2 * slice_spacing);  per particle scl_w*w_j
!
! in which the rest energy and the wavenumber have cancelled. Per particle
! part = sqrt(faw2(x,y))*scl_w*w_j/gamma, deposited as (sin theta + i cos theta)*part
! with the bilinear weights, added times 2 in real space after the transform pair.
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
type (fel_coherent_struct), optional :: coh    ! Coherent source (fel-physics.md sec-coherent-source):
real(rp), optional :: kappa                    !   present together when und%source_model says so.

real(rp) xks, dgrid, scl_w, part, theta, beta, wx, wy, gam, p_mc, p0_mc, s_t, c_t
complex(rp) cpart
complex(rp), allocatable :: crsource(:,:)    ! Per-call source accumulator: thread safe.
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

! The kernel is read-only here: a rebuild would race with concurrent slices. The caller
! initializes it serially (fel_field_kernel_init). A mismatch is a caller bug.

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

  ! The coherent-Gaussian source (fel-physics.md sec-coherent-source): the slice's whole
  ! phasor S, deposited as one analytic Gaussian at the phasor's charge centroid
  ! with covariance kappa^2 * (second-moment matrix). S carries the PHYSICAL shot
  ! noise through B(s), which the Fawley loading pins to <|B|^2> N_lambda = 1. The
  ! discrete sum is normalized to exactly S (the deposit's own contract, SUM crsource
  ! = SUM cpart), so edge truncation loses nothing silently.

  call coherent_gaussian_source ()

elseif (scl_w /= 0) then
  do ip = 1, sl%n
    call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
    if (.not. on_grid) cycle

    p_mc = p0_mc * (1 + sl%pz(ip))
    gam = sqrt(p_mc**2 + 1)
    beta = p_mc / gam

    ! The fundamental ponderomotive phase (xks1 = the FUNDAMENTAL wavenumber, passed by
    ! the caller so a harmonic field never reconstructs it from its own wavelength's
    ! rounding), scaled to harm*theta for a harmonic source: FieldSolver.cpp's
    ! theta = harm * particle.theta. harm = 1 leaves the phase untouched.

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

! Propagate and add: FFT, multiply the cached exp(K2 delz), inverse FFT, normalize,
! plus 2*crsource in real space.

call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_forward$, err);  if (err) return
wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) * fel_kernels(ik)%exp_k2
call wavefront_fft2 (wf%Ex(:,:,ifld), wf_fft_backward$, err);  if (err) return

if (allocated(wf%Ey)) then

  ! Two live polarizations: the element's source lands as pol * src on the pair
  ! (the exact dual of the kick's conj(pol).E read), and Ey diffracts identically.

  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * und%pol(1) * crsource
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_forward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) * fel_kernels(ik)%exp_k2
  call wavefront_fft2 (wf%Ey(:,:,ifld), wf_fft_backward$, err);  if (err) return
  wf%Ey(:,:,ifld) = wf%Ey(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * und%pol(2) * crsource
else
  wf%Ex(:,:,ifld) = wf%Ex(:,:,ifld) / real(ngrid*ngrid, rp) + 2 * crsource
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine coherent_gaussian_source ()
!
! Routine to deposit the coherent-Gaussian source (fel-physics.md sec-coherent-source):
! the slice's bunching phasor as a Gaussian of the phasor-weighted second moments
! scaled by kappa, normalized so the deposited sum equals the particle deposit's.
!-

subroutine coherent_gaussian_source ()

real(rp) det, a11, a12, a22, xg, yg, dx_c, dy_c, q, gsum, gridmax
real(rp), allocatable :: gk(:,:)
complex(rp) s_scaled
integer i, j

! Covariance kappa^2 * [[sxx, sxy],[sxy, syy]]. Its inverse for the quad form.

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
! Subroutine fel_field_kernel_init (ngrid, dgrid, ks, dz)
!
! Routine to (re)build the K2 kernel when the grid or the step changes, exactly as
! FieldSolverFFT::init builds it: dk = twopi/(ngrid*dgrid), integer offsets from the grid
! center, fftshift index mapping, K2 = -i*(dx^2+dy^2)*dk^2/(2*ks). Also built: the step
! propagator exp(K2*dz). K2 and dz are constant across each caller's whole
! parallel loop, so the ngrid^2 complex exponentials are paid once per (grid, ks, dz)
! rather than once per slice per step.
!
! MUST be called serially: fel_track_und_step and fel_track_interlude_genesis call it
! before their parallel slice loops. fel_field_step only reads the kernel and errors on a
! mismatch rather than rebuilding, so that nothing writes module state concurrently.
!
! Input:
!   ngrid -- integer: Grid points per side.
!   dgrid -- real(rp): Grid half-width [m].
!   ks    -- real(rp): Radiation wavenumber [1/m].
!   dz    -- real(rp): Step length [m].
!
! Output:
!   None directly: the module kernel cache (fel_kernels) gains an entry, and the
!   FFTW plans are warmed serially (the parallel loops then only execute).
!-

subroutine fel_field_kernel_init (ngrid, dgrid, ks, dz)

integer ngrid
real(rp) dgrid, ks, dz
real(rp) dk, shift, dx, dy
type (fel_kernel_struct), allocatable :: grow(:)
type (fel_kernel_struct), pointer :: kn
integer ix, iy, iix, iiy, ik
logical err

!

! Warm every thread's FFT plan cache from this serial context. No FFTW planner then
! runs concurrently with transform execution in the parallel slice loops that follow
! (wavefront_fft2_plan_threads' note). Runs even on a kernel-cache hit: the thread team
! may have grown since the plans were made. A planning failure reports there and then
! errors again, catchably, at the first transform.

call wavefront_fft2_plan_threads (ngrid, ngrid, err)

! One cache entry per wavelength: find this ks's entry (rebuilding it when the grid or
! the step changed, exactly the single-entry behavior each wavelength saw before), or
! append a new one. A single-field run lives entirely in entry 1.

if (.not. allocated(fel_kernels)) allocate (fel_kernels(0))

do ik = 1, size(fel_kernels)
  if (fel_kernels(ik)%ks /= ks) cycle
  if (fel_kernels(ik)%ngrid == ngrid .and. fel_kernels(ik)%dgrid == dgrid .and. &
      fel_kernels(ik)%dz == dz) return
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
shift = -0.5_rp * (ngrid - 1)

do iy = 0, ngrid-1
  dy = iy + shift
  do ix = 0, ngrid-1
    dx = ix + shift
    iiy = mod(iy + (ngrid+1)/2, ngrid)
    iix = mod(ix + (ngrid+1)/2, ngrid)
    kn%k2(iix+1, iiy+1) = cmplx(0.0_rp, -(dx*dx + dy*dy) * dk * dk / 2 / ks, rp)
  enddo
enddo

kn%exp_k2 = exp(kn%k2 * dz)

kn%ngrid = ngrid
kn%dgrid = dgrid
kn%ks = ks
kn%dz = dz

end subroutine fel_field_kernel_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_kernel_index (ngrid, dgrid, ks, dz) result (ik)
!
! Look up the cached kernel entry matching all four keys exactly. A miss returns 0.
! Read-only, so callable from inside the parallel slice loops.
!
! Input:
!   ngrid -- integer: Grid points per side.
!   dgrid -- real(rp): Grid half-width [m].
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
! Pure free-space diffraction of field slice ifld over the step dz: FFT, multiply
! the cached exp(K2 dz), inverse FFT, normalize. No source, no coupling factor, no
! undulator knowledge: this is the field half the unaveraged mode (fel-physics.md
! sec:unaveraged) shares with the averaged solver. The caller deposits its own source.
! The kernel cache must have been initialized serially (fel_field_kernel_init) for
! this grid, wavelength and step. A mismatch errors, exactly as fel_field_step.
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
! Subroutine fel_field_diag (wf, ifld, power, on_axis_intensity)
!
! Routine to compute the field diagnostics of field slice ifld from a field in V/m, the
! definitions matching Genesis's DiagField::calc expressed in physical units:
! power = sum|E|^2 * dgrid^2/(2*Z0) [W] (cell intensity |E|^2/(2*Z0) times cell area),
! on-axis intensity = |E(center)|^2/(2*Z0) [W/m^2]. Accumulation order: y outer, x inner.
!
! ifld is the raw record index. To report the field at time window position is, pass
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

! With a live second polarization, power and intensity are TOTALS (|Ex|^2 + |Ey|^2).
! Single-polarization lines never take this branch, keeping them bit-identical.

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
