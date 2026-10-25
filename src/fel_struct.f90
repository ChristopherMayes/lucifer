!+
! Module fel_struct
!
! Input and run-state structures of the FEL tracker, laid out as Tao lays out tao_struct.
! Every user-visible input struct is in this file, with its defaults in the declarations.
! One assembled-state struct, fel_run_struct, is passed explicitly by the library.
!
! The input structs have no allocatable or pointer components. So an input file sets
! them directly through the input file, as in Tao:
!
!   [run]
!     lattice = "line.bmad"
!     run.out_root = "run1"
!     bmad.radiation_damping_on = true
!     wake.radius = 2.5e-3
!     space_charge.nz = 2
!   /
!   [beam]
!     beam.n_particle = 8192
!     shot_noise = true
!   /
!   [field]
!     field.lambda0 = 1e-10
!     field.grid_n_pts = 255
!   /
!
! bmad_com and space_charge_com are Bmad's own globals, exposed directly as Tao's
! &tao_params exposes them. The parsing is in fel_input_mod. An embedding program may
! use it or fill the structs itself.
!
! Unlike Tao, there is no module-level singleton. Tao's super-universe s serves its
! command loop. A library wants explicit state, and an oscillator or scan driver wants
! several passes over one state. So fel_run_struct is passed as an argument everywhere.
!-

module fel_struct

use fel_track_mod
use fel_unaveraged_mod
use fel_import_mod
use fel_stats_mod

implicit none

! The sin^2 ramp of an unaveraged element, in periods at each end, where the input file leaves
! run.unaveraged_ramp_periods unset or sets it to 0 (fel-physics.md sec-unaveraged).

integer, parameter :: fel_ramp_periods_default$ = 2

!+
! Struct fel_slicing_struct
!
! Longitudinal discretization of the run, set in [run] as window.. The beam's
! slices and the field's z samples are one partition. So the window belongs to neither
! the beam nor the field and sits beside them.
!
! Slippage is an index rotation of the slice ring, with no interpolation. The rotation is
! exact only when the slice spacing is a whole number of wavelengths. n_wavelength is
! that whole number, and it is carried as an integer to every consumer. Dividing a
! spacing in metres by a wavelength gives 12 to within the last bit, and not 12. The
! device's exact bucket arithmetic refuses that.
!
! The window is set as a length or as a slice count, never both. current sets a flat
! current directly, as Genesis's &beam current does. With a window this is a flat
! time-dependent run. With no window it is the steady state of one slice, so there is no
! separate steady-state switch. Unset, the current comes from the bunch's z structure.
!-

type fel_slicing_struct
  real(rp) :: window_length = 0      ! Time window [m]. 0 = from the bunch, or one slice.
  integer :: n_wavelength = 0        ! Wavelengths per slice (Genesis's sample). An integer:
                                     !   slippage rotates the ring by one. 0 derives it from
                                     !   the gain. It is 1 in the steady state, where the
                                     !   spacing sets the current.
  integer :: n_slice = 0             ! 0 = from window_length. Both set is refused.
  real(rp) :: current = 0            ! Flat current [A]. 0 = from the bunch's z structure.
end type

!+
! Struct fel_global_struct
!
! Run-level switches and names, set in [run] as run.. The analog of
! tao_global_struct. Everything here concerns one run of the tracker. The physics is
! described by the lattice, beam_init and wavefront_init.
!-

type fel_global_struct
  character(400) :: out_root = 'fel_track'   ! Output file root.
  character(16) :: interlude_model = 'bmad'  ! 'bmad' (the seam) or 'genesis' (transcribed).
  ! Transverse transport inside averaged FEL elements. 'bmad' is Bmad's own periodic-wiggler
  ! kernel, the production model. 'genesis' is Genesis4's maps verbatim, for validation
  ! only: the comparisons against Genesis4 need its transport. So it is a run switch rather
  ! than a tracking method, and no production run sets it.
  character(16) :: transport_model = 'bmad'
  ! The two numbers of the unaveraged model. They were per-element attributes registered
  ! by this program, and every other Bmad program then refused to parse the lattice. So
  ! they are run switches, the same for every element. Averaged and unaveraged elements
  ! still mix per element through Bmad's fel_method attribute.
  ! Steps per undulator period. Below 10 the push does not converge (fel-physics.md
  ! sec-unaveraged). A smaller number is refused rather than rounded up.
  integer :: unaveraged_steps_per_period = 20
  ! Periods over which the field ramps at each end. 0 keeps the default of 2 and -1 asks
  ! for a hard edge. A hard edge reintroduces the handoff phase jump, so it is never the
  ! result of leaving the setting at 0.
  integer :: unaveraged_ramp_periods = fel_ramp_periods_default$
  ! The FEL source model (fel-physics.md sec-coherent-source). 'deposit' is the standard
  ! per-particle scatter. 'coherent' is the SIMPLEX-hybrid coherent-Gaussian source
  ! (Tanaka, PRAB 27, 030703 (2024)). It drops the spatially incoherent artifact, keeps
  ! the physical shot noise in the slice bunch factor B(s), and makes the transverse shape
  ! a guarded Gaussian.
  character(16) :: source_model = 'deposit'
  ! The angular filter on the source term (fel-physics.md sec-convergence), transcribed
  ! from Genesis4's source_filter, added in release 4.6.12 with its FFT solver. The
  ! transformed source is multiplied by a sigmoid in normalized transverse spatial
  ! frequency before it reaches the field. This suppresses the wide-angle emission of the
  ! point-like beamlets, and the field's own propagation is untouched. On by default
  ! because what it removes is emission no physical beam radiates (the four-machine
  ! results are in fel-physics.md sec-source-filter). Turn it off to compare against a
  ! code that has no such filter.
  logical :: source_filter = .true.
  ! The sigmoid's edge as an angle [rad]. Unset, the run derives it from the beam and the
  ! gain. It is the larger of four mode diffraction angles and the angle at which the
  ! resonant wavelength red-shifts by rho. The run prints both and which one won.
  ! Cutting into the mode loses real radiation, while leaving artifact in only weakens
  ! the filter. So the default errs wide. Setting the angle bypasses
  ! source_filter_tolerance.
  real(rp) :: source_filter_angle = 0
  ! The largest fraction of source intensity the filter may remove anywhere inside the
  ! protected angle, the larger of the two derived angles. The edge is placed outside
  ! that angle by as much as the sigmoid's width needs to hold this. A half-amplitude
  ! edge on the protected angle passes a quarter of the intensity there. So 0.75 puts
  ! the half-amplitude point on the protected angle, to the bit. Strictly between 0 and 1.
  real(rp) :: source_filter_tolerance = 0.75_rp
  ! The sigmoid's edge in Genesis4's own units, xcut and ycut: the shifted grid index
  ! over ngrid. So 1 is twice the Nyquist index, and the angle is xcut lambda/dx. For
  ! validation only: they place the edge where Genesis4 places it. An angle is
  ! independent of the grid and these are not. Setting them with the angle is refused.
  real(rp) :: source_filter_xcut = 0
  real(rp) :: source_filter_ycut = 0
  ! The sigmoid's width as a fraction of its edge, Genesis4's sigmoid. Small is a sharp
  ! edge. Genesis4's own default of 1 leaves the sigmoid at 0.73 on axis. That attenuates
  ! the coherent source as much as the wide angles and delays saturation
  ! (fel-physics.md sec-source-filter).
  real(rp) :: source_filter_width = 0.05_rp
  ! Self-test of the filter check, as fp32_mutate is for the FP32 check: apply the
  ! sigmoid to the propagated field instead of to the source. Both suppress wide angles,
  ! so the run still completes and looks reasonable. Only the comparison against
  ! Genesis4 says which one the code did.
  logical :: source_filter_mutate = .false.
  ! The tracking window as element locators (lat_ele_locator syntax). The names are Tao's
  ! (tao_beam_init carries track_start and track_end), and the role is Genesis's zstop.
  ! Blank is the whole line. The schedule (slippage, autophasing, break geometry) is always
  ! built on the full lattice. So a windowed run composes exactly with the full run:
  ! [1,k] and then [k+1,end] from its dumps reproduce [1,end] bit for bit.
  character(60) :: track_start = '', track_end = ''
  ! What beam.file and field.file are for. Off, they initialize a fresh run from sliced
  ! inputs, prepared or ramped. Only the standard records are read, phi0 and the slippage
  ! residual start at zero, and a checkpoint group in either file is ignored.
  ! On, they continue a run. Both files must be frames of one eligible boundary and carry
  ! the checkpoint group (fel_write_checkpoint_group). The reader rebuilds phi0, every z
  ! and the residual from them. track_start must open the walk at the element after the
  ! one the checkpoint completed. The mode is never inferred from what the files carry.
  logical :: continuation = .false.
  ! The comb, with the name and meaning of Bmad's bunch_track_struct%ds_save and Tao's
  ! comb_ds_save: the minimum z advance between per-record stats rows.
  !   < 0: No per-record rows. Element ends, dumps and the final state remain.
  !   = 0: A row at every record position. This is the default here, since the per-record
  !        arrays are the primary stats output. Tao's default is -1, with the same meaning.
  !   > 0: A row when z has advanced comb_ds_save past the last row. Element ends always.
  real(rp) :: comb_ds_save = 0
  ! The interlude step. Tao's comb has two halves. save_a_bunch_step selects among the
  ! positions a tracker already reaches. tao_lattice_calc_mod's beam-track loop cuts an
  ! element into n_slice pieces with element_slice_iterator, where the tracker would
  ! otherwise cross it in one map. comb_ds_save above is the first half and this is the
  ! second. It is a setting of its own for two reasons. The comb defaults to zero here,
  ! so it names no length to cut by. And the comb must not change what the tracking does.
  ! The default of zero tracks an interlude whole.
  !   <= 0: One piece, the element.
  !   >  0: Pieces of about this length [m]. Stats rows and frames then land inside a
  !         drift as they already land inside a wiggler.
  real(rp) :: interlude_ds_step = 0
  character(60) :: dump_beam_at(40) = ''     ! Element locators for mid-run beam dumps.
  character(60) :: dump_field_at(40) = ''    ! Element locators for mid-run field dumps.
  ! Write the beam and the field at every comb position as a frame series. This gives a
  ! visualization z resolution inside the undulators and not only at their ends. A frame
  ! and its stats row share an index (docs/reading-output.md).
  logical :: dump_at_comb = .false.
  ! The slices a frame carries, in window order. The default pair is the whole window.
  ! A range keeps part of the bunch: a sample = 1 run has thousands of slices, and a
  ! whole-window field frame is gigabytes. The element-end, initial and final dumps are
  ! restart points and are unaffected.
  integer :: dump_slice_first = 1
  integer :: dump_slice_last = -1            ! -1 is the last slice of the window.
  ! Write the field's own reductions per record into the stats file. These are the
  ! projections a picture is drawn from, computed once rather than from raw frames.
  logical :: dump_reduced = .false.
  integer :: ran_seed = 12345                ! The one RNG seed (generation, import, noise).
  logical :: write_diag = .false.            ! The Genesis-comparison text diag file (large).
  logical :: write_initial = .false.         ! Dump the initial state before tracking.
  logical :: load_only = .false.             ! Build the initial state, dump it, stop.
  logical :: keep_escaped_field = .false.    ! Keep the escaped-slice bank file.
  logical :: migrate = .false.               ! Slice migration (fel-physics.md sec-migration).
  logical :: migrate_check = .false.         ! Migration's bunching-invariance instrument.
  logical :: reference_run = .false.         ! No FEL interaction: Bmad tracks everything.
  ! Provenance detail in the meta/ group of stats.h5. On, meta/ also records the login
  ! name and the working directory. That is useful in a lab notebook and a leak outside
  ! one. Off by default because a stats file travels: with a paper, to a collaborator,
  ! beside a figure. Genesis always records them.
  logical :: record_environment = .false.
  ! The FP32 lockstep check (fel_fp32_mod): 'off', 'lockstep' or 'freerun'. With
  ! 'lockstep' the FP32 twin is rebuilt from FP64 every step. With 'freerun' its
  ! longitudinal state carries, so the error compounds and its rate shows. For checks
  ! only. No production run sets it.
  character(16) :: fp32_check = 'off'
  logical :: fp32_mutate = .false.     ! The check's self-test: coarsen the residual.
  ! The device backend (fel_device_mod): 'off', or 'metal' for the Apple Silicon backend
  ! on a build that carries it. Alone, it runs the averaged FEL step resident on the
  ! device. With fp32_check set, the device takes the lockstep twin's role and the FP64
  ! run is untouched. Anything the backend does not cover is refused at setup or first
  ! use. It is never run on the CPU instead.
  character(16) :: device = 'off'
  ! Per-pass device timing (lucifer_device.h). The backend's own busy seconds cover one
  ! command buffer, which holds every step between two host touches. So they cannot give
  ! the time of one pass. On, each pass of a step is encoded into its own compute encoder
  ! with a timestamp attachment, which stage-boundary counter sampling can read. The run
  ! then reports the seconds per pass. The extra encoder boundaries slow the run, so this
  ! is off by default and timings are taken with it off (docs/performance.md perf-device).
  logical :: device_timing = .false.
  integer :: device_dep_mutate = 0     ! Shift the deposit's fixed-point scale, in bits.
end type

!+
! Struct wavefront_init_struct
!
! The radiation starting condition ([field]), the analog of beam_init_struct.
! A wavefront is an optical object: a transverse grid, a wavelength and a longitudinal
! spacing, and nothing else longitudinal. The time window that sets that spacing belongs
! to the FEL interaction and lives in fel_slicing_struct. harmonics requests the field
! set. A field.file import overrides the seed.
!-

type wavefront_init_struct
  real(rp) :: lambda0 = 0            ! Resonant wavelength [m]. Required for generation.
  ! The transverse grid (Genesis ngrid and dgrid). Zero derives it from the beam at
  ! setup, with fel_cells_per_sigma$ and fel_widths_per_sigma$ below. One may be set and
  ! the other derived against it.
  integer :: grid_n_pts = 0          ! Transverse grid points per side.
  real(rp) :: grid_half_width = 0    ! Transverse half width [m].
  real(rp) :: seed_power = 0         ! Gaussian seed power [W]. 0 = dark start.
  real(rp) :: seed_waist_size = 0    ! Seed intensity 1/e^2 radius [m].
  character(1) :: seed_polarization = 'x'   ! 'x' or 'y'.
  ! The field set: harmonics(1) must be 1 (the fundamental). Further entries are
  ! harmonic numbers in increasing order, 0 = unused.
  integer :: harmonics(9) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
end type

! The transverse grid derived from the beam when the input file sets none
! (docs/startup-noise.md, Recommendations). Cells of a seventh of the rms beam size
! resolve the beam and the mode. Finer cells raise the wide-angle emission of the point
! beamlets without moving the mode power. A half width of nine beam sizes contains the
! mode, and above that the half width does not enter. fel-physics.md sec-field gives
! both constants on three machines a hundred times apart in wavelength.

real(rp), parameter :: fel_cells_per_sigma$ = 7
real(rp), parameter :: fel_widths_per_sigma$ = 9

!+
! Struct fel_chamber_wake_init_struct
!
! The chamber-wake description ([run] wake.), with Genesis &wake names. A
! plain mirror of the configuration fields of fel_wake_struct. That struct carries
! allocatable state, which an input mirror need not. fel_setup copies these in.
!-

type fel_chamber_wake_init_struct
  logical :: on = .false.
  ! Which chamber-wake implementation runs. 'genesis' is the transcribed solver at slice
  ! granularity. It convolves the kernels with the weighted slice currents. A second
  ! implementation would be a new value of this setting.
  character(16) :: model = 'genesis'
  real(rp) :: loss = 0             ! External loss [eV/m].
  real(rp) :: radius = 2.5e-3_rp   ! Chamber radius, or half gap if flat [m].
  real(rp) :: conductivity = 0     ! DC conductivity [1/(Ohm m)]. 0: no resistive wake.
  real(rp) :: relaxation = 0       ! AC relaxation distance c*tau [m].
  logical :: roundpipe = .true.    ! Round chamber. False: flat (parallel plates).
  character(8) :: material = ''    ! 'CU' or 'AL' shortcut for conductivity+relaxation.
  real(rp) :: gap = 0              ! Undulator gap [m]. 0: no geometric wake.
  real(rp) :: lgap = 1             ! Period of the gaps [m].
  real(rp) :: hrough = 0           ! Roughness amplitude [m]. 0: no roughness wake.
  real(rp) :: lrough = 1           ! Roughness period [m].
  ! Check instrument: export the transcribed kernels for building z_long tables.
  character(400) :: write_kernels = ''
end type

! Space charge needs no init mirror. fel_space_charge_struct is pure scalars and is read
! directly in [run] as space_charge..

!+
! Struct fel_beam_init_param_struct
!
! The beam-side scalars of [beam]. They sit beside Bmad's beam_init, which
! describes the bunch, and imp, the resampler. They are the source and output files,
! the quiet-start settings and the check instruments the benchmark sets.
!-

type fel_beam_init_param_struct
  character(400) :: beam_file = ''        ! An FEL beam already in slices (with field.file(1)).
  character(400) :: write_genesis_dist = ''  ! Write the bunch as a Genesis4 distribution file.
  character(400) :: write_openpmd_file = ''  ! Write the bunch as openPMD-beamphysics.
  ! How a bunch becomes slices (fel-physics.md sec-loading). The bunch is beam_init's,
  ! generated or read through beam.position_file. It is binned by arrival time into
  ! the slices of slicing.
  ! "sample" draws the same number of beamlets in every slice from what fell into it,
  ! weighted by the slice charge, so the current follows the bunch. This is Genesis's
  ! &beam with &profile and its importdistribution in one. beam.n_particle is the
  ! count per slice after phase copies.
  ! "keep" keeps every real particle as a beamlet of beamlet_size copies at
  ! weight/beamlet_size. The counts follow the charge, so the correlations of a
  ! start-to-end bunch survive whole.
  character(16) :: load_mode = 'sample'
  ! The two switches of the load. quiet_start makes beamlets of copies with phases spread
  ! over 2 pi about each particle's own. The load then has no bunching at any harmonic
  ! below beamlet_size. shot_noise imposes the physical level on whatever load there is,
  ! by independent phasors per group. With both off the beam is taken as it came.
  logical :: quiet_start = .true.
  integer :: beamlet_size = 8             ! Copies per beamlet. A loading parameter only.
  logical :: shot_noise = .false.         ! Impose physical shot noise on the load.
  ! Check instruments, set by the benchmark. They are not physics inputs.
  logical :: split_weights = .false.      ! Coincident w/3 + 2w/3 copies after loading.
  logical :: swap_beam_xy = .false.       ! Swap (x,px) <-> (y,py) after generation.
  logical :: gen_test_weights = .false.   ! Alternate beamlet weights 0.25x/1.75x.
  logical :: resample_split_weights = .false.  ! Split-weight copies before the resample.
end type

!+
! Struct fel_run_struct
!
! The assembled state of one run: the parsed lattice, the beam, the field set, the
! collective-effects state, the setup schedule and the diagnostics state. It also holds
! a copy of every input struct, so the resolved inputs are one object. The Meta/
! provenance echo and any embedding caller read them from here. Passed explicitly, with
! no singleton.
!-

type fel_run_struct
  ! The resolved inputs.
  type (fel_global_struct) :: global
  type (wavefront_init_struct) :: winit
  type (fel_slicing_struct) :: slicing
  type (fel_chamber_wake_init_struct) :: chamber_wake
  type (fel_space_charge_input_struct) :: space_charge
  type (beam_init_struct) :: beam_init
  type (fel_resample_param_struct) :: resample
  type (fel_beam_init_param_struct) :: bparam
  character(400) :: lat_file = ''
  character(400) :: field_file(9) = ''
  ! The assembled state.
  type (lat_struct) :: lat
  type (fel_beam_struct) :: fbeam
  type (fel_field_struct), allocatable :: ffield(:)
  type (fel_collective_struct) :: coll
  type (fel_fp32_struct) :: fp32
  type (fel_device_struct) :: dev
  type (fel_stats_struct) :: stats
  type (fel_unavg_struct) :: ustate
  ! The schedule (fel_setup): per tracked element.
  type (fel_und_struct), allocatable :: und_of(:)
  ! The angle at which the stats split the field power, whether or not the filter is on.
  ! It is the larger of four mode diffraction angles and the rho angle, derived from the
  ! beam and the first FEL element (fel-physics.md sec-source-filter). Inside it lies what
  ! can couple to the mode. Outside it lies the wide-angle emission of the point beamlets.
  real(rp) :: split_angle = 0
  character(24) :: split_origin = ''       ! Where split_angle came from, for the report.
  integer, allocatable :: fel_mode(:), fel_spp(:)
  real(rp), allocatable :: fel_ramp(:)
  real(rp), allocatable :: ele_slip(:)     ! Slippage after each element's last step [wavelengths].
  real(rp), allocatable :: fel_zoff(:)     ! The z_offset off-phase knob per element [m].
  real(rp), allocatable :: light_corr(:)   ! Chord-vs-arc correction on a break's last element [m].
  logical, allocatable :: is_fel(:)
  ! Space charge per element: ele%space_charge_method = slice and Bmad's master switch
  ! bmad.csr_and_space_charge_on, resolved once at setup for the walk to read.
  logical, allocatable :: sc_here(:)
  logical, allocatable :: dump_beam_here(:), dump_field_here(:)
  integer :: i_start = 1, i_end = 0        ! The resolved tracking window [elements].
  type (fel_checkpoint_struct) :: ckpt     ! Where a continuation's checkpoint was written.
  integer :: dump_is1 = 1, dump_is2 = 0    ! The resolved frame slice range [slices].
  ! Run facts.
  real(rp) :: gamma0 = 0                   ! From the lattice e_tot.
  real(rp) :: gamma0_ref = 0               ! fel_gamma0(fbeam), the walk's reference.
  real(rp) :: phase_rate = 0               ! 2pi/(2 gamma0_ref^2 lambda) [rad/m].
  real(rp) :: ks = 0                       ! 2pi/lambda.
  real(rp) :: z_now = 0
  integer :: nslice = 0, n_harm = 1
  logical :: two_pol = .false., any_unavg = .false., timerun = .false.
  ! Running counters and ledger terms.
  real(rp) :: u_spont_cum = 0, e_rad_cum = 0
  real(rp) :: charge_dropped_tot = 0, b_dev_max = 0
  integer :: n_moved_tot = 0
  ! Diagnostics state.
  type (fel_slice_diag_struct), allocatable :: bdiag_arr(:)
  real(rp), allocatable :: fpow_arr(:), fonax_arr(:)
  real(rp), allocatable :: e_rad_slice(:), rad_kick(:,:)
  ! The interlude pieces, counted by the record-count precompute that replays the walk.
  ! The pre-run header reports them without a second loop or a second rule.
  integer :: n_int_piece = 0            ! Pieces the interludes will take.
  integer :: n_int_ele = 0              ! Interlude elements that take more than one.
  integer :: nrec_stats = 0, nend_stats = 0
  real(rp) :: z_last_rec = -1e30_rp        ! The comb's last-row position.
  ! Escaped-field bank state, one slot per field.
  integer, allocatable :: n_banked(:)
  real(rp), allocatable :: bank_z(:,:), bank_pms(:,:,:)
  integer(hid_t), allocatable :: esc_id(:)
  ! Whole-window wake scratch (element sr wakes across the window).
  type (bunch_struct) :: wake_bunch
  real(rp), allocatable :: wake_beta0(:)
  ! Open output units (0 = not open).
  integer :: iu_diag = 0, iu_ledger = 0, iu_wake = 0
end type

contains

!------------------------------------------------------------------------------
!+
! Function fel_si_str (value, unit) result (str)
!
! Routine to format a value for display. The SI prefix is chosen so the mantissa lies in
! [1, 1000), with three decimals. The result is right-justified to a fixed width, so a
! column of values lines up (4.230 kW above 105.000 GW). This is display only. Stdout is
! for humans, and the files carry full precision (docs/user-guide.md).
!
! Values outside the prefix range print in es11.3 with the bare unit. Exact zero prints
! as 0.000 with the bare unit.
!
! Input:
!   value -- real(rp): The value, in the unit's own base (W, J, m, ...).
!   unit  -- character(*): The unit symbol, appended after the prefix.
!
! Output:
!   str   -- character(14): The formatted value, right-justified. Callers wanting it
!              inline use trim(adjustl(...)).
!-

function fel_si_str (value, unit) result (str)

real(rp) value, av, mant
character(*) unit
character(14) str
character(14) num
integer ip

! The table spans yocto to peta. That covers everything this program reports: a wake
! run's pulse energy is attojoules and a saturated pulse is millijoules. Index 9 is the
! bare unit.

character(1), parameter :: pfx(14) = &
      ['y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', ' ', 'k', 'M', 'G', 'T', 'P']

!

av = abs(value)
if (av == 0) then
  write (num, '(f8.3, 1x, a)') 0.0_rp, trim(unit)
  str = adjustr(num)
  return
endif

ip = 9 + floor(log10(av) / 3.0_rp)
if (ip < 1 .or. ip > 14) then
  write (num, '(es11.3, 1x, a)') value, trim(unit)
else
  mant = value / 10.0_rp**(3 * (ip - 9))
  write (num, '(f8.3, 1x, 2a)') mant, pfx(ip), trim(unit)
endif
str = adjustr(num)

end function fel_si_str

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_comb_take (comb_ds_save, z, z_last, at_end) result (take)
!
! Routine to apply the comb rule. The rule is Bmad's bunch_track_struct%ds_save, with
! save_a_bunch_step's guards and Tao's comb_ds_save note "< 0 => No comb calculated":
!   comb < 0: Element ends only.
!   comb = 0: A row at every record position.
!   comb > 0: A row when z has advanced comb_ds_save past the last row.
!
! Note: There is one difference from Bmad. An element end is always a row, whatever the
! comb. This lets the stats file carry one record axis with a boolean mask
! (coords/at_element_end), instead of a second axis and a copy of every element-end
! quantity. Where Bmad's comb < 0 drops the comb, here the element ends remain. They are
! the positions the evaluated bunch_params live on.
!
! z_last is updated when the row is taken. The walk calls this routine live. The setup's
! nrec precompute replays it with the same z arithmetic, so the stats arrays are sized
! exactly in every mode.
!
! Input:
!   comb_ds_save -- real(rp): The comb setting (see above).
!   z            -- real(rp): Current position [m].
!   z_last       -- real(rp): Position of the last row taken [m].
!   at_end       -- logical: True at an element end.
!
! Output:
!   z_last       -- real(rp): Updated when the row is taken.
!   take         -- logical: True when a stats row is due.
!-

function fel_comb_take (comb_ds_save, z, z_last, at_end) result (take)

real(rp) comb_ds_save, z, z_last
logical at_end, take

!

if (at_end) then
  take = .true.
elseif (comb_ds_save < 0) then
  take = .false.
else
  take = (z >= z_last + comb_ds_save)
endif
if (take) z_last = z

end function fel_comb_take

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_interlude_pieces (ds_step, length) result (n_piece)
!
! Routine to return the number of pieces an interlude element is tracked in. The walk
! cuts the element this many ways. The setup's record-count precompute replays the same
! count, so the stats arrays stay exactly sized.
!
! The expression is Tao's own, from the beam-track loop of tao_lattice_calc_mod,
! including the 1.01 factor. The factor keeps an element of exactly n whole steps from
! rounding down to n-1 pieces.
!
! Input:
!   ds_step -- real(rp): run.interlude_ds_step. Zero or negative gives one piece.
!   length  -- real(rp): The element's length [m].
!
! Output:
!   n_piece -- integer: Pieces, at least one.
!-

function fel_interlude_pieces (ds_step, length) result (n_piece)

real(rp) ds_step, length
integer n_piece

!

n_piece = 1
if (ds_step <= 0 .or. length <= 0) return
n_piece = max(1, int(1.01_rp * length / ds_step))

end function fel_interlude_pieces

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_interlude_wake (ele) result (wake_src)
!
! Routine to return the element whose short-range wake acts on an interlude. The walk
! routes a wake-carrying element through the whole-window path and refuses to cut it into
! pieces. The setup's record-count precompute replays both decisions. The two must agree
! or the stats arrays are mis-sized.
!
! Bmad's own resolver, pointer_to_wake_ele, finds the wake on the element or on exactly
! one slave of a split lord. With bmad.sr_wakes_on off, track1_bunch applies no
! wake, so no element carries one here. A zero-length element whose wake scales with
! its length (Bmad's default) kicks by exactly zero and carries none either. A
! long-range wake has no scale_with_length and acts at zero length, so it is kept.
!
! Input:
!   ele      -- ele_struct: The element.
!
! Output:
!   wake_src -- ele_struct, pointer: The element whose wake acts, or null.
!-

function fel_interlude_wake (ele) result (wake_src)

type (ele_struct), target :: ele
type (ele_struct), pointer :: wake_src

!

wake_src => null()
if (.not. bmad_com%sr_wakes_on) return
wake_src => pointer_to_wake_ele(ele)
if (.not. associated(wake_src)) return

if (ele%value(l$) == 0 .and. wake_src%wake%sr%scale_with_length) then
  if (.not. bmad_com%lr_wakes_on) then
    wake_src => null()
  elseif (.not. allocated(wake_src%wake%lr%mode)) then
    wake_src => null()
  elseif (size(wake_src%wake%lr%mode) == 0) then
    wake_src => null()
  endif
endif

end function fel_interlude_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_element_inert (ele, wake_src) result (inert)
!
! Routine to say whether the walk skips an element. The walk and the setup's
! record-count precompute both use it.
!
! Only a zero-length element that cannot act is skipped. A thin kicker, a multipole, a
! patch, a zero-length cavity or any other element whose key carries a map is tracked
! whatever its length. Bmad tracks it, and skipping it would drop its kick with no
! message. A marker, drift, pipe, instrument or monitor of zero length is skipped when
! it has no kick, no multipole, no aperture and no wake that can act. Bmad tracks such
! an element as the identity.
!
! Input:
!   ele      -- ele_struct: The element.
!   wake_src -- ele_struct, pointer: fel_interlude_wake's answer for it.
!
! Output:
!   inert    -- logical: True when the walk skips the element.
!-

function fel_element_inert (ele, wake_src) result (inert)

type (ele_struct) ele
type (ele_struct), pointer :: wake_src
logical inert

!

inert = .false.
if (ele%value(l$) /= 0) return
if (associated(wake_src)) return

select case (ele%key)
case (marker$, beginning_ele$, null_ele$, fiducial$, floor_shift$, fork$, photon_fork$, &
      drift$, pipe$, instrument$, monitor$, detector$)
case default
  return
end select

if (ele_has_nonzero_kick(ele)) return
if (associated(ele%a_pole) .or. associated(ele%a_pole_elec)) return
if (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0) return
if (ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) return

inert = .true.

end function fel_element_inert

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_refuse_two_pol_combinations (run, err_flag)
!
! Routine to refuse the combinations that two live polarizations are not validated with.
! run%two_pol is decided twice. It is set at lattice setup from the seed's polarization
! and the element tilts. It is set again when a field file carries both planes. The
! refusals are made at both points by this routine.
!
! Input:
!   run      -- fel_run_struct: Run state with %two_pol, %n_harm and %global set.
!
! Output:
!   err_flag -- logical: Set True when a combination is refused. False otherwise.
!-

subroutine fel_refuse_two_pol_combinations (run, err_flag)

type (fel_run_struct) run
logical err_flag
character(*), parameter :: r_name = 'fel_refuse_two_pol_combinations'

!

err_flag = .false.
if (.not. run%two_pol) return

if (run%n_harm > 1) then
  call out_io (s_error$, r_name, 'HARMONIC FIELDS WITH TWO LIVE POLARIZATIONS ARE NOT VALIDATED', &
                                 'TOGETHER YET; RUN ONE OR THE OTHER.')
  err_flag = .true.;  return
endif
if (run%global%source_model == 'coherent') then
  call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH TWO LIVE POLARIZATIONS IS NOT IN V1.')
  err_flag = .true.;  return
endif

end subroutine fel_refuse_two_pol_combinations

end module fel_struct
