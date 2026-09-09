!+
! Module fel_struct
!
! The input and run-state structures of the FEL tracker, laid out the way Tao lays out
! tao_struct: every user-visible input struct in one file, defaults in the declarations,
! and one assembled-state struct (fel_run_struct) that the library passes explicitly.
!
! The input structs are namelist-safe (no allocatable or pointer components), so an
! input file sets them directly, Tao-style:
!
!   &fel_params
!     lat_file = "line.bmad"
!     global%out_root = "run1"
!     bmad_com%radiation_damping_on = T
!     chamber_wake%radius = 2.5e-3
!     space_charge%nz = 2
!   /
!   &fel_beam_init
!     beam_init%n_particle = 8192
!     shot_noise = T
!   /
!   &fel_wavefront_init
!     wavefront_init%lambda0 = 1e-10
!     wavefront_init%grid_n_pts = 255
!   /
!
! (bmad_com and space_charge_com are Bmad's own globals, exposed directly as Tao's
! &tao_params exposes them. The parsing lives in fel_input_mod, which an embedding
! program may reuse or skip by filling the structs itself.)
!
! One deliberate deviation from Tao: no module-level singleton. Tao's super-universe s
! serves its command loop. A library wants explicit state, and an oscillator or scan
! driver wants several passes over one state. fel_run_struct passes as an argument
! everywhere.
!-

module fel_struct

use fel_track_mod
use fel_unaveraged_mod
use fel_import_mod
use fel_stats_mod

implicit none

!+
! Struct fel_slicing_struct
!
! The longitudinal discretization of the run, exposed in &fel_params as slicing%. One
! window, one definition: the beam's slices and the field's z samples are the same
! partition, so it belongs to neither the beam nor the field and sits beside them.
!
! Slippage is an index rotation of the slice ring, exact and free of interpolation,
! which holds only when the spacing is a whole number of wavelengths. n_wavelength is
! that whole number, and it is carried as an integer from here to every consumer:
! recovering it by dividing a spacing in metres by a wavelength gives 12 to the last
! bit rather than 12, which the device's exact bucket arithmetic then refuses
! (FINDINGS 7.52).
!
! The window is stated as a length or as a slice count, never both. current states a
! flat current directly, which is Genesis's &beam current: with a window it is a flat
! time-dependent run, and with none it is the steady state of one slice, so no separate
! steady-state switch exists. Unset, the current comes from the bunch's own z structure.
!-

type fel_slicing_struct
  real(rp) :: window_length = 0      ! Time window [m]. 0 = from the bunch, or one slice.
  integer :: n_wavelength = 0        ! Wavelengths per slice (Genesis's sample). Integer by
                                     !   construction: slippage rotates the ring by one.
                                     !   0 derives it from the gain, and is 1 in the steady
                                     !   state, where the spacing states the current.
  integer :: n_slice = 0             ! 0 = from window_length. Both set is refused.
  real(rp) :: current = 0            ! Flat current [A]. 0 = from the bunch's z structure.
end type

!+
! Struct fel_global_struct
!
! The run-level switches and names, exposed in &fel_params as global%... (the
! tao_global_struct analog). Everything here is about one run of the tracker. The
! physics description lives in the lattice, beam_init and wavefront_init.
!-

type fel_global_struct
  character(400) :: out_root = 'fel_track'   ! Output file root.
  character(16) :: interlude_model = 'bmad'  ! 'bmad' (the seam) or 'genesis' (transcribed).
  ! Transverse transport inside averaged FEL elements: 'bmad' is Bmad's own periodic-wiggler
  ! kernel, the production model, and 'genesis' is Genesis4's maps verbatim. The second is
  ! validation-internal: the comparison tiers need transcription-level transport, so it is a
  ! run switch rather than a tracking method, and no production run sets it.
  character(16) :: transport_model = 'bmad'
  ! The unaveraged model's two numbers. They were per-element attributes registered by
  ! this program, which every other Bmad program then refused to parse, so they are run
  ! switches selected uniformly. Averaged and unaveraged elements still mix per element,
  ! through the fel_method attribute Bmad carries.
  ! Steps per undulator period. Below 10 the push does not converge (fel-physics.md
  ! sec-unaveraged), so a smaller number is refused rather than rounded up.
  integer :: unaveraged_steps_per_period = 20
  ! Periods over which the field ramps at each end. 0 is not a hard edge: a silent hard
  ! edge reintroduces the handoff phase jump by omission, so the hard edge is asked for
  ! with -1 and 0 keeps the default of 2.
  integer :: unaveraged_ramp_periods = 2
  ! The FEL source model (fel-physics.md sec-coherent-source). 'deposit' is the standard
  ! per-particle scatter (bit-for-bit unchanged). 'coherent' is the
  ! SIMPLEX-hybrid coherent-Gaussian source (Tanaka, PRAB 27, 030703 (2024)): the
  ! spatially incoherent artifact is dropped, the slice bunch factor B(s) keeps the
  ! physical shot noise, and the transverse shape is a guarded Gaussian.
  character(16) :: source_model = 'deposit'
  ! The angular filter on the source term (fel-physics.md sec-convergence), transcribed
  ! from Genesis4's source_filter, which release 4.6.12 added with its FFT solver. The
  ! transformed source is multiplied by a sigmoid in normalized transverse spatial
  ! frequency before it reaches the field, so the wide-angle emission of the point-like
  ! beamlets is suppressed while the field itself propagates untouched. On by default:
  ! it is measured on four machines, and what it removes is emission no physical beam
  ! radiates. Turn it off to compare against a code that has no such filter, which is
  ! what the comparison tiers do.
  logical :: source_filter = .true.
  ! The sigmoid's edge as an angle [rad]. Unset, the run derives it from the beam and the
  ! gain: the larger of four mode diffraction angles and the angle at which the resonant
  ! wavelength red-shifts by rho, printed with both and with which one won. Cutting into
  ! the mode loses real radiation while leaving artifact in only weakens the filter, so
  ! the default errs wide. Setting the angle bypasses source_filter_tolerance below.
  real(rp) :: source_filter_angle = 0
  ! The largest fraction of source intensity the filter may take anywhere inside the
  ! protected angle, which is the larger of the two derived angles. The edge is then
  ! placed outside that angle by however much the sigmoid's width needs to hold the
  ! transmission. A half-amplitude edge sitting on the protected angle passes a quarter
  ! of the intensity there, so 0.75 is the loss the placement used before this input
  ! existed implies, and it reproduces that edge to the bit. Strictly between 0 and 1.
  real(rp) :: source_filter_tolerance = 0.75_rp
  ! The sigmoid's edge in Genesis4's own units, xcut and ycut: the shifted grid index over
  ! ngrid, so 1 is twice the Nyquist index and the angle xcut lambda/dx. Validation-internal:
  ! they place the edge where Genesis4 places it, and a run that sets them and the angle
  ! together is refused. An angle is independent of the grid and these are not.
  real(rp) :: source_filter_xcut = 0
  real(rp) :: source_filter_ycut = 0
  ! The sigmoid's width as a fraction of its edge, Genesis4's sigmoid. Small is a sharp
  ! edge. Genesis4's own default of 1 leaves the sigmoid at 0.73 on axis, which attenuates
  ! the coherent source as much as the wide angles and measurably delays saturation.
  real(rp) :: source_filter_width = 0.05
  ! The filter check's self-test, the fp32_mutate pattern: apply the sigmoid to the
  ! propagated field instead of to the source. Both suppress wide angles, so a run still
  ! completes and still looks reasonable, and only the comparison against Genesis4 says
  ! which one the code did.
  logical :: source_filter_mutate = .false.
  ! The tracking window, Tao's names (tao_beam_init carries track_start/track_end) and
  ! Genesis's zstop parity: element locators (lat_ele_locator syntax). Blank = the whole
  ! line. The schedule (slippage, autophasing, break geometry) is always built on the
  ! Full lattice, so a windowed run composes exactly with the full run: [1,k] then
  ! [k+1,end] from its dumps reproduces [1,end] bit for bit.
  character(60) :: track_start = '', track_end = ''
  ! The comb (Bmad's bunch_track_struct%ds_save name and semantics, Tao's
  ! comb_ds_save): the minimum z advance between per-record stats rows.
  !   < 0: no per-record rows at all (element ends, dumps and the final state remain).
  !   = 0: a row at every record position (the default here: the per-record arrays are
  !        the primary stats contract, where Tao's comb is an optional extra. The one
  !        deliberate deviation from Tao's -1 default, with identical semantics).
  !   > 0: a row when z has advanced comb_ds_save past the last row, and element ends always.
  real(rp) :: comb_ds_save = 0
  character(60) :: dump_beam_at(40) = ''     ! Element locators for mid-run beam dumps.
  character(60) :: dump_field_at(40) = ''    ! Element locators for mid-run field dumps.
  ! Write the beam and the field at every comb position as a frame series, for a
  ! visualization that wants z resolution inside the undulators rather than at their
  ! ends. The comb already schedules the stats rows, so a frame and its row share an
  ! index (doc/reading-output.md).
  logical :: dump_at_comb = .false.
  ! The slices a frame carries, in window order. The default pair is the whole window,
  ! and a range cuts the series where a view wants only part of the bunch: a sample = 1
  ! run has thousands of slices and a whole-window field frame is gigabytes. The
  ! element-end, initial and final dumps are unaffected, being restart points.
  integer :: dump_slice_first = 1
  integer :: dump_slice_last = -1            ! -1 is the last slice of the window.
  ! The field's own reductions, per record, into the stats file: the projections a
  ! picture is drawn from, written once rather than recomputed from raw frames.
  logical :: dump_reduced = .false.
  integer :: ran_seed = 12345                ! The one RNG seed (generation, import, noise).
  logical :: write_diag = .false.            ! The Genesis-comparison text diag file (large).
  logical :: write_initial = .false.         ! Dump the initial state before tracking.
  logical :: load_only = .false.             ! Build the initial state, dump it, stop.
  logical :: keep_escaped_field = .false.    ! Keep the escaped-slice bank file.
  logical :: migrate = .false.               ! Slice migration (fel-physics.md sec-migration).
  logical :: migrate_check = .false.         ! Migration's bunching-invariance instrument.
  logical :: reference_run = .false.         ! No FEL interaction: Bmad tracks everything.
  ! Provenance detail in stats.h5's meta/ group. Off by default because a stats file is
  ! meant to travel: attached to a paper, mailed to a collaborator, posted beside a
  ! figure. On, meta/ also records the user name and the working directory, which is
  ! useful in a lab notebook and a leak outside one. Genesis records them always. Parity
  ! is not a reason to leak.
  logical :: record_environment = .false.
  ! The FP32 lockstep instrument (fel_fp32_mod): 'off', 'lockstep' (the FP32 twin
  ! rebuilt from FP64 every step) or 'freerun' (its longitudinal state carries, so the
  ! compounding rate is measured). A check instrument, not a production mode.
  character(16) :: fp32_check = 'off'
  logical :: fp32_mutate = .false.     ! The instrument's self-test: coarsen the residual.
  ! The device backend (fel_device_mod): 'off', or 'metal' for the Apple Silicon
  ! backend on a build that carries it. Alone it runs the averaged FEL step resident
  ! on the device; combined with fp32_check the device takes the lockstep twin's role
  ! and the FP64 run is untouched. Everything the backend does not cover is refused
  ! at setup or first use, never quietly run on the CPU instead.
  character(16) :: device = 'off'
end type

!+
! Struct wavefront_init_struct
!
! The radiation starting condition, the beam_init_struct analog (&fel_wavefront_init).
! A wavefront is an optical object: a transverse grid, a wavelength, a longitudinal
! spacing and nothing else longitudinal. The time window that sets that spacing is the
! FEL interaction's, not the field's, and lives in fel_slicing_struct. harmonics requests
! the field set. field_file imports override the seed.
!-

type wavefront_init_struct
  real(rp) :: lambda0 = 0            ! Resonant wavelength [m]. Required for generation.
  ! The transverse grid (Genesis ngrid and dgrid). Zero means derive it from the beam
  ! at setup, by fel_cells_per_sigma$ and fel_widths_per_sigma$ below. One may be set
  ! and the other derived against it.
  integer :: grid_n_pts = 0          ! Transverse grid points per side.
  real(rp) :: grid_half_width = 0    ! Transverse half width [m].
  real(rp) :: seed_power = 0         ! Gaussian seed power [W]. 0 = dark start.
  real(rp) :: seed_waist_size = 0    ! Seed intensity 1/e^2 radius [m].
  character(1) :: seed_polarization = 'x'   ! 'x' or 'y'.
  ! The field set: harmonics(1) must be 1 (the fundamental). Further entries are
  ! harmonic numbers in increasing order, 0 = unused.
  integer :: harmonics(9) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
end type

! The transverse grid derived from the beam when the deck states none
! (doc/startup-noise.md, Recommendations). Cells of a seventh of the rms beam size
! resolve the beam and the mode, and finer cells raise the wide-angle emission of the
! point beamlets without moving the mode power. A half width of nine beam sizes contains
! the mode, and the half width does not enter above that. Both were measured on three
! machines a hundred times apart in wavelength (FINDINGS 7.54).

real(rp), parameter :: fel_cells_per_sigma$ = 7
real(rp), parameter :: fel_widths_per_sigma$ = 9

!+
! Struct fel_chamber_wake_init_struct
!
! The chamber-wake description (&fel_params chamber_wake%...), Genesis &wake names. A plain
! mirror of fel_wake_struct's configuration fields (that struct carries allocatable
! state, which a namelist object cannot). fel_setup copies these in.
!-

type fel_chamber_wake_init_struct
  logical :: on = .false.
  ! Which chamber-wake implementation runs. 'genesis' is the transcribed solver at
  ! slice granularity, convolving the kernels with the weighted slice currents. The
  ! field exists so a second implementation arrives as a value, not a rework.
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

! (Space charge needs no init mirror: fel_space_charge_struct is already pure scalars and
! reads directly in &fel_params as space_charge%... .)

!+
! Struct fel_beam_init_param_struct
!
! The beam-side scalars of &fel_beam_init that sit beside Bmad's beam_init (which
! describes the bunch itself) and imp (the resampler): source/output files, the
! quiet-start knobs, and the check instruments the validation harness sets.
!-

type fel_beam_init_param_struct
  character(400) :: beam_file = ''        ! An FEL beam already in slices (with field_file(1)).
  character(400) :: write_genesis_dist = ''  ! Write the bunch as a Genesis4 distribution file.
  character(400) :: write_openpmd_file = ''  ! Write the bunch as openPMD-beamphysics.
  ! How a bunch becomes slices (fel-physics.md sec-loading). The bunch is beam_init's,
  ! generated or read through beam_init%position_file, and is binned by arrival time into
  ! slicing's slices. "sample" draws the same number of beamlets in every slice from what
  ! fell into it, weights from the slice charge, so the current follows the bunch: Genesis's
  ! &beam with &profile and its importdistribution in one, with beam_init%n_particle the
  ! count per slice after phase copies. "keep" keeps every real particle as a beamlet of
  ! beamlet_size copies at weight/beamlet_size, counts following the charge, so a
  ! start-to-end bunch's correlations survive whole.
  character(16) :: load_mode = 'sample'
  ! The two switches of the load. quiet_start makes it quiet: beamlets of copies with phases
  ! spread over 2 pi about each particle's own, so the load has no bunching at any harmonic
  ! below beamlet_size. shot_noise then imposes the physical level on whatever load there
  ! is, by independent phasors per group. Off both, the beam is taken as it came.
  logical :: quiet_start = .true.
  integer :: beamlet_size = 8             ! Copies per beamlet. A loading parameter only.
  logical :: shot_noise = .false.         ! Impose physical shot noise on the load.
  ! Check instruments (the validation harness's knobs, not physics inputs):
  logical :: split_weights = .false.      ! Coincident w/3 + 2w/3 copies after loading.
  logical :: swap_beam_xy = .false.       ! Swap (x,px) <-> (y,py) after generation.
  logical :: gen_test_weights = .false.   ! Alternate beamlet weights 0.25x/1.75x.
  logical :: resample_split_weights = .false.  ! Split-weight copies before the resample.
end type

!+
! Struct fel_run_struct
!
! The assembled state of one run: the parsed lattice, the beam, the field set, the
! collective-effects state, the schedule the setup pass computed, the diagnostics
! state, and copies of every input struct (so the resolved inputs are one object:
! the Meta/ provenance echo and any embedding caller read them from here). Passed
! explicitly, no singleton.
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
  ! The angle at which the stats split the field power, whether or not the filter is on:
  ! the larger of four mode diffraction angles and the rho angle, derived from the beam
  ! and the first FEL element (fel-physics.md sec-source-filter). Inside it lies what can
  ! couple to the mode, outside it the wide-angle emission of the point beamlets.
  real(rp) :: split_angle = 0
  character(24) :: split_origin = ''       ! Where split_angle came from, for the report.
  integer, allocatable :: fel_mode(:), fel_spp(:)
  real(rp), allocatable :: fel_ramp(:)
  real(rp), allocatable :: ele_slip(:)     ! Slippage after each element's last step [wavelengths].
  real(rp), allocatable :: fel_zoff(:)     ! The z_offset off-phase knob per element [m].
  real(rp), allocatable :: light_corr(:)   ! Chord-vs-arc correction on a break's last element [m].
  logical, allocatable :: is_fel(:)
  ! Space charge, per element: ele%space_charge_method = slice, and Bmad's master switch
  ! bmad_com%csr_and_space_charge_on, resolved once at setup so the walk reads a fact.
  logical, allocatable :: sc_here(:)
  logical, allocatable :: dump_beam_here(:), dump_field_here(:)
  integer :: i_start = 1, i_end = 0        ! The resolved tracking window [elements].
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
! Routine to format a value for a human: an SI prefix chosen so the mantissa lands in
! [1, 1000), three decimals, and the whole thing right-justified to a fixed width so a
! column of them lines up (4.230 kW and 105.000 GW under each other). This is display
! only -- stdout is for humans and the files carry full precision (doc/user-guide.md).
!
! Values outside the prefix range, and exact zero, fall back to es10.3 with the bare
! unit rather than inventing a prefix.
!
! Input:
!   value -- real(rp): The value, in the unit's own base (W, J, m, ...).
!   unit  -- character(*): The unit symbol, appended after the prefix.
!
! Output:
!   str   -- character(14): The formatted value, right-justified so a column of them
!              lines up. Callers wanting it inline use trim(adjustl(...)).
!-

function fel_si_str (value, unit) result (str)

real(rp) value, av, mant
character(*) unit
character(14) str
character(14) num
integer ip

! The table spans yocto to peta, which covers everything this program reports: a wake
! run's pulse energy is attojoules and a saturated pulse is millijoules. Index 9 is the
! bare unit. Outside the table, and for exact zero, no prefix is invented.

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
! The comb rule, in one place: Bmad's bunch_track_struct%ds_save semantics, with one
! deliberate difference (save_a_bunch_step's guards, and Tao's comb_ds_save note
! "< 0 => No comb calculated"):
!   comb < 0: element ends only;
!   comb = 0: a row at every record position;
!   comb > 0: a row when z has advanced comb_ds_save past the last row.
! An element end is always a row, whatever the comb. That is the difference, and it is
! what lets the stats file carry one record axis with a boolean mask (coords/
! at_element_end) instead of a second axis and a duplicated copy of every element-end
! quantity. Bmad's comb < 0 drops the comb, and here that leaves the element ends,
! which are the positions the evaluated bunch_params live on.
! z_last updates when the row is taken. The walk consults this rule live and the
! setup's nrec precompute replays it with the same z arithmetic, so the stats
! arrays are exact-sized in every mode.
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

end module fel_struct
