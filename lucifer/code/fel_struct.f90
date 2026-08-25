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
!     wake%radius = 2.5e-3
!     sc%nz = 2
!   /
!   &fel_beam_init
!     beam_init%n_particle = 8192
!     shotnoise = T
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
! Struct fel_global_struct
!
! The run-level switches and names, exposed in &fel_params as global%... (the
! tao_global_struct analog). Everything here is about ONE run of the tracker. The
! physics description lives in the lattice, beam_init and wavefront_init.
!-

type fel_global_struct
  character(400) :: out_root = 'fel_track'   ! Output file root.
  character(16) :: interlude_model = 'bmad'  ! 'bmad' (the seam) or 'genesis' (transcribed).
  ! The FEL source model (manual sec:coherent-source). 'deposit' is the standard
  ! per-particle scatter (the referee, bit-for-bit unchanged). 'coherent' is the
  ! SIMPLEX-hybrid coherent-Gaussian source (Tanaka, PRAB 27, 030703 (2024)): the
  ! spatially incoherent artifact is dropped, the slice bunch factor B(s) keeps the
  ! physical shot noise, and the transverse shape is a guarded Gaussian.
  character(16) :: source_model = 'deposit'
  ! The tracking window, Tao's names (tao_beam_init carries track_start/track_end) and
  ! Genesis's zstop parity: element locators (lat_ele_locator syntax). Blank = the whole
  ! line. The schedule (slippage, autophasing, break geometry) is always built on the
  ! FULL lattice, so a windowed run composes exactly with the full run: [1,k] then
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
  integer :: ran_seed = 12345                ! The one RNG seed (generation, import, noise).
  logical :: write_diag = .false.            ! The Genesis-comparison text diag file (large).
  logical :: write_initial = .false.         ! Dump the initial state before tracking.
  logical :: load_only = .false.             ! Build the initial state, dump it, stop.
  logical :: keep_escaped_field = .false.    ! Keep the escaped-slice bank file.
  logical :: migrate = .false.               ! Slice migration (manual sec:migration).
  logical :: migrate_check = .false.         ! Migration's bunching-invariance instrument.
  logical :: reference_run = .false.         ! No FEL interaction: Bmad tracks everything.
end type

!+
! Struct wavefront_init_struct
!
! The radiation starting condition, the beam_init_struct analog (&fel_wavefront_init).
! The field record IS the time window, so the window lives here: window_length and
! window_sample set the slice count and spacing for the field AND the generated beam
! (one window, one definition). harmonics requests the field set. field_file imports
! override the seed.
!-

type wavefront_init_struct
  real(rp) :: lambda0 = 0            ! Resonant wavelength [m]. Required for generation.
  real(rp) :: window_length = 0      ! Time window [m]. 0 = derived from the bunch.
  integer :: window_sample = 1       ! Slice spacing in wavelengths (Genesis's sample).
  integer :: grid_n_pts = 255        ! Transverse grid points per side (Genesis ngrid).
  real(rp) :: grid_half_width = 0    ! Transverse half width [m] (Genesis dgrid).
  real(rp) :: seed_power = 0         ! Gaussian seed power [W]. 0 = dark start.
  real(rp) :: seed_waist_size = 0    ! Seed intensity 1/e^2 radius [m].
  character(1) :: seed_polarization = 'x'   ! 'x' or 'y'.
  ! The field set: harmonics(1) must be 1 (the fundamental). Further entries are
  ! harmonic numbers in increasing order, 0 = unused.
  integer :: harmonics(9) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
end type

!+
! Struct fel_wake_init_struct
!
! The chamber-wake description (&fel_params wake%...), Genesis &wake names. A plain
! mirror of fel_wake_struct's configuration fields (that struct carries allocatable
! state, which a namelist object cannot). fel_setup copies these in.
!-

type fel_wake_init_struct
  logical :: on = .false.
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

! (Space charge needs no init mirror: fel_efield_struct is already pure scalars and
! reads directly in &fel_params as sc%... .)

!+
! Struct fel_beam_init_param_struct
!
! The beam-side scalars of &fel_beam_init that sit beside Bmad's beam_init (which
! describes the bunch itself) and imp (the resampler): source/output files, the
! quiet-start knobs, and the check instruments the validation harness sets.
!-

type fel_beam_init_param_struct
  character(400) :: beam_file = ''        ! Genesis .par.h5 dump (with field_file(1)).
  character(400) :: dist_file = ''        ! openPMD-beamphysics particle file to import.
  character(400) :: write_dist_file = ''  ! Write the bunch as a Genesis DISTRIBUTION file.
  character(400) :: write_opmd_file = ''  ! Write the bunch as openPMD-beamphysics.
  logical :: use_beam_init = .false.      ! Generate the bunch from beam_init, then import.
  integer :: nbins = 8                    ! Beamlet size of the quiet start.
  logical :: shotnoise = .false.          ! Physical (Fawley) shot noise on the phases.
  ! Check instruments (the validation harness's knobs, not physics inputs):
  logical :: split_weights = .false.      ! Coincident w/3 + 2w/3 copies after loading.
  logical :: swap_beam_xy = .false.       ! Swap (x,px) <-> (y,py) after generation.
  logical :: gen_test_weights = .false.   ! Alternate beamlet weights 0.25x/1.75x.
  logical :: imp_split_weights = .false.  ! Split-weight copies BEFORE the import resample.
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
  type (fel_wake_init_struct) :: wake_init
  type (fel_efield_struct) :: sc_init
  type (beam_init_struct) :: beam_init
  type (fel_import_param_struct) :: imp
  type (fel_beam_init_param_struct) :: bparam
  character(400) :: lat_file = ''
  character(400) :: field_file(9) = ''
  ! The assembled state.
  type (lat_struct) :: lat
  type (fel_beam_struct) :: fbeam
  type (fel_field_struct), allocatable :: ffield(:)
  type (fel_collective_struct) :: coll
  type (fel_stats_struct) :: stats
  type (fel_unavg_struct) :: ustate
  ! The schedule (fel_setup): per tracked element.
  type (fel_und_struct), allocatable :: und_of(:)
  integer, allocatable :: fel_mode(:), fel_spp(:)
  real(rp), allocatable :: fel_ramp(:)
  real(rp), allocatable :: ele_slip(:)     ! Slippage after each element's last step [wavelengths].
  real(rp), allocatable :: fel_zoff(:)     ! The z_offset off-phase knob per element [m].
  real(rp), allocatable :: light_corr(:)   ! Chord-vs-arc correction on a break's last element [m].
  logical, allocatable :: is_fel(:)
  logical, allocatable :: dump_beam_here(:), dump_field_here(:)
  integer :: i_start = 1, i_end = 0        ! The resolved tracking window [elements].
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
! only -- stdout is for humans and the files carry full precision (manual sec:program).
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
! THE COMB RULE, in one place: Bmad's bunch_track_struct%ds_save semantics verbatim
! (save_a_bunch_step's guards, and Tao's comb_ds_save note "< 0 => No comb calculated"):
!   comb < 0: no per-record rows at all (element ends, dumps and finals remain
!             through their own arrays);
!   comb = 0: a row at every record position;
!   comb > 0: a row when z has advanced comb_ds_save past the last row, and always
!             at an element end.
! z_last updates when the row is taken. The walk consults this rule live and the
! setup's nrec precompute REPLAYS it with the same z arithmetic, so the stats
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

if (comb_ds_save < 0) then
  take = .false.
elseif (at_end) then
  take = .true.
else
  take = (z >= z_last + comb_ds_save)
endif
if (take) z_last = z

end function fel_comb_take

end module fel_struct
