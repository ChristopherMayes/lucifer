!+
! Program fel_track_test
!
! FEL tracker validated against Genesis 1.3 Version 4. The PHYSICS -- coordinates and
! conventions, the FEL step, the field solver, slippage, loading, import, migration,
! collective effects, and each piece's Genesis provenance and measured validation
! level -- lives in the manual, bsim/fel/doc/fel-physics.tex. Measured numbers and how
! to run the checks: bsim/fel/README.md. This header documents the inputs.
!
! The program walks a Bmad lattice and applies the seam of the design (manual
! sec:element, sec:seam): FEL segments are real Bmad wiggler/undulator elements with
! tracking_method = custom, their FEL parameters read from lattice attributes (aw from
! b_max/l_period, helicity from field_calc, step from ds_step -- assertions enforced by
! name at setup) and stepped with the transcribed Genesis physics of fel_track_mod;
! every other element tracks each slice's bunch with Bmad's track1_bunch (the packed
! arrays ARE Bmad coordinates, so conversion is a plain copy) and drifts the field by
! wavefront_drift, with one advance of the reference phase phi0 per element.
!
! Time dependence follows from the starting state alone: a multi-slice window makes a
! time-dependent run (slippage active), a single slice the steady state, no separate
! switch -- the same rule as Genesis. The slippage schedule is precomputed over the
! lattice, transcribing Lattice::calcSlippage (manual sec:slippage, including the
! drift autophasing and its unguarded end-of-lattice fixup), and applied after each
! step's field solve, before its diagnostics -- Gencore's step order. The field record
! rotates rather than moves; everything reading it in time order goes through
! fel_field_index.
!
! The starting state is a pair of Genesis dumps (&write of beam and field), so both
! codes track from bitwise-identical initial conditions; or self-generated; or an
! imported distribution (below). Diagnostics matching Genesis's definitions are
! recorded at the same z positions Genesis records them: once at the start and once
! after every integration step, one step per interlude element -- one row per slice
! per record, in time-window order.
!
! Input is a namelist file:
!
!   &fel_track_params
!     lat_file = "aramis.bmad"                 ! Bmad lattice.
!     beam_file = "Aramis-initial.par.h5"      ! Genesis particle dump to start from.
!     field_file = "Aramis-initial.fld.h5"     ! Genesis field dump to start from.
!     out_root = "fel_td"                      ! Prefix for the output files.
!     interlude_model = "bmad"                 ! "bmad" (the seam, default) or "genesis".
!     split_weights = F                        ! Weight-invariance test mode; see below.
!
! The FEL tracking mode and unaveraged parameters are per-element LATTICE attributes
! (registered by this program; usable on any wiggler/undulator, class-settable as
! wiggler::*[attr] = ...). Use NAMED values, defined as one-line lattice variables --
! fel_transcribed = -1, fel_averaged = 0, fel_unaveraged = 1 -- matching the code's
! fel_transcribed$/fel_averaged$/fel_unaveraged$ parameters (fel_track_mod):
!
!     fel_tracking          ! unset/0 = averaged, the bmad_standard wiggler kernel's
!                           !   transverse maps -- BMAD'S OWN KERNEL IS THE DEFAULT.
!                           ! 1 = the unaveraged verification mode (full Newton-Lorentz
!                           !   quiver, no fc/faw; manual sec:unaveraged; the run
!                           !   writes <out_root>.ledger.txt).
!                           ! -1 = averaged with the transcribed-Genesis transverse
!                           !   maps: VALIDATION-INTERNAL (the Genesis tiers require
!                           !   transcription-level transport; no production lattice
!                           !   writes it). Differences priced in the README.
!     fel_steps_per_period  ! Unaveraged substeps per period. unset/0 -> 20;
!                           !   below 10 refused (MINERVA's floor).
!     fel_ramp_periods      ! sin^2 entry/exit ramp length [periods]. unset/0 -> 2;
!                           !   a TRUE hard edge (test configuration) is the explicit
!                           !   sentinel -1 -- silence never means hard edge.
!     write_initial = F                        ! Also dump the initial state (Genesis format).
!     dump_beam_at = "UND3", "quadrupole::*"   ! Dump the beam (Genesis .par.h5 format) at the
!                                              !   end of the named elements (Bmad locator
!                                              !   syntax, class::name allowed; an entry
!                                              !   matching nothing is refused by name).
!     dump_field_at = "UND3"                   ! Same for the field (Genesis .fld.h5, unrotated).
!     keep_escaped_field = F                   ! Bank the field slices slippage transmits out of
!                                              !   the window (<out_root>-escaped.fld.h5, with
!                                              !   wavefront_params and z_transmit per slice) and
!                                              !   reconstruct the FULL PULSE at the exit plane
!                                              !   (<out_root>-pulse.fld.h5) by free-space
!                                              !   propagation at finalize. Manual sec:stats.
!     migrate = F                              ! Slice migration (see below).
!     migrate_check = F                        ! Verify phase continuity at each migration.
!   &end
!
! migrate = T moves particles between slices when their ponderomotive phase leaves the
! slice window (fel_migrate_slices, manual sec:migration), called serially after every
! element. OFF BY DEFAULT, deliberately: the Genesis-comparison tiers run against
! Genesis WITHOUT one4one, which never migrates, so migration would be a physics-model
! difference inside a transcription-level comparison. Dropped charge is counted and
! reported. migrate_check = T additionally verifies exact phase continuity at every
! migration and reports the worst deviation.
!
! Alternatively, leave beam_file and field_file blank and the program generates its own
! starting condition -- a quiet-start beam and a Gaussian seed field -- making a
! self-contained single run with no Genesis anywhere (see bsim/fel/examples). The BEAM
! is described by Bmad's standard beam_init_struct (the same block the import path
! uses: one bulk-bunch description, two generation methods -- the import resamples
! real particles from it, the quiet-start loader evaluates it analytically per slice).
! The Twiss is always the lattice's. Honored beam_init fields:
!
!     beam_init%n_particle    ! Macroparticles PER SLICE; a positive multiple of nbins.
!                             !   (The import path reads it as bunch particles instead.)
!     beam_init%a_norm_emit   ! Normalized emittances [m rad] (a_emit/b_emit refused:
!     beam_init%b_norm_emit   !   normalized only, Bmad's preferred form).
!     beam_init%sig_pz        ! Fractional momentum spread dP/P0. The Gaussian gamma
!                             !   spread is delta_gamma = beta0*p0_mc*sig_pz. (sig_e is
!                             !   deprecated Bmad-wide and does not exist here.)
!     beam_init%bunch_charge  ! Charge [C]. The current is DERIVED, never input.
!     beam_init%sig_z         ! Bunch length [m], with distribution_type(3):
!     beam_init%distribution_type(3)  ! "RAN_GAUSS" (default): Gaussian current profile
!                             !   I(s) = Q*c/(sqrt(2pi)*sig_z) * exp(-s^2/2 sig_z^2) at
!                             !   the slice centers, bunch centered in the window.
!                             !   sig_z = 0 is the STEADY STATE: one slice holding the
!                             !   whole charge, I = Q*c/slice_spacing; refused by name
!                             !   for time-dependent windows.
!                             ! "GRID" (Bmad's uniform): flat I = Q*c/(x_max - x_min)
!                             !   over grid(3)%x_min..x_max (the z extent).
!
! EVERY other beam_init field that is set is REFUSED BY NAME (see
! check_beam_init_contract) -- a standard structure that silently drops fields would
! be worse than a custom one. Remaining generation knobs, with the Genesis4
! &field / &time names they map to:
!
!     lambda0 = 1e-10          ! Radiation wavelength [m], REQUIRED (with dumps it comes
!                              !   from the file). Deliberately not defaulted from the
!                              !   lattice resonance: the first undulator may be off.
!     nbins = 8                ! Beamlet size of the quiet start (quiet below nbins).
!     seed_power = 5e3         ! Seed power [W] (Genesis &field power); 0 = dark start.
!     seed_waist_size = 30e-6  ! Seed 1/e^2 intensity radius w0 [m] (&field waist_size).
!     grid_n_pts = 255         ! Transverse grid points per side (&field ngrid).
!     grid_half_width = 2e-4   ! Grid half width [m] (&field dgrid).
!     ran_seed = 12345         ! Random seed (Bmad's ran_seed_put), for reproducibility.
!     window_length = 0        ! Time window [m] (&time slen). 0: derived from the bunch
!                              !   (+-4 sig_z Gaussian; the grid extent flat; one slice
!                              !   steady state). Override it for slippage headroom; a
!                              !   window that clips the bunch warns with numbers.
!     window_sample = 1        ! Slice spacing / lambda0 (integer, &time sample).
!     shotnoise = F            ! Impose physical shot noise (&time shotnoise;
!                              !   time-dependent windows only, Genesis's dotime rule).
!     gen_test_weights = F     ! Validation knob: alternate beamlet weights 0.25x/1.75x
!                              !   (charge preserving, uniform within each beamlet) to
!                              !   exercise the weighted-noise paths. Not physics input.
!     load_only = F            ! Generate, write <out_root>-initial dumps, exit without
!                              !   tracking. For the shot-noise statistical check.
!
! Element wakes: elements carrying Bmad sr_wake definitions -- pseudomodes or a
! tabular z_long -- act across the WHOLE time window via slice concatenation (manual
! sec:seamwake: conventions, ds_wake, the mid-element wiggler kick, refusals).
! write_wake_kernels = "<file>" exports the transcribed wake kernels for building
! matching z_long tables (see the README's seam-wake section and examples/bmad_wake).
!
! Third way in: import a particle DISTRIBUTION -- an arbitrary bunch, resampled into
! slices by the transcribed Genesis importdistribution method (fel_import_mod; manual
! sec:import). The bunch comes from the SAME beam_init block (use_beam_init = T:
! init_beam_distribution generates it, honoring everything Bmad honors -- the
! honored-fields contract above applies to the quiet-start generator only) or from an
! openPMD-beamphysics file (dist_file). window_sample, ran_seed and the seed field are
! shared with the generator; one seed governs generation, resampling and noise. Knobs,
! named after &importdistribution's where one exists:
!
!     imp%slicewidth = 0.01    ! Sampling window / bunch length (Genesis's slicewidth).
!     imp%npart = 8192         ! Macroparticles per slice after resampling.
!     imp%nbins = 4            ! Beamlet size (quiet-load bins) of the resample.
!     imp%nslice = 0           ! 0: round(bunch_length/spacing), Genesis's rule.
!                              ! (Genesis's match/center are NOT ported: a Bmad lattice
!                              !   carries its Twiss and beam_init generates matched
!                              !   bunches already -- see fel_import_mod's header.)
!     use_beam_init = F        ! Generate the bunch from the beam_init%... block.
!     dist_file = ""           ! Or read an openPMD-beamphysics file.
!     write_dist_file = ""     ! Write the bunch as a Genesis &importdistribution
!                              !   input (t/p/x/xp/y/yp + charge, t = -tau/c) -- the
!                              !   shared file of the cross-code checks.
!     write_opmd_file = ""     ! Write the bunch as openPMD-beamphysics.
!     imp_split_weights = F    ! Check knob: coincident w/3 + 2w/3 copies before import.
!
! The quiet start and the weighted Fawley shot noise (shotnoise = T), including
! the noise-level algebra and the N_eff refusal guard, are the manual's sec:loading;
! the loader warns where Genesis silently clamps beamlets with fewer real electrons
! than macroparticles.
!
! interlude_model selects how the field-free elements are handled: "bmad" is the seam
! (track1_bunch, exact theta mapping, wavefront_drift); "genesis" uses the transcribed
! Genesis interlude step everywhere, which prices what the seam changes -- the manual's
! sec:interlude and sec:seam. The slippage schedule is identical in both models.
!
! split_weights = T replaces each imported particle by two copies at identical
! coordinates carrying 1/3 and 2/3 of its weight; every collective observable must be
! identical to the unsplit run. This checks the weighted paths, which no Genesis
! comparison can (Genesis dumps carry no weights).
!
! Outputs: <out_root>.diag.txt (one row per slice per record: z, slice, field and beam
! diagnostics), <out_root>.stats.h5 (the production statistics file, manual sec:stats:
! per-record per-slice beam moments named as bunch_params_struct components, per-record
! per-slice wavefront_params, and the evaluated calc_bunch_params at element ends;
! fixed Bmad units), <out_root>-final.fld.h5 and <out_root>-final.par.h5 (Genesis-format dumps
! of the end state, for field-by-field comparison; the field dump is unrotated to time
! order first, as writeFieldHDF5 does).
!-

program fel_track_test

use fel_track_mod
use fel_unaveraged_mod
use fel_import_mod
use wavefront_hdf5_mod
use fel_stats_mod
use wavefront_openpmd_mod
use beam_mod
use wake_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: wake_src   ! The element whose wake applies at ele, resolved
                                         ! through lords (pointer_to_wake_ele) -- a wake on
                                         ! a superimposed/split element lives on the LORD
                                         ! and ele%wake is null on its slaves.
type (fel_beam_struct), target :: fbeam
! The radiation field SET (manual sec:field-set): ffield(1) is the fundamental,
! further entries are harmonic fields (namelist harmonics). The wf/slip/bank pointers
! alias entry 1, so every fundamental-only statement below reads exactly as it did
! when the walk carried one field -- the bit-for-bit keystone.
type (fel_field_struct), allocatable, target :: ffield(:)
type (wavefront_struct), pointer :: wf => null()
type (fel_slip_struct), pointer :: slip => null()
type (fel_bank_struct), pointer :: bank => null()
type (bunch_struct) bunch
type (bunch_struct) wake_bunch          ! Whole-window bunch for element sr wakes.
real(rp), allocatable :: wake_beta0(:)  ! Entry betas for the exact concat/split inverse.
type (fel_und_struct) und
type (fel_slice_diag_struct) bdiag

real(rp) :: gamma0 = 0                  ! DERIVED from the lattice e_tot.
logical :: split_weights = .false.
logical :: write_initial = .false.
logical :: migrate = .false., migrate_check = .false.
! Bmad's own synchrotron radiation, for elements Bmad tracks (interludes and any
! wiggler whose tracking_method is not custom). The FEL path has its own radiation
! physics -- these knobs do not touch it. Used by check_spontaneous.py to measure
! Bmad's independent spontaneous-loss implementation against the FEL modes and the
! analytic rate; off by default, as in Bmad.
logical :: radiation_damping = .false., radiation_fluctuations = .false.
! A deliberate run with NO FEL interaction: Bmad tracks every element, the field just
! drifts. The reference leg of check_spontaneous.py, which measures Bmad's own
! radiation physics through the same wiggler the FEL modes use.
logical :: reference_run = .false.
! Two-polarization radiation (manual sec:field vector convention): Ey goes LIVE when
! any FEL element is tilted or the seed is y-polarized; single-polarization lines
! never allocate it and stay bit-identical.
character(1) :: seed_polarization = 'x'
logical :: two_pol = .false.
! Check instrument (check_two_polarization.py's rotation identity): swap the beam's
! transverse planes after generation, (x,px) <-> (y,py). A y-planar line fed the
! swapped beam must reproduce the x-planar line fed the original, exactly -- the RNG
! draws its planes sequentially, so the generated beam itself is never swap-symmetric.
logical :: swap_beam_xy = .false.
! field_file: entry 1 starts the FUNDAMENTAL (Genesis dump or openPMD wavefront,
! auto-detected by signature); entries 2+ are per-harmonic imports, each matched to
! the field-set entry whose photon energy the file carries -- no match is refused.
character(400) :: lat_file = '', beam_file = '', out_root = 'fel_track'
character(400) :: field_file(9) = ''
character(16) :: interlude_model = 'bmad'

! Collective effects (deliverable 8), Genesis &wake / &efield names with wake_/sc_
! prefixes. All off by default; see the header.
logical :: wake_on = .false.
real(rp) :: wake_loss = 0, wake_radius = 2.5e-3_rp, wake_conductivity = 0, wake_relaxation = 0
real(rp) :: wake_gap = 0, wake_lgap = 1, wake_hrough = 0, wake_lrough = 1
logical :: wake_roundpipe = .true.
character(8) :: wake_material = ''
real(rp) :: sc_rmax = 0
integer :: sc_ngrid = 100, sc_nz = 0, sc_nphi = 0
logical :: sc_longrange = .false.

! The FEL tracking mode and unaveraged parameters are per-element LATTICE attributes
! (fel_tracking / fel_steps_per_period / fel_ramp_periods; see the header), read at
! setup into these arrays. any_unavg is derived: does ANY element run unaveraged?
type (fel_unavg_struct) ustate
integer, allocatable :: fel_mode(:), fel_spp(:)
real(rp), allocatable :: fel_ramp(:)
real(rp) :: dE_step, dU_step, u_spont_cum = 0
! Spontaneous radiation inside FEL elements (bmad_com's GLOBAL switches; manual
! sec:core spontaneous paragraph): the actual drawn energy the beam radiated away,
! the ledger's E_radiated column, plus the per-record scratch.
real(rp) :: e_rad_cum = 0
real(rp), allocatable :: e_rad_slice(:), rad_kick(:,:)
integer :: iu_ledger = 0
logical any_unavg

real(rp) :: lambda0 = 0                  ! Generation parameters; see the header.
! The generated beam is described by beam_init (honored fields in the header table);
! these are the remaining generation knobs, Genesis4 mappings in the header.
real(rp) tw_beta_x, tw_alpha_x, tw_beta_y, tw_alpha_y   ! From the lattice beginning
                                                        ! element -- NOT namelist input:
                                                        ! the lattice is the one Twiss
                                                        ! authority (as for the import).
real(rp) :: seed_power = 0, seed_waist_size = 0, grid_half_width = 0
real(rp) :: window_length = 0
integer :: nbins = 8, grid_n_pts = 255, ran_seed = 12345
integer :: window_sample = 1
logical :: shotnoise = .false., gen_test_weights = .false., load_only = .false.
integer npart_gen                        ! From beam_init%n_particle (per slice here).
real(rp) delgam_gen                      ! beta0*p0_mc*beam_init%sig_pz.

! Distribution import (deliverable 10): a bunch_struct -- generated natively from
! Bmad's beam_init_struct, or read from an openPMD-beamphysics file -- resampled into
! FEL slices by the transcribed Genesis method (fel_import_mod). The knobs mirror
! &importdistribution's names in the imp block; window_sample/ran_seed and the seed
! field are shared with the quiet-start generator (one field generator, one seed).
type (beam_init_struct) :: beam_init     ! Bmad's native bunch description (&beam_init).
type (fel_import_param_struct) :: imp
logical :: use_beam_init = .false.       ! Generate the bunch from beam_init.
character(400) :: dist_file = ''         ! Or read it from an openPMD-beamphysics file.
character(400) :: write_dist_file = ''   ! Write the bunch as a Genesis DISTRIBUTION
                                         ! file (t/p/x/xp/y/yp + charge, t = -tau/c),
                                         ! the shared input of the cross-code checks.
character(400) :: write_opmd_file = ''   ! Write the bunch as openPMD-beamphysics
                                         ! (hdf5_write_beam), the dist_file round trip.
logical :: imp_split_weights = .false.   ! Check knob: coincident w/3 + 2w/3 copies
                                         ! BEFORE import; RNG-free outputs must not move.
character(400) :: write_wake_kernels = ''  ! Write the deliverable-8 wake kernels
                                           ! (s, wakeres, wakegeo, wakerou; eV/(m e-))
                                           ! for the seam-wake cross-validation check.

real(rp), allocatable :: ele_slip(:)     ! Slippage applied after each element's last step [wavelengths].
type (fel_und_struct), allocatable :: und_of(:)   ! Per-element FEL parameters, from lattice attributes.
logical, allocatable :: is_fel(:)                 ! Which tracked elements are FEL segments.
real(rp) z_now, ks, qf, und_slip_step, Lz, gamma0_ref
real(rp) charge_dropped_tot, b_dev_max
integer ie, is, ih, istep, n_arg, iu_diag, iu_nml, nslice, prev_ie, n_moved_tot, iu_wake
logical err, timerun

type (fel_collective_struct) coll

! Diagnostics (manual sec:stats): the stats file, dump-at element lists, and the
! escaped-field bank (transmitted slices streamed to file; the full pulse is
! reconstructed at finalize by free-space propagation -- transmitted light never
! re-interacts, so it is fixed information).
type (fel_stats_struct) stats
type (fel_slice_diag_struct), allocatable :: bdiag_arr(:)  ! Diag rows, filled by the stats loop.
real(rp), allocatable :: fpow_arr(:), fonax_arr(:)         ! fel_field_diag per slice, same loop.
type (wavefront_struct) wf1     ! One-slice scratch for banked-slice propagation.
character(60) :: dump_beam_at(40) = '', dump_field_at(40) = ''
logical :: keep_escaped_field = .false.
logical, allocatable :: dump_beam_here(:), dump_field_here(:)
integer nrec_stats, nend_stats
integer(8) prog_count0, prog_count_last, prog_rate    ! Wall clock for progress lines.
real(rp) lat_length
! Escaped-field stream state, one per field object (the fundamental's file keeps its
! pre-harmonic name; a harmonic's carries -h<h>).
integer, allocatable :: n_banked(:)
real(rp), allocatable :: bank_z(:,:), bank_pms(:,:,:)  ! (slot, field) / (25, slot, field).
integer(hid_t), allocatable :: esc_id(:)

! The field set requested by the namelist: harmonics(1) must be 1 (the fundamental);
! further entries are harmonic numbers in increasing order, 0 = unused. Harmonic
! fields share the fundamental's grid and start dark, growing from the bunching the
! fundamental drives (or filled by an openPMD import whose photonEnergy matches).
integer :: harmonics(9) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
integer :: n_harm = 1
! Radiation dump format: 'genesis' (the default; one polarization per file, the
! validation tiers' format), 'openpmd' (EXT_Wavefront; both polarizations as
! components of one mesh record, one file per harmonic), or 'both'.
character(8) :: wavefront_format = 'genesis'

character(400) param_file
character(*), parameter :: r_name = 'fel_track_test'

namelist / fel_track_params / lat_file, beam_file, field_file, out_root, &
                           interlude_model, &
                           split_weights, write_initial, lambda0, nbins, &
                           seed_power, seed_waist_size, grid_n_pts, &
                           grid_half_width, ran_seed, window_length, window_sample, shotnoise, &
                           gen_test_weights, load_only, migrate, migrate_check, &
                           wake_on, wake_loss, wake_radius, wake_conductivity, wake_relaxation, &
                           wake_roundpipe, wake_material, wake_gap, wake_lgap, wake_hrough, &
                           wake_lrough, sc_rmax, sc_ngrid, sc_nz, sc_nphi, sc_longrange, &
                           beam_init, imp, use_beam_init, dist_file, write_dist_file, &
                           write_opmd_file, imp_split_weights, write_wake_kernels, &
                           dump_beam_at, dump_field_at, keep_escaped_field, &
                           radiation_damping, radiation_fluctuations, reference_run, &
                           seed_polarization, swap_beam_xy, harmonics, wavefront_format

! Read parameters.

n_arg = command_argument_count()
if (n_arg /= 1) then
  print '(a)', 'Usage: fel_track_test <param_file>'
  stop 1
endif
call get_command_argument (1, param_file)

open (newunit = iu_nml, file = param_file, status = 'old', action = 'read')
read (iu_nml, nml = fel_track_params)
close (iu_nml)

! (The FEL tracking mode is a lattice attribute; its refusals live in
! setup_fel_elements, and the unaveraged-vs-collective refusal follows setup, once
! any_unavg is known.)

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  print '(a)', 'fel_track_test: interlude_model must be "bmad" or "genesis", got: ' // trim(interlude_model)
  stop 1
endif

bmad_com%radiation_damping_on = radiation_damping
bmad_com%radiation_fluctuations_on = radiation_fluctuations

if (radiation_fluctuations .and. migrate) then
  print '(a)', 'fel_track_test: radiation_fluctuations draws one kick per BEAMLET (the quiet start'
  print '(a)', '  cancels per beamlet), and slice migration scrambles beamlet grouping. Pick one.'
  stop 1
endif

! Read the lattice and the starting state: a pair of Genesis dumps (the shared-start
! benchmark methodology), or a self-generated steady-state condition when both file
! names are blank.
!
! FEL elements carry tracking_method = custom, and Bmad's bookkeeping (the reference
! time/energy pass inside bmad_parser, any track1 at the seam) resolves custom tracking
! through track1_custom_ptr: point it at the standard periodic-wiggler kernel, so the
! element behaves as the plain Bmad wiggler it is everywhere EXCEPT inside this
! driver's own FEL walk. In particular the reference time acquires the resonant
! undulation delay of brief 7.5 from Bmad's own code, not from anything written here.

! Both hooks are needed: mat6_calc_method resolves to custom too (auto follows the
! tracking method), and make_mat6 calls through a null make_mat6_custom_ptr otherwise.

track1_custom_ptr => fel_ele_as_wiggler
make_mat6_custom_ptr => fel_mat6_as_wiggler

! The FEL mode and unaveraged parameters live on the LATTICE (Stage A of the CUSTOM
! retirement; manual sec:element), registered program-side so no lattice declares
! them. The same slot index serves wigglers and undulators.

call set_custom_attribute_name ('WIGGLER::FEL_TRACKING', err, 1)
if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_TRACKING', err, 1)
if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_STEPS_PER_PERIOD', err, 2)
if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_STEPS_PER_PERIOD', err, 2)
if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_RAMP_PERIODS', err, 3)
if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_RAMP_PERIODS', err, 3)
if (err) then
  print '(a)', 'fel_track_test: could not register the FEL lattice attributes.'
  stop 1
endif

! err_flag matters: bmad_parser reports attribute errors (e.g. a wake on an element
! type that cannot carry one) and RETURNS -- without the check the run continues on a
! partial lattice, found when an example's drift wakes silently never attached.

call bmad_parser (lat_file, lat, err_flag = err)
if (err) then
  print '(a)', 'fel_track_test: lattice parse errors (above); refusing to run on a partial lattice.'
  stop 1
endif
branch => lat%branch(0)

gamma0 = branch%ele(0)%value(e_tot$) / m_electron
print '(a, f0.6, a)', 'fel_track_test: gamma0 = ', gamma0, ' (from the lattice e_tot).'

call setup_fel_elements ()   ! Needs only the parsed lattice; sets two_pol for the
                             ! beam/field construction below.

! The field set (manual sec:field-set). Validate the harmonics request, then allocate
! the set and point the fundamental aliases at entry 1 BEFORE any field construction:
! every construction path below builds the fundamental through wf.

n_harm = count(harmonics /= 0)
if (harmonics(1) /= 1) then
  print '(a)', 'fel_track_test: harmonics(1) must be 1 -- the FUNDAMENTAL anchors the phase,'
  print '(a)', '  the phi0 advance and the slippage schedule; harmonic fields ride on it.'
  stop 1
endif
do ih = 2, 9
  if (harmonics(ih) == 0) then
    if (any(harmonics(ih:) /= 0)) then
      print '(a)', 'fel_track_test: harmonics must be a gap-free increasing list (0 padding at the end).'
      stop 1
    endif
    exit
  endif
  if (harmonics(ih) <= harmonics(ih-1)) then
    print '(a)', 'fel_track_test: harmonics must be strictly increasing (no duplicates).'
    stop 1
  endif
enddo
if (n_harm > 1 .and. any(fel_mode == fel_unaveraged$ .and. is_fel)) then
  print '(a)', 'fel_track_test: harmonic FIELDS with an UNAVERAGED element are not implemented'
  print '(a)', '  (the unaveraged mode carries the fundamental envelope only; its harmonic'
  print '(a)', '  couplings are validated through the particle spectra, not a carried field).'
  stop 1
endif
if (n_harm > 1 .and. two_pol) then
  print '(a)', 'fel_track_test: harmonic fields with TWO LIVE POLARIZATIONS are not validated'
  print '(a)', '  together yet; run one or the other.'
  stop 1
endif
select case (wavefront_format)
case ('genesis', 'openpmd', 'both')
case default
  print '(a)', 'fel_track_test: wavefront_format must be genesis, openpmd or both, not "' // &
                trim(wavefront_format) // '".'
  stop 1
end select

allocate (ffield(n_harm))
do ih = 1, n_harm
  ffield(ih)%harm = harmonics(ih)
enddo
allocate (n_banked(n_harm), esc_id(n_harm))
n_banked = 0
esc_id = 0
wf => ffield(1)%wf
slip => ffield(1)%slip
bank => ffield(1)%bank

! ONE reference energy, and the lattice is it: gamma0 = e_tot/m_e c^2 from the lattice
! header, never a namelist input. There used to be a namelist gamma0 for Genesis-deck
! symmetry; the first external user fed it a hand-rounded value against a round lattice
! e_tot, the two disagreed at 1.4e-9, and the run died mid-tracking on the seam's
! backstop p0c check with raw numbers -- the FEL physics ran on one reference while
! Bmad's momenta were normalized by the other. Two specifications of one truth is the
! defect; the redundant one was removed (the deliverable-9 rule: parameters live on
! the lattice).

if ((beam_file == '') .neqv. (field_file(1) == '')) then
  print '(a)', 'fel_track_test: give both beam_file and field_file, or neither (to generate).'
  stop 1
endif
if (field_file(1) == '' .and. any(field_file(2:) /= '')) then
  print '(a)', 'fel_track_test: harmonic field files need the fundamental in field_file(1).'
  stop 1
endif
if (beam_file /= '' .and. (dist_file /= '' .or. use_beam_init)) then
  print '(a)', 'fel_track_test: dump files and a distribution import are mutually exclusive.'
  stop 1
endif
if (dist_file /= '' .and. use_beam_init) then
  print '(a)', 'fel_track_test: give dist_file or use_beam_init, not both.'
  stop 1
endif

if (beam_file /= '') then
  call fel_read_genesis4_beam (fbeam, beam_file, gamma0, err)
  if (err) stop 1
  if (wavefront_file_is_openpmd(field_file(1))) then
    call read_openpmd_into_field (field_file(1), 1)
  else
    call wavefront_read_genesis4 (wf, field_file(1), err)
    if (err) stop 1
  endif
  if (two_pol .and. .not. allocated(wf%Ey)) then
    allocate (wf%Ey(size(wf%Ex,1), size(wf%Ex,2), size(wf%Ex,3)))
    wf%Ey = 0
  endif
elseif (dist_file /= '' .or. use_beam_init) then
  call import_initial_state ()
else
  call generate_initial_state ()
endif

! Harmonic fields: the fundamental's grid and window, its wavelength / h, dark. An
! openPMD import (field_file entries 2+) fills the entry whose photon energy it
! carries; the fundamental's grid must match, per the same one-window rule the
! fundamental import obeys.

do ih = 2, n_harm
  call wavefront_init (ffield(ih)%wf, size(wf%Ex,1), size(wf%Ex,2), size(wf%Ex,3), &
                       wf%dx, wf%dy, wf%dz, wf%wavelength / ffield(ih)%harm, 'x', wf%ref_position)
enddo

do is = 2, 9
  if (field_file(is) == '') cycle
  if (.not. wavefront_file_is_openpmd(field_file(is))) then
    print '(a)', 'fel_track_test: harmonic field files must be openPMD EXT_Wavefront'
    print '(a)', '  (the Genesis format carries no photon energy to match on): ' // trim(field_file(is))
    stop 1
  endif
  call import_harmonic_field (field_file(is))
enddo

if (split_weights) call do_split_weights (fbeam)
nslice = size(fbeam%slice)
ks = twopi / wf%wavelength

! The beam and field must describe the same time window: one field slice per beam slice,
! at the same wavelength. Checked, never assumed (FINDINGS.md section 5).

if (size(wf%Ex, 3) /= nslice) then
  print '(2(a, i0))', 'fel_track_test: beam has ', nslice, ' slices but the field has ', size(wf%Ex, 3)
  stop 1
endif
if (abs(wf%wavelength - fbeam%wavelength) > 1e-12_rp * fbeam%wavelength) then
  print '(a, 2es20.12)', 'fel_track_test: beam and field disagree on the wavelength: ', &
                         fbeam%wavelength, wf%wavelength
  stop 1
endif

! Collective effects (deliverable 8): configure, build the wake kernels, hoist the
! convolution once (Genesis's behavior; recomputed at the migration stride when
! migration can change the currents). The per-slice eloss is written to
! <out_root>.wake.txt so the energy-bookkeeping check can check the applied loss exactly.

coll%efield%on = (sc_nz >= 1 .or. sc_longrange)
coll%efield%rmax = sc_rmax
coll%efield%ngrid = sc_ngrid
coll%efield%nz = sc_nz
coll%efield%nphi = sc_nphi
coll%efield%longrange = sc_longrange

if (coll%efield%on .and. interlude_model == 'bmad') then
  print '(a)', 'fel_track_test: NOTE space charge acts inside undulators and genesis-model'
  print '(a)', '  interludes only; the Bmad seam''s own space charge is deliverable 9''s domain.'
endif

coll%wake%on = wake_on
coll%wake%loss = wake_loss
coll%wake%radius = wake_radius
coll%wake%conductivity = wake_conductivity
coll%wake%relaxation = wake_relaxation
coll%wake%roundpipe = wake_roundpipe
coll%wake%material = wake_material
coll%wake%gap = wake_gap
coll%wake%lgap = wake_lgap
coll%wake%hrough = wake_hrough
coll%wake%lrough = wake_lrough

if (wake_on) then
  call fel_wake_init (coll%wake, nslice, nint(fbeam%slice_spacing / fbeam%wavelength), &
                      fbeam%wavelength, err)
  if (err) stop 1
  call fel_wake_update (coll%wake, fbeam)
  open (newunit = iu_wake, file = trim(out_root) // '.wake.txt', action = 'write')
  call write_wake_block (0.0_rp)

  ! The single-particle kernels, for the deliverable-11 cross-validation: the same
  ! physical wake fed to Bmad's z_long machinery must reproduce these kernels'
  ! convolution. NOTE the s = 0 entries carry the Bane self-slice half factor
  ! (fel_wake_init halves them); a plain W(z) table wants the unhalved value.

  if (write_wake_kernels /= '') then
    block
      integer iu_k, i_k
      open (newunit = iu_k, file = trim(write_wake_kernels), action = 'write')
      write (iu_k, '(a)') '# s [m]   wakeres   wakegeo   wakerou   [eV/(m electron)]; s=0 rows are HALVED (Bane self-slice)'
      do i_k = 1, coll%wake%ns
        write (iu_k, '(4es24.15e3)') coll%wake%ds * (i_k-1), coll%wake%wakeres(i_k), &
                                     coll%wake%wakegeo(i_k), coll%wake%wakerou(i_k)
      enddo
      close (iu_k)
      print '(a)', 'fel_track_test: wrote wake kernels: ' // trim(write_wake_kernels)
    end block
  endif
endif

if (write_initial .or. load_only) then
  call wavefront_write_genesis4 (wf, trim(out_root) // '-initial.fld.h5', err, 'x')
  if (err) stop 1
  call fel_write_genesis4_beam (fbeam, trim(out_root) // '-initial.par.h5', err)
  if (err) stop 1
endif

if (load_only) then
  print '(a)', 'fel_track_test: load_only set; initial state written, no tracking.'
  stop 0
endif

! Time dependence follows from the dumps: more than one slice makes a time-dependent run
! with slippage active; one slice is the steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
do ih = 1, n_harm       ! Same values for every field: one window, lockstep rotation
                        ! in fundamental-wavelength units (Genesis's one Control::sample).
  ffield(ih)%slip%timerun = timerun
  ffield(ih)%slip%sample = fbeam%slice_spacing / fbeam%wavelength
enddo
gamma0_ref = fel_gamma0(fbeam)

! FEL segments are real Bmad wiggler/undulator elements carrying
! tracking_method = custom -- Bmad's own semantics for "the program supplies the
! tracking", which this driver does. Their FEL parameters come from LATTICE ATTRIBUTES,
! not the namelist (deliverable 9): aw (rms, Genesis's convention) derives from b_max
! and l_period through K = c*b_max/(k_u*m_e c^2) -- reference-energy independent -- with
! aw = K for a helical device and K/sqrt(2) for a planar one; helicity from field_calc;
! Genesis's natural-focusing split kx/ky from the helicity defaults (0.5/0.5 helical,
! 0/1 planar; Bmad's kx roll-off attribute is not yet mapped and must be zero). The 7.5
! assertions are enforced here: a wiggler with zero b_max or l_period would silently get
! factor = 0 in Bmad's own kernel (no resonance, no error), and a fieldmap field_calc
! gets osc_amplitude without focusing -- both are refused by name.
! (setup_fel_elements itself runs right after the parse, before the beam and field are
! built: two_pol -- does any element tilt? -- must be known when the field is made.)

call check_wake_window ()
any_unavg = any(fel_mode == fel_unaveraged$ .and. is_fel)

! The unaveraged mode is a verification mode (fel-physics.tex sec:unaveraged): the
! collective terms are not wired into its step, and a mixed line would apply them in
! some segments and silently drop them in others. Refuse by name.

if (any_unavg .and. (wake_on .or. sc_nz >= 1 .or. sc_longrange)) then
  print '(a)', 'fel_track_test: wakes/space charge are NOT wired into the unaveraged mode'
  print '(a)', '  (a verification mode; see fel-physics.tex sec:unaveraged). Turn them off.'
  stop 1
endif

! The rest of the schedule: drift autophasing. Interludes accumulate Lz; the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is UNGUARDED in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is --
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run. (Citations kept AT THE LINES:
! this quirk's exactness matters here, at the call site; manual sec:slippage.)

allocate (ele_slip(branch%n_ele_track))
ele_slip = 0
Lz = 0
prev_ie = 0

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%value(l$) == 0) cycle
  if (is_fel(ie)) then
    if (Lz > 0 .and. prev_ie > 0) then
      ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
      Lz = 0
    endif
  else
    Lz = Lz + ele%value(l$)
  endif
  prev_ie = ie
enddo
if (prev_ie > 0) then
  ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
endif

! Diagnostics file, one row per slice per record at Genesis's record positions, slices in
! time-window order.

open (newunit = iu_diag, file = trim(out_root) // '.diag.txt', action = 'write')
write (iu_diag, '(a, i0)') '# nslice = ', nslice
write (iu_diag, '(a, es22.14)') '# slice_spacing = ', fbeam%slice_spacing
write (iu_diag, '(a)') '#         z            slice        power         on_axis_intensity        bunching        ' // &
      'bunching_phase        mean_energy         sigma_energy          sigma_x               sigma_y' // &
      '               current               n_eff'

! The unaveraged energy ledger (fel-physics.tex sec:unaveraged): one row per record
! step inside FEL segments -- total weighted beam energy, total window field energy,
! and the kick-side energy change of the step. The ledger check holds
! d(E_beam + U_field) to its measured floor.

if (any_unavg) then
  open (newunit = iu_ledger, file = trim(out_root) // '.ledger.txt', action = 'write')
  write (iu_ledger, '(a)') '#         z          E_beam_rel [J]          U_field [J]           dE_kick [J]          U_escaped [J]          U_spont [J]         E_radiated [J]'
endif

n_moved_tot = 0
charge_dropped_tot = 0
b_dev_max = 0

! Diagnostics setup (manual sec:stats): resolve the dump-at element lists through
! Bmad's own locator (class::name syntax comes for free; an entry matching nothing is
! refused by name), precompute the EXACT record and element-end counts by replaying
! the walk's skip rule, and size the stats arrays.

if (swap_beam_xy) then
  do is = 1, nslice
    call swap_arrays (fbeam%slice(is)%x, fbeam%slice(is)%y)
    call swap_arrays (fbeam%slice(is)%px, fbeam%slice(is)%py)
  enddo
endif

call setup_diagnostics ()

! Progress goes to stdout, throttled by wall clock (slow modes print a steady trickle,
! fast runs just the element boundaries); scripts that redirect stdout get it as their
! log. The numbers come from the stats row just taken -- no extra computation.

call system_clock (prog_count0, prog_rate)
prog_count_last = prog_count0
lat_length = branch%ele(branch%n_ele_track)%s

z_now = 0
call take_stats_record (.true.)   ! Evaluates the diag instrument too; the writer prints.
call write_diag_rows()            ! Initial record, matching Genesis's diag before the first step.

! Walk the lattice.

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record: Genesis's unrolled lattice has no counterpart for them -- UNLESS a wake
  ! applies there (a zero-length wake element is a standard Bmad idiom), in which case
  ! the element takes the interlude wake path below. wake_src resolves the wake through
  ! lords: a wake on a superimposed or split element lives on the LORD, and ele%wake is
  ! null on its slaves (pointer_to_wake_ele; the resolution also picks exactly ONE
  ! slave of a split lord -- the one containing the lord's midpoint -- so a split wake
  ! applies once, Bmad's own convention). Checking ele%wake directly was the deliverable
  ! 11 hole: lord wakes fell through to the per-slice path, where Bmad applied them
  ! within single slices and noted every zero-charge filler.

  wake_src => pointer_to_wake_ele(ele)
  if (ele%value(l$) == 0 .and. .not. associated(wake_src)) cycle

  if (is_fel(ie)) then

    ! FEL segment: Genesis's unroll in the element's OWN step -- Bmad's standard
    ! ds_step/num_steps attributes, whose bookkeeper computes exactly Genesis's
    ! num_steps = round(l/ds_step) (attribute_bookkeeper.f90; there is no namelist
    ! step size, the same rule as every other parameter). Equal steps; slippage after
    ! each step's field solve; any end-of-lattice fixup lands on the last step.

    und = und_of(ie)
    und%bmad_transport = (fel_mode(ie) == fel_averaged$)   ! The default: bmad_standard's kernel.
    und_slip_step = (1 + und%aw**2) / (2 * gamma0_ref**2 * wf%wavelength)
    und%nstep = max(1, nint(ele%value(num_steps$)))
    und%dz = ele%value(l$) / und%nstep

    if (fel_mode(ie) == fel_unaveraged$) then
      ! The concatenated wake kick would meet the quiver-carrying chart mid-segment;
      ! nothing in this mode needs element wakes, so refuse rather than approximate.
      if (associated(wake_src)) then
        print '(a)', 'fel_track_test: element sr wakes are not supported in the unaveraged mode,'
        print '(a)', '  at element: ' // trim(ele%name)
        stop 1
      endif
      call fel_unavg_setup (und, ustate, ele%value(l$), und%dz, fel_spp(ie), fel_ramp(ie), err)
      if (err) stop 1
    endif

    do istep = 1, und%nstep
      if (fel_mode(ie) == fel_unaveraged$) then
        call fel_unavg_step (und, ustate, fbeam, wf, slip, und%dz, istep == 1, &
                             istep == und%nstep, dE_step, dU_step, err)
        u_spont_cum = u_spont_cum + dU_step
      else
        call fel_track_und_step (und, fbeam, ffield, coll, err)
      endif
      if (err) stop 1
      call apply_radiation ()

      ! Element sr wake, Bmad's once-per-passage convention mirrored: one kick at the
      ! step nearest mid-element, scaled to the full element length (scale_with_length
      ! uses ele's l), applied across the WHOLE window (deliverable 11). Direct kick,
      ! no transport -- this walk owns transport inside wigglers.

      if (associated(wake_src) .and. istep == (und%nstep + 1)/2) call apply_bmad_wake_kick (wake_src)

      z_now = z_now + und%dz
      if (istep == und%nstep) then
        call apply_slippage_banked (und%dz * und_slip_step + ele_slip(ie))
      else
        call apply_slippage_banked (und%dz * und_slip_step)
      endif
      if (istep == und%nstep) call do_migrate ()
      call take_stats_record (istep == und%nstep)
      if (fel_mode(ie) == fel_unaveraged$) call write_ledger_row ()
      call write_diag_rows()
      call progress_line (istep == und%nstep, istep, und%nstep)
    enddo
    call end_of_element ()

  elseif (interlude_model == 'bmad') then

    ! The seam: Bmad tracks each slice's bunch (coordinate copies in and out),
    ! wavefront_drift moves the field (every slice, rotation-invariant), and the common
    ! phase phi0 advances by the reference rate with Genesis's drift surrogate
    ! ks/(2*gamma0^2) as the reference wavenumber.
    !
    ! This slice loop is deliberately SERIAL: track1_bunch parallelizes over particles
    ! internally (track1_bunch_hom, "$OMP parallel do if (thread_safe)", on by default
    ! via global_com%mp_threading_is_safe), so the threads are already busy inside each
    ! call, and parallelizing here as well would nest. The FEL step's parallelism over
    ! slices lives in fel_track_mod.

    if (associated(wake_src)) then

      ! Wake-carrying interlude: ALL slices as one bunch in global window coordinates,
      ! through Bmad's own track1_bunch, which applies the sr wake at ds_wake with the
      ! whole window visible head to tail (deliverable 11). The per-slice path below is
      ! untouched for everything else, keeping its numerics bit-identical.

      call fel_concat_slices (fbeam, ele, wake_bunch, wake_beta0, err)
      if (err) stop 1
      call track1_bunch (wake_bunch, ele, err)
      if (err) then
        print '(2a)', 'fel_track_test: tracking error in element ', trim(ele%name)
        stop 1
      endif
      call fel_split_slices (wake_bunch, ele, fbeam, wake_beta0, .false., err)
      if (err) stop 1

    else
      do is = 1, nslice
        call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele, bunch, err)
        if (err) stop 1
        call track1_bunch (bunch, ele, err)
        if (err) then
          print '(2a)', 'fel_track_test: tracking error in element ', trim(ele%name)
          stop 1
        endif
        call fel_bunch_to_slice (bunch, ele, fbeam%slice(is), err)
        if (err) stop 1
      enddo
    endif

    fbeam%phi0 = fbeam%phi0 + ele%value(l$) * &
                    fel_phi0_rate(ks, ks * 0.5_rp / gamma0_ref**2, fel_p0_mc(fbeam))

    do ih = 1, n_harm      ! Each field diffracts at its own wavelength.
    call wavefront_drift (ffield(ih)%wf, ele%value(l$), err)
    if (err) exit
  enddo
    if (err) stop 1

    ! The chamber does not end where the undulator does: the wake's energy loss applies
    ! through seam interludes too, as one kick of the element's length (Genesis applies
    ! it every step, and an interlude is one step).

    call fel_wake_apply (coll%wake, fbeam, ele%value(l$))

    z_now = z_now + ele%value(l$)
    call apply_slippage_banked (ele_slip(ie))

    call do_migrate ()
    call take_stats_record (.true.)
    call write_diag_rows()
    call progress_line (.true., 1, 1)
    call end_of_element ()

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, ffield, coll, err)
    if (associated(wake_src)) call apply_bmad_wake_kick (wake_src)
    if (err) stop 1

    z_now = z_now + ele%value(l$)
    call apply_slippage_banked (ele_slip(ie))

    call do_migrate ()
    call take_stats_record (.true.)
    call write_diag_rows()
    call progress_line (.true., 1, 1)
    call end_of_element ()
  endif
enddo

close (iu_diag)
if (any_unavg) close (iu_ledger)
if (wake_on) close (iu_wake)

if (migrate) then
  print '(a, i0, a, es12.4, a)', 'fel_track_test: migration moved ', n_moved_tot, &
        ' particles; dropped charge ', charge_dropped_tot, ' C off the window ends.'
  if (migrate_check) then
    print '(a, es10.2)', '  worst whole-beam bunching deviation across migrations: ', b_dev_max
  endif
endif

! Final dumps. The field records are unrotated to time order first -- time window
! position is holds record slice 1 + mod(is-1+first, nslice) -- which is what
! Genesis's field writer does on the fly (manual sec:slippage).

call fel_write_genesis4_beam (fbeam, trim(out_root) // '-final.par.h5', err)
if (err) stop 1

print '(a)', 'fel_track_test done.'
print '(a)', '  ' // trim(out_root) // '.diag.txt'
call dump_field_set (trim(out_root) // '-final')
print '(a)', '  ' // trim(out_root) // '-final.par.h5'

call finalize_diagnostics ()

call wavefront_fft_free()

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Generate the starting condition from the beam_init description: quiet-start beam
! slices (one, or a time window of them), optional physical shot noise generalized to
! weights, and a Gaussian seed field (or a dark start at seed_power = 0). See the program
! header for the parameter list, the noise algorithm and its provenance, and the N_eff
! guard. The no-noise single-slice path is arithmetic-identical to the deliverable-4
! loader -- same draw order, same operations -- which the bit-identity anchors rely on.
!-

subroutine generate_initial_state ()

type (fel_slice_struct), pointer :: sl
real(rp) p0_mc, ks_l, eg_x, eg_y, u, v, x, xp, y, yp, gam, p_mc, beta, pz, theta0
real(rp) dx_grid, w_part, e0, xg, yg, wsum, w2sum, n_lambda, n_eff, floor_b2, target_b2
real(rp) phi, an, nbl, br, bi
real(rp), allocatable :: theta_work(:), beta_work(:), kick(:), cur_gen(:)
real(rp) nl_min, nl_max, neff_min, neff_max, floor_max
real(rp) spacing_gen, zlen_gen, s_i
integer ib, im, ip, mbase, ix, iy, is_g, nslice_gen, ih, nharm, n_clamp
logical flat_z
character(*), parameter :: r_name = 'fel_track_test'

!

if (lambda0 <= 0) then
  print '(a)', 'fel_track_test: generation needs lambda0 > 0 (required by decision; not defaulted'
  print '(a)', '  from the lattice resonance -- the first undulator may be off).'
  stop 1
endif

call check_beam_init_contract ()

npart_gen = beam_init%n_particle
if (npart_gen < 1 .or. nbins < 1 .or. mod(npart_gen, nbins) /= 0) then
  print '(a)', 'fel_track_test: beam_init%n_particle (macroparticles PER SLICE here) must be a'
  print '(a)', '  positive multiple of nbins.'
  stop 1
endif
if (beam_init%a_norm_emit <= 0 .or. beam_init%b_norm_emit <= 0) then
  print '(a)', 'fel_track_test: beam_init%a_norm_emit and %b_norm_emit must be positive.'
  stop 1
endif
if (beam_init%bunch_charge <= 0) then
  print '(a)', 'fel_track_test: the current is DERIVED from the description; beam_init%bunch_charge'
  print '(a)', '  must be positive (there is no current parameter).'
  stop 1
endif

! The Twiss is the LATTICE's (one specification of one truth, as with e_tot and the
! import path's init_beam_distribution): read the beginning element, refuse by name
! when a lattice carries none.

tw_beta_x = branch%ele(0)%a%beta;  tw_alpha_x = branch%ele(0)%a%alpha
tw_beta_y = branch%ele(0)%b%beta;  tw_alpha_y = branch%ele(0)%b%alpha
if (tw_beta_x <= 0 .or. tw_beta_y <= 0) then
  print '(a)', 'fel_track_test: the lattice carries no beginning Twiss (beginning[beta_a], etc.);'
  print '(a)', '  the generated quiet start is matched to the lattice, so the lattice must say.'
  stop 1
endif
if (beam_init%sig_pz < 0 .or. seed_power < 0 .or. grid_n_pts < 3 .or. grid_half_width <= 0 .or. &
    window_sample < 1) then
  print '(a)', 'fel_track_test: check beam_init%sig_pz, seed_power, grid_n_pts, grid_half_width,'
  print '(a)', '  window_sample.'
  stop 1
endif
if (seed_power > 0 .and. seed_waist_size <= 0) then
  print '(a)', 'fel_track_test: seed_waist_size must be positive when seed_power > 0.'
  stop 1
endif

! The window and the per-slice current derive from the beam_init description (manual
! sec:loading): one bulk bunch, evaluated analytically at the slice centers. The
! default window covers the described bunch (as the import derives its window from
! real particles); window_length overrides it for slippage headroom and warns when it
! clips the bunch. sig_z = 0 is the steady state -- the whole charge in one slice
! window -- and is refused by name for time-dependent windows.

flat_z = .false.
select case (trim(beam_init%distribution_type(3)))
case ('', 'RAN_GAUSS', 'ran_gauss', 'Ran_Gauss')
case ('GRID', 'grid', 'Grid')
  flat_z = .true.
case default
  print '(a)', 'fel_track_test: beam_init%distribution_type(3) must be RAN_GAUSS (Gaussian bunch)'
  print '(a)', '  or GRID (flat, Bmad''s uniform) for the quiet-start generator, got: ' // &
        trim(beam_init%distribution_type(3))
  stop 1
end select

spacing_gen = window_sample * lambda0

if (flat_z) then
  zlen_gen = beam_init%grid(3)%x_max - beam_init%grid(3)%x_min
  if (zlen_gen <= 0) then
    print '(a)', 'fel_track_test: a GRID (flat) z-plane needs beam_init%grid(3)%x_min < %x_max.'
    stop 1
  endif
elseif (beam_init%sig_z > 0) then
  zlen_gen = 8 * beam_init%sig_z          ! +-4 sigma covers the described bunch.
else
  zlen_gen = 0                            ! Steady state.
endif

if (window_length > 0) then
  if (zlen_gen == 0 .and. window_length > 1.5_rp * spacing_gen) then
    print '(a)', 'fel_track_test: sig_z = 0 (the steady-state description) is invalid for a'
    print '(a)', '  time-dependent window. Give the bunch a length (sig_z, or a GRID extent).'
    stop 1
  endif
  if (window_length < zlen_gen) then
    print '(a, es10.3, a, es10.3, a)', 'fel_track_test: WARNING: window_length = ', window_length, &
          ' m CLIPS the described bunch (', zlen_gen, ' m).'
  endif
  nslice_gen = max(1, nint(window_length / spacing_gen))
else
  nslice_gen = max(1, nint(zlen_gen / spacing_gen))
endif

if (shotnoise .and. nslice_gen < 2) then
  print '(a)', 'fel_track_test: shotnoise needs a time-dependent window, the same rule as Genesis.'
  stop 1
endif

mbase = npart_gen / nbins
if (gen_test_weights .and. mod(mbase, 2) /= 0) then
  print '(a)', 'fel_track_test: gen_test_weights needs an even number of beamlets.'
  stop 1
endif

p0_mc = sqrt(gamma0**2 - 1)
delgam_gen = (p0_mc**2 / gamma0) * beam_init%sig_pz    ! beta0*p0_mc*sig_pz.

fbeam%p0c = p0_mc * m_electron
fbeam%phi0 = 0
fbeam%wavelength = lambda0
fbeam%slice_spacing = spacing_gen
fbeam%s0 = 0
fbeam%nbins = nbins
fbeam%one4one = .false.

if (allocated(fbeam%slice)) deallocate(fbeam%slice)
allocate (fbeam%slice(nslice_gen), cur_gen(nslice_gen))

! The derived per-slice current: flat Q*c/extent inside the GRID extent; Gaussian
! profile at the slice centers, bunch centered in the window; steady state = the
! whole charge in the one slice window, I = Q*c/spacing.

if (flat_z) then
  cur_gen = 0
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (nslice_gen - 1) * spacing_gen / 2
    if (abs(s_i) <= zlen_gen / 2) cur_gen(is_g) = beam_init%bunch_charge * c_light / zlen_gen
  enddo
elseif (zlen_gen > 0) then
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (nslice_gen - 1) * spacing_gen / 2
    cur_gen(is_g) = beam_init%bunch_charge * c_light / (sqrt(twopi) * beam_init%sig_z) * &
                    exp(-s_i**2 / (2 * beam_init%sig_z**2))
  enddo
else
  cur_gen(1) = beam_init%bunch_charge * c_light / spacing_gen
endif

call ran_seed_put (ran_seed)

ks_l = twopi / lambda0
eg_x = beam_init%a_norm_emit / p0_mc  ! Normalized emittance to geometric.
eg_y = beam_init%b_norm_emit / p0_mc

allocate (theta_work(npart_gen), beta_work(npart_gen))
n_clamp = 0
nl_min = huge(1.0_rp); nl_max = 0; neff_min = huge(1.0_rp); neff_max = 0; floor_max = 0

do is_g = 1, nslice_gen
  sl => fbeam%slice(is_g)
  call fel_slice_reallocate (sl, npart_gen)
  sl%n = npart_gen
  w_part = cur_gen(is_g) * fbeam%slice_spacing / (c_light * npart_gen)

  ! Quiet start: mbase base samples, each replicated at nbins equally spaced
  ! ponderomotive phases (theta0 spread on a uniform grid within one beamlet spacing),
  ! so bunching harmonics below nbins vanish to roundoff. Weights and coordinates
  ! follow fel_read_genesis4_beam: z = beta*theta/ks with phi0 = 0,
  ! weight = I*slice_spacing/(c*npart). theta and beta are held in work arrays so noise
  ! can kick the phases before the z conversion.

  ip = 0
  do ib = 1, mbase
    call ran_gauss (u);  call ran_gauss (v)
    x  = sqrt(eg_x * tw_beta_x) * u
    xp = sqrt(eg_x / tw_beta_x) * (v - tw_alpha_x * u)
    call ran_gauss (u);  call ran_gauss (v)
    y  = sqrt(eg_y * tw_beta_y) * u
    yp = sqrt(eg_y / tw_beta_y) * (v - tw_alpha_y * u)

    call ran_gauss (u)
    gam = gamma0 + delgam_gen * u
    p_mc = sqrt(gam**2 - 1)
    beta = p_mc / gam
    pz = (p_mc - p0_mc) / p0_mc

    theta0 = (ib - 0.5_rp) * twopi / (nbins * mbase)

    do im = 0, nbins - 1
      ip = ip + 1
      theta_work(ip) = theta0 + im * twopi / nbins
      beta_work(ip) = beta
      sl%x(ip) = x;   sl%px(ip) = xp
      sl%y(ip) = y;   sl%py(ip) = yp
      sl%pz(ip) = pz
      sl%weight(ip) = w_part
    enddo
  enddo

  ! Validation knob: alternate beamlet weights 0.25x/1.75x, charge preserving, uniform
  ! within each beamlet so the quiet cancellation is untouched. Exercises every
  ! weighted-noise path (the asymmetry is strong enough that using a slice-uniform
  ! electron count where the beamlet's charge belongs mis-sets <|b|^2> by 56 percent,
  ! far outside the statistical check); not a physics input.

  if (gen_test_weights) then
    do ib = 1, mbase
      sl%weight((ib-1)*nbins+1 : ib*nbins) = &
              sl%weight((ib-1)*nbins+1 : ib*nbins) * (1 + 0.75_rp * (-1)**ib)
    enddo
  endif

  ! Bookkeeping the brief's 6.2 demands: real electrons N_lambda = charge/e, effective
  ! macroparticle number N_eff = (sum w)^2/sum w^2, both per slice.

  wsum = sum(sl%weight(1:npart_gen))
  w2sum = sum(sl%weight(1:npart_gen)**2)
  n_lambda = wsum / e_charge
  n_eff = wsum**2 / w2sum
  nl_min = min(nl_min, n_lambda);  nl_max = max(nl_max, n_lambda)
  neff_min = min(neff_min, n_eff); neff_max = max(neff_max, n_eff)

  ! Zero-current slices (Gaussian tails, outside a flat extent) carry no noise --
  ! Genesis's own zero-current skip, shared with the import.

  if (shotnoise .and. wsum > 0) then

    ! The N_eff guard: measure the pre-noise quiet floor. A representation whose floor
    ! is not far below the target 1/N_lambda cannot carry physical noise -- imposing on
    ! top would give a silently wrong startup level. The sweep covers EVERY harmonic the
    ! beamlet structure can resolve (1..nbins-1), not just the imposed ones: an
    ! unquiet weight pattern can park its floor on a harmonic the imposition never
    ! touches (an alternating within-beamlet pattern lands exactly on nbins/2, found
    ! by the guard's own mutation test) and still corrupt the dynamics through the
    ! nonlinear phase evolution.

    target_b2 = 1 / n_lambda
    floor_b2 = 0
    do ih = 1, nbins - 1
      br = 0; bi = 0
      do ip = 1, npart_gen
        br = br + sl%weight(ip) * cos(ih * theta_work(ip))
        bi = bi + sl%weight(ip) * sin(ih * theta_work(ip))
      enddo
      floor_b2 = max(floor_b2, (br**2 + bi**2) / wsum**2)
    enddo
    floor_max = max(floor_max, floor_b2 * n_lambda)

    if (floor_b2 > 0.01_rp * target_b2) then
      print '(a, i0, a)',      'fel_track_test: slice ', is_g, ': the quiet-start floor is not far below the'
      print '(a)',             '  physical shot-noise level -- this representation cannot carry the requested noise.'
      print '(a, es10.2, a, es10.2)', '  max_h |b(h)|^2 = ', floor_b2, '  vs target 1/N_lambda = ', target_b2
      print '(a, es10.2, a, es10.2)', '  N_eff = ', n_eff, '  N_lambda = ', n_lambda
      stop 1
    endif

    ! Fawley-style shot noise: fel_fawley_noise (fel_beam_mod), the ShotNoise
    ! transcription generalized to weights, shared with the distribution import so the
    ! two paths stay one implementation. Draw order is unchanged from when this block
    ! lived inline here (two ran_uniform per harmonic per beamlet, Genesis's loops).

    call fel_fawley_noise (theta_work(1:npart_gen), sl%weight(1:npart_gen), npart_gen, nbins, n_clamp)
  endif

  ! To the stored chart: z = beta*theta/ks with phi0 = 0, beta of the base sample.

  do ip = 1, npart_gen
    sl%z(ip) = beta_work(ip) * theta_work(ip) / ks_l
  enddo
enddo

deallocate (theta_work, beta_work)

if (shotnoise) then
  print '(a, i0, a)',        'fel_track_test: shot noise imposed on ', nslice_gen, ' slices.'
  print '(2(a, es10.3))',    '  N_lambda per slice: ', nl_min, ' to ', nl_max
  print '(2(a, es10.3))',    '  N_eff per slice:    ', neff_min, ' to ', neff_max
  print '(a, es10.2)',       '  worst quiet floor, |b|^2 * N_lambda: ', floor_max
  if (n_clamp > 0) then
    print '(a, i0, a)',      '  WARNING: ', n_clamp, ' beamlet draws had fewer than one real electron'
    print '(a)',             '  (nbl clamped to 1, as Genesis does silently). The noise level in those'
    print '(a)',             '  beamlets is not physical; use fewer macroparticles or more charge.'
  endif
endif

call generate_seed_field (nslice_gen)

end subroutine generate_initial_state

!------------------------------------------------------------------------------

subroutine check_beam_init_contract ()

! The quiet-start generator honors the beam_init fields in the header table and
! REFUSES BY NAME every other field that is set -- a standard structure that silently
! dropped fields would be worse than a custom one. (The import path is exempt:
! init_beam_distribution honors everything Bmad honors.) renorm_center/renorm_sigma,
! random_engine defaults and n_bunch = 0/1 are generation details with no analytic
! counterpart and are accepted at their defaults only.

character(60) bad

!

bad = ''
if (beam_init%position_file /= '')                          bad = 'position_file'
if (beam_init%a_emit /= 0 .or. beam_init%b_emit /= 0)       bad = 'a_emit/b_emit (use a_norm_emit/b_norm_emit)'
if (beam_init%dPz_dz /= 0)                                  bad = 'dPz_dz'
if (any(beam_init%center /= 0))                             bad = 'center'
if (any(beam_init%spin /= 0))                               bad = 'spin'
if (any(beam_init%center_jitter /= 0))                      bad = 'center_jitter'
if (any(beam_init%emit_jitter /= 0))                        bad = 'emit_jitter'
if (beam_init%sig_z_jitter /= 0)                            bad = 'sig_z_jitter'
if (beam_init%sig_pz_jitter /= 0)                           bad = 'sig_pz_jitter'
if (beam_init%t_offset /= 0)                                bad = 't_offset'
if (beam_init%dt_bunch /= 0)                                bad = 'dt_bunch'
if (beam_init%n_bunch > 1)                                  bad = 'n_bunch'
if (beam_init%ix_turn /= 0)                                 bad = 'ix_turn'
if (beam_init%full_6D_coupling_calc)                        bad = 'full_6D_coupling_calc'
if (beam_init%use_particle_start)                           bad = 'use_particle_start'
if (beam_init%use_t_coords)                                 bad = 'use_t_coords'
if (beam_init%file_name /= '')                              bad = 'file_name'
if (beam_init%random_engine /= '' .and. beam_init%random_engine /= 'pseudo') bad = 'random_engine'
if (beam_init%random_gauss_converter /= '' .and. beam_init%random_gauss_converter /= 'ziggurat') &
                                                            bad = 'random_gauss_converter'
if (beam_init%random_sigma_cutoff /= -1)                    bad = 'random_sigma_cutoff'
if (beam_init%species /= '' .and. beam_init%species /= 'electron') bad = 'species (electron only)'
if (trim(beam_init%distribution_type(1)) /= '' .and. trim(beam_init%distribution_type(1)) /= 'RAN_GAUSS' &
    .and. trim(beam_init%distribution_type(1)) /= 'ran_gauss') bad = 'distribution_type(1) (transverse: RAN_GAUSS only)'
if (trim(beam_init%distribution_type(2)) /= '' .and. trim(beam_init%distribution_type(2)) /= 'RAN_GAUSS' &
    .and. trim(beam_init%distribution_type(2)) /= 'ran_gauss') bad = 'distribution_type(2) (transverse: RAN_GAUSS only)'

if (bad /= '') then
  print '(a)', 'fel_track_test: beam_init%' // trim(bad) // ' is set but NOT honored by the'
  print '(a)', '  quiet-start generator (see the honored-fields table in the program header).'
  print '(a)', '  Refusing rather than silently ignoring it.'
  stop 1
endif

end subroutine check_beam_init_contract

!------------------------------------------------------------------------------

subroutine generate_seed_field (nslice_f)

! The field: a Gaussian seed at its waist in every slice, E = E0*exp(-r^2/w0^2),
! intensity 1/e^2 radius w0, integrating to seed_power; seed_power = 0 is a dark start.
! Grid convention matches Genesis's dgrid: ngrid points spanning +-dgrid,
! dx = 2*dgrid/(ngrid-1), center on axis. Shared by the built-in generator and the
! distribution import (both make their own beam, neither brings a field).

integer nslice_f, ix, iy, is_g
real(rp) dx_grid, e0, xg, yg

!

if (grid_n_pts < 3 .or. grid_half_width <= 0) then
  print '(a)', 'fel_track_test: check grid_n_pts and grid_half_width.'
  stop 1
endif
if (seed_power > 0 .and. seed_waist_size <= 0) then
  print '(a)', 'fel_track_test: seed_waist_size must be positive when seed_power > 0.'
  stop 1
endif

dx_grid = 2 * grid_half_width / (grid_n_pts - 1)
call wavefront_init (wf, grid_n_pts, grid_n_pts, nslice_f, dx_grid, dx_grid, &
                     fbeam%slice_spacing, lambda0, 'x', 0.0_rp)

if (two_pol .and. .not. allocated(wf%Ey)) then
  allocate (wf%Ey(grid_n_pts, grid_n_pts, nslice_f))
  wf%Ey = 0
endif

if (seed_power > 0 .and. seed_polarization == 'y') then
  e0 = sqrt(4 * (mu_0_vac * c_light) * seed_power / (pi * seed_waist_size**2))
  do iy = 1, grid_n_pts
    yg = (iy - 1) * dx_grid - grid_half_width
    do ix = 1, grid_n_pts
      xg = (ix - 1) * dx_grid - grid_half_width
      wf%Ey(ix, iy, 1) = e0 * exp(-(xg**2 + yg**2) / seed_waist_size**2)
    enddo
  enddo
  do is_g = 2, nslice_f
    wf%Ey(:, :, is_g) = wf%Ey(:, :, 1)
  enddo
elseif (seed_power > 0) then
  e0 = sqrt(4 * (mu_0_vac * c_light) * seed_power / (pi * seed_waist_size**2))
  do iy = 1, grid_n_pts
    yg = (iy - 1) * dx_grid - grid_half_width
    do ix = 1, grid_n_pts
      xg = (ix - 1) * dx_grid - grid_half_width
      wf%Ex(ix, iy, 1) = e0 * exp(-(xg**2 + yg**2) / seed_waist_size**2)
    enddo
  enddo
  do is_g = 2, nslice_f
    wf%Ex(:, :, is_g) = wf%Ex(:, :, 1)
  enddo
endif

end subroutine generate_seed_field

!------------------------------------------------------------------------------

subroutine import_initial_state ()

! Deliverable 10: a bunch_struct -- generated from Bmad's beam_init_struct (the native
! equivalent of Genesis's &beam description) or read from an openPMD-beamphysics file
! -- resampled into FEL slices by the transcribed Genesis importdistribution method
! (fel_import_mod, where the algorithm and its provenance are documented). The seed
! field comes from the same generator as the built-in loader. The RNG-free outputs the
! exactness checks read -- the analysis moments and the per-slice current profile --
! are printed at full precision.

type (beam_struct), target :: beam_b
type (bunch_struct), pointer :: bp
real(rp) moments(11)
integer is_g, ip_g, n0
logical err_i

!

if (lambda0 <= 0) then
  print '(a)', 'fel_track_test: import needs lambda0 > 0.'
  stop 1
endif
if (window_sample < 1) then
  print '(a)', 'fel_track_test: window_sample must be a positive integer (Genesis''s sample).'
  stop 1
endif

! One seed governs the whole import: the bunch generation, the resampler's draws and
! the shot noise. Seeding AFTER generation was the first mutation this path caught in
! development -- every run then imports a different bunch, and the split-weight and
! thread-determinism checks both fail on what looks like resampler noise.

call ran_seed_put (ran_seed)

if (use_beam_init) then
  if (beam_init%n_particle < 1) then
    print '(a)', 'fel_track_test: beam_init%n_particle must be positive.'
    stop 1
  endif
  beam_init%n_bunch = 1
  call init_beam_distribution (branch%ele(0), lat%param, beam_init, beam_b, err_i)
  if (err_i) stop 1
  print '(a, i0, a)', 'fel_track_test: generated ', size(beam_b%bunch(1)%particle), &
                      ' particles from beam_init.'
else
  call hdf5_read_beam (dist_file, beam_b, err_i, branch%ele(0))
  if (err_i) stop 1
  print '(a, i0, a)', 'fel_track_test: read ', size(beam_b%bunch(1)%particle), &
                      ' particles from: ' // trim(dist_file)
endif

bp => beam_b%bunch(1)

! Check knob: coincident split-weight copies before anything downstream sees the bunch.
! The current profile (weighted sums) and the analysis moments (unweighted, over
! coincident copies) must then be bit-identical to the unsplit run.

if (imp_split_weights) then
  n0 = size(bp%particle)
  call reallocate_bunch (bp, 2*n0, save = .true.)
  do ip_g = 1, n0
    bp%particle(n0+ip_g) = bp%particle(ip_g)
    bp%particle(n0+ip_g)%charge = 2 * bp%particle(ip_g)%charge / 3
    bp%particle(ip_g)%charge = bp%particle(ip_g)%charge / 3
  enddo
endif

if (write_dist_file /= '') then
  call fel_write_genesis4_distribution (bp, write_dist_file, err_i)
  if (err_i) stop 1
  print '(a)', 'fel_track_test: wrote Genesis distribution file: ' // trim(write_dist_file)
endif

if (write_opmd_file /= '') then
  call hdf5_write_beam (write_opmd_file, beam_b%bunch(1:1), .false., err_i, lat)
  if (err_i) stop 1
  print '(a)', 'fel_track_test: wrote openPMD-beamphysics file: ' // trim(write_opmd_file)
endif

! imp%npart and imp%nbins come from the imp block directly (the resample's own knobs;
! beam_init%n_particle is the BUNCH particle count on this path).
call fel_import_bunch (bp, gamma0, lambda0, window_sample * lambda0, imp, fbeam, err_i, moments)
if (err_i) stop 1

print '(a, i0, a, i0, a)', 'fel_track_test: imported into ', size(fbeam%slice), &
                           ' slices of ', imp%npart, ' particles.'
print '(a, 11es24.15e3)', 'import moments (gavg xavg pxavg yavg pyavg ex ey bx by ax ay):', moments
do is_g = 1, size(fbeam%slice)
  print '(a, i0, a, es24.15e3)', 'import current ', is_g, ': ', &
        c_light * sum(fbeam%slice(is_g)%weight(1:fbeam%slice(is_g)%n)) / fbeam%slice_spacing
enddo

call generate_seed_field (size(fbeam%slice))

end subroutine import_initial_state

!------------------------------------------------------------------------------

subroutine do_split_weights (beam)

! Replace each particle by two coincident copies with weights w/3 and 2w/3. The order --
! all first copies, then all second copies -- keeps the original particles' storage
! order, which keeps the RK4 arithmetic per copy identical to the unsplit run.

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sp
integer is, ip, n0

do is = 1, size(beam%slice)
  sp => beam%slice(is)
  n0 = sp%n
  call fel_slice_reallocate (sp, 2*n0)
  do ip = 1, n0
    sp%x(n0+ip) = sp%x(ip);  sp%px(n0+ip) = sp%px(ip)
    sp%y(n0+ip) = sp%y(ip);  sp%py(n0+ip) = sp%py(ip)
    sp%z(n0+ip) = sp%z(ip);  sp%pz(n0+ip) = sp%pz(ip)
    sp%weight(n0+ip) = 2 * sp%weight(ip) / 3
    sp%weight(ip) = sp%weight(ip) / 3
  enddo
  sp%n = 2*n0
enddo

end subroutine do_split_weights

!------------------------------------------------------------------------------

subroutine apply_bmad_wake_kick (wake_ele)

! One whole-window application of an element's Bmad sr wake, as a pure kick
! (deliverable 11): concatenate the slices into global window coordinates, let Bmad's
! own machinery order and kick (track1_sr_wake: pseudomode accumulation head to tail,
! z_long binned FFT), split back holding theta -- Genesis's convention for wake energy
! loss, the same z rescale fel_wake_apply_slice does -- so the phase every deposition
! sees is continuous through the kick. Used inside wigglers (mid-element) and after
! genesis-model interludes; Bmad-model interludes instead go through track1_bunch,
! where the wake applies at ds_wake in Bmad's own chart.

type (ele_struct) wake_ele
logical err_w

!

call fel_concat_slices (fbeam, wake_ele, wake_bunch, wake_beta0, err_w)
if (err_w) stop 1
call order_particles_in_z (wake_bunch)
call track1_sr_wake (wake_bunch, wake_ele)
call fel_split_slices (wake_bunch, wake_ele, fbeam, wake_beta0, .true., err_w)
if (err_w) stop 1

end subroutine apply_bmad_wake_kick

!------------------------------------------------------------------------------

subroutine setup_fel_elements ()

! Recognize FEL segments -- wiggler/undulator elements with tracking_method = custom,
! Bmad's own semantics for program-supplied tracking, which this driver is -- and derive
! their FEL parameters from lattice attributes, enforcing the brief's 7.5 assertions.
! The stored k1x/k1y wiggler attributes are deliberately NOT read (7.5: their helical
! sign disagrees with the tracking locals; nothing here cross-uses them).

type (ele_struct), pointer :: w
integer je
real(rp) kw, kk, rv
logical err_a

allocate (is_fel(branch%n_ele_track), und_of(branch%n_ele_track))
allocate (fel_mode(branch%n_ele_track), fel_spp(branch%n_ele_track), fel_ramp(branch%n_ele_track))
is_fel = .false.
fel_mode = 0;  fel_spp = 0;  fel_ramp = 0

do je = 1, branch%n_ele_track
  w => branch%ele(je)
  if (.not. (w%key == wiggler$ .or. w%key == undulator$)) cycle
  if (w%tracking_method /= custom$) cycle

  ! The FEL mode and unaveraged parameters, from the element's own attributes.
  ! fel_tracking: unset/0 = averaged with the bmad_standard kernel's transverse maps
  ! (Bmad's own kernel is the default); 1 = unaveraged; -1 = averaged with the
  ! transcribed-Genesis maps (validation-internal: the Genesis tiers require
  ! transcription-level transport; no production lattice writes it).

  rv = value_of_attribute(w, 'FEL_TRACKING', err_a)
  if (err_a) stop 1
  fel_mode(je) = nint(rv)
  if (abs(rv - fel_mode(je)) > 1e-9_rp .or. fel_mode(je) < fel_transcribed$ .or. &
      fel_mode(je) > fel_unaveraged$) then
    print '(a)', 'fel_track_test: fel_tracking must be -1 (transcribed maps, validation-internal),'
    print '(a)', '  0/unset (averaged, bmad_standard kernel maps) or 1 (unaveraged), at element: ' // trim(w%name)
    stop 1
  endif

  rv = value_of_attribute(w, 'FEL_STEPS_PER_PERIOD', err_a)
  if (err_a) stop 1
  fel_spp(je) = nint(rv)
  if (fel_spp(je) == 0) fel_spp(je) = 20
  if (fel_spp(je) < 10) then
    print '(a)', 'fel_track_test: fel_steps_per_period is below the floor of 10 (MINERVA''s envelope),'
    print '(a)', '  at element: ' // trim(w%name)
    stop 1
  endif

  ! fel_ramp_periods: an attribute's unset value is 0, and a silent hard edge would
  ! reintroduce the K/gamma handoff hazard by omission -- so unset/0 means the default
  ! of 2 periods, and a TRUE hard edge (the mutation/test configuration) must be asked
  ! for by name with the explicit sentinel -1.

  rv = value_of_attribute(w, 'FEL_RAMP_PERIODS', err_a)
  if (err_a) stop 1
  if (rv == 0) then
    fel_ramp(je) = 2
  elseif (rv == -1) then
    fel_ramp(je) = 0
  elseif (rv < 0) then
    print '(a)', 'fel_track_test: fel_ramp_periods must be positive, 0/unset (default 2), or the'
    print '(a)', '  hard-edge test sentinel -1, at element: ' // trim(w%name)
    stop 1
  else
    fel_ramp(je) = rv
  endif

  ! The 7.5 assertions live in fel_assert_wiggler_sane, ONE authority called from the
  ! track1/mat6 hooks (where they fire first, during the parse) and again here. Keeping
  ! a second copy inline was tried and rejected: redundant assertions mask the removal
  ! of either copy, which defeats mutation testing of the refusal checks.

  call fel_assert_wiggler_sane (w)

  is_fel(je) = .true.
  kw = twopi / w%value(l_period$)

  ! aw (rms, Genesis's convention) from the peak field:
  ! K = c*b_max/(k_u * m_e c^2), exactly and independent of the reference energy;
  ! helical aw = K, planar aw = K/sqrt(2). Focusing split: Genesis's defaults by
  ! helicity, scaled by ku^2 as Genesis's unroll does (manual sec:element).

  kk = c_light * w%value(b_max$) / (kw * m_electron)

  und_of(je)%ku = kw
  und_of(je)%helical = (w%field_calc == helical_model$)
  if (und_of(je)%helical) then
    und_of(je)%aw = kk
    und_of(je)%kx = 0.5_rp * kw**2
    und_of(je)%ky = 0.5_rp * kw**2
  else
    und_of(je)%aw = kk / sqrt(2.0_rp)
    und_of(je)%kx = 0
    und_of(je)%ky = kw**2
  endif

  ! Tilt: the wiggle-plane rotation (planar only -- a tilted helical is a no-op that
  ! reads as confusion, refused; the transcribed-Genesis maps know no tilt, refused).
  ! The polarization 2-vector on (Ex, Ey): planar (cos t, sin t); helical (1,-i)/sqrt2.

  und_of(je)%tilt = w%value(tilt_tot$)
  if (und_of(je)%tilt /= 0) then
    if (und_of(je)%helical) then
      print '(2a)', 'fel_track_test: tilt on a HELICAL FEL element is a rotation of a ', &
                    'circularly symmetric field -- a no-op that reads as a mistake: ' // trim(w%name)
      stop 1
    endif
    if (fel_mode(je) == fel_transcribed$) then
      print '(2a)', 'fel_track_test: the transcribed-Genesis maps (fel_tracking = -1) know ', &
                    'no tilt (Genesis has none); use the default maps on: ' // trim(w%name)
      stop 1
    endif
  endif
  und_of(je)%cos_t = cos(und_of(je)%tilt)
  und_of(je)%sin_t = sin(und_of(je)%tilt)
  if (und_of(je)%helical) then
    und_of(je)%pol = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.0_rp, -1.0_rp, rp)] / sqrt(2.0_rp)
  else
    und_of(je)%pol = [cmplx(und_of(je)%cos_t, 0.0_rp, rp), cmplx(und_of(je)%sin_t, 0.0_rp, rp)]
  endif
enddo

two_pol = seed_polarization == 'y'
do je = 1, branch%n_ele_track
  if (is_fel(je) .and. und_of(je)%sin_t /= 0) two_pol = .true.
enddo
if (seed_polarization /= 'x' .and. seed_polarization /= 'y') then
  print '(a)', 'fel_track_test: seed_polarization must be "x" or "y".'
  stop 1
endif

if (.not. any(is_fel) .and. .not. reference_run) then
  print '(a)', 'fel_track_test: the lattice has no FEL elements (wiggler/undulator with tracking_method = custom).'
  print '(a)', '  Set reference_run = T for a deliberate no-FEL run (Bmad tracks everything).'
  stop 1
endif



end subroutine setup_fel_elements

!------------------------------------------------------------------------------
! Element sr wakes act across the WHOLE window (deliverable 11): all slices
! concatenate into one bunch in global window coordinates and Bmad's wake machinery
! applies unmodified. What is checked here, by name: lr (multi-bunch) wakes are not
! supported; a pseudomode wake whose z_max is shorter than the window would have Bmad
! kill the bunch mid-run; a z_long table narrower than the window would overflow its
! binning grid the same way. Runs AFTER the beam is built (the window length is the
! subject); setup_fel_elements runs before it (two_pol must precede the field).

subroutine check_wake_window ()

type (ele_struct), pointer :: w
integer je

do je = 1, branch%n_ele_track
  w => pointer_to_wake_ele(branch%ele(je))
  if (.not. associated(w)) cycle
  if (allocated(w%wake%lr%mode)) then
    if (size(w%wake%lr%mode) > 0) then
      print '(2a)', 'fel_track_test: lr (multi-bunch) wakes are not supported; ', &
                    'remove them from: ' // trim(w%name)
      stop 1
    endif
  endif
  if (w%wake%sr%z_max > 0 .and. size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_max) then
    print '(2a)', 'fel_track_test: the time window is longer than this element''s sr wake ', &
                  'z_max can handle: ' // trim(w%name)
    stop 1
  endif
  if (w%wake%sr%z_long%dz > 0 .and. &
      size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_long%z0) then
    print '(2a)', 'fel_track_test: the time window is longer than this element''s z_long ', &
                  'wake table extent z0: ' // trim(w%name)
    stop 1
  endif
enddo

end subroutine check_wake_window

!------------------------------------------------------------------------------

subroutine do_migrate ()

! Slice migration at the per-element stride, serial, between the parallel regions (the
! thread check stays untouched). Called AFTER z_now is advanced, so per-event drop
! reports carry the z of the diagnostic record they precede -- the conservation
! timeline reconstructs exactly from the log. With migrate_check, the whole-beam
! weighted phasor S = sum(w e^{i theta}) must satisfy S_before = S_after + S_dropped to
! rounding: every mover's phase shifts by an exact multiple of 2*pi*sample and a drop
! removes exactly its own term, so any deviation beyond rounding is a bookkeeping bug
! (wrong z adjustment, weight not moved), not statistics.

real(rp) chd, sb_re, sb_im, sa_re, sa_im, d_re, d_im, wsum
integer nm

if (.not. migrate) return

if (migrate_check) call whole_beam_phasor (sb_re, sb_im, wsum)

call fel_migrate_slices (fbeam, ks, nm, chd, d_re, d_im, err)
if (err) stop 1

n_moved_tot = n_moved_tot + nm
charge_dropped_tot = charge_dropped_tot + chd

! Migration changes the current profile, which the wake convolution was hoisted on
! (brief 4.3's premise predates migration): recompute at this stride. Every recompute
! appends a z-stamped block to <out_root>.wake.txt, so "the wake followed the currents"
! is a structural fact a check can parse without reimplementing the convolution.

if (nm > 0 .and. coll%wake%on) then
  call fel_wake_update (coll%wake, fbeam)
  call write_wake_block (z_now)
endif
if (chd > 0) then
  print '(a, es22.14, a, es22.14, a)', 'fel_track_test: migration dropped ', chd, &
                                       ' C off the window ends at z = ', z_now, ' m.'
endif

if (migrate_check .and. (nm > 0 .or. chd > 0)) then
  call whole_beam_phasor (sa_re, sa_im, wsum)
  if (wsum > 0) then
    b_dev_max = max(b_dev_max, sqrt((sb_re - sa_re - d_re)**2 + (sb_im - sa_im - d_im)**2) / wsum)
  endif
endif

end subroutine do_migrate

!------------------------------------------------------------------------------

subroutine whole_beam_phasor (s_re, s_im, wsum)

! The whole-beam weighted phasor sum(w e^{i theta}) and total weight, all slices.
! Exactly conserved across migration (moves shift phases by 2*pi*sample multiples;
! drops are accounted separately), which is what migrate_check verifies.

real(rp) s_re, s_im, wsum, theta, beta, p0_mc, w
integer is, ip

s_re = 0; s_im = 0; wsum = 0
p0_mc = fel_p0_mc(fbeam)

do is = 1, size(fbeam%slice)
  do ip = 1, fbeam%slice(is)%n
    w = fbeam%slice(is)%weight(ip)
    beta = fel_beta_of(p0_mc, fbeam%slice(is)%pz(ip))
    theta = fbeam%phi0 + ks * fbeam%slice(is)%z(ip) / beta
    s_re = s_re + w * cos(theta)
    s_im = s_im + w * sin(theta)
    wsum = wsum + w
  enddo
enddo

end subroutine whole_beam_phasor

!------------------------------------------------------------------------------

subroutine write_wake_block (z)

! One block of per-slice eloss, z-stamped. Written at the hoisted update and at every
! migration-stride recompute; the energy-bookkeeping and stale-wake checks parse these.

real(rp) z
integer is_w

write (iu_wake, '(a, es22.14)') '# z = ', z
do is_w = 1, nslice
  write (iu_wake, '(i8, es24.16)') is_w, coll%wake%eloss(is_w)
enddo

end subroutine write_wake_block

!------------------------------------------------------------------------------

subroutine write_diag_rows ()

! One row per slice, slices in time-window order: beam slice is against field slice
! fel_field_index(slip, is, nslice), the unrotation of manual sec:slippage. The values
! are the ones the stats loop just evaluated with the SAME fel_field_diag and
! fel_slice_diag calls, slice-parallel (each slice's arithmetic identical to the old
! serial sweep, so this file is bit-for-bit what it always was); this routine only
! prints. take_stats_record must have run for this record first.

integer is

do is = 1, nslice
  write (iu_diag, '(es24.16, i8, 10es24.16)') z_now, is, fpow_arr(is), fonax_arr(is), &
        bdiag_arr(is)%bunching, bdiag_arr(is)%bunching_phase, bdiag_arr(is)%mean_energy, &
        bdiag_arr(is)%sigma_energy, bdiag_arr(is)%sigma_x, bdiag_arr(is)%sigma_y, &
        bdiag_arr(is)%current, bdiag_arr(is)%n_eff
enddo

end subroutine write_diag_rows

!------------------------------------------------------------------------------
subroutine write_ledger_row ()

! One energy-ledger row (unaveraged mode): beam energy RELATIVE to the reference,
! sum(w*(gamma-gamma0))*me [C*eV = J] -- relative so the per-record change is not
! differenced off a large baseline at its own summation-rounding floor (the
! FINDINGS 4.8 lesson) -- total window field energy sum(P_is)*slice_spacing/c [J],
! and the kick-side change dE_step returned by fel_unavg_step. The last column is the
! cumulative energy transmitted out of the window by slippage (banked at the zero fill
! in fel_apply_slippage), the cumulative spontaneous deposit energy sum|dE_src|^2 (the
! one field-energy term the kick/deposit duality does not charge to the beam), and the
! cumulative energy the beam radiated away under bmad_com's radiation switches (the
! ACTUAL drawn sums, not expectations, so closure stays exact; zero with the switches
! off). In a time-dependent run the window is an open system, and the EXACTLY closing
! quantity is E_beam + U_field + U_escaped - U_spont + E_radiated. Wakes would be a
! second, unbanked exit channel; this ledger only exists where they are refused.

real(rp) e_beam, u_field, p_mc_l, g0_l
integer is, ip

! The field power per slice comes from the stats loop's fel_field_diag evaluation
! (same call, same unrotation, same summation order over slices -- bit-identical);
! take_stats_record must have run for this record first.

e_beam = 0
u_field = 0
g0_l = fel_gamma0(fbeam)
do is = 1, nslice
  do ip = 1, fbeam%slice(is)%n
    p_mc_l = fel_p0_mc(fbeam) * (1 + fbeam%slice(is)%pz(ip))
    e_beam = e_beam + fbeam%slice(is)%weight(ip) * (sqrt(p_mc_l**2 + 1) - g0_l) * m_electron
  enddo
  u_field = u_field + fpow_arr(is) * fbeam%slice_spacing / c_light
enddo

write (iu_ledger, '(7es24.16)') z_now, e_beam, u_field, dE_step, slip%u_escaped, u_spont_cum, e_rad_cum

end subroutine write_ledger_row

!------------------------------------------------------------------------------
! Spontaneous radiation inside FEL elements, honoring Bmad's GLOBAL switches
! bmad_com%radiation_damping_on / %radiation_fluctuations_on -- the same switches every
! Bmad tracking path honors (interludes get theirs through track1; this covers the
! custom-tracked FEL step, whose radiation reaction cannot emerge from Newton-Lorentz:
! FINDINGS 7.27). Per record:
!
!   damping:      dgamma_j = -(2/3) r_e gamma_j^2 ku^2 aw^2 * INT g^2 ds
!                 (each particle's own gamma; the envelope integral over the record, so
!                 the unaveraged ramps radiate by their actual strength; averaged g = 1)
!   fluctuations: the standard undulator quantum-diffusion form -- variance
!                 1.015e-27 * ku^3 aw^2 F(aw) gamma0^4 per meter with Genesis's
!                 Incoherent.cpp F(aw) fits (the Saldin closed form) -- drawn GAUSSIAN.
!                 Genesis draws uniform scaled by sqrt(3) to reach this SAME variance;
!                 the sqrt(3) normalizes the uniform draw and must not appear with a
!                 Gaussian (the physical limit, a merit choice). ONE DRAW PER BEAMLET,
!                 exactly as
!                 Genesis: independent per-particle kicks would break the quiet start's
!                 per-beamlet harmonic cancellation. Draws happen SERIALLY in fixed
!                 slice order from the one seeded stream, so results are independent of
!                 thread count; only the application parallelizes.
!
! The actual drawn energy accumulates into e_rad_cum (the ledger's E_radiated column):
! drawn sums, not expectations, so the TD closure stays exact. In the unaveraged mode
! this double-counts the grid-captured band (~3% at the reference grid, bounded live by
! check_spontaneous.py) -- accepted, smaller than the angular-distribution estimate's
! own uncertainty.

subroutine apply_radiation ()

real(rp) c_loss, fform, sig, gam, p_mc_r, dgam, intg2, s0, g_env, gp_env
integer isl, ipr, ibl, nb, nbl_max
!

if (.not. (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) return

! The envelope integral over this record: the substep-grid midpoint sum for the
! unaveraged mode (matching its own integration grid), dz exactly for the averaged.

if (fel_mode(ie) == fel_unaveraged$) then
  intg2 = 0
  do ipr = 1, ustate%nsub
    s0 = (ustate%s - und%dz) + (ipr - 0.5_rp) * ustate%dsub
    g_env = fel_unavg_envelope(ustate, s0, gp_env)
    intg2 = intg2 + g_env**2 * ustate%dsub
  enddo
else
  intg2 = und%dz
endif

c_loss = (2.0_rp / 3.0_rp) * r_e * und%ku**2 * und%aw**2
if (und%helical) then
  fform = 1.42_rp * und%aw + 1 / (1 + 1.5_rp * und%aw + 0.95_rp * und%aw**2)
else
  fform = 1.697_rp * und%aw + 1 / (1 + 1.88_rp * und%aw + 0.8_rp * und%aw**2)
endif
sig = 0
if (bmad_com%radiation_fluctuations_on) then
  sig = sqrt(1.015e-27_rp * und%ku**3 * und%aw**2 * fform * gamma0_ref**4 * intg2)
endif

nb = max(1, fbeam%nbins)
if (.not. allocated(e_rad_slice)) allocate (e_rad_slice(nslice))
e_rad_slice = 0

if (bmad_com%radiation_fluctuations_on) then
  nbl_max = 0
  do isl = 1, nslice
    nbl_max = max(nbl_max, (fbeam%slice(isl)%n + nb - 1) / nb)
  enddo
  if (.not. allocated(rad_kick)) allocate (rad_kick(nbl_max, nslice))
  do isl = 1, nslice
    do ibl = 1, (fbeam%slice(isl)%n + nb - 1) / nb
      call ran_gauss (rad_kick(ibl, isl))
    enddo
  enddo
endif

!$OMP parallel do private(ipr, p_mc_r, gam, dgam)
do isl = 1, nslice
  do ipr = 1, fbeam%slice(isl)%n
    p_mc_r = fel_p0_mc(fbeam) * (1 + fbeam%slice(isl)%pz(ipr))
    gam = sqrt(p_mc_r**2 + 1)
    dgam = 0
    if (bmad_com%radiation_damping_on) dgam = dgam - c_loss * gam**2 * intg2
    if (bmad_com%radiation_fluctuations_on) dgam = dgam + sig * rad_kick((ipr-1)/nb + 1, isl)
    e_rad_slice(isl) = e_rad_slice(isl) - fbeam%slice(isl)%weight(ipr) * dgam * m_electron
    gam = gam + dgam
    fbeam%slice(isl)%pz(ipr) = (sqrt(gam**2 - 1) - fel_p0_mc(fbeam)) / fel_p0_mc(fbeam)
  enddo
enddo
!$OMP end parallel do

e_rad_cum = e_rad_cum + sum(e_rad_slice)

end subroutine apply_radiation

!------------------------------------------------------------------------------
subroutine swap_arrays (a, b)

real(rp) a(:), b(:), tmp(size(a))

tmp = a;  a = b;  b = tmp

end subroutine swap_arrays

!------------------------------------------------------------------------------
! One progress line to stdout: where the walk is and what the light and beam are
! doing, so the slow modes (the unaveraged mode runs ~30x the averaged) show signs of
! life. Element boundaries always print; inside elements a wall-clock throttle (2 s)
! keeps fast runs quiet. All numbers are read from the stats row just taken.

subroutine progress_line (at_element_end, i_step, n_step)

logical at_element_end
integer i_step, n_step
integer(8) now
real(rp) elapsed

call system_clock (now)
if (.not. at_element_end .and. real(now - prog_count_last, rp) / prog_rate < 2.0_rp) return
prog_count_last = now
elapsed = real(now - prog_count0, rp) / prog_rate

print '(a, f5.1, a, f8.3, a, i0, a, i0, 3a, i0, a, i0, a, es9.2, a, es9.2, a, f8.5, a, i0, a)', &
      'progress: ', 100 * z_now / lat_length, '%  z = ', z_now, ' m  ele ', ie, '/', &
      branch%n_ele_track, ' ', trim(ele%name), '  step ', i_step, '/', n_step, &
      '  P = ', sum(stats%f_power(:, stats%irec)), ' W  U = ', sum(stats%f_energy(:, stats%irec)), &
      ' J  <|b|> = ', sum(stats%bunching(:, stats%irec)) / nslice, '  t = ', nint(elapsed), ' s'

end subroutine progress_line

!------------------------------------------------------------------------------
! Diagnostics (manual sec:stats). setup_diagnostics resolves the dump-at lists through
! Bmad's own lat_ele_locator (class::name syntax for free; an entry matching nothing is
! refused by name) and precomputes the EXACT record and element-end counts by replaying
! the walk's skip rule -- the stats arrays are sized once, never grown.

subroutine setup_diagnostics ()

type (ele_pointer_struct), allocatable :: eles(:)
integer i, j, n_loc
logical derr

!

allocate (dump_beam_here(0:branch%n_ele_track), dump_field_here(0:branch%n_ele_track))
dump_beam_here = .false.;  dump_field_here = .false.

do i = 1, size(dump_beam_at)
  if (dump_beam_at(i) == '') cycle
  call lat_ele_locator (dump_beam_at(i), lat, eles, n_loc, derr)
  if (derr .or. n_loc == 0) then
    print '(2a)', 'fel_track_test: dump_beam_at entry matches no element: ', trim(dump_beam_at(i))
    stop 1
  endif
  do j = 1, n_loc
    if (eles(j)%ele%ix_ele >= 1 .and. eles(j)%ele%ix_ele <= branch%n_ele_track) &
                            dump_beam_here(eles(j)%ele%ix_ele) = .true.
  enddo
enddo

do i = 1, size(dump_field_at)
  if (dump_field_at(i) == '') cycle
  call lat_ele_locator (dump_field_at(i), lat, eles, n_loc, derr)
  if (derr .or. n_loc == 0) then
    print '(2a)', 'fel_track_test: dump_field_at entry matches no element: ', trim(dump_field_at(i))
    stop 1
  endif
  do j = 1, n_loc
    if (eles(j)%ele%ix_ele >= 1 .and. eles(j)%ele%ix_ele <= branch%n_ele_track) &
                            dump_field_here(eles(j)%ele%ix_ele) = .true.
  enddo
enddo

nrec_stats = 1
nend_stats = 0
do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  wake_src => pointer_to_wake_ele(ele)
  if (ele%value(l$) == 0 .and. .not. associated(wake_src)) cycle
  nend_stats = nend_stats + 1
  if (is_fel(i)) then
    nrec_stats = nrec_stats + max(1, nint(ele%value(num_steps$)))
  else
    nrec_stats = nrec_stats + 1
  endif
enddo

call fel_stats_init (stats, nrec_stats, nend_stats, nslice, fbeam%p0c, two_pol, harmonics(2:n_harm))
allocate (bdiag_arr(nslice), fpow_arr(nslice), fonax_arr(nslice))

end subroutine setup_diagnostics

!------------------------------------------------------------------------------
subroutine take_stats_record (with_angles)

logical with_angles, serr

call fel_stats_record (stats, fbeam, ffield, z_now, with_angles, bdiag_arr, fpow_arr, fonax_arr, serr)
if (serr) stop 1

end subroutine take_stats_record

!------------------------------------------------------------------------------
! End of element: the evaluated bunch_params (the Tao end-of-element pattern) and any
! requested dumps. Mid-run field dumps unrotate exactly as the final dump does -- the
! rotation is a gauge, and re-zeroing slip%first keeps the mapping consistent.

subroutine end_of_element ()

logical eerr
character(500) fname

!

call fel_stats_element_end (stats, fbeam, ele, z_now, eerr)
if (eerr) stop 1

if (dump_beam_here(ie)) then
  write (fname, '(2a, i0, 3a)') trim(out_root), '-at', ie, '-', trim(ele%name), '.par.h5'
  call fel_write_genesis4_beam (fbeam, trim(fname), eerr)
  if (eerr) stop 1
endif

if (dump_field_here(ie)) then
  write (fname, '(2a, i0, 2a)') trim(out_root), '-at', ie, '-', trim(ele%name)
  call dump_field_set (trim(fname))
endif

end subroutine end_of_element

!------------------------------------------------------------------------------
! Write the whole field set at the given filename prefix, honoring wavefront_format:
! Genesis dumps (one polarization per file: -x/-y when Ey is live) and/or openPMD
! EXT_Wavefront (both polarizations as components of ONE mesh record). The
! fundamental keeps the pre-harmonic names; a harmonic's files carry -h<h>. Every
! record is unrotated to time order first (each field owns its rotation state; they
! move in lockstep, but the cshift must run per record).

subroutine dump_field_set (prefix)

character(*) prefix
integer ihh
logical eerr
character(8) hsuf

!

do ihh = 1, n_harm
  if (ffield(ihh)%slip%first /= 0) then
    ffield(ihh)%wf%Ex = cshift(ffield(ihh)%wf%Ex, shift = ffield(ihh)%slip%first, dim = 3)
    if (allocated(ffield(ihh)%wf%Ey)) &
        ffield(ihh)%wf%Ey = cshift(ffield(ihh)%wf%Ey, shift = ffield(ihh)%slip%first, dim = 3)
    ffield(ihh)%slip%first = 0
  endif

  hsuf = ''
  if (ffield(ihh)%harm /= 1) write (hsuf, '(a, i0)') '-h', ffield(ihh)%harm

  if (wavefront_format /= 'openpmd') then         ! genesis or both
    if (allocated(ffield(ihh)%wf%Ey)) then        ! One component per Genesis file.
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '-x.fld.h5', eerr, 'x')
      if (eerr) stop 1
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '-y.fld.h5', eerr, 'y')
      if (eerr) stop 1
      print '(a)', '  ' // prefix // trim(hsuf) // '-{x,y}.fld.h5'
    else
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '.fld.h5', eerr, 'x')
      if (eerr) stop 1
      print '(a)', '  ' // prefix // trim(hsuf) // '.fld.h5'
    endif
  endif

  if (wavefront_format /= 'genesis') then         ! openpmd or both
    call wavefront_write_openpmd (ffield(ihh)%wf, prefix // trim(hsuf) // '.wf.h5', z_now, eerr)
    if (eerr) stop 1
    print '(a)', '  ' // prefix // trim(hsuf) // '.wf.h5'
  endif
enddo

end subroutine dump_field_set

!------------------------------------------------------------------------------
! Read an openPMD wavefront into field-set entry ihh (the fundamental import path).
! The photon energy must be the fundamental's -- a file carrying a harmonic in
! field_file(1) is refused by name.

subroutine read_openpmd_into_field (fname, ihh)

character(*) fname
integer ihh
real(rp) e_photon
logical rerr

!

call wavefront_read_openpmd (ffield(ihh)%wf, fname, rerr, e_photon)
if (rerr) stop 1
if (abs(ffield(ihh)%wf%wavelength - fbeam%wavelength) > 1e-6_rp * fbeam%wavelength) then
  print '(a)', 'fel_track_test: the openPMD file in field_file(1) does not carry the FUNDAMENTAL:'
  print '(a, 2es20.12)', '  its photonEnergy wavelength vs the beam: ', &
                         ffield(ihh)%wf%wavelength, fbeam%wavelength
  stop 1
endif
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = fbeam%wavelength   ! One wavelength authority (the 1e-12
                                               ! beam/field consistency check's spirit).

end subroutine read_openpmd_into_field

!------------------------------------------------------------------------------
! Import a harmonic field: match the file's photon energy to the field-set entry
! carrying that harmonic (no match is refused by name), require the fundamental's
! grid and window, and keep the walk's wavelength convention (fundamental / h).

subroutine import_harmonic_field (fname)

character(*) fname
type (wavefront_struct) wtmp
real(rp) e_photon, e1
integer ihh
logical rerr

!

call wavefront_read_openpmd (wtmp, fname, rerr, e_photon)
if (rerr) stop 1
e1 = h_planck * c_light / fbeam%wavelength * e_charge

do ihh = 2, n_harm
  if (abs(e_photon - ffield(ihh)%harm * e1) <= 1e-6_rp * ffield(ihh)%harm * e1) exit
enddo
if (ihh > n_harm) then
  print '(a)', 'fel_track_test: the photonEnergy of ' // trim(fname)
  print '(a, es13.5, a)', '  (', e_photon, ' J) matches NO field of this run''s harmonics list.'
  stop 1
endif

if (size(wtmp%Ex,1) /= size(wf%Ex,1) .or. size(wtmp%Ex,3) /= size(wf%Ex,3) .or. &
    abs(wtmp%dx - wf%dx) > 1e-12_rp * wf%dx) then
  print '(a)', 'fel_track_test: harmonic import ' // trim(fname)
  print '(a)', '  does not match the fundamental''s grid and window (one time window, one grid).'
  stop 1
endif

call move_alloc (wtmp%Ex, ffield(ihh)%wf%Ex)
if (allocated(wtmp%Ey)) call move_alloc (wtmp%Ey, ffield(ihh)%wf%Ey)
ffield(ihh)%wf%dx = wtmp%dx;  ffield(ihh)%wf%dy = wtmp%dy
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = wf%wavelength / ffield(ihh)%harm

end subroutine import_harmonic_field

!------------------------------------------------------------------------------
subroutine apply_slippage_banked (slippage)

real(rp) slippage
integer ihh

! Every field of the set slips by the same amount in fundamental-wavelength units
! (one window, lockstep rotation; the harm argument only fixes the escape bank's
! slice light-time for a harmonic field).

if (keep_escaped_field) then
  do ihh = 1, n_harm
    call fel_apply_slippage (ffield(ihh)%slip, ffield(ihh)%wf, slippage, ffield(ihh)%bank, &
                             ffield(ihh)%harm)
    call drain_bank (ihh)
  enddo
else
  do ihh = 1, n_harm
    call fel_apply_slippage (ffield(ihh)%slip, ffield(ihh)%wf, slippage, harm = ffield(ihh)%harm)
  enddo
endif

end subroutine apply_slippage_banked

!------------------------------------------------------------------------------
! Drain the slippage bank: each transmitted slice streams to <out_root>-escaped.fld.h5
! as it leaves (peak memory a handful of grid planes), with its wavefront_params (full
! 4x4 at bank time -- exactly what the analytic free-space propagation of pulse
! statistics needs) and its transmission z. Genesis field-file conventions (dfl units)
! so existing tooling reads it; root datasets land at finalize.

subroutine drain_bank (ihh)

type (wavefront_params_struct) pms
type (fel_bank_struct), pointer :: bnk
type (wavefront_struct), pointer :: wfl
integer(hid_t) g_id
real(rp), allocatable :: work(:), gz(:,:), gp(:,:,:)
real(rp) dfl_scale
integer ihh, k, nx, h5e
logical berr
character(20) gname

!

bnk => ffield(ihh)%bank
wfl => ffield(ihh)%wf
if (bnk%n == 0) return
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
allocate (work(nx*nx))

if (esc_id(ihh) == 0) then
  call hdf5_open_file (trim(escaped_file_name(ihh)), 'WRITE', esc_id(ihh), berr)
  if (berr) stop 1
endif

do k = 1, bnk%n
  n_banked(ihh) = n_banked(ihh) + 1
  if (.not. allocated(bank_z)) then
    allocate (bank_z(1024, n_harm), bank_pms(25, 1024, n_harm))
  elseif (n_banked(ihh) > size(bank_z, 1)) then
    call move_alloc (bank_z, gz);  call move_alloc (bank_pms, gp)
    allocate (bank_z(2*size(gz,1), n_harm), bank_pms(25, 2*size(gz,1), n_harm))
    bank_z(1:size(gz,1), :) = gz;  bank_pms(:, 1:size(gz,1), :) = gp
    deallocate (gz, gp)
  endif

  call wavefront_params_of_plane (bnk%plane(:,:,k), wfl%dx, wfl%wavelength, &
                                  fbeam%slice_spacing, pms, .true., berr)
  if (berr) stop 1
  pms%s = z_now
  bank_z(n_banked(ihh), ihh) = z_now
  bank_pms(:, n_banked(ihh), ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  if (two_pol) then             ! The y component's params add to the banked energy.
    call wavefront_params_of_plane (bnk%plane_y(:,:,k), wfl%dx, wfl%wavelength, &
                                    fbeam%slice_spacing, pms, .true., berr)
    if (berr) stop 1
  endif

  write (gname, '(a, i0.6)') 'slice', n_banked(ihh)
  call H5Gcreate_f (esc_id(ihh), trim(gname), g_id, h5e)
  if (h5e < 0) stop 1
  work = dfl_scale * reshape(real(bnk%plane(:,:,k), rp), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-real', work, berr);  if (berr) stop 1
  work = dfl_scale * reshape(aimag(bnk%plane(:,:,k)), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-imag', work, berr);  if (berr) stop 1
  if (two_pol) then
    work = dfl_scale * reshape(real(bnk%plane_y(:,:,k), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real-y', work, berr);  if (berr) stop 1
    work = dfl_scale * reshape(aimag(bnk%plane_y(:,:,k)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag-y', work, berr);  if (berr) stop 1
    call hdf5_write_dataset_real (g_id, 'wavefront_params_y', &
            [pms%centroid, reshape(pms%sigma, [16]), pms%energy, pms%power, &
             pms%on_axis_intensity, pms%emit_x, pms%emit_y], berr);  if (berr) stop 1
  endif
  call hdf5_write_dataset_real (g_id, 'z_transmit', [z_now], berr);  if (berr) stop 1
  call hdf5_write_dataset_real (g_id, 'wavefront_params', bank_pms(:, n_banked(ihh), ihh), berr);  if (berr) stop 1
  call H5Gclose_f (g_id, h5e)
enddo

end subroutine drain_bank

!------------------------------------------------------------------------------
! Escaped-file name of field ihh: the fundamental keeps its pre-harmonic name.

function escaped_file_name (ihh) result (fname)

integer ihh
character(420) fname

if (ffield(ihh)%harm == 1) then
  fname = trim(out_root) // '-escaped.fld.h5'
else
  write (fname, '(2a, i0, a)') trim(out_root), '-escaped-h', ffield(ihh)%harm, '.fld.h5'
endif

end function escaped_file_name

!------------------------------------------------------------------------------
! Finalize: the stats file always; with keep_escaped_field also the escaped file's root
! datasets and the FULL PULSE at the exit plane -- each banked slice read back, free-
! space propagated over z_end - z_transmit (transmitted light is fixed information;
! undulator vacuum is free space for light with no beam under it), and concatenated
! above the live window: earliest-transmitted light is furthest ahead, so banked slice
! k lands at pulse index nslice + (n_banked - k + 1). The caller has already unrotated
! the live window (the final-dump cshift).

subroutine finalize_diagnostics ()

integer nx, h5e, ihh
logical ferr
character(420) fname_h

!

call fel_stats_write (stats, trim(out_root) // '.stats.h5', ferr)
if (ferr) stop 1
print '(a)', '  ' // trim(out_root) // '.stats.h5'

if (.not. keep_escaped_field) return

do ihh = 1, n_harm
  nx = size(ffield(ihh)%wf%Ex, 1)

  if (esc_id(ihh) /= 0) then
    call hdf5_write_dataset_int  (esc_id(ihh), 'gridpoints',   [nx],                        ferr)
    call hdf5_write_dataset_real (esc_id(ihh), 'gridsize',     [ffield(ihh)%wf%dx],         ferr)
    call hdf5_write_dataset_real (esc_id(ihh), 'refposition',  [ffield(ihh)%wf%ref_position], ferr)
    call hdf5_write_dataset_real (esc_id(ihh), 'wavelength',   [ffield(ihh)%wf%wavelength], ferr)
    call hdf5_write_dataset_int  (esc_id(ihh), 'slicecount',   [n_banked(ihh)],             ferr)
    call hdf5_write_dataset_real (esc_id(ihh), 'slicespacing', [fbeam%slice_spacing],       ferr)
    call H5Fclose_f (esc_id(ihh), h5e)
    esc_id(ihh) = 0
    print '(a)', '  ' // trim(escaped_file_name(ihh))
  endif

  ! The full pulse at the exit plane; with two live polarizations, one file per
  ! component (Genesis's format holds one). Harmonic pulses carry -h<h>.

  if (two_pol) then
    call write_pulse_file (trim(out_root) // '-pulse-x.fld.h5', .false., ihh)
    call write_pulse_file (trim(out_root) // '-pulse-y.fld.h5', .true., ihh)
    print '(a)', '  ' // trim(out_root) // '-pulse-{x,y}.fld.h5'
  elseif (ffield(ihh)%harm == 1) then
    call write_pulse_file (trim(out_root) // '-pulse.fld.h5', .false., ihh)
    print '(a)', '  ' // trim(out_root) // '-pulse.fld.h5'
  else
    write (fname_h, '(2a, i0, a)') trim(out_root), '-pulse-h', ffield(ihh)%harm, '.fld.h5'
    call write_pulse_file (trim(fname_h), .false., ihh)
    print '(a)', '  ' // trim(fname_h)
  endif
enddo

end subroutine finalize_diagnostics

!------------------------------------------------------------------------------
! One component of the full exit-plane pulse (Genesis's format holds one per file):
! the live window's slices, then each banked slice read back from the escaped file
! and free-space propagated over z_end - z_transmit. Shared by the -x and -y files.

subroutine write_pulse_file (fname, use_y, ihh)

character(*) fname
logical use_y
integer ihh
type (wavefront_struct), pointer :: wfl
integer(hid_t) p_id, e_id, g_id
real(rp), allocatable :: work(:), re_w(:), im_w(:)
real(rp) dfl_scale
integer k, is_f, nx, h5e
logical ferr
character(24) gname, dset_r, dset_i

!

wfl => ffield(ihh)%wf
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
dset_r = 'field-real';  dset_i = 'field-imag'
if (use_y) then
  dset_r = 'field-real-y';  dset_i = 'field-imag-y'
endif

call hdf5_open_file (trim(fname), 'WRITE', p_id, ferr);  if (ferr) stop 1
call hdf5_write_dataset_int  (p_id, 'gridpoints',   [nx],                       ferr)
call hdf5_write_dataset_real (p_id, 'gridsize',     [wfl%dx],                   ferr)
call hdf5_write_dataset_real (p_id, 'refposition',  [wfl%ref_position],         ferr)
call hdf5_write_dataset_real (p_id, 'wavelength',   [wfl%wavelength],           ferr)
call hdf5_write_dataset_int  (p_id, 'slicecount',   [nslice + n_banked(ihh)],   ferr)
call hdf5_write_dataset_real (p_id, 'slicespacing', [fbeam%slice_spacing],      ferr)

allocate (work(nx*nx), re_w(nx*nx), im_w(nx*nx))

do is_f = 1, nslice
  write (gname, '(a, i0.6)') 'slice', is_f
  call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
  if (use_y) then
    work = dfl_scale * reshape(real(wfl%Ey(:,:,is_f), rp), [nx*nx])
  else
    work = dfl_scale * reshape(real(wfl%Ex(:,:,is_f), rp), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) stop 1
  if (use_y) then
    work = dfl_scale * reshape(aimag(wfl%Ey(:,:,is_f)), [nx*nx])
  else
    work = dfl_scale * reshape(aimag(wfl%Ex(:,:,is_f)), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) stop 1
  call H5Gclose_f (g_id, h5e)
enddo

if (n_banked(ihh) > 0) then
  if (.not. allocated(wf1%Ex)) then
    allocate (wf1%Ex(nx, nx, 1))
    wf1%dx = wfl%dx;  wf1%dy = wfl%dy;  wf1%dz = fbeam%slice_spacing
  endif
  wf1%wavelength = wfl%wavelength      ! Per field: banked light drifts at ITS wavelength.

  call hdf5_open_file (trim(escaped_file_name(ihh)), 'READ', e_id, ferr);  if (ferr) stop 1
  do k = 1, n_banked(ihh)
    write (gname, '(a, i0.6)') 'slice', k
    g_id = hdf5_open_group (e_id, trim(gname), ferr, .true.);  if (ferr) stop 1
    call hdf5_read_dataset_real (g_id, trim(dset_r), re_w, ferr, trim(gname));  if (ferr) stop 1
    call hdf5_read_dataset_real (g_id, trim(dset_i), im_w, ferr, trim(gname));  if (ferr) stop 1
    call H5Gclose_f (g_id, h5e)
    wf1%Ex(:,:,1) = reshape(cmplx(re_w, im_w, wf_rp), [nx, nx]) / dfl_scale

    call wavefront_drift (wf1, z_now - bank_z(k, ihh), ferr);  if (ferr) stop 1

    write (gname, '(a, i0.6)') 'slice', nslice + (n_banked(ihh) - k + 1)
    call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
    work = dfl_scale * reshape(real(wf1%Ex(:,:,1), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) stop 1
    work = dfl_scale * reshape(aimag(wf1%Ex(:,:,1)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) stop 1
    call H5Gclose_f (g_id, h5e)
  enddo
  call H5Fclose_f (e_id, h5e)
endif

call H5Fclose_f (p_id, h5e)

end subroutine write_pulse_file

end program fel_track_test
