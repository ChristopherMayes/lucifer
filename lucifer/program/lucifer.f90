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
! wiggler::*[attr] = ...):
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
! diagnostics), <out_root>-final.fld.h5 and <out_root>-final.par.h5 (Genesis-format dumps
! of the end state, for field-by-field comparison; the field dump is unrotated to time
! order first, as writeFieldHDF5 does).
!-

program fel_track_test

use fel_track_mod
use fel_unaveraged_mod
use fel_import_mod
use wavefront_hdf5_mod
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
type (wavefront_struct) wf
type (bunch_struct) bunch
type (bunch_struct) wake_bunch          ! Whole-window bunch for element sr wakes.
real(rp), allocatable :: wake_beta0(:)  ! Entry betas for the exact concat/split inverse.
type (fel_und_struct) und
type (fel_slip_struct) slip
type (fel_slice_diag_struct) bdiag

real(rp) :: gamma0 = 0                  ! DERIVED from the lattice e_tot.
logical :: split_weights = .false.
logical :: write_initial = .false.
logical :: migrate = .false., migrate_check = .false.
character(400) :: lat_file = '', beam_file = '', field_file = '', out_root = 'fel_track'
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
real(rp) dE_step
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
integer ie, is, istep, n_arg, iu_diag, iu_nml, nslice, prev_ie, n_moved_tot, iu_wake
logical err, timerun

type (fel_collective_struct) coll

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
                           write_opmd_file, imp_split_weights, write_wake_kernels

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

! ONE reference energy, and the lattice is it: gamma0 = e_tot/m_e c^2 from the lattice
! header, never a namelist input. There used to be a namelist gamma0 for Genesis-deck
! symmetry; the first external user fed it a hand-rounded value against a round lattice
! e_tot, the two disagreed at 1.4e-9, and the run died mid-tracking on the seam's
! backstop p0c check with raw numbers -- the FEL physics ran on one reference while
! Bmad's momenta were normalized by the other. Two specifications of one truth is the
! defect; the redundant one was removed (the deliverable-9 rule: parameters live on
! the lattice).

if ((beam_file == '') .neqv. (field_file == '')) then
  print '(a)', 'fel_track_test: give both beam_file and field_file, or neither (to generate).'
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
  call wavefront_read_genesis4 (wf, field_file, err)
  if (err) stop 1
elseif (dist_file /= '' .or. use_beam_init) then
  call import_initial_state ()
else
  call generate_initial_state ()
endif

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
slip%timerun = timerun
slip%sample = fbeam%slice_spacing / fbeam%wavelength
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

call setup_fel_elements ()
any_unavg = any(fel_mode == 1 .and. is_fel)

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
  write (iu_ledger, '(a)') '#         z          E_beam_rel [J]          U_field [J]           dE_kick [J]'
endif

n_moved_tot = 0
charge_dropped_tot = 0
b_dev_max = 0

z_now = 0
call write_diag_rows()     ! Initial record, matching Genesis's diag before the first step.

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
    und%bmad_transport = (fel_mode(ie) == 0)     ! The default: bmad_standard's kernel.
    und_slip_step = (1 + und%aw**2) / (2 * gamma0_ref**2 * wf%wavelength)
    und%nstep = max(1, nint(ele%value(num_steps$)))
    und%dz = ele%value(l$) / und%nstep

    if (fel_mode(ie) == 1) then
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
      if (fel_mode(ie) == 1) then
        call fel_unavg_step (und, ustate, fbeam, wf, slip, und%dz, istep == 1, &
                             istep == und%nstep, dE_step, err)
      else
        call fel_track_und_step (und, fbeam, wf, slip, coll, err)
      endif
      if (err) stop 1

      ! Element sr wake, Bmad's once-per-passage convention mirrored: one kick at the
      ! step nearest mid-element, scaled to the full element length (scale_with_length
      ! uses ele's l), applied across the WHOLE window (deliverable 11). Direct kick,
      ! no transport -- this walk owns transport inside wigglers.

      if (associated(wake_src) .and. istep == (und%nstep + 1)/2) call apply_bmad_wake_kick (wake_src)

      if (istep == und%nstep) then
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step + ele_slip(ie))
      else
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step)
      endif
      z_now = z_now + und%dz
      if (fel_mode(ie) == 1) call write_ledger_row ()
      if (istep == und%nstep) call do_migrate ()
      call write_diag_rows()
    enddo

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

    call wavefront_drift (wf, ele%value(l$), err)
    if (err) stop 1

    ! The chamber does not end where the undulator does: the wake's energy loss applies
    ! through seam interludes too, as one kick of the element's length (Genesis applies
    ! it every step, and an interlude is one step).

    call fel_wake_apply (coll%wake, fbeam, ele%value(l$))

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call do_migrate ()
    call write_diag_rows()

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, wf, slip, coll, err)
    if (associated(wake_src)) call apply_bmad_wake_kick (wake_src)
    if (err) stop 1

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call do_migrate ()
    call write_diag_rows()
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

! Final dumps in Genesis format. The field record is unrotated to time order first --
! time window position is holds record slice 1 + mod(is-1+first, nslice) -- which is
! what Genesis's field writer does on the fly (manual sec:slippage).

if (slip%first /= 0) then
  wf%Ex = cshift(wf%Ex, shift = slip%first, dim = 3)
  slip%first = 0
endif

call wavefront_write_genesis4 (wf, trim(out_root) // '-final.fld.h5', err, 'x')
if (err) stop 1

call fel_write_genesis4_beam (fbeam, trim(out_root) // '-final.par.h5', err)
if (err) stop 1

print '(a)', 'fel_track_test done.'
print '(a)', '  ' // trim(out_root) // '.diag.txt'
print '(a)', '  ' // trim(out_root) // '-final.fld.h5'
print '(a)', '  ' // trim(out_root) // '-final.par.h5'

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

if (seed_power > 0) then
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
  if (abs(rv - fel_mode(je)) > 1e-9_rp .or. fel_mode(je) < -1 .or. fel_mode(je) > 1) then
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
enddo

if (.not. any(is_fel)) then
  print '(a)', 'fel_track_test: the lattice has no FEL elements (wiggler/undulator with tracking_method = custom).'
  stop 1
endif

! Element sr wakes act across the WHOLE window (deliverable 11): all slices
! concatenate into one bunch in global window coordinates and Bmad's wake machinery
! applies unmodified. What is checked here, by name: lr (multi-bunch) wakes are not
! supported; a pseudomode wake whose z_max is shorter than the window would have Bmad
! kill the bunch mid-run; a z_long table narrower than the window would overflow its
! binning grid the same way.

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

end subroutine setup_fel_elements

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
! fel_field_index(slip, is, nslice), the unrotation of manual sec:slippage.

real(rp) power, on_axis
integer is

do is = 1, nslice
  call fel_field_diag (wf, fel_field_index(slip, is, nslice), power, on_axis)
  call fel_slice_diag (fbeam, fbeam%slice(is), ks, bdiag)

  write (iu_diag, '(es24.16, i8, 10es24.16)') z_now, is, power, on_axis, bdiag%bunching, &
        bdiag%bunching_phase, bdiag%mean_energy, bdiag%sigma_energy, bdiag%sigma_x, bdiag%sigma_y, &
        bdiag%current, bdiag%n_eff
enddo

end subroutine write_diag_rows

!------------------------------------------------------------------------------
subroutine write_ledger_row ()

! One energy-ledger row (unaveraged mode): beam energy RELATIVE to the reference,
! sum(w*(gamma-gamma0))*me [C*eV = J] -- relative so the per-record change is not
! differenced off a large baseline at its own summation-rounding floor (the
! FINDINGS 4.8 lesson) -- total window field energy sum(P_is)*slice_spacing/c [J],
! and the kick-side change dE_step returned by fel_unavg_step.

real(rp) e_beam, u_field, power, on_axis, p_mc_l, g0_l
integer is, ip

e_beam = 0
u_field = 0
g0_l = fel_gamma0(fbeam)
do is = 1, nslice
  do ip = 1, fbeam%slice(is)%n
    p_mc_l = fel_p0_mc(fbeam) * (1 + fbeam%slice(is)%pz(ip))
    e_beam = e_beam + fbeam%slice(is)%weight(ip) * (sqrt(p_mc_l**2 + 1) - g0_l) * m_electron
  enddo
  call fel_field_diag (wf, fel_field_index(slip, is, nslice), power, on_axis)
  u_field = u_field + power * fbeam%slice_spacing / c_light
enddo

write (iu_ledger, '(4es24.16)') z_now, e_beam, u_field, dE_step

end subroutine write_ledger_row

end program fel_track_test
