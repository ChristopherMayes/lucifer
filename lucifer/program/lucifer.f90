!+
! Program lucifer
!
! FEL tracker validated against Genesis 1.3 Version 4. The PHYSICS (coordinates and
! conventions, the FEL step, the field solver, slippage, loading, import, migration,
! collective effects, and each piece's Genesis provenance and measured validation
! level) lives in the manual, lucifer/doc/fel-physics.tex. The measured numbers and
! how to run the checks are in lucifer/README.md. This header documents the inputs.
!
! The program walks a Bmad lattice and applies the seam of the design (manual
! sec:element, sec:seam). FEL segments are real Bmad wiggler/undulator elements with
! tracking_method = custom, their FEL parameters read from lattice attributes (aw from
! b_max/l_period, helicity from field_calc, step from ds_step, asserted by name at
! setup), stepped with the transcribed Genesis physics of fel_track_mod. Every other
! element tracks each slice's bunch with Bmad's track1_bunch (the packed arrays ARE
! Bmad coordinates, so conversion is a plain copy) and drifts the field by
! wavefront_drift, with one advance of the reference phase phi0 per element.
!
! Time dependence follows from the starting state alone: a multi-slice window makes a
! time-dependent run (slippage active), a single slice the steady state, no separate
! switch (the same rule as Genesis). The slippage schedule is precomputed over the
! lattice, transcribing Lattice::calcSlippage (manual sec:slippage, including the
! drift autophasing and its unguarded end-of-lattice fixup), and applied after each
! step's field solve, before its diagnostics -- Gencore's step order. The field record
! rotates rather than moves. Everything reading it in time order goes through
! fel_field_index.
!
! The starting state can be a pair of Genesis dumps (&write of beam and field), so
! both codes track from bitwise-identical initial conditions. It can also be
! self-generated, or an imported distribution (below). Diagnostics matching Genesis's
! definitions are recorded at the same z positions Genesis records them: once at the
! start and once after every integration step, one step per interlude element. Each
! record is one row per slice, in time-window order.
!
! Input is a namelist file of THREE GROUPS, laid out the way Tao lays out its init file
! (manual sec:program). The structs live in fel_struct, the parsing in fel_input_mod,
! both library. &fel_params carries the run: the lattice, the global%... switches,
! Bmad's own bmad_com and space_charge_com set directly (Tao's &tao_params pattern),
! and the wake%/sc% collective descriptions. &fel_beam_init carries the beam:
! Bmad's beam_init%... description, the imp%... resampler, source/output files and the
! beam-side check knobs. &fel_wavefront_init carries the radiation: the
! wavefront_init%... starting condition (the beam_init analog: the field record IS
! the time window, so the window lives here) and field_file imports. The retired flat
! &fel_track_params group is refused by name, each parameter mapped to its new home.
!
!   &fel_params
!     lat_file = "aramis.bmad"                 ! Bmad lattice.
!     global%out_root = "fel_td"               ! Prefix for the output files.
!     global%interlude_model = "bmad"          ! "bmad" (the seam, default) or "genesis".
!   /
!   &fel_beam_init
!     beam_file = "Aramis-initial.par.h5"      ! Genesis particle dump to start from.
!     split_weights = F                        ! Weight-invariance test mode (see below).
!   /
!   &fel_wavefront_init
!     field_file = "Aramis-initial.fld.h5"     ! Genesis field dump to start from.
!   /
!
! The FEL tracking mode and unaveraged parameters are per-element LATTICE attributes
! (registered by this program, usable on any wiggler/undulator, class-settable as
! wiggler::*[attr] = ...). Use NAMED values, defined as one-line lattice variables
! (fel_transcribed = -1, fel_averaged = 0, fel_unaveraged = 1), matching the code's
! fel_transcribed$/fel_averaged$/fel_unaveraged$ parameters (fel_track_mod):
!
!     fel_tracking          ! unset/0 = averaged, the bmad_standard wiggler kernel's
!                           !   transverse maps -- BMAD'S OWN KERNEL IS THE DEFAULT.
!                           ! 1 = the unaveraged verification mode: full Newton-Lorentz
!                           !   quiver, no fc/faw (manual sec:unaveraged). The run
!                           !   writes <out_root>.ledger.txt.
!                           ! -1 = averaged with the transcribed-Genesis transverse
!                           !   maps: VALIDATION-INTERNAL (the Genesis tiers require
!                           !   transcription-level transport, so no production lattice
!                           !   writes it). Differences priced in the README.
!     fel_steps_per_period  ! Unaveraged substeps per period. unset/0 -> 20.
!                           !   Below 10 refused (MINERVA's floor).
!     fel_ramp_periods      ! sin^2 entry/exit ramp length [periods]. unset/0 -> 2.
!                           !   A TRUE hard edge (test configuration) is the explicit
!                           !   sentinel -1: silence never means hard edge.
! The remaining global%... switches of &fel_params:
!
!     global%write_initial = F                 ! Also dump the initial state (Genesis format).
!     global%dump_beam_at = "UND3", "quadrupole::*"  ! Dump the beam (Genesis .par.h5 format)
!                                              !   at the end of the named elements (Bmad
!                                              !   locator syntax, class::name allowed).
!                                              !   An entry matching nothing is refused by name.
!     global%dump_field_at = "UND3"            ! Same for the field (Genesis .fld.h5, unrotated).
!     global%keep_escaped_field = F            ! Bank the field slices slippage transmits out of
!                                              !   the window (<out_root>-escaped.fld.h5, with
!                                              !   wavefront_params and z_transmit per slice) and
!                                              !   reconstruct the FULL PULSE at the exit plane
!                                              !   (<out_root>-pulse.fld.h5) by free-space
!                                              !   propagation at finalize. Manual sec:stats.
!     global%migrate = F                       ! Slice migration (see below).
!     global%migrate_check = F                 ! Verify phase continuity at each migration.
!     global%ran_seed = 12345                  ! The one RNG seed (generation, import, noise).
!     global%load_only = F                     ! Build the initial state, dump it, exit
!                                              !   without tracking (shot-noise checks).
!     global%reference_run = F                 ! NO FEL interaction: Bmad tracks everything.
!
! global%migrate = T moves particles between slices when their ponderomotive phase leaves the
! slice window (fel_migrate_slices, manual sec:migration), called serially after every
! element. OFF BY DEFAULT, deliberately: the Genesis-comparison tiers run against
! Genesis WITHOUT one4one, which never migrates. Migration would be a physics-model
! difference inside a transcription-level comparison. Dropped charge is counted and
! reported. migrate_check = T also verifies exact phase continuity at every
! migration and reports the worst deviation.
!
! Alternatively, leave beam_file and field_file blank and the program generates its own
! starting condition (a quiet-start beam and a Gaussian seed field), making a
! self-contained single run with no Genesis anywhere (see lucifer/examples). The BEAM
! is described by Bmad's standard beam_init_struct, the same block the import path
! uses: one bulk-bunch description, two generation methods. The import resamples
! real particles from it, the quiet-start loader evaluates it analytically per slice.
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
! be worse than a custom one. The radiation starting condition is
! &fel_wavefront_init's wavefront_init%... (the beam_init analog), with the Genesis4
! &field / &time names each field maps to:
!
!     wavefront_init%lambda0 = 1e-10          ! Radiation wavelength [m], REQUIRED (with
!                              !   dumps it comes from the file). Deliberately not
!                              !   defaulted from the lattice resonance: the first
!                              !   undulator may be off.
!     wavefront_init%seed_power = 5e3         ! Seed power [W] (&field power). 0 = dark start.
!     wavefront_init%seed_waist_size = 30e-6  ! Seed 1/e^2 intensity radius w0 [m]
!                              !   (&field waist_size).
!     wavefront_init%seed_polarization = 'x'  ! 'x' or 'y'.
!     wavefront_init%grid_n_pts = 255         ! Transverse grid points per side (&field ngrid).
!     wavefront_init%grid_half_width = 2e-4   ! Grid half width [m] (&field dgrid).
!     wavefront_init%window_length = 0        ! Time window [m] (&time slen). 0: derived from
!                              !   the bunch (+-4 sig_z Gaussian, the grid extent flat,
!                              !   one slice steady state). Override it for slippage
!                              !   headroom. A window that clips the bunch warns.
!     wavefront_init%window_sample = 1        ! Slice spacing / lambda0 (integer, &time sample).
!     wavefront_init%harmonics = 1, 3         ! The field set (manual sec:field-set).
!
! The beam-side generation knobs of &fel_beam_init:
!
!     nbins = 8                ! Beamlet size of the quiet start (quiet below nbins).
!     shotnoise = F            ! Impose physical shot noise (&time shotnoise).
!                              !   Time-dependent windows only, Genesis's dotime rule.
!     gen_test_weights = F     ! Validation knob: alternate beamlet weights 0.25x/1.75x
!                              !   (charge preserving, uniform within each beamlet) to
!                              !   exercise the weighted-noise paths. Not physics input.
!
! Element wakes: elements carrying Bmad sr_wake definitions (pseudomodes or a
! tabular z_long) act across the WHOLE time window via slice concatenation (manual
! sec:seamwake: conventions, ds_wake, the mid-element wiggler kick, refusals).
! write_wake_kernels = "<file>" exports the transcribed wake kernels for building
! matching z_long tables (see the README's seam-wake section and examples/bmad_wake).
!
! Third way in: import a particle DISTRIBUTION, an arbitrary bunch resampled into
! slices by the transcribed Genesis importdistribution method (fel_import_mod, manual
! sec:import). The bunch comes from the SAME beam_init block (use_beam_init = T:
! init_beam_distribution generates it, honoring everything Bmad honors) or from an
! openPMD-beamphysics file (dist_file). Note: the honored-fields contract above
! applies to the quiet-start generator only. window_sample, ran_seed and the seed
! field are shared with the generator: one seed governs generation, resampling and
! noise. Knobs, named after &importdistribution's where one exists:
!
!     imp%slicewidth = 0.01    ! Sampling window / bunch length (Genesis's slicewidth).
!     imp%npart = 8192         ! Macroparticles per slice after resampling.
!     imp%nbins = 4            ! Beamlet size (quiet-load bins) of the resample.
!     imp%nslice = 0           ! 0: round(bunch_length/spacing), Genesis's rule.
!                              ! (Genesis's match/center are NOT ported: a Bmad lattice
!                              !   carries its Twiss and beam_init generates matched
!                              !   bunches already. Also see fel_import_mod's header.)
!     use_beam_init = F        ! Generate the bunch from the beam_init%... block.
!     dist_file = ""           ! Or read an openPMD-beamphysics file.
!     write_dist_file = ""     ! Write the bunch as a Genesis &importdistribution
!                              !   input (t/p/x/xp/y/yp + charge, t = -tau/c), the
!                              !   shared file of the cross-code checks.
!     write_opmd_file = ""     ! Write the bunch as openPMD-beamphysics.
!     imp_split_weights = F    ! Check knob: coincident w/3 + 2w/3 copies before import.
!
! The quiet start and the weighted Fawley shot noise (shotnoise = T), including
! the noise-level algebra and the N_eff refusal guard, are the manual's sec:loading.
! The loader warns where Genesis silently clamps beamlets with fewer real electrons
! than macroparticles.
!
! interlude_model selects how the field-free elements are handled. "bmad" is the seam
! (track1_bunch, exact theta mapping, wavefront_drift). "genesis" uses the transcribed
! Genesis interlude step everywhere, which prices what the seam changes (the manual's
! sec:interlude and sec:seam). The slippage schedule is identical in both models.
!
! split_weights = T replaces each imported particle by two copies at identical
! coordinates carrying 1/3 and 2/3 of its weight. Every collective observable must be
! identical to the unsplit run. This checks the weighted paths, which no Genesis
! comparison can (Genesis dumps carry no weights).
!
! Outputs: <out_root>.diag.txt, ONLY with global%write_diag = T (one row per slice per
! record: z, slice, field and beam diagnostics). <out_root>.stats.h5 is the production
! statistics file (manual sec:stats): per-record per-slice beam moments named as
! bunch_params_struct components, per-record per-slice wavefront_params, and the
! evaluated calc_bunch_params at element ends, in fixed Bmad units.
! <out_root>-final.fld.h5 and <out_root>-final.par.h5 are Genesis-format dumps of the
! end state, for field-by-field comparison. The field dump is unrotated to time order
! first, as writeFieldHDF5 does.
!-

program lucifer

use fel_struct
use fel_input_mod
use fel_setup_mod
use fel_init_mod
use fel_track_line_mod
use fel_io_mod

implicit none

! The driver is read-parse-call (manual sec:program): the namelist layer fills the
! input structs, the library builds and walks the run, and every library error
! RETURNS here. This is the one place that stops. The check instruments
! (split_weights, swap_beam_xy, gen_test_weights, imp_split_weights,
! write_wake_kernels) ride along in the input structs and act inside init/setup.

type (fel_run_struct), target :: run
integer n_arg, iu_k, i_k
logical err
character(400) param_file
character(*), parameter :: r_name = 'lucifer'

!

n_arg = command_argument_count()
if (n_arg /= 1) then
  print '(a)', 'Usage: lucifer <param_file>'
  stop 1
endif
call get_command_argument (1, param_file)

call fel_read_input (param_file, run, err)
if (err) stop 1

call fel_setup_lattice (run, err)
if (err) stop 1

call fel_init_beam (run, err)
if (err) stop 1

call fel_init_wavefront (run, err)
if (err) stop 1

call fel_setup_schedule (run, err)
if (err) stop 1

call fel_write_header (run)

! Check instrument: export the transcribed single-particle wake kernels for
! cross-validation against Genesis (manual sec:wakes). The kernels are built by
! fel_setup_schedule's fel_wake_init. Note: The s = 0 entries carry the Bane
! self-slice half factor.

if (run%wake_init%write_kernels /= '' .and. run%coll%wake%on) then
  open (newunit = iu_k, file = trim(run%wake_init%write_kernels), action = 'write')
  write (iu_k, '(a)') '# s [m]   wakeres   wakegeo   wakerou   [eV/(m electron)]; s=0 rows are HALVED (Bane self-slice)'
  do i_k = 1, run%coll%wake%ns
    write (iu_k, '(4es24.15e3)') run%coll%wake%ds * (i_k-1), run%coll%wake%wakeres(i_k), &
                                 run%coll%wake%wakegeo(i_k), run%coll%wake%wakerou(i_k)
  enddo
  close (iu_k)
  call out_io (s_info$, r_name, 'Wrote wake kernels: ' // trim(run%wake_init%write_kernels))
endif

if (run%global%write_initial .or. run%global%load_only) then
  call wavefront_write_genesis4 (run%ffield(1)%wf, trim(run%global%out_root) // '-initial.fld.h5', err, 'x')
  if (err) stop 1
  call fel_write_genesis4_beam (run%fbeam, trim(run%global%out_root) // '-initial.par.h5', err)
  if (err) stop 1
endif

if (run%global%load_only) then
  call out_io (s_info$, r_name, 'load_only set; initial state written, no tracking.')
  stop 0
endif

call track_fel_line (run, err)
if (err) stop 1

! Final dumps. The field records are unrotated to time order first (position is of
! the time window holds record slice 1 + mod(is-1+first, nslice)), which is what
! Genesis's field writer does on the fly (manual sec:slippage).

call fel_write_genesis4_beam (run%fbeam, trim(run%global%out_root) // '-final.par.h5', err)
if (err) stop 1

call fel_dump_field_set (run, trim(run%global%out_root) // '-final', err)
if (err) stop 1

call fel_finalize_diagnostics (run, err)
if (err) stop 1

call fel_write_footer (run)

call wavefront_fft_free()

end program lucifer
