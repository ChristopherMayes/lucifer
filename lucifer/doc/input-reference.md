# Input reference

The normative reference for every namelist parameter Lucifer honors: its default, its meaning, and what refuses it. The program reads three namelist groups from one input file, each setting Bmad or program structures directly, which is Tao's `&tao_params` pattern. Defaults live in the struct declarations, so this document and the code agree by construction.

How to build and run is [`user-guide.md`](user-guide.md). What the outputs hold is [`reading-output.md`](reading-output.md). The measured levels are [`validation.md`](validation.md). The physics is the manual, [`fel-physics.tex`](fel-physics.tex).

Every resolved input is also written into the stats file's `params/` group, one subgroup per structure, so a finished run states its own configuration with every default made explicit.

The program walks a Bmad lattice and applies the seam of the design (manual: the FEL-element and seam sections). FEL segments are real Bmad wiggler/undulator elements with
tracking_method = custom, their FEL parameters read from lattice attributes (aw from
b_max/l_period, helicity from field_calc, step from ds_step, asserted by name at
setup), stepped with the transcribed Genesis physics of fel_track_mod. Every other
element tracks each slice's bunch with Bmad's track1_bunch (the packed arrays ARE
Bmad coordinates, so conversion is a plain copy) and drifts the field by
wavefront_drift, with one advance of the reference phase phi0 per element.

Time dependence follows from the starting state alone: a multi-slice window makes a
time-dependent run (slippage active), a single slice the steady state, no separate
switch (the same rule as Genesis). The slippage schedule is precomputed over the
lattice, transcribing Lattice::calcSlippage (manual: the slippage section, including the
drift autophasing and its unguarded end-of-lattice fixup), and applied after each
step's field solve, before its diagnostics -- Gencore's step order. The field record
rotates rather than moves. Everything reading it in time order goes through
fel_field_index.

The starting state can be a pair of Genesis dumps (&write of beam and field), so
both codes track from bitwise-identical initial conditions. It can also be
self-generated, or an imported distribution (below). Diagnostics matching Genesis's
definitions are recorded at the same z positions Genesis records them: once at the
start and once after every integration step, one step per interlude element. Each
record is one row per slice, in time-window order.

Input is a namelist file of three groups, laid out the way Tao lays out its init file
(doc/user-guide.md). The structs live in fel_struct, the parsing in fel_input_mod,
both library. &fel_params carries the run: the lattice, the global%... switches,
Bmad's own bmad_com and space_charge_com set directly (Tao's &tao_params pattern),
and the wake%/sc% collective descriptions. &fel_beam_init carries the beam:
Bmad's beam_init%... description, the imp%... resampler, source/output files and the
beam-side check knobs. &fel_wavefront_init carries the radiation: the
wavefront_init%... starting condition (the beam_init analog: the field record IS
the time window, so the window lives here) and field_file imports. The retired flat
&fel_track_params group is refused by name, each parameter mapped to its new home.

```
  &fel_params
    lat_file = "aramis.bmad"                 ! Bmad lattice.
    global%out_root = "fel_td"               ! Prefix for the output files.
    global%interlude_model = "bmad"          ! "bmad" (the seam, default) or "genesis".
  /
  &fel_beam_init
    beam_file = "Aramis-initial.beam.h5"     ! openPMD particle dump to start from. A
                                             !   file that is not openPMD is refused by
                                             !   name, with the conversion command:
                                             !   tests/scripts/convert_genesis.py.
    nbins = 8                                ! Beamlet size. No dump format carries it.
    split_weights = F                        ! Weight-invariance test mode (see below).
  /
  &fel_wavefront_init
    field_file = "Aramis-initial.wf.h5"      ! openPMD EXT_Wavefront field dump.
    wavefront_init%lambda0 = 1e-10           ! Required with a beam dump: the file
                                             !   carries the slice partition and not the
                                             !   radiation it was sliced on.
    wavefront_init%window_sample = 3         ! Slice spacing in wavelengths, likewise.
  /
```

The FEL tracking mode and unaveraged parameters are per-element LATTICE attributes
(registered by this program, usable on any wiggler/undulator, class-settable as
wiggler::*[attr] = ...). Use NAMED values, defined as one-line lattice variables
(fel_transcribed = -1, fel_averaged = 0, fel_unaveraged = 1), matching the code's
fel_transcribed$/fel_averaged$/fel_unaveraged$ parameters (fel_track_mod):

```
    fel_tracking          ! unset/0 = averaged, the bmad_standard wiggler kernel's
                          !   transverse maps -- BMAD'S own kernel is the default.
                          ! 1 = the unaveraged mode: full Newton-Lorentz
                          !   quiver, no fc/faw (manual: the unaveraged-mode section). The run
                          !   writes <out_root>.ledger.txt.
                          ! -1 = averaged with the transcribed-Genesis transverse
                          !   maps: VALIDATION-INTERNAL (the Genesis tiers require
                          !   transcription-level transport, so no production lattice
                          !   writes it). Differences priced in the README.
    fel_steps_per_period  ! Unaveraged substeps per period. unset/0 -> 20.
                          !   Below 10 refused (MINERVA's floor).
    fel_ramp_periods      ! sin^2 entry/exit ramp length [periods]. unset/0 -> 2.
                          !   A TRUE hard edge (test configuration) is the explicit
                          !   sentinel -1: silence never means hard edge.
```

The remaining global%... switches of &fel_params:

```
    global%write_initial = F                 ! Also dump the initial state.
    global%dump_beam_at = "UND3", "quadrupole::*"  ! Dump the beam (.beam.h5) at the end
                                             !   of the named elements (Bmad locator
                                             !   syntax, class::name allowed). An entry
                                             !   matching nothing is refused by name.
    global%dump_field_at = "UND3"            ! Same for the field (.wf.h5, unrotated).
    global%keep_escaped_field = F            ! Bank the field slices slippage transmits out of
                                             !   the window (<out_root>-escaped.fld.h5, with
                                             !   wavefront_params and z_transmit per slice) and
                                             !   reconstruct the full pulse at the exit plane
                                             !   (<out_root>-pulse.fld.h5) by free-space
                                             !   propagation at finalize. Manual: the diagnostic-output section.
    global%migrate = F                       ! Slice migration (see below).
    global%migrate_check = F                 ! Verify phase continuity at each migration.
    global%ran_seed = 12345                  ! The one RNG seed (generation, import, noise).
    global%load_only = F                     ! Build the initial state, dump it, exit
                                             !   without tracking (shot-noise checks).
    global%reference_run = F                 ! No FEL interaction: Bmad tracks everything.
    global%record_environment = F            ! Also record the user name and working
                                             !   directory in stats.h5's meta/ group. OFF
                                             !   by default: a stats file is meant to
                                             !   travel, and those identify a person and a
                                             !   machine. Manual: the provenance section.
```

global%migrate = T moves particles between slices when their ponderomotive phase leaves the
slice window (fel_migrate_slices, manual: the slice-migration section), called serially after every
element. off by default, deliberately: the Genesis-comparison tiers run against
Genesis WITHOUT one4one, which never migrates. Migration would be a physics-model
difference inside a transcription-level comparison. Dropped charge is counted and
reported. migrate_check = T also verifies exact phase continuity at every
migration and reports the worst deviation.

Alternatively, leave beam_file and field_file blank and the program generates its own
starting condition (a quiet-start beam and a Gaussian seed field), making a
self-contained single run with no Genesis anywhere (see lucifer/examples). The BEAM
is described by Bmad's standard beam_init_struct, the same block the import path
uses: one bulk-bunch description, two generation methods. The import resamples
real particles from it, the quiet-start loader evaluates it analytically per slice.
The Twiss is always the lattice's. Honored beam_init fields:

```
    beam_init%n_particle    ! Macroparticles per slice; a positive multiple of nbins.
                            !   (The import path reads it as bunch particles instead.)
    beam_init%a_norm_emit   ! Normalized emittances [m rad] (a_emit/b_emit refused:
    beam_init%b_norm_emit   !   normalized only, Bmad's preferred form).
    beam_init%sig_pz        ! Fractional momentum spread dP/P0. The Gaussian gamma
                            !   spread is delta_gamma = beta0*p0_mc*sig_pz. (sig_e is
                            !   deprecated Bmad-wide and does not exist here.)
    beam_init%bunch_charge  ! Charge [C]. The current is DERIVED, never input.
    beam_init%sig_z         ! Bunch length [m], with distribution_type(3):
    beam_init%distribution_type(3)  ! "RAN_GAUSS" (default): Gaussian current profile
                            !   I(s) = Q*c/(sqrt(2pi)*sig_z) * exp(-s^2/2 sig_z^2) at
                            !   the slice centers, bunch centered in the window.
                            !   sig_z = 0 is the steady state: one slice holding the
                            !   whole charge, I = Q*c/slice_spacing; refused by name
                            !   for time-dependent windows.
                            ! "GRID" (Bmad's uniform): flat I = Q*c/(x_max - x_min)
                            !   over grid(3)%x_min..x_max (the z extent).
```

EVERY other beam_init field that is set is refused by name (see
check_beam_init_contract) -- a standard structure that silently drops fields would
be worse than a custom one. The radiation starting condition is
&fel_wavefront_init's wavefront_init%... (the beam_init analog), with the Genesis4
&field / &time names each field maps to:

```
    wavefront_init%lambda0 = 1e-10          ! Radiation wavelength [m], REQUIRED (with
                             !   dumps it comes from the file). Deliberately not
                             !   defaulted from the lattice resonance: the first
                             !   undulator may be off.
    wavefront_init%seed_power = 5e3         ! Seed power [W] (&field power). 0 = dark start.
    wavefront_init%seed_waist_size = 30e-6  ! Seed 1/e^2 intensity radius w0 [m]
                             !   (&field waist_size).
    wavefront_init%seed_polarization = 'x'  ! 'x' or 'y'.
    wavefront_init%grid_n_pts = 255         ! Transverse grid points per side (&field ngrid).
    wavefront_init%grid_half_width = 2e-4   ! Grid half width [m] (&field dgrid).
    wavefront_init%window_length = 0        ! Time window [m] (&time slen). 0: derived from
                             !   the bunch (+-4 sig_z Gaussian, the grid extent flat,
                             !   one slice steady state). Override it for slippage
                             !   headroom. A window that clips the bunch warns.
    wavefront_init%window_sample = 1        ! Slice spacing / lambda0 (integer, &time sample).
    wavefront_init%harmonics = 1, 3         ! The field set (manual: harmonic radiation).
```

The beam-side generation knobs of &fel_beam_init:

```
    nbins = 8                ! Beamlet size of the quiet start (quiet below nbins).
    shotnoise = F            ! Impose physical shot noise (&time shotnoise).
                             !   Time-dependent windows only, Genesis's dotime rule.
    gen_test_weights = F     ! Validation knob: alternate beamlet weights 0.25x/1.75x
                             !   (charge preserving, uniform within each beamlet) to
                             !   exercise the weighted-noise paths. Not physics input.
```

Element wakes: elements carrying Bmad sr_wake definitions (pseudomodes or a
tabular z_long) act across the WHOLE time window via slice concatenation (manual: the
Bmad-element-wakes section, for conventions, ds_wake, the mid-element wiggler kick, refusals).
write_wake_kernels = "<file>" exports the transcribed wake kernels for building
matching z_long tables (see the README's seam-wake section and examples/bmad_wake).

Third way in: import a particle DISTRIBUTION, an arbitrary bunch resampled into
slices by the transcribed Genesis importdistribution method (fel_import_mod, manual: the
distribution-import section). The bunch comes from the SAME beam_init block (use_beam_init = T:
init_beam_distribution generates it, honoring everything Bmad honors) or from an
openPMD-beamphysics file (dist_file). Note: the honored-fields contract above
applies to the quiet-start generator only. window_sample, ran_seed and the seed
field are shared with the generator: one seed governs generation, resampling and
noise. Knobs, named after &importdistribution's where one exists:

```
    imp%slicewidth = 0.01    ! Sampling window / bunch length (Genesis's slicewidth).
    imp%npart = 8192         ! Macroparticles per slice after resampling.
    imp%nbins = 4            ! Beamlet size (quiet-load bins) of the resample.
    imp%nslice = 0           ! 0: round(bunch_length/spacing), Genesis's rule.
                             ! (Genesis's match/center are NOT ported: a Bmad lattice
                             !   carries its Twiss and beam_init generates matched
                             !   bunches already. Also see fel_import_mod's header.)
    use_beam_init = F        ! Generate the bunch from the beam_init%... block.
    dist_file = ""           ! Or read an openPMD-beamphysics file.
    write_dist_file = ""     ! Write the bunch as a Genesis &importdistribution
                             !   input (t/p/x/xp/y/yp + charge, t = -tau/c), the
                             !   shared file of the cross-code checks.
    write_opmd_file = ""     ! Write the bunch as openPMD-beamphysics.
    imp_split_weights = F    ! Check knob: coincident w/3 + 2w/3 copies before import.
```

The quiet start and the weighted Fawley shot noise (shotnoise = T), including
the noise-level algebra and the N_eff refusal guard, are the manual's loading section.
The loader warns where Genesis silently clamps beamlets with fewer real electrons
than macroparticles.

interlude_model selects how the field-free elements are handled. "bmad" is the seam
(track1_bunch, exact theta mapping, wavefront_drift). "genesis" uses the transcribed
Genesis interlude step everywhere, which prices what the seam changes (the manual's
interlude and seam sections). The slippage schedule is identical in both models.

split_weights = T replaces each imported particle by two copies at identical
coordinates carrying 1/3 and 2/3 of its weight. Every collective observable must be
identical to the unsplit run. This checks the weighted paths, which no Genesis
comparison can (Genesis dumps carry no weights).

Outputs: <out_root>.diag.txt, ONLY with global%write_diag = T (one row per slice per
record: z, slice, field and beam diagnostics). <out_root>.stats.h5 is the production
statistics file (manual: the diagnostic-output section): per-record per-slice beam moments named as
bunch_params_struct components, per-record per-slice wavefront_params, and the
evaluated calc_bunch_params at element ends, in fixed Bmad units.
The end state is dumped as openPMD, the one format this program speaks:
<out_root>-final.beam.h5 and -final.wf.h5 are openPMD, <out_root>-final.par.h5 and
-final.fld.h5 are Genesis format, for field-by-field comparison against Genesis. The
field dump is unrotated to time order first, as writeFieldHDF5 does.
