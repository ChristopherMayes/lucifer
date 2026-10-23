# User guide: building, describing a run, and running it

This page holds everything needed to get from a checkout to a tracked FEL run. The parameter-by-parameter reference is [`input-reference.md`](input-reference.md). What the output files hold and how to read them is [`reading-output.md`](reading-output.md). The measured levels and how to reproduce them are [`validation.md`](validation.md). The physics is the manual, [`fel-physics.md`](fel-physics.md).

## Building

Lucifer is a Bmad user program, built against a Bmad distribution rather than as part of one. Its repository holds Lucifer alone, and a Bmad checkout sits beside it, named by `BMAD_DIR`:

```
git clone https://github.com/ChristopherMayes/lucifer.git
git clone https://github.com/bmad-sim/bmad-ecosystem.git bmad
git -C bmad checkout "$(cat lucifer/lucifer/BMAD_RELEASE)"
cd lucifer
BMAD_DIR=../bmad BUILD_PRODUCTION=N lucifer/util/conda_compile   # debug: Bmad, then Lucifer
```

The Bmad it needs is a `bmad-sim/bmad-ecosystem` release at or after the one named in `lucifer/BMAD_RELEASE`, currently `20261018-0`. That is the first release carrying the two Bmad changes Lucifer needs: the `fel_method` element attribute and openPMD particle patches in the beam files. An older Bmad refuses the lattices and reads a multi-slice beam file as one bunch. `lucifer/util/bmad_release_check <bmad root>` compares a checkout's `bmad_version_date` with that minimum, the date and the revision as numbers. The build script runs it before it builds anything, and an older checkout is refused. The tracker prints the Bmad version it was built against at startup and records it in the stats file.

Its `CMakeLists.txt` includes Bmad's `Master.cmake`, which compiles it against the distribution's libraries and modules. The binaries land in `debug/` and `production/` beside the Lucifer source, which is the repository root. One script builds Bmad and then Lucifer, with the conda preferences the recorded levels were measured under:

```
conda env create -n bmad-build -f <bmad>/.github/bmad-build-env.yaml   # once
BUILD_PRODUCTION=N lucifer/util/conda_compile      # debug, the build to develop against
lucifer/util/conda_compile                         # production
BUILD_BMAD=N lucifer/util/conda_compile            # Lucifer alone, Bmad already built
```

The Bmad source is `BMAD_DIR`, which a Lucifer clone needs since its parent holds no Bmad checkout. A `BMAD_DIR` naming a directory with no `util/dist_source_me` is passed over with a note, since `dist_source_me` uses that name as a loop variable and leaves it at the checkout's `util/` once sourced. By hand, the same build is `source <bmad>/util/dist_source_me` from the Bmad root and then `mk` or `mkd` in the Lucifer directory.

Running the validation harness needs one thing beyond the binary: a Python environment with numpy, h5py, pytest and openPMD-beamphysics 0.16.3 or later, described by `lucifer/wavefront/tests/environment.yml`. The scripts use whatever version that environment holds and print it. A development version of openPMD-beamphysics is installed into the environment as well. The harness takes four locations, each defaulting to the layout above: the Lucifer source, which is the scripts' own tree, the Bmad source with its regression suite (`BMAD_DIR`), and the two binary directories (`LUCIFER_DEBUG_BIN`, `LUCIFER_PRODUCTION_BIN`). A relative value is taken against the directory a script was started in. The benchmark records the build its `--build` option declares, debug or production, and never reads it off a directory name. The exact commands the harness runs are in [`validation.md`](validation.md).

## The two FEL methods

The FEL method is Bmad's own `fel_method` attribute, set on the element as any attribute is, and averaged and unaveraged elements mix freely in one line.

| `tracking_method` | what it is for |
|---|---|
| `averaged` | The wiggle-averaged (KMR) model on Bmad's own `bmad_standard` kernel maps. The production method. |
| `unaveraged` | Direct RK4 integration through the analytic undulator field: no averaging, no resonance approximation, fc and JJ nowhere in its inputs. A production method whose ~30x cost per step buys full quiver dynamics, energy accounting the beam actually pays, polarization-agnostic coupling, and arbitrary harmonic content in the current. It is also an independent check on the averaged method, since the two share no approximation. |

A lattice writes them as it writes any attribute, `fel_method = unaveraged`, and a class statement sets them too, `wiggler::*[FEL_METHOD] = unaveraged`. The averaged method's maps have one more option behind them, `global%transport_model`, described in [`input-reference.md`](input-reference.md).

## Running

A run is one input file, which names the lattice:

```
lucifer lucifer.in
```

The name is an argument and any name works. `lucifer.in` is the convention the examples follow, so that a directory reads the same way twice, and it is not a default the program looks for. Running `lucifer` with no argument prints that much, the Bmad version it was built against and the device backend this build carries.

The file vocabulary of a run has four extensions:

| Extension | Written by | What it is |
|---|---|---|
| `.in` | you | The input deck: `&fel_params`, `&fel_beam_init` and `&fel_wavefront_init`, naming one lattice |
| `.bmad` | you | A Bmad lattice, named by the deck's `lat_file`. Variants call a shared one and override it |
| `.h5` | the run | HDF5, openPMD throughout: the statistics file, and the beam and field dumps |
| `.txt` | the run | The named streams, one per subject: diag, ledger, import, fp32, migration, wake |

What each written file holds is the output section below, and what to do with it is [`reading-output.md`](reading-output.md).

The examples are the fastest way in. Each is a directory of real input files, a README with the numbers measured on it, and nothing else to fetch:

```
cd lucifer/examples/steady_state && lucifer lucifer.in
```

The one that needs Genesis4, and the only one with a runner, is the three-way saturation comparison:

```
cd lucifer/examples/saturation_demo && ./run.sh
```

The validation harness is a separate thing, described with its commands in [`validation.md`](validation.md).

## Files

| Path | Contents |
|---|---|
| `lucifer/doc/fel-physics.md` | **The physics manual**: equations, conventions, Genesis4 provenance and validation pointers, one section per subsystem |
| `lucifer/code/fel_beam_mod.f90` | Packed particle slices in Bmad coordinates plus per-particle weight, openPMD `.beam.h5` dump read/write through Bmad's own beam I/O, copy-only `coord_struct` conversion, weighted beam diagnostics with `N_eff` |
| `lucifer/code/fel_track_mod.f90` | The transcribed FEL step: transverse push with natural focusing, RK4 ponderomotive advance, source deposition, FFT field solve; the rotating-record slippage machinery (`fel_slip_struct`, `fel_apply_slippage`, `fel_field_index`); plus the transcribed Genesis4 interlude model |
| `lucifer/program/lucifer.f90` | The tracker: walks a Bmad lattice, FEL steps inside wiggler/undulator elements tracked by an FEL method (parameters from the lattice attributes, described in the FEL element section), seam everywhere else, slippage schedule transcribed from Genesis4; loads the beam from a `beam_init` bunch and generates the seed field when no dumps are named |
| `lucifer/examples/` | Self-contained single-command examples, one directory per feature, each with its own README and measured numbers. Its index says which example shows which feature |
| `lucifer/tests/scripts/check_shot_noise.py` | Statistical check: `<\|b(h)\|^2> = 1/N_lambda` over many seeds, uniform and nonuniform weights, in keep mode, and on copies and deposit cells grouped by the loader, plus the stats file's mode split against a field dump |
| `lucifer/tests/scripts/check_load.py` | Load-path checks: keep mode exact against the bunch, split-weight invariance, the in-cone startup power of a Gaussian bunch, the refusals of noise on unquiet loads |
| `lucifer/tests/scripts/bunchfile.py` | Reads and writes the openPMD-beamphysics bunch files `beam_init%position_file` reads, for the checks that build a bunch outside the loader |
| `lucifer/tests/scripts/check_sase_startup.py` | Cross-code check: SASE startup power, our loader against Genesis4's, independent RNGs |
| `lucifer/tests/scripts/check_migration.py` | Migration checks: charge conservation under heavy migration, exact phase continuity, window residency, no-op bit identity, and migration under a full carrier turn of `z_offset` |
| `lucifer/code/fel_import_mod.f90` | The resampler: a bunch_struct resampled into FEL slices by Genesis4's importdistribution method, the sample mode from a file, plus the Genesis4-distribution-file writer. Also see the resampler in [`validation.md`](validation.md) |
| `lucifer/tests/scripts/check_seam_wake.py` | Seam-wake checks: closed-form pseudomode, lord resolution on a split pipe, exact causality cross-checked against the transcribed wake model, the z_long kernel against a numpy convolution, split-weight invariance, thread determinism |
| `lucifer/tests/scripts/check_import.py` | Resampler checks, exact: the current profile against Genesis4 on the same file, the generated moments, split-weight invariance, the openPMD round trip, the zero-charge refusal, thread determinism and the quiet-loaded current. Statistical: Twiss recovery, startup power and the imported current against the analytic Gaussian |
| `lucifer/code/fel_collective_mod.f90` | Wakes and space charge at Genesis4's granularity: the numerical resistive-wall impedance (Bane-Stupakov, a separable future Bmad port), geometric and roughness kernels, the causal convolution, the per-slice eloss application, and the short/long-range space-charge solvers behind a swappable interface |
| `lucifer/tests/genesis4/collective/` | Genesis4 decks: the collective tiers, importing the shared time-dependent dumps |
| `lucifer/tests/scripts/check_collective.py` | Collective checks: exact wake energy bookkeeping with sigma_energy invariance, stale-wake structure under migration, a stride that only drops, the charge representation under splitting, the long-range kernel at zero width, the short-range solve on a collapsed slice |
| `lucifer/code/fel_unaveraged_mod.f90` | The unaveraged mode: full Newton-Lorentz quiver through the analytic undulator field with sin² end ramps, radiation kick + coupling-free source, energy ledger; fc and faw appear nowhere in it |
| `lucifer/tests/scripts/check_unaveraged.py` | Unaveraged checks: energy ledger, ballistic/handoff, fc measured vs closed form (planar, helical, h=3), step-size convergence, priced gain-curve comparison, the mixed line, thread determinism, the time-dependent ledger and the single-precision twin |
| `lucifer/tests/bmad/unavg_probe_*.bmad` | The paired coupling probes (12/20 periods, planar and helical) |
| `lucifer/tests/run_fel_benchmark.sh` | The whole validation, one command |
| `lucifer/tests/run_perf_benchmark.sh` | Performance head-to-head vs Genesis4, serial and at the machine's performance-core count. Also see the parallelism record in [`validation.md`](validation.md) |
| `lucifer/tests/scripts/compare_fel.py` | Comparison: four steady-state tiers and six time-dependent tiers against Genesis4, plus the split-weight invariance check |
| `lucifer/tests/scripts/plot_fel_compare.py` | Visual companion to `compare_fel.py`: overlays one tier's Bmad and Genesis4 curves (power, checked relative difference, per-slice exit power, bunching) from the diag file and the `.out.h5` |
| `lucifer/tests/genesis4/time_dependent/Aramis-td-sase.in` | Genesis4 deck: pure SASE, the time-dependent window with the seed removed (`power = 0`), writing its own shared dumps |
| `lucifer/tests/genesis4/steady_state/Aramis-ss.in`, `genesis4/Aramis.lat` | Genesis4 deck: Benchmark1-SASE steady state, modified as documented in the deck header |
| `lucifer/tests/genesis4/steady_state/Aramis-1seg.in`, `genesis4/Aramis-1seg.lat` | Genesis4 deck: one undulator segment, importing the same dumps |
| `lucifer/tests/genesis4/time_dependent/Aramis-td.in`, `Aramis-td-1seg.in` | Genesis4 decks: the time-dependent pair, 32 slices with shot noise |
| `lucifer/tests/bmad/aramis.bmad`, `aramis_1seg.bmad` | The Bmad lattices: real wiggler elements, `b_max` encoding aw = 0.84853 exactly in Bmad's constants |
| `lucifer/tests/genesis4/sweep/Aramis-td-s12.in`, `run_delz_sweep.sh` | The coarse-step measurement: one shared dump at `sample = 12`, tracker runs at several `ds_step` values (wrapper lattices overriding that attribute) |

## Architecture

The physics is in the manual: [](fel-physics.md#sec-core), [](fel-physics.md#sec-chart), [](fel-physics.md#sec-field) and [](fel-physics.md#sec-slippage).

Inside FEL elements (wigglers tracked by an FEL method, described in the FEL element section), `fel_track_und_step` advances the coupled system in steps of the element's `ds_step`. That is Bmad's standard step attribute, which is Genesis4's `delz` living on the lattice like every other parameter, and the bookkeeper's `num_steps = round(l/ds_step)` is exactly Genesis4's unroll. The order is Genesis4's: transverse half step, RK4 advance of (theta, gamma) with the field gathered once per step, transverse half step, then source deposition and the `exp(K2 dz)` field solve. Bmad tracking is never used inside, since `symp_lie_bmad` resolves the wiggle motion the period-averaged map assumes away.

Everywhere else (the seam) the packed slice converts to `coord_struct`s by plain copies, `track1_bunch` tracks them, and the field goes through `wavefront_drift`.

**Coordinates.** The packed arrays store Bmad's `(x, px/p0, y, py/p0, z, pz)` plus a per-particle weight (macroparticle charge in Coulombs, mapping to `coord_struct%charge`). The ponderomotive phase is derived and not stored. Genesis4's per-particle theta splits into a common reference advance (one scalar per beam, `phi0`, advanced once per step) and the particle-specific lag carried by Bmad's z:

```
theta_j = phi0 - ks*tau_j,    tau_j = -z_j/beta_j = c*(t_j - t_ref)
```

This reference-offset formulation is the FP32-safe one, and it removes the wrap hazard outright: z does not wrap, so slice migration has no theta-wrap-plus-index update to get wrong. The longitudinal RK4 still runs in Genesis4's (theta, gamma) chart as per-step working variables. theta <-> z is an affine map, under which RK4 is exactly invariant, so the verbatim transcription survives. The per-step gamma <-> pz round trip costs about 1 ulp of gamma per step. With Genesis4's constants transcribed, that cost puts tier1 at 2.8e-11 where the Genesis4-coordinate representation sat at 2.5e-12. The particles live in the packed arrays at all times except inside `track1_bunch`. `coord_struct` never enters the FEL step loop.

**Weights.** Every reduction is weighted: the source deposition scales per particle as `c*w_j/slice_spacing` (Genesis4's `current/N` for uniform weights), bunching is `|sum w e^(i theta)|/sum w`, and `N_eff = (sum w)^2/sum w^2` is a per-slice diagnostic. The Genesis4 dump format carries no weights, so imports are uniform and a Genesis4 comparison can only exercise the uniform case. The split-weight check below covers the rest.

**Units and constants.** The field is in V/m throughout (the `wavefront_struct` convention) and all constants are Bmad's (`m_electron`, `mu_0_vac*c_light`). In those units power is `sum|E|^2*dA/(2*Z0)` and the coupling is `fc*conj(E)/(sqrt(2)*m_electron)`. Genesis4's impedance constant is truncated, 376.73 against the exact 376.7303..., so the two codes differ by 8.3e-7 relative before any physics runs. That difference is the floor of every Genesis4 comparison, and the checks are sized to it rather than to round-off.

**Time dependence and slippage.** A multi-slice starting dump makes a time-dependent run, a single-slice dump the steady state. There is no separate switch, the same rule as Genesis4, whose imports carry the time window. Beam slice `is` couples to field slice `1 + mod(is-1+first, nslice)`. The field record is a circular buffer over the wavefront's slice index, rotated by slippage rather than moved, as Genesis4 does. Slippage accumulates at `dz*(1+aw^2)/(2*gamma0^2*lambda)` wavelengths per undulator step and rotates the record one slice whenever the accumulation exceeds `0.8*sample`. That is Genesis4's rule reduced to one shared-memory node, where the MPI ring exchange is the identity. The slice rotating out at the head of the time window is discarded and re-enters zeroed at the tail. Drifts autophase `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths onto the last interlude element before each undulator. The end-of-lattice fixup is unguarded in Genesis4 (`+1` even with no trailing drift at all). A transcription with a guard Genesis4 does not have loses one final rotation, a difference of 0.84 in the final field. Everything reading the record in time order (the per-slice diagnostics, the final field dump) unrotates through `fel_field_index`, as Genesis4 does on the way out. Beam slices never rotate.

## The FEL element: parameters live on the lattice

The physics is in the manual, [](fel-physics.md#sec-element).

An FEL segment is a real Bmad `wiggler` (or `undulator`) element carrying `fel_method = averaged` or `unaveraged`, Bmad's own attribute for the two FEL methods, which this driver supplies through the same hook Bmad uses for `custom`. Recognition is by key and FEL method, never by name. A wiggler carrying `tracking_method = custom` is some other program's element, and this driver leaves it to the seam. There are no per-undulator namelist parameters. The FEL parameters derive from the element attributes:

- `aw` (rms, Genesis4's convention) from the peak field: `K = c*b_max/(k_u*m_e c^2)`, exactly and independent of the reference energy. Helical undulators have `aw = K`, and planar undulators have `aw = K/sqrt(2)`. The benchmark lattices write `b_max` as an expression in Bmad's own constants so it encodes `aw = 0.84853` exactly, and the round trip still differs by 1 ulp, measured below.
- Helicity from `field_calc` (`helical_model` / `planar_model`), never from the stored `k1x`/`k1y` wiggler attributes, whose helical sign convention disagrees with Bmad's own tracking locals.
- Genesis4's natural-focusing split from the helicity defaults (`kx = ky = 0.5*k_u^2` helical, `0`/`k_u^2` planar). Bmad's `kx` roll-off attribute is not yet mapped onto that split and must be zero.

Outside the driver's own FEL walk the element is just a periodic wiggler, tracked by Bmad's standard kernel, because its `tracking_method` says `bmad_standard` and the FEL method rides alongside. The reference time acquires the resonant undulation delay from Bmad's own code and the transfer-matrix bookkeeping works, with nothing wired in. The driver still supplies `track1_custom_ptr` and `make_mat6_custom_ptr`, both delegating to `track_a_wiggler`, for a wiggler whose `tracking_method` is `custom` and so names some other program's tracking. Both are needed when that happens. `mat6_calc_method` resolves to custom alongside the tracking method, and Bmad calls through a null `make_mat6_custom_ptr`, a jump to address zero, if a program sets only the tracking hook.

Three assertions guard silent failures where the driver claims the element, in the setup pass that reads `fel_method`. They refuse a missing `b_max`, where Bmad's own kernel would silently give `osc_amplitude = 0` (no field, no resonance, no error). They refuse a missing `l_period`, the same silence. And they refuse a fieldmap `field_calc`, which segfaults the parse-time reference tracking if allowed through. Any later is too late, since the missing-`b_max` lattice otherwise dies downstream on an unrelated generation message. The assertions live in one routine (`fel_assert_wiggler_sane`) and nowhere else, since a second copy would mask the removal of either under mutation testing. The same routine also refuses a super_slave undulator, a nonzero roll-off attribute, a tilt on a helical element and an integration step equal to Bmad's fallback for an unset one, which a wiggler with a period never carries. All three refusals are permanent harness checks.

The element-driven full time-dependent line is anchored against the namelist-driven reference: it reproduces the earlier namelist-driven run to a max relative difference of 3.4e-12, zero at z = 0 and growing with gain. That is the amplification of the 1-ulp `aw` difference between deriving `aw` through the `b_max` attribute round trip and reading it from input.

The FEL mode is the element's `fel_method`, class-settable as `wiggler::*[FEL_METHOD] = unaveraged`, so averaged and unaveraged segments mix freely in one line. `averaged` uses the transverse maps of Bmad's own `bmad_standard` periodic-wiggler kernel, flattened per `ds_step` (`track_a_wiggler`'s matrix with the octupole-like end kicks, chromatic via `p0/p`). `unaveraged` is a production method whose cost buys the full quiver dynamics (see the manual's unaveraged-mode section).

Behind the averaged method there is one more choice, and it is a run switch rather than a method. `global%transport_model = "genesis"` replaces those maps with the transcribed-Genesis4 focusing (matrix from `aw`, `kx`, `ky`, chromatic via `gammaz`). It is validation-internal, since the Genesis4 tiers need transcription-level transport, and no production run sets it. The two averaged models are priced, measured over the full time-dependent line (32 slices, 90 records). Power differs by 5.0e-5 max (exit total power 3.0e-7), on-axis intensity 7.3e-4, spot sizes 6.7e-6 and wrapped bunching phase 1.3e-2 rad max, for +3.8% runtime. That is period-averaged Genesis4 focusing against Bmad's end-field treatment. The unaveraged mode's two numbers are run switches for the same reason: `global%unaveraged_steps_per_period` (20, below 10 refused) and `global%unaveraged_ramp_periods` (2). A true hard edge, the test configuration, is the explicit sentinel -1, because a silent hard edge would reintroduce the K/gamma handoff hazard.

`examples/taper/` uses this per-element mixing. It is the same line with the last two cells' undulators a second element definition with `b_max` 0.4% lower, bit-identical to the untapered run until the taper starts and 12.7x its exit power after (see `examples/README.md`).

## Program structure

The tracker is laid out the way Tao is laid out. The input structs are in `fel_struct`, with defaults in the declarations. The namelist layer is quarantined in `fel_input_mod`, three groups in one file: `&fel_params` with `global%...`, `slicing%`, `bmad_com`, `space_charge_com`, `chamber_wake%` and `space_charge%`, `&fel_beam_init` with `beam_init%...`, `resample%...` and files, and `&fel_wavefront_init` with `wavefront_init%...` and `field_file`. The work is in library modules (`fel_setup_mod` / `fel_init_mod` / `fel_track_line_mod` / `fel_io_mod`) over one explicit `fel_run_struct`. Nothing in the library stops. Errors return and the driver decides, and `track_fel_line` is re-entrant: twice in one process is bit-identical to two processes, and `tests/lucifer_smoke_test.f90` drives the library with no namelist anywhere. A flat `&fel_track_params` group is refused, and the message names the three groups that replace it. `global%comb_ds_save` is the stats comb, with Bmad's `ds_save` semantics and one deliberate difference. An element end is always a record, whatever the comb, which is what lets the file carry one record axis and a mask instead of a second axis. The default, 0, is a row at every record. `global%track_start`/`track_end` bound the walk over a schedule always built on the full lattice, so windowed runs compose exactly. `global%continuation = T` declares that `beam_file` and `field_file` are the two frames of one checkpoint and the run continues from it, `track_start` opening at the element after the one the checkpoint completed. Off, the same files initialize a fresh run from their standard records ([reading-output](reading-output.md)). `stats.h5` carries `meta/` provenance (the resolved input echo, which lattice, a timestamp and the Bmad version) as datasets. `global%record_environment` adds the user name and working directory for a lab notebook, off by default because a stats file travels.

## Output: stdout is for humans, files are for programs

The program's terminal output is formatted for reading and carries no contract beyond that: units chosen for the eye, columns aligned, wording free to improve. Parsing it is discouraged. Everything a program should read is written to a file, at full precision:

| File | What |
|---|---|
| `<out_root>.stats.h5` | The run: per-record beam, field and Twiss data, the lattice table, and the axes and parameters to read them by. Self-describing (see [`reading-output.md`](reading-output.md)) |
| `<out_root>-final.beam.h5`, `-final.wf.h5` | Final beam and field, openPMD, the only dump format this code writes. The beam file carries the per-particle weight, which no Genesis4 dump can. At an eligible boundary both carry the `lucifer` checkpoint group a continuation reads (`global%continuation = T`). `lucifer/tests/scripts/convert_genesis.py` converts either kind to Genesis4 conventions, for feeding Genesis4 |
| `<out_root>.diag.txt` | The per-record Genesis4-comparison instrument, one row per slice per record |
| `<out_root>.ledger.txt` | The unaveraged energy ledger, one row per record step |
| `<out_root>.import.txt` | The load record: the per-slice current profile, with the bunch's analysis moments from the resampler or the per-slice counts and first moments from the keep mode. Written at load time, so it exists under `load_only = T` |
| `<out_root>.migration.txt` | Slice migration: one row per event (s, particles moved, charge dropped, phase-continuity residual, the smallest phase distance any examined particle had to a slice boundary) plus the run summary |

This is Tao's division of labour. `show ele` is for people, and `show value` and pytao are for precision. It is why the validation suite reads files rather than scraping the screen.

A run prints a framed block naming the configuration, a progress table with a wall time per row, and a completion block listing what was written. The lines below are the steady-state example's, written from its log by the examples runner with the progress rows left out, since their wall times belong to the machine that ran it:

```{include} generated/sample-output.md
```

Four choices in the progress table, which the sample leaves out, are deliberate. Values carry an SI prefix chosen per value, so one column holds a kilowatt startup and a megawatt saturation without losing the early digits a fixed unit would flatten. The element name is last, unpadded, so every numeric column is fixed however long a name gets: no truncation rule, no name shoving the numbers rightward. The step column is blank for one-step elements, which is most of a real lattice. And the distance column is `s`, Bmad's arc length along the reference orbit. Through a bending break the arc exceeds the chord the light takes, which is what the light-path correction accounts for.

The Radiation line reports the grid the run built, which is the file's when the field came from one and the derived one when the deck stated none ([](input-reference.md#param-wavefront-grid)).

The Work, Estimate and Load lines are the pre-run. The sample above leaves the Estimate line out, since it names a machine and a date. Work counts what the run will do, and every number in it is known before the first step. The numbers are the integration steps over the line, the macroparticles in the window, their product, and the transverse grid points the field solve touches per step over the window and every member of the field set. Estimate turns that into wall clock at rates measured on one machine, which the line names, from the particle push and the field solve ([](performance.md)). It is an estimate and it is checked: the Cost line in the completion block prints the clock beside it on every run. Load says where the macroparticle count stands against the load the source filter's convergence was measured at, 1024 macroparticles in 128 beamlets. It also says what the beamlet size buys, which is the harmonics the shot noise is imposed at, `(beamlet_size - 1)/2`. A harmonic field above that is refused when shot noise is on, with the beamlet size that would carry it.

The Frames line appears when `global%dump_at_comb` is on, and states how many frames the run will write and roughly how many bytes they take. A comb of zero is a frame at every integration step, so read the number before the run rather than after the disk fills ([](input-reference.md#param-global-dump-at-comb), [](reading-output.md)).

Setting `global%load_only` stops the run after this block. That is the pre-run on its own: everything the deck derives and everything it will cost, with the initial state written and nothing tracked ([](input-reference.md#param-global-load-only)).

The Split block is the run's own convergence report. The field power is separated into the part inside an angle of about four mode diffraction angles and the part outside it. The inside part is what couples to the mode. The outside part is the emission of macroparticles that occupy one grid point each. That second part is an artifact of the representation. It grows as the cells shrink and falls as the beamlet count rises, and at common settings it is most of the total power of an unseeded run. The angular filter that removes it is on by default, and what it costs and what load it wants is [](fel-physics.md#sec-convergence). The block reports both at the record where the mode peaked and at the last record. The run then warns when the outside part reaches a tenth of the inside part and confirms convergence below a hundredth. On the steady-state example with the filter on, 0.2 percent of the 792.7 MW exit power is that emission and 0.08 percent of the 1.56 GW at saturation. The same numbers are in the statistics file under `run/`, and the per-record per-slice arrays are `field/total/power_inside_angle` beside `field/total/power`.

Below the two records the block names the thin slices, by index and position, up to eight of them with a count of the rest. A thin slice carries charge and has a beamlet count under 128, or a power outside the split angle over a tenth of the power inside it. A window ratio is a mean over the window, so a thin tail disappears into it, and a beam imported with irregular counts per slice is where that happens. The remedy the line names is `load_mode = "sample"`, which draws the same number of beamlets in every slice ([](input-reference.md#param-beam-load-mode)).

The physics columns are filled whatever the comb setting. With per-record rows they come from the current record, and with `comb_ds_save < 0` from the element-end row, which is filled just before the row that reads it. Between element ends in that mode the last element-end values are carried forward rather than blanked. A mid-element row exists to answer "is it alive and roughly where".

There is one deliberate exception. A check has to recognize a refusal, so the checks match the ALL-CAPS text of a refusal together with a nonzero exit status. Those texts are a supported contract and change only together with their checks, in the same commit. Nothing else printed to the terminal is.
