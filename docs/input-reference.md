# Input reference

This page is the normative reference for every key of the deck: its default, its meaning, and what refuses it. A run is one TOML file. Its top level names the lattice, and its tables set the program's structures and two of Bmad's directly, which is Tao's `&tao_params` pattern in TOML. The key tables on this page are generated from the schema the program writes with `lucifer --schema`, so the names, the types, the defaults and the bounds are the program's own, and the prose between them says what the tables do not: what each key is for and what refuses it. The schema is `docs/generated/schema.json`, and an editor that reads JSON Schema validates a deck against it as it is typed.

| Table | Carries |
|---|---|
| `lattice` | The one top-level key, the Bmad lattice file |
| [`[run]`](#group-fel-params) | The run switches: the output root, the models, the window of the walk, the dumps, the seed and the check instruments |
| [`[window]`](#param-slicing) | The longitudinal discretization: the time window and the slice spacing |
| [`[beam]`](#group-fel-beam-init) | The bunch description, the load and its switches, the files, and `[beam.resample]`, the resampler |
| [`[field]`](#group-fel-wavefront-init) | The radiation starting condition and the field files |
| [`[wake]`](#chamber-wakes) | The chamber wake, with Genesis4's `&wake` names |
| [`[space_charge]`](#space-charge) | The space-charge solver |
| [`[bmad]`](#bmads-own-switches) | The Bmad switches the tracker honors, `bmad_com` and `space_charge_com` |

Per-element settings are [lattice attributes](#lattice-attributes) rather than deck keys. How to build and run is [](user-guide.md). What the outputs hold is [](reading-output.md). The measured levels are [](validation.md). The physics is the manual, [](fel-physics.md). Every resolved input is echoed into the stats file's `meta/input_echo` as a complete deck, defaults made explicit, and recorded key by key under `params/`.

A key the schema does not name is refused with the keys of its table, so a misspelling cannot pass as a default. A key in the wrong table is refused with the table it belongs in. A value of the wrong type or outside its bounds is refused with what was expected. A deck is read whole before anything runs.

(group-fel-params)=
## The lattice and `[run]`

```toml
lattice = "aramis.bmad"

[run]
out_root = "fel_td"
```

```{include} generated/input-keys.md
:start-after: <!-- begin keys-top -->
:end-before: <!-- end keys-top -->
```

```{include} generated/input-keys.md
:start-after: <!-- begin keys-run -->
:end-before: <!-- end keys-run -->
```

### The lattice and the run


(param-global-interlude-model)=
**`run.interlude_model`** selects how the field-free elements are handled. `"bmad"` is the seam: Bmad's own `track1_bunch`, an exact theta mapping, and `wavefront_drift` for the field. `"genesis"` uses the transcribed interlude step everywhere, which measures what the seam changes. The slippage schedule is identical in both models. The physics is in the manual's interlude and seam sections.

(param-global-transport-model)=
**`run.transport_model`** selects the transverse maps inside an averaged FEL element. `"bmad"`, the default and the production model, is Bmad's own periodic-wiggler kernel with its end-field treatment and chromaticity through `p0/p`. `"genesis"` is Genesis4's maps verbatim, with the focusing matrix built from `aw`, `kx`, `ky` and chromaticity through `gammaz`. The second is validation-internal. A transcription comparison needs transcription-level transport, so the comparison tiers set it and no production run does. The two are compared in the manual's element section. The unaveraged mode integrates the field rather than applying a map, so the switch does not reach it.

(param-global-unaveraged-steps-per-period)=
**`run.unaveraged_steps_per_period`** is 20 by default. Below 10 is refused. That floor is where the coupling-factor convergence check stands: [](validation.md) tabulates the factor at 10, 20 and 30 steps per period. The setting reaches unaveraged elements only. Averaged elements take their step from `ds_step`.

(param-global-unaveraged-ramp-periods)=
**`run.unaveraged_ramp_periods`** is 2 by default, and 0 takes the default too. A true hard edge is a test configuration and has the explicit sentinel `-1`, so silence never means hard edge. The pre-run block's Unaveraged line states the ramp the run took and where it came from. A ramp pair longer than the segment is refused. It is a run switch, the same for every element. A per-element attribute this program registered made every other Bmad program refuse to parse the lattice. Averaged and unaveraged elements still mix per element, through `fel_method`.

(param-global-source-model)=
**`run.source_model`** is `"deposit"`, the standard per-particle scatter, or `"coherent"`, the coherent-Gaussian retrieval. In the coherent retrieval the spatially incoherent part of the source is dropped, the slice bunch factor keeps the physical shot noise, and the transverse shape is a guarded Gaussian. A profile that is measurably not Gaussian is refused. See [](fel-physics.md#sec-coherent-source).

(param-global-source-filter)=
**`run.source_filter`** filters the source term at wide transverse angles. Each beamlet is a single transverse point on the deposition grid, and the source has no angular dependence. So a beamlet radiates into every transverse wavenumber the grid carries, well outside the undulator's central cone. On, the transformed source is multiplied by a sigmoid in normalized spatial frequency before it is added to the field. The field itself is never multiplied: a seed, an imported field and everything already radiated propagate exactly as they do with the filter off. `source_filter_angle` places the sigmoid's half-height, in radians. At zero, the run derives it as the larger of two angles. The first is four diffraction angles of the fundamental mode, `4*lambda0/(2*pi*sigma)`, with `sigma` the rms transverse size of the loaded beam. The second is the angle at which the resonant wavelength red-shifts by the Pierce parameter, `sqrt(2*rho*lambda0/l_period)`. Inside the first lies the radiation that can couple to the mode. Outside the second, emission is beyond the bandwidth of the interaction. Their ratio goes as `6.8/sqrt(z_R/L_g)`, measured to 5 percent over four machines. So the mode angle is the larger below a `z_R/L_g` of about 46, and neither alone serves for every machine. The larger is taken because cutting into the mode loses real radiation while leaving artifact in only weakens the filter. That larger angle is the one protected. `source_filter_tolerance` places the edge relative to it: the tolerance is the largest fraction of source intensity the filter may take anywhere inside the protected angle. A sigmoid whose half-amplitude point sits on that angle passes a quarter of the intensity there. So the default of 0.75 puts the half-amplitude point on the protected angle, to the bit. A smaller tolerance moves the edge out by `1/(1 - width*log(t/(1-t)))` with `t = sqrt(1 - tolerance)`, a factor of 1.12 at 0.2 and 1.36 at 0.01. That holds the transmission across the protected angle at the cost of leaving more wide-angle source in. The tolerance must lie strictly between 0 and 1. A width too large to reach the transmission asked for is refused, with the width that would reach it. Setting `source_filter_angle` or the grid-relative cuts places the edge directly and bypasses the tolerance. Setup says which happened. The statistics split stays on the protected angle whatever the tolerance, so a scan of it compares against one acceptance. Both angles, their ratio and which one won are printed at setup, with a warning when the ratio leaves 0.3 to 3. `source_filter_width` is the sigmoid's width as a fraction of that edge, so the sigmoid on axis is `1/(1 + exp(-1/width))`. The default of 0.05 leaves the axis untouched. A width that puts the axis below 0.99 is refused, with the width that would pass. `source_filter_xcut` and `source_filter_ycut` place the edge in Genesis4's own units instead, the shifted grid index over the grid size. So 1 is the angle `lambda0/dx`, twice the Nyquist angle and off the grid, and the value moves with the grid. They exist so the transcription check can place the edge where Genesis4 places it ([](validation.md#val-source-filter)). They lift the width refusal with them. Setting them together with the angle is refused. The filter is on by default, since what it removes is emission no physical beam radiates and it is the cheaper route to a converged answer ([](fel-physics.md#sec-convergence)). Turn it off to compare against a code that carries no such filter or against a level recorded without one, which is what the comparison tiers do. A run with it off is bit-for-bit a run with no filter in the code. The filter reaches the averaged elements only. The unaveraged mode deposits from the resolved motion and builds no transformed source to multiply. So a mixed line filters its averaged elements and leaves the rest, and a line at setup says so. `source_model = "coherent"` turns the filter off for the run and says so at setup, since the coherent source is an analytic Gaussian that carries no wide-angle content to remove. Neither is refused, because the filter is the default and a deck that never asked for it must still run. A non-positive cut or width is refused rather than quietly turning the filter off. What the filter is for, what it costs and what load it wants is [](fel-physics.md#sec-convergence). The measurements behind it are [](startup-noise.md). `source_filter_mutate` applies the sigmoid to the propagated field instead of to the source, which is how the benchmark proves the check can fail.

(param-global-load-only)=
**`run.load_only`** is the pre-run. The run parses the lattice, derives everything it derives, builds the beam and the field, prints the header and stops. It writes the initial state as `<out_root>-initial.beam.h5` and `<out_root>-initial.wf.h5` through the same writers the final dumps use, so the tracker reads them straight back. The header it prints is the one a full run prints. It states the grid, the slice spacing and the step with their origins, the work the run will do, and the wall clock that work is estimated to cost on the machine the estimate names. It also states where the macroparticle count stands against the load the source filter's convergence was measured at ([](user-guide.md), [](fel-physics.md#sec-convergence)). A deck's cost and load can therefore be read before any of the tracking is run. `run.write_initial` writes the same two files and then tracks.

(param-global-device-timing)=
**`run.device_timing`** times each pass of a device step separately. The backend's own busy seconds come from one command buffer. A command buffer holds every step between two host touches, so those seconds say what a step cost and not what a pass cost. On, each pass is encoded into its own compute encoder carrying a timestamp counter attachment, still inside the one command buffer the step batching rests on. The run reports seconds and encoder counts for the transverse map, the push, the source clear, the deposit, the filter and the solve. The passes of one kind sum into one line, so the transverse map reports two encoders a step and the filter and the solve four each. Stage-boundary counter sampling is what an Apple GPU carries. It times an encoder rather than a dispatch, which is why a pass gets an encoder of its own. A device that samples nothing at an encoder boundary refuses rather than reporting zeros. The instrument therefore pays encoder boundaries the production path does not, measured at 1.2 to 2.6 times the device busy seconds. So it is off by default, and a timing is taken with it off ([](performance.md#perf-device)). The sample buffer holds 2048 passes, about 157 steps. A command buffer needing more is committed early, with a warning saying how often, since that is a batching change the instrument caused.

**`run.device_dep_mutate`** shifts the deposit's fixed-point scale by that many bits and changes nothing else. That is how `check_device` measures what the quantum contributes to the deposit's error rather than bounding it away ([](validation.md#val-device)). It is zero in every run but that check. A shift down coarsens the quantum and is refused once it would leave the fixed point no finer than the float accumulation it replaced.

(param-global-continuation)=
**`run.continuation`** declares what `beam.file` and `field.file` are for, and nothing is inferred from the files. Off, the files initialize a fresh run from sliced inputs. Only the standard openPMD records are read, the reference phase and the slippage residual start at zero, and a checkpoint group present in either file is ignored. A frame written inside an undulator is refused as a start. On, the files continue an interrupted run. Both must be the two frames of one checkpoint, written at an eligible boundary and carrying the group `lucifer`. `track_start` must open the walk at the element right after the one the checkpoint completed. The reader refuses anything less before tracking, naming the condition. The group's contents and the refusals are in [](reading-output.md).

(param-global-track-start)=
**`run.track_start`** and **`run.track_end`** bound the walk, using Bmad's element-locator syntax. A locator that matches no element, or more than one, is refused. The schedule (slippage, autophasing, break geometry) is always built on the full lattice. So a windowed run composes exactly with the full one: the first span followed by the second, started from the first's dumps, reproduces the one-shot run. The measured level is in [](validation.md#val-the-programs-own-identities).

(param-global-fp32-check)=
**`run.fp32_check`** turns on the FP32 lockstep instrument, a check rather than a production mode. An FP32 twin of the FEL advance steps beside the FP64 path from a shared state, and the per-quantity divergence goes to `<out_root>.fp32.txt` with the worst case in the footer. `"lockstep"` rebuilds the FP32 state from FP64 every step, so a wrong formula shows as a jump. `"freerun"` carries the FP32 longitudinal state and the FP32 field across steps, a complete single-precision run beside the FP64 one. So the compounding rate and the end-to-end observable divergence are measured. It covers a single-slice window and refuses more, since the twin keeps no slippage rotation of its own. The twin covers the fundamental-only, collective-free advance of either mode, each segment by the twin that mirrors it. It refuses anything else (harmonics, two polarizations, the coherent source, wakes, space charge). The unaveraged twin is a lockstep only, since `"freerun"` compounds a single-precision state across steps and that twin carries none. So `"freerun"` on a lattice with an unaveraged segment is refused. With the device in the twin's role the harmonic and two-polarization refusals lift, since the device twin covers the set. A run whose FP32 residual cannot resolve the per-step phase increment is refused by the runtime check rather than reported as a clean number. The FP64 physics is untouched: the instrumented run's outputs are byte-identical to the uninstrumented run's. The measured levels are in [](validation.md#val-fp32-lockstep). `fp32_mutate` coarsens the FP32 residual by eight mantissa bits so the benchmark can prove the check fails when the path is wrong. The instrument needs the single-precision FFTW (`libfftw3f`), since the field twin's transform is FFTW's own `fftwf` interface. A build without that library refuses `fp32_check` rather than transforming in another precision, and the build reports which way it went at configure time.

(param-global-device)=
**`run.device`** selects a device backend for the averaged FEL step. `"metal"` runs it on an Apple Silicon GPU in single precision. The beam and the field upload at each FEL element's entry and stay resident through the element. Every integration step is one command buffer (transverse maps, the longitudinal push in a fixed-point phase chart, the source deposit and the FFT field solve). The state returns to the host at the stats comb's positions and the element end. Between those positions the host arrays are stale by design, so `comb_ds_save` is also the readback knob. Combined with `fp32_check` the device instead takes the lockstep twin's role. The FP64 run is untouched and the device's per-step divergence goes to the same `.fp32.txt` stream, which is how the backend is judged ([](validation.md#val-device)). The field set rides the device whole: every harmonic member and both polarization planes, though not the two together, which is refused for the CPU as well. Slice migration runs with the device. The pass reads the beam back at the element's last step and re-slices it on the host, and the next element uploads the result. A slice that outgrows the device's particle capacity grows it, which the run reports. The unaveraged mode rides the device too, with three kernels of its own for the quiver-resolving push, the kick with its deposit and the ledger's spontaneous term, and the same solve underneath ([](validation.md#val-device-unaveraged)). Everything else the kernels do not cover is refused at setup or first use. That is the coherent source, wakes, space charge, spontaneous radiation, the transcribed transport model and the filter's self-test. It is also two live polarizations inside the unaveraged mode, the escaped-field bank, and a transverse grid that is not a power of two from 64 to 1024. A build without the backend refuses the knob itself. An unsupported configuration stops the run and never quietly takes the CPU path. The backend is built where the toolchain can carry it, which is macOS with a Clang-family Objective-C++ compiler rather than macOS alone, since it is ARC-managed Objective-C++ against the Metal framework. Every other build takes the refusing stub and says so at configure time. Running `lucifer` with no arguments names what the build in hand carries, either the device it found or the reason there is none. The question needs no run and no deck to answer.

(param-global-migrate)=
**`run.migrate`** is off by default. The reason is the comparison rather than the physics. The tiers that compare against Genesis 1.3 Version 4 (Genesis4) run against a code that never migrates, so migration inside a transcription-level comparison would be a model difference. Dropped charge is counted and reported per event. With `migrate_check = true` the run also verifies exact phase continuity at every migration. See [](fel-physics.md#sec-migration).

### Output switches


(param-global-dump-beam-at)=
**`run.dump_beam_at`** and **`run.dump_field_at`** name elements through Bmad's own locator, so `class::name` syntax works. An entry matching no element is refused. Dumps are openPMD, and the field dump is unrotated into time order first.

(param-global-dump-at-comb)=
**`run.dump_at_comb`** writes the beam and the field set at every position the stats comb takes a row, as `<out_root>-<record>.beam.h5` and `<out_root>-<record>.wf.h5` with the record zero padded to six digits. That is openPMD's fileBased series. The index is the stats record's, so a frame and the row describing it are read on one axis. The writers are the ones the element-end dumps use, so a frame taken at an element end is the same file that element's `dump_beam_at` entry would write. `comb_ds_save` therefore sets the frame spacing. A comb of zero is a frame at every integration step, which on a long line is tens of gigabytes. The header's `Frames` line states the count and the bytes before the run tracks. A frame inside an unaveraged segment is written in the chart that mode tracks in, with the undulator quiver in `px` and `py`, and the frame's `felMethod` says so. The averaged physics and the Bmad seam still refuse that chart, since they would read it as the averaged one. The beam dump of `dump_beam_at` keeps that refusal because it writes no `felMethod` to distinguish them. That refusal is never reached: `dump_beam_at` names elements through Bmad's locator and so lands on element boundaries, where the ramps have handed the averaged convention back.

Every frame carries where it was taken: `sPosition`, `elementName`, `elementIndex`, `phi0`, and inside an FEL element `felMethod` with `aw`, `ku`, `tilt` and `helical`. The averaged mode integrates the quiver away, so a reader that wants the physical orbit rebuilds it from those ([](reading-output.md)). A frame at the entry face carries a position and a name and no undulator.

Writing frames does not change the run. The field's records are rotated to time order to be written and rotated back, so a run with frames is dataset-identical to the same run without them.

**`run.dump_slice_first`** and **`run.dump_slice_last`** cut the series to a range of the window, in window order, with `-1` as the last meaning the last slice. Frames then carry the particles and the field of those slices only. The beam file's patch count is the range's length, and both files state `sliceFirst` and `sliceLast`. A run at one wavelength per slice has thousands of slices, and a whole-window field frame is gigabytes. So a view that wants the core of the bunch asks for it here rather than reading it out afterwards. A range outside the window, or with the first past the last, is refused. Nothing else is affected: `dump_beam_at`, `dump_field_at` and the initial and final dumps stay whole, because a run restarts from those. A frame from a sub-window, or from inside an element, is not a restart point.

(param-global-dump-reduced)=
**`run.dump_reduced`** writes the field's intensity projections into the statistics file, one set per member and plane of the field set, at every record. They are the transverse intensity summed over the window, the intensity against slice and x summed over y, the same against slice and y, and the complex field on axis per slice. They are what a picture of a run is drawn from, summed once from the field the tracker already holds rather than recomputed from a raw field frame. They integrate back to the record's `power` exactly. The group is [](BMAD-STATS-EXT-FEL.md)'s F16 and the axes it needs, `coords/grid_x` and `coords/grid_y`, appear with it. It is off by default because the transverse projection is the grid's own size at every record. The header's `Reduced` line states the bytes before the run tracks.

(param-global-keep-escaped-field)=
**`run.keep_escaped_field`** banks each field slice that slippage transmits beyond the window, with its `wavefront_params` and transmission position. At finalize it propagates each to the exit plane to write the whole pulse. Field that has left the window never re-interacts, so it is fixed information. See [](fel-physics.md#sec-stats).

(param-global-comb-ds-save)=
**`run.comb_ds_save`** is Bmad's `bunch_track_struct%ds_save` name and semantics. Negative keeps no per-record rows at all, while element ends, dumps and the final state remain. Zero, the default here, records every position the walk reaches. Positive records a row once z has advanced that far past the last one, with element ends always kept. The default deviates from Tao's, deliberately: the per-record arrays are this program's primary statistics contract where Tao's comb is an optional extra.

Note: this is the selecting half of Tao's comb and not the whole of it. `save_a_bunch_step` selects among the positions a tracker already reaches, and that is what this setting does. Tao's beam-track loop also cuts an element into pieces where the tracker would otherwise cross it in one map, and here that is `interlude_ds_step` below. An FEL element reaches one position per `ds_step`, so this setting alone gives a fine comb inside a wiggler. An interlude is one tracking call over the whole element, so this setting alone gives one row there, at its end, whatever it is set to.

**`run.interlude_ds_step`** cuts an interlude element into pieces of about this length, so statistics rows and openPMD frames land along a drift or a quadrupole the way they already land along a wiggler. Zero, the default, tracks the element whole, which is what the recorded validation digits are taken at. The pieces come from Bmad's own `element_slice_iterator`, the routine Tao's comb uses, so an element's fringes stay at its own ends rather than appearing at every cut. The count is Tao's expression, `max(1, int(1.01 * L / ds))`, so an element shorter than the step is one piece. The run's header states the total before it tracks. Each piece tracks the beam, drifts the field, advances the common phase and takes its share of the element's slippage, so the element's exit is the element's exit however many pieces crossed it ([](validation.md#val-interlude-steps)). A chicane's chord-versus-arc correction is a property of the whole break and lands on its last element's last piece. Two configurations decline the cut. An element carrying a Bmad short-range wake is tracked whole, since that wake is a once-per-passage kick of the element's length and would otherwise become one kick a piece. Tao makes the same judgment where CSR or space charge would see a cut. `interlude_model = "genesis"` refuses the setting outright, that model being a transcription whose step is Genesis's own element.

(param-global-record-environment)=
**`run.record_environment`** is off by default because a statistics file is meant to travel, attached to a paper or mailed to a collaborator, and the user name and working directory identify a person and a machine. The timestamp and the Bmad version identify the run without them. See [](fel-physics.md#sec-meta).

### `[bmad]`: Bmad's own switches

The `[bmad]` table sets the `bmad_com` and `space_charge_com` switches the tracker honors, as Tao sets them, and their defaults are Bmad's. Two matter most here.

```{include} generated/input-keys.md
:start-after: <!-- begin keys-bmad -->
:end-before: <!-- end keys-bmad -->
```

Both are off by default, matching Genesis4's `&sponrad`. Fluctuations with `run.migrate = true` are refused: the quiet start cancels per beamlet, and migration scrambles the grouping. See [](fel-physics.md#sec-eom).

### `[wake]`: the chamber wake

The names are Genesis4's `&wake` names. Bmad element wakes are a separate mechanism, described under [element wakes](#element-wakes) below.

```{include} generated/input-keys.md
:start-after: <!-- begin keys-wake -->
:end-before: <!-- end keys-wake -->
```

(param-write-wake-kernels)=
**`wake.write_kernels`** names a file, and setup writes the transcribed kernels into it after they are built. The export builds matching `z_long` tables, and `examples/bmad_wake` uses one.

(param-chamber-wake-model)=
**`wake.model`** names the implementation. `"genesis"` is the only value accepted: the transcribed solver, which convolves its kernels with the weighted slice currents and produces one energy loss per slice. Anything else is refused rather than treated as the transcribed solver by default. The field exists so that a second implementation can arrive as a value here rather than as a rework.

The three kernels, their numerical impedance and the causal convolution are in [](fel-physics.md#sec-wakes).

### `[space_charge]`

Whether space charge acts is not in `[space_charge]`. It is the element's own `space_charge_method`, described with the other [lattice attributes](#lattice-attributes), and this group holds only the solver's numbers. A run with no element asking for it never reaches them.

```{include} generated/input-keys.md
:start-after: <!-- begin keys-space_charge -->
:end-before: <!-- end keys-space_charge -->
```

(param-space-charge-model)=
**`space_charge.model`** names the implementation. `"genesis"` is the only value accepted: the transcribed solver, which works per slice on a radial grid over azimuthal modes and longitudinal harmonics. Anything else is refused. Note: this family is Lucifer's own, and `space_charge_com` in the same group is Bmad's global structure, a separate thing.

The default is expected to stay `"genesis"` when Bmad's own slice solver becomes the second value. The reason is the physics rather than the order they arrived in. The transcribed solver carries two terms that matter inside an undulator and that Bmad's slice model does not have. The first is the space charge of the microbunching itself, solved per longitudinal harmonic of the ponderomotive phase. The second is the longitudinal Lorentz factor, which at aw = 0.85 rms is worth a factor of 1.7 in the field's own scaling. What Bmad's model adds that this one lacks, a transverse defocusing kick, falls as the inverse cube of gamma and is an injector term rather than an undulator-line one. The transcribed path is also the measured one, at the level [](validation.md) records for the space-charge tier. A second value would take the default on a measurement.

`space_charge.ngrid` needs at least two radial points, and `space_charge.nphi`, `space_charge.nz` and `space_charge.rmax` cannot be negative. Each is refused at setup. `rmax` is a floor under the radial scale. Zero lets each slice's own extent set it, which is enough for any slice with transverse size. A slice whose charge sits at one transverse point has none, and the solve is refused there, naming `rmax` as the scale to state.

If neither term is configured, `nz = 0` with `longrange = false`, an element asking for `slice` is refused while Bmad's master switch is on. The solve would run at its full cost and return an exact zero, and the refusal asks which of the two terms was meant.

See [](fel-physics.md#sec-spacecharge).

(element-wakes)=
### Element wakes

Elements carrying Bmad `sr_wake` definitions, either pseudomodes or a tabular `z_long`, act across the whole time window through slice concatenation. The conventions, the step size, the mid-element wiggler kick and the refusals are in [](fel-physics.md#sec-seamwake). Wakes and space charge are both refused in the unaveraged mode, which does not wire them into its step.

(group-fel-beam-init)=
## `[beam]`

```toml
[beam]
n_particle = 2048
a_norm_emit = 4e-7
b_norm_emit = 4e-7
sig_z = 1.2e-9
sig_pz = 8.8e-5
bunch_charge = 3.0e-14
shot_noise = true
```

```{include} generated/input-keys.md
:start-after: <!-- begin keys-beam -->
:end-before: <!-- end keys-beam -->
```

One path loads the beam: a bunch from `beam.`, generated or read from a file, binned into the slices and loaded by `load_mode`, with `quiet_start` and `shot_noise` acting on every load ([](fel-physics.md#sec-loading)). `beam.file` starts from a dump instead.

### Files and sources


(param-beam-beam-file)=
**`beam.file`** reads openPMD only. A file that is not openPMD is refused, and the message carries the conversion command. `tests/scripts/convert_genesis.py` converts in both directions, and [](genesis4.md) describes the exchange. A dump named in `beam.position_file` is refused, since its slices would collapse onto one another when read as a bunch.

### The bunch description

This is Bmad's standard `beam_init_struct`. With `position_file` blank the bunch is generated from the description. In the sample mode the loader evaluates the description analytically per slice without making a bunch at all. The Twiss is always the lattice's.


(param-beam-init-contract)=
The analytic loader honors the fields above and refuses every other `beam_init` field that is set. A standard structure that silently dropped fields would be worse than a custom one. So the honored set is the contract, and `check_beam_init_contract` enforces it. A bunch that Bmad makes, in keep mode or from a file, honors everything `init_beam_distribution` honors.

(param-beam-sig-z)=
**`beam.sig_z`** with `"RAN_GAUSS"` gives a Gaussian current profile evaluated at the slice centers, with the bunch centered in the window. A zero length is the steady state, one slice holding the whole charge, and it is refused for a time-dependent window. With `"GRID"` the profile is flat over the z extent of `grid(3)`. Bmad's `a_emit` and `b_emit` are refused: normalized emittances only, which is Bmad's preferred form. `sig_e` is deprecated Bmad-wide and does not exist here.

### The load


(param-beam-load-mode)=
**`load_mode`** is what becomes of the particles in a slice. In the sample mode a generated bunch takes the analytic loader and a file takes the resampler ([](fel-physics.md#sec-import)). `beam.n_particle` is the count per slice either way. In the keep mode each particle becomes a beamlet of copies at weight over `beamlet_size` sharing its coordinates, so every moment of the bunch is the load's, slice by slice. The window is then the bunch's extent unless `window.` states one. The method is in [](fel-physics.md#sec-loading).

(param-beam-shot-noise)=
**`shot_noise`** imposes the weighted Fawley loading on the groups the load has: the beamlets when the loader made them, else particles sharing their transverse coordinates, else the occupants of a deposit cell, and the message says which. The loader measures the quiet floor first and refuses noise on a load above it with the value in the message, so noise is never counted twice. The algebra, the rule and the four switch combinations are in [](fel-physics.md#sec-noise). The loader warns where Genesis4 silently clamps groups holding fewer real electrons than macroparticles.

### `[beam.resample]`: the resampler

The names are Genesis4's `&importdistribution` names where an equivalent exists. `window.n_wavelength`, `ran_seed` and the seed field are shared with the generator: one seed governs generation, resampling and noise.

```{include} generated/input-keys.md
:start-after: <!-- begin keys-beam.resample -->
:end-before: <!-- end keys-beam.resample -->
```

On the user's path the counts are `beam.n_particle` and `beamlet_size`. Genesis4's `match` and `center` are not ported: a Bmad lattice carries its Twiss and `beam_init` generates matched bunches already. The method is in [](fel-physics.md#sec-import).

### Check instruments

These are not physics input. The benchmark sets them.


(param-beam-split-weights)=
**`split_weights`** exists because no Genesis4 comparison can test the weighted paths: its dumps carry no weights, so every cross-code comparison sees the uniform case. Every collective observable must be identical to the unsplit run, which makes a bug like using one particle's weight for all visible. The measured level is the `weight_split` tier in [](validation.md).

(param-beam-swap-beam-xy)=
**`swap_beam_xy`** feeds the rotation identity of the two-polarization checks: an all-y line fed the swapped beam must reproduce the all-x line.

(group-fel-wavefront-init)=
## `[field]` and `[window]`

```toml
[field]
lambda0 = 1e-10
seed_power = 5e3
seed_waist_size = 30e-6
grid_half_width = 2e-4

[window]
n_wavelength = 3
```

```{include} generated/input-keys.md
:start-after: <!-- begin keys-field -->
:end-before: <!-- end keys-field -->
```


| `field.grid_n_pts` | `0` | Transverse grid points per side. Zero derives it |
| `field.grid_half_width` | `0` | Transverse grid half width [m]. Zero derives it |
| `field.seed_power` | `0` | Gaussian seed power [W]. Zero is a dark start |
| `field.seed_waist_size` | `0` | Seed intensity 1/e^2 radius [m] |
| `field.seed_polarization` | `"x"` | `"x"` or `"y"`. A field from `field.file` decides its own planes |
| `field.harmonics` | `1` | The field set. The first entry must be the fundamental |

(param-wavefront-grid)=
**`field.grid_n_pts`** and **`field.grid_half_width`** are derived from the beam when the deck leaves them at zero, and both are printed with their origin at setup. The rms beam size comes from the emittances the deck states and the matched Twiss the lattice states, averaged over the FEL elements by length. The half width is nine of those, which contains the mode, and the cells are a seventh of one, which resolves the beam and the mode. Both constants were measured on three machines a hundred times apart in wavelength ([](startup-noise.md#sn-recommendations)). The two together fix the point count at 127 before rounding, since the beam size cancels. That count is then rounded up to the next power of two, 128, so that one grid serves the CPU and the device field solver alike and a derived deck moves between the backends unchanged. Rounding up only ever refines, more points across one half width being smaller cells, and setup prints the rounded count beside the count it came from.

One may be set and the other derived against it. A stated half width keeps the cell rule, so the point count follows from it. A stated point count keeps the containment rule, so the half width is nine beam sizes and the cells are whatever the count makes them. Nothing is derived without an emittance to derive from: a run that loads its beam from a dump states its own grid, or takes the field's from `field.file`. A field from `field.file` is read with its own `unitSI` and `gridUnitSI` applied, and its `y` component must be shaped as its `x` is. The source filter's angle is then converted on the file's grid rather than on any grid the deck states.

(param-wavefront-lambda0)=
**`field.lambda0`** is required. It is deliberately not defaulted from the lattice resonance, since the first undulator may be detuned. Starting from a beam dump it is required too. The file carries the slice partition and not the radiation wavelength it was sliced on, and `window.n_wavelength` is required for the same reason.

(param-slicing)=
```{include} generated/input-keys.md
:start-after: <!-- begin keys-window -->
:end-before: <!-- end keys-window -->
```

**`[window]`** is the longitudinal discretization of the run, and it belongs to neither the beam nor the field: the beam's slices and the field's longitudinal samples are one partition, stated once. A wavefront carries only the spacing that follows, since a wavefront travelling through mirrors and gratings has no slices.

`window.n_wavelength` is the spacing in wavelengths, Genesis4's `sample`. It is an integer because slippage rotates the field's slice ring by one index per slice, exactly and without interpolation, which holds only for a whole number. It is carried as an integer from the deck to every consumer rather than recovered by dividing a spacing in metres by a wavelength. A negative value is refused. Zero derives it as the whole number of wavelengths in a quarter of the cooperation length $\lambda/(4\pi\rho)$, four slices per cooperation length, which is where FLASH1's pulse energy converged. The Pierce parameter is taken from the beam the deck describes and the first FEL element. In the steady state it stays one, because the whole charge sits in the one slice and the spacing there states the current rather than a resolution. A beam or field read from a dump is the one place the spacing is recovered from a length, since no dump format records the integer. A spacing that is not a whole number of wavelengths is refused there.

`window.window_length` when zero holds every particle of the bunch and the line's slippage, one wavelength per undulator period, ahead of the window head. Radiation that slips forward then has somewhere to go before it leaves. A drawn bunch spans its last particle to its first. A Gaussian has no last particle, so it spans the length beyond which a slice would hold less than one electron of the charge. That rule scales with the charge where a fixed number of sigmas does not. A bunch too weak to put one electron in its peak slice spans eight sigmas. A flat bunch takes its grid extent, and the steady state is one slice. A stated window says its own headroom and gets none added, and one that clips the bunch warns. `window.n_slice` states the same window as a count instead, and setting both is refused.

`window.current` states a flat current directly, which is Genesis4's `&beam current`. The bunch's z structure is then not read. With a window it is a flat time-dependent run, and with none it is the steady state of one slice, so no separate steady-state switch exists. Unset, the current comes from `beam.bunch_charge` and the bunch's own z distribution. Setting both is refused.

(param-wavefront-harmonics)=
**`field.harmonics`** is a gap-free increasing list whose first entry must be `1`: the fundamental anchors the optical phase, the reference advance and the slippage schedule, and the harmonics ride on it. Anything else is refused. Nothing is ever summed across harmonics. See [](fel-physics.md#sec-field-set).

(lattice-attributes)=
## Lattice attributes

An FEL segment is a wiggler or undulator carrying Bmad's `fel_method` attribute, and the method is what selects the physics:

| `fel_method` | Meaning |
|---|---|
| `off` | The default. Not an FEL segment, and this program leaves it to the seam |
| `averaged` | The wiggle-averaged model. The production workhorse |
| `unaveraged` | Direct integration through the undulator field, with no period averaging |

It is Bmad's own attribute rather than anything this program registers, so `show ele` prints it, a written lattice keeps it, the parser refuses a misspelling, and every other Bmad program loads the lattice. It is class-settable as any attribute is, `wiggler::*[FEL_METHOD] = unaveraged`, and the two settings mix freely in one line. `tracking_method` is a separate axis and names how one particle crosses the element, as it does everywhere in Bmad. An FEL segment leaves it at `bmad_standard` and is the plain periodic wiggler its field attributes describe outside the FEL walk.

Space charge is per element too, through Bmad's own `space_charge_method` attribute rather than anything this program registers:

| `space_charge_method` | Meaning |
|---|---|
| `off` | The default. No space charge in this element |
| `slice` | The slice-binned longitudinal solve, with the FEL slices as the bins |

`fft_3d` and `cathode_fft_3d` are refused on an FEL element, since their solvers want a three-dimensional grid this walk does not build. On an element that is not an FEL element, any method other than `off`, and any `csr_method` other than `off`, is refused while the master switch is on and the interlude model is `bmad`. The seam tracks one slice at a time through Bmad's `track1_bunch`, whose collective paths need a centroid orbit this program does not carry. Bmad would return the bunch untracked with no error. The transcribed Genesis interlude applies the slice solve itself, so with `interlude_model = "genesis"` a `slice` interlude is tracked. Bmad's master switch applies as it does everywhere else: `bmad.csr_and_space_charge_on` must also be true, and when elements ask for `slice` while it is false the run says so and tracks without space charge. Bmad's `space_charge_com` `n_bin`, set in the lattice file is ignored, because the slices are the bins. Inside the Bmad seam the same attribute drives Bmad's own machinery, so one lattice reads the same way in every Bmad program.

The unaveraged mode's two numbers are run switches, [](#param-global-unaveraged-steps-per-period) and [](#param-global-unaveraged-ramp-periods), rather than per-element attributes.

(attr-fel-tracking)=
**`unaveraged`** is a full Newton-Lorentz quiver with no period averaging and no coupling factor anywhere in its inputs. The run writes an energy ledger beside its other outputs. **`averaged`** is the wiggle-averaged model, whose transverse maps are chosen by [](#param-global-transport-model). See [](fel-physics.md#sec-unaveraged).

The FEL parameters themselves come from the element: `aw` from `b_max` and `l_period`, the helicity from `field_calc`, and the step from `ds_step`. Each is asserted at setup. See [](fel-physics.md#sec-element).

## Outputs

Every output file, and how to read it, is [](reading-output.md). In brief: `<out_root>.stats.h5` is the production statistics file, and `<out_root>.diag.txt` is the text comparison instrument written only under `write_diag`. The end state is dumped as openPMD. The escaped-field and pulse files are written beside it in Genesis4's field format, and [](genesis4.md) covers the exchange.
