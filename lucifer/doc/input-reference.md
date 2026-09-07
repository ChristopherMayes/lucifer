# Input reference

The normative reference for every namelist parameter Lucifer honors: its default, its meaning, and what refuses it. A run is one namelist file of three groups, each setting Bmad or program structures directly, which is Tao's `&tao_params` pattern. Defaults live in the struct declarations, and a harness check holds this page to them in both directions, so a stated default cannot drift from the real one and a new parameter cannot go undocumented.

| Group | Carries |
|---|---|
| [`&fel_params`](#group-fel-params) | The run: the lattice, the `global%` switches, Bmad's `bmad_com` and `space_charge_com`, and the `chamber_wake%` and `space_charge%` collective descriptions |
| [`&fel_beam_init`](#group-fel-beam-init) | The beam: Bmad's `beam_init%` description, the load mode and its two switches, the `resample%` resampler, files, and the check instruments |
| [`&fel_wavefront_init`](#group-fel-wavefront-init) | The radiation: the `wavefront_init%` starting condition and `field_file`. The field record is the time window, so the window lives here |

Per-element settings are [lattice attributes](#lattice-attributes) rather than namelist parameters. How to build and run is [](user-guide.md). What the outputs hold is [](reading-output.md). The measured levels are [](validation.md). The physics is the manual, [](fel-physics.md). Every resolved input is written into the statistics file's `params/` group, so a finished run states its own configuration with every default explicit.

A flat `&fel_track_params` group is refused, with each parameter mapped to the group that now carries it.

(group-fel-params)=
## &fel_params

```
&fel_params
  lat_file = "aramis.bmad"
  global%out_root = "fel_td"
/
```

### The lattice and the run

| Parameter | Default | Meaning |
|---|---|---|
| `lat_file` | `""` | The Bmad lattice file |
| `global%out_root` | `"fel_track"` | Prefix for every output file |
| `global%interlude_model` | `"bmad"` | How field-free elements are tracked: `"bmad"` (the seam) or `"genesis"` (transcribed) |
| `global%transport_model` | `"bmad"` | Transverse transport inside averaged FEL elements: `"bmad"` (Bmad's own kernel) or `"genesis"` (transcribed, validation-internal) |
| `global%unaveraged_steps_per_period` | `20` | Integration steps per undulator period in unaveraged elements ([](#param-global-unaveraged-steps-per-period)) |
| `global%unaveraged_ramp_periods` | `2` | Length of the sin^2 entry and exit ramps, in periods ([](#param-global-unaveraged-ramp-periods)) |
| `slicing%window_length` | `0` | Time window [m]. Zero holds every particle and the line's slippage, or is one slice ([](#param-slicing)) |
| `slicing%n_wavelength` | `0` | Slice spacing in wavelengths, an integer. Zero derives it from the gain ([](#param-slicing)) |
| `slicing%n_slice` | `0` | The window as a slice count instead. Both it and `window_length` is refused |
| `slicing%current` | `0` | A flat current [A], which needs no z description. Zero takes the current from the bunch |
| `global%source_model` | `"deposit"` | The FEL source: `"deposit"` per particle, or `"coherent"` for the coherent retrieval |
| `global%source_filter` | `T` | Filter the source term at wide transverse angles ([](#param-global-source-filter)) |
| `global%source_filter_angle` | `0` | The filter's sigmoid edge as an angle [rad]. 0 derives it from the beam and the gain |
| `global%source_filter_xcut` | `0` | Validation-internal: the edge in x in units of half the grid's Nyquist frequency |
| `global%source_filter_ycut` | `0` | Validation-internal: the edge in y, same units |
| `global%source_filter_width` | `0.05` | The filter's sigmoid width as a fraction of its edge. Small is a sharp edge |
| `global%source_filter_mutate` | `F` | The filter check's self-test: filter the field instead of the source |
| `global%track_start` | `""` | Element locator bounding the walk below. Blank is the whole line |
| `global%track_end` | `""` | Element locator bounding the walk above. Blank is the whole line |
| `global%ran_seed` | `12345` | The one random seed, governing generation, resampling and noise |
| `global%reference_run` | `F` | Permit a lattice with no FEL element, so Bmad tracks every element |
| `global%load_only` | `F` | The pre-run: build the initial state, print the header, dump it, exit without tracking ([](#param-global-load-only)) |
| `global%migrate` | `F` | Move particles between slices when their ponderomotive phase leaves the slice |
| `global%migrate_check` | `F` | Verify phase continuity at every migration and report the worst deviation |
| `global%fp32_check` | `"off"` | Step a single-precision twin beside the FP64 path and record its divergence: `"lockstep"` or `"freerun"` |
| `global%fp32_mutate` | `F` | The instrument's self-test: coarsen the FP32 residual so the recorded level must move |
| `global%device` | `"off"` | Run the averaged FEL step on a device backend: `"metal"` on a build that carries it |

(param-global-interlude-model)=
**`global%interlude_model`** selects how the field-free elements are handled. `"bmad"` is the seam: Bmad's own `track1_bunch`, an exact theta mapping, and `wavefront_drift` for the field. `"genesis"` uses the transcribed interlude step everywhere, which prices what the seam changes. The slippage schedule is identical in both models. The physics is in the manual's interlude and seam sections.

(param-global-transport-model)=
**`global%transport_model`** selects the transverse maps inside an averaged FEL element. `"bmad"`, the default and the production model, is Bmad's own periodic-wiggler kernel with its end-field treatment and chromaticity through `p0/p`. `"genesis"` is Genesis4's maps verbatim, with the focusing matrix built from `aw`, `kx`, `ky` and chromaticity through `gammaz`. The second is validation-internal: a transcription comparison needs transcription-level transport, so the comparison tiers set it and no production run does. The two are priced against each other in the manual's element section, and the unaveraged mode integrates the field rather than applying a map, so the switch does not reach it.

(param-global-unaveraged-steps-per-period)=
**`global%unaveraged_steps_per_period`** is 20 by default. Below 10 is refused, the floor our own coupling-factor convergence supports: [](validation.md) tabulates it at 10, 20 and 30 steps per period. It reaches unaveraged elements only, and averaged elements take their step from `ds_step` as before.

(param-global-unaveraged-ramp-periods)=
**`global%unaveraged_ramp_periods`** is 2 by default. A true hard edge is a test configuration and has the explicit sentinel `-1`, so silence never means hard edge. A ramp pair longer than the segment is refused. It was a per-element attribute this program registered, which every other Bmad program then refused to parse, so it is a run switch selected uniformly. Averaged and unaveraged elements still mix per element, through `fel_method`.

(param-global-source-model)=
**`global%source_model`** is `"deposit"`, the standard per-particle scatter, or `"coherent"`, the coherent-Gaussian retrieval in which the spatially incoherent part of the source is dropped, the slice bunch factor keeps the physical shot noise, and the transverse shape is a guarded Gaussian. A profile that is measurably not Gaussian is refused. See [](fel-physics.md#sec-coherent-source).

(param-global-source-filter)=
**`global%source_filter`** filters the source term at wide transverse angles. Each beamlet is a single transverse point on the deposition grid and the source has no angular dependence, so a beamlet radiates into every transverse wavenumber the grid carries, well outside the undulator's central cone. On, the transformed source is multiplied by a sigmoid in normalized spatial frequency before it is added to the field, and the field itself is never multiplied: a seed, an imported field and everything already radiated propagate exactly as they do with the filter off. `source_filter_angle` places the sigmoid's half-height, in radians. Left at zero the run derives it as the larger of four diffraction angles of the fundamental mode, `4*lambda0/(2*pi*sigma)` with `sigma` the rms transverse size of the loaded beam, and the angle at which the resonant wavelength red-shifts by the Pierce parameter, `sqrt(2*rho*lambda0/l_period)`. Inside the first lies the radiation that can couple to the mode; outside the second, emission is beyond the bandwidth of the interaction. Their ratio goes as `6.8/sqrt(z_R/L_g)`, measured to 5 percent over four machines, so the mode angle is the larger below a `z_R/L_g` of about 46 and neither alone serves for every machine. The larger is taken because cutting into the mode loses real radiation while leaving artifact in only weakens the filter. Both angles, their ratio and which one won are printed at setup, with a warning when the ratio leaves 0.3 to 3. `source_filter_width` is the sigmoid's width as a fraction of that edge, so the sigmoid on axis is `1/(1 + exp(-1/width))`: the default of 0.05 leaves the axis untouched, and a width that puts the axis below 0.99 is refused. `source_filter_xcut` and `source_filter_ycut` place the edge in Genesis4's own units instead, the shifted grid index over the grid size, so 1 is the angle `lambda0/dx`, twice the Nyquist angle and off the grid, and the value moves with the grid. They exist so the transcription check can place the edge where Genesis4 places it ([](validation.md#val-source-filter)), they lift the width guard with them, and setting them together with the angle is refused. On by default, since what it removes is emission no physical beam radiates and it is the cheaper route to a converged answer ([](fel-physics.md#sec-convergence)). Turn it off to compare against a code that carries no such filter or against a level recorded without one, which is what the comparison tiers do. A run with it off is bit-for-bit the run before the filter existed. The filter reaches the averaged elements only, the unaveraged mode depositing from the resolved motion and building no transformed source to multiply, so a mixed line filters its averaged elements and leaves the rest, with a line at setup saying so. `source_model = "coherent"` turns the filter off for the run and says so at setup, since the coherent source is an analytic Gaussian that carries no wide-angle content to remove. Neither is refused, because the filter is now the default and a deck that never asked for it must still run. A non-positive cut or width is refused rather than quietly turning the filter off. What the filter is for, what it costs and what load it wants is [](fel-physics.md#sec-convergence), and the measurements behind it are [](startup-noise.md). `source_filter_mutate` applies the sigmoid to the propagated field instead of to the source, which is how the harness proves the check can fail.

(param-global-load-only)=
**`global%load_only`** is the pre-run. The run parses the lattice, derives everything it derives, builds the beam and the field, prints the header and stops, writing the initial state as `<out_root>-initial.beam.h5` and `<out_root>-initial.wf.h5` through the same writers the final dumps use, so the tracker reads them straight back. The header it prints is the one a full run prints, which is the point: the grid, the slice spacing and the step with their origins, the work the run will do, the wall clock that work is estimated to cost on the machine the estimate names, and where the macroparticle count stands against the load the source filter's convergence was measured at ([](user-guide.md), [](fel-physics.md#sec-convergence)). A deck can therefore be priced and its load read before any of the tracking is paid for. `global%write_initial` writes the same two files and then tracks.

(param-global-track-start)=
**`global%track_start`** and **`global%track_end`** bound the walk, using Bmad's element-locator syntax. The schedule (slippage, autophasing, break geometry) is always built on the full lattice, so a windowed run composes exactly with the full one: the first span followed by the second, started from the first's dumps, reproduces the one-shot run. The measured level is in [](validation.md#val-the-programs-own-identities).

(param-global-fp32-check)=
**`global%fp32_check`** arms the FP32 lockstep instrument, a check rather than a production mode: an FP32 twin of the averaged FEL advance steps beside the FP64 path from a shared state, and the per-quantity divergence goes to `<out_root>.fp32.txt` with the worst case in the footer. `"lockstep"` rebuilds the FP32 state from FP64 every step, so a wrong formula shows as a jump. `"freerun"` carries the FP32 longitudinal state and the FP32 field across steps, a complete single-precision run beside the FP64 one, so the compounding rate and the end-to-end observable divergence are measured; it covers a single-slice window and refuses more, since the twin keeps no slippage rotation of its own. The twin covers the fundamental-only, collective-free averaged advance and refuses anything else (harmonics, two polarizations, the coherent source, wakes, space charge, the unaveraged mode). With the device in the twin's role the harmonic and two-polarization refusals lift, since the device twin covers the set. A run whose FP32 residual cannot resolve the per-step phase increment is refused by the runtime guard rather than reported as a clean number. The FP64 physics is untouched: the instrumented run's outputs are byte-identical to the uninstrumented run's. The measured levels are in [](validation.md#val-fp32-lockstep). `fp32_mutate` coarsens the FP32 residual by eight mantissa bits so the harness can prove the check fails when the path is wrong. The instrument needs the single-precision FFTW (`libfftw3f`), since the field twin's transform is FFTW's own `fftwf` interface; a build without that library refuses `fp32_check` rather than transforming in another precision, and the build reports which way it went at configure time.

(param-global-device)=
**`global%device`** selects a device backend for the averaged FEL step. `"metal"` runs it on an Apple Silicon GPU in single precision: the beam and the field upload at each FEL element's entry and stay resident through the element, every integration step is one command buffer (transverse maps, the longitudinal push in a fixed-point phase chart, the source deposit and the FFT field solve), and the state returns to the host at the stats comb's positions and the element end. Between those positions the host arrays are stale by design, so `comb_ds_save` is also the readback knob. Combined with `fp32_check` the device instead takes the lockstep twin's role: the FP64 run is untouched and the device's per-step divergence goes to the same `.fp32.txt` stream, which is how the backend is judged ([](validation.md#val-device)). The field set rides the device whole: every harmonic member and both polarization planes, though not the two together, which is refused for the CPU as well. Slice migration runs with the device: the pass reads the beam back at the element's last step and re-slices it on the host, the next element uploads the result, and a slice that outgrows the device's particle capacity grows it, which the run reports. Everything else the kernels do not cover is refused at setup or first use (the coherent source, wakes, space charge, spontaneous radiation, the unaveraged mode, the escaped-field bank, a transverse grid that is not a power of two from 64 to 1024), and a build without the backend refuses the knob itself: an unsupported configuration stops the run and never quietly takes the CPU path. The backend is built where the toolchain can carry it, which is macOS with a Clang-family Objective-C++ compiler rather than macOS alone, since it is ARC-managed Objective-C++ against the Metal framework. Every other build takes the refusing stub and says so at configure time. Running `lucifer` with no arguments names what the build in hand carries, either the device it found or the reason there is none, so the question needs no run and no deck to answer.

(param-global-migrate)=
**`global%migrate`** is off by default, and the reason is the comparison rather than the physics: the tiers that compare against Genesis 1.3 Version 4 (Genesis4) run against a code that never migrates, so migration inside a transcription-level comparison would be a model difference. Dropped charge is counted and reported per event. With `migrate_check = T` the run also verifies exact phase continuity at every migration. See [](fel-physics.md#sec-migration).

### Output switches

| Parameter | Default | Meaning |
|---|---|---|
| `global%write_diag` | `F` | Write the text comparison diagnostics, one row per slice per record. Large |
| `global%write_initial` | `F` | Also dump the initial state |
| `global%dump_beam_at` | `""` | Element locators for mid-run beam dumps |
| `global%dump_field_at` | `""` | Element locators for mid-run field dumps |
| `global%keep_escaped_field` | `F` | Bank the field slices slippage carries out of the window, and rebuild the full pulse at exit |
| `global%comb_ds_save` | `0` | Minimum z advance between per-record statistics rows |
| `global%record_environment` | `F` | Also record the user name and working directory in the statistics file |

(param-global-dump-beam-at)=
**`global%dump_beam_at`** and **`global%dump_field_at`** name elements through Bmad's own locator, so `class::name` syntax works. An entry matching no element is refused. Dumps are openPMD, and the field dump is unrotated into time order first.

(param-global-keep-escaped-field)=
**`global%keep_escaped_field`** banks each field slice that slippage transmits beyond the window, with its `wavefront_params` and transmission position, and at finalize propagates each to the exit plane to write the whole pulse. Field that has left the window never re-interacts, so it is fixed information. See [](fel-physics.md#sec-stats).

(param-global-comb-ds-save)=
**`global%comb_ds_save`** is Bmad's `bunch_track_struct%ds_save` name and semantics. Negative keeps no per-record rows at all, while element ends, dumps and the final state remain. Zero, the default here, records every position. Positive records a row once z has advanced that far past the last one, with element ends always kept. The default deviates from Tao's, deliberately: the per-record arrays are this program's primary statistics contract where Tao's comb is an optional extra.

(param-global-record-environment)=
**`global%record_environment`** is off by default because a statistics file is meant to travel, attached to a paper or mailed to a collaborator, and the user name and working directory identify a person and a machine. The timestamp and the Bmad version identify the run without them. See [](fel-physics.md#sec-meta).

### Bmad's own structures

`bmad_com` and `space_charge_com` are set directly in this group, as Tao sets them, and their defaults are Bmad's. Two matter most here.

| Parameter | Meaning |
|---|---|
| `bmad_com%radiation_damping_on` | Spontaneous energy loss in the FEL step and through Bmad's elements |
| `bmad_com%radiation_fluctuations_on` | Quantum diffusion, one draw per beamlet |

Both are off by default, matching Genesis4's `&sponrad`. Fluctuations with `global%migrate = T` are refused: the quiet start cancels per beamlet, and migration scrambles the grouping. See [](fel-physics.md#sec-eom).

### Chamber wakes: `chamber_wake%`

Genesis4's `&wake` names. Bmad element wakes are a separate mechanism, described under [element wakes](#element-wakes) below.

| Parameter | Default | Meaning |
|---|---|---|
| `chamber_wake%on` | `F` | Enable the chamber wake |
| `chamber_wake%model` | `"genesis"` | Which implementation runs. `"genesis"` is the only accepted value |
| `chamber_wake%loss` | `0` | External loss [eV/m] |
| `chamber_wake%radius` | `2.5e-3` | Chamber radius, or half gap if flat [m] |
| `chamber_wake%conductivity` | `0` | DC conductivity [1/(Ohm m)]. Zero means no resistive wake |
| `chamber_wake%relaxation` | `0` | AC relaxation distance c*tau [m] |
| `chamber_wake%roundpipe` | `T` | Round chamber. `F` is flat, parallel plates |
| `chamber_wake%material` | `""` | `"CU"` or `"AL"` shortcut for conductivity and relaxation |
| `chamber_wake%gap` | `0` | Undulator gap [m]. Zero means no geometric wake |
| `chamber_wake%lgap` | `1` | Period of the gaps [m] |
| `chamber_wake%hrough` | `0` | Roughness amplitude [m]. Zero means no roughness wake |
| `chamber_wake%lrough` | `1` | Roughness period [m] |
| `chamber_wake%write_kernels` | `""` | Export the transcribed kernels to this file, for building `z_long` tables |

(param-write-wake-kernels)=
**`chamber_wake%write_kernels`** is a bare name in this group rather than a `chamber_wake%` component, and it is the only name that works: it lands in `chamber_wake%write_kernels`, which the parser assigns unconditionally, so a value written as `chamber_wake%write_kernels` is overwritten. The export builds matching `z_long` tables, and `examples/bmad_wake` uses one.

(param-chamber-wake-model)=
**`chamber_wake%model`** names the implementation, and `"genesis"` is the only value accepted today: the transcribed solver, which convolves its kernels with the weighted slice currents and produces one energy loss per slice. Anything else is refused rather than treated as the transcribed solver by default. The field exists so that a second implementation can arrive as a value here rather than as a rework.

The three kernels, their numerical impedance and the causal convolution are in [](fel-physics.md#sec-wakes).

### Space charge: `space_charge%`

Whether space charge acts is not in this group. It is the element's own
`space_charge_method`, described with the other [lattice attributes](#lattice-attributes),
and this group holds only the solver's numbers. A run with no element asking for it never
reaches them.

| Parameter | Default | Meaning |
|---|---|---|
| `space_charge%model` | `"genesis"` | Which implementation runs. `"genesis"` is the only accepted value |
| `space_charge%rmax` | `0` | Radial grid extent scale [m]. Grows adaptively |
| `space_charge%ngrid` | `100` | Radial grid points |
| `space_charge%nz` | `0` | Longitudinal harmonics. Zero disables the short-range solve |
| `space_charge%nphi` | `0` | Azimuthal modes, m over -nphi to nphi |
| `space_charge%longrange` | `F` | The whole-window long-range term |

(param-space-charge-model)=
**`space_charge%model`** names the implementation, and `"genesis"` is the only value accepted today: the transcribed solver, which works per slice on a radial grid over azimuthal modes and longitudinal harmonics. Anything else is refused. Note that this family is Lucifer's own, and `space_charge_com` in the same group is Bmad's global structure, a separate thing.

The default is expected to stay `"genesis"` when Bmad's own slice solver becomes the second value, and the reason is the physics rather than the order they arrived in. The transcribed solver carries two terms that matter inside an undulator and that Bmad's slice model does not have: the space charge of the microbunching itself, solved per longitudinal harmonic of the ponderomotive phase, and the longitudinal Lorentz factor, which at aw = 0.85 rms is worth a factor of 1.7 in the field's own scaling. What Bmad's model adds that this one lacks, a transverse defocusing kick, falls as the inverse cube of gamma and is an injector term rather than an undulator-line one. The transcribed path is also the measured one, at the level [](validation.md) records for the space-charge tier. A second value earns the default by measurement, not by being newer.

If neither term is configured, `nz = 0` with `longrange = F`, an element asking for `slice` is refused: the solve would cost its full price and return an exact zero, and which of the two terms was meant is worth asking.

See [](fel-physics.md#sec-spacecharge).

(element-wakes)=
### Element wakes

Elements carrying Bmad `sr_wake` definitions, either pseudomodes or a tabular `z_long`, act across the whole time window through slice concatenation. The conventions, the step size, the mid-element wiggler kick and the refusals are in [](fel-physics.md#sec-seamwake). Wakes and space charge are both refused in the unaveraged mode, which does not wire them into its step.

(group-fel-beam-init)=
## &fel_beam_init

```
&fel_beam_init
  beam_init%n_particle = 2048
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beam_init%sig_z = 1.2e-9
  beam_init%sig_pz = 8.8e-5
  beam_init%bunch_charge = 3.0e-14
  shot_noise = T
/
```

One path loads the beam: a bunch from `beam_init%`, generated or read from a file, binned
into the slices and loaded by `load_mode`, with `quiet_start` and `shot_noise` acting on
every load ([](fel-physics.md#sec-loading)). `beam_file` starts from a dump instead.

### Files and sources

| Parameter | Default | Meaning |
|---|---|---|
| `beam_file` | `""` | openPMD particle dump to start from, slices as they are. Blank loads a bunch instead |
| `write_genesis_dist` | `""` | Write the bunch as a Genesis4 `&importdistribution` input |
| `write_openpmd_file` | `""` | Write the bunch as openPMD-beamphysics |

(param-beam-beam-file)=
**`beam_file`** reads openPMD only. A file that is not openPMD is refused, and the message carries the conversion command. `tests/scripts/convert_genesis.py` converts in both directions, and [](genesis4.md) describes the exchange. A dump named in `beam_init%position_file` is refused, since its slices would collapse onto one another when read as a bunch.

### The bunch description: `beam_init%`

Bmad's standard `beam_init_struct`. With `position_file` blank the bunch is generated from the description, and in the sample mode the loader evaluates the description analytically per slice without making a bunch at all. The Twiss is always the lattice's.

| Parameter | Meaning |
|---|---|
| `beam_init%position_file` | openPMD-beamphysics bunch to load. Read whole, then sliced by `load_mode` |
| `beam_init%n_particle` | Sample mode: macroparticles per slice after the phase copies, a positive multiple of `beamlet_size`. Keep mode: the particles of the generated bunch |
| `beam_init%a_norm_emit`, `beam_init%b_norm_emit` | Normalized emittances [m rad] |
| `beam_init%sig_pz` | Fractional momentum spread dP/P0 |
| `beam_init%bunch_charge` | Charge [C]. The current is derived from it, never input |
| `beam_init%sig_z` | Bunch length [m], read with `distribution_type(3)` |
| `beam_init%distribution_type(3)` | `"RAN_GAUSS"` for a Gaussian current profile, `"GRID"` for Bmad's uniform one |

(param-beam-init-contract)=
The analytic loader honors the fields above and refuses every other `beam_init` field that is set. A standard structure that silently dropped fields would be worse than a custom one, so the honored set is the contract, and `check_beam_init_contract` enforces it. A bunch that Bmad makes, in keep mode or from a file, honors everything `init_beam_distribution` honors.

(param-beam-sig-z)=
**`beam_init%sig_z`** with `"RAN_GAUSS"` gives a Gaussian current profile evaluated at the slice centers, with the bunch centered in the window. A zero length is the steady state, one slice holding the whole charge, and it is refused for a time-dependent window. With `"GRID"` the profile is flat over the z extent of `grid(3)`. `beam_init%a_emit` and `b_emit` are refused: normalized emittances only, which is Bmad's preferred form. `sig_e` is deprecated Bmad-wide and does not exist here.

### The load

| Parameter | Default | Meaning |
|---|---|---|
| `load_mode` | `"sample"` | `"sample"`: the same count in every slice, weights from the slice charge. `"keep"`: every particle kept, counts following the charge |
| `quiet_start` | `T` | Load beamlets of `beamlet_size` phase copies, quiet below that harmonic |
| `beamlet_size` | `8` | Copies per beamlet, and the harmonics the noise resolves |
| `shot_noise` | `F` | Impose the physical noise on the load's groups. Time-dependent windows only |

(param-beam-load-mode)=
**`load_mode`** is what becomes of the particles in a slice. In the sample mode a generated bunch takes the analytic loader and a file takes the resampler ([](fel-physics.md#sec-import)); `beam_init%n_particle` is the count per slice either way. In the keep mode each particle becomes a beamlet of copies at weight over `beamlet_size` sharing its coordinates, so every moment of the bunch is the load's, slice by slice, and the window is the bunch's extent unless `slicing%` states one. The method is in [](fel-physics.md#sec-loading).

(param-beam-shot-noise)=
**`shot_noise`** imposes the weighted Fawley loading on the groups the load has: the beamlets when the loader made them, else particles sharing their transverse coordinates, else the occupants of a deposit cell, and the message says which. The loader measures the quiet floor first and refuses noise on a load above it with the value in the message, so noise is never counted twice. The algebra, the rule and the four switch combinations are in [](fel-physics.md#sec-noise). The loader warns where Genesis4 silently clamps groups holding fewer real electrons than macroparticles.

### The resampler: `resample%`

Named after Genesis4's `&importdistribution` where an equivalent exists. `slicing%n_wavelength`, `ran_seed` and the seed field are shared with the generator: one seed governs generation, resampling and noise.

| Parameter | Default | Meaning |
|---|---|---|
| `resample%slice_width` | `0.01` | Sampling window over bunch length |
| `resample%n_slice` | `0` | Slice count. Zero derives it from the bunch length and the spacing |
| `resample%use_beam_init` | `F` | Validation route: resample a bunch Bmad generated from `beam_init%` instead of taking the analytic loader |
| `resample%n_particle_per_slice` | `8192` | Macroparticles per slice on the validation route |
| `resample%beamlet_size` | `4` | Beamlet size on the validation route |

On the user's path the counts are `beam_init%n_particle` and `beamlet_size`. Genesis4's `match` and `center` are not ported: a Bmad lattice carries its Twiss and `beam_init` generates matched bunches already. The method is in [](fel-physics.md#sec-import).

### Check instruments

Not physics input. The validation harness sets these.

| Parameter | Default | Meaning |
|---|---|---|
| `split_weights` | `F` | Replace each particle by coincident copies carrying 1/3 and 2/3 of its weight |
| `resample_split_weights` | `F` | The same split, applied to the bunch before it is sliced |
| `gen_test_weights` | `F` | Alternate beamlet weights of 0.25x and 1.75x, charge preserving, to exercise the weighted-noise paths |
| `swap_beam_xy` | `F` | Swap (x, px) with (y, py) after generation |

(param-beam-split-weights)=
**`split_weights`** exists because no Genesis4 comparison can test the weighted paths: its dumps carry no weights, so every cross-code comparison sees the uniform case. Every collective observable must be identical to the unsplit run, which makes a bug like using one particle's weight for all visible. The measured level is the `weight_split` tier in [](validation.md).

(param-beam-swap-beam-xy)=
**`swap_beam_xy`** feeds the rotation identity of the two-polarization checks: an all-y line fed the swapped beam must reproduce the all-x line.

(group-fel-wavefront-init)=
## &fel_wavefront_init

```
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  slicing%n_wavelength = 3
  wavefront_init%seed_power = 5e3
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_half_width = 2e-4
/
```

| Parameter | Default | Meaning |
|---|---|---|
| `field_file` | `""` | openPMD EXT_Wavefront field dump to start from |
| `wavefront_init%lambda0` | `0` | Radiation wavelength [m]. Required for generation |

| `wavefront_init%grid_n_pts` | `0` | Transverse grid points per side. Zero derives it |
| `wavefront_init%grid_half_width` | `0` | Transverse grid half width [m]. Zero derives it |
| `wavefront_init%seed_power` | `0` | Gaussian seed power [W]. Zero is a dark start |
| `wavefront_init%seed_waist_size` | `0` | Seed intensity 1/e^2 radius [m] |
| `wavefront_init%seed_polarization` | `"x"` | `"x"` or `"y"` |
| `wavefront_init%harmonics` | `1` | The field set. The first entry must be the fundamental |

(param-wavefront-grid)=
**`wavefront_init%grid_n_pts`** and **`wavefront_init%grid_half_width`** are derived from the beam when the deck leaves them at zero, and both are printed with their origin at setup. The rms beam size comes from the emittances the deck states and the matched Twiss the lattice states, averaged over the FEL elements by length. The half width is nine of those, which contains the mode, and the cells are a seventh of one, which resolves the beam and the mode. Both constants were measured on three machines a hundred times apart in wavelength ([](startup-noise.md#sn-recommendations)), and the two together fix the point count at 127, since the beam size cancels. On the device the count is rounded up to the power of two its field solver takes.

One may be set and the other derived against it. A stated half width keeps the cell rule, so the point count follows from it. A stated point count keeps the containment rule, so the half width is nine beam sizes and the cells are whatever the count makes them. Nothing is derived without an emittance to derive from: a run that loads its beam from a dump states its own grid, or takes the field's from `field_file`.

(param-wavefront-lambda0)=
**`wavefront_init%lambda0`** is required, and deliberately not defaulted from the lattice resonance, since the first undulator may be detuned. Starting from a beam dump it is required too: the file carries the slice partition but not the radiation wavelength it was sliced on, and `slicing%n_wavelength` is required for the same reason.

(param-slicing)=
**`slicing%`** is the longitudinal discretization of the run, and it belongs to neither the beam nor the field: the beam's slices and the field's longitudinal samples are one partition, stated once. A wavefront carries only the spacing that follows, since a wavefront travelling through mirrors and gratings has no slices.

`slicing%n_wavelength` is the spacing in wavelengths, Genesis4's `sample`. It is an integer because slippage rotates the field's slice ring by one index per slice, exactly and without interpolation, which holds only for a whole number, and it is carried as an integer from the deck to every consumer rather than recovered by dividing a spacing in metres by a wavelength. Zero derives it as the whole number of wavelengths in a quarter of the cooperation length $\lambda/(4\pi\rho)$, four slices per cooperation length, which is where FLASH1's pulse energy converged, with the Pierce parameter taken from the beam the deck describes and the first FEL element. In the steady state it stays one, because the whole charge sits in the one slice and the spacing there states the current rather than a resolution. A beam or field read from a dump is the one place it is recovered, since no dump format records it, and a spacing that is not a whole number of wavelengths is refused there.

`slicing%window_length` when zero holds every particle of the bunch and the line's slippage, one wavelength per undulator period, ahead of the window head, so radiation that slips forward has somewhere to go before it leaves. A drawn bunch spans its last particle to its first. A Gaussian has no last particle, so it spans the length beyond which a slice would hold less than one electron of the charge, a rule that scales with the charge where a fixed number of sigmas does not. A flat bunch takes its grid extent, and the steady state is one slice. A stated window says its own headroom and gets none added, and one that clips the bunch warns. `slicing%n_slice` states the same window as a count instead, and setting both is refused.

`slicing%current` states a flat current directly, which is Genesis4's `&beam current`: the bunch's z structure is then not read, so with a window it is a flat time-dependent run and with none it is the steady state of one slice, and no separate steady-state switch exists. Unset, the current comes from `beam_init%bunch_charge` and the bunch's own z distribution. Setting both is refused.

(param-wavefront-harmonics)=
**`wavefront_init%harmonics`** is a gap-free increasing list whose first entry must be `1`: the fundamental anchors the optical phase, the reference advance and the slippage schedule, and the harmonics ride on it. Anything else is refused. Nothing is ever summed across harmonics. See [](fel-physics.md#sec-field-set).

(lattice-attributes)=
## Lattice attributes

An FEL segment is a wiggler or undulator carrying Bmad's `fel_method` attribute, and the method is what selects the physics:

| `fel_method` | Meaning |
|---|---|
| `off` | The default. Not an FEL segment, and this program leaves it to the seam |
| `averaged` | The wiggle-averaged model. The production workhorse |
| `unaveraged` | Direct integration through the undulator field, with no period averaging |

It is Bmad's own attribute rather than anything this program registers, so `show ele` prints it, a written lattice keeps it, the parser refuses a misspelling, and every other Bmad program loads the lattice. It is class-settable as any attribute is, `wiggler::*[FEL_METHOD] = unaveraged`, and the two settings mix freely in one line. `tracking_method` is a separate axis and names how one particle crosses the element, as it does everywhere in Bmad: an FEL segment leaves it at `bmad_standard` and is the plain periodic wiggler its field attributes describe outside the FEL walk.

Space charge is per element too, through Bmad's own `space_charge_method` attribute rather than anything this program registers:

| `space_charge_method` | Meaning |
|---|---|
| `off` | The default. No space charge in this element |
| `slice` | The slice-binned longitudinal solve, with the FEL slices as the bins |

`fft_3d` and `cathode_fft_3d` are refused on an FEL element, since their solvers want a three-dimensional grid this walk does not build. Bmad's master switch applies as it does everywhere else: `bmad_com%csr_and_space_charge_on` must also be true, and when elements ask for `slice` while it is false the run says so and tracks without space charge. `space_charge_com%n_bin` is ignored, because the slices are the bins. Inside the Bmad seam the same attribute drives Bmad's own machinery, so one lattice reads the same way in every Bmad program.

The unaveraged mode's two numbers are run switches, [](#param-global-unaveraged-steps-per-period) and [](#param-global-unaveraged-ramp-periods), rather than per-element attributes.

(attr-fel-tracking)=
**`unaveraged`** is a full Newton-Lorentz quiver with no period averaging and no coupling factor anywhere in its inputs, and the run writes an energy ledger beside its other outputs. **`averaged`** is the wiggle-averaged model, whose transverse maps are chosen by [](#param-global-transport-model). See [](fel-physics.md#sec-unaveraged).

The FEL parameters themselves come from the element: `aw` from `b_max` and `l_period`, the helicity from `field_calc`, and the step from `ds_step`. Each is asserted at setup. See [](fel-physics.md#sec-element).

## Outputs

Every output file, and how to read it, is [](reading-output.md). In brief: `<out_root>.stats.h5` is the production statistics file, `<out_root>.diag.txt` is the text comparison instrument written only under `write_diag`, and the end state is dumped as openPMD, the one format this program writes. Files in Genesis4 format are written beside them for field-by-field comparison, and [](genesis4.md) covers the exchange.
