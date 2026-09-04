# Input reference

The normative reference for every namelist parameter Lucifer honors: its default, its meaning, and what refuses it. A run is one namelist file of three groups, each setting Bmad or program structures directly, which is Tao's `&tao_params` pattern. Defaults live in the struct declarations, and a harness check holds this page to them in both directions, so a stated default cannot drift from the real one and a new parameter cannot go undocumented.

| Group | Carries |
|---|---|
| [`&fel_params`](#group-fel-params) | The run: the lattice, the `global%` switches, Bmad's `bmad_com` and `space_charge_com`, and the `chamber_wake%` and `space_charge%` collective descriptions |
| [`&fel_beam_init`](#group-fel-beam-init) | The beam: Bmad's `beam_init%` description, the `resample%` resampler, source and output files, and the check instruments |
| [`&fel_wavefront_init`](#group-fel-wavefront-init) | The radiation: the `wavefront_init%` starting condition and `field_file`. The field record is the time window, so the window lives here |

Per-element settings are [lattice attributes](#lattice-attributes) rather than namelist parameters. How to build and run is [](user-guide.md). What the outputs hold is [](reading-output.md). The measured levels are [](validation.md). The physics is the manual, [](fel-physics.md). Every resolved input is written into the statistics file's `params/` group, so a finished run states its own configuration with every default explicit.

A flat `&fel_track_params` group is refused by name, with each parameter mapped to the group that now carries it.

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
| `global%source_model` | `"deposit"` | The FEL source: `"deposit"` per particle, or `"coherent"` for the coherent retrieval |
| `global%track_start` | `""` | Element locator bounding the walk below. Blank is the whole line |
| `global%track_end` | `""` | Element locator bounding the walk above. Blank is the whole line |
| `global%ran_seed` | `12345` | The one random seed, governing generation, resampling and noise |
| `global%reference_run` | `F` | Permit a lattice with no FEL element, so Bmad tracks every element |
| `global%load_only` | `F` | Build the initial state, dump it, exit without tracking |
| `global%migrate` | `F` | Move particles between slices when their ponderomotive phase leaves the slice |
| `global%migrate_check` | `F` | Verify phase continuity at every migration and report the worst deviation |
| `global%fp32_check` | `"off"` | Step a single-precision twin beside the FP64 path and record its divergence: `"lockstep"` or `"freerun"` |
| `global%fp32_mutate` | `F` | The instrument's self-test: coarsen the FP32 residual so the recorded level must move |
| `global%device` | `"off"` | Run the averaged FEL step on a device backend: `"metal"` on a build that carries it |

(param-global-interlude-model)=
**`global%interlude_model`** selects how the field-free elements are handled. `"bmad"` is the seam: Bmad's own `track1_bunch`, an exact theta mapping, and `wavefront_drift` for the field. `"genesis"` uses the transcribed interlude step everywhere, which prices what the seam changes. The slippage schedule is identical in both models. The physics is in the manual's interlude and seam sections.

(param-global-transport-model)=
**`global%transport_model`** selects the transverse maps inside an averaged FEL element. `"bmad"`, the default and the production model, is Bmad's own periodic-wiggler kernel with its end-field treatment and chromaticity through `p0/p`. `"genesis"` is Genesis4's maps verbatim, with the focusing matrix built from `aw`, `kx`, `ky` and chromaticity through `gammaz`. The second is validation-internal: a transcription comparison needs transcription-level transport, so the comparison tiers set it and no production run does. The two are priced against each other in the manual's element section, and the unaveraged mode integrates the field rather than applying a map, so the switch does not reach it.

(param-global-source-model)=
**`global%source_model`** is `"deposit"`, the standard per-particle scatter, or `"coherent"`, the coherent-Gaussian retrieval in which the spatially incoherent part of the source is dropped, the slice bunch factor keeps the physical shot noise, and the transverse shape is a guarded Gaussian. A profile that is measurably not Gaussian is refused by name. See [](fel-physics.md#sec-coherent-source).

(param-global-track-start)=
**`global%track_start`** and **`global%track_end`** bound the walk, using Bmad's element-locator syntax. The schedule (slippage, autophasing, break geometry) is always built on the full lattice, so a windowed run composes exactly with the full one: the first span followed by the second, started from the first's dumps, reproduces the one-shot run. The measured level is in [](validation.md#val-the-programs-own-identities).

(param-global-fp32-check)=
**`global%fp32_check`** arms the FP32 lockstep instrument, a check rather than a production mode: an FP32 twin of the averaged FEL advance steps beside the FP64 path from a shared state, and the per-quantity divergence goes to `<out_root>.fp32.txt` with the worst case in the footer. `"lockstep"` rebuilds the FP32 state from FP64 every step, so a wrong formula shows as a jump. `"freerun"` carries the FP32 longitudinal state and the FP32 field across steps, a complete single-precision run beside the FP64 one, so the compounding rate and the end-to-end observable divergence are measured; it covers a single-slice window and refuses more, since the twin keeps no slippage rotation of its own. The twin covers the fundamental-only, collective-free averaged advance and refuses anything else by name (harmonics, two polarizations, the coherent source, wakes, space charge, the unaveraged mode). A run whose FP32 residual cannot resolve the per-step phase increment is refused by the runtime guard rather than reported as a clean number. The FP64 physics is untouched: the instrumented run's outputs are byte-identical to the uninstrumented run's. The measured levels are in [](validation.md#val-fp32-lockstep). `fp32_mutate` coarsens the FP32 residual by eight mantissa bits so the harness can prove the check fails when the path is wrong. The instrument needs the single-precision FFTW (`libfftw3f`), since the field twin's transform is FFTW's own `fftwf` interface; a build without that library refuses `fp32_check` by name rather than transforming in another precision, and the build reports which way it went at configure time.

(param-global-device)=
**`global%device`** selects a device backend for the averaged FEL step. `"metal"` runs it on an Apple Silicon GPU in single precision: the beam and the field upload at each FEL element's entry and stay resident through the element, every integration step is one command buffer (transverse maps, the longitudinal push in a fixed-point phase chart, the source deposit and the FFT field solve), and the state returns to the host at the stats comb's positions and the element end. Between those positions the host arrays are stale by design, so `comb_ds_save` is also the readback knob. Combined with `fp32_check` the device instead takes the lockstep twin's role: the FP64 run is untouched and the device's per-step divergence goes to the same `.fp32.txt` stream, which is how the backend is judged ([](validation.md#val-device)). Everything the kernels do not cover is refused by name at setup or first use (harmonics, two polarizations, the coherent source, wakes, space charge, spontaneous radiation, migration, the unaveraged mode, the escaped-field bank, a transverse grid that is not a power of two from 64 to 1024), and a build without the backend refuses the knob itself: an unsupported configuration stops the run and never quietly takes the CPU path. The backend is built where the toolchain can carry it, which is macOS with a Clang-family Objective-C++ compiler rather than macOS alone, since it is ARC-managed Objective-C++ against the Metal framework. Every other build takes the refusing stub and says so at configure time. Running `lucifer` with no arguments names what the build in hand carries, either the device it found or the reason there is none, so the question needs no run and no deck to answer.

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
**`global%dump_beam_at`** and **`global%dump_field_at`** name elements through Bmad's own locator, so `class::name` syntax works. An entry matching no element is refused by name. Dumps are openPMD, and the field dump is unrotated into time order first.

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

Both are off by default, matching Genesis4's `&sponrad`. Fluctuations with `global%migrate = T` are refused by name: the quiet start cancels per beamlet, and migration scrambles the grouping. See [](fel-physics.md#sec-eom).

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
**`chamber_wake%model`** names the implementation, and `"genesis"` is the only value accepted today: the transcribed solver, which convolves its kernels with the weighted slice currents and produces one energy loss per slice. Anything else is refused by name rather than treated as the transcribed solver by default. The field exists so that a second implementation can arrive as a value here rather than as a rework.

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
**`space_charge%model`** names the implementation, and `"genesis"` is the only value accepted today: the transcribed solver, which works per slice on a radial grid over azimuthal modes and longitudinal harmonics. Anything else is refused by name. Note that this family is Lucifer's own, and `space_charge_com` in the same group is Bmad's global structure, a separate thing.

The default is expected to stay `"genesis"` when Bmad's own slice solver becomes the second value, and the reason is the physics rather than the order they arrived in. The transcribed solver carries two terms that matter inside an undulator and that Bmad's slice model does not have: the space charge of the microbunching itself, solved per longitudinal harmonic of the ponderomotive phase, and the longitudinal Lorentz factor, which at aw = 0.85 rms is worth a factor of 1.7 in the field's own scaling. What Bmad's model adds that this one lacks, a transverse defocusing kick, falls as the inverse cube of gamma and is an injector term rather than an undulator-line one. The transcribed path is also the measured one, at the level [](validation.md) records for the space-charge tier. A second value earns the default by measurement, not by being newer.

If neither term is configured, `nz = 0` with `longrange = F`, an element asking for `slice` is refused by name: the solve would cost its full price and return an exact zero, and which of the two terms was meant is worth asking.

See [](fel-physics.md#sec-spacecharge).

(element-wakes)=
### Element wakes

Elements carrying Bmad `sr_wake` definitions, either pseudomodes or a tabular `z_long`, act across the whole time window through slice concatenation. The conventions, the step size, the mid-element wiggler kick and the refusals are in [](fel-physics.md#sec-seamwake). Wakes and space charge are both refused by name in the unaveraged mode, which does not wire them into its step.

(group-fel-beam-init)=
## &fel_beam_init

```
&fel_beam_init
  beam_file = "Aramis-initial.beam.h5"
  nbins = 8
/
```

There are three ways to obtain a beam: start from a dump, generate one from a description, or import and resample a distribution.

### Files and sources

| Parameter | Default | Meaning |
|---|---|---|
| `beam_file` | `""` | openPMD particle dump to start from. Blank generates a beam instead |
| `dist_file` | `""` | openPMD-beamphysics file to import and resample |
| `use_beam_init` | `F` | Import path: generate the bunch from the `beam_init%` block instead of a file |
| `write_genesis_dist` | `""` | Write the bunch as a Genesis4 `&importdistribution` input |
| `write_openpmd_file` | `""` | Write the bunch as openPMD-beamphysics |

(param-beam-beam-file)=
**`beam_file`** reads openPMD only. A file that is not openPMD is refused by name, and the message carries the conversion command. `tests/scripts/convert_genesis.py` converts in both directions, and [](genesis4.md) describes the exchange.

### The bunch description: `beam_init%`

Bmad's standard `beam_init_struct`, the same block both generation paths read. One bulk-bunch description, two methods: the quiet-start loader evaluates it analytically per slice, and the import resamples real particles from it. The Twiss is always the lattice's.

| Parameter | Meaning |
|---|---|
| `beam_init%n_particle` | Macroparticles per slice, a positive multiple of `beamlet_size`. The import path reads it as bunch particles |
| `beam_init%a_norm_emit`, `beam_init%b_norm_emit` | Normalized emittances [m rad] |
| `beam_init%sig_pz` | Fractional momentum spread dP/P0 |
| `beam_init%bunch_charge` | Charge [C]. The current is derived from it, never input |
| `beam_init%sig_z` | Bunch length [m], read with `distribution_type(3)` |
| `beam_init%distribution_type(3)` | `"RAN_GAUSS"` for a Gaussian current profile, `"GRID"` for Bmad's uniform one |

(param-beam-init-contract)=
Every other `beam_init` field that is set is refused by name. A standard structure that silently dropped fields would be worse than a custom one, so the honored set above is the contract, and `check_beam_init_contract` enforces it. Note that this contract covers the quiet-start generator: the import path honors everything Bmad honors, since `init_beam_distribution` generates the bunch.

(param-beam-sig-z)=
**`beam_init%sig_z`** with `"RAN_GAUSS"` gives a Gaussian current profile evaluated at the slice centers, with the bunch centered in the window. A zero length is the steady state, one slice holding the whole charge, and it is refused by name for a time-dependent window. With `"GRID"` the profile is flat over the z extent of `grid(3)`. `beam_init%a_emit` and `b_emit` are refused: normalized emittances only, which is Bmad's preferred form. `sig_e` is deprecated Bmad-wide and does not exist here.

### The quiet start

| Parameter | Default | Meaning |
|---|---|---|
| `beamlet_size` | `8` | Beamlet size of the quiet start. The load is quiet below this |
| `shot_noise` | `F` | Impose physical shot noise. Time-dependent windows only |

(param-beam-shot-noise)=
**`shot_noise`** imposes the weighted Fawley loading. The noise-level algebra and the effective-count refusal guard are in [](fel-physics.md#sec-noise). The loader warns where Genesis4 silently clamps beamlets holding fewer real electrons than macroparticles.

### The resampler: `resample%`

Named after Genesis4's `&importdistribution` where an equivalent exists. `window_sample`, `ran_seed` and the seed field are shared with the generator: one seed governs generation, resampling and noise.

| Parameter | Default | Meaning |
|---|---|---|
| `resample%slice_width` | `0.01` | Sampling window over bunch length |
| `resample%n_particle_per_slice` | `8192` | Macroparticles per slice after resampling |
| `resample%beamlet_size` | `4` | Beamlet size of the resample |
| `resample%n_slice` | `0` | Slice count. Zero derives it from the bunch length and the spacing |

Genesis4's `match` and `center` are not ported: a Bmad lattice carries its Twiss and `beam_init` generates matched bunches already. The method is in [](fel-physics.md#sec-import).

### Check instruments

Not physics input. The validation harness sets these.

| Parameter | Default | Meaning |
|---|---|---|
| `split_weights` | `F` | Replace each particle by coincident copies carrying 1/3 and 2/3 of its weight |
| `resample_split_weights` | `F` | The same split, applied before the import |
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
  wavefront_init%window_sample = 3
  wavefront_init%seed_power = 5e3
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_half_width = 2e-4
/
```

| Parameter | Default | Meaning |
|---|---|---|
| `field_file` | `""` | openPMD EXT_Wavefront field dump to start from |
| `wavefront_init%lambda0` | `0` | Radiation wavelength [m]. Required for generation |
| `wavefront_init%window_length` | `0` | Time window [m]. Zero derives it from the bunch |
| `wavefront_init%window_sample` | `1` | Slice spacing in wavelengths, an integer |
| `wavefront_init%grid_n_pts` | `255` | Transverse grid points per side |
| `wavefront_init%grid_half_width` | `0` | Transverse grid half width [m] |
| `wavefront_init%seed_power` | `0` | Gaussian seed power [W]. Zero is a dark start |
| `wavefront_init%seed_waist_size` | `0` | Seed intensity 1/e^2 radius [m] |
| `wavefront_init%seed_polarization` | `"x"` | `"x"` or `"y"` |
| `wavefront_init%harmonics` | `1` | The field set. The first entry must be the fundamental |

(param-wavefront-lambda0)=
**`wavefront_init%lambda0`** is required, and deliberately not defaulted from the lattice resonance, since the first undulator may be detuned. Starting from a beam dump it is required too: the file carries the slice partition but not the radiation wavelength it was sliced on, and `window_sample` is required for the same reason.

(param-wavefront-window-length)=
**`wavefront_init%window_length`** derives from the bunch when zero, four sigma either side of a Gaussian or the grid extent when flat, and one slice in the steady state. Override it for slippage headroom. A window that clips the bunch warns.

(param-wavefront-harmonics)=
**`wavefront_init%harmonics`** is a gap-free increasing list whose first entry must be `1`: the fundamental anchors the optical phase, the reference advance and the slippage schedule, and the harmonics ride on it. Anything else is refused by name. Nothing is ever summed across harmonics. See [](fel-physics.md#sec-field-set).

(lattice-attributes)=
## Lattice attributes

An FEL segment is a wiggler or undulator whose `tracking_method` is one of Bmad's two FEL methods, and the method is what selects the physics:

| `tracking_method` | Meaning |
|---|---|
| `fel_averaged` | The wiggle-averaged model. The production workhorse |
| `fel_unaveraged` | Direct integration through the undulator field, with no period averaging |

Both are Bmad's own named methods rather than anything this program registers, so `show ele` prints them, a written lattice keeps them, and the parser refuses a misspelling. They are class-settable as any attribute is, `wiggler::*[TRACKING_METHOD] = fel_unaveraged`, and they mix freely in one line. A wiggler tracked any other way is not an FEL segment: `tracking_method = custom` means some other program's tracking and this program leaves it to the seam.

Space charge is per element too, through Bmad's own `space_charge_method` attribute rather than anything this program registers:

| `space_charge_method` | Meaning |
|---|---|
| `off` | The default. No space charge in this element |
| `slice` | The slice-binned longitudinal solve, with the FEL slices as the bins |

`fft_3d` and `cathode_fft_3d` are refused by name on an FEL element, since their solvers want a three-dimensional grid this walk does not build. Bmad's master switch applies as it does everywhere else: `bmad_com%csr_and_space_charge_on` must also be true, and when elements ask for `slice` while it is false the run says so and tracks without space charge. `space_charge_com%n_bin` is ignored, because the slices are the bins. Inside the Bmad seam the same attribute drives Bmad's own machinery, so one lattice reads the same way in every Bmad program.

The unaveraged mode's two numbers are per-element attributes, registered by this program and usable on any wiggler or undulator:

| Attribute | Default | Meaning |
|---|---|---|
| `fel_steps_per_period` | unset, meaning 20 | Unaveraged substeps per undulator period |
| `fel_ramp_periods` | unset, meaning 2 | Length of the sin^2 entry and exit ramps, in periods |

(attr-fel-tracking)=
**`fel_unaveraged`** is a full Newton-Lorentz quiver with no period averaging and no coupling factor anywhere in its inputs, and the run writes an energy ledger beside its other outputs. **`fel_averaged`** is the wiggle-averaged model, whose transverse maps are chosen by [](#param-global-transport-model). See [](fel-physics.md#sec-unaveraged).

(attr-fel-steps-per-period)=
**`fel_steps_per_period`** defaults to 20 when unset. Below 10 is refused by name, the floor our own coupling-factor convergence supports: [](validation.md) tabulates it at 10, 20 and 30 steps per period.

(attr-fel-ramp-periods)=
**`fel_ramp_periods`** defaults to 2 when unset. A true hard edge is a test configuration, and it has an explicit sentinel of `-1`, so silence never means hard edge. A ramp pair longer than the segment is refused by name.

The FEL parameters themselves come from the element: `aw` from `b_max` and `l_period`, the helicity from `field_calc`, and the step from `ds_step`. Each is asserted by name at setup. See [](fel-physics.md#sec-element).

## Outputs

Every output file, and how to read it, is [](reading-output.md). In brief: `<out_root>.stats.h5` is the production statistics file, `<out_root>.diag.txt` is the text comparison instrument written only under `write_diag`, and the end state is dumped as openPMD, the one format this program writes. Files in Genesis4 format are written beside them for field-by-field comparison, and [](genesis4.md) covers the exchange.
