# User guide: building, describing a run, and running it

Everything needed to get from a checkout to a tracked FEL run. The parameter-by-parameter reference is [`input-reference.md`](input-reference.md). What the output files hold and how to read them is [`reading-output.md`](reading-output.md). The measured levels and how to reproduce them are [`validation.md`](validation.md). The physics is the manual, [`fel-physics.md`](fel-physics.md).

## Building

Lucifer builds as part of Bmad's normal distribution build. No step here is specific to
it: the binary lands in `production/bin/lucifer`, or `debug/bin/lucifer` for a debug
build, alongside every other Bmad program.

Running the validation harness needs one thing beyond the binary, a Python environment
with numpy, h5py and pytest, described by `lucifer/wavefront/tests/environment.yml`. The
exact commands the harness runs are in [`validation.md`](validation.md).

## The three tracking methods

Selected per element by the `fel_tracking` lattice attribute, and they mix freely in one line.

| `fel_tracking` | method | what it is for |
|---|---|---|
| unset or `0` | averaged | The default. The wiggle-averaged (KMR) model on Bmad's own `bmad_standard` kernel maps. The production workhorse. |
| `1` | unaveraged | Direct RK4 integration through the analytic undulator field: no averaging, no resonance approximation, fc and JJ nowhere in its inputs. A production method whose ~30x cost per step buys full quiver dynamics, energy accounting the beam actually pays, polarization-agnostic coupling, and arbitrary harmonic content in the current. It is also the tree's referee, since it shares no approximation with the averaged path. |
| `-1` | transcribed Genesis4 | Genesis4's own transverse maps, verbatim. Validation-internal: the Genesis4 comparison tiers select it through wrapper lattices, and no production lattice writes it. |

Lattices name these with one-line variables rather than raw numbers, so `fel_unaveraged = 1` then `wiggler::*[FEL_TRACKING] = fel_unaveraged`.

## Running

A run is one input file, which names the lattice:

```
lucifer lucifer.in
```

The name is an argument and any name works. `lucifer.in` is the convention the
examples follow, so that a directory reads the same way twice, and it is not a default
the program looks for. Running `lucifer` with no argument prints that much and the
version it was built from.

The file vocabulary of a run, in four extensions:

| Extension | Written by | What it is |
|---|---|---|
| `.in` | you | The input deck: `&fel_params`, `&fel_beam_init` and `&fel_wavefront_init`, naming one lattice |
| `.bmad` | you | A Bmad lattice, named by the deck's `lat_file`. Variants call a shared one and override it |
| `.h5` | the run | HDF5, openPMD throughout: the statistics file, and the beam and field dumps |
| `.txt` | the run | The named streams, one per subject: diag, ledger, import, migration, wake |

What each written file holds is the output section below, and what to do with it is
[`reading-output.md`](reading-output.md).

The examples are the fastest way in. Each is a directory of real input files, a README
with the numbers measured on it, and nothing else to fetch:

```
cd lucifer/examples/steady_state && lucifer lucifer.in
```

The one that needs Genesis4, and the only one with a runner, is the three-way
saturation comparison:

```
cd lucifer/examples/saturation_demo && ./run.sh
```

The validation harness is a separate thing, described with its commands in
[`validation.md`](validation.md).

## Files

| Path | Contents |
|---|---|
| `lucifer/doc/fel-physics.md` | **The physics manual**: equations, conventions, Genesis4 provenance and validation pointers, one section per subsystem |
| `lucifer/code/fel_beam_mod.f90` | Packed particle slices in Bmad coordinates plus per-particle weight, openPMD `.beam.h5` dump read/write through Bmad's own beam I/O, copy-only `coord_struct` conversion, weighted beam diagnostics with `N_eff` |
| `lucifer/code/fel_track_mod.f90` | The transcribed FEL step: transverse push with natural focusing, RK4 ponderomotive advance, source deposition, FFT field solve; the rotating-record slippage machinery (`fel_slip_struct`, `fel_apply_slippage`, `fel_field_index`); plus the transcribed Genesis4 interlude model |
| `lucifer/program/lucifer.f90` | The tracker: walks a Bmad lattice, FEL steps inside wiggler/undulator elements with `tracking_method = custom` (parameters from the lattice attributes, described in the FEL element section), seam everywhere else, slippage schedule transcribed from Genesis4; generates its own quiet-start beam and seed field when no dumps are named |
| `lucifer/examples/` | Self-contained single-command examples, one directory per feature, each with its own README and measured numbers. Its index says which example shows which feature |
| `lucifer/tests/scripts/check_shot_noise.py` | Statistical check: `<\|b(h)\|^2> = 1/N_lambda` over many seeds, uniform and nonuniform weights |
| `lucifer/tests/scripts/check_sase_startup.py` | Cross-code check: SASE startup power, our loader against Genesis4's, independent RNGs |
| `lucifer/tests/scripts/check_migration.py` | Migration checks: charge conservation under heavy migration, exact phase continuity, window residency, no-op bit identity |
| `lucifer/code/fel_import_mod.f90` | The distribution import: a bunch_struct resampled into FEL slices by Genesis4's importdistribution method, plus the Genesis4-distribution-file writer. Also see the distribution import in [`validation.md`](validation.md) |
| `lucifer/tests/scripts/check_seam_wake.py` | Seam-wake checks: closed-form pseudomode, exact causality with the d8 direction cross-check, z_long kernel cross-validation, split-weight, thread determinism |
| `lucifer/tests/scripts/check_import.py` | Import checks: exact current profile vs Genesis4 on the same file, match exactness, split-weight invariance, openPMD round trip, thread determinism; statistical Twiss recovery and startup power |
| `lucifer/code/fel_collective_mod.f90` | Wakes and space charge at Genesis4's granularity: the numerical resistive-wall impedance (Bane-Stupakov, a separable future Bmad port), geometric and roughness kernels, the causal convolution, the per-slice eloss application, and the short/long-range space-charge solvers behind a swappable interface |
| `lucifer/tests/genesis4/collective/` | Genesis4 decks: the collective tiers, importing the shared TD dumps |
| `lucifer/tests/scripts/check_collective.py` | Collective checks: exact wake energy bookkeeping, sigma_energy invariance, stale-wake structure under migration |
| `lucifer/code/fel_unaveraged_mod.f90` | The unaveraged mode: full Newton-Lorentz quiver through the analytic undulator field with sin² end ramps, radiation kick + coupling-free source, energy ledger; no fc/faw anywhere (grep-checked) |
| `lucifer/tests/scripts/check_unaveraged.py` | Unaveraged checks: energy ledger, ballistic/handoff, fc measured vs closed form (planar, helical, h=3), step-size convergence, priced gain-curve comparison |
| `lucifer/tests/bmad/unavg_probe_*.bmad` | The paired coupling probes (12/20 periods, planar and helical) |
| `lucifer/tests/run_fel_benchmark.sh` | The whole validation, one command |
| `lucifer/tests/run_perf_benchmark.sh` | Performance head-to-head vs Genesis4, serial and at the machine's performance-core count. Also see the parallelism record in [`validation.md`](validation.md) |
| `lucifer/tests/scripts/compare_fel.py` | Comparison: three steady-state tiers plus five time-dependent tiers against Genesis4, plus the split-weight invariance check |
| `lucifer/tests/scripts/plot_fel_compare.py` | Visual companion to `compare_fel.py`: overlays one tier's Bmad and Genesis4 curves (power, checked relative difference, per-slice exit power, bunching) from the diag file and the `.out.h5` |
| `lucifer/tests/genesis4/time_dependent/Aramis-td-sase.in` | Genesis4 deck: pure SASE, the TD window with the seed removed (`power = 0`), writing its own shared dumps |
| `lucifer/tests/genesis4/steady_state/Aramis-ss.in`, `genesis4/Aramis.lat` | Genesis4 deck: Benchmark1-SASE steady state, modified as documented in the deck header |
| `lucifer/tests/genesis4/steady_state/Aramis-1seg.in`, `genesis4/Aramis-1seg.lat` | Genesis4 deck: one undulator segment, importing the same dumps |
| `lucifer/tests/genesis4/time_dependent/Aramis-td.in`, `Aramis-td-1seg.in` | Genesis4 decks: the time-dependent pair, 32 slices with shot noise |
| `lucifer/tests/bmad/aramis.bmad`, `aramis_1seg.bmad` | The Bmad lattices: real wiggler elements, `b_max` encoding aw = 0.84853 exactly in Bmad's constants |
| `lucifer/tests/genesis4/sweep/Aramis-td-s12.in`, `run_delz_sweep.sh` | The coarse-step measurement: one shared dump at `sample = 12`, tracker runs at several `ds_step` values (wrapper lattices overriding the element attribute) |

## Architecture

(Physics: manual [](fel-physics.md#sec-core), [](fel-physics.md#sec-chart), [](fel-physics.md#sec-field), [](fel-physics.md#sec-slippage).)

Inside FEL elements (wigglers with `tracking_method = custom`, described in the FEL
element section), `fel_track_und_step` advances the coupled system in steps of the
element's `ds_step` (Bmad's standard step attribute, which is Genesis4's `delz` living on
the lattice like every other parameter, and the bookkeeper's `num_steps =
round(l/ds_step)` is exactly Genesis4's unroll), in Genesis4's exact order: transverse half step, RK4 advance of (theta, gamma) with
the field gathered once per step, transverse half step, then source deposition and the
`exp(K2 dz)` field solve. Bmad tracking is never used inside (the rule:
`symp_lie_bmad` resolves the wiggle motion the period-averaged map assumes away).

Everywhere else (the seam) the packed slice converts to `coord_struct`s by plain
copies, `track1_bunch` tracks them, and the field goes through `wavefront_drift`.

**Coordinates.** The packed arrays store Bmad's `(x, px/p0, y, py/p0, z, pz)` plus a
per-particle weight (macroparticle charge in Coulombs, mapping to `coord_struct%charge`).
The ponderomotive phase is derived, not stored: Genesis4's per-particle theta splits into
a common reference advance (one scalar per beam, `phi0`, advanced once per step) and
the particle-specific lag carried by Bmad's z:

```
theta_j = phi0 - ks*tau_j,    tau_j = -z_j/beta_j = c*(t_j - t_ref)
```

This is the reference-offset formulation identified as the
FP32-safe one, and it removes the wrap hazard outright: z does not wrap, so slice
migration has no theta-wrap-plus-index update to get wrong. The longitudinal RK4 still
runs in Genesis4's (theta, gamma) chart as per-step working variables. theta <-> z is an
affine map, under which RK4 is exactly invariant, so the verbatim transcription survives.
The per-step gamma <-> pz round trip costs ~1 ulp of gamma per step, which is what moved
tier1 from 2.5e-12 to 2.8e-11 when this representation replaced the Genesis4-coordinate
one. The particles live in the packed arrays at all times except inside `track1_bunch`.
`coord_struct` never enters the FEL step loop.

**Weights.** Every reduction is weighted: the source deposition scales per particle as
`c*w_j/slice_spacing` (Genesis4's `current/N` for uniform weights), bunching is
`|sum w e^(i theta)|/sum w`, and `N_eff = (sum w)^2/sum w^2` is a per-slice diagnostic.
The Genesis4 dump format carries no weights, so imports are uniform and a Genesis4
comparison can only exercise the uniform case. The split-weight check below covers the
rest.

**Units and constants.** The field is in V/m throughout (the `wavefront_struct`
convention) and all constants are Bmad's (`m_electron`, `mu_0_vac*c_light`). Expressed
that way the formulas are simpler than Genesis4's internal-unit originals (power is
`sum|E|^2*dA/(2*Z0)`, the coupling is `fc*conj(E)/(sqrt(2)*m_electron)`). Genesis4's
impedance constant is truncated, 376.73 against the exact 376.7303..., so the two codes
differ by 8.3e-7 relative before any physics runs. That difference is the floor of every
Genesis4 comparison, and the checks are sized to it rather than to round-off.

**Time dependence and slippage.** A multi-slice starting dump makes a time-dependent run,
a single-slice dump the steady state. There is no separate switch, the same rule as
Genesis4, whose imports carry the time window. Beam slice `is` couples to field slice
`1 + mod(is-1+first, nslice)`: the field record is a circular buffer over the wavefront's
slice index, rotated by slippage rather than moved, as Genesis4 does. Slippage accumulates at
`dz*(1+aw^2)/(2*gamma0^2*lambda)` wavelengths per undulator step and rotates the record
one slice whenever the accumulation exceeds `0.8*sample`, reduced to one shared-memory
node, where the MPI ring exchange is the identity. The
slice rotating out at the head of the time window is discarded and re-enters zeroed at
the tail. Drifts autophase `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths onto the last
interlude element before each undulator, and the end-of-lattice
fixup is **unguarded in Genesis4** (`+1` even with no trailing drift at all), which
costs one final rotation if transcribed with a guard Genesis4 does not have (found at
0.84 of the final field). Everything reading the record in time order
(the per-slice diagnostics, the final field dump) unrotates through `fel_field_index`,
as Genesis4 does on the way out. Beam slices never rotate.

## The FEL element: parameters live on the lattice

(Physics: manual [](fel-physics.md#sec-element).)

An FEL segment is a real Bmad `wiggler` (or `undulator`) element carrying
`tracking_method = custom`, Bmad's own semantics for "the program supplies the
tracking", which this driver does. Recognition is by key and tracking method, never by
name. There are no per-undulator namelist parameters. The FEL parameters derive from
the element attributes:

- `aw` (rms, Genesis4's convention) from the peak field: `K = c*b_max/(k_u*m_e c^2)`,
  exactly and independent of the reference energy. Helical `aw = K`, planar
  `aw = K/sqrt(2)`. The benchmark lattices write `b_max` as an expression in Bmad's own
  constants so it encodes `aw = 0.84853` exactly, and the round trip still differs by
  1 ulp, measured below.
- Helicity from `field_calc` (`helical_model` / `planar_model`), never from the stored
  `k1x`/`k1y` wiggler attributes, whose helical sign convention disagrees with Bmad's
  own tracking locals.
- Genesis4's natural-focusing split from the helicity defaults (`kx = ky = 0.5*k_u^2`
  helical, `0`/`k_u^2` planar). Bmad's `kx` roll-off
  attribute is not yet mapped onto that split and must be zero.

Outside the driver's own FEL walk the element is just a periodic wiggler, tracked by
Bmad's standard kernel through two hooks wired before `bmad_parser`:
`track1_custom_ptr` and `make_mat6_custom_ptr` both delegate to `track_a_wiggler`, so
the reference time acquires the resonant undulation delay from Bmad's own code
and the transfer-matrix bookkeeping works. Both hooks are load-bearing:
`mat6_calc_method` resolves to custom alongside the tracking method, and Bmad calls
through a null `make_mat6_custom_ptr` (a jump to address zero) if a program sets only
the tracking hook.

The 7.5 assertions are enforced at the element's first touch (the reference pass
inside `bmad_parser`, through those hooks) and refuse by name: a missing `b_max`
(Bmad's own kernel would silently give `osc_amplitude = 0`: no field, no resonance, no
error), a missing `l_period` (same silence), and a fieldmap `field_calc` (which
segfaults the parse-time reference tracking if allowed through). Enforcing them any
later is provably too late. The missing-`b_max` lattice otherwise dies downstream on
an unrelated generation message. The assertions live in one routine
(`fel_assert_wiggler_sane`). A second inline copy was tried and removed because
redundant assertions mask the removal of either copy under mutation testing. All three
refusals are permanent harness checks.

**Anchored against the namelist-driven reference:** the element-driven full TD line
reproduces the last namelist-driven build's run to a max relative difference of
3.4e-12, zero at z = 0 and growing with gain: exactly the amplification of the 1-ulp
`aw` difference from deriving `aw` through the `b_max` attribute round trip rather than
reading it from input.

**The FEL tracking mode is a per-element lattice attribute** (`fel_tracking`,
registered by the driver, class-settable as `wiggler::*[fel_tracking] = ...`), so
averaged and unaveraged segments mix freely in one line. Unset/0, the default, is
averaged with the transverse maps of Bmad's own `bmad_standard` periodic-wiggler
kernel, flattened per `ds_step` (`track_a_wiggler`'s matrix with the octupole-like
end kicks, chromatic via `p0/p`). `1` is the unaveraged mode, a production method whose cost buys the full quiver dynamics (see the manual's unaveraged-mode section). `-1` is
averaged with the transcribed-Genesis4 focusing (matrix from `aw`, `kx`, `ky`,
chromatic via `gammaz`), and it is validation-internal: the Genesis4 tiers require
transcription-level transport and select it via the `*_val.bmad` wrapper lattices. No
production lattice writes it. The two averaged models are priced, measured over the
full time-dependent line (32 slices, 90 records): power differs by 5.0e-5 max (exit
total power 3.0e-7), on-axis intensity 7.3e-4, spot sizes 6.7e-6, wrapped bunching
phase 1.3e-2 rad max, for +3.8% runtime. That is period-averaged Genesis4 focusing
against Bmad's end-field treatment. The unaveraged parameters are attributes too:
`fel_steps_per_period` (unset becomes 20, below 10 refused) and `fel_ramp_periods`
(unset becomes 2). A true hard edge, the test configuration, is the explicit sentinel
-1, because an attribute's unset value is 0 and a silent hard edge would reintroduce
the K/gamma handoff hazard.

The examples directory exercises the heterogeneity this buys: `examples/taper/` is the same
line with the last two cells' undulators a second element definition with `b_max` 0.4%
lower: bit-identical to the untapered run until the taper starts, 12.7x its exit
power after (see `examples/README.md`).

## Program structure

The tracker is laid out the way Tao is laid out: the input structs in `fel_struct`
(defaults in the declarations), the namelist layer quarantined in `fel_input_mod`
(three groups in one file -- `&fel_params` with `global%...`, `bmad_com`,
`space_charge_com`, `chamber_wake%`, `space_charge%`. `&fel_beam_init` with `beam_init%...`, `resample%...`
and files. `&fel_wavefront_init` with `wavefront_init%...` and `field_file`), and the
work in library modules (`fel_setup_mod` / `fel_init_mod` / `fel_track_line_mod` /
`fel_io_mod`) over one explicit `fel_run_struct`. Nothing in the library stops --
errors return and the driver decides, and `track_fel_line` is re-entrant (twice in
one process is bit-identical to two processes. `tests/lucifer_smoke_test.f90` drives the
library with no namelist anywhere). A flat `&fel_track_params` group is refused by
name, with each parameter's current home in the message. `global%comb_ds_save` is the stats
comb (Bmad's `ds_save` semantics with one deliberate difference: an element end is
always a record, whatever the comb, which is what lets the file carry one record axis
and a mask instead of a second axis. Default 0 = every record).
`global%track_start`/`track_end` bound the walk over a schedule always
built on the full lattice (windowed runs compose exactly). `stats.h5` carries `meta/`
provenance (the resolved input echo, which lattice, a timestamp and the Bmad version) as
datasets, and `global%record_environment` adds the user name and working directory for a
lab notebook, off by default because a stats file travels.

## Output: stdout is for humans, files are for programs

The program's terminal output is formatted for reading and carries no contract beyond
that: units chosen for the eye, columns aligned, wording free to improve. **Parsing it
is discouraged.** Everything a program should read is written to a file, at full
precision:

| File | What |
|---|---|
| `<out_root>.stats.h5` | The run: per-record beam, field and Twiss data, the lattice table, and the axes and parameters to read them by. Self-describing (see [`reading-output.md`](reading-output.md)) |
| `<out_root>-final.beam.h5`, `-final.wf.h5` | Final beam and field, openPMD, the only dump format this code writes. The beam file carries the per-particle weight, which no Genesis4 dump can. `lucifer/tests/scripts/convert_genesis.py` converts either kind to Genesis4 conventions, for feeding Genesis4 |
| `<out_root>.diag.txt` | The per-record Genesis4-comparison instrument, one row per slice per record |
| `<out_root>.ledger.txt` | The unaveraged energy ledger, one row per record step |
| `<out_root>.import.txt` | The distribution import's analysis moments and per-slice current profile. Written at import time, so it exists under `load_only = T` |
| `<out_root>.migration.txt` | Slice migration: one row per event (s, particles moved, charge dropped, phase-continuity residual) plus the run summary |

This is Tao's division of labour. `show ele` is for people, `show value` and pytao are
for precision. It is why the validation suite reads files rather than scraping the
screen.

What a run looks like: a framed block naming the configuration, a progress table, and a
completion block listing what was written.

```
================================================================================
 Lucifer -- FEL tracking in Bmad, Bmad version 20260810-0
--------------------------------------------------------------------------------
 Lattice     aramis.bmad
             49 elements, 57.000 m, 12 FEL segments
 Beam        1 slice x 8192 particles, gamma0 = 11357.82
 Radiation   lambda0 = 100.000 pm, slice spacing 100.000 pm, 1 field(s), grid half width 200.000 um
 Switches    sr wakes F, space charge F, radiation damping F
 Output      out_root "steady_state", threads 12
================================================================================
     %        s       ele       step        power       energy    <|b|>  elapsed  element
    7.0     3.990      1/49      89/89     4.230 kW     1.411 fJ   0.0006    0:00  UND
   15.3     8.740      5/49      89/89    34.562 kW    11.529 fJ   0.0017    0:00  UND
  100.0    57.000     48/49              761.499 MW   254.009 pJ   0.1854    0:06  D2
--------------------------------------------------------------------------------
 Done        57.000 m, 48 element ends
 Exit        power 761.499 MW, pulse energy 254.009 pJ, <|b|> 0.1854
 Wrote       steady_state-final.beam.h5                      710.760 kB
             steady_state-final.wf.h5                        1.048 MB
             steady_state.stats.h5                           841.220 kB
================================================================================
```

Four choices in that table are deliberate. Values carry an **SI prefix** chosen per
value, so one column holds both `4.230 kW` at startup and `761.499 MW` at saturation
without losing the early digits a fixed unit would flatten. The **element name is
last**, unpadded, so every numeric column is fixed however long a name gets: no
truncation rule, no name shoving the numbers rightward. The **step column is blank for
one-step elements**, which is most of a real lattice. And the distance column is **`s`**,
Bmad's arc length along the reference orbit: through a bending break the arc exceeds the
chord the light takes, which is exactly what the light-path correction accounts for.

The physics columns are filled whatever the comb setting: with per-record rows they come
from the current record, and with `comb_ds_save < 0` from the element-end row, which is
filled just before the row that reads it. Between element ends in that mode the last
element-end values are carried forward rather than blanked. A mid-element row exists to
answer "is it alive and roughly where".

There is exactly **one** deliberate exception. A refusal has to be recognizable *by
name*, so the checks match the ALL-CAPS text of a refusal together with a nonzero exit
status. Those texts are a supported contract and change only together with their
checks, in the same commit. Nothing else printed to the terminal is.

