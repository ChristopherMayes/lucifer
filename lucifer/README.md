# FEL tracker, validated against Genesis 1.3 Version 4

An FEL tracker whose physics is transcribed from Genesis 1.3 Version 4 (GPL permits
transcription), embedded in Bmad's lattice machinery by the seam of the design brief's
section 4.1, and validated against Genesis over its `benchmark/Benchmark1-SASE`
configuration from bitwise-identical starting states.

**The physics reference is the manual, [`doc/fel-physics.tex`](doc/fel-physics.tex)**
(`make` in `doc/` builds the PDF): every equation the code integrates, each subsystem's
Genesis provenance, and which check pins it at what measured level. This README carries
the measured numbers, the validation methodology, and how to run everything; where it
touches physics it cites the manual's sections (`sec:core`, `sec:slippage`, ...).

No harmonics beyond the coupling formula. Everything else the brief's section 10 sequences through step 10 is in — including the FEL element proper (deliverable 9): undulator segments are real Bmad wiggler elements with `tracking_method = custom`, their FEL parameters derived from lattice attributes, with the brief's 7.5 assertions enforced by name; see the FEL element section. Per-particle weights are carried from day one (brief section 5): the
packed arrays store one, every reduction uses it, and the split-weight check below tests
the nonuniform case that no Genesis comparison can reach. The FEL step is OpenMP-parallel
over slices (deliverable 5), with thread-count independence -- bit-identical results at
any thread count -- as a checked property; see the parallelism section. The built-in
loader imposes physical shot noise generalized to weights (deliverable 6), checked
statistically and cross-checked against Genesis's independent loader; see the shot-noise
section.

## Files

| Path | Contents |
|---|---|
| `bsim/fel/doc/fel-physics.tex` | **The physics manual**: equations, conventions, Genesis provenance and validation pointers, one section per subsystem |
| `bsim/modules/fel_beam_mod.f90` | Packed particle slices in Bmad coordinates plus per-particle weight, Genesis `.par.h5` dump read/write (converting), copy-only `coord_struct` conversion, weighted beam diagnostics with `N_eff` |
| `bsim/modules/fel_track_mod.f90` | The transcribed FEL step: transverse push with natural focusing, RK4 ponderomotive advance, source deposition, FFT field solve; the rotating-record slippage machinery (`fel_slip_struct`, `fel_apply_slippage`, `fel_field_index`); plus the transcribed Genesis interlude model |
| `bsim/fel/fel_track_test.f90` | The tracker: walks a Bmad lattice, FEL steps inside wiggler/undulator elements with `tracking_method = custom` (parameters from the lattice attributes; see the FEL element section), seam everywhere else, slippage schedule transcribed from `Lattice::calcSlippage`; generates its own quiet-start beam and seed field when no dumps are named |
| `bsim/fel/examples/` | Self-contained single-command examples (no Genesis, no dump files): a seeded steady-state run and a pure-SASE time-dependent run of the benchmark line, with plotting script and README |
| `bsim/fel/tests/scripts/check_shot_noise.py` | Statistical check: `<\|b(h)\|^2> = 1/N_lambda` over many seeds, uniform and nonuniform weights |
| `bsim/fel/tests/scripts/check_sase_startup.py` | Cross-code check: SASE startup power, our loader against Genesis's, independent RNGs |
| `bsim/fel/tests/scripts/check_migration.py` | Migration checks: charge conservation under heavy migration, exact phase continuity, window residency, no-op bit identity |
| `bsim/modules/fel_import_mod.f90` | The distribution import (brief 10 step 10): a bunch_struct resampled into FEL slices by Genesis's importdistribution method, transcribed from SDDSBeam.cpp, plus the Genesis-distribution-file writer; see the distribution-import section |
| `bsim/fel/tests/scripts/check_seam_wake.py` | Seam-wake checks: closed-form pseudomode, exact causality with the d8 direction cross-check, z_long kernel cross-validation, split-weight, thread determinism |
| `bsim/fel/tests/scripts/check_import.py` | Import checks: exact current profile vs Genesis on the same file, match exactness, split-weight invariance, openPMD round trip, thread determinism; statistical Twiss recovery and startup power |
| `bsim/modules/fel_collective_mod.f90` | Wakes and space charge at Genesis's granularity: the numerical resistive-wall impedance (Bane-Stupakov, a separable future Bmad port), geometric and roughness kernels, the causal convolution, the per-slice eloss application, and the short/long-range space-charge solvers behind a swappable interface |
| `bsim/fel/tests/genesis4/collective/` | Genesis decks: the collective tiers, importing the shared TD dumps |
| `bsim/fel/tests/scripts/check_collective.py` | Collective checks: exact wake energy bookkeeping, sigma_energy invariance, stale-wake structure under migration |
| `bsim/modules/fel_unaveraged_mod.f90` | The unaveraged verification mode (brief 6.6, manual `sec:unaveraged`): full Newton-Lorentz quiver through the analytic undulator field with sin² end ramps, radiation kick + coupling-free source, energy ledger; no fc/faw anywhere (grep-checked) |
| `bsim/fel/tests/scripts/check_unaveraged.py` | Unaveraged checks: energy ledger, ballistic/handoff, fc measured vs closed form (planar, helical, h=3), step-size convergence, priced gain-curve comparison |
| `bsim/fel/tests/bmad/unavg_probe_*.bmad` | The paired coupling probes (12/20 periods, planar and helical) |
| `bsim/fel/tests/run_fel_benchmark.sh` | The whole validation, one command |
| `bsim/fel/tests/run_perf_benchmark.sh` | Performance head-to-head vs Genesis4, serial and at the machine's performance-core count; see the parallelism section |
| `bsim/fel/tests/scripts/compare_fel.py` | Comparison: three steady-state tiers plus five time-dependent tiers against Genesis, plus the split-weight invariance check |
| `bsim/fel/tests/scripts/plot_fel_compare.py` | Visual companion to `compare_fel.py`: overlays one tier's Bmad and Genesis curves (power, checked relative difference, per-slice exit power, bunching) from the diag file and the `.out.h5` |
| `bsim/fel/tests/genesis4/time_dependent/Aramis-td-sase.in` | Genesis deck: pure SASE — the TD window with the seed removed (`power = 0`), writing its own shared dumps |
| `bsim/fel/tests/genesis4/steady_state/Aramis-ss.in`, `genesis4/Aramis.lat` | Genesis deck: Benchmark1-SASE steady state, modified as documented in the deck header |
| `bsim/fel/tests/genesis4/steady_state/Aramis-1seg.in`, `genesis4/Aramis-1seg.lat` | Genesis deck: one undulator segment, importing the same dumps |
| `bsim/fel/tests/genesis4/time_dependent/Aramis-td.in`, `Aramis-td-1seg.in` | Genesis decks: the time-dependent pair, 32 slices with shot noise |
| `bsim/fel/tests/bmad/aramis.bmad`, `aramis_1seg.bmad` | The Bmad lattices: real wiggler elements, `b_max` encoding aw = 0.84853 exactly in Bmad's constants |
| `bsim/fel/tests/genesis4/sweep/Aramis-td-s12.in`, `run_delz_sweep.sh` | The coarse-step measurement (brief 8.3): one shared dump at `sample = 12`, tracker runs at several `ds_step` values (wrapper lattices overriding the element attribute) |

## Running

```
cd <bmad-ecosystem>
BUILD_PRODUCTION=N ./util/conda_compile                      # builds fel_track_test
./bsim/fel/tests/run_fel_benchmark.sh [--genesis <path to genesis4>]
```

The harness runs Genesis six times (full line and single segment, steady state and time
dependent, plus the collective tiers), the Bmad tracker for every tier and check, and
prints the largest relative difference of each check. It fails loudly if the genesis4 binary is missing; there is no comparison
without it, so there is nothing to skip to. Genesis must be built with FFTW, since the
benchmark runs with `fft_fieldsolver=true` (the Bmad tracker transcribes the FFT solver;
Genesis's default ADI solver is out of scope).

## Architecture

(Physics: manual `sec:core`, `sec:chart`, `sec:field`, `sec:slippage`.)

Inside FEL elements (wigglers with `tracking_method = custom`; see the FEL element
section), `fel_track_und_step` advances the coupled system in steps of the element's
`ds_step` (Bmad's standard step attribute -- Genesis's `delz`, living on the lattice
like every other parameter; the bookkeeper's `num_steps = round(l/ds_step)` is exactly
Genesis's unroll), in Genesis's exact order: transverse half step, RK4 advance of (theta, gamma) with
the field gathered once per step, transverse half step, then source deposition and the
`exp(K2 dz)` field solve. Bmad tracking is never used inside (the brief's rule:
`symp_lie_bmad` resolves the wiggle motion the period-averaged map assumes away).

Everywhere else -- the seam -- the packed slice converts to `coord_struct`s by plain
copies, `track1_bunch` tracks them, and the field goes through `wavefront_drift`.

**Coordinates.** The packed arrays store Bmad's `(x, px/p0, y, py/p0, z, pz)` plus a
per-particle weight (macroparticle charge in Coulombs, mapping to `coord_struct%charge`).
The ponderomotive phase is derived, not stored: Genesis's per-particle theta splits into
a common reference advance -- one scalar per beam, `phi0`, advanced once per step -- and
the particle-specific lag carried by Bmad's z:

```
theta_j = phi0 - ks*tau_j,    tau_j = -z_j/beta_j = c*(t_j - t_ref)
```

This is the reference-offset formulation the design brief's section 8 identifies as the
FP32-safe one, and it removes the brief's 6.4 hazard outright: z does not wrap, so slice
migration has no theta-wrap-plus-index update to get wrong. The longitudinal RK4 still
runs in Genesis's (theta, gamma) chart as per-step working variables -- theta <-> z is an
affine map, under which RK4 is exactly invariant, so the verbatim transcription survives;
the per-step gamma <-> pz round trip costs ~1 ulp of gamma per step, which is what moved
tier1 from 2.5e-12 to 2.8e-11 when this representation replaced the Genesis-coordinate
one. The particles live in the packed arrays at all times except inside `track1_bunch`;
`coord_struct` never enters the FEL step loop.

**Weights.** Every reduction is weighted: the source deposition scales per particle as
`c*w_j/slice_spacing` (Genesis's `current/N` for uniform weights), bunching is
`|sum w e^(i theta)|/sum w`, and `N_eff = (sum w)^2/sum w^2` is a per-slice diagnostic.
The Genesis dump format carries no weights, so imports are uniform and a Genesis
comparison can only exercise the uniform case; the split-weight check below covers the
rest.

**Units and constants.** The field is in V/m throughout — the `wavefront_struct`
convention — and all constants are Bmad's (`m_electron`, `mu_0_vac*c_light`); expressed
that way the formulas are simpler than Genesis's internal-unit originals (power is
`sum|E|^2*dA/(2*Z0)`, the coupling is `fc*conj(E)/(sqrt(2)*m_electron)`). During
development the tracker instead transcribed Genesis's internal units and constants —
including its truncated impedance, 376.73 against the exact 376.7303... — and agreed at
transcription level: tier1 2.8e-11, tier2_genesis 5.9e-8, on commit 236dc372f. With that
validation banked, the code moved to Bmad constants by decision; the 8.3e-7 relative
impedance difference is now the floor of every Genesis comparison, and the checks below
are sized to it.

**Time dependence and slippage.** A multi-slice starting dump makes a time-dependent run,
a single-slice dump the steady state -- no separate switch, the same rule as Genesis,
whose imports carry the time window. Beam slice `is` couples to field slice
`1 + mod(is-1+first, nslice)`: the field record is a circular buffer over the wavefront's
slice index, rotated by slippage rather than moved (Genesis's `Field::first`;
`BeamSolver.cpp:66`, `FieldSolverFFT.cpp:54`). Slippage accumulates at
`dz*(1+aw^2)/(2*gamma0^2*lambda)` wavelengths per undulator step and rotates the record
one slice whenever the accumulation exceeds `0.8*sample` (`Control::applySlippage`,
reduced to one shared-memory node -- the MPI ring exchange is the identity); the slice
rotating out at the head of the time window is discarded and re-enters zeroed at the
tail. Drifts autophase `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths onto the last
interlude element before each undulator (`Lattice::calcSlippage`), and the end-of-lattice
fixup is **unguarded in Genesis** -- `+1` even with no trailing drift at all -- which
costs one final rotation if transcribed with a guard Genesis does not have (FINDINGS.md
7.1; found at 0.84 of the final field). Everything reading the record in time order --
the per-slice diagnostics, the final field dump -- unrotates through `fel_field_index`,
as Genesis does at `writeFieldHDF5.cpp:86` and `Diagnostic.cpp:852`. Beam slices never
rotate.

## The FEL element (brief 7.5): parameters live on the lattice

(Physics: manual `sec:element`.)

An FEL segment is a real Bmad `wiggler` (or `undulator`) element carrying
`tracking_method = custom` — Bmad's own semantics for "the program supplies the
tracking", which this driver does. Recognition is by key and tracking method, never by
name; there are no per-undulator namelist parameters. The FEL parameters derive from
the element attributes:

- `aw` (rms, Genesis's convention) from the peak field: `K = c*b_max/(k_u*m_e c^2)`,
  exactly and independent of the reference energy; helical `aw = K`, planar
  `aw = K/sqrt(2)`. The benchmark lattices write `b_max` as an expression in Bmad's own
  constants so it encodes `aw = 0.84853` exactly — and the round trip still differs by
  1 ulp, measured below.
- Helicity from `field_calc` (`helical_model` / `planar_model`), never from the stored
  `k1x`/`k1y` wiggler attributes, whose helical sign convention disagrees with Bmad's
  own tracking locals (brief 7.5; nothing here cross-uses them).
- Genesis's natural-focusing split from the helicity defaults (`kx = ky = 0.5*k_u^2`
  helical, `0`/`k_u^2` planar, LatticeParser.cpp:328-333). Bmad's `kx` roll-off
  attribute is not yet mapped onto that split and must be zero.

Outside the driver's own FEL walk the element is just a periodic wiggler, tracked by
Bmad's standard kernel through two hooks wired before `bmad_parser`:
`track1_custom_ptr` and `make_mat6_custom_ptr` both delegate to `track_a_wiggler`, so
the reference time acquires the resonant undulation delay from Bmad's own code
(brief 7.5) and the transfer-matrix bookkeeping works. Both hooks are load-bearing:
`mat6_calc_method` resolves to custom alongside the tracking method, and Bmad calls
through a null `make_mat6_custom_ptr` (a jump to address zero) if a program sets only
the tracking hook.

The 7.5 assertions are enforced at the element's first touch — the reference pass
inside `bmad_parser`, through those hooks — and refuse BY NAME: a missing `b_max`
(Bmad's own kernel would silently give `osc_amplitude = 0`: no field, no resonance, no
error), a missing `l_period` (same silence), and a fieldmap `field_calc` (which
segfaults the parse-time reference tracking if allowed through). Enforcing them any
later is provably too late — the missing-`b_max` lattice otherwise dies downstream on
an unrelated generation message. The assertions live in ONE routine
(`fel_assert_wiggler_sane`); a second inline copy was tried and removed because
redundant assertions mask the removal of either copy under mutation testing. All three
refusals are permanent harness checks.

**Anchored against the namelist-driven reference:** the element-driven full TD line
reproduces the last namelist-driven build's run to a max relative difference of
3.4e-12, zero at z = 0 and growing with gain — exactly the amplification of the 1-ulp
`aw` difference from deriving `aw` through the `b_max` attribute round trip rather than
reading it from input.

**The FEL tracking mode is a per-element lattice attribute** (`fel_tracking`,
registered by the driver; class-settable as `wiggler::*[fel_tracking] = ...`), so
averaged and unaveraged segments mix freely in one line. Unset/0 — THE DEFAULT — is
averaged with the transverse maps of Bmad's own `bmad_standard` periodic-wiggler
kernel, flattened per `ds_step` (`track_a_wiggler`'s matrix with the octupole-like
end kicks, chromatic via `p0/p`). `1` is the unaveraged verification mode. `-1` is
averaged with the transcribed-Genesis focusing (matrix from `aw`, `kx`, `ky`,
chromatic via `gammaz`) — VALIDATION-INTERNAL: the Genesis tiers require
transcription-level transport and select it via the `*_val.bmad` wrapper lattices; no
production lattice writes it. The two averaged models are priced, measured over the
full time-dependent line (32 slices, 90 records): power differs by 5.0e-5 max (exit
total power 3.0e-7), on-axis intensity 7.3e-4, spot sizes 6.7e-6, wrapped bunching
phase 1.3e-2 rad max, for +3.8% runtime — period-averaged Genesis focusing vs Bmad's
end-field treatment. The unaveraged parameters are attributes too:
`fel_steps_per_period` (unset → 20; below 10 refused) and `fel_ramp_periods`
(unset → 2; a true hard edge — the test configuration — is the explicit sentinel -1,
because an attribute's unset value is 0 and a silent hard edge would reintroduce the
K/gamma handoff hazard).

The examples directory exercises the heterogeneity this buys: `examples/taper/` is the same
line with the last two cells' undulators a second element definition with `b_max` 0.4%
lower — bit-identical to the untapered run until the taper starts, 12.7x its exit
power after (see `examples/README.md`).

## Validation, from one command

Both codes start from the same Genesis `&write` dumps, so the initial state is bitwise
identical and no loader is reproduced. Genesis records diagnostics once at the start and
once per integration step, with each interlude element being a single step; the tracker
records at the same z positions (per slice, in time-window order, for the time-dependent
tiers). Measured, on the numbers this tree was developed against:

| Tier | What runs | Largest relative difference |
|---|---|---|
| `tier1` | One undulator segment: the FEL core alone | **1.8e-6** (the impedance-constant floor; was 2.8e-11 with Genesis's constants transcribed) |
| `tier1_unavg` | The same segment and dumps, tracked by the UNAVERAGED mode against Genesis's averaged run | **6.9e-2** — a priced model difference (sin² ramps vs hard edges, no averaging, integrator structure), dominated by the final-field phase; the power curve agrees at 6.1e-3 and per-particle gamma at 3.7e-6, and theta shows a CONSTANT ~6.6 rad ramp-phase offset with only 2.8e-3 rms about it |
| `tier2_genesis` | Full 6-FODO line, interludes via the transcribed Genesis model | **1.8e-5** (constants floor through full gain; was 5.9e-8) |
| `tier2_bmad` | Full line, interludes via the Bmad seam | **5.0e-2** (power curve 1.3e-2) -- a measured model difference, see below |
| `weight_split` | tier1 rerun with every particle split into coincident w/3 + 2w/3 copies, against the unsplit run | **3.6e-13** (Fortran vs Fortran; constants-independent) |
| `td1` | One undulator segment, 32 slices: FEL core plus slippage (accumulation, threshold, rotation, zero fill, end-of-lattice autophasing) | **8.5e-7** (constants floor) |
| `td2_genesis` | Full line time dependent, transcribed Genesis interludes: adds the drift autophasing schedule | **2.4e-6** (constants floor) |
| `td2_bmad` | Full line time dependent through the Bmad seam | **4.1e-2** -- the tier2_bmad transport model difference with slippage interleaved |
| `tdsase` | Full line, pure SASE: dark start (`power = 0`), shot noise on, both codes tracking the identical noisy beam and identical zero field from shared dumps -- the deterministic startup-from-noise comparison the seeded tiers and the statistical startup check leave uncovered | **2.3e-6** (constants floor; exit total power agrees at 1.9e-6) |
| `tdsc` | One segment TD, space charge on (short-range harmonics nz=2/nphi=1 plus long range) | **2.4e-4** (the epsilon_0-truncation floor of Genesis's longRange, 8.85e-12) |
| `tdwk` | One segment TD, all three wake kernels on (numerical-impedance resistive wall, gap, roughness) | **8.7e-7** (the impedance floor) |

Particle ordering is preserved by both codes (no sorting happens without one4one), so the
final dumps compare particle by particle, not just statistically, in every tier.

All outputs land in the benchmark's work directory (`--work-dir <path>`; without it
a temporary directory is used, removed on success and kept on failure). To *see* any
tier rather than check it, run from a kept work directory:
`scripts/plot_fel_compare.py <tier>.diag.txt <GenesisRoot>.out.h5` overlays the two codes' power
and bunching curves, plots the checked elementwise relative power difference along the
line, and the per-slice exit power — e.g. `plot_fel_compare.py tdsase.diag.txt
AramisTDSASE.out.h5`. Add `--fld <tier>-final.fld.h5 <GenesisRoot>-final.fld.h5` for a
second figure overlaying the FINAL FIELD (on-axis amplitude and unwrapped-phase
lineouts, plus the transverse map of the complex difference) — the panel that shows
what a phase-dominated tier number is made of: for `tier1_unavg`, an amplitude overlay
indistinguishable by eye, a −0.03 rad on-axis phase offset, and the 6.9e-2 difference
localized at the beam center, i.e. the freshly radiated model-dependent part. The slice count comes from the Genesis file, so steady-state and
time-dependent tiers both work.

The `weight_split` check is Fortran against itself and exists because Genesis cannot test
the weighted paths: its dump format has no weights, so every cross-code comparison sees
the uniform case, where a bug like using one particle's weight for all is invisible.
Collective observables are linear in the weights, so the split run must reproduce the
unsplit one to round-off; the weight-dropping mutation fails this check at 2.4e-1.

### The tier2_bmad difference is a transport model difference, located and priced

The divergence localizes to the quadrupoles: after the first quad the bunching phase steps
by 2.1e-3 rad while the transverse beam sizes still agree to 2.5e-11. The mechanism:
Genesis advances theta through an interlude element as a single step whose path-length
term `(px^2+py^2)/(2 gamma^2)` is sampled once, at mid element (transverse half step, then
the theta step, then the second half; `TrackBeam::track` + the `BeamSolver` drift case).
Bmad's z advance integrates the same term exactly through the quad map. In a quad px
changes substantially across the element, so midpoint sampling differs from the integral
at the ~10 percent level of the ~0.01 rad px^2 contribution -- the observed 2e-3 rad per
quad. Through exponential gain and twelve quads this compounds to the 1.3e-2 power
difference at saturation.

The proof is `tier2_genesis`: the identical build with only the interlude model swapped to
the transcribed Genesis step collapses the difference by six orders of magnitude. Nothing
else changes between those runs, so nothing else contributes at that level. Where the two
transport models differ, Bmad's is the better one (the exact integral); the benchmark's
job is to price the difference against Genesis, not to prefer Genesis's answer.

Residual budget for `tier2_genesis`'s own 1e-8: rounding differences in the interlude
theta advance. Genesis evaluates `ks*(1 - 1/beta_z)` per particle -- a ~4e-9 cancellation
carrying ~1e-16 absolute rounding -- inside its RK4 bookkeeping; the transcription
evaluates the same expression but sums the step differently (the RK4 collapses exactly
when the slope is theta independent, and is collapsed). Multiplied by `ks*L` this is
microradians of per-particle phase noise per interlude, amplified through gain. The
per-particle theta medians tell the same story: 6.4e-14 rad for the typical particle, with
a chaotic separatrix tail (see below).

### theta is reported, not checked

The final per-particle theta difference is printed with max, rms and median but does not
check the comparison. theta is a bucket phase: at saturation neighboring trajectories
separate exponentially, so the worst-particle difference measures Lyapunov amplification,
not implementation quality (the brief's 9.1 warning about chaotic growth, met in
practice: `tier2_genesis` has median 6.4e-14 rad against max 2.0e-5). Its collective
effect is checked, through the bunching curve and the final field.

### The harness bites

Verified by mutation, the FINDINGS.md 4.1 discipline: dropping the factor 2 on the source
term fails tier1 at 3.5e-1; dropping the conjugation in the field gather fails at 2.8e-2;
and using one particle's weight for every deposition, invisible to every Genesis-based
tier, fails the split-weight check at 2.4e-1. One sensitivity was knowingly traded away
with the move to Bmad constants: the near-degenerate replacement of `sqrt(faw2)` by `faw`
in the deposition, a 1.5e-10-level effect that the transcription-era 1e-10 check caught,
now sits below the 2e-6 constants floor and is not detectable by the Genesis comparison.
Recovering that class of sensitivity is a job for Fortran-vs-Fortran regression baselines
(the weight_split pattern), not for tighter Genesis checks.

The time-dependent checks were mutation-tested the same way: a slippage rotation that
never fires fails td1 at 1.0e10 (elementwise per-slice power); zeroing the wrong slice at
rotation fails td1 at 1.0e9; and dropping the drift autophasing fails td2_genesis at
1.9e8. The fourth sensitivity was demonstrated live rather than by mutation: guarding the
end-of-lattice autophasing the way the mid-lattice case is guarded -- one wavelength short
on the final step -- failed td1 at 0.84 of the final field during development, which is
how the unguarded `+1` in Genesis (FINDINGS.md 7.1) was found. The elementwise per-slice
power comparison is what makes these loud: a one-slice misalignment puts finite power
against near-vacuum slices, so the failure signature is orders of magnitude, not percent.

The element-parameter path was mutation-tested the same way: a spurious `1/sqrt(2)` on
the helical aw (the rms-convention error, exactly the mistake a translator would make)
kills the gain outright and fails tier1 at 1.0 relative; deriving helicity from the
wrong attribute (element key instead of `field_calc`) fails identically; and removing
the `b_max` assertion is caught by the refusal check — the lattice still dies, but on an
unrelated downstream message instead of by name, which the check's grep rejects.

The thread-independence check bites too: reintroducing a shared source accumulator across
slices (the exact state of the code before deliverable 5) puts the 1-thread and 8-thread
runs apart by 7.0 relative in power -- while the mutated 1-thread run is IDENTICAL to the
pristine one. That is the defining property of this bug class: invisible to every
single-threaded check, including all seven Genesis tiers, and caught only by comparing
across thread counts.

## Parallelism (brief 4.3): OpenMP over slices, bit-identical by construction

Slices are independent within an integration step: the only cross-slice operations are
slippage (an index rotation, applied serially between steps) and the per-beam `phi0`
advance (one scalar, computed before the loop). The slice loops of `fel_track_und_step`
and `fel_track_interlude_genesis` are therefore plain `parallel do` regions, and the work
of this deliverable was making the per-slice step safe to run concurrently:

- The FFTW plan cache in `wavefront_mod` is **threadprivate**: each thread owns its plans
  and its aligned work buffer (the change the cache's design note anticipated). Plan
  *creation* stays inside a named critical section -- FFTW's planner is globally
  serialized -- while plan *execution* touches only the calling thread's buffer.
- The field-solve kernel (`fel_k2`) is built **once, serially**, by
  `fel_field_kernel_init` before any parallel loop, and is read-only ever after;
  `fel_field_step` errors on a mismatch rather than rebuilding, so nothing writes module
  state inside a parallel region. The source accumulator is a local of `fel_field_step`,
  one per invocation.
- The Bmad-seam interlude loop stays **serial** at the slice level: `track1_bunch`
  parallelizes over particles internally (`track1_bunch_hom`, on by default via
  `global_com%mp_threading_is_safe`), so the threads are already busy inside each call
  and an outer parallel loop would nest.

Because each slice's arithmetic is identical regardless of which thread runs it -- no
cross-slice reductions exist in the step loop -- the result is **bit-identical across
thread counts**, and the harness checks that: the time-dependent single-segment
configuration reruns with `OMP_NUM_THREADS=8` against the 1-thread run, requiring the
diag file byte-equal and every dump dataset exactly equal. (Whole-file `cmp` of HDF5
would false-alarm on object-header timestamps; the comparison is per dataset.)

Measured scaling, full 6-FODO line, 32 slices x 2048 particles (Apple Silicon, debug
build):

| Threads | Wall time | Speedup | Efficiency |
|---|---|---|---|
| 1 | 129.6 s | 1.00 | -- |
| 2 | 74.4 s | 1.74 | 87% |
| 4 | 46.4 s | 2.79 | 70% |
| 8 | 32.6 s | 3.97 | 50% |

It saturates near 4x at 8 threads. The serial fraction per step is real: the per-slice
diagnostics reduction (32 grids of 255^2 every record) and the slippage rotation run
serially between the parallel regions, two parallel regions are spawned per integration
step, and on this machine 8 threads includes efficiency cores. With 32 slices there is
also little schedule slack -- 4 slices per thread at 8 threads. Production-size runs
(hundreds to thousands of slices) have more parallel work per serial byte, so this is
the floor of the scaling, not its ceiling.

### Head to head against Genesis at 12 workers

One command reproduces this measurement on any machine:

```
./bsim/fel/tests/run_perf_benchmark.sh
```

It detects the performance-core count (`--workers N` to override), sizes the window to
4 slices per worker so Genesis does not pad, finds the MPI launcher that matches the
Genesis binary's own linkage (an OpenMPI `mpirun` aborts an MPICH-linked genesis4 on
sight), times all four runs with one external clock, and refuses to print a table
unless the two codes' answers agree at the documented seam level first.

Measured on a 48-slice run (slice count divisible by the worker count, so Genesis does
not pad its window), full line, 2048 particles per slice, identical starting dumps,
production/release builds on both sides, M3 Max with 12 performance cores; the answers
agree at the documented seam level (4.0e-2) before any timing is quoted:

| | Serial | 12 workers | Parallel speedup |
|---|---|---|---|
| Genesis 1.3 v4 (MPI) | 152.7 s | 35.8 s (12 ranks) | 4.3x |
| `fel_track_test` (OpenMP) | 123.6 s | 19.4 s (12 threads) | 6.4x |

The tracker is 1.2x faster serial and 1.85x faster at 12 workers. The parallel gap is
the design brief's section 4.3 bet showing up in a measurement: Genesis's 12-rank run
spent 68 s of system time on the per-step MPI slippage ring exchange and diagnostics
gather, where this code's slippage is an index rotation in shared memory and costs
nothing to communicate.

## Shot noise under weights (brief 6.2): the loader's noise is physical, and checked

(Physics: manual `sec:loading`.)

A slice of current I represents `N_lambda = I*slice_spacing/(e*c)` real electrons, and
physical shot noise means `<|b(h)|^2> = 1/N_lambda` per harmonic. The built-in loader
imposes it Fawley style, transcribed from Genesis's `ShotNoise::applyShotNoise` and
generalized to per-particle weights in one substitution: the per-beamlet electron count
in the amplitude is the beamlet's REAL charge over e, not Genesis's slice-uniform
`ne/mpart` (identical for uniform weights). The algebra then gives `1/N_lambda` for any
cross-beamlet weight distribution with no correction factors (FINDINGS.md 7.6). Weights
must stay uniform within a beamlet -- the quiet cancellation is per beamlet. Genesis's
silent `nbl < 1` clamp is kept but counted and warned.

The N_eff discipline (brief 6.2's trap): the loader reports per-slice `N_lambda` and
`N_eff = (sum w)^2/sum w^2`, and refuses to impose noise on a slice whose pre-noise
quiet floor `max_h |b(h)|^2` is not far below the target -- swept over every harmonic
the beamlet structure can resolve (`1..nbins-1`), not just the imposed ones, because an
unquiet weight pattern can park its floor on a harmonic the imposition never touches
(FINDINGS.md 7.7, found by this guard's own mutation test).

Two permanent checks in the harness:

- `check_shot_noise.py` (self-referenced, FINDINGS 6.9): many-seed loading-only runs,
  `<|b(h)|^2>*N_lambda` against 1 within 5/sqrt(n) for harmonics 1-3, uniform AND
  0.25x/1.75x alternating beamlet weights. Measured 1.03 in both modes.
- `check_sase_startup.py` (cross-code): dark start, one segment, each code generating
  its own noisy beam -- fully independent loaders and RNGs -- and the mean SASE startup
  power must agree. Measured ln ratio -0.003 (0.3 percent) over 6 seeds x 32 slices.

Mutations bite: amplitudes from macroparticle count fail the statistical check at 13-20x;
a slice-uniform electron count (weights ignored) passes uniform and fails the nonuniform
mode at 1.52 (theory 1.5625); the within-beamlet weight mutation makes the guard refuse
at 2.3e-1 against a 5.3e-5 target.

These two checks are statistical by necessity (independent RNGs). The `tdsase` tier is
their deterministic complement: Genesis generates the noisy beam, writes it, and both
codes track the identical realization dark through the full line, so startup-from-noise
is also compared elementwise like any other tier.

## Distribution import (brief 10 step 10): a bunch_struct, resampled Genesis's way

(Physics: manual `sec:import`.)

The `importdistribution` equivalent: an arbitrary bunch -- arbitrary times, arbitrary
weights -- resampled into the evenly spaced, equal-population slices the FEL step
wants, by Genesis's own method, transcribed from `SDDSBeam.cpp` (`fel_import_mod`; the
class name is historical, it reads plain HDF5). The bunch comes in two ways: generated
natively from Bmad's `beam_init_struct` via `init_beam_distribution` (a `&beam_init`
namelist block -- Bmad's equivalent of Genesis's `&beam` description), or read from an
openPMD-beamphysics file via `hdf5_read_beam`. The driver can write any bunch back out
as a Genesis distribution file (`t/p/x/xp/y/yp` + charge, with `t = -tau/c` so
Genesis's `s = -c*t` reproduces this port's window position exactly) and as
openPMD-beamphysics (`hdf5_write_beam`).

Per slice: every particle inside a sampling window `dslen = slicewidth*bunch_length`
(default 0.01 of the bunch, deliberately much wider than a slice) is a candidate, and
the slice current comes from the same window -- Genesis's `count*dQ*c/dslen`,
generalized to the weight sum `c*sum(w)/dslen`, identical for uniform weights and the
only weighted generalization made. Candidates are brought to `npart/nbins` beamlet
seeds by random deletion or by Genesis's phase-space interpolation (normalize 5D to
unit rms; nearest ORIGINAL neighbor under a metric whose per-coordinate weights are
fresh random draws; child at midpoint plus uniform[-1,1] times the difference); theta
is refilled over one beamlet spacing, mirrored into nbins bins, and the deliverable-6
Fawley loader imposes shot noise with `ne = round(I*lambda*sample/(e*c))` -- shared
code (`fel_fawley_noise`), so the generator and the import stay one implementation.
Genesis's `match`/`center` transforms are NOT ported, by decision: they exist because
Genesis lattices carry no optics, so an imported bunch must be rematched by hand. A
Bmad lattice carries its Twiss and `init_beam_distribution` generates bunches matched
to it already -- the transform would be a second way to say what the lattice says.
For the same reason there is no namelist `gamma0` anywhere in the driver anymore: the
reference energy derives from the lattice's `e_tot`, after the first external user fed
a hand-rounded gamma0 against a round lattice energy and the run died mid-tracking on
the seam's backstop p0c check (two specifications of one truth was the defect;
FINDINGS.md 7.19). Facts pinned by reading: the `align*` parameters are parsed but
never used in v4, and the `shotnoise` flag is read but never consulted (the import
applies noise unconditionally, skipping only zero-current slices) -- both kept as
Genesis has them, neither transcribed as functional. one4one is out of scope: weights
supersede it.

Validation (`scripts/check_import.py`, in the harness) splits along the RNG boundary
-- exact checks where no random number enters, statistical only where one does:

| Check | Kind | Measured |
|---|---|---|
| per-slice current profile vs Genesis importing the SAME file | exact | **8.2e-13** of peak (24 slices) |
| the generated bunch carries the specified emittance | exact | ex, ey within 3e-5 of spec |
| split-weight invariance (coincident w/3 + 2w/3 copies) | exact | moments 1.8e-15, currents 7.4e-14 |
| openPMD round trip (write_opmd_file -> dist_file) | exact | moments 9.4e-19, currents 0 |
| thread determinism (1 vs 8 threads, same seed) | exact | byte-identical diag |
| slice Twiss/emittance recover the spec (mean, central slices) | statistical | beta 0.8%, alpha 1.0%, emit 1.4% |
| dark-start startup power vs Genesis, independent resampling RNGs | statistical | ln ratio +0.023 (check 0.30) |

Mutations bite, each on the check built for it: normalizing the current by the slice
spacing instead of `dslen` fails the exact current check at 6.5e-1; skipping the shot
noise fails the startup check at ln ratio -57 (a dead-quiet start); collapsing the
beamlet mirroring fails the startup check at ln ratio +4.7. One planned mutation --
refilling theta over 2pi instead of 2pi/nbins -- turned out to be an EQUIVALENT
MUTANT: under the beamlet mirroring, a uniform seed over the full turn is uniform
modulo one beamlet spacing, the quiet cancellation is untouched, and the checks
correctly pass it (FINDINGS.md 7.16). It is a convention, not a defect class; the
load-bearing neighbor (the mirroring itself) is what gets mutation-tested. (A fourth
mutation, the match transform's slope/momentum order, retired with the match
transform itself -- see below.)

Recorded improvement path: with per-particle weights the resampling is OPTIONAL -- a
direct weighted import (every bunch particle a macroparticle in its slice, no deletion,
no interpolation) has no Genesis counterpart outside one4one and is likely the better
default once validated; Genesis's O(n^2) randomly-reweighted nearest-neighbor
interpolation is the first thing worth replacing. `examples/import/` is the
self-contained demonstration: a beam_init bunch, matched by the import, tracked dark
through the full line.

## Slice migration under weights (brief 6.4)

(Physics: manual `sec:migration`.)

Genesis permits migration only under one4one, because with uniform weighting the charge
a mover carries cannot be expressed; per-particle weights dissolve the problem, and this
port's coordinates dissolve the other half: z is continuous, the slice index is derived,
and a one-slice move adjusts z by exactly `beta*slice_spacing` -- a phase shift of
exactly `2*pi*sample`, so for integer `sample` the ponderomotive phase is continuous
across the move to rounding. No wrap protocol exists to get wrong (the 4.2 decision
paying off a third time). `fel_migrate_slices` is the weighted generalization of
`Sorting::localSort`: same criterion (`atar = floor(theta/slen)`), same swap-with-last
removal, same rescan semantics; the MPI `globalSort` has no counterpart here by design
(brief 4.3). Particles leaving the window ends are dropped WITH THEIR CHARGE COUNTED and
reported per event -- Genesis discards them silently at the world edges
(`Sorting.cpp:194-195`). Migration is OFF BY DEFAULT (`migrate = T` enables): the
Genesis-comparison tiers run against Genesis without one4one, which never migrates, so
enabling it inside a transcription-level comparison would be a model difference. The
pass runs serially at a per-element stride between the parallel regions, so the
thread-independence check is untouched. Per-slice `current` and `n_eff` are appended to
the diag columns -- N_eff drifts once particles migrate (brief 6.4) and is monitored,
not assumed.

Four permanent checks (`check_migration.py`, self-referenced per FINDINGS 6.9): charge
conservation under heavy migration (a 60-m_ec² energy-spread beam, 74k moves, in-window
charge plus reported drops equals initial charge at every record; measured 1.5e-14),
exact phase continuity (the whole-beam weighted phasor obeys
`S_before = S_after + S_dropped`; measured 6.9e-15), window residency (every surviving
particle's phase inside its slice window in the final dump -- the routine's
postcondition), and no-op bit identity (a frozen-phase run with migration enabled
reports zero moves and reproduces the disabled run byte for byte).

Mutations bite, each caught by the instrument built for it: charge left behind on the
move fails conservation at 8.9e-1; z unadjusted on the move SURVIVES conservation and
continuity -- the mover cascade evaporates through accounted window-end drops, a balance
sheet that balances while the beam disappears -- and is caught by window residency (17
survivors outside); removing the high-side bounds check dies on the constructed
off-the-end mover with a bounds trap, never touching memory beyond the arrays
(FINDINGS.md 7.8).

## Wakes and space charge (brief 10 step 8): Genesis's granularity, checked

(Physics: manual `sec:wakes`, `sec:spacecharge`.)

Inside undulators and Genesis-model interludes, the collective terms are transcribed at
Genesis's granularity (`fel_collective_mod`): wakes as a per-slice energy-loss rate --
the three single-particle kernels (resistive wall via the NUMERICAL impedance of Bane &
Stupakov SLAC-PUB-10707, with AC conductivity and round/flat geometry; the undulator gap
wake, convolved with dI/ds; surface roughness via the complex-q contour) superposed and
causally convolved with the window's current profile, applied between the longitudinal
advance and the second transverse half step, theta held fixed through the kick -- and
longitudinal space charge as the per-particle `ez` in the pendulum equation, from the
per-slice radial-harmonic tridiagonal solves plus the whole-window long-range term, both
weighted (`c*w_j/slice_spacing` where Genesis has `current/npart`). The convolution is
hoisted once, as Genesis's is, and recomputed at the migration stride when migration can
change the currents; every recompute appends a z-stamped eloss block to
`<out_root>.wake.txt`, making "the wake followed the currents" a parseable fact.

Two decisions recorded here as much as in the code: the numerical impedance is a clean,
SEPARABLE routine (`fel_resistive_wall_wake`) because that computation is a future port
target into Bmad proper as a wake source; and the space-charge solver sits behind an
interface a Bmad-slice implementation can later fill -- Bmad's slice method is suspected
the better model long-term, Genesis's is transcribed now for consistency, and the
comparison between the two is an explicit future task.

Genesis's collective code carries more truncated constants than its FEL core, and each
tier's check is sized to the floor of the terms it enables:

| Term | Genesis constant | Floor | Measured tier |
|---|---|---|---|
| resistive/geometric kernels | `vacimp = 376.73` | 8.3e-7 | tdwk 8.7e-7 |
| roughness coefficient | `e = 1.6e-19`, `eps0 = 8.854e-12` | 1.4e-3 of that kernel | (inside tdwk) |
| long-range space charge | `eps0 = 8.85e-12` | 4.7e-4 | tdsc 2.4e-4 |

Self-referenced checks (`check_collective.py`): on a cold dark beam the wake is the only
energy channel, and every record's `d<gamma>` must equal `eloss*dz/m_electron` exactly
(measured 8.6e-11 against 4e-4 kicks, in gamma units at the time) with the energy
spread invariant under the uniform
kicks (4.9e-13, after the diagnostics moved to a two-pass variance -- the one-pass form
hid sigma-scale cancellation noise for five deliverables because nothing ever moved the
mean; FINDINGS.md 7.9); and under heavy migration the eloss blocks must multiply and
change (49 blocks measured).

Mutations bite, each on its named check: reversing the convolution's causality (wake
collected from trailing charge) fails tdwk at 7.0e-3 against its 8.7e-7 pristine;
flipping the sign of `ez` fails tdsc at 9.1e-1; removing the migration-stride recompute
leaves one eloss block and fails the stale-wake structural check.

## Bmad element wakes across the whole bunch (brief 10 step 11)

(Physics and conventions: manual `sec:seamwake`.)

Elements carrying Bmad `sr_wake` definitions -- pseudomodes and the tabular `z_long`
(binning + FFT) -- act across the WHOLE time window. For wake elements only, the seam
concatenates all slices into one bunch in global window coordinates,

```
z_global = z_local + beta * (islice-1) * slice_spacing
```

-- higher slice index is the window head (larger Bmad z). The formula is the slice-
migration invariant run backward (a mover's z shifts by exactly -atar*beta*spacing),
and the direction is triple-pinned: by that invariant, by the deliverable-8
convolution (eloss collects `current(is+i)`, the wake trailing its source), and
empirically by the causality check. The deliverable-11 goal guessed the opposite sign;
the checks corrected it. Interlude elements pass through Bmad's own `track1_bunch`
(wake applied at `ds_wake`, Bmad's once-per-passage convention); FEL wigglers carrying
`sr_wake` get one whole-window kick at the step nearest mid-element via
`track1_sr_wake` directly -- a pure kick, no transport, with z rescaled by
`beta_new/beta_old` to hold theta (the same convention as `wake_on`'s energy kick).
Bmad's wake machinery is used AS IS: the seam supplies global z and charges, Bmad
supplies the physics (`ix_z(1)` is the bunch head at largest `vec(5)`;
`wake_mod.f90`). Refused by name: lr (multi-bunch) wakes; a pseudomode `z_max` or a
`z_long` table extent `z0` shorter than the window (Bmad would kill the bunch
mid-run); a Bmad drift cannot carry a wake at all (use a pipe -- and the driver now
stops on ANY `bmad_parser` error rather than running a partial lattice, found when an
example's drift wakes silently never attached). Wakes are resolved through LORDS
(`pointer_to_wake_ele`): a wake on a superimposed or split element lives on the lord
and `ele%wake` is null on its slaves -- checking `ele%wake` directly was the first
shipped version's hole, found by a user's lattice whose lord wakes fell through to the
per-slice path (the zero-charge INFO spam returning was the symptom). Zero-length wake
elements are kept in the walk for the same reason (a wake on a marker-like element is
a standard Bmad idiom); zero-length elements without wakes are skipped as before.

Checks (`scripts/check_seam_wake.py`, self-referenced per FINDINGS 6.9, every wake
measurement an A-B difference against a bit-identical no-wake run on a one-step
wiggler so the FEL evolution cancels exactly):

| Check | Measured |
|---|---|
| constant-pseudomode closed form (W = amp, self = W(0)/2), per-slice means | **6.2e-10** |
| causality: kick ahead of ALL charge (must be exactly zero) | **0.0** bitwise |
| d8 direction cross-check: `wake_on` marks the same affected mask | agrees, 0 ahead |
| z_long vs first-principles particle convolution of the same table | **3.4e-8** |
| resolved-beam z_long vs `wake_on`, same kernel, per-slice bound derived from Genesis's half-slice head deficit | 0.55 of bound |
| split-weight invariance of the kick profile | 1.8e-10 |
| thread determinism with wake elements | byte-identical |
| lord resolution: a wake on a superimposition-split element applies exactly once | **1.8e-10** closed form |

The kernel bridge: `write_wake_kernels` exports the deliverable-8 Bane-Stupakov
kernels (eV/(m electron), s = 0 rows carrying the Bane self-slice half factor), and a
`z_long` table built from them (sign-flipped -- the d8 kernels are stored as SIGNED
energy loss, Bmad's table is positive-decelerating -- unhalved at s = 0, causal side
z < 0, zero-padded past the window) gives the SAME physical wake through two
independent implementations. Measured on the full 96-slice SASE line with a 0.5 mm
copper chamber on every element (`examples/bmad_wake/`): exit mean energy drop
-2.324 (Bmad z_long, once per element) vs -2.308 m_e c^2 (`wake_on`, per step), the
interior per-slice profiles agreeing to **0.7%**. On a sub-slice-width bunch the two
diverge by design: `wake_on` convolves slice-density (Genesis's model, linearly
interpolated with a zero pad past the head -- half a slice of charge missing at the
window head), `z_long` bins actual particles; the derived per-slice bound
0.5/(slices ahead + 0.5) brackets the measured differences at half its size.

Mutations bite: flipping the head/tail direction fails the closed form at 22 and
causality at 2.8e-4; dropping the slice offset (all slices coincide) fails both;
unwiring the charge fails the closed form at exactly 1.0 (kicks vanish).

## The unaveraged verification mode (brief 6.6, step 13): fc measured, not assumed

(Physics, conventions, and provenance: manual `sec:unaveraged`. MINERVA is the
production existence proof — `minerva-code-analysis.md` — and a statistical
~1e-3-class reference only, never a bit check.)

`fel_tracking = 1` — a per-element lattice attribute, so unaveraged segments mix with
averaged ones in a single line — integrates the particles through the undulator's
REAL field: the full Newton-Lorentz quiver, RK4 at `fel_steps_per_period` (default
20; MINERVA's envelope), with the radiation field as a Strang-split kick and a source
built from the actual quiver current. Nothing from the averaged coupling path appears
in it (the harness greps `fel_und_coupling|faw` out of the module), so the averaged
mode's inputs become measurements. Entry/exit are sin² amplitude ramps
(`fel_ramp_periods`, default 2) with continuous slope, so the quiver vanishes at the
segment ends where the averaged and unaveraged momentum conventions coincide; the
beam carries a `quiver_in_px` convention flag that every averaged/seam entry asserts
— in a mixed line those handoffs happen at real internal boundaries, and the sandwich
check exercises them (ledger conserved on the middle segment at 3.3e-4 of turnover,
mixed-vs-averaged exit price 2.9e-2 ln, a wake on the unaveraged segment refused by
name). Parallel over slices with the averaged path's guarantees (results are
bit-identical across thread counts; the harness checks it); wakes/space
charge/element wakes are refused by name. Each run writes `<out_root>.ledger.txt`.

Measured (`check_unaveraged.py`, in the harness — self-referenced or closed-form):

| Check | Measured | Check level |
|---|---|---|
| energy ledger: max d(E_beam+U_field) over field-energy turnover | **1.0e-5** | 1e-4 |
| ledger internal (kick-side vs realized beam change) | 1.0e-5 | 1e-4 |
| ballistic dark run: max dgamma (B does no work) | **exactly 0** | 1e-12 |
| ramp handoff: emittance ratio − 1 / orbit shift / mean-px shift | 9.8e-11 / 1.9e-9 m / 3.6e-14 | 1e-6 / 1e-7 / 1e-6 |
| fc planar (h=1) vs closed-form JJ: 0.75051 vs 0.75095 | **5.9e-4** | 5e-3 |
| fc helical (h=1) vs aw: 0.84803 vs 0.84853 | **5.9e-4** | 5e-3 |
| fc planar h=3 vs closed form: 0.21287 vs 0.21302 | **7.1e-4** | 5e-3 |
| convergence fc(10/20/30 steps per period) | 0.750582 / 0.750505 / 0.750502 | 5e-4 on 30 vs 20 |
| full-segment gain curve, unaveraged vs the averaged default (bmad_standard maps), ln ratio at exit | **2.7e-3** | 0.2 (priced) |

Two merit choices over the references' conventions, both measured: the source
deposits with `1/u_s` (not Genesis's averaged `1/gamma`), making kick and source
exact energy duals — the ledger tightened 65× to the gamma-pz round-trip floor, while
the period-averaged limit moves only ~5e-9; and the magnetic push is explicit RK4
because fourth order is what makes fc measurable at 6e-4 with 20 steps/period — its
non-symplecticity is priced at gamma exact / emittance ≤ 3.3e-6 over the longest
benchmark segment (266 periods, 5320 steps), with the ballistic check standing watch
should production-length unaveraged runs ever appear.

The JJ Bessel factor and its h=3 counterpart EMERGE from the raw dynamics at 6e-4 —
the first independent check on `fel_und_coupling`'s closed forms, and on the
harmonic-load rule (the beamlet quiet start cancels every harmonic below `nbins`, so
`nbins = 8` is quiet at h=3 and no quadrature load is needed for this measurement).
The h-probe is just the same undulator with the field at `lambda1/h`: the mode is
harmonic-agnostic.

Mutations bite, each on its named check: flipping the E·v kick sign fails the ledger
at 2.0 (energy created) and the gain comparison at 0.69 — while the magnitude-blind
fc checks PASS it, which is why the ledger is check zero; a hard-edge entry
(`fel_ramp_periods = -1`, the explicit sentinel) fails the orbit-handoff check at 3.2e-5 m (19 sigma of the
probe beam; note the exit momentum re-absorbs the quiver at integer-period lengths,
so the ORBIT, not the exit mean px, is the reliable instrument — FINDINGS.md 7.24);
skipping the exit handoff flag is refused by name at the first seam element.

## The saturation demo: one practical SASE case, three trackers, one clock

```
bsim/fel/examples/saturation_demo/run.sh
```

Every input of the demo is a real file in `examples/saturation_demo/` (Genesis decks,
Bmad namelists, the two-line unaveraged wrapper lattice) -- read them, edit them,
rerun; the script is a thin clock-and-check runner. The expert's question — does this work for a practical case? — answered on Genesis's
own Benchmark1-SASE configuration run to saturation: the full 57 m 6-FODO Aramis line,
dark start, growth from shot noise alone, 96 slices × 2048 particles, all three
trackers fed IDENTICAL initial dumps and each given the machine's full performance-core
count. Wall times come from one external clock; the exit answers must agree at the
documented levels before any timing or report is produced. The run ends with a
multi-page PDF summary (`tests/scripts/report_fel_saturation.py`), regenerated from the
run's own files: a cover carrying the timing and agreement tables, then gain curves,
pulse structure, beam evolution, and the energy accounting -- every figure with the
paragraph that explains how to read it. The energy page shows per tracker the energy
the beam gave against where it is now (still in the window vs slipped out forward),
and the unaveraged ledger closing exactly along z.

Measured (M3 Max, 12 performance cores, production builds both sides):

| | wall | exit total power | vs Genesis |
|---|---|---|---|
| Genesis 1.3 v4, 12 MPI ranks | 38.0 s | 3.381 GW (4.52 GW peak at 56.8 m) | — |
| Bmad averaged (`bmad_standard` default), 12 threads | 30.2 s | 3.380 GW | **rel 4.9e-4** |
| Bmad unaveraged (`fel_tracking = fel_unaveraged`), 12 threads | 1143.2 s | 6.25 GW | ln ratio +0.62 |

The averaged mode tracks Genesis through eight decades of z and three of power to
**4.9e-4** at saturation — the ~4e-2 seam-transport difference the benchmark tiers
price is invisible here because saturation self-limits the power. It is also 1.26x
faster than Genesis at equal cores, each code computing its own in-run diagnostics —
and ours are the full 6x6/4x4 moment sets of the stats file where Genesis's are
scalar columns (the diag/stats fusion in the diagnostics section is where that speed
came from). Spontaneous radiation is OFF in the demo on both sides, matching Genesis's
`&sponrad` default; turning `radiation_damping`/`radiation_fluctuations` on costs the
beam ~1.3e-4 of its energy over the line (~0.5 rho of accumulated detuning -- the size
that moves hard-X-ray saturation, and the reason the switches exist).

The unaveraged mode — an independent integrator with fc/JJ nowhere in its inputs,
paying its documented ~32x cost — reproduces the startup (coherent shot-noise
radiation matches both codes to ~8%), the gain curve shape, and the saturation
location (56.2 vs 56.8 m), and rides ~2%/m above the KMR codes through the
exponential regime, while its beam gives up ~14x more energy. The energy difference is
understood and measured (FINDINGS.md 7.27): both models emit spontaneous shot-noise
radiation of the same magnitude, but the averaged/KMR model does not DEBIT the beam
for it (its step adds 2S to the field and kicks with E, so the 4|S|^2 part of the
field energy is created — measured factor 134 on a dark segment), while the unaveraged
mode conserves energy by construction and therefore pays. The unaveraged mode's
captured spontaneous loss agrees with the analytic rate (2/3)r_e gamma^2 ku^2 aw^2
restricted to the grid's angular acceptance to 8% — the same formula Genesis's own
optional &sponrad module uses. Neither model yet carries the ~90% of spontaneous
power radiated outside the grid acceptance (the named follow-on). The energy panels
are where all of this is visible (the budget panel:
Genesis and the averaged mode both keep ~72% of the beam's energy in the window at
exit, the unaveraged mode keeps 10% -- its beam gave 6.6e-8 J against their 4.7e-9,
the difference radiated and slipped out forward; also the faster energy-spread
growth), and "slipped out forward" is bookkept, not asserted: the unaveraged ledger
banks the energy of every slice the slippage zero-fill discards (U_escaped) and the
deposit's own |src|^2 (U_spont, the one term the kick/deposit duality does not charge
to the beam -- physically the substep's spontaneous emission), so the time-dependent
books close EXACTLY: E_beam + U_window + U_escaped - U_spont conserved to 2.9e-3 of
turnover over the whole demo (8.0e-6 on the harness configuration, where it is a
standing check at 1e-3; the demo's 6.6e-8 J = 7.9e-9 held + 5.9e-8 escaped - 1.1e-9
spontaneous credit). The figure's dotted curve is that closure drawn on the budget
panel. Wakes would be a second, unbanked exit channel from the beam -- the ledger
exists only in the unaveraged mode, where wakes and space charge are refused by name. Measured independent of
particle count (1024/2048/4096) and steps-per-period (20/40): physics, not statistics
or resolution. A dark segment with real shot noise isolates it: same in-window noise
power in both models, 20x the beam-side energy cost in the unaveraged one. The absolute
calibration of that channel (the in-band, in-grid-acceptance fraction of undulator
radiation) is the named future check; until then the demo prices the difference at
|ln| <= 1.0 at exit, measured 0.615.

Two more unaveraged end effects were found and dispatched on the way to this figure:
the ramps' slippage deficit (3.3 rad of optical phase per end — compensated exactly by
the built-in handoff phase jump, FINDINGS.md 7.26; tier1_unavg's theta median fell
from 6.6 rad to 6.4e-2) and their reduced coupling length (~2% ln per segment at the
default 2-period ramps, real field physics, priced and left visible).

The report is banked at `fel-benchmark-plots/sat-demo-report.pdf` in the project root.

## Diagnostic output: the stats file, dumps at elements, the escaped-field bank

The production statistics live in `<out_root>.stats.h5` (manual sec:stats), in FIXED
Bmad units (m, rad, eV, s, C, J, W -- units attributes are documentation, never
load-bearing): per-record, per-slice arrays in the Genesis4 visualization layout
(record index first under h5py), with beam datasets named exactly as
`bunch_params_struct` components. Those per-record datasets are SUFFICIENT statistics
-- `scripts/bunch_params_from_stats.py` reconstructs a bunch_params dict from any
(record, slice), and at element ends (where the tracker stores Bmad's own
`calc_bunch_params`, whole-window and per slice, the Tao end-of-element pattern) the
harness holds the reconstruction against the stored values. The field side stores
`wavefront_params_struct` per slice: `centroid(4)` = (x, theta_x, y, theta_y),
`sigma(4,4)` Wigner moments, energy, power, on_axis_intensity, emit_x/y = sqrt(det)
(= M^2 lambda/4pi), with `angle_moments_valid` marking where the FFT-costed theta
rows were filled (element ends and bank time -- the `twiss_valid` pattern). Pulse
values are pooled downstream; the file stays raw. `dump_beam_at` / `dump_field_at`
dump Genesis-format files at named elements (Bmad locator syntax; unknown names
refused). `keep_escaped_field` banks every slice slippage transmits out of the window
(`-escaped.fld.h5`, with per-slice wavefront_params and z_transmit) and reconstructs
the FULL PULSE at the exit plane (`-pulse.fld.h5`) by free-space propagation at
finalize -- transmitted light is fixed information and never re-interacts, so
whole-pulse statistics use the ABCD map on the banked moment matrices and never
propagate numerically.

Measured (check_diagnostics.py, in the harness -- cross-identities, not references):

| Check | Measured | Check level |
|---|---|---|
| bunch_params reconstruction vs stored calc_bunch_params | **4.6e-8** | 1e-6 |
| banked slice energies vs the ledger's U_escaped | 6.7e-16 | 1e-12 |
| analytic (moment-map) vs FFT-propagated rms at exit, 267 slices | **8.4e-3** | 2e-2 |
| pooled pulse sigma, analytic vs numerical routes | 3.4e-3 | 2e-2 |
| stats/escaped/pulse dataset-identical, 1 vs 8 threads | exact | -- |
| unknown dump element | refused by name | -- |

diag.txt is untouched (the Genesis-comparison instrument): every benchmark tier
reproduces bit for bit -- including through the diag/stats FUSION: the per-record
stats loop also evaluates the diag instrument (the identical fel_field_diag and
fel_slice_diag calls per slice, so each slice's arithmetic is unchanged), and the
diag writer only prints. That retired the formerly SERIAL per-record diag sweeps
(all 96 field planes and 96 particle slices, every record), which were worth more
than the whole stats machinery costs: measured on the saturation demo's averaged run
(96 slices x 2048 particles, 255^2 grid, 12 threads), the pre-stats baseline was
36.5 s and the run WITH full statistics is ~32 s -- the diagnostics deliverable made
the tracker 12% faster, net. Per-slice element-end twiss is evaluated through Bmad's
own calc_emittances_and_twiss_from_sigma_matrix FED FROM the already-computed
per-record moments (an element end always coincides with its last record), not by
re-summing particles.
Known scaling limit, named for the follow-on: the stats accumulate in memory and write
once (demo: 64 MB); a tens-of-thousands-of-slices hard-X-ray window wants chunked
incremental writes instead.

## Two polarizations: vector radiation, tilt honored, the crossed undulator

The radiation carries (Ex, Ey) when any FEL element is TILTED -- `UNDY: UNDX,
tilt = pi/2` is a y-planar undulator, standard Bmad, no new attribute -- or the seed
is y-polarized (`seed_polarization = 'y'`); otherwise Ey is never allocated and the
single-component path runs untouched (every tier bit-for-bit, the compatibility
keystone). Kick and deposit act through each element's polarization 2-vector (planar:
(cos t, sin t); helical: (1,-i)/sqrt2); the unaveraged mode needs no polarization
code at all -- its real per-particle currents work against and deposit into their own
components. Dumps hold one polarization per Genesis-format file (-final-{x,y});
stats.h5's field group carries totals plus the x component, with a field/y/ group.

Measured (check_two_polarization.py, in the harness -- symmetries and physics, not
reference files):

| Check | Measured | Check level |
|---|---|---|
| rotation identity, averaged (all-y line + swapped beam == all-x line) | **4.5e-15** | 1e-8 |
| rotation identity, unaveraged | 4.1e-6 | 1e-5 (the vector path carries the betatron-current radiation the scalar convention omits -- a refinement, priced) |
| beam-size swap identity, both modes | ~1e-16 | 1e-8 |
| crossed undulator: x-field gain isolation through the y set | **1.6e-14** (avg), 1.7e-6 (unavg) | 5e-2 |
| crossed undulator: the y field lights up from carried-over bunching | Py/Px ~ 1e-2 | floor 5e-3 |
| averaged vs unaveraged crossed y-power, ln | 0.25 | 0.5 (priced) |
| unaveraged TD ledger over both components | 9.6e-5 | 1e-3 |
| 1 vs 8 threads, crossed TD | byte-identical | -- |
| helical re-anchor (vector vs scalar path, shot-noise power) | **7.2e-15** | 1e-12 |
| tilt on helical / tilt with transcribed maps | refused by name | -- |

The rotation identity uses the driver's `swap_beam_xy` check instrument (the RNG
draws its planes sequentially, so the generated beam itself is never swap-symmetric).
The elliptical follow-on (cartesian_map-derived coupling, APPLE-II, the Ming-Xie
ladder, B0's em_field_calc unification) builds on this foundation.

## Spontaneous emission: the two FEL modes against Bmad's own radiation

Bmad-only, no Genesis (the averaged mode's Genesis agreement is settled by the tiers).
One short helical wiggler (`tests/bmad/spont_probe.bmad`, 40 periods, benchmark
parameters), the same beam, tracked four ways -- `check_spontaneous.py` in the harness:

| | Δγ over 0.6 m | vs analytic |
|---|---|---|
| analytic `(2/3) r_e γ² k_u² a_w²` (= the coefficient in Genesis's `Incoherent.cpp`) | 0.018369 | — |
| **Bmad `runge_kutta` + `radiation_damping`** (the same wiggler, Bmad's own tracking) | 0.018371 | **1.0e-4** |
| Bmad `runge_kutta`, radiation off | exactly 0 | the integrator conserves |
| averaged FEL mode | 1.6e-5 | 8.9e-4 — nothing |
| unaveraged FEL mode | 6.0e-4 | 3.3% |

Bmad's radiation damping reproduces the analytic rate to **1e-4**, so the reference is
independent and not in doubt. The two FEL modes then say something specific:

The **averaged (KMR) mode does not charge the beam at all**. Its step adds `2S` to the
field while the particles are kicked by `E`, so the `4|S|²` part of the field-energy
increment is created rather than taken from the beam. That is the model, not a defect
-- Genesis carries an optional `&sponrad` module (`doLoss`/`doSpread`, off by default)
precisely to add the missing loss and diffusion by hand. The check pins it as a known
zero, so a future change that starts debiting the beam cannot pass unnoticed.

The **unaveraged mode conserves energy by construction**, so its beam pays -- but only
for the radiation the grid can hold. An SVEA grid represents angles to the FFT Nyquist
θ_max = λ/2dx, and the evidence that the captured 3.3% really is acceptance-limited
undulator radiation is its SCALING: varying only the acceptance (ngrid 127/255/511 at
fixed box, θ_max = 1.6/3.2/6.4e-5 rad) moves the captured loss 0.84% -> 3.28% -> 11.4%,
a measured 13.7x against a predicted 10.5x across a 16x range in captured solid angle.
The absolute normalization sits ~3x below a dipole-limit estimate of the angular
distribution, which is that estimate's own accuracy at a_w ~ 1 -- so the test checks
the shape tightly and the magnitude loosely, which is the honest split.

RESOLVED: both FEL modes now honor Bmad's GLOBAL switches
`bmad_com%radiation_damping_on` / `%radiation_fluctuations_on` (set from the driver's
`radiation_damping` / `radiation_fluctuations` namelist conveniences; interludes always
honored them through track1). Measured, same instrument:

| with the switches on | measured | check level |
|---|---|---|
| damping: averaged vs analytic | **8.9e-4** | 5e-3 |
| damping: unaveraged vs the ramp+capture composite (0.9703 of analytic: the explicit term integrates the ramp envelope, the grid-captured self-field adds on top) | **2.1e-6** | 2e-2 |
| fluctuations: averaged / unaveraged sigma growth vs the Saldin form | 1.3e-2 / 1.8e-2 | 5e-2 |
| fluctuations: FEL form vs Bmad runge_kutta+fluctuations, ln | 0.13 | 0.25 (the references' own F-convention spread) |
| 1 vs 8 threads with fluctuations on | byte-identical | (serial per-beamlet draws in fixed slice order) |
| unaveraged TD ledger with radiation on: E_beam + U + U_esc − U_spont + E_rad | 4.6e-5 of turnover | 1e-3 |

Fluctuation kicks are ONE DRAW PER BEAMLET, exactly as Genesis: independent
per-particle kicks would break the quiet start's per-beamlet harmonic cancellation
(and fluctuations + slice migration is refused by name for the same reason). Genesis
reaches the same variance with uniform x sqrt(3) draws; ours are Gaussian -- the
physical limit. With the switches off (the default) every tier and check is unchanged
bit for bit. FINDINGS.md 7.27 has the full trail, including the measurement that first
exposed the gap (a dark segment where the averaged model's field gained 134x what its
beam paid) and the measurement-design lesson (cold beam for fluctuation growth: with a
real energy spread the sigma^2 differencing drowns in cross-covariance sampling
noise).

## The coarse-step measurement (brief 8.3)

(Summarized in manual `sec:numerics`.)

SIMPLEX's reference case integrates twelve undulator periods per step with the slice
spacing matched so slippage is one slice per step -- an order of magnitude fewer steps
than one-step-per-period. Whether that economy transfers to this integrator is the
brief's 8.3 question, measured by `run_delz_sweep.sh`: Genesis generates ONE
time-dependent initial state (32 slices of spacing `12*lambda0`, shot noise on), and the
tracker runs the full line from that same dump at `ds_step` of 1, 2, 3, 6 and 12 periods
(two-line wrapper lattices overriding the element attribute), so
every run shares one shot-noise realization and the differences are pure integration
error. Total power at the twelve undulator-segment exits, against the one-period run:

| `ds_step` | max over exits | at saturation |
|---|---|---|
| 2 periods | 1.3e-1 | 1.5e-2 |
| 3 periods | 2.5e-1 | 2.8e-2 |
| 6 periods | 5.1e-1 | 2.1e-2 |
| 12 periods | 7.0e-1 | 2.6e-1 |

The max-over-exits error lives in the exponential-gain region, where a step-size error is
a quasi-systematic gain-length shift; fitted there, convergence is roughly first order
(pairwise p of 1.7, 1.1, 0.5). The saturation error is oscillatory (post-saturation power
oscillates, so a small phase shift moves the sampled value) and fits no clean order.

The answer to 8.3: saturation power holds to ~3% up to six periods per step, and the
twelve-period matched configuration misses it by 26% -- SIMPLEX's step-size economy is
tied to its semianalytic field advance and does not transfer to this Genesis-style
integrator as-is. `ds_step` of two to three periods is the operating point here; six periods
is defensible when only saturation power matters.

## Facts about Genesis this work pinned down

- Outside undulators Genesis does not subdivide into `delz` steps: each interlude element
  is one integration step of the element's full length (`Lattice::unrollLattice` pushes
  one entry per non-undulator layout segment). The step count over the benchmark is
  12*89 undulator steps + 36 interlude steps = 1104.
- In steady state `slippage` and `phaseshift` are no-ops: `Control::applySlippage`
  returns immediately when not time dependent, and the phaseshift array is all zero
  without phase shifter elements.
- Time dependence follows from `&time` or from importing a multi-slice dump
  (`readBeamHDF5.cpp:68-77` reconstructs the window from `slicespacing` and `slicecount`,
  time on by default), with `nslice = round(slen/(sample*lambda0))` (`GenTime.cpp:70`).
- The end-of-lattice autophasing is unguarded: the last step always gets
  `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths of slippage, `+1` even with no trailing
  drift (`Lattice.cpp:191-193`; FINDINGS.md 7.1).
- The field record's rotation never appears in Genesis's outputs: `writeFieldHDF5` and
  `DiagField::calc` both unrotate on the fly, so `.out.h5` per-slice arrays and `.fld.h5`
  dumps are in time-window order, aligned with beam-slice indexing (FINDINGS.md 7.3).
- A helical undulator defaults to `kx = ky = 0.5` (LatticeParser.cpp:329), scaled by
  `ku^2` in the unroll.
- The `&importbeam` / `&importfield` namelists take full filenames and make the
  shared-start methodology possible; `&write` before `&track` produces the dumps.
- Genesis's internal field unit: `dfl [sqrt(W)] = u * dgrid * eev / (ks * sqrt(vacimp))`
  (writeFieldHDF5.cpp:70), with `vacimp = 376.73` truncated and `eev = 510998.95069`
  (identical to Bmad's `m_electron`). Equivalently `u = E*ks/(sqrt(2)*eev)` with E in
  V/m, the relation used to derive this tracker's physical-unit formulas.
- The quad transport is chromatic through per-particle `gammaz` with `foc^2 =
  k1*gamma0/gammaz`; Bmad's equivalent scaling is `k1*p0/p`. The two differ by ~4e-9
  relative at this energy (1 - beta0), far below the path-length-term difference.
