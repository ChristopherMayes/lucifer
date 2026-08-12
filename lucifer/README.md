# FEL tracker, validated against Genesis 1.3 Version 4

Deliverables 3 and 4 of the FEL port: an FEL tracker whose physics is transcribed from
Genesis 1.3 Version 4 (GPL permits transcription), embedded in Bmad's lattice machinery by
the seam of the design brief's section 4.1, and validated against Genesis over its
`benchmark/Benchmark1-SASE` configuration from bitwise-identical starting states --
single-slice steady state (deliverable 3) and multi-slice time dependence with slippage
(deliverable 4).

No space charge, no wakes, no harmonics, no slice migration. Those are later
deliverables. Per-particle weights are carried from day one (brief section 5): the
packed arrays store one, every reduction uses it, and the split-weight check below tests
the nonuniform case that no Genesis comparison can reach. The FEL step is OpenMP-parallel
over slices (deliverable 5), with thread-count independence -- bit-identical results at
any thread count -- as a gated property; see the parallelism section. The built-in
loader imposes physical shot noise generalized to weights (deliverable 6), gated
statistically and cross-checked against Genesis's independent loader; see the shot-noise
section.

## Files

| Path | Contents |
|---|---|
| `bsim/modules/fel_beam_mod.f90` | Packed particle slices in Bmad coordinates plus per-particle weight, Genesis `.par.h5` dump read/write (converting), copy-only `coord_struct` conversion, weighted beam diagnostics with `N_eff` |
| `bsim/modules/fel_track_mod.f90` | The transcribed FEL step: transverse push with natural focusing, RK4 ponderomotive advance, source deposition, FFT field solve; the rotating-record slippage machinery (`fel_slip_struct`, `fel_apply_slippage`, `fel_field_index`); plus the transcribed Genesis interlude model |
| `bsim/fel/fel_track_test.f90` | The tracker: walks a Bmad lattice, FEL steps in undulator segments, seam everywhere else, slippage schedule transcribed from `Lattice::calcSlippage`; generates its own quiet-start beam and seed field when no dumps are named |
| `bsim/fel/examples/` | Self-contained single-command examples (no Genesis, no dump files): a seeded steady-state run and a pure-SASE time-dependent run of the benchmark line, with plotting script and README |
| `bsim/fel/tests/check_shot_noise.py` | Statistical gate: `<\|b(h)\|^2> = 1/N_lambda` over many seeds, uniform and nonuniform weights |
| `bsim/fel/tests/check_sase_startup.py` | Cross-code gate: SASE startup power, our loader against Genesis's, independent RNGs |
| `bsim/fel/tests/run_fel_benchmark.sh` | The whole validation, one command |
| `bsim/fel/tests/compare_fel.py` | Comparison: three steady-state tiers plus three time-dependent tiers against Genesis, plus the split-weight invariance check |
| `bsim/fel/tests/Aramis-ss.in`, `Aramis.lat` | Genesis deck: Benchmark1-SASE steady state, modified as documented in the deck header |
| `bsim/fel/tests/Aramis-1seg.in`, `Aramis-1seg.lat` | Genesis deck: one undulator segment, importing the same dumps |
| `bsim/fel/tests/Aramis-td.in`, `Aramis-td-1seg.in` | Genesis decks: the time-dependent pair, 32 slices with shot noise |
| `bsim/fel/tests/aramis.bmad`, `aramis_1seg.bmad` | The Bmad lattices |
| `bsim/fel/tests/Aramis-td-s12.in`, `run_delz_sweep.sh` | The coarse-step measurement (brief 8.3): one shared dump at `sample = 12`, tracker runs at several `delz` |

## Running

```
cd <bmad-ecosystem>
BUILD_PRODUCTION=N ./util/conda_compile                      # builds fel_track_test
./bsim/fel/tests/run_fel_benchmark.sh [--genesis <path to genesis4>]
```

The harness runs Genesis four times (full line and single segment, steady state and time
dependent), the Bmad tracker seven times, and prints the largest relative difference of
each check. It fails loudly if the genesis4 binary is missing; there is no comparison
without it, so there is nothing to skip to. Genesis must be built with FFTW, since the
benchmark runs with `fft_fieldsolver=true` (the Bmad tracker transcribes the FFT solver;
Genesis's default ADI solver is out of scope).

## Architecture

Inside elements named `UND*`, `fel_track_und_step` advances the coupled system in steps of
`delz`, in Genesis's exact order: transverse half step, RK4 advance of (theta, gamma) with
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
impedance difference is now the floor of every Genesis comparison, and the gates below
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

Undulator segments are marked by name (`UND*`) and their FEL parameters come from the
namelist, not from the lattice file. A real FEL element type with its own tracking method
is a later deliverable; this is the smallest scheme that exercises the seam.

## Validation: seven checks, from one command

Both codes start from the same Genesis `&write` dumps, so the initial state is bitwise
identical and no loader is reproduced. Genesis records diagnostics once at the start and
once per integration step, with each interlude element being a single step; the tracker
records at the same z positions (per slice, in time-window order, for the time-dependent
tiers). Measured, on the numbers this tree was developed against:

| Tier | What runs | Largest relative difference |
|---|---|---|
| `tier1` | One undulator segment: the FEL core alone | **1.8e-6** (the impedance-constant floor; was 2.8e-11 with Genesis's constants transcribed) |
| `tier2_genesis` | Full 6-FODO line, interludes via the transcribed Genesis model | **1.8e-5** (constants floor through full gain; was 5.9e-8) |
| `tier2_bmad` | Full line, interludes via the Bmad seam | **5.0e-2** (power curve 1.3e-2) -- a measured model difference, see below |
| `weight_split` | tier1 rerun with every particle split into coincident w/3 + 2w/3 copies, against the unsplit run | **3.6e-13** (Fortran vs Fortran; constants-independent) |
| `td1` | One undulator segment, 32 slices: FEL core plus slippage (accumulation, threshold, rotation, zero fill, end-of-lattice autophasing) | **8.5e-7** (constants floor) |
| `td2_genesis` | Full line time dependent, transcribed Genesis interludes: adds the drift autophasing schedule | **2.4e-6** (constants floor) |
| `td2_bmad` | Full line time dependent through the Bmad seam | **4.1e-2** -- the tier2_bmad transport model difference with slippage interleaved |

Particle ordering is preserved by both codes (no sorting happens without one4one), so the
final dumps compare particle by particle, not just statistically, in every tier.

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

### theta is reported, not gated

The final per-particle theta difference is printed with max, rms and median but does not
gate the comparison. theta is a bucket phase: at saturation neighboring trajectories
separate exponentially, so the worst-particle difference measures Lyapunov amplification,
not implementation quality (the brief's 9.1 warning about chaotic growth, met in
practice: `tier2_genesis` has median 6.4e-14 rad against max 2.0e-5). Its collective
effect is gated, through the bunching curve and the final field.

### The harness bites

Verified by mutation, the FINDINGS.md 4.1 discipline: dropping the factor 2 on the source
term fails tier1 at 3.5e-1; dropping the conjugation in the field gather fails at 2.8e-2;
and using one particle's weight for every deposition, invisible to every Genesis-based
tier, fails the split-weight check at 2.4e-1. One sensitivity was knowingly traded away
with the move to Bmad constants: the near-degenerate replacement of `sqrt(faw2)` by `faw`
in the deposition, a 1.5e-10-level effect that the transcription-era 1e-10 gate caught,
now sits below the 2e-6 constants floor and is not detectable by the Genesis comparison.
Recovering that class of sensitivity is a job for Fortran-vs-Fortran regression baselines
(the weight_split pattern), not for tighter Genesis gates.

The time-dependent gates were mutation-tested the same way: a slippage rotation that
never fires fails td1 at 1.0e10 (elementwise per-slice power); zeroing the wrong slice at
rotation fails td1 at 1.0e9; and dropping the drift autophasing fails td2_genesis at
1.9e8. The fourth sensitivity was demonstrated live rather than by mutation: guarding the
end-of-lattice autophasing the way the mid-lattice case is guarded -- one wavelength short
on the final step -- failed td1 at 0.84 of the final field during development, which is
how the unguarded `+1` in Genesis (FINDINGS.md 7.1) was found. The elementwise per-slice
power comparison is what makes these loud: a one-slice misalignment puts finite power
against near-vacuum slices, so the failure signature is orders of magnitude, not percent.

The thread-independence gate bites too: reintroducing a shared source accumulator across
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
thread counts**, and the harness gates that: the time-dependent single-segment
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

## Shot noise under weights (brief 6.2): the loader's noise is physical, and gated

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

Two permanent gates in the harness:

- `check_shot_noise.py` (self-referenced, FINDINGS 6.9): many-seed loading-only runs,
  `<|b(h)|^2>*N_lambda` against 1 within 5/sqrt(n) for harmonics 1-3, uniform AND
  0.25x/1.75x alternating beamlet weights. Measured 1.03 in both modes.
- `check_sase_startup.py` (cross-code): dark start, one segment, each code generating
  its own noisy beam -- fully independent loaders and RNGs -- and the mean SASE startup
  power must agree. Measured ln ratio -0.003 (0.3 percent) over 6 seeds x 32 slices.

Mutations bite: amplitudes from macroparticle count fail the statistical gate at 13-20x;
a slice-uniform electron count (weights ignored) passes uniform and fails the nonuniform
mode at 1.52 (theory 1.5625); the within-beamlet weight mutation makes the guard refuse
at 2.3e-1 against a 5.3e-5 target.

## The coarse-step measurement (brief 8.3)

SIMPLEX's reference case integrates twelve undulator periods per step with the slice
spacing matched so slippage is one slice per step -- an order of magnitude fewer steps
than one-step-per-period. Whether that economy transfers to this integrator is the
brief's 8.3 question, measured by `run_delz_sweep.sh`: Genesis generates ONE
time-dependent initial state (32 slices of spacing `12*lambda0`, shot noise on), and the
tracker runs the full line from that same dump at `delz` of 1, 2, 3, 6 and 12 periods, so
every run shares one shot-noise realization and the differences are pure integration
error. Total power at the twelve undulator-segment exits, against the one-period run:

| `delz` | max over exits | at saturation |
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
integrator as-is. `delz` of two to three periods is the operating point here; six periods
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
