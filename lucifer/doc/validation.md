# Validation: what is checked, how, and at what measured level

Lucifer's physics is transcribed from Genesis 1.3 Version 4 (Genesis4) and validated against it from bitwise-identical starting states. This document is the measured record: the rule every commit obeys, the tier table and its recorded digits, the check sections that run on every validation pass, and the subsystem-by-subsystem story of what each measurement pins.

The physics itself is the manual, [`fel-physics.md`](fel-physics.md): every equation the code integrates, each subsystem's Genesis4 provenance, and which check pins it. This document carries the numbers, and links to the manual section where a section here touches physics.

(val-the-keystone-rule)=
## The keystone rule

Every commit is validated before it lands, and the tiers land on their recorded digits. A moved digit is a bug, not a new baseline. From the `bmad-ecosystem` root:

```
BUILD_PRODUCTION=N ./util/conda_compile      # debug
./util/conda_compile                         # production

# Everything below is independent, so start it all and wait. Each benchmark pass
# needs its own --work-dir once they run together.
./lucifer/tests/run_fel_benchmark.sh --results /tmp/fel-debug.txt --work-dir /tmp/wd-dbg &
./lucifer/tests/run_fel_benchmark.sh --exe $PWD/production/bin/lucifer \
        --results /tmp/fel-prod.txt --work-dir /tmp/wd-prd &
(cd regression_tests && pytest test_fortran.py --bmad-bin=$PWD/../debug/bin) &
./lucifer/wavefront/tests/run_validation.sh &
./lucifer/examples/run_examples.sh --no-figures &
wait
```

The five run at once because they share only a source tree they read: separate work directories, separate output roots, and no shared state anywhere else (checked: nothing persists FFTW wisdom or writes outside its own directory). Measured on an M3 Max, the whole thing takes 7.7 minutes with the Genesis references cached and 9.2 without, against 25 minutes when every step ran in sequence and the references were regenerated twice per run. The references are cached because they are a pure function of the reference binary and the decks that make them, under a key naming both and the pinned version besides, and each run says whether it hit or missed. Deleting `~/.cache/lucifer/genesis-refs` is always safe and costs one cold run.

| step | cached | cold |
|---|---|---|
| both benchmark passes, regression, wavefront and examples, concurrent | 458 s | 551 s |
| the regeneration below | 1 s | 1 s |
| total | 7.7 min | 9.2 min |
| the same work in sequence, before this arrangement | 25 min | 25 min |

Every section runs in every keystone. Nothing is behind a flag, and there is no shorter mode to reach for, which is deliberate: a cheap run that checks less is the thing a keystone exists to prevent.

Then regenerate the documentation that is generated, and require no diff. Both halves are the check, and running the diff alone is a trap: it then asks only whether someone hand-edited a generated file, and a page that no longer describes the code passes it. Two pages drifted for several commits under exactly that mistake, one of them missing a whole module (FINDINGS 7.43), so treat these four commands as one step.

```
python3 lucifer/tests/scripts/report_validation.py \
        --debug /tmp/fel-debug.txt --production /tmp/fel-prod.txt \
        --out lucifer/doc/generated/validation-measured.md
python3 lucifer/tests/scripts/report_api.py \
        --code lucifer/code --code lucifer/program --out lucifer/doc/generated/api.md
python3 lucifer/tests/scripts/report_examples.py \
        --examples lucifer/examples --out lucifer/doc/generated/examples
git diff --exit-code -- 'lucifer/doc/generated/*.md' 'lucifer/doc/generated/examples/*.md'
```

The measured levels in this document are written by the harness that measured them, so
a moved digit is a failing command rather than a discrepancy someone has to notice while
reading. The benchmark regenerates the example pages itself, in its `examples` section,
so the command above is the same work made explicit.

The diff names the Markdown. The figures in `doc/generated/examples/` are committed
beside the pages and deliberately excluded: a plotting-library upgrade rewrites every
one of them byte for byte without changing any physics, and a check that fails on that
is a check that gets switched off. `examples/run_examples.sh` writes them, and the
examples check asserts that each page's figure exists rather than what it contains.

The regression suite reports 52 passed and 3 skipped. Build in the `bmad-build` environment. The harness runs its Python in `bmad-fel-validate`, which also carries the `genesis4` the comparison runs against, from conda-forge and pinned in `lucifer/wavefront/tests/environment.yml`. The harness takes that one rather than searching PATH, because the levels recorded here belong to the build that produced them. `--genesis <path>` or `$GENESIS4` names a different one, which the version check then has to accept: it refuses whatever does not report the recorded version, however it was named.

The levels here were measured against the Genesis4 4.6.15 release, and the harness refuses a reference reporting any other version by name rather than running against it. That matters because the failure it prevents is silent: a different reference moves the recorded digits while every tolerance still passes. Two properties of 4.6.15 are load-bearing. It carries the CODATA electron rest energy, which is also Bmad's `m_electron`, where releases to v4.6.14 carried a value 2.14e-7 above it and loosened the transcription tiers by several percent of their own value ([](genesis4.md) states the constant and the size). And it seeds each stochastic stream from the global slice index, which this project contributed upstream, so noise realizations differ from every earlier release. The digits below were re-recorded on that migration, and which of them moved is stated in the changelog with the reason. The regeneration check refuses to absorb a reference change either way: it shows up as a non-empty diff on the generated table rather than as quietly different digits. Only environments this project created: an unrelated `devel` environment on PATH once caused an HDF5 mismatch that cost hours, and the same class of mistake is why the reference's own variant is chosen to match the HDF5 its environment already runs.

Debug and production binaries are never bit-comparable to each other. Compare like builds only.

(val-second-toolchain)=
### The same checks on a second toolchain

Every level in this document was measured on one machine, an M3 Max running macOS with the conda-forge gfortran. A second toolchain finds a class of defect the first one hides, and this project has paid for an instance: five wake cases failed on Linux because a steady-state run read past the end of a one-element vector, and the macOS allocator had been answering with a number that happened to work. `.github/workflows/lucifer-keystone.yml` builds the debug tree under Ubuntu with its system gfortran and runs the benchmark harness and the wavefront validation there. The debug build is the point, since `-fbounds-check` is what turns that read into a failure rather than a wrong number.

That run asserts the tolerances and not the digits. A different FFTW and HDF5 give different last digits by construction, so the levels recorded here stay a property of the machine that recorded them. Its first green run is the demonstration: all eleven tiers passed, four of them on digits identical to the ones above and seven moved in their last places. It takes 37.5 minutes on a four-core runner, of which the harness is 30, so it runs on demand and weekly rather than per push, and the distribution's `ci.yml` covers the build on every push.

The device section is the one section allowed to skip, and it skips there by name: a Linux build carries the refusing stub, so it has no backend to judge. The skip keys on what `lucifer` answers when run with no arguments, which names the backend the build carries and the machine under it, so a machine that has a usable backend cannot skip. A results file recording that skip cannot write the generated table, which `report_validation.py` refuses by name: that table lists every section that ran, and a run one section short would undercount them.

(val-crossidentities-not-reference-files)=
## Cross-identities, not reference files

Almost nothing here is checked against a stored expected-output file. A reference file records what the code did once. An identity records what must be true of any correct implementation, and it stays valid when the code is rewritten underneath it. So the checks are conservation laws, closed forms, independent routes to the same number, invariance under a change that must not matter (thread count, weight splitting, a no-op), and refusals by name.

A check that has never failed on a real defect is untested, so several here carry their own mutation record: what was broken deliberately, and how loudly the check noticed. And a measurement without a stated tolerance is not a check, so every number below is paired with the level it is held at.

(val-validation-from-one-command)=
## Validation, from one command

The harness runs Genesis4 six times (full line and single segment, steady state and time
dependent, plus the collective tiers), converts each chain's initial dumps to openPMD,
which is the only format the tracker reads, runs the Bmad tracker for every tier and
check, and prints the largest relative difference of each check. It fails loudly if the
genesis4 binary is missing, and likewise without an openPMD-beamphysics checkout
(`--beamphysics <path>`, a sibling of bmad-ecosystem by default), which is what performs
the conversion. There is no comparison without either, so there is nothing to skip to. Genesis4 must be built with FFTW, since the
benchmark runs with `fft_fieldsolver=true` (the Bmad tracker transcribes the FFT solver,
and Genesis4's default ADI solver is out of scope).

Both codes start from the same Genesis4 `&write` dumps, so the initial state is bitwise
identical and no loader is reproduced. Genesis4 records diagnostics once at the start and
once per integration step, with each interlude element being a single step. The tracker
records at the same z positions (per slice, in time-window order, for the time-dependent
tiers). Measured, on the numbers this tree was developed against:

```{include} generated/validation-measured.md
```

What each level means, which the digits above cannot say on their own:

| Tier | Attribution |
|---|---|
| `tier1` | The impedance-constant floor. With Genesis4's own constants transcribed it was 2.8e-11, which is what makes the floor an attribution rather than a guess |
| `tier1_unavg` | A priced model difference (sin² ramps against hard edges, no averaging, integrator structure), dominated by the final-field phase. The power curve agrees at 6.1e-3 and per-particle gamma at 3.7e-6, and theta shows a constant ~6.6 rad ramp-phase offset with only 2.8e-3 rms about it |
| `tier2_genesis` | The constants floor through full gain. With Genesis4's constants it was 5.9e-8 |
| `tier2_bmad` | A measured transport model difference, located and priced below. The power curve is 1.3e-2 |
| `weight_split` | Fortran against itself, so constants-independent |
| `td1` | The constants floor |
| `td2_genesis` | The constants floor |
| `td2_bmad` | The tier2_bmad transport model difference with slippage interleaved |
| `tdsase` | The constants floor. Exit total power agrees at 1.9e-6 |
| `tdsc` | The epsilon_0-truncation floor of Genesis4's longRange, 8.85e-12 |
| `tdwk` | The impedance floor |

Particle ordering is preserved by both codes (no sorting happens without one4one), so the
final dumps compare particle by particle, not just statistically, in every tier.

All outputs land in the benchmark's work directory (`--work-dir <path>`, and without it
a temporary directory is used, removed on success and kept on failure). To *see* any
tier rather than check it, run from a kept work directory:
`scripts/plot_fel_compare.py <tier>.diag.txt <GenesisRoot>.out.h5` overlays the two codes' power
and bunching curves, plots the checked elementwise relative power difference along the
line, and the per-slice exit power, e.g. `plot_fel_compare.py tdsase.diag.txt
AramisTDSASE.out.h5`. Add `--fld <tier>-final.wf.h5 <GenesisRoot>-final.fld.h5` for a
second figure overlaying the final field (on-axis amplitude and unwrapped-phase
lineouts, plus the transverse map of the complex difference). That panel shows
what a phase-dominated tier number is made of: for `tier1_unavg`, an amplitude overlay
indistinguishable by eye, a −0.03 rad on-axis phase offset, and the 6.9e-2 difference
localized at the beam center, i.e. the freshly radiated model-dependent part. The slice count comes from the Genesis4 file, so steady-state and
time-dependent tiers both work.

The `weight_split` check is Fortran against itself and exists because Genesis4 cannot test
the weighted paths: its dump format has no weights, so every cross-code comparison sees
the uniform case, where a bug like using one particle's weight for all is invisible.
Collective observables are linear in the weights, so the split run must reproduce the
unsplit one to round-off. The weight-dropping mutation fails this check at 2.4e-1.

### The tier2_bmad difference is a transport model difference, located and priced

The divergence localizes to the quadrupoles: after the first quad the bunching phase steps
by 2.1e-3 rad while the transverse beam sizes still agree to 2.5e-11. The mechanism:
Genesis4 advances theta through an interlude element as a single step whose path-length
term `(px^2+py^2)/(2 gamma^2)` is sampled once, at mid element (transverse half step, then
the theta step, then the second half, in `TrackBeam::track` + the `BeamSolver` drift case).
Bmad's z advance integrates the same term exactly through the quad map. In a quad px
changes substantially across the element, so midpoint sampling differs from the integral
at the ~10 percent level of the ~0.01 rad px^2 contribution: the observed 2e-3 rad per
quad. Through exponential gain and twelve quads this compounds to the 1.3e-2 power
difference at saturation.

The proof is `tier2_genesis`: the identical build with only the interlude model swapped to
the transcribed Genesis4 step collapses the difference by six orders of magnitude. Nothing
else changes between those runs, so nothing else contributes at that level. Where the two
transport models differ, Bmad's is the better one (the exact integral). The benchmark's
job is to price the difference against Genesis4, not to prefer Genesis4's answer.

Residual budget for `tier2_genesis`'s own 1e-8: rounding differences in the interlude
theta advance. Genesis4 evaluates `ks*(1 - 1/beta_z)` per particle (a ~4e-9 cancellation
carrying ~1e-16 absolute rounding) inside its RK4 bookkeeping. The transcription
evaluates the same expression but sums the step differently (the RK4 collapses exactly
when the slope is theta independent, and is collapsed). Multiplied by `ks*L` this is
microradians of per-particle phase noise per interlude, amplified through gain. The
per-particle theta medians tell the same story: 6.4e-14 rad for the typical particle, with
a chaotic separatrix tail (see below).

### theta is reported, not checked

The final per-particle theta difference is printed with max, rms and median but does not
check the comparison. theta is a bucket phase: at saturation neighboring trajectories
separate exponentially, so the worst-particle difference measures Lyapunov amplification,
not implementation quality (the warning about chaotic growth, met in
practice: `tier2_genesis` has median 6.4e-14 rad against max 2.0e-5). Its collective
effect is checked, through the bunching curve and the final field.

### The harness bites

Verified by mutation. Dropping the factor 2 on the source
term fails tier1 at 3.5e-1. Dropping the conjugation in the field gather fails at 2.8e-2.
Using one particle's weight for every deposition, invisible to every Genesis4-based
tier, fails the split-weight check at 2.4e-1. One sensitivity was knowingly traded away
with the move to Bmad constants: the near-degenerate replacement of `sqrt(faw2)` by `faw`
in the deposition, a 1.5e-10-level effect that the transcription-era 1e-10 check caught,
now sits below the 2e-6 constants floor and is not detectable by the Genesis4 comparison.
Recovering that class of sensitivity is a job for Fortran-vs-Fortran regression baselines
(the weight_split pattern), not for tighter Genesis4 checks.

The time-dependent checks were mutation-tested the same way. A slippage rotation that
never fires fails td1 at 1.0e10 (elementwise per-slice power). Zeroing the wrong slice at
rotation fails td1 at 1.0e9. Dropping the drift autophasing fails td2_genesis at
1.9e8. The fourth sensitivity was demonstrated live rather than by mutation: guarding the
end-of-lattice autophasing the way the mid-lattice case is guarded (one wavelength short
on the final step) failed td1 at 0.84 of the final field during development, which is
how the unguarded `+1` in Genesis4 was found. The elementwise per-slice
power comparison is what makes these loud: a one-slice misalignment puts finite power
against near-vacuum slices, so the failure signature is orders of magnitude, not percent.

The element-parameter path was mutation-tested the same way. A spurious `1/sqrt(2)` on
the helical aw (the rms-convention error, exactly the mistake a translator would make)
kills the gain outright and fails tier1 at 1.0 relative. Deriving helicity from the
wrong attribute (element key instead of `field_calc`) fails identically. Removing
the `b_max` assertion is caught by the refusal check: the lattice still dies, but on an
unrelated downstream message instead of by name, which the check's grep rejects.

The thread-independence check bites too: reintroducing a shared source accumulator across
slices (the exact state of the code before the parallel step landed) puts the 1-thread and 8-thread
runs apart by 7.0 relative in power, while the mutated 1-thread run is identical to the
pristine one. That is the defining property of this bug class: invisible to every
single-threaded check, including all seven Genesis4 tiers, and caught only by comparing
across thread counts.

(val-shot-noise-under-weights)=
## Shot noise under weights: the loader's noise is physical, and checked

(Physics: manual [](fel-physics.md#sec-loading).)

A slice of current I represents `N_lambda = I*slice_spacing/(e*c)` real electrons, and
physical shot noise means `<|b(h)|^2> = 1/N_lambda` per harmonic. The built-in loader
imposes it Fawley style, transcribed from Genesis4's `ShotNoise::applyShotNoise` and
generalized to per-particle weights in one substitution: the per-beamlet electron count
in the amplitude is the beamlet's real charge over e, not Genesis4's slice-uniform
`ne/mpart` (identical for uniform weights). The algebra then gives `1/N_lambda` for any
cross-beamlet weight distribution with no correction factors. Weights
must stay uniform within a beamlet. The quiet cancellation is per beamlet. Genesis4's
silent `nbl < 1` clamp is kept but counted and warned.

The N_eff discipline: the loader reports per-slice `N_lambda` and
`N_eff = (sum w)^2/sum w^2`, and refuses to impose noise on a slice whose pre-noise
quiet floor `max_h |b(h)|^2` is not far below the target: swept over every harmonic
the beamlet structure can resolve (`1..nbins-1`), not just the imposed ones, because an
unquiet weight pattern can park its floor on a harmonic the imposition never touches
(found by this guard's own mutation test).

Two permanent checks in the harness:

- `check_shot_noise.py` (self-referenced): many-seed loading-only runs,
  `<|b(h)|^2>*N_lambda` against 1 within 5/sqrt(n) for harmonics 1-3, uniform AND
  0.25x/1.75x alternating beamlet weights. Measured 1.03 in both modes.
- `check_sase_startup.py` (cross-code): dark start, one segment, each code generating
  its own noisy beam (fully independent loaders and RNGs), and the mean SASE startup
  power must agree. Measured ln ratio -0.003 (0.3 percent) over 6 seeds x 32 slices.

Mutations bite. Amplitudes from macroparticle count fail the statistical check at 13-20x.
A slice-uniform electron count (weights ignored) passes uniform and fails the nonuniform
mode at 1.52 (theory 1.5625). The within-beamlet weight mutation makes the guard refuse
at 2.3e-1 against a 5.3e-5 target.

These two checks are statistical by necessity (independent RNGs). The `tdsase` tier is
their deterministic complement: Genesis4 generates the noisy beam, writes it, and both
codes track the identical realization dark through the full line, so startup-from-noise
is also compared elementwise like any other tier.

The loader's own current profile reproduces the derived-current identity exactly, **9.9e-15** of peak, and a dropped $\sqrt{2\pi}$ fails that at 1.5. Importing the same Gaussian description agrees with it at rms **4.9e-2** of peak, which is a statistical comparison over 50k particles against a truncated-tail fit rather than an identity.

(val-distribution-import-a-bunchstruct)=
## Distribution import: a bunch_struct, resampled into slices

(Physics: manual [](fel-physics.md#sec-import).)

The `importdistribution` equivalent: an arbitrary bunch (arbitrary times, arbitrary
weights) resampled into the evenly spaced, equal-population slices the FEL step
wants, by Genesis4's own method, transcribed from `SDDSBeam.cpp` (`fel_import_mod`, the
class name is historical, it reads plain HDF5). The bunch comes in two ways: generated
natively from Bmad's `beam_init_struct` via `init_beam_distribution` (a `&beam_init`
namelist block, Bmad's equivalent of Genesis4's `&beam` description), or read from an
openPMD-beamphysics file via `hdf5_read_beam`. The driver can write any bunch back out
as a Genesis4 distribution file (`t/p/x/xp/y/yp` + charge, with `t = -tau/c` so
Genesis4's `s = -c*t` reproduces this port's window position exactly) and as
openPMD-beamphysics (`hdf5_write_beam`).

Per slice: every particle inside a sampling window `dslen = slicewidth*bunch_length`
(default 0.01 of the bunch, deliberately much wider than a slice) is a candidate, and
the slice current comes from the same window -- Genesis4's `count*dQ*c/dslen`,
generalized to the weight sum `c*sum(w)/dslen`, identical for uniform weights and the
only weighted generalization made. Candidates are brought to `npart/nbins` beamlet
seeds by random deletion or by Genesis4's phase-space interpolation (normalize 5D to
unit rms. Nearest original neighbor under a metric whose per-coordinate weights are
fresh random draws. Child at midpoint plus uniform[-1,1] times the difference). Theta
is refilled over one beamlet spacing, mirrored into nbins bins, and the Fawley
loader imposes shot noise with `ne = round(I*lambda*sample/(e*c))`, shared
code (`fel_fawley_noise`), so the generator and the import stay one implementation.
Genesis4's `match`/`center` transforms are NOT ported, by decision: they exist because
Genesis4 lattices carry no optics, so an imported bunch must be rematched by hand. A
Bmad lattice carries its Twiss and `init_beam_distribution` generates bunches matched
to it already -- the transform would be a second way to say what the lattice says.
For the same reason there is no namelist `gamma0` anywhere in the driver anymore: the
reference energy derives from the lattice's `e_tot`, after the first external user fed
a hand-rounded gamma0 against a round lattice energy and the run died mid-tracking on
the seam's backstop p0c check (two specifications of one truth was the defect,
one specification of one truth). Facts pinned by reading: the `align*` parameters are parsed but
never used in v4, and the `shot_noise` flag is read but never consulted (the import
applies noise unconditionally, skipping only zero-current slices) -- both kept as
Genesis4 has them, neither transcribed as functional. one4one is out of scope: weights
supersede it.

Validation (`scripts/check_import.py`, in the harness) splits along the RNG boundary
-- exact checks where no random number enters, statistical only where one does:

| Check | Kind | Measured |
|---|---|---|
| per-slice current profile vs Genesis4 importing the same file | exact | **8.2e-13** of peak (24 slices) |
| the generated bunch carries the specified emittance | exact | ex, ey within 3e-5 of spec |
| split-weight invariance (coincident w/3 + 2w/3 copies) | exact | moments 1.8e-15, currents 7.4e-14 |
| openPMD round trip (write_openpmd_file -> dist_file) | exact | moments 9.4e-19, currents 0 |
| thread determinism (1 vs 8 threads, same seed) | exact | byte-identical diag |
| slice Twiss/emittance recover the spec (mean, central slices) | statistical | beta 0.8%, alpha 1.0%, emit 1.4% |
| dark-start startup power vs Genesis4, independent resampling RNGs | statistical | ln ratio +0.023 (check 0.30) |

Mutations bite, each on the check built for it: normalizing the current by the slice
spacing instead of `dslen` fails the exact current check at 6.5e-1. Skipping the shot
noise fails the startup check at ln ratio -57 (a dead-quiet start). Collapsing the
beamlet mirroring fails the startup check at ln ratio +4.7. One planned mutation --
refilling theta over 2pi instead of 2pi/nbins, turned out to be an equivalent
mutant: under the beamlet mirroring, a uniform seed over the full turn is uniform
modulo one beamlet spacing, the quiet cancellation is untouched, and the checks
correctly pass it. It is a convention, not a defect class. The
load-bearing neighbor (the mirroring itself) is what gets mutation-tested.

Recorded improvement path: with per-particle weights the resampling is optional -- a
direct weighted import (every bunch particle a macroparticle in its slice, no deletion,
no interpolation) has no Genesis4 counterpart outside one4one and is likely the better
default once validated. Genesis4's O(n^2) randomly-reweighted nearest-neighbor
interpolation is the first thing worth replacing. `examples/import/` is the
self-contained demonstration: a beam_init bunch, matched by the import, tracked dark
through the full line.

(val-slice-migration-under-weights)=
## Slice migration under weights

(Physics: manual [](fel-physics.md#sec-migration).)

Genesis4 permits migration only under one4one, because with uniform weighting the charge
a mover carries cannot be expressed. Per-particle weights dissolve the problem, and this
port's coordinates dissolve the other half: z is continuous, the slice index is derived,
and a one-slice move adjusts z by exactly `beta*slice_spacing` -- a phase shift of
exactly `2*pi*sample`, so for integer `sample` the ponderomotive phase is continuous
across the move to rounding. No wrap protocol exists to get wrong (the 4.2 decision
paying off a third time). `fel_migrate_slices` is the weighted generalization of
`Sorting::localSort`: same criterion (`atar = floor(theta/slen)`), same swap-with-last
removal, same rescan semantics. The MPI `globalSort` has no counterpart here by design
Particles leaving the window ends are dropped with their charge counted and
reported per event. Genesis4 discards them silently at the world edges
(`Sorting.cpp:194-195`). Migration is off by default (`migrate = T` enables): the
Genesis4-comparison tiers run against Genesis4 without one4one, which never migrates, so
enabling it inside a transcription-level comparison would be a model difference. The
pass runs serially at a per-element stride between the parallel regions, so the
thread-independence check is untouched. Per-slice `current` and `n_eff` are appended to
the diag columns. N_eff drifts once particles migrate and is monitored,
not assumed.

Four permanent checks (`check_migration.py`, self-referenced): charge
conservation under heavy migration (a 60-m_ec² energy-spread beam, 74k moves, in-window
charge plus reported drops equals initial charge at every record, measured 1.5e-14),
exact phase continuity (the whole-beam weighted phasor obeys
`S_before = S_after + S_dropped`, measured 6.9e-15), window residency (every surviving
particle's phase inside its slice window in the final dump -- the routine's
postcondition), and no-op bit identity (a frozen-phase run with migration enabled
reports zero moves and reproduces the disabled run byte for byte).

Mutations bite, each caught by the instrument built for it: charge left behind on the
move fails conservation at 8.9e-1. Leaving z unadjusted on the move survives conservation and
continuity, because the mover cascade evaporates through accounted window-end drops, a
balance sheet that balances while the beam disappears. It is caught by window residency (17
survivors outside). Removing the high-side bounds check dies on the constructed
off-the-end mover with a bounds trap, never touching memory beyond the arrays


(val-wakes-and-space-charge)=
## Wakes and space charge: kernels, convolution and the solvers

(Physics: manual [](fel-physics.md#sec-wakes), [](fel-physics.md#sec-spacecharge).)

Inside undulators and Genesis4-model interludes, the collective terms are transcribed at
Genesis4's granularity (`fel_collective_mod`): wakes as a per-slice energy-loss rate --
the three single-particle kernels (resistive wall via the numerical impedance of Bane &
Stupakov SLAC-PUB-10707, with AC conductivity and round/flat geometry. The undulator gap
wake, convolved with dI/ds, and surface roughness via the complex-q contour) superposed and
causally convolved with the window's current profile, applied between the longitudinal
advance and the second transverse half step, theta held fixed through the kick, and
longitudinal space charge as the per-particle `ez` in the pendulum equation, from the
per-slice radial-harmonic tridiagonal solves plus the whole-window long-range term, both
weighted (`c*w_j/slice_spacing` where Genesis4 has `current/npart`). The convolution is
hoisted once, as Genesis4's is, and recomputed at the migration stride when migration can
change the currents. Every recompute appends a z-stamped eloss block to
`<out_root>.wake.txt`, making "the wake followed the currents" a parseable fact.

Two decisions recorded here as much as in the code: the numerical impedance is a clean,
separable routine (`fel_resistive_wall_wake`) because that computation is a future port
target into Bmad proper as a wake source, and the space-charge solver sits behind an
interface a Bmad-slice implementation can later fill. Bmad's slice method is suspected
the better model long-term, Genesis4's is transcribed now for consistency, and the
comparison between the two is an explicit future task.

Genesis4's collective code carries more truncated constants than its FEL core, and each
tier's check is sized to the floor of the terms it enables:

| Term | Genesis4 constant | Floor | Measured tier |
|---|---|---|---|
| resistive/geometric kernels | `vacimp = 376.73` | 8.3e-7 | tdwk 8.7e-7 |
| roughness coefficient | `e = 1.6e-19`, `eps0 = 8.854e-12` | 1.4e-3 of that kernel | (inside tdwk) |
| long-range space charge | `eps0 = 8.85e-12` | 4.7e-4 | tdsc 2.4e-4 |

Self-referenced checks (`check_collective.py`): on a cold dark beam the wake is the only
energy channel, and every record's `d<gamma>` must equal `eloss*dz/m_electron` exactly
(measured 8.6e-11 against 4e-4 kicks, in gamma units) with the energy
spread invariant under the uniform
kicks (4.9e-13, after the diagnostics moved to a two-pass variance, since the one-pass form
hid sigma-scale cancellation noise for a long time because nothing ever moved the
mean), and under heavy migration the eloss blocks must multiply and
change (49 blocks measured).

Mutations bite, each on its named check: reversing the convolution's causality (wake
collected from trailing charge) fails tdwk at 7.0e-3 against its 8.7e-7 pristine.
flipping the sign of `ez` fails tdsc at 9.1e-1. Removing the migration-stride recompute
leaves one eloss block and fails the stale-wake structural check.

(val-bmad-element-wakes-across)=
## Bmad element wakes across the whole bunch

(Physics and conventions: manual [](fel-physics.md#sec-seamwake).)

Elements carrying Bmad `sr_wake` definitions (pseudomodes and the tabular `z_long`,
binning plus FFT) act across the whole time window. For wake elements only, the seam
concatenates all slices into one bunch in global window coordinates,

```
z_global = z_local + beta * (islice-1) * slice_spacing
```

-- higher slice index is the window head (larger Bmad z). The formula is the slice-
migration invariant run backward (a mover's z shifts by exactly -atar*beta*spacing),
and the direction is triple-pinned: by that invariant, by the wake
convolution (eloss collects `current(is+i)`, the wake trailing its source), and
empirically by the causality check. The sign was guessed the other way when the seam wake was
specified, and the checks corrected it. Interlude elements pass through Bmad's own `track1_bunch`
(wake applied at `ds_wake`, Bmad's once-per-passage convention). FEL wigglers carrying
`sr_wake` get one whole-window kick at the step nearest mid-element via
`track1_sr_wake` directly -- a pure kick, no transport, with z rescaled by
`beta_new/beta_old` to hold theta (the same convention as `wake_on`'s energy kick).
Bmad's wake machinery is used as is: the seam supplies global z and charges, Bmad
supplies the physics (`ix_z(1)` is the bunch head at largest `vec(5)`,
`wake_mod.f90`). Refused by name: lr (multi-bunch) wakes. A pseudomode `z_max` or a
`z_long` table extent `z0` shorter than the window (Bmad would kill the bunch
mid-run). A Bmad drift cannot carry a wake at all (use a pipe, and the driver now
stops on any `bmad_parser` error rather than running a partial lattice, found when an
example's drift wakes silently never attached). Wakes are resolved through lords
(`pointer_to_wake_ele`): a wake on a superimposed or split element lives on the lord
and `ele%wake` is null on its slaves. Checking `ele%wake` directly was the first
shipped version's hole, found by a user's lattice whose lord wakes fell through to the
per-slice path (the zero-charge INFO spam returning was the symptom). Zero-length wake
elements are kept in the walk for the same reason (a wake on a marker-like element is
a standard Bmad idiom). Zero-length elements without wakes are skipped as before.

Checks (`scripts/check_seam_wake.py`, self-referenced, every wake
measurement an A-B difference against a bit-identical no-wake run on a one-step
wiggler so the FEL evolution cancels exactly):

| Check | Measured |
|---|---|
| constant-pseudomode closed form (W = amp, self = W(0)/2), per-slice means | **6.2e-10** |
| causality: kick ahead of all charge (must be exactly zero) | **0.0** bitwise |
| d8 direction cross-check: `wake_on` marks the same affected mask | agrees, 0 ahead |
| z_long vs first-principles particle convolution of the same table | **3.4e-8** |
| resolved-beam z_long vs `wake_on`, same kernel, per-slice bound derived from Genesis4's half-slice head deficit | 0.55 of bound |
| split-weight invariance of the kick profile | 1.8e-10 |
| thread determinism with wake elements | byte-identical |
| lord resolution: a wake on a superimposition-split element applies exactly once | **1.8e-10** closed form |

The kernel bridge: `chamber_wake%write_kernels` exports the Bane-Stupakov
kernels (eV/(m electron), s = 0 rows carrying the Bane self-slice half factor), and a
`z_long` table built from them (sign-flipped, since the d8 kernels are stored as signed
energy loss where Bmad's table is positive-decelerating, unhalved at s = 0, causal side
z < 0, zero-padded past the window) gives the same physical wake through two
independent implementations. Measured on the full 96-slice SASE line with a 0.5 mm
copper chamber on every element (`examples/bmad_wake/`): exit mean energy drop
-2.324 (Bmad z_long, once per element) vs -2.308 m_e c^2 (`wake_on`, per step), the
interior per-slice profiles agreeing to **0.7%**. On a sub-slice-width bunch the two
diverge by design: `wake_on` convolves slice-density (Genesis4's model, linearly
interpolated with a zero pad past the head -- half a slice of charge missing at the
window head), `z_long` bins actual particles. The derived per-slice bound
0.5/(slices ahead + 0.5) brackets the measured differences at half its size.

Mutations bite: flipping the head/tail direction fails the closed form at 22 and
causality at 2.8e-4. Dropping the slice offset (all slices coincide) fails both,
unwiring the charge fails the closed form at exactly 1.0 (kicks vanish).

Two bookkeeping identities ride with the seam wake. Energy bookkeeping, $d\langle E\rangle = \mathcal{E}\,\delta z$, holds to **4.4e-5** eV against 203 eV kicks. The energy spread is invariant under a uniform kick to **2.5e-7** eV, since a rigid shift moves the centroid and must leave $\sigma_E$ alone.

(val-the-unaveraged-mode-fc)=
## The unaveraged mode: fc measured, not assumed

(Physics, conventions, and provenance: manual [](fel-physics.md#sec-unaveraged). MINERVA (Freund and
van der Slot) is the published existence proof for this physics. Only its published
work was used: nothing here is compared against a MINERVA run.)

`tracking_method = fel_unaveraged` is per element, so unaveraged segments mix
with averaged ones in a single line. It integrates the particles through the undulator's
real field: the full Newton-Lorentz quiver, RK4 at `fel_steps_per_period` (default
20, floor 10 on the convergence measured below), with the radiation field as a Strang-split kick and a source
built from the actual quiver current. Nothing from the averaged coupling path appears
in it (the harness greps `fel_und_coupling|faw` out of the module), so the averaged
mode's inputs become measurements. Entry/exit are sin² amplitude ramps
(`fel_ramp_periods`, default 2) with continuous slope, so the quiver vanishes at the
segment ends where the averaged and unaveraged momentum conventions coincide. The
beam carries a `quiver_in_px` convention flag that every averaged/seam entry asserts.
In a mixed line those handoffs happen at real internal boundaries, and the sandwich
check exercises them (ledger conserved on the middle segment at 3.3e-4 of turnover,
mixed-vs-averaged exit price 2.9e-2 ln, a wake on the unaveraged segment refused by
name). Parallel over slices with the averaged path's guarantees (results are
bit-identical across thread counts. The harness checks it). Wakes/space
charge/element wakes are refused by name. Each run writes `<out_root>.ledger.txt`.

Measured (`check_unaveraged.py`, in the harness, self-referenced or closed-form):

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
deposits with `1/u_s` (not Genesis4's averaged `1/gamma`), making kick and source
exact energy duals: the ledger tightened 65× to the gamma-pz round-trip floor, while
the period-averaged limit moves only ~5e-9, and the magnetic push is explicit RK4
because fourth order is what makes fc measurable at 6e-4 with 20 steps/period. Its
non-symplecticity is priced at gamma exact / emittance ≤ 3.3e-6 over the longest
benchmark segment (266 periods, 5320 steps), with the ballistic check standing watch
should production-length unaveraged runs ever appear.

The JJ Bessel factor and its h=3 counterpart emerge from the raw dynamics at 6e-4,
the first independent check on `fel_und_coupling`'s closed forms, and on the
harmonic-load rule (the beamlet quiet start cancels every harmonic below `beamlet_size`, so
`nbins = 8` is quiet at h=3 and no quadrature load is needed for this measurement).
The h-probe is just the same undulator with the field at `lambda1/h`: the mode is
harmonic-agnostic.

Mutations bite, each on its named check: flipping the E·v kick sign fails the ledger
at 2.0 (energy created) and the gain comparison at 0.69, while the magnitude-blind
fc checks pass it, which is why the ledger is check zero. A hard-edge entry
(`fel_ramp_periods = -1`, the explicit sentinel) fails the orbit-handoff check at 3.2e-5 m (19 sigma of the
probe beam. Note the exit momentum re-absorbs the quiver at integer-period lengths,
so the orbit, not the exit mean px, is the reliable instrument).
skipping the exit handoff flag is refused by name at the first seam element.

(val-two-polarizations-vector-radiation)=
## Two polarizations: vector radiation, tilt honored, the crossed undulator

The radiation carries (Ex, Ey) when any FEL element is tilted (`UNDY: UNDX,
tilt = pi/2` is a y-planar undulator, standard Bmad, no new attribute), or when the seed
is y-polarized (`seed_polarization = 'y'`). Otherwise Ey is never allocated and the
single-component path runs untouched (every tier bit-for-bit, the compatibility
keystone). Kick and deposit act through each element's polarization 2-vector (planar:
(cos t, sin t), helical (1,-i)/sqrt2). The unaveraged mode needs no polarization
code at all, since its real per-particle currents work against and deposit into their own
components. One openPMD dump holds both polarizations as components x and y of its
mesh record, and the stats file carries field/total, field/x and field/y under the same
dataset names.

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
| 1 vs 8 threads, crossed TD | byte-identical | (exactness, no level) |
| helical re-anchor (vector vs scalar path, shot-noise power) | **7.2e-15** | 1e-12 |
| tilt on helical / tilt with transcribed maps | refused by name | (refusal, no level) |

The rotation identity uses the driver's `swap_beam_xy` check instrument (the RNG
draws its planes sequentially, so the generated beam itself is never swap-symmetric).
The elliptical follow-on (cartesian_map-derived coupling, APPLE-II, the Ming-Xie
ladder, B0's em_field_calc unification) builds on this foundation.

The helical re-anchor, which rotates the reference polarization basis and must return the same field, closes at **7e-15**.

(val-harmonic-fields-and-the)=
## Harmonic fields and the openPMD wavefront

The walk carries an ordered set of radiation fields, Genesis4's `vector<Field*>`
shape, with the fundamental always entry 1. A harmonic field h (namelist
`harmonics = 1, 3, ...`) lives at wavelength lambda_1/h and enters the physics in
exactly three places: its coupling fc(h) (the Bessel factor, zero for every
harmonic of a helical device), the harmonic ponderomotive phase h*theta in kick and
deposit, and its own diffraction kernel and escape accounting. In V/m no other
factor appears. Genesis4's per-field ks factors are its internal unit conversion,
absorbed once. Harmonic fields start dark and grow from the beam's bunching
(nonlinear harmonic generation: `examples/harmonics/`, P3 climbing nine decades to
5.1 kW under a 200 MW fundamental in 4 m), or import. A single-entry set is the
pre-harmonic walk bit for bit (the two-polarization overlay discipline again).
Harmonics with two live polarizations, and with an unaveraged element, are refused
by name, so unvalidated combinations refuse rather than guess.

Dumps speak one format, openPMD, in both directions. A field dump is an openPMD
EXT_Wavefront file (`.wf.h5`), both polarizations as complex components of one mesh
record (h5py reads them as complex128 natively), one file per harmonic with
photonEnergy identifying it, and a harmonic's file carrying `-h<h>`. The standard
document is authoritative: the harness verifies every required attribute against its
text independently of the writer. A file that is not openPMD is refused by name on
import, with the conversion command in the message, and
`lucifer/tests/scripts/convert_genesis.py` converts either kind of file in either
direction. The Python side of the round trip is openPMD-beamphysics's own
`Wavefront.from_openpmd`/`write_openpmd`, exercised against the Fortran writer and
reader every harness run.

One format rather than a choice of formats is a decision, and its cost is one
conversion step in the benchmark: the Genesis4 reference chains write their own dumps
and the harness converts them once per chain before the tiers read them. What that
buys is that no writer in the physics code has to hold a second chart, no deck names a
format, and every dump carries the per-particle weight.

Measured (check_harmonics.py, the harness's harmonics section):

| check | level |
|---|---|
| planar SS tier vs Genesis4, fundamental power | 5.3e-8 |
| planar SS tier vs Genesis4, dark harmonic-3 growth | 1.3e-4 |
| one-step dark deposit P3/P1 vs the exact Bessel sum | 3.3e-16 |
| the converter's Genesis4 view of a dump, complex values | 1.5e-16 |
| Fortran openPMD read-back, re-dump | dataset-identical |
| Python Wavefront class energy (its Genesis4 path) | 1.2e-12 |
| Python openPMD read + round trip | 1.4e-16 / exact |
| Fortran reads the Python writer's file | dataset-identical |
| 1 vs 8 threads (TD harmonic run) | byte-identical |
| six refusals (anchor, unavg, two-pol, frequency domain, no-match, format) | by name |

(val-phasing-between-segments-measured)=
## Phasing between segments: measured, then built

Between undulator segments the default is Genesis4's behavior, verified by
measurement before anything was written: scanning an inter-segment gap by fractions
of 2 gamma^2 lambda (25.8 mm per turn at the benchmark) leaves the bunching phase
entering the next segment flat on Genesis4 and on this code alike, since segments are
autophased, and the fractional drift slip never reaches the coupling. The
deliberate off-phase knob is the wiggler's own `z_offset` (standard Bmad
misalignment, anchored at the nominal position): a displacement delta shifts the
entry phase by exactly -2 pi delta/(2 gamma^2 lambda), and the same scan run as
Genesis4's own PHASESHIFTER element (which needs finite length to register, since a
zero-length one silently does nothing) reproduces our curve at 1.9e-8, and the
exit power against phi at 9.3e-9, so the phase reaches the physics identically, not
just the bookkeeping (both codes tracking the same Genesis4-written initial dumps,
on independently loaded beams the same comparison sits at 1.4e-4, which measures
the two loaders, not the phasing). No new element, no new attribute.

With `bmad_com[absolute_time_tracking] = T` in the lattice (honored per element
through Bmad's own resolver), phasing follows the real spacings instead -- the
recirculating-linac discipline: a wrong-length break visibly detunes the next
segment, and the absolute-mode gap scan reproduces the relative-mode z_offset scan
point by point. Chicanes work in both modes: the beam detours the closed bump via
the seam while the radiation drifts the chord between undulator faces (from
ele%floor, derived never entered). The arc-minus-chord delay becomes whole-slice
window rotations plus, in absolute mode, its carrier phase. examples/chicane/
flips a second segment between 1.32x gain and 0.94x absorption on 0.43 urad of
bend angle -- half a wavelength of geometric delay.

Measured (check_phasing.py, the harness's phasing section):

| check | level |
|---|---|
| re-anchor baseline: gap scan flat (rad span) | 1.5e-3 |
| z_offset knob vs -2pi delta/(2 gamma^2 lambda) | 1.5e-3 |
| cross-mode identity (absolute gap == relative knob) | 1.2e-3 |
| our knob curve vs Genesis4's PHASESHIFTER curve (shared dumps) | 1.9e-8 |
| exit power vs phi: our knob vs Genesis4's shifter (shared dumps) | 9.3e-9 |
| chicane, relative mode: geometric fraction dropped | 6.3e-6 |
| chicane, absolute mode vs independent 2D geometry | 6.7e-4 |
| 2D trace vs the small-angle closed form theta^2(2L_b/3 + L_d) | 6.4e-8 |
| unaveraged ledger closure across the chicane | 4.0e-6 |
| TD chicane: extra banked slices == floor(delay/lambda) at 3.5, 6.5, 10.5 lambda | exact |
| TD chicane 1 vs 8 threads | byte-identical |
| four refusals (open bump, genesis-model bend, oversize z_offset, first-element z_offset) | by name |

(val-spontaneous-emission-the-two)=
## Spontaneous emission: the two FEL modes against Bmad's own radiation

Bmad-only, no Genesis4 (the averaged mode's Genesis4 agreement is settled by the tiers).
One short helical wiggler (`tests/bmad/spont_probe.bmad`, 40 periods, benchmark
parameters), the same beam, tracked four ways -- `check_spontaneous.py` in the harness:

| | Δγ over 0.6 m | vs analytic |
|---|---|---|
| analytic `(2/3) r_e γ² k_u² a_w²` (= the coefficient in Genesis4's `Incoherent.cpp`) | 0.018369 | the reference |
| **Bmad `runge_kutta` + `radiation_damping`** (the same wiggler, Bmad's own tracking) | 0.018371 | **1.0e-4** |
| Bmad `runge_kutta`, radiation off | exactly 0 | the integrator conserves |
| averaged FEL mode | 1.6e-5 | 8.9e-4 of it, nothing |
| unaveraged FEL mode | 6.0e-4 | 3.3% |

Bmad's radiation damping reproduces the analytic rate to **1e-4**, so the reference is
independent and not in doubt. The two FEL modes then say something specific:

The **averaged (KMR) mode does not charge the beam at all**. Its step adds `2S` to the
field while the particles are kicked by `E`, so the `4|S|²` part of the field-energy
increment is created rather than taken from the beam. That is the model, not a defect
-- Genesis4 carries an optional `&sponrad` module (`doLoss`/`doSpread`, off by default)
precisely to add the missing loss and diffusion by hand. The check pins it as a known
zero, so a future change that starts debiting the beam cannot pass unnoticed.

The **unaveraged mode conserves energy by construction**, so its beam pays, but only
for the radiation the grid can hold. An SVEA grid represents angles to the FFT Nyquist
θ_max = λ/2dx, and the evidence that the captured 3.3% really is acceptance-limited
undulator radiation is its scaling: varying only the acceptance (ngrid 127/255/511 at
fixed box, θ_max = 1.6/3.2/6.4e-5 rad) moves the captured loss 0.84% -> 3.28% -> 11.4%,
a measured 13.7x against a predicted 10.5x across a 16x range in captured solid angle.
The absolute normalization sits ~3x below a dipole-limit estimate of the angular
distribution, which is that estimate's own accuracy at a_w ~ 1 -- so the test checks
the shape tightly and the magnitude loosely, which is the honest split.

Resolved: both FEL modes now honor Bmad's global switches
`bmad_com%radiation_damping_on` / `%radiation_fluctuations_on` (set directly in
`&fel_params`, which exposes `bmad_com` the way Tao's `&tao_params` does. Interludes
always honored them through track1). Measured, same instrument:

| with the switches on | measured | check level |
|---|---|---|
| damping: averaged vs analytic | **8.9e-4** | 5e-3 |
| damping: unaveraged vs the ramp+capture composite (0.9703 of analytic: the explicit term integrates the ramp envelope, the grid-captured self-field adds on top) | **2.1e-6** | 2e-2 |
| fluctuations: averaged / unaveraged sigma growth vs the Saldin form | 1.3e-2 / 1.8e-2 | 5e-2 |
| fluctuations: FEL form vs Bmad runge_kutta+fluctuations, ln | 0.13 | 0.25 (the references' own F-convention spread) |
| 1 vs 8 threads with fluctuations on | byte-identical | (serial per-beamlet draws in fixed slice order) |
| unaveraged TD ledger with radiation on: E_beam + U + U_esc − U_spont + E_rad | 4.6e-5 of turnover | 1e-3 |

Fluctuation kicks are one draw per beamlet, exactly as Genesis4: independent
per-particle kicks would break the quiet start's per-beamlet harmonic cancellation
(and fluctuations + slice migration is refused by name for the same reason). Genesis4
reaches the same variance with uniform x sqrt(3) draws. Ours are Gaussian. The
physical limit. With the switches off (the default) every tier and check is unchanged
bit for bit. The full trail includes the measurement that first
exposed the gap (a dark segment where the averaged model's field gained 134x what its
beam paid) and the measurement-design lesson (cold beam for fluctuation growth: with a
real energy spread the sigma^2 differencing drowns in cross-covariance sampling
noise).

(val-the-coherent-source-simplex)=
## The coherent source (SIMPLEX hybrid)

`global%source_model = "coherent"` swaps the per-particle source deposit for Tanaka's
coherent retrieval (PRAB 27, 030703 (2024), implemented from the paper): the spatially
incoherent part of the bunching, whose per-cell sampling noise makes under-populated
runs overestimate gain, is dropped from the source. The coherent part deposits as an
analytic Gaussian carrying the slice's exact source phasor, centered and tilted by
phasor-weighted moments (this port's extension) with Tanaka's kappa width fit. The
field, gather, and every diagnostic are untouched. `"deposit"` stays the default and
the reference the coherent model is measured against.

Measured (check_coherent.py, the harness's coherent-source section): at M = 128/slice the
plain deposit fakes ln P by +0.42 on a curve that truly absorbs, the coherent source
stays at 0.048 -- a 64x particle reduction at the model's own bias (1.9e-2 at large
M). Guarded by name: per-slice Gaussianity vs sampling significance (weighted by charge),
and refusals for unaveraged/harmonics/two-polarization and for dark starts -- measured
~175x startup deficit: spontaneous spatially-incoherent emission dominates SASE
startup and the coherent model drops it, so seeded runs only. Also priced once:
SIMPLEX's coarse stepping (12 periods/step) costs 2.6e-3 in |ln P| here (taper and
harmonics untested at coarse steps).

(val-the-programs-own-identities)=
## The program's own identities

Three invariances of the driver, none of them physics, each one a way a rewrite could break the tracker silently.

The library never stops. Every error returns through `err_flag` and the program decides, so `tests/lucifer_smoke_test.f90` drives the whole library with no namelist anywhere and reproduces a namelist run dataset-identically. `track_fel_line` is re-entrant, and twice in one process is bit-identical to two processes.

A windowed run composes with the full one. `global%track_start`/`track_end` bound the walk, but the schedule (slippage, autophasing, break geometry, including the end-of-lattice fixup) is always built on the full lattice. So `[start, D]` followed by `[after D, end]` from its dumps reproduces the one-shot run to the dump format's own round-trip floor, measured **3e-13** and held at 1e-10, with the walk itself bit-for-bit and a windowed run's finals equal to the full run's mid-line dumps exactly.

Refusal texts are a supported contract. stdout is otherwise for humans and carries no contract at all, which is why the suite reads files rather than scraping the screen. The one exception is deliberate: a refusal must be recognizable by name, so the suite matches the capitalized text of a refusal together with a nonzero exit status. Those texts therefore change only with their checks, in the same commit.

A fourth: an element-driven line reproduces the namelist-driven reference bit-identically at $z=0$, growing only to **3.4e-12** over the line, which is the one-ulp round trip of the $a_w$ lattice attribute. The FEL element's assertion checks refuse each malformed lattice by name.

(val-parallelism-openmp-over-slices)=
## Parallelism: OpenMP over slices, bit-identical by construction

Slices are independent within an integration step: the only cross-slice operations are
slippage (an index rotation, applied serially between steps) and the per-beam `phi0`
advance (one scalar, computed before the loop). The slice loops of `fel_track_und_step`
and `fel_track_interlude_genesis` are therefore plain `parallel do` regions, and the work
of the parallel step was making the per-slice step safe to run concurrently:

- The FFTW plan cache in `wavefront_mod` is **threadprivate**: each thread owns its plans
  and its aligned work buffer (the change the cache's design note anticipated). Plan
  *creation* stays inside a named critical section (FFTW's planner is globally
  serialized), while plan *execution* touches only the calling thread's buffer.
- The field-solve kernel (`fel_k2`) is built **once, serially**, by
  `fel_field_kernel_init` before any parallel loop, and is read-only ever after.
  `fel_field_step` errors on a mismatch rather than rebuilding, so nothing writes module
  state inside a parallel region. The source accumulator is a local of `fel_field_step`,
  one per invocation.
- The Bmad-seam interlude loop stays **serial** at the slice level: `track1_bunch`
  parallelizes over particles internally (`track1_bunch_hom`, on by default via
  `global_com%mp_threading_is_safe`), so the threads are already busy inside each call
  and an outer parallel loop would nest.

Because each slice's arithmetic is identical regardless of which thread runs it (no
cross-slice reductions exist in the step loop), the result is **bit-identical across
thread counts**, and the harness checks that: the time-dependent single-segment
configuration reruns with `OMP_NUM_THREADS=8` against the 1-thread run, requiring the
diag file byte-equal and every dump dataset exactly equal. (Whole-file `cmp` of HDF5
would false-alarm on object-header timestamps. The comparison is per dataset.)

Measured scaling, full 6-FODO line, 32 slices x 2048 particles (Apple Silicon, debug
build):

| Threads | Wall time | Speedup | Efficiency |
|---|---|---|---|
| 1 | 129.6 s | 1.00 | the reference |
| 2 | 74.4 s | 1.74 | 87% |
| 4 | 46.4 s | 2.79 | 70% |
| 8 | 32.6 s | 3.97 | 50% |

It saturates near 4x at 8 threads. Two parallel regions are spawned per integration
step, and on this machine 8 threads includes efficiency cores. With 32 slices there is
also little schedule slack: 4 slices per thread at 8 threads. Production-size runs
(hundreds to thousands of slices) have more parallel work per serial byte, so this is
the floor of the scaling rather than its ceiling.

[](performance.md#perf-thread-scaling) measures the same curve at 96 slices in a
production build and gets 6.76x at 8 threads and 9.16x at 12, which holds that claim. It
also prices what the serial fraction is made of. The slippage rotation is 0.4% of the
walk, the per-slice diagnostics reduction runs slice-parallel, and every phase this file
can name totals 0.6%, so the rest is the cost of entering a parallel region.

### Head to head against Genesis4 at 12 workers

One command reproduces this measurement on any machine:

```
./lucifer/tests/run_perf_benchmark.sh
```

It detects the performance-core count (`--workers N` to override), sizes the window to
4 slices per worker so Genesis4 does not pad, finds the MPI launcher that matches the
Genesis4 binary's own linkage (an OpenMPI `mpirun` aborts an MPICH-linked genesis4 on
sight), times all four runs with one external clock, and refuses to print a table
unless the two codes' answers agree at the documented seam level first.

Measured on a 48-slice run (slice count divisible by the worker count, so Genesis4 does
not pad its window), full line, 2048 particles per slice, identical starting dumps,
production/release builds on both sides, M3 Max with 12 performance cores. The answers
agree at the documented seam level (4.0e-2) before any timing is quoted:

| | Serial | 12 workers | Parallel speedup |
|---|---|---|---|
| Genesis4 (MPI) | 152.7 s | 35.8 s (12 ranks) | 4.3x |
| `lucifer` (OpenMP) | 123.6 s | 19.4 s (12 threads) | 6.4x |

The tracker is 1.2x faster serial and 1.85x faster at 12 workers. The parallel gap is
the shared-memory bet showing up in a measurement: Genesis4's 12-rank run
spent 68 s of system time on the per-step MPI slippage ring exchange and diagnostics
gather, where this code's slippage is an index rotation in shared memory and costs
nothing to communicate.

(val-the-particlepath-cost-measured)=
## The particle-path cost: measured, two levers adopted

A real 42 m case (131 slices) at 2048 particles/slice ran at parity with Genesis4
(96 vs 94 s, 12 threads vs 12 MPI ranks, same machine back-to-back). At 8192/slice it
did not (164 vs 99 s). The per-element time stamps and an in-wiggler profile put the
gap in one place: the FEL step's particle path (the RK4 + gather), whose cost
quadrupled with the particles while the per-slice field FFTs (~68 s) stayed constant.
[](performance.md#perf-the-particles-against-the-grid) measures that crossover as a
table, and [](performance.md#perf-the-phase-profile) is the standing phase profile every
run's own footer now prints.

`tests/lucifer_advance_bench.f90` (built as `lucifer_advance_bench`) times the path serially
at fixed state -- measured (min-of-5, this machine): full fel_advance **152
ns/particle-step**, the RK4+ODE alone 87.5, the four sin/cos pairs the ODE evaluates
48. A clang -O3 twin of the same verbatim RK loop runs 68 ns, and its sincos loop
13.7 ns. The gap is codegen and libm rather than physics (macOS gfortran has no vectorized
libm and emits scalar calls).

Adopted, as one named value change (each pays nothing alone, together full
fel_advance drops to ~108-125 ns and the real case 164 -> 143 s):

- per-file `-O3` on `fel_track_mod` (CMakeLists `set_source_files_properties`), and
- the paired-sincos shim (`code/fel_sincos.c`): one libm call for the (sin, cos)
  pair at the ODE's three call sites.

The change is 1-ulp-level: gfortran's sin/cos intrinsics differ from libm's sincos
by one ulp on ~2e-6 of arguments (measured: 73 mismatches in a 44M-point sweep of
the theta domain. The C-side test of __sincos against libm sin/cos was clean, so
the divergence is gfortran's intrinsic lowering), and -O3 re-contracts floating
arithmetic. Per the tier contract, the full benchmark re-ran and every tier holds
its Genesis4 tolerance. Four recorded values moved at the ulp-amplified level
(tier1_unavg 6.934012e-02 -> 6.933979e-02, tier2_genesis 1.772017e-05 ->
1.771895e-05, tdsase 2.292890e-06 -> 2.292906e-06, weight_split 7.917158e-13 ->
3.508953e-13). The other seven are unchanged.

The interludes (the seam) now track SLICE-PARALLEL where they can: the per-slice
track1_bunch loop runs one slice per thread (coarse granularity instead of
per-particle inner threading), BIT-FOR-BIT at any thread count against the serial
baseline (checked at 1 and 8 threads on tier2 and td2). Two boundaries, both
deliberate: radiation fluctuations keep the serial loop (they draw from the one
shared RNG stream inside track1, whose order must stay fixed), and wake-carrying
elements were always whole-window concatenations (manual [](fel-physics.md#sec-seamwake)), so lattices whose
interludes all carry wakes see no change. One hard-won implementation note: the
scratch bunch is BLOCK-LOCAL inside the loop body, not an OMP `private` variable --
gfortran left a `private` copy of the allocatable-component derived type improperly
initialized (a deterministic 2x bunching shift, present even at one thread), while a
block-local declaration gets Fortran's own default initialization.

Also measured and worth knowing: the production (-O2) and debug (-O0) trees were
never bit-identical to each other (FP contraction). The bit-for-bit keystone is a
WITHIN-tree contract, and the harness runs the debug tree. And the remaining gap to
Genesis4 on the 8192 case (143 vs 99 s) now sits in the per-slice field FFTs and the
wake-path interludes, with the clang-twin numbers bounding what further particle-path
codegen work could still buy (~68 ns vs our ~125 per particle-step).

(val-fp32-lockstep)=
## The FP32 particle path, priced by lockstep

Apple GPUs have no FP64, so a single-precision form of the averaged FEL advance exists as a measured object before any kernel: an FP32 twin steps beside the FP64 path from a shared state at `fel_advance`'s own sequence point, and `global%fp32_check` records the per-quantity divergence (`fel_fp32_mod`, [](input-reference.md#param-global-fp32-check)). The twin carries the three reformulations FP32 forces, each documented in the module header with the quantum that forces it: energy as an offset from the reference (the FP32 quantum of absolute gamma is 1.35e-3 against a per-step change of 3.0e-4), phase as a residual with the base rotator evaluated in FP64 once per slice, and the detuning as a difference of like small quantities (the FP64 ODE's square-root argument is 1 minus ~1.3e-8, which is exactly 1 in FP32 and erases the phase equation).

The longitudinal residual needs a moving reference, and the instrument's own guard found that. The split is $z = z_{ref} + dz$ with the reference one FP64 number per slice; on the instrument's first run the reference was static, the beam's common $z$ drift accumulated into the residual, and at step 514 of 1068 the guard refused the run at 28 ulps against its floor of 32. A secular term in the residual is silent death by quantization, so the reference now follows the slice's FP64 mean each step, one host number per slice, and in freerun mode the persistent residuals rebase across the move with `fel_fp32_renorm`, the same operation migration needs (round trip measured at 1 ulp).

The field solve has its own twin: the deposit scattered into an FP32 source grid, FFTW's single-precision transform pair, and the propagator rounded from the FP64 kernel, applied to the twin's own field record each step. Two rows join the stream, the source grid and the post-solve field, both against the FP64 post-solve field's norm. In freerun the FP32 field carries across steps, fed by the twin's own deposit and read by its own gather, so the twin is a complete single-precision run beside the FP64 one, longitudinal and field state alike (the transverse maps stay FP64 and enter rounded). A freerun window is one slice by refusal: the twin keeps no slippage rotation of its own. That transform is FFTW's own single-precision interface, and the build detects it rather than assuming it: the conda toolchain carries `libfftw3f` and the off-site distribution's from-source FFTW is double precision only, so a build without the library refuses `fp32_check` by name instead of transforming in a precision other than the one it reports.

Measured levels, worst over the run per quantity (M3 Max, production build; debug agrees):

| case | pz (relative) | theta [rad] | phasor (of full charge) | source | field | guard [ulp] |
|---|---|---|---|---|---|---|
| steady state, 2048 particles, 57 m through saturation | 4.2e-7 | 4.7e-5 | 7.2e-6 | 7.2e-5 | 7.2e-5 | 1.9e3 |
| SASE, 96 slices x 2048, full line through saturation | 4.1e-7 | 5.7e-5 | 2.8e-6 | 4.1e-5 | 4.1e-5 | 8.2e2 |
| freerun (compounding), steady state | 1.8e0 | 3.3e1 | 2.1e-1 | 4.0e-2 | 1.3e0 | 2.2e3 |

The lockstep rows are the per-step cost of FP32: about 4e-7 on the energy chart, 5e-5 rad on the phase, under 1e-5 on the source term and under 1e-4 on the field step, flat across a run because the state refreshes each step. The source and field rows agree to three digits, so the deposit's FP32 phase noise dominates the step and the single-precision transform contributes beneath it, consistent with published transform-only figures of 2.2e-7 to 4.0e-7. The transverse rows measure representation only (~6e-8, FP32 epsilon), since the transverse maps stay FP64.

The end-to-end number, the one a device port is judged against:

| case | exit power (relative) | exit bunching (absolute) |
|---|---|---|
| freerun, one 3.99 m segment, exponential gain | 1.0e-3 | 2.9e-8 |
| freerun, 57 m steady state, deep saturation | 2.7e-1 | 3.6e-2 |

The gain-regime figure sits at the order of the published backends' end-to-end SASE comparison (6.7e-4). The deep-saturation figure is a different animal and reads as one: past saturation the power oscillates with the synchrotron rotation, an FP32 trajectory ends at a visibly different phase of that oscillation, and a 27% instantaneous power difference at one z is what phase decorrelation looks like there, with the bunching magnitude off by 3.6e-2 on 0.185. A device port validated per step by this instrument's lockstep rows is not exposed to that compounding; a device run trusted end to end through deep saturation is, in FP32, and now the number is on record rather than discovered later.

| check (harness section `fp32-lockstep`) | level |
|---|---|
| lockstep steady and TD levels within the recorded ceilings | pz 2e-6, theta 3e-4, phasor 3e-5, source 1e-4, field 1e-4 |
| residual guard healthy | >= 256 ulp (measured 5.7e3 to 6.8e3 on the check cases) |
| bucket renormalization round trip | <= 1.5 ulp (measured 1.0) |
| FP64 untouched: diag byte-identical, instrument on vs off | exact |
| the mutation moves the theta level and the source row | >= 10x and >= 5x (measured 12.5x, 29x) |
| freerun compounds past lockstep on the phasor | measured 14x |
| freerun end-to-end power divergence on the gain segment | <= 1e-2 (measured 6.2e-4) |
| instrument stream at 1 vs 8 threads | byte-identical |
| wake with the instrument on | refused by name |

The instrument is read-only on the physics by construction and by check: nothing outside it reads the FP32 state, and the instrumented run's FP64 outputs are byte-identical to the uninstrumented run's. Configurations the twin does not cover (harmonics, two polarizations, the coherent source, wakes, space charge, the unaveraged mode) are refused by name rather than half-measured.

(val-device)=
## The Metal backend, judged by the instrument

`global%device = "metal"` runs the averaged FEL step on an Apple Silicon GPU ([](input-reference.md#param-global-device)). The seam is one C interface (`device/lucifer_device.h`) behind `iso_c_binding`, one Objective-C++ file behind that, and a stub that refuses by name on every other platform; no Metal type or call appears anywhere else, so a second backend is a new implementation of the same interface plus one CMake branch. The design is this project's own prior work, the Genesis 1.3 v4 backends on branch `gpu/metal-engine` (commit 4919b01, unmerged upstream): the transform kernels, the one-command-buffer step and the buffer-sync rule transcribe `MetalEngine.mm` from that branch, while the physics kernels transcribe Lucifer's own `fel_fp32_mod` twin and the Bmad transverse map, so the arithmetic the device runs is the arithmetic the lockstep instrument already priced.

One deliberate divergence from the reference backend, stated where the citation is: the longitudinal state is not an absolute FP32 theta but a 64-bit fixed-point phase, in ticks of $2\pi/2^{32}$ off a static FP64 per-slice reference. The low 32 bits are the phase modulo one radiation period at a uniform 1.5e-9 rad, the high bits count whole periods, and a bucket crossing is exact integer arithmetic -- the migration operation a later landing needs. The uniform quantum is what lets the reference stay static for a whole element where the FP32 residual of [](#val-fp32-lockstep) needed a moving one. Bucket wraps are modular arithmetic and are asserted exactly, on the device's own arithmetic at every setup: a bucket shift and its return must be bit-exact and the extracted phase must never see the shift, a statement no tolerance is allowed to soften.

The backend is judged by the lockstep machinery with the device in the twin's role: `device = "metal"` with `fp32_check = "lockstep"` leaves the FP64 run untouched (diag byte-identical, on against off) and the device advances every step from the shared pre-step state, whole -- its own transverse maps, push, deposit and field solve -- with the rows in the same `.fp32.txt` stream. The ceilings are three times the recorded CPU-twin levels above, absorbing that Metal's precise-math sincos is not libm and its fused ops round differently; the transverse rows get their own recorded levels, since the CPU twin refreshed those from FP64 and the device runs its own FP32 maps. A jump against the refreshed reference is a kernel bug by definition. Measured levels, worst over the run (M3 Max, both builds agree; the check cases of `check_device.py`, grid 64):

| case | x..py (relative) | pz (relative) | theta [rad] | phasor | source | field | guard [tick] |
|---|---|---|---|---|---|---|---|
| steady state, 1024 particles | 2.7e-7 | 1.4e-7 | 6.6e-6 | 7.3e-8 | 1.5e-5 | 1.5e-5 | 8.8e6 |
| SASE window, 8 slices x 1024, slippage live | 2.8e-7 | 1.4e-7 | 6.8e-6 | 1.0e-7 | 6.4e-6 | 6.4e-6 | 6.8e6 |

The guard column is the median per-step phase increment in ticks of the fixed-point quantum; its floor plays the same silent-z role as the residual guard above. The mutation hook reaches the kernels themselves: `fp32_mutate` perturbs the detuning resonance inside the shader by one part in $2^{12}$ and coarsens the uploaded phase, and the theta and source levels move by 710x and 52x, so the ceilings are demonstrated to catch a wrong kernel constant, not just a wrong upload.

End to end the two roles corroborate each other. The freerun twin (state resident and compounding, single-slice window) lands at 5.1e-4 on exit power; the production run (`fp32_check` off, the device is the run, beam and field resident through each element with readbacks only at the stats comb's positions) lands at 5.1e-4 against the CPU run's diag rows on the same deck -- the freerun-measured order, and the same number twice by different routes. Wall times against the CPU are in [](performance.md#perf-device), with the dispatch floor stated.

| check (harness section `device`) | level |
|---|---|
| device lockstep steady and TD levels within the recorded ceilings | x..py 1e-6, pz 1.3e-6, theta 1.5e-4, phasor 2.2e-5, source 2.2e-4, field 2.2e-4 |
| phase-increment guard healthy | >= 1e6 ticks (measured 6.8e6 to 8.8e6) |
| bucket wraps exact on device arithmetic | asserted bit-exact, never a tolerance |
| FP64 untouched: diag byte-identical, device twin on vs off, steady and TD | exact |
| the perturbed kernel constant moves the recorded levels | theta >= 10x, source >= 5x (measured 710x, 52x) |
| device freerun compounds past lockstep on the phasor | measured 14x |
| freerun and production exit power vs CPU, steady and TD | <= 1.6e-3 (measured 5.1e-4, 5.1e-4, 1.5e-5) |
| wakes, migration, radiation, harmonics, the unaveraged mode, an unknown backend, an unsupported grid | each refused by name |

One property is inherited from the reference backends and documented rather than fought: the deposit accumulates with device atomic adds whose ordering is not fixed, so two runs of the same step differ in the source's last bit or two. No device output is therefore asserted byte-identical against another device run; the read-only proofs compare FP64 outputs only, and the ceilings absorb the flutter. Everything the kernels do not cover is refused by name at setup or first use, and a build without the backend refuses the knob itself with the stub's own message. Which build that is comes from detection rather than from the platform: the backend needs macOS and a Clang-family Objective-C++ compiler, since it is ARC-managed Objective-C++ against the Metal framework, and a macOS build whose Objective-C++ compiler is a GNU one takes the stub like any other. The configure output names the backend it built, and a capability-free build was measured to produce byte-identical physics to a full one.

(val-the-coarsestep-measurement)=
## The coarse-step measurement

(Summarized in manual [](fel-physics.md#sec-numerics).)

SIMPLEX's reference case integrates twelve undulator periods per step with the slice
spacing matched so slippage is one slice per step -- an order of magnitude fewer steps
than one-step-per-period. Whether that economy transfers to this integrator is the
question, measured by `run_delz_sweep.sh`: Genesis4 generates one
time-dependent initial state (32 slices of spacing `12*lambda0`, shot noise on), and the
tracker runs the full line from that same dump at `ds_step` of 1, 2, 3, 6 and 12 periods
(two-line wrapper lattices overriding `ds_step`), so
every run shares one shot-noise realization and the differences are pure integration
error. That sharing is what the measurement rests on, and it holds whichever realization
the reference produces, so the numbers below stand as measured against the reference of
their day (a Genesis4 predating the 4.6.15 seeding change). A rerun on the current
reference draws a different realization and would move them in their own right, which is
a property of the comparison rather than of the integrator. Total power at the twelve
undulator-segment exits, against the one-period run:

| `ds_step` | max over exits | at saturation |
|---|---|---|
| 2 periods | 1.3e-1 | 1.5e-2 |
| 3 periods | 2.5e-1 | 2.8e-2 |
| 6 periods | 5.1e-1 | 2.1e-2 |
| 12 periods | 7.0e-1 | 2.6e-1 |

The max-over-exits error lives in the exponential-gain region, where a step-size error is
a quasi-systematic gain-length shift. Fitted there, convergence is roughly first order
(pairwise p of 1.7, 1.1, 0.5). The saturation error is oscillatory (post-saturation power
oscillates, so a small phase shift moves the sampled value) and fits no clean order.

The answer to 8.3: saturation power holds to ~3% up to six periods per step, and the
twelve-period matched configuration misses it by 26%. SIMPLEX's step-size economy is
tied to its semianalytic field advance and does not transfer to this Genesis4-style
integrator as-is. `ds_step` of two to three periods is the operating point here. Six periods
is defensible when only saturation power matters.

(val-the-saturation-demo)=
## The saturation demo

The practical end-to-end case: Genesis4's own Benchmark1-SASE configuration run to saturation over the full 57 m Aramis line, dark start, 96 slices, tracked three ways from identical dumps against one external clock. Its inputs, its run script and its measured result live with the demo itself, in [`../examples/saturation_demo/README.md`](../examples/saturation_demo/README.md). The headline: the averaged mode matches Genesis4 at 4.9e-4 through saturation, and the unaveraged mode rides +0.6 ln on the shot-noise radiation channel it physically resolves. On that configuration and machine (M3 Max, 12 performance cores, production builds on both sides) it also ran 1.26x faster at equal cores. That is one case on one machine, not a general claim about either code.

Properties of Genesis4 that this transcription had to establish by reading its source are recorded in [](genesis4.md), which is the one page that tracks Genesis4 as it changes.
