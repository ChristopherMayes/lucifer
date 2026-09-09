---
title: Performance
short_title: Performance
---

Where a run spends its time, measured. Every number here is machine-local and carries the machine and the build that produced it, which is why this page is hand written and stays out of `generated/`.

The instrument is `code/fel_timer_mod.f90`. It accumulates wall clock per phase at region boundaries and the footer of every run prints the table. The phases partition the walk, so the fractions sum to it and the remainder has its own row named `unaccounted`. Measured below, that remainder is 0.02% of the walk, so the partition is the walk.

Timers cannot see inside a parallel region. The deposit and its FFT interleave per slice, and so do the field gather and the RK4, so a clock call between them would measure the loop it perturbs. That split comes from a sampling profiler and appears in its own section.

The tables below are also the calibration of the cost estimate every run prints in its header. `code/fel_cost_mod.f90` holds the fitted rates, one for the particle push and two for the field solve, and `tests/scripts/fit_cost_rates.py` fits them to the logs the runs below wrote. Over the thirteen runs of this page the worst residual on the walk is 4.0 percent. The estimate is for this machine and says so on the line, and the run's footer prints the measured walk beside it. Re-measuring this page means re-running that script and pasting its five constants back.

(perf-reproducing)=
## Reproducing this

```
./lucifer/tests/run_perf_benchmark.sh --phases
```

No genesis4 and no MPI: the phases mode compares the code against itself, so it starts from the tracker's own shot-noise quiet start rather than an imported dump. It runs the profile at two slice counts, sweeps the thread count over the smaller one, and prints both tables. `--npart` sets the particles per slice and `--big-slices` the larger case.

The sampling split, on macOS:

```
OMP_NUM_THREADS=12 production/bin/lucifer <deck>.in &
sample $(pgrep -n lucifer) 15 1 -file prof.txt -mayDie
```

Read the `Sort by top of stack` section, which is exclusive cost per symbol. The call graph is per thread and fragments the same information across dozens of stacks.

The vectorization audit:

```
./lucifer/tests/scripts/vec_audit.sh
```

It takes the compile command from the production build's own makefile, so the report always describes the binary that ships. Build production first.

(perf-the-phase-profile)=
## The phase profile

Apple M3 Max, 12 performance cores, production build, full 6-FODO Benchmark1-SASE line, 2048 particles per slice, `ngrid` 255, 12 threads, the source filter at its default. The `was` columns are the same runs measured before that default changed, and they are given wherever the move exceeds the run-to-run spread below.

| phase | 96 slices | was | 504 slices | was |
|---|---|---|---|---|
| field solve (deposit and FFT) | 21.687 s, 71.5% | 15.435 s, 64.5% | 118.060 s, 71.7% | 82.072 s, 64.4% |
| particle push (transverse maps, RK4) | 4.735 s, 15.6% | 4.485 s, 18.8% | 24.509 s, 14.9% | 23.608 s, 18.5% |
| `= FEL step`, the sum of the phases inside it | 26.468 s, 87.3% | 19.966 s, 83.5% | 142.645 s, 86.6% | 105.731 s, 83.0% |
| stats and diag | 3.140 s, 10.4% | 3.201 s, 13.4% | 18.559 s, 11.3% | 18.223 s, 14.3% |
| field drift through the breaks | 0.490 s, 1.6% | 0.510 s, 2.1% | 2.597 s, 1.6% | 2.556 s, 2.0% |
| seam interlude | 0.096 s, 0.3% | 0.103 s, 0.4% | 0.578 s, 0.4% | 0.574 s, 0.5% |
| slippage | 0.104 s, 0.3% | 0.102 s, 0.4% | 0.106 s, 0.1% | 0.109 s, 0.1% |
| undulator prep (plans, kernel cache) | 0.045 s, 0.1% | 0.046 s, 0.2% | 0.075 s, 0.0% | 0.051 s, 0.0% |
| element end | 0.032 s, 0.1% | 0.032 s, 0.1% | 0.159 s, 0.1% | 0.170 s, 0.1% |
| unaccounted | 0.004 s, 0.0% | 0.003 s, 0.0% | 0.007 s, 0.0% | 0.006 s, 0.0% |
| walk | 30.333 s | 23.917 s | 164.650 s | 127.370 s |

One phase moved and it moved the walk with it. The field solve is up 40 percent at 96 slices and 44 percent at 504, and everything else is inside the spread. The source filter is why: it adds a forward transform of the source per field per step, where an unfiltered solve adds the source in real space and transforms twice rather than three times ([](fel-physics.md#sec-source-filter)). Measured directly on the 96-slice configuration, the same deck with the filter off and on:

| | filter off | filter on | ratio |
|---|---|---|---|
| field solve | 15.011 s | 21.339 s | 1.42 |
| particle push | 4.484 s | 4.441 s | 0.99 |
| walk | 23.392 s | 29.646 s | 1.27 |

Two runs each, mean. The filter costs 27 percent of this walk and 42 percent of its field solve, and it leaves the push alone. The 8 to 14 percent [](fel-physics.md#sec-convergence) quotes is the cost at the derived grid on a shorter line, where the transform is a smaller share of a run; this line at `ngrid` 255 is where the transform dominates, so it is where the filter costs most. The filter-off walk of 23.392 s reproduces the 23.917 s recorded before the default changed, so nothing else in the code moved between the two measurements.

The `= FEL step` row is derived rather than measured, which its label says: it is the sum of the undulator prep, the space-charge profile, the particle push and the field solve, so it is the share of the walk the undulator segments own. It is left out of the leaf sum, and the other rows still add to the walk minus the remainder.

The shares are the result and the third decimal is noise. Two 96-slice runs back to back give walks of 29.625 s and 29.667 s, and a third taken in a different batch gives 30.333 s, so the run-to-run spread is about 2% on a walk and a few tenths of a point on a share.

The shares hold across a factor of 5.25 in slice count. The run row against the walk row is everything outside the walk: the parse, the beam and field build, the slippage schedule, the final dumps and the stats file. That is 0.2 s at 96 slices and 0.8 s at 504.

Slippage is 0.4% of the walk and 0.1% at the larger size. [](validation.md#val-parallelism-openmp-over-slices) named the slippage rotation as part of the serial cost that capped scaling, which the measurement does not support. The field solve and the particle push are where the time is, and both are slice-parallel.

(perf-the-particles-against-the-grid)=
## The particles against the grid

The field solve and the particle push scale differently, so their ratio is a property of the configuration and not of the code. Same machine and line, 96 slices, 12 threads.

Two runs at each count, mean.

| particles per slice | particle push | field solve | walk |
|---|---|---|---|
| 512 | 1.178 s, 4.6% | 21.182 s, 83.2% | 25.446 s |
| 2048 | 4.441 s, 15.0% | 21.339 s, 72.0% | 29.646 s |
| 8192 | 16.722 s, 36.2% | 22.934 s, 49.6% | 46.224 s |

The push is linear in the particles (1.178, 4.441, 16.722 for 512, 2048, 8192) and the field solve is nearly flat (21.182, 21.339, 22.934). Three FFTs of $255^2$ points per slice per step do not care how many particles deposited into the grid. The crossover has moved past 8192 particles per slice at `ngrid` 255, since the filter raised the flat part and left the push where it was.

That flatness also splits the field solve without a profiler. Fitting the 512 and 8192 points to $C + d \cdot n$ gives $C = 21.06$ s and $d = 2.28 \times 10^{-4}$ s per particle. The slope is the deposit and it is unchanged, which it should be: the filter multiplies a transformed source and does not touch how the source is built. At 2048 particles per slice the fit predicts 21.06 s of transform and propagator plus 0.47 s of deposit, so 21.53 s against the 21.339 s measured, a 0.9% residual on a two-point fit. The transform and propagator are 98.7% of the solve by that fit against 95.6% by the sampling profile below, where the two agreed to half a point before the filter. The fit's deposit term is now 2% of a larger solve, so the two-point $C$ carries the solve's own 1.5% spread almost whole.

[](validation.md#val-the-particlepath-cost-measured) measured a real 131-slice case where the particle path dominated at 8192 particles per slice while the field FFTs stayed constant. That is this table, and the two agree.

(perf-thread-scaling)=
## Thread scaling

96 slices, 2048 particles per slice, full line, same machine and build. The implied serial fraction inverts Amdahl's law, $f = (n/S - 1)/(n - 1)$, and bounds this configuration rather than the code.

| threads | walk | speedup | efficiency | implied serial |
|---|---|---|---|---|
| 1 | 290.8 s | 1.00 | the reference | |
| 2 | 149.7 s | 1.94 | 97% | 3.0% |
| 4 | 76.5 s | 3.80 | 95% | 1.7% |
| 8 | 41.6 s | 7.00 | 87% | 2.0% |
| 12 | 30.3 s | 9.59 | 80% | 2.3% |

The scaling improved with the filter on, from 9.16x to 9.59x at twelve threads and from 76% efficiency to 80%. The filter's work is a transform per slice, which is the parallel part, so adding it lengthens the parallel region and leaves the serial part alone. The estimate takes 2.1% as the serial fraction, fitted over the whole sweep.

[](validation.md#val-parallelism-openmp-over-slices) records 3.97x at 8 threads on 32 slices in a debug build and says that is the floor of the scaling rather than its ceiling. At 96 slices in a production build the same 8 threads give 6.76x, so the claim holds.

The nameable serial phases account for 0.6% of the walk, well under the 2.8% the speedup implies. The rest is not serial code. Two parallel regions are spawned per integration step (1068 steps, so about 4300 region entries), and the sampling profile below puts 10.0% of its samples in the thread wait, which at 76% efficiency is the barrier and not a lock.

A steady-state run gets nothing from threads. The parallelism is over slices, and a steady-state window has one: the unaveraged example measures 21.94 s at 1 thread and 21.76 s at 12.

(perf-the-sampling-split)=
## The sampling split, averaged mode

96 slices, 2048 particles per slice, 12 threads, the filter at its default, 15 s of samples at 1 ms, 120732 samples. Exclusive cost, categorized by symbol. The `was` column is the same measurement before the default changed.

| | share | was |
|---|---|---|
| FFT transform (libfftw3, `wavefront_fft2`) | 68.8% | 61.3% |
| libm sin and cos | 6.5% | 8.5% |
| diagnostics | 5.5% | 4.6% |
| thread wait | 5.1% | 10.0% |
| deposit and propagator multiply | 4.9% | 5.7% |
| transverse maps (`quad_mat2_calc`, and its `cexp`) | 4.9% | 5.9% |
| longitudinal RK4 and ODE | 3.7% | 3.6% |
| slippage and drift | 0.1% | 0.1% |
| everything else | 0.6% | 0.4% |

The transform took the filter's third pass and rose 7.5 points, which is the whole of what moved. The thread wait halved, and the thread scan above agrees: a longer parallel region against the same serial part is better efficiency, not worse.

An earlier estimate put the field solve near 10% of a run, and at these parameters it is 71.5% by the phase timer and 68.8% in the transform alone by samples. That estimate came from a configuration with many more particles per slice, where the table above shows the balance reversing.

`quad_mat2_calc` is Bmad's own quadrupole map, and its cost is a complex exponential per particle: `cexp` and the `exp` beneath it are the whole 5.9%. The libm sin and cos are the deposit's phase factor and the ODE's, through `fel_sincos`.

(perf-the-sampling-split-unaveraged)=
## The sampling split, unaveraged mode

One slice, 16384 particles, `ngrid` 129, serial (the mode's parallelism is over slices and there is one). The phase table gives the unaveraged step 99.7% of the walk, so this is the mode. Both columns are 20 s of samples at 1 ms, and the totals differ because the hoisted run finishes sooner inside that window.

| | before the hoist | after |
|---|---|---|
| libm transcendentals | 52.1% | 15.3% |
| substep RK4 (`unavg_push`, `unavg_ode`) | 17.6% | 35.3% |
| undulator field, this code's own arithmetic | 8.9% | 27.7% |
| FFT transform and `fel_field_diffract` | 13.0% | 17.8% |
| ramp envelope (`fel_unavg_envelope`) | 4.8% | 0.0% |
| deposit, this code's own arithmetic | 1.2% | 1.9% |
| everything else | 2.4% | 1.8% |
| total samples | 13302 | 7964 |

Half of this mode was libm, and none of it went through `fel_sincos`. `fel_unavg_bfield` called `cos(und%ku * s)` and `sin(und%ku * s)` as two separate intrinsics and `fel_unavg_envelope` called two more, all four per particle per RK stage, when the argument depends only on the substep position and is therefore the same for every particle in the loop. The four s-dependent factors now arrive as arguments, evaluated once per stage, and `fel_unavg_bfield` no longer takes `s` at all: the s-independence is structural rather than a comment. FINDINGS 7.37 records the defect.

The share row hides how much moved, because the denominator moved too. In absolute samples libm fell from 6936 to 1221, which is 82% of the transcendental work gone, and the step's wall clock fell by the amounts in the next table.

Read the before column's label carefully, since the P1 table got it wrong. It said "libm sin and cos" where it measured every libm transcendental the run entered, `exp` and `cexp` and the shared reduction helper included, and only 13.9% of the run was `sin` and `cos` proper. The label was wrong and the number was right for what it summed. The hoist recovers most of it either way, because the envelope's own trig was per-particle for the same reason.

(perf-the-hoist-measured)=
## The hoist, measured

Same machine and build, one thread, median of five runs at 2048 particles and three at 16384. The change is bit-identical by construction and checked as such: the diag and ledger files of the unaveraged example are byte-identical before and after at 1 and 8 threads.

| configuration | before | after | change |
|---|---|---|---|
| 2048 particles per slice | 5.293 s | 4.507 s | -14.8% |
| 16384 particles per slice | 21.842 s | 15.806 s | -27.6% |

The gain grows with the particle count because what it removes is per-particle and what remains at fixed count is the per-substep FFT.

The next lever after the hoist was the loop order, and it landed. The substep RK4 and the undulator field arithmetic were 63% of the mode's samples, four RK stages inside each particle of a `do ip` loop. The push now runs stage-major: per stage, one loop over the slice's particles into per-particle stage arrays, then one loop for the RK4 combination. Per particle the operations, values and order are unchanged, so the results are byte-identical, checked the same way as the hoist.

| configuration | before | after | change |
|---|---|---|---|
| 2048 particles per slice | 4.577 s | 4.217 s | -7.9% |
| 16384 particles per slice | 16.073 s | 13.189 s | -17.9% |

The win is locality and call overhead, and deliberately not vectorization. The stage loops still refuse to vectorize because each iteration calls `unavg_ode`, and inlining the field into the stage loop was tried and refused on the identity check: the compiler fuses the inlined arithmetic differently than it fused the called form, the results differ in the last bits, and a moved bit is a moved digit somewhere downstream. So the call stays, and the vectorizer's verdict stands recorded rather than fought. What vectorization of this loop now needs is a reduction-free arithmetic core whose bits are the contract rather than an accident of fusion, which is the FP32 lockstep work's territory, where the reference is measured divergence rather than identity.

After the interchange the mode's samples split 37.5% substep RK4, 21.4% FFT, 18.8% undulator field arithmetic, 18.7% the deposit's own transcendentals, 2.2% deposit weights.

(perf-the-sincos-ceiling)=
## The sincos ceiling

An inline (sin, cos) pair with no libm call, Cody-Waite reduction and the fdlibm kernel polynomials, swapped into `fel_sincos.c` for the measurement and reverted. A first probe build recorded the argument range of that call site: $[-159.8, 19.0]$ on the 96-slice averaged configuration and $[-127.4, 17.1]$ on the unaveraged example. Both sit well inside the range a two-term Cody-Waite reduction handles, and the probe run reproduces the reference exit power, pulse energy and bunching to every printed digit. The physics is the same, so the timing difference is the call.

| case | reference | probe | change |
|---|---|---|---|
| 96 slices, 2048 particles, 12 threads | 23.917 s | 23.158 s | -3.2% |
| 96 slices, 8192 particles, 12 threads | 41.375 s | 38.197 s | -7.7% |
| unaveraged, 16384 particles, serial | 21.937 s | 21.848 s | -0.4% |

The averaged rows bound the libm call from below at 3.2% of a run at 2048 particles per slice, against the 8.5% the sampler attributes to sin and cos: the polynomial is cheaper than the call and not free. The unaveraged row was the finding. That mode never reaches this file, and its own transcendentals were the undulator field's and the ramp envelope's intrinsics, which the hoist above removed instead.

Neither probe is a candidate implementation, and the second one is refused on its own audit rather than on judgement. `tests/scripts/sincos_audit.c` sweeps it against libm:

| range | worst error | arguments differing |
|---|---|---|
| the measured range, $[-159.8, 19.0]$ | 6 ulp | 34% |
| $\lvert \theta \rvert < 10^4$ | 6 ulp | 34% |
| $\lvert \theta \rvert < 10^5$ | 10 ulp | 34% |
| $\lvert \theta \rvert < 10^6$ | 24 ulp | 34% |
| $\lvert \theta \rvert < 10^8$ | 2681 ulp | 34% |

The shim it would replace was admitted at one ulp on $2 \times 10^{-6}$ of arguments, 73 mismatches in a 44M-point sweep ([](validation.md#val-the-particlepath-cost-measured)). This candidate is six ulp on a third of them, which is five orders of magnitude more arguments and six times the error, so it cannot be recorded as a 1-ulp change. The error is the two-term reduction rather than the kernels, and a reduction accurate enough to reach one ulp costs arithmetic that eats a 3.2% gain. The range is also not a constant: $\theta$ carries the common phase accumulating along the line, so a longer line moves toward the degradation with nothing to detect it.

The version of this that could pay is a vectorized pair inside a loop that vectorizes, where one reduction amortizes over a whole vector. The averaged path has no such loop, for the reasons in the audit section below. FINDINGS 7.38 records the refusal.

(perf-the-vectorization-audit)=
## The vectorization audit

gfortran 13.3.0 at the production flags, which are `-O2 -ftree-vectorize -march=armv8.3-a` with per-file `-O3` on `fel_track_mod`.

| file | vectorized | not vectorized |
|---|---|---|
| `fel_track_mod.f90` | 26 | 158 |
| `fel_unaveraged_mod.f90` | 12 | 33 |
| `fel_collective_mod.f90` | 5 | 49 |

Not one of the hot per-particle loops is among the vectorized. Every one of them is refused for the same reason, `control flow in loop`:

| loop | routine | what blocks it |
|---|---|---|
| `fel_track_mod.f90:1141` | `fel_transverse_track` | `fel_apply_focus`, called per particle, branches inside |
| `fel_track_mod.f90:1213` | `fel_transverse_track_bmad` | the same, plus `quad_mat2_calc` |
| `fel_track_mod.f90:1434` | `fel_advance` | `faw`, `fel_sincos` and `fel_runge_kutta`, called per particle |
| `fel_track_mod.f90:1494` | `fel_advance`, the multi-harmonic twin | the same |
| `fel_track_mod.f90:2099` | `fel_field_step`, the deposit | `if (.not. on_grid) cycle`, and the report also names `statement clobbers memory: fel_grid_weights` |

The aggregate reasons, most common first, are `vectorization is not profitable` (56 in `fel_track_mod`), `control flow in loop` (30) and `complicated access pattern` (19). The first is mostly short fixed-length loops over the 6-vector, which the compiler is right about. The second covers the loops in the table.

One pattern accounts for all five. A per-particle loop calls a routine that the compiler does not inline, the routine carries a branch, and the branch ends the loop's chance of vectorizing. The deposit's `on_grid` guard is the clearest case: it is a mask written as a `cycle`.

(perf-what-the-numbers-say-to-do)=
## What the numbers say to do

Ordered by measured share, for the averaged mode at 2048 particles per slice and `ngrid` 255.

The FFT is 61% and it is in FFTW, not in this code. What this code controls is how many transforms it asks for: two per slice per field per step, at the grid the deck sets. `ngrid` is a deck parameter and $255^2$ is 65025 points against 2048 particles, so a convergence study on `ngrid` is worth more here than any change to the deposit.

The libm transcendentals are 8.5% of the averaged mode. They were half the unaveraged mode and the hoist took 82% of that, so what remains there is the RK4 itself. Bmad's `cexp` in the transverse map is a further 5.9% of the averaged mode, and a quadrupole map has a real closed form.

The per-particle loops do not vectorize, and the blocker is uniform: an un-inlined call with a branch in it. Three of the five in the table above resist for reasons worth stating rather than fixing.

The deposit's `on_grid` guard cannot become a mask. `fel_grid_weights` leaves `ix` and `iy` unassigned when it reports off-grid, so a masked body would index `crsource` with an undefined subscript. Making that safe means clamping the indices, which is no longer the same program.

`fel_apply_focus` is called only by `fel_transverse_track`, which runs only under `transport_model = "genesis"`. That default is `bmad`, no example sets it, and only the comparison tiers select it, so restructuring that loop speeds up a validation-internal path. `faw`, on the production path, carries one branch on `und%sin_t` which is a property of the element and already loop-invariant.

`fel_advance`'s loops call `fel_sincos` and `fel_runge_kutta` per particle, so `!$omp simd` cannot apply to them until the sincos itself is vectorizable.

Thread scaling at production slice counts is 9.16x on 12 cores and the nameable serial cost is 0.6% of the walk, so there is little left to win by removing serial work.

The single-precision question is measured rather than argued. `global%fp32_check` steps an FP32 twin of the averaged advance beside the FP64 path and records the divergence per quantity, with the reformulations FP32 forces (offset energy, residual phase, difference-form detuning) and a runtime guard on the one failure that is silent, a residual too coarse for the per-step increment. The measured levels and the guard's own first catch are in [](validation.md#val-fp32-lockstep), and the field solve has its own twin there: FP32 deposit, single-precision transform pair, rounded propagator, with an end-to-end freerun number for the complete single-precision run (1.0e-3 on exit power over a gain segment, and 2.7e-1 in deep saturation where phase decorrelation dominates). An FP32 six-vector is 24 bytes per particle against 48, which is the bandwidth argument for a device path. Whether an FP32 CPU mode ever ships is a decision against those recorded levels, not taken here.

(perf-estimate)=
## The estimate against the clock

Every run prints an estimated walk in its header and the measured walk beside it in its footer. The rates are fitted to the tables above and are this machine's. Measured here on the eleven comparison tiers and four examples, each run on its own rather than beside the other keystone jobs, since a machine running five jobs at once measures the load and not the deck.

| run | filter | measured | estimate | estimate/measured |
|---|---|---|---|---|
| tier1 | off | 0.297 s | 0.394 s | 1.33 |
| tier1_unavg | off | 12.6 s | 23.6 s | 1.87 |
| tier2_bmad | off | 3.8 s | 4.7 s | 1.24 |
| tier2_genesis | off | 3.7 s | 4.7 s | 1.27 |
| td1 | off | 0.671 s | 0.845 s | 1.26 |
| td2_bmad | off | 8.6 s | 10.1 s | 1.17 |
| td2_genesis | off | 8.8 s | 10.1 s | 1.15 |
| tdsase | off | 8.9 s | 10.1 s | 1.13 |
| tdsc | off | 0.738 s | 0.845 s | 1.15 |
| tdwk | off | 0.691 s | 0.845 s | 1.22 |
| weight_split | off | 0.453 s | 0.576 s | 1.27 |
| `sase` | on | 31.5 s | 30.4 s | 0.97 |
| `steady_state` | on | 5.0 s | 4.7 s | 0.94 |
| `taper` | on | 5.1 s | 4.7 s | 0.92 |
| `flash1` | on | 22.6 s | 18.5 s | 0.82 |

Every row is inside a factor of two, which is what the estimate is held to. The filter column is what separates the two blocks. The examples take the default, and three of the four now read between 0.92 and 0.97 where all four read 0.80 against the previous rates. Every tier turns the filter off, since each compares against a code that carries none ([](validation.md)), so their field solve is the cheaper two-transform one while the rate pricing it is the filtered default's. That is the 1.13 to 1.33 the tier block sits at, and it is a property of the configuration rather than of the machine.

`tier1_unavg` is the highest at 1.87 and is high for a second reason: the unaveraged mode resolves the quiver, so it integrates and deposits on twenty substeps per period, and the count carries those substeps while the rates it multiplies are the averaged path's. `flash1` is the lowest at 0.82, on 351 slices at `ngrid` 129, which is the furthest configuration here from the one the rates were fitted on.

(perf-frames)=
## The frame series, measured

What a frame costs, on the same machine, so that the layout questions a series raises are decided against numbers rather than against expectations. The deck is the `sase` example, 96 slices of 2048 macroparticles at `ngrid` 255, one frame every 20 m.

| | size | per unit |
|---|---|---|
| field frame, 96 x 255 x 255 complex128 | 95.3 MB | 1.02 kB per grid point |
| beam frame, 196608 macroparticles | 12.8 MB | 68.2 B per macroparticle |

The beam frame carries no overhead to trim. Its 68 bytes are six coordinates, the time the reference phase folds into, the time offset that places the particle's slice in the bunch, the reference momentum and the label, and nothing else: the uniform weight and the zero longitudinal position are constant records that cost nothing per particle. The field is seven times the beam and is where a series' bytes are.

Compression does not pay on this data. Chunked one slice at a time with the shuffle filter:

| | size | ratio | read |
|---|---|---|---|
| stored, uncompressed | 95.3 MB | | 0.009 s |
| gzip level 1 | 73.2 MB | 1.30x | 0.135 s |
| gzip level 4 | 72.6 MB | 1.31x | 0.139 s |
| float32 | 47.6 MB | 2.00x | |

A saturated FEL field is not sparse. The wide-angle emission of the point beamlets fills the grid, and a float64 mantissa of a noise-like field is incompressible, so gzip buys 1.3 times for sixteen times the read. Storing the field as float32 halves it exactly and costs 5.9e-8 of relative accuracy, which is far below anything a picture or a moment needs and far above single-precision tracking's own floor. Neither is taken here: the measurement is what a later decision is made against.

Writing frames dominates a run at a movie-grade comb, and the device readback does not. One undulator segment, 15 slices of 4096 at `ngrid` 256 on the Metal backend, a comb of 0.1 m and 31 records:

| | frames off | frames on |
|---|---|---|
| walk | 0.134 s | 0.582 s |
| dumps | 0.000 s | 0.402 s, 69.1% |
| stats and diag | 0.027 s | 0.025 s |
| device step | 0.001 s | 0.002 s |

The readback is about 0.10 s either way, since it happens at every comb row whether or not a frame is written, and it is the unaccounted remainder in both columns. With frames on it is a sixth of the walk and the file writing is two thirds. So a series that costs too much is not answered by moving the per-slice reductions onto the device, which would save readback the run is not spending its time in. It is answered by writing fewer slices, which `global%dump_slice_first` and `dump_slice_last` do, or fewer records, which the comb does.

(perf-device)=
## The device, measured against the CPU

The Metal backend ([](validation.md#val-device)) run against the CPU path on the same decks, one 9 m Aramis segment of 89 steps at `ngrid` 256, `comb_ds_save = -1` so both sides pay stats at the element end only. Walk seconds, best of three, M3 Max (12 performance cores for the CPU, production builds). The device-busy column is the backend's own command-buffer timestamps; the difference from its wall clock is the host side of each step, dominated by the dispatch floor.

| case | CPU, 12 threads | device wall | device busy | ratio |
|---|---|---|---|---|
| steady state, 1 x 8192 particles | 0.360 s | 0.049 s | 0.015 s | 7.3x |
| SASE window, 24 x 4096, slippage live | 0.873 s | 0.107 s | 0.048 s | 8.2x |
| SASE window, 96 x 8192 | 3.92 s | 0.316 s | 0.196 s | 12.4x |
| dispatch floor: 1 x 1024, ngrid 64 | 0.032 s | 0.037 s | 0.009 s | 0.9x |

Where a step's device seconds go, per pass. `global%device_timing` encodes each pass into
its own compute encoder carrying a timestamp attachment, which is what stage-boundary
counter sampling can measure on an Apple GPU: dispatch-boundary sampling would time each
dispatch inside one encoder and the M3 Max reports it absent. The full Aramis line, 12
segments and 1068 steps, one slice, `ngrid` 256, the filter on and `comb_ds_save = -1`, so
the readbacks fall at element ends and each element's steps batch into one command buffer.
Best complete run of five, never a per-pass minimum, with the spread over the five beside
it. The transverse map runs twice a step, and the filter and the solve are four dispatches
each, so those lines sum more encoders than the others.

| particles a slice | walk | busy | transverse | push | zero | deposit | filter | solve | the six | spread |
|---|---|---|---|---|---|---|---|---|---|---|
| 1024 | 0.953 s | 0.311 s | 0.0278 s | 0.0273 s | 0.0143 s | 0.0201 s | 0.0831 s | 0.0871 s | 0.2597 s | 12.5% |
| 8192 | 0.875 s | 0.256 s | 0.0197 s | 0.0184 s | 0.0116 s | 0.0142 s | 0.0577 s | 0.0605 s | 0.1821 s | 27.7% |
| 131072 | 1.800 s | 0.312 s | 0.0482 s | 0.0423 s | 0.0147 s | 0.0825 s | 0.0590 s | 0.0614 s | 0.3081 s | 7.1% |

The instrument is not free, and its price is what says how far the table above can be
read. The same three loads with the timing off against on, best of five:

| particles a slice | walk off | walk on | busy off | busy on | busy |
|---|---|---|---|---|---|
| 1024 | 0.611 s | 0.912 s | 0.1120 s | 0.2920 s | 2.61x |
| 8192 | 0.683 s | 0.882 s | 0.1240 s | 0.2580 s | 2.08x |
| 131072 | 1.700 s | 1.800 s | 0.2740 s | 0.3170 s | 1.16x |

So the absolute seconds in the per-pass table are the instrument's rather than a
production run's, and at the two lighter loads the encoder boundaries cost about as much
as the work. The shares survive that better than the seconds do, with one bias to keep in
mind: a pass bears one boundary per encoder, so the deposit's one encoder a step carries
less of the overhead than the filter's four or the solve's four, and the deposit's share
is therefore a floor. Read that way the deposit is 7.7 percent of the six at 1024 and 8192
particles a slice and 26.8 percent at 131072, where the overhead is smallest and the
figure is firmest. The transform, filter and solve together are two thirds of a step at
production loads, which is the same conclusion the load scan of [](#perf-device) reaches
from the outside: the step is dominated by the grid and not by the particles until the
load is very high.

What the fixed-point deposit costs. The same three loads and the same deck as the two
tables above, timing off on both sides so the instrument's own boundaries are out of it,
best of five. The float deposit is the run before the accumulator changed.

| particles a slice | busy, float | busy, fixed point | busy | walk, float | walk, fixed point | walk |
|---|---|---|---|---|---|---|
| 1024 | 0.1120 s | 0.1390 s | 1.24x | 0.611 s | 0.649 s | 1.06x |
| 8192 | 0.1240 s | 0.1640 s | 1.32x | 0.683 s | 0.733 s | 1.07x |
| 131072 | 0.2740 s | 0.3550 s | 1.30x | 1.700 s | 1.800 s | 1.06x |

A quarter to a third more device time, and six to seven percent more wall clock, because
the step is dominated by the host's dispatch at these sizes rather than by the device.
Three passes pay it. The deposit does two atomic adds a component where it did one, since
the 64-bit accumulator is two 32-bit words. The source clear writes twice the words. And
the pass that converts reads four integers where it read two floats, which is the first of
the filter's source passes with the filter on and the solve's last pass with it off. What
that buys is in [](validation.md#val-device): two runs of one deck identical on every
array of the statistics file, where 54 of 127 differed before.

The field set, on the planar segment of the harmonics check (3.96 m, 88 steps) with the same 96 x 8192 window at `ngrid` 256 and `comb_ds_save = -1`, since a helical segment couples only the fundamental and the set has nothing to carry there:

| case | CPU, 12 threads | device wall | device busy | ratio |
|---|---|---|---|---|
| fundamental only | 3.344 s | 0.269 s | 0.185 s | 12.4x |
| harmonics 1, 3 | 5.249 s | 0.503 s | 0.365 s | 10.4x |
| two polarizations (tilt 0.3) | 4.502 s | 0.418 s | 0.296 s | 10.8x |

A second member costs the device 87% more wall clock and the CPU 57% more. A second plane costs the device 55% more and the CPU 35% more. The device pays the extra planes' transforms in full and the gather and deposit loops over members per particle, where the CPU's field solve is already the smaller share, so the ratio narrows by a fifth and no more. Production build, min of two runs, walk time.

The readback schedule is the other knob: the rows above read back at the element end only, and `comb_ds_save = 0` (a stats row every step) pays a beam-and-field readback per step, the same trade the reference manual prices as `output_step`. Measured at that finest comb on the 96 x 8192 case, with the readback's conversion loops first serial and then parallel over slices:

| comb = 0, walk | host stats | readback | per stats row |
|---|---|---|---|
| serial readback: 1.231 s | 0.425 s | 0.409 s | ~10 ms |
| parallel readback: 0.897 s | 0.412 s | 0.181 s | ~7 ms |

For calibration, the same deck is 0.303 s at `comb_ds_save = -1` and 4.41 s on the CPU at comb 0, so the finest comb costs the device path a factor of three and the CPU path a tenth. The parallelization recovered the conversion arithmetic and no more, and the thread scan says why: the readback share is 0.371 s at one thread, 0.195 at four and 0.185 at twelve, a memory-bandwidth plateau. What remains is the bytes, roughly 70 MB of FP32 transfers plus their FP64 stores per record, about 2 ms on this machine, and no thread count reduces bytes. That number is the honest input to the device-side moments question: computing the per-slice stats reductions on the GPU would replace the megabytes with some fifty floats per slice, which is the only way under the plateau, at the cost of FP32 accumulation moving the stats digits. Whether that trade is ever taken is decided against this table, and the analysis side (Twiss, emittance, the element-end bunch parameters) stays on the CPU in every variant.
