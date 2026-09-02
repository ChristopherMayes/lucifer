---
title: Performance
short_title: Performance
---

Where a run spends its time, measured. Every number here is machine-local and carries the machine and the build that produced it, which is why this page is hand written and stays out of `generated/`.

The instrument is `code/fel_timer_mod.f90`. It accumulates wall clock per phase at region boundaries and the footer of every run prints the table. The phases partition the walk, so the fractions sum to it and the remainder has its own row named `unaccounted`. Measured below, that remainder is 0.02% of the walk, so the partition is the walk.

Timers cannot see inside a parallel region. The deposit and its FFT interleave per slice, and so do the field gather and the RK4, so a clock call between them would measure the loop it perturbs. That split comes from a sampling profiler and appears in its own section.

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

Apple M3 Max, 12 performance cores, production build, full 6-FODO Benchmark1-SASE line, 2048 particles per slice, `ngrid` 255, 12 threads.

| phase | 96 slices | 504 slices |
|---|---|---|
| field solve (deposit and FFT) | 15.435 s, 64.5% | 82.072 s, 64.4% |
| particle push (transverse maps, RK4) | 4.485 s, 18.8% | 23.608 s, 18.5% |
| `= FEL step`, the sum of the phases inside it | 19.966 s, 83.5% | 105.731 s, 83.0% |
| stats and diag | 3.201 s, 13.4% | 18.223 s, 14.3% |
| field drift through the breaks | 0.510 s, 2.1% | 2.556 s, 2.0% |
| seam interlude | 0.103 s, 0.4% | 0.574 s, 0.5% |
| slippage | 0.102 s, 0.4% | 0.109 s, 0.1% |
| undulator prep (plans, kernel cache) | 0.046 s, 0.2% | 0.051 s, 0.0% |
| element end | 0.032 s, 0.1% | 0.170 s, 0.1% |
| unaccounted | 0.003 s, 0.0% | 0.006 s, 0.0% |
| walk | 23.917 s | 127.370 s |

The `= FEL step` row is derived rather than measured, which its label says: it is the sum of the undulator prep, the space-charge profile, the particle push and the field solve, so it is the share of the walk the undulator segments own. It is left out of the leaf sum, and the other rows still add to the walk minus the remainder.

The shares are the result and the third decimal is noise. Repeating the 96-slice run gives a walk of 23.583 s against 23.917 s and a field solve of 64.0% against 64.5%, so the run-to-run spread is about 1.4% on the total and half a point on a share.

The shares hold across a factor of 5.25 in slice count. The run row against the walk row is everything outside the walk: the parse, the beam and field build, the slippage schedule, the final dumps and the stats file. That is 0.2 s at 96 slices and 0.8 s at 504.

Slippage is 0.4% of the walk and 0.1% at the larger size. [](validation.md#val-parallelism-openmp-over-slices) named the slippage rotation as part of the serial cost that capped scaling, which the measurement does not support. The field solve and the particle push are where the time is, and both are slice-parallel.

(perf-the-particles-against-the-grid)=
## The particles against the grid

The field solve and the particle push scale differently, so their ratio is a property of the configuration and not of the code. Same machine and line, 96 slices, 12 threads.

| particles per slice | particle push | field solve | walk |
|---|---|---|---|
| 512 | 1.217 s, 6.3% | 14.878 s, 77.3% | 19.250 s |
| 2048 | 4.485 s, 18.8% | 15.435 s, 64.5% | 23.917 s |
| 8192 | 17.619 s, 42.6% | 16.677 s, 40.3% | 41.375 s |

The push is linear in the particles (1.217, 4.485, 17.619 for 512, 2048, 8192) and the field solve is nearly flat (14.878, 15.435, 16.677). Two FFTs of $255^2$ points per slice per step do not care how many particles deposited into the grid. The crossover is near 8192 particles per slice at `ngrid` 255.

That flatness also splits the field solve without a profiler. Fitting the 512 and 8192 points to $C + d \cdot n$ gives $C = 14.76$ s and $d = 2.34 \times 10^{-4}$ s per particle. At 2048 particles per slice that predicts 14.76 s of transform and propagator plus 0.48 s of deposit, so 15.24 s against the 15.435 s measured, a 1.3% residual on a two-point fit. The transform and propagator are 95.6% of the solve by that fit. The sampling profile below gives 95.2% from independent evidence, which is the same number.

[](validation.md#val-the-particlepath-cost-measured) measured a real 131-slice case where the particle path dominated at 8192 particles per slice while the field FFTs stayed constant. That is this table, and the two agree.

(perf-thread-scaling)=
## Thread scaling

96 slices, 2048 particles per slice, full line, same machine and build. The implied serial fraction inverts Amdahl's law, $f = (n/S - 1)/(n - 1)$, and bounds this configuration rather than the code.

| threads | walk | speedup | efficiency | implied serial |
|---|---|---|---|---|
| 1 | 219.1 s | 1.00 | the reference | |
| 2 | 113.6 s | 1.93 | 96% | 3.7% |
| 4 | 58.5 s | 3.74 | 94% | 2.3% |
| 8 | 32.4 s | 6.76 | 84% | 2.6% |
| 12 | 23.9 s | 9.16 | 76% | 2.8% |

[](validation.md#val-parallelism-openmp-over-slices) records 3.97x at 8 threads on 32 slices in a debug build and says that is the floor of the scaling rather than its ceiling. At 96 slices in a production build the same 8 threads give 6.76x, so the claim holds.

The nameable serial phases account for 0.6% of the walk, well under the 2.8% the speedup implies. The rest is not serial code. Two parallel regions are spawned per integration step (1068 steps, so about 4300 region entries), and the sampling profile below puts 10.0% of its samples in the thread wait, which at 76% efficiency is the barrier and not a lock.

A steady-state run gets nothing from threads. The parallelism is over slices, and a steady-state window has one: the unaveraged example measures 21.94 s at 1 thread and 21.76 s at 12.

(perf-the-sampling-split)=
## The sampling split, averaged mode

96 slices, 2048 particles per slice, 12 threads, 15 s of samples at 1 ms, 116236 samples. Exclusive cost, categorized by symbol.

| | share |
|---|---|
| FFT transform (libfftw3, `wavefront_fft2`) | 61.3% |
| thread wait | 10.0% |
| libm sin and cos | 8.5% |
| transverse maps (`quad_mat2_calc`, and its `cexp`) | 5.9% |
| deposit and propagator multiply | 5.7% |
| diagnostics | 4.6% |
| longitudinal RK4 and ODE | 3.6% |
| slippage and drift | 0.1% |
| everything else | 0.4% |

An earlier estimate put the field solve near 10% of a run, and at these parameters it is 64.4% by the phase timer and 61.3% in the transform alone by samples. That estimate came from a configuration with many more particles per slice, where the table above shows the balance reversing.

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
