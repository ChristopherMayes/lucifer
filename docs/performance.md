---
title: Performance
short_title: Performance
---

This page is the measured record of where a run spends its time. Every number here is machine-local and carries the machine and the build that produced it. That is why this page is hand written and stays out of `generated/`.

The instrument is `src/fel_timer_mod.f90`. It accumulates wall clock per phase at region boundaries, and the footer of every run prints the table. The phases partition the walk, so the fractions sum to it and the remainder has its own row named `unaccounted`. That remainder, measured below, is 0.02% of the walk, so the partition is the walk.

Timers cannot see inside a parallel region. The deposit and its FFT interleave per slice, and so do the field gather and the RK4, so a clock call between them would measure the loop it perturbs. That split comes from a sampling profiler and appears in its own section.

The tables below are also the calibration of the cost estimate every run prints in its header. `src/fel_cost_mod.f90` holds the fitted rates, one for the particle push and two for the field solve, and `tests/scripts/fit_cost_rates.py` fits them to the logs the runs below wrote. Over the sixteen runs of this page the worst residual on the walk is 16.5 percent, on the four- and eight-thread rows of the thread sweep, where one serial fraction cannot carry a sweep whose efficiency is 92 percent at four threads and 94 at twelve. Every other row is inside 10 percent. The estimate is for this machine and says so on the line, and the run's footer prints the measured walk beside it. Re-measuring this page means re-running that script and pasting its five constants back.

Two toolchains sit side by side on this page, and every table carries one line naming the toolchain, the machine and the day it was measured. The tables re-measured on 2026-10-24 were taken under gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits of [](validation.md) belong to since 2026-10-23. The tables kept from September 2026 were taken under gfortran 13.3 and the clang of 2025-06, and no later measurement has replaced them. Where a kept table and a re-measured one are compared, the comparison crosses the toolchain and the text says so.

(perf-reproducing)=
## Reproducing this

```
./tests/run_perf_benchmark.sh --phases
```

The phases mode needs no genesis4 and no MPI. It compares the code against itself, so it starts from the tracker's own shot-noise quiet start rather than an imported dump. It runs the profile at two slice counts, sweeps the thread count over the smaller one, and prints both tables. `--npart` sets the particles per slice, `--big-slices` the larger case, `--ngrid` the transverse grid, 256 by default, and `--no-sweep` leaves the sweep out. The particles-against-the-grid table is the same command at `--npart 512` and `--npart 8192`, and the transform-length table is `--grid-scan "128 151 255 256 257"`, which runs the small configuration at each grid three times and takes the best walk. The kept tables that state `ngrid` 255 were measured before the benchmark moved to power-of-two grids, and on this toolchain the two grids are within 0.3 percent of each other ([](#perf-the-transform-length)).

The device table runs as:

```
./tests/perf/device_table.sh
```

Its four input files are the `*.in` files beside it. The estimate-against-the-clock table reads the footers of runs made one at a time:

```
./tests/scripts/estimate_vs_clock.py tier1=<log> tier1_unavg=<log> ...
```

The sampling split is taken on macOS with:

```
OMP_NUM_THREADS=12 build-release/bin/lucifer <input>.toml &
sample $(pgrep -n lucifer) 15 1 -file prof.txt -mayDie
```

The averaged input file is the 96-slice one the phases mode writes, and the unaveraged one is `tests/perf/unavg_sample.in`, run at one thread and sampled for 20 s.

Read the `Sort by top of stack` section, which is exclusive cost per symbol. The call graph is per thread and fragments the same information across dozens of stacks.

The vectorization audit runs as:

```
./tests/scripts/vec_audit.sh
```

It takes the compile command from the production build's own makefile, so the report always describes the binary that ships. Build production first.

(perf-the-phase-profile)=
## The phase profile

The machine and configuration are an Apple M3 Max with 12 performance cores, a production build, the full 6-FODO Benchmark1-SASE line, 2048 particles per slice, `ngrid` 256, 12 threads and the source filter at its default. Each column is the mean of two runs. The `2026-09-07` columns are the same configuration at `ngrid` 255 under gfortran 13.3, the record this one replaces.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| phase | 96 slices | 2026-09-07 | 504 slices | 2026-09-07 |
|---|---|---|---|---|
| field solve (deposit and FFT) | 25.152 s, 74.9% | 21.687 s, 71.5% | 124.885 s, 75.4% | 118.060 s, 71.7% |
| particle push (transverse maps, RK4) | 4.199 s, 12.5% | 4.735 s, 15.6% | 21.031 s, 12.7% | 24.509 s, 14.9% |
| `= FEL step`, the sum of the phases inside it | 29.406 s, 87.6% | 26.468 s, 87.3% | 145.972 s, 88.2% | 142.645 s, 86.6% |
| stats and diag | 3.326 s, 9.9% | 3.140 s, 10.4% | 16.083 s, 9.7% | 18.559 s, 11.3% |
| field drift through the breaks | 0.619 s, 1.8% | 0.490 s, 1.6% | 2.783 s, 1.7% | 2.597 s, 1.6% |
| seam interlude | 0.096 s, 0.3% | 0.096 s, 0.3% | 0.469 s, 0.3% | 0.578 s, 0.4% |
| slippage | 0.095 s, 0.3% | 0.104 s, 0.3% | 0.092 s, 0.1% | 0.106 s, 0.1% |
| undulator prep (plans, kernel cache) | 0.056 s, 0.2% | 0.045 s, 0.1% | 0.056 s, 0.0% | 0.075 s, 0.0% |
| wake | 0.000 s, 0.0% | | 0.000 s, 0.0% | |
| element end | 0.030 s, 0.1% | 0.032 s, 0.1% | 0.142 s, 0.1% | 0.159 s, 0.1% |
| unaccounted | 0.003 s, 0.0% | 0.004 s, 0.0% | 0.006 s, 0.0% | 0.007 s, 0.0% |
| walk | 33.574 s | 30.333 s | 165.546 s | 164.650 s |

Across the toolchain the push is faster and the solve is slower. The push, which is this code's own arithmetic through gfortran, takes 11 percent less at 96 slices and 14 percent less at 504. The field solve, whose transform is FFTW's and not the compiler's, takes 16 percent more at 96 slices and 6 percent more at 504, and the grid moved from 255 to 256 between the two records, which the transform-length table below puts at 0.3 percent on this toolchain. The 96-slice walk is up 11 percent and the 504-slice walk is unchanged within the spread. The shares hold: the solve is three quarters of the walk at both sizes. The transform-length table below names the move: the June-2024 fftw under today's binary reproduces the September walk at `ngrid` 255 to half a percent, so the solve's extra seconds are the grid the benchmark moved to and the library's price for a power of two, and nothing in the code.

The source filter, which has been the default since 2026-09-07, adds a forward transform of the source per field per step, where an unfiltered solve adds the source in real space and transforms twice rather than three times ([](fel-physics.md#sec-source-filter)). The same input file, measured on 2026-09-07 on the 96-slice configuration with the filter off and on, gives:

*Measured 2026-09-07: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| | filter off | filter on | ratio |
|---|---|---|---|
| field solve | 15.011 s | 21.339 s | 1.42 |
| particle push | 4.484 s | 4.441 s | 0.99 |
| walk | 23.392 s | 29.646 s | 1.27 |

Each column is the mean of two runs. The filter costs 27 percent of this walk and 42 percent of its field solve, and it leaves the push alone. The 8 to 14 percent [](fel-physics.md#sec-convergence) quotes is the cost at the derived grid on a shorter line, where the transform is a smaller share of a run. This line at `ngrid` 255 is where the transform dominates, so it is where the filter costs most. The filter-off walk of 23.392 s reproduces the 23.917 s of the `was` column, so nothing else in the code moved between the two measurements.

The `= FEL step` row is derived rather than measured, as its label says. It is the sum of the undulator prep, the space-charge profile, the particle push and the field solve, so it is the share of the walk the undulator segments own. It is left out of the leaf sum, and the other rows add to the walk minus the remainder.

The shares are the result and the third decimal is noise. The two 96-slice runs behind the table give walks of 33.182 s and 33.965 s, and the two 504-slice runs 164.300 s and 166.791 s. The run-to-run spread is therefore about 2% on a walk and a few tenths of a point on a share. The phase timers' own cost, two `system_clock` calls per region, is below that spread.

The shares hold across a factor of 5.25 in slice count. The run row against the walk row is everything outside the walk: the parse, the beam and field build, the slippage schedule, the final dumps and the stats file. That is 0.2 s at 96 slices and 0.9 s at 504.

Slippage is 0.4% of the walk and 0.1% at the larger size. [](validation.md#val-parallelism-openmp-over-slices) lists the slippage rotation among the serial operations between steps, and it is too small to cap the scaling. The field solve and the particle push are where the time is, and both are slice-parallel.

(perf-the-transform-length)=
## The transform length

Every input file in the benchmark runs a power-of-two transverse grid, so one input file runs on the CPU and on the Metal field solver unchanged, and the derived count rounds up to one ([](input-reference.md#param-wavefront-grid)). What that costs on the CPU is measured here. The input file is the phases input file above, 96 slices of 2048 macroparticles over the full line, 12 threads, production build, varying only `grid_n_pts`, three runs at each grid and the best taken. The last column is the walk over the grid points, and at the small grids the particle work, which does not move with the grid, is most of it.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| ngrid | factors | points | walk | seconds a point |
|---|---|---|---|---|
| 128 | 2^7 | 16384 | 12.563 s | 7.67e-4 |
| 151 | prime | 22801 | 22.575 s | 9.90e-4 |
| 255 | 3 x 5 x 17 | 65025 | 32.418 s | 4.99e-4 |
| 256 | 2^8 | 65536 | 32.308 s | 4.93e-4 |
| 257 | prime | 66049 | 53.614 s | 8.12e-4 |

The two large grids cost the same on this toolchain: 256 is 0.3 percent faster than 255 for 0.8 percent more points. On gfortran 13.3 and the FFTW beside it, measured on 2026-09-07 on a one-slice input file, the power of two was the slower of the two at 1.17x, which the previous version of this table recorded and the `2026-09-07` columns above still reflect. The library is the difference, since the transform is FFTW's and the compiler never sees it, and the test that says so is cheap: the same binary with the June-2024 build of fftw 3.3.10 placed first on its library path, the build the September environment held, runs this input file in 30.481 s at 255 and 32.952 s at 256, best of three, where today's fftw 3.3.11 gives 32.418 and 32.308. The old library reproduces the September walk at 255 to half a percent, and it is slower than the new one at 256 by 2 percent, so the 11 percent between the two phase-profile records is the grid the benchmark moved to, 255 to 256, through a library that paid 8 percent for the power of two, and not the compiler. A prime is worse than either on both toolchains, 1.66x against 256 here, which is the transform falling back to its general algorithm, and 151 costs 1.8x what 128 does for 1.4x the points. What the grid convention buys is an input file that moves between the backends, and on this toolchain it costs nothing.

## The particles against the grid

The field solve and the particle push scale differently, so their ratio is a property of the configuration and not of the code. The machine and line are the same, at 96 slices and 12 threads.

Each row is the mean of two runs at that count.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| particles per slice | particle push | field solve | walk |
|---|---|---|---|
| 512 | 0.996 s, 3.5% | 24.049 s, 84.8% | 28.366 s |
| 2048 | 4.199 s, 12.5% | 25.152 s, 74.9% | 33.574 s |
| 8192 | 15.871 s, 32.7% | 25.702 s, 53.0% | 48.492 s |

The push is linear in the particles (0.996, 4.199, 15.871 for 512, 2048, 8192) and the field solve is nearly flat (24.049, 25.152, 25.702). Three FFTs of $256^2$ points per slice per step cost the same whatever the particle count deposited into the grid. The crossover sits past 8192 particles per slice at `ngrid` 256, where the push is a third of the walk and the solve half.

That flatness also splits the field solve without a profiler. Fitting the 512 and 8192 points to $C + d \cdot n$ gives $C = 23.94$ s and $d = 2.15 \times 10^{-4}$ s per particle. The slope is the deposit, and it is within 6 percent of the 2.28e-4 of 2026-09-07, as it should be, since the deposit is this code's own arithmetic and the filter does not touch how the source is built. At 2048 particles per slice the fit predicts 23.94 s of transform and propagator plus 0.44 s of deposit, 24.38 s against the 25.152 s measured, a 3 percent residual on a two-point fit. The transform and propagator are 98.2 percent of the solve by that fit against 98.6 percent by the sampling profile below, the transform's 73.9 percent of the walk over the solve's 74.9. The fit's deposit term is 2 percent of the solve, so the two-point $C$ carries the solve's own spread almost whole.

[](validation.md#val-the-particlepath-cost-measured) measured a real 131-slice case where the particle path dominated at 8192 particles per slice while the field FFTs stayed constant. That is this table, and the two agree.

(perf-thread-scaling)=
## Thread scaling

The configuration is 96 slices of 2048 particles per slice through the full line, on the same machine and build. The implied serial fraction inverts Amdahl's law, $f = (n/S - 1)/(n - 1)$, and bounds this configuration rather than the code.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| threads | walk | speedup | efficiency | implied serial |
|---|---|---|---|---|
| 1 | 372.3 s | 1.00 | the reference | |
| 2 | 195.7 s | 1.90 | 95% | 5.1% |
| 4 | 101.6 s | 3.66 | 92% | 3.0% |
| 8 | 52.8 s | 7.05 | 88% | 1.9% |
| 12 | 33.2 s | 11.22 | 94% | 0.6% |

The twelve-thread speedup is 11.22x against the 9.59x of 2026-09-07 under gfortran 13.3, and the one-thread walk is 372.3 s against 290.8. The two moves are one: the one-thread walk grew 28 percent across the toolchain and the twelve-thread walk 9.5 percent, so the parallel part grew more than the serial part did, and the sweep's efficiency at twelve threads rose from 80 to 94 percent. The implied serial fraction is not one number on this toolchain, 5.1 percent at two threads and 0.6 at twelve, which is why one Amdahl constant carries the sweep to 16.5 percent at worst. The estimate takes 1.05 percent, fitted over the whole sweep.

[](validation.md#val-parallelism-openmp-over-slices) records 3.97x at 8 threads on 32 slices in a debug build and says that is the floor of the scaling rather than its ceiling. At 96 slices in a production build the same 8 threads give 7.05x, so the claim holds.

The nameable serial phases account for 0.6% of the walk. The rest is the cost of the parallel regions themselves. Two are spawned per integration step (1068 steps, so about 4300 region entries). The sampling profile below puts 3.7 percent of its samples in the thread wait, which at 94 percent efficiency is the barrier and not a lock.

A steady-state run gets nothing from threads. The parallelism is over slices, and a steady-state window has one: the unaveraged example measured 21.94 s at 1 thread and 21.76 s at 12 on 2026-09-01, before the hoist and the interchange took the mode's step to 60 percent of itself.

The field-free drift is parallel over slices for the same reason. A time-dependent line makes one `wavefront_drift` call per field-free element. With that call's slice loop serial, a 590-slice, 256-point case took about 0.5 s per element on one core, the largest serial block of the walk.

(perf-the-sampling-split)=
## The sampling split, averaged mode

The configuration is 96 slices of 2048 particles per slice on 12 threads with the filter at its default, sampled for 15 s at 1 ms, 92446 samples. The table is exclusive cost, categorized by symbol. The `2026-09-07` column is the same measurement under gfortran 13.3, and its `was` column the same with the filter off.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| | share | 2026-09-07 | was |
|---|---|---|---|
| FFT transform (libfftw3, `wavefront_fft2`) | 73.9% | 68.8% | 61.3% |
| diagnostics | 5.3% | 5.5% | 4.6% |
| transverse maps (`quad_mat2_calc`, and its `cexp`) | 4.2% | 4.9% | 5.9% |
| longitudinal RK4 and ODE | 4.0% | 3.7% | 3.6% |
| thread wait | 3.7% | 5.1% | 10.0% |
| libm sin and cos | 3.0% | 6.5% | 8.5% |
| libm, internals the sampler could not name | 2.8% | | |
| deposit and propagator multiply | 2.8% | 4.9% | 5.7% |
| slippage and drift | 0.0% | 0.1% | 0.1% |
| everything else | 0.3% | 0.6% | 0.4% |

The transform is up 5.1 points across the toolchain and everything this code compiles is down, which is the phase table's split read by symbol: the push's share fell and the solve's rose. The sampler attributes 2.8 percent to addresses inside libm it cannot name, which the earlier toolchain's libm did not show and which belong with the sin, cos and exp lines. The thread wait is 3.7 percent against 5.1, and the thread sweep above agrees. On 2026-09-07 the transform took the filter's third pass and rose 7.5 points over the `was` column, which was the whole of what moved then.

An estimate of the field solve near 10% of a run holds only at many more particles per slice, where the table above shows the balance reversing. At these parameters it is 74.9% by the phase timer and 73.9% in the transform alone by samples.

`quad_mat2_calc` is Bmad's own quadrupole map, and its cost is a complex exponential per particle. `cexp` and the `exp` beneath it are the whole 5.9%. The libm sin and cos are the deposit's phase factor and the ODE's, through `fel_sincos`.

(perf-the-sampling-split-unaveraged)=
## The sampling split, unaveraged mode

The configuration is one slice of 16384 particles at `ngrid` 129, run serially (the mode's parallelism is over slices and there is one). The phase table gives the unaveraged step 99.7% of the walk, so this is the mode. Both columns are 20 s of samples at 1 ms, and the totals differ because the hoisted run finishes sooner inside that window.

*Measured 2026-09-02 for the first two columns under gfortran 13.3, and 2026-10-24 for the third under gfortran 15.3, Apple M3 Max, production build, one thread.*

| | before the hoist | after | 2026-10-24 |
|---|---|---|---|
| libm transcendentals | 52.1% | 15.3% | 18.7% |
| substep RK4 (`unavg_push_all`, `unavg_ode`) | 17.6% | 35.3% | 32.7% |
| undulator field, this code's own arithmetic | 8.9% | 27.7% | 15.0% |
| FFT transform and `fel_field_diffract` | 13.0% | 17.8% | 21.9% |
| ramp envelope (`fel_und_envelope`) | 4.8% | 0.0% | 0.0% |
| deposit, this code's own arithmetic | 1.2% | 1.9% | 2.1% |
| everything else | 2.4% | 1.8% | 9.6% |
| total samples | 13302 | 7964 | 8085 |

The third column is the same input file, `tests/perf/unavg_sample.in`, under gfortran 15.3 on 2026-10-24, with the interchange of 2026-09-02 in the binary as well as the hoist. Its `everything else` holds 7.7 points the sampler attributes to the step's own parallel-region body, `fel_unavg_step`, where the newer compiler inlines what the older one called, and the undulator field's line is lower by the same amount. The walk was 12.2 s against the 21.94 s of 2026-09-01.

Half of this mode was libm, and none of it went through `fel_sincos`. `fel_unavg_bfield` called `cos(und%ku * s)` and `sin(und%ku * s)` as two separate intrinsics, and `fel_und_envelope` called two more, all four per particle per RK stage. The argument depends only on the substep position and is therefore the same for every particle in the loop. The four s-dependent factors now arrive as arguments, evaluated once per stage, and `fel_unavg_bfield` no longer takes `s` at all.

The share row hides how much moved, because the denominator moved too. In absolute samples libm fell from 6936 to 1221, which is 82% of the transcendental work gone, and the step's wall clock fell by the amounts in the next table.

The before column's label covers every libm transcendental the run entered, `exp` and `cexp` and the shared reduction helper included. Only 13.9% of the run was `sin` and `cos` proper. The hoist recovers most of it either way, because the envelope's own trig was per particle for the same reason.

(perf-the-hoist-measured)=
## The hoist, measured

The machine and build are the same, on one thread, and each entry is the median of five runs at 2048 particles and three at 16384. The change is bit-identical by construction and checked as such: the diag and ledger files of the unaveraged example are byte-identical before and after at 1 and 8 threads.

*Measured 2026-09-02: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| configuration | before | after | change |
|---|---|---|---|
| 2048 particles per slice | 5.293 s | 4.507 s | -14.8% |
| 16384 particles per slice | 21.842 s | 15.806 s | -27.6% |

The gain grows with the particle count because what it removes is per-particle and what remains at fixed count is the per-substep FFT.

The loop order is the next lever after the hoist. The substep RK4 and the undulator field arithmetic were 63% of the mode's samples, four RK stages inside each particle of a `do ip` loop. The push runs stage-major: per stage, one loop over the slice's particles into per-particle stage arrays, then one loop for the RK4 combination. Per particle the operations, values and order are unchanged, so the results are byte-identical, checked the same way as the hoist.

*Measured 2026-09-02: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| configuration | before | after | change |
|---|---|---|---|
| 2048 particles per slice | 4.577 s | 4.217 s | -7.9% |
| 16384 particles per slice | 16.073 s | 13.189 s | -17.9% |

The gain is locality and call overhead, and not vectorization. The stage loops do not vectorize, because each iteration calls `unavg_ode`. Inlining the field into the stage loop fails the identity check. The compiler fuses the inlined arithmetic differently from the called form, the results differ in the last bits, and a moved bit is a moved digit somewhere downstream. So the call stays. Vectorizing this loop needs an arithmetic core whose bits are fixed by its formulation rather than by fusion. That is the territory of the FP32 twin, where the reference is measured divergence rather than identity.

After the interchange the mode's samples split 37.5% substep RK4, 21.4% FFT, 18.8% undulator field arithmetic, 18.7% the deposit's own transcendentals, 2.2% deposit weights.

(perf-the-sincos-ceiling)=
## The sincos ceiling

The probe is an inline (sin, cos) pair with no libm call, Cody-Waite reduction and the fdlibm kernel polynomials, swapped into `fel_sincos.c` for the measurement and reverted. A first probe build recorded the argument range of that call site: $[-159.8, 19.0]$ on the 96-slice averaged configuration and $[-127.4, 17.1]$ on the unaveraged example. Both sit well inside the range a two-term Cody-Waite reduction handles. The probe run reproduces the reference exit power, pulse energy and bunching to every printed digit. The physics is the same, so the timing difference is the call.

*Measured 2026-09-02: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| case | reference | probe | change |
|---|---|---|---|
| 96 slices, 2048 particles, 12 threads | 23.917 s | 23.158 s | -3.2% |
| 96 slices, 8192 particles, 12 threads | 41.375 s | 38.197 s | -7.7% |
| unaveraged, 16384 particles, serial | 21.937 s | 21.848 s | -0.4% |

The averaged rows bound the libm call from below at 3.2% of a run at 2048 particles per slice, against the 8.5% the sampler attributes to sin and cos. The polynomial is cheaper than the call and not free. The unaveraged row barely moves. That mode never reaches this file. Its own transcendentals were the undulator field's and the ramp envelope's intrinsics, which the hoist above removed instead.

Neither probe build is a candidate implementation, and the audit rules the inline pair out. `tests/scripts/sincos_audit.c` sweeps it against libm:

*Measured 2026-09-02: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| range | worst error | arguments differing |
|---|---|---|
| the measured range, $[-159.8, 19.0]$ | 6 ulp | 34% |
| $\lvert \theta \rvert < 10^4$ | 6 ulp | 34% |
| $\lvert \theta \rvert < 10^5$ | 10 ulp | 34% |
| $\lvert \theta \rvert < 10^6$ | 24 ulp | 34% |
| $\lvert \theta \rvert < 10^8$ | 2681 ulp | 34% |

The shim it would replace was admitted at one ulp on $2 \times 10^{-6}$ of arguments, 73 mismatches in a 44M-point sweep ([](validation.md#val-the-particlepath-cost-measured)). This candidate is six ulp on a third of them, five orders of magnitude more arguments and six times the error, so it cannot be recorded as a 1-ulp change. The error is the two-term reduction rather than the kernels. A reduction accurate enough to reach one ulp costs arithmetic that eats a 3.2% gain. The range is also not a constant. $\theta$ carries the common phase accumulating along the line, so a longer line moves toward the degradation with nothing to detect it.

The version of this that could pay is a vectorized pair inside a loop that vectorizes, where one reduction amortizes over a whole vector. The averaged path has no such loop, for the reasons in the audit section below.

(perf-the-vectorization-audit)=
## The vectorization audit

The compiler is gfortran 15.3.0 from conda-forge at the production flags, which are `-O2` with per-file `-O3` on `fel_track_mod`, read from the production build's own `flags.make`. At `-O2` gfortran runs its vectorizer at the cheap cost model, and no `-march` is set, so the code runs on any arm64.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| file | vectorized | not vectorized |
|---|---|---|
| `fel_track_mod.f90` | 47 | 251 |
| `fel_unaveraged_mod.f90` | 12 | 77 |
| `fel_collective_mod.f90` | 1 | 49 |

None of the hot per-particle loops is among the vectorized. The compiler rejects every one of them for the same reason, `unsupported control flow in loop`. The audit's report gives each loop's line, which moves with every edit to the file, so the table names the loop by its routine:

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| loop | routine | what blocks it |
|---|---|---|
| the particle loop | `fel_transverse_track` | `fel_apply_focus`, called per particle, branches inside |
| the particle loop | `fel_transverse_track_bmad` | the same, plus `quad_mat2_calc` |
| the particle loop | `fel_advance` | `faw`, `fel_sincos` and `fel_runge_kutta`, called per particle |
| the particle loop of the multi-harmonic twin | `fel_advance` | the same |
| the deposit's particle loop | `fel_field_step` | `if (.not. on_grid) cycle`, and the report also names `statement clobbers memory: fel_grid_weights` |

The aggregate reasons, most common first, are `vectorization is not profitable` (112 in `fel_track_mod`), `complicated access pattern` (33) and `unsupported control flow in loop` (31). The first is mostly short fixed-length loops over the 6-vector, which the compiler is right about. The second covers the loops in the table.

One pattern accounts for all five. A per-particle loop calls a routine that the compiler does not inline, the routine carries a branch, and the branch ends the loop's chance of vectorizing. The deposit's `on_grid` guard is the clearest case: it is a mask written as a `cycle`.

(perf-what-the-numbers-say-to-do)=
## What the numbers say to do

The items below are ordered by measured share, for the averaged mode at 2048 particles per slice and `ngrid` 255, and read from the `was` columns above, which were measured with the filter off.

The FFT is 61% of the run, and it is in FFTW rather than in this code. What this code controls is how many transforms it asks for: two per slice per field per step without the filter and three with it, at the grid the input file sets. `ngrid` is an input file parameter, and $255^2$ is 65025 points against 2048 particles, so a convergence study on `ngrid` is worth more here than any change to the deposit.

The libm transcendentals are 8.5% of the averaged mode. They were half the unaveraged mode and the hoist took 82% of that, so what remains there is the RK4 itself. Bmad's `cexp` in the transverse map is a further 5.9% of the averaged mode, and a quadrupole map has a real closed form.

The per-particle loops do not vectorize, and the blocker is uniform: an un-inlined call with a branch in it. Three of the five in the table above have reasons that stand.

The deposit's `on_grid` guard cannot become a mask. `fel_grid_weights` leaves `ix` and `iy` unassigned when it reports off-grid, so a masked body would index `crsource` with an undefined subscript. Making that safe means clamping the indices, which is no longer the same program.

`fel_apply_focus` is called by `fel_transverse_track`, which runs only under `transport_model = "genesis"`, and by the transcribed interlude, which runs only under `interlude_model = "genesis"`. Both default to `bmad`, no example sets either, and only the comparison tiers select them, so restructuring that loop speeds up a validation-internal path. `faw`, on the production path, carries one branch on `und%sin_t`, which is a property of the element and already loop-invariant.

`fel_advance`'s loops call `fel_sincos` and `fel_runge_kutta` per particle, so `!$omp simd` cannot apply to them until the sincos itself is vectorizable.

Thread scaling at production slice counts is 9.16x on 12 cores with the filter off, and the nameable serial cost is 0.6% of the walk, so there is little to win by removing serial work.

The single-precision question is measured. `run.fp32_check` steps an FP32 twin of the averaged advance beside the FP64 path and records the divergence per quantity. The twin carries the reformulations FP32 forces (offset energy, residual phase, difference-form detuning) and a runtime guard on the one failure that is silent, a residual too coarse for the per-step increment. The measured levels are in [](validation.md#val-fp32-lockstep). The field solve has its own twin there: FP32 deposit, single-precision transform pair, rounded propagator. The end-to-end freerun number for the complete single-precision run is 1.0e-3 on exit power over a gain segment, and 2.7e-1 in deep saturation where phase decorrelation dominates. An FP32 six-vector is 24 bytes per particle against 48, which is the bandwidth argument for a device path. Whether an FP32 CPU mode ships is a decision against those recorded levels.

The unaveraged mode has its own twin, and it costs more than the averaged one because of the substep. It repeats the whole advance in single precision on every substep of a record step, sixty of them on the benchmark segment, and it carries a field record of its own through a single-precision transform pair on each. On that segment at 2048 macroparticles and grid 128 the walk is 6.1 s with `run.fp32_check = "lockstep"` and 3.5 s without. The instrumented run takes 1.75 times the walk of the same run without it. The input file is one slice, where the mode's parallelism is over slices, so the figure is a serial one. The instrument is a measuring device rather than a mode, and no production run turns it on. [](validation.md#val-unaveraged-fp32) records what it measures.

(perf-estimate)=
## The estimate against the clock

Every run prints an estimated walk in its header and the measured walk beside it in its footer. The rates are fitted to the tables above and are this machine's. The estimate is measured here on the eleven comparison tiers and four examples, each run on its own rather than beside the other keystone jobs. A machine running five jobs at once measures the load and not the input file.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| run | filter | measured | estimate | estimate/measured |
|---|---|---|---|---|
| tier1 | off | 0.325 s | 0.43 s | 1.32 |
| tier1_unavg | off | 13.3 s | 25.8 s | 1.94 |
| tier2_bmad | off | 4.1 s | 5.2 s | 1.27 |
| tier2_genesis | off | 4.0 s | 5.2 s | 1.30 |
| td1 | off | 0.72 s | 0.867 s | 1.20 |
| td2_bmad | off | 9.1 s | 10.4 s | 1.14 |
| td2_genesis | off | 9.1 s | 10.4 s | 1.14 |
| tdsase | off | 9.2 s | 10.4 s | 1.13 |
| tdsc | off | 0.768 s | 0.867 s | 1.13 |
| tdwk | off | 0.76 s | 0.867 s | 1.14 |
| weight_split | off | 0.502 s | 0.615 s | 1.23 |
| `sase` | on | 32.4 s | 31.4 s | 0.97 |
| `steady_state` | on | 6.0 s | 5.2 s | 0.87 |
| `taper` | on | 6.0 s | 5.2 s | 0.87 |
| `flash1` | on | 19.6 s | 17.7 s | 0.90 |

Every row is inside a factor of two, and the table is `tests/scripts/estimate_vs_clock.py` over the fifteen footers. The rates are fitted with the source filter on, which the measurement of 2026-09-07 put at 58 percent of the transform rate with the push and deposit rates unchanged. The work outside the walk (parsing, building the beam and field, the final dumps) is 1 percent of the walk on these runs, and the estimate leaves it out. The filter column is what separates the two blocks. The examples take the default and read between 0.87 and 0.97. Every tier turns the filter off, since each compares against a code that carries none ([](validation.md)). Their field solve is therefore the cheaper two-transform one, while the rate costing it is the filtered default's. That is the 1.13 to 1.32 the tier block sits at, a property of the configuration rather than of the machine. On 2026-09-07 the same fifteen rows read 0.82 to 1.87.

`tier1_unavg` is the highest at 1.94, for a second reason. The unaveraged mode resolves the quiver, so it integrates and deposits on twenty substeps per period. The count carries those substeps while the rates it multiplies are the averaged path's. `steady_state` and `taper` are the lowest at 0.87, one slice of 8192 particles over the full line, where the estimate's push rate, fitted at 2048 particles a slice, runs ahead of a push that is a third of the walk.

(perf-frames)=
## The frame series, measured

This section measures what a frame costs, on the same machine. The input file is the `sase` example, 96 slices of 2048 macroparticles at `ngrid` 255, one frame every 20 m.

*Measured 2026-09-09: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| | size | per unit |
|---|---|---|
| field frame, 96 x 255 x 255 complex128 | 95.3 MB | 1.02 kB per grid point |
| beam frame, 196608 macroparticles | 12.8 MB | 68.2 B per macroparticle |

The beam frame carries no overhead to trim. Its 68 bytes are six coordinates, the time the reference phase folds into, the time offset that places the particle's slice in the bunch, the reference momentum and the label. The uniform weight and the zero longitudinal position are constant records that cost nothing per particle. The field is seven times the beam and is where a series' bytes are.

Compression does not pay on this data. The field is chunked one slice at a time with the shuffle filter:

*Measured 2026-09-09: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| | size | ratio | read |
|---|---|---|---|
| stored, uncompressed | 95.3 MB | | 0.009 s |
| gzip level 1 | 73.2 MB | 1.30x | 0.135 s |
| gzip level 4 | 72.6 MB | 1.31x | 0.139 s |
| float32 | 47.6 MB | 2.00x | |

A saturated FEL field is not sparse. The wide-angle emission of the point beamlets fills the grid, and a float64 mantissa of a noise-like field is incompressible, so gzip buys 1.3 times for sixteen times the read. Storing the field as float32 halves it exactly and costs 5.9e-8 of relative accuracy. That is far below anything a picture or a moment needs and far above single-precision tracking's own floor. Neither is taken. The measurement is what a later decision is made against.

Writing frames dominates a run at a movie-grade comb, and the device readback does not. The input file is one undulator segment, 15 slices of 4096 at `ngrid` 256 on the Metal backend, with a comb of 0.1 m and 31 records:

*Measured 2026-09-09: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| | frames off | frames on |
|---|---|---|
| walk | 0.134 s | 0.582 s |
| dumps | 0.000 s | 0.402 s, 69.1% |
| stats and diag | 0.027 s | 0.025 s |
| device step | 0.001 s | 0.002 s |

The readback is about 0.10 s either way, since it happens at every comb row whether or not a frame is written. It is the unaccounted remainder in both columns. With frames on it is a sixth of the walk and the file writing is two thirds. Moving the per-slice reductions onto the device would save readback the run is not spending its time in. A series that costs too much is answered by writing fewer slices, which `run.dump_slice_first` and `dump_slice_last` do, or fewer records, which the comb does.

(perf-interlude-steps)=
## The interlude's pieces, measured

`run.interlude_ds_step` cuts an interlude element into pieces so that stats rows and frames land along it ([](validation.md#val-interlude-steps)). What that costs is the field's drift rather than the beam's tracking, since a piece drifts every slice of the window where it tracks one bunch per slice. The line is three undulators of 0.6 m with two pipes of 0.6 m between them, with 32 slices of 512 macroparticles at `ngrid` 128, `comb_ds_save = 0.02`, 12 threads and a production build, best of three.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| `interlude_ds_step` | pieces an element | walk | rows |
|---|---|---|---|
| 0, the default | 1 | 0.285 s | 2016 |
| 0.06 | 10 | 0.325 s | 2592 |
| 0.02 | 30 | 0.427 s | 3680 |

The run's own phase table at the two ends of that scan says where the time goes:

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| phase | whole | cut at 0.02 |
|---|---|---|
| field drift | 0.004 s over 2 regions | 0.101 s over 60 |
| seam interlude | 0.001 s over 2 | 0.020 s over 60 |
| stats and diag | 0.035 s | 0.055 s |
| the FEL step | 0.240 s | 0.244 s |

The FEL step is untouched, as it must be. The drift is 1.3% of the walk whole and 24% cut thirty ways. It dominates because it is one transform pair per slice per piece, where the seam's own tracking is one bunch per slice. The price is set by the window's slice count and the pieces per element. An input file that wants frames through a break pays for the field's diffraction rather than for the beam's transport. Rows cost the rest: 2016 to 3680 of them moves stats and diag by 0.020 s.

Bmad's global `bmad.sr_wakes_on` chooses the interlude path for every element of a line at once. With it on, an element carrying a wake takes the serial whole-window path, and with it off every interlude takes the slice-parallel loop. On a 103-element line with 74 wake-carrying pipes at 12 threads, that one switch decides the path of the whole lattice.

(perf-device)=
## The device, measured against the CPU

The Metal backend ([](validation.md#val-device)) is run against the CPU path on the same input files: one 9 m Aramis segment of 89 steps at `ngrid` 256, with `comb_ds_save = -1` so both sides pay stats at the element end only. The table gives walk seconds, best of three, on the M3 Max (12 performance cores for the CPU, production builds). The device-busy column is the backend's own command-buffer timestamps. The difference from its wall clock is the host side of each step, dominated by the dispatch floor.

*Measured 2026-10-24: Apple M3 Max, 12 performance cores, production build, gfortran 15.3 and clang 21.1 from conda-forge, the toolchain the recorded digits belong to.*

| case | CPU, 12 threads | device wall | device busy | ratio | 2026-09-04 |
|---|---|---|---|---|---|
| steady state, 1 x 8192 particles | 0.436 s | 0.059 s | 0.018 s | 7.3x | 7.3x |
| SASE window, 24 x 4096, slippage live | 0.758 s | 0.124 s | 0.061 s | 6.1x | 8.2x |
| SASE window, 96 x 8192 | 3.509 s | 0.383 s | 0.290 s | 9.1x | 12.4x |
| dispatch floor: 1 x 1024, ngrid 64 | 0.032 s | 0.045 s | 0.014 s | 0.7x | 0.9x |

The input files are `tests/perf/*.in` and `tests/perf/device_table.sh` prints the table. The last column is the ratio of 2026-09-04, under gfortran 13.3 and the clang of 2025-06, when the CPU walks were 0.360, 0.873, 3.92 and 0.032 s and the device walls 0.049, 0.107, 0.316 and 0.037. The device walls are longer on the current toolchain, by 20 to 48 percent, and the busy seconds with them, which is the host side and the Metal library of the newer clang rather than the shaders, which are compiled at run time by the same Metal compiler. The source filter is on in these runs and was not the default then.

This table says where a step's device seconds go, per pass. `run.device_timing` encodes each pass into its own compute encoder carrying a timestamp attachment, which is what stage-boundary counter sampling can measure on an Apple GPU. Dispatch-boundary sampling would time each dispatch inside one encoder, and the M3 Max reports it absent. The input file is the full Aramis line, 12 segments and 1068 steps, one slice, `ngrid` 256, with the filter on and `comb_ds_save = -1`. The readbacks therefore fall at element ends, and each element's steps batch into one command buffer. Each row is the best complete run of five, never a per-pass minimum, with the spread over the five beside it. The transverse map runs twice a step, and the filter and the solve are four dispatches each, so those lines sum more encoders than the others.

*Measured 2026-09-09: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| particles a slice | walk | busy | transverse | push | zero | deposit | filter | solve | the six | spread |
|---|---|---|---|---|---|---|---|---|---|---|
| 1024 | 0.953 s | 0.311 s | 0.0278 s | 0.0273 s | 0.0143 s | 0.0201 s | 0.0831 s | 0.0871 s | 0.2597 s | 12.5% |
| 8192 | 0.875 s | 0.256 s | 0.0197 s | 0.0184 s | 0.0116 s | 0.0142 s | 0.0577 s | 0.0605 s | 0.1821 s | 27.7% |
| 131072 | 1.800 s | 0.312 s | 0.0482 s | 0.0423 s | 0.0147 s | 0.0825 s | 0.0590 s | 0.0614 s | 0.3081 s | 7.1% |

The instrument is not free, and its price says how far the table above can be read. The same three loads are run with the timing off against on, best of five:

*Measured 2026-09-09: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| particles a slice | walk off | walk on | busy off | busy on | busy |
|---|---|---|---|---|---|
| 1024 | 0.611 s | 0.912 s | 0.1120 s | 0.2920 s | 2.61x |
| 8192 | 0.683 s | 0.882 s | 0.1240 s | 0.2580 s | 2.08x |
| 131072 | 1.700 s | 1.800 s | 0.2740 s | 0.3170 s | 1.16x |

The absolute seconds in the per-pass table are therefore the instrument's rather than a production run's. At the two lighter loads the encoder boundaries cost about as much as the work. The shares survive that better than the seconds do, with one bias. A pass bears one boundary per encoder, so the deposit's one encoder a step carries less of the overhead than the filter's four or the solve's four, and the deposit's share is therefore a floor. On that reading the deposit is 7.7 percent of the six at 1024 and 8192 particles a slice and 26.8 percent at 131072, where the overhead is smallest and the figure is firmest. The transform, filter and solve together are two thirds of a step at production loads. That is the conclusion the load scan of [](#perf-device) reaches from the outside: the step is dominated by the grid and not by the particles until the load is very high.

The fixed-point deposit's cost is measured on the same three loads and the same input file as the two tables above, timing off on both sides so the instrument's own boundaries are out of it, best of five. The float deposit is the build before the accumulator changed.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| particles a slice | busy, float | busy, fixed point | busy | walk, float | walk, fixed point | walk |
|---|---|---|---|---|---|---|
| 1024 | 0.1120 s | 0.1390 s | 1.24x | 0.611 s | 0.649 s | 1.06x |
| 8192 | 0.1240 s | 0.1640 s | 1.32x | 0.683 s | 0.733 s | 1.07x |
| 131072 | 0.2740 s | 0.3550 s | 1.30x | 1.700 s | 1.800 s | 1.06x |

The fixed point costs a quarter to a third more device time on those input files, and six to seven percent more wall clock. A single slice leaves the step dominated by the host's dispatch rather than by the device, so the backend is busy for only about a sixth of the walk there. That is the cheap end of the range and not a general figure, so the same comparison is made below on a window that keeps the device busy.

This table says where the device time goes after the change, on the input file and conditions of the per-pass table above, so that the passes can be compared rather than inferred.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| particles a slice | walk | busy | transverse | push | zero | deposit | filter | solve | the six | spread |
|---|---|---|---|---|---|---|---|---|---|---|
| 1024 | 0.820 s | 0.227 s | 0.0157 s | 0.0162 s | 0.0150 s | 0.0162 s | 0.0605 s | 0.0507 s | 0.1743 s | 26.4% |
| 8192 | 0.891 s | 0.236 s | 0.0188 s | 0.0176 s | 0.0187 s | 0.0187 s | 0.0664 s | 0.0571 s | 0.1974 s | 3.8% |
| 131072 | 1.898 s | 0.356 s | 0.0452 s | 0.0399 s | 0.0230 s | 0.1347 s | 0.0692 s | 0.0593 s | 0.3713 s | 3.4% |

The two tables were taken at different times, so the passes the change does not touch say whether they can be compared at all. At 8192 and 131072 the transverse, push and solve lines sit at 0.94 to 0.97 of the earlier ones, close enough to read the rest. At 1024 they sit at 0.56 to 0.59 and its spread is 26 percent, so that row measures nothing. On the two rows that carry, the clear is 1.6x and the deposit 1.3x at 8192 and 1.6x at 131072. Every other pass is unmoved. The deposit does two atomic adds a component where it did one, the clear writes twice the words, and the pass that converts reads four integers where it read two floats.

What it costs where the device is the bottleneck is measured on one Aramis segment, 96 slices at `ngrid` 256 with shot noise and `comb_ds_save = -1`. That window leaves the backend busy for about 68 percent of the walk instead of a sixth. The two builds are the one before the accumulator changed and the one after, timing off on both, best of five after a discarded warm-up.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| slices x particles | walk, float | walk, fixed point | walk | busy, float | busy, fixed point | busy |
|---|---|---|---|---|---|---|
| 96 x 1024 | 0.329 s | 0.385 s | 1.17x | 0.226 s | 0.283 s | 1.25x |
| 96 x 8192 | 0.357 s | 0.466 s | 1.31x | 0.244 s | 0.337 s | 1.38x |

That table has the clear as it first stood, one word a thread. On this window the clear was the largest single pass at 0.0863 s of 0.328 s busy, a quarter of all device time, since it writes 96 slices of 256 by 256 accumulators every step. It is a memory fill, so what it costs is the width of a store. Four words a thread halves it, to 0.0427 s. The same widening on the float deposit is worth nothing, 0.99x and 1.01x on the walk, because that clear writes half the words and was never the bottleneck. The saving belongs to the fixed point's own bytes and is credited there.

With the clear widened on both sides and the deposit's gamma guard in the kernel, six pairs are run alternating between the two builds so that drift falls on both:

*Measured 2026-09-12: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| slices x particles | walk, float | walk, fixed point | walk | busy, float | busy, fixed point | busy |
|---|---|---|---|---|---|---|
| 96 x 1024 | 0.318 s | 0.348 s | 1.09x | 0.217 s | 0.239 s | 1.10x |
| 96 x 8192 | 0.366 s | 0.405 s | 1.11x | 0.248 s | 0.293 s | 1.18x |

The wall-clock price is six to seven percent where the host's dispatch dominates and nine to eleven percent where the device does. The share of the walk the backend holds sets which end an input file sits at. The kernel's gamma guard is inside those numbers and costs nothing the walk can resolve, one comparison and one branch a particle.

What the price buys is in [](validation.md#val-device): two runs of one input file identical on every array of the statistics file, where 54 of 127 differed before.

(perf-device-unaveraged)=
## The unaveraged mode on the device

A second polarization plane costs one more plane: every field plane transforms, and the kick gathers and deposits into one more. On the one-segment unaveraged input file at 89 record steps, forcing the second plane raises the device's busy time from 0.250 s to 0.295 s and the walk from 0.325 s to 0.372 s. That is 18% of the device and 14% of the walk. The kernels are the same kernels with a branch, so a one-plane run pays none of it.

The mode's parallelism on the CPU is over slices and nothing else, so what the device is worth depends on how many slices an input file has. Both regimes are measured, since they are different questions, and only the single-slice one shows the parallelism the CPU leaves unused. The input file is one 4 m Aramis segment in the unaveraged mode, 89 record steps of 60 substeps, `ngrid` 128, `comb_ds_save` at its default, production builds, on the M3 Max. Each entry is the best of five after a discarded warm-up, with the spread over the five beside it.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| slices x particles | CPU | device | ratio | device busy | busy share of the walk |
|---|---|---|---|---|---|
| 32 x 512, 12 threads | 8.508 s, 3.8% | 0.941 s, 13.1% | 9.0x | 0.833 s | 89% |
| 32 x 4096, 12 threads | 16.662 s, 3.0% | 1.573 s, 3.1% | 10.6x | 1.308 s | 83% |
| 1 x 2048, one core of twelve | 3.687 s, 1.5% | 0.365 s, 3.6% | 10.1x | 0.299 s | 82% |
| 1 x 16384, one core of twelve | 12.616 s, 1.0% | 0.455 s, 2.9% | 27.7x | 0.335 s | 74% |

The single-slice rows are against a CPU using one core of twelve, which is what that input file gets on the CPU. The last row is where the difference is largest: the device carries sixteen times the particles for 25% more time while the CPU pays the full factor. The many-slice rows are against twelve CPU threads and land near ten. The backend is busy for three quarters to nine tenths of the walk in every row. This mode is therefore device-bound rather than dispatch-bound, which the averaged path's single-slice input files are not.

This table says where a step's device seconds go, per pass, with `run.device_timing = true`. The input file is the 32-slice one above, under the same conditions, in one complete run. The half push runs twice a substep and the solve is four dispatches, so those lines carry more encoders than the others, and the record step's energy baseline adds one dispatch to the quiver line.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| particles a slice | busy | solve | zero | quiver | kick | spont | the five |
|---|---|---|---|---|---|---|---|
| 512 | 1.568 s | 0.7011 s | 0.1895 s | 0.1521 s | 0.1407 s | 0.1702 s | 1.354 s |
| 4096 | 1.722 s | 0.5098 s | 0.1475 s | 0.2990 s | 0.4431 s | 0.1248 s | 1.524 s |

The kernels this mode adds are the quiver, kick and spont lines: 34% of the five at 512 particles a slice and 57% at 4096, against the solve and the clear at 66% and 43%. That is the same crossover the CPU's own sampling split shows ([](#perf-the-sampling-split-unaveraged)). There the transform is 21% of the mode at 16384 particles a slice, and the mode is field-solve dominated at low load. The ledger's spontaneous reduction is 8 to 13% of the device's time, which is what banking that term costs.

The instrument's own price applies here as it does above: encoder boundaries the production path does not pay. This mode encodes eight dispatches a substep, so a segment runs past 42,000 of them and the sample buffer is committed early many times over. The seconds in that table are the instrument's. The shares are what carries.

What the fixed-point deposit costs is re-measured here. This mode deposits twenty times an undulator period where the averaged path deposits once a record step, so the question is whether the higher rate changes the answer. The comparison is the accumulator narrowed to one 32-bit atomic a component and a two-word clear, the float accumulator's own traffic, against the fixed point as it ships. The input files are the same, timing off on both, best of five.

*Measured 2026-09-10: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| slices x particles | walk, narrow | walk, fixed point | walk | busy, narrow | busy, fixed point | busy |
|---|---|---|---|---|---|---|
| 32 x 512 | 0.897 s | 0.901 s | 1.00x | 0.776 s | 0.819 s | 1.06x |
| 32 x 4096 | 1.421 s | 1.598 s | 1.12x | 1.140 s | 1.327 s | 1.16x |
| 1 x 16384 | 0.427 s | 0.488 s | 1.14x | 0.300 s | 0.353 s | 1.18x |

The price is zero to fourteen percent on the walk and six to eighteen on device time, the range the averaged path measured for the same change. The load per deposit sets the price rather than the rate. An input file with few particles a slice pays almost nothing, because it is paying for the clear rather than the atomics. What the price buys is in [](validation.md#val-device-unaveraged): two device runs of one unaveraged input file agree on their diag and ledger files bit for bit, at 4 threads against 8.

The field set is measured on the planar segment of the harmonics check (3.96 m, 88 steps), since a helical segment couples only the fundamental and the set has nothing to carry there. The window is the same 96 x 8192 at `ngrid` 256 and `comb_ds_save = -1`:

*Measured 2026-09-12: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| case | CPU, 12 threads | device wall | device busy | ratio |
|---|---|---|---|---|
| fundamental only | 3.344 s | 0.269 s | 0.185 s | 12.4x |
| harmonics 1, 3 | 5.249 s | 0.503 s | 0.365 s | 10.4x |
| two polarizations (tilt 0.3) | 4.502 s | 0.418 s | 0.296 s | 10.8x |

A second member costs the device 87% more wall clock and the CPU 57% more. A second plane costs the device 55% more and the CPU 35% more. The device pays the extra planes' transforms in full and the gather and deposit loops over members per particle. The CPU's field solve is already the smaller share, so the ratio narrows by a fifth and no more. The rows are walk time in a production build, min of two runs.

The readback schedule is the other setting. The rows above read back at the element end only. `comb_ds_save = 0` (a stats row every step) pays a beam-and-field readback per step, the same trade the reference manual describes for `output_step`. The 96 x 8192 case is measured at that finest comb, with the readback's conversion loops first serial and then parallel over slices:

*Measured 2026-09-12: Apple M3 Max, production build, gfortran 13.3 and the clang of 2025-06, kept as measured.*

| comb = 0, walk | host stats | readback | per stats row |
|---|---|---|---|
| serial readback: 1.231 s | 0.425 s | 0.409 s | ~10 ms |
| parallel readback: 0.897 s | 0.412 s | 0.181 s | ~7 ms |

For calibration, the same input file is 0.303 s at `comb_ds_save = -1` and 4.41 s on the CPU at comb 0. The finest comb therefore costs the device path a factor of three and the CPU path a tenth. The parallelization recovered the conversion arithmetic and no more, and the thread scan says why. The readback share is 0.371 s at one thread, 0.195 at four and 0.185 at twelve, a memory-bandwidth plateau. What remains is the bytes, roughly 70 MB of FP32 transfers plus their FP64 stores per record, about 2 ms on this machine, and no thread count reduces bytes. That number is the input to the device-side moments question. Computing the per-slice stats reductions on the GPU would replace the megabytes with some fifty floats per slice, the only way under the plateau, at the cost of FP32 accumulation moving the stats digits. Whether that trade is taken is decided against this table. The analysis side (Twiss, emittance, the element-end bunch parameters) stays on the CPU in every variant.
