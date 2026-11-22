---
title: "The electrons mode: one particle for every electron"
short_title: Electrons mode
---

The [load](glossary.md#gl-load) of a run ordinarily represents the electrons of a [slice](glossary.md#gl-slice) by a much smaller number of [macroparticles](glossary.md#gl-macroparticle), on which the [shot noise](glossary.md#gl-shot-noise) of the beam is imposed ([](shot-noise.md)). With `beam.load_mode = "electrons"` the load instead contains one particle for every electron of the bunch. Every particle has coordinates and a phase of its own, there are no [beamlets](glossary.md#gl-beamlet), and the shot noise is a property of the load.

The page has two parts. The first describes the mode: the load ([](#el-load)), the keys that the mode does and does not use ([](#el-keys)), the number of particles and the memory that they occupy ([](#el-size)), the run time ([](#el-time)) and the limits of the mode ([](#el-scope)). It follows from the method and from the sizes of the arrays of the program, and it contains no result of a run. The second part, [](#el-measured), records what was measured on the examples of the repository, and its numbers belong to those examples and to the machines named there.

(el-load)=
## The load

The charge $Q$ of the [window](glossary.md#gl-window) is counted as the nearest whole number of electrons, $N = \mathrm{nint}(Q/e)$, and every particle has the weight $Q/N$, which differs from $e$ by at most one part in $2N$. The loaded charge is therefore $Q$ exactly.

The numbers of electrons in the slices are one multinomial draw. The slices are taken in order, and slice $i$, whose charge $Q_i$ follows from the current profile of the bunch ([](fel-physics.md#sec-loading)), receives a binomial share of the electrons not yet assigned, with a probability equal to $Q_i$ divided by the charge not yet loaded. The numbers $N_i$ sum to $N$, and each fluctuates about its mean as the number of electrons in a slice of a real bunch does. Each electron is then assigned four transverse coordinates from the Gaussian distribution matched to the Twiss parameters of the lattice, an energy from the Gaussian distribution of width `sig_pz`, and a phase distributed uniformly over the $2\pi\,\texttt{n\_wavelength}$ of its slice. The current is taken as uniform within a slice, as it is in the other modes.

No two particles share a coordinate. For $N_i$ independent phases

$$
  \Bigl\langle \Bigl|\sum_{j=1}^{N_i} e^{-ik\theta_j}\Bigr|^2 \Bigr\rangle = N_i
$$ (eq-electrons-noise)

at every wavenumber $k$ for which $k\,\texttt{n\_wavelength}$ is a nonzero integer. The bunching of a slice consequently has the mean square $1/N_i$ at the fundamental and at every harmonic, which is the shot noise of $N_i$ electrons ([](shot-noise.md#shot-real)). No noise is imposed on the load.

Each slice takes its random numbers from a state of its own, which is seeded by `run.ran_seed` and the index of the slice. The slices are loaded in threads, and the load does not depend on the number of threads.

### A bunch from a file

With `beam.position_file` the particles of the bunch are retained at their own coordinates, as in the keep mode without the [quiet start](glossary.md#gl-quiet-start). The file is accepted only if each of its particles is one electron: the charges are equal, and their sum is as many electrons as there are particles. Any other file is refused.

(el-keys)=
## The keys of the mode

The keys of the beamlets and of the imposed noise have no meaning in this mode. A value other than the default in `beam.n_particle`, `quiet_start`, `beamlet_size`, `beamlet_min_electrons` or `shot_noise_method` is refused, and `shot_noise` is not read. The file `<out_root>.resolved.toml` of a run states every key at its default, and is accepted as the input file of the same run.

The transverse grid, the step and the window are derived as in the other modes ([](input-reference.md)).

| key | in the electrons mode |
|---|---|
| `beam.load_mode` | `"electrons"` |
| `beam.bunch_charge` | Determines the number of particles, $6.24\times10^{6}$ per pC. |
| `run.ran_seed` | Seeds the numbers of electrons in the slices and the state of every slice. |
| `run.write_final_beam` | False unless it is set, since the beam file at the exit would contain every electron. |
| `run.migrate` | Moves an electron to the neighboring slice when its phase leaves its own ([](#el-migration)). |

(el-migration)=
### Migration

A particle of this load is one electron, so that a move to the neighboring slice separates no beamlet, and the radiation that [migration](glossary.md#gl-migration) adds to a load of beamlets ([](fel-physics.md#sec-migration)) does not arise. The pass runs in threads. It leaves in every slice the particles that a serial scan of the slices in order leaves there, in the same order, so that a result does not depend on the number of threads.

(el-size)=
## The number of particles and their memory

The number of particles is the charge of the window divided by the elementary charge, $N = Q/e$, which is $6.24\times10^{6}$ per pC. It does not depend on the grid, on the number of slices or on any other key.

On the host a particle occupies 64 bytes: six coordinates and a weight in FP64, and a label of 8 bytes that follows the particle. A slice is allocated with a margin of one sixteenth for the particles that migrate into it, so that the particles of a run occupy 68 bytes per electron.

On the [device](glossary.md#gl-device), `run.device = "metal"`, a particle occupies 32 bytes: five coordinates and a weight in [FP32](glossary.md#gl-fp32), and a phase held as an 8-byte integer. The device allocates every slice at the number of particles of the largest slice, so that its memory is 32 bytes times the number of slices times that number. This equals 32 bytes per electron for a uniform current and exceeds it otherwise. The host retains its own particles during a run on the device, together with 8 bytes per particle for the phase of the previous step. The memory of the device is therefore required in addition to that of the host and does not replace it.

| charge in the window | particles | host, FP64 | device, FP32 |
|---|---|---|---|
| 0.1 pC | $6.24\times10^{5}$ | 0.0424 GB | 0.0200 GB |
| 1 pC | $6.24\times10^{6}$ | 0.424 GB | 0.200 GB |
| 10 pC | $6.24\times10^{7}$ | 4.24 GB | 2.00 GB |
| 100 pC | $6.24\times10^{8}$ | 42.4 GB | 20.0 GB |
| 1000 pC | $6.24\times10^{9}$ | 424 GB | 200 GB |

:::{figure} figures/electrons-memory.svg
:label: fig-electrons-memory
:width: 80%

The memory of the particles as a function of the charge in the window, at 68 bytes per electron on the host and at 32 bytes per electron on the device for a uniform current. The upper axis is the number of particles. The broken line is at $2^{31}-1$ particles, 344 pC, above which a beam file and an element wake are refused. `docs/figures/electrons_memory_figure.py` draws the figure and prints the table.
:::

Two further quantities of memory increase with the number of particles. In an [interlude](glossary.md#gl-interlude) each thread converts one slice at a time to a Bmad bunch of 224 bytes per particle, and retains the bunch of the largest slice that it has converted. The memory of the interludes is therefore at most 224 bytes times the number of threads times the number of particles of the largest slice. The beam file at the exit, where it is requested, is written through a Bmad bunch of 224 bytes per particle that contains the whole window.

A Bmad bunch is indexed by a default integer. A beam file, and an element wake of Bmad, which requires the whole window as one bunch ([](fel-physics.md#sec-seamwake)), are therefore refused above $2^{31}-1$ particles, which is 344 pC. The tracking itself counts particles in 8-byte integers and has no such limit.

The memory of the field, of the record of the run and of the [escaped field](glossary.md#gl-bank) is that of the other modes and does not depend on the number of particles.

(el-time)=
## The run time

The particle push and the statistics pass over every electron at every step, and every interlude and every migration pass does so once. The field solve is that of the other modes: its operation count depends on the grid and on the number of slices, and not on the number of particles. The measured division of the [walk](glossary.md#gl-walk) among these phases is given for one example in [](#el-measured).

(el-scope)=
## Scope

The [validation suite](glossary.md#gl-validation-suite) tests the mode in its `electron-load` section ([](validation.md#val-electron-load)): the number and the weights of the particles, the noise of [](#eq-electrons-noise) at the harmonics and between them, the independence of the slices, the refusals, the migration pass, and the pulse energy after one undulator as a ratio to that of the weight load.

A run is on one node, and its memory is that of [](#el-size).

With `bmad_com%radiation_fluctuations_on` the slices of an interlude are tracked one after another, so that the random numbers are taken in a fixed order.

The validation suite runs the mode on the CPU alone. No test covers it under `run.device = "metal"`.

(el-measured)=
## Measurements on the examples

The sections above describe the mode and hold for any run. This section records measurements. Every number in it was obtained on an example of the repository, on the machine and with the number of threads stated, and is revised when the example or the program changes.

### Memory and run time of `lcls_electrons`

The example [`lcls_electrons`](generated/examples/lcls_electrons.md) is the window of the example `lcls`, 30.8 pC, with 192,386,391 electrons in 300 slices over 132 m. By [](#el-size) its particles occupy 13.1 GB. The run reached 16.9 GB on 12 threads and 31.3 GB on 128, the difference being the bunches of the interludes. It is run by hand and by no test.

| phase | Apple M3 Max, 12 threads | two AMD EPYC 7763, 128 threads |
|---|---|---|
| particle push | 4330 s | 974 s |
| statistics | 1253 s | 217 s |
| field solve | 746 s | 195 s |
| interludes | 375 s | 133 s |
| walk | 6704 s | 1519 s |

The two runs gave the same exit power, pulse energy and bunching to every digit that the program prints. The conditions of the measurement, and the time of an interlude by part as a function of the number of threads, are given in [](performance.md#perf-electrons).

With `bmad_com%radiation_fluctuations_on`, where the slices of an interlude are tracked one after another, an interlude took 4.1 times as long on twelve threads.

### The migration pass

On $9.6\times10^{6}$ electrons in 96 slices, 8 passes that moved $1.9\times10^{6}$ particles took 0.054 s on twelve threads, 0.13% of the walk.

(el-compare)=
### Comparison with the weight load

The weight load of the sample mode represents the same shot noise with a number of macroparticles that does not depend on the charge ([](shot-noise.md#shot-weight)). The two loads were compared at the charge of a short window over many seeds, with a seed field, and at the full charge of a machine.

#### Start-up from noise

In the window of the example `sase`, which contains $1.8\times10^{6}$ electrons in 96 slices, the two loads were run over 48 seeds each on a 128-point grid, the weight load with 2048 macroparticles per slice. The uncertainty of each ratio is the standard error.

| position | mean pulse energy, electrons over weight load |
|---|---|
| 4 m | $1.00 \pm 0.03$ |
| 12 m | $1.00 \pm 0.05$ |
| 20 m | $0.95 \pm 0.06$ |
| 30 m | $0.87 \pm 0.09$ |
| 57 m, the exit | $0.83 \pm 0.11$ |

The pulse energy of a single seed at the exit has a standard deviation of 0.6 to 0.7 of its mean in either load. The first 16 of the seeds gave $0.63 \pm 0.13$ at the exit, and the other 32 gave $0.95 \pm 0.16$. The start-up level is therefore the same to between 3% and 5%, and the ratio at the exit is known to 11%.

#### Gain with a seed

The gain was compared without that scatter by seeding the same window with 100 MW, a power beside which the noise is small. The pulse energy of the electrons load is then within 3% of that of the weight load along the whole line and within 1.3% at the exit. It is within 2% of that of a load with no noise, at 2048 and at 32768 macroparticles per slice.

#### The full charge of LCLS

The example `lcls_electrons` is the window of the example `lcls` with its $1.9\times10^{8}$ electrons. The example `lcls` represents the same window with 2048 macroparticles per slice. Each run is one seed.

| quantity | electrons | weight load |
|---|---|---|
| particles | 192,386,391 | 614,400 |
| gain length | 3.537 m | 3.549 m |
| saturation point | 65.02 m | 65.20 m |
| largest power of a slice at saturation | 36.9 GW | 37.0 GW |
| pulse energy at the exit | 278 µJ | 285 µJ |
| power outside the split angle divided by the power inside it | $2.3\times10^{-5}$ | $6.1\times10^{-4}$ |

The gain length, the saturation point and the power at saturation agree to better than 1%, and the pulse energy at the exit to 2.5%. The power radiated at wide angles is 27 times smaller for the electrons, which are 313 times as numerous in a slice as the macroparticles of the weight load.
