---
title: The spectral field propagator
short_title: Spectral propagator
---

The spectral field propagator, `field.propagator = "spectral"`, transforms the whole window in $(x, y, \tau)$, so each frequency of the field diffracts at its own wavenumber and the slippage is a continuous phase ramp with no rotation and no threshold. The slice propagator stays the default and every tier of the keystone takes it. This page is what the spectral propagator is, what it is measured to do against the slice propagator, Genesis4 and an exact solution, and what it costs. Its check in the keystone is the `spectral` section of the validation suite ([](validation.md#val-spectral)). In one paragraph: the steady state is reproduced to rounding, the SASE tier with the source filter on agrees with the slice propagator to 8.7e-3 of the exit power, and the difference is the slippage's form and not the kernel's, since a seeded chirped pulse through Aramis shows the same agreement at every bandwidth from 0.5 to 5 percent. Where the slice propagator is wrong in closed form, a chirped pulse diffracting in vacuum, the spectral propagator is right to 1e-9. What the spectral propagator is for is a field sampled coarser than the beam's slices along $\tau$, which only a continuous slippage can carry, and which is the next step.

(spectral-what-changes)=
## What changes and what does not

The slice propagator ([](fel-physics.md#sec-field), [](fel-physics.md#sec-slippage)) holds the field as a stack of transverse grids, one a slice, and advances every slice over a step by the paraxial kernel $\exp(-i\,\delta z\, k_\perp^2 / 2k_0)$ at the one reference wavenumber $k_0 = 2\pi/\lambda_s$. The slippage rotates the stack by whole slices once $0.8$ of a slice has accumulated, and the light that leaves the head is banked.

The spectral propagator, `src/fel_spectral_mod.f90`, transforms a member's whole window in $(x, y, \tau)$, where $\tau$ is the window coordinate, $\tau = (i - 1)\,\Delta$ for window slice $i$ and slice spacing $\Delta$, so $\tau$ grows toward the head. Two things change.

- **The kernel is frequency dependent.** The envelope rides the carrier $\exp(i(k_0 z - \omega_0 t))$, the convention the kernel's sign fixes, so a mode $\exp(i\kappa\tau)$ of the envelope has the carrier's frequency plus $c\kappa$ and the wavenumber $k(\kappa) = k_0 + \kappa$. Its kernel is $\exp(-i\,\delta z\, k_\perp^2 / 2k(\kappa))$. A pulse that is short or chirped carries a spread of $\kappa$, and one $k_0$ for all of it is wrong at relative order $\delta\omega/\omega$.
- **The slippage is a phase ramp.** A shift of the field toward the head by $\delta$ is $\exp(-i\kappa\delta)$ in the spectral domain, exact and continuous. A step may slip a fraction of a slice, and the field is never a fraction of a slice behind the beam, where the rotation leaves it up to a slice behind between crossings of its threshold.

The push, the deposit, the source filter and the averaged equations are untouched. A step runs in the slice propagator's order: the kernel over the step, then the deposit of every slice's source into the record, then the slippage, which moves the field and the fresh source together. The deposit is `fel_spectral_deposit`, the per-particle scatter of `fel_field_step` kept apart from it so the slice propagator stays as it was, with the coherent source model refused.

The kernel factors: $\exp(-i\,\delta z\,(m_x^2 + m_y^2)\,\delta k^2 / 2k(\kappa))$ is the product of a table in $(m_x, \kappa)$ and one in $(m_y, \kappa)$, each `grid_n_pts` by the transform length, so the whole kernel costs one complex multiply a point and the tables are rebuilt at every call. The transverse transforms are the per-slice ones of `wavefront_fft2`, parallel over slices, and the transform along $\tau$ runs on each line of the transverse spectrum, parallel over the second index, through a plan cache in the same form, one plan pair and one buffer a thread.

In the steady state the window is one slice, the transform along $\tau$ is the identity, and the two propagators agree to rounding: on the `tier1` tier the power curves differ by $4.8\times10^{-15}$.

(spectral-keys)=
## The keys

| Key | Default | What it does |
|---|---|---|
| `field.propagator` | `slices` | `spectral` selects the spectral propagator. |
| `field.n_wavelength` | 0 | The field's plane along $\tau$ in wavelengths, a whole multiple of `window.n_wavelength`. 0 takes the window's slice under `slices` and derives a coarser plane under `spectral` ([](#spectral-coarse)). |

The spectral propagator covers the averaged mode on the CPU in double precision, for the fundamental alone when the field is coarser than the slices. The unaveraged mode, the device, the FP32 check and the coherent source model each carry a propagator or a deposit of their own, and each is refused with the key, the message naming which. The header prints a `Propagator` line that states the kernel, the ramp and the pad in force, and a `Field` line with the planes when they are coarser than the slices.

(spectral-coarse)=
## The field's own sampling along tau

The window is the beam's partition, `window.n_wavelength` wavelengths a slice, and under the slice propagator the field has one plane a slice and no other sampling. Under the spectral propagator the field may be sampled coarser, `field.n_wavelength` wavelengths a plane, a whole multiple of the slice ([](fel-physics.md#sec-window)). A SASE field at saturation is band limited to a few $\rho$ about the carrier, so a plane every fraction of the cooperation length $\lambda_s / 4\pi\rho$ carries it, while the beam's slice stays where the loader put it, and a step then pays for the field at that sampling and for the beam at its slice. The continuous slippage is what makes it possible: a step slips a fraction of a plane, which the ramp applies exactly and the whole-slice rotation could not.

The field's record holds the window at that sampling, `fel_field_planes` planes of it, with one more past the last slice when the last slice falls between two planes. A beam slice sits between two planes and reads them by linear weights along $\tau$ (`fel_field_sample`), $1 - w$ on the plane at or below it and $w$ on the one above, so a slice on a plane reads that plane alone, which is every slice at one plane a slice, where nothing changes to the bit. The deposit splits each slice's source between the same two planes by the same weights and divides a plane's sum by the slices a plane spans, so the deposit is the readback's adjoint under each grid's own measure, a uniform source over the slices gives every plane that one slice's source, and the source filter runs once on each plane's sum (`fel_spectral_deposit_coarse`), one transform pair a plane a step and not one a slice. The stats' field rows stay on the beam's axis, each slice's row taken from the plane it reads, so every reader and comparison works as before, and the field's own planes are what the dumps and the pulse file carry, their spacing the field's. A frame of a slice range writes the planes those slices read. The bank banks a plane each time the slippage crosses a whole plane, and the header's `Bank` line counts planes. A field file at another spacing is refused unless it is dark, when its transverse grid alone is taken and the field starts dark on its own planes, so a SASE input file whose field file fixes the grid runs at any sampling. The escaped light's planes are the field's spacing apart in the escaped and pulse files, and `lucifer.field` reads that spacing as it reads the slice's.

The sampling is derived at setup when `field.n_wavelength` is unset: a tenth of the cooperation length at the $\rho$ of the beam as built (`fel_pierce_of_beam`: the slice with the most charge, its rms size and the first undulator's coupling), rounded to whole slices, printed with its origin, and cut so the window keeps 8 planes. A stated value is held to a whole multiple of the slice, to at least 8 planes over the window, and below the cliff the sweep found, a quarter of the cooperation length, where the run is refused with the cooperation length and the cliff in the message. Harmonic fields and `run.dump_reduced` are refused with a field coarser than the slices.

(spectral-sweep)=
### The sweep

Four lines from noise under the spectral propagator, each at 1, 2, 5, 10, 20 and 40 slices a plane, against its own run at one plane a slice: the SASE tier's beam over 1280 slices of 3 wavelengths through Aramis, which is the seeded case's line and bunch with no seed, and the `sase`, `flash1` and `lcls` examples as they ship, 96 slices of 3 wavelengths, 351 of 12 and 300 of 68. Each row gives the field's planes, the window's and the pulse's energy at the exit, the peak slice power, the exit power per beam slice against the finest run over its peak, and the walk, at 12 threads on the release build, gfortran 15.3, Apple M3 Max, 2026-11-03, with the cliff open. The plane's length is also given as the fraction of the cooperation length $\lambda_s / 4\pi\rho$ at the beam's $\rho$: 155 wavelengths at 5.1e-4 on Aramis, 40.5 at 2.0e-3 on FLASH1 and 135 at 5.9e-4 on LCLS. The spectrum's rms width was also taken and is not a measure here: at one plane a slice the window's Nyquist band is a sixth of $k_0$ and the shot noise's floor across it sets the width, which then falls with the band as the plane coarsens, while the spectrum's centroid holds at 1e-3 below the carrier on every Aramis row.

| Line | Slices a plane | Plane [wavelengths] | Of $L_c$ | Planes | Window energy | Pulse energy | Peak power | Exit power against the finest | Walk |
|---|---|---|---|---|---|---|---|---|---|
| Aramis, 1280 slices | 1 | 3 | 0.019 | 1280 | 2.776 µJ | 4.151 µJ | 8.14 GW | 0 | 692 s |
| | 2 | 6 | 0.039 | 641 | 2.762 µJ | 4.128 µJ | 8.09 GW | 1.5e-2 | 500 s |
| | 5 | 15 | 0.097 | 257 | 2.762 µJ | 4.095 µJ | 8.04 GW | 3.6e-2 | 195 s |
| | 10 | 30 | 0.19 | 129 | 2.714 µJ | 4.018 µJ | 7.78 GW | 7.0e-2 | 169 s |
| | 20 | 60 | 0.39 | 65 | 2.477 µJ | 3.691 µJ | 6.45 GW | 2.6e-1 | 121 s |
| | 40 | 120 | 0.77 | 33 | 1.905 µJ | 2.889 µJ | 4.92 GW | 5.1e-1 | 105 s |
| `sase`, 96 slices | 1 | 3 | 0.019 | 96 | 34.2 pJ | 851 pJ | 1.95 MW | 0 | 97 s |
| | 2 | 6 | 0.039 | 49 | 31.2 pJ | 803 pJ | 1.75 MW | 1.2e-1 | 50 s |
| | 5 | 15 | 0.097 | 20 | 36.6 pJ | 937 pJ | 2.19 MW | 1.5e-1 | 28 s |
| | 10 | 30 | 0.19 | 11 | 30.5 pJ | 724 pJ | 1.71 MW | 1.2e-1 | 24 s |
| | 20, 40 | | | | refused: fewer than 8 planes | | | | |
| `flash1`, 351 slices | 1 | 12 | 0.30 | 351 | 47.2 µJ | 47.2 µJ | 4.59 GW | 0 | 55 s |
| | 2 | 24 | 0.59 | 176 | 42.9 µJ | 42.9 µJ | 4.21 GW | 1.8e-1 | 26 s |
| | 5 | 60 | 1.5 | 71 | 33.4 µJ | 33.4 µJ | 2.67 GW | 6.1e-1 | 18 s |
| | 10 | 120 | 3.0 | 36 | 28.2 µJ | 28.2 µJ | 1.02 GW | 9.2e-1 | 16 s |
| | 20 | 240 | 5.9 | 19 | 24.4 µJ | 24.4 µJ | 1.10 GW | 9.8e-1 | 16 s |
| `lcls`, 300 slices | 1 | 68 | 0.50 | 300 | 288 µJ | 315 µJ | 121 GW | 0 | 118 s |
| | 2 | 136 | 1.0 | 151 | 155 µJ | 171 µJ | 79 GW | 9.4e-1 | 129 s |
| | 5 | 340 | 2.5 | 61 | 84 µJ | 91 µJ | 16 GW | 9.2e-1 | 55 s |
| | 10 | 680 | 5.0 | 31 | 78 µJ | 79 µJ | 16 GW | 9.4e-1 | 53 s |

The long Aramis window is the measurement, since its 1280 slices hold enough spikes for the energy to be a mean and its slices are fine enough to place the plane from a fiftieth of the cooperation length to most of it. There the energy holds to 0.5 percent up to a tenth of the cooperation length, falls 2 percent at a fifth, 11 percent at 0.4 and 31 percent at 0.8, and the exit power per slice differs by 1.5, 3.6 and 7 percent of the peak at the first three samplings, the saturated light's own sensitivity to any change of the operations, which the step's halving showed at 9.5 percent on the seeded line ([](#spectral-seeded)). The walk falls from 692 s to 195 s at a tenth of the cooperation length, 3.5 times, and the push then holds what is left. The other three lines agree with that curve where they reach it: FLASH1 loses 9 percent at 0.6 of the cooperation length and LCLS 46 percent at 1.0, and their own slices, 12 and 68 wavelengths, already sit at 0.3 and 0.5 of it, so a coarser field than their slice costs output and the derived sampling keeps their slice. The `sase` example's window of 96 slices holds two cooperation lengths, so its energy moves 10 percent either way between samplings, a realization's own spread and not a trend, which is why a short window is no measure of the sampling.

From that curve the derived sampling is a tenth of the cooperation length, `fel_field_coop_fraction$`, where the energy holds to 0.5 percent and the field costs a fifth of its fine cost, and the cliff is a quarter, `fel_field_coop_cliff$`, past which the energy falls by percent and a stated `field.n_wavelength` is refused. On the SASE tier's window of 32 slices the derived sampling is cut to 4 slices a plane so the window keeps 9 planes, and the exit power per slice there differs from one plane a slice by 8.0e-2 at saturation, which the `spectral` section of the validation suite holds below 1.5e-1 ([](validation.md#val-spectral)). The cost estimate counts the field's planes and not the beam's slices, so the header's estimate falls with the sampling ([](performance.md)).

(spectral-bank)=
## The head and the bank

A transform along $\tau$ is periodic, so light the ramp moves past the head would reappear at the tail. The window is followed along $\tau$ by a pad of 4 free slots and a $\cos^2$ sponge of 12 slots that falls to zero at the pad's end, and every transform runs over the window and the pad together. Light the ramp carries past the head enters the pad, diffracts there under the kernel as escaped light does in vacuum, and is attenuated by the sponge after every slip pass, so the field is smooth at the periodic seam and nothing reaches the tail.

The escaped-field bank ([](fel-physics.md#sec-slippage)) takes its planes from the pad's first slot. The slippage accumulates in the ramp's own count, and each time it crosses a whole slice the pass is split there in the transform domain, the field at that slot is read by a direct sum over the modes, and the plane goes to the bank as the rotation's transmitted slice does, its energy added to the escaped energy. Between two crossings the ramp moves the field by exactly one slot, which the transform does exactly, so the bank receives the window's samples one slice apart, each once, and the escaped and pulse files are assembled as before, with the same readers. The sponge removes light the bank already holds and counts nothing. One pass slips at most one slice, so light walks through every slot of the sponge and reaches the seam below 1e-5 of its amplitude, and an interlude's slippage is taken in such passes. A pass of several slices was measured to jump light into the sponge's last slots, where the next pass's interpolation carried 2e-3 of the exit peak power round to the tail slice of the SASE tier's window. A dump carries the window and not the pad, and a continuation starts with an empty pad: every plane that had entered it was banked when it crossed. Backward slippage is refused with the key, as it is with the bank.

Three other forms were tried on the way and are retired. The headroom alone, with the run refused once the head slice held more than a stated share of the window's energy, refused the SASE tier's flat window at the first step, and with the refusal off it let the light pass through the beam again and again until the window held 475 times the slice run's energy. A $\cos^2$ sponge on the window's own head slices, applied before each shift, held the head at 4.5e-4 of the exit power but removed the light instead of banking it, so the pulse file held none of it. Zero padding beyond the largest shift one pass applies, with what crossed into it discarded, banked nothing either and rang where the discarded light's edge met the window. The whole-slice rotation under the spectral kernel is retired with them, since the ramp is what the spectral propagator is for.

(spectral-measurements)=
## The measurements

```
python tests/scripts/spectral_field.py --exe build-release/bin/lucifer \
        --work-dir <dir> --validation-dir <a completed validation pass's work directory>
```

The script runs the first four measurements below by default and prints them, and `seeded` and `chirped` name the two FEL cases, which take hours at 12 threads. The vacuum test needs the tracker alone. The others take the tier inputs, the imported dumps and the Genesis4 references from a validation pass's work directory, which they read and never write. Every number on this page was measured on the release build, gfortran 15.3 from conda-forge, on an Apple M3 Max with 12 cores, and the thread count of each table is stated with it. The vacuum test, the null test, the bank and the cost were measured on 2026-11-03 on the propagator as this page describes it. The two FEL cases were measured the same day on the propagator as it first stood, with zero padding at the head in place of the bank, and are kept as measured: the padding and the bank differ at the head of a window alone, and the seeded window has the line's slippage of headroom.

(spectral-vacuum)=
### The vacuum test

A Gaussian pulse with a linear chirp, $\exp(-r^2/w_0^2)\exp(-\xi^2/4\sigma^2 + i b\,\xi^2)$ with $\xi$ the distance from the pulse center, drifts through a pipe of two Rayleigh ranges with no beam, under both propagators, against the exact paraxial solution. Each frequency of the pulse is a Gaussian beam with its own Rayleigh range $k w_0^2/2$, and the exact solution superposes them over the pulse's analytic spectrum, on a wavenumber grid eight times finer than the window carries, read at the slices. The wavelength is 1 nm, the window 1024 slices of one wavelength, the pulse's rms length 50 slices, the waist 20 µm, the grid 128 points of 3.15 µm over a 400 µm box, and the pipe 2.513 m. The error is the largest difference over the window, divided by the exact solution's largest amplitude.

| Relative rms bandwidth | Chirp $b$ [1/m²] | Slice propagator | Spectral propagator |
|---|---|---|---|
| 1.6e-3 (no chirp) | 0 | 1.22e-3 | 1.10e-9 |
| 5.0e-3 | 2.98e14 | 3.83e-3 | 1.12e-9 |
| 1.0e-2 | 6.20e14 | 7.67e-3 | 1.21e-9 |
| 2.0e-2 | 1.25e15 | 1.54e-2 | 1.65e-9 |
| 5.0e-2 | 3.14e15 | 3.88e-2 | 1.28e-8 |

The slice propagator's error is 0.77 of the relative bandwidth at every row, which is the spatio-temporal coupling one $k_0$ cannot carry: over two Rayleigh ranges the red and the blue ends of the pulse diffract to different sizes and phases. The spectral propagator sits at the rounding of its transform pairs. This is the one place the slice propagator is wrong in closed form, and it is wrong at the bandwidth, not at the bandwidth squared. A SASE pulse at saturation carries a bandwidth of order $\rho$, so the effect on the tiers is at that order and no larger, which the null test below measures.

(spectral-null)=
### The null test

The steady-state tier and the SASE tier of the validation suite ([](validation.md#val-validation-from-one-command)) rerun under `spectral`, against Genesis4 at the tier tolerances as `compare_fel.py` compares them, and against the validation pass's own slice run through `lucifer.compare`. The SASE tier is the full line from a dark start, 32 slices of 3 wavelengths, 1068 steps of 3 periods, and the line slips 3552 wavelengths, 1184 slices, so every undulator step slips 0.9963 of a slice and the interludes whole wavelengths. The tiers run with the source filter off, since Genesis4 is compared with it off. 8 threads.

| Tier | Propagator | Against Genesis4 | Tolerance | Against the slice run: exit power, bunching, energy |
|---|---|---|---|---|
| `tier1` | slices, the recorded level | 1.825899e-06 | 1e-4 | |
| `tier1` | spectral | 1.825899e-06 | 1e-4 | 4.8e-15, 2.3e-14, 0 |
| `tdsase` | slices, the recorded level | 2.349843e-06 | 1e-3 | |
| `tdsase` | spectral | 8.26e+06 | 1e-3 | 0.124, 0.212, 0.516 |

The steady state is reproduced to rounding, as [](#spectral-what-changes) says it must be. The SASE tier as it runs is not reproduced within its tolerance: the spectral propagator differs from the slice run by 12 percent of the exit power, and the number against Genesis4 is the elementwise power at records where the slice run holds no light yet, floored at 1e-9 of the peak. Two runs on the propagator as it first stood placed the difference. With the kernel's frequency dependence removed by a temporary edit, $k(\kappa) = k_0$ for every mode, the spectral path reproduced the slice run to 7e-10 in power over the 1068 steps, so the transform along $\tau$, its order with the deposit and the slippage are right, and the whole difference is the kernel's $\kappa$ dependence. The script does not carry that edit. With the source filter on, the production default ([](fel-physics.md#sec-source-filter)), both propagators rerun against the slice run with it on:

| Propagator, source filter on | Exit power | Bunching | Energy |
|---|---|---|---|
| spectral | 8.7e-3 | 2.1e-3 | 6.0e-3 |

The filter removes the wide-angle emission of the point beamlets, and with it 93 percent of the difference. The per-particle deposit is white in $\kappa$ up to the slice Nyquist, a sixth of $k_0$ at 3 wavelengths a slice, and white in $k_\perp$ up to the grid's. The frequency-dependent kernel gives a ray at angle $\theta$ the group delay $\delta z\,\theta^2/2$ a step, 0.2 wavelengths a step at the grid's corner, which moves the wide-angle noise backward through the window where the slice propagator keeps every slice's field in its slice. The gain band's own modes sit at angles near 1e-6 rad and within $\rho$ of $k_0$, where the delay is 2e-4 wavelengths a step and the kernel differs at 5e-4, and on the propagator as it first stood the whole-slice rotation under the spectral kernel showed that 1.6e-3. The one-dimensional $\rho$ of the tier is 5.4e-4 (3000 A, $K = 1.2$, 17 µm rms size, 15 mm period). The rest, 7e-3 of the exit power, is the slippage's form, continuous against a whole slice at a time, 13 times $\rho$, which grows toward the head of the window where the most light has slipped in: from 5e-4 of the peak on the tail's first lit slices to 8.7e-3 on the head's. The chirped seed below places that form's difference at 2.6e-3 over one segment with the seed the only light.

(spectral-bank-measured)=
### The bank

The SASE tier's window is a flat bunch of 32 slices under a line that slips 1184, so light reaches the head from the first step on and every slice, the head among them, holds a thirty-second of the noise at once. Both propagators ran on it with the escaped field kept and the source filter off, as the tier runs, at 8 threads. The table gives the exit energies the footer states, the escaped energy and the window's energy at the exit, their sum, what the pulse file holds, and the head and tail slices' shares of the exit power.

| Propagator | Window energy at the exit | Escaped energy | Sum | Pulse file | Head slice | Tail slice |
|---|---|---|---|---|---|---|
| slices | 1.166e-10 J | 8.026e-9 J | 8.142e-9 J | 1216 planes, 8.142e-9 J | 9.6e-2 | 0 |
| spectral | 1.354e-10 J | 8.700e-9 J | 8.835e-9 J | 1216 planes, 8.835e-9 J | 9.3e-2 | 2e-33 |

Both bank the 1216 planes the header states, the pulse file holds the sum to 1e-5, and the tail stays dark under the ramp. The spectral run's pulse carries 8.5 percent more energy than the slice run's, the kernel's difference of [](#spectral-null) with the filter off. The pulse check of the validation suite closes the same ledger on a Gaussian pulse through a pipe to 7e-11 under the rotation and 1.2e-6 under the ramp ([](validation.md#val-spectral)): the frequency-dependent kernel moves light along $\tau$ across the head on its own, by the group delay of the wide angles, and the bank samples the ramp's crossings alone.

The beam's energy change over the walk, each slice's live charge times its mean energy change, is 1e-10 J of either sign under the two propagators, 80 times smaller than the field's energy. The averaged equations charge the beam for the field at the particles and not for the incoherent emission the deposit adds, so on a run from noise the field's energy is not the beam's loss and no closure can be read off an averaged SASE run. The unaveraged mode's ledger is the closure statement ([](fel-physics.md#sec-unaveraged)), and the spectral propagator does not reach that mode.

(spectral-cost)=
### The cost

The field phase of the SASE tier's walk, the `field solve` row of the footer's timing table, under both propagators, so the row holds the deposit, the kernel pass over the window and its pad with transforms of 48 along $\tau$, and the slip passes. The slice propagator's row holds the deposit and its transform pair. 255 by 255 grid, 32 slices and a pad of 16, 1068 steps, one run each on an otherwise idle machine, 2026-11-03.

| Threads | Propagator | Field solve | A step | Walk |
|---|---|---|---|---|
| 1 | slices | 54.6 s | 51.1 ms | 77.2 s |
| 1 | spectral | 120.2 s | 112.5 ms | 146.9 s |
| 8 | slices | 7.46 s | 6.99 ms | 11.0 s |
| 8 | spectral | 17.8 s | 16.6 ms | 21.8 s |

The spectral propagator costs 2.2 times the slice propagator's field phase at one thread and 2.4 at eight. Three things are added on this window. The transverse transform pair runs on the pad's 16 planes beside the window's 32, half again the count, which on a window of a thousand slices is 2 percent. Each step takes two passes along $\tau$, the kernel's and the slip's, each a gather of every line of the window and pad across the stack's stride, a transform pair of 48 and a scatter, with the sponge and a plane's direct sum in the slip pass. And an interlude's slippage of up to 17 slices takes a pass a slice, 400 passes over the line's 24 interludes beside the steps' 2136. The header's cost estimate carries the passes along $\tau$ as one term fitted to this table through `fel_cost_mod` ([](performance.md)).

(spectral-seeded)=
### A chirped seed through Aramis

The vacuum test's pulse with gain. The seed is the chirped Gaussian at 1 Å, rms length 15 wavelengths, waist 30 µm, 1 MW at its peak, ten times the start-up noise light a slice holds at the segment's exit, on the tier's grid, 30 slices from the tail of a window of 1280 slices of 3 wavelengths, which holds the line's slippage of 1184 slices so the seed stays with the beam to the exit. The beam is the tier's, 3000 A flat over the window at 1024 particles a slice, dumped by a run that tracks a marker alone and loaded with the seed. Relative rms bandwidths 0.5, 1, 2 and 5 percent, where 0.53 percent is the unchirped pulse's own, through the first segment and through the full line, under both propagators with the source filter on, at 12 threads. The segment also runs on the same bunch loaded without its shot noise, so the seed is the only light.

The slice propagator leaves its field a fraction of a slice behind the ramp's: the segment slips 88.67 slices, of which it rotates 88 and records the rest as the residual its field file carries, 2.0003 wavelengths, and the line slips 1184.00. The comparison moves the slice run's field by that residual along $\tau$ through the transform before it differences, since a chirped field two thirds of a slice out of step differs by its phase slope times the offset and by nothing physical. The field difference is the largest over the lit slices and four on each side, over the spectral run's peak, and the raw column is the whole window as the runs left it. Each slice's wavelength comes from the phase advance to its neighbor on the slices above a quarter of the peak power, the seed's core.

| Bandwidth | Light | Offset [slices] | Field, raw | Field, aligned | Wavelength | Peak power, slices / spectral |
|---|---|---|---|---|---|---|
| 0.53 % | seed alone, segment | 0.667 | 5.7e-2 | 2.9e-3 | 1.9e-5 | 0.994 / 0.994 MW |
| 1 % | seed alone, segment | 0.667 | 1.1e-1 | 2.3e-3 | 2.0e-5 | 0.979 / 0.979 MW |
| 2 % | seed alone, segment | 0.667 | 2.1e-1 | 3.0e-3 | 2.1e-5 | 0.987 / 0.987 MW |
| 5 % | seed alone, segment | 0.667 | 5.2e-1 | 6.5e-3 | 8.5e-5 | 0.992 / 0.993 MW |
| 0.53 % | seed and noise, segment | 0.667 | 1.2e-1 | 8.1e-2 | 1.6e-3 | 1.010 / 1.012 MW |
| 1 % | seed and noise, segment | 0.667 | 1.7e-1 | 8.2e-2 | 1.5e-3 | 0.996 / 0.999 MW |
| 2 % | seed and noise, segment | 0.667 | 2.8e-1 | 8.3e-2 | 2.1e-3 | 1.003 / 1.006 MW |
| 5 % | seed and noise, segment | 0.667 | 5.1e-1 | 6.5e-2 | 2.7e-3 | 1.006 / 1.008 MW |
| 0.53 % | seed and noise, line | 0.001 | 1.9e-2 | 1.9e-2 | 4.7e-4 | 9.58 / 9.54 GW |
| 1 % | seed and noise, line | 0.001 | 2.0e-2 | 2.0e-2 | 3.9e-4 | 9.70 / 9.67 GW |
| 2 % | seed and noise, line | 0.001 | 2.1e-2 | 1.8e-2 | 4.0e-4 | 9.32 / 9.26 GW |
| 5 % | seed and noise, line | 0.001 | 1.7e-2 | 1.6e-2 | 4.1e-4 | 8.91 / 8.84 GW |

With the seed alone, the propagators differ by 2.3e-3 to 3.0e-3 of the peak field at the three lower bandwidths and 6.5e-3 at 5 percent, and the peak powers agree to 1e-3. The floor is the slippage's form and not the kernel's: the spectral kernel with the whole-slice rotation differed from the slice propagator by 6.5e-4 at 0.53 percent and from the ramp by 2.6e-3, so what the beam sees over the segment is the field's position against it, up to a slice behind the ramp's between rotations. The part that grows with the bandwidth, 3.9e-3 above that floor at 5 percent, is the vacuum test's error scaled to this pulse, $0.77 \times 0.05 \times L / 2 z_R$ with the waist's Rayleigh range 28.3 m and the segment 3.99 m, 2.7e-3. With the beam's shot noise on, every segment row sits at 6.5e-2 to 8.3e-2 and the seed's own error is under it: the noise-only slices of the same runs differ by 8e-2 of the seed's peak and by 13 percent of their own rms, since the two kernels diffract the start-up light's wide angles differently before the filter and the gain narrow it. Through the line, where the seed has swept the window and the light is at saturation, the propagators differ by 1.6e-2 to 2.0e-2 of the peak field, 4e-4 in wavelength and 0.4 to 0.7 percent in peak power, with no trend in the bandwidth.

The reference's convergence: the 2 percent seed with the undulator step cut from three periods to one. Through the segment, the spectral run moves by 3.4e-2 of the peak field in the seed's slices, 5.1e-4 in wavelength and 1.9 percent in peak power, and its noise-only slices by 6 percent of their rms. Through the line it moves by 9.5e-2 in the lit slices, 7.9e-4 in wavelength and 2.4 percent in peak power, and by 0.80 of the peak in the tail slices, which hold 0.3 percent of the energy and are past saturation from 35 m on. The step therefore moves a seeded Aramis run more than the propagator does at every bandwidth, by 5 times through the line, and the 1e-3 the vacuum test's reference holds to is not available with a beam in the window.

A 2 percent seed through the full line is not wrong under the slice propagator at a level a measurement could see: the two propagators' exit fields differ by 1.8e-2 of the peak, their peak powers by 0.6 percent and their wavelengths by 4e-4, with no trend in the bandwidth from 0.53 to 5 percent, and halving the step moves the spectral run itself by 9.5e-2 and 2.4 percent. The gain guiding holds the mode at the beam's size through the line, so the vacuum test's growth over the Rayleigh range does not accumulate, and the bandwidth error the vacuum test isolates is below the floor the whole-slice slippage sets against the beam, 2.6e-3 over a segment.

(spectral-chirped)=
### An energy-chirped beam from noise

The same beam with a linear energy chirp along the window, 0.25, 0.5 and 1 percent head to tail, the head above the tail, and the unchirped bunch beside them, from noise through the full line under both propagators with the source filter on, at 12 threads. The bunch is the dump's with every momentum component and the momentum deviation record scaled by the same factor, so the angles stay, loaded as the dump is, one patch a slice. A chirp $\epsilon$ in energy chirps the resonant wavelength by $-2\epsilon$, the head shorter, and the slope of a line fitted to the exit wavelengths over the slices above a thousandth of the peak power, as a fraction of the mean over the window, is placed against it. The window is 1.28 fs, so 1 percent across it moves the resonance by 1.5 $\rho$ over a cooperation length of 53 slices.

| Chirp | Field | Wavelength | Lit slices | Peak power, slices / spectral | Exit energy | Slope, slices / spectral | $-2\epsilon$ |
|---|---|---|---|---|---|---|---|
| 0 | 2.1e-2 | 4.5e-4 | 514 | 8.22 / 8.13 GW | 2.77 µJ | +0.0020 / +0.0020 | 0 |
| 0.25 % | 2.5e-2 | 1.9e-4 | 274 | 43.2 / 43.3 GW | 7.98 µJ | -0.0011 / -0.0011 | -0.0050 |
| 0.5 % | 3.9e-2 | 7.5e-4 | 528 | 2.34 / 2.32 GW | 0.77 µJ | -0.0068 / -0.0067 | -0.0100 |
| 1 % | 6.7e-2 | 1.6e-3 | 368 | 0.103 / 0.101 GW | 0.028 µJ | -0.0147 / -0.0147 | -0.0200 |

The chirp decides the output and the propagator does not. At 0.25 percent the light slipping forward meets electrons of rising energy, which is what a taper does after saturation, and the exit energy is 2.9 times the unchirped run's with a single spike of 43 GW. At 0.5 and 1 percent the resonance moves across a cooperation length by more than the gain bandwidth and the exit energy falls to 0.28 and 0.010 of the unchirped. The exit wavelengths slope at 0.2, 0.7 and 0.7 of the resonance's, above the +0.002 the unchirped SASE carries. The two propagators agree on every one of those numbers: the peak powers to 0.3 to 2 percent, the slopes to 1e-4 and the slice wavelengths to 2e-4 to 1.6e-3. The field difference grows with the chirp, from 2.1e-2 to 6.7e-2 of the peak, which is the saturated light's own sensitivity: the unchirped run's 2.1e-2 is the null test's filter-on level, and halving the step moves a seeded run by 9.5e-2.

(spectral-next)=
## What is next

The measurements say the kernel's frequency dependence does not reach the output of an Aramis run through the gain, and that the one difference the beam sees is the slippage's form: continuous under the ramp, a whole slice at a time under the rotation. That is what lets the field be sampled coarser than the beam's slices along $\tau$, at the cooperation length instead of at the slice, which is the spectral propagator's purpose and the next step. The spectral propagator keeps the envelope and its carrier, so the SVEA stands and the bandwidth the window carries is bounded by the slice spacing at $\lambda_s / 2\Delta$. A field without an envelope, Puffin's model, is a different propagator and is not this one.
