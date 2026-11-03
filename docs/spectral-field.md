---
title: The spectral field propagator
short_title: Spectral propagator
---

This page is the record of a research prototype, `field.propagator = "spectral"`, and of what it measured. The slice propagator stays the default, every run of the keystone takes it, and the prototype is outside the keystone. The question the prototype asks is whether a field envelope propagated by one transform over the whole window, in $(x, y, \tau)$, reproduces the SASE tier within its tolerance with today's particle push and deposit untouched, and what that costs. The answer is below. The steady state is reproduced to rounding. The SASE tier as it is configured, with the source filter off, is not reproduced within its tolerance, and the difference is the wide-angle shot-noise light, which the frequency-dependent kernel delays as physics says it should and which the source filter removes, after which the two propagators agree to 1.6e-3 under whole-slice slippage and the continuous slippage adds a difference of a few $\rho$. The vacuum test shows the one place the slice propagator is wrong in closed form.

(spectral-what-changes)=
## What changes and what does not

The slice propagator ([](fel-physics.md#sec-field), [](fel-physics.md#sec-slippage)) holds the field as a stack of transverse grids, one a slice, and advances every slice over a step by the paraxial kernel $\exp(-i\,\delta z\, k_\perp^2 / 2k_0)$ at the one reference wavenumber $k_0 = 2\pi/\lambda_s$. The slippage rotates the stack by whole slices once $0.8$ of a slice has accumulated, and the light that leaves the head is banked.

The spectral propagator, `src/fel_spectral_mod.f90`, transforms a member's whole window in $(x, y, \tau)$, where $\tau$ is the window coordinate, $\tau = (i - 1)\,\Delta$ for window slice $i$ and slice spacing $\Delta$, so $\tau$ grows toward the head. Two things change.

- **The kernel is frequency dependent.** The envelope rides the carrier $\exp(i(k_0 z - \omega_0 t))$, the convention the kernel's sign fixes, so a mode $\exp(i\kappa\tau)$ of the envelope has the carrier's frequency plus $c\kappa$ and the wavenumber $k(\kappa) = k_0 + \kappa$. Its kernel is $\exp(-i\,\delta z\, k_\perp^2 / 2k(\kappa))$. A pulse that is short or chirped carries a spread of $\kappa$, and one $k_0$ for all of it is wrong at relative order $\delta\omega/\omega$.
- **The slippage is a phase ramp.** A shift of the field toward the head by $\delta$ is $\exp(-i\kappa\delta)$ in the spectral domain, exact and continuous, with no threshold and no rotation. Under `field.spectral_slippage = "slices"` the ramp is left out and the whole-slice rotation of the slice propagator runs after the pass, bank and all.

The push, the deposit, the source filter and the averaged equations are untouched. A step runs in the slice propagator's order: the kernel over the step, then the deposit of every slice's source into the record, then the slippage, which moves the field and the fresh source together. The deposit is `fel_spectral_deposit`, the per-particle scatter of `fel_field_step` kept apart from it so the slice propagator stays as it was, with the coherent source model refused.

The kernel factors: $\exp(-i\,\delta z\,(m_x^2 + m_y^2)\,\delta k^2 / 2k(\kappa))$ is the product of a table in $(m_x, \kappa)$ and one in $(m_y, \kappa)$, each `grid_n_pts` by the transform length, so the whole kernel costs one complex multiply a point and the tables are rebuilt at every call. The transverse transforms are the per-slice ones of `wavefront_fft2`, parallel over slices, and the transform along $\tau$ runs on each line of the transverse spectrum, parallel over the second index, through a plan cache in the same form, one plan pair and one buffer a thread.

In the steady state the window is one slice, the transform along $\tau$ is the identity, and the two propagators agree to rounding: on `examples/steady_state` the power curves differ by $2\times10^{-13}$.

(spectral-keys)=
## The keys

| Key | Default | What it does |
|---|---|---|
| `field.propagator` | `slices` | `spectral` selects the prototype. |
| `field.spectral_slippage` | `ramp` | The slippage as the continuous phase ramp, or `slices`, the whole-slice rotation. |
| `field.spectral_boundary` | `headroom` | What meets light the ramp moves to the window's head: `headroom`, `sponge` or `pad` ([](#spectral-boundary)). |
| `field.spectral_sponge_slices` | 0 | Slices of the sponge. 0 derives an eighth of the window, at least 4. |
| `field.spectral_head_tolerance` | 1e-3 | Under `headroom`, the head slice's share of the window's energy above which the run is refused. |

The prototype covers the averaged mode on the CPU in double precision. The unaveraged mode, the device, the FP32 check and the coherent source model each carry a propagator of their own, and each is refused with the key. The header prints a `Propagator` line that states the form in force.

(spectral-boundary)=
## The boundary at the head

A transform along $\tau$ is periodic, so light the ramp moves past the head would reappear at the tail. A derived window holds the line's slippage ahead of the bunch ([](fel-physics.md#sec-window)) and nothing reaches the tail. A stated window shorter than the line's slippage needs something at the head, and three forms are carried.

- **`headroom`.** Nothing is done, and after every pass the head slice's share of the window's energy is held below `field.spectral_head_tolerance`. Above it the run is refused, the message giving the position and the share, and the footer states the largest share a completed run saw and where.
- **`sponge`.** A $\cos^2$ weight falls from just under 1 on its first slice to 0 on the head slice, over `field.spectral_sponge_slices` slices, and the field is multiplied by it before each shift. The ramp then takes the slippage in chunks of at most one slice with the sponge before each, so what a chunk carries past the head is the sponge's zero on the head slice and nothing else. An interlude's whole slippage in one pass, ten slices on Aramis, would carry nine unsponged slices round to the tail, which the first form of the sponge did. The energy the sponge removes joins the member's escaped energy, and the light itself is gone.
- **`pad`.** The transform along $\tau$ is zero padded by eight slices beyond the largest shift one pass applies, an element's end slippage with its last step's. What a pass moves into the padding is removed, and its energy joins the escaped energy the same way. Under `spectral_slippage = "slices"` this is the hybrid: the rotation and today's bank carry the slippage and keep the light, and the padding meets the kernel's own spread alone, which is below one wavelength a step.

Under whole-slice slippage the rotation zero-fills the slice it transmits and the bank keeps it, so nothing wraps, and the headroom form checks nothing.

(spectral-measurements)=
## The measurements

```
python tests/scripts/spectral_field.py --exe build-release/bin/lucifer \
        --work-dir <dir> --validation-dir <a completed validation pass's work directory>
```

The script runs the first four measurements below by default and prints them, and `seeded` and `chirped` name the two FEL cases, which take hours at 12 threads. The vacuum test needs the tracker alone. The others take the tier inputs, the imported dumps and the Genesis4 references from a validation pass's work directory, which they read and never write. Every number on this page was measured on 2026-11-03 on the release build, gfortran 15.3 from conda-forge, on an Apple M3 Max with 12 cores, and the thread count of each table is stated with it.

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

The steady-state tier and the SASE tier of the validation suite ([](validation.md#val-validation-from-one-command)) rerun under `spectral`, against Genesis4 at the tier tolerances as `compare_fel.py` compares them, and against the validation pass's own slice run through `lucifer.compare`. The SASE tier is the full line from a dark start, 32 slices of 3 wavelengths, 1068 steps of 3 periods, and the line slips 3552 wavelengths, 1184 slices, so every undulator step slips exactly one slice and the ramp and the rotation differ only at the interludes, which slip whole wavelengths and not whole slices. The tiers run with the source filter off, since Genesis4 is compared with it off. Under the ramp the SASE window needs a boundary, and the padding is the one that discards what crosses as the rotation's zero fill does.

| Tier | Propagator | Against Genesis4 | Tolerance | Against the slice run: exit power, bunching, energy |
|---|---|---|---|---|
| `tier1` | slices, the recorded level | 1.825899e-06 | 1e-4 | |
| `tier1` | spectral, rotation | 1.825899e-06 | 1e-4 | 4.8e-15, 2.3e-14, 0 |
| `tier1` | spectral, ramp | 1.825899e-06 | 1e-4 | 4.8e-15, 2.3e-14, 0 |
| `tdsase` | slices, the recorded level | 2.349843e-06 | 1e-3 | |
| `tdsase` | spectral, rotation | 2.25e+07 | 1e-3 | 0.129, 0.167, 0.295 |
| `tdsase` | spectral, ramp, padding | 7.57e+06 | 1e-3 | 0.123, 0.199, 0.349 |

The steady state is reproduced to rounding, as [](#spectral-what-changes) says it must be. The SASE tier is not reproduced within its tolerance: the spectral propagator differs from the slice run by 13 percent of the exit power under either slippage, and the number against Genesis4 is the elementwise power at records where the slice run holds no light yet, floored at 1e-9 of the peak. Two further runs place the difference. With the kernel's frequency dependence removed by a temporary edit, $k(\kappa) = k_0$ for every mode, the spectral path reproduces the slice run to 7e-10 in power over the 1068 steps, so the transform along $\tau$, its order with the deposit and the rotation are right, and the whole difference is the kernel's $\kappa$ dependence. The script does not carry that edit. With the source filter on, the production default ([](fel-physics.md#sec-source-filter)), both propagators rerun against the slice run with it on:

| Propagator, source filter on | Exit power | Bunching | Energy |
|---|---|---|---|
| spectral, rotation | 1.6e-3 | 2.0e-4 | 4.1e-4 |
| spectral, ramp, padding | 9.2e-3 | 2.0e-3 | 6.1e-3 |

The filter removes the wide-angle emission of the point beamlets, and with it 99 percent of the difference. The per-particle deposit is white in $\kappa$ up to the slice Nyquist, a sixth of $k_0$ at 3 wavelengths a slice, and white in $k_\perp$ up to the grid's. The frequency-dependent kernel gives a ray at angle $\theta$ the group delay $\delta z\,\theta^2/2$ a step, 0.2 wavelengths a step at the grid's corner, which moves the wide-angle noise backward through the window where the slice propagator keeps every slice's field in its slice. The first step agrees exactly and the second differs by 4 percent in one slice, the delay acting on the first step's noise. The gain band's own modes sit at angles near 1e-6 rad and within $\rho$ of $k_0$, where the delay is 2e-4 wavelengths a step and the kernel differs at 5e-4, and that is the 1e-3 the filtered rotation run shows. The one-dimensional $\rho$ of the tier is 5.4e-4 (3000 A, $K = 1.2$, 17 µm rms size, 15 mm period). The continuous slippage adds 7.6e-3 of the exit power and 1.8e-3 of the bunching over the rotation, 14 and 3 times $\rho$, which is the difference the design expected by construction, since Genesis4's digits are measured with whole-slice rotation.

(spectral-boundary-measured)=
### The boundary forms

The SASE tier's window is a flat bunch of 32 slices under a line that slips 1184, so light reaches the head from the first step on and every slice, the head among them, holds a thirty-second of the noise at once. Each form ran on it with the bank on and the source filter off, as the tier runs, at 4 threads. The table gives the exit energies the footer states, the escaped energy and the window's energy at the exit, their sum, what the pulse file holds, and the head and tail slices' shares of the exit power.

| Form | Window energy at the exit | Escaped energy | Sum | Pulse file | Head slice | Tail slice |
|---|---|---|---|---|---|---|
| the slice propagator | 1.166e-10 J | 8.026e-9 J | 8.142e-9 J | 1216 planes, 8.142e-9 J | 9.6e-2 | 0 |
| spectral, ramp, `headroom` | refused at the first step: the head slice holds 3.2e-2 of the window's energy | | | | | |
| spectral, ramp, `headroom`, tolerance 2 | 5.538e-8 J | 0 | 5.538e-8 J | none, nothing left the window | 3.1e-2 | 3.1e-2 |
| spectral, ramp, `sponge` of 4 slices | 1.105e-10 J | 7.731e-9 J | 7.842e-9 J | none, the light is gone | 4.5e-4 | 9e-34 |
| spectral, rotation, `pad` of 18 slices | 1.347e-10 J | 8.293e-9 J | 8.428e-9 J | 1216 planes, 8.428e-9 J | 9.3e-2 | 0 |

What wraps under the headroom alone is everything: with no boundary the window is periodic, the light the beam emits passes through the beam again and again, every slice settles at a thirty-second of the window's energy, and the window holds 475 times the energy of the slice run. The refusal catches it at the first step on this window, since the headroom form is for a window with headroom. The sponge holds the head at 4.5e-4 of the exit power and the tail at nothing, and removes 7.73e-9 J against the 8.03e-9 J the rotation banks, so the light that escapes is counted to 4 percent and the pulse file holds none of it. The hybrid keeps the light, 1216 banked planes as the slice propagator banks, and its pulse carries 3.5 percent more energy than the slice run's, the kernel's difference of [](#spectral-null).

The beam's energy loss over the walk, each slice's live charge times its mean energy change, is 9.9e-11 J on the slice run and of that order under every form, 82 times smaller than the field's energy. The averaged equations charge the beam for the field at the particles and not for the incoherent emission the deposit adds, so on a run from noise the field's energy is not the beam's loss and no closure can be read off an averaged SASE run. The unaveraged mode's ledger is the closure statement ([](fel-physics.md#sec-unaveraged)), and the prototype does not reach that mode.

(spectral-cost)=
### The cost

The field phase of the SASE tier's walk, the `field solve` row of the footer's timing table, under both propagators, with the ramp and the padding of 18 slices under `spectral`, so the row holds the deposit, the kernel pass over transforms of 50 along $\tau$ and the ramp pass. The slice propagator's row holds the deposit and its transform pair. 255 by 255 grid, 32 slices, 1068 steps, one run each on an otherwise idle machine.

| Threads | Propagator | Field solve | A step | Walk |
|---|---|---|---|---|
| 1 | slices | 54.5 s | 51.1 ms | 77.2 s |
| 1 | spectral | 98.5 s | 92.2 ms | 122.9 s |
| 8 | slices | 7.45 s | 6.98 ms | 11.0 s |
| 8 | spectral | 14.3 s | 13.4 ms | 18.0 s |

The spectral propagator costs 1.8 times the slice propagator's field phase at one thread and 1.9 at eight. The transverse transform pair is the same count, and the rest is the two passes along $\tau$, each a gather of every line of the window across the stack's stride, a transform pair of 50 and a scatter, and the second pass only under the ramp.

(spectral-seeded)=
### A chirped seed through Aramis

The vacuum test's pulse with gain. The seed is the chirped Gaussian at 1 Å, rms length 15 wavelengths, waist 30 µm, 1 MW at its peak, ten times the start-up noise light a slice holds at the segment's exit, on the tier's grid, 30 slices from the tail of a window of 1280 slices of 3 wavelengths, which holds the line's slippage of 1184 slices so the seed stays with the beam to the exit. The beam is the tier's, 3000 A flat over the window at 1024 particles a slice, dumped by a run that tracks a marker alone and loaded with the seed. Relative rms bandwidths 0.5, 1, 2 and 5 percent, where 0.53 percent is the unchirped pulse's own, through the first segment and through the full line, under both propagators with the source filter on, at 12 threads. The segment also runs on the same bunch loaded without its shot noise, so the seed is the only light.

The slice propagator leaves its field a fraction of a slice behind the ramp's: the segment slips 88.67 slices, of which it rotates 88 and records the rest as the residual its field file carries, 2.0003 wavelengths, and the line slips 1184.00. The comparison moves the slice run's field by that residual along $	au$ through the transform before it differences, since a chirped field two thirds of a slice out of step differs by its phase slope times the offset and by nothing physical. The field difference is the largest over the lit slices and four on each side, over the spectral run's peak, and the raw column is the whole window as the runs left it. Each slice's wavelength comes from the phase advance to its neighbor on the slices above a quarter of the peak power, the seed's core.

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

With the seed alone, the propagators differ by 2.3e-3 to 3.0e-3 of the peak field at the three lower bandwidths and 6.5e-3 at 5 percent, and the peak powers agree to 1e-3. The floor is the slippage's form and not the kernel's: the spectral kernel with the whole-slice rotation differs from the slice propagator by 6.5e-4 at 0.53 percent and from the ramp by 2.6e-3, so what the beam sees over the segment is the field's position against it, up to a slice behind the ramp's between rotations. The part that grows with the bandwidth, 3.9e-3 above that floor at 5 percent, is the vacuum test's error scaled to this pulse, $0.77 \times 0.05 \times L / 2 z_R$ with the waist's Rayleigh range 28.3 m and the segment 3.99 m, 2.7e-3. With the beam's shot noise on, every segment row sits at 6.5e-2 to 8.3e-2 and the seed's own error is under it: the noise-only slices of the same runs differ by 8e-2 of the seed's peak and by 13 percent of their own rms, since the two kernels diffract the start-up light's wide angles differently before the filter and the gain narrow it. Through the line, where the seed has swept the window and the light is at saturation, the propagators differ by 1.6e-2 to 2.0e-2 of the peak field, 4e-4 in wavelength and 0.4 to 0.7 percent in peak power, with no trend in the bandwidth.

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

The vacuum test says where the prototype pays: a chirped or short pulse, where the slice propagator errs at the relative bandwidth. The two FEL cases say that error does not reach the output of an Aramis run through the gain: the seed's bandwidth error is under the floor the whole-slice slippage sets, that floor is under the start-up noise light's difference, and both are under what the step moves. Where the prototype differs from the slice propagator in a way the beam can see is the slippage, continuous under the ramp and a whole slice at a time under the rotation, which is what lets the field be sampled coarser than the beam's slices along $\tau$. The prototype keeps the envelope and its carrier, so the SVEA stands and the bandwidth the window carries is bounded by the slice spacing at $\lambda_s / 2\Delta$. A field without an envelope, Puffin's model, is the step after that one and is not this prototype.
