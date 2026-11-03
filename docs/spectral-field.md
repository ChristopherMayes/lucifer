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

The script runs the four measurements below and prints them. The vacuum test needs the tracker alone. The other three take the tier inputs, the imported dumps and the Genesis4 references from a validation pass's work directory, which they read and never write. Every number on this page was measured on 2026-11-03 on the release build, gfortran 15.3 from conda-forge, on an Apple M3 Max with 12 cores, and the thread count of each table is stated with it.

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

(spectral-next)=
## What is next

The vacuum test says where the prototype pays: a chirped or short pulse, where the slice propagator errs at the relative bandwidth. A chirped FEL case is the next measurement. The prototype keeps the envelope and its carrier, so the SVEA stands and the bandwidth the window carries is bounded by the slice spacing at $\lambda_s / 2\Delta$. A field without an envelope, Puffin's model, is the step after that one and is not this prototype.
