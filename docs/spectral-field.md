---
title: The spectral field propagator
short_title: Spectral propagator
---

The spectral field propagator, `field.propagator = "spectral"`, transforms the whole window in $(x, y, \tau)$. Each frequency of the field then diffracts at its own wavenumber, and the slippage is a continuous phase ramp. The slice propagator stays the default, and every tier of the keystone takes it. What the spectral propagator is for is a field sampled coarser than the beam's slices, `field.n_wavelength`, which only a continuous slippage can carry.

Every number on this page comes from one record, `docs/generated/spectral/measured.json`, which `tests/scripts/spectral_field.py all` writes, and `tests/scripts/report_spectral.py` draws the figures and writes the tables from it ([](#spectral-measurements)). The keystone regenerates them and requires no diff. The mode's own check in the keystone is the `spectral` section of the validation suite ([](validation.md#val-spectral)).

(spectral-when)=
## When to select it

Select the spectral propagator for a long window. Under it the field is sampled at a tenth of the cooperation length instead of at every beam slice. The transforms along the window then run over the field's planes alone. On the SASE tier's beam over 1280 slices through the Aramis line that sampling holds the pulse energy to 1.3 percent of the run at one plane a slice. The walk at 12 threads falls from 690 s to 198 s ([](#spectral-sweep)). The derived sampling needs no key.

Select it also for a seeded pulse whose bandwidth matters. The slice propagator diffracts every frequency at one wavenumber, which is wrong at first order in the relative bandwidth ([](#spectral-vacuum)). Through an FEL's gain that error does not reach the output of a chirped Aramis run at any bandwidth measured ([](#spectral-seeded)), so this reason is the smaller one.

Keep the slice propagator for a short window, where the cooperation length spans the whole window and the field keeps one plane a slice either way, and for every configuration in [](#spectral-limits). On a window of one slice the two propagators agree to rounding.

(spectral-keys)=
## The keys

| Key | Default | What it does |
|---|---|---|
| `field.propagator` | `slices` | `spectral` selects the spectral propagator. |
| `field.n_wavelength` | 0 | The field's plane along $\tau$ in wavelengths, a whole multiple of `window.n_wavelength`. 0 takes the window's slice under `slices` and derives a coarser plane under `spectral` ([](#spectral-coarse)). |

The header prints a `Propagator` line that states the kernel, the ramp and the pad in force, and a `Field` line with the planes when they are coarser than the slices. A stated `field.n_wavelength` that is not a whole multiple of the slice is refused. So is one coarser than the slice and above a quarter of the cooperation length, with the cooperation length and the cliff in the message. A stated value equal to the window's slice runs on any line, as the derived one does. On a line whose slice is itself above the cliff, the `flash1` and `lcls` examples among them, the header's `Field sampling` line says where the slice sits.

(spectral-what-changes)=
## How it works

The slice propagator ([](fel-physics.md#sec-field), [](fel-physics.md#sec-slippage)) holds the field as a stack of transverse grids, one a slice. It advances every slice over a step by the paraxial kernel $\exp(-i\,\delta z\, k_\perp^2 / 2k_0)$ at the one reference wavenumber $k_0 = 2\pi/\lambda_s$. The slippage rotates the stack by whole slices once $0.8$ of a slice has accumulated, and the light that leaves the head is banked.

The spectral propagator, `src/fel_spectral_mod.f90`, transforms a member's whole window in $(x, y, \tau)$. Here $\tau$ is the window coordinate, $\tau = (i - 1)\,\Delta$ for window slice $i$ and slice spacing $\Delta$, so $\tau$ grows toward the head. Two things change.

The kernel is frequency dependent. The envelope rides the carrier $\exp(i(k_0 z - \omega_0 t))$, so a mode $\exp(i\kappa\tau)$ of the envelope has the carrier's frequency plus $c\kappa$ and the wavenumber $k(\kappa) = k_0 + \kappa$. Its kernel is $\exp(-i\,\delta z\, k_\perp^2 / 2k(\kappa))$. A pulse that is short or chirped carries a spread of $\kappa$, and one $k_0$ for all of it is wrong at relative order $\delta\omega/\omega$. The kernel factors into a table in $(m_x, \kappa)$ and one in $(m_y, \kappa)$, so it costs one complex multiply a point.

The slippage is a phase ramp. A shift of the field toward the head by $\delta$ is $\exp(-i\kappa\delta)$ in the spectral domain, exact and continuous. A step may slip a fraction of a slice. The field is never a fraction of a slice behind the beam, where the rotation leaves it up to a slice behind between crossings of its threshold.

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-diagram-slippage -->
:end-before: <!-- end fig-sp-diagram-slippage -->
```

The push, the deposit, the source filter and the averaged equations are untouched. A step runs in the slice propagator's order: the kernel over the step, then the deposit of every slice's source into the record, then the slippage. The slippage moves the field and the fresh source together. In the steady state the window is one slice, the transform along $\tau$ is the identity, and the two propagators agree to rounding.

(spectral-bank)=
### The head and the bank

A transform along $\tau$ is periodic, so light the ramp moves past the head would reappear at the tail. The window is followed along $\tau$ by a pad of 4 free slots and a $\cos^2$ sponge of 12 slots that falls to zero at the pad's end. Every transform runs over the window and the pad together. Light the ramp carries past the head enters the pad, diffracts there as escaped light does in vacuum, and is attenuated by the sponge after every slip pass. The field is smooth at the periodic seam, and nothing reaches the tail.

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-diagram-pad -->
:end-before: <!-- end fig-sp-diagram-pad -->
```

The escaped-field bank ([](fel-physics.md#sec-slippage)) takes its planes from the pad's first slot. The slippage accumulates in the ramp's own count. Each time it crosses a whole plane, the pass is split there in the transform domain. The field at the pad's first slot is read by a direct sum over the modes, and the plane goes to the bank as the rotation's transmitted slice does. Between two crossings the ramp moves the field by exactly one slot. So the bank receives the window's samples one plane apart, each once, and the escaped and pulse files are assembled as before with the same readers. One pass slips at most one slice, so light walks through every slot of the sponge and reaches the seam below 1e-5 of its amplitude. An interlude's slippage is taken in such passes. A dump carries the window and not the pad, and a continuation starts with an empty pad, since every plane that had entered it was banked when it crossed. Backward slippage is refused, as it is with the bank.

(spectral-coarse)=
### The field's own sampling along tau

The window is the beam's partition, `window.n_wavelength` wavelengths a slice. Under the slice propagator the field has one plane a slice and no other sampling. Under the spectral propagator the field may be sampled coarser, `field.n_wavelength` wavelengths a plane, a whole multiple of the slice ([](fel-physics.md#sec-window)). A SASE field at saturation is band limited to a few $\rho$ about the carrier. A plane every fraction of the cooperation length $\lambda_s / 4\pi\rho$ carries it, while the beam's slice stays where the loader put it. A step then pays for the field at that sampling and for the beam at its slice. The continuous slippage is what makes this possible: a step slips a fraction of a plane, which the ramp applies exactly and the whole-slice rotation could not.

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-diagram-planes -->
:end-before: <!-- end fig-sp-diagram-planes -->
```

The field's record holds the window at that sampling, `fel_field_planes` planes of it, with one more past the last slice when the last slice falls between two planes. A beam slice reads the two planes around it by linear weights along $\tau$ (`fel_field_sample`). A slice on a plane reads that plane alone, which is every slice at one plane a slice, where nothing changes to the bit. The deposit splits each slice's source between the same two planes by the same weights and divides a plane's sum by the slices a plane spans. The deposit is then the readback's adjoint under each grid's own measure. The source filter runs once on each plane's sum (`fel_spectral_deposit_coarse`). The stats' field rows stay on the beam's axis, each slice's row taken from the plane it reads, so every reader and comparison works as before. The dumps, the frames and the pulse file carry the field's own planes at their spacing, and `lucifer.field` reads that spacing as it reads the slice's. The bank banks a plane each time the slippage crosses a whole plane, and the header's `Bank` line counts planes. A field file at another spacing is refused unless it is dark. A dark file gives its transverse grid alone, and the field starts dark on its own planes, so a SASE input file whose field file fixes the grid runs at any sampling.

The sampling is derived at setup when `field.n_wavelength` is unset: a tenth of the cooperation length at the $\rho$ of the beam as built, rounded to whole slices and cut so the window keeps 8 planes. The $\rho$ comes from the slice with the most charge, its rms size and the first undulator's coupling (`fel_pierce_of_beam`), and the header prints the sampling with its origin. A stated value is held to a whole multiple of the slice, to at least 8 planes over the window, and, when it is coarser than the slice, below the cliff the sweep found, a quarter of the cooperation length ([](#spectral-sweep)). The window's slice itself is never refused, since on a line whose slice is above the cliff it is what the derivation returns.

(spectral-measurements)=
## What it is measured to do

```
python tests/scripts/spectral_field.py --exe build-release/bin/lucifer \
        --work-dir <dir> --validation-dir <a completed validation pass's work directory> \
        --threads 12 --cost-threads 1,8 --toolchain "<the compiler>" \
        --json docs/generated/spectral/measured.json all
```

The script runs the eight measurements below and writes every number and curve into the record as each finishes. It takes the tier inputs, the imported dumps and the Genesis4 references from a validation pass's work directory, which it reads and never writes. The vacuum test needs the tracker alone. The two FEL cases and the sweep take hours at 12 threads, and the cost rows need an otherwise idle machine. So the record is re-measured by hand as [](performance.md) is, and never in the keystone. The tracker refuses the sweep's rows coarser than the slice above the cliff, so they run on a second tracker built with `fel_field_coop_cliff$` raised, given as `--exe-open`. Each line's first row, one plane a slice, runs on the tracker itself. The record's first rows on `flash1` and `lcls` were measured on the second tracker and are marked so. The cliff's constant moves nothing but the refusal, so their numbers are the tracker's own. Each measurement can be run alone by its name, and `--reuse` re-analyzes a finished set from its files. `report_spectral.py` then writes the tables and draws the figures, and the keystone runs it and diffs the tables and the captions.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-provenance -->
:end-before: <!-- end tab-sp-provenance -->
```

(spectral-vacuum)=
### The vacuum test

A Gaussian pulse with a linear chirp, $\exp(-r^2/w_0^2)\exp(-\xi^2/4\sigma^2 + i b\,\xi^2)$ with $\xi$ the distance from the pulse center, drifts through a pipe of two Rayleigh ranges with no beam, under both propagators, against the exact paraxial solution. Each frequency of the pulse is a Gaussian beam with its own Rayleigh range $k w_0^2/2$, and the exact solution superposes them over the pulse's analytic spectrum, on a wavenumber grid eight times finer than the window carries, read at the slices.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-vacuum -->
:end-before: <!-- end tab-sp-vacuum -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-vacuum -->
:end-before: <!-- end fig-sp-vacuum -->
```

This is the one place the slice propagator is wrong in closed form. Over two Rayleigh ranges the red and the blue ends of the pulse diffract to different sizes and phases, and the error is at the bandwidth and not at the bandwidth squared. A SASE pulse at saturation carries a bandwidth of order $\rho$, so the effect on the tiers is at that order and no larger, which the null test measures.

(spectral-null)=
### The null test

The steady-state tier and the SASE tier of the validation suite ([](validation.md#val-validation-from-one-command)) rerun under `spectral`, against Genesis4 at the tier tolerances as `compare_fel.py` compares them, and against the validation pass's own slice run through `lucifer.compare`. The SASE tier is the full line from a dark start, 32 slices of 3 wavelengths, 1068 steps of 3 periods, and the line slips 3552 wavelengths, 1184 slices, so every undulator step slips 0.9963 of a slice and the interludes whole wavelengths. The tiers run with the source filter off, since Genesis4 is compared with it off.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-null -->
:end-before: <!-- end tab-sp-null -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-null -->
:end-before: <!-- end fig-sp-null -->
```

The steady state is reproduced to rounding. The SASE tier as it runs is not reproduced within its tolerance. The number against Genesis4 is the elementwise power at records where the slice run holds no light yet, floored at 1e-9 of the peak, and the difference from the slice run is the kernel's $\kappa$ dependence: with it removed by a temporary edit, $k(\kappa) = k_0$ for every mode, the spectral path reproduced the slice run to 7e-10 in power over the 1068 steps, so the transform along $\tau$, its order with the deposit and the slippage are right. The script does not carry that edit.

The source filter, the production default ([](fel-physics.md#sec-source-filter)), removes most of that difference. The per-particle deposit is white in $\kappa$ up to the slice Nyquist, a sixth of $k_0$ at 3 wavelengths a slice, and white in $k_\perp$ up to the grid's. The frequency-dependent kernel gives a ray at angle $\theta$ the group delay $\delta z\,\theta^2/2$ a step, 0.2 wavelengths a step at the grid's corner, which moves the wide-angle noise backward through the window where the slice propagator keeps every slice's field in its slice. The gain band's own modes sit at angles near 1e-6 rad and within $\rho$ of $k_0$, where the delay is 2e-4 wavelengths a step. What is left with the filter on is the slippage's form, continuous against a whole slice at a time, which grows toward the head of the window where the most light has slipped in. The chirped seed below places that form's difference at a few 1e-3 over one segment with the seed the only light.

(spectral-bank-measured)=
### The bank

The SASE tier's window is a flat bunch of 32 slices under a line that slips 1184, so light reaches the head from the first step on, and every slice, the head among them, holds a thirty-second of the noise at once. Both propagators run on it with the escaped field kept and the source filter off, as the tier runs.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-bank -->
:end-before: <!-- end tab-sp-bank -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-bank -->
:end-before: <!-- end fig-sp-bank -->
```

The pulse check of the validation suite closes the same ledger on a Gaussian pulse through a pipe ([](validation.md#val-spectral)). The ramp's ledger closes less tightly than the rotation's there, since the frequency-dependent kernel moves light along $\tau$ across the head on its own, by the group delay of the wide angles, and the bank samples the ramp's crossings alone. The beam's energy change over the walk, each slice's live charge times its mean energy change, is far below the field's energy under either propagator. The averaged equations charge the beam for the field at the particles and not for the incoherent emission the deposit adds, so on a run from noise the field's energy is not the beam's loss and no closure can be read off an averaged SASE run. The unaveraged mode's ledger is the closure statement ([](fel-physics.md#sec-unaveraged)).

(spectral-seeded)=
### A chirped seed through Aramis

The vacuum test's pulse with gain. The seed is the chirped Gaussian at 1 Å on the tier's grid, 30 slices from the tail of a window of 1280 slices of 3 wavelengths, which holds the line's slippage of 1184 slices so the seed stays with the beam to the exit. The beam is the tier's, 3000 A flat over the window at 1024 particles a slice, dumped by a run that tracks a marker alone and loaded with the seed. Relative rms bandwidths 0.5, 1, 2 and 5 percent, where 0.53 percent is the unchirped pulse's own, run through the first segment and through the full line under both propagators with the source filter on. The segment also runs on the same bunch loaded without its shot noise, so the seed is the only light.

The slice propagator leaves its field a fraction of a slice behind the ramp's: the segment slips 88.67 slices, of which it rotates 88 and records the rest as the residual its field file carries. The comparison moves the slice run's field by that residual along $\tau$ through the transform before it differences, since a chirped field two thirds of a slice out of step differs by its phase slope times the offset and by nothing physical. The field difference is the largest over the lit slices and four on each side, over the spectral run's peak, and the raw column is the whole window as the runs left it. Each slice's wavelength comes from the phase advance to its neighbor on the slices above a quarter of the peak power, the seed's core.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-seeded -->
:end-before: <!-- end tab-sp-seeded -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-seeded -->
:end-before: <!-- end fig-sp-seeded -->
```

With the seed alone, the floor is the slippage's form and not the kernel's: the spectral kernel with the whole-slice rotation, a form since retired, differed from the slice propagator by 6.5e-4 at 0.53 percent and from the ramp by 2.6e-3, so what the beam sees over the segment is the field's position against it. The part that grows with the bandwidth is the vacuum test's error scaled to this pulse, $0.77 \times B \times L / 2 z_R$ with the waist's Rayleigh range 28.3 m and the segment 3.99 m. With the beam's shot noise on, the seed's own error is under the start-up light's, since the two kernels diffract that light's wide angles differently before the filter and the gain narrow it. Through the line, where the seed has swept the window and the light is at saturation, the propagators differ by percent of the peak field, with no trend in the bandwidth.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-convergence -->
:end-before: <!-- end tab-sp-convergence -->
```

Through the line, cutting the step from three periods to one moves the 2 percent seed's run by 5.8e-2 of the peak field, more than the propagators' 1.6e-2 to 2.0e-2 at any bandwidth. Over the segment with the shot noise on the step moves it by 3.4e-2 and the propagator by 8.3e-2, since the start-up light's wide angles are what the two kernels treat differently. The gain guiding holds the mode at the beam's size through the line, so the vacuum test's growth over the Rayleigh range does not accumulate, and the bandwidth error the vacuum test isolates is below the floor the whole-slice slippage sets against the beam.

(spectral-chirped)=
### An energy-chirped beam from noise

The same beam with a linear energy chirp along the window, 0.25, 0.5 and 1 percent head to tail, the head above the tail, and the unchirped bunch beside them, from noise through the full line under both propagators with the source filter on. The bunch is the dump's with every momentum component and the momentum deviation record scaled by the same factor, so the angles stay. A chirp $\epsilon$ in energy chirps the resonant wavelength by $-2\epsilon$, the head shorter. The window is 1.28 fs, so 1 percent across it moves the resonance by 1.6 $\rho$ over a cooperation length of 52 slices.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-chirped -->
:end-before: <!-- end tab-sp-chirped -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-chirped -->
:end-before: <!-- end fig-sp-chirped -->
```

The chirp decides the output and the propagator does not. At 0.25 percent the light slipping forward meets electrons of rising energy, which is what a taper does after saturation, and the exit energy rises with a single spike. At 0.5 and 1 percent the resonance moves across a cooperation length by more than the gain bandwidth and the exit energy falls. The two propagators agree on every one of those numbers. The field difference grows with the chirp, which is the saturated light's own sensitivity to any change of the operations: the unchirped run's 2.1e-2 is the seeded line's level, and cutting the step to a third moves the seeded line by 5.8e-2.

(spectral-sweep)=
### The sweep

Four lines from noise under the spectral propagator, each at 1, 2, 5, 10, 20 and 40 slices a plane, against its own run at one plane a slice: the SASE tier's beam over 1280 slices of 3 wavelengths through Aramis, which is the seeded case's line and bunch with no seed, and the `sase`, `flash1` and `lcls` examples as they ship. The spectrum's rms width is also taken and is not a measure here: at one plane a slice the window's Nyquist band is a sixth of $k_0$ and the shot noise's floor across it sets the width, which then falls with the band as the plane coarsens.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-sweep -->
:end-before: <!-- end tab-sp-sweep -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-sweep -->
:end-before: <!-- end fig-sp-sweep -->
```

The long Aramis window is the measurement. Its 1280 slices hold enough spikes for the energy to be a mean, and its slices are fine enough to place the plane from a fiftieth of the cooperation length to most of it. The pulse energy holds to 1.3 percent up to a tenth of the cooperation length and then falls, and the exit power per slice differs by a few percent of the peak at the first samplings, the saturated light's own sensitivity. The push then holds what is left of the walk. The other three lines agree with that curve where they reach it. FLASH1 and LCLS lose output at the first coarsening, and their own slices, 12 and 68 wavelengths, already sit at 0.3 and 0.5 of their cooperation lengths, so the derived sampling keeps their slice. The `sase` example's window of 96 slices holds two cooperation lengths, so its energy moves by 10 percent either way between samplings, a realization's own spread and not a trend, which is why a short window is no measure of the sampling.

From that curve the derived sampling is a tenth of the cooperation length, `fel_field_coop_fraction$`, and the cliff is a quarter, `fel_field_coop_cliff$`, past which the energy falls by percent and a stated `field.n_wavelength` coarser than the slice is refused. On the SASE tier's window of 32 slices the derived sampling is cut to 4 slices a plane so the window keeps 9 planes, and the `spectral` section of the validation suite holds that run against one plane a slice ([](validation.md#val-spectral)).

(spectral-unaveraged)=
### The unaveraged mode

The unaveraged mode ([](fel-physics.md#sec-unaveraged)) advances a record step as substeps that resolve the quiver, each slice diffracting at $k_0$ after its own substep, and the walk slips the record by rotation after the step. Under the spectral propagator the substeps run as they do, and after them the kernel's remainder at $k_0 + \kappa$ and the ramp run once over the window and its pad (`fel_spectral_correct`, `fel_spectral_slip`), the remainder being the kernel with its $k_0$ part divided out, which commutes with the substeps' diffraction. The sandwich probe measures it: an averaged segment, a pipe, an unaveraged segment, a pipe and an averaged segment, over 200 slices of one wavelength from shot noise, so the unaveraged segment's light reaches the trailing segment, under both propagators, with the all-averaged twin under both as the level to read against.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-unavg -->
:end-before: <!-- end tab-sp-unavg -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-unavg -->
:end-before: <!-- end fig-sp-unavg -->
```

The mode under the ramp sits where the averaged line sits between the propagators. A second form was measured on the same probe and not taken: the kernel at $k_0 + \kappa$ and the ramp at every substep, the whole window transformed twenty times a record step. It departed from the slice propagator's run by 1.3e-1 of the peak field and 0.8 percent of the exit power, from the form above by 1.1e-1, and cost 3.5 times the walk, 6.4 s against 1.8 s under the slice propagator, with its ledger at 2.4e-7. The field holds one plane a slice, so a ramp by a twentieth of a slice at a time moves the band-limited interpolant of that sampling, which is the departure, and the mode's own emission is white to the slice's Nyquist and feels it most.

The ledger needs one thing under the ramp. Light between the window's head and its next crossing sits in the pad, and so does light the kernel's own transport along $\tau$ carries past the head, in neither the window's column nor the bank's, up to a plane's energy at a time. On the eight-slice window of the validation suite that is 12 percent of the field at a record, and it is recovered at the next crossing, so a ledger on the bank's count reads 1e-2 while nothing is lost, which a seeded window without charge shows conserving its energy to 2e-8. Under the spectral propagator the ledger's field column therefore counts the pad's planes with the window's and its escaped column the energy the sponge has removed for good, both read from the field before and after the step's passes, since the kernel and the ramp are unitary on the window and pad together. The ledger then closes as the table shows.

(spectral-cost)=
## What it costs

Three things are added to a step. The transverse transform pair runs on the pad's 16 planes beside the window's, half again the count on the tier's 32 and 2 percent on a window of a thousand. Each step takes two passes along $\tau$, the kernel's and the slip's, each a gather of every line of the window and pad across the stack's stride, a transform pair and a scatter, with the sponge and a plane's direct sum in the slip pass. And an interlude's slippage of up to 17 slices takes a pass a slice, 400 passes over the Aramis line's 24 interludes beside the steps' 2136.

```{include} generated/spectral-measured.md
:start-after: <!-- begin tab-sp-cost -->
:end-before: <!-- end tab-sp-cost -->
```

```{include} generated/spectral-measured.md
:start-after: <!-- begin fig-sp-cost -->
:end-before: <!-- end fig-sp-cost -->
```

The header's cost estimate counts what the mode runs ([](performance.md)): the pad's planes in the transverse transform at the solve's rate, and the passes along $\tau$, two a step and the interludes' own, at `fel_cost_tau$`, the rate the cost table's one-thread rows give. Under a coarse field the counts are the field's planes and not the beam's slices, so the estimate falls with the sampling. The walk falls with it as the sweep's last column shows.

(spectral-limits)=
## What it does not cover

- The unaveraged mode on the device, whose unaveraged kernels step the slice propagator (`fel_device_mod`), and the FP32 check with the unaveraged mode under the propagator, since the device would be its twin.
- The FP32 check without the device, whose CPU twin steps the slice propagator in single precision. With `run.device = "metal"` the device takes the twin's role ([](validation.md#val-device-spectral)).
- The coherent source model, whose deposit is `fel_field_step`'s own.
- Harmonic fields and `run.dump_reduced` with a field coarser than the slices. The coarse field is carried for the fundamental alone, and the field reductions are taken plane by plane of the record and not on a slice's interpolated plane.
- A field coarser than the slices in a line with an unaveraged segment, whose source resolves structure below a slice. The message names the slice as the sampling to state.
- Backward slippage, since the pad that feeds the bank is at the window's head.

Each is refused at setup with the reason in the message. The spectral propagator keeps the envelope and its carrier, so the SVEA stands and the bandwidth the window carries is bounded by the slice spacing at $\lambda_s / 2\Delta$. A field without an envelope, Puffin's model, is a different propagator and is not this one.

The device carries the propagator ([](validation.md#val-device-spectral)). Its step's command buffer runs the kernel pass along $\tau$ over the window and its pad on the resident record, the slip passes follow it with the plane at the pad's first slot read at each crossing and downloaded for the bank, and the coarse field's push and deposit read and write the two planes about a slice by the weights above. Harmonic members and two polarization planes ride under it at one plane a slice. The device's own limits are two. The window and its pad are transformed along $\tau$ as one length in the threadgroup memory, 2048 planes on an M3 Max, and a longer one is refused with both counts. The resident record is counted against the device's working set before anything is allocated, and one that does not fit is refused with the bytes it needs and the bytes the device offers. The header's `Device` line states the resident bytes beside the bank's. The GPU time of the transform along $\tau$, before and after its stages were compiled into the kernel, is in [](performance.md#perf-device-tau).

(spectral-next)=
## What is next

The unaveraged mode on the device under the propagator, which needs the transform along $\tau$ inside the mode's own command buffer. The transform's GPU time per pass is in [](performance.md#perf-device-tau). On a window above 1024 planes a threadgroup transforms one line, and the measurement there found two lines no faster.
