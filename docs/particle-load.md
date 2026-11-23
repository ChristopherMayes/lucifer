---
title: "The particle load"
short_title: Particle load
tags: [explanation]
---

The [load](glossary.md#gl-load) is the set of [macroparticles](glossary.md#gl-macroparticle) with which a run begins. It is generated once, before the first step, from a bunch: either the description of a bunch in the `beam_init` structure of Bmad, or a file of particles. The bunch is divided by arrival time into the [slices](glossary.md#gl-slice) of the [window](glossary.md#gl-window), and the key `beam.load_mode` determines how the particles of a slice represent the charge that arrived in it.

The page describes the method: the three modes and the two switches ([](#pl-modes)), the sequence of operations that generates a load and the modes in turn ([](#sec-loading)), the [quiet start](glossary.md#gl-quiet-start) ([](#sec-quiet)), the [shot noise](glossary.md#gl-shot-noise) imposed on a load ([](#sec-noise)), and the [resampler](glossary.md#gl-resampler) of a bunch from a file ([](#sec-import)). The page follows from the method and from the source of the program, and it contains no result of a run. Its last section, [](#pl-evidence), states where the measurements of the load are recorded.

Three pages continue this one. [](shot-noise.md) derives in closed form the noise of each load. [](electrons-mode.md) describes the mode in which every electron is a particle. [](migration.md) describes the [migration](glossary.md#gl-migration) of a particle from one slice to another during the run, whose effect depends on the load.

(pl-modes)=
## The three modes and the two switches

The key `beam.load_mode` has three values. They differ in the number of particles of a slice and in the weight of a particle, which is the charge that it represents.

| `beam.load_mode` | particles in a slice | weight of a particle |
|---|---|---|
| `"sample"`, the default | `beam.n_particle`, or fewer in a slice of few electrons ([](#pl-fewer)) | the charge of the slice divided by the number of its particles, before the noise is imposed |
| `"keep"` | the particles of the bunch that arrived in the slice, each as $b$ copies under the quiet start | the weight of the particle of the bunch, divided by $b$ under the quiet start |
| `"electrons"` | the electrons of the slice | $Q/N$ for a window of charge $Q$ and $N$ electrons, which differs from $e$ by at most one part in $2N$ |

Here $b$ is the number of macroparticles in a [beamlet](glossary.md#gl-beamlet), the key `beamlet_size`.

Two switches act on a load of the first two modes. The switch `quiet_start` arranges the macroparticles in beamlets, whose bunching vanishes ([](#sec-quiet)). The switch `shot_noise` imposes on the load the shot noise of the electrons that it represents ([](#sec-noise)). Their four combinations are the following.

| `quiet_start` | `shot_noise` | the load |
|---|---|---|
| true | true | Quiet beamlets with the shot noise of the beam imposed on them. This is the load of a SASE run. |
| true | false | Quiet beamlets and no noise. This is the load of a seeded run, or of a run with the coherent source. |
| false | true | The groups of the load as they are, with the shot noise of the beam imposed on them by the first-order kick, which the input file must select. This is the load of a bunch that was made quiet by other means. |
| false | false | The bunch as it is, with the noise of its macroparticles. |

In the electrons mode neither switch is used: the particles are independent, and the shot noise is a property of the load.

The key `beam.file` is separate from all of the above. It names a dump of slices, which is loaded as it is, together with the field of `field.file`.

(sec-loading)=
## Loading

A single sequence of operations generates every load. The bunch is obtained from the `beam_init` structure of Bmad, which either describes it or, through `beam.position_file`, names a file of particles in the openPMD-beamphysics format. The bunch is divided by arrival time into the slices of [](#sec-window), and `beam.load_mode` determines what is done with the particles that arrived in a slice. The switches `quiet_start` and `shot_noise` then act on the result.

The arrival coordinate is $\tau$ of Eq. [](#eq-kinematics) of [](#sec-chart), shifted so that the earliest particle is at zero, and the spacing of the slices is $\Delta s = \texttt{n\_wavelength}\, \lambda_0$. The Bmad routine `init_beam_distribution` generates the bunch from either source, so that a file is read in full. That routine applies `beam.n_particle` as a limit on the number of particles read from a file. The limit is suppressed, since in Lucifer the key is the number of macroparticles in a slice.

### The sample mode

With `load_mode = "sample"`, the default, a slice contains `beam.n_particle` macroparticles, which includes the phase copies of the quiet start. The weight of a macroparticle is the charge of its slice divided by that number, so that the current of the load follows that of the bunch. A slice of few electrons contains fewer macroparticles, and a slice with no charge contains none, as [](#pl-fewer) states. The two loads of Genesis4, `&beam` with a `&profile` and `importdistribution`, both correspond to this mode.

(pl-described)=
#### A bunch from its description

When no file is named, the loader evaluates `beam_init` analytically for each slice, and no bunch is generated. In each slice it generates $m = n_{\mathrm{part}}/b$ transverse samples from the distribution matched to the Twiss parameters of the lattice, where $n_{\mathrm{part}}$ is `beam.n_particle` and $b$ is `beamlet_size`, and assigns one energy to each. It uses the fields `n_particle`, which is the number per slice, `a_norm_emit` and `b_norm_emit`, `sig_pz`, where $\delta\gamma = \beta_0 p_0 \,\texttt{sig\_pz}/mc$, `bunch_charge`, `sig_z` and `distribution_type(3)`. A value in any other field of the structure is refused, since a field that was set and then ignored would describe a bunch other than the one that is loaded. The Twiss parameters are always those of the lattice ([](#sec-element)). The current is derived from the bunch and is not an input:

$$
  \begin{aligned}
    I(s_i) &= \frac{Q c}{\sqrt{2\pi}\,\sigma_z}\, e^{-s_i^2/2\sigma_z^2}
      && \text{(Gaussian)},\\[4pt]
    I &= \frac{Q c}{x_{\max}-x_{\min}}
      && \text{(GRID: Bmad's uniform, flat)},
  \end{aligned}
$$ (eq-derivedcurrent)

The current is evaluated at the centers of the slices, with the bunch centered in the window. Where the input file gives `window.current`, that uniform current is used in its place. The value $\sigma_z = 0$ denotes the steady state, in which the whole charge is in the one slice and $I = Qc/\Delta s$, and it is refused for a time-dependent window.

When the window is not stated, it covers the described bunch. For a Gaussian bunch the window extends to the position at which a slice would contain less than one electron: its length is $2\sigma_z\sqrt{2\ln N_0}$, where $N_0$ is the number of electrons in the central slice. A bunch whose central slice contains less than one electron has no such length, and its window is $8\sigma_z$ long. For a uniform bunch the window is the extent of the grid. The key `window.window_length` replaces this length, for example to provide room for [slippage](glossary.md#gl-slippage), and a warning states how much of the bunch a shorter window excludes.

(pl-file)=
#### A bunch from a file

With `beam.position_file` the bunch is resampled slice by slice ([](#sec-import)). The macroparticles of a slice are generated from the particles of the bunch within the sampling window of the slice, and the current of the slice is obtained from the same window.

### The keep mode

With `load_mode = "keep"` every particle of the bunch is retained. Under the quiet start each particle becomes a beamlet of $b$ copies of weight $w/b$, which share its five other coordinates. Without the quiet start it remains one particle of its whole weight $w$. The number of macroparticles in a slice therefore follows the charge: it is the number of particles of the bunch that arrived in the slice, multiplied by $b$ under the quiet start. The window is $\lceil \tau_{\mathrm{total}}/\Delta s \rceil$ slices from the earliest particle, unless the input file states `window.n_slice` or `window.window_length`. The number is rounded up so that no particle is outside the window, whereas the resampler, which follows Genesis4, rounds to the nearest integer ([](#sec-import)). In a stated window the bunch is centered, and a warning gives the number of particles outside it. No random number is used before the noise is imposed, so that a given bunch always produces the same load.

The copies are placed at $\theta_0 + 2\pi j/b$, $j = 0,\ldots,b-1$, where $\theta_0 = k_s (s - s_i) \bmod 2\pi\,\texttt{n\_wavelength}$ is the phase of the particle within its slice, $s_i$ being here the position at which the slice begins, and $z = \beta\theta/k_s$ is the coordinate that is stored ([](#sec-chart)). The copies therefore occupy the one wavelength that begins at the phase of the particle, and those of a particle near the end of its slice extend beyond that end by less than one wavelength. The transverse momenta are exact, $P_x/(mc)$ divided by $p_0/(mc)$, whereas the resampler uses $\gamma x'$ ([](#sec-import)). Every charge-weighted moment of the five coordinates other than the phase is consequently a moment of the load, slice by slice. This property is the purpose of the mode, which is intended for a bunch obtained from a simulation of the accelerator. The current of each slice and its five first moments are written to `<out_root>.import.txt` at full precision, and the checks of this exactness read them there ([](validation.md#val-the-load-path)).

(sec-electrons)=
### The electrons mode

With `load_mode = "electrons"` the load contains one particle for every electron of the bunch. The charge $Q$ of the window is counted as the nearest whole number of electrons, $N = \mathrm{nint}(Q/e)$, and every particle has the weight $Q/N$. The numbers of electrons in the slices are one multinomial draw, and each electron is assigned transverse coordinates, an energy and a phase of its own. There are no beamlets. The bunching of a slice of $N_i$ electrons has the mean square $1/N_i$ at every harmonic, so that the shot noise is a property of the load and nothing is imposed on it. The keys of the beamlets and of the imposed noise are not used, and a value other than the default in any of them is refused.

The particles of a run occupy 68 bytes per electron, 0.42 GB per pC of charge in the window. The load, its keys, the number of particles and their memory, the run time, and the measurements on the examples, among them the comparison with the weight load, are the subject of [](electrons-mode.md).

(sec-quiet)=
## Quiet start

With `quiet_start = true`, the default, the load consists of beamlets. A beamlet is one transverse sample and one energy at $b$ phases, $\theta = \theta_0 + j\,2\pi/b$, so that its contribution to the bunching vanishes, to rounding, at every harmonic below $b$ ([](#fig-shot-phasors-quiet)). Before the noise is imposed, the weights within a beamlet are equal, since the cancellation is that of the single beamlet. The averaged step depends on $z$ only through $\theta$, which is periodic with the period $2\pi$. A beamlet that is displaced as a whole therefore remains quiet at every harmonic, whereas a beamlet whose members are displaced singly does not. The [quiet floor](glossary.md#gl-quiet-floor) of [](#sec-noise) measures the second case.

### A slice of several wavelengths

A slice is `window.n_wavelength` wavelengths long and contains the charge of all of them, and its macroparticles are distributed over all of its periods (`fel_slice_fill`, [](#fig-shot-slice-fill)). The $N = b\,n_b$ macroparticles of a slice of $n_b$ beamlets lie on one uniform grid over the slice, at $(j + \tfrac12)\,n/N$ periods for $j = 0, \ldots, N-1$, where $n$ is the number of wavelengths of the slice. A beamlet consists of the $b$ grid points whose phases within a period differ by a multiple of $2\pi/b$. It cancels at every harmonic as before, and its points lie in several periods of the slice. The charge in any interval of the slice is then uniform to within one macroparticle, the weights are equal, and a dump of the beam has no structure at the length of the slice: the sum of $e^{-ikz}$ over a quiet slice vanishes, to rounding, at the multiples $m k_s/n$ with $m = 1, \ldots, N-1$. In the keep mode a particle is placed at its own phase within its slice.

For a slice of one wavelength the arrangement places the particles as a load within one period does. The averaged step is insensitive to the period in which a particle lies, whereas a dump that is tracked further by another code is not.

### A load without the quiet start

With `quiet_start = false` no beamlets are loaded. The analytic generator assigns to every particle its own five coordinates and a phase distributed uniformly. The resampler retains the coordinates of the copies and assigns an independent phase to each. The keep mode retains the phase and the whole weight of each particle. Each of these loads has the noise of its macroparticles at every harmonic, as a bunch of that number of electrons would have. A bunch that was made quiet by other means is loaded in this way, and `shot_noise` then identifies its groups ([](#pl-groups)).

### The number of beamlets and the size of a beamlet

When `beam.n_particle` is not set in the sample mode, the load has 512 beamlets per slice, each of the beamlet size. The default is a number of beamlets since the result depends on the number of beamlets in a slice and not on the number of macroparticles. The value 512 was chosen by measurement, which [](shot-noise.md#shot-count) records.

The key `beamlet_size` is a quantity of the load. When it is not set, it is derived from the harmonics that the [field set](glossary.md#gl-field-set) contains, as twice the highest harmonic plus 2: 4 for the fundamental alone, 8 with the third harmonic and 12 with the fifth. In the sample mode it is the smallest divisor of `beam.n_particle` that is not less than this value. The unaveraged mode uses the same size. For evenly spaced phases the harmonic $b - h$ is not distinguished from the harmonic $h$, so that a smaller beamlet gives a harmonic of the field set the bunching of a lower one. The measurements of this effect are in [](shot-noise.md#shot-size).

During a run the size is used once, in the number of independent transverse samples of the coherent source ([](#sec-coherent-source)). That number is the number of distinct transverse positions in the slice, which is exact for copies, and a dump is therefore loaded without the key.

(sec-noise)=
## The imposed shot noise

With `shot_noise = true` the shot noise of the beam is imposed on the load by the method that `beam.shot_noise_method` names. The weight method applies to the beamlets of this loader, and the first-order kick applies to any load.

### The weight method

The default, `"weight"`, applies to the beamlets that this loader generates. No phase is changed. Each macroparticle is assigned a whole number of electrons, distributed as the electrons of a beam are distributed among the macroparticles, with the total of the bunch fixed: the numbers are one multinomial draw of the electrons of the bunch over its macroparticles. The bunching of a slice is then a sum over independent electrons, and $\langle|b(h)|^2\rangle = 1/N_\lambda$ holds exactly at every harmonic below `beamlet_size` and for any number of electrons, where $N_\lambda$ is the number of electrons of the slice ([](shot-noise.md#shot-weight)). The charge of a slice fluctuates about its mean as that of a beam does, and the loaded charge is exactly that of the bunch.

In the keep mode under the quiet start, a beamlet is one particle of the bunch with its copies, and it retains its own charge. A particle of $N$ electrons, with $N \ge 1$, is represented by $\lfloor N \rfloor$ electrons or by one more, with the probability for which the mean of the reciprocal of that number is $1/N$, and the electrons are assigned to the copies with equal probability. A particle of less than one electron is represented by one electron ([](shot-noise.md#shot-weight)).

A load with `quiet_start = false` has no beamlets of this loader, in the sample mode and in the keep mode alike. The weight method refuses such a load, and the message names the two remedies, which are the quiet start and the first-order kick.

### The first-order kick

The other value, `"first_order"`, displaces the phases by the algorithm of Fawley ([](references.md#ref-fawley)), generalized to weights, and applies to any load. Each group of macroparticles that share their transverse coordinates is assigned one phasor for each harmonic, whose amplitude is determined by the number of electrons of the group, which is its charge divided by $e$:

$$
  \begin{aligned}
    \theta_{gj} &\mathrel{-}= \sum_{h=0}^{n_h-1} a_{gh}\,
       \sin\!\bigl((h{+}1)\theta_{gj} + \phi_{gh}\bigr),\\[4pt]
    a_{gh} &= \frac{2}{h+1}\sqrt{\frac{-\ln u_{gh}}{n_{g}}},
      \qquad n_{g} = \frac{w_g}{e},
  \end{aligned}
$$ (eq-fawley)

where $n_h = \lfloor (b-1)/2 \rfloor$, $u_{gh}$ is uniform in $(0,1)$, $\phi_{gh}$ is uniform in $[0,2\pi)$, and the [kicks](glossary.md#gl-kick-noise) are evaluated from the unperturbed phases. To first order in the kick this gives $\langle|b(h)|^2\rangle = 1/N_\lambda$ for any distribution of the weight among the groups. The kick is small when a group contains many electrons. [](shot-noise.md) derives in closed form the noise that a kick imposes at any number of electrons. For the kick of Eq. [](#eq-fawley) the noise power at the fundamental is low by approximately $2.4/n_g$ where a group contains many electrons, which is 0.95 of the shot-noise power at 50 electrons per group. The closed form of [](shot-noise.md#shot-kick) gives 0.16 at one electron.

This is the kick of Genesis4, and the checks that compare a load with that of Genesis4 select it. A group of less than one electron is kicked as one electron, and the report of the load lists the slices in which this occurred ([](shot-noise.md#shot-clamp)).

(pl-groups)=
### The groups

The groups are identified by the first of three rules that applies, and the rule that was used is recorded as `setup/noise_groups` in the [stats file](glossary.md#gl-statistics-file):

1. The beamlets, when this loader generated them.
2. Otherwise the particles whose five coordinates other than the phase are equal bit for bit. Copies generated by another loader have this property, and no other particles do.
3. Otherwise the particles in each [deposit](glossary.md#gl-deposit) cell of the field grid, which is the finest structure that the source term resolves.

### The quiet floor

Before the noise is imposed, the loader evaluates the quiet floor of each slice, $\max_h |b(h)|^2 N_\lambda$ over $h = 1,\ldots,b-1$, and refuses `shot_noise` for a slice in which it exceeds $0.01$, with the value in the message. The noise of a load that already has its own would otherwise be counted twice. A generated load without the quiet start has a value of the order of $N_\lambda/N_{\mathrm{part}}$, which is far above the threshold. A bunch that was made quiet outside the loader has a value at the level of rounding and is accepted.

(pl-fewer)=
### Fewer beamlets in a slice of little charge

Under the first-order kick a beamlet that represents few electrons receives less than the shot noise of those electrons, and [](shot-noise.md) derives the level as a function of the number of electrons per beamlet. For this reason, under the first-order kick and the quiet start, the sample mode does not load the same number of beamlets in every slice. A slice of $N_\lambda$ electrons is loaded with at most $N_\lambda/N_{\min}$ beamlets, where $N_{\min}$ is `beam.beamlet_min_electrons`, which is one electron by default, and with one beamlet when it contains fewer than $N_{\min}$ electrons. The charge of the slice is unchanged. The generator produces that number of base samples, and the resampler reduces its candidates to that number of seeds. [](shot-noise.md#shot-fewer) states the rule and the measurements on which the default was chosen.

The weight method does not use the key. Under it a beamlet in which no macroparticle is assigned an electron is omitted from the load, so that a slice of few electrons contains approximately one beamlet for each electron ([](shot-noise.md#shot-weight)). Under either method, and without the noise, a slice with no charge contains no particles in the sample mode. With the key at 0 and the first-order kick every slice contains `beam.n_particle` macroparticles, which a comparison with Genesis4 requires. The case of a beamlet of less than one electron, which [](shot-noise.md#shot-clamp) treats, remains only in a slice of less than one electron, in the keep mode, and with the key at 0.

:::{dropdown} Genesis4 source
`fel_fawley_noise` is transcribed from [`ShotNoise::applyShotNoise`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L34) and generalized by $n_{bl} = n_e/m \to w_g/e$ (identical for uniform weights). Genesis4 sets a group under one electron to one electron silently ([line 47](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L47)). Lucifer does the same and reports the slices it happened in. Genesis4 gives every slice the same beamlets, and `fel_beamlet_count` has no Genesis4 original. The generator, the resampler and the keep mode share the code (`fel_impose_noise`). Genesis4's `importdistribution` applied the noise whether or not its `shotnoise` flag said so. Here the two switches decide. The grouping of an arbitrary bunch (`fel_group_noise`), the weight load (`fel_weight_noise`) and the closed forms of [](shot-noise.md) have no Genesis4 original.
:::

The measured levels and how they are checked are in [](validation.md#val-shot-noise-under-weights), which compares the load with the closed forms.

(sec-import)=
## The resampler

In the sample mode a bunch from a file is resampled into the slices (`fel_import_mod`). The centers of the slices are at $(i-1)\Delta s$, and $n_{\mathrm{slice}} = \mathrm{nint}(\tau_{\mathrm{total}}/\Delta s)$ when the number is not given. For each slice, every particle within a sampling window of length $d s_{\mathrm{len}} = \texttt{slice\_width}\cdot\tau_{\mathrm{total}}$ centered on the slice, with strict inequalities at its ends, is a candidate, and the current of the slice is obtained from the same window,

$$
  I_i = \frac{c\sum_{j\in\mathrm{window}} w_j}{d s_{\mathrm{len}}},
$$ (eq-importcurrent)

the sum being over the weights of the candidates.

The set of candidates is reduced to $n_{\mathrm{part}}/b$ seeds by random deletion, or is enlarged by an interpolation in phase space. The interpolation normalizes the five coordinates $(\gamma,x,y,\hat p_x, \hat p_y)$ to zero mean and unit rms value, selects a parent at random, and finds the nearest of the original particles under a metric in which the weight of each coordinate is a new random number. The child is placed at the midpoint of the two, displaced in each coordinate by the difference multiplied by a number uniform in $[-1,1]$. The phase $\theta$ is then assigned uniformly over one spacing of the beamlet and repeated in $b$ copies, and the noise of [](#sec-noise) is imposed with $n_e = \mathrm{nint}(I\lambda_0\,\texttt{n\_wavelength}/(e c))$.

The quantities $p_x$ and $p_y$ of the file are slopes, and the conversion from slope to momentum, $\hat p = x'\gamma$, is made when a particle is copied into the set of candidates. No matching or centering transformation is applied, since a Bmad lattice contains its optics and `init_beam_distribution` generates matched bunches. A bunch of zero charge is refused.

The keys `beam.resample.n_particle_per_slice` and `beam.resample.beamlet_size` belong to the validation route of the resampler, `beam.resample.use_beam_init`. On that route the resampler is given a bunch that Bmad generated, so that it can be compared, on a known bunch, with Genesis4 importing the same file. On the ordinary route the two numbers are `beam.n_particle` and `beamlet_size`, which the generator also uses.

:::{dropdown} Genesis4 source
`fel_import_bunch` is transcribed from [`SDDSBeam::init`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/SDDSBeam.cpp#L58) (window, current, ordering), [`removeParticles`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/SDDSBeam.cpp#L606), [`addParticles`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/SDDSBeam.cpp#L471) and [`distance`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/SDDSBeam.cpp#L590). Slice centers are from [`Time::getPosition`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Main/GenTime.cpp#L98). The slice current is the weighted generalization of Genesis4's count$\cdot dQ$, and $\hat p = \gamma x'$ is kept as Genesis4 has it. Genesis4's `match` and `center` transforms have no counterpart here. Its `align*` inputs are parsed but unused in Genesis4, and are not transcribed. The keep mode and the generated bunch have no Genesis4 original.
:::

The measured levels and how they are checked are in [](validation.md#val-distribution-import-a-bunchstruct).

(pl-evidence)=
## Where the load is measured

This page contains no result of a run. The measurements that concern the load are recorded on the pages of their subjects.

| subject | record |
|---|---|
| the number of beamlets in a slice and the size of a beamlet | [](shot-noise.md#shot-count), [](shot-noise.md#shot-size) |
| the rule of fewer beamlets in a slice of little charge | [](shot-noise.md#shot-fewer) |
| the imposed noise compared with the closed forms | [](shot-noise.md#shot-measured) |
| the electrons mode: its memory, its run time and its comparison with the weight load | [](electrons-mode.md#el-measured) |
| the checks of the keep mode, of the resampler and of the imposed noise | [](validation.md#val-the-load-path), [](validation.md#val-distribution-import-a-bunchstruct), [](validation.md#val-shot-noise-under-weights) |
| the radiation that migration adds to a load of beamlets | [](migration.md#mig-measured) |
