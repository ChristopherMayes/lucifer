---
title: "The shot-noise load: what a beamlet radiates"
short_title: Shot-noise load
---

A SASE run starts from the [shot noise](glossary.md#gl-shot-noise) of the electron beam: $N$ electrons at random phases have a mean square bunching of $1/N$. Lucifer makes the [macroparticles](glossary.md#gl-macroparticle) of every [slice](glossary.md#gl-slice) before tracking begins. That is the [load](glossary.md#gl-load), and the code that makes it is the loader. The load is made quiet first, in [beamlets](glossary.md#gl-beamlet) of evenly spaced phases that have no bunching ([quiet start](glossary.md#gl-quiet-start)). The loader then imposes the shot noise, by `beam.shot_noise_method`:

- `"weight"`, the default, gives every macroparticle a whole number of electrons and moves no phase ([](#shot-weight)). The noise power it loads is the shot noise at every harmonic below the beamlet size, at any number of electrons, and at the wavelengths between the harmonics. The bunch's charge is held exactly.
- `"first_order"` is Fawley's algorithm ([](references.md#ref-fawley)), which shifts every phase by a small random amount, the [kick](glossary.md#gl-kick-noise) ([](#shot-kick)). It is Genesis4's load. It loads less than the shot noise where a beamlet stands for few electrons, 0.95 of it at 50 electrons a beamlet and 0.16 at one, and away from the harmonics it kicks it loads from 0.004 of the shot noise to 4 times it.

Take the default. `"first_order"` is for a comparison with Genesis4, and for a bunch whose beamlets another code made, which the weight load does not take ([](#shot-scope)).

:::{figure} generated/shot-noise/bunch.svg
:name: fig-shot-bunch
:alt: Two plots. Top: the noise power of each slice of a Gaussian bunch against its place along the bunch on a logarithmic axis of eleven decades. The points of the weight load lie on the line of the real beam down to a tenth of an electron. The points of the first-order kick fall below the line in the tails, by three decades at five rms lengths. Bottom: the bunching spectrum of the whole bunch against the wavenumber over the resonant one, from zero to 8.6, on a logarithmic axis. The weight load is level at one. The first-order kick falls to zero below the fundamental and is three to four above the fourth harmonic. The load with no noise is a row of narrow lines at every third of the resonant wavenumber, the first at one and the rest falling, and nothing between them.
A Gaussian bunch of 384 slices, each three wavelengths long, with 50 electrons a beamlet in its peak slice, loaded at 48 seeds under each method. Top: the mean of $|S_1|^2$ in each slice against the slice's place, where $S_1$ is the sum of $e^{-i\theta}$ over the slice's electrons ([](#eq-shot-s)). The line is the real beam's value, the slice's electrons $N_\lambda$. Bottom: the mean of $|S(k)|^2$ over the bunch's electrons $N$, from the particles' own places along the bunch, against the wavenumber over the resonant one, at the [window](glossary.md#gl-window)'s own frequencies. The real beam's value is 1 at every wavenumber, the black line. The gray curve is the load with no noise, whose lines are the steps of the current from slice to slice ([](#shot-spectrum)).
:::

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-bands -->
:end-before: <!-- end tab-shot-bands -->
```

The load itself, its switches and its grouping rules are in [](fel-physics.md#sec-noise).

(shot-counts)=
## Macroparticles and electrons

Four counts enter, and two of them are easily confused.

A slice is one longitudinal sample of the beam, a whole number of radiation wavelengths long, and its macroparticles stand over all of those wavelengths ([](fel-physics.md#sec-quiet)). A slice of current $I$ and length $\Delta$ holds

$$
  N_\lambda = \frac{I\,\Delta}{e\,c}
$$ (eq-shot-nlambda)

real electrons. The load represents the slice by $n_p$ macroparticles, `beam.n_particle`. Each macroparticle carries a weight, the charge it stands for.

The macroparticles are loaded in beamlets. A beamlet is $M$ macroparticles, `beamlet_size`, 8 by default. The macroparticles of one beamlet share their transverse coordinates and their energy. Their phases are spaced evenly over one radiation period,

$$
  \theta_j = \theta_0 + \frac{2\pi j}{M}, \qquad j = 0 \ldots M-1 ,
$$ (eq-shot-quiet)

each plus a whole number of $2\pi$ where a slice is several wavelengths long, which no harmonic sees. The sum of $e^{-ih\theta_j}$ over $j$ is zero unless $h$ is a multiple of $M$, so the beamlet has no bunching at any harmonic below $M$. A slice of $n_p/M$ such beamlets radiates nothing.

The count that sets the noise of a beamlet is the number of real electrons it stands for,

$$
  N = \frac{M\,w}{e} ,
$$ (eq-shot-n)

with $w$ the weight of one of its macroparticles in the quiet load. With uniform weights $N = N_\lambda M / n_p$. $N$ is a count of electrons. It is *not* the number of macroparticles in the beamlet, which is $M$. At a given current $N$ rises with the slice length, `window.n_wavelength`, and falls with `beam.n_particle`.

A slice one wavelength long at 1 Å and 3 kA holds $N_\lambda = 6246$ electrons. With 1024 macroparticles in beamlets of 8 the slice has 128 beamlets, and each beamlet stands for $N = 48.8$ electrons. A slice in the tail of the same bunch at $10^{-3}$ of the current holds 6 electrons, 0.049 of an electron to each of 128 beamlets.

(shot-real)=
## The shot noise of N electrons

The bunching of a set of $N$ electrons at harmonic $h$ is

$$
  b_h = \frac{1}{N} \sum_{j=1}^{N} e^{-i h \theta_j} ,
$$ (eq-shot-b)

with $\theta_j$ the ponderomotive phase of electron $j$. The electrons of a real beam arrive independently, so the phases are independent and uniform over $2\pi$. The cross terms of $|b_h|^2$ then average to zero and the $N$ diagonal terms remain:

$$
  \langle |b_h|^2 \rangle = \frac{1}{N} .
$$ (eq-shot-real)

This holds at every harmonic and at every wavelength between them. For a slice it reads $\langle|b_h|^2\rangle = 1/N_\lambda$, which the start-up power of a SASE run is proportional to.

A set of macroparticles has its own noise, $1/n_p$ for random phases, which is larger than the real noise by $N_\lambda/n_p$. The quiet start removes that noise and the load puts the real level back. The page's measure of a load is the ratio $R_h(N)$ of the noise power it imposes to the shot-noise power of the real beam.

(shot-weight)=
## The noise through the weights

The phases stay the quiet beamlet's, and the weights carry the noise. This is the load of McNeil, Poole and Robb ([](references.md#ref-mcneil)), put on a beamlet as Hwang and Qiang put it ([](references.md#ref-hwang-qiang-weight)).

A macroparticle of a beamlet of $N$ electrons stands for $N/M$ of them. Let it stand for a whole number $c_j$ of electrons instead, drawn as the electrons of a real beam fall: independently, each into any of the $M$ macroparticles with the same probability. Independent counts of that kind are Poisson numbers of mean and variance $N/M$. The beamlet's source at harmonic $h$ is the sum over its electrons,

$$
  S_h = \sum_{j=0}^{M-1} c_j\, e^{-ih\theta_j} .
$$ (eq-shot-s)

The mean of $S_h$ is zero for the evenly spaced phases at every $h$ that is not a multiple of $M$. Its mean square is the sum of the variances of the counts:

$$
  \langle |S_h|^2 \rangle = N,
  \qquad \langle |b_h|^2 \rangle = \frac{\langle |S_h|^2 \rangle}{N^2} = \frac{1}{N} .
$$ (eq-shot-weight)

This is [](#eq-shot-real), and no step of it is an approximation. The bunching is linear in the counts, so there is no order of a kick to stop at. It holds at any $N$, whole or not, and at every harmonic below $M$: 1 to 7 for the default beamlet. At one electron a beamlet most macroparticles draw no electron and a few draw one. The load is then the electrons themselves, each at one of $M$ phases.

### The bunch's charge, held

Independent counts leave the bunch's charge to chance. The loader holds it. A bunch of charge $Q$ is counted as $n$ electrons, the whole number nearest $Q/e$ and at least one. The counts of all the macroparticles of the bunch are then one multinomial draw of $n$ electrons, each macroparticle with a probability $p_j$ in proportion to the weight the quiet load gave it. The variance of $c_j$ is $n p_j (1 - p_j)$ and the covariance of two counts is $-n p_j p_k$, so for any set of macroparticles, a slice or a beamlet,

$$
  \langle |S_h|^2 \rangle = n \sum_j p_j \;-\; n\, \Bigl| \sum_j p_j\, e^{-ih\theta_j} \Bigr|^2 .
$$ (eq-shot-multinomial)

The second term is the square of the mean source, which is zero beamlet by beamlet. What is left is the mean count of the set. [](#eq-shot-weight) therefore holds with the bunch's charge held, for a slice of any charge and for any weights across the beamlets.

The charge of one counted electron is $Q/n$, and the loaded charge is $Q$ exactly. The noise power goes as the square of that charge times $n$, so

$$
  R = \frac{Q/e}{n}
$$ (eq-shot-weight-ratio)

at every harmonic and in every slice. $R$ is within $1/(2n)$ of 1, which is $10^{-9}$ for the $5\times10^8$ electrons of FLASH1's window. A bunch of less than half an electron is counted as one.

A slice's count is a binomial number of mean $N_\lambda$ and variance $N_\lambda (1 - N_\lambda/n)$. The charge of a slice thus varies from shot to shot by $1/\sqrt{N_\lambda}$ of itself, as the real beam's does: 0.7% for the 18740 electrons of a slice of `examples/sase`. That variation is the noise the load carries at long wavelengths, below the fundamental. The current profile of the load is the bunch's in the mean.

`fel_weight_noise` takes the draw in order. A slice takes a binomial share of the electrons not yet given, and each macroparticle a binomial share of what its slice has left. `fel_ran_binomial` is exact in distribution: Hörmann's transformed rejection ([](references.md#ref-hormann)) where ten or more successes are expected, and a sum of the probabilities upward from zero successes below. One `ran_seed` can be a different shot under another compiler or build, since a probability that differs in its last digit can change one count and with it every count after.

A beamlet none of whose copies drew an electron is not loaded, so a slice of few electrons loads about one beamlet for each. A copy that drew none, in a beamlet that holds charge, stays at zero weight, so a loaded slice is still whole beamlets. A slice that drew no electron has no particles.

### The keep mode

In the keep mode a beamlet is a particle of the bunch and its copies, and every particle keeps its own charge. A particle of $N$ electrons is counted as $k = \lfloor N \rfloor$ electrons or as $k + 1$, and the charge of one counted electron is the particle's charge over its count. Its noise power goes as its charge squared over the count. The count is $k$ with the probability

$$
  P(k) = \frac{k\,(k + 1 - N)}{N} ,
$$ (eq-shot-kept)

which makes the mean of one over the count $1/N$. The ratio is then 1 for every particle of one electron or more. The count falls on the $M$ copies with equal probability, so at least one copy holds charge. A particle under one electron is counted as one. One of its copies holds its whole charge, and the ratio is $N$.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-methods -->
:end-before: <!-- end tab-shot-methods -->
```

:::{figure} generated/shot-noise/weight.svg
:name: fig-shot-weight
:alt: The ratio of imposed to physical noise power at the fundamental against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis. The first-order kick rises from near zero below one electron to 1 above one hundred. The weight load is a level line at 1, with measured points at five electron counts lying on it.
The ratio $R_1(N)$ at the fundamental against the electrons a beamlet holds, under the two methods. The points are the weight load's levels at harmonics 1 to 7 measured on loaded beams by `check_shot_noise.py`, with one standard deviation, set side by side at each count.
:::

(shot-spectrum)=
### Along the bunch and across wavelength

:::{figure} generated/shot-noise/slices.svg
:name: fig-shot-slices
:alt: The noise power of each slice of a Gaussian bunch against its place along the bunch in nanometers, on a logarithmic axis from 1e-7 to 1e4. The weight load's points lie on the real beam's line. The first-order kick's points fall away from it in the tails.
The top plot of [](#fig-shot-bunch): 384 slices of three wavelengths, 48 seeds. A slice of the weight load reads its electrons down to the slices that drew one electron in 48 shots. Under the first-order kick a slice of the tails reads its electrons times the ratio of [](#shot-clamp).
:::

The weight load has the real beam's noise in every slice that drew an electron. The first-order kick is below it by 10 at four rms lengths from the peak and by $10^3$ at five.

:::{figure} generated/shot-noise/spectrum.svg
:name: fig-shot-spectrum
:alt: The bunching spectrum of the whole bunch from zero to 8.6 resonant wavenumbers on a logarithmic axis. The weight load fluctuates about 1 at every wavenumber. The first-order kick rises from zero to 1 at the fundamental, dips to 0.04 below the second harmonic and is near 4 above the fifth. The load with no noise is narrow lines at the thirds of the resonant wavenumber.
The bottom plot of [](#fig-shot-bunch), at the window's own frequencies with no binning. The fluctuation about the mean is that of 48 shots.
:::

:::{figure} generated/shot-noise/spectrum-low.svg
:name: fig-shot-spectrum-low
:alt: The same spectrum from zero to 2.2 resonant wavenumbers on a linear axis, with the mean in bins of 0.05 drawn over the unbinned curve. The weight load's mean stays at 1. The first-order kick's rises from zero to 1 at the fundamental, falls to 0.04 at 1.33 and peaks near 2 at 1.7.
The spectrum below 2.2 resonant wavenumbers on a linear axis. The marked curves are the mean in bins of 0.05, over the unbinned curves.
:::

The weight load's spectrum is white: the table under [](#fig-shot-bunch) gives its mean in six bands. The first-order kick's is not. Below the fundamental it falls to zero, because the charge of a slice does not vary: 0.004 of the shot noise below $0.3\,k_r$. It is 0.87 in the band about the fundamental, 0.04 at $1.33\,k_r$, 2 at $1.7\,k_r$ and 3 to 4.4 above the fourth harmonic, from the kick's higher orders.

With no noise the spectrum is below $10^{-14}$ between the multiples of $k_r/3$ and zero at them, since a slice's macroparticles stand on one even grid over the slice ([](fel-physics.md#sec-quiet)). Beside multiple $m$ of $k_r/3$ it has a pair of lines $0.003\,k_r$ apart, 1.2 high at $m = 1$ and falling as $1/m^2$, 0.13 at the fundamental. They are the steps of the current from one slice to the next, a slice being three wavelengths of one current, and the two loads with noise carry them too.

`tests/scripts/shot_noise_monte_carlo.py` draws the counts with no run of the program, in three forms: independent Poisson counts on one beamlet, a slice of 16 beamlets with its total held, and the keep mode's particle. At 150, 18, 4, 3.76, 1 and 0.25 electrons a beamlet, with $10^8$ draws of the first and third and $1.25\times10^7$ of the second, the means at harmonics 1 to 7 agree with the forms within 2.7 standard errors of the mean at all 126 points (2026-11-11). The standard error is $1.4\times10^{-4}$ to $2.8\times10^{-4}$ and the largest difference is $5\times10^{-4}$. At 3.76 electrons a beamlet the slice of 16 beamlets holds 60.16 electrons counted as 60, and its mean is 1.0024 to 1.0031 where [](#eq-shot-weight-ratio) gives 1.0027.

| load | electrons a beamlet | loaded power over shot-noise power, `"first_order"` | `"weight"` | start-up power, `"weight"` over `"first_order"` | the levels' ratio | exit pulse energy, `"first_order"` | `"weight"` | run time, `"first_order"` | `"weight"` |
|---|---|---|---|---|---|---|---|---|---|
| import, kept bunch | 3.76 | 0.5452 | 1.0000 | 1.85 ± 0.18 | 1.83 | 8.10 ± 0.27 pJ | 14.82 ± 0.67 pJ | 53.9 s | 53.0 s |
| sase | 73.2 | 0.9673 | 1.0000 | 1.11 ± 0.09 | 1.034 | 3.40 ± 0.91 nJ | 2.41 ± 0.55 nJ | 40.2 s | 40.0 s |
| FLASH1 | 16700 at the peak | 0.9995 | 1.0000 | 1.035 ± 0.024 | 1.0005 | 48.4 ± 0.9 µJ | 48.9 ± 1.9 µJ | 19.1 s | 19.0 s |
| LCLS | 1250 at the peak | 0.9981 | 1.0000 | 1.011 ± 0.011 | 1.0019 | 278.9 ± 1.2 µJ | 277.4 ± 1.7 µJ | 63.1 s | 63.1 s |

The table is four loads tracked at `run.ran_seed` 1 to 6 under each method on 2026-11-11, as means over the six shots with their standard errors, the start-up power summed over the window one gain length into the line and the run time a median on 12 threads. The weight load starts each load at the level the closed forms give within the spread of six shots and in the same run time, with the loaded charge the input file's to rounding, so it is the default.

(shot-kick)=
## The first-order kick

`beam.shot_noise_method = "first_order"` moves every phase of a beamlet by

$$
  \theta_j \;\to\; \theta_j - \sum_{k=1}^{H} a_k \sin(k\theta_j + \phi_k),
  \qquad a_k = \frac{z_k}{k}, \qquad z_k = 2\sqrt{\frac{x_k}{N}} ,
$$ (eq-shot-kick)

with $H = (M-1)/2$ harmonics in integer division, 3 for a beamlet of 8. The sines are taken at the unperturbed phases. Each harmonic of each beamlet has its own two draws. The phase $\phi_k$ is uniform over $2\pi$. The number $x_k = -\ln u_k$, with $u_k$ uniform on $(0,1)$, is exponentially distributed with mean 1. The draw asks the beamlet for a bunching of $\sqrt{x_k/N}$ at harmonic $k$, whose square has the mean $1/N$ and the exponential distribution of real shot noise. $z_k$ is the size of the kick at its own harmonic. This is [](fel-physics.md#eq-fawley) for one beamlet.

:::{figure} generated/shot-noise/beamlet.svg
:name: fig-shot-beamlet
:alt: Two unit circles. On the left, eight points evenly spaced, with zero bunching at the center. On the right, the same eight points displaced along the circle toward one side by a sinusoidal kick, with a bunching phasor drawn from the center.
A beamlet of 8 on the unit circle, each point at $e^{-i\theta_j}$. Left: the quiet beamlet, whose points sum to zero. Right: the same beamlet after a kick at the fundamental of size $z = 0.7$ rad, drawn large. The open circles are the quiet phases. The mean of the points is the bunching phasor of length $J_1(z)$.
:::

### One harmonic

Take the kick at the fundamental alone and write $z = z_1$. The bunching of the beamlet after the kick is

$$
  b_1 = \frac{1}{M} \sum_{j} e^{-i\theta_j}\, e^{\,i z \sin(\theta_j + \phi)}
      = \sum_{n} J_n(z)\, e^{in\phi}\; \frac{1}{M} \sum_{j} e^{\,i(n-1)\theta_j} ,
$$ (eq-shot-ja)

by the Jacobi-Anger expansion. The even spacing makes the sum over $j$ zero unless $n - 1$ is a multiple of $M$. The term $n = 1$ survives, and the next are of order $z^{M-1}$: $1.7\times10^{-4}$ at $z = 2$ for a beamlet of 8. So

$$
  |b_1| = J_1(z) .
$$ (eq-shot-j1)

For a small kick $J_1(z) \approx z/2 = \sqrt{x/N}$, the bunching asked for. That is the algorithm's design: the kick is sized so that its first-order bunching has the shot-noise mean. The kick is small only when $N$ is large. Its typical size is $2/\sqrt{N}$: 0.16 rad at 150 electrons, 0.47 rad at 18 and 2 rad at one. $J_1(z)$ falls below $z/2$ as $z$ grows and has its maximum of 0.58 at $z = 1.84$.

:::{figure} generated/shot-noise/kick.svg
:name: fig-shot-kick
:alt: The Bessel function J1 of z and the straight line z over 2 against the kick z from 0 to 4 radians, with three shaded distributions of the kick for 150, 18 and 1 electrons a beamlet. The two curves part near z of 0.5.
The bunching $|b_1|$ a beamlet takes from a kick of size $z$, $J_1(z)$, against the first-order value $z/2$. The shaded curves are the distributions of the kicks the load draws at 150, 18 and 1 electrons a beamlet, each scaled to the same height.
:::

The mean over the draw of $x$ is Weber's second exponential integral ([](references.md#ref-watson)),

$$
  \langle |b_1|^2 \rangle = \int_0^\infty e^{-x}\, J_1^2\!\bigl(2\sqrt{x/N}\bigr)\, dx
     = e^{-2/N}\, I_1(2/N) ,
  \qquad R(N) = N\, e^{-2/N}\, I_1(2/N) = 1 - \frac{2}{N} + \frac{5}{2N^2} - \ldots
$$ (eq-shot-r1)

with $I_1$ the modified Bessel function. The same holds for any single harmonic kicked alone.

### The harmonics together

The kick of [](#eq-shot-kick) sums $H$ harmonics, and each kick also changes the bunching at the other harmonics. The mean is found from the pairs of macroparticles:

$$
  |b_h|^2 = \frac{1}{M^2} \sum_{j,l} e^{-ih(\theta_j - \theta_l)}
     \prod_{k=1}^{H} \exp\!\Bigl( i h a_k \bigl[ \sin(k\theta_j + \phi_k) - \sin(k\theta_l + \phi_k) \bigr] \Bigr) .
$$ (eq-shot-pairs)

The difference of the two sines is $2 \sin(k\delta/2) \cos(\psi_k)$, with $\delta = \theta_j - \theta_l$ and $\psi_k$ uniform. The mean of $e^{i\zeta\cos\psi}$ over a uniform $\psi$ is $J_0(\zeta)$, so the mean over $\phi_k$ of the factor $k$ is $J_0(2\sqrt{\beta_k x_k})$ with $\beta_k = (2h^2/k^2 N)(1 - \cos k\delta)$. The mean over the exponential draw follows from the series of $J_0$ term by term, Weber's first exponential integral:

$$
  \int_0^\infty e^{-x} J_0\!\bigl(2\sqrt{\beta x}\bigr)\, dx = e^{-\beta} .
$$ (eq-shot-weber0)

The draws of different harmonics are independent, and the pair $(j, l)$ enters only through $\delta = 2\pi(j-l)/M$, so

$$
  \langle |b_h|^2 \rangle = \frac{1}{M} \sum_{j=0}^{M-1} \cos\frac{2\pi j h}{M}\;
    \exp\!\Bigl( -\frac{2h^2}{N} \sum_{k=1}^{H} \frac{1 - \cos(2\pi j k/M)}{k^2} \Bigr) ,
  \qquad R_h(N) = N\, \langle |b_h|^2 \rangle .
$$ (eq-shot-rh)

This is exact for the kick of [](#eq-shot-kick). It is a sum of $M$ terms, and `fel_shot_ratio` and `tests/scripts/shot_noise_ratio.py` evaluate it. With one harmonic and dense phases it is [](#eq-shot-r1). A harmonic above $H$ is not kicked and has a bunching of order $1/N^2$. To second order the default beamlet has

$$
  R_1 \approx 1 - \frac{22}{9N}, \qquad R_2 \approx 1 - \frac{1.01}{N}, \qquad R_3 \approx 1 - \frac{2.00}{N} .
$$ (eq-shot-second)

The fundamental loses more than the $2/N$ of [](#eq-shot-r1), because the kicks at the second and third harmonic also spread its phases.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-ratio -->
:end-before: <!-- end tab-shot-ratio -->
```

:::{figure} generated/shot-noise/ratio.svg
:name: fig-shot-ratio
:alt: The ratio of imposed to physical noise power against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis. The curves rise from near zero below one electron through about 0.5 at four electrons to 1 above one hundred. Measured points with error bars lie on the curves.
The ratio $R_h(N)$ of [](#eq-shot-rh) for the default beamlet at the fundamental and the two harmonics above it. The dashed curve is one harmonic kicked alone, [](#eq-shot-r1), and the dotted curve $1 - 2/N$. The points are the levels measured on loaded beams, with one standard deviation, from `check_shot_noise.py`. In the shaded band a beamlet holds less than one electron ([](#shot-clamp)).
:::

The ratio has no threshold. The deficit at the fundamental is close to $2.4/N$ above 20 electrons a beamlet. The worked slice of [](#shot-counts), at 48.8 electrons a beamlet, starts with 0.951 of the shot-noise power.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-floor -->
:end-before: <!-- end tab-shot-floor -->
```

A slice whose beamlets hold different counts has the charge-weighted mean of their ratios,

$$
  \langle |b_h|^2 \rangle\, N_\lambda = \frac{\sum_g N_g\, R_h(N_g)}{\sum_g N_g} ,
$$ (eq-shot-slice)

with $N_g$ the electrons of beamlet $g$, since the beamlets are kicked independently.

(shot-clamp)=
### Below one electron

The kick at $N < 1$ would be larger than at one electron, and the loader sets $N$ to 1 in the kick when it is below 1. The beamlet's bunching is then the bunching at one electron, whatever its charge, and

$$
  R_h(N) = N\, R_h(1), \qquad N < 1 .
$$ (eq-shot-clamp)

$R_1(1)$ is 0.161 for the default beamlet. The tail slice of [](#shot-counts), 6 electrons over 128 beamlets, would start with 0.008 of its shot-noise power: $\langle|b|^2\rangle = 0.0013$ where the real slice has $1/6$.

Three loads have such a beamlet: a slice that holds less than one electron in all, which has one beamlet under [](#eq-shot-fewer), the keep mode, and the sample mode with `beam.beamlet_min_electrons = 0`. The load's report counts the slices of the first kind. It warns of the slices that hold one electron or more and have a group under one, with their level by [](#eq-shot-slice):

```
<count> slices have a noise group of less than one real electron.
They hold <share> of the charge that carries noise, and the noise power loaded
in them is <ratio> of their shot-noise power at the fundamental.
```

The loader also takes a kick $a_k$ above $2\pi$ modulo $2\pi$. At one electron that is $e^{-\pi^2} = 5\times10^{-5}$ of the draws at the fundamental. The closed form leaves it out.

:::{dropdown} Genesis4 source
`fel_fawley_noise` is transcribed from [`ShotNoise::applyShotNoise`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L34), with the electron count of the kick taken from the group's own charge. Genesis4 sets $N$ to 1 below one electron without a message ([line 47](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L47)), which the warning's text refers to. The closed forms of this page and the weight load, `fel_weight_noise`, have no Genesis4 original.
:::

(shot-fewer)=
### Fewer beamlets in a slice of little charge

Under the first-order kick with the quiet start, the sample mode loads a slice of $N_\lambda$ electrons with

$$
  n_b = \min\!\left(\frac{n_p}{M},\ \max\!\left(1,\ \left\lfloor \frac{N_\lambda}{N_{\min}} \right\rfloor\right)\right)
$$ (eq-shot-fewer)

beamlets, where $N_{\min}$ is `beam.beamlet_min_electrons`, 1 by default (`fel_beamlet_count`). The slice keeps its charge: the weight of a macroparticle is the slice charge over $n_b M$. A slice of $N_{\min}$ electrons or more then has no beamlet under $N_{\min}$ electrons, and a slice of fewer has one beamlet. The generator draws $n_b$ base samples for the slice, and the [resampler](glossary.md#gl-resampler) reduces the slice's candidates to $n_b$ seeds by the random deletion that reduces them to $n_p/M$. The keep mode keeps every particle of the bunch, and the weight load does not read the key. In every load of the sample mode a slice with no charge has no particles, unless the key is 0.

At the default the tail slice of [](#shot-counts) has 6 beamlets of one electron and starts with 0.161 of its shot-noise power, where 128 beamlets gave it 0.008. A larger $N_{\min}$ raises the level in the slices it acts on, to 0.56 at 4 electrons, 0.87 at 18 and 0.95 at 50, and acts on more slices: $N_{\min} = 18$ takes beamlets from the worked core slice at 2784 macroparticles a slice or more, and $N_{\min} = 50$ at 1000. FLASH1, a Gaussian bunch of 351 slices of 512 beamlets, was loaded and then tracked at six seeds at each value on 2026-11-10:

| $N_{\min}$ | slices with fewer beamlets | their share of the charge | macroparticles | loaded power over shot-noise power, window | exit pulse energy, six seeds |
|---|---|---|---|---|---|
| 0 | 0 | 0 | 1437696 | 0.9995 | 48.2 ± 2.7 µJ |
| 1 | 146 | 1.0e-5 | 880832 | 0.9995 | 49.9 ± 1.9 µJ |
| 4 | 162 | 4.7e-5 | 821904 | 0.9995 | 48.8 ± 2.0 µJ |
| 18 | 180 | 2.3e-4 | 752544 | 0.9995 | 48.2 ± 1.5 µJ |
| 50 | 192 | 6.2e-4 | 701472 | 0.9996 | 45.4 ± 1.5 µJ |

The values differ by less than the shots do, and a Gaussian bunch through one undulator at twelve seeds starts with 0.99 ± 0.06 of the flat window's power at 0 and 0.98 ± 0.04 at 1 ([](validation.md#val-shot-noise-under-weights)). The measurements therefore do not choose, and the default is the one value that takes a beamlet from a slice only when the slice has no electron for it.

A slice with fewer beamlets has fewer samples of the transverse distribution and of the energy spread. The header's `Load` line states the beamlets of a slice that was not reduced and the fewest any slice with particles has. The [device](glossary.md#gl-device) and the [FP32](glossary.md#gl-fp32) [twin](glossary.md#gl-twin) take each slice's own count, and the twin's residual guard is read on the slices of 256 particles or more ([](performance.md)). A load with a reduced slice is a different shot from the load with the key at 0 at one `ran_seed`, and a load with none is the same bits. `beam.beamlet_min_electrons = 0` gives every slice $n_p/M$ beamlets, as in Genesis4, and the [tiers](glossary.md#gl-tier) that compare a load with Genesis4's set it.

(shot-measured)=
## The load, measured against the closed forms

`check_shot_noise.py` loads a uniform bunch of 256 slices with 15 seeds and reads the bunching of every slice from the particle dump. It runs in every pass of the [validation suite](glossary.md#gl-validation-suite), and its levels are in [](validation.md#val-levels-shot-noise).

- The weight load, at 146, 18, 4 and 1 electrons a beamlet and at a quarter of an electron, is held to [](#eq-shot-weight-ratio) at each of harmonics 1 to 7 within five standard deviations of the mean. The loaded charge is held to the input file's to $10^{-12}$, every weight to a whole number of the counted electron's charge, and the variance of the slices' counts to the binomial variance. A slice that drew no electron must have no particles, and the keep mode every particle of the bunch at its own charge.
- The first-order kick, at 146, 50, 18 and 4 electrons a beamlet and at a quarter of an electron, with uniform weights and with beamlet weights alternating between 0.25 and 1.75 of the mean, is held at harmonics 1 to 3 to [](#eq-shot-slice), and the level the loader prints to the four decimals printed. A Gaussian bunch is held to the charge and the beamlets of [](#eq-shot-fewer), and the keep mode and two loads whose beamlets were made outside the loader to the same form.
- The fill of a slice of 3 and of 12 wavelengths ([](fel-physics.md#sec-quiet)). With no noise the sum of $w\,e^{-i(m/n)\theta}$ over a slice, at the $m$ that are not multiples of the slice's $n$ wavelengths, is held to $10^{-12}$ of the slice's charge. Under the weight load its square is held to the shot noise, and the first-order kick's is recorded: 0.36 of the shot noise at 3 wavelengths and 0.46 at 12.

The points of [](#fig-shot-ratio) and [](#fig-shot-weight) are these levels. The record behind [](#fig-shot-bunch) is written by `tests/scripts/shot_noise_spectrum.py`, which refuses a record whose weight load is more than 0.05 from 1 in any band.

(shot-scope)=
## Scope

The weight load needs the beamlets this loader makes, in either load mode. A load without the quiet start is refused under `"weight"`, with `"first_order"` in the message.

The closed form [](#eq-shot-rh) is for a beamlet of evenly spaced phases kicked by `fel_fawley_noise`. That covers the beamlets Lucifer's loader makes in both load modes, and evenly spaced copies a bunch file brings from another code's loader, which `fel_group_noise` kicks by the same rule. Where the groups are the occupants of a deposit cell, the dense form, [](#eq-shot-rh) at large $M$, applies with $N$ the electrons of the whole group. A group whose phases are not copies at even spacing is outside this page.

The page is about the load. What the FEL does with the noise, and how the start-up power depends on the transverse grid, is in [](startup-noise.md).
