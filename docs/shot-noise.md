---
title: "The shot-noise load: what a beamlet radiates"
short_title: Shot-noise load
---

A SASE run starts from the shot noise of the electron beam. The loader imposes that noise on a quiet load by Fawley's algorithm ([](references.md#ref-fawley)), as Genesis4 does in `ShotNoise::applyShotNoise`. The algorithm is exact to first order in the size of its phase kick. This page derives what the kick imposes to every order, in closed form, as a function of the number of real electrons one beamlet stands for.

The result is one number per harmonic, the ratio $R_h(N)$ of the noise power the load imposes to the shot-noise power of the real beam. The ratio is $1 - 2/N$ to first order for one harmonic, and it is 1 when a beamlet holds many electrons. It is 0.95 at 50 electrons a beamlet and 0.16 at one electron. The loader is held to the closed form by `check_shot_noise.py` at five electron counts ([](#shot-measured)).

The page uses the FEL quantities of [](fel-physics.md) and the words of the [glossary](glossary.md). The load itself, its switches and its grouping rules are in [](fel-physics.md#sec-noise).

(shot-counts)=
## Macroparticles and electrons

Four counts enter, and two of them are easily confused.

A slice is one longitudinal sample of the beam, a whole number of radiation wavelengths long. A slice of current $I$ and length $\Delta$ holds

$$
  N_\lambda = \frac{I\,\Delta}{e\,c}
$$ (eq-shot-nlambda)

real electrons. The loader represents the slice by $n_p$ macroparticles, `beam.n_particle`. Each macroparticle carries a weight, the charge it stands for.

The macroparticles are loaded in beamlets. A beamlet is $M$ macroparticles, `beamlet_size`, 8 by default. The macroparticles of one beamlet share their transverse coordinates and their energy. Their phases are spaced evenly over one radiation period. A slice holds $n_p/M$ beamlets.

The count that sets the noise is the number of real electrons one beamlet stands for,

$$
  N = \frac{M\,w}{e} ,
$$ (eq-shot-n)

with $w$ the weight of one of its macroparticles. This is `nbl` in `fel_fawley_noise`. With uniform weights $N = N_\lambda M / n_p$. $N$ is a count of electrons. It is *not* the number of macroparticles in the beamlet, which is $M$.

For example, take a slice one wavelength long at 1 Å and 3 kA. The slice holds $N_\lambda = 6246$ electrons. With 1024 macroparticles in beamlets of 8 the slice has 128 beamlets, and each beamlet stands for $N = 48.8$ electrons. A slice in the tail of the same bunch at $10^{-3}$ of the current holds 6 electrons, and each of its 128 beamlets stands for 0.049 of an electron.

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

This holds at every harmonic. For a slice it reads $\langle|b_h|^2\rangle = 1/N_\lambda$, which is the level the start-up power of a SASE run is proportional to. One electron has $|b_h| = 1$ at every harmonic.

A set of macroparticles has its own noise, $1/n_p$ for random phases, which is larger than the real noise by $N_\lambda/n_p$. The quiet start removes that noise and the kick puts the real level back.

(shot-kick)=
## The quiet beamlet and the kick

The $M$ phases of a quiet beamlet are

$$
  \theta_j = \theta_0 + \frac{2\pi j}{M}, \qquad j = 0 \ldots M-1 .
$$ (eq-shot-quiet)

The sum of $e^{-ih\theta_j}$ over $j$ is zero unless $h$ is a multiple of $M$. So the beamlet has no bunching at any harmonic below $M$, whatever $\theta_0$ is. A slice of such beamlets radiates nothing.

The load then moves every phase of the beamlet by

$$
  \theta_j \;\to\; \theta_j - \sum_{k=1}^{H} a_k \sin(k\theta_j + \phi_k),
  \qquad a_k = \frac{2}{k} \sqrt{\frac{x_k}{N}} ,
$$ (eq-shot-kick)

with $H = (M-1)/2$ harmonics in integer division, 3 for a beamlet of 8. Each harmonic of each beamlet has its own two draws. The phase $\phi_k$ is uniform over $2\pi$. The number $x_k = -\ln u_k$, with $u_k$ uniform on $(0,1)$, is exponentially distributed with mean 1. The sines are taken at the unperturbed phases. This is [](fel-physics.md#eq-fawley) for one beamlet.

:::{figure} generated/shot-noise/beamlet.svg
:name: fig-shot-beamlet
:alt: Two unit circles. On the left, eight points evenly spaced, with zero bunching at the center. On the right, the same eight points displaced along the circle toward one side by a sinusoidal kick, with a bunching phasor drawn from the center.
A beamlet of 8 on the unit circle, each point at $e^{-i\theta_j}$. Left: the quiet beamlet, whose points sum to zero. Right: the same beamlet after a kick at the fundamental of size $z = 0.7$ rad, drawn large. The open circles are the quiet phases. The points crowd toward one side, and their mean is the bunching phasor of length $J_1(z)$.
:::

(shot-one)=
## One harmonic

Take the kick at one harmonic first, the fundamental, and write $z = a_1 = 2\sqrt{x/N}$. The bunching of the beamlet at the fundamental after the kick is

$$
  b_1 = \frac{1}{M} \sum_{j} e^{-i\theta_j}\, e^{\,i z \sin(\theta_j + \phi)} .
$$ (eq-shot-b1)

The Jacobi-Anger expansion, $e^{iz\sin\psi} = \sum_n J_n(z)\, e^{in\psi}$, turns this into

$$
  b_1 = \sum_{n} J_n(z)\, e^{in\phi}\; \frac{1}{M} \sum_{j} e^{\,i(n-1)\theta_j} .
$$ (eq-shot-ja)

The even spacing makes the sum over $j$ zero unless $n - 1$ is a multiple of $M$. The term $n = 1$ survives, and the next terms are $n = 1 \pm M$, which are of order $z^{M-1}$. For a beamlet of 8 the larger of them is $1.7\times10^{-4}$ at $z = 2$ and $1.5\times10^{-6}$ at $z = 1$. So

$$
  |b_1| = J_1(z), \qquad z = 2\sqrt{x/N} .
$$ (eq-shot-j1)

For a small kick $J_1(z) \approx z/2 = \sqrt{x/N}$. The mean of $x$ is 1, so $\langle|b_1|^2\rangle = 1/N$, which is [](#eq-shot-real). That is the algorithm's design: the kick is sized so that its first-order bunching has the shot-noise mean, and the exponential draw gives $|b_1|^2$ the exponential distribution of real shot noise.

The kick is small only when $N$ is large. Its typical size is $2/\sqrt{N}$: 0.16 rad at 150 electrons, 0.47 rad at 18 and 2 rad at one. $J_1(z)$ falls below $z/2$ as $z$ grows and has its maximum of 0.58 at $z = 1.84$.

:::{figure} generated/shot-noise/kick.svg
:name: fig-shot-kick
:alt: The Bessel function J1 of z and the straight line z over 2 against the kick z from 0 to 4 radians, with three shaded distributions of the kick for 150, 18 and 1 electrons a beamlet. The two curves part near z of 0.5.
The bunching a beamlet takes from a kick of size $z$, $J_1(z)$, against the first-order value $z/2$. The shaded curves are the distributions of the kicks the load draws at 150, 18 and 1 electrons a beamlet, each scaled to the same height. At 150 electrons the kicks lie where the two curves agree. At one electron most kicks lie where $J_1$ is well below $z/2$.
:::

The mean over the draw of $x$ is Weber's second exponential integral ([](references.md#ref-watson)),

$$
  \langle |b_1|^2 \rangle = \int_0^\infty e^{-x}\, J_1^2\!\bigl(2\sqrt{x/N}\bigr)\, dx
     = e^{-2/N}\, I_1(2/N) ,
$$ (eq-shot-weber)

with $I_1$ the modified Bessel function. The ratio to the shot-noise value is

$$
  R(N) = N\, e^{-2/N}\, I_1(2/N) = 1 - \frac{2}{N} + \frac{5}{2N^2} - \ldots
$$ (eq-shot-r1)

The same holds for any single harmonic $h$ kicked alone, since the kick enters $b_h$ as $z = h\,a_h = 2\sqrt{x/N}$.

(shot-all)=
## The harmonics together

The loader kicks $H$ harmonics together, and each kick also changes the bunching at the other harmonics. The mean is found without the expansion, from the pairs of macroparticles.

From [](#eq-shot-b) and [](#eq-shot-kick),

$$
  |b_h|^2 = \frac{1}{M^2} \sum_{j,l} e^{-ih(\theta_j - \theta_l)}
     \prod_{k=1}^{H} \exp\!\Bigl( i h a_k \bigl[ \sin(k\theta_j + \phi_k) - \sin(k\theta_l + \phi_k) \bigr] \Bigr) .
$$ (eq-shot-pairs)

The difference of the two sines is $2 \sin(k\delta/2) \cos(\psi_k)$, with $\delta = \theta_j - \theta_l$ and $\psi_k$ a phase that is uniform because $\phi_k$ is. The mean of $e^{i\zeta\cos\psi}$ over a uniform $\psi$ is $J_0(\zeta)$. So the mean over $\phi_k$ of the factor $k$ is

$$
  J_0\!\Bigl( 2\sqrt{\beta_k x_k} \Bigr), \qquad
  \beta_k = \frac{4h^2}{k^2 N} \sin^2\frac{k\delta}{2} = \frac{2h^2}{k^2 N}\,(1 - \cos k\delta) .
$$ (eq-shot-j0)

The mean over the exponential draw follows from the series of $J_0$ term by term, and is Weber's first exponential integral:

$$
  \int_0^\infty e^{-x} J_0\!\bigl(2\sqrt{\beta x}\bigr)\, dx
    = \sum_{m} \frac{(-\beta)^m}{(m!)^2} \int_0^\infty e^{-x} x^m\, dx = e^{-\beta} .
$$ (eq-shot-weber0)

The draws of different harmonics are independent, so the product of the means is the mean of the product. The pair $(j, l)$ enters only through $\delta = 2\pi(j-l)/M$, so the double sum is $M$ times a single sum:

$$
  \langle |b_h|^2 \rangle = \frac{1}{M} \sum_{j=0}^{M-1} \cos\frac{2\pi j h}{M}\;
    \exp\!\Bigl( -\frac{2h^2}{N} \sum_{k=1}^{H} \frac{1 - \cos(2\pi j k/M)}{k^2} \Bigr) .
$$ (eq-shot-closed)

This is exact for the kick of [](#eq-shot-kick). It is a sum of $M$ terms, and `tests/scripts/shot_noise_ratio.py` evaluates it. The ratio the tables and figures below give is

$$
  R_h(N) = N\, \langle |b_h|^2 \rangle .
$$ (eq-shot-rh)

Three limits check the form. With one harmonic and dense phases the sum over $j$ is an integral over $\delta$, and the integral is $e^{-2/N} I_1(2/N)$ by the integral representation of $I_1$, which is [](#eq-shot-weber). To first order in $1/N$ only the term $k = h$ survives the sum over $j$, and it gives $1/N$ for every kicked harmonic. A harmonic above $H$ is not kicked and has no first-order term. Its bunching is of order $1/N^2$.

To second order the loader's beamlet of 8 with 3 harmonics has

$$
  R_1 \approx 1 - \frac{22}{9N}, \qquad R_2 \approx 1 - \frac{1.01}{N}, \qquad R_3 \approx 1 - \frac{2.00}{N} .
$$ (eq-shot-second)

The fundamental loses more than the $2/N$ of [](#eq-shot-r1), because the kicks at the second and third harmonic also spread its phases.

(shot-table)=
## The ratio against the electron count

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-ratio -->
:end-before: <!-- end tab-shot-ratio -->
```

:::{figure} generated/shot-noise/ratio.svg
:name: fig-shot-ratio
:alt: The ratio of imposed to physical noise power against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis. The curves rise from near zero below one electron through about 0.5 at four electrons to 1 above one hundred. Measured points with error bars lie on the curves.
The ratio $R_h(N)$ of [](#eq-shot-rh) for the loader's beamlet of 8, at the fundamental and the two harmonics it kicks. The dashed curve is one harmonic kicked alone, [](#eq-shot-r1), and the dotted curve its first-order form. The points are the loader's own levels with one standard deviation, from `check_shot_noise.py` ([](#shot-measured)). In the shaded band a beamlet holds less than one electron ([](#shot-clamp)).
:::

The ratio has no threshold. The deficit at the fundamental is close to $2.4/N$ above 20 electrons a beamlet. The worked slice of [](#shot-counts), at 48.8 electrons a beamlet, starts with 0.951 of the shot-noise power at the fundamental.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-floor -->
:end-before: <!-- end tab-shot-floor -->
```

$N$ follows from the slice and the load by [](#eq-shot-n). At a given current it rises with the slice length, `window.n_wavelength`, and falls with the macroparticles per slice, `beam.n_particle`.

A slice whose beamlets hold different counts, as with nonuniform weights, has the charge-weighted mean of their ratios:

$$
  \langle |b_h|^2 \rangle\, N_\lambda = \frac{\sum_g N_g\, R_h(N_g)}{\sum_g N_g} ,
$$ (eq-shot-slice)

with $N_g$ the electrons of beamlet $g$. The beamlets are kicked independently, so their bunchings add in the mean square with the squares of their charges.

(shot-clamp)=
## Below one electron

A beamlet may stand for less than one electron. The tail slice of [](#shot-counts) is such a load: 6 electrons over 128 beamlets. The kick of [](#eq-shot-kick) at $N < 1$ would be larger than at one electron. Genesis4 sets $N$ to 1 in the kick when it is below 1, and the loader does the same, counts the draws and warns:

```
<count> beamlet draws had fewer than one real electron
(nbl clamped to 1, as Genesis does silently). The noise level in those
beamlets is not physical; use fewer macroparticles or more charge.
```

The beamlet's bunching is then the bunching at one electron, $\langle|b_h|^2\rangle = R_h(1)$, whatever its charge. A slice of such beamlets has

$$
  R_h(N) = N\, R_h(1), \qquad N < 1 .
$$ (eq-shot-clamp)

$R_1(1)$ is 0.161 for the loader's beamlet. The tail slice above, at 0.049 of an electron a beamlet, starts with 0.008 of its shot-noise power. A real slice of 6 electrons has $\langle|b|^2\rangle = 1/6$. The loaded slice has $0.161/128 = 0.0013$.

The loader also takes a kick $a_k$ above $2\pi$ modulo $2\pi$. At one electron that is $e^{-\pi^2} = 5\times10^{-5}$ of the draws at the fundamental, and fewer at every larger count. The closed form leaves it out.

(shot-measured)=
## The loader against the closed form

`check_shot_noise.py` loads a uniform bunch of 256 slices with 15 seeds and reads the bunching of every slice at harmonics 1 to 3 from the particle dump. It runs at 146, 50, 18 and 4 electrons a beamlet and at a quarter of an electron, each with uniform weights and with beamlet weights alternating between 0.25 and 1.75 of the mean. The mean of $|b_h|^2 N_\lambda$ over the samples is held to [](#eq-shot-slice) within five standard deviations of the mean. The check runs in every pass of the validation suite, and its levels are in [](validation.md#val-levels-shot-noise). They are the points of [](#fig-shot-ratio).

(shot-scope)=
## Scope

The closed form [](#eq-shot-closed) is for a beamlet of evenly spaced phases kicked by `fel_fawley_noise`. That covers the loader's own beamlets in both load modes. It also covers copies made by another loader, which `fel_group_noise` kicks by the same rule when the copies are evenly spaced.

Where the groups are the occupants of a deposit cell, one group holds the copies of several particles at different base phases. The terms of [](#eq-shot-ja) that do not depend on the base phase are the same, so the dense form, [](#eq-shot-closed) at large $M$, applies with $N$ the electrons of the whole group. A group whose phases are not copies at even spacing is outside this page.

The page is about the load. What the FEL does with the noise, and how the start-up power depends on the transverse grid, is in [](startup-noise.md).
