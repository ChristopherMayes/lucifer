---
title: "The shot-noise load: what a beamlet radiates"
short_title: Shot-noise load
---

A SASE run starts from the [shot noise](glossary.md#gl-shot-noise) of the electron beam: $N$ electrons at random phases have a mean square bunching of $1/N$. Lucifer makes the [macroparticles](glossary.md#gl-macroparticle) of every [slice](glossary.md#gl-slice) before tracking begins. That is the [load](glossary.md#gl-load), and the code that makes it is the loader. The load is made quiet first, in [beamlets](glossary.md#gl-beamlet) of evenly spaced phases that have no bunching ([quiet start](glossary.md#gl-quiet-start)). The loader then imposes the shot noise in one of three ways, `beam.shot_noise_method`. Two of them are Fawley's algorithm ([](references.md#ref-fawley)), which shifts every phase by a small random amount, the [kick](glossary.md#gl-kick-noise). The algorithm sizes the kick to first order in its size. This page derives what a kick imposes to every order, in closed form, as a function of the number of real electrons one beamlet stands for, and then a kick sized from the bunching it gives. The third way moves no phase and imposes the noise through the weights, which is exact at every electron count and is the default.

The result is one number per harmonic, the ratio $R_h(N)$ of the noise power the load imposes to the shot-noise power of the real beam. Under the first-order kick the ratio is $1 - 2/N$ to first order for one harmonic. It is 0.95 at 50 electrons a beamlet and 0.16 at one electron. Lucifer sizes the kick from the Bessel function $J_1$, which is the bunching a kick gives ([](#shot-bessel)). The ratio at the fundamental is then $1 - 4/(9N)$ to first order: 0.99 at 50 electrons a beamlet, 0.91 at 10 and 0.22 at one. This is `beam.shot_noise_method = "bessel"`. `"first_order"` gives the first-order kick, which is Genesis4's. `"weight"`, the default, gives every macroparticle a whole number of electrons with the bunch's charge held ([](#shot-weight)). Its ratio is 1 at every harmonic below the beamlet size, at any number of electrons. The load is held to the three closed forms by `check_shot_noise.py` ([](#shot-measured)). Under the kicks a slice of little charge is loaded with fewer beamlets, so that none stands for less than one electron ([](#shot-fewer)).

The page uses the FEL quantities of [](fel-physics.md) and the words of the [glossary](glossary.md). The load itself, its switches and its grouping rules are in [](fel-physics.md#sec-noise).

(shot-counts)=
## Macroparticles and electrons

Four counts enter, and two of them are easily confused.

A slice is one longitudinal sample of the beam, a whole number of radiation wavelengths long, and its macroparticles stand over all of those wavelengths ([](fel-physics.md#sec-quiet)). A slice of current $I$ and length $\Delta$ holds

$$
  N_\lambda = \frac{I\,\Delta}{e\,c}
$$ (eq-shot-nlambda)

real electrons. The load represents the slice by $n_p$ macroparticles, `beam.n_particle`. Each macroparticle carries a weight, the charge it stands for.

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
  \theta_j = \theta_0 + \frac{2\pi j}{M}, \qquad j = 0 \ldots M-1 ,
$$ (eq-shot-quiet)

each plus a whole number of $2\pi$ where a slice is several wavelengths long, which no harmonic sees.

The sum of $e^{-ih\theta_j}$ over $j$ is zero unless $h$ is a multiple of $M$. So the beamlet has no bunching at any harmonic below $M$, whatever $\theta_0$ is. A slice of such beamlets radiates nothing.

The load then moves every phase of the beamlet by

$$
  \theta_j \;\to\; \theta_j - \sum_{k=1}^{H} a_k \sin(k\theta_j + \phi_k),
  \qquad a_k = \frac{z_k}{k} ,
$$ (eq-shot-kick)

with $H = (M-1)/2$ harmonics in integer division, 3 for a beamlet of 8. The sines are taken at the unperturbed phases. $z_k$ is the size of the kick at its own harmonic, since the kick enters the phase $k\theta_j$ as $k\,a_k$.

Each harmonic of each beamlet has its own two draws. The phase $\phi_k$ is uniform over $2\pi$. The number $x_k = -\ln u_k$, with $u_k$ uniform on $(0,1)$, is exponentially distributed with mean 1. The draw asks the beamlet for a bunching of $\sqrt{x_k/N}$ at harmonic $k$. Its square has the mean $1/N$ and the exponential distribution of real shot noise.

The two kicks of this page differ in how $z_k$ follows from the bunching asked for. Fawley's algorithm takes

$$
  z_k = 2\sqrt{\frac{x_k}{N}} ,
$$ (eq-shot-first)

the first-order kick. This is [](fel-physics.md#eq-fawley) for one beamlet, and [](#shot-one) to [](#shot-clamp) are about it. The kick sized from $J_1$ is in [](#shot-bessel).

:::{figure} generated/shot-noise/beamlet.svg
:name: fig-shot-beamlet
:alt: Two unit circles. On the left, eight points evenly spaced, with zero bunching at the center. On the right, the same eight points displaced along the circle toward one side by a sinusoidal kick, with a bunching phasor drawn from the center.
A beamlet of 8 on the unit circle, each point at $e^{-i\theta_j}$. Left: the quiet beamlet, whose points sum to zero. Right: the same beamlet after a kick at the fundamental of size $z = 0.7$ rad, drawn large. The open circles are the quiet phases. The points crowd toward one side, and their mean is the bunching phasor of length $J_1(z)$.
:::

(shot-one)=
## One harmonic

Take the kick at one harmonic first, the fundamental, and write $z = z_1$. The bunching of the beamlet at the fundamental after the kick is

$$
  b_1 = \frac{1}{M} \sum_{j} e^{-i\theta_j}\, e^{\,i z \sin(\theta_j + \phi)} .
$$ (eq-shot-b1)

The Jacobi-Anger expansion, $e^{iz\sin\psi} = \sum_n J_n(z)\, e^{in\psi}$, turns this into

$$
  b_1 = \sum_{n} J_n(z)\, e^{in\phi}\; \frac{1}{M} \sum_{j} e^{\,i(n-1)\theta_j} .
$$ (eq-shot-ja)

The even spacing makes the sum over $j$ zero unless $n - 1$ is a multiple of $M$. The term $n = 1$ survives, and the next terms are $n = 1 \pm M$, which are of order $z^{M-1}$. For a beamlet of 8 the larger of them is $1.7\times10^{-4}$ at $z = 2$ and $1.5\times10^{-6}$ at $z = 1$. So

$$
  |b_1| = J_1(z) .
$$ (eq-shot-j1)

This holds for a kick of any size. For a small kick $J_1(z) \approx z/2$, and the first-order kick of [](#eq-shot-first) has $z/2 = \sqrt{x/N}$, the bunching asked for. The mean of $x$ is 1, so $\langle|b_1|^2\rangle = 1/N$, which is [](#eq-shot-real). That is the algorithm's design: the kick is sized so that its first-order bunching has the shot-noise mean.

The first-order kick is small only when $N$ is large. Its typical size is $2/\sqrt{N}$: 0.16 rad at 150 electrons, 0.47 rad at 18 and 2 rad at one. $J_1(z)$ falls below $z/2$ as $z$ grows and has its maximum of 0.58 at $z = 1.84$.

:::{figure} generated/shot-noise/kick.svg
:name: fig-shot-kick
:alt: The Bessel function J1 of z and the straight line z over 2 against the kick z from 0 to 4 radians, with three shaded distributions of the kick for 150, 18 and 1 electrons a beamlet. The two curves part near z of 0.5.
The bunching a beamlet takes from a kick of size $z$, $J_1(z)$, against the first-order value $z/2$. The shaded curves are the distributions of the first-order kicks the load draws at 150, 18 and 1 electrons a beamlet, each scaled to the same height. At 150 electrons the kicks lie where the two curves agree. At one electron most kicks lie where $J_1$ is well below $z/2$.
:::

The mean over the draw of $x$, with $z = 2\sqrt{x/N}$, is Weber's second exponential integral ([](references.md#ref-watson)),

$$
  \langle |b_1|^2 \rangle = \int_0^\infty e^{-x}\, J_1^2\!\bigl(2\sqrt{x/N}\bigr)\, dx
     = e^{-2/N}\, I_1(2/N) ,
$$ (eq-shot-weber)

with $I_1$ the modified Bessel function. The ratio to the shot-noise value is

$$
  R(N) = N\, e^{-2/N}\, I_1(2/N) = 1 - \frac{2}{N} + \frac{5}{2N^2} - \ldots
$$ (eq-shot-r1)

The same holds for any single harmonic $h$ kicked alone, since the kick enters $b_h$ as $z_h = h\,a_h = 2\sqrt{x/N}$.

(shot-all)=
## The harmonics together

The kick of [](#eq-shot-kick) sums $H$ harmonics, and each kick also changes the bunching at the other harmonics. The mean is found without the expansion, from the pairs of macroparticles.

From [](#eq-shot-b) and [](#eq-shot-kick),

$$
  |b_h|^2 = \frac{1}{M^2} \sum_{j,l} e^{-ih(\theta_j - \theta_l)}
     \prod_{k=1}^{H} \exp\!\Bigl( i h a_k \bigl[ \sin(k\theta_j + \phi_k) - \sin(k\theta_l + \phi_k) \bigr] \Bigr) .
$$ (eq-shot-pairs)

The difference of the two sines is $2 \sin(k\delta/2) \cos(\psi_k)$, with $\delta = \theta_j - \theta_l$ and $\psi_k$ a phase that is uniform because $\phi_k$ is. The mean of $e^{i\zeta\cos\psi}$ over a uniform $\psi$ is $J_0(\zeta)$. So the mean over $\phi_k$ of the factor $k$ is

$$
  J_0( q_k z_k ), \qquad q_k = \frac{2h}{k} \left| \sin\frac{k\delta}{2} \right| .
$$ (eq-shot-j0)

Write $g(q)$ for the mean of $J_0(q z)$ over the draw that sizes $z$. The draws of different harmonics are independent, so the product of the means is the mean of the product. The pair $(j, l)$ enters only through $\delta = 2\pi(j-l)/M$, so the double sum is $M$ times a single sum:

$$
  \langle |b_h|^2 \rangle = \frac{1}{M} \sum_{j=0}^{M-1} \cos\frac{2\pi j h}{M}\;
    \prod_{k=1}^{H} g\!\left( \frac{2h}{k} \left| \sin\frac{\pi j k}{M} \right| \right) .
$$ (eq-shot-general)

This holds for any kick whose size follows from its own draw of $x$. For the first-order kick $q z = 2\sqrt{\beta x}$ with $\beta = q^2/N$. The mean over the exponential draw follows from the series of $J_0$ term by term, and is Weber's first exponential integral:

$$
  g(q) = \int_0^\infty e^{-x} J_0\!\bigl(2\sqrt{\beta x}\bigr)\, dx
    = \sum_{m} \frac{(-\beta)^m}{(m!)^2} \int_0^\infty e^{-x} x^m\, dx = e^{-\beta} .
$$ (eq-shot-weber0)

With $q^2 = (2h^2/k^2)(1 - \cos k\delta)$ the product is one exponential:

$$
  \langle |b_h|^2 \rangle = \frac{1}{M} \sum_{j=0}^{M-1} \cos\frac{2\pi j h}{M}\;
    \exp\!\Bigl( -\frac{2h^2}{N} \sum_{k=1}^{H} \frac{1 - \cos(2\pi j k/M)}{k^2} \Bigr) .
$$ (eq-shot-closed)

This is exact for the first-order kick. It is a sum of $M$ terms, and `tests/scripts/shot_noise_ratio.py` evaluates it. The ratio the tables and figures below give is

$$
  R_h(N) = N\, \langle |b_h|^2 \rangle .
$$ (eq-shot-rh)

Three limits check the form. With one harmonic and dense phases the sum over $j$ is an integral over $\delta$, and the integral is $e^{-2/N} I_1(2/N)$ by the integral representation of $I_1$, which is [](#eq-shot-weber). To first order in $1/N$ only the term $k = h$ survives the sum over $j$, and it gives $1/N$ for every kicked harmonic. A harmonic above $H$ is not kicked and has no first-order term. Its bunching is of order $1/N^2$.

To second order the default beamlet, $M = 8$ with $H = 3$ harmonics kicked, has

$$
  R_1 \approx 1 - \frac{22}{9N}, \qquad R_2 \approx 1 - \frac{1.01}{N}, \qquad R_3 \approx 1 - \frac{2.00}{N} .
$$ (eq-shot-second)

The fundamental loses more than the $2/N$ of [](#eq-shot-r1), because the kicks at the second and third harmonic also spread its phases.

(shot-table)=
## The first-order kick against the electron count

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-ratio -->
:end-before: <!-- end tab-shot-ratio -->
```

:::{figure} generated/shot-noise/ratio.svg
:name: fig-shot-ratio
:alt: The ratio of imposed to physical noise power against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis. The curves rise from near zero below one electron through about 0.5 at four electrons to 1 above one hundred. Measured points with error bars lie on the curves.
The ratio $R_h(N)$ of [](#eq-shot-rh) under the first-order kick for the default beamlet, $M = 8$ with 3 harmonics kicked, at the fundamental and the two harmonics above it. The dashed curve is one harmonic kicked alone, [](#eq-shot-r1), and the dotted curve its first-order form. The points are the levels measured on loaded beams, with one standard deviation, from `check_shot_noise.py` ([](#shot-measured)). In the shaded band a beamlet holds less than one electron ([](#shot-clamp)).
:::

The ratio has no threshold. The deficit at the fundamental is close to $2.4/N$ above 20 electrons a beamlet. The worked slice of [](#shot-counts), at 48.8 electrons a beamlet, starts with 0.951 of the shot-noise power at the fundamental under the first-order kick.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-floor -->
:end-before: <!-- end tab-shot-floor -->
```

The table's last three columns are the kick of [](#shot-bessel). $N$ follows from the slice and the load by [](#eq-shot-n). At a given current it rises with the slice length, `window.n_wavelength`, and falls with the macroparticles per slice, `beam.n_particle`.

A slice whose beamlets hold different counts, as with nonuniform weights, has the charge-weighted mean of their ratios:

$$
  \langle |b_h|^2 \rangle\, N_\lambda = \frac{\sum_g N_g\, R_h(N_g)}{\sum_g N_g} ,
$$ (eq-shot-slice)

with $N_g$ the electrons of beamlet $g$. The beamlets are kicked independently, so their bunchings add in the mean square with the squares of their charges.

(shot-clamp)=
## Below one electron

A beamlet may stand for less than one electron. The tail slice of [](#shot-counts) is such a load when it is given the beamlets of every other slice: 6 electrons over 128 beamlets. The first-order kick at $N < 1$ would be larger than at one electron. Under that kick the loader sets $N$ to 1 when it is below 1.

The beamlet's bunching is then the bunching at one electron, $\langle|b_h|^2\rangle = R_h(1)$, whatever its charge. A slice of such beamlets has

$$
  R_h(N) = N\, R_h(1), \qquad N < 1 .
$$ (eq-shot-clamp)

$R_1(1)$ is 0.161 for the default beamlet. The tail slice above, at 0.049 of an electron a beamlet, starts with 0.008 of its shot-noise power. A real slice of 6 electrons has $\langle|b|^2\rangle = 1/6$. The loaded slice has $0.161/128 = 0.0013$.

The kick sized from $J_1$ sets nothing to 1. Its largest kick bounds it, and [](#eq-shot-general) holds at every $N$ ([](#shot-bessel)). The same tail slice then starts with 0.012 of its shot-noise power.

Three loads have such a beamlet. The first is a slice that holds less than one electron in all, which has one beamlet under any rule. The second is the keep mode, where a group is a particle of the bunch and its copies, and the bunch's own weights decide. The third is the sample mode with `beam.beamlet_min_electrons = 0`, which gives every slice the same beamlets. The load's report counts the slices of the first kind. It warns of the slices that hold one electron or more and have a group under one, with their share of the charge and the noise power they are loaded with by [](#eq-shot-slice):

```
<count> slices have a noise group of less than one real electron.
They hold <share> of the charge that carries noise, and the noise power loaded
in them is <ratio> of their shot-noise power at the fundamental.
```

The warning is the same under both kicks, and the level it states is the level of the kick the load took.

The loader also takes a first-order kick $a_k$ above $2\pi$ modulo $2\pi$. At one electron that is $e^{-\pi^2} = 5\times10^{-5}$ of the draws at the fundamental, and fewer at every larger count. The closed form leaves it out.

:::{dropdown} Genesis4 source
`fel_fawley_noise` is transcribed from [`ShotNoise::applyShotNoise`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L34), with the electron count of the kick taken from the group's own charge. Genesis4 sets $N$ to 1 below one electron without a message ([line 47](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Loading/ShotNoise.cpp#L47)), which the warning's text refers to. The closed forms of this page are derived here and have no Genesis4 original. Neither has the kick sized from $J_1$, `fel_bessel_kick`.
:::

(shot-bessel)=
## The kick sized from $J_1$

The deficit of the first-order kick comes from one step. The kick is sized as if a kick of size $z$ gave the bunching $z/2$, and it gives $J_1(z)$ ([](#eq-shot-j1)). Lucifer sizes the kick from $J_1$. For each draw it takes the $z_k$ between 0 and $z_{\max}$ with

$$
  J_1(z_k) = \sqrt{\frac{x_k}{N}} .
$$ (eq-shot-bessel)

`fel_bessel_kick` finds the root by Newton's method, to rounding. The draws $u_k$ and $\phi_k$ and their order are those of the first-order kick. The two loads of one `ran_seed` are therefore the same shot, and this one has the larger kicks. The key is `beam.shot_noise_method`, `"bessel"` for this kick and `"first_order"` for [](#eq-shot-first). Fawley's algorithm ([](references.md#ref-fawley)) is the first-order kick, and sizing it from $J_1$ is this program's choice.

$J_1$ has its first maximum, $b_{\max} = 0.5819$, at $z_{\max} = 1.8412$. No kick gives a beamlet more bunching than $b_{\max}$ at the kicked harmonic. A draw with $x_k$ above $b_{\max}^2 N = 0.3386\,N$ has no root. The loader holds such a draw at $z_{\max}$ and counts it, and the load's report states the count with the share of the charge in the slices that have one. The share of the draws held is $e^{-0.3386 N}$: 0.2% at 18 electrons a beamlet, 3.4% at 10 and 71% at one.

With one harmonic kicked the beamlet's bunching is the bunching asked for, or $b_{\max}$ where the draw is held. The ratio is the mean of $\min(x,\, b_{\max}^2 N)$ over the exponential draw:

$$
  R(N) = 1 - e^{-b_{\max}^2 N} .
$$ (eq-shot-bessel-one)

It is 0.966 at 10 electrons a beamlet and 0.287 at one, where the first-order kick has 0.823 and 0.215. It holds below one electron as well, where it is close to $0.34\,N$.

The three kicks of the default beamlet act on one another as in [](#shot-all), and [](#eq-shot-general) gives the mean with

$$
  g(q) = \int_0^{z_{\max}} 2N\, J_1(z)\, J_1'(z)\, e^{-N J_1^2(z)}\, J_0(q z)\, dz
         \;+\; e^{-b_{\max}^2 N}\, J_0(q\, z_{\max}) .
$$ (eq-shot-g)

The first term is the draws that have a root, written as an integral over the kick's size through $x = N J_1^2(z)$. The second is the draws held. The integral has no closed form. `fel_shot_ratio` and `shot_noise_ratio.py` take it by Gauss-Legendre quadrature, 16 points in each of 4 panels, which a rule of three times the points reproduces to rounding.

To first order in $1/N$ the default beamlet has

$$
  R_1 \approx 1 - \frac{4}{9N}, \qquad R_2 \approx 1 + \frac{0.99}{N}, \qquad R_3 \approx 1 + \frac{0.00}{N} .
$$ (eq-shot-bessel-second)

Each is the first-order kick's coefficient of [](#eq-shot-second) less the $2/N$ that sizing from $J_1$ removes. What is left is the action of the other harmonics' kicks. The kicks at the second and third harmonic spread the fundamental's phases, which lowers it by $13/(18N)$, and the products of two kicks add $5/(18N)$ to it. The second harmonic gains about $1/N$ more than it loses, mostly from the square of the fundamental's kick, and so starts above the shot noise. The third harmonic is within 3% above 18 electrons a beamlet.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-bessel -->
:end-before: <!-- end tab-shot-bessel -->
```

:::{figure} generated/shot-noise/bessel.svg
:name: fig-shot-bessel
:alt: The ratio of imposed to physical noise power against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis, under the kick sized from J1. The fundamental's curve rises from near zero below one electron to 0.91 at ten electrons and to 1 above fifty. The first-order kick's fundamental, dashed, lies below it everywhere. Measured points with error bars lie on the curves.
The ratio $R_h(N)$ under the kick sized from $J_1$ for the default beamlet, at the fundamental and the two harmonics above it. The dashed curve is the first-order kick's fundamental, from [](#fig-shot-ratio). The dotted curve is one harmonic kicked alone, [](#eq-shot-bessel-one). The points are the levels measured on loaded beams, with one standard deviation, from `check_shot_noise.py` ([](#shot-measured)).
:::

The worked slice of [](#shot-counts), at 48.8 electrons a beamlet, starts with 0.991 of the shot-noise power at the fundamental, where the first-order kick gives it 0.951. A beamlet of one electron starts with 0.219, against 0.161.

The kick is sized one harmonic at a time. Sizing the three together, so that the beamlet had the three bunchings asked for at once, would remove what [](#eq-shot-bessel-second) leaves wherever such a kick exists. `tests/scripts/shot_noise_monte_carlo.py --joint` looks for it by Newton's method on a beamlet's own phases, continued from zero kick. It found the kick for every one of 2000 beamlets at 150 electrons a beamlet, for 90% at 50, for 42% at 18 and for 16% at 10 (2026-11-11). A beamlet of 8 phases kicked at three harmonics cannot take every triple of bunchings, even small ones. The kick of [](#eq-shot-bessel) has a root for every draw below $b_{\max}$, and what it leaves is known, so it is the kick the loader takes.

`tests/scripts/shot_noise_monte_carlo.py` kicks beamlets as the loader does, with no run of the program, and sums their bunching. At 150, 50, 18, 10, 4 and 1 electrons a beamlet, with $4\times10^8$ beamlets at each, the mean of $|b_h|^2 N$ at the three kicked harmonics agrees with [](#eq-shot-general) within two standard errors of the mean at all 18 points (2026-11-11). The standard error is $5\times10^{-5}$ or less and the largest difference is $9\times10^{-5}$. The share of the draws held agrees with $e^{-0.3386 N}$ to $10^{-5}$.

Three machines were tracked under each kick on 2026-11-11: `examples/sase`, `examples/flash1` and `examples/lcls`, each at `run.ran_seed` 1 to 6. The two runs of a seed are one shot, so their ratio is taken seed by seed, and the table gives its mean and standard error. The start-up power is the power summed over the [window](glossary.md#gl-window) one gain length into the line. The loaded power is the level the load's report states for the window, [](#eq-shot-slice) summed over the slices.

| machine | electrons a beamlet | loaded power over shot-noise power, first order | from $J_1$ | their ratio | start-up power, $J_1$ over first order | exit pulse energy, first order | from $J_1$ | $J_1$ over first order |
|---|---|---|---|---|---|---|---|---|
| sase | 73.2 | 0.9673 | 0.9938 | 1.0274 | 1.0266 ± 0.0012 | 3.40 ± 0.91 nJ | 3.48 ± 0.92 nJ | 1.026 ± 0.003 |
| FLASH1 | 16700 at the peak | 0.9995 | 0.9999 | 1.0004 | 1.000368 ± 0.000004 | 48.40 ± 0.90 µJ | 48.40 ± 0.90 µJ | 1.00002 ± 0.00001 |
| LCLS | 1250 at the peak | 0.9981 | 0.9996 | 1.0015 | 1.001585 ± 0.000007 | 278.9 ± 1.2 µJ | 279.0 ± 1.2 µJ | 1.0001 ± 0.0001 |

The start-up power follows the closed forms on all three machines. FLASH1 and LCLS saturate, and the kick moves their exit pulse energy by 2e-5 and 1e-4 of itself, where the standard error over the six seeds is 1.9% and 0.4%. The sase example does not saturate, and its exit pulse energy is higher by 2.6%, as its noise is. A run with many electrons a beamlet is therefore changed little by the kick. The kick matters where a beamlet holds tens of electrons or fewer: short wavelengths, low currents and large `beam.n_particle`.

(shot-weight)=
## The noise through the weights

The two kicks move the phases and leave every macroparticle of a beamlet at one weight. The shot noise can be imposed the other way round. The phases stay the quiet beamlet's, and the weights carry the noise. This is the load of McNeil, Poole and Robb ([](references.md#ref-mcneil)), put on a beamlet as Hwang and Qiang put it ([](references.md#ref-hwang-qiang-weight)). It is `beam.shot_noise_method = "weight"`.

A macroparticle of a beamlet of $N$ electrons stands for $N/M$ of them. Let it stand for a whole number $c_j$ of electrons instead, drawn as the electrons of a real beam fall: independently, each into any of the $M$ macroparticles with the same probability. Independent counts of that kind are Poisson numbers of mean $N/M$ and variance $N/M$. The beamlet's source at harmonic $h$ is the sum over its electrons,

$$
  S_h = \sum_{j=0}^{M-1} c_j\, e^{-ih\theta_j} .
$$ (eq-shot-s)

The mean of $S_h$ is $(N/M)\sum_j e^{-ih\theta_j}$, which is zero for the evenly spaced phases at every $h$ that is not a multiple of $M$. Its mean square is the sum of the variances of the counts:

$$
  \langle |S_h|^2 \rangle = \sum_j \langle (c_j - N/M)^2 \rangle = N,
  \qquad \langle |b_h|^2 \rangle = \frac{\langle |S_h|^2 \rangle}{N^2} = \frac{1}{N} .
$$ (eq-shot-weight)

This is [](#eq-shot-real), and no step of it is an approximation. The bunching is linear in the counts, so there is no order of a kick to stop at and no largest bunching. It holds at any $N$, whole or not, and at every harmonic below $M$: 1 to 7 for the default beamlet, where the kicks reach 3. At one electron a beamlet most macroparticles draw no electron and a few draw one. The load is then the electrons themselves, each at one of $M$ phases.

### The bunch's charge, held

Independent counts leave the bunch's charge to chance. The loader holds it. A bunch of charge $Q$ is counted as $n$ electrons, the whole number nearest $Q/e$ and at least one. The counts of all the macroparticles of the bunch are then one multinomial draw of $n$ electrons, each macroparticle with a probability $p_j$ in proportion to the weight the quiet load gave it. The counts sum to $n$ in every draw. They are no longer independent: the variance of $c_j$ is $n p_j (1 - p_j)$ and the covariance of two counts is $-n p_j p_k$. For any set of macroparticles, a slice or a beamlet,

$$
  \langle |S_h|^2 \rangle = n \sum_j p_j \;-\; n\, \Bigl| \sum_j p_j\, e^{-ih\theta_j} \Bigr|^2 .
$$ (eq-shot-multinomial)

The second term is the square of the mean source. The probabilities are equal within a beamlet and the beamlet's phases are evenly spaced, so it is zero beamlet by beamlet. What is left is $n \sum_j p_j$, the mean count of the set. [](#eq-shot-weight) therefore holds with the bunch's charge held, for a slice of any charge and for any weights across the beamlets.

The charge of one counted electron is $Q/n$, which is $e$ only where $Q/e$ is whole. The loaded charge is $Q$ exactly. The noise power goes as the square of that charge times $n$, so

$$
  R = \frac{Q/e}{n}
$$ (eq-shot-weight-ratio)

at every harmonic and in every slice. $R$ is within $1/(2n)$ of 1. That is $10^{-9}$ for the $5\times10^8$ electrons of FLASH1's window, and it is the whole cost of a count that is not whole. A bunch of less than half an electron is counted as one.

A slice's count is a binomial number: $n$ trials, each with the slice's share of the charge as its probability. Its mean is the slice's $N_\lambda$ and its variance is $N_\lambda (1 - N_\lambda/n)$. The charge of a slice thus varies from shot to shot by $1/\sqrt{N_\lambda}$ of itself, as the real beam's does: 0.7% for the 18740 electrons of a slice of `examples/sase`. The current profile of the load is the bunch's in the mean.

`fel_weight_noise` takes the draw in order. A slice takes a binomial share of the electrons not yet given, by its share of the charge not yet loaded. Each macroparticle then takes a binomial share of what its slice has left. `fel_ran_binomial` is exact in distribution: it is Hörmann's transformed rejection ([](references.md#ref-hormann)) where ten or more successes are expected, and a sum of the probabilities upward from zero successes below.

A beamlet none of whose copies drew an electron is not loaded. A copy that drew none, in a beamlet that holds charge, stays at zero weight, so a loaded slice is still whole beamlets. A slice that drew no electron has no particles. `beam.beamlet_min_electrons` is a rule for the kicks and this load does not read it: a slice of few electrons loads few beamlets with no rule.

### The keep mode

In the keep mode a beamlet is a particle of the bunch and its copies, and every particle keeps its own charge. A particle of $N$ electrons is counted as $k = \lfloor N \rfloor$ electrons or as $k + 1$, and the charge of one counted electron is the particle's charge over its count. Its noise power goes as its charge squared over the count. The count is $k$ with the probability

$$
  P(k) = \frac{k\,(k + 1 - N)}{N} ,
$$ (eq-shot-kept)

which makes the mean of one over the count $1/N$. The ratio is then 1 for every particle of one electron or more. The count falls on the $M$ copies with equal probability, so at least one copy holds charge and the others may hold none. A particle under one electron is counted as one. One of its copies holds its whole charge, and the ratio is $N$.

```{include} generated/shot-noise.md
:start-after: <!-- begin tab-shot-three -->
:end-before: <!-- end tab-shot-three -->
```

:::{figure} generated/shot-noise/weight.svg
:name: fig-shot-weight
:alt: The ratio of imposed to physical noise power at the fundamental against the electrons a beamlet holds, from 0.1 to 1000 on a logarithmic axis. The first-order kick and the kick from J1 rise from near zero below one electron to 1 above one hundred. The weight load is a level line at 1, with measured points at five electron counts lying on it.
The ratio at the fundamental under the three values of `beam.shot_noise_method`. The points are the weight load's levels at harmonics 1 to 7 measured on loaded beams by `check_shot_noise.py`, with one standard deviation, set side by side at each electron count.
:::

### What else differs

The weight load carries shot noise at harmonics 4 to 7 of the default beamlet, which the kicks leave quiet. On $M$ phases the harmonics $h$ and $M - h$ of one beamlet have the same magnitude, since $e^{-i(M-h)\theta_j}$ is the conjugate of $e^{-ih\theta_j}$ up to the beamlet's base phase. The beamlets of a slice stand at different base phases, so the two harmonics of the slice are independent sums.

The weights of a beamlet's copies differ, and at few electrons some are zero. Every pass over the particles takes each particle's own weight. The radiated power, the [FP32](glossary.md#gl-fp32) [twin](glossary.md#gl-twin) and the [device](glossary.md#gl-device) are held on such loads by the [validation suite](glossary.md#gl-validation-suite).

A weight load and a kick load of one `ran_seed` are different shots, since the weight load draws its counts and the kicks draw phases and sizes.

`tests/scripts/shot_noise_monte_carlo.py --weight` draws the counts with no run of the program, in three forms: independent Poisson counts on one beamlet, a slice of 16 beamlets with its total held, and the keep mode's particle. At 150, 18, 4, 3.76, 1 and 0.25 electrons a beamlet, with $10^8$ draws of the first and third and $1.25\times10^7$ of the second, the means at harmonics 1 to 7 agree with the forms within 2.7 standard errors of the mean at all 126 points (2026-11-11). The standard error is $1.4\times10^{-4}$ to $2.8\times10^{-4}$ and the largest difference is $5\times10^{-4}$. At 3.76 electrons a beamlet the slice of 16 beamlets holds 60.16 electrons counted as 60, and its mean is 1.0024 to 1.0031 where [](#eq-shot-weight-ratio) gives 1.0027. A keep-mode particle under one electron reads its electron count exactly, with no spread.

Four loads were tracked under the weight load and under the kick from $J_1$ on 2026-11-11: the kept bunch of `examples/import`, `examples/sase`, `examples/flash1` and `examples/lcls`, each at `run.ran_seed` 1 to 6. The two loads of a seed are different shots, so each column is a mean over the six seeds with its standard error, and the ratio is the ratio of the means. The start-up power is the power summed over the window one gain length into the line.

| load | electrons a beamlet | loaded power over shot-noise power, $J_1$ | weights | start-up power, weights over $J_1$ | the forms' ratio | exit pulse energy, $J_1$ | weights | macroparticles, $J_1$ | weights | run time, $J_1$ | weights |
|---|---|---|---|---|---|---|---|---|---|---|---|
| import, kept bunch | 3.76 | 0.6201 | 1.0000 | 1.62 ± 0.16 | 1.61 | 9.14 ± 0.27 pJ | 14.82 ± 0.67 pJ | 400000 | 400000 | 54.0 s | 52.9 s |
| sase | 73.2 | 0.9938 | 1.0000 | 1.08 ± 0.08 | 1.006 | 3.48 ± 0.92 nJ | 2.41 ± 0.55 nJ | 196608 | 196608 | 39.7 s | 39.8 s |
| FLASH1 | 16700 at the peak | 0.9999 | 1.0000 | 1.034 ± 0.024 | 1.0001 | 48.4 ± 0.9 µJ | 48.9 ± 1.9 µJ | 880832 | 862161 | 18.9 s | 18.9 s |
| LCLS | 1250 at the peak | 0.9996 | 1.0000 | 1.009 ± 0.011 | 1.0004 | 278.9 ± 1.2 µJ | 277.4 ± 1.7 µJ | 1228800 | 1228800 | 62.7 s | 62.2 s |

The run times are medians, on 12 threads. The start-up power agrees with the forms on all four within its spread. The kept bunch, at 3.76 electrons a particle, starts 1.6 times higher under the weights, which is the 0.62 the kick from $J_1$ loads there. The three machines with many electrons a beamlet start and exit as they do under the kick, within the spread of six shots. The macroparticles tracked and the run time are the same. FLASH1 tracks 2% fewer macroparticles, since its 82 slices of the least charge draw no electron. The loaded charge is the input file's to rounding in every run.

The weight load is therefore the default. It loads the shot noise at every electron count and at every harmonic below the beamlet size, it holds the bunch's charge, and it takes the same run time. The kick from $J_1$ is within 1% of it above 50 electrons a beamlet and stays in the program, as does the first-order kick, which the checks that compare a load with Genesis4 set.

(shot-fewer)=
## Fewer beamlets in a slice of little charge

Under the two kicks the sample mode does not give every slice $n_p/M$ beamlets. A slice of $N_\lambda$ electrons is loaded with

$$
  n_b = \min\!\left(\frac{n_p}{M},\ \max\!\left(1,\ \left\lfloor \frac{N_\lambda}{N_{\min}} \right\rfloor\right)\right)
$$ (eq-shot-fewer)

beamlets, where $N_{\min}$ is `beam.beamlet_min_electrons`. The slice keeps its charge: the weight of a macroparticle is the slice charge over $n_b M$. A slice of $N_{\min}$ electrons or more then has no beamlet under $N_{\min}$ electrons. A slice of fewer has one beamlet. A slice with no charge has none and is loaded with no particles, as the keep mode loads a slice that no particle of the bunch falls in. The routine is `fel_beamlet_count`. The generator draws $n_b$ base samples for the slice. The [resampler](glossary.md#gl-resampler) reduces the slice's candidates to $n_b$ seeds by the random deletion that reduces them to $n_p/M$. The rule is about the noise, so it acts where `shot_noise` and the quiet start are both set. The slice with no charge is the exception: it has no particles in every load of the sample mode, with the noise or without it, unless the key is 0. The keep mode is outside it, since that mode keeps every particle of the bunch.

The default is $N_{\min} = 1$: a slice is loaded with no more beamlets than it has electrons. The tail slice of [](#shot-counts), 6 electrons, has 6 beamlets of one electron and starts with 0.219 of its shot-noise power, 0.161 under the first-order kick. Over 128 beamlets it starts with 0.012 and 0.008.

The measurements of this section were taken with the first-order kick, before the kick sized from $J_1$ was in the program, and the levels quoted with them are that kick's. A larger $N_{\min}$ raises the level in the slices it acts on, to 0.56 at 4 electrons, 0.87 at 18 and 0.95 at 50. It also acts on more slices. The worked core slice of [](#shot-counts) holds 6246 electrons. $N_{\min} = 18$ takes beamlets from it at 348 beamlets or more, 2784 macroparticles a slice, and $N_{\min} = 50$ at 125 beamlets or more, which the slice's own 128 already exceed. At $N_{\min} = 1$ a slice loses beamlets only when it was given more of them than it has electrons. Three machines were loaded with `run.load_only` at each value on 2026-11-10. The sase example, at 0.1 nm, holds 73.2 electrons a beamlet in every slice and LCLS 1250 at its peak, and neither has a slice under 50 electrons a beamlet. FLASH1 is a Gaussian bunch in a [window](glossary.md#gl-window) with slippage headroom, 351 slices of 512 beamlets:

| $N_{\min}$ | slices with fewer beamlets | their share of the charge | macroparticles | loaded power over shot-noise power, window |
|---|---|---|---|---|
| 0 | 0 | 0 | 1437696 | 0.9995 |
| 1 | 146 | 1.0e-5 | 880832 | 0.9995 |
| 4 | 162 | 4.7e-5 | 821904 | 0.9995 |
| 18 | 180 | 2.3e-4 | 752544 | 0.9995 |
| 50 | 192 | 6.2e-4 | 701472 | 0.9996 |

FLASH1 was then tracked at six seeds for each value. The exit pulse energy, as the mean and its standard error, is 48.2 ± 2.7 µJ at 0, 49.9 ± 1.9 µJ at 1, 48.8 ± 2.0 µJ at 4, 48.2 ± 1.5 µJ at 18 and 45.4 ± 1.5 µJ at 50. The values differ by less than the shots do, since the slices the rule acts on hold 1e-5 to 6e-4 of the charge. The run at 1 tracks 0.61 of the macroparticles of the run at 0 and took 19 s where that run took 23 s. The sase example is the same bits at every value up to 50. A Gaussian bunch of 3 kA at its peak was tracked through one undulator at twelve seeds, the load of [](validation.md#val-shot-noise-under-weights) with 52 of its 85 slices reduced at the default. Its start-up power inside the mode's angle is 0.99 ± 0.06 of the flat window's at 0 and 0.98 ± 0.04 at 1. The measurements therefore do not choose among the values. The default is the one value that takes a beamlet from a slice only when the slice has no electron for it.

A slice with fewer beamlets has fewer samples of the transverse distribution and of the energy spread: $n_b$ draws where the other slices have $n_p/M$. A slice of one beamlet is a point in the transverse plane with one energy. At the default the slice has one sample for each electron it holds, which is what the real slice has. The header's `Load` line states the beamlets of a slice that was not reduced and the fewest any slice with particles has, and a line under it the slices loaded with none. The [device](glossary.md#gl-device) and the [FP32](glossary.md#gl-fp32) [twin](glossary.md#gl-twin) take each slice's own count, and the twin's residual guard is read on the slices of 256 particles or more ([](performance.md)). The random numbers of the slices after the first reduced slice are not those of a load with the key at 0, so the two are different shots at one `ran_seed`. A load in which no slice is reduced is the same bits either way.

With `beam.beamlet_min_electrons = 0` every slice has $n_p/M$ beamlets, as in Genesis4. The [tiers](glossary.md#gl-tier) of the [validation suite](glossary.md#gl-validation-suite) that compare a load with Genesis4's set it.

(shot-measured)=
## The load, measured against the closed form

`check_shot_noise.py` loads a uniform bunch of 256 slices with 15 seeds and reads the bunching of every slice at harmonics 1 to 3 from the particle dump. Under the first-order kick it runs at 146, 50, 18 and 4 electrons a beamlet and at a quarter of an electron, each with uniform weights and with beamlet weights alternating between 0.25 and 1.75 of the mean. The mean of $|b_h|^2 N_\lambda$ over the samples is held to [](#eq-shot-slice) within five standard deviations of the mean. These loads set `beam.beamlet_min_electrons = 0`, so that a beamlet holds the count stated.

Under the kick sized from $J_1$ it runs at 146, 50, 18, 10 and 1 electrons a beamlet with uniform weights, held to that kick's form the same way. Five standard deviations of such a mean are 8% at one harmonic, which is more than the two kicks differ by above 18 electrons. The two kicks take the same draws, so a seed's two loads are one shot kicked twice, and the difference of their statistics slice by slice has a far smaller spread than either. The mean difference is held to the difference of the two closed forms within five standard deviations of its own mean. At 146 electrons a beamlet it reads 0.0140 against 0.0139 with a bound of 0.0008. The level the loader prints for each of these loads, from `fel_shot_ratio`, is held to `shot_noise_ratio.py`'s to the four decimals printed.

The weight load runs at 146, 18, 4 and 1 electrons a beamlet and at a quarter of an electron, with uniform weights. Its statistic is held to [](#eq-shot-weight-ratio) at each of harmonics 1 to 7. The loaded charge is held to the input file's to $10^{-12}$, every weight to a whole number of the counted electron's charge, and the variance of the slices' counts to the binomial variance within five standard deviations. A Gaussian bunch whose last slices hold a hundredth of an electron is loaded the same way: 30% of its slices draw no electron and have no particles, the bunch is a whole number of electrons, and the statistic of the other slices is held to 1. In the keep mode every particle of the bunch must be present at the charge the kick load gives it.

A slice of several wavelengths is loaded over all its periods ([](fel-physics.md#sec-quiet)). A Gaussian bunch in slices of 3 and of 12 wavelengths is loaded with no noise, under the weight load and under the first-order kick. With no noise the sum of $w\,e^{-i(m/n)\theta}$ over a slice, at the $m$ that are not multiples of the slice's $n$ wavelengths, is held to $10^{-12}$ of the slice's charge. Under the weight load its square is held to the shot noise within five standard deviations, and the electrons in each period of a slice to the binomial spread. The first-order kick's is recorded: 0.36 of the shot noise at 3 wavelengths and 0.46 at 12.

A Gaussian bunch whose last slices hold two electrons is then loaded under the kicks at the default $N_{\min}$ and at $N_{\min} = 18$, and at the default under both kicks. Each slice is held to the charge it has with the key at 0 and to the beamlets of [](#eq-shot-fewer), and the slices with fewer beamlets to the closed form at the electrons their beamlets hold. A flat bunch over the middle half of the window is loaded both ways: a slice with no charge must have no particles at the default and the input file's count at 0. The keep mode and the two loads whose beamlets were made outside the loader run under both kicks as well, with the same paired difference. The check runs in every pass of the [validation suite](glossary.md#gl-validation-suite), and its levels are in [](validation.md#val-levels-shot-noise). They are the points of [](#fig-shot-ratio), [](#fig-shot-bessel) and [](#fig-shot-weight).

(shot-scope)=
## Scope

The weight load needs the beamlets this loader makes, in either load mode. A load without the quiet start is refused under `"weight"`, with `"bessel"` in the message. The closed forms [](#eq-shot-general) and [](#eq-shot-closed) are for a beamlet of evenly spaced phases kicked by `fel_fawley_noise`. That covers the beamlets Lucifer's loader makes in both load modes. It also covers copies a bunch file brings from another code's loader, which `fel_group_noise` kicks by the same rule when the copies are evenly spaced.

Where the groups are the occupants of a deposit cell, one group holds the copies of several particles at different base phases. The terms of [](#eq-shot-ja) that do not depend on the base phase are the same, so the dense form, [](#eq-shot-general) at large $M$, applies with $N$ the electrons of the whole group. Under the kick sized from $J_1$ a group of more than 64 phases is taken as 64, which differs at the 63rd order of the kick. A group whose phases are not copies at even spacing is outside this page.

The page is about the load. What the FEL does with the noise, and how the start-up power depends on the transverse grid, is in [](startup-noise.md).
