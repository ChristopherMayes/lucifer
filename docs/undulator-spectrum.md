---
title: "Undulator radiation: the spectrum of one undulator and its formulas"
short_title: Undulator spectrum
---

An electron that crosses an undulator of $N$ periods radiates a line at each harmonic of a fundamental frequency, and the shape of the line, its level and its dependence on the angle of observation are known in closed form. A beam whose current is far below that required for gain radiates the incoherent sum of these lines. Lucifer computes the same radiation from the source and field equations of [](fel-physics.md), which were written for the amplified field and not for the spontaneous one. This page compares the two for the fundamental of two undulators: that of FLASH1, which is planar, and the Delta afterburner of the LCLS, which is helical. The reference in every comparison is a formula. No other code is involved.

The comparison is made with the default [propagator](glossary.md#gl-propagator) and the averaged model unless stated otherwise, and gives four results.

- On the axis the spectrum is the line $\mathrm{sinc}^2(\pi N\nu)$ of undulator radiation, in width, in the positions of its zeros and in shape ([](#fig-und-line), [](#tab-und-line)).
- The level on the axis, in joules per unit angular frequency and per steradian for each electron, is equal to the closed form for the planar and for the helical undulator ([](#und-level)).
- Away from the axis the line is displaced to lower frequency by an amount that is equal to the red shift of undulator radiation to first order in the square of the angle. With `field.propagator = "spectral"` the displacement is equal to the red shift itself ([](#und-angle)).
- The level away from the axis is that of the source alone. The dependence of undulator radiation on the angle of observation is absent from the averaged model and from the unaveraged model ([](#und-factor)).

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-line -->
:end-before: <!-- end fig-und-line -->
```

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin tab-und-line -->
:end-before: <!-- end tab-und-line -->
```

(und-formulas)=
## The radiation of one electron

Let an electron of Lorentz factor $\gamma$ cross $N$ periods of an undulator of period $\lambda_u$, and let $K = eB\lambda_u/(2\pi mc)$ be the deflection parameter, where $B$ is the peak field of a planar undulator and the constant magnitude of the field of a helical one. The mean square of the normalized vector potential is $a_w^2 = K^2/2$ for the planar undulator and $a_w^2 = K^2$ for the helical one. The mean longitudinal velocity of the electron corresponds to the Lorentz factor $\gamma_z$, where

$$
  \gamma_z^2 = \frac{\gamma^2}{1 + a_w^2} .
$$ (eq-und-gammaz)

Far from the undulator the energy radiated per unit angular frequency and per unit solid angle in the direction $\mathbf n$ is ([](references.md#ref-jackson), Eq. 14.67)

$$
  \frac{d^2W}{d\omega\,d\Omega} = \frac{e^2\omega^2}{16\pi^3\epsilon_0 c}
  \left| \int \mathbf n\times(\mathbf n\times\boldsymbol\beta)\; e^{\,i\omega(t - \mathbf n\cdot\mathbf r/c)}\,dt \right|^2 ,
$$ (eq-und-jackson)

where $\mathbf r(t)$ is the position of the electron and $c\boldsymbol\beta$ its velocity. The motion repeats in every period. The integral over $N$ periods is therefore the integral over one period multiplied by a sum of $N$ phase factors, and near the fundamental

$$
  \frac{d^2W}{d\omega\,d\Omega} = W_1(\theta)\, N^2\, \mathrm{sinc}^2(\pi N\nu), \qquad
  \nu = \frac{\omega}{\omega_1(\theta)} - 1 ,
$$ (eq-und-line)

where $\mathrm{sinc}\,u = \sin u/u$. The line has zeros at $\nu = m/N$ for every nonzero integer $m$ and a full width at half maximum of $0.8859/N$. Its center is the fundamental frequency in the direction of observation,

$$
  \omega_1(\theta) = \frac{\omega_1(0)}{1 + \gamma_z^2\theta^2}, \qquad
  \omega_1(0) = \frac{4\pi c\,\gamma_z^2}{\lambda_u} ,
$$ (eq-und-resonance)

where $\theta$ is the angle between $\mathbf n$ and the axis. The radiation at a larger angle has a lower frequency. The quantity $W_1(\theta)$ is Eq. [](#eq-und-jackson) evaluated over one period at $\omega_1(\theta)$. On the axis it has the closed forms ([](references.md#ref-kim))

$$
  W_1(0) = \frac{e^2\gamma^2}{4\pi\epsilon_0 c}\,\frac{K^2\,[JJ]^2}{(1 + K^2/2)^2}, \qquad
  [JJ] = J_0(\xi) - J_1(\xi), \quad \xi = \frac{K^2}{4 + 2K^2}
$$ (eq-und-planar)

for the planar undulator and

$$
  W_1(0) = \frac{e^2\gamma^2}{4\pi\epsilon_0 c}\,\frac{2K^2}{(1 + K^2)^2}
$$ (eq-und-helical)

for the helical one, in which the two polarizations are summed. The factor $[JJ]$ is the coupling factor of the planar undulator, which appears in the averaged model as Eq. [](fel-physics.md#eq-fc). Away from the axis $W_1(\theta)$ is evaluated numerically from Eq. [](#eq-und-jackson) on the ideal trajectory (`tests/scripts/undulator_radiation.py`). The numerical value on the axis is equal to Eqs. [](#eq-und-planar) and [](#eq-und-helical) to six figures.

A transverse grid resolves the transverse wavenumber $k_\perp$ of the field and not the angle. A component of wavenumber $k$ that propagates at the angle $\theta$ has $k_\perp = k\theta$ in the paraxial approximation. With $k_0 = \omega_1(0)/c$ and

$$
  x = \gamma_z^2 \left(\frac{k_\perp}{k_0}\right)^2 ,
$$ (eq-und-x)

Eq. [](#eq-und-resonance) is equivalent to

$$
  \nu\,(1 + \nu) = -x, \qquad \nu = \frac{\omega_1}{\omega_1(0)} - 1 ,
$$ (eq-und-nux)

so that $\nu = -x - x^2$ to second order in $x$. This form of the resonance is used in [](#und-angle). The fundamental has no component with $x > 1/4$.

(und-undulators)=
## The two undulators

:::{table} The undulators of the comparison. The FLASH1 undulator is that of `examples/flash1`, of which one segment is used and no quadrupole. The Delta undulator is taken in its circular mode, with the beam energy and the photon energy of [](references.md#ref-delta), and its field is the one for which the fundamental is at that photon energy.
:name: tab-und-undulators

| | FLASH1 | Delta |
|---|---|---|
| type | planar | helical |
| period $\lambda_u$ | 27.3 mm | 32 mm |
| periods $N$ | 165, and 55 | 100 |
| field $B$ | 0.47 T | 0.800 T |
| $K$ | 1.198 | 2.391 |
| beam energy | 668.5 MeV | 4 GeV |
| fundamental | 13.70 nm | 1.7537 nm (707 eV) |
| lattice | `tests/bmad/und_flash1.bmad` | `tests/bmad/und_delta.bmad` |
:::

(und-emitter)=
## One slice as the emitter

The radiation of a beam is the sum of the radiation of its electrons, and in the averaged model the sum is taken over each [slice](glossary.md#gl-slice) in turn. Let particle $j$ of a slice represent $w_j$ electrons and have the phase $\theta_j$, and let

$$
  S = \sum_j w_j\, e^{-i\theta_j}
$$ (eq-und-s)

be the sum over the particles of the slice, as in [](shot-noise.md#eq-shot-s). In each step the slice adds to the [plane](glossary.md#gl-plane-field) of the field that it overlaps a source proportional to the complex conjugate of $S$ (Eq. [](fel-physics.md#eq-source)). With a step of one period and slices of one wavelength, the field then advances by one plane with respect to the beam ([slippage](glossary.md#gl-slippage), [](fel-physics.md#sec-slippage)). A slice therefore radiates into a new plane in every period, and after $N$ periods it has radiated $N$ planes, which lie ahead of it in the [pulse](glossary.md#gl-pulse). Each plane has received the source of one period and has then propagated in free space to the exit.

Free propagation changes the transverse profile of a plane and leaves its integral over $x$ and $y$ unchanged, because the integral is the component of zero transverse wavenumber, on which the kernel of Eq. [](fel-physics.md#eq-k2) is 1. Let $E_0(p)$ be this integral for the plane at position $p$ of the pulse. If the current is far below that required for gain, the field does not act back on the particles, and

$$
  E_0(p_j) = \kappa\, S_j^* ,
$$ (eq-und-plane)

where $S_j$ is the sum of Eq. [](#eq-und-s) at period $j$, $p_j$ is the plane into which the slice radiated at that period, and $\kappa$ is one constant. The spectrum on the axis is consequently the squared magnitude of the Fourier sum of $S_j^*$ over the periods. For a slice in which $S$ does not change along the undulator, the pulse on the axis consists of $N$ equal planes, and its spectrum is the line of Eq. [](#eq-und-line) multiplied by $|S|^2$.

One loaded slice thus gives the line in a single run, with no average over seeds, although its $S$ is one draw of the shot noise. The runs of this page load one slice of 62 electrons in a longer [window](glossary.md#gl-window), with a transverse size below 10 nm and a relative energy spread of $10^{-6}$, for which the change of $S$ along the undulator is negligible. In the run of the FLASH1 undulator the planes of the pulse differ from $\kappa S_j^*$ by at most $1.2 \times 10^{-9}$ of the largest. The effect of a beam of finite emittance and energy spread is the subject of [](#und-beam).

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-pulse -->
:end-before: <!-- end fig-und-pulse -->
```

(und-on-axis)=
## The on-axis spectrum

The spectrum is computed from the field of the pulse file, which holds the complex amplitude $E(x, y, p)$ that the run propagated to the exit of the undulator. The field is summed over $x$ and $y$ in each plane, which gives $E_0(p)$, and the Fourier sum of $E_0(p)$ over the planes of the pulse is evaluated at each relative frequency $\nu$. The squared magnitude of this sum, divided by its value at the peak, is the quantity that [](#fig-und-line) shows. It is the transform of Eq. [](#eq-und-transform) at zero transverse wavenumber. The sums $S_j$ of the [statistics file](glossary.md#gl-statistics-file) do not enter the figure. They are used in the comparison of Eq. [](#eq-und-plane) and in the normalization of [](#und-level) only.

The emitter of these runs is a slice whose $S$ does not change along the undulator, so that the pulse on the axis consists of $N$ equal planes ([](#fig-und-pulse)). The shape of the spectrum then follows from three properties of the run: the field advances by exactly one plane in each period, every period deposits the same source, and no part of the pulse is lost at the ends of the undulator or of the window. The agreement of [](#fig-und-line) with the line establishes these three properties. It is not an independent test of the shape of the line, because $N$ equal planes can have no other spectrum. The spectra of pulses whose planes are not equal are those of [](#und-beam) and [](#und-unaveraged).

The pulse of a run is sampled once per slice, which is once per wavelength here. The spectrum of $N$ equal samples is

$$
  \left( \frac{\sin \pi N\nu}{N \sin \pi\nu} \right)^2 = \mathrm{sinc}^2(\pi N\nu)\left( 1 + \frac{(\pi\nu)^2}{3} + \ldots \right) ,
$$ (eq-und-sampled)

which has the zeros of the line and differs from it by $(\pi\nu)^2/3$ of its value. The difference is largest in the side lobes of the shortest undulator and is the whole of the difference between the runs and the line in [](#fig-und-line): it is $1.1 \times 10^{-4}$ of the peak at 55 periods and $1.2 \times 10^{-5}$ at 165.

The width of the line is equal to $0.8859/N$ to four figures for the three undulators ([](#tab-und-line)), and the ratio of the width at 55 periods to that at 165 is 3.0004. The level of the line at its zeros is less than $10^{-7}$ of the peak. The center of the line is at the fundamental of the undulator as Eq. [](#eq-und-resonance) gives it from the field and the energy that the lattice states, to $10^{-4}$ of the line width.

(und-level)=
## The on-axis spectral intensity

The field of a run is the complex amplitude $E(x, y, p)$ of a wave of wavenumber $k_0$, in volts per meter, on a grid of cells $\Delta x\,\Delta y$ and planes that are one slice apart, which is the time $\Delta t = \lambda/c$ for slices of one wavelength. The intensity of the wave is $\epsilon_0 c\,|E|^2/2$, and the energy of the pulse is the sum of the intensity over the cells and the planes multiplied by $\Delta x\,\Delta y\,\Delta t$. Let

$$
  \tilde E(\mathbf k_\perp, \nu) = \sum_{x, y, p} E(x, y, p)\; e^{-i\mathbf k_\perp\cdot\mathbf x}\, e^{-2\pi i\nu p}\;\Delta x\,\Delta y\,\Delta t
$$ (eq-und-transform)

be the transform of the field in the transverse plane and along the pulse. By Parseval's theorem the energy of the pulse is the integral of $\epsilon_0 c\,|\tilde E|^2/2$ over $d^2k_\perp\,d\omega/(2\pi)^3$. A component of transverse wavenumber $\mathbf k_\perp$ propagates in the direction $\mathbf k_\perp/k$, so that $d^2k_\perp = k^2\,d\Omega$ and

$$
  \frac{d^2W}{d\omega\,d\Omega} = \frac{\epsilon_0 c}{2}\,\frac{k^2}{(2\pi)^3}\,\bigl|\tilde E(\mathbf k_\perp, \nu)\bigr|^2 .
$$ (eq-und-level)

On the axis $\tilde E(0, \nu)$ is the Fourier sum of $E_0(p)$ of [](#und-emitter) multiplied by $\Delta t$, and by Eq. [](#eq-und-plane) the level at the peak of the line is to be compared with $W_1(0)$ multiplied by the squared magnitude of the sum of $S_j$ over the periods. The values of $S_j$ are taken from the [statistics file](glossary.md#gl-statistics-file), which records the bunching and the charge of every slice at every step.

The ratio of the two is 1.0000 for the planar undulator at 165 and at 55 periods and for the helical undulator ([](#tab-und-line)). The source of the averaged model, Eq. [](fel-physics.md#eq-source), thus contains the coupling factor, the wavelength and the energy in the combination of Eq. [](#eq-und-planar), although it was transcribed from a code for the amplified field. For the helical undulator the single complex field of the run represents the circularly polarized wave, and its level is equal to Eq. [](#eq-und-helical), in which the two linear polarizations are summed.

(und-angle)=
## The red shift with angle

The component of the field at the transverse wavenumber $k_\perp$ is multiplied in each step of length $\lambda_u$ by the kernel of Eq. [](fel-physics.md#eq-k2), which is the phase factor $\exp(-i k_\perp^2\lambda_u/2k_0)$. A plane that was radiated $m$ periods before the exit has therefore acquired $m$ times this phase, and the planes of the pulse, which are one wavelength apart, differ in phase from one to the next by

$$
  \frac{k_\perp^2\,\lambda_u}{2k_0} = 2\pi\,\gamma_z^2\left(\frac{k_\perp}{k_0}\right)^2 = 2\pi x ,
$$

where $\lambda = \lambda_u/2\gamma_z^2$ has been used. A phase that advances by $2\pi x$ per wavelength along the pulse is a displacement of the line by

$$
  \nu = -x .
$$ (eq-und-envelope)

The envelope equation of the default propagator therefore displaces the line at each transverse wavenumber to lower frequency, by an amount that is equal to the red shift of Eq. [](#eq-und-nux) to first order in $x$. The two differ by $x^2$ because the kernel takes the wavenumber of the carrier, $k_0$, in its denominator for every frequency of the pulse. This is the approximation of the slowly varying envelope, and the kernel is that of Genesis4.

[](#fig-und-ridge) shows the frequency of the peak at each transverse wavenumber. The runs lie on Eq. [](#eq-und-envelope) to $3 \times 10^{-4}$ of the line width over the whole range of the grid, in both transverse directions and for both undulators. At $x = 0.14$ the resonance is 5.2 line widths below the run for 165 periods of the FLASH1 undulator.

The [spectral propagator](spectral-field.md), `field.propagator = "spectral"`, transforms the field along the pulse as well and applies to each frequency the kernel of its own wavenumber, $k_0(1 + \nu)$. The phase per plane is then $2\pi x/(1 + \nu)$, the displacement satisfies $\nu(1 + \nu) = -x$, and this is Eq. [](#eq-und-nux). The run with the spectral propagator, of 55 periods of the FLASH1 undulator, lies on the resonance to 0.03 of the line width for $x \le 0.08$, where the line $-x$ is 0.46 of the line width away.

The difference between the two propagators is small at the angles of an amplified mode. At the edge of the [source filter](glossary.md#gl-source-filter) for the beam of FLASH1, $x = 0.007$ ([](#und-filter)), the term $x^2$ is $5 \times 10^{-5}$, which is 0.009 of the line width at 165 periods.

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-ridge -->
:end-before: <!-- end fig-und-ridge -->
```

(und-factor)=
## The off-axis spectral intensity

The source of Eq. [](fel-physics.md#eq-source) is deposited at the transverse position of each particle and does not depend on the direction in which the radiation leaves. The level of the line at a transverse wavenumber is therefore the level on the axis multiplied by the squared magnitude of the transverse transform of the deposited source. For the beam of these runs, whose size is far below a cell of the grid, that transform is the response of the bilinear deposit to a point, which decreases toward the Nyquist wavenumber of the grid.

To separate the model from the grid, the level of each run along the ridge is divided by the spectrum of what the slice deposits in one step, which is taken from the planes of the same pulse. The quotient is 1 at every wavenumber, to $3 \times 10^{-7}$ for the FLASH1 undulator and to $10^{-7}$ for the Delta undulator ([](#fig-und-angular)). The averaged model applies no factor that depends on the angle.

Undulator radiation has such a factor. Expressed per unit area of transverse wavenumber, as the runs are, it is $W_1(\theta)/k^2$ at $\theta = k_\perp/k$. For the planar undulator at $x = 0.14$ it is 0.34 of its value on the axis in the plane of the oscillation and 1.05 in the perpendicular plane. For the helical undulator at $x = 0.2$ it is 0.38. A run consequently contains more radiation at wide angles in the plane of the oscillation than an undulator emits, by these factors. This is in accord with the observation of Fawley ([](references.md#ref-fawley-2003)) that a code of this kind reproduces the intensity on the axis and the power within the central cone of the spontaneous emission, and that its total spontaneous power depends on the grid. The radiation at wide angles and its effect on a run that starts from noise are the subject of [](startup-noise.md), and the source filter is the remedy that Lucifer applies.

The unaveraged model was run on 55 periods of the FLASH1 undulator. Its level along the ridge is within 0.01 of the spectrum of its own source, so that it applies no angular factor either. Its ridge lies on Eq. [](#eq-und-envelope) to 0.02 of the line width, since it uses the same propagator.

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-angular -->
:end-before: <!-- end fig-und-angular -->
```

(und-filter)=
## The source filter

The [source filter](glossary.md#gl-source-filter), `run.source_filter`, removes the source above a transverse wavenumber that the run derives from the size of the beam ([](fel-physics.md#sec-source-filter)). It has no effect on a beam of no transverse size, for which the derived edge lies beyond the grid. The filter is therefore examined on one slice with the normalized emittance of FLASH1, 1.5 µm, for which the rms beam size is 107 µm and the edge is 83.5 µrad. The shot noise of such a beam radiates into every wavenumber of the grid, since each [beamlet](glossary.md#gl-beamlet) has a phase of its own.

[](#fig-und-filter) shows the energy at each wavenumber with the filter on, divided by that of the same beam with the filter off. The ratio is 1 below the edge and vanishes above it. The filter does not change the line on the axis. At the edge the angular factor of undulator radiation is 0.95 in the plane of the oscillation, so that the radiation which the filter retains is that for which the absence of the factor from the model is a correction of at most 5 percent.

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-filter -->
:end-before: <!-- end fig-und-filter -->
```

(und-beam)=
## Emittance and energy spread

An electron whose energy differs from the reference by the relative amount $\delta$ radiates a line that is displaced by $2\delta$, and an electron with a transverse angle or a betatron oscillation radiates a line that is displaced to lower frequency. In a run each beamlet has an energy and transverse coordinates of its own. Its phase with respect to the field then advances at a rate of its own, and the sum $S$ of the slice changes along the undulator. Equation [](#eq-und-plane) continues to hold: with the emittance of FLASH1 the planes differ from $\kappa S_j^*$ by at most $2 \times 10^{-4}$ of the largest, and with its relative energy spread, $6 \times 10^{-4}$, by at most $7 \times 10^{-4}$. These differences are of the order of the variation over the particles of the coupling of Eq. [](fel-physics.md#eq-source), which depends on the energy and on the transverse position of a particle, and which the bunching recorded in the statistics file does not contain. The level at the peak remains equal to $W_1(0)$ multiplied by the squared Fourier sum of $S_j$, to 0.2 percent.

The line of one slice is then the sum of the lines of its beamlets, each displaced by its own amount and multiplied by the complex amplitude that the shot noise assigned to the beamlet. The result is one draw of a random spectrum. Its peak is displaced by a fraction of the line width, in a direction and by an amount that differ from one seed to the next, and its zeros are filled ([](#fig-und-beam)). The spectrum of a bunch is the average over many such slices, which is the line convolved with the distribution of the displacements.

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-beam -->
:end-before: <!-- end fig-und-beam -->
```

(und-window)=
## The window and the slices

The radiation of a slice leaves the window at its head when the window is shorter than the [slippage](glossary.md#gl-slippage), and the planes that leave are kept in the [bank](glossary.md#gl-bank), from which the pulse file is assembled. In a window of 65 wavelengths, against a slippage of 165, the line of the FLASH1 undulator is equal to that of the window of 257 wavelengths in every digit that the check records.

With slices of three wavelengths the field advances by one plane in three periods, each plane of the pulse receives the source of three periods, and the pulse is sampled every three wavelengths. The width of the line is then 1.0001 of $0.8859/N$ and its level 1.0000 of the closed form.

(und-unaveraged)=
## The unaveraged model

The unaveraged model, `fel_method = unaveraged`, advances the particles through the period in [substeps](glossary.md#gl-substep) and deposits the source at each of them ([](fel-physics.md#sec-unaveraged)). The pulse of one slice is then not cut off between two planes: its first two and its last two planes are lower than the others ([](#fig-und-pulse)). For 55 periods the width of the line is 1.043 of $0.8859/N$, its level is 0.939 of the closed form, and it differs from $\mathrm{sinc}^2$ by 0.038 of the peak ([](#tab-und-line)). Its spectrum is the continuous colored curve of [](#fig-und-line). The minima of that curve are at 1.04 to 1.05 of the zeros of the line, as for a pulse of 52.9 equal planes, which is the sum of the magnitudes of the planes of [](#fig-und-pulse) in units of the largest. The curve is not symmetric about its peak: its first side lobe is 0.043 of the peak below the fundamental and 0.051 above it, where that of the line is 0.047. The cause of this asymmetry was not determined. These differences are those of a pulse whose ends are rounded over two planes in 55, and they decrease with $N$: a run at 165 periods, which is not part of the check, gave a level of 0.973.

(und-bunch)=
## A bunch

A bunch of 1024 loaded slices, each of about 62 electrons with shot noise, was run through 55 periods of the FLASH1 undulator with eight seeds. Every plane of the pulse then contains the radiation of up to 55 slices. Equation [](#eq-und-plane) holds for the sum: in each run the planes differ from $\kappa$ multiplied by the sum of $S_j^*$ over the slices and periods that radiated into them by at most $1.2 \times 10^{-6}$ of the largest. The spectrum of a run is therefore that of the loaded shot noise, and the comparison with the formula reduces to the statement that the mean of $|S|^2$ in a slice is equal to the number of its electrons, which is the subject of [](shot-noise.md).

The mean spectrum of the eight seeds is shown in [](#fig-und-bunch) with the line of one electron multiplied by the number of electrons. The energy radiated on the axis, summed over the line, is 0.97 of that of one electron multiplied by the number of electrons, with a standard deviation of 0.12 over the seeds.

```{include} generated/undulator-spectrum.md
:start-after: <!-- begin fig-und-bunch -->
:end-before: <!-- end fig-und-bunch -->
```

(und-measured)=
## How the comparison is checked

Every number of this page is measured by `tests/scripts/check_undulator_spectrum.py`, which is the `undulator-spectrum` section of the [validation suite](glossary.md#gl-validation-suite) and runs on both builds. Each value is held to a bound, and the values and bounds are listed in [](validation.md#val-levels-undulator-spectrum). The figures and their captions are drawn by `tests/scripts/report_undulator_spectrum.py` from the curves of the release build's run, at every [keystone](glossary.md#gl-keystone), and a caption whose numbers change fails the keystone.

(und-scope)=
## What is not reproduced

- The dependence of the level on the angle of observation, $W_1(\theta)$, is absent from both models ([](#und-factor)).
- With the default propagator the red shift is correct to first order in $x$. The term of second order requires `field.propagator = "spectral"` ([](#und-angle)).
- The harmonics above the fundamental were not examined.
- The polarization away from the axis was not examined. A planar run has one transverse component of the field, and undulator radiation has a second component away from the planes of symmetry.
- The Metal [device](glossary.md#gl-device) was not run.
