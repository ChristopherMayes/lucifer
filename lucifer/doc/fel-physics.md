---
title: "Lucifer, the Bmad FEL Tracker: Physics Manual"
short_title: Physics manual
---

(sec-intro)=
## Introduction

This manual states the physics and the numerical conventions of the Bmad FEL tracker.
The code references these sections by label. The validation harness
(`lucifer/tests/`) proves each statement against Genesis 1.3 Version 4 (Genesis4), against
closed forms, or against independent implementations, and the measured levels live in
`lucifer/doc/validation.md`. Each section ends with a *Provenance* note (the Genesis4
routine the physics was transcribed from, moved here from the source comments) and a
*Validation* note (which check pins it). The physics was transcribed from Genesis4
(GPL permits transcription), validated at transcription level, and then re-expressed in
Bmad's units and conventions. The transcription-era agreement is banked in `validation.md`.

Throughout: $m_e c^2$ is Bmad's `m_electron` (eV), $c$ is `c_light`,
$Z_0 = \mu_0 c$ (`mu_0_vac*c_light`), $e$ is `e_charge`, and
$\varepsilon_0$ is `eps_0_vac`. Genesis4 carries its own truncated constants
(EG $Z_0 = 376.73$ against $376.7303\ldots$). The resulting $8.3\times10^{-7}$
relative difference is the accepted floor of every Genesis4 comparison
([](#sec-numerics)).

(sec-coords)=
## Coordinates and conventions

(sec-chart)=
### The packed chart

Particles are stored in packed structure-of-arrays slices (`fel_slice_struct`),
in Bmad's phase-space coordinates plus a weight:

$$
  \bigl(x,\; p_x/p_0,\; y,\; p_y/p_0,\; z,\; p_z\bigr), \qquad
  w = \text{macroparticle charge [C]},
$$

with $z = -\beta c\,(t - t_{\mathrm{ref}})$ and $p_z = (P - p_0)/p_0$. The single stored
normalization reference is $p_0 c$ in eV (`fel_beam_struct%p0c`). Everything
kinematic derives from it exactly (`fel_gamma_of`, `fel_beta_of`):

$$
  \hat P = \frac{p_0}{m_e c}(1+p_z), \qquad
  \gamma = \sqrt{\hat P^2 + 1}, \qquad
  \beta = \hat P/\gamma, \qquad
  \tau \equiv -z/\beta = c\,(t - t_{\mathrm{ref}}).
$$ (eq-kinematics)


The ponderomotive phase is *derived, never stored*. Genesis4's per-particle
$\theta$ splits into a common reference advance (one scalar per beam, $\varphi_0$,
advanced once per step) and the particle lag carried by Bmad's $z$:

$$
  \theta_j \;=\; \varphi_0 + k_s\, z_j/\beta_j \;=\; \varphi_0 - k_s\,\tau_j ,
  \qquad k_s = 2\pi/\lambda_s .
$$ (eq-theta)

This reference-offset split is the FP32-safe formulation and removes
the wrap hazard outright: $z$ does not wrap, so slice migration
([](#sec-migration)) needs no $\theta$-wrap protocol. The common advance rate over
an element with reference wavenumber $k_u^{\mathrm{like}}$ is computed
cancellation-free (`fel_phi0_rate`):

$$
  \begin{aligned}
    \frac{d\varphi_0}{ds} &= k_u^{\mathrm{like}} + k_s\Bigl(1 - \frac{1}{\beta_0}\Bigr),\\[4pt]
    1 - \frac{1}{\beta_0} &= \frac{-1}{\beta_0\,\gamma_{0b}^2\,(1+\beta_0)},
      \qquad \gamma_{0b}^2 = \hat P_0^2 + 1 .
  \end{aligned}
$$ (eq-phi0rate)

Inside undulators $k_u^{\mathrm{like}} = k_u$. Over field-free interludes Genesis4's
drift surrogate $k_u^{\mathrm{like}} = k_s/(2\gamma_0^2)$ is used
([](#sec-interlude)).

(sec-window)=
### The time window

A time-dependent beam is a set of slices spaced by
$\Delta = \texttt{n\_wavelength}\cdot\lambda_s$, the deck's `slicing%n_wavelength`.

A deck that states no spacing has it derived from the gain. The cooperation length
$\ell_c = \lambda_s/(4\pi\rho)$ is the distance radiation slips forward in one gain
length, so it is the scale over which the beam's slices are not independent, and the
spacing is the whole number of wavelengths in a quarter of it: four slices per
cooperation length, the resolution at which FLASH1's pulse energy converged. The Pierce
parameter comes from the beam the deck describes and the first FEL element, the same
expression the source filter's edge uses, computed before the beam is sliced because the
spacing has to be known to slice it. Two descriptions are not resolutions and keep one
wavelength: the steady state, where the whole charge sits in the one slice so the spacing
states the current, and a run whose beam is loaded rather than described, which has
nothing to derive from. The run prints which of the three it took.

The default window holds every particle of the bunch and the line's slippage, one
wavelength per undulator period, ahead of the window head. Radiation slips toward higher
slice index, so the added length is at the high-index end and the bunch sits low in the
window it was given. A Gaussian has no last particle, so its extent is the length beyond
which a slice would hold less than one electron of the charge, which scales with the
charge where a fixed number of sigmas does not. A flat bunch takes its grid extent. A
loaded bunch spans its own particles and gets no slippage added, because a loaded window
is data: it came from the run or the code that wrote the dump, and widening it would put
this program's beam in a different window from the one a comparison starts in. A stated
window says its own headroom, centers the bunch it is given, and warns when it clips it.
The window belongs to neither the beam nor the field: the beam's slices and the field's
longitudinal samples are one partition, so `slicing` states it once and both read it. A
wavefront carries only the spacing $\Delta$ that follows, since a wavefront travelling
through mirrors and gratings has no slices. `n_wavelength` is an integer from the deck to
every consumer, never recovered by dividing $\Delta$ by $\lambda_s$, which returns 12 to
the last bit rather than 12 and which the device's exact bucket arithmetic then refuses.
The one place it is recovered is a beam or field read from a dump, since no dump format
records it, and there `fel_n_wavelength` refuses a spacing that is not whole. *Higher slice
index is the window head* (earlier arrival, larger Bmad $z$). The global window
position of a particle in slice $i$ is

$$
  z_{\mathrm{global}} = z_{\mathrm{local}} + \beta\,(i-1)\,\Delta ,
$$ (eq-zglobal)

the invariant that slice migration preserves ([](#sec-migration)) and that the
whole-window wake concatenation uses directly ([](#sec-seamwake)). A single-slice
beam is the steady state: slippage is a no-op and the window degenerates to one
$\lambda_s$-periodic slice.

$\Delta$ is a light-travel distance, not a separation in Bmad's $z$, and the difference
is not a convention. The $\beta$ in Eq. [](#eq-zglobal) is the particle's, and with
$z = -\beta c\,(t - t_{\mathrm{ref}})$ it cancels at a grid point:

$$
  t - t_{\mathrm{ref}} = -(i-1)\,\Delta/c ,
$$ (eq-tslice)

so the grid is exactly uniform in arrival time and in $ct$, and only $\beta$-dependently
uniform in $z$. That is why the concatenation stores every particle's entry $\beta$
rather than one number. It is also what makes slippage an exact integer shift: the light
advances one $\lambda_s$ per undulator period, at the rate
$(1+a_w^2)/(2\gamma^2\lambda_s)$ per unit length, so one slice is exactly `n_wavelength`
periods of slippage and the field record rotates by one index with no interpolation
([](#sec-slippage)). Uniformity in $ct$ is the reason that works. The file therefore
publishes the two exact coordinates, $t$ and $ct$, and Bmad's $z$ only at the reference
$\beta$, which is the one of the three that needs a reference at all
([](#sec-stats)).

(sec-units)=
### Units

The radiation field is stored in V/m (the `wavefront_struct` convention). There
is no internal unit system. Beam energies in diagnostics are total energies
$E = \gamma\, m_e c^2$ in eV, which is Bmad's convention. Here *energy never
means* $\gamma$. Genesis4's `energy` output *is* $\gamma$, and the comparison
scripts convert. Genesis4's internal field unit $u$ relates to the physical field by
$u = E\,k_s/(\sqrt2\, m_e c^2)$, the identity under which every formula in this manual
was converted from its internal-unit original. The converted forms are simpler
([](#sec-eom), [](#sec-field)).

:::{admonition} Provenance
:class: note
The chart and window conventions are this port's own.
Genesis4 stores $(\theta,\gamma)$ per particle with per-slice bounded $\theta$. The
unit identity composes Genesis4's dump scale (`writeFieldHDF5.cpp`) with
$E = \mathit{dfl}\,\sqrt{2 Z_0}/\mathit{dgrid}$.
:::
Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-element)=
## The FEL element

An FEL segment is a Bmad `wiggler`/`undulator` whose `fel_method` is not `off`,
recognized by key and method, never by name. Its FEL
parameters derive from lattice attributes:

$$
  K = \frac{c\,B_{\max}}{k_u\, m_e c^2/e c} = \frac{c\,B_{\max}}{k_u}\frac{e}{m_e c},
$$ (eq-awmap)

$$
  a_w = \begin{cases}
    K          & \text{helical (\texttt{field\_calc = helical\_model})},\\[2pt]
    K/\sqrt2   & \text{planar (\texttt{planar\_model})}.
  \end{cases}
$$

Here $a_w$ is rms, Genesis4's convention. In code $K = c\,\texttt{b\_max}/(k_u\,
\texttt{m\_electron})$ with `b_max` in Tesla and `m_electron` in eV. The
natural-focusing split follows the helicity defaults, scaled by $k_u^2$:
helical $k_x = k_y = \tfrac12 k_u^2$. Planar $k_x = 0$, $k_y = k_u^2$. The integration
step is the element's own `ds_step`/`num_steps`. The bookkeeper's
$n_{\mathrm{step}} = \max(1, \mathrm{nint}(L/\texttt{ds\_step}))$ is exactly Genesis4's
unroll.

An element that states neither has its step derived from the gain, as the slice spacing
is ([](#sec-window)): a twentieth of the one-dimensional gain length
$L_{1D} = \lambda_u/(4\pi\sqrt{3}\rho)$, rounded down to whole periods and at least one.
Twenty steps per gain length is where the Aramis step sweep puts the saturated power
within 3 percent, and one period per step against three landed inside the seed spread on
FLASH1 and reproduced LCLS's gain length to two decimals. What it replaces is not silence
but Bmad's own default for a wiggler, $\lambda_u/20$, which resolves the wiggle this
model has already averaged over: the FLASH-scale probe ran at 3300 steps a segment
against FLASH1's 55 for that reason. The unset state is that exact quotient, since the
bookkeeper leaves nothing else to read, so an element whose step is $\lambda_u/20$ by the
deck's own choice is derived over. Any other stated step is kept.

The mode is the element's own `fel_method` attribute, `averaged` or `unaveraged`,
class-settable as `wiggler::*[FEL_METHOD] = ...`, so averaged and unaveraged segments mix
freely in one line. `averaged` is the averaged model on the `bmad_standard` kernel's
transverse maps ([](#sec-undbmad)) and `unaveraged` is [](#sec-unaveraged). Bmad carries
the attribute, as it carries `space_charge_method`, so a lattice, the parser and `show
ele` say the same thing, a misspelling is refused rather than read as an integer, and
every other Bmad program loads the lattice. The element's `tracking_method` is a separate
axis naming how one particle crosses it, and an FEL segment leaves it at `bmad_standard`.
The transcribed maps are still selectable, as the run switch
`global%transport_model = "genesis"`, since they are chosen for a whole comparison rather
than per element. The unaveraged mode's own two numbers are run switches for the same
reason, `global%unaveraged_steps_per_period` (20, below 10 refused) and
`global%unaveraged_ramp_periods` (2). A silent hard edge would reintroduce the
$K/\gamma$ handoff hazard, so a true hard edge (the test configuration) must be asked for
explicitly with the sentinel $-1$. Silence never means hard edge. The ramps carry two priced end
effects ([](#sec-unaveraged)): their slippage deficit is compensated exactly by
the mode's built-in end-of-segment phase jump, and their reduced coupling length
costs $\sim$2% in $\ln P$ per 266-period benchmark segment at the default 2
periods. That is real field physics of adiabatic ends, scaling with the ramp length,
deliberately left visible. `tilt` is honored on planar elements as the polarization specification
([](#sec-vector)): the wiggle plane, the focusing, the coupling roll-off and the
unaveraged field all rotate with it, and the radiation grows a second component.
Refused at parse time: missing `b_max` or `l_period`, a
fieldmap `field_calc`, nonzero `kx` (unmapped roll-off), an unset
`ds_step`, `tilt` on a helical element (a rotation of a circularly
symmetric field, a no-op that reads as a mistake) and `tilt` with the
transcribed-Genesis4 maps (Genesis4 has no tilt). Outside the FEL walk the element is a plain periodic wiggler: the
`track1_custom` and `make_mat6_custom` hooks delegate to
`track_a_wiggler`, so reference time and transfer matrices come from Bmad's own
kernel.

:::{admonition} Provenance
:class: note
Helicity defaults from Genesis4's `LatticeParser` ($k_x{=}k_y{=}0.5$ helical,
$0/1$ planar, scaled by $k_u^2$ in the unroll); the step rule from
`Lattice::unrollLattice`.
:::
Measured levels and how they are checked: [](validation.md#val-the-programs-own-identities).

(sec-core)=
## The FEL core

(sec-step)=
### Step composition

One integration step of length $\delta z$ inside an undulator advances, in order
(`fel_track_und_step`): (i) transverse half step, (ii) longitudinal RK4 with
the field gathered once and held, (iii) the wake's energy decrement
([](#sec-wakes)), (iv) transverse half step, (v) source deposition and field
propagation ([](#sec-field)), (vi) slippage, applied by the caller
([](#sec-slippage)). $\varphi_0$ advances once per step by Eq. [](#eq-phi0rate)
with $k_u^{\mathrm{like}} = k_u$. Beam slice $i$ couples to field record slice
$1 + \mathrm{mod}(i-1+\texttt{first},\, n_{\mathrm{slice}})$ ([](#sec-slippage)).

(sec-coupling)=
### Coupling and field roll-off

The coupling factor at harmonic $h$ (`fel_und_coupling`):

$$
  f_c(h) = \begin{cases}
    a_w
      & \text{helical, } h = 1; \quad 0 \text{ otherwise},\\[4pt]
    a_w\,(-1)^{(h-1)/2}\Bigl[J_{\frac{h-1}{2}}(\xi) - J_{\frac{h+1}{2}}(\xi)\Bigr]
      & \text{planar, odd } h; \quad 0 \text{ even},
  \end{cases}
$$ (eq-fc)

with $\xi = (h/2)\,a_w^2/(1+a_w^2)$. Only $h=1$ is used (harmonics beyond the
coupling formula are out of scope). The transverse roll-off of the undulator field
about its offset $(a_x, a_y)$:

$$
  \begin{aligned}
    f_{aw}(x,y)   &= 1 + \tfrac12\!\left[k_x (x-a_x)^2 + k_y (y-a_y)^2\right],\\[2pt]
    f_{aw}^2(x,y) &= 1 + k_x (x-a_x)^2 + k_y (y-a_y)^2 ,
  \end{aligned}
$$ (eq-faw)

the first-order form in the particle gather, the squared form in the deposition.

(sec-eom)=
### Longitudinal equations of motion

Per particle, with $\hat p_{x,y} = \gamma\beta_{x,y}$ (in code
$p_{x}\cdot p_0/m_ec$),

$$
  b_\perp = 1 + \hat p_x^2 + \hat p_y^2 + a_w^2 f_{aw}^2(x,y),
$$ (eq-btpar)

and the complex drive from the local field $E$ (bilinear interpolation from the
particle's cell, zero off-grid),

$$
  r = \frac{f_c(1)}{\sqrt2\, m_e c^2}\, f_{aw}(x,y)\, E^{*} \quad [\mathrm{m}^{-1}],
  \qquad \hat c = r\, e^{-i\theta},
$$ (eq-rpart)

the equations of motion integrated by RK4 (`fel_ode`) are

$$
\begin{aligned}
  \beta_\parallel &= \sqrt{1 - \frac{b_\perp - (2/k_s)\,\mathrm{Re}\,\hat c}{\gamma^2}},\\
  \frac{d\theta}{dz} &= k_s\Bigl(1 - \frac{1}{\beta_\parallel}\Bigr) + k_u,\\
  \frac{d\gamma}{dz} &= \frac{\mathrm{Im}\,\hat c}{\beta_\parallel\,\gamma} - e_z,
\end{aligned}
$$ (eq-dgamma)

with $e_z$ the space-charge term ([](#sec-spacecharge)) in units of $m_e c^2$ per
meter, held fixed through the stages along with $r$, $\hat p_{x,y}$ and $f_{aw}$. The
integrator (`fel_runge_kutta`) is the classical fourth-order Runge--Kutta in
Genesis4's incremental bookkeeping, kept verbatim so the audited transcription
survives.

*Chart round trip.* $(\theta,\gamma)$ are derived at step entry from
Eqs. [](#eq-kinematics)--[](#eq-theta) and written back at exit:
$\hat P = \sqrt{\gamma^2-1}$, $p_z = (\hat P m_e c - p_0)/p_0$ (exact once $\hat P$ is
formed), $z = -\beta\,(\varphi_0^{\mathrm{new}} - \theta)/k_s$ with the updated
$\beta$. The $\theta \leftrightarrow z$ map is affine with a state-independent shift,
under which RK4 is exactly invariant. The $\gamma \leftrightarrow p_z$ round trip
costs $\sim$1 ulp of $\gamma$ per step ([](#sec-numerics)).

:::{admonition} Provenance
:class: note
`BeamSolver::advance` (gather, $r$, stage order), `BeamSolver::RungeKutta`
(verbatim), `BeamSolver::ODE`; `Undulator::fc`, `faw`, `faw2`;
grid weights from `Field::getLLGridpoint`. The V/m forms follow from the unit
identity of [](#sec-units), under which the rest energy and $k_s$ cancel out of
$r$.
:::
Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-transverse)=
## Transverse motion

(sec-natfocus)=
### Natural focusing inside undulators

The default in-undulator transverse transport (`fel_transverse_track`) is a
drift-plus-natural-focusing map over each half step. With the per-particle
longitudinal factor in stored normalization,

$$
  \hat\gamma_z = \gamma_z\frac{m_ec}{p_0} = \sqrt{(1+p_z)^2 - (p_x/p_0)^2 - (p_y/p_0)^2
   - (a_w m_ec/p_0)^2},
$$ (eq-gzhat)

the effective focusing strengths are

$$
  \hat q_{x,y} = \frac{k_{x,y}\, a_w^2}{\gamma_0\,\beta_{\parallel 0}}\frac{m_ec}{p_0},
  \qquad
  \beta_{\parallel 0} = \sqrt{1 - \frac{1 + a_w^2}{\gamma_0^2}},
$$ (eq-qnat)

and each plane maps by (`fel_apply_focus`, $\Omega = \sqrt{\hat q/\hat\gamma_z}\,
\delta z$)

$$
  \begin{pmatrix} x \\ p_x/p_0 \end{pmatrix} \mapsto
  \begin{pmatrix} \cos\Omega & \dfrac{\sin\Omega}{\sqrt{\hat q \hat\gamma_z}} \\[6pt]
    -\sqrt{\hat q\hat\gamma_z}\,\sin\Omega & \cos\Omega \end{pmatrix}
  \begin{pmatrix} x \\ p_x/p_0 \end{pmatrix},
$$ (eq-focusmap)

with $\cosh/\sinh$ for $\hat q<0$ and the plain drift $x \mathrel{+}= p_x\,\delta z/
(p_0\hat\gamma_z)$ for $\hat q = 0$. Chromaticity is through $\hat\gamma_z$
(Genesis4's $\gamma_z$ scaling).

(sec-interlude)=
### The Genesis4 interlude model

Field-free elements (drifts, quadrupoles) tracked "Genesis4's way" are one
integration step of the full element length: transverse half (quad strength
$\hat q = k_1\gamma_0\, m_ec/p_0$, Genesis4's chromatic convention $k_1\gamma_0/
\gamma_z$ against Bmad's $k_1 p_0/P$), then the $\theta$ advance with the path-length
term sampled once at mid element,

$$
  \Delta\theta = L\left[k_s\Bigl(1-\frac{1}{\beta_{\parallel}}\Bigr)
   + \frac{k_s}{2\gamma_0^2}\right],\qquad
  \beta_\parallel = \sqrt{1 - \frac{1 + \hat p_x^2 + \hat p_y^2}{\gamma^2}},
$$ (eq-interludetheta)

then the wake decrement and the second transverse half. With space charge off the RK4
collapses exactly (the slope is $\theta$-independent and $\gamma$ constant) and the
collapsed step is kept for bit-identity. With space charge on the full RK4 runs with
$r=0$ and the per-particle $e_z$. This model exists to *price* what the Bmad seam
changes: the seam integrates the path-length term exactly through the quad map where
this samples it once. That is a measured $\sim\!2\times10^{-3}$ rad of bunching phase
per quadrupole, compounding to the documented tier2_bmad difference.

(sec-undbmad)=
### The Bmad-kernel default

The default transverse maps (`transport_model = "bmad"`, [](#sec-element)) are
Bmad's own: the `bmad_standard` periodic-wiggler kernel, flattened to per-step
granularity. The transcribed maps of [](#sec-natfocus) remain selectable as
`transport_model = "genesis"`, validation-internal. The Genesis4 comparisons require
transcription-level transport. No production lattice writes it. The kernel, per
particle: $k_1 = -\tfrac12 g_{\max}^2/
(1+p_z)^2$ with $g_{\max} = K k_u\,m_ec/p_0$ ($K$ from Eq. [](#eq-awmap)), quadrupole
bodies via `quad_mat2_calc`, and half-octupole edge kicks
$\delta p_y = k_3 L\,(1{+}p_z)\, k_u^2\, y^3/3$ (and $x$ for helical) with
$k_3 L = 2\,\delta z\, k_1$. The kernel's own $z$ bookkeeping is deliberately absent.
The phase evolution lives entirely in [](#sec-eom), and applying it here too
would double count. Model differences from [](#sec-natfocus) (chromaticity
$1/(1+p_z)^2$ vs $1/\hat\gamma_z$, the octupole edges) are priced, not argued: power
$5\times10^{-5}$, intensity $7.3\times10^{-4}$, $+3.8\%$ runtime over the full TD
line.

:::{admonition} Provenance
:class: note
[](#sec-natfocus)--[](#sec-interlude): `TrackBeam::track`,
`applyDrift`/`applyFQuad`/`applyDQuad`, and `BeamSolver::advance`'s
drift path with the division order kept. [](#sec-undbmad):
`track_a_wiggler.f90` tracking locals (never the stored `k1x`/`k1y`
attributes, whose helical sign disagrees).
:::
Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-field)=
## The field solver

Each field slice advances by split-step Fourier propagation with source addition
(`fel_field_step`). The paraxial kernel, built once per grid
(`fel_field_kernel_init`) with $\Delta k = 2\pi/(N\,\mathit{dgrid})$ and
fftshift index mapping,

$$
  K_2(k_x,k_y) = -\,i\,\frac{k_x^2 + k_y^2}{2 k_s}, \qquad
  \tilde E \mapsto \tilde E\, e^{K_2\,\delta z},
$$ (eq-k2)

followed by the inverse transform and the source added in real space with the factor
2. The source deposits per particle with bilinear weights at
$(\sin\theta + i\cos\theta)\cdot p_j$, where

$$
  p_j = \sqrt{f_{aw}^2(x_j,y_j)}\;\frac{s_w\, w_j}{\gamma_j}, \qquad
  s_w = \frac{f_c(1)\, Z_0\, \sqrt2\, c\;\delta z}{4\,\mathit{dgrid}^2\,\Delta s} ,
$$ (eq-source)

the weighted V/m form in which the rest energy and $k_s$ have cancelled. Genesis4's
$\mathit{current}/n_{\mathrm{part}}$ is the uniform-weight special case of
$c\,w_j/\Delta s$. Field diagnostics: power $= \sum_{\mathrm{cells}}|E|^2\,
\mathit{dgrid}^2/(2Z_0)$, and on-axis intensity is $|E(0,0)|^2/(2Z_0)$. The window field
energy plotted by the tools is $U = \sum_i P_i\,\Delta s/c$.

The grid itself is derived from the beam when the deck states neither `grid_n_pts` nor
`grid_half_width`, and both are printed with their origin at setup. The rms beam size
$\sigma$ is taken from the emittances the deck states and the matched Twiss the lattice
states, averaged over the FEL elements by length, as the quadratic mean of the two
planes. The half width is then $9\sigma$ and the cell size $\sigma/7$:

$$
  \mathit{dgrid} = 9\sigma, \qquad dx = \sigma/7, \qquad
  N = 2\lceil \mathit{dgrid}/dx \rceil + 1 = 127 .
$$ (eq-gridrule)

The point count is 127 whatever the beam is, since $\sigma$ cancels between the two
rules, and it is rounded up from there to 128. Cells of a seventh of the beam resolve the beam and the mode, and finer cells
raise the wide-angle emission of [](#sec-source-filter) without moving the mode power.
Nine beam sizes contain the mode, and above that the half width does not enter. Both
constants were measured on three machines a hundred times apart in wavelength
([](startup-noise.md#sn-recommendations)). The rounding to a power of two is what lets
one derived deck run on the CPU and on the device at one grid, and it only ever refines,
more points across one half width being smaller cells. One of the two may be set and the other is
then derived against it: a stated half width keeps the cell rule and a stated point
count keeps the containment rule. Nothing is derived without an emittance to derive
from, so a run that starts from a dump states its own grid or takes the field's.

:::{admonition} Provenance
:class: note
`FieldSolverFFT::advance`, `FFT`, `init` (unfiltered path);
diagnostics per `DiagField::calc`. The internal-unit source scale
$f_c\,\mathit{vacimp}\,I\,k_s\,\delta z/(4\,\mathit{eev}\,n\,\mathit{dgrid}^2)$
converts to Eq. [](#eq-source) by [](#sec-units).
:::
Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-convergence)=
## Convergence of an unseeded run

An unseeded run starts the field at zero and lets the beam's shot noise seed it, and the
power it reports depends on two numbers with no physical counterpart: the transverse cell
size and the macroparticles per slice. This section is the whole of it. The mechanism and
its scaling are [](startup-noise.md), the equations of the filter that answers it are
[](#sec-source-filter), and what follows is what a user has to decide.

**What the artifact is.** Each beamlet is one transverse point on the deposition grid and
the source has no angular dependence, so a beamlet radiates into every angle the grid
carries. A real undulator radiates the resonant wavelength only within its central cone.
The difference is the coherent emission of each beamlet with itself, it goes as one over
the beamlet count, it grows as the cells shrink, and at grids and loads in common use it
is most of the power an unseeded run reports. The power inside the mode carries none of
that dependence.

**How it depends on the cell size.** The deposit alone radiates as one over the cell
size squared, since the grid carries that many transverse modes and a point beamlet
feeds them all alike. The power a run reports rises more slowly than that once the
grid's widest angle passes a knee, and the knee is set by two lengths the deck fixes
without naming them. Between deposits the solver propagates the field, and a mode at
angle $\theta$ slips against the axial wave by $k_s\theta^2\Delta z/2$ a step, which is
the paraxial form of the angular red shift. The field record slips one slice past the
beam every $s\lambda_u$ of undulator, $s$ being `slicing%n_wavelength`, so a field slice
receives the deposits of one static beam slice for that long and of an independent one
after. A wide-angle mode therefore adds coherently over one slippage of the record and
incoherently across them, and the modes beyond

$$\theta_t = \sqrt{\lambda / (s\,\lambda_u)}$$

count a fraction $\Delta z/(s\lambda_u)$ of what the modes below it count. On the Aramis
benchmark at a spacing of 12 wavelengths that angle is 23.6 µrad, the exponent of the
reported power is 2.25, 1.42 and 1.08 over cells of 6.35 to 0.78 µm, and the derived sum
gives 2.26, 1.47 and 0.90 with nothing fitted. The half width does not enter, the spacing
moves the knee as $s^{-1/2}$, and a step of one slice rotation returns the exponent to 2
([](startup-noise.md#fig-sn-cause)). The mechanism is physical and the cutoff is what is
missing: a real electron dephases at the same rate, and the red shift leaves it no power
at that angle in the undulator's band. The model gets the coherence length of the
wide-angle light right and its spectrum wrong, which is why the remedy is a cutoff on the
source and not a change to the propagation.

**What the filter does and what it costs.** It multiplies the transformed source by a
sigmoid beyond a derived edge, so the beam adds nothing at angles the undulator does not
radiate into, while the field already present propagates untouched. Measured on four
machines, two of them real, at a converged load it removes about a hundred times the
wide-angle power while moving the mode power a few percent and the saturation point not at
all. What it removes with the artifact is the rim of the physical cone, 14 percent of the
in-cone startup power where the derived edge sits well inside the cone and 28 to 30 percent
where the edge and the cone nearly coincide. Saturation absorbed that on every machine
measured. In wall clock it costs 8 to 14 percent, the source gaining a transform pair.

**Why it is on.** Because what it removes is emission no physical beam produces, and
because it is the cheaper route to a converged answer. The four configurations below were
measured on one machine at its derived grid, one seed each.

```{table} Four configurations on the Aramis benchmark at the derived grid of 3.15 um cells, a 300-slice window, one seed. The wide-angle columns are the power outside the mode over the power inside it, at the mode's saturation and at the exit.
:name: tab-convergence

| macroparticles | beamlet size | beamlets | filter | wall clock | mode power at exit | wide/mode at saturation | wide/mode at exit | saturation |
|---|---|---|---|---|---|---|---|---|
| 1024 | 8 | 128 | on | 7 s | 1.64 GW | 0.000 | 0.001 | 32.9 m |
| 4096 | 8 | 512 | off | 9 s | 0.52 GW | 0.062 | 0.641 | 32.9 m |
| 8192 | 4 | 2048 | off | 12 s | 1.67 GW | 0.020 | 0.069 | 32.9 m |
| 8192 | 4 | 2048 | on | 13 s | 1.73 GW | 0.000 | 0.000 | 32.9 m |
```

The first and the third are both converged and the first costs 60 percent of the third.
The third is the reference code's own default load, which reaches convergence without any
filter, and the fourth shows the filter changing that answer by 4 percent, inside the
fluctuation of the process. The second is the configuration to avoid: 512 beamlets without
the filter leaves two thirds of the exit power outside the mode.

**When to turn it off.** One case: a comparison against a code that carries no such
filter, or against a level recorded without one. That is what the comparison tiers do, and
their decks say so. A run with the filter off is otherwise identical arithmetic, so the
recorded levels are unmoved.

**The load, once the filter is on.** At or above 1024 macroparticles per slice in beamlets
of eight, which is 128 beamlets, the total power is converged by the criterion: the run's
own report reads under a hundredth. Below about 128 beamlets the power inside the mode
rises in the exponential regime by about a factor of two, which is each beamlet's field
acting back on itself, and the filter is not measured to remove it. With the filter off
there is no formula that carries between machines, since the estimate fitted on the first
machine was 40 times low on the third, so the test is to double the macroparticles at the
same beamlet size and see whether the mode power moves, or to read the run's report.

**Why the filter is not applied per slice.** An imported beam can carry very different
particle counts from slice to slice, and the artifact is largest where the count is
smallest, so filtering only the thin slices is a natural thought. It is refused. The
filter is physics that holds for every slice, and its cost on a well-populated slice has
no measured consequence, so there is nothing to save. Slippage carries the field from
slice to slice, so unfiltered light from the core would pass through every filtered tail
slice regardless. And an irregular load is what the resampler exists for: sample mode
redraws every slice at one beamlet count from the imported distribution, which gives every
slice the same statistics and leaves the filter one job.

**What the run says about its own load.** The header states the macroparticles in the
thinnest slice, the beamlet size and the beamlet count, and then where that stands: at the
load the filter's convergence was measured on, or under it, or, with the filter off, that
no estimate carries between machines and the two tests are a run at twice the load and the
report in the footer. It also states what the beamlet size buys, which is the harmonics the
shot noise is imposed at, `(beamlet_size - 1)/2`. A harmonic field above that would start
from a load with no noise at its own frequency, so it is refused with the beamlet size that
would carry it.

The footer then names the thin slices by index and position: those carrying charge whose
beamlet count is under 128, or whose power outside the split angle is over a tenth of the
power inside it, at the record where the mode peaked. A window ratio is a mean over the
window, so a thin tail averages away inside it, and an imported beam with irregular counts
is where a thin tail comes from. Setting `global%load_only` prints the header and stops,
which prices a deck and reads its load before any tracking is paid for
([](user-guide.md)).

(sec-source-filter)=
## The source filter

Each beamlet is a single transverse point on the deposition grid, and Eq. [](#eq-source)
has no angular dependence, so a beamlet radiates into every transverse wavenumber the grid
carries up to its Nyquist wavenumber. A real undulator of $N_w$ periods emits at the
resonant wavelength only within the central cone of half angle
$\theta_{con} = \sqrt{1+K^2}/(\gamma\sqrt{N_w})$, and emission at larger angles is
red-shifted out of the bandwidth. The difference is an artifact of the representation, and
at grids and macroparticle counts in common use it is most of the power an unseeded run
reports. What that means for a deck is [](#sec-convergence); the equations are here.

The source filter suppresses it. With `global%source_filter` on, the source of Eq.
[](#eq-source) gets a forward transform of its own and is multiplied there by a sigmoid in
normalized transverse spatial frequency,

$$
  \sigma(k_x, k_y) = \frac{1}{1 + \exp\!\big[(\rho_k - 1)/w\big]}, \qquad
  \rho_k = \sqrt{\Big(\frac{k_x}{k_{\mathrm{Ny}}\, c_x}\Big)^2 +
                  \Big(\frac{k_y}{k_{\mathrm{Ny}}\, c_y}\Big)^2} ,
$$ (eq-sigmoid)

in which $k_{\mathrm{Ny}} = \pi/\mathit{dgrid}$ is the grid's Nyquist wavenumber, $c_x$
and $c_y$ place the sigmoid's half-height per plane in units of half that frequency, and
$w$ is its width as a fraction of the edge, so the sigmoid on axis is
$1/(1 + e^{-1/w})$.

The edge is set as an angle, since it is a property of the beam and the gain and not of
the mesh. Genesis normalizes the shifted grid index by the grid size, so $c_x = c_y = 1$
is twice the Nyquist index, the angle $\lambda/\mathit{dgrid}$, off the grid, and
`source_filter_angle` divided by that is what the kernel builds from. Left unset, the
angle is

$$
  \theta_{edge} = \max\!\left(\frac{4\lambda}{2\pi\sigma},\;
                              \sqrt{\frac{2\rho\lambda}{\lambda_u}}\right) ,
$$ (eq-filter-edge)

the larger of four diffraction angles of the fundamental mode, with $\sigma$ the rms
transverse size of the loaded beam, and the angle at which the resonant wavelength
red-shifts by $\rho$. Inside the first lies the radiation that can couple to the mode.
Outside the second, emission is beyond the bandwidth of the interaction. Their ratio goes
as $1/\sqrt{z_R/L_g}$ with a coefficient of 6.8 measured to 5 percent over four machines,
so the mode angle is the larger below $z_R/L_g \approx 46$ and the $\rho$ angle above it.
None of the four machines reaches that, and the measured ratios run from 1.14 to 1.86, so
neither angle serves alone and the larger is taken. The larger is taken because
cutting into the mode loses real radiation while leaving artifact in only weakens the
filter. $\rho$ is the one-dimensional Pierce parameter in Genesis's form, from the
undulator's coupling and the peak current of the loaded beam. Both angles, their ratio
and which one won are printed at setup, with a warning when the ratio leaves 0.3 to 3.

That larger angle $\theta_p$ is the one the filter protects, and where the edge goes
relative to it is a separate question. An edge on $\theta_p$ puts the sigmoid's
half-amplitude point there, so the source keeps a quarter of its intensity at the angle
the run just derived as physical. `global%source_filter_tolerance` names the largest
fraction of source intensity the filter may take anywhere inside $\theta_p$, and the edge
is placed to hold it,

$$
  \theta_{edge} = \frac{\theta_p}{1 - w \ln\!\big[t/(1-t)\big]}, \qquad
  t = \sqrt{1 - \epsilon} ,
$$ (eq-filter-tolerance)

with $\epsilon$ the tolerance. The sigmoid passes amplitude $t$ at $\theta_p$ and
intensity $t^2$ there, so the tolerance fixes $t$. At $\epsilon = 3/4$ the logarithm is
zero, the margin is one and the edge is $\theta_p$ itself, which is the placement in use
before the tolerance existed and is its default. The 0.75 is the loss that placement
implies rather than a number chosen for it. Smaller tolerances move the edge out by 1.12
at 0.2 and 1.36 at 0.01, holding the transmission across $\theta_p$ at the cost of leaving
more wide-angle source in.

The logarithm is taken as $\ln t + \ln(1+t) - \ln\epsilon$, which is the same number since
$\epsilon = (1-t)(1+t)$ and which never forms $1 - t$. Formed directly, $1 - t$ is exactly
zero for any tolerance below about $2.2\times10^{-16}$, where $1 - \epsilon$ rounds to
one, and at a width of 0.01 such a tolerance still has a margin of 0.60 and is an edge the
run can build. The three logarithms also cancel to the bit at $\epsilon = 3/4$, so the
default edge is bitwise the edge derived before the input existed. A tolerance outside the
open interval from 0 to 1 is refused, checked before anything takes a logarithm: both ends
are singular, and at one the margin runs to positive infinity and the edge collapses onto
the axis, which a test on the margin's sign would take. A width too large for the
tolerance leaves the margin at or below zero and is refused with the width that would meet
it. Setting the angle, or the grid-relative cuts, places the edge directly and bypasses the
tolerance, and setup says which happened.

The filtered source is then added to the field in
Fourier space, before the one inverse transform:

$$
  \tilde E \mapsto \tilde E\, e^{K_2\,\delta z} + 2\,\sigma\,\tilde s .
$$ (eq-filtered-add)

The field itself is never multiplied by $\sigma$. Only what the beam adds is filtered, so
a seed, an imported field and everything already radiated propagate exactly as they do with
the filter off. The switch is on by default ([](#sec-convergence)), and a run with it off
takes the arithmetic of Eq. [](#eq-k2) unchanged, which is why a level recorded without it
is unmoved by a deck that says so.

The unaveraged mode ([](#sec-unaveraged)) deposits from the resolved motion and builds no
transformed source, and the coherent source ([](#sec-coherent-source)) has already replaced
the per-particle deposit by an analytic Gaussian that carries no wide-angle content. Neither
is refused, since the filter is the default and a deck that never asked for it must still
run. The filter is a property of each element, so it reaches the averaged elements and
leaves the unaveraged ones, and a mixed line filters the averaged part of itself. The
coherent source turns the filter off for the whole run. Both cases say so at setup, and the
harness measures them rather than reading those lines: an unaveraged line and a coherent
source give the same exit power with the switch either way, and a mixed line of two averaged
segments around one unaveraged segment moves by 1.2e-3. A non-positive cut or width is refused:
Genesis4 turns its own filter off there, and a run that asked for the filter and silently
did not get it is worse than one that stops. So is a width that leaves the sigmoid below
0.99 on axis, since such a filter attenuates the coherent source as much as the wide
angles: Genesis4's own default of 1 puts the axis at 0.73, and measured, that halves the
physical in-cone startup power and delays saturation by five to nine metres
([](startup-noise.md)). The grid-relative cuts lift that guard, since placing the edge
where Genesis4 places it means taking its width too.

The protected angle serves the statistics whether or not the filter is on. Every record
that takes the field angle moments reports `power_inside_angle`, the power within
$\theta_p$ of the axis, from the transform those moments already pay for. It does not
follow the edge out as the tolerance falls, so a scan of the tolerance is measured
against one acceptance. Its
complement to the total is the wide-angle emission of the point beamlets, which is what a
SASE run has to separate from the mode before quoting a power.

:::{admonition} Provenance
:class: note
`FieldSolverFFT::init` (the sigmoid table, lines 129-146), `FieldSolverFFT::FFT` (the
filtered branch, lines 74-100) and `FieldSolverFFT::initSourceFilter` (lines 158-170).
Genesis4 added the filter in release 4.6.12 with its FFT solver, describing it in the
changelog as suppressing strongly diffracting field components from a rough distribution of
the source term. At its defaults the edge lies beyond the grid, so its default filter is
the soft roll-off of the sigmoid alone. The exponent of Eq. [](#eq-sigmoid) is clamped at 700 here, which Genesis
does not do: it is unreachable at Genesis4's own cuts and overflows for a cut small enough
to carry the grid's corner far outside the edge, where $\sigma$ is already zero in double
precision.
:::
Measured levels and how they are checked: [](validation.md#val-source-filter).

(sec-slippage)=
## Time dependence and slippage

The field record is a circular buffer over the wavefront's slice index, rotated by
slippage rather than moved. Beam slice $i$ couples to record slice

$$
  \mathrm{ifld}(i) = 1 + \mathrm{mod}(i - 1 + \texttt{first},\, n_{\mathrm{slice}}),
$$ (eq-fieldindex)

and everything reading the record in time order (diagnostics, dumps) unrotates through
the same map. Slippage accumulates in units of $\lambda_s$ in a running total $S$
(Genesis4's `accuslip`). Each undulator step adds

$$
  \delta_{\mathrm{slip}} = \frac{\delta z\,(1 + a_w^2)}{2\gamma_0^2\,\lambda_s},
$$ (eq-slipstep)

and whenever $|S| > 0.8\cdot\texttt{n\_wavelength}$ the record rotates one slice and
$S$ drops by `n_wavelength`: the slice coupled to the window *head* is zeroed and
re-enters at the tail. Radiation leaves at the head and fresh vacuum enters behind the
bunch (a non-periodic fill). Backward slippage mirrors the direction.

Interludes rotate the record by no step of their own. The particles still advance
through them. The Genesis4 interlude model samples the path-length term of
Eq. [](#eq-interludetheta) once per element, and the Bmad seam tracks each slice
through the element's own map, so the momentum dependence of the path length,
$R_{56}$ and its transverse partners, reaches every particle's $z$ and through
Eq. [](#eq-theta) its $\theta$, while the reference phase follows the beam through
`fel_phi0_rate`. What the interludes contribute to the record is their lengths $L_z$,
accumulated until an undulator follows, when the last interlude element receives
$\lfloor L_z/(2\gamma_0^2\lambda_s)\rfloor + 1$ wavelengths of autophasing, whole
rotations of the record. The end-of-lattice fixup is *unguarded*: the final element
always receives the same $+1$, even with no trailing interlude at all. That is a
Genesis4 quirk kept deliberately, since guarding it leaves the record one rotation
short.

:::{admonition} Provenance
:class: note
`Control::applySlippage` reduced to one shared-memory node (the MPI ring exchange
is the identity); the schedule from `Lattice::calcSlippage`, including the
unguarded fixup (found live: the guarded version failed td1 at 0.84 of the final
field).
:::
Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-loading)=
## Loading

One path loads the beam. A bunch comes from Bmad's `beam_init` structure, generated from
its description or read through `beam_init%position_file`, which reads openPMD-beamphysics.
It is binned by arrival time into the slices of [](#sec-window), and `load_mode` says what
becomes of the particles that fell into a slice. Two switches then act on every load:
`quiet_start` makes it quiet ([](#sec-quiet)) and `shot_noise` imposes the physical noise
on whatever load there is ([](#sec-noise)). `beam_file` keeps its one meaning, a dump of
slices loaded as they are with `field_file`.

The arrival coordinate is $\tau$ of Eq. [](#eq-kinematics), shifted so the earliest
particle sits at zero, and the slice spacing is $\Delta s = \texttt{n\_wavelength}\,
\lambda_0$. Bmad's `init_beam_distribution` makes the bunch on both sources, so a file is
read whole: the reader's own use of `beam_init%n_particle` as a cap on the particles read is
suppressed, since here the count means something else.

### The sample mode

`load_mode = "sample"`, the default, loads the same number of macroparticles into every
slice, `beam_init%n_particle` counted after the phase copies, with weights from the slice's
charge, so the current follows the bunch. Genesis4's `&beam` with a `&profile` and its
`importdistribution` are this one mode.

*From the description.* With no file the loader evaluates `beam_init` analytically per
slice and no bunch is ever made. It honors `n_particle` (per slice),
`a_norm_emit`/`b_norm_emit`, `sig_pz` ($\delta\gamma = \beta_0 p_0 \,\texttt{sig\_pz}/mc$),
`bunch_charge`, `sig_z` and `distribution_type(3)`, and refuses every other field that is
set. A standard structure that silently dropped fields would be worse than a custom one.
The Twiss is always the lattice's ([](#sec-element)'s one-truth rule), and the current is
derived, never input:

$$
  \begin{aligned}
    I(s_i) &= \frac{Q c}{\sqrt{2\pi}\,\sigma_z}\, e^{-s_i^2/2\sigma_z^2}
      && \text{(Gaussian)},\\[4pt]
    I &= \frac{Q c}{x_{\max}-x_{\min}}
      && \text{(GRID: Bmad's uniform, flat)},
  \end{aligned}
$$ (eq-derivedcurrent)

evaluated at the slice centers with the bunch centered in the window, or `slicing%current`
where the deck states a flat current directly. $\sigma_z = 0$ is the steady state (the whole
charge in one slice window, $I = Qc/\Delta s$) and is refused for time-dependent windows.
The default window covers the described bunch ($\pm4\sigma_z$ Gaussian, the grid extent
flat). `slicing%window_length` overrides it for slippage headroom and warns with numbers
when it clips the bunch. Each slice draws $m = n_{\mathrm{part}}/b$ transverse samples
matched to the lattice Twiss, $b$ the `beamlet_size`, and one energy each.

*From a file.* With `beam_init%position_file` the bunch is resampled slice by slice by
Genesis4's method, [](#sec-import), which draws the slice's macroparticles from the
particles that fell into its sampling window and takes the current from the same window.

### The keep mode

`load_mode = "keep"` keeps every real particle. Each becomes a beamlet of $b$ copies at
weight $w/b$ sharing its five other coordinates, or stays one particle at its whole weight
without the quiet start, and the counts follow the charge: a slice holds as many
macroparticles as the bunch put there. The window is
$\lceil \tau_{\mathrm{total}}/\Delta s \rceil$ slices from the earliest particle, or
`slicing%n_slice` or `slicing%window_length` where the deck states one, the bunch centered
and the particles outside counted in a warning. No random number is drawn before the
noise, so the same bunch loads the same way every time.

The copies stand at $\theta_0 + 2\pi j/b$, $j = 0..b-1$, about the particle's own phase
$\theta_0 = k_s (s - s_i) \bmod 2\pi$ inside its slice, and $z = \beta\theta/k_s$ is the
stored chart's. The transverse momenta are exact, $P_x/(mc)$ over $p_0/(mc)$, where the
resampler carries Genesis4's $\gamma x'$ instead as transcription fidelity. Every
charge-weighted moment of the bunch is therefore the load's, slice by slice, which is
what a bunch from an accelerator simulation is loaded for. The per-slice current and the
five first moments go to `<out_root>.import.txt` at full precision, where the exactness
checks read them ([](validation.md#val-the-load-path)).

(sec-quiet)=
### Quiet start

`quiet_start = T`, the default, loads beamlets: one transverse sample, one energy, $b$
phases $\theta = \theta_0 + j\,2\pi/b$, so every bunching harmonic below $b$ cancels to
roundoff. Weights are uniform within a beamlet, since the cancellation is per beamlet. The
averaged step sees $z$ only through $\theta$, periodic in $2\pi$, so a beamlet shifted whole
stays quiet at every harmonic and particles shifted singly do not, which is what the floor
of [](#sec-noise) measures.

`quiet_start = F` loads no beamlets. The analytic generator draws every particle its own
five coordinates and its own uniform phase. The resampler keeps the copies' coordinates and
gives every phase its own draw. The keep mode keeps the particle's phase and its whole
weight. Each is a load with macroparticle noise at every harmonic, as a real bunch's would
be, and a user who quietened a bunch by other means loads it this way and lets
`shot_noise` group it.

The `beamlet_size` is a loading quantity. Its one use at run time, the coherent source's
count of independent transverse samples ([](#sec-coherent-source)), is a count of the slice's
distinct transverse positions, exact for copies, so a dump loads without it.

(sec-noise)=
### Weighted shot noise

`shot_noise = T` imposes Fawley's noise, generalized to weights, on whatever load there is.
Every group of macroparticles that shares its transverse coordinates gets one phasor per
harmonic, drawn on the group's real electron count, its charge over $e$:

$$
  \begin{aligned}
    \theta_{gj} &\mathrel{-}= \sum_{h=0}^{n_h-1} a_{gh}\,
       \sin\!\bigl((h{+}1)\theta_{gj} + \phi_{gh}\bigr),\\[4pt]
    a_{gh} &= \frac{2}{h+1}\sqrt{\frac{-\ln u_{gh}}{n_{g}}},
      \qquad n_{g} = \frac{w_g}{e},
  \end{aligned}
$$ (eq-fawley)

with $n_h = (b-1)/2$, $\phi_{gh}$ uniform in $[0,2\pi)$, and kicks accumulated from the
unperturbed phases. This makes $\langle|b(h)|^2\rangle = 1/N_\lambda$ exact for any
distribution of weight across the groups, $N_\lambda$ the slice's real electron count, and
it is what [](validation.md#val-shot-noise-under-weights) measures. Genesis4's silent
$n_{g}<1$ clamp is kept but counted and warned.

The groups are found by rule, and the rule used is printed once per load:

1. the beamlets, when this loader made them;
2. else the particles sharing their five transverse coordinates to the bit, which is what
   copies made by another loader look like and what nothing else does;
3. else the occupants of each deposit cell of the field grid, the finest structure the
   source term can see.

Before imposing, the loader measures the quiet floor, $\max_h |b(h)|^2 N_\lambda$ over
$h = 1..b-1$, and refuses `shot_noise` on a slice above $0.01$ with the value in the
message: a load that already carries its noise would count it twice. A generated load
without the quiet start sits at $N_\lambda/N_{\mathrm{part}}$, far above the floor, and a
bunch quietened outside the loader sits at roundoff and passes. The four combinations:

| `quiet_start` | `shot_noise` | The load |
|---|---|---|
| T | T | Quiet beamlets with the physical noise: the SASE start |
| T | F | Quiet beamlets and nothing else: a seeded run, or a coherent-source run |
| F | T | The load's own groups with the physical noise: a bunch quietened elsewhere |
| F | F | The bunch as it is, macroparticle noise included |

:::{admonition} Provenance
:class: note
`ShotNoise::applyShotNoise`, generalized by $n_{bl} = n_e/m \to w_g/e$ (identical for
uniform weights). Shared code between the generator, the resampler and the keep mode
(`fel_impose_noise`). Genesis4's `importdistribution` applied the noise whether or not its
`shotnoise` flag said so, a quirk this port carried until the two switches became explicit.
:::
Measured levels and how they are checked: [](validation.md#val-shot-noise-under-weights).

(sec-import)=
## The resampler

The sample mode from a file resamples the bunch into slices by Genesis4's method
(`fel_import_mod`). Slice centers sit at $(i-1)\Delta s$, with
$n_{\mathrm{slice}} = \mathrm{nint}(\tau_{\mathrm{total}}/\Delta s)$ when not given.
Per slice: every particle inside a sampling window $d s_{\mathrm{len}} =
\texttt{slice\_width}\cdot\tau_{\mathrm{total}}$ centered on the slice (strict
inequalities) is a candidate, and the slice current comes from the same window,

$$
  I_i = \frac{c\sum_{j\in\mathrm{window}} w_j}{d s_{\mathrm{len}}},
$$ (eq-importcurrent)

the weighted generalization of Genesis4's count$\cdot dQ$. The candidate set reduces to
$n_{\mathrm{part}}/b$ seeds by random deletion, or grows by Genesis4's
phase-space interpolation: normalize the five coordinates $(\gamma,x,y,\hat p_x,
\hat p_y)$ to zero mean and unit rms, pick a random parent, find its nearest
*original* neighbor under a metric whose per-coordinate weights are fresh random
draws, and place the child at the midpoint plus $\mathrm{uniform}[-1,1]$ times the
difference per coordinate. $\theta$ is refilled uniformly over one beamlet spacing,
mirrored into $b$ copies, and [](#sec-noise) imposes the noise with
$n_e = \mathrm{nint}(I\lambda_0\,\texttt{n\_wavelength}/(e c))$. The file's $p_x/p_y$ are
*slopes*. The slope-to-momentum conversion $\hat p = x'\gamma$ happens at the
copy into the candidate set. Genesis4's `match`/`center` transforms are not
ported (a Bmad lattice carries its optics, and `init_beam_distribution`
generates matched bunches), and Genesis4's `align*` inputs are parsed but dead in v4 and
not transcribed as functional. A zero-charge bunch is refused.

The resampler's own counts, `resample%n_particle_per_slice` and `resample%beamlet_size`,
serve its validation route, `resample%use_beam_init`, which resamples a bunch Bmad generated
so the transcription can be checked on a known bunch against Genesis4 importing the same
file. On the user's path the counts are `beam_init%n_particle` and `beamlet_size`, the same
two numbers the generator reads.

:::{admonition} Provenance
:class: note
`SDDSBeam::init` (window, current, ordering), `removeParticles`,
`addParticles`, `distance`; slice centers from `Time::getPosition`.
:::
Measured levels and how they are checked: [](validation.md#val-distribution-import-a-bunchstruct).

(sec-pardump)=
### Particle dumps

A run writes its beam as openPMD (`.beam.h5`) through Bmad's own
`hdf5_write_beam`, and reads the same. There is no second format and no knob to
select one: a `beam_file` that is not openPMD is refused, with the
conversion command in the message, and
`lucifer/tests/scripts/convert_genesis.py` converts a Genesis4 `.par.h5` of
[](#sec-units)'s conventions either way.

The reason is charge. openPMD's macro-charge IS the per-particle weight $w_j$, so a
weighted beam survives. Genesis4 `.par.h5` carries one current per slice,
Eq. [](#eq-importcurrent) read backwards, so a reader must divide it out uniformly, and
per-particle weights are this port's day-one difference from Genesis4. The converter
refuses to write a nonuniform-weight beam in that format rather than writing one that
would silently read back uniform.

The slice partition is `particlePatches`, the standard's own partition of a species
record: one patch per slice, in window order, an empty slice as a patch of no particles.
The patch count is therefore the window, and the file states nothing else about it. What
openPMD has no place for comes from the deck: `lambda0`, `slicing%n_wavelength` and
`beamlet_size`. Reading a dump with no `lambda0` is refused rather than defaulted,
since a wrong wavelength rescales every phase in the run. `one4one` is not stored
either: the flag asserts that every macroparticle carries one electron, which is what the
weights say.

The reference phase is folded into the file's time coordinate. Neither format has
anywhere to put $\varphi_0$, so every reader restarts it at zero, and a dump that wrote
the lag $z_j$ alone would come back with every $\theta_j$ short by $\varphi_0$: a
different phase against the same dumped field, which is a change of state and not of
bookkeeping. So a dump writes $\beta_j \theta_j / k_s$, making the file's time
$-\theta_j/(k_s c)$ and a $\varphi_0 = 0$ reader exact. Genesis4 stores $\theta$ itself and
its reader does the same fold, so the two formats agree on what a dump means.

The chart, not the file, is what could cost a digit: openPMD stores absolute momenta and
a time, so $\hat p_x$, $\hat p_y$ and $z$ pass through $P/p_0$ and
$-\beta c\,\delta t$. A check that differences energies between two runs should read
$\hat p_z$, which the file states, rather than $\gamma$, which costs an ulp of $\gamma$
against a difference many orders smaller.

:::{admonition} Provenance
:class: note
`hdf5_write_beam` and `hdf5_read_beam` (Bmad's own openPMD particle I/O,
macro-charge as `weight`, `particlePatches` added for the window);
`fel_slice_to_bunch` and `fel_bunch_to_slice` for the chart, the same
pair the seam uses every step, with `fold_phi0` set for a dump alone.
:::
Measured levels and how they are checked: [](validation.md#val-distribution-import-a-bunchstruct).

(sec-migration)=
## Slice migration

With per-particle weights, migration needs no one4one: a mover carries its own charge.
The criterion is Genesis4's in this chart: with the derived $\theta$ of
Eq. [](#eq-theta) and the slice window $[0, 2\pi\,\texttt{n\_wavelength})$,
$a_{\mathrm{tar}} = \lfloor\theta/(2\pi\,\texttt{n\_wavelength})\rfloor$ is the relative
destination (positive $\theta$ drift moves toward higher index, the head). A mover's
$z$ shifts by exactly $-a_{\mathrm{tar}}\,\beta\,\Delta s$, changing $\theta$ by
$-a_{\mathrm{tar}}\cdot2\pi\,\texttt{n\_wavelength}$: for integer `n_wavelength` the phase every
deposition sees is continuous across the move to rounding, at every harmonic $h$ too
since $h\cdot2\pi\,\texttt{n\_wavelength}$ stays a multiple of $2\pi$, and
Eq. [](#eq-zglobal) is invariant. Removal is swap-with-last with rescan. Particles
whose destination lies beyond the window are dropped *with their charge counted*
(Genesis4 discards silently). Off by default: the Genesis4 tiers compare against
non-one4one Genesis4, which never migrates.

`sample > 1` is the intended case, not a tolerated one. In a sampled window the slices
are the only represented longitudinal positions, so the slice spacing is the only move
the discretization can express, and a particle drifting within the $2\pi\,\texttt{n\_wavelength}$
window has nowhere finer to go. The sub-window position is invisible to the per-slice
current until the threshold, and that error is below the slice spacing, the window's own
resolution. Genesis4's `localSort` uses the same sample-scaled window. What migration
does require is *integer* `n_wavelength`, since the phase continuity above holds only then,
and a non-integer value is refused at the first migration pass.

:::{admonition} Provenance
:class: note
`Sorting::localSort` (criterion, swap-with-last, rescan semantics); the silent
world-edge discard replaced by counted drops.
:::
Measured levels and how they are checked: [](validation.md#val-slice-migration-under-weights).

(sec-collective)=
## Collective effects

(sec-wakes)=
### Wake kernels

Wakes act as a per-slice energy-loss rate built from three single-particle kernels at
wavelength resolution over the window, $w(s_i)$, $s_i = (i-1)\lambda_s$, in eV per
meter per electron (negative = loss).

*Resistive wall* (`fel_resistive_wall_wake`), the numerical impedance of
Bane & Stupakov (SLAC-PUB-10707) with AC (Drude) conductivity and Leontovich surface
impedance:

$$
\begin{gathered}
  \sigma(k) = \frac{\sigma_0}{1 - i k c\tau}, \qquad
  \zeta(k) = (1-i)\sqrt{\frac{k}{2\sigma(k) Z_0}}, \notag\\[2pt]
  Z_{\mathrm{round}}(k) = \frac{Z_0/(2\pi a)}{1/\zeta - i k a/2},\\[2pt]
  Z_{\mathrm{flat}}(k) = \int_0^{X}\!
     \frac{Z_0/(2\pi a)\; dx}{\cosh x\,[\cosh x/\zeta - i k a\,\mathrm{sinhc}\,x]},\\[2pt]
  w(s) = -e\,\frac{2c}{\pi}\int_0^{k_{\max}} \mathrm{Re}\,Z(k)\cos(k s)\, dk,
\end{gathered}
$$ (eq-bszk)

with Genesis4's exact numerics kept: $k \in [0, 100/s_0]$ on 1000 trapezoid intervals,
$s_0 = (2a^2/(Z_0\sigma_0))^{1/3}$, flat integral on $[0,15]$ with 20000 points. This
routine is deliberately separable (a future Bmad-proper wake source, and
[](#sec-seamwake) already feeds it to Bmad's machinery). *Geometric* (undulator gap, convolved with
$dI/ds$): $w(s) = -\frac{Z_0 c\, e}{\pi^2 a\, L_{\mathrm{gap}}}\,2\sqrt{g/2}\,\sqrt{s}$
($\times0.956$ flat). *Roughness*: $w(s) = \frac{r_r}{\pi}\frac{4}{a^2}
\frac{e}{4\pi\varepsilon_0}\,\mathrm{Re}\!\int e^{-iq\tau} K(q)\,dq$ over Genesis4's
four-segment complex-$q$ contour, $r_r = \pi^3 h^2 a/L_r^3$, $\tau = 2\pi s/L_r$. The
$s{=}0$ bin of every kernel carries half weight (the self-loading theorem).

*Convolution and application* (`fel_wake_update`, `fel_wake_apply`):
slice currents interpolate to $\lambda_s$ resolution with a zero pad past the head
(the trapezoidal density model whose half-slice head deficit is the derived bound of
the seam-wake cross-validation), convert to electrons per bin, and each slice sums
causally *toward the head* (a trailing slice collects the wake of the charge
ahead), averaged over the `n_wavelength` sub-steps:

$$
  \begin{aligned}
    \mathcal{E}_i = \frac{1}{\texttt{n\_wavelength}}
      \sum_{j=0}^{\texttt{n\_wavelength}-1} \sum_{k\ge0} \Bigl[
        &N_{i,j+k}\,\bigl(w_{\mathrm{res}}+w_{\mathrm{rou}}\bigr)_k \\
        &+ N'_{i,j+k}\, w_{\mathrm{geo},k} \Bigr]
      + \mathcal{E}_{\mathrm{ext}} .
  \end{aligned}
$$ (eq-eloss)

The kick is uniform within a slice: every particle loses $\Delta\gamma =
\mathcal{E}_i\,\delta z/(m_e c^2)$ per step, in every element (the chamber does not
end where the undulator does), with $z$ rescaled by $\beta_{\mathrm{new}}/
\beta_{\mathrm{old}}$ to hold $\theta$ fixed through the kick. The convolution is
hoisted while currents cannot change and recomputed at the migration stride.

(sec-spacecharge)=
### Space charge

*Short range* (`fel_shortrange_ez`), per slice: center and radially bin the
particles, then for azimuthal modes $m = -n_\phi..n_\phi$ and longitudinal harmonics
$l = 1..n_z$ solve the radial tridiagonal system of the harmonic potential
(cell volumes $V_1 = \pi\,dr^2$, $V_i = \pi dr^2(2i-1)$, off-diagonals
$2\pi(i-1)$. $\log$ ring terms $-2\pi m^2\ln\frac{i}{i-1}$. Outer boundary
$-2\pi N$), with the weighted source per particle

$$
  \frac{Z_0}{m_ec^2}\,\frac{c\, w_j}{\Delta s}\,\frac{1}{k_s}\;
  e^{-im\phi_j}\, e^{-il\theta_j},
$$ (eq-scsource)

and the operator scale $-\gamma_z^2/k_s^2$ ($\gamma_z^2 = \gamma_0^2/(1+a_w^2)$ inside
undulators, $\gamma_0^2$ in drifts). Each particle receives
$e_z \mathrel{+}= 2\,\mathrm{Re}\bigl[e^{im\phi}e^{il\theta} u_m^l(r)\bigr]$, in
$m_ec^2$ per meter, entering Eq. [](#eq-dgamma). *Long range*
(`fel_longrange_esc`), whole window, one node:

$$
  E_i = \frac{\Delta s}{2\pi c\,\varepsilon_0}\sum_{j\ne i}
    \mathrm{sgn}(i-j)\,
    \Bigl[1 - \frac{|d_{ij}|}{\sqrt{d_{ij}^2 + \sigma_{x,j}\sigma_{y,j}}}\Bigr]
    \frac{I_j}{\sigma_{x,j}\sigma_{y,j}},
$$ (eq-longrange)

$$
  d_{ij} = (j-i)\,\Delta s\,\gamma_z,
$$

where the size scale is the *product* $\sigma_x\sigma_y$ (an effective area, not
a variance sum: transcribed wrong once, caught at $1.7\times10^{-1}$). The
per-particle ODE term is $e_z^{\mathrm{short}}(j) - E_i/(m_ec^2)$.

*Where it acts.* Each element's own `space_charge_method` decides, and Bmad's
`csr_and_space_charge_on` gates the lattice as a whole, exactly as in any Bmad program.
An FEL element set to `slice` gets the two terms above, with the FEL slices as the bins,
so `space_charge_com%n_bin` is not consulted. A genesis-model interlude set to `slice`
gets them through the transcribed interlude step. A Bmad-seam interlude carries neither.
It tracks one slice's bunch at a time through `track1_bunch`, whose own space-charge and
CSR paths want the centroid orbit of the whole beam, which this walk does not carry, and
Bmad would return the bunch untracked with no error set. An interlude element that sets
`space_charge_method` or `csr_method` while the master switch is on is therefore refused
at setup. Nothing in an FEL element is Bmad's, and the seam's collective terms are the
chamber wake and the element wake.

A Bmad-slice implementation of the terms above, selected by `space_charge%model`, remains
a named follow-on. The default stays the transcribed solver on the physics: it carries the
microbunching harmonics ($l \ge 1$ above) and the $\gamma_z$ of Eq. [](#eq-scsource)'s
operator scale, which Bmad's slice model has no counterpart for, while what that model
adds is a transverse kick falling as $\gamma^{-3}$.

(sec-seamwake)=
### Bmad element wakes across the window

Elements carrying Bmad `sr_wake` definitions (pseudomodes, or tabular
`z_long` with binning + FFT) act across the whole window. For wake elements only,
the slices concatenate into one bunch in global coordinates (Eq. [](#eq-zglobal)) and
Bmad's machinery applies unmodified: interludes through `track1_bunch` (wake
at `ds_wake`, once per passage). FEL wigglers through one whole-window
`track1_sr_wake` kick at the step nearest mid-element, $z$ rescaled to hold
$\theta$ as in [](#sec-wakes). Bmad's conventions, pinned: `ix_z(1)` is the
bunch head at largest $z$. The pseudomode wake is $W(\delta z) = A\,e^{\Gamma\delta z}
\sin(2\pi\phi + k\,\delta z)$ with self-kick $W(0)/2$. The `z_long` table is
positive-decelerating, causal side $\delta z<0$ (the [](#sec-wakes) kernels are
stored as signed loss, so flip the sign when bridging). Wakes on superimposed/split
elements live on the *lord* and resolve through `pointer_to_wake_ele`.
`sr%z_max` and the table extent must cover the window (refused
otherwise). Lr wakes are refused. `chamber_wake%write_kernels` exports the
[](#sec-wakes) kernels for building matching tables.

:::{admonition} Provenance
:class: note
[](#sec-wakes): `Wake::singleWakeResistive`, `singleWakeGeometric`,
`singleWakeRoughness` + `KernelRoughness` + `TrapIntegrateRoughness`,
`Collective::update` (convolution, self-half), `Collective::apply`.
[](#sec-spacecharge): `EFieldSolver::shortRange`, `analyseBeam`,
`constructLaplaceOperator`, `tridiag`, `longRange`; the ODE entry is
`BeamSolver.cpp`'s `ez = getEField + eloss`. [](#sec-seamwake): Bmad's
`wake_mod.f90` used as-is; the concatenation is this port's.
:::
Measured levels and how they are checked: [](validation.md#val-bmad-element-wakes-across).

(sec-seam)=
## The seam

Outside FEL elements, each slice converts to a Bmad bunch by plain copies (the packed
chart *is* Bmad's, and the element's $p_0c$ is asserted against the beam's),
`track1_bunch` tracks it, and the field propagates by `wavefront_drift`.
The only phase bookkeeping is one $\varphi_0$ advance per element
(Eq. [](#eq-phi0rate) with the drift surrogate $k_s/(2\gamma_0^2)$). The reference
energy is the lattice's ($\gamma_0 = E_{\mathrm{tot}}/m_ec^2$ from the beginning
element, and there is no namelist duplicate), and every FEL parameter lives on the
lattice ([](#sec-element)), one specification of one truth.

Measured levels and how they are checked: [](validation.md#val-validation-from-one-command).

(sec-unaveraged)=
## The unaveraged mode

This mode (`fel_unaveraged_mod`, selected per element by
`fel_method = unaveraged`, [](#sec-element)) integrates the
particles through the undulator's *real* field with the radiation as a
co-evolving kick. There is no period averaging and no resonance approximation. It
keeps the grid field and the Lorentz force where MINERVA (the published existence
proof of this physics) uses modal fields. The cost is $\sim$30$\times$ per step,
priced and not hidden. What it buys is physics the averaged map cannot reach: the
full quiver dynamics, the energy accounting the beam actually pays, coupling that is
polarization-agnostic, and arbitrary harmonic content.

The mode is also an independent check on the averaged path. Every input of the averaged map (the coupling
$f_c(h)$, the roll-off $f_{aw}$, the entry/exit behavior) is an output of the
underlying unaveraged dynamics, so here those inputs become measurements. The two
paths share no approximation. Because the kick/deposit pair resolves each particle's full quiver
current, a beam carrying *physical shot noise* radiates it continuously into the
SVEA band and pays for it (self-interaction back-reaction). That is spontaneous-emission
physics the period-averaged model does not track, since it injects shot noise once, at
load time. On the saturation demo this lifts the unaveraged SASE curve $\sim$2%/m above
the averaged codes and deepens the beam's energy loss. The effect is independent of
particle count and steps per period (physics, not statistics or resolution), and
quiet-start probes are blind to it by design (the beamlet loading cancels shot noise).
Its absolute calibration against the analytic undulator spectrum is the named future
check. Like the averaged path it is parallel over slices with
bit-identical results at any thread count (per-slice private state, fixed-order
energy reduction. The harness runs a 1-thread/8-thread byte comparison). Nothing
from the averaged coupling path may
appear in it (a harness grep enforces this). The mode measures its agreement with
Eq. [](#eq-fc), and does *not* take it as an input.

*The step.* Per substep $\delta$ ($\lambda_u/20$ by default, with a floor of 10 that
this mode's own convergence measurement supports and the published 10 to 30 envelope
corroborates), a Strang split: half magnetic push, radiation kick and source
deposit at the midpoint, half magnetic push, then the source landed on the record the kick read and the pair carried through the shared pure diffraction
(Eq. [](#eq-k2)) and the $+2\cdot$source convention of [](#sec-field). The
magnetic push is classical RK4 on the exact $z$-ODEs in kinetic variables
$u = \gamma\beta$,

$$
  \frac{du_x}{ds} = b_y - \frac{u_y b_z}{u_s},\quad
  \frac{du_y}{ds} = -b_x + \frac{u_x b_z}{u_s},\quad
  \frac{d\tau}{ds} = \frac{\gamma}{u_s} - \frac{1}{\beta_0},\quad
  u_s = \sqrt{\gamma^2 - 1 - u_\perp^2},
$$ (eq-unavgode)

with $\gamma$ exact ($B$ does no work, so $\gamma$ changes only in the kick), through
$b = \nabla\times a$, $a = eA/m_ec$:
planar $a_x = a_0\, g(s) \cos(k_u s)\cosh(k_u y)$, $a_0 = \sqrt2\,a_w$. Helical
$(a_x, a_y) = a_0\, g(s)(\cos, \sin)(k_u s)\,(1 + k_u^2 r^2/4)$, $a_0 = a_w$. These
are near-axis models whose ponderomotive-average focusing reproduces
[](#sec-natfocus)'s split exactly. The envelope $g(s)$ is a $\sin^2$ ramp over
$n_{\mathrm{ramp}}$ periods at each end, continuous in amplitude *and* slope
with the $g'$ terms retained in $b$: the quiver builds adiabatically and vanishes at
the segment ends, which is where the averaged and unaveraged momentum conventions
coincide. The beam carries a convention flag (`quiver_in_px`), set inside the
region and asserted by every averaged-physics and seam entry: a hard-edge handoff
injects $K/\gamma$ of spurious transverse momentum.

*Kick and source.* With the optical phase $\Psi = (\varphi_0 - k_u s) - k_s\tau$
(the reference $\varphi_0$ advances at the averaged rate of [](#sec-chart), so
every diagnostic stays comparable, and the $k_u s$ term is undone because the particles
now actually ride $\cos k_u s$), the work phasor $W = -i\hat E\, e^{i\Psi}$ and the
polarization-basis current $\hat\jmath = u_x$ (planar), $(u_x - iu_y)/\sqrt2$
(helical):

$$
  \frac{d\gamma}{ds} = -\frac{\mathrm{Re}[W \bar{\hat\jmath}\,]}{u_s\, m_e},
  \qquad
  \mathrm{src} \mathrel{+}= i\,e^{-i\Psi}\,\hat\jmath\;
     \frac{Z_0\, c\,\delta}{2\,\Delta_{\mathrm{grid}}^2\,\Delta s}\,\frac{w}{u_s},
$$ (eq-unavgkick)

the averaged deposit of Eq. [](#eq-source) with the coupling removed and the actual
quiver current in its place. The $1/u_s$ (where the averaged solver carries
Genesis4's $1/\gamma$) is a merit choice: kick and source are exact duals of one
wave equation, and using the same $u_s$ in both (same operands, same bilinear
weights, and the source landing on the record the kick read before the unitary
diffraction, $E' = D(E + 2S)$) makes the coherent energy exchange
cancel in $E+U$ *identically*, leaving only physical spontaneous emission and
rounding in the ledger. The order matters: with the source added after the
diffraction, $E' = DE + 2S$, the beam is charged the cross term against $E$ while the
record gains it against $DE$, and the difference $4\,\mathrm{Re}\langle (D-I)E, 2S\rangle$
is first order in the substep's diffraction phase, 1.2 percent of the turnover at 20
substeps a period on a seed of 100 um waist and halving with the substep. The stored
record is the substep's midpoint field, the same Strang split read at the other half
step. Measured, this tightened the ledger 65$\times$ over the
$1/\gamma$ convention. The period-averaged limit shifts by
$O(1-\beta_\parallel) \sim 5\times10^{-9}$. Period-averaging the pair reproduces
[](#sec-coupling): the $JJ$ Bessel factor emerges from the figure-8.

*The energy ledger* (check zero, written per record step to
`<out_root>.ledger.txt`): conservation of

$$
  E_{\mathrm{beam}} + U + U_{\mathrm{esc}} - U_{\mathrm{sp}} + E_{\mathrm{rad}},
$$ (eq-ledger)

$$
  E_{\mathrm{beam}} = \sum_j w_j\,(\gamma_j - \gamma_0)\, m_e, \qquad
  U = \sum_{\mathrm{slices}} P\,\Delta s/c,
$$

with $E_{\mathrm{beam}}$ stored *relative* to the reference so per-record
changes are not differenced off a large baseline at its own summation-rounding floor
(the two-pass-variance lesson of [](#sec-numerics)), and $P$ the power of
[](#sec-field). The kick-side change is written alongside as an internal
cross-check. The last two columns make the time-dependent ledger close
*exactly*: $U_{\mathrm{esc}}$ is the energy transmitted out of the window by
slippage, banked at each zero fill (the window is an open system, so light escapes
forward), and $U_{\mathrm{sp}}$ is the cumulative deposit energy
$\sum|\Delta E_{\mathrm{src}}|^2$, the one field-energy term the kick/deposit
duality does not charge to the beam (physically, the substep's spontaneous emission,
as in Genesis4 and MINERVA). Both vanish in steady state. $E_{\mathrm{rad}}$ is the energy the beam radiated away
under `bmad_com`'s radiation switches (actual drawn sums rather than expectations,
and zero with the switches off, their default). Measured closure:
$8\times10^{-6}$ of turnover on the time-dependent harness configuration (checked at
$10^{-3}$), $1.0\times10^{-5}$ steady state. The closure premise (radiation is the
beam's only exit channel) is enforced by this mode's wake and space-charge
refusals. A wake would drain beam energy to the environment with no ledger column,
which is one more reason those refusals exist. The physical field of the
scalar envelope: planar $E_x = \mathrm{Re}[-i\hat E e^{i\Psi}]$. Helical
$(E_x, E_y) = (\mathrm{Re}[-i\hat E e^{i\Psi}], \mathrm{Re}[\hat E e^{i\Psi}])/\sqrt2$,
the pair $(1, +i)$ that rotates with the electron.
Both give intensity $|\hat E|^2/2Z_0$, so [](#sec-field)'s power diagnostic is
mode-independent.

*Spontaneous emission.* The classical rate is
$d\gamma/ds = \tfrac23 r_e \gamma^2 k_u^2 a_w^2$ ($a_w$ rms, both polarizations). Bmad's
own `radiation_damping` through the same wiggler reproduces it to $10^{-4}$, and it
is the coefficient Genesis4 uses. The two FEL modes stand in
very different relations to it. The averaged (KMR) mode *emits* the radiation but
never debits the beam: its step adds $2S$ to the field while particles are kicked by
$E$, so the $4|S|^2$ part of the field energy is created (measured on a dark segment:
the field gains $134\times$ what the beam pays). Genesis4's optional `&sponrad`
module exists to add that loss by hand and is off by default. The unaveraged mode
conserves energy by construction ([](#sec-unaveraged)'s $1/u_s$ deposit), so its
beam pays for the fraction of the emission the grid can hold: angles up to the FFT
Nyquist $\theta_{\max} = \lambda/2\,dx$. Measured at the benchmark parameters that is
3.3% of the analytic rate, and the evidence that this is acceptance-limited undulator
radiation is its scaling: over a $16\times$ range in captured solid angle the loss
moves $0.84\% \to 11.4\%$, a measured $13.7\times$ against a predicted $10.5\times$.
Both modes now honor Bmad's global switches
`bmad_com%radiation_damping_on` and
`%radiation_fluctuations_on`, the same switches every Bmad tracking path
honors. Interludes get theirs through `track1`, and this covers the custom-tracked
FEL step. Damping applies $d\gamma_j = -\tfrac23 r_e \gamma_j^2 k_u^2 a_w^2
\int g^2\,ds$ per record (each particle's own $\gamma$, and the envelope integral, so
the unaveraged ramps radiate by their actual strength). Measured: the averaged mode
lands on the analytic rate at $9\times10^{-4}$, and the unaveraged mode on the
ramp+capture composite $(1 - \tfrac{5}{4}l_{\mathrm{ramp}}/L) + f_{\mathrm{cap}}$
at $2\times10^{-6}$ (its grid-captured self-field work adds to the explicit term).
Fluctuations apply the standard quantum-diffusion variance
$1.015\times10^{-27} k_u^3 a_w^2 F(a_w)\gamma_0^4$ per meter (Genesis4's fits of the
Saldin form, see [](references.md)), one draw per
beamlet exactly as Genesis4 (independent per-particle kicks would break the quiet
start's per-beamlet cancellation), drawn serially in fixed slice order so thread
count stays invisible. Genesis4 reaches the same variance with uniform$\times\sqrt3$
draws where ours are Gaussian, the physical limit. Measured: both FEL modes sit on the
analytic form at $\le 2\%$. Bmad's own `runge_kutta`+fluctuations through the
same wiggler sits 11% below it, which is the two references' $F$-convention spread,
recorded as a cross-check level, not absorbed. The radiated energy is the ledger's
$E_{\mathrm{rad}}$ term (actual drawn sums, so closure stays exact). With the
switches off (the default) every model is unchanged bit for bit.

*Measuring $f_c$.* Paired probes (12 and 20 periods, identical 2-period ramps,
small-signal seed): the difference of the energy-modulation phasors
$F = (2/N)\sum_j \Delta\gamma_j e^{i\theta_j^0}$ isolates the flat region (ramps
and their detuning cancel exactly), and

$$
  f_c^{\mathrm{meas}} =
  \frac{|F_B - F_A|\;\beta_0\gamma_0\sqrt2\, m_e}{|\hat E_0|\,\Delta L\,
        \mathrm{sinc}(\delta_h \Delta L/2)},
$$ (eq-fcmeas)

$$
  \delta_h = h k_u - k_{s,h}\frac{1+a_w^2}{2\gamma_0^2},
$$

with $|\hat E_0| = \sqrt{4 Z_0 P/\pi w_0^2}$ the Gaussian seed's on-axis envelope.
The harmonic probe is the same undulator with the field at $\lambda_1/h$ (the mode
is harmonic-agnostic), and the load is quiet there because the beamlet start
cancels every harmonic below `beamlet_size` ([](#sec-quiet)). A quadrature load
(MINERVA's) is therefore not needed for this measurement.

:::{admonition} Provenance
:class: note
The Strang split, the handoff hazard and its convention flag are this port's own.
MINERVA (Freund and van der Slot) is the published existence proof that unaveraged
FEL simulation is practical, and the published practices this mode shares with it
are the 10 to 30 steps per period envelope, the energy ledger, the $\sin^2$ ramp
treatment and kinematic phase everywhere. Only the published work was used: MINERVA's
source was not read, no MINERVA code was transcribed, and no run of it was compared
against. MINERVA precedent decided nothing here on its own. The magnetic
push is explicit RK4 on merit: fourth order is what makes $f_c$ measurable at
$6\times10^{-4}$ with 20 steps/period (second-order symplectic needs $\sim$100 to
match), and the non-symplecticity is priced by measurement: $\gamma$ exact,
emittance drift $\le 3.3\times10^{-6}$ over the longest benchmark segment
(266 periods, 5320 steps).
:::
Measured levels and how they are checked: [](validation.md#val-the-unaveraged-mode-fc).

(sec-vector)=
### Two polarizations

The radiation carries the envelope pair $(\hat E_x, \hat E_y)$ when any FEL element
is tilted (or the seed is y-polarized). Otherwise $\hat E_y$ is never allocated and
the single-component path runs with its arithmetic untouched (every tier reproduces
bit for bit). Diffraction is polarization-diagonal. Power and intensity are totals
$|\hat E_x|^2 + |\hat E_y|^2$. Each averaged element carries a polarization 2-vector
$\hat p$: planar with tilt $t$ is $(\cos t, \sin t)$, helical is $(1, +i)/\sqrt2$.
Its kick reads $E_{\mathrm{eff}} = \hat p^\dagger\!\cdot\!(\hat E_x,\hat E_y)$
and its deposit writes $\hat p\,S$, which reproduces the scalar convention exactly
when one polarization is live. The helical handedness is not a choice: the potential
$a = a_0(\cos k_u s, \sin k_u s)$ turns the velocity counterclockwise along $s$, as
Bmad's `helical_model` field does, and with a component's physical field
$\mathrm{Re}[-i\hat E_c e^{i\Psi}]$ and the carrier decreasing in time, the pair that
rotates with the electron is $(1, +i)$. The unaveraged mode's Lorentz force selects it
with no vector at all, and a seed of the other handedness leaves its beam unmodulated. The unaveraged mode needs no polarization vector at
all: the instantaneous kinetic momenta are real, and each component works against and
deposits into its own field component (the tilted element's field is the untilted
analytic potential evaluated in the rotated frame). The scalar helical envelope is
the co-rotating combination $(\hat E_x - i\hat E_y)/\sqrt2$. Measured on dark
shot-noise runs the two representations agree at $7\times10^{-15}$.

Measured levels and how they are checked: [](validation.md#val-two-polarizations-vector-radiation).

(sec-field-set)=
### The field set: harmonic radiation

The walk carries an ordered set of radiation fields (`fel_field_struct`:
harmonic number, wavefront, slippage state, escape bank, Genesis4's
`vector<Field*>`), with the fundamental always entry 1: the ponderomotive
phase, the `phi0` advance and the slippage schedule are all defined against it.
A harmonic field $h$ carries wavelength $\lambda_1/h$ and enters the physics in
exactly three places: its coupling $f_c(h)$ (the Bessel factor of
[](#sec-coupling), zero for every harmonic of a helical device), the harmonic
ponderomotive phase $h\theta$ -- live per Runge--Kutta stage in the kick and scale
the deposit phase -- and its own diffraction
kernel and escape accounting at $\lambda_1/h$. In the V/m field
convention no other factor appears: Genesis4's per-field $k_{s,h}$ factors in kick
and deposit are exactly its per-field internal unit conversion, which this
formulation absorbed once ([](#sec-units)). Every field shares the time window,
so all slippage states advance in fundamental-wavelength units and the records
rotate in lockstep on one slippage sample. Harmonic fields start dark and
grow from the bunching the beam carries, or import. A single-entry set is the
pre-harmonic walk, bit for bit: the same overlay discipline as
[](#sec-vector). The two are not yet validated together: harmonic fields
with two live polarizations, and with an unaveraged element (the unaveraged mode
carries the fundamental envelope only), are refused.

Validation (`check_harmonics.py`, ninth harness section): a planar
steady-state segment at the benchmark's $\gamma$, $\lambda$ and rms $a_w$, both
codes tracking the same Genesis4-written start with a dark third-harmonic field --
fundamental power to $5.3\times10^{-8}$, harmonic growth to $1.3\times10^{-4}$ of
its curve maximum. A one-step dark restart from a hard-bunched beam, whose exit
$P_3/P_1$ the exact deposit sum (Bessel $f_c$, $h\theta$, bilinear weights, from
the dumped particles) reproduces at $3.3\times10^{-16}$. Thread byte-identity on a
time-dependent harmonic run. Refusals matched on their messages. `examples/harmonics/` is the
runnable instance.

:::{admonition} Provenance
:class: note
The harmonic convention is Genesis4's: the harmonic ponderomotive phase enters the
kick per Runge-Kutta stage (`BeamSolver::ODE`, its $\sum_i r_i e^{-ih_i\theta}$) and
scales the deposit phase (`FieldSolver.cpp`, its $\theta_h = h\,\theta$), and every
harmonic record rotates on the fundamental's slippage sample (`Control::sample`).
:::

(sec-openpmd)=
### The openPMD wavefront

Radiation dumps speak one format, the openPMD `Wavefront` extension
(`EXT_Wavefront.md`, standard branch `upcoming-2.0.0`), written as
`.wf.h5` with a harmonic's file carrying `-h<h>`. The standard document is
authoritative: the harness verifies the written attributes against its text
independently of the writer. A `field_file` that is not openPMD is refused by
name, with the conversion command in the message, and
`lucifer/tests/scripts/convert_genesis.py` converts a Genesis4 field dump
either way. Per-harmonic imports are matched to the field whose photon energy they
carry, and a file matching none is refused.

The mapping, one decision per row where the extension's text under-determines the
file (each carried upstream into openPMD-beamphysics's own openPMD I/O):

| what | how it is written |
|---|---|
| one mesh record `electricField` | complex compound $\{r,i\}$ components `x`, `y` (V/m), both polarizations in one file. `z` never written (paraxial) |
| required attributes | on the mesh record (the extension's heading. Its body says "series", resolved for the record: `photonEnergy` is per field) |
| `photonEnergy` | $h_{\rm Planck} c/\lambda$ in joules, identifies the harmonic |
| `temporalDomain`, `spatialDomain` | `time` (field in V/m) and `r` only. The rest refused |
| slice axis | a mesh axis, since slices are simultaneous and the iteration is not time: stored order $(z,y,x)$, declared by `axisLabels`, zero-copy for Fortran and numpy alike. One iteration per file |
| harmonics | one file per harmonic (the record name is fixed by the extension) |

Validation: every required attribute and `unitDimension` checked against the
spec text in h5py. The converter's Genesis4 view of a dump agrees complex-value-wise at
$1.5\times10^{-16}$. The Fortran reader's re-dump is dataset-identical to the file it
read. The openPMD-beamphysics `Wavefront` class
agrees on the field energy through its own Genesis4 path at $1.2\times10^{-12}$.
Its `from_openpmd` reads the harmonic file to the same complex values,
`write_openpmd` round-trips exactly, and the Fortran reader takes what that
writer produces and re-dumps it dataset-identical. All of it runs every harness
pass, so the two implementations of this extension stay in step.

(sec-phasing)=
### Phasing between segments

What the beam-field phase does between undulator segments, measured before it was
built (every claim below is a scan in `check_phasing.py`). In the relative
mode (the default, `absolute_time_tracking` false) segments are
autophased: scanning an inter-segment gap by fractions of $2\gamma^2\lambda$ (one
full turn of drift slip, 25.8 mm at the benchmark) leaves the
bunching phase entering the next segment flat, measured on Genesis4 and on this code
alike, both interlude models. The mechanism is the co-moving reference convention:
`fel_phi0_rate` tracks the beam through a drift, and the window rotations are
whole wavelengths, carrier-neutral. The fractional beam-vs-light slip
$k_s L/(2\gamma^2)$ never reaches the coupling. (Genesis4's `calcSlippage`
carries commented-out code that would have applied that fraction to $\theta$.
The active code re-anchors, as the measurement shows.)

The deliberate off-phase knob is the wiggler's own `z_offset`: standard
Bmad misalignment, girder-composed, no new element or attribute. The autophase
anchor is the element's nominal reference position, so a physical displacement
$\delta$ adds the real, never-re-anchored entry phase
$\Delta\phi = -2\pi\,\delta/(2\gamma^2\lambda)$, returned at the element's last step
(the downstream break is shorter by the same $\delta$, and downstream elements stay
anchored). It is returned before the slice migration that step runs, since migration
reads each particle's phase against its slice's window and an assignment made against
the displaced phase left every mover outside its new slice once the nominal phase came
back. Measured
against the analytic slope, and against Genesis4's own
`PHASESHIFTER` element scanned over $\phi = 2\pi\,\delta/(2\gamma^2\lambda)$:
the two bunching-phase curves agree at $1.9\times10^{-8}$ and the exit power against
$\phi$ at $9.3\times10^{-9}$, both codes tracking the same Genesis4-written initial
dumps, so the phase reaches the physics identically, not just the bookkeeping. (On
independently loaded beams the same power comparison sits at $1.4\times10^{-4}$,
which measures the two loaders. Hence the shared start.) The sign is Genesis4's
"delay goes backwards", $\theta \to \theta - \phi$, applied persistently.
Note: Genesis4's phase shifter must have finite length to register at all. A
zero-length one silently does nothing.

In the absolute mode (`bmad_com[absolute_time_tracking] = T` in the
lattice, honored per element through Bmad's own resolver, never a namelist knob)
the walk keeps the real carrier phase of every break:
$\phi_0 \mathrel{-}= 2\pi L/(2\gamma^2\lambda)$ per interlude element, no floors
(whole turns wrap). Inter-segment phasing then follows the real spacings, the
recirculating-linac discipline: a wrong-length break visibly detunes the next
segment, and moving an element moves its phase. The absolute-mode gap scan
reproduces the relative-mode `z_offset` scan point by point (cross-mode
identity).

Chicanes: a break whose elements bend the reference (sbends, patches) detours the
beam while the radiation goes straight, so the light's path is the chord between
the flanking undulator faces (computed from `ele%floor`, never entered by
hand) while `vec(5)` is measured against the reference arc. The
arc-minus-chord delay is charged as whole-wavelength window rotations on the
break's last element (Genesis4's chicane semantics: "always autophasing") plus,
in absolute mode, its carrier phase. The light drift is shortened by the same
delay. Relative mode drops the geometric fraction (measured flat under an angle
scan). Absolute mode ramps at exactly the slope an independent 2D geometric
computation of $d(\mathrm{arc}-\mathrm{chord})/d\alpha$ predicts. Only a closed
bump keeps the light on the next undulator's axis: floor orientations and the
chord direction are asserted, anything else refused, as are bends and
patches under the genesis-model interludes (Genesis4's drift/quad set cannot
represent them) and a `z_offset` that exceeds its upstream break. The
unaveraged energy ledger closes across a chicane sandwich. The window rotations the
delay buys are checked time-dependently (steady state never exercises them): a
chicane banks exactly $\lfloor \mathrm{delay}/\lambda \rfloor$ more escaped slices
than its arc-matched straight twin, at three delays, and the run is byte-identical
at 1 and 8 threads. Delays for that check are chosen mid-interval: a target of
exactly $3.000\lambda$ lands on the floor boundary, where the $6\times10^{-8}$
difference between the exact trace and the small-angle form decides between 2 and 3.
The boundary was located at $2.99 \to 2$, $3.01 \to 3$, which independently
confirms the walk's floor-derived delay against the trace at that scale.

(sec-stats)=
## Diagnostic output

This section is what the tracker computes at each record and how each quantity is
defined. Where those quantities land in the file is not: the layout is
`bmad-stats` 1.0 with the `fel` extension, defined normatively by
`doc/BMAD-STATS-SPEC.md` with `doc/BMAD-STATS-EXT-FEL.md`, and described for
a reader by `doc/reading-output.md`. The file is self-describing, so nothing here
restates it.

*Beam statistics.* Per record and per slice the tracker accumulates the
sufficient statistics of `bunch_params_struct`, named as that struct names
them: `centroid`, the `sigma` covariance matrix at its natural rank,
`charge_live`, `n_particle_live`, `t` and `sigma_t`. The
moments are computed two-pass, for the reason [](#sec-numerics) gives. Beside them
are the per-coordinate extremes `rel_max` and `rel_min`, taken relative to
the centroid in Bmad's own convention, so the envelope a plot wants is the centroid
plus the extreme. Those are order statistics, which no moment can reconstruct, and
they are accumulated in the same sweep. `current` and `energy` in eV follow,
since every consumer would otherwise re-derive them.

*Twiss.* `beam/slice_twiss/` and `beam/bunch/` hold the fully
evaluated struct from Bmad's own `calc_bunch_params`, per slice and for the
whole window. The projected planes follow from the sigma matrix exactly. The normal
modes follow from an eigendecomposition, Wolski's eigen-emittances. The two are kept
apart because they are two decompositions of one beam: an eigen-emittance is not a
projected emittance, and carrying both on one axis would invite an average across
them. Every twiss quantity is a pure function of the stored moments, and the file says
so, so a reader re-deriving them knows what is independent information. Any (record,
slice) reconstructs a `bunch_params_struct`
(`scripts/bunch_params_from_stats.py`, transcribing Bmad's
`projected_twiss_calc` and Wolski's eigen-emittances), and the harness holds
that reconstruction against the stored twiss.

*Slice positions.* The slice axis carries three positions. `ct_slice` and
`t_slice` are exact and free of $\beta$. `z_slice` is Bmad's $z$ at the
reference $\beta$ and says so (Eq. [](#eq-tslice) and [](#sec-window)). All three
publish `@head_direction`, the migration invariant ([](#sec-migration)) that
the high slice index is the window head. Without it a per-slice profile cannot be
drawn, let alone compared with another code's. Stating it on one of a set would invite
the reader to assume the set agrees.

*One sweep for two instruments.* The per-record stats loop also evaluates the
`diag.txt` instrument, the identical `fel_field_diag` and
`fel_slice_diag` calls per slice, so each slice's arithmetic is unchanged and
the diag writer only prints. Every benchmark tier reproduces `diag.txt` bit for
bit. That removed the serial per-record diag sweeps, which cost more than the
whole statistics machinery does, and the measured saving is recorded in
`doc/reading-output.md`.

*Wavefront params.* The field analog (`wavefront_params_struct`,
`wavefront_mod`) describes one field slice: `centroid(4)` $= (x,
\theta_x, y, \theta_y)$ intensity-weighted, `sigma(4,4)` Wigner second moments,
`energy`, `power`, `on_axis_intensity`, `emit_x/y`
$= \sqrt{\det}$ of a plane's $2\times2$ block ($= M^2\lambda/4\pi$), and
`angle_moments_valid`, which follows the `twiss_valid` pattern: spatial moments
come from $|E|^2$ every record, the $\theta$ rows need spectral sums (one forward, two
inverse FFTs) and fill at element ends and bank time only. Not computed is NaN, never a
zero that reads as an answer: an empty slice has no moments and says so, and a NaN
propagates through whatever a consumer does with it where a sentinel would quietly plot.
Pulse-level values are pooled downstream (energy-weighted mean of slice sigmas plus the
variance of slice centroids, in the analysis scripts). The file stays raw.

Nothing is ever summed across harmonics: a detector separates colors. The polarization
components and the harmonics each carry their own copy of the same quantities, and the
sum over the live polarizations of one wavelength is written whether one polarization
is live or two, so no dataset's meaning depends on what else the file holds.

*Dumps at elements.* `dump_beam_at` / `dump_field_at` name
elements through Bmad's own locator (`class::name` syntax for free). An entry
matching nothing is refused. openPMD, like every other dump.

*The escaped-field bank* (`keep_escaped_field`). Field that slippage has
transmitted beyond the window is fixed information: slippage is one-directional, so it
never re-interacts and evolves by free space alone. Each transmitted slice streams to
`<out_root>-escaped.fld.h5` at the zero fill (peak memory a handful of grid
planes. One banked slice per `n_wavelength`$\cdot\lambda$ of slippage, inheriting the
window's decimation) with its `wavefront_params` and transmission $z$. These two
diagnostic files keep Genesis4 field conventions rather than openPMD, since the
per-slice `wavefront_params` and transmission $z$ beside each plane have no home
in the wavefront extension. They are diagnostics: nothing reads them back as state. At
finalize each is propagated to the exit plane by the `exp(K2\,dz)` kernel and
the full pulse is written (`<out_root>-pulse.fld.h5`, live window + banked, the
earliest-transmitted light furthest ahead). Whole-pulse statistics never propagate
numerically: free-space propagation is an ABCD map on `sigma(4,4)`, so
$\sigma_x^2(z) = \sigma_x^2 + 2z\,\sigma_{x\theta} + z^2\sigma_\theta^2$ per
banked slice, exact.

Measured levels and how they are checked: [](reading-output.md).

(sec-pool)=
### The whole-window row

The element-end row for the whole window is assembled from the per-slice
moments by the pooled-covariance identity,

$$
  m = \frac{\sum_s w_s m_s}{\sum_s w_s}, \qquad
  S = \frac{\sum_s w_s \left[ S_s + (m_s - m)(m_s - m)^{\mathsf T} \right]}{\sum_s w_s},
$$ (eq-pool)

never by concatenating the window's particles into one bunch. The between-group term
is what makes [](#eq-pool) exact rather than an average of covariances. Each slice
enters moved from its local $z$ chart to the global window chart. Because
`fel_concat_slices` places a particle at $z_{\rm global} = z_{\rm local} +
\beta (i_s - 1)\Delta$ with the particle's own $\beta$, that move is a shear of
$(z, p_z)$ once $\beta$ is linearized about the slice's mean $p_z$, and the pool
applies it as $S \to J S J^{\mathsf T}$ with the single off-diagonal
$J_{56} = (i_s - 1)\Delta \, {\rm d}\beta/{\rm d}p_z$. The residual is second order
in that expansion.

Measured levels and how they are checked: [](reading-output.md).

(sec-meta)=
### Provenance

`stats.h5` carries a `meta/` group: the resolved input echo, the
top-level lattice file's base name and text, how many files the parser opened, an ISO
timestamp and the Bmad version. What it holds and how it is stored are in
`doc/reading-output.md`. Two things about it are design rather than layout.

It is not a reproducibility record, and it does not claim to be.
`meta/lattice_source` is the top-level file only: Bmad's `call, file =`
pulls in more, and every wrapper lattice in this tree records a call statement while the
lattice it calls is absent. What the file offers for reproduction is the
`lattice/` table, every tracked element with the values the physics used, beside
the input echo. Serialization is not an option:
`write_bmad_lattice_file` inlines a `grid_field` as ASCII under
`one_file$` and writes sibling binary files otherwise, so no `output_form`
is both complete and bounded, and Bmad has no HDF5 lattice format to borrow.
`lat%creation_hash` is no substitute either, hashing inode and size, which makes
it machine specific rather than a fingerprint of content. The complete record wants an
upstream change: `bmad_parser` builds the full list of files it opened, uses it
for that hash, and deallocates it before returning.

A stats file is meant to travel, attached to a paper, mailed to a collaborator, posted
beside a figure, so nothing in it identifies a person by default. The timestamp and the
Bmad version identify the run. The user name and working directory go in only under
`global%record_environment`, and the lattice file is recorded as a base name.
What a user types into the namelist is echoed as typed, so an absolute path in the input
is still an absolute path in the file: relative paths are the user's half of that
contract. Genesis4 records user and cwd unconditionally. Parity is not a reason to leak.

(sec-numerics)=
## Numerical practices

*Two-pass variance.* All second moments subtract the mean before squaring: for
$\gamma \sim 10^4$ with $\sigma \sim 10^{-3}$ the one-pass
$\langle v^2\rangle - \langle v\rangle^2$ form loses the entire $\sigma$ to
cancellation, which stayed hidden until a uniform wake kick first moved the mean.

*Constants floors.* Comparisons against Genesis4 bottom out at the mismatch of
its truncated constants, per term: $Z_0 = 376.73$ ($8.3\times10^{-7}$, the FEL core
and resistive/geometric wakes), $\varepsilon_0 = 8.85\times10^{-12}$ in the
long-range space charge ($4.7\times10^{-4}$), $e = 1.6\times10^{-19}$ and
$\varepsilon_0 = 8.854\times10^{-12}$ in the roughness coefficient
($1.4\times10^{-3}$ of that kernel). Each tier's tolerance is sized to the floor of
the terms it enables.

*Step size.* Measured convergence (one shared 32-slice dump at
$\texttt{n\_wavelength}=12$, `ds_step` swept): roughly first order in the exponential
gain region. Saturation power holds to $\sim$3% up to six periods per step and
misses by 26% at twelve. Two to three periods per step is the operating point.
SIMPLEX's twelve-period economy belongs to its semianalytic solver and does not
transfer.

*Ulp accounting.* The $\gamma \leftrightarrow p_z$ round trip costs $\sim$1 ulp
of $\gamma$ per step (it moved tier1 from $2.5\times10^{-12}$ to
$2.8\times10^{-11}$ at transcription level). A lattice-attribute round trip
($a_w$ through `b_max`) costs 1 ulp amplified to $3.4\times10^{-12}$ over the
line. Moving verbatim code between compilation units shifts arithmetic 1--2 ulp
(gain-amplified to $1.2\times10^{-4}$ of bunching phase). Bit-identity anchors hold
within one compilation unit. Across moves, anchor at ulp tolerance and record the
drift.

*Thread safety.* Cross-slice state changes only at serial points: the FFTW plan
cache is threadprivate, the $K_2$ kernel and its step propagator $e^{K_2\delta z}$
build serially before parallel loops and are read-only inside them, source
accumulators are per call, and slippage/migration/wake updates run at barriers. Each slice's arithmetic is independent of which thread runs
it, so results are bit-identical across thread counts. That is a checked property, and
the only instrument that catches shared-accumulator races (invisible to every
single-threaded check).

Measured levels and how they are checked: [](validation.md).

(sec-coherent-source)=
## The coherent source (SIMPLEX hybrid)

`global%source_model = 'coherent'` replaces the per-particle source deposit
with Tanaka's coherent retrieval (PRAB **27**, 030703 (2024), implemented from
the paper, since SIMPLEX's own source is unlicensed): the slice bunch factor splits
$b = b_C + b_S$, the spatially incoherent $b_S$ -- whose radiation diffracts away,
and whose per-cell sampling noise is the artifact that makes under-populated runs
overestimate gain -- is dropped, and $b_C$ deposits as one analytic Gaussian carrying
the slice's exact source phasor $S$ (the deposit's own normalization contract,
$\sum \mathrm{crsource} = S$). The width is $\kappa\,\sigma$ with $\kappa$ from the
order-0/1 Laguerre-Gauss sums (the paper's Eq. 26), charge-weighted across the window
so near-empty edge slices cannot poison the fit. The centering and the full
$2\times2$ second-moment matrix (tilt included) are this port's extension, letting
offset and mismatched (but transversely Gaussian) beams through. Everything else
is untouched: the grid, the FFT diffraction, slippage, banking, the stats file, and
the per-particle gather.

Measured on the seeded 6\,GeV lethargy configuration (the regime where the artifact
bites hardest): the plain deposit at $M=128$ per slice fakes $\ln P$ by $+0.42$
where the $M=8192$ reference absorbs. The coherent source at the same $M=128$ stays
within $0.05$, a demonstrated $64\times$ particle reduction at the model's own bias,
which is $1.9\times10^{-2}$ in $|\ln P|$ at large $M$. Guards, each with its own message: a
per-slice excess-kurtosis test against its sampling significance
($\sqrt{24/m_\mathrm{ind}}$), charge-weighted so sparse edge slices are exempt. Refusals
for the unaveraged mode, harmonics (even harmonics are invalid in the method itself),
two polarizations, and -- measured, a $\sim$175$\times$ startup deficit -- dark
starts: spontaneous, spatially incoherent emission dominates SASE startup and is
exactly what the coherent model drops, so seeded runs only. The related SIMPLEX
coarse-stepping observation was priced separately: 12 undulator periods per step
costs $2.6\times10^{-3}$ in $|\ln P|$ on the same configuration (taper and harmonics
not yet tested at coarse steps).
