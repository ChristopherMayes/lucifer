# Source and diffraction explorer

Browser prototype for sampled transverse sources and paraxial diffraction.
The calculation deposits one source contribution and propagates its field.
The calculation contains no FEL gain or particle tracking.

## Run

Use Node.js 22.12 or later. Node.js 26.6 was used for validation.
Run these commands from this directory:

```sh
npm ci
npm run dev -- --host 127.0.0.1
```

Vite prints the local URL. Vite selects another port when the default port is occupied.
The browser uses JavaScript and a Web Worker. Python is not required.
Fonts and the numerical libraries are served locally.

```sh
npm test
npm run build
npx playwright install chromium
npm run test:browser
```

The numerical tests check complex FFT normalization, CIC charge conservation,
Gaussian diffraction, filter inversion, and random-phase variance.
The browser tests check rendered pixels, controls, exports, configuration links,
and layouts from 320 to 1920 pixels. The browser tests use a separate server on
port 4175. That port must be free.

`npm run build` writes `dist/`. Serve that directory with an HTTP server.
Opening the HTML through a `file:` URL does not load the worker.
This prototype is independent of the Fortran build and benchmark harness.

## Source and field

The transverse sample count is $N$. Each sample has equal charge weight $1/N$.
The transverse positions are independent Gaussian draws with rms width $\sigma$
on each axis. The source phase of sample $j$ is $\theta_j$.
Equal phases use $\theta_j=0$. Random phases are independent uniform draws on
$[0,2\pi)$.

Cloud-in-cell deposition distributes each complex weight
$\exp(-i\theta_j)/N$ among four mesh nodes. Let $D_{ab}$ denote the deposited
sum at node $(a,b)$, where $a$ and $b$ are transverse grid indices.
The mesh spacing is $\Delta x=W/M$. The window width is $W$, and the mesh
has $M$ nodes on each axis. Node $a$ lies at $(a-M/2)\Delta x$.
The mesh sizes are 128, 256, 512, and 1024.

The illustrative field scale is $E_*=1\,\mathrm{V/m}$:

$$
E_{ab}(0)=iE_*\frac{2\pi\sigma^2}{\Delta x^2}D_{ab}.
$$

The smooth source uses the sampled Gaussian amplitude instead of particles.
The smooth source is normalized to the same deposited charge sum.
Its field approaches $iE_*\exp[-(x^2+y^2)/(2\sigma^2)]$ for a large window.
The normalization is illustrative. The normalization does not predict FEL power.

The physical-electron option uses
$N_e=\operatorname{round}[I\lambda/(ec)]$ samples for one longitudinal wavelength.
Here $I$ is beam current, $\lambda$ is radiation wavelength, $e$ is the positive
elementary charge, and $c$ is the speed of light. At 1 kA and 10 nm,
$N_e=208194$. Calculations above 2,000,000 samples are refused.

The seed fixes positions across mesh and phase changes. Increasing the sample
count preserves the existing prefix of positions. Phase draws use a separate
random stream. Particle markers show at most 4096 positions.

Independent random macroparticle phases give an on-axis mean intensity of $1/N$
relative to the coherent smooth on-axis intensity. Actual electron shot noise
has $1/N_e$ in the same normalization. This macroparticle option does not implement
Fawley loading or a quiet start. A single realization need not equal either mean.

## Drift and diagnostics

The radiation wavenumber is $k=2\pi/\lambda$. The transverse Fourier components
are $k_x$ and $k_y$, and $k_\perp^2=k_x^2+k_y^2$.
Free propagation through distance $s$ applies

$$
\widetilde E(k_x,k_y;s)=\widetilde E(k_x,k_y;0)
\exp[-ik_\perp^2s/(2k)].
$$

The tilde denotes a transverse Fourier transform. FFT.js supplies the complex
one-dimensional transforms. Row and column transforms produce the two-dimensional
transform. Each inverse transform includes its normalization.
Changing the drift distance changes Fourier phase. Angular intensity and total
power remain constant to numerical precision.

The displayed fluence is $\mathcal F=|E|^2\lambda/(2Z_0c)$ in J/m$^2$.
The vacuum impedance is $Z_0$. The pulse duration used here is $\lambda/c$.
Projected fluence integrates $\mathcal F$ over $y$ and has units J/m.
The total power is $P=\sum_{ab}|E_{ab}|^2\Delta x^2/(2Z_0)$.

The two real-space maps can share a color scale. Their shared maximum is the
larger peak of the two fields. Angular intensity uses the unfiltered smooth
on-axis intensity as its fixed reference. The angular spectrum is an annular
mean with logarithmic bins. The angular spectrum excludes the zero-frequency node.
The $1/N$ line is the point-sample noise level before CIC attenuation.

## Source filter

The filter is initially off. When enabled, the default intensity-loss tolerance
is $\varepsilon=0.75$ and the fractional transition width is $w=0.05$.
The filter multiplies the deposited source spectrum once. Changing the drift
distance does not apply the filter again.

The angular radius is $\alpha=k_\perp/k$. The mode angle is
$\alpha_{\rm mode}=4/(\sigma k)$. The Pierce angle is
$\alpha_\rho=\sqrt{2\rho\lambda/\lambda_u}$, where $\lambda_u=0.03\,\mathrm m$
is the undulator period. The fixed planar-undulator amplitude is $a_w=1$.
The relativistic factor is $\gamma=\sqrt{\lambda_u/\lambda}$.
The undulator wavenumber is $k_u=2\pi/\lambda_u$.
The Alfven current used here is $I_A=17045.090\,\mathrm A$.
The planar coupling factor is $J=J_0(1/4)-J_1(1/4)$, where $J_0$ and $J_1$
are Bessel functions of the first kind. The Pierce parameter is

$$
\rho=\left[\frac{I}{I_A}\frac{J^2}{8\gamma^3\sigma^2k_u^2}\right]^{1/3}.
$$

The protected angle is $\alpha_p=\max(\alpha_{\rm mode},\alpha_\rho)$.
The amplitude transmission at angular radius $\alpha$ is

$$
H(\alpha)=\frac{1}{1+\exp[(\alpha/\alpha_e-1)/w]}.
$$

The filter edge $\alpha_e$ follows from $H(\alpha_p)^2=1-\varepsilon$.
The numerically stable inversion uses $t=\sqrt{1-\varepsilon}$:

$$
\alpha_e=\frac{\alpha_p}
 {1-w[\log t+\log(1+t)-\log\varepsilon]}.
$$

The inversion matches `fel_filter_margin` in [fel_track_mod.f90](../code/fel_track_mod.f90).
The default 75% tolerance gives $\alpha_e=\alpha_p$.
At 1% tolerance and $w=0.05$, $\alpha_e/\alpha_p=1.359691661$.
A nonpositive denominator is refused. The on-axis amplitude must exceed 0.99.
The intensity-loss tolerance does not specify FEL gain accuracy.

The protected disk remains fixed during tolerance scans. Power inside this disk
is $P_{\rm in}$. All remaining represented modes, including square-grid corners,
contribute to $P_{\rm out}$. The diagnostics report $P_{\rm out}/P_{\rm in}$ and
$P_{\rm out}/(P_{\rm in}+P_{\rm out})$. Source retention is filtered power divided
by unfiltered power.

## Limits

The angular band uses the prescribed $\sigma$. Lucifer and the teaching notebook
use measured particle beam sizes for their corresponding setup calculations.
The browser uses a separate seeded random generator and an even periodic mesh.
The browser does not reproduce the notebook's odd-grid particle images bit for bit.

The FFT boundary is periodic. The boundary diagnostic measures the power in the
outer 5% strips on each side of the propagated mesh. A small boundary value does
not establish convergence or exclude earlier wraparound. Compare larger windows
at fixed spacing and finer meshes at fixed window.

The angular interval is $\lambda/W$. The filter-resolution diagnostic reports the
10%-to-90% amplitude transition width in these intervals. Fewer than four intervals
produces a warning. This warning is a resolution indicator, not an error bound.

Calculation time measures worker wall-clock elapsed time. Calculation time excludes plot
rendering. Double precision is used for calculation, and single precision is used
for the displayed heatmaps. Configuration links store model parameters in the URL
fragment. Downloads contain parameters, scalar diagnostics, and annular spectra.