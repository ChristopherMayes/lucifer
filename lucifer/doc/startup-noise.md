---
title: "Grid and particle-count dependence of SASE power"
short_title: SASE convergence
---

Self-amplified spontaneous emission (SASE) is simulated by starting the radiation field at zero and letting the shot noise of the electron beam seed the instability. Such a simulation is commonly called unseeded. The field power an unseeded simulation reports depends on two numerical parameters that have no physical counterpart: the transverse cell size of the field grid and the number of macroparticles per slice. The dependence is large. For the beam and undulator line described below, with a 29 nm bunch and 1024 macroparticles per slice, the exit power is 0.25 GW with 6.35 µm cells and 37 GW with 0.78 µm cells. At either cell size, four times the macroparticles reduces the exit power by about a factor of six.

The problem matters because the power of a SASE simulation is usually the quantity of interest, and a power that moves by two orders of magnitude with the discretization cannot be quoted without a convergence statement. It also matters because the standard comparison between codes cannot detect it: two codes that discretize the transverse plane and the particle distribution in the same way agree with each other at every grid and every particle count, whatever the physical answer is.

This page reports measurements that separate the physical power from the numerical part. The power in the radiation mode of the FEL is the same for every grid and every particle count within the statistical fluctuation of SASE, it starts at the physical spontaneous emission of the beam, and it saturates where theory predicts. The dependence lies entirely in emission at angles outside the mode, which is the radiation of the point-like macroparticle beamlets into the angles the grid can represent. That emission scales inversely with the number of beamlets and increases as the cells shrink. It is what the total field power of an unseeded simulation measures at the grids and particle counts in common use, including those of the examples in this documentation. The page then gives a criterion for the grid, the time window and the particle count, and names the two changes to the field solver and the particle deposition that would remove the artifact.

The analysis was carried out with the script [startup_noise.py](../tests/scripts/startup_noise.py), which runs the simulations, draws the figures and writes every number quoted here to [startup-noise.json](generated/startup-noise/startup-noise.json). Both files are part of this documentation:

```
python3 startup_noise.py --exe <lucifer> --genesis <genesis4> --pyrepo <openPMD-beamphysics> \
        --examples <lucifer/examples> --latdir <lucifer/tests> --workdir <dir> --out <dir>
```

## The model

The particle loader groups the macroparticles of a slice into beamlets of eight that share their transverse coordinates and energy and differ only in phase, with the phases spaced uniformly over one radiation period. This is the quiet start: a beamlet has no bunching at any harmonic until the loader imposes it. The loader then perturbs the phases so that the rms bunching of each slice at each harmonic equals the shot-noise value $1/\sqrt{N_\lambda}$ for the $N_\lambda$ electrons the slice represents, following Fawley's algorithm ([](references.md#ref-fawley), [](fel-physics.md#sec-noise)). For the beam of [](#tab-sn-beam) a slice of twelve wavelengths represents $7.5\times10^{4}$ electrons and the rms bunching is $3.7\times10^{-3}$ ([](#tab-sn-derived)).

Each beamlet occupies a single transverse point on the deposition grid. The paraxial field solver carries every transverse wavenumber up to the Nyquist wavenumber of the grid, and the source term has no angular dependence, so a point beamlet radiates its noise at the resonant wavelength into every angle the grid can represent. A real undulator of $N_w$ periods emits at the resonant wavelength only within the central cone of half angle $\theta_{con} = \sqrt{1+K^2}/(\gamma\sqrt{N_w})$. Emission at an angle $\theta$ is red-shifted by the relative amount $\gamma^2\theta^2/(1+K^2)$, so emission outside the cone is outside the bandwidth $1/N_w$ of the undulator and, at the angles a fine grid represents, far outside the bandwidth $\rho$ of the FEL. The model applies no such cutoff. Spontaneous emission enters the field in no other way. The switch `radiation_fluctuations` acts on the particle energies and never on the field.

Once the beam is bunched, each beamlet radiates coherently with its own bunching factor, again into every angle the grid can represent. A real beam radiates coherently only within the mode of the FEL, because the bunching is correlated across the beam only there. The wide-angle emission of $N_b$ beamlets of charge $Q_b$ and bunching $b$ is $N_b (Q_b b)^2$ times a geometric factor set by the grid, while the coherent emission is $(N_b Q_b b)^2$ times the coupling to the mode, so the ratio of the two is of order $1/N_b$ times a factor that depends on the cell size. This part of the power is an artifact of the macroparticle representation. The theory used for comparison is Ming Xie's fitting formula for the gain length and the effective shot-noise power and saturation estimates of Saldin, Schneidmiller and Yurkov ([](references.md#ref-sase-theory)).

## Method

The beam and the undulator line are those of the Aramis benchmark distributed with Genesis4, which the examples in this documentation also use ([](#tab-sn-beam), [](#tab-sn-lattice)). The line is six FODO cells, each holding two helical undulator segments of 266 periods, two quadrupoles and four drifts, for twelve segments and 57 m in all. The simulation settings are those of [](#tab-sn-simulation), and the quantities derived from them that the page uses are in [](#tab-sn-derived).

```{table} Electron beam parameters. The beam is matched to the periodic focusing of the cell.
:name: tab-sn-beam

| quantity | symbol | value | unit |
|---|---|---|---|
| beam energy | $E$ | 5.80 | GeV |
| Lorentz factor | $\gamma$ | 11358 | |
| peak current, uniform along the window | $I$ | 3.0 | kA |
| normalized emittance, both planes | $\varepsilon_n$ | 0.40 | µm |
| rms relative energy spread | $\sigma_\gamma/\gamma$ | $8.8\times10^{-5}$ | |
| rms energy spread | $\sigma_\gamma$ | 1.0 | |
| matched beta function, horizontal, at the cell entrance | $\beta_x$ | 8.54 | m |
| matched beta function, vertical, at the cell entrance | $\beta_y$ | 17.4 | m |
| rms beam size at the mean beta function | $\sigma_x$ | 21 | µm |
```

```{table} Undulator line. Each FODO cell holds two undulator segments with a drift, a quadrupole and a drift after each.
:name: tab-sn-lattice

| quantity | symbol | value | unit |
|---|---|---|---|
| undulator type | | helical | |
| undulator period | $\lambda_u$ | 15.0 | mm |
| rms undulator parameter | $a_w$ | 0.849 | |
| segment length | | 3.99 | m |
| periods per segment | $N_w$ | 266 | |
| number of segments | | 12 | |
| drift after each segment | | 0.44 | m |
| quadrupole length | | 0.08 | m |
| quadrupole strength, alternating in sign | $k_1$ | 2.0 | m$^{-2}$ |
| drift after each quadrupole | | 0.24 | m |
| cell length | | 9.5 | m |
| line length | | 57 | m |
| resonant wavelength | $\lambda$ | 1.00 | A |
```

```{table} Simulation settings.
:name: tab-sn-simulation

| quantity | value | unit |
|---|---|---|
| integration step | 45 | mm |
| periods per step | 3 | |
| field grid half width | 200 | µm |
| field grid points per side | 64, 128, 256, 512 | |
| cell size | 6.35, 3.15, 1.57, 0.78 | µm |
| macroparticles per slice | 1024, 4096, 16384, 65536 | |
| beamlet size | 8 | |
| slice spacing, long window | 12 | wavelengths |
| slices, long window | 300 | |
| slices, doubled line | 600 | |
```

```{table} Derived quantities.
:name: tab-sn-derived

| quantity | symbol | value | unit |
|---|---|---|---|
| Pierce parameter | $\rho$ | $5.1\times10^{-4}$ | |
| cooperation length | $\lambda/(4\pi\rho)$ | 16 | nm |
| slippage over the line | | 320 | nm |
| central cone of one segment | $\theta_{con}$ | 7.1 | µrad |
| angle at which the resonant wavelength shifts by $\rho$ | | 2.6 | µrad |
| diffraction angle of the mode | $\lambda/(2\pi\sigma_x)$ | 0.74 | µrad |
| electrons per slice of twelve wavelengths | $N_\lambda$ | $7.5\times10^{4}$ | |
| rms shot-noise bunching of such a slice | $1/\sqrt{N_\lambda}$ | $3.7\times10^{-3}$ | |
```

The time window holds 300 slices at a spacing of twelve wavelengths, 360 nm in all. It is longer than the slippage accumulated over the line, so the interior of the window behaves as a long bunch. The sixty slices nearest the tail lag by the cooperation length, and the thirty slices nearest the head collect the field that slips forward out of the interior ([](#fig-sn-profile)). Every power quoted on this page is therefore the mean over slices 80 to 230 of the per-slice power.

```{figure} generated/startup-noise/window-profile.png
:name: fig-sn-profile

Power per slice against the slice index, from the tail to the head of the window, at the ends of undulators 1, 2, 4, 8 and 12, for 1.57 µm cells and 4096 macroparticles per slice. The shaded band marks the interior slices 80 to 230 over which the powers on this page are averaged.
```

The field is written to file at the ends of undulators 1, 2, 4, 8 and 12, at z of 4.0, 8.7, 18.2, 37.2 and 56.2 m. The far field of each record is its two-dimensional Fourier transform, and the power within an angular radius of 3 µrad is separated from the power outside it. The cut is four diffraction angles of the mode and close to the angle at which the resonant wavelength shifts by $\rho$, so it contains the radiation that can take part in the interaction. At saturation the split is the same for cuts of 2 and 5 µrad. Early in the line the far field is flat in angle, so a cut of radius $\theta$ contains a share of the wide-angle emission proportional to $\theta^2$, and the coherent part becomes identifiable once it rises above that share, from about z = 10 m.

The simulations use the GPU backend, which computes in single precision. One simulation on the CPU in double precision, at 6.35 µm and 1024 particles, differs from the corresponding GPU simulation by $1.5\times10^{-3}$ in the interior power and by $1.9\times10^{-3}$ in the power within the cut, which is the level at which the two backends agree in general ([](validation.md#val-device)).

The script computes the theoretical quantities from the parameters of [](#tab-sn-beam) and [](#tab-sn-lattice). The one-dimensional gain length is 1.34 m and Ming Xie's three-dimensional power gain length is 1.79 m of undulator. The effective shot-noise power is 1.9 kW per slice, the number of electrons per coherence volume is $N_c = 1.5\times10^{6}$, the saturation length is 30.5 m, and the saturated power is 0.5 GW by the estimate of Saldin, Schneidmiller and Yurkov and 8.1 GW by Ming Xie's. A seeded steady-state simulation on the same line, with a 5 kW seed and 8192 macroparticles, has a power e-folding length of 2.25 m along the line, which is 1.9 m of undulator once the drift spaces are excluded, within 6 percent of the fit.

## Emission outside the mode

```{figure} generated/startup-noise/power-outside-and-inside-the-mode.png
:name: fig-sn-split

Power per slice outside (left) and inside (right) an angular radius of 3 µrad against z. Top: 1024 macroparticles per slice at four cell sizes. Bottom: 1.57 µm cells at four particle counts. The dashed line is the effective shot-noise power amplified as $\exp(z/L_g)/9$ with Ming Xie's gain length, capped at his saturated power.
```

At the end of the first segment the power outside the cut is independent of the particle count and increases as the cells shrink, between $1/dx$ and $1/dx^2$ ([](#tab-sn-floor), [](#fig-sn-split), [](#fig-sn-scaling)). The number of transverse modes the grid carries scales as $1/dx^2$, and the linear interpolation of the deposition suppresses the modes nearest the Nyquist wavenumber, which is why the measured exponent falls short of two.

```{table} Power per slice outside 3 µrad at the end of the first undulator, z = 3.99 m, in MW.
:name: tab-sn-floor

| cell size (µm) | 1024 particles | 4096 | 16384 | 65536 |
|---|---|---|---|---|
| 6.35 | 0.37 | 0.36 | 0.36 | 0.29 |
| 3.15 | 1.7 | 1.7 | 1.7 | 1.3 |
| 1.57 | 4.7 | 4.5 | 4.3 | 3.5 |
| 0.78 | 9.9 | 9.2 | 8.7 | 7.1 |
```

At the same point the power within the cut is 0.09 to 0.13 MW per slice for every grid and every particle count. The spontaneous power of the beam within the central cone, eq. (1) of Saldin, Schneidmiller and Yurkov, is 0.71 MW and independent of the number of periods, of which the share within 3 µrad of the 7.1 µrad cone of one segment is 0.13 MW. The loaded noise therefore radiates the physical power at the angles where a real beam radiates. What the model adds is the emission outside the cone.

```{figure} generated/startup-noise/wide-angle-vs-cell-size-and-particle-count.png
:name: fig-sn-scaling

Left: power per slice outside 3 µrad at z = 3.99 m against the cell size, with $1/dx^2$ (dashed) and $1/dx$ (dotted) for comparison. Right: power per slice inside (solid) and outside (dashed) 3 µrad at z = 18.2 m against the particle count, for four cell sizes.
```

Beyond the first segments the power outside the cut grows with the bunching, and it then depends on the particle count. At 1.57 µm and z = 56 m it is 4.2 GW, 0.96 GW, 0.26 GW and 0.096 GW per slice for 128, 512, 2048 and 8192 beamlets. Varying the beamlet size at a fixed particle count, and the particle count at a fixed beamlet count, separates the two ([](#tab-sn-beamlets), [](#fig-sn-beamlets)).

```{table} Power per slice outside 3 µrad with the beamlet size and the particle count varied separately, 1.57 µm cells, in MW.
:name: tab-sn-beamlets

| particles per slice | beamlet size | beamlets | z = 37 m | z = 56 m |
|---|---|---|---|---|
| 4096 | 4 | 1024 | 200 | 520 |
| 4096 | 8 | 512 | 410 | 960 |
| 4096 | 16 | 256 | 1000 | 1900 |
| 2048 | 4 | 512 | 420 | 1100 |
| 8192 | 16 | 512 | 450 | 890 |
```

```{figure} generated/startup-noise/beamlets-vs-particles.png
:name: fig-sn-beamlets

Power per slice inside (top) and outside (bottom) 3 µrad against z, at 1.57 µm cells. Left: 4096 macroparticles per slice with beamlet sizes of 4, 8 and 16. Right: 512 beamlets with 2048 and 8192 macroparticles per slice.
```

The power outside the cut is inversely proportional to the beamlet count and independent of the particle count at a fixed beamlet count. This is the coherent emission of each beamlet with itself.

## Emission inside the mode

The power within the cut follows the same curve for every grid and every particle count ([](#fig-sn-split), right), within the statistical fluctuation of the interior mean, which is about 25 percent for 150 slices and a coherence length of about 10 slices. At z = 37.2 m it is 1.8 to 3.9 GW per slice in fifteen of the sixteen simulations. The saturation point of the doubled line falls between 28 and 33 m in every simulation, against the estimate of 30.5 m. The saturated power lies between the two estimates of 0.5 and 8.1 GW. The interior bunching factor at 37 m is 0.20 to 0.30 in the same fifteen simulations.

Two departures from convergence are measured. In the exponential regime the power within the cut rises as the beamlet count falls: 3.5, 5.3 and 9.8 MW per slice at z = 18.2 m for 1024, 512 and 256 beamlets at 4096 particles, and 11 MW against 6.9 MW at 1.57 µm for 128 against 8192 beamlets. This is consistent with the field of a beamlet acting back on the beamlet itself, and the effect is a factor of about two in the exponential regime at 128 beamlets per slice, the count the examples use. The sixteenth simulation, with 0.78 µm cells and 1024 particles, is the one in which the wide-angle emission reaches 6.8 GW per slice at 37 m, above the coherent saturated power. That emission is drawn from the beam: the bunching factor there is 0.08 and the power within the cut 0.49 GW, a fifth of the other simulations.

## Short and long bunches

```{figure} generated/startup-noise/line-and-doubled.png
:name: fig-sn-lines

Left: total window power against z for a window of 96 slices at three wavelengths, at two cell sizes and two particle counts. Right: power per interior slice against z for a window of 600 slices at twelve wavelengths through a line of twice the length, with the saturated power estimates of Ming Xie (dashed) and of Saldin, Schneidmiller and Yurkov (dotted), and the latter's saturation length (vertical).
```

The SASE examples in this documentation use a window of 96 slices at a spacing of three wavelengths, a bunch 29 nm long, against a cooperation length of 16 nm and a slippage over the line of 320 nm. The field leaves such a bunch within a segment and a half, so the window power is the emission of the bunched beam into the field passing through it, and most of that emission is the wide-angle part ([](#tab-sn-window), [](#fig-sn-lines), left).

```{table} Total window power at the exit of the line for a window of 96 slices at three wavelengths, in GW.
:name: tab-sn-window

| cell size (µm) | 1024 particles | 65536 particles |
|---|---|---|
| 3.15 | 1.18 | 0.047 |
| 1.57 | 7.3 | 0.17 |
```

The two simulations with 65536 particles still differ by a factor of 3.6, which is the ratio of the wide-angle emission between the two cell sizes. The window power of this bunch is therefore the wide-angle emission at every particle count measured. The powers recorded for the SASE examples in this documentation, 3.0 GW at 255 grid points and 2048 particles and 3.13 GW at 151 points and 1024 particles, are of this kind. The ratios those examples measure, slice migration enabled against disabled and a wakefield enabled against disabled, are taken at one grid and one particle count and remain valid as ratios.

A window of 600 slices through a line of twice the length gives the long-bunch result ([](#fig-sn-lines), right). With 65536 particles the interior power at 3.15 and 1.57 µm agrees: a peak of 2.24 and 2.28 GW per slice at z = 37 m, falling to 0.24 and 0.31 GW at 114 m as the saturated field loses coherence. With 1024 particles the same simulations saturate at the same position but hold 1.3 and 2.8 GW at 114 m, which is the accumulated wide-angle emission. Saturation is insensitive to the startup level and to the saturation position. It does not remove the particle-count dependence of the total power reported after saturation.

## Comparison with Genesis4

Genesis4 discretizes the transverse plane and the particle distribution in the same way as Lucifer, and the FEL physics of Lucifer's averaged tracking is transcribed from it ([](validation.md)). The two codes therefore share the dependence described here, and their agreement cannot detect it. The comparison uses an unseeded simulation of the same line with a window of 32 slices spanning 9.6 nm and 2048 macroparticles per slice, at 151, 255 and 511 grid points, with both of Genesis4's field solvers ([](#tab-sn-genesis)). Genesis4 4.6.15 generates the noisy beam and writes its initial particle and field distributions, and Lucifer is started from those same distributions. The two codes agree on the exit power to $1.9\times10^{-6}$ at all three grids, with a worst relative difference over all records of $2.1\times10^{-6}$, while the power itself changes by a factor of 36 across the three grids. That bunch is 0.6 cooperation lengths long, so its power is wide-angle emission by the analysis above. Genesis4's alternating-direction implicit solver gives less wide-angle power than its FFT solver at the finest grid, because its finite-difference form of the transverse Laplacian damps the wavenumbers nearest the Nyquist wavenumber, which the FFT solver propagates exactly.

```{table} Exit window power of Genesis4 4.6.15 for an unseeded simulation of a 32-slice window on the line of the table above, with its FFT solver and its alternating-direction implicit (ADI) solver, in MW, and the relative difference of Lucifer's exit power from the FFT result when started from the same distributions.
:name: tab-sn-genesis

| grid points | cell size (µm) | FFT solver | ADI solver | Lucifer against FFT |
|---|---|---|---|---|
| 151 | 2.67 | 22.6 | 22.6 | $1.9\times10^{-6}$ |
| 255 | 1.57 | 117 | 120 | $1.9\times10^{-6}$ |
| 511 | 0.78 | 816 | 548 | $1.9\times10^{-6}$ |
```

## Recommendations

**Distinguish the mode power from the total.** Write the field to file at the position of interest, take its two-dimensional Fourier transform, and compare the power within an angular radius of about four diffraction angles of the mode, $4\lambda/(2\pi\sigma_x)$, with the power outside it. A converged simulation has the outside part below about a tenth of the inside part at saturation. Without a field dump, repeat the simulation with twice the macroparticles at the same beamlet size. The mode power agrees within the SASE fluctuation, about 25 percent for a window of 150 interior slices, while wide-angle emission halves.

**Estimate the wide-angle share in advance.** At saturation the ratio of the wide-angle emission to the mode power measured here is
$$\frac{P_{wide}}{P_{mode}} \approx \frac{100}{N_b}\left(\frac{1.57\,\mathrm{µm}}{dx}\right)^2,$$
within a factor of 1.5 over cell sizes of 0.78 to 6.35 µm and beamlet counts of 128 to 2048, where $N_b$ is the number of beamlets per slice and $dx$ the cell size. The prefactor is for the beam of [](#tab-sn-beam) and scales with the bunching factor squared. Keeping the ratio below a tenth requires $N_b > 1000\,(1.57\,\mathrm{µm}/dx)^2$, which is 4096 macroparticles per slice at a beamlet size of 8 for 3.15 µm cells and 16384 for 1.57 µm cells.

**Choose the cell size from the beam, not from the field.** Cells of about $\sigma_x/7$, 3 µm here, resolve the beam and the mode. Finer cells increase the wide-angle emission as $1/dx^2$ and leave the mode power unchanged.

**Choose the window from the cooperation length.** A long-bunch result requires a window several cooperation lengths long, 16 nm here, and the slices within a few cooperation lengths of the tail are not representative. A shorter window is a short-bunch problem, and its total power was dominated by wide-angle emission at every particle count measured here.

**Check the bunching factor.** In the converged simulations the interior bunching factor at saturation is 0.20 to 0.30. A markedly lower value with a high total power indicates that wide-angle emission is draining the beam.

## Possible code improvements

Two changes to the code would remove the artifact rather than manage it, and neither has been measured. A transverse deposition kernel wider than one cell would suppress the wavenumbers a point beamlet radiates into. An angular filter on the source term at the central cone would remove the emission a real undulator does not produce at the resonant wavelength. Either would change the field the two codes compute from the same distributions, so it would be measured against the coherent power of this page rather than against Genesis4.
