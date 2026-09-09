---
title: "Grid and particle-count dependence of SASE power"
short_title: SASE convergence
---

Self-amplified spontaneous emission (SASE) is simulated by starting the radiation field at zero and letting the shot noise of the electron beam seed the instability. Such a simulation is commonly called unseeded. The field power an unseeded simulation reports depends on two numerical parameters that have no physical counterpart: the transverse cell size of the field grid and the number of macroparticles per slice. The dependence is large. For the beam and undulator line described below, with a 29 nm bunch and 1024 macroparticles per slice, the exit power is 0.25 GW with 6.35 µm cells and 37 GW with 0.78 µm cells. At either cell size, four times the macroparticles reduces the exit power by about a factor of six.

The problem matters because the power of a SASE simulation is usually the quantity of interest, and a power that moves by two orders of magnitude with the discretization cannot be quoted without a convergence statement. It also matters because the standard comparison between codes cannot detect it: two codes that discretize the transverse plane and the particle distribution in the same way agree with each other at every grid and every particle count, whatever the physical answer is.

The measurements reported here separate the physical power from the numerical part by resolving the far field of the radiation into the part inside the mode of the FEL and the part outside it. The power inside the mode is the same for every grid and every particle count within the statistical fluctuation of SASE. It starts at the physical spontaneous emission of the beam and saturates where theory predicts. The power outside the mode carries the whole dependence. It is the radiation of the point-like macroparticle beamlets into the angles the grid can represent, it scales inversely with the number of beamlets, and it increases as the cells shrink. At the grids and particle counts in common use, including those of the examples in this documentation, the total field power of an unseeded simulation is mostly this part. Two of these findings are known. That the power inside the mode is what amplifies, and that two codes agreeing with each other can both depend on the macroparticle count, were established by Tanaka in 2024, and the agreement of the power inside the cone with the spontaneous emission is the known behaviour of the loading. The measurements here add the far-field split against z, the separation of the beamlet count from the particle count, and the dependence on the cell size ([](#prior-work)).

The sections below describe the model and the prior work, then the method, present the two parts of the emission in turn, compare with Genesis4, and close with a criterion for choosing the grid, the time window and the particle count, and with the known remedies.

The analysis was carried out with the script [startup_noise.py](../tests/scripts/startup_noise.py), which runs the simulations, draws the figures and writes every number quoted here to [startup-noise.json](generated/startup-noise/startup-noise.json). Both files are part of this documentation:

```
python3 startup_noise.py --exe <lucifer> --genesis <genesis4> --pyrepo <openPMD-beamphysics> \
        --examples <lucifer/examples> --latdir <lucifer/tests> --workdir <dir> --out <dir>
```

## The model

The particle loader groups the macroparticles of a slice into beamlets of eight that share their transverse coordinates and energy and differ only in phase, with the phases spaced uniformly over one radiation period. This is the quiet start: a beamlet has no bunching at any harmonic until the loader imposes it. The loader then perturbs the phases so that the rms bunching of each slice at each harmonic equals the shot-noise value $1/\sqrt{N_\lambda}$ for the $N_\lambda$ electrons the slice represents, following Fawley's algorithm ([](references.md#ref-fawley), [](fel-physics.md#sec-noise)). For the beam of [](#tab-sn-beam) a slice of twelve wavelengths represents $7.5\times10^{4}$ electrons and the rms bunching is $3.7\times10^{-3}$ ([](#tab-sn-derived)).

Each beamlet occupies a single transverse point on the deposition grid. The paraxial field solver carries every transverse wavenumber up to the Nyquist wavenumber of the grid, and the source term has no angular dependence, so a point beamlet radiates its noise at the resonant wavelength into every angle the grid can represent. A real undulator of $N_w$ periods emits at the resonant wavelength only within the central cone of half angle $\theta_{con} = \sqrt{1+K^2}/(\gamma\sqrt{N_w})$. Emission at an angle $\theta$ is red-shifted by the relative amount $\gamma^2\theta^2/(1+K^2)$, so emission outside the cone is outside the bandwidth $1/N_w$ of the undulator and, at the angles a fine grid represents, far outside the bandwidth $\rho$ of the FEL. The model applies no such cutoff. Spontaneous emission enters the field in no other way. The switch `radiation_fluctuations` acts on the particle energies and never on the field.

Once the beam is bunched, each beamlet radiates coherently with its own bunching factor, again into every angle the grid can represent. A real beam radiates coherently only within the mode of the FEL, because the bunching is correlated across the beam only there. The wide-angle emission of $N_b$ beamlets of charge $Q_b$ and bunching $b$ is $N_b (Q_b b)^2$ times a geometric factor set by the grid, while the coherent emission is $(N_b Q_b b)^2$ times the coupling to the mode, so the ratio of the two is of order $1/N_b$ times a factor that depends on the cell size. This part of the power is an artifact of the macroparticle representation. The theory used for comparison is Ming Xie's fitting formula for the gain length and the effective shot-noise power and saturation estimates of Saldin, Schneidmiller and Yurkov ([](references.md#ref-sase-theory)).

(prior-work)=
## Prior work

The loading of shot noise on macroparticles was solved in one dimension by Penman and McNeil ([](references.md#ref-penman-mcneil)), and Fawley gave the multidimensional algorithm on beamlets that Genesis4 and this code use, with the higher harmonics loaded to the same statistics and the spontaneous emission of the loaded beam as its test ([](references.md#ref-fawley)). McNeil, Poole and Robb extended the loading to unaveraged models ([](references.md#ref-mcneil-2003)). None of these addresses the transverse plane beyond the beamlet sharing one position. What a grid code does with the spontaneous emission it loads was examined by Fawley and by Huang and Kim ([](references.md#ref-fawley-2003), [](references.md#ref-huang-kim-2003)). A transverse grid represents an artificially limited number of modes. The on-axis far-field intensity and the power within the central cone agree with the analytic spontaneous emission, and the total spontaneous power depends on the mode content of the code, growing by more than an order of magnitude from an axisymmetric code to a three-dimensional one while the two agree in the exponential regime. The agreement inside the cone reported on this page is therefore the known behaviour of the loading. The disagreement outside the cone is the known dependence on the mode content, here in the direction of too much rather than too little, because the grid carries the resonant wavelength to angles where the real emission is red-shifted.

That a macroparticle radiates as a single charge many times the electron's, and so radiates anomalously strongly, was stated by Litvinenko, who proposed paired clones of opposite charge to separate the induced from the spontaneous part ([](references.md#ref-litvinenko)). Pausch and co-workers give the mechanism exactly for the radiation diagnostics of particle-in-cell codes: a point-like macroparticle radiates coherently with its whole weight, the Fourier transform of its shape function bounds that emission above the inverse cell size, and a form factor over the electrons it represents restores the incoherent part ([](references.md#ref-pausch)). Andriyash, Lehe and Malka note that the spontaneous emission of a simulation with few macroparticles can differ greatly from the physical one ([](references.md#ref-andriyash)). Tanaka established the point for FEL grid codes ([](references.md#ref-tanaka)): SIMPLEX and GENESIS agree with each other while both overestimate the gain at low macroparticle counts, the local bunching factor holds a spatially incoherent component whose radiation diverges by diffraction and takes no part in the amplification, and dropping that component and adding the physical spontaneous emission back analytically converges with far fewer macroparticles. His grid was fixed. The measurements here add the separation of the beamlet count from the particle count, the dependence on the cell size, the scaling law, and the far-field split against z. The diagnosis and the shared dependence of two codes are his.

Remedies of three kinds exist. Genesis4 carries an angular filter on the source term since release 4.6.12, `source_filter`, a sigmoid in normalized spatial frequency applied to the transformed source term, described in its changelog as suppressing strongly diffracting field components from a rough distribution of the source term ([](references.md#ref-genesis4)). MINERVA represents the field in 20 to 30 Gaussian modes and requires at least 8000 particles per slice, more for transverse resolution or harmonics ([](references.md#ref-minerva)). A field of a few tens of smooth modes carries no wide-angle content. Tanaka's coherent retrieval keeps a few Laguerre-Gauss orders of the bunching and is the extreme filter. Hwang and Qiang treat the artificial noise of particle migration across the longitudinal mesh, the corresponding problem in the other direction ([](references.md#ref-hwang-qiang)). No printed guidance on the cell size was found beyond the Genesis4 manual's advice to use an odd number of grid points for a point on axis. Here 255 and 256 points agree.

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

At the same point the power within the cut is 0.09 to 0.13 MW per slice for every grid and every particle count. The spontaneous power of the beam within the central cone, eq. (1) of Saldin, Schneidmiller and Yurkov, is 0.71 MW and independent of the number of periods, of which the share within 3 µrad of the 7.1 µrad cone of one segment is 0.13 MW. The loaded noise therefore radiates the physical power at the angles where a real beam radiates. What the model adds is the emission outside the cone. This comparison holds where the first segment is short in gain lengths, 2.2 of them here. Over four machines the measured ratio tracks that length, from 0.81 at 1.2 gain lengths to 2.6 at 5.8 ([](#tab-sn-startup-four)).

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

The two simulations with 65536 particles still differ by a factor of 3.6, which is the ratio of the wide-angle emission between the two cell sizes. The window power of this bunch is therefore the wide-angle emission at every particle count measured. The powers the SASE examples report with the filter off are of this kind: 3.04 GW for the `sase` example at 256 grid points and 2048 particles, and 1.83 GW for the `migration` example at 128 points and 1024 particles. The ratios those examples measure, slice migration enabled against disabled and a wakefield enabled against disabled, are taken at one grid and one particle count and remain valid as ratios.

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

(sn-recommendations)=
## Recommendations

**Distinguish the mode power from the total.** Compare the power within an angular radius of about four diffraction angles of the mode, $4\lambda/(2\pi\sigma_x)$, with the power outside it. A converged simulation has the outside part below about a tenth of the inside part at saturation. The tracker measures this and reports it: every element end that takes the field angle moments stores the power inside that radius beside the total, and the run closes by printing the two at the record where the mode peaked and at the last record, warning above a tenth and confirming below a hundredth. The same numbers are in the statistics file, in `run/` for those two records and per record and per slice in `field/total/power_inside_angle`. Without that, repeat the simulation with twice the macroparticles at the same beamlet size. The mode power agrees within the SASE fluctuation, about 25 percent for a window of 150 interior slices, while wide-angle emission halves.

**Estimate the wide-angle share in advance, then measure it.** At saturation the ratio of the wide-angle emission to the mode power measured on the first machine is
$$\frac{P_{wide}}{P_{mode}} \approx \frac{100}{N_b}\left(\frac{1.57\,\mu\mathrm{m}}{dx}\right)^2,$$
within a factor of 1.5 over cell sizes of 0.78 to 6.35 µm and beamlet counts of 128 to 2048, where $N_b$ is the number of beamlets per slice and $dx$ the cell size. The prefactor is for the beam of [](#tab-sn-beam) and scales with the bunching factor squared. It does not carry to another machine: on FLASH1 at 13.7 nm the same expression is 40 times low over cell sizes of 3.9 to 31 µm. The form is the guide for choosing a load before a run, and the number that decides the question is the one the run reports. Keeping the ratio below a tenth requires $N_b > 1000\,(1.57\,\mu\mathrm{m}/dx)^2$, which is 4096 macroparticles per slice at a beamlet size of 8 for 3.15 µm cells and 16384 for 1.57 µm cells.

**Filter the source, which is now the default.** On 2026-09-07 `global%source_filter` became true by default, on the measurements of this page: four machines, two of them real, at angle ratios of 1.14 to 1.86. At a converged load it removes about a hundred times the wide-angle power while moving the mode power a few percent and the saturation point not at all, and it makes 1024 macroparticles per slice reach the answer of the unfiltered 4096 at 60 percent of the wall clock. Turn it off to compare against a code that carries no such filter, which is what the comparison tiers do. What it costs, and the four configurations measured side by side, is [](fel-physics.md#sec-convergence).

**Choose the cell size from the beam, not from the field.** Cells of about $\sigma_x/7$, 3 µm here, resolve the beam and the mode. Finer cells increase the wide-angle emission as $1/dx^2$ and leave the mode power unchanged. A half width of about nine beam sizes contains the mode, and above that the half width does not enter. The tracker derives both from the beam when the deck states neither, which fixes the point count at 127 before it is rounded up to 128, since the beam size cancels between the two rules ([](fel-physics.md#sec-field)). On the Aramis benchmark that is 128 points over 192 µm against the 256 over 200 µm these examples state, and it lowers the ratio of wide-angle to mode power at the exit of the `sase` example from 24 to 5.0.

**Choose the window from the cooperation length.** A long-bunch result requires a window several cooperation lengths long, 16 nm here, and the slices within a few cooperation lengths of the tail are not representative. A shorter window is a short-bunch problem, and its total power was dominated by wide-angle emission at every particle count measured here.

**Check the bunching factor.** In the converged simulations the interior bunching factor at saturation is 0.20 to 0.30. A markedly lower value with a high total power indicates that wide-angle emission is draining the beam.

## The source filter, measured

Of the four remedies below, one is now in this code and measured against the split this page defines. `global%source_filter` multiplies the transformed source by a sigmoid in normalized transverse spatial frequency, transcribed from Genesis4's own filter and agreeing with it at 2.2e-6 ([](fel-physics.md#sec-source-filter), [](validation.md#val-source-filter)). Two things set it: where the sigmoid's edge sits, and how sharply it falls. Both matter, and the second decides whether the filter helps or harms.

The measurement is the sweep of this page at 1.57 µm cells, with the filter's edge at the central cone of one segment, 7.1 µrad, and at the 3 µrad cut this page splits at, each at two widths. Genesis4's own default width of 1 is in the same units as the edge, which puts the sigmoid at 0.73 on axis: at that width the filter attenuates the coherent source as much as the wide angles. The sharp rows use a width of 0.05 of the edge. Genesis4's default edge, `xcut = 1`, lies at twice the Nyquist angle, off the grid, so its default filter is the soft roll-off alone.

The derived default rows set nothing. The run computes the edge from the beam and the gain as the larger of four mode diffraction angles, 2.97 µrad here, and the angle at which the resonant wavelength red-shifts by ρ, 2.61 µrad ([](fel-physics.md#sec-source-filter)), at the default width. It lands on the 3 µrad cut this page chose by hand, and the two rows are the same filter to the digits shown.

```{table} The source filter against the split, 1.57 µm cells, 300-slice window. The wide-angle column is the factor by which the power outside 3 µrad falls at the exit. The in-cone column is the power inside 3 µrad at the end of the first undulator, where the physical spontaneous emission within that angle is 0.13 MW. The mode column is the change in the power inside 3 µrad at the exit, against a SASE fluctuation of about 25 percent.
:name: tab-sn-filter

| particles | edge | width | wide-angle at exit | in-cone at z = 4 m (MW) | mode at exit | bunching at 37 m | saturation (m) |
|---|---|---|---|---|---|---|---|
| 1024 | none | | 1 | 0.133 | | | 28.2 |
| 1024 | derived, 2.97 µrad | 0.05 | 1500x | 0.111 | +102% | +28% | 32.9 |
| 1024 | 7.1 µrad | 1 | 12x | 0.058 | +24% | -40% | 37.7 |
| 1024 | 7.1 µrad | 0.05 | 10x | 0.129 | +69% | +17% | 32.9 |
| 1024 | 3 µrad | 1 | 85x | 0.043 | +56% | -53% | 37.7 |
| 1024 | 3 µrad | 0.05 | 1300x | 0.112 | +102% | +28% | 32.9 |
| 4096 | none | | 1 | 0.123 | | | 32.9 |
| 4096 | derived, 2.97 µrad | 0.05 | 1100x | 0.106 | +19% | +5% | 32.9 |
| 4096 | 7.1 µrad | 1 | 13x | 0.055 | -6% | -63% | 37.7 |
| 4096 | 7.1 µrad | 0.05 | 9x | 0.122 | +12% | 0% | 32.9 |
| 4096 | 3 µrad | 1 | 83x | 0.042 | -7% | -70% | 42.4 |
| 4096 | 3 µrad | 0.05 | 860x | 0.107 | +19% | +5% | 32.9 |
```

```{figure} generated/startup-noise/source-filter.png
:name: fig-sn-filter

Power per slice inside 3 µrad (circles) and outside it (crosses) against z, with the filter off and at four settings, for 1024 and 4096 macroparticles per slice at 1.57 µm cells.
```

The wide-angle column is the filter doing what a filter does. Free propagation conserves angle, so once the source is cut beyond the edge nothing arrives there, and a sharp edge on the 3 µrad cut leaves 0.9 MW per slice outside it at the exit where the unfiltered run had 0.96 GW. That column says the filter works. The other three columns say what it costs, and they are the measurement.

At Genesis4's default width the filter halves the power inside the cone at the end of the first undulator, from 0.13 MW to 0.04 or 0.06 MW, where 0.13 MW is the physical spontaneous emission the beam radiates into that angle. It removes half the startup seed along with the artifact, and the gain arrives late: the bunching at 37 m falls by a third to two thirds and saturation moves from 28 and 33 m out to 38 and 42 m. The power at the exit is then no better a measure of the machine than before, since the run has not finished saturating.

At a sharp edge on the 3 µrad cut, the derived default, the in-cone startup power is 0.107 to 0.112 MW, which is 14 percent under the unfiltered 0.123 to 0.133: the sigmoid's fall, 0.15 µrad wide, sits inside the measurement cut and trims the edge of the physical cone. Saturation absorbs it. At 4096 particles the bunching at 37 m is unchanged to 5 percent, the saturation point does not move, and the power inside the mode at the exit rises 19 percent, which is the mode keeping power that used to diffract away and is inside the fluctuation of the process. At 1024 particles the same filter raises the mode power by a factor of two and the bunching at 37 m by 28 percent, and saturation moves from 28.2 to 32.9 m. This page's earlier sections found that at 128 beamlets the wide-angle emission drains the beam. The filter removes the drain, and the beam at 1024 particles then saturates where the beam at 4096 does. A sharp edge at the cone, 7.1 µrad, removes 9 to 10 times the wide-angle power and leaves a tenth of it, since the emission between the cut and the cone is the larger share.

With the default filter on, the outside part is 0.2 percent of the inside part at the exit, and the total power is a converged quantity by the criterion of [](#sn-recommendations). The filter costs 14 percent of the wall clock at 4096 macroparticles, 31 s against 27 s, since the source gains a transform pair it did not need.

Two cautions carry from this. The width is not a detail: Genesis4's default leaves the sigmoid at 0.73 on axis, and a filter meant to remove only the wide angles wants a width well below the edge it cuts at, which is why the code refuses a width that puts the axis below 0.99. And an edge near the mode trims the physical seed by a measured 14 percent that saturation hides. A filter removes what it removes and restores nothing, and the seed a soft filter takes out is exactly what Tanaka's remedy adds back analytically ([](references.md#ref-tanaka)). This measurement is the case for that add-back rather than against it.

## The default on a second machine

A default verified on one machine has been verified against that machine's coincidences. On Aramis the two angles the default is built from, four mode diffraction angles and the ρ angle, are 2.96 and 2.61 µrad, and both sit in the wide band between the mode and the grid's Nyquist angle of 32 µrad, so an edge four times too wide would have looked as good as the right one. It did, for a day: the first version of the conversion from an angle to Genesis4's `xcut` was a factor of four off, and every edge quoted above was measured at four times its stated angle before a second machine exposed it (FINDINGS 7.51).

The second machine is a probe at FLASH's scales, `tests/bmad/flash.bmad`: 10 nm from an 803 MeV beam through twelve helical segments of 165 periods at 27.3 mm, 1.5 kA, 1.5 µm normalized emittance, a beam five times the size of Aramis's and a wavelength a hundred times longer. It is not FLASH, and no number in it comes from a FLASH publication. FLASH1 itself is the third machine below. Its quadrupoles hold a mean beta function of 10 m over the cell, and the Twiss it states is that cell's matched periodic solution. It stated a beta of 10 m with zero alpha until 2026-09-06, which is not that solution, so every number in this section was measured again on the matched beam (FINDINGS 7.57). Its rms beam size is 98.7 µm and its ratio of Rayleigh length to gain length is 2.2 times smaller than Aramis's, so the two angles come out in a different ratio: the mode angle is 64.5 µrad and the ρ angle 34.7, ratio 1.86 against Aramis's 1.14, and the mode angle wins again. The central cone of one segment is 66.6 µrad, so on this machine the default edge and the physical cone nearly coincide. The grid is the beam's, 128 points at 15.7 µm cells, and the load is 1024 macroparticles per slice in a 300-slice window. Its lattice states no integration step, so the step is derived from the gain at one undulator period, 165 steps a segment. It ran at 3300 until then, Bmad's own default for a wiggler being a twentieth of a period, and each run takes 8 to 9 s on the device rather than 67 to 85. It took 35 minutes on twelve CPU threads until the slice spacing became an integer the deck states rather than a real the code recovered by division, which the device's exact bucket arithmetic refused at 12 to the last bit ([](fel-physics.md#sec-window)).

```{table} The derived default on the second machine. The split is at the derived edge, 64.5 µrad, and the bunching is read at the end of the eighth undulator, z = 43 m.
:name: tab-sn-filter-flash

| filter | in-cone at z = 4.5 m (kW) | inside at exit (GW) | outside at exit | bunching at 43 m | saturation (m) |
|---|---|---|---|---|---|
| none | 6.59 | 1.74 | 625 MW | 0.124 | 21.6 |
| derived default, 64.5 µrad, width 0.05 | 4.73 | 2.12 | 3.9 MW | 0.133 | 21.6 |
```

```{figure} generated/startup-noise/source-filter-flash.png
:name: fig-sn-filter-flash

Power per slice inside 64.5 µrad (circles) and outside it (crosses) against z on the second machine, with the filter off and at its derived default, 1024 macroparticles per slice at 15.7 µm cells.
```

The default holds. The power inside the mode at the exit rises 22 percent, saturation does not move, and the power outside the edge falls from 625 to 3.9 MW per slice, 160 times. The default removes physical emission along with the artifact, as it does on Aramis and more of it: the in-cone startup power falls from 6.59 to 4.73 kW, 28 percent, because the edge sits on the physical cone and the sigmoid's fall, 3.2 µrad wide, takes the cone's rim with it. Saturation absorbs it, as it did the 14 percent on Aramis, and the bunching at 43 m is 7 percent higher with the filter than without. On both machines the mode angle won, so the ρ angle's case, where four mode angles fall inside the bandwidth of the interaction, is not yet measured: it needs $z_R/L_g$ above about 50, a hard X-ray line with a small beam. Until it is, the default is measured over ratios of 1.1 to 1.9, and the warning the code prints outside 0.3 to 3 is where that measurement's authority ends.

## The default on a third machine, a real one

The third machine is FLASH1 at DESY, the line that lased at 13.7 nm in 2006, at its published parameters: six planar segments of 165 periods at 27.3 mm and 0.47 T, a 668 MeV beam of 2.5 kA and 1.5 µm normalized emittance through a doublet lattice holding the beta function at 10 m. The lattice, the beam and the sources are [`examples/flash1`](../examples/flash1/README.md). It is planar where the other two machines are helical, so its undulator couples to the resonance through the Bessel factor 0.885 rather than through unity, and that factor enters the Pierce parameter and the spontaneous power alike.

Its rms beam size is 107 µm, so the grid is 128 points over a 1 mm half width at 15.7 µm cells, again the beam's own $\sigma_x/7$. The two angles the default is built from are 81.6 µrad, four mode diffraction angles, and 44.4 µrad, the ρ angle, a ratio of 1.84 against the second machine's 1.86 and the first's 1.14. The mode angle wins for the third time. The run derives its edge from the beam it loaded rather than from the nominal one, which puts it at 81 µrad, again beside the central cone of one segment, 78.0 µrad.

```{table} The derived default on the third machine. The split is at the derived edge, 81 µrad, 1024 macroparticles per slice at 15.7 µm cells, and the bunching is read at the end of the fourth undulator, z = 20.4 m.
:name: tab-sn-filter-flash1

| filter | in-cone at z = 4.5 m (kW) | inside at exit (GW) | outside at exit | bunching at 20.4 m | saturation (m) |
|---|---|---|---|---|---|
| none | 8.93 | 2.07 | 754 MW | 0.134 | 20.7 |
| derived default, 81 µrad, width 0.05 | 6.22 | 2.75 | 6.1 MW | 0.177 | 20.7 |
```

```{figure} generated/startup-noise/source-filter-flash1.png
:name: fig-sn-filter-flash1

Power per slice inside 81 µrad (circles) and outside it (crosses) against z on FLASH1, with the filter off and at its derived default, 1024 macroparticles per slice at 15.7 µm cells.
```

The default holds a third time. The power outside the edge at the exit falls from 754 to 6.1 MW per slice, 123 times, the power inside the mode rises 33 percent, the bunching at 20.4 m rises 10 percent, and the saturation point does not move. What it removes with the artifact is again the rim of the physical cone: the in-cone startup power falls from 8.93 to 6.22 kW, 30 percent, against 28 percent on the second machine and for the same reason, that the sigmoid's fall sits on the cone. The fourth machine puts the range these ratios cover on a footing ([](#tab-sn-filter-four)).

The startup power is the one measurement that does not transfer. On the first machine the power within the cut at the end of the first segment equals the spontaneous emission of the beam within that cut, which is what says the loading is right. On FLASH1 the measured 8.93 kW is 2.6 times the 3.40 kW the same formula gives, and the excess is not the discretization: it is the same 2.5 to 3.4 times at every cell size from 3.9 to 31 µm and at every load from 1024 to 65536 macroparticles per slice. The reason is the gain length. FLASH1's first segment is 4.50 m against a power gain length of 0.78 m, which is 5.8 gain lengths, where the Aramis benchmark's first segment is 2.2 and the fourth machine's is 1.2 ([](#tab-sn-startup-four)). The in-cone power grows by a factor of 70 between the first and second dumps on FLASH1 and by 3.2 on Aramis over the same ratio in z, so by the end of the first FLASH1 segment the radiation has been amplified and is no longer the spontaneous emission alone. A startup measurement needs a dump within the first gain length, not at the first segment's end.

## The default on a fourth machine, and the range it is measured over

The fourth machine is LCLS at SLAC, the line that lased at 1.5 Angstrom in April 2009, at its published parameters: thirty-three planar segments of 114 periods at 3.0 cm and 1.25 T, a 13.6 GeV beam of 3.0 kA and 0.4 mm mrad through a doublet lattice holding the beta function at 30 m over 132 m. The lattice, the beam and the sources are [`examples/lcls`](../examples/lcls/README.md). Its rms beam size is 21.3 µm, within a tenth of the Aramis benchmark's, so it tests the default at the first machine's transverse scale and at a gain length two and a half times longer.

```{table} The derived default on the four machines. The split is at each machine's derived edge. The wide-angle column is the factor by which the power outside that edge falls at the exit, the mode column the change in the power inside it, and the bunching is read at the last dump before saturation.
:name: tab-sn-filter-four

| machine | λ | σₓ (µm) | z_R/L_g | mode angle | ρ angle | ratio | load | wide-angle | mode | bunching | saturation (m) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Aramis | 0.1 nm | 21.4 | 32.1 | 2.98 µrad | 2.62 µrad | 1.14 | 1024 | 1546x | +102% | +28% | 28.2 to 32.9 |
| Aramis | | | | | | | 4096 | 1064x | +19% | +5% | 32.9 to 32.9 |
| FLASH probe | 10 nm | 98.7 | 14.4 | 64.5 µrad | 34.7 µrad | 1.86 | 1024 | 160x | +22% | +7% | 21.6 to 21.6 |
| FLASH1 | 13.7 nm | 106.9 | 13.5 | 81.6 µrad | 44.4 µrad | 1.84 | 1024 | 123x | +33% | +10% | 20.7 to 20.7 |
| LCLS | 0.151 nm | 21.3 | 13.6 | 4.52 µrad | 2.43 µrad | 1.86 | 1024 | 112x | +19% | +17% | 51.8 to 55.7 |
| LCLS | | | | | | | 4096 | 105x | +4% | +3% | 55.8 to 55.8 |
```

The default holds a fourth time, and at the converged load it is close to free: at 4096 macroparticles per slice on LCLS it removes 105 times the wide-angle power while moving the mode power 4 percent, the bunching 3 percent and the saturation point not at all. Where the load is not converged the filter also removes what the wide-angle emission was taking from the beam, which is the +102 percent on Aramis at 1024 and the +19 percent on LCLS at the same count.

**The range the default is measured over.** The two angles are four mode diffraction angles $4\lambda/(2\pi\sigma_x)$ and the ρ angle $\sqrt{2\rho\lambda/\lambda_u}$, and their ratio goes as $1/\sqrt{z_R/L_g}$. Over the four machines that scaling holds to 5 percent, the four coefficients being 6.45, 6.75, 6.88 and 7.05, with

$$
  \frac{4\lambda/(2\pi\sigma_x)}{\sqrt{2\rho\lambda/\lambda_u}} \;\approx\; \frac{6.8}{\sqrt{z_R/L_g}} ,
$$ (eq-anglerule)

so the mode angle stops winning only above $z_R/L_g \approx 46$. None of the four reaches it, and two of them are real machines at their published focusing. A beam that did would be some ten diffraction-limited emittances wide, which a hard X-ray line is built to avoid: LCLS at 1.5 Angstrom sits at 13.6 and the Aramis benchmark at 0.1 Angstrom, the widest of the four, at 32.1. The default is therefore measured over ratios of 1.14 to 1.86, and that is the range machines occupy rather than a limitation of the measurement. The setup warning outside 0.3 to 3 is left as it is, since it bounds a wider band than any of the four.

**The in-cone startup, and what decides whether it can be measured.** [](#emission-inside-the-mode) reports that the power inside the cut at the end of the first segment equals the spontaneous emission of the beam inside that cut, which is what says the loading is right, and the third machine read 2.6 times the analytic value instead. The fourth machine settles the reason. The comparison is clean only where the first segment is short in gain lengths, and over the four the measured ratio tracks that length exactly:

```{table} The power inside the reported cut at the end of the first segment, against eq. (1) of Saldin, Schneidmiller and Yurkov for the same cut, at the most converged load of each machine.
:name: tab-sn-startup-four

| machine | first segment in power gain lengths | measured | analytic | ratio |
|---|---|---|---|---|
| LCLS | 1.24 | 0.100 MW | 0.123 MW | 0.81 |
| Aramis | 2.23 | 0.123 MW | 0.128 MW | 0.96 |
| FLASH probe | 5.29 | 6.59 kW | 3.58 kW | 1.84 |
| FLASH1 | 5.81 | 8.93 kW | 3.40 kW | 2.63 |
```

Below about two gain lengths the first dump is still spontaneous emission and the comparison measures the loading. Above about five it measures amplified light, and the excess is the gain, not a defect in the loading. A startup measurement on a machine with a short segment relative to its gain length needs a dump inside the first gain length rather than at the segment's end.

## Where the edge sits inside the protected angle

The four sections above put the edge on the larger of the two derived angles, which places the sigmoid's half-amplitude point on it and passes a quarter of the source intensity there. `global%source_filter_tolerance` names the largest source-intensity loss allowed anywhere inside that angle and moves the edge out until the transmission holds, and its default of 0.75 is the loss the older placement implies rather than a number anyone chose ([](fel-physics.md#sec-source-filter)). The reported split stays on the protected angle whatever the tolerance, so every row below is measured against one acceptance. The sweep is experiment `h` of `tests/scripts/startup_noise.py`.

```{table} The passband tolerance on the four machines, grid 256 on Aramis and 128 on the others, 1024 macroparticles per slice, one seed a row, each machine's rows against its own 0.75 row. The accepted power is the power inside the reported cut at the end of the first undulator, which on the two lines whose first segment is long in gain lengths is amplified light rather than the loading ([](#tab-sn-startup-four)).
:name: tab-sn-tolerance

| machine | tolerance | accepted at the first segment | wide/mode at the exit | gain length | saturation | pulse energy |
|---|---|---|---|---|---|---|
| Aramis | 0.75 | 16.4 MW | 1.66e-3 | 2.40 m | 32.9 m | 0.974 uJ |
| | 0.2 | +9.3% | x6.2 | +1.6% | 32.9 m | -2.0% |
| | 0.1 | +10.6% | x9.8 | -2.0% | 32.9 m | -2.1% |
| | 0.01 | +12.0% | x28.4 | -0.5% | 32.9 m | -1.8% |
| FLASH probe | 0.75 | 0.706 MW | 1.84e-3 | 1.02 m | 21.6 m | 112 uJ |
| | 0.2 | +11.0% | x6.0 | +0.8% | 21.6 m | -3.1% |
| | 0.1 | +13.7% | x9.3 | +1.2% | 21.6 m | -3.2% |
| | 0.01 | +19.4% | x24.9 | +2.1% | 21.6 m | -3.8% |
| FLASH1 | 0.75 | 0.935 MW | 2.23e-3 | 1.01 m | 20.7 m | 226 uJ |
| | 0.2 | +12.5% | x4.9 | +1.0% | 20.7 m | +0.8% |
| | 0.1 | +15.7% | x7.5 | +1.3% | 20.7 m | +1.2% |
| | 0.01 | +21.9% | x19.2 | +2.1% | 20.7 m | +1.4% |
| LCLS | 0.75 | 13.1 MW | 1.98e-3 | 3.64 m | 55.7 m | 145 uJ |
| | 0.2 | +10.1% | x6.1 | +1.5% | 55.7 m | +1.6% |
| | 0.1 | +11.6% | x9.5 | +2.0% | 55.7 m | +2.2% |
| | 0.01 | +13.1% | x26.3 | +3.9% | 55.7 m | +2.9% |
```

**What the sweep says.** The accepted power at the first segment rises with a tighter tolerance on every machine, by 9 to 13 percent at 0.2 and by 12 to 22 percent at 0.01. The wide-angle power over the same range rises faster, by 4.9 to 6.2 times at 0.2 and by 19 to 28 times at 0.01. The saturation point does not move on any machine at any tolerance. The fitted gain length stays inside 4 percent and the pulse energy inside 4 percent. Over this range and at this load the tolerance reaches the startup seed and leaves the gain alone, and what it buys in seed saturates while what it costs in wide-angle power does not.

**What it does not settle.** One seed a row. The interior mean of a 300-slice window fluctuates by about 25 percent on this page's own measurement, so a 2 percent move in pulse energy is inside the noise of one realization and only the trends that hold across all four machines carry. The fitted gain length comes from the records a decade above the starting power and below saturation, which on the six-segment lines is a handful of element ends. Choosing a default needs paired ensembles over seeds, priced against the seed saturating while the retained wide-angle power does not. The default stays at 0.75.

## The other remedies

Three remain, and none has been measured against the power inside the mode. Tanaka's coherent retrieval keeps the lowest Laguerre-Gauss orders of the bunching and adds the physical spontaneous emission back analytically ([](references.md#ref-tanaka)); this code carries it as `global%source_model = "coherent"`, which refuses a dark start for the reason the measurement above illustrates. Litvinenko's clones separate the induced from the spontaneous radiation by pairing each macroparticle with one of opposite charge ([](references.md#ref-litvinenko)). A deposition kernel of fixed physical width, wider than one cell, is the real-space form of the form factor of Pausch and co-workers ([](references.md#ref-pausch)). Any of the three changes the field the two codes compute from the same distributions, so each would be a switch, off by default, and measured against the power inside the mode rather than against Genesis4.
