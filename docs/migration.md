---
title: "Slice migration"
short_title: Slice migration
tags: [explanation, record]
---

[Slice migration](glossary.md#gl-migration) is an operation of the run and not of the [load](glossary.md#gl-load). The key `run.migrate` selects it, and it is off by default. Its effect depends on the load. A load of independent particles remains a load of independent particles after a transfer, whereas the members of a [beamlet](glossary.md#gl-beamlet) are separated and no longer cancel.

The page has two parts. The first describes the method: the transfer of a particle from one [slice](glossary.md#gl-slice) to another ([](#sec-migration)), and the effect of the transfer on a load of beamlets ([](#mig-beamlets)). This part follows from the method and from the source of the program, and it contains no result of a run. The second part, [](#mig-measured), records what was measured, and its numbers belong to the configurations named there. The loads are described in [](particle-load.md).

(sec-migration)=
## The transfer of a particle

Migration transfers a particle to another slice when its phase leaves the interval of its own slice. Because every particle has its own weight, the transfer does not require that a [macroparticle](glossary.md#gl-macroparticle) be one electron, the condition that Genesis4 names one4one: a particle that moves takes its charge with it. With $\theta$ of Eq. [](#eq-theta) of [](#sec-chart) and the interval $[0, 2\pi\,\texttt{n\_wavelength})$ of a slice, $a_{\mathrm{tar}} = \lfloor\theta/(2\pi\,\texttt{n\_wavelength})\rfloor$ is the number of slices by which the particle moves. A positive drift of $\theta$ is toward the higher index, which is the head. The coordinate $z$ of the particle is changed by exactly $-a_{\mathrm{tar}}\,\beta\,\Delta s$, which changes $\theta$ by $-a_{\mathrm{tar}}\cdot2\pi\,\texttt{n\_wavelength}$. For an integer `n_wavelength` the phase that enters every deposition is continuous through the move, to rounding. It is continuous at every harmonic $h$ as well, since $h\cdot2\pi\,\texttt{n\_wavelength}$ remains a multiple of $2\pi$, and Eq. [](#eq-zglobal) is invariant. Migration therefore requires an integer `n_wavelength`, and a value that is not an integer is refused at the first migration pass.

A particle is removed from its slice by exchange with the last particle of the slice, and the position is then examined again. A particle whose destination is beyond the [window](glossary.md#gl-window) is dropped, and its charge is counted.

Within a slice of several wavelengths a particle is not moved from one wavelength to the next, and nothing is lost by this. In such a window the slices are the only longitudinal positions that are represented. The spacing of the slices is therefore the only displacement that the discretization can express, and a particle that drifts within the interval $2\pi\,\texttt{n\_wavelength}$ of its slice has no finer position to which it could be moved. The position within the slice does not affect the current of the slice until the end of the interval is reached, and the error is less than the spacing of the slices, which is the resolution of the window.

The pass runs in threads. It leaves in every slice the particles that a serial scan of the slices in order leaves there, in the same order, so that a result does not depend on the number of threads. The serial scan is retained as `fel_migrate_slices_serial`, and under `run.migrate_check` both are run at every pass and their slices are compared bit for bit.

Migration is off by default, for two reasons. With the default load, which is a load of beamlets, it adds a radiation that is not that of the beam ([](#mig-beamlets)). The Genesis4 [tiers](glossary.md#gl-tier) are compared with Genesis4 runs without one4one, in which no particle migrates.

:::{dropdown} Genesis4 source
`fel_migrate_slices` is the weighted generalization of [`Sorting::localSort`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L74): its criterion in Lucifer's chart, with the same sample-scaled window, swap-with-last removal and rescan. Genesis4 discards the particles that leave at the world edges silently ([lines 194 to 195](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L194-L195)). Here they are dropped with their charge counted. Genesis4 migrates only under one4one, and Lucifer's weights lift that condition.
:::

The measured levels and how they are checked are in [](validation.md#val-slice-migration-under-weights).

(mig-beamlets)=
## Migration of a load of beamlets

Migration moves each macroparticle singly. It therefore separates the members of a beamlet, and the two slices that then contain members of one beamlet are no longer quiet. The radiation that results is not that of the beam. At a given size of a beamlet its power is inversely proportional to the number $N$ of macroparticles per slice. The fraction of the beamlets that a pass separates is a property of the beam and not of the load, so that the number of separated beamlets is proportional to $N$. The members of a separated beamlet that remain in one slice contribute to the bunching of that slice a phasor whose magnitude is of the order of $1/N$, and the phasors of different beamlets are independent, so that the mean square of the bunching is proportional to $N \cdot N^{-2} = 1/N$. The coefficient depends on the size of a beamlet through the sum over the members that remain together, and this argument does not determine it. With a few hundred to several thousand macroparticles per slice the power was found to lie between one percent of the SASE power and nearly the whole of it in the first part of a line, and the amplified radiation exceeds it further along ([](#pl-migration-radiation)). The power of a run with migration on and `shot_noise = false` is this radiation alone, so that such a run determines it for a given input file.

Genesis4 migrates only when each macroparticle is one electron, in which case there is no beamlet to separate, and Hwang and Qiang move the members of a beamlet together ([](references.md#ref-hwang-qiang)). The electrons mode ([](#sec-electrons)) is the first of these cases: its particles are independent before and after a move. The pulse energy of that mode with migration on does not differ from its value with migration off by more than the uncertainty of the comparison ([](#pl-migration-electrons)).

What follows for a run with migration is stated here in words, and the numbers are in [](#pl-migration-radiation). A load of independent particles, which is the electrons mode, is the load for which migration is intended. With the beamlets of the sample mode, the one scan of the number of macroparticles per slice ends at 8192, where the whole range of the radiation that migration adds is still above one percent of the SASE power. The inverse proportionality, continued beyond the scan, places one percent near ten thousand macroparticles per slice for the lower end of the range and near twenty thousand for the upper end. At a given number the radiation depends on the configuration and on the position along the line. In the keep mode under the [quiet start](glossary.md#gl-quiet-start) the number of macroparticles per slice is that of the bunch and is not chosen by a key of the load ([](#mig-keep)). In the one configuration measured, the radiation that migration adds to that mode was between two percent and nine percent of the SASE energy ([](#mig-keep-measured)). The steps of a run with migration are in [](#ug-migration).

(mig-keep)=
### The copies of the keep mode

In the keep mode under the quiet start the $b$ copies of a particle are placed within the wavelength that contains the particle ([](#sec-loading)), and no copy is stored outside the slice of its particle. Every load thus places its beamlets within their slices. The sample mode distributes the members of a beamlet over the periods of its slice, the keep mode without the quiet start has no copies, and the electrons mode has no beamlets. With migration on, a beamlet of any load is therefore separated only by the drift of its members.

The members of a beamlet of the keep mode have the same energy and the same transverse coordinates. Where the interaction with the field is weak they drift by the same phase $\delta$ between two passes, and a pass transfers the members that lie within $|\delta|$ of the end of the slice toward which they drift. For $|\delta| < 2\pi$ these are members of the beamlets of the wavelength at that end, which are the fraction $1/\texttt{n\_wavelength}$ of the beamlets where the current is uniform over the slice. Such a beamlet loses $b|\delta|/2\pi$ of its members on average, so that a pass transfers the fraction $|\delta|/(2\pi\,\texttt{n\_wavelength})$ of the macroparticles. The beamlet loses at least one member, and is separated, when $2\pi/b \le |\delta| < 2\pi$. In a slice of one wavelength every beamlet extends over the whole slice, and a drift of that size separates every beamlet at one pass.

The weight of a separated member is the charge of a particle of the bunch divided by $b$, and the number of macroparticles per slice is $b$ times the number of particles of the bunch that arrived in the slice. The radiation of [](#mig-beamlets) is therefore lowered by a bunch of more particles and not by a key of the load. In the one configuration measured it was a few percent of the SASE energy after two undulators, and it rose along the line to nearly a tenth of the SASE energy, in slices of three wavelengths and in slices of one ([](#mig-keep-measured)).

(mig-measured)=
## Measurements

The sections above describe the method and hold for any run. This section records measurements. The first two records were made before the manual separated its measurements from its explanations, and they are kept as they were recorded, each with the day on which it entered the manual. A field that one of them does not state, such as the version of the program or the number of [random seeds](glossary.md#gl-seed-random), was not recorded with it. The third record, [](#mig-keep-measured), states the fields of its measurement.

(pl-migration-radiation)=
### The radiation that migration adds to a load of beamlets

Recorded on 2026-11-12, with the program at commit `31be058`. The radiation was measured with a load of beamlets of 4 macroparticles that has neither shot noise nor a [seed field](glossary.md#gl-seed-field), which radiates $10^{-18}$ to $10^{-15}$ W with migration off. With migration on, its power was compared with the SASE power of the same configuration, which is the mean of 8 random seeds with migration off.

The first six rows are of the example `sase`. The rows LCLS and FLASH1 were recorded under the names of the two FEL lines, whose configurations in the program are the examples `lcls` and `flash1`. The input files of those runs were not recorded.

The manual first held the rows of 1024 and of 8192 macroparticles per slice of the example `sase`, the rows of LCLS and the rows of FLASH1 at 8 m and at 14 m, with no commit and no position for the example `sase`. The commit, the beamlet size, the numbers of random seeds, the positions, the four other rows of the example `sase` and the row of FLASH1 at 12 m were added on 2026-11-22 and 2026-11-23 from the notes made on the day of the measurement. Every fraction is as it was recorded on that day.

Each range of the example `sase` is over the positions from 5 to 40 m. The position at which a range has its upper end was not recorded, so that the dependence on position within a range is not known.

| configuration | macroparticles per slice | position | power with migration on, as a fraction of the SASE power |
|---|---|---|---|
| the example `sase` | 256 | 5 to 40 m | 0.38 to 0.85 |
| the example `sase` | 512 | 5 to 40 m | 0.19 to 0.40 |
| the example `sase` | 1024 | 5 to 40 m | 0.09 to 0.21 |
| the example `sase` | 2048 | 5 to 40 m | 0.05 to 0.10 |
| the example `sase` | 4096 | 5 to 40 m | 0.02 to 0.05 |
| the example `sase` | 8192 | 5 to 40 m | 0.011 to 0.027 |
| LCLS | 4096 | 10 to 20 m | 0.06 |
| LCLS | 4096 | 40 m | 0.01 |
| FLASH1 | 4096 | 8 m | 0.38 |
| FLASH1 | 4096 | 12 m | 0.11 |
| FLASH1 | 4096 | 14 m | 0.05 |

In the six rows of the example `sase` the lower end of the range multiplied by the number of macroparticles per slice lies between 82 and 102, and the upper end between 205 and 221. The rows are therefore consistent with a power that is inversely proportional to the number of macroparticles, as [](#mig-beamlets) derives. At 512 macroparticles per slice beamlets of 2 gave the fraction of beamlets of 4. Beamlets of 8 gave less in this example, 0.03 to 0.045 at 2048 macroparticles per slice and 0.005 to 0.012 at 8192, and the same fractions as beamlets of 4 for LCLS and for FLASH1. The argument of [](#mig-beamlets) gives the dependence on the number of macroparticles at a fixed size of a beamlet and does not determine the dependence on that size. The lower fraction of beamlets of 8 in the example `sase`, which LCLS and FLASH1 do not show, is not explained. In the rows of LCLS and of FLASH1 the fraction falls along the line, as the amplified radiation exceeds the radiation that migration adds.

With migration on, the SASE power of FLASH1 over the same 8 random seeds is 1.30 of its value with migration off at 8 m, 1.07 at 12 m and 1.01 at 18 m, and the pulse energy at the exit is $0.966 \pm 0.015$ of its value with migration off. The pulse energy at the exit is thus lower with migration on by more than twice its uncertainty, on the same random seeds, whereas the radiation that migration adds would raise it. The charge that was dropped beyond the window was not recorded, and the difference is not explained. The power of a FLASH1 run with a seed field and with migration on, over 4 random seeds, differs from that with migration off by 1.5% at 6 m and by less than 0.7% from 10 to 18 m. The commit of the program was not recorded for the runs with a seed field.

(pl-migration-electrons)=
### Migration of independent particles

Recorded on 2026-11-21. In the window of the example `sase`, over 16 random seeds, the pulse energy with migration on was divided by the pulse energy with migration off, for the electrons mode and for the weight load of the sample mode with 2048 macroparticles per slice.

| position | electrons mode | weight load, 2048 macroparticles per slice |
|---|---|---|
| 4 m | $1.00 \pm 0.05$ | $1.00 \pm 0.06$ |
| 20 m | $1.02 \pm 0.10$ | $1.06 \pm 0.14$ |
| the exit | $1.08 \pm 0.25$ | $1.08 \pm 0.20$ |

The record does not state whether an uncertainty is a standard deviation or a standard error of the mean. The pulse energy of a single random seed at the exit has a standard deviation of 0.6 to 0.7 of its mean ([](electrons-mode.md#el-compare)), which gives a standard error near 0.2 for the ratio of two means of 16 random seeds, and that is the uncertainty of the last row.

The electrons mode shows no difference outside the uncertainty, and the weight load shows none either. In [](#pl-migration-radiation) the radiation that migration adds to a load of 2048 macroparticles per slice with beamlets of 4 is 0.05 to 0.10 of the SASE power. The uncertainty of the ratio at 20 m and at the exit does not resolve a fraction of that size, so that this record neither shows that radiation nor excludes it. The beamlet size of the weight load was not recorded.

(mig-keep-measured)=
### The radiation that migration adds to the keep mode under the quiet start

Measured on 2026-11-24 on the example `import`, with the release build of Lucifer at commit `b6392a4` and with the same source in which `keep_initial_state` places the copies as [](#mig-keep) states. Both were built with gfortran 15.3.0 against Bmad 20261122-0 and were run on four threads. At commit `b6392a4` copy $j$ of a particle of phase $\theta_0$ was placed at $\theta_0 + 2\pi j/b$ and was not returned to the wavelength of the particle. A copy of a particle in the last wavelength of a slice was then stored beyond the end of the slice, and the first pass transferred it. The two placements are called the former and the present one in this record.

The input file of the example loads a bunch of 50000 particles of 3.76 electrons each as 400000 macroparticles in beamlets of 8, with no seed field. Three runs were made of it, each with the keys stated and no other change: the quiet load with migration on, `beam.shot_noise=false run.migrate=true`, the quiet load with migration off, `beam.shot_noise=false run.migrate=false`, and the SASE run, `run.migrate=false`. The first and the third were made at the eight values 12345 to 12352 of `run.ran_seed`, which generates the bunch as well as the noise, and the second at 12345. Each was made in slices of three wavelengths, which is the `window.n_wavelength` of the example, and in slices of one wavelength. Over the eight random seeds the window is 32 to 36 slices of three wavelengths and 94 to 108 slices of one.

The step of the example is three undulator periods, `ds_step = 0.045` of the element `UND` of `examples/aramis.bmad`, at which the field advances by one slice of three wavelengths. In slices of one wavelength the field advances three slices at such a step, and the mean SASE pulse energy at the exit was then 42.6 pJ, three times that of the slices of three wavelengths. The runs in slices of one wavelength were therefore made with `ds_step = 0.015` in that lattice, a step of one period, and every number below for that slice length is from those runs.

A pulse energy at a position is the energy of the field in the window together with the energy that has left the window before that position. The mean SASE pulse energy at the exit is 14.3 pJ in slices of three wavelengths and 14.9 pJ in slices of one. At a step of three periods the SASE runs of the two placements agree at every record of the line to $6\times10^{-11}$ of the largest pulse energy. The SASE runs at a step of one period were made with the former placement alone. The quiet load with migration off has a pulse energy at the exit of $6\times10^{-35}$ J in slices of three wavelengths and of $3\times10^{-34}$ J in slices of one, with either placement.

The first table gives the number of macroparticles that one pass transferred, as its range over the eight random seeds.

| wavelengths per slice | placement | the pass at the end of the first undulator | a pass at the end of a later undulator |
|---|---|---|---|
| 3 | former | 45356 to 46382 | 33636 to 34909 |
| 3 | present | 33645 to 34504 | 33401 to 34695 |
| 1 | former | 136238 to 137582 | 100186 to 100734 |
| 1 | present | 100291 to 100681 | 100168 to 100737 |

With the present placement the first pass transfers what a later pass transfers. That number is one twelfth of the macroparticles in slices of three wavelengths and one quarter in slices of one, which by [](#mig-keep) is a drift of a quarter of a wavelength in one undulator at either slice length.

The second table gives the pulse energy of the quiet load with migration on. The columns of the two placements are that energy divided by the SASE pulse energy at the same position and with the same random seed, as the mean over the eight random seeds and its standard error. The last column is the energy with the present placement divided by that with the former, seed by seed.

| wavelengths per slice | end of undulator | position | former placement | present placement | present divided by former |
|---|---|---|---|---|---|
| 3 | 2 | 8.74 m | $0.0273 \pm 0.0010$ | $0.0198 \pm 0.0008$ | $0.727 \pm 0.003$ |
| 3 | 3 | 13.49 m | $0.0357 \pm 0.0010$ | $0.0315 \pm 0.0009$ | $0.882 \pm 0.002$ |
| 3 | 5 | 22.99 m | $0.0635 \pm 0.0020$ | $0.0559 \pm 0.0017$ | $0.881 \pm 0.005$ |
| 3 | 8 | 37.24 m | $0.0876 \pm 0.0030$ | $0.0850 \pm 0.0029$ | $0.970 \pm 0.006$ |
| 3 | 12 | 56.24 m | $0.0902 \pm 0.0031$ | $0.0907 \pm 0.0032$ | $1.005 \pm 0.007$ |
| 1 | 2 | 8.74 m | $0.0276 \pm 0.0028$ | $0.0203 \pm 0.0020$ | $0.736 \pm 0.003$ |
| 1 | 3 | 13.49 m | $0.0335 \pm 0.0026$ | $0.0298 \pm 0.0023$ | $0.892 \pm 0.004$ |
| 1 | 5 | 22.99 m | $0.0554 \pm 0.0023$ | $0.0494 \pm 0.0020$ | $0.892 \pm 0.005$ |
| 1 | 8 | 37.24 m | $0.0787 \pm 0.0034$ | $0.0774 \pm 0.0033$ | $0.985 \pm 0.007$ |
| 1 | 12 | 56.24 m | $0.0832 \pm 0.0036$ | $0.0846 \pm 0.0034$ | $1.018 \pm 0.006$ |

The file `docs/figures/migration-keep.tsv` holds these quantities at the end of every undulator from the second to the twelfth, with the mean pulse energies.

The present placement lowers the energy at the end of the second undulator by 27% in slices of three wavelengths and by 26% in slices of one. The difference decreases along the line. At the exit the two placements do not differ in slices of three wavelengths, and the present one is higher by $1.8 \pm 0.6$ percent in slices of one, which is not explained. The radiation that migration adds to this load is therefore, from the third undulator onward, mainly that of the beamlets that the drift has separated. At the exit it is 9% of the SASE pulse energy in slices of three wavelengths and 8% in slices of one.
