---
title: "Slice migration"
short_title: Slice migration
tags: [explanation, record]
---

[Slice migration](glossary.md#gl-migration) is an operation of the run and not of the [load](glossary.md#gl-load). The key `run.migrate` selects it, and it is off by default. Its effect depends on the load: a load of independent particles is unchanged in kind by it, whereas the members of a [beamlet](glossary.md#gl-beamlet) are separated.

The page has two parts. The first describes the method: the transfer of a particle from one [slice](glossary.md#gl-slice) to another ([](#sec-migration)), and the effect of the transfer on a load of beamlets ([](#mig-beamlets)). This part follows from the method and from the source of the program, and it contains no result of a run. The second part, [](#mig-measured), records what was measured, and its numbers belong to the configurations named there. The loads are described in [](particle-load.md).

(sec-migration)=
## The transfer of a particle

Migration transfers a particle to another slice when its phase leaves the interval of its own slice. Because every particle has its own weight, the transfer does not require that a [macroparticle](glossary.md#gl-macroparticle) be one electron, the condition that Genesis4 names one4one: a particle that moves takes its charge with it. With $\theta$ of Eq. [](#eq-theta) of [](#sec-chart) and the interval $[0, 2\pi\,\texttt{n\_wavelength})$ of a slice, $a_{\mathrm{tar}} = \lfloor\theta/(2\pi\,\texttt{n\_wavelength})\rfloor$ is the number of slices by which the particle moves. A positive drift of $\theta$ is toward the higher index, which is the head. The coordinate $z$ of the particle is changed by exactly $-a_{\mathrm{tar}}\,\beta\,\Delta s$, which changes $\theta$ by $-a_{\mathrm{tar}}\cdot2\pi\,\texttt{n\_wavelength}$. For an integer `n_wavelength` the phase that enters every deposition is continuous through the move, to rounding. It is continuous at every harmonic $h$ as well, since $h\cdot2\pi\,\texttt{n\_wavelength}$ remains a multiple of $2\pi$, and Eq. [](#eq-zglobal) is invariant. Migration therefore requires an integer `n_wavelength`, and a value that is not an integer is refused at the first migration pass.

A particle is removed from its slice by exchange with the last particle of the slice, and the position is then examined again. A particle whose destination is beyond the [window](glossary.md#gl-window) is dropped, and its charge is counted. Migration is off by default. With the default load, which is a load of beamlets, it adds a radiation that is not that of the beam ([](#mig-beamlets)), and the Genesis4 [tiers](glossary.md#gl-tier) are compared with Genesis4 runs without one4one, in which no particle migrates.

Within a slice of several wavelengths a particle is not moved from one wavelength to the next, and nothing is lost by this. In such a window the slices are the only longitudinal positions that are represented. The spacing of the slices is therefore the only displacement that the discretization can express, and a particle that drifts within the interval $2\pi\,\texttt{n\_wavelength}$ of its slice has no finer position to which it could be moved. The position within the slice does not affect the current of the slice until the end of the interval is reached, and the error is less than the spacing of the slices, which is the resolution of the window.

The pass runs in threads. It leaves in every slice the particles that a serial scan of the slices in order leaves there, in the same order, so that a result does not depend on the number of threads. The serial scan is retained as `fel_migrate_slices_serial`, and under `run.migrate_check` both are run at every pass and their slices are compared bit for bit.

:::{dropdown} Genesis4 source
`fel_migrate_slices` is the weighted generalization of [`Sorting::localSort`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L74): its criterion in Lucifer's chart, with the same sample-scaled window, swap-with-last removal and rescan. Genesis4 discards the particles that leave at the world edges silently ([lines 194 to 195](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L194-L195)). Here they are dropped with their charge counted. Genesis4 migrates only under one4one, and Lucifer's weights lift that condition.
:::

The measured levels and how they are checked are in [](validation.md#val-slice-migration-under-weights).

(mig-beamlets)=
## Migration of a load of beamlets

Migration moves each macroparticle singly. It therefore separates the members of a beamlet, and the two slices that then contain members of one beamlet are no longer quiet. The radiation that results is not that of the beam. At a given size of a beamlet its power is inversely proportional to the number $N$ of macroparticles per slice. The fraction of the beamlets that a pass separates is a property of the beam and not of the load, so that the number of separated beamlets is proportional to $N$. The members of a separated beamlet that remain in one slice contribute to the bunching of that slice a phasor whose magnitude is of the order of $1/N$, and the phasors of different beamlets are independent, so that the mean square of the bunching is proportional to $N \cdot N^{-2} = 1/N$. The coefficient depends on the size of a beamlet through the sum over the members that remain together, and this argument does not determine it. With a few hundred to several thousand macroparticles per slice the power was found to lie between one percent of the SASE power and nearly the whole of it in the first part of a line, and the amplified radiation exceeds it further along ([](#pl-migration-radiation)). The power of a run with migration on and `shot_noise = false` is this radiation alone, so that such a run determines it for a given input file.

Genesis4 migrates only when each macroparticle is one electron, in which case there is no beamlet to separate, and Hwang and Qiang move the members of a beamlet together ([](references.md#ref-hwang-qiang)). The electrons mode ([](#sec-electrons)) is the first of these cases: its particles are independent before and after a move. The pulse energy of that mode with migration on does not differ from its value with migration off by more than the uncertainty of the comparison ([](#pl-migration-electrons)).

What follows for a run with migration is stated here in words, and the numbers are in [](#pl-migration-radiation). A load of independent particles, which is the electrons mode, is the load for which migration is intended. With a load of beamlets, in the one configuration in which the number of macroparticles was scanned, the radiation that migration adds falls below one percent of the SASE power only at about ten thousand macroparticles per slice or more. At a given number it depends on the configuration and on the position along the line. The steps of a run with migration are in [](#ug-migration).

(mig-keep)=
### The copies of the keep mode

In the keep mode under the [quiet start](glossary.md#gl-quiet-start) the $b$ copies of a particle are placed over one wavelength from the phase $\theta_0$ of the particle ([](#sec-loading)). A copy at $\theta_0 + 2\pi j/b \ge 2\pi\,\texttt{n\_wavelength}$ is loaded beyond the end of its slice. The criterion of [](#sec-migration) transfers it to the next slice at the first pass, which follows the first element of the line, unless the drift within that element has returned it to its slice. The beamlet of that particle is then separated although no particle has drifted. A particle has a copy beyond the end when $\theta_0$ lies within $2\pi(b-1)/b$ of the end of its slice. For phases distributed uniformly over the slice this is the fraction $(1 - 1/b)/\texttt{n\_wavelength}$ of the beamlets, which is seven eighths of them for a beamlet of 8 in a slice of one wavelength.

With migration on, the keep mode under the quiet start therefore adds the radiation of separated beamlets from the first pass onward, and the weight of a separated member is the charge of a particle of the bunch divided by $b$. The sample mode places every beamlet within its slice, so that its beamlets are separated by drift alone. The keep mode without the quiet start has no copies, and the electrons mode has no beamlets.

This separation follows from the placement of the copies and not from the transfer. The averaged step depends on $z$ only through $\theta$ modulo $2\pi$, so that the copies of a particle could be placed within the wavelength that contains the particle, with the same bunching of the beamlet at every harmonic, and no copy would then lie beyond the end of its slice. The loader does not at present place them so. The radiation that the keep mode then adds under migration has not been measured.

(mig-measured)=
## Measurements

The sections above describe the method and hold for any run. This section records measurements. They were made before the manual separated its measurements from its explanations, and they are kept as they were recorded, each with the day on which it entered the manual. A field that a record does not state, such as the version of the program or the number of [random seeds](glossary.md#gl-seed-random), was not recorded with it.

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
