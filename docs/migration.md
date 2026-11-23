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

A particle is removed from its slice by exchange with the last particle of the slice, and the position is then examined again. A particle whose destination is beyond the [window](glossary.md#gl-window) is dropped, and its charge is counted. Migration is off by default, since the Genesis4 [tiers](glossary.md#gl-tier) are compared with Genesis4 runs without one4one, in which no particle migrates.

Within a slice of several wavelengths a particle is not moved from one wavelength to the next, and nothing is lost by this. In such a window the slices are the only longitudinal positions that are represented. The spacing of the slices is therefore the only displacement that the discretization can express, and a particle that drifts within the interval $2\pi\,\texttt{n\_wavelength}$ of its slice has no finer position to which it could be moved. The position within the slice does not affect the current of the slice until the end of the interval is reached, and the error is less than the spacing of the slices, which is the resolution of the window.

The pass runs in threads. It leaves in every slice the particles that a serial scan of the slices in order leaves there, in the same order, so that a result does not depend on the number of threads. The serial scan is retained as `fel_migrate_slices_serial`, and under `run.migrate_check` both are run at every pass and their slices are compared bit for bit.

(mig-beamlets)=
## Migration of a load of beamlets

Migration moves each macroparticle singly. It therefore separates the members of a beamlet, and the two slices that then contain members of one beamlet are no longer quiet. The radiation that results is not that of the beam. Its power is inversely proportional to the number $N$ of macroparticles per slice. The fraction of the beamlets that a pass separates is a property of the beam and not of the load, so that the number of separated members is proportional to $N$. Each contributes to the bunching of its slice a phasor of magnitude $1/N$, and the phases of the members of different beamlets are independent, so that the mean square of the bunching is proportional to $N \cdot N^{-2} = 1/N$. With one thousand to eight thousand macroparticles per slice the power was found to lie between one percent and a few tens of percent of the SASE power in the first part of a line, and the amplified radiation exceeds it further along ([](#pl-migration-radiation)). The power of a run with migration on and `shot_noise = false` is this radiation alone, so that such a run determines it for a given input file.

Genesis4 migrates only when each macroparticle is one electron, in which case there is no beamlet to separate, and Hwang and Qiang move the members of a beamlet together ([](references.md#ref-hwang-qiang)). The electrons mode ([](#sec-electrons)) is the first of these cases: its particles are independent before and after a move. The pulse energy of that mode with migration on does not differ from its value with migration off by more than the spread among seeds ([](#pl-migration-electrons)).

:::{dropdown} Genesis4 source
`fel_migrate_slices` is the weighted generalization of [`Sorting::localSort`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L74): its criterion in Lucifer's chart, with the same sample-scaled window, swap-with-last removal and rescan. Genesis4 discards the particles that leave at the world edges silently ([lines 194 to 195](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/Sorting.cpp#L194-L195)). Here they are dropped with their charge counted. Genesis4 migrates only under one4one, and Lucifer's weights lift that condition.
:::

The measured levels and how they are checked are in [](validation.md#val-slice-migration-under-weights).

(mig-measured)=
## Measurements

The sections above describe the method and hold for any run. This section records measurements. They were made before the manual separated its measurements from its explanations, and they are kept as they were recorded, each with the day on which it entered the manual. A field that a record does not state, such as the version of the program or the seeds, was not recorded with it.

(pl-migration-radiation)=
### The radiation that migration adds to a load of beamlets

Recorded on 2026-11-12, with the program at commit `31be058`. The radiation was measured with a load of beamlets of 4 macroparticles that has neither shot noise nor a seed, which radiates $10^{-18}$ to $10^{-15}$ W with migration off. With migration on, its power was compared with the SASE power of the same configuration, which is the mean of 8 seeds with migration off.

The first two rows are of the example `sase`. The rows LCLS and FLASH1 were recorded under the names of the two FEL lines, whose configurations in the program are the examples `lcls` and `flash1`. The input files of those runs were not recorded.

| configuration | macroparticles per slice | position | power with migration on, as a fraction of the SASE power |
|---|---|---|---|
| the example `sase` | 1024 | 5 to 40 m | 0.09 to 0.21 |
| the example `sase` | 8192 | 5 to 40 m | 0.011 to 0.027 |
| LCLS | 4096 | 10 to 20 m | 0.06 |
| LCLS | 4096 | 40 m | 0.01 |
| FLASH1 | 4096 | 8 m | 0.38 |
| FLASH1 | 4096 | 14 m | 0.05 |

The two rows of the example `sase` are consistent with a power that is inversely proportional to the number of macroparticles, as [](#sec-migration) derives. In the rows of LCLS and of FLASH1 the fraction falls along the line, as the amplified radiation exceeds the radiation that migration adds.

With migration on, the pulse energy of FLASH1 at the exit, over the same 8 seeds, is $0.966 \pm 0.015$ of its value with migration off. The power of a seeded FLASH1 run with migration on, over 4 seeds, differs from that with migration off by 1.5% at 6 m and by less than 0.7% from 10 to 18 m. The commit of the program was not recorded for the seeded runs.

(pl-migration-electrons)=
### Migration of independent particles

Recorded on 2026-11-21. In the window of the example `sase`, over 16 seeds, the pulse energy with migration on was divided by the pulse energy with migration off, for the electrons mode and for the weight load of the sample mode with 2048 macroparticles per slice.

| position | electrons mode | weight load, 2048 macroparticles per slice |
|---|---|---|
| 4 m | $1.00 \pm 0.05$ | $1.00 \pm 0.06$ |
| 20 m | $1.02 \pm 0.10$ | $1.06 \pm 0.14$ |
| the exit | $1.08 \pm 0.25$ | $1.08 \pm 0.20$ |
