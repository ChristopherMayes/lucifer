---
title: References
short_title: References
---

This page lists what this code implements, and where each thing came from. Every entry names the source and the routine that carries it, so a reader can go from a paper to the transcription and back. The physics manual repeats these citations in place, section by section, with the equations beside them.

## The codes this one is measured against

(ref-genesis4)=
**Genesis 1.3 Version 4 (Genesis4)**, S. Reiche. The averaged path's physics was transcribed from it, routine by routine: the equations of motion and their integrator, the coupling factors, the transform field solver and its [source filter](glossary.md#gl-source-filter), the slippage schedule, Fawley's [load](glossary.md#gl-load), the distribution [resampler](glossary.md#gl-resampler), [slice migration](glossary.md#gl-migration)y.md#gl-slice) migration, the wakes and the space-charge solvers. The Genesis4 [tiers](glossary.md#gl-tier) compare against it from shared dumps. Genesis4 is GPLv3, which permits the transcription, and the per-routine citations in the source are kept for that reason. The manual's collapsed Genesis4 source blocks ([](fel-physics.md)) link each original at release 4.6.16. Not transcribed, and with no Genesis4 original: the unaveraged method, the spectral field propagator ([](spectral-field.md)), the [device](glossary.md#gl-device) step, the numerics derived from the beam, the load of an arbitrary bunch with its per-particle weights, the escaped-field [bank](glossary.md#gl-bank) and the pulse, and the closed forms of [](shot-noise.md). The lattice and the transport outside the FEL step are Bmad's. The original code is described in S. Reiche, *GENESIS 1.3: a fully 3D time-dependent FEL simulation code*, [Nucl. Instr. and Meth. A **429** (1999) 243](<https://doi.org/10.1016/S0168-9002(99)00114-X>). The current code, its [manual](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/manual/MAIN_INPUT.md) and its [changelog](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/CHANGELOG.md) are on GitHub, linked here at release 4.6.16. Release 4.6.12 added an FFT field solver and the `source_filter` option, a sigmoid filter on the transformed source term. See [](startup-noise.md#prior-work).

**Bmad**, D. Sagan. The lattice, the tracking of every non-FEL element, the wake machinery, the beam I/O and the units and phase-space conventions are Bmad's. D. Sagan, *Bmad: A relativistic charged particle simulation library*, [Nucl. Instr. and Meth. A **558** (2006) 356](https://doi.org/10.1016/j.nima.2005.11.001).

(ref-minerva)=
**MINERVA**, H. P. Freund and P. J. M. van der Slot. H. P. Freund, P. J. M. van der Slot, D. L. A. G. Grimminck, I. D. Setija and P. Falgari, [New J. Phys. **19**, 023020 (2017)](https://doi.org/10.1088/1367-2630/aa59f1). The published existence proof that unaveraged FEL dynamics is practical for production. It is also the published source of the substeps-per-period envelope that this mode's own convergence measurement is read against. Lucifer differs by keeping the grid field and the Lorentz force where MINERVA evaluates modal fields. Only the published work informed this port. MINERVA's source was not read, and no run of it was compared against. See [](fel-physics.md#sec-unaveraged).

**SIMPLEX**, T. Tanaka. Its source is unlicensed and unread. The coherent source is implemented from the published paper alone.

**openPMD-beamphysics** is a library the [validation suite](glossary.md#gl-validation-suite) calls: it converts the Genesis4 reference dumps, and the [wavefront](glossary.md#gl-wavefront) round trip goes through its openPMD I/O in both directions.

## Papers the code implements

(ref-tanaka)=
**Coherent source.** T. Tanaka, [Phys. Rev. Accel. Beams **27**, 030703 (2024)](https://doi.org/10.1103/PhysRevAccelBeams.27.030703), [arXiv:2310.20197](https://arxiv.org/abs/2310.20197). The [slice](glossary.md#gl-slice) bunch factor splits into coherent and spatially incoherent parts, and the incoherent part's radiation diffracts away. It is implemented from the paper in `fel_track_mod` (`run.source_model = 'coherent'`), which refuses the configurations the model does not cover. See [](fel-physics.md#sec-coherent-source).

(ref-bane-stupakov)=
**Resistive-wall wakes.** K. L. F. Bane and G. Stupakov, [SLAC-PUB-10707](https://www.slac.stanford.edu/pubs/slacpubs/10500/slac-pub-10707.pdf). The numerical impedance with AC (Drude) conductivity and Leontovich surface impedance, round and flat geometry. It is computed in `fel_collective_mod` (`fel_resistive_wall_wake`) on the integration grids Genesis4 uses. See [](fel-physics.md#sec-wakes).

(ref-fawley)=
**Shot-noise loading.** W. M. Fawley's algorithm for loading shot-noise microbunching in FEL simulation codes, [Phys. Rev. ST Accel. Beams **5**, 070701 (2002)](https://doi.org/10.1103/PhysRevSTAB.5.070701), in the one-substitution weighted generalization that makes it correct for per-particle weights. `fel_fawley_noise`, checked statistically against the closed form of its kick under uniform and nonuniform weights. See [](fel-physics.md#sec-noise) and [](shot-noise.md).

(ref-mcneil)=
**Shot noise through macroparticle charge.** B. W. J. McNeil, M. W. Poole and G. R. M. Robb, [Phys. Rev. ST Accel. Beams **6**, 070701 (2003)](https://doi.org/10.1103/PhysRevSTAB.6.070701), and their conference paper, [Proc. PAC 2003, p. 950](https://accelconf.web.cern.ch/p03/PAPERS/MPPB059.PDF). A [macroparticle](glossary.md#gl-macroparticle) stands for a Poisson number of electrons, so the statistics of the load are those of the electrons' arrival and not of an averaged bunching. `fel_weight_noise` was written from the conference paper and from the restatement below, with the bunch's charge held. See [](shot-noise.md#shot-weight).

(ref-hwang-qiang-weight)=
**The weight load on a beamlet.** K. Hwang and J. Qiang, [Phys. Rev. Accel. Beams **21**, 120702 (2018)](https://doi.org/10.1103/PhysRevAccelBeams.21.120702), Section III. Reviews Fawley's phase perturbation and McNeil, Poole and Robb's charge-weight perturbation on a [beamlet](glossary.md#gl-beamlet) of evenly spaced phases, and shows the second has the shot-noise bunching as an equality. It divides a beamlet into twice as many particles as the highest harmonic of the simulation, the size [](shot-noise.md#shot-size) measures to be the smallest that is usable. The same paper is cited below for its migration scheme.

(ref-hormann)=
**Binomial random numbers.** W. Hörmann, "The generation of binomial random variates", J. Stat. Comput. Simul. **46**, 101 (1993), algorithm BTRS. `fel_ran_binomial`, which draws the weight load's counts.

(ref-watson)=
**Weber's exponential integrals.** G. N. Watson, *A Treatise on the Theory of Bessel Functions*, 2nd edition (Cambridge University Press, 1944), sections 13.3 and 13.31. The two integrals of Bessel functions against a Gaussian that give the mean of the shot-noise kick in closed form. See [](shot-noise.md).

**Incoherent (quantum) diffusion.** The standard undulator quantum-diffusion variance in the Saldin form, as fitted in Genesis4's [`Incoherent.cpp`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Core/Incoherent.cpp#L56), with one draw per [beamlet](glossary.md#gl-beamlet) so the [quiet start](glossary.md#gl-quiet-start)'s cancellation survives. Applied only when Bmad's own `radiation_fluctuations` switch is on. See [](fel-physics.md#sec-eom).

(ref-sase-theory)=
**SASE startup and gain estimates.** Ming Xie's fitting formula for the 3D power gain length, [Nucl. Instrum. Methods A **445**, 59 (2000)](<https://doi.org/10.1016/S0168-9002(00)00114-5>). E. L. Saldin, E. A. Schneidmiller and M. V. Yurkov's effective shot-noise power, saturation length and efficiency, and the incoherent undulator power of a beam in its central cone, [New J. Phys. **12**, 035010 (2010)](https://doi.org/10.1088/1367-2630/12/3/035010), [arXiv:0912.4161](https://arxiv.org/abs/0912.4161), eqs. 1, 5 and 18. They are computed from the input's parameters by the script of [](startup-noise.md) and placed beside the measurements there.

(ref-wolski)=
**Normal-mode emittances.** A. Wolski's approach to general coupled linear optics, [Phys. Rev. ST Accel. Beams **9**, 024001 (2006)](https://doi.org/10.1103/PhysRevSTAB.9.024001), supplies the eigen-emittances stored beside the projected Twiss planes. Bmad's `calc_bunch_params` evaluates them, and `tests/scripts/bunch_params_from_stats.py` re-derives them independently as a check. See [](fel-physics.md#sec-stats).

## Papers the SASE convergence page cites

These are the sources of [](startup-noise.md#prior-work), each read for that page.

(ref-penman-mcneil)=
**Shot noise in one dimension.** C. Penman and B. W. J. McNeil, [Opt. Commun. **90**, 82 (1992)](<https://doi.org/10.1016/0030-4018(92)90333-M>). The shot noise to impose on the electron phases as the initial condition of a simulation, from the Poisson statistics of the real electrons.

(ref-mcneil-2003)=
**Shot noise in unaveraged models.** B. W. J. McNeil, M. W. Poole and G. R. M. Robb, [Phys. Rev. ST Accel. Beams **6**, 070701 (2003)](https://doi.org/10.1103/PhysRevSTAB.6.070701). The loading extended to unaveraged one-dimensional simulation, where it also models coherent spontaneous emission.

(ref-fawley-2003)=
**Spontaneous emission in the codes.** W. M. Fawley, [Nucl. Instrum. Methods A **507**, 19 (2003)](<https://doi.org/10.1016/S0168-9002(03)00827-1>). A transverse grid represents an artificially limited number of modes. GINGER reproduces the on-axis far-field intensity and the power within the central cone of the spontaneous emission, while the total spontaneous power depends on the bandpass and the mode content of the code.

(ref-huang-kim-2003)=
**Spontaneous and amplified radiation.** Z. Huang and K.-J. Kim, [Nucl. Instrum. Methods A **507**, 65 (2003)](<https://doi.org/10.1016/S0168-9002(03)00839-8>). The transition from spontaneous emission to exponential gain in the paraxial wave equation, compared with GINGER and GENESIS.

(ref-ssy-2000)=
**Diffraction in the SASE FEL.** E. L. Saldin, E. A. Schneidmiller and M. V. Yurkov, [Opt. Commun. **186**, 185 (2000)](<https://doi.org/10.1016/S0030-4018(00)01061-0>). The analytical theory of the linear regime: power, spectrum, angular distribution, correlation functions and transverse coherence.

(ref-litvinenko)=
**Macroparticle spontaneous radiation.** V. N. Litvinenko, [arXiv:1504.01430](https://arxiv.org/abs/1504.01430) (2015, presented at FEL 2002). [Macroparticles](glossary.md#gl-macroparticle) of charge far above the electron's radiate anomalously strongly. Paired clones of opposite charge separate the induced radiation from the spontaneous part.

(ref-pausch)=
**Form factors for macroparticles.** R. Pausch, A. Debus, A. Huebl, U. Schramm, K. Steiniger, R. Widera and M. Bussmann, [Nucl. Instrum. Methods A **909**, 419 (2018)](https://doi.org/10.1016/j.nima.2018.02.020), [arXiv:1802.03972](https://arxiv.org/abs/1802.03972). A point-like macroparticle radiates coherently with its whole weight. The transform of its shape function bounds the emission, and a form factor over the electrons it represents restores the incoherent part.

(ref-andriyash)=
**A spectral unaveraged algorithm.** I. A. Andriyash, R. Lehe and V. Malka, [J. Comput. Phys. **282**, 397 (2015)](https://doi.org/10.1016/j.jcp.2014.11.026), [arXiv:1405.1972](https://arxiv.org/abs/1405.1972). Notes that the spontaneous emission of a simulation with few macroparticles can differ greatly from the physical one.

(ref-hwang-qiang)=
**Artificial noise from mesh migration.** K. Hwang and J. Qiang, [Phys. Rev. Accel. Beams **21**, 120702 (2018)](https://doi.org/10.1103/PhysRevAccelBeams.21.120702). An improved shot-noise model under period averaging, suppressing the artificial noise that particle migration across the longitudinal mesh produces.

## Standards

**openPMD**, with the beamphysics extension for particle dumps and EXT_Wavefront for field dumps. The only dump format this code writes. Note: openPMD's `unitSI` is load-bearing, which is the opposite of the [statistics file](glossary.md#gl-statistics-file)'s convention, so the two never appear in one file. See [](fel-physics.md#sec-openpmd) and [](BMAD-STATS-SPEC.md).

**bmad-stats**, defined by this project in [](BMAD-STATS-SPEC.md) with the FEL extension in [](BMAD-STATS-EXT-FEL.md). The statistics file is deliberately not FEL-specific.

## On the citations in the source

Genesis4 and Bmad are both GPLv3. Every Fortran routine transcribed from Genesis4 cites the Genesis4 routine or file it came from in its header, and those citations are kept deliberately. They record provenance and let a reader check a transcription against its original, and they are the license's own trail. The manual carries the same citations in the collapsed Genesis4 source blocks of [](fel-physics.md), each linked to the file at the reference release.
