---
title: References
short_title: References
---

What this code implements, and where each thing came from. Every entry names the source
and the routine that carries it, so a reader can go from a paper to the transcription
and back. The physics manual repeats these citations in place, section by section, with
the equations beside them.

## The codes this one is measured against

**Genesis 1.3 Version 4 (Genesis4)**, S. Reiche. The FEL physics of the averaged path is
transcribed from it, routine by routine, and the benchmark tiers compare against it from
shared dumps. Genesis4 is GPLv3, which is what permits the transcription, and the
per-routine citations in the source are kept for that reason. The original code is
described in S. Reiche, *GENESIS 1.3: a fully 3D time-dependent FEL simulation code*,
[Nucl. Instr. and Meth. A **429** (1999) 243](https://doi.org/10.1016/S0168-9002(99)00114-X).

**Bmad**, D. Sagan. The lattice, the tracking of every non-FEL element, the wake
machinery, the beam I/O and the units and phase-space conventions are Bmad's.
D. Sagan, *Bmad: A relativistic charged particle simulation library*, [Nucl. Instr. and
Meth. A **558** (2006) 356](https://doi.org/10.1016/j.nima.2005.11.001).

**MINERVA**, H. P. Freund and P. J. M. van der Slot. The published existence proof that
unaveraged FEL dynamics is practical for production, and the published source of the
substeps-per-period envelope this mode's own convergence measurement is read against.
Lucifer differs by keeping the grid field and the Lorentz force where MINERVA evaluates
modal fields. Only the published work informed this port: MINERVA's source was not
read, and no run of it was compared against. See [](fel-physics.md#sec-unaveraged).

**SIMPLEX**, T. Tanaka. Its source is unlicensed and unread. The coherent source is
implemented from the published paper alone, not from that code.

**openPMD-beamphysics** is a library this code calls: the harness converts the Genesis4
reference dumps through it, and the wavefront round trip goes through its own openPMD
I/O in both directions.

## Papers the code implements

(ref-tanaka)=
**Coherent source.** T. Tanaka, [Phys. Rev. Accel. Beams **27**, 030703 (2024)](https://doi.org/10.1103/PhysRevAccelBeams.27.030703),
[arXiv:2310.20197](https://arxiv.org/abs/2310.20197). The slice bunch factor splits into coherent and spatially incoherent
parts, and the incoherent part's radiation diffracts away. Implemented from the paper in
`fel_track_mod` (`global%source_model = 'coherent'`), which refuses the configurations
the model does not cover. See [](fel-physics.md#sec-coherent-source).

(ref-bane-stupakov)=
**Resistive-wall wakes.** K. L. F. Bane and G. Stupakov, [SLAC-PUB-10707](https://www.slac.stanford.edu/pubs/slacpubs/10500/slac-pub-10707.pdf). The numerical
impedance with AC (Drude) conductivity and Leontovich surface impedance, round and flat
geometry. Transcribed in `fel_collective_mod` (`fel_resistive_wall_wake`), with
Genesis4's exact numerics kept. See [](fel-physics.md#sec-wakes).

(ref-fawley)=
**Shot-noise loading.** W. M. Fawley's algorithm for loading shot-noise microbunching in
FEL simulation codes, [Phys. Rev. ST Accel. Beams **5**, 070701 (2002)](https://doi.org/10.1103/PhysRevSTAB.5.070701), in the one-substitution weighted generalization that makes it
correct for per-particle weights. `fel_fawley_noise`, checked statistically against
$\langle|b(h)|^2\rangle = 1/N_\lambda$ under uniform and nonuniform weights. See
[](fel-physics.md#sec-noise).

**Incoherent (quantum) diffusion.** The standard undulator quantum-diffusion variance in
the Saldin form, as fitted in Genesis4's `Incoherent.cpp`, with one draw per beamlet so
the quiet start's cancellation survives. Applied only when Bmad's own
`radiation_fluctuations` switch is on. See [](fel-physics.md#sec-eom).

(ref-sase-theory)=
**SASE startup and gain estimates.** Ming Xie's fitting formula for the 3D power gain length, [Nucl. Instrum. Methods A **445**, 59 (2000)](https://doi.org/10.1016/S0168-9002(00)00114-5), and E. L. Saldin, E. A. Schneidmiller and M. V. Yurkov's effective shot-noise power, saturation length and efficiency, and the incoherent undulator power of a beam in its central cone, [New J. Phys. **12**, 035010 (2010)](https://doi.org/10.1088/1367-2630/12/3/035010), [arXiv:0912.4161](https://arxiv.org/abs/0912.4161), eqs. 1, 5 and 18. Computed from the input's parameters by the script of [](startup-noise.md) and placed beside the measurements there.

(ref-wolski)=
**Normal-mode emittances.** A. Wolski's approach to general coupled linear optics,
[Phys. Rev. ST Accel. Beams **9**, 024001 (2006)](https://doi.org/10.1103/PhysRevSTAB.9.024001),
supplies the eigen-emittances stored beside the projected Twiss planes. Bmad's
`calc_bunch_params` evaluates them, and
`tests/scripts/bunch_params_from_stats.py` re-derives them independently as a check. See
[](fel-physics.md#sec-stats).

## Standards

**openPMD**, with the **beamphysics** extension for particle dumps and **EXT_Wavefront**
for field dumps. The only dump format this code writes. Note that openPMD's `unitSI` is
load-bearing, which is the opposite of the statistics file's convention, so the two never
appear in one file. See [](fel-physics.md#sec-openpmd) and [](BMAD-STATS-SPEC.md).

**bmad-stats**, defined by this project in [](BMAD-STATS-SPEC.md) with the FEL extension
in [](BMAD-STATS-EXT-FEL.md). The statistics file is deliberately not FEL-specific.

## On the citations in the source

Genesis4 and Bmad are both GPLv3. The per-routine transcription citations in the Fortran
(`Field::track`, `Beam::track`, `Lattice::calcSlippage`, `SDDSBeam.cpp` and the rest) are
kept deliberately: they record provenance, they let a reader check a transcription
against its original, and they are the license's own trail.
