---
title: References
short_title: References
---

What this code implements, and where each thing came from. Every entry names the source
and the routine that carries it, so a reader can go from a paper to the transcription
and back. The physics manual repeats these citations in place, section by section, with
the equations beside them.

## The codes this one is measured against

**Genesis 1.3 Version 4**, S. Reiche. The FEL physics of the averaged path is
transcribed from it, routine by routine, and the benchmark tiers compare against it from
shared dumps. Genesis is GPLv3, which is what permits the transcription, and the
per-routine citations in the source are kept for that reason. The original code is
described in S. Reiche, *GENESIS 1.3: a fully 3D time-dependent FEL simulation code*,
Nucl. Instr. and Meth. A **429** (1999) 243.

**Bmad**, D. Sagan. The lattice, the tracking of every non-FEL element, the wake
machinery, the beam I/O and the units and phase-space conventions are Bmad's.
D. Sagan, *Bmad: A relativistic charged particle simulation library*, Nucl. Instr. and
Meth. A **558** (2006) 356.

**MINERVA**, H. P. Freund and P. J. M. van der Slot. The production existence proof for
unaveraged FEL dynamics, and the source of the substeps-per-period envelope the
unaveraged mode is sized against. Lucifer differs by keeping the grid field and the
Lorentz force where MINERVA evaluates modal fields. See
[](fel-physics.md#sec-unaveraged).

**SIMPLEX**, T. Tanaka. Its source is unlicensed and unread. The coherent source is
implemented from the published paper alone, not from that code.

**Puffin** and **openPMD-beamphysics** are read-only references in this project's
working tree, used for conventions and for file conversion respectively.

## Papers the code implements

**Coherent source.** T. Tanaka, Phys. Rev. Accel. Beams **27**, 030703 (2024),
arXiv:2310.20197. The slice bunch factor splits into coherent and spatially incoherent
parts, and the incoherent part's radiation diffracts away. Implemented from the paper in
`fel_track_mod` (`global%source_model = 'coherent'`), with its guardrails refusing by
name. See [](fel-physics.md#sec-coherent-source).

**Resistive-wall wakes.** K. L. F. Bane and G. Stupakov, SLAC-PUB-10707. The numerical
impedance with AC (Drude) conductivity and Leontovich surface impedance, round and flat
geometry. Transcribed in `fel_collective_mod` (`fel_resistive_wall_wake`), with
Genesis's exact numerics kept. See [](fel-physics.md#sec-wakes).

**Shot-noise loading.** W. M. Fawley's algorithm for loading shot-noise microbunching in
FEL simulation codes, in the one-substitution weighted generalization that makes it
correct for per-particle weights. `fel_fawley_noise`, checked statistically against
$\langle|b(h)|^2\rangle = 1/N_\lambda$ under uniform and nonuniform weights. See
[](fel-physics.md#sec-noise).

**Incoherent (quantum) diffusion.** The standard undulator quantum-diffusion variance in
the Saldin form, as fitted in Genesis's `Incoherent.cpp`, with one draw per beamlet so
the quiet start's cancellation survives. Applied only when Bmad's own
`radiation_fluctuations` switch is on. See [](fel-physics.md#sec-eom).

**Normal-mode emittances.** A. Wolski's approach to general coupled linear optics
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

Genesis and Bmad are both GPLv3. The per-routine transcription citations in the Fortran
(`Field::track`, `Beam::track`, `Lattice::calcSlippage`, `SDDSBeam.cpp` and the rest) are
kept deliberately: they record provenance, they let a reader check a transcription
against its original, and they are the license's own trail.
