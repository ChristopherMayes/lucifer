---
title: "Lucifer: an FEL tracker inside Bmad"
short_title: Introduction
---

Lucifer (*lux ferre*, light-bringer) tracks a free-electron laser inside Bmad. Its
physics is transcribed from Genesis 1.3 Version 4, embedded in Bmad's lattice
machinery, and validated against Genesis from bitwise-identical starting states.

An FEL segment is a real Bmad wiggler element with `tracking_method = custom`, so a
Lucifer lattice is a Bmad lattice. Quadrupoles, chicanes, phase shifters, wakes and
apertures are Bmad's own, the FEL parameters are lattice attributes rather than a
parallel namelist, and every non-FEL element tracks each slice's bunch through Bmad's
`track1_bunch`. A run is one input file naming a lattice.

## Three tracking methods

Selected per element by the `fel_tracking` lattice attribute, and they mix freely in
one line.

| `fel_tracking` | method | what it is for |
|---|---|---|
| unset or `0` | averaged | The default. The wiggle-averaged model on Bmad's own `bmad_standard` kernel maps. The production workhorse, and faster than Genesis at equal cores with richer in-run diagnostics. |
| `1` | unaveraged | Direct integration through the analytic undulator field, with no period averaging and no resonance approximation. A production method whose cost per step buys the full quiver dynamics, the energy accounting the beam actually pays, polarization-agnostic coupling and arbitrary harmonic content. It is also the tree's referee, sharing no approximation with the averaged path. |
| `-1` | transcribed Genesis | Genesis's transverse maps verbatim. Validation-internal: the comparison tiers select it through wrapper lattices, and no production lattice writes it. |

## Where to start

| If you want to | Read |
|---|---|
| build it, describe a run, and run it | [](user-guide.md) |
| know what a namelist parameter does | [](input-reference.md) |
| read an output file you were handed | [](reading-output.md) |
| know the physics and the conventions | [](fel-physics.md) |
| know what is checked and at what level | [](validation.md) |
| write or read the statistics format | [](BMAD-STATS-SPEC.md) |
| see what changed | [](changelog.md) |
| find a paper this code implements | [](references.md) |

Runnable cases live in `lucifer/examples/`, each a directory of real input files with
its own README and a runner.

## What is in, and what is not

In: the FEL element, per-particle weights throughout, OpenMP over slices with
bit-identical results at any thread count, physical shot noise under weights,
distribution import, slice migration, Genesis's collective effects and Bmad element
wakes across the whole bunch, the unaveraged mode, two polarizations, harmonic fields,
openPMD dumps in both directions, phasing between segments, spontaneous emission
honoring Bmad's global switches, and a self-describing statistics file.

Not in: simultaneous harmonic fields in the unaveraged mode, elliptical polarization
beyond the tilt-honored planar and helical limits, MPI (deliberate, since the
shared-memory design is measured faster at equal cores), one-to-one particle tracking,
undulator field errors, and GPU support.

Validation runs on eleven comparison tiers and thirteen check sections, on both debug
and production builds, before every commit. Nine of the eleven tiers are transcription
checks that agree with Genesis at the floor set by its truncated impedance constants.
The other two are priced model differences rather than defects. The levels, the
attributions and the tier table are in [](validation.md).
