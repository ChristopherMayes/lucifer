---
title: "Lucifer: an FEL tracker inside Bmad"
short_title: Introduction
---

Lucifer (*lux ferre*, light-bringer) tracks a free-electron laser inside Bmad. Its
physics is transcribed from Genesis 1.3 Version 4 (Genesis4), embedded in Bmad's lattice
machinery, and validated against Genesis4 from bitwise-identical starting states.

An FEL segment is a real Bmad wiggler element tracked by a Bmad tracking method, so a
Lucifer lattice is a Bmad lattice. Quadrupoles, chicanes, phase shifters, wakes and
apertures are Bmad's own, the FEL parameters are lattice attributes rather than a
parallel namelist, and every non-FEL element tracks each slice's bunch through Bmad's
`track1_bunch`. A run is one input file naming a lattice.

## Features

- **The FEL element.** An FEL segment is a Bmad wiggler, so its parameters are lattice attributes and one lattice serves tracking, optics and layout.
- **Two tracking methods**, averaged and unaveraged, set per element like any Bmad tracking method and mixing freely in one line. The [next section](#intro-two-methods) says what each is for.
- **Per-particle weights** throughout: physical shot noise, collective effects and diagnostics are all weight-correct.
- **Time dependence** with an exact integer slippage shift, and slice migration when a particle's ponderomotive phase leaves its slice.
- **Collective effects**: resistive-wall, geometric and roughness wakes, short-range and long-range space charge, and Bmad element wakes applied across the whole time window.
- **Two polarizations** and harmonic field sets, with tilt honored on planar and helical undulators.
- **Undulator tapering**, as a line whose segments are different elements rather than a per-segment input.
- **Phasing between segments**, following the real geometry in absolute time tracking and autophased in the default relative mode.
- **Spontaneous emission** gated on Bmad's own `radiation_damping` and `radiation_fluctuations` switches.
- **Shared-memory parallelism** over slices, bit-identical at any thread count.
- **A GPU backend** for the averaged method on Apple Silicon (`global%device = "metal"`): the beam and the field stay resident on the device through each FEL element, measured between 7x and 12x against twelve CPU cores and judged per step against the FP64 path ([](validation.md#val-device)), harmonic field sets, both polarizations and slice migration included. Configurations its kernels do not cover are refused, never run quietly on the CPU instead.
- **openPMD** particle and field dumps in both directions, and a self-describing statistics file described by [](BMAD-STATS-SPEC.md).
- **Distribution import**, resampling a `bunch_struct` into FEL slices.
- **A coherent source model**, Tanaka's retrieval, for transversely coherent beams.

(intro-two-methods)=
## The two tracking methods

Bmad's own named methods, set on the element as any tracking method is.

| `tracking_method` | what it is for |
|---|---|
| `fel_averaged` | The wiggle-averaged (KMR) model on Bmad's own `bmad_standard` kernel maps. The production method. |
| `fel_unaveraged` | Direct integration through the analytic undulator field, with no period averaging and no resonance approximation. Its cost per step buys the full quiver dynamics, the energy accounting the beam actually pays, polarization-agnostic coupling and arbitrary harmonic content. It is also an independent check on the averaged method, since the two share no approximation. |

The averaged method's transverse maps have a third option behind them, the transcribed
Genesis4 maps, which the comparison tiers select through [](input-reference.md)'s
`global%transport_model`. It is validation-internal and no production run sets it.

## Where to start

| If you want to | Read |
|---|---|
| build it, describe a run, and run it | [](user-guide.md) |
| know what a namelist parameter does | [](input-reference.md) |
| read an output file you were handed | [](reading-output.md) |
| know the physics and the conventions | [](fel-physics.md) |
| know what is checked and at what level | [](validation.md) |
| know where a run spends its time | [](performance.md) |
| run something and see a figure | [](generated/examples/examples.md) |
| write or read the statistics format | [](BMAD-STATS-SPEC.md) |
| see what changed | [](changelog.md) |
| find a paper this code implements | [](references.md) |

Runnable cases live in `lucifer/examples/`, each a directory of real input files with
its own README and the numbers measured on it. Every one of them has a page here, with
its figure and its input files: [](generated/examples/examples.md) carries the table of
which example shows which feature.

## Known missing features

- Simultaneous harmonic fields in the unaveraged mode.
- Elliptical polarization beyond the tilt-honored planar and helical limits.
- One-to-one particle tracking, and the particle sorting it implies.
- Undulator field errors.
- On the GPU: collective effects and the unaveraged mode, each refused. A CUDA backend, which the device interface is shaped for.
- MPI. This one is a design decision rather than a gap: the shared-memory design is what this program is, and the case for adding MPI would have to be made on measurement.

Validation runs on eleven comparison tiers and twenty-five check sections, on both
debug and production builds, before every commit. Seven of the eleven tiers are
transcription checks that agree with Genesis4 at the floor set by its truncated
impedance constants, three are priced model differences rather than defects, and one
proves nonuniform weights invisible against the uniform run. The levels, the
attributions and the tier table are in [](validation.md).
