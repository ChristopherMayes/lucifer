# Lucifer: an FEL tracker, validated against Genesis 1.3 Version 4 (Genesis4)

Lucifer (*lux ferre*, light-bringer) is a free-electron-laser tracker inside Bmad. Its physics is transcribed from Genesis4 (GPL permits transcription), embedded in Bmad's lattice machinery by the seam, and validated against Genesis4 over its `benchmark/Benchmark1-SASE` configuration from bitwise-identical starting states.

FEL segments are real Bmad wiggler elements with `tracking_method = custom`, so a Lucifer lattice is a Bmad lattice: quadrupoles, chicanes, phase shifters, wakes and apertures are Bmad's, and the FEL parameters are read from lattice attributes rather than a parallel namelist. Every other element tracks each slice's bunch with Bmad's own `track1_bunch`.

## Documentation

| document | what it answers |
|---|---|
| [`doc/user-guide.md`](doc/user-guide.md) | How do I build it, describe a run, and run it? |
| [`doc/input-reference.md`](doc/input-reference.md) | What does this namelist parameter do, and what refuses it? |
| [`doc/reading-output.md`](doc/reading-output.md) | I have an output file. What is in it and how do I read it? |
| [`doc/validation.md`](doc/validation.md) | What is checked, how, and at what measured level? |
| [`doc/fel-physics.md`](doc/fel-physics.md) | What does it compute, and why is that right? |
| [`doc/BMAD-STATS-SPEC.md`](doc/BMAD-STATS-SPEC.md) | The statistics file format, normatively, with [`doc/BMAD-STATS-EXT-FEL.md`](doc/BMAD-STATS-EXT-FEL.md) |
| [`examples/`](examples) | Runnable cases, each a directory of real input files with its own README |
| [`doc/changelog.md`](doc/changelog.md) | What changed on this branch, newest first |

The documents render as one site. With `mystmd` available, `myst build --html` in
`doc/` produces it and `myst start` serves it locally. Each page also reads on its own
as Markdown.

## The three tracking methods

Chosen per element with the `fel_tracking` lattice attribute, and they mix freely in one line.

| `fel_tracking` | method | what it is for |
|---|---|---|
| unset or `0` | averaged | The default. The wiggle-averaged (KMR) model on Bmad's own kernel maps: the production workhorse. |
| `1` | unaveraged | Direct integration through the analytic undulator field, with no averaging and no resonance approximation. A production method whose ~30x cost buys full quiver dynamics, energy accounting the beam actually pays, polarization-agnostic coupling, and arbitrary harmonic content. Also the tree's referee, sharing no approximation with the averaged path. |
| `-1` | transcribed Genesis4 | Genesis4's transverse maps verbatim. Validation-internal: the comparison tiers select it through wrapper lattices, and no production lattice writes it. |

## Building and running

Lucifer builds as part of Bmad's normal distribution build, and the binary lands in
`production/bin/lucifer` (or `debug/bin/lucifer` for a debug build) alongside every
other Bmad program.

A run is one input file naming a lattice:

```
lucifer lucifer.in
```

`lucifer.in` is the name the examples use by convention. It is an argument rather than
a default, and any name works.

The examples are the fastest way in. Each is a directory of real inputs with its own
README and the numbers measured on it:

```
cd examples/steady_state && lucifer lucifer.in
```

See [`doc/user-guide.md`](doc/user-guide.md) for the input file's three namelist groups, the lattice attributes, and what each output file is.

## Where it stands

Validated against Genesis4 on eleven tiers plus eighteen check sections, on both debug and production builds, before every commit. A moved digit is treated as a bug rather than a new baseline.

Nine of the eleven tiers are transcription checks and agree with Genesis4 at the floor set by its truncated impedance constants. The other two are priced model differences, not defects: the Bmad seam's interlude transport, and the unaveraged mode against Genesis4's averaged one. Every level, every attribution and the full tier table live in [`doc/validation.md`](doc/validation.md), which is their one home, so the numbers here are a click away rather than a copy that can drift.

In: the FEL element, per-particle weights throughout, OpenMP over slices with bit-identical results at any thread count, physical shot noise under weights, distribution import, slice migration, Genesis4's collective effects and Bmad element wakes across the whole bunch, the unaveraged mode, two polarizations, harmonic fields, openPMD dumps in both directions, phasing between segments, spontaneous emission honoring Bmad's global switches, and a self-describing statistics file.

Not in: simultaneous harmonic fields in the unaveraged mode, elliptical polarization beyond the tilt-honored planar and helical limits, MPI (deliberate: the shared-memory design is measured faster at equal cores), one-to-one particle tracking, undulator field errors, and GPU support.
