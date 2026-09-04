# Lucifer: an FEL tracker, validated against Genesis 1.3 Version 4 (Genesis4)

Lucifer (*lux ferre*, light-bringer) is a free-electron-laser tracker inside Bmad. Its physics is transcribed from Genesis4 (GPL permits transcription), embedded in Bmad's lattice machinery by the seam, and validated against Genesis4 over its `benchmark/Benchmark1-SASE` configuration from bitwise-identical starting states.

FEL segments are real Bmad wiggler elements tracked by a Bmad FEL method, so a Lucifer lattice is a Bmad lattice: quadrupoles, chicanes, phase shifters, wakes and apertures are Bmad's, and the FEL parameters are read from lattice attributes rather than a parallel namelist. Every other element tracks each slice's bunch with Bmad's own `track1_bunch`.

## Documentation

| document | what it answers |
|---|---|
| [`doc/user-guide.md`](doc/user-guide.md) | How do I build it, describe a run, and run it? |
| [`doc/input-reference.md`](doc/input-reference.md) | What does this namelist parameter do, and what refuses it? |
| [`doc/reading-output.md`](doc/reading-output.md) | I have an output file. What is in it and how do I read it? |
| [`doc/validation.md`](doc/validation.md) | What is checked, how, and at what measured level? |
| [`doc/performance.md`](doc/performance.md) | Where does a run spend its time, on what machine? |
| [`doc/fel-physics.md`](doc/fel-physics.md) | What does it compute, and why is that right? |
| [`doc/BMAD-STATS-SPEC.md`](doc/BMAD-STATS-SPEC.md) | The statistics file format, normatively, with [`doc/BMAD-STATS-EXT-FEL.md`](doc/BMAD-STATS-EXT-FEL.md) |
| [`examples/`](examples) | Runnable cases, each a directory of real input files with its own README |
| [`doc/changelog.md`](doc/changelog.md) | What changed on this branch, newest first |

The documents render as one site. With `mystmd` available, `myst build --html` in
`doc/` produces it and `myst start` serves it locally. Each page also reads on its own
as Markdown.

## The two tracking methods

Bmad's own named methods, set on the element as any tracking method is, and they mix freely in one line.

| `tracking_method` | what it is for |
|---|---|
| `fel_averaged` | The wiggle-averaged (KMR) model on Bmad's own kernel maps: the production workhorse. |
| `fel_unaveraged` | Direct integration through the analytic undulator field, with no averaging and no resonance approximation. A production method whose ~30x cost buys full quiver dynamics, energy accounting the beam actually pays, polarization-agnostic coupling, and arbitrary harmonic content. Also an independent check on the averaged path, since the two share no approximation. |

The averaged method's transverse maps have one more option behind them, the transcribed Genesis4 maps, selected for a whole run by `global%transport_model`. It is validation-internal and no production run sets it.

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

Validated against Genesis4 on eleven tiers plus nineteen check sections, on both debug and production builds, before every commit. A moved digit is treated as a bug rather than a new baseline.

Nine of the eleven tiers are transcription checks and agree with Genesis4 at the floor set by its truncated impedance constants. The other two are priced model differences, not defects: the Bmad seam's interlude transport, and the unaveraged mode against Genesis4's averaged one. Every level, every attribution and the full tier table live in [`doc/validation.md`](doc/validation.md), which is their one home, so the numbers here are a click away rather than a copy that can drift.

What the program does, and what it does not do yet, are the two lists in
[`doc/index.md`](doc/index.md). They are kept in one place because the examples are
checked against them: every declared feature names the example that shows it, and the
harness refuses a feature that has lost its example.

## License

Lucifer is part of the Bmad distribution and is distributed on the distribution's terms.
Its physics is transcribed from Genesis 1.3 Version 4, which is licensed under the GNU
General Public License version 3, and every transcribed routine carries a citation to
the Genesis4 file and lines it came from. Those citations are kept deliberately.
