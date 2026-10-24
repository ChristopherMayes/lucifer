# Lucifer

Lucifer (*lux ferre*, light-bringer) is a free-electron-laser tracker inside Bmad. Its physics is transcribed from Genesis 1.3 Version 4 (Genesis4), embedded in Bmad's lattice machinery, and validated against Genesis4 over its `benchmark/Benchmark1-SASE` configuration from bitwise-identical starting states. The documentation is published at <https://christophermayes.github.io/lucifer/>.

FEL segments are real Bmad wiggler elements tracked by a Bmad FEL method, so a Lucifer lattice is a Bmad lattice: quadrupoles, chicanes, phase shifters, wakes and apertures are Bmad's, and the FEL parameters are read from lattice attributes rather than a parallel namelist. Every other element tracks each slice's bunch with Bmad's own `track1_bunch`.

Lucifer is under development and unreleased. It has no version number, and its inputs, outputs and interfaces may change without notice.

## Layout

| directory | what it holds |
|---|---|
| [`src/`](src) | The library: the FEL modules, the wavefront propagator and, under `src/device/`, the Metal backend and its stub |
| [`programs/`](programs) | The main programs: `lucifer`, the smoke test, the advance benchmark and the wavefront drift test |
| [`docs/`](docs) | The documentation, a MyST site that also reads as Markdown |
| [`examples/`](examples) | Runnable cases, each a directory of real input files with its own README |
| [`tests/`](tests) | The keystone suite, the benchmark harness, its check scripts and, under `tests/wavefront/`, the wavefront validation |
| [`util/`](util) | The build script, the release check and the locations rule the scripts share |
| [`cmake/`](cmake) | `FindBmad.cmake`, which finds the Bmad tree to build against |

`BMAD_RELEASE` names the Bmad release Lucifer needs, a minimum. `searchf.namelist` is the index Bmad's `getf` reads.

## Building

Lucifer is a Bmad user program and builds with its own CMake against an installed Bmad, named by `LUCIFER_BMAD_ROOT`. The first route is conda-forge's `bmad` package, which `environment.yml` at the root names with the compilers and cmake:

```
conda env create -n lucifer-build -f environment.yml
conda activate lucifer-build
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DLUCIFER_BMAD_ROOT=$CONDA_PREFIX
cmake --build build-release --parallel
```

The binaries land in `build-release/bin`. `LUCIFER_BMAD_ROOT` can as well name a distribution build's `production/` or `debug/` root, and `util/conda_compile` builds that distribution from a checkout `BMAD_DIR` names and then Lucifer against it, the route the recorded digits come from. The Bmad is a release at or after the one in [`BMAD_RELEASE`](BMAD_RELEASE), and the configure refuses an older one. The whole of it is in [`docs/user-guide.md`](docs/user-guide.md#building).

## Documentation

| document | what it answers |
|---|---|
| [`docs/user-guide.md`](docs/user-guide.md) | How do I build it, describe a run, and run it? |
| [`docs/input-reference.md`](docs/input-reference.md) | What does this namelist parameter do, and what refuses it? |
| [`docs/reading-output.md`](docs/reading-output.md) | I have an output file. What is in it and how do I read it? |
| [`docs/validation.md`](docs/validation.md) | What is checked, how, and at what measured level? |
| [`docs/performance.md`](docs/performance.md) | Where does a run spend its time, on what machine? |
| [`docs/fel-physics.md`](docs/fel-physics.md) | What does it compute, and why is that right? |
| [`docs/BMAD-STATS-SPEC.md`](docs/BMAD-STATS-SPEC.md) | The statistics file format, normatively, with [`docs/BMAD-STATS-EXT-FEL.md`](docs/BMAD-STATS-EXT-FEL.md) |
| [`docs/startup-noise.md`](docs/startup-noise.md) | How much of an unseeded run's power is the FEL, measured on four machines |
| [`docs/changelog.md`](docs/changelog.md) | What changed, newest first |

The documents render as one site. With `mystmd` available, `myst build --html` in `docs/` produces it and `myst start` serves it locally.

## The two FEL methods

The FEL method is Bmad's own `fel_method` attribute, set on the element as any attribute is, and averaged and unaveraged elements mix freely in one line.

| `fel_method` | what it is for |
|---|---|
| `averaged` | The wiggle-averaged (KMR) model on Bmad's own kernel maps: the production workhorse. |
| `unaveraged` | Direct integration through the analytic undulator field, with no averaging and no resonance approximation. A production method whose ~30x cost buys full quiver dynamics, energy accounting the beam actually pays, polarization-agnostic coupling, and arbitrary harmonic content. Also an independent check on the averaged path, since the two share no approximation. |

The averaged method's transverse maps have one more option behind them, the transcribed Genesis4 maps, selected for a whole run by `global%transport_model`. It is validation-internal and no production run sets it.

## Running

A run is one input file naming a lattice:

```
lucifer lucifer.in
```

`lucifer.in` is the name the examples use by convention. It is an argument rather than a default, and any name works.

The examples are the fastest way in. Each is a directory of real inputs with its own README and the numbers measured on it:

```
cd examples/steady_state && ../../build-release/bin/lucifer lucifer.in
```

See [`docs/user-guide.md`](docs/user-guide.md) for the input file's three namelist groups, the lattice attributes, and what each output file is.

## Where it stands

Validated against Genesis4 on eleven tiers plus twenty-eight check sections, on both debug and production builds, before every commit. A moved digit is treated as a bug rather than a new baseline.

Nine of the eleven tiers are transcription checks and agree with Genesis4 at the floor set by its truncated impedance constants. The other two are priced model differences, not defects: the Bmad seam's interlude transport, and the unaveraged mode against Genesis4's averaged one. Every level, every attribution and the full tier table live in [`docs/validation.md`](docs/validation.md), which is their one home, so the numbers here are a click away rather than a copy that can drift.

What the program does, and what it does not do yet, are the two lists in [`docs/index.md`](docs/index.md). They are kept in one place because the examples are checked against them: every declared feature names the example that shows it, and the harness refuses a feature that has lost its example.

## License

Copyright 2026 the Lucifer authors.

Lucifer is licensed under the GNU General Public License version 3 only (SPDX: `GPL-3.0-only`). The license text is in [`LICENSE`](LICENSE).

The routines transcribed from Genesis 1.3 Version 4 are derivative works of it, modified from 2026-07-31 on, and are distributed under the same license. Genesis 1.3 Version 4 carries the notice

```
Copyright 2010, 2011 Paul Scherrer Institut. All rights reserved.
```

and is licensed under the GNU General Public License version 3. Every transcribed routine cites the Genesis4 file and lines it came from, and those citations are kept.
