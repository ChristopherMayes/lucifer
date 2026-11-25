# FEL examples

The examples are self-contained runs of the FEL tracker, one command each, with no Genesis 1.3 Version 4 (Genesis4) and no dump files. The one exception is `saturation_demo/`, the three-way comparison, which needs the genesis4 binary. Every directory holds its own README. It says what the example demonstrates and what to look for in the plot. The quantities measured on each run are in a table on the example's page, written by the runner. `lcls_electrons/` is run by hand, since it takes over an hour, and its README states what its run measured.

The validation suite, which compares against Genesis4, is in `tests/`. The physics itself, with equations, conventions and provenance, is the manual, [`docs/fel-physics.md`](../docs/fel-physics.md).

## Running an example

    cd examples/<example>
    ../../build-release/bin/lucifer lucifer.toml        # or debug/bin
    python -m lucifer.plot <out_root>.stats.h5            # needs pip install -e "python[plot]"

Every input file is named `lucifer.toml`, and variants carry a descriptive suffix: `lucifer_detuned.toml`, `lucifer_transcribed.toml`. That is a convention for reading a directory. The program looks for no default name: the file name is an argument, and any name works.

An input file carries no comment beyond the one line that says why a setting departs from its default. Each line is an input key documented in [the input reference](../docs/input-reference.md). What an input file is for is said in this file and in the directory's own README, and an input file is short enough to paste into a terminal.

## Shared pieces

`aramis.bmad` is the Aramis line, which most examples call and which the next section states. The examples that need a different line carry their own lattice.

`python -m lucifer.plot` draws the summary figure of any stats file, ten panels against s with the lattice layout under them, an eleventh when a second polarization is live and a twelfth when the run carries harmonic fields. Each panel is one call of `lucifer.plot.plot`, which [the Python page](../docs/python.md) describes, and what the panels mean is [reading an output file](../docs/reading-output.md).

## The Aramis line

`aramis.bmad` is the Bmad translation of `Aramis.lat`, the lattice of the Benchmark1-SASE case that is distributed with Genesis 1.3 Version 4 (Genesis4). The manual calls it the Aramis line, after that file. Its parameters are those of the Benchmark1-SASE case and are not taken from a publication of a facility.

| the line | |
|---|---|
| undulator | helical, period 15 mm, rms undulator parameter 0.84853 |
| undulator segment | 3.99 m, 266 periods |
| segments | 12, two in each of 6 FODO cells of 9.5 m |
| after each segment | a drift of 0.44 m, a quadrupole of 0.08 m with k1 = 2 or -2 m^-2, a drift of 0.24 m |
| length | 57 m |
| reference energy | 11357.82 electron rest energies, 5.80 GeV |
| matched Twiss functions at the entrance | beta 8.54 m and alpha -0.703 in x, beta 17.4 m and alpha 1.40 in y |
| resonant wavelength | 0.1 nm |

The examples send one beam through this line, and the page of an example states where its beam differs.

| the beam | |
|---|---|
| energy | 5.80 GeV |
| current | 3.0 kA |
| normalized emittance | 0.4 µm in both planes |
| rms relative energy spread | 8.8e-5, which is one unit of the Lorentz factor |
| rms size at the entrance | 17 µm in x and 25 µm in y, matched to the cell |
| FEL parameter of the one-dimensional theory | 5.1e-4 |
| gain length of the one-dimensional theory | 1.34 m |

A time-dependent example states its window and its current profile, and a steady-state example tracks one slice of this current. The page of each example that `run_examples.sh` runs opens with a table of its case, the lattice and the beam, which is written from the statistics file of the run.

## The examples

Times are wall clock for a production build on twelve cores. A debug build is several times slower.

| Example | What it is | Line and beam | Time |
|---|---|---|---|
| [`steady_state/`](steady_state/) | Seeded single-slice gain curve, and mid-run openPMD dumps | Aramis at 0.1 nm, 5.80 GeV, 3.0 kA, steady state | about five seconds |
| [`frames/`](frames/) | A frame series along z over a slice range, with the field's reductions | two segments of Aramis, 3.0 kA uniform, 8 slices | under a second |
| [`taper/`](taper/) | The same line with a two-stage undulator taper, past saturation | Aramis with its last four segments tapered, 3.0 kA, steady state | about five seconds |
| [`sase/`](sase/) | Pure SASE: a time window of slices, dark start, physical shot noise, slippage | Aramis, 3.0 kA uniform, 96 slices | about half a minute |
| [`flash1/`](flash1/) | FLASH1 at DESY at 13.7 nm from its published parameters, against its measured gain length and pulse energy | FLASH1 at 13.7 nm, 668 MeV, 2.5 kA at the peak of a Gaussian bunch | about twenty seconds |
| [`lcls/`](lcls/) | LCLS at SLAC at 1.5 Angstrom from its published parameters, against its measured gain length and saturated power | LCLS at 0.151 nm, 13.6 GeV, 3.0 kA uniform | about a minute and a half |
| [`sase_wake/`](sase_wake/) | The SASE run through a narrow chamber: resistive-wall, gap and roughness wakes | Aramis, the beam and the window of `sase` | about half a minute |
| [`bmad_wake/`](bmad_wake/) | The same resistive-wall kernel through Bmad's own `z_long` machinery, against the transcribed model | Aramis, the beam and the window of `sase` | about half a minute each |
| [`space_charge/`](space_charge/) | Long-range space charge on an attosecond bunch, and the short-range harmonics on a bunched beam | Aramis, a Gaussian bunch of 4 nm rms with 3.0 kA at its peak, and one steady-state run | a few seconds each |
| [`migration/`](migration/) | Slice migration: 149k moves, the charge accounted, and what leaving it off costs | Aramis, the beam and the window of `sase` | about fifteen seconds each |
| [`spontaneous/`](spontaneous/) | Undulator radiation under Bmad's own damping and fluctuation switches | Aramis, 3.0 kA, steady state | about five seconds each |
| [`unaveraged/`](unaveraged/) | One segment with no period averaging, beside its averaged twin | one segment of Aramis, 3.0 kA, steady state | about five seconds |
| [`mixed_line/`](mixed_line/) | One unaveraged segment among eleven averaged ones, in one line | Aramis, 3.0 kA, steady state | a few seconds |
| [`import/`](import/) | A `beam_init` bunch loaded in the keep mode and in the sample mode, tracked dark, and the openPMD round trip | Aramis, a Gaussian bunch of 1.2 nm rms with 3.0 kA at its peak | under a minute and under ten seconds |
| [`crossed_undulator/`](crossed_undulator/) | Two polarizations: an x-planar set bunches, its quarter-turn twin radiates orthogonally | two planar undulators of 0.6 m with the period and the rms parameter of Aramis, 5.80 GeV, 3.0 kA, steady state | about a second |
| [`harmonics/`](harmonics/) | Harmonic lasing: a dark third harmonic grows from the fundamental's bunching | one planar segment of 3.96 m with the period and the rms parameter of Aramis, 5.80 GeV, 3.0 kA, steady state | about a second |
| [`chicane/`](chicane/) | A four-bend chicane between segments in absolute time: half a wavelength of geometry flips gain to absorption | two helical segments of 0.99 m with the undulator of Aramis and a chicane between them, 5.80 GeV, 3.0 kA, steady state | about a second each |
| [`coherent_source/`](coherent_source/) | Tanaka's coherent-source retrieval against the per-particle deposit | Aramis, 3.0 kA, steady state | a couple of seconds each |
| [`saturation_demo/`](saturation_demo/) | The full 57 m case to saturation, three trackers from identical dumps, one clock. Needs genesis4 | Aramis, 3.0 kA uniform, 96 slices | about twenty-five minutes |
| [`lcls_electrons/`](lcls_electrons/) | The LCLS window of [`lcls/`](lcls/) with one particle for every electron, $1.9\times10^{8}$ of them. Needs 17 GB | LCLS, the beam and the window of `lcls` | over an hour on twelve cores |

## What each feature is exemplified by

The table has one row per feature declared in [the introduction](../docs/index.md).

| Feature | Shown by |
|---|---|
| The FEL element: parameters as lattice attributes | [`steady_state/`](steady_state/), and [`taper/`](taper/) for a heterogeneous line |
| Two FEL methods, Bmad's own fel_method attribute, mixable per element in one line | [`mixed_line/`](mixed_line/), [`unaveraged/`](unaveraged/) |
| Per-particle weights throughout | [`sase/`](sase/) for weighted shot noise, [`migration/`](migration/) for weighted migration |
| Time dependence with an exact integer slippage shift | [`sase/`](sase/) |
| Slice migration | [`migration/`](migration/) |
| Collective effects: chamber wakes | [`sase_wake/`](sase_wake/) |
| Collective effects: Bmad element wakes across the window | [`bmad_wake/`](bmad_wake/) |
| Collective effects: space charge, short and long range | [`space_charge/`](space_charge/) |
| Two polarizations, tilt honored | [`crossed_undulator/`](crossed_undulator/) |
| Harmonic field sets | [`harmonics/`](harmonics/) |
| Spontaneous emission on Bmad's switches | [`spontaneous/`](spontaneous/) |
| Undulator tapering | [`taper/`](taper/) |
| Phasing between segments | [`chicane/`](chicane/) |
| Shared-memory parallelism, bit-identical at any thread count | Demonstrated inline, below |
| A GPU backend for the averaged method | [`steady_state/`](steady_state/), [`sase/`](sase/), [`harmonics/`](harmonics/), [`crossed_undulator/`](crossed_undulator/) and [`migration/`](migration/) run on it with `run.device = "metal"` and `grid_n_pts = 256`, on an Apple Silicon build (the device field solver takes powers of two, and refusing 255 it names 256). Measured levels and wall times are in the validation and performance pages |
| openPMD dumps in both directions | [`steady_state/`](steady_state/) writes mid-run, [`import/`](import/) round-trips a bunch through a file, [`frames/`](frames/) writes a series along z, [`saturation_demo/`](saturation_demo/) reads Genesis4's |
| A self-describing statistics file | Every example writes one |
| One load path | [`import/`](import/), and [`lcls_electrons/`](lcls_electrons/) for the electrons mode |
| A coherent source model | [`coherent_source/`](coherent_source/) |

Thread independence is a few commands rather than a directory. The threads divide the slices, so it takes a time-dependent example. Add `run.write_diag = true` to `migration/lucifer.toml` and run it twice:

    cd migration
    for n in 1 12; do
      OMP_NUM_THREADS=$n ../../build-release/bin/lucifer lucifer.toml > t$n.log
      mv migration.diag.txt t$n.diag.txt
    done
    cmp t1.diag.txt t12.diag.txt

On this configuration, 1 thread against 12, the two diag files are identical. Compare the text file rather than the stats file. HDF5 object headers carry timestamps, so a whole-file compare would false-alarm on them. The validation suite runs the same check on the td1 tier at every commit.

## Planned examples

Each of these is a feature with no example yet, and the reason it has none.

| Planned | Waiting on |
|---|---|
| Pitched and misaligned undulators | The feature |
| Simultaneous harmonics in the unaveraged mode | The feature |
| Elliptical polarization beyond the tilt-honored limits | The feature |
| Undulator field errors | The feature |
| One-to-one particle tracking | The feature |
| `wake.model = 'bmad'` | The model, which is a named follow-on |
| `space_charge.model = 'bmad_slice'` | The model, which is a named follow-on |

## What a run writes

`<out_root>.stats.h5` is the statistics file, and [reading an output file](../docs/reading-output.md) covers it and the plot panels. `<out_root>.diag.txt` is the Genesis4-comparison text file, written on `run.write_diag = true` and large. Several examples write one more file named for what it carries: `.migration.txt`, `.wake.txt`, `.ledger.txt`, `.import.txt`. A run from which light slips out of the window writes `<out_root>-escaped.wf.h5`, that light, and `<out_root>-pulse.wf.h5`, the whole pulse at the exit, both openPMD wavefront files, unless `run.keep_escaped_field = false`.

Run outputs are not committed. `.gitignore` in this directory covers them by extension.
