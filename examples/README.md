# FEL examples

The examples are self-contained runs of the FEL tracker, one command each, with no Genesis 1.3 Version 4 (Genesis4) and no dump files. The one exception is `saturation_demo/`, the three-way comparison, which needs the genesis4 binary. Every directory holds its own README. It says what the example demonstrates and what to look for in the plot. The quantities measured on each run are in a table on the example's page, written by the runner.

The validation benchmark against Genesis4 is in `tests/`. The physics itself, with equations, conventions and provenance, is the manual, [`docs/fel-physics.md`](../docs/fel-physics.md).

## Running an example

    cd examples/<example>
    ../../build-release/bin/lucifer lucifer.toml        # or debug/bin
    python ../plot_fel.py <out_root>.stats.h5             # needs h5py + matplotlib

Every input file is named `lucifer.toml`, and variants carry a descriptive suffix: `lucifer_averaged.in`, `lucifer_detuned.in`. That is a convention for reading a directory. The program looks for no default name: the file name is an argument, and any name works.

An input file carries no comment beyond the one line that says why a setting departs from its default. Each line is a namelist parameter documented in [the input reference](../docs/input-reference.md). What an input file is for is said in this file and in the directory's own README, and an input file is short enough to paste into a terminal.

## Shared pieces

`aramis.bmad` is the Benchmark1-SASE line: 6 FODO cells of 3.99 m helical undulator segments, aw = 0.85 rms, 15 mm period, 5.8 GeV, resonant at 1 Angstrom. Most examples call it, and the ones that need a different line carry their own lattice.

`plot_fel.py` draws ten panels against s from any stats file, with an eleventh when a second polarization is live and a twelfth when the run carries harmonic fields. What the panels mean is [reading an output file](../docs/reading-output.md).

## The examples

Times are wall clock for a production build on twelve cores. A debug build is several times slower.

| Example | What it is | Time |
|---|---|---|
| [`steady_state/`](steady_state/) | Seeded single-slice gain curve, and mid-run openPMD dumps | about five seconds |
| [`frames/`](frames/) | A frame series along z over a slice range, with the field's reductions | under a second |
| [`taper/`](taper/) | The same line with a two-stage undulator taper, past saturation | about five seconds |
| [`sase/`](sase/) | Pure SASE: a time window of slices, dark start, physical shot noise, slippage | about half a minute |
| [`flash1/`](flash1/) | FLASH1 at DESY at 13.7 nm from its published parameters, against its measured gain length and pulse energy | about twenty seconds |
| [`lcls/`](lcls/) | LCLS at SLAC at 1.5 Angstrom from its published parameters, against its measured gain length and saturated power | about a minute and a half |
| [`sase_wake/`](sase_wake/) | The SASE run through a narrow chamber: resistive-wall, gap and roughness wakes | about half a minute |
| [`bmad_wake/`](bmad_wake/) | The same resistive-wall kernel through Bmad's own `z_long` machinery, against the transcribed model | about half a minute each |
| [`space_charge/`](space_charge/) | Long-range space charge on an attosecond bunch, and the short-range harmonics on a bunched beam | a few seconds each |
| [`migration/`](migration/) | Slice migration: 149k moves, the charge accounted, and what leaving it off costs | about fifteen seconds each |
| [`spontaneous/`](spontaneous/) | Undulator radiation under Bmad's own damping and fluctuation switches | about five seconds each |
| [`unaveraged/`](unaveraged/) | One segment with no period averaging, beside its averaged twin | about five seconds |
| [`mixed_line/`](mixed_line/) | One unaveraged segment among eleven averaged ones, in one line | a few seconds |
| [`import/`](import/) | A `beam_init` bunch loaded in the keep mode and in the sample mode, tracked dark, and the openPMD round trip | under a minute and under ten seconds |
| [`crossed_undulator/`](crossed_undulator/) | Two polarizations: an x-planar set bunches, its quarter-turn twin radiates orthogonally | about a second |
| [`harmonics/`](harmonics/) | Harmonic lasing: a dark third harmonic grows from the fundamental's bunching | about a second |
| [`chicane/`](chicane/) | A four-bend chicane between segments in absolute time: half a wavelength of geometry flips gain to absorption | about a second each |
| [`coherent_source/`](coherent_source/) | Tanaka's coherent-source retrieval against the per-particle deposit | a couple of seconds each |
| [`saturation_demo/`](saturation_demo/) | The full 57 m case to saturation, three trackers from identical dumps, one clock. Needs genesis4 | about twenty-five minutes |

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
| One load path | [`import/`](import/) |
| A coherent source model | [`coherent_source/`](coherent_source/) |

Thread independence is a few commands rather than a directory. The threads divide the slices, so it takes a time-dependent example. Add `run.write_diag = true` to `migration/lucifer.toml` and run it twice:

    cd migration
    for n in 1 12; do
      OMP_NUM_THREADS=$n ../../build-release/bin/lucifer lucifer.toml > t$n.log
      mv migration.diag.txt t$n.diag.txt
    done
    cmp t1.diag.txt t12.diag.txt

On this configuration, 1 thread against 12, the two diag files are identical. Compare the text file rather than the stats file. HDF5 object headers carry timestamps, so a whole-file compare would false-alarm on them. The benchmark runs the same check on the td1 tier at every commit.

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
| The escaped-field bank and the reconstructed pulse | A small-window case. On the `sase/` configuration those two files come to more than a gigabyte each |

## What a run writes

`<out_root>.stats.h5` is the statistics file, and [reading an output file](../docs/reading-output.md) covers it and the plot panels. `<out_root>.diag.txt` is the Genesis4-comparison text file, written on `run.write_diag = true` and large. Several examples write one more file named for what it carries: `.migration.txt`, `.wake.txt`, `.ledger.txt`, `.import.txt`.

Run outputs are not committed. `.gitignore` in this directory covers them by extension.
