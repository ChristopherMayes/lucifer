# FEL examples

Self-contained runs of the FEL tracker, one command each, no Genesis 1.3 Version 4
(Genesis4) and no dump files. The one exception is `saturation_demo/`, the three-way
comparison, which needs the genesis4 binary. Every directory holds its own README with
what the example demonstrates, the numbers measured on it, and what to look for in the
plot.

For the validation benchmark against Genesis4 see `lucifer/tests/`. That is where the
physics is proven. The physics itself, with equations, conventions and provenance, is
the manual, [`lucifer/doc/fel-physics.md`](../doc/fel-physics.md).

## Running an example

    cd lucifer/examples/<example>
    ../../../production/bin/lucifer lucifer.in        # or debug/bin
    python ../plot_fel.py <out_root>.stats.h5             # needs h5py + matplotlib

Every input deck is named `lucifer.in`, and variants carry a descriptive suffix:
`lucifer_averaged.in`, `lucifer_detuned.in`. That is a convention for reading a
directory, not a default the program looks for. The file name is an argument, and any
name works.

The decks carry no comments. Each line is a namelist parameter documented in
[the input reference](../doc/input-reference.md), the story is this file and the
directory's own README, and a deck is short enough to paste into a terminal.

## Shared pieces

`aramis.bmad` is the Benchmark1-SASE line: 6 FODO cells of 3.99 m helical undulator
segments, aw = 0.85 rms, 15 mm period, 5.8 GeV, resonant at 1 Angstrom. Most examples
call it, and the ones that need a different line carry their own lattice.

`plot_fel.py` draws ten panels against s from any stats file, with an eleventh when a
second polarization is live and a twelfth per harmonic. What the panels mean is
[reading an output file](../doc/reading-output.md).

## The examples

Times are wall clock for a production build on twelve cores. A debug build is several
times slower.

| Example | What it is | Time |
|---|---|---|
| [`steady_state/`](steady_state/) | Seeded single-slice gain curve, and mid-run openPMD dumps | 5 s |
| [`taper/`](taper/) | The same line with a two-stage undulator taper, past saturation | 5 s |
| [`sase/`](sase/) | Pure SASE: 96 slices, dark start, physical shot noise, slippage | 25 s |
| [`sase_wake/`](sase_wake/) | The SASE run through a narrow chamber: resistive-wall, gap and roughness wakes | 25 s |
| [`bmad_wake/`](bmad_wake/) | The same resistive-wall kernel through Bmad's own `z_long` machinery, against the transcribed model | 25 s each |
| [`space_charge/`](space_charge/) | Long-range space charge on an attosecond bunch, and the short-range harmonics on a bunched beam | 7 s each |
| [`migration/`](migration/) | Slice migration: 149k moves, the charge accounted, and what leaving it off costs | 15 s each |
| [`spontaneous/`](spontaneous/) | Undulator radiation under Bmad's own damping and fluctuation switches | 5 s each |
| [`unaveraged/`](unaveraged/) | One segment with no period averaging, beside its averaged twin | 5 s |
| [`mixed_line/`](mixed_line/) | One unaveraged segment among eleven averaged ones, in one line | 7 s |
| [`import/`](import/) | A `beam_init` bunch resampled into slices and tracked dark, and the openPMD round trip | 10 s each |
| [`crossed_undulator/`](crossed_undulator/) | Two polarizations: an x-planar set bunches, its quarter-turn twin radiates orthogonally | 1 s |
| [`harmonics/`](harmonics/) | Harmonic lasing: a dark third harmonic grows from the fundamental's bunching | 1 s |
| [`chicane/`](chicane/) | A four-bend chicane between segments in absolute time: half a wavelength of geometry flips gain to absorption | 1 s each |
| [`coherent_source/`](coherent_source/) | Tanaka's coherent-source retrieval against the per-particle deposit | 2 s each |
| [`saturation_demo/`](saturation_demo/) | The full 57 m case to saturation, three trackers from identical dumps, one clock. Needs genesis4 | 25 min |

## What each feature is exemplified by

One row per feature declared in [the introduction](../doc/index.md).

| Feature | Shown by |
|---|---|
| The FEL element: parameters as lattice attributes | [`steady_state/`](steady_state/), and [`taper/`](taper/) for a heterogeneous line |
| Two tracking methods, Bmad's own named methods, mixable per element in one line | [`mixed_line/`](mixed_line/), [`unaveraged/`](unaveraged/) |
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
| A GPU backend for the averaged method | [`steady_state/`](steady_state/) and [`sase/`](sase/) run on it with `global%device = "metal"` and `grid_n_pts = 256`, on an Apple Silicon build (the device field solver takes powers of two, and refusing 255 it names 256). Measured levels and wall times are in the validation and performance pages |
| openPMD dumps in both directions | [`steady_state/`](steady_state/) writes mid-run, [`import/`](import/) round-trips a bunch through a file, [`saturation_demo/`](saturation_demo/) reads Genesis4's |
| A self-describing statistics file | Every example writes one |
| Distribution import | [`import/`](import/) |
| A coherent source model | [`coherent_source/`](coherent_source/) |

Thread independence is a few commands rather than a directory. The threads divide the
slices, so it takes a time-dependent example. Add `global%write_diag = T` to
`migration/lucifer.in` and run it twice:

    cd migration
    for n in 1 12; do
      OMP_NUM_THREADS=$n ../../../production/bin/lucifer lucifer.in > t$n.log
      mv migration.diag.txt t$n.diag.txt
    done
    cmp t1.diag.txt t12.diag.txt

Measured on this configuration: 96 slices, 1 thread against 12, the two diag files
identical at 28,960,122 bytes. Compare the text file rather than the stats file, since
HDF5 object headers carry timestamps and a whole-file compare would false-alarm on
them. The harness runs the same check on the td1 tier at every commit.

## Planned examples

Each of these is a feature with no example yet, and the reason it has none.

| Planned | Waiting on |
|---|---|
| Pitched and misaligned undulators | The feature |
| Simultaneous harmonics in the unaveraged mode | The feature |
| Elliptical polarization beyond the tilt-honored limits | The feature |
| Undulator field errors | The feature |
| One-to-one particle tracking | The feature |
| `chamber_wake%model = 'bmad'` | The model, which is a named follow-on |
| `space_charge%model = 'bmad_slice'` | The model, which is a named follow-on |
| The escaped-field bank and the reconstructed pulse | A small-window case. On the `sase/` configuration those two files come to 1.2 GB and 1.3 GB |

## What a run writes

`<out_root>.stats.h5` is the statistics file, and
[reading an output file](../doc/reading-output.md) covers it and the plot panels.
`<out_root>.diag.txt` is the Genesis4-comparison text file, written on
`global%write_diag = T` and large. Several examples write one more file named for what
it carries: `.migration.txt`, `.wake.txt`, `.ledger.txt`.

Run outputs are ignored by git, not committed. `.gitignore` in this directory covers
them by extension.
