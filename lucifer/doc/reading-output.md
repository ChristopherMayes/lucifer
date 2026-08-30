# Reading the output

For someone handed a Lucifer output file who never runs the program. Every output is self-describing, so the tooling here hard-codes almost nothing and neither should yours.

| file | what it holds |
|---|---|
| `<out_root>.stats.h5` | The run: per-record beam and field statistics, the lattice table, the resolved inputs, provenance, and the axes to read them by. A `bmad-stats` 1.0 file, defined by [`BMAD-STATS-SPEC.md`](BMAD-STATS-SPEC.md) and [`BMAD-STATS-EXT-FEL.md`](BMAD-STATS-EXT-FEL.md) |
| `<out_root>-final.beam.h5`, `-at<n>-<ele>.beam.h5` | Particle dumps, openPMD-beamphysics |
| `<out_root>-final.wf.h5`, `-at<n>-<ele>.wf.h5` | Field dumps, openPMD EXT_Wavefront |
| `<out_root>-escaped.fld.h5`, `-pulse.fld.h5` | The escaped-field bank and the reconstructed full pulse, Genesis field conventions |
| `<out_root>.ledger.txt` | The unaveraged mode's energy ledger, one row per record |
| `<out_root>.diag.txt` | The Genesis-comparison text diagnostics, written only under `write_diag` |

## The Python side

`tests/scripts/read_stats.py` is the one reader everything in the tree uses. Because the file describes itself, it hard-codes no dataset names, so a file that grows a polarization or a harmonic needs no change in it.

```python
from read_stats import read_stats

with read_stats("run.stats.h5") as st:
    st.s                          # (nz,) path length, a variable on the record axis
    st.record                     # (nz,) the record axis itself
    st["beam/slice/current"]      # (nz, ns) amperes, as the file states them
    st.coord("plane")             # ('x', 'y', 'z')
    st.dim_coords("beam/slice/sigma")   # the coordinate of each dimension, in order
    st.units["beam/slice/t"]      # 's'
    st.run["p0c"]                 # what the run produced
    st.params["global"]["ran_seed"]     # what the user set
    st.ele_name                   # (nz,) element name per record, gathered through ix_ele
```

`tests/scripts/validate_bmad_stats.py` is the format's conformance checker. It knows no program and no extension, needs only numpy and h5py, and exits nonzero on any failure:

```
python3 validate_bmad_stats.py run.stats.h5
```

`examples/plot_fel.py` draws ten panels against `s` straight from a stats file, with an eleventh when a second polarization is live and a twelfth per harmonic:

```
python3 examples/plot_fel.py run.stats.h5 -o run.png
```

Units in these files are documentation, never a factor to apply. The values are already SI and eV. That is the opposite of openPMD's `unitSI`, which is why the two never appear in one file.

## The stats file describes itself

The production statistics live in `<out_root>.stats.h5`: a `bmad-stats` 1.0 file with
the `fel` extension (manual: the diagnostic-output section), the planned reset from the development format
`lucifer-stats` 2.x. The layout is deliberately not FEL-specific, and the file says
what it holds. Every dataset carries `@unit`, `@long_name`, `@description` and
`@axes`, and every name in `@axes` resolves to a `coords/` dataset, trailing label axes
included, so a reader needs no table of names and never infers a dimension from its
length: `scripts/read_stats.py` is the one reader everything in the tree uses, and it
hard-codes nothing. The acceptance test is `scripts/validate_bmad_stats.py`, the
format's own conformance checker, program-blind and run by the harness, which must
report zero failures. Units are fixed Bmad units (m, rad, eV, s, C, J, W) and the
attributes are documentation, never load-bearing, which is the opposite of openPMD's
`unitSI` and deliberately not mixed with it in one file. The version is refused by name
rather than negotiated.

| group | what |
|---|---|
| `coords/` | every axis once: `record` (the record axis) with `s`, `ix_ele` (`@indexes='ele'`) and `at_element_end` (`@selects='element_end'`) as variables on it, `element_end` and `s_element_end`, the `slice` axis with `ct_slice`, `t_slice` and `z_slice` on it, the `ele` axis, and the label axes `bmad`, `bmad_col`, `bmad_t`, `wavefront`, `wavefront_col`, `plane`, `mode` |
| `params/` | the input tree: one subgroup per honored input struct (`global`, `beam_init`, `bmad_com`, `space_charge_com`, `wake`, `sc`, `wavefront_init`, ...), each with `@struct`, every honored component resolved after defaults |
| `run/` | what the run produced: `p0c`, the species, `slice_spacing`, the axis lengths as bookkeeping cross-checks |
| `beam/slice/` | per-record sufficient statistics, `bunch_params_struct` names, `sigma` at natural rank (nz, ns, 6, 6), the envelope extremes `rel_max`/`rel_min` over `bmad_t`, plus derived `current` and `energy` in eV |
| `beam/slice_twiss/`, `beam/bunch/` | Bmad's own `calc_bunch_params`, per slice and whole window, on the element-end axis, the nine twiss quantities in `twiss/` over the `plane` axis and in `modes/` over the `mode` axis |
| `field/total/`, `field/x/`, `field/y/`, `field/harm<h>/` | one wavelength's total and each component, all with the same dataset names |
| `lattice/` | one row per element on the `ele` axis, for layout plots |
| `meta/` | which lattice and which input, as datasets. Not the lattice, and not a person |

Six rules earn their keep. **The record number is the axis**: `s` rides along as a
variable, because `s` repeats wherever two records land on one plane, which is what a
zero-length element applying a wake kick does, and a selection on a repeating index
answers silently wrong. It is `s`, Bmad's name for position along the lattice, and `z` is
left to mean the phase-space coordinate. **One record axis**: an element end is always
a record, so `at_element_end` is a boolean mask where a duplicated copy of every
element-end quantity used to be. That mask is information only the writer has, and no
reader can recover afterwards which record was the end. **Every axis has a coordinate**,
label axes included, and things that must not be confused are named apart: a square
matrix's two sides (`bmad`, `bmad_col`) so that selecting one entry needs no rule, and
the projected twiss planes (`plane`) from the normal modes (`mode`), because an
eigen-emittance is not a projected emittance and one axis carrying both invites an
average across them. The nine twiss quantities sit on those axes rather than in six
groups per set for two reasons. A group named `z` cannot sit beside a `z` coordinate,
which xarray and netCDF both refuse, and one array per quantity is 18 datasets where
six groups per set were 108. They are subgroups (`twiss/`, `modes/`) because one of the
nine is named `sigma`, as `bunch_params_struct` names it, and the covariance matrix
beside them is `sigma` too. **The slice axis exists**, and its
coordinates are exact: the slice number is the axis, with `ct_slice` and `t_slice` on it
free of any beta and `z_slice` marked as Bmad's z at the reference beta. The slice grid is
uniform in time, not in z, which is what makes slippage an exact integer shift of the
field record. All three carry `@head_direction`, publishing the migration invariant that
the high slice index is the window head, without which no per-slice profile can be
trusted, let alone overlaid on another code's. **No dataset's meaning depends on what
else the file holds**: `field/total/power` is the sum over live polarizations whether one
or two are live, `field/x/power` is always the x component, and `@components`,
`@harmonics` and `@derived_from` say which children are which so that a sum over them
cannot double-count. **Not computed is NaN**, not a zero that reads as an answer: the
theta moments away from element ends, an empty slice's moments, the twiss of a degenerate
slice.

The per-record beam datasets are sufficient statistics, so
`scripts/bunch_params_from_stats.py` reconstructs a bunch_params dict from any
(record, slice), and the harness holds that reconstruction against the twiss the
tracker stored from Bmad's own `calc_bunch_params`. The field side stores
`wavefront_params_struct` per slice: `centroid(4)` = (x, theta_x, y, theta_y),
`sigma(4,4)` Wigner moments, energy, power, on_axis_intensity, emit_x/y = sqrt(det)
(= M^2 lambda/4pi), with `angle_moments_valid` marking where the FFT-costed theta
rows were filled (element ends and bank time -- the `twiss_valid` pattern). Pulse
values are pooled downstream. The file stays raw.

`lattice/` is what a layout plot needs and what the file did without for too long: one
row per tracked element with element 0 included, on the `ele` axis that `coords/ix_ele`
indexes so a join is a gather, carrying `name`, `key`, `s_start`, `s_end`, `l`, `ds_step`, `is_fel`,
`fel_tracking`, `b_max`, `aw` as the physics used it, `l_period`, `ku`, `helical`,
`k1`, `tilt` and `z_offset`. Genesis writes per-step arrays. A table plus the existing
join says the same thing without a second copy of the record axis, and says what
Genesis cannot: signed quad strengths with a length, wake-carrying pipes, and the
tracking mode per element. It is not a lattice serialization, and nothing else in
the file is one either: see below. `dump_beam_at` / `dump_field_at`
dump openPMD files at named elements (Bmad locator syntax, unknown names
refused). `keep_escaped_field` banks every slice slippage transmits out of the window
(`-escaped.fld.h5`, with per-slice wavefront_params and z_transmit, the one place that
keeps Genesis field conventions since those two records have no home in the wavefront
extension) and reconstructs the full pulse at the exit plane (`-pulse.fld.h5`) by free-space propagation at
finalize, because transmitted light is fixed information and never re-interacts, so
whole-pulse statistics use the ABCD map on the banked moment matrices and never
propagate numerically.

`meta/` says which lattice and which input, and deliberately not more. Everything in it
is a dataset rather than an attribute, because HDF5 caps one attribute at 64 kB (the
largest that writes here is 65495 bytes, where a scalar string dataset took 3 MB) while
the echoed namelist is 12 kB and a real lattice text 37 kB, and the failure path was a
warning: provenance whose failure mode is silent absence must not sit on a resource
limit. The cost is that `meta/` no longer sits outside dataset-level identity
comparisons for free, since `input_echo` carries `out_root`, so the harness excludes it
by name instead. `meta/lattice_source` is the top-level
lattice file only, and says so, because Bmad's `call, file =` pulls in more and every
wrapper lattice here recorded a call statement while the lattice it called was absent.
`n_lattice_files` reports how many files the parser opened, so one means the text is the
whole story. Reproduction rests on the `lattice/` table and the input echo, not on this
group. And nothing here identifies a person by default: the timestamp and the Bmad
version identify the run, `lattice_file` is a base name, and the user and working
directory go in only under `global%record_environment`. Note that `input_echo` echoes
file names as the user typed them, so relative paths are the user's half of that.

Measured (check_diagnostics.py, in the harness -- cross-identities, not references):

| Check | Measured | Check level |
|---|---|---|
| bunch_params reconstruction vs stored calc_bunch_params | **4.6e-8** | 1e-6 |
| banked slice energies vs the ledger's U_escaped | 6.7e-16 | 1e-12 |
| analytic (moment-map) vs FFT-propagated rms at exit, 267 slices | **8.4e-3** | 2e-2 |
| pooled pulse sigma, analytic vs numerical routes | 3.4e-3 | 2e-2 |
| stats/escaped/pulse dataset-identical, 1 vs 8 threads | exact | (exactness, no level) |
| unknown dump element | refused by name | (refusal, no level) |

diag.txt is untouched (the Genesis-comparison instrument): every benchmark tier
reproduces bit for bit, including through the diag/stats fusion: the per-record
stats loop also evaluates the diag instrument (the identical fel_field_diag and
fel_slice_diag calls per slice, so each slice's arithmetic is unchanged), and the
diag writer only prints. That retired the formerly serial per-record diag sweeps
(all 96 field planes and 96 particle slices, every record), which were worth more
than the whole stats machinery costs: measured on the saturation demo's averaged run
(96 slices x 2048 particles, 255^2 grid, 12 threads), the pre-stats baseline was
36.5 s and the run with full statistics is ~32 s. The diagnostics deliverable made
the tracker 12% faster, net. Per-slice element-end twiss is evaluated through Bmad's
own calc_emittances_and_twiss_from_sigma_matrix fed from the already-computed
per-record moments (an element end always coincides with its last record), not by
re-summing particles.

The whole-window element-end row is assembled the same way, from the per-slice
moments, by the pooled-covariance identity, with each slice first moved from its local
z chart to the global window chart. That move is not a plain offset, since
`fel_concat_slices` places a particle at z_global = z_local + beta*(is-1)*spacing with
the particle's own beta, so within a slice z depends on pz. The identity, the shear
that carries the chart change, and why the between-group term makes the pool exact
rather than an average of covariances are in the manual's whole-window-row section.
What matters for reading the file is that the row is a pool of the slice moments and
never a re-sum of particles. This replaced a concatenation of every
particle in the window into one bunch plus Bmad's full 6D moments and Twiss on it, at
every element end -- 110 million particle visits on a 131-slice x 8192 case, all on
one thread. Measured against that particle sum when it landed: 4.0e-12 worst relative
over every whole-window quantity on the diagnostics config, 5.0e-11 over 48 element
ends of the 96-slice SASE example. Cost: the SASE example went 28.2 -> 26.7 s, and a
131-slice x 8192, 103-element case went 137.7 -> 126.9 s with core utilization
931% -> 1048% (the serial block it removed was 16.5% of that run's wall clock).
Across every configuration the harness runs, the element-end rows moved by 1e-15 to
1e-13 on most configs, with the moments themselves (centroid, sigma) worst at 2.7e-9
and a tail reaching 1.1e-7 confined to the normal modes' `eta`/`etap`, the dispersion
parameters that are ratios of near-cancelling small terms, so their relative agreement
says more about their conditioning than about either computation. The same caution
applies to whole beams: on a degenerate one (a window resampled far shorter than its
bunch, where the transverse moments nearly vanish) the z cross terms agree only to
4e-6 relative while agreeing to 1e-29 absolute -- physically the same number. That is
why the check measures on a physical configuration, and why it compares the whole
matrix rather than chasing the relative error of an individual near-zero entry.
The `element_end/` group is SELF-SUFFICIENT: beam moments and Twiss (whole window and
per slice) plus the radiation power, energy, on-axis intensity and bunching per slice.
That is what makes `comb_ds_save < 0` usable rather than a trap. Bmad's comb semantics
are kept verbatim ("< 0 => no comb calculated"), so that mode writes no per-record
rows at all, and the element-end row is then the only record of the run, which is why
it carries the field as well as the beam. With records present those field values are
the record's own (copied, so the datasets agree exactly). With no records they are
evaluated at the element end by the same routines, angle moments not needed. Measured
on the 96-slice SASE example: the two paths agree bit-for-bit.

Measured cost of the modes on a 131-slice x 8192, 103-element case (12 threads):
`comb_ds_save = -1` 124 s and 9.0 MB (element ends only, no per-record rows),
a large positive value 133 s and 15.7 MB (a record at each element end, since element
ends always record when the comb is not negative). `= 0.1` 132 s and 35.4 MB. The
per-record sweep is not the cost driver -- 86 records and 338 records time the same.

Known scaling limit, named for the follow-on: the stats accumulate in memory and write
once (demo: 64 MB). A tens-of-thousands-of-slices hard-X-ray window wants chunked
incremental writes instead.

## Particle dumps: openPMD carries the weights, Genesis .par cannot

A Genesis `.par.h5` holds one current per slice. A writer sends
`c*sum(w)/slice_spacing` and a reader divides it back out uniformly, so a beam whose
particles carry different weights comes back uniform. Per-particle weights are this
port's day-one difference from Genesis, so that format cannot hold this code's state,
and the dump is openPMD (`.beam.h5`): Bmad's own `hdf5_write_beam`, where openPMD's
macro-charge IS the per-particle weight. Converting such a beam to a Genesis `.par.h5`
is refused by name by the converter, with the slice, the weight spread and the total
charge in the message.

The slice partition is `particlePatches`, the standard's own partition of a species
record: one patch per slice, in window order, and an empty slice is a patch of no
particles. So the patch count IS the window, and the file needs no attributes of this
code's invention to describe it. What openPMD has no place for comes from the deck
instead: the wavelength, the slice spacing (`lambda0` and `window_sample`) and the
beamlet size (`nbins`). Reading a dump without `lambda0` is refused by name rather than
defaulted, since a wrong wavelength rescales every phase in the run. `one4one` needs no
storage at all: the flag asserts that every macroparticle carries one electron, which is
what the weights say.

The reference phase is folded into the file's time coordinate. The chart splits a
particle's ponderomotive phase into a per-beam reference and a per-particle lag,
`theta_j = phi0 + ks z_j / beta_j`, and no dump format has anywhere to put `phi0`, so
every reader restarts it at zero. A dump therefore writes the lag the whole phase
implies, which makes the file's time `-theta_j/(ks c)` and a restart exact. This is not
bookkeeping: the beam's phase against the field's phase is what the next segment's gain
is made of, and without the fold a mid-line restart lands 2.1e-2 away on the
windowed-composition check. Genesis stores `theta` itself and its reader does the same
fold, and the converter maps a Genesis `theta` to the same time, so the two formats mean
the same thing.

Measured (check_beam_format.py, the harness's beam-format section):

| check | level |
|---|---|
| openPMD file write, read, write | dataset-identical |
| x, y, px, py, gamma, current through a restart, in Genesis's chart | exact |
| theta through a restart, absolute | exact |
| split weights stored per particle, and bit-identical on read | exact |
| per-slice counts restored with an empty slice in the window | identical |
| one patch per slice, empty ones included | identical |
| four refusals (nonuniform weights to Genesis, not openPMD, patch count, no charge) | by name |

The chart, not the file, is what could cost a digit here: openPMD stores absolute
momenta and a time where this code keeps `px`, `py` and a lag, so those pass through
`P/p0` and `-beta*c*dt` on the way out and back. On this configuration every column
comes back exact, so the levels above are bounds on what the chart could cost rather
than expected values.
