# FEL examples

Self-contained runs of the FEL tracker: one command each, no Genesis, no dump files
(one exception: `saturation_demo/`, the three-way comparison, needs the genesis4 binary). For
the validation benchmark against Genesis see `lucifer/tests/`. That is where the physics
is proven. The physics itself (equations, conventions, provenance) is the
manual, `lucifer/doc/fel-physics.tex`.

One directory per example, each holding a `run.nml`. Shared pieces live at this level:
`aramis.bmad` (the Benchmark1-SASE line: 6 FODO cells of 3.99 m helical undulator
segments, aw = 0.85 rms, 15 mm period, 5.8 GeV, resonant at 1 Angstrom) and
`plot_fel.py`. Every example runs the same way, and its outputs land in its own
directory:

```
cd lucifer/examples/<example>
../../../../debug/bin/lucifer run.nml        # or production/bin
python ../plot_fel.py <example>.stats.h5            # needs h5py + matplotlib; writes <example>.png
```

| Example | What it is | Time |
|---|---|---|
| `steady_state/` | Seeded single-slice gain curve | ~1 min |
| `taper/` | The same, with a two-stage undulator taper | ~1 min |
| `sase/` | Pure SASE: 96 slices, dark start, shot noise | ~90 s |
| `sase_wake/` | The SASE run plus resistive-wall/gap/roughness wakes | ~100 s |
| `import/` | A beam_init bunch resampled into slices (Genesis's importdistribution method), tracked dark | ~1 min |
| `bmad_wake/` | The SASE run with the chamber wake via Bmad's z_long machinery on every element | ~2 min |
| `unaveraged/` | One seeded segment with no period averaging: the real quiver, the coupling as an outcome, plus its averaged twin for the overlay | ~1 min |
| `crossed_undulator/` | The two-polarization afterburner: an x-planar set bunches, its quarter-turn twin radiates orthogonally from that bunching | ~30 s |
| `harmonics/` | Harmonic lasing (the field set): a dark third harmonic grows from the fundamental's bunching on a planar segment, with openPMD wavefront output | ~2 min |
| `chicane/` | A four-bend chicane between segments in absolute-time mode: half a wavelength of geometric delay flips the second segment between amplifying and absorbing | ~2 min |
| `saturation_demo/` | The one exception to "no Genesis": the full 57 m SASE case to saturation, three trackers (Genesis4 MPI, Bmad averaged, Bmad unaveraged) from identical dumps, one clock (its own `run.sh`, every input a real file in the directory) | ~25 min |

With no dump files named in the namelist, the program generates its own starting state
(quiet-start beam, and a Gaussian seed where `wavefront_init%seed_power > 0`). The beam is described
by Bmad's standard `beam_init` block (one bulk-bunch description shared with the
import path). The current is derived, never input: a Gaussian profile from
`bunch_charge` + `sig_z`, a flat bunch via Bmad's `distribution_type(3) = "GRID"`, and
`sig_z = 0` meaning steady state (the whole charge in one slice window). Any set
`beam_init` field the quiet start does not honor is refused by name. The header of
`lucifer.f90` documents every parameter and the honored-fields table. The undulator segments in the lattices
are real Bmad wiggler elements with `tracking_method = custom`, and their FEL parameters
live on the lattice. The attribute-to-parameter map and the parse-time refusals are the
manual's `sec:element`. There are no per-undulator namelist parameters.

`<example>.stats.h5` is the production statistics file (manual sec:stats): fixed Bmad
units, per-record per-slice arrays with beam datasets named as `bunch_params_struct`
components (plus `bunching`), `wavefront_params` for the field, and the evaluated
`calc_bunch_params` at element ends. The plot is ten panels against z: radiation power
and window field energy (log and linear), bunching, beam energy change and rms spread
(MeV, Bmad's convention), rms beam and field sizes
(the field sizes from the wavefront sigma(4,4) -- watch gain guiding), and beam
normalized emittances beside the field emittance sqrt(det) = M^2 lambda/4pi at element
ends. Two more panels appear when the run calls for them: a polarization split when the
run carries two live polarizations (field/y in the stats file), and per-harmonic power
when it carries harmonic fields (field/harm<h>). Note: Energies are eV, never gamma.
`<example>.diag.txt` remains the Genesis-comparison instrument
(same columns as always) and is what the benchmark harness reads.

The sections below cover the single-command examples. `crossed_undulator/`,
`harmonics/`, `chicane/` and `saturation_demo/` carry their own READMEs with their
measured tables.

## steady_state

A seeded, steady-state (single-slice) FEL: the benchmark line with a 3 kA slice and a
5 kW Gaussian seed at its waist. The quiet start is exact (initial bunching 4e-17, so
the FEL grows from the seed, not sampling noise). Measured: 7.5 m power gain length,
saturation at 1.6 GW around z = 38 m, and the beam sizes start on the matched values.

## taper

The same run over `taper/taper.bmad`: identical to `aramis.bmad` for the first four FODO
cells, but the last two cells' undulators are a second element definition, `UND2`, with
`b_max` (hence aw) 0.4% lower. This is what driving the FEL from lattice attributes
buys: a heterogeneous line is just different elements, with no per-segment program
input. Measured against `steady_state` (same seed, same start): the two gain curves are
bit-identical until the taper begins at z = 31.92 m. The untapered line saturates at
1.6 GW and falls back to 0.76 GW at z = 57 m as particles rotate in the bucket. The
step-down taper re-matches the resonance to the decelerated beam and the power still
climbs at the exit: 9.6 GW at z = 57 m, 12.7x the untapered exit power (5.9x its
saturation peak).

## import

A bunch described by Bmad's `beam_init_struct` (the native equivalent of Genesis's
`&beam`), generated, resampled into FEL slices by the transcribed Genesis
`importdistribution` method (`fel_import_mod`), and tracked dark through the full
line: SASE from an imported bunch (the resampling method is the manual's `sec:import`).
The bunch is a Gaussian test bunch sized to the FEL
window's economics (sigma_z = 1.2 nm, 30 fC, 3 kA peak). It is labeled a test bunch for
a reason: physical bunches are micrometers and need thousands of slices. The time
window derives from the bunch itself, so the diag file's per-slice current is the
Gaussian profile. The lattice is the whole optics specification. The reference energy
comes from its `e_tot`, and `init_beam_distribution` generates the bunch matched to the
Twiss in its beginning statement. Genesis's `match` transform is not ported, since a
Bmad lattice already says what it would say. Set `write_dist_file` to hand the
identical bunch to Genesis's `&importdistribution`, or `write_opmd_file` for
openPMD-beamphysics. `dist_file` reads openPMD back in place of `use_beam_init`.

## sase

Pure SASE with nothing external at all: the loader generates a 96-slice time window
(spacing 3 wavelengths), imposes physical shot noise (weighted Fawley loading), starts
the field dark, and the FEL grows from its own noise through the full line with slippage
active. The loader prints the per-slice electron count N_lambda, N_eff, and the quiet
floor it verified before imposing.

Measured on this input (seed 12345): startup power settles near 4 MW per slice after the
first segment, total power reaches 3.0 GW at z = 57 m with a per-slice spread of 0.91
(the SASE fluctuation), and the induced energy spread grows from 0.99 to 1.14 m_e c^2.
(These numbers moved when the transverse-map default became Bmad's own kernel: in a
dark start, a tiny transport change re-rolls the effective noise realization.) The
plot shows the physics directly. The total-power sawtooth is radiation slipping out of
the head of the finite window at each drift while fresh vacuum enters at the tail (a
real effect of any finite time window, identical in Genesis). Deep saturation of every
slice needs a window longer than the total slippage. The per-slice spaghetti in the
power panel is the slippage cascade itself.

## bmad_wake

The same chamber-wake physics as `sase_wake`, delivered through Bmad's own wake
machinery instead of the transcribed Genesis model (conventions: manual
`sec:seamwake`): `ztable.wake` is the Bane-Stupakov resistive-wall kernel for a 0.5 mm
copper chamber (exported by `write_wake_kernels`, sign-flipped to Bmad's
positive-decelerating convention, self-slice unhalved, causal side z < 0, padded past
the window), attached as an `sr_wake` `z_long` table to every element of the line. The
undulators apply it once per element at mid-element across the whole 96-slice
window, the quads and pipes through Bmad's own `track1_bunch`. (The drift slots are
pipes here: a Bmad drift cannot carry a wake.) Measured against the identical kernel
in the `wake_on` model: exit mean energy drop -2.324 vs -2.308 m_e c^2, interior
per-slice profiles agreeing to 0.7% -- one physical wake, two independent
implementations, two application granularities. Regenerate the table by running any
wake_on configuration with `write_wake_kernels = "kern.txt"` and applying the recipe
in `wake_lattice.bmad`'s header.

## unaveraged

The FEL with no period averaging (manual `sec:unaveraged`): one seeded steady-state
segment (`seg1.bmad`, 266 periods of the benchmark undulator, carrying
`fel_tracking = 1` as its own lattice attribute). The particles ride the real helical
field, quiver and all, at 20 integration substeps per period, with sin² entry/exit ramps
and the radiation as a co-evolving kick. Nothing in this path knows the coupling factor
`fc`. The energy exchange is just what the Lorentz force does. `run_averaged.nml` is the
identical configuration through the averaged default (a wrapper overrides the
element's attribute). The two gain curves are two independent formulations of the same
physics, and they agree to 7.7e-4 ln at the segment exit (measured here, closer to
Bmad's own kernel than the 3.8e-3 the transcribed maps gave). The harness pins the
coupling itself at ~6e-4 with dedicated probes, see the main README. Measured cost of
not averaging: 12 s vs 0.4 s for the averaged twin, ~30× on this config. The run also
writes `unaveraged.ledger.txt`: beam energy (relative to the reference) and window field
energy per record, whose sum is conserved (manual `eq:ledger`, checked at 1e-4 of the
turnover under the harness's strong-exchange probe). One segment sits in the lethargy
regime. The seed diffracts and the exponential growth is only beginning at the exit, so
the point of the plot is the overlay, not the gain.

## sase_wake

The same SASE run through a deliberately narrow copper chamber (0.5 mm radius), plus the
undulator gap wake and a 100 nm rough surface. The chamber is a tuned demonstration
case: at 5.8 GeV and 1 Angstrom a normal chamber's wake is small, and this exists to
make the physics visible. The run writes the per-slice energy-loss rate to
`sase_wake.wake.txt` (measured: 1.9 to 121 keV/m across the window, the head slices
losing least). The wake is causal, and the resistive-wall numerical impedance of Bane &
Stupakov sets the shape. Measured effect: the mean energy drops 8.3 m_e c^2 (about
4.2 MeV) over the 57 m line while the SASE still reaches 3.9 GW. The drop is clearly
visible in the energy panel against the 1 m_e c^2 initial spread.
