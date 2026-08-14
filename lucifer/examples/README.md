# FEL examples

Self-contained runs of the FEL tracker: one command each, no Genesis, no dump files. For
the validation benchmark against Genesis — which is where the physics is proven — see
`bsim/fel/tests/`.

One directory per example, each holding a `run.nml`. Shared pieces live at this level:
`aramis.bmad` (the Benchmark1-SASE line: 6 FODO cells of 3.99 m helical undulator
segments, aw = 0.85 rms, 15 mm period, 5.8 GeV, resonant at 1 Angstrom) and
`plot_fel.py`. Every example runs the same way, and its outputs land in its own
directory:

```
cd bsim/fel/examples/<example>
../../../../debug/bin/fel_track_test run.nml        # or production/bin
python ../plot_fel.py <example>.diag.txt            # needs matplotlib; writes <example>.diag.png
```

| Example | What it is | Time |
|---|---|---|
| `steady_state/` | Seeded single-slice gain curve | ~1 min |
| `taper/` | The same, with a two-stage undulator taper | ~1 min |
| `sase/` | Pure SASE: 96 slices, dark start, shot noise | ~90 s |
| `sase_wake/` | The SASE run plus resistive-wall/gap/roughness wakes | ~100 s |
| `import/` | A beam_init bunch resampled into slices (Genesis's importdistribution method), tracked dark | ~1 min |

With no dump files named in the namelist, the program generates its own starting state
(quiet-start beam, and a Gaussian seed where `gen_power > 0`); the header of
`fel_track_test.f90` documents every parameter. The undulator segments in the lattices
are real Bmad wiggler elements with `tracking_method = custom`, and their FEL parameters
live on the lattice: aw derives from `b_max` and `l_period` (helical
aw = c·b_max/(k_u·m_e c²), rms convention; a planar device divides by √2), helicity from
`field_calc`. There are no per-undulator namelist parameters. A lattice whose FEL
element is missing `b_max` or `l_period`, or uses a fieldmap `field_calc`, is refused by
name at parse time.

`<example>.diag.txt` is the gain curve: one row per slice per integration step with
columns `z, slice, power, on_axis_intensity, bunching, bunching_phase, mean_gamma,
sigma_gamma, sigma_x, sigma_y`. The plot is four panels against z: radiation power
(log), bunching, energy change and rms spread, and the transverse rms sizes showing the
FODO betatron oscillation. On time-dependent files each thin gray line is a slice, the
bold line the total power or slice average, and the slippage echelon is directly visible
in the per-slice power. `<example>-final.fld.h5` / `-final.par.h5` are the end state in
Genesis dump format (readable by `openPMD-beamphysics`, h5py, or Genesis itself); add
`write_initial = T` to also dump the generated start — useful for handing the identical
initial condition to Genesis (`&importbeam` / `&importfield`), which is how the loader
was validated (Genesis imported the generated dumps and agreed at 1.5e-5, the constants
floor of every Genesis comparison here).

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
bit-identical until the taper begins at z = 31.92 m; the untapered line saturates at
1.6 GW and falls back to 0.76 GW at z = 57 m as particles rotate in the bucket, while
the step-down taper re-matches the resonance to the decelerated beam and the power still
climbs at the exit — 9.6 GW at z = 57 m, 12.7x the untapered exit power (6x its
saturation peak).

## import

A bunch described by Bmad's `beam_init_struct` -- the native equivalent of Genesis's
`&beam` -- generated, resampled into FEL slices by the transcribed Genesis
`importdistribution` method (`fel_import_mod`), and tracked dark through the full
line: SASE from an imported bunch. The bunch is a Gaussian TEST bunch sized to the FEL
window's economics (sigma_z = 1.2 nm, 30 fC, 3 kA peak -- honestly labeled: physical
bunches are micrometers and need thousands of slices). The time window derives from
the bunch itself, so the diag file's per-slice current is the Gaussian profile; the
lattice is the whole optics specification -- the reference energy comes from its
`e_tot` and `init_beam_distribution` generates the bunch matched to the Twiss in its
beginning statement (Genesis's `match` transform is not ported; a Bmad lattice already
says what it would say). Set `write_dist_file` to hand the identical bunch to
Genesis's `&importdistribution`, or `write_opmd_file` for openPMD-beamphysics;
`dist_file` reads openPMD back in place of `use_beam_init`.

## sase

Pure SASE with nothing external at all: the loader generates a 96-slice time window
(spacing 3 wavelengths), imposes physical shot noise (weighted Fawley loading — the
loader prints the per-slice electron count N_lambda, N_eff, and the quiet floor it
verified before imposing), starts the field dark, and the FEL grows from its own noise
through the full line with slippage active.

Measured on this input (seed 12345): startup power settles near 4 MW per slice after the
first segment, total power reaches 4.0 GW at z = 57 m with a per-slice spread of 0.79 —
the SASE fluctuation — and the induced energy spread grows from 1.0 to 1.15 m_e c^2. The
plot shows the physics directly: the total-power sawtooth is radiation slipping out of
the head of the finite window at each drift while fresh vacuum enters at the tail (a
real effect of any finite time window, identical in Genesis; deep saturation of every
slice needs a window longer than the total slippage), and the per-slice spaghetti in the
power panel is the slippage cascade itself.

## sase_wake

The same SASE run through a deliberately NARROW copper chamber — 0.5 mm radius, a tuned
demonstration case, honestly labeled: at 5.8 GeV and 1 Angstrom a normal chamber's wake
is small, and this exists to make the physics visible — plus the undulator gap wake and
a 100 nm rough surface. The run writes the per-slice energy-loss rate to
`sase_wake.wake.txt` (measured: 1.9 to 121 keV/m across the window, the head slices
losing least — the wake is causal, and the resistive-wall numerical impedance of Bane &
Stupakov sets the shape). Measured effect: the mean energy drops 8.3 m_e c^2 (about
4.2 MeV) over the 57 m line — clearly visible in the energy panel against the 1 m_e c^2
initial spread — while the SASE still reaches 3.9 GW.
