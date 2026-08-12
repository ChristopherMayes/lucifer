# FEL examples

Self-contained runs of the FEL tracker: one command, no Genesis, no dump files. For the
validation benchmark against Genesis — which is where the physics is proven — see
`bsim/fel/tests/`.

## Steady state

```
cd bsim/fel/examples
../../../debug/bin/fel_track_test steady_state.nml     # or production/bin
```

This simulates a seeded, steady-state (single-slice) FEL: the Benchmark1-SASE
configuration — a 6-cell FODO line of 3.99 m helical undulator segments (aw = 0.85 rms,
15 mm period) at 5.8 GeV, resonant at 1 Angstrom, 3 kA, with a 5 kW Gaussian seed. With
no dump files named in the namelist, the program generates its own starting state: a
quiet-start beam (beamlet loading, so the initial bunching is zero to roundoff and the
FEL grows from the seed, not from sampling noise) and the seed field at its waist. The
run takes about a minute and prints three output files.

`steady_state.diag.txt` is the gain curve: one row per integration step with columns
`z, slice, power, on_axis_intensity, bunching, bunching_phase, mean_gamma, sigma_gamma,
sigma_x, sigma_y` (slice is always 1 here). Measured on this input: the seed enters at
exactly 5 kW with bunching 4e-17 (the quiet start is exact), power grows exponentially
with a 7.5 m power gain length and saturates at 1.6 GW around z = 38 m, and the beam
sizes start on the matched values. To plot (needs matplotlib; the bmad-fel-validate
environment has it):

```
python plot_fel.py steady_state.diag.txt        # writes steady_state.diag.png
```

Four panels against z: radiation power (log), bunching, energy change and rms spread
(both in units of m_e c^2 on one axis), and the transverse rms sizes showing the FODO
betatron oscillation. The same script reads time-dependent diag files -- give it one
from the benchmark's td tiers and each thin gray line is a slice, with the bold line the
total power or the slice average; the slippage echelon is directly visible in the
per-slice power.

`steady_state-final.fld.h5` and `steady_state-final.par.h5` are the end state in Genesis
dump format (readable by `openPMD-beamphysics`'s `Wavefront.from_genesis4`, h5py, or
Genesis itself). Add `write_initial = T` to the namelist to also dump the generated
starting state — useful for handing the identical initial condition to Genesis
(`&importbeam` / `&importfield`). That is how this loader was validated: Genesis imported
the generated dumps and tracked the full line, agreeing with this tracker (transcribed
interlude model) at 1.5e-5 — the constants floor of every Genesis comparison here.

Everything is set in `steady_state.nml`; the header of `fel_track_test.f90` documents
every parameter. The undulator segments are the elements named `UND*` in `aramis.bmad`,
with the FEL parameters (`und_aw`, `und_lambdau`, ...) coming from the namelist — a real
FEL element type carrying them on the lattice is a later deliverable.

## SASE, time dependent

```
../../../debug/bin/fel_track_test sase.nml     # ~90 s (OpenMP over slices)
python plot_fel.py sase.diag.txt
```

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

## SASE with wakefields

```
../../../debug/bin/fel_track_test sase_wake.nml     # ~100 s
python plot_fel.py sase_wake.diag.txt
```

The same SASE run through a deliberately NARROW copper chamber — 0.5 mm radius, a tuned
demonstration case, honestly labeled: at 5.8 GeV and 1 Angstrom a normal chamber's wake
is small, and this exists to make the physics visible — plus the undulator gap wake and
a 100 nm rough surface. The loader writes the per-slice energy-loss rate to
`sase_wake.wake.txt` (measured: 1.9 to 121 keV/m across the window, the head slices
losing least — the wake is causal, and the resistive-wall numerical impedance of Bane &
Stupakov sets the shape). Measured effect: the mean energy drops 8.3 m_e c^2 (about
4.2 MeV) over the 57 m line — clearly visible in the energy panel against the 1 m_e c^2
initial spread — while the SASE still reaches 3.9 GW.
