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
sizes start on the matched values. To plot:

```python
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("steady_state.diag.txt")
plt.semilogy(d[:, 0], d[:, 2]); plt.xlabel("z [m]"); plt.ylabel("power [W]"); plt.show()
```

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
