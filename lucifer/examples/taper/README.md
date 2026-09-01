# A post-saturation taper: heterogeneity as different elements

One command, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py taper.stats.h5

The `steady_state` run over `taper.bmad`, which is `../aramis.bmad` for the first
four FODO cells and then switches to a second undulator definition, `UND2`, whose
`b_max` and therefore `aw` is 0.4% lower. A step-down taper re-matches the
resonance to the beam the FEL has already decelerated, so the power keeps climbing
past the point where the untapered line saturates and turns over.

This is what driving the FEL from lattice attributes buys. A heterogeneous line is
a line of elements that differ, each carrying its own `b_max`, `l_period` and
`field_calc`, and the tracker derives each segment's `aw` from the element it is
standing in. There is no per-segment namelist input and nothing to keep in step
with the lattice.

Measured against `../steady_state` (same seed, same starting state):

| line | P at 57 m | behavior |
|---|---|---|
| `../aramis.bmad` (untapered) | 761.5 MW | saturates at 1.62 GW at z = 37.24 m, then falls back |
| `taper.bmad` (0.4% step at z = 38.0 m) | 9.64 GW | still climbing at the exit |

The step is at z = 38.0 m, where the fourth FODO cell ends. The two gain curves are
bit-identical until the first record inside the first `UND2`, at z = 38.045 m, which
is the check that the two runs differ in the taper and in nothing else. The exit
power is 12.66x the untapered exit and 5.95x the untapered saturation peak.

Runs in ~4 s.
