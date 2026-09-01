# Harmonic lasing: the field set

One command, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py harmonics.stats.h5

A strongly seeded steady-state run on a planar segment carrying two radiation
fields (namelist `harmonics = 1, 3`): the fundamental at 1 Angstrom and a third
harmonic at 1/3 Angstrom that starts dark and grows from the bunching the
fundamental drives -- nonlinear harmonic generation, at three times the
fundamental's growth rate (P3 rides |b3|^2, and b3 grows as b1^3).

Measured by this example (200 MW seed, 3.96 m planar, rms aw = 0.84853):

| z (m) | P1 (W) | P3 (W) |
|---|---|---|
| 0.00 | 2.000e8 | 0 |
| 1.98 | 1.975e8 | 7.3e-1 |
| 3.01 | 1.905e8 | 1.7e2 |
| 3.96 (exit) | 1.847e8 | 5.1e3 |

`plot_fel.py`'s harmonic panel shows both curves on one log scale (banked as
`fel-benchmark-plots/harmonic-lasing.png`). The harmonic's full wavefront_params
live under `field/harm3/` in the stats file. Its dump carries `-h3` in the name,
`harmonics-final-h3.wf.h5`, an openPMD EXT_Wavefront file with the photon energy
identifying the harmonic.

A helical undulator would give exactly nothing here: its coupling fc(h) vanishes
for every harmonic (manual: the harmonic-radiation section), which is why this example is planar.
Validated against Genesis 1.3 Version 4 (Genesis4) running the same configuration (fundamental 5.3e-8,
harmonic growth 1.3e-4) and against the exact Bessel deposit sum (3.3e-16) by the
harness's harmonic check section. Runs in ~1 s.
