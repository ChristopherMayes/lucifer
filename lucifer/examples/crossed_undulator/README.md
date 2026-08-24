# Crossed undulator: the two-polarization afterburner

An x-planar undulator (0.6 m) seeds and microbunches the beam. The same undulator
rotated a quarter turn (`UNDY: UNDX, tilt = pi/2` -- tilt is the polarization spec,
standard Bmad) follows after a short gap. Microbunching is longitudinal and
polarization-blind, so the y set radiates orthogonally polarized light seeded by the
bunching the x set built, while the x field passes through gaining nothing (the
harness holds that isolation at 1.6e-14, manual sec:field vector convention).

    ../../../../production/bin/lucifer run.nml
    python ../plot_fel.py crossed.stats.h5

Measured by this example (1 MW x seed, 0.6 m x-set, 0.2 m gap, 0.6 m y-set):

| z (m) | P_x (W) | P_y (W) |
|---|---|---|
| 0.00 (entry) | 1.0000e6 | 0 |
| 0.60 (x-set exit) | 9.9982e5 | 0 |
| 1.00 (inside the y-set) | 9.9982e5 | 1.63e1 |
| 1.26 (exit) | 9.9982e5 | 1.35e2 |

P_x is flat through the y-set to five digits, while P_y climbs from nothing on the
microbunching the x-set left behind. The tilted undulator cannot amplify the
orthogonal component. It only diffracts it. `plot_fel.py` shows exactly this in its
polarization panel (the last one, present only when a run carries two components).
The figure is banked as `fel-benchmark-plots/crossed-afterburner.png`.

Outputs: `crossed-final-x.fld.h5` / `crossed-final-y.fld.h5` (one polarization per
file, Genesis's format), and `crossed.stats.h5` with the field group carrying totals
plus the x component and a `field/y/` group. Runs in ~30 s.
