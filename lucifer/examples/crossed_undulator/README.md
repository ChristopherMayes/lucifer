# Crossed undulator: the two-polarization afterburner

An x-planar undulator (0.6 m) seeds and microbunches the beam. The same undulator rotated a quarter turn, `UNDY: UNDX, tilt = pi/2`, follows after a short gap. The tilt is standard Bmad, and it is what states the polarization. Microbunching is longitudinal and blind to polarization, so the y set radiates orthogonally polarized light seeded by the bunching the x set built. The x field passes through the y set gaining nothing, and the harness holds that isolation at 1.6e-14 (the manual's field-solver vector convention).

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py crossed.stats.h5

Measured by this example (1 MW x seed, 0.6 m x-set, 0.2 m gap, 0.6 m y-set):

| z (m) | P_x (W) | P_y (W) |
|---|---|---|
| 0.00 (entry) | 1.0000e6 | 0 |
| 0.60 (x-set exit) | 9.9982e5 | 0 |
| 1.00 (inside the y-set) | 9.9982e5 | 1.63e1 |
| 1.26 (exit) | 9.9982e5 | 1.35e2 |

P_x is flat through the y-set to five digits, while P_y climbs from nothing on the microbunching the x-set left behind. The tilted undulator cannot amplify the orthogonal component. It only diffracts it. `plot_fel.py` shows this in its polarization panel, which is present only when a run carries two components, and that figure is this example's page in the documentation.

Outputs: `crossed-final.wf.h5`, an openPMD EXT_Wavefront file carrying both polarizations as the x and y components of one mesh record, and `crossed.stats.h5`, whose field group carries the total, the x component and a `field/y/` group. Runs in ~1 s.
