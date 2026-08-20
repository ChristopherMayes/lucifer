# Crossed undulator: the two-polarization afterburner

An x-planar undulator (0.6 m) seeds and microbunches the beam; the SAME undulator
rotated a quarter turn (`UNDY: UNDX, tilt = pi/2` -- tilt is the polarization spec,
standard Bmad) follows after a short gap. Microbunching is longitudinal and
polarization-blind, so the y set radiates ORTHOGONALLY polarized light seeded by the
bunching the x set built, while the x field passes through gaining nothing (the
harness holds that isolation at 1.6e-14; manual sec:field vector convention).

    ../../../../production/bin/fel_track_test run.nml
    python ../plot_fel.py crossed.stats.h5

Outputs: `crossed-final-x.fld.h5` / `crossed-final-y.fld.h5` (one polarization per
file, Genesis's format), and `crossed.stats.h5` with the field group carrying totals
plus the x component and a `field/y/` group. Runs in ~30 s.
