# Crossed undulator: the two-polarization afterburner

An x-planar undulator (0.6 m) seeds and microbunches the beam. The same undulator rotated a quarter turn, `UNDY: UNDX, tilt = pi/2`, follows after a short gap. The tilt is standard Bmad, and it is what states the polarization. Microbunching is longitudinal and blind to polarization, so the y set radiates orthogonally polarized light seeded by the bunching the x set built. The x field passes through the y set gaining nothing, and the validation suite holds that isolation at the level in [](../../docs/validation.md#val-levels-two-polarization) (the manual's field-solver vector convention).

    ../../build-release/bin/lucifer lucifer.toml
    python -m lucifer.plot crossed.stats.h5

With a 1 MW x seed, a 0.6 m x-set, a 0.2 m gap and a 0.6 m y-set, P_x is flat through the y-set to five digits, while P_y climbs from nothing on the microbunching the x-set left behind. The measured quantities are in the table below, the total exit power and the power inside the split angle among them, and the plot carries P_y along z. The tilted undulator cannot amplify the orthogonal component. It only diffracts it. The summary figure shows this in its polarization panel, which is present only when a run carries two components, and that figure is this example's page in the documentation.

The outputs are `crossed-final.wf.h5`, an openPMD EXT_Wavefront file carrying both polarizations as the x and y components of one mesh record, and `crossed.stats.h5`, whose field group carries the total, the x component and a `field/y/` group. The run takes about a second.
