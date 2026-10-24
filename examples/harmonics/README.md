# Harmonic lasing: the field set

The example runs in one command, with Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py harmonics.stats.h5

This is a strongly seeded steady-state run on a planar segment carrying two radiation fields (namelist `harmonics = 1, 3`): the fundamental at 1 Angstrom and a third harmonic at 1/3 Angstrom. The harmonic starts dark and grows from the bunching the fundamental drives. That is nonlinear harmonic generation, at three times the fundamental's growth rate, since P3 rides |b3|^2 and b3 grows as b1^3.

The seed is 200 MW on 3.96 m of planar undulator at rms aw = 0.84853. Over the segment the fundamental gives up a small fraction of its seed power to the beam, and the third harmonic climbs from nothing through several decades. The measured quantities are in the table below. The exit power there is the fundamental's, and the harmonic's power along z is in the plot.

`plot_fel.py`'s harmonic panel shows both curves on one log scale, and that figure is this example's page in the documentation. The harmonic's full wavefront_params live under `field/harm3/` in the stats file. Its dump carries `-h3` in the name, `harmonics-final-h3.wf.h5`, an openPMD EXT_Wavefront file with the photon energy identifying the harmonic.

A helical undulator would give nothing here, since its coupling fc(h) vanishes for every harmonic (the manual's harmonic-radiation section). That is why this example is planar. The deck turns the source filter off, as the comparison against Genesis 1.3 Version 4 (Genesis4) requires. The harness's harmonics section runs the same configuration against Genesis4, on the fundamental's power and on the harmonic's growth, and against the exact Bessel deposit sum, at the levels in [](../../docs/validation.md#val-levels-harmonics).

The run takes about a second.
