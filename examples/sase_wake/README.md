# Chamber wakes on the SASE run

The example runs in one command, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    python ../plot_fel.py sase_wake.stats.h5

This is the `../sase` run through a deliberately narrow copper chamber of 0.5 mm radius, with the undulator gap wake and a 100 nm rough surface added. All three wake families are transcribed from Genesis 1.3 Version 4 (Genesis4) and act as a per-slice energy-loss rate. The resistive wall acts through the numerical impedance of Bane and Stupakov, the gap wake is convolved with dI/ds, and roughness acts through the complex-q contour (the manual's [wakes section](../../docs/fel-physics.md)).

The chamber is a demonstration case. At 5.8 GeV and 1 Angstrom a normal chamber's wake is small, and this radius makes the physics visible in one run.

On this input the run writes the per-slice energy-loss rate to `sase_wake.wake.txt` at 1.94 to 121.5 keV/m across the window. The head slices lose least, because the wake is causal and the head has little charge ahead of it. The mean energy drops over the 57 m line by the energy change in the table, and the drop varies across the window by the span the table gives at the exit. The drop is visible in the energy panel against the initial spread. The SASE still lases, to an exit power below the `../sase` run's. The measured quantities are in the table below on the example's page. The power of an unseeded run depends on the grid and the particle count of its input ([SASE convergence](../../docs/startup-noise.md)).

The same resistive-wall kernel applied through Bmad's own wake machinery instead of the transcribed model is `../bmad_wake`, which carries both runs and their agreement. What a loss this size does to the ponderomotive phases is `../migration`.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is about a percent of the power inside it, so the exit power is the mode. With the filter off the same input file reports far more power, nearly all of it that emission. The wake's effect is measured on the beam and on the energy panel, which the split does not touch.

The run takes about half a minute.
