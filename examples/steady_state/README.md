# A seeded steady-state FEL: the first example to run

The example runs in one command, with Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py steady_state.stats.h5

This is a seeded, single-slice FEL on the benchmark line (`../aramis.bmad`): a 3 kA slice and a 5 kW Gaussian seed at its waist are tracked through 6 FODO cells of helical undulator. Nothing external is read. The program generates the quiet-start beam and the seed field from the namelist, and the lattice supplies everything else.

`sig_z = 0` makes the run steady state. The whole charge sits in the one slice window, so the current is `I = Q*c/slice_spacing`, and `bunch_charge` in the deck encodes 3 kA exactly. The current is always derived this way and never entered ([the input reference](../../doc/input-reference.md)).

On this input the quiet start is exact to the initial bunching in the table, so the FEL grows from the seed and not from sampling noise. The power grows exponentially from the seed, saturates past the middle of the line, and falls back at the exit as particles rotate in the bucket. The gain length, the peak power with its position, and the exit power are in the table below on the example's page. The beam sizes start on the lattice's matched values, since `init_beam_distribution` generates the bunch from the Twiss in the lattice's `beginning` statement.

The two `dump_*_at` lines write the beam and the field at the third undulator, `steady_state-at9-UND.beam.h5` and `-at9-UND.wf.h5`, both openPMD. `at9` is the element index the locator resolved. Mid-run dumps are how a run is handed to another program or restarted from the middle, and a single slice keeps them under 1.5 MB here. The same two lines on the `../sase` window write about a hundred megabytes.

Read [reading an output file](../../doc/reading-output.md) for the stats file and the ten plot panels. The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../doc/startup-noise.md)). The table gives the power outside the split angle over the power inside it, at the mode peak and at the exit. Both ratios are well under a percent, so the exit power is the mode. With the filter off the same deck reports more power at the exit and at the peak, and nearly half of that exit power is this emission.

The run takes about five seconds.
