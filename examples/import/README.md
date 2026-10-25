# Loading a bunch: the keep mode and the sample mode

The example runs in two commands, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    ../../build-release/bin/lucifer lucifer_from_file.in
    python ../plot_fel.py import.stats.h5

The bunch is described by Bmad's `beam_init_struct`, the native equivalent of what Genesis 1.3 Version 4 (Genesis4)'s `&beam` describes. It is loaded into FEL slices by each of the two `load_mode` settings and tracked dark through the full line: SASE from a real bunch (the manual's [loading section](../../docs/fel-physics.md)).

The lattice is the whole optics specification. The reference energy comes from its `e_tot`, and the bunch is generated matched to the Twiss in its `beginning` statement, which is the FODO line's periodic solution. There is no `gamma0` knob and no `match` transform, because both would be a second way of saying what the lattice already says.

The bunch is a Gaussian test bunch sized to what the FEL window can afford: `sig_z = 1.2 nm` with 30 fC gives the benchmark 3 kA peak in a window a few hundred wavelengths long. A physical accelerator bunch is micrometers long and carries hundreds of pC. It needs thousands of slices and belongs to a production run rather than to an example.

`lucifer.toml` generates 50000 particles and keeps every one. `load_mode = "keep"` makes each particle a beamlet of eight phase copies at an eighth of its charge. The slices then hold as many macroparticles as the bunch put there, a handful in the tails and tens of thousands at the peak, 400000 in all. Every charge-weighted moment of the bunch is the load's, and the window is the bunch's own extent. `quiet_start = true` spreads the copies over 2 pi about each particle's phase. `shot_noise = true` imposes the physical noise on the beamlets after measuring that the load is quiet, 2.6e-29 in |b|^2 N_lambda before the noise goes on. The input file also writes the bunch with `write_openpmd_file` as `import_bunch.h5`, 2.8 MB of openPMD.

`lucifer_from_file.in` names that file in `beam.position_file` and loads it in the sample mode. `beam.n_particle = 2048` macroparticles go into every slice, resampled by Genesis4's `importdistribution` method from the particles that fell into each slice's sampling window, with the current taken from the same window. The same slices come out, since both modes bin the same bunch on the same spacing. `write_genesis_dist` writes the same bunch as a Genesis4 distribution file for its `&importdistribution`.

On these inputs both modes exit far below saturation, a startup from a real bunch rather than a saturated FEL. The measured quantities are in the table below, the slice count and the peak slice's macroparticle count among them. The two exit powers differ, and the two runs are different SASE realizations of the same physics. The keep mode draws no random number before the noise and the sample mode draws for every resampled particle, so the noise falls on different phasors. In a dark start a small change in the load re-rolls the realization, as the per-slice spread of `../sase` shows. The keep mode costs what its 400000 macroparticles cost. A production run with a bunch from an accelerator simulation would keep fewer copies or sample.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is a few percent of the power inside it in both modes, and the table carries the ratio. With the filter off the same input files reported far more power, nearly all of it that emission, which is why the exit powers here are the mode. Both modes load the same bunch on the same grid, so the comparison between them is unaffected.
