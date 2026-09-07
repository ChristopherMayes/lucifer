# Loading a bunch: the keep mode and the sample mode

Two commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_from_file.in
    python ../plot_fel.py import.stats.h5

A bunch described by Bmad's `beam_init_struct`, which is the native equivalent of what
Genesis 1.3 Version 4 (Genesis4)'s `&beam` describes, loaded into FEL slices by each of the
two `load_mode` settings and tracked dark through the full line: SASE from a real bunch
(the manual's [loading section](../../doc/fel-physics.md)).

The lattice is the whole optics specification. The reference energy comes from its
`e_tot`, and the bunch is generated matched to the Twiss in its `beginning`
statement, which is the FODO line's periodic solution. There is no `gamma0` knob and
no `match` transform, because both would be a second way of saying what the lattice
already says.

The bunch is a Gaussian test bunch sized to the FEL window's economics, and it is
labeled a test bunch for a reason: `sig_z = 1.2 nm` with 30 fC gives the benchmark
3 kA peak in a window a few hundred wavelengths long. A physical accelerator bunch
is micrometers and hundreds of pC, needs thousands of slices, and belongs to a
production run rather than to an example.

`lucifer.in` generates 50000 particles and keeps every one: `load_mode = "keep"` makes
each particle a beamlet of eight phase copies at an eighth of its charge, so the slices
hold as many macroparticles as the bunch put there, 16 in the tails and 40320 at the
peak, 400000 in all over 34 slices, and every charge-weighted moment of the bunch is the
load's. The window is the bunch's own extent. `quiet_start = T` spreads the copies over
2 pi about each particle's phase, and `shot_noise = T` imposes the physical noise on the
beamlets, after measuring that the load is quiet, 2.6e-29 in |b|^2 N_lambda before the
noise goes on. The deck also writes the bunch with `write_openpmd_file` as
`import_bunch.h5`, 2.8 MB of openPMD.

`lucifer_from_file.in` names that file in `beam_init%position_file` and loads it in the
sample mode: `beam_init%n_particle = 2048` macroparticles in every slice, resampled from
the particles that fell into each slice's sampling window by Genesis4's
`importdistribution` method, with the current taken from the same window. The same
34 slices come out, since both modes bin the same bunch on the same spacing. `write_genesis_dist`
is the third direction, handing the identical bunch to Genesis4's `&importdistribution`.

Measured on these inputs: the keep mode exits at 42.2 kW after 57 m in 45 s, and the
sample mode at 84.3 kW in 12 s. Both are a startup from a real bunch rather than a
saturated FEL, and the two are different SASE realizations of the same physics: the keep
mode draws no random number before the noise and the sample mode draws for every
resampled particle, so the noise falls on different phasors. `../sase` describes the same
sensitivity: in a dark start, a small change re-rolls the realization. The keep mode
costs what its 400000 macroparticles cost. A production run with a bunch from an
accelerator simulation would keep fewer copies or sample.

The run reports the split between the mode and the wide-angle emission of the point beamlets, which is the emission of macroparticles that occupy one grid point each and is no part of what a real beam radiates ([SASE convergence](../../doc/startup-noise.md)). At the exit the power outside the split angle is 3 percent of the power inside it in both modes. With the filter off the same decks reported 3.5 and 11.2 MW, of which 99 percent was that emission, which is why the exit powers here are eighty times lower and are the mode. Both modes load the same bunch on the same grid, so the comparison between them is unaffected.
