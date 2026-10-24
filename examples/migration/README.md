# Slice migration: particles that outgrow their slice

The example runs in two commands, with Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_no_migrate.in
    python ../plot_fel.py migration.stats.h5

A particle's ponderomotive phase drifts, and once it leaves the slice window the particle belongs to a neighboring slice. `global%migrate = T` moves it there. This run is the `../sase` configuration at 1024 particles per slice with migration on, and `lucifer_no_migrate.in` is the same run with it off.

No artificial driver is needed to make this happen. The phase advance a particle accumulates depends on its energy and on its betatron angles, both of which vary across the beam, so a real beam in a real line migrates on its own. Adding `../sase_wake`'s chamber block to this deck raises both the moves and the dropped charge. A wake that costs the tail slices 121 keV/m and the head slices 1.9 keV/m (`../sase_wake`) is a systematic phase walk on top of the beam's own spread.

On this input the run moves particles between slices at 48 element ends and drops a small fraction of the beam's 2.88e-13 C off the window ends. `migration.migration.txt` carries the totals. The measured quantities are in the table below, the exit power and the exit bunching with migration on and off among them.

The exit powers are those of this input's grid and particle count. The power of an unseeded run depends on both ([SASE convergence](../../docs/startup-noise.md)). Each run is also one SASE realization on a 96-slice window. Over four seeds the ratio of the two exit powers falls on either side of one, so the exit power of a single pair does not measure what migration is worth. The exit bunching is the higher with migration on in every one of those seeds. The bunching is a mean over 96 slices where the exit power follows the largest spike, which is why one of the two separates the runs and the other is swamped.

`migration.migration.txt` carries one row per event: the s position, the count, the charge dropped, the phasor deviation, and the smallest phase distance any particle at that event had to a slice boundary. The totals are at the end. The rows are what a plot of migration would be drawn from.

The per-slice current shows the difference most directly. With `migrate = F` it cannot change. All 96 slices sit at exactly 3000 A at the exit, the value they were loaded with, because no particle ever leaves the slice it started in. With `migrate = T` the profile evolves across the window. The drain is at the high-index end of the window, where the last two slices are the emptiest. In the plot, the emptiest slice is the one gray line whose bunching climbs far above every other slice's, to the largest slice bunching the table carries. That is a small-population artifact and not physics. A slice with a small fraction of its charge left has a noisy bunching estimate, and `n_eff` in the diag columns is what says so.

Migration changes the physics and not only the bookkeeping. With migration off, a particle whose phase has left its slice keeps radiating into the field of the slice it started in, and the bunching each slice sees is smeared by contributions that belong elsewhere. On this configuration that shows in the exit bunching, which the table carries for both runs. The default is off because Genesis 1.3 Version 4 (Genesis4) does not migrate without one4one, and the comparison tiers would otherwise measure a model difference rather than a transcription ([validation](../../docs/validation.md), the slice-migration section).

`migrate_check = T` is the check and not the feature. It verifies that the whole-beam weighted phasor obeys `S_before = S_after + S_dropped` across every move, and the worst deviation at the end of `migration.migration.txt` is that residual. Migration under weights generalizes a Genesis4 method that requires one4one, and the check is what validates the generalization.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is a few percent of the power inside it in both runs, so both exit powers are the mode. With the filter off the same decks reported far more power, most of it that emission.

Each run takes about a quarter of a minute.
