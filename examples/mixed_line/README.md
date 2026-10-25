# Two FEL methods in one line

The example runs in two commands, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    ../../build-release/bin/lucifer lucifer_averaged.in
    python ../plot_fel.py mixed_line.stats.h5

The FEL method is an element attribute, so a line whose segments are tracked differently is a line of elements that differ. `mixed.bmad` calls `../aramis.bmad` for its parameters, defines `UNDU` as the same undulator carrying `fel_method = unaveraged`, and puts one `UNDU` in the third FODO cell. Eleven segments run the averaged default and the fifth runs unaveraged, integrating the real helical field with no period averaging (the manual's [unaveraged-mode section](../../docs/fel-physics.md)).

`lucifer_averaged.in` is the all-averaged twin, which overrides that one element back to the default. The two runs differ in one segment's method and in nothing else.

On this input, a seeded steady state, the two runs are bit-identical up to the first record inside `UNDU`, at z = 19.045 m, and differ from there. The measured quantities are in the table below, the exit power of each run and the step count that the one unaveraged segment adds among them.

In the plot the log-power panel runs straight through the boundary at 19 to 23 m, with no step and no kink where the method changes, and saturation arrives at the same place on both lines. A difference of a few percent is not visible on a six-decade axis, which is why the table is the check and the plot is the sanity check.

The two curves are bit-identical up to the method boundary, which is the check that the mixing itself is what differs. Through the unaveraged segment they separate by a few percent in ln P, and they stay within that at the exit. Two formulations of the same physics, sharing no approximation, hand the same beam and field back and forth across an element boundary.

The cost is where the difference is spent. The one unaveraged segment takes most of the run's wall time at 20 substeps per period, and the eleven averaged segments share the rest. Per-element selection puts the expensive method in the segment where it is wanted and the cheap one everywhere else.

Collective effects are refused in a line with any unaveraged segment, so this example carries none, and the benchmark checks the refusal.

The benchmark runs the same configuration as its sandwich check. There the mixed line must agree with the all-averaged line at the exit, the energy ledger rows must appear only inside the unaveraged segment, and a wake on that segment must be refused ([validation](../../docs/validation.md)). This directory is the runnable version of that check.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is several times the power inside it in both runs, so most of each exit power is that emission. The table carries the exit power inside and outside the split for both runs, so the difference between them can be read in the mode as well as in the total. Both input files state `run.source_filter = false`, so the two FEL methods compare like for like. The filter reaches the averaged elements only and would filter eleven segments of one line and twelve of the other.

The run takes a few seconds, and the averaged twin about a second.
