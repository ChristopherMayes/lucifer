# Pure SASE: the time-dependent run with nothing external

The example runs in one command, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    python -m lucifer.plot sase.stats.h5

This is self-amplified spontaneous emission from nothing but the electrons. The loader generates a time window of slices spaced an integer number of wavelengths apart, imposes physical shot noise by weighted Fawley loading, and starts the field dark (`seed_power = 0`). The FEL then grows through the full line with slippage active. Before imposing the noise the loader prints the per-slice electron count `N_lambda`, the effective count `N_eff`, and the quiet floor it verified. The noise level is therefore a reported number.

The exit power of a run like this one is set by the transverse cell size and the macroparticle count as much as by the physics, since each beamlet radiates into every angle the grid carries. The SASE convergence page measures that and says how to choose the grid, the window and the load. `run.source_filter = true` suppresses the wide-angle part, with its edge derived from the beam and the gain and printed at setup. Every run also reports `power_inside_angle` beside the total power at each element end, so the separation is in the stats file without a dump.

The beam is a flat coasting bunch in Bmad's own vocabulary: `distribution_type(3) = "GRID"` over the window extent. The current is derived as `I = Q*c/extent`, so `bunch_charge` encodes 3 kA exactly.

On this input (`ran_seed` at its 12345 default) the startup power settles after the first segment, grows through the line, and fluctuates from slice to slice at the exit, which is the SASE fluctuation. The slice energy spread grows by a few percent over the line. The measured quantities are in the table below on the example's page: the power per slice after the first element, the exit power, its rms over mean across the slices, and the energy spread at the start and the exit. With the source filter off the same input file reports nearly two orders of magnitude more power at the exit, almost all of it the wide-angle emission of the point beamlets ([the convergence section](../../docs/fel-physics.md#sec-convergence)).

Two features of the plot are physics and not artifacts. The total-power sawtooth is radiation slipping out of the head of a finite time window at each drift while fresh vacuum enters at the tail. Genesis 1.3 Version 4 (Genesis4) shows the same sawtooth, as does any finite window. Deep saturation of every slice needs a window longer than the total slippage. The per-slice lines in the power panel are the slippage cascade itself.

`run.keep_escaped_field = true` banks what slips out and reconstructs the whole pulse at the exit. It is off here for two reasons. On this configuration the bank and the reconstructed pulse come to more than a gigabyte each. And those two files are the one place the tracker still writes Genesis4 field conventions rather than openPMD ([reading an output file](../../docs/reading-output.md)).

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). The filter is on by default and this input file runs with it. At the exit the power outside the split angle is a few percent of the power inside it, so the exit power is the mode. With the filter off the same input file reports far more power, nearly all of it that emission. Four times the macroparticles is not enough on this window without the filter. Leaving the grid out of the input file derives a coarser grid with half the points of this one. Two grids a factor of two apart in cell size then agree to within ten percent once the artifact is gone.

The run takes about half a minute.
