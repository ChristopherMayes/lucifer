# A post-saturation taper: heterogeneity as different elements

The example runs in one command, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    python ../plot_fel.py taper.stats.h5

This is the `steady_state` run over `taper.bmad`. That lattice is `../aramis.bmad` for the first four FODO cells, and the last two cells use a second undulator definition, `UND2`, whose `b_max`, and so `aw`, is 0.4% lower. A step-down taper re-matches the resonance to the beam the FEL has already decelerated, so the power keeps climbing past the point where the untapered line saturates and turns over.

The FEL parameters are lattice attributes, so a heterogeneous line is a line of elements that differ. Each element carries its own `b_max`, `l_period` and `field_calc`, and the tracker derives each segment's `aw` from the element it is standing in. There is no per-segment namelist input and nothing to keep in step with the lattice.

Against `../steady_state` with the same seed and the same starting state, the untapered `../aramis.bmad` saturates and then falls back, and `taper.bmad` is still climbing at the exit. The exit power and the peak power of each line are in the table below on its example's page. On the tapered line the peak sits at the exit.

The step is at z = 38.0 m, where the fourth FODO cell ends. The two gain curves are bit-identical until the first record inside the first `UND2`, at z = 38.045 m. That is the check that the two runs differ in the taper and in nothing else. The exit power is an order of magnitude above the untapered exit and several times the untapered saturation peak.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). The table gives the power outside the split angle over the power inside it at the exit. On the tapered line it is well under a percent, so the exit power is the mode and the taper's gain is the mode's. With the filter off the same input file reports more power, the excess being that emission.

The run takes about four seconds.
