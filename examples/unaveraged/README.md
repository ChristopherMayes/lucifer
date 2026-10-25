# The FEL with no period averaging

The example runs in two commands, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    ../../build-release/bin/lucifer lucifer_averaged.in
    python ../plot_fel.py unaveraged.stats.h5

This is one seeded steady-state segment tracked with no period averaging and no resonance approximation (the manual's [unaveraged-mode section](../../docs/fel-physics.md)). The particles ride the real helical field, quiver and all, at 20 integration substeps per period with sin^2 entry and exit ramps, and the radiation is a co-evolving kick. Nothing in this path knows the coupling factor `fc`: the energy exchange is what the Lorentz force does.

The mode is `seg1.bmad`'s own `fel_method`, `unaveraged`. The substep count and the ramp length are global switches, `run.unaveraged_steps_per_period` and `run.unaveraged_ramp_periods`, at their defaults here. `lucifer_averaged.in` runs the identical configuration through the averaged default, with a lattice that calls the same file and overrides that one attribute.

On this input the two gain curves agree in ln P at the segment exit to under one part in a thousand, and the table below carries both exit powers. Those are two independent formulations of the same physics, sharing no approximation, and the agreement is what pins the averaged mode's coupling factor. `lucifer_averaged.in` states `run.source_filter = false`, since the unaveraged mode builds no source to filter and the two must compare like for like. The benchmark measures the factor itself to about 6e-4 with dedicated probes ([validation](../../docs/validation.md)).

The run also writes `unaveraged.ledger.txt`: beam energy relative to the reference and window field energy per record, whose sum is conserved. The benchmark checks the closure at 1e-4 of the turnover under a strong-exchange probe.

One segment sits in the lethargy regime, so the plot shows the two curves overlaid rather than a gain curve. The seed diffracts and exponential growth is only beginning at the exit. For the same comparison over a full line, with one segment unaveraged among eleven averaged ones, see `../mixed_line`.

The run takes about five seconds, against under a second for the averaged twin.
