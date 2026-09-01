# The FEL with no period averaging

Two commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_averaged.in
    python ../plot_fel.py unaveraged.stats.h5

One seeded steady-state segment tracked with no period averaging and no resonance
approximation (the manual's [unaveraged-mode section](../../doc/fel-physics.md)). The particles ride the
real helical field, quiver and all, at 20 integration substeps per period with sin^2
entry and exit ramps, and the radiation is a co-evolving kick. Nothing in this path
knows the coupling factor `fc`: the energy exchange is what the Lorentz force does.

The mode is `seg1.bmad`'s own tracking method, `fel_unaveraged`, and the substep count and
ramp length are element attributes too, `fel_steps_per_period` and `fel_ramp_periods`.
`lucifer_averaged.in` runs the identical configuration through the averaged default by
calling the same lattice and overriding that one attribute.

Measured on this input: the two gain curves agree to 7.6e-4 in ln P at the segment
exit, 4.421 kW against 4.418 kW. That is two independent formulations of the same
physics, sharing no approximation, and the agreement is what pins the averaged mode's
coupling factor. The harness measures the factor itself to about 6e-4 with dedicated
probes ([validation](../../doc/validation.md)).

The run also writes `unaveraged.ledger.txt`: beam energy relative to the reference and
window field energy per record, whose sum is conserved. The harness checks the closure
at 1e-4 of the turnover under a strong-exchange probe.

One segment sits in the lethargy regime, so the plot shows the two curves overlaid
rather than a gain curve. The seed diffracts and exponential growth is only beginning
at the exit. For the same comparison over a full line, with one segment unaveraged
among eleven averaged ones, see `../mixed_line`.

Runs in ~5 s, against under a second for the averaged twin.
