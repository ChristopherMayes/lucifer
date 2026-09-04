# Two tracking methods in one line

Two commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_averaged.in
    python ../plot_fel.py mixed_line.stats.h5

The tracking method is an element attribute, so a line whose segments are tracked
differently is a line of elements that differ. `mixed.bmad` calls `../aramis.bmad`
for its parameters, defines `UNDU` as the same undulator carrying
`tracking_method = fel_unaveraged`, and puts one `UNDU` in the third FODO cell. Eleven
segments run the averaged default and the fifth runs unaveraged, integrating the real
helical field with no period averaging (the manual's [unaveraged-mode section](../../doc/fel-physics.md)).

`lucifer_averaged.in` is the all-averaged twin, which overrides that one element back
to the default. The two runs differ in one segment's method and in nothing else.

Measured on this input (seeded steady state, 129-point grid):

| | value |
|---|---|
| first difference between the runs | z = 19.045 m, the first record inside `UNDU` |
| largest difference along the line | 2.68e-2 in ln P, at z = 24.7 m |
| difference at the exit | -7.39e-3 in ln P (1.757 GW against 1.770 GW) |

In the plot the log-power panel runs straight through the boundary at 19 to 23 m, with
no step and no kink where the method changes, and saturation arrives at 37 m as it does
on the all-averaged line. A 2.7% difference is not visible on a six-decade axis, which
is why the table above is the actual check and the plot is the sanity check.

The two curves are bit-identical up to the method boundary, which is the check that
the mixing itself is what differs. Through the unaveraged segment they separate by
2.7%, and by the exit they have converged back to 0.74%: two formulations of the same
physics, sharing no approximation, handing the same beam and field back and forth
across an element boundary.

The cost is where the difference is spent. The one unaveraged segment takes about 5 s
of the 6 s run at 20 substeps per period, and the eleven averaged segments share the
rest. That is what per-element selection buys: the expensive method in the segment
where it is wanted, the cheap one everywhere else.

Collective effects are refused in a line with any unaveraged segment, so this example
carries none, and the refusal itself is checked by the harness.

The harness runs the same configuration as its sandwich check, where the mixed line must
agree with the all-averaged line at the exit, the energy ledger rows must appear only
inside the unaveraged segment, and a wake on that segment must be refused
([validation](../../doc/validation.md)). This directory is the runnable version of that
check.

Runs in ~7 s, and the averaged twin in ~1 s.
