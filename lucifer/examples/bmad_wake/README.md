# The same wake through Bmad's own machinery

Two commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_transcribed.in
    python ../plot_fel.py bmad_wake.stats.h5

The chamber-wake physics of `../sase_wake`, delivered through Bmad's own wake
machinery instead of the transcribed Genesis 1.3 Version 4 (Genesis4) model.
`ztable.wake` is the Bane and Stupakov resistive-wall kernel for a 0.5 mm copper
chamber, attached as an `sr_wake` `z_long` table to every element of
`wake_lattice.bmad`. The undulators apply it once per element at mid-element across
the whole 96-slice window, and the quadrupoles and pipes apply it through Bmad's own
`track1_bunch`. The drift slots are pipes here, because a Bmad drift cannot carry a
wake.

The table was exported by `chamber_wake%write_kernels` and then converted to Bmad's
conventions: sign flipped to Bmad's positive-decelerating sense, the self-slice term
unhalved, the causal side at z < 0, and padding past the window. The recipe is in
`wake_lattice.bmad`'s header, and any `chamber_wake%on` configuration will
regenerate the kernel (the manual's [Bmad element wakes section](../../doc/fel-physics.md)).

`lucifer_transcribed.in` is the comparison partner: the same beam and the same window
through the transcribed model with the same chamber, so the two runs differ in the
implementation and in nothing else. `../sase_wake` is not the partner, because it adds
the gap and roughness terms that this table has no counterpart for.

| implementation | exit mean energy change | granularity |
|---|---|---|
| transcribed `chamber_wake` | -2.319 m_e c^2 | per integration step |
| Bmad `sr_wake` `z_long` table | -2.335 m_e c^2 | once per element |

The two means differ by 0.69%. Per slice across the window interior the agreement is
0.80% on average and 1.73% at worst. One physical wake, two independent
implementations, two application granularities, and what separates them is the
granularity rather than the physics.

Runs in ~25 s.
