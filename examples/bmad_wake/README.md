# The same wake through Bmad's own machinery

The example runs in two commands, with Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_transcribed.in
    python ../plot_fel.py bmad_wake.stats.h5

This is the resistive-wall wake of `../sase_wake`, delivered through Bmad's own wake machinery instead of the transcribed Genesis 1.3 Version 4 (Genesis4) model. `ztable.wake` is the Bane and Stupakov resistive-wall kernel for a 0.5 mm copper chamber, attached as an `sr_wake` `z_long` table to every element of `wake_lattice.bmad`. The undulators apply it once per element at mid-element across the whole window. The quadrupoles and pipes apply it through Bmad's own `track1_bunch`. The drift slots are pipes here, because a Bmad drift cannot carry a wake.

The table was exported by `chamber_wake%write_kernels` and then converted to Bmad's conventions. The sign is flipped to Bmad's positive-decelerating sense, the self-slice term is unhalved, the causal side is at z < 0, and the table is padded past the window. The conversion is described in the header of `wake_lattice.bmad`, and any `chamber_wake%on` run regenerates the exported kernel (the manual's [Bmad element wakes section](../../docs/fel-physics.md)).

`lucifer_transcribed.in` is the comparison partner: the same beam and the same window through the transcribed model with the same chamber, so the two runs differ in the implementation and in nothing else. `../sase_wake` is not the partner, because it adds the gap and roughness terms that this table has no counterpart for.

The transcribed `chamber_wake` applies the wake per integration step, and the Bmad `sr_wake` `z_long` table applies it once per element. The exit mean energy change of each deck is in the table below on the example's page, and the two differ by well under a percent. Per slice across the window interior the agreement is 0.80% on average and 1.73% at worst. The two implementations apply one physical wake at two granularities, and the granularity is what separates them.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is a few percent of the power inside it in both runs, so the exit powers are the mode. With the filter off both decks reported far more power, nearly all of it that emission. The comparison between the two wake implementations is made on the beam's energy loss, which the split does not touch.

Each run takes about half a minute.
