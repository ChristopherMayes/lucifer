# Chamber wakes on the SASE run

One command, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py sase_wake.stats.h5

The `../sase` run through a deliberately narrow copper chamber of 0.5 mm radius,
plus the undulator gap wake and a 100 nm rough surface. All three wake families are
transcribed from Genesis 1.3 Version 4 (Genesis4) and act as a per-slice energy-loss
rate: the resistive wall through the numerical impedance of Bane and Stupakov, the
gap wake convolved with dI/ds, and roughness through the complex-q contour
(the manual's [wakes section](../../doc/fel-physics.md)).

The chamber is a tuned demonstration case, labeled as one. At 5.8 GeV and
1 Angstrom a normal chamber's wake is small, and this radius exists to make the
physics visible in one run.

Measured on this input: the run writes the per-slice energy-loss rate to
`sase_wake.wake.txt` at 1.94 to 121.5 keV/m across the window, the head slices losing
least because the wake is causal and the head has little charge ahead of it. The mean
energy drops 8.29 m_e c^2, about 4.24 MeV, over the 57 m line, and the per-slice drop
runs from 0.25 to 13.56 m_e c^2 across the window. That is clearly visible in the
energy panel against the 1 m_e c^2 initial spread. The SASE still reaches 2.94 GW,
against 3.02 GW for the same run with no wake. Both powers are at this deck's grid and
particle count, and a dark start's power depends on both
([startup noise](../../doc/startup-noise.md)).

The same resistive-wall kernel applied through Bmad's own wake machinery instead of
the transcribed model is `../bmad_wake`, which carries both runs and their agreement.
What a loss this size does to the ponderomotive phases is `../migration`.

Runs in ~25 s.
