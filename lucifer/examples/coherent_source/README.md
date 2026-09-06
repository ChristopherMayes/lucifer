# The coherent source: converged results at a fraction of the particles

Three commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_low_m.in
    ../../../production/bin/lucifer lucifer_reference.in
    python ../plot_fel.py coherent.stats.h5

The use case (manual: the coherent-source section): runs whose beams are transversely
Gaussian (idealized machines, parameter scans, seeded amplifier studies) can
trade the per-particle source deposit for Tanaka's coherent retrieval
(`global%source_model = "coherent"`, PRAB 27, 030703 (2024)) and converge at one
to two orders of magnitude fewer macroparticles. The reason: with few particles
the per-cell source is spiky, and that sampling noise always overestimates the
gain. The coherent source carries the slice's exact bunching phasor on a fitted
Gaussian instead, so the artifact never enters.

Three runs of the same seeded steady-state case (the full Aramis line, 5.2
decades of gain, ../steady_state's configuration):

| run | source, particles/slice | exit power | vs reference | wall (12 threads) |
|---|---|---|---|---|
| `lucifer_reference.in` | deposit, 8192 | 7.62e+08 W | -- | 4.3 s |
| `lucifer_low_m.in` | deposit, 512 | 5.82e+09 W | 7.6x high (ln +2.03) | 2.0 s |
| `lucifer.in` | coherent, 512 | 7.46e+08 W | ln 0.021 | 2.2 s |

The middle row is the trap this feature exists to remove: cutting particles
without the coherent source multiplies the predicted power by 7.6 on
this case. Nothing about the power curve says so, and the run's own convergence report
is what does: the power outside the split angle is 8.5 times the power inside it at the
exit on that row, against 0.83 on the reference row and 0.050 with the coherent source
([SASE convergence](../../doc/startup-noise.md)). The power the middle row adds is the
wide-angle emission of 64 beamlets, which the coherent source does not have because it
never puts a macroparticle on a grid point. The last row is the same
particle count giving the converged answer. (Per-slice cost dominates real
time-dependent runs, where the same particle reduction pays proportionally.)

The guardrails are part of the feature (all refusals, the manual's coherent-source section). A per-slice Gaussianity test is sized against its own
sampling significance. A genuinely structured profile refuses. An offset,
mismatched or tilted Gaussian beam passes -- the source centers and tilts with
the beam's phasor-weighted moments. Harmonics, two polarizations and the
unaveraged mode are out of scope in v1. Dark starts are refused, at a measured
~175x startup deficit: SASE grows from spontaneous, spatially-incoherent
emission, which is exactly what the coherent model drops. Seed the field, as
here, or use the default deposit.

Plot any run with ../plot_fel.py <root>.stats.h5.
