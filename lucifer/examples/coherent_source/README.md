# The coherent source: converged results at a fraction of the particles

The use case (manual sec:coherent-source): runs whose beams are transversely
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
| `run_reference.nml` | deposit, 8192 | 7.62e+08 W | -- | 4.3 s |
| `run_low_m.nml` | deposit, 512 | 5.82e+09 W | 7.6x high (ln +2.03) | 2.0 s |
| `run.nml` | coherent, 512 | 7.46e+08 W | ln 0.021 | 2.2 s |

The middle row is the trap this feature exists to remove: cutting particles
without the coherent source silently multiplies the predicted power by 7.6 on
this case, with nothing anywhere reporting a problem. The last row is the same
particle count giving the converged answer. (Per-slice cost dominates real
time-dependent runs, where the same particle reduction pays proportionally.)

The guardrails are part of the feature (all refusals by name, manual
sec:coherent-source). A per-slice Gaussianity test is sized against its own
sampling significance. A genuinely structured profile refuses. An offset,
mismatched or tilted Gaussian beam passes -- the source centers and tilts with
the beam's phasor-weighted moments. Harmonics, two polarizations and the
unaveraged mode are out of scope in v1. Dark starts are refused, at a measured
~175x startup deficit: SASE grows from spontaneous, spatially-incoherent
emission, which is exactly what the coherent model drops. Seed the field, as
here, or use the default deposit.

Plot any run with ../plot_fel.py <root>.stats.h5.
