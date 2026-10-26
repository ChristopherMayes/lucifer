# The coherent source: converged results at a fraction of the particles

The example runs in three commands, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    ../../build-release/bin/lucifer lucifer_low_m.toml
    ../../build-release/bin/lucifer lucifer_reference.toml
    python -m lucifer.plot coherent.stats.h5

A run whose beam is transversely Gaussian (an idealized machine, a parameter scan, a seeded amplifier study) can trade the per-particle source deposit for Tanaka's coherent retrieval, `run.source_model = "coherent"` (PRAB 27, 030703 (2024), the manual's coherent-source section). It then converges at one to two orders of magnitude fewer macroparticles. With few particles the per-cell source is spiky, and that sampling noise always overestimates the gain. The coherent source carries the slice's exact bunching phasor on a fitted Gaussian instead, so the artifact never enters.

Three runs of the same seeded steady-state case (the full Aramis line, 5.2 decades of gain, ../steady_state's configuration) are compared. `lucifer_reference.toml` deposits 8192 particles per slice. `lucifer_low_m.toml` deposits 512. `lucifer.toml` runs the coherent source on the same 512. The measured quantities are in the table below. The low-M deposit comes out high against the reference, and the coherent source at the same count lands close to it. Each run takes a few seconds.

The low-M deposit is the artifact the coherent source removes, and the source filter, on by default, removes most of it on its own. Cutting particles without the coherent source raises the predicted power on this case, and with the filter off the rise was far larger. The run's own convergence report gives the power outside the split angle over the power inside it, and the table carries that ratio at the exit for each run. The coherent source builds no transformed source for the filter to act on, so the filter turns itself off for that run and says so. Its own ratio at the exit is the largest of the three ([the convergence section](../../docs/fel-physics.md#sec-convergence)). The power the low-M deposit adds is the wide-angle emission of 64 beamlets. The coherent source does not have it, since it never puts a macroparticle on a grid point. `lucifer.toml` is the same particle count giving the converged answer. Per-slice cost dominates a time-dependent run, where the same particle reduction pays in proportion.

The refusals are part of the feature (the manual's coherent-source section). A per-slice Gaussianity test is sized against its own sampling significance, and a profile with real structure is refused. An offset, mismatched or tilted Gaussian beam passes, since the source centers and tilts with the beam's phasor-weighted moments. Harmonic fields, two live polarizations and the unaveraged mode are refused in v1. A dark start is refused, since the coherent model drops the spontaneous, spatially incoherent emission that SASE grows from, and its startup comes out far too low. Seed the field, as here, or use the default deposit.

Plot any run with python -m lucifer.plot <root>.stats.h5.
