# The saturation demo: one practical SASE case, three trackers, one clock

This is the Benchmark1-SASE configuration of Genesis 1.3 Version 4 (Genesis4), run to saturation: the full 57 m 6-FODO Aramis line, a dark start, growth from shot noise alone, 96 slices x 2048 particles. It is tracked three ways from identical initial dumps, each on the machine's full performance-core count. Every input is a file in this directory:

| file | what it is |
|---|---|
| `Aramis.lat` | the Genesis4 lattice (6 FODO cells, 12 undulator segments) |
| `aramis.bmad` | the Bmad translation (real wigglers, the manual's FEL-element section) |
| `sat_unavg.bmad` | two-line wrapper selecting `fel_method = unaveraged` |
| `sat-prep.in` | Genesis4 deck that writes the shared initial dumps (no tracking) |
| `sat-genesis.in` | the timed Genesis4 deck (same seed and ranks as the prep) |
| `lucifer-avg.in` | Bmad averaged mode (the `bmad_standard` default) |
| `lucifer-unavg.in` | Bmad unaveraged mode (~32x cost). The point is the same answer from raw dynamics |
| `check_agreement.py` | exit powers must agree at the documented levels before timings mean anything |
| `run.sh` | thin runner: prep, three timed runs (one external clock), check, PDF report |

```
./run.sh                     # outputs land in ./output; summary: output/sat-demo-report.pdf
```

The report (`tests/scripts/report_fel_saturation.py`) is a multi-page PDF built from the run's own files: a cover with the timing and agreement tables, then the gain curves, the pulse structure, the beam evolution and the energy accounting. Each figure carries the paragraph that says how to read it. The energy page shows, per tracker, the energy the beam gave against where it is now, in the window or slipped out forward, and the unaveraged ledger closing along z.

## The measured result

The grid here is `ngrid` 255, and odd where the rest of the examples run powers of two. This demo's headline is an agreement with Genesis4, and Genesis4 propagates the transverse field on half-integer wavenumbers when `ngrid` is even, so a comparison this tight has to run where the reference is right (`lucifer/doc/validation.md`). The grid moves to 256 once that is fixed upstream.

On an M3 Max with 12 performance cores and production builds on both sides, the three trackers give:

| | wall | exit total power | vs Genesis4 |
|---|---|---|---|
| Genesis4, 12 MPI ranks | 38.0 s | 3.381 GW (4.52 GW peak at 56.8 m) | the reference |
| Bmad averaged (`bmad_standard` default), 12 threads | 30.2 s | 3.380 GW | rel 4.9e-4 |
| Bmad unaveraged (`fel_method = unaveraged`), 12 threads | 1143.2 s | 6.25 GW | ln ratio +0.62 |

The averaged mode tracks Genesis4 through eight decades of z and three of power to 4.9e-4 at saturation. The seam-transport difference of about 4e-2 that the benchmark tiers measure is invisible here, because saturation self-limits the power. On this configuration and machine the averaged mode also ran 1.26x faster than Genesis4 at equal cores. Each code computes its own in-run diagnostics, and Lucifer's are the full 6x6/4x4 moment sets of the stats file where Genesis4's are scalar columns (the diagnostics record in `lucifer/doc/validation.md`). Spontaneous radiation is off in the demo on both sides, matching Genesis4's `&sponrad` default. Turning `radiation_damping`/`radiation_fluctuations` on costs the beam about 1.3e-4 of its energy over the line, about 0.5 rho of accumulated detuning. That is the size that moves hard-X-ray saturation, and the reason the switches exist.

The unaveraged mode is an independent integrator with fc/JJ nowhere in its inputs, at its documented cost of about 32x. It reproduces the startup, where the coherent shot-noise radiation matches both codes to about 8%, the shape of the gain curve, and the saturation location, 56.2 against 56.8 m. It rides about 2%/m above the KMR codes through the exponential regime, and its beam gives up about 14x more energy.

The energy difference is measured. Both models emit spontaneous shot-noise radiation of the same magnitude. The averaged (KMR) model does not debit the beam for it. Its step adds 2S to the field and kicks with E, so the 4|S|^2 part of the field energy is created, a measured factor of 134 on a dark segment. The unaveraged mode conserves energy by construction and therefore pays. Its captured spontaneous loss agrees to 8% with the analytic rate (2/3)r_e gamma^2 ku^2 aw^2 restricted to the grid's angular acceptance, the formula Genesis4's own optional &sponrad module uses. Neither model carries the roughly 90% of spontaneous power radiated outside the grid acceptance, which is a named follow-on.

The energy panels show all of this. In the budget panel Genesis4 and the averaged mode both keep about 72% of the beam's energy in the window at the exit, and the unaveraged mode keeps 10%. Its beam gave 6.6e-8 J against their 4.7e-9, and the difference radiated and slipped out forward. The energy-spread panel shows its faster growth. The slipped-out energy is bookkept. The unaveraged ledger banks the energy of every slice the slippage zero-fill discards (U_escaped) and the deposit's own |src|^2 (U_spont, the one term the kick/deposit duality does not charge to the beam, physically the substep's spontaneous emission). The time-dependent books then close: E_beam + U_window + U_escaped - U_spont is conserved to 2.9e-3 of turnover over the whole demo, and to 8.0e-6 on the harness configuration, where the check level is 1e-3. The demo's 6.6e-8 J = 7.9e-9 held + 5.9e-8 escaped - 1.1e-9 spontaneous credit. The figure's dotted curve is that closure drawn on the budget panel.

A wake would be a second exit channel from the beam that the ledger does not bank. The ledger exists only in the unaveraged mode, where wakes and space charge are refused. The energy difference is independent of the particle count (1024/2048/4096) and of the steps per period (20/40), so it is physics and neither statistics nor resolution. A dark segment with real shot noise isolates it: the same in-window noise power in both models, and 20x the beam-side energy cost in the unaveraged one. The absolute calibration of that channel, the in-band and in-grid-acceptance fraction of undulator radiation, is a named future check. The demo holds the difference at |ln| <= 1.0 at the exit, measured 0.615.

The unaveraged ramps carry two end effects. Their slippage deficit, 3.3 rad of optical phase per end, is compensated exactly by the built-in handoff phase jump, which took the theta median of the unaveraged steady-state tier (tier1_unavg) from 6.6 rad to 6.4e-2. Their reduced coupling length costs about 2% in ln power per segment at the default 2-period ramps. That is the field physics of an adiabatic end, measured and left visible.

The runner writes the report to `output/sat-demo-report.pdf`.
