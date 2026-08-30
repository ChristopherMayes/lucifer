# The saturation demo: one practical SASE case, three trackers, one clock

Genesis's own Benchmark1-SASE configuration run to saturation (the full 57 m 6-FODO
Aramis line, dark start, growth from shot noise alone, 96 slices x 2048 particles),
tracked three ways from identical initial dumps, each on the machine's full
performance-core count. Every input is a real file in this directory:

| file | what it is |
|---|---|
| `Aramis.lat` | the Genesis lattice (6 FODO cells, 12 undulator segments) |
| `aramis.bmad` | the Bmad translation (real wigglers, manual sec:element) |
| `sat_unavg.bmad` | two-line wrapper selecting `fel_tracking = fel_unaveraged` |
| `sat-prep.in` | Genesis deck that writes the shared initial dumps (no tracking) |
| `sat-genesis.in` | the timed Genesis deck (same seed and ranks as the prep) |
| `sat-avg.nml` | Bmad averaged mode (the `bmad_standard` default) |
| `sat-unavg.nml` | Bmad unaveraged mode (~32x cost). The point is the same answer from raw dynamics |
| `check_agreement.py` | exit powers must agree at the documented levels before timings mean anything |
| `run.sh` | thin runner: prep, three timed runs (one external clock), check, PDF report |

```
./run.sh                     # outputs land in ./output; summary: output/sat-demo-report.pdf
```

The report (`tests/scripts/report_fel_saturation.py`) is a multi-page PDF regenerated
from the run's own files: a cover with the timing and agreement tables, then gain
curves, pulse structure, beam evolution and the energy accounting -- each figure with
the paragraph that explains how to read it. The energy page shows, per tracker, the
energy the beam gave against where it is now (window vs slipped-out-forward), and the
unaveraged ledger closing exactly along z.

## The measured result

Measured (M3 Max, 12 performance cores, production builds both sides):

| | wall | exit total power | vs Genesis |
|---|---|---|---|
| Genesis 1.3 v4, 12 MPI ranks | 38.0 s | 3.381 GW (4.52 GW peak at 56.8 m) | the reference |
| Bmad averaged (`bmad_standard` default), 12 threads | 30.2 s | 3.380 GW | **rel 4.9e-4** |
| Bmad unaveraged (`fel_tracking = fel_unaveraged`), 12 threads | 1143.2 s | 6.25 GW | ln ratio +0.62 |

The averaged mode tracks Genesis through eight decades of z and three of power to
**4.9e-4** at saturation. The ~4e-2 seam-transport difference the benchmark tiers
price is invisible here because saturation self-limits the power. It is also 1.26x
faster than Genesis at equal cores, each code computing its own in-run diagnostics,
and ours are the full 6x6/4x4 moment sets of the stats file where Genesis's are
scalar columns (the diagnostics record in `lucifer/doc/validation.md` is where that speed
came from). Spontaneous radiation is off in the demo on both sides, matching Genesis's
`&sponrad` default. Turning `radiation_damping`/`radiation_fluctuations` on costs the
beam ~1.3e-4 of its energy over the line (~0.5 rho of accumulated detuning -- the size
that moves hard-X-ray saturation, and the reason the switches exist).

The unaveraged mode is an independent integrator with fc/JJ nowhere in its inputs,
paying its documented ~32x cost. It reproduces the startup (coherent shot-noise
radiation matches both codes to ~8%), the gain curve shape, and the saturation
location (56.2 vs 56.8 m), and rides ~2%/m above the KMR codes through the
exponential regime, while its beam gives up ~14x more energy. The energy difference is
understood and measured: both models emit spontaneous shot-noise
radiation of the same magnitude, but the averaged/KMR model does not debit the beam
for it (its step adds 2S to the field and kicks with E, so the 4|S|^2 part of the
field energy is created, measured factor 134 on a dark segment), while the unaveraged
mode conserves energy by construction and therefore pays. The unaveraged mode's
captured spontaneous loss agrees with the analytic rate (2/3)r_e gamma^2 ku^2 aw^2
restricted to the grid's angular acceptance to 8%, the same formula Genesis's own
optional &sponrad module uses. Neither model yet carries the ~90% of spontaneous
power radiated outside the grid acceptance (the named follow-on). The energy panels
are where all of this is visible (the budget panel:
Genesis and the averaged mode both keep ~72% of the beam's energy in the window at
exit, the unaveraged mode keeps 10%. Its beam gave 6.6e-8 J against their 4.7e-9,
the difference radiated and slipped out forward. Also the faster energy-spread
growth), and "slipped out forward" is bookkept, not asserted: the unaveraged ledger
banks the energy of every slice the slippage zero-fill discards (U_escaped) and the
deposit's own |src|^2 (U_spont, the one term the kick/deposit duality does not charge
to the beam -- physically the substep's spontaneous emission), so the time-dependent
books close exactly: E_beam + U_window + U_escaped - U_spont conserved to 2.9e-3 of
turnover over the whole demo (8.0e-6 on the harness configuration, where it is a
standing check at 1e-3. The demo's 6.6e-8 J = 7.9e-9 held + 5.9e-8 escaped - 1.1e-9
spontaneous credit). The figure's dotted curve is that closure drawn on the budget
panel. Wakes would be a second, unbanked exit channel from the beam -- and the ledger
exists only in the unaveraged mode, where wakes and space charge are refused by name. Measured independent of
particle count (1024/2048/4096) and steps-per-period (20/40): physics, not statistics
or resolution. A dark segment with real shot noise isolates it: same in-window noise
power in both models, 20x the beam-side energy cost in the unaveraged one. The absolute
calibration of that channel (the in-band, in-grid-acceptance fraction of undulator
radiation) is the named future check. Until then the demo prices the difference at
|ln| <= 1.0 at exit, measured 0.615.

Two more unaveraged end effects were found and dispatched on the way to this figure:
the ramps' slippage deficit (3.3 rad of optical phase per end, compensated exactly by
the built-in handoff phase jump, and tier1_unavg's theta median fell
from 6.6 rad to 6.4e-2) and their reduced coupling length (~2% ln per segment at the
default 2-period ramps, real field physics, priced and left visible).

The report is banked at `fel-benchmark-plots/sat-demo-report.pdf` in the project root.

