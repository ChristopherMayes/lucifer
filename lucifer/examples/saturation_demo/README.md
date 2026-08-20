# The saturation demo: one practical SASE case, three trackers, one clock

Genesis's own Benchmark1-SASE configuration run to saturation -- the full 57 m 6-FODO
Aramis line, dark start, growth from shot noise alone, 96 slices x 2048 particles --
tracked three ways from IDENTICAL initial dumps, each on the machine's full
performance-core count. Every input is a real file in this directory:

| file | what it is |
|---|---|
| `Aramis.lat` | the Genesis lattice (6 FODO cells, 12 undulator segments) |
| `aramis.bmad` | the Bmad translation (real wigglers; manual sec:element) |
| `sat_unavg.bmad` | two-line wrapper selecting `fel_tracking = fel_unaveraged` |
| `sat-prep.in` | Genesis deck that writes the shared initial dumps (no tracking) |
| `sat-genesis.in` | the timed Genesis deck (same seed and ranks as the prep) |
| `sat-avg.nml` | Bmad averaged mode (the `bmad_standard` default) |
| `sat-unavg.nml` | Bmad unaveraged mode (~32x cost; the point is the same answer from raw dynamics) |
| `check_agreement.py` | exit powers must agree at the documented levels before timings mean anything |
| `run.sh` | thin runner: prep, three timed runs (one external clock), check, figure |

```
./run.sh                     # outputs land in ./output, figure output/sat-demo.png
```

Measured results, attribution of every difference (the averaged mode matches Genesis
at 4.9e-4 through saturation; the unaveraged mode rides +0.6 ln on the shot-noise
radiation channel it physically resolves, FINDINGS 7.27), and the banked figure live
in the "The saturation demo" section of `bsim/fel/README.md` and
`fel-benchmark-plots/sat-demo.png` at the project root.
