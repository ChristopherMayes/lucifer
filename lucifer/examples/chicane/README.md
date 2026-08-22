# A chicane between segments: phasing from real geometry

Two commands, Bmad only:

    ../../../../production/bin/lucifer run.nml
    ../../../../production/bin/lucifer run_detuned.nml

Two undulator segments with a four-bend closed-bump chicane in the break, in
ABSOLUTE TIME TRACKING mode (`bmad_com[absolute_time_tracking] = T` in the
lattice): the inter-segment phase follows the real geometry. The beam detours the
bump while the radiation drifts the chord between the undulator faces (computed
from `ele%floor`, never entered by hand); the arc-minus-chord delay -- 583
wavelengths at the reference angle -- re-injects the light onto the bunched beam
at whatever carrier phase the bend angle dictates.

The two lattices differ by 0.43 microradians of bend angle: HALF A WAVELENGTH of
delay. Measured (200 MW seed, two 0.99 m helical segments):

| lattice | P at segment-1 exit | P at line exit | segment-2 ratio |
|---|---|---|---|
| `chicane.bmad` (in phase) | 1.998e8 W | 2.642e8 W | 1.32 |
| `detuned.bmad` (+lambda/2 of delay) | 1.998e8 W | 1.872e8 W | 0.94 |

The second segment amplifies or absorbs on a half-wavelength of geometry -- the
re-phasing knob a real machine turns by trimming its chicane. In the default
RELATIVE mode the same two lattices produce identical output (the geometric
fraction is autophased away, Genesis's chicane semantics); the deliberate
off-phase knob there is the wiggler's own `z_offset` (manual sec:phasing).

Validated by the harness's phasing section: the re-anchor baseline flat, the
z_offset knob on the analytic slope, the knob curve equal to Genesis's own
PHASESHIFTER scan at 6.0e-6, the absolute-mode chicane ramp on the independent
geometric prediction at 6.8e-4, the unaveraged ledger closing across the chicane
at 4.0e-6, and four refusals by name. Runs in ~1 min each.
