# A chicane between segments: phasing from real geometry

The example runs in two commands, with Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_detuned.in

The line is two undulator segments with a four-bend closed-bump chicane in the break, in absolute time tracking mode (`bmad_com[absolute_time_tracking] = T` in the lattice). The inter-segment phase then follows the real geometry. The beam detours the bump while the radiation drifts the chord between the undulator faces, which is computed from `ele%floor` and never entered by hand. The arc-minus-chord delay (583 wavelengths at the reference angle) re-injects the light onto the bunched beam at whatever carrier phase the bend angle dictates.

The two lattices differ by 0.43 microradians of bend angle: half a wavelength of delay. With a 200 MW seed and two 0.99 m helical segments, the two lattices give:

| lattice | P at segment-1 exit | P at line exit | segment-2 ratio |
|---|---|---|---|
| `chicane.bmad` (in phase) | 1.998e8 W | 2.642e8 W | 1.32 |
| `detuned.bmad` (+lambda/2 of delay) | 1.998e8 W | 1.872e8 W | 0.94 |

The second segment amplifies or absorbs on a half-wavelength of geometry. A real machine turns that knob by trimming its chicane. In the default relative mode the same two lattices produce identical output, since the geometric fraction is autophased away, which is the chicane semantics of Genesis 1.3 Version 4 (Genesis4). The deliberate off-phase knob in that mode is the wiggler's own `z_offset` (the manual's phasing section).

The harness's phasing section checks the same physics. The re-anchor baseline is flat, and the z_offset knob sits on the analytic slope. The knob curve equals Genesis4's own PHASESHIFTER scan at 6.0e-6. The absolute-mode chicane ramp sits on the independent geometric prediction at 6.8e-4, and the unaveraged ledger closes across the chicane at 4.0e-6. Four refusals are asserted.

Each run takes ~1 s.
