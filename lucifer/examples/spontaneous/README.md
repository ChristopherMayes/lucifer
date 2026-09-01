# Spontaneous emission on Bmad's own switches

Two commands, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    ../../../production/bin/lucifer lucifer_damping_only.in
    python ../plot_fel.py spontaneous.stats.h5

Undulator radiation the FEL model does not account for, applied inside FEL elements
under Bmad's own global switches. `bmad_com%radiation_damping_on` subtracts the
classical rate from every particle, and `bmad_com%radiation_fluctuations_on` adds the
quantum kick, drawn once per beamlet because the quiet start cancels within a beamlet.
There is no separate FEL-side switch, and there is no separate FEL-side default:
these are the same two switches Tao exposes, off unless the input turns them on
([the input reference](../../doc/input-reference.md)).

Both decks are the `../steady_state` deck plus the switches, so that run is the
control and the difference is the radiation.

The plot shows it in the beam-energy panel, where the mean falls in a straight line
from the first meter. Over the first 20 m the FEL has amplified the seed to only a few
MW and has extracted nothing measurable, so that slope is the radiation and nothing
else. Fitting it there is the measurement:

| run | d\<E\>/ds per metre of undulator | against the analytic rate |
|---|---|---|
| `lucifer_damping_only.in` | -0.03071 m_e c^2 | 100.3% |
| `lucifer.in`, both switches | -0.03103 m_e c^2 | 101.4% |
| `../steady_state`, neither switch | -0.00008 m_e c^2 | 0.25% |

The classical rate `(2/3) r_e gamma^2 k_u^2 a_w^2` is 0.030616 m_e c^2 per metre of
undulator at these parameters, and the damping run reproduces it to 0.3%. The last row
is the control: with both switches off the same fit finds a quarter of a percent of that
rate, which is the FEL's own start-up extraction rather than radiation. Over the whole
line the analytic loss is 1.47 m_e c^2 across 47.88 m of undulator, and by then the FEL
is extracting far more than radiation is, so the pre-gain slope is where the two can be
told apart.

The harness makes the same comparison without any FEL model, tracking the wiggler with
Bmad's own `runge_kutta`, and reproduces the analytic rate to 1.0e-4
([validation](../../doc/validation.md), the spontaneous-emission section). That check
also prices what the averaged FEL mode does not do: its field gains `2S` while the
particles are kicked by `E`, so the `|S|^2` part of the field energy is created rather
than taken from the beam. Genesis 1.3 Version 4 (Genesis4) carries an optional
`&sponrad` module for the same reason, off by default.

Fluctuations barely move the exit numbers here because at 1 Angstrom the quantum kick
is small against the 1 m_e c^2 the beam starts with and the 4.5 m_e c^2 the FEL itself
drives by saturation. They are refused together with slice migration, since one draw
per beamlet and a pass that regroups beamlets cannot both be right.

Runs in ~5 s each.
