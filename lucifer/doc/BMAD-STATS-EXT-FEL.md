# BMAD-STATS-EXT-FEL: the FEL extension

Version 1.0, in force for this writer, with open questions pending. Extends BMAD-STATS-SPEC.md for time-dependent FEL tracking output. Declared by `'fel'` in the root `@extensions` and `@fel_version` on the root. Everything here satisfies the core contract (spec R41): a generic core reader loads an FEL file without knowing this document exists. Rules are numbered F1, F2, ... and cite the core rules they instantiate.

Lucifer is the reference writer. The escaped-field and pulse field files are companion files, not part of this extension: they carry Genesis field conventions for reasons recorded in Lucifer's manual, and nothing in a bmad-stats file points into them normatively.

## 1. The slice axis

- **F1.** The time window is a `slice` axis (core R16): `coords/slice = 0..n-1`, the slice number, with three positions as variables on it. The index is the axis because the positions are derived and one of them is reference-dependent.
- **F2.** `coords/ct_slice = i * slice_spacing` (0-based), the light-travel distance of each slice ahead of the reference. Exact and free of beta, and the coordinate slippage counts in: one slice is exactly `window_sample` radiation wavelengths of slippage, which is what lets the field record rotate by one whole index with no interpolation.
- **F3.** `coords/t_slice = -ct_slice / c`, the arrival time relative to the reference. Exact: the slice grid is uniform in time, because the migration invariant carries each particle's own beta and `z = -beta c (t - t_ref)` cancels it at a grid point. A writer whose `t_slice` carries any beta does not conform.
- **F4.** `coords/z_slice = beta0 * ct_slice`, Bmad's z at the reference beta, with `beta0` derived from `run/p0c`. Its description MUST say it is reference-dependent: a particle's own offset uses its own beta.
- **F5.** All three position variables MUST carry `@head_direction`, value `'+index'` or `'-index'`, stating which end of the slice index is the window head. Stated on all three because a convention stated on one of a set invites the reader to assume the set agrees.
- **F6.** `coords/slice` SHOULD carry `@plot_against` (core R17) naming one of the three.

## 2. The record axis, FEL reading

- **F7.** The core `record` axis carries the per-record beam and field data. `coords/s` (path length, Bmad's s) may legitimately repeat at a zero-length element that applies a wake kick, which is the case core R16 exists for. A repeat inside one element is a writer defect, and the writer's harness checks it.

## 3. Groups

- **F8.** `beam/slice/`, kind `per_slice`: the per-record sufficient statistics, at every record, over `(record, slice)` leading axes. Datasets: `centroid`, `sigma` (over `bmad`, `bmad_col`), `charge_live`, `n_particle_live`, `t`, `sigma_t`, `bunching`, `bunching_phase`, `current`, `energy` (eV, never gamma), `sigma_energy`, `rel_max`, `rel_min`. The moments, weights, bunching and extremes are the primary layer (core R43). `current`, `energy` and `sigma_energy` are derived conveniences and declare `@derived_from` (core R38). `rel_max` and `rel_min` are the per-coordinate extremes relative to the centroid, Bmad's own convention, over the seven-entry label axis `bmad_t` (the six phase-space names plus `t`): order statistics, primary by definition since no moment reconstructs them, and what an envelope plot needs (envelope = centroid + rel). An empty slice's extremes are NaN (core R28).
- **F9.** `beam/slice_twiss/`, kind `per_slice`: the evaluated bunch_params per slice on the `element_end` axis, with `twiss/` over `plane` and `modes/` over `mode` (core R32). Its `centroid` and `sigma` are the masked records' own moments, and `twiss/` and `modes/` derive from them (core R38, R44): the mode labels are eigenvector-identified, so a re-deriver matches them as a set. `beam/bunch/` is itself derived, pooled from `beam/slice` by the covariance identity, and says so.
- **F10.** `field/`, kind `field`: the radiation, per wavelength and polarization. It MUST declare `@components` (core R39) and, when harmonics are present, `@harmonics`, an integer array.
- **F11.** One group per polarization component, kind `component` (`field/x`, `field/y`), each holding: `power`, `energy`, `on_axis_intensity`, `centroid` and `sigma` over the `wavefront` and `wavefront_col` label axes (`x, theta_x, y, theta_y`, with `coords/wavefront_unit`), `emit_x`, `emit_y`, and `angle_moments_valid` (a core R25 boolean, since the theta moments cost FFTs and are NaN where not computed, core R28). `emit_x` and `emit_y` are pure functions of `sigma` and declare `@derived_from` (core R38).
- **F12.** `field/total/`, kind `derived` with `@derived_from` (core R38): the sum over the live polarizations of one wavelength, always written, whether one component is live or two, so no reader asks what else the file holds before it knows what power means (core R29).
- **F13.** One group per harmonic, kind `harmonic` (`field/harm<h>`), carrying `@harmonic`, a true scalar integer (core R24), and its own components and `total/` under the same dataset names. Nothing is ever summed across harmonics: a detector separates colors.

## 4. Parameters

- **F14.** `params/` gains the program's window and seed structs (for Lucifer: `lambda0`, `window_sample`, `nbins`, the grid, the seed knobs), under core R34.
- **F15.** `run/` gains `slice_spacing` and `n_slice`. `slice_spacing`'s description MUST state that it is a light-travel distance, c times the slice time separation, with the z separation being beta times it, per particle. It is derived (`window_sample * lambda0`), which is why it is in `run/` and not `params/`.

## 5. Kinds and axes this extension adds

Kinds: `per_slice`, `field`, `component`, `harmonic`, `beam`, `projected`, `twiss`, `modes`. Axes: `slice`, `wavefront`, `wavefront_col`, `bmad_t`. Attributes: `@head_direction`, `@harmonic`, `@harmonics`, `@fel_version`.

## 6. Writer-harness identities

Not checkable by the generic validator, required of a conforming writer's own tests: `t_slice == -ct_slice/c` exactly, `z_slice == beta0 * ct_slice`, `field/total` equal to the sum of its `@derived_from`, the angle-moment validity flags matching where the FFTs actually ran, the position entries of `rel_max` and `rel_min` matching extremes over a particle dump at the same plane, and every repeated `coords/s` straddling an element boundary.

## 7. Open questions

- **OF1.** `emit_x`/`emit_y` as two datasets against one `emit` over a two-entry axis. Kept as two on the grounds that a two-entry axis buys nothing; a standard may want the uniformity anyway.
- **OF2.** Whether the companion field files (escaped, pulse) should eventually get an extension of their own rather than Genesis conventions, tabled with the escaped-field openPMD question.
