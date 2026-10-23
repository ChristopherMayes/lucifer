# BMAD-STATS-EXT-FEL: the FEL extension

This is version 1.0, in force for this writer, with open questions pending. It extends BMAD-STATS-SPEC.md for time-dependent FEL tracking output. It is declared by `'fel'` in the root `@extensions` and by `@fel_version` on the root. Everything here satisfies the core contract (spec R41): a generic core reader loads an FEL file without knowing this document exists. Rules are numbered F1, F2, ... and cite the core rules they instantiate.

Lucifer is the reference writer. The escaped-field and pulse field files are companion files and not part of this extension. They carry Genesis 1.3 Version 4 (Genesis4) field conventions for reasons recorded in Lucifer's manual, and nothing in a bmad-stats file points into them normatively.

## 1. The slice axis

- **F1.** The time window is a `slice` axis (core R16): `coords/slice = 0..n-1`, the slice number, with three positions as variables on it. The index is the axis because the positions are derived and one of them is reference-dependent.
- **F2.** `coords/ct_slice = i * slice_spacing` (0-based) is the light-travel distance of each slice ahead of the reference. It is exact and free of beta, and it is the coordinate slippage counts in: one slice is exactly `window_sample` radiation wavelengths of slippage, which is what lets the field record rotate by one whole index with no interpolation.
- **F3.** `coords/t_slice = -ct_slice / c` is the arrival time relative to the reference. It is exact: the slice grid is uniform in time, because the migration invariant carries each particle's own beta and `z = -beta c (t - t_ref)` cancels it at a grid point. A writer whose `t_slice` carries any beta does not conform.
- **F4.** `coords/z_slice = beta0 * ct_slice` is Bmad's z at the reference beta, with `beta0` derived from `run/p0c`. Its description MUST say it is reference-dependent: a particle's own offset uses its own beta.
- **F5.** All three position variables MUST carry `@head_direction`, value `'+index'` or `'-index'`, stating which end of the slice index is the window head. It is stated on all three, because a convention stated on one of a set invites the reader to assume the set agrees.
- **F6.** `coords/slice` SHOULD carry `@plot_against` (core R17) naming one of the three.

## 2. The record axis, FEL reading

- **F7.** The core `record` axis carries the per-record beam and field data. `coords/s` (path length, Bmad's s) may legitimately repeat at a zero-length element that acts, a wake kick or a thin kicker, which is the case core R16 exists for. A repeat inside one element is a writer defect, and the writer's harness checks it.

## 3. Groups

- **F8.** `beam/slice/`, kind `per_slice`, holds the per-record sufficient statistics, at every record, over `(record, slice)` leading axes. Its datasets are `centroid`, `sigma` (over `bmad`, `bmad_col`), `charge_live`, `n_particle_live`, `t`, `sigma_t`, `bunching`, `bunching_phase`, `current`, `energy` (eV, never gamma), `sigma_energy`, `rel_max`, `rel_min`. The moments, weights, bunching and extremes are the primary layer (core R43). `current`, `energy` and `sigma_energy` are derived conveniences and declare `@derived_from` (core R38). `rel_max` and `rel_min` are the per-coordinate extremes relative to the centroid, Bmad's own convention, over the seven-entry label axis `bmad_t` (the six phase-space names plus `t`). They are order statistics, primary by definition since no moment reconstructs them, and what an envelope plot needs (envelope = centroid + rel). An empty slice's extremes are NaN (core R28).
- **F9.** `beam/slice_twiss/`, kind `per_slice`, holds the evaluated bunch_params per slice on the `element_end` axis, with `twiss/` over `plane` and `modes/` over `mode` (core R32). Its `centroid` and `sigma` are the masked records' own moments, and `twiss/` and `modes/` derive from them (core R38, R44): the mode labels are eigenvector-identified, so a re-deriver matches them as a set. `beam/bunch/` is itself derived, pooled from `beam/slice` by the covariance identity, and says so.
- **F10.** `field/`, kind `field`, holds the radiation, per wavelength and polarization. It MUST declare `@components` (core R39) and, when harmonics are present, `@harmonics`, an integer array.
- **F11.** Each polarization component has a group of kind `component` (`field/x`, `field/y`), holding `power`, `energy`, `on_axis_intensity`, `centroid` and `sigma` over the `wavefront` and `wavefront_col` label axes (`x, theta_x, y, theta_y`, with `coords/wavefront_unit`), `emit_x`, `emit_y`, and `angle_moments_valid`. The last is a core R25 boolean, since the theta moments cost FFTs and are NaN where not computed (core R28). `emit_x` and `emit_y` are pure functions of `sigma` and declare `@derived_from` (core R38).
- **F12.** `field/total/`, kind `derived` with `@derived_from` (core R38), is the sum over the live polarizations of one wavelength. It is always written, whether one component is live or two. So no reader asks what else the file holds before it knows what power means (core R29).
- **F12a.** `field/total/power_inside_angle` and `field/total/split_angle` are the power of each slice within `split_angle` of the axis, summed over both polarization planes, and that angle. The angle is a run-level scalar, the protected angle the filter's edge is placed against, whether or not the filter is on. The power rides the transform the theta moments take, so it is NaN at every record where they were not (core R28). Its complement to `power` is the wide-angle emission of the point macroparticle beamlets.
- **F13.** Each harmonic has a group of kind `harmonic` (`field/harm<h>`), carrying `@harmonic`, a true scalar integer (core R24), and its own components and `total/` under the same dataset names. Nothing is ever summed across harmonics: a detector separates colors.
- **F16.** `field/<component>/reduced/`, kind `projected`, is OPTIONAL and holds intensity projections of that component, each summed over the axis it drops, and the field on axis. Its datasets are `xy_intensity` over `(record, grid_y, grid_x)`, `slice_x_intensity` over `(record, slice, grid_x)`, `slice_y_intensity` over `(record, slice, grid_y)`, and `on_axis_field_re` and `on_axis_field_im` over `(record, slice)` for the field on axis. Intensity is `|E|^2 / (2 mu_0 c)`, the same quantity as `on_axis_intensity` beside it, so the projections integrate to `power`. A conforming writer's `xy_intensity` summed over the transverse grid times the cell area equals that record's total `power`, and `slice_x_intensity` does the same per slice. The group exists so that a reader can draw a run without opening a field dump, and it is optional because the transverse projection is the grid's own size at every record. `coords/grid_x` and `coords/grid_y`, the transverse positions in metres measured from the axis, are present with it.

## 4. Parameters

- **F14.** `params/` gains the program's slicing and seed structs (for Lucifer: `lambda0`, `n_wavelength`, `beamlet_size`, the grid, the seed knobs), under core R34.
- **F15.** `run/` gains `slice_spacing` and `n_slice`. `slice_spacing`'s description MUST state that it is a light-travel distance, c times the slice time separation, with the z separation being beta times it, per particle. It is derived (`window_sample * lambda0`), which is why it is in `run/` and not `params/`.

## 5. Kinds and axes this extension adds

The kinds are `per_slice`, `field`, `component`, `harmonic`, `beam`, `projected`, `twiss`, `modes`. The axes are `slice`, `wavefront`, `wavefront_col`, `bmad_t`, `grid_x`, `grid_y`. The attributes are `@head_direction`, `@harmonic`, `@harmonics`, `@fel_version`.

## 6. Writer-harness identities

These identities are not checkable by the generic validator and are required of a conforming writer's own tests. `t_slice == -ct_slice/c` holds exactly, and so does `z_slice == beta0 * ct_slice`. The `reduced/` projections integrate to `power` per record and per slice where they are written. `field/total` equals the sum of its `@derived_from`. The angle-moment validity flags match where the FFTs actually ran. The position entries of `rel_max` and `rel_min` match extremes over a particle dump at the same plane. Every repeated `coords/s` straddles an element boundary.

## 7. Open questions

- **OF1.** `emit_x`/`emit_y` could be two datasets or one `emit` over a two-entry axis. They are kept as two on the grounds that a two-entry axis buys nothing. A standard may want the uniformity anyway.
- **OF2.** The companion field files (escaped, pulse) may eventually want an extension of their own rather than Genesis4 conventions. That question is tabled with the escaped-field openPMD question.
