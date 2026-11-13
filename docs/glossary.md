---
title: Glossary
short_title: Glossary
---

The words this documentation uses for the pieces of a run and for the checks on it, each with the page that defines it. A page links a word here where it first uses it. Four words carry two meanings each, and both meanings are listed. A `record` is a row of the statistics file in one place and the field's stored array in another. A `plane` is a transverse field sample in one place and a position along the lattice in another. A `kick` is the phase kick of the shot-noise load in one place and the radiation's force on a particle in another. A `ramp` is an undulator's end field in one place and the spectral propagator's slippage phase in another.

## The beam and the time window

:::{figure} figures/window.svg
:name: fig-window
:alt: A row of sixteen slice boxes in window order, the bunch over the low-index slices, the slippage headroom at the head, a slab of four slices, and one slice enlarged to show its wavelengths, its beamlets and its transverse grid.
The window, in window order. Slice 1 is the tail and the highest index the head, which arrives first. The bunch sits low in the window and the slippage headroom sits ahead of it. Light slips toward the head and leaves there into the bank, and vacuum enters at the tail. A slab is the slices a frame writes. One slice, enlarged, is `n_wavelength` wavelengths long and holds its macroparticles in beamlets and one transverse grid of the field. The counts drawn are for the picture: `docs/figures/window_figure.py` draws it.
:::

(gl-wavelength)=
**wavelength.** The radiation wavelength $\lambda_s$, `field.lambda0`. Every longitudinal length below is a whole number of it.

(gl-slice)=
**slice.** One longitudinal sample of the run. The beam's particles are partitioned into slices, and the field holds one transverse grid per slice. Slices are spaced by $\Delta = $ `window.n_wavelength` $\cdot\lambda_s$, so a slice is `n_wavelength` wavelengths long ([](fel-physics.md#sec-window)). Its particles stand over all of those wavelengths ([](fel-physics.md#sec-quiet)). A steady-state run is one slice.

(gl-slice-spacing)=
**slice spacing, `n_wavelength`.** The whole number of wavelengths one slice spans. It is an integer end to end, since slippage moves the field by exactly one slice per `n_wavelength` undulator periods ([](fel-physics.md#sec-slippage)).

(gl-window)=
**window.** The set of slices the run carries, `window.n_slice` of them, fixed in the beam's frame. Slice 1 is the tail and the highest index is the head, which arrives first. Light slips toward the head and leaves the window there. A default window holds every particle and, ahead of the head, the slippage the line applies ([](fel-physics.md#sec-window)). The window is *not* the pulse: see `pulse` below.

(gl-head)=
**head, tail.** The high-index end and the low-index end of the window. Radiation leaves at the head, and fresh vacuum enters at the tail.

(gl-window-order)=
**window order, time order.** Slices in the order of their index, 1 to `n_slice`. Every file is written in this order, and so is every statistics array over slices.

(gl-bunch)=
**bunch, beam.** The particles of every slice together. Bmad's `beam_init` generates the bunch or reads it from a file, and it is binned by arrival time into the slices ([](#gl-load), [](fel-physics.md#sec-loading)).

(gl-macroparticle)=
**macroparticle, weight.** A tracked particle standing for `weight` electrons. The weights are kept per particle, so a slice with little charge has light particles rather than few.

(gl-beamlet)=
**beamlet.** A group of `beamlet_size` macroparticles with one transverse sample and one energy, their phases spaced evenly over $2\pi$ so that every bunching harmonic below the group's size cancels. That is the quiet start ([](#gl-quiet-start)), and the shot noise is imposed per beamlet ([](fel-physics.md#sec-quiet), [](fel-physics.md#sec-noise)). The size is derived from the highest harmonic carried, twice it plus 2, where `beam.beamlet_size` is not set ([](shot-noise.md#shot-size)). A slice has `beam.n_particle / beamlet_size` beamlets in the sample mode. Under `shot_noise` a slice of little charge has fewer, so that none stands for less than `beam.beamlet_min_electrons` electrons ([](shot-noise.md#shot-fewer)).

## The load

(gl-load)=
**load, loader.** Loading makes the macroparticles of every slice from the bunch, before the walk. The load is what results: each slice's macroparticles with their weights. The loader is the part of Lucifer that does it, `fel_import_bunch` in `src/fel_import_mod.f90` with the noise routines of `src/fel_beam_mod.f90`. It bins the bunch by arrival time into the slices. `beam.load_mode = "sample"` then draws `beam.n_particle` macroparticles in every slice, and `"keep"` keeps the bunch's own particles ([](fel-physics.md#sec-loading)).

(gl-resampler)=
**resampler.** The sample mode's loader for a bunch read from a file. It draws each slice's macroparticles from the bunch's particles near that slice, at one count a slice, with weights that follow the current. Its keys are the table `[beam.resample]` ([](fel-physics.md#sec-import)).

(gl-quiet-start)=
**quiet start, quiet load.** A load made of beamlets, `beam.quiet_start = true`, the default. A beamlet's phases are evenly spaced, so a quiet load has no bunching at any harmonic below the beamlet size and radiates nothing coherent until noise is imposed ([](fel-physics.md#sec-quiet)).

(gl-shot-noise)=
**shot noise.** The bunching $N$ electrons at random phases have, $\langle|b|^2\rangle = 1/N$, which is what a SASE run starts from. `beam.shot_noise = true` imposes it on a quiet load, through the weights by default (`fel_weight_noise`) or by the kick of Fawley's algorithm (`fel_fawley_noise`), as `beam.shot_noise_method` says ([](fel-physics.md#sec-noise), [](shot-noise.md)).

(gl-kick-noise)=
**kick, first meaning.** The phase shift the shot-noise load gives each macroparticle of a beamlet, a sum of sines of the phase with a random amplitude and a random phase for each harmonic ([](shot-noise.md#shot-kick)). `beam.shot_noise_method = "first_order"` kicks the phases. The default imposes the noise through the weights and kicks nothing ([](shot-noise.md#shot-weight)).

(gl-quiet-floor)=
**quiet floor.** The largest bunching a load has before the noise is imposed, $\max_h |b(h)|^2 N_\lambda$ over the harmonics below the beamlet size. The loader measures it and refuses `shot_noise` on a load above 0.01, since that load already carries noise ([](fel-physics.md#sec-noise)).

(gl-migration)=
**slice migration.** Moving a particle to the slice its phase has drifted into, `run.migrate`, in `fel_migrate_slices`. It is off by default ([](fel-physics.md#sec-migration)).

## The field

(gl-grid)=
**grid.** The transverse grid of one slice's field, `grid_n_pts` points a side at `grid_half_width` from the axis, derived from the beam when the input file states neither ([](fel-physics.md#sec-field)).

(gl-plane-field)=
**plane, first meaning.** One slice's field: one transverse grid of complex values, one for each polarization. The banked planes of the escaped field are planes in this sense.

(gl-plane-lattice)=
**plane, second meaning.** A position along the lattice, as in the entry plane and the exit plane of an element, or the plane a frame is named for. Inside an unaveraged segment the stored field sits half a substep from the plane the particles stand on, and a frame is carried back onto its own plane for the write ([](fel-physics.md#sec-unaveraged)).

(gl-polarization)=
**polarization, component.** The field's `x` and `y` components. A planar undulator with no tilt radiates one, and a tilted or helical one radiates both ([](fel-physics.md#sec-vector)).

(gl-field-set)=
**field set, harmonic, member.** The fields a run carries, one member per harmonic, the fundamental first. Each member is a wavefront with its own slices and polarizations, and the statistics file has a group per member ([](fel-physics.md#sec-field-set)).

(gl-wavefront)=
**wavefront.** The field structure: the planes of one member over the slices, with the grid spacing, the slice spacing and the wavelength. It is what the openPMD field files hold ([](fel-physics.md#sec-openpmd)).

(gl-record-field)=
**record, field meaning.** The array in which a member's planes are stored. It is a circular buffer: slippage rotates it by changing the index `first` and never moves a plane. Window slice $i$ lives at record index $1 + \mathrm{mod}(i - 1 + \texttt{first}, n_{\mathrm{slice}})$ ([](fel-physics.md#sec-slippage)). Every reader of the field in window order unrotates through that map, and every file is written unrotated.

(gl-rotation-seam)=
**rotation seam.** The place in the record where window order wraps. The slice at record index $n_{\mathrm{slice}}$ is followed, in window order, by the slice at record index 1. A range of slices can straddle it. A reader that takes the record in its stored order there reads two separate runs of slices.

(gl-slippage)=
**slippage, slippage residual.** Light advances one wavelength per undulator period, so the field moves toward the head by `n_wavelength` wavelengths for every slice spacing. The run accumulates the slippage in wavelengths and rotates the record one slice when the total passes 0.8 of `n_wavelength`. The part not yet rotated is the residual, written into every field file as `slippageResidual` ([](fel-physics.md#sec-slippage)).

(gl-bank)=
**escaped field, bank.** The light slippage carries out of the window at its head. The slice that leaves is banked, plane by plane, with the position along the lattice at which it left. The bank is on by default, `run.keep_escaped_field` ([](reading-output.md)).

(gl-pulse)=
**pulse.** The whole radiation pulse at the exit: the window's planes and the banked planes behind them, each banked plane drifted in free space to the exit position. It is the file `<out_root>-pulse.wf.h5`. The window energy, the escaped energy and the pulse energy are three different numbers in the statistics file ([](reading-output.md)).

(gl-push)=
**push.** The particles' part of a step: the transverse maps and the longitudinal equations of motion, integrated with the field held fixed ([](fel-physics.md#sec-step)).

(gl-deposit)=
**deposit, source.** The source is the current each slice's particles drive the field with. The deposit puts it on the slice's transverse grid, each beamlet at one point of it, before the field is advanced ([](fel-physics.md#sec-field)).

(gl-solve)=
**field solve.** The field's part of a step: the deposited source is added and every plane is diffracted over the step by a pair of transforms ([](fel-physics.md#sec-field)).

(gl-kick-force)=
**kick, second meaning.** In the unaveraged method, the force of the radiation field on a particle over a substep. The kick and the deposit are a pair, and the energy the kick takes from the beam is what the deposit gives the field ([](fel-physics.md#sec-unaveraged)). Bmad's `hkick` and `vkick` attributes are a third thing, an element's steering.

(gl-source-filter)=
**source filter.** A factor on the transformed source that falls to zero at wide transverse angles, `run.source_filter`, on by default. A beamlet at one grid point radiates into every angle the grid carries, where a real undulator radiates into its central cone, and the filter removes the part outside ([](fel-physics.md#sec-source-filter)).

(gl-propagator)=
**slice propagator, spectral propagator.** The two ways the field is advanced and slipped, `field.propagator`. The slice propagator, the default, diffracts each slice's plane at the one reference wavenumber and applies the slippage by rotating the record a whole slice at a time. The spectral propagator transforms the window along $\tau$ as well, diffracts each frequency at its own wavenumber and applies the slippage continuously ([](spectral-field.md)).

(gl-kernel)=
**kernel.** The factor a transformed field is multiplied by to diffract it over a step, $\exp(-i\,\delta z\,k_\perp^2/2k)$. The wake pages use the word for a wake function, and the device pages for a routine compiled for the GPU.

(gl-ramp-phase)=
**ramp, first meaning.** Under the spectral propagator, the phase $\exp(-i\kappa\delta)$ that slips the field toward the head by $\delta$, which need not be a whole slice ([](spectral-field.md#spectral-what-changes)).

(gl-pad)=
**pad, sponge.** Under the spectral propagator, 16 planes that follow the window's head along $\tau$: 4 free ones, then 12 that a $\cos^2$ factor, the sponge, brings to zero. Light the ramp carries past the head enters the pad and is banked or absorbed there, so none of it comes round to the tail ([](spectral-field.md#spectral-bank)).

## The lattice and the walk

(gl-line)=
**line, element.** Bmad's lattice and the elements that make it up. Every element is Bmad's own, and the FEL ones are wigglers with an attribute.

(gl-fel-segment)=
**FEL segment, undulator.** A Bmad `wiggler` or `undulator` element whose `fel_method` is `averaged` or `unaveraged`. The FEL parameters derive from its lattice attributes ([](fel-physics.md#sec-element)).

(gl-averaged)=
**averaged, unaveraged.** The two FEL methods, set per element. The averaged method integrates the wiggle-averaged equations on Bmad's own transverse maps. The unaveraged method integrates through the analytic undulator field with no period averaging ([](index.md#intro-two-methods), [](fel-physics.md#sec-unaveraged)).

(gl-interlude)=
**interlude.** An element between FEL segments: a drift, a quadrupole, a phase shifter, a chicane. Each slice's bunch is tracked through it by Bmad's `track1_bunch`, and the field drifts across it in free space ([](fel-physics.md#sec-seam)). `run.interlude_model = "genesis"` takes Genesis4's transcribed maps instead, for the comparison tiers ([](#gl-tier), [](fel-physics.md#sec-interlude)).

(gl-seam)=
**seam.** Where the FEL walk hands each slice to Bmad and takes it back, outside the FEL elements. A slice becomes a Bmad bunch by plain copies, `track1_bunch` tracks it, the field drifts by `wavefront_drift`, and the reference phase advances once per element ([](fel-physics.md#sec-seam)).

(gl-step)=
**step.** The FEL integration step of a segment, the element's `ds_step`, a whole number of undulator periods. Each step pushes the particles, deposits the source and advances the field, and accumulates the slippage ([](fel-physics.md#sec-step)).

(gl-substep)=
**substep.** The unaveraged method's own step, `run.unaveraged_steps_per_period` of them to an undulator period, 20 by default. A record step of the averaged kind holds many substeps ([](fel-physics.md#sec-unaveraged)).

(gl-ramp-end)=
**ramp, second meaning.** The periods at each end of an unaveraged segment over which the undulator field rises from zero and falls back, `run.unaveraged_ramp_periods`, 2 by default. A segment with none is a hard edge ([](fel-physics.md#sec-unaveraged)).

(gl-walk)=
**walk.** One pass down the line, from `run.track_start` to `run.track_end`, element by element and step by step. The footer's timing table partitions the walk into its phases, and the header's estimate is an estimate of the walk ([](#gl-estimate), [](performance.md)).

(gl-comb)=
**comb.** The positions along the line where the run takes a statistics row, every `run.comb_ds_save` metres, with an element end always one of them. The frame series rides the comb ([](user-guide.md)).

(gl-record-stats)=
**record, statistics meaning.** One row of the statistics file: one position on the comb, with every per-slice and whole-window quantity at that position. The record number is the file's axis, and `s` rides along as a variable ([](reading-output.md), [](BMAD-STATS-SPEC.md)).

## The files

(gl-statistics-file)=
**statistics file, stats file.** `<out_root>.stats.h5`, the run's history: one record per comb position, every diagnostic along the record axis and the slice axis, with the lattice table and the run's provenance ([](reading-output.md)).

(gl-dump)=
**dump.** An openPMD beam file or field file written at an element's end, `run.dump_beam_at` and `run.dump_field_at`, or at the exit, `-final`. A dump holds the whole window and is a point a later run can start from ([](reading-output.md)).

(gl-frame)=
**frame, frame series.** A beam file and a field file written at every comb position when `run.dump_at_comb` is on. Frame $n$ and statistics record $n$ are the same position, so the two plot on one axis. A frame inside an element holds the coordinates the method produced there, and a frame at an element's end is that element's dump ([](reading-output.md)).

(gl-slab)=
**slab, slice range.** The slices a frame writes when `run.dump_slice_first` and `run.dump_slice_last` cut it to a range of the window. The field file carries those slices alone, placed where the whole window would put them, and the frame costs those slices and not the window ([](reading-output.md), [](performance.md#perf-frames-slab)). The whole window is the slab when the range is not set.

(gl-checkpoint)=
**checkpoint, continuation.** A dump or element-end frame at an eligible boundary carries the group `lucifer/checkpoint`, the state a continuation needs beyond the standard records. `run.continuation = true` declares that the beam and field files given are such a pair and continues the run from them to the bit ([](reading-output.md)).

(gl-ledger)=
**ledger.** The unaveraged mode's energy accounting, one row per record, `<out_root>.ledger.txt`: the beam's energy change, the field's, the escaped and spontaneous terms, and the closure ([](fel-physics.md#sec-unaveraged)).

(gl-header)=
**header, footer.** What a run prints before and after the walk. The header states the derived numbers with their origins, the work and the estimate. The footer states the timing by phase, the exit energies and every file the run wrote ([](user-guide.md)).

## The device and single precision

(gl-device)=
**device, backend, host.** The device is the GPU a run puts the FEL step on, `run.device`. The backend is the code that drives it, and the one backend is `"metal"`, for an Apple Silicon GPU. The host is the CPU side of the same run. With no device the whole run is on the CPU ([](input-reference.md#param-global-device)).

(gl-fp32)=
**FP32, FP64.** Single and double precision, 32-bit and 64-bit floating point. The CPU path is FP64. The Metal device has FP32 alone, so the FEL step on it is written in forms that hold in FP32 ([](validation.md#val-fp32-lockstep)).

(gl-twin)=
**twin, lockstep, freerun.** The twin is an FP32 copy of the FEL step that `run.fp32_check` runs beside the FP64 step on the CPU, in `src/fel_fp32_mod.f90`. It changes nothing in the run and records how far the two differ, quantity by quantity, in `<out_root>.fp32.txt`. In `"lockstep"` the twin starts every step from the FP64 state, so the difference is that of one step. In `"freerun"` it keeps its own state, so the difference accumulates. The pages also call this check the instrument ([](input-reference.md#param-global-fp32-check)).

## The checks

(gl-validation-suite)=
**validation suite, check section.** The validation suite is `tests/run_fel_validation.sh` with its check scripts. It runs the tiers and then the check sections, each section a script that holds one part of the program to identities, closed forms and refusals ([](validation.md#val-what-is-tested)).

(gl-tier)=
**tier.** One comparison of a Lucifer run with a Genesis4 run from the same starting dumps, named `tier1`, `td1`, `tdsase` and so on. There are eleven. A comparison tier is one that selects a transcribed Genesis4 model to make the comparison ([](validation.md#val-validation-from-one-command)).

(gl-level)=
**level, recorded digit.** A level is a number a check measured: a largest relative difference, a ratio, a closure. The levels are written on the generated pages as recorded digits, and a run of the checks must reproduce them ([](validation.md#val-the-keystone-rule)).

(gl-keystone)=
**keystone.** The run every commit passes before it is merged: the validation suite on the debug build and on the release build, Bmad's regression suite, the wavefront validation and the examples, and then the generated pages regenerated and required to equal the committed ones ([](validation.md#val-what-is-tested)).

(gl-estimate)=
**estimate.** The run time the header predicts before the walk, from the counts of the run's operations and constants fitted to measured runs, in `src/fel_cost_mod.f90` ([](performance.md#perf-estimate)).
