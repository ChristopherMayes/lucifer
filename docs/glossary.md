---
title: Glossary
short_title: Glossary
---

The words this documentation uses for the pieces of a run, each with the page that defines it. Two words carry two meanings each, and both meanings are listed. A `record` is a row of the statistics file in one place and the field's stored array in another. A `plane` is a transverse field sample in one place and a position along the lattice in another.

## The beam and the time window

**wavelength.** The radiation wavelength $\lambda_s$, `field.lambda0`. Every longitudinal length below is a whole number of it.

**slice.** One longitudinal sample of the run. The beam's particles are partitioned into slices, and the field holds one transverse grid per slice. Slices are spaced by $\Delta = $ `window.n_wavelength` $\cdot\lambda_s$, so a slice is `n_wavelength` wavelengths long ([](fel-physics.md#sec-window)). A steady-state run is one slice.

**slice spacing, `n_wavelength`.** The whole number of wavelengths one slice spans. It is an integer end to end, since slippage moves the field by exactly one slice per `n_wavelength` undulator periods ([](fel-physics.md#sec-slippage)).

**window.** The set of slices the run carries, `window.n_slice` of them, fixed in the beam's frame. Slice 1 is the tail and the highest index is the head, which arrives first. Light slips toward the head and leaves the window there. A default window holds every particle and, ahead of the head, the slippage the line applies ([](fel-physics.md#sec-window)). The window is *not* the pulse: see `pulse` below.

**head, tail.** The high-index end and the low-index end of the window. Radiation leaves at the head, and fresh vacuum enters at the tail.

**window order, time order.** Slices in the order of their index, 1 to `n_slice`. Every file is written in this order, and so is every statistics array over slices.

**bunch, beam.** The particles of every slice together. Bmad's `beam_init` generates the bunch or reads it from a file, and the loader slices it ([](fel-physics.md#sec-loading)).

**macroparticle, weight.** A tracked particle standing for `weight` electrons. The weights are kept per particle, so a slice with little charge has light particles rather than few.

**beamlet.** A group of `beamlet_size` macroparticles with one transverse sample and one energy, their phases spaced evenly over $2\pi$ so that every bunching harmonic below the group's size cancels. That is the quiet start, and the shot noise is imposed per beamlet ([](fel-physics.md#sec-quiet), [](fel-physics.md#sec-noise)).

## The field

**grid.** The transverse grid of one slice's field, `grid_n_pts` points a side at `grid_half_width` from the axis, derived from the beam when the input file states neither ([](fel-physics.md#sec-field)).

**plane, first meaning.** One slice's field: one transverse grid of complex values, one for each polarization. The banked planes of the escaped field are planes in this sense.

**plane, second meaning.** A position along the lattice, as in the entry plane and the exit plane of an element, or the plane a frame is named for. Inside an unaveraged segment the stored field sits half a substep from the plane the particles stand on, and a frame is carried back onto its own plane for the write ([](fel-physics.md#sec-unaveraged)).

**polarization, component.** The field's `x` and `y` components. A planar undulator with no tilt radiates one, and a tilted or helical one radiates both ([](fel-physics.md#sec-vector)).

**field set, harmonic, member.** The fields a run carries, one member per harmonic, the fundamental first. Each member is a wavefront with its own slices and polarizations, and the statistics file has a group per member ([](fel-physics.md#sec-field-set)).

**wavefront.** The field structure: the planes of one member over the slices, with the grid spacing, the slice spacing and the wavelength. It is what the openPMD field files hold ([](fel-physics.md#sec-openpmd)).

**record, field meaning.** The array in which a member's planes are stored. It is a circular buffer: slippage rotates it by changing the index `first` and never moves a plane. Window slice $i$ lives at record index $1 + \mathrm{mod}(i - 1 + \texttt{first}, n_{\mathrm{slice}})$ ([](fel-physics.md#sec-slippage)). Every reader of the field in window order unrotates through that map, and every file is written unrotated.

**rotation seam.** The place in the record where window order wraps. The slice at record index $n_{\mathrm{slice}}$ is followed, in window order, by the slice at record index 1. A range of slices can straddle it. A reader that takes the record in its stored order there reads two separate runs of slices.

**slippage, slippage residual.** Light advances one wavelength per undulator period, so the field moves toward the head by `n_wavelength` wavelengths for every slice spacing. The run accumulates the slippage in wavelengths and rotates the record one slice when the total passes 0.8 of `n_wavelength`. The part not yet rotated is the residual, written into every field file as `slippageResidual` ([](fel-physics.md#sec-slippage)).

**escaped field, bank.** The light slippage carries out of the window at its head. The slice that leaves is banked, plane by plane, with the position along the lattice at which it left. The bank is on by default, `run.keep_escaped_field` ([](reading-output.md)).

**pulse.** The whole radiation pulse at the exit: the window's planes and the banked planes behind them, each banked plane drifted in free space to the exit position. It is the file `<out_root>-pulse.wf.h5`. The window energy, the escaped energy and the pulse energy are three different numbers in the statistics file ([](reading-output.md)).

## The lattice and the walk

**line, element.** Bmad's lattice and the elements that make it up. Every element is Bmad's own, and the FEL ones are wigglers with an attribute.

**FEL segment, undulator.** A Bmad `wiggler` or `undulator` element whose `fel_method` is `averaged` or `unaveraged`. The FEL parameters derive from its lattice attributes ([](fel-physics.md#sec-element)).

**averaged, unaveraged.** The two FEL methods, set per element. The averaged method integrates the wiggle-averaged equations on Bmad's own transverse maps. The unaveraged method integrates through the analytic undulator field with no period averaging ([](index.md#intro-two-methods), [](fel-physics.md#sec-unaveraged)).

**interlude.** An element between FEL segments: a drift, a quadrupole, a phase shifter, a chicane. Each slice's bunch is tracked through it by Bmad's `track1_bunch`, and the field drifts across it in free space ([](fel-physics.md#sec-seam)). `run.interlude_model = "genesis"` takes Genesis4's transcribed maps instead, for the comparison tiers ([](fel-physics.md#sec-interlude)).

**seam.** Where the FEL walk hands each slice to Bmad and takes it back, outside the FEL elements. A slice becomes a Bmad bunch by plain copies, `track1_bunch` tracks it, the field drifts by `wavefront_drift`, and the reference phase advances once per element ([](fel-physics.md#sec-seam)).

**step.** The FEL integration step of a segment, the element's `ds_step`, a whole number of undulator periods. Each step pushes the particles, deposits the source and advances the field, and accumulates the slippage ([](fel-physics.md#sec-step)).

**substep.** The unaveraged method's own step, `run.unaveraged_steps_per_period` of them to an undulator period, 20 by default. A record step of the averaged kind holds many substeps ([](fel-physics.md#sec-unaveraged)).

**walk.** One pass down the line, from `run.track_start` to `run.track_end`, element by element and step by step. The footer's timing table partitions the walk into its phases, and the header's estimate is an estimate of the walk ([](performance.md)).

**comb.** The positions along the line where the run takes a statistics row, every `run.comb_ds_save` metres, with an element end always one of them. The frame series rides the comb ([](user-guide.md)).

**record, statistics meaning.** One row of the statistics file: one position on the comb, with every per-slice and whole-window quantity at that position. The record number is the file's axis, and `s` rides along as a variable ([](reading-output.md), [](BMAD-STATS-SPEC.md)).

## The files

**statistics file.** `<out_root>.stats.h5`, the run's history: one record per comb position, every diagnostic along the record axis and the slice axis, with the lattice table and the run's provenance ([](reading-output.md)).

**dump.** An openPMD beam file or field file written at an element's end, `run.dump_beam_at` and `run.dump_field_at`, or at the exit, `-final`. A dump holds the whole window and is a point a later run can start from ([](reading-output.md)).

**frame, frame series.** A beam file and a field file written at every comb position when `run.dump_at_comb` is on. Frame $n$ and statistics record $n$ are the same position, so the two plot on one axis. A frame inside an element holds the coordinates the method produced there, and a frame at an element's end is that element's dump ([](reading-output.md)).

**slab, slice range.** The slices a frame writes when `run.dump_slice_first` and `run.dump_slice_last` cut it to a range of the window. The field file carries those slices alone, placed where the whole window would put them, and the frame costs those slices and not the window ([](reading-output.md), [](performance.md#perf-frames-slab)). The whole window is the slab when the range is not set.

**checkpoint, continuation.** A dump or element-end frame at an eligible boundary carries the group `lucifer/checkpoint`, the state a continuation needs beyond the standard records. `run.continuation = true` declares that the beam and field files given are such a pair and continues the run from them to the bit ([](reading-output.md)).

**ledger.** The unaveraged mode's energy accounting, one row per record, `<out_root>.ledger.txt`: the beam's energy change, the field's, the escaped and spontaneous terms, and the closure ([](fel-physics.md#sec-unaveraged)).

**header, footer.** What a run prints before and after the walk. The header states the derived numbers with their origins, the work and the estimate. The footer states the timing by phase, the exit energies and every file the run wrote ([](user-guide.md)).
