---
title: Lucifer changelog
short_title: Changelog
---

# Lucifer changelog

Development history of the FEL tracker on the `lucifer-dev` branch, newest first.
This is the branch's own record. Bmad's `changelog.md` carries what a merge changes,
and it is written at the merge.

- 2026-09-11 Added: the run header projects what the device's FP32 transform pair will cost a deck. The pair
  loses about 0.17 of an FP32 quantum of the field's energy a butterfly stage, in one direction, and how
  much a deck loses is set by how often the pair runs, once a record step in an averaged element and nsub
  times in an unaveraged one. The header counts the pairs a field record sees over the line and states the
  projected loss, 6.5e-4 on the one-segment unaveraged deck and 1.1e-5 on the averaged one. The check holds
  the count to the deck's and the projection within a factor of two of what the same run measured.
- 2026-09-11 Added: the frame series with the device resident in the unaveraged mode is checked. That
  combination became supported when the Metal backend took that mode and nothing in the harness held it,
  which is the shape a producer and a consumer written to the same feature take when neither is exercised.
  A device run writes the same nine beam and nine field frames the CPU run writes at the same s exactly, a
  frame taken inside the segment names its chart through felMethod, and its mean px carries the undulator
  quiver at 4.2e5 eV/c where the averaged chart sits near 1e2. The comparison is the mean and not the
  spread, an rms being blind to the common offset the quiver is, and the device's mean tracks the CPU's to
  1.9e-7.

- 2026-09-10 Added: global%interlude_ds_step cuts an interlude element into pieces, so statistics rows and
  openPMD frames land along a drift or a quadrupole the way they already land along a wiggler. A frame
  series had particles through an undulator and one frame through the break after it, and no comb setting
  changed that. comb_ds_save is the selecting half of Tao's comb, save_a_bunch_step's, which chooses among
  the positions a tracker already reaches. Tao's beam-track loop carries a second half that cuts an element
  the tracker would otherwise cross in one map, and this is it. An FEL element reaches one position per
  ds_step, so the comb was already fine-grained inside a wiggler. An interlude was one track1_bunch call
  over the whole element, so it reached one position, at its end.

  Each piece tracks the beam through its own sliced element, drifts the field by its own length, advances
  the common phase over it and takes its share of the element's slippage, which is an accumulator in whole
  slices and therefore sums to the element's own. The pieces come from Bmad's element_slice_iterator, the
  routine Tao's comb uses, so an element's fringes stay at its ends rather than appearing at every cut. On a
  wiggler, quadrupole, pipe, wiggler line a cut run's exit power agrees with the whole-element run's at
  1.2e-13 and its exit sigma_x at 6.0e-16. Zero, the default, tracks the element whole: a deck naming zero
  reproduces a deck that never names it byte for byte in the diag stream and on every array of the
  statistics file but the input echo and the timestamp, which is what leaves the recorded tier digits alone.

  Two configurations decline the cut, both in the shape Tao declines where CSR or space charge would see
  one. An element carrying a Bmad short-range wake is tracked whole, that wake being a once-per-passage kick
  of the element's length. interlude_model = "genesis" refuses the setting, that model being a transcription
  whose step is Genesis's own element. The chamber wake needs no refusal: it is an energy loss proportional
  to the length applied additively in gamma, so the pieces sum to the element's kick.

  The price is the field's drift rather than the beam's tracking, since a piece drifts every slice of the
  window where it tracks one bunch per slice. On 32 slices of 512 macroparticles at ngrid 128 the walk goes
  from 0.285 s whole to 0.325 s at ten pieces an element and 0.427 s at thirty, and the drift goes from 1.3
  percent of the walk to 24 percent while the FEL step is unmoved. doc/performance.md carries the table.

- 2026-09-10 Changed: doc/input-reference.md said comb_ds_save carries Bmad's bunch_track_struct%ds_save
  name and semantics. It carries the selecting half of them. The entry now says which half and points at
  interlude_ds_step for the other, and says what the setting does and does not do inside an interlude.

- 2026-09-10 Added: the Metal backend runs the unaveraged mode. device = "metal" refused that mode until now,
  so the only way to resolve the undulator quiver was the CPU path, whose parallelism is over slices and gives
  a steady-state deck one core of twelve. Three kernels are the mode's own, the quiver-resolving Strang push,
  the radiation kick with its deposit from the resolved motion, and the ledger's spontaneous reduction. The
  rest of a substep is the kernels the averaged path already had: the four-pass solve, the accumulator's clear
  and the fixed-point scatter. On one Aramis segment the device runs 9.0x a twelve-thread CPU at 32 slices of
  512 particles and 10.6x at 4096, and 10.1x and 27.7x on a single slice where the CPU has no parallelism to
  use. The backend is busy for three quarters to nine tenths of the walk, so the mode is device-bound.

  The lag rides the 64-bit tick accumulator rather than a single-precision residual. dtau/ds depends on gamma,
  ux and uy alone, so the push integrates it from zero over a substep and adds the increment exactly, and no
  arithmetic differences the increment against a lag that grew. The phase row reads 3.4e-6 rad against the CPU
  twin's 9.9e-6 and the phasor row 6.5e-8 against 7.7e-8. The three per-substep quantities that do not depend
  on the particle are built in FP64 on the host and uploaded once a substep: the envelope and the undulator
  phase at the four RK stage positions, the optical carrier's base rotator, and the propagator at the substep.

  The instrument judges it. The device takes the twin's role in the unaveraged lockstep, the nine rows land
  inside recorded ceilings on a steady deck and an 8-slice window, the FP64 diag and ledger are byte identical
  with the twin on against off, and the mutation takes the phase row from 3.4e-6 to 9.8e-5 and the field row
  from 6.4e-6 to 7.7e-4. Two device runs of one unaveraged deck agree on their diag and ledger bit for bit at
  4 threads against 8: the deposit accumulates in fixed point and the ledger's spontaneous term reduces with
  one threadgroup owning one slot and no atomic at all. Two live polarizations inside the mode is refused,
  which the CPU carries and the device does not.

  One measurement sets those ceilings and it predates the mode. An FP32 forward transform, propagator multiply
  and inverse transform lose about 1.3e-7 of the field's energy every time they run, systematically, where the
  FP64 pair drifts 1.6e-12 over 5340 applications. The averaged mode runs the pair once a record step and a
  12-segment line loses 1.4e-4 of the field energy, under every ceiling it records. This mode runs it nsub
  times, so an 89-step segment of 60 substeps loses 7.2e-4, which is the whole of a weakly seeded deck's
  1.7e-3 exit-power divergence from the CPU. The unaveraged energy ledger therefore does not close on the
  device at the 4.2e-6 of turnover it closes at on the CPU, and doc/validation.md says so with the numbers.

- 2026-09-10 Fixed: the device deposit's residual-coarsening mutation hook returned its argument unchanged.
  The step was written as fabs(d) * 3.05175781e-5, a fixed fraction of the value, and that constant is exactly
  2^-15, so d/step is exactly 32768 for every d and rounding it to the nearest integer is the identity. The
  Fortran it transcribes uses 256 * spacing(d), whose step comes from the exponent alone and really quantizes.
  check_device's mutation assertions passed throughout on the other hook, the coarsening of the uploaded phase
  accumulator, which covers the lockstep path but not freerun's deposit. The shader now takes the step from
  frexp and ldexp. The mutated levels move further than they did and the straight levels are unchanged.

- 2026-09-10 Added: the unaveraged twin carries a ledger row and a field row, and the stream's source and
  field columns hold them. The ledger row is the energy the twin's own kicks took from the beam against the
  FP64 step's, scaled by the energy the step moved rather than by the energy it netted. Over a record step
  the gains and the losses very nearly cancel, so the net is a small difference of large exchanges, and
  against it the ratio reaches 3, which says something about the cancellation and nothing about precision.
  Scaled by the turnover the row reads 1.56e-5. The field row is the twin's own single-precision record,
  advanced on every substep by a single-precision transform pair against a rounded propagator with the
  step's source added, compared with the FP64 field in the L2 norm. Sixty substeps stand behind one row on
  the benchmark segment and it reads 6.0e-6. check_unaveraged holds the two under 1e-4 and 5e-5, and
  doc/validation.md carries the levels beside the six that were already there.

  The field row runs the single-precision transform inside the slice loop, so fel_fp32_mod's FFTW plan cache
  moved from routine-local save variables to module scope under threadprivate, with a warm-up that fills
  every thread's cache before the loops start. That is wavefront_mod's arrangement and its reason: the axis
  of parallelism is the slice and not the transform, and FFTW's new-array execute rule makes a plan on a
  differently aligned array undefined, so each thread owns its plans and its aligned buffer. The build holds
  a named lock, the planner and the allocator carrying global state where the executor carries none. Four
  threads planning at once segfault inside the library, and a single-threaded run never shows it.

  The instrument is 1.75 times the walk of the same run without it on that segment, 6.1 s against 3.5 s at
  one slice, which doc/performance.md records. Nothing in a production run turns it on.

- 2026-09-10 Added: the unaveraged advance has a single-precision twin, so fp32_check no longer refuses that
  mode. Apple GPUs carry no FP64 and this tree's rule is that a single-precision form of an advance exists
  with its divergence from FP64 measured before a kernel is written for it. The averaged advance has had one
  since fel_fp32_mod; the unaveraged advance had none, so nothing could have judged a device port of it.

  The twin advances every substep of a record step from the state that step received and is compared where
  the record step ends, a substep being internal to the integrator where a record step is where the FP64 path
  and a port would synchronize. Four quantities single precision destroys there, each measured on a real
  mid-segment state: the slippage rate gamma/u_s - 1/beta0, whose two terms are both one to within 1.3e-8 and
  whose difference is 0.016 of the quantum of one, so the naive form returns exactly zero and the slippage
  disappears; gamma, whose per-substep change is a thousandth of its quantum; the lag, resolved at 4.4e5
  quanta on one slice and 0.84 on a 351-slice window; and the phase carried whole, at 2.7e-5 rad by a
  segment's end. The replacements are the offset charts the averaged twin already uses, a base phase held in
  FP64 one number a slice a substep, and for the slippage the exact identity that moves the cancellation into
  a difference of two small like quantities. u_s needs no reformulation, which measurement settled against
  expectation: the terms it loses are 6.7e-9 of the result.

  Two are mutation records rather than arguments. The naive slippage takes the phase row from 9.9e-6 rad to
  1.7 rad, which is comparable to pi and leaves the bunching phase meaningless, and an absolute gamma takes
  the energy row from 2.1e-6 to 2.0e-4. check_unaveraged asserts five rows, the guard on the lag residual,
  and that the FP64 diag and ledger are byte identical with the twin on against off. The field is not
  reached: this mode diffracts sixty times a record step where the averaged one diffracts once, so its
  single-precision field accumulation is a separate quantity and the stream's source and field columns are
  zero in the mode to say so. freerun stays refused, that mode compounding a state this twin does not carry.

- 2026-09-09 Fixed: an unaveraged run can write a frame series. global%dump_at_comb died at the first comb
  position inside a segment, so the mode could produce no particle frames at all and dump_at_comb = F was
  the only way to finish, which writes the final beam and nothing else. The chart assertion sat at the top
  of fel_slice_to_bunch, which converts a slice into a Bmad bunch, and quiver_in_px is true for the whole
  interior of an unaveraged segment. That refused every caller, the frame writer among them.

  The conversion is the same six numbers whichever chart the beam is in, so the assertion belongs to what
  the caller does next rather than to the conversion. A caller handing the bunch to Bmad's tracking needs
  the averaged convention and gets the refusal by default. The frame writer states that its file names the
  chart, since fel_frame_attributes writes felMethod beside the element and the undulator, and is allowed
  the chart the beam is in. The beam dump of fel_dump_beam writes no felMethod and keeps the refusal, which
  costs nothing: dump_beam_at resolves through Bmad's locator and lands on element boundaries, where the
  ramps have handed back the averaged convention.

  Measured on the one-segment unaveraged example at a 0.5 m comb: nine frames where there were none, each
  carrying felMethod Unaveraged, and the quiver is in them. The mean px swings over 2.9e5, -5.4e4, -3.6e5
  and -4.2e5 eV/c across successive frames, which is the aw/gamma amplitude of 4.3e5 the chart implies,
  where the averaged run holds -6.5e2 throughout. The rms is unmoved at 1.46e4 either way, the quiver being
  a common offset at a position along the undulator rather than a spread.

- 2026-09-09 Fixed: the device deposit's overflow bound is now checked where it was assumed, and its clear
  costs a quarter of what it did. The bound took the reference gamma over a fixed floor, where a
  contribution carries w/gamma and so needs the lowest gamma the run will ever deposit at: setup measures
  that over the loaded beam, every readback checks the beam has not fallen under an eighth of it, and a
  particle that has refuses the run. The roll-off took the box half width with no axis offset and no
  rotation, where the kernel offsets and rotates before forming it, which was short by up to four and a
  half in the roll-off term. And the scale is now required to be finite and normal in the kernel's single
  precision: the headroom refusal beside it cannot see that, the scale being chosen against the same bound
  the headroom is measured against, so a bound below about 5e-20 V/m or above 4e56 would have deposited
  through an infinity while the headroom read a healthy 37 bits. No deck reaches either edge, the tested
  ones sitting thirty orders inside the nearer.

  What the quantum costs is now measured. global%device_dep_mutate shifts the scale by whole bits and
  moves nothing else, so the source row's composition can be read off: on the cancelling-phase deck it is
  1.4673e-05 at the derived scale and the same to every digit with the quantum sixteen million times
  coarser, first moves at 26 bits coarser and reaches 1.05x at 29, one short of the refusal. The level the
  precision cases record is therefore phase and FP32 arithmetic, not quantization.

  The clear writes four words a thread where it wrote one. It is a memory fill, it was the largest single
  pass on a many-slice window at 0.0863 s of 0.328 s busy, and widening the store halves it. The same
  widening on the float deposit is worth nothing, 0.99x and 1.01x, since that clear writes half the words
  and was never the bottleneck, so the saving belongs to the fixed point's own bytes. On a 96-slice window
  at ngrid 256 the deposit's price falls from 1.17x to 1.09x of the float walk at 1024 particles a slice
  and from 1.31x to 1.11x at 8192, measured six pairs alternating between the two builds. Output is
  unchanged: two runs still agree on all 127 arrays and every recorded precision level holds to its
  printed digits.

  The gamma the bound assumes is checked in the kernel rather than at a readback. A readback sees an
  element's last step, so a particle can cross the floor and return before anything looks, and a breach
  found there is found after that element's deposits have used the bound. The kernel tests every particle
  before its contribution is converted, drops one that fails, and records a fault the host reads at the
  next drain and stops the run on, ahead of any stats row or dump drawn from the step. It costs nothing
  the walk can resolve. Guarding there also covers the twin, which is never resident and whose readback a
  host check never reaches.

- 2026-09-09 Fixed: the field solver's transverse propagator carried the wrong wavenumbers on an even grid.
  The table was built as an offset index, which lands on whole multiples of the fundamental only when the
  point count is odd. On an even count every mode was propagated half a step out and the constant mode at
  half a step rather than at zero, an error that grows with the mode. The kernel now builds the wavenumbers
  the transform carries. On an odd grid the two tables are the same one, so the change is bit-identical
  there, verified by building both and running one deck in one directory, and all eleven tier digits are
  unmoved. On an even grid it is a correction of 1.4e-1 peak normalized on the 57 m line at 256 points.
  Genesis4 builds the same table the same way, so the transcription was faithful and the defect is upstream,
  and the two codes agreed with each other at 1.4e-6 on an even grid while both were wrong. What found it
  was the tier comparing the Bmad seam against the transcribed interlude, the seam moving the field with
  wavefront_drift, whose wavenumbers were right at both parities.

- 2026-09-09 Changed: every deck in the tests and examples runs a power-of-two transverse grid, and the
  derived count rounds up to one instead of only doing so for the device. A deck now runs on the CPU and on
  the Metal field solver unchanged: the LCLS example at 300 slices of 4096 macroparticles and grid 128 lands
  at 7.407 TW against 7.404 TW on exit power, a relative 4.1e-4 inside the device's 5.1e-4 band, where before
  it derived 127 points on the CPU and 128 on the device and the two could not be compared. Six replacements
  in check_device that existed only to patch an odd deck to an even one are gone with it. The cost is CPU
  time and it is priced in doc/performance.md: 256 points take 1.17x the walk of 255 on the Aramis line,
  FFTW taking 3 x 5 x 17 better than it takes 2^8. Four places hold an odd grid on purpose and say why,
  the five Genesis reference decks, the harmonic and phasing checks, the saturation demo, and the smoke
  test with its namelist twin, which keeps the odd-grid path under test.

- 2026-09-09 Fixed: two device runs of one deck now agree bit for bit. The source deposit accumulated with
  atomic float adds, four corners a macroparticle and two a corner, and float addition is neither associative
  nor commutative, so the order threads reached a cell set the bits. On the time-dependent window with
  shot noise, 54 of the 127 arrays of the statistics file differed between two runs, the field power by 2.0e-7
  and the exit power by 1.6e-8. That sat far inside the 5.1e-4 band against the CPU and was never an accuracy
  problem, but it cost every device output the right to be compared exactly, and check_device measured the
  migration no-op against the device's own spread between two identical runs for want of one.

  The deposit accumulates in fixed point instead. Integers do not depend on arrival order, so the answer does
  not either. Sixty-four bits are carried as two 32-bit words because this hardware has no 64-bit atomic of
  any kind, atomic_ulong existing as a type with store, load, fetch_add, fetch_min, fetch_or and
  compare_exchange all invalid for it: every low addend is unsigned, so the number of carries is a property of
  the sum rather than of the order. Thirty-two bits were tried and rejected on measurement, a scale that keeps
  the worst case inside an int leaving a quantum that loses to the float accumulation by a factor of seven
  where sixty-four bits beat it by the same factor. The scale is a power of two, so it and its reciprocal are
  exact and the conversion costs one rounding where the float deposit paid one per contribution. It is bounded
  by every macroparticle in one cell in phase at a gamma floored an eighth of the reference. The run's charge
  gives the bound at setup, the element's source scale and the grid's roll-off complete it at the first step,
  and the run prints the scale with its origin and its headroom there: 2^24 ticks per volt per metre against
  a bound of 2.7e11 V/m on the Aramis segment, the quantum 37 bits below the FP32 spacing of that bound. A
  headroom too small to beat the float accumulation refuses the run before anything converts.

  The conversion adds no dispatch: it rides the first of the filter's four source passes with the filter on,
  and the solve's last pass when it is off, since nothing transforms the source then. The instrument's source
  readback decodes the accumulator for the same reason. Accuracy is unchanged, measured against the CPU's FP64
  deposit on cancelling phases, unequal weights, charge in few cells and the largest load: 1.46e-5, 1.46e-5,
  6.95e-6 and 2.51e-6 against the float deposit's 1.46e-5, 1.47e-5, 6.89e-6 and 2.51e-6. The cost depends on how much of a
  walk the backend holds. On a single slice, where the host's dispatch dominates, it is a quarter to a third
  more device time and six to seven percent more wall clock. On a 96-slice window at ngrid 256, where the
  backend is busy for about 68 percent of the walk, it is a quarter to two fifths more device time and 17 to
  31 percent more wall clock. Both are in doc/performance.md, with the per-pass table repeated after the
  change: the clear and the deposit carry all of it and every other pass is unmoved.
  check_device gained the two-run comparison over every array, on both filter settings, and its no-op check
  is equality with its floor deleted.

- 2026-09-09 Added: global%device_timing times each pass of a device step separately, so the deposit's share
  step is a measured number rather than an inference. The backend's busy seconds come from one command buffer,
  and a command buffer holds every step between two host touches, so they price a step, not a pass. On, each
  pass is encoded into its own compute encoder carrying a timestamp counter attachment, still inside the one
  command buffer the step batching rests on, and the run reports seconds and encoder counts for the transverse
  map, the push, the source clear, the deposit, the filter and the solve. Stage-boundary counter sampling is
  what an Apple GPU carries and it times an encoder rather than a dispatch, which is why a pass gets its own;
  the M3 Max reports dispatch-boundary sampling absent, and a device that samples nothing at an encoder
  boundary refuses rather than reporting zeros.

  Measured on the full Aramis line, 12 segments and 1068 steps, one slice, ngrid 256, the filter on, best
  complete run of five: the deposit is 7.7 percent of the six passes at 1024 and at 8192 particles a slice and
  26.8 percent at 131072, and the filter and the solve together are two thirds of a step at production loads.
  The instrument costs 2.61, 2.08 and 1.16 times the device busy seconds at those three loads, so its absolute
  seconds are its own rather than a run's and doc/performance.md says so beside the table. It is off
  by default for that reason. The sample buffer holds 2048 passes, about 157 steps, and a command buffer
  needing more is committed early with a warning naming how often, since that batching was the instrument's.

- 2026-09-08 Added: global%source_filter_tolerance places the source filter's edge relative to the angle it
  protects rather than on that angle. The derived edge was the larger of the mode and rho angles, which puts
  the sigmoid's half-amplitude point on it and so keeps a quarter of the source intensity at the angle setup
  had just derived as physical. The tolerance is the largest fraction of source intensity the filter may take
  anywhere inside the protected angle, and the edge moves out by 1/(1 - width*log(t/(1-t))) with
  t = sqrt(1 - tolerance) to hold it. doc/startup-noise.md records what the old placement costs the in-cone
  startup power, 14 percent where the derived edge sits well inside the cone and 28 to 30 percent where the
  two nearly coincide, and the on-axis width guard cannot see any of it because it constrains one angle.

  The default of 0.75 is the loss the old placement implies rather than a number chosen for it, and its margin
  is exactly one, so the derived edge is bitwise what it was and no recorded digit moves. At a width of 0.05 a
  tolerance of 0.2 moves the edge out by 1.12 and 0.01 by 1.36. The statistics split stays on the protected
  angle, so power_inside_angle answers one question across a tolerance scan. A stated source_filter_angle or
  the grid-relative cuts place the edge directly and bypass the tolerance, and setup says which happened.

  Swept on the four machines of doc/startup-noise.md at 0.75, 0.2, 0.1 and 0.01, one seed a row, the accepted
  power at the first segment rises by 9 to 13 percent at 0.2 and by 12 to 22 percent at 0.01, while the
  wide-angle power rises by 4.9 to 6.2 times and then by 19 to 28 times. The saturation point does not move
  on any machine at any tolerance, and the fitted gain length and the pulse energy stay inside 4 percent. The
  sweep is experiment h of tests/scripts/startup_noise.py. The default is unchanged, since one seed a row
  cannot select one.

  A tolerance outside the open interval from 0 to 1 is refused before anything takes a logarithm, since at one
  the margin runs to positive infinity and the edge collapses onto the axis, which a test on the margin's sign
  would take. A width too large for the tolerance is refused with the width that would meet it. The logarithm
  is taken as log t + log(1+t) - log(tolerance), which never forms 1 - t: that difference is exactly zero for
  any tolerance below about 2.2e-16, and at a width of 0.01 such a tolerance still has a margin of 0.60 and is
  an edge the run can build. check_source_filter gained thirteen checks, and three of them fail against a run
  that reads the input and derives the old edge anyway.

- 2026-09-07 Fixed: a beam dump places its slices in the bunch, so a reader that knows openPMD and
  nothing about this program gets the whole thing. The per-particle time record was the lag inside a
  particle's own slice and timeOffset was a constant reference, so the slice's own position lived only in
  particlePatches. Measured with ParticleGroup on a four-slice frame: a bunch 1.00 slices long, every
  slice on top of the others. The slice's offset now rides p%t and so timeOffset, which openPMD sums with
  time, and the same frame reads 4.00 slices long. vec(5) is untouched, so the time record is byte for
  byte what it was and the windowed composition still restarts mid-line at 3.7e-14. The two numbers stay
  apart on purpose: the lag is a part in 1e8 of the arrival time and a single global time would lose it,
  entirely so in single precision. The cost is eight bytes a macroparticle, 60.2 to 68.2 measured, and
  doc/performance.md says so. convert_genesis.py carries the placement too, so a converted file and a
  native one mean the same thing.

  The harness restarted a run mid-line from its own dumps at 3.9e-14 throughout and could not have caught
  this, since this program wrote and read both ends and re-supplied the offset from the patches. The file
  comparison beside it collected only shaped datasets, so a constant record was invisible and timeOffset
  was the one record it never compared. check_beam_format now expands constant records and gains a
  section that reads a dump with openPMD-beamphysics, which shares none of this program's conventions:
  the slices sit one spacing apart, the bunch spans its window, the reader sees every particle and the
  whole charge, and a frame cut to a slice range places those slices exactly where the whole window does.
  The round trip now also names the three records a restart may move, the element, the path length and
  the arrival time, and holds the placement across it. The placement lives inside an absolute time, so a
  double carries it to 6.5e-6 of a slice spacing and particlePatches stays the exact partition
  (FINDINGS 7.64).

- 2026-09-07 Added: the frame series can be cut to a slice range, and the field's reductions are
  written per record. global%dump_slice_first and dump_slice_last name the slices a frame carries, in
  window order, with -1 as the last meaning the last slice. The beam file's patch count is then the
  range and the field file carries that many slices, both stating sliceFirst and sliceLast. A run at one
  wavelength per slice has thousands of slices and a whole-window field frame at ngrid 255 is 95 MB, so
  a view that wants the core of the bunch asks for it rather than reading it out afterwards. The dumps a
  run restarts from stay whole, and a range outside the window is refused.

  global%dump_reduced writes the field's intensity projections into the statistics file at every record,
  one set per member and plane: the transverse intensity summed over the window, the intensity against
  slice and x summed over y, the same against slice and y, and the complex field on axis. They are
  BMAD-STATS-EXT-FEL's F16, with coords/grid_x and coords/grid_y as their axes, and they integrate back
  to the record's power to 5e-15. Two hundred times smaller than the frames they replace, so a picture of
  a run needs no field frame.

  Every frame now states frameFormat, sliceFirst and sliceLast, the floor position and angles at its own
  s rather than the element's end, and the iteration time as Bmad's reference time there, since openPMD
  orders a series by time and has no notion of s. The Frames line prices the beam, the field and the
  reductions with the range's slice count, at 60 bytes a macroparticle measured rather than assumed.

  Three identities are checked, and they pin the phase convention, the time-order rotation and the slice
  indexing end to end. Per-slice bunching from a frame's own particles agrees with the statistics row to
  4.2e-15, per-slice power from a raw field frame to 7.4e-16, and the projections integrate to the row's
  power to 5.5e-15. doc/performance.md gains the frame measurements. A field frame compresses 1.30 times under
  gzip with shuffle for sixteen times the read, and stores as float32 at half the size for 5.9e-8 of
  relative error, neither of which is taken here. A beam frame carries no overhead to trim. At a comb of
  a tenth of a metre with frames on, writing them is 69 percent of the walk and the device readback is a
  sixth, so a series that costs too much is answered by writing fewer slices rather than by moving
  reductions onto the device. examples/frames runs the whole of it in 0.4 s.

- 2026-09-07 Added: global%dump_at_comb writes the beam and the field at every comb position as an
  openPMD fileBased series, <out_root>-<record>.beam.h5 and <out_root>-<record>.wf.h5 with the record
  zero padded to six digits. The index is the stats record's, so a frame and the row describing it are
  read on one axis, and comb_ds_save sets the frame spacing. The writers are the element-end dumps' own,
  so a frame taken at an element end is byte for byte the file that element's dump_beam_at entry writes.
  The header's Frames line states the count and the bytes before the run tracks, since a comb of zero on
  a long line is tens of gigabytes. Writing frames does not change the run: the field's records are
  rotated to time order to be written and rotated back, so a run with frames is dataset-identical to the
  same run without them, which the harness measures.

  Every macroparticle now carries a label, assigned at load, unique over the window, and unchanged by
  migration, so a series reads as trajectories rather than as unrelated snapshots. It is written as
  openPMD's own id record through coord_struct%ix_user, which Bmad's beam writer and reader gained, and
  a beam read from a dump keeps the labels the file carried. On a migrating deck 419 of 2048 labels
  changed slice between the first frame and the last and none appeared that the first frame did not
  carry.

  Both files say where they were taken: sPosition, elementName, elementIndex and phi0, and inside an FEL
  element felMethod with aw, ku, tilt and helical. The averaged mode integrates the quiver away, so a
  reader that wants the physical orbit rebuilds it from aw, ku and s, which doc/reading-output.md states
  along with which mode it applies to. The diagnostics section gains nine checks over the series.

- 2026-09-07 Changed: doc/performance.md is re-measured with the source filter at its default, and the
  cost estimate's rates are re-fitted to it. One phase moved. The field solve is up 40 percent at 96
  slices and 44 percent at 504, and the walk with it, from 23.917 to 30.333 s and from 127.370 to
  164.650 s. The filter is the whole of it: measured on one deck both ways, it costs 42 percent of the
  field solve, 27 percent of the walk, and nothing of the particle push, which is the arithmetic of a
  third transform per field per step against two. The filter-off walk of 23.392 s reproduces the 23.917 s
  recorded before the default changed, so nothing else moved between the two measurements. Thread scaling
  improved, 9.59x at twelve threads against 9.16x and 80 percent efficiency against 76, since the filter
  lengthens the parallel region and leaves the serial part alone. The sampling split puts the transform at
  68.8 percent of the run against 61.3, and the thread wait at 5.1 against 10.0. The estimate's five
  constants are re-fitted by tests/scripts/fit_cost_rates.py, which reads the logs the runs wrote, so the
  next re-measurement is a command: the transform rate rises 58 percent, the push and the deposit stay
  where they were, the serial fraction falls from 2.8 to 2.1 percent and the walk share rises from 0.83
  to 0.878, with a worst residual of 4.0 percent over thirteen runs. Against the clock the four examples
  now read 0.82 to 0.97 where all four read 0.80, and the eleven tiers read 1.13 to 1.87 because every
  one of them turns the filter off to compare against a code that carries none, so their solve is the
  cheaper two-transform one priced at the filtered default's rate. Every row is inside the factor of two
  the estimate is held to. The phases mode of run_perf_benchmark.sh states slicing%n_wavelength and checks
  the slice count the run built: an unstated spacing is derived from the beam and the gain, which on this
  line is 38 wavelengths rather than 3, so the profile had been measuring 8 slices where it reports 96
  (FINDINGS 7.63).

- 2026-09-07 Added: a run states its cost and its convergence standing before it tracks. The header gains
  three lines. Work counts the integration steps over the line, the macroparticles in the window, their
  product in particle-steps and the transverse grid points the field solve touches per step, all of them
  known before the first step. Estimate turns that into wall clock from rates fitted to the eight
  configurations of doc/performance.md, one for the particle push and two for the field solve, which the
  fit reproduces to within 2 percent. The line names the machine the rates are for, and the footer prints
  the measured walk beside the estimate on every run. Load states the macroparticles in the thinnest
  slice, the beamlet size and the beamlet count, then where that stands against the 1024 macroparticles
  in 128 beamlets the source filter's convergence was measured at, and what the beamlet size buys, which
  is the harmonics the shot noise is imposed at. A harmonic field above (beamlet_size - 1)/2 is refused
  with the beamlet size that would carry it, since it would otherwise start from a load with no noise at
  its own frequency. The convergence footer names the thin slices by index and position, those carrying
  charge whose beamlet count is under 128 or whose power outside the split angle is over a tenth of the
  power inside it, up to eight with a count of the rest, and names sample mode as the remedy: on the
  import example that is 12 slices of 34. global%load_only prints the header and stops, which is the
  pre-run, and the user guide and the input reference describe it as one. Against the comparison tiers
  and the examples the estimate reads 20 to 40 percent low, inside the factor of two it is held to, and
  part of that is the source filter, which became the default after the performance page was measured.

- 2026-09-07 Changed: global%source_filter is on by default. The filter is measured on four machines, two
  of them real, over angle ratios of 1.14 to 1.86, and at a converged load it removes about a hundred
  times the wide-angle power while moving the mode power a few percent and the saturation point not at
  all. On the Aramis benchmark at the derived grid, 1024 macroparticles per slice with the filter and
  8192 in beamlets of four without it are both converged, at 7 and 12 s, and the filter's own cost is 8
  percent of wall clock there and 14 at the fine grid. fel-physics.md gains sec-convergence, one section
  carrying the whole story, and every other mention points to it. The filter reaches the elements that
  build a source to filter: an unaveraged element is left alone and the coherent source turns it off
  for the run, each with a line at setup, where both were refused before, since a deck that asks for
  the unaveraged mode has not asked for the filter. Nothing recorded moves: the tier template and every
  deck a check writes state the filter, off where a level was recorded without it and on where the
  filter is what is tested, 25 templates in 16 scripts, and the library twin states it in Fortran as
  it states the transport. Three of those were silent off cases that the new default would have turned
  on, in the device check, the load check and the startup-noise instrument, and each would have passed
  while measuring something else. The examples take the default and their exit lines move where the
  filter reaches. The unseeded SASE decks fall by up to a factor of ninety, sase from 3.02 GW to 33.6 MW
  and import from 3.5 MW to 42 kW, since at 128 beamlets on a 96-slice window almost all of the
  reported power was the artifact, and their READMEs say so. The seeded and steady-state decks move by
  one to five percent. The two real machines move inside their seed spread, FLASH1 from 46.35 to 46.72
  uJ and LCLS from 255.3 to 252.8 uJ, and their four-seed tables stand as measured. Three decks state
  the filter off because they compare against Genesis4, which carries none, and three because they
  compare the two FEL methods like for like. The coherent-source example's trap shrinks with the
  filter on, a low load costing 39 percent rather than a factor of 7.6. The migration example measures
  its comparison over four seeds instead of one: the exit power ratio between migration on and off runs
  from 0.30 to 1.54 and settles nothing, while the exit bunching is higher with migration in all four
  seeds, by a factor of 1.41 +- 0.30. The source-filter section gains a check on where the filter
  reaches, run rather than read off the line the run prints, and it found a coherent-source run arming
  the per-element filter with no edge set and moving its exit power by 3.3e-9 (FINDINGS 7.60).

- 2026-09-07 Fixed: the FP32 lockstep instrument allocated its field record on first use inside the
  per-slice loop, which carries an OpenMP parallel do, so two threads could pass the allocation test
  together and both allocate the same array. A bounds-checked build stops there, on about one run in
  five of the eight-slice deck under load, and a one-slice window never raced at all. The record spans
  the window rather than a slice, so fel_fp32_field_prep now allocates it before the loop and the
  per-slice routine writes only its own page. The FP32 check prints a failed run's stderr beside its
  stdout, since a compiler runtime error is written there and two keystones reported a run that
  stopped after its banner with nothing to read (FINDINGS 7.61).

- 2026-09-06 Added: the slice spacing and the integration step are derived from the gain, and the
  default window holds the line's slippage. A deck that states no slicing%n_wavelength gets the whole
  number of wavelengths in a quarter of the cooperation length, four slices per cooperation length,
  which is where FLASH1's pulse energy converged. On FLASH1's own beam that derives 10 wavelengths
  against the 12 its example states by hand, and over four seeds the two give 50.0 +- 4.2 and
  50.8 +- 6.8 uJ, 1.4 percent apart inside a seed spread of 8 to 14 percent. Two descriptions are not
  resolutions and keep one wavelength: the steady state, where the whole charge sits in the one slice
  so the spacing states the current, and a run whose beam is loaded rather than described. An element
  that states neither ds_step nor num_steps gets a twentieth of the one-dimensional gain length in
  whole periods. What that replaces is not silence but Bmad's own default for a wiggler, a twentieth
  of a period, which resolves the wiggle the averaged model has already averaged over: the FLASH-scale
  probe of doc/startup-noise.md ran at 3300 steps a segment against FLASH1's 55 for that reason, and
  at the derived 165 its two instrument runs take 8 and 9 s rather than 67 and 85. Against half the
  derived step at matched seeds the pulse energy moves 1.8 and 2.1 percent. The default window holds
  every particle of the bunch and one wavelength per undulator period ahead of the window head, so
  radiation that slips forward has somewhere to go. A Gaussian has no last particle, so it spans the
  length beyond which a slice would hold less than one electron of the charge, a rule that scales with
  the charge where a fixed number of sigmas does not. A loaded bunch spans its own particles and gets
  no slippage added, a loaded window being data rather than a choice. No example moves, since every
  time-dependent deck states its own window and spacing. The second machine's numbers are re-recorded
  at the derived step (FINDINGS 7.59).

- 2026-09-06 Changed: an FEL segment is a wiggler or undulator carrying Bmad's fel_method attribute,
  which replaces the two FEL tracking methods. tracking_method names how one particle crosses an
  element, and multiparticle physics riding on it has its own per-element attribute in Bmad, as
  space_charge_method and csr_method do. Encoding the FEL step as a tracking method made every
  lattice here unloadable in Tao: the released Tao refused the switch at parse and this tree's Tao
  parsed it and then died in track1 on an unset custom pointer. All twenty-five lattices now carry
  fel_method and leave tracking_method at bmad_standard, and a new benchmark section runs Tao over
  every committed lattice and requires a zero exit with no error line. Twenty-one load. The
  unaveraged mode's two numbers stop being attributes this program registered, which the released
  Tao also refused, and become global%unaveraged_steps_per_period, 20 and refused below 10, and
  global%unaveraged_ramp_periods, 2 with -1 the explicit hard edge. Averaged and unaveraged elements
  still mix per element. Two refusal checks changed with the encoding, since an ordinary Bmad
  wiggler meets Bmad's own sanity checks first: a missing l_period is now Bmad's refusal naming the
  attribute, and the field-map case carries a real cartesian map, because Bmad accepts only the two
  analytic field models without one and a map-free element would be refused for lacking a map rather
  than for having one. All eleven tier digits are unchanged, since only recognition and plumbing
  moved.

- 2026-09-06 Fixed: tests/bmad/flash.bmad holds the beta function its header claims. The probe stated
  beta = 10 m with zero alpha, which is not the matched periodic solution of its own cell, and at
  quadrupole strengths of plus and minus 2.5 m^-2 the cell held a mean beta of 8.19 m rather than the
  10 m the header named. The beam beat through the line at a mismatch of 2.3, its rms size swinging
  from 27 to 256 um. The strengths are now plus and minus 1.185 m^-2, which hold a mean of 10.0 m in
  both planes with beta between 6.74 and 13.91 m, and the beginning statement is that cell's matched
  solution. A helical undulator focuses both planes alike and the cell's two halves differ only in
  which quadrupole they carry, so one strength sets both means. The run confirms the match: the beam
  returns to the stated Twiss at every cell end and its rms size holds 78 to 116 um. Nothing else in
  the tree reads this lattice, so no recorded digit moves. doc/startup-noise.md re-records the second
  machine on the matched beam. Its geometry barely moves, the mode angle going from 65.2 to 64.5 urad
  and the ratio of the two angles from 1.87 to 1.86, while its dynamics move a great deal: the power
  inside the mode at the exit goes from 1.23 to 1.88 GW, the bunching at 43 m from 0.062 to 0.127, and
  the saturation point from 27.1 to 21.6 m. The filter's default holds as it did, removing 171 times
  the wide-angle power rather than 165 and raising the mode power 19 percent rather than 6. The
  four-machine tables take the new rows, the range of ratios the default is measured over becomes 1.14
  to 1.86, and the coefficient of the angle rule is unchanged at 6.8 with its threshold at 46
  (FINDINGS 7.57).

- 2026-09-06 Changed: the branch takes Bmad main at 20260904-1, 158 commits since the last common
  point. The two FEL tracking methods, the wiggler-averaged and the unaveraged, and the openPMD
  wavefront work in the HDF5 layer merged without conflict, and the parser hook that gives an FEL
  element its custom pointers merged with them. Nothing under the wiggler, the undulator, the Twiss
  or the field solver moved upstream, and all eleven tier digits and the 27 check sections land where
  they were recorded. The regression suite reports 53 passed and 3 skipped, one more than before:
  upstream added exact_bend_edge_test, which passes. One upstream change reaches a path this branch
  runs, the multipole conversion, which now recomputes rather than reading its cache, and no lattice
  here carries a multipole.

- 2026-09-06 Added: examples/lcls, the LCLS undulator line at SLAC at 1.5 Angstrom from its published
  parameters, and with it the fourth machine the filter's derived default is measured on. The lattice
  is thirty-three planar segments of 114 periods at 3.0 cm and 1.25 T with a quadrupole in each break,
  and the beam is a flat 3.0 kA slab at 0.4 mm mrad and 13.6 GeV, all from P. Emma, PAC09 TH3PBI01,
  and Emma et al., Nature Photonics 4, 641 (2010). The field and K agree to 0.04 percent and the
  resonance that follows is 1.5099 Angstrom against the published 1.5. The paper gives the maximum
  quadrupole gradient and not the operating one, so the strengths are chosen to hold the published
  30 m mean beta, and they land at 96 percent of that maximum. Over four seeds the run gives a power
  gain length of 3.81 +- 0.04 m against the measured 3.3 m and the design 4.5 m, saturation at
  64.7 +- 1.9 m against about 60 m, and a peak power of 14.1 +- 1.0 GW in one slice against 15 GW,
  which at the published 75 fs pulse duration is 1.06 mJ against 1.1 mJ. The deck states no grid and
  takes the derived one, 127 points over 191 um. The page gains the fourth machine's section and two
  tables over the four machines. The default holds a fourth time and at a converged load is nearly
  free, removing 105 times the wide-angle power while moving the mode power 4 percent and the
  saturation point not at all. The ratio of the two angles the default is built from goes as
  1/sqrt(z_R/L_g) to 4 percent over the four, so the rho angle wins only above z_R/L_g of about 46,
  which none of the four reaches and no hard X-ray line is built to: the default is measured over the
  range machines occupy. The condition on the in-cone startup comparison is settled with it, since
  over the four the measured ratio rises with the first segment's length in gain lengths, from 0.81 at
  1.2 to 2.6 at 5.8 (FINDINGS 7.56). tests/scripts/startup_noise.py gains --machine lcls and a
  per-machine long window, since twelve wavelengths per slice do not reach the slippage at 1.5
  Angstrom.

- 2026-09-06 Added: every run reports the split between the mode and the wide-angle emission of the
  point macroparticle beamlets, and derives its transverse grid from the beam when the deck states
  none. The split costs nothing new: each element end that takes the field angle moments already
  stores the power inside `field/total/split_angle` beside the total, and the footer reads it back
  at the record where the mode peaked and at the last record. Above a tenth the run warns that a
  total power quoted from it is mostly representation and names the two remedies, more
  macroparticles per slice at the same `beamlet_size` or `global%source_filter = T`, and below a
  hundredth it says the total is converged. The same numbers go to `run/` in the statistics file.
  Eleven of the sixteen examples exceed a tenth, the dark starts worst at 21 to 86, and every
  README that quotes a power now says what share of it is that emission. Where the comparison an
  example exists for survives the split, the README states it in the mode as well: migration is
  worth 1.75 there against 1.73 in the total, and the mixed line costs 3.6e-2 in ln P against
  7.4e-3. `wavefront_init%grid_n_pts` and `%grid_half_width` now default to zero and derive from
  the rms beam size, cells of a seventh of it and a half width of nine of it, which fixes the point
  count at 127 since the beam size cancels, rounded up to a power of two on the device. Both are
  printed with their origin. The derivation reproduces the grids chosen by hand on the FLASH probe
  and on FLASH1, and on the Aramis benchmark it gives 127 points over 192 um against the examples'
  255 over 200, which lowers the wide-angle ratio of the `sase` deck from 21 to 5.4 (FINDINGS 7.55).
  The header's Radiation line now reports the grid the run built rather than the one the deck
  stated, which a run starting from a field file did not have.

- 2026-09-06 Added: examples/flash1, the FLASH1 undulator line at DESY at 13.7 nm from its published
  parameters, and with it a third machine for every numerical default. The lattice is six planar
  segments of 165 periods at 27.3 mm and 0.47 T with a quadrupole doublet in each intersection, and
  the beam is 2.5 kA over a 30 fs spike at 1.5 mm mrad, all from Ackermann et al., Nature Photonics
  1, 336 (2007) and Schreiber and Faatz, High Power Laser Science and Engineering 3, e20 (2015). The
  two published numbers for the undulator disagree by 2.7 percent, since 0.47 T at 27.3 mm gives
  K = 1.19807 against the stated 1.23, so the field is taken and the resonance at 13.7 nm sets
  668.494 MeV against the machine's stated 700. The lattice header marks that and the three numbers
  the file chose, which are the whole-period segment length, the intersection length and the doublet
  strengths. Over four seeds the run gives a field gain length of 2.05 +- 0.16 m against the measured
  2.5 +- 0.3 m, an exit pulse energy of 50.3 +- 3.1 uJ against 40 uJ characterized and 70 uJ
  demonstrated, and the end of exponential growth at 20.9 +- 1.1 m inside the 27 m line. The README
  carries one row per numerical setting with the Aramis value, the value this machine needs and the
  measurement behind it. Two of the SASE convergence page's numbers do not carry to it: the estimate
  of the wide-angle share is 40 times low, and the agreement of the in-cone startup power with the
  spontaneous emission holds only where the first segment is short in gain lengths, 2.2 on Aramis
  against 5.8 here (FINDINGS 7.54). The page gains the third machine's table, tests/scripts/startup_noise.py
  gains --machine flash1 and the undulator coupling factor that a planar device needs, and
  tests/bmad/flash.bmad says in its header that it is a probe at FLASH's scales and not FLASH.

- 2026-09-06 Changed: one path loads the beam. A bunch comes from `beam_init`, generated or
  read through `beam_init%position_file`, and `load_mode` says what becomes of the particles in
  a slice: `"sample"` loads `beam_init%n_particle` macroparticles into every slice with weights
  from the slice's charge, and `"keep"` keeps every particle as a beamlet of `beamlet_size`
  copies at weight over `beamlet_size`, so every charge-weighted moment of the bunch is the
  load's, measured at 1.5e-11 over the slices. `dist_file` and `use_beam_init` are retired and
  refused with the names that replaced them, `beam_init%position_file` and
  `resample%use_beam_init`, the latter the resampler's validation route. Two switches act on
  every load. `quiet_start`, default true, makes beamlets of copies with their phases spread
  over 2 pi. `shot_noise` imposes Fawley's noise on whatever load there is, by one phasor per
  group of macroparticles sharing their transverse coordinates: the beamlets when the loader
  made them, else particles sharing their coordinates exactly, else the occupants of a deposit
  cell, and the message says which. The loader measures the quiet floor first and refuses noise
  above 0.01 in |b|^2 N_lambda with the value in the message, so a bunch that carries its own
  noise is never given a second dose: the unquiet generator sits at 0.28 and a quiet bunch made
  elsewhere at 3e-28. The resampler honors `shot_noise` where it applied the noise unconditionally
  before, so its decks now say `shot_noise = T`. `beamlet_size` is a loading quantity: the coherent
  source counts the slice's distinct transverse positions instead, so a dump loads without it. A
  generated Gaussian bunch in sample mode radiates 1.385e5 W per 3 kA slice inside 3 urad at the
  end of the first Aramis undulator, against the 1.33e5 W of the flat window's interior. The
  `import` example moves to the path, the keystone's eleven tier digits hold since the tiers load
  dumps, and check_load.py joins the harness.

- 2026-09-05 Changed: the time window moved out of `wavefront_init` into `slicing`, a struct of
  its own in `&fel_params`, and the wavelengths per slice are carried as the integer they are.
  A wavefront is an optical object with a longitudinal spacing and no slices; the rule that the
  spacing is a whole number of wavelengths belongs to the FEL interaction, since it is what lets
  the field record rotate by one index with no interpolation. `slicing%n_wavelength` replaces
  `wavefront_init%window_sample` and `slicing%window_length` replaces
  `wavefront_init%window_length`, both retired names refused with the name that replaced them.
  `slicing%n_slice` states the window as a count instead, and `slicing%current` states a flat
  current directly, which needs no z description: with a window it is a flat time-dependent run
  and with none the steady state of one slice, so a steady-state deck no longer computes a
  bunch charge and a z grid to say what a current says. The integer reaches every consumer from
  the deck rather than being recovered by dividing a spacing in metres by a wavelength, which
  returns 12 to the last bit rather than 12 and which the device refused: a 10 nm deck that ran
  35 minutes on twelve CPU threads now runs on the device. A beam or field from a dump is the one
  place the integer is recovered, and a spacing that is not whole is refused there.

- 2026-09-05 Added: doc/startup-noise.md, the dependence of SASE power on the transverse cell
  size and the macroparticle count, with tests/scripts/startup_noise.py producing every number and
  figure on demand. The field power of an unseeded run has two parts. The power inside the mode, taken from
  the far field of dumps at five undulator ends inside 3 urad, is the same for cells of 6.35 to
  0.78 µm and loads of 1024 to 65536 particles per slice within the SASE fluctuation, 1.8e9 to
  3.9e9 W per interior slice at z = 37 m, saturating at 28 to 33 m beside Saldin, Schneidmiller and
  Yurkov's 30.5 m, and at the end of the first segment it equals the beam's physical spontaneous
  power inside that cone, 1.3e5 W per slice. The power outside the mode is the point beamlets'
  own emission: independent of the particle count at first and between 1/dx and 1/dx^2 in the cell size, then
  growing with the bunching as one over the beamlet count and independent of the particle count
  at a fixed beamlet count. At 1.57 µm cells and 1024 particles it is 4.2e9 W per slice at the
  exit against 8e8 W inside the mode. The examples' 96-slice window is 1.9 cooperation lengths
  long and its power is that emission at every particle count measured, 7.3 GW at 1024 particles and
  0.17 GW at 65536 at 1.57 µm, so the SASE examples' recorded powers are numbers of that kind
  and their ratios stand. Genesis4 4.6.15 on the SASE tier input reports 2.26e7, 1.17e8 and
  8.16e8 W at 151, 255 and 511 points, and Lucifer from its dumps lands on it at 1.9e-6 at each.
  The page states what the codes model, what converged, and how to choose the grid, the window
  and the particle count. The index gains a row, validation.md's SASE tier paragraph a sentence, and the
  sase, migration and sase_wake examples one sentence each. No setting or recorded number moved.

- 2026-09-04 Added: slice migration runs with the device. The pass is fel_migrate_slices
  unchanged, serial on the host at an FEL element's last step, and it touches no field, so no
  kernel is involved: with the beam resident at that step the pass reads the beam and the field
  set back first and then releases residency through a new fel_device_release, since the host is
  now authoritative and the field came back with it. The comb readback and the element end that
  follow find nothing resident, so the readback moves rather than doubles, and the next element
  uploads the re-sliced beam. The device's particle buffers are a rectangle sized at setup to the
  largest fill, and migration grows fills, so the seam gains luc_dev_resize_particles: growth
  only, no shader recompile since the kernels take the count at dispatch, called from the element
  entry and the twin's staging when the largest fill exceeds the capacity, and reported in the
  run's output. The device's refusal of migration leaves fel_setup_mod. Fluctuations with
  migration stays refused, the CPU's own rule. Measured on the debug build with check_migration's
  decks at grid 64: the heavy-migration run bites past 73,000 moves on the device with charge
  conservation at 1.6e-14, phase continuity at 5.2e-15 and every survivor inside its window, its
  exit power 4.0e-4 from the CPU's, and the capacity grows under it. The lockstep instrument on
  the migrating deck, the twin re-staging every step across the changed fills, holds x..py at
  2.9e-7, theta at 7.4e-6, the phasor at 2.0e-7 and source and field at 7.9e-6, with the FP64
  diag byte-identical, twin on against off. The no-op is self-referenced, since the deposit's
  atomics rule out byte identity on the device: zero moves, and the final beam 1.7e-14 slice
  spacings from the migrate = F device run. One instrument property surfaced and is recorded in
  validation.md: on a quiet-start deck without shot noise every slippage feeds the tail slice a
  field at FP64 roundoff, and the source and field rows, normalized by the slice's own field,
  read 0.4 there on a deposit whose end-to-end power agreed at 4e-4, so the check's deck carries
  a seed and shot noise. examples/migration runs on the device at grid 256, agreeing with the CPU
  at 2.8e-4 on exit power, 148,654 moves against 148,652 with the dropped charge identical, the
  capacity grown from 1024 to 1215 per slice along the line.

- 2026-09-04 Added: harmonic field sets and two polarizations run on the device. The resident
  field is now the run's whole set, every harmonic member and both polarization planes, stored
  member-major with one propagator per member. The step's constants carry each member's harmonic
  number, coupling fc(h) and deposit scale plus the element's polarization pair. The push gathers
  every member at h*theta into every RK stage, fel_ode_multi's structure in the FP32 reformulation,
  the deposit writes every member's source at exp(-i h theta), and the transform runs over every
  plane with its member's propagator, the polarization factors applied at the field add as
  fel_field_step writes them. Every member rides the one fixed-point phase at h times it and every
  plane slips together. With one member and one plane every kernel is the single-field arithmetic
  it was, and the eleven tier digits are byte-identical on both builds. Harmonics together with
  two polarizations stays refused, for the device as for the CPU. The field-set types moved to
  fel_field_mod so the device seam can take the set as ff(:): passing ffield(:)%wf, a component
  section of a pointer array whose elements hold allocatable arrays, made the compiler pack a
  temporary and deep-free it on return, a heap abort with no Fortran message (FINDINGS 7.46).
  The device twin judges every member: phasor, source and field rows are the worst over members,
  the footer adds each member's own worst, and every member's source and field rows are
  normalized by the fundamental's field, since a quiet-start harmonic's own norm is FP64 roundoff
  on the CPU (FINDINGS 7.47). Measured on the debug build: harmonic lockstep phasor_h3 6.6e-8
  steady and 2.0e-7 time dependent, three times the fundamental's phase error as h = 3 predicts;
  production P3 against the CPU 6.0e-4 on a strongly bunched beam and 3.4e-6 on a shot-noise
  window. The one-step Bessel identity P3/P1 is 4.4e-6 in FP32 against 2.7e-14 on the CPU.
  Against one time-dependent Genesis reference with shot noise at grid 64, imported by CPU and
  device alike, the device lands at 1.5e-5 on P1 and 2.3e-5 on P3 where the CPU lands at 8.4e-7
  and 1.1e-6. The crossed undulator on the device: Px 1.4e-4 and Py 4.3e-3 against the CPU steady,
  4.1e-6 both with slippage live, x-field isolation through the y set 4.3e-6, Py/Px 8.8e-3 above
  the 5e-3 floor. One limitation is recorded rather than hidden: the device resolves harmonic
  bunching to about 2e-7 of the charge, and a quiet-start dark harmonic whose true bunching is
  1e-8 radiates that floor instead, 137 times the CPU's P3 on the planar steady-state tier deck
  while P1 holds its band. On the production build, the 96 x 8192 SASE window at ngrid 256 on the
  planar segment runs 12.4x faster than twelve CPU cores with the fundamental alone, 10.4x with
  harmonics 1 and 3, and 10.8x with two polarizations. examples/harmonics and
  examples/crossed_undulator run on the device at grid 256, agreeing with the CPU at 4.7e-4 on P3
  and 9.5e-4 on Py.

- 2026-09-03 Added: the checks run on a second platform, in CI. `lucifer-keystone.yml` builds the
  debug tree on Ubuntu against a conda environment and runs the benchmark harness and the wavefront
  validation there, on demand and on a weekly schedule. The job takes 37.5 minutes on a
  four-core runner, of which the harness is 30 and the build 6.4, against the local keystone's 7.7
  minutes on sixteen cores with its references cached. The defect class it reaches is one platform
  hiding what another shows, and this project has paid for an instance: five wake cases failed on
  Linux because a steady-state run read past the end of a one-element vector, and the macOS
  allocator had been answering with a number that happened to work. The debug build is the point,
  since `-fbounds-check` turns that read into a failure rather than a wrong number. The run asserts
  the tolerances and not the digits, because a different FFTW and HDF5 give different last digits by
  construction. It varies the operating system, the C library, the allocator and every dependency
  build, and not the compiler project, which is conda-forge gfortran on both sides. Two changes
  made it possible. `lucifer` run with no arguments now names
  the device backend the build carries, either the device it found or the reason there is none, and
  the harness reads that line to decide whether the device section can run, so the skip keys on the
  build and the machine rather than on the platform and a machine that has a backend cannot skip.
  And `report_validation.py` refuses a results file whose sections did not all run, since the table
  it writes lists every section that ran and a file recording the skip is one short. The reference
  pin is what unblocked the rest: `genesis4` 4.6.15 is on conda-forge for linux-64, so a runner
  installs the same reference the levels were recorded against. The first green run is the
  demonstration of what the job does and does not hold: all eleven tiers passed their tolerances,
  four of them on digits identical to the recorded ones and seven moved in their last places, and
  the device section skipped with its reason recorded, on a build whose banner reported no backend.

- 2026-09-03 Changed: the everyday keystone runs in 7.7 minutes rather than 25, with nothing behind
  a flag and every section still running in every run. Three measured changes, in order of what they
  bought. The Genesis references are cached: they are a pure function of the reference binary and the
  decks that make them, so they sit under a content-hash key naming both and the pinned version
  besides, and each run reports hit or miss. A set is staged privately and moved into place with its
  completion marker written last, since the two build passes now run at once and a reader must never
  find a half-filled directory. Everything independent runs concurrently, which is the two benchmark
  passes in separate work directories plus the regression suite, the wavefront validation and the
  examples: they share only a source tree they read, with nothing persisting FFTW wisdom or writing
  outside its own directory. And the spontaneous section's acceptance scan moved an octave down to
  ngrid (63, 127, 255), keeping 255 as the reference point so its fraction check is untouched and the
  shape test still spans the same 16x solid-angle range; that section goes from 114 s to 77 s on the
  debug tree, and the scaling test is better centred than before (measured against predicted 1.07,
  where the old bracket's top point sat at 1.31 because the dipole model it compares to is furthest
  from valid there). The per-job thread count and the job pool's width moved together, from eight
  threads two wide to four threads four wide, which is the same subscription with more overlap. Every
  tier digit is unchanged on both builds, and the deliberately single-threaded job in that section
  was left alone on purpose: it also feeds the energy-ledger closure, whose recorded level is a grid
  quantity, so a cheaper deck there would move a measured number rather than only the clock.

- 2026-09-03 Changed: the reference is the Genesis4 4.6.15 release, and six of the eleven tiers are
  re-recorded. 4.6.15 is the first release carrying the CODATA electron rest energy, which is what
  the retired `$GENESIS4` variable existed to pin, so the reference is now `genesis4` from
  conda-forge in `bmad-fel-validate` and no export is needed. It also carries this project's
  rank-independent seeding, and that moves noise realizations: `ShotNoise` is re-keyed from the
  global slice index rather than drawn continuously per rank, `Incoherent` holds one generator per
  slice, and `RandomU::set`, which had an empty body, now re-keys as the surrounding code always
  intended. The digits move exactly where that reaches, which is the attribution: the four
  steady-state tiers (`shotnoise = 0`) and `weight_split` (Fortran against Fortran) are
  bit-identical, and the six loading from a shot-noise deck moved. td1 8.467690e-07 ->
  8.362073e-07, td2_genesis 2.398226e-06 -> 2.376734e-06, td2_bmad 4.127587e-02 -> 3.442186e-02,
  tdsase 2.292496e-06 -> 2.349627e-06, tdsc 2.440477e-04 -> 2.868371e-04, tdwk 8.708129e-07 ->
  8.917019e-07, all on the debug tree, with tier1 1.825899e-06, tier1_unavg 6.934613e-02,
  tier2_genesis 1.771890e-05, tier2_bmad 5.001254e-02 and weight_split 3.532394e-13 unchanged. The
  production tree lands at the same digits except where the build reaches: tier1_unavg
  6.934017e-02, td2_genesis 2.376736e-06, tdsase 2.349843e-06, weight_split 3.536001e-13. Every
  tolerance is unchanged and every level stayed in its own order, between 0.83 and 1.18 of itself.
  The harness now refuses a reference reporting any other version, printing what it found and
  wanted, since the old failure mode was a silent re-baseline.

- 2026-09-03 Changed: the device readback runs parallel over slices. At the finest stats comb
  (`comb_ds_save = 0`) every row paid ~4.3 ms of serial FP32-to-FP64 conversion on the 96 x 8192
  case, half the row's cost. The conversions are elementwise per slice, so the loops now run under
  OMP with block-local staging, after one serial drain of the device: the seam's transfers are then
  pure disjoint copies, a contract `lucifer_device.h` now states. Measured: the fine-comb walk goes
  from 1.231 s to 0.897 s, and the readback share from 0.409 s to 0.181 s. The share did not return
  to noise, and the thread scan says why: it plateaus by four threads at the memory-bandwidth floor
  of the bytes themselves, ~70 MB of transfers per record. That finding is recorded in
  `doc/performance.md` as the input to the separate device-side moments decision, which is the only
  route under the plateau and is not taken here. Bit-identical by construction (no accumulation,
  per-slice writes), and every tier digit holds.

- 2026-09-03 Fixed: two build-time capabilities are detected rather than assumed, so builds whose
  toolchains lack them compile and refuse the knob instead of failing to build. The single-precision
  FFTW the FP32 field twin transforms with is present in the conda environment and absent from the
  off-site distribution, whose FFTW is built from source in double precision only, where an
  unconditional `-lfftw3f` does not link; the build now looks for the library and `fp32_check`
  refuses without it, since transforming in another precision would measure something other
  than what it reports. The Metal backend is ARC-managed Objective-C++, so it needs a Clang-family
  Objective-C++ compiler and not merely macOS: a build that hands `OBJCXX` to a GNU compiler
  rejected `-fobjc-arc` outright, and such a build now takes the refusing stub like any other.
  Configure output names both outcomes. A capability-free build was measured to produce
  byte-identical physics, and every tier digit holds.

- 2026-09-02 Added: the Metal device backend. `global%device = "metal"` runs the averaged FEL step
  on an Apple Silicon GPU: the beam and the field upload at each FEL element's entry and stay
  resident through it, every integration step is one command buffer (transverse maps, the RK4 push,
  the deposit and the FFT field solve, all FP32), and the state returns at the stats comb's
  positions and the element end. The seam is one C interface behind `iso_c_binding` with one
  Objective-C++ file behind it and a refusing stub on every other platform, shaped as the GPUEngine
  design of this project's own prior Genesis 1.3 v4 backends (branch `gpu/metal-engine`, 4919b01,
  unmerged upstream); the transform transcribes that backend and the physics kernels transcribe
  `fel_fp32_mod`'s priced reformulations. One divergence, stated at the citation: the longitudinal
  state is a 64-bit fixed-point phase in ticks of two pi over 2^32 off a static FP64 per-slice
  reference, so bucket wraps are exact integer arithmetic and are asserted exactly on the device at
  every setup. Acceptance is the lockstep instrument with the device in the twin's role: with
  `fp32_check` set the FP64 run is untouched and the device's rows land in the same stream, inside
  ceilings at three times the recorded CPU-twin levels (measured: theta 6.8e-6 rad, source 1.5e-5,
  guard 6.8e6 ticks), with `fp32_mutate` perturbing a kernel constant so a wrong shader moves a
  recorded level (710x on theta). End to end the freerun twin and the resident production run land
  at the same 5.1e-4 on exit power against the CPU. Everything the kernels do not cover is refused:
  harmonics, two polarizations, the coherent source, wakes, space charge, radiation,
  migration, the unaveraged mode, the escaped-field bank, and any grid that is not a power of two
  from 64 to 1024 (the message names the nearest size). Measured on an M3 Max against 12 CPU cores:
  7.3x steady, 8.2x to 12.4x time-dependent, with the 0.31 ms per-step dispatch floor stated and the
  deck below it left to the CPU. The harness gains a device section of thirty-six assertions, and
  every tier digit holds with the knob off.

- 2026-09-02 Added: the FP32 field solve joins the lockstep instrument, and the end-to-end single-
  precision run has its number. The field twin scatters the deposit into an FP32 source grid,
  transforms it with FFTW's single-precision interface (FFTW_ESTIMATE, the determinism decision the
  FP64 solver already records), multiplies the propagator rounded from the FP64 kernel, and two rows
  join `.fp32.txt`: the source grid and the post-solve field, both against the FP64 field's norm.
  Measured through saturation the field step costs under 1e-4, with the deposit's phase noise
  dominating and the transform beneath it. In freerun the FP32 field carries across steps, fed and
  read by the twin itself, so the end-to-end divergence of the exit observables is recorded: 1.0e-3
  on power over a gain segment (the published backends' comparable figure is 6.7e-4) and 2.7e-1 in
  deep saturation, where an FP32 trajectory ends at a different phase of the synchrotron oscillation
  and the number says so. Freerun covers a single-slice window and refuses more, since the twin
  keeps no slippage rotation. The harness section grows to twenty-one assertions including field-row
  falsifiability, and every tier digit holds.

- 2026-09-02 Added: the FP32 lockstep instrument. `global%fp32_check = "lockstep"` steps a single-
  precision twin of the averaged FEL advance beside the FP64 path from a shared state and records
  per-quantity divergence to `<out_root>.fp32.txt` ("freerun" carries the FP32 state so the
  compounding rate is measured too). The twin carries the reformulations FP32 forces: energy as an
  offset (the FP32 quantum of absolute gamma exceeds the per-step change), phase as a residual off
  an FP64 per-slice reference with the base rotator evaluated once per slice, and the ODE detuning
  as a difference of like small quantities (the FP64 form's square-root argument rounds to exactly 1
  in FP32). A runtime guard refuses a run whose residual cannot resolve the per-step increment, the
  silent failure mode, and it fired on the instrument's own first run: a static reference let the
  beam's secular z drift grow the residual until step 514 of 1068, so the reference now follows the
  slice mean. Measured on the way to saturation, FP32 costs about 4e-7 on the energy chart, 5e-5 rad
  on the phase and under 1e-5 on the source phasor per step. The FP64 path is untouched: the
  instrumented run's outputs are byte-identical to the uninstrumented run's, and every tier digit
  holds. Harness section `fp32-lockstep` records the levels and proves the check can fail.
  doc/validation.md carries the tables.
- 2026-09-02 Changed: the unaveraged push runs stage-major. Each RK4 stage is one loop over the
  slice's particles into per-particle stage arrays, with the combination a last loop, instead of
  four stages walked inside each particle through a call tree. Per particle the arithmetic is the
  same operations in the same order, so every result is byte-identical, and the unaveraged step is
  7.9% faster at 2048 particles per slice and 17.9% at 16384. Inlining the field into the stage loop
  was tried for vectorization and refused: the fused arithmetic differs from the called form in the
  last bits, and bit-identity is the contract. The stage loops therefore still carry a call and
  still do not vectorize, which doc/performance.md records with the reason.

- 2026-09-02 Added: `tests/scripts/sincos_audit.c` records why `code/fel_sincos.c` still calls libm.
  Replacing that call with a Cody-Waite reduction and the fdlibm kernel polynomials, which is what a
  hand-written double-precision sincos looks like, lands 6 ulp from libm on 34% of arguments where
  the shim it would replace was admitted at one ulp on 2e-6 of them. It also degrades with the
  argument, 10 ulp at 1e5 and 2681 at 1e8, and the argument carries the common phase accumulating
  along the line, so a longer line moves toward that with nothing to detect it. Measured gain 3.2%
  of an averaged run at 2048 particles per slice. Refused, and the audit is committed so the next
  reader starts from the number. FINDINGS 7.38.

- 2026-09-02 Changed: the unaveraged mode's undulator field no longer evaluates transcendentals per
  particle. `fel_unavg_bfield` took a position `s` and computed the ramp envelope, its slope,
  `cos(ku s)` and `sin(ku s)` from it, once per particle per RK4 stage, when all four depend only on
  the substep position and are therefore the same for every particle in the loop. They now arrive as
  arguments built once per push by `unavg_field_quartet`, and the routine no longer takes `s`. The
  unaveraged step is 14.8% faster at 2048 particles per slice and 27.6% at 16384, libm samples fall
  from 6936 to 1221, and every result is bit-identical: the stage positions are the expressions
  `unavg_push` evaluated inline, in the same order, so every value is the same bits it was. FINDINGS
  7.37.

- 2026-09-02 Fixed: `doc/performance.md` labelled a profile category "libm sin and cos" where it had
  summed every libm transcendental, `exp` and `cexp` and the shared reduction helper included. `sin`
  and `cos` proper were 13.9% of the unaveraged mode rather than the 53.6% the table implied. The
  category is now named for what it sums and the sub-split is stated.

- 2026-09-01 Added: every run's footer prints where its time went. `code/fel_timer_mod.f90`
  accumulates wall clock per phase at region boundaries, always on, and the phases partition the
  walk so the fractions sum to it and the remainder is a row named `unaccounted` (measured, 0.02%).
  `tests/run_perf_benchmark.sh --phases` runs the profile at two slice counts and sweeps the thread
  count, needing neither genesis4 nor MPI. `tests/scripts/vec_audit.sh` reports what the vectorizer
  did to the hot files, taking its compile command from the production build's own makefile.
  `doc/performance.md` carries the measured tables, attributed to the machine and the build that
  produced them. No physics changed and no recorded digit moved.

- 2026-09-01 Changed: `doc/validation.md` no longer names the slippage rotation as part of the
  serial cost that caps thread scaling. Measured, slippage is 0.4% of the walk at 96 slices and 0.1%
  at 504, and the per-slice diagnostics reduction it named beside slippage runs slice-parallel. Its
  claim that the 32-slice curve is the floor of the scaling holds: at 96 slices in a production
  build, 8 threads give 6.76x against the recorded 3.97x, and 12 threads give 9.16x.

- 2026-08-31 Changed: the last of the prose-style debt in the Fortran. 203 single-word
  stress capitals in comments are lowercased per prose 1.4, and the 12 remaining
  `deliverable N` citations are gone, so nothing in the tree cites a document outside the
  repository. Untouched by design: `out_io` message text, which 6.3 sanctions as capitals
  and which refusal checks match; `!$OMP` directives, which are not prose; Bmad attribute
  registration strings; acronyms and units; and `SUM` where it transliterates a sum.

- 2026-08-31 Fixed: no machine-local paths in the harness. The genesis4 binary, the Python
  interpreter and the openPMD-beamphysics checkout were defaulted to paths under one
  developer's home directory, in nine files. Each now resolves through an explicit flag,
  then an environment variable (`$GENESIS4`, `$LUCIFER_PYTHON`, `$OPENPMD_BEAMPHYSICS`),
  then a portable search, and the refusal says how to fix it. `convert_genesis.py`
  lost its home-directory default and states its import contract: installed, or a checkout
  named explicitly.

- 2026-08-31 Added: `tests/scripts/genesis_constants.py` reports which electron rest energy
  a genesis4 binary carries, and the harness prints it beside the binary's path. Genesis
  releases to v4.6.14 use a value 2.14e-7 above CODATA; master carries the CODATA value,
  which is also Bmad's `m_electron`. Comparing against a release loosens the transcription
  tiers by several percent of their own value while still passing every tolerance, so the
  reference is now identified rather than assumed.

- 2026-08-31 Fixed: two check scripts ran genesis4 with a stripped environment, which a
  conda-provided MPI build cannot initialize in. They inherit the environment now.

- 2026-08-31 Changed: The documentation cites published work rather than code this port
  never used. MINERVA's papers informed the unaveraged mode, its source was not read, and
  no run of it was compared against, so the claim that it served as a statistical
  reference is gone along with a parenthetical naming one of its source variables. The
  ten-substep floor rests on this tree's own coupling-factor convergence, measured at ten,
  twenty and thirty steps per period, and both refusal texts say so. Puffin is removed: it
  was named once and nothing here traces to it. openPMD-beamphysics is described as what
  it is, a library this code calls.

- 2026-08-31 Fixed: three datasets under `params/resample` were still written as `npart`,
  `nslice` and `slicewidth`. The rename had changed the Fortran accessors but not the
  strings naming the datasets, so a stats file read back in the old vocabulary. The
  conformance check gained a third assertion, mutation tested like the others: a dataset
  written under `params/<family>` must name a component of that family's struct.

- 2026-08-31 Changed: The input names say what they are. `imp%` is `resample%`, with
  `n_particle_per_slice` (which is what it always was), `n_slice`, `slice_width` and
  `beamlet_size`. `wake%` is `chamber_wake%` and `sc%` is `space_charge%`, so neither
  collides with Bmad's own element wakes or `space_charge_com`. Bare `nbins` is
  `beamlet_size`, `shotnoise` is `shot_noise`, `write_opmd_file` is
  `write_openpmd_file`, `write_dist_file` is `write_genesis_dist`, and
  `imp_split_weights` is `resample_split_weights`. The stats file's `params/` subgroups
  follow, so a file reads with the vocabulary its input was written in. Genesis4's own
  names live in the translation table of `genesis4.md`.

- 2026-08-31 Added: `chamber_wake%model` and `space_charge%model`, both `"genesis"`, both
  refusing any other value. Each names the transcribed solver that runs at slice
  granularity, so a second implementation can arrive as a value rather than a rework.

- 2026-08-31 Fixed: `chamber_wake%write_kernels` works as written. The bare
  `write_wake_kernels` and the unconditional assignment that silently overwrote the
  struct path are both gone.

- 2026-08-31 Removed: the 104-line mapping of retired flat-group parameters to their new
  homes. The group name is still refused, since Fortran ignores an unknown namelist group
  in silence and a stale deck must fail loudly.

- 2026-08-31 Changed: The input reference reads as a reference. One section per namelist
  group with a subsection per family, each opening with a table of parameter, default and
  one-line meaning, and the long descriptions as anchored paragraphs below so another page
  can link a single parameter. The comment-laden namelist fences are gone, each group
  keeping one short uncommented example.

- 2026-08-31 Added: `tests/scripts/check_input_reference.py` holds the page to the struct
  declarations in both directions: every stated default must match the code, and every
  component of the six fel-owned input structs must be named on the page. It found 23
  parameters the reference had never documented, the whole of `wake%` and `sc%` among them,
  and one name that does not exist. Wired into the harness as its own check section.

- 2026-08-30 Added: The measured levels and the API reference are generated. The benchmark
  writes a results file under `--results`, `tests/scripts/report_validation.py` turns the
  debug and production results into `doc/generated/validation-measured.md`, and
  `validation.md` includes it in place of a hand-maintained tier table.
  `tests/scripts/report_api.py` reads the 217 `!+ ... !-` headers the source already
  carries, the convention `util/getf` reads, into `doc/generated/api.md`. Both are
  deterministic, so the keystone regenerates them and requires an empty diff: a moved
  digit is now a failing command.

- 2026-08-30 Changed: Development history leaves the user documentation. The user guide's
  units paragraph stated the convention and the comparison floor instead of narrating how
  the code got there, and "retired", "formerly" and "used to be" go with it. The mutation
  records in `validation.md` stay: what was broken deliberately and how loudly a check
  noticed is evidence about the check. Incidental Genesis4 source references leave the
  user guide and the input reference, and the manual's move into its Provenance notes.
  Every provenance citation is kept.

- 2026-08-30 Added: `genesis4.md`, one page for everything about Genesis 1.3 Version 4
  (Genesis4). Shared physics, what each code has that the other does not, a settings
  translation table, the conventions that differ, and file exchange. It pins the Genesis4
  release it was checked against, so it is the only page to revisit when Genesis4 changes.
  The "Facts about Genesis this work pinned down" section moved here from `validation.md`.

- 2026-08-30 Changed: The documentation presents Lucifer as its own code. Genesis4 is named
  in full at first mention on each page and `Genesis4` after, in place of five spellings.
  The introduction replaces its "In:"/"Not in:" paragraphs with Features and Known missing
  features. The claim of being faster than Genesis4 leaves the introduction, the user guide
  and the README: it was one configuration on one machine, and it now appears only in the
  two places that carry the measurement, each naming the machine.

- 2026-08-30 Changed: The physics manual no longer carries measured levels. Its 17
  Validation admonitions are removed and each section links to the measured record instead.
  The 12 Provenance admonitions stay, since transcription provenance is not comparison.
  Seven levels that lived only in the manual moved to `validation.md` and
  `reading-output.md` with their stories.

- 2026-08-30 Changed: The LaTeX manual is removed. `fel-physics.md` is the manual now, and
  `fel-physics.tex`, its `Makefile` and every reference to either are gone. The Markdown
  was kept beside the LaTeX for one commit so the two could be compared, and the
  comparison is done.

- 2026-08-30 Added: An introduction page and a references page. `index.md` is the site's
  front page: what Lucifer is, the three tracking methods, where to start for each kind of
  question, and what is in and out. `references.md` collects what the code implements and
  where each thing came from, every entry naming both the source and the routine that
  carries it. The manual keeps its citations in place beside the equations.

- 2026-08-30 Changed: The unaveraged mode is documented as a peer of the averaged one, not
  as a probe. The manual section retitles from "The unaveraged verification mode" and its
  opening states what the cost per step buys: the full quiver dynamics, the energy
  accounting the beam actually pays, polarization-agnostic coupling and arbitrary harmonic
  content. Its role as an independent check on the averaged path is still there, as a second
  paragraph rather than as the mode's name. The physics of the section is unchanged, and
  "verification mode" is swept from all eleven places it survived.

- 2026-08-30 Changed: The manual sheds the file format and the program. `fel-physics.tex`
  goes from 1580 lines to 1384. The diagnostic-output section was 213 lines of layout,
  attribute vocabulary and reader rules that `BMAD-STATS-SPEC.md` and `reading-output.md`
  own normatively; it is 117 lines of what the tracker computes and how each quantity is
  defined. The program, comb and tracking-window sections are dissolved into
  `user-guide.md`, which already held that content, and the fourteen pointers that named
  the deleted section now name the user guide. Two things filed under the program were
  physics and stayed: the pooled-covariance identity that builds the whole-window row is
  now a subsection of the diagnostic output, and the coherent source is a top-level
  section. Number conservation over the restructure: 103 distinct numeric values in the
  old manual, none lost.

- 2026-08-30 Changed: The README is a front door. It was 1467 lines and 26 sections, and
  it is now 60: what Lucifer is, the three tracking methods, how to build and run, where
  it stands, and a table of links. Its content moved rather than being rewritten, into
  four new documents under `lucifer/doc/`. `user-guide.md` holds building, the source
  map, the architecture, the FEL element's lattice attributes, the program structure and
  the output inventory. `reading-output.md` holds the statistics file, the particle dumps
  and the Python side for a reader who never runs the program. `input-reference.md` is
  the normative parameter reference. `validation.md` holds the keystone rule, the tier
  table, the check sections and the measured record subsystem by subsystem. The
  saturation demo's measured story moved to the demo's own directory, where its inputs
  already lived. Every scientific-notation value in the old README was checked into its
  new home: none was lost.

- 2026-08-30 Changed: The tier digits have one home, `doc/validation.md`. The README now
  states the shape of the validation and links, rather than carrying a second copy of the
  numbers that could drift from the first.

- 2026-08-30 Changed: `program/lucifer.f90`'s input documentation went from 234 lines to
  36: a group summary and a pointer to `doc/input-reference.md`, which is now normative.
  The refusal text that named the honored-fields table retargets in the same commit,
  since a refusal must be recognizable from its message and this one told the user where to look.

- 2026-08-30 Added: The format specifications join the tree as
  `lucifer/doc/BMAD-STATS-SPEC.md` and `lucifer/doc/BMAD-STATS-EXT-FEL.md`. The tree
  already shipped files claiming `@file_format = 'bmad-stats'`, a reader refusing other
  versions, and a validator citing rule numbers, so the definition of those bytes
  belongs beside them rather than in a working document a repository reader cannot see.
  The manual's stats section names the specification as the normative home, and the
  validator's docstring cites it by its path in the tree.

- 2026-08-30 Changed: Committed prose cites committed artifacts only. Every
  `FINDINGS.md n.m` and design-brief reference is gone from the manual, the README, the
  code comments and the scripts, with the load-bearing lesson inlined where it earned its
  line and the pointer dropped otherwise. `changelog.md` is append-only history, so past
  entries stand as written. The manual's citation-convention paragraph, which announced
  that those references pointed outside the repository, goes with them.

- 2026-08-30 Changed: Committed prose no longer shouts. Multi-word capital phrases used
  for emphasis are a house-style violation and are swept from the manual, the README, the
  scripts' docstrings and check names, and the `@description`, `@long_name` and
  `@units_note` strings the writers emit into every statistics file. What stays capital is
  what the style guide sanctions: `out_io` message text, the refusal strings the checks
  match against it, machine-parsed banners such as the wavefront suite's
  `LARGEST RELATIVE DIFFERENCE`, acronyms and code identifiers.

- 2026-08-29 Changed: The statistics file is `bmad-stats` 1.0 with the `fel`
  extension, the PLANNED VERSION RESET from the development format `lucifer-stats` 2.x.
  The layout was never FEL-specific, so the general contract now carries a general
  name: the root states `@file_format`, `@file_format_version`, `@writer`,
  `@extensions` and `@kinds` (an array), and the reader refuses an unknown version by
  name. The acceptance test is the format's own conformance checker,
  `tests/scripts/validate_bmad_stats.py`, program-blind, numpy and h5py only, run by
  the harness on the diagnostics file and required to report zero failures. The join
  key and the mask are machine-readable: `coords/ix_ele` carries `@indexes = 'ele'`
  and `coords/at_element_end` carries `@selects = 'element_end'`, the last two
  relationships a reader had to learn from prose.

- 2026-08-29 Changed: `params/` is the INPUT TREE: one subgroup per input structure
  the program honors (`global`, `beam_init` as the quiet start's honored set,
  `beam_param`, `imp`, `wavefront_init`, `wake`, `sc`, and `bmad_com` and
  `space_charge_com` whole), each carrying `@struct` naming the Fortran type, every
  honored component resolved after defaults as a true HDF5 scalar. Two runs diff by
  their inputs and no reader needs a defaults table. `out_root` alone is left out: it
  is the run's own name, and as data it made two otherwise identical runs compare
  different, which the thread-identity check caught immediately. What the run PRODUCED
  is the new `run/` group: `p0c`, the species, `slice_spacing` and the axis lengths.
  Everything that is a pure function of other datasets declares `@derived_from`: the
  twiss and modes groups, `current`, `energy`, `sigma_energy`, the field emittances,
  and `beam/bunch`, which is pooled from the slice moments.

- 2026-08-29 Added: The per-slice envelope extremes, at every record. `rel_max` and
  `rel_min` are bunch_params_struct's per-coordinate extremes RELATIVE TO THE CENTROID
  over the new seven-entry `bmad_t` axis (the six phase-space names plus t), order
  statistics no moment can reconstruct, accumulated in the per-record sweep that
  already visits every particle. NaN for an empty slice. `plot_fel`'s beam-size panel
  draws the envelope band, centroid plus rel, which is what they are for. Checked
  against numpy extremes over a `dump_beam_at` file's particles at the same plane:
  the position entries match EXACTLY (same particles, and IEEE subtraction of the
  stored centroid is deterministic), the momentum entries at 2.9e-16 across the dump's
  unit round trip, and the t entry is the z entry through -dz/(beta0 c) exactly.

- 2026-08-28 Fixed: `coords/t_slice` carried a spurious factor of `beta0`, and the
  statistics file is at `@file_format_version` 2.3. The slice grid is uniform in TIME,
  which is provable from the tracker's own code: the migration invariant carries each
  PARTICLE's beta, and with `z = -beta*c*(t - t_ref)` that beta cancels at a grid point,
  leaving an arrival-time separation of `slice_spacing/c` with no beta in it. So
  `t_slice` is `-ct_slice/c` exactly, where it read `-s_slice/(beta0*c)` and was wrong by
  3.9e-9. Nothing read `coords/` back, which is why nothing caught it, and no tier digit
  can move. FINDINGS 7.32.

- 2026-08-28 Changed: The slice axis is the slice NUMBER, with three positions as
  variables on it, which is the treatment `coords/record` already had. `ct_slice` is the
  light-travel distance ahead of the reference and `t_slice` the arrival time, both exact
  and free of beta. `z_slice` is Bmad's z AT THE REFERENCE beta and says so in its
  description, because a particle's own offset uses its own beta, which is why the
  slice-to-bunch concatenation stores every entry beta rather than one number.
  `slice_spacing` is documented as a light-travel distance: one slice is exactly
  `window_sample` wavelengths of slippage, which is what makes the field record's rotation
  an integer index shift with no interpolation, and uniformity in `ct` is the reason it
  works.

- 2026-08-28 Changed: `coords/z` is `coords/s`. It holds path length along the lattice,
  which Bmad calls s, while `coords/s_element_end` held the same quantity under an `s`
  name with a `@long_name` that said "z at element end", and `lattice/s_start` and `s_end`
  said s all along. So one quantity had two names, and z was already taken twice over as
  the fifth entry of `coords/bmad` and the third of `coords/plane`. z now means only the
  phase-space coordinate and the longitudinal twiss plane.

- 2026-08-28 Changed: A zero-length element whose wake CANNOT act is skipped like any
  other zero-length element. Bmad's `scale_with_length` defaults true and `wake_mod`
  scales the kick by the element length, so at zero length it is identically zero. This is
  not refused, since a wake assigned over an element range or a class lands on
  zero-length members as a matter of course. But honoring one cost the serial interlude
  path, a record, an element end and a repeated `coords/s`, all for a kick of zero.
  Measured: such a run is now dataset-identical to one whose zero-length pipes carry no
  wake at all. A long-range wake has no `scale_with_length` and would act, so it keeps
  the element.

- 2026-08-28 Added: The harness tracks a lattice with zero-length wake pipes in BOTH
  polarities. The check that every repeated `coords/s` straddles an element boundary had
  only ever seen zero repeats, so it was untested. It now meets a real duplicate, at the
  plane where a kicking zero-length pipe follows an undulator.

- 2026-08-28 Fixed: Attribute shape expresses arity. `@unit_power` and a harmonic group's
  `@harmonic` are true HDF5 scalars, where the high-level Fortran layer had made them
  shape-(1,) arrays that a reader had to unwrap. `@components` and `@derived_from` are
  arrays of length one or more, so a one-component file parses exactly like a
  two-component one. `units_note` also states that a coordinate variable may repeat, and
  names the case.

- 2026-08-27 Fixed: `stats.h5`'s `meta/` group was three ways wrong at once, and is
  rebuilt at `@file_format_version` 2.2 (FINDINGS 7.31). **The 64 kB attribute cap.**
  HDF5 caps a single attribute at 64 kB (measured: the largest that writes is 65495
  bytes, where a scalar string dataset took 3 MB), the echoed namelist is already 12 kB
  and a real lattice text 37 kB, and a failure only warned. A lattice under twice a real
  one's size would have written a file whose provenance was SILENTLY absent. Every text
  in `meta/` is a scalar string dataset now. The cost is that `meta/` no longer sits
  outside dataset-level identity comparisons for free, since `input_echo` carries
  `out_root`, so the harness excludes it by name with the reason on the line.
  **Completeness.** `file_text` reads ONE file, so `lattice_text` never was the
  reproducibility record it claimed to be: every wrapper lattice in the tree recorded a
  call statement while the lattice it called was absent, 91 bytes of the diagnostics
  wrapper and 485 of the tier wrapper against the 2569-byte lattice. It is now
  `lattice_source`, described as the top-level file only, beside a new
  `n_lattice_files` from the parser's own tally, and the claim is withdrawn from the
  manual and the README: reproduction rests on the `lattice/` table and the input echo.
  Serializing the lattice was examined and rejected, `write_bmad_lattice_file` inlining a
  `grid_field` as ASCII under `one_file$` and writing sibling binary files otherwise, so
  no `output_form` is both complete and bounded. **Privacy.** A stats file is meant to
  travel, so `user` and `cwd` leave the default file behind the new
  `global%record_environment`, and `lattice_file` records a base name. Genesis records
  user and cwd unconditionally. Parity is not a reason to leak.

- 2026-08-27 Changed: The twiss planes and the normal modes are SEPARATE axes.
  `coords/plane` keeps the projected x, y and z, the new `coords/mode` takes a, b and c,
  and `beam/bunch/` and `beam/slice_twiss/` each hold nine datasets in `twiss/` and nine
  in `modes/`. One axis carrying both was one axis carrying two decompositions of one
  beam: `beta` read 16.65, 8.99, 2.5e-6 beside 11.67, 6.47, 2.5e-6, a mean over the axis
  was meaningless, and a reader plotting all planes got six curves where it wanted three.
  `coords/mode` also states that its labels are eigenvector-identified rather than
  magnitude-sorted, which is why the harness compares mode emittances as a set.

- 2026-08-27 Added: The per-entry units of a centroid and a sigma live on their axis.
  `coords/bmad_unit` and `coords/wavefront_unit` are variables on those axes, and a
  dataset over one of them carries `@unit_of_axis` and `@unit_power` (1 for a centroid, 2
  for a second moment) beside the human `@unit` string, which is a comma list nothing can
  parse. Also: a root `@kinds` enumerating the group vocabulary, `@dtype_hint = bool` on
  every int8 flag so the boolean convention is something a dataset says rather than a
  rule to pattern-match, and `@plot_against` on `coords/record` and
  `coords/element_end` naming `z` and `s_element_end`, which is information only the
  writer has.

- 2026-08-27 Added: Four provenance checks in `check_diagnostics.py`. No attribute
  anywhere near the 64 kB cap. Every text in `meta/` a dataset. `n_lattice_files`
  reporting the wrapper lattice's second file, which is the case that was broken.
  And no machine-local value in a default run, with `global%record_environment`
  restoring them.

- 2026-08-26 Changed: The statistics file's axis vocabulary is complete, at
  `@file_format_version` 2.1, so a reader guesses nothing. EVERY NAME IN `@axes` NOW
  RESOLVES TO A `coords/` DATASET, the trailing label axes included: `bmad` and
  `bmad_col` for the six phase-space coordinates, `wavefront` and `wavefront_col` for the
  four field moments, `plane` for the six twiss planes. The two sides of a square matrix
  are named apart on purpose, so selecting the (x, pz) entry needs no rule rather than a
  dedupe. THE RECORD NUMBER IS THE AXIS, `coords/record`, with `z` demoted to a variable
  on it: `z` repeats wherever two records land on one plane, and a selection on a
  repeating index answers silently wrong. The element-end axis gains its own coordinates,
  `coords/element_end` and `coords/s_element_end`, and the lattice table's axis is now
  `ele` with `coords/ele` beside it, which ends the collision where `ix_ele` named two
  axes of different length. `coords/ix_ele` stays the per-record join key. Every dataset
  also carries `@long_name`, since `@description` is a sentence and an axis label wants
  three words, and every group carries `@kind` and `@description`. `params/` holds true
  HDF5 scalars rather than shape-(1,) arrays. `t_slice` carries `@head_direction` too,
  and a `sigma` matrix carries `@unit_of_axis` and `@unit_power` beside the human unit
  string, which no reader could parse. The version is a development marker and is refused
  outright, with no compatibility machinery for older files. It resets to 1.0 at the first
  external release.

- 2026-08-26 Changed: The six twiss planes of `beam/slice_twiss/` and `beam/bunch/` are
  an AXIS rather than six groups. One `twiss/` subgroup per set holds nine datasets over
  `coords/plane`, so 108 datasets became 18 and no group is named after a coordinate. A
  group named `z` beside a `z` coordinate is something xarray and netCDF both refuse. The
  names stay `bunch_params_struct`'s, as labels now, so the mapping to Bmad's struct is
  exact and machine-readable. They sit in a subgroup because one of them is `sigma`, and
  the covariance matrix beside them is `sigma` too.

- 2026-08-26 Added: `field/@components` and `field/@harmonics` name what the children of
  `field/` are, and `total/@derived_from` names what it sums, so a reader adding up the
  children cannot take the always-written derived sibling for a component and
  double-count. `@components` was specified when 2.0 landed and never written.

- 2026-08-26 Added: Two structural checks in `check_diagnostics.py`. The acceptance test
  is a GENERIC LOAD: label every dimension of every dataset from `@axes` alone, failing
  if a name does not resolve to a coordinate, if a dimension has no name, if a length
  disagrees with its coordinate, or if one dataset names an axis twice. Nothing in it
  knows a dataset name, which is the property being checked. The second holds the record
  axis: `z` non-decreasing, and every repeated `z` straddling an element boundary, since
  a repeat inside one element would be a defect in the walk that demoting `z` to a
  variable would otherwise hide. The same check runs on the three-element line in
  `check_program.py`, where every comb comparison matches rows by `z`.

- 2026-08-25 Changed: The statistics file `<out_root>.stats.h5` describes itself, at
  `@file_format_version` 2.0. Every dataset carries `@unit`, `@description` and `@axes`,
  the last naming the `coords/` datasets its dimensions run over, so a reader needs no
  table of names: `lucifer/tests/scripts/read_stats.py` is the one reader the tree uses
  and it hard-codes nothing. Units stay DOCUMENTATION, never a factor to apply. Five
  groups replace the old flat layout: `coords/` holds every axis once (including the
  SLICE axis, which the file never had, with `@head_direction` publishing which end of
  the index is the window head), `params/` holds every scalar as data so nothing is
  scraped out of the echoed namelist, `beam/slice/` the per-record sufficient statistics
  with `sigma` at its natural `(nz, ns, 6, 6)` rank, `beam/slice_twiss/` and
  `beam/bunch/` the evaluated Bmad bunch_params on the element-end grid, and `field/`
  one group per component and per harmonic, each with its own always-written `total/`.
  No dataset's meaning now depends on what else the file holds: `field/power` used to
  become a sum when a second polarization was live. Not-computed is NaN rather than a
  zero that reads as an answer. `element_end/` is gone: an element end is always a
  record, and `coords/at_element_end` marks it, which removes a duplicated copy of every
  element-end quantity. `charge_tot`, `n_particle_tot` and `beam/s` are gone too, being
  other datasets under second names. The provenance group is `meta/`, lower case with the
  rest.

- 2026-08-25 Added: The statistics file carries a `lattice/` table, one row per tracked
  element indexed by `coords/ix_ele`, so a layout plot needs nothing but the file: name,
  key, s_start, s_end, l, ds_step, is_fel, fel_tracking, b_max, aw as the physics used
  it, l_period, ku, helical, k1, tilt, z_offset. Genesis writes its lattice as per-step
  arrays; a table joined through the element index says the same thing without a second
  copy of the record axis. `beam/slice/` also gains `current`, `energy` in eV and
  `sigma_energy`, which every consumer used to re-derive.

- 2026-08-25 Changed: An element end is always a stats record, whatever
  `global%comb_ds_save` says. Bmad's comb semantics drop the comb entirely at a negative
  value; here that leaves the element ends, because the file now carries one record axis
  with a mask rather than a second axis. A `comb < 0` run therefore writes a file whose
  every record is an element end, at the positions an every-record run puts them.

- 2026-08-25 Fixed: The first stats record of a run said it sat in an uninitialized
  element rather than at the entry face. The debug build zeroed the stack, which is the
  right answer by accident, so only the production build showed it, and only the new join
  check between `coords/ix_ele` and the `lattice/` table could see it at all.

- 2026-08-25 Fixed: `hdf5_write_dataset_int_rank0` and `hdf5_write_dataset_real_rank0`
  wrote an uninitialized local to the file and then copied it back over the CALLER's
  variable, a reader's body in a writer. Neither had a caller in the tree until the stats
  file grew a group of scalars, and the symptom was memory corruption in the caller
  rather than a wrong number in the file. Both now write the value they were given, and
  `intent(in)` makes the direction structural.

- 2026-08-25 Changed: Lucifer reads and writes openPMD and nothing else. A particle dump
  is `<out_root>-final.beam.h5` and a field dump `<out_root>-final.wf.h5`, and the format
  knobs `beam_formats` and `wavefront_formats` are gone with no alias. A file that is not
  openPMD is refused on import, with the conversion command in the message.
  `lucifer/tests/scripts/convert_genesis.py` converts particles and fields between the
  Genesis format and openPMD in either direction, and the validation harness converts the
  Genesis reference dumps once per chain at its boundary, so both codes still start from
  the same state (measured: the field bit-identical, the particles at 3e-16 steady state
  and 8e-15 over 32 slices). The slice partition is now the standard's own
  `particlePatches`, one patch per slice with an empty slice as a patch of no particles,
  so the seven `fel*` root attributes are gone: the patch count is the window, `one4one`
  is what the weights say, and the wavelength, the slice spacing and the beamlet size come
  from the deck. Reading a dump with no `lambda0` is refused rather than defaulted, since a
  wrong wavelength rescales every phase in the run.

- 2026-08-25 Fixed: A particle dump now carries the whole ponderomotive phase. The chart
  splits it into a per-beam reference phase and a per-particle lag, and no dump format has
  anywhere to put the reference, so every reader restarts it at zero. The writer folds the
  reference into the lag it writes, which makes the file's time coordinate
  `-theta/(ks c)` and a restart exact. Without the fold a mid-line restart placed the beam
  at a different phase against the same dumped field, which the windowed-composition check
  measured at 2.1e-2 and now measures at 3.7e-14.

- 2026-08-25 Fixed: `hdf5_write_beam` writes its `particlePatches` records as datasets
  even where every patch holds the same particle count. Bmad's dataset writer collapses an
  all-equal array to a constant-value group, which is a legal openPMD record component and
  not a legal patch list, so a single-patch file, or any file whose patches held equal
  counts, came back with no partition at all. `hdf5_read_beam` names that form rather than
  reading zero particles from it.

- 2026-08-25 Changed: The eleven benchmark tiers are re-recorded. Both codes still start
  from the same state, but the tracker now reads it through a conversion, which costs a
  multiply and a divide: the initial field is bit-identical and the initial particles move
  in their last digits (3e-16 steady state, 8e-15 over 32 slices). Each tier's digits then
  move by that times its own sensitivity, measured here by perturbing the input by 1e-15
  and rerunning: 2.6e4 for the averaged FEL core, and saturating at ~1e-5 for the
  unaveraged mode, whose final-field phase is chaotic and whose recorded digits are
  therefore only reproducible on a bit-identical input. tier1 1.825901e-06 ->
  1.825899e-06, tier1_unavg 6.933979e-02 -> 6.934613e-02, tier2_genesis 1.771895e-05 ->
  1.771890e-05, tier2_bmad 5.001254e-02, td1 8.467690e-07, td2_genesis 2.398226e-06,
  td2_bmad 4.127587e-02, tdsase 2.292906e-06 -> 2.292496e-06, tdsc 2.440477e-04, tdwk
  8.708129e-07, weight_split 3.508953e-13 -> 3.532394e-13, all on the debug tree. The
  production tree lands at the same digits except where the build reaches: tier1_unavg
  6.934017e-02, td2_genesis 2.397983e-06, tdwk 8.708128e-07, weight_split 3.536001e-13.
  Every tolerance is unchanged.

- 2026-08-25 Changed: A Lucifer run's particle dumps are openPMD by default
  (`<out_root>-final.beam.h5`), so a beam with per-particle weights now survives a dump
  and a restart. Genesis `.par.h5` holds one current per slice, and writing a
  nonuniform-weight beam to it is refused rather than silently returning a
  uniform beam on read. Both dump kinds now take a LIST of formats instead of an enum,
  `beam_formats` and `wavefront_formats`, so a third code costs a token rather than
  another combination. `wavefront_format = 'both'` is retired with no alias, and
  `beam_file` accepts either format by signature.

- 2026-08-24 Changed: Running Lucifer's validation harness now needs an
  openPMD-beamphysics checkout carrying `beamphysics/wavefront/openpmd.py`
  (`../openPMD-beamphysics` by default), and the harmonics section refuses
  without it. The Python side of the wavefront round trip was a patch carried in
  `lucifer/openpmd/`. It has landed upstream, so the patch is removed and the checks
  exercise `Wavefront.from_openpmd` and `write_openpmd` directly. The Fortran reader
  is now also checked against Python-written files.

- 2026-08-23 Changed: Lucifer's terminal output is formatted for humans: a framed configuration
  header, a progress table with SI-prefixed values and fixed numeric columns, and a completion
  block listing the files written with their sizes. The progress row now carries power, energy
  and bunching in every comb mode. Do not parse stdout. Machine-readable output goes to files,
  with the ALL-CAPS refusal texts the one documented exception.

- 2026-08-23 Added: Lucifer writes the distribution-import moments and per-slice current profile
  to `<out_root>.import.txt`, and one row per slice-migration event to
  `<out_root>.migration.txt`. Both data streams previously went to stdout.

- 2026-08-23 Added: Radiation power, energy, on-axis intensity and bunching per slice are now in
  Lucifer's `element_end/` stats group, so a run with `comb_ds_save < 0` (no per-record rows
  kept) still has these quantities at element ends. Bmad's comb semantics are unchanged.

- 2026-08-23 Changed: A 131-slice Lucifer run finishes in 126.9 s instead of 137.7 s, with
  utilization up from 931% to 1048%. The element-end whole-window bunch statistics are now
  assembled from the per-slice moments by the pooled-covariance identity instead of
  concatenating every particle in the time window into one bunch and running the full 6D
  moments and Twiss on it. That removes 110 million single-threaded particle visits, 16.5% of
  the run's wall clock. Agreement with the particle sum measured at 4.0e-12 / 5.0e-11 on two
  physical configurations.

- 2026-08-22 Fixed: `util/searchf.py` (getf/listf/create_searchf_namelist): the
  interface-end regex was `end\s+ interface` (a stray space requiring two whitespace
  characters), so a bare `interface` block made the scanner swallow the rest of the file.
  Every routine after such a block was missing from `searchf.namelist` and invisible to
  getf/listf wherever an index file existed. Directories with tracked index files may want
  to regenerate them with `util/create_searchf_namelist`.

- 2026-08-22 Changed: Lucifer messages appear as routine-tagged, severity-tagged blocks, and the
  old `fel_track_test:` stdout prefixes are gone. The messages go through `out_io` (Bmad's
  standard message system). Data lines that scripts parse (import moments/currents, progress,
  file listings) stay bare and full-precision.

- 2026-08-22 Added: Lucifer, an FEL tracker validated against Genesis 1.3 Version 4, as a
  top-level program directory (`lucifer/`, executable `lucifer`, library `liblucifer`).
  Time-dependent SASE and seeded tracking with slippage on Bmad lattices (wiggler/undulator
  elements), wakes and space charge, distribution import, harmonic fields, an unaveraged
  verification mode, openPMD wavefront I/O, and a validation harness
  (`lucifer/tests/run_fel_benchmark.sh`). Physics manual in `lucifer/doc/fel-physics.tex`.
