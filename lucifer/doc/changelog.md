---
title: Lucifer changelog
short_title: Changelog
---

# Lucifer changelog

Development history of the FEL tracker on the `lucifer-dev` branch, newest first.
This is the branch's own record. Bmad's `changelog.md` carries what a merge changes,
and it is written at the merge.

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
