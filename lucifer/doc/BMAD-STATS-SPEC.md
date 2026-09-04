# bmad-stats: a self-describing HDF5 layout for accelerator statistics

Version 1.0, in force for this writer, with open questions pending. This document is normative for any writer claiming `@file_format = 'bmad-stats'`. It is versioned by its own number, independently of any program. Lucifer is the reference writer, as `bmad-stats 1.0` plus the `fel` extension (BMAD-STATS-EXT-FEL.md).

Rules are numbered R1, R2, ... so markup can point at them, and the numbers are append-only: a rule keeps its number forever, so late additions sit out of sequence (R43 and R44 live in section 12). MUST, SHOULD and MAY are used as in RFC 2119. Open questions are collected at the end as O1, O2, ...

## 1. Purpose and scope

This standard covers the statistics and diagnostics a tracking program produces about a run: moments and Twiss parameters along a machine, per-element tables, resolved input parameters, derived run scalars, and provenance. Candidate writers beyond Lucifer are Tao's bunch tracking, beam_track, long_term_tracking, and anything else that today invents an ASCII column format per program.

It deliberately does *not* cover particle dumps or field meshes. Those are openPMD's job, and the two standards carry opposite unit conventions (R21), so one file never carries both.

The design goal in one sentence: a reader guesses nothing. Every dimension, unit, meaning and relationship is stated in the file, and the machine-checkable part of that claim is enforced by a generic validator that knows no program (section 14).

## 2. Definitions

The datasets of a file share dimensions: many run over the same records, the same slices, the same six phase-space entries. Everything in this standard hangs off naming those shared dimensions and stating, once each, what lies along them.

- **axis**: a named dimension that datasets share (`record`, `slice`, `bmad`). Every dimension of every dataset names its axis in `@axes` (R11).
- **coordinate**: the one dataset that defines an axis by listing what lies along it: the record numbers, the slice numbers, the six phase-space names. Axes and coordinates are one to one (R15).
- **variable on an axis**: a dataset giving one more per-entry quantity of an axis without defining it: the path length at each record, the lattice element holding each record, the unit of each phase-space entry. An axis has one coordinate and any number of variables.
- **`coords/`**: the group holding every coordinate and every variable, once each (section 6). The name is netCDF's, not Bmad's: see the terminology note below.
- **label axis**: an axis whose coordinate is a string array (`bmad`, `plane`), used for the trailing dimensions of vectors and matrices so no reader guesses a dimension from its length.
- **join key**: an integer variable whose values index another axis (R18).
- **mask**: a boolean variable selecting another axis's entries, in order (R19).
- **data tree**: a group holding one `coords/` and the data groups whose `@axes` resolve against it (section 4).
- **kind**: the one-word class of a group, from the root's `@kinds` vocabulary.
- **extension**: a named, separately documented addition to this core (section 13).

Terminology note. "Coordinate" is used here in netCDF's coordinate-variable sense, three decades of prior art for exactly this concept and what xarray exposes as `coords`. It *never* means phase-space coordinates in this document: phase-space values live in data groups (`beam/centroid`), and the phase-space names are the labels of the `bmad` axis. O7 records the alternative of renaming the group `axes/` to avoid the domain collision entirely.

## 3. File identity

- **R1.** The root group MUST carry `@file_format = 'bmad-stats'`.
- **R2.** The root MUST carry `@file_format_version = 'MAJOR.MINOR'`. A reader SHOULD refuse a version it does not know rather than guess, and say which version it found.
- **R3.** The root MUST carry `@writer`, the program name and version (`'lucifer 1.0'`).
- **R4.** A file using extensions MUST list them in `@extensions`, a string array. A file using none MUST omit the attribute: absence means none, because a zero-length string array cannot pass safely through HDF5's Fortran layer (see R11).
- **R5.** The root MUST carry `@kinds`, a string array enumerating every `@kind` value the file uses, so a reader that sorts groups never discovers the vocabulary by inspection.
- **R6.** The root SHOULD carry `@units_note`, one human paragraph stating R21 and the file's conventions in prose.

## 4. The data tree

- **R7.** The root group is itself a data tree. Every `@axes` name in a tree resolves against that tree's own `coords/`.
- **R8.** A program tracking several branches or universes puts each additional tree under `branches/<name>/`, a complete data tree carrying `@universe` and `@branch` identity attributes. The root tree is the default branch. (O2.)

## 5. The self-description contract

- **R9.** Every dataset MUST carry four attributes: `@unit` (R21), `@long_name` (a few words, for an axis label or table heading), `@description` (one sentence, for a reader who has never seen the file), and `@axes` (R11).
- **R10.** Every group MUST carry `@kind` and `@description`.
- **R11.** `@axes` is a comma-separated string of axis names in the stored (h5py, row-major) order, or the single word `none` for a scalar. Every name in it MUST resolve to an axis (R15) of the same data tree. This is a deliberate exception to R24: HDF5's Fortran string layer overruns a zero-length value into the caller's stack, so the scalar case needs a sentinel, and a reader pays one split.
- **R12.** A dataset's rank MUST equal the number of names in its `@axes`, and each dimension's length MUST equal its axis's coordinate length.
- **R13.** A dataset MUST NOT name one axis twice. A square matrix uses a paired second axis carrying the same labels under a distinct name (`bmad`, `bmad_col`), so selecting one entry needs no dedupe rule.

## 6. Axes and coordinates: the coords/ group

`coords/` is a data tree's table of dimensions. It holds no physics results. It holds what the results are indexed by: one coordinate per axis plus the variables attached to each axis, once each, so a reader has one place to look and every name in every `@axes` has one place to resolve. From a real FEL file, the three roles marked:

```
coords/record          0 1 2 ...          axis      the record number
coords/s               [m]                variable  path length at each record (repeats, R16)
coords/ix_ele                             variable  join key, @indexes = 'ele'
coords/at_element_end                     variable  mask, @selects = 'element_end'
coords/element_end                        axis      the element ends
coords/s_element_end   [m]                variable  path length at each end
coords/slice           0 1 2 ...          axis      the slice number
coords/ct_slice        [m]                variable  light-travel distance of each slice
coords/ele             0 1 2 ...          axis      the lattice table's rows
coords/bmad            x px y py z pz     axis      labels of a phase-space vector
coords/bmad_unit       m 1 m 1 m 1        variable  unit of each bmad entry (R23)
coords/plane           x y z              axis      the projected twiss planes
```

- **R14.** Each data tree MUST hold `coords/`, kind `axis`, containing every axis and every axis-attached variable, once each.
- **R15.** An axis is a `coords/` dataset whose own `@axes` names itself (`coords/record` has `@axes = 'record'`). A `coords/` dataset whose `@axes` names a different axis is a variable on that axis.
- **R16.** Where a physical coordinate can repeat, the index MUST be the axis and the positions variables on it. A repeated index answers a selection silently wrong, and which repeats are legitimate is information only the writer has.
- **R17.** An index axis SHOULD carry `@plot_against`, naming the variable a viewer puts on the abscissa in place of the index.
- **R18.** A join key MUST declare `@indexes = '<axis>'` on the variable that carries it (`coords/ix_ele` indexes the `ele` axis). A join is then a gather, machine-readably.
- **R19.** A mask MUST declare `@selects = '<axis>'` (`coords/at_element_end` selects the `element_end` axis, in order).
- **R20.** A label axis's coordinate is a fixed-length string array. When its entries carry different units, `coords/<axis>_unit` holds them, one per entry (R23).

## 7. Units

- **R21.** Values are stored in fixed units: m, rad, s, eV, C, A, J, W, T, and `'1'` for a dimensionless quantity. `@unit` is documentation. A reader MUST NOT scale by it. This is the opposite of openPMD's `unitSI`, which is a factor to apply, and the two conventions MUST NOT appear in one file.
- **R22.** Counts and ratios carry `@unit = '1'`.
- **R23.** A dataset whose entries carry different units (a phase-space vector, a sigma matrix) MUST declare `@unit_of_axis` naming the label axis whose `coords/<axis>_unit` carries the per-entry units, and `@unit_power`, a true scalar integer: 1 for a first moment, 2 for a second. The human `@unit` string stays for display and is not parseable by design.

## 8. Types, shape, arity

- **R24.** Shape expresses arity. A quantity that is one value is a true HDF5 scalar (scalar dataspace). A quantity that is a list is an array, length one included, so a one-component file parses exactly like a two-component one. This applies to datasets and attributes alike, with R11's stated exception.
- **R25.** A boolean is int8 with `@unit = '1'` and `@dtype_hint = 'bool'`. HDF5 has no boolean type, and a flag stored as a float is a lie about what it is.
- **R26.** An enumeration is written as its name, never an internal integer code: element class via `key_name`, tracking methods, particle species and state likewise. A file outlives any code table.
- **R27.** Text of provenance size is a scalar string dataset, never an attribute. A single HDF5 attribute stops writing near 64 kB (measured: 65495 bytes is the largest that writes), and for provenance the failure mode is silent absence, the worst a provenance field has.

## 9. Missing data, and independence of meaning

- **R28.** Not computed is NaN, never a zero that reads as an answer. Where NaN cannot serve (integer data), a boolean validity dataset accompanies it.
- **R29.** No dataset's meaning may depend on what else the file holds. A quantity that would change meaning when a sibling appears is instead written under a derived group that names its inputs (R38).

## 10. Structs

- **R30.** A group mirroring a Fortran derived type MUST carry `@struct = '<type_name>'` and use the type's component names verbatim, so the mapping to the code is exact and machine-readable. It maps names, not completeness: a writer stores the components it computed (R43).
- **R31.** A quantity varying over a set is one dataset over a label axis, never N datasets with suffixed names.
- **R32.** Two decompositions of one thing are two axes (projected planes `x, y, z` against normal modes `a, b, c`), because a mean across them is meaningless and one axis invites it.

## 11. Core groups

- **R33.** `coords/`: section 6.
- **R34.** `params/`, kind `input`: one subgroup per input structure the program honors, each also kind `input` and carrying `@struct` (R30), holding the resolved values of every honored component, after defaults. Components the program refuses are not written: recording a knob that did nothing would be a lie about the run. A reader never needs a defaults table, and two runs diff by their inputs alone. For a Bmad program that is typically `params/bmad_com/`, `params/space_charge_com/`, `params/beam_init/`, and the program's own globals.
- **R35.** `run/`, kind `run`: the scalars a run produces, separated from `params/` because a reader must know what the user controlled apart from what came out: the reference momentum, species, axis lengths kept as bookkeeping cross-checks. Not `derived` in R38's sense, since their inputs are the particles and the lattice rather than other datasets here.
- **R36.** `lattice/`, kind `table`: when the program tracks a lattice, one row per tracked element on an `ele` axis, element 0 included. Minimum columns `name`, `key`, `s_start`, `s_end`. A writer SHOULD add what a layout plot needs. Per-element attributes MAY be one dataset per non-default attribute. It is *not* a lattice serialization and MUST NOT claim to be one. (O4.)
- **R37.** `meta/`, kind `provenance`, at the root: datasets only (R27), an ISO 8601 timestamp, and the writer's version. It MUST NOT record user identity or machine-local paths unless the run explicitly opts in, because a stats file is meant to travel, and it MUST NOT claim completeness it lacks: a top-level input text says it is the top-level text.

## 12. Derived data and sufficiency

Within this file, the primary data are what the program computed from particles: the moments (centroid, sigma matrix, time moments), the order statistics (per-coordinate minima and maxima, which no moment can reconstruct), and the weights and counts. Everything Twiss is a pure function of those: the projected planes from the sigma matrix exactly, the normal modes through an eigendecomposition. Bmad separates the two steps already (`calc_emittances_and_twiss_from_sigma_matrix`), and the reference writer feeds that routine from its own stored moments.

- **R38.** Data that is a pure function of other data in the tree declares `@derived_from`, a string array of the sibling names it derives from. A group may declare it once for everything it holds, keeping its own kind (`twiss`), and a group whose whole identity is the derivation (`field/total`) has kind `derived`, which then requires the declaration. Two readers depend on the declaration: one summing children must not double count, and one re-deriving must know what is independent information. Scalar context from `params/` and `run/` (a spacing, a reference momentum) is ambient and is not named: `@derived_from` lists the axis-carrying inputs, and the formula belongs in `@description`.
- **R39.** A parent whose children are summable components MUST declare `@components`, a string array, and MAY declare further child lists (`@harmonics`).
- **R43.** Sufficiency layers. A writer MUST store the primary statistics it computed and MAY store derived conveniences beside them, marked per R38. Derived data MAY live on a sparser axis than its inputs (element ends against every record), because the primary layer makes any point reconstructible offline. A writer that computes only moments conforms.
- **R44.** Stored derived data is the reference implementation's answer, which is the reason to store it at all: the normal-mode labels a, b, c are identified by eigenvector structure, a convention no reader should reimplement by guesswork, and independent eigensolvers agree only to their common floor, where the projected planes are exact algebra. The writer's harness SHOULD verify the round trip, re-deriving from the stored primary layer and comparing against the stored derived layer, matching mode sets rather than mode labels.

## 13. Extensions

- **R40.** An extension is named in `@extensions`, documented in its own `BMAD-STATS-EXT-<NAME>.md`, and carries its version in a root attribute the extension defines (`@fel_version`). (O3.)
- **R41.** An extension adds groups, axes, kinds, `params/` structs and documented attributes. It MUST NOT alter core semantics, and everything it adds MUST satisfy R9 through R32, so a generic core load never breaks on an extension it has never heard of. That property is the point of the whole contract.

## 14. Conformance

- **R42.** A file conforms when `validate_bmad_stats.py` reports no MUST failure. The validator is part of this standard, knows no program and no extension, and checks everything generically checkable. Some MUSTs are writer-side honesty the validator cannot see: that every join key and mask is declared (it can only verify declared ones resolve), that R26's strings really are the code's enums, that R28's NaN means not-computed, and that R34's trees are complete. Those belong in each writing program's own test harness, beside the physics checks.

## 15. Relation to openPMD

Complementary, not competing. Particles and field meshes go to openPMD, whose `unitSI` convention is load-bearing there and banned here (R21). The base-plus- extensions structure of this document is deliberately parallel to openPMD's, because this ecosystem's readers already know that model.

## 16. Worked example

```
/                          @file_format @file_format_version @writer @extensions @kinds
  meta/                    @kind=provenance   input_echo  lattice_source  timestamp ...
  coords/                  @kind=axis
    record  s(@axes=record)  ix_ele(@indexes=ele)  at_element_end(@selects=element_end)
    element_end  s_element_end  ele  bmad  bmad_col  bmad_unit  bmad_t  plane  mode
  params/                  @kind=input
    bmad_com/              @struct=bmad_common_struct     one scalar per component
    beam_init/             @struct=beam_init_struct
  run/                     @kind=run          p0c  n_record  species ...
  lattice/                 @kind=table        name  key  s_start  s_end  k1 ...
  beam/                    @kind=beam         @struct=bunch_params_struct
    centroid(record,bmad)  sigma(record,bmad,bmad_col)  rel_max  rel_min    primary
    twiss/(plane)  modes/(mode)               @derived_from=['centroid','sigma']
```

## 17. Open questions

- **O1.** `@axes` as a comma string with a `none` sentinel (R11) or as a string array, accepting the Fortran zero-length hazard on scalars. Current draft: string, for the stack-corruption reason, one split per reader.
- **O2.** Branch layout: root-as-default-tree plus `branches/` (current draft, zero churn for single-branch writers) against everything under `branches/` with a `main` link (one uniform shape, two loops for nobody).
- **O3.** Extension versioning: a root attribute per extension (draft) or versioned entries in `@extensions` (`'fel-1.0'`).
- **O4.** Whether `lattice/` is core or belongs to an accelerator extension. Draft: core, since every candidate writer tracks a lattice.
- **O5.** Whether Bmad's name arrays (`key_name`, species names) become normative vocabularies of this standard or stay writer-defined strings.
- **O6.** Whether extension-defined kinds should be namespaced (`fel:field`) or flat with first-come naming (draft: flat, the vocabulary is small).
- **O7.** The group name. `coords/` follows netCDF and xarray, the right prior art for the concept, but collides with this domain's reading of "coordinates" as phase space. `axes/` has no collision and matches the `@axes` attribute it resolves, at the cost of holding variables that are not axes. Draft: `coords/`, with the section 2 terminology note. The rename is mechanical if the collision confuses in practice.
- **O8.** Name-based `@axes` against HDF5's own Dimension Scales (H5DS), the format's native mechanism for exactly this and the base netCDF-4 is built on. Name-based resolution is one string attribute in Fortran, greppable in h5dump output, and trivial for the validator, where scales are object references. Scales would buy dimension display in HDF5 tools and, through h5netcdf, named dimensions in xarray with no custom loader. Draft: name-based is normative and a writer MAY attach dimension scales as well. Revisit if native xarray opening proves worth the machinery.
