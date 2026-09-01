---
title: Lucifer and Genesis 1.3 Version 4 (Genesis4)
short_title: Genesis4
---

Lucifer's averaged physics was transcribed from Genesis4, and
the benchmark tiers compare against it. This page is for a user who knows one code and
wants the other: what physics the two share, what each has that the other does not, how
a setting translates, and how a file moves between them.

**Checked against Genesis4 release v4.6.14.** This page is the only one that tracks
Genesis4, so it is the only one to revisit when Genesis4 changes. Everything else in
this documentation describes Lucifer on its own terms.

## Shared physics

These are the same model in both codes, transcribed routine by routine, and the
comparison tiers hold them to the level [](validation.md) records.

| Method | Shared |
|---|---|
| Field solver | The SVEA paraxial wave equation on a transverse grid, split-step with an FFT diffraction kernel and a $+2\cdot$source deposit |
| Averaged FEL step | The wiggle-averaged (KMR) model: transverse push with natural focusing, RK4 advance of the ponderomotive phase and energy, per-particle source deposition |
| Slippage | An exact integer rotation of the field record, one wavelength per undulator period, on a slice grid uniform in arrival time |
| Shot noise | Fawley's loading algorithm, here in a weighted generalization |
| Wakes | The Bane and Stupakov numerical resistive-wall impedance, plus geometric and roughness kernels, with Genesis4's exact numerics |
| Space charge | Short-range harmonic and long-range solvers at Genesis4's granularity |
| Distribution import | The `importdistribution` resampling method, transcribed from `SDDSBeam.cpp` |
| Spontaneous emission | The same undulator quantum-diffusion form, though Lucifer gates it on Bmad's own switches (see below) |

## What each code has that the other does not

Lucifer has, and Genesis4 does not:

- **Per-particle weights** throughout, including weighted shot noise and weighted diagnostics. No Genesis4 dump format carries a weight.
- **The unaveraged mode**, integrating the real undulator field with no period averaging and no resonance approximation.
- **Bmad element wakes** applied across the whole time window, and every non-FEL element tracked by Bmad itself: real quadrupoles, chicanes, patches, apertures and collimators.
- **Slice migration**, moving particles between slices when their ponderomotive phase leaves the slice window.
- **A lattice**, in Bmad's sense. FEL parameters are element attributes, so one lattice serves tracking, optics and layout.
- **A self-describing statistics file**, described by [](BMAD-STATS-SPEC.md).

Genesis4 has, and Lucifer does not:

- **MPI**, where Lucifer is shared-memory OpenMP over slices.
- **One-to-one tracking** (`one4one`) and the particle sorting it implies.
- **Undulator field errors**, and taper and error models Lucifer has not transcribed.
- **The ADI field solver**. Lucifer transcribes the FFT solver only, so a comparison run sets `fft_fieldsolver = true`.
- **Its own semianalytic and scan machinery**, none of which is in scope here.

## Translating settings

Genesis4 is driven by namelists over its own lattice format. Lucifer is driven by three
namelists over a Bmad lattice, with the FEL parameters living on the elements. The
parameter reference is [](input-reference.md).

| Genesis4 | Lucifer | Note |
|---|---|---|
| `&setup lambda0` | `wavefront_init%lambda0` | Required when starting from a beam dump, since the dump carries the slice partition and not the wavelength it was sliced on |
| `&setup delz` | the element's `ds_step` | Bmad's own step attribute. The step count `round(l/ds_step)` matches Genesis4's unroll |
| `&setup rootname` | `global%out_root` | |
| `&time sample` | `wavefront_init%window_sample` | Slice spacing in wavelengths, an integer in both |
| `&time slen` | the number of slices in the starting state | Lucifer takes the window from the beam it is given rather than a length |
| `&beam npart` | `beam_init%n_particle` | Per slice in both. Lucifer requires a positive multiple of `beamlet_size` |
| `&beam nbins` | `beamlet_size` | Beamlet size of the quiet start. No dump format carries it, so it is named on both sides |
| `&importdistribution npart` | `resample%n_particle_per_slice` | Per slice in both, and Lucifer's name says so |
| `&importdistribution nbins` | `resample%beamlet_size` | Beamlet size of the resample |
| `&importdistribution nslice` | `resample%n_slice` | Slice count. Zero derives it from the bunch length |
| `&track zstop` | `global%track_end` | Lucifer's bound is an element locator, not a distance. `global%track_start` has no Genesis4 equivalent |
| `&track fft_fieldsolver` | (always) | Lucifer transcribes the FFT solver only |
| `&importdistribution slicewidth` | `resample%slice_width` | Same meaning, sampling window over bunch length |
| `&importdistribution` (match, center) | not ported | A Bmad lattice matches the beam, so these are the lattice's job |
| `&sponrad` | `bmad_com%radiation_damping_on`, `bmad_com%radiation_fluctuations_on` | Lucifer exposes Bmad's own switches directly, as Tao does, rather than an FEL-local flag. Both default off, matching Genesis4's `&sponrad` default |
| `&write beam`, `&write field` | `beam_file`, `field_file` on the input side; openPMD dumps on the output side | See file exchange below |
| `&importbeam`, `&importfield` | `beam_file`, `field_file` | Lucifer reads openPMD only, so convert first |
| (no equivalent) | `global%migrate` | Slice migration. Off by default, since the comparison tiers run against a code that never migrates |
| (no equivalent) | `tracking_method = fel_averaged` or `fel_unaveraged` | Selects the method per element, and `global%transport_model` selects the averaged method's maps |

## Conventions that differ

- **Energy.** Lucifer reports energy as $E = \gamma m_e c^2$ in eV, which is Bmad's convention. Genesis4's `energy` output is $\gamma$. The comparison scripts convert.
- **`s` and `z`.** In Bmad, `s` is position along the lattice and `z` is the longitudinal phase-space coordinate. The statistics file uses both words in exactly that sense. Genesis4 uses `s` for the position within the time window, which is closest to Lucifer's `z_slice`.
- **Units in files.** openPMD's `unitSI` is a factor a reader must apply. The `bmad-stats` file's `@unit` is documentation, and its values are already SI and eV. The two conventions never appear in one file.
- **Field amplitude.** Genesis4 stores an internal field unit; Lucifer stores V/m. The identity between them, taken from `writeFieldHDF5.cpp`, is the one every physical-unit formula here was derived through.
- **Particle ordering.** Neither code sorts without `one4one`, so final dumps compare particle by particle rather than only statistically.

## Moving files between the codes

`lucifer/tests/scripts/convert_genesis.py` holds all of the format knowledge, in both
directions:

```
python3 convert_genesis.py <in> <out>
```

Fields convert both ways and Genesis4 particles convert inward through
openPMD-beamphysics, which owns those conversions. The fourth direction, openPMD
particles out to a Genesis4 `.par.h5`, is implemented here because no upstream
implementation exists. That is what Genesis4's `&importbeam` reads, so it is how
Genesis4 restarts from a state Lucifer produced.

A file that is not openPMD is refused by name, with the conversion command in the
message.

## Facts about Genesis4 this work pinned down

These are properties of Genesis4 that the transcription had to establish by reading the
source, and they are recorded because they are not stated in its documentation. Each is
true of the release named at the top of this page.

- Outside undulators Genesis4 does not subdivide into `delz` steps: each interlude element
  is one integration step of the element's full length (`Lattice::unrollLattice` pushes
  one entry per non-undulator layout segment). The step count over the benchmark is
  12*89 undulator steps + 36 interlude steps = 1104.
- In steady state `slippage` and `phaseshift` are no-ops: `Control::applySlippage`
  returns immediately when not time dependent, and the phaseshift array is all zero
  without phase shifter elements.
- Time dependence follows from `&time` or from importing a multi-slice dump
  (`readBeamHDF5.cpp:68-77` reconstructs the window from `slicespacing` and `slicecount`,
  time on by default), with `nslice = round(slen/(sample*lambda0))` (`GenTime.cpp:70`).
- The end-of-lattice autophasing is unguarded: the last step always gets
  `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths of slippage, `+1` even with no trailing
  drift (`Lattice.cpp:191-193`).
- The field record's rotation never appears in Genesis4's outputs: `writeFieldHDF5` and
  `DiagField::calc` both unrotate on the fly, so `.out.h5` per-slice arrays and `.fld.h5`
  dumps are in time-window order, aligned with beam-slice indexing.
- A helical undulator defaults to `kx = ky = 0.5` (`LatticeParser.cpp:329`), scaled by
  `ku^2` in the unroll.
- The `&importbeam` / `&importfield` namelists take full filenames and make the
  shared-start methodology possible. `&write` before `&track` produces the dumps.
- Genesis4's internal field unit: `dfl [sqrt(W)] = u * dgrid * eev / (ks * sqrt(vacimp))`
  (`writeFieldHDF5.cpp:70`), with `vacimp = 376.73` truncated and `eev = 510998.95069`
  (identical to Bmad's `m_electron`). Equivalently `u = E*ks/(sqrt(2)*eev)` with E in
  V/m, the relation used to derive this tracker's physical-unit formulas.
- Releases up to v4.6.14 use `eev = 510999.06` for the electron rest energy, which is
  2.14e-7 above the CODATA 2022 value. Genesis master carries 510998.95069, which is also
  Bmad's `m_electron`, so a comparison against master has no electron-mass discrepancy at
  all while one against a release has this one. The field unit is linear in `eev`, so the
  difference reaches the field-normalization integrals as 4.28e-7 and loosens the
  transcription tiers by several percent of their own value without failing any of them.
  The levels recorded in [](validation.md) belong to a build that carries the CODATA value.
- The quad transport is chromatic through per-particle `gammaz` with `foc^2 =
  k1*gamma0/gammaz`. Bmad's equivalent scaling is `k1*p0/p`. The two differ by ~4e-9
  relative at this energy (1 - beta0), far below the path-length-term difference.

Genesis4's truncated impedance constants set the floor of every transcription
comparison, which is why nine of the eleven tiers agree at that floor rather than at
round-off. [](validation.md) says which tiers those are and what the remaining two
price.
