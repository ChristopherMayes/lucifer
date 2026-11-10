---
title: Lucifer and Genesis 1.3 Version 4 (Genesis4)
short_title: Genesis4
---

Lucifer's averaged physics was transcribed from Genesis4, and the Genesis4 tiers compare against it. This page is for a user who knows one code and wants the other. It says what physics the two share, what each has that the other does not, how a setting translates, and how a file moves between them.

This page is checked against Genesis4 release v4.6.16, which is also the reference the comparison tiers run against. This page is the only one that tracks Genesis4, so it is the only one to revisit when Genesis4 changes. Everything else in this documentation describes Lucifer on its own terms.

## Shared physics

These are the same model in both codes, transcribed routine by routine, and the comparison tiers hold them to the level [](validation.md) records.

| Method | Shared |
|---|---|
| Field solver | The SVEA paraxial wave equation on a transverse grid, split-step with an FFT diffraction kernel and a $+2\cdot$source deposit |
| Averaged FEL step | The wiggle-averaged (KMR) model: transverse push with natural focusing, RK4 advance of the ponderomotive phase and energy, per-particle source deposition |
| Slippage | An exact integer rotation of the field record, one wavelength per undulator period, on a slice grid uniform in arrival time |
| Shot noise | Fawley's loading algorithm, here in a weighted generalization |
| Wakes | The Bane and Stupakov numerical resistive-wall impedance, plus geometric and roughness kernels, with Genesis4's exact numerics |
| Space charge | Short-range harmonic and long-range solvers at Genesis4's granularity |
| Distribution import | The `importdistribution` resampling method, transcribed from [`SDDSBeam.cpp`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/SDDSBeam.cpp) |
| Spontaneous emission | The same undulator quantum-diffusion form, though Lucifer turns it on with Bmad's own switches (see below) |

## What each code has that the other does not

Lucifer has, and Genesis4 does not:

- **Per-particle weights** throughout, including weighted shot noise and weighted diagnostics. No Genesis4 dump format carries a weight.
- **The unaveraged mode**, integrating the real undulator field with no period averaging and no resonance approximation.
- **Bmad element wakes** applied across the whole time window, and every non-FEL element tracked by Bmad itself: real quadrupoles, chicanes, patches, apertures and collimators.
- **Slice migration**, moving particles between slices when their ponderomotive phase leaves the slice window. With `load_mode = "keep"`, which keeps every real particle at its own weight, migration is what Genesis4's one-to-one tracking and its sorting do ([](fel-physics.md#sec-loading), [](fel-physics.md#sec-migration)).
- **A lattice**, in Bmad's sense. FEL parameters are element attributes, so one lattice serves tracking, optics and layout.
- **A self-describing statistics file**, described by [](BMAD-STATS-SPEC.md).

Genesis4 has, and Lucifer does not:

- **MPI**, where Lucifer is shared-memory OpenMP over slices.
- **Undulator field errors**, the field-error model of its `&lattice` namelist. A taper is a per-element attribute of the Bmad lattice and needs no model.
- **The ADI field solver**. Lucifer transcribes the FFT solver only, so a comparison run sets `fft_fieldsolver = true`.
- **Its own semianalytic and scan machinery**, none of which is in scope here.

## Translating settings

Genesis4 is driven by namelists over its own lattice format. Lucifer is driven by three namelists over a Bmad lattice, with the FEL parameters living on the elements. The parameter reference is [](input-reference.md).

| Genesis4 | Lucifer | Note |
|---|---|---|
| `&setup lambda0` | `field.lambda0` | Required when starting from a beam dump, since the dump carries the slice partition and not the wavelength it was sliced on |
| `&setup delz` | the element's `ds_step` | Bmad's own step attribute. The step count `round(l/ds_step)` matches Genesis4's unroll |
| `&setup rootname` | `run.out_root` | |
| `&time sample` | `window.n_wavelength` | Slice spacing in wavelengths, an integer in both |
| `&time slen` | the number of slices in the starting state | Lucifer takes the window from the beam it is given rather than a length |
| `&beam npart` | `beam.n_particle` | Per slice in both. Lucifer requires a positive multiple of `beamlet_size` |
| `&beam nbins` | `beamlet_size` | Beamlet size of the quiet start. No dump format carries it, so it is named on both sides |
| `&importdistribution npart` | `beam.resample.n_particle_per_slice` | Per slice in both, and Lucifer's name says so |
| `&importdistribution nbins` | `beam.resample.beamlet_size` | Beamlet size of the resample |
| `&importdistribution nslice` | `beam.resample.n_slice` | Slice count. Zero derives it from the bunch length |
| `&track zstop` | `run.track_end` | Lucifer's bound is an element locator, not a distance. `run.track_start` has no Genesis4 equivalent |
| `&track fft_fieldsolver` | (always) | Lucifer transcribes the FFT solver only |
| `&importdistribution slicewidth` | `beam.resample.slice_width` | Same meaning, sampling window over bunch length |
| `&importdistribution` (match, center) | no counterpart | A Bmad lattice matches the beam, so these are the lattice's job |
| `&sponrad` | `bmad.radiation_damping_on`, `bmad.radiation_fluctuations_on` | Lucifer exposes Bmad's own switches directly, as Tao does, rather than an FEL-local flag. Both default off, matching Genesis4's `&sponrad` default |
| `&write beam`, `&write field` | `beam.file`, `field.file` on the input side; openPMD dumps on the output side | See file exchange below |
| `&track beamdumpstep`, `fielddumpstep` | `run.dump_at_comb` with `run.comb_ds_save` | Genesis4 strides in integration steps and separately for the beam and the field. Here one comb sets the cadence in metres and both files are written on it, so a frame and its statistics row share an index. A field that is wanted less often is read less often rather than written less often |
| `&importbeam`, `&importfield` | `beam.file`, `field.file` | Lucifer reads openPMD only, so convert first |
| `&lattice`'s `ax`, `ay` on an undulator | the element's `x_offset`, `y_offset` | The undulator's transverse offset, Bmad's own misalignment attributes |
| `&lattice`'s `kx`, `ky` on an undulator | the element's `kx` | Bmad's planar roll-off wavenumber: Genesis4's pair is `-kx^2` and `kx^2 + k_u^2` in units of `k_u^2`, and the helicity defaults otherwise |
| `corrector` | the element's `hkick`, `vkick`, or a `kicker` element | A kick on the FEL element is carried as Bmad's tracking carries it. A `kicker` between elements is an interlude |
| (no equivalent) | `run.migrate` | Slice migration. Off by default, since the comparison tiers run against a code that never migrates |
| (no equivalent) | `fel_method = averaged` or `unaveraged` | Selects the method per element, and `run.transport_model` selects the averaged method's maps |

## Conventions that differ

- **Energy.** Lucifer reports energy as $E = \gamma m_e c^2$ in eV, which is Bmad's convention. Genesis4's `energy` output is $\gamma$. The comparison scripts convert.
- **`s` and `z`.** In Bmad, `s` is position along the lattice and `z` is the longitudinal phase-space coordinate. The statistics file uses both words in exactly that sense. Genesis4 uses `s` for the position within the time window, which is closest to Lucifer's `z_slice`.
- **Units in files.** openPMD's `unitSI` is a factor a reader must apply. The `bmad-stats` file's `@unit` is documentation, and its values are already SI and eV. The two conventions never appear in one file.
- **Field amplitude.** Genesis4 stores an internal field unit. Lucifer stores V/m. Every physical-unit formula here was derived through the identity between the two, taken from [`writeFieldHDF5.cpp`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/writeFieldHDF5.cpp#L70).
- **Particle ordering.** Neither code sorts without `one4one`, so final dumps compare particle by particle rather than only statistically.

## Moving files between the codes

`tests/scripts/convert_genesis.py` holds all of the format knowledge, in both directions:

```
python3 convert_genesis.py <in> <out>
```

Fields convert both ways. Genesis4 particles convert inward through openPMD-beamphysics, which owns those conversions. The fourth direction, openPMD particles out to a Genesis4 `.par.h5`, is implemented here, since no upstream implementation exists. Genesis4's `&importbeam` reads that file, so this is how Genesis4 restarts from a state Lucifer produced.

A file that is not openPMD is refused, with the conversion command in the message.

## Facts about Genesis4 this work pinned down

The transcription had to establish these properties of Genesis4 by reading its source, and they are recorded here because its documentation does not state them. Each is true of the release named at the top of this page.

- Outside undulators Genesis4 does not subdivide into `delz` steps. Each interlude element is one integration step of the element's full length, since [`Lattice::unrollLattice`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Lattice/Lattice.cpp#L340) pushes one entry per non-undulator layout segment. The step count over the Benchmark1-SASE line is 12*89 undulator steps + 36 interlude steps = 1104.
- In steady state `slippage` and `phaseshift` are no-ops. [`Control::applySlippage`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Core/Control.cpp#L158-L159) returns immediately when the run is not time dependent, and the phaseshift array is all zero without phase shifter elements.
- Time dependence follows from `&time` or from importing a multi-slice dump. [`readBeamHDF5.cpp:68-77`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/readBeamHDF5.cpp#L68-L77) reconstructs the window from `slicespacing` and `slicecount`, with time on by default, and `nslice = round(slen/(sample*lambda0))` ([`GenTime.cpp:70-71`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Main/GenTime.cpp#L70-L71)).
- The end-of-lattice autophasing is unguarded. The last step always gets `floor(Lz/(2*gamma0^2*lambda)) + 1` wavelengths of slippage, the `+1` included when there is no trailing drift ([`Lattice.cpp:191-193`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Lattice/Lattice.cpp#L191-L193)).
- The field record's rotation never appears in Genesis4's outputs. [`writeFieldHDF5.cpp`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/writeFieldHDF5.cpp#L86) and [`DiagField::getValues`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Core/Diagnostic.cpp#L755) both unrotate on the fly, so the per-slice arrays of `.out.h5` and the `.fld.h5` dumps are in time-window order, aligned with the beam-slice indexing.
- A helical undulator defaults to `kx = ky = 0.5` ([`LatticeParser.cpp:329-331`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Lattice/LatticeParser.cpp#L329-L331)), scaled by `ku^2` in the unroll.
- The `&importbeam` / `&importfield` namelists take full filenames and make a shared start possible. `&write` before `&track` produces the dumps.
- Genesis4's internal field unit is `dfl [sqrt(W)] = u * dgrid * eev / (ks * sqrt(vacimp))` ([`writeFieldHDF5.cpp:70`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/IO/writeFieldHDF5.cpp#L70)), with `vacimp = 376.73` truncated and `eev = 510998.95069`, which is Bmad's `m_electron`. Equivalently, `u = E*ks/(sqrt(2)*eev)` holds with E in V/m, and that relation derived this tracker's physical-unit formulas.
- Releases up to v4.6.14 use `eev = 510999.06` for the electron rest energy, which is 2.14e-7 above the CODATA 2022 value. From v4.6.15 the value is 510998.95069, which is also Bmad's `m_electron`. A comparison against the current release has no electron-mass discrepancy, and one against an older release has this one. The field unit is linear in `eev`, so the difference reaches the field-normalization integrals as 4.28e-7. It loosens the transcription tiers by several percent of their own value without failing any of them. The levels recorded in [](validation.md) belong to a reference that carries the CODATA value.
- From v4.6.15 every stochastic stream is a pure function of the user seed, the global slice index and a per-stream identifier, mixed with the splitmix64 finaliser. A slice's noise therefore no longer depends on which MPI rank owns it. This project contributed that change. Two consequences follow for a comparison. A result is reproducible across rank counts, where before a fixed seed gave different physics at different core counts and the per-slice bunching could differ by 90 per cent. And the realizations differ from every earlier release, because `ShotNoise` is re-keyed per slice instead of drawn continuously, and [`RandomU::set`](https://github.com/svenreiche/Genesis-1.3-Version4/blob/v4.6.16/src/Util/RandomU.cpp#L31), which had an empty body, now re-keys as the surrounding code always intended. A level recorded against an earlier release compares against a different noise realization. It is neither tighter nor looser.
- From v4.6.16 `FieldSolverFFT` propagates the transverse field on the wavenumbers the transform carries at either grid parity. Releases to v4.6.15 built the table as an offset index, which lands on whole multiples of the fundamental only when `ngrid` is odd, so on an even grid every mode was propagated half a wavenumber out and the constant mode at half a wavenumber rather than at zero, with the error growing with the mode. This project contributed that change too. Until it, every comparison against Genesis4 ran an odd grid, and since it they run the power-of-two grids the rest of the tree runs ([](validation.md#val-the-transverse-grid)). A comparison against v4.6.15 or earlier on an even grid measures that defect, 6.4e-5 of the power on a short planar segment at 64 points and 1.4e-3 on the harmonic tier at 256.
- The quad transport is chromatic through per-particle `gammaz` with `foc^2 = k1*gamma0/gammaz`. Bmad's equivalent scaling is `k1*p0/p`. The two differ by ~4e-9 relative at this energy (1 - beta0), far below the path-length-term difference.

Genesis4's truncated impedance constants set the floor of every transcription comparison. Nine of the eleven tiers agree at that floor and not at round-off. [](validation.md) says which tiers those are and what the remaining two cost.
