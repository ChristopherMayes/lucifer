# EXT_Wavefront: clarifications from a first cross-code implementation

Notes for the openPMD-standard maintainers (branch `upcoming-2.0.0`,
`EXT_Wavefront.md`), from implementing the extension independently in Fortran (the
Bmad FEL tracker, `bsim/modules/wavefront_openpmd_mod.f90`) and Python (the carried
`wavefront_openpmd_patch.py` beside this file, targeting openPMD-beamphysics'
`Wavefront` class), with the two implementations validated against each other every
harness run. Each item is a place where the extension's text under-determines the
file, resolved here one way; the extension would benefit from saying it explicitly.

1. **Attribute placement contradicts itself.** The heading says "Additional
   attributes on the `mesh record` named `electricField`"; the sentence under it
   says "On the `series` object, set the following attributes". These are different
   HDF5 objects. We resolved in the RECORD's favor: `photonEnergy` is a property of
   one field, and a series-level attribute would forbid a file (or a future
   multi-record layout) from ever carrying more than one. Recommend deleting the
   "series object" sentence.

2. **Say which axis the slice/time axis is.** For a wavefront the openPMD iteration
   cannot be the time sample: slices of one pulse are simultaneous, and propagation
   steps are what an iteration naturally is. The slice axis must therefore be a MESH
   axis, declared by `axisLabels`. We store `(z, y, x)` order (slice axis first in
   the file's C-order view, which is also the natural zero-copy layout for both
   Fortran and numpy writers) and read whatever `axisLabels` declare. The extension
   states none of this and every implementer will guess differently.

3. **The `z` component of `electricField`.** The extension lists components
   `x/y/z`. A paraxial code has no `z` component; absent components are ordinary
   openPMD, but the extension should say whether `z` is required (we say no).

4. **Harmonics / multiple wavelengths.** The record name is fixed
   (`electricField`), so one file holds one `photonEnergy`. Multi-color output
   (FEL harmonics) therefore needs one file per harmonic. If that is intended,
   say so; if not, the extension needs a naming convention for additional records
   (`electricField3`, ...).

5. **`photonEnergy` units.** The attribute's `unitDimension` is given (energy) but
   attributes carry no `unitSI`. We write JOULES (SI). Recommend stating the unit.

6. **Complex storage.** The 2.0 base standard's `complexX` types map naturally to
   an HDF5 compound with members `r` and `i` -- which h5py reads as `complex128`
   natively (no manual reassembly). Worth an "advice to implementors" note.

7. **`temporalDomain = 'frequency'` unit text.** The inverse-length unitDimension
   given for the frequency domain reads as a Genesis-family convention
   (sqrt(J/eV)/m); a sentence deriving it would prevent wrong implementations.
