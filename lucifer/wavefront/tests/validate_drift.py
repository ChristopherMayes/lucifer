#!/usr/bin/env python3
"""
Python half of the Fortran against Python free-space drift comparison.

Run through run_validation.sh, which builds the Fortran side, sets PYTHONPATH to an
openPMD-beamphysics checkout and calls this script twice:

    validate_drift.py write   <input.h5>  --nx ... --ny ... --nz ... --wavelength ...
    validate_drift.py compare <input.h5> <fortran_out.h5> --z ...

The first call builds a test wavefront and writes it as an openPMD EXT_Wavefront file.
The second reads the same input, drifts it with openPMD-beamphysics' drift_wavefront,
reads the file the Fortran program produced, and reports the largest relative difference
between the two.

The input is deliberately not a centered symmetric Gaussian. A symmetric field hides a
whole class of mistakes: a transposed transverse grid, x and y swapped somewhere in the
HDF5 index ordering, or a sign error in the kernel's negative wavenumbers all leave a
symmetric field unchanged. The test field is therefore offset from the axis, elliptical,
tilted, and different in every slice.

The input carries BOTH transverse polarizations, with a different field in each. One
openPMD file holds both as components of the electricField record, which is what the
Genesis4 format could not do, so the two-component path through the Fortran reader,
the propagator and the Fortran writer is compared against Python here rather than only
against itself. Giving Ey its own centroid, size and tilt is what makes a swap of the
two components visible.

Two conventions the Fortran side fixes, and this script matches so that the whole loop
is information preserving. The transverse grid is centered on the axis. The slice axis
starts at zero, so the head of the radiation window is z = 0 rather than the window
being centered. Both reach the file as gridGlobalOffset, which openPMD defines as the
position of the first cell.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from dataclasses import replace

import h5py
import numpy as np

from beamphysics.wavefront import Wavefront
from beamphysics.wavefront.propagators import drift_wavefront

# EXT_Wavefront states a photon energy rather than a wavelength, and the two codes convert
# between them with different Planck constants. Bmad's h_planck is 4.135667696e-15 eV s,
# truncated from the exact post-2019 value of h/e, which is 4.135667696923859e-15 and is
# what scipy carries. So the same file read by the two readers gives wavelengths differing
# by 2.2339e-10 of each other.
#
# The drift kernel's phase is z * lambda * (kx^2 + ky^2) / (4 pi), which reaches 20 rad at
# the corner of the default grid, so that difference alone puts the drifted fields 7.5e-11
# apart: a hundred times the tolerance, and nothing to do with either propagator. The
# Python side is therefore moved onto Bmad's wavelength before it propagates, which leaves
# the comparison measuring the propagators and the file I/O. The constants themselves are
# then visible in exactly one place, the photonEnergy attribute, where the file comparison
# allows this much and no more.
BMAD_WAVELENGTH_RATIO = 4.135667696e-15 / 4.135667696923859e-15
PLANCK_TOL = 1.0e-9


def build_wavefront(nx: int, ny: int, nz: int, dx: float, dy: float, dz: float,
                    wavelength: float, s_position: float) -> Wavefront:
    """
    Build the test wavefront: off-axis tilted elliptical Gaussians whose centroid, size
    and phase all vary from slice to slice, one in each transverse polarization.

    Asymmetry is the point. Every deliberate feature below breaks a symmetry that would
    otherwise let a bug through:

      - x0 != y0 and sigma_x != sigma_y  break the x/y interchange symmetry, so a
        transposed slice or a swapped dx/dy shows up.
      - the xy cross term tilts the ellipse, so a transpose is visible even if the two
        sizes were equal.
      - the linear phase gives the field a nonzero mean transverse wavenumber, so the
        centroid moves during the drift and the sign of the kernel matters.
      - the quadratic phase makes the field converge in x and diverge in y, so the two
        transverse axes evolve differently.
      - varying everything with slice index means a mistake in the slice axis, or a
        slice written into the wrong group, cannot cancel.
      - Ey is not Ex, nor a transpose of it, so interchanging the two components is
        visible even where the grid is square.

    s_position is the lattice position of the dump plane, which the file carries as the
    extension's zCoordinate. It is nonzero here so that a Fortran side failing to
    advance it over the drift cannot pass.
    """
    x = (np.arange(nx) - (nx - 1) / 2) * dx
    y = (np.arange(ny) - (ny - 1) / 2) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")

    k0 = 2 * np.pi / wavelength
    w = 8 * dx  # Reference transverse scale, comfortably resolved and well inside the grid.

    def gaussian_pulse(x0_scale, y0_scale, sx_scale, sy_scale, cross, tilt_x, tilt_y,
                       curv_x, curv_y, slice_phase, amp0):
        """One polarization: a tilted ellipse whose parameters walk along the pulse."""
        field = np.empty((nx, ny, nz), dtype=np.complex128)

        for iz in range(nz):
            f = iz / max(nz - 1, 1)  # 0 .. 1 across the pulse.

            x0 = (x0_scale + 0.4 * f) * w  # Off axis, and moving with slice.
            y0 = y0_scale * w
            sx = (sx_scale + 0.5 * f) * w
            sy = sy_scale * w  # Elliptical.

            u = (X - x0) / sx
            v = (Y - y0) / sy

            # Amplitude: tilted ellipse (the u*v term), amplitude varying along the pulse.
            amp = (amp0 + f) * np.exp(-(u**2 + v**2 + cross * u * v) / 2)

            # Phase: a linear part that steers the beam, and a quadratic part with
            # opposite signs in x and y so the two axes focus and defocus respectively.
            tilt = k0 * (tilt_x * (X - x0) + tilt_y * (Y - y0))
            curv = k0 * ((X - x0) ** 2 / (2 * curv_x) - (Y - y0) ** 2 / (2 * curv_y))

            field[:, :, iz] = amp * np.exp(1j * (tilt + curv + slice_phase * iz))

        return field

    Ex = gaussian_pulse(0.9, -1.3, 1.0, 1.7, 0.6, 2.5e-5, -1.1e-5, 4.0, 7.0, 0.37, 0.4)
    Ey = gaussian_pulse(-1.1, 0.7, 1.4, 0.9, -0.5, -1.7e-5, 3.1e-5, 9.0, 3.0, -0.53, 0.7)

    # Scale to a physically plausible field amplitude in V/m, the same factor on both
    # components so that their relative size is the one built above.
    scale = 1e9 / max(np.abs(Ex).max(), np.abs(Ey).max())

    # zmid puts the first slice at z = 0: the window's head, not its center.
    return Wavefront(Ex=scale * Ex, Ey=scale * Ey, dx=dx, dy=dy, dz=dz,
                     wavelength=wavelength, zmid=(nz - 1) * dz / 2,
                     s_position=s_position)


def file_layout(path: str) -> dict[str, tuple]:
    """
    Walk an HDF5 file and return {name: (kind, shape)} for every object in it, where kind is
    'group' or the NumPy dtype kind character of a dataset.

    The dtype is recorded in full so that width differences can be reported, but see
    compare_file for which differences are treated as failures.
    """
    layout: dict[str, tuple] = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Group):
            layout[name] = ("group", None)
        else:
            layout[name] = (str(obj.dtype), obj.shape)

    with h5py.File(path, "r") as h5:
        h5.visititems(visit)

    return layout


def file_attributes(path: str) -> dict[tuple, object]:
    """
    Every attribute in the file, as {(object name, attribute name): value}, with the
    root group under the name '/'.

    An openPMD file states nearly all of its metadata in attributes rather than in
    datasets, so this is where a comparison of two writers has to look.
    """
    attrs: dict[tuple, object] = {}

    def collect(name, obj):
        for key, value in obj.attrs.items():
            attrs[(name, key)] = value

    with h5py.File(path, "r") as h5:
        collect("/", h5)
        h5.visititems(collect)

    return attrs


def _text(value):
    """A string attribute as text, whether it arrived as bytes or as str, scalar or array."""
    flat = np.atleast_1d(value)
    out = [v.decode() if isinstance(v, bytes) else str(v) for v in flat]
    return out[0] if len(out) == 1 else tuple(out)


def _attrs_agree(got, want, tol):
    """
    Whether two attribute values carry the same content.

    Two differences are content preserving and are accepted here rather than reported.
    Bmad's HDF5 helpers write every scalar attribute as a length-1 array, where h5py
    writes a true scalar, so both sides are flattened before comparison; every openPMD
    reader involved, ours and openPMD-beamphysics', reads through np.atleast_1d. And a
    string reaches the file as fixed length from Fortran and as variable length from
    h5py, which HDF5 hands back as the same text either way.
    """
    got_arr, want_arr = np.atleast_1d(got), np.atleast_1d(want)

    if want_arr.dtype.kind in "SUO" or got_arr.dtype.kind in "SUO":
        return _text(got) == _text(want)

    if got_arr.shape != want_arr.shape:
        return False

    scale = np.maximum(np.abs(want_arr), 1.0)
    return bool(np.all(np.abs(got_arr - want_arr) <= tol * scale))


def compare_file(fortran_file: str, reference_file: str) -> bool:
    """
    Require the Fortran output file to hold the same objects, shapes and attributes as a
    file written by Wavefront.write_openpmd.

    This is not redundant with the field comparison. openPMD states the grid spacing, the
    grid offset, the axis order, the photon energy and the position along the lattice in
    attributes, and a reader that finds them missing or wrong cannot reconstruct the field
    it is looking at however exact the numbers in the datasets are. The field comparison
    reads both files with the same reader, which papers over exactly that: a wrong
    gridSpacing on both sides cancels. So the metadata is checked against an independent
    writer here, name by name and value by value.

    A dtype difference is reported unless _TYPE_DIFF_ALLOWED explains it.
    """
    got, want = file_layout(fortran_file), file_layout(reference_file)
    got_attrs, want_attrs = file_attributes(fortran_file), file_attributes(reference_file)

    ok = True

    missing = sorted(set(want) - set(got))
    extra = sorted(set(got) - set(want))

    if missing:
        print(f"  FAIL: objects the reference has and the Fortran output does not: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")
        ok = False
    if extra:
        print(f"  FAIL: objects the Fortran output has and the reference does not: "
              f"{', '.join(extra[:8])}{' ...' if len(extra) > 8 else ''}")
        ok = False

    noted = []
    for name in sorted(set(got) & set(want)):
        g_type, g_shape = got[name]
        w_type, w_shape = want[name]

        if g_shape != w_shape:
            print(f"  FAIL: {name}: Fortran shape {g_shape}, reference shape {w_shape}")
            ok = False
            continue

        if g_type == w_type:
            continue

        if name in _TYPE_DIFF_ALLOWED:
            noted.append(f"{name}: {g_type} against {w_type} ({_TYPE_DIFF_ALLOWED[name]})")
        else:
            print(f"  FAIL: {name}: Fortran dtype {g_type}, reference dtype {w_type}")
            ok = False

    # Attributes: the same names on the same objects, carrying the same values.

    missing_attrs = sorted(set(want_attrs) - set(got_attrs))
    extra_attrs = sorted(set(got_attrs) - set(want_attrs))

    for where, what in ((missing_attrs, "the reference has and the Fortran output does not"),
                        (extra_attrs, "the Fortran output has and the reference does not")):
        named = [f"{obj}:{key}" for obj, key in where if not _skipped(obj, key)]
        if named:
            print(f"  FAIL: attributes {what}: {', '.join(named)}")
            ok = False

    for key in sorted(set(got_attrs) & set(want_attrs)):
        if _skipped(*key):
            continue
        tol = _ATTR_TOL.get(key[1], (1.0e-14, ""))
        if not _attrs_agree(got_attrs[key], want_attrs[key], tol[0]):
            print(f"  FAIL: attribute {key[0]}:{key[1]}: Fortran {got_attrs[key]!r}, "
                  f"reference {want_attrs[key]!r}, relative tolerance {tol[0]:.1e}"
                  + (f" [{tol[1]}]" if tol[1] else ""))
            ok = False

    if ok:
        n_comp = sum(1 for k in want if k.startswith(MESH_PATH + "/"))
        n_attr = sum(1 for k in want_attrs if not _skipped(*k))
        print(f"  {len(want)} objects match by name and shape, including {n_comp} "
              f"polarization components of the electricField record")
        print(f"  {n_attr} attributes match by name and value, over the series root, the "
              f"iteration, the mesh record and the components")
    for note in noted:
        print(f"  note: {note}")

    return ok


# Dataset types that may differ between the Fortran writer and Wavefront.write_openpmd
# without that being an incompatibility. The complex compound is the only one: HDF5 has
# no native complex type, so both writers build a two-field compound, and h5py presents
# either as complex128 on read while reporting the field names it was built from.
_TYPE_DIFF_ALLOWED: dict[str, str] = {}

# Attributes not compared, and why. All three are written by Bmad's shared openPMD dataset
# writer, which stamps the unit record of every dataset it creates, and none is written by
# openPMD-beamphysics on a mesh component. They carry no part of the field or its geometry,
# an openPMD reader passes over an attribute it does not know, and the required unitSI sits
# beside them and IS compared.
#
#   localName       the descriptive name Bmad gives the dataset, here 'x' or 'y'
#   unitSymbol      'V/m', the same statement unitDimension makes in exponents
#   unitDimension   the component's dimensions. The record carries the authoritative
#                   copy, which is compared; this one repeats it.
#
# The allowance is on a component and nowhere else, so the record's own unitDimension is
# still compared name by name and value by value.
MESH_PATH = "data/1/meshes/electricField"

_COMPONENT_EXTRA = {
    "localName": "written by Bmad's dataset writer, not by openPMD-beamphysics",
    "unitSymbol": "written by Bmad's dataset writer, not by openPMD-beamphysics",
    "unitDimension": "on a component: Bmad's dataset writer repeats the record's value",
}


def _skipped(obj_name, attr_name):
    """Whether this attribute is one of the documented allowances above."""
    return obj_name.startswith(MESH_PATH + "/") and attr_name in _COMPONENT_EXTRA

# Attributes compared to something other than round-off, with the cause. photonEnergy is
# the only one: it is where the two Planck constants land, since each side writes the energy
# its own wavelength implies. See BMAD_WAVELENGTH_RATIO.
_ATTR_TOL = {
    "photonEnergy": (PLANCK_TOL, "Bmad and scipy carry different Planck constants"),
}


def parse_check_lines(path: str) -> dict[str, float]:
    """
    Pull the 'CHECK <name> <value>' lines out of the Fortran program's stdout.
    """
    values: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) == 3 and fields[0] == "CHECK":
                values[fields[1]] = float(fields[2])
    return values


def report(name: str, a: np.ndarray, b: np.ndarray) -> float:
    """
    Compare two complex arrays and print the difference measures. Returns the peak
    normalized difference, that is, max|a-b| divided by max|b|.

    Peak normalization rather than elementwise relative error is the meaningful measure
    here. The field spans many orders of magnitude across the grid, and the round-off
    floor of an FFT based propagator is set by the largest value in the transform, not by
    the local one; an elementwise ratio in the far tail measures nothing but that floor
    divided by a number near zero.
    """
    if a is None or b is None:
        print(f"  {name}: absent on one side, python {b is not None}, fortran {a is not None}")
        return np.inf

    scale = np.abs(b).max()
    if scale == 0:
        print(f"  {name}: reference is identically zero, nothing to compare")
        return np.inf

    abs_diff = np.abs(a - b)
    peak_norm = abs_diff.max() / scale

    # Elementwise relative difference, restricted to points carrying real amplitude, as a
    # cross-check that the agreement is not confined to the peak.
    mask = np.abs(b) > 1e-6 * scale
    elem_rel = (abs_diff[mask] / np.abs(b)[mask]).max() if mask.any() else np.nan

    print(f"  {name}:")
    print(f"    max |fortran - python|            = {abs_diff.max():.6e}")
    print(f"    max |python|                      = {scale:.6e}")
    print(f"    max |diff| / max |python|         = {peak_norm:.6e}   <-- largest relative difference")
    print(f"    max elementwise rel diff (>1e-6)  = {elem_rel:.6e}   ({mask.sum()} of {b.size} points)")

    return peak_norm


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pw = sub.add_parser("write", help="Build the test wavefront and write it as openPMD")
    pw.add_argument("input_file")
    pw.add_argument("--nx", type=int, default=64)
    pw.add_argument("--ny", type=int, default=None, help="Default: nx")
    pw.add_argument("--nz", type=int, default=5)
    pw.add_argument("--dx", type=float, default=2.0e-6)
    pw.add_argument("--dy", type=float, default=None, help="Default: dx")
    pw.add_argument("--dz", type=float, default=1.0e-9)
    pw.add_argument("--wavelength", type=float, default=1.0e-9)
    pw.add_argument("--s-position", type=float, default=3.75,
                    help="Lattice position of the dump plane [m], the extension's "
                         "zCoordinate. The Fortran side must advance it by the drift.")

    pc = sub.add_parser("compare", help="Drift in Python and compare against the Fortran output")
    pc.add_argument("input_file")
    pc.add_argument("fortran_file")
    pc.add_argument("--z", type=float, required=True)
    pc.add_argument("--tolerance", type=float, default=1.0e-12,
                    help="Maximum acceptable peak normalized difference")
    pc.add_argument("--fortran-log", default=None,
                    help="File holding the Fortran program's stdout. Its 'CHECK <name> <value>' "
                         "lines are recomputed here and verified.")
    pc.add_argument("--scalar-tolerance", type=float, default=1.0e-11,
                    help="Maximum acceptable relative difference for the CHECK scalars")

    args = p.parse_args()

    if args.command == "write":
        # openPMD puts no constraint on the grid: nx may differ from ny and dx from dy,
        # which the Genesis4 format forbade.
        w = build_wavefront(nx=args.nx, ny=args.ny if args.ny else args.nx, nz=args.nz,
                            dx=args.dx, dy=args.dy if args.dy else args.dx, dz=args.dz,
                            wavelength=args.wavelength, s_position=args.s_position)
        w.write_openpmd(args.input_file)
        print(f"Wrote test wavefront to {args.input_file}")
        print(f"  shape {w.shape}, dx {w.dx:.6e} m, dy {w.dy:.6e} m, dz {w.dz:.6e} m, "
              f"wavelength {w.wavelength:.6e} m")
        print(f"  both polarizations present, s_position {w.s_position:.6e} m, "
              f"first slice at z = {w.zmin:.6e} m")
        print(f"  energy {w.energy:.6e} J, sigma_x {w.sigma_x:.6e} m, sigma_y {w.sigma_y:.6e} m")
        return 0

    # Compare.

    print("Python drift_wavefront reference")

    # Read the same file the Fortran program read, so that the two start from bitwise
    # identical input rather than from two separate constructions of the same field.
    w_file = Wavefront.from_openpmd(args.input_file)
    w_in = replace(w_file, wavelength=w_file.wavelength * BMAD_WAVELENGTH_RATIO,
                   attrs=w_file.attrs.copy())
    w_py = drift_wavefront(w_in, args.z)

    w_f = Wavefront.from_openpmd(args.fortran_file)

    print(f"  input:   {args.input_file}")
    print(f"  fortran: {args.fortran_file}")
    print(f"  drift:   {args.z:.6e} m")
    print(f"  shape:   python {w_py.shape}, fortran {w_f.shape}")
    print(f"  wavelength: {w_file.wavelength:.15e} m as the file's photonEnergy reads here, "
          f"{w_in.wavelength:.15e} m as it reads in Bmad")

    if w_py.shape != w_f.shape:
        print(f"FAIL: shape mismatch, python {w_py.shape} vs fortran {w_f.shape}")
        return 1

    # The grid and the photon energy must come back out of the Fortran program exactly as
    # they went in, so these compare its output against its input rather than against the
    # Python wavefront, whose wavelength was moved onto Bmad's constant above.
    for name, val_in, val_f in (("dx", w_file.dx, w_f.dx),
                                ("dy", w_file.dy, w_f.dy),
                                ("dz", w_file.dz, w_f.dz),
                                ("wavelength", w_file.wavelength, w_f.wavelength)):
        if val_in != val_f:
            print(f"FAIL: {name} mismatch, input file {val_in!r} vs fortran {val_f!r}")
            return 1

    # s_position is the lattice position of the plane the field sits on, so a drift must
    # advance it. Nothing in the field arrays can see that, and a reader that trusts the
    # file would place the dump at the wrong place along the line.
    if w_py.s_position != w_f.s_position:
        print(f"FAIL: s_position mismatch, python {w_py.s_position!r} vs "
              f"fortran {w_f.s_position!r}")
        return 1

    print(f"  s_position advanced by the drift: {w_file.s_position:.6e} m -> "
          f"{w_f.s_position:.6e} m, matching Python")

    print()
    print("File content, against a file written by Wavefront.write_openpmd")
    with tempfile.TemporaryDirectory() as tmp:
        ref_file = os.path.join(tmp, "python_reference.wf.h5")
        w_py.write_openpmd(ref_file)
        structure_ok = compare_file(args.fortran_file, ref_file)

    print()
    print("Field comparison, drifted")
    peak_norm = max(report("Ex", w_f.Ex, w_py.Ex), report("Ey", w_f.Ey, w_py.Ey))

    # An undrifted round trip, to separate an I/O error from a propagator error. If this
    # one is large but the drifted one is not, or vice versa, that says which half is wrong.
    print()
    print("Field comparison, input file read back (isolates HDF5 round trip from the drift)")
    w_again = Wavefront.from_openpmd(args.input_file)
    report("Ex", w_again.Ex, w_in.Ex)
    report("Ey", w_again.Ey, w_in.Ey)

    print()
    print(f"Energy: python {w_py.energy:.12e} J, fortran output file {w_f.energy:.12e} J, "
          f"relative difference {abs(w_f.energy - w_py.energy) / w_py.energy:.6e}")

    # The Fortran program's own values for the mirrored derived quantities, recomputed here
    # from the file it wrote. Without this the Fortran implementations of wavefront_energy,
    # wavefront_transverse_moments and the coordinate vectors underneath them would be
    # printed but never checked against anything.

    scalars_ok = True
    if args.fortran_log:
        print()
        print("Derived quantities computed by Fortran, recomputed here from its output file")

        # One of these cannot reach the field comparison's round-off floor, for a reason that
        # is not a defect in either implementation. It gets a tolerance sized to its cause,
        # rather than one loose tolerance that would hide the next real problem.
        #
        # energy. A sum over every grid point of both polarizations, which is 40960 of them
        # at the default size. NumPy sums pairwise, with error growing as log(N) * eps. The
        # Fortran sum intrinsic sums sequentially, with error growing as N * eps, and N * eps
        # here is 9.1e-12. The two also differ in association: Fortran forms the intensity in
        # W/m^2 and divides the sum by c, where Python forms the energy density directly.
        #
        # photon_energy is held to the same tolerance as the rest even though the two codes
        # convert it with different Planck constants. The file states the energy itself, and
        # each side recovers it with the constant it went in with, so the difference cancels
        # here. It shows up where it belongs: in the photonEnergy attribute of a file the
        # two writers wrote from wavelengths of their own.
        expected = {
            "photon_energy": (w_f.photon_energy, args.scalar_tolerance, ""),
            "energy":        (w_f.energy, 1.0e-10,
                              "sequential against pairwise summation over all grid points"),
            "mean_x":        (w_f.mean_x, args.scalar_tolerance, ""),
            "mean_y":        (w_f.mean_y, args.scalar_tolerance, ""),
            "sigma_x":       (w_f.sigma_x, args.scalar_tolerance, ""),
            "sigma_y":       (w_f.sigma_y, args.scalar_tolerance, ""),
        }
        reported = parse_check_lines(args.fortran_log)

        missing = sorted(set(expected) - set(reported))
        if missing:
            print(f"  FAIL: the Fortran log reported no value for: {', '.join(missing)}")
            scalars_ok = False

        for name, (want, tol, note) in expected.items():
            if name not in reported:
                continue
            got = reported[name]
            denom = abs(want) if want != 0 else 1.0
            rel = abs(got - want) / denom
            status = "ok" if rel <= tol else "FAIL"
            if status == "FAIL":
                scalars_ok = False
            print(f"  {name:<14} fortran {got:+.12e}   python {want:+.12e}   "
                  f"rel {rel:.3e}  tol {tol:.1e}  {status}"
                  + (f"   [{note}]" if note else ""))

    print()
    print(f"LARGEST RELATIVE DIFFERENCE: {peak_norm:.6e}   (tolerance {args.tolerance:.1e})")

    if not np.isfinite(peak_norm) or peak_norm > args.tolerance or not scalars_ok or not structure_ok:
        print("FAIL")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
