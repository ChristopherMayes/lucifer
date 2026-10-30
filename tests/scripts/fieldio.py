#!/usr/bin/env python3
"""
The check scripts' name for the package's field readers, with the Genesis writer.

read_field, is_openpmd and components are lucifer.field's, imported here, so the check
scripts keep the name they have always used. The tracker writes openPMD EXT_Wavefront
(.wf.h5) and Genesis writes its own .fld.h5, and the checks compare the two codes slice
by slice, so the reading of both is in one place, the package, and a user comparing a
run against Genesis reads the same arrays the checks read. The constants are the
package's too, so write_genesis below inverts read_field with the same Z0.

What comes back from read_field: the field as (nslice, ny, nx) complex V/m, with the grid
spacings, the wavelength, the coordinate of the first slice within the pulse and, for
openPMD, the plane's position along the lattice. lucifer.field's docstring gives the unit
conversion of the Genesis side.
"""

from __future__ import annotations

import h5py
import numpy as np

from lucifer.field import C_LIGHT, E_CHARGE, H_PLANCK_EVS, MESH_PATH, MU0_C  # noqa: F401
from lucifer.field import components, is_openpmd, read_field  # noqa: F401


def write_genesis(path, field):
    """
    Write one polarization of a field, as read_field returns it, in Genesis's own layout:
    the root datasets gridpoints, gridsize, refposition (zmin), wavelength, slicecount and
    slicespacing, and a group per slice with field-real and field-imag in the dfl scale.
    The inverse of the Genesis side of read_field, with the same constants, so a field
    written here and read back agrees to rounding. A check that hands Genesis a field
    the program wrote writes the copy with this rather than with a converter whose
    constants differ from Bmad's, since Genesis scales the field by the wavelength the
    file states on import and by its own on export.
    """
    u, n = field["u"], field["u"].shape[1]
    if field["u"].shape[2] != n or field["dx"] != field["dy"]:
        raise ValueError("Genesis's layout holds a square grid of one spacing")
    dfl = u * field["dx"] / np.sqrt(2 * MU0_C)
    with h5py.File(path, "w") as h5:
        h5["gridpoints"] = np.array([n], dtype=np.int32)
        h5["gridsize"] = np.array([field["dx"]])
        h5["refposition"] = np.array([field["zmin"]])
        h5["wavelength"] = np.array([field["wavelength"]])
        h5["slicecount"] = np.array([u.shape[0]], dtype=np.int32)
        h5["slicespacing"] = np.array([field["dz"]])
        for i in range(u.shape[0]):
            g = h5.create_group(f"slice{i + 1:06d}")
            g["field-real"] = np.ascontiguousarray(dfl[i].real).reshape(-1)
            g["field-imag"] = np.ascontiguousarray(dfl[i].imag).reshape(-1)


def _read_openpmd(path, component):
    with h5py.File(path) as h5:
        m = h5[MESH_PATH]
        components = tuple(name for name in ("x", "y") if name in m)
        if component not in components:
            raise ValueError(f"{path.name} carries polarization {components} and does not "
                             f"carry {component!r}")
        # Stored (z, y, x), which is what a caller wants slice by slice, so no transpose.
        u = np.asarray(m[component][...], dtype=complex)
        spacing = np.atleast_1d(m.attrs["gridSpacing"]).astype(float)
        offset = np.atleast_1d(m.attrs["gridGlobalOffset"]).astype(float)
        e_photon = float(np.atleast_1d(m.attrs["photonEnergy"])[0])
        return dict(u=u, dx=spacing[2], dy=spacing[1], dz=spacing[0],
                    wavelength=H_PLANCK_EVS * C_LIGHT * E_CHARGE / e_photon,
                    zmin=offset[0],
                    s_position=float(np.atleast_1d(m.attrs["zCoordinate"])[0]),
                    components=components)


def field_power(u, dx, dy):
    """Radiated power per slice [W] from the field in V/m: |E|^2 dx dy / (2 Z0) summed over the grid."""
    axes = tuple(range(1, u.ndim))
    return (np.abs(u) ** 2).sum(axis=axes) * dx * dy / (2 * MU0_C)
