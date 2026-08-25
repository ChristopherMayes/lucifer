#!/usr/bin/env python3
"""
Read a radiation field dump into per-slice complex arrays, in either format.

The companion of beamio for the other kind of dump. The tracker writes openPMD
EXT_Wavefront (.wf.h5) and Genesis writes its own .fld.h5, and the checks compare the two
codes slice by slice, so the reading of both belongs in one place rather than in every
script.

What comes back: the field as (nslice, ny, nx) complex V/m, with the grid spacings, the
wavelength and the dump's position along the lattice. Both formats are read with h5py
alone, which keeps this an independent reader of our own writer: nothing here shares a
line with the Fortran that produced the file.

UNITS ARE THE ONE REAL DIFFERENCE. openPMD carries the electric field in V/m, the SI
quantity, while Genesis carries a dimensionless amplitude dfl in sqrt(W) with

    E = dfl * sqrt(2 * Z0) / dx,      Z0 = mu_0 * c

which scales with the grid spacing, so it is not a change of units alone. This module
converts the Genesis side, leaving every comparison in V/m. Any common factor cancels in a
relative difference, so the choice of which side to convert cannot move a measurement; the
choice of Z0 can, at the 1e-9 level of the two codes' impedances, and Bmad's is used here.

A POLARIZATION NOTE. One openPMD file holds both transverse components as components x and
y of the electricField record, where Genesis holds one component per file. Callers say
which they want, and asking for a component the file does not carry is an error rather than
a silent zero.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np

C_LIGHT = 2.99792458e8
MU0_C = 1.25663706127e-6 * C_LIGHT       # Bmad's mu_0_vac (2018 CODATA), not 4 pi e-7.
H_PLANCK_EVS = 4.135667696e-15           # Bmad's h_planck [eV s].
E_CHARGE = 1.602176634e-19

MESH_PATH = "data/1/meshes/electricField"


def is_openpmd(path):
    """True for an openPMD file, by the root attribute the standard requires."""
    with h5py.File(path) as h5:
        return "openPMD" in h5.attrs


def components(path):
    """The polarizations a field dump carries, ('x',) or ('x', 'y')."""
    with h5py.File(path) as h5:
        if "openPMD" not in h5.attrs:
            return ("x",)
        m = h5[MESH_PATH]
        return tuple(name for name in ("x", "y") if name in m)


def read_field(path, component="x", nslice=None):
    """
    A field dump as a dict: u (nslice, ny, nx) complex V/m, dx, dy, dz [m], wavelength [m],
    s_position [m], and components, the polarizations the file carries.

    nslice, when given, is required to match the file's slice count. A caller that knows
    what it asked the tracker for gets the mismatch named here rather than a silently
    truncated comparison.
    """
    path = pathlib.Path(path)
    out = _read_openpmd(path, component) if is_openpmd(path) else _read_genesis(path, component)

    if nslice is not None and out["u"].shape[0] != nslice:
        raise ValueError(f"{path.name} holds {out['u'].shape[0]} slices, caller expected "
                         f"{nslice}")
    return out


def _read_genesis(path, component):
    if component != "x":
        raise ValueError(f"{path.name} is a Genesis field dump, which carries one "
                         f"polarization; component {component!r} was asked for")
    with h5py.File(path) as h5:
        n = int(h5["gridpoints"][0])
        dx = float(h5["gridsize"][0])
        nslice = int(h5["slicecount"][0])
        u = np.empty((nslice, n, n), dtype=complex)
        for i in range(nslice):
            g = h5[f"slice{i + 1:06d}"]
            u[i] = (g["field-real"][:].reshape(n, n)
                    + 1j * g["field-imag"][:].reshape(n, n))
        return dict(u=u * np.sqrt(2 * MU0_C) / dx, dx=dx, dy=dx,
                    dz=float(h5["slicespacing"][0]),
                    wavelength=float(h5["wavelength"][0]),
                    s_position=float(h5["refposition"][0]), components=("x",))


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
        e_photon = float(np.atleast_1d(m.attrs["photonEnergy"])[0])
        return dict(u=u, dx=spacing[2], dy=spacing[1], dz=spacing[0],
                    wavelength=H_PLANCK_EVS * C_LIGHT * E_CHARGE / e_photon,
                    s_position=float(np.atleast_1d(m.attrs["zCoordinate"])[0]),
                    components=components)


def field_power(u, dx, dy):
    """Radiated power per slice [W] from the field in V/m: |E|^2 dx dy summed over the grid."""
    axes = tuple(range(1, u.ndim))
    return (np.abs(u) ** 2).sum(axis=axes) * dx * dy / (2 * MU0_C)
