"""
Read a field dump, and the escaped light a run writes beside its final dump.

``read_field`` reads one polarization of a field dump into per-slice complex arrays, in
V/m, from either format: the openPMD EXT_Wavefront file the program writes
(``.wf.h5``) or the ``.fld.h5`` of Genesis 1.3 Version 4. Both are read with h5py alone.
Units are the one real difference between the two. openPMD carries the electric field in
V/m, the SI quantity, while Genesis carries a dimensionless amplitude dfl in sqrt(W) with

    E = dfl * sqrt(2 * Z0) / dx,      Z0 = mu_0 * c

which scales with the grid spacing, so it is not a change of units alone. The Genesis
side is converted, leaving every comparison in V/m. Any common factor cancels in a
relative difference, so the choice of which side to convert cannot move a result. The
choice of Z0 can, at the 1e-9 level of the two codes' impedances, and Bmad's is used
here. One openPMD file holds both transverse components as components x and y of the
electricField record, where Genesis holds one component per file. Callers say which they
want, and asking for a component the file does not carry is an error rather than a
silent zero.

``<out_root>-pulse.wf.h5`` is the whole pulse at the exit, the window's slices and the
light slippage carried out of the window propagated to the exit plane. A run keeps it
whenever light left the window. ``<out_root>-escaped.wf.h5`` is the bank itself, one
plane per slice as it was when it left, and a run writes it only with
``run.write_escaped``. Both are openPMD EXT_Wavefront files, so ``Wavefront.from_openpmd``
of openpmd-beamphysics reads either whole, with both polarizations as its Ex and Ey.

Each file carries beside its mesh what openPMD has no name for, in a group of its own,
``lucifer/pulse`` or ``lucifer/escaped``: ``zTransmit``, the beamline position at which
each banked plane left the window, and ``wavefrontParams``, each banked plane's moments
as the program computed them then, one row per banked plane in the mesh's order. The two
groups hold the same rows. In the escaped file every plane is banked. In the pulse the
banked planes follow the window's, so row i is plane ``nz - n + i``. ``read_pulse`` and
``read_escaped`` return the Wavefront with those two arrays, and the y component's params
where the run carried one.

No one beamline position holds the escaped planes, so the escaped file's
``zCoordinate``, and the Wavefront's ``s_position``, is NaN. Each plane's own is
``zTransmit``. The pulse file's is the exit plane.
"""

from __future__ import annotations

import dataclasses
import pathlib

import h5py
import numpy as np

__all__ = ["read_field", "is_openpmd", "components", "Escaped", "Pulse", "read_escaped", "read_pulse"]

C_LIGHT = 2.99792458e8
MU0_C = 1.25663706127e-6 * C_LIGHT       # Bmad's mu_0_vac (2018 CODATA), not 4 pi e-7.
H_PLANCK_EVS = 4.135667696e-15           # Bmad's h_planck [eV s].
E_CHARGE = 1.602176634e-19

MESH_PATH = "data/1/meshes/electricField"

ESCAPED_GROUP = "lucifer/escaped"
PULSE_GROUP = "lucifer/pulse"

# The 25 columns of wavefrontParams: the program's wavefront_params_struct flattened.
# The centroid over (x, theta_x, y, theta_y) [m, rad], the 4x4 second-moment matrix over
# the same axes column by column, then the energy [J], the power [W], the on-axis
# intensity [W/m^2] and the two emittances [m rad]. The file states them too, as the
# dataset's columns attribute, and the readers hold the two to agree.
PARAM_NAMES = (
    "centroid_x", "centroid_theta_x", "centroid_y", "centroid_theta_y",
    "sigma_11", "sigma_21", "sigma_31", "sigma_41", "sigma_12", "sigma_22", "sigma_32", "sigma_42",
    "sigma_13", "sigma_23", "sigma_33", "sigma_43", "sigma_14", "sigma_24", "sigma_34", "sigma_44",
    "energy", "power", "on_axis_intensity", "emit_x", "emit_y",
)


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
    zmin [m], s_position [m] or None, and components, the polarizations the file carries.
    Slices are in time-window order, as both codes write them. zmin is the coordinate of
    the first slice within the pulse, Genesis's refposition or openPMD's gridGlobalOffset
    along z, and s_position the plane's position along the lattice, openPMD's zCoordinate,
    which the Genesis format has no record of, so it comes back None there.

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
                    zmin=float(h5["refposition"][0]), s_position=None, components=("x",))


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


@dataclasses.dataclass
class Escaped:
    """The escaped-field bank of one run: the planes as a Wavefront, in increasing
    position along the pulse, and per plane where it left the window and its moments
    then. ``params[:, PARAM_NAMES.index("energy")]`` is each plane's energy."""

    wavefront: object
    """beamphysics.wavefront.Wavefront: the planes, (nx, ny, n) in V/m, s_position NaN."""
    z_transmit: np.ndarray
    """(n,) the beamline position at which each plane left the window [m]."""
    params: np.ndarray
    """(n, 25) each plane's wavefront params at banking, the columns PARAM_NAMES."""
    params_y: np.ndarray | None = None
    """(n, 25) the y component's params, when the run carried a second polarization."""

    def _energy(self):
        e = self.params[:, PARAM_NAMES.index("energy")].sum()
        if self.params_y is not None:
            e += self.params_y[:, PARAM_NAMES.index("energy")].sum()
        return e

    def __repr__(self):
        n = len(self.z_transmit)
        span = f"{self.z_transmit.min():.4g} to {self.z_transmit.max():.4g} m" if n else "none"
        return f"Escaped({n} planes, left at {span}, {self._energy():.4g} J)"


@dataclasses.dataclass(repr=False)
class Pulse(Escaped):
    """The whole pulse of one run at the exit: the window's planes and then the banked
    ones as a Wavefront, in increasing position along the pulse, and per banked plane
    where it left the window and its moments then. The banked planes are
    ``wavefront.Ex[:, :, n_window:]``, and their rows are those of the escaped file."""

    n_window: int = 0
    """The planes from the window, which come first in the mesh."""

    def __repr__(self):
        n = len(self.z_transmit)
        span = f"{self.z_transmit.min():.4g} to {self.z_transmit.max():.4g} m" if n else "none"
        return (f"Pulse({self.n_window} planes of the window and {n} banked, left at {span}, "
                f"{self._energy():.4g} J when they left)")


def _names(dataset):
    return tuple(v.decode() if isinstance(v, bytes) else str(v) for v in dataset.attrs["columns"])


def _read(path, group, what):
    """The Wavefront of ``path`` and the three arrays of its ``group``, refused as the
    readers state."""
    from beamphysics.wavefront import Wavefront

    path = pathlib.Path(path)
    wavefront = Wavefront.from_openpmd(str(path))
    with h5py.File(path) as h5:
        if group not in h5:
            raise ValueError(f"{path.name} has no {group} group, so it is not {what}")
        g = h5[group]
        z_transmit = np.asarray(g["zTransmit"][...], dtype=float)
        params = np.asarray(g["wavefrontParams"][...], dtype=float)
        names = _names(g["wavefrontParams"])
        params_y = None
        if "wavefrontParamsY" in g:
            params_y = np.asarray(g["wavefrontParamsY"][...], dtype=float)
    if names != PARAM_NAMES:
        raise ValueError(f"{path.name} names its params columns {names}, and this reader "
                         f"knows {PARAM_NAMES}")
    n = len(z_transmit)
    if z_transmit.ndim != 1 or n > wavefront.nz or params.shape != (n, 25) or \
            (params_y is not None and params_y.shape != (n, 25)):
        raise ValueError(f"{path.name} holds {wavefront.nz} planes, zTransmit of shape "
                         f"{z_transmit.shape} and wavefrontParams of shape {params.shape}")
    return wavefront, z_transmit, params, params_y


def read_escaped(path):
    """
    The escaped-field file at ``path`` as an ``Escaped``: its planes as a Wavefront
    through ``Wavefront.from_openpmd``, with ``zTransmit`` and ``wavefrontParams`` of
    its ``lucifer/escaped`` group. A file without that group is refused, and so is one
    whose params columns are not ``PARAM_NAMES`` or whose arrays do not match the mesh's
    plane count. A run writes the file only with ``run.write_escaped``. Needs
    openpmd-beamphysics.
    """
    wavefront, z_transmit, params, params_y = _read(path, ESCAPED_GROUP, "an escaped-field file")
    if len(z_transmit) != wavefront.nz:
        raise ValueError(f"{pathlib.Path(path).name} holds {wavefront.nz} planes and "
                         f"{len(z_transmit)} rows of zTransmit")
    return Escaped(wavefront, z_transmit, params, params_y)


def read_pulse(path):
    """
    The pulse file at ``path`` as a ``Pulse``: its planes as a Wavefront through
    ``Wavefront.from_openpmd``, with ``zTransmit`` and ``wavefrontParams`` of its
    ``lucifer/pulse`` group, one row per banked plane, and the count of the window's
    planes before them. A file without that group is refused, and so is one whose params
    columns are not ``PARAM_NAMES`` or that has more rows than planes. Needs
    openpmd-beamphysics.
    """
    wavefront, z_transmit, params, params_y = _read(path, PULSE_GROUP, "a pulse file")
    return Pulse(wavefront, z_transmit, params, params_y, n_window=wavefront.nz - len(z_transmit))
