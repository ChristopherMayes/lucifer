"""
Read the escaped-field bank a run writes beside its final dump.

``<out_root>-escaped.wf.h5`` holds the light slippage carried out of the window, one
plane per slice as it was when it left, and ``<out_root>-pulse.wf.h5`` the whole pulse
at the exit, the window's slices and the escaped ones propagated to the exit plane.
Both are openPMD EXT_Wavefront files, so ``Wavefront.from_openpmd`` of
openpmd-beamphysics reads either whole, with both polarizations as its Ex and Ey. The
escaped file carries beside its mesh what openPMD has no name for, in the group
``lucifer/escaped``: ``zTransmit``, the beamline position at which each plane left the
window, and ``wavefrontParams``, each plane's moments as the program computed them
then, one row per plane in the mesh's order. ``read_escaped`` returns the Wavefront
with those two arrays, and the y component's params where the run carried one.

No one beamline position holds the escaped planes, so the file's ``zCoordinate``, and
the Wavefront's ``s_position``, is NaN. Each plane's own is ``zTransmit``. The pulse
file's is the exit plane.
"""

from __future__ import annotations

import dataclasses
import pathlib

import h5py
import numpy as np

__all__ = ["Escaped", "read_escaped"]

ESCAPED_GROUP = "lucifer/escaped"

# The 25 columns of wavefrontParams: the program's wavefront_params_struct flattened.
# The centroid over (x, theta_x, y, theta_y) [m, rad], the 4x4 second-moment matrix over
# the same axes column by column, then the energy [J], the power [W], the on-axis
# intensity [W/m^2] and the two emittances [m rad]. The file states them too, as the
# dataset's columns attribute, and read_escaped holds the two to agree.
PARAM_NAMES = (
    "centroid_x", "centroid_theta_x", "centroid_y", "centroid_theta_y",
    "sigma_11", "sigma_21", "sigma_31", "sigma_41", "sigma_12", "sigma_22", "sigma_32", "sigma_42",
    "sigma_13", "sigma_23", "sigma_33", "sigma_43", "sigma_14", "sigma_24", "sigma_34", "sigma_44",
    "energy", "power", "on_axis_intensity", "emit_x", "emit_y",
)


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

    def __repr__(self):
        n = len(self.z_transmit)
        e = self.params[:, PARAM_NAMES.index("energy")].sum()
        if self.params_y is not None:
            e += self.params_y[:, PARAM_NAMES.index("energy")].sum()
        span = f"{self.z_transmit.min():.4g} to {self.z_transmit.max():.4g} m" if n else "none"
        return f"Escaped({n} planes, left at {span}, {e:.4g} J)"


def _names(dataset):
    return tuple(v.decode() if isinstance(v, bytes) else str(v) for v in dataset.attrs["columns"])


def read_escaped(path):
    """
    The escaped-field file at ``path`` as an ``Escaped``: its planes as a Wavefront
    through ``Wavefront.from_openpmd``, with ``zTransmit`` and ``wavefrontParams`` of
    its ``lucifer/escaped`` group. A file without that group is refused, and so is one
    whose params columns are not ``PARAM_NAMES`` or whose arrays do not match the mesh's
    plane count. Needs openpmd-beamphysics.
    """
    from beamphysics.wavefront import Wavefront

    path = pathlib.Path(path)
    wavefront = Wavefront.from_openpmd(str(path))
    with h5py.File(path) as h5:
        if ESCAPED_GROUP not in h5:
            raise ValueError(f"{path.name} has no {ESCAPED_GROUP} group, so it is not an "
                             f"escaped-field file")
        g = h5[ESCAPED_GROUP]
        z_transmit = np.asarray(g["zTransmit"][...], dtype=float)
        params = np.asarray(g["wavefrontParams"][...], dtype=float)
        names = _names(g["wavefrontParams"])
        params_y = None
        if "wavefrontParamsY" in g:
            params_y = np.asarray(g["wavefrontParamsY"][...], dtype=float)
    if names != PARAM_NAMES:
        raise ValueError(f"{path.name} names its params columns {names}, and this reader "
                         f"knows {PARAM_NAMES}")
    n = wavefront.nz
    if z_transmit.shape != (n,) or params.shape != (n, 25) or \
            (params_y is not None and params_y.shape != (n, 25)):
        raise ValueError(f"{path.name} holds {n} planes, zTransmit of shape {z_transmit.shape} "
                         f"and wavefrontParams of shape {params.shape}")
    return Escaped(wavefront, z_transmit, params, params_y)
