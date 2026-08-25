#!/usr/bin/env python3
"""
Read a Lucifer particle dump into per-slice arrays, in either format.

The two formats carry the same beam in different charts, so the checks should not each
learn both. Genesis .par.h5 stores (x, y, px, py, theta, gamma) per slice with one
current, and openPMD .beam.h5 stores position, momentum, time and a per-particle
weight, one bunch per nonempty slice, with the window in root attributes.

What comes back per slice: x, y [m], px, py [gamma*beta, Genesis's units, which
openPMD reaches through the file's own p0c], gamma, theta [rad], weight [C],
current [A], and n. Empty slices come back with n = 0 and empty arrays.

ONE CAVEAT, and it is the reason theta is documented rather than assumed. A reader
cannot recover the reference phase phi0 from either format: the Genesis path folds it
into theta, which a reader restarts at zero, and the openPMD path never stores it. So
theta here is defined with phi0 = 0, exactly as a restart would see it, and differs
from the writing run's theta by a constant. Every statistic these checks take of theta
is invariant under that constant (|b(h)| multiplies by a unit phasor), which is why the
formats can be compared at all.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np

C_LIGHT = 2.99792458e8
M_ELECTRON = 0.51099895069e6      # eV, Bmad's value
E_CHARGE = 1.602176634e-19


def _record(group, name, n):
    """A dataset, or the constant a pseudo-dataset stands for, as an array of length n."""
    obj = group[name]
    if isinstance(obj, h5py.Dataset) and obj.shape:
        return np.asarray(obj[...], dtype=float)
    return np.full(n, float(np.atleast_1d(obj.attrs["value"])[0]))


def _attr(h5, name):
    return np.atleast_1d(h5.attrs[name])[0]


def is_openpmd(path):
    """True for an openPMD file, by the root attribute the standard requires."""
    with h5py.File(path) as h5:
        return "openPMD" in h5.attrs


def read_slices(path):
    """Per-slice arrays from a Lucifer particle dump, either format."""
    path = pathlib.Path(path)
    return _read_openpmd(path) if is_openpmd(path) else _read_genesis(path)


def _read_genesis(path):
    out = []
    with h5py.File(path) as h5:
        nslice = int(h5["slicecount"][0])
        spacing = float(h5["slicespacing"][0])
        for i in range(1, nslice + 1):
            g = h5[f"slice{i:06d}"]
            n = g["gamma"].shape[0]
            current = float(g["current"][0])
            # The format has no per-particle weight: it carries one current, which the
            # reader divides out uniformly. That is the whole reason for .beam.h5.
            w = np.full(n, current * spacing / (C_LIGHT * n)) if n else np.zeros(0)
            out.append(dict(n=n, x=g["x"][...], y=g["y"][...],
                            px=g["px"][...], py=g["py"][...],
                            gamma=g["gamma"][...], theta=g["theta"][...],
                            weight=w, current=current))
    return out


def _read_openpmd(path):
    with h5py.File(path) as h5:
        nslice = int(_attr(h5, "felNumSlice"))
        spacing = float(_attr(h5, "felSliceSpacing"))
        ks = 2 * np.pi / float(_attr(h5, "felWavelength"))
        slice_ix = [int(i) for i in np.atleast_1d(h5.attrs["felSliceIndex"])]
        empty = dict(n=0, x=np.zeros(0), y=np.zeros(0), px=np.zeros(0), py=np.zeros(0),
                     gamma=np.zeros(0), theta=np.zeros(0), weight=np.zeros(0), current=0.0)
        out = [dict(empty) for _ in range(nslice)]

        for k, name in enumerate(sorted(h5["data"].keys())):
            g = h5[f"data/{name}/particles/electron"]
            n = int(np.atleast_1d(g.attrs["numParticles"])[0])
            p0c = _record(g, "totalMomentumOffset", n)
            pz = _record(g, "totalMomentum", n) / p0c        # written as pz * p0c
            p_mc = p0c * (1 + pz) / M_ELECTRON
            w = _record(g, "weight", n)
            # time is t - t_ref = -z/(beta c), so theta = ks*z/beta = -ks*c*time with
            # phi0 = 0. See the module docstring on the constant.
            theta = -ks * C_LIGHT * _record(g, "time", n)
            out[slice_ix[k] - 1] = dict(
                n=n, x=_record(g, "position/x", n), y=_record(g, "position/y", n),
                px=_record(g, "momentum/x", n) / M_ELECTRON,
                py=_record(g, "momentum/y", n) / M_ELECTRON,
                gamma=np.sqrt(p_mc ** 2 + 1), theta=theta, weight=w,
                current=float(w.sum()) * C_LIGHT / spacing)
    return out
