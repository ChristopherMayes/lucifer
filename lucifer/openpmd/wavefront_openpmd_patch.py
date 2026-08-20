"""openPMD EXT_Wavefront I/O for the openPMD-beamphysics Wavefront class -- a carried
patch, written to be grafted into beamphysics/wavefront/wavefront.py (which today has
only the Genesis4 pair, from_genesis4/write_genesis4). Until it lands upstream it is
importable standalone; the Bmad FEL harness (check_harmonics.py) exercises it against
the Fortran writer every run.

Layout (the standard document EXT_Wavefront.md, branch upcoming-2.0.0, is
authoritative; UPSTREAM-NOTES.md in this directory records the clarifications the
implementation needed):

  /                             openPMD, openPMDextension "Wavefront", basePath
                                "/data/%T/", meshesPath "meshes/", iterationEncoding
  /data/1/meshes/electricField  the mesh record; required attributes ON THE RECORD:
                                photonEnergy [J], temporalDomain 'time',
                                spatialDomain 'r', zCoordinate [m], plus the base
                                standard's geometry/axisLabels/gridSpacing/
                                gridGlobalOffset/gridUnitSI/unitDimension/timeOffset
  .../electricField/x, y        complex compound {r, i} datasets in V/m -- h5py maps
                                them to complex128 natively. Both transverse
                                polarizations in ONE file; z never written (paraxial).

The stored dataset order is (z, y, x) -- axisLabels declare it -- so axis 0 is the
slice axis in the file's C-order view while the class's (nx, ny, nz) convention is
recovered by the transpose below. One file per harmonic (the record name is fixed by
the extension, and photonEnergy is per record).
"""

from __future__ import annotations

import h5py
import numpy as np

H_PLANCK = 6.62607015e-34      # [J s]
C_LIGHT = 2.99792458e8

_MESH = "data/1/meshes/electricField"


def write_openpmd(wavefront, file, z_coordinate=0.0):
    """Write a Wavefront as an openPMD EXT_Wavefront file (module docstring layout).

    Parameters
    ----------
    wavefront : Wavefront
        Real-space wavefront; Ex and/or Ey in V/m, shape (nx, ny, nz).
    file : str or h5py.File
    z_coordinate : float
        Beamline position of the dump plane [m] (the extension's zCoordinate).
    """
    wf = wavefront
    if isinstance(file, h5py.File):
        return _write(wf, file, z_coordinate)
    with h5py.File(file, "w") as h5:
        return _write(wf, h5, z_coordinate)


def _write(wf, h5, z_coordinate):
    h5.attrs["openPMD"] = np.bytes_("2.0.0")
    h5.attrs["openPMDextension"] = np.bytes_("Wavefront")
    h5.attrs["basePath"] = np.bytes_("/data/%T/")
    h5.attrs["meshesPath"] = np.bytes_("meshes/")
    h5.attrs["iterationEncoding"] = np.bytes_("groupBased")
    h5.attrs["iterationFormat"] = np.bytes_("/data/%T/")

    it = h5.create_group("data/1")
    it.attrs["time"] = 0.0
    it.attrs["dt"] = 0.0
    it.attrs["timeUnitSI"] = 1.0

    m = h5.create_group(_MESH)
    some = wf.Ex if wf.Ex is not None else wf.Ey
    nx, ny, _ = some.shape
    m.attrs["geometry"] = np.bytes_("cartesian")
    m.attrs["axisLabels"] = np.array([b"z", b"y", b"x"])
    m.attrs["gridSpacing"] = np.array([wf.dz, wf.dy, wf.dx])
    m.attrs["gridGlobalOffset"] = np.array([0.0, -(ny - 1) * wf.dy / 2, -(nx - 1) * wf.dx / 2])
    m.attrs["gridUnitSI"] = 1.0
    m.attrs["unitDimension"] = np.array([1.0, 1.0, -3.0, -1.0, 0.0, 0.0, 0.0])
    m.attrs["timeOffset"] = 0.0
    m.attrs["photonEnergy"] = H_PLANCK * C_LIGHT / wf.wavelength
    m.attrs["temporalDomain"] = np.bytes_("time")
    m.attrs["spatialDomain"] = np.bytes_("r")
    m.attrs["zCoordinate"] = float(z_coordinate)

    for name, comp in (("x", wf.Ex), ("y", wf.Ey)):
        if comp is None:
            continue
        d = m.create_dataset(name, data=np.ascontiguousarray(comp.transpose(2, 1, 0)))
        d.attrs["unitSI"] = 1.0
        d.attrs["position"] = np.array([0.0, 0.0, 0.0])


def from_openpmd(file, wavefront_cls=None):
    """Read an openPMD EXT_Wavefront file into a Wavefront.

    Refuses, by name, what the class cannot represent: a frequency-domain or k-space
    field, and any missing required attribute. Returns wavefront_cls(...) if given,
    else a dict of the constructor arguments (so this module stays importable without
    the class for harness use).
    """
    with h5py.File(file, "r") as h5:
        if _MESH not in h5:
            raise ValueError(f"no {_MESH} mesh record: not an EXT_Wavefront file")
        m = h5[_MESH]

        for req in ("temporalDomain", "spatialDomain", "photonEnergy", "zCoordinate",
                    "axisLabels", "gridSpacing"):
            if req not in m.attrs:
                raise ValueError(f"required EXT_Wavefront attribute {req} missing")

        def scalar(a):
            v = m.attrs[a]
            return v[0] if np.ndim(v) else v

        def text(a):
            v = scalar(a)
            return v.decode() if isinstance(v, bytes) else str(v)

        if text("temporalDomain") != "time":
            raise ValueError(f"temporalDomain {text('temporalDomain')!r}: only the "
                             "time domain (field in V/m) is implemented")
        if text("spatialDomain") != "r":
            raise ValueError(f"spatialDomain {text('spatialDomain')!r}: only "
                             "cartesian r space is implemented")
        labels = [v.decode() if isinstance(v, bytes) else str(v) for v in m.attrs["axisLabels"]]
        if labels != ["z", "y", "x"]:
            raise ValueError(f"axisLabels {labels}: only the (z, y, x) stored order "
                             "is implemented")

        spacing = np.asarray(m.attrs["gridSpacing"], dtype=float)
        kwargs = dict(
            dz=float(spacing[0]), dy=float(spacing[1]), dx=float(spacing[2]),
            wavelength=H_PLANCK * C_LIGHT / float(scalar("photonEnergy")),
            Ex=None, Ey=None,
        )
        for name, key in (("x", "Ex"), ("y", "Ey")):
            if name in m:
                kwargs[key] = np.ascontiguousarray(m[name][:].transpose(2, 1, 0))

    if wavefront_cls is None:
        return kwargs
    return wavefront_cls(**kwargs)
