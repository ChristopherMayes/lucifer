"""
Read and write the openPMD-beamphysics bunch files Bmad's beam_init%position_file reads.

A bunch file is one iteration, one species, every particle against the file's one time
reference: what Bmad's hdf5_write_beam writes and hdf5_read_beam reads. The tracker's own
beam dumps are a different thing, one particle patch per slice, and are not handled here.

Writing starts from a Bmad-written template, so every attribute Bmad expects is already
there. The per-particle records are replaced and the constant records are re-shaped.
"""

import shutil

import h5py
import numpy as np

M_ELECTRON = 0.51099895e6      # eV, Bmad's m_electron
C_LIGHT = 2.99792458e8
EV_C = 1.602176634e-19 / C_LIGHT   # kg m/s per eV/c, the momentum records' unitSI

RECORDS = ("position/x", "position/y", "momentum/x", "momentum/y", "momentum/z",
           "totalMomentum", "time", "weight")


def species_group(f):
    """The one species of the one iteration."""
    d = f["data"]
    it = sorted(d)
    if len(it) != 1:
        raise ValueError(f"{f.filename}: {len(it)} iterations; a bunch file has one")
    parts = d[it[0]]["particles"]
    sp = list(parts)
    if len(sp) != 1:
        raise ValueError(f"{f.filename}: {len(sp)} species; a bunch file has one")
    return parts[sp[0]]


def record(g, name):
    """A per-particle record as an array, constant records expanded."""
    o = g[name]
    n = int(np.atleast_1d(g.attrs["numParticles"])[0])
    if isinstance(o, h5py.Dataset):
        v = o[...]
    else:
        v = np.full(n, float(np.atleast_1d(o.attrs["value"])[0]))
    return v * float(np.atleast_1d(o.attrs.get("unitSI", 1.0))[0])


def read_bunch(path):
    """Every record in SI: x, y [m]; px, py, pz [kg m/s, EV_C per eV/c]; t [s, against the
    file's reference]; w [C]. p0c is the reference momentum in eV/c."""
    with h5py.File(path) as f:
        g = species_group(f)
        out = dict(x=record(g, "position/x"), y=record(g, "position/y"),
                   px=record(g, "momentum/x"), py=record(g, "momentum/y"),
                   pz=record(g, "momentum/z"), t=record(g, "time"), w=record(g, "weight"))
        out["p0c"] = float(np.atleast_1d(g["totalMomentumOffset"].attrs["value"])[0])
    return out


def write_bunch(template, path, b):
    """
    Write b (the dict read_bunch returns, any length) as a bunch file shaped like template.

    Every per-particle record becomes a dataset. Every constant record keeps its value and
    takes the new count in its shape attribute. The species attributes follow the weights.
    """
    shutil.copy(template, path)
    n = len(b["x"])
    with h5py.File(path, "r+") as f:
        g = species_group(f)
        values = {"position/x": b["x"], "position/y": b["y"], "momentum/x": b["px"],
                  "momentum/y": b["py"], "momentum/z": b["pz"], "time": b["t"], "weight": b["w"],
                  "totalMomentum": np.sqrt(b["px"] ** 2 + b["py"] ** 2 + b["pz"] ** 2)}
        for name in RECORDS:
            old = g[name]
            attrs = {k: old.attrs[k] for k in old.attrs if k not in ("shape", "value", "minValue", "maxValue")}
            del g[name]
            v = np.asarray(values[name], dtype=np.float64) / float(np.atleast_1d(attrs.get("unitSI", 1.0))[0])
            ds = g.create_dataset(name, data=v)
            for k, a in attrs.items():
                ds.attrs[k] = a
            ds.attrs["minValue"] = np.array([v.min()])
            ds.attrs["maxValue"] = np.array([v.max()])

        def reshape(name, o):
            if isinstance(o, h5py.Group) and "value" in o.attrs:
                o.attrs["shape"] = np.array([n], dtype=np.int32)
        g.visititems(reshape)
        g.attrs["numParticles"] = np.array([n], dtype=np.int32)
        g.attrs["totalCharge"] = np.array([float(b["w"].sum())])
        g.attrs["chargeLive"] = np.array([float(b["w"].sum())])


def arrival(b):
    """The loader's longitudinal coordinate s = c (t - t_ref) [m], the import's tau chart."""
    return C_LIGHT * b["t"]


def slice_index(s, spacing):
    """The keep loader's binning: ceiling(extent / spacing) slices from the earliest particle."""
    smin = s.min()
    ttotal = s.max() - smin
    nslice = max(1, int(np.ceil(ttotal / spacing - 1e-9)))
    isl = np.floor((s - smin) / spacing).astype(int) + 1
    isl[(isl == nslice + 1) & (s - smin <= nslice * spacing)] = nslice
    return isl, nslice, smin


def copies(b, n_copy, wavelength, spacing, jitter=0.0, shift=None, rng=None):
    """
    Every particle as n_copy copies at weight / n_copy, phases spread over one wavelength
    and kept inside the particle's own slice, so the bunch is quiet below n_copy the way
    the loader's own beamlets are, but made outside it. jitter perturbs the transverse
    coordinates of every copy by that many metres, so no two share them exactly. shift
    moves every copy singly by a uniform fraction of a wavelength, which is the mutation
    that breaks the cancellation.
    """
    rng = rng or np.random.default_rng(1)
    s = arrival(b)
    isl, nslice, smin = slice_index(s, spacing)
    n = len(s)
    out = {k: np.repeat(b[k], n_copy) for k in ("x", "y", "px", "py", "pz", "w")}
    out["w"] = out["w"] / n_copy
    out["p0c"] = b["p0c"]
    m = np.tile(np.arange(n_copy), n)
    sc = np.repeat(s, n_copy) + m * wavelength / n_copy
    if shift is not None:
        sc = sc + rng.uniform(0, 1, sc.size) * wavelength * shift
    lo = smin + (np.repeat(isl, n_copy) - 1) * spacing
    hi = lo + spacing
    over = sc >= hi
    sc[over] -= wavelength * np.ceil((sc[over] - hi[over]) / wavelength + 1e-12)
    sc = np.minimum(sc, np.nextafter(hi, -np.inf))
    out["t"] = sc / C_LIGHT
    if jitter:
        out["x"] = out["x"] + rng.uniform(-jitter, jitter, out["x"].size)
        out["y"] = out["y"] + rng.uniform(-jitter, jitter, out["y"].size)
    return out
