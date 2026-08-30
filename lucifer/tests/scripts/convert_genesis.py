#!/usr/bin/env python3
"""
Convert particle and field files between Genesis 1.3 v4 and openPMD.

The tracker speaks openPMD only. Genesis conversion lives here so the format knowledge
sits in one place, outside the physics.

Four directions. Fields both ways and Genesis particles inward are openPMD-beamphysics
calls: that library owns those conversions and this script does not repeat them. The
fourth, openPMD particles out to a Genesis .par.h5, has no upstream implementation and
is written here. It is what Genesis's &importbeam reads, so it is how Genesis restarts
from a state this tracker produced.

Particle files carry an FEL time window, which openPMD expresses as particlePatches: one
patch per slice, in window order, and an empty slice is a patch of no particles. The
Genesis format expresses the same thing as one group per slice with a current. Neither
carries the radiation wavelength for the beam, so a Genesis .par.h5 written here takes
it from the source file and one read here reports it for the caller to pass on.

Usage:
  convert_genesis.py to-openpmd <genesis file> <openpmd file> [--wavelength <m>]
  convert_genesis.py to-genesis <openpmd file> <genesis file> --wavelength <m>
                                                              [--sample <n>]
  convert_genesis.py round-trip <genesis file>      # convert out and back, compare

The kind of file (particles or field) is detected from its contents, so the same two
commands cover both. Exit 0 only on success.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import h5py
import numpy as np

C_LIGHT = 2.99792458e8
M_ELECTRON = 0.51099895069e6      # eV, Bmad's value
E_CHARGE = 1.602176634e-19


# ----------------------------------------------------------------------------------
# What kind of file is this


def is_openpmd(path):
    with h5py.File(path) as h5:
        return "openPMD" in h5.attrs


def is_particles(path):
    """Genesis particle dumps have slice groups with a gamma record; fields have none."""
    with h5py.File(path) as h5:
        if "openPMD" in h5.attrs:
            it = h5["data"][sorted(h5["data"].keys())[0]]
            return "particles" in it
        return "slicecount" in h5 and "slice000001/gamma" in h5


# ----------------------------------------------------------------------------------
# Genesis particles in


def read_genesis_par(path, pyrepo=None):
    """Per-slice ParticleGroups from a Genesis .par.h5, empty slices as None.

    The physics conversion is openPMD-beamphysics's: theta and gamma to a time and a
    momentum, the slice current to per-particle weights.
    """
    if pyrepo:
        sys.path.insert(0, str(pyrepo))
    from beamphysics import ParticleGroup

    with h5py.File(path) as h5:
        nslice = int(h5["slicecount"][0])
        scalars = dict(wavelength=float(h5["slicelength"][0]),
                       spacing=float(h5["slicespacing"][0]),
                       refposition=float(h5["refposition"][0]),
                       beamletsize=int(h5["beamletsize"][0]),
                       one4one=int(h5["one4one"][0]))
        counts = [h5[f"slice{i:06d}/gamma"].shape[0] for i in range(1, nslice + 1)]

    groups = []
    for i, n in enumerate(counts, start=1):
        # from_genesis4 numbers slices from 1 and refuses one with no particles, so an
        # empty slice is a zero-count patch here rather than a call.
        groups.append(None if n == 0 else ParticleGroup.from_genesis4(path, slices=[i]))
    return groups, scalars


# ----------------------------------------------------------------------------------
# openPMD particles out, in the layout Bmad's hdf5_read_beam reads


def _pmd_dataset(group, name, values, unit_si, unit_dim, local_name):
    d = group.create_dataset(name, data=np.asarray(values, dtype=float))
    d.attrs["unitSI"] = np.array([unit_si])
    d.attrs["unitDimension"] = np.array(unit_dim, dtype=float)
    d.attrs["unitSymbol"] = np.bytes_(b"")
    d.attrs["localName"] = np.bytes_(local_name.encode())
    return d


UNIT_M = (1.0, 0, 0, 0, 0, 0, 0)
UNIT_EV_C = (1.0, 1.0, -1.0, 0, 0, 0, 0)
UNIT_SEC = (0, 0, 1.0, 0, 0, 0, 0)
UNIT_1 = (0, 0, 0, 0, 0, 0, 0)

# unitSI per record, matching Bmad's pmd_unit_struct table: a reader multiplies by it.
SI_M = 1.0
SI_EV_C = E_CHARGE / C_LIGHT
SI_SEC = 1.0
SI_1 = 1.0


def write_openpmd_beam(path, groups, spacing):
    """One species record, one particlePatch per slice, in the Bmad reader's layout.

    spacing [m] is the slice spacing, needed because a ParticleGroup read from a Genesis
    dump carries a global z: slice i sits at (i-1)*spacing within the window. A patch
    holds its slice's own coordinate, the same local z the tracker keeps, so the slice
    offset comes back off here.
    """
    n_pat = [0 if g is None else g.n_particle for g in groups]
    n_tot = sum(n_pat)
    if n_tot == 0:
        raise ValueError("every slice is empty, so there is no beam to write")

    def cat(key, offset_per_slice=False):
        out = np.empty(n_tot)
        i = 0
        for k, g in enumerate(groups):
            if g is None:
                continue
            v = np.asarray(getattr(g, key), dtype=float)
            if offset_per_slice:
                v = v - k * spacing
            out[i:i + v.size] = v
            i += v.size
        return out

    x, y = cat("x"), cat("y")
    px, py, pz_tot = cat("px"), cat("py"), cat("pz")     # eV/c, absolute
    w = cat("weight")

    # The longitudinal coordinate. A ParticleGroup from a Genesis dump puts it in z, with
    # t the same for every particle: the dump is one plane. Bmad's beam file is the other
    # way round, position/z zero and the coordinate in time, and its reader rebuilds
    # vec(5) = -beta*c*time. Since a Genesis z is theta*lambda/(2 pi), time = -z/c makes
    # vec(5) come out as beta*theta/ks, which is exactly what reading the Genesis file
    # directly would have given.
    t = -cat("z", offset_per_slice=True) / C_LIGHT
    t_ref = 0.0

    with h5py.File(path, "w") as h5:
        h5.attrs["dataType"] = np.bytes_(b"openPMD")
        h5.attrs["openPMD"] = np.bytes_(b"2.0.0")
        h5.attrs["openPMDextension"] = np.bytes_(b"BeamPhysics;SpeciesType")
        h5.attrs["basePath"] = np.bytes_(b"/data/%T/")
        h5.attrs["particlesPath"] = np.bytes_(b"particles/")
        h5.attrs["software"] = np.bytes_(b"convert_genesis.py")
        h5.attrs["softwareVersion"] = np.bytes_(b"1.0")

        sp = h5.create_group("data/1/particles/electron")
        sp.attrs["speciesType"] = np.bytes_(b"electron")
        sp.attrs["numParticles"] = np.array([n_tot])
        sp.attrs["totalCharge"] = np.array([float(w.sum())])
        sp.attrs["chargeLive"] = np.array([float(w.sum())])
        sp.attrs["chargeUnitSI"] = np.array([1.0])

        pos = sp.create_group("position")
        _pmd_dataset(pos, "x", x, 1.0, UNIT_M, "x")
        _pmd_dataset(pos, "y", y, 1.0, UNIT_M, "y")
        _pmd_dataset(pos, "z", np.zeros(n_tot), 1.0, UNIT_M, "z")

        mom = sp.create_group("momentum")
        _pmd_dataset(mom, "x", px, SI_EV_C, UNIT_EV_C, "px * p0c")
        _pmd_dataset(mom, "y", py, SI_EV_C, UNIT_EV_C, "py * p0c")
        _pmd_dataset(mom, "z", pz_tot, SI_EV_C, UNIT_EV_C, "ps * p0c")

        # No totalMomentum and no totalMomentumOffset. Without them Bmad's reader takes
        # the momentum from the three absolute components and references it to the
        # lattice's p0c, which is the run's one reference. Writing a reference here would
        # invent a second one, and a beam's mean momentum is not the lattice's: doing that
        # shifted gamma by 5.1e-4 in the first version of this converter.

        _pmd_dataset(sp, "time", t - t_ref, 1.0, UNIT_SEC, "t - t_ref")
        _pmd_dataset(sp, "timeOffset", np.full(n_tot, t_ref), 1.0, UNIT_SEC, "t_ref")
        _pmd_dataset(sp, "weight", w, 1.0, UNIT_1, "macro-charge")
        _pmd_dataset(sp, "sPosition", np.zeros(n_tot), 1.0, UNIT_M, "s")

        st = sp.create_dataset("particleStatus", data=np.ones(n_tot, dtype=np.int32))
        st.attrs["unitSI"] = np.array([1.0])
        st.attrs["unitDimension"] = np.array(UNIT_1, dtype=float)
        st.attrs["unitSymbol"] = np.bytes_(b"")

        # Bmad's reader takes these through its openPMD dataset reader, which expects the
        # unit attributes every openPMD record carries.
        pp = sp.create_group("particlePatches")
        off = np.concatenate([[0], np.cumsum(n_pat)[:-1]]).astype(np.int32)
        for name_, vals in (("numParticles", np.asarray(n_pat, dtype=np.int32)),
                            ("numParticlesOffset", off)):
            d = pp.create_dataset(name_, data=vals)
            d.attrs["unitSI"] = np.array([1.0])
            d.attrs["unitDimension"] = np.array(UNIT_1, dtype=float)
            d.attrs["unitSymbol"] = np.bytes_(b"")
            d.attrs["localName"] = np.bytes_(name_.encode())
        lo, hi = pp.create_group("offset"), pp.create_group("extent")
        for comp, vals in (("x", x), ("y", y), ("z", np.zeros(n_tot))):
            b_lo, b_ex = [], []
            for k, n in enumerate(n_pat):
                seg = vals[off[k]:off[k] + n]
                b_lo.append(float(seg.min()) if n else 0.0)
                b_ex.append(float(seg.max() - seg.min()) if n else 0.0)
            _pmd_dataset(lo, comp, b_lo, 1.0, UNIT_M, comp)
            _pmd_dataset(hi, comp, b_ex, 1.0, UNIT_M, comp)


# ----------------------------------------------------------------------------------
# Genesis particles out, the direction openPMD-beamphysics does not have


def write_genesis_par(path, slices, wavelength, spacing, refposition=0.0,
                      beamletsize=1, one4one=None):
    """A Genesis .par.h5 from per-slice arrays, the format &importbeam reads.

    slices is a list of dicts with x, y, px, py, gamma, theta and weight per slice, the
    shape beamio.read_slices returns. one4one defaults to what the weights say: the flag
    asserts that every macroparticle carries one electron.

    A beam whose weights differ within a slice is refused by name. This format stores one
    current per slice, so a read-back would return a uniform beam, and per-particle weights
    are the whole reason the tracker's own format is openPMD. Weights differing between
    slices are fine: that is what a current profile is.
    """
    for i, s in enumerate(slices, start=1):
        if s["n"] < 2:
            continue
        w = np.asarray(s["weight"], dtype=float)
        if w.max() - w.min() <= 1e-12 * w.max():
            continue
        raise ValueError(
            f"GENESIS FORMAT CANNOT CARRY PER-PARTICLE WEIGHTS, AND THIS BEAM HAS "
            f"NONUNIFORM WEIGHTS.\n"
            f"  slice {i} spreads {w.min():.4e} to {w.max():.4e} C over {s['n']} "
            f"particles, total {w.sum():.4e} C.\n"
            f"  The format stores one current per slice, so a read-back would return a "
            f"uniform beam.")

    if one4one is None:
        w = np.concatenate([s["weight"] for s in slices if s["n"]])
        one4one = bool(np.all(np.abs(w - E_CHARGE) <= 1e-9 * E_CHARGE))

    with h5py.File(path, "w") as h5:
        h5.create_dataset("slicelength", data=np.array([wavelength]))
        h5.create_dataset("slicespacing", data=np.array([spacing]))
        h5.create_dataset("refposition", data=np.array([refposition]))
        h5.create_dataset("beamletsize", data=np.array([beamletsize], dtype=np.int32))
        h5.create_dataset("slicecount", data=np.array([len(slices)], dtype=np.int32))
        h5.create_dataset("one4one", data=np.array([1 if one4one else 0], dtype=np.int32))
        for i, s in enumerate(slices, start=1):
            g = h5.create_group(f"slice{i:06d}")
            n = s["n"]
            cur = float(s["weight"].sum()) * C_LIGHT / spacing if n else 0.0
            g.create_dataset("current", data=np.array([cur]))
            for key in ("x", "y", "px", "py", "gamma", "theta"):
                g.create_dataset(key, data=np.asarray(s[key], dtype=float))


# ----------------------------------------------------------------------------------
# Fields, both ways, entirely upstream


def convert_field(src, dst, pyrepo=None):
    if pyrepo:
        sys.path.insert(0, str(pyrepo))
    from beamphysics.wavefront import Wavefront

    if is_openpmd(src):
        Wavefront.from_openpmd(src).write_genesis4(str(dst))
    else:
        Wavefront.from_genesis4(str(src)).write_openpmd(str(dst))


# ----------------------------------------------------------------------------------

DEFAULT_PYREPO = str(pathlib.Path.home() / "Code/GitHub/bmad-fel/openPMD-beamphysics")


def to_openpmd(src, dst, pyrepo=DEFAULT_PYREPO):
    """Convert one Genesis file, particles or field, to openPMD. The importable form of
    the to-openpmd command, for the check scripts that run Genesis themselves."""
    src, dst = pathlib.Path(src), pathlib.Path(dst)
    if is_particles(src):
        groups, scalars = read_genesis_par(src, pyrepo)
        write_openpmd_beam(dst, groups, scalars["spacing"])
    else:
        convert_field(src, dst, pyrepo)
    return dst


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("command", choices=("to-openpmd", "to-genesis", "round-trip"))
    p.add_argument("src")
    p.add_argument("dst", nargs="?")
    p.add_argument("--wavelength", type=float, default=0.0)
    p.add_argument("--sample", type=int, default=1)
    p.add_argument("--pyrepo", default=DEFAULT_PYREPO)
    a = p.parse_args()
    src = pathlib.Path(a.src)

    if a.command == "round-trip":
        return round_trip(src, a)

    dst = pathlib.Path(a.dst)
    if not is_particles(src):
        convert_field(src, dst, a.pyrepo)
        print(f"field: {src.name} -> {dst.name}")
        return 0

    if a.command == "to-openpmd":
        groups, scalars = read_genesis_par(src, a.pyrepo)
        write_openpmd_beam(dst, groups, scalars["spacing"])
        print(f"particles: {src.name} -> {dst.name}, {len(groups)} patches")
        return 0

    # to-genesis
    if a.wavelength <= 0:
        print("FAIL: --wavelength is required to write a Genesis particle file: the "
              "format states it and an openPMD beam does not carry it.")
        return 1
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import beamio
    slices = beamio.read_slices(src, a.wavelength, a.sample * a.wavelength)
    try:
        write_genesis_par(dst, slices, a.wavelength, a.sample * a.wavelength)
    except ValueError as exc:
        # A refusal, not a crash: the message names the slice and the spread.
        print(f"FAIL: {exc}")
        return 1
    print(f"particles: {src.name} -> {dst.name}, {len(slices)} slices")
    return 0


def round_trip(src, a):
    """Convert out and back, and report how far the numbers moved."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        if not is_particles(src):
            convert_field(src, d / "mid.wf.h5", a.pyrepo)
            convert_field(d / "mid.wf.h5", d / "back.fld.h5", a.pyrepo)
            worst = 0.0
            with h5py.File(src) as A, h5py.File(d / "back.fld.h5") as B:
                names, back = [], []
                A.visit(lambda n: names.append(n) if isinstance(A[n], h5py.Dataset) else None)
                B.visit(lambda n: back.append(n) if isinstance(B[n], h5py.Dataset) else None)
                # A Genesis-written file carries blocks of its own beside the field: a
                # Meta/Version group and precomputed intensity projections. They are not
                # part of the format's field content, so they are named here and skipped
                # rather than counted as a difference.
                extra = sorted(set(names) - set(back))
                if extra:
                    print(f"  outside the field format, not compared: {', '.join(extra)}")
                for n in sorted(set(names) & set(back)):
                    u, v = A[n][...], B[n][...]
                    scale = np.max(np.abs(u)) or 1.0
                    worst = max(worst, float(np.max(np.abs(u - v))) / scale)
                    print(f"  {n:34s} rel {float(np.max(np.abs(u - v))) / scale:.3e}")
            print(f"field round trip worst relative difference: {worst:.3e}")
            return 0 if worst < 1e-12 else 1

        with h5py.File(src) as h5:
            lam = float(h5["slicelength"][0])
            spacing = float(h5["slicespacing"][0])
        groups, scalars = read_genesis_par(src, a.pyrepo)
        write_openpmd_beam(d / "mid.beam.h5", groups, spacing)
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
        import beamio
        write_genesis_par(d / "back.par.h5", beamio.read_slices(d / "mid.beam.h5", lam, spacing),
                          lam, spacing, scalars["refposition"], scalars["beamletsize"])
        worst = 0.0
        with h5py.File(src) as A, h5py.File(d / "back.par.h5") as B:
            ns = int(A["slicecount"][0])
            for i in range(1, ns + 1):
                for key in ("x", "y", "px", "py", "gamma", "theta", "current"):
                    u, v = A[f"slice{i:06d}/{key}"][...], B[f"slice{i:06d}/{key}"][...]
                    if u.size == 0:
                        continue
                    scale = np.max(np.abs(u)) or 1.0
                    rel = float(np.max(np.abs(u - v))) / scale
                    worst = max(worst, rel)
        print(f"particle round trip worst relative difference: {worst:.3e}")
        return 0 if worst < 1e-12 else 1


if __name__ == "__main__":
    sys.exit(main())
