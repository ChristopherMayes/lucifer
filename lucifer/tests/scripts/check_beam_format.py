#!/usr/bin/env python3
"""
Checks for the particle dump format (fel-physics.md sec-import). Most claims here are a
round trip against the beam that was written. A round trip proves that this program agrees
with itself, which is not the same as proving the file says what it means: a convention
that is wrong and self-consistent survives every one of them, and one did (FINDINGS 7.64).
Section 8 therefore reads a dump with openPMD-beamphysics, which knows the standard and
nothing about this program.

1. The file round trip is exact. A run writes .beam.h5, a second run reads it back and
   writes it again, and every dataset of the two files must be bit-identical. This is the
   claim the format exists for: the file is a faithful container.

2. The state round trip. The file stores absolute momenta and a time where the tracker
   keeps packed (px, py) and a lag, so those pass through P/p0 and -beta*c*dt on the way
   out and back, while x and y are untouched. The state is read out here through the
   converter, in Genesis's chart, which is the only other view of the beam there is now:
   theta and gamma instead of a time and a momentum. Every column comes back exact on this
   configuration, so the levels are bounds on what the chart could cost rather than
   expected values: which coordinate pays is fixed, whether it pays is configuration
   dependent.

3. The phase survives. A dump carries the particle lag, and a reader restarts the
   reference phase phi0 at zero, so the writer folds phi0 into the lag. theta must
   therefore come back with NO offset at all, constant or otherwise. Without the fold the
   restarted beam sits at a different phase against the field, which is a real change of
   state and not a bookkeeping one: 2.1e-2 on the windowed-composition check.

4. Weights survive, which the Genesis format cannot do. The split-weight instrument makes
   every particle two coincident copies of w/3 and 2w/3, so the file's weight record is a
   real dataset of two distinct values. It must come back bit-identical. Converting that
   beam to a Genesis .par.h5 must be refused, since that format carries one
   current per slice and a read-back would silently return a uniform beam.

5. An empty slice survives. Every slice is a particlePatch, in window order, and an empty
   slice is a patch of no particles, so the patch list IS the window. A heavy migration run
   (sig_pz = 0.15) empties a slice; the restored beam must have the same per-slice counts,
   empty slice included, and the file must hold one patch per slice.

6. Refusals: a beam file that is not openPMD (the message names the
   converter), a file whose patch count disagrees with the window the deck states, and a
   file that carries no charge.

Usage: check_beam_format.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys

import h5py
import numpy as np

import beamio
import convert_genesis
from nml import to_groups

FAILED = False

LAMBDA0 = 1e-10
SAMPLE = 3
SPACING = SAMPLE * LAMBDA0

BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis.bmad"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 1024
  beam_init%bunch_charge = 4.803322970853e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -2.4e-9
  beam_init%grid(3)%x_max = 2.4e-9
  beam_init%sig_pz = {sig_pz}
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = 64
  grid_half_width = 2e-4
  beamlet_size = 8
  window_length = 4.8e-9
  n_wavelength = 3
  ran_seed = 777
  write_diag = T
{extra}&end
"""


def run(exe, wd, root, text, expect_fail=False):
    (wd / (root + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), root + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
    if expect_fail:
        return r.returncode, r.stdout
    if r.returncode != 0:
        print(f"FAIL: {root}.nml exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    return r.returncode, r.stdout


def check(label, value, level):
    global FAILED
    ok = value <= level
    FAILED = FAILED or not ok
    print(f"--- {label}: {value:.3e} (check {level:.0e}) {'ok' if ok else '** FAIL **'}")


def refused(label, code, out, phrase):
    global FAILED
    ok = code != 0 and phrase in out
    FAILED = FAILED or not ok
    print(f"--- refusal {label}: {'ok' if ok else '** FAIL **'}")


def datasets(fn):
    """Every record of an openPMD particle file, by path, constant ones expanded.

    A constant record is a group carrying value and shape rather than a dataset, and
    collecting only shaped datasets makes those invisible. timeOffset was one until it
    began carrying the slice placement, so this round trip compared everything in the
    file except the one record it had no eyes for.
    """
    out = {}

    def take(name, obj):
        if isinstance(obj, h5py.Dataset) and obj.shape:
            out[name] = obj[...]
        elif isinstance(obj, h5py.Group) and "value" in obj.attrs and "shape" in obj.attrs:
            n = int(np.ravel(obj.attrs["shape"])[0])
            out[name] = np.full(n, np.ravel(obj.attrs["value"])[0])

    with h5py.File(fn) as h5:
        h5.visititems(take)
    return out


def weights(fn):
    """The weight record, whether stored as a dataset or as a constant pseudo-dataset."""
    with h5py.File(fn) as h5:
        w = h5["data/00001/particles/electron/weight"]
        if isinstance(w, h5py.Dataset) and w.shape:
            return w[...], "dataset"
        return np.full(int(w.attrs["shape"][0]), float(w.attrs["value"][0])), "constant"


def as_genesis(wd, src, dst):
    """The same beam in Genesis's chart, through the converter: the second view of the
    state a check can compare against."""
    slices = beamio.read_slices(wd / src, LAMBDA0, SPACING)
    convert_genesis.write_genesis_par(wd / dst, slices, LAMBDA0, SPACING)
    return wd / dst


def par_field(fn, key, islice=1):
    with h5py.File(fn) as h5:
        return h5[f"slice{islice:06d}/{key}"][...]


def slice_counts(fn):
    return [sl["n"] for sl in beamio.read_slices(fn, LAMBDA0, SPACING)]


def patch_count(fn):
    with h5py.File(fn) as h5:
        g = h5["data/00001/particles/electron/particlePatches/numParticles"]
        return np.atleast_1d(g[...]).size


def restart(root, source, field, extra=""):
    """A load_only run that reads `source` and re-dumps its initial state. The driver
    requires a field beside a beam, so the source run's own field dump comes along."""
    return BASE.format(root=root, sig_pz="5.282703940115e-03",
                       extra=f'  beam_file = "{source}"\n  field_file = "{field}"\n'
                             f"  load_only = T\n  write_initial = T\n{extra}")


def main():
    global FAILED
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--pyrepo", help="openPMD-beamphysics checkout. "
                   "Default: $OPENPMD_BEAMPHYSICS, or the installed package.")
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(a.exe).resolve()

    # ------------------------------------------------------------------
    # A tracked beam, dumped, read back and dumped again.

    print("== the file round trip ==")
    run(exe, wd, "bf1", BASE.format(root="bf1", sig_pz="5.282703940115e-03", extra=""))
    run(exe, wd, "bfp", restart("bfp", "bf1-final.beam.h5", "bf1-final.wf.h5"))

    a1, a2 = datasets(wd / "bf1-final.beam.h5"), datasets(wd / "bfp-initial.beam.h5")

    # Three records say where the beam is rather than what it is, and a restart that
    # re-emits it at the start of the line rather than at the plane it was dumped from is
    # entitled to move all three: the element it sits in, the path length there, and the
    # arrival time. They were invisible to this check until constant records were
    # expanded, which is how timeOffset came to carry the slice placement unexamined.
    #
    # What must survive is the placement, which is timeOffset measured against its own
    # head. It survives to a few ulps of the absolute time rather than bit for bit, since
    # a double at 1.9e-7 s resolves 4e-23 s and a slice step here is 1e-18 s.
    WHERE = ["data/00001/particles/electron/" + k
             for k in ("timeOffset", "sPosition", "elementIndex")]
    same = sorted(a1) == sorted(a2) and \
           all(np.array_equal(a1[k], a2[k]) for k in a1 if k not in WHERE)
    check("file datasets bit-identical (write, read, write)", 0.0 if same else 1.0, 0.5)

    toff = "data/00001/particles/electron/timeOffset"
    pa, pb = a1[toff] - a1[toff].max(), a2[toff] - a2[toff].max()
    check("the slice placement survives the round trip",
          float(np.max(np.abs(pa - pb))) / (SPACING / 2.99792458e8), 1e-4)

    # ------------------------------------------------------------------
    # The field file read as an import. Three things a file written elsewhere may do
    # that this writer never does: state its samples in a unit other than V/m with
    # unitSI carrying the factor, state its grid with gridUnitSI carrying the factor,
    # and carry a y component shaped unlike x. The first two must read as the file
    # this writer would have written for the same field, and the third is refused
    # before the buffer is written past.

    print("== the field file as an import ==")

    def mesh(h5):
        """The one mesh record of a field file, whichever iteration holds it."""
        return h5["data"][sorted(h5["data"])[0]]["meshes/electricField"]

    def field_arrays(fn):
        """Every dataset of a field file, by path."""
        out = {}
        with h5py.File(fn) as h5:
            h5.visititems(lambda n, o: out.__setitem__(n, o[...])
                          if isinstance(o, h5py.Dataset) else None)
        return out

    def scaled_copy(src, dst, sample_factor, unit, grid_factor):
        shutil.copy(src, dst)
        with h5py.File(dst, "r+") as h5:
            m = mesh(h5)
            for comp in ("x", "y"):
                if comp in m:
                    m[comp][...] = m[comp][...] * sample_factor
                    m[comp].attrs["unitSI"] = np.array([unit])
            m.attrs["gridSpacing"] = np.ravel(m.attrs["gridSpacing"]) / grid_factor
            m.attrs["gridUnitSI"] = np.ravel(m.attrs["gridUnitSI"]) * grid_factor

    plain = field_arrays(wd / "bfp-initial.wf.h5")
    scaled_copy(wd / "bf1-final.wf.h5", wd / "unit_half.wf.h5", 2.0, 0.5, 1.0)
    run(exe, wd, "bfun", restart("bfun", "bf1-final.beam.h5", "unit_half.wf.h5"))
    got = field_arrays(wd / "bfun-initial.wf.h5")
    same = sorted(got) == sorted(plain) and all(np.array_equal(got[k], plain[k]) for k in plain)
    check("unitSI = 0.5 on doubled samples reads as the plain file (0 = yes)",
          0.0 if same else 1.0, 0.5)

    scaled_copy(wd / "bf1-final.wf.h5", wd / "grid_unit.wf.h5", 1.0, 1.0, 2.0)
    run(exe, wd, "bfgu", restart("bfgu", "bf1-final.beam.h5", "grid_unit.wf.h5"))
    got = field_arrays(wd / "bfgu-initial.wf.h5")
    with h5py.File(wd / "bfgu-initial.wf.h5") as h5:
        sp_got = np.ravel(mesh(h5).attrs["gridSpacing"])
    with h5py.File(wd / "bfp-initial.wf.h5") as h5:
        sp_plain = np.ravel(mesh(h5).attrs["gridSpacing"])
    same = sorted(got) == sorted(plain) and all(np.array_equal(got[k], plain[k]) for k in plain)
    check("gridUnitSI = 2 on a halved spacing reads as the plain file (0 = yes)",
          0.0 if (same and np.array_equal(sp_got, sp_plain)) else 1.0, 0.5)

    shutil.copy(wd / "bf1-final.wf.h5", wd / "ey_wide.wf.h5")
    with h5py.File(wd / "ey_wide.wf.h5", "r+") as h5:
        m = mesh(h5)
        x = m["x"]
        wide = np.concatenate([x[...], x[...]], axis=2)
        y = m.create_dataset("y", data=wide)
        for k, v in x.attrs.items():
            y.attrs[k] = v
    code, out = run(exe, wd, "bfey", restart("bfey", "bf1-final.beam.h5", "ey_wide.wf.h5"),
                    expect_fail=True)
    refused("a y component shaped unlike x", code, out, "COMPONENT y IS SHAPED")

    # The source filter on an imported field. Its angle is converted to a cutoff on the
    # grid, and the conversion once came from the grid the deck stated, which a deck that
    # imports its field need not state: the cutoff was then infinite and every mode
    # passed while the header named the angle asked for. The conversion comes from the
    # field's own spacing now, so a deck without the grid keys filters as one with them.
    nogrid = BASE.replace("  grid_n_pts = 64\n", "").replace("  grid_half_width = 2e-4\n", "")
    load = ('  beam_file = "bf1-final.beam.h5"\n  field_file = "bf1-final.wf.h5"\n'
            '  global%track_end = "UND##1"\n  source_filter = T\n')
    for root, deck, ang in (("bff6", nogrid, "1e-6"), ("bff5", nogrid, "1e-5"), ("bffg", BASE, "1e-5")):
        run(exe, wd, root, deck.replace("  source_filter = F\n", "").format(
            root=root, sig_pz="5.282703940115e-03",
            extra=load + f"  global%source_filter_angle = {ang}\n"))
    def x_field(fn):
        with h5py.File(fn) as h5:
            return mesh(h5)["x"][...]
    f6, f5, fg = (x_field(wd / f"{r}-final.wf.h5") for r in ("bff6", "bff5", "bffg"))
    scale = float(np.max(np.abs(f5)))
    check("the filter angle acts on an imported field with no grid in the deck (relative move)",
          1.0 / max(float(np.max(np.abs(f6 - f5))) / scale, 1e-300) * 1e-6, 1.0)
    check("and the deck with the grid keys filters as the deck without them",
          float(np.max(np.abs(fg - f5))) / scale, 1e-9)

    # ------------------------------------------------------------------
    print("== the state round trip, seen in the other chart ==")

    g0 = as_genesis(wd, "bf1-final.beam.h5", "bf1-conv.par.h5")
    g1 = as_genesis(wd, "bfp-initial.beam.h5", "bfp-conv.par.h5")

    for key, level in (("x", 1e-16), ("y", 1e-16), ("px", 1e-15), ("py", 1e-15),
                       ("gamma", 1e-15), ("current", 1e-16)):
        u, v = par_field(g0, key), par_field(g1, key)
        scale = np.max(np.abs(u)) or 1.0
        check(f"restart, {key}", float(np.max(np.abs(u - v))) / scale, level)

    # theta with the reference phase folded in: no offset, not even a constant one.
    u, v = par_field(g0, "theta"), par_field(g1, "theta")
    check("restart, theta absolute [rad]", float(np.max(np.abs(u - v))), 1e-9)
    check("restart, theta offset spread [rad]", float(np.ptp(u - v)), 1e-9)

    # ------------------------------------------------------------------
    print("== weights, which only openPMD carries ==")

    run(exe, wd, "bfw", BASE.format(root="bfw", sig_pz="5.282703940115e-03",
        extra="  split_weights = T\n"))
    run(exe, wd, "bfwr", restart("bfwr", "bfw-final.beam.h5", "bfw-final.wf.h5"))

    w0, kind0 = weights(wd / "bfw-final.beam.h5")
    w1, kind1 = weights(wd / "bfwr-initial.beam.h5")
    nd = np.unique(w0).size
    print(f"    written as a {kind0} of {w0.size} weights, {nd} distinct, "
          f"total {w0.sum():.6e} C")
    check("split weights are stored per particle",
          0.0 if (kind0 == "dataset" and nd == 2) else 1.0, 0.5)
    check("split weights bit-identical through the round trip",
          0.0 if np.array_equal(w0, w1) else 1.0, 0.5)

    r = subprocess.run([sys.executable,
                        str(pathlib.Path(__file__).resolve().parent / "convert_genesis.py"),
                        "to-genesis", "bfw-final.beam.h5", "bfw-conv.par.h5",
                        "--wavelength", str(LAMBDA0), "--sample", str(SAMPLE)],
                       cwd=wd, capture_output=True, text=True)
    refused("nonuniform weights to Genesis format", r.returncode, r.stdout,
            "GENESIS FORMAT CANNOT CARRY PER-PARTICLE WEIGHTS")

    # ------------------------------------------------------------------
    print("== an empty slice survives ==")

    run(exe, wd, "bfe", BASE.format(root="bfe", sig_pz="1.5e-01",
        extra="  migrate = T\n"))
    run(exe, wd, "bfer", restart("bfer", "bfe-final.beam.h5", "bfe-final.wf.h5"))

    c0 = slice_counts(wd / "bfe-final.beam.h5")
    c1 = slice_counts(wd / "bfer-initial.beam.h5")
    n_empty = sum(1 for c in c0 if c == 0)
    print(f"    {len(c0)} slices, {n_empty} empty, counts {c0}")
    check("the configuration actually empties a slice", 0.0 if n_empty > 0 else 1.0, 0.5)
    check("per-slice counts restored, empty slice included", 0.0 if c0 == c1 else 1.0, 0.5)
    check("one patch per slice, empty ones included",
          0.0 if patch_count(wd / "bfe-final.beam.h5") == len(c0) else 1.0, 0.5)

    # ------------------------------------------------------------------
    print("== refusals ==")

    # Not openPMD: a Genesis dump, which the message must name the converter for.
    convert_genesis.write_genesis_par(wd / "asgenesis.par.h5",
                                      beamio.read_slices(wd / "bf1-final.beam.h5",
                                                         LAMBDA0, SPACING),
                                      LAMBDA0, SPACING)
    code, out = run(exe, wd, "bfng", restart("bfng", "asgenesis.par.h5", "bf1-final.wf.h5"),
                    expect_fail=True)
    refused("a beam file that is not openPMD", code, out, "BEAM FILE IS NOT openPMD")

    # A window the file does not have: the deck states one slice fewer than the file holds.
    code, out = run(exe, wd, "bfnw",
                    restart("bfnw", "bf1-final.beam.h5", "bf1-final.wf.h5",
                            extra="  window_length = 4.5e-9\n"), expect_fail=True)
    refused("patch count against the deck's window", code, out, "PARTICLE PATCHES BUT THE DECK")

    # A file with no charge: every weight zeroed, which would track and radiate nothing.
    shutil.copy(wd / "bf1-final.beam.h5", wd / "noq.beam.h5")
    with h5py.File(wd / "noq.beam.h5", "r+") as h5:
        g = h5["data/00001/particles/electron"]
        w = g["weight"]
        if isinstance(w, h5py.Dataset) and w.shape:
            w[...] = 0.0
        else:
            w.attrs["value"] = np.array([0.0])
        g.attrs["totalCharge"] = np.array([0.0])
    code, out = run(exe, wd, "bfnq", restart("bfnq", "noq.beam.h5", "bf1-final.wf.h5"),
                    expect_fail=True)
    refused("a file that carries no charge", code, out, "FILE CARRIES NO CHARGE")

    # 8. The file read by a reader that shares none of this program's conventions.
    # particlePatches partitions the window and a general reader does not walk it, so the
    # placement of a slice in the bunch has to be in the particle data. It is: openPMD's
    # time and timeOffset sum to a particle's time, the lag rides time and the slice's own
    # offset rides timeOffset. Dropped, every slice lands on one and the bunch reads one
    # slice long (FINDINGS 7.64).
    #
    # The placement is checked on the per-patch offset and not on the spread of t. With
    # migration off a particle whose phase leaves its slice keeps going, so the lag is not
    # bounded by the slice width and the spread of t is a property of the physics. The
    # offset is a property of the format.
    #
    # The level is 1e-4 and it is a floor of the standard's own arithmetic, not of this
    # writer. timeOffset is an absolute time, and the ratio of a beamline's arrival time
    # to one slice spacing is 3e10 here, so a double leaves about five digits for the
    # slice structure: the ulp at t_ref is 6.5e-6 of a slice step. That is ample for
    # placing a picture and is not the exact partition, which stays particlePatches.

    print("== the file as openPMD, read by a reader that knows nothing of this program ==")
    convert_genesis._beamphysics(a.pyrepo)
    from beamphysics import ParticleGroup

    C_LIGHT = 2.99792458e8

    def patch_offsets(path):
        """Each patch's own time offset, which is what places its slice in the bunch."""
        with h5py.File(path) as h5:
            g = h5["data/00001/particles/electron"]
            n = g["particlePatches/numParticles"][()]
            off = g["particlePatches/numParticlesOffset"][()]
            to = g["timeOffset"]
            n_tot = g["id"].shape[0]
            tov = to[()] if isinstance(to, h5py.Dataset) else \
                  np.full(n_tot, float(np.ravel(to.attrs["value"])[0]))
        return np.array([tov[off[k]] for k in range(len(n)) if n[k] > 0]), n, n_tot

    tos, n_patch, n_file = patch_offsets(wd / "bf1-final.beam.h5")
    P = ParticleGroup(str(wd / "bf1-final.beam.h5"))

    # Consecutive slices are one spacing apart in light-travel distance, toward the head.
    step = -C_LIGHT * np.diff(tos)
    print(f"    slice step {step.mean():.6e} m over {len(tos)} patches, spacing {SPACING:.6e} m")
    check("consecutive slices sit one spacing apart in the file",
          float(np.max(np.abs(step - SPACING))) / SPACING, 1e-4)
    span = C_LIGHT * (tos.max() - tos.min())
    check("the bunch spans its window and not one slice",
          abs(span - (len(tos) - 1) * SPACING) / ((len(tos) - 1) * SPACING), 1e-4)
    check("the reader sees every particle", abs(len(P) - n_file), 0.5)
    with h5py.File(wd / "bf1-final.beam.h5") as h5:
        w = h5["data/00001/particles/electron/weight"]
        q_file = float(np.sum(w[()])) if isinstance(w, h5py.Dataset) \
                 else n_file * float(np.ravel(w.attrs["value"])[0])
    check("the reader sees the whole charge",
          abs(P.charge - q_file) / max(abs(q_file), 1e-30), 1e-12)

    # A frame cut to a range places the slices it carries at their true window positions.
    # Two runs at one comb, one whole and one cut, so the same frame is compared with and
    # without the range and the reference time is the same on both sides.

    frame_extra = "  dump_at_comb = T\n  comb_ds_save = 2.0\n"
    run(exe, wd, "bffull", BASE.format(root="bffull", sig_pz="8.8045e-5", extra=frame_extra))
    run(exe, wd, "bfrng", BASE.format(root="bfrng", sig_pz="8.8045e-5",
        extra=frame_extra + "  dump_slice_first = 4\n  dump_slice_last = 9\n"))
    full = sorted(wd.glob("bffull-[0-9]*.beam.h5"))
    rng = sorted(wd.glob("bfrng-[0-9]*.beam.h5"))
    if not full or not rng:
        check("a frame cut to a range places its slices where the whole window does", 1.0, 0.5)
    else:
        to_f = patch_offsets(full[-1])[0]
        to_r, nr, _ = patch_offsets(rng[-1])
        print(f"    range: {len(nr)} patches against {len(to_f)} in the whole window")
        check("a range carries the slices it asked for",
              abs(len(nr) - 6), 0.5)
        # Patch k of the range is window slice k + 3, so it must carry that slice's offset.
        worst = max(abs(C_LIGHT * (to_r[k] - to_f[k + 3])) / SPACING for k in range(len(to_r)))
        check("a frame cut to a range places its slices where the whole window does",
              worst, 1e-4)

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
