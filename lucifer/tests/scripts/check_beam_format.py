#!/usr/bin/env python3
"""
Checks for the particle dump formats (manual sec:import). Self-referenced: both formats
are this port's own writers, so every claim here is a round trip against the beam that
was written, never against another code.

1. THE openPMD FILE ROUND TRIP is exact. A run writes .beam.h5, a second run reads it
   back and writes it again, and every dataset of the two files must be bit-identical.
   This is the claim the format exists for: the file is a faithful container.

2. THE PACKED ARRAYS come back to a measured level, not bit for bit, and the reason is
   the chart. openPMD stores absolute momenta and a time, so a packed (px, py) makes a
   round trip through P/p0 and z through -beta*c*dt. Genesis format stores (theta,
   gamma), so it loses its ulp there instead: the two formats are mirror images and
   neither is the identity in floating point. Measured with both re-dumped Genesis
   format from the same beam, so the comparison is like for like:

     openPMD:  x, y and gamma exact;  px, py at 1e-16
     Genesis:  x, y and px, py exact;  gamma and theta at 1e-16

   Which coordinate pays is fixed by the chart. Whether it pays on a given beam is
   configuration-dependent, so the levels are bounds rather than expected values.

   theta additionally differs by a CONSTANT: a reader sets phi0 = 0 because the particle
   lag lives in z and the reference phase is free. That the offset is constant to
   rounding is the check, since a non-constant part would mean z did not survive.

3. WEIGHTS SURVIVE, which Genesis format cannot do. The split-weight instrument makes
   every particle two coincident copies of w/3 and 2w/3, so the file's weight record
   becomes a real dataset of two distinct values. It must come back bit-identical. The
   same beam written Genesis format must be REFUSED BY NAME, since that format carries
   one current per slice and a read-back would silently return a uniform beam.

4. AN EMPTY SLICE SURVIVES. hdf5_write_beam takes the species from p(1), which a
   zero-particle bunch does not have, so the writer skips empty slices and the window
   rides in root attributes instead of being implied by the bunch count. A heavy
   migration run (sig_pz = 0.15) empties a slice; the restored beam must have the same
   per-slice counts, empty slice included, and the file must hold one bunch fewer than
   the window with felSliceIndex skipping it.

5. REFUSALS, each by name: an unknown format token, and an openPMD particle file with no
   FEL window attributes (a bunch, not a sliced window, so it belongs on dist_file).

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

from nml import to_groups

FAILED = False

BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis.bmad"
  out_root = "{root}"
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
  grid_n_pts = 65
  grid_half_width = 2e-4
  nbins = 8
  window_length = 4.8e-9
  window_sample = 3
  ran_seed = 777
  write_diag = T
  beam_formats = {formats}
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
    """Every shaped dataset of an openPMD particle file, by path."""
    out = {}
    with h5py.File(fn) as h5:
        h5.visititems(lambda n, o: out.__setitem__(n, o[...])
                      if isinstance(o, h5py.Dataset) and o.shape else None)
    return out


def weights(fn):
    """The weight record, whether stored as a dataset or as a constant pseudo-dataset."""
    with h5py.File(fn) as h5:
        w = h5["data/00001/particles/electron/weight"]
        if isinstance(w, h5py.Dataset) and w.shape:
            return w[...], "dataset"
        return np.full(int(w.attrs["shape"][0]), float(w.attrs["value"][0])), "constant"


def slice_counts(fn):
    with h5py.File(fn) as h5:
        ns = int(h5["slicecount"][0])
        return [h5[f"slice{i:06d}/gamma"].shape[0] for i in range(1, ns + 1)]


def par_field(fn, key, islice=1):
    with h5py.File(fn) as h5:
        return h5[f"slice{islice:06d}/{key}"][...]


def restart(root, source, field, formats, extra=""):
    """A load_only run that reads `source` and re-dumps its initial state. The driver
    requires a field beside a beam, so the source run's own field dump comes along."""
    return BASE.format(root=root, sig_pz="5.282703940115e-03", formats=formats,
                       extra=f'  beam_file = "{source}"\n  field_file = "{field}"\n'
                             f"  load_only = T\n  write_initial = T\n{extra}")


def main():
    global FAILED
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(a.exe).resolve()

    # ------------------------------------------------------------------
    # A tracked beam, written in both formats, then read back from each.

    print("== openPMD round trip ==")
    run(exe, wd, "bf1", BASE.format(root="bf1", sig_pz="5.282703940115e-03",
        formats="'openpmd', 'genesis'", extra=""))
    run(exe, wd, "bfp", restart("bfp", "bf1-final.beam.h5", "bf1-final.fld.h5", "'openpmd', 'genesis'"))
    run(exe, wd, "bfg", restart("bfg", "bf1-final.par.h5", "bf1-final.fld.h5", "'openpmd', 'genesis'"))

    a1, a2 = datasets(wd / "bf1-final.beam.h5"), datasets(wd / "bfp-initial.beam.h5")
    same = sorted(a1) == sorted(a2) and all(np.array_equal(a1[k], a2[k]) for k in a1)
    check("file datasets bit-identical (write, read, write)", 0.0 if same else 1.0, 0.5)

    # The packed content, both restarts re-dumped Genesis format so the comparison is
    # like for like. x and y pass through untouched by either chart.
    for key, level in (("x", 1e-16), ("y", 1e-16), ("px", 1e-15), ("py", 1e-15), ("gamma", 1e-15)):
        u = par_field(wd / "bf1-final.par.h5", key)
        v = par_field(wd / "bfp-initial.par.h5", key)
        scale = np.max(np.abs(u)) or 1.0
        check(f"openPMD restart, {key}", float(np.max(np.abs(u - v))) / scale, level)

    for key, level in (("x", 1e-16), ("y", 1e-16), ("px", 1e-16), ("py", 1e-16), ("gamma", 1e-15)):
        u = par_field(wd / "bf1-final.par.h5", key)
        v = par_field(wd / "bfg-initial.par.h5", key)
        scale = np.max(np.abs(u)) or 1.0
        check(f"Genesis restart, {key}", float(np.max(np.abs(u - v))) / scale, level)

    # theta carries the reference phase, which a reader is free to restart at zero. The
    # offset must be constant: a spread in it would mean z did not survive.
    u = par_field(wd / "bf1-final.par.h5", "theta")
    v = par_field(wd / "bfp-initial.par.h5", "theta")
    check("openPMD restart, theta offset is constant [rad]", float(np.ptp(u - v)), 1e-9)

    # ------------------------------------------------------------------
    print("== weights, which only openPMD carries ==")

    run(exe, wd, "bfw", BASE.format(root="bfw", sig_pz="5.282703940115e-03",
        formats="'openpmd'", extra="  split_weights = T\n"))
    run(exe, wd, "bfwr", restart("bfwr", "bfw-final.beam.h5", "bfw-final.fld.h5", "'openpmd'"))

    w0, kind0 = weights(wd / "bfw-final.beam.h5")
    w1, kind1 = weights(wd / "bfwr-initial.beam.h5")
    nd = np.unique(w0).size
    print(f"    written as a {kind0} of {w0.size} weights, {nd} distinct, "
          f"total {w0.sum():.6e} C")
    check("split weights are stored per particle", 0.0 if (kind0 == "dataset" and nd == 2) else 1.0, 0.5)
    check("split weights bit-identical through the round trip",
          0.0 if np.array_equal(w0, w1) else 1.0, 0.5)

    code, out = run(exe, wd, "bfwg", BASE.format(root="bfwg", sig_pz="5.282703940115e-03",
                    formats="'genesis'", extra="  split_weights = T\n"), expect_fail=True)
    refused("nonuniform weights to Genesis format", code, out,
            "GENESIS FORMAT CANNOT CARRY PER-PARTICLE WEIGHTS")

    # ------------------------------------------------------------------
    print("== an empty slice survives ==")

    run(exe, wd, "bfe", BASE.format(root="bfe", sig_pz="1.5e-01", formats="'openpmd', 'genesis'",
        extra="  migrate = T\n"))
    run(exe, wd, "bfer", restart("bfer", "bfe-final.beam.h5", "bfe-final.fld.h5", "'genesis'"))

    c0, c1 = slice_counts(wd / "bfe-final.par.h5"), slice_counts(wd / "bfer-initial.par.h5")
    n_empty = sum(1 for c in c0 if c == 0)
    print(f"    {len(c0)} slices, {n_empty} empty, counts {c0}")
    check("the configuration actually empties a slice", 0.0 if n_empty > 0 else 1.0, 0.5)
    check("per-slice counts restored, empty slice included", 0.0 if c0 == c1 else 1.0, 0.5)

    with h5py.File(wd / "bfe-final.beam.h5") as h5:
        nb = len(h5["data"])
        ix = [int(i) for i in np.atleast_1d(h5.attrs["felSliceIndex"])]
        nslice = int(np.atleast_1d(h5.attrs["felNumSlice"])[0])
    expect = [i + 1 for i, c in enumerate(c0) if c > 0]
    check("one bunch per nonempty slice, felSliceIndex says which",
          0.0 if (nb == len(expect) and ix == expect and nslice == len(c0)) else 1.0, 0.5)

    # ------------------------------------------------------------------
    print("== refusals ==")

    code, out = run(exe, wd, "bfx", BASE.format(root="bfx", sig_pz="5.282703940115e-03",
                    formats="'openpmd', 'sdds'", extra=""), expect_fail=True)
    refused("unknown format token", code, out, "HAS AN UNKNOWN FORMAT")

    # An openPMD beam file that is a bunch rather than a sliced window: the same file
    # with its FEL attributes removed. It cannot be placed into slices, so it must be
    # refused with the import path named.
    shutil.copy(wd / "bf1-final.beam.h5", wd / "nofel.beam.h5")
    with h5py.File(wd / "nofel.beam.h5", "r+") as h5:
        for key in [k for k in h5.attrs if k.startswith("fel")]:
            del h5.attrs[key]
    code, out = run(exe, wd, "bfbr", restart("bfbr", "nofel.beam.h5", "bf1-final.fld.h5", "'openpmd'"), expect_fail=True)
    refused("openPMD file with no FEL window attributes", code, out,
            "CARRIES NO FEL WINDOW ATTRIBUTES")

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
