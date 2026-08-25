#!/usr/bin/env python3
"""
Checks for slice migration under weights (deliverable 7, self-referenced per FINDINGS
6.9: Genesis migrates only under one4one, so weighted migration has no reference).

Three checks:

1. Conservation. A beam with artificially large energy spread (gen_delgam = 60) tracked
   through the full line migrates heavily and bleeds charge off the window ends. At
   every diagnostic record, in-window charge (sum of the per-slice current column times
   spacing/c) plus the charge dropped so far (per-event log lines, full precision) must
   equal the initial charge to 1e-10 relative. The run must actually bite: >10000 moves
   and nonzero drops are asserted.

2. Phase continuity. The same run carries migrate_check = T: at every migration the
   whole-beam weighted phasor must satisfy S_before = S_after + S_dropped -- every
   mover's phase shifts by an exact multiple of 2*pi*sample, every drop removes exactly
   its own term. The reported worst deviation must be < 1e-10 (measured: rounding,
   ~7e-15).

3. No-op. A frozen-phase configuration (delgam ~ 0, negligible emittance, dark, no
   noise) run with migrate = T must report zero moves and reproduce the migrate = F run
   bit for bit (diag byte-equal, dumps dataset-equal): migration must not fire
   spuriously and its inactive presence must change nothing.

Usage: check_migration.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis.bmad and aramis_1seg.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

import beamio

from nml import to_groups

C_LIGHT = 2.99792458e8

BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = {npart}
  beam_init%bunch_charge = {q}
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -{half}
  beam_init%grid(3)%x_max = {half}
  beam_init%sig_pz = {sig_pz}
  beam_init%a_norm_emit = {emit}
  beam_init%b_norm_emit = {emit}
  seed_power = 0
  grid_n_pts = 65
  grid_half_width = 2e-4
  nbins = 8
  window_length = {slen}
  window_sample = 3
  ran_seed = 777
  migrate = {mig}
  migrate_check = T
  write_diag = T
&end
"""


def run(exe, wd, nml_name, text):
    (wd / nml_name).write_text(to_groups(text))
    r = subprocess.run([str(exe), nml_name], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {nml_name} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    return r.stdout


def read_migration_file(wd, root):
    """Per-event migration rows and the summary from <out_root>.migration.txt. Program
    stdout is for humans and is deliberately not parsed (manual sec:program).
    Returns (moved, worst_bunching_deviation, [(s, charge_dropped), ...])."""
    moved = 0
    bdev = 0.0
    drops = []
    for line in (pathlib.Path(wd) / (root + ".migration.txt")).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if f[0] == "moved":
            moved = int(f[1])
        elif f[0] == "worst_bunching_deviation":
            bdev = float(f[1])
        elif f[0] == "charge_dropped_total":
            pass
        else:
            drops.append((float(f[0]), float(f[2])))
    return moved, bdev, drops


def in_window_charge(diag_file):
    d = np.loadtxt(diag_file)
    ns = int(d[:, 1].max())
    d = d.reshape(-1, ns, d.shape[1])
    return d[:, 0, 0], d[:, :, 10].sum(axis=1) * 3e-10 / C_LIGHT


def dumps_equal(fa, fb):
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names = []
        a.visit(lambda n: names.append(n) if isinstance(a[n], h5py.Dataset) else None)
        return all(np.array_equal(a[n][...], b[n][...]) for n in names)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    exe = pathlib.Path(a.exe).resolve()
    ok = True

    # 1 + 2: conservation and phase continuity under heavy migration.
    out = run(exe, wd, "migc.nml", BASE.format(lat="aramis.bmad", root="migc",
              sig_pz="5.282703940115e-03", emit="4e-7", npart=1024, slen="4.8e-9", q="4.803322970853e-14", half="2.4e-9", mig="T"))
    moved, bdev, drops = read_migration_file(wd, "migc")
    z, q_win = in_window_charge(wd / "migc.diag.txt")
    q0 = q_win[0]
    worst = max(abs(q_win[i] + sum(q for zd, q in drops if zd <= zr + 1e-9) - q0) / q0
                for i, zr in enumerate(z))
    tot_drop = sum(q for _, q in drops)
    c_ok = moved > 10000 and tot_drop > 0 and worst < 1e-10
    ok = ok and c_ok
    print(f"--- migration conservation: {moved} moves, {tot_drop:.3e} C dropped, "
          f"worst violation {worst:.2e} (tol 1e-10)  {'ok' if c_ok else 'FAIL'}")
    p_ok = bdev < 1e-10
    ok = ok and p_ok
    print(f"--- migration phase continuity: worst phasor deviation {bdev:.2e} "
          f"(tol 1e-10)  {'ok' if p_ok else 'FAIL'}")

    # Window residency: the routine's postcondition. The final dump is written right
    # after the last migration, so every surviving particle's theta must lie in its slice
    # window [0, 2*pi*sample). An unadjusted z on the move survives conservation (the
    # cascade drops are accounted) but not this.
    #
    # The window is the same [0, slen) for every slice, since theta is the phase inside a
    # slice and a move shifts z by exactly one spacing. The dump carries the lag rather
    # than the phase, and the writer folds the reference phase into it, so theta comes
    # back here as the absolute phase the tracker held.
    slen = 2 * np.pi * 3
    n_out = 0
    for sl in beamio.read_slices(wd / "migc-final.beam.h5", 1e-10, 3e-10):
        th = sl["theta"]
        n_out += int(np.sum((th < -1e-9) | (th >= slen + 1e-9)))
    w_ok = n_out == 0
    ok = ok and w_ok
    print(f"--- migration window residency: {n_out} particles outside their window "
          f"in the final dump  {'ok' if w_ok else 'FAIL'}")

    # 3: no-op bit identity on a frozen-phase beam.
    for mig, root in (("F", "mignf"), ("T", "mignt")):
        out = run(exe, wd, f"{root}.nml", BASE.format(lat="aramis_1seg.bmad", root=root,
                  sig_pz="8.804506566858e-08", emit="1e-13", npart=256, slen="1.2e-9", q="1.200830742713e-14", half="6e-10", mig=mig))
    moved = read_migration_file(wd, "mignt")[0]
    diag_eq = (wd / "mignf.diag.txt").read_bytes() == (wd / "mignt.diag.txt").read_bytes()
    d_eq = all(dumps_equal(wd / f"mignf-final.{s}.h5", wd / f"mignt-final.{s}.h5")
               for s in ("wf", "beam"))
    n_ok = moved == 0 and diag_eq and d_eq
    ok = ok and n_ok
    print(f"--- migration no-op: {moved} moves, diag byte-equal {diag_eq}, "
          f"dumps dataset-equal {d_eq}  {'ok' if n_ok else 'FAIL'}")

    print("migration checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
