#!/usr/bin/env python3
"""
Self-referenced gates for the collective terms (deliverable 8; the Genesis-comparison
tiers live in compare_fel.py).

1. Energy bookkeeping. A dark, quiet, cold beam (no seed, no shot noise, negligible
   energy spread) tracked with all wake kernels on has exactly one energy channel: the
   applied per-slice eloss. At every consecutive record pair, the measured
   d<gamma> must equal eloss(slice)*dz/m_electron to the diagnostic's print resolution,
   and sigma_gamma must stay constant (the kick is uniform within a slice).

2. Stale-wake structure. With migration on and a large energy spread, the current
   profile changes, and the wake convolution must follow: <out_root>.wake.txt must hold
   MORE THAN ONE z-stamped eloss block, and the blocks must differ. Removing the
   migration-stride recompute (the stale-wake mutation) leaves one block and fails
   loudly. The checker parses the driver's record of recomputes rather than
   reimplementing the convolution (FINDINGS 4.3).

Usage: check_collective.py --exe <fel_track_test> --workdir <dir>
The workdir must hold aramis_1seg.bmad and aramis.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

import numpy as np

M_ELECTRON = 0.51099895069e6

BASE = """&fel_track_params
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  gen_current = 3000
  gen_delgam = {delgam}
  gen_ex = 4e-7, gen_ey = 4e-7
  gen_beta_x = 8.53711,  gen_alpha_x = -0.703306
  gen_beta_y = 17.3899,  gen_alpha_y = 1.40348
  gen_power = 0
  gen_ngrid = 65
  gen_dgrid = 2e-4
  gen_npart = 512
  gen_nbins = 8
  gen_slen = 2.4e-9
  gen_sample = 3
  gen_seed = 555
  wake_on = T
  wake_radius = 2.5e-3
  wake_conductivity = 5.813e7
  wake_relaxation = 8.1e-6
  wake_gap = 0.5e-3
  wake_lgap = 0.015
  wake_hrough = 100e-9
  wake_lrough = 100e-6
{extra}&end
"""


def run(exe, wd, name, text):
    (wd / name).write_text(text)
    r = subprocess.run([str(exe), name], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    return r.stdout


def wake_blocks(fn):
    """[(z, eloss array)] per z-stamped block."""
    blocks, z, vals = [], None, []
    for line in pathlib.Path(fn).read_text().splitlines():
        m = re.match(r"# z =\s+(\S+)", line)
        if m:
            if z is not None:
                blocks.append((z, np.array(vals)))
            z, vals = float(m.group(1)), []
        elif line.strip() and not line.startswith("#"):
            vals.append(float(line.split()[1]))
    if z is not None:
        blocks.append((z, np.array(vals)))
    return blocks


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    exe = pathlib.Path(a.exe).resolve()
    ok = True

    # 1. Energy bookkeeping on a cold dark beam, one segment.
    run(exe, wd, "colle.nml", BASE.format(lat="aramis_1seg.bmad", root="colle",
                                          delgam="0.001", extra=""))
    eloss = wake_blocks(wd / "colle.wake.txt")[0][1]
    d = np.loadtxt(wd / "colle.diag.txt")
    ns = int(d[:, 1].max())
    d = d.reshape(-1, ns, d.shape[1])
    z, mg, sg = d[:, 0, 0], d[:, :, 6], d[:, :, 7]
    dz = np.diff(z)
    dg_meas = np.diff(mg, axis=0)
    dg_exp = eloss[None, :] * dz[:, None] / M_ELECTRON
    err = np.abs(dg_meas - dg_exp).max()
    scale = np.abs(dg_exp).max()
    e_ok = err < 1e-6 * scale + 1e-7
    ok = ok and e_ok
    print(f"--- collective energy bookkeeping: max |d<gamma> - eloss*dz/mc^2| = {err:.2e}"
          f" against kicks of {scale:.2e}  {'ok' if e_ok else 'FAIL'}")
    s_dev = np.abs(np.diff(sg, axis=0)).max()
    s_ok = s_dev < 1e-8
    ok = ok and s_ok
    print(f"--- collective sigma_gamma invariance under uniform kicks: {s_dev:.2e}  "
          f"{'ok' if s_ok else 'FAIL'}")

    # 2. Stale-wake structure: heavy migration must force recomputes that change eloss.
    run(exe, wd, "collm.nml", BASE.format(lat="aramis.bmad", root="collm",
        delgam="60.0", extra="  migrate = T\n"))
    blocks = wake_blocks(wd / "collm.wake.txt")
    changed = len(blocks) > 1 and any(
        not np.array_equal(blocks[0][1], b[1]) for _, b in [(None, blk) for _, blk in blocks[1:]])
    m_ok = len(blocks) > 1 and changed
    ok = ok and m_ok
    print(f"--- stale-wake structure: {len(blocks)} eloss blocks under heavy migration, "
          f"changed = {changed}  {'ok' if m_ok else 'FAIL'}")

    print("collective gates:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
