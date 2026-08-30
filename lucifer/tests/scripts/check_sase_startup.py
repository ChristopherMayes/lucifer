#!/usr/bin/env python3
"""
Cross-code validation of the shot-noise level (deliverable 6): a dark start (no seed
field), one undulator segment, each code generating its own noisy beam -- Genesis with
its loader and RNG, lucifer with its weighted Fawley loader -- so the two noise
implementations are fully independent. SASE startup power is proportional to the
imposed <|b|^2>, and the tracking itself is validated elsewhere at 1e-6, so agreement of
the mean startup power validates the noise level and nothing else.

Statistics: per side, mean over slices and seeds of the final-record slice power. Slice
powers are near-exponentially distributed, so with n samples the mean carries a relative
sigma of about 1/sqrt(n); the check is |ln(P_bmad/P_genesis)| < 0.30, roughly 3 sigma at
the default 6 seeds x 32 slices per side and far below the factor-type errors a wrong
noise normalization produces.

Usage: check_sase_startup.py --exe <lucifer> --genesis <genesis4> --workdir <dir> [--seeds N]
The workdir must hold Aramis-1seg.lat and aramis_1seg.bmad. Exit 0 on pass.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from nml import to_groups

import pool

GENESIS_DECK = """&setup
rootname=SASE{seed}
lattice=Aramis-1seg.lat
beamline=SEG1
lambda0=1e-10
gamma0=11357.82
delz=0.045000
shotnoise=1
npart = 2048
nbins = 8
seed = {seed}
beam_global_stat = true
field_global_stat = true
&end

&time
slen = 9.6e-9
sample = 3
&end

&field
power=0
dgrid=2.000000e-04
ngrid=255
waist_size=30e-6
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
&end

&track
fft_fieldsolver = true
&end
"""

BMAD_NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis_1seg.bmad"
  out_root = "bsase{seed}"
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 9.606645941707e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4.800000e-09
  beam_init%grid(3)%x_max = 4.800000e-09
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = 255
  grid_half_width = 2e-4
  nbins = 8
  window_length = 9.6e-9
  window_sample = 3
  shotnoise = T
  ran_seed = {seed}
  write_diag = T
&end
"""


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--genesis", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--seeds", type=int, default=6)
    a = p.parse_args()

    wd = pathlib.Path(a.workdir)
    exe = pathlib.Path(a.exe).resolve()
    env = dict(os.environ, FI_PROVIDER="tcp", OMP_NUM_THREADS="8")

    # Every seed's two runs are independent processes. The pool runs them all.
    def gen_one(seed):
        (wd / f"gsase{seed}.in").write_text(GENESIS_DECK.format(seed=seed))
        r = subprocess.run([a.genesis, f"gsase{seed}.in"], cwd=wd, env=env,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: Genesis seed {seed} exited {r.returncode}");  sys.exit(1)

    def bmad_one(seed):
        (wd / f"bsase{seed}.nml").write_text(to_groups(BMAD_NML.format(seed=seed)))
        r = subprocess.run([str(exe), f"bsase{seed}.nml"], cwd=wd, env=env,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: lucifer seed {seed} exited {r.returncode}:\n{r.stdout[-1500:]}")
            sys.exit(1)

    seeds = [100000 + 991 * s for s in range(1, a.seeds + 1)]
    pool.run_all([lambda sd=sd: gen_one(sd) for sd in seeds]
                 + [lambda sd=sd: bmad_one(sd) for sd in seeds], threads_per_job=8)

    gp, bp = [], []
    for seed in seeds:
        with h5py.File(wd / f"SASE{seed}.out.h5") as h5:
            gp.extend(h5["Field/power"][-1, :])
        d = np.loadtxt(wd / f"bsase{seed}.diag.txt")
        nslice = int(d[:, 1].max())
        bp.extend(d.reshape(-1, nslice, d.shape[1])[-1, :, 2])

    gp, bp = np.asarray(gp), np.asarray(bp)
    ratio = bp.mean() / gp.mean()
    n = min(len(gp), len(bp))
    bound = 0.30
    ok = abs(np.log(ratio)) < bound
    print(f"--- SASE startup, dark start, one segment, {a.seeds} seeds x {n // a.seeds} slices/side")
    print(f"    Genesis  <P>: {gp.mean():.4e} W   (slice spread {gp.std()/gp.mean():.2f})")
    print(f"    Bmad     <P>: {bp.mean():.4e} W   (slice spread {bp.std()/bp.mean():.2f})")
    print(f"    ln ratio: {np.log(ratio):+.3f}  (bound +-{bound})  {'ok' if ok else 'FAIL'}")
    print("sase startup cross-check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
