#!/usr/bin/env python3
"""
Statistical check for the weighted shot-noise loader (deliverable 6: a
feature Genesis cannot represent is tested against its own statistics, not a reference).

Physics under test: after quiet loading plus Fawley-style noise, each slice's bunching
satisfies <|b(h)|^2> = 1/N_lambda per imposed harmonic, with N_lambda the slice's real
electron count (charge/e) -- for uniform AND nonuniform per-particle weights.

Method: run lucifer with load_only = T over many seeds, in both weight modes
(uniform, and gen_test_weights = T which alternates beamlet weights 0.25x/1.75x at
constant charge). Read each .beam.h5, compute the charge-weighted |b(h)|^2 per slice for
the imposed harmonics, and test the scaled mean m = <|b(h)|^2 * N_lambda> against 1.
b is a sum of many independent beamlet contributions, so |b|^2*N_lambda is Exp(1) to
excellent approximation and the mean over n samples has sigma = 1/sqrt(n); the check is
|m - 1| < 5/sqrt(n) per weight mode, plus a looser per-harmonic check (n/3 samples).

Usage: check_shot_noise.py --exe <lucifer> --workdir <dir> [--seeds N] [--lat <bmad file>]
Exit 0 only if every check passes.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

import beamio
from nml import to_groups

E_CHARGE = 1.602176634e-19

# The window the deck below states. A dump carries the slice partition and not the
# radiation it was sliced on, so the reader is told.
LAMBDA0 = 1e-10
SPACING = 3 * LAMBDA0

NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 1024
  beam_init%bunch_charge = 4.803322970853e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -2.400000e-09
  beam_init%grid(3)%x_max = 2.400000e-09
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = 33
  grid_half_width = 2e-4
  beamlet_size = 8
  window_length = 4.8e-9
  window_sample = 3
  shot_noise = T
  gen_test_weights = {testw}
  ran_seed = {seed}
  load_only = T
&end
"""


def b2_samples(dump_file, harmonics):
    """Charge-weighted |b(h)|^2 * N_lambda for every slice and harmonic.

    The dump is openPMD, so the weights are the ones the loader wrote. An earlier
    version read Genesis format, which carries one current per slice, and had to
    reconstruct the alternating test pattern to weight the sum. Reading the weights
    tests the loader instead of assuming it.

    |b(h)| is invariant under a constant shift of theta, which is what a dump's missing
    reference phase amounts to, so this measurement is unaffected by it."""
    out = []
    for sl in beamio.read_slices(dump_file, LAMBDA0, SPACING):
        if sl["n"] == 0:
            continue
        n_lambda = sl["weight"].sum() / E_CHARGE
        w, theta = sl["weight"], sl["theta"]
        for h in harmonics:
            b = np.sum(w * np.exp(-1j * h * theta)) / np.sum(w)
            out.append(abs(b) ** 2 * n_lambda)
    return out


def run_mode(exe, lat, workdir, seeds, test_weights):
    samples = {h: [] for h in (1, 2, 3)}
    for seed in range(1, seeds + 1):
        root = f"sn_{'w' if test_weights else 'u'}_{seed}"
        nml = workdir / f"{root}.nml"
        nml.write_text(to_groups(NML.format(lat=lat, root=root, seed=1000 + 7 * seed,
                                             testw="T" if test_weights else "F")))
        r = subprocess.run([str(exe), nml.name], cwd=workdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        vals = b2_samples(workdir / f"{root}-initial.beam.h5", (1, 2, 3))
        for k, h in enumerate((1, 2, 3)):
            samples[h].extend(vals[k::3])
    return samples


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--lat", default="aramis_1seg.bmad")
    p.add_argument("--seeds", type=int, default=25)
    global PARSED
    PARSED = p.parse_args()

    workdir = pathlib.Path(PARSED.workdir)
    exe = pathlib.Path(PARSED.exe).resolve()

    ok = True
    for test_weights, label in ((False, "uniform weights"), (True, "nonuniform weights (0.25x/1.75x)")):
        samples = run_mode(exe, PARSED.lat, workdir, PARSED.seeds, test_weights)
        allv = np.concatenate([samples[h] for h in (1, 2, 3)])
        n = len(allv)
        m = allv.mean()
        bound = 5 / np.sqrt(n)
        good = abs(m - 1) < bound
        ok = ok and good
        print(f"--- shot noise, {label}: <|b(h)|^2 * N_lambda> = {m:.4f} "
              f"(target 1 +- {bound:.3f}, {n} samples)  {'ok' if good else 'FAIL'}")
        for h in (1, 2, 3):
            v = np.asarray(samples[h])
            mh = v.mean()
            bh = 5 / np.sqrt(len(v))
            goodh = abs(mh - 1) < bh
            ok = ok and goodh
            print(f"      harmonic {h}: {mh:.4f} (+- {bh:.3f})  {'ok' if goodh else 'FAIL'}")

    print("shot-noise statistical check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
