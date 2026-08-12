#!/usr/bin/env python3
"""
Statistical gate for the weighted shot-noise loader (deliverable 6, FINDINGS 6.9: a
feature Genesis cannot represent is tested against its own statistics, not a reference).

Physics under test: after quiet loading plus Fawley-style noise, each slice's bunching
satisfies <|b(h)|^2> = 1/N_lambda per imposed harmonic, with N_lambda the slice's real
electron count (charge/e) -- for uniform AND nonuniform per-particle weights.

Method: run fel_track_test with load_only = T over many seeds, in both weight modes
(uniform, and gen_test_weights = T which alternates beamlet weights 0.25x/1.75x at
constant charge). Read each .par.h5, compute the charge-weighted |b(h)|^2 per slice for
the imposed harmonics, and test the scaled mean m = <|b(h)|^2 * N_lambda> against 1.
b is a sum of many independent beamlet contributions, so |b|^2*N_lambda is Exp(1) to
excellent approximation and the mean over n samples has sigma = 1/sqrt(n); the gate is
|m - 1| < 5/sqrt(n) per weight mode, plus a looser per-harmonic check (n/3 samples).

Usage: check_shot_noise.py --exe <fel_track_test> --workdir <dir> [--seeds N] [--lat <bmad file>]
Exit 0 only if every check passes.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

E_CHARGE = 1.602176634e-19

NML = """&fel_track_params
  lat_file = "{lat}"
  out_root = "{root}"
  gamma0 = 11357.82
  lambda0 = 1e-10
  delz = 0.045
  gen_current = 3000
  gen_delgam = 1.0
  gen_ex = 4e-7, gen_ey = 4e-7
  gen_beta_x = 8.53711,  gen_alpha_x = -0.703306
  gen_beta_y = 17.3899,  gen_alpha_y = 1.40348
  gen_power = 0
  gen_ngrid = 33
  gen_dgrid = 2e-4
  gen_npart = 1024
  gen_nbins = 8
  gen_slen = 4.8e-9
  gen_sample = 3
  gen_shotnoise = T
  gen_test_weights = {testw}
  gen_seed = {seed}
  load_only = T
&end
"""


def b2_samples(par_file, harmonics):
    """Charge-weighted |b(h)|^2 * N_lambda for every slice and harmonic."""
    out = []
    with h5py.File(par_file) as h5:
        nslice = int(h5["slicecount"][0])
        spacing = float(h5["slicespacing"][0])
        for i in range(nslice):
            g = h5[f"slice{i+1:06d}"]
            theta = g["theta"][:]
            current = float(g["current"][0])
            n_lambda = current * spacing / (E_CHARGE * 2.99792458e8)
            # The dump format carries no weights; reconstruct the loader's test pattern:
            # uniform within beamlets of nbins=8, alternating 0.25/1.75 across beamlets
            # when the test mode is on. The generator wrote current = c*sum(w)/spacing,
            # so relative weights are all the statistic needs.
            w = np.ones(len(theta))
            if PARSED.reconstruct_weights:
                nb = 8
                scale = np.where(np.arange(len(theta)) // nb % 2 == 0, 0.25, 1.75)
                w = scale
            for h in harmonics:
                b = np.sum(w * np.exp(-1j * h * theta)) / np.sum(w)
                out.append(abs(b) ** 2 * n_lambda)
    return out


def run_mode(exe, lat, workdir, seeds, test_weights):
    samples = {h: [] for h in (1, 2, 3)}
    for seed in range(1, seeds + 1):
        root = f"sn_{'w' if test_weights else 'u'}_{seed}"
        nml = workdir / f"{root}.nml"
        nml.write_text(NML.format(lat=lat, root=root, seed=1000 + 7 * seed,
                                  testw="T" if test_weights else "F"))
        r = subprocess.run([str(exe), nml.name], cwd=workdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        PARSED.reconstruct_weights = test_weights
        vals = b2_samples(workdir / f"{root}-initial.par.h5", (1, 2, 3))
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
    PARSED.reconstruct_weights = False

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

    print("shot-noise statistical gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
