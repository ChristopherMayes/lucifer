#!/usr/bin/env python3
"""
Load a Gaussian bunch under each shot-noise method and record its noise along the bunch
and across wavelength, for the figures and the band table of docs/shot-noise.md.

The bunch is N_SLICE slices of N_WAVELENGTH wavelengths in a window of twelve rms lengths,
PEAK electrons a beamlet in its peak slice, so its last slices hold 1e-5 of an electron.
It is loaded with run.load_only under beam.shot_noise_method = "weight", under
"first_order", with shot_noise = false, and with quiet_start = false and
shot_noise = false, at --seeds seeds each. From every particle
dump two things are taken, with S the sum of w exp(-i phase) over particles and w a
particle's electrons:

    slices    |S_1|^2 of each slice at the fundamental, the mean over the seeds, beside
              the slice's mean electrons
    spectrum  |S(k)|^2 over the bunch's electrons, from the particles' own places along
              the bunch, at the window's own frequencies from K_MIN to K_MAX resonant
              wavenumbers, the mean over the seeds

A real beam has |S_1|^2 equal to a slice's electrons and a spectrum of 1. The record
carries the mean spectrum in BANDS. A record whose weight load is outside BOUND of 1 in
any band is refused, and nothing is written.

The loads and the sums take a few minutes on a dozen cores, one process a seed, so the
script runs on demand and the
record is committed. shot_noise_figures.py draws from it and report_shot_noise.py writes
its band table, both in the keystone's regeneration.

Usage: shot_noise_spectrum.py --exe <lucifer> --workdir <dir> [--seeds 48] [--jobs N]
                              --out docs/generated/shot-noise/spectrum.json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import pathlib
import subprocess
import sys

import numpy as np

import lucifer.beam as beamio
from lucifer.flat import to_toml

E_CHARGE = 1.602176634e-19
LAMBDA0 = 1e-10
N_WAVELENGTH = 3
N_SLICE = 384
N_PARTICLE = 256
BEAMLET_SIZE = 8
PEAK = 50.0
K_MIN, K_MAX = 0.01, 8.6
BANDS = ((0.01, 0.3), (0.4, 0.6), (0.9, 1.1), (1.4, 1.6), (4.0, 5.0), (8.0, 8.6))
BOUND = 0.05

SPACING = N_WAVELENGTH * LAMBDA0
SIG_Z = N_SLICE * SPACING / 12
CHARGE = PEAK * (N_PARTICLE / BEAMLET_SIZE) * E_CHARGE * np.sqrt(2 * np.pi) * SIG_Z / SPACING

METHODS = {"weight": '  shot_noise = T\n  shot_noise_method = "weight"\n',
           "first_order": '  shot_noise = T\n  shot_noise_method = "first_order"\n',
           "quiet": "  shot_noise = F\n",
           "loud": "  quiet_start = F\n  shot_noise = F\n"}

TEXT = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = {lambda0}
  beam.n_particle = {n_particle}
  beam.bunch_charge = {charge:.12e}
  beam.sig_z = {sig_z:.12e}
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = 32
  grid_half_width = 2e-4
  beamlet_size = {beamlet_size}
  window_length = {window:.6e}
  n_wavelength = {n_wavelength}
{noise}  ran_seed = {seed}
  load_only = T
&end
"""


def one_seed(args):
    """One load and its sums: the spectrum, the slices' |S_1|^2 and the slices' electrons."""
    exe, workdir, lat, name, noise, seed, k = args
    root = f"spec_{name}_{seed}"
    (workdir / f"{root}.toml").write_text(to_toml(TEXT.format(
        lat=lat, root=root, lambda0=LAMBDA0, n_particle=N_PARTICLE, charge=CHARGE, sig_z=SIG_Z,
        beamlet_size=BEAMLET_SIZE, window=N_SLICE * SPACING, n_wavelength=N_WAVELENGTH, noise=noise,
        seed=100 + seed)))
    r = subprocess.run([exe, f"{root}.toml"], cwd=workdir, capture_output=True, text=True,
                       env={**os.environ, "OMP_NUM_THREADS": "1"})
    if r.returncode != 0:
        print(f"FAIL: the load {root} exited {r.returncode}:\n{r.stdout[-2000:]}", flush=True)
        return None
    spectrum = np.zeros(len(k))
    slices = np.zeros(N_SLICE)
    electrons = np.zeros(N_SLICE)
    phase, weight = [], []
    for i, sl in enumerate(beamio.read_slices(workdir / f"{root}-initial.beam.h5", LAMBDA0, SPACING)):
        n = sl["n"]
        if n == 0:
            continue
        w = sl["weight"][:n] / E_CHARGE
        theta = sl["theta"][:n]
        slices[i] = abs(np.sum(w * np.exp(-1j * theta))) ** 2
        electrons[i] = w.sum()
        phase.append(theta + 2 * np.pi * N_WAVELENGTH * i)
        weight.append(w)
    phase, weight = np.concatenate(phase), np.concatenate(weight)
    for j in range(0, len(k), 128):
        spectrum[j:j + 128] = np.abs(np.exp(-1j * np.outer(k[j:j + 128], phase)) @ weight) ** 2 / weight.sum()
    for f in workdir.glob(f"{root}[-.]*"):
        f.unlink()
    return spectrum, slices, electrons


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True, help="A directory that holds the lattice")
    p.add_argument("--lat", default="aramis_1seg.bmad")
    p.add_argument("--seeds", type=int, default=48)
    p.add_argument("--jobs", type=int, default=os.cpu_count(), help="Loads run at once")
    p.add_argument("--out", required=True, help="The record to write")
    a = p.parse_args()
    workdir = pathlib.Path(a.workdir).resolve()
    exe = str(pathlib.Path(a.exe).resolve())

    # The window's own frequencies, in resonant wavenumbers.
    n_window = N_SLICE * N_WAVELENGTH
    k = np.arange(int(np.ceil(K_MIN * n_window)), int(K_MAX * n_window) + 1) / n_window
    record = {"lambda0": LAMBDA0, "n_wavelength": N_WAVELENGTH, "n_slice": N_SLICE, "n_particle": N_PARTICLE,
              "beamlet_size": BEAMLET_SIZE, "peak_electrons_a_beamlet": PEAK, "sig_z": SIG_Z, "seeds": a.seeds,
              "bound": BOUND, "k": [float(f"{x:.6g}") for x in k], "methods": {}}
    for name, noise in METHODS.items():
        spectrum = np.zeros(len(k))
        slices = np.zeros(N_SLICE)
        electrons = np.zeros(N_SLICE)
        # One process a seed. The sums are added in the order of the seeds, so the record
        # does not depend on the number of processes.
        with concurrent.futures.ProcessPoolExecutor(a.jobs) as pool:
            for one in pool.map(one_seed, [(exe, workdir, a.lat, name, noise, seed, k)
                                           for seed in range(1, a.seeds + 1)]):
                if one is None:
                    return 1
                spectrum += one[0]
                slices += one[1]
                electrons += one[2]
        spectrum /= a.seeds
        bands = [float(spectrum[(k >= lo) & (k < hi)].mean()) for lo, hi in BANDS]
        record["methods"][name] = {"spectrum": [float(f"{x:.5g}") for x in spectrum],
                                   "slice_power": [float(f"{x:.5g}") for x in slices / a.seeds],
                                   "slice_electrons": [float(f"{x:.5g}") for x in electrons / a.seeds],
                                   "bands": [float(f"{x:.5g}") for x in bands]}
        print(f"{name}: bands " + " ".join(f"{lo:g}-{hi:g}: {b:.4g}" for (lo, hi), b in zip(BANDS, bands)), flush=True)
    record["bands"] = [list(b) for b in BANDS]

    worst = max(abs(b - 1) for b in record["methods"]["weight"]["bands"])
    if worst > BOUND:
        print(f"REFUSED: the weight load's spectrum is {worst:.3f} from 1 in a band, over the bound of {BOUND}. "
              f"No record is written.")
        return 1
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    print(f"wrote {a.out}: {len(k)} wavenumbers, {a.seeds} seeds, the weight load within {worst:.3f} of 1 in every band")
    return 0


if __name__ == "__main__":
    sys.exit(main())
