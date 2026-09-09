#!/usr/bin/env python3
"""
Statistical check for the weighted shot-noise loader (a
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

Fawley's kick is linear in the noise amplitude only while a group holds many electrons,
so every deck here keeps about 150 electrons per group, as the generator's has. At 18 the
scaled mean sits at 0.88 and the higher harmonics far lower, in every mode alike, which is
the algorithm's own limit (Genesis4 has the same one) and not a load defect.

The same statistic then runs on the loads that do not draw their own beamlets
(fel-physics.md sec-noise): load_mode = "keep" from a generated bunch, where the loader
makes the beamlets; a bunch of copies made outside the loader and loaded without the quiet
start, where the noise goes on the particles sharing their transverse coordinates; and the
same copies with their coordinates jittered so none coincide, where it goes on the
occupants of each deposit cell. The loader's message must name the grouping it used.

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

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fieldio  # noqa: E402

import beamio
import bunchfile
from nml import to_groups

E_CHARGE = 1.602176634e-19

# The window the deck below states. A dump carries the slice partition and not the
# radiation it was sliced on, so the reader is told.
LAMBDA0 = 1e-10
SPACING = 3 * LAMBDA0

NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
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
  grid_n_pts = 32
  grid_half_width = 2e-4
  beamlet_size = 8
  window_length = 4.8e-9
  n_wavelength = 3
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


KEEP_NML = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%source_filter = F
  global%ran_seed = {seed}
  global%load_only = T
  slicing%n_wavelength = 3
/
&fel_beam_init
  load_mode = "keep"
  quiet_start = {quiet}
  shot_noise = T
  beamlet_size = 8
{source}/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 0
  wavefront_init%grid_n_pts = 32
  wavefront_init%grid_half_width = 2e-4
/
"""

KEEP_BUNCH = """  beam_init%n_particle = 2048
  beam_init%bunch_charge = 4.803322970853e-14
  beam_init%sig_z = 1.2e-9
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
"""


def run_keep(exe, lat, workdir, seeds, label, source, quiet, want):
    """The statistic on a load the loader did not make the beamlets of, or made in keep mode."""
    samples = {h: [] for h in (1, 2, 3)}
    for seed in range(1, seeds + 1):
        root = f"sn_{label}_{seed}"
        (workdir / f"{root}.nml").write_text(KEEP_NML.format(lat=lat, root=root, seed=2000 + 7 * seed,
                                                              quiet=quiet, source=source))
        r = subprocess.run([str(exe), f"{root}.nml"], cwd=workdir, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        if want not in r.stdout:
            print(f"FAIL: {root} did not report the noise on {want!r}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        vals = b2_samples(workdir / f"{root}-initial.beam.h5", (1, 2, 3))
        for k, h in enumerate((1, 2, 3)):
            samples[h].extend(vals[k::3])
    return samples


def judge(samples, label):
    """The scaled mean against 1, over all harmonics and per harmonic."""
    ok = True
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
    return ok


def other_loads(exe, lat, workdir, seeds):
    """Keep mode, then copies and cells made outside the loader from one written bunch."""
    ok = True
    src = KEEP_BUNCH + '  write_openpmd_file = "sn_bunch.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, seeds, "keep", src, "T", "the beamlets the loader made"),
               "keep mode, beamlets of a generated bunch") and ok

    b = bunchfile.read_bunch(workdir / "sn_bunch.h5")
    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_copies.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING))
    src = '  beam_init%position_file = "sn_copies.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, max(3, seeds // 3), "copies", src, "F",
                        "particles sharing their transverse coordinates"),
               "copies made outside the loader") and ok

    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_cells.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING, jitter=1e-12))
    src = '  beam_init%position_file = "sn_cells.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, max(3, seeds // 3), "cells", src, "F",
                        "the occupants of each deposit cell"),
               "the occupants of each deposit cell") and ok
    return ok


SPLIT_NML = """  lat_file = "{lat}"
  out_root = "sn_split"
  source_filter = F
  dump_field_at = "END"
  ran_seed = 4321
  lambda0 = 1e-10
  beam_init%n_particle = 1024
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
&end
"""


def split_check(exe, lat, workdir):
    """
    The stats file's mode split against the same quantity taken from a field dump.

    Every record that takes the field angle moments also reports the power within
    split_angle of the axis, from the transform those moments already pay for
    (fel-physics.md sec-source-filter). A dump's own far field is the independent route to
    the same number, and Parseval makes them equal, so a difference is a bug in the
    masking or the normalization rather than a tolerance. One seeded steady-state segment
    is enough: the split does not care what the field is, only how it is resolved.
    """
    root = "sn_split"
    (workdir / f"{root}.nml").write_text(to_groups(SPLIT_NML.format(lat=lat)))
    r = subprocess.run([str(exe), f"{root}.nml"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: split run exited {r.returncode}:\n{r.stdout[-2000:]}")
        return False

    dumps = sorted(workdir.glob(f"{root}-at*.wf.h5")) + sorted(workdir.glob(f"{root}-final.wf.h5"))
    if not dumps:
        print("FAIL: the split run wrote no field dump")
        return False

    with h5py.File(workdir / f"{root}.stats.h5") as h:
        g = h["field/total"]
        angle = float(np.ravel(g["split_angle"][()])[0])
        p_in = np.atleast_1d(g["power_inside_angle"][-1])

    f = fieldio.read_field(str(dumps[-1]))
    u, dx = f["u"], f["dx"]
    n = u.shape[-1]
    kx = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    theta = np.sqrt(kx[None, :] ** 2 + kx[:, None] ** 2) / (2 * np.pi / f["wavelength"])
    spec = np.abs(np.fft.fft2(u, axes=(1, 2))) ** 2 / (n * n)
    inside = (spec * (theta <= angle)).sum(axis=(1, 2)) * dx * dx / (2 * fieldio.MU0_C)

    live = np.isfinite(p_in) & (inside > 0)
    if not live.any():
        print("FAIL: no slice carried both a split and a dump")
        return False
    worst = float(np.max(np.abs(p_in[live] - inside[live]) / inside[live]))
    good = worst <= 1e-6
    print(f"--- mode split against the dump's far field: worst {worst:.3e} over "
          f"{int(live.sum())} slices, angle {angle:.3e} rad (check <= 1e-6)  "
          f"{'ok' if good else 'FAIL'}")
    return good


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
        ok = judge(run_mode(exe, PARSED.lat, workdir, PARSED.seeds, test_weights), label) and ok

    ok = other_loads(exe, PARSED.lat, workdir, PARSED.seeds) and ok
    ok = split_check(exe, PARSED.lat, workdir) and ok

    print("shot-noise statistical check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
