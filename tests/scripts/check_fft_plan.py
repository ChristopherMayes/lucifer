#!/usr/bin/env python3
"""
Checks of the transverse transform's plan (run.fft_plan, docs/input-reference.md).

The transform of one plane is FFTW's, and FFTW can make more than one plan for it. Two
plans give outputs that differ in the last bit, so a run's digits belong to its plan.
`wisdom`, the default, makes the measured plan the source records
(src/wavefront_wisdom_mod.f90) and makes it alike in every run. `estimate` makes FFTW's
estimated plan, also alike in every run. `measure` measures a plan in the run itself,
which can differ from one run to the next. Where the recorded text holds no plan for a
transform, or the FFTW of the build does not accept the text, `wisdom` makes the
estimated plan, the run says so on a warning line, and the stats file and the resolved
input file name `estimate` as the plan in use.

One time-dependent probe, eight slices on a 128-point grid, is run:

  1. Under `wisdom` twice. Every array of the two stats files is the same.
  2. Under `estimate`. Where the plan in use under `wisdom` is the recorded one, the two
     runs differ in at least one array, which shows the recorded plan is the one that
     ran, and the window energy of the two agrees to RTOL. Where it is the
     estimated one, as on a platform the text was not measured on, every array of the
     `wisdom` run is the `estimate` run's.
  3. Under `measure`, held to `estimate` in the window energy at RTOL. The difference
     is printed and not recorded, since the plan of this run can differ from the next.
  4. On a 127-point grid under `wisdom`, a size the recorded text does not hold. The
     plan in use is `estimate`, the warning line is printed, the resolved input file
     names `estimate`, and every array is that of the same run under `estimate`.
  5. At one thread and at eight under `wisdom`. Every array is the same.

The script also prints the field solve's seconds on a 256-point grid under `estimate`
and under `measure`. They are printed for a reader of the log and are not judged.

Run by the validation suite. Exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from checkrec import record

from lucifer.flat import to_toml

FAILED = False
# Two plans of one transform differ in the last bits. The window energy of the probe
# under the recorded plan and under the estimated one differs by 9.6e-16 of itself.
RTOL = 1e-10
SKIP = ("params", "meta", "timing", "setup/fft_plan")

FLAT = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "crossed_probe.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam.n_particle = 2048
  beam.bunch_charge = 8.0e-15
  beam.z_distribution = "GRID"
  beam.z_min = -4e-10
  beam.z_max = 4e-10
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = {grid}
  grid_half_width = 2e-4
  window_length = {window}
  n_wavelength = 1
  shot_noise = T
  ran_seed = 777
  keep_escaped_field = F
  fft_plan = "{plan}"
&end
"""


def check(name, value, bound, kind="max"):
    """Print and record one check of a number against its bound."""
    global FAILED
    ok = {"max": value <= bound, "min": value >= bound, "equal": value == bound}[kind]
    print(f"--- {name}: {value:.3e} ({kind} {bound:.0e}) {'ok' if ok else '** FAIL **'}")
    record(name, value, bound, kind, ok)
    FAILED = FAILED or not ok


def check_true(name, ok, note=""):
    """Print and record one check that has no number."""
    global FAILED
    print(f"--- {name}: {note} {'ok' if ok else '** FAIL **'}")
    record(name, None, None, "bool", ok)
    FAILED = FAILED or not ok


def run(exe, wd, root, plan, grid=128, window="8e-10", threads="8"):
    """Run the probe under one plan and return the program's printout."""
    (wd / f"{root}.toml").write_text(to_toml(FLAT.format(root=root, plan=plan, grid=grid, window=window)))
    r = subprocess.run([str(exe), f"{root}.toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
        sys.exit(1)
    return r.stdout


def arrays(wd, root):
    """Every dataset of a run's stats file, by path, as raw bytes."""
    out = {}
    with h5py.File(wd / f"{root}.stats.h5") as h:
        def walk(name, obj):
            if isinstance(obj, h5py.Dataset) and not name.startswith(SKIP):
                out[name] = np.asarray(obj[()]).tobytes()
        h.visititems(walk)
    return out


def differing(wd, a, b):
    """The number of stats arrays that differ between two runs."""
    x, y = arrays(wd, a), arrays(wd, b)
    return len(set(x) ^ set(y)) + sum(1 for k in x if k in y and x[k] != y[k])


def plan_in_use(wd, root):
    """The plan the run's stats file names."""
    with h5py.File(wd / f"{root}.stats.h5") as h:
        v = h["setup/fft_plan"][()]
    return v.decode() if isinstance(v, bytes) else str(v)


def energy_difference(wd, a, b):
    """The largest relative difference of the window energy along the line."""
    with h5py.File(wd / f"{a}.stats.h5") as h:
        x = np.asarray(h["field/total/energy"])
    with h5py.File(wd / f"{b}.stats.h5") as h:
        y = np.asarray(h["field/total/energy"])
    return float(np.max(np.abs(x - y)) / np.max(np.abs(y)))


def solve_seconds(wd, root):
    """The seconds the run's field solve took."""
    with h5py.File(wd / f"{root}.stats.h5") as h:
        return float(np.asarray(h["timing/field_solve"]).ravel()[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    exe = pathlib.Path(args.exe).resolve()
    wd = pathlib.Path(args.workdir).resolve()

    # 1. The default plan, twice.
    run(exe, wd, "fp_w1", "wisdom")
    run(exe, wd, "fp_w2", "wisdom")
    check("wisdom twice: stats arrays that differ", differing(wd, "fp_w1", "fp_w2"), 0, "equal")
    used = plan_in_use(wd, "fp_w1")
    print(f"--- the plan in use under wisdom on a 128-point grid: {used}")
    check_true("the plan in use is wisdom or estimate", used in ("wisdom", "estimate"), used)

    # 2. The estimated plan.
    run(exe, wd, "fp_e", "estimate")
    check_true("estimate: the stats file names estimate", plan_in_use(wd, "fp_e") == "estimate")
    if used == "wisdom":
        check("wisdom against estimate: stats arrays that differ", differing(wd, "fp_w1", "fp_e"), 1, "min")
        check("wisdom against estimate: window energy", energy_difference(wd, "fp_w1", "fp_e"), RTOL)
    else:
        check("wisdom without a plan: stats arrays that differ from estimate",
              differing(wd, "fp_w1", "fp_e"), 0, "equal")

    # 3. A plan measured in the run.
    run(exe, wd, "fp_m", "measure")
    check_true("measure: the stats file names measure", plan_in_use(wd, "fp_m") == "measure")
    d_m = energy_difference(wd, "fp_m", "fp_e")
    check_true(f"measure against estimate: window energy within {RTOL:.0e}", d_m <= RTOL, f"{d_m:.3e}")

    # 4. A size the recorded text does not hold.
    out = run(exe, wd, "fp_w127", "wisdom", grid=127)
    run(exe, wd, "fp_e127", "estimate", grid=127)
    check_true("127-point grid under wisdom: the plan in use is estimate",
               plan_in_use(wd, "fp_w127") == "estimate")
    check_true("127-point grid under wisdom: the warning line",
               "! fft_plan = wisdom has no plan for a 127 by 127 transform" in out)
    check_true("127-point grid under wisdom: the resolved input file names estimate",
               'fft_plan = "estimate"   # derived' in (wd / "fp_w127.resolved.toml").read_text())
    check("127-point grid: stats arrays that differ from estimate",
          differing(wd, "fp_w127", "fp_e127"), 0, "equal")

    # 5. One thread and eight.
    run(exe, wd, "fp_w1t", "wisdom", threads="1")
    check("wisdom at 1 and 8 threads: stats arrays that differ", differing(wd, "fp_w1", "fp_w1t"), 0, "equal")

    # The field solve's seconds on a 256-point grid, for a reader of the log.
    for plan in ("estimate", "measure", "wisdom"):
        run(exe, wd, f"fp_t_{plan}", plan, grid=256, window="16e-10")
        print(f"    field solve on a 256-point grid, {plan}: {solve_seconds(wd, f'fp_t_{plan}'):.3f} s "
              f"(plan in use {plan_in_use(wd, f'fp_t_{plan}')}, not judged)")

    if FAILED:
        print("FAIL: fft-plan checks")
        sys.exit(1)
    print("PASS: fft-plan checks")


if __name__ == "__main__":
    main()
