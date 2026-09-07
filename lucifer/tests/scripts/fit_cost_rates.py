#!/usr/bin/env python3
"""
Fit the cost estimate's rates to a set of measured runs (doc/performance.md).

code/fel_cost_mod.f90 estimates a run's walk from five constants. They are measured
rather than chosen, and the page they are measured from says that re-measuring it means
re-fitting them. This script is that fit, so the constants are never hand-derived twice.

The model is the module's: one rate for the particle push, linear in particle-steps; two
for the field solve, a deposit linear in particle-steps and a transform going as
ngrid^2 log2(ngrid) per slice per step; a serial fraction carrying the rates to any
thread count through Amdahl's law; and the share of the walk the push and the solve own.

Every run's configuration is read from its own header and its phase seconds from its own
Timing block, so the input is the logs the runs already wrote. Rates come out per thread.

Usage:

  fit_cost_rates.py <log> [<log> ...]

The set needs at least two particle counts at one slice count for the solve to separate,
and a thread sweep at one configuration for the serial fraction. Runs whose deck turned
the source filter off are reported and left out of the fit, since the estimate is for the
default.
"""

import argparse
import math
import re
import sys


def read_run(path):
    """One run's configuration and phase seconds, from the log it wrote."""
    txt = open(path, errors="ignore").read()
    r = {"log": path}

    m = re.search(r"^ Beam\s+(\d+) slices? x (\d+) particles", txt, re.M)
    if not m:
        m = re.search(r"^ Beam\s+(\d+) slices?, (\d+) particles", txt, re.M)
    if not m:
        return None
    r["nslice"], r["npart"] = int(m.group(1)), int(m.group(2))

    m = re.search(r"grid (\d+) points", txt)
    if not m:
        return None
    r["ngrid"] = int(m.group(1))

    m = re.search(r"^ Work\s+(\d+) FEL steps", txt, re.M)
    if not m:
        return None
    r["nstep"] = int(m.group(1))

    m = re.search(r"threads (\d+)", txt)
    r["threads"] = int(m.group(1)) if m else 1

    # A deck that turned the filter off is a different configuration: the solve then
    # runs one transform fewer per field per step. The run says so in its own warning.
    r["filtered"] = "global%source_filter = F" not in txt

    for key, pat in (("push", r"^\s+particle push\s+\d+\s+([\d.]+)"),
                     ("solve", r"^\s+field solve\s+\d+\s+([\d.]+)"),
                     ("walk", r"^\s+walk\s+\d+\s+([\d.]+)")):
        m = re.search(pat, txt, re.M)
        if not m:
            return None
        r[key] = float(m.group(1))

    r["pstep"] = r["nslice"] * r["npart"] * r["nstep"]
    r["gwork"] = r["nslice"] * r["ngrid"] ** 2 * math.log2(r["ngrid"]) * r["nstep"]
    return r


def speedup(f, n):
    return 1 / (f + (1 - f) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    a = ap.parse_args()

    runs = []
    for p in a.logs:
        r = read_run(p)
        if r is None:
            print(f"  skipped, no header or no Timing block: {p}")
            continue
        runs.append(r)
    if not runs:
        sys.exit("no runs read")

    off = [r for r in runs if not r["filtered"]]
    for r in off:
        print(f"  left out, the filter is off: {r['log']}")
    runs = [r for r in runs if r["filtered"]]
    if len(runs) < 2:
        sys.exit("need at least two filtered runs")

    # The serial fraction first, from whichever configuration was run at more than one
    # thread count. Least squares on the walk against Amdahl's law, over a grid: five
    # points and one parameter do not need anything cleverer.

    groups = {}
    for r in runs:
        groups.setdefault((r["nslice"], r["npart"], r["ngrid"], r["nstep"]), []).append(r)
    sweep = max(groups.values(), key=lambda g: len({r["threads"] for r in g}))
    if len({r["threads"] for r in sweep}) < 3:
        sys.exit("no thread sweep in the set: need one configuration at three thread counts")
    base = min(sweep, key=lambda r: r["threads"])
    best_f, best_err = None, None
    f = 0.0
    while f < 0.2:
        err = sum((math.log(r["walk"] / (base["walk"] * base["threads"] /
                                         speedup(f, r["threads"]) / base["threads"])))**2
                  for r in sweep if r["threads"] > 0)
        if best_err is None or err < best_err:
            best_f, best_err = f, err
        f += 0.0005
    serial = best_f

    # Everything else is fitted at one thread, so each run's phases are lifted by the
    # speedup its own thread count had.

    for r in runs:
        n_eff = max(1, min(r["threads"], r["nslice"]))
        r["lift"] = speedup(serial, n_eff)
        r["push1"] = r["push"] * r["lift"]
        r["solve1"] = r["solve"] * r["lift"]

    # The push, one rate over every run.
    push_rate = sum(r["push1"] for r in runs) / sum(r["pstep"] for r in runs)

    # The solve, two rates by least squares on solve1 = dep * pstep + fft * gwork.
    sxx = sum(r["pstep"] ** 2 for r in runs)
    syy = sum(r["gwork"] ** 2 for r in runs)
    sxy = sum(r["pstep"] * r["gwork"] for r in runs)
    sxb = sum(r["pstep"] * r["solve1"] for r in runs)
    syb = sum(r["gwork"] * r["solve1"] for r in runs)
    det = sxx * syy - sxy * sxy
    if det == 0:
        sys.exit("the solve cannot separate: every run has the same particles per slice")
    dep_rate = (sxb * syy - syb * sxy) / det
    fft_rate = (syb * sxx - sxb * sxy) / det

    # The walk share, the push and the solve over the walk.
    share = sum(r["push"] + r["solve"] for r in runs) / sum(r["walk"] for r in runs)

    print()
    print(f"{'run':34} {'threads':>7} {'walk':>8} {'model':>8} {'ratio':>6}")
    worst = 0.0
    for r in sorted(runs, key=lambda r: (r["nslice"], r["npart"], r["threads"])):
        t1 = push_rate * r["pstep"] + dep_rate * r["pstep"] + fft_rate * r["gwork"]
        est = t1 / r["lift"] / share
        ratio = est / r["walk"]
        worst = max(worst, abs(math.log(ratio)))
        tag = f"{r['nslice']}x{r['npart']}"
        print(f"{tag:34} {r['threads']:7d} {r['walk']:8.2f} {est:8.2f} {ratio:6.3f}")
    print(f"\nworst residual over the set: {100 * (math.exp(worst) - 1):.1f} percent")
    print()
    print("real(rp), parameter :: fel_cost_push$ = %.3g_rp" % push_rate)
    print("real(rp), parameter :: fel_cost_deposit$ = %.3g_rp" % dep_rate)
    print("real(rp), parameter :: fel_cost_solve$ = %.3g_rp" % fft_rate)
    print("real(rp), parameter :: fel_cost_serial$ = %.3g_rp" % serial)
    print("real(rp), parameter :: fel_cost_walk_share$ = %.3g_rp" % share)


if __name__ == "__main__":
    sys.exit(main())
