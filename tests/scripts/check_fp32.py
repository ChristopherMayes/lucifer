#!/usr/bin/env python3
"""Checks of the FP32 lockstep instrument (performance.md, and validation.md at
val-fp32-lockstep). The divergence of the single-precision particle path from the FP64
path is a recorded number, the instrument leaves the FP64 physics alone, and the check
can fail.

1. Lockstep levels. A steady one-segment run with fp32_check = "lockstep" must land its
   worst divergence per quantity inside the recorded ceilings. This holds for the
   particle rows and for the field twin's source and post-solve rows. The residual
   guard must stay at or above 256 ulps, and the bucket renormalization round trip must
   stay within 1.5 ulps.
2. Read-only. The same input file with the instrument off must reproduce the instrumented
   run's diag file byte for byte. The twin observes the FP64 path and never steers it.
3. Falsifiable. fp32_mutate = T coarsens the FP32 residual by eight mantissa bits. The
   theta level must then grow by 10 times and the field twin's source row by 5 times,
   so a wrong twin on either side fails the check.
4. Compounding and the end-to-end number. Freerun carries the FP32 longitudinal state
   and the FP32 field across steps. The phasor must diverge past lockstep, and the
   end-to-end power divergence on this gain-regime segment must lie in (0, 1e-2]. That
   number is the figure a device port is judged against.
5. Time dependence. An 8-slice shot-noise window holds the same ceilings. The
   instrument's stream is byte-identical at 1 and 8 threads, since its work is per slice.
6. Refused. The twin mirrors the advance at the fundamental with no collective effects,
   so a wake run with the instrument on must be refused rather than measured wrongly.

Usage: check_fp32.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis_1seg.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

from checkrec import record
from lucifer.flat import to_toml as to_groups

FAILED = False

BASE = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "aramis_1seg.bmad"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 1024
  beam.bunch_charge = 1.000692285594e-15
  beam.sig_z = 0
  beam.sig_pz = 8.804506566858e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  ran_seed = 777
  write_diag = T
{extra}&end
"""

TD_EXTRA = """  beam.bunch_charge = 1.601107656951e-14
  beam.sig_z = 0
  beam.z_distribution = "GRID"
  beam.z_min = -1.200000e-09
  beam.z_max = 1.200000e-09
  seed_power = 0
  shot_noise = T
  beamlet_size = 8
  window_length = 2.4e-9
  n_wavelength = 3
"""

# The recorded ceilings are three times the levels of the cases below, which
# validation.md (val-fp32-lockstep) records. The production and debug builds give the
# same levels, since the twin's arithmetic is its own. A level past its ceiling means the
# error of the FP32 path changed.

CEIL = {"pz": 2e-6, "theta": 3e-4, "phasor": 3e-5, "source": 1e-4, "field": 1e-4}
GUARD_FLOOR = 256.0
RENORM_CEIL = 1.5


def ok(label, value, check, good, num=None, bound=None, kind="bool"):
    global FAILED
    tag = "ok" if good else "FAIL"
    print(f"--- {label}: {value} (check {check})  {tag}")

    # The printed value and check are strings. The record takes the number the caller
    # compared and the bound it compared against, and kind "bool" where there is none.

    record(label, num, bound, kind, good)
    if not good:
        FAILED = True


def run(exe, wd, name, text, threads="4", expect_fail=False):
    (wd / name).write_text(to_groups(text))
    r = subprocess.run([exe, name], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if expect_fail:
        return r
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:")
        print(r.stdout[-3000:])

        # A runtime error from the compiler goes to stderr and not to out_io. A run that
        # stops there prints nothing on stdout after its banner, so the tail of stderr is
        # printed too.

        if r.stderr:
            print(f"--- {name} stderr:")
            print(r.stderr[-3000:])
        sys.exit(1)
    return r


def summary(wd, root):
    vals = {}
    for line in (wd / f"{root}.fp32.txt").read_text().splitlines():
        if line.startswith("worst_"):
            k, v = line.split()
            vals[k.removeprefix("worst_")] = float(v)
        elif line.startswith("guard_ulp_min"):
            vals["guard"] = float(line.split()[1])
        elif line.startswith("endtoend_power_rel"):
            vals["e2e_power"] = float(line.split()[1])
        elif line.startswith("endtoend_bunching_abs"):
            vals["e2e_bunching"] = float(line.split()[1])
        elif line.startswith("# renorm_roundtrip_ulp"):
            vals["renorm"] = float(line.split()[2])
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)

    # 1. Steady lockstep: the recorded levels, particle and field rows alike.
    run(args.exe, wd, "fp32_ss.in", BASE.format(root="fp32ss",
        extra='  run.fp32_check = "lockstep"\n'))
    s = summary(wd, "fp32ss")
    for q in ("pz", "theta", "phasor", "source", "field"):
        ok(f"lockstep steady worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.0e}", s[q] <= CEIL[q],
           num=s[q], bound=CEIL[q], kind="max")
    ok("lockstep steady guard ulps", f"{s['guard']:.3e}", f">= {GUARD_FLOOR:.0f}",
       s["guard"] >= GUARD_FLOOR,
       num=s["guard"], bound=GUARD_FLOOR, kind="min")
    ok("bucket renormalization round trip [ulp]", f"{s['renorm']:.2f}", f"<= {RENORM_CEIL}",
       s["renorm"] <= RENORM_CEIL,
       num=s["renorm"], bound=RENORM_CEIL, kind="max")

    # 2. Read-only: the instrumented run's FP64 output is the uninstrumented run's.
    run(args.exe, wd, "fp32_off.in", BASE.format(root="fp32off", extra=""))
    same = (wd / "fp32ss.diag.txt").read_bytes() == (wd / "fp32off.diag.txt").read_bytes()
    ok("FP64 diag byte-identical, instrument on vs off", same, "True", same)

    # 3. Falsifiable: the mutation must move the theta level and the field twin's source row.
    run(args.exe, wd, "fp32_mut.in", BASE.format(root="fp32mut",
        extra='  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "fp32mut")
    ok("mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"],
       num=m["theta"], bound=10 * s["theta"], kind="min")
    ok("mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"],
       num=m["source"], bound=5 * s["source"], kind="min")

    # 4. Compounding: freerun must diverge past lockstep, and its end-to-end power
    # divergence on this gain-regime segment must lie in (0, 1e-2]. That number is the
    # figure a device port is judged against.
    run(args.exe, wd, "fp32_fr.in", BASE.format(root="fp32fr",
        extra='  run.fp32_check = "freerun"\n'))
    f = summary(wd, "fp32fr")
    ok("freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"],
       num=f["phasor"], bound=s["phasor"], kind="min")
    ok("freerun end-to-end power divergence recorded", f"{f['e2e_power']:.3e}",
       "<= 1e-2 on the gain segment", 0 < f["e2e_power"] <= 1e-2,
       num=f["e2e_power"], bound=[0, 1e-2], kind="range")

    # 5. Time dependence, and thread identity of the instrument's own stream.
    run(args.exe, wd, "fp32_td.in", BASE.format(root="fp32td",
        extra=TD_EXTRA + '  run.fp32_check = "lockstep"\n'), threads="8")
    t8 = (wd / "fp32td.fp32.txt").read_bytes()
    t = summary(wd, "fp32td")
    for q in ("pz", "theta", "phasor", "source", "field"):
        ok(f"lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.0e}", t[q] <= CEIL[q],
           num=t[q], bound=CEIL[q], kind="max")
    ok("lockstep TD guard ulps", f"{t['guard']:.3e}", f">= {GUARD_FLOOR:.0f}",
       t["guard"] >= GUARD_FLOOR,
       num=t["guard"], bound=GUARD_FLOOR, kind="min")
    run(args.exe, wd, "fp32_td.in", BASE.format(root="fp32td",
        extra=TD_EXTRA + '  run.fp32_check = "lockstep"\n'), threads="1")
    same = t8 == (wd / "fp32td.fp32.txt").read_bytes()
    ok("instrument stream byte-identical, 1 vs 8 threads", same, "True", same)

    # 6. A wake run lies outside the twin's coverage and must be refused.
    r = run(args.exe, wd, "fp32_ref.in", BASE.format(root="fp32ref",
        extra='  run.fp32_check = "lockstep"\n  wake_on = T\n  wake_radius = 2.5e-3\n'
              '  wake_conductivity = 5.813e7\n  wake_relaxation = 8.1e-6\n'),
        expect_fail=True)
    refused = r.returncode != 0 and "FP32_CHECK DOES NOT COVER WAKES" in r.stdout
    ok("wake with the instrument refused", refused, "True", refused)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
