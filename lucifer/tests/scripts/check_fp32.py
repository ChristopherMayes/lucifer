#!/usr/bin/env python3
"""FP32 lockstep checks (doc/performance.md): the single-precision particle path's
divergence from the FP64 reference is a recorded number, the instrument cannot touch
the FP64 physics, and the check can fail.

1. Lockstep levels. A steady one-segment run with fp32_check = "lockstep" must land
   its worst per-quantity divergence inside the recorded ceilings, particle rows and
   the field twin's source and post-solve rows alike, its residual guard must stay in
   healthy ulps, and the bucket-renormalization round trip must be exact to one ulp.
2. Read-only. The same deck with the instrument off must reproduce the instrumented
   run's diag file byte for byte: the twin observes the FP64 path and never steers it.
3. Falsifiable. fp32_mutate = T coarsens the FP32 residual by eight mantissa bits, and
   the theta level and the field twin's source row must both move, so a wrong twin on
   either side cannot hide behind a check that always passes.
4. Compounding and the end-to-end number. freerun carries the FP32 longitudinal state
   and the FP32 field across steps: the phasor must diverge past lockstep, and the
   recorded end-to-end power divergence on this gain-regime segment is the yardstick a
   device port will be judged against.
5. Time dependence. An 8-slice shot-noise window holds the same ceilings, and the
   instrument's stream is byte-identical at 1 and 8 threads (per-slice work only).
6. Refused by name. The twin mirrors the fundamental-only collective-free advance, so
   a wake run with the instrument on must refuse, not measure wrongly in silence.

Usage: check_fp32.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis_1seg.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

from nml import to_groups

FAILED = False

BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 1024
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 65
  grid_half_width = 2e-4
  ran_seed = 777
  write_diag = T
{extra}&end
"""

TD_EXTRA = """  beam_init%bunch_charge = 1.601107656951e-14
  beam_init%sig_z = 0
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.200000e-09
  beam_init%grid(3)%x_max = 1.200000e-09
  seed_power = 0
  shot_noise = T
  beamlet_size = 8
  window_length = 2.4e-9
  window_sample = 3
"""

# The recorded ceilings, set at three times the measured levels of the cases below
# (production and debug builds measure alike: the twin's arithmetic is its own).
# Moving past one is a real change in the FP32 path's error, not noise.

CEIL = {"pz": 2e-6, "theta": 3e-4, "phasor": 3e-5, "source": 1e-4, "field": 1e-4}
GUARD_FLOOR = 256.0
RENORM_CEIL = 1.5


def ok(label, value, check, good):
    global FAILED
    tag = "ok" if good else "FAIL"
    print(f"--- {label}: {value} (check {check})  {tag}")
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
        extra='  global%fp32_check = "lockstep"\n'))
    s = summary(wd, "fp32ss")
    for q in ("pz", "theta", "phasor", "source", "field"):
        ok(f"lockstep steady worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.0e}", s[q] <= CEIL[q])
    ok("lockstep steady guard ulps", f"{s['guard']:.3e}", f">= {GUARD_FLOOR:.0f}",
       s["guard"] >= GUARD_FLOOR)
    ok("bucket renormalization round trip [ulp]", f"{s['renorm']:.2f}", f"<= {RENORM_CEIL}",
       s["renorm"] <= RENORM_CEIL)

    # 2. Read-only: the instrumented run's FP64 output is the uninstrumented run's.
    run(args.exe, wd, "fp32_off.in", BASE.format(root="fp32off", extra=""))
    same = (wd / "fp32ss.diag.txt").read_bytes() == (wd / "fp32off.diag.txt").read_bytes()
    ok("FP64 diag byte-identical, instrument on vs off", same, "True", same)

    # 3. Falsifiable: the mutation must move the theta level and the field rows.
    run(args.exe, wd, "fp32_mut.in", BASE.format(root="fp32mut",
        extra='  global%fp32_check = "lockstep"\n  global%fp32_mutate = T\n'))
    m = summary(wd, "fp32mut")
    ok("mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"])
    ok("mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"])

    # 4. Compounding: freerun past lockstep, and the end-to-end number exists. The
    # freerun end-to-end on this gain-regime segment is the device-port yardstick.
    run(args.exe, wd, "fp32_fr.in", BASE.format(root="fp32fr",
        extra='  global%fp32_check = "freerun"\n'))
    f = summary(wd, "fp32fr")
    ok("freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"])
    ok("freerun end-to-end power divergence recorded", f"{f['e2e_power']:.3e}",
       "<= 1e-2 on the gain segment", 0 < f["e2e_power"] <= 1e-2)

    # 5. Time dependence, and thread identity of the instrument's own stream.
    run(args.exe, wd, "fp32_td.in", BASE.format(root="fp32td",
        extra=TD_EXTRA + '  global%fp32_check = "lockstep"\n'), threads="8")
    t8 = (wd / "fp32td.fp32.txt").read_bytes()
    t = summary(wd, "fp32td")
    for q in ("pz", "theta", "phasor", "source", "field"):
        ok(f"lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.0e}", t[q] <= CEIL[q])
    ok("lockstep TD guard ulps", f"{t['guard']:.3e}", f">= {GUARD_FLOOR:.0f}",
       t["guard"] >= GUARD_FLOOR)
    run(args.exe, wd, "fp32_td.in", BASE.format(root="fp32td",
        extra=TD_EXTRA + '  global%fp32_check = "lockstep"\n'), threads="1")
    same = t8 == (wd / "fp32td.fp32.txt").read_bytes()
    ok("instrument stream byte-identical, 1 vs 8 threads", same, "True", same)

    # 6. Refused by name outside the twin's coverage.
    r = run(args.exe, wd, "fp32_ref.in", BASE.format(root="fp32ref",
        extra='  global%fp32_check = "lockstep"\n  wake_on = T\n  wake_radius = 2.5e-3\n'
              '  wake_conductivity = 5.813e7\n  wake_relaxation = 8.1e-6\n'),
        expect_fail=True)
    refused = r.returncode != 0 and "FP32_CHECK DOES NOT COVER WAKES" in r.stdout
    ok("wake with the instrument refused by name", refused, "True", refused)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
