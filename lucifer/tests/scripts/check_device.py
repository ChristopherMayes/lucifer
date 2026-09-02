#!/usr/bin/env python3
"""Device backend checks (doc/validation.md): the Metal backend's divergence from the
FP64 path is a recorded number judged by the lockstep instrument with the device in
the twin's role, the instrument cannot touch the FP64 physics, the check can fail,
and everything the backend does not cover is refused by name.

1. Lockstep levels. A steady one-segment run with device = "metal" and fp32_check =
   "lockstep" must land its worst per-quantity divergence inside the recorded device
   ceilings, set at three times the CPU twin's recorded levels where those exist and
   at three times the device's own first measurement for the transverse rows (the CPU
   twin refreshed those from FP64; the device runs its own FP32 maps). The stream
   must carry the exact-wrap verdict: bucket shifts on the fixed-point phase are
   modular arithmetic and are asserted exactly on the device, never to a tolerance.
2. Read-only. The instrumented run must reproduce the uninstrumented run's diag file
   byte for byte: the device twin observes the FP64 path and never steers it.
3. Falsifiable. fp32_mutate = T perturbs a kernel constant (the detuning resonance,
   inside the shader) and coarsens the uploaded phase, and the theta and source rows
   must both move, so a wrong kernel cannot hide behind a check that always passes.
4. End to end. freerun with the device carries the resident state across steps: the
   phasor must compound past lockstep, and the exit-power divergence lands at the
   freerun-measured order. The resident production run (fp32_check off) must land in
   the same band against the CPU run's diag rows, so the instrumented twin and the
   production path corroborate each other.
5. Time dependence. An 8-slice shot-noise window with slippage rotating the resident
   record holds the same ceilings, read-only proof and production band.
6. Refused by name. Wakes, migration, spontaneous radiation, harmonic field sets,
   the unaveraged mode, an unsupported grid (the message names the nearest supported
   size) and an unknown backend name each stop the run.

No device output is asserted byte-identical against another device run: the deposit
accumulates with atomic adds whose ordering is not fixed, so two runs differ in the
source's last bit or two, the reference backends' own documented behavior. The
ceilings absorb it; the read-only proof compares FP64 outputs only.

Usage: check_device.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis_1seg.bmad and aramis_1seg_unavg.bmad. Exit 0 only if
all pass.
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
  grid_n_pts = 64
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

DEV = '  global%device = "metal"\n'

# The recorded ceilings. pz through field are three times the CPU twin's recorded
# levels (doc/validation.md's lockstep table), absorbing that the device's sincos is
# not libm and its fused ops round differently; x through py are three times the
# device's own first recorded levels (2.8e-7 worst, steady and TD, both builds),
# since the CPU twin's transverse rows priced representation only and the device
# runs its own FP32 maps. Moving past one is a real change in the kernels' error.

CEIL = {"x": 1e-6, "px": 1e-6, "y": 1e-6, "py": 1e-6,
        "pz": 1.3e-6, "theta": 1.5e-4, "phasor": 2.2e-5, "source": 2.2e-4, "field": 2.2e-4}

# The guard is the median per-step phase increment in ticks of 2 pi / 2^32; the
# recorded minimum is 6.8e6 and the floor sits a factor of several under it. The
# silent failure is this number collapsing toward zero.

GUARD_FLOOR = 1e6

# The end-to-end band: the freerun-measured order on this gain segment (recorded
# 5.1e-4 on exit power, both builds), ceiling at three times that. The production
# run must land in the same band, or the two roles disagree and that is a finding.

E2E_CEIL = 1.6e-3


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
        elif line.startswith("# device_wrap_exact"):
            vals["wrap"] = line.split()[2]
    return vals


def diag_power(wd, root, nslice):
    rows = [l.split() for l in (wd / f"{root}.diag.txt").read_text().splitlines()
            if not l.startswith("#")]
    return sum(float(r[2]) for r in rows[-nslice:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)

    # 1. Steady lockstep: the device in the twin's role, the recorded ceilings, the
    # exact-wrap verdict from the device's own arithmetic.
    run(args.exe, wd, "dev_ss.in", BASE.format(root="devss",
        extra=DEV + '  global%fp32_check = "lockstep"\n'))
    s = summary(wd, "devss")
    for q in CEIL:
        ok(f"device lockstep steady worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q])
    ok("device lockstep steady guard [ticks]", f"{s['guard']:.3e}", f">= {GUARD_FLOOR:.0e}",
       s["guard"] >= GUARD_FLOOR)
    ok("bucket wraps asserted exact on device", s.get("wrap", "absent"), "T", s.get("wrap") == "T")

    # 2. Read-only: the instrumented run's FP64 output is the uninstrumented run's.
    run(args.exe, wd, "dev_off.in", BASE.format(root="devoff", extra=""))
    same = (wd / "devss.diag.txt").read_bytes() == (wd / "devoff.diag.txt").read_bytes()
    ok("FP64 diag byte-identical, device twin on vs off", same, "True", same)

    # 3. Falsifiable: the perturbed kernel constant and the coarsened upload must
    # move the recorded levels.
    run(args.exe, wd, "dev_mut.in", BASE.format(root="devmut",
        extra=DEV + '  global%fp32_check = "lockstep"\n  global%fp32_mutate = T\n'))
    m = summary(wd, "devmut")
    ok("mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"])
    ok("mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"])

    # 4. End to end, both roles. The freerun twin carries the resident state and
    # compounds; the production run replaces the CPU step outright. Both land in
    # the recorded band and therefore corroborate each other.
    run(args.exe, wd, "dev_fr.in", BASE.format(root="devfr",
        extra=DEV + '  global%fp32_check = "freerun"\n'))
    f = summary(wd, "devfr")
    ok("device freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"])
    ok("device freerun end-to-end power", f"{f['e2e_power']:.3e}",
       f"in (0, {E2E_CEIL:.1e}]", 0 < f["e2e_power"] <= E2E_CEIL)

    run(args.exe, wd, "dev_pr.in", BASE.format(root="devpr", extra=DEV))
    pc = diag_power(wd, "devoff", 1)
    pd = diag_power(wd, "devpr", 1)
    rel = abs(pd - pc) / pc
    ok("device production exit power vs CPU", f"{rel:.3e}",
       f"in (0, {E2E_CEIL:.1e}]", 0 < rel <= E2E_CEIL)

    # 5. Time dependence: slippage rotates the resident record. Same ceilings,
    # same read-only proof, same production band, at 8 threads.
    run(args.exe, wd, "dev_td.in", BASE.format(root="devtd",
        extra=TD_EXTRA + DEV + '  global%fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devtd")
    for q in CEIL:
        ok(f"device lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q])
    ok("device lockstep TD guard [ticks]", f"{t['guard']:.3e}", f">= {GUARD_FLOOR:.0e}",
       t["guard"] >= GUARD_FLOOR)
    run(args.exe, wd, "dev_tdoff.in", BASE.format(root="devtdoff", extra=TD_EXTRA), threads="8")
    same = (wd / "devtd.diag.txt").read_bytes() == (wd / "devtdoff.diag.txt").read_bytes()
    ok("TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "dev_tdpr.in", BASE.format(root="devtdpr", extra=TD_EXTRA + DEV), threads="8")
    pc = diag_power(wd, "devtdoff", 8)
    pd = diag_power(wd, "devtdpr", 8)
    rel = abs(pd - pc) / pc
    ok("device production TD window power vs CPU", f"{rel:.3e}",
       f"<= {E2E_CEIL:.1e}", rel <= E2E_CEIL)

    # 6. Refused by name, never a quiet CPU fallback.
    refusals = [
        ("dev_rw.in", DEV + '  wake_on = T\n  wake_radius = 2.5e-3\n'
                            '  wake_conductivity = 5.813e7\n  wake_relaxation = 8.1e-6\n',
         "DOES NOT COVER WAKES"),
        ("dev_rm.in", DEV + '  global%migrate = T\n', "DOES NOT COVER SLICE MIGRATION"),
        ("dev_rr.in", DEV + '  bmad_com%radiation_damping_on = T\n',
         "DOES NOT COVER SPONTANEOUS RADIATION"),
        ("dev_rh.in", DEV + '  harmonics(2) = 3\n', "DOES NOT COVER HARMONIC FIELD SETS"),
        ("dev_rn.in", '  global%device = "cuda"\n', "UNRECOGNIZED DEVICE"),
    ]
    for name, extra, msg in refusals:
        r = run(args.exe, wd, name, BASE.format(root=name.removesuffix(".in"), extra=extra),
                expect_fail=True)
        refused = r.returncode != 0 and msg in r.stdout
        ok(f"refused by name: {msg.lower()}", refused, "True", refused)

    r = run(args.exe, wd, "dev_rg.in",
            BASE.format(root="dev_rg", extra=DEV).replace("grid_n_pts = 64", "grid_n_pts = 65"),
            expect_fail=True)
    refused = r.returncode != 0 and "nearest supported size is 64" in r.stdout
    ok("refused by name: unsupported grid, nearest size named", refused, "True", refused)

    r = run(args.exe, wd, "dev_ru.in",
            BASE.format(root="dev_ru", extra=DEV).replace("aramis_1seg.bmad", "aramis_1seg_unavg.bmad"),
            expect_fail=True)
    refused = r.returncode != 0 and "DOES NOT COVER THE UNAVERAGED MODE" in r.stdout
    ok("refused by name: the unaveraged mode", refused, "True", refused)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
