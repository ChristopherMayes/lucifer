#!/usr/bin/env python3
"""Device backend checks (doc/validation.md): the Metal backend's divergence from the
FP64 path is a recorded number judged by the lockstep instrument with the device in
the twin's role, the instrument cannot touch the FP64 physics, the check can fail,
and everything the backend does not cover is refused.

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
6. Refused. Wakes, spontaneous radiation, harmonics together with two live
   polarizations, the unaveraged mode, an unsupported grid (the message names the
   nearest supported size) and an unknown backend name each stop the run.
7. The field set. Harmonic members ride the device (planar segment, harmonics 1 and 3)
   and so do two polarization planes (the crossed undulator of check_two_polarization),
   each judged the same four ways: the lockstep rows inside the recorded ceilings with
   the per-member footer lines the set adds, the read-only proof, the mutation moving
   the levels, and the production run against the CPU on the fundamental and on the
   harmonic or the second plane. The harmonic's own identity check is the one-step
   dark deposit's P3/P1 against the Bessel closed form from the dumped particles
   (check_harmonics' identity, its helper shared), which the device holds in FP32. A
   time-dependent Genesis run with shot noise and a third-harmonic field is imported
   on both the CPU and the device at grid 64, the same particles and the same field,
   so both are judged against one reference. One measurement is recorded rather than
   asserted: a quiet-start dark harmonic whose true bunching sits at ~1e-8 of the
   charge (the planar steady-state tier deck) radiates the device's FP32 phase floor
   instead, ~1e-7 of the charge in the phasor row, and its P3 lands two decades above
   the CPU's. The device resolves harmonic bunching down to that floor and not below,
   which the phasor row states per member. The harmonic checks here sit above it.
8. Slice migration, from check_migration's own decks at the device's grid. Migration
   needs no kernel: at an element's last step the walk reads the beam back and releases
   residency before re-slicing, so conservation and phase continuity hold on the device
   exactly as on the CPU, and the check proves the order is right. The production run
   lands in a band against the CPU on the heavy-migration window. The no-op check cannot
   be byte identity on the device (the deposit's atomics), so it is self-referenced:
   zero moves, and the migrate = T run within three times the flutter between two
   migrate = F device runs. Lockstep on the migrating window holds the ceilings and the
   read-only proof, the twin re-staging each step across the changed fills. A slice
   outgrowing the setup rectangle exercises the seam's growth-only resize, which the
   log records.

No device output is asserted byte-identical against another device run: the deposit
accumulates with atomic adds whose ordering is not fixed, so two runs differ in the
source's last bit or two, the reference backends' own documented behavior. The
ceilings absorb it; the read-only proof compares FP64 outputs only.

Usage: check_device.py --exe <lucifer> --workdir <dir> --latdir <tests/bmad>
                       --genesis <genesis4> --pyrepo <openPMD-beamphysics>
The workdir must hold aramis_1seg.bmad and aramis_1seg_unavg.bmad. Exit 0 only if
all pass.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import subprocess
import sys

import h5py
import numpy as np

import beamio
import check_harmonics as ch
import check_migration as cm
import check_two_polarization as tp
import convert_genesis
import fieldio
from nml import to_groups
from read_stats import read_stats

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

# The field set's own levels, each at three times its first measurement (debug build,
# this machine, the check's own first run). The harmonic phasor row is the device's
# harmonic bunching floor, measured 6.6e-8 steady and 2.0e-7 time dependent, three
# times the fundamental's phase error as h = 3 predicts. The harmonic's source and
# field rows are normalized by the fundamental's field and measured 3.5e-9 steady and
# 2.7e-6 time dependent, the latter on a shot-noise window where the third harmonic is
# a quarter of the fundamental's power. The harmonic end-to-end band is the
# third-harmonic power against the CPU on a strongly bunched beam (measured 6.0e-4, the
# fundamental's phase error tripled and squared into a power), and the Bessel identity
# is the one-step P3/P1 in FP32 (measured 4.4e-6 against 2.7e-14 on the CPU). The
# isolation level is |ln(Px/P_drift)| through the y set on the device (measured 4.3e-6;
# the CPU's is 1.6e-15, so this is the FP32 price of "the y set only diffracts Ex").
# The second plane's own band is wider than the fundamental's because Py is a
# hundredth of Px on this line and its relative error carries that ratio (measured
# 4.3e-3). Against the one Genesis reference the device lands at 1.5e-5 on P1 and
# 2.3e-5 on P3 where the CPU lands at 8.4e-7 and 1.1e-6, both inside these bands.

CEIL_H = {"phasor_h3": 6e-7, "source_h3": 8e-6, "field_h3": 8e-6}
H_E2E_CEIL = 1.8e-3
BESSEL_CEIL = 1.3e-5
ISO_CEIL = 1.3e-5
PY_E2E_CEIL = 1.3e-2
PYPX_FLOOR = 5e-3           # check_two_polarization's afterburner floor.

# Migration: the heavy-migration window's exit power against the CPU, at three times
# the first measurement, and the no-op flutter floor: two migrate = F device runs may
# be bit-identical on a small dark deck, and a zero flutter tripled is zero, so the
# floor is a level in slice spacings under which FP32 cannot express a difference.
MIG_E2E_CEIL = 3e-3
NOOP_FLOOR = 1e-12

# The planar segment of check_harmonics, since fc(3) is alive there where the Aramis
# segment is helical and couples only the fundamental, at the device's grid.
PLANAR_LAT = ch.BMAD_LAT

HARM_BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 8192
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e9
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  ran_seed = 777
  harmonics = 1, 3
  write_diag = T
{extra}&end
"""


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
    ap.add_argument("--latdir", required=True, help="tests/bmad, for crossed_probe.bmad")
    ap.add_argument("--genesis", required=True)
    ap.add_argument("--pyrepo", help="openPMD-beamphysics checkout for the dump conversion")
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    exe = pathlib.Path(args.exe).resolve()

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

    # 6. Refused, never a quiet CPU fallback.
    refusals = [
        ("dev_rw.in", DEV + '  wake_on = T\n  wake_radius = 2.5e-3\n'
                            '  wake_conductivity = 5.813e7\n  wake_relaxation = 8.1e-6\n',
         "DOES NOT COVER WAKES"),
        ("dev_rr.in", DEV + '  bmad_com%radiation_damping_on = T\n',
         "DOES NOT COVER SPONTANEOUS RADIATION"),
        ("dev_rn.in", '  global%device = "cuda"\n', "UNRECOGNIZED DEVICE"),
    ]
    for name, extra, msg in refusals:
        r = run(args.exe, wd, name, BASE.format(root=name.removesuffix(".in"), extra=extra),
                expect_fail=True)
        refused = r.returncode != 0 and msg in r.stdout
        ok(f"refused: {msg.lower()}", refused, "True", refused)

    r = run(args.exe, wd, "dev_rg.in",
            BASE.format(root="dev_rg", extra=DEV).replace("grid_n_pts = 64", "grid_n_pts = 65"),
            expect_fail=True)
    refused = r.returncode != 0 and "nearest supported size is 64" in r.stdout
    ok("refused: unsupported grid, nearest size named", refused, "True", refused)

    r = run(args.exe, wd, "dev_ru.in",
            BASE.format(root="dev_ru", extra=DEV).replace("aramis_1seg.bmad", "aramis_1seg_unavg.bmad"),
            expect_fail=True)
    refused = r.returncode != 0 and "DOES NOT COVER THE UNAVERAGED MODE" in r.stdout
    ok("refused: the unaveraged mode", refused, "True", refused)

    field_set(args, wd, exe)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0


def field_set(args, wd, exe):
    """The field set on the device: harmonic members, then two polarization planes."""
    (wd / "dev_planar.bmad").write_text(PLANAR_LAT.format(length="3.96"))
    (wd / "dev_planar_mid.bmad").write_text(PLANAR_LAT.format(length="1.98"))
    (wd / "dev_planar_short.bmad").write_text(PLANAR_LAT.format(length="0.045"))

    # 7a. Harmonics on a strongly bunched beam: check_harmonics' buncher (1 GW over
    # 1.98 m drives the pendulum nonlinear, so b3 is real and far above the device's
    # floor). Lockstep rows and the per-member footer, read-only, mutation, freerun,
    # and the production run against the CPU on P1 and P3.
    print("== the field set: harmonics 1 and 3 on the device ==")
    hb = HARM_BASE.replace('lat_file = "{lat}"', 'lat_file = "dev_planar_mid.bmad"')
    run(args.exe, wd, "devh_ss.in", hb.format(root="devhss", extra=DEV + '  global%fp32_check = "lockstep"\n'))
    s = summary(wd, "devhss")
    for q in CEIL:
        ok(f"harmonic lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q])
    for q in CEIL_H:
        ok(f"harmonic lockstep worst_{q} (per-member footer)", f"{s.get(q, float('nan')):.3e}",
           f"<= {CEIL_H[q]:.1e}", q in s and s[q] <= CEIL_H[q])
    ok("harmonic footer carries the fundamental too", "phasor_h1" in s, "True", "phasor_h1" in s)
    run(args.exe, wd, "devh_off.in", hb.format(root="devhoff", extra=""))
    same = (wd / "devhss.diag.txt").read_bytes() == (wd / "devhoff.diag.txt").read_bytes()
    ok("harmonic FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devh_mut.in", hb.format(root="devhmut",
        extra=DEV + '  global%fp32_check = "lockstep"\n  global%fp32_mutate = T\n'))
    m = summary(wd, "devhmut")
    ok("harmonic mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"])
    ok("harmonic mutation moves worst_source_h3", f"{m['source_h3']:.3e} vs {s['source_h3']:.3e}",
       ">= 5x", m["source_h3"] >= 5 * s["source_h3"])
    run(args.exe, wd, "devh_fr.in", hb.format(root="devhfr", extra=DEV + '  global%fp32_check = "freerun"\n'))
    f = summary(wd, "devhfr")
    ok("harmonic freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"])
    run(args.exe, wd, "devh_pr.in", hb.format(root="devhpr", extra=DEV))
    c1, c3 = harm_powers(wd, "devhoff")
    d1, d3 = harm_powers(wd, "devhpr")
    ok("harmonic production P1 vs CPU (bunched beam)", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL)
    ok("harmonic production P3 vs CPU (bunched beam)", f"{abs(d3-c3)/c3:.3e}", f"<= {H_E2E_CEIL:.1e}",
       abs(d3 - c3) / c3 <= H_E2E_CEIL)

    # 7b. The Bessel identity in FP32: a one-step dark restart from the CPU buncher's
    # exit beam, P3/P1 against the closed form from the dumped particles.
    shutil.copy(wd / "devhoff-final.wf.h5", wd / "devh_dark.wf.h5")
    with h5py.File(wd / "devh_dark.wf.h5", "r+") as h5:
        h5[fieldio.MESH_PATH + "/x"][...] = 0.0
    dep = ch.NML_IMPORT.format(lat="dev_planar_short.bmad", root="devhdep", beam="devhoff-final.beam.h5",
                               field="devh_dark.wf.h5", extra='  transport_model = "bmad"\n' + DEV)
    (wd / "devh_dep.in").write_text(to_groups(dep))
    r = subprocess.run([str(exe), "devh_dep.in"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: devh_dep.in exited {r.returncode}:\n{r.stdout[-3000:]}")
        sys.exit(1)
    expect = ch.bessel_p3_over_p1(wd / "devhdep-final.beam.h5", wd / "devh_dark.wf.h5")
    p1, p3 = harm_powers(wd, "devhdep")
    dv = abs(p3 / p1 / expect - 1)
    ok("device one-step deposit P3/P1 vs the Bessel closed form", f"{dv:.3e}", f"<= {BESSEL_CEIL:.1e}",
       dv <= BESSEL_CEIL)

    # 7c. Time dependence with the set: the harmonics check's TD deck (8 slices, shot
    # noise, dark third harmonic) at the device's grid, slippage rotating every plane.
    # Shot noise puts b3 at ~1/sqrt(N), far above the floor, so P3 is a real comparison.
    td = ch.NML_TD.replace("grid_n_pts = 63", "grid_n_pts = 64").replace('transport_model = "genesis"',
                                                                            'transport_model = "bmad"')
    run(args.exe, wd, "devh_td.in", td.format(lat="dev_planar.bmad", root="devhtd",
        extra=DEV + '  global%fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devhtd")
    for q in CEIL:
        ok(f"harmonic lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q])
    for q in CEIL_H:
        ok(f"harmonic lockstep TD worst_{q}", f"{t.get(q, float('nan')):.3e}", f"<= {CEIL_H[q]:.1e}",
           q in t and t[q] <= CEIL_H[q])
    run(args.exe, wd, "devh_tdoff.in", td.format(lat="dev_planar.bmad", root="devhtdoff", extra=""), threads="8")
    same = (wd / "devhtd.diag.txt").read_bytes() == (wd / "devhtdoff.diag.txt").read_bytes()
    ok("harmonic TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devh_tdpr.in", td.format(lat="dev_planar.bmad", root="devhtdpr", extra=DEV), threads="8")
    c1, c3 = harm_powers(wd, "devhtdoff", window=True)
    d1, d3 = harm_powers(wd, "devhtdpr", window=True)
    ok("harmonic production TD window P1 vs CPU", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL)
    ok("harmonic production TD window P3 vs CPU (shot-noise b3)", f"{abs(d3-c3)/c3:.3e}",
       f"<= {H_E2E_CEIL:.1e}", abs(d3 - c3) / c3 <= H_E2E_CEIL)

    # 7d. One Genesis reference for both: check_harmonics' planar deck made time
    # dependent (8 slices, shot noise on, so the third harmonic's bunching is real on
    # every side) at the device's grid, its starting state imported by the CPU and the
    # device alike. The CPU's level is the anchor the device's is read against.
    genesis_harmonics(args, wd, exe)

    # 7e. The quiet-start floor, recorded: the planar SS tier deck at grid 64, dark
    # third harmonic, no shot noise. The fundamental holds its band. The harmonic's
    # true bunching is ~1e-8 of the charge, below the device's FP32 phase floor, so
    # its P3 is the floor radiating and is printed for the record, not asserted.
    qs = ch.NML_TD.replace("grid_n_pts = 63", "grid_n_pts = 64").replace('transport_model = "genesis"',
                                                                            'transport_model = "bmad"')
    qs = tp.ss(qs).replace("beam_init%bunch_charge = 8.0e-15", "beam_init%bunch_charge = 1.000692285594e-15")
    qs = qs.replace("seed_power = 1e4", "seed_power = 5e3")
    run(args.exe, wd, "devh_qs.in", qs.format(lat="dev_planar.bmad", root="devhqs", extra=""))
    run(args.exe, wd, "devh_qsdev.in", qs.format(lat="dev_planar.bmad", root="devhqsdev", extra=DEV))
    c1, c3 = harm_powers(wd, "devhqs")
    d1, d3 = harm_powers(wd, "devhqsdev")
    ok("quiet-start dark harmonic: fundamental P1 vs CPU", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL)
    print(f"--- quiet-start dark harmonic: device P3 / CPU P3 = {d3/c3:.3e} "
          f"(CPU {c3:.3e} W, device {d3:.3e} W): the FP32 floor radiating, recorded")

    # 7f. Two polarizations: the crossed undulator of check_two_polarization, SS and
    # TD, on its own deck at the device's grid.
    print("== the field set: two polarization planes on the device ==")
    (wd / "crossed_probe.bmad").write_bytes((pathlib.Path(args.latdir) / "crossed_probe.bmad").read_bytes())
    (wd / "devp_c.bmad").write_text("call, file = crossed_probe.bmad\nuse, CROSSED\n")
    (wd / "devp_xd.bmad").write_text("call, file = crossed_probe.bmad\nDEQ: pipe, l = 0.60\n"
                                     "XDRIFT: line = (UNDX, D1, DEQ)\nuse, XDRIFT\n")
    pn = tp.ss(tp.NML).replace("grid_n_pts = 63", "grid_n_pts = 64")
    run(args.exe, wd, "devp_ss.in", pn.format(lat="devp_c.bmad", root="devpss",
        extra=DEV + '  global%fp32_check = "lockstep"\n'))
    s = summary(wd, "devpss")
    for q in CEIL:
        ok(f"two-plane lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q])
    run(args.exe, wd, "devp_off.in", pn.format(lat="devp_c.bmad", root="devpoff", extra=""))
    same = (wd / "devpss.diag.txt").read_bytes() == (wd / "devpoff.diag.txt").read_bytes()
    ok("two-plane FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devp_mut.in", pn.format(lat="devp_c.bmad", root="devpmut",
        extra=DEV + '  global%fp32_check = "lockstep"\n  global%fp32_mutate = T\n'))
    m = summary(wd, "devpmut")
    ok("two-plane mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"])
    ok("two-plane mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"])
    run(args.exe, wd, "devp_pr.in", pn.format(lat="devp_c.bmad", root="devppr", extra=DEV))
    run(args.exe, wd, "devp_dr.in", pn.format(lat="devp_xd.bmad", root="devpdr", extra=DEV))
    px_c = tp.dump_power(wd, "devpoff-final.wf.h5", "x")
    py_c = tp.dump_power(wd, "devpoff-final.wf.h5", "y")
    px_d = tp.dump_power(wd, "devppr-final.wf.h5", "x")
    py_d = tp.dump_power(wd, "devppr-final.wf.h5", "y")
    p_drift = tp.dump_power(wd, "devpdr-final.wf.h5", "x")
    ok("crossed production Px vs CPU", f"{abs(px_d-px_c)/px_c:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(px_d - px_c) / px_c <= E2E_CEIL)
    ok("crossed production Py vs CPU", f"{abs(py_d-py_c)/py_c:.3e}", f"<= {PY_E2E_CEIL:.1e}",
       abs(py_d - py_c) / py_c <= PY_E2E_CEIL)
    iso = abs(np.log(px_d / p_drift))
    ok("crossed device: x-field isolation through the y set, |ln(Px/P_drift)|", f"{iso:.3e}",
       f"<= {ISO_CEIL:.1e}", iso <= ISO_CEIL)
    ok("crossed device: afterburner floor, Py/Px", f"{py_d/px_d:.3e}", f">= {PYPX_FLOOR:.0e}",
       py_d / px_d >= PYPX_FLOOR)
    tdn = tp.NML.replace("grid_n_pts = 63", "grid_n_pts = 64")
    run(args.exe, wd, "devp_td.in", tdn.format(lat="devp_c.bmad", root="devptd",
        extra=DEV + '  global%fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devptd")
    for q in CEIL:
        ok(f"two-plane lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q])
    run(args.exe, wd, "devp_tdoff.in", tdn.format(lat="devp_c.bmad", root="devptdoff", extra=""), threads="8")
    same = (wd / "devptd.diag.txt").read_bytes() == (wd / "devptdoff.diag.txt").read_bytes()
    ok("two-plane TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devp_tdpr.in", tdn.format(lat="devp_c.bmad", root="devptdpr", extra=DEV), threads="8")
    for comp, band in (("x", E2E_CEIL), ("y", E2E_CEIL)):
        pc = tp.dump_power(wd, "devptdoff-final.wf.h5", comp)
        pd = tp.dump_power(wd, "devptdpr-final.wf.h5", comp)
        ok(f"crossed TD production P{comp} vs CPU (slippage live)", f"{abs(pd-pc)/pc:.3e}",
           f"<= {band:.1e}", abs(pd - pc) / pc <= band)

    migration(args, wd, exe)

    # 7g. The combination is refused for the device as for the CPU.
    r = run(args.exe, wd, "devp_rh.in", pn.format(lat="devp_c.bmad", root="devprh",
            extra=DEV + '  harmonics = 1, 3\n'), expect_fail=True)
    refused = r.returncode != 0 and "TWO LIVE POLARIZATIONS" in r.stdout
    ok("refused: harmonics together with two live polarizations", refused, "True", refused)


def migration(args, wd, exe):
    """Slice migration with the device, on check_migration's decks at grid 64."""
    print("== slice migration with the device ==")
    base = cm.BASE.replace("grid_n_pts = 65", "grid_n_pts = 64").replace("&end\n", "{extra}&end\n")
    # check_migration's heavy deck is a dark quiet start, which suits its accounting and
    # not the instrument: a slice whose field sits at FP64 roundoff normalizes the source
    # and field rows to nothing (the harmonic floor's lesson, doc/validation.md), and
    # every slippage feeds the tail slice exactly such a field. A seed fills the window
    # and shot noise fills each fresh tail slice, as the TD lockstep deck above does, and
    # the 60 m_e c^2 energy spread migrates just as heavily.
    heavy_base = base.replace("seed_power = 0\n",
                              "seed_power = 1e4\n  seed_waist_size = 30e-6\n  shot_noise = T\n")
    heavy = dict(lat="aramis.bmad", sig_pz="5.282703940115e-03", emit="4e-7", npart=1024,
                 slen="4.8e-9", q="4.803322970853e-14", half="2.4e-9")
    frozen = dict(lat="aramis_1seg.bmad", sig_pz="8.804506566858e-08", emit="1e-13", npart=256,
                  slen="1.2e-9", q="1.200830742713e-14", half="6e-10")

    # 8a. Heavy migration on the device: conservation, phase continuity and window
    # residency by check_migration's own readers, the production band against the CPU,
    # and the capacity growth in the log.
    run(args.exe, wd, "devm_c.in", heavy_base.format(root="devmc", mig="T", extra="", **heavy), threads="8")
    r = run(args.exe, wd, "devm_d.in", heavy_base.format(root="devmd", mig="T", extra=DEV, **heavy), threads="8")
    moved, bdev, drops = cm.read_migration_file(wd, "devmd")
    z, q_win = cm.in_window_charge(wd / "devmd.diag.txt")
    q0 = q_win[0]
    worst = max(abs(q_win[i] + sum(q for zd, q in drops if zd <= zr + 1e-9) - q0) / q0
                for i, zr in enumerate(z))
    ok("migration on the device bites", f"{moved} moves, {sum(q for _, q in drops):.3e} C dropped",
       "> 10000 moves and drops", moved > 10000 and sum(q for _, q in drops) > 0)
    ok("migration conservation on the device, worst violation", f"{worst:.2e}", "< 1e-10", worst < 1e-10)
    ok("migration phase continuity on the device, worst deviation", f"{bdev:.2e}", "< 1e-10", bdev < 1e-10)
    slen = 2 * np.pi * 3
    n_out = sum(int(np.sum((sl["theta"] < -1e-9) | (sl["theta"] >= slen + 1e-9)))
                for sl in beamio.read_slices(wd / "devmd-final.beam.h5", 1e-10, 3e-10))
    ok("migration window residency on the device", f"{n_out} outside", "0", n_out == 0)
    nslice = 16
    pc = diag_power(wd, "devmc", nslice)
    pd = diag_power(wd, "devmd", nslice)
    ok("heavy-migration window power, device vs CPU", f"{abs(pd-pc)/pc:.3e}", f"<= {MIG_E2E_CEIL:.1e}",
       abs(pd - pc) / pc <= MIG_E2E_CEIL)
    grown = "particle capacity grown" in r.stdout
    ok("a migrated slice outgrew the rectangle and the seam grew it", grown, "True", grown)

    # 8b. Lockstep across migration: the twin re-stages every step over the changed
    # fills, inside the ceilings, and the FP64 path is untouched.
    run(args.exe, wd, "devm_l.in", heavy_base.format(root="devml", mig="T",
        extra=DEV + '  global%fp32_check = "lockstep"\n', **heavy), threads="8")
    s = summary(wd, "devml")
    for q in CEIL:
        ok(f"migrating lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q])
    same = (wd / "devml.diag.txt").read_bytes() == (wd / "devmc.diag.txt").read_bytes()
    ok("migrating FP64 diag byte-identical, device twin on vs off", same, "True", same)

    # 8c. The no-op, self-referenced against the device's own flutter. The frozen deck
    # is dark, so its power is FP32 noise and says nothing; the read-back beam does. With
    # zero moves the migrate = T run differs from migrate = F only in where the readback
    # happened, so the final beam dumps must agree to within the flutter two migrate = F
    # device runs show between themselves, measured on z in units of the slice spacing.
    run(args.exe, wd, "devm_f1.in", base.format(root="devmf1", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_f2.in", base.format(root="devmf2", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_t.in", base.format(root="devmt", mig="T", extra=DEV, **frozen), threads="8")
    moved = cm.read_migration_file(wd, "devmt")[0]
    z1 = beam_z(wd / "devmf1-final.beam.h5")
    z2 = beam_z(wd / "devmf2-final.beam.h5")
    zt = beam_z(wd / "devmt-final.beam.h5")
    flutter = float(np.max(np.abs(z1 - z2))) / 3e-10
    tol = max(3 * flutter, NOOP_FLOOR)
    dev_t = float(np.max(np.abs(zt - z1))) / 3e-10
    ok("migration no-op on the device: moves", moved, "0", moved == 0)
    ok("migration no-op on the device: migrate T vs F final beam z [spacings]", f"{dev_t:.2e}",
       f"<= {tol:.1e} (3 x flutter {flutter:.1e})", dev_t <= tol)


def beam_z(path):
    """Every particle's longitudinal position c*t [m] from an openPMD beam dump, in the
    file's own order. Bmad's dumps hold z as a constant and the coordinate in time."""
    with h5py.File(path) as h5:
        rec = h5["data"][list(h5["data"].keys())[0]]["particles"]
        name = list(rec.keys())[0]
        return 2.99792458e8 * np.asarray(rec[name]["time"][...], dtype=float)


def harm_powers(wd, root, window=False):
    """Exit power of the fundamental and the third harmonic from the stats file, the
    window's sum when time dependent."""
    with read_stats(wd / f"{root}.stats.h5") as st:
        p1 = st["field/total/power"][-1]
        p3 = st["field/harm3/total/power"][-1]
    if window:
        return float(np.sum(p1)), float(np.sum(p3))
    return float(p1[0]), float(p3[0])


def genesis_harmonics(args, wd, exe):
    """check_harmonics' planar deck made time dependent, with shot noise, at grid 64: one
    Genesis reference, its starting state imported by the CPU and the device."""
    (wd / "devg_planar.lat").write_text(ch.GENESIS_LAT)
    deck = ch.GENESIS_DECK.replace("ngrid=255", "ngrid=64").replace("shotnoise=0\n", "shotnoise=true\n")
    deck = deck.replace("lattice=planar.lat", "lattice=devg_planar.lat").replace("rootname=H3", "rootname=DEVG")
    deck = deck.replace("&field\npower=5e3", "&time\nslen=8e-10\nsample=1\n&end\n\n&field\npower=5e3")
    deck = deck.replace("H3-initial", "DEVG-initial").replace("H3-final", "DEVG-final")
    (wd / "DEVG.in").write_text(deck)
    r = subprocess.run([args.genesis, "DEVG.in"], cwd=wd, capture_output=True, text=True,
                       env=dict(os.environ, FI_PROVIDER="tcp"))
    if r.returncode != 0:
        print(f"FAIL: genesis exited {r.returncode}:\n{r.stdout[-2000:]}\n{r.stderr[-500:]}")
        sys.exit(1)
    for src, dst in (("DEVG-initial.par.h5", "DEVG-initial.beam.h5"),
                     ("DEVG-initial.fld.h5", "DEVG-initial.wf.h5")):
        convert_genesis.to_openpmd(wd / src, wd / dst, args.pyrepo)
    imp = ch.NML_IMPORT.replace('transport_model = "genesis"', 'transport_model = "bmad"')
    curves = {}
    for tag, extra in (("cpu", ""), ("dev", DEV)):
        name = f"devg_{tag}"
        (wd / (name + ".nml")).write_text(to_groups(imp.format(lat="dev_planar.bmad", root=name,
            beam="DEVG-initial.beam.h5", field="DEVG-initial.wf.h5", extra=extra)))
        r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                           env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
        if r.returncode != 0:
            print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}")
            sys.exit(1)
        with read_stats(wd / f"{name}.stats.h5") as st:
            curves[tag] = (st["field/total/power"][:].sum(axis=1), st["field/harm3/total/power"][:].sum(axis=1))
    with h5py.File(wd / "DEVG.out.h5") as h5:
        p1g = h5["Field/power"][:].sum(axis=1)
        p3g = h5["Field3/power"][:].sum(axis=1)

    def level(p1b, p3b):
        n = min(len(p1g), len(p1b))
        lo = n // 8
        return (float(np.max(np.abs(p1b[:n] - p1g[:n])) / np.max(p1g[:n])),
                float(np.max(np.abs(p3b[lo:n] - p3g[lo:n])) / np.max(p3g[:n])))

    c1, c3 = level(*curves["cpu"])
    d1, d3 = level(*curves["dev"])
    print(f"--- Genesis TD harmonic reference at grid 64: CPU P1 {c1:.3e} P3 {c3:.3e} (the anchor)")
    ok("device vs Genesis, TD planar with shot noise, P1", f"{d1:.3e}", f"<= {E2E_CEIL:.1e}", d1 <= E2E_CEIL)
    ok("device vs Genesis, TD planar with shot noise, P3", f"{d3:.3e}", f"<= {H_E2E_CEIL:.1e}", d3 <= H_E2E_CEIL)


if __name__ == "__main__":
    sys.exit(main())
