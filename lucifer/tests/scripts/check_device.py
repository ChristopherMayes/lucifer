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
   polarizations, two live polarizations inside the unaveraged mode, an unsupported
   grid (the message names the nearest supported size) and an unknown backend name
   each stop the run.
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
9. The unaveraged mode, whose quiver-resolving push, kick and deposit and ledger
   reduction are the kernels this mode adds and whose field solve is the one above.
   Judged the same ways, with the transform pair's own energy loss measured beside
   them because it is what sets its ceilings: that mode runs the pair nsub times a
   record step where the averaged mode runs it once. The frame series rides the comb
   there too, so a frame is written from arrays the readback refreshed, and what it
   must carry is the chart: the undulator quiver sits in px as a common offset, which
   the mean sees and an rms does not (FINDINGS 7.68). See unaveraged() below.

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
import re
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
  source_filter = F
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
  n_wavelength = 3
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
# the first measurement. The no-op flutter floor that sat here is gone with the flutter:
# the fixed-point deposit made the no-op comparison exact, so it needs no level at all.
MIG_E2E_CEIL = 3e-3

# The source filter's device-against-CPU ceiling. The device filters the source by four
# FP32 passes of its own where the CPU multiplies inside the field's FP64 transform pair,
# so the two differ by more than the plain solve does. Measured 1.7e-4 end to end.
SF_E2E_CEIL = 1.0e-3

# The unaveraged mode's own ceilings, each three times its first measurement (this
# machine, both builds, the check's own decks). They are wider than the averaged ones
# and the reason is arithmetic rather than the port: this mode runs the transform pair
# nsub times a record step where the averaged mode runs it once, and an FP32 pair loses
# about 1.3e-7 of the field's energy every time it runs (doc/validation.md's unaveraged
# device section). The transverse rows are the FP32 push's own, so they sit an order
# above the averaged mode's, whose transverse maps are one map a step.
#
# The source column carries this mode's ledger term rather than a source grid, scaled
# by the energy the record step turned over. On a dark-start shot-noise window a slice
# can turn over almost nothing, which widens that row and nothing else.

CEIL_U = {"x": 3.5e-6, "px": 2.5e-5, "y": 3.0e-6, "py": 2.1e-5, "pz": 6.5e-6,
          "theta": 1.1e-5, "phasor": 2.5e-7, "source": 5.0e-3, "field": 2.0e-5}
GUARD_FLOOR_U = 1e9
U_E2E_CEIL = 5.0e-3

# The frame series with the device resident in this mode. The mean px a frame carries is
# the quiver chart's own, and the device's differs from the CPU's by the FP32 round trip
# of the chart rather than by the transform loss that sets the band above: measured
# 1.9e-7, and this sits well under that band so a chart error cannot hide in it.

U_FRAME_CEIL = 1.0e-5

# The end-to-end band against the CPU's own unaveraged run. It is set by the transform
# pair's loss below rather than by anything the port chose: on the steady deck, whose
# seed barely grows over one segment, the exit power sits 1.7e-3 under the CPU's, which
# is 5340 pairs of that loss carried into the power. The time-dependent window measures
# 1.5e-6 because its power comes from a beam that is radiating rather than from a seed
# the arithmetic is slowly eating.

# What one FP32 transform pair does to the field's energy, with no beam to feed it.
# Measured -1.29e-7 a pair on the averaged path and -1.34e-7 on the unaveraged one,
# and it is what sets the ceilings above. The band is wide enough to carry a grid
# change and narrow enough that a pair which stopped being nearly unitary would show.

PAIR_LOSS = (5e-8, 4e-7)
# The header's projection of that loss against what the same run measures. The rate it
# projects from is one machine's fit to four grids, so a factor of two either way.
PAIR_RE = re.compile(r"FP32 transform pair (\d+) times a field record \((\d+) unaveraged "
                     r"substeps\), projected field-energy loss ([-+.0-9Ee]+) a pair")
PROJ_BAND = (0.5, 2.0)

# The planar segment of check_harmonics, since fc(3) is alive there where the Aramis
# segment is helical and couples only the fundamental, at the device's grid.
PLANAR_LAT = ch.BMAD_LAT

HARM_BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
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

    field_set(args, wd, exe)
    source_filter(args, wd, exe)
    unaveraged(args, wd, exe)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0



def unaveraged(args, wd, exe):
    """
    The unaveraged mode on the device (doc/validation.md's unaveraged device section).

    The mode resolves the undulator quiver, so a record step is nsub Strang substeps of
    half push, radiation kick with its deposit, half push and the slice's own diffract
    and source add. Three kernels are the mode's own and the rest of the step is the
    kernels the averaged path already had, the four-pass solve above all.

    The instrument judges it the way it judges the averaged path: the device takes the
    twin's role in the unaveraged lockstep, the nine rows land inside recorded ceilings,
    the FP64 stream is byte identical with the twin on against off, and the mutation
    moves the rows the mutation reaches. What the mode's own twin makes the columns mean
    is in check_unaveraged.py: rows 1 to 4 are the quiver chart and are a worst
    per-particle difference, the source column carries the ledger's kick-side term, and
    the field column carries the twin's own record.

    The transform pair's energy loss is measured here rather than inferred, because it
    is what sets these ceilings: a run with no charge leaves the repeated diffract as the
    only thing acting on the field, and what the field then loses is the pair's own
    arithmetic.
    """
    print("== the unaveraged mode on the device ==")
    latdir = pathlib.Path(args.latdir)
    for f in ("aramis_1seg.bmad", "aramis_1seg_unavg.bmad", "crossed_probe.bmad"):
        (wd / f).write_bytes((latdir / f).read_bytes())
    (wd / "dvu_cross.bmad").write_text(
        "call, file = crossed_probe.bmad\nuse, CROSSED\nwiggler::*[FEL_METHOD] = unaveraged\n")
    lat = "aramis_1seg_unavg.bmad"
    base = BASE.replace("aramis_1seg.bmad", lat)

    # 1. Steady lockstep, the device in the twin's role.
    run(args.exe, wd, "dvu_ss.in", base.format(root="dvuss",
        extra=DEV + '  global%fp32_check = "lockstep"\n'))
    s = summary(wd, "dvuss")
    for q in CEIL_U:
        ok(f"unaveraged lockstep steady worst_{q}", f"{s[q]:.3e}",
           f"<= {CEIL_U[q]:.1e}", s[q] <= CEIL_U[q])
    ok("unaveraged lockstep steady guard [ticks]", f"{s['guard']:.3e}",
       f">= {GUARD_FLOOR_U:.0e}", s["guard"] >= GUARD_FLOOR_U)

    # 2. Read-only. The FP64 run here carries neither the device nor the instrument, so
    # what is compared is the CPU's own unaveraged step with the twin beside it and
    # without.
    run(args.exe, wd, "dvu_off.in", base.format(root="dvuoff", extra=""))
    same = (wd / "dvuss.diag.txt").read_bytes() == (wd / "dvuoff.diag.txt").read_bytes()
    ok("unaveraged FP64 diag byte-identical, device twin on vs off", same, "True", same)
    same = (wd / "dvuss.ledger.txt").read_bytes() == (wd / "dvuoff.ledger.txt").read_bytes()
    ok("unaveraged FP64 ledger byte-identical, device twin on vs off", same, "True", same)

    # 3. Falsifiable. The hook coarsens the uploaded lag accumulator by 65536 ticks
    # (9.6e-5 rad) and the residual angle the kick reads, and the rows those reach must
    # move. The transverse rows are not among them: gamma changes only in the kick, so a
    # phase perturbation reaches x and px only at second order.
    run(args.exe, wd, "dvu_mut.in", base.format(root="dvumut",
        extra=DEV + '  global%fp32_check = "lockstep"\n  global%fp32_mutate = T\n'))
    m = summary(wd, "dvumut")
    for q, factor in (("theta", 10), ("phasor", 50), ("field", 20)):
        ok(f"unaveraged mutation moves worst_{q}", f"{m[q]:.3e} vs {s[q]:.3e}",
           f">= {factor}x", m[q] >= factor * s[q])

    # 4. The production role, against the CPU's own unaveraged run.
    run(args.exe, wd, "dvu_pr.in", base.format(root="dvupr", extra=DEV))
    pc = diag_power(wd, "dvuoff", 1)
    pd = diag_power(wd, "dvupr", 1)
    rel = abs(pd - pc) / pc
    ok("unaveraged production power vs CPU", f"{rel:.3e}", f"<= {U_E2E_CEIL:.1e}",
       rel <= U_E2E_CEIL)

    # 4b. The frame series with the device resident in this mode. The comb's positions
    # are where the resident state comes back to the host, so a frame is written from
    # arrays the device refreshed rather than from stale ones. The chart is the thing to
    # check: a frame taken mid-segment carries the undulator quiver in px, and a
    # statistic that cannot see a common offset would pass with the quiver missing
    # (FINDINGS 7.68), so the mean is compared and not the spread.
    frame_extra = '  global%comb_ds_save = 0.5\n  global%dump_at_comb = T\n'
    run(args.exe, wd, "dvu_fr.in", base.format(root="dvufr", extra=DEV + frame_extra))
    run(args.exe, wd, "dvu_froff.in", base.format(root="dvufroff", extra=frame_extra))
    fd = sorted(wd.glob("dvufr-[0-9]*.beam.h5"))
    fc = sorted(wd.glob("dvufroff-[0-9]*.beam.h5"))
    ok("unaveraged frames on the device, count against the CPU",
       f"{len(fd)} against {len(fc)}", "equal and nonzero",
       len(fd) == len(fc) and len(fd) > 0)
    ok("unaveraged field frames on the device, count against the CPU",
       f"{len(sorted(wd.glob('dvufr-[0-9]*.wf.h5')))} against {len(fd)}", "equal",
       len(sorted(wd.glob("dvufr-[0-9]*.wf.h5"))) == len(fd))

    if fd and len(fd) == len(fc):
        def read(path):
            with h5py.File(path) as h:
                it = list(h["data"])[0]
                m = h.attrs.get("felMethod")
                m = np.ravel(m)[0] if m is not None else None
                if isinstance(m, bytes):
                    m = m.decode()
                return (float(np.ravel(h.attrs["sPosition"])[0]), m,
                        float(np.mean(h[f"data/{it}/particles/electron/momentum/x"])))

        d = [read(f) for f in fd]
        c = [read(f) for f in fc]
        sdev = max(abs(a[0] - b[0]) for a, b in zip(d, c))
        ok("unaveraged frames on the device land at the CPU's own s", f"{sdev:.3e}",
           "<= 1e-9", sdev <= 1e-9)

        # A frame taken outside an FEL element carries no felMethod, which the frame at
        # the line's start is. Inside one, this mode's frames must name the chart they
        # hold, since the averaged physics and the Bmad seam read the two differently.
        inside = [a for a in d if a[1] is not None]
        m = sorted({a[1].lower() for a in inside})
        ok("unaveraged frames on the device name the chart they hold",
           f"{m} over {len(inside)} of {len(d)} frames", "['unaveraged'], some inside",
           m == ["unaveraged"] and len(inside) > 0)

        # The quiver is a common offset of order aw/gamma, so the mean carries it and an
        # rms does not. It reaches 1e5 eV/c mid-segment against 1e2 in the averaged chart.
        quiver = max(abs(a[2]) for a in inside) if inside else 0.0
        ok("unaveraged frames on the device carry the quiver in px", f"{quiver:.3e}",
           "> 1e4", quiver > 1e4)
        scale = max(max(abs(b[2]) for b in c), 1e-300)
        rel = max(abs(a[2] - b[2]) for a, b in zip(d, c)) / scale
        ok("unaveraged frame px mean, device against CPU", f"{rel:.3e}",
           f"<= {U_FRAME_CEIL:.1e}", rel <= U_FRAME_CEIL)

    # 5. A time-dependent window: slippage rotates the resident record between record
    # steps and the ledger's banked terms cross the seam every record step.
    run(args.exe, wd, "dvu_td.in", base.format(root="dvutd",
        extra=DEV + TD_EXTRA + '  global%fp32_check = "lockstep"\n'))
    t = summary(wd, "dvutd")
    for q in CEIL_U:
        ok(f"unaveraged lockstep TD worst_{q}", f"{t[q]:.3e}",
           f"<= {CEIL_U[q]:.1e}", t[q] <= CEIL_U[q])
    run(args.exe, wd, "dvu_tdoff.in", base.format(root="dvutdoff", extra=TD_EXTRA))
    same = (wd / "dvutd.diag.txt").read_bytes() == (wd / "dvutdoff.diag.txt").read_bytes()
    ok("unaveraged TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "dvu_tdpr.in", base.format(root="dvutdpr", extra=DEV + TD_EXTRA))
    pc = diag_power(wd, "dvutdoff", 8)
    pd = diag_power(wd, "dvutdpr", 8)
    rel = abs(pd - pc) / pc
    ok("unaveraged production TD window power vs CPU", f"{rel:.3e}",
       f"<= {U_E2E_CEIL:.1e}", rel <= U_E2E_CEIL)

    # 6. Two device runs of one unaveraged deck agree bit for bit. The deposit
    # accumulates in fixed point and the ledger's spontaneous term reduces without an
    # atomic, so no arrival order reaches either.
    run(args.exe, wd, "dvu_r1.in", base.format(root="dvur1", extra=DEV + TD_EXTRA))
    run(args.exe, wd, "dvu_r2.in", base.format(root="dvur2", extra=DEV + TD_EXTRA),
        threads="8")
    same = ((wd / "dvur1.diag.txt").read_bytes() == (wd / "dvur2.diag.txt").read_bytes()
            and (wd / "dvur1.ledger.txt").read_bytes() == (wd / "dvur2.ledger.txt").read_bytes())
    ok("two unaveraged device runs, diag and ledger identical", same, "True", same)

    # 7. The transform pair's own energy loss, with no charge to feed the field. The
    # CPU's FP64 pair holds the energy to 1.6e-12 over the same 5340 applications. The
    # run header counts the pairs and projects the loss from the recorded rate, and the
    # count must be the deck's and the projection must sit near what the run measured.
    for root, extra, npair, nsub in (("dvu_darku", DEV, 89 * 60, 89 * 60),
                                     ("dvu_darka", DEV, 89, 0)):
        text = base.format(root=root, extra=extra)
        if root.endswith("a"):
            text = text.replace(lat, "aramis_1seg.bmad")
        r = run(args.exe, wd, f"{root}.in",
                text.replace("bunch_charge = 1.000692285594e-15", "bunch_charge = 1e-30"))
        p0 = diag_power(wd, root, 1)
        rows = [l.split() for l in (wd / f"{root}.diag.txt").read_text().splitlines()
                if not l.startswith("#")]
        loss = (float(rows[0][2]) - p0) / float(rows[0][2]) / npair
        mode = "unaveraged" if root.endswith("u") else "averaged"
        good = PAIR_LOSS[0] <= loss <= PAIR_LOSS[1]
        ok(f"FP32 transform pair energy loss, {mode} deck", f"{loss:.3e} a pair",
           f"in [{PAIR_LOSS[0]:.0e}, {PAIR_LOSS[1]:.0e}]", good)
        m = PAIR_RE.search(r.stdout)
        counts = (int(m.group(1)), int(m.group(2))) if m else None
        ok(f"header counts the transform pairs, {mode} deck", counts, f"== {(npair, nsub)}",
           counts == (npair, nsub))
        ratio = float(m.group(3)) / loss if m else float("nan")
        ok(f"header projection over measured loss, {mode} deck", f"{ratio:.3f}",
           f"in [{PROJ_BAND[0]}, {PROJ_BAND[1]}]", PROJ_BAND[0] <= ratio <= PROJ_BAND[1])

    # 8. Refused. Two live polarizations in this mode is the one configuration the CPU
    # carries and the device does not.
    r = run(args.exe, wd, "dvu_rp.in",
            base.format(root="dvurp", extra=DEV).replace(
                'lat_file = "' + lat + '"', 'lat_file = "dvu_cross.bmad"'),
            expect_fail=True)
    refused = r.returncode != 0 and "TWO POLARIZATIONS IN THE UNAVERAGED MODE" in r.stdout
    ok("refused: two polarizations in the unaveraged mode", refused, "True", refused)

def source_filter(args, wd, exe):
    """
    The source filter on the device (fel-physics.md sec-source-filter).

    The CPU adds the filtered source inside the field's own transform pair. The device's
    fused solve adds the source in real space in its last pass, so the source is filtered
    in place first by four passes of its own and the solve is left untouched. The two
    reach the same quantity by a different route, and this is what says so. The
    transcription itself is checked against Genesis4 on the CPU
    (check_source_filter.py), which the device cannot join: Genesis4's reference chain is
    a 255-point grid and the Metal solver takes powers of two.

    The mutation hook is the CPU check's own, so the device refuses it rather than
    carrying a second wrong solve.
    """
    # A 256-point grid rather than the 64 the rest of this file uses. The filter removes
    # emission at angles the grid carries beyond the edge, and a coarse grid carries
    # almost none: at 6.3 um cells the Nyquist angle is 7.9 urad against an edge near 3,
    # so there is nothing to remove and the check could not see a filter that did nothing.
    nslice = 8
    for root, filt, extra in (("dev_sf", "  source_filter = T\n", DEV),
                              ("dev_sfoff", "  source_filter = F\n", DEV),
                              ("cpu_sf", "  source_filter = T\n", "")):
        text = BASE.format(root=root, extra=extra + TD_EXTRA + filt)
        run(args.exe, wd, f"{root}.in", text.replace("grid_n_pts = 64", "grid_n_pts = 256"))

    pd = diag_power(wd, "dev_sf", nslice)
    pc = diag_power(wd, "cpu_sf", nslice)
    pn = diag_power(wd, "dev_sfoff", nslice)
    ok("source filter, device vs CPU window power", f"{abs(pd - pc) / pc:.3e}",
       f"<= {SF_E2E_CEIL:.1e}", abs(pd - pc) / pc <= SF_E2E_CEIL)
    ok("source filter on the device removes wide-angle power", f"{pn / pd:.2f}x", "> 1.5",
       pn / pd > 1.5)

    r = run(args.exe, wd, "dev_sfm.in",
            BASE.format(root="dev_sfm", extra=DEV + TD_EXTRA +
                        "  source_filter = T\n  source_filter_mutate = T\n"),
            expect_fail=True)  # refused before the grid matters
    refused = r.returncode != 0 and "DOES NOT COVER SOURCE_FILTER_MUTATE" in r.stdout
    ok("refused: source_filter_mutate on the device", refused, "True", refused)


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
    pn = tp.ss(tp.NML)
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
    tdn = tp.NML
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
    reproducible(args, wd, exe)
    precision(args, wd, exe)
    quantization(args, wd, exe)

    # 7g. The combination is refused for the device as for the CPU.
    r = run(args.exe, wd, "devp_rh.in", pn.format(lat="devp_c.bmad", root="devprh",
            extra=DEV + '  harmonics = 1, 3\n'), expect_fail=True)
    refused = r.returncode != 0 and "TWO LIVE POLARIZATIONS" in r.stdout
    ok("refused: harmonics together with two live polarizations", refused, "True", refused)


# Three things cannot match and none of them is physics: the resolved parameters and the
# echoed input carry the out_root, which differs by construction, and the timestamp is
# the wall clock. Everything else is compared, the field and beam records among them,
# and the Bmad versions and the lattice source are compared and do match. The timestamp
# was found by this check rather than reasoned about: the first pair of runs landed in
# one second and the second pair did not.
REPRO_SKIP = ("params", "meta/input_echo", "meta/timestamp")


# The four configurations the recorded ceilings do not isolate, each measured as the
# lockstep source row: the device's deposit against the same deposit in FP64. The
# fixed-point accumulator has to be no worse than the float one it replaced on every
# one, and the numbers beside them are that float deposit measured on these same decks.
# Cancelling phases come first because a quiet start leaves the slice phasor at
# roundoff, which is where a coarse quantum would show before anywhere else (7.47).

PRECISION_CASES = (
    ("cancelling phases", "  beamlet_size = 8\n", 1.464e-05),
    ("unequal weights", "  beamlet_size = 8\n  split_weights = T\n", 1.465e-05),
    ("charge in few cells", "  beam_init%a_norm_emit = 4e-9\n"
                            "  beam_init%b_norm_emit = 4e-9\n", 6.889e-06),
    ("the largest load", "  beam_init%n_particle = 32768\n", 2.507e-06),
    # The load doc/performance.md prices the deposit at. It carries no float-deposit
    # number beside it, the float deposit never having been measured here and not being
    # measurable now, so it is held to the recorded ceiling and to the headroom instead.
    # Its purpose is that the accuracy cases reach the load the cost table uses.
    ("the performance load", "  beam_init%n_particle = 131072\n", None),
)

# The float deposit is not reproducible, so its numbers above carry their own scatter and
# the comparison allows five percent over them. Measured against it the fixed-point
# accumulator came in at 1.462e-05, 1.462e-05, 6.948e-06 and 2.506e-06: the same levels,
# two of them slightly better and one nine parts in a thousand worse.
PRECISION_MARGIN = 1.05

# The quantum has to sit this far below the FP32 spacing of the deposit's own bound, in
# bits, on every deck the precision cases run. The scale is derived per run and the run
# prints the headroom it reached, so this is read back rather than assumed. The setup
# refusal requires only 8, which is where fixed point stops beating the float accumulation
# it replaced. 24 is the margin the derivation actually delivers, and
# holding the checks to it means a derivation that quietly lost precision would be caught
# here rather than passing on a combined source row.
PRECISION_BITS_MIN = 24

# What the quantum costs, measured rather than bounded. global%device_dep_mutate shifts
# the derived scale by whole bits and touches nothing else, so coarsening the quantum by a
# known factor and watching the source row is a direct separation of the quantization from
# the phase and arithmetic error beside it. Measured on the cancelling-phase deck: the row
# is 1.4673e-05 at the derived scale and the same to every digit 24 bits coarser, first
# moves at 26 (1.0100x), and reaches 1.0505x at 29, which is the last shift the headroom
# refusal allows. So the 8-bit floor sits where quantization starts to cost a few percent,
# and at the 37 bits the derivation delivers it costs nothing that the row can see.
QUANT_SHIFT = -24          # 16 million times the quantum, and still nothing
QUANT_UNMOVED = 1.001      # the row may not move by even a part in a thousand there
QUANT_REFUSED = -30        # one bit past the floor, where the run must be refused
BITS_RE = re.compile(r"(\d+) bits below the FP32 spacing")


def precision(args, wd, exe):
    """
    What the fixed-point deposit costs in accuracy, which is nothing.

    The lockstep source row is the device's deposit against the same deposit carried in
    FP64 on the host, so it prices the accumulator directly rather than through a power.
    Each case is asserted twice: under the recorded source ceiling, and no worse than the
    float accumulator measured on the same deck, which is the number carried beside it.

    A fixed-point accumulator can be reproducible and still be worse, by taking a scale
    that keeps the bound but quantizes away what matters. Thirty-two bits would have done
    exactly that here, losing to the float deposit by a factor of seven. These four say
    the sixty-four-bit accumulator did not.

    The source row is a combined number, carrying the FP32 arithmetic and the phase error
    along with the quantization, so on its own it cannot say which of the three moved. The
    quantization is separated here by asserting the property it depends on directly, per
    deck: the run derives its own scale and prints how far the quantum sits below the FP32
    spacing of the bound it was chosen against, and that headroom is read back and required
    to be at least PRECISION_BITS_MIN. A quantum that far under the spacing cannot be what
    a source row of 1e-5 is made of. The two together are the separation: the headroom
    bounds the quantization per deck, and the source row prices everything at once.
    """
    print("== 10. the deposit's accuracy against the CPU's FP64 deposit ==")
    for name, extra, before in PRECISION_CASES:
        root = "devprec" + name.split()[0][:4]
        r = run(exe, wd, f"{root}.in", BASE.format(
            root=root, extra=DEV + '  global%fp32_check = "lockstep"\n' + extra))
        if r.returncode != 0:
            ok(f"deposit accuracy, {name}", "run failed", "exit 0", False)
            continue
        v = summary(wd, root)["source"]
        if before is None:
            ok(f"deposit accuracy, {name}", f"{v:.3e}",
               f"<= {CEIL['source']:.1e}", v <= CEIL["source"])
        else:
            ok(f"deposit accuracy, {name}", f"{v:.3e}",
               f"<= {CEIL['source']:.1e} and <= the float deposit's {before:.2e}",
               v <= CEIL["source"] and v <= before * PRECISION_MARGIN)
        m = BITS_RE.search(r.stdout)
        bits = int(m.group(1)) if m else -1
        ok(f"deposit quantum below the FP32 spacing of its bound, {name}",
           f"{bits} bits" if m else "the run printed no headroom",
           f">= {PRECISION_BITS_MIN} bits", bits >= PRECISION_BITS_MIN)


def quantization(args, wd, exe):
    """
    What the quantum contributes to the measured deposit error, which is nothing.

    The source row prices the phase error, the FP32 arithmetic and the quantization
    together, so on its own it cannot say which of the three it is made of. Shifting the
    deposit's scale moves one of them and leaves the other two alone. The row does not
    move at all with a quantum 16 million times coarser, so the level the precision cases
    record is not quantization, and the accumulator could be far coarser than it is before
    the deposit noticed.

    The other end is checked too. One bit past the headroom floor the run is refused, so
    the guard that keeps the quantum out of the answer is demonstrated rather than assumed.
    """
    print("== 11. what the deposit's quantum contributes ==")
    base = None
    for shift in (0, QUANT_SHIFT):
        root = f"devq{abs(shift)}"
        extra = DEV + '  global%fp32_check = "lockstep"\n  beamlet_size = 8\n'
        if shift:
            extra += f"  global%device_dep_mutate = {shift}\n"
        r = run(exe, wd, f"{root}.in", BASE.format(root=root, extra=extra))
        if r.returncode != 0:
            ok(f"quantum shifted {shift} bits", "run failed", "exit 0", False)
            return
        v = summary(wd, root)["source"]
        if shift == 0:
            base = v
            continue
        ok(f"the source row with the quantum {abs(shift)} bits coarser",
           f"{v:.4e} against {base:.4e}, {v/base:.4f}x",
           f"unmoved to {QUANT_UNMOVED:.3f}x, so the row is not quantization",
           v <= base * QUANT_UNMOVED)

    root = "devqref"
    extra = (DEV + '  global%fp32_check = "lockstep"\n  beamlet_size = 8\n'
             f"  global%device_dep_mutate = {QUANT_REFUSED}\n")
    r = run(exe, wd, f"{root}.in", BASE.format(root=root, extra=extra), expect_fail=True)
    refused = r.returncode != 0 and "CANNOT CARRY ITS BOUND AND ITS PRECISION" in r.stdout
    ok(f"refused: a quantum {abs(QUANT_REFUSED)} bits coarser, one past the floor",
       refused, "True", refused)


def stats_arrays(path, skip=REPRO_SKIP):
    """Every dataset of a statistics file, by path, as raw bytes."""
    out = {}
    with h5py.File(path) as h:
        def walk(name, obj):
            if isinstance(obj, h5py.Dataset) and not name.startswith(skip):
                out[name] = np.asarray(obj[()]).tobytes()
        h.visititems(walk)
    return out


def reproducible(args, wd, exe):
    """
    Two runs of one deck, bit for bit.

    The deposit accumulates in fixed point precisely so that this holds: integer
    addition is associative and commutative where float addition is neither, so the
    order threads reach a cell cannot reach the answer. Before that it could, and this
    deck is the one that showed it: the time-dependent window with shot noise, whose
    eight slices of 1024 macroparticles contend for cells on every step. Run against the
    float deposit, 54 of these 127 arrays differ with the filter on and 54 with it off,
    the field power by 2.0e-07 and 1.5e-07 and the exit power by 1.6e-08 and 2.2e-08. The
    level below is exact equality, asserted rather than toleranced.

    A small dark deck is no test: two runs of the float deposit could already agree on
    one, which is why the no-op floor beneath section 8c existed. This deck carries
    charge, noise and slippage.

    Both filter settings run. They read the source at different points, since with the
    filter on the first of its four source passes converts the accumulator and with it
    off nothing transforms the source and the solve's last pass converts as it reads.
    """
    print("== 9. two runs of one deck, bit for bit ==")
    for on, what in (("T", "the filter on"), ("F", "the filter off")):
        roots = []
        for i in (1, 2):
            root = f"devrep{on}{i}"
            extra = TD_EXTRA + DEV + f"  source_filter = {on}\n"
            r = run(exe, wd, f"{root}.in", BASE.format(root=root, extra=extra), threads="8")
            if r.returncode != 0:
                print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
                ok(f"two runs bit for bit, {what}", "run failed", "exit 0", False)
                return
            roots.append(root)
        a = stats_arrays(wd / f"{roots[0]}.stats.h5")
        b = stats_arrays(wd / f"{roots[1]}.stats.h5")
        same_keys = sorted(a) == sorted(b)
        differ = [k for k in sorted(a) if same_keys and a[k] != b[k]]
        ok(f"two runs bit for bit, {what}: {len(a)} arrays compared",
           "identical" if same_keys and not differ else f"{len(differ)} differ: {differ[:3]}",
           "every array identical", same_keys and not differ)


def migration(args, wd, exe):
    """Slice migration with the device, on check_migration's decks at grid 64."""
    print("== slice migration with the device ==")
    base = cm.BASE.replace("&end\n", "{extra}&end\n")
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

    # 8c. The no-op, exact. The frozen deck is dark, so its power is FP32 noise and says
    # nothing. The read-back beam is what carries the check. With zero moves the migrate = T run differs from
    # migrate = F only in where the readback happened, so the final beam dumps must agree
    # exactly. They once could not: the float deposit made two migrate = F runs of this
    # deck differ between themselves, and the level was three times that flutter with a
    # floor under it for the decks where the flutter came out zero. The fixed-point
    # deposit removed the flutter, so both comparisons below are equality.
    run(args.exe, wd, "devm_f1.in", base.format(root="devmf1", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_f2.in", base.format(root="devmf2", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_t.in", base.format(root="devmt", mig="T", extra=DEV, **frozen), threads="8")
    moved = cm.read_migration_file(wd, "devmt")[0]
    z1 = beam_z(wd / "devmf1-final.beam.h5")
    z2 = beam_z(wd / "devmf2-final.beam.h5")
    zt = beam_z(wd / "devmt-final.beam.h5")
    ok("migration no-op on the device: moves", moved, "0", moved == 0)
    ok("two migrate = F device runs, final beam z", "identical" if np.array_equal(z1, z2)
       else f"{float(np.max(np.abs(z1 - z2))) / 3e-10:.2e} spacings",
       "identical", np.array_equal(z1, z2))
    ok("migration no-op on the device: migrate T vs F final beam z",
       "identical" if np.array_equal(zt, z1)
       else f"{float(np.max(np.abs(zt - z1))) / 3e-10:.2e} spacings",
       "identical", np.array_equal(zt, z1))


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
