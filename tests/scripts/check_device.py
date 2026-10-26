#!/usr/bin/env python3
"""Device backend checks (validation.md, val-device and val-device-unaveraged). The
divergence of the Metal backend from the FP64 path is a recorded number, judged by the
lockstep instrument with the device in the twin's role. The instrument leaves the FP64
physics alone, the check can fail, and everything the backend does not cover is refused.

1. Lockstep levels. A steady one-segment run with device = "metal" and fp32_check =
   "lockstep" must land its worst divergence per quantity inside the recorded device
   ceilings. These are three times the CPU twin's recorded levels where those exist,
   and three times the device's own first measurement for the transverse rows, since
   the CPU twin refreshed those from FP64 and the device runs its own FP32 maps. The
   stream must carry the exact-wrap verdict. Bucket shifts on the fixed-point phase are
   modular arithmetic and are asserted exactly on the device, never to a tolerance.
2. Read-only. The instrumented run must reproduce the uninstrumented run's diag file
   byte for byte. The device twin observes the FP64 path and never steers it.
3. Falsifiable. fp32_mutate = T perturbs a kernel constant (the detuning resonance,
   inside the shader) and coarsens the uploaded phase. The theta row must then grow by
   10 times and the source row by 5 times, so a wrong kernel fails the check.
4. End to end. freerun with the device carries the resident state across steps. The
   phasor must compound past lockstep, and the exit-power divergence must lie in
   (0, 1.6e-3]. The resident production run (fp32_check off) must land in the same
   band against the CPU run's diag rows, so the instrumented twin and the production
   path corroborate each other.
5. Time dependence. An 8-slice shot-noise window with slippage rotating the resident
   record holds the same ceilings, read-only proof and production band.
6. Refused. Wakes, spontaneous radiation, harmonics together with two live
   polarizations, source_filter_mutate, an unsupported grid (the message names the
   nearest supported size), an unknown backend name and keep_escaped_field set in the
   input file each stop the run. The bank is on by default, so a device run that does
   not set it runs without one, says so, and writes neither field file, while its stats
   file still records the escaped and pulse energies.
7. The field set. Harmonic members ride the device (planar segment, harmonics 1 and 3)
   and so do two polarization planes (the crossed undulator of check_two_polarization).
   Each is judged the same four ways: the lockstep rows inside the recorded ceilings
   with the per-member footer lines the set adds, the read-only proof, the mutation
   moving the levels, and the production run against the CPU on the fundamental and on
   the harmonic or the second plane. The harmonic's own identity check is the one-step
   dark deposit's P3/P1 against the Bessel closed form from the dumped particles
   (check_harmonics' identity, its helper shared), which the device holds in FP32. A
   time-dependent Genesis run with shot noise and a third-harmonic field is imported
   on both the CPU and the device at grid 64, the same particles and the same field,
   so both are judged against one reference. One figure is printed and not asserted.
   A quiet-start dark harmonic whose true bunching sits at ~1e-8 of the charge (the
   planar steady-state tier input file) radiates the device's FP32 phase floor instead,
   ~1e-7 of the charge in the phasor row, and its P3 lands two decades above the
   CPU's. The device resolves harmonic bunching down to that floor and not below,
   which the phasor row states per member. The harmonic checks here sit above it.
8. Slice migration, from check_migration's own input files at the device's grid. Migration
   needs no kernel. At an element's last step the walk reads the beam back and releases
   residency before re-slicing, so conservation and phase continuity hold on the device
   exactly as on the CPU, and the check shows the order is right. The production run
   lands in a band against the CPU on the heavy-migration window. The no-op check is
   exact: zero moves, two migrate = F device runs with identical final beams, and the
   migrate = T run identical to them. Lockstep on the migrating window holds the
   ceilings and the read-only proof, the twin re-staging each step across the changed
   fills. A slice outgrowing the setup rectangle exercises the seam's growth-only
   resize, which the log records.
9. Two device runs of one input file agree bit for bit, every array of the statistics file,
   with the source filter on and off.
10. The deposit's accuracy. The lockstep source row of the fixed-point deposit must stay
   under the source ceiling and no worse than the float deposit it replaced, and its
   quantum must sit at least 24 bits below the FP32 spacing of its bound.
11. The deposit's quantum. A quantum 24 bits coarser must leave the source row unmoved,
   and one 30 bits coarser must be refused.
12. The source filter. The device and the CPU must agree on the window power with the
   filter on, and the filter must remove wide-angle power on the device.
13. The unaveraged mode, whose quiver-resolving push, kick and deposit and ledger
   reduction are the kernels this mode adds and whose field solve is the one above.
   It is judged the same ways, with the transform pair's own energy loss measured
   beside them because that loss sets its ceilings. The mode runs the pair nsub times a
   record step where the averaged mode runs it once. The frame series rides the comb
   there too, so a frame is written from arrays the readback refreshed. What a frame
   must carry is the chart: the undulator quiver sits in px as a common offset, which
   the mean sees and an rms does not. See unaveraged() below.

The deposit accumulates in fixed point, so two device runs of one input file agree bit for
bit, and checks 8, 9 and the unaveraged section assert that. Against the CPU the
device is held to the ceilings, and the read-only proof compares FP64 outputs only.

Usage: check_device.py --exe <lucifer> --workdir <dir> --latdir <tests/bmad>
                       --genesis <genesis4>
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

import lucifer.beam as beamio
import check_harmonics as ch
import check_migration as cm
import check_two_polarization as tp
import convert_genesis
import fieldio
from checkrec import record
from lucifer.flat import to_toml as to_groups
from lucifer.stats import read_stats

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

DEV = '  run.device = "metal"\n'

# The recorded ceilings. pz through field are three times the CPU twin's recorded
# levels (validation.md, val-fp32-lockstep), absorbing that the device's sincos is not
# libm and its fused ops round differently. x through py are three times the device's
# own first recorded levels (2.8e-7 worst, steady and time dependent, both builds),
# since the CPU twin's transverse rows measured representation only and the device
# runs its own FP32 maps. A level past its ceiling means the kernels' error changed.

CEIL = {"x": 1e-6, "px": 1e-6, "y": 1e-6, "py": 1e-6,
        "pz": 1.3e-6, "theta": 1.5e-4, "phasor": 2.2e-5, "source": 2.2e-4, "field": 2.2e-4}

# The guard is the median per-step phase increment in ticks of 2 pi / 2^32. The
# recorded minimum is 6.8e6, and the floor sits a factor of several under it. The
# silent failure is this number collapsing toward zero.

GUARD_FLOOR = 1e6

# The end-to-end band. Freerun on this gain segment gives 5.1e-4 on exit power, both
# builds, and the ceiling is three times that. The production run must land in the
# same band, or the two roles disagree and the check fails.

E2E_CEIL = 1.6e-3

# The field set's own levels, each at three times its first measurement on the debug
# build (validation.md, val-device). The harmonic phasor row is the device's harmonic
# bunching floor, 6.6e-8 steady and 2.0e-7 time dependent, three times the
# fundamental's phase error as h = 3 predicts. The harmonic's source and field rows are
# normalized by the fundamental's field and read 3.5e-9 steady and 2.7e-6 time
# dependent, the latter on a shot-noise window where the third harmonic is a quarter
# of the fundamental's power. The harmonic end-to-end band is the third-harmonic power
# against the CPU on a strongly bunched beam, 6.0e-4, which is the fundamental's phase
# error tripled and squared into a power. The Bessel identity is the one-step P3/P1 in
# FP32, 4.4e-6 against 2.7e-14 on the CPU. The isolation level is |ln(Px/P_drift)|
# through the y set on the device, 4.3e-6. The CPU's is 1.6e-15, so 4.3e-6 is what
# FP32 costs the statement that the y set only diffracts Ex. The second plane's own
# band is wider than the fundamental's because Py is a hundredth of Px on this line
# and its relative error carries that ratio, 4.3e-3. Against the one Genesis reference
# the device lands at 1.5e-5 on P1 and 2.3e-5 on P3 where the CPU lands at 8.4e-7 and
# 1.1e-6, both inside these bands.

CEIL_H = {"phasor_h3": 6e-7, "source_h3": 8e-6, "field_h3": 8e-6}
H_E2E_CEIL = 1.8e-3
BESSEL_CEIL = 1.3e-5
ISO_CEIL = 1.3e-5
PY_E2E_CEIL = 1.3e-2
PYPX_FLOOR = 5e-3           # check_two_polarization's afterburner floor.

# Migration: the heavy-migration window's exit power against the CPU, at three times
# the first measurement. The fixed-point deposit makes the no-op comparison exact, so
# the no-op needs no level.
MIG_E2E_CEIL = 3e-3

# The source filter's device-against-CPU ceiling. The device filters the source by four
# FP32 passes of its own where the CPU multiplies inside the field's FP64 transform pair,
# so the two differ by more than the plain solve does. The difference is 1.7e-4 end to
# end.
SF_E2E_CEIL = 1.0e-3

# The unaveraged mode's own ceilings, each three times its first measurement (both
# builds, the check's own input files). They are wider than the averaged ones, and the reason
# is arithmetic rather than the port. This mode runs the transform pair nsub times a
# record step where the averaged mode runs it once, and an FP32 pair loses about 1.3e-7
# of the field's energy every time it runs (validation.md, val-device-unaveraged). The
# transverse rows are the FP32 push's own, so they sit an order above the averaged
# mode's, whose transverse maps are one map a step.
#
# The source column carries this mode's ledger term rather than a source grid, scaled
# by the energy the record step turned over. On a dark-start shot-noise window a slice
# can turn over almost nothing, which widens that row and nothing else.

CEIL_U = {"x": 3.5e-6, "px": 2.5e-5, "y": 3.0e-6, "py": 2.1e-5, "pz": 6.5e-6,
          "theta": 1.1e-5, "phasor": 2.5e-7, "source": 5.0e-3, "field": 2.0e-5}
GUARD_FLOOR_U = 1e9
# The crossed probe's own guard floor. The guard is the median phase increment per
# record step in ticks, so it belongs to the input file and not to the mode. This line's
# elements are 0.6 m where the steady input file's are 4 m, and the guard reads 1.1e8.
GUARD_FLOOR_2P = 1e7
# The two-plane rows all sit inside the one-plane ceilings above. On the crossed
# unaveraged input file they read x 1.9e-6, px 2.5e-6, y 2.2e-6, py 3.0e-6, pz 5.0e-7, theta
# 2.0e-6, phasor 1.6e-8, source 2.2e-4, field 2.0e-6. The field row is the whole record,
# both planes summed, and the source row is the ledger's kick-side term.
U_E2E_CEIL = 5.0e-3

# The frame series with the device resident in this mode. The mean px a frame carries is
# the quiver chart's own. The device's differs from the CPU's by the FP32 round trip of
# the chart, 1.9e-7, rather than by the transform loss that sets the band above. This
# ceiling sits well under that band, so a chart error cannot hide in it.

U_FRAME_CEIL = 1.0e-5

# The end-to-end band against the CPU's own unaveraged run, U_E2E_CEIL above. The
# transform pair's loss below sets it, rather than anything the port chose. On the
# steady input file, whose seed barely grows over one segment, the exit power sits 1.7e-3
# under the CPU's, which is 5340 pairs of that loss carried into the power. The
# time-dependent window reads 1.5e-6, because its power comes from a beam that is
# radiating rather than from a seed that each pair reduces.

# What one FP32 transform pair does to the field's energy, with no beam to feed it.
# The loss is -1.29e-7 a pair on the averaged path and -1.34e-7 on the unaveraged one,
# and it sets the ceilings above. The band is wide enough to carry a grid change and
# narrow enough that a pair which stopped being nearly unitary would show.

PAIR_LOSS = (5e-8, 4e-7)
# The footer's line on what slippage transmitted out of the window, per member.
ESCAPED_RE = re.compile(r"Escaped\s+h1\s+([-+.0-9Ee]+) J,\s+h3\s+([-+.0-9Ee]+) J")
# The header's projection of that loss against what the same run measures. The rate it
# projects from is one machine's fit to four grids, so a factor of two either way.
PAIR_RE = re.compile(r"FP32 transform pair (\d+) times a field record \((\d+) unaveraged "
                     r"substeps\), projected field-energy loss ([-+.0-9Ee]+) a pair")
PROJ_BAND = (0.5, 2.0)

# The planar segment of check_harmonics, since fc(3) is alive there where the Aramis
# segment is helical and couples only the fundamental, at the device's grid.
PLANAR_LAT = ch.BMAD_LAT

HARM_BASE = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 8192
  beam.bunch_charge = 1.000692285594e-15
  beam.sig_z = 0
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
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
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    exe = pathlib.Path(args.exe).resolve()

    # 1. Steady lockstep: the device in the twin's role, the recorded ceilings, the
    # exact-wrap verdict from the device's own arithmetic.
    run(args.exe, wd, "dev_ss.in", BASE.format(root="devss",
        extra=DEV + '  run.fp32_check = "lockstep"\n'))
    s = summary(wd, "devss")
    for q in CEIL:
        ok(f"device lockstep steady worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q],
           num=s[q], bound=CEIL[q], kind="max")
    ok("device lockstep steady guard [ticks]", f"{s['guard']:.3e}", f">= {GUARD_FLOOR:.0e}",
       s["guard"] >= GUARD_FLOOR,
       num=s["guard"], bound=GUARD_FLOOR, kind="min")
    ok("bucket wraps asserted exact on device", s.get("wrap", "absent"), "T", s.get("wrap") == "T")

    # 2. Read-only: the instrumented run's FP64 output is the uninstrumented run's.
    run(args.exe, wd, "dev_off.in", BASE.format(root="devoff", extra=""))
    same = (wd / "devss.diag.txt").read_bytes() == (wd / "devoff.diag.txt").read_bytes()
    ok("FP64 diag byte-identical, device twin on vs off", same, "True", same)

    # 3. Falsifiable: the perturbed kernel constant and the coarsened upload must
    # move the recorded levels.
    run(args.exe, wd, "dev_mut.in", BASE.format(root="devmut",
        extra=DEV + '  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "devmut")
    ok("mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"],
       num=m["theta"], bound=10 * s["theta"], kind="min")
    ok("mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"],
       num=m["source"], bound=5 * s["source"], kind="min")

    # 4. End to end, both roles. The freerun twin carries the resident state and
    # compounds. The production run replaces the CPU step outright. Both must land in
    # the recorded band, so they corroborate each other.
    run(args.exe, wd, "dev_fr.in", BASE.format(root="devfr",
        extra=DEV + '  run.fp32_check = "freerun"\n'))
    f = summary(wd, "devfr")
    ok("device freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"],
       num=f["phasor"], bound=s["phasor"], kind="min")
    ok("device freerun end-to-end power", f"{f['e2e_power']:.3e}",
       f"in (0, {E2E_CEIL:.1e}]", 0 < f["e2e_power"] <= E2E_CEIL,
       num=f["e2e_power"], bound=[0, E2E_CEIL], kind="range")

    run(args.exe, wd, "dev_pr.in", BASE.format(root="devpr", extra=DEV))
    pc = diag_power(wd, "devoff", 1)
    pd = diag_power(wd, "devpr", 1)
    rel = abs(pd - pc) / pc
    ok("device production exit power vs CPU", f"{rel:.3e}",
       f"in (0, {E2E_CEIL:.1e}]", 0 < rel <= E2E_CEIL,
       num=rel, bound=[0, E2E_CEIL], kind="range")

    # 5. Time dependence: slippage rotates the resident record. Same ceilings,
    # same read-only proof, same production band, at 8 threads.
    run(args.exe, wd, "dev_td.in", BASE.format(root="devtd",
        extra=TD_EXTRA + DEV + '  run.fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devtd")
    for q in CEIL:
        ok(f"device lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q],
           num=t[q], bound=CEIL[q], kind="max")
    ok("device lockstep TD guard [ticks]", f"{t['guard']:.3e}", f">= {GUARD_FLOOR:.0e}",
       t["guard"] >= GUARD_FLOOR,
       num=t["guard"], bound=GUARD_FLOOR, kind="min")
    run(args.exe, wd, "dev_tdoff.in", BASE.format(root="devtdoff", extra=TD_EXTRA), threads="8")
    same = (wd / "devtd.diag.txt").read_bytes() == (wd / "devtdoff.diag.txt").read_bytes()
    ok("TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "dev_tdpr.in", BASE.format(root="devtdpr", extra=TD_EXTRA + DEV), threads="8")
    pc = diag_power(wd, "devtdoff", 8)
    pd = diag_power(wd, "devtdpr", 8)
    rel = abs(pd - pc) / pc
    ok("device production TD window power vs CPU", f"{rel:.3e}",
       f"<= {E2E_CEIL:.1e}", rel <= E2E_CEIL,
       num=rel, bound=E2E_CEIL, kind="max")

    # 6. Refused, never a quiet CPU fallback.
    refusals = [
        ("dev_rw.in", DEV + '  wake_on = T\n  wake_radius = 2.5e-3\n'
                            '  wake_conductivity = 5.813e7\n  wake_relaxation = 8.1e-6\n',
         "DOES NOT COVER WAKES"),
        ("dev_rr.in", DEV + '  bmad.radiation_damping_on = T\n',
         "DOES NOT COVER SPONTANEOUS RADIATION"),
        ("dev_rn.in", '  run.device = "cuda"\n', "UNRECOGNIZED DEVICE"),
        ("dev_rb.in", DEV + '  keep_escaped_field = T\n', "DOES NOT COVER THE ESCAPED-FIELD BANK"),
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

    r = run(args.exe, wd, "dev_bank.in", BASE.format(root="dev_bank", extra=DEV))
    with h5py.File(wd / "dev_bank.stats.h5") as h5:
        recorded = "field/total/escaped_energy" in h5 and "field/total/pulse_energy" in h5
    dropped = ("SO THIS RUN KEEPS NONE" in r.stdout and not list(wd.glob("dev_bank-*.fld.h5"))
               and recorded)
    ok("bank at its default: the device runs without it, says so, writes no field file",
       dropped, "True", dropped)

    field_set(args, wd, exe)
    source_filter(args, wd, exe)
    unaveraged(args, wd, exe)

    print("PASS" if not FAILED else "FAIL")
    return 1 if FAILED else 0



def unaveraged(args, wd, exe):
    """
    The unaveraged mode on the device (validation.md, val-device-unaveraged).

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
    for f in ("aramis_1seg.bmad", "aramis_1seg_unavg.bmad", "crossed_probe.bmad",
              "spont_probe.bmad"):
        (wd / f).write_bytes((latdir / f).read_bytes())
    (wd / "dvu_cross.bmad").write_text(
        "call, file = crossed_probe.bmad\nuse, CROSSED\nwiggler::*[FEL_METHOD] = unaveraged\n")
    lat = "aramis_1seg_unavg.bmad"
    base = BASE.replace("aramis_1seg.bmad", lat)

    # 1. Steady lockstep, the device in the twin's role.
    run(args.exe, wd, "dvu_ss.in", base.format(root="dvuss",
        extra=DEV + '  run.fp32_check = "lockstep"\n'))
    s = summary(wd, "dvuss")
    for q in CEIL_U:
        ok(f"unaveraged lockstep steady worst_{q}", f"{s[q]:.3e}",
           f"<= {CEIL_U[q]:.1e}", s[q] <= CEIL_U[q],
           num=s[q], bound=CEIL_U[q], kind="max")
    ok("unaveraged lockstep steady guard [ticks]", f"{s['guard']:.3e}",
       f">= {GUARD_FLOOR_U:.0e}", s["guard"] >= GUARD_FLOOR_U,
       num=s["guard"], bound=GUARD_FLOOR_U, kind="min")

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
        extra=DEV + '  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "dvumut")
    for q, factor in (("theta", 10), ("phasor", 50), ("field", 20)):
        ok(f"unaveraged mutation moves worst_{q}", f"{m[q]:.3e} vs {s[q]:.3e}",
           f">= {factor}x", m[q] >= factor * s[q],
           num=m[q], bound=factor * s[q], kind="min")

    # 4. The production role, against the CPU's own unaveraged run.
    run(args.exe, wd, "dvu_pr.in", base.format(root="dvupr", extra=DEV))
    pc = diag_power(wd, "dvuoff", 1)
    pd = diag_power(wd, "dvupr", 1)
    rel = abs(pd - pc) / pc
    ok("unaveraged production power vs CPU", f"{rel:.3e}", f"<= {U_E2E_CEIL:.1e}",
       rel <= U_E2E_CEIL,
       num=rel, bound=U_E2E_CEIL, kind="max")

    # 4b. The frame series with the device resident in this mode. The comb's positions
    # are where the resident state comes back to the host, so a frame is written from
    # arrays the device refreshed rather than from stale ones. The chart is the thing to
    # check: a frame taken mid-segment carries the undulator quiver in px, and a
    # statistic that cannot see a common offset would pass with the quiver missing, so
    # the mean is compared and not the spread.
    frame_extra = '  run.comb_ds_save = 0.5\n  run.dump_at_comb = T\n'
    run(args.exe, wd, "dvu_fr.in", base.format(root="dvufr", extra=DEV + frame_extra))
    run(args.exe, wd, "dvu_froff.in", base.format(root="dvufroff", extra=frame_extra))
    fd = sorted(wd.glob("dvufr-[0-9]*.beam.h5"))
    fc = sorted(wd.glob("dvufroff-[0-9]*.beam.h5"))
    ok("unaveraged frames on the device, count against the CPU",
       f"{len(fd)} against {len(fc)}", "equal and nonzero",
       len(fd) == len(fc) and len(fd) > 0,
       num=len(fd), bound=len(fc), kind="equal")
    ok("unaveraged field frames on the device, count against the CPU",
       f"{len(sorted(wd.glob('dvufr-[0-9]*.wf.h5')))} against {len(fd)}", "equal",
       len(sorted(wd.glob("dvufr-[0-9]*.wf.h5"))) == len(fd),
       num=len(sorted(wd.glob("dvufr-[0-9]*.wf.h5"))), bound=len(fd), kind="equal")

    # The frames steer nothing. Inside the unaveraged segment the frame carries the
    # plane's field where the record holds the substep's midpoint field, so the writer
    # diffracts, and it does so on a copy: the same run with the same comb and no frames
    # ends in the same beam and field to the bit, on the device and on the CPU.
    comb_only = '  run.comb_ds_save = 0.5\n'
    run(args.exe, wd, "dvu_frnf.in", base.format(root="dvufrnf", extra=DEV + comb_only))
    run(args.exe, wd, "dvu_froffnf.in", base.format(root="dvufroffnf", extra=comb_only))
    for tag, a, b in (("device", "dvufr", "dvufrnf"), ("CPU", "dvufroff", "dvufroffnf")):
        same = all(dumps_identical(wd / f"{a}-final.{k}.h5", wd / f"{b}-final.{k}.h5")
                   for k in ("beam", "wf"))
        ok(f"unaveraged frames on against off, final beam and field identical ({tag})",
           same, "True", same)

    if fd and len(fd) == len(fc):
        def read(path):
            with h5py.File(path) as h:
                it = list(h["data"])[0]
                m = h["lucifer/frame"].attrs.get("felMethod")
                m = np.ravel(m)[0] if m is not None else None
                if isinstance(m, bytes):
                    m = m.decode()
                return (float(np.ravel(h["lucifer/frame"].attrs["sPosition"])[0]), m,
                        float(np.mean(h[f"data/{it}/particles/electron/momentum/x"])))

        d = [read(f) for f in fd]
        c = [read(f) for f in fc]
        sdev = max(abs(a[0] - b[0]) for a, b in zip(d, c))
        ok("unaveraged frames on the device land at the CPU's own s", f"{sdev:.3e}",
           "<= 1e-9", sdev <= 1e-9,
           num=sdev, bound=1e-9, kind="max")

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
           "> 1e4", quiver > 1e4,
           num=quiver, bound=1e4, kind="min")
        scale = max(max(abs(b[2]) for b in c), 1e-300)
        rel = max(abs(a[2] - b[2]) for a, b in zip(d, c)) / scale
        ok("unaveraged frame px mean, device against CPU", f"{rel:.3e}",
           f"<= {U_FRAME_CEIL:.1e}", rel <= U_FRAME_CEIL,
           num=rel, bound=U_FRAME_CEIL, kind="max")

    # 5. A time-dependent window: slippage rotates the resident record between record
    # steps and the ledger's banked terms cross the seam every record step.
    run(args.exe, wd, "dvu_td.in", base.format(root="dvutd",
        extra=DEV + TD_EXTRA + '  run.fp32_check = "lockstep"\n'))
    t = summary(wd, "dvutd")
    for q in CEIL_U:
        ok(f"unaveraged lockstep TD worst_{q}", f"{t[q]:.3e}",
           f"<= {CEIL_U[q]:.1e}", t[q] <= CEIL_U[q],
           num=t[q], bound=CEIL_U[q], kind="max")
    run(args.exe, wd, "dvu_tdoff.in", base.format(root="dvutdoff", extra=TD_EXTRA))
    same = (wd / "dvutd.diag.txt").read_bytes() == (wd / "dvutdoff.diag.txt").read_bytes()
    ok("unaveraged TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "dvu_tdpr.in", base.format(root="dvutdpr", extra=DEV + TD_EXTRA))
    pc = diag_power(wd, "dvutdoff", 8)
    pd = diag_power(wd, "dvutdpr", 8)
    rel = abs(pd - pc) / pc
    ok("unaveraged production TD window power vs CPU", f"{rel:.3e}",
       f"<= {U_E2E_CEIL:.1e}", rel <= U_E2E_CEIL,
       num=rel, bound=U_E2E_CEIL, kind="max")

    # 6. Two device runs of one unaveraged input file agree bit for bit. The deposit
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
    # count must be the input file's and the projection must sit near what the run measured.
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
        ok(f"FP32 transform pair energy loss, {mode} input file", f"{loss:.3e} a pair",
           f"in [{PAIR_LOSS[0]:.0e}, {PAIR_LOSS[1]:.0e}]", good,
           num=loss, bound=list(PAIR_LOSS), kind="range")
        m = PAIR_RE.search(r.stdout)
        counts = (int(m.group(1)), int(m.group(2))) if m else None
        ok(f"header counts the transform pairs, {mode} input file", counts, f"== {(npair, nsub)}",
           counts == (npair, nsub),
           num=counts, bound=[npair, nsub], kind="equal")
        ratio = float(m.group(3)) / loss if m else float("nan")
        ok(f"header projection over measured loss, {mode} input file", f"{ratio:.3f}",
           f"in [{PROJ_BAND[0]}, {PROJ_BAND[1]}]", PROJ_BAND[0] <= ratio <= PROJ_BAND[1],
           num=ratio, bound=list(PROJ_BAND), kind="range")

    # 8. Two live polarizations. The mode needs no polarization vector: the two real
    # kinetic momenta work against and deposit into their own field components, so the
    # device's kick reads both planes and fills both source planes, and each plane meets
    # the diffraction with its own source. The crossed line is the input file for it, an x set
    # then a y set, so both planes carry power by the end.
    xb = base.replace('lattice = "' + lat + '"', 'lattice = "dvu_cross.bmad"')

    run(args.exe, wd, "dvu_2p.in", xb.format(root="dvu2p",
        extra=DEV + '  run.fp32_check = "lockstep"\n'))
    t = summary(wd, "dvu2p")
    for q in CEIL_U:
        ok(f"two-plane unaveraged lockstep worst_{q}", f"{t[q]:.3e}", f"<= {CEIL_U[q]:.1e}",
           t[q] <= CEIL_U[q],
           num=t[q], bound=CEIL_U[q], kind="max")
    ok("two-plane unaveraged lockstep guard [ticks]", f"{t['guard']:.3e}",
       f">= {GUARD_FLOOR_2P:.0e}", t["guard"] >= GUARD_FLOOR_2P,
       num=t["guard"], bound=GUARD_FLOOR_2P, kind="min")

    run(args.exe, wd, "dvu_2poff.in", xb.format(root="dvu2poff", extra=""))
    same = (wd / "dvu2p.diag.txt").read_bytes() == (wd / "dvu2poff.diag.txt").read_bytes()
    ok("two-plane unaveraged FP64 diag byte-identical, device twin on vs off", same, "True", same)
    same = (wd / "dvu2p.ledger.txt").read_bytes() == (wd / "dvu2poff.ledger.txt").read_bytes()
    ok("two-plane unaveraged FP64 ledger byte-identical, device twin on vs off", same, "True", same)

    run(args.exe, wd, "dvu_2pm.in", xb.format(root="dvu2pm",
        extra=DEV + '  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "dvu2pm")
    for q, factor in (("theta", 10), ("phasor", 50), ("field", 20)):
        ok(f"two-plane unaveraged mutation moves worst_{q}", f"{m[q]:.3e} vs {t[q]:.3e}",
           f">= {factor}x", m[q] >= factor * t[q],
           num=m[q], bound=factor * t[q], kind="min")

    # The production role. Each plane is compared against its own CPU power, since a
    # total would hide a plane that went dark. A propagator that asks for a member the
    # kernel buffer does not hold leaves its plane dark while the first plane and the
    # total stay close to the CPU's. The y plane of this line is the afterburner's
    # transfer and carries 1.3e-4 of the x plane. Its own relative agreement, 3.0e-2, is
    # what single precision does to a weak coherent transfer, where a plane that went
    # dark would read 1.
    run(args.exe, wd, "dvu_2ppr.in", xb.format(root="dvu2ppr", extra=DEV))
    for comp, ceil in (("x", U_E2E_CEIL), ("y", 1.0e-1)):
        pc = tp.dump_power(wd, "dvu2poff-final.wf.h5", comp)
        pd = tp.dump_power(wd, "dvu2ppr-final.wf.h5", comp)
        rel = abs(pd - pc) / max(pc, 1e-300)
        ok(f"two-plane unaveraged production power {comp} vs CPU", f"{rel:.3e}",
           f"<= {ceil:.1e}", pc > 0 and rel <= ceil,
           num=rel, bound=ceil, kind="max")

    # And with both planes: frames inside the segment carry both planes back onto the
    # plane, on a copy, and the run with frames ends where the run without them does.
    run(args.exe, wd, "dvu_2pfr.in", xb.format(root="dvu2pfr",
        extra=DEV + '  run.comb_ds_save = 0.5\n  run.dump_at_comb = T\n'))
    run(args.exe, wd, "dvu_2pfrnf.in", xb.format(root="dvu2pfrnf",
        extra=DEV + '  run.comb_ds_save = 0.5\n'))
    same = all(dumps_identical(wd / f"dvu2pfr-final.{k}.h5", wd / f"dvu2pfrnf-final.{k}.h5")
               for k in ("beam", "wf"))
    ok("two-plane unaveraged frames on against off, final beam and field identical (device)",
       same, "True", same)

    # The transform pair's loss on two planes, with no charge to feed either. The loss is
    # each plane's own, so the pair costs the set the same fraction it costs one plane,
    # and the band is the one section 7 holds.
    r = run(args.exe, wd, "dvu_d2.in",
            xb.format(root="dvud2", extra=DEV + "  run.write_initial = T\n").replace(
                "bunch_charge = 1.000692285594e-15", "bunch_charge = 1e-30"))
    p0 = sum(tp.dump_power(wd, "dvud2-initial.wf.h5", c) for c in ("x", "y"))
    p1 = sum(tp.dump_power(wd, "dvud2-final.wf.h5", c) for c in ("x", "y"))
    npair = int(PAIR_RE.search(r.stdout).group(1))   # this line's own count, not the steady input file's
    loss = (p0 - p1) / p0 / npair
    good = PAIR_LOSS[0] <= loss <= PAIR_LOSS[1]
    ok("FP32 transform pair energy loss, two-plane input file", f"{loss:.3e} a pair",
       f"in [{PAIR_LOSS[0]:.0e}, {PAIR_LOSS[1]:.0e}]", good,
       num=loss, bound=list(PAIR_LOSS), kind="range")

    # 9. The handedness the mode selects, on the device. A helical element radiates one
    # circular polarization and couples to it, and which one is the electron's own. A
    # total power cannot tell the two apart, so the control is seeded:
    # a cold beam driven by a uniform seed as the scalar envelope and as each circular
    # pair of the same total intensity, and the rms energy modulation says which pair the
    # device couples to. It must be the CPU's answer, and the CPU's is check_two_polarization's.
    unavg_handedness(args, wd, exe)


def unavg_handedness(args, wd, exe):
    """The seeded handedness control of check_two_polarization, run on the device."""
    (wd / "dvu_hel.bmad").write_text(
        "call, file = spont_probe.bmad\nwiggler::*[FEL_METHOD] = unaveraged\n")
    cold = tp.ss(tp.NML).replace("seed_waist_size = 30e-6", "seed_waist_size = 1e-2").replace(
        "seed_power = 1e4", "seed_power = 1e9").replace(
        "beam.sig_pz = 8.804506566858e-05", "beam.sig_pz = 1e-12").replace(
        "beam.bunch_charge = 8.0e-15", "beam.bunch_charge = 8.0e-18")
    run(args.exe, wd, "dvh_src.in", cold.format(lat="dvu_hel.bmad", root="dvhsrc",
        extra="  write_initial = T\n  load_only = T\n"))

    def circular(src, dst, sign):
        shutil.copy(wd / src, wd / dst)
        with h5py.File(wd / dst, "r+") as h5:
            m = h5["data"][sorted(h5["data"])[0]]["meshes/electricField"]
            x = m["x"][...] / np.sqrt(2.0)
            m["x"][...] = x
            y = m.create_dataset("y", data=sign * 1j * x)
            for k, v in m["x"].attrs.items():
                y.attrs[k] = v

    circular("dvhsrc-initial.wf.h5", "dvh_plus.wf.h5", +1)
    circular("dvhsrc-initial.wf.h5", "dvh_minus.wf.h5", -1)
    sig = {}
    for tag, ff in (("scalar", "dvhsrc-initial.wf.h5"), ("plus", "dvh_plus.wf.h5"),
                    ("minus", "dvh_minus.wf.h5")):
        run(args.exe, wd, f"dvh_{tag}.in", cold.format(lat="dvu_hel.bmad", root=f"dvh{tag}",
            extra=f'  beam.file = "dvhsrc-initial.beam.h5"\n  field.file = "{ff}"\n'
                  '  load_mode = "keep"\n' + DEV))
        sig[tag] = float(tp.diag(wd, f"dvh{tag}")[-1, 0, 7])
    rel = abs(sig["plus"] - sig["scalar"]) / sig["scalar"]
    ok("handedness on the device: (1, +i) reproduces the scalar seed's energy modulation",
       f"{rel:.3e} [{sig['plus']:.4f} vs {sig['scalar']:.4f} eV]", "<= 1.0e-5", rel <= 1.0e-5,
       num=rel, bound=1.0e-5, kind="max")
    rel = sig["minus"] / sig["scalar"]
    ok("handedness on the device: (1, -i) leaves the beam unmodulated",
       f"{rel:.3e} [{sig['minus']:.3e} eV]", "<= 1.0e-2", rel <= 1.0e-2,
       num=rel, bound=1.0e-2, kind="max")

def source_filter(args, wd, exe):
    """
    The source filter on the device (fel-physics.md sec-source-filter).

    The CPU adds the filtered source inside the field's own transform pair. The device's
    fused solve adds the source in real space in its last pass, so the source is filtered
    in place first by four passes of its own and the solve is left untouched. The two
    reach the same quantity by a different route, and this check compares them. The
    transcription itself is checked against Genesis4 on the CPU
    (check_source_filter.py), which the device cannot join: Genesis4's reference chain is
    a 255-point grid and the Metal solver takes powers of two.

    The mutation hook is the CPU check's own, so the device refuses it rather than
    carrying a second wrong solve.
    """
    # A 256-point grid rather than the 64 the rest of this file uses. The filter removes
    # emission at angles the grid carries beyond the edge, and a coarse grid carries
    # almost none. At 6.3 um cells the Nyquist angle is 7.9 urad against an edge near 3,
    # so there is nothing to remove, and the check could not see a filter that did nothing.
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
       f"<= {SF_E2E_CEIL:.1e}", abs(pd - pc) / pc <= SF_E2E_CEIL,
       num=abs(pd - pc) / pc, bound=SF_E2E_CEIL, kind="max")
    ok("source filter on the device removes wide-angle power", f"{pn / pd:.2f}x", "> 1.5",
       pn / pd > 1.5,
       num=pn / pd, bound=1.5, kind="min")

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
    hb = HARM_BASE.replace('lattice = "{lat}"', 'lattice = "dev_planar_mid.bmad"')
    run(args.exe, wd, "devh_ss.in", hb.format(root="devhss", extra=DEV + '  run.fp32_check = "lockstep"\n'))
    s = summary(wd, "devhss")
    for q in CEIL:
        ok(f"harmonic lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q],
           num=s[q], bound=CEIL[q], kind="max")
    for q in CEIL_H:
        ok(f"harmonic lockstep worst_{q} (per-member footer)", f"{s.get(q, float('nan')):.3e}",
           f"<= {CEIL_H[q]:.1e}", q in s and s[q] <= CEIL_H[q],
           num=s.get(q), bound=CEIL_H[q], kind="max")
    ok("harmonic footer carries the fundamental too", "phasor_h1" in s, "True", "phasor_h1" in s)
    run(args.exe, wd, "devh_off.in", hb.format(root="devhoff", extra=""))
    same = (wd / "devhss.diag.txt").read_bytes() == (wd / "devhoff.diag.txt").read_bytes()
    ok("harmonic FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devh_mut.in", hb.format(root="devhmut",
        extra=DEV + '  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "devhmut")
    ok("harmonic mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"],
       num=m["theta"], bound=10 * s["theta"], kind="min")
    ok("harmonic mutation moves worst_source_h3", f"{m['source_h3']:.3e} vs {s['source_h3']:.3e}",
       ">= 5x", m["source_h3"] >= 5 * s["source_h3"],
       num=m["source_h3"], bound=5 * s["source_h3"], kind="min")
    run(args.exe, wd, "devh_fr.in", hb.format(root="devhfr", extra=DEV + '  run.fp32_check = "freerun"\n'))
    f = summary(wd, "devhfr")
    ok("harmonic freerun compounds past lockstep (phasor)", f"{f['phasor']:.3e} vs {s['phasor']:.3e}",
       "> 1x", f["phasor"] > s["phasor"],
       num=f["phasor"], bound=s["phasor"], kind="min")
    run(args.exe, wd, "devh_pr.in", hb.format(root="devhpr", extra=DEV))
    c1, c3 = harm_powers(wd, "devhoff")
    d1, d3 = harm_powers(wd, "devhpr")
    ok("harmonic production P1 vs CPU (bunched beam)", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL,
       num=abs(d1 - c1) / c1, bound=E2E_CEIL, kind="max")
    ok("harmonic production P3 vs CPU (bunched beam)", f"{abs(d3-c3)/c3:.3e}", f"<= {H_E2E_CEIL:.1e}",
       abs(d3 - c3) / c3 <= H_E2E_CEIL,
       num=abs(d3 - c3) / c3, bound=H_E2E_CEIL, kind="max")

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
       dv <= BESSEL_CEIL,
       num=dv, bound=BESSEL_CEIL, kind="max")

    # 7c. Time dependence with the set: the harmonics check's time-dependent input file (8
    # slices, shot noise, dark third harmonic) at the device's grid, slippage rotating
    # every plane. Shot noise puts b3 at ~1/sqrt(N), far above the floor, so P3 is a real
    # comparison.
    td = ch.NML_TD.replace("grid_n_pts = 63", "grid_n_pts = 64").replace('transport_model = "genesis"',
                                                                            'transport_model = "bmad"')
    run(args.exe, wd, "devh_td.in", td.format(lat="dev_planar.bmad", root="devhtd",
        extra=DEV + '  run.fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devhtd")
    for q in CEIL:
        ok(f"harmonic lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q],
           num=t[q], bound=CEIL[q], kind="max")
    for q in CEIL_H:
        ok(f"harmonic lockstep TD worst_{q}", f"{t.get(q, float('nan')):.3e}", f"<= {CEIL_H[q]:.1e}",
           q in t and t[q] <= CEIL_H[q],
           num=t.get(q), bound=CEIL_H[q], kind="max")
    r_off = run(args.exe, wd, "devh_tdoff.in", td.format(lat="dev_planar.bmad", root="devhtdoff", extra=""), threads="8")
    same = (wd / "devhtd.diag.txt").read_bytes() == (wd / "devhtdoff.diag.txt").read_bytes()
    ok("harmonic TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    r_pr = run(args.exe, wd, "devh_tdpr.in", td.format(lat="dev_planar.bmad", root="devhtdpr", extra=DEV), threads="8")

    # The energy slippage transmitted out of the window, per member. Bookkeeping that
    # formed the transmitted slice's light time from the harmonic's wavelength alone
    # would give a third harmonic's escaped energy as a third of what left the record,
    # while the surviving field and the window power stayed right (validation.md,
    # val-device).
    e_off, e_pr = ESCAPED_RE.search(r_off.stdout), ESCAPED_RE.search(r_pr.stdout)
    for name, m in (("CPU", e_off), ("device", e_pr)):
        ok(f"harmonic TD footer states the escaped energy ({name})", bool(m), "True", bool(m))
    if e_off and e_pr:
        for grp, label, ceil in ((1, "h1", E2E_CEIL), (2, "h3", H_E2E_CEIL)):
            c, d = float(e_off.group(grp)), float(e_pr.group(grp))
            ok(f"harmonic production TD escaped energy {label} vs CPU",
               f"{abs(d-c)/c:.3e} [{d:.4e} vs {c:.4e} J]", f"<= {ceil:.1e}",
               c > 0 and abs(d - c) / c <= ceil,
               num=abs(d - c) / c, bound=ceil, kind="max")
    c1, c3 = harm_powers(wd, "devhtdoff", window=True)
    d1, d3 = harm_powers(wd, "devhtdpr", window=True)
    ok("harmonic production TD window P1 vs CPU", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL,
       num=abs(d1 - c1) / c1, bound=E2E_CEIL, kind="max")
    ok("harmonic production TD window P3 vs CPU (shot-noise b3)", f"{abs(d3-c3)/c3:.3e}",
       f"<= {H_E2E_CEIL:.1e}", abs(d3 - c3) / c3 <= H_E2E_CEIL,
       num=abs(d3 - c3) / c3, bound=H_E2E_CEIL, kind="max")

    # 7d. One Genesis reference for both: check_harmonics' planar input file made time
    # dependent (8 slices, shot noise on, so the third harmonic's bunching is real on
    # every side) at the device's grid, its starting state imported by the CPU and the
    # device alike. The CPU's level is the anchor the device's is read against.
    genesis_harmonics(args, wd, exe)

    # 7e. The quiet-start floor, printed: the planar steady-state tier input file at grid 64,
    # dark third harmonic, no shot noise. The fundamental holds its band. The harmonic's
    # true bunching is ~1e-8 of the charge, below the device's FP32 phase floor, so its
    # P3 is the floor radiating and is printed, not asserted.
    qs = ch.NML_TD.replace("grid_n_pts = 63", "grid_n_pts = 64").replace('transport_model = "genesis"',
                                                                            'transport_model = "bmad"')
    qs = tp.ss(qs).replace("beam.bunch_charge = 8.0e-15", "beam.bunch_charge = 1.000692285594e-15")
    qs = qs.replace("seed_power = 1e4", "seed_power = 5e3")
    run(args.exe, wd, "devh_qs.in", qs.format(lat="dev_planar.bmad", root="devhqs", extra=""))
    run(args.exe, wd, "devh_qsdev.in", qs.format(lat="dev_planar.bmad", root="devhqsdev", extra=DEV))
    c1, c3 = harm_powers(wd, "devhqs")
    d1, d3 = harm_powers(wd, "devhqsdev")
    ok("quiet-start dark harmonic: fundamental P1 vs CPU", f"{abs(d1-c1)/c1:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(d1 - c1) / c1 <= E2E_CEIL,
       num=abs(d1 - c1) / c1, bound=E2E_CEIL, kind="max")
    print(f"--- quiet-start dark harmonic: device P3 / CPU P3 = {d3/c3:.3e} "
          f"(CPU {c3:.3e} W, device {d3:.3e} W): the FP32 floor radiating, recorded")

    # 7f. Two polarizations: the crossed undulator of check_two_polarization, steady
    # state and time dependent, on its own input file at the device's grid.
    print("== the field set: two polarization planes on the device ==")
    (wd / "crossed_probe.bmad").write_bytes((pathlib.Path(args.latdir) / "crossed_probe.bmad").read_bytes())
    (wd / "devp_c.bmad").write_text("call, file = crossed_probe.bmad\nuse, CROSSED\n")
    (wd / "devp_xd.bmad").write_text("call, file = crossed_probe.bmad\nDEQ: pipe, l = 0.60\n"
                                     "XDRIFT: line = (UNDX, D1, DEQ)\nuse, XDRIFT\n")
    pn = tp.ss(tp.NML)
    run(args.exe, wd, "devp_ss.in", pn.format(lat="devp_c.bmad", root="devpss",
        extra=DEV + '  run.fp32_check = "lockstep"\n'))
    s = summary(wd, "devpss")
    for q in CEIL:
        ok(f"two-plane lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q],
           num=s[q], bound=CEIL[q], kind="max")
    run(args.exe, wd, "devp_off.in", pn.format(lat="devp_c.bmad", root="devpoff", extra=""))
    same = (wd / "devpss.diag.txt").read_bytes() == (wd / "devpoff.diag.txt").read_bytes()
    ok("two-plane FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devp_mut.in", pn.format(lat="devp_c.bmad", root="devpmut",
        extra=DEV + '  run.fp32_check = "lockstep"\n  run.fp32_mutate = T\n'))
    m = summary(wd, "devpmut")
    ok("two-plane mutation moves worst_theta", f"{m['theta']:.3e} vs {s['theta']:.3e}",
       ">= 10x", m["theta"] >= 10 * s["theta"],
       num=m["theta"], bound=10 * s["theta"], kind="min")
    ok("two-plane mutation moves worst_source", f"{m['source']:.3e} vs {s['source']:.3e}",
       ">= 5x", m["source"] >= 5 * s["source"],
       num=m["source"], bound=5 * s["source"], kind="min")
    run(args.exe, wd, "devp_pr.in", pn.format(lat="devp_c.bmad", root="devppr", extra=DEV))
    run(args.exe, wd, "devp_dr.in", pn.format(lat="devp_xd.bmad", root="devpdr", extra=DEV))
    px_c = tp.dump_power(wd, "devpoff-final.wf.h5", "x")
    py_c = tp.dump_power(wd, "devpoff-final.wf.h5", "y")
    px_d = tp.dump_power(wd, "devppr-final.wf.h5", "x")
    py_d = tp.dump_power(wd, "devppr-final.wf.h5", "y")
    p_drift = tp.dump_power(wd, "devpdr-final.wf.h5", "x")
    ok("crossed production Px vs CPU", f"{abs(px_d-px_c)/px_c:.3e}", f"<= {E2E_CEIL:.1e}",
       abs(px_d - px_c) / px_c <= E2E_CEIL,
       num=abs(px_d - px_c) / px_c, bound=E2E_CEIL, kind="max")
    ok("crossed production Py vs CPU", f"{abs(py_d-py_c)/py_c:.3e}", f"<= {PY_E2E_CEIL:.1e}",
       abs(py_d - py_c) / py_c <= PY_E2E_CEIL,
       num=abs(py_d - py_c) / py_c, bound=PY_E2E_CEIL, kind="max")
    iso = abs(np.log(px_d / p_drift))
    ok("crossed device: x-field isolation through the y set, |ln(Px/P_drift)|", f"{iso:.3e}",
       f"<= {ISO_CEIL:.1e}", iso <= ISO_CEIL,
       num=iso, bound=ISO_CEIL, kind="max")
    ok("crossed device: afterburner floor, Py/Px", f"{py_d/px_d:.3e}", f">= {PYPX_FLOOR:.0e}",
       py_d / px_d >= PYPX_FLOOR,
       num=py_d / px_d, bound=PYPX_FLOOR, kind="min")
    tdn = tp.NML
    run(args.exe, wd, "devp_td.in", tdn.format(lat="devp_c.bmad", root="devptd",
        extra=DEV + '  run.fp32_check = "lockstep"\n'), threads="8")
    t = summary(wd, "devptd")
    for q in CEIL:
        ok(f"two-plane lockstep TD worst_{q}", f"{t[q]:.3e}", f"<= {CEIL[q]:.1e}", t[q] <= CEIL[q],
           num=t[q], bound=CEIL[q], kind="max")
    run(args.exe, wd, "devp_tdoff.in", tdn.format(lat="devp_c.bmad", root="devptdoff", extra=""), threads="8")
    same = (wd / "devptd.diag.txt").read_bytes() == (wd / "devptdoff.diag.txt").read_bytes()
    ok("two-plane TD FP64 diag byte-identical, device twin on vs off", same, "True", same)
    run(args.exe, wd, "devp_tdpr.in", tdn.format(lat="devp_c.bmad", root="devptdpr", extra=DEV), threads="8")
    for comp, band in (("x", E2E_CEIL), ("y", E2E_CEIL)):
        pc = tp.dump_power(wd, "devptdoff-final.wf.h5", comp)
        pd = tp.dump_power(wd, "devptdpr-final.wf.h5", comp)
        ok(f"crossed TD production P{comp} vs CPU (slippage live)", f"{abs(pd-pc)/pc:.3e}",
           f"<= {band:.1e}", abs(pd - pc) / pc <= band,
           num=abs(pd - pc) / pc, bound=band, kind="max")

    migration(args, wd, exe)
    reproducible(args, wd, exe)
    precision(args, wd, exe)
    quantization(args, wd, exe)

    # 7g. The combination is refused for the device as for the CPU.
    r = run(args.exe, wd, "devp_rh.in", pn.format(lat="devp_c.bmad", root="devprh",
            extra=DEV + '  harmonics = 1, 3\n'), expect_fail=True)
    refused = r.returncode != 0 and "TWO LIVE POLARIZATIONS" in r.stdout
    ok("refused: harmonics together with two live polarizations", refused, "True", refused)


# Three things cannot match and none of them is physics. The resolved parameters and
# the echoed input carry the out_root, which differs by construction, and the timestamp
# is the wall clock, which differs whenever two runs land in different seconds.
# Everything else is compared, the field and beam records among them, and the Bmad
# versions and the lattice source are compared and do match.
REPRO_SKIP = ("params", "meta/input_echo", "meta/timestamp")


# The four configurations the recorded ceilings do not isolate, each read as the
# lockstep source row: the device's deposit against the same deposit in FP64. The
# fixed-point accumulator has to be no worse than the float one it replaced on every
# one, and the numbers beside them are the float deposit's on these same input files.
# Cancelling phases come first because a quiet start leaves the slice phasor at
# roundoff, which is where a coarse quantum would show before anywhere else.

PRECISION_CASES = (
    ("cancelling phases", "  beamlet_size = 8\n", 1.464e-05),
    ("unequal weights", "  beamlet_size = 8\n  split_weights = T\n", 1.465e-05),
    ("charge in few cells", "  beam.a_norm_emit = 4e-9\n"
                            "  beam.b_norm_emit = 4e-9\n", 6.889e-06),
    ("the largest load", "  beam.n_particle = 32768\n", 2.507e-06),
    # The load at which performance.md gives the deposit's cost. It carries no
    # float-deposit number, since the float deposit was not run at this load and no
    # longer exists to run, so it is held to the recorded ceiling and to the headroom
    # instead. It makes the accuracy cases reach the load the cost table uses.
    ("the performance load", "  beam.n_particle = 131072\n", None),
)

# The float deposit is not reproducible, so its numbers above carry their own scatter and
# the comparison allows five percent over them. The fixed-point accumulator gives
# 1.462e-05, 1.462e-05, 6.948e-06 and 2.506e-06: the same levels, two of them slightly
# better and one nine parts in a thousand worse (validation.md, val-device).
PRECISION_MARGIN = 1.05

# The quantum has to sit this far below the FP32 spacing of the deposit's own bound, in
# bits, on every input file the precision cases run. The scale is derived per run and the run
# prints the headroom it reached, so this is read back rather than assumed. The setup
# refusal requires only 8, which is where fixed point stops beating the float
# accumulation it replaced. The derivation delivers 37 bits. Holding the checks to 24
# means a derivation that quietly lost precision is caught here rather than passing on
# a combined source row.
PRECISION_BITS_MIN = 24

# What the quantum costs, measured rather than bounded (validation.md, val-device).
# run.device_dep_mutate shifts the derived scale by whole bits and touches nothing
# else, so coarsening the quantum by a known factor and watching the source row separates
# the quantization from the phase and arithmetic error beside it. On the cancelling-phase
# input file the row is 1.4673e-05 at the derived scale and the same to every digit 24 bits
# coarser. It first moves at 26 (1.0100x), and reaches 1.0505x at 29, which is the last
# shift the headroom refusal allows. So the 8-bit floor sits where quantization starts to
# cost a few percent, and at the 37 bits the derivation delivers it costs nothing that
# the row can see.
QUANT_SHIFT = -24          # 16 million times the quantum, where the row does not move
QUANT_UNMOVED = 1.001      # the row may not move by even a part in a thousand there
QUANT_REFUSED = -30        # one bit past the floor, where the run must be refused
BITS_RE = re.compile(r"(\d+) bits below the FP32 spacing")


def precision(args, wd, exe):
    """
    What the fixed-point deposit costs in accuracy, which is nothing.

    The lockstep source row is the device's deposit against the same deposit carried in
    FP64 on the host, so it measures the accumulator directly rather than through a
    power. Each case is asserted twice: under the recorded source ceiling, and no worse
    than the float accumulator on the same input file, which is the number carried beside it.

    A fixed-point accumulator can be reproducible and still be worse, by taking a scale
    that keeps the bound but quantizes away what matters. Thirty-two bits do exactly that
    here and lose to the float deposit by a factor of seven. These four cases show that
    the sixty-four-bit accumulator does not.

    The source row is a combined number, carrying the FP32 arithmetic and the phase error
    along with the quantization, so on its own it cannot say which of the three moved. The
    quantization is separated here by asserting the property it depends on directly, per
    input file: the run derives its own scale and prints how far the quantum sits below the FP32
    spacing of the bound it was chosen against, and that headroom is read back and required
    to be at least PRECISION_BITS_MIN. A quantum that far under the spacing cannot be what
    a source row of 1e-5 is made of. The two together are the separation: the headroom
    bounds the quantization per input file, and the source row measures everything at once.
    """
    print("== 10. the deposit's accuracy against the CPU's FP64 deposit ==")
    for name, extra, before in PRECISION_CASES:
        root = "devprec" + name.split()[0][:4]
        r = run(exe, wd, f"{root}.in", BASE.format(
            root=root, extra=DEV + '  run.fp32_check = "lockstep"\n' + extra))
        if r.returncode != 0:
            ok(f"deposit accuracy, {name}", "run failed", "exit 0", False)
            continue
        v = summary(wd, root)["source"]
        if before is None:
            ok(f"deposit accuracy, {name}", f"{v:.3e}",
               f"<= {CEIL['source']:.1e}", v <= CEIL["source"],
               num=v, bound=CEIL["source"], kind="max")
        else:
            ok(f"deposit accuracy, {name}", f"{v:.3e}",
               f"<= {CEIL['source']:.1e} and <= the float deposit's {before:.2e}",
               v <= CEIL["source"] and v <= before * PRECISION_MARGIN,
               num=v, bound=min(CEIL["source"], before * PRECISION_MARGIN), kind="max")
        m = BITS_RE.search(r.stdout)
        bits = int(m.group(1)) if m else -1
        ok(f"deposit quantum below the FP32 spacing of its bound, {name}",
           f"{bits} bits" if m else "the run printed no headroom",
           f">= {PRECISION_BITS_MIN} bits", bits >= PRECISION_BITS_MIN,
           num=bits, bound=PRECISION_BITS_MIN, kind="min")


def quantization(args, wd, exe):
    """
    What the quantum contributes to the measured deposit error, which is nothing.

    The source row carries the phase error, the FP32 arithmetic and the quantization
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
        extra = DEV + '  run.fp32_check = "lockstep"\n  beamlet_size = 8\n'
        if shift:
            extra += f"  run.device_dep_mutate = {shift}\n"
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
           v <= base * QUANT_UNMOVED,
           num=v, bound=base * QUANT_UNMOVED, kind="max")

    root = "devqref"
    extra = (DEV + '  run.fp32_check = "lockstep"\n  beamlet_size = 8\n'
             f"  run.device_dep_mutate = {QUANT_REFUSED}\n")
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
    Two runs of one input file, bit for bit.

    The deposit accumulates in fixed point so that this holds. Integer addition is
    associative and commutative where float addition is neither, so the order in which
    threads reach a cell cannot reach the answer. This input file is the time-dependent window
    with shot noise, whose eight slices of 1024 macroparticles contend for cells on every
    step. With a float deposit, 54 of these 127 arrays differ with the filter on and 54
    with it off, the field power by 2.0e-07 and 1.5e-07 and the exit power by 1.6e-08 and
    2.2e-08 (validation.md, val-device). The level below is exact equality, asserted
    rather than toleranced.

    A small dark input file is no test, since two runs of a float deposit can already agree
    on one. This input file carries charge, noise and slippage.

    Both filter settings run. They read the source at different points, since with the
    filter on the first of its four source passes converts the accumulator and with it
    off nothing transforms the source and the solve's last pass converts as it reads.
    """
    print("== 9. two runs of one input file, bit for bit ==")
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
           "every array identical", same_keys and not differ,
           num=len(differ) if same_keys else None, bound=0, kind="max")


def migration(args, wd, exe):
    """Slice migration with the device, on check_migration's input files at grid 64."""
    print("== slice migration with the device ==")
    base = cm.BASE.replace("&end\n", "{extra}&end\n")
    # check_migration's heavy input file is a dark quiet start, which suits its accounting and
    # not the instrument. A slice whose field sits at FP64 roundoff normalizes the source
    # and field rows to nothing (validation.md, val-device), and every slippage feeds the
    # tail slice exactly such a field. A seed fills the window and shot noise fills each
    # fresh tail slice, as the time-dependent lockstep input file above does, and the 60 m_e c^2
    # energy spread migrates just as heavily.
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
       "> 10000 moves and drops", moved > 10000 and sum(q for _, q in drops) > 0,
       num=moved, bound=10000, kind="min")
    ok("migration conservation on the device, worst violation", f"{worst:.2e}", "< 1e-10", worst < 1e-10,
       num=worst, bound=1e-10, kind="max")
    ok("migration phase continuity on the device, worst deviation", f"{bdev:.2e}", "< 1e-10", bdev < 1e-10,
       num=bdev, bound=1e-10, kind="max")
    slen = 2 * np.pi * 3
    n_out = sum(int(np.sum((sl["theta"] < -1e-9) | (sl["theta"] >= slen + 1e-9)))
                for sl in beamio.read_slices(wd / "devmd-final.beam.h5", 1e-10, 3e-10))
    ok("migration window residency on the device", f"{n_out} outside", "0", n_out == 0,
       num=n_out, bound=0, kind="max")
    nslice = 16
    pc = diag_power(wd, "devmc", nslice)
    pd = diag_power(wd, "devmd", nslice)
    ok("heavy-migration window power, device vs CPU", f"{abs(pd-pc)/pc:.3e}", f"<= {MIG_E2E_CEIL:.1e}",
       abs(pd - pc) / pc <= MIG_E2E_CEIL,
       num=abs(pd - pc) / pc, bound=MIG_E2E_CEIL, kind="max")
    grown = "particle capacity grown" in r.stdout
    ok("a migrated slice outgrew the rectangle and the seam grew it", grown, "True", grown)

    # 8b. Lockstep across migration: the twin re-stages every step over the changed
    # fills, inside the ceilings, and the FP64 path is untouched.
    run(args.exe, wd, "devm_l.in", heavy_base.format(root="devml", mig="T",
        extra=DEV + '  run.fp32_check = "lockstep"\n', **heavy), threads="8")
    s = summary(wd, "devml")
    for q in CEIL:
        ok(f"migrating lockstep worst_{q}", f"{s[q]:.3e}", f"<= {CEIL[q]:.1e}", s[q] <= CEIL[q],
           num=s[q], bound=CEIL[q], kind="max")
    same = (wd / "devml.diag.txt").read_bytes() == (wd / "devmc.diag.txt").read_bytes()
    ok("migrating FP64 diag byte-identical, device twin on vs off", same, "True", same)

    # 8c. The no-op, exact. The frozen input file is dark, so its power is FP32 noise and says
    # nothing. The read-back beam is what carries the check. With zero moves the
    # migrate = T run differs from migrate = F only in where the readback happened, so the
    # final beam dumps must agree exactly. The fixed-point deposit leaves two migrate = F
    # runs of this input file identical, so both comparisons below are equality.
    run(args.exe, wd, "devm_f1.in", base.format(root="devmf1", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_f2.in", base.format(root="devmf2", mig="F", extra=DEV, **frozen), threads="8")
    run(args.exe, wd, "devm_t.in", base.format(root="devmt", mig="T", extra=DEV, **frozen), threads="8")
    moved = cm.read_migration_file(wd, "devmt")[0]
    z1 = beam_z(wd / "devmf1-final.beam.h5")
    z2 = beam_z(wd / "devmf2-final.beam.h5")
    zt = beam_z(wd / "devmt-final.beam.h5")
    ok("migration no-op on the device: moves", moved, "0", moved == 0,
       num=moved, bound=0, kind="max")
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


def dumps_identical(fa, fb):
    """Whether two openPMD dumps hold the same datasets to the bit, meta excluded."""
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names = []
        a.visititems(lambda n, o: names.append(n)
                     if isinstance(o, h5py.Dataset) and not n.startswith("meta/") else None)
        return all(n in b and np.array_equal(a[n][()], b[n][()]) for n in names)


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
    """check_harmonics' planar input file made time dependent, with shot noise, at grid 64: one
    Genesis reference, its starting state imported by the CPU and the device."""
    (wd / "devg_planar.lat").write_text(ch.GENESIS_LAT)
    genesis_input = ch.GENESIS_INPUT.replace("ngrid=255", "ngrid=64").replace("shotnoise=0\n", "shotnoise=true\n")
    genesis_input = genesis_input.replace("lattice=planar.lat", "lattice=devg_planar.lat").replace("rootname=H3", "rootname=DEVG")
    genesis_input = genesis_input.replace("&field\npower=5e3", "&time\nslen=8e-10\nsample=1\n&end\n\n&field\npower=5e3")
    genesis_input = genesis_input.replace("H3-initial", "DEVG-initial").replace("H3-final", "DEVG-final")
    (wd / "DEVG.in").write_text(genesis_input)
    r = subprocess.run([args.genesis, "DEVG.in"], cwd=wd, capture_output=True, text=True,
                       env=dict(os.environ, FI_PROVIDER="tcp"))
    if r.returncode != 0:
        print(f"FAIL: genesis exited {r.returncode}:\n{r.stdout[-2000:]}\n{r.stderr[-500:]}")
        sys.exit(1)
    for src, dst in (("DEVG-initial.par.h5", "DEVG-initial.beam.h5"),
                     ("DEVG-initial.fld.h5", "DEVG-initial.wf.h5")):
        convert_genesis.to_openpmd(wd / src, wd / dst)
    imp = ch.NML_IMPORT.replace('transport_model = "genesis"', 'transport_model = "bmad"')
    curves = {}
    for tag, extra in (("cpu", ""), ("dev", DEV)):
        name = f"devg_{tag}"
        (wd / (name + ".toml")).write_text(to_groups(imp.format(lat="dev_planar.bmad", root=name,
            beam="DEVG-initial.beam.h5", field="DEVG-initial.wf.h5", extra=extra)))
        r = subprocess.run([str(exe), name + ".toml"], cwd=wd, capture_output=True, text=True,
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
    ok("device vs Genesis, TD planar with shot noise, P1", f"{d1:.3e}", f"<= {E2E_CEIL:.1e}", d1 <= E2E_CEIL,
       num=d1, bound=E2E_CEIL, kind="max")
    ok("device vs Genesis, TD planar with shot noise, P3", f"{d3:.3e}", f"<= {H_E2E_CEIL:.1e}", d3 <= H_E2E_CEIL,
       num=d3, bound=H_E2E_CEIL, kind="max")


if __name__ == "__main__":
    sys.exit(main())
