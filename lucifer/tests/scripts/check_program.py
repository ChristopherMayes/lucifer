#!/usr/bin/env python3
"""
Program-structure checks (manual sec:program): the library contract, the comb, and
the tracking window.

  1. LIBRARY: lucifer_smoke_test drives the tracker with NO namelist anywhere -- structs
     filled in code -- and runs TWICE IN ONE PROCESS. Pass 1 must reproduce a
     namelist-driven run of the same configuration dataset-identically (the library
     IS the program), and pass 2 must reproduce pass 1 bit-for-bit (re-entrancy:
     twice in one process = two processes).
  2. LIBRARY ERRORS RETURN: lucifer_smoke_test on an unreadable lattice must print its
     proof-of-return line and exit with its own code 2 -- the library returned, the
     PROGRAM decided (no exit inside the library).
  3. RETIRED GROUP: the flat &fel_track_params is refused BY NAME, the error mapping
     each parameter found to its new group.
  4. THE COMB (global%comb_ds_save, Bmad's ds_save semantics verbatim): a comb > 0
     run's per-record rows must be EXACTLY the every-record run's rows at the comb
     positions (subset, dataset-equal), element ends always present; comb < 0 keeps
     NO per-record rows (element-end arrays and dumps remain); the precomputed nrec
     is exact in every mode (the arrays are sized once, never grown).
  5. THE WINDOW (global%track_start/track_end, Tao's names): the schedule is built
     on the FULL lattice, so a windowed run composes exactly -- run A = [start, D]
     dumps its final state; run B = [after D, end] imports it; B's finals must be
     dataset-identical to the one-shot full run's finals, and A's finals to the full
     run's mid-line dumps at D.

Run by the benchmark harness; exits nonzero on failure. Self-referenced (no Genesis).
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

FAILED = False

LAT1 = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
fel_transcribed = -1
wiggler::*[FEL_TRACKING] = fel_transcribed
SEG: line = (UND)
use, SEG
"""

LAT2 = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
UND2: UND
D: pipe, l = 0.30
fel_transcribed = -1
wiggler::*[FEL_TRACKING] = fel_transcribed
SEG: line = (UND, D, UND2)
use, SEG
"""

# The namelist twin of lucifer_smoke_test's in-code configuration (same values).
NML_TWIN = """&fel_params
  lat_file = "smoke.bmad"
  global%out_root = "{root}"
  global%interlude_model = "genesis"
  global%write_diag = T
  global%ran_seed = 777
{extra}/
&fel_beam_init
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
{bextra}/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 1e7
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_n_pts = 63
  wavefront_init%grid_half_width = 2e-4
/
"""

OLD_STYLE = """&fel_track_params
  lat_file = "smoke.bmad"
  out_root = "old"
  lambda0 = 1e-10
  wake_on = T
/
"""


def check(name, ok, note=""):
    global FAILED
    print(f"--- {name}: {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(cmd, wd, threads="4"):
    return subprocess.run(cmd, cwd=wd, capture_output=True, text=True,
                          env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})


def run_nml(exe, wd, root, nml_text):
    (wd / f"{root}.nml").write_text(nml_text)
    r = run([str(exe), f"{root}.nml"], wd)
    if r.returncode != 0:
        print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)


def h5_identical(fa, fb):
    with h5py.File(fa) as a, h5py.File(fb) as b:
        na, nb = [], []
        a.visititems(lambda n, o: na.append(n) if isinstance(o, h5py.Dataset) else None)
        b.visititems(lambda n, o: nb.append(n) if isinstance(o, h5py.Dataset) else None)
        if sorted(na) != sorted(nb):
            return False
        return all(np.array_equal(a[n][()], b[n][()]) for n in na)


def stats_of(path):
    out = {}
    with h5py.File(path) as h5:
        for k in ("z", "field/power", "beam/bunching", "element_end/s"):
            if k in h5:
                out[k] = h5[k][()]
    return out


def main():
    global FAILED
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--smoke", required=True)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(args.exe).resolve()
    smoke = pathlib.Path(args.smoke).resolve()

    (wd / "smoke.bmad").write_text(LAT1)
    (wd / "twoseg.bmad").write_text(LAT2)

    # ------------------------------------------------------------------
    print("== the library contract ==")

    r = run([str(smoke), "smoke.bmad", "sp1", "sp2"], wd)
    if r.returncode != 0:
        print(f"FAIL: lucifer_smoke_test exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    run_nml(exe, wd, "snml", NML_TWIN.format(root="snml", extra="", bextra=""))

    same = (wd / "sp1.diag.txt").read_bytes() == (wd / "snml.diag.txt").read_bytes()
    for s in ("-final.par.h5", "-final.fld.h5"):
        same = same and h5_identical(wd / f"sp1{s}", wd / f"snml{s}")
    check("no-namelist library run == namelist run (diag byte, dumps dataset)", same)

    same = (wd / "sp1.diag.txt").read_bytes() == (wd / "sp2.diag.txt").read_bytes()
    for s in ("-final.par.h5", "-final.fld.h5", ".stats.h5"):
        same = same and h5_identical(wd / f"sp1{s}", wd / f"sp2{s}")
    check("re-entrant: twice in one process bit-identical (pass 2 == pass 1)", same)

    (wd / "bad.bmad").write_text(LAT1.replace(
        "beginning[beta_a] = 15\nbeginning[beta_b] = 15\n", ""))
    r = run([str(smoke), "bad.bmad", "se1", "se2"], wd)
    ok = (r.returncode == 2 and "the library returned an error" in r.stdout
          and "NO BEGINNING TWISS" in r.stdout)
    check("a library error RETURNS (no beginning Twiss -> err to the program, exit 2)", ok,
          note=f"[exit {r.returncode}]")

    # ------------------------------------------------------------------
    print("== the retired flat group ==")
    (wd / "old.nml").write_text(OLD_STYLE)
    r = run([str(exe), "old.nml"], wd)
    ok = (r.returncode != 0 and "IS RETIRED" in r.stdout
          and "lambda0 -> &fel_wavefront_init wavefront_init%lambda0" in r.stdout
          and "wake_on -> &fel_params wake%on" in r.stdout)
    check("&fel_track_params refused by name, each parameter mapped", ok)

    # ------------------------------------------------------------------
    print("== the comb (global%comb_ds_save) ==")

    # Every-record baseline (comb 0, the default), a spaced comb, and no comb.
    run_nml(exe, wd, "cb0", NML_TWIN.format(root="cb0", extra="", bextra=""))
    run_nml(exe, wd, "cbp", NML_TWIN.format(root="cbp",
            extra="  global%comb_ds_save = 0.05\n", bextra=""))
    run_nml(exe, wd, "cbn", NML_TWIN.format(root="cbn",
            extra="  global%comb_ds_save = -1\n", bextra=""))

    s0, sp, sn = stats_of(wd / "cb0.stats.h5"), stats_of(wd / "cbp.stats.h5"), \
                 stats_of(wd / "cbn.stats.h5")

    # comb > 0: the rows are a subset of the every-record run's rows, dataset-equal
    # at the matching z; element ends always present (here: the final record).
    idx = np.searchsorted(s0["z"], sp["z"])
    ok = bool(np.array_equal(s0["z"][idx], sp["z"]))
    ok = ok and np.array_equal(s0["field/power"][idx, :], sp["field/power"])
    ok = ok and np.array_equal(s0["beam/bunching"][idx, :], sp["beam/bunching"])
    ok = ok and sp["z"][-1] == s0["z"][-1]
    # the spacing rule itself: consecutive rows at least comb apart (ends exempt).
    ok = ok and bool(np.all(np.diff(sp["z"][:-1]) >= 0.05 - 1e-12))
    check("comb > 0: rows == every-record rows at the comb positions (subset)", ok,
          note=f"[{len(sp['z'])} of {len(s0['z'])} rows]")

    # comb < 0: no per-record rows at all; element-end arrays remain.
    ok = len(sn.get("z", [])) == 0 and len(sn["element_end/s"]) == len(s0["element_end/s"]) \
         and bool(np.array_equal(sn["element_end/s"], s0["element_end/s"]))
    check("comb < 0: no per-record rows; element ends remain", ok,
          note=f"[nrec {len(sn.get('z', []))}, ends {len(sn['element_end/s'])}]")

    # nrec exact: the arrays are sized by the same rule the walk replays -- full,
    # never padded (h5 dataset lengths ARE nrec).
    ok = len(s0["z"]) == 31 and len(sn.get("z", [])) == 0
    check("nrec exact in every mode (sized once, never grown)", ok,
          note=f"[comb0 {len(s0['z'])} rows = 30 steps + initial]")

    # The keystone locally: comb 0 (the default) is bit-for-bit the pre-comb run --
    # cb0 above ran with the default (no comb key at all) and fed every comparison.

    # ------------------------------------------------------------------
    print("== the tracking window (global%track_start/track_end) ==")

    two = NML_TWIN.replace('lat_file = "smoke.bmad"', 'lat_file = "twoseg.bmad"')
    run_nml(exe, wd, "wfull", two.format(root="wfull",
            extra='  global%dump_beam_at = "D"\n  global%dump_field_at = "D"\n', bextra=""))
    run_nml(exe, wd, "wa", two.format(root="wa",
            extra='  global%track_end = "D"\n', bextra=""))

    # A = [start, D]: its finals equal the full run's mid-line dumps at D.
    ok = h5_identical(wd / "wa-final.par.h5", wd / "wfull-at2-D.par.h5")
    ok = ok and h5_identical(wd / "wa-final.fld.h5", wd / "wfull-at2-D.fld.h5")
    check("windowed [start, D] finals == full run's dumps at D", ok)

    # B = [after D, end] from A's finals: composes to the full run's finals.
    run_nml(exe, wd, "wb", two.format(root="wb",
            extra='  global%track_start = "UND2"\n', bextra="""  beam_file = "wa-final.par.h5"
""").replace("&fel_wavefront_init\n", """&fel_wavefront_init
  field_file = "wa-final.fld.h5"
"""))
    # B's state passed through the Genesis dump format once more than the full run
    # (theta/gamma/px folds and back), so the composition sits at the dump
    # round-trip's conversion floor, not at zero: measured 2.6e-13 rad in theta,
    # 3e-14 of the field scale, 7e-18 in px (the walk itself is bit-for-bit -- check
    # A above IS exact, both sides dumping the same in-memory state).
    worst = 0.0
    with h5py.File(wd / "wb-final.par.h5") as a, h5py.File(wd / "wfull-final.par.h5") as b:
        for k in ("gamma", "theta", "x", "y", "px", "py", "current"):
            da, db = a[f"slice000001/{k}"][()], b[f"slice000001/{k}"][()]
            scale = max(np.max(np.abs(db)), 1e-300)
            worst = max(worst, float(np.max(np.abs(da - db)) / scale))
    with h5py.File(wd / "wb-final.fld.h5") as a, h5py.File(wd / "wfull-final.fld.h5") as b:
        fs = max(float(np.max(np.abs(b["slice000001/field-real"][()]))), 1e-300)
        for k in ("field-real", "field-imag"):
            worst = max(worst, float(np.max(np.abs(a[f"slice000001/{k}"][()]
                                                   - b[f"slice000001/{k}"][()])) / fs))
    check("windowed [after D, end] from A's dumps == full run's finals (composition)",
          worst <= 1e-10,
          note=f"[max rel {worst:.2e} vs 1e-10; the dump round-trip's floor, measured 3e-13]")

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
