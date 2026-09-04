#!/usr/bin/env python3
"""
Program-structure checks (doc/user-guide.md): the library contract, the comb, and
the tracking window.

  1. Library: lucifer_smoke_test drives the tracker with NO namelist anywhere -- structs
     filled in code -- and runs twice in one process. Pass 1 must reproduce a
     namelist-driven run of the same configuration dataset-identically (the library
     IS the program), and pass 2 must reproduce pass 1 bit-for-bit (re-entrancy:
     twice in one process = two processes).
  2. Library errors return: lucifer_smoke_test on an unreadable lattice must print its
     proof-of-return line and exit with its own code 2 -- the library returned, the
     program decided (no exit inside the library).
  3. Retired group: the flat &fel_track_params is refused. Fortran ignores an
     unknown namelist group in silence, so without the refusal a stale deck would run
     on defaults.
  4. The comb (global%comb_ds_save, Bmad's ds_save semantics verbatim): a comb > 0
     run's per-record rows must be exactly the every-record run's rows at the comb
     positions (subset, dataset-equal), element ends always present; comb < 0 keeps
     NO per-record rows (element-end arrays and dumps remain); the precomputed nrec
     is exact in every mode (the arrays are sized once, never grown).
  5. The window (global%track_start/track_end, Tao's names): the schedule is built
     on the full lattice, so a windowed run composes exactly -- run A = [start, D]
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

import beamio
import fieldio
from read_stats import read_stats, same_data

FAILED = False

LAMBDA0 = 1e-10          # the wavelength both decks below state

LAT1 = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = fel_averaged, ds_step = 0.015
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
      tracking_method = fel_averaged, ds_step = 0.015
UND2: UND
D: pipe, l = 0.30
SEG: line = (UND, D, UND2)
use, SEG
"""

# The namelist twin of lucifer_smoke_test's in-code configuration (same values).
NML_TWIN = """&fel_params
  lat_file = "smoke.bmad"
  global%out_root = "{root}"
  global%interlude_model = "genesis"
  global%transport_model = "genesis"
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
    # meta/ is excluded, deliberately. Provenance is datasets rather than attributes,
    # for HDF5's 64 kB attribute cap (fel-physics.md sec-meta), and meta/timestamp differs
    # between any two runs by construction. Nothing in meta/ is physics. Before the
    # move, the exclusion existed only by the accident of being attributes.
    with h5py.File(fa) as a, h5py.File(fb) as b:
        na, nb = [], []
        keep = lambda n, o: isinstance(o, h5py.Dataset) and not n.startswith("meta/")
        a.visititems(lambda n, o: na.append(n) if keep(n, o) else None)
        b.visititems(lambda n, o: nb.append(n) if keep(n, o) else None)
        if sorted(na) != sorted(nb):
            return False
        return all(same_data(a[n][()], b[n][()]) for n in na)


def stats_of(path):
    """The comb-relevant content of a stats file: path length, two per-record
    quantities, and the element-end mask that selects the ends out of the record axis."""
    with read_stats(path) as st:
        return {"s": st.s, "power": st["field/total/power"],
                "bunching": st["beam/slice/bunching"], "at_end": st.at_end,
                "ix_ele": st.ix_ele}


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
    for s in ("-final.beam.h5", "-final.wf.h5"):
        same = same and h5_identical(wd / f"sp1{s}", wd / f"snml{s}")
    check("no-namelist library run == namelist run (diag byte, dumps dataset)", same)

    same = (wd / "sp1.diag.txt").read_bytes() == (wd / "sp2.diag.txt").read_bytes()
    for s in ("-final.beam.h5", "-final.wf.h5", ".stats.h5"):
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
    ok = (r.returncode != 0 and "IS NOT AN INPUT GROUP" in r.stdout
          and "&fel_params" in r.stdout and "&fel_beam_init" in r.stdout
          and "&fel_wavefront_init" in r.stdout)
    check("&fel_track_params refused, naming the three groups", ok)

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

    # First, what the comparison below rests on. The record number is the axis and s is
    # a variable on it (fel-physics.md sec-stats), so matching rows BY s is only legitimate while
    # s does not repeat, and it can: a zero-length element that applies a wake kick sits
    # at the plane the element before it ended on. Every repeat must therefore straddle
    # an element boundary. A repeat inside one element would be a defect in the walk, and
    # every s-keyed comparison in this file would then silently pick the wrong row.
    ds = np.diff(s0["s"])
    dup = np.flatnonzero(ds == 0)
    ok = bool(np.all(ds >= 0))
    ok = ok and not any(s0["ix_ele"][i] == s0["ix_ele"][i + 1] for i in dup)
    check("the record axis: s non-decreasing, every repeat at an element boundary", ok,
          note=f"[{len(dup)} repeats in {len(s0['s'])} rows over 3 elements]")

    # comb > 0: the rows are a subset of the every-record run's rows, dataset-equal
    # at the matching s. Element ends always present (here: the final record).
    idx = np.searchsorted(s0["s"], sp["s"])
    ok = bool(np.array_equal(s0["s"][idx], sp["s"]))
    ok = ok and np.array_equal(s0["power"][idx, :], sp["power"])
    ok = ok and np.array_equal(s0["bunching"][idx, :], sp["bunching"])
    ok = ok and sp["s"][-1] == s0["s"][-1]
    # the spacing rule itself: consecutive rows at least comb apart (ends exempt).
    ok = ok and bool(np.all(np.diff(sp["s"][:-1]) >= 0.05 - 1e-12))
    check("comb > 0: rows == every-record rows at the comb positions (subset)", ok,
          note=f"[{len(sp['s'])} of {len(s0['s'])} rows]")

    # comb < 0: the element ends, and nothing else. Bmad's comb semantics drop the comb
    # there; this tracker always keeps the element ends, because the stats file carries
    # one record axis and marks the ends inside it (fel-physics.md sec-stats). So a comb < 0 run
    # is a file whose every record is an element end, at the same positions the
    # every-record run put them.
    ok = bool(np.all(sn["at_end"])) and len(sn["s"]) == int(s0["at_end"].sum())
    ok = ok and bool(np.array_equal(sn["s"], s0["s"][s0["at_end"]]))
    ok = ok and bool(np.array_equal(sn["power"], s0["power"][s0["at_end"]]))
    check("comb < 0: the rows are exactly the element ends", ok,
          note=f"[{len(sn['s'])} rows, all element ends]")

    # nrec exact: the arrays are sized by the same rule the walk replays -- full,
    # never padded (h5 dataset lengths are nrec).
    ok = len(s0["s"]) == 31 and len(sn["s"]) == int(s0["at_end"].sum())
    check("nrec exact in every mode (sized once, never grown)", ok,
          note=f"[comb0 {len(s0['s'])} rows = 30 steps + initial, "
               f"comb<0 {len(sn['s'])} ends]")

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
    ok = h5_identical(wd / "wa-final.beam.h5", wd / "wfull-at2-D.beam.h5")
    ok = ok and h5_identical(wd / "wa-final.wf.h5", wd / "wfull-at2-D.wf.h5")
    check("windowed [start, D] finals == full run's dumps at D", ok)

    # B = [after D, end] from A's finals: composes to the full run's finals.
    run_nml(exe, wd, "wb", two.format(root="wb",
            extra='  global%track_start = "UND2"\n', bextra="""  beam_file = "wa-final.beam.h5"
""").replace("&fel_wavefront_init\n", """&fel_wavefront_init
  field_file = "wa-final.wf.h5"
"""))
    # B's state passed through the dump format once more than the full run (a pack and an
    # unpack of every coordinate), so the composition sits at the dump round-trip's
    # conversion floor, not at zero: measured 2.6e-13 rad in theta, 3e-14 of the field
    # scale, 7e-18 in px (the walk itself is bit-for-bit -- check A above IS exact, both
    # sides dumping the same in-memory state).
    #
    # theta is the sharpest column here. A dump carries the particle lag and the reader
    # restarts the reference phase at zero, so a restart reproduces the absolute phase only
    # if the writer folded the reference in. Nothing else in this check can see that, and
    # the beam's phase against the field's is what the next segment's gain is made of.
    worst = 0.0
    pa = beamio.read_slices(wd / "wb-final.beam.h5", LAMBDA0, LAMBDA0)[0]
    pb = beamio.read_slices(wd / "wfull-final.beam.h5", LAMBDA0, LAMBDA0)[0]
    for k in ("gamma", "theta", "x", "y", "px", "py", "weight"):
        scale = max(float(np.max(np.abs(pb[k]))), 1e-300)
        worst = max(worst, float(np.max(np.abs(pa[k] - pb[k])) / scale))
    fa = fieldio.read_field(wd / "wb-final.wf.h5")["u"]
    fb = fieldio.read_field(wd / "wfull-final.wf.h5")["u"]
    worst = max(worst, float(np.max(np.abs(fa - fb))) / max(float(np.max(np.abs(fb))), 1e-300))
    check("windowed [after D, end] from A's dumps == full run's finals (composition)",
          worst <= 1e-10,
          note=f"[max rel {worst:.2e} vs 1e-10; the dump round-trip's floor, measured 3e-13]")

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
