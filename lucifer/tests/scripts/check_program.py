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
  6. Committed prose cites committed artifacts: no file this program ships references
     the design repository's numbered findings, which a reader of this tree cannot follow.
  7. A restart across migration and the chamber wake: a beam that migrates and drops
     through three undulators, checkpointed after two events with a residual. Every later
     decision clears its recorded margin by 100x the phase error the continuation has
     propagated, and then ids, populations and membership are identical at every frame,
     events match by position and moves, the wake table is rebuilt from the surviving
     charge and refreshed at every later event, and moves, dropped charge and escaped
     energy since the checkpoint agree at full precision. A separate case shifts the
     restarted phases past the smallest margin and shows one decision flip.
  8. A restart across a Bmad element wake: checkpoints immediately before and after the
     passage that applies the kick, which the check reads from the runs rather than
     assuming, once on the element itself and once on a wake a superimposed marker has
     split onto a lord. Each is judged against a wake-off history run from the same
     initial state, and the wake's own effect is the scale: a kick omitted or applied
     twice would move it by order one.
  9. A continuation across an unaveraged segment boundary: the entry and exit
     conversions are scratch of the segment, so an element-boundary dump carries the
     averaged chart and the plane field, and continuations from both boundaries hold at
     the floor. The one before the segment holds because the checkpoint group carries
     phi0 and every z as the tracker holds them, rebuilt bit for bit before any advance,
     where the folded records alone left the entry state bit-identical yet the
     continuation amplified to 6.5e-8. Frames from inside the segment, from inside an
     averaged undulator and from a sliced endpoint inside one carry no group and are
     refused, as are an incomplete or inconsistent group and a continuation point that
     does not follow the checkpoint. The same checkpoint loaded as an initialization
     folds the records and starts the residual at zero.
 10. A continuation across a slippage residual: the window check above is steady state,
     where slippage is a no-op, so a time-dependent pair follows it. The continuous run
     dumps at a pipe whose accumulated slip leaves half a slice of residual, the dump's
     slippageResidual attribute carries it, and the continuation must reproduce the
     continuous run's frames, field slice by slice and particles by id, to the same
     floor. A continuation from the dump with the attribute stripped is refused, and so
     is an impossible residual.

Run by the benchmark harness; exits nonzero on failure. Self-referenced (no Genesis).
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

import beamio
import fieldio
from read_stats import read_stats, same_data

FAILED = False

LAMBDA0 = 1e-10          # the wavelength both decks below state

# What the citation check below walks: the text this program ships.
SOURCE_SUFFIXES = {".f90", ".mm", ".h", ".c", ".py", ".sh", ".md", ".bmad", ".lat", ".in",
                   ".nml", ".yml", ".yaml", ".txt"}

LAT1 = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      fel_method = averaged, ds_step = 0.015
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
      fel_method = averaged, ds_step = 0.015
UND2: UND
D: pipe, l = 0.30
SEG: line = (UND, D, UND2)
use, SEG
"""

# The namelist twin of lucifer_smoke_test's in-code configuration (same values).
# grid_n_pts is 63 here on purpose, and this pair is the only place in the harness that is
# not a power of two. Every other deck runs a power-of-two grid, which is what the Metal
# field solver takes, so a deck moves between the backends unchanged. This pair holds the
# odd-grid path on the CPU, which has to keep working: a user writing a Genesis-style odd
# count gets one, and an imported Genesis field can carry one. The value has to match
# lucifer_smoke_test.f90, since the check below compares the two runs byte for byte.

NML_TWIN = """&fel_params
  lat_file = "smoke.bmad"
  global%out_root = "{root}"
  global%source_filter = F
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
  source_filter = F
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
    # A refusal that quotes a malformed file can carry bytes that are not text, and the
    # check reads the message rather than dying on it.
    return subprocess.run(cmd, cwd=wd, capture_output=True, text=True, errors="replace",
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
            extra='  global%continuation = T\n  global%track_start = "UND2"\n', bextra="""  beam_file = "wa-final.beam.h5"
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

    # ------------------------------------------------------------------
    print("== a restart across a slippage residual (time dependent) ==")

    # The composition check above runs steady state, where slippage is a no-op, so it
    # cannot see the one piece of the field's state a dump did not carry: the slip the
    # record has accumulated and not yet rotated, fel_slip_struct%accuslip. A restart
    # that started it at zero rotated on a different schedule from the run it continued.
    # Two undulators with a pipe between, four wavelengths a slice, a
    # deterministic seed ramped slice by slice so a rotation in the wrong place shows,
    # and the beam and field both imported so no loader noise enters. The residual at
    # the boundary is read from the dump's own slippageResidual attribute. Two pipes: one
    # leaves half a slice of residual, the other leaves a residual near zero. A
    # continuation from a dump stripped of the attribute is refused: the residual is part
    # of the field's state, and a continuation takes nothing from a file that lacks it.

    rs_nml = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%source_filter = F
  global%interlude_model = "genesis"
  global%transport_model = "genesis"
  global%write_diag = T
  global%ran_seed = 777
  slicing%n_wavelength = 4
{extra}/
&fel_beam_init
{beam}/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
{field}/
"""
    gen_beam = """  beam_init%n_particle = 1024
  beam_init%bunch_charge = 2.401661485427e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.2e-9
  beam_init%grid(3)%x_max =  1.2e-9
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
"""
    gen_field = """  wavefront_init%seed_power = 1e7
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_n_pts = 63
  wavefront_init%grid_half_width = 2e-4
"""
    def rs_lat(pipe):
        return LAT2.replace("l = 0.30", f"l = {pipe}").replace("D: pipe", "D: pipe")

    def frames_by_s(root):
        out = {}
        for f in sorted(wd.glob(f"{root}-0*.wf.h5")):
            with h5py.File(f) as h5:
                out[round(float(np.atleast_1d(h5.attrs["sPosition"])[0]), 9)] = f
        return out

    def beam_by_id(path):
        with h5py.File(path) as h5:
            name = sorted(h5["data"].keys())[0]
            sp = h5[f"data/{name}/particles"]
            g = sp[sorted(sp.keys())[0]]
            n = int(np.atleast_1d(g.attrs["numParticles"])[0])
            rec = lambda k: g[k][...] if g[k].shape else np.full(n, g[k][()])
            ids = rec("id")
            order = np.argsort(ids)
            pops = np.atleast_1d(g["particlePatches/numParticles"][...])
            off = np.atleast_1d(g["particlePatches/numParticlesOffset"][...])
            slice_of = np.zeros(n, int)
            for k, (o0, c) in enumerate(zip(off, pops)):
                slice_of[o0:o0 + c] = k
            return dict(id=ids[order], w=rec("weight")[order], t=rec("time")[order],
                        pz=rec("totalMomentum")[order], sl=slice_of[order], pops=pops)

    def compare_runs(full, other):
        """Worst relative differences over the frames the two runs share, and the s of the
        first frame whose field differs by more than the composition floor."""
        ff, fo = frames_by_s(full), frames_by_s(other)
        common = sorted(set(ff) & set(fo))
        worst = dict(field=0.0, theta=0.0, pz=0.0, weight=0.0)
        ids_ok = pops_ok = True
        first = None
        for s in common:
            uf = fieldio.read_field(ff[s])["u"]
            uo = fieldio.read_field(fo[s])["u"]
            d = float(np.max(np.abs(uf - uo))) / max(float(np.max(np.abs(uf))), 1e-300)
            worst["field"] = max(worst["field"], d)
            if first is None and d > 1e-10:
                first = s
            bf = beam_by_id(str(ff[s]).replace(".wf.h5", ".beam.h5"))
            bo = beam_by_id(str(fo[s]).replace(".wf.h5", ".beam.h5"))
            ids_ok = ids_ok and np.array_equal(bf["id"], bo["id"])
            pops_ok = pops_ok and np.array_equal(bf["pops"], bo["pops"])
            for k, key in (("t", "theta"), ("pz", "pz"), ("w", "weight")):
                scale = max(float(np.max(np.abs(bf[k]))), 1e-300)
                worst[key] = max(worst[key], float(np.max(np.abs(bf[k] - bo[k]))) / scale)
        with read_stats(wd / f"{full}.stats.h5") as a, read_stats(wd / f"{other}.stats.h5") as b:
            ia = {round(float(s), 9): i for i, s in enumerate(a.s)}
            pa, pb = a["field/total/power"], b["field/total/power"]
            dp = 0.0
            for j, s in enumerate(b.s):
                i = ia.get(round(float(s), 9))
                if i is not None:
                    dp = max(dp, float(np.max(np.abs(np.asarray(pa[i]) - np.asarray(pb[j]))))
                             / max(float(np.max(np.abs(pa[i]))), 1e-300))
        worst["power"] = dp
        return worst, len(common), ids_ok, pops_ok, first

    def escaped_of(root):
        m = re.findall(r"Escaped\s+h1\s+([0-9.Ee+-]+)", (wd / f"{root}.out").read_text())
        return float(m[0])

    def rs_run(root, nml_text, expect_fail=False, fragment=""):
        (wd / f"{root}.nml").write_text(nml_text)
        r = run([str(exe), f"{root}.nml"], wd)
        (wd / f"{root}.out").write_text(r.stdout + r.stderr)
        if expect_fail:
            ok = r.returncode != 0 and fragment in r.stdout + r.stderr
            return ok
        if r.returncode != 0:
            print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        return True

    def residual_of(path):
        with h5py.File(path) as h5:
            return float(np.atleast_1d(h5.attrs["slippageResidual"])[0]) if "slippageResidual" in h5.attrs else None

    # The starting state: one seeded run through the first undulator, its final dumps
    # taken as the imported beam and, ramped, the imported field.
    (wd / "rs.bmad").write_text(rs_lat(0.30))
    rs_run("rsprep", rs_nml.format(lat="rs.bmad", root="rsprep",
           extra='  slicing%window_length = 4.8e-9\n  global%track_end = "UND"\n',
           beam=gen_beam, field=gen_field))
    shutil.copy(wd / "rsprep-final.wf.h5", wd / "rsramp.wf.h5")
    with h5py.File(wd / "rsramp.wf.h5", "r+") as h5:
        m = h5[fieldio.MESH_PATH]
        u = m["x"][...]
        m["x"][...] = u * (1.0 + 0.5 * np.arange(u.shape[0]) / u.shape[0])[:, None, None]
        if "slippageResidual" in h5.attrs:
            del h5.attrs["slippageResidual"]      # the seed starts every run from zero
    imp = dict(beam='  beam_file = "rsprep-final.beam.h5"\n', field='  field_file = "rsramp.wf.h5"\n')

    for pipe, tag, kind in ((0.30, "rsh", "half a slice"), (0.35, "rsz", "near zero")):
        lat = f"rs_{tag}.bmad"
        (wd / lat).write_text(rs_lat(pipe))
        full, rest, bare, toD = f"{tag}full", f"{tag}rest", f"{tag}bare", f"{tag}toD"
        rs_run(full, rs_nml.format(lat=lat, root=full,
               extra='  global%dump_beam_at = "D"\n  global%dump_field_at = "D"\n  global%dump_at_comb = T\n', **imp))
        rs_run(toD, rs_nml.format(lat=lat, root=toD, extra='  global%track_end = "D"\n', **imp))
        dump_b, dump_f = f"{full}-at2-D.beam.h5", f"{full}-at2-D.wf.h5"
        resid = residual_of(wd / dump_f)
        # The continuation carrying the residual, and one from the same dump with the
        # residual stripped, which is what a file written before the attribute existed is.
        rs_run(rest, rs_nml.format(lat=lat, root=rest,
               extra='  global%continuation = T\n  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{dump_b}"\n', field=f'  field_file = "{dump_f}"\n'))
        shutil.copy(wd / dump_f, wd / f"{tag}bare.wf.h5")
        with h5py.File(wd / f"{tag}bare.wf.h5", "r+") as h5:
            del h5.attrs["slippageResidual"]
        stripped = rs_run(bare, rs_nml.format(lat=lat, root=bare,
               extra='  global%continuation = T\n  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{dump_b}"\n', field=f'  field_file = "{tag}bare.wf.h5"\n'),
               expect_fail=True, fragment="NO slippageResidual")

        w_rest, n_rest, ids_r, pops_r, first_r = compare_runs(full, rest)
        worst_rest = max(w_rest.values())
        if kind == "half a slice":
            check(f"the dump at the pipe carries a residual of half a slice ({pipe} m pipe)",
                  resid is not None and 0.5 < abs(resid) < 3.2,
                  note=f"[slippageResidual {resid} wavelengths of 4 a slice]")
        else:
            check(f"the dump at the pipe carries a residual near zero ({pipe} m pipe)",
                  resid is not None and abs(resid) < 1e-3,
                  note=f"[slippageResidual {resid} wavelengths, pipe slip is whole]")
        check(f"continuation carrying the residual == continuous, {kind}: field, phase, energy, weight, "
              f"power at {n_rest} frames", worst_rest <= 1e-10 and ids_r and pops_r,
              note=f"[worst rel {worst_rest:.2e} vs 1e-10; ids and populations identical: "
                   f"{ids_r and pops_r}; field {w_rest['field']:.1e} theta {w_rest['theta']:.1e} "
                   f"pz {w_rest['pz']:.1e} power {w_rest['power']:.1e}]")
        check(f"a continuation from the dump with the residual stripped is refused, {kind}, and the "
              "message names slippageResidual", stripped)
        e_full, e_toD, e_rest = escaped_of(full), escaped_of(toD), escaped_of(rest)
        check(f"escaped energy: the continuous run's increment past the pipe == the restart's total, "
              f"{kind}", abs((e_full - e_toD) - e_rest) <= 1e-5 * e_full,
              note=f"[continuous {e_full:.6e} J less {e_toD:.6e} J at the pipe = "
                   f"{e_full - e_toD:.3e} J; restart {e_rest:.3e} J; within the footer's 1e-5]")

    # A present but impossible residual is refused, never read as absent: the magnitude
    # can never reach n_wavelength, since the threshold is 0.8 of it.
    shutil.copy(wd / "rshfull-at2-D.wf.h5", wd / "rsbad.wf.h5")
    with h5py.File(wd / "rsbad.wf.h5", "r+") as h5:
        h5.attrs["slippageResidual"] = np.array([7.0])
    ok = rs_run("rsbad", rs_nml.format(lat="rs_rsh.bmad", root="rsbad",
                extra='  global%track_start = "UND2"\n',
                beam='  beam_file = "rshfull-at2-D.beam.h5"\n', field='  field_file = "rsbad.wf.h5"\n'),
                expect_fail=True, fragment="slippageResidual")
    check("a residual of 7 wavelengths at 4 a slice is refused, and the message names slippageResidual", ok)

    # ------------------------------------------------------------------
    print("== a restart across migration, then with the chamber wake ==")

    # The residual case above migrates nothing. Here the beam migrates and drops through
    # three undulators, the checkpoint is the first pipe's end, after two events, and the
    # continuation holds three more. Migration derives a destination from the whole
    # phase, phi0 + ks z/beta, floored over the window's 2 pi n_wavelength, and the dump
    # folds exactly that phase into the file's time, so a decision carries across a
    # restart unless a particle sits nearer a boundary than the error in its
    # reconstructed phase. The migration file records that margin per event, over every
    # particle examined and the dropped ones included, which no frame after the event can
    # show. The margin is measured against the phase error the continuation has actually
    # propagated to the frame before each event, read by id, and not against the file's
    # bound alone. With the margin cleared, exact identity is required: one particle in
    # another slice changes the current every wake is built from, and one drop changes the
    # source charge. A separate case shifts every restarted phase by a little less and a
    # little more than the smallest recorded margin: the first flips nothing, the second
    # flips a decision, which is what the margin means, and it establishes nothing about
    # equivalence.

    mw_lat = LAT2.replace("SEG: line = (UND, D, UND2)",
                          "UND3: UND\nD2: pipe, l = 0.30\nSEG: line = (UND, D, UND2, D2, UND3)")
    (wd / "mw.bmad").write_text(mw_lat)
    mw_beam = """  beam_init%n_particle = 1024
  beam_init%bunch_charge = 2.401661485427e-14
  beam_init%distribution_type(3) = "RAN_GAUSS"
  beam_init%sig_z = 1.0e-9
  beam_init%sig_pz = 5.282703940115e-03
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
"""
    wake_extra = """  chamber_wake%on = T
  chamber_wake%radius = 2.5e-3
  chamber_wake%conductivity = 5.813e7
  chamber_wake%relaxation = 8.1e-6
  chamber_wake%gap = 0.5e-3
  chamber_wake%lgap = 0.015
  chamber_wake%hrough = 100e-9
  chamber_wake%lrough = 100e-6
"""
    mig_extra = "  global%migrate = T\n  global%migrate_check = T\n"
    KS_C = 2 * np.pi / LAMBDA0 * 2.99792458e8      # phase per second of the file's time

    def mig_events(root):
        """[(s, moved, dropped charge, margin)] per event, and the summary totals."""
        rows, moved, dropped = [], None, None
        for ln in (wd / f"{root}.migration.txt").read_text().splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            if f[0] == "moved":
                moved = int(f[1])
            elif f[0] == "charge_dropped_total":
                dropped = float(f[1])
            elif f[0] == "worst_bunching_deviation":
                pass
            else:
                rows.append((float(f[0]), int(f[1]), float(f[2]), float(f[4])))
        return rows, moved, dropped

    def wake_table(root):
        """[(z, eloss array)] per z-stamped block of <root>.wake.txt."""
        blocks, z, vals = [], None, []
        for ln in (wd / f"{root}.wake.txt").read_text().splitlines():
            m = re.match(r"# z =\s+(\S+)", ln)
            if m:
                if z is not None:
                    blocks.append((z, np.array(vals)))
                z, vals = float(m.group(1)), []
            elif ln.strip() and not ln.startswith("#"):
                vals.append(float(ln.split()[1]))
        if z is not None:
            blocks.append((z, np.array(vals)))
        return blocks

    def identity_and_floor(full, other):
        """Exact identity of ids, populations and membership over the shared frames, the
        field's worst relative difference, and the largest absolute phase difference [rad]
        by id at each shared s."""
        ff, fo = frames_by_s(full), frames_by_s(other)
        common = sorted(set(ff) & set(fo))
        exact = True
        w_field = w_pz = 0.0
        dtheta = {}
        for s in common:
            bf = beam_by_id(str(ff[s]).replace(".wf.h5", ".beam.h5"))
            bo = beam_by_id(str(fo[s]).replace(".wf.h5", ".beam.h5"))
            same = (np.array_equal(bf["id"], bo["id"]) and np.array_equal(bf["pops"], bo["pops"])
                    and np.array_equal(bf["sl"], bo["sl"]))
            exact = exact and same
            uf = fieldio.read_field(ff[s])["u"]
            uo = fieldio.read_field(fo[s])["u"]
            w_field = max(w_field, float(np.max(np.abs(uf - uo))) / max(float(np.max(np.abs(uf))), 1e-300))
            if same:
                dtheta[s] = float(np.max(np.abs(bf["t"] - bo["t"]))) * KS_C
                w_pz = max(w_pz, float(np.max(np.abs(bf["pz"] - bo["pz"]))) / max(float(np.max(np.abs(bf["pz"]))), 1e-300))
        return exact, w_field, w_pz, dtheta, common

    (wd / "mw.bmad").write_text(mw_lat)
    rs_run("mwprep", rs_nml.format(lat="mw.bmad", root="mwprep",
           extra=mig_extra + '  slicing%window_length = 4.8e-9\n  global%track_end = "UND"\n',
           beam=mw_beam, field=gen_field))
    shutil.copy(wd / "mwprep-final.wf.h5", wd / "mwramp.wf.h5")
    with h5py.File(wd / "mwramp.wf.h5", "r+") as h5:
        m = h5[fieldio.MESH_PATH]
        u = m["x"][...]
        m["x"][...] = u * (1.0 + 0.5 * np.arange(u.shape[0]) / u.shape[0])[:, None, None]
        if "slippageResidual" in h5.attrs:
            del h5.attrs["slippageResidual"]
    mw_imp = dict(beam='  beam_file = "mwprep-final.beam.h5"\n', field='  field_file = "mwramp.wf.h5"\n')
    s_ckpt = 0.75

    for tag, wake, what in (("mwm", "", "migration alone"), ("mww", wake_extra, "with the chamber wake")):
        full, toD, rest = f"{tag}full", f"{tag}toD", f"{tag}rest"
        rs_run(full, rs_nml.format(lat="mw.bmad", root=full,
               extra=mig_extra + wake + '  global%dump_beam_at = "D"\n  global%dump_field_at = "D"\n  global%dump_at_comb = T\n', **mw_imp))
        rs_run(toD, rs_nml.format(lat="mw.bmad", root=toD, extra=mig_extra + wake + '  global%track_end = "D"\n', **mw_imp))
        rs_run(rest, rs_nml.format(lat="mw.bmad", root=rest,
               extra=mig_extra + wake + '  global%continuation = T\n  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{full}-at2-D.beam.h5"\n', field=f'  field_file = "{full}-at2-D.wf.h5"\n'))
        resid = residual_of(wd / f"{full}-at2-D.wf.h5")
        ev_full, mv_full, dr_full = mig_events(full)
        ev_toD, mv_toD, dr_toD = mig_events(toD)
        ev_rest, mv_rest, dr_rest = mig_events(rest)
        before = [e for e in ev_full if e[0] <= s_ckpt + 1e-9]
        after = [e for e in ev_full if e[0] > s_ckpt + 1e-9]
        check(f"{what}: the checkpoint follows actual moves and drops and carries a residual",
              len(before) >= 1 and sum(e[1] for e in before) > 0 and sum(e[2] for e in before) > 0
              and resid is not None and 0.5 < abs(resid) < 3.2,
              note=f"[{sum(e[1] for e in before)} moves and {sum(e[2] for e in before):.3e} C dropped before "
                   f"s = {s_ckpt}; slippageResidual {resid}]")

        exact, w_field, w_pz, dtheta, common = identity_and_floor(full, rest)
        # The margin at every later event against the phase error propagated to the frame
        # before it. The event count is the element ends of the continuation, so no
        # decision went unrecorded.
        n_ends_after = 3
        margins_ok = len(ev_rest) == n_ends_after and len(after) == n_ends_after
        worst_ratio = float("inf")
        for s_ev, _, _, mg in ev_rest:
            prev = max((x for x in common if x < s_ev - 1e-9), default=None)
            err = dtheta.get(prev, float("nan"))
            ratio = mg / err if err and err > 0 else float("inf")
            worst_ratio = min(worst_ratio, ratio)
            margins_ok = margins_ok and err == err and mg > 100 * err
        check(f"{what}: every later decision clears its margin by 100x the propagated phase error",
              margins_ok,
              note=f"[{len(ev_rest)} events after the checkpoint, margins "
                   f"{', '.join(f'{e[3]:.1e}' for e in ev_rest)} rad against propagated errors "
                   f"{', '.join(f'{dtheta.get(max((x for x in common if x < e[0] - 1e-9), default=None), float(chr(110)+chr(97)+chr(110))):.1e}' for e in ev_rest)} rad; smallest ratio {worst_ratio:.1e}]")
        check(f"{what}: identity of ids, populations and membership at every frame, field and energy at the floor",
              exact and w_field <= 1e-10 and w_pz <= 1e-10,
              note=f"[{len(common)} frames; identical {exact}; field {w_field:.2e}, energy {w_pz:.2e}, "
                   f"phase {max(dtheta.values()) if dtheta else float(chr(110)+chr(97)+chr(110)):.2e} rad vs 1e-10]")
        events_ok = ([(round(e[0], 9), e[1]) for e in after] == [(round(e[0], 9), e[1]) for e in ev_rest]
                     and all(abs(a[2] - r[2]) <= 1e-10 * max(abs(a[2]), 1e-300) for a, r in zip(after, ev_rest)))
        check(f"{what}: every later migration event matches, position and moves exactly, dropped charge at the floor",
              events_ok, note=f"[{[(round(e[0], 3), e[1]) for e in ev_rest]}]")
        e_full, e_toD, e_rest = escaped_of(full), escaped_of(toD), escaped_of(rest)
        inc = e_full - e_toD
        check(f"{what}: since the checkpoint, moves exact, dropped charge and escaped energy at the floor",
              mv_full - mv_toD == mv_rest
              and abs((dr_full - dr_toD) - dr_rest) <= 1e-10 * max(abs(dr_rest), 1e-300)
              and abs(inc - e_rest) <= 1e-10 * max(abs(inc), 1e-300),
              note=f"[moves {mv_full - mv_toD} vs {mv_rest}; dropped {dr_full - dr_toD:.6e} vs {dr_rest:.6e} C; "
                   f"escaped {inc:.6e} vs {e_rest:.6e} J, differing by {abs(inc - e_rest):.1e}]")

        if wake:
            bf_, br_ = wake_table(full), wake_table(rest)
            in_force = [b for b in bf_ if b[0] <= s_ckpt + 1e-9][-1]
            d0 = float(np.max(np.abs(in_force[1] - br_[0][1]))) / max(float(np.max(np.abs(in_force[1]))), 1e-300)
            later_f = [b for b in bf_ if b[0] > s_ckpt + 1e-9]
            later_r = br_[1:]
            pos_ok = [round(b[0], 9) for b in later_f] == [round(b[0], 9) for b in later_r]
            d_later = max((float(np.max(np.abs(a[1] - r[1]))) / max(float(np.max(np.abs(a[1]))), 1e-300)
                           for a, r in zip(later_f, later_r)), default=0.0)
            check("the restart rebuilds the wake table in force at the checkpoint from the surviving charge",
                  d0 <= 1e-10, note=f"[rel {d0:.2e} vs 1e-10; block at z = {in_force[0]:.3f} m]")
            check("the wake refreshes after every later move or drop: events by count and position, then block by block",
                  pos_ok and len(later_f) == n_ends_after and d_later <= 1e-10,
                  note=f"[{len(later_r)} refreshes at z = {[round(b[0], 3) for b in later_r]}; worst rel {d_later:.2e}]")

    # Boundary sensitivity, on the migration-alone case: shift every restarted phase by a
    # little less and a little more than the smallest margin any later event recorded. The
    # shift goes into the records' time and into the group's phi0 alike, so the checkpoint
    # stays consistent with its records and the shift is the whole of what changes.
    ev_rest, _, _ = mig_events("mwmrest")
    m_min = min(e[3] for e in ev_rest)
    s_min = [e[0] for e in ev_rest if e[3] == m_min][0]
    outcome = {}
    for side, factor in (("below", 0.9), ("above", 1.1)):
        delta = factor * m_min
        root = f"mwsens_{side}"
        shutil.copy(wd / "mwmfull-at2-D.beam.h5", wd / f"{root}.beam.h5")
        with h5py.File(wd / f"{root}.beam.h5", "r+") as h5:
            name = sorted(h5["data"].keys())[0]
            g = h5[f"data/{name}/particles"]
            g = g[sorted(g.keys())[0]]
            g["time"][...] = g["time"][...] - delta / KS_C
            ck = h5["lucifer"]
            ck.attrs.modify("phi0", np.atleast_1d(ck.attrs["phi0"]) + delta)
        rs_run(root, rs_nml.format(lat="mw.bmad", root=root,
               extra=mig_extra + '  global%continuation = T\n  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{root}.beam.h5"\n', field='  field_file = "mwmfull-at2-D.wf.h5"\n'))
        exact, w_field, _, _, common = identity_and_floor("mwmfull", root)
        ev, _, _ = mig_events(root)
        outcome[side] = (exact, w_field, [(round(e[0], 3), e[1]) for e in ev])
    check("boundary sensitivity: a phase shift below the smallest recorded margin flips no decision, one above flips one",
          outcome["below"][0] and not outcome["above"][0],
          note=f"[margin {m_min:.3e} rad at s = {s_min:.3f}; below: identical {outcome['below'][0]}, "
               f"field {outcome['below'][1]:.1e}, events {outcome['below'][2]}; above: identical {outcome['above'][0]}, "
               f"field {outcome['above'][1]:.1e}, events {outcome['above'][2]}. This case establishes no equivalence]")

    # ------------------------------------------------------------------
    print("== a restart across a Bmad element wake ==")

    # An element wake acts once per passage through the whole window: the concatenation
    # adds each particle's slice offset to its longitudinal coordinate at its own entry
    # beta, Bmad kicks that bunch, and the split subtracts the offsets back. All of it is
    # scratch of one passage, so a checkpoint at an element end carries what the next
    # passage needs. What a restart has to get right is the placement: a checkpoint before
    # the wake element must apply the kick, and one after it must not apply it again, the
    # dump already holding it. Both checkpoints are judged against a wake-off history run
    # from the same initial state and never from a wake-on dump, which holds the kick that
    # no continuation can undo, and the wake's own effect is the scale the agreement is
    # measured on: a kick omitted or applied twice moves that effect by order one, where
    # the beam's total energy would hide it.

    ew_wake = (", sr_wake = {amp_scale = 1, scale_with_length = T, "
               "longitudinal = {1e17, 0, 0, 0.25, none}}")

    def ew_lat(wake, lord):
        mk = "\nMK: marker, superimpose, ref = PW" if lord else ""
        return (LAT2.replace("D: pipe, l = 0.30", f"PW: pipe, l = 0.30{wake}{mk}")
                    .replace("SEG: line = (UND, D, UND2)", "SEG: line = (UND, PW, UND2)"))

    def ew_frames(root):
        """{s: (field file, element name)} over a run's frame series."""
        out = {}
        for f in sorted(wd.glob(f"{root}-0*.wf.h5")):
            with h5py.File(f) as h5:
                name = np.atleast_1d(h5.attrs["elementName"])[0]
                out[round(float(np.atleast_1d(h5.attrs["sPosition"])[0]), 9)] = (
                    f, name.decode() if isinstance(name, bytes) else str(name))
        return out

    ew_beam = """  beam_init%n_particle = 1024
  beam_init%bunch_charge = 2.401661485427e-14
  beam_init%distribution_type(3) = "RAN_GAUSS"
  beam_init%sig_z = 1.0e-9
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
"""
    # One starting state for both cases: the run to the first undulator's end is the same
    # line either way, and its field is ramped slice by slice as the residual case does.
    (wd / "ew_n.bmad").write_text(ew_lat("", False))
    rs_run("ewprep", rs_nml.format(lat="ew_n.bmad", root="ewprep",
           extra='  slicing%window_length = 4.8e-9\n  global%track_end = "UND"\n',
           beam=ew_beam, field=gen_field))
    shutil.copy(wd / "ewprep-final.wf.h5", wd / "ewramp.wf.h5")
    with h5py.File(wd / "ewramp.wf.h5", "r+") as h5:
        m = h5[fieldio.MESH_PATH]
        u = m["x"][...]
        m["x"][...] = u * (1.0 + 0.5 * np.arange(u.shape[0]) / u.shape[0])[:, None, None]
        if "slippageResidual" in h5.attrs:
            del h5.attrs["slippageResidual"]
    ew_imp = dict(beam='  beam_file = "ewprep-final.beam.h5"\n', field='  field_file = "ewramp.wf.h5"\n')

    for lord, what in ((False, "on the element itself"), (True, "resolved through a split lord")):
        tag = "lw" if lord else "ew"
        (wd / f"{tag}_w.bmad").write_text(ew_lat(ew_wake, lord))
        (wd / f"{tag}_n.bmad").write_text(ew_lat("", lord))
        dumps = '  global%dump_beam_at = "UND,PW*"\n  global%dump_field_at = "UND,PW*"\n  global%dump_at_comb = T\n'
        rs_run(f"{tag}on", rs_nml.format(lat=f"{tag}_w.bmad", root=f"{tag}on", extra=dumps, **ew_imp))
        rs_run(f"{tag}off", rs_nml.format(lat=f"{tag}_n.bmad", root=f"{tag}off", extra=dumps, **ew_imp))

        # Which passage applies the kick is read from the runs rather than assumed: the
        # first frame whose beam differs from the wake-off history names the element,
        # which for a split lord is the slave holding the lord's midpoint.
        fon, foff = ew_frames(f"{tag}on"), ew_frames(f"{tag}off")
        shared = sorted(set(fon) & set(foff))
        applied_at = None
        for s in shared:
            bo = beam_by_id(str(fon[s][0]).replace(".wf.h5", ".beam.h5"))
            bf = beam_by_id(str(foff[s][0]).replace(".wf.h5", ".beam.h5"))
            d = float(np.max(np.abs(bo["pz"] - bf["pz"]))) / max(float(np.max(np.abs(bf["pz"]))), 1e-300)
            if d > 1e-12:
                applied_at = (s, fon[s][1], d)
                break
        # The continuation opens at the element right after the one the checkpoint completed,
        # which for the split lord is the superimposed marker between its two slaves.
        after = {"PW": "UND2", "PW#1": "MK"}.get(applied_at[1] if applied_at else "", "UND2")
        check(f"the wake {what} acts once, and the passage that applies it is read from the run",
              applied_at is not None and applied_at[2] > 1e-6,
              note=f"[first change at s = {applied_at[0]} m in {applied_at[1]}, "
                   f"relative energy change {applied_at[2]:.2e}; the continuation restarts at {after}]")
        apply_ele = applied_at[1]

        def ew_restart(root, src, start, wake):
            text = rs_nml.format(lat=f"{tag}_{'w' if wake else 'n'}.bmad", root=root,
                                 extra=f'  global%continuation = T\n  global%track_start = "{start}"\n  global%dump_at_comb = T\n',
                                 beam=f'  beam_file = "{src}.beam.h5"\n',
                                 field=f'  field_file = "{src}.wf.h5"\n')
            rs_run(root, text)

        ew_restart(f"{tag}rA", f"{tag}on-at1-UND", apply_ele, True)
        ew_restart(f"{tag}rB", f"{tag}on-at2-{apply_ele}", after, True)
        ew_restart(f"{tag}rOff", f"{tag}off-at1-UND", apply_ele, False)

        for label, root in ((f"before the wake, restarting at {apply_ele}", f"{tag}rA"),
                            (f"after the wake, restarting at {after}", f"{tag}rB")):
            fr = ew_frames(root)
            common = sorted(set(fon) & set(fr))
            exact = True
            w_field = w_pz = w_theta = 0.0
            effect = w_effect = 0.0
            pz_scale = 1.0
            for s in common:
                bo = beam_by_id(str(fon[s][0]).replace(".wf.h5", ".beam.h5"))
                br = beam_by_id(str(fr[s][0]).replace(".wf.h5", ".beam.h5"))
                bf = beam_by_id(str(foff[s][0]).replace(".wf.h5", ".beam.h5"))
                same = (np.array_equal(bo["id"], br["id"]) and np.array_equal(bo["pops"], br["pops"])
                        and np.array_equal(bo["sl"], br["sl"]))
                exact = exact and same
                uo = fieldio.read_field(fon[s][0])["u"]
                ur = fieldio.read_field(fr[s][0])["u"]
                w_field = max(w_field, float(np.max(np.abs(uo - ur))) / max(float(np.max(np.abs(uo))), 1e-300))
                if same:
                    scale = max(float(np.max(np.abs(bo["pz"]))), 1e-300)
                    w_pz = max(w_pz, float(np.max(np.abs(bo["pz"] - br["pz"]))) / scale)
                    w_theta = max(w_theta, float(np.max(np.abs(bo["t"] - br["t"]))) * KS_C)
                    # Both in energy units, and divided only at the end: at a frame
                    # before the kick the effect is zero, and a per-frame ratio there
                    # says nothing.
                    d_cont = bo["pz"] - bf["pz"]
                    d_rest = br["pz"] - bf["pz"]
                    effect = max(effect, float(np.max(np.abs(d_cont))))
                    w_effect = max(w_effect, float(np.max(np.abs(d_cont - d_rest))))
                    pz_scale = scale
            rel_effect = effect / pz_scale
            reproduced = w_effect / max(effect, 1e-300)
            check(f"a checkpoint {label}: the kick is neither omitted nor applied twice",
                  exact and rel_effect > 1e-6 and reproduced <= 1e-6 and w_pz <= 1e-10
                  and w_field <= 1e-10 and w_theta <= 1e-10,
                  note=f"[{len(common)} frames; identical {exact}; the wake moves the energies by "
                       f"{rel_effect:.2e} and the restart reproduces that to {reproduced:.2e} of it; "
                       f"field {w_field:.1e}, energy {w_pz:.1e}, phase {w_theta:.1e} rad vs 1e-10]")

        fr = ew_frames(f"{tag}rOff")
        common = sorted(set(foff) & set(fr))
        w_off = 0.0
        for s in common:
            bf = beam_by_id(str(foff[s][0]).replace(".wf.h5", ".beam.h5"))
            bo = beam_by_id(str(fr[s][0]).replace(".wf.h5", ".beam.h5"))
            w_off = max(w_off, float(np.max(np.abs(bf["pz"] - bo["pz"]))) / max(float(np.max(np.abs(bf["pz"]))), 1e-300))
        check(f"the wake-off history restarted from its own checkpoint reproduces itself ({what})",
              w_off <= 1e-10, note=f"[{len(common)} frames; worst energy {w_off:.2e} vs 1e-10]")



    # ------------------------------------------------------------------
    print("== a continuation across an unaveraged segment boundary ==")

    # An unaveraged segment converts on entry and exit: it asserts the averaged chart,
    # marks that px is now the kinetic momentum, jumps the ramp phase and carries the
    # entry field a half substep to the midpoint the record advances, and undoes all of
    # it at the last step. None of that outlives the segment, so an element-boundary dump
    # holds the field on its plane and the particles in the averaged chart, and a
    # continuation there runs the conversions once. Both boundaries are taken, before and
    # after the segment, and both hold at the floor. The one before it did not until the
    # checkpoint group existed: the standard records fold phi0 into each particle's time,
    # the reader rebuilt an equal phase at another magnitude, off by phi0 times the machine
    # epsilon, and the segment's gain amplified that to 6.5e-8 from an entry state
    # bit-identical in every record (doc/validation.md). The group carries phi0 and each z
    # as the tracker holds them, and the continuation rebuilds the split to the bit before
    # any conversion or advance, which is checked first, so a later disagreement could only
    # be amplification. A frame from inside the segment, from inside an averaged undulator
    # or from a sliced endpoint inside one carries no group and is refused, as is every
    # incomplete or inconsistent group and a continuation point that does not follow the
    # checkpoint. The same checkpoint loaded without the switch is an initialization: the
    # records are folded and the residual starts at zero.

    uv_lat = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 8.53711
beginning[alpha_a] = -0.703306
beginning[beta_b] = 17.3899
beginning[alpha_b] = 1.40348
UNDA: wiggler, l = 0.6, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      fel_method = averaged, ds_step = 0.015
UNDB: UNDA, fel_method = unaveraged
P1: pipe, l = 0.21
SEG: line = (UNDA, P1, UNDB, P1, UNDA)
use, SEG
"""
    (wd / "uvsand.bmad").write_text(uv_lat)
    # The same line with the first undulator cut in two by a superimposed marker, so a dump
    # at the cut is a sliced endpoint inside a physical undulator.
    (wd / "uvsplit.bmad").write_text(uv_lat.replace("UNDB: UNDA, fel_method = unaveraged\n",
                                                    "UNDB: UNDA, fel_method = unaveraged\nUNDS: UNDA\nMK: marker, superimpose, ref = UNDS\n")
                                           .replace("SEG: line = (UNDA, P1, UNDB, P1, UNDA)", "SEG: line = (UNDS, P1, UNDB, P1, UNDA)"))
    uv_beam = """  beam_init%n_particle = 1024
  beam_init%bunch_charge = 2.401661485427e-14
  beam_init%distribution_type(3) = "RAN_GAUSS"
  beam_init%sig_z = 1.0e-9
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
"""
    CONT = '  global%continuation = T\n'

    def uv_frames(root):
        out = {}
        for f in sorted(wd.glob(f"{root}-0*.wf.h5")):
            with h5py.File(f) as h5:
                out[round(float(np.atleast_1d(h5.attrs["sPosition"])[0]), 9)] = f
        return out

    def raw_field(path):
        return fieldio.read_field(path)["u"]

    def raw_beam(path):
        with h5py.File(path) as h5:
            n0 = sorted(h5["data"].keys())[0]
            g = h5[f"data/{n0}/particles"]
            g = g[sorted(g.keys())[0]]
            n = int(np.atleast_1d(g.attrs["numParticles"])[0])
            rec = lambda k: g[k][...] if g[k].shape else np.full(n, g[k][()])
            ids = rec("id")
            o = np.argsort(ids)
            return {"time": rec("time")[o], "id": ids[o], "pz": rec("totalMomentum")[o],
                    "x": rec("position/x")[o], "y": rec("position/y")[o]}

    def group_of(path):
        """The checkpoint group's members, or None where the frame carries none."""
        with h5py.File(path) as h5:
            g = h5.get("lucifer")
            if g is None:
                return None
            return {"phi0": float(np.atleast_1d(g.attrs["phi0"])[0]),
                    "s": float(np.atleast_1d(g.attrs["sPosition"])[0]),
                    "z": g["z"][...], "id": g["id"][...], "counts": g["sliceCount"][...]}

    def floor_of(full, other):
        """Shared frames by position: identity of ids, populations and membership, and the
        worst field, energy and phase differences by id."""
        fon, fr = uv_frames(full), uv_frames(other)
        common = sorted(set(fon) & set(fr))
        exact = True
        w_field = w_pz = w_theta = 0.0
        for sc in common:
            bo = beam_by_id(str(fon[sc]).replace(".wf.h5", ".beam.h5"))
            br = beam_by_id(str(fr[sc]).replace(".wf.h5", ".beam.h5"))
            same = (np.array_equal(bo["id"], br["id"]) and np.array_equal(bo["pops"], br["pops"])
                    and np.array_equal(bo["sl"], br["sl"]))
            exact = exact and same
            w_field = max(w_field, float(np.max(np.abs(raw_field(fon[sc]) - raw_field(fr[sc])))) / max(float(np.max(np.abs(raw_field(fon[sc])))), 1e-300))
            if same:
                w_pz = max(w_pz, float(np.max(np.abs(bo["pz"] - br["pz"]))) / max(float(np.max(np.abs(bo["pz"]))), 1e-300))
                w_theta = max(w_theta, float(np.max(np.abs(bo["t"] - br["t"]))) * KS_C)
        return len(common), exact, w_field, w_pz, w_theta

    # The prepared inputs, an initialization as before: the beam from one run's final
    # dump and its field ramped slice by slice with the residual stripped, neither a
    # checkpoint pair.
    rs_run("uvprep", rs_nml.format(lat="uvsand.bmad", root="uvprep",
           extra='  slicing%window_length = 4.8e-9\n  global%track_end = "UNDA##1"\n',
           beam=uv_beam, field=gen_field))
    shutil.copy(wd / "uvprep-final.wf.h5", wd / "uvramp.wf.h5")
    with h5py.File(wd / "uvramp.wf.h5", "r+") as h5:
        m = h5[fieldio.MESH_PATH]
        u = m["x"][...]
        m["x"][...] = u * (1.0 + 0.5 * np.arange(u.shape[0]) / u.shape[0])[:, None, None]
        if "slippageResidual" in h5.attrs:
            del h5.attrs["slippageResidual"]
    uv_imp = dict(beam='  beam_file = "uvprep-final.beam.h5"\n', field='  field_file = "uvramp.wf.h5"\n')

    rs_run("uvfull", rs_nml.format(lat="uvsand.bmad", root="uvfull",
           extra='  global%dump_beam_at = "P1,UNDB"\n  global%dump_field_at = "P1,UNDB"\n  global%dump_at_comb = T\n', **uv_imp))
    r_before = residual_of(wd / "uvfull-at2-P1.wf.h5")
    r_after = residual_of(wd / "uvfull-at3-UNDB.wf.h5")
    check("both unaveraged-boundary checkpoints carry a residual off a whole-slice boundary, and the "
          "prepared inputs initialized as before",
          r_before is not None and 0.5 < abs(r_before) < 3.2 and r_after is not None and 0.5 < abs(r_after) < 3.2,
          note=f"[before {r_before}, after {r_after}, of 4 a slice]")

    fon = uv_frames("uvfull")
    ckb = group_of(wd / "uvfull-at2-P1.beam.h5")
    cka = group_of(wd / "uvfull-at3-UNDB.beam.h5")
    interior_b = max((sc for sc in fon if 0.81 + 1e-9 < sc < 1.41 - 1e-9), default=None)
    interior_a = max((sc for sc in fon if 1e-9 < sc < 0.6 - 1e-9), default=None)
    n_grp = sum(1 for sc in fon if group_of(str(fon[sc]).replace(".wf.h5", ".beam.h5")) is not None)
    check("the checkpoint group rides the boundary dumps and the element-end frames, and no interior frame",
          ckb is not None and cka is not None and n_grp == 6
          and group_of(str(fon[interior_b]).replace(".wf.h5", ".beam.h5")) is None
          and group_of(str(fon[interior_a]).replace(".wf.h5", ".beam.h5")) is None,
          note=f"[{n_grp} of {len(fon)} frames carry it: the entry face and five element ends; "
               f"none at s = {interior_a} (averaged) or s = {interior_b} (unaveraged)]")

    rs_run("uvrafter", rs_nml.format(lat="uvsand.bmad", root="uvrafter",
           extra=CONT + '  global%track_start = "P1##2"\n  global%dump_at_comb = T\n',
           beam='  beam_file = "uvfull-at3-UNDB.beam.h5"\n', field='  field_file = "uvfull-at3-UNDB.wf.h5"\n'))
    n_a, exact_a, wf_a, wp_a, wt_a = floor_of("uvfull", "uvrafter")
    check("continuation after the segment: field and particles at the floor through the averaged segment",
          exact_a and wf_a <= 1e-10 and wp_a <= 1e-10 and wt_a <= 1e-10,
          note=f"[{n_a} frames; identical {exact_a}; field {wf_a:.1e}, energy {wp_a:.1e}, phase {wt_a:.1e} rad vs 1e-10]")

    # Before the segment. The internal state first: the continuation's first frame is
    # written before any entrance conversion or advance, and its group and records are
    # the checkpoint's to the bit, so whatever follows is the segment's doing.
    rs_run("uvrbefore", rs_nml.format(lat="uvsand.bmad", root="uvrbefore",
           extra=CONT + '  global%track_start = "UNDB"\n  global%dump_at_comb = T\n',
           beam='  beam_file = "uvfull-at2-P1.beam.h5"\n', field='  field_file = "uvfull-at2-P1.wf.h5"\n'))
    frb = uv_frames("uvrbefore")
    entry = round(0.81, 9)
    g0 = group_of(str(frb[entry]).replace(".wf.h5", ".beam.h5"))
    same_state = (g0 is not None and g0["phi0"] == ckb["phi0"] and np.array_equal(g0["z"], ckb["z"])
                  and np.array_equal(g0["id"], ckb["id"]) and np.array_equal(g0["counts"], ckb["counts"]))
    same_resid = residual_of(frb[entry]) == r_before
    be = raw_beam(str(wd / "uvfull-at2-P1.beam.h5"))
    bc = raw_beam(str(frb[entry]).replace(".wf.h5", ".beam.h5"))
    same_records = all(np.array_equal(be[k], bc[k]) for k in ("time", "id", "pz", "x", "y"))
    same_field = float(np.max(np.abs(raw_field(wd / "uvfull-at2-P1.wf.h5") - raw_field(frb[entry])))) == 0.0
    check("continuation before the segment: phi0, every z, the residual, the records and the field are the "
          "checkpoint's bit for bit before any conversion or advance",
          same_state and same_resid and same_records and same_field,
          note=f"[phi0 {ckb['phi0']:.9f} rad over {len(ckb['z'])} particles; group {same_state}, residual {same_resid}, "
               f"records {same_records}, field {same_field}]")
    n_b, exact_b, wf_b, wp_b, wt_b = floor_of("uvfull", "uvrbefore")
    check("continuation before the segment: field and particles at the floor through the unaveraged segment "
          "and the averaged one after it",
          exact_b and wf_b <= 1e-10 and wp_b <= 1e-10 and wt_b <= 1e-10,
          note=f"[{n_b} frames; identical {exact_b}; field {wf_b:.1e}, energy {wp_b:.1e}, phase {wt_b:.1e} rad vs 1e-10; "
               f"the folded records alone reached 6.5e-8 here]")

    # The same checkpoint loaded the other way, an initialization, stopped after its
    # initial dump: the records are folded, phi0 and the residual start at zero.
    rs_run("uvinit", rs_nml.format(lat="uvsand.bmad", root="uvinit",
           extra='  global%track_start = "UNDB"\n  global%load_only = T\n  global%write_initial = T\n',
           beam='  beam_file = "uvfull-at2-P1.beam.h5"\n', field='  field_file = "uvfull-at2-P1.wf.h5"\n'))
    gi = group_of(wd / "uvinit-initial.beam.h5")
    ri = residual_of(wd / "uvinit-initial.wf.h5")
    dz_init = float(np.max(np.abs(gi["z"] - ckb["z"]))) if gi is not None else -1.0
    check("one checkpoint loaded both ways: the initialization folds the records with phi0 and the residual at "
          "zero, the continuation restores both exactly, so the declaration controls the loading",
          gi is not None and gi["phi0"] == 0.0 and dz_init > 0 and ri == 0.0 and same_state and same_resid,
          note=f"[initialization: phi0 0, every z moved by up to {dz_init:.2e} m, residual 0; "
               f"continuation: phi0 {ckb['phi0']:.6f} rad, residual {r_before} wavelengths]")

    # Refusals, each naming its condition.
    def refused(root, beam, field, start, fragment, lat="uvsand.bmad"):
        return rs_run(root, rs_nml.format(lat=lat, root=root,
                      extra=CONT + f'  global%track_start = "{start}"\n',
                      beam=f'  beam_file = "{beam}"\n', field=f'  field_file = "{field}"\n'),
                      expect_fail=True, fragment=fragment)

    def variant(name, src, edit):
        shutil.copy(wd / src, wd / name)
        with h5py.File(wd / name, "r+") as h5:
            edit(h5)
        return name

    fb = str(fon[interior_b]).replace(".wf.h5", ".beam.h5")
    check("a continuation from a frame inside the unaveraged segment is refused, no group",
          refused("uvr_inu", pathlib.Path(fb).name, fon[interior_b].name, "P1##2", "NO CHECKPOINT GROUP"))
    fa = str(fon[interior_a]).replace(".wf.h5", ".beam.h5")
    check("a continuation from a frame inside an averaged undulator is refused, no group",
          refused("uvr_ina", pathlib.Path(fa).name, fon[interior_a].name, "P1##1", "NO CHECKPOINT GROUP"))

    # A superimposed element cannot cut an undulator this tracker tracks: the slaves refer
    # to their lord for the field model, and the FEL element assertion refuses the lattice
    # at setup. So no frame from a sliced endpoint inside an undulator exists to continue
    # from, and fel_checkpoint_eligible's guard against that boundary waits for the day one
    # does.
    split_refused = rs_run("uvsplit", rs_nml.format(lat="uvsplit.bmad", root="uvsplit",
           extra='  global%track_end = "UNDS#2"\n  global%dump_beam_at = "UNDS#1"\n  global%dump_field_at = "UNDS#1"\n', **uv_imp),
           expect_fail=True, fragment="SUPER_SLAVE")
    check("an undulator cut by a superimposed marker is refused at setup, so no sliced endpoint inside one is "
          "ever dumped", split_refused,
          note="[UNDS cut at its middle: its slaves refer to the lord for the field model, and the message names the slave]")

    def del_z(h5): del h5["lucifer/z"]
    def del_group(h5): del h5["lucifer"]
    def del_resid(h5): del h5.attrs["slippageResidual"]
    def permute(h5): h5["lucifer/z"][...] = h5["lucifer/z"][...][::-1]
    def bad_version(h5): h5["lucifer"].attrs["checkpointFormat"] = np.bytes_("lucifer-checkpoint 9.9")
    check("a checkpoint with a member removed is refused",
          refused("uvr_member", variant("uvck_noz.beam.h5", "uvfull-at2-P1.beam.h5", del_z), "uvfull-at2-P1.wf.h5",
                  "UNDB", "CHECKPOINT MEMBER MISSING: z"))
    check("a checkpoint with the group removed is refused",
          refused("uvr_group", variant("uvck_nogrp.beam.h5", "uvfull-at2-P1.beam.h5", del_group), "uvfull-at2-P1.wf.h5",
                  "UNDB", "NO CHECKPOINT GROUP"))
    check("a checkpoint whose field file lacks the residual is refused",
          refused("uvr_resid", "uvfull-at2-P1.beam.h5", variant("uvck_nores.wf.h5", "uvfull-at2-P1.wf.h5", del_resid),
                  "UNDB", "NO slippageResidual"))
    check("a beam frame and a field frame from two positions are refused",
          refused("uvr_pos", "uvfull-at2-P1.beam.h5", "uvfull-at3-UNDB.wf.h5", "UNDB", "CHECKPOINT POSITIONS DISAGREE"))
    check("a checkpoint whose coordinates are permuted against its records is refused",
          refused("uvr_perm", variant("uvck_perm.beam.h5", "uvfull-at2-P1.beam.h5", permute), "uvfull-at2-P1.wf.h5",
                  "UNDB", "CHECKPOINT COORDINATES DISAGREE"))
    check("a checkpoint of a version this tracker does not know is refused",
          refused("uvr_ver", variant("uvck_ver.beam.h5", "uvfull-at2-P1.beam.h5", bad_version), "uvfull-at2-P1.wf.h5",
                  "UNDB", "CHECKPOINT VERSION NOT KNOWN"))
    check("a continuation point one element past the checkpoint's is refused",
          refused("uvr_start", "uvfull-at2-P1.beam.h5", "uvfull-at2-P1.wf.h5", "P1##2", "CONTINUATION POINT DOES NOT FOLLOW"))

    # ------------------------------------------------------------------
    print("== committed prose cites committed artifacts ==")

    # The numbered findings live in the design repository and do not ship with this
    # program, so a reference to one dangles for every reader of this tree. A comment
    # states the fact, or points at the page here that records it. The generated pages
    # follow the headers they are built from, so they are not walked.
    lucifer = pathlib.Path(__file__).resolve().parents[2]
    skip = {"generated", "_build", "__pycache__", ".pytest_cache"}
    # Spelled in halves so this file is not its own offender. Naming the file as an
    # exception instead would blind the check to a real citation written here.
    token = "FIND" + "INGS"
    offenders = []
    for f in sorted(lucifer.rglob("*")):
        if not f.is_file() or f.suffix not in SOURCE_SUFFIXES:
            continue
        if skip & set(f.relative_to(lucifer).parts):
            continue
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        for i, line in enumerate(text.split("\n"), 1):
            if token in line:
                offenders.append(f"{f.relative_to(lucifer.parent)}:{i}")
    check("no committed file cites the design repository's findings", not offenders,
          note=f"[{len(offenders)} reference(s)"
               + (": " + ", ".join(offenders[:6]) if offenders else "") + "]")

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
