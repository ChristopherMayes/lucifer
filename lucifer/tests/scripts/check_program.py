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
  8. A restart across a slippage residual: the window check above is steady state,
     where slippage is a no-op, so a time-dependent pair follows it. The continuous run
     dumps at a pipe whose accumulated slip leaves half a slice of residual, the dump's
     slippageResidual attribute carries it, and the restart must reproduce the
     continuous run's frames, field slice by slice and particles by id, to the same
     floor. A restart from the dump with the attribute stripped rotates on its own
     schedule and differs from the second step of the second undulator, a pipe leaving a
     residual near zero makes that difference vanish, and an impossible residual is
     refused.

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
    # leaves half a slice of residual, the other leaves a residual near zero, and a
    # restart from a dump stripped of the attribute (a legacy file) tells the two apart,
    # which is the demonstration that the residual is the cause.

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
        full, rest, legacy, toD = f"{tag}full", f"{tag}rest", f"{tag}leg", f"{tag}toD"
        rs_run(full, rs_nml.format(lat=lat, root=full,
               extra='  global%dump_beam_at = "D"\n  global%dump_field_at = "D"\n  global%dump_at_comb = T\n', **imp))
        rs_run(toD, rs_nml.format(lat=lat, root=toD, extra='  global%track_end = "D"\n', **imp))
        dump_b, dump_f = f"{full}-at2-D.beam.h5", f"{full}-at2-D.wf.h5"
        resid = residual_of(wd / dump_f)
        # The restart carrying the residual, and the legacy restart with it stripped.
        rs_run(rest, rs_nml.format(lat=lat, root=rest,
               extra='  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{dump_b}"\n', field=f'  field_file = "{dump_f}"\n'))
        shutil.copy(wd / dump_f, wd / f"{tag}legacy.wf.h5")
        with h5py.File(wd / f"{tag}legacy.wf.h5", "r+") as h5:
            del h5.attrs["slippageResidual"]
        rs_run(legacy, rs_nml.format(lat=lat, root=legacy,
               extra='  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
               beam=f'  beam_file = "{dump_b}"\n', field=f'  field_file = "{tag}legacy.wf.h5"\n'))

        w_rest, n_rest, ids_r, pops_r, first_r = compare_runs(full, rest)
        w_leg, n_leg, ids_l, pops_l, first_l = compare_runs(full, legacy)
        worst_rest = max(w_rest.values())
        if kind == "half a slice":
            check(f"the dump at the pipe carries a residual of half a slice ({pipe} m pipe)",
                  resid is not None and 0.5 < abs(resid) < 3.2,
                  note=f"[slippageResidual {resid} wavelengths of 4 a slice]")
        else:
            check(f"the dump at the pipe carries a residual near zero ({pipe} m pipe)",
                  resid is not None and abs(resid) < 1e-3,
                  note=f"[slippageResidual {resid} wavelengths, pipe slip is whole]")
        check(f"restart carrying the residual == continuous, {kind}: field, phase, energy, weight, "
              f"power at {n_rest} frames", worst_rest <= 1e-10 and ids_r and pops_r,
              note=f"[worst rel {worst_rest:.2e} vs 1e-10; ids and populations identical: "
                   f"{ids_r and pops_r}; field {w_rest['field']:.1e} theta {w_rest['theta']:.1e} "
                   f"pz {w_rest['pz']:.1e} power {w_rest['power']:.1e}]")
        if kind == "half a slice":
            check("a legacy restart, the residual stripped, rotates on its own schedule: the field "
                  "differs from the second step of the second undulator",
                  first_l is not None and abs(first_l - (0.45 + pipe + 0.030)) < 1e-6
                  and w_leg["field"] > 1e-2,
                  note=f"[first divergence at s = {first_l} m, expected {0.45 + pipe + 0.030:.3f}; "
                       f"worst field {w_leg['field']:.2e}, theta {w_leg['theta']:.1e}, pz {w_leg['pz']:.1e}]")
        else:
            check("a legacy restart at a near-zero residual agrees, which is why the steady-state "
                  "and whole-slice checks never saw the loss",
                  max(w_leg.values()) <= 1e-10 and ids_l and pops_l,
                  note=f"[worst rel {max(w_leg.values()):.2e} vs 1e-10]")
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
               extra=mig_extra + wake + '  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
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
    # little less and a little more than the smallest margin any later event recorded.
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
        rs_run(root, rs_nml.format(lat="mw.bmad", root=root,
               extra=mig_extra + '  global%track_start = "UND2"\n  global%dump_at_comb = T\n',
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
