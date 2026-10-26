#!/usr/bin/env python3
"""
Checks of the diagnostic output (fel-physics.md, sec-stats): the stats file, the dump_at
elements and the escaped-field bank. Each is held by an identity between two routes
rather than by a reference file.

  0. The stats file against the bmad-stats standard's own validator, and the file's
     description of itself: units, axes, the record axis, provenance, slice
     coordinates, envelope extremes and the element join.
  1. bunch_params reconstruction. bunch_params_from_stats.py applied to the per-record
     sufficient statistics must reproduce the calc_bunch_params values the tracker
     stored at element ends. This shows the per-record datasets are sufficient.
  2. The banked energy equals the ledger's U_escaped column. The slices are the same and
     the bookkeeping paths are independent: wavefront_params at drain time against the
     zero-fill energy sum.
  3. Analytic against numerical propagation. Each banked slice's rms size at the exit,
     from the free-space ABCD map on its moment matrix at bank time, must match the rms
     of the FFT-propagated slice in the pulse file.
  4. Pulse pooling. The pooled whole-pulse sigma from per-slice params must match the
     moments computed directly from the concatenated pulse file.
  5. Thread invariance. Every dataset of the stats, escaped and pulse files is identical
     at 1 and 8 threads. The comparison is by dataset, since HDF5 headers embed
     creation times. Also zero-length wake elements of both polarities.
  6. Refusal. A dump_beam_at entry that matches no element is refused.
  7. The pre-run. The header's work count, cost estimate and load standing, the footer's
     thin-slice list on an input file loaded under the measured floor, load_only stopping
     after the header, and the refusal of a harmonic at which the load carries no shot
     noise.
  8. The frame series (run.dump_at_comb). One frame per stats record, a frame at an
     element end identical to that element's own dump, labels unique in every frame and
     following migration, the frame's attributes against the lattice, frames landing
     inside an element, and a run with frames identical to the same run without.
  9. A frame against the stats row it sits on. Per-slice bunching and transverse
     centroid from the frame's own particles, per-slice power from the field frame, and
     the reduced projections integrating to the row's power. Also a slice range
     cutting both files.
 10. The interlude's own steps (run.interlude_ds_step). The knob at zero reproduces
     the run that never named it byte for byte. Cutting an element reproduces tracking
     it whole. Rows and frames land inside a quadrupole and a pipe at the s the rows
     carry. The stats file is exact-sized against the rows written, and the header
     states the pieces before it tracks. A wake-carrying element is not cut, space
     charge on a seam interlude is refused, a zero-length kicker kicks, and the
     transcribed Genesis interlude refuses the knob.
 11. refposition in the two Genesis-format files means what Genesis means: the in-bunch
     coordinate of slice 1 in Genesis's order. That is 0 for the pulse file and the
     window length for the escaped file, unmoved by a drift along the lattice. The
     escaped file's slices are numbered tail to head. openPMD-beamphysics reads the
     pulse at zmin 0 with the live slices equal to the final field. Genesis itself
     imports the pulse and re-dumps it slice by slice with the same refposition. Genesis
     is an external reader, and a round trip through this program's own reader cannot
     stand in for it.

Run by the benchmark. Exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

from checkrec import record
from lucifer.flat import to_toml
to_groups = to_toml

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bunch_params_from_stats import bunch_params_at, pool_wavefront, M_ELECTRON  # noqa: E402
from lucifer.stats import read_stats, same_data  # noqa: E402
import lucifer.beam as beamio
from lucifer.beam import read_slices  # noqa: E402
import fieldio  # noqa: E402

import validate_bmad_stats  # noqa: E402  The standard's own conformance checker.

FAILED = False

NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 512
  beam.bunch_charge = 8.0e-15
  beam.z_distribution = "GRID"
  beam.z_min = -4e-10
  beam.z_max = 4e-10
  beam.sig_pz = 8.8045e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 8e-10
  n_wavelength = 1
  shot_noise = T
  ran_seed = 999
  keep_escaped_field = T
  dump_beam_at = "UND"
  dump_field_at = "UND"
{extra}&end
"""

WRAP = """call, file = aramis_1seg.bmad
wiggler::*[FEL_METHOD] = unaveraged
"""

# The frame-series input file. A short segment, a coarse comb so the frame count stays small,
# and migration on so labels have to follow particles between slices.

FRAME_NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam.n_particle = 256
  beam.bunch_charge = 8.0e-15
  beam.z_distribution = "GRID"
  beam.z_min = -4e-10
  beam.z_max = 4e-10
  beam.sig_pz = 8.8045e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 8e-10
  n_wavelength = 1
  shot_noise = T
  ran_seed = 999
{extra}&end
"""

# The pre-run input file. It is averaged, seeded so it reaches saturation in one segment at
# little cost, and loaded at 512 macroparticles in beamlets of 8. That is 64 beamlets,
# half the floor at which the source filter's convergence was measured, so the load
# standing and the thin-slice list both have something to report.

PRE_NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam.n_particle = 512
  beam.bunch_charge = 8.0e-15
  beam.z_distribution = "GRID"
  beam.z_min = -4e-10
  beam.z_max = 4e-10
  beam.sig_pz = 8.8045e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 8e-10
  n_wavelength = 1
  shot_noise = T
  ran_seed = 999
{extra}&end
"""

# Zero-length wake elements of both polarities, for check 5a. A zero-length element
# carrying an sr wake is a standard Bmad idiom: a wake assigned over an element range or
# a class lands on zero-length members. What differs is whether the kick can act. With
# scale_with_length = F the kick acts, and the element must be tracked, take a record,
# and repeat coords/s at the plane where the element before it ended. With T the kick is
# identically zero (wake_mod scales by l$), so the element is skipped and the file must
# come out as though the wake were absent. WKT_NONE replaces the T pipes with plain ones
# for that comparison.
ZL_WAKE = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
{wkt}
WKF: pipe, l = 0, sr_wake = {{amp_scale = 1, scale_with_length = F,
  longitudinal = {{1e12, 0, 0, 0.25, none}}}}
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      fel_method = averaged, ds_step = 0.015
SEG: line = (WKT, UND, WKF, WKT)
use, SEG
"""

WKT_WAKE = """WKT: pipe, l = 0, sr_wake = {amp_scale = 1, scale_with_length = T,
  longitudinal = {1e12, 0, 0, 0.25, none}}"""

WKT_NONE = """WKT: pipe, l = 0"""


def particle_frames(wd, root):
    """The particle frame of every record of a series, in record order."""
    files = list(wd.glob(f"{root}-[0-9]*.beam.h5"))
    return sorted(files, key=lambda p: int(p.name.split("-")[-1].split(".")[0]))


def h5_identical(fa, fb):
    """Whether two HDF5 files hold the same data, with meta/ and params/ excluded.

    meta/ holds provenance as datasets (fel-physics.md, sec-meta), since HDF5 caps an
    attribute at 64 kB. input_echo carries out_root and timestamp carries the clock, so
    any two runs differ in meta/. params/ records every key of the input file, the frame and
    dump switches among them, so two runs that differ only in what they write differ
    there and nowhere else. Nothing in either is physics. Only datasets are compared, so
    attributes do not enter.

    same_data counts NaN as equal to NaN, which the stats file needs. See lucifer.stats.
    """
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names_a, names_b = [], []
        keep = lambda n, o: isinstance(o, h5py.Dataset) and not n.startswith(("meta/", "params/"))
        a.visititems(lambda n, o: names_a.append(n) if keep(n, o) else None)
        b.visititems(lambda n, o: names_b.append(n) if keep(n, o) else None)
        if sorted(names_a) != sorted(names_b):
            return False
        return all(same_data(a[n][()], b[n][()]) for n in names_a)


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    record(name, value, tol, "max", ok)
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".toml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-1000:]}")
        sys.exit(1)
    return r



# The interlude's own steps (run.interlude_ds_step). Tao's comb has two halves.
# comb_ds_save is the first: save_a_bunch_step selects among positions a tracker already
# reaches. The interlude step is the second, tao_lattice_calc_mod's n_slice branch,
# which cuts an element that would otherwise be crossed in one map. The lattice is a
# wiggler, a quadrupole whose fringes are live, a pipe and a wiggler. A cut element of
# each kind is exercised, and the wigglers keep their own steps either way.

INT_LAT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 8.53711
beginning[alpha_a] = -0.703306
beginning[beta_b] = 17.3899
beginning[alpha_b] = 1.40348
UND: wiggler, l = 0.3, l_period = 0.015, field_calc = helical_model,
     b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light,
     fel_method = averaged, ds_step = 0.015
QF: quadrupole, l = 0.20, k1 = 8.0, fringe_type = full
P1: pipe, l = 0.40
SEG: line = (UND, QF, P1, UND)
use, SEG
"""

INT_NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 256
  beam.bunch_charge = 1.000692285594e-15
  beam.sig_z = 0
  beam.sig_pz = 8.804506566858e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  ran_seed = 999
  comb_ds_save = 0.02
  write_diag = T
{extra}&end
"""


def interludes(exe, wd):
    """
    Stats rows and frames along an interlude, and what the knob guarantees.

    The knob is off by default, and off must be the path of an input file that never names it.
    The first check is byte identity against such an input file. The second is that cutting an
    element reproduces tracking it whole. The pieces come from Bmad's
    element_slice_iterator, which keeps an element's fringes at its own ends rather than
    giving every piece a pair. On FEL beam sizes the fringe difference is below notice.
    The hand-split lattice, whose every piece carries a fringe pair, checks that here
    at 1e-9.
    """
    print("--- the interlude's own steps:")
    (wd / "intlat.bmad").write_text(INT_LAT)
    hand = INT_LAT.replace("QF: quadrupole, l = 0.20, k1 = 8.0, fringe_type = full",
                           "QFS: quadrupole, l = 0.02, k1 = 8.0, fringe_type = full\n"
                           "QF: line = (10*QFS)")
    hand = hand.replace("P1: pipe, l = 0.40", "P1S: pipe, l = 0.04\nP1: line = (10*P1S)")
    (wd / "intlath.bmad").write_text(hand)

    run(exe, wd, "int_w", INT_NML.format(lat="intlat.bmad", root="int_w", extra=""))
    run(exe, wd, "int_off", INT_NML.format(lat="intlat.bmad", root="int_off",
        extra="  run.interlude_ds_step = 0\n"))
    same = (wd / "int_w.diag.txt").read_bytes() == (wd / "int_off.diag.txt").read_bytes()
    check("interlude: the knob at zero is the run that never named it (0 = yes)",
          0.0 if same else 1.0, 0.5)

    run(exe, wd, "int_p", INT_NML.format(lat="intlat.bmad", root="int_p",
        extra="  run.interlude_ds_step = 0.04\n"))
    run(exe, wd, "int_h", INT_NML.format(lat="intlath.bmad", root="int_h", extra=""))
    w = np.loadtxt(wd / "int_w.diag.txt")
    p = np.loadtxt(wd / "int_p.diag.txt")
    h = np.loadtxt(wd / "int_h.diag.txt")

    # The pieces are tracked, so they must not change the answer. The element's exit is
    # the same however many pieces cross it.
    check("interlude: cut into pieces, exit power against the whole element",
          abs(p[-1, 2] - w[-1, 2]) / w[-1, 2], 1e-9)
    check("interlude: cut into pieces, exit sigma_x against the whole element",
          abs(p[-1, 8] - w[-1, 8]) / w[-1, 8], 1e-9)
    check("interlude: the hand-split lattice agrees too, so the fringe difference the "
          "iterator buys is below notice here",
          abs(h[-1, 2] - w[-1, 2]) / w[-1, 2], 1e-9)

    # With the pieces on, rows land inside the quadrupole and inside the pipe. Tracked
    # whole, neither has a row inside. The element spans are 0.3 to 0.5 m and 0.5 to
    # 0.9 m on this line.
    def inside(rows, z0, z1):
        z = rows[:, 0]
        return int(((z > z0 + 1e-9) & (z < z1 - 1e-9)).sum())
    check("interlude: rows strictly inside the quadrupole, whole (0 expected)",
          abs(inside(w, 0.3, 0.5) - 0), 0.5)
    check("interlude: rows strictly inside the quadrupole, cut (4 expected)",
          abs(inside(p, 0.3, 0.5) - 4), 0.5, note=f"[{inside(p, 0.3, 0.5)}]")
    check("interlude: rows strictly inside the pipe, cut (9 expected)",
          abs(inside(p, 0.5, 0.9) - 9), 0.5, note=f"[{inside(p, 0.5, 0.9)}]")

    # The record-count precompute replays the walk, so the stats file is exact-sized. A
    # mismatch would mean the two disagree about where a row falls.
    with h5py.File(wd / "int_p.stats.h5") as f:
        nrec = f["coords/s"].shape[0]
        step = float(np.ravel(f["params/run/interlude_ds_step"])[0])
    check("interlude: the stats file holds exactly the rows the diag carries",
          abs(nrec - len(p)), 0.5, note=f"[{nrec} against {len(p)}]")
    check("interlude: the stats file records the step it ran at", abs(step - 0.04), 1e-12)

    # Frames along the interlude, which is what the knob is for.
    run(exe, wd, "int_f", INT_NML.format(lat="intlat.bmad", root="int_f",
        extra="  run.interlude_ds_step = 0.04\n  run.dump_at_comb = T\n"))
    frames = particle_frames(wd, "int_f")
    check("interlude: one frame per stats row with the pieces on",
          abs(len(frames) - len(p)), 0.5, note=f"[{len(frames)} frames, {len(p)} rows]")
    zs = []
    for fr in frames:
        with h5py.File(fr) as f:
            zs.append(float(np.ravel(f["lucifer/frame"].attrs["sPosition"])[0]))
    zs = np.array(zs)
    check("interlude: the frames carry the rows' own s", np.abs(zs - p[:, 0]).max(), 1e-9)
    nin = int(((zs > 0.5 + 1e-9) & (zs < 0.9 - 1e-9)).sum())
    check("interlude: frames strictly inside the pipe (9 expected)", abs(nin - 9), 0.5,
          note=f"[{nin}]")

    # The pre-run header states the pieces before it tracks, and says nothing when the
    # element is tracked whole.
    r = run(exe, wd, "int_hdr", INT_NML.format(lat="intlat.bmad", root="int_hdr",
        extra="  run.interlude_ds_step = 0.04\n  run.load_only = T\n"))
    m = re.search(r"Interlude   steps of .*: (\d+) pieces, (\d+) element", r.stdout)
    check("interlude: the header states the pieces before it tracks",
          0.0 if (m and int(m.group(1)) == 15 and int(m.group(2)) == 2) else 1.0, 0.5,
          note=f"[{m.group(0) if m else 'absent'}]")
    r = run(exe, wd, "int_hdr0", INT_NML.format(lat="intlat.bmad", root="int_hdr0",
        extra="  run.load_only = T\n"))
    check("interlude: and says nothing where nothing is cut",
          0.0 if "Interlude" not in r.stdout else 1.0, 0.5)

    # The refusals. Bmad applies an element's short-range wake once per passage, as one
    # kick for the element's length, so a cut element would take one kick per piece. The
    # walk does not cut such an element, as Tao does not cut where CSR or space charge
    # would see the cut. The transcribed Genesis interlude is validated by transcription
    # fidelity and steps by the element, so the knob is refused outright with that model.
    # The wake below has a table that is zero everywhere. The element carries a wake,
    # which is what the walk keys on, and the physics is the no-wake run's, so only the
    # rows may differ.
    (wd / "intlatw.bmad").write_text(INT_LAT.replace(
        "P1: pipe, l = 0.40",
        "P1: pipe, l = 0.40, sr_wake = {scale_with_length = T, z_long = "
        "{position_dependence = none, w = {-1e-6 0.0, 0.0 0.0, 1e-6 0.0,}}}"))
    run(exe, wd, "int_wk", INT_NML.format(lat="intlatw.bmad", root="int_wk",
        extra="  run.interlude_ds_step = 0.04\n"))
    wk = np.loadtxt(wd / "int_wk.diag.txt")
    check("interlude: a wake-carrying element is not cut", abs(inside(wk, 0.5, 0.9) - 0), 0.5,
          note=f"[{inside(wk, 0.5, 0.9)} rows inside it]")
    check("interlude: and the elements without a wake still are",
          abs(inside(wk, 0.3, 0.5) - 4), 0.5, note=f"[{inside(wk, 0.3, 0.5)}]")
    check("interlude: the zero wake left the physics alone, so only the rows differ",
          abs(wk[-1, 2] - p[-1, 2]) / p[-1, 2], 1e-9)

    (wd / "int_gen.toml").write_text(to_groups(INT_NML.format(lat="intlat.bmad", root="int_gen",
        extra='  run.interlude_ds_step = 0.04\n  run.interlude_model = "genesis"\n')))
    rg = subprocess.run([str(exe), "int_gen.toml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = rg.returncode != 0 and "TRANSCRIPTION OF GENESIS" in (rg.stdout + rg.stderr)
    check("interlude: refused with the transcribed Genesis interlude (0 = yes)",
          0.0 if refused else 1.0, 0.5)

    # The switch off. With bmad.sr_wakes_on off the walk applies no wake and cuts the
    # element like any other. The setup's record-count precompute must make the same
    # decision. If the precompute counts the element whole while the walk cuts it, the
    # run stops at the stats writer with more records than the precomputed count.
    wo = run(exe, wd, "int_wko", INT_NML.format(lat="intlatw.bmad", root="int_wko",
        extra="  run.interlude_ds_step = 0.04\n  bmad.sr_wakes_on = F\n"))
    wko = np.loadtxt(wd / "int_wko.diag.txt")
    check("interlude: with sr_wakes_on off the wake element is cut like any other",
          abs(inside(wko, 0.5, 0.9) - 9), 0.5, note=f"[{inside(wko, 0.5, 0.9)} rows inside it]")
    check("interlude: and its physics is the whole-element run's",
          abs(wko[-1, 2] - p[-1, 2]) / p[-1, 2], 1e-9)

    # Space charge on a Bmad-seam interlude is refused. The seam tracks one slice's bunch
    # at a time through track1_bunch, whose space-charge and CSR paths need the centroid
    # orbit of the whole beam. Without it Bmad prints that the centroid must be supplied
    # and returns the bunch untracked with no error set. A drift carrying the attribute
    # then loses its transport while the run reports success. The input file configures a
    # space-charge term, so the refusal reached is this one and not the one for a solver
    # with nothing to solve.
    (wd / "intlat_sc.bmad").write_text(INT_LAT.replace(
        "P1: pipe, l = 0.40", "P1: pipe, l = 0.40, space_charge_method = slice"))
    (wd / "int_sc.toml").write_text(to_groups(INT_NML.format(lat="intlat_sc.bmad", root="int_sc",
        extra="  bmad.csr_and_space_charge_on = T\n  space_charge.nz = 1\n"
              "  space_charge.rmax = 250e-6\n")))
    rs = subprocess.run([str(exe), "int_sc.toml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = rs.returncode != 0 and "NOT AN FEL ELEMENT" in (rs.stdout + rs.stderr)
    check("interlude: space charge on a seam interlude is refused (0 = yes)",
          0.0 if refused else 1.0, 0.5)
    check("interlude: and Bmad's own centroid message never printed (0 = yes)",
          0.0 if "CENTROID MUST BE SUPPLIED" not in (rs.stdout + rs.stderr) else 1.0, 0.5)

    # A zero-length element that can act is tracked. The walk skips a zero-length element
    # that Bmad tracks as the identity. A walk that skipped every zero-length element
    # without a wake would drop a thin kicker's kick with no message. The same kick
    # through a kicker of finite length is the control. A zero-length pipe in the
    # kicker's place gives the beam's own centroid, which the kick is taken against. The
    # pipe is the one element of the three that the walk skips, so its stats file is one
    # record shorter. A kicker is tracked whatever its kick, since its key carries a map.
    for tag, hk in (("thin", "HK: hkicker, l = 0, kick = 1e-5"),
                    ("long", "HK: hkicker, l = 1e-6, kick = 1e-5"),
                    ("none", "HK: pipe, l = 0")):
        (wd / f"intlat_{tag}.bmad").write_text(INT_LAT.replace(
            "P1: pipe, l = 0.40\n", f"P1: pipe, l = 0.40\n{hk}\n").replace(
            "SEG: line = (UND, QF, P1, UND)", "SEG: line = (UND, QF, HK, P1, UND)"))
        run(exe, wd, f"int_{tag}", INT_NML.format(lat=f"intlat_{tag}.bmad", root=f"int_{tag}",
            extra=""))
    px = {}
    nrec = {}
    for tag in ("thin", "long", "none"):
        with read_stats(wd / f"int_{tag}.stats.h5") as st:
            px[tag] = float(st["beam/slice/centroid"][-1, 0, 1])
            nrec[tag] = len(st.s)
    kick = px["thin"] - px["none"]
    check("thin kicker: a zero-length kicker kicks (1e-5 expected)", abs(kick - 1e-5) / 1e-5, 1e-3,
          note=f"[{kick:.4e}]")
    check("thin kicker: the same kick through a finite length agrees",
          abs(px["thin"] - px["long"]) / abs(kick), 1e-3)
    check("thin kicker: it takes a record, and the zero-length pipe is still skipped",
          abs((nrec["thin"] - nrec["none"]) - 1), 0.5, note=f"[{nrec['thin']} vs {nrec['none']}]")

    # A frame's particle records say where the frame is: sPosition is the frame's s and
    # timeOffset is the iteration's time. A writer that set the coords at the containing
    # element's upstream face would give a frame inside a piece the element's start and
    # entry time, while the file's own attributes gave the frame's position.
    worst_s, worst_t = 0.0, 0.0
    for fr in frames:
        with h5py.File(fr) as f:
            it = f["data"][sorted(f["data"])[0]]
            el = it["particles/electron"]
            def value(rec):
                return float(np.ravel(rec.attrs["value"] if "value" in rec.attrs else rec[()])[0])
            worst_s = max(worst_s, abs(value(el["sPosition"]) - float(np.ravel(f["lucifer/frame"].attrs["sPosition"])[0])))
            worst_t = max(worst_t, abs(value(el["timeOffset"]) - float(np.ravel(it.attrs["time"])[0])))
    check("frames: the particle sPosition is the frame's own [m]", worst_s, 1e-12)
    check("frames: the particle timeOffset is the iteration's time [s]", worst_t, 1e-18)

# refposition in the two Genesis-format files. Genesis's reader places slice i of a field
# file at refposition + (i-1) slicespacing (readFieldHDF5.cpp), so refposition is the
# in-bunch coordinate of slice 1, and the files this program writes say that: 0 for the
# pulse, whose slice 1 is the live window's tail, and the window length for the escaped
# file, whose slice 1 is the slot just past the head. The run below has an interlude, so
# the field is drifted along the lattice before the files are written. The check holds
# refposition unmoved by that drift, so that it never carries a lattice position.

REF_LAT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model,
     b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light,
     fel_method = averaged, ds_step = 0.015
P1: pipe, l = 0.40
SEG: line = (UND, P1, UND)
use, SEG
"""

GEN_LAT_REF = """UND: undulator = { lambdau=0.015000, nwig=30, aw=0.84853, helical=True};
SEG: line={UND};
"""

# Genesis imports the pulse file and writes it back with no tracking. &importfield sets
# the time window from the file (its refposition, slicecount and slicespacing), and
# &write dumps every slice of that window.
GEN_IMPORT_REF = """&setup
rootname=gref
lattice=ref.lat
beamline=SEG
lambda0=1e-10
gamma0=11357.82
delz=0.015
shotnoise=0
nbins = 8
&end

&importfield
file = {pulse}
&end

&write
field = GREF-copy
&end
"""


def grid_parity(exe, wd, examples):
    """The parity cost of a power-of-two grid, on the steady-state example's input file. The
    input file runs at its 256 points and at 255 with the cell held, the half width scaled by
    254/255, and the exit power's relative difference is recorded. On an even
    grid the axis falls between four cells, so a beam on the axis reads a centred field
    low by dx^2 over twice the waist squared through the bilinear gather, and the deposit
    returns the same deficit. An odd grid puts a cell on the axis. The example's dumps are
    left out, since the exit power is what the check reads, and the lattice is the
    example's own (docs/validation.md, the transverse grid)."""
    input_text = (examples / "steady_state" / "lucifer.toml").read_text()
    shutil.copy(examples / "aramis.bmad", wd / "aramis.bmad")
    input_text = input_text.replace('lattice = "../aramis.bmad"', 'lattice = "aramis.bmad"')
    input_text = "\n".join(line for line in input_text.splitlines()
                     if "dump_beam_at" not in line and "dump_field_at" not in line) + "\n"
    for needed in ('out_root = "steady_state"', "grid_n_pts = 256", "grid_half_width = 2e-4"):
        if needed not in input_text:
            print(f"FAIL: the steady-state example's input file no longer states {needed}")
            sys.exit(1)
    power = {}
    for n in (256, 255):
        text = input_text.replace('out_root = "steady_state"', f'out_root = "gp{n}"')
        if n == 255:
            text = text.replace("grid_n_pts = 256", "grid_n_pts = 255")
            text = text.replace("grid_half_width = 2e-4", f"grid_half_width = {2e-4 * 254 / 255:.12e}")
        (wd / f"gp{n}.toml").write_text(to_toml(text))
        r = subprocess.run([str(exe), f"gp{n}.toml"], cwd=wd, capture_output=True, text=True,
                           env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
        if r.returncode != 0:
            print(f"FAIL: gp{n} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
            sys.exit(1)
        with h5py.File(wd / f"gp{n}.stats.h5") as h5:
            power[n] = float(np.sum(h5["field/total/power"][-1]))
    check("grid parity: the steady-state example's exit power at 255 against 256 points, |rel|",
          abs(power[255] / power[256] - 1), 3e-2, note="(the cell held)")


def refposition_checks(exe, genesis, wd):
    """refposition in the pulse and escaped files, held against Genesis's own reader."""
    from beamphysics.wavefront import Wavefront

    print("--- refposition in the Genesis-format files:")
    (wd / "ref_int.bmad").write_text(REF_LAT)
    run(exe, wd, "dgi", NML.format(lat="ref_int.bmad", root="dgi", extra=""))

    def roots(root):
        with h5py.File(wd / f"{root}-pulse.fld.h5") as h5:
            p_ref = float(h5["refposition"][0])
            p_n = int(h5["slicecount"][0])
            dz = float(h5["slicespacing"][0])
        with h5py.File(wd / f"{root}-escaped.fld.h5") as h5:
            b_ref = float(h5["refposition"][0])
            b_n = int(h5["slicecount"][0])
            zt = np.array([h5[f"slice{k:06d}/z_transmit"][0] for k in range(1, b_n + 1)])
        return p_ref, p_n, dz, b_ref, b_n, zt

    p_ref, p_n, dz, b_ref, b_n, zt = roots("dgi")
    p_ref0, _, _, b_ref0, _, _ = roots("dg")
    nlive = p_n - b_n
    check("refposition: the pulse file's is 0 [m]", abs(p_ref), 1e-18)
    check("refposition: the escaped file's is the window length n_slice * slicespacing [m]",
          abs(b_ref - nlive * dz), 1e-18, note=f"[{b_ref:.3e}, {nlive} live slices]")
    check("refposition: both unmoved by the interlude's drift [m]",
          max(abs(p_ref - p_ref0), abs(b_ref - b_ref0)), 1e-18)
    # Tail to head: the escaped file's slice 1 left last, and the furthest ahead left first.
    check("refposition: the escaped file's z_transmit falls with its slice number (0 = yes)",
          0.0 if np.all(np.diff(zt) <= 0) and zt[0] > zt[-1] else 1.0, 0.5,
          note=f"[{b_n} banked, z_transmit {zt[0]:.3f} to {zt[-1]:.3f} m]")

    # openPMD-beamphysics reads the pulse as Genesis would: its first slice is at zmin 0,
    # and the live slices are the final openPMD field, slice by slice. The two agree to
    # 1.6e-16, since its Z0 and Bmad's agree and the dfl scale round-trips in double.
    w = Wavefront.from_genesis4(str(wd / "dgi-pulse.fld.h5"))
    fin = fieldio.read_field(wd / "dgi-final.wf.h5")
    live = np.transpose(w.Ex[:, :, :nlive], (2, 1, 0))
    rel = float(np.max(np.abs(live - fin["u"])) / np.max(np.abs(fin["u"])))
    check("refposition: Wavefront.from_genesis4 puts the pulse's first slice at zmin 0 [m]",
          abs(float(w.zmin)), 1e-18)
    check("refposition: the pulse's live slices are the final field, slice by slice, rel", rel, 1e-8)
    # Its banked slices are the escaped file's in order. Free-space propagation conserves
    # each slice's energy and the slices differ in it, so energy in order is the test.
    bank = fieldio.read_field(wd / "dgi-escaped.fld.h5")
    e_b = fieldio.field_power(bank["u"], bank["dx"], bank["dy"])
    e_w = fieldio.field_power(np.transpose(w.Ex[:, :, nlive:], (2, 1, 0)), w.dx, w.dy)
    lit = e_b > 1e-6 * e_b.max()
    check("refposition: pulse slice n_slice + i carries the energy of escaped slice i, rel",
          float(np.max(np.abs(e_w[lit] / e_b[lit] - 1))), 1e-8, note=f"[{lit.sum()} of {b_n} slices lit]")

    # Genesis itself: import the pulse and write it back. The re-dump places every slice
    # where the pulse said, with the pulse's own refposition.
    (wd / "ref.lat").write_text(GEN_LAT_REF)
    (wd / "gref.in").write_text(GEN_IMPORT_REF.format(pulse="dgi-pulse.fld.h5"))
    r = subprocess.run([str(genesis), "gref.in"], cwd=wd, capture_output=True, text=True,
                       env=dict(os.environ, FI_PROVIDER="tcp"))
    if r.returncode != 0:
        print(f"FAIL: genesis import of the pulse file:\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")
        sys.exit(1)
    pul = fieldio.read_field(wd / "dgi-pulse.fld.h5")
    cp = fieldio.read_field(wd / "GREF-copy.fld.h5")
    check("refposition: Genesis imports the pulse and re-dumps the same slice count",
          abs(cp["u"].shape[0] - pul["u"].shape[0]), 0.5, note=f"[{pul['u'].shape[0]} slices]")
    if cp["u"].shape == pul["u"].shape:
        check("refposition: Genesis's re-dump is the pulse, slice by slice, rel",
              float(np.max(np.abs(cp["u"] - pul["u"])) / np.max(np.abs(pul["u"]))), 1e-12)
    check("refposition: and carries the pulse's own refposition [m]", abs(cp["zmin"] - pul["zmin"]), 1e-18)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--latdir", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--genesis", required=True, help="genesis4 binary, for the import of the pulse file")
    ap.add_argument("--examples", required=True, help="the examples directory, for the grid parity check")
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    latdir = pathlib.Path(args.latdir).resolve()
    wd = pathlib.Path(args.workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "aramis_1seg.bmad").write_bytes((latdir / "aramis_1seg.bmad").read_bytes())
    (wd / "dg_wrap.bmad").write_text(WRAP)

    run(exe, wd, "dg", NML.format(lat="dg_wrap.bmad", root="dg", extra=""))

    # 0. The acceptance test is the standard's own validator, which knows no program. A
    # bmad-stats file conforms when the validator reports zero MUST failures, and
    # everything the validator checks comes from the spec alone. The named checks below
    # hold what it cannot see: the physics identities, and the declarations only this
    # writer knows it owes.
    rep = validate_bmad_stats.validate(wd / "dg.stats.h5")
    for level, where, msg in rep.rows:
        if level == "MUST":
            print(f"    {level} {where}: {msg}")
    check("bmad-stats: the generic validator reports zero MUST failures", rep.n_must, 0.5,
          note=f"[{sum(1 for r in rep.rows if r[0] != 'INFO')} findings in "
               f"{len(rep.rows)} rows]")

    with read_stats(wd / "dg.stats.h5") as st:
        nslice = int(st.run["n_slice"])

        # 0a. The file describes itself. A walk of every dataset checks that each says
        # what it is and what its dimensions run over. A generic reader depends on this,
        # so it is checked generically rather than name by name.
        bad = []
        for name, unit, label, axes, descrip in st.units_table():
            if unit == '' or descrip == '' or label == '':
                bad.append(f"{name}: no unit, label or description")
        for name, kind in st.kind.items():
            if kind == '' or st.group_description[name] == '':
                bad.append(f"{name}/: group carries no kind or description")
        check("stats: every dataset carries unit, label, description and axes", len(bad), 0.5,
              note=f"[{len(st.units_table())} datasets, {len(st.kind)} groups" +
                   (f"; first problem {bad[0]}]" if bad else "]"))

        # 0a2. The acceptance test is a generic load: label every dimension of every
        # dataset from @axes alone. The load fails if a name does not resolve to a
        # coordinate, if a dimension has no name, if a length disagrees with its
        # coordinate, or if one dataset names an axis twice. Any of these would force a
        # reader to guess from a length or to dedupe by rule. Nothing here knows a
        # dataset name.
        bad = []
        for name, unit, label, axes, descrip in st.units_table():
            shape = st[name].shape
            if len(axes) != len(shape):
                bad.append(f"{name}: {len(shape)} dimensions, {len(axes)} named")
                continue
            if len(set(axes)) != len(axes):
                bad.append(f"{name}: names an axis twice ({','.join(axes)})")
            for k, ax in enumerate(axes):
                if ax not in st.axis_names:
                    bad.append(f"{name}: @axes names {ax}, which is not an axis")
                elif len(st.coord(ax)) != shape[k]:
                    bad.append(f"{name}: dimension {k} is {shape[k]} long, {ax} is "
                               f"{len(st.coord(ax))}")
        check("stats: every dimension of every dataset labels itself", len(bad), 0.5,
              note=f"[{len(st.axis_names)} axes: {', '.join(st.axis_names)}" +
                   (f"; first problem {bad[0]}]" if bad else "]"))

        # 0a3. The record number is the axis, and s is a variable on it. s must be
        # non-decreasing. Where s repeats, the two records must sit in different
        # elements, which is what a zero-length element applying a wake kick gives: it
        # records at the plane where the element before it ended. Two records at one
        # plane inside one element would be a defect in the walk that this axis choice
        # hides. The zero-length wake lattice of check 5a gives this a real duplicate.
        ds = np.diff(st.s)
        dup = np.flatnonzero(ds == 0)
        unexplained = [int(i) for i in dup if st.ix_ele[i] == st.ix_ele[i + 1]]
        check("stats: s is non-decreasing along the record axis", int(np.sum(ds < 0)), 0.5)
        check("stats: every repeated s straddles an element boundary", len(unexplained), 0.5,
              note=f"[{len(dup)} repeats in {len(st.s)} records]")

        # 0a4. Provenance is held in datasets. HDF5 caps one attribute at 64 kB, and
        # 65495 bytes is the largest that writes. The echoed namelist is 12 kB and a real
        # lattice text 37 kB, so as attributes they would sit at half the cap, and a
        # larger text would go missing with only a warning. Every text in meta/ is a
        # dataset, and nothing in the file depends on the cap.
        meta = st.meta
        big = []

        def note_attrs(name, obj):
            for key, val in obj.attrs.items():
                if isinstance(val, bytes) and len(val) > 60000:
                    big.append(f"{name}@{key} is {len(val)} bytes")

        with h5py.File(wd / "dg.stats.h5") as h5:
            h5.visititems(note_attrs)
            note_attrs("/", h5)
            n_meta_attrs = len([k for k in h5["meta"].attrs if k not in ("kind", "description")])
        check("stats: no attribute is near HDF5's 64 kB cap", len(big), 0.5,
              note=f"[{'; '.join(big) if big else 'largest well under'}]")
        check("stats: meta/ holds its texts as datasets", n_meta_attrs, 0.5,
              note=f"[{len(meta)} datasets, input_echo {len(meta['input_echo'])} bytes]")

        # 0a5. The file does not overstate itself. dg tracks a wrapper lattice, so the
        # parser opens two files and lattice_source holds only the outer one.
        # n_lattice_files records the count, and this run is the case that needs it: a
        # reader of lattice_source alone would think it had the whole lattice.
        check("stats: n_lattice_files reports the wrapper's second file",
              abs(int(meta["n_lattice_files"]) - 2), 0.5,
              note=f"[{meta['n_lattice_files']} files, source "
                   f"{len(meta['lattice_source'])} bytes]")

        # 0a6. A stats file travels. Nothing identifies a person or a machine unless the
        # run asked for it, and a path is a basename.
        leaks = [k for k in ("user", "cwd") if k in meta]
        leaks += ["lattice_file has a directory"] if "/" in meta["lattice_file"] else []
        check("stats: the default file carries no machine-local values", len(leaks), 0.5,
              note=f"[{', '.join(leaks) if leaks else 'clean'}]")

        # 0a7. Two slice coordinates are exact and one is not. The grid is uniform in
        # time: the migration invariant carries each particle's own beta, and
        # z = -beta*c*(t-t_ref), so at a grid point the beta cancels and the arrival-time
        # separation is ct_slice/c with no beta anywhere. t_slice must therefore be
        # -ct_slice/c exactly. A spurious factor beta0 would put it off by 3.9e-9. z_slice
        # needs a reference, so it equals beta0*ct_slice only to the rounding of two ways
        # of forming beta0.
        c_light = 299792458.0
        p0c = float(st.run["p0c"])
        beta0 = p0c / np.sqrt(p0c**2 + M_ELECTRON**2)
        ct = st.ct_slice
        check("stats: t_slice is -ct_slice/c exactly (abs)",
              float(np.max(np.abs(st.t_slice + ct / c_light))), 1e-30)
        rel = np.abs(st.z_slice[1:] / (beta0 * ct[1:]) - 1)
        check("stats: z_slice is beta0*ct_slice", float(np.max(rel)), 1e-15,
              note=f"[beta0 = 1 - {1 - beta0:.3e}]")
        check("stats: the slice axis is the index", abs(len(st.slice) - nslice), 0.5,
              note=f"[{nslice} slices, positions ct_slice, t_slice, z_slice on it]")

        # 0a8. The envelope data are the particles'. rel_max and rel_min are order
        # statistics relative to the stored centroid. The dump at the UND end holds the
        # same particles the record saw, so the position entries must match to the bit:
        # the file stores x verbatim, and IEEE subtraction of the same centroid is
        # deterministic. The momentum entries cross the dump's unit round trip
        # (px*p0c on write, /p0c on read), so they get an ulp-scale tolerance.
        irec = int(np.flatnonzero(st.at_end)[0])
        rel_hi = st["beam/slice/rel_max"][irec]
        rel_lo = st["beam/slice/rel_min"][irec]
        cen_r = st["beam/slice/centroid"][irec]
        dump_sl = read_slices(wd / "dg-at1-UND.beam.h5", wavelength=1e-10)
        p0_mc_l = float(st.run["p0c"]) / M_ELECTRON
        worst_pos, worst_mom = 0.0, 0.0
        for isl, sd in enumerate(dump_sl):
            if sd["n"] == 0:
                continue
            for j, arr in ((0, sd["x"]), (2, sd["y"])):
                worst_pos = max(worst_pos,
                                abs(np.max(arr - cen_r[isl, j]) - rel_hi[isl, j]),
                                abs(np.min(arr - cen_r[isl, j]) - rel_lo[isl, j]))
            for j, arr in ((1, sd["px"] / p0_mc_l), (3, sd["py"] / p0_mc_l),
                           (5, sd["pz"])):
                sc = max(abs(rel_hi[isl, j]), abs(rel_lo[isl, j]), 1e-30)
                worst_mom = max(worst_mom,
                                abs(np.max(arr - cen_r[isl, j]) - rel_hi[isl, j]) / sc,
                                abs(np.min(arr - cen_r[isl, j]) - rel_lo[isl, j]) / sc)
        check("envelope: position extremes vs the dump's particles, exactly (abs)",
              worst_pos, 1e-30, note=f"[{len(dump_sl)} slices, stored centroid]")
        check("envelope: momentum extremes across the dump's unit round trip",
              worst_mom, 1e-12)

        # And the time entry is the z entry through the exact map t - <t> =
        # -(z - <z>)/(beta0 c): the largest time offset is the smallest z one.
        beta0_l = p0_mc_l / np.sqrt(p0_mc_l**2 + 1)
        c_l = 299792458.0
        rmx = st["beam/slice/rel_max"]
        rmn = st["beam/slice/rel_min"]
        fin = np.isfinite(rmx[..., 6])
        worst_t = float(np.max(np.abs(rmx[..., 6][fin] +
                                      rmn[..., 4][fin] / (beta0_l * c_l))))
        worst_t = max(worst_t, float(np.max(np.abs(rmn[..., 6][fin] +
                                                   rmx[..., 4][fin] / (beta0_l * c_l)))))
        check("envelope: the t entry is the z entry through -dz/(beta0 c), exactly (abs)",
              worst_t, 1e-30)

        # 0b. The join. at_element_end selects the element-end rows, and every record
        # sits inside the element its ix_ele names. An off-by-one in either would break
        # every layout plot and no physics check would see it.
        n_end = int(st.at_end.sum())
        check("stats: at_element_end selects exactly the element ends",
              abs(n_end - int(st.run["n_element_end"])), 0.5,
              note=f"[{n_end} of {len(st.s)} records]")
        check("stats: element-end arrays are aligned with the mask",
              abs(st["beam/bunch/centroid"].shape[0] - n_end), 0.5)
        s1, s2 = st["lattice/s_start"], st["lattice/s_end"]
        ix = st.ix_ele
        outside = int(np.sum((st.s < s1[ix] - 1e-9) | (st.s > s2[ix] + 1e-9)))
        check("stats: every record sits inside the element ix_ele names", outside, 0.5,
              note=f"[{len(st.s)} records against {len(s1)} lattice rows]")

        # 1. Reconstruction at the element end, which is the last record. The entry is an
        # index into the group's own axis, looked up rather than assumed. That lookup is
        # all that stands between a reader and a silent transposition. The projected
        # planes and the normal modes are separate groups over separate axes, because an
        # eigen-emittance and a projected emittance are different quantities.
        plane = {name: i for i, name in enumerate(st.coord("plane"))}
        mode = {name: i for i, name in enumerate(st.coord("mode"))}
        worst = 0.0
        for isl in range(nslice):
            bp = bunch_params_at(st, -1, isl)
            for m in ("x", "y", "z"):
                for pn in ("beta", "alpha", "emit", "norm_emit", "sigma", "sigma_p", "eta", "etap"):
                    stored = st[f"beam/slice_twiss/twiss/{pn}"][0, isl, plane[m]]
                    worst = max(worst, abs(bp[m][pn] - stored) / max(abs(stored), 1e-30))
            # Bmad's mat_eigen labels modes by eigenvector structure rather than by
            # magnitude, which coords/mode says, so this compares them as a sorted set.
            modes = sorted(st["beam/slice_twiss/modes/emit"][0, isl, mode[m]]
                           for m in ("a", "b", "c"))
            mine = sorted(bp[m]["emit"] for m in ("a", "b", "c"))
            for s0, m0 in zip(modes, mine):
                worst = max(worst, abs(m0 - s0) / max(abs(s0), 1e-30))
        # numpy's eig and mat_eigen differ by 4.6e-8. The projected planes agree exactly.
        check("stats: bunch_params reconstruction vs stored calc_bunch_params", worst, 1e-6)

        # 1b. The whole-window row is assembled from the per-slice moments by the
        # pooled-covariance identity rather than summed over particles. The identity is
        # re-implemented here, independently, with the local-to-global z map: a shift of
        # beta*(is-1)*spacing, plus the pz shear the map carries because beta is the
        # particle's own. Against the particle sum the row agrees to 4.0e-12 on this
        # configuration and to 5.0e-11 on the 96-slice SASE example.
        spacing = float(st.run["slice_spacing"])
        p0_mc = float(st.run["p0c"]) / M_ELECTRON
        worst_pool = 0.0
        worst_nobg = 0.0
        for ie in range(st["beam/bunch/sigma"].shape[0]):
            cen = st["beam/slice_twiss/centroid"][ie]
            sig = st["beam/slice_twiss/sigma"][ie]
            wsl = st["beam/slice_twiss/charge_live"][ie]
            wtot = wsl.sum()
            if wtot <= 0:
                continue
            pmc = p0_mc * (1.0 + cen[:, 5])
            beta = pmc / np.sqrt(pmc**2 + 1)
            ell = np.arange(len(wsl)) * spacing
            shear = ell * p0_mc / np.sqrt(pmc**2 + 1)**3
            cg = cen.copy()
            cg[:, 4] = cg[:, 4] + beta * ell
            m = (wsl[:, None] * cg).sum(0) / wtot
            pooled = np.zeros((6, 6))
            nobg = np.zeros((6, 6))
            for isl in range(len(wsl)):
                s6 = sig[isl].copy()
                k = shear[isl]
                s0 = sig[isl]
                s6[4, 4] = s0[4, 4] + 2 * k * s0[4, 5] + k * k * s0[5, 5]
                for j in range(6):
                    if j == 4:
                        continue
                    s6[4, j] = s0[4, j] + k * s0[5, j]
                    s6[j, 4] = s6[4, j]
                d = cg[isl] - m
                pooled += wsl[isl] * (s6 + np.outer(d, d))
                nobg += wsl[isl] * s6
            pooled /= wtot
            nobg /= wtot
            st_c = st["beam/bunch/centroid"][ie]
            st_s = st["beam/bunch/sigma"][ie]
            den = np.maximum(np.abs(st_s), np.abs(pooled))
            worst_pool = max(worst_pool, float(np.max(np.where(den > 0, np.abs(st_s - pooled) / np.where(den > 0, den, 1), 0.0))))
            dc = np.maximum(np.abs(st_c), np.abs(m))
            worst_pool = max(worst_pool, float(np.max(np.where(dc > 0, np.abs(st_c - m) / np.where(dc > 0, dc, 1), 0.0))))
            den = np.maximum(np.abs(st_s), np.abs(nobg))
            worst_nobg = max(worst_nobg, float(np.max(np.where(den > 0, np.abs(st_s - nobg) / np.where(den > 0, den, 1), 0.0))))
        check("stats: whole-window row vs the pooled-covariance identity", worst_pool, 1e-9)
        # Dropping the between-group term (m_s - m)(m_s - m)^T, which makes the identity
        # exact, must be caught and not absorbed. Without the term sigma(5,5) collapses
        # from the window length squared to a slice length squared.
        check("stats: the between-group term is load-bearing (0 = confirmed)",
              0.0 if worst_nobg > 1e-3 else 1.0, 0.5,
              note=f"[dropping it moves the row by {worst_nobg:.2e}]")

        # 1c. The element-end twiss rests on the record's own moments. At an element end
        # beam/slice_twiss/centroid and sigma must equal beam/slice/centroid and sigma at
        # the record the mask selects. They agree exactly, since one is written from the
        # other.
        worst_ef = float(np.max(np.abs(st["beam/slice/centroid"][st.at_end] -
                                      st["beam/slice_twiss/centroid"])))
        worst_ef = max(worst_ef, float(np.max(np.abs(st["beam/slice/sigma"][st.at_end] -
                                                    st["beam/slice_twiss/sigma"]))))
        check("stats: the element-end moments ARE the record's (abs)", worst_ef, 1e-30)

        # The field's theta moments are computed where the file says they are: (nz, ns)
        # per slice, at every element end that has field, plus the initial record. The
        # walk takes the initial record with angles so the starting state is complete.
        # An empty slice has no moments to compute and says so.
        valid = st["field/x/angle_moments_valid"].astype(bool)
        has_field = st["field/x/power"] > 0
        rows = np.zeros(len(st.s), bool)
        rows[0] = True
        rows |= st.at_end
        want = has_field & rows[:, None]
        check("stats: angle moments valid exactly where they were computed",
              int(np.sum(valid != want)), 0.5,
              note=f"[{int(valid.sum())} of {valid.size} slice-records]")

        # Everywhere else the moments are NaN rather than a zero that reads as an answer.
        # NaN lets a consumer tell an empty slice's missing moments from a real value.
        cen = st["field/x/centroid"]
        check("stats: theta moments are NaN where they were not computed",
              0.0 if np.all(np.isnan(cen[..., 1][~valid])) else 1.0, 0.5)
        check("stats: and finite where they were",
              0.0 if np.all(np.isfinite(cen[..., 1][valid])) else 1.0, 0.5)

        f_cen = st["field/x/centroid"][-1]
        f_sig = st["field/x/sigma"][-1]
        f_en = st["field/total/energy"][-1]

    # 2. Banked energy == ledger U_escaped.
    led = np.loadtxt(wd / "dg.ledger.txt")
    with h5py.File(wd / "dg-escaped.fld.h5") as h5:
        nb = int(h5["slicecount"][0])
        pms = np.array([h5[f"slice{k:06d}/wavefront_params"][:] for k in range(1, nb + 1)])
        zt = np.array([h5[f"slice{k:06d}/z_transmit"][0] for k in range(1, nb + 1)])
    e_banked = pms[:, 20].sum()
    check("bank: banked slice energies vs ledger U_escaped, rel",
          abs(e_banked / led[-1, 4] - 1), 1e-12)

    # 2b. The pulse record. The window is fixed and light leaves it at the head, so the
    # window's energy alone falls as the pulse slips out; escaped_energy is the light that
    # left and pulse_energy the whole pulse. Each record's pulse is the window's energy
    # plus the escaped light, the escaped light at the exit is the bank's and the
    # ledger's, and with the gain still growing over this one seeded segment the pulse
    # never falls where the window does.
    with read_stats(wd / "dg.stats.h5") as st:
        e_win = np.nansum(st["field/total/energy"], axis=1)
        e_esc = st["field/total/escaped_energy"]
        e_pul = st["field/total/pulse_energy"]
    check("pulse: pulse_energy vs window energy + escaped_energy, every record, max rel",
          float(np.max(np.abs(e_pul - (e_win + e_esc)) / np.maximum(e_pul, 1e-300))), 1e-12)
    check("pulse: escaped_energy at the exit vs the banked slices' energies, rel",
          abs(e_esc[-1] / e_banked - 1), 1e-12)
    check("pulse: escaped_energy at the exit vs the ledger's U_escaped, rel",
          abs(e_esc[-1] / led[-1, 4] - 1), 1e-12)
    falls = int(np.sum(np.diff(e_win) < 0))
    check("pulse: the largest fall of pulse_energy between records, rel (the window falls)",
          float(max(0.0, -np.min(np.diff(e_pul) / np.maximum(e_pul[:-1], 1e-300)))), 1e-12,
          note=f"[the window's energy falls at {falls} of {len(e_win) - 1} steps]")
    if falls == 0:
        check("pulse: the window's energy falls somewhere, so the check above is not empty", 1.0, 0.5)

    # A window with room ahead of the bunch for all the slippage the run applies keeps the
    # light the bunch makes: escaped_energy is zero at every record. The one segment slips the
    # light 267 wavelengths, its 266 periods and the schedule's end-of-lattice +1, and the
    # window gives 600 slices of 1e-10 m to an 8e-10 m bunch. The run is unseeded, since a
    # seed fills the whole window, the empty slices ahead of the bunch included, and that
    # seed light does slip out.
    run(exe, wd, "dgw", NML.format(lat="aramis_1seg.bmad", root="dgw", extra="").replace(
        "window_length = 8e-10", "window_length = 6e-8").replace("seed_power = 1e4", "seed_power = 0"))
    with read_stats(wd / "dgw.stats.h5") as st:
        esc_w = float(np.max(st["field/total/escaped_energy"]))
        n_slice_w = len(st.coord("slice"))
    check("pulse: a window holding the slippage loses no light, escaped_energy max [J]", esc_w, 0.0,
          note=f"[{n_slice_w} slices]")

    # 3. Analytic (moment-map) vs numerical (FFT) propagation to the exit plane.
    with h5py.File(wd / "dg-pulse.fld.h5") as h5:
        n_all = int(h5["slicecount"][0])
        nx = int(h5["gridpoints"][0])
        dx = float(h5["gridsize"][0])
        nlive = n_all - nb
        ax = (np.arange(nx) - 0.5 * (nx - 1)) * dx
        sxx_num, syy_num, w_num = [], [], []
        cx_num, cy_num = [], []
        for j in range(nlive + 1, n_all + 1):
            g = h5[f"slice{j:06d}"]
            fr = g["field-real"][:].reshape(nx, nx)   # x fastest: columns = x.
            fi = g["field-imag"][:].reshape(nx, nx)
            I = fr**2 + fi**2
            w = I.sum()
            w_num.append(w)
            if w == 0:
                sxx_num.append(0); syy_num.append(0); cx_num.append(0); cy_num.append(0)
                continue
            cx = (I * ax[None, :]).sum() / w
            cy = (I * ax[:, None]).sum() / w
            cx_num.append(cx); cy_num.append(cy)
            sxx_num.append((I * (ax[None, :] - cx) ** 2).sum() / w)
            syy_num.append((I * (ax[:, None] - cy) ** 2).sum() / w)
    # Pulse index nlive+i holds escaped slice i: both files are in Genesis's order,
    # tail to head, and so is zt, read from the escaped file above.
    sxx_num = np.array(sxx_num)
    syy_num = np.array(syy_num)
    w_num = np.array(w_num)

    z_end = led[-1, 0]  # Ledger rows end at the segment exit for this one-element line.
    dz = z_end - zt
    # Flattened 4x4 column-major: (i,j) -> 4*(j-1)+(i-1).
    s11, s21, s22 = pms[:, 4 + 0], pms[:, 4 + 1], pms[:, 4 + 5]
    s33, s43, s44 = pms[:, 4 + 10], pms[:, 4 + 11], pms[:, 4 + 15]
    sxx_ana = s11 + 2 * dz * s21 + dz**2 * s22
    syy_ana = s33 + 2 * dz * s43 + dz**2 * s44
    m = w_num > 1e-6 * w_num.max()
    rel = np.abs(np.sqrt(sxx_num[m] / sxx_ana[m]) - 1)
    rel = np.maximum(rel, np.abs(np.sqrt(syy_num[m] / syy_ana[m]) - 1))
    # The largest difference is 8.4e-3, from the grid-edge wrap of the FFT route on
    # divergent noise slices.
    check("bank: analytic vs FFT-propagated rms at exit, max rel", rel.max(), 2e-2)

    # 4. Pulse pooling: the pooled per-slice params against pulse moments computed
    # directly. The banked side is the params propagated analytically to the exit. The
    # live side is the last stats record, where the angles are valid. sigma_xx of the
    # pooled result is compared against pooling the pulse file's numerically computed
    # per-slice moments.
    cen_b = np.stack([pms[:, 0] + dz * pms[:, 1], pms[:, 1],
                      pms[:, 2] + dz * pms[:, 3], pms[:, 3]], axis=1)
    sig_b = np.zeros((nb, 4, 4))
    sig_b[:, 0, 0] = sxx_ana;  sig_b[:, 2, 2] = syy_ana
    en_b = pms[:, 20]
    cen_l = f_cen;  sig_l = f_sig;  en_l = f_en
    _, pooled_ana, _ = pool_wavefront(np.vstack([cen_l, cen_b]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_b.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_b]))
    # Numerical pooling over the same banked slices, from the pulse file.
    cen_n = np.stack([np.array(cx_num), np.zeros(nb), np.array(cy_num), np.zeros(nb)], axis=1)
    sig_n = np.zeros((nb, 4, 4))
    sig_n[:, 0, 0] = sxx_num;  sig_n[:, 2, 2] = syy_num
    en_n = w_num * (en_b.sum() / max(w_num.sum(), 1e-300))   # Same total, per-slice proportional.
    _, pooled_num, _ = pool_wavefront(np.vstack([cen_l, cen_n]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_n.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_n]))
    check("pulse: pooled sigma_xx analytic vs numerical routes, rel",
          abs(pooled_ana[0, 0] / pooled_num[0, 0] - 1), 2e-2)

    # 5. Thread invariance. Every dataset of the three files is identical at 1 and 8
    # threads. The comparison is by dataset rather than by raw bytes, since HDF5 object
    # headers embed creation times and two files with bit-identical data differ as byte
    # streams. The invariance holds for the data and not for the container metadata.
    run(exe, wd, "dg1", NML.format(lat="dg_wrap.bmad", root="dg1", extra=""), threads="1")

    same = all(h5_identical(wd / f"dg{s}", wd / f"dg1{s}")
               for s in (".stats.h5", "-escaped.fld.h5", "-pulse.fld.h5"))
    check("thread invariance: stats/escaped/pulse data identical 1 vs 8 (0 = yes)",
          0.0 if same else 1.0, 0.5)

    # 5a. Zero-length wake elements of both polarities. This lattice gives the
    # repeated-s check of 0a3 a real repeat. On a lattice without one that check sees
    # zero repeats and tests nothing. WKF can kick, so it is tracked and its record
    # lands on the plane where UND ended, repeating coords/s. The two WKT pipes cannot
    # kick, so they are skipped and the file must be identical to one whose zero-length
    # pipes carry no wake at all.
    (wd / "zlw.bmad").write_text(ZL_WAKE.format(wkt=WKT_WAKE))
    (wd / "zlwn.bmad").write_text(ZL_WAKE.format(wkt=WKT_NONE))
    run(exe, wd, "zlw", NML.format(lat="zlw.bmad", root="zlw", extra=""), threads="4")
    run(exe, wd, "zlwn", NML.format(lat="zlwn.bmad", root="zlwn", extra=""), threads="4")

    with read_stats(wd / "zlw.stats.h5") as st:
        ds = np.diff(st.s)
        dup = np.flatnonzero(ds == 0)
        unexplained = [int(i) for i in dup if st.ix_ele[i] == st.ix_ele[i + 1]]
        names = st.ele_name
        check("zero-length wake: coords/s repeats where one CAN kick", 1.0 / max(len(dup), 1) - 1.0,
              0.5, note=f"[{len(dup)} repeat(s), at " +
                        ", ".join(f"{names[i]}->{names[i+1]}" for i in dup) + "]")
        check("zero-length wake: and every repeat straddles an element boundary",
              len(unexplained), 0.5)
        check("zero-length wake: s is still non-decreasing", int(np.sum(ds < 0)), 0.5)

    same = h5_identical(wd / "zlw.stats.h5", wd / "zlwn.stats.h5")
    check("zero-length wake: one that CANNOT kick leaves the file unchanged (0 = yes)",
          0.0 if same else 1.0, 0.5)

    # 5b. The environment switch: opt-in, so the fields appear only when asked for.
    run(exe, wd, "dgenv", NML.format(lat="dg_wrap.bmad", root="dgenv",
        extra="  run.record_environment = T\n"), threads="2")
    with read_stats(wd / "dgenv.stats.h5") as st:
        got = sorted(k for k in st.meta if k in ("user", "cwd"))
    check("provenance: run.record_environment restores user and cwd",
          abs(len(got) - 2), 0.5, note=f"[{', '.join(got) if got else 'neither'}]")

    # 6. Refusal of a dump list entry that matches nothing.
    (wd / "dg_bad.toml").write_text(to_groups(NML.format(lat="dg_wrap.bmad", root="dg_bad",
                                             extra='  dump_beam_at(2) = "NO_SUCH_ELEMENT"\n')))
    r = subprocess.run([str(exe), "dg_bad.toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = r.returncode != 0 and "NO_SUCH_ELEMENT" in (r.stdout + r.stderr)
    check("refusal: unknown dump_beam_at element refused (0 = yes)",
          0.0 if refused else 1.0, 0.5)

    # 7. The pre-run (docs/user-guide.md). The input file loads 512 macroparticles per slice in
    # beamlets of 8, which is 64 beamlets and half the measured floor, so the load
    # standing and the thin-slice list both have something to report. The count is
    # arithmetic and is checked as arithmetic. The estimate is one machine's, and is
    # checked only for being present and positive and for naming the machine.

    r = run(exe, wd, "pre", PRE_NML.format(root="pre", extra=""), threads="4")
    out = r.stdout
    work = next((l for l in out.splitlines() if l.startswith(" Work ")), "")
    m = re.search(r"(\d+) FEL steps, (\d+) macroparticles, ([\d.E+]+) particle-steps, "
                  r"([\d.E+]+) grid points per step", work)
    check("pre-run: the header states the work", 0.0 if m else 1.0, 0.5,
          note=f"[{work.strip()}]")
    if m:
        nstep, npart, pstep = int(m.group(1)), int(m.group(2)), float(m.group(3))
        # The line carries three significant figures, so the identity is checked to
        # what it prints and not to what it computed.
        check("pre-run: particle-steps is the product of the two counts",
              abs(pstep - nstep * npart) / max(pstep, 1.0), 1e-2)

    est = next((l for l in out.splitlines() if l.startswith(" Estimate ")), "")
    m = re.search(r"Estimate\s+([\d.]+) s of walk", est)
    check("pre-run: the header estimates the walk and names the machine",
          0.0 if (m and float(m.group(1)) > 0 and "M3 Max" in est) else 1.0, 0.5,
          note=f"[{est.strip()}]")

    check("pre-run: the load standing reads under the measured floor",
          0.0 if "Under the measured load" in out else 1.0, 0.5)
    check("pre-run: the beamlet size is stated by the harmonics it carries",
          0.0 if "Shot noise on harmonics 1 to 3" in out else 1.0, 0.5)

    thin = [l for l in out.splitlines() if l.strip().startswith("slice ")]
    check("pre-run: the footer lists the thin slices", 0.0 if thin else 1.0, 0.5,
          note=f"[{len(thin)} listed]")
    check("pre-run: and names the remedy",
          0.0 if 'load_mode = "sample"' in out else 1.0, 0.5)
    cost = next((l for l in out.splitlines() if l.startswith(" Cost ")), "")
    check("pre-run: the footer prints the clock beside the estimate",
          0.0 if ("measured" in cost and "estimated" in cost) else 1.0, 0.5,
          note=f"[{cost.strip()}]")

    # load_only prints the same header and stops: the pre-run is the header without the
    # tracking, so the two headers have to agree line for line.

    r2 = run(exe, wd, "pre2", PRE_NML.format(root="pre2", extra="  load_only = T\n"),
             threads="4")
    head = lambda s, root: [l.replace(root, "R") for l in s.splitlines()
                            if l.startswith((" Work ", " Estimate ", " Load "))]
    check("pre-run: load_only gives the same header as the run does",
          0.0 if head(out, "pre") == head(r2.stdout, "pre2") else 1.0, 0.5)
    check("pre-run: load_only tracks nothing", 0.0 if " Cost " not in r2.stdout else 1.0, 0.5)

    # The refusal: a harmonic above (beamlet_size - 1)/2 has no shot noise at its own
    # frequency, so a dark-start run there would report a number with no startup behind
    # it. The message names the beamlet size that would carry it.

    (wd / "pre_h.toml").write_text(to_groups(PRE_NML.format(root="pre_h",
                                  extra="  harmonics = 1, 5\n")))
    r3 = subprocess.run([str(exe), "pre_h.toml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    said = r3.stdout + r3.stderr
    refused = r3.returncode != 0 and "BEAMLET_SIZE TO 11" in said
    check("refusal: harmonic 5 on beamlets of 8 refused, naming the size that carries it",
          0.0 if refused else 1.0, 0.5)

    # 8. The frame series (docs/reading-output.md). The input file migrates, so labels have to
    # survive particles changing slice, and its comb is coarse enough to keep the frame
    # count small. The element-end dump is asked for by name so the two writers can be
    # compared on one position.

    fr_extra = ('  run.dump_at_comb = T\n  run.comb_ds_save = 1.0\n'
                '  run.migrate = T\n  run.dump_beam_at = "UND"\n')
    r = run(exe, wd, "fr", FRAME_NML.format(root="fr", extra=fr_extra), threads="4")
    frames = sorted(wd.glob("fr-[0-9]*.beam.h5"))
    wframes = sorted(wd.glob("fr-[0-9]*.wf.h5"))
    with read_stats(wd / "fr.stats.h5") as st:
        nrec = len(st.s)
    check("frames: one beam frame per stats record", abs(len(frames) - nrec), 0.5,
          note=f"[{len(frames)} frames, {nrec} records]")
    check("frames: one field frame per stats record", abs(len(wframes) - nrec), 0.5)

    # The cost line counts them before the run, from the schedule alone.
    m = re.search(r"^ Frames\s+(\d+) at the comb", r.stdout, re.M)
    check("frames: the header's count matches the files written",
          abs(int(m.group(1)) - len(frames)) if m else 1.0, 0.5,
          note=f"[header {m.group(1) if m else 'absent'}]")

    # A frame at an element end and that element's own dump are one writer at one state.
    end_dump = wd / "fr-at1-UND.beam.h5"
    same = any(h5_identical(f, end_dump) for f in frames) if end_dump.exists() else False
    check("frames: a frame at an element end equals that element's own dump (0 = yes)",
          0.0 if same else 1.0, 0.5)

    # Labels: unique inside a frame, and no label appears that was not there before.
    def labels(path):
        with h5py.File(path) as h:
            g = h[f"data/{list(h['data'])[0]}/particles/electron"]
            ids = g["id"][()]
            npart = g["particlePatches/numParticles"][()]
            slot, k = {}, 0
            for islice, cnt in enumerate(npart):
                for _ in range(int(cnt)):
                    slot[int(ids[k])] = islice
                    k += 1
            return ids, slot

    dup = 0
    for f in frames:
        ids, _ = labels(f)
        dup += len(ids) - len(set(ids.tolist()))
    check("frames: every label appears once in every frame", dup, 0.5)

    _, first = labels(frames[0])
    _, last = labels(frames[-1])
    moved = sum(1 for i in first if i in last and first[i] != last[i])
    check("frames: no label appears that the first frame did not carry",
          len(set(last) - set(first)), 0.5,
          note=f"[{moved} labels changed slice, {len(set(first) - set(last))} left the window]")
    check("frames: migration moved labels between slices, so the check has something to see",
          0.0 if moved > 0 else 1.0, 0.5)

    # The attributes say where the frame was taken. The element the run reports for a
    # record is the element the frame must name.
    with read_stats(wd / "fr.stats.h5") as st:
        names = [n.strip() for n in st.ele_name]
        s_rec = np.asarray(st.s).copy()
    bad_attr = 0
    for i, f in enumerate(frames):
        with h5py.File(f) as h:
            got = h["lucifer/frame"].attrs["elementName"]
            got = got.decode() if isinstance(got, bytes) else str(got)
            spos = float(np.ravel(h["lucifer/frame"].attrs["sPosition"])[0])
        if got.strip() != names[i] or abs(spos - s_rec[i]) > 1e-9:
            bad_attr += 1
    check("frames: each frame names the element and s the stats row carries", bad_attr, 0.5)

    # An FEL frame carries the undulator a reader rebuilds the wiggle from.
    with h5py.File(frames[-1]) as h:
        has_und = all(k in h["lucifer/frame"].attrs for k in ("aw", "ku", "helical", "tilt", "felMethod"))
    check("frames: an FEL frame carries aw, ku and the method (0 = yes)",
          0.0 if has_und else 1.0, 0.5)

    # A run that writes frames must equal the same run without them. Writing a frame
    # must not change the state it records, which is why the field's rotation is put
    # back.
    run(exe, wd, "nofr", FRAME_NML.format(root="nofr",
        extra='  run.comb_ds_save = 1.0\n  run.migrate = T\n'), threads="4")
    check("frames: the run is identical to the same run without them (0 = yes)",
          0.0 if h5_identical(wd / "fr.stats.h5", wd / "nofr.stats.h5") else 1.0, 0.5)

    # The unaveraged method, written and named. A frame inside an unaveraged segment is
    # taken where px carries the resolved undulator oscillation. The writer allows that
    # chart because the frame names its method. A chart assertion inside the
    # slice-to-bunch conversion would refuse the whole segment interior, and an
    # unaveraged run could then finish only with dump_at_comb = F.

    (wd / "dg_unavg.bmad").write_text(WRAP)
    ua_extra = ('  run.dump_at_comb = T\n  run.comb_ds_save = 1.0\n'
                '  run.device = "off"\n')
    # Run it without run()'s own exit, so that a mode which cannot write a frame is a
    # named failing check here rather than an abort with the reason buried in a log.
    ua_text = (FRAME_NML.format(root="ua", extra=ua_extra)
               .replace('lattice = "aramis_1seg.bmad"', 'lattice = "dg_unavg.bmad"'))
    (wd / "ua.toml").write_text(to_groups(ua_text))
    r = subprocess.run([str(exe), "ua.toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    ua_frames = sorted(wd.glob("ua-[0-9]*.beam.h5"))
    check("frames: an unaveraged run writes them at all (exit 0)", float(r.returncode), 0.5,
          note="" if r.returncode == 0 else f"[{r.stdout.strip().splitlines()[-1][:70]}]")
    check("frames: the unaveraged run wrote a series",
          0.0 if len(ua_frames) > 1 else 1.0, 0.5, note=f"[{len(ua_frames)} frames]")

    if len(ua_frames) > 1:
        def method_of(f):
            with h5py.File(f) as h:
                m = h["lucifer/frame"].attrs.get("felMethod", b"")
                return m.decode() if isinstance(m, bytes) else str(m)
        bad_name = sum(1 for f in ua_frames[1:] if method_of(f) != "Unaveraged")
        check("frames: every unaveraged frame names its method", bad_name, 0.5)

    # Restarting from a frame needs no check of its own. A frame at an element end is
    # dataset-identical to that element's dump, which the check above holds, and
    # check_program holds the windowed composition for those dumps. A frame taken inside
    # an element is not a restart point at all, since run.track_start names elements
    # and not positions.
    #
    # What needs a check is that frames land inside an element, which is the reason for
    # the series. A comb that fired only at element ends would pass every check above
    # and give a movie with one frame per undulator.
    inside = 0
    for i, f in enumerate(frames):
        with h5py.File(f) as h:
            spos = float(np.ravel(h["lucifer/frame"].attrs["sPosition"])[0])
            nm = h["lucifer/frame"].attrs["elementName"]
            nm = (nm.decode() if isinstance(nm, bytes) else str(nm)).strip()
        if nm == "UND" and 0.0 < spos < 3.99:
            inside += 1
    check("frames: the series lands inside the undulator, not only at its ends",
          0.0 if inside >= 2 else 1.0, 0.5, note=f"[{inside} of {len(frames)} inside UND]")

    # 9. A frame against its stats row. These three identities pin the phase convention,
    # the time-order rotation and the slice indexing end to end: a picture drawn from a
    # frame and a number read from the row have to be the same run.

    rr_extra = ('  run.dump_at_comb = T\n  run.comb_ds_save = 1.0\n'
                '  run.dump_reduced = T\n')
    run(exe, wd, "id", FRAME_NML.format(root="id", extra=rr_extra), threads="4")
    idf = particle_frames(wd, "id")
    idw = sorted(wd.glob("id-[0-9]*.wf.h5"))

    with read_stats(wd / "id.stats.h5") as st:
        b_row = np.asarray(st["beam/slice/bunching"]).copy()
        c_row = np.asarray(st["beam/slice/centroid"]).copy()
        p_row = np.asarray(st["field/total/power"]).copy()
        spacing = float(np.ravel(st.run["slice_spacing"])[0])
        p0_mc_id = float(st.run["p0c"]) / M_ELECTRON

    # Bunching per slice, from the frame's own particles, with no conversion between the
    # frame and the row: the frame holds the coordinates the tracker held, so a reduction
    # over it is the row's reduction.
    worst_b = 0.0
    for i, f in enumerate(idf):
        sl = read_slices(f, wavelength=1e-10, spacing=spacing)
        for isl, s in enumerate(sl):
            if s["n"] == 0:
                continue
            w = np.asarray(s["weight"])
            th = np.asarray(s["theta"])
            b = abs((w * np.exp(1j * th)).sum() / w.sum())
            worst_b = max(worst_b, abs(b - b_row[i, isl]))
    check("frame vs row: per-slice bunching from the frame's particles", worst_b, 2e-12,
          note=f"[{len(idf)} frames]")

    # The transverse centroid per slice, the frame's weighted mean against the row's in the
    # chart's own units, x and y in meters and px and py over p0. The frame holds the
    # coordinates the tracker held, inside the undulator included, so nothing is converted
    # and the two agree at the rounding of a weighted mean.
    worst_c = 0.0
    for i, f in enumerate(idf):
        sl = read_slices(f, wavelength=1e-10, spacing=spacing)
        for isl, s in enumerate(sl):
            if s["n"] == 0:
                continue
            w = np.asarray(s["weight"])
            got = np.array([np.average(s["x"], weights=w), np.average(s["px"], weights=w) / p0_mc_id,
                            np.average(s["y"], weights=w), np.average(s["py"], weights=w) / p0_mc_id])
            worst_c = max(worst_c, float(np.max(np.abs(got - c_row[i, isl, :4]))))
    check("frame vs row: per-slice transverse centroid from the frame's particles, no conversion",
          worst_c, 1e-14, note=f"[{len(idf)} frames, x and y in m, px and py over p0]")

    # Power per slice, from the raw field frame.
    Z0 = 1.25663706127e-6 * 299792458.0   # mu_0 * c, Bmad's own constants.
    worst_p = 0.0
    for i, f in enumerate(idw):
        with h5py.File(f) as h:
            it = list(h["data"])[0]
            m = h[f"data/{it}/meshes/electricField"]
            E = m["x"][()]
            dx = float(np.ravel(m.attrs["gridSpacing"])[2])
            dy = float(np.ravel(m.attrs["gridSpacing"])[1])
        e2 = E["r"] ** 2 + E["i"] ** 2 if E.dtype.names else np.abs(E) ** 2
        pw = e2.sum(axis=(1, 2)) * dx * dy / (2 * Z0)
        rel = np.abs(pw - p_row[i]) / np.maximum(p_row[i], 1e-30)
        worst_p = max(worst_p, float(np.max(rel)))
    check("frame vs row: per-slice power from the field frame", worst_p, 1e-10)

    # The reduced projections integrate to the row's power, the identity that F16 of
    # BMAD-STATS-EXT-FEL.md states.
    with h5py.File(wd / "id.stats.h5") as h:
        gx = h["coords/grid_x"][()]
        gy = h["coords/grid_y"][()]
        xy = h["field/x/reduced/xy_intensity"][()]
        sx = h["field/x/reduced/slice_x_intensity"][()]
        pwx = h["field/x/power"][()]
    dxg = float(gx[1] - gx[0]);  dyg = float(gy[1] - gy[0])
    tot = xy.sum(axis=(1, 2)) * dxg * dyg
    rel_xy = np.abs(tot - pwx.sum(axis=1)) / np.maximum(pwx.sum(axis=1), 1e-30)
    check("frame vs row: the transverse projection integrates to the record's power",
          float(np.max(rel_xy)), 1e-10)
    per = sx.sum(axis=2) * dxg * dyg
    rel_sx = np.abs(per - pwx) / np.maximum(pwx, 1e-30)
    check("frame vs row: the slice projection integrates to each slice's power",
          float(np.max(rel_sx)), 1e-10)

    # A slice range cuts both files and says which slices it carried.
    run(exe, wd, "rg", FRAME_NML.format(root="rg",
        extra='  run.dump_at_comb = T\n  run.comb_ds_save = 1.0\n'
              '  run.dump_slice_first = 3\n  run.dump_slice_last = 5\n'), threads="4")
    rgf = sorted(wd.glob("rg-[0-9]*.beam.h5"))
    rgw = sorted(wd.glob("rg-[0-9]*.wf.h5"))
    with h5py.File(rgf[-1]) as h:
        it = list(h["data"])[0]
        npatch = h[f"data/{it}/particles/electron/particlePatches/numParticles"].shape[0]
        s1 = int(np.ravel(h["lucifer/frame"].attrs["sliceFirst"])[0])
        s2 = int(np.ravel(h["lucifer/frame"].attrs["sliceLast"])[0])
    with h5py.File(rgw[-1]) as h:
        it = list(h["data"])[0]
        nsf = h[f"data/{it}/meshes/electricField/x"].shape[0]
        off_cut = np.ravel(h[f"data/{it}/meshes/electricField"].attrs["gridGlobalOffset"])
        dz_cut = float(np.ravel(h[f"data/{it}/meshes/electricField"].attrs["gridSpacing"])[0])
    with h5py.File(wframes[-1]) as h:
        it = list(h["data"])[0]
        off_all = np.ravel(h[f"data/{it}/meshes/electricField"].attrs["gridGlobalOffset"])
    check("range: the beam frame's patch count is the range", abs(npatch - 3), 0.5,
          note=f"[slices {s1} to {s2}, {npatch} patches]")
    check("range: the field frame carries the range's slices", abs(nsf - 3), 0.5)
    # The cut field keeps its place. A mesh that started at zero whatever the range would
    # put the cut field at the head of the window in a standard reader, while the cut
    # beam's timeOffset kept each slice where it was.
    check("range: the cut field's mesh starts at its first slice's z [m]",
          abs(off_cut[0] - (off_all[0] + (s1 - 1) * dz_cut)), 1e-18,
          note=f"[{off_cut[0]:.3e} against {(s1 - 1) * dz_cut:.3e}]")
    check("range: and its transverse origin is the whole window's [m]",
          float(np.max(np.abs(off_cut[1:] - off_all[1:]))), 1e-18)

    # The range is refused where it leaves the window.
    (wd / "rgbad.toml").write_text(to_groups(FRAME_NML.format(root="rgbad",
        extra='  run.dump_slice_first = 4\n  run.dump_slice_last = 99\n')))
    r4 = subprocess.run([str(exe), "rgbad.toml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    check("refusal: a slice range outside the window refused (0 = yes)",
          0.0 if (r4.returncode != 0 and "INSIDE THE WINDOW" in (r4.stdout + r4.stderr)) else 1.0, 0.5)

    interludes(exe, wd)
    refposition_checks(exe, pathlib.Path(args.genesis).resolve(), wd)
    grid_parity(exe, wd, pathlib.Path(args.examples).resolve())

    if FAILED:
        print("diagnostic checks: FAIL")
        sys.exit(1)
    print("diagnostic checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
