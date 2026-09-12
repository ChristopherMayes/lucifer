#!/usr/bin/env python3
"""
Diagnostic-output checks (fel-physics.md sec-stats): the stats file, dump-at elements, and the
escaped-field bank, held by cross-identities rather than reference files --

  1. bunch_params reconstruction: bunch_params_from_stats.py applied to the per-record
     sufficient statistics must reproduce the calc_bunch_params values the tracker
     stored at element ends (proves the per-record datasets are sufficient).
  2. banked energy == the ledger's U_escaped column: same slices, independent
     bookkeeping paths (drain-time wavefront_params vs the zero-fill energy sum).
  3. analytic vs numerical propagation: each banked slice's rms size at the exit from
     the free-space ABCD map on its bank-time moment matrix must match the rms of the
     FFT-propagated slice in the pulse file (two independent routes).
  4. pulse pooling: the pooled whole-pulse sigma from per-slice params must match the
     directly computed moments of the concatenated pulse file.
  5. thread invariance: every dataset of stats.h5, escaped and pulse files identical
     at 1 vs 8 threads (dataset-level: HDF5 headers embed creation times).
  6. refusal: a dump_beam_at entry matching no element is refused.
  7. the pre-run: the header's work count, cost estimate and load standing, the footer's
     thin-slice list on a deck loaded under the measured floor, load_only stopping after
     the header, and the refusal of a harmonic the load carries no shot noise at.
  8. the frame series (global%dump_at_comb): one frame per stats record, a frame at an
     element end identical to that element's own dump, labels unique in every frame and
     following migration, the frame's attributes against the lattice, frames landing
     inside an element, and a run with frames identical to the same run without.
  9. a frame against the stats row it sits on: per-slice bunching from the frame's own
     particles, per-slice power from the field frame, and the reduced projections
     integrating to the row's power. Plus a slice range cutting both files.
 10. the interlude's own steps (global%interlude_ds_step): the knob at zero reproduces
     the run that never named it byte for byte, cutting an element reproduces tracking
     it whole, rows and frames land inside a quadrupole and a pipe at the s the rows
     carry, the stats file is exact-sized against the rows written, the header states
     the pieces before it tracks, a wake-carrying element is not cut, and the
     transcribed Genesis interlude refuses the knob.

Run by the benchmark harness; exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

from nml import to_groups

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bunch_params_from_stats import bunch_params_at, pool_wavefront, M_ELECTRON  # noqa: E402
from read_stats import read_stats, same_data  # noqa: E402
from beamio import read_slices  # noqa: E402

import validate_bmad_stats  # noqa: E402  The standard's own conformance checker.

FAILED = False

NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
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

# The pre-run deck. Averaged, seeded so it reaches saturation in one segment and costs
# little, and loaded at 512 macroparticles in beamlets of 8. That is 64 beamlets, half
# the floor the source filter's convergence was measured at, so the load standing and
# the thin-slice list both have something to report.

# The frame-series deck. A short segment, a coarse comb so the frame count stays small,
# and migration on so labels have to follow particles between slices.

FRAME_NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 256
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
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

PRE_NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
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

# Zero-length wake elements, both polarities, for check 7. A zero-length element carrying
# an sr wake is a standard Bmad idiom and legitimate: a wake assigned over an element
# range or a class lands on zero-length members as a matter of course. What differs is
# whether the kick can act. scale_with_length = F means it does, and the element must be
# tracked, take a record, and repeat coords/s at the plane the element before it ended
# on. T means the kick is identically zero (wake_mod scales by l$), so the element is
# skipped and the file must come out as though the wake were not there at all. WKT_NONE
# replaces the T pipes with plain ones to make that comparison.
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


def h5_identical(fa, fb):
    """Whether two HDF5 files hold the same data, meta/ excluded.

    meta/ is excluded, deliberately. Provenance moved from attributes to datasets for
    HDF5's 64 kB attribute cap (fel-physics.md sec-meta), and input_echo carries out_root while
    timestamp carries the clock, so any two runs differ there by construction. Nothing in
    meta/ is physics. Before the move the exclusion existed only by the accident of being
    attributes, which is not something a reader of this check could have reasoned about.

    same_data counts NaN as equal to NaN, which the stats file needs: see read_stats.
    """
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names_a, names_b = [], []
        keep = lambda n, o: isinstance(o, h5py.Dataset) and not n.startswith("meta/")
        a.visititems(lambda n, o: names_a.append(n) if keep(n, o) else None)
        b.visititems(lambda n, o: names_b.append(n) if keep(n, o) else None)
        if sorted(names_a) != sorted(names_b):
            return False
        return all(same_data(a[n][()], b[n][()]) for n in names_a)


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-1000:]}")
        sys.exit(1)
    return r



# The interlude's own steps (global%interlude_ds_step). Tao's comb has two halves and
# comb_ds_save is the first: save_a_bunch_step selects among positions a tracker already
# reaches. This is the second, tao_lattice_calc_mod's n_slice branch, which cuts an
# element that would otherwise be crossed in one map. The lattice is a wiggler, a
# quadrupole whose fringes are live, a pipe and a wiggler, so a cut element of each kind
# is exercised and the wigglers keep their own steps either way.

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

INT_NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 256
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
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
    stats rows and frames along an interlude, and the guarantees the knob carries.

    The knob is off by default and off has to be the path every run took before it
    existed, so the first check is byte identity against a deck that never names it.
    The second is that cutting an element reproduces tracking it whole: the pieces come
    from Bmad's element_slice_iterator, which is what keeps an element's fringes at its
    own ends rather than giving every piece a pair. That construction is why the
    iterator is used; on FEL beam sizes the difference it buys is below notice, and the
    hand-split comparison here records that rather than asserting it.
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
        extra="  global%interlude_ds_step = 0\n"))
    same = (wd / "int_w.diag.txt").read_bytes() == (wd / "int_off.diag.txt").read_bytes()
    check("interlude: the knob at zero is the run that never named it (0 = yes)",
          0.0 if same else 1.0, 0.5)

    run(exe, wd, "int_p", INT_NML.format(lat="intlat.bmad", root="int_p",
        extra="  global%interlude_ds_step = 0.04\n"))
    run(exe, wd, "int_h", INT_NML.format(lat="intlath.bmad", root="int_h", extra=""))
    w = np.loadtxt(wd / "int_w.diag.txt")
    p = np.loadtxt(wd / "int_p.diag.txt")
    h = np.loadtxt(wd / "int_h.diag.txt")

    # The pieces are tracking, not bookkeeping, so what they must not do is change the
    # answer: the element's exit is the element's exit however many pieces crossed it.
    check("interlude: cut into pieces, exit power against the whole element",
          abs(p[-1, 2] - w[-1, 2]) / w[-1, 2], 1e-9)
    check("interlude: cut into pieces, exit sigma_x against the whole element",
          abs(p[-1, 8] - w[-1, 8]) / w[-1, 8], 1e-9)
    check("interlude: the hand-split lattice agrees too, so the fringe difference the "
          "iterator buys is below notice here",
          abs(h[-1, 2] - w[-1, 2]) / w[-1, 2], 1e-9)

    # Rows land inside the quadrupole and inside the pipe, where none did before. The
    # element spans are 0.3 to 0.5 and 0.5 to 0.9 on this line.
    def inside(rows, z0, z1):
        z = rows[:, 0]
        return int(((z > z0 + 1e-9) & (z < z1 - 1e-9)).sum())
    check("interlude: rows strictly inside the quadrupole, whole (0 expected)",
          abs(inside(w, 0.3, 0.5) - 0), 0.5)
    check("interlude: rows strictly inside the quadrupole, cut (4 expected)",
          abs(inside(p, 0.3, 0.5) - 4), 0.5, note=f"[{inside(p, 0.3, 0.5)}]")
    check("interlude: rows strictly inside the pipe, cut (9 expected)",
          abs(inside(p, 0.5, 0.9) - 9), 0.5, note=f"[{inside(p, 0.5, 0.9)}]")

    # The record-count precompute replays the walk, so the stats file is exact-sized:
    # a mismatch would mean the two disagree about where a row falls.
    with h5py.File(wd / "int_p.stats.h5") as f:
        nrec = f["coords/s"].shape[0]
        step = float(np.ravel(f["params/global/interlude_ds_step"])[0])
    check("interlude: the stats file holds exactly the rows the diag carries",
          abs(nrec - len(p)), 0.5, note=f"[{nrec} against {len(p)}]")
    check("interlude: the stats file records the step it ran at", abs(step - 0.04), 1e-12)

    # Frames along the interlude, which is what the knob is for.
    run(exe, wd, "int_f", INT_NML.format(lat="intlat.bmad", root="int_f",
        extra="  global%interlude_ds_step = 0.04\n  global%dump_at_comb = T\n"))
    frames = sorted(wd.glob("int_f-[0-9]*.beam.h5"))
    check("interlude: one frame per stats row with the pieces on",
          abs(len(frames) - len(p)), 0.5, note=f"[{len(frames)} frames, {len(p)} rows]")
    zs = []
    for fr in frames:
        with h5py.File(fr) as f:
            zs.append(float(np.ravel(f.attrs["sPosition"])[0]))
    zs = np.array(zs)
    check("interlude: the frames carry the rows' own s", np.abs(zs - p[:, 0]).max(), 1e-9)
    nin = int(((zs > 0.5 + 1e-9) & (zs < 0.9 - 1e-9)).sum())
    check("interlude: frames strictly inside the pipe (9 expected)", abs(nin - 9), 0.5,
          note=f"[{nin}]")

    # The pre-run header states the pieces before it tracks, and says nothing when the
    # element is tracked whole.
    r = run(exe, wd, "int_hdr", INT_NML.format(lat="intlat.bmad", root="int_hdr",
        extra="  global%interlude_ds_step = 0.04\n  global%load_only = T\n"))
    m = re.search(r"Interlude   steps of .*: (\d+) pieces, (\d+) element", r.stdout)
    check("interlude: the header states the pieces before it tracks",
          0.0 if (m and int(m.group(1)) == 15 and int(m.group(2)) == 2) else 1.0, 0.5,
          note=f"[{m.group(0) if m else 'absent'}]")
    r = run(exe, wd, "int_hdr0", INT_NML.format(lat="intlat.bmad", root="int_hdr0",
        extra="  global%load_only = T\n"))
    check("interlude: and says nothing where nothing is cut",
          0.0 if "Interlude" not in r.stdout else 1.0, 0.5)

    # Refused. The element short-range wake is Bmad's own once-per-passage kick of the
    # element's length, so a cut element would take one kick a piece: the walk declines
    # to cut there, as Tao declines where CSR or space charge would see the cut. The
    # transcribed Genesis interlude is validated by transcription fidelity and its step
    # is the element, so the knob is refused outright with that model.
    # A wake whose table is zero everywhere: the element carries one, which is what the
    # walk keys on, and the physics is the no-wake run's, so the rows are the only thing
    # that may differ.
    (wd / "intlatw.bmad").write_text(INT_LAT.replace(
        "P1: pipe, l = 0.40",
        "P1: pipe, l = 0.40, sr_wake = {scale_with_length = T, z_long = "
        "{position_dependence = none, w = {-1e-6 0.0, 0.0 0.0, 1e-6 0.0,}}}"))
    run(exe, wd, "int_wk", INT_NML.format(lat="intlatw.bmad", root="int_wk",
        extra="  global%interlude_ds_step = 0.04\n"))
    wk = np.loadtxt(wd / "int_wk.diag.txt")
    check("interlude: a wake-carrying element is not cut", abs(inside(wk, 0.5, 0.9) - 0), 0.5,
          note=f"[{inside(wk, 0.5, 0.9)} rows inside it]")
    check("interlude: and the elements without a wake still are",
          abs(inside(wk, 0.3, 0.5) - 4), 0.5, note=f"[{inside(wk, 0.3, 0.5)}]")
    check("interlude: the zero wake left the physics alone, so only the rows differ",
          abs(wk[-1, 2] - p[-1, 2]) / p[-1, 2], 1e-9)

    (wd / "int_gen.nml").write_text(to_groups(INT_NML.format(lat="intlat.bmad", root="int_gen",
        extra='  global%interlude_ds_step = 0.04\n  global%interlude_model = "genesis"\n')))
    rg = subprocess.run([str(exe), "int_gen.nml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = rg.returncode != 0 and "TRANSCRIPTION OF GENESIS" in (rg.stdout + rg.stderr)
    check("interlude: refused with the transcribed Genesis interlude (0 = yes)",
          0.0 if refused else 1.0, 0.5)

    # The switch off. With bmad_com%sr_wakes_on off the walk applies no wake and cuts the
    # element like any other, and the setup's record-count precompute must make the same
    # decision: it once counted the element whole while the walk cut it, and the run
    # stopped at the stats writer with more records than the precomputed count.
    wo = run(exe, wd, "int_wko", INT_NML.format(lat="intlatw.bmad", root="int_wko",
        extra="  global%interlude_ds_step = 0.04\n  bmad_com%sr_wakes_on = F\n"))
    wko = np.loadtxt(wd / "int_wko.diag.txt")
    check("interlude: with sr_wakes_on off the wake element is cut like any other",
          abs(inside(wko, 0.5, 0.9) - 9), 0.5, note=f"[{inside(wko, 0.5, 0.9)} rows inside it]")
    check("interlude: and its physics is the whole-element run's",
          abs(wko[-1, 2] - p[-1, 2]) / p[-1, 2], 1e-9)

    # Space charge on a Bmad-seam interlude is refused. The seam tracks one slice's bunch
    # at a time through track1_bunch, whose own space-charge and CSR paths want the centroid
    # orbit of the whole beam, and without it Bmad printed that it must be supplied and
    # returned the bunch untracked with no error set: a 0.1 m drift carrying the attribute
    # lost its transport and the run reported success. The deck configures a term so the
    # refusal reached is this one and not the one for a solver with nothing to solve.
    (wd / "intlat_sc.bmad").write_text(INT_LAT.replace(
        "P1: pipe, l = 0.40", "P1: pipe, l = 0.40, space_charge_method = slice"))
    (wd / "int_sc.nml").write_text(to_groups(INT_NML.format(lat="intlat_sc.bmad", root="int_sc",
        extra="  bmad_com%csr_and_space_charge_on = T\n  space_charge%nz = 1\n"
              "  space_charge%rmax = 250e-6\n")))
    rs = subprocess.run([str(exe), "int_sc.nml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = rs.returncode != 0 and "NOT AN FEL ELEMENT" in (rs.stdout + rs.stderr)
    check("interlude: space charge on a seam interlude is refused (0 = yes)",
          0.0 if refused else 1.0, 0.5)
    check("interlude: and Bmad's own centroid message never printed (0 = yes)",
          0.0 if "CENTROID MUST BE SUPPLIED" not in (rs.stdout + rs.stderr) else 1.0, 0.5)

    # A zero-length element that can act is tracked. The walk skips a zero-length element
    # Bmad tracks as the identity, and once skipped every zero-length element without a
    # wake, a thin kicker among them, whose kick then vanished with no message. The same
    # kick through a kicker of finite length is the control, and a zero-length pipe in the
    # kicker's place says how large the kick is against the beam's own centroid and is
    # the one element of the three the walk still skips, so its stats file is one record
    # shorter: a kicker is tracked whatever its kick, since its key carries a map.
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

    # A frame's particle records say where the frame is. The beam writer initialized the
    # coords at the containing element's upstream face whatever the frame's position, so
    # a frame inside a piece carried the element's start in sPosition and its entry time
    # in timeOffset while the file's own attributes said the frame's position.
    worst_s, worst_t = 0.0, 0.0
    for fr in frames:
        with h5py.File(fr) as f:
            it = f["data"][sorted(f["data"])[0]]
            el = it["particles/electron"]
            def value(rec):
                return float(np.ravel(rec.attrs["value"] if "value" in rec.attrs else rec[()])[0])
            worst_s = max(worst_s, abs(value(el["sPosition"]) - float(np.ravel(f.attrs["sPosition"])[0])))
            worst_t = max(worst_t, abs(value(el["timeOffset"]) - float(np.ravel(it.attrs["time"])[0])))
    check("frames: the particle sPosition is the frame's own [m]", worst_s, 1e-12)
    check("frames: the particle timeOffset is the iteration's time [s]", worst_t, 1e-18)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--latdir", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    latdir = pathlib.Path(args.latdir).resolve()
    wd = pathlib.Path(args.workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "aramis_1seg.bmad").write_bytes((latdir / "aramis_1seg.bmad").read_bytes())
    (wd / "dg_wrap.bmad").write_text(WRAP)

    run(exe, wd, "dg", NML.format(lat="dg_wrap.bmad", root="dg", extra=""))

    # 0. The acceptance test is the standard's own validator, which knows no program:
    # a bmad-stats file conforms when it reports zero MUST failures, and everything it
    # checks is checked from the spec alone. The named checks below hold what it cannot
    # see: the physics identities and the declarations only this writer knows it owes.
    rep = validate_bmad_stats.validate(wd / "dg.stats.h5")
    for level, where, msg in rep.rows:
        if level == "MUST":
            print(f"    {level} {where}: {msg}")
    check("bmad-stats: the generic validator reports zero MUST failures", rep.n_must, 0.5,
          note=f"[{sum(1 for r in rep.rows if r[0] != 'INFO')} findings in "
               f"{len(rep.rows)} rows]")

    with read_stats(wd / "dg.stats.h5") as st:
        nslice = int(st.run["n_slice"])

        # 0a. The file describes itself: a walk of every dataset, checking that it says
        # what it is and what its dimensions run over. This is what makes a generic
        # reader possible, so it is checked generically rather than name by name.
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
        # dataset from @axes alone. It fails if a name does not resolve to a coordinate,
        # if a dimension has no name, if a length disagrees with its coordinate, or if
        # one dataset names an axis twice, which is what would force a reader to guess
        # from a length or to dedupe by rule. Nothing here knows a dataset name.
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
        # non-decreasing, and where it repeats the two records must sit in different
        # elements, which is what a zero-length element applying a wake kick does: it
        # records at the plane the element before it ended on. Two records at one plane
        # inside one element would be a defect in the walk that this axis choice hides.
        # The zero-length wake lattice below is what makes this fire on real duplicates.
        ds = np.diff(st.s)
        dup = np.flatnonzero(ds == 0)
        unexplained = [int(i) for i in dup if st.ix_ele[i] == st.ix_ele[i + 1]]
        check("stats: s is non-decreasing along the record axis", int(np.sum(ds < 0)), 0.5)
        check("stats: every repeated s straddles an element boundary", len(unexplained), 0.5,
              note=f"[{len(dup)} repeats in {len(st.s)} records]")

        # 0a4. Provenance is data, not an attribute. HDF5 caps one attribute at 64 kB
        # (measured: 65495 bytes is the largest that writes), and the echoed namelist is
        # already 12 kB with a real lattice text at 37 kB, so meta/ as attributes was at
        # half the cap with a warning for a failure path: the provenance would have gone
        # missing silently. Every text in meta/ is a dataset now, and nothing anywhere in
        # the file leans on the cap.
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

        # 0a5. And it does not overstate itself. dg tracks a wrapper lattice, so the
        # parser opened two files and lattice_source holds only the outer one. Recording
        # the count is what keeps that honest, and this run is the case that proves it:
        # a reader of lattice_source alone would think it had the lattice.
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

        # 0a7. The slice coordinates are exact, and one of them is not. The grid is
        # uniform in time: the migration invariant carries each particle's own beta, and
        # z = -beta*c*(t-t_ref), so at a grid point the beta cancels and the arrival-time
        # separation is ct_slice/c with no beta anywhere. t_slice must therefore be
        # -ct_slice/c exactly, which it was not until 2.3 (it carried a spurious beta0,
        # 3.9e-9). z_slice is the one that needs a reference, so it is only equal to
        # beta0*ct_slice to the rounding of two ways of forming beta0.
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

        # 0a8. The envelope data are the particles'. rel_max/rel_min are order
        # statistics relative to the stored centroid, and the dump at the UND end holds
        # the same particles the record saw, so the position entries must match to the
        # bit: the file stores x verbatim, and IEEE subtraction of the same centroid is
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

        # 1. Reconstruction at the element end (= the last record). The entry is an index
        # into the group's own axis, looked up rather than assumed: that lookup is the
        # only thing standing between a reader and a silent transposition. The projected
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
        # Measured 4.6e-8 (numpy eig vs mat_eigen). The projected planes agree exactly.
        check("stats: bunch_params reconstruction vs stored calc_bunch_params", worst, 1e-6)

        # 1b. The whole-window row is assembled from the per-slice moments by the
        # pooled-covariance identity rather than summed over particles. Re-implement the
        # identity here, independently, including the local-to-global z map: a shift of
        # beta*(is-1)*spacing plus the pz shear that map carries because beta is the
        # particle's own. Measured against the particle sum when this landed: 4.0e-12 on
        # this config, 5.0e-11 on the 96-slice SASE example.
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
        # The check has teeth: dropping the between-group term (m_s - m)(m_s - m)^T --
        # the term that makes the identity exact -- must be caught, not absorbed. It
        # collapses sigma(5,5) from the window length squared to a slice length squared.
        check("stats: the between-group term is load-bearing (0 = confirmed)",
              0.0 if worst_nobg > 1e-3 else 1.0, 0.5,
              note=f"[dropping it moves the row by {worst_nobg:.2e}]")

        # 1c. The moments the element-end twiss rests on are the record's, not a second
        # copy: beam/slice_twiss/centroid at an element end must be beam/slice/centroid
        # at the record the mask selects. Exactly, since one is written from the other.
        worst_ef = float(np.max(np.abs(st["beam/slice/centroid"][st.at_end] -
                                      st["beam/slice_twiss/centroid"])))
        worst_ef = max(worst_ef, float(np.max(np.abs(st["beam/slice/sigma"][st.at_end] -
                                                    st["beam/slice_twiss/sigma"]))))
        check("stats: the element-end moments ARE the record's (abs)", worst_ef, 1e-30)

        # And the field's theta moments are computed where the file says they are:
        # (nz, ns) per slice, at every element end that has field, plus the initial
        # record (the walk takes that one with angles so the starting state is
        # complete). An empty slice has no moments to compute and says so.
        valid = st["field/x/angle_moments_valid"].astype(bool)
        has_field = st["field/x/power"] > 0
        rows = np.zeros(len(st.s), bool)
        rows[0] = True
        rows |= st.at_end
        want = has_field & rows[:, None]
        check("stats: angle moments valid exactly where they were computed",
              int(np.sum(valid != want)), 0.5,
              note=f"[{int(valid.sum())} of {valid.size} slice-records]")

        # NaN, not a zero that reads as an answer, everywhere else. This is what lets a
        # consumer tell an empty slice's missing moments from a real measurement.
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
    # Pulse index nlive+j holds banked slice nb-j+1: reverse into transmission order.
    sxx_num = np.array(sxx_num)[::-1]
    syy_num = np.array(syy_num)[::-1]
    w_num = np.array(w_num)[::-1]

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
    # Measured max 8.4e-3 (grid-edge wrap of the FFT route on divergent noise slices).
    check("bank: analytic vs FFT-propagated rms at exit, max rel", rel.max(), 2e-2)

    # 4. Pulse pooling: pooled per-slice params vs directly computed pulse moments.
    #    Banked side: params propagated analytically to the exit. Live side: the last
    #    stats record (angles valid there). Compare sigma_xx of the pooled result
    #    against pooling the pulse file's numerically computed per-slice moments.
    cen_b = np.stack([pms[:, 0] + dz * pms[:, 1], pms[:, 1],
                      pms[:, 2] + dz * pms[:, 3], pms[:, 3]], axis=1)
    sig_b = np.zeros((nb, 4, 4))
    sig_b[:, 0, 0] = sxx_ana;  sig_b[:, 2, 2] = syy_ana
    en_b = pms[:, 20]
    cen_l = f_cen;  sig_l = f_sig;  en_l = f_en
    _, pooled_ana, _ = pool_wavefront(np.vstack([cen_l, cen_b]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_b.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_b]))
    # Numerical pooling over the same banked slices from the pulse file:
    cen_n = np.stack([np.array(cx_num)[::-1], np.zeros(nb), np.array(cy_num)[::-1], np.zeros(nb)], axis=1)
    sig_n = np.zeros((nb, 4, 4))
    sig_n[:, 0, 0] = sxx_num;  sig_n[:, 2, 2] = syy_num
    en_n = w_num * (en_b.sum() / max(w_num.sum(), 1e-300))   # Same total, per-slice proportional.
    _, pooled_num, _ = pool_wavefront(np.vstack([cen_l, cen_n]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_n.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_n]))
    check("pulse: pooled sigma_xx analytic vs numerical routes, rel",
          abs(pooled_ana[0, 0] / pooled_num[0, 0] - 1), 2e-2)

    # 5. Thread invariance: every dataset of all three files identical at 1 vs 8
    #    threads. Dataset-level, not raw bytes: HDF5 object headers embed creation
    #    times, so two files with bit-identical data differ as byte streams -- the
    #    data is the invariance claim, the container metadata is not.
    run(exe, wd, "dg1", NML.format(lat="dg_wrap.bmad", root="dg1", extra=""), threads="1")

    same = all(h5_identical(wd / f"dg{s}", wd / f"dg1{s}")
               for s in (".stats.h5", "-escaped.fld.h5", "-pulse.fld.h5"))
    check("thread invariance: stats/escaped/pulse data identical 1 vs 8 (1 = yes)",
          0.0 if same else 1.0, 0.5)

    # 5a. Zero-length wake elements, both polarities. This is the lattice the repeated-s
    # check was written for, and until it existed that check had only ever seen zero
    # repeats, which is to say it was untested. WKF can kick, so it is tracked and its
    # record lands on the plane UND ended at, repeating coords/s. The two WKT pipes
    # cannot kick, so they are skipped and the file must be identical to one whose
    # zero-length pipes carry no wake at all.
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
        extra="  global%record_environment = T\n"), threads="2")
    with read_stats(wd / "dgenv.stats.h5") as st:
        got = sorted(k for k in st.meta if k in ("user", "cwd"))
    check("provenance: global%record_environment restores user and cwd",
          abs(len(got) - 2), 0.5, note=f"[{', '.join(got) if got else 'neither'}]")

    # 6. Refusal: dump list entry matching nothing.
    (wd / "dg_bad.nml").write_text(to_groups(NML.format(lat="dg_wrap.bmad", root="dg_bad",
                                             extra='  dump_beam_at(2) = "NO_SUCH_ELEMENT"\n')))
    r = subprocess.run([str(exe), "dg_bad.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = r.returncode != 0 and "NO_SUCH_ELEMENT" in (r.stdout + r.stderr)
    check("refusal: unknown dump_beam_at element refused (1 = yes)",
          0.0 if refused else 1.0, 0.5)

    # 7. The pre-run (doc/user-guide.md). The deck loads 512 macroparticles per slice in
    # beamlets of 8, which is 64 beamlets and half the measured floor, so the load
    # standing and the thin-slice list both have something to report. The count is
    # arithmetic and is checked as arithmetic; the estimate is one machine's and is
    # checked only for being present and positive.

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

    (wd / "pre_h.nml").write_text(to_groups(PRE_NML.format(root="pre_h",
                                  extra="  harmonics = 1, 5\n")))
    r3 = subprocess.run([str(exe), "pre_h.nml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    said = r3.stdout + r3.stderr
    refused = r3.returncode != 0 and "BEAMLET_SIZE TO 11" in said
    check("refusal: harmonic 5 on beamlets of 8 refused, naming the size that carries it",
          0.0 if refused else 1.0, 0.5)

    # 8. The frame series (doc/reading-output.md). The deck migrates, so labels have to
    # survive particles changing slice, and its comb is coarse enough to keep the frame
    # count small. The element-end dump is asked for by name so the two writers can be
    # compared on one position.

    fr_extra = ('  global%dump_at_comb = T\n  global%comb_ds_save = 1.0\n'
                '  global%migrate = T\n  global%dump_beam_at = "UND"\n')
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
            got = h.attrs["elementName"]
            got = got.decode() if isinstance(got, bytes) else str(got)
            spos = float(np.ravel(h.attrs["sPosition"])[0])
        if got.strip() != names[i] or abs(spos - s_rec[i]) > 1e-9:
            bad_attr += 1
    check("frames: each frame names the element and s the stats row carries", bad_attr, 0.5)

    # An FEL frame carries the undulator a reader rebuilds the wiggle from.
    with h5py.File(frames[-1]) as h:
        has_und = all(k in h.attrs for k in ("aw", "ku", "helical", "tilt", "felMethod"))
    check("frames: an FEL frame carries aw, ku and the method (0 = yes)",
          0.0 if has_und else 1.0, 0.5)

    # A run that writes frames must be the run that does not. The instrument may not
    # steer what it observes, which is why the field's rotation is put back.
    run(exe, wd, "nofr", FRAME_NML.format(root="nofr",
        extra='  global%comb_ds_save = 1.0\n  global%migrate = T\n'), threads="4")
    check("frames: the run is identical to the same run without them (0 = yes)",
          0.0 if h5_identical(wd / "fr.stats.h5", wd / "nofr.stats.h5") else 1.0, 0.5)

    # The unaveraged chart, written and named. A frame inside an unaveraged segment is
    # taken where px still carries the undulator quiver, and the writer is allowed that
    # chart because the frame names it. The mode wrote no frames at all before: the chart
    # assertion sat inside the slice-to-bunch conversion, so the whole segment interior
    # refused, and dump_at_comb = F was the only way an unaveraged run could finish.
    #
    # The quiver is a common offset at a position along the undulator and not a spread,
    # so it is the mean of px that carries it and not the rms. Its size is aw/gamma of
    # the reference momentum, which on this deck is about 4e5 eV/c, where the averaged
    # chart holds a few hundred. Asserting the swing is what proves the file holds the
    # chart it claims rather than a quietly averaged one.

    (wd / "dg_unavg.bmad").write_text(WRAP)
    ua_extra = ('  global%dump_at_comb = T\n  global%comb_ds_save = 1.0\n'
                '  global%device = "off"\n')
    # Run it without run()'s own exit, so that a mode which cannot write a frame is a
    # named failing check here rather than an abort with the reason buried in a log.
    ua_text = (FRAME_NML.format(root="ua", extra=ua_extra)
               .replace('lat_file = "aramis_1seg.bmad"', 'lat_file = "dg_unavg.bmad"'))
    (wd / "ua.nml").write_text(to_groups(ua_text))
    r = subprocess.run([str(exe), "ua.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    ua_frames = sorted(wd.glob("ua-[0-9]*.beam.h5"))
    check("frames: an unaveraged run writes them at all (exit 0)", float(r.returncode), 0.5,
          note="" if r.returncode == 0 else f"[{r.stdout.strip().splitlines()[-1][:70]}]")
    check("frames: the unaveraged run wrote a series",
          0.0 if len(ua_frames) > 1 else 1.0, 0.5, note=f"[{len(ua_frames)} frames]")

    if len(ua_frames) > 1:
        def px_mean(f):
            with h5py.File(f) as h:
                m = h.attrs.get("felMethod", b"")
                m = m.decode() if isinstance(m, bytes) else str(m)
                g = h["data"]; it = g[list(g.keys())[0]]
                pg = it["particles"]; b0 = pg[list(pg.keys())[0]]
                return m, float(np.asarray(b0["momentum"]["x"][()]).mean())
        named = [px_mean(f) for f in ua_frames[1:]]
        avg = [px_mean(f) for f in frames[1:]]
        bad_name = sum(1 for m, _ in named if m != "Unaveraged")
        check("frames: every unaveraged frame names its chart", bad_name, 0.5)
        swing_u = max(abs(v) for _, v in named)
        swing_a = max(abs(v) for _, v in avg)
        check("frames: the unaveraged chart carries the quiver in the mean of px",
              0.0 if swing_u > 20 * swing_a else 1.0, 0.5,
              note=f"[unaveraged {swing_u:.2e}, averaged {swing_a:.2e} eV/c]")

    # Restarting from a frame needs no check of its own. A frame at an element end is
    # dataset-identical to that element's dump, which the check above measures, and
    # check_program holds the windowed composition for those dumps. A frame taken inside
    # an element is not a restart point at all, since global%track_start names elements
    # and not positions.
    #
    # What does need measuring is that frames land inside an element, which is the whole
    # reason for the series: a comb that only fired at element ends would pass every
    # check above and give a movie with one frame per undulator.
    inside = 0
    for i, f in enumerate(frames):
        with h5py.File(f) as h:
            spos = float(np.ravel(h.attrs["sPosition"])[0])
            nm = h.attrs["elementName"]
            nm = (nm.decode() if isinstance(nm, bytes) else str(nm)).strip()
        if nm == "UND" and 0.0 < spos < 3.99:
            inside += 1
    check("frames: the series lands inside the undulator, not only at its ends",
          0.0 if inside >= 2 else 1.0, 0.5, note=f"[{inside} of {len(frames)} inside UND]")

    # 9. A frame against its stats row. These three identities pin the phase convention,
    # the time-order rotation and the slice indexing end to end: a picture drawn from a
    # frame and a number read from the row have to be the same run.

    rr_extra = ('  global%dump_at_comb = T\n  global%comb_ds_save = 1.0\n'
                '  global%dump_reduced = T\n')
    run(exe, wd, "id", FRAME_NML.format(root="id", extra=rr_extra), threads="4")
    idf = sorted(wd.glob("id-[0-9]*.beam.h5"))
    idw = sorted(wd.glob("id-[0-9]*.wf.h5"))

    with read_stats(wd / "id.stats.h5") as st:
        b_row = np.asarray(st["beam/slice/bunching"]).copy()
        p_row = np.asarray(st["field/total/power"]).copy()
        spacing = float(np.ravel(st.run["slice_spacing"])[0])

    # Bunching per slice, from the frame's own particles. read_slices does the documented
    # reconstruction, theta = -ks c t, which is where the reference phase the writer
    # folded into the time coordinate comes back.
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

    # The reduced projections integrate to the row's power, which is F16's own identity.
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
        extra='  global%dump_at_comb = T\n  global%comb_ds_save = 1.0\n'
              '  global%dump_slice_first = 3\n  global%dump_slice_last = 5\n'), threads="4")
    rgf = sorted(wd.glob("rg-[0-9]*.beam.h5"))
    rgw = sorted(wd.glob("rg-[0-9]*.wf.h5"))
    with h5py.File(rgf[-1]) as h:
        it = list(h["data"])[0]
        npatch = h[f"data/{it}/particles/electron/particlePatches/numParticles"].shape[0]
        s1 = int(np.ravel(h.attrs["sliceFirst"])[0])
        s2 = int(np.ravel(h.attrs["sliceLast"])[0])
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
    # The cut field keeps its place. Its mesh once started at zero whatever the range, so a
    # standard reader put the cut field at the head of the window while the cut beam's
    # timeOffset kept each slice where it was.
    check("range: the cut field's mesh starts at its first slice's z [m]",
          abs(off_cut[0] - (off_all[0] + (s1 - 1) * dz_cut)), 1e-18,
          note=f"[{off_cut[0]:.3e} against {(s1 - 1) * dz_cut:.3e}]")
    check("range: and its transverse origin is the whole window's [m]",
          float(np.max(np.abs(off_cut[1:] - off_all[1:]))), 1e-18)

    # The range is refused where it leaves the window.
    (wd / "rgbad.nml").write_text(to_groups(FRAME_NML.format(root="rgbad",
        extra='  global%dump_slice_first = 4\n  global%dump_slice_last = 99\n')))
    r4 = subprocess.run([str(exe), "rgbad.nml"], cwd=wd, capture_output=True, text=True,
                        env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    check("refusal: a slice range outside the window refused (1 = yes)",
          0.0 if (r4.returncode != 0 and "INSIDE THE WINDOW" in (r4.stdout + r4.stderr)) else 1.0, 0.5)

    interludes(exe, wd)

    if FAILED:
        print("diagnostic checks: FAIL")
        sys.exit(1)
    print("diagnostic checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
