#!/usr/bin/env python3
"""
Statistical check of the weighted shot-noise loader against the closed form of its kick.

Genesis cannot represent nonuniform weights, so there is no reference to compare to.

Physics under test: after quiet loading plus Fawley-style noise, a beamlet standing for
N real electrons has <|b(h)|^2> = R_h(N)/N at every imposed harmonic, where R_h is the
ratio shot_noise_ratio.py computes and docs/shot-noise.md derives. R_h is 1 - 2/N to
first order and goes to 1 as N grows. A slice's <|b(h)|^2> * N_lambda is then the
charge-weighted mean of R_h over its beamlets, N_lambda the slice's real electron count.
This holds for uniform and for nonuniform per-particle weights.

Method: run lucifer with load_only = T over many seeds, in both weight modes. The modes
are uniform, and gen_test_weights = T, which alternates beamlet weights of 0.25x and
1.75x at constant charge. Each .beam.h5 gives the charge-weighted |b(h)|^2 per slice for
the imposed harmonics, and the expected value from the weights it carries, over 256 slices. The scaled
mean m = <|b(h)|^2 * N_lambda> is tested against the mean expected value r. b is a sum
of many independent beamlet contributions, so |b|^2*N_lambda is exponential with mean r
to excellent approximation, and the mean over n samples has sigma = r/sqrt(n). The check
is |m - r| < 5 r/sqrt(n) per weight mode, plus a looser per-harmonic check on n/3 samples.

The bunch charge is set so that a beamlet holds about 146, 50, 18 and 4 electrons, where
the fundamental's ratio is 0.98, 0.95, 0.87 and 0.56. A last charge leaves a quarter of an
electron to a beamlet. The loader then draws the kick for one electron and warns, and
the ratio is a quarter of the ratio at one electron, 0.040 at the fundamental. These
loads set beam.beamlet_min_electrons = 0, which gives every slice the same beamlets, so
that the count a beamlet holds is the one stated.

A Gaussian bunch is then loaded with the key at its default and at 18 electrons. Its
peak slice holds 50 electrons a beamlet and its last slices two electrons in all, so
the slices away from the peak are loaded with fewer beamlets. Each such slice must
carry the charge it carries with the key at 0, must have the beamlets
min(n_particle/beamlet_size, max(1, floor(N/n_min))) for its N electrons, and its
statistic must agree with the closed form for the beamlets it has. The mean is held
within five standard deviations of the mean, taken from the samples, since a slice of
one beamlet is not a sum of many. A uniform load at 146 electrons a beamlet, where no
slice has fewer, must be the same bits with the key at its default as with the key at 0.

The same statistic then runs on the loads that do not draw their own beamlets
(fel-physics.md, sec-noise):
- load_mode = "keep" from a generated bunch, where the loader makes the beamlets.
- A bunch of copies made outside the loader and loaded without the quiet start. The
  noise goes on the particles sharing their transverse coordinates.
- The same copies with their coordinates jittered so none coincide. The noise goes on
  the occupants of each deposit cell.
The loader's message must name the grouping it used.

Last, the stats file's mode split, the power within split_angle of the axis, is checked
against the same quantity taken from a field dump to 1e-6.

Usage: check_shot_noise.py --exe <lucifer> --workdir <dir> [--seeds N] [--lat <bmad file>]
Exit 0 only if every check passes.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fieldio  # noqa: E402

import lucifer.beam as beamio
import bunchfile
import shot_noise_ratio
from checkrec import record
from lucifer.flat import to_toml
to_groups = to_toml

E_CHARGE = 1.602176634e-19

# The loads below share these. The cell rule reads the grid.
BEAMLET_SIZE = 8
GRID_N_PTS = 32
GRID_HALF_WIDTH = 2e-4

# The generated load: a uniform bunch over N_SLICE slices of N_PARTICLE / BEAMLET_SIZE
# beamlets. The slice count sets the samples a seed gives, 256 slices and 15 seeds a
# bound of 4.7% of the closed form over the three harmonics.
N_SLICE = 256
N_PARTICLE = 1024

# The electrons a beamlet holds in each scan of the generated load, the first the level
# every other load here is run at. The last is under one electron.
ELECTRONS = 146.4
SCAN = (ELECTRONS, 50.0, 18.0, 4.0, 0.25)

# The window the input file below states. A dump carries the slice partition and not the
# radiation it was sliced on, so the reader is told.
LAMBDA0 = 1e-10
SPACING = 3 * LAMBDA0

NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = {n_particle}
  beam.bunch_charge = {charge:.12e}
  beam.z_distribution = "GRID"
  beam.z_min = {z_min:.6e}
  beam.z_max = {z_max:.6e}
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = {grid_n_pts}
  grid_half_width = {grid_half_width}
  beamlet_size = {beamlet_size}
  window_length = {window:.6e}
  n_wavelength = 3
  shot_noise = T
  beam.beamlet_min_electrons = 0
  gen_test_weights = {testw}
  ran_seed = {seed}
  load_only = T
&end
"""

# The Gaussian bunch: the window is eight sigmas, the peak slice holds PEAK electrons a
# beamlet at the input file's count, and the last slices 3.4e-4 of the peak's charge.
PEAK = 50.0
GAUSS_NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = {n_particle}
  beam.bunch_charge = {charge:.12e}
  beam.sig_z = {sig_z:.12e}
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = {grid_n_pts}
  grid_half_width = {grid_half_width}
  beamlet_size = {beamlet_size}
  window_length = {window:.6e}
  n_wavelength = 3
  shot_noise = T
{limit}  ran_seed = {seed}
  load_only = T
&end
"""


def cell_groups(sl):
    """The deposit cell of every particle, by the loader's own rule (fel_impose_noise)."""
    dx = 2 * GRID_HALF_WIDTH / (GRID_N_PTS - 1)
    ix = np.clip(((sl["x"] + GRID_HALF_WIDTH) / dx).astype(int), 0, GRID_N_PTS - 1)
    iy = np.clip(((sl["y"] + GRID_HALF_WIDTH) / dx).astype(int), 0, GRID_N_PTS - 1)
    return iy * GRID_N_PTS + ix


def b2_samples(dump_file, harmonics, cells=False):
    """Charge-weighted |b(h)|^2 * N_lambda for every slice and harmonic, each with the
    value the closed form expects of that slice.

    The dump is openPMD, so the weights are the ones the loader wrote. Reading the
    weights tests the loader instead of assuming it. A Genesis-format dump carries one
    current per slice and would need the alternating test pattern reconstructed.

    A particle's group holds BEAMLET_SIZE times its weight over e electrons, which is the
    count the loader draws the kick on. With cells set the group is the deposit cell's
    occupants, whose copies stand at many base phases, and the expected value is the
    closed form's for dense phases.

    |b(h)| is invariant under a constant shift of theta, which is what a dump's missing
    reference phase amounts to. This statistic is therefore unaffected by it."""
    out, want = [], []
    for sl in beamio.read_slices(dump_file, LAMBDA0, SPACING):
        if sl["n"] == 0:
            continue
        n_lambda = sl["weight"].sum() / E_CHARGE
        w, theta = sl["weight"], sl["theta"]
        if cells:
            _, inverse = np.unique(cell_groups(sl), return_inverse=True)
            electrons = np.bincount(inverse, weights=w) / E_CHARGE
            m = 256
        else:
            electrons = BEAMLET_SIZE * w[::BEAMLET_SIZE] / E_CHARGE
            m = BEAMLET_SIZE
        for h in harmonics:
            b = np.sum(w * np.exp(-1j * h * theta)) / np.sum(w)
            out.append(abs(b) ** 2 * n_lambda)
            want.append(shot_noise_ratio.slice_ratio(electrons, h, m, (BEAMLET_SIZE - 1) // 2))
    return out, want


def run_mode(exe, lat, workdir, seeds, test_weights, electrons):
    """The generated load over the seeds, at the given electrons a beamlet."""
    samples = {h: ([], []) for h in (1, 2, 3)}
    tag = f"{'w' if test_weights else 'u'}{electrons:g}".replace(".", "p")
    clamped = False
    for seed in range(1, seeds + 1):
        root = f"sn_{tag}_{seed}"
        nml = workdir / f"{root}.toml"
        nml.write_text(to_groups(NML.format(lat=lat, root=root, seed=1000 + 7 * seed,
                                             n_particle=N_PARTICLE, window=N_SLICE * SPACING,
                                             z_min=-N_SLICE * SPACING / 2, z_max=N_SLICE * SPACING / 2,
                                             charge=electrons * E_CHARGE * N_SLICE * N_PARTICLE / BEAMLET_SIZE,
                                             grid_n_pts=GRID_N_PTS, grid_half_width=GRID_HALF_WIDTH,
                                             beamlet_size=BEAMLET_SIZE,
                                             testw="T" if test_weights else "F")))
        r = subprocess.run([str(exe), nml.name], cwd=workdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        clamped = clamped or "noise group of less than one real electron" in r.stdout
        dump = workdir / f"{root}-initial.beam.h5"
        vals, want = b2_samples(dump, (1, 2, 3))
        # The scan runs 150 seeds at 25 MB of particles and 4 MB of field each.
        dump.unlink()
        (workdir / f"{root}-initial.wf.h5").unlink(missing_ok=True)
        for k, h in enumerate((1, 2, 3)):
            samples[h][0].extend(vals[k::3])
            samples[h][1].extend(want[k::3])
    return samples, clamped


def run_load(exe, workdir, root, text):
    """One load_only run from flat keys. Returns its output and its particle dump's path."""
    (workdir / f"{root}.toml").write_text(to_groups(text))
    r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    (workdir / f"{root}-initial.wf.h5").unlink(missing_ok=True)
    return r.stdout, workdir / f"{root}-initial.beam.h5"


def gauss_text(lat, root, seed, limit):
    sig_z = N_SLICE * SPACING / 8
    charge = PEAK * (N_PARTICLE / BEAMLET_SIZE) * E_CHARGE * np.sqrt(2 * np.pi) * sig_z / SPACING
    return GAUSS_NML.format(lat=lat, root=root, seed=seed, n_particle=N_PARTICLE, charge=charge,
                            sig_z=sig_z, window=N_SLICE * SPACING, grid_n_pts=GRID_N_PTS,
                            grid_half_width=GRID_HALF_WIDTH, beamlet_size=BEAMLET_SIZE,
                            limit="" if limit is None else f"  beam.beamlet_min_electrons = {limit:g}\n")


def fewer_beamlets(exe, lat, workdir, seeds, limit, label):
    """The Gaussian bunch with the fewest electrons a beamlet may hold at limit, None
    being the key left at its default, which the load is then held to as one electron."""
    n_min = 1.0 if limit is None else limit
    mbase = N_PARTICLE // BEAMLET_SIZE
    names = ("one beamlet, under the limit", "fewer beamlets", "the input file's beamlets")
    vals = {k: [] for k in names}
    want = {k: [] for k in names}
    worst_charge, count_ok, floor_ok, said_ok, n_fewer = 0.0, True, True, True, 0
    for seed in range(1, seeds + 1):
        tag = "d" if limit is None else f"{limit:g}"
        _, ref = run_load(exe, workdir, f"sn_g0_{tag}_{seed}", gauss_text(lat, f"sn_g0_{tag}_{seed}", 3000 + 7 * seed, 0))
        out, dump = run_load(exe, workdir, f"sn_g_{tag}_{seed}", gauss_text(lat, f"sn_g_{tag}_{seed}", 3000 + 7 * seed, limit))
        fewer_seed = 0
        for sl0, sl in zip(beamio.read_slices(ref, LAMBDA0, SPACING), beamio.read_slices(dump, LAMBDA0, SPACING)):
            q0, q = sl0["weight"].sum(), sl["weight"].sum()
            n_e = q0 / E_CHARGE
            worst_charge = max(worst_charge, abs(q - q0) / q0)
            beamlets = sl["n"] // BEAMLET_SIZE
            # A slice within rounding of a whole multiple of the limit may take either count.
            near = abs(n_e / n_min - round(n_e / n_min)) < 1e-9 * max(1.0, n_e / n_min)
            if not near and beamlets != min(mbase, max(1, int(np.floor(n_e / n_min)))):
                count_ok = False
            per = n_e / beamlets
            if n_e >= n_min and per < n_min * (1 - 1e-12):
                floor_ok = False
            key = names[2] if beamlets == mbase else names[1] if n_e >= n_min else names[0]
            fewer_seed += beamlets < mbase
            w, theta = sl["weight"], sl["theta"]
            for h in (1, 2, 3):
                b = np.sum(w * np.exp(-1j * h * theta)) / np.sum(w)
                vals[key].append(abs(b) ** 2 * n_e)
                want[key].append(shot_noise_ratio.slice_ratio(np.full(beamlets, per), h, BEAMLET_SIZE,
                                                              (BEAMLET_SIZE - 1) // 2))
        n_fewer += fewer_seed
        said_ok = said_ok and f": {fewer_seed} slices have fewer beamlets for it" in out
        ref.unlink()
        dump.unlink()

    ok = True
    for name, good, what in ((f"{label}: each slice's charge is the charge with the key at 0", worst_charge < 1e-12,
                              f"worst {worst_charge:.1e}"),
                             (f"{label}: each slice has min(n/M, max(1, floor(N/n_min))) beamlets", count_ok, ""),
                             (f"{label}: no beamlet is under the limit where the slice holds it", floor_ok, ""),
                             (f"{label}: the loader states the slices with fewer beamlets", said_ok,
                              f"{n_fewer} over the seeds")):
        print(f"--- {name}  {what}  {'ok' if good else 'FAIL'}")
        record(f"shot noise, {name}", None, None, "bool", good)
        ok = ok and good
    for name in names:
        if not vals[name]:
            continue
        v, r = np.asarray(vals[name]), np.asarray(want[name])
        m, rm = v.mean(), r.mean()
        bound = 5 * v.std() / np.sqrt(len(v))
        good = abs(m - rm) < bound
        ok = ok and good
        print(f"--- shot noise, {label}, {name}: <|b(h)|^2 * N_lambda> = {m:.4f} "
              f"(closed form {rm:.4f} +- {bound:.3f}, {len(v)} samples)  {'ok' if good else 'FAIL'}")
        record(f"shot noise, {label}, {name}", m, [rm - bound, rm + bound], "range", good)
    return ok


def same_bits(exe, lat, workdir):
    """A load in which no slice has fewer beamlets, with the key at its default and at 0."""
    def text(root, limit):
        base = NML.format(lat=lat, root=root, seed=1234, n_particle=N_PARTICLE, window=N_SLICE * SPACING,
                          z_min=-N_SLICE * SPACING / 2, z_max=N_SLICE * SPACING / 2,
                          charge=ELECTRONS * E_CHARGE * N_SLICE * N_PARTICLE / BEAMLET_SIZE,
                          grid_n_pts=GRID_N_PTS, grid_half_width=GRID_HALF_WIDTH,
                          beamlet_size=BEAMLET_SIZE, testw="F")
        return base if limit else base.replace("  beam.beamlet_min_electrons = 0\n", "")
    _, a = run_load(exe, workdir, "sn_bits_0", text("sn_bits_0", True))
    _, b = run_load(exe, workdir, "sn_bits_d", text("sn_bits_d", False))
    same = True
    for s0, s1 in zip(beamio.read_slices(a, LAMBDA0, SPACING), beamio.read_slices(b, LAMBDA0, SPACING)):
        for k in ("x", "y", "px", "py", "pz", "theta", "weight"):
            same = same and np.array_equal(s0[k], s1[k])
    a.unlink()
    b.unlink()
    print(f"--- a load with no slice of fewer beamlets is the same bits with the key at its default "
          f"as at 0  {'ok' if same else 'FAIL'}")
    record("shot noise, no slice with fewer beamlets: the same bits as with beam.beamlet_min_electrons = 0",
           None, None, "bool", same)
    return same


KEEP_NML = """[run]
  lattice = "{lat}"
  run.out_root = "{root}"
  run.source_filter = F
  run.ran_seed = {seed}
  run.load_only = T
  window.n_wavelength = 3
/
[beam]
  load_mode = "keep"
  quiet_start = {quiet}
  shot_noise = T
  beamlet_size = {beamlet_size}
{source}/
[field]
  field.lambda0 = 1e-10
  field.seed_power = 0
  field.grid_n_pts = {grid_n_pts}
  field.grid_half_width = {grid_half_width}
/
"""

KEEP_BUNCH = """  beam.n_particle = 2048
  beam.bunch_charge = 4.803322970853e-14
  beam.sig_z = 1.2e-9
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
"""


def run_keep(exe, lat, workdir, seeds, label, source, quiet, want, cells=False):
    """The statistic on a load the loader did not make the beamlets of, or made in keep mode."""
    samples = {h: ([], []) for h in (1, 2, 3)}
    for seed in range(1, seeds + 1):
        root = f"sn_{label}_{seed}"
        (workdir / f"{root}.toml").write_text(to_toml(KEEP_NML.format(lat=lat, root=root, seed=2000 + 7 * seed,
                                                              quiet=quiet, source=source,
                                                              grid_n_pts=GRID_N_PTS,
                                                              grid_half_width=GRID_HALF_WIDTH,
                                                              beamlet_size=BEAMLET_SIZE)))
        r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        if want not in r.stdout:
            record(f"shot noise, {label} load: the noise reported on {want}", None, None, "bool", False)
            print(f"FAIL: {root} did not report the noise on {want!r}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        vals, expect = b2_samples(workdir / f"{root}-initial.beam.h5", (1, 2, 3), cells)
        for k, h in enumerate((1, 2, 3)):
            samples[h][0].extend(vals[k::3])
            samples[h][1].extend(expect[k::3])

    # Every seed's message named the grouping, or the run exited above.
    record(f"shot noise, {label} load: the noise reported on {want}", None, None, "bool", True)
    return samples


def judge(samples, label):
    """The scaled mean against the closed form's, over all harmonics and per harmonic."""
    ok = True
    allv = np.concatenate([samples[h][0] for h in (1, 2, 3)])
    r = float(np.mean(np.concatenate([samples[h][1] for h in (1, 2, 3)])))
    n = len(allv)
    m = allv.mean()
    bound = 5 * r / np.sqrt(n)
    good = abs(m - r) < bound
    ok = ok and good
    print(f"--- shot noise, {label}: <|b(h)|^2 * N_lambda> = {m:.4f} "
          f"(closed form {r:.4f} +- {bound:.3f}, {n} samples)  {'ok' if good else 'FAIL'}")
    record(f"shot noise, {label}", m, [r - bound, r + bound], "range", good)
    for h in (1, 2, 3):
        v = np.asarray(samples[h][0])
        rh = float(np.mean(samples[h][1]))
        mh = v.mean()
        bh = 5 * rh / np.sqrt(len(v))
        goodh = abs(mh - rh) < bh
        ok = ok and goodh
        print(f"      harmonic {h}: {mh:.4f} (closed form {rh:.4f} +- {bh:.3f})  {'ok' if goodh else 'FAIL'}")
        record(f"shot noise, {label}, harmonic {h}", mh, [rh - bh, rh + bh], "range", goodh)
    return ok


def other_loads(exe, lat, workdir, seeds):
    """Keep mode, then copies and cells made outside the loader from one written bunch."""
    ok = True
    src = KEEP_BUNCH + '  write_openpmd_file = "sn_bunch.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, seeds, "keep", src, "T", "the beamlets the loader made"),
               "keep mode, beamlets of a generated bunch") and ok

    b = bunchfile.read_bunch(workdir / "sn_bunch.h5")
    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_copies.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING))
    src = '  beam.position_file = "sn_copies.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, max(3, seeds // 3), "copies", src, "F",
                        "particles sharing their transverse coordinates"),
               "copies made outside the loader") and ok

    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_cells.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING, jitter=1e-12))
    src = '  beam.position_file = "sn_cells.h5"\n'
    ok = judge(run_keep(exe, lat, workdir, max(3, seeds // 3), "cells", src, "F",
                        "the occupants of each deposit cell", cells=True),
               "the occupants of each deposit cell") and ok
    return ok


SPLIT_NML = """  lattice = "{lat}"
  out_root = "sn_split"
  source_filter = F
  dump_field_at = "END"
  ran_seed = 4321
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
&end
"""


def split_check(exe, lat, workdir):
    """
    The stats file's mode split against the same quantity taken from a field dump.

    Every record that takes the field angle moments also reports the power within
    split_angle of the axis, from the transform those moments already compute
    (fel-physics.md, sec-source-filter). A dump's own far field is an independent route
    to the same number, and Parseval makes the two equal. A difference is a bug in the
    masking or the normalization, and the 1e-6 check allows rounding only. One seeded
    steady-state segment is enough. The split depends on how the field is resolved and
    not on what the field is.
    """
    root = "sn_split"
    (workdir / f"{root}.toml").write_text(to_groups(SPLIT_NML.format(lat=lat)))
    r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: split run exited {r.returncode}:\n{r.stdout[-2000:]}")
        return False

    dumps = sorted(workdir.glob(f"{root}-at*.wf.h5")) + sorted(workdir.glob(f"{root}-final.wf.h5"))
    if not dumps:
        print("FAIL: the split run wrote no field dump")
        return False

    with h5py.File(workdir / f"{root}.stats.h5") as h:
        g = h["field/total"]
        angle = float(np.ravel(g["split_angle"][()])[0])
        p_in = np.atleast_1d(g["power_inside_angle"][-1])

    f = fieldio.read_field(str(dumps[-1]))
    u, dx = f["u"], f["dx"]
    n = u.shape[-1]
    kx = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    theta = np.sqrt(kx[None, :] ** 2 + kx[:, None] ** 2) / (2 * np.pi / f["wavelength"])
    spec = np.abs(np.fft.fft2(u, axes=(1, 2))) ** 2 / (n * n)
    inside = (spec * (theta <= angle)).sum(axis=(1, 2)) * dx * dx / (2 * fieldio.MU0_C)

    live = np.isfinite(p_in) & (inside > 0)
    if not live.any():
        print("FAIL: no slice carried both a split and a dump")
        return False
    worst = float(np.max(np.abs(p_in[live] - inside[live]) / inside[live]))
    good = worst <= 1e-6
    print(f"--- mode split against the dump's far field: worst {worst:.3e} over "
          f"{int(live.sum())} slices, angle {angle:.3e} rad (check <= 1e-6)  "
          f"{'ok' if good else 'FAIL'}")
    record("mode split against the dump's far field", worst, 1e-6, "max", good)
    return good


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--lat", default="aramis_1seg.bmad")
    p.add_argument("--seeds", type=int, default=25)
    global PARSED
    PARSED = p.parse_args()

    workdir = pathlib.Path(PARSED.workdir)
    exe = pathlib.Path(PARSED.exe).resolve()

    ok = True
    for electrons in SCAN:
        for test_weights, label in ((False, "uniform weights"), (True, "nonuniform weights (0.25x/1.75x)")):
            samples, clamped = run_mode(exe, PARSED.lat, workdir, PARSED.seeds, test_weights, electrons)
            if electrons != ELECTRONS:
                label += f", {electrons:g} electrons a beamlet"
            ok = judge(samples, label) and ok

            # The loader says so when a beamlet holds less than one electron, and only then.
            # The nonuniform weights at 4 electrons leave the light beamlets one electron to
            # rounding, where either answer is right.
            lightest = electrons * (0.25 if test_weights else 1)
            under = lightest < 1
            if clamped != under and abs(lightest - 1) > 0.01:
                print(f"FAIL: {label}: the loader {'warned' if clamped else 'did not warn'} of a "
                      f"beamlet under one electron")
                ok = False
            if under and abs(lightest - 1) > 0.01:
                record(f"shot noise, {label}: the loader warns of beamlets under one electron",
                       None, None, "bool", clamped)

    ok = fewer_beamlets(exe, PARSED.lat, workdir, PARSED.seeds, None,
                        "Gaussian bunch, beam.beamlet_min_electrons at its default") and ok
    ok = fewer_beamlets(exe, PARSED.lat, workdir, PARSED.seeds, 18.0,
                        "Gaussian bunch, beamlets of 18 electrons or more") and ok
    ok = same_bits(exe, PARSED.lat, workdir) and ok

    ok = other_loads(exe, PARSED.lat, workdir, PARSED.seeds) and ok
    ok = split_check(exe, PARSED.lat, workdir) and ok

    print("shot-noise statistical check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
