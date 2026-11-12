#!/usr/bin/env python3
"""
Statistical check of the weighted shot-noise loader against the closed form of its kick.

Genesis cannot represent nonuniform weights, so there is no reference to compare to.

Physics under test: after quiet loading plus Fawley-style noise, a beamlet standing for
N real electrons has <|b(h)|^2> = R_h(N)/N at every imposed harmonic, where R_h is the
ratio shot_noise_ratio.py computes and docs/shot-noise.md derives for the kick the input
file names, beam.shot_noise_method. Under the first-order kick R_h is 1 - 2/N to first
order and goes to 1 as N grows. A slice's <|b(h)|^2> * N_lambda is then the
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
that the count a beamlet holds is the one stated, and the first-order kick. The level
the loader prints for each load, which fel_shot_ratio computes, is held to
shot_noise_ratio.py's to the four decimals it prints.

A Gaussian bunch is then loaded with the key at its default and at 18 electrons. Its
peak slice holds 50 electrons a beamlet and its last slices two electrons in all, so
the slices away from the peak are loaded with fewer beamlets. Each such slice must
carry the charge it carries with the key at 0, must have the beamlets
min(n_particle/beamlet_size, max(1, floor(N/n_min))) for its N electrons, and its
statistic must agree with the closed form for the beamlets it has. The mean is held
within five standard deviations of the mean, taken from the samples, since a slice of
one beamlet is not a sum of many. A uniform load at 146 electrons a beamlet, where no
slice has fewer, must be the same bits with the key at its default as with the key at 0.

The weight load, beam.shot_noise_method = "weight", counts a whole number of electrons
into every macroparticle and moves no phase. It runs at 146, 18, 4 and 1 electrons a
beamlet and at a quarter of an electron, and its statistic is held at harmonics 1 to 7 to
the bunch's electrons over their nearest whole number, which is 1 to a part in 1e4 here.
The loaded charge is held to the input file's to rounding, every weight to a whole number
of the counted electron's charge, and the variance of the slices' counts to the binomial
variance within five standard deviations. A Gaussian bunch whose last slices hold a
hundredth of an electron is loaded the same way: the charge is held, the slices that drew
no electron have no particles, and the statistic of the others is held to 1. In the keep
mode every particle must be present at the charge the kick load gives it.

A slice of several wavelengths is loaded over all its periods. A Gaussian bunch with
slices of 3 and of 12 wavelengths is loaded under the weight load, the first-order kick
and no noise, and every slice's sum of w exp(-i (m/n) theta) is taken at the m that are
not multiples of n, the slice's wavelengths: the wavenumbers at which a load that stood
its particles in one period had lines. With no noise the sum over the slice's charge is
held to rounding. Under the weight load its square over the slice's electrons is held to
1 within five standard deviations. The charge in each period of a slice is held to the
binomial spread at the load. In the dump at the end of a tracked segment, with migration
off and on, it is held to be no wider than that of the slice's macroparticles placed at
random.

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
import re
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

# The loader's own statement of its level at the fundamental, from fel_shot_ratio, and
# the worst difference from shot_noise_ratio.py over the scans of the first-order kick.
# The line prints four decimals.
STATED = re.compile(r"noise power loaded over the shot-noise power, fundamental:\s+([0-9.]+)")
stated_worst = {"first_order": 0.0}

# The electrons a beamlet holds in each scan of the weight load, and the harmonics its
# statistic is held at: every harmonic below the beamlet size.
WEIGHT_SCAN = (ELECTRONS, 18.0, 4.0, 1.0, 0.25)
WEIGHT_HARMONICS = tuple(range(1, BEAMLET_SIZE))

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
  beam.shot_noise_method = "{kick}"
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
  beam.shot_noise_method = "{kick}"
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
    closed form's for dense phases. The load is under the first-order kick.

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
                                             beamlet_size=BEAMLET_SIZE, kick="first_order",
                                             testw="T" if test_weights else "F")))
        r = subprocess.run([str(exe), nml.name], cwd=workdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        clamped = clamped or "noise group of less than one real electron" in r.stdout
        dump = workdir / f"{root}-initial.beam.h5"
        vals, want = b2_samples(dump, (1, 2, 3))
        said = STATED.search(r.stdout)
        if not said:
            print(f"FAIL: loader run {root} did not state its noise level:\n{r.stdout[-2000:]}")
            sys.exit(1)
        stated_worst["first_order"] = max(stated_worst["first_order"],
                                          abs(float(said.group(1)) - float(np.mean(want[0::3]))))
        # The scan runs 150 seeds at 25 MB of particles and 4 MB of field each.
        dump.unlink()
        (workdir / f"{root}-initial.wf.h5").unlink(missing_ok=True)
        for k, h in enumerate((1, 2, 3)):
            samples[h][0].extend(vals[k::3])
            samples[h][1].extend(want[k::3])
    return samples, clamped


def weight_statistics(dump, charge):
    """One weight-load dump: the statistic of every slice that has particles at
    WEIGHT_HARMONICS, the slices' electron counts, the relative error of the loaded charge
    against charge, and whether every weight is a whole number of the counted electron's
    charge. That charge is the dump's total over its nearest whole number of electrons."""
    slices = list(beamio.read_slices(dump, LAMBDA0, SPACING))
    total = float(sum(sl["weight"][: sl["n"]].sum() for sl in slices))
    unit = total / max(1.0, np.rint(total / E_CHARGE))
    vals = {h: [] for h in WEIGHT_HARMONICS}
    counts, whole = [], True
    for sl in slices:
        w = sl["weight"][: sl["n"]]
        counts.append(w.sum() / unit)
        if sl["n"] == 0 or w.sum() == 0:
            continue
        whole = whole and bool(np.all(np.abs(w / unit - np.rint(w / unit)) < 1e-6))
        whole = whole and sl["n"] % BEAMLET_SIZE == 0
        for h in WEIGHT_HARMONICS:
            b = np.sum(w * np.exp(-1j * h * sl["theta"][: sl["n"]])) / w.sum()
            vals[h].append(abs(b) ** 2 * w.sum() / E_CHARGE)
    return vals, np.asarray(counts), abs(total - charge) / charge, whole


def weight_scan(exe, lat, workdir, seeds):
    """The weight load of the uniform bunch over the seeds at each count of WEIGHT_SCAN."""
    ok = True
    worst_charge, whole_ok, stated = 0.0, True, 0.0
    for electrons in WEIGHT_SCAN:
        charge = electrons * E_CHARGE * N_SLICE * N_PARTICLE / BEAMLET_SIZE
        want = float(shot_noise_ratio.ratio_weight_total(charge / E_CHARGE))
        vals = {h: [] for h in WEIGHT_HARMONICS}
        counts = []
        tag = f"{electrons:g}".replace(".", "p")
        for seed in range(1, seeds + 1):
            root = f"sn_wt{tag}_{seed}"
            out, dump = run_load(exe, workdir, root, NML.format(
                lat=lat, root=root, seed=1000 + 7 * seed, n_particle=N_PARTICLE, window=N_SLICE * SPACING,
                z_min=-N_SLICE * SPACING / 2, z_max=N_SLICE * SPACING / 2, charge=charge, grid_n_pts=GRID_N_PTS,
                grid_half_width=GRID_HALF_WIDTH, beamlet_size=BEAMLET_SIZE, kick="weight", testw="F"))
            v, c, dq, whole = weight_statistics(dump, charge)
            dump.unlink()
            said = STATED.search(out)
            stated = max(stated, abs(float(said.group(1)) - want) if said else 1.0)
            worst_charge = max(worst_charge, dq)
            whole_ok = whole_ok and whole
            counts.append(c)
            for h in WEIGHT_HARMONICS:
                vals[h].extend(v[h])

        count = f"{electrons:g} electron{'' if electrons == 1 else 's'} a beamlet"
        allv = np.concatenate([vals[h] for h in WEIGHT_HARMONICS])
        bound = 5 * want / np.sqrt(len(allv))
        good = abs(allv.mean() - want) < bound
        ok = ok and good
        print(f"--- shot noise, weight load, uniform weights, {count}: <|b(h)|^2 * N_lambda> = {allv.mean():.4f} "
              f"(closed form {want:.4f} +- {bound:.3f}, {len(allv)} samples)  {'ok' if good else 'FAIL'}")
        record(f"shot noise, weight load, uniform weights, {count}", allv.mean(), [want - bound, want + bound],
               "range", good)
        for h in WEIGHT_HARMONICS:
            v = np.asarray(vals[h])
            bh = 5 * want / np.sqrt(len(v))
            goodh = abs(v.mean() - want) < bh
            ok = ok and goodh
            print(f"      harmonic {h}: {v.mean():.4f} (closed form {want:.4f} +- {bh:.3f})  {'ok' if goodh else 'FAIL'}")
            record(f"shot noise, weight load, uniform weights, {count}, harmonic {h}", v.mean(),
                   [want - bh, want + bh], "range", goodh)

        # A slice's count is binomial: the bunch's electrons, each in this slice with probability 1/N_SLICE.
        c = np.concatenate(counts)
        n_total = np.rint(charge / E_CHARGE)
        variance = n_total * (1 - 1 / N_SLICE) / N_SLICE
        ratio = float(np.mean((c - n_total / N_SLICE) ** 2) / variance)
        bound = 5 * np.sqrt(2 / len(c))
        good = abs(ratio - 1) < bound
        ok = ok and good
        print(f"      the variance of the slices' electron counts over the binomial variance: {ratio:.4f} "
              f"(1 +- {bound:.3f}, {len(c)} slices)  {'ok' if good else 'FAIL'}")
        record(f"shot noise, weight load, {count}: the variance of the slices' counts over the binomial variance",
               ratio, [1 - bound, 1 + bound], "range", good)

    for name, value, limit in (("the loaded charge against the input file's, worst over the scans", worst_charge, 1e-12),
                               ("the level the loader states against the closed form", stated, 1e-4)):
        good = value <= limit
        ok = ok and good
        print(f"--- weight load, {name}: {value:.1e} (check <= {limit:g})  {'ok' if good else 'FAIL'}")
        record(f"shot noise, weight load, {name}", value, limit, "max", good)
    print(f"--- weight load, every weight is a whole number of counted electrons, in whole beamlets  "
          f"{'ok' if whole_ok else 'FAIL'}")
    record("shot noise, weight load, every weight is a whole number of counted electrons, in whole beamlets",
           None, None, "bool", whole_ok)
    return ok and whole_ok


def weight_gauss(exe, lat, workdir, seeds):
    """The Gaussian bunch under the weight load. Its last slices hold two electrons at PEAK
    electrons a beamlet in the peak slice, so the peak is lowered until they hold a
    hundredth of an electron and most of them draw none."""
    sig_z = N_SLICE * SPACING / 8
    scale = 0.005
    charge = scale * PEAK * (N_PARTICLE / BEAMLET_SIZE) * E_CHARGE * np.sqrt(2 * np.pi) * sig_z / SPACING
    vals, worst_charge, whole_ok, n_empty, n_slices = [], 0.0, True, 0, 0
    for seed in range(1, seeds + 1):
        root = f"sn_wg_{seed}"
        text = gauss_text(lat, root, 3000 + 7 * seed, None, "weight")
        line = [x for x in text.splitlines() if "beam.bunch_charge" in x][0]
        _, dump = run_load(exe, workdir, root, text.replace(line, f"  beam.bunch_charge = {charge:.12e}"))
        v, c, dq, whole = weight_statistics(dump, charge)
        # The window is eight sigmas, so the charge loaded is the bunch's less its tails outside.
        worst_charge = max(worst_charge, abs(c.sum() - np.rint(c.sum())))
        whole_ok = whole_ok and whole
        n_empty += int(np.count_nonzero(c == 0))
        n_slices += len(c)
        for sl in beamio.read_slices(dump, LAMBDA0, SPACING):
            whole_ok = whole_ok and (sl["n"] == 0) == (sl["weight"][: sl["n"]].sum() == 0)
        dump.unlink()
        vals.extend(np.concatenate([v[h] for h in WEIGHT_HARMONICS]))
    v = np.asarray(vals)
    bound = 5 * v.std() / np.sqrt(len(v))
    good = abs(v.mean() - 1) < bound
    print(f"--- shot noise, weight load, Gaussian bunch, the slices with particles: <|b(h)|^2 * N_lambda> = "
          f"{v.mean():.4f} (closed form 1.0000 +- {bound:.3f}, {len(v)} samples)  {'ok' if good else 'FAIL'}")
    record("shot noise, weight load, Gaussian bunch, the slices with particles", v.mean(), [1 - bound, 1 + bound],
           "range", good)
    empty_ok = n_empty > 0 and whole_ok
    print(f"--- weight load, Gaussian bunch: {n_empty} of {n_slices} slices drew no electron and have no particles, "
          f"and the bunch is a whole number of counted electrons to {worst_charge:.1e}  "
          f"{'ok' if empty_ok and worst_charge < 1e-6 else 'FAIL'}")
    record("shot noise, weight load, Gaussian bunch: a slice that drew no electron has no particles", None, None,
           "bool", empty_ok)
    record("shot noise, weight load, Gaussian bunch: the loaded electrons are a whole number", worst_charge, 1e-6,
           "max", worst_charge < 1e-6)
    return good and empty_ok and worst_charge < 1e-6


def weight_kept(exe, lat, workdir, seeds):
    """Keep mode under the weight load, against the first-order kick's load of the same
    bunch: every particle is present, and its copies hold the charge the kick load's hold.
    The statistic is held to the keep mode's form, 1 at these 146 electrons a particle."""
    src = KEEP_BUNCH
    vals, worst, same_count = [], 0.0, True
    for seed in range(1, seeds + 1):
        dumps = {}
        for kick in ("first_order", "weight"):
            root = f"sn_keepw_{kick[0]}_{seed}"
            (workdir / f"{root}.toml").write_text(to_toml(KEEP_NML.format(
                lat=lat, root=root, seed=2000 + 7 * seed, quiet="T", source=src, kick=kick, grid_n_pts=GRID_N_PTS,
                grid_half_width=GRID_HALF_WIDTH, beamlet_size=BEAMLET_SIZE)))
            r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
                sys.exit(1)
            dumps[kick] = workdir / f"{root}-initial.beam.h5"
        for s0, s1 in zip(beamio.read_slices(dumps["first_order"], LAMBDA0, SPACING),
                          beamio.read_slices(dumps["weight"], LAMBDA0, SPACING)):
            same_count = same_count and s0["n"] == s1["n"]
            if s0["n"] == 0 or s0["n"] != s1["n"]:
                continue
            q0 = s0["weight"].reshape(-1, BEAMLET_SIZE).sum(axis=1)
            q1 = s1["weight"].reshape(-1, BEAMLET_SIZE).sum(axis=1)
            worst = max(worst, float(np.max(np.abs(q1 - q0) / q0)))
            w, theta = s1["weight"], s1["theta"]
            for h in WEIGHT_HARMONICS:
                b = np.sum(w * np.exp(-1j * h * theta)) / w.sum()
                vals.append(abs(b) ** 2 * w.sum() / E_CHARGE)
        for d in dumps.values():
            d.unlink()
    v = np.asarray(vals)
    bound = 5 / np.sqrt(len(v))
    good = abs(v.mean() - 1) < bound
    print(f"--- shot noise, weight load, keep mode: <|b(h)|^2 * N_lambda> = {v.mean():.4f} "
          f"(closed form 1.0000 +- {bound:.3f}, {len(v)} samples)  {'ok' if good else 'FAIL'}")
    record("shot noise, weight load, keep mode", v.mean(), [1 - bound, 1 + bound], "range", good)
    held = same_count and worst < 1e-12
    print(f"--- weight load, keep mode: every particle is present at its own charge, worst {worst:.1e}  "
          f"{'ok' if held else 'FAIL'}")
    record("shot noise, weight load, keep mode: every particle of the bunch is present at its own charge",
           worst, 1e-12, "max", held)
    return good and held


def fill_statistics(dump, n_period):
    """One dump of slices n_period wavelengths long: every slice's |sum w exp(-i (m/n) theta)|
    over its charge at the m below 2 n that are not multiples of n, the same squared over
    the slice's electrons, and each period's count of electrons less its share, squared,
    over the binomial variance of that share: the electrons' own, and the variance of the
    slice's macroparticles placed at random with the counts they carry."""
    spacing = n_period * LAMBDA0
    slices = list(beamio.read_slices(dump, LAMBDA0, spacing))
    total = float(sum(sl["weight"][: sl["n"]].sum() for sl in slices))
    unit = total / max(1.0, np.rint(total / E_CHARGE))
    orders = [m for m in range(1, 2 * n_period) if m % n_period]
    bunching, power, spread, placed = [], [], [], []
    for sl in slices:
        n = sl["n"]
        w = sl["weight"][:n]
        if n == 0 or w.sum() == 0:
            continue
        theta = sl["theta"][:n]
        s = np.array([np.sum(w * np.exp(-1j * m * theta / n_period)) for m in orders])
        bunching.append(np.abs(s).max() / w.sum())
        power.extend(np.abs(s) ** 2 / (w.sum() * unit))
        electrons = w.sum() / unit
        period = np.floor(theta / (2 * np.pi)).astype(int) % n_period
        counts = np.bincount(period, weights=w / unit, minlength=n_period)
        if electrons >= 8 * n_period:
            share = (1 - 1 / n_period) / n_period
            spread.extend((counts - electrons / n_period) ** 2 / (electrons * share))
            placed.extend((counts - electrons / n_period) ** 2 / (np.sum((w / unit) ** 2) * share))
    return np.asarray(bunching), np.asarray(power), np.asarray(spread), np.asarray(placed)


def slice_fill(exe, lat, workdir, seeds):
    """A Gaussian bunch in slices of 3 and of 12 wavelengths, loaded over the slices' periods."""
    ok = True
    sig_z = N_SLICE * SPACING / 8
    charge = PEAK * (N_PARTICLE / BEAMLET_SIZE) * E_CHARGE * np.sqrt(2 * np.pi) * sig_z / SPACING

    def text(root, seed, n_period, noise, track=False):
        t = gauss_text(lat, root, seed, None, "weight")
        t = t.replace("  n_wavelength = 3\n", f"  n_wavelength = {n_period}\n")
        t = t.replace('  shot_noise = T\n  beam.shot_noise_method = "weight"\n', noise)
        if n_period != 3:
            scale = n_period / 3
            for key, value in (("beam.sig_z", sig_z * scale), ("window_length", N_SLICE * SPACING * scale),
                               ("beam.bunch_charge", charge * scale)):
                line = [x for x in t.splitlines() if x.strip().startswith(key + " =")][0]
                t = t.replace(line, f"  {key} = {value:.12e}")
        return t.replace("  load_only = T\n", "") if track else t

    methods = (("weight", '  shot_noise = T\n  beam.shot_noise_method = "weight"\n'),
               ("first_order", '  shot_noise = T\n  beam.shot_noise_method = "first_order"\n'),
               ("no noise", "  shot_noise = F\n"))
    for n_period in (3, 12):
        for name, noise in methods:
            bunching, power, spread = [], [], []
            for seed in range(1, (seeds if name == "weight" else 2) + 1):
                root = f"sn_fill{n_period}_{name[0]}_{seed}"
                _, dump = run_load(exe, workdir, root, text(root, 5000 + 7 * seed, n_period, noise))
                b, p, c, _ = fill_statistics(dump, n_period)
                dump.unlink()
                bunching.extend(b)
                power.extend(p)
                spread.extend(c)
            label = f"slices of {n_period} wavelengths, {name}"
            if name == "no noise":
                worst = float(np.max(bunching))
                good = worst <= 1e-12
                print(f"--- {label}: worst |sum w exp(-i (m/n) theta)| over a slice's charge, m not a multiple of n: "
                      f"{worst:.1e} (check <= 1e-12)  {'ok' if good else 'FAIL'}")
                record(f"shot noise, {label}: bunching of a slice at the multiples of k_r/n that are not harmonics",
                       worst, 1e-12, "max", good)
            elif name == "weight":
                v = np.asarray(power)
                bound = 5 * v.std() / np.sqrt(len(v))
                good = abs(v.mean() - 1) < bound
                print(f"--- {label}: noise power at the multiples of k_r/n that are not harmonics over the shot "
                      f"noise: {v.mean():.4f} (1 +- {bound:.3f}, {len(v)} samples)  {'ok' if good else 'FAIL'}")
                record(f"shot noise, {label}: noise power at the multiples of k_r/n that are not harmonics",
                       v.mean(), [1 - bound, 1 + bound], "range", good)
                c = np.asarray(spread)
                cb = 5 * np.sqrt(2 / len(c))
                goodc = abs(c.mean() - 1) < cb
                print(f"      the spread of the electrons in a period of a slice over the binomial spread: "
                      f"{c.mean():.4f} (1 +- {cb:.3f}, {len(c)} periods)  {'ok' if goodc else 'FAIL'}")
                record(f"shot noise, {label}: the spread of a period's electrons over the binomial spread",
                       c.mean(), [1 - cb, 1 + cb], "range", goodc)
                good = good and goodc
            else:
                v = np.asarray(power)
                good = True
                print(f"--- {label}: noise power at the multiples of k_r/n that are not harmonics over the shot "
                      f"noise: {v.mean():.4f} (recorded, not judged)")
                record(f"shot noise, {label}: noise power at the multiples of k_r/n that are not harmonics",
                       v.mean(), 0.0, "min", True)
            ok = ok and good

    # The dump at the end of a tracked segment, with migration off and on. The particles
    # have slipped by a fraction of a period, each with the electrons it carries, so a
    # period's count is no longer the electrons' binomial one. It is held to be no wider
    # than the count of the same macroparticles placed at random. A slice whose particles
    # stood in one period reads a third of its macroparticle count here, hundreds.
    for migrate in ("F", "T"):
        root = f"sn_fillt_{migrate}"
        t = text(root, 5007, 3, methods[0][1], track=True) .replace("&end", f"  migrate = {migrate}\n&end")
        (workdir / f"{root}.toml").write_text(to_groups(t))
        r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: the tracked run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        _, _, _, c = fill_statistics(workdir / f"{root}-final.beam.h5", 3)
        for f in workdir.glob(f"{root}-*.h5"):
            f.unlink()
        cb = 1 + 5 * np.sqrt(2 / len(c))
        good = c.mean() <= cb
        ok = ok and good
        name = "with migration" if migrate == "T" else "without migration"
        print(f"--- slices of 3 wavelengths after a tracked segment, {name}: the spread of the charge in a "
              f"period over that of its macroparticles placed at random: {c.mean():.4f} (check <= {cb:.3f}, "
              f"{len(c)} periods)  {'ok' if good else 'FAIL'}")
        record(f"shot noise, slices of 3 wavelengths after a tracked segment, {name}: the spread of a period's "
               f"charge over that of its macroparticles placed at random", c.mean(), cb, "max", good)
    return ok


def run_load(exe, workdir, root, text):
    """One load_only run from flat keys. Returns its output and its particle dump's path."""
    (workdir / f"{root}.toml").write_text(to_groups(text))
    r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    (workdir / f"{root}-initial.wf.h5").unlink(missing_ok=True)
    return r.stdout, workdir / f"{root}-initial.beam.h5"


def gauss_text(lat, root, seed, limit, kick):
    sig_z = N_SLICE * SPACING / 8
    charge = PEAK * (N_PARTICLE / BEAMLET_SIZE) * E_CHARGE * np.sqrt(2 * np.pi) * sig_z / SPACING
    return GAUSS_NML.format(lat=lat, root=root, seed=seed, n_particle=N_PARTICLE, charge=charge,
                            sig_z=sig_z, window=N_SLICE * SPACING, grid_n_pts=GRID_N_PTS,
                            grid_half_width=GRID_HALF_WIDTH, beamlet_size=BEAMLET_SIZE, kick=kick,
                            limit="" if limit is None else f"  beam.beamlet_min_electrons = {limit:g}\n")


def fewer_beamlets(exe, lat, workdir, seeds, limit, label):
    """The Gaussian bunch with the fewest electrons a beamlet may hold at limit, None
    being the key left at its default, which the load is then held to as one electron.
    The load is under the first-order kick, which the rule is for."""
    kick = "first_order"
    n_min = 1.0 if limit is None else limit
    mbase = N_PARTICLE // BEAMLET_SIZE
    names = ("one beamlet, under the limit", "fewer beamlets", "the input file's beamlets")
    vals = {k: [] for k in names}
    want = {k: [] for k in names}
    worst_charge, count_ok, floor_ok, said_ok, n_fewer = 0.0, True, True, True, 0
    for seed in range(1, seeds + 1):
        tag = "d" if limit is None else f"{limit:g}"
        _, ref = run_load(exe, workdir, f"sn_g0_{tag}_{seed}",
                          gauss_text(lat, f"sn_g0_{tag}_{seed}", 3000 + 7 * seed, 0, kick))
        out, dump = run_load(exe, workdir, f"sn_g_{tag}_{seed}",
                             gauss_text(lat, f"sn_g_{tag}_{seed}", 3000 + 7 * seed, limit, kick))
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
                          beamlet_size=BEAMLET_SIZE, testw="F", kick="first_order")
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


def empty_slices(exe, lat, workdir):
    """A flat bunch over the middle half of the window, with the key at its default and at 0.
    At the default a slice with no charge has no particles, and at 0 it has the input file's
    count at zero weight. A slice with charge is the same either way."""
    def text(root, limit):
        base = NML.format(lat=lat, root=root, seed=4321, n_particle=N_PARTICLE, window=N_SLICE * SPACING,
                          z_min=-N_SLICE * SPACING / 4, z_max=N_SLICE * SPACING / 4,
                          charge=ELECTRONS * E_CHARGE * N_SLICE * N_PARTICLE / BEAMLET_SIZE,
                          grid_n_pts=GRID_N_PTS, grid_half_width=GRID_HALF_WIDTH,
                          beamlet_size=BEAMLET_SIZE, testw="F", kick="first_order")
        return base if limit else base.replace("  beam.beamlet_min_electrons = 0\n", "")
    _, a = run_load(exe, workdir, "sn_empty_0", text("sn_empty_0", True))
    _, b = run_load(exe, workdir, "sn_empty_d", text("sn_empty_d", False))
    n_empty, none_ok, same_ok, worst = 0, True, True, 0.0
    for s0, s1 in zip(beamio.read_slices(a, LAMBDA0, SPACING), beamio.read_slices(b, LAMBDA0, SPACING)):
        q0 = s0["weight"].sum()
        if q0 == 0:
            n_empty += 1
            none_ok = none_ok and s0["n"] == N_PARTICLE and s1["n"] == 0
        else:
            same_ok = same_ok and s1["n"] == s0["n"]
            worst = max(worst, abs(s1["weight"].sum() - q0) / q0)
    a.unlink()
    b.unlink()
    none_ok = none_ok and n_empty > 0
    same_ok = same_ok and worst < 1e-12
    print(f"--- a slice with no charge has no particles: {n_empty} such slices  {'ok' if none_ok else 'FAIL'}")
    record("shot noise, a slice with no charge is loaded with no particles", None, None, "bool", none_ok)
    print(f"--- the slices with charge beside them keep their particles and their charge, "
          f"worst {worst:.1e}  {'ok' if same_ok else 'FAIL'}")
    record("shot noise, the slices with charge beside slices with none keep their particles and charge",
           None, None, "bool", same_ok)
    return none_ok and same_ok


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
  shot_noise_method = "{kick}"
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
                                                              quiet=quiet, source=source, kick="first_order",
                                                              grid_n_pts=GRID_N_PTS,
                                                              grid_half_width=GRID_HALF_WIDTH,
                                                              beamlet_size=BEAMLET_SIZE)))
        r = subprocess.run([str(exe), f"{root}.toml"], cwd=workdir, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: loader run {root} exited {r.returncode}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        said = f"shot noise, {label} load: the noise reported on {want}"
        if want not in r.stdout:
            record(said, None, None, "bool", False)
            print(f"FAIL: {root} did not report the noise on {want!r}:\n{r.stdout[-2000:]}")
            sys.exit(1)
        vals, expect = b2_samples(workdir / f"{root}-initial.beam.h5", (1, 2, 3), cells)
        for k, h in enumerate((1, 2, 3)):
            samples[h][0].extend(vals[k::3])
            samples[h][1].extend(expect[k::3])

    # Every seed's message named the grouping, or the run exited above.
    record(said, None, None, "bool", True)
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
    """Keep mode, then copies and cells made outside the loader from one written bunch,
    under the first-order kick."""
    ok = True

    def both(label, name, src, quiet, want, n_seed, cells=False):
        return judge(run_keep(exe, lat, workdir, n_seed, label, src, quiet, want, cells), name)

    src = KEEP_BUNCH + '  write_openpmd_file = "sn_bunch.h5"\n'
    ok = both("keep", "keep mode, beamlets of a generated bunch", src, "T", "the beamlets the loader made", seeds) and ok

    b = bunchfile.read_bunch(workdir / "sn_bunch.h5")
    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_copies.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING))
    src = '  beam.position_file = "sn_copies.h5"\n'
    ok = both("copies", "copies made outside the loader", src, "F",
              "particles sharing their transverse coordinates", max(3, seeds // 3)) and ok

    bunchfile.write_bunch(workdir / "sn_bunch.h5", workdir / "sn_cells.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING, jitter=1e-12))
    src = '  beam.position_file = "sn_cells.h5"\n'
    ok = both("cells", "the occupants of each deposit cell", src, "F",
              "the occupants of each deposit cell", max(3, seeds // 3), cells=True) and ok
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

    good = stated_worst["first_order"] <= 1e-4
    ok = ok and good
    print(f"--- the level the loader states for the first-order kick against the closed form: worst "
          f"{stated_worst['first_order']:.1e} over the scans (check <= 1e-4)  {'ok' if good else 'FAIL'}")
    record("shot noise, the level the loader states for the first-order kick against the closed form",
           stated_worst["first_order"], 1e-4, "max", good)

    ok = weight_scan(exe, PARSED.lat, workdir, PARSED.seeds) and ok
    ok = weight_gauss(exe, PARSED.lat, workdir, PARSED.seeds) and ok
    ok = weight_kept(exe, PARSED.lat, workdir, max(3, PARSED.seeds // 3)) and ok
    ok = slice_fill(exe, PARSED.lat, workdir, max(3, PARSED.seeds // 3)) and ok

    ok = fewer_beamlets(exe, PARSED.lat, workdir, PARSED.seeds, None,
                        "Gaussian bunch, beam.beamlet_min_electrons at its default") and ok
    ok = fewer_beamlets(exe, PARSED.lat, workdir, PARSED.seeds, 18.0,
                        "Gaussian bunch, beamlets of 18 electrons or more") and ok
    ok = same_bits(exe, PARSED.lat, workdir) and ok
    ok = empty_slices(exe, PARSED.lat, workdir) and ok

    ok = other_loads(exe, PARSED.lat, workdir, PARSED.seeds) and ok
    ok = split_check(exe, PARSED.lat, workdir) and ok

    print("shot-noise statistical check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
