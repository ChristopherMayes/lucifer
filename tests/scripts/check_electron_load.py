#!/usr/bin/env python3
"""
Checks of beam.load_mode = "electrons", the load of one particle for every electron.

No reference code enters. The load is held to what independent electrons are: a count in
each slice drawn with the charge of the window held, coordinates from the described
distributions, and a bunching whose mean square is one over the count at every
wavenumber. Genesis4's one4one load is the same physics from another generator, and a
comparison with it could hold only these same statistics.

1. The load. A flat bunch of 256 slices of 3 wavelengths, about 400 electrons a slice, is
   loaded with run.load_only over 16 seeds.
   - The particles are the nearest whole number to the window's charge over e, every
     weight is the window's charge over that number, and their sum is the window's charge.
   - S(k) = sum over a slice's electrons of exp(-i k theta) has <|S|^2> / N = 1 for N
     independent electrons. It is held at k = 1 to 7 and at k = 1/3 and 2/3, the
     wavenumbers a slice of 3 wavelengths has below its first harmonic, each within five
     standard deviations of the mean of an exponential variate.
   - S of neighboring slices is uncorrelated: the mean of Re(S_i conj(S_i+1)) over
     sqrt(N_i N_i+1) is held to zero. Each slice has a random state of its own, and this
     is the test that the states are independent.
   - The variance of the slices' counts is the multinomial's, N p (1 - p).
   - The rms of x, y and pz are those of the matched beam and of beam.sig_pz.
   - The load at 1 thread and at 8 is the same file, dataset for dataset.

2. The keys that do not apply. A value other than the default in beam.quiet_start,
   beam.beamlet_size, beam.n_particle, beam.beamlet_min_electrons or
   beam.shot_noise_method is refused with the key named. The run's own
   <out_root>.resolved.toml, which states every key at its default, loads the same beam.

3. A bunch file. A bunch of single electrons, written by the keep mode with
   beam.write_openpmd_file, is loaded with every particle kept at its own coordinates:
   the same dump as the keep mode's with no copies and no noise. The same bunch at twice
   the charge, whose particles are two electrons each, is refused.

4. Tracking with migration. A beam with a large energy spread is tracked through one
   undulator with run.migrate and run.migrate_check. The check compares, at every
   migration, the slices of fel_migrate_slices with those of fel_migrate_slices_serial
   to the bit, and the run stops where they differ. The run must move particles, leave
   each inside its slice's window, and give the same final field at 1 thread and at 8.

5. The start-up power. Over one undulator the pulse energy of this load and that of the
   weight load of the sample mode, each a mean over 8 seeds of a window of 256 slices,
   agree within five standard errors of their ratio.

Usage: check_electron_load.py --exe <lucifer> --workdir <dir> [--seeds N]
The workdir must hold aramis_1seg.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

import lucifer.beam as beamio
from checkrec import record
from lucifer.flat import to_toml

E_CHARGE = 1.602176634e-19
M_ELECTRON = 0.51099895e6
LAMBDA0 = 1e-10
SAMPLE = 3
SPACING = SAMPLE * LAMBDA0
GAMMA0 = 11357.82
EMIT = 4e-7
SIG_PZ = 8.804506566858e-5
BETA_X, BETA_Y = 8.53711, 17.3899      # aramis_1seg.bmad's beginning Twiss.

# The statistics load: slices, and the electrons of the whole window.
N_SLICE = 256
N_ELECTRON = 102400
WAVENUMBERS = (1, 2, 3, 4, 5, 6, 7, 1 / 3, 2 / 3)

BASE = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "aramis_1seg.bmad"
  out_root = "{root}"
  ran_seed = {seed}
  lambda0 = 1e-10
  n_wavelength = 3
  seed_power = 0
  grid_n_pts = {ngrid}
  grid_half_width = 2e-4
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beam.sig_pz = {sig_pz}
{extra}&end
"""

FLAT = """  beam.bunch_charge = {q:.15e}
  beam.z_distribution = "GRID"
  beam.z_min = -{half:.15e}
  beam.z_max = {half:.15e}
  window_length = {wlen:.15e}
"""


def flat(n_slice, n_electron):
    """The keys of a flat bunch that fills a window of n_slice slices."""
    wlen = n_slice * SPACING
    return FLAT.format(q=n_electron * E_CHARGE, half=wlen / 2, wlen=wlen)


def run(exe, wd, root, extra, seed=777, ngrid=32, sig_pz=SIG_PZ, threads="4", args=()):
    """One run. Returns (returncode, output)."""
    (wd / f"{root}.toml").write_text(to_toml(BASE.format(root=root, seed=seed, ngrid=ngrid,
                                                         sig_pz=sig_pz, extra=extra)))
    r = subprocess.run([str(exe), f"{root}.toml", *args], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    (wd / f"{root}.log").write_text(r.stdout + r.stderr)
    return r.returncode, r.stdout + r.stderr


def must_run(exe, wd, root, extra, **kw):
    rc, out = run(exe, wd, root, extra, **kw)
    if rc != 0:
        print(f"FAIL: {root} exited {rc}:\n{out[-2000:]}")
        sys.exit(1)
    return out


def dumps_equal(fa, fb):
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names = []
        a.visit(lambda n: names.append(n) if isinstance(a[n], h5py.Dataset) else None)
        return all(n in b and a[n].shape == b[n].shape and np.array_equal(a[n][...], b[n][...])
                   for n in names)


def held(name, value, bound, kind="max"):
    """Judge, print and record one value against its bound."""
    if kind == "max":
        good = value <= bound
    elif kind == "min":
        good = value >= bound
    else:
        good = value == bound
    sign = {"max": "<=", "min": ">=", "equal": "=="}[kind]
    print(f"  {name}: {value:.6g} (check {sign} {bound:.6g})  {'ok' if good else 'FAIL'}")
    record(name, value, bound, kind, bool(good))
    return bool(good)


def flag(name, good):
    print(f"  {name}: {'ok' if good else 'FAIL'}")
    record(name, None, None, "bool", bool(good))
    return bool(good)


def the_load(exe, wd, seeds):
    ok = True
    load = flat(N_SLICE, N_ELECTRON) + '  beam.load_mode = "electrons"\n  load_only = T\n'
    q_window = N_ELECTRON * E_CHARGE
    s2 = {k: [] for k in WAVENUMBERS}
    cross, counts, xs, ys, pzs = [], [], [], [], []
    n_worst = w_worst = q_worst = 0.0
    for i in range(seeds):
        root = f"el_load{i}"
        must_run(exe, wd, root, load, seed=4000 + i)
        sls = beamio.read_slices(wd / f"{root}-initial.beam.h5", LAMBDA0, SPACING)
        n = np.array([sl["n"] for sl in sls])
        w = np.concatenate([sl["weight"] for sl in sls])
        n_worst = max(n_worst, abs(int(n.sum()) - N_ELECTRON))
        w_worst = max(w_worst, float(np.max(np.abs(w / (q_window / N_ELECTRON) - 1))))
        q_worst = max(q_worst, abs(float(w.sum()) / q_window - 1))
        counts.append(n)
        s1 = np.array([np.sum(np.exp(-1j * sl["theta"])) for sl in sls])
        cross.append(np.real(s1[:-1] * np.conj(s1[1:])) / np.sqrt(n[:-1] * n[1:]))
        for k in WAVENUMBERS:
            s2[k].append([abs(np.sum(np.exp(-1j * k * sl["theta"]))) ** 2 / sl["n"] for sl in sls])
        xs.append(np.concatenate([sl["x"] for sl in sls]))
        ys.append(np.concatenate([sl["y"] for sl in sls]))
        pzs.append(np.concatenate([sl["pz"] for sl in sls]))

    print(f"--- the electrons load: {seeds} seeds of {N_SLICE} slices, {N_ELECTRON} electrons")
    ok &= held("electrons load: particles less the window's electrons", n_worst, 0, "equal")
    ok &= held("electrons load: weight against the window's charge over its electrons", w_worst, 1e-12)
    ok &= held("electrons load: loaded charge against the window's", q_worst, 1e-12)

    n_sample = seeds * N_SLICE
    worst = max(abs(float(np.mean(s2[k])) - 1) for k in WAVENUMBERS)
    ok &= held(f"electrons load: <|S(k)|^2> / N against 1, worst of {len(WAVENUMBERS)} wavenumbers",
               worst, 5 / np.sqrt(n_sample))
    c = abs(float(np.mean(np.concatenate(cross))))
    ok &= held("electrons load: correlation of S between neighboring slices", c,
               5 * np.sqrt(0.5 / (seeds * (N_SLICE - 1))))
    var = float(np.mean([np.var(n, ddof=1) for n in counts]))
    binom = N_ELECTRON * (1 / N_SLICE) * (1 - 1 / N_SLICE)
    ok &= held("electrons load: variance of the slice counts against the multinomial's",
               abs(var / binom - 1), 5 * np.sqrt(2 / (n_sample - seeds)))

    p0_mc = np.sqrt(GAMMA0 ** 2 - 1)
    n_all = seeds * N_ELECTRON
    for name, v, sigma in (("x", xs, np.sqrt(EMIT / p0_mc * BETA_X)),
                           ("y", ys, np.sqrt(EMIT / p0_mc * BETA_Y)), ("pz", pzs, SIG_PZ)):
        rms = float(np.std(np.concatenate(v)))
        ok &= held(f"electrons load: rms of {name} against the described beam's",
                   abs(rms / sigma - 1), 5 / np.sqrt(2 * n_all))

    for threads in ("1", "8"):
        must_run(exe, wd, f"el_thr{threads}", load, seed=4000, threads=threads)
    ok &= flag("electrons load: the same file at 1 thread and at 8",
               dumps_equal(wd / "el_thr1-initial.beam.h5", wd / "el_thr8-initial.beam.h5"))
    return ok


def the_keys(exe, wd):
    ok = True
    load = flat(64, 6400) + '  beam.load_mode = "electrons"\n  load_only = T\n'
    print("--- the keys the electrons load does not read")
    for key, value in (("quiet_start", "F"), ("beamlet_size", "8"), ("n_particle", "512"),
                       ("beamlet_min_electrons", "0"), ("shot_noise_method", '"first_order"')):
        rc, out = run(exe, wd, f"el_key_{key}", load + f"  beam.{key} = {value}\n")
        ok &= flag(f"electrons load: beam.{key} = {value.strip(chr(34))} is refused with the key named",
                   rc != 0 and f"BEAM.{key.upper()} IS SET" in out)

    must_run(exe, wd, "el_res", load)
    r = subprocess.run([str(exe), "el_res.resolved.toml", "run.out_root=el_res2"], cwd=wd,
                       capture_output=True, text=True, env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    same = r.returncode == 0 and dumps_equal(wd / "el_res-initial.beam.h5", wd / "el_res2-initial.beam.h5")
    ok &= flag("electrons load: the resolved input file loads the same beam", same)
    return ok


def a_bunch_file(exe, wd):
    ok = True
    n = 20000
    bunch = ("  beam.n_particle = {n}\n  beam.sig_z = 3e-9\n  beam.bunch_charge = {q:.15e}\n"
             "  beam.quiet_start = F\n  load_only = T\n")
    print("--- a bunch file of single electrons")
    for tag, per in (("one", 1), ("two", 2)):
        must_run(exe, wd, f"el_mk_{tag}", bunch.format(n=n, q=per * n * E_CHARGE) +
                 f'  beam.load_mode = "keep"\n  beam.write_openpmd_file = "el_bunch_{tag}.h5"\n')
    must_run(exe, wd, "el_file", '  beam.load_mode = "electrons"\n  beam.position_file = "el_bunch_one.h5"\n'
             "  load_only = T\n")
    sls = beamio.read_slices(wd / "el_file-initial.beam.h5", LAMBDA0, SPACING)
    kept = int(sum(sl["n"] for sl in sls))
    ok &= held("electrons from a file: particles kept", kept, n, "equal")
    w = np.concatenate([sl["weight"] for sl in sls])
    ok &= held("electrons from a file: weight against e", float(np.max(np.abs(w / E_CHARGE - 1))), 1e-12)
    must_run(exe, wd, "el_file_keep", '  beam.load_mode = "keep"\n  beam.position_file = "el_bunch_one.h5"\n'
             "  beam.quiet_start = F\n  load_only = T\n")
    ok &= flag("electrons from a file: the keep mode's dump of the same file with no copies",
               dumps_equal(wd / "el_file-initial.beam.h5", wd / "el_file_keep-initial.beam.h5"))
    rc, out = run(exe, wd, "el_file2", '  beam.load_mode = "electrons"\n'
                  '  beam.position_file = "el_bunch_two.h5"\n  load_only = T\n')
    ok &= flag("electrons from a file: particles of two electrons are refused",
               rc != 0 and "NOT" in out and "SINGLE ELECTRONS" in out)
    return ok


def with_migration(exe, wd):
    ok = True
    extra = flat(32, 32 * 2000) + '  beam.load_mode = "electrons"\n  migrate = T\n  migrate_check = T\n' \
        "  write_final_beam = T\n"
    print("--- tracking the electrons load with migration")
    for threads in ("1", "8"):
        must_run(exe, wd, f"el_mig{threads}", extra, sig_pz="5.282703940115e-3", threads=threads)
    moved = 0
    for line in (wd / "el_mig8.migration.txt").read_text().splitlines():
        if line.startswith("moved "):
            moved = int(line.split()[1])
    ok &= held("electrons with migration: moves, each compared with the serial routine", moved, 10000, "min")
    slen = 2 * np.pi * SAMPLE
    n_out = 0
    for sl in beamio.read_slices(wd / "el_mig8-final.beam.h5", LAMBDA0, SPACING):
        n_out += int(np.sum((sl["theta"] < -1e-9) | (sl["theta"] >= slen + 1e-9)))
    ok &= held("electrons with migration: particles outside their slice's window", n_out, 0, "equal")
    same = all(dumps_equal(wd / f"el_mig1-final.{s}.h5", wd / f"el_mig8-final.{s}.h5") for s in ("wf", "beam"))
    ok &= flag("electrons with migration: the same final beam and field at 1 thread and at 8", same)
    return ok


def pulse_energy(wd, root):
    with h5py.File(wd / f"{root}.stats.h5") as h:
        return float(h["run/exit_pulse_energy"][()])


def the_startup(exe, wd):
    seeds, n_slice, per_slice = STARTUP_SEEDS, 256, 200
    bunch = flat(n_slice, n_slice * per_slice)
    energy = {"electrons": [], "weight": []}
    for i in range(seeds):
        must_run(exe, wd, f"el_su_e{i}", bunch + '  beam.load_mode = "electrons"\n', seed=5000 + i)
        must_run(exe, wd, f"el_su_w{i}", bunch + "  beam.shot_noise = T\n  beam.n_particle = 512\n",
                 seed=5000 + i)
        energy["electrons"].append(pulse_energy(wd, f"el_su_e{i}"))
        energy["weight"].append(pulse_energy(wd, f"el_su_w{i}"))
    e, w = np.array(energy["electrons"]), np.array(energy["weight"])
    ratio = float(e.mean() / w.mean())
    sem = float(np.sqrt(e.var(ddof=1) / e.mean() ** 2 + w.var(ddof=1) / w.mean() ** 2) / np.sqrt(seeds))
    print(f"--- the start-up pulse energy over one undulator, {seeds} seeds: electrons "
          f"{e.mean():.4e} J, weight load {w.mean():.4e} J, relative scatter of the ratio {sem:.3f}")
    return held("electrons against the weight load: pulse energy after one undulator, ratio less 1",
                abs(ratio - 1), STARTUP_BOUND)


# The seeds of the start-up comparison and the bound on its ratio. Twelve seeds of windows
# of 512 and of 1024 slices gave ratios of 0.956 and 0.997 with standard errors of 0.031
# and 0.025, and a seed's pulse energy scatters by 0.06 to 0.09 of its mean at 512 slices.
# At 256 slices and 8 seeds the standard error of the ratio is 0.05 to 0.06, and the
# bound is five of them.
STARTUP_SEEDS = 8
STARTUP_BOUND = 0.3


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--seeds", type=int, default=16)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    exe = pathlib.Path(a.exe).resolve()

    ok = the_load(exe, wd, a.seeds)
    ok &= the_keys(exe, wd)
    ok &= a_bunch_file(exe, wd)
    ok &= with_migration(exe, wd)
    ok &= the_startup(exe, wd)

    print("electron-load checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
