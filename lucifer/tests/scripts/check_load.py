#!/usr/bin/env python3
"""
Check the load path (fel-physics.md sec-loading).

Four things, each on the check built for it:

  keep exactness   load_mode = "keep" hands every particle of the bunch through as a
                   beamlet of copies. The per-slice currents and charge-weighted moments
                   the loader records are recomputed here from the bunch file with the same
                   binning, and must agree to roundoff. The dump the load writes must carry
                   the bunch's charge.
  split weights    the same load with every particle split into coincident copies at a
                   third and two thirds of its charge must record the same currents and
                   moments to roundoff.
  in-cone startup  a generated Gaussian bunch in sample mode radiates, inside 3 urad at
                   the end of the first undulator, the physical spontaneous power of the
                   beam within that angle: 0.133 MW per slice at 3 kA on 1.57 um cells
                   (doc/startup-noise.md, the interior mean of a flat window). A field
                   record collects the emission of every beam slice it slipped past, so
                   for a bunch the comparable quantity is the sum over a window long
                   enough that nothing escapes: the in-cone power summed over the records,
                   over the summed current in units of 3 kA, read from the stats split.
  refusals         noise asked for on a load that is not quiet is refused with the floor
                   in the message: the generator without the quiet start, and a bunch of
                   copies whose phases were shifted singly. The unshifted copies load, and
                   the message says the noise went on particles sharing their coordinates.

Usage:

  check_load.py --exe <lucifer> --workdir <dir>
"""

import argparse
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import beamio  # noqa: E402
import bunchfile  # noqa: E402

LAMBDA0 = 1e-10
SAMPLE = 3
SPACING = SAMPLE * LAMBDA0
E_CHARGE = 1.602176634e-19

DECK = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%ran_seed = {seed}
{pextra}  slicing%n_wavelength = {sample}
/
&fel_beam_init
{bextra}/
&fel_wavefront_init
  wavefront_init%lambda0 = {lam:.6e}
  wavefront_init%seed_power = 0
  wavefront_init%grid_n_pts = {ngrid}
  wavefront_init%grid_half_width = 2e-4
/
"""

BUNCH = """  beam_init%n_particle = 40000
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beam_init%sig_z = 3e-10
  beam_init%sig_pz = 8.8e-5
  beam_init%bunch_charge = 1.5e-14
"""


def run(exe, wd, root, bextra, pextra="", lat="aramis_1seg.bmad", seed=777, ngrid=64,
        sample=SAMPLE, threads="4"):
    """One run. Returns (returncode, stdout)."""
    (wd / f"{root}.nml").write_text(DECK.format(lat=lat, root=root, seed=seed, pextra=pextra,
                                                  bextra=bextra, lam=LAMBDA0, ngrid=ngrid, sample=sample))
    r = subprocess.run([exe, f"{root}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    (wd / f"{root}.log").write_text(r.stdout + r.stderr)
    return r.returncode, r.stdout + r.stderr


def load_record(wd, root):
    """The loader's per-slice current, count and five charge-weighted moments."""
    rows = []
    for line in (wd / f"{root}.import.txt").read_text().splitlines():
        if line.startswith("slice "):
            f = line.split()
            rows.append([float(f[2]), int(f[3])] + [float(v) for v in f[4:9]])
    return np.array(rows)


def bunch_record(b, n_copy):
    """The same record recomputed from the bunch: the loader's binning and its chart."""
    s = bunchfile.arrival(b)
    isl, nslice, _ = bunchfile.slice_index(s, SPACING)
    p0_mc = b["p0c"] / bunchfile.M_ELECTRON            # The lattice's reference, the loader's too.
    px = b["px"] / bunchfile.EV_C / (bunchfile.M_ELECTRON * p0_mc)
    py = b["py"] / bunchfile.EV_C / (bunchfile.M_ELECTRON * p0_mc)
    p_mc = np.sqrt(b["px"] ** 2 + b["py"] ** 2 + b["pz"] ** 2) / bunchfile.EV_C / bunchfile.M_ELECTRON
    pz = (p_mc - p0_mc) / p0_mc
    rows = []
    for i in range(1, nslice + 1):
        m = isl == i
        w = b["w"][m]
        cur = bunchfile.C_LIGHT * w.sum() / SPACING
        if w.sum() > 0:
            mom = [np.sum(w * q[m]) / w.sum() for q in (b["x"], px, b["y"], py, pz)]
        else:
            mom = [0.0] * 5
        rows.append([cur, int(m.sum()) * n_copy] + mom)
    return np.array(rows)


def compare(a, b, tol, what):
    """Rows of the two records against each other, relative to each column's scale."""
    if a.shape != b.shape:
        print(f"FAIL: {what}: {a.shape[0]} slices against {b.shape[0]}")
        return False
    scale = np.maximum(np.abs(b).max(axis=0), 1e-300)
    worst = float(np.max(np.abs(a - b) / scale))
    good = worst <= tol
    print(f"  {what}: worst relative difference {worst:.3e} over {a.shape[0]} slices "
          f"(check <= {tol:.0e})  {'ok' if good else 'FAIL'}")
    return good


def keep_exactness(exe, wd):
    ok = True
    bextra = BUNCH + '  load_mode = "keep"\n  quiet_start = T\n  shot_noise = F\n' \
                     '  write_openpmd_file = "load_bunch.h5"\n'
    rc, out = run(exe, wd, "ldkeep", bextra, pextra="  global%load_only = T\n")
    if rc != 0:
        print(f"FAIL: the keep load exited {rc}:\n{out[-1500:]}")
        return False
    b = bunchfile.read_bunch(wd / "load_bunch.h5")
    rec = load_record(wd, "ldkeep")
    ok = compare(rec, bunch_record(b, 8), 1e-10, "keep mode against the bunch, currents, counts and moments") and ok

    # The dump the load wrote: the bunch's charge, in eight copies per particle.
    q_dump = sum(sl["weight"].sum() for sl in beamio.read_slices(wd / "ldkeep-initial.beam.h5", LAMBDA0, SPACING))
    n_dump = sum(sl["n"] for sl in beamio.read_slices(wd / "ldkeep-initial.beam.h5", LAMBDA0, SPACING))
    dq = abs(q_dump - b["w"].sum()) / b["w"].sum()
    good = dq <= 1e-12 and n_dump == 8 * len(b["x"])
    print(f"  the dump carries the bunch's charge: relative difference {dq:.3e}, "
          f"{n_dump} macroparticles for {len(b['x'])} particles  {'ok' if good else 'FAIL'}")
    ok = ok and good

    # Split weights: a third and two thirds, coincident, must record the same load.
    rc, out = run(exe, wd, "ldsplit", bextra.replace('  write_openpmd_file = "load_bunch.h5"\n', '')
                  + "  resample_split_weights = T\n", pextra="  global%load_only = T\n")
    if rc != 0:
        print(f"FAIL: the split keep load exited {rc}:\n{out[-1500:]}")
        return False
    rec_s = load_record(wd, "ldsplit")
    rec_s[:, 1] = rec[:, 1]                       # Twice the particles by construction.
    ok = compare(rec_s, rec, 1e-11, "split weights against the unsplit load") and ok
    return ok


# The in-cone measurement. The page's deck is a flat 3 kA over 300 slices of 12
# wavelengths on a 255-point grid of 2e-4 m half width, 1024 macroparticles per slice in
# beamlets of 8, and its statistic is the mean over the interior records. Here the current
# is a Gaussian, the bunch's own, at the same peak: sigma_z of 5.5e-9 m, 37 slices of 12
# wavelengths across 8 sigma, in a window with 24 slices of headroom on each side, since
# the field slips 22 slices over the first undulator (266 periods) and every record that
# started on the bunch must still be in the window at its end.

INCONE_TARGET = 0.133e6
INCONE_CUT = 3e-6
INCONE_SAMPLE = 12
SIG_Z = 5.5e-9
I_PEAK = 3000.0
Q_INCONE = I_PEAK * np.sqrt(2 * np.pi) * SIG_Z / bunchfile.C_LIGHT
WINDOW = 8 * SIG_Z + 2 * 24 * INCONE_SAMPLE * LAMBDA0


def incone(exe, wd):
    bextra = f"""  beam_init%n_particle = 1024
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beam_init%sig_z = {SIG_Z:.6e}
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%bunch_charge = {Q_INCONE:.12e}
  beamlet_size = 8
  quiet_start = T
  shot_noise = T
"""
    pextra = (f'  global%track_end = "UND##1"\n  global%source_filter = F\n'
              f'  global%source_filter_angle = {INCONE_CUT:.1e}\n'
              f'  slicing%window_length = {WINDOW:.6e}\n')
    rc, out = run(exe, wd, "ldcone", bextra, pextra=pextra, lat="aramis.bmad", ngrid=256,
                  sample=INCONE_SAMPLE)
    if rc != 0:
        print(f"FAIL: the in-cone run exited {rc}:\n{out[-1500:]}")
        return False
    with h5py.File(wd / "ldcone.stats.h5") as h:
        angle = float(np.ravel(h["field/total/split_angle"][()])[0])
        p_in = np.asarray(h["field/total/power_inside_angle"])
        p_tot = np.asarray(h["field/total/power"])
        cur = np.asarray(h["beam/slice/current"])[0]
    last = np.where(np.isfinite(p_in).any(axis=1))[0][-1]
    edge = max(p_tot[last, 0], p_tot[last, -1]) / p_tot[last].max()
    per_3ka = float(np.nansum(p_in[last]) / np.sum(cur / I_PEAK))
    rel = abs(per_3ka / INCONE_TARGET - 1)
    good = rel <= 0.35 and abs(angle - INCONE_CUT) <= 1e-9 and edge < 1e-3
    print(f"  in-cone startup, Gaussian bunch in sample mode: {per_3ka:.4e} W per 3 kA slice inside "
          f"{angle:.1e} rad at the end of the first undulator, summed over {p_in.shape[1]} records "
          f"(edge records carry {edge:.1e} of the peak), against {INCONE_TARGET:.3e} W: "
          f"relative {rel:.3f} (check <= 0.35)  {'ok' if good else 'FAIL'}")
    return good


def refusals(exe, wd):
    ok = True
    floor_re = re.compile(r"MAX_H \|B\(H\)\|\^2 \* N_LAMBDA =\s*([0-9.E+-]+)")

    # The generator without the quiet start carries the load's own noise already.
    rc, out = run(exe, wd, "ldloud", BUNCH.replace("40000", "512") + "  quiet_start = F\n  shot_noise = T\n",
                  pextra="  global%load_only = T\n")
    m = floor_re.search(out)
    good = rc != 0 and "NOT QUIET ENOUGH" in out and m is not None
    print(f"  refused, noise on the unquiet generator: {'ok' if good else 'MISSED'}"
          + (f", floor |b|^2 N_lambda = {float(m.group(1)):.2e}" if m else ""))
    ok = ok and good

    # A bunch of copies made outside the loader loads with the noise on the copies.
    b = bunchfile.read_bunch(wd / "load_bunch.h5")
    bunchfile.write_bunch(wd / "load_bunch.h5", wd / "load_copies.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING))
    rc, out = run(exe, wd, "ldcopies", '  beam_init%position_file = "load_copies.h5"\n'
                  '  load_mode = "keep"\n  quiet_start = F\n  shot_noise = T\n',
                  pextra="  global%load_only = T\n")
    good = rc == 0 and "particles sharing their transverse coordinates" in out
    fl = re.search(r"worst quiet floor before imposing, \|b\|\^2 \* N_lambda:\s*([0-9.E+-]+)", out)
    print(f"  copies made outside the loader take the noise on the copies: {'ok' if good else 'FAIL'}"
          + (f", floor {float(fl.group(1)):.2e}" if fl else ""))
    ok = ok and good

    # The same copies shifted singly are not quiet, and the message says by how much.
    bunchfile.write_bunch(wd / "load_bunch.h5", wd / "load_shifted.h5",
                          bunchfile.copies(b, 8, LAMBDA0, SPACING, shift=1.0))
    rc, out = run(exe, wd, "ldshift", '  beam_init%position_file = "load_shifted.h5"\n'
                  '  load_mode = "keep"\n  quiet_start = F\n  shot_noise = T\n',
                  pextra="  global%load_only = T\n")
    m = floor_re.search(out)
    good = rc != 0 and "NOT QUIET ENOUGH" in out and m is not None
    print(f"  refused, copies shifted singly: {'ok' if good else 'MISSED'}"
          + (f", floor |b|^2 N_lambda = {float(m.group(1)):.2e}" if m else ""))
    ok = ok and good

    # A slice dump is not a bunch.
    rc, out = run(exe, wd, "lddump", '  beam_init%position_file = "ldkeep-initial.beam.h5"\n'
                  '  load_mode = "keep"\n', pextra="  global%load_only = T\n")
    good = rc != 0 and "BEAM DUMP IN SLICES" in out
    print(f"  refused, a slice dump through beam_init%position_file: {'ok' if good else 'MISSED'}")
    ok = ok and good
    return ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir).resolve()
    exe = str(pathlib.Path(a.exe).resolve())

    ok = keep_exactness(exe, wd)
    ok = refusals(exe, wd) and ok
    ok = incone(exe, wd) and ok
    print("  load-path checks: PASS" if ok else "  load-path checks: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
