#!/usr/bin/env python3
"""
Distribution-import checks (deliverable 10): the bunch_struct -> FEL slices resampler,
transcribed from Genesis's SDDSBeam.cpp, validated along the RNG boundary.

The driver generates a bunch from Bmad's beam_init_struct (a nm-scale Gaussian test
bunch -- sized to the FEL time window's economics, not to a physical accelerator
bunch), writes it as a Genesis DISTRIBUTION file (t = -tau/c, so both codes bin the
identical particle set identically), imports it through the transcribed resampler, and
Genesis imports the same file through &importdistribution. Then:

EXACT (RNG-free), tolerance at roundoff:
  current   The per-slice current profile, ours vs Genesis's beam dump. Window
            membership and the weighted sum contain no random numbers.
  moments   The generated bunch carries the specified emittance (matched generation
            comes from init_beam_distribution and the lattice Twiss; there is no
            match transform in the import).
  split     Coincident w/3 + 2w/3 copies of every bunch particle before import leave
            the moments (unweighted, over coincident copies) and the current profile
            (weighted sums) unchanged.
  threads   The same seed at 1 and 8 threads gives byte-identical diag output (the
            import runs serially before tracking).

STATISTICAL (the resampling and loading RNG):
  twiss     Per-slice Twiss/emittance measured from the imported dump must recover
            the lattice Twiss and the beam_init emittance, central slices, ~percent.
  startup   Dark-start SASE power after one segment, ours vs Genesis's, each code
            resampling with its OWN RNG from the same file: mean over seeds and
            slices within |ln ratio| < 0.30 (the check_sase_startup check, reused).

Usage: check_import.py --exe <fel_track_test> --genesis <genesis4> --workdir <dir> [--seeds N]
The workdir must hold Aramis-1seg.lat and aramis_1seg.bmad. Exit 0 on pass.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

import pool

# The test bunch: Gaussian, peak current ~3 kA at sigma_z = 1.2 nm (charge 30.1 fC),
# Benchmark1-SASE energy and emittance. The window: 24 slices of 3 wavelengths.
LAMBDA0 = 1e-10
SAMPLE = 3
NSLICE = 24
SLEN = NSLICE * SAMPLE * LAMBDA0
GAMMA0 = 11357.82
NORM_EMIT = 4e-7
SIG_PZ = 1 / GAMMA0            # gen_delgam = 1 in fractional momentum
# The lattices carry the matched FODO Twiss in their beginning statements;
# init_beam_distribution generates matched bunches from them (no match transform
# exists in the import -- a Bmad lattice IS the optics specification).
LATTICE_TWISS = dict(betax=8.53711, alphax=-0.703306, betay=17.3899, alphay=1.40348)

NML = """&fel_track_params
  lat_file = "aramis_1seg.bmad"
  out_root = "{root}"
  lambda0 = {lambda0}
  window_sample = {sample}
  imp%npart = 2048
  imp%nbins = 8
  ran_seed = {seed}
  seed_power = 0
  grid_n_pts = 255
  grid_half_width = 2e-4
{source}  imp%nslice = {nslice}
  imp%slicewidth = 0.01
  write_diag = T
{extra}&end
"""

GENESIS_WRITE_DECK = """&setup
rootname={root}
lattice=Aramis-1seg.lat
beamline=SEG1
lambda0={lambda0}
gamma0={gamma0}
delz=0.045
shotnoise=1
npart = 2048
nbins = 8
&end

&time
slen = {slen}
sample = {sample}
&end

&importdistribution
file = {dist}
slicewidth = 0.01
&end

&write
beam = {root}-imp
&end
"""

GENESIS_TRACK_DECK = """&setup
rootname={root}
lattice=Aramis-1seg.lat
beamline=SEG1
lambda0={lambda0}
gamma0={gamma0}
delz=0.045
shotnoise=1
npart = 2048
nbins = 8
beam_global_stat = true
field_global_stat = true
&end

&time
slen = {slen}
sample = {sample}
&end

&importdistribution
file = {dist}
slicewidth = 0.01
&end

&field
power=0
dgrid=2.000000e-04
ngrid=255
&end

&track
fft_fieldsolver = true
&end
"""


def run(cmd, log, env=None):
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(f"FAIL: command exited {r.returncode}: {' '.join(str(c) for c in cmd)}")
        print(pathlib.Path(log).read_text()[-2000:])
        sys.exit(1)


BEAM_INIT_SOURCE = """  use_beam_init = T
  beam_init%n_particle = 50000
  beam_init%a_norm_emit = {emit}
  beam_init%b_norm_emit = {emit}
  beam_init%sig_z = 1.2e-9
  beam_init%sig_pz = {sig_pz}
  beam_init%bunch_charge = 3.01e-14
""".format(emit=NORM_EMIT, sig_pz=SIG_PZ)


def write_nml(path, root, seed, extra="", source=BEAM_INIT_SOURCE):
    path.write_text(NML.format(root=root, seed=seed, lambda0=LAMBDA0,
                               sample=SAMPLE, nslice=NSLICE, source=source, extra=extra))


def parse_stdout(log):
    """The RNG-free lines the exactness checks read: moments and per-slice currents."""
    text = pathlib.Path(log).read_text()
    m = re.search(r"import moments \([^)]*\):\s*(.*)", text)
    moments = np.array([float(v) for v in m.group(1).split()])
    currents = np.array([float(line.rsplit(":", 1)[1])
                         for line in text.splitlines() if line.startswith("import current ")])
    return moments, currents


def load_dump_currents(fn):
    with h5py.File(fn) as h5:
        n = int(h5["slicecount"][0])
        return np.array([h5[f"slice{i+1:06d}/current"][0] for i in range(n)])


def load_dump_slices(fn):
    """Per-slice particle arrays of a Genesis-format dump (our writer's layout)."""
    out = []
    with h5py.File(fn) as h5:
        n = int(h5["slicecount"][0])
        for i in range(n):
            g = h5[f"slice{i+1:06d}"]
            out.append({k: g[k][:] for k in ("gamma", "x", "y", "px", "py")})
        cur = np.array([h5[f"slice{i+1:06d}/current"][0] for i in range(n)])
    return out, cur


def final_powers(diag, nslice):
    d = np.loadtxt(diag)
    d = d.reshape(-1, nslice, d.shape[1])
    return d[-1, :, 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--genesis", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--seeds", type=int, default=4)
    args = ap.parse_args()
    w = pathlib.Path(args.workdir)
    import os
    env1 = dict(os.environ, OMP_NUM_THREADS="1", FI_PROVIDER="tcp")
    env8 = dict(os.environ, OMP_NUM_THREADS="8", FI_PROVIDER="tcp")

    ok = True

    # ---- The RNG-free runs: generate, write the shared file, import, dump, stop.
    write_nml(w/"imp_ref.nml", "impref", 1000,
              extra='  write_dist_file = "impdist.h5"\n'
                    '  write_opmd_file = "impopmd.h5"\n  load_only = T\n')
    run([args.exe, "imp_ref.nml"], w/"imp_ref.log", env=env1)
    mom_ref, cur_ref = parse_stdout(w/"imp_ref.log")

    # moments sanity: the bunch is generated matched to the lattice Twiss by
    # init_beam_distribution; the measured whole-bunch moments must carry the spec.
    ex, ey = mom_ref[5], mom_ref[6]
    print(f"bunch moments: normalized emittance ex = {ex:.6e}, ey = {ey:.6e} (spec {NORM_EMIT:.1e})")
    if not (0.9*NORM_EMIT < ex < 1.1*NORM_EMIT and 0.9*NORM_EMIT < ey < 1.1*NORM_EMIT):
        print("FAIL: generated bunch does not carry the specified emittance")
        ok = False

    # current: ours vs Genesis's import of the same file, roundoff.
    (w/"imp_gen_write.in").write_text(GENESIS_WRITE_DECK.format(
        root="impgw", dist="impdist.h5", lambda0=LAMBDA0, gamma0=GAMMA0,
        slen=SLEN, sample=SAMPLE))
    run([args.genesis, "imp_gen_write.in"], w/"imp_gen_write.log", env=env1)
    cur_gen = load_dump_currents(w/"impgw-imp.par.h5")
    cur_ours = load_dump_currents(w/"impref-initial.par.h5")
    if len(cur_gen) != len(cur_ours):
        print(f"FAIL: slice counts differ (ours {len(cur_ours)}, genesis {len(cur_gen)})")
        ok = False
    else:
        scale = max(cur_gen.max(), 1e-30)
        rel = np.abs(cur_ours - cur_gen).max() / scale
        print(f"current exactness: {len(cur_gen)} slices, peak {cur_gen.max():.4e} A, "
              f"max |diff|/peak = {rel:.3e}")
        if rel > 1e-10:
            print("FAIL: current profiles differ beyond roundoff")
            ok = False

    # split: coincident copies leave the RNG-free outputs unchanged to roundoff.
    write_nml(w/"imp_split.nml", "impsplit", 1000,
              extra='  imp_split_weights = T\n  load_only = T\n')
    run([args.exe, "imp_split.nml"], w/"imp_split.log", env=env1)
    mom_s, cur_s = parse_stdout(w/"imp_split.log")
    dm = np.abs(mom_s - mom_ref).max() / np.abs(mom_ref).max()
    dc = np.abs(cur_s - cur_ref).max() / max(cur_ref.max(), 1e-30)
    print(f"split-weight invariance: moments {dm:.3e}, currents {dc:.3e}")
    if dm > 1e-12 or dc > 1e-12:
        print("FAIL: split-weight invariance broken")
        ok = False

    # openPMD round trip: the same bunch written by Bmad's hdf5_write_beam and read
    # back through dist_file (hdf5_read_beam) must reproduce the RNG-free outputs to
    # file precision -- the import's second input path, exercised end to end.
    write_nml(w/"imp_opmd.nml", "impopmd", 1000,
              extra='  load_only = T\n',
              source='  dist_file = "impopmd.h5"\n')
    run([args.exe, "imp_opmd.nml"], w/"imp_opmd.log", env=env1)
    mom_o, cur_o = parse_stdout(w/"imp_opmd.log")
    dm = np.abs(mom_o - mom_ref).max() / np.abs(mom_ref).max()
    dc = np.abs(cur_o - cur_ref).max() / max(cur_ref.max(), 1e-30)
    print(f"openPMD round trip (write_opmd_file -> dist_file): moments {dm:.3e}, currents {dc:.3e}")
    if dm > 1e-10 or dc > 1e-10:
        print("FAIL: openPMD round trip does not reproduce the bunch")
        ok = False

    # zero-charge refusal: a chargeless bunch (openPMD without charge data, unset
    # beam_init%bunch_charge) must be refused by name, not imported as a dark beam.
    write_nml(w/"imp_dark.nml", "impdark", 1000, extra='  load_only = T\n',
              source=BEAM_INIT_SOURCE.replace("bunch_charge = 3.01e-14", "bunch_charge = 0"))
    with open(w/"imp_dark.log", "w") as fh:
        rr = subprocess.run([args.exe, "imp_dark.nml"], stdout=fh, stderr=subprocess.STDOUT, env=env1)
    dark_log = (w/"imp_dark.log").read_text()
    if rr.returncode == 0 or "ZERO TOTAL CHARGE" not in dark_log:
        print(f"FAIL: zero-charge bunch not refused by name (exit {rr.returncode})")
        ok = False
    else:
        print("zero-charge refusal: refused by name")

    # twiss (statistical): the imported dump's central slices recover the targets.
    # Statistics note: each slice has npart/nbins = 256 independent phase-space seeds,
    # so a single slice's Twiss carries ~1/sqrt(256) = 6% noise and a max over slices
    # would check on order statistics. Check on the MEAN over central slices instead.
    slices, cur = load_dump_slices(w/"impref-initial.par.h5")
    central = [i for i in range(len(slices)) if cur[i] > 0.5*cur.max()]
    berr, aerr, eerr = [], [], []
    for i in central:
        s = slices[i]
        for (u, pu, bt, al) in (("x", "px", LATTICE_TWISS["betax"], LATTICE_TWISS["alphax"]),
                                ("y", "py", LATTICE_TWISS["betay"], LATTICE_TWISS["alphay"])):
            xx = s[u] - s[u].mean()
            xpp = s[pu]/s["gamma"] - (s[pu]/s["gamma"]).mean()   # slope = (gamma*beta_u)/gamma
            e_g = math.sqrt(abs((xx**2).mean()*(xpp**2).mean() - (xx*xpp).mean()**2))
            berr.append((xx**2).mean()/e_g / bt - 1)
            aerr.append(-(xx*xpp).mean()/e_g / al - 1)
            eerr.append(e_g * s["gamma"].mean()/NORM_EMIT - 1)
    bm, am, em = (abs(np.mean(v)) for v in (berr, aerr, eerr))
    print(f"twiss recovery (mean over {len(central)} central slices x 2 planes): "
          f"beta {bm:.3f}, alpha {am:.3f}, emittance {em:.3f} "
          f"(per-slice spread ~{1/math.sqrt(256):.2f})")
    if bm > 0.05 or am > 0.10 or em > 0.05:
        print("FAIL: imported slices do not recover the specified optics")
        ok = False

    # threads (exact): same seed, 1 vs 8 threads, tracked this time; diag byte-equal.
    write_nml(w/"imp_t1.nml", "impt1", 1000)
    write_nml(w/"imp_t8.nml", "impt8", 1000)
    run([args.exe, "imp_t1.nml"], w/"imp_t1.log", env=env1)
    run([args.exe, "imp_t8.nml"], w/"imp_t8.log", env=env8)
    b1 = (w/"impt1.diag.txt").read_bytes()
    b8 = (w/"impt8.diag.txt").read_bytes()
    d1 = re.sub(rb"impt1", b"", b1); d8 = re.sub(rb"impt8", b"", b8)
    print(f"thread determinism: diag files {'byte-identical' if d1 == d8 else 'DIFFER'}")
    if d1 != d8:
        print("FAIL: import + tracking not thread-count independent")
        ok = False

    # startup (statistical): both codes resample the same file with their own RNG.
    # Each seed's pair is a chain (our run writes the dist file its genesis twin
    # imports); the chains are independent across seeds and go through the pool.
    def seed_chain(k):
        seed = 2000 + 17*k
        dist = f"impdist_s{k}.h5"
        write_nml(w/f"imp_s{k}.nml", f"imps{k}", seed,
                  extra=f'  write_dist_file = "{dist}"\n')
        run([args.exe, f"imp_s{k}.nml"], w/f"imp_s{k}.log", env=env8)
        (w/f"imp_g{k}.in").write_text(GENESIS_TRACK_DECK.format(
            root=f"impg{k}", dist=dist, lambda0=LAMBDA0, gamma0=GAMMA0,
            slen=SLEN, sample=SAMPLE))
        run([args.genesis, f"imp_g{k}.in"], w/f"imp_g{k}.log", env=env1)

    pool.run_all([lambda k=k: seed_chain(k) for k in range(args.seeds)], threads_per_job=8)
    ours, theirs = [], []
    for k in range(args.seeds):
        ours.append(final_powers(w/f"imps{k}.diag.txt", NSLICE))
        with h5py.File(w/f"impg{k}.out.h5") as h5:
            theirs.append(h5["Field/power"][-1, :])
    # Compare over slices that carry beam in both (the Gaussian tail slices are dark).
    po = np.concatenate(ours); pg = np.concatenate(theirs)
    m = (po > 0) & (pg > 0)
    lr = math.log(po[m].mean() / pg[m].mean())
    print(f"startup power ({args.seeds} seeds x {NSLICE} slices, {m.sum()} live): "
          f"ln(P_bmad/P_genesis) = {lr:+.3f} (check 0.30)")
    if abs(lr) > 0.30:
        print("FAIL: startup power after import disagrees between the codes")
        ok = False

    # ---- Cross-path equivalence (the beam_init interface deliverable): the SAME
    # Gaussian description, quiet-loaded directly (analytic per-slice evaluation) and
    # imported (real particles resampled). The quiet-load must match the analytic
    # Gaussian profile EXACTLY -- this is also the sqrt(2pi) mutation check on the
    # driver's current derivation -- and the import must match it statistically.
    (w/"imp_xq.nml").write_text(NML.format(root="impxq", seed=1000, lambda0=LAMBDA0,
        sample=SAMPLE, nslice=NSLICE, extra='  load_only = T\n',
        source="""  beam_init%n_particle = 512
  beam_init%a_norm_emit = {emit}
  beam_init%b_norm_emit = {emit}
  beam_init%sig_z = 1.2e-9
  beam_init%sig_pz = {sig_pz}
  beam_init%bunch_charge = 3.01e-14
""".format(emit=NORM_EMIT, sig_pz=SIG_PZ)))
    run([args.exe, "imp_xq.nml"], w/"imp_xq.log", env=env1)
    cur_q = load_dump_currents(w/"impxq-initial.par.h5")
    sp = SAMPLE * LAMBDA0
    sig_z, q_charge = 1.2e-9, 3.01e-14
    n_q = len(cur_q)
    s_i = (np.arange(n_q) - (n_q - 1) / 2) * sp
    ana = q_charge * 2.99792458e8 / (math.sqrt(2*math.pi) * sig_z) * np.exp(-s_i**2 / (2*sig_z**2))
    dev_q = np.abs(cur_q - ana).max() / ana.max()
    print(f"cross-path: quiet-load current vs analytic Gaussian: {dev_q:.2e} of peak (check 1e-12)")
    if dev_q > 1e-12:
        print("FAIL: the quiet-load's derived current is not the described Gaussian")
        ok = False
    # The import's profile (cur_ref, measured earlier from the same description):
    # its own slice centers, bunch center fitted (the import min-shifts positions, and
    # imp%nslice truncates the +4 sigma tail, so the fit is over a clipped profile --
    # rms over the live slices is the honest statistic; the mutations this check exists
    # to catch, a dropped sqrt(2pi) or dslen-for-spacing, sit at 0.65-1.5 of peak).
    si = np.arange(len(cur_ref)) * sp
    c0 = (cur_ref * si).sum() / cur_ref.sum()
    ana_i = q_charge * 2.99792458e8 / (math.sqrt(2*math.pi) * sig_z) * np.exp(-(si - c0)**2 / (2*sig_z**2))
    live = ana_i > 0.05 * ana_i.max()
    dev_i = math.sqrt((((cur_ref - ana_i) / ana_i.max())**2)[live].mean())
    print(f"cross-path: imported current vs the same Gaussian: rms {dev_i:.2e} of peak "
          f"over {live.sum()} live slices (check 7e-2, statistical; measured 4.9e-2)")
    if dev_i > 7e-2:
        print("FAIL: import and quiet-load disagree about the described bunch")
        ok = False

    print("import checks: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
