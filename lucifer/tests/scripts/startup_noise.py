#!/usr/bin/env python3
"""SASE startup against the transverse grid and the macroparticle load (doc/startup-noise.md).

A SASE run started dark has two sources of power in the field. The bunching noise the
quiet start loads radiates, and the FEL instability amplifies the part of that radiation
inside the guided mode. Each beamlet is a transverse point on the deposit grid, and the
paraxial solver radiates it into every transverse wavenumber the grid carries, so the
radiated power depends on the cell size. The measurements here separate the two parts
and follow each against the cell size dx, the particle count and the beamlet count.

The experiments, each a set of runs and one figure:

  a. floor      A window of 300 slices at twelve wavelengths, longer than the line's
                whole slippage, through the Aramis line with field dumps at five
                undulator ends. The far field of each dump by FFT, power inside an
                angular radius of a few mode diffraction angles against power outside it,
                versus z, at four cell sizes and four particle counts. Powers are per
                slice, averaged over the window's interior (slices 80 to 230): the tail
                sixty slices lag by the cooperation length and the head thirty pile up.
  b. beamlets   The beamlet size varied at a fixed particle count, so the beamlet count
                moves alone, and the particle count varied at a fixed beamlet count.
  c. line       The examples' own 96-slice window at three wavelengths through the full
                line at the sweep's cell sizes and loads, beside the examples' numbers.
  d. doubled    The long window, 600 slices, through a line of twice the length, for the
                saturated power and the saturation point against startup.
  e. theory     Ming Xie's gain length and Saldin, Schneidmiller and Yurkov's effective
                shot-noise power and saturation estimates, from the deck's parameters.
  f. genesis    Genesis4 on the harness's SASE tier deck at three grids with both of its
                field solvers, and Lucifer started from each of its dumps.

Lucifer runs on the device by default (global%device = "metal"), with one CPU run at the
sweep's smallest grid as the cross-check. Genesis runs on the CPU. Every run is cached
in the work directory by its deck text, so a rerun costs only what changed. The figures
and a JSON file of every number go to --out. This script is run on demand and never by
the benchmark harness.

Usage:
  startup_noise.py --exe <lucifer> --genesis <genesis4> --pyrepo <openPMD-beamphysics>
                   --examples <lucifer/examples> --latdir <lucifer/tests> --workdir <dir>
                   --out <doc/generated/startup-noise> [--device metal|off] [--only a,b,...]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fieldio  # noqa: E402
from read_stats import Stats  # noqa: E402

C_LIGHT = 2.99792458e8
E_CHARGE = 1.602176634e-19
M_E_EV = 0.51099895e6
EPS0 = 8.8541878128e-12
I_ALFVEN = 17045.0

LAMBDA0 = 1e-10
HALF_WIDTH = 2e-4
CURRENT = 3000.0
SIG_PZ = 8.804506566858e-5
NORM_EMIT = 4e-7
GAMMA0 = 11357.82
AW = 0.84853
LAMBDA_U = 0.015
SEG_LENGTH = 3.99
BETA_A = 8.53711
BETA_B = 17.3899
FILL = SEG_LENGTH / (SEG_LENGTH + 0.44 + 0.08 + 0.24)   # undulator length per meter of line

THETA_CUTS = (2e-6, 3e-6, 5e-6)   # rad, the far-field radii. 3e-6 is the one reported.
GRIDS = (64, 128, 256, 512)
NPARTS = (1024, 4096, 16384, 65536)
LONG = (300, 12)          # (slices, sample) the window longer than the line's slippage
DOUBLED = (600, 12)       # the same for the line twice
EXAMPLE = (96, 3)         # the examples' own window
INTERIOR = slice(80, 230)
DUMP_ELES = ("UND##1", "UND##2", "UND##4", "UND##8", "UND##12")
UND_END_RECORDS = {1: 0, 2: 4, 4: 12, 8: 28, 12: 44}   # record index of each undulator end

DECK = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%comb_ds_save = -1
{extra}/

&fel_beam_init
  beam_init%n_particle = {npart}
  beam_init%bunch_charge = {charge:.12e}
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -{half:.6e}
  beam_init%grid(3)%x_max = {half:.6e}
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  shot_noise = T
  beamlet_size = {beamlet}
/

&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 0
  wavefront_init%grid_n_pts = {ngrid}
  wavefront_init%grid_half_width = 2e-4
  wavefront_init%window_length = {slen:.6e}
  wavefront_init%window_sample = {sample}
/
"""

GENESIS_DECK = """&setup
rootname={root}
lattice=Aramis.lat
beamline=ARAMIS
lambda0=1e-10
gamma0=11357.82
delz=0.045000
shotnoise=1
npart = 2048
nbins = 8
beam_global_stat = true
field_global_stat = true
&end

&lattice
zmatch=9.5
&end

&time
slen = 9.6e-9
sample = 3
&end

&field
power=0
dgrid=2.000000e-04
ngrid={ngrid}
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
&end

&write
field = {root}-initial
beam = {root}-initial
&end

&track
fft_fieldsolver = {fft}
&end
"""

LUCIFER_FROM_GENESIS = """&fel_params
  lat_file = "aramis.bmad"
  global%out_root = "{root}"
  global%interlude_model = "genesis"
  global%transport_model = "genesis"
  global%write_diag = T
/
&fel_beam_init
  beam_file = "{gen}-initial.beam.h5"
  beamlet_size = 8
/
&fel_wavefront_init
  field_file = "{gen}-initial.wf.h5"
  wavefront_init%lambda0 = 1e-10
  wavefront_init%window_sample = 3
/
"""


# ---------------------------------------------------------------------------
# Running and caching

def cell_size(ngrid):
    return 2 * HALF_WIDTH / (ngrid - 1)


def deck_text(lat, root, ngrid, npart, beamlet, window, device, dumps=()):
    nslice, sample = window
    slen = nslice * sample * LAMBDA0
    charge = CURRENT * slen / C_LIGHT
    extra = ""
    if device != "off":
        extra += f'  global%device = "{device}"\n'
    if dumps:
        extra += "  global%dump_field_at = " + ", ".join(f'"{d}"' for d in dumps) + "\n"
    return DECK.format(lat=lat, root=root, ngrid=ngrid, npart=npart, beamlet=beamlet,
                       charge=charge, half=slen / 2, slen=slen, sample=sample, extra=extra)


class Runner:
    """Runs decks in a work directory, skipping any whose text and outputs are already
    there. The stamp is the deck's hash, so an edited deck reruns and an unchanged one
    does not."""

    def __init__(self, exe, workdir, threads):
        self.exe = str(exe)
        self.wd = pathlib.Path(workdir)
        self.wd.mkdir(parents=True, exist_ok=True)
        self.threads = str(threads)

    def run(self, name, text, exe=None, env_threads=None):
        deck = self.wd / f"{name}.in"
        stamp = self.wd / f"{name}.done"
        digest = hashlib.sha256(text.encode()).hexdigest()
        if stamp.exists() and stamp.read_text().strip() == digest:
            return
        deck.write_text(text)
        env = {"OMP_NUM_THREADS": env_threads or self.threads, "PATH": os.environ["PATH"],
               "HOME": os.environ.get("HOME", "")}
        with open(self.wd / f"{name}.log", "w") as log:
            r = subprocess.run([exe or self.exe, deck.name], cwd=self.wd, stdout=log,
                               stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            print(f"FAIL: {name} exited {r.returncode}, see {self.wd / (name + '.log')}")
            sys.exit(1)
        stamp.write_text(digest)
        print(f"  ran {name}")


# ---------------------------------------------------------------------------
# Reading

def element_end_power(wd, root):
    """(z, per-slice power, per-slice |bunching|) at every record the stats file kept,
    (nrec,), (nrec, nslice) and (nrec, nslice)."""
    with Stats(wd / f"{root}.stats.h5") as st:
        z = np.asarray(st.s, dtype=float)
        P = np.asarray(st["field/total/power"], dtype=float)
        b = np.asarray(st["beam/slice/bunching"], dtype=float)
    return z, P, b


def interior(P):
    """Mean per-slice power over the window's interior, (nrec,) from (nrec, nslice)."""
    return P[:, INTERIOR].mean(axis=1)


def far_field_split(path, cuts):
    """Power inside and outside each angular radius, summed over slices, from one field
    dump. The far field is the FFT of the field over the grid, and the angle of a
    transverse wavenumber is k_perp / k. Parseval keeps the total equal to the
    near-field power, which is checked here to 1e-9."""
    f = fieldio.read_field(path)
    u, dx, dy = f["u"], f["dx"], f["dy"]
    nslice, ny, nx = u.shape
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=dy)
    k = 2 * np.pi / f["wavelength"]
    theta = np.sqrt(kx[None, :] ** 2 + ky[:, None] ** 2) / k
    spec = np.abs(np.fft.fft2(u, axes=(1, 2))) ** 2 / (nx * ny)
    scale = dx * dy / (2 * fieldio.MU0_C)
    total = spec.sum(axis=(1, 2)) * scale
    near = fieldio.field_power(u, dx, dy)
    if np.abs(total.sum() - near.sum()) > 1e-9 * max(near.sum(), 1e-300):
        raise RuntimeError(f"{path}: Parseval fails, {total.sum()} vs {near.sum()}")
    out = {"total": float(total[INTERIOR].mean()), "total_window": float(total.sum()),
           "z": float(f["s_position"]), "theta_nyquist": float(np.pi / dx / k)}
    for c in cuts:
        inside = (spec * (theta <= c)).sum(axis=(1, 2)) * scale
        out[f"in_{c:.0e}"] = float(inside[INTERIOR].mean())
        out[f"out_{c:.0e}"] = float((total - inside)[INTERIOR].mean())
    return out


def dumps_of(wd, root, nele=len(DUMP_ELES)):
    """The field dumps of a run in element order, whatever element indices they carry."""
    files = sorted(wd.glob(f"{root}-at*-UND.wf.h5"),
                   key=lambda p: int(p.name.split("-at")[1].split("-")[0]))
    if len(files) != nele:
        raise RuntimeError(f"{root}: expected {nele} field dumps, found {len(files)}")
    return files


def analyze_dumps(wd, root, cache):
    """Far-field splits of every dump of a run, cached as JSON so the dumps can go."""
    if cache.exists():
        return json.loads(cache.read_text())
    rows = [far_field_split(p, THETA_CUTS) for p in dumps_of(wd, root)]
    cache.write_text(json.dumps(rows))
    for p in dumps_of(wd, root):
        p.unlink()
    return rows


# ---------------------------------------------------------------------------
# Theory

def pierce():
    """The 1D Pierce parameter for the deck, Genesis's own form with the rms aw and the
    helical coupling fc = 1, and the transverse size from the mean matched beta."""
    beta = 0.5 * (BETA_A + BETA_B)
    emit = NORM_EMIT / GAMMA0
    sigma = math.sqrt(emit * beta)
    ku = 2 * math.pi / LAMBDA_U
    rho3 = (CURRENT / I_ALFVEN) * AW ** 2 / (8 * GAMMA0 ** 3 * sigma ** 2 * ku ** 2)
    return rho3 ** (1 / 3), sigma, emit, beta


def ming_xie():
    """Ming Xie's fit for the 3D power gain length (Nucl. Instrum. Methods A 445, 59
    (2000)), and the saturation power estimate 1.6 rho (L1D/Lg)^2 P_beam that goes with it."""
    rho, sigma, emit, beta = pierce()
    l1d = LAMBDA_U / (4 * math.pi * math.sqrt(3) * rho)
    lr = 4 * math.pi * sigma ** 2 / LAMBDA0
    eta_d = l1d / lr
    eta_e = (l1d / beta) * (4 * math.pi * emit / LAMBDA0)
    eta_g = 4 * math.pi * (l1d / LAMBDA_U) * SIG_PZ
    a = [None, 0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7,
         1.9, 1140.0, 2.2, 2.9, 3.2]
    lam = (a[1] * eta_d ** a[2] + a[3] * eta_e ** a[4] + a[5] * eta_g ** a[6]
           + a[7] * eta_e ** a[8] * eta_g ** a[9] + a[10] * eta_d ** a[11] * eta_g ** a[12]
           + a[13] * eta_d ** a[14] * eta_e ** a[15]
           + a[16] * eta_d ** a[17] * eta_e ** a[18] * eta_g ** a[19])
    lg = l1d * (1 + lam)
    p_beam = CURRENT * GAMMA0 * M_E_EV
    p_sat = 1.6 * rho * (l1d / lg) ** 2 * p_beam
    return {"rho": rho, "sigma_x": sigma, "L1D": l1d, "eta_d": eta_d, "eta_eps": eta_e,
            "eta_gamma": eta_g, "Lg": lg, "P_beam": p_beam, "P_sat_xie": p_sat}


def ssy(lg):
    """Saldin, Schneidmiller and Yurkov (New J. Phys. 12, 035010 (2010)). The effective
    shot-noise power that seeds the exponential regime, P(z) = P_eff exp(z/Lg) / 9 per
    slice, in the 1D form 6 sqrt(pi) rho^2 P_beam / (N_lambda sqrt(ln(N_lambda/rho))).
    The number of electrons per coherence volume N_c = I N_g lambda / (e c) with N_g the
    field gain length in periods. The saturation length 0.6 L_g,field ln N_c and the
    efficiency 0.17 / eps_hat with eps_hat = 2 pi eps / lambda, their eq. (18). The
    incoherent undulator power of the beam in the central cone, their eq. (1), for a
    helical undulator with A_JJ = 1."""
    rho, sigma, emit, beta = pierce()
    p_beam = CURRENT * GAMMA0 * M_E_EV
    n_lambda = CURRENT * LAMBDA0 / (E_CHARGE * C_LIGHT)
    p_eff = 6 * math.sqrt(math.pi) * rho ** 2 * p_beam / (n_lambda * math.sqrt(math.log(n_lambda / rho)))
    lg_field = 2 * lg
    n_c = CURRENT * (lg_field / LAMBDA_U) * LAMBDA0 / (E_CHARGE * C_LIGHT)
    l_sat = 0.6 * lg_field * math.log(n_c)
    eps_hat = 2 * math.pi * emit / LAMBDA0
    rho_bar = LAMBDA_U / (4 * math.pi * math.sqrt(3) * lg)
    p_sat = 0.17 / eps_hat * rho_bar * p_beam
    k2 = AW ** 2
    w_incoh = math.pi * E_CHARGE * CURRENT / (EPS0 * LAMBDA0) * k2 / (1 + k2)
    theta_cone_seg = math.sqrt(1 + k2) / (GAMMA0 * math.sqrt(SEG_LENGTH / LAMBDA_U))
    return {"N_lambda": n_lambda, "P_eff": p_eff, "N_c": n_c, "L_sat_ssy": l_sat,
            "eps_hat": eps_hat, "P_sat_ssy": p_sat, "W_incoh_cone": w_incoh,
            "theta_cone_one_segment": theta_cone_seg,
            "W_incoh_in_3urad_one_segment": w_incoh * (3e-6 / theta_cone_seg) ** 2,
            "theta_mode": LAMBDA0 / (2 * math.pi * sigma)}


# ---------------------------------------------------------------------------
# Experiments

def exp_floor(rn, args, lat, results):
    print("== a. floor and gain, the long window, dumps at five undulator ends ==")
    res = {}
    for ng in GRIDS:
        for npart in NPARTS:
            root = f"floor_g{ng}_n{npart}"
            rn.run(root, deck_text(lat, root, ng, npart, 8, LONG, args.device, DUMP_ELES))
            rows = analyze_dumps(rn.wd, root, rn.wd / f"{root}.farfield.json")
            z, P, b = element_end_power(rn.wd, root)
            res[root] = {"ngrid": ng, "dx": cell_size(ng), "npart": npart, "beamlets": npart // 8,
                         "dumps": rows, "z": z.tolist(), "P": interior(P).tolist(),
                         "b": interior(b).tolist()}
            if root == "floor_g256_n4096":
                res[root]["profile"] = {str(k): P[i].tolist() for k, i in UND_END_RECORDS.items()}
    # The CPU cross-check at the smallest grid and load.
    root = "floor_g64_n1024_cpu"
    rn.run(root, deck_text(lat, root, 64, 1024, 8, LONG, "off", DUMP_ELES), env_threads=args.cpu_threads)
    rows = analyze_dumps(rn.wd, root, rn.wd / f"{root}.farfield.json")
    dev = res["floor_g64_n1024"]["dumps"]
    rel = max(abs(a["total"] - b["total"]) / b["total"] for a, b in zip(dev, rows))
    rel_in = max(abs(a["in_3e-06"] - b["in_3e-06"]) / b["in_3e-06"] for a, b in zip(dev, rows))
    res["cpu_cross_check"] = {"worst_total_rel": rel, "worst_inside_rel": rel_in, "dumps": rows}
    print(f"  device vs CPU at grid 64, 1024 particles: total {rel:.2e}, inside 3 urad {rel_in:.2e}")
    results["floor"] = res


def exp_beamlets(rn, args, lat, results):
    print("== b. beamlets against particles, grid 256 ==")
    res = {}
    cases = [(4096, 4), (4096, 8), (4096, 16), (2048, 4), (8192, 16)]
    for npart, bl in cases:
        root = f"beamlet_n{npart}_b{bl}"
        if (npart, bl) == (4096, 8):
            src = "floor_g256_n4096"
            res[root] = {"npart": npart, "beamlet": bl, "beamlets": npart // bl,
                         "dumps": results["floor"][src]["dumps"], "z": results["floor"][src]["z"],
                         "P": results["floor"][src]["P"], "b": results["floor"][src]["b"]}
            continue
        rn.run(root, deck_text(lat, root, 256, npart, bl, LONG, args.device, DUMP_ELES))
        rows = analyze_dumps(rn.wd, root, rn.wd / f"{root}.farfield.json")
        z, P, b = element_end_power(rn.wd, root)
        res[root] = {"npart": npart, "beamlet": bl, "beamlets": npart // bl, "dumps": rows,
                     "z": z.tolist(), "P": interior(P).tolist(), "b": interior(b).tolist()}
    results["beamlets"] = res


def saturation_point(z, P):
    """The z of the first record where the power's growth rate over one FODO cell falls
    under a tenth of its peak, and the power there. NaN if the run never turned over."""
    lnP = np.log(np.maximum(P, 1e-300))
    rate = np.gradient(lnP, z)
    i0 = int(np.argmax(rate))
    for i in range(i0, len(z)):
        if rate[i] < 0.1 * rate[i0]:
            return float(z[i]), float(P[i])
    return float("nan"), float("nan")


def exp_line(rn, args, lat, lat2, results):
    print("== c. the examples' window through the line, d. the long window through the line twice ==")
    res = {}
    for tag, lt, window in (("line", lat, EXAMPLE), ("doubled", lat2, DOUBLED)):
        for ng in (128, 256):
            for npart in (1024, 65536):
                root = f"{tag}_g{ng}_n{npart}"
                rn.run(root, deck_text(lt, root, ng, npart, 8, window, args.device))
                z, P, b = element_end_power(rn.wd, root)
                curve = P.sum(axis=1) if tag == "line" else interior(P)
                zs, ps = saturation_point(z, curve)
                res[root] = {"ngrid": ng, "dx": cell_size(ng), "npart": npart, "z": z.tolist(),
                             "P": curve.tolist(), "b": interior(b).tolist(), "P_exit": float(curve[-1]),
                             "z_sat": zs, "P_sat": ps}
                what = "window power" if tag == "line" else "interior slice power"
                print(f"  {root}: exit {what} {curve[-1]:.3e} W, saturation at z = {zs:.1f} m with {ps:.3e} W")
    results["line"] = res


def exp_genesis(rn, args, results):
    print("== f. Genesis4 on the SASE tier deck at three grids, both solvers ==")
    latdir = pathlib.Path(args.latdir)
    for name in ("genesis4/Aramis.lat", "bmad/aramis.bmad"):
        (rn.wd / pathlib.Path(name).name).write_bytes((latdir / name).read_bytes())
    convert = pathlib.Path(__file__).with_name("convert_genesis.py")
    cases = [(ng, fft) for ng in (151, 255, 511) for fft in (True, False)]

    def one(case):
        ng, fft = case
        root = f"G{ng}{'fft' if fft else 'adi'}"
        text = GENESIS_DECK.format(root=root, ngrid=ng, fft="true" if fft else "false")
        rn.run(root, text, exe=args.genesis, env_threads="1")
        return root

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(cases)) as ex:
        roots = list(ex.map(one, cases))
    res = {}
    for (ng, fft), root in zip(cases, roots):
        with h5py.File(rn.wd / f"{root}.out.h5") as h:
            gP = h["Field/power"][...].sum(axis=1)
        entry = {"ngrid": ng, "fft": fft, "genesis_exit": float(gP[-1]), "genesis_max": float(gP.max())}
        if fft:
            for pair in ((f"{root}-initial.par.h5", f"{root}-initial.beam.h5"),
                         (f"{root}-initial.fld.h5", f"{root}-initial.wf.h5")):
                if not (rn.wd / pair[1]).exists():
                    subprocess.run([sys.executable, str(convert), "to-openpmd", *pair,
                                    "--pyrepo", args.pyrepo], cwd=rn.wd, check=True,
                                   stdout=subprocess.DEVNULL)
            lroot = f"L{ng}"
            rn.run(lroot, LUCIFER_FROM_GENESIS.format(root=lroot, gen=root), env_threads=args.cpu_threads)
            d = np.loadtxt(rn.wd / f"{lroot}.diag.txt", comments="#")
            zs = np.unique(d[:, 0])
            lP = np.array([d[d[:, 0] == zz, 2].sum() for zz in zs])
            m = min(len(lP), len(gP))
            rel = np.abs(lP[:m] - gP[:m]) / np.maximum(gP[:m], 1e-300)
            entry.update({"lucifer_exit": float(lP[-1]), "exit_rel": float(abs(lP[-1] - gP[-1]) / gP[-1]),
                          "max_rel_over_records": float(rel[1:].max())})
        res[root] = entry
        print(f"  {root}: Genesis exit {entry['genesis_exit']:.3e} W" +
              (f", Lucifer {entry['lucifer_exit']:.3e} W, rel {entry['exit_rel']:.1e}, "
               f"max over records {entry['max_rel_over_records']:.1e}" if fft else ""))
    results["genesis"] = res


def rounded(obj, digits=5):
    """The same structure with every float at `digits` significant figures."""
    if isinstance(obj, float):
        return float(f"{obj:.{digits}g}") if math.isfinite(obj) else obj
    if isinstance(obj, dict):
        return {k: rounded(v, digits) for k, v in obj.items()}
    if isinstance(obj, list):
        return [rounded(v, digits) for v in obj]
    return obj


# ---------------------------------------------------------------------------
# Figures

def figures(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
    cut = "3e-06"
    th = results["theory"]

    # 1. The floor and the gain, versus z.
    fl = results["floor"]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
    for ng in GRIDS:
        r = fl[f"floor_g{ng}_n1024"]
        z = [d["z"] for d in r["dumps"]]
        ax[0, 0].semilogy(z, [d[f"out_{cut}"] for d in r["dumps"]], "o-", label=f"{r['dx']*1e6:.2f} µm")
        ax[0, 1].semilogy(z, [d[f"in_{cut}"] for d in r["dumps"]], "o-", label=f"{r['dx']*1e6:.2f} µm")
    for npart in NPARTS:
        r = fl[f"floor_g256_n{npart}"]
        z = [d["z"] for d in r["dumps"]]
        ax[1, 0].semilogy(z, [d[f"out_{cut}"] for d in r["dumps"]], "o-", label=f"{npart} particles")
        ax[1, 1].semilogy(z, [d[f"in_{cut}"] for d in r["dumps"]], "o-", label=f"{npart} particles")
    zz = np.linspace(0, 57, 50)
    for a in (ax[0, 1], ax[1, 1]):
        a.semilogy(zz, np.minimum(th["P_eff"] / 9 * np.exp(zz * FILL / th["Lg"]), th["P_sat_xie"]), "k--", lw=1,
                   label=r"$P_\mathrm{eff}\,e^{z/L_g}/9$, capped at $P_\mathrm{sat}$ (Xie)")
    ax[0, 0].set_title("1024 particles per slice, cell size varied")
    ax[0, 1].set_title("1024 particles per slice, cell size varied")
    ax[1, 0].set_title("1.57 µm cells, particle count varied")
    ax[1, 1].set_title("1.57 µm cells, particle count varied")
    for a in ax.flat:
        a.legend(fontsize=8)
    for a in ax[:, 0]:
        a.set_ylabel("power per slice outside 3 µrad (W)")
    for a in ax[:, 1]:
        a.set_ylabel("power per slice inside 3 µrad (W)")
    for a in ax[1]:
        a.set_xlabel("z (m)")
    fig.tight_layout()
    fig.savefig(out / "power-outside-and-inside-the-mode.png", dpi=150)
    plt.close(fig)

    # 2. The floor against the cell size, the gain against the load, at the first dump.
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    dxs = [fl[f"floor_g{ng}_n1024"]["dx"] for ng in GRIDS]
    for npart, mk in zip(NPARTS, "osd^"):
        ax[0].loglog(np.array(dxs) * 1e6, [fl[f"floor_g{ng}_n{npart}"]["dumps"][0][f"out_{cut}"] for ng in GRIDS],
                     mk, label=f"{npart} particles")
    ref = fl["floor_g256_n1024"]["dumps"][0][f"out_{cut}"]
    ax[0].loglog(np.array(dxs) * 1e6, ref * (cell_size(256) / np.array(dxs)) ** 2, "k--", lw=1, label=r"$1/dx^2$")
    ax[0].loglog(np.array(dxs) * 1e6, ref * (cell_size(256) / np.array(dxs)), "k:", lw=1, label=r"$1/dx$")
    ax[0].set_xlabel("cell size (µm)")
    ax[0].set_ylabel("power per slice outside 3 µrad at z = 3.99 m (W)")
    ax[0].legend(fontsize=8)
    for ng, mk in zip(GRIDS, "osd^"):
        ax[1].loglog(NPARTS, [fl[f"floor_g{ng}_n{n}"]["dumps"][2][f"in_{cut}"] for n in NPARTS], mk + "-",
                     label=f"{cell_size(ng)*1e6:.2f} µm, inside")
        ax[1].loglog(NPARTS, [fl[f"floor_g{ng}_n{n}"]["dumps"][2][f"out_{cut}"] for n in NPARTS], mk + "--",
                     label=f"{cell_size(ng)*1e6:.2f} µm, outside")
    ax[1].set_xlabel("particles per slice")
    ax[1].set_ylabel("power per slice at z = 18.2 m (W)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "wide-angle-vs-cell-size-and-particle-count.png", dpi=150)
    plt.close(fig)

    # 3. Beamlets against particles.
    bl = results["beamlets"]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
    for root, r in bl.items():
        z = [d["z"] for d in r["dumps"]]
        col = 0 if r["npart"] == 4096 else 1
        lab = f"{r['npart']} particles, beamlet {r['beamlet']}, {r['beamlets']} beamlets"
        ax[0, col].semilogy(z, [d[f"in_{cut}"] for d in r["dumps"]], "o-", label=lab)
        ax[1, col].semilogy(z, [d[f"out_{cut}"] for d in r["dumps"]], "o-", label=lab)
    ax[0, 0].set_title("4096 particles per slice, beamlet size varied")
    ax[0, 1].set_title("512 beamlets, particle count varied")
    ax[1, 0].set_title("4096 particles per slice, beamlet size varied")
    ax[1, 1].set_title("512 beamlets, particle count varied")
    for a in ax.flat:
        a.legend(fontsize=8)
    for a in ax[0]:
        a.set_ylabel("power per slice inside 3 µrad (W)")
    for a in ax[1]:
        a.set_ylabel("power per slice outside 3 µrad (W)")
        a.set_xlabel("z (m)")
    fig.tight_layout()
    fig.savefig(out / "beamlets-vs-particles.png", dpi=150)
    plt.close(fig)

    # 4. The full and doubled lines.
    ln = results["line"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for a, tag in zip(ax, ("line", "doubled")):
        for ng in (128, 256):
            for npart in (1024, 65536):
                r = ln[f"{tag}_g{ng}_n{npart}"]
                a.semilogy(r["z"], r["P"], "-", label=f"{r['dx']*1e6:.2f} µm, {npart} particles")
        if tag == "doubled":
            a.axhline(th["P_sat_xie"], color="k", lw=0.8, ls="--", label=r"$P_\mathrm{sat}$, Xie")
            a.axhline(th["P_sat_ssy"], color="k", lw=0.8, ls=":", label=r"$P_\mathrm{sat}$, Saldin, Schneidmiller and Yurkov")
            a.axvline(th["L_sat_ssy"], color="gray", lw=0.8)
        a.set_xlabel("z (m)")
        a.legend(fontsize=8)
    ax[0].set_ylabel("window power (W)")
    ax[1].set_ylabel("power per interior slice (W)")
    ax[0].set_title("96 slices at three wavelengths, the line")
    ax[1].set_title("600 slices at twelve wavelengths, the line twice")
    fig.tight_layout()
    fig.savefig(out / "line-and-doubled.png", dpi=150)
    plt.close(fig)

    # 5. The window profile: per-slice power against the slice index at five undulator ends.
    r = fl["floor_g256_n4096"]
    fig, ax = plt.subplots(figsize=(8, 4))
    for k, prof in r["profile"].items():
        ax.semilogy(prof, lw=1, label=f"undulator {k} end")
    ax.axvspan(INTERIOR.start, INTERIOR.stop, color="gray", alpha=0.15, label="interior")
    ax.set_xlabel("slice index, tail to head")
    ax.set_ylabel("power per slice (W)")
    ax.set_title("1.57 µm cells, 4096 particles per slice")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "window-profile.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--genesis", required=True)
    ap.add_argument("--pyrepo", required=True)
    ap.add_argument("--examples", required=True, help="lucifer/examples, for aramis.bmad")
    ap.add_argument("--latdir", required=True, help="lucifer/tests, for the Genesis lattice and deck")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", default="metal")
    ap.add_argument("--cpu-threads", default="12")
    ap.add_argument("--only", default="a,b,c,d,e,f")
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rn = Runner(args.exe, args.workdir, "4")
    lat = "aramis.bmad"
    (rn.wd / lat).write_bytes((pathlib.Path(args.examples) / lat).read_bytes())
    lat2 = "aramis_twice.bmad"
    (rn.wd / lat2).write_text("call, file = aramis.bmad\nARAMIS2: line = (ARAMIS, ARAMIS)\nuse, ARAMIS2\n")

    results_file = out / "startup-noise.json"
    results = json.loads(results_file.read_text()) if results_file.exists() else {}
    want = set(args.only.split(","))

    th = ming_xie()
    th.update(ssy(th["Lg"]))
    results["theory"] = th
    print("== e. theory for the deck ==")
    for k, v in th.items():
        print(f"  {k:24s} {v:.4e}")

    if "a" in want:
        exp_floor(rn, args, lat, results)
    if "b" in want:
        exp_beamlets(rn, args, lat, results)
    if "c" in want or "d" in want:
        exp_line(rn, args, lat, lat2, results)
    if "f" in want:
        exp_genesis(rn, args, results)
    results_file.write_text(json.dumps(rounded(results), indent=1))
    if all(k in results for k in ("floor", "beamlets", "line")):
        figures(results, out)
    print(f"wrote {results_file}")


if __name__ == "__main__":
    main()
