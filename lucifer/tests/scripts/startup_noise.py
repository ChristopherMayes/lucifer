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
  g. filter     The source filter (fel-physics.md sec-source-filter) at grid 256 and two
                loads, with its sigmoid's edge at the central cone and at the 3 urad cut,
                against the same runs with the filter off. What the filter does to the
                power inside the mode is what this measures.

Lucifer runs on the device by default (global%device = "metal"), with one CPU run at the
sweep's smallest grid as the cross-check. Genesis runs on the CPU. Every run is cached
in the work directory by its deck text, so a rerun costs only what changed. The figures
and a JSON file of every number go to --out. This script is run on demand and never by
the benchmark harness.

Usage:
  startup_noise.py --exe <lucifer> --genesis <genesis4> --pyrepo <openPMD-beamphysics>
                   --examples <lucifer/examples> --latdir <lucifer/tests> --workdir <dir>
                   --out <doc/generated/startup-noise> [--device metal|off] [--only a,b,...]
                   [--machine aramis|flash|flash1]

The second and third machines take experiments e and g: the line of experiments c and d
is the first machine's own, and the sweep of a and b is priced against what its answer
adds. Their results and figures carry the machine's name.
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
import time

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

# Four machines. The Aramis benchmark is the page's own. The second is a probe at FLASH's
# scales, a hundred times the wavelength on a beam five times the size, whose z_R/L_g is
# 2.6 times smaller. The filter's default is built from two angles whose ratio goes as
# 1/sqrt(z_R/L_g), so a second machine is what says whether the default is general or a
# property of the first (fel-physics.md sec-source-filter). The third and fourth are real
# machines at their published parameters, FLASH1 at DESY at 13.7 nm (examples/flash1) and
# LCLS at SLAC at 1.5 Angstrom (examples/lcls). Both are planar where the other two are
# helical, so their coupling factor is not one.
#
# coupling is the undulator's own coupling to the resonant wavelength, one for a helical
# device and the Bessel factor JJ for a planar one. It multiplies aw wherever the
# interaction reads it and never where the resonance does.
#
# dumps and und_end_record place the far-field dumps: five undulator ends, and the record
# index of each in a stats file written with comb_ds_save = -1, where a record is an
# element end. window is the long window of experiments a, b and g, (slices, wavelengths
# per slice), which must be longer than the line's whole slippage. At 1.5 Angstrom the
# slippage is 3762 wavelengths, so the first machine's twelve per slice would not reach it.

MACHINES = {
    "aramis": dict(lat="aramis.bmad", lat2="aramis_x2.bmad", lambda0=1e-10, half_width=2e-4,
                   current=3000.0, sig_pz=8.804506566858e-5, norm_emit=4e-7, gamma0=11357.82,
                   aw=0.84853, coupling=1.0, lambda_u=0.015, seg_length=3.99,
                   beta_a=8.53711, beta_b=17.3899,
                   gaps=(0.44, 0.08, 0.24), filter_grid=256, cuts=(2e-6, 3e-6, 5e-6),
                   window=(300, 12),
                   dumps=("UND##1", "UND##2", "UND##4", "UND##8", "UND##12"),
                   und_end_record={1: 0, 2: 4, 4: 12, 8: 28, 12: 44}),
    "flash": dict(lat="flash.bmad", lat2=None, lambda0=1e-8, half_width=1e-3,
                  current=1500.0, sig_pz=1.0e-4, norm_emit=1.5e-6, gamma0=1572.0,
                  aw=0.9, coupling=1.0, lambda_u=0.0273, seg_length=4.5045,
                  beta_a=7.02965, beta_b=13.39578,
                  gaps=(0.6, 0.1, 0.3), filter_grid=128, cuts=(4e-5, 6.45e-5, 1e-4),
                  window=(300, 12),
                  dumps=("UND##1", "UND##2", "UND##4", "UND##8", "UND##12"),
                  und_end_record={1: 0, 2: 4, 4: 12, 8: 28, 12: 44}),
    "flash1": dict(lat="flash1.bmad", lat2=None, lambda0=1.37e-8, half_width=1e-3,
                   current=2500.0, sig_pz=6.0e-4, norm_emit=1.5e-6, gamma0=1308.2103,
                   aw=0.847162, coupling=0.885232, lambda_u=0.0273, seg_length=4.5045,
                   beta_a=10.1189, beta_b=9.8031,
                   gaps=(0.25, 0.10, 0.25), filter_grid=128, cuts=(5e-5, 8.1e-5, 1.2e-4),
                   window=(300, 12),
                   dumps=("UND##1", "UND##2", "UND##3", "UND##4", "UND##6"),
                   und_end_record={1: 0, 2: 6, 3: 12, 4: 18, 6: 30}),
    # LCLS: the grid derives to 127 points over 191 um, so the sweep runs at the device's
    # 128. Dumps at the ends of undulators 1, 6, 12, 20 and 30 span the startup, the
    # exponential regime and saturation, which on this line sits near 60 m. A cell is six
    # segments and twelve breaks, so undulator n ends at record 4(n-1) + floor((n-1)/3).
    "lcls": dict(lat="lcls.bmad", lat2=None, lambda0=1.50992e-10, half_width=1.914e-4,
                 current=3000.0, sig_pz=1.3e-4, norm_emit=4.0e-7, gamma0=26614.536,
                 aw=2.475918, coupling=0.744325, lambda_u=0.03, seg_length=3.42,
                 beta_a=26.830, beta_b=33.341,
                 gaps=(0.135, 0.20, 0.135), filter_grid=128, cuts=(3e-6, 4.5e-6, 7e-6),
                 window=(300, 68),
                 dumps=("UND##1", "UND##6", "UND##12", "UND##20", "UND##30"),
                 und_end_record={1: 0, 6: 20, 12: 44, 20: 76, 30: 116}),
}

LAMBDA0 = 1e-10
FC = 1.0                  # undulator coupling: 1 helical, the Bessel factor JJ planar
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


def select_machine(name):
    """Point the module's deck and theory constants at one machine."""
    m = MACHINES[name]
    g = globals()
    g["LAMBDA0"] = m["lambda0"];      g["HALF_WIDTH"] = m["half_width"]
    g["CURRENT"] = m["current"];      g["SIG_PZ"] = m["sig_pz"]
    g["NORM_EMIT"] = m["norm_emit"];  g["GAMMA0"] = m["gamma0"]
    g["AW"] = m["aw"];                g["LAMBDA_U"] = m["lambda_u"]
    g["SEG_LENGTH"] = m["seg_length"]
    g["BETA_A"] = m["beta_a"];        g["BETA_B"] = m["beta_b"]
    g["FC"] = m["coupling"]
    g["FILL"] = m["seg_length"] / (m["seg_length"] + sum(m["gaps"]))
    g["DUMP_ELES"] = m["dumps"]
    g["UND_END_RECORDS"] = m["und_end_record"]
    g["LONG"] = m["window"]
    # The middle cut is the one reported: on the first machine the 3 urad of the page, on
    # the second the edge the code derives for it, 65 urad, since a cut inside the mode
    # would report the mode's own shape rather than the split.
    g["THETA_CUTS"] = m["cuts"]
    g["REPORT_CUT"] = f"{m['cuts'][1]:.0e}"
    return m

THETA_CUTS = (2e-6, 3e-6, 5e-6)   # rad, the far-field radii. 3e-6 is the one reported.
REPORT_CUT = "3e-06"              # the key of the reported radius, f"{cut:.0e}"


def cut_key(kind):
    """The far-field key for the reported radius, 'in' or 'out'."""
    return f"{kind}_{REPORT_CUT}"
GRIDS = (64, 128, 256, 512)
NPARTS = (1024, 4096, 16384, 65536)
LONG = (300, 12)          # (slices, sample), per machine, longer than the line's slippage
DOUBLED = (600, 12)       # the same for the line twice
EXAMPLE = (96, 3)         # the examples' own window
INTERIOR = slice(80, 230)
DUMP_ELES = ("UND##1", "UND##2", "UND##4", "UND##8", "UND##12")
UND_END_RECORDS = {1: 0, 2: 4, 4: 12, 8: 28, 12: 44}   # record index of each undulator end
SHARP_WIDTH = 0.05        # the source filter's sigmoid width where the edge is meant to be sharp

DECK = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%comb_ds_save = -1
  slicing%window_length = {slen:.6e}
  slicing%n_wavelength = {sample}
{extra}/

&fel_beam_init
  beam_init%n_particle = {npart}
  beam_init%bunch_charge = {charge:.12e}
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -{zhalf:.6e}
  beam_init%grid(3)%x_max = {zhalf:.6e}
  beam_init%sig_pz = {sigpz:.12e}
  beam_init%a_norm_emit = {emit:.6e}
  beam_init%b_norm_emit = {emit:.6e}
  shot_noise = T
  beamlet_size = {beamlet}
/

&fel_wavefront_init
  wavefront_init%lambda0 = {lam0:.12e}
  wavefront_init%seed_power = 0
  wavefront_init%grid_n_pts = {ngrid}
  wavefront_init%grid_half_width = {half:.6e}
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
  slicing%n_wavelength = 3
/
&fel_beam_init
  beam_file = "{gen}-initial.beam.h5"
  beamlet_size = 8
/
&fel_wavefront_init
  field_file = "{gen}-initial.wf.h5"
  wavefront_init%lambda0 = 1e-10
/
"""


# ---------------------------------------------------------------------------
# Running and caching

def cell_size(ngrid):
    return 2 * HALF_WIDTH / (ngrid - 1)


def filter_xcut(ngrid, theta):
    """The source filter's xcut that puts its sigmoid's edge at the angle theta.

    Genesis normalizes the shifted grid index by ngrid, so xcut = 1 is the angle
    lambda/dx, twice the Nyquist angle (fel-physics.md sec-source-filter)."""
    return theta * cell_size(ngrid) / LAMBDA0


def deck_text(lat, root, ngrid, npart, beamlet, window, device, dumps=(), xcut=None, width=1.0,
              tolerance=None, half=None, more=""):
    nslice, sample = window
    slen = nslice * sample * LAMBDA0
    charge = CURRENT * slen / C_LIGHT
    extra = more
    if device != "off":
        extra += f'  global%device = "{device}"\n'
    if xcut is None:
        extra += '  global%source_filter = F\n'      # Silence now means on, and this page's
                                                      # unfiltered rows are the comparison.
    elif xcut == "default":
        extra += '  global%source_filter = T\n'
    else:
        extra += ('  global%source_filter = T\n'
                  f'  global%source_filter_xcut = {xcut:.9f}\n'
                  f'  global%source_filter_ycut = {xcut:.9f}\n'
                  f'  global%source_filter_width = {width:.9f}\n')
    if tolerance is not None:
        extra += f'  global%source_filter_tolerance = {tolerance:.9f}\n'
    if dumps:
        extra += "  global%dump_field_at = " + ", ".join(f'"{d}"' for d in dumps) + "\n"
    return DECK.format(lat=lat, root=root, ngrid=ngrid, npart=npart, beamlet=beamlet,
                       charge=charge, zhalf=slen / 2, slen=slen, sample=sample, extra=extra,
                       lam0=LAMBDA0, sigpz=SIG_PZ, emit=NORM_EMIT,
                       half=HALF_WIDTH if half is None else half)


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
    """The 1D Pierce parameter for the deck, Genesis's own form with the rms aw times the
    undulator's coupling, and the transverse size from the mean matched beta."""
    beta = 0.5 * (BETA_A + BETA_B)
    emit = NORM_EMIT / GAMMA0
    sigma = math.sqrt(emit * beta)
    ku = 2 * math.pi / LAMBDA_U
    rho3 = (CURRENT / I_ALFVEN) * (AW * FC) ** 2 / (8 * GAMMA0 ** 3 * sigma ** 2 * ku ** 2)
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
    incoherent undulator power of the beam in the central cone, their eq. (1), carrying
    the undulator's coupling squared, which is one for a helical device. The cone is the
    resonance's, so its own aw has no coupling in it. The share of that power inside the
    reported cut goes as the square of the angle and is capped at the whole cone, since
    outside the cone the emission is off resonance and the square is not its measure."""
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
    w_incoh = math.pi * E_CHARGE * CURRENT / (EPS0 * LAMBDA0) * k2 * FC ** 2 / (1 + k2)
    theta_cone_seg = math.sqrt(1 + k2) / (GAMMA0 * math.sqrt(SEG_LENGTH / LAMBDA_U))
    cut = THETA_CUTS[1]
    return {"N_lambda": n_lambda, "P_eff": p_eff, "N_c": n_c, "L_sat_ssy": l_sat,
            "eps_hat": eps_hat, "P_sat_ssy": p_sat, "W_incoh_cone": w_incoh,
            "theta_cone_one_segment": theta_cone_seg,
            "W_incoh_in_cut_one_segment": w_incoh * min(1.0, (cut / theta_cone_seg) ** 2),
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
    rel_in = max(abs(a[cut_key("in")] - b[cut_key("in")]) / b[cut_key("in")] for a, b in zip(dev, rows))
    res["cpu_cross_check"] = {"worst_total_rel": rel, "worst_inside_rel": rel_in, "dumps": rows}
    print(f"  device vs CPU at grid 64, 1024 particles: total {rel:.2e}, inside 3 urad {rel_in:.2e}")
    results["floor"] = res


def exp_filter(rn, args, lat, results):
    """
    The source filter against the split this page defines (fel-physics.md
    sec-source-filter). The filter suppresses the source at wide transverse angles, and
    what has never been measured is what that does to the power inside the mode, which is
    the part that amplifies. Grid 256, the cell size the page's scaling law is written
    at, and the two loads that bracket the examples. Two edge positions: the central cone
    of one segment, which is where the physical emission stops, and the 3 urad cut this
    page splits at, which is inside the cone and inside the coherent mode.
    """
    # The grid comes from the beam, as the page's own criterion says: cells of about
    # sigma/7. The second machine's beam is five times the first's, so its grid is coarser.
    ng = MACHINES[args.machine]["filter_grid"]
    print(f"== g. the source filter against the mode power, grid {ng} ==")
    th = results["theory"]
    # The edge the code derives, from the nominal beam. It takes the larger of four mode
    # diffraction angles and the angle at which the resonance red-shifts by rho, and the
    # run derives its own from the beam it loaded, which differs by the load's statistics.
    mode_edge, rho_edge = 4 * th["theta_mode"], math.sqrt(2 * th["rho"] * LAMBDA0 / LAMBDA_U)
    print(f"  derived edge {max(mode_edge, rho_edge):.3e} rad: four mode angles {mode_edge:.3e}, "
          f"rho angle {rho_edge:.3e}, ratio {max(mode_edge, rho_edge) / min(mode_edge, rho_edge):.2f}")
    edges = {"cone": th["theta_cone_one_segment"], "cut": THETA_CUTS[1]}   # the reported cut
    res = {"ngrid": ng, "dx": cell_size(ng), "edges": {}, "machine": args.machine}
    for name, theta in edges.items():
        res["edges"][name] = {"theta": theta, "xcut": filter_xcut(ng, theta)}

    # Genesis4's default width is 1 in the same normalized units as the edge, which is a
    # very soft roll-off: the sigmoid is 1/(1 + exp(-1/w)) = 0.73 on axis, so it attenuates
    # the coherent source as well as the wide angles. The sharp pair separates the filter's
    # angular selectivity from that softness.
    cases = [("off", None, 1.0), ("default", "default", 0.0)]
    if args.machine == "aramis":
        # The edge and width sweep is the first machine's. The second answers one
        # question, whether the derived default holds there, and its runs are on the CPU
        # at about 48 minutes each on four threads, so it takes the two cases that answer it.
        for k, v in res["edges"].items():
            cases.append((k, v["xcut"], 1.0))
            cases.append((k + "_sharp", v["xcut"], SHARP_WIDTH))
    res["sharp_width"] = SHARP_WIDTH

    # The second machine answers one question at one load. The first and the fourth carry
    # both, since the load is what the wide-angle share scales with.
    loads = (1024,) if args.machine in ("flash", "flash1") else (1024, 4096)
    for npart in loads:
        for name, xcut, width in cases:
            # The unfiltered case is experiment a's own deck, so the runner skips it when
            # that experiment has already run in this work directory.
            root = f"floor_g{ng}_n{npart}" if name == "off" else f"filt_{name}_n{npart}"
            t0 = time.time()
            rn.run(root, deck_text(lat, root, ng, npart, 8, LONG, args.device, DUMP_ELES,
                                   xcut=xcut, width=width))
            wall = time.time() - t0
            rows = analyze_dumps(rn.wd, root, rn.wd / f"{root}.farfield.json")
            z, P, b = element_end_power(rn.wd, root)
            zs, Ps = saturation_point(z, interior(P))
            res[f"n{npart}_{name}"] = {
                "npart": npart, "filter": name, "xcut": xcut, "width": width, "dumps": rows,
                "z": z.tolist(), "P": interior(P).tolist(), "b": interior(b).tolist(),
                "z_sat": zs, "P_sat": Ps, "wall_s": wall,
            }

    # The three numbers the page needs, at the last dump and at the one before it, which
    # is the last still inside the exponential regime on all three machines.
    late = sorted(UND_END_RECORDS)[-2]
    for npart in loads:
        base = res[f"n{npart}_off"]
        for name, _, _ in cases[1:]:
            r = res[f"n{npart}_{name}"]
            for tag, i in (("late", 3), ("exit", 4)):
                o0, o1 = base["dumps"][i][cut_key("out")], r["dumps"][i][cut_key("out")]
                i0, i1 = base["dumps"][i][cut_key("in")], r["dumps"][i][cut_key("in")]
                r[f"wide_factor_{tag}"] = o0 / o1 if o1 > 0 else float("inf")
                r[f"mode_rel_{tag}"] = (i1 - i0) / i0
            r["b_z_late"] = r["z"][UND_END_RECORDS[late]]
            r["b_rel_late"] = (r["b"][UND_END_RECORDS[late]] - base["b"][UND_END_RECORDS[late]]) \
                / base["b"][UND_END_RECORDS[late]]
            print(f"  {npart:5d} particles, edge at the {name}: wide-angle down "
                  f"{r['wide_factor_exit']:.1f}x at the exit, mode power "
                  f"{r['mode_rel_exit']:+.0%}, bunching at {r['b_z_late']:.0f} m "
                  f"{r['b_rel_late']:+.0%}, "
                  f"saturation {r['z_sat']:.1f} m against {base['z_sat']:.1f} m")
    results["filter"] = res


def fit_gain_length(z, power, z_sat):
    """
    Power gain length fitted through the exponential regime, in undulator metres.

    The fit runs from ten times the starting power up to two thirds of the way to
    saturation, which keeps the startup transient and the roll-over out of it. FILL
    converts the slope along the line into a gain length along undulator, which is the
    convention Ming Xie's estimate and this page's figures use.
    """
    p = np.asarray(power, dtype=float)
    live = np.isfinite(p) & (p > 0)
    if live.sum() < 4 or not z_sat:
        return None
    zz, pp = np.asarray(z, dtype=float)[live], p[live]

    # Every record a decade above the starting power and below saturation. A band set from
    # the saturated power instead collapses onto one metre on the two six-segment lines,
    # where the records are element ends and few, and it returned a gain length of 60 m
    # against a line that saturates at 21.6 m.

    band = (pp > 10 * pp[0]) & (zz < z_sat)
    if band.sum() < 4:
        return None
    slope = np.polyfit(zz[band], np.log(pp[band]), 1)[0]
    return float(FILL / slope) if slope > 0 else None


def accepted_power(wd, root, record):
    """Power inside the reported split angle at one record, summed over interior slices."""
    with h5py.File(wd / f"{root}.stats.h5") as h:
        g = h["field/total"]
        p_in = np.asarray(g["power_inside_angle"], dtype=float)
    row = p_in[record]
    return float(np.nansum(row[INTERIOR])) if np.isfinite(row).any() else None


def exp_tolerance(rn, args, lat, results):
    """
    Where the source filter's edge sits relative to the angle it protects
    (fel-physics.md sec-source-filter).

    The edge used to be the protected angle itself, which puts the sigmoid's
    half-amplitude point on it and takes three quarters of the source intensity there.
    global%source_filter_tolerance names the largest loss allowed inside that angle and
    moves the edge out to hold it. This sweeps it while the reported split stays on the
    protected angle, so every row is measured against one acceptance.

    One seed a row. That is an exploratory sweep: it says whether the tolerance reaches
    the gain at all and how far, and it cannot select a default, which needs paired
    ensembles over the seed spread this page measures at about 25 percent.
    """
    ng = MACHINES[args.machine]["filter_grid"]
    tols = (0.75, 0.2, 0.1, 0.01)
    print(f"== h. the filter's passband tolerance, grid {ng} ==")
    res = {"ngrid": ng, "machine": args.machine, "tolerances": list(tols), "rows": {}}
    first_end = UND_END_RECORDS[min(UND_END_RECORDS)]

    for tol in tols:
        root = f"tol{str(tol).replace('.', 'p')}_n1024"
        t0 = time.time()
        rn.run(root, deck_text(lat, root, ng, 1024, 8, LONG, args.device, DUMP_ELES,
                               xcut="default", tolerance=tol))
        wall = time.time() - t0
        rows = analyze_dumps(rn.wd, root, rn.wd / f"{root}.farfield.json")
        z, P, b = element_end_power(rn.wd, root)
        Pi = interior(P)
        zs, Ps = saturation_point(z, Pi)
        slen = LONG[1] * LAMBDA0
        row = {
            "tolerance": tol, "wall_s": wall,
            "accepted_first_segment_W": accepted_power(rn.wd, root, first_end),
            "z_sat": zs, "P_sat": Ps,
            "Lg_fit": fit_gain_length(z, Pi, zs),
            "pulse_energy_exit_J": float(np.nansum(P[-1][INTERIOR]) * slen / C_LIGHT),
        }
        for tag, i in (("late", 3), ("exit", 4)):
            row[f"wide_over_mode_{tag}"] = (rows[i][cut_key("out")] / rows[i][cut_key("in")]
                                            if rows[i][cut_key("in")] > 0 else float("inf"))
        res["rows"][str(tol)] = row
        def g(v, fmt=".4g"):
            return "n/a" if v is None else format(v, fmt)
        print(f"  tolerance {tol:<5}: accepted at the first segment "
              f"{g(row['accepted_first_segment_W'])} W, wide/mode "
              f"{g(row['wide_over_mode_late'], '.3g')} late and "
              f"{g(row['wide_over_mode_exit'], '.3g')} at the exit, gain length "
              f"{g(row['Lg_fit'])} m, saturation {g(row['z_sat'], '.1f')} m, pulse energy "
              f"{g(row['pulse_energy_exit_J'])} J, {wall:.0f} s")

    base = res["rows"]["0.75"]
    for tol in tols[1:]:
        r = res["rows"][str(tol)]
        for k in ("accepted_first_segment_W", "pulse_energy_exit_J", "Lg_fit"):
            if base[k] and r[k]:
                r[k + "_rel"] = (r[k] - base[k]) / base[k]
        r["wide_over_mode_exit_factor"] = (r["wide_over_mode_exit"] / base["wide_over_mode_exit"]
                                           if base["wide_over_mode_exit"] > 0 else float("inf"))
        print(f"  against 0.75, tolerance {tol}: accepted {r.get('accepted_first_segment_W_rel', 0):+.1%}, "
              f"pulse energy {r.get('pulse_energy_exit_J_rel', 0):+.1%}, gain length "
              f"{r.get('Lg_fit_rel', 0):+.2%}, wide/mode x{r['wide_over_mode_exit_factor']:.2f}, "
              f"saturation {r['z_sat']:.1f} m against {base['z_sat']:.1f} m")
    results["tolerance"] = res


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
    # The growth rate dips to zero in every drift, and a run whose startup is far below
    # its saturated power has a peak rate high enough for such a dip to pass the test in
    # the exponential regime. Saturation is where the growth stops near the top, so the
    # power there has to be within a decade of the run's maximum.
    pmax = float(np.max(P))
    for i in range(i0, len(z)):
        if rate[i] < 0.1 * rate[i0] and P[i] > 0.1 * pmax:
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



# ---------------------------------------------------------------------------
# i. The cell-size exponent, derived and measured to its cause
#
# The source's exponent is 2 (deposit-demo.ipynb) and the reported floor's falls from
# 2 toward 1 as the cells shrink (doc/startup-noise.md's floor table). The derivation
# in the design record (docs/analyses/cell-size-exponent.md) says why: the field
# record slips one slice every L_ref = sample * lambda_u of undulator, so a field slice
# sees the same static source for N_ref = L_ref / dz steps and an independent one
# after, a mode at k_perp dephases by phi = k_perp^2 dz / (2 k_s) between deposits, and
# a block of N deposits adds as sin^2(N phi/2)/sin^2(phi/2), N^2 while N phi is small
# and N once it is not. The knee is theta_t = sqrt(lambda / (sample lambda_u)). The
# same sum is evaluated here beside every measurement, with nothing fitted, so the
# JSON carries the prediction and the number in one place.

CAUSE_GRIDS = (64, 128, 256, 512)
CAUSE_DS = (0.0225, 0.045, 0.09, 0.18)      # m, N_ref of 8, 4, 2 and 1 at sample 12
CAUSE_SAMPLES = (3, 6, 12, 24, 48)          # wavelengths a slice, theta_t of 47 to 12 urad
CAUSE_WIDTHS = (2e-4, 4e-4, 8e-4)           # m, at cells near 3.15 and 1.57 um
CAUSE_STAIR_LENGTH = 0.585                  # m, thirteen steps, three slice rotations
CAUSE_PLANE_COMB = 0.27                     # m, a frame every six steps


def cause_blocks(dz, sample, z):
    """The deposit blocks a field slice sees, in steps, replaying the tracker's own
    rotation rule: the record turns when the accumulated slip passes 0.8 of a slice."""
    nstep = int(round(SEG_LENGTH / dz))
    dz_eff = SEG_LENGTH / nstep
    slip = dz_eff / LAMBDA_U
    n = int(round(z / dz_eff))
    acc, cur, out = 0.0, 0, []
    for _ in range(n):
        cur += 1
        acc += slip
        if abs(acc) > 0.8 * sample:
            acc -= sample
            out.append(cur)
            cur = 0
    return out + ([cur] if cur else []), dz_eff


def cause_model(ngrid, half, dz, sample, z, cut, shape="cic", sigma_xp=None, nxi=33):
    """The wide-angle power outside the cut at plane z, up to a constant common to every
    grid: the sum over the grid's own modes of the shape's power transform times the
    block factor, summed over the blocks. Returns (power, cell size)."""
    if sigma_xp is None:
        sigma_xp = math.sqrt(NORM_EMIT / GAMMA0 / (0.5 * (BETA_A + BETA_B)))
    ks = 2 * math.pi / LAMBDA0
    blocks, dz_eff = cause_blocks(dz, sample, z)
    dx = 2 * half / (ngrid - 1)
    k1 = 2 * np.pi * np.fft.fftfreq(ngrid, d=dx)
    kx, ky = np.meshgrid(k1, k1)
    k2 = kx ** 2 + ky ** 2
    theta = np.sqrt(k2) / ks
    shp = (np.sinc(kx * dx / (2 * np.pi)) * np.sinc(ky * dx / (2 * np.pi))) ** (2 if shape == "cic" else 1)
    phi = k2 * dz_eff / (2 * ks)
    xi = np.linspace(-3, 3, nxi)
    wxi = np.exp(-xi ** 2 / 2)
    wxi /= wxi.sum()
    g = np.zeros_like(k2)
    for nb in blocks:
        m = np.arange(nb)
        for x, w in zip(xi, wxi):
            ph = (phi - np.sqrt(k2) * sigma_xp * x * dz_eff).ravel()
            g += w * np.abs(np.exp(-1j * np.outer(m, ph)).sum(axis=0)).reshape(k2.shape) ** 2
    return float((shp ** 2 * g * (theta > cut)).sum()), dx


def cause_exponents(powers):
    """The exponent in 1/dx over each halving of the cell, from the powers at successive
    grids, coarse to fine."""
    p = np.asarray(powers, dtype=float)
    return (np.log2(p[1:] / p[:-1])).tolist()


def cause_frames(wd, root, cut):
    """Far-field splits of a run's comb frames in record order, cached as JSON beside the
    run so the frames can go and a second pass reproduces the numbers without them."""
    cache = wd / f"{root}.farfield.json"
    if cache.exists():
        return json.loads(cache.read_text())
    files = sorted(wd.glob(f"{root}-[0-9]*.wf.h5"),
                   key=lambda p: int(p.name.split("-")[-1].split(".")[0]))
    if not files:
        raise RuntimeError(f"{root}: no comb frames and no cache")
    rows = [far_field_split(p, (cut,)) for p in files]
    cache.write_text(json.dumps(rows))
    for p in files:
        p.unlink()
    for p in wd.glob(f"{root}-[0-9]*.beam.h5"):
        p.unlink()
    return rows


def cause_fit(results):
    """The scaling law's cell-size factor, re-fitted on the data it was built from: the
    ratio of the power outside the cut to the power inside it at the saturation dump of
    the floor runs, over the law's own range of 128 to 2048 beamlets. Two factors are
    fitted, the (1.57 um/dx)^2 the law carries and the floor's own factor as this block
    measures it, and the residual range of each is what says which to keep. The run at
    0.78 um cells and 128 beamlets is set aside from the fit and reported beside it: the
    artifact drains that beam and the ratio there is a property of the drain."""
    fl = results["floor"]
    c = results["cause"]
    cut = c["cut"]
    ko, ki = f"out_{cut:.0e}", f"in_{cut:.0e}"
    base = {r["ngrid"]: r for r in c["base"]["rows"]}
    ref = base[256]
    factors = {"law": {ng: (ref["dx"] / base[ng]["dx"]) ** 2 for ng in base},
               "floor": {ng: base[ng]["out"] / ref["out"] for ng in base}}
    idump = 3    # the eighth undulator's end, the saturation plane the law is stated at
    pts = [(ng, n // 8, fl[f"floor_g{ng}_n{n}"]["dumps"][idump])
           for ng in CAUSE_GRIDS for n in NPARTS]
    z_sat = pts[0][2]["z"]
    drained = [(512, 128)]
    out = {"z": z_sat, "beamlets_fitted": [128, 512, 2048], "set_aside": drained,
           "factor": {k: {str(ng): v[ng] for ng in CAUSE_GRIDS} for k, v in factors.items()},
           "ratio": {str(ng): {str(nb): p[ko] / p[ki] for g, nb, p in pts if g == ng} for ng in CAUSE_GRIDS},
           "saturation_exponent": {}, "fit": {}, "criterion": {}}
    for nb in (128, 512, 2048, 8192):
        o = [p[ko] for g, b, p in pts if b == nb]
        out["saturation_exponent"][str(nb)] = cause_exponents(o)
    sel = [q for q in pts if q[1] <= 2048 and (q[0], q[1]) not in drained]
    for name, F in factors.items():
        a = math.exp(np.mean([math.log(p[ko] / p[ki] * nb / F[ng]) for ng, nb, p in sel]))
        res = [p[ko] / p[ki] / (a * F[ng] / nb) for ng, nb, p in sel]
        out["fit"][name] = {"prefactor": a, "residual_low": min(res), "residual_high": max(res)}
        # The criterion N_b > 1000 F(dx), and the ratio it lands on at each cell, read from the
        # measured ratios by interpolation in log N_b.
        crit = {}
        for ng in CAUSE_GRIDS:
            nbs = np.array([128, 512, 2048, 8192], dtype=float)
            rat = np.array([out["ratio"][str(ng)][str(int(b))] for b in nbs])
            need = 1000 * F[ng]
            crit[str(ng)] = {"beamlets": need,
                             "ratio_there": float(np.exp(np.interp(math.log(need), np.log(nbs), np.log(rat))))}
        out["criterion"][name] = crit
    c["fit"] = out
    print(f"  law re-fitted at z = {z_sat:.1f} m over {out['beamlets_fitted']} beamlets:")
    for name in factors:
        f = out["fit"][name]
        print(f"    {name:6s} factor: prefactor {f['prefactor']:.0f}, residual {f['residual_low']:.2f} to {f['residual_high']:.2f}")


def exp_cause(rn, args, results):
    print("== i. the cell-size exponent, one segment, 1024 particles ==")
    cut = THETA_CUTS[1]
    key_out, key_in = f"out_{cut:.0e}", f"in_{cut:.0e}"
    lat1 = "aramis_1seg.bmad"
    src = pathlib.Path(args.latdir) / "bmad" / lat1
    (rn.wd / lat1).write_bytes(src.read_bytes())
    npart, beamlet = 1024, 8
    nslice = LONG[0]
    res = {"cut": cut, "npart": npart, "nslice": nslice, "seg_length": SEG_LENGTH,
           "lambda_u": LAMBDA_U, "ds_step_lattice": 0.045}

    def one(root, lat, ngrid, sample, half, dz, more="", dumps=("UND",), device=None, threads=None):
        text = deck_text(lat, root, ngrid, npart, beamlet, (nslice, sample), device or args.device,
                         dumps=dumps, half=half, more=more)
        rn.run(root, text, env_threads=threads)
        row = None
        if dumps:
            cache = rn.wd / f"{root}.farfield.json"
            if cache.exists():
                row = json.loads(cache.read_text())
            else:
                row = far_field_split(dumps_of(rn.wd, root, nele=1)[0], THETA_CUTS)
                cache.write_text(json.dumps(row))
        pred, dx = cause_model(ngrid, half, dz, sample, SEG_LENGTH, cut)
        out = {"ngrid": ngrid, "dx": dx, "sample": sample, "half": half, "dz": dz, "model": pred}
        if row is not None:
            out.update({"out": row[key_out], "in": row[key_in], "theta_nyquist": row["theta_nyquist"]})
        return out

    def sweep(rows):
        """Measured and predicted exponents over the halvings, coarse to fine."""
        return {"rows": rows, "exponent": cause_exponents([r["out"] for r in rows]),
                "exponent_model": cause_exponents([r["model"] for r in rows])}

    def say(label, sw):
        print(f"  {label:14s} exponents", [f"{e:.2f}" for e in sw["exponent"]],
              "derived", [f"{e:.2f}" for e in sw["exponent_model"]])

    # The base: the lattice's own step and the sweep's own spacing and window.
    base = [one(f"cause_base_g{ng}", lat1, ng, LONG[1], HALF_WIDTH, 0.045) for ng in CAUSE_GRIDS]
    res["base"] = sweep(base)
    say("base", res["base"])

    # The step, at the spacing and window of the base. A wrapper sets the element's step.
    res["step"] = {}
    for dz in CAUSE_DS:
        lat = lat1
        tag = int(round(dz * 1e4))
        if dz != 0.045:
            lat = f"cause_ds{tag}.bmad"
            (rn.wd / lat).write_text(f"call, file = {lat1}\nUND[ds_step] = {dz}\n")
        rows = [one(f"cause_ds{tag}_g{ng}", lat, ng, LONG[1], HALF_WIDTH, dz) for ng in CAUSE_GRIDS]
        res["step"][str(dz)] = sweep(rows)
        say(f"dz {dz:.4f}", res["step"][str(dz)])

    # The slice spacing, at the step and window of the base.
    res["spacing"] = {}
    for smp in CAUSE_SAMPLES:
        rows = [one(f"cause_s{smp}_g{ng}", lat1, ng, smp, HALF_WIDTH, 0.045) for ng in CAUSE_GRIDS]
        res["spacing"][str(smp)] = sweep(rows)
        say(f"sample {smp}", res["spacing"][str(smp)])

    # The half width, at cells near 3.15 and 1.57 um. The device takes powers of two, so
    # the cell moves by 0.4 percent between widths; both grids of a pair move together and
    # the exponent over the pair is what is compared. The 1024-point grid is the CPU's.
    res["width"] = {}
    for half in CAUSE_WIDTHS:
        grids = tuple(1 << int(round(math.log2(2 * half / d + 1))) for d in (3.15e-6, 1.569e-6))
        cpu = grids[1] > 512
        rows = [one(f"cause_w{int(half * 1e6)}_g{ng}", lat1, ng, LONG[1], half, 0.045,
                    device="off" if cpu else None, threads=args.cpu_threads if cpu else None)
                for ng in grids]
        res["width"][str(half)] = sweep(rows)
        say(f"half {half * 1e6:.0f} um", res["width"][str(half)])

    # The plane: frames along the segment, the wide-angle power against z.
    res["plane"] = {}
    for ng in (128, 256):
        root = f"cause_plane_g{ng}"
        one(root, lat1, ng, LONG[1], HALF_WIDTH, 0.045, dumps=(),
            more=f"  global%dump_at_comb = T\n  global%comb_ds_save = {CAUSE_PLANE_COMB}\n")
        rows = cause_frames(rn.wd, root, cut)
        zs = [r["z"] for r in rows]
        pred = [cause_model(ng, HALF_WIDTH, 0.045, LONG[1], z, cut)[0] if z > 0 else 0.0 for z in zs]
        res["plane"][str(ng)] = {"z": zs, "out": [r[key_out] for r in rows],
                                 "in": [r[key_in] for r in rows], "model": pred}
    # The staircase: every step over three slice rotations of a short element.
    lat = "cause_stair.bmad"
    (rn.wd / lat).write_text(f"call, file = {lat1}\nUND[l] = {CAUSE_STAIR_LENGTH}\n")
    root = "cause_stair_g256"
    one(root, lat, 256, LONG[1], HALF_WIDTH, 0.045, dumps=(),
        more="  global%dump_at_comb = T\n  global%comb_ds_save = 0\n")
    rows = cause_frames(rn.wd, root, cut)
    res["stair"] = {"z": [r["z"] for r in rows], "out": [r[key_out] for r in rows]}

    # The CPU cross-check at one grid of the base.
    cpu = one("cause_base_g128_cpu", lat1, 128, LONG[1], HALF_WIDTH, 0.045, device="off",
              threads=args.cpu_threads)
    rel = abs(cpu["out"] - base[1]["out"]) / base[1]["out"]
    res["cpu_cross_check"] = {"rel_out": rel}
    print(f"  device vs CPU at grid 128, wide-angle power: {rel:.2e}")
    results["cause"] = res
    if "floor" in results:
        cause_fit(results)


def cause_figure(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
    c = results["cause"]
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))

    rows = c["base"]["rows"]
    dx = np.array([r["dx"] for r in rows]) * 1e6
    meas = np.array([r["out"] for r in rows])
    mod = np.array([r["model"] for r in rows]) * (meas[2] / rows[2]["model"])
    a = ax[0, 0]
    a.loglog(dx, meas, "o", label="measured")
    a.loglog(dx, mod, "-", label="derived, scaled at 1.57 um")
    a.loglog(dx, meas[2] * (dx[2] / dx) ** 2, "--", color="gray", label="$1/dx^2$")
    a.set_xlabel("cell size (um)")
    a.set_ylabel("power per slice outside 3 urad (W)")
    a.set_title("the floor at z = 3.99 m")
    a.legend(fontsize=8)

    def expo_panel(a, key, xs, xlabel, title, which):
        a.semilogx(xs, [c[key][str(x)]["exponent"][which] for x in xs], "o", label="measured")
        a.semilogx(xs, [c[key][str(x)]["exponent_model"][which] for x in xs], "-", label="derived")
        a.axhline(2, color="gray", ls="--", lw=0.8)
        a.set_xticks(list(xs))
        a.set_xticklabels([f"{x:g}" for x in xs])
        a.minorticks_off()
        a.set_xlabel(xlabel)
        a.set_ylabel("exponent over one halving")
        a.set_title(title)
        a.set_ylim(0, 2.6)
        a.legend(fontsize=8)

    expo_panel(ax[0, 1], "spacing", CAUSE_SAMPLES, "slice spacing (wavelengths)",
               "the slice spacing moves the knee, 3.15 to 1.57 um", 1)
    expo_panel(ax[0, 2], "step", CAUSE_DS, "ds_step (m)", "the step at fixed spacing, 1.57 to 0.78 um", 2)
    a = ax[1, 0]
    ws = [w * 1e6 for w in CAUSE_WIDTHS]
    a.semilogx(ws, [c["width"][str(w)]["exponent"][0] for w in CAUSE_WIDTHS], "o", label="measured")
    a.semilogx(ws, [c["width"][str(w)]["exponent_model"][0] for w in CAUSE_WIDTHS], "-", label="derived")
    a.axhline(2, color="gray", ls="--", lw=0.8)
    a.set_xticks(ws)
    a.set_xticklabels([f"{w:g}" for w in ws])
    a.minorticks_off()
    a.set_xlabel("half width (um)")
    a.set_ylabel("exponent, 3.15 to 1.57 um")
    a.set_title("the window does not enter")
    a.set_ylim(0, 2.6)
    a.legend(fontsize=8)

    a = ax[1, 1]
    for ng, mk in ((128, "o"), (256, "s")):
        p = c["plane"][str(ng)]
        z = np.array(p["z"])
        o = np.array(p["out"])
        m = np.array(p["model"])
        keep = z > 0
        a.plot(z[keep], o[keep] / o[keep][-1], mk, label=f"grid {ng}, measured")
        a.plot(z[keep], m[keep] / m[keep][-1], "-", lw=0.8, label=f"grid {ng}, derived")
    a.plot([0, c["seg_length"]], [0, 1], "--", color="gray", label="linear in z")
    a.set_xlabel("z (m)")
    a.set_ylabel("power outside 3 urad, relative to the end")
    a.set_title("the plane: blocks add in power")
    a.legend(fontsize=7)

    a = ax[1, 2]
    st = c["stair"]
    z = np.array(st["z"])
    o = np.array(st["out"])
    keep = z > 0
    dz_eff = CAUSE_STAIR_LENGTH / round(CAUSE_STAIR_LENGTH / c["ds_step_lattice"])
    a.plot(z[keep] / dz_eff, o[keep], "o-")
    for k in (4, 8, 12):
        a.axvline(k, color="gray", ls=":", lw=0.8)
    a.set_xlabel("step")
    a.set_ylabel("power per slice outside 3 urad (W)")
    a.set_title("the first three slice rotations, grid 256")
    fig.tight_layout()
    fig.savefig(out / "cell-size-cause.png", dpi=130)
    plt.close(fig)


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

def filter_figure(results, out):
    """The source filter against the split, one panel per load."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
    fi = results["filter"]
    cases = [("off", "no filter", "k", "-"),
             ("cone", "edge at the cone, width 1", "tab:orange", "-"),
             ("cut", "edge at 3 urad, width 1", "tab:red", "-"),
             ("cone_sharp", f"edge at the cone, width {fi['sharp_width']}", "tab:green", "--"),
             ("cut_sharp", f"edge at 3 urad, width {fi['sharp_width']}", "tab:blue", "--")]
    loads = [n for n in (1024, 4096) if f"n{n}_off" in fi]
    fig, ax = plt.subplots(1, len(loads), figsize=(5 * len(loads), 4), sharey=True, squeeze=False)
    ax = ax[0]
    for j, npart in enumerate(loads):
        for name, label, c, ls in cases:
            if f"n{npart}_{name}" not in fi:
                continue          # the second machine runs the default alone
            r = fi[f"n{npart}_{name}"]
            z = [d["z"] for d in r["dumps"]]
            ax[j].semilogy(z, [d[cut_key("in")] for d in r["dumps"]], ls, color=c, marker="o",
                           ms=3, label=label if j == 0 else None)
            ax[j].semilogy(z, [d[cut_key("out")] for d in r["dumps"]], ls, color=c, marker="x",
                           ms=4, alpha=0.45)
        ax[j].set_title(f"{npart} macroparticles per slice")
        ax[j].set_xlabel("z [m]")
    ax[0].set_ylabel("power per slice [W]")
    ax[0].legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Inside {float(REPORT_CUT) * 1e6:.3g} urad (circles) and outside it (crosses), "
                 "with and without the source filter", fontsize=10)
    fig.tight_layout()
    # The second machine's figure carries its name, so the two never overwrite each other.
    suffix = "" if fi.get("machine", "aramis") == "aramis" else f"-{fi['machine']}"
    fig.savefig(pathlib.Path(out) / f"source-filter{suffix}.png", dpi=130)
    plt.close(fig)


def figures(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
    cut = REPORT_CUT
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
    ap.add_argument("--only", default="a,b,c,d,e,f,g,h,i")
    ap.add_argument("--machine", default="aramis", choices=sorted(MACHINES))
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rn = Runner(args.exe, args.workdir, "4" if args.device != "off" else args.cpu_threads)
    mach = select_machine(args.machine)
    lat = mach["lat"]
    ex = pathlib.Path(args.examples)
    for src in (ex / lat, ex / args.machine / lat, pathlib.Path(args.latdir) / "bmad" / lat):
        if src.exists():
            break
    else:
        raise SystemExit(f"{lat}: not found under {ex} or {args.latdir}/bmad")
    (rn.wd / lat).write_bytes(src.read_bytes())
    lat2 = "aramis_twice.bmad"
    if mach["lat2"]:
        (rn.wd / lat2).write_text("call, file = aramis.bmad\nARAMIS2: line = (ARAMIS, ARAMIS)\nuse, ARAMIS2\n")

    # The second machine writes beside the first, never over it.
    suffix = "" if args.machine == "aramis" else f"-{args.machine}"
    results_file = out / f"startup-noise{suffix}.json"
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
    if ("c" in want or "d" in want) and mach["lat2"]:
        exp_line(rn, args, lat, lat2, results)
    if "f" in want:
        exp_genesis(rn, args, results)
    if "g" in want:
        exp_filter(rn, args, lat, results)
    if "h" in want:
        exp_tolerance(rn, args, lat, results)
    if "i" in want and args.machine == "aramis":
        exp_cause(rn, args, results)
    results_file.write_text(json.dumps(rounded(results), indent=1))
    if "cause" in results:
        cause_figure(results, out)
    if all(k in results for k in ("floor", "beamlets", "line")):
        figures(results, out)
    if "filter" in results:
        filter_figure(results, out)
    print(f"wrote {results_file}")


if __name__ == "__main__":
    main()
