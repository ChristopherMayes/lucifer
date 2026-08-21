#!/usr/bin/env python3
"""
Phasing checks (manual sec:phasing): how the beam-field phase behaves between
undulator segments, held by symmetries, closed forms and Genesis itself.

  1. RE-ANCHOR BASELINE (relative mode, the default): scanning an inter-segment gap
     by fractions of 2 gamma^2 lambda leaves the bunching phase entering the next
     segment FLAT -- the transcription of Genesis's drift autophasing (measured flat
     on Genesis itself the same way; the residual is a slow adiabatic trend).
  2. THE Z_OFFSET KNOB: displacing the downstream wiggler by delta (standard Bmad
     misalignment, anchor at the nominal position) shifts that phase by EXACTLY
     -2 pi delta / (2 gamma^2 lambda) -- the analytic drift slip rate, no fit.
  3. CROSS-MODE IDENTITY: in absolute mode (bmad_com[absolute_time_tracking] = T in
     the lattice, honored through Bmad's own resolver) the same phase arrives via the
     real gap length instead: the absolute-mode gap scan must reproduce the relative-
     mode knob scan point by point.
  4. PHASE-SHIFTER PARITY vs GENESIS4: Genesis scans PHASESHIFTER phi (which needs
     FINITE LENGTH to register -- a zero-length one silently does nothing); we scan
     z_offset with delta = phi * 2 gamma^2 lambda / (2 pi). Same physical scan; the
     bunching-phase curves must agree point by point. This anchors the SIGN
     convention ("delay goes backwards") against Genesis's own element.
  5. CHICANE: a four-bend closed bump between segments (Bmad seam tracks the beam;
     the radiation drifts the CHORD from ele%floor). Relative mode drops the
     geometric fraction (Genesis's "chicane is always autophasing"); absolute mode
     ramps with the bend angle at the slope an INDEPENDENT geometric computation of
     d(arc - chord)/d(angle) predicts. The unaveraged ledger closes across a
     chicane sandwich.
  6. REFUSALS by name: a non-closed-bump break; a bend under the genesis-model
     interludes; a z_offset exceeding its upstream break; a z_offset on the first
     element (no break to displace into).

Run by the benchmark harness; exits nonzero on failure. Needs genesis4 for section 4.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

FAILED = False

TWO_G2L = 2 * 11357.82**2 * 1e-10       # 2 gamma^2 lambda: one full turn of gap phase.
LAM = 1e-10

TOL_FLAT = 5e-3          # measured 1.5e-3 span: the re-anchor baseline's adiabatic trend.
TOL_KNOB = 5e-3          # measured 1.5e-3: knob ramp vs the exact analytic slope.
TOL_XMODE = 5e-3         # measured 1.2e-3: absolute gap scan vs relative knob scan.
TOL_PARITY = 1e-4        # measured 6.0e-6: our knob curve vs Genesis's PHASESHIFTER curve.
TOL_CHIC_SLOPE = 2e-3    # measured 6.8e-4: absolute chicane ramp vs the geometric slope
                         #   (583 wavelengths of delay at the base angle).
TOL_LEDGER = 1e-4        # measured 4.0e-6: unaveraged ledger closure across the chicane.

LAT2SEG = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
{absline}
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
UND2: UND, z_offset = {zoff}
D: pipe, l = {gap}
SEG: line = (UND, D, UND2)
use, SEG
{modeline}
"""

LATCHIC = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
{absline}
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
ANG = {ang}
B1: sbend, l = 0.05, g =  ANG / 0.05
B2: sbend, l = 0.05, g = -ANG / 0.05
B3: sbend, l = 0.05, g = -ANG / 0.05
B4: sbend, l = 0.05, g =  ANG / 0.05
DD: pipe, l = 0.025
SEG: line = (UND, DD, {mid}, DD, UND)
use, SEG
{modeline}
"""
CHIC_MID = "B1, DD, B2, DD, B3, DD, B4"
OPEN_MID = "B1, DD, B2, DD, B3"          # Three bends: NOT a closed bump.

TRANSCRIBED = "fel_transcribed = -1\nwiggler::*[FEL_TRACKING] = fel_transcribed"
UNAVG = "fel_unaveraged = 1\nwiggler::*[FEL_TRACKING] = fel_unaveraged"

NML = """&fel_track_params
  lat_file = "{tag}.bmad"
  out_root = "{tag}"
  lambda0 = 1e-10
  interlude_model = '{imodel}'
  beam_init%n_particle = 8192
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 1e7
  seed_waist_size = 30e-6
  grid_n_pts = 151
  grid_half_width = 2e-4
  ran_seed = 777
/
"""

GEN_DECK = """&setup
rootname={root}
lattice={lat}
beamline=SEG
lambda0=1e-10
gamma0=11357.82
delz=0.015
shotnoise=0
nbins = 8
field_global_stat = true
&end

&field
power=1e7
dgrid=2.000000e-04
ngrid=151
waist_size=30e-6
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
betax=15
betay=15
&end

&track
fft_fieldsolver = true
&end
"""

GEN_LAT_PS = """UND: UNDULATOR = {{ lambdau=0.015000, nwig=30, aw=0.84853, helical=True}};
D1: DRIFT = {{ l = 0.1425 }};
PS: PHASESHIFTER = {{ l = 0.015, phi = {phi} }};
D2: DRIFT = {{ l = 0.1425 }};
SEG: LINE={{UND,D1,PS,D2,UND}};
"""


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def wrap(dp):
    return (np.asarray(dp) + np.pi) % (2 * np.pi) - np.pi


def run(exe, wd, tag, lat_text, imodel="bmad", extra=""):
    (wd / f"{tag}.bmad").write_text(lat_text)
    nml = NML.format(tag=tag, imodel=imodel)
    if extra:
        nml = nml.replace("/\n", extra + "/\n")
    (wd / f"{tag}.nml").write_text(nml)
    r = subprocess.run([str(exe), f"{tag}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {tag} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)


def refuse(exe, wd, tag, lat_text, fragment, imodel="bmad"):
    (wd / f"{tag}.bmad").write_text(lat_text)
    (wd / f"{tag}.nml").write_text(NML.format(tag=tag, imodel=imodel))
    r = subprocess.run([str(exe), f"{tag}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    return r.returncode != 0 and fragment in (r.stdout + r.stderr)


def bphase(wd, tag, z_probe):
    d = np.loadtxt(wd / f"{tag}.diag.txt").reshape(-1, 1, 12)
    z = d[:, 0, 0]
    return float(d[np.searchsorted(z, z_probe), 0, 5])


def chicane_delay(ang):
    """arc - chord of the four-bend closed bump, EXACT 2D geometry, independently of
    the walk's floor arithmetic: bends of arc length lb and bend angle +a,-a,-a,+a
    with drifts ld between (and the outer DD pieces adding straight length only)."""
    lb, ld = 0.05, 0.025
    # Trace the reference: (x, z, angle); bend of angle a: chord lb*sin(a/2)/(a/2) at
    # heading (angle_in + a/2); drift: ld at current heading.
    x = z = th = 0.0
    arc = 0.0
    for piece, a in (("b", ang), ("d", 0), ("b", -ang), ("d", 0), ("b", -ang), ("d", 0), ("b", ang)):
        if piece == "b":
            c = lb * np.sin(a / 2) / (a / 2)
            x += c * np.sin(th + a / 2)
            z += c * np.cos(th + a / 2)
            th += a
            arc += lb
        else:
            x += ld * np.sin(th)
            z += ld * np.cos(th)
            arc += ld
    # The outer DD pieces are straight and cancel in arc - chord.
    chord = np.hypot(x, z)
    return arc - chord


def main():
    global FAILED
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--genesis", required=True)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(args.exe).resolve()

    NSCAN = 6
    ABS = "bmad_com[absolute_time_tracking] = T"

    # ------------------------------------------------------------------
    print("== re-anchor baseline + z_offset knob + cross-mode identity ==")
    base_bp, knob_bp, xmode_bp = [], [], []
    for k in range(NSCAN):
        f = k / 8
        run(exe, wd, f"pb{k}", LAT2SEG.format(absline="", zoff="0",
            gap=repr(0.30 + f * TWO_G2L), modeline=TRANSCRIBED), imodel="genesis")
        base_bp.append(bphase(wd, f"pb{k}", 0.76))
        run(exe, wd, f"pk{k}", LAT2SEG.format(absline="", zoff=repr(f * TWO_G2L),
            gap="0.30", modeline=TRANSCRIBED), imodel="genesis")
        knob_bp.append(bphase(wd, f"pk{k}", 0.76))
        run(exe, wd, f"px{k}", LAT2SEG.format(absline=ABS, zoff="0",
            gap=repr(0.30 + f * TWO_G2L), modeline=TRANSCRIBED), imodel="genesis")
        xmode_bp.append(bphase(wd, f"px{k}", 0.76))

    d_base = wrap(np.array(base_bp) - base_bp[0])
    check("re-anchor baseline: gap scan flat (rad, span)", float(np.ptp(d_base)), TOL_FLAT)

    d_knob = wrap(np.array(knob_bp) - knob_bp[0])
    expect = wrap(-2 * np.pi * np.arange(NSCAN) / 8)
    resid = float(np.max(np.abs(wrap(d_knob - expect - d_base))))   # trend-corrected
    check("z_offset knob vs -2pi delta/(2 gamma^2 lambda)", resid, TOL_KNOB)

    d_x = wrap(np.array(xmode_bp) - xmode_bp[0])
    check("cross-mode identity: absolute gap scan == knob scan",
          float(np.max(np.abs(wrap(d_x - d_knob)))), TOL_XMODE)

    # ------------------------------------------------------------------
    print("== phase-shifter parity vs Genesis4 ==")
    gen_bp = []
    for k in range(NSCAN):
        (wd / f"ps{k}.lat").write_text(GEN_LAT_PS.format(phi=repr(2 * np.pi * k / 8)))
        (wd / f"ps{k}.in").write_text(GEN_DECK.format(root=f"PS{k}", lat=f"ps{k}.lat"))
        r = subprocess.run([args.genesis, f"ps{k}.in"], cwd=wd, capture_output=True, text=True,
                           env={"PATH": "/usr/bin:/bin", "FI_PROVIDER": "tcp"})
        if r.returncode != 0:
            print(f"FAIL: genesis ps{k}:\n{r.stdout[-800:]}")
            sys.exit(1)
        with h5py.File(wd / f"PS{k}.out.h5") as h5:
            z = h5["Lattice/zplot"][:].ravel()
            aw = h5["Lattice/aw"][:].ravel()
            bp = h5["Beam/bunchingphase"][:].ravel()
        n = min(len(z), len(aw), len(bp))
        i2 = np.where((aw[:n] > 0) & (z[:n] > 0.5))[0][0] + 1
        gen_bp.append(float(bp[i2]))
    d_gen = wrap(np.array(gen_bp) - gen_bp[0])
    check("our z_offset curve == Genesis PHASESHIFTER curve",
          float(np.max(np.abs(wrap(d_knob - d_gen)))), TOL_PARITY,
          note="[delta = phi 2 gamma^2 lambda / 2pi]")

    # ------------------------------------------------------------------
    print("== chicane: relative flat, absolute at the geometric slope, ledger ==")
    a0 = 1e-3
    da = LAM / 8 / (2 * a0 * 0.058)     # ~lam/8 of delay per step at the true L_eff scale
    rel_bp, abs_bp, delays = [], [], []
    for k in range(NSCAN):
        a = a0 + k * da
        run(exe, wd, f"cr{k}", LATCHIC.format(absline="", ang=repr(a), mid=CHIC_MID,
                                              modeline=TRANSCRIBED))
        rel_bp.append(bphase(wd, f"cr{k}", 0.83))
        run(exe, wd, f"ca{k}", LATCHIC.format(absline=ABS, ang=repr(a), mid=CHIC_MID,
                                              modeline=TRANSCRIBED))
        abs_bp.append(bphase(wd, f"ca{k}", 0.83))
        delays.append(chicane_delay(a))

    check("chicane, relative mode: geometric fraction dropped (flat, rad span)",
          float(np.ptp(wrap(np.array(rel_bp) - rel_bp[0]))), TOL_FLAT)

    d_abs = wrap(np.array(abs_bp) - abs_bp[0])
    d_pred = wrap(-2 * np.pi * (np.array(delays) - delays[0]) / LAM)
    check("chicane, absolute mode vs independent geometry (rad)",
          float(np.max(np.abs(wrap(d_abs - d_pred)))), TOL_CHIC_SLOPE,
          note=f"[delay(a0) = {delays[0]:.3e} m = {delays[0]/LAM:.0f} wavelengths]")

    # The unaveraged ledger closes across a chicane sandwich (energy bookkeeping
    # survives the seam detour; the ledger's columns cover the unaveraged segments).
    run(exe, wd, "cled", LATCHIC.format(absline="", ang=repr(a0), mid=CHIC_MID,
                                        modeline=UNAVG))
    led = np.loadtxt(wd / "cled.ledger.txt")
    closure = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turnover = np.max(np.abs(led[:, 1])) + np.max(np.abs(led[:, 2]))
    check("unaveraged ledger closure across the chicane",
          float(np.max(np.abs(closure - closure[0])) / turnover), TOL_LEDGER)

    # ------------------------------------------------------------------
    print("== refusals ==")
    ok = refuse(exe, wd, "rf_open", LATCHIC.format(absline="", ang="1e-3", mid=OPEN_MID,
                                                   modeline=TRANSCRIBED), "NOT A CLOSED BUMP")
    print(f"--- refusal non-closed-bump break: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    ok = refuse(exe, wd, "rf_genb", LATCHIC.format(absline="", ang="1e-3", mid=CHIC_MID,
                                                   modeline=TRANSCRIBED),
                "GENESIS-MODEL interlude", imodel="genesis")
    print(f"--- refusal bend under genesis-model interludes: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    ok = refuse(exe, wd, "rf_zbig", LAT2SEG.format(absline="", zoff="0.35", gap="0.30",
                                                   modeline=TRANSCRIBED),
                "exceeds its upstream break")
    print(f"--- refusal z_offset exceeding the break: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    first = LAT2SEG.format(absline="", zoff="0", gap="0.30", modeline=TRANSCRIBED)
    first = first.replace("SEG: line = (UND, D, UND2)", "SEG: line = (UND2, D, UND)")
    first = first.replace("UND2: UND, z_offset = 0", "UND2: UND, z_offset = 0.01")
    ok = refuse(exe, wd, "rf_zfirst", first, "no upstream break")
    print(f"--- refusal z_offset on the first element: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
