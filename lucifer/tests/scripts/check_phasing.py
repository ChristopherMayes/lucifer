#!/usr/bin/env python3
"""
Phasing checks (fel-physics.md sec-phasing): how the beam-field phase behaves between
undulator segments, held by symmetries, closed forms and Genesis itself.

  1. Re-anchor baseline (relative mode, the default): scanning an inter-segment gap
     by fractions of 2 gamma^2 lambda leaves the bunching phase entering the next
     segment flat -- the transcription of Genesis's drift autophasing (measured flat
     on Genesis itself the same way; the residual is a slow adiabatic trend).
  2. The Z_OFFSET knob: displacing the downstream wiggler by delta (standard Bmad
     misalignment, anchor at the nominal position) shifts that phase by exactly
     -2 pi delta / (2 gamma^2 lambda) -- the analytic drift slip rate, no fit.
  3. Cross-mode identity: in absolute mode (bmad_com[absolute_time_tracking] = T in
     the lattice, honored through Bmad's own resolver) the same phase arrives via the
     real gap length instead: the absolute-mode gap scan must reproduce the relative-
     mode knob scan point by point.
  4. Phase-shifter parity vs Genesis4: Genesis scans PHASESHIFTER phi (which needs
     finite length to register -- a zero-length one silently does nothing); we scan
     z_offset with delta = phi * 2 gamma^2 lambda / (2 pi). Same physical scan; the
     bunching-phase curves must agree point by point. This anchors the sign
     convention ("delay goes backwards") against Genesis's own element.
  5. Chicane: a four-bend closed bump between segments (Bmad seam tracks the beam;
     the radiation drifts the chord from ele%floor). Relative mode drops the
     geometric fraction (Genesis's "chicane is always autophasing"); absolute mode
     ramps with the bend angle at the slope an independent geometric computation of
     d(arc - chord)/d(angle) predicts -- itself cross-checked against the textbook
     small-angle path lengthening theta^2 (2 L_bend/3 + L_drift). The unaveraged
     ledger closes across a chicane sandwich.
  6. Time-dependent chicane: the window rotations the geometric delay buys (the
     steady-state sections above never exercise them -- slippage is a no-op with one
     slice). A delay of a few wavelengths must bank exactly floor(delay/lambda) more
     escaped slices than the straight-line twin of the same arc length, and the run
     must be byte-identical at 1 and 8 threads.
  7. Refusals by name: a non-closed-bump break; a bend under the genesis-model
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

import convert_genesis
from nml import to_groups

import pool

FAILED = False

TWO_G2L = 2 * 11357.82**2 * 1e-10       # 2 gamma^2 lambda: one full turn of gap phase.
LAM = 1e-10

TOL_FLAT = 5e-3          # measured 1.5e-3 span: the re-anchor baseline's adiabatic trend.
TOL_KNOB = 5e-3          # measured 1.5e-3: knob ramp vs the exact analytic slope.
TOL_XMODE = 5e-3         # measured 1.2e-3: absolute gap scan vs relative knob scan.
TOL_PARITY = 1e-6        # measured 1.9e-8 from shared dumps: knob curve vs PHASESHIFTER.
TOL_CHIC_SLOPE = 2e-3    # measured 6.8e-4: absolute chicane ramp vs the geometric slope
                         #   (583 wavelengths of delay at the base angle).
TOL_LEDGER = 1e-4        # measured 4.0e-6: unaveraged ledger closure across the chicane.
TOL_CF = 1e-6            # measured 6.4e-8: the 2D trace vs the small-angle closed form.
TOL_POWER = 1e-6         # measured 9.3e-9 from shared dumps: exit power vs phi, our knob
                         #   vs Genesis's shifter (1.4e-4 on independently loaded beams --
                         #   that measures the loaders, so this tier shares the start).

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

NML = """! flat keys; routed into the three groups by nml.to_groups
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
  write_diag = T
/
"""

LATSTRAIGHT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
{absline}
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
DD: pipe, l = 0.025
STR: pipe, l = 0.275               ! ARC-MATCHED to the bump (5*0.025 + 4*0.05 = 0.325
                                   ! total with the flanking DDs), so the two runs share
                                   ! the autophase term and only the geometric delay differs.
SEG: line = (UND, DD, STR, DD, UND)
use, SEG
{modeline}
"""

NML_TD = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{tag}.bmad"
  out_root = "{tag}"
  lambda0 = 1e-10
  interlude_model = 'bmad'
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 1.6e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -8e-10
  beam_init%grid(3)%x_max = 8e-10
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 1e6
  seed_waist_size = 30e-6
  grid_n_pts = 63
  grid_half_width = 2e-4
  window_length = 1.6e-9
  window_sample = 1
  ran_seed = 777
  keep_escaped_field = T
  write_diag = T
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

GEN_PREP = """&setup
rootname=psp
lattice=ps0.lat
beamline=SEG
lambda0=1e-10
gamma0=11357.82
delz=0.015
shotnoise=0
nbins = 8
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

&write
field = PSP-initial
beam = PSP-initial
&end
"""

GEN_IMPORT = """&setup
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

&importfield
file = PSP-initial.fld.h5
&end

&importbeam
file = PSP-initial.par.h5
&end

&track
fft_fieldsolver = true
&end
"""

# lambda0 and nbins are the deck's: an openPMD beam file carries the slice partition and
# not the radiation it was sliced on. Both match GEN_PREP, which writes the shared dumps.
NML_IMP = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{tag}.bmad"
  out_root = "{tag}"
  beam_file = "PSP-initial.beam.h5"
  field_file = "PSP-initial.wf.h5"
  lambda0 = 1e-10
  nbins = 8
  interlude_model = '{imodel}'
  write_diag = T
/
"""

GEN_LAT_PS = """UND: undulator = {{ lambdau=0.015000, nwig=30, aw=0.84853, helical=True}};
D1: drift = {{ l = 0.1425 }};
PS: PHASESHIFTER = {{ l = 0.015, phi = {phi} }};
D2: drift = {{ l = 0.1425 }};
SEG: line={{UND,D1,PS,D2,UND}};
"""


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def wrap(dp):
    return (np.asarray(dp) + np.pi) % (2 * np.pi) - np.pi


def run(exe, wd, tag, lat_text, imodel="bmad", extra="", nml_t=None, threads="4"):
    (wd / f"{tag}.bmad").write_text(lat_text)
    nml = (nml_t or NML).format(tag=tag, imodel=imodel)
    if extra:
        nml = nml.replace("/\n", extra + "/\n")
    (wd / f"{tag}.nml").write_text(to_groups(nml))
    r = subprocess.run([str(exe), f"{tag}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {tag} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)


def refuse(exe, wd, tag, lat_text, fragment, imodel="bmad"):
    (wd / f"{tag}.bmad").write_text(lat_text)
    (wd / f"{tag}.nml").write_text(to_groups(NML.format(tag=tag, imodel=imodel)))
    r = subprocess.run([str(exe), f"{tag}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    return r.returncode != 0 and fragment in (r.stdout + r.stderr)


def bphase(wd, tag, z_probe):
    d = np.loadtxt(wd / f"{tag}.diag.txt").reshape(-1, 1, 12)
    z = d[:, 0, 0]
    return float(d[np.searchsorted(z, z_probe), 0, 5])


def exit_power(wd, tag):
    d = np.loadtxt(wd / f"{tag}.diag.txt").reshape(-1, 1, 12)
    return float(d[-1, 0, 2])


def escaped_count(wd, tag):
    with h5py.File(wd / f"{tag}-escaped.fld.h5") as h5:
        return int(h5["slicecount"][0])


def chicane_delay(ang):
    """arc - chord of the four-bend closed bump, exact 2D geometry, independently of
    the walk's floor arithmetic: bends of arc length lb and bend angle +a,-a,-a,+a
    with drifts ld between (and the outer DD pieces adding straight length only)."""
    lb, ld = 0.05, 0.025
    # Trace the reference: (x, z, angle). Bend of angle a: chord lb*sin(a/2)/(a/2) at
    # heading (angle_in + a/2). Drift: ld at current heading.
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
    ap.add_argument("--pyrepo", default=convert_genesis.DEFAULT_PYREPO,
                    help="openPMD-beamphysics checkout, for the dump conversion")
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(args.exe).resolve()

    NSCAN = 6
    ABS = "bmad_com[absolute_time_tracking] = T"

    # ------------------------------------------------------------------
    print("== re-anchor baseline + z_offset knob + cross-mode identity ==")

    # The scan points are independent processes. The pool runs them a few at a
    # time (same configurations, same thread counts, same tolerances).
    jobs = []
    for k in range(NSCAN):
        f = k / 8
        jobs.append(lambda k=k, f=f: run(exe, wd, f"pb{k}", LAT2SEG.format(absline="",
            zoff="0", gap=f"{0.30 + f * TWO_G2L:.12e}", modeline=TRANSCRIBED), imodel="genesis"))
        jobs.append(lambda k=k, f=f: run(exe, wd, f"pk{k}", LAT2SEG.format(absline="",
            zoff=f"{f * TWO_G2L:.12e}", gap="0.30", modeline=TRANSCRIBED), imodel="genesis"))
        jobs.append(lambda k=k, f=f: run(exe, wd, f"px{k}", LAT2SEG.format(absline=ABS,
            zoff="0", gap=f"{0.30 + f * TWO_G2L:.12e}", modeline=TRANSCRIBED), imodel="genesis"))
    pool.run_all(jobs, threads_per_job=4)
    base_bp = [bphase(wd, f"pb{k}", 0.76) for k in range(NSCAN)]
    knob_bp = [bphase(wd, f"pk{k}", 0.76) for k in range(NSCAN)]
    knob_pw = [exit_power(wd, f"pk{k}") for k in range(NSCAN)]
    xmode_bp = [bphase(wd, f"px{k}", 0.76) for k in range(NSCAN)]

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
    print("== phase-shifter parity vs Genesis4 (SHARED initial dumps) ==")

    # Both codes start from the same particles and field: Genesis writes the dumps
    # once (the loaders are independently written, so comparing power on
    # independently generated beams would measure the loaders, not the phasing --
    # measured at 1.4e-4 that way, against a 2.2e-2 power swing).

    (wd / "ps0.lat").write_text(GEN_LAT_PS.format(phi="0"))
    (wd / "psp.in").write_text(GEN_PREP)
    r = subprocess.run([args.genesis, "psp.in"], cwd=wd, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "FI_PROVIDER": "tcp"})
    if r.returncode != 0:
        print(f"FAIL: genesis prep run:\n{r.stdout[-800:]}")
        sys.exit(1)

    # The tracker reads openPMD only, so the shared dumps convert once here. Genesis keeps
    # reading its own, so both codes still start from the same particles and field.
    for src, dst in (("PSP-initial.par.h5", "PSP-initial.beam.h5"),
                     ("PSP-initial.fld.h5", "PSP-initial.wf.h5")):
        convert_genesis.to_openpmd(wd / src, wd / dst, args.pyrepo)

    def gen_ps(k):
        (wd / f"ps{k}.lat").write_text(GEN_LAT_PS.format(phi=f"{2 * np.pi * k / 8:.12e}"))
        (wd / f"ps{k}.in").write_text(GEN_IMPORT.format(root=f"PS{k}", lat=f"ps{k}.lat"))
        r = subprocess.run([args.genesis, f"ps{k}.in"], cwd=wd, capture_output=True, text=True,
                           env={"PATH": "/usr/bin:/bin", "FI_PROVIDER": "tcp"})
        if r.returncode != 0:
            print(f"FAIL: genesis ps{k}:\n{r.stdout[-800:]}")
            sys.exit(1)

    jobs = []
    for k in range(NSCAN):
        jobs.append(lambda k=k: run(exe, wd, f"pi{k}", LAT2SEG.format(absline="",
            zoff=f"{k / 8 * TWO_G2L:.12e}", gap="0.30", modeline=TRANSCRIBED),
            imodel="genesis", nml_t=NML_IMP))
        jobs.append(lambda k=k: gen_ps(k))
    pool.run_all(jobs, threads_per_job=4)

    knob_bp2, knob_pw2 = [], []
    gen_bp, gen_pw = [], []
    for k in range(NSCAN):
        knob_bp2.append(bphase(wd, f"pi{k}", 0.76))
        knob_pw2.append(exit_power(wd, f"pi{k}"))
        with h5py.File(wd / f"PS{k}.out.h5") as h5:
            z = h5["Lattice/zplot"][:].ravel()
            aw = h5["Lattice/aw"][:].ravel()
            bp = h5["Beam/bunchingphase"][:].ravel()
            pw = h5["Field/power"][:].ravel()
        n = min(len(z), len(aw), len(bp))
        i2 = np.where((aw[:n] > 0) & (z[:n] > 0.5))[0][0] + 1
        gen_bp.append(float(bp[i2]))
        gen_pw.append(float(pw[-1]))
    d_gen = wrap(np.array(gen_bp) - gen_bp[0])
    d_knob2 = wrap(np.array(knob_bp2) - knob_bp2[0])
    check("our z_offset curve == Genesis PHASESHIFTER curve",
          float(np.max(np.abs(wrap(d_knob2 - d_gen)))), TOL_PARITY,
          note="[delta = phi 2 gamma^2 lambda / 2pi]")

    # The phase must reach the physics identically, not just the bookkeeping: exit
    # power against phi, code vs code, from the shared start.

    gp = np.array(gen_pw);  kp = np.array(knob_pw2)
    check("exit power vs phi: our knob == Genesis shifter (max rel)",
          float(np.max(np.abs(kp - gp) / gp)), TOL_POWER,
          note=f"[power swing over the scan {np.ptp(gp)/gp.mean():.2e}]")

    # ------------------------------------------------------------------
    print("== chicane: relative flat, absolute at the geometric slope, ledger ==")
    a0 = 1e-3
    da = LAM / 8 / (2 * a0 * 0.058)     # ~lam/8 of delay per step at the true L_eff scale
    jobs = []
    for k in range(NSCAN):
        a = a0 + k * da
        jobs.append(lambda k=k, a=a: run(exe, wd, f"cr{k}", LATCHIC.format(absline="",
            ang=f"{a:.12e}", mid=CHIC_MID, modeline=TRANSCRIBED)))
        jobs.append(lambda k=k, a=a: run(exe, wd, f"ca{k}", LATCHIC.format(absline=ABS,
            ang=f"{a:.12e}", mid=CHIC_MID, modeline=TRANSCRIBED)))
    jobs.append(lambda: run(exe, wd, "cled", LATCHIC.format(absline="", ang=f"{a0:.12e}",
        mid=CHIC_MID, modeline=UNAVG)))
    pool.run_all(jobs, threads_per_job=4)
    rel_bp = [bphase(wd, f"cr{k}", 0.83) for k in range(NSCAN)]
    abs_bp = [bphase(wd, f"ca{k}", 0.83) for k in range(NSCAN)]
    delays = [chicane_delay(a0 + k * da) for k in range(NSCAN)]

    check("chicane, relative mode: geometric fraction dropped (flat, rad span)",
          float(np.ptp(wrap(np.array(rel_bp) - rel_bp[0]))), TOL_FLAT)

    d_abs = wrap(np.array(abs_bp) - abs_bp[0])
    d_pred = wrap(-2 * np.pi * (np.array(delays) - delays[0]) / LAM)
    check("chicane, absolute mode vs independent geometry (rad)",
          float(np.max(np.abs(wrap(d_abs - d_pred)))), TOL_CHIC_SLOPE,
          note=f"[delay(a0) = {delays[0]:.3e} m = {delays[0]/LAM:.0f} wavelengths]")

    # The 2D trace itself against the textbook small-angle path lengthening
    # theta^2 (2 L_bend/3 + L_drift): two independent derivations of the same delay.

    cf = np.array([a0 + k * da for k in range(NSCAN)]) ** 2 * (2 * 0.05 / 3 + 0.025)
    check("2D trace vs the small-angle closed form", float(np.max(np.abs(np.array(delays) / cf - 1))),
          TOL_CF)

    # The unaveraged ledger closes across a chicane sandwich (energy bookkeeping
    # survives the seam detour, and the ledger's columns cover the unaveraged segments).
    # The cled run itself went through the pool above.
    led = np.loadtxt(wd / "cled.ledger.txt")
    closure = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turnover = np.max(np.abs(led[:, 1])) + np.max(np.abs(led[:, 2]))
    check("unaveraged ledger closure across the chicane",
          float(np.max(np.abs(closure - closure[0])) / turnover), TOL_LEDGER)

    # ------------------------------------------------------------------
    # The window rotations the geometric delay buys. Steady state never exercises
    # them (slippage is a no-op with one slice), so this runs time dependent with a
    # bend angle tuned for a few wavelengths of delay: the chicane must bank exactly
    # floor(delay/lambda) more escaped slices than its straight-line twin of the same
    # arc length, and the run must be thread-invariant.

    print("== time-dependent chicane: window rotations and threads ==")

    # Delays chosen MID-INTERVAL (x.5 wavelengths), never on an integer: a target of
    # exactly 3.000 lambda sits on the floor boundary, where the 6e-8 difference
    # between the exact trace and the small-angle form decides between 2 and 3 (the
    # boundary was located: 2.99 -> 2, 3.01 -> 3, which is itself an independent
    # confirmation of the walk's floor geometry against the trace at this scale).

    a_td = float(np.sqrt(3.5 * LAM / (2 * 0.05 / 3 + 0.025)))
    jobs = [lambda: run(exe, wd, "tds", LATSTRAIGHT.format(absline="", modeline=TRANSCRIBED),
                        nml_t=NML_TD),
            lambda: run(exe, wd, "tdc", LATCHIC.format(absline="", ang=f"{a_td:.12e}",
                        mid=CHIC_MID, modeline=TRANSCRIBED), nml_t=NML_TD),
            lambda: run(exe, wd, "tdc8", LATCHIC.format(absline="", ang=f"{a_td:.12e}",
                        mid=CHIC_MID, modeline=TRANSCRIBED), nml_t=NML_TD, threads="8")]
    for target in (3.5, 6.5, 10.5):
        a = float(np.sqrt(target * LAM / (2 * 0.05 / 3 + 0.025)))
        jobs.append(lambda a=a, tag=f"tdc{int(target)}": run(exe, wd, tag,
            LATCHIC.format(absline="", ang=f"{a:.12e}", mid=CHIC_MID,
                           modeline=TRANSCRIBED), nml_t=NML_TD))
    pool.run_all(jobs, threads_per_job=4)

    n_str = escaped_count(wd, "tds")
    worst, detail = 0.0, []
    for target in (3.5, 6.5, 10.5):
        a = float(np.sqrt(target * LAM / (2 * 0.05 / 3 + 0.025)))
        n_expect = int(np.floor(chicane_delay(a) / LAM))
        got = escaped_count(wd, f"tdc{int(target)}") - n_str
        worst = max(worst, abs(got - n_expect))
        detail.append(f"{n_expect}:{got}")
    check("extra banked slices == floor(delay/lambda), three delays", float(worst), 0.5,
          note="[expected:got " + " ".join(detail) + "]")
    same = (wd / "tdc.diag.txt").read_bytes() == (wd / "tdc8.diag.txt").read_bytes()
    with h5py.File(wd / "tdc-final.wf.h5") as a, h5py.File(wd / "tdc8-final.wf.h5") as b:
        names = []
        a.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
        same = same and bool(names)
        for n in names:
            same = same and bool(np.array_equal(a[n][...], b[n][...]))
    check("TD chicane 1 vs 8 threads byte/dataset-identical", 0.0 if same else 1.0, 0.5)

    # ------------------------------------------------------------------
    print("== refusals ==")
    first = LAT2SEG.format(absline="", zoff="0", gap="0.30", modeline=TRANSCRIBED)
    first = first.replace("SEG: line = (UND, D, UND2)", "SEG: line = (UND2, D, UND)")
    first = first.replace("UND2: UND, z_offset = 0", "UND2: UND, z_offset = 0.01")
    oks = pool.run_all([
        lambda: refuse(exe, wd, "rf_open", LATCHIC.format(absline="", ang="1e-3",
                       mid=OPEN_MID, modeline=TRANSCRIBED), "NOT A CLOSED BUMP"),
        lambda: refuse(exe, wd, "rf_genb", LATCHIC.format(absline="", ang="1e-3",
                       mid=CHIC_MID, modeline=TRANSCRIBED),
                       "GENESIS-MODEL INTERLUDE", imodel="genesis"),
        lambda: refuse(exe, wd, "rf_zbig", LAT2SEG.format(absline="", zoff="0.35",
                       gap="0.30", modeline=TRANSCRIBED), "EXCEEDS ITS UPSTREAM BREAK"),
        lambda: refuse(exe, wd, "rf_zfirst", first, "NO UPSTREAM BREAK"),
    ], threads_per_job=4)
    for name, ok in zip(("non-closed-bump break", "bend under genesis-model interludes",
                         "z_offset exceeding the break", "z_offset on the first element"), oks):
        print(f"--- refusal {name}: {'ok' if ok else '** FAIL **'}")
        FAILED = FAILED or not ok

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
