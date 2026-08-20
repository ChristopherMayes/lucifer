#!/usr/bin/env python3
"""
Spontaneous emission: the FEL modes against Bmad's own radiation physics and against
the analytic undulator-radiation rate. Bmad-only -- Genesis is not involved, because
the averaged mode's agreement with Genesis is already established by the tiers.

ONE short helical wiggler (spont_probe.bmad, 40 periods, benchmark parameters) tracked
four ways from the same beam:

  1. Bmad's own runge_kutta with radiation damping (spont_probe_rk.bmad: the same
     wiggler with tracking_method = runge_kutta, so the FEL walk hands it to
     track1_bunch -- no FEL model, no grid, no SVEA band). This is the independent
     implementation, and it must reproduce the analytic rate.
  2. The same, radiation OFF: the integrator must lose exactly nothing.
  3. The averaged FEL mode: loses essentially NOTHING, by design. The KMR/SVEA step
     adds 2S to the field while kicking particles with E, so the |S|^2 part of the
     field energy is created rather than taken from the beam. This is a documented
     model property, not a defect -- Genesis carries an optional &sponrad module to
     add the missing loss by hand, off by default. The check pins it as a known zero.
  4. The unaveraged FEL mode: it conserves energy by construction (the /u_s deposit
     makes kick and source exact duals), so its beam DOES pay -- but only for the
     radiation the grid can hold. An SVEA grid represents angles up to the FFT
     Nyquist, theta_max = lambda/(2 dx), so only a few percent of the emission is
     captured, and THAT is what the beam is charged.

The analytic rate is the classical result dgamma/ds = (2/3) r_e gamma^2 ku^2 aw^2
(aw rms, both polarizations) -- identical to the coefficient in Genesis's own
Incoherent.cpp, 1.88e-15 * (ku gamma aw)^2.

The physics check on (4) is the SCALING: vary only the grid's angular acceptance
(ngrid at fixed box) and the captured loss must track the acceptance the way undulator
radiation does. The absolute normalization is compared against a dipole-limit estimate
of the angular distribution, which is only good to a factor of a few at aw ~ 1 -- so
the shape is checked tightly and the magnitude loosely, which is the honest split.

Run by the benchmark harness; exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

import numpy as np

M_ELECTRON = 0.51099895069e6
R_E = 2.8179403262e-15
GAMMA0 = 11357.82
LAMBDA_U = 0.015
AW = 0.84853
L_UND = 0.60          # spont_probe.bmad
DGRID = 2e-4          # grid half width used below
NGRID_REF = 255

FAILED = False

NML = """&fel_track_params
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 2.4017e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.2e-9
  beam_init%grid(3)%x_max = 1.2e-9
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 0
  seed_waist_size = 30e-6
  grid_n_pts = {ngrid}
  grid_half_width = 2e-4
  window_length = 2.4e-9
  window_sample = 1
  shotnoise = T
  ran_seed = 4242
{extra}&end
"""

WRAP = """call, file = spont_probe.bmad
fel_unaveraged = 1
wiggler::*[FEL_TRACKING] = fel_unaveraged
"""


def check(name, value, lo, hi, note=""):
    """Two-sided: lo <= value <= hi."""
    global FAILED
    ok = lo <= value <= hi
    print(f"--- {name}: {value:.4e} (check {lo:.2e} .. {hi:.2e}) {note} "
          f"{'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text):
    (wd / (name + ".nml")).write_text(text)
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
        sys.exit(1)


def loss(wd, root):
    """-Delta<gamma> over the wiggler, slice-averaged, from the diag file."""
    fn = wd / f"{root}.diag.txt"
    ns = int(re.search(r"nslice = (\d+)", fn.open().read(400)).group(1))
    d = np.loadtxt(fn).reshape(-1, ns, 12)
    return -(d[-1, :, 6] - d[0, :, 6]).mean() / M_ELECTRON


def dipole_fraction(u):
    """Fraction of a dipole-limit angular distribution (1+u^2)/(1+u)^4 inside u,
    u = gamma^2 theta^2 / (1 + aw^2). Closed form; exact for aw << 1, indicative
    at aw ~ 1 (which is why only the SHAPE is checked tightly)."""
    x = 1.0 + u
    return 1.5 * (2 / 3 - 1 / x + 1 / x**2 - (2 / 3) / x**3)


def acceptance(ngrid):
    dx = 2 * DGRID / (ngrid - 1)
    theta = 1e-10 / (2 * dx)                     # FFT Nyquist angle
    return theta, (GAMMA0 * theta) ** 2 / (1 + AW**2)


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
    for lat in ("spont_probe.bmad", "spont_probe_rk.bmad"):
        (wd / lat).write_bytes((latdir / lat).read_bytes())
    (wd / "sp_uv.bmad").write_text(WRAP)

    ku = 2 * np.pi / LAMBDA_U
    analytic = (2 / 3) * R_E * GAMMA0**2 * ku**2 * AW**2 * L_UND
    print(f"--- analytic spontaneous loss over {L_UND} m: {analytic:.6f} in gamma "
          f"({analytic/L_UND:.5f}/m; = Genesis Incoherent.cpp doLoss)")

    # 1-2. Bmad's own tracking through the same wiggler, radiation on and off.
    run(exe, wd, "sp_rk", NML.format(lat="spont_probe_rk.bmad", root="sp_rk",
                                     ngrid=NGRID_REF,
                                     extra="  radiation_damping = T\n  reference_run = T\n"))
    run(exe, wd, "sp_rk0", NML.format(lat="spont_probe_rk.bmad", root="sp_rk0",
                                      ngrid=NGRID_REF, extra="  reference_run = T\n"))
    l_rk, l_rk0 = loss(wd, "sp_rk"), loss(wd, "sp_rk0")
    check("Bmad runge_kutta + radiation vs analytic, |ratio - 1|",
          abs(l_rk / analytic - 1), 0.0, 5e-3,
          note="(independent implementation of the same physics)")
    check("Bmad runge_kutta, radiation OFF: |d gamma| (integrator conserves)",
          abs(l_rk0), 0.0, 1e-9)

    # 3. The averaged mode: a known, documented zero.
    run(exe, wd, "sp_avg", NML.format(lat="spont_probe.bmad", root="sp_avg",
                                      ngrid=NGRID_REF, extra=""))
    f_avg = abs(loss(wd, "sp_avg")) / analytic
    check("averaged mode: beam debit / analytic (KMR does NOT charge the beam)",
          f_avg, 0.0, 2e-2, note="(model property; Genesis has &sponrad for it)")

    # 4. The unaveraged mode: grid-acceptance-limited, and it must SCALE that way.
    print("--- unaveraged mode, grid angular-acceptance scan (box fixed):")
    meas, pred = {}, {}
    for ngrid in (127, NGRID_REF, 511):
        root = f"sp_uv{ngrid}"
        run(exe, wd, root, NML.format(lat="sp_uv.bmad", root=root, ngrid=ngrid, extra=""))
        theta, u = acceptance(ngrid)
        meas[ngrid] = loss(wd, root)
        pred[ngrid] = dipole_fraction(u) * analytic
        print(f"      ngrid {ngrid:4d}: theta_max {theta:.3e} rad (gamma*theta "
              f"{GAMMA0*theta:.3f}), captured {meas[ngrid]:.3e} = "
              f"{100*meas[ngrid]/analytic:5.2f}% of analytic "
              f"(dipole estimate {100*pred[ngrid]/analytic:5.2f}%)")

    check("unaveraged: captured fraction of the analytic rate at the reference grid",
          meas[NGRID_REF] / analytic, 1e-2, 8e-2,
          note="(only what the SVEA grid can hold)")
    # The shape test: a 16x range in captured solid angle.
    r_meas = meas[511] / meas[127]
    r_pred = pred[511] / pred[127]
    check("unaveraged: acceptance SCALING, measured/predicted ratio over 16x solid angle",
          r_meas / r_pred, 0.5, 2.0,
          note=f"(measured {r_meas:.2f}x vs predicted {r_pred:.2f}x)")

    if FAILED:
        print("SPONTANEOUS CHECKS: FAIL")
        sys.exit(1)
    print("spontaneous checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
