#!/usr/bin/env python3
"""
Checks for the unaveraged mode (fel-physics.md sec-unaveraged). Self-referenced by design where the point is self-consistency (the
energy ledger, ballistic conservation) and closed-form where the point is measuring
the averaged mode's inputs (the coupling factor fc).

1. Energy ledger (check zero). A strong-seed helical probe: at every record,
   E_beam + U_field must be conserved -- the radiation kick and the source deposit are
   independent transcriptions of one wave equation, and only their consistency makes
   this hold. Check: max |d(E+U)| over the cumulative field-energy turnover.

2. Ballistic: a dark run's magnetic push does no work, exactly (gamma only changes in
   the radiation kick), and the emittance survives the RK4 push through the ramps.
   This also checks the handoff: with the sin^2 ramps the quiver vanishes at the
   segment ends, so the exit emittance equals the entry emittance; a hard-edge entry
   (fel_ramp_periods = -1, the explicit test sentinel) fails the orbit instrument loudly.

3. fc measured, both limits. Paired probes (12 and 20 periods, identical 2-period
   ramps): the difference of the two energy-modulation phasors
       F = (2/N) sum dgamma_j exp(+i theta0_j)
   isolates the flat region (ramps and their detuning cancel exactly), and
       fc_meas = |F_B - F_A| * beta0*gamma0*sqrt(2)*m_e / (|E0| * dL * sinc)
   must match the closed forms: planar fc = aw*(J0(xi)-J1(xi)), xi = aw^2/(2(1+aw^2));
   helical fc = aw. |E0| = sqrt(4 Z0 P / (pi w0^2)) is the Gaussian seed's on-axis
   envelope; sinc corrects the (tiny) off-resonance detuning.

4. h = 3: the same planar pair run with lambda0 = lambda1/3 measures the third-harmonic
   coupling against fc3 = aw*|J1(xi3)-J2(xi3)|, xi3 = (3/2) aw^2/(1+aw^2). The load is
   quiet at h = 3 because the beamlet quiet start cancels every harmonic below nbins
   (nbins = 8 here) -- the recorded reason a quadrature load is not needed for this.

5. Convergence: fc at 10/20/30 steps per period, tabulated (MINERVA's envelope).

6. Gain curve: the benchmark single segment, seeded steady state, unaveraged vs
   averaged from the same generated start. The exit ln power ratio prices the
   integrator-structure difference (ramps + split + quiver diagnostics)
   style: measured and bounded, not litigated.

Usage: check_unaveraged.py --exe <lucifer> --latdir <tests/bmad> --workdir <dir>
"""

from __future__ import annotations

import argparse
import math
import pathlib
import subprocess
import sys

import numpy as np

import beamio
from nml import to_groups
from scipy.special import jv

C_LIGHT = 2.99792458e8
M_ELECTRON = 510998.95069  # eV, Bmad's value
Z0 = 1.25663706127e-6 * C_LIGHT  # Bmad's mu_0_vac (2018 CODATA) times c

AW = 0.84853
LAMBDA_W = 0.015
GAMMA0 = 11357.82
L_A, L_B = 0.18, 0.30          # probe lengths [m]
N_RAMP = 2.0                   # ramp periods (each end)
DL_FLAT = L_B - L_A            # the isolated flat length [m]
LAMBDA1 = 1e-10
SEED_P = 1e3                   # small-signal seed power [W]
SEED_W0 = 4e-4                 # 1/e^2 intensity radius [m]

PROBE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = {lam}
  beam_init%n_particle = 2048
  beam_init%bunch_charge = {q}
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-07
  beam_init%a_norm_emit = 4e-9
  beam_init%b_norm_emit = 4e-9
  nbins = 8
  seed_power = {power}
  seed_waist_size = {w0}
  grid_n_pts = 129
  grid_half_width = 2e-3
  ran_seed = 4242
  write_initial = T
  write_diag = T
&end
"""

def probe_nml(wd, **kw):
    """probe with the steady-state charge derived (I = Q*c/spacing, spacing = lam) and
    the unaveraged mode selected by a wrapper lattice (attributes, not namelist)."""
    kw.setdefault("q", f"{3000 * float(kw['lam']) / 2.99792458e8:.12e}")
    kw["lat"] = unavg_wrapper(wd, kw["lat"], kw.pop("spp"), kw.pop("ramp"))
    return PROBE.format(**kw)

GAIN = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 129
  grid_half_width = 2e-4
  ran_seed = 4242
  write_diag = T
&end
"""

TDID = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 63
  grid_half_width = 2e-4
  window_length = 8e-10
  window_sample = 1
  shotnoise = T
  ran_seed = 999
  write_diag = T
&end
"""

def unavg_wrapper(wd, base, spp, ramp):
    """A wrapper lattice selecting the unaveraged mode with per-run parameters --
    the delz-sweep pattern: the mode and its knobs are lattice attributes."""
    name = f"w_{base.replace('.bmad','')}_s{spp}_r{str(ramp).replace('-','m').replace('.','p')}.bmad"
    (wd / name).write_text(
        f"call, file = {base}\n"
        f"fel_unaveraged = 1\n"
        f"wiggler::*[FEL_TRACKING] = fel_unaveraged\n"
        f"wiggler::*[FEL_STEPS_PER_PERIOD] = {spp}\n"
        f"wiggler::*[FEL_RAMP_PERIODS] = {ramp}\n")
    return name

FAILED = False


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="4"):
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-1000:]}")
        sys.exit(1)


def run_expect_refusal(exe, wd, name, text, fragment):
    """The run must fail, and by name."""
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    return r.returncode != 0 and fragment in r.stdout


SAND_WAKE_LAT = """call, file = unavg_sandwich.bmad
UNDW: undb, sr_wake = {amp_scale = 1, scale_with_length = T,
  longitudinal = {1e14, 0, 0, 0.25, none}}
SEGW: line = (UNDA, P1, UNDW, P1, UNDA)
use, SEGW
"""


def read_par(path, lam):
    """A dump -> dict of concatenated per-slice arrays (these probes are single slice).

    Read through beamio, which is told the run's wavelength: an openPMD beam file states
    the slice partition and not the radiation it was sliced on. These probes are steady
    state, so the slice spacing is one wavelength."""
    slices = beamio.read_slices(path, lam, lam)
    return {q: np.concatenate([sl[q] for sl in slices])
            for q in ("gamma", "theta", "x", "y", "px", "py")}


def phasor(root, wd, lam):
    """F = (2/N) sum dgamma * exp(+i theta0) between the initial and final dumps.

    theta0 comes from the initial dump, where the reference phase phi0 is still zero, so
    the phase here is the absolute one the tracker used."""
    p0 = read_par(wd / f"{root}-initial.beam.h5", lam)
    p1 = read_par(wd / f"{root}-final.beam.h5", lam)
    dg = p1["gamma"] - p0["gamma"]
    return 2.0 * np.mean(dg * np.exp(1j * p0["theta"]))


def fc_measured(fa, fb, lam, h):
    ks = 2 * math.pi / lam
    ku = 2 * math.pi / LAMBDA_W
    beta0 = math.sqrt(GAMMA0**2 - 1) / GAMMA0
    e0 = math.sqrt(4 * Z0 * SEED_P / (math.pi * SEED_W0**2))
    delta = h * ku - ks * (1 + AW**2) / (2 * GAMMA0**2)   # h-th resonance detuning [1/m]
    sinc = abs(math.sin(delta * DL_FLAT / 2) / (delta * DL_FLAT / 2)) if delta != 0 else 1.0
    return abs(fb - fa) * beta0 * GAMMA0 * math.sqrt(2) * M_ELECTRON / (e0 * DL_FLAT * sinc)


def fc_closed(h):
    xi = (h / 2) * AW**2 / (1 + AW**2)
    h0, h1 = (h - 1) // 2, (h + 1) // 2
    return AW * abs(jv(h0, xi) - jv(h1, xi))


def emit(p):
    x, px = p["x"], p["px"]
    vx = x - x.mean(); vp = px - px.mean()
    return math.sqrt(max((vx**2).mean() * (vp**2).mean() - (vx * vp).mean()**2, 0.0))


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
    for lat in ("unavg_probe_planar_a.bmad", "unavg_probe_planar_b.bmad",
                "unavg_probe_helical_a.bmad", "unavg_probe_helical_b.bmad",
                "aramis_1seg.bmad", "unavg_sandwich.bmad"):
        (wd / lat).write_bytes((latdir / lat).read_bytes())

    # 1. Energy ledger: strong-seed helical long probe (real turnover per record).
    run(exe, wd, "uv_ledger", probe_nml(wd, lat="unavg_probe_helical_b.bmad", root="uv_ledger",
        spp=20, ramp=N_RAMP, lam=LAMBDA1, power=1e8, w0=SEED_W0))
    led = np.loadtxt(wd / "uv_ledger.ledger.txt")
    etot = led[:, 1] + led[:, 2]
    turnover = np.abs(np.diff(led[:, 2])).sum()  # cumulative field-energy turnover |dU|
    ledger_dev = np.abs(etot - etot[0]).max() / max(turnover, 1e-300)
    # The /u_s deposit makes kick and source exact energy duals, so this closes at the
    # gamma<->pz round-trip floor (measured 1.0e-5, it was 6.7e-4 with Genesis's
    # averaged /gamma convention in the deposit).
    check("energy ledger: max|d(E_beam+U_field)| / sum|dU|", ledger_dev, 1e-4)
    # Internal consistency: the kick-side column must equal the realized beam change.
    # The floor is the per-record gamma <-> pz round trip (~1 ulp of gamma per particle
    # per record, summed). A missing or double-counted kick would sit at O(1).
    dE_beam = np.diff(led[:, 1])
    kick_col = led[1:, 3]
    check("ledger internal: |dE_kick - dE_beam| / max|dE|",
          np.abs(kick_col - dE_beam).max() / max(np.abs(dE_beam).max(), 1e-300), 1e-4)

    # 2. Ballistic dark run: B does no work. Ramps hand the emittance back.
    run(exe, wd, "uv_dark", probe_nml(wd, lat="unavg_probe_planar_b.bmad", root="uv_dark",
        spp=20, ramp=N_RAMP, lam=LAMBDA1, power=0.0, w0=SEED_W0))
    d0 = read_par(wd / "uv_dark-initial.beam.h5", LAMBDA1)
    d1 = read_par(wd / "uv_dark-final.beam.h5", LAMBDA1)
    check("ballistic: max|dgamma| (B does no work)",
          float(np.abs(d1["gamma"] - d0["gamma"]).max() / GAMMA0), 1e-12)
    check("ballistic: |emit_x out/in - 1| (ramp handoff)",
          abs(emit(d1) / emit(d0) - 1), 1e-6)
    # The coherent-quiver instruments. A hard-edge entry (fel_ramp_periods = -1, the
    # explicit test sentinel) starts the quiver about the wrong DC (pi = -a0 instead of 0), which
    # integrates to a centroid displacement ~a0*L/gamma -- measured 3.2e-5 m against
    # 1.9e-9 pristine, 19 sigma of this beam -- and a non-integer-period hard exit
    # would also leave the coherent quiver ~a0 in <px>. Both watched.
    check("handoff: |<x>_out - <x>_in| (dark) [m]",
          abs(float(d1["x"].mean() - d0["x"].mean())), 1e-7)
    check("handoff: |<px>_out - <px>_in| (dark, gamma*beta units)",
          abs(float(d1["px"].mean() - d0["px"].mean())), 1e-6)

    # 3. fc in both limits, differential pair measurement.
    results = {}
    for pol, lat_a, lat_b, lam, h in (
            ("planar", "unavg_probe_planar_a.bmad", "unavg_probe_planar_b.bmad", LAMBDA1, 1),
            ("helical", "unavg_probe_helical_a.bmad", "unavg_probe_helical_b.bmad", LAMBDA1, 1),
            ("planar_h3", "unavg_probe_planar_a.bmad", "unavg_probe_planar_b.bmad", LAMBDA1 / 3, 3)):
        for tag, lat in (("a", lat_a), ("b", lat_b)):
            run(exe, wd, f"uv_{pol}_{tag}", probe_nml(wd, lat=lat, root=f"uv_{pol}_{tag}",
                spp=20, ramp=N_RAMP, lam=lam, power=SEED_P, w0=SEED_W0))
        fa = phasor(f"uv_{pol}_a", wd, lam)
        fb = phasor(f"uv_{pol}_b", wd, lam)
        fm = fc_measured(fa, fb, lam, h)
        fx = AW if (pol == "helical") else fc_closed(h)
        results[pol] = (fm, fx)
        check(f"fc measured vs closed form, {pol} (h={h}): |ratio-1|  [{fm:.5f} vs {fx:.5f}]",
              abs(fm / fx - 1), 5e-3)

    # 5. Convergence over MINERVA's envelope (planar pair at 10/20/30 steps/period).
    print("--- step-size convergence (planar fc, steps/period):")
    fcs = {}
    for spp in (10, 30):
        for tag, lat in (("a", "unavg_probe_planar_a.bmad"), ("b", "unavg_probe_planar_b.bmad")):
            run(exe, wd, f"uv_cv{spp}_{tag}", probe_nml(wd, lat=lat, root=f"uv_cv{spp}_{tag}",
                spp=spp, ramp=N_RAMP, lam=LAMBDA1, power=SEED_P, w0=SEED_W0))
        fcs[spp] = fc_measured(phasor(f"uv_cv{spp}_a", wd, LAMBDA1),
                               phasor(f"uv_cv{spp}_b", wd, LAMBDA1), LAMBDA1, 1)
    fcs[20] = results["planar"][0]
    for spp in (10, 20, 30):
        print(f"      {spp:3d} steps/period: fc = {fcs[spp]:.6f}")
    check("convergence: |fc(30)/fc(20) - 1|", abs(fcs[30] / fcs[20] - 1), 5e-4)

    # 6. Gain curve: benchmark segment, unaveraged vs averaged, same start.
    run(exe, wd, "uv_gain_unavg", GAIN.format(root="uv_gain_unavg",
        lat=unavg_wrapper(wd, "aramis_1seg.bmad", 20, 2)))
    run(exe, wd, "uv_gain_avg", GAIN.format(root="uv_gain_avg", lat="aramis_1seg.bmad"))
    pu = np.loadtxt(wd / "uv_gain_unavg.diag.txt")[:, 2]
    pa = np.loadtxt(wd / "uv_gain_avg.diag.txt")[:, 2]
    lnr = abs(math.log(pu[-1] / pa[-1]))
    check(f"gain curve: |ln(P_unavg/P_avg)| at segment exit  [{pu[-1]:.4e} vs {pa[-1]:.4e}]",
          lnr, 0.2, note="(priced integrator-structure difference)")

    # 7. Mixed line (Stage A): averaged / unaveraged / averaged sandwich with pipe
    # interludes. Completing at all exercises the convention-flag asserts at real
    # internal boundaries. The ledger must be confined to and conserved over the
    # unaveraged segment. A wake on that segment must refuse by name, and the exit
    # power is priced against the all-averaged twin.
    (wd / "sandwich_avg.bmad").write_text(
        "call, file = unavg_sandwich.bmad\nfel_averaged = 0\n"
        "UNDB[FEL_TRACKING] = fel_averaged\n")
    run(exe, wd, "uv_sand", GAIN.format(root="uv_sand", lat="unavg_sandwich.bmad"))
    run(exe, wd, "uv_sand_avg", GAIN.format(root="uv_sand_avg", lat="sandwich_avg.bmad"))
    led = np.loadtxt(wd / "uv_sand.ledger.txt")
    confined = float(((led[:, 0] > 0.81 - 1e-9) & (led[:, 0] < 1.41 + 1e-9)).all())
    check("sandwich: ledger rows confined to the unaveraged segment (1 = yes)",
          1.0 - confined, 0.5)
    etot = led[:, 1] + led[:, 2]
    turn = np.abs(np.diff(led[:, 2])).sum()
    check("sandwich: ledger max|d(E+U)| / sum|dU| on the middle segment",
          np.abs(etot - etot[0]).max() / max(turn, 1e-300), 1e-3)
    ps = np.loadtxt(wd / "uv_sand.diag.txt")[:, 2]
    pa2 = np.loadtxt(wd / "uv_sand_avg.diag.txt")[:, 2]
    check(f"sandwich: |ln(P_mixed/P_averaged)| at exit  [{ps[-1]:.4e} vs {pa2[-1]:.4e}]",
          abs(math.log(ps[-1] / pa2[-1])), 5e-2, note="(one segment's ramp+mode price)")
    (wd / "sandwich_wake.bmad").write_text(SAND_WAKE_LAT)
    refused = run_expect_refusal(exe, wd, "uv_sandw",
        GAIN.format(root="uv_sandw", lat="sandwich_wake.bmad"),
        "ELEMENT SR WAKES ARE NOT SUPPORTED IN THE UNAVERAGED MODE")
    check("sandwich: wake on the unaveraged segment refused by name (1 = yes)",
          0.0 if refused else 1.0, 0.5)

    # 8. Thread invariance: the parallel slice loop must be invisible. A multi-slice
    # time-dependent unaveraged run (slippage active, so slices genuinely interleave
    # through the field) at 1 thread and at 8 threads must produce byte-identical
    # diagnostics AND ledger -- the same guarantee the averaged path carries
    # (per-slice private state, fixed-order energy reduction).
    wl = unavg_wrapper(wd, "aramis_1seg.bmad", 20, 2)
    run(exe, wd, "uv_tid1", TDID.format(root="uv_tid1", lat=wl), threads="1")
    run(exe, wd, "uv_tid8", TDID.format(root="uv_tid8", lat=wl), threads="8")
    same = all((wd / f"uv_tid1{s}").read_bytes() == (wd / f"uv_tid8{s}").read_bytes()
               for s in (".diag.txt", ".ledger.txt"))
    check("thread invariance: 1-thread vs 8-thread TD run byte-identical (1 = yes)",
          0.0 if same else 1.0, 0.5)

    # 9. The TIME-DEPENDENT ledger closure. The window is an open system -- slippage
    # transmits the head slice's light out of the simulation -- and the deposit's own
    # |dE_src|^2 is the one field-energy term the kick/deposit duality does not charge
    # to the beam (physically: the substep's spontaneous emission). Both are banked as
    # ledger columns, so the closing quantity is exact:
    #     E_beam + U_window + U_escaped - U_spont = const.
    # Wakes would be a second, unbanked exit channel from the beam. They are refused in
    # this mode, which is what entitles this check to exist.
    led = np.loadtxt(wd / "uv_tid1.ledger.txt")
    etot = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turn = np.abs(np.diff(led[:, 2])).sum() + abs(led[-1, 4]) + abs(led[-1, 5]) + abs(led[-1, 6])
    check("TD ledger: max|d(E_beam+U_window+U_escaped-U_spont+E_rad)| / turnover",
          np.abs(etot - etot[0]).max() / max(turn, 1e-300), 1e-3)

    if FAILED:
        print("unaveraged checks: FAIL")
        sys.exit(1)
    print("unaveraged checks: PASS")


if __name__ == "__main__":
    main()
