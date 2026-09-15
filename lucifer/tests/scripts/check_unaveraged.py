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
   (unaveraged_ramp_periods = -1, the explicit sentinel) fails the orbit instrument loudly.

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
   (beamlet_size = 8 here) -- the recorded reason a quadrature load is not needed for this.

5. Convergence: fc at 10/20/30 steps per period, tabulated. This is what the floor of 10 rests on.

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
import shutil
import subprocess
import sys

import h5py
import numpy as np

import beamio
import fieldio
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
  source_filter = F
  lambda0 = {lam}
  beam_init%n_particle = 2048
  beam_init%bunch_charge = {q}
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-07
  beam_init%a_norm_emit = 4e-9
  beam_init%b_norm_emit = 4e-9
  beamlet_size = 8
  seed_power = {power}
  seed_waist_size = {w0}
  grid_n_pts = 128
  grid_half_width = 2e-3
  ran_seed = 4242
  unaveraged_steps_per_period = {spp}
  unaveraged_ramp_periods = {ramp}
  write_initial = T
  write_diag = T
&end
"""

def probe_nml(wd, **kw):
    """probe with the steady-state charge derived (I = Q*c/spacing, spacing = lam), the
    unaveraged mode selected on the lattice and its two numbers stated by the run."""
    kw.setdefault("q", f"{3000 * float(kw['lam']) / 2.99792458e8:.12e}")
    kw["ramp"] = int(kw["ramp"])
    kw["lat"] = unavg_wrapper(wd, kw["lat"])
    return PROBE.format(**kw)

GAIN = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 5e3
  seed_waist_size = 30e-6
  grid_n_pts = 128
  grid_half_width = 2e-4
  ran_seed = 4242
  write_diag = T
&end
"""

TDID = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 8e-10
  n_wavelength = 1
  shot_noise = T
  ran_seed = 999
  write_diag = T
&end
"""

def unavg_wrapper(wd, base):
    """A wrapper lattice selecting the unaveraged mode. The mode is per element, so it is
    a lattice attribute; its two numbers are uniform, so the run states them."""
    name = f"w_{base.replace('.bmad','')}.bmad"
    (wd / name).write_text(
        f"call, file = {base}\n"
        f"wiggler::*[FEL_METHOD] = unaveraged\n")
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
    """The run must fail, with its own message."""
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


def dumps_identical(fa, fb):
    """Whether two openPMD dumps hold the same datasets to the bit, meta excluded."""
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names = []
        a.visititems(lambda n, o: names.append(n)
                     if isinstance(o, h5py.Dataset) and not n.startswith("meta/") else None)
        return all(n in b and np.array_equal(a[n][()], b[n][()]) for n in names)


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


Q_L = 0.12            # the quiver probe's undulator: 8 periods, 2-period ramps at each end
Q_NSTEP = 8           # records a period, so a frame series samples 8 phases of the quiver
Q_RAMP = 2
P0_MC = math.sqrt(GAMMA0**2 - 1)

QUIVER = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 128
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 0
  grid_n_pts = 32
  grid_half_width = 2e-4
  ran_seed = 4242
  unaveraged_steps_per_period = 40
  unaveraged_ramp_periods = 2
  dump_at_comb = {frames}
  dump_orbit = T
  write_initial = T
&end
"""


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

    # 1b. The ledger on a field that diffracts. The broad seed above diffracts little in a
    #     substep, and a source added after the diffraction carried a cross term the kick
    #     had not charged, first order in the substep's diffraction phase: 1.2 percent of
    #     the turnover at 20 substeps a period on a 100 um waist, halving with the
    #     substep. The source now lands on the record the kick read, so this closes at the
    #     same floor as the broad seed at both resolutions, and no trend with the substep.
    for spp in (20, 40):
        root = f"uv_ledgd{spp}"
        run(exe, wd, root, probe_nml(wd, lat="unavg_probe_helical_b.bmad", root=root,
            spp=spp, ramp=N_RAMP, lam=LAMBDA1, power=1e8, w0=1e-4))
        led = np.loadtxt(wd / f"{root}.ledger.txt")
        etot = led[:, 1] + led[:, 2]
        turnover = np.abs(np.diff(led[:, 2])).sum()
        check(f"energy ledger, 100 um seed, {spp} substeps a period: max|d(E_beam+U_field)| / sum|dU|",
              np.abs(etot - etot[0]).max() / max(turnover, 1e-300), 1e-4)

    # 1c. The order of the split, measured. The lattice, the beam and the transverse grid
    #     are held and the substep is refined, and the error of each resolution against the
    #     finest is read for the field, as a complex array, and for the particles' energies.
    #     The split is second order, and it measured first order for a diffracting seed:
    #     the record the step advances is the substep's midpoint field, and it was named
    #     for the plane half a substep behind it, so the particles worked against a field
    #     offset by a half-substep diffraction and the exit plane was handed the same
    #     offset. With the half-step carries at the segment's ends the
    #     observed order between successive resolutions must sit near two.
    print("--- the split's order, fixed grid, 100 um seed:")
    spps = (20, 40, 80, 160)
    fld, gam = {}, {}
    for spp in spps:
        root = f"uv_ord{spp}"
        run(exe, wd, root, probe_nml(wd, lat="unavg_probe_helical_b.bmad", root=root,
            spp=spp, ramp=N_RAMP, lam=LAMBDA1, power=1e8, w0=1e-4))
        fld[spp] = fieldio.read_field(wd / f"{root}-final.wf.h5", "x")["u"]
        gam[spp] = read_par(wd / f"{root}-final.beam.h5", LAMBDA1)["gamma"]
    ref_f, ref_g = fld[spps[-1]], gam[spps[-1]]
    e_f = {s: float(np.linalg.norm(fld[s] - ref_f) / np.linalg.norm(ref_f)) for s in spps[:-1]}
    e_g = {s: float(np.max(np.abs(gam[s] - ref_g)) * M_ELECTRON) for s in spps[:-1]}
    for s in spps[:-1]:
        print(f"      {s:3d} substeps/period: field {e_f[s]:.3e}, energy {e_g[s]:.3e} eV")
    p_f = math.log2(e_f[20] / e_f[40])
    p_g = math.log2(e_g[20] / e_g[40])
    check("order: field, observed between 20 and 40 substeps a period (2 expected, floor 1.7)",
          1.7 / p_f, 1.0, note=f"[p = {p_f:.2f}]")
    check("order: particle energy, observed between 20 and 40 (2 expected, floor 1.7)",
          1.7 / p_g, 1.0, note=f"[p = {p_g:.2f}]")

    # 2. Ballistic dark run: B does no work. Ramps hand the emittance back.
    run(exe, wd, "uv_dark", probe_nml(wd, lat="unavg_probe_planar_b.bmad", root="uv_dark",
        spp=20, ramp=N_RAMP, lam=LAMBDA1, power=0.0, w0=SEED_W0))
    d0 = read_par(wd / "uv_dark-initial.beam.h5", LAMBDA1)
    d1 = read_par(wd / "uv_dark-final.beam.h5", LAMBDA1)
    check("ballistic: max|dgamma| (B does no work)",
          float(np.abs(d1["gamma"] - d0["gamma"]).max() / GAMMA0), 1e-12)
    check("ballistic: |emit_x out/in - 1| (ramp handoff)",
          abs(emit(d1) / emit(d0) - 1), 1e-6)
    # The coherent-quiver instruments. A hard-edge entry (unaveraged_ramp_periods = -1, the
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

    # 5. Convergence at 10/20/30 steps/period (planar pair), which the floor rests on.
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
        lat=unavg_wrapper(wd, "aramis_1seg.bmad")))
    run(exe, wd, "uv_gain_avg", GAIN.format(root="uv_gain_avg", lat="aramis_1seg.bmad"))
    pu = np.loadtxt(wd / "uv_gain_unavg.diag.txt")[:, 2]
    pa = np.loadtxt(wd / "uv_gain_avg.diag.txt")[:, 2]
    lnr = abs(math.log(pu[-1] / pa[-1]))
    check(f"gain curve: |ln(P_unavg/P_avg)| at segment exit  [{pu[-1]:.4e} vs {pa[-1]:.4e}]",
          lnr, 0.2, note="(priced integrator-structure difference)")

    # 7. Mixed line (Stage A): averaged / unaveraged / averaged sandwich with pipe
    # interludes. Completing at all exercises the convention-flag asserts at real
    # internal boundaries. The ledger must be confined to and conserved over the
    # unaveraged segment. A wake on that segment must refuse, and the exit
    # power is priced against the all-averaged twin.
    (wd / "sandwich_avg.bmad").write_text(
        "call, file = unavg_sandwich.bmad\n"
        "UNDB[FEL_METHOD] = averaged\n")
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
    # A diagnostic must not steer the run it observes, and a frame written inside the
    # unaveraged segment is the one that could: the record there is the substep's midpoint
    # field and the frame carries the plane's, so the writer diffracts. It does so on a
    # copy. With frames every 0.1 m through the mixed line, six of them inside the
    # unaveraged segment, the final beam and field must be the frames-off run's to the bit.
    run(exe, wd, "uv_sandfr", GAIN.format(root="uv_sandfr", lat="unavg_sandwich.bmad").replace(
        "  write_diag = T\n", "  write_diag = T\n  global%dump_at_comb = T\n  global%comb_ds_save = 0.1\n"))
    inside = 0
    for fr in sorted(wd.glob("uv_sandfr-[0-9]*.wf.h5")):
        with h5py.File(fr) as h5:
            z = float(np.ravel(h5.attrs["sPosition"])[0])
        inside += int(0.81 + 1e-9 < z < 1.41 - 1e-9)
    same = all(dumps_identical(wd / f"uv_sandfr-final.{k}.h5", wd / f"uv_sand-final.{k}.h5")
               for k in ("beam", "wf"))
    check(f"sandwich: frames on against off, final beam and field identical to the bit "
          f"({inside} frames inside the unaveraged segment) (0 = yes)", 0.0 if same else 1.0, 0.5)
    check("sandwich: frames landed inside the unaveraged segment (5 expected)", abs(inside - 5), 0.5)

    ps = np.loadtxt(wd / "uv_sand.diag.txt")[:, 2]
    pa2 = np.loadtxt(wd / "uv_sand_avg.diag.txt")[:, 2]
    check(f"sandwich: |ln(P_mixed/P_averaged)| at exit  [{ps[-1]:.4e} vs {pa2[-1]:.4e}]",
          abs(math.log(ps[-1] / pa2[-1])), 5e-2, note="(one segment's ramp+mode price)")
    (wd / "sandwich_wake.bmad").write_text(SAND_WAKE_LAT)
    refused = run_expect_refusal(exe, wd, "uv_sandw",
        GAIN.format(root="uv_sandw", lat="sandwich_wake.bmad"),
        "ELEMENT SR WAKES ARE NOT SUPPORTED IN THE UNAVERAGED MODE")
    check("sandwich: wake on the unaveraged segment refused (1 = yes)",
          0.0 if refused else 1.0, 0.5)

    # 8. Thread invariance: the parallel slice loop must be invisible. A multi-slice
    # time-dependent unaveraged run (slippage active, so slices genuinely interleave
    # through the field) at 1 thread and at 8 threads must produce byte-identical
    # diagnostics AND ledger -- the same guarantee the averaged path carries
    # (per-slice private state, fixed-order energy reduction).
    wl = unavg_wrapper(wd, "aramis_1seg.bmad")
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

    # 8. The single-precision twin (doc/validation.md, and this mode's section of
    # fel_unaveraged_mod). Apple GPUs carry no FP64, so a device port of this mode needs
    # a single-precision form of the advance with its divergence from FP64 measured
    # first. The twin is that measurement, and the levels below are what it reads.
    #
    # Four rows are asserted. The transverse pair and the energy price the reformulated
    # charts. theta prices the phase, and it is the row the mode lives on: the slippage
    # rate is a difference of two numbers that are both one to within 1.3e-8, so the
    # naive single-precision form returns exactly zero and the phase is destroyed.
    # Replacing the identity with that form takes this row from 1.3e-5 rad to 1.7 rad,
    # which is the mutation record for it.
    #
    # The instrument may not steer what it observes, so the FP64 stream and the ledger
    # are required byte-identical with the twin on against off. That is the same demand
    # the averaged instrument carries and it is the one that would catch a twin writing
    # into the run's own state.

    print("--- the single-precision twin of the unaveraged advance:")
    run(exe, wd, "uv_tw", GAIN.format(root="uv_tw", lat=unavg_wrapper(wd, "aramis_1seg.bmad"))
        .replace("&end", '  global%fp32_check = "lockstep"\n&end'))
    rows = [l.split() for l in (wd / "uv_tw.fp32.txt").read_text().splitlines()
            if l and not l.startswith("#") and l.split()[0].isdigit()]
    if not rows:
        check("twin: the stream carries per-step rows", 1.0, 0.5, note="[none written]")
    else:
        cols = np.array([[float(v) for v in r] for r in rows])
        worst = cols.max(axis=0)
        print(f"      {len(rows)} steps instrumented, worst over the run:")
        for i, nm in ((1, "x"), (2, "px"), (3, "y"), (4, "py"), (5, "pz"), (6, "theta"),
                      (7, "phasor"), (8, "ledger"), (9, "field")):
            print(f"        {nm:<7} {worst[i]:.3e}")
        check("twin: worst transverse row", max(worst[1], worst[3]), 1e-5)
        check("twin: worst transverse momentum row", max(worst[2], worst[4]), 1e-4)
        check("twin: worst energy row", worst[5], 1e-4)
        check("twin: worst phase row [rad]", worst[6], 1e-3,
              note="(the slippage identity holds it here; the naive form reads 1.7)")
        check("twin: worst phasor row", worst[7], 1e-5)
        # The energy the twin's kicks moved against the FP64 step's, scaled by the
        # energy the step actually moved rather than by what it netted: over a record
        # step the gains and losses very nearly cancel, and against that net the ratio
        # reaches 3, which says something about the cancellation and nothing about
        # single precision.
        check("twin: worst ledger row, its energy against the FP64 step's", worst[8], 1e-4)
        # The field row is this mode's and not the averaged instrument's: sixty substeps
        # of transform pair, rounded propagator and source add stand behind one row where
        # the averaged twin's field crosses one. The row is the twin's record against the
        # FP64 field in the L2 norm, so the tolerance is the deeper accumulation's.
        check("twin: worst field row, after every substep of the record step", worst[9], 5e-5)

    run(exe, wd, "uv_twoff", GAIN.format(root="uv_twoff",
        lat=unavg_wrapper(wd, "aramis_1seg.bmad")))
    same_diag = (wd / "uv_tw.diag.txt").read_bytes() == (wd / "uv_twoff.diag.txt").read_bytes()
    same_led = (wd / "uv_tw.ledger.txt").read_bytes() == (wd / "uv_twoff.ledger.txt").read_bytes()
    check("twin: the FP64 stream is identical with it on against off (0 = yes)",
          0.0 if same_diag else 1.0, 0.5)
    check("twin: the ledger is identical with it on against off (0 = yes)",
          0.0 if same_led else 1.0, 0.5)

    # freerun compounds a single-precision state across steps and this twin carries none,
    # so it is refused rather than reported as something it is not.
    run_expect_refusal(exe, wd, "uv_twfree", GAIN.format(root="uv_twfree",
        lat=unavg_wrapper(wd, "aramis_1seg.bmad"))
        .replace("&end", '  global%fp32_check = "freerun"\n&end'),
        "DOES NOT COVER THE UNAVERAGED MODE")


    # 8. The quiver an exported frame carries. An openPMD momentum record is the
    # instantaneous kinetic momentum, and the averaged map stores the guiding centre, so
    # the writer keeps the two apart: the guiding centre goes to a frame's .gc.h5, and the
    # .beam.h5 inside an undulator is written only with global%dump_orbit, its records
    # the orbit fel_restore_quiver rebuilds (fel-physics.md sec-quiverdump). These runs
    # ask for it, and this is that reconstruction's validation. This mode is the reference: the same device tracked
    # here resolves the orbit itself. The runs are dark and share a beam, so the two differ
    # by the averaging alone, and the comb takes 8 frames a period so the series samples 8
    # phases of the quiver. The residual that indicts the conversion is the part that
    # oscillates with the undulator phase: the smooth drift of the two integrators is not
    # it. What the tolerance is for is a wrong constant, not the averaging order: a missing
    # sqrt(2) in the planar amplitude leaves 0.29 of the quiver behind and a quarter-period
    # phase origin 1.0 of it, where the worst measured here is 1e-4 of it.

    def quiver_lat(name, helical, tilt):
        """The device tracked averaged, and the wrapper that tracks the same one here."""
        peak = "" if helical else "sqrt(2) * "
        tstr = f", tilt = {tilt}" if tilt else ""
        (wd / f"{name}.bmad").write_text(f"""no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = {GAMMA0} * m_electron
beginning[beta_a] = 8.53711
beginning[alpha_a] = -0.703306
beginning[beta_b] = 17.3899
beginning[alpha_b] = 1.40348
QU: wiggler, l = {Q_L}, l_period = {LAMBDA_W}, &
    field_calc = {'helical_model' if helical else 'planar_model'}, &
    b_max = {peak}{AW} * (twopi / {LAMBDA_W}) * m_electron / c_light, &
    fel_method = averaged, ds_step = {LAMBDA_W / Q_NSTEP}{tstr}
QLINE: line = (QU)
use, QLINE
""")
        return f"{name}.bmad", unavg_wrapper(wd, f"{name}.bmad")

    def q_frames(root):
        """{s: the frame's records, particles in id order} over a run's series."""
        out = {}
        for f in sorted(wd.glob(f"{root}-0*.beam.h5")):
            with h5py.File(f) as h5:
                s_f = round(float(np.atleast_1d(h5.attrs["sPosition"])[0]), 9)
                phi0 = float(np.atleast_1d(h5.attrs["phi0"])[0])
                n0 = sorted(h5["data"].keys())[0]
                gp = h5[f"data/{n0}/particles"]
                ids = gp[sorted(gp.keys())[0]]["id"][...]
            par = beamio.read_slices(f, LAMBDA1, LAMBDA1)[0]
            o = np.argsort(ids)
            out[s_f] = {k: par[k][o] for k in ("x", "y", "px", "py", "gamma", "theta", "weight")}
            out[s_f]["phi0"] = phi0
        return out

    def oscillating(fa, fu, key, ss):
        """The residual's part that oscillates with the undulator phase, over one period at
        a time, and the quiver's own amplitude there. Both in the record's own units."""
        d = np.array([fa[s][key] - fu[s][key] for s in ss])
        q = np.array([fu[s][key] for s in ss])
        worst = amp = 0.0
        for i in range(0, len(ss) - Q_NSTEP + 1, Q_NSTEP):
            w = d[i:i + Q_NSTEP]
            worst = max(worst, float(np.max(np.abs(w - w.mean(axis=0)))))
            qq = q[i:i + Q_NSTEP]
            amp = max(amp, float(np.max(qq.max(axis=0) - qq.min(axis=0)) / 2))
        return worst, amp

    l_ramp = Q_RAMP * LAMBDA_W
    for name, helical, tilt in (("planar", False, 0.0), ("helical", True, 0.0),
                                ("planar tilted", False, 0.4)):
        tag = "uv_q" + name.split()[0][0] + ("t" if tilt else "")
        base, wrap = quiver_lat(f"q_{tag}", helical, tilt)
        run(exe, wd, tag + "a", QUIVER.format(lat=base, root=tag + "a", frames="T"))
        run(exe, wd, tag + "u", QUIVER.format(lat=wrap, root=tag + "u", frames="T"))
        check(f"quiver ({name}): the two runs start from the same beam (0 = yes)",
              0.0 if dumps_identical(wd / f"{tag}a-initial.beam.h5",
                                     wd / f"{tag}u-initial.beam.h5") else 1.0, 0.5)

        fa, fu = q_frames(tag + "a"), q_frames(tag + "u")
        shared = sorted(set(fa) & set(fu))
        body = [s for s in shared if l_ramp + 1e-9 < s < Q_L - l_ramp - 1e-9]
        ramp = [s for s in shared if 1e-9 < s < l_ramp - 1e-9]
        worst_rel, note = 0.0, []
        for where, ss, keys in (("body", body, ("x", "px", "theta")), ("ramp", ramp, ("x", "px"))):
            for key in keys:
                res, amp = oscillating(fa, fu, key, ss)
                worst_rel = max(worst_rel, res / amp)
                note.append(f"{where} {key} {res:.1e} of {amp:.1e}")
        check(f"quiver ({name}): the restored orbit is this mode's, to a fraction of the quiver",
              worst_rel, 1e-3, note="[" + ", ".join(note) + " (m, u, rad)]")

        # Inside a ramp the two devices differ in phase as well: this mode jumps the ramp's
        # phase at the segment's first step and the deficit then accrues through the ramp,
        # so the difference opens at the entrance and closes by the ramp's end. It is not
        # the orbit's, and the conversion leaves it alone (fel-physics.md sec-quiverdump).
        d_ramp = [float(np.max(np.abs(fa[s_f]["theta"] - fu[s_f]["theta"]))) for s_f in ramp]
        d_body = float(np.max(np.abs(fa[body[0]]["theta"] - fu[body[0]]["theta"])))
        check(f"quiver ({name}): the ramp's phase transient opens at the entrance and closes by its end [rad]",
              d_body, 1e-2, note=f"[{max(d_ramp):.2f} rad at the entrance, {d_ramp[-1]:.4f} at the last ramp "
                                 f"frame, {d_body:.4f} through the body, on a quiver of "
                                 f"{oscillating(fa, fu, 'theta', body)[1]:.2f} rad]")

        # The inverse, against the run's own record of the state it tracked: the stats
        # centroid is taken from the live beam, before the writer's copy is converted.
        with h5py.File(wd / f"{tag}a.stats.h5") as h5:
            s_rec, cen = h5["coords/s"][...], h5["beam/slice/centroid"][...]
        back = {k: 0.0 for k in ("x", "px", "y", "py", "z", "pz")}
        for f in sorted(wd.glob(f"{tag}a-0*.beam.h5")):
            with h5py.File(f) as h5:
                s_f = round(float(np.atleast_1d(h5.attrs["sPosition"])[0]), 9)
                phi0 = float(np.atleast_1d(h5.attrs["phi0"])[0])
            if s_f not in body and s_f not in ramp:
                continue
            gc = beamio.unquiver(f, beamio.read_slices(f, LAMBDA1, LAMBDA1), LAMBDA1)[0]
            wgt, gam = gc["weight"], gc["gamma"]
            beta = np.sqrt(1 - 1 / gam**2)
            p_mc = np.sqrt(gam**2 - 1)
            got = {"x": gc["x"], "y": gc["y"], "px": gc["px"] / P0_MC, "py": gc["py"] / P0_MC,
                   "z": beta * (gc["theta"] - phi0) * LAMBDA1 / (2 * math.pi),
                   "pz": (p_mc - P0_MC) / P0_MC}
            irec = int(np.argmin(np.abs(s_rec - s_f)))
            for i, k in enumerate(("x", "px", "y", "py", "z", "pz")):
                back[k] = max(back[k], abs(float(np.average(got[k], weights=wgt)) - float(cen[irec, 0, i])))
        _, amp_x = oscillating(fa, fu, "x", body)
        # beamio.unquiver is the page's inverse, reading the device's numbers from the
        # frame's own attributes. A planar roll-off depends on the wiggle frame's y, which
        # the quiver does not move, and a helical one moves with r, which the inverse's
        # second pass answers, so all three land on the rounding floor.
        check(f"quiver ({name}): the page's inverse returns the guiding centre the run tracked [m]",
              back["x"], 1e-16,
              note=f"[x {back['x']:.1e} m on a quiver of {amp_x:.1e} m, px {back['px']:.1e}, "
                   f"z {back['z']:.1e} m, pz {back['pz']:.1e}]")

        if not helical and not tilt:
            # The same frame's .gc.h5 is the guiding centre with no inverse at all, and the
            # inverse of the exported .beam.h5 lands on it. A frame of the 1.0 vintage may
            # hold either representation, so the inverse refuses to guess at one.
            f_o = sorted(wd.glob(f"{tag}a-0*.beam.h5"))[len(body) // 2 + Q_NSTEP * 2]
            f_g = pathlib.Path(str(f_o).replace(".beam.h5", ".gc.h5"))
            gc = beamio.read_slices(f_g, LAMBDA1, LAMBDA1)[0]
            inv = beamio.unquiver(f_o, beamio.read_slices(f_o, LAMBDA1, LAMBDA1), LAMBDA1)[0]
            check("quiver: the export's inverse lands on the frame's own guiding-centre file [m]",
                  float(np.max(np.abs(inv["x"] - gc["x"]))), 1e-16,
                  note=f"[x {np.max(np.abs(inv['x'] - gc['x'])):.1e} m, px {np.max(np.abs(inv['px'] - gc['px'])):.1e}, "
                       f"theta {np.max(np.abs(inv['theta'] - gc['theta'])):.1e} rad; unquiver leaves the .gc.h5 "
                       f"itself alone: {beamio.unquiver(f_g, [gc], LAMBDA1)[0] is gc}]")
            old = wd / "uv_q_old.beam.h5"
            shutil.copy(f_o, old)
            with h5py.File(old, "r+") as h5:
                h5.attrs["frameFormat"] = np.bytes_("lucifer-frames 1.0")
            try:
                beamio.unquiver(old, beamio.read_slices(old, LAMBDA1, LAMBDA1), LAMBDA1)
                guessed = 1.0
            except ValueError:
                guessed = 0.0
            stated = beamio.unquiver(old, beamio.read_slices(old, LAMBDA1, LAMBDA1), LAMBDA1,
                                     representation="orbit")[0]
            check("quiver: an interior frame of the 1.0 vintage is not guessed at, and is inverted once its "
                  "representation is stated (0 = yes)",
                  guessed + float(np.max(np.abs(stated["x"] - inv["x"])) > 0), 0.5)

            # The conversion is the writer's own copy: a run with the series on holds the
            # state a run without it holds, to the bit.
            run(exe, wd, tag + "n", QUIVER.format(lat=base, root=tag + "n", frames="F"))
            check("quiver: the frame series does not move the run it observes (0 = yes)",
                  0.0 if dumps_identical(wd / f"{tag}a-final.beam.h5",
                                         wd / f"{tag}n-final.beam.h5") else 1.0, 0.5)

    if FAILED:
        print("unaveraged checks: FAIL")
        sys.exit(1)
    print("unaveraged checks: PASS")


if __name__ == "__main__":
    main()
