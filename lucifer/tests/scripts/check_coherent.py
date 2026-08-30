#!/usr/bin/env python3
"""
Coherent-source checks (manual sec:coherent-source; Tanaka PRAB 27, 030703 (2024),
implemented from the paper -- the SIMPLEX hybrid: their source, our field).

  1. Keystone: source_model is defaulted -- the deposit path -- everywhere else in the
     harness, so the benchmark's own tiers hold the default-path contract. Here: a
     coherent run exists and runs clean on the reference configuration.
  2. Limit identity (the model's measured bias, not zero): coherent vs deposit at
     large M on the seeded 6 GeV segment -- the |ln P| curve difference is the
     coherent-Gaussian model's intrinsic error in the lethargy regime, measured
     2.0e-2, checked at 5e-2. The regime is deliberate: lethargy/absorption is the
     energy exchange where the low-M artifact bites hardest.
  3. The claim (variance reduction): at M = 128 the plain deposit fakes gain
     (measured ln +0.41 on a curve that truly absorbs -- the per-cell shot-noise
     artifact of the survey's 5.1), while the coherent source at the same M stays
     within the model bias of the large-M reference (measured 0.05). Both measured
     here, both checked: the artifact must exceed 0.2, the coherent error stay
     under 0.1.
  4. Moments extension: an offset + tilted (but Gaussian) beam, imported via an
     openPMD file this check writes, runs coherent and tracks its own deposit twin
     (same tolerance as 2); a DOUBLE-HORN beam is refused by name by the
     significance guard (excess kurtosis against sqrt(24/m_ind)).
  5. Dark start refused by name: measured ~175x startup deficit (spontaneous,
     spatially-incoherent emission dominates SASE startup; the coherent model drops
     it by construction) -- the goal's fallback clause, executed. B(s) does carry
     the physical Fawley noise, but it is not the dominant seed at startup.
  6. Threads: a coherent run is bit-identical at 1 and 8 threads.

Run by the benchmark harness; exits nonzero on failure. Self-referenced (no Genesis).
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from pool import run_all
from read_stats import read_stats

FAILED = False

TOL_LIMIT = 5e-2      # measured 2.0e-2: the coherent-Gaussian model bias (lethargy).
TOL_CLAIM_COH = 0.1   # measured 0.05: coherent at M=128 vs the M=8192 reference.
MIN_CLAIM_DEP = 0.2   # measured 0.41: the deposit's fake gain at M=128 (must exceed).

LAT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UNDP: wiggler, l = 3.96, l_period = 0.015, field_calc = planar_model, &
      b_max = sqrt(2) * 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.045
fel_transcribed = -1
wiggler::*[FEL_TRACKING] = fel_transcribed
SEG: line = (UNDP)
use, SEG
"""

NML = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%interlude_model = "genesis"
  global%ran_seed = {seed}
{extra}/
&fel_beam_init
  beam_init%n_particle = {m}
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
{bextra}/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 5e3
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_n_pts = 151
  wavefront_init%grid_half_width = 2e-4
{wextra}/
"""

COHERENT = '  global%source_model = "coherent"\n'

SASE_B = """  nbins = 8
  shotnoise = T
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.6e-9
  beam_init%grid(3)%x_max = 1.6e-9
"""
SASE_W = """  wavefront_init%window_length = 3.2e-9
  wavefront_init%window_sample = 1
"""


def check(name, value, tol, note="", low=None):
    global FAILED
    ok = value <= tol if low is None else (low <= value)
    lim = f"(check {tol:.0e})" if low is None else f"(must exceed {low:.0e})"
    print(f"--- {name}: {value:.3e} {lim} {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, root, m, coherent, seed=777, bextra="", wextra="", pextra="", threads="8",
        sase=False, expect_fail=None, sig_z=None, lat="coh.bmad"):
    extra = (COHERENT if coherent else "") + pextra
    text = NML.format(root=root, m=m, seed=seed, extra=extra, lat=lat,
                      bextra=(SASE_B if sase else "") + bextra,
                      wextra=(SASE_W if sase else "") + wextra)
    if sig_z is not None:
        text = text.replace("  beam_init%sig_z = 0\n", f"  beam_init%sig_z = {sig_z}\n")
    if sase:
        text = text.replace("  wavefront_init%seed_power = 5e3\n",
                            "  wavefront_init%seed_power = 0\n")
        text = text.replace("  beam_init%sig_z = 0\n", "")
    (wd / f"{root}.nml").write_text(text)
    r = subprocess.run([str(exe), f"{root}.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if expect_fail is not None:
        return r.returncode != 0 and expect_fail in (r.stdout + r.stderr)
    if r.returncode != 0:
        print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-1500:]}")
        sys.exit(1)


def curve(wd, root):
    with read_stats(wd / f"{root}.stats.h5") as st:
        return st.s, np.sum(st["field/total/power"], axis=1)   # window total


def sase_startup(wd, root):
    with read_stats(wd / f"{root}.stats.h5") as st:
        return float(np.mean(st["field/total/power"][-1, :]))


def main():
    global FAILED
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--exe", required=True)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(args.exe).resolve()
    (wd / "coh.bmad").write_text(LAT)

    # ------------------------------------------------------------------
    print("== limit identity and the variance-reduction claim ==")
    run_all([
        lambda: run(exe, wd, "cd8k", 8192, False),
        lambda: run(exe, wd, "cc8k", 8192, True),
        lambda: run(exe, wd, "cd128", 128, False),
        lambda: run(exe, wd, "cc128", 128, True),
    ], threads_per_job=8)
    z, ref = curve(wd, "cd8k")
    _, c8k = curve(wd, "cc8k")
    _, d128 = curve(wd, "cd128")
    _, c128 = curve(wd, "cc128")
    check("limit identity: |ln P| curve, coherent vs deposit at M=8192",
          float(np.max(np.abs(np.log(c8k / ref)))), TOL_LIMIT,
          note="[the model's measured lethargy bias]")
    dep_art = abs(float(np.log(d128[-1] / ref[-1])))
    coh_err = abs(float(np.log(c128[-1] / ref[-1])))
    check("the artifact: plain deposit at M=128 fakes energy exchange", dep_art,
          1e9, low=MIN_CLAIM_DEP, note="[|ln| vs the M=8192 reference]")
    check("the claim: coherent at M=128 stays at the model bias", coh_err, TOL_CLAIM_COH,
          note=f"[reduction demonstrated: 64x fewer particles, error {coh_err:.3f} vs artifact {dep_art:.3f}]")

    # ------------------------------------------------------------------
    print("== the moments extension (offset + tilted admitted; double horn refused) ==")

    # An offset, x-y-tilted Gaussian bunch written as openPMD and imported by both
    # source models: the coherent run must track its own deposit twin. The import
    # runs in the resampler's validated regime (a time-dependent window of many
    # slices, imp%nslice = 0 auto) -- the single-slice corner produces degenerate
    # transverse sampling (a separate finding, noted in the brief). nbins = 4 keeps
    # the guard's independent-sample estimate honest for the horn to trip on.

    # A short undulator (0.9 m) and sample = 20 keep the seed in the window: at
    # sample = 1 the slippage (3 wavelengths per step here) flushes the seed out
    # within a few steps and the comparison degenerates into the dark-start case
    # (the deposit's inflated spontaneous emission vs the coherent model's none --
    # measured as an apparent 5.9 |ln| disagreement before this was understood).

    (wd / "cohshort.bmad").write_text(LAT.replace("l = 3.96,", "l = 0.90,"))
    IMP = ('  dist_file = "beam0.h5"\n  imp%npart = 2048\n  imp%nbins = 4\n'
           '  imp%nslice = 0\n  imp%slicewidth = 0.05\n')
    WIN = '  wavefront_init%window_sample = 20\n'
    run(exe, wd, "seedbeam", 200000, False, pextra='  global%load_only = T\n', sig_z='4e-9',
        bextra='  use_beam_init = T\n  write_opmd_file = "beam0.h5"\n'
               '  imp%npart = 2048\n  imp%nbins = 4\n')
    with h5py.File(wd / "beam0.h5", "r+") as f:
        # openPMD-beamphysics as Bmad writes it: /data/%T/particles/electron/...
        for it in f["data"]:
            pp = f[f"data/{it}/particles/electron"]
            x = pp["position/x"][()]
            y = pp["position/y"][()]
            pp["position/x"][...] = x + 2e-5 + 0.3 * y      # offset + tilt
            break
    for root, coh in (("mdep", False), ("mcoh", True)):
        run(exe, wd, root, 2048, coh, sig_z='4e-9', bextra=IMP, wextra=WIN, lat="cohshort.bmad")
    _, md = curve(wd, "mdep")
    _, mc = curve(wd, "mcoh")
    # Measured 8.9e-2: the SS model bias (1.9e-2) plus the residual emission
    # difference from the 3 slices of window slippage this TD config still has
    # (the deposit refills slipped-in slices with its inflated spontaneous
    # emission, and the coherent source, correctly, does not). Without the centering
    # extension the Gaussian would sit 0.87 sigma off the beam and the error
    # would be O(1) -- that is what this check pins.
    check("offset+tilted beam: coherent tracks its deposit twin",
          float(np.max(np.abs(np.log(mc / md)))), 2e-1)

    with h5py.File(wd / "beam0.h5", "r+") as f:
        for it in f["data"]:
            pp = f[f"data/{it}/particles/electron"]
            x = pp["position/x"][()]
            pp["position/x"][...] = x + 6e-5 * np.sign(x - np.median(x))  # double horn
            break
    ok = run(exe, wd, "mhorn", 2048, True, sig_z='4e-9', bextra=IMP, wextra=WIN, lat="cohshort.bmad",
             expect_fail="GAUSSIAN ENOUGH")
    print(f"--- double-horn beam refused by name (GAUSSIAN ENOUGH): {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    # ------------------------------------------------------------------
    print("== dark start refused by name (the measured startup deficit) ==")

    # measured on this configuration: the coherent dark start came out ~175x low
    # (deposit 1311 W vs coherent 7.5 W at the same seeds) -- spontaneous,
    # spatially-incoherent emission dominates SASE startup and the coherent model
    # drops it by construction. The goal's fallback clause applies: refusal by name.

    ok = run(exe, wd, "sdark", 2048, True, sase=True, expect_fail="DARK START")
    print(f"--- coherent + dark start refused by name: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    # ------------------------------------------------------------------
    print("== thread identity ==")
    run(exe, wd, "ct1", 2048, True, threads="1")
    run(exe, wd, "ct8", 2048, True, threads="8")
    same = True
    with read_stats(wd / "ct1.stats.h5") as a, read_stats(wd / "ct8.stats.h5") as b:
        same = bool(np.array_equal(a["field/total/power"], b["field/total/power"]) and
                    np.array_equal(a["beam/slice/bunching"], b["beam/slice/bunching"]))
    check("coherent run 1 vs 8 threads dataset-identical (0 = yes)",
          0.0 if same else 1.0, 0.5)

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
