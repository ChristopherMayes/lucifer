#!/usr/bin/env python3
"""
Checks of the coherent source (fel-physics.md, sec-coherent-source), implemented from
Tanaka, PRAB 27, 030703 (2024): SIMPLEX's coherent source driving Lucifer's field.

  1. Every other check in the harness leaves source_model at its default, the deposit,
     so the benchmark's tiers hold that path. Here a coherent run on the reference
     configuration runs clean.
  2. Limit identity. At large M on the seeded 6 GeV segment, the coherent and deposit
     |ln P| curves differ by the coherent Gaussian's own error in the lethargy regime,
     the regime where the low-M artifact is largest. The bias is 2.0e-2, checked at 5e-2.
  3. Variance reduction. At M = 128 the plain deposit fakes gain, ln P of +0.41 on a
     curve that absorbs, while the coherent source stays within the model bias of the
     large-M reference, at 0.05. The artifact must exceed 0.2 and the coherent error
     stay under 0.1. validation.md (val-the-coherent-source-simplex) records both.
  4. Moments extension. An offset and tilted Gaussian beam, written here as openPMD and
     imported, runs coherent and tracks its own deposit twin at the tolerance of 2. A
     double-horn beam is refused by the significance guard, the excess kurtosis
     against sqrt(24/m_ind).
  5. A dark start is refused. The coherent model drops the spontaneous, spatially
     incoherent emission that dominates SASE startup, and its startup comes out about
     175 times low. B(s) carries the physical Fawley noise, but that noise is not the
     dominant seed at startup.
  6. A coherent run is bit-identical at 1 and 8 threads.

Run by the benchmark harness. Exits nonzero on failure. No Genesis run is involved.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from checkrec import record
from pool import run_all
from read_stats import read_stats

FAILED = False

TOL_LIMIT = 5e-2      # The coherent Gaussian's bias in lethargy is 2.0e-2.
TOL_CLAIM_COH = 0.1   # Coherent at M = 128 against the M = 8192 reference is 0.05.
MIN_CLAIM_DEP = 0.2   # The deposit's fake gain at M = 128 is 0.41 and must exceed this.

LAT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
UNDP: wiggler, l = 3.96, l_period = 0.015, field_calc = planar_model, &
      b_max = sqrt(2) * 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      fel_method = averaged, ds_step = 0.045
SEG: line = (UNDP)
use, SEG
"""

NML = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%source_filter = F
  global%interlude_model = "genesis"
  global%transport_model = "genesis"
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
  wavefront_init%grid_n_pts = 128
  wavefront_init%grid_half_width = 2e-4
{wextra}/
"""

COHERENT = '  global%source_model = "coherent"\n'

SASE_B = """  beamlet_size = 8
  shot_noise = T
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.6e-9
  beam_init%grid(3)%x_max = 1.6e-9
"""
SASE_W = """  slicing%window_length = 3.2e-9
  slicing%n_wavelength = 1
"""


def check(name, value, tol, note="", low=None):
    global FAILED
    ok = value <= tol if low is None else (low <= value)
    lim = f"(check {tol:.0e})" if low is None else f"(must exceed {low:.0e})"
    print(f"--- {name}: {value:.3e} {lim} {note} {'ok' if ok else '** FAIL **'}")
    if low is None:
        record(name, value, tol, "max", ok)
    else:
        record(name, value, low, "min", ok)
    if not ok:
        FAILED = True


def run(exe, wd, root, m, coherent, seed=777, bextra="", wextra="", pextra="", threads="8",
        sase=False, expect_fail=None, sig_z=None, lat="coh.bmad"):
    extra = (COHERENT if coherent else "") + (SASE_W if sase else "") + pextra
    text = NML.format(root=root, m=m, seed=seed, extra=extra, lat=lat,
                      bextra=(SASE_B if sase else "") + bextra,
                      wextra=wextra)
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

    # An offset Gaussian bunch, tilted in x-y, is written as openPMD and imported by both
    # source models, and the coherent run must track its own deposit twin. The import
    # runs where the resampler is validated, a time-dependent window of many slices
    # with resample%n_slice = 0 so the count is derived. A single slice samples the
    # beam degenerately in x and y. beamlet_size = 4 keeps the guard's estimate of
    # independent samples honest, so the double horn trips it.

    # A short undulator, 0.9 m, and sample = 20 keep the seed in the window. At
    # sample = 1 the slippage, 3 wavelengths a step here, flushes the seed out within a
    # few steps and the comparison becomes the dark-start case: the deposit's inflated
    # spontaneous emission against none from the coherent model, an apparent |ln P|
    # disagreement of 5.9.

    (wd / "cohshort.bmad").write_text(LAT.replace("l = 3.96,", "l = 0.90,"))
    IMP = ('  beam_init%position_file = "beam0.h5"\n  beam_init%n_particle = 2048\n  beamlet_size = 4\n'
           '  shot_noise = T\n  resample%n_slice = 0\n  resample%slice_width = 0.05\n')
    WIN = '  slicing%n_wavelength = 20\n'
    run(exe, wd, "seedbeam", 200000, False, pextra='  global%load_only = T\n', sig_z='4e-9',
        bextra='  resample%use_beam_init = T\n  write_openpmd_file = "beam0.h5"\n'
               '  resample%n_particle_per_slice = 2048\n  resample%beamlet_size = 4\n')
    with h5py.File(wd / "beam0.h5", "r+") as f:
        # openPMD-beamphysics as Bmad writes it: /data/%T/particles/electron/...
        for it in f["data"]:
            pp = f[f"data/{it}/particles/electron"]
            x = pp["position/x"][()]
            y = pp["position/y"][()]
            pp["position/x"][...] = x + 2e-5 + 0.3 * y      # offset + tilt
            break
    for root, coh in (("mdep", False), ("mcoh", True)):
        run(exe, wd, root, 2048, coh, sig_z='4e-9', bextra=IMP, pextra=WIN, lat="cohshort.bmad")
    _, md = curve(wd, "mdep")
    _, mc = curve(wd, "mcoh")
    # The twins differ by 8.9e-2: the steady-state model bias of 1.9e-2 plus the
    # emission over the 3 slices this window still slips. The deposit refills a
    # slipped-in slice with its inflated spontaneous emission and the coherent source
    # does not. Without the centering of the moments extension the Gaussian would sit
    # 0.87 sigma off the beam and the error would be of order 1, which this check pins.
    check("offset+tilted beam: coherent tracks its deposit twin",
          float(np.max(np.abs(np.log(mc / md)))), 2e-1)

    with h5py.File(wd / "beam0.h5", "r+") as f:
        for it in f["data"]:
            pp = f[f"data/{it}/particles/electron"]
            x = pp["position/x"][()]
            pp["position/x"][...] = x + 6e-5 * np.sign(x - np.median(x))  # double horn
            break
    ok = run(exe, wd, "mhorn", 2048, True, sig_z='4e-9', bextra=IMP, pextra=WIN, lat="cohshort.bmad",
             expect_fail="GAUSSIAN ENOUGH")
    print(f"--- double-horn beam refused (GAUSSIAN ENOUGH): {'ok' if ok else '** FAIL **'}")
    record("double-horn beam refused (GAUSSIAN ENOUGH)", None, None, "bool", ok)
    FAILED = FAILED or not ok

    # ------------------------------------------------------------------
    print("== dark start refused (the measured startup deficit) ==")

    # On this configuration a coherent dark start comes out about 175 times low, 7.5 W
    # against 1311 W from the deposit at the same seeds. The coherent model drops the
    # spontaneous, spatially incoherent emission that dominates SASE startup, so a dark
    # start is refused rather than reported as a startup the model cannot represent.

    ok = run(exe, wd, "sdark", 2048, True, sase=True, expect_fail="DARK START")
    print(f"--- coherent + dark start refused: {'ok' if ok else '** FAIL **'}")
    record("coherent + dark start refused", None, None, "bool", ok)
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
