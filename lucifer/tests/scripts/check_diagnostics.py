#!/usr/bin/env python3
"""
Diagnostic-output checks (manual sec:stats): the stats file, dump-at elements, and the
escaped-field bank, held by cross-identities rather than reference files --

  1. bunch_params reconstruction: bunch_params_from_stats.py applied to the per-record
     sufficient statistics must reproduce the calc_bunch_params values the tracker
     stored at element ends (proves the per-record datasets ARE sufficient).
  2. banked energy == the ledger's U_escaped column: same slices, independent
     bookkeeping paths (drain-time wavefront_params vs the zero-fill energy sum).
  3. analytic vs numerical propagation: each banked slice's rms size at the exit from
     the free-space ABCD map on its bank-time moment matrix must match the rms of the
     FFT-propagated slice in the pulse file (two independent routes).
  4. pulse pooling: the pooled whole-pulse sigma from per-slice params must match the
     directly computed moments of the concatenated pulse file.
  5. thread invariance: every dataset of stats.h5, escaped and pulse files identical
     at 1 vs 8 threads (dataset-level: HDF5 headers embed creation times).
  6. refusal: a dump_beam_at entry matching no element is refused by name.

Run by the benchmark harness; exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bunch_params_from_stats import bunch_params_at, pool_wavefront  # noqa: E402

FAILED = False

NML = """&fel_track_params
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
  keep_escaped_field = T
  dump_beam_at = "UND"
  dump_field_at = "UND"
{extra}&end
"""

WRAP = """call, file = aramis_1seg.bmad
fel_unaveraged = 1
wiggler::*[FEL_TRACKING] = fel_unaveraged
"""


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".nml")).write_text(text)
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-1000:]}")
        sys.exit(1)


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
    (wd / "aramis_1seg.bmad").write_bytes((latdir / "aramis_1seg.bmad").read_bytes())
    (wd / "dg_wrap.bmad").write_text(WRAP)

    run(exe, wd, "dg", NML.format(lat="dg_wrap.bmad", root="dg", extra=""))

    with h5py.File(wd / "dg.stats.h5") as h5:
        nslice = h5["beam/centroid"].shape[1]
        ee = h5["element_end"]

        # 1. Reconstruction at the element end (= the last record).
        worst = 0.0
        for isl in range(nslice):
            bp = bunch_params_at(h5, -1, isl)
            for m in ("x", "y", "z"):
                for pn in ("beta", "alpha", "emit", "norm_emit", "sigma", "sigma_p", "eta", "etap"):
                    stored = ee[f"slice/{m}/{pn}"][0, isl]
                    worst = max(worst, abs(bp[m][pn] - stored) / max(abs(stored), 1e-30))
            st = sorted(ee[f"slice/{m}/emit"][0, isl] for m in ("a", "b", "c"))
            mine = sorted(bp[m]["emit"] for m in ("a", "b", "c"))
            for s0, m0 in zip(st, mine):
                worst = max(worst, abs(m0 - s0) / max(abs(s0), 1e-30))
        # Measured 4.6e-8 (numpy eig vs mat_eigen); the projected planes agree exactly.
        check("stats: bunch_params reconstruction vs stored calc_bunch_params", worst, 1e-6)

        f_cen = h5["field/centroid"][-1]
        f_sig = h5["field/sigma"][-1]
        f_en = h5["field/energy"][-1]

    # 2. Banked energy == ledger U_escaped.
    led = np.loadtxt(wd / "dg.ledger.txt")
    with h5py.File(wd / "dg-escaped.fld.h5") as h5:
        nb = int(h5["slicecount"][0])
        pms = np.array([h5[f"slice{k:06d}/wavefront_params"][:] for k in range(1, nb + 1)])
        zt = np.array([h5[f"slice{k:06d}/z_transmit"][0] for k in range(1, nb + 1)])
    e_banked = pms[:, 20].sum()
    check("bank: banked slice energies vs ledger U_escaped, rel",
          abs(e_banked / led[-1, 4] - 1), 1e-12)

    # 3. Analytic (moment-map) vs numerical (FFT) propagation to the exit plane.
    with h5py.File(wd / "dg-pulse.fld.h5") as h5:
        n_all = int(h5["slicecount"][0])
        nx = int(h5["gridpoints"][0])
        dx = float(h5["gridsize"][0])
        nlive = n_all - nb
        ax = (np.arange(nx) - 0.5 * (nx - 1)) * dx
        sxx_num, syy_num, w_num = [], [], []
        cx_num, cy_num = [], []
        for j in range(nlive + 1, n_all + 1):
            g = h5[f"slice{j:06d}"]
            fr = g["field-real"][:].reshape(nx, nx)   # x fastest: columns = x.
            fi = g["field-imag"][:].reshape(nx, nx)
            I = fr**2 + fi**2
            w = I.sum()
            w_num.append(w)
            if w == 0:
                sxx_num.append(0); syy_num.append(0); cx_num.append(0); cy_num.append(0)
                continue
            cx = (I * ax[None, :]).sum() / w
            cy = (I * ax[:, None]).sum() / w
            cx_num.append(cx); cy_num.append(cy)
            sxx_num.append((I * (ax[None, :] - cx) ** 2).sum() / w)
            syy_num.append((I * (ax[:, None] - cy) ** 2).sum() / w)
    # Pulse index nlive+j holds banked slice nb-j+1: reverse into transmission order.
    sxx_num = np.array(sxx_num)[::-1]
    syy_num = np.array(syy_num)[::-1]
    w_num = np.array(w_num)[::-1]

    z_end = led[-1, 0]  # Ledger rows end at the segment exit for this one-element line.
    dz = z_end - zt
    # Flattened 4x4 column-major: (i,j) -> 4*(j-1)+(i-1).
    s11, s21, s22 = pms[:, 4 + 0], pms[:, 4 + 1], pms[:, 4 + 5]
    s33, s43, s44 = pms[:, 4 + 10], pms[:, 4 + 11], pms[:, 4 + 15]
    sxx_ana = s11 + 2 * dz * s21 + dz**2 * s22
    syy_ana = s33 + 2 * dz * s43 + dz**2 * s44
    m = w_num > 1e-6 * w_num.max()
    rel = np.abs(np.sqrt(sxx_num[m] / sxx_ana[m]) - 1)
    rel = np.maximum(rel, np.abs(np.sqrt(syy_num[m] / syy_ana[m]) - 1))
    # Measured max 8.4e-3 (grid-edge wrap of the FFT route on divergent noise slices).
    check("bank: analytic vs FFT-propagated rms at exit, max rel", rel.max(), 2e-2)

    # 4. Pulse pooling: pooled per-slice params vs directly computed pulse moments.
    #    Banked side: params propagated analytically to the exit; live side: the last
    #    stats record (angles valid there). Compare sigma_xx of the pooled result
    #    against pooling the pulse file's numerically computed per-slice moments.
    cen_b = np.stack([pms[:, 0] + dz * pms[:, 1], pms[:, 1],
                      pms[:, 2] + dz * pms[:, 3], pms[:, 3]], axis=1)
    sig_b = np.zeros((nb, 4, 4))
    sig_b[:, 0, 0] = sxx_ana;  sig_b[:, 2, 2] = syy_ana
    en_b = pms[:, 20]
    cen_l = f_cen;  sig_l = f_sig.reshape(-1, 4, 4);  en_l = f_en
    _, pooled_ana, _ = pool_wavefront(np.vstack([cen_l, cen_b]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_b.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_b]))
    # Numerical pooling over the SAME banked slices from the pulse file:
    cen_n = np.stack([np.array(cx_num)[::-1], np.zeros(nb), np.array(cy_num)[::-1], np.zeros(nb)], axis=1)
    sig_n = np.zeros((nb, 4, 4))
    sig_n[:, 0, 0] = sxx_num;  sig_n[:, 2, 2] = syy_num
    en_n = w_num * (en_b.sum() / max(w_num.sum(), 1e-300))   # Same total, per-slice proportional.
    _, pooled_num, _ = pool_wavefront(np.vstack([cen_l, cen_n]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_n.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_n]))
    check("pulse: pooled sigma_xx analytic vs numerical routes, rel",
          abs(pooled_ana[0, 0] / pooled_num[0, 0] - 1), 2e-2)

    # 5. Thread invariance: every dataset of all three files identical at 1 vs 8
    #    threads. Dataset-level, not raw bytes: HDF5 object headers embed creation
    #    times, so two files with bit-identical DATA differ as byte streams -- the
    #    data is the invariance claim, the container metadata is not.
    run(exe, wd, "dg1", NML.format(lat="dg_wrap.bmad", root="dg1", extra=""), threads="1")

    def h5_identical(fa, fb):
        with h5py.File(fa) as a, h5py.File(fb) as b:
            names_a, names_b = [], []
            a.visititems(lambda n, o: names_a.append(n) if isinstance(o, h5py.Dataset) else None)
            b.visititems(lambda n, o: names_b.append(n) if isinstance(o, h5py.Dataset) else None)
            if sorted(names_a) != sorted(names_b):
                return False
            return all(np.array_equal(a[n][()], b[n][()]) for n in names_a)

    same = all(h5_identical(wd / f"dg{s}", wd / f"dg1{s}")
               for s in (".stats.h5", "-escaped.fld.h5", "-pulse.fld.h5"))
    check("thread invariance: stats/escaped/pulse data identical 1 vs 8 (1 = yes)",
          0.0 if same else 1.0, 0.5)

    # 6. Refusal: dump list entry matching nothing, refused by name.
    (wd / "dg_bad.nml").write_text(NML.format(lat="dg_wrap.bmad", root="dg_bad",
                                              extra='  dump_beam_at(2) = "NO_SUCH_ELEMENT"\n'))
    r = subprocess.run([str(exe), "dg_bad.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = r.returncode != 0 and "NO_SUCH_ELEMENT" in (r.stdout + r.stderr)
    check("refusal: unknown dump_beam_at element refused by name (1 = yes)",
          0.0 if refused else 1.0, 0.5)

    if FAILED:
        print("DIAGNOSTIC CHECKS: FAIL")
        sys.exit(1)
    print("diagnostic checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
