#!/usr/bin/env python3
"""
The saturation-demo figure: one SASE case run to saturation three ways -- Genesis 1.3
Version 4, the Bmad FEL tracker's averaged mode, and its unaveraged mode -- from
IDENTICAL initial dumps, overlaid with comparative statistics and wall-clock timing.
This is the figure run_saturation_demo.sh ends with; see that script for the
methodology (identical particles and field, one external clock, answers verified
before any timing is quoted).

Six panels:
  1. Total power vs z (log): the gain curve through saturation, all three codes,
     with the wall-clock annotation. Agreement here is the headline.
  2. Relative total-power difference vs z, each Bmad mode against Genesis (log).
  3. Power per slice at the exit: the SASE spike structure across the time window.
  4. Slice-averaged bunching |b| vs z.
  5. Mean energy change <gamma> - gamma(0) vs z: the energy the light took.
  6. Slice-averaged rms energy spread vs z.

Usage:
  plot_fel_saturation.py <genesis .out.h5> <avg diag.txt> <unavg diag.txt>
      [--times G,A,U] [--workers N] [-o out.png]

--times takes the three wall-clock seconds (Genesis, averaged, unaveraged); --workers
the core count, for the annotation only.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_fel import load_fortran_diag, load_genesis_out  # noqa: E402

M_ELECTRON = 0.51099895069e6

STYLE = {  # identity by entity, fixed order; gray/orange/blue survive CVD checks
    "genesis": dict(color="0.25", lw=2.6, label="Genesis4 (MPI)"),
    "avg": dict(color="tab:orange", lw=1.4, label="Bmad averaged"),
    "unavg": dict(color="tab:blue", lw=1.4, label="Bmad unaveraged"),
}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("genesis_out")
    p.add_argument("avg_diag")
    p.add_argument("unavg_diag")
    p.add_argument("--times", default=None,
                   help="wall seconds as G,A,U (Genesis, averaged, unaveraged)")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("-o", "--out", default="sat-demo.png")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with h5py.File(args.genesis_out) as h5:
        nslice = h5["Field/power"].shape[1]
        espread_g = h5["Beam/energyspread"][:]     # gamma units, (nrec, nslice)

    g = load_genesis_out(args.genesis_out, nslice)
    a = load_fortran_diag(args.avg_diag, nslice)
    u = load_fortran_diag(args.unavg_diag, nslice)
    # sigma_energy is diag column 7 (eV); load it directly for the spread panel.
    sig_a = np.loadtxt(args.avg_diag).reshape(-1, nslice, 12)[:, :, 7] / M_ELECTRON
    sig_u = np.loadtxt(args.unavg_diag).reshape(-1, nslice, 12)[:, :, 7] / M_ELECTRON

    n = min(len(g["z"]), len(a["z"]), len(u["z"]))
    z = g["z"][:n]
    total = {k: d["power"][:n].sum(axis=1) for k, d in (("genesis", g), ("avg", a), ("unavg", u))}
    bunch = {k: d["bunching"][:n].mean(axis=1) for k, d in (("genesis", g), ("avg", a), ("unavg", u))}
    # Flat-current window: the plain slice mean IS the bunch mean.
    dgam = {k: d["energy"][:n].mean(axis=1) - d["energy"][0].mean()
            for k, d in (("genesis", g), ("avg", a), ("unavg", u))}
    spread = {"genesis": espread_g[:n].mean(axis=1),
              "avg": sig_a[:n].mean(axis=1), "unavg": sig_u[:n].mean(axis=1)}

    fig, ax = plt.subplots(2, 3, figsize=(15, 8.6), constrained_layout=True)
    fig.suptitle(f"SASE to saturation, one initial condition, three trackers "
                 f"({nslice} slices × {z[-1]:.1f} m)")

    A = ax[0, 0]
    for k in STYLE:
        A.semilogy(z, np.maximum(total[k], 1e-30), **STYLE[k])
    pos = total["genesis"][total["genesis"] > 0]
    if len(pos):
        A.set_ylim(bottom=max(pos.min() * 0.3, total["genesis"].max() * 1e-10),
                   top=total["genesis"].max() * 4)
    A.set_xlabel("z (m)"); A.set_ylabel("total power (W)")
    A.legend(loc="lower right")
    if args.times:
        tg, ta, tu = [float(x) for x in args.times.split(",")]
        wtxt = f" on {args.workers} cores" if args.workers else ""
        A.text(0.03, 0.97,
               f"wall clock{wtxt}:\n"
               f"Genesis4 {tg:.0f} s (MPI)\n"
               f"Bmad averaged {ta:.0f} s ({tg/ta:.1f}× vs G4)\n"
               f"Bmad unaveraged {tu:.0f} s ({tu/ta:.0f}× averaged)",
               transform=A.transAxes, va="top", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    A = ax[0, 1]
    den = np.maximum(total["genesis"], 1e-9 * total["genesis"].max())
    for k in ("avg", "unavg"):
        A.semilogy(z, np.maximum(np.abs(total[k] - total["genesis"]) / den, 1e-12),
                   color=STYLE[k]["color"], lw=1.4, label=STYLE[k]["label"])
    A.set_xlabel("z (m)"); A.set_ylabel("|P/P_Genesis − 1|, total power")
    A.legend()

    A = ax[0, 2]
    s = np.arange(1, nslice + 1)
    for k, d in (("genesis", g), ("avg", a), ("unavg", u)):
        A.step(s, d["power"][n - 1], where="mid", **STYLE[k])
    A.set_xlabel("slice (time-window order)")
    A.set_ylabel(f"power at z = {z[-1]:.1f} m (W)")
    A.legend()

    A = ax[1, 0]
    for k in STYLE:
        A.plot(z, bunch[k], **STYLE[k])
    A.set_xlabel("z (m)"); A.set_ylabel("bunching |b| (slice mean)")
    A.legend()

    A = ax[1, 1]
    for k in STYLE:
        A.plot(z, dgam[k], **STYLE[k])
    A.set_xlabel("z (m)"); A.set_ylabel("⟨γ⟩ − γ(0)")
    A.legend()

    A = ax[1, 2]
    for k in STYLE:
        A.plot(z, spread[k], **STYLE[k])
    A.set_xlabel("z (m)"); A.set_ylabel("rms energy spread σ_γ (slice mean)")
    A.legend()

    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")

    import math
    isat = int(np.argmax(total["genesis"]))
    print(f"  Genesis saturation:  P = {total['genesis'][isat]:.3e} W at z = {z[isat]:.1f} m; "
          f"gain from first record {total['genesis'][isat]/total['genesis'][1]:.1e}")
    for k in ("avg", "unavg"):
        print(f"  {STYLE[k]['label']:16s} exit P = {total[k][-1]:.3e} W "
              f"(Genesis {total['genesis'][-1]:.3e}); |ln ratio| = "
              f"{abs(math.log(total[k][-1] / total['genesis'][-1])):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
