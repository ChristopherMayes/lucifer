#!/usr/bin/env python3
"""
Overlay the Bmad FEL tracker and Genesis 1.3 Version 4 outputs of one benchmark tier:
the visual companion to compare_fel.py's numbers. Reads the tracker's diag file and
Genesis's .out.h5 (the slice count comes from the Genesis file, so steady-state and
time-dependent tiers both work unannounced) and writes a four-panel figure:

  1. Total power vs z, both codes overlaid (log). Agreement here is the headline.
  2. The elementwise relative power difference vs z (max and median over slices) --
     the same quantity compare_fel.py gates, resolved along the line.
  3. Final-record power per slice, both codes: the SASE spectrum across the time
     window (or two dots for steady state).
  4. Slice-averaged bunching vs z, both codes overlaid.

Usage:
  plot_fel_compare.py <bmad diag.txt> <genesis .out.h5> [-o out.png]

e.g., from a benchmark work directory:
  plot_fel_compare.py tdsase.diag.txt AramisTDSASE.out.h5
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_fel import load_fortran_diag, load_genesis_out  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bmad_diag", help="fel_track_test diag file (<out_root>.diag.txt)")
    p.add_argument("genesis_out", help="Genesis .out.h5 of the same configuration")
    p.add_argument("-o", "--out", default=None,
                   help="Output figure (default: <bmad_diag stem>.compare.png)")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with h5py.File(args.genesis_out) as h5:
        nslice = h5["Field/power"].shape[1]

    f = load_fortran_diag(args.bmad_diag, nslice)
    g = load_genesis_out(args.genesis_out, nslice)

    n = min(len(f["z"]), len(g["z"]))
    if len(f["z"]) != len(g["z"]):
        print(f"note: record counts differ (bmad {len(f['z'])}, genesis {len(g['z'])}); "
              f"comparing the first {n}")
    z = g["z"][:n]

    if nslice > 1:
        fp, gp = f["power"][:n], g["power"][:n]
        total_f, total_g = fp.sum(axis=1), gp.sum(axis=1)
        bunch_f, bunch_g = f["bunching"][:n].mean(axis=1), g["bunching"][:n].mean(axis=1)
    else:
        fp, gp = f["power"][:n, None], g["power"][:n, None]
        total_f, total_g = f["power"][:n], g["power"][:n]
        bunch_f, bunch_g = f["bunching"][:n], g["bunching"][:n]

    # The gated quantity: elementwise relative power, denominator floored exactly as in
    # compare_fel.py so vacuum slices cannot divide by zero.
    den = np.maximum(np.abs(gp), 1e-9 * np.abs(gp).max())
    rel = np.abs(fp - gp) / den

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    name = pathlib.Path(args.bmad_diag).stem.removesuffix(".diag")
    fig.suptitle(f"{name}: Bmad FEL tracker vs Genesis4  "
                 f"({nslice} slice{'s' if nslice > 1 else ''})")

    a = ax[0, 0]
    a.semilogy(z, np.maximum(total_g, 1e-30), color="0.25", lw=2.4, label="Genesis4")
    a.semilogy(z, np.maximum(total_f, 1e-30), color="tab:orange", lw=1.2, label="Bmad")
    # A dark start's z = 0 record is exactly zero; don't let it stretch the axis.
    pos = total_g[total_g > 0]
    if len(pos):
        a.set_ylim(bottom=max(pos.min() * 0.3, total_g.max() * 1e-14),
                   top=total_g.max() * 3)
    a.set_xlabel("z [m]"); a.set_ylabel("total power [W]")
    a.set_title("Power along the line")
    a.legend()

    a = ax[0, 1]
    a.semilogy(z, rel.max(axis=1), color="tab:red", lw=1.2, label="max over slices")
    if nslice > 1:
        a.semilogy(z, np.median(rel, axis=1), color="tab:red", lw=1.2, ls="--",
                   alpha=0.6, label="median")
    a.set_xlabel("z [m]"); a.set_ylabel("|P_bmad - P_genesis| / P_genesis")
    a.set_title("Relative power difference (the gated quantity)")
    a.legend()

    a = ax[1, 0]
    s = np.arange(1, nslice + 1)
    a.step(s, gp[-1], where="mid", color="0.25", lw=2.4, label="Genesis4")
    a.step(s, fp[-1], where="mid", color="tab:orange", lw=1.2, label="Bmad")
    a.set_xlabel("slice (time-window order)"); a.set_ylabel("power [W]")
    a.set_title(f"Per-slice power at exit (z = {z[-1]:.2f} m)")
    a.legend()

    a = ax[1, 1]
    a.plot(z, bunch_g, color="0.25", lw=2.4, label="Genesis4")
    a.plot(z, bunch_f, color="tab:orange", lw=1.2, label="Bmad")
    a.set_xlabel("z [m]"); a.set_ylabel("bunching |b|" + (" (slice mean)" if nslice > 1 else ""))
    a.set_title("Bunching along the line")
    a.legend()

    out = args.out or f"{name}.compare.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    print(f"  exit total power: bmad {total_f[-1]:.4e} W, genesis {total_g[-1]:.4e} W, "
          f"rel {abs(total_f[-1]-total_g[-1])/total_g[-1]:.2e}")
    print(f"  elementwise relative power, whole run: max {rel.max():.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
