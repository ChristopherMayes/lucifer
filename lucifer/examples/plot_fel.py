#!/usr/bin/env python3
"""
Plot the standard diagnostics of a fel_track_test run from its .diag.txt file.

    python plot_fel.py steady_state.diag.txt            # writes steady_state.png
    python plot_fel.py run.diag.txt -o gain.png --show

Seven panels against z: radiation power (log and linear), FIELD ENERGY in the
window (log and linear,
joules -- the honest growth curve for SASE, since per-slice power fluctuates as
radiation slips through and out of the window while the window energy grows
smoothly), bunching, beam energy change and rms spread (MeV -- Bmad's convention:
energy means eV, never gamma), and transverse rms beam sizes. Works on
steady-state and time-dependent files alike: with more than one slice, thin gray
lines show every slice and the bold line the total (power, energy) or the plain
slice average (everything else). The field-energy scaling U = P*slice_spacing/c
reads slice_spacing from the diag header.

Needs numpy and matplotlib (the bmad-fel-validate environment has both).
"""

from __future__ import annotations

import argparse
import pathlib

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fixed series colors (validated categorical palette, assigned in order, never
# cycled); per-slice context lines are neutral gray so identity stays with the
# bold aggregate.
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#b0afa9"
INK, INK2 = "#0b0b0b", "#52514e"


def load(fn):
    """Return z (nrec,), per-quantity arrays (nrec, nslice) in time order, nslice,
    and slice_spacing [m] from the header (None for pre-header files)."""
    spacing = None
    with open(fn) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "slice_spacing" in line:
                spacing = float(line.split("=")[1])
    if spacing is None:
        raise SystemExit(f"{fn}: no slice_spacing header; rerun with a current fel_track_test")
    d = np.loadtxt(fn)
    if d.ndim == 1:
        d = d[None, :]
    nslice = int(d[:, 1].max())
    # An interrupted run leaves a partial trailing record; plot the complete ones.
    n_full = (d.shape[0] // nslice) * nslice
    if n_full != d.shape[0]:
        print(f"note: file ends mid-record (interrupted run?); "
              f"plotting {d.shape[0]//nslice} complete records, dropping "
              f"{d.shape[0]-n_full} trailing rows")
        d = d[:n_full]
    d = d.reshape(-1, nslice, d.shape[1])
    q = {name: d[:, :, i] for i, name in enumerate(
        ("z", "slice", "power", "on_axis", "bunching", "bunching_phase",
         "mean_energy", "sigma_energy", "sigma_x", "sigma_y"))}
    return d[:, 0, 0], q, nslice, spacing


def panel_series(ax, z, a2d, color, label=None, reduce="mean"):
    """One quantity: thin gray per-slice lines (if several) plus a bold aggregate."""
    nslice = a2d.shape[1]
    if nslice > 1:
        ax.plot(z, a2d, color=GRAY, lw=0.6, alpha=0.6, zorder=1)
        agg = a2d.sum(axis=1) if reduce == "sum" else a2d.mean(axis=1)
    else:
        agg = a2d[:, 0]
    ax.plot(z, agg, color=color, lw=1.8, label=label, zorder=3)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("diag", help="fel_track_test .diag.txt file")
    p.add_argument("-o", "--output", help="Output image (default: <diag>.png)")
    p.add_argument("--show", action="store_true", help="Also open a window")
    args = p.parse_args()

    z, q, nslice, spacing = load(args.diag)
    out = args.output or str(pathlib.Path(args.diag).with_suffix(".png"))

    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": INK2, "grid.alpha": 0.18, "grid.linewidth": 0.6,
        "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "font.size": 10,
    })

    fig, axd = plt.subplot_mosaic([["power", "energy"], ["power_lin", "energy_lin"],
                                   ["bunching", "gamma"], ["size", "size"]],
                                  figsize=(9.5, 12), sharex=True)
    ax_p, ax_u, ax_pl, ax_ul, ax_b, ax_g, ax_s = (axd[k] for k in
        ("power", "energy", "power_lin", "energy_lin", "bunching", "gamma", "size"))

    # Radiation power. The one panel where the multi-slice aggregate is a sum:
    # total power is what a detector sees.
    label = "total" if nslice > 1 else None
    panel_series(ax_p, z, q["power"], BLUE, label=label, reduce="sum")
    ax_p.set_yscale("log")
    ax_p.set_ylabel("power (W)")
    ax_p.set_title("Radiation power", loc="left")
    if nslice > 1:
        ax_p.legend(frameon=False)

    # Field energy in the window: U_i = P_i * slice_spacing / c per slice, summed.
    # This is the curve to read for SASE growth: power per slice churns as radiation
    # slips forward through the window (and out of its head), while the window energy
    # integrates it honestly.
    dt = spacing / 2.99792458e8
    label = "window total" if nslice > 1 else None
    panel_series(ax_u, z, q["power"] * dt, BLUE, label=label, reduce="sum")
    ax_u.set_yscale("log")
    ax_u.set_ylabel("field energy (J)")
    ax_u.set_title("Field energy", loc="left")
    if nslice > 1:
        ax_u.legend(frameon=False)

    # The same two quantities on LINEAR axes: the log panels show the exponential
    # gain regime; the linear ones show where the energy actually is (saturation and
    # the post-saturation behavior are nearly invisible on a log axis).
    label = "total" if nslice > 1 else None
    panel_series(ax_pl, z, q["power"], BLUE, label=label, reduce="sum")
    ax_pl.set_ylabel("power (W)")
    ax_pl.set_title("Radiation power (linear)", loc="left")

    label = "window total" if nslice > 1 else None
    panel_series(ax_ul, z, q["power"] * dt, BLUE, label=label, reduce="sum")
    ax_ul.set_ylabel("field energy (J)")
    ax_ul.set_title("Field energy (linear)", loc="left")

    # Bunching factor.
    label = "slice average" if nslice > 1 else None
    panel_series(ax_b, z, q["bunching"], BLUE, label=label)
    ax_b.set_ylabel("|b|")
    ax_b.set_title("Bunching", loc="left")
    if nslice > 1:
        ax_b.legend(frameon=False)

    # Two-series legends sit above their panel, clear of the data (betatron
    # oscillations fill the beam-size panel top to bottom).
    above = dict(frameon=False, ncol=2, loc="lower right",
                 bbox_to_anchor=(1, 0.99), borderaxespad=0)

    # Beam energy: change of the mean and the rms spread share the unit MeV
    # (Bmad's convention -- energy is eV, never gamma), so one axis carries both.
    de = (q["mean_energy"] - q["mean_energy"][0, :]) / 1e6
    panel_series(ax_g, z, de, BLUE, label=r"$\Delta\langle E\rangle$")
    panel_series(ax_g, z, q["sigma_energy"] / 1e6, ORANGE, label=r"$\sigma_E$")
    ax_g.set_ylabel("(MeV)")
    ax_g.set_title("Beam energy", loc="left")
    ax_g.legend(**above)

    # Transverse rms sizes.
    panel_series(ax_s, z, 1e6 * q["sigma_x"], BLUE, label=r"$\sigma_x$")
    panel_series(ax_s, z, 1e6 * q["sigma_y"], ORANGE, label=r"$\sigma_y$")
    ax_s.set_ylabel(r"rms size ($\mu$m)")
    ax_s.set_xlabel("z (m)")
    ax_s.set_title("Beam size", loc="left")
    ax_s.legend(**above)

    title = pathlib.Path(args.diag).name
    if nslice > 1:
        title += f"   ({nslice} slices; thin lines are individual slices)"
    fig.suptitle(title, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    fig.savefig(out, dpi=160)
    print(f"wrote {out}")
    if args.show:
        matplotlib.use("TkAgg", force=False)
        plt.show()


if __name__ == "__main__":
    main()
