#!/usr/bin/env python3
"""
Plot the standard diagnostics of a lucifer run from its stats file.

    python plot_fel.py run.stats.h5              # writes run.png
    python plot_fel.py run.stats.h5 -o gain.png --show

The stats file (manual sec:stats) carries per-record, per-slice beam MOMENTS (named as
bunch_params_struct components) and wavefront_params, in fixed Bmad units -- everything
below derives from it, no text files involved. Ten panels against z:

  radiation power and window field energy (log and linear -- the log pair shows the
  exponential gain regime, the linear pair shows where the energy actually is);
  bunching |b|; beam energy change and rms spread (MeV -- Bmad's convention: energy is
  eV, never gamma); rms BEAM sizes and rms FIELD sizes (the field sizes come from the
  wavefront_params sigma(4,4) -- watch gain guiding pull the light onto the beam);
  beam normalized emittances (projected, dispersion removed, the bunch_params
  convention) and the field "emittance" sqrt(det sigma_plane) = M^2 lambda/4pi at the
  records where the FFT-costed angle moments were taken (element ends).

Works on steady-state and time-dependent runs alike: with more than one slice, thin
gray lines show every slice and the bold line the total (power, energy) or the slice
average (everything else). When the run carries TWO POLARIZATIONS (any tilted FEL
element; the stats file then has a field/y group) an eleventh panel splits the power
by polarization -- on a crossed line that panel IS the afterburner story. When it
carries HARMONIC FIELDS (namelist harmonics; field/harm<h> groups) a panel plots each
harmonic's power against the fundamental's.

Needs numpy, h5py and matplotlib (the bmad-fel-validate environment has all three).
"""

from __future__ import annotations

import argparse
import pathlib

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

M_ELECTRON = 0.51099895069e6   # Bmad's m_electron [eV]

# Fixed series colors (validated categorical palette, assigned in order, never
# cycled). Per-slice context lines are neutral gray so identity stays with the
# bold aggregate.
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#b0afa9"
INK, INK2 = "#0b0b0b", "#52514e"


def load(fn):
    """Everything the panels need, all (nrec, nslice) except z and p0c."""
    with h5py.File(fn) as h5:
        q = {
            "z": h5["z"][:], "p0c": float(h5["p0c"][0]),
            "centroid": h5["beam/centroid"][:],          # (nrec, nslice, 6)
            "sigma": h5["beam/sigma"][:],                # (nrec, nslice, 36)
            "bunching": h5["beam/bunching"][:],
            "power": h5["field/power"][:],
            "energy": h5["field/energy"][:],
            "f_sigma": h5["field/sigma"][:],             # (nrec, nslice, 16)
            "f_emit_x": h5["field/emit_x"][:],
            "f_emit_y": h5["field/emit_y"][:],
            "f_valid": h5["field/angle_moments_valid"][:].astype(bool),
        }
        # Two live polarizations: field/power is the TOTAL and field/* the x
        # component. Field/y/ carries the second. Absent on single-polarization runs.
        if "field/y" in h5:
            q["power_y"] = h5["field/y/power"][:]
            q["power_x"] = q["power"] - q["power_y"]
        # Harmonic fields: field/harm<h>/ groups, one per harmonic beyond the
        # fundamental. field/power stays the FUNDAMENTAL's -- harmonics are separate
        # wavelengths, never summed with it.
        q["harmonics"] = {}
        for k in h5["field"]:
            if k.startswith("harm"):
                q["harmonics"][int(k[4:])] = h5[f"field/{k}/power"][:]
    return q


def beam_energy_ev(q):
    """Mean total energy per slice [eV] from the stored chart: E = sqrt(P^2 + m^2)."""
    p_hat = q["p0c"] * (1.0 + q["centroid"][:, :, 5])
    return np.sqrt(p_hat**2 + M_ELECTRON**2)


def norm_emit(q, i, j):
    """Projected normalized emittance per (record, slice), dispersion removed --
    the bunch_params convention (Bmad's projected_twiss_calc)."""
    s = q["sigma"].reshape(q["sigma"].shape[0], -1, 6, 6)
    s66 = s[:, :, 5, 5]
    with np.errstate(divide="ignore", invalid="ignore"):
        x2 = s[:, :, i, i] - np.where(s66 > 0, s[:, :, i, 5] ** 2 / s66, 0)
        xp = s[:, :, i, j] - np.where(s66 > 0, s[:, :, i, 5] * s[:, :, j, 5] / s66, 0)
        p2 = s[:, :, j, j] - np.where(s66 > 0, s[:, :, j, 5] ** 2 / s66, 0)
    emit = np.sqrt(np.maximum(0.0, x2 * p2 - xp**2))
    f_emit = q["p0c"] * (1.0 + q["centroid"][:, :, 5]) / M_ELECTRON
    return f_emit * emit


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
    p.add_argument("stats", help="lucifer .stats.h5 file")
    p.add_argument("-o", "--output", help="Output image (default: <stats stem>.png)")
    p.add_argument("--show", action="store_true", help="Also open a window")
    args = p.parse_args()

    if args.stats.endswith(".diag.txt"):
        raise SystemExit("plot_fel.py reads the stats file now; pass <out_root>.stats.h5 "
                         "(the diag file remains the Genesis-comparison instrument).")

    q = load(args.stats)
    z = q["z"]
    nslice = q["power"].shape[1]
    out = args.output or str(pathlib.Path(args.stats).name.removesuffix(".stats.h5") + ".png")

    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": INK2, "grid.alpha": 0.18, "grid.linewidth": 0.6,
        "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "font.size": 10,
    })

    rows = [["power", "energy"], ["power_lin", "energy_lin"], ["bunching", "gamma"],
            ["bsize", "fsize"], ["bemit", "femit"]]
    two_pol = "power_y" in q
    if two_pol:
        rows.append(["pol", "pol"])
    if q["harmonics"]:
        rows.append(["harm", "harm"])
    fig, axd = plt.subplot_mosaic(rows, figsize=(9.5, 14 + 2.6 * (two_pol + bool(q["harmonics"]))),
                                  sharex=True)

    # Radiation power and window field energy, log then linear. The multi-slice
    # aggregate is a sum: total power is what a detector sees, and the window energy
    # is the honest SASE growth curve (per-slice power churns as light slips through
    # and out of the window, while the energy integrates it).
    for key, quant, ylab in (("power", "power", "radiation power (W)"),
                             ("energy", "energy", "field energy (J)")):
        label = "total" if nslice > 1 else None
        panel_series(axd[key], z, q[quant], BLUE, label=label, reduce="sum")
        axd[key].set_yscale("log")
        axd[key].set_ylabel(ylab)
        if nslice > 1:
            axd[key].legend(frameon=False)
        panel_series(axd[key + "_lin"], z, q[quant], BLUE, reduce="sum")
        axd[key + "_lin"].set_ylabel(ylab)

    # Bunching factor.
    label = "slice average" if nslice > 1 else None
    panel_series(axd["bunching"], z, q["bunching"], BLUE, label=label)
    axd["bunching"].set_ylabel("bunching |b|")
    if nslice > 1:
        axd["bunching"].legend(frameon=False)

    above = dict(frameon=False, ncol=2, loc="lower right",
                 bbox_to_anchor=(1, 0.99), borderaxespad=0)

    # Beam energy: change of the mean and the rms spread share the unit MeV.
    e_mean = beam_energy_ev(q)
    sig = q["sigma"].reshape(len(z), nslice, 6, 6)
    beta = q["p0c"] * (1 + q["centroid"][:, :, 5]) / e_mean
    sig_e = beta * q["p0c"] * np.sqrt(np.maximum(0.0, sig[:, :, 5, 5]))
    panel_series(axd["gamma"], z, (e_mean - e_mean[0]) / 1e6, BLUE, label=r"$\Delta\langle E\rangle$")
    panel_series(axd["gamma"], z, sig_e / 1e6, ORANGE, label=r"$\sigma_E$")
    axd["gamma"].set_ylabel("beam energy (MeV)")
    axd["gamma"].legend(**above)

    # Transverse rms sizes: beam from the 6x6, FIELD from the wavefront 4x4 --
    # gain guiding is the field-size curve bending toward the beam-size curve.
    panel_series(axd["bsize"], z, 1e6 * np.sqrt(np.maximum(0, sig[:, :, 0, 0])), BLUE, label=r"$\sigma_x$")
    panel_series(axd["bsize"], z, 1e6 * np.sqrt(np.maximum(0, sig[:, :, 2, 2])), ORANGE, label=r"$\sigma_y$")
    axd["bsize"].set_ylabel(r"rms beam size ($\mu$m)")
    axd["bsize"].legend(**above)

    fsig = q["f_sigma"].reshape(len(z), nslice, 4, 4)
    panel_series(axd["fsize"], z, 1e6 * np.sqrt(np.maximum(0, fsig[:, :, 0, 0])), BLUE, label=r"$\sigma_x$")
    panel_series(axd["fsize"], z, 1e6 * np.sqrt(np.maximum(0, fsig[:, :, 2, 2])), ORANGE, label=r"$\sigma_y$")
    axd["fsize"].set_ylabel(r"rms field size ($\mu$m)")
    axd["fsize"].legend(**above)

    # Emittances: beam normalized (projected, dispersion removed). Field
    # sqrt(det sigma_plane) at the records where angle moments exist (element ends),
    # against the diffraction limit lambda/4pi.
    panel_series(axd["bemit"], z, 1e6 * norm_emit(q, 0, 1), BLUE, label=r"$\gamma\epsilon_x$")
    panel_series(axd["bemit"], z, 1e6 * norm_emit(q, 2, 3), ORANGE, label=r"$\gamma\epsilon_y$")
    axd["bemit"].set_ylabel(r"norm. emittance ($\mu$m)")
    axd["bemit"].legend(**above)

    valid = q["f_valid"].any(axis=1)
    if valid.any():
        zv = z[valid]
        ex = np.where(q["f_valid"][valid], q["f_emit_x"][valid], np.nan)
        ey = np.where(q["f_valid"][valid], q["f_emit_y"][valid], np.nan)
        axd["femit"].plot(zv, 1e12 * np.nanmean(ex, axis=1), "o-", color=BLUE, ms=4, label=r"$\epsilon_x$")
        axd["femit"].plot(zv, 1e12 * np.nanmean(ey, axis=1), "o-", color=ORANGE, ms=4, label=r"$\epsilon_y$")
    axd["femit"].set_ylabel(r"field emit (pm$\,$rad)")
    axd["femit"].legend(**above)

    axd["bemit"].set_xlabel("z (m)")
    axd["femit"].set_xlabel("z (m)")

    # The polarization split (two live components only): x and y power against z. On a
    # crossed line this is the afterburner -- the x set amplifies Ex and bunches the
    # beam, then the tilted set radiates Ey from that bunching while Ex only diffracts.
    if two_pol:
        A = axd["pol"]
        A.semilogy(z, np.maximum(q["power_x"][:len(z)].sum(axis=1), 1e-30),
                   color=BLUE, lw=1.8, label="x polarization")
        A.semilogy(z, np.maximum(q["power_y"][:len(z)].sum(axis=1), 1e-30),
                   color=ORANGE, lw=1.8, label="y polarization")
        A.set_xlabel("z (m)"); A.set_ylabel("radiation power by polarization (W)")
        A.legend(frameon=False)
        axd["bemit"].set_xlabel("")
        axd["femit"].set_xlabel("")

    # Harmonic power (harmonic field-set runs): each harmonic against the
    # fundamental, on the fundamental's log scale -- harmonic lasing reads directly
    # off the vertical gaps.
    if q["harmonics"]:
        A = axd["harm"]
        A.semilogy(z, np.maximum(q["power"][:len(z)].sum(axis=1), 1e-30),
                   color=BLUE, lw=1.8, label="fundamental")
        for hh, ph in sorted(q["harmonics"].items()):
            A.semilogy(z, np.maximum(ph[:len(z)].sum(axis=1), 1e-30),
                       color=ORANGE, lw=1.8, label=f"harmonic {hh}")
        A.set_xlabel("z (m)"); A.set_ylabel("radiation power by harmonic (W)")
        A.legend(frameon=False)
        axd["bemit"].set_xlabel("")
        axd["femit"].set_xlabel("")

    title = pathlib.Path(args.stats).name
    if nslice > 1:
        title += f"   ({nslice} slices; thin lines are individual slices)"
    fig.suptitle(title, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    fig.savefig(out, dpi=160)
    print(f"wrote {out}")
    if args.show:
        matplotlib.use("TkAgg", force=False)
        plt.show()


if __name__ == "__main__":
    main()
