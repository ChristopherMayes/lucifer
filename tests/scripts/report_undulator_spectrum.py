#!/usr/bin/env python3
"""
Draw the figures of docs/undulator-spectrum.md and write the page of figures, captions and
the table that it includes.

  report_undulator_spectrum.py --work-dir <kept work directory> --out docs/generated/undulator-spectrum
                               --page docs/generated/undulator-spectrum.md

The curves are those check_undulator_spectrum.py wrote to undulator_spectrum.json in its
section's directory of a kept work directory of the validation suite. The formulas drawn
under them are tests/scripts/undulator_radiation.py's.

    line.svg      the line on the axis of the three undulators over sinc^2
    pulse.svg     the zero-angle component of each plane of the pulse, in the averaged
                  mode and in the unaveraged one
    ridge.svg     the frequency of the ridge at each transverse wavenumber over -x and over the
                  resonance of undulator radiation, with the spectral propagator's beside it
    angular.svg   the level along the ridge, the model's and the formula's
    filter.svg    the energy at each angle with run.source_filter on, over the filter off
    beam.svg      the line of one slice with FLASH1's emittance and with its energy spread
    bunch.svg     the mean spectrum of a bunch of 1024 slices over the line of one electron

Each figure is written through lucifer.plot.save, so the same curves give the same bytes.
The caption of each states the largest differences at two significant figures, computed
from the curves the figure draws. The keystone's regeneration runs this from the
production pass's kept work directory and diffs the page with the other generated pages,
so a curve that moves fails the keystone through its caption, and the figures stay out of
the diff.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import undulator_radiation as ur  # noqa: E402
from check_undulator_spectrum import FLASH, resonance  # noqa: E402

SIZE = (6.4, 3.8)


def two(x):
    """A number at two significant figures, as the captions print it."""
    return f"{x:.2g}"


def undulator(und, n_period):
    return ur.Undulator(und["e_tot"], und["b_max"], und["l_period"], n_period, und["helical"])


def block(name, lines):
    """One block of the page, between the markers docs/undulator-spectrum.md includes by."""
    return [f"<!-- begin {name} -->", *lines, f"<!-- end {name} -->", ""]


def figure(name, alt, caption):
    return block(f"fig-und-{name}", [f":::{{figure}} generated/undulator-spectrum/{name}.svg",
                                     f":name: fig-und-{name}", f":alt: {alt}", "", caption, ":::"])


def draw_line(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    u = np.linspace(-4, 4, 1601)
    ax.semilogy(u, np.maximum(ur.line(1, u), 1e-12), color=lp.INK, lw=1.0, label=r"$\mathrm{sinc}^2(\pi N\nu)$")
    worst = {}
    for k, (key, label) in enumerate((("us_flash", "FLASH1, $N = 165$"), ("us_flash55", "FLASH1, $N = 55$"),
                                      ("us_delta", "Delta, $N = 100$"))):
        r = rec[key]
        x = np.array(r["nu"]) * r["n_period"]
        y = np.array(r["line"])
        ax.semilogy(x[k * 3::16], np.maximum(y[k * 3::16], 1e-12), "o", ms=3.2, color=lp.COLORS[k], label=label)
        worst[key] = float(np.abs(y - ur.line(r["n_period"], np.array(r["nu"]))).max())
    # The unaveraged mode's pulse does not consist of equal planes, so its line is drawn
    # whole: its minima are 4% outside the zeros, which sparse points would not show.
    r = rec["us_unavg"]
    y = np.array(r["line"])
    ax.semilogy(np.array(r["nu"]) * r["n_period"], np.maximum(y, 1e-12), lw=0.9, color=lp.COLORS[3],
                label="FLASH1, $N = 55$, unaveraged")
    worst["us_unavg"] = float(np.abs(y - ur.line(r["n_period"], np.array(r["nu"]))).max())
    ax.set_ylim(1e-7, 2)
    ax.set_xlabel(r"$N\nu$")
    ax.set_ylabel(r"$W(\nu)/W(0)$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    return worst


def draw_pulse(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    for k, (key, label) in enumerate((("us_flash55", '`fel_method = averaged`'), ("us_unavg", '`fel_method = unaveraged`'))):
        r = rec[key]
        y = np.array(r["planes"])
        p = np.arange(len(y)) - r["first_plane"]
        sel = (p >= -4) & (p <= r["n_period"] + 4)
        ax.plot(p[sel], y[sel], "o-" if k else "s-", ms=3, lw=0.8, color=lp.COLORS[k], label=label.strip("`"))
    ax.set_xlabel(r"$p$")
    ax.set_ylabel(r"$|E_0(p)|/\max|E_0|$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    un = np.array(rec["us_unavg"]["planes"])
    lit = un[un > 1e-6]
    return float(lit[0]), float(lit[1])


RIDGES = (("us_flash", "ridge_x", "FLASH1, along $x$", 165), ("us_flash", "ridge_y", "FLASH1, along $y$", 165),
          ("us_delta", "ridge_x", "Delta", 100), ("us_unavg", "ridge_x", "FLASH1, unaveraged, along $x$", 55))
X_LABEL = r"$x = \gamma_z^2 (k_\perp/k_0)^2$"


def draw_ridge(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    x = np.linspace(0, 0.21, 200)
    ax.plot(x, -x, color=lp.INK, lw=1.0, label=r"$-x$")
    ax.plot(x, resonance(x), color=lp.INK2, lw=1.0, ls="--", label=r"$\nu(1 + \nu) = -x$")
    worst = {}
    for k, (key, axis, label, n) in enumerate(RIDGES):
        rows = np.array(rec[key][axis])
        nu = rows[:, 1] - rows[0, 1]
        ax.plot(rows[::2, 0], nu[::2], "osD^"[k], ms=3.4, color=lp.COLORS[k], label=label)
        worst[(key, axis)] = float(np.abs(nu + rows[:, 0]).max() * n / ur.LINE_WIDTH)
    rows = np.array(rec["us_spectral"]["ridge_x"])
    nu = rows[:, 1] - rows[0, 1]
    ax.plot(rows[::2, 0], nu[::2], "v", ms=3.6, color=lp.COLORS[4], label='FLASH1, 55 periods, `field.propagator = "spectral"`'.replace("`", ""))
    worst["spectral"] = float(np.abs(nu - resonance(rows[:, 0])).max() * 55 / ur.LINE_WIDTH)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(r"$\nu$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    return worst


def draw_angular(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    last = {}
    for k, (key, axis, label, n) in enumerate(RIDGES):
        rows = np.array(rec[key][axis])
        ax.plot(rows[::2, 0], (rows[:, 2] / rows[:, 3])[::2], "osD^"[k], ms=3.4, color=lp.COLORS[k], label=label)
        if key != "us_unavg":
            ax.plot(rows[:, 0], rows[:, 4], color=lp.COLORS[k], lw=1.0)
        last[(key, axis)] = (float(rows[-1, 0]), float(rows[-1, 4]), float(rows[-1, 2] / rows[-1, 3]))
    ax.set_ylim(0, 1.15)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(r"$w(k_\perp)/w(0)$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    return last


def draw_filter(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    r = rec["us_filter"]
    rows = np.array(r["kept_x"])
    sel = rows[:, 0] < 3
    ax.plot(rows[sel, 0], rows[sel, 1], "o", ms=3.4, color=lp.COLORS[0])
    ax.axvline(1.0, color=lp.INK2, lw=0.8, ls="--")
    ax.set_xlabel(r"$\theta/\theta_f$")
    ax.set_ylabel(r"$W_\mathrm{on}(\theta)/W_\mathrm{off}(\theta)$")
    lp.save(fig, path)
    plt.close(fig)
    both = np.vstack([np.array(r["kept_x"]), np.array(r["kept_y"])])
    u = undulator(FLASH, 165)
    return (float(np.abs(both[both[:, 0] < 0.8, 1] - 1).max()), float(both[both[:, 0] > 1.2, 1].max()),
            u.angular_factor(r["edge"], 0.0), u.angular_factor(0.0, r["edge"]))


def draw_beam(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    u = np.linspace(-4, 4, 1601)
    ax.plot(u, ur.line(1, u), color=lp.INK, lw=1.0, label=r"$\mathrm{sinc}^2(\pi N\nu)$")
    shift = {}
    for k, (key, label) in enumerate((("us_size", "`beam.a_norm_emit = 1.5e-6`"), ("us_spread", "`beam.sig_pz = 6e-4`"))):
        r = rec[key]
        move = (r["center"] - r["detuning"]) * r["n_period"]
        ax.plot(np.array(r["nu"]) * r["n_period"] + move, r["line"], color=lp.COLORS[k], lw=1.2, label=label.strip("`"))
        shift[key] = (move / ur.LINE_WIDTH, r["width"])
    ax.set_xlim(-4, 4)
    ax.set_xlabel(r"$N\nu$")
    ax.set_ylabel(r"$W(\nu)/\max W$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    return shift


def draw_bunch(plt, lp, path, rec):
    fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
    r = rec["bunch"]
    n = r["n_period"]
    nu = np.array(r["nu"])
    ax.plot(nu * n, ur.sampled_line(n, nu), color=lp.INK, lw=1.0, label="one electron, times the electrons")
    ax.plot(nu * n, r["line"], "o", ms=3, color=lp.COLORS[0], label=f"mean of {r['seeds']} seeds")
    ax.set_xlabel(r"$N\nu$")
    ax.set_ylabel(r"$W(\nu)/(N_e N^2 W_1(0))$")
    ax.legend(frameon=False, fontsize=8)
    lp.save(fig, path)
    plt.close(fig)
    e = np.array(r["energy"])
    return float(e.mean()), float(e.std())


def table(rec):
    """The table of the three undulators' lines."""
    head = ["```{table} The line on the axis of each undulator, from one loaded slice with no beam size and no "
            "energy spread. The width is the full width at half maximum in units of $0.8859/N$, the level is "
            "$W(0)$ divided by $N^2 W_1(0)\\,|S|^2$, and the last two columns are the largest values of "
            "$|W(\\nu)/W(0) - \\mathrm{sinc}^2(\\pi N\\nu)|$ for $|N\\nu| \\le 4$ and of $W(m/N)/W(0)$ for "
            "$m = 1, 2, 3$.",
            ":name: tab-und-line", "",
            "| undulator | $N$ | width | level | difference from $\\mathrm{sinc}^2$ | at the zeros |", "|---|---|---|---|---|---|"]
    rows = []
    for key, name in (("us_flash", "FLASH1, planar"), ("us_flash55", "FLASH1, planar"), ("us_delta", "Delta, helical"),
                      ("us_unavg", "FLASH1, planar, `fel_method = unaveraged`")):
        r = rec[key]
        nu = np.array(r["nu"])
        y = np.array(r["line"])
        n = r["n_period"]
        diff = float(np.abs(y - ur.line(n, nu)).max())
        zeros = max(float(np.interp(m / n, nu, y)) for m in (1, 2, 3))
        rows.append(f"| {name} | {n} | {r['width']:.4f} | {r['level']:.4f} | {two(diff)} | {two(zeros)} |")
    return block("tab-und-line", head + rows + ["```"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--work-dir", required=True, help="a kept work directory of run_fel_validation.sh")
    ap.add_argument("--out", required=True, help="where the figures are written")
    ap.add_argument("--page", required=True, help="the Markdown page of figures, captions and the table")
    args = ap.parse_args()
    record = pathlib.Path(args.work_dir) / "sections" / "undulator-spectrum" / "undulator_spectrum.json"
    if not record.is_file():
        raise SystemExit(f"{record}: no record. The undulator-spectrum section of the validation suite writes it.")
    rec = json.loads(record.read_text())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import lucifer.plot as lp

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with lp._style():
        line = draw_line(plt, lp, out / "line.svg", rec)
        ends = draw_pulse(plt, lp, out / "pulse.svg", rec)
        ridge = draw_ridge(plt, lp, out / "ridge.svg", rec)
        angular = draw_angular(plt, lp, out / "angular.svg", rec)
        kept, cut, edge_x, edge_y = draw_filter(plt, lp, out / "filter.svg", rec)
        shift = draw_beam(plt, lp, out / "beam.svg", rec)
        mean, std = draw_bunch(plt, lp, out / "bunch.svg", rec)

    f, d, un = rec["us_flash"], rec["us_delta"], rec["us_unavg"]
    fx, fy, dx = angular[("us_flash", "ridge_x")], angular[("us_flash", "ridge_y")], angular[("us_delta", "ridge_x")]
    ux = angular[("us_unavg", "ridge_x")]
    b = rec["bunch"]
    flash_last = np.array(f["ridge_x"])[-1]
    delta_last = np.array(d["ridge_x"])[-1]
    lines = ["<!-- Generated by tests/scripts/report_undulator_spectrum.py from the undulator-spectrum section's",
             "     record in a validation pass's work directory. Do not edit. The figures beside it are redrawn",
             "     by the same run, and each caption's numbers are computed from the curves its figure draws. -->", ""]
    lines += figure("line", "The line on the axis of three undulators on a logarithmic axis, as a function of the "
                    "relative frequency times the number of periods, from minus four to four. The points of all "
                    "three undulators lie on the curve of sinc squared, through its central lobe and three side "
                    "lobes on each side. A fourth curve, of the unaveraged mode, follows the central lobe and has "
                    "its minima displaced outward from the zeros of sinc squared by four percent, and its side lobes "
                    "are higher at positive frequency than at negative frequency.",
                    "The energy radiated on the axis per unit frequency and solid angle, divided by its value at the "
                    "peak, as a function of $N\\nu$, where $\\nu = \\omega/\\omega_1 - 1$ and $N$ is the number of "
                    "periods. Each point is the squared magnitude of the Fourier sum, over the planes of the pulse, of the "
                    "field of one run summed over $x$ and $y$ in each plane. The curve is $\\mathrm{sinc}^2(\\pi N\\nu)$. "
                    "The pulse of each run consists of $N$ equal planes, of which this curve is the spectrum. "
                    f"The largest difference from the curve is {two(line['us_flash'])} for FLASH1 at 165 periods, "
                    f"{two(line['us_flash55'])} at 55 periods and {two(line['us_delta'])} for the Delta undulator. "
                    "The difference is that of a pulse sampled once per wavelength, $(\\pi\\nu)^2/3$ of the line "
                    "([](#und-on-axis)). The continuous colored curve is the run of `fel_method = unaveraged` at 55 periods, "
                    "whose pulse does not consist of equal planes ([](#fig-und-pulse)). It differs from "
                    f"$\\mathrm{{sinc}}^2$ by {two(line['us_unavg'])} of the peak: its minima lie outside the zeros of the line, and its side lobes are higher "
                    "above the fundamental than below it ([](#und-unaveraged)).")
    lines += figure("pulse", "The magnitude of the zero-angle component of each plane of the pulse as a function of "
                    "the plane's position, for a 55-period undulator. In the averaged mode it is one value over 55 "
                    "planes and zero outside them. In the unaveraged mode the first and the last two planes are lower.",
                    "The magnitude of the zero-angle component of each plane of the pulse, $E_0(p)$, divided by its "
                    "largest value, as a function of the position $p$ of the plane in wavelengths from the plane "
                    "that the slice radiated last, for 55 periods of the FLASH1 undulator. In the averaged mode the "
                    "pulse is 55 equal planes. In the unaveraged mode its first two planes are "
                    f"{two(ends[0])} and {two(ends[1])} of the others, and its last two the same.")
    lines += figure("ridge", "The relative frequency of the ridge as a function of x, the square of the transverse "
                    "wavenumber in units of the carrier wavenumber over gamma z, from zero to 0.2. The points of both "
                    "undulators and of both planes lie on the straight line minus x. The dashed curve of the "
                    "resonance of undulator radiation lies below the line, by 0.03 at x of 0.14. The points of the "
                    "spectral propagator lie on the dashed curve.",
                    "The relative frequency $\\nu$ at which the spectrum at the transverse wavenumber $k_\\perp$ is "
                    "largest, as a function of $x$. The solid line is $-x$, which the envelope equation gives, and "
                    "the dashed curve is the resonance of undulator radiation, $\\nu(1 + \\nu) = -x$. The largest "
                    "distance of the points from $-x$, in units of the line width $0.8859/N$, is "
                    f"{two(ridge[('us_flash', 'ridge_x')])} for FLASH1 along $x$, {two(ridge[('us_flash', 'ridge_y')])} "
                    f"along $y$, {two(ridge[('us_delta', 'ridge_x')])} for the Delta undulator and "
                    f"{two(ridge[('us_unavg', 'ridge_x')])} in the unaveraged mode. At the widest wavenumber of FLASH1, "
                    f"$x = {flash_last[0]:.3f}$, the resonance is at $\\nu = {resonance(flash_last[0]):.4f}$ and the run at "
                    f"${flash_last[1] - np.array(f['ridge_x'])[0, 1]:.4f}$, and at that of the Delta undulator, "
                    f"$x = {delta_last[0]:.3f}$, the two are ${resonance(delta_last[0]):.4f}$ and "
                    f"${delta_last[1] - np.array(d['ridge_x'])[0, 1]:.4f}$. With `field.propagator = \"spectral\"`, at 55 periods, the "
                    f"largest distance of the points from the resonance is {two(ridge['spectral'])} of the line width.")
    lines += figure("angular", "The level along the ridge divided by the level on the axis, as a function of x. The "
                    "points of the runs are at 1 for every x. The curves of the formula for FLASH1 fall along x to "
                    "0.34 at x of 0.14 and rise along y to 1.05, and that for the Delta undulator falls to 0.38 at "
                    "x of 0.2.",
                    "The energy per unit frequency and per unit area of transverse wavenumber at the peak of the "
                    "line, $w(k_\\perp)$, divided by its value on the axis, as a function of $x$. The curves are those "
                    "of undulator radiation, $W_1(\\theta)/k^2$ at $\\theta = k_\\perp/k$, in the colors of the points. "
                    "The points are the runs, each divided by the spectrum of the source that the slice deposits in "
                    "one step, which removes the response of the grid to a point. The runs give 1 at every "
                    "wavenumber: the model applies no angular factor. At the widest wavenumber the formula gives "
                    f"{fx[1]:.3f} for FLASH1 along $x$ and {fy[1]:.3f} along $y$, at $x = {fx[0]:.3f}$, and "
                    f"{dx[1]:.3f} for the Delta undulator at $x = {dx[0]:.3f}$. The unaveraged mode gives "
                    f"{ux[2]:.3f} at $x = {ux[0]:.3f}$.")
    lines += figure("filter", "The energy at each angle with the source filter on, divided by that with the filter "
                    "off, as a function of the angle in units of the filter's edge. It is 1 below the edge and "
                    "falls to zero across it.",
                    "The energy radiated at the angle $\\theta$ with `run.source_filter = true`, divided by that of "
                    "the same beam with the filter off, along $x$, as a function of $\\theta/\\theta_f$, where $\\theta = k_\\perp/k_0$ and $\\theta_f$ is "
                    f"the edge that the run derived, {rec['us_filter']['edge'] * 1e6:.1f} µrad, "
                    f"at which $x = {rec['us_filter']['x_edge']:.4f}$. The beam has the emittance of FLASH1. Inside "
                    f"0.8 of the edge the ratio differs from 1 by at most {two(kept)}, and outside 1.2 of the edge "
                    f"it is at most {two(cut)}, along $x$ and along $y$ alike. At the edge the factor $W_1(\\theta)/W_1(0)$ of "
                    f"undulator radiation is {edge_x:.3f} along $x$ and {edge_y:.3f} along $y$.")
    lines += figure("beam", "Two lines of one slice as a function of the relative frequency times the number of "
                    "periods, with the curve of sinc squared. The line of the beam with emittance and the line of "
                    "the beam with energy spread are each displaced from the curve by a fraction of its width, and "
                    "their first zeros are filled in.",
                    "The line on the axis of one slice of 62 electrons in 165 periods of the FLASH1 undulator, with "
                    "the normalized emittance of FLASH1, 1.5 µm, and with its relative energy spread, $6 \\times 10^{-4}$, "
                    "each divided by its largest value, and the curve $\\mathrm{sinc}^2(\\pi N\\nu)$. Each is one "
                    "draw of the shot noise. The peak is displaced by "
                    f"{two(shift['us_size'][0])} of the line width with the emittance and by {two(shift['us_spread'][0])} "
                    f"with the energy spread, and the width is {shift['us_size'][1]:.3f} and {shift['us_spread'][1]:.3f} "
                    "of $0.8859/N$.")
    lines += figure("bunch", "The mean spectrum on the axis of a bunch of 1024 slices as a function of the relative "
                    "frequency times the number of periods. The points scatter about the line of one electron "
                    "multiplied by the number of electrons.",
                    f"The spectrum on the axis of a bunch of {b['slices']} slices, {b['electrons']:.0f} electrons, in 55 "
                    f"periods of the FLASH1 undulator: the mean of {b['seeds']} seeds, divided by $N_e N^2 W_1(0)$. "
                    "The curve is the line of one electron. The energy on the axis, which is the sum over the "
                    f"line, divided by that of one electron times $N_e$, has a mean of {mean:.3f} over the seeds and "
                    f"a standard deviation of {std:.3f}.")
    lines += table(rec)
    pathlib.Path(args.page).write_text("\n".join(lines))
    print(f"wrote {args.page} and 7 figures in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
