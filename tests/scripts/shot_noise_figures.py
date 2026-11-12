#!/usr/bin/env python3
"""
Draw the figures of docs/shot-noise.md.

    beamlet.svg   a beamlet's phases on the unit circle before and after the kick, with
                  the bunching phasor the kick leaves
    kick.svg      J1(z) against z/2, under the distribution of the kick z at 150, 18 and
                  1 electrons a beamlet
    ratio.svg     what the first-order kick delivers of the shot-noise power against the
                  electrons a beamlet holds, with the levels check_shot_noise.py measured
                  on loaded beams
    weight.svg    the fundamental under the weight load and under the first-order kick,
                  with the weight load's measured levels at harmonics 1 to 7
    slices.svg    the noise power of each slice of a loaded Gaussian bunch against its
                  place along the bunch, under each method, beside the real beam's
    spectrum.svg  the spectrum of that bunch from the particles' own places, from zero to
                  past the eighth harmonic
    spectrum-low.svg  the same below 2.2 resonant wavenumbers on a linear axis, with its
                  mean in bins
    bunch.svg     slices.svg over spectrum.svg in one figure, which opens the page

The curves are arithmetic on shot_noise_ratio.py. The measured points are read from the
levels page, docs/generated/validation-measured.md, which the keystone regenerates before
this script runs. The last four are drawn from spectrum.json in the output directory,
the record shot_noise_spectrum.py writes on demand. Every figure is written through
lucifer.plot.save, so the same numbers give the same bytes.

No figure has a title, the axes carry the page's symbols with their units, and the methods
are named by the values of beam.shot_noise_method.

Usage: shot_noise_figures.py --levels docs/generated/validation-measured.md
                             --out docs/generated/shot-noise
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import shot_noise_ratio as snr  # noqa: E402

M = 8
HARMONICS = (1, 2, 3)

# The kick drawn in beamlet.svg: one harmonic, z = 0.7 rad, which is the size of a typical
# kick at eight electrons a beamlet and is drawn large so the displacement is visible.
Z_DRAWN = 0.7
PHI_DRAWN = 0.6


def j1(z):
    """J1 from its integral over one period, summed by the trapezoid rule."""
    t = 2 * np.pi * np.arange(512) / 512
    return np.cos(t - np.multiply.outer(z, np.sin(t))).mean(axis=-1)


def measured(levels):
    """The measured levels at the fundamental and the two harmonics, uniform weights.

    Returns {harmonic: [(electrons, value, sigma), ...]}. A row's range is the closed
    form with five standard deviations either side."""
    row = re.compile(r"^\| shot noise, uniform weights(?:, ([0-9.]+) electrons a beamlet)?, "
                     r"harmonic ([123]) \| (\S+) \| (\S+) \| (\S+) to (\S+) \|")
    out = {h: [] for h in HARMONICS}
    for line in pathlib.Path(levels).read_text().splitlines():
        m = row.match(line)
        if not m:
            continue
        electrons = float(m.group(1)) if m.group(1) else 146.4
        value = float(m.group(3) if m.group(4) == "same" else m.group(4))
        sigma = (float(m.group(6)) - float(m.group(5))) / 10
        out[int(m.group(2))].append((electrons, value, sigma))
    if not all(out.values()):
        raise SystemExit(f"{levels}: no shot-noise rows for the uniform weights. Run the validation "
                         "suite and report_validation.py first.")
    return out


def draw_beamlet(plt, lp, path):
    theta = 2 * np.pi * np.arange(M) / M + 0.25
    kicked = theta - Z_DRAWN * np.sin(theta + PHI_DRAWN)
    b = np.exp(-1j * kicked).mean()
    circle = np.linspace(0, 2 * np.pi, 361)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7), layout="constrained")
    for ax, ph, title in ((axes[0], theta, "quiet"), (axes[1], kicked, "after the kick")):
        ax.plot(np.cos(circle), np.sin(circle), color=lp.GRAY, lw=1)
        ax.plot([-1.15, 1.15], [0, 0], color=lp.GRAY, lw=0.6)
        ax.plot([0, 0], [-1.15, 1.15], color=lp.GRAY, lw=0.6)
        ax.set_aspect("equal")
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.axis("off")
        ax.set_title(title, fontsize=10)
        ax.plot(np.cos(ph), -np.sin(ph), "o", color=lp.COLORS[0], ms=7, zorder=3)
    ax = axes[1]
    ax.plot(np.cos(theta), -np.sin(theta), "o", mfc="none", mec=lp.GRAY, ms=7, zorder=2)
    for t0, t1 in zip(theta, kicked):
        if abs(t1 - t0) < 0.12:
            continue
        arc = np.linspace(t0, t1 - 0.06 * np.sign(t1 - t0), 30)
        ax.plot(1.13 * np.cos(arc[:-1]), -1.13 * np.sin(arc[:-1]), color=lp.INK2, lw=0.9)
        ax.annotate("", xy=(1.13 * np.cos(arc[-1]), -1.13 * np.sin(arc[-1])),
                    xytext=(1.13 * np.cos(arc[-3]), -1.13 * np.sin(arc[-3])),
                    arrowprops=dict(arrowstyle="-|>", color=lp.INK2, lw=0.9, shrinkA=0, shrinkB=0))
    ax.annotate("", xy=(b.real, b.imag), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=lp.COLORS[1], lw=2, shrinkA=0, shrinkB=0))
    ax.text(b.real + 0.06, b.imag + 0.1, r"$|b| = J_1(z)$", color=lp.COLORS[1],
            ha="center", va="bottom")
    axes[0].plot(0, 0, "o", color=lp.COLORS[1], ms=4)
    axes[0].text(0.08, 0.08, r"$b = 0$", color=lp.COLORS[1])
    lp.save(fig, path)
    plt.close(fig)
    return abs(b)


def draw_kick(plt, lp, path):
    z = np.linspace(0, 4, 801)
    fig, ax = plt.subplots(figsize=(6.4, 3.8), layout="constrained")
    ax.plot(z, z / 2, color=lp.INK2, ls="--", lw=1.2, label=r"$z/2$, the first-order kick")
    ax.plot(z, j1(z), color=lp.COLORS[0], lw=2, label=r"$J_1(z)$, what a kick gives")
    for n, color in ((150, lp.COLORS[2]), (18, lp.COLORS[3]), (1, lp.COLORS[4])):
        p = n * z / 2 * np.exp(-n * z**2 / 4)
        ax.fill_between(z, 0, 0.55 * p / p.max(), color=color, alpha=0.2, lw=0,
                        label=f"$N = {n}$")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 0.75)
    ax.set_xlabel("$z$ (rad)")
    ax.set_ylabel("$|b_1|$")
    ax.legend(frameon=False, loc="upper right")
    lp.save(fig, path)
    plt.close(fig)


def draw_ratio(plt, lp, path, points):
    n = np.logspace(-1, 3, 801)
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    ax.axvspan(0.1, 1, color=lp.GRAY, alpha=0.18, lw=0)
    big = n >= 2
    ax.plot(n[big], 1 - 2 / n[big], color=lp.INK2, ls=":", lw=1.2, label=r"$1 - 2/N$")
    ax.plot(n, snr.ratio_one_harmonic(n), color=lp.INK2, ls="--", lw=1.2, label="one harmonic")
    for h, color in zip(HARMONICS, lp.COLORS):
        ax.plot(n, snr.ratio(n, h, M), color=color, lw=1.8, label=f"$h = {h}$")
        e, v, s = np.array(points[h]).T
        ax.errorbar(e, v, yerr=s, fmt="o", color=color, ms=4, lw=1, capsize=2, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("$N$")
    ax.set_ylabel("$R_h(N)$")
    ax.legend(frameon=False, loc="lower right")
    lp.save(fig, path)
    plt.close(fig)


def measured_weight(levels):
    """The weight load's measured levels: {harmonic: [(electrons, value, sigma), ...]}."""
    row = re.compile(r"^\| shot noise, weight load, uniform weights, ([0-9.]+) electrons? a beamlet, "
                     r"harmonic ([1-7]) \| (\S+) \| (\S+) \| (\S+) to (\S+) \|")
    out = {}
    for line in pathlib.Path(levels).read_text().splitlines():
        m = row.match(line)
        if not m:
            continue
        value = float(m.group(3) if m.group(4) == "same" else m.group(4))
        out.setdefault(int(m.group(2)), []).append((float(m.group(1)), value, (float(m.group(6)) - float(m.group(5))) / 10))
    if not out:
        raise SystemExit(f"{levels}: no shot-noise rows for the weight load. Run the validation suite and "
                         "report_validation.py first.")
    return out


def draw_weight(plt, lp, path, points):
    n = np.logspace(-1, 3, 401)
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    ax.plot(n, snr.ratio(n, 1, M), color=lp.COLORS[2], ls="--", lw=1.4, label='"first_order"')
    ax.plot(n, np.ones_like(n), color=lp.COLORS[1], lw=1.8, label='"weight"')
    for h, pts in sorted(points.items()):
        e, v, s = np.array(pts).T
        ax.errorbar(e * (1 + 0.035 * (h - 4)), v, yerr=s, fmt="o", color=lp.COLORS[1], ms=3, lw=0.8, capsize=1.5,
                    zorder=4, label='"weight", measured at $h = 1$ to 7' if h == 1 else None)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$R_1(N)$")
    ax.legend(frameon=False, loc="lower right")
    lp.save(fig, path)
    plt.close(fig)


# The methods of the record, in drawing order: the key in the record, the value of
# beam.shot_noise_method as the input file takes it, the color and the marker.
def _methods(lp):
    return (("first_order", '"first_order"', lp.COLORS[2], "^"), ("weight", '"weight"', lp.COLORS[1], "o"))


# The load with neither the quiet start nor the shot noise, which has the noise of its
# macroparticles. The gray curve of the spectrum is the quiet start with no shot noise.
def _LOUD(lp):
    return ("loud", "quiet_start = false", lp.COLORS[3], "s")


def _above(ax, ncol):
    """The legend above the axes, where no curve crosses it."""
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=ncol)


def _slices(ax, lp, rec):
    """The noise power of each slice at the fundamental against its place along the bunch."""
    n_slice, spacing = rec["n_slice"], rec["n_wavelength"] * rec["lambda0"]
    z = (np.arange(n_slice) - (n_slice - 1) / 2) * spacing * 1e9
    peak = rec["peak_electrons_a_beamlet"] * rec["n_particle"] / rec["beamlet_size"]
    zz = np.linspace(z[0], z[-1], 600)
    ax.plot(zz, peak * np.exp(-(zz * 1e-9) ** 2 / (2 * rec["sig_z"] ** 2)), color=lp.INK, lw=1.2, label=r"$N_\lambda$")
    for key, label, color, mark in (_LOUD(lp),) + _methods(lp):
        s = np.array(rec["methods"][key]["slice_power"])
        m = s > 0
        ax.plot(z[m], s[m], mark, color=color, ms=2.6, mew=0, label=label)
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 1e5)
    ax.set_xlim(z[0], z[-1])
    ax.set_xlabel("$z$ (nm)")
    ax.set_ylabel(r"$\langle |S_1|^2 \rangle$")


def _spectrum(ax, lp, rec, k_max=None, binned=False):
    """The spectrum of the whole bunch, on a logarithmic axis, or below k_max on a linear one."""
    k = np.array(rec["k"])
    m = np.ones(len(k), bool) if k_max is None else k <= k_max
    if not binned:
        ax.plot(k[m], np.array(rec["methods"]["quiet"]["spectrum"])[m], color=lp.GRAY, lw=0.5,
                label="shot_noise = false")
        key, label, color, _ = _LOUD(lp)
        ax.plot(k[m], np.array(rec["methods"][key]["spectrum"])[m], color=color, lw=0.5, label=label)
    for key, label, color, mark in _methods(lp):
        sp = np.array(rec["methods"][key]["spectrum"])
        if binned:
            ax.plot(k[m], sp[m], color=color, lw=0.4, alpha=0.3)
            edges = np.arange(0, k_max + 1e-9, 0.05)
            mid = 0.5 * (edges[1:] + edges[:-1])
            mean = np.array([sp[(k >= a) & (k < b)].mean() for a, b in zip(edges[:-1], edges[1:])])
            ax.plot(mid, mean, mark + "-", color=color, ms=3.5, lw=1.2, mew=0, label=label)
        else:
            ax.plot(k[m], sp[m], color=color, lw=0.5, label=label)
    ax.axhline(1, color=lp.INK, lw=1.0)
    if binned:
        ax.set_xlim(0, k_max)
        ax.set_ylim(0, 2.0)
    else:
        ax.set_yscale("log")
        ax.set_xlim(0, k[-1])
        ax.set_ylim(1e-6, 1e2)
    ax.set_xlabel("$k / k_r$")
    ax.set_ylabel(r"$\langle |S(k)|^2 \rangle / N$")


def draw_bunch(plt, lp, out, rec):
    """slices.svg, spectrum.svg, spectrum-low.svg and bunch.svg from the record."""
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    _slices(ax, lp, rec)
    _above(ax, 4)
    lp.save(fig, out / "slices.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.4), layout="constrained")
    _spectrum(ax, lp, rec)
    _above(ax, 4)
    lp.save(fig, out / "spectrum.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.4), layout="constrained")
    _spectrum(ax, lp, rec, k_max=2.2, binned=True)
    _above(ax, 2)
    lp.save(fig, out / "spectrum-low.svg")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(6.8, 8.2), layout="constrained")
    _slices(axes[0], lp, rec)
    _above(axes[0], 4)
    _spectrum(axes[1], lp, rec)
    _above(axes[1], 4)
    lp.save(fig, out / "bunch.svg")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--levels", required=True, help="docs/generated/validation-measured.md")
    p.add_argument("--out", required=True, help="The directory the figures are written to")
    a = p.parse_args()
    points = measured(a.levels)
    points_weight = measured_weight(a.levels)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import lucifer.plot as lp

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    record = out / "spectrum.json"
    if not record.is_file():
        raise SystemExit(f"{record}: no record of the loaded bunch. Run shot_noise_spectrum.py first.")
    rec = json.loads(record.read_text())
    with lp._style():
        b = draw_beamlet(plt, lp, out / "beamlet.svg")
        draw_kick(plt, lp, out / "kick.svg")
        draw_ratio(plt, lp, out / "ratio.svg", points)
        draw_weight(plt, lp, out / "weight.svg", points_weight)
        draw_bunch(plt, lp, out, rec)
    print(f"wrote {out}: 8 figures, the drawn kick of {Z_DRAWN:g} rad leaves |b| = {b:.4f}, "
          f"J1 = {float(j1(Z_DRAWN)):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
