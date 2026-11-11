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
    bessel.svg    the same for the kick sized from J1, beside the first-order kick's
                  fundamental

The curves are arithmetic on shot_noise_ratio.py. The measured points are read from the
levels page, docs/generated/validation-measured.md, which the keystone regenerates before
this script runs. Every figure is written through lucifer.plot.save, so the same numbers
give the same bytes.

Usage: shot_noise_figures.py --levels docs/generated/validation-measured.md
                             --out docs/generated/shot-noise
"""

from __future__ import annotations

import argparse
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


def measured(levels, kick="first_order"):
    """The measured levels at the fundamental and the two harmonics, uniform weights, under
    the given kick.

    Returns {harmonic: [(electrons, value, sigma), ...]}. A row's range is the closed
    form with five standard deviations either side."""
    lead = "kick from J1, " if kick == "bessel" else ""
    row = re.compile(r"^\| shot noise, " + lead + r"uniform weights(?:, ([0-9.]+) electrons? a beamlet)?, "
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
        raise SystemExit(f"{levels}: no shot-noise rows for the uniform weights, kick {kick}. Run the validation "
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
                        label=f"the kicks drawn at $N = {n}$")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 0.75)
    ax.set_xlabel(r"kick $z = 2\sqrt{x/N}$ (rad)")
    ax.set_ylabel(r"bunching $|b|$ of the beamlet")
    ax.legend(frameon=False, loc="upper right")
    lp.save(fig, path)
    plt.close(fig)


def draw_ratio(plt, lp, path, points):
    n = np.logspace(-1, 3, 801)
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    ax.axvspan(0.1, 1, color=lp.GRAY, alpha=0.18, lw=0)
    ax.text(0.32, 1.04, "kick drawn for\none electron", ha="center", va="top", fontsize=9, color=lp.INK2)
    big = n >= 2
    ax.plot(n[big], 1 - 2 / n[big], color=lp.INK2, ls=":", lw=1.2, label=r"$1 - 2/N$")
    ax.plot(n, snr.ratio_one_harmonic(n, "first_order"), color=lp.INK2, ls="--", lw=1.2, label="one harmonic kicked")
    for h, color in zip(HARMONICS, lp.COLORS):
        ax.plot(n, snr.ratio(n, h, M, kick="first_order"), color=color, lw=1.8, label=f"default beamlet, $h = {h}$")
        e, v, s = np.array(points[h]).T
        ax.errorbar(e, v, yerr=s, fmt="o", color=color, ms=4, lw=1, capsize=2, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("electrons one beamlet holds, $N$")
    ax.set_ylabel(r"$\langle|b_h|^2\rangle\,N$")
    ax.legend(frameon=False, loc="lower right")
    lp.save(fig, path)
    plt.close(fig)


def draw_bessel(plt, lp, path, points):
    n = np.logspace(-1, 3, 401)
    fig, ax = plt.subplots(figsize=(6.8, 4.2), layout="constrained")
    ax.axhline(1, color=lp.GRAY, lw=0.8)
    ax.plot(n, snr.ratio(n, 1, M, kick="first_order"), color=lp.INK2, ls="--", lw=1.4,
            label="first-order kick, $h = 1$")
    ax.plot(n, snr.ratio_one_harmonic(n), color=lp.INK2, ls=":", lw=1.2,
            label=r"kick from $J_1$, one harmonic kicked")
    for h, color in zip(HARMONICS, lp.COLORS):
        ax.plot(n, snr.ratio(n, h, M), color=color, lw=1.8, label=f"kick from $J_1$, $h = {h}$")
        e, v, s = np.array(points[h]).T
        ax.errorbar(e, v, yerr=s, fmt="o", color=color, ms=4, lw=1, capsize=2, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("electrons one beamlet holds, $N$")
    ax.set_ylabel(r"$\langle|b_h|^2\rangle\,N$")
    ax.legend(frameon=False, loc="lower right")
    lp.save(fig, path)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--levels", required=True, help="docs/generated/validation-measured.md")
    p.add_argument("--out", required=True, help="The directory the figures are written to")
    a = p.parse_args()
    points = measured(a.levels)
    points_bessel = measured(a.levels, "bessel")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import lucifer.plot as lp

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with lp._style():
        b = draw_beamlet(plt, lp, out / "beamlet.svg")
        draw_kick(plt, lp, out / "kick.svg")
        draw_ratio(plt, lp, out / "ratio.svg", points)
        draw_bessel(plt, lp, out / "bessel.svg", points_bessel)
    print(f"wrote {out}: 4 figures, the drawn kick of {Z_DRAWN:g} rad leaves |b| = {b:.4f}, "
          f"J1 = {float(j1(Z_DRAWN)):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
