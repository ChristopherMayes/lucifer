#!/usr/bin/env python3
"""
Write the tables of docs/shot-noise.md from the closed form of the shot-noise load.

The tables are arithmetic on shot_noise_ratio.py and take no run. The keystone
regenerates the file and requires no difference, so the page quotes no ratio by hand.

    tab-shot-ratio   what the first-order kick delivers of the shot-noise power against
                     the electrons a beamlet holds, for one harmonic and for the default
                     beamlet's three
    tab-shot-bessel  the same for the kick sized from J1, with the share of the draws held
                     at the largest bunching a kick gives
    tab-shot-floor   the fewest electrons a beamlet may hold for the load to be within a
                     stated fraction of the shot noise, under each kick

Usage: report_shot_noise.py --out docs/generated/shot-noise.md
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import shot_noise_ratio as snr  # noqa: E402

# The default beamlet, beamlet_size = 8, and the harmonics fel_fawley_noise kicks on it.
M = 8
HARMONICS = (1, 2, 3)

# The electron counts the page quotes. 146.4 is the count check_shot_noise.py loads at.
COUNTS = (0.25, 0.5, 1, 2, 4, 10, 18, 50, 100, 146.4, 150, 1000)
DEFICITS = (0.10, 0.05, 0.02, 0.01)


def block(label, caption, header, rows):
    L = [f"<!-- begin {label} -->", f"```{{table}} {caption}", f":name: {label}", ""]
    L.append("| " + " | ".join(header) + " |")
    L.append("|" + "---|" * len(header))
    for r in rows:
        L.append("| " + " | ".join(r) + " |")
    L += ["```", f"<!-- end {label} -->", ""]
    return L


def table_ratio():
    rows = []
    for n in COUNTS:
        first = f"{1 - 2 / n:.3f}" if n >= 4 else ""
        rows.append([f"{n:g}", first, f"{float(snr.ratio_one_harmonic(n, 'first_order')):.3f}"]
                    + [f"{float(snr.ratio(n, h, M, kick='first_order')):.3f}" for h in HARMONICS])
    return block("tab-shot-ratio",
                 "The mean square bunching the first-order kick imposes, over the shot-noise value, against "
                 "the electrons one beamlet holds. The first two columns kick one harmonic. The "
                 f"last three are the default beamlet's, {M} macroparticles with {len(HARMONICS)} harmonics "
                 "kicked together. Below one electron the kick is the one for one electron.",
                 ["electrons a beamlet", "$1 - 2/N$", "one harmonic"]
                 + [f"default beamlet, $h = {h}$" for h in HARMONICS], rows)


def table_bessel():
    rows = []
    for n in COUNTS:
        rows.append([f"{n:g}", f"{np.exp(-snr.B_MAX**2 * n):.3f}", f"{float(snr.ratio_one_harmonic(n)):.3f}"]
                    + [f"{float(snr.ratio(n, h, M)):.3f}" for h in HARMONICS])
    return block("tab-shot-bessel",
                 "The mean square bunching the kick sized from $J_1$ imposes, over the shot-noise value, against "
                 "the electrons one beamlet holds. The second column is the share of the draws that ask for more "
                 "bunching than a kick gives and are held at the largest. The third kicks one harmonic. The "
                 f"last three are the default beamlet's, {M} macroparticles with {len(HARMONICS)} harmonics "
                 "kicked together.",
                 ["electrons a beamlet", "draws held", "one harmonic"]
                 + [f"default beamlet, $h = {h}$" for h in HARMONICS], rows)


def floor(fraction, h, kick):
    """The fewest electrons at and above which the default beamlet's ratio stays within the
    fraction of 1. The last count outside it is found on 100 counts a decade from 1 to
    1e4, and the crossing above it by bisection."""
    n = np.logspace(0, 4, 401)
    out = np.nonzero(np.abs(snr.ratio(n, h, M, kick=kick) - 1) > fraction)[0]
    if not len(out):
        return n[0]
    lo, hi = n[out[-1]], n[out[-1] + 1]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if abs(float(snr.ratio(mid, h, M, kick=kick)) - 1) > fraction:
            lo = mid
        else:
            hi = mid
    return hi


def table_floor():
    rows = [[f"{100 * d:g}%"] + [f"{int(np.ceil(floor(d, h, kick)))}" for kick in ("first_order", "bessel")
                                 for h in HARMONICS] for d in DEFICITS]
    return block("tab-shot-floor",
                 "The fewest electrons a beamlet may hold for the imposed noise power to be "
                 f"within the stated fraction of the shot-noise power, for the default beamlet of {M} macroparticles "
                 "under each kick.",
                 ["within"] + [f"{name}, $h = {h}$" for name in ("first order", "from $J_1$") for h in HARMONICS], rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="The Markdown file to write")
    a = p.parse_args()
    L = ["<!-- Generated by tests/scripts/report_shot_noise.py from shot_noise_ratio.py. "
         "Do not edit. -->", ""]
    L += table_ratio()
    L += table_bessel()
    L += table_floor()
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L))
    print(f"wrote {a.out}: {sum(1 for x in L if x.startswith('<!-- begin'))} tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
