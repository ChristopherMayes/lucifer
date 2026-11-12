#!/usr/bin/env python3
"""
Write the tables of docs/shot-noise.md from the closed form of the shot-noise load.

The tables are arithmetic on shot_noise_ratio.py and take no run. The keystone
regenerates the file and requires no difference, so the page quotes no ratio by hand.

    tab-shot-methods  the fundamental under the two values of beam.shot_noise_method
    tab-shot-ratio    what the first-order kick delivers of the shot-noise power against the
                      electrons a beamlet holds, for one harmonic and for the default
                      beamlet's three
    tab-shot-floor    the fewest electrons a beamlet may hold under that kick for a stated
                      deficit

Usage: report_shot_noise.py --out docs/generated/shot-noise.md
"""

from __future__ import annotations

import argparse
import json
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
        rows.append([f"{n:g}", first, f"{float(snr.ratio_one_harmonic(n)):.3f}"]
                    + [f"{float(snr.ratio(n, h, M)):.3f}" for h in HARMONICS])
    return block("tab-shot-ratio",
                 "The ratio of the mean square bunching imposed by the first-order kick to the shot-noise value, "
                 "as a function of the number of electrons per beamlet. The first two columns are for a kick at a "
                 f"single harmonic. The last three are for the default beamlet of {M} macroparticles, with kicks at "
                 f"{len(HARMONICS)} harmonics together. Below one electron the kick is that for one electron.",
                 ["electrons per beamlet", "$1 - 2/N$", "one harmonic"]
                 + [f"default beamlet, $h = {h}$" for h in HARMONICS], rows)


def table_methods():
    rows = []
    for n in COUNTS:
        rows.append([f"{n:g}", "1.000", f"{float(snr.ratio_weight_kept(n)):.3f}", f"{float(snr.ratio(n, 1, M)):.3f}"])
    return block("tab-shot-methods",
                 "The ratio of the noise power of the load at the fundamental to the shot-noise power, as a function "
                 "of the number of electrons per beamlet, for the two values of `beam.shot_noise_method`. The columns "
                 "for the weight load apply at every harmonic below the beamlet size. In the keep mode a particle "
                 "of less than one electron is represented by one electron.",
                 ["electrons per beamlet", '`"weight"`, sample mode', '`"weight"`, keep mode', '`"first_order"`'],
                 rows)


def floor(deficit, h):
    """The fewest electrons at and above which the default beamlet's ratio stays within the deficit."""
    n = np.logspace(0, 4, 40001)
    low = np.nonzero(snr.ratio(n, h, M) < 1 - deficit)[0]
    return n[low[-1] + 1] if len(low) else n[0]


def table_floor():
    rows = [[f"{100 * d:g}%"] + [f"{int(np.ceil(floor(d, h)))}" for h in HARMONICS] for d in DEFICITS]
    return block("tab-shot-floor",
                 "The smallest number of electrons per beamlet for which the noise power imposed by the first-order "
                 f"kick is within the stated fraction of the shot-noise power, for the default beamlet of {M} "
                 "macroparticles.",
                 ["deficit"] + [f"$h = {h}$" for h in HARMONICS], rows)


def table_bands(record):
    """The mean spectrum in the record's bands, from shot_noise_spectrum.py's committed record."""
    names = (("weight", '`"weight"`'), ("first_order", '`"first_order"`'), ("quiet", "`shot_noise = false`"),
             ("loud", "`quiet_start = false`, `shot_noise = false`"))
    rows = []
    for i, (lo, hi) in enumerate(record["bands"]):
        rows.append([f"{lo:g} to {hi:g}"] + [f"{record['methods'][key]['bands'][i]:.3g}" for key, _ in names])
    return block("tab-shot-bands",
                 f"The mean of $|S(k)|^2/N$ for the bunch of [](#fig-shot-bunch) in six bands of the wavenumber, "
                 f"in units of the resonant wavenumber, with {record['seeds']} seeds. The value for the real beam "
                 f"is 1. The script `shot_noise_spectrum.py` requires the weight load to be within "
                 f"{record['bound']:g} of 1 in every band.",
                 ["$k/k_r$"] + [label for _, label in names], rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="The Markdown file to write")
    p.add_argument("--record", help="shot_noise_spectrum.py's record, by default shot-noise/spectrum.json beside --out")
    a = p.parse_args()
    out = pathlib.Path(a.out)
    record = pathlib.Path(a.record) if a.record else out.parent / "shot-noise" / "spectrum.json"
    if not record.is_file():
        print(f"report_shot_noise.py: no record at {record}. Run tests/scripts/shot_noise_spectrum.py.", file=sys.stderr)
        return 1
    L = ["<!-- Generated by tests/scripts/report_shot_noise.py from shot_noise_ratio.py. "
         "Do not edit. -->", ""]
    L += table_bands(json.loads(record.read_text()))
    L += table_methods()
    L += table_ratio()
    L += table_floor()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L))
    print(f"wrote {a.out}: {sum(1 for x in L if x.startswith('<!-- begin'))} tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
