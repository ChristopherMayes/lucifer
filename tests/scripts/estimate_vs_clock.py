#!/usr/bin/env python3
"""
Print the estimate-against-the-clock table of docs/performance.md from run logs.

  estimate_vs_clock.py <name>=<log> [<name>=<log> ...]

Every run prints an estimated walk in its header and the measured walk beside it in its
footer, on one line:

   Cost        walk 11.4 s measured, 8.2 s estimated on the CPU path for <machine>.

This script reads that line from each log, and whether the source filter was on, from
the header's load block, and prints one Markdown row per run: the name given, the filter
state, the measured walk, the estimate and their ratio. The runs are the eleven comparison
tiers and four examples, each run on its own at the machine's thread count, since a run
beside the keystone's other jobs reads a longer walk than the same run alone.
"""
import re
import sys

COST = re.compile(r"^ Cost\s+walk\s+([\d.]+) s measured, ([\d.]+) s estimated", re.M)
OFF = re.compile(r"The source filter is off")


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    print("| run | filter | measured | estimate | estimate/measured |")
    print("|---|---|---|---|---|")
    bad = 0
    for arg in argv:
        name, _, path = arg.partition("=")
        txt = open(path, errors="ignore").read()
        m = COST.search(txt)
        if not m:
            print(f"| {name} | | no Cost line in {path} | | |")
            bad += 1
            continue
        measured, estimate = float(m.group(1)), float(m.group(2))
        filt = "off" if OFF.search(txt) else "on"
        print(f"| {name} | {filt} | {measured:g} s | {estimate:g} s | {estimate / measured:.2f} |")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
