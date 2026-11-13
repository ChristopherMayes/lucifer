#!/usr/bin/env python3
"""
Print the estimate-against-the-clock table of docs/performance.md from the stats files of runs.

  estimate_vs_clock.py <name>=<stats file> [<name>=<stats file> ...]

Every run records an estimate of its walk before it tracks and the measured walk after,
as timing/estimate_walk and timing/walk of its stats file, and whether the source filter
was on as params/run/source_filter. This script prints one Markdown row per run: the name
given, the filter state, the measured walk, the estimate and their ratio. The runs are the
eleven comparison tiers and four examples, each run on its own at the machine's thread
count, since a run beside the keystone's other jobs reads a longer walk than the same run
alone.
"""
import sys

import h5py


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    print("| run | filter | measured | estimate | estimate/measured |")
    print("|---|---|---|---|---|")
    bad = 0
    for arg in argv:
        name, _, path = arg.partition("=")
        try:
            with h5py.File(path, "r") as f:
                measured, estimate = float(f["timing/walk"][()]), float(f["timing/estimate_walk"][()])
                filt = "on" if bool(f["params/run/source_filter"][()]) else "off"
        except (OSError, KeyError):
            print(f"| {name} | | no timing/ in {path} | | |")
            bad += 1
            continue
        print(f"| {name} | {filt} | {measured:.3g} s | {estimate:.3g} s | {estimate / measured:.2f} |")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
