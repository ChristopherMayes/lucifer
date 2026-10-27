#!/usr/bin/env python3
"""
Check the saturation demo's answers before any timing or figure is read. The exit total
powers must agree at the documented levels (this directory's README, "The measured
result").

The averaged mode against Genesis carries the seam-transport difference, about 4e-2 on
the Benchmark1-SASE line and 4.9e-4 at saturation, where the power self-limits. The check
holds that difference to a relative 0.15. The unaveraged mode also carries the
shot-noise radiation channel that it resolves and the averaged model does not track.
That channel puts the unaveraged exit power +0.6 in ln above Genesis, independent of
the particle count and of the steps per period. The check holds |ln| <= 1.0, so an
order-of-magnitude disagreement still fails the demo.

Usage: check_agreement.py <genesis .out.h5> <avg diag> <unavg diag>
"""

import math
import sys

import h5py
import numpy as np


def total_exit(fn, nslice):
    d = np.loadtxt(fn)
    return d.reshape(-1, nslice, d.shape[1])[-1, :, 2].sum()


def main():
    gen, avg, unavg = sys.argv[1:4]
    with h5py.File(gen) as h5:
        p_gen = h5["Field/power"][-1, :].sum()
        nslice = h5["Field/power"].shape[1]
    p_avg = total_exit(avg, nslice)
    p_unavg = total_exit(unavg, nslice)
    rel = abs(p_avg - p_gen) / p_gen
    lnr = abs(math.log(p_unavg / p_gen))
    print(f"check: exit total power  Genesis {p_gen:.4e} W")
    print(f"check: Bmad averaged     {p_avg:.4e} W, rel {rel:.2e} (must be <= 0.15)")
    print(f"check: Bmad unaveraged   {p_unavg:.4e} W, |ln ratio| {lnr:.3f} (must be <= 1.0)")
    if rel > 0.15 or lnr > 1.0:
        print("FAIL: disagreement beyond the documented levels; timings are meaningless.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
