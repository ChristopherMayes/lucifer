#!/usr/bin/env python3
"""
The bytes a run states for the escaped-field bank before it tracks, against the files it kept.

A run whose light leaves the window records in its stats file, before it tracks, how many
slices leave (setup/n_slice_leaving) and, for each member of the field set, the bytes of
the files it keeps (setup/bank_bytes_h<h>): the pulse, and with run.write_escaped the
escaped file too. The banked planes stream into those files and into nothing else, so the
bytes kept are also the peak on disk. stated() reads those values and kept() sums the
files of one member on disk.

    python3 bank_bytes.py <stats file> <dir>

holds every member a stats file states to its files within TOLERANCE, and a stats file
that states no bank to a directory with no pulse or escaped file of its out_root. It
prints one line per member and exits nonzero on a failure. The examples job runs it on
every example.
"""

from __future__ import annotations

import pathlib
import re
import sys

import h5py

# The estimate counts the planes, 16 bytes a grid point a component, and the rows of the
# private group. HDF5's own headers and chunk index are the rest, tens of kilobytes. That
# is 2.5e-3 of the two-polarization check's 24 MB, the largest share any check measures,
# and 6.6e-5 of the sase example's 1.34 GB pulse.
TOLERANCE = 1e-2

MEMBER = re.compile(r"^bank_bytes_h(\d+)$")


def stated(stats_file):
    """What a stats file states of the bank: (slices, {harmonic: (kept, peak, both)}),
    both True when the run keeps the escaped file too. (0, {}) when it states no bank."""
    members = {}
    with h5py.File(stats_file, "r") as f:
        setup = f["setup"]
        both = bool(f["params/run/write_escaped"][()])
        for name in setup:
            m = MEMBER.match(name)
            if m:
                kept_b = float(setup[name][()])
                members[int(m.group(1))] = (kept_b, kept_b, both)
        n = int(setup["n_slice_leaving"][()]) if members else 0
    return n, members


def bank_file(wd, root, harm, what):
    """The pulse or escaped file of one member, as fel_bank_file_name names it."""
    suffix = "" if harm == 1 else f"-h{harm}"
    return pathlib.Path(wd) / f"{root}-{what}{suffix}.wf.h5"


def kept(wd, root, harm):
    """The bytes of the pulse and the escaped file of one member on disk, each 0 when
    the file is not there."""
    return tuple(bank_file(wd, root, harm, w).stat().st_size if bank_file(wd, root, harm, w).exists()
                 else 0 for w in ("pulse", "escaped"))


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    log, wd = pathlib.Path(argv[1]), pathlib.Path(argv[2])
    root = log.name.removesuffix(".stats.h5")
    n, members = stated(log)
    failed = False
    if not members:
        left = sorted(p.name for w in ("pulse", "escaped") for p in wd.glob(f"{root}-{w}*.wf.h5"))
        ok = not left
        failed = not ok
        print(f"bank_bytes: {log.name} states no bank, files {left or 'none'}  {'ok' if ok else 'FAIL'}")
    for h, (kept_b, _, both) in sorted(members.items()):
        pulse_b, esc_b = kept(wd, root, h)
        on_disk = pulse_b + esc_b
        rel = abs(kept_b / on_disk - 1) if on_disk else float("inf")
        ok = rel <= TOLERANCE and (esc_b > 0) == both
        failed |= not ok
        print(f"bank_bytes: {log.name} h{h}, {n} slices: stated {kept_b:.4e} bytes kept, "
              f"{on_disk:.6e} on disk, rel {rel:.1e} (check {TOLERANCE:.0e})  {'ok' if ok else 'FAIL'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
