#!/usr/bin/env python3
"""
The bytes the header states for the escaped-field bank, against the files the run kept.

A run whose light leaves the window states before it tracks how many slices leave and,
for each member of the field set, the bytes of the files it keeps and the peak while the
escaped file and the pulse are both on disk (write_bank_bytes in src/fel_io_mod.f90). A
run that keeps the escaped file (run.write_escaped) states one number, the two files
together, which is also its peak. stated() reads those lines and kept() sums the files of
one member on disk.

    python3 bank_bytes.py <log> <dir>

holds every member a log states to its files within TOLERANCE, and a log that states no
bank to a directory with no pulse or escaped file of its out_root. It prints one line per
member and exits nonzero on a failure. The examples job runs it on every example.
"""

from __future__ import annotations

import pathlib
import re
import sys

# The estimate counts the planes, 16 bytes a grid point a component, and the rows of the
# private group. HDF5's own headers and chunk index are the rest, tens of kilobytes. That
# is 1.5e-3 of the two-polarization check's 24 MB, the largest share any check measures,
# and 6e-6 of the sase example's 1.34 GB pulse.
TOLERANCE = 1e-2

OUT_ROOT = re.compile(r'^ Output\s+out_root "([^"]*)"', re.M)
SLICES = re.compile(r"^ Bank\s+(\d+) slices leave the window", re.M)
MEMBER = re.compile(r"^\s+h(\d+): the (pulse file|pulse and escaped files) about (\S+) bytes kept, "
                    r"(?:(\S+) at the peak with the escaped file|which is the peak)", re.M)


def stated(text):
    """The header's bank lines: (slices, {harmonic: (kept, peak, both)}), both True when
    the run keeps the escaped file too. (0, {}) when the header states no bank."""
    m = SLICES.search(text)
    if not m:
        return 0, {}
    members = {}
    for h, what, kept_s, peak_s in MEMBER.findall(text):
        both = what != "pulse file"
        kept_b = float(kept_s)
        members[int(h)] = (kept_b, kept_b if both else float(peak_s), both)
    return int(m.group(1)), members


def out_root(text):
    """The out_root the header names."""
    m = OUT_ROOT.search(text)
    if not m:
        raise ValueError("the log has no Output line naming its out_root")
    return m.group(1)


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
    text = log.read_text(errors="replace")
    root = out_root(text)
    n, members = stated(text)
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
