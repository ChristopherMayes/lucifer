#!/usr/bin/env python3
"""
Report which electron rest energy a genesis4 binary was compiled with.

Genesis releases up to v4.6.14 use eev = 510999.06, which is 2.14e-7 above the CODATA
value. Master carries 510998.95069 (CODATA 2022), which is also Bmad's m_electron. The
field unit is linear in eev, so the difference reaches the comparison as 4.28e-7 in
intensity and loosens the transcription tiers by several percent of their own value.

That is invisible in a run: every tier still passes its tolerance. What it moves is the
recorded level, so the reference binary is worth identifying before a level is recorded
against it. The constant is read straight out of the binary, since it is compiled in and
the banner does not carry it.

Usage:

  genesis_constants.py <path to genesis4>

Prints codata, legacy or unknown, and exits zero unless the binary cannot be read.
"""

import struct
import sys
from pathlib import Path

CODATA = 510998.95069      # CODATA 2022, and Bmad's m_electron
LEGACY = 510999.06         # Genesis releases up to v4.6.14


def classify(path):
    try:
        blob = Path(path).read_bytes()
    except OSError as exc:
        return "unknown", f"cannot read {path}: {exc}"
    # Both byte orders, so the answer does not depend on the host.
    def present(value):
        return any(blob.count(struct.pack(fmt, value)) for fmt in ("<d", ">d"))
    has_codata, has_legacy = present(CODATA), present(LEGACY)
    if has_codata and not has_legacy:
        return "codata", f"eev = {CODATA} (CODATA 2022, matches Bmad's m_electron)"
    if has_legacy and not has_codata:
        return "legacy", (f"eev = {LEGACY}, which is {abs(LEGACY - CODATA) / CODATA:.2e} "
                          "above CODATA. Transcription levels read looser against this "
                          "build than against one built from master.")
    if has_codata and has_legacy:
        return "unknown", "both constants appear in the binary"
    return "unknown", "neither constant found, so the layout is not what this expects"


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__.strip().split("\n")[0])
    verdict, detail = classify(sys.argv[1])
    print(f"{verdict}: {detail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
