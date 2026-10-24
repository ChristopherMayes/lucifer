"""Record what a check measured, beside the line it prints.

A check script prints one line per check: the name, the value it measured, the bound it
held the value to, and ok or FAIL. The benchmark reads only the exit status, so the
measured values reached the documentation by hand. record() writes the same facts as
one line of the results file, and the documentation is generated from them.

    from checkrec import record
    record(name, value, bound, kind, passed)

name    the check's printed name, unique within its section, since the generator pairs
        the debug and production records by it.
value   the number measured, or None when the check has no number (a message matched).
bound   the tolerance or bound the value was held to: a number, or [lo, hi] for a range,
        or None.
kind    how the value met the bound: "max" (value <= bound), "min" (value >= bound),
        "range" (lo <= value <= hi), "equal" (value == bound), or "bool" (no number).
passed  whether the check passed, as the script itself judged it.

The line is written to the file named by LUCIFER_CHECK_RECORD, which run_fel_benchmark.sh
sets for each section, with the section's name from LUCIFER_CHECK_SECTION, as

    check|{"section": ..., "name": ..., "value": ..., "bound": ..., "kind": ..., "pass": ...}

The benchmark appends each section's file to its results file beside the section's
outcome. With LUCIFER_CHECK_RECORD unset, as in a script run by hand, record() does
nothing, so a script behaves the same inside and outside the benchmark.
"""

from __future__ import annotations

import json
import math
import os

KINDS = ("max", "min", "range", "equal", "bool")


def _number(x):
    """A JSON number for a value, None for no value. NaN and infinities become strings,
    since JSON has no spelling for them and a reader must still see them."""
    if x is None or isinstance(x, bool):
        return x
    if isinstance(x, (list, tuple)):
        return [_number(v) for v in x]
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    if math.isnan(f) or math.isinf(f):
        return repr(f)
    if isinstance(x, int):
        return x
    return f


def record(name: str, value, bound=None, kind: str = "max", passed: bool = True) -> None:
    """Append one check record to the file LUCIFER_CHECK_RECORD names, if set."""
    path = os.environ.get("LUCIFER_CHECK_RECORD")
    if not path:
        return
    if kind not in KINDS:
        raise ValueError(f"record kind {kind!r} is not one of {KINDS}")
    row = {"section": os.environ.get("LUCIFER_CHECK_SECTION", ""), "name": str(name),
           "value": _number(value), "bound": _number(bound), "kind": kind, "pass": bool(passed)}
    with open(path, "a") as fh:
        fh.write("check|" + json.dumps(row) + "\n")
