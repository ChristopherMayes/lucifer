#!/usr/bin/env python3
"""
Check doc/input-reference.md against the struct declarations that define the defaults.

The input reference states a default for every parameter Lucifer owns. Those defaults
live in the Fortran declarations, so a hand-written table can drift from them. This
check holds the two together in both directions:

  1. Every fel-owned parameter the page names must exist in its struct, with the
     default the page states.
  2. Every component of those structs must be named somewhere on the page, so the
     reference cannot go quietly incomplete when a parameter is added.

Bmad-owned structures (beam_init, bmad_com, space_charge_com) are exempt. Their
defaults are Bmad's to state, and the page documents only which fields Lucifer honors.

Usage:

  check_input_reference.py [--doc <file>] [--code <dir>]
"""

import argparse
import re
import sys
from pathlib import Path

# struct name -> (source file, the prefix the page writes it under)
STRUCTS = {
    "fel_global_struct":           ("fel_struct.f90",          "global%"),
    "wavefront_init_struct":       ("fel_struct.f90",          "wavefront_init%"),
    "fel_chamber_wake_init_struct": ("fel_struct.f90",         "chamber_wake%"),
    "fel_beam_init_param_struct":  ("fel_struct.f90",          ""),
    "fel_resample_param_struct":   ("fel_import_mod.f90",      "resample%"),
    "fel_space_charge_struct":     ("fel_collective_mod.f90",  "space_charge%"),
}

DECL = re.compile(r"""^\s*(?:type\s*\([^)]*\)|character\([^)]*\)|real\(rp\)|integer|logical)
                       \s*(?:,\s*allocatable)?\s*::\s*(.*?)\s*(?:!.*)?$""", re.X)


def norm_default(v):
    """Fortran literal -> the form a document would write."""
    v = v.strip().rstrip(",").strip()
    v = re.sub(r"_rp\b", "", v)
    if v in (".false.", ".true."):
        return {".false.": "F", ".true.": "T"}[v]
    # An array literal whose tail is all zeros: the leading value is the default a
    # user writes, since "harmonics = 1" sets the first entry and leaves the rest.
    am = re.match(r"\[\s*([^,\]]+?)\s*(?:,\s*0\s*)*\]$", v)
    if am:
        return am.group(1).strip()
    if v.startswith("'") and v.endswith("'"):
        inner = v[1:-1]
        return '""' if inner == "" else f'"{inner}"'
    return v


def parse_struct(text, name):
    """{component: normalized default} for one type block."""
    m = re.search(r"^type " + name + r"\b(.*?)^end type", text, re.S | re.M)
    if not m:
        raise SystemExit(f"struct {name} not found")
    out = {}
    for line in m.group(1).split("\n"):
        d = DECL.match(line)
        if not d:
            continue
        # One line may declare several: "track_start = '', track_end = ''"
        for part in re.split(r",(?![^(\[]*[)\]])", d.group(1)):
            part = part.strip()
            if not part:
                continue
            mm = re.match(r"(\w+)\s*(?:\([^)]*\))?\s*(?:=\s*(.+))?$", part)
            if not mm:
                continue
            comp, val = mm.group(1), mm.group(2)
            if val is None:
                # a continued array literal from the previous component, not a new one
                if re.match(r"^[\d.\s\[\]'-]+$", part):
                    continue
                out[comp] = None
            else:
                out[comp] = norm_default(val)
    return out


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--doc", default=str(here / "../../doc/input-reference.md"))
    p.add_argument("--code", default=str(here / "../../code"))
    args = p.parse_args()

    doc = Path(args.doc).read_text()
    code = Path(args.code)
    failed = []
    checked = 0

    # The page states parameters in table rows: | `name` | `default` | meaning |
    # Split on the pipes rather than pattern-matching the cells: a cell holding
    # anything unexpected must land in the comparison, not slip past the regex.
    rows = {}
    for line in doc.split("\n"):
        if not line.strip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        name = cells[0]
        if not (name.startswith("`") and name.endswith("`")):
            continue
        rows[name.strip("`").strip()] = cells[1].strip().strip("`").strip()

    for struct, (src, prefix) in STRUCTS.items():
        decl = parse_struct((code / src).read_text(), struct)
        for comp, default in decl.items():
            name = prefix + comp
            checked += 1
            # (1) named on the page at all
            if name not in rows and not re.search(r"`" + re.escape(name) + r"`", doc):
                failed.append(f"{struct}: {name} is not named in the reference")
                continue
            # (2) if it has a table row, the stated default must match
            if name in rows and default is not None:
                stated = rows[name]
                if stated in ("", "-", "(none)"):
                    continue
                if stated != default:
                    failed.append(f"{struct}: {name} default is {default} in code, "
                                  f"{stated} on the page")

    # Parameters the page states under a fel-owned prefix that no struct declares.
    known = {prefix + c for s, (src, prefix) in STRUCTS.items()
             for c in parse_struct((code / src).read_text(), s)}
    for name in rows:
        if "%" not in name:
            continue
        fam = name.split("%")[0] + "%"
        if fam in {p for _, p in STRUCTS.values()} and name not in known:
            failed.append(f"{name} is on the page but no struct declares it")

    for f in failed:
        print(f"  FAIL: {f}")
    print(f"  input reference: {checked} fel-owned parameters checked, "
          f"{len(failed)} problems")
    print("checks: " + ("FAIL" if failed else "PASS"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
