#!/usr/bin/env python3
"""
Generate doc/generated/validation-measured.md from benchmark results files.

The levels in the documentation are written by the harness that produces them, and not
by hand. Run the benchmark on each build with --results, then point this script at both
files.

The output is deterministic: the same results in give the same bytes out. Nothing time
dependent, machine dependent or path dependent is written, because the keystone suite
regenerates this file and requires the diff to be empty. A timestamp would fail that
check on every run.

Usage:

  report_validation.py --debug <results> --production <results> --out <file.md>

The results files are written by run_fel_benchmark.sh --results, one record per line:

  build|<flavor>
  section|<name>|<outcome>[|<reason>]
  tier|<key>|<level>|<outcome>|<description>
  check|<json>

A check line is one check a section's script made, written by checkrec.record() with
the section, the check's name, the value it measured, the bound it held the value to,
how (max, min, range, equal or bool) and whether it passed. The two builds' records are
paired by section and name. A name recorded twice in one build, or recorded in one
build and not the other, is refused, since a table that silently dropped a check would
say less than it claims.

A results file with a section recorded as anything but pass is refused. The device
section skips on a build with no usable backend, as on a Linux CI runner. A table
written from such a file would name one section fewer while claiming to list them all.
So the generated table comes from a full local run on both builds.
"""

import argparse
import json
import sys
from pathlib import Path


def read_results(path):
    """Parse one results file into {'build':…, 'sections':[…], 'tiers':[…], 'checks':[…]}."""
    out = {"build": None, "sections": [], "tiers": [], "checks": []}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        kind, _, rest = line.partition("|")
        if kind == "check":
            try:
                row = json.loads(rest)
            except json.JSONDecodeError as exc:
                raise SystemExit(f"{path}: a check record is not JSON: {rest[:80]} ({exc})")
            for key in ("section", "name", "value", "bound", "kind", "pass"):
                if key not in row:
                    raise SystemExit(f"{path}: a check record lacks {key!r}: {rest[:80]}")
            out["checks"].append(row)
            continue
        f = rest.split("|")
        if kind == "build":
            out["build"] = f[0]
        elif kind == "section":
            out["sections"].append((f[0], f[1]))
        elif kind == "tier":
            out["tiers"].append({"key": f[0], "level": f[1], "outcome": f[2],
                                 "desc": f[3] if len(f) > 3 else ""})
        else:
            raise SystemExit(f"{path}: unknown record kind {kind!r}")
    if out["build"] is None:
        raise SystemExit(f"{path}: no build record")
    skipped = [n for n, outcome in out["sections"] if outcome != "pass"]
    if skipped:
        raise SystemExit(
            f"{path}: this run skipped {', '.join(skipped)}, so it records fewer sections "
            "than a full\n  run and the table would undercount them. Write the table from a "
            "full run on both builds.")
    seen = set()
    for row in out["checks"]:
        key = (row["section"], row["name"])
        if key in seen:
            raise SystemExit(f"{path}: the check {row['name']!r} is recorded twice in section "
                             f"{row['section']!r}, so the two builds cannot be paired by name")
        seen.add(key)
    return out


def fmt(x):
    """A value as the table prints it: five significant digits for a float, the integer
    itself, yes or no for a truth value, and nothing for no value. Five digits is what a
    check's level reproduces between two runs of one build. The tier levels keep their
    seven, which they have held run after run."""
    if x is None:
        return ""
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.4e}"
    if isinstance(x, list):
        return "[" + ", ".join(fmt(v) for v in x) + "]"
    return str(x)


def fmt_bound(row):
    """The bound a check held its value to, with the sense of the comparison."""
    kind, bound = row["kind"], row["bound"]
    if kind == "max":
        return f"<= {fmt(bound)}"
    if kind == "min":
        return f">= {fmt(bound)}"
    if kind == "range" and isinstance(bound, list) and len(bound) == 2:
        return f"{fmt(bound[0])} to {fmt(bound[1])}"
    if kind == "equal":
        return f"= {fmt(bound)}"
    if kind == "bool":
        return "holds"
    return fmt(bound)


def check_tables(dbg, prod):
    """The per-section tables of checks, debug and production side by side."""
    pmap = {(r["section"], r["name"]): r for r in prod["checks"]}
    dkeys = {(r["section"], r["name"]) for r in dbg["checks"]}
    missing = [k for k in pmap if k not in dkeys] + [k for k in dkeys if k not in pmap]
    if missing:
        raise SystemExit("checks recorded in one build and not the other:\n  "
                         + "\n  ".join(f"{s}: {n}" for s, n in sorted(missing)))
    sections = [s for s, _ in dbg["sections"]]
    by_section = {s: [r for r in dbg["checks"] if r["section"] == s] for s in sections}
    L = []
    n = 0
    for s in sections:
        rows = by_section[s]
        if not rows:
            continue
        # A label per table, so a page can point at one section's levels.
        L.append(f"(val-levels-{s})=")
        L.append(f"### `{s}`")
        L.append("")
        L.append("| Check | Debug | Production | Bound |")
        L.append("|---|---|---|---|")
        for r in rows:
            p = pmap[(s, r["name"])]
            dval, pval = fmt(r["value"]), fmt(p["value"])
            mark = "" if r["pass"] and p["pass"] else " (FAILING)"
            dbound, pbound = fmt_bound(r), fmt_bound(p)
            bound = dbound if dbound == pbound else f"{dbound} (debug), {pbound} (production)"
            name = r["name"].replace("|", "\\|")
            L.append(f"| {name} | {dval}{mark} | {'same' if pval == dval else pval} | {bound} |")
            n += 1
        L.append("")
    return L, n


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--debug", required=True, help="Results file from the debug build")
    p.add_argument("--production", required=True, help="Results file from the production build")
    p.add_argument("--out", required=True, help="Markdown file to write")
    args = p.parse_args()

    dbg = read_results(args.debug)
    prod = read_results(args.production)

    for got, want, path in ((dbg["build"], "debug", args.debug),
                            (prod["build"], "production", args.production)):
        if got != want:
            raise SystemExit(f"{path}: expected a {want} build, got {got!r}")

    dkeys = [t["key"] for t in dbg["tiers"]]
    pkeys = [t["key"] for t in prod["tiers"]]
    if dkeys != pkeys:
        raise SystemExit(f"tier sets differ between builds:\n  debug {dkeys}\n  production {pkeys}")

    plevel = {t["key"]: t["level"] for t in prod["tiers"]}
    poutcome = {t["key"]: t["outcome"] for t in prod["tiers"]}

    L = []
    L.append("<!-- Generated by tests/scripts/report_validation.py. Do not edit. -->")
    L.append("")
    L.append("The tiers are given as measured on the build that ran. A digit here moves only when the")
    L.append("physics moves, and the keystone regenerates this file and requires no diff.")
    L.append("")
    L.append("| Tier | What runs | Debug | Production |")
    L.append("|---|---|---|---|")
    for t in dbg["tiers"]:
        k = t["key"]
        mark = "" if t["outcome"] == "pass" and poutcome[k] == "pass" else " (FAILING)"
        same = t["level"] == plevel[k]
        prod_cell = "same" if same else plevel[k]
        L.append(f"| `{k}` | {t['desc']} | {t['level']}{mark} | {prod_cell} |")
    L.append("")
    L.append('Production levels read "same" where the two builds agree to every printed digit.')
    L.append("Debug and production binaries are never bit-comparable to each other, so the")
    L.append("columns are two independent records rather than a comparison.")
    L.append("")
    L.append("## Check sections")
    L.append("")
    L.append(f"Every section below ran and passed on both builds, {len(dbg['sections'])} in all.")
    L.append("A section that fails stops the run, so an entry here means it passed.")
    L.append("")
    dsec = [s for s, _ in dbg["sections"]]
    psec = [s for s, _ in prod["sections"]]
    if dsec != psec:
        raise SystemExit(f"section sets differ between builds:\n  debug {dsec}\n  production {psec}")
    for s in dsec:
        L.append(f"- `{s}`")
    L.append("")

    tables, n_checks = check_tables(dbg, prod)
    if tables:
        L.append("## Measured levels of the checks")
        L.append("")
        L.append("Each table below is one check section as its script recorded it: the value each check")
        L.append("measured on each build and the bound it held the value to. The scripts write these")
        L.append('records as they run, and the keystone regenerates them. Production reads "same"')
        L.append("where the two builds agree to every printed digit.")
        L.append("")
        L.extend(tables)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(L))
    print(f"wrote {args.out}: {len(dbg['tiers'])} tiers, {len(dsec)} check sections, {n_checks} checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
