#!/usr/bin/env python3
"""
Check that an edit changed comments and nothing else, and count what the comments say.

  check_comments_only.py --base <rev> <file>...   # code against <rev>, comments ignored
  check_comments_only.py --counts <file>...       # comment counts, one row per file

The first form strips the comments from each file as it is at <rev> and as it is in
the working tree, drops blank lines and trailing blanks, and requires the two to be
identical, so no statement, declaration or continuation can move. Comments are what the
language says they are:

  Fortran   a line whose first nonblank character is !, except an OpenMP sentinel (!$),
            and the rest of a line from an ! outside a quoted string. A string
            continued onto the next line with & carries its quote state there.
  C         /* ... */ across lines, and // to the end of a line, outside strings.
  CMake     # to the end of a line, outside a quoted string.

The second form counts, per file and in total: lines, comment lines, comment lines
that narrate (a date, the user, a keystone, a tier, a check section, a mutation, a
commit, a price, a measurement), comment lines holding a semicolon, and lines over
110 columns that carry a comment. A semicolon in a Fortran statement and a long
statement are code and are not counted.

Exit status is nonzero when a file's code differs from <rev>.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

WIDTH = 110

NARRATION = {
    "date": r"\b20\d\d-\d\d(-\d\d)?\b",
    "the user": r"\bthe user\b",
    "keystone": r"\bkeystone\b",
    "tier": r"\btiers?\b",
    "check section": r"\bcheck sections?\b",
    "mutation": r"\bmutat",
    "commit": r"\bcommit(s|ted)?\b",
    "price": r"\bpric(e|ed|es|ing)\b",
    "measured": r"\bmeasur",
}


def kind(path: str) -> str:
    p = pathlib.Path(path)
    if p.suffix.lower() in (".f90", ".f", ".f95", ".f03"):
        return "fortran"
    if p.suffix.lower() in (".c", ".h", ".cpp", ".mm"):
        return "c"
    return "cmake"


def split_fortran(text: str):
    """(code lines, comment strings) of a Fortran source."""
    code, comments = [], []
    quote = None
    for line in text.splitlines():
        stripped = line.lstrip()
        if quote is None and stripped.startswith("!") and not stripped.startswith("!$"):
            comments.append(stripped[1:])
            continue
        out = []
        i = 0
        while i < len(line):
            ch = line[i]
            if quote:
                out.append(ch)
                if ch == quote:
                    quote = None
            elif ch in ("'", '"'):
                quote = ch
                out.append(ch)
            elif ch == "!" and not line[i:].startswith("!$"):
                comments.append(line[i + 1:])
                break
            else:
                out.append(ch)
            i += 1
        # A string that runs off the end of a line continues only when the line ends in &.
        if quote and not "".join(out).rstrip().endswith("&"):
            quote = None
        code.append("".join(out))
    return code, comments


def split_c(text: str):
    code, comments = [], []
    out, i, quote, n = [], 0, None, len(text)
    while i < n:
        ch = text[i]
        if quote:
            out.append(ch)
            if ch == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
            out.append(ch)
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            comments.extend(text[i + 2:j - 2].splitlines() or [""])
            out.extend("\n" * text[i:j].count("\n"))
            i = j
            continue
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            comments.append(text[i + 2:j])
            i = j
            continue
        else:
            out.append(ch)
        i += 1
    return "".join(out).splitlines(), comments


def split_cmake(text: str):
    code, comments = [], []
    for line in text.splitlines():
        out, quote = [], False
        for i, ch in enumerate(line):
            if ch == '"':
                quote = not quote
            if ch == "#" and not quote:
                comments.append(line[i + 1:])
                break
            out.append(ch)
        code.append("".join(out))
    return code, comments


def split(path: str, text: str):
    return {"fortran": split_fortran, "c": split_c, "cmake": split_cmake}[kind(path)](text)


def normalized(path: str, text: str):
    code, _ = split(path, text)
    return [line.rstrip() for line in code if line.strip()]


def at_rev(rev: str, path: str) -> str:
    top = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True,
                         text=True, check=True).stdout.strip()
    rel = str(pathlib.Path(path).resolve().relative_to(top))
    return subprocess.run(["git", "-C", top, "show", f"{rev}:{rel}"], capture_output=True,
                          text=True, check=True).stdout


def compare(rev: str, paths) -> int:
    bad = 0
    for path in paths:
        before = normalized(path, at_rev(rev, path))
        after = normalized(path, pathlib.Path(path).read_text())
        if before == after:
            print(f"ok    {path}: code identical to {rev}, {len(after)} code lines")
            continue
        bad += 1
        print(f"FAIL  {path}: code differs from {rev}")
        for k, (a, b) in enumerate(zip(before, after)):
            if a != b:
                print(f"      first difference at code line {k + 1}:\n        was: {a}\n        now: {b}")
                break
        else:
            print(f"      code line count {len(before)} became {len(after)}")
    return bad


def counts(paths) -> None:
    keys = list(NARRATION)
    head = f"{'file':40s} {'lines':>6s} {'comm':>5s} {'%':>3s} " + " ".join(
        f"{k[:5]:>5s}" for k in keys) + f" {'semi':>5s} {'>110':>5s}"
    print(head)
    tot = dict.fromkeys(["lines", "comm", "semi", "wide", *keys], 0)
    for path in paths:
        text = pathlib.Path(path).read_text()
        lines = text.splitlines()
        _, comments = split(path, text)
        # A long statement is code. Only a line that carries a comment counts as wide.
        wide = sum(1 for line in lines if len(line) > WIDTH and split(path, line)[1])
        row = {"lines": len(lines), "comm": len(comments),
               "semi": sum(1 for c in comments if ";" in c), "wide": wide}
        for k, pat in NARRATION.items():
            row[k] = sum(1 for c in comments if re.search(pat, c, re.I))
        for k in tot:
            tot[k] += row[k]
        pct = 100 * row["comm"] // max(row["lines"], 1)
        print(f"{pathlib.Path(path).name:40s} {row['lines']:6d} {row['comm']:5d} {pct:3d} "
              + " ".join(f"{row[k]:5d}" for k in keys) + f" {row['semi']:5d} {row['wide']:5d}")
    pct = 100 * tot["comm"] // max(tot["lines"], 1)
    print(f"{'total':40s} {tot['lines']:6d} {tot['comm']:5d} {pct:3d} "
          + " ".join(f"{tot[k]:5d}" for k in keys) + f" {tot['semi']:5d} {tot['wide']:5d}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--base", help="Revision whose code the files must still have")
    ap.add_argument("--counts", action="store_true", help="Count the comments instead")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    if a.counts:
        counts(a.files)
        return 0
    if not a.base:
        ap.error("give --base <rev> or --counts")
    return 1 if compare(a.base, a.files) else 0


if __name__ == "__main__":
    sys.exit(main())
