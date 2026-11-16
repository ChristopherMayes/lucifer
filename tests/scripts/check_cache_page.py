#!/usr/bin/env python3
"""
Check that the user guide's table of the cache directory agrees with the code.

  check_cache_page.py [--root <repository root>]

The program and the validation suite keep a few files outside any run, in the directory
lucifer of the user's cache directory, and docs/user-guide.md lists them under the label
ug-cache: each file, what writes and reads it, and the key or environment variable that
places it elsewhere. A file added to that directory in the code and not on the page is a
file a user finds and cannot account for. This check reads both.

From the code, src/, programs/, tests/, examples/, util/ and python/, it takes every
name written after the directory: `.cache}/lucifer/<name>` and `.cache/lucifer/<name>`
in a script, and '/lucifer/<name>' in a Fortran string, where the name may be the
opening of a longer one that the program completes, as fft-wisdom- is. From the page it
takes the first column of the table and every LUCIFER_ variable the table names.

Three failures, each printed:

  a name the code places in the directory that opens no file of the table,
  a file of the table that no name in the code opens,
  a variable of the table that the code does not read.

It needs no build, so it runs in the docs workflow and in the examples job.
"""
import argparse
import re
import sys
from pathlib import Path

CODE_DIRS = ("src", "programs", "tests", "examples", "util", "python")
SUFFIXES = {".f90", ".sh", ".py", ".mm", ".c", ".h"}
IN_SCRIPT = re.compile(r"\.cache\}?/lucifer/([A-Za-z0-9._-]+)")
IN_FORTRAN = re.compile(r"'/lucifer/([A-Za-z0-9._-]+)'")
VARIABLE = re.compile(r"\bLUCIFER_[A-Z0-9_]+\b")


def code_names(root):
    """The names the code writes after the cache directory, with where each was found."""
    found = {}
    me = Path(__file__).resolve()
    for d in CODE_DIRS:
        for path in sorted((root / d).rglob("*")):
            if not path.is_file() or path.resolve() == me:
                continue
            if path.suffix not in SUFFIXES and path.name not in ("conda_compile", "bmad_release_check"):
                continue
            try:
                text = path.read_text()
            except UnicodeDecodeError:
                continue
            for n, line in enumerate(text.splitlines(), 1):
                for pattern in (IN_SCRIPT, IN_FORTRAN):
                    for m in pattern.finditer(line):
                        found.setdefault(m.group(1), f"{path.relative_to(root)}:{n}")
    return found


def code_text(root):
    """Every line of the code as one string, for the variables' lookup."""
    parts = []
    for d in CODE_DIRS:
        for path in sorted((root / d).rglob("*")):
            if path.is_file() and (path.suffix in SUFFIXES or path.name in ("conda_compile",)):
                try:
                    parts.append(path.read_text())
                except UnicodeDecodeError:
                    pass
    return "\n".join(parts)


def page_table(root):
    """The files of the page's table, and the LUCIFER_ variables its rows name."""
    text = (root / "docs" / "user-guide.md").read_text()
    start = text.index("(ug-cache)=")
    end = text.index("\n## ", text.index("\n## ", start) + 1)
    files, variables = [], set()
    for line in text[start:end].splitlines():
        if not line.startswith("| `"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        files.append(cells[0].strip("`").rstrip("/"))
        variables.update(VARIABLE.findall(line))
    return files, variables


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]))
    args = ap.parse_args()
    root = Path(args.root).resolve()

    names = code_names(root)
    files, variables = page_table(root)
    text = code_text(root)
    bad = 0
    for name, where in sorted(names.items()):
        if not any(f.startswith(name.rstrip("/")) for f in files):
            print(f"{where}: the code places '{name}' in the cache directory, and the table of "
                  "docs/user-guide.md (ug-cache) has no such file")
            bad += 1
    for f in files:
        if not any(f.startswith(name.rstrip("/")) for name in names):
            print(f"docs/user-guide.md: the table lists '{f}', and no code places it in the cache directory")
            bad += 1
    for v in sorted(variables):
        if v not in text:
            print(f"docs/user-guide.md: the table names {v}, and the code does not read it")
            bad += 1
    print(f"check_cache_page: {len(names)} names in the code, {len(files)} files and "
          f"{len(variables)} variables in the table, {bad} failures")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
