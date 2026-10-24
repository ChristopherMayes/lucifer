#!/usr/bin/env python3
"""
Check that no committed text carries a word this repository has retired.

  check_retired_words.py [--root <repository root>]

The files are the ones git tracks with a Markdown, Python, shell, Fortran, YAML or
CMake suffix, less docs/changelog.md, whose older entries are history and keep their
words, and docs/generated/, which is written from the Fortran headers the check reads
at their source. Each line carrying a retired word is printed as file:line: the line,
and the check fails if any is.

The words, and why each is retired. They are this port's own coinages or testing
jargon, absent from the Bmad corpus the prose is measured against, and each had spread
through the pages and comments before it was noticed:

  harness          testing jargon for tests/run_fel_benchmark.sh and its check scripts,
                   which are the benchmark, and for tests/run_perf_benchmark.sh, which
                   is the performance benchmark.
  gate, gates      process jargon for a check, a test, a tolerance or a refusal.
  refused by name  this port's coinage; the run is refused, and the message says why.

A retired word in a quoted title or a URL is still a hit, so the list stays short and
the words plain.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

RETIRED = [
    ("harness", re.compile(r"\bharness(es|'s)?\b", re.I)),
    ("gate", re.compile(r"\bgates?\b", re.I)),
    ("refused by name", re.compile(r"\brefused by name\b", re.I)),
]
SUFFIXES = (".md", ".py", ".sh", ".f90", ".yml", ".yaml", ".cmake", ".txt")
SKIP = ("docs/changelog.md", "docs/generated/")


def tracked(root):
    out = subprocess.run(["git", "ls-files"], cwd=root, check=True,
                         capture_output=True, text=True).stdout.split("\n")
    return [root / f for f in out if f.endswith(SUFFIXES) and not f.startswith(SKIP)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    args = ap.parse_args()
    root = Path(args.root).resolve()

    hits = {}
    n_files = 0
    for path in tracked(root):
        if path.resolve() == Path(__file__).resolve():
            continue
        n_files += 1
        for n, line in enumerate(path.read_text(encoding="utf-8", errors="replace").split("\n"), 1):
            for word, pat in RETIRED:
                if pat.search(line):
                    hits.setdefault(word, []).append((path.relative_to(root), n, line.strip()))
    total = 0
    for word, _ in RETIRED:
        for rel, n, text in hits.get(word, []):
            print(f"{rel}:{n}: {text[:120]}")
        total += len(hits.get(word, []))
    counts = ", ".join(f"{word} {len(hits.get(word, []))}" for word, _ in RETIRED)
    print(f"check_retired_words: {n_files} files, {total} lines with a retired word ({counts})")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
