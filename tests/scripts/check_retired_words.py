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

  harness          testing jargon for tests/run_fel_validation.sh and its check scripts,
                   which are the validation suite, and for tests/run_perf_benchmark.sh, which
                   is the performance benchmark.
  gate, gates      process jargon for a check, a test, a tolerance or a refusal.
  refused by name  this port's coinage; the run is refused, and the message says why.
  global% and the other namelist prefixes, the &fel_ group names and lucifer.in
                   the namelist input, replaced by the TOML input file on 2026-10-24,
                   whose keys are table.name and whose examples' files are lucifer.toml.
  deck, decks      the punched-card word for an input file, which the Bmad corpus never
                   uses: Bmad says input file, and master input file for a program's own,
                   and Genesis4 says main input file. The file a run reads is the input
                   file, and its keys are input keys. An identifier carrying the word is
                   a hit too, make_deck or Deck, since a name is read as a word.

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
    # The namelist deck's names, gone since 2026-10-24: the deck is TOML, and its keys are
    # table.name, run.out_root for what was global%out_root.
    # Not inside a Fortran component chain, run%global%out_root, which is code.
    ("namelist prefix", re.compile(r"(?<![%\w])(global|beam_init|wavefront_init|slicing|chamber_wake|resample)%")),
    ("namelist group", re.compile(r"&fel_(params|beam_init|wavefront_init|track_params)\b")),
    ("lucifer.in", re.compile(r"\blucifer\.in\b")),
    # The word alone or inside an identifier: deck, Deck, DECK, make_deck, deck_path,
    # GENESIS_DECK, DeckOps. Not inside a longer word, so bedeck and decker would pass.
    ("deck", re.compile(r"(?:\b|_|(?<=[a-z]))[Dd][Ee][Cc][Kk]s?(?=\b|_|[A-Z])|\bDECKS?(?=\b|_)")),
]
SUFFIXES = (".md", ".py", ".sh", ".f90", ".yml", ".yaml", ".cmake", ".txt", ".c", ".h", ".mm", ".in",
            ".toml", ".tsv", ".json", ".bmad", ".lat")
SKIP = ("docs/changelog.md", "docs/generated/")


def fortran_prose(line):
    """The comment and the string literals of a Fortran line, the code left out."""
    if re.match(r"\s*!", line):
        return line
    out, i, n = [], 0, len(line)
    while i < n:
        ch = line[i]
        if ch == "'":
            j = i + 1
            while j < n:
                if line[j] == "'":
                    if j + 1 < n and line[j + 1] == "'":
                        j += 2
                        continue
                    break
                j += 1
            out.append(line[i + 1:j])
            i = j + 1
        elif ch == "!":
            out.append(line[i:])
            break
        else:
            i += 1
    return " ".join(out)


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
        # The converter of flat text knows the old namelist names on purpose: it routes them
        # to the new. Every other word is checked there as anywhere.
        converter = path.relative_to(root).as_posix() in ("python/lucifer/deck.py", "python/lucifer/flat.py")
        n_files += 1
        fortran = path.suffix == ".f90"
        for n, line in enumerate(path.read_text(encoding="utf-8", errors="replace").split("\n"), 1):
            for word, pat in RETIRED:
                if converter and word.startswith("namelist"):
                    continue
                # In Fortran, a retired deck name is prose: a comment or a message. A local
                # variable named like the old structure, beam_init%n_particle in a routine
                # that takes beam_init, is code and stays.
                text = fortran_prose(line) if fortran and word.startswith("namelist") else line
                if pat.search(text):
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
