#!/usr/bin/env python3
"""
Examples conformance: the feature list, the matrix and the directories must agree.

doc/index.md declares what the program does and examples/README.md says which example
shows each declared feature. Neither statement is checkable on its own, and a feature
that quietly loses its example reads exactly like one that never had it. Six assertions:

  1. Every Features bullet is represented by a matrix row. The row's name carries the
     bullet's own bold lead phrase, so the two files use one vocabulary.
  2. No matrix row invents a feature: a row either carries a bullet's lead phrase or
     names something the Features prose says.
  3. Every directory a matrix row points at exists. A row with no link at all is
     allowed only in the two documented prose forms below.
  4. Every example directory carries a README and a deck, with lucifer.in the primary
     name wherever a reader rather than a runner picks the deck.
  5. Every deck fenced into a generated page is byte-identical to the file on disk,
     which is what makes the pages copy-paste safe.
  6. Every generated page's figure exists.

Nothing here runs the program, so it is also the part of the harness that a machine with
no Fortran build can run.

Usage:

  check_examples.py [--lucifer-dir <dir>]

Exits nonzero on any problem.
"""

import argparse
import re
import sys
from pathlib import Path

# A Shown-by cell with no directory link is allowed only where the feature genuinely has
# no directory. Both forms are spelled out here so that a third one cannot appear by
# accident: adding to this list is a decision, not an oversight.
PROSE_CELLS = ("Demonstrated inline, below", "Every example writes one")

FENCE = re.compile(r"^### `([^`]+)`\n\n```[a-z]*\n(.*?)\n```", re.S | re.M)
FIGURE = re.compile(r"^!\[[^\]]*\]\(([^)]+\.png)\)", re.M)

# The generated index page, whose source is examples/README.md rather than one
# example directory. report_examples.py names it, and the two must agree.
INDEX_NAME = "examples"
LINKED_DIR = re.compile(r"\]\(([A-Za-z0-9_]+)/\)")


def section(text, start, end, where):
    if start not in text:
        raise SystemExit(f"{where}: no {start!r} heading")
    body = text.split(start, 1)[1]
    return body.split(end, 1)[0] if end and end in body else body


def lead_phrases(features, where):
    """The bold lead of every Features bullet, which is the name the matrix must use."""
    out = []
    for line in features.splitlines():
        if not line.startswith("- "):
            continue
        m = re.search(r"\*\*(.+?)\*\*", line)
        if not m:
            raise SystemExit(f"{where}: Features bullet without a bold lead phrase: "
                             f"{line[:60]!r}")
        out.append(m.group(1).strip().rstrip(".").rstrip(","))
    if not out:
        raise SystemExit(f"{where}: no Features bullets")
    return out


def rows(table, where, want_cols=2):
    """(name, cell) for each body row of a pipe table."""
    out = []
    for line in table.splitlines():
        s = line.strip()
        if not s.startswith("|") or s.startswith("|---"):
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if len(cells) != want_cols:
            raise SystemExit(f"{where}: row has {len(cells)} cells, expected {want_cols}: "
                             f"{s[:70]!r}")
        if cells[0] in ("Feature", "Planned", "Example"):
            continue
        out.append((cells[0], cells[1]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lucifer-dir", default=str(Path(__file__).resolve().parents[2]),
                   help="The lucifer directory (default: two levels above this script)")
    args = p.parse_args()

    luc = Path(args.lucifer_dir)
    ex = luc / "examples"
    index_md = luc / "doc" / "index.md"
    matrix_md = ex / "README.md"
    pages = luc / "doc" / "generated" / "examples"

    problems = []
    idx_text = index_md.read_text()
    mtx_text = matrix_md.read_text()

    features = section(idx_text, "## Features", "## Known missing features", index_md)
    leads = lead_phrases(features, index_md)
    matrix = section(mtx_text, "## What each feature is exemplified by", "## Planned examples",
                     matrix_md)
    planned = section(mtx_text, "## Planned examples", "## What a run writes", matrix_md)

    matrix_rows = rows(matrix, f"{matrix_md} matrix")
    planned_rows = rows(planned, f"{matrix_md} planned")

    dirs = sorted(d.name for d in ex.iterdir() if d.is_dir())

    # 1. Every bullet is represented.
    def norm(s):
        return s.lower()
    for lead in leads:
        if not any(norm(lead) in norm(name) for name, _ in matrix_rows):
            problems.append(f"Features bullet {lead!r} has no matrix row naming it")

    # 2. No row invents a feature.
    feat_prose = norm(re.sub(r"\*\*", "", features))
    for name, _ in matrix_rows:
        if any(norm(lead) in norm(name) for lead in leads):
            continue
        if norm(name) in feat_prose:
            continue
        problems.append(f"matrix row {name!r} names nothing the Features list declares")

    # 3. Every row points somewhere real.
    for name, cell in matrix_rows:
        linked = LINKED_DIR.findall(cell)
        if not linked:
            if cell not in PROSE_CELLS:
                problems.append(f"matrix row {name!r} links no example directory and is "
                                f"not one of the documented prose forms: {cell!r}")
            continue
        for d in linked:
            if d not in dirs:
                problems.append(f"matrix row {name!r} points at examples/{d}/, which does "
                                "not exist")

    # 3b. A planned feature must not already have a directory, and must say what it waits on.
    for name, cell in planned_rows:
        if not cell:
            problems.append(f"planned row {name!r} does not say what it waits on")
        for d in LINKED_DIR.findall(name) + LINKED_DIR.findall(cell):
            if d in dirs:
                problems.append(f"planned row {name!r} names examples/{d}/, which exists: "
                                "a planned example with a directory is a contradiction")

    # 4. Every directory is a complete example. A directory run by hand carries
    # lucifer.in as its primary deck. One with its own runner may name its decks by
    # role instead, since the runner rather than the reader chooses them.
    for d in dirs:
        if not (ex / d / "README.md").is_file():
            problems.append(f"examples/{d}/ has no README.md")
        decks = sorted(f.name for f in (ex / d).glob("lucifer*.in"))
        if not decks:
            problems.append(f"examples/{d}/ has no lucifer*.in deck")
        elif "lucifer.in" not in decks and not (ex / d / "run.sh").is_file():
            problems.append(f"examples/{d}/ has no lucifer.in and no run.sh to name its "
                            f"decks for it: {decks}")

    # 5 and 6. The generated pages quote the real files, and their figures exist.
    n_fenced = 0
    if not pages.is_dir():
        problems.append(f"{pages} does not exist: run report_examples.py")
    else:
        for page in sorted(pages.glob("*.md")):
            name = page.stem
            text = page.read_text()
            src_dir = ex if name == INDEX_NAME else ex / name
            if name != INDEX_NAME:
                if not (src_dir).is_dir():
                    problems.append(f"generated page {page.name} has no example directory")
                    continue
                # A page requires a figure exactly when it references one, so the two
                # statements cannot disagree. A page with no figure is a deliberate
                # choice report_examples.py records, not an omission to detect here.
                for fig_name in FIGURE.findall(text):
                    if not (pages / fig_name).is_file():
                        problems.append(f"generated page {page.name} references "
                                        f"{fig_name}, which does not exist: run "
                                        "run_examples.sh")
            for fname, body in FENCE.findall(text):
                src = src_dir / fname
                if not src.is_file():
                    problems.append(f"{page.name} fences {fname}, which is not a file in "
                                    f"{src_dir}")
                    continue
                if body.rstrip("\n") != src.read_text().rstrip("\n"):
                    problems.append(f"{page.name} fences {fname} but the text differs from "
                                    "the file: regenerate with report_examples.py")
                n_fenced += 1

    print(f"  examples: {len(dirs)} directories, {len(leads)} declared features, "
          f"{len(matrix_rows)} matrix rows, {len(planned_rows)} planned, "
          f"{n_fenced} fenced files verified")
    for msg in problems:
        print(f"  PROBLEM: {msg}")
    if problems:
        print(f"examples conformance: FAIL ({len(problems)} problem(s))")
        return 1
    print("checks: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
