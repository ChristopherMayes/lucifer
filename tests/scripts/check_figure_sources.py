#!/usr/bin/env python3
"""
Check that every image a page shows has a script that draws it, and that an image drawn
by hand-run script is not older than what it was drawn from.

  check_figure_sources.py [--root <repository root>] [--list]
  check_figure_sources.py --update [<image> ...]

Every figure of the manual is drawn by a committed script, and none is drawn by hand. The
images are committed, and no test compares them, since the fonts and the plotting library
differ from one computer to the next and rewrite an image with nothing in it changed.
This check compares what an image was drawn from.

docs/figures/sources.tsv lists the images in three classes. Its columns are separated by
tabs: the image, the class, the script, the files the script reads, and a hash.

  keystone   The keystone's regeneration redraws the image at every run. The image may be
             a pattern. No hash is kept.
  examples   A full run of examples/run_examples.sh redraws the image, and that run is
             made whenever a recorded digit or a plot routine changes. The image may be a
             pattern. No hash is kept.
  hand       A script that is run by hand draws the image. The line names the image, the
             script and every file the script reads to draw it, separated by commas, or a
             hyphen where it reads none. The hash is the SHA-256 of the script and of
             those files.

Four things fail the check, each printed as file: the reason.

  A page shows an image that no line names. The pages are the Markdown files under docs/
  that git tracks, and an image is the target of a figure directive or of ![text](target).
  An image inside a notebook's stored output is not read.

  A line of the class hand whose hash differs from the hash of its script and its files as
  they stand. The script or a file it reads was changed after the image was drawn. The
  image is redrawn, and --update then rewrites the hash of its line.

  A line that names a script, a file or, in the class hand, an image that is not in the
  tree.

  A line that is not in the form above.

--update rewrites the hash of each line of the class hand, or of the images named, and
changes nothing else. It is run after the image is redrawn and looked at. Running it
without redrawing records a hash for an image that was not drawn from those files.

--list prints each image a page shows with its class and its script.

Standard library only, so the docs job runs it without an environment.
"""
import argparse
import fnmatch
import hashlib
import re
import subprocess
import sys
from pathlib import Path

SOURCES = "docs/figures/sources.tsv"
CLASSES = ("keystone", "examples", "hand")
FIGURE = re.compile(r"\{figure\}\s+(\S+)")
IMAGE = re.compile(r"!\[[^\]]*\]\((?:<([^>]*)>|([^)\s]+))(?:\s+\"[^\"]*\")?\)")
SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")


def read_sources(root):
    """(lines of the file, rows). A row is [line index, image, class, script, reads, hash]."""
    lines = (root / SOURCES).read_text(encoding="utf-8").split("\n")
    rows, bad = [], []
    for i, line in enumerate(lines):
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 5 or cols[1] not in CLASSES:
            bad.append(f"{SOURCES}:{i + 1}: not five columns with a class among {', '.join(CLASSES)}")
            continue
        reads = [] if cols[3] == "-" else [r for r in cols[3].split(",") if r]
        rows.append([i, cols[0], cols[1], cols[2], reads, cols[4]])
    return lines, rows, bad


def digest(root, script, reads):
    """SHA-256 over the script and the files it reads, each under its own path."""
    h = hashlib.sha256()
    for rel in [script, *reads]:
        h.update(rel.encode("utf-8") + b"\0")
        h.update((root / rel).read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def shown_images(root):
    """{image path relative to the root: the first page that shows it}."""
    out = subprocess.run(["git", "ls-files", "--", "docs/*.md"], cwd=root, check=True,
                         capture_output=True, text=True).stdout.split("\n")
    shown = {}
    for f in sorted(x for x in out if x):
        page = root / f
        text = page.read_text(encoding="utf-8")
        targets = FIGURE.findall(text) + [a or b for a, b in IMAGE.findall(text)]
        for target in targets:
            if SCHEME.match(target) or target.startswith("#"):
                continue
            target = target.split("#", 1)[0].split("?", 1)[0]
            # A generated page that another page includes names its images from docs/, and
            # every other page names them from its own directory.
            for base in (page.parent, root / "docs"):
                path = (base / target).resolve()
                if path.exists():
                    break
            else:
                path = (page.parent / target).resolve()
            try:
                rel = path.relative_to(root).as_posix()
            except ValueError:
                rel = target
            shown.setdefault(rel, f)
    return shown


def row_of(rows, image):
    for row in rows:
        if row[1] == image or (row[2] != "hand" and fnmatch.fnmatchcase(image, row[1])):
            return row
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2], type=Path)
    ap.add_argument("--list", action="store_true", help="print each image a page shows and stop")
    ap.add_argument("--update", nargs="*", metavar="image", default=None,
                    help="rewrite the hash of every line of the class hand, or of these images")
    args = ap.parse_args()
    root = args.root.resolve()
    lines, rows, problems = read_sources(root)

    if args.update is not None:
        if problems:
            print("\n".join(problems))
            return 1
        names = set(args.update)
        unknown = names - {row[1] for row in rows if row[2] == "hand"}
        if unknown:
            print(f"check_figure_sources: no line of the class hand names {', '.join(sorted(unknown))}")
            return 1
        n = 0
        for i, image, cls, script, reads, old in rows:
            if cls != "hand" or (names and image not in names):
                continue
            new = digest(root, script, reads)
            if new != old:
                cols = lines[i].split("\t")
                cols[4] = new
                lines[i] = "\t".join(cols)
                n += 1
        (root / SOURCES).write_text("\n".join(lines), encoding="utf-8")
        print(f"check_figure_sources: {n} hashes rewritten in {SOURCES}")
        return 0

    shown = shown_images(root)
    if args.list:
        for image, page in sorted(shown.items()):
            row = row_of(rows, image)
            print(f"{row[2] if row else '-':9s} {image}  {row[3] if row else ''}")
        return 0

    for i, image, cls, script, reads, old in rows:
        where = f"{SOURCES}:{i + 1}"
        missing = [f for f in [script, *reads] if not (root / f).is_file()]
        if missing:
            problems.append(f"{where}: not in the tree: {', '.join(missing)}")
            continue
        if cls != "hand":
            continue
        if not (root / image).is_file():
            problems.append(f"{where}: the image {image} is not in the tree")
        elif digest(root, script, reads) != old:
            problems.append(f"{image}: {script} or a file it reads was changed after the image was "
                            "drawn. Redraw the image, then run check_figure_sources.py --update "
                            f"{image}")
    count = {c: 0 for c in CLASSES}
    for image, page in sorted(shown.items()):
        row = row_of(rows, image)
        if row is None:
            problems.append(f"{page}: shows {image}, which no line of {SOURCES} names. Every "
                            "figure is drawn by a committed script")
        else:
            count[row[2]] += 1
    if problems:
        print("\n".join(problems))
        print(f"check_figure_sources: {len(problems)} failures")
        return 1
    print(f"check_figure_sources: {len(shown)} images shown by the pages, each with its script (" +
          ", ".join(f"{c} {count[c]}" for c in CLASSES) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())
