#!/usr/bin/env python3
"""
Check that the manual names no Genesis4 routine or source file without a link.

  check_genesis_names.py [--root <repository root>]

The files are the Markdown files git tracks, less docs/changelog.md, whose entries are
history, and docs/generated/, which is written from the Fortran headers and keeps their
citations. Two things fail the check, each printed as file:line: the text.

  A C++ name of the form Class::method, or a file name ending in .cpp, that is not the
  text of a link and is not inside a block opened by ":::{dropdown} Genesis4 source".
  The pages state what Lucifer does. A Genesis4 routine is named where a reader can
  open it.

  A link into the Genesis4 repository that does not name the tag of the reference
  release. The release is the genesis4 version of tests/environment-validate.yml, and
  a line number is only right at one release.

The CMake target Lucifer::lucifer_lib is not a Genesis4 name and passes.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

SKIP = ("docs/changelog.md", "docs/generated/")
REPOSITORY = "github.com/svenreiche/Genesis-1.3-Version4/"
NAME = re.compile(r"\b[A-Z][A-Za-z0-9]*::[A-Za-z_][A-Za-z0-9_]*|\b[A-Za-z0-9_]+\.cpp\b")
LINK = re.compile(r"\[([^\]]*)\]\(<?([^)\s>]*)>?\)")
SOURCE_LINK = re.compile(re.escape(REPOSITORY) + r"(?:blob|tree)/([^/)\s]+)/")
ALLOWED = ("Lucifer::lucifer_lib",)
BLOCK = re.compile(r"^(:{3,})\{dropdown\}\s+Genesis4 source\s*$")


def reference_tag(root):
    """The tag of the reference release, from the validation environment's genesis4 entry."""
    text = (root / "tests/environment-validate.yml").read_text()
    m = re.search(r"^\s*-\s*genesis4\s*=+\s*([0-9][0-9.]*)", text, re.M)
    if not m:
        raise SystemExit("check_genesis_names: no genesis4 version in tests/environment-validate.yml")
    return "v" + m.group(1)


def tracked_markdown(root):
    out = subprocess.run(["git", "ls-files", "--", "*.md"], cwd=root, check=True,
                         capture_output=True, text=True).stdout.split("\n")
    return [root / f for f in out if f and not f.startswith(SKIP)]


def check_file(path, tag):
    """(line number, reason, text) for every failure in one file."""
    found = []
    fence = None
    for n, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        m = BLOCK.match(line)
        if fence is None and m:
            fence = m.group(1)
            continue
        if fence is not None and line.strip() == fence:
            fence = None
            continue
        for ref in SOURCE_LINK.findall(line):
            if ref != tag:
                found.append((n, f"link at {ref}, the reference release is {tag}", line.strip()))
        if fence is not None:
            continue
        # A name that is the text of a link, or part of its target, is linked.
        bare = LINK.sub(" ", line)
        for name in NAME.findall(bare):
            if name not in ALLOWED:
                found.append((n, f"{name} with no link", line.strip()))
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    args = ap.parse_args()
    root = Path(args.root).resolve()
    tag = reference_tag(root)

    total = 0
    files = tracked_markdown(root)
    for path in files:
        for n, reason, text in check_file(path, tag):
            print(f"{path.relative_to(root)}:{n}: {reason}: {text[:110]}")
            total += 1
    print(f"check_genesis_names: {len(files)} files, {total} Genesis4 names without a link "
          f"or links off {tag}")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
