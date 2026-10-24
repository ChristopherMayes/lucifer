#!/usr/bin/env python3
"""
Check that every relative link in the committed Markdown files resolves to a file or a
directory in the tree.

  check_links.py [--root <repository root>]

The files are the ones git tracks with a .md suffix, less the generated pages under
docs/generated, whose links report_examples.py retargets when it writes them. Each
inline link or image, [text](target) or ![text](target), is read outside fenced code
blocks and inline code spans. A target with a scheme, http:, https: or mailto:, an
absolute path, or an anchor alone is not a relative link and is passed over. The file
part of a target is what is checked: the anchor after # and any query after ? are
removed, percent-encoding is decoded, and the rest is resolved against the file's own
directory. A link whose target does not exist is reported as file:line: target, and the
check fails if any is.

The check exists because the site can be whole while the tree is not: the example pages
are generated from the examples' READMEs with their links retargeted, so a README whose
link points at a directory that no longer exists reads correctly on the site and is
dead on GitHub. On 2026-10-24 forty-one links pointed at doc/, the pages' directory
before it became docs/, and nothing had said so.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

# A target is bare, or in angle brackets when it holds spaces or parentheses, as a DOI does.
LINK = re.compile(r"!?\[[^\]]*\]\((?:<([^>]*)>|([^)\s]+))(?:\s+\"[^\"]*\")?\)")
FENCE = re.compile(r"^\s*(```|~~~)")
CODE_SPAN = re.compile(r"`[^`]*`")
SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")


def tracked_markdown(root):
    out = subprocess.run(["git", "ls-files", "--", "*.md"], cwd=root, check=True,
                         capture_output=True, text=True).stdout.split("\n")
    return [root / f for f in out if f and not f.startswith("docs/generated/")]


def links_of(path):
    """(line number, target) for each inline link outside code."""
    found = []
    fenced = False
    for n, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        if FENCE.match(line):
            fenced = not fenced
            continue
        if fenced:
            continue
        for bracketed, bare in LINK.findall(CODE_SPAN.sub("", line)):
            found.append((n, bracketed or bare))
    return found


def is_relative(target):
    if SCHEME.match(target) or target.startswith("/") or target.startswith("#"):
        return False
    return True


def resolve(path, target):
    """The file part of the target, against the file's own directory."""
    file_part = target.split("#", 1)[0].split("?", 1)[0]
    if not file_part:
        return path
    return (path.parent / unquote(file_part)).resolve()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    args = ap.parse_args()
    root = Path(args.root).resolve()

    n_files = n_links = 0
    dead = []
    for md in tracked_markdown(root):
        n_files += 1
        for n, target in links_of(md):
            if not is_relative(target):
                continue
            n_links += 1
            if not resolve(md, target).exists():
                dead.append((md.relative_to(root), n, target))

    for rel, n, target in dead:
        print(f"{rel}:{n}: {target}")
    print(f"check_links: {n_links} relative links in {n_files} Markdown files, {len(dead)} dead")
    return 1 if dead else 0


if __name__ == "__main__":
    sys.exit(main())
