#!/usr/bin/env python3
"""
Check that every page of the site's contents declares its kind.

  check_page_kinds.py [--root <repository root>] [--list]

A page of the manual is one of five kinds, named by what its reader is doing:

  tutorial      a first run of the program, followed step by step
  how-to        the steps of a task the reader already has
  reference     one fact to look up: a key, a routine, a file layout, a term
  explanation   the physics or the numerical method, which changes when the method does
  record        what was measured, on which example, computer and day

The kind is declared in the front matter of the page, under the key tags, which holds
kinds and nothing else:

  ---
  tags: [explanation, record]
  ---

A page is one kind, or it is an explanation or a how-to guide followed by its record, and
its tags are then those two in that order. The key is tags since MyST accepts it. A key
of another name is reported as an extra key at every build of the site.

The pages are the files that docs/myst.yml lists under its toc. A Markdown page declares
its kind in the front matter that opens the file. A notebook declares it in a front
matter block at the head of its first Markdown cell, which MyST reads in the same way. A
generated page is written with its tags by its generator. Each failure is printed as
file: the reason.

--list prints every page with its tags.

Standard library only, so the docs job runs it without an environment.
"""
import argparse
import json
import re
import sys
from pathlib import Path

KINDS = ("tutorial", "how-to", "reference", "explanation", "record")
PAIRS = (("explanation", "record"), ("how-to", "record"))
TOC_FILE = re.compile(r"^\s*-\s*file:\s*(\S+)\s*$")
TAGS_FLOW = re.compile(r"^tags:\s*\[(.*)\]\s*$")
TAGS_BLOCK = re.compile(r"^tags:\s*$")
BLOCK_ITEM = re.compile(r"^\s+-\s*(\S.*?)\s*$")


def contents(root):
    """The files of the toc of docs/myst.yml, in its order."""
    files = []
    in_toc = False
    for line in (root / "docs" / "myst.yml").read_text(encoding="utf-8").split("\n"):
        if re.match(r"^\s*toc:\s*$", line):
            in_toc = True
            continue
        if in_toc and re.match(r"^\S", line):
            break
        m = TOC_FILE.match(line)
        if in_toc and m:
            files.append(m.group(1))
    return files


def head_text(path):
    """The text that may open with front matter: the file, or a notebook's first Markdown cell."""
    if path.suffix != ".ipynb":
        return path.read_text(encoding="utf-8")
    for cell in json.loads(path.read_text(encoding="utf-8")).get("cells", []):
        if cell.get("cell_type") == "markdown":
            src = cell.get("source", "")
            return "".join(src) if isinstance(src, list) else src
    return ""


def tags_of(text):
    """The values of tags in the front matter, or None where the page has no tags."""
    lines = text.split("\n")
    if not lines or lines[0].strip() != "---":
        return None
    try:
        end = next(i for i in range(1, len(lines)) if lines[i].strip() == "---")
    except StopIteration:
        return None
    block = lines[1:end]
    for i, line in enumerate(block):
        m = TAGS_FLOW.match(line)
        if m:
            return [t.strip().strip("\"'") for t in m.group(1).split(",") if t.strip()]
        if TAGS_BLOCK.match(line):
            found = []
            for item in block[i + 1:]:
                mi = BLOCK_ITEM.match(item)
                if not mi:
                    break
                found.append(mi.group(1).strip("\"'"))
            return found
    return None


def problem(tags):
    """Why the tags do not declare a kind, or None where they do."""
    if tags is None:
        return "no tags in its front matter"
    if not tags:
        return "tags is empty"
    other = [t for t in tags if t not in KINDS]
    if other:
        return f"tags holds {', '.join(other)}, which is not one of {', '.join(KINDS)}"
    if len(tags) == 1:
        return None
    if tuple(tags) in PAIRS:
        return None
    return (f"tags is {', '.join(tags)}: a page is one kind, or explanation or how-to "
            "followed by record")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2], type=Path)
    ap.add_argument("--list", action="store_true", help="print every page with its tags and stop")
    args = ap.parse_args()
    root = args.root.resolve()
    files = contents(root)
    if not files:
        print("check_page_kinds: no file found under the toc of docs/myst.yml")
        return 1
    bad = 0
    count = {k: 0 for k in KINDS}
    for f in files:
        path = root / "docs" / f
        if not path.exists():
            print(f"docs/{f}: listed in docs/myst.yml and not in the tree")
            bad += 1
            continue
        tags = tags_of(head_text(path))
        if args.list:
            print(f"{', '.join(tags) if tags else '-':22s} docs/{f}")
            continue
        why = problem(tags)
        if why:
            print(f"docs/{f}: {why}")
            bad += 1
        else:
            for t in tags:
                count[t] += 1
    if args.list:
        return 0
    if bad:
        print(f"check_page_kinds: {bad} of {len(files)} pages do not declare a kind. Give the page "
              "tags in its front matter, or have its generator write them.")
        return 1
    print(f"check_page_kinds: {len(files)} pages declare a kind (" +
          ", ".join(f"{k} {count[k]}" for k in KINDS) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())
