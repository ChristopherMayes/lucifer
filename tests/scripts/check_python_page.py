#!/usr/bin/env python3
"""
Check that every Python example on docs/python.md is code a test runs.

  check_python_page.py [--root <repository root>]

Each ```python block on the page follows a MyST comment line naming it,
"% example: <name>", and the code between "# example: <name>" and
"# end example: <name>" in one of SOURCES is the same text once the region's common
indentation is removed. pytest runs python/tests/test_page.py's regions, the keystone
runs check_program.py's, which call the program, and the keystone's examples job runs
python_page_figures.py's on the sase example's run, saving the figure each draws. So an
example that stops running fails a test, and one edited on the page alone fails here.

Refused: a Python block with no name, a name no file holds or two files hold, a block
that differs from its region, a region no block on the page shows, an image the page
shows that does not exist, and a figure python_page_figures.py saves that the page does
not show.

Standard library only, so the docs job runs it without an environment.
"""
import argparse
import re
import sys
import textwrap
from pathlib import Path

PAGE = "docs/python.md"
SOURCES = ("python/tests/test_page.py", "tests/scripts/check_program.py",
           "tests/scripts/python_page_figures.py")
FIGURES = "tests/scripts/python_page_figures.py"
IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+\.svg)\)")
SAVED = re.compile(r'\bsave\("([\w-]+)"\)')
NAME = re.compile(r"^% example: ([\w-]+)\s*$")
START = re.compile(r"^\s*# example: ([\w-]+)\s*$")
END = re.compile(r"^\s*# end example: ([\w-]+)\s*$")


def page_blocks(text):
    """[(name or None, line number, code)] for every ```python block of the page."""
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        if lines[i].strip() == "```python":
            j = i + 1
            while j < len(lines) and lines[j].strip() != "```":
                j += 1
            k = i - 1
            while k >= 0 and not lines[k].strip():
                k -= 1
            m = NAME.match(lines[k]) if k >= 0 else None
            out.append((m.group(1) if m else None, i + 1, "\n".join(lines[i + 1:j])))
            i = j + 1
        else:
            i += 1
    return out


def regions(path, text):
    """{name: (file, code)} for every marked region of one source file."""
    out, lines, i = {}, text.split("\n"), 0
    while i < len(lines):
        m = START.match(lines[i])
        if not m:
            i += 1
            continue
        name, j = m.group(1), i + 1
        while j < len(lines) and not (END.match(lines[j]) and END.match(lines[j]).group(1) == name):
            j += 1
        if j == len(lines):
            raise SystemExit(f"{path}:{i + 1}: example {name} has no end marker")
        out[name] = (path, textwrap.dedent("\n".join(lines[i + 1:j])).strip("\n"))
        i = j + 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    root = Path(ap.parse_args().root).resolve()

    held = {}
    failed = []
    for src in SOURCES:
        for name, where in regions(src, (root / src).read_text()).items():
            if name in held:
                failed.append(f"example {name} is in both {held[name][0]} and {src}")
            held[name] = where
    shown = set()
    for name, line, code in page_blocks((root / PAGE).read_text()):
        if name is None:
            failed.append(f"{PAGE}:{line}: a Python block with no '% example: <name>' line above it")
            continue
        shown.add(name)
        if name not in held:
            failed.append(f"{PAGE}:{line}: example {name} is in none of {', '.join(SOURCES)}")
        elif code.strip("\n") != held[name][1]:
            failed.append(f"{PAGE}:{line}: example {name} differs from its region in {held[name][0]}")
    for name in sorted(set(held) - shown):
        failed.append(f"{held[name][0]}: example {name} is shown on no block of {PAGE}")
    # The figures: every image the page shows exists, and every figure the script saves
    # is shown.
    page_text = (root / PAGE).read_text()
    images = IMAGE.findall(page_text)
    for img in images:
        if not (root / PAGE).parent.joinpath(img).is_file():
            failed.append(f"{PAGE} shows {img}, which does not exist")
    for name in SAVED.findall((root / FIGURES).read_text()):
        if f"generated/python/{name}.svg" not in images:
            failed.append(f"{FIGURES} saves {name}.svg, which {PAGE} does not show")
    if failed:
        print("\n".join(failed))
        print(f"check_python_page: FAIL, {len(failed)} problems")
        return 1
    print(f"check_python_page: {len(shown)} examples and {len(images)} figures on {PAGE}, "
          f"each the code a test runs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
