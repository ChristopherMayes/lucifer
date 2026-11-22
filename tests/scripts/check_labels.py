#!/usr/bin/env python3
"""
Check that every reference to a label of the manual names a label that exists.

  check_labels.py [--root <repository root>] [--list]

A section, an equation, a figure or a table of the manual is referred to by its label.
The site resolves a label wherever its section stands, so a section can move to another
page and keep its references. Nothing else holds a reference to its label: the link
check reads the file of a link and not its anchor, and the site reports a reference that
does not resolve as a warning and builds.

The pages are the files in the toc of docs/myst.yml and the generated pages those
include. A notebook is read by its Markdown cells. A label is one of

  (name)=              before a heading or a paragraph
  :label: name         an option of a directive
  :name: name          the same, in its older spelling
  $$ (name)            the label of an equation, at the closing of its block

Five things fail the check, each printed as file:line: the reason.

  A label defined twice in the manual. A reference to it would be ambiguous, and a
  section could not move to a page that holds its namesake.

  A reference by label, [text](#name), that names no label of the manual.

  A link by file and anchor, [text](page.md#name), whose page is a page of the manual
  that does not define the label, in its own text or in a generated page it includes.

  The same link in README.md or in the README of an example, whose links keep the file
  name since the forge does not resolve a label.

  A comment of the source or of a script that cites a section label, sec-<name>, that
  the manual does not define. The files are the Fortran, C and Objective-C of src/ and
  programs/, and the Python, shell, TOML and Bmad files of tests/, python/, examples/ and
  docs/figures/.

References are read outside fenced code blocks and inline code. docs/changelog.md is
history: its labels count, and its references are not checked, since an entry cites the
page as it was on its day.

--list prints every label with the page that defines it.

Standard library only, so the docs job runs it without an environment.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HISTORY = ("docs/changelog.md",)
TOC_FILE = re.compile(r"^\s*-\s*file:\s*(\S+)\s*$")
INCLUDE = re.compile(r"^\s*(?:```|:::)\{include\}\s+(\S+)")
LABEL_LINE = re.compile(r"^\(([A-Za-z0-9_:.-]+)\)=\s*$")
LABEL_OPTION = re.compile(r"^\s*:(?:label|name):\s*([A-Za-z0-9_:.-]+)\s*$")
LABEL_EQUATION = re.compile(r"\$\$\s*\(([A-Za-z0-9_:.-]+)\)\s*$")
FENCE = re.compile(r"^\s*(```|~~~)(.*)$")
CODE_SPAN = re.compile(r"`[^`]*`")
LINK = re.compile(r"\]\((?:<([^>]*)>|([^)\s]+))(?:\s+\"[^\"]*\")?\)")
SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
SECTION = re.compile(r"(?<![A-Za-z0-9_-])sec-[a-z][a-z0-9-]*[a-z0-9]")
SOURCE_GLOBS = ("src/*.f90", "src/*.c", "src/*.in", "src/device/*", "programs/*.f90",
                "tests/*.sh", "tests/*.py", "tests/scripts/*.py", "tests/scripts/*.sh",
                "tests/**/*.toml", "tests/**/*.bmad", "python/**/*.py", "examples/*.sh",
                "examples/**/*.toml", "examples/**/*.bmad", "docs/figures/*.py")


def contents(root):
    """The files of the toc of docs/myst.yml, as paths from the root."""
    files, in_toc = [], False
    for line in (root / "docs" / "myst.yml").read_text(encoding="utf-8").split("\n"):
        if re.match(r"^\s*toc:\s*$", line):
            in_toc = True
            continue
        if in_toc and re.match(r"^\S", line):
            break
        m = TOC_FILE.match(line)
        if in_toc and m:
            files.append("docs/" + m.group(1))
    return files


def lines_of(path):
    """The Markdown of a page as (line number, line). A notebook gives its Markdown cells."""
    if path.suffix != ".ipynb":
        return list(enumerate(path.read_text(encoding="utf-8").split("\n"), 1))
    out = []
    for n, cell in enumerate(json.loads(path.read_text(encoding="utf-8")).get("cells", []), 1):
        if cell.get("cell_type") == "markdown":
            src = cell.get("source", "")
            for line in ("".join(src) if isinstance(src, list) else src).split("\n"):
                out.append((n, line))
    return out


def scan(path):
    """(labels, references, includes) of one file.

    A label is (line, name). A reference is (line, file part, name), the file part
    empty for a reference by label alone. An include is the path as written."""
    labels, refs, includes = [], [], []
    fence = None
    for n, line in lines_of(path):
        m = FENCE.match(line)
        if m:
            if fence is None:
                inc = INCLUDE.match(line)
                if inc:
                    includes.append(inc.group(1))
                # A directive's body is Markdown and is read. A plain code block is not.
                fence = "directive" if m.group(2).lstrip().startswith("{") else "code"
            else:
                fence = None
            continue
        if fence == "code":
            continue
        for rx in (LABEL_LINE, LABEL_OPTION, LABEL_EQUATION):
            lm = rx.search(line)
            if lm:
                labels.append((n, lm.group(1)))
                break
        inc = INCLUDE.match(line)
        if inc:
            includes.append(inc.group(1))
        for bracketed, bare in LINK.findall(CODE_SPAN.sub("", line)):
            target = bracketed or bare
            if SCHEME.match(target) or "#" not in target:
                continue
            file_part, name = target.split("#", 1)
            refs.append((n, file_part.split("?", 1)[0], name))
    return labels, refs, includes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2], type=Path)
    ap.add_argument("--list", action="store_true", help="print every label with its page and stop")
    args = ap.parse_args()
    root = args.root.resolve()

    # Each page with the files it includes, and every file read once.
    pages = contents(root)
    scanned, members = {}, {}
    for page in pages:
        todo, mine = [page], []
        while todo:
            f = todo.pop()
            if f in mine or not (root / f).is_file():
                continue
            mine.append(f)
            if f not in scanned:
                scanned[f] = scan(root / f)
            for inc in scanned[f][2]:
                todo.append((Path(f).parent / inc).as_posix())
        members[page] = mine

    problems = []
    defined = {}                                  # name -> (file, line) of its first definition
    for f, (labels, _, _) in scanned.items():
        for n, name in labels:
            if name in defined and defined[name] != (f, n):
                problems.append(f"{f}:{n}: the label {name} is also defined at "
                                f"{defined[name][0]}:{defined[name][1]}")
            defined.setdefault(name, (f, n))
    page_labels = {page: {name for f in mine for _, name in scanned[f][0]}
                   for page, mine in members.items()}
    if args.list:
        for name, (f, n) in sorted(defined.items()):
            print(f"{name:44s} {f}:{n}")
        return 0

    n_refs = 0
    for f, (_, refs, _) in scanned.items():
        if f in HISTORY:
            continue
        for n, file_part, name in refs:
            n_refs += 1
            if not file_part:
                if name not in defined:
                    problems.append(f"{f}:{n}: #{name} names no label of the manual")
                continue
            target = (Path(f).parent / file_part).as_posix()
            target = Path(*[p for p in Path(target).parts]).as_posix()
            target = str(Path(root / target).resolve().relative_to(root)) if (root / target).exists() else target
            if target in page_labels and name not in page_labels[target]:
                where = f", which is defined in {defined[name][0]}" if name in defined else ""
                problems.append(f"{f}:{n}: {file_part}#{name}: {target} defines no label {name}{where}")

    # The READMEs, whose links keep the file name.
    tracked = subprocess.run(["git", "ls-files", "--", "README.md", "examples/*/README.md"], cwd=root,
                             check=True, capture_output=True, text=True).stdout.split()
    for f in tracked:
        for n, file_part, name in scan(root / f)[1]:
            if not file_part:
                continue
            target = (root / f).parent / file_part
            if not target.exists():
                continue
            target = target.resolve().relative_to(root).as_posix()
            if target in page_labels:
                n_refs += 1
                if name not in page_labels[target]:
                    where = f", which is defined in {defined[name][0]}" if name in defined else ""
                    problems.append(f"{f}:{n}: {file_part}#{name}: {target} defines no label {name}{where}")

    # The sections that the source and the scripts cite.
    n_cited = 0
    files = subprocess.run(["git", "ls-files", "--", *SOURCE_GLOBS], cwd=root, check=True,
                           capture_output=True, text=True).stdout.split("\n")
    for f in sorted(x for x in files if x):
        try:
            text = (root / f).read_text(encoding="utf-8")
        except (UnicodeDecodeError, IsADirectoryError):
            continue
        for n, line in enumerate(text.split("\n"), 1):
            for name in SECTION.findall(line):
                n_cited += 1
                if name not in defined:
                    problems.append(f"{f}:{n}: cites {name}, which the manual does not define")

    if problems:
        print("\n".join(problems))
        print(f"check_labels: {len(problems)} failures")
        return 1
    print(f"check_labels: {len(defined)} labels in {len(scanned)} files of the manual, {n_refs} "
          f"references and {n_cited} citations of a section in the source, each to a label that exists")
    return 0


if __name__ == "__main__":
    sys.exit(main())
