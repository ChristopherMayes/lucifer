#!/usr/bin/env python3
"""
Check that an edit changed comments and nothing else, and count what the comments say.

  check_comments_only.py --base <rev> <file>...   # code against <rev>, comments ignored
  check_comments_only.py --counts <file>...       # comment counts, one row per file

The first form strips the comments from each file as it is at <rev> and as it is in
the working tree, drops blank lines and trailing blanks, and requires the two to be
identical, so no statement, declaration or continuation can move. Comments are what the
language says they are:

  Fortran   a line whose first nonblank character is !, except an OpenMP sentinel (!$),
            and the rest of a line from an ! outside a quoted string. A string
            continued onto the next line with & carries its quote state there.
  C         /* ... */ across lines, and // to the end of a line, outside strings.
  CMake     # to the end of a line, outside a quoted string.
  Shell     # that begins a word, to the end of the line, outside single and double
            quotes and outside a ${...} expansion, as the shell reads it. So $#, ${#a[@]},
            a quoted "#" and x=#y are code. The shebang is code. The body of a here
            document is the input of a command and is code, comments and all, up to its
            delimiter. A quote left open at the end of a line carries over to the next.
            A .sh file is shell, and so is a file with no suffix whose first line is a
            shebang naming a shell.
  Python    # comments as the tokenizer finds them, and docstrings: the first
            statement of a module, class or function when it is a string. The code
            compared is the syntax tree with the docstrings removed, so a change to any
            statement, literal or message fails.

A docstring is also a value at run time. A script that hands a line of its docstring
to argparse, as __doc__.split("\n")[1] does, must still find a line there that is not
blank, and each such line is printed so it can be read. A shell script's header is
too, where --help prints a range of the script's own lines with sed -n 'A,Bp'. Lines 1
to B+1 must keep the comment rows they had at <rev>, so the help prints the same span
of the header and stops where it stopped, and the range's first and last lines are
printed.

The second form counts, per file and in total: lines, comment lines, comment lines
that narrate (a date, the user, a keystone, a tier, a check section, a mutation, a
commit, a price, a measurement), comment lines holding a semicolon, and lines over
110 columns that carry a comment. A semicolon in a Fortran statement and a long
statement are code and are not counted.

Exit status is nonzero when a file's code differs from <rev>.
"""

from __future__ import annotations

import argparse
import ast
import io
import pathlib
import re
import subprocess
import sys
import tokenize

WIDTH = 110

NARRATION = {
    "date": r"\b20\d\d-\d\d(-\d\d)?\b",
    "the user": r"\bthe user\b",
    "keystone": r"\bkeystone\b",
    "tier": r"\btiers?\b",
    "check section": r"\bcheck sections?\b",
    "mutation": r"\bmutat",
    "commit": r"\bcommit(s|ted)?\b",
    "price": r"\bpric(e|ed|es|ing)\b",
    "measured": r"\bmeasur",
}


def kind(path: str) -> str:
    p = pathlib.Path(path)
    if p.suffix.lower() in (".f90", ".f", ".f95", ".f03"):
        return "fortran"
    if p.suffix.lower() in (".c", ".h", ".cpp", ".mm"):
        return "c"
    if p.suffix.lower() == ".py":
        return "python"
    if p.suffix.lower() == ".sh":
        return "shell"
    if p.suffix == "" and p.is_file():
        with open(p, errors="replace") as f:
            first = f.readline()
        if first.startswith("#!") and re.search(r"\b(ba|z|k|da)?sh\b", first):
            return "shell"
    return "cmake"


def split_fortran(text: str):
    """(code lines, comment strings) of a Fortran source."""
    code, comments = [], []
    quote = None
    for line in text.splitlines():
        stripped = line.lstrip()
        if quote is None and stripped.startswith("!") and not stripped.startswith("!$"):
            comments.append(stripped[1:])
            continue
        out = []
        i = 0
        while i < len(line):
            ch = line[i]
            if quote:
                out.append(ch)
                if ch == quote:
                    quote = None
            elif ch in ("'", '"'):
                quote = ch
                out.append(ch)
            elif ch == "!" and not line[i:].startswith("!$"):
                comments.append(line[i + 1:])
                break
            else:
                out.append(ch)
            i += 1
        # A string that runs off the end of a line continues only when the line ends in &.
        if quote and not "".join(out).rstrip().endswith("&"):
            quote = None
        code.append("".join(out))
    return code, comments


def split_c(text: str):
    comments = []
    out, i, quote, n = [], 0, None, len(text)
    while i < n:
        ch = text[i]
        if quote:
            out.append(ch)
            if ch == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
            out.append(ch)
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            comments.extend(text[i + 2:j - 2].splitlines() or [""])
            out.extend("\n" * text[i:j].count("\n"))
            i = j
            continue
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            comments.append(text[i + 2:j])
            i = j
            continue
        else:
            out.append(ch)
        i += 1
    return "".join(out).splitlines(), comments


def split_cmake(text: str):
    code, comments = [], []
    for line in text.splitlines():
        out, quote = [], False
        for i, ch in enumerate(line):
            if ch == '"':
                quote = not quote
            if ch == "#" and not quote:
                comments.append(line[i + 1:])
                break
            out.append(ch)
        code.append("".join(out))
    return code, comments


# A here document's operator and delimiter: << or <<-, then the word, bare or quoted. A
# here-string, <<<, is not one.
HEREDOC = re.compile(r"<<(-?)[ \t]*(['\"]?)([A-Za-z_][A-Za-z0-9_]*)\2")


def shell_parse(text: str):
    """(code lines, comment strings, comment rows from 1) of a shell script.

    The scan keeps a stack of the contexts it is inside: a single or double quote, a
    command substitution $(...) with its parenthesis depth, and a ${...} expansion. A
    command substitution opens a new context even inside double quotes, so its own
    quotes nest there. The stack carries over from one line to the next."""
    code, comments, rows = [], [], set()
    stack, body = [], []
    for n, line in enumerate(text.splitlines(), 1):
        if body:
            code.append(line)
            delim, tabs = body[0]
            if (line.lstrip("\t") if tabs else line) == delim:
                body.pop(0)
            continue
        if n == 1 and line.startswith("#!"):
            code.append(line)
            continue
        out, i, opened = [], 0, []
        while i < len(line):
            ch = line[i]
            top = stack[-1][0] if stack else None
            if top == "'":
                out.append(ch)
                if ch == "'":
                    stack.pop()
                i += 1
                continue
            if ch == "\\" and i + 1 < len(line):
                out.append(line[i:i + 2])
                i += 2
                continue
            if line.startswith("$(", i) and top != "{":
                stack.append(["(", 1])
                out.append("$(")
                i += 2
                continue
            if line.startswith("${", i):
                stack.append(["{", 0])
                out.append("${")
                i += 2
                continue
            if top == '"':
                out.append(ch)
                if ch == '"':
                    stack.pop()
                i += 1
                continue
            if ch in ("'", '"'):
                stack.append([ch, 0])
            elif ch == "}" and top == "{":
                stack.pop()
            elif ch == "(" and top == "(":
                stack[-1][1] += 1
            elif ch == ")" and top == "(":
                stack[-1][1] -= 1
                if stack[-1][1] == 0:
                    stack.pop()
            elif line.startswith("<<", i) and not line.startswith("<<<", i) \
                    and (i == 0 or line[i - 1] != "<") and top in (None, "("):
                m = HEREDOC.match(line, i)
                if m:
                    opened.append((m.group(3), m.group(1) == "-"))
                    out.append(m.group(0))
                    i = m.end()
                    continue
            elif ch == "#" and top in (None, "(") and (i == 0 or line[i - 1] in " \t;|&()"):
                comments.append(line[i + 1:])
                rows.add(n)
                break
            out.append(ch)
            i += 1
        code.append("".join(out))
        body.extend(opened)
    return code, comments, rows


HELP_RANGE = re.compile(r"sed -n '(\d+),(\d+)p' \"\$\{BASH_SOURCE\[0\]\}\"")


def help_ranges(before: str, after: str):
    """(a, b, rows before, rows after, the range's lines after) for each help range."""
    out = []
    for m in HELP_RANGE.finditer(after):
        a, b = int(m.group(1)), int(m.group(2))
        rows_b = {r for r in shell_parse(before)[2] if r <= b + 1}
        rows_a = {r for r in shell_parse(after)[2] if r <= b + 1}
        out.append((a, b, rows_b, rows_a, after.splitlines()[a - 1:b]))
    return out


def split_shell(text: str):
    code, comments, _ = shell_parse(text)
    return code, comments


DOCSTRING_OWNERS = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def docstring_nodes(tree):
    """The string statements that are docstrings, outermost first."""
    found = []
    for node in ast.walk(tree):
        if isinstance(node, DOCSTRING_OWNERS) and node.body:
            first = node.body[0]
            if (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                    and isinstance(first.value.value, str)):
                found.append(first)
    return found


def python_code(text: str):
    """The syntax tree of a Python source with its docstrings removed, one line per node."""
    tree = ast.parse(text)
    drop = {id(n) for n in docstring_nodes(tree)}
    for node in ast.walk(tree):
        if isinstance(node, DOCSTRING_OWNERS):
            node.body = [n for n in node.body if id(n) not in drop] or [ast.Pass()]
    return ast.dump(tree, indent=1, include_attributes=False).splitlines()


def split_python(text: str):
    comments = [tok.string[1:] for tok in tokenize.generate_tokens(io.StringIO(text).readline)
                if tok.type == tokenize.COMMENT]
    for node in docstring_nodes(ast.parse(text)):
        comments.extend(node.value.value.splitlines() or [""])
    return python_code(text), comments


def python_comment_rows(text: str):
    """Row numbers, from 1, of the lines holding a comment or part of a docstring."""
    rows = {tok.start[0] for tok in tokenize.generate_tokens(io.StringIO(text).readline)
            if tok.type == tokenize.COMMENT}
    for node in docstring_nodes(ast.parse(text)):
        rows.update(range(node.lineno, node.end_lineno + 1))
    return rows


# The ways a script takes one line of its docstring, and the line each takes.
DOC_USES = {
    r'__doc__\.split\("\\n"\)\[(\d+)\]': lambda d, k: d.split("\n")[k],
    r"__doc__\.splitlines\(\)\[(\d+)\]": lambda d, k: d.splitlines()[k],
    r'__doc__\.strip\(\)\.split\("\\n"\)\[(\d+)\]': lambda d, k: d.strip().split("\n")[k],
}


def doc_lines_used(path: str, text: str):
    """(expression, line or None) for each docstring line the script reads at run time."""
    doc = ast.get_docstring(ast.parse(text), clean=False)
    used = []
    for pat, take in DOC_USES.items():
        for m in re.finditer(pat, text):
            try:
                line = take(doc, int(m.group(1))) if doc is not None else None
            except IndexError:
                line = None
            used.append((m.group(0), line))
    return used


def split(path: str, text: str):
    return {"fortran": split_fortran, "c": split_c, "cmake": split_cmake,
            "python": split_python, "shell": split_shell}[kind(path)](text)


def normalized(path: str, text: str):
    code, _ = split(path, text)
    return [line.rstrip() for line in code if line.strip()]


def at_rev(rev: str, path: str) -> str:
    top = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True,
                         text=True, check=True).stdout.strip()
    rel = str(pathlib.Path(path).resolve().relative_to(top))
    return subprocess.run(["git", "-C", top, "show", f"{rev}:{rel}"], capture_output=True,
                          text=True, check=True).stdout


def compare(rev: str, paths) -> int:
    bad = 0
    for path in paths:
        text = pathlib.Path(path).read_text()
        before = normalized(path, at_rev(rev, path))
        after = normalized(path, text)
        if kind(path) == "python":
            blank = False
            for expr, line in doc_lines_used(path, text):
                if line is None or not line.strip():
                    blank = True
                    print(f"FAIL  {path}: {expr} is blank or missing")
                else:
                    print(f"      {path}: {expr} is {line.strip()!r}")
            if blank:
                bad += 1
                continue
        if kind(path) == "shell":
            moved = False
            for a, b, rows_b, rows_a, shown in help_ranges(at_rev(rev, path), text):
                if rows_b != rows_a:
                    moved = True
                    print(f"FAIL  {path}: --help prints lines {a} to {b}, and the comment rows "
                          f"up to {b + 1} moved")
                else:
                    print(f"      {path}: --help prints lines {a} to {b}, from "
                          f"{shown[0][:60]!r} to {shown[-1][:60]!r}")
            if moved:
                bad += 1
                continue
        if before == after:
            unit = "syntax tree lines" if kind(path) == "python" else "code lines"
            print(f"ok    {path}: code identical to {rev}, {len(after)} {unit}")
            continue
        bad += 1
        print(f"FAIL  {path}: code differs from {rev}")
        for k, (a, b) in enumerate(zip(before, after)):
            if a != b:
                print(f"      first difference at code line {k + 1}:\n        was: {a}\n        now: {b}")
                break
        else:
            print(f"      code line count {len(before)} became {len(after)}")
    return bad


def counts(paths) -> None:
    keys = list(NARRATION)
    head = f"{'file':40s} {'lines':>6s} {'comm':>5s} {'%':>3s} " + " ".join(
        f"{k[:5]:>5s}" for k in keys) + f" {'semi':>5s} {'>110':>5s}"
    print(head)
    tot = dict.fromkeys(["lines", "comm", "semi", "wide", *keys], 0)
    for path in paths:
        text = pathlib.Path(path).read_text()
        lines = text.splitlines()
        _, comments = split(path, text)
        # A long statement is code. Only a line that carries a comment counts as wide.
        if kind(path) in ("python", "shell"):
            rows = python_comment_rows(text) if kind(path) == "python" else shell_parse(text)[2]
            wide = sum(1 for k, line in enumerate(lines, 1) if len(line) > WIDTH and k in rows)
        else:
            wide = sum(1 for line in lines if len(line) > WIDTH and split(path, line)[1])
        row = {"lines": len(lines), "comm": len(comments),
               "semi": sum(1 for c in comments if ";" in c), "wide": wide}
        for k, pat in NARRATION.items():
            row[k] = sum(1 for c in comments if re.search(pat, c, re.I))
        for k in tot:
            tot[k] += row[k]
        pct = 100 * row["comm"] // max(row["lines"], 1)
        print(f"{pathlib.Path(path).name:40s} {row['lines']:6d} {row['comm']:5d} {pct:3d} "
              + " ".join(f"{row[k]:5d}" for k in keys) + f" {row['semi']:5d} {row['wide']:5d}")
    pct = 100 * tot["comm"] // max(tot["lines"], 1)
    print(f"{'total':40s} {tot['lines']:6d} {tot['comm']:5d} {pct:3d} "
          + " ".join(f"{tot[k]:5d}" for k in keys) + f" {tot['semi']:5d} {tot['wide']:5d}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--base", help="Revision whose code the files must still have")
    ap.add_argument("--counts", action="store_true", help="Count the comments instead")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    if a.counts:
        counts(a.files)
        return 0
    if not a.base:
        ap.error("give --base <rev> or --counts")
    return 1 if compare(a.base, a.files) else 0


if __name__ == "__main__":
    sys.exit(main())
