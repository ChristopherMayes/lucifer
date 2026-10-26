#!/usr/bin/env python3
"""
Check that the notebooks under docs/notebooks read as practical examples.

  check_notebooks.py [--root <repository root>]

A notebook shows the package in use one step at a time, with the output of each step
stored beneath it. tests/scripts/run_notebooks.py executes them, and the site renders
what they store. Refused, each with the notebook and the cell:

  - a code cell over six lines;
  - a code cell that calls print, since a notebook shows a value rather than prints it;
  - a code cell whose last statement is an expression, not ended by a semicolon, with
    no stored output, which is a notebook not executed since that cell was written;
  - an output holding an error;
  - execution counts that are not 1, 2, 3 and on in cell order;
  - a markdown cell over three sentences, headings aside.

Standard library only, so the docs job runs it without an environment.
"""
import argparse
import ast
import json
import re
import sys
from pathlib import Path

MAX_CODE_LINES = 6
MAX_SENTENCES = 3
PRINT = re.compile(r"(?<![\w.])print\s*\(")
SENTENCE_END = re.compile(r"[.!?](?=\s|$)")


def source(cell):
    src = cell.get("source", "")
    return "".join(src) if isinstance(src, list) else src


def returns_value(code):
    """Whether the cell's last statement is an expression IPython would display."""
    if code.rstrip().endswith(";"):
        return False
    try:
        tree = ast.parse(code)
    except SyntaxError:
        return False
    return bool(tree.body) and isinstance(tree.body[-1], ast.Expr)


def sentences(text):
    prose = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    # Code spans and numbers with a decimal point do not end a sentence.
    prose = re.sub(r"`[^`]*`", "x", prose)
    prose = re.sub(r"\d\.\d", "0", prose)
    return len(SENTENCE_END.findall(prose))


def check(path):
    problems = []
    book = json.loads(path.read_text())
    count = 0
    for i, cell in enumerate(book.get("cells", []), 1):
        text = source(cell)
        where = f"{path.name} cell {i}"
        if cell.get("cell_type") == "markdown":
            n = sentences(text)
            if n > MAX_SENTENCES:
                problems.append(f"{where}: a markdown cell of {n} sentences, over {MAX_SENTENCES}")
            continue
        if cell.get("cell_type") != "code":
            continue
        lines = [line for line in text.splitlines() if line.strip()]
        if len(lines) > MAX_CODE_LINES:
            problems.append(f"{where}: {len(lines)} lines of code, over {MAX_CODE_LINES}")
        if PRINT.search(text):
            problems.append(f"{where}: calls print; show the value instead")
        outputs = cell.get("outputs", [])
        if any(o.get("output_type") == "error" for o in outputs):
            problems.append(f"{where}: an output holds an error")
        if returns_value(text) and not outputs:
            problems.append(f"{where}: returns a value and stores no output; execute the notebook")
        count += 1
        if cell.get("execution_count") != count:
            problems.append(f"{where}: execution count {cell.get('execution_count')}, expected {count}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    root = Path(ap.parse_args().root).resolve()
    books = sorted((root / "docs" / "notebooks").glob("*.ipynb"))
    problems = [p for b in books for p in check(b)]
    if not books:
        problems.append("no notebooks under docs/notebooks")
    if problems:
        print("\n".join(problems))
        print(f"check_notebooks: FAIL, {len(problems)} problems")
        return 1
    print(f"check_notebooks: {len(books)} notebooks, every cell short, shown rather than printed, "
          f"and executed in order")
    return 0


if __name__ == "__main__":
    sys.exit(main())
