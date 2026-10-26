#!/usr/bin/env python3
"""
Execute the notebooks under docs/notebooks and write them with their outputs.

  run_notebooks.py --out <dir> [--exe <lucifer>]

Each notebook runs on the steady-state example from a scratch copy of examples/
steady_state with the shared lattice beside it, so the notebooks write nothing into
examples/ and their relative paths read as they would for a user who copied the
directory. The stored outputs are dropped first, the notebook is executed with nbclient
on the python3 kernel, and the executed notebook is written to <out> under its own name.
A cell that raises stops the script with the notebook, the cell's source and the error.
LUCIFER_EXE is set from --exe, else kept as the environment has it.

examples/run_examples.sh runs this after the examples, into docs/notebooks with figures
on and into its scratch directory under --no-figures, so every keystone executes the
notebooks while their stored outputs, like the committed figures, stay out of the
regeneration diff.
"""
import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
NOTEBOOKS = ROOT / "docs" / "notebooks"
EXAMPLE = "steady_state"
TIMEOUT = 600


def scratch_copy(dest):
    """A copy of the example's tracked files, with the shared lattice its input file names
    by a relative path, at the same relative place."""
    files = subprocess.run(["git", "ls-files", f"examples/{EXAMPLE}", "examples/aramis.bmad"],
                           cwd=ROOT, check=True, capture_output=True, text=True).stdout.split()
    for f in files:
        target = dest / pathlib.Path(f).relative_to("examples")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / f, target)
    return dest / EXAMPLE


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--out", required=True, help="where the executed notebooks are written")
    ap.add_argument("--exe", help="the lucifer executable, set as LUCIFER_EXE")
    args = ap.parse_args()
    import nbformat
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    if args.exe:
        os.environ["LUCIFER_EXE"] = str(pathlib.Path(args.exe).resolve())
    out = pathlib.Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    books = sorted(NOTEBOOKS.glob("*.ipynb"))
    if not books:
        print(f"run_notebooks: no notebooks under {NOTEBOOKS}")
        return 1
    for path in books:
        book = nbformat.read(path, as_version=4)
        for cell in book.cells:
            if cell.cell_type == "code":
                cell.outputs, cell.execution_count = [], None
                cell.metadata.pop("execution", None)
        with tempfile.TemporaryDirectory(prefix="lucifer_notebook_") as tmp:
            cwd = scratch_copy(pathlib.Path(tmp))
            client = NotebookClient(book, timeout=TIMEOUT, kernel_name="python3", record_timing=False,
                                    resources={"metadata": {"path": str(cwd)}})
            try:
                client.execute()
            except CellExecutionError as err:
                print(f"run_notebooks: {path.name} failed:\n{err}")
                return 1
        nbformat.write(book, out / path.name)
        n = sum(c.cell_type == "code" for c in book.cells)
        print(f"run_notebooks: {path.name}, {n} code cells executed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
