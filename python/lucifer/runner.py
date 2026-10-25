"""Run the program on an input file and hand back where its output went.

    from lucifer import Input, run
    inp = Input(lattice="aramis.bmad")
    inp["beam.n_particle"] = 4096
    result = run(inp, cwd="scan7", threads=12)
    result.stats["field/total/power"]

    for n in (1024, 2048, 4096):
        run("lucifer.toml", overrides=[f"beam.n_particle={n}", f"run.out_root=n{n}"])

The executable is ``exe``, else ``$LUCIFER_EXE``, else ``lucifer`` on PATH. A nonzero
exit raises ``LuciferError`` carrying the program's message, so a refusal reads as the
program wrote it.
"""
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
from dataclasses import dataclass, field


class LuciferError(RuntimeError):
    """The program exited nonzero. ``message`` is its error text and ``output`` all of it."""

    def __init__(self, message: str, output: str = "", returncode: int = 1):
        super().__init__(message)
        self.message = message
        self.output = output
        self.returncode = returncode


def executable(exe=None) -> pathlib.Path:
    """The program: ``exe``, ``$LUCIFER_EXE`` or ``lucifer`` on PATH."""
    cand = exe or os.environ.get("LUCIFER_EXE") or shutil.which("lucifer")
    if not cand:
        raise LuciferError("no lucifer executable: pass exe=, set LUCIFER_EXE or put lucifer on PATH")
    p = pathlib.Path(cand)
    if not p.is_file():
        raise LuciferError(f"no lucifer executable at {p}")
    return p.resolve()


@dataclass
class Result:
    """One run: its exit code and output, the input file it ran, and the files it wrote."""
    returncode: int
    output: str
    input_file: pathlib.Path
    cwd: pathlib.Path
    out_root: str
    overrides: list[str] = field(default_factory=list)

    @property
    def stats_path(self) -> pathlib.Path:
        """``<out_root>.stats.h5`` in the run's directory."""
        return self.cwd / f"{self.out_root}.stats.h5"

    @property
    def stats(self):
        """The stats file, opened through ``lucifer.stats.read_stats``."""
        from lucifer.stats import read_stats
        return read_stats(self.stats_path)

    @property
    def files(self) -> list[pathlib.Path]:
        """Every file the footer says the run wrote."""
        return [self.cwd / f for f in WROTE.findall(self.output)]


OUT_ROOT = re.compile(r'^ Output\s+out_root "([^"]*)"', re.M)
WROTE = re.compile(r"^(?: Wrote\s+|\s{13})(\S+\.(?:h5|txt))", re.M)
ERROR_BLOCK = re.compile(r"^\[ERROR[^\n]*\n((?:    [^\n]*\n?)+)", re.M)


def run(input_file, cwd=".", exe=None, threads=None, overrides=(), name="lucifer.toml") -> Result:
    """Run the program on ``input_file``, a ``Input`` written as ``name`` in ``cwd`` or the path of
    an input file relative to ``cwd``, with ``threads`` as OMP_NUM_THREADS and the overrides
    as ``table.key=value`` arguments. Raises ``LuciferError`` on a nonzero exit."""
    cwd = pathlib.Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    if hasattr(input_file, "toml"):
        input_path = input_file.write(cwd / name)
        out_root = input_file.run.out_root
    else:
        input_path = cwd / input_file
        out_root = _out_root_of(input_path)
    env = dict(os.environ)
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(int(threads))
    cmd = [str(executable(exe)), str(input_path.relative_to(cwd) if input_path.is_relative_to(cwd) else input_path),
           *overrides]
    r = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True, errors="replace")
    output = r.stdout + r.stderr
    m = OUT_ROOT.search(output)
    if m:
        out_root = m.group(1)
    if r.returncode != 0:
        blocks = ERROR_BLOCK.findall(output)
        message = "\n".join(b.strip() for b in blocks) or output.strip()[-2000:]
        raise LuciferError(message, output, r.returncode)
    return Result(r.returncode, output, input_path, cwd.resolve(), out_root, list(overrides))


def _out_root_of(input_path: pathlib.Path) -> str:
    import tomllib
    try:
        data = tomllib.loads(input_path.read_text())
    except (OSError, tomllib.TOMLDecodeError):
        return "fel_td"
    return str(data.get("run", {}).get("out_root", "fel_td"))
