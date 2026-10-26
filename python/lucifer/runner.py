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

__all__ = ["run", "Result", "LuciferError", "executable"]


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


@dataclass(repr=False)
class Result:
    """One run: its exit code and output, the input file it ran, and the files it wrote.
    ``stats_file`` names the stats file where it is not ``<out_root>.stats.h5`` in the run's
    directory, as for a run restored from an archive."""
    returncode: int
    output: str
    input_file: pathlib.Path
    cwd: pathlib.Path
    out_root: str
    overrides: list[str] = field(default_factory=list)
    stats_file: pathlib.Path | None = None
    _stats: object = field(default=None, init=False, repr=False, compare=False)

    @property
    def stats_path(self) -> pathlib.Path:
        """The stats file: ``stats_file``, or ``<out_root>.stats.h5`` in the run's directory."""
        if self.stats_file is not None:
            return pathlib.Path(self.stats_file)
        return self.cwd / f"{self.out_root}.stats.h5"

    def __repr__(self):
        return (f"Result(out_root={self.out_root!r}, returncode={self.returncode}, "
                f"stats={str(self.stats_path)!r})")

    def _repr_html_(self):
        from lucifer._display import table
        rows = [("exit", self.returncode), ("out_root", self.out_root), ("directory", self.cwd),
                ("stats file", self.stats_path.name)]
        rows += [("overrides", " ".join(self.overrides))] if self.overrides else []
        rows += [("files", ", ".join(f.name for f in self.files))] if self.output else []
        return table("Result", rows)

    @property
    def stats(self):
        """The stats file, opened through ``lucifer.stats.read_stats`` once and kept open
        until close(). A run cannot overwrite a stats file a reader holds open, so
        Lucifer.run() closes it before running again."""
        from lucifer.stats import read_stats
        if self._stats is None or not self._stats._h5:
            self._stats = read_stats(self.stats_path)
        return self._stats

    def close(self):
        """Close the stats file if ``stats`` opened it."""
        if self._stats is not None:
            self._stats.close()
            self._stats = None

    def plot(self, **kwargs):
        """Plot curves of the run's stats file against s: lucifer.plot.plot, which names
        its arguments. Needs the extra lucifer[plot]."""
        from lucifer.plot import plot
        return plot(self.stats_path, **kwargs)

    @property
    def files(self) -> list[pathlib.Path]:
        """Every file the footer says the run wrote."""
        return [self.cwd / f for f in WROTE.findall(self.output)]


OUT_ROOT = re.compile(r'^ Output\s+out_root "([^"]*)"', re.M)
WROTE = re.compile(r"^(?: Wrote\s+|\s{13})(\S+\.(?:h5|txt))", re.M)
ERROR_BLOCK = re.compile(r"^\[ERROR[^\n]*\n((?:    [^\n]*\n?)+)", re.M)


def run(input_file, cwd=".", exe=None, threads=None, overrides=(), name="lucifer.toml",
        timeout=None) -> Result:
    """Run the program on ``input_file``, a ``Input`` written as ``name`` in ``cwd`` or the path of
    an input file relative to ``cwd``, with ``threads`` as OMP_NUM_THREADS and the overrides
    as ``table.key=value`` arguments. ``timeout`` stops the run after that many seconds.
    Raises ``LuciferError`` on a nonzero exit or a timeout."""
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
    try:
        r = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True, errors="replace",
                           timeout=timeout)
    except subprocess.TimeoutExpired as err:
        out = (err.stdout or b"").decode(errors="replace") if isinstance(err.stdout, bytes) else (err.stdout or "")
        raise LuciferError(f"the run passed its timeout of {timeout} s", out, -1) from None
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
        data = {}
    from lucifer.models import Run
    return str(data.get("run", {}).get("out_root", Run().out_root))
