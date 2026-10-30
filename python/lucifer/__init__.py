"""The Python side of Lucifer.

lucifer  Lucifer, one object to build, run, read and plot a run, with lume-base's names.
models   the input file as typed dataclasses generated from the schema: Input, Run, Beam, ...
flat     the flat key = value text of the check scripts as a TOML input file.
schema   the JSON Schema the program writes with lucifer --schema, carried here.
runner   runs the program on an input file and hands back where its output went.
dumps    a run's particle and field dumps, as files and as objects read on first use.
stats    reads the stats file.
beam     reads a particle dump, openPMD or Genesis, into per-slice arrays.
field    reads a field dump, openPMD or Genesis, and the pulse and the escaped-field bank.
compare  reads a Lucifer run and a reference run, Genesis4's or another Lucifer run, into
         one form and takes their differences.
plot     plots a stats file against s with the lattice layout, its summary figure, and a
         run over a reference run with the difference beneath each curve.

Nothing here holds a default the program does not, generates anything but the TOML the
program reads, or checks what the program refuses itself. The program is the source of
truth, and this package makes input files and reads results.

The names below are reached on first use rather than imported here, so that
python -m lucifer.generate can write models.py before it exists.
"""
from __future__ import annotations

import importlib
from importlib.metadata import PackageNotFoundError, version

# The version setuptools-scm gave the installed package, from the git tag the program's
# version comes from. An editable install keeps the version it was installed at.
try:
    __version__ = version("lucifer")
except PackageNotFoundError:
    __version__ = "0.0.0"

_FROM = {name: "lucifer.models" for name in
         ("Input", "Run", "Window", "Beam", "Resample", "Field", "Wake", "SpaceCharge", "Bmad")}
_FROM.update({name: "lucifer.runner" for name in ("LuciferError", "Result", "run")})
_FROM["Lucifer"] = "lucifer.lucifer"

__all__ = sorted(_FROM)


def __getattr__(name):
    if name in _FROM:
        return getattr(importlib.import_module(_FROM[name]), name)
    raise AttributeError(f"module 'lucifer' has no attribute {name!r}")
