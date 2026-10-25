"""The Python side of Lucifer.

models   the deck as typed dataclasses generated from the schema: Deck, Run, Beam, ...
deck     the flat key = value text of the check scripts as a TOML deck.
schema   the JSON Schema the program writes with lucifer --schema, carried here.
runner   runs the program on a deck and hands back where its output went.
stats    reads the stats file.
beam     reads a particle dump, openPMD or Genesis, into per-slice arrays.

Nothing here holds a default the program does not, generates anything but the TOML the
program reads, or checks what the program refuses itself. The program is the source of
truth, and this package makes decks and reads results.

The names below are reached on first use rather than imported here, so that
python -m lucifer.generate can write models.py before it exists.
"""
from __future__ import annotations

import importlib

_FROM = {name: "lucifer.models" for name in
         ("Deck", "Run", "Window", "Beam", "Resample", "Field", "Wake", "SpaceCharge", "Bmad")}
_FROM.update({name: "lucifer.runner" for name in ("LuciferError", "Result", "run")})

__all__ = sorted(_FROM)


def __getattr__(name):
    if name in _FROM:
        return getattr(importlib.import_module(_FROM[name]), name)
    raise AttributeError(f"module 'lucifer' has no attribute {name!r}")
