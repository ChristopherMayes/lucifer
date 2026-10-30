"""Write ``models.py`` from the schema: one dataclass per table and ``Input``.

    python -m lucifer.generate

The keystone runs this after ``lucifer --schema`` and requires no diff, so the models
cannot disagree with the program. Every field carries the program's type and default
and a docstring of the schema's description, which ends with the unit where the key
has one, and the bounds. A key with a set of names, an enum in the schema, is typed as a
Literal of the set, which a type checker reads and nothing enforces at run time. The
dataclasses hold data alone: what a ``Input`` does is in ``lucifer.inputops``, which the
generated ``Input`` inherits.
"""
from __future__ import annotations

import keyword
import pathlib

from lucifer import schema as schema_mod

__all__ = ["render", "main"]

CLASS = {"run": "Run", "window": "Window", "beam": "Beam", "beam.resample": "Resample",
         "field": "Field", "wake": "Wake", "space_charge": "SpaceCharge", "bmad": "Bmad"}
PY_TYPE = {"string": "str", "number": "float", "integer": "int", "boolean": "bool"}


def py_type(prop: dict) -> str:
    if "enum" in prop:
        return "Literal[" + ", ".join(repr(v) for v in prop["enum"]) + "]"
    if prop["type"] == "array":
        return f"list[{PY_TYPE[prop['items']['type']]}]"
    return PY_TYPE[prop["type"]]


def py_default(prop: dict) -> str:
    v = prop["default"]
    if prop["type"] == "array":
        return f"_field(default_factory=lambda: {v!r})"
    if prop["type"] == "number":
        return repr(float(v))
    return repr(v)


def bounds_text(prop: dict) -> str:
    parts = []
    if "minimum" in prop:
        parts.append(f"from {prop['minimum']}")
    if "exclusiveMinimum" in prop:
        parts.append(f"above {prop['exclusiveMinimum']}")
    if "maximum" in prop:
        parts.append(f"to {prop['maximum']}")
    if "exclusiveMaximum" in prop:
        parts.append(f"below {prop['exclusiveMaximum']}")
    return " ".join(parts)


def field_lines(name: str, prop: dict) -> list[str]:
    assert name.isidentifier() and not keyword.iskeyword(name), name
    doc = prop["description"].rstrip()
    b = bounds_text(prop)
    if b:
        doc += f" Bounds: {b}."
    return [f"    {name}: {py_type(prop)} = {py_default(prop)}",
            f'    """{doc}"""']


def render(schema: dict) -> str:
    """The text of models.py for ``schema``, the document ``lucifer.schema.load`` returns."""
    keys = schema_mod.keys(schema)
    out = [
        '"""The input file as typed models, generated from the schema. Do not edit.',
        "",
        f"Written by python -m lucifer.generate from schema.json, schema version "
        f"{schema['x-schema-version']}. A dataclass per table of the input file, and Input with",
        "the lattice and one of each. The defaults are the program's, read from the schema,",
        "and nothing here checks a bound or a set of names: the program refuses a value outside",
        "one. The field named field is the [field] table, so dataclasses.field is imported as",
        "_field.",
        '"""',
        "from __future__ import annotations",
        "",
        "from dataclasses import dataclass, field as _field",
        "from typing import Literal",
        "",
        "from lucifer.inputops import InputOps",
        "",
        "__all__ = [" + ", ".join(f'"{c}"' for c in ["Input", *(CLASS[t] for t in schema_mod.TABLES)]) + "]",
        "",
        "",
    ]
    # The tables, the nested one before the table that holds it.
    order = ["run", "window", "beam.resample", "beam", "field", "wake", "space_charge", "bmad"]
    for table in order:
        out += ["@dataclass", f"class {CLASS[table]}:", f'    """The [{table}] table of the input file."""']
        for name, prop in keys[table].items():
            out += field_lines(name, prop)
        if table == "beam":
            out += ["    resample: Resample = _field(default_factory=Resample)",
                    '    """The [beam.resample] table, the resampler."""']
        out += ["", ""]
    # repr=False, so the representation is InputOps's: the keys set away from the defaults.
    out += ["@dataclass(repr=False)", "class Input(InputOps):", '    """One input file: the lattice and the tables."""']
    for name, prop in keys[""].items():
        out += field_lines(name, prop)
    for table in schema_mod.TABLES:
        if "." in table:
            continue
        out += [f"    {table}: {CLASS[table]} = _field(default_factory={CLASS[table]})",
                f'    """The [{table}] table."""']
    out += ["", "", "TABLES = (" + ", ".join(f'("{t}", {CLASS[t]})' for t in schema_mod.TABLES) + ")",
            '"""The tables in the input file\'s order, with their classes."""', ""]
    return "\n".join(out)


def main(argv=None) -> int:
    """python -m lucifer.generate: write models.py beside this module from schema.json."""
    text = render(schema_mod.load())
    pathlib.Path(__file__).with_name("models.py").write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
