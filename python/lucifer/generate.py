"""Write ``models.py`` from the schema: one dataclass per table and ``Deck``.

    python -m lucifer.generate

The keystone runs this after ``lucifer --schema`` and requires no diff, so the models
cannot disagree with the program. Every field carries the program's type and default
and a docstring of the schema's description, which ends with the unit where the key
has one, and the bounds. The dataclasses hold data alone: what a ``Deck`` does is in
``lucifer.deckops``, which the generated ``Deck`` inherits.
"""
from __future__ import annotations

import keyword
import pathlib

from lucifer import schema as schema_mod

CLASS = {"run": "Run", "window": "Window", "beam": "Beam", "beam.resample": "Resample",
         "field": "Field", "wake": "Wake", "space_charge": "SpaceCharge", "bmad": "Bmad"}
PY_TYPE = {"string": "str", "number": "float", "integer": "int", "boolean": "bool"}


def py_type(prop: dict) -> str:
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
    keys = schema_mod.keys(schema)
    out = [
        '"""The deck as typed models, generated from the schema. Do not edit.',
        "",
        f"Written by python -m lucifer.generate from schema.json, schema version "
        f"{schema['x-schema-version']}. A dataclass per table of the deck, and Deck with",
        "the lattice and one of each. The defaults are the program's, read from the schema,",
        "and nothing here checks a bound: the program refuses a value outside one. The",
        "field named field is the [field] table, so dataclasses.field is imported as _field.",
        '"""',
        "from __future__ import annotations",
        "",
        "from dataclasses import dataclass, field as _field",
        "",
        "from lucifer.deckops import DeckOps",
        "",
        "",
    ]
    # The tables, the nested one before the table that holds it.
    order = ["run", "window", "beam.resample", "beam", "field", "wake", "space_charge", "bmad"]
    for table in order:
        out += ["@dataclass", f"class {CLASS[table]}:", f'    """The [{table}] table of the deck."""']
        for name, prop in keys[table].items():
            out += field_lines(name, prop)
        if table == "beam":
            out += ["    resample: Resample = _field(default_factory=Resample)",
                    '    """The [beam.resample] table, the resampler."""']
        out += ["", ""]
    out += ["@dataclass", "class Deck(DeckOps):", '    """One deck: the lattice and the tables."""']
    for name, prop in keys[""].items():
        out += field_lines(name, prop)
    for table in schema_mod.TABLES:
        if "." in table:
            continue
        out += [f"    {table}: {CLASS[table]} = _field(default_factory={CLASS[table]})",
                f'    """The [{table}] table."""']
    out += ["", "", "TABLES = (" + ", ".join(f'("{t}", {CLASS[t]})' for t in schema_mod.TABLES) + ")",
            '"""The tables in the deck\'s order, with their classes."""', ""]
    return "\n".join(out)


def main(argv=None) -> int:
    text = render(schema_mod.load())
    pathlib.Path(__file__).with_name("models.py").write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
