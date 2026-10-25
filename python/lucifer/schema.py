"""The schema the program writes with ``lucifer --schema``, carried by the package.

``schema.json`` beside this module is the one committed copy. The keystone writes it from
the production executable and requires no diff, so the keys, the types, the defaults and
the bounds here are the program's own. ``python -m lucifer.schema`` prints its path, for
an editor that validates a deck against JSON Schema.
"""
from __future__ import annotations

import json
import pathlib

TABLES = ("run", "window", "beam", "beam.resample", "field", "wake", "space_charge", "bmad")


def path() -> pathlib.Path:
    """The schema file."""
    return pathlib.Path(__file__).with_name("schema.json")


def load() -> dict:
    """The schema as the JSON Schema document the program wrote."""
    return json.loads(path().read_text())


def keys(schema: dict | None = None) -> dict[str, dict[str, dict]]:
    """``{table: {name: property}}``, the top level under ``""`` and the nested table
    ``beam.resample`` flattened to its dotted name, in the schema's order."""
    schema = load() if schema is None else schema
    out: dict[str, dict[str, dict]] = {"": {}}
    for name, prop in schema["properties"].items():
        if prop.get("type") != "object":
            out[""][name] = prop
            continue
        out[name] = {}
        for k, kp in prop["properties"].items():
            if kp.get("type") == "object":
                out[f"{name}.{k}"] = dict(kp["properties"])
            else:
                out[name][k] = kp
    return out


def main(argv=None) -> int:
    print(path())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
