#!/usr/bin/env python3
"""
Check the input reference against the schema the program writes.

The deck's keys, their types, defaults and bounds are the program's, written by
lucifer --schema into python/lucifer/schema.json and rendered into
docs/generated/input-keys.md by report_input_keys.py, which docs/input-reference.md
includes. What this check holds together:

  1. python/lucifer/deck.py's routing table names exactly the schema's keys, table by
     table, so a deck the package writes is one the program reads, and the generated
     models python/lucifer/models.py carry exactly those keys with the schema's defaults.
  2. Every key of the schema appears in docs/generated/input-keys.md.
  3. Every dotted key a page names in code font, `beam.n_particle` or `[wake]`, is a key
     or a table of the schema, so no page describes a key the program does not read.
     The pages are the committed Markdown under docs/ and the READMEs, less the
     changelog, whose entries are history.

Standard library only, so the docs job runs it without an environment.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

DOTTED = re.compile(r"`((?:run|window|beam\.resample|beam|field|wake|space_charge|bmad)\.[a-z_]+)`")
# A file named in code font, `run.sh` or `bmad.mod`, reads like a dotted key and is not one.
FILE_SUFFIXES = {"sh", "py", "md", "mod", "h5", "toml", "in", "bmad", "lat", "txt", "json", "yml", "cmake"}
TABLE = re.compile(r"`\[((?:run|window|beam\.resample|beam|field|wake|space_charge|bmad|[a-z_.]+))\]`")


def schema_keys(schema):
    keys = {"": set()}
    for name, p in schema["properties"].items():
        if p.get("type") != "object":
            keys[""].add(name)
            continue
        keys[name] = set()
        for k, kp in p["properties"].items():
            if kp.get("type") == "object":
                keys[f"{name}.{k}"] = set(kp["properties"])
            else:
                keys[name].add(k)
    return keys


def schema_table(schema, table):
    props = schema["properties"]
    for part in table.split("."):
        props = props[part]["properties"]
    return props


def main():
    failed = []
    schema = json.loads((ROOT / "python/lucifer/schema.json").read_text())
    keys = schema_keys(schema)

    # (1) The package's routing table, and the models' fields and defaults.
    from lucifer.deck import KEYS
    from lucifer import models
    for table, names in keys.items():
        cls = models.Deck if table == "" else dict(models.TABLES)[table]
        obj = cls()
        have = {k for k in vars(obj) if not hasattr(vars(obj)[k], "__dataclass_fields__")}
        if have != names:
            failed.append(f"python/lucifer/models.py {cls.__name__} differs from the schema: "
                          f"missing {sorted(names - have)}, extra {sorted(have - names)}")
        props = schema["properties"] if table == "" else schema_table(schema, table)
        for k in names & have:
            if getattr(obj, k) != props[k]["default"]:
                failed.append(f"python/lucifer/models.py {cls.__name__}.{k} defaults to "
                              f"{getattr(obj, k)!r}, the schema to {props[k]['default']!r}")
    for table, names in keys.items():
        mine = set(KEYS.get(table, []))
        if mine != names:
            failed.append(f"python/lucifer/deck.py KEYS[{table!r}] differs from the schema: "
                          f"missing {sorted(names - mine)}, extra {sorted(mine - names)}")
    for table in KEYS:
        if table not in keys:
            failed.append(f"python/lucifer/deck.py names a table the schema does not: {table!r}")

    # (2) The generated key tables.
    page = (ROOT / "docs/generated/input-keys.md").read_text()
    for table, names in keys.items():
        for k in names:
            if f"| `{k}` |" not in page:
                failed.append(f"docs/generated/input-keys.md does not list {table}.{k}")

    # (3) Dotted keys on the pages.
    tracked = subprocess.run(["git", "ls-files", "--", "*.md"], cwd=ROOT, capture_output=True,
                             text=True).stdout.split()
    tables = set(keys)
    for f in tracked:
        if f.startswith("docs/generated/") or f == "docs/changelog.md":
            continue
        text = (ROOT / f).read_text()
        for m in DOTTED.finditer(text):
            table, _, name = m.group(1).rpartition(".")
            if name in FILE_SUFFIXES:
                continue
            if table not in keys or name not in keys[table]:
                failed.append(f"{f} names `{m.group(1)}`, which the schema does not have")
        for m in TABLE.finditer(text):
            if m.group(1) not in tables - {""}:
                failed.append(f"{f} names `[{m.group(1)}]`, which is not a table of the deck")

    n = sum(len(v) for v in keys.values())
    if failed:
        print("\n".join(failed))
        print(f"check_input_reference: FAIL, {len(failed)} problems over {n} keys")
        return 1
    print(f"check_input_reference: {n} keys in {len(keys)} tables, the package, the key tables and the pages agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
