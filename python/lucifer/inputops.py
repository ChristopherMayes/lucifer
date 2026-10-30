"""What a ``Input`` does: render itself as TOML, read one, and reach a key by its dotted
name. ``models.Input`` inherits this, so the generated file holds data alone.

The TOML written here is the program's: a key is written when its value differs from
the program's default, every key with ``explicit=True``, and nothing is checked that
the program refuses itself.
"""
from __future__ import annotations

import dataclasses
import pathlib
import tomllib

__all__ = ["InputOps", "toml_value"]


def toml_value(v) -> str:
    """One value as TOML text."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    raise TypeError(f"no TOML for {type(v).__name__}: {v!r}")


def _default(f: dataclasses.Field):
    if f.default is not dataclasses.MISSING:
        return f.default
    return f.default_factory()


def _is_table(f: dataclasses.Field) -> bool:
    # A table's annotation is its class name, Run or Resample. A key with a set of names is
    # annotated Literal[...], which also starts upper case and is not a name.
    return dataclasses.is_dataclass(f.type) or (isinstance(f.type, str) and f.type.isidentifier()
                                                and f.type[:1].isupper())


class InputOps:
    """Mixin for the generated ``Input``."""

    # ------------------------------------------------------------------
    # TOML out.

    def toml(self, explicit: bool = False) -> str:
        """The input file as TOML: the keys set away from the program's defaults, or every
        key with ``explicit``. The lattice is always written."""
        lines = []
        for f in dataclasses.fields(self):
            if _is_table(f):
                continue
            v = getattr(self, f.name)
            if explicit or f.name == "lattice" or v != _default(f):
                lines.append(f"{f.name} = {toml_value(v)}")
        for f in dataclasses.fields(self):
            if _is_table(f):
                lines += _table_lines(f.name, getattr(self, f.name), explicit)
        return "\n".join(lines) + "\n"

    def write(self, path, explicit: bool = False) -> pathlib.Path:
        """Write the input file at ``path`` and return it."""
        path = pathlib.Path(path)
        path.write_text(self.toml(explicit))
        return path

    # ------------------------------------------------------------------
    # TOML in.

    @classmethod
    def from_dict(cls, data: dict):
        """An input file from the mapping a TOML parser returns. An unknown key is Python's own
        TypeError from the dataclass, and a bound is the program's."""
        kw = {}
        for f in dataclasses.fields(cls):
            if f.name not in data:
                continue
            kw[f.name] = _table_from(f, data[f.name]) if _is_table(f) else data[f.name]
        extra = set(data) - {f.name for f in dataclasses.fields(cls)}
        if extra:
            raise TypeError(f"the input file has no key {sorted(extra)}")
        return cls(**kw)

    @classmethod
    def loads(cls, text: str):
        """An input file from TOML text."""
        return cls.from_dict(tomllib.loads(text))

    @classmethod
    def read(cls, path):
        """An input file from a TOML file."""
        return cls.loads(pathlib.Path(path).read_text())

    @classmethod
    def from_flat(cls, text: str):
        """An input file from flat ``key = value`` text, through ``lucifer.flat.to_toml``."""
        from lucifer.flat import to_toml
        return cls.loads(to_toml(text))

    # ------------------------------------------------------------------
    # What it shows: the lattice and the keys set away from the defaults.

    def __repr__(self):
        keys = ", ".join(f"{k}={v!r}" for k, v in self.items())
        return f"{type(self).__name__}({keys})"

    def _repr_html_(self):
        from lucifer._display import short, table
        return table(type(self).__name__, [(k, short(v)) for k, v in self.items()])

    # ------------------------------------------------------------------
    # Dotted access: inp["beam.n_particle"].

    def _resolve(self, key: str):
        obj = self
        parts = key.split(".")
        for p in parts[:-1]:
            obj = getattr(obj, p)
        if not hasattr(obj, parts[-1]):
            raise KeyError(key)
        return obj, parts[-1]

    def __getitem__(self, key: str):
        """The value of a key by its dotted name, ``inp["beam.n_particle"]``. A name the
        input file does not have raises KeyError."""
        obj, name = self._resolve(key)
        return getattr(obj, name)

    def __setitem__(self, key: str, value) -> None:
        """Set a key by its dotted name, ``inp["beam.n_particle"] = 4096``."""
        obj, name = self._resolve(key)
        setattr(obj, name, value)

    def items(self, explicit: bool = False):
        """``(dotted key, value)`` pairs, the keys set away from the defaults or every key."""
        out = []
        for f in dataclasses.fields(self):
            if _is_table(f):
                out += _table_items(f.name, getattr(self, f.name), explicit)
            elif explicit or f.name == "lattice" or getattr(self, f.name) != _default(f):
                out.append((f.name, getattr(self, f.name)))
        return out


def _table_from(f: dataclasses.Field, data: dict):
    cls = f.default_factory
    kw = {}
    for sf in dataclasses.fields(cls):
        if sf.name in data:
            kw[sf.name] = _table_from(sf, data[sf.name]) if _is_table(sf) else data[sf.name]
    extra = set(data) - {sf.name for sf in dataclasses.fields(cls)}
    if extra:
        raise TypeError(f"[{f.name}] has no key {sorted(extra)}")
    return cls(**kw)


def _table_items(table: str, obj, explicit: bool) -> list:
    out = []
    for f in dataclasses.fields(obj):
        v = getattr(obj, f.name)
        if _is_table(f):
            out += _table_items(f"{table}.{f.name}", v, explicit)
        elif explicit or v != _default(f):
            out.append((f"{table}.{f.name}", v))
    return out


def _table_lines(table: str, obj, explicit: bool) -> list[str]:
    own = [(f.name, getattr(obj, f.name)) for f in dataclasses.fields(obj)
           if not _is_table(f) and (explicit or getattr(obj, f.name) != _default(f))]
    lines = []
    if own:
        lines += ["", f"[{table}]"] + [f"{k} = {toml_value(v)}" for k, v in own]
    for f in dataclasses.fields(obj):
        if _is_table(f):
            lines += _table_lines(f"{table}.{f.name}", getattr(obj, f.name), explicit)
    return lines
