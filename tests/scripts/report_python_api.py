#!/usr/bin/env python3
"""
Generate docs/generated/python-api.md from the Python package's own docstrings.

The Fortran reference, docs/generated/api.md, is written by report_api.py from the
routine headers. This is the same page for python/lucifer: every public module, class
and function, with its signature from inspect and its docstring as written, so the
reference says what the code says. A module's public names are its __all__. Nothing
here is written by hand, and the keystone suite and the docs workflow regenerate the
page and require an empty diff.

The generated models are listed by name and not field by field. Their fields are the
input keys, which docs/generated/input-keys.md lists with their meaning, so the page
points there rather than repeating 106 docstrings.

Refusals, so the page cannot fall behind the package silently:
  - a module of the package that ORDER does not name, or a name in ORDER with no module;
  - a module without __all__, or a name in it the module does not define;
  - a public module, class, function, method or property without a docstring.

Deterministic by construction: modules in ORDER, names in __all__ order, members in
definition order, no timestamp and no path outside the tree.

Usage:

  report_python_api.py --out <file.md>
"""

import argparse
import importlib
import inspect
import pkgutil
import re
import sys
from pathlib import Path

import lucifer

# The modules in reading order: what a user builds, runs and reads first, the tools last.
ORDER = ("models", "inputops", "runner", "stats", "beam", "plot", "flat", "schema", "generate")
# The generated models are listed by name: their fields are on the input reference.
MODELS = "models"


class Refusal(Exception):
    pass


def slug(*parts):
    return "pyapi-" + re.sub(r"[^a-z0-9]+", "-", "-".join(parts).lower()).strip("-")


def own_doc(obj):
    """The docstring the object itself carries, cleaned, and never one inherited from a
    base class, which inspect.getdoc would supply."""
    d = obj.__doc__ if not isinstance(obj, property) else obj.fget.__doc__
    return inspect.cleandoc(d) if d else ""


def doc_of(obj, where):
    d = own_doc(obj)
    if not d:
        raise Refusal(f"{where} has no docstring")
    return d


def signature(obj, bound=False):
    """The signature as Python prints it, annotations evaluated. With bound, the first
    parameter, self or cls, is left out, as Python's own documentation leaves it out."""
    try:
        sig = inspect.signature(obj, eval_str=True)
    except (TypeError, ValueError):
        return ""
    if bound:
        sig = sig.replace(parameters=list(sig.parameters.values())[1:])
    if inspect.isclass(obj):
        sig = sig.replace(return_annotation=inspect.Signature.empty)
    return str(sig)


def fenced(text):
    return ["```", text, "```", ""]


def members(cls, qual):
    """(name, kind, object) for the members a class defines itself, in definition order.
    A public name, and a special method that carries its own docstring."""
    out = []
    for name, raw in vars(cls).items():
        if name == "__init__":
            continue
        if isinstance(raw, property):
            kind, obj = "Property", raw
        elif isinstance(raw, classmethod):
            kind, obj = "Class method", raw.__func__
        elif isinstance(raw, staticmethod):
            kind, obj = "Static method", raw.__func__
        elif inspect.isfunction(raw):
            kind, obj = "Method", raw
        else:
            continue
        if name.startswith("_"):
            if not (name.startswith("__") and name.endswith("__") and own_doc(obj)):
                continue
        out.append((name, kind, obj))
    return out


def class_lines(mod, name, cls, short):
    qual = f"{mod}.{name}"
    L = [f"({slug(mod, name)})=", f"### `{name}`", ""]
    bases = [b for b in cls.__bases__ if b is not object]
    sig = "" if short else signature(cls)
    L.append(f"*Class* `{name}{sig}`")
    L.append("")
    if bases:
        L.append("*Bases* " + ", ".join(f"`{b.__module__}.{b.__qualname__}`"
                                         if b.__module__.startswith("lucifer") else f"`{b.__qualname__}`"
                                         for b in bases))
        L.append("")
    L += fenced(doc_of(cls, qual))
    if short:
        return L
    for mname, kind, obj in members(cls, qual):
        L.append(f"({slug(mod, name, mname)})=")
        L.append(f"#### `{name}.{mname}`")
        L.append("")
        L.append(f"*{kind}*" + ("" if kind == "Property" else f" `{mname}{signature(obj, bound=kind != 'Static method')}`"))
        L.append("")
        L += fenced(doc_of(obj, f"{qual}.{mname}"))
    return L


def module_lines(mname):
    mod = importlib.import_module(f"lucifer.{mname}")
    full = mod.__name__
    if not hasattr(mod, "__all__"):
        raise Refusal(f"{full} has no __all__")
    L = [f"({slug(full)})=", f"## `{full}`", "", "*Module*", ""]
    L += fenced(doc_of(mod, full))
    if mname == MODELS:
        L.append(f"The table classes' fields are the input keys of their tables, with the program's "
                 f"defaults, and [](../input-reference.md) lists each with its meaning, so they are "
                 f"not repeated here. `Input` takes its methods from `lucifer.inputops.InputOps`.")
        L.append("")
    n = 0
    for name in mod.__all__:
        if not hasattr(mod, name):
            raise Refusal(f"{full}.__all__ names {name}, which {full} does not define")
        obj = getattr(mod, name)
        if inspect.isclass(obj):
            L += class_lines(full, name, obj, short=(mname == MODELS))
        elif inspect.isfunction(obj):
            L += [f"({slug(full, name)})=", f"### `{name}`", "",
                  f"*Function* `{name}{signature(obj)}`", ""]
            L += fenced(doc_of(obj, f"{full}.{name}"))
        else:
            raise Refusal(f"{full}.__all__ names {name}, which is neither a class nor a function")
        n += 1
    return L, n


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip(),
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    found = {m.name for m in pkgutil.iter_modules(lucifer.__path__)}
    try:
        if found != set(ORDER):
            raise Refusal(f"the package's modules and ORDER differ: not in ORDER "
                          f"{sorted(found - set(ORDER))}, no such module {sorted(set(ORDER) - found)}")
        L = ["---", "title: Python reference", "short_title: Python reference", "---", "",
             "<!-- Generated by tests/scripts/report_python_api.py from the package's docstrings. "
             "Do not edit. -->", "",
             "This page lists every public module, class and function of the Python package "
             "`python/lucifer`, with its signature and its docstring as written. It is generated, "
             "so it says what the code says. A module's public names are its `__all__`. What the "
             "package is for, with examples, is [](../python.md).", ""]
        total = 0
        for m in ORDER:
            lines, n = module_lines(m)
            L += lines
            total += n
    except Refusal as err:
        print(f"report_python_api: REFUSED: {err}", file=sys.stderr)
        return 1
    text = "\n".join(L).rstrip("\n") + "\n"
    if re.search(r"0x[0-9a-f]{6,}", text):
        print("report_python_api: REFUSED: a signature carries a memory address", file=sys.stderr)
        return 1
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(text)
    print(f"wrote {args.out}: {total} entries from {len(ORDER)} modules")
    return 0


if __name__ == "__main__":
    sys.exit(main())
