#!/usr/bin/env python3
"""
Check that the code compares every name-valued input key against the names in its set and
no other, and that every name in the set is compared somewhere or is the key's default.

  check_allowed_values.py [--root <repository root>] [--list]

--list prints, for each key, every name a comparison holds and where the first one is.

A string key that takes one of a fixed set of names carries the set in the key table,
src/fel_schema_mod.f90, as the allowed argument of its add call, and fel_schema_set
refuses any other value as the input file is read. The routines that act on the value then
compare it against literal names. This check reads the sets from the table and the
component each key fills from the store in fel_schema_set. It follows each component
through the copies the code makes of it, an assignment to another variable or an argument
into a routine's dummy, and collects every comparison against a literal: an ==, a /= or a
case of a select on it. trim, upcase, downcase and adjustl around the value are seen
through. A component chain is followed by its last two parts, so a copy of a whole
structure, beam_init for run%beam_init, is followed too. A plain name is followed inside
the routine that made it a copy, a contained routine counting as part of its host.

Two failures, each printed with the file, the line and the key:

  a compared name outside the set, a value the program would act on that the input file
  can never hold. A set that matches in any case is compared upcased, so its names are
  compared upper case.
  a name in the set that no comparison holds and that is not the key's default (read from
  python/lucifer/schema.json), a value the input file can hold that the program ignores.

It needs no build, so it runs in the docs workflow and in the examples job.
"""
import argparse
import json
import re
import sys
from pathlib import Path

WRAPPER = re.compile(r"^(trim|upcase|downcase|adjustl)\((.*)\)$")
LITERAL = r"'((?:[^']|'')*)'"
CHAIN = r"[A-Za-z_]\w*(?:\(\d+\))?(?:%\w+(?:\(\d+\))?)*"
OPERAND = r"((?:trim|upcase|downcase|adjustl)\(\s*" + CHAIN + r"\s*\)|" + CHAIN + ")"
COMPARE_L = re.compile(OPERAND + r"\s*(?:==|/=)\s*" + LITERAL, re.I)
COMPARE_R = re.compile(LITERAL + r"\s*(?:==|/=)\s*" + OPERAND, re.I)
ASSIGN = re.compile(r"^\s*(" + CHAIN + r")\s*=(?!=)\s*(.+)$")
ROUTINE = re.compile(r"^\s*(?:(?:real|integer|logical|complex|character|type)\s*(?:\([^)]*\))?\s+)?"
                     r"(?:(?:recursive|pure|elemental|impure)\s+)*(subroutine|function)\s+(\w+)\s*(?:\(([^)]*)\))?",
                     re.I)
END_ROUTINE = re.compile(r"^\s*end\s+(subroutine|function)\b", re.I)
CALLED = re.compile(r"\b(\w+)\s*\(")


def norm(expr):
    """An operand without blanks, lower case, with the wrappers taken off."""
    e = re.sub(r"\s+", "", expr).lower()
    m = WRAPPER.match(e)
    while m:
        e = m.group(2)
        m = WRAPPER.match(e)
    return e


def tail(e):
    """A chain's last two parts, which a copy of the structure that holds it still has."""
    parts = e.split("%")
    return "%".join(parts[-2:]) if len(parts) > 1 else e


def statements(path):
    """(first line number, text) of each statement of a Fortran file, with continuations
    joined, comments removed and semicolon-separated statements apart."""
    out, buf, first = [], "", None
    for n, raw in enumerate(path.read_text(errors="replace").split("\n"), 1):
        kept, inq = [], False
        for ch in raw:
            if ch == "'":
                inq = not inq
            if ch == "!" and not inq:
                break
            kept.append(ch)
        text = "".join(kept).rstrip()
        if not text.strip():
            continue
        if first is None:
            first = n
        lead = text.lstrip()
        if buf and lead.startswith("&"):
            lead = lead[1:]
        buf = buf + " " + lead if buf else text
        if buf.endswith("&"):
            buf = buf[:-1]
            continue
        cur, inq = "", False
        for ch in buf:
            if ch == "'":
                inq = not inq
            if ch == ";" and not inq:
                if cur.strip():
                    out.append((first, cur))
                cur = ""
            else:
                cur += ch
        if cur.strip():
            out.append((first, cur))
        buf, first = "", None
    return out


def split_args(text):
    """The top-level comma-separated arguments of an argument list."""
    args, depth, cur, inq = [], 0, "", False
    for ch in text:
        if ch == "'":
            inq = not inq
        if not inq and ch in "()":
            depth += 1 if ch == "(" else -1
        if not inq and ch == "," and depth == 0:
            args.append(cur.strip())
            cur = ""
            continue
        cur += ch
    if cur.strip():
        args.append(cur.strip())
    return args


def calls(stmt, known):
    """(routine, arguments) for each call in a statement of a routine in known."""
    out = []
    for m in CALLED.finditer(stmt):
        name = m.group(1).lower()
        if name not in known:
            continue
        depth, i = 1, m.end()
        while i < len(stmt) and depth:
            depth += {"(": 1, ")": -1}.get(stmt[i], 0)
            i += 1
        out.append((name, split_args(stmt[m.end():i - 1])))
    return out


def table_sets(schema_src):
    """{table.name: (names, any_case)} from the add calls of the key table."""
    sets = {}
    for n, s in statements(schema_src):
        m = re.match(r"^\s*call add \('([\w.]*)', '(\w+)',(.*)\)\s*$", s)
        if not m:
            continue
        a = re.search(r"allowed = \[character\(\d+\) ::(.*?)\]", m.group(3))
        if a:
            names = [x.replace("''", "'") for x in re.findall(LITERAL, a.group(1))]
            sets[f"{m.group(1)}.{m.group(2)}"] = (names, "any_case = .true." in m.group(3))
    return sets


def store_components(schema_src):
    """{table.name: component} from the store of fel_schema_set."""
    comp = {}
    for m in re.finditer(r"case \('([\w.]+)'\);\s*call put \(([\w%()]+)\)", schema_src.read_text()):
        comp[m.group(1)] = norm(m.group(2))
    return comp


def schema_defaults(schema_json):
    """{table.name: default} from the JSON Schema, the lattice as .lattice."""
    out = {}

    def walk(props, prefix):
        for k, v in props.items():
            if v.get("type") == "object":
                walk(v["properties"], f"{prefix}{k}.")
            else:
                out[f"{prefix}{k}" if prefix else f".{k}"] = v.get("default")
    walk(json.loads(schema_json.read_text())["properties"], "")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]),
                    help="the repository root (default: this script's own tree)")
    ap.add_argument("--list", action="store_true", help="print each key's compared names and where")
    opts = ap.parse_args()
    root = Path(opts.root).resolve()
    schema_src = root / "src" / "fel_schema_mod.f90"
    sets = table_sets(schema_src)
    comps = store_components(schema_src)
    defaults = schema_defaults(root / "python" / "lucifer" / "schema.json")

    # The scopes: each top-level routine of each source file, with the statements of its
    # contained routines. routines maps a name to its dummies and the scope that holds it.
    files = sorted((root / "src").glob("*.f90")) + sorted((root / "programs").glob("*.f90"))
    scopes, routines = [], {}
    for f in files:
        stmts = statements(f)
        stack = []
        for idx, (n, s) in enumerate(stmts):
            if END_ROUTINE.match(s):
                name, start, dummies = stack.pop()
                if not stack:
                    scopes.append((f, stmts[start:idx + 1]))
                    sid = len(scopes) - 1
                    for nm, dm in pending:
                        routines[nm] = (dm, sid)
                continue
            m = ROUTINE.match(s)
            if m:
                if not stack:
                    pending = []
                dummies = [d.strip().lower() for d in (m.group(3) or "").split(",") if d.strip()]
                stack.append((m.group(2).lower(), idx, dummies))
                pending.append((m.group(2).lower(), dummies))
    calls_in = [[calls(s, routines) for n, s in stmts] for f, stmts in scopes]

    # Follow each key to a fixed point. held holds ("*", chain tail) and (scope, name).
    def holds(held, sid, expr):
        e = norm(expr)
        return ("*", tail(e)) in held if "%" in e else (sid, e) in held

    tracked = {key: {("*", tail(c))} for key, c in comps.items() if key in sets}
    for key, held in tracked.items():
        changed = True
        while changed:
            changed = False
            for sid, (f, stmts) in enumerate(scopes):
                for (n, s), cl in zip(stmts, calls_in[sid]):
                    new = []
                    m = ASSIGN.match(s)
                    if m and holds(held, sid, m.group(2)):
                        lhs = norm(m.group(1))
                        new.append(("*", tail(lhs)) if "%" in lhs else (sid, lhs))
                    for name, args in cl:
                        dummies, rsid = routines[name]
                        for i, a in enumerate(args):
                            kw = re.match(r"^(\w+)\s*=(?!=)\s*(.+)$", a)
                            dummy, val = (kw.group(1).lower(), kw.group(2)) if kw else (
                                dummies[i] if i < len(dummies) else None, a)
                            if dummy and holds(held, sid, val):
                                new.append((rsid, dummy))
                    for x in new:
                        if x not in held:
                            held.add(x)
                            changed = True

    # The comparisons.
    failures, compared = [], {key: {} for key in tracked}
    for key, held in tracked.items():
        names, any_case = sets[key]
        for sid, (f, stmts) in enumerate(scopes):
            rel = f.relative_to(root)
            on_select = []
            for n, s in stmts:
                low = s.strip().lower()
                m = re.match(r"^select\s+case\s*\((.*)\)\s*$", low)
                if m:
                    on_select.append(holds(held, sid, m.group(1)))
                    continue
                if re.match(r"^end\s*select\b", low):
                    if on_select:
                        on_select.pop()
                    continue
                lits = []
                m = re.match(r"^case\s*\((.*)\)\s*$", s.strip(), re.I)
                if m and on_select and on_select[-1]:
                    lits += [x.replace("''", "'") for x in re.findall(LITERAL, m.group(1))]
                for c in COMPARE_L.finditer(s):
                    if holds(held, sid, c.group(1)):
                        lits.append(c.group(2).replace("''", "'"))
                for c in COMPARE_R.finditer(s):
                    if holds(held, sid, c.group(2)):
                        lits.append(c.group(1).replace("''", "'"))
                for lit in lits:
                    compared[key].setdefault(lit, f"{rel}:{n}")
                    if lit not in names:
                        why = ", compared upper case since the set matches in any case" if any_case else ""
                        failures.append(f"{rel}:{n}: {key}: \"{lit}\" is compared and is not in the set "
                                        + ", ".join(f'"{x}"' for x in names) + why)

    for key, (names, any_case) in sorted(sets.items()):
        if key not in tracked:
            failures.append(f"src/fel_schema_mod.f90: {key}: has a set and the store fills no component")
            continue
        dflt = defaults.get(key)
        for name in names:
            is_default = dflt is not None and (name.upper() == str(dflt).upper() if any_case else name == dflt)
            if name not in compared[key] and not is_default:
                failures.append(f"src/fel_schema_mod.f90: {key}: \"{name}\" is in the set, no comparison "
                                f"holds it and it is not the default")

    if opts.list:
        for key in sorted(compared):
            print(f"{key}: " + ", ".join(f'"{x}" {w}' for x, w in compared[key].items()))
    for line in failures:
        print(line)
    print(f"check_allowed_values: {len(sets)} keys with a set, "
          f"{sum(len(v) for v in compared.values())} compared names, {len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
