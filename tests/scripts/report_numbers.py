#!/usr/bin/env python3
"""
Classify every number on the documentation pages, and say which ones are typed by hand.

  report_numbers.py [--unmatched <tsv>] [--list <class>]

Every numeric token in the prose, the table cells and the code spans of the pages is
read and classified by matching its value against the sources a number can come from:

  default    a default of an input structure, read from the Fortran struct files
  standard   a constant of a derivation or a format version: the parameter constants of
             fel_struct.f90 and the committed list tests/number_standards.tsv
  tolerance  a tolerance or bound a check script holds a value to: its module-level
             constants, the literal bound in a check() call, and the bounds of the
             generated tables
  measured   a value a generated page carries, from the validation suite's records, or a value of
             the examples' record, examples/measured.tsv
  example    a setting of an input file under examples/ or tests/

A token that matches no source is classified by hand in tests/number_classes.tsv, which
this tool reads first, so a hand classification stands and the rest is rerunnable. The
hand table may also say citation, for a year, a volume or an identifier in a reference,
and count, for a small integer that counts things in a sentence.

Three rules classify what no source matches, and the hand table overrides them. A
small integer, 0 to 12, with no decimal point or exponent is a count, since such a
number counts things far more often than it sets them. A year, and every number on
references.md, is a citation. On a record page, validation.md, performance.md,
startup-noise.md and the example READMEs, a number that matches no source is a
measurement typed from a run. On the manual pages such a number is unmatched and goes
to the hand table.

The report prints counts per page and class, then the measured numbers that no generated
page carries, which are the numbers still typed by hand and the work of the goals that
generate them. --unmatched writes the unmatched tokens as rows for the hand table, and
--list prints every token of one class with its page and sentence.
"""

from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from check_prose_only import DATE, NUMBER  # noqa: E402

LUCIFER = pathlib.Path(__file__).resolve().parents[2]
CLASSES = ("default", "standard", "tolerance", "measured", "example")
HAND_ONLY = ("citation", "count", "other")
SMALL = 12
RECORD_PAGES = ("validation.md", "performance.md", "startup-noise.md")


def tokens_of(path: pathlib.Path):
    """(token, where, line number, context) for every numeric token of a page, with
    where one of prose, cell or code. Math, fences, labels and headings are skipped."""
    text = re.sub(r"\A---\n.*?\n---\n", "", path.read_text(), flags=re.S)
    out = []
    fence = math = False
    for k, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if s.startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        if s.startswith("$$"):
            math = not (s.count("$$") == 2 and len(s) > 2) if not math else False
            continue
        if math or s.startswith((":::", "#", "<!--")) or re.match(r"^\([\w.-]+\)=\s*$", s):
            continue
        if re.match(r"^\s*:[\w-]+:", s):
            continue
        is_row = s.startswith("|")
        if is_row and re.match(r"^\|[\s:|-]+\|?\s*$", s):
            continue
        for code in re.findall(r"`[^`]+`", line):
            for n in re.findall(NUMBER, code):
                out.append((n, "code", k, code))
        rest = re.sub(r"`[^`]+`", " ", line)
        rest = re.sub(r"\$[^$]+\$", " ", rest)
        rest = re.sub(r"\]\([^)]*\)", "]", rest)
        rest = re.sub(DATE, " ", rest)
        if is_row:
            for cell in rest.strip("|").split("|"):
                for n in re.findall(NUMBER, cell):
                    out.append((n, "cell", k, " ".join(cell.split())))
        else:
            for m in re.finditer(NUMBER, rest):
                a, b = max(0, m.start() - 60), min(len(rest), m.end() + 60)
                out.append((m.group(0), "prose", k, " ".join(rest[a:b].split())))
    return out


def as_float(token: str):
    try:
        return float(token)
    except ValueError:
        return None


def small_integer(token: str) -> bool:
    return re.fullmatch(r"[-+]?\d+", token) is not None and abs(int(token)) <= SMALL


class Sources:
    """The values each class can match, as {value: [names]}."""

    def __init__(self):
        self.values = {c: collections.defaultdict(list) for c in CLASSES}

    def add(self, cls, token, name):
        v = as_float(str(token).replace("_rp", "").replace("d", "e").replace("D", "e"))
        if v is not None:
            self.values[cls][v].append(name)

    def read_struct(self, path: pathlib.Path):
        """Defaults, and the parameter constants as standards, of a Fortran struct file."""
        for line in path.read_text().splitlines():
            m = re.match(r"^\s*(?:real\(rp\)|integer|logical)\s*::\s*(\w+)\s*=\s*([-+0-9.eEdD_rp]+)", line)
            if m:
                self.add("default", m.group(2), f"{path.name}:{m.group(1)}")
            m = re.match(r"^\s*(?:real\(rp\)|integer)\s*,\s*parameter\s*::\s*(\w+\$?)\s*=\s*([-+0-9.eEdD_rp]+)", line)
            if m:
                self.add("standard", m.group(2), f"{path.name}:{m.group(1)}")

    def read_tsv(self, path: pathlib.Path, cls: str):
        if not path.exists():
            return
        for row in csv.DictReader(path.open(), delimiter="\t"):
            if row.get("value"):
                self.add(cls, row["value"], row.get("name", ""))

    def read_scripts(self, paths):
        """Module-level numeric constants and the literal bound of each check() call."""
        for p in paths:
            text = p.read_text()
            for m in re.finditer(r"^([A-Z][A-Z0-9_]+)\s*=\s*([-+0-9.eE]+)\s*(?:#|$)", text, re.M):
                self.add("tolerance", m.group(2), f"{p.name}:{m.group(1)}")
            for m in re.finditer(r"\bcheck\((?:[^()]|\([^()]*\))*?,\s*([-+0-9.eE]+)\s*[,)]", text):
                self.add("tolerance", m.group(1), f"{p.name}:check")

    def read_generated(self, path: pathlib.Path):
        """Measured values from a generated page's table cells, and bounds as tolerances."""
        if not path.exists():
            return
        for line in path.read_text().splitlines():
            if not line.startswith("|") or re.match(r"^\|[\s:|-]+\|?\s*$", line):
                continue
            cells = [c.strip() for c in line.strip("|").split("|")]
            for i, cell in enumerate(cells):
                cls = "tolerance" if cells[0] == "Check" or (len(cells) == 4 and i == 3 and ("<=" in cell or ">=" in cell or " to " in cell or cell.startswith("="))) else "measured"
                for n in re.findall(NUMBER, cell):
                    self.add(cls, n, f"{path.name}:{cells[0][:30]}")

    def read_record(self, path: pathlib.Path):
        """Measured values from the examples' record, one row per input file and quantity."""
        if not path.exists():
            return
        for row in csv.DictReader(path.open(), delimiter="\t"):
            if row.get("value") not in (None, "", "yes", "no"):
                self.add("measured", row["value"], f"measured.tsv:{row['directory']}/{row['input']}:{row['quantity']}")

    def read_inputs(self, roots):
        """Settings of the input files, and of the namelists the check scripts carry as text."""
        for root in roots:
            files = list(root.rglob("*.in")) + list(root.rglob("*.bmad")) + list(root.rglob("*.lat"))
            files += list(root.rglob("check_*.py"))
            for p in files:
                for m in re.finditer(r"^\s*([\w%()\[\]:*]+)\s*=\s*([-+0-9.eEdD]+)\s*(?:[,!;]|$)", p.read_text(), re.M):
                    self.add("example", m.group(2), f"{p.parent.name}/{p.name}:{m.group(1)}")

    def classify(self, token: str):
        if small_integer(token):
            return None, None
        v = as_float(token)
        if v is None:
            return None, None
        digits = significant_digits(token)
        for cls in CLASSES:
            names = self.values[cls].get(v)
            if names:
                return cls, names[0]
            # A page quotes a value to fewer digits than its source carries, so the source
            # rounded to the token's digits must equal the token. Two digits or more, since
            # one digit matches too much.
            if digits >= 2:
                for src, nm in self.values[cls].items():
                    if src != 0 and float(f"{src:.{digits - 1}e}") == float(f"{v:.{digits - 1}e}"):
                        return cls, nm[0]
        return None, None


def significant_digits(token: str) -> int:
    mant = re.split(r"[eE]", token)[0].lstrip("+-").replace(".", "").lstrip("0")
    return len(mant)


def rule_class(page: str, token: str):
    """The class a rule gives a token no source matched, or None for the hand table."""
    if small_integer(token):
        return "count", "rule: a small integer"
    if page == "references.md" or re.fullmatch(r"(19|20)\d\d", token):
        return "citation", "rule: a year or a reference entry"
    if page in RECORD_PAGES or page.endswith("/README.md"):
        return "measured", "rule: a number on a record page that no source gives"
    return None, None


def read_hand(path: pathlib.Path):
    """{(page, token, context key): (class, note)} from the hand table."""
    hand = {}
    if not path.exists():
        return hand
    for row in csv.DictReader(path.open(), delimiter="\t"):
        hand[(row["page"], row["token"], row["context"])] = (row["class"], row.get("note", ""))
    return hand


def hand_lookup(hand, page, token, context):
    for (hp, ht, hc), v in hand.items():
        if hp == page and ht == token and hc in context:
            return v
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--doc", default=str(LUCIFER / "docs"))
    ap.add_argument("--examples", default=str(LUCIFER / "examples"))
    ap.add_argument("--tests", default=str(LUCIFER / "tests"))
    ap.add_argument("--code", default=str(LUCIFER / "src"))
    ap.add_argument("--classes", default=str(LUCIFER / "tests" / "number_classes.tsv"))
    ap.add_argument("--standards", default=str(LUCIFER / "tests" / "number_standards.tsv"))
    ap.add_argument("--unmatched", help="Write the unmatched tokens here as hand-table rows")
    ap.add_argument("--list", help="Print every token of this class")
    a = ap.parse_args()

    doc, examples, tests, code = (pathlib.Path(x) for x in (a.doc, a.examples, a.tests, a.code))
    pages = sorted(p for p in doc.glob("*.md") if p.name != "changelog.md")
    pages += sorted(examples.glob("*/README.md")) + [examples / "README.md"]

    src = Sources()
    for f in ("fel_struct.f90", "fel_import_mod.f90", "fel_collective_mod.f90"):
        if (code / f).exists():
            src.read_struct(code / f)
    src.read_tsv(pathlib.Path(a.standards), "standard")
    src.read_scripts(sorted((tests / "scripts").glob("check_*.py")))
    src.read_generated(doc / "generated" / "validation-measured.md")
    if (doc / "generated" / "startup-noise-measured.md").exists():
        src.read_generated(doc / "generated" / "startup-noise-measured.md")
    src.read_record(examples / "measured.tsv")
    src.read_inputs([examples, tests / "bmad", tests / "genesis4"])
    hand = read_hand(pathlib.Path(a.classes))

    rows = []
    for page in pages:
        pname = page.parent.name + "/README.md" if page.name == "README.md" else page.name
        for token, where, k, context in tokens_of(page):
            h = hand_lookup(hand, pname, token, context)
            if h:
                cls, name = h[0], "hand: " + h[1]
            else:
                cls, name = src.classify(token)
                if cls is None:
                    cls, name = rule_class(pname, token)
            rows.append((pname, token, where, k, context, cls or "unmatched", name or ""))

    classes = list(CLASSES) + list(HAND_ONLY) + ["unmatched"]
    print(f"{'page':28s}" + "".join(f"{c[:9]:>10s}" for c in classes) + f"{'total':>8s}")
    tot = collections.Counter()
    for pname in dict.fromkeys(r[0] for r in rows):
        c = collections.Counter(r[5] for r in rows if r[0] == pname)
        tot.update(c)
        print(f"{pname[:28]:28s}" + "".join(f"{c[k]:10d}" for k in classes) + f"{sum(c.values()):8d}")
    print(f"{'total':28s}" + "".join(f"{tot[k]:10d}" for k in classes) + f"{sum(tot.values()):8d}")

    by_hand = [r for r in rows if r[5] == "measured" and not r[6].startswith("validation-measured")]
    print(f"\nmeasured numbers no generated page carries: {len(by_hand)}")
    per_page = collections.Counter(r[0] for r in by_hand)
    for pname, n in per_page.most_common():
        print(f"  {pname:28s} {n}")

    if a.list:
        print()
        for r in rows:
            if r[5] == a.list:
                print(f"{r[0]}:{r[3]}\t{r[1]}\t{r[2]}\t{r[6]}\t{r[4][:100]}")
    if a.unmatched:
        with open(a.unmatched, "w") as fh:
            fh.write("page\ttoken\tcontext\tclass\tnote\n")
            for r in rows:
                if r[5] == "unmatched":
                    fh.write(f"{r[0]}\t{r[1]}\t{r[4][:80]}\t\t\n")
        print(f"\nwrote {sum(1 for r in rows if r[5] == 'unmatched')} unmatched rows to {a.unmatched}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
