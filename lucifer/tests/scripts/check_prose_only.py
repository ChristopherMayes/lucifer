#!/usr/bin/env python3
"""
Check that an edit to a Markdown page changed prose and nothing else, and count what the prose does.

  check_prose_only.py --base <rev> [--accept <file>] <page>...  # the page against <rev>
  check_prose_only.py --unwrapped <rev> <page>...                # text identical, lines aside
  check_prose_only.py --counts <page>...                         # prose counts, one row per page

The first form takes each page as it is at <rev> and as it is in the working tree, and
pulls out of each what a prose edit must leave alone: fenced blocks, directive openers
and their options, labels, headings, table cells, inline code, math, link targets, dates
and numbers. Each is a multiset, so a sentence can move within its page, and the two
sides must agree in every one. A difference is forgiven only by a line in the accept
file that names it:

  number - 64            the token 64 left the page
  number + 256           the token 256 arrived
  cell + Five permanent checks
  # a comment

so a changed number stands in the accept file beside the code that shows it, and a
changed table cell beside its reason. The check fails on any difference not named.

The second form requires the whitespace-normalized text of each page identical to the
text at <rev>, which is what a reflow of the lines must leave.

The third form counts, per page and in total: words, sentences, the median sentence
length, the 90th percentile, sentences over 20 and over 40 words, prose semicolons,
em-dashes, " -- ", dates, sentences naming the user or a goal, and hard-wrapped lines,
which are prose lines of 100 columns or less whose next line continues the paragraph.

Code blocks, indented blocks, tables, headings, labels and directive lines are not
prose and are not measured in the third form. The sentence counts take each line as one
paragraph or list item, so a page is reflowed to that form first, which the second
form proves.
"""

from __future__ import annotations

import argparse
import collections
import pathlib
import re
import statistics
import subprocess
import sys

WIDTH = 100

DATE = r"\b20\d\d-\d\d-\d\d\b"
NUMBER = r"(?<![\w.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![\w.]\w)"
STRUCTURAL = ("|", "#", ":::", "```", "(")


def at_rev(rev: str, path: str) -> str:
    top = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True,
                         text=True, check=True).stdout.strip()
    rel = str(pathlib.Path(path).resolve().relative_to(top))
    return subprocess.run(["git", "-C", top, "show", f"{rev}:{rel}"], capture_output=True,
                          text=True, check=True).stdout


def split_page(text: str):
    """(frozen, prose lines): what a prose edit must keep, and the lines it may change.

    frozen is a dict of category to multiset. prose is the list of lines outside fences,
    tables, headings, labels and directive lines, with the front matter removed.
    """
    frozen = {k: collections.Counter() for k in
              ("fence", "directive", "label", "heading", "cell", "code", "math", "link",
               "date", "number")}
    text = re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.S)
    prose = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        s = line.strip()
        if s.startswith("```"):
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith("```"):
                j += 1
            frozen["fence"]["\n".join(lines[i:j + 1])] += 1
            i = j + 1
            continue
        if s.startswith("$$"):
            if s.count("$$") == 2 and len(s) > 2:
                j = i
            else:
                j = i + 1
                while j < len(lines) and not (lines[j].strip().startswith("$$")
                                              or lines[j].strip().endswith("$$")):
                    j += 1
            frozen["math"][" ".join(" ".join(lines[i:j + 1]).split())] += 1
            i = j + 1
            continue
        if s.startswith(":::"):
            frozen["directive"][s] += 1
            # Options follow the opener as :key: value lines.
            j = i + 1
            while j < len(lines) and re.match(r"^\s*:[\w-]+:", lines[j]):
                frozen["directive"][lines[j].strip()] += 1
                j += 1
            i = j
            continue
        if re.match(r"^\([\w.-]+\)=\s*$", s):
            frozen["label"][s] += 1
        elif s.startswith("#"):
            frozen["heading"][s] += 1
        elif s.startswith("|"):
            if not re.match(r"^\|[\s:|-]+\|?\s*$", s):
                for cell in s.strip("|").split("|"):
                    frozen["cell"][" ".join(cell.split())] += 1
        else:
            prose.append(line)
        i += 1
    # A span may cross one line break on a page whose paragraphs are wrapped. Spans are
    # compared with their whitespace normalized, so a reflow does not change them.
    joined = "\n".join(prose)
    for code in re.findall(r"`[^`\n]*(?:\n[^`\n]*)?`", joined):
        frozen["code"][" ".join(code.split())] += 1
    rest = re.sub(r"`[^`\n]*(?:\n[^`\n]*)?`", " ", joined)
    for m in re.findall(r"\$[^$\n]+(?:\n[^$\n]+)?\$", rest):
        frozen["math"][" ".join(m.split())] += 1
    rest = re.sub(r"\$[^$\n]+(?:\n[^$\n]+)?\$", " ", rest)
    for target in re.findall(r"\]\(([^)\n]*(?:\n[^)\n]*)?)\)", rest):
        frozen["link"][" ".join(target.split())] += 1
    rest = re.sub(r"\]\([^)\n]*(?:\n[^)\n]*)?\)", "]", rest)
    for d in re.findall(DATE, rest):
        frozen["date"][d] += 1
    rest = re.sub(DATE, " ", rest)
    for n in re.findall(NUMBER, rest):
        frozen["number"][n] += 1
    return frozen, prose


def read_accept(path: str | None):
    accepted = collections.Counter()
    if path is None:
        return accepted
    for raw in pathlib.Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # The token may be empty, as a removed table row's empty cell is.
        m = re.match(r"^(\w+)\s+([-+])(?:\s+(.*))?$", line)
        if not m:
            raise SystemExit(f"{path}: cannot read accept line: {raw}")
        accepted[(m.group(1), m.group(2), " ".join((m.group(3) or "").split()))] += 1
    return accepted


def compare(rev: str, paths, accept: str | None) -> int:
    accepted = read_accept(accept)
    used = collections.Counter()
    bad = 0
    for path in paths:
        before, _ = split_page(at_rev(rev, path))
        after, _ = split_page(pathlib.Path(path).read_text())
        problems = []
        for cat in before:
            for sign, diff in (("-", before[cat] - after[cat]), ("+", after[cat] - before[cat])):
                for token, n in sorted(diff.items()):
                    key = (cat, sign, " ".join(token.split()))
                    forgiven = min(n, accepted[key] - used[key])
                    used[key] += forgiven
                    if n > forgiven:
                        problems.append(f"{cat} {sign} {token}" + (f"  (x{n - forgiven})" if n - forgiven > 1 else ""))
        if problems:
            bad += 1
            print(f"FAIL  {path}: {len(problems)} difference(s) from {rev} not in the accept file")
            for p in problems[:60]:
                print("      " + p.replace("\n", "\n        "))
            if len(problems) > 60:
                print(f"      ... and {len(problems) - 60} more")
        else:
            kept = sum(sum(c.values()) for c in after.values())
            print(f"ok    {path}: {kept} frozen items agree with {rev}")
    unused = {k: v - used[k] for k, v in accepted.items() if v > used[k]}
    for (cat, sign, token), n in sorted(unused.items()):
        print(f"note  accept line not needed: {cat} {sign} {token}")
    return bad


def unwrapped(rev: str, paths) -> int:
    bad = 0
    for path in paths:
        a = " ".join(at_rev(rev, path).split())
        b = " ".join(pathlib.Path(path).read_text().split())
        if a == b:
            print(f"ok    {path}: text identical to {rev}, lines aside")
        else:
            bad += 1
            k = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
            print(f"FAIL  {path}: text differs from {rev} near: ...{b[max(0, k - 40):k + 60]}...")
    return bad


ABBREVIATIONS = ("e.g", "i.e", "vs", "cf", "et", "al", "Fig", "Figs", "Eq", "Eqs", "Sec", "No",
                 "approx", "Dr", "St", "Phys", "Rev", "Accel", "Beams", "Nucl", "Instrum", "Meth",
                 "Methods", "Lett", "Opt", "Commun", "Proc", "Conf", "Vol", "pp", "ed", "eds", "Ref",
                 "Refs", "Appl", "Sci", "Rep", "Mod", "Part", "Univ", "Lab", "Inst", "Tech", "Trans",
                 "Soc", "Am", "Ann", "Int", "Eng", "Comput", "Res", "Suppl", "Ser", "J", "Chap",
                 "Instr", "Technol", "ST", "Stab")


def sentences(prose_lines):
    """The sentences of the prose, each as its text. A paragraph or list item is measured
    on its own, so a list never reads as one sentence, and a sentence may begin with a
    lowercase name such as openPMD."""
    out = []
    for line in prose_lines:
        # A span stands in as one word that no abbreviation or initial rule can match.
        t = re.sub(r"`[^`]*`", "span", line)
        t = re.sub(r"\$[^$]*\$", "span", t)
        t = re.sub(r"\[\]\([^)]*\)", "span", t)
        t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", t)
        t = re.sub(r"^\s*([-*]|\d+\.)\s+", "", t)
        t = t.replace("*", "")
        t = re.sub(r"(?<=\d)\.(?=\d)", "N", t)
        t = re.sub(r"\b(" + "|".join(re.escape(a) for a in ABBREVIATIONS) + r")\.", r"\1", t)
        t = re.sub(r"\b([A-Z])\.(?=\s)", r"\1", t)
        out += [s for s in re.split(r"(?<=[.!?])\s+", t) if len(s.split()) > 1]
    return out


def counts(paths) -> None:
    keys = ["words", "sents", "median", "p90", ">20", ">40", "semi", "emdash", "dbldash",
            "date", "user", "goal", "wrapped"]
    print(f"{'page':28s}" + "".join(f"{k:>8s}" for k in keys))
    tot = collections.Counter()
    all_lens = []
    for path in paths:
        _, prose = split_page(pathlib.Path(path).read_text())
        # An indented block is a command or a listing, and is neither measured nor wrapped.
        prose = [line for line in prose if not line.startswith("    ")]
        sents = sentences(prose)
        lens = [len(s.split()) for s in sents]
        all_lens += lens
        text = " ".join(prose)
        wrapped = sum(1 for a, b in zip(prose, prose[1:])
                      if a.strip() and b.strip() and len(a) <= WIDTH
                      and not b.lstrip().startswith(("-", "*", ":")) and not re.match(r"^\s*\d+\.", b))
        row = {"words": len(text.split()), "sents": len(sents),
               "median": statistics.median(lens) if lens else 0,
               "p90": sorted(lens)[int(0.9 * len(lens))] if lens else 0,
               ">20": sum(1 for n in lens if n > 20), ">40": sum(1 for n in lens if n > 40),
               "semi": len(re.findall(r"[a-z)]; \S", re.sub(r"`[^`]*`|\$[^$]*\$", "", text))),
               "emdash": text.count("—") + text.count("–"),
               "dbldash": text.count(" -- "),
               "date": len(re.findall(DATE, text)),
               "user": sum(1 for s in sents if re.search(r"\bthe user\b", s)),
               "goal": sum(1 for s in sents if re.search(r"\bgoal\b", s, re.I)),
               "wrapped": wrapped}
        for k in keys:
            if k not in ("median", "p90"):
                tot[k] += row[k]
        print(f"{pathlib.Path(path).name[:28]:28s}" + "".join(f"{row[k]:>8}" for k in keys))
    tot["median"] = statistics.median(all_lens) if all_lens else 0
    tot["p90"] = sorted(all_lens)[int(0.9 * len(all_lens))] if all_lens else 0
    print(f"{'total':28s}" + "".join(f"{tot[k]:>8}" for k in keys))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--base", help="Revision whose frozen items the pages must still carry")
    ap.add_argument("--accept", help="File of the differences from --base that are allowed")
    ap.add_argument("--unwrapped", metavar="REV", help="Require the text identical to REV, lines aside")
    ap.add_argument("--counts", action="store_true", help="Count the prose instead")
    ap.add_argument("pages", nargs="+")
    a = ap.parse_args()
    if a.counts:
        counts(a.pages)
        return 0
    if a.unwrapped:
        return 1 if unwrapped(a.unwrapped, a.pages) else 0
    if not a.base:
        ap.error("give --base <rev>, --unwrapped <rev> or --counts")
    return 1 if compare(a.base, a.pages, a.accept) else 0


if __name__ == "__main__":
    sys.exit(main())
