#!/usr/bin/env python3
"""
Check that a page links the glossary where it first uses one of the glossary's words.

  check_glossary.py [--root <repository root>] [--list]

The words are read from docs/glossary.md. An entry there is a label line, (gl-<name>)=,
followed by a paragraph that opens with the words in bold, **load, loader.** for
example. The words of one entry share its label. A word in code voice and the
qualifiers "first meaning", "second meaning", "field meaning" and "statistics meaning"
are not words to look for.

The pages are README.md and the Markdown files of docs/ that git tracks, less
docs/glossary.md itself, docs/changelog.md, whose entries are history, and the two
statistics specifications, which are documents of their own. The examples' READMEs
are not pages of the manual until report_examples.py writes them under
docs/generated/, and they are not read. Two things fail the check, each printed as
file:line: the text.

  A page uses a word of an entry and no link to that entry stands at or before the
  first use. The link is [text](glossary.md#gl-<name>), and the first use may be the
  link's own text. A word inside the text of a link to another entry, slice in
  [slice migration](glossary.md#gl-migration), is that entry's word and not a use.
  Text in a heading, in bold, in code, in
  mathematics, in a link's target and in a directive's options is not a use.

  A link names a label the glossary does not have.

The words in PLAIN are defined in the glossary and are not held to the first rule. Each
is also an ordinary word of English or of accelerator physics, record or step or beam,
and a page uses it in that sense before it uses it in the glossary's. A link at such a
use would point a reader at the wrong meaning. NO_PLURAL lists the words whose plural
is usually a verb, as in "the run banks the plane".

--list prints every word with its label and whether it is held to the first rule.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

GLOSSARY = "docs/glossary.md"
SKIP = (GLOSSARY, "docs/changelog.md", "docs/BMAD-STATS-SPEC.md", "docs/BMAD-STATS-EXT-FEL.md")
QUALIFIERS = ("first meaning", "second meaning", "field meaning", "statistics meaning")
PLAIN = {
    "wavelength", "head", "tail", "window order", "time order", "bunch", "beam", "weight",
    "grid", "plane", "polarization", "component", "harmonic", "member", "record",
    "slippage", "pulse", "line", "element", "undulator", "averaged", "unaveraged", "step",
    "dump", "frame", "slice range", "continuation", "header", "footer",
    "load", "shot noise", "kick", "push", "deposit", "source", "field solve", "kernel",
    "ramp", "sponge", "host", "FP64", "freerun", "check section", "level",
    "recorded digit", "estimate",
}
NO_PLURAL = {"bank", "walk", "seam", "pad", "twin", "window", "slice", "comb", "slab"}

LABEL = re.compile(r"^\(gl-([a-z0-9-]+)\)=\s*$")
HEAD = re.compile(r"^\*\*(.+?)\.\*\*")
GLOSSARY_LINK = re.compile(r"\[([^\]]*)\]\((?:docs/)?glossary\.md#gl-([a-z0-9-]+)\)")


def entries(root):
    """[(label, [words])] in the glossary's order."""
    out = []
    label = None
    for line in (root / GLOSSARY).read_text(encoding="utf-8").split("\n"):
        m = LABEL.match(line)
        if m:
            label = m.group(1)
            continue
        m = HEAD.match(line)
        if m:
            if label is None:
                raise SystemExit(f"check_glossary: no (gl-...)= label before the entry {m.group(0)}")
            words = [w.strip() for w in m.group(1).split(",")]
            out.append((label, [w for w in words if w not in QUALIFIERS and "`" not in w]))
        if line.strip():
            label = None
    return out


def word_pattern(word):
    """The word and its plural, not inside a longer word, a hyphenated compound or code."""
    body = re.escape(word)
    plural = "" if word in NO_PLURAL else "(?:e?s)?"
    flags = 0 if word.isupper() or word[:2].isupper() else re.I
    return re.compile(r"(?<![\w`/.#-])" + body + plural + r"(?![\w`/-])", flags)


def blank(text, pattern, group=0, flags=0):
    """The text with every match of the pattern, or one group of it, turned to spaces."""
    def sub(m):
        a, b = m.span(group)
        s = m.group(0)
        a -= m.start()
        b -= m.start()
        return s[:a] + re.sub(r"[^\n]", " ", s[a:b]) + s[b:]
    return re.sub(pattern, sub, text, flags=flags)


def prose(text):
    """The text with everything that is not prose turned to spaces, offsets kept.

    Links to the glossary are found in the original text, so a link's target is blanked
    here with every other target."""
    t = text
    t = blank(t, r"\A---\n.*?\n---\n", flags=re.S)                  # front matter
    t = blank(t, r"^(`{3,}|~{3,})[^\n]*\n.*?^\1[ \t]*$", flags=re.S | re.M)   # fenced blocks
    t = blank(t, r"<!--.*?-->", flags=re.S)
    t = blank(t, r"\$\$.*?\$\$", flags=re.S)
    t = blank(t, r"`[^`\n]*`")
    t = blank(t, r"\$[^$\n]*\$")
    t = blank(t, r"\*\*[^*\n]+\*\*")                                # bold: a list item's head
    t = blank(t, r"^#{1,6} [^\n]*$", flags=re.M)                    # headings
    t = blank(t, r"^\([A-Za-z0-9_-]+\)=[ \t]*$", flags=re.M)        # labels
    t = blank(t, r"^:{3,}[^\n]*$", flags=re.M)                      # directive lines
    t = blank(t, r"^:[a-z-]+:[^\n]*$", flags=re.M)                  # directive options
    t = blank(t, r"\]\((?:<[^>\n]*>|[^)\s]*)(?:\s+\"[^\"\n]*\")?\)")  # link targets
    return t


def tracked_pages(root):
    out = subprocess.run(["git", "ls-files", "--", "README.md", "docs/*.md"], cwd=root, check=True,
                         capture_output=True, text=True).stdout.split("\n")
    return [root / f for f in out if f and f.count("/") <= 1 and not f.startswith(SKIP)]


def first_uses(text, held):
    """{label: (offset, matched text)} for the first use of each held entry in the page."""
    t = prose(text)
    others = [(m.start(1), m.end(1), m.group(2)) for m in GLOSSARY_LINK.finditer(text)]
    found = {}
    for label, words in held:
        for w in words:
            for m in word_pattern(w).finditer(t):
                # "slice" inside the link text of "slice migration" is that entry's word.
                if any(a <= m.start() < b and other != label for a, b, other in others):
                    continue
                if label not in found or m.start() < found[label][0]:
                    found[label] = (m.start(), m.group(0))
                break
    return found


def check_page(text, held, labels):
    """([(offset, message)], {label: (offset, text)} of the uses with no link)."""
    bad = []
    linked = {}
    for m in GLOSSARY_LINK.finditer(text):
        if m.group(2) not in labels:
            bad.append((m.start(), f"no glossary entry gl-{m.group(2)}"))
        linked.setdefault(m.group(2), m.start())
    missing = {}
    for label, (at, word) in first_uses(text, held).items():
        if label in linked and linked[label] <= at:
            continue
        missing[label] = (at, word)
        bad.append((at, f"'{word}' is used before any link to glossary.md#gl-{label}"))
    return sorted(bad), missing


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2], type=Path)
    ap.add_argument("--list", action="store_true", help="print the words and stop")
    args = ap.parse_args()
    root = args.root.resolve()
    all_entries = entries(root)
    labels = {label for label, _ in all_entries}
    held = [(label, [w for w in words if w not in PLAIN]) for label, words in all_entries]
    held = [(label, words) for label, words in held if words]
    if args.list:
        for label, words in all_entries:
            for w in words:
                print(f"{'plain' if w in PLAIN else 'held ':5s}  gl-{label:18s} {w}")
        return 0
    n_bad = 0
    pages = tracked_pages(root)
    for path in pages:
        text = path.read_text(encoding="utf-8")
        for at, message in check_page(text, held, labels)[0]:
            line = text.count("\n", 0, at) + 1
            print(f"{path.relative_to(root)}:{line}: {message}")
            n_bad += 1
    n_words = sum(len(w) for _, w in held)
    if n_bad:
        print(f"check_glossary: {n_bad} failures in {len(pages)} pages. Link the word's first use on "
              "the page to its glossary entry, or reword the sentence.")
        return 1
    print(f"check_glossary: {len(pages)} pages, {n_words} words of {len(held)} entries held to a "
          f"link at their first use, {len(all_entries)} entries in the glossary")
    return 0


if __name__ == "__main__":
    sys.exit(main())
