#!/usr/bin/env python3
"""
Generate docs/generated/examples/ from the example directories themselves.

Each page carries the example's own README prose, the figure the runner saved, the
case of the example (its lattice and its beam, a table before the prose), the
quantities the runner recorded for each of its input files (examples/measured.tsv, one table
per example with a column per input file), and every input file it runs on, fenced verbatim
from the file on disk. The table of the case is written from the same record, in the
units a reader uses, with one column where the input files of an example do not differ. The embedded text is read
rather than transcribed, so an input file shown in the documentation cannot drift from the input file
that runs: editing the input file and regenerating is the only way to change the page.

The output is deterministic. Nothing time dependent, machine dependent or path dependent
is written, because the keystone suite regenerates these files and requires the diff to
be empty. Figures are referenced, never read. A figure's bytes vary with the plotting
library, so run_examples.sh writes them and check_examples.py asserts only that each one
exists.

Usage:

  report_examples.py --examples <dir> --out <dir>

Every relative link in a README is rewritten for the page's new depth, and any link form
this script does not recognize is an error rather than a guess. A wrong rewrite would
leave a dead link on the page with no message.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

# Input files shown on a page, in the order they appear. A run's own input file first,
# then the reference code's input files, then the lattices.
INPUT_GLOBS = ("lucifer.toml", "lucifer[-_]*.toml", "sat-*.toml", "sat-*.in")
LATTICE_GLOBS = ("*.bmad", "*.lat")

# An embedded file is meant to be read on a page. Anything larger is data, not input,
# and belongs in the repository only (examples/bmad_wake/ztable.wake is the case).
MAX_EMBED_BYTES = 4000

# run_examples.sh does not run these, so their pages carry no figure. saturation_demo
# needs the genesis4 binary and its own clock, and it has run.sh and a report of its own.
# lcls_electrons tracks every electron of its window and takes over an hour.
NO_FIGURE = ("saturation_demo", "lcls_electrons")

# The index page's own name. The site serves a page under its file stem, so this cannot
# be index.md: that collides with the site's own index and is served as /index-1.
INDEX_PAGE = "examples.md"

LANG = {".toml": "toml", ".in": "fortran", ".bmad": "fortran", ".lat": "text"}

LINK = re.compile(r"\]\(([^)]+)\)")


def rewrite_links(text, names, where):
    """Retarget a README's relative links for a page in docs/generated/examples/."""
    def sub(m):
        t = m.group(1)
        if t.startswith(("http://", "https://", "#")):
            return m.group(0)
        # A doc page, from a per-directory README (../../docs/x.md) or the index (../docs/x.md),
        # with an optional section anchor on the end (../../docs/x.md#sec-y).
        d = re.fullmatch(r"\.\./(?:\.\./)?docs/([A-Za-z0-9._-]+\.md)(#[A-Za-z0-9._-]+)?", t)
        if d:
            return f"](../../{d.group(1)}{d.group(2) or ''})"
        # A section of the index of the examples, from a per-directory README. The index
        # page carries a label for each of its sections, which label_sections() writes.
        i = re.fullmatch(r"\.\./README\.md#([a-z0-9-]+)", t)
        if i:
            return f"]({INDEX_PAGE}#ex-{i.group(1)})"
        # A sibling example, which the index names as a directory.
        s = re.fullmatch(r"([A-Za-z0-9_]+)/", t)
        if s and s.group(1) in names:
            return f"]({s.group(1)}.md)"
        raise SystemExit(f"{where}: unrecognized link target {t!r}. Teach "
                         "report_examples.py this form rather than letting it guess.")
    return LINK.sub(sub, text)


def split_title(text, where):
    """Return (title line, body) for a README whose first heading is its title."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("# "):
            return line[2:].strip(), "\n".join(lines[i + 1:]).lstrip("\n")
    raise SystemExit(f"{where}: no '# ' title line")


def gather(d):
    """The input files of one example, deduplicated, in a fixed order."""
    out = []
    for globs in (INPUT_GLOBS, LATTICE_GLOBS):
        found = []
        for g in globs:
            found += [p for p in d.glob(g) if p.is_file()]
        for p in sorted(set(found), key=lambda p: (p.name != "lucifer.toml", p.name)):
            if p not in out:
                out.append(p)
    return out


def front_matter(title, short_title, source, tags):
    """Frontmatter, the generated-by note, and nothing else.

    The tags are the kind of the page, which check_page_kinds.py holds every page of the
    site's contents to. The page of an example is the guide of its README followed by what
    its run measured, and the page that lists the examples is a reference.

    The title lives here rather than as a body heading. The site derives a page's title
    from its first heading when there is no frontmatter, and it removes that heading only
    when the heading leads the file. Every page here opens with the note below, so a body
    heading would be lifted into the title and left in place, and the reader would see it
    twice. YAML also needs the quoting: most of these titles carry a colon."""
    if '"' in title or "\\" in title:
        raise SystemExit(f"{source}: title needs escaping this script does not do: {title!r}")
    return ["---", f'title: "{title}"', f'short_title: "{short_title}"',
            f"tags: [{tags}]", "---", "",
            f"<!-- Generated by {Path(__file__).name} from {source}. Do not edit. -->", ""]


def label_sections(text):
    """The index page's text with a label before each of its sections, (ex-<slug>)=, where
    the slug is the heading as a forge writes its anchor: lower case, spaces as hyphens.
    A README links a section as ../README.md#<slug>, which rewrite_links() retargets."""
    out = []
    for line in text.split("\n"):
        if line.startswith("## "):
            slug = re.sub(r"[^a-z0-9 -]", "", line[3:].strip().lower()).replace(" ", "-")
            out += [f"(ex-{slug})=", line]
        else:
            out.append(line)
    return "\n".join(out)


def embed(page, path, label=None):
    size = path.stat().st_size
    if size > MAX_EMBED_BYTES:
        raise SystemExit(f"{path}: {size} bytes is too large to embed on a page. Add it "
                         "to the exclusions or shorten it.")
    page.append(f"### `{label or path.name}`")
    page.append("")
    page.append("```" + LANG.get(path.suffix, "text"))
    page.append(path.read_text().rstrip("\n"))
    page.append("```")
    page.append("")


ORDER = ["exit_power", "exit_power_per_slice", "exit_power_slice_spread", "exit_window_energy",
         "exit_pulse_energy",
         "peak_power", "peak_power_s", "mode_peak_inside",
         "mode_peak_s", "exit_power_inside", "exit_power_outside", "exit_power_inside_s", "wide_angle_ratio",
         "wide_angle_ratio_exit", "split_angle", "gain_length", "power_after_first_element",
         "power_after_first_element_per_slice",
         "initial_bunching", "exit_bunching", "exit_bunching_max", "bunching_max", "bunching_max_s",
         "bunching_peak", "slice_power_max_at_bunching_max", "exit_power_slice_max", "initial_energy",
         "exit_energy", "energy_change", "initial_sigma_energy", "exit_sigma_energy",
         "exit_energy_span", "n_slice",
         "n_record", "slice_spacing", "n_particle", "grid_n_pts", "grid_cell", "grid_derived",
         "grid_derived_from", "fel_steps", "filter_mode_angle", "filter_rho_angle", "filter_ratio",
         "filter_edge"]
WORDS = {"exit_power": "exit power", "exit_power_per_slice": "exit power per slice",
         "exit_power_slice_spread": "exit power, rms over mean across slices",
         "exit_power_inside_s": "s of the last record with the split",
         "power_after_first_element_per_slice": "power per slice at the first element end",
         "exit_window_energy": "exit window energy",
         "exit_pulse_energy": "exit pulse energy, the window's and the light slipped out of it",
         "peak_power": "peak power", "peak_power_s": "s of the peak power",
         "mode_peak_inside": "mode peak, power inside the split", "mode_peak_s": "s of the mode peak",
         "exit_power_inside": "exit power inside the split", "exit_power_outside": "exit power outside the split",
         "wide_angle_ratio": "outside over inside, the run's maximum",
         "wide_angle_ratio_exit": "outside over inside at the exit", "split_angle": "split angle",
         "gain_length": "gain length, the shortest e-folding length over five consecutive records",
         "power_after_first_element": "power at the first element end",
         "exit_bunching": "exit bunching, mean |b| over slices", "exit_bunching_max": "exit bunching, the largest slice",
         "bunching_peak": "bunching at the mode peak", "initial_bunching": "initial bunching, mean |b| over slices",
         "bunching_max": "bunching, the maximum along the line of the mean |b|",
         "bunching_max_s": "s of the bunching maximum", "exit_power_slice_max": "exit power of the largest slice",
         "slice_power_max_at_bunching_max": "power of the largest slice at the bunching maximum",
         "initial_energy": "mean slice energy at the start", "exit_energy": "mean slice energy at the exit",
         "energy_change": "mean slice energy, exit minus start",
         "exit_energy_span": "span of the slice energies at the exit, head to tail",
         "initial_sigma_energy": "slice energy spread at the start", "exit_sigma_energy": "slice energy spread at the exit",
         "n_slice": "slices", "n_record": "records", "slice_spacing": "slice spacing",
         "n_particle": "macroparticles per slice", "grid_n_pts": "grid points", "grid_cell": "grid cell",
         "grid_derived": "grid derived by setup", "grid_derived_from": "grid rounded up from",
         "n_particle_derived": "macroparticles per slice derived by setup", "beamlet_size": "beamlet size",
         "beamlet_size_derived": "beamlet size derived by setup",
         "fel_steps": "FEL steps", "filter_mode_angle": "source filter, mode angle",
         "filter_rho_angle": "source filter, rho angle", "filter_ratio": "source filter, ratio of the angles",
         "filter_edge": "source filter, edge"}


def read_record(path):
    """{directory: {input: {quantity: (value, unit)}}} from examples/measured.tsv."""
    rec = {}
    if not path.exists():
        return rec
    for row in csv.DictReader(path.open(), delimiter="\t"):
        rec.setdefault(row["directory"], {}).setdefault(row["input"], {})[row["quantity"]] = (row["value"], row["unit"])
    return rec


# The quantities of the record that the table of the case shows and the table of measured
# quantities leaves out. initial_energy, initial_sigma_energy, n_slice and slice_spacing
# are in both.
CASE = ("wavelength", "current", "current_uniform", "bunch_charge", "bunch_length_rms",
        "norm_emit_x", "norm_emit_y", "initial_sigma_x", "initial_sigma_y", "initial_beta_x",
        "initial_beta_y", "initial_alpha_x", "initial_alpha_y", "initial_sigma_pz", "n_undulator",
        "n_undulator_helical",
        "undulator_period", "undulator_period_min", "undulator_aw", "undulator_aw_min",
        "undulator_length", "undulator_length_min", "n_quadrupole", "quadrupole_k1_min",
        "quadrupole_k1_max", "line_length", "tracked_length", "rho_1d", "gain_length_1d")


def sig(x, n=3):
    """x at n significant figures in plain decimal, 0.100 and 13.7 and 132."""
    if x == 0:
        return "0"
    d = n - 1 - math.floor(math.log10(abs(x)))
    y = round(x, d)
    d = n - 1 - math.floor(math.log10(abs(y)))       # 9.996 rounds to 10.0, one digit fewer
    return f"{y:.{max(d, 0)}f}"


def scaled(x, units):
    """x in the first unit of (factor, name) pairs, largest first, in which it is 1 or more."""
    for factor, name in units:
        if abs(x) >= factor:
            return f"{sig(x / factor)} {name}"
    factor, name = units[-1]
    return f"{sig(x / factor)} {name}"


def power_of_ten(x):
    mant, expo = f"{x:.1e}".split("e")
    return f"${mant}\\times10^{{{int(expo)}}}$"


METRE, MILLI, MICRO, NANO = ((1.0, "m"),), ((1e-3, "mm"),), ((1e-6, "µm"),), ((1e-9, "nm"),)
SHORT = ((1e-6, "µm"), (1e-9, "nm"))


def case_cells(q, lattice):
    """The rows of the table of the case for one input file, as (quantity, cell) pairs, from
    its rows q of the record, {quantity: (value, unit)}, and the lattice file it names."""
    def v(name):
        return float(q[name][0]) if name in q else None

    def pair(a, b, units):
        if v(a) is None or v(b) is None:
            return None
        first, second = scaled(v(a), units), scaled(v(b), units)
        return f"{first.rsplit(' ', 1)[0]} and {second}" if first.rsplit(" ", 1)[1] == second.rsplit(" ", 1)[1] \
            else f"{first} and {second}"

    def spread(name, fmt):
        lo, hi = v(name + "_min"), v(name)
        return fmt(hi) if lo is None else f"{fmt(lo)} to {fmt(hi)}"

    out = [("lattice file", f"`{lattice}`" if lattice else None)]
    out.append(("radiation wavelength", scaled(v("wavelength"), NANO) if "wavelength" in q else None))
    out.append(("beam energy", scaled(v("initial_energy"), ((1e9, "GeV"), (1e6, "MeV"))) if "initial_energy" in q else None))
    if "current" in q:
        kind = "uniform along the window" if q["current_uniform"][0] == "yes" else "at the peak of the bunch"
        if "bunch_charge" not in q:
            kind = "steady state"
        out.append(("current", f"{scaled(v('current'), ((1e3, 'kA'), (1.0, 'A')))}, {kind}"))
    out.append(("charge in the window", scaled(v("bunch_charge"), ((1e-9, "nC"), (1e-12, "pC"), (1e-15, "fC")))
                if "bunch_charge" in q else None))
    out.append(("rms length of the current profile", scaled(v("bunch_length_rms"), SHORT) if "bunch_length_rms" in q else None))
    out.append(("normalized emittance, x and y", pair("norm_emit_x", "norm_emit_y", MICRO)))
    if "initial_sigma_pz" in q:
        out.append(("rms relative momentum spread of a slice",
                    power_of_ten(v("initial_sigma_pz")) if v("initial_sigma_pz") > 0 else "0"))
    out.append(("rms beam size at the entrance, x and y", pair("initial_sigma_x", "initial_sigma_y", MICRO)))
    out.append(("beta function at the entrance, x and y", pair("initial_beta_x", "initial_beta_y", METRE)))
    if "initial_alpha_x" in q and "initial_alpha_y" in q:
        out.append(("alpha function at the entrance, x and y", f"${v('initial_alpha_x'):.2f}$ and ${v('initial_alpha_y'):.2f}$"))
    if "n_undulator" in q:
        n, hel = int(v("n_undulator")), int(v("n_undulator_helical"))
        kind = "helical" if hel == n else "planar" if hel == 0 else f"{hel} helical and {n - hel} planar"
        out.append(("undulator", f"{kind}, period {spread('undulator_period', lambda x: scaled(x, MILLI))}, "
                                 f"rms parameter $a_w$ = {spread('undulator_aw', lambda x: f'{x:.4f}')}"))
        length = spread("undulator_length", lambda x: scaled(x, METRE))
        out.append(("undulator segments", f"{n} of {length}" if n > 1 else f"one of {length}"))
    if "n_quadrupole" in q:
        out.append(("quadrupoles", f"{int(v('n_quadrupole'))}, $k_1$ from ${v('quadrupole_k1_min'):.3g}$ "
                                   f"to ${v('quadrupole_k1_max'):.3g}$ m$^{{-2}}$"))
    out.append(("length of the line", scaled(v("line_length"), METRE) if "line_length" in q else None))
    if "tracked_length" in q and v("tracked_length") < v("line_length") * (1 - 1e-6):
        out.append(("length tracked, from the entrance", scaled(v("tracked_length"), METRE)))
    if "n_slice" in q and "slice_spacing" in q and "wavelength" in q:
        n, per = int(v("n_slice")), v("slice_spacing") / v("wavelength")
        # The record holds five significant digits of the spacing and of the wavelength.
        per = f"{round(per)}" if abs(per - round(per)) < 1e-3 * per else sig(per)
        out.append(("window", f"{n} slices of {per} wavelengths" if n > 1 else f"one slice of {per} wavelengths"
                    if per != "1" else "one slice of one wavelength"))
    if "rho_1d" in q:
        out.append(("FEL parameter $\\rho$ of the one-dimensional theory", power_of_ten(v("rho_1d"))))
        out.append(("gain length of the one-dimensional theory", scaled(v("gain_length_1d"), METRE)))
    return [(name, cell) for name, cell in out if cell is not None]


def lattice_of(path):
    """The lattice file an input file names, from the directory above examples/."""
    m = re.search(r'^lattice\s*=\s*"([^"]+)"', path.read_text(), re.M)
    if not m:
        return ""
    parts = []
    for part in (path.parent.name, *m.group(1).split("/")):
        parts = parts[:-1] if part == ".." else parts + [part]
    return "/".join(["examples", *parts])


def case_table(page, d, inputs):
    """The table of the case of one example: a row per quantity and a column per input
    file, or one column where the input files do not differ."""
    names = sorted(inputs, key=lambda n: (n != "lucifer.toml", n))
    cols = {n: dict(case_cells(inputs[n], lattice_of(d / n))) for n in names}
    order = list(dict.fromkeys(k for n in names for k in cols[n]))
    if not order:
        return
    same = all(cols[n] == cols[names[0]] for n in names)
    shown = names[:1] if same else names
    head = "value" if same else None
    page += [f":::{{table}} The case of this example: its lattice and its beam, from the statistics file of "
             f"{'its run' if len(names) == 1 else 'each of its runs'}. The emittance, the beam size and the "
             "Twiss functions are those of the slices at the first record, averaged over the slices by "
             "their current. They are moments of the macroparticles that were loaded, and so differ from the "
             "values of the input file by the statistics of the load." + (" The input files of the example do not differ in these quantities." if same and len(names) > 1 else ""),
             f":label: tab-case-{d.name.replace('_', '-')}", "",
             "| quantity | " + " | ".join(head or f"`{n}`" for n in shown) + " |",
             "|---|" + "---|" * len(shown)]
    for k in order:
        page.append(f"| {k} | " + " | ".join(cols[n].get(k, "") for n in shown) + " |")
    page += [":::", ""]


def measured_table(page, inputs):
    """One table for an example: a row per quantity, a column per input file, and the unit."""
    names = sorted(inputs, key=lambda d: (d != "lucifer.toml", d))
    quantities = [q for q in ORDER if any(q in inputs[d] for d in names)]
    quantities += sorted(q for d in names for q in inputs[d] if q not in ORDER and q not in CASE)
    page += ["## Measured on this run", "",
             "The runner wrote these from each run's stats file and log, and the keystone "
             "regenerates them. A float carries five significant digits.", "",
             "| Quantity | " + " | ".join(f"`{d}`" for d in names) + " | Unit |",
             "|---|" + "---|" * len(names) + "---|"]
    for q in dict.fromkeys(quantities):
        cells = [inputs[d].get(q, ("", ""))[0] for d in names]
        unit = next((inputs[d][q][1] for d in names if q in inputs[d]), "")
        word = WORDS.get(q) or (f"exit power, the {q[11:]} component" if q.startswith("exit_power_") else q)
        page.append(f"| {word} | " + " | ".join(cells) + f" | {unit} |")
    page.append("")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--examples", required=True, help="The examples directory")
    p.add_argument("--out", required=True, help="Directory to write the pages into")
    args = p.parse_args()

    ex = Path(args.examples)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    dirs = sorted(d for d in ex.iterdir() if d.is_dir() and (d / "README.md").is_file())
    names = {d.name for d in dirs}
    if not dirs:
        raise SystemExit(f"{ex}: no example directory carries a README.md")

    written = []
    record = read_record(ex / "measured.tsv")

    # The index page, from the examples index, plus the lattice every example shares.
    idx = ex / "README.md"
    title, body = split_title(idx.read_text(), idx)
    page = front_matter(title, "All examples", f"examples/{idx.name}", "reference")
    page += [label_sections(rewrite_links(body, names, idx)).rstrip("\n"), "",
            "## Where the numbers come from", "",
            "Every table headed \"Measured on this run\" on these pages is written by the runner, "
            "`examples/run_examples.sh`, from the stats file and the log of each run it makes, into "
            "`examples/measured.tsv`. The keystone regenerates that file and these pages and "
            "requires no diff, so a number on a page is the number the current binary gives. The table "
            "of the case that opens each page, the lattice and the beam of the example, is written from "
            "the same file.", "",
            "## The shared lattice", "",
            "Most examples call this line. The ones that need a different lattice carry "
            "their own, shown on their own page.", ""]
    embed(page, ex / "aramis.bmad")
    (out / INDEX_PAGE).write_text("\n".join(page).rstrip("\n") + "\n")
    written.append((INDEX_PAGE[:-3], 1))

    for d in dirs:
        readme = d / "README.md"
        title, body = split_title(readme.read_text(), readme)
        files = gather(d)
        if not files:
            raise SystemExit(f"{d}: no input file to embed")

        page = front_matter(title, d.name, f"examples/{d.name}/", "how-to, record")
        # The figure the runner saved, referenced beside the page it belongs to.
        if d.name not in NO_FIGURE:
            page += [f"![{d.name}]({d.name}.svg)", ""]
        # The commands below carry paths relative to the example's own directory, which a
        # reader of the repository knows and a reader of this page has not been told.
        page += [f"Everything on this page lives in and runs from "
                 f"`examples/{d.name}/`.", ""]
        if d.name in record:
            case_table(page, d, record[d.name])
        page += [rewrite_links(body, names, readme).rstrip("\n"), ""]
        if d.name in record:
            measured_table(page, record[d.name])
        page += ["## The input files", "",
                 "The input file, its variants, and any lattice they name, as they are on disk.", ""]
        for f in files:
            embed(page, f)
        (out / f"{d.name}.md").write_text("\n".join(page).rstrip("\n") + "\n")
        written.append((d.name, len(files)))

    print(f"wrote {len(written)} pages to {out}")
    for name, n in written:
        print(f"  {name:20s} {n} file(s) embedded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
