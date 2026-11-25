#!/usr/bin/env python3
"""
Draw docs/figures/migration-radiation.svg and write docs/generated/migration-radiation.md, the
figure and the table of the radiation that slice migration adds to a load of beamlets.

  python docs/figures/migration_radiation_figure.py [--no-figure]
         [--tsv <file>] [--out <image>] [--page <page>]

Both are made from docs/figures/migration-radiation.tsv, which holds the measurements, and
the script holds none. A row of that file is the power of a run with migration on and a
load without noise, as a fraction of the SASE power of its configuration, at one number N
of macroparticles per slice.

The figure draws that fraction against N on logarithmic axes. A range over positions is a
bar and a single value is a point with its position beside it. The two lines are C / N.
Each C is fitted here: the geometric mean, over the rows of the example sase, which runs the
Aramis line, with beamlets of FIT_SIZE, of one end of the range multiplied by N. A row is
named by its FEL line and its radiation wavelength, which the file states. The rows that are not of that
example and that size are drawn beside the bar of their N and not on it, displaced by the
factor DODGE, so that marks at one N do not cover each other.

The page holds the table of the rows and the sentence that states the two coefficients,
each between the markers that docs/migration.md includes it by.

The keystone's regeneration runs the script with --no-figure and compares the page with
the committed one, so that a row of the file that changes fails it. The figure is drawn by
hand, without that option, and is committed with its line of docs/figures/sources.tsv,
which holds the hash of this script, of the file of measurements and of the style module.
"""
import argparse
import math
from pathlib import Path

FIT_CONFIGURATION = "sase"      # the rows the two coefficients are fitted to
FIT_SIZE = 4
DODGE = 1.16                    # the horizontal factor between marks at one N
SIZE = (6.4, 4.2)


def read_rows(tsv):
    """The rows of the file as dictionaries. A fraction and a position keep the text they
    were recorded with, which the table prints, beside the number the figure draws."""
    lines = [x for x in Path(tsv).read_text(encoding="utf-8").split("\n") if x and not x.startswith("#")]
    head = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        r = dict(zip(head, line.split("\t")))
        if len(r) != len(head):
            raise SystemExit(f"{tsv}: a row without {len(head)} columns: {line}")
        r["n"] = int(r["n_per_slice"])
        r["size"] = int(r["beamlet_size"])
        r["low"], r["high"] = float(r["fraction_low"]), float(r["fraction_high"])
        rows.append(r)
    return rows


def fit_rows(rows):
    return [r for r in rows if r["configuration"] == FIT_CONFIGURATION and r["size"] == FIT_SIZE]


def coefficient(rows, end):
    """(C, least product, largest product) of one end of the range, low or high: C is the
    geometric mean of that end multiplied by N over the fitted rows."""
    products = [r[end] * r["n"] for r in fit_rows(rows)]
    c = math.exp(sum(math.log(p) for p in products) / len(products))
    return c, min(products), max(products)


def case(r):
    """The FEL line of a row with its radiation wavelength, as the legend and the table name it."""
    return f"{r['line']} at {r['wavelength_nm']} nm"


def position(r):
    if r["s_from_m"] == r["s_to_m"]:
        return f"{r['s_from_m']} m"
    return f"{r['s_from_m']} to {r['s_to_m']} m"


def block(name, lines):
    """One block of the page, between the markers docs/migration.md includes it by."""
    return [f"<!-- begin {name} -->", *lines, f"<!-- end {name} -->", ""]


def page(rows):
    table = ["| configuration | beamlet size | macroparticles per slice | position | "
             "power with migration on, as a fraction of the SASE power |", "|---|---|---|---|---|"]
    for r in rows:
        fraction = r["fraction_low"] if r["fraction_low"] == r["fraction_high"] else \
            f"{r['fraction_low']} to {r['fraction_high']}"
        name = case(r) + (f", the example `{r['configuration']}`" if r["configuration"] != r["line"] else "")
        table.append(f"| {name} | {r['size']} | {r['n']} | "
                     f"{position(r)} | {fraction} |")
    (c_low, low_min, low_max), (c_high, high_min, high_max) = coefficient(rows, "low"), coefficient(rows, "high")
    n_fit = len(fit_rows(rows))
    fit = [f"The two lines of [](#fig-mig-radiation) are $C/N$ with $C = {c_low:.0f}$ and $C = {c_high:.0f}$. "
           f"The first coefficient is the geometric mean, over the {n_fit} rows of the "
           f"{fit_rows(rows)[0]['line']} line with beamlets of {FIT_SIZE}, of the lower end of the range multiplied by "
           f"the number of macroparticles per slice, and the second is that of the upper end. "
           f"The {n_fit} products lie between {low_min:.0f} and {low_max:.0f} for the lower end "
           f"and between {high_min:.0f} and {high_max:.0f} for the upper end."]
    out = ["<!-- Generated by docs/figures/migration_radiation_figure.py from "
           "docs/figures/migration-radiation.tsv. Do not edit. -->", ""]
    out += block("tab-mig-radiation", table)
    out += block("fit-mig-radiation", fit)
    return "\n".join(out)


def draw(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import lucifer.plot as lp

    c_low, c_high = coefficient(rows, "low")[0], coefficient(rows, "high")[0]
    counts = sorted({r["n"] for r in rows})
    with lp._style():
        fig, ax = plt.subplots(figsize=SIZE, layout="constrained")
        span = [counts[0] / 1.3, counts[-1] * 1.45]
        ax.plot(span, [c_high / x for x in span], color=lp.INK2, lw=1.0, label=f"${c_high:.0f}/N$")
        ax.plot(span, [c_low / x for x in span], color=lp.INK2, lw=1.0, ls="--", label=f"${c_low:.0f}/N$")

        def bars(size, factor, color):
            sel = [r for r in rows if r["configuration"] == FIT_CONFIGURATION and r["size"] == size]
            ax.vlines([r["n"] * factor for r in sel], [r["low"] for r in sel], [r["high"] for r in sel],
                      color=color, lw=5.0, label=f"{case(sel[0])}, beamlet_size = {size}")

        bars(FIT_SIZE, 1.0, lp.COLORS[0])
        for size in sorted({r["size"] for r in rows} - {FIT_SIZE}):
            bars(size, DODGE, lp.COLORS[1])

        # A configuration with single values: points on one side of the bar of their N,
        # each with its position on the side away from the bar.
        others = [c for c in dict.fromkeys(r["configuration"] for r in rows) if c != FIT_CONFIGURATION]
        for k, name in enumerate(others):
            sel = [r for r in rows if r["configuration"] == name]
            factor = DODGE if k == 0 else 1 / DODGE
            ax.plot([r["n"] * factor for r in sel], [r["low"] for r in sel], "sDo^"[k], ms=5.0,
                    color=lp.COLORS[2 + k], label=case(sel[0]), ls="none")
            for r in sel:
                ax.annotate(position(r), (r["n"] * factor, r["low"]), xytext=(7 if k == 0 else -7, 0),
                            textcoords="offset points", ha="left" if k == 0 else "right", va="center",
                            fontsize=8, color=lp.INK2)

        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlim(*span)
        ax.set_xticks(counts, [str(n) for n in counts])
        ax.minorticks_off()
        ax.set_xlabel("$N$")
        ax.set_ylabel(r"$P_\mathrm{m} / P_\mathrm{SASE}$")
        ax.legend(frameon=False, loc="lower left", fontsize=9)
        lp.save(fig, out)
        plt.close(fig)


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tsv", default=str(here / "migration-radiation.tsv"))
    ap.add_argument("--out", default=str(here / "migration-radiation.svg"))
    ap.add_argument("--page", default=str(here.parent / "generated" / "migration-radiation.md"))
    ap.add_argument("--no-figure", action="store_true", help="write the page and draw nothing")
    args = ap.parse_args()
    rows = read_rows(args.tsv)
    Path(args.page).write_text(page(rows), encoding="utf-8")
    if not args.no_figure:
        draw(rows, args.out)


if __name__ == "__main__":
    main()
