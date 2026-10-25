#!/usr/bin/env python3
"""
Draw the figures of docs/python.md's Plotting section from the sase example's run.

  python_page_figures.py --stats-dir examples/sase --out docs/generated/python

Each region between "# example: <name>" and "# end example: <name>" is the code of one
block on the page, as a user types it, with the run's own file name, and
tests/scripts/check_python_page.py holds the page to it. After each region the figure it
drew is saved as <out>/<name>.png at one fixed resolution and checked: its labels carry
the file's units with their prefix, and the layout is drawn. A figure that lacks what
the page says it shows exits nonzero.

examples/run_examples.sh runs this after every example has run, into
docs/generated/python with figures on and into its scratch directory under
--no-figures, so the keystone's examples job runs the page's code on the real run while
the committed figures stay out of the regeneration diff, as the examples' figures do.
"""
import argparse
import contextlib
import os
import pathlib
import sys

DPI = 90


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--stats-dir", required=True, help="the sase example's directory, holding sase.stats.h5")
    ap.add_argument("--out", required=True, help="where the figures are written")
    args = ap.parse_args()
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib.pyplot as plt
    out = pathlib.Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    failed = []

    def save(name):
        fig = plt.gcf()
        fig.savefig(out / f"{name}.png", dpi=DPI)
        return fig

    def expect(name, ok, what):
        print(f"--- {name}: {what} {'ok' if ok else '** FAIL **'}")
        if not ok:
            failed.append(name)

    with contextlib.chdir(args.stats_dir):
        # example: plot-power
        from lucifer.plot import plot

        plot("sase.stats.h5", y="field/total/power", y2="beam/slice/bunching", yscale="log")
        # end example: plot-power
        fig = save("power")
        main_ax, layout_ax, twin_ax, k1_ax = fig.axes
        expect("power", main_ax.get_ylabel() == "$P$ (MW)" and twin_ax.get_ylabel() == "$|b|$"
               and len(layout_ax.patches) == 12 and len(k1_ax.patches) == 12,
               f"[{main_ax.get_ylabel()} | {twin_ax.get_ylabel()}; {len(layout_ax.patches)} undulators "
               f"and {len(k1_ax.patches)} quadrupoles drawn]")

        # example: plot-sizes
        plot("sase.stats.h5", y=["beam_size_x", "field_size_x"])
        # end example: plot-sizes
        fig = save("sizes")
        labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        expect("sizes", fig.axes[0].get_ylabel() == r"$\sigma_{ x }, \sigma_{x,\mathrm{field}}$ (µm)"
               and labels == [r"$\sigma_{ x }$", r"$\sigma_{x,\mathrm{field}}$"],
               f"[{fig.axes[0].get_ylabel()}; legend {labels}]")

        # example: plot-components
        plot("sase.stats.h5", y=["beam/bunch/twiss/norm_emit[x]", "beam/slice_twiss/twiss/norm_emit[x]"],
             labels=["projected", "slice average"])
        # end example: plot-components
        fig = save("emittance")
        expect("emittance", fig.axes[0].get_ylabel() == r"$\epsilon_{n, x}$ (nm rad)",
               f"[{fig.axes[0].get_ylabel()}]")

        # example: summary
        from lucifer.plot import summary

        summary("sase.stats.h5")
        # end example: summary
        fig = save("summary")
        expect("summary", len(fig.axes) == 14,
               f"[{len(fig.axes)} axes, ten panels and two layouts with their k1 axes]")
        plt.close("all")

    print(f"python_page_figures: {4 - len(failed)} of 4 figures written to {out}"
          + (f", failed: {', '.join(failed)}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
