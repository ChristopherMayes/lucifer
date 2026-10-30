#!/usr/bin/env python3
"""
Draw each Genesis4 tier's figure from a kept work directory of the validation suite, and
write the page of figures and captions that docs/validation.md includes.

  report_tier_figures.py --work-dir <kept work directory> --results <fel-production.txt>
                         --out docs/generated/validation --page docs/generated/validation-figures.md
                         [--included-by docs/validation.md]

The page is written to be included by another, docs/validation.md, and MyST resolves an
included page's relative paths against the page that includes it, so the figure paths are
written relative to --included-by's directory.

Each figure is lucifer.plot.compare's, the run over the reference with the difference
beneath each curve, drawn from the tier's stats file and the Genesis4 .out.h5 it was
compared against, with the two final field dumps for the field row. weight_split is
Lucifer over Lucifer, tier1 with its weights split over tier1. The caption of each
figure states what the tier runs, what each block shows, the largest difference in each
block at two significant figures, the tier's level from the results file against its
tolerance from compare_fel.py's table, and which quantity sets the level: the power, the
bunching, the complex field over the whole exit grid, or the final particles'
coordinates, which the figure does not draw. The numbers come from lucifer.compare, the
same arrays the figure draws, so a caption cannot disagree with its figure.

The figures are SVG, so they stay sharp on a high-density display and under the reader's
zoom, at 167 to 219 kB each, near the 143 to 183 kB of a 110 dpi PNG of the same figure.
matplotlib writes a date and random element ids into an SVG unless told otherwise, so
the date is left out and the ids are salted with a fixed string: the same arrays give
the same bytes, and a keystone that redraws the figures leaves the tree as it found it.
Each figure takes the site's col-body-outset-right class, which widens it 3rem past the
text column to the right. The wider col-page and col-page-right reach under the left
navigation and under the sticky outline at the right of the book theme.

The keystone's regeneration runs this from the production pass's kept work directory and
diffs the page with the other generated pages, so a tier whose curves move fails the
keystone through its caption. The figures stay out of the diff, as the example figures
do, since a plotting library upgrade rewrites every one of them with no physics moving.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fel import TOLERANCE  # noqa: E402

# Each tier: its run's out_root, its reference's root, the reference's kind, and what it
# runs, for the caption's first sentence.
TIERS = (
    ("tier1", "tier1", "Aramis1seg", "genesis4",
     "one 4 m undulator segment of the Aramis line in steady state from a 5 kW seed, with no "
     "interlude, so the FEL core alone: the push, the deposition and the field solve"),
    ("tier1_unavg", "tier1u", "Aramis1seg", "genesis4",
     "the same segment with the unaveraged dynamics, whose model differs from Genesis4's in its "
     "sin² end ramps, its lack of period averaging and its integrator, so the tier holds a model "
     "difference and shows a drift"),
    ("tier2_genesis", "tier2g", "Aramis", "genesis4",
     "the full Aramis line, six FODO cells over 57 m, in steady state from a 5 kW seed through "
     "saturation, with the interludes between the segments transcribed from Genesis4"),
    ("tier2_bmad", "tier2", "Aramis", "genesis4",
     "the full line with Bmad tracking the interludes. Genesis4 samples the quadrupoles' path-length "
     "term once at mid element and Bmad integrates it, a transport model difference the page locates "
     "and prices below"),
    ("td1", "td1", "AramisTD1seg", "genesis4",
     "one segment time dependent, 32 slices of three wavelengths from a 5 kW seed: the FEL core with "
     "slippage"),
    ("td2_genesis", "td2g", "AramisTD", "genesis4",
     "the full line time dependent, 32 slices from a 5 kW seed, with the interludes transcribed"),
    ("td2_bmad", "td2", "AramisTD", "genesis4",
     "the full line time dependent with Bmad tracking the interludes, the tier2_bmad model difference "
     "with slippage interleaved"),
    ("tdsase", "tdsase", "AramisTDSASE", "genesis4",
     "the full line from a dark start with shot noise on, pure SASE over 32 slices, both codes "
     "tracking the identical noisy beam from shared dumps"),
    ("tdsc", "tdsc", "AramisTDSC", "genesis4",
     "one segment time dependent with space charge on, short range and long range"),
    ("tdwk", "tdwk", "AramisTDWK", "genesis4",
     "one segment time dependent with every wake kernel on"),
    ("weight_split", "tier1s", "tier1", "lucifer",
     "tier1 run again with the particles' weights split nonuniformly, drawn over tier1 itself, "
     "Lucifer over Lucifer. A collective observable is linear in the weights, so the two must agree "
     "to round-off, which Genesis4's dumps cannot test since they carry no weights"),
)
# The steady-state tiers compare their final fields over the whole grid. weight_split
# has no particle comparison.
NO_PARTICLES = ("weight_split",)


def sci(x, digits):
    """x in scientific notation with digits significant figures and a plain exponent,
    1.77e-5 rather than 1.77e-05, and 0 for zero."""
    if x == 0:
        return "0"
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m}e{int(e)}"


def read_levels(path):
    """{tier: level} from a results file of run_fel_validation.sh --results."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("tier|"):
            _, key, level, *_ = line.split("|")
            out[key] = float(level)
    return out


def caption(key, what, run, ref, diffs, level, tolerance, flevel, has_particles, image):
    td = run.n_slice > 1
    s_exit = f"{run.s[-1]:.0f} m"
    parts = [f"{key}: {what}. {ref.name} is the solid line and {run.name} the dashed line over it, "
             f"and δ is {run.name}'s value less {ref.name}'s."]
    if td:
        parts.append("(a) The energy in the window on a linear and a log axis, and under it the power "
                     "difference at each record over the reference's power, the one over the slices of "
                     "largest magnitude.")
        parts.append("(b) The slice-mean bunching, and the largest slice difference over the peak bunching.")
        parts.append("(c) The slice-mean energy change from the first record, and the largest slice "
                     "difference over the largest change.")
        parts.append(f"(d) The power across the window at the exit, s = {s_exit}, slice by slice, and its "
                     "difference over the peak exit power.")
        d_key, d_what = "exit_power", "the exit power"
    else:
        parts.append("(a) The power on a linear and a log axis, and δP/P at each record.")
        parts.append("(b) The bunching, and its difference over the peak bunching.")
        parts.append("(c) The mean energy change from the first record, and its difference over the "
                     "largest change.")
        parts.append(f"(d) The field amplitude along the row through the center of the exit plane at "
                     f"s = {s_exit}, and its difference over the peak amplitude.")
        d_key, d_what = "field", "the field row"
    parts.append(f"The largest differences are {sci(diffs['power'].largest, 2)} in power, "
                 f"{sci(diffs['bunching'].largest, 2)} in bunching, {sci(diffs['energy'].largest, 2)} in "
                 f"energy and {sci(diffs[d_key].largest, 2)} in {d_what}.")
    # Which quantity the level is: the one of compare_fel.py's that equals it.
    candidates = [("the power", diffs["power"].largest), ("the bunching", diffs["bunching"].largest),
                  (f"the complex field at the exit over the whole grid, {sci(flevel, 3)} of its peak", flevel)]
    setter = next((name for name, v in candidates if abs(v - level) <= 1e-3 * level), None)
    if setter is None:
        if not has_particles:
            raise SystemExit(f"{key}: the level {level:.6e} equals none of the power, bunching or field "
                             f"differences, {[v for _, v in candidates]}, and the tier compares no particles")
        setter = "the final particles' coordinates, which the figure does not draw"
    parts.append(f"The tier's level is {sci(level, 3)} against a tolerance of {sci(tolerance, 1)}, set by "
                 f"{setter}.")
    parts.append(f"[Full size]({image}).")
    return " ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--work-dir", required=True, help="a kept work directory of run_fel_validation.sh")
    ap.add_argument("--results", required=True, help="the pass's results file, for each tier's level")
    ap.add_argument("--out", required=True, help="where the figures are written, <out>/<tier>.svg")
    ap.add_argument("--page", required=True, help="the Markdown page of figures and captions to write")
    ap.add_argument("--included-by", default="docs/validation.md",
                    help="the page that includes --page, whose directory the figure paths are relative to")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # The SVG's element ids from a fixed salt, so a figure is the same bytes on every run.
    matplotlib.rcParams["svg.hashsalt"] = "lucifer-tier-figures"
    from lucifer import compare as cmp
    from lucifer.plot import compare

    w = Path(args.work_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    link_dir = Path(args.included_by).resolve().parent
    levels = read_levels(args.results)
    missing = [k for k, *_ in TIERS if k not in levels]
    if missing:
        raise SystemExit(f"{args.results} records no level for {', '.join(missing)}")

    lines = ["<!-- Generated by tests/scripts/report_tier_figures.py through lucifer.plot.compare from a",
             "     validation pass's work directory. Do not edit. The figures beside it are redrawn by the",
             "     same run; the captions' numbers are computed from the arrays the figures draw. -->", ""]
    for key, root, ref_root, kind, what in TIERS:
        run = cmp.read_curves(w / f"{root}.stats.h5")
        if kind == "genesis4":
            ref = cmp.read_genesis_out(w / f"{ref_root}.out.h5")
            run.name = "Lucifer"
            fields = (w / f"{root}-final.wf.h5", w / f"{ref_root}-final.fld.h5")
        else:
            ref = cmp.read_curves(w / f"{ref_root}.stats.h5")
            run.name, ref.name = f"{ref_root}, weights split", ref_root
            fields = (w / f"{root}-final.wf.h5", w / f"{ref_root}-final.wf.h5")
        reason = cmp.aligned(run, ref)
        if reason:
            raise SystemExit(f"{key}: the two runs cannot be differenced: {reason}")
        diffs = cmp.differences(run, ref, fields=fields if run.n_slice == 1 else None)
        flevel = cmp.field_level(fields, run.n_slice)
        fig = compare(run, ref, fields=fields, labels=(run.name, ref.name), return_figure=True)
        image = out / f"{key}.svg"
        fig.savefig(image, metadata={"Date": None})
        plt.close(fig)
        rel = image.resolve().relative_to(link_dir).as_posix()
        text = caption(key, what, run, ref, diffs, levels[key], TOLERANCE[key], flevel,
                       key not in NO_PARTICLES, rel)
        lines += [f":::{{figure}} {rel}", f":label: fig-{key.replace('_', '-')}", ":class: col-body-outset-right", "",
                  text, ":::", ""]
        print(f"{key}: {image.name}, level {sci(levels[key], 3)}")
    Path(args.page).write_text("\n".join(lines))
    print(f"wrote {args.page}: {len(TIERS)} figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
