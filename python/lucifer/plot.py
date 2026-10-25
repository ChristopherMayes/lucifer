"""Plots of a stats file against s, with the lattice layout under them.

    from lucifer.plot import plot, summary
    plot("run.stats.h5", y="field/total/power", yscale="log")
    plot("run.stats.h5", y=["beam_size_x", "beam_size_y"], y2="beam/slice/bunching")
    summary("run.stats.h5")

    python -m lucifer.plot run.stats.h5 -o run.png

The file describes itself, so a plot needs no table of names: each dataset carries its
@long_name, its @unit and its @axes, and the lattice group carries each element's
extent and strength. A curve is named in one of three ways:

  a dataset path            field/total/power
  a path and components     beam/bunch/twiss/norm_emit[x], beam/slice/sigma[x,x]
  a derived quantity        beam_size_x, one of the names in DERIVED below

A component is named by the labels of the dataset's own axes, coords/plane, coords/bmad
or coords/wavefront, and its unit comes from the file through Stats.unit_of. A dataset on
the record and slice axes is drawn as a thin gray line per slice with the window total
bold for watts, joules and coulombs and the slice mean bold for everything else. A
dataset on the element-end axis is drawn as points at the element ends. The values are
scaled by the SI prefix beamphysics.units chooses where the unit takes one, and a unit
with a power or a quotient keeps its stated form, since a prefix on m^2 would misstate
the scale.

The derived quantities are the ones the summary figure needs and the file states only
inside a matrix, each with its formula:

  beam_size_x, beam_size_y      sqrt(beam/slice/sigma[x,x]), and [y,y]
  field_size_x, field_size_y    sqrt(field/x/sigma[x,x]), and [y,y]
  norm_emit_x, norm_emit_y      the projected normalized emittance per slice with the
                                dispersion removed, Bmad's projected_twiss_calc, from
                                beam/slice/sigma, beam/slice/centroid and run/p0c
  energy_change                 beam/slice/energy less its value at the first record
  field_emit_x, field_emit_y    field/x/emit_x and emit_y where field/x/angle_moments_valid

matplotlib and beamphysics.units are imported inside the functions, so importing this
module costs nothing beyond numpy. They come with the extra lucifer[plot].
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys
import warnings

import numpy as np

__all__ = ["plot", "summary", "layout", "PlotError", "main"]

M_ELECTRON = 0.51099895069e6   # Bmad's m_electron [eV]

# Series colors from a categorical palette, assigned in order. Per-slice context lines are
# gray, so the color stays with the bold aggregate.
COLORS = ("#2a78d6", "#eb6834", "#2b9f6a", "#8a5cc2", "#c23b4a", "#8c6d3f")
GRAY, INK, INK2 = "#b0afa9", "#0b0b0b", "#52514e"
# The figures' style, applied with matplotlib's rc_context around each figure, so the
# user's own settings are never changed.
STYLE = {
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": INK2, "grid.alpha": 0.18, "grid.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "font.size": 10,
}
# The units summed over the window rather than averaged.
EXTENSIVE = ("W", "J", "C")
# The axes a curve runs along. Every other axis of a dataset is a label axis, whose
# entries are named as components.
DATA_AXES = ("record", "slice", "element_end")


class PlotError(ValueError):
    """A curve the file cannot give: a dataset it does not have, a component its axes do
    not label, or curves on one axis with different units. The message names what the
    file does have."""


def _sigma_size(st, path, i):
    sig = st[path]
    return np.sqrt(np.maximum(0.0, sig[:, :, i, i]))


def _norm_emit(st, i):
    s = st["beam/slice/sigma"]
    s66 = s[:, :, 5, 5]
    j = i + 1
    with np.errstate(divide="ignore", invalid="ignore"):
        x2 = s[:, :, i, i] - np.where(s66 > 0, s[:, :, i, 5] ** 2 / s66, 0)
        xp = s[:, :, i, j] - np.where(s66 > 0, s[:, :, i, 5] * s[:, :, j, 5] / s66, 0)
        p2 = s[:, :, j, j] - np.where(s66 > 0, s[:, :, j, 5] ** 2 / s66, 0)
    emit = np.sqrt(np.maximum(0.0, x2 * p2 - xp ** 2))
    return float(st.run["p0c"]) * (1.0 + st["beam/slice/centroid"][:, :, 5]) / M_ELECTRON * emit


def _field_emit(st, plane):
    valid = st["field/x/angle_moments_valid"].astype(bool)
    return np.where(valid, st[f"field/x/emit_{plane}"], np.nan)


def _energy_change(st):
    e = st["beam/slice/energy"]
    return e - e[0]


# name: (quantity, component, unit, function of the Stats). Every one runs on the record
# and slice axes.
DERIVED = {
    "beam_size_x": ("rms beam size", "x", "m", lambda st: _sigma_size(st, "beam/slice/sigma", 0)),
    "beam_size_y": ("rms beam size", "y", "m", lambda st: _sigma_size(st, "beam/slice/sigma", 2)),
    "field_size_x": ("rms field size", "x", "m", lambda st: _sigma_size(st, "field/x/sigma", 0)),
    "field_size_y": ("rms field size", "y", "m", lambda st: _sigma_size(st, "field/x/sigma", 2)),
    "norm_emit_x": ("normalized emittance", "x", "m rad", lambda st: _norm_emit(st, 0)),
    "norm_emit_y": ("normalized emittance", "y", "m rad", lambda st: _norm_emit(st, 2)),
    "energy_change": ("energy change", None, "eV", _energy_change),
    "field_emit_x": ("field emittance", "x", "m rad", lambda st: _field_emit(st, "x")),
    "field_emit_y": ("field emittance", "y", "m rad", lambda st: _field_emit(st, "y")),
}


class _Series:
    def __init__(self, name, quantity, component, unit, s, values, reduce, at_ends):
        self.name, self.quantity, self.component, self.unit = name, quantity, component, unit
        self.s, self.values, self.reduce, self.at_ends = s, values, reduce, at_ends


def _open(stats):
    """(Stats, whether this call opened it)."""
    from lucifer.stats import Stats, read_stats
    if isinstance(stats, Stats):
        return stats, False
    return read_stats(stats), True


def _component_unit(units, power, idx):
    if len(units) == 1:
        return units[0] if power == 1 else f"{units[0]}^{power}"
    picked = [units[i] for i in idx]
    if power == 1:
        return picked[0]
    if len(set(picked)) == 1:
        return "1" if picked[0] == "1" else f"{picked[0]}^{len(picked)}"
    named = [u for u in picked if u != "1"]
    return " ".join(named) if named else "1"


def _refuse_missing(st, path):
    group = path.rsplit("/", 1)[0] if "/" in path else ""
    if group and group in st:
        have = [k for k in st._h5[group] if hasattr(st._h5[f"{group}/{k}"], "shape")]
        raise PlotError(f"{st.path.name} has no dataset {path}. {group} holds: {', '.join(have)}.")
    if "/" not in path:
        raise PlotError(f"{path} is neither a dataset path nor a derived quantity. "
                        f"The derived quantities are: {', '.join(DERIVED)}.")
    raise PlotError(f"{st.path.name} has no group {group}. Its groups are: "
                    f"{', '.join(k for k in st._h5 if not hasattr(st._h5[k], 'shape'))}.")


def _series(st, key, x, reduce):
    if key in DERIVED:
        quantity, component, unit, fn = DERIVED[key]
        values = fn(st)
        axes, path = ("record", "slice"), None
    else:
        m = re.fullmatch(r"\s*([\w/]+)\s*(?:\[([\w, ]+)\])?\s*", key)
        if not m:
            raise PlotError(f"{key} is not a dataset path, path[component] or derived quantity")
        path, sel = m.group(1), [c.strip() for c in (m.group(2) or "").split(",") if c.strip()]
        if path not in st.units:
            _refuse_missing(st, path)
        axes = st.axes[path]
        data = [a for a in axes if a in DATA_AXES]
        labels = [a for a in axes if a not in DATA_AXES]
        if not data:
            raise PlotError(f"{path} has no record or element-end axis to plot against s")
        if len(sel) != len(labels):
            named = "; ".join(f"{a}: {', '.join(st.coord(a))}" for a in labels)
            raise PlotError(f"{path} runs over {', '.join(labels) or 'no label axis'}, so it takes "
                            f"{len(labels)} component(s) in brackets, {key} names {len(sel)}"
                            + (f". The labels are {named}." if labels else "."))
        idx = []
        for a, c in zip(labels, sel):
            names = st.coord(a)
            if c not in names:
                raise PlotError(f"{path}'s axis {a} has no component {c}; its labels are {', '.join(names)}")
            idx.append(names.index(c))
        values = np.asarray(st[path], dtype=float)
        for i in idx:
            values = np.take(values, i, axis=len(data))
        units, power = st.unit_of(path)
        unit = _component_unit(units, power, idx)
        quantity = st.long_name[path] or path
        component = ",".join(sel) or None
        axes = tuple(data)
    at_ends = axes[0] == "element_end"
    if x == "s":
        s = st.s_end if at_ends else st.s
    elif x == "record":
        if at_ends:
            raise PlotError(f"{key} runs over element ends and has no record axis")
        s = st.record
    else:
        raise PlotError(f"x is s or record, not {x}")
    if values.ndim == 1:
        values = values[:, None]
    if reduce is None:
        reduce = "sum" if unit in EXTENSIVE else "mean"
    return _Series(key, quantity, component, unit, np.asarray(s, dtype=float), values, reduce, at_ends)


def _scale(unit, arrays):
    """(factor, prefix) for a unit and the values drawn in it."""
    if unit == "1" or not re.fullmatch(r"[A-Za-z]+( rad)?", unit):
        return 1.0, ""
    from beamphysics.units import nice_scale_prefix
    vals = np.concatenate([np.abs(a[np.isfinite(a)]).ravel() for a in arrays] or [np.zeros(0)])
    if vals.size == 0 or vals.max() == 0:
        return 1.0, ""
    return nice_scale_prefix(float(vals.max()))


def _unit_label(prefix, unit):
    return "" if unit == "1" else f" ({prefix}{unit})"


def _draw(ax, series, first_color, xfactor, yscale, include_legend, slices=True, named=False,
          labels=None):
    """Draw the series on one axes in one unit, and return the axis label. With slices,
    a curve over several slices draws a thin line per slice under its aggregate. With
    named, the legend names each curve by its quantity, as it must when two axes share
    one legend."""
    units = sorted({se.unit for se in series})
    if len(units) > 1:
        raise PlotError("the curves on one axis must share a unit: "
                        + ", ".join(f"{se.name} [{se.unit}]" for se in series))
    factor, prefix = _scale(units[0], [se.values for se in series])
    quantities = list(dict.fromkeys(se.quantity for se in series))
    for k, se in enumerate(series):
        color = COLORS[(first_color + k) % len(COLORS)]
        v, s = se.values / factor, se.s / xfactor
        live = np.isfinite(v).any(axis=1)
        sparse = se.at_ends or live.sum() < 0.5 * len(live)
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore", RuntimeWarning)
            if v.shape[1] > 1:
                # A curve with values at a few records alone draws its mean as points;
                # per-slice lines through those records would read as streaks.
                if slices and not sparse:
                    ax.plot(s, v, color=GRAY, lw=0.6, alpha=0.6, zorder=1)
                agg = np.nansum(v, axis=1) if se.reduce == "sum" else np.nanmean(v, axis=1)
            else:
                agg = v[:, 0]
        agg = np.where(live, agg, np.nan)
        ok = np.isfinite(agg)
        same = len(quantities) == 1 and len(series) > 1 and not any(x.component for x in series)
        if labels is not None:
            label = labels[k]
        elif same:
            # Curves of one quantity with no component, a power and a harmonic's power, are
            # told apart by their paths.
            label = se.name
        elif len(quantities) == 1 and not named:
            label = se.component or ("total" if se.reduce == "sum" else "slice mean")
            if se.component is None and len(series) == 1 and v.shape[1] == 1:
                label = se.quantity
        else:
            label = se.quantity + (f" {se.component}" if se.component else "")
        ax.plot(s[ok], agg[ok], "o-" if sparse else "-", ms=3, lw=1.8 if not sparse else 1.2,
                color=color, label=label, zorder=3)
    ax.set_yscale(yscale)
    if yscale == "log":
        # A log axis labels its decades. One that spans less than a decade would show no
        # decade, or crowd its minor ticks with digits, so it takes a few evenly spaced
        # ticks in plain numbers instead.
        from matplotlib.ticker import MaxNLocator, NullFormatter, ScalarFormatter
        ax.yaxis.set_minor_formatter(NullFormatter())
        lo, hi = ax.dataLim.y0, ax.dataLim.y1
        if lo > 0 and hi / lo < 10:
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    multi = len(series) > 1 or any(se.values.shape[1] > 1 for se in series)
    # The axis carries the first quantity and the unit, and the legend names each curve.
    return quantities[0] + _unit_label(prefix, units[0]), factor, multi


def layout(stats, ax, xfactor=1.0, xlabel=True):
    """Draw the lattice layout on ``ax`` from the stats file's lattice group: an FEL
    element as a box up from the axis scaled by its aw, an element with a k1 as a box
    scaled by k1 with its sign, and the line itself as a rule. ``xfactor`` divides s, so
    the layout matches a plot drawn with a prefix on s. Nothing is modeled: the extents
    and strengths are the file's."""
    from matplotlib.patches import Rectangle
    st, close = _open(stats)
    try:
        s0, s1 = st["lattice/s_start"] / xfactor, st["lattice/s_end"] / xfactor
        aw, k1 = st["lattice/aw"], st["lattice/k1"]
        fel = st["lattice/is_fel"].astype(bool)
        aw_max = max(float(np.max(np.abs(aw[fel]))) if fel.any() else 0.0, 1e-300)
        quad = (~fel) & (k1 != 0)
        k_max = max(float(np.max(np.abs(k1[quad]))) if quad.any() else 0.0, 1e-300)
        ax.axhline(0, color=INK2, lw=0.8, zorder=1)
        for a, b, isf, w, k in zip(s0, s1, fel, aw, k1):
            if b <= a:
                continue
            if isf:
                ax.add_patch(Rectangle((a, 0), b - a, abs(w) / aw_max, color=COLORS[0], alpha=0.8, lw=0))
            elif k != 0:
                ax.add_patch(Rectangle((a, 0), b - a, 0.6 * k / k_max, color=COLORS[1], alpha=0.8, lw=0))
        ax.set_ylim(-0.75, 1.15)
        ax.set_yticks([])
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)
        ax.grid(False)
        if xlabel:
            prefix = {1.0: ""}.get(xfactor)
            ax.set_xlabel("s" + (" (m)" if prefix == "" else ""))
    finally:
        if close:
            st.close()


def plot(stats, y="field/total/power", x="s", y2=(), xlim=None, ylim=None, ylim2=None,
         yscale="linear", yscale2="linear", include_layout=True, include_legend=True,
         return_figure=False, ax=None, reduce=None, labels=None):
    """Plot curves of a stats file against s, the form of LUME-Genesis's plot().

    ``stats`` is a Stats or a path. ``y`` and ``y2`` are a curve or a list of curves for
    the left and the right axis, each a dataset path, a path with components in
    brackets, or a derived quantity, as the module docstring describes. ``x`` is "s" or
    "record". The limits are in the file's units, before any prefix. ``reduce`` is
    "sum" or "mean" across the slices, and by default the window sum for watts, joules
    and coulombs and the mean otherwise. With ``include_layout`` the lattice layout is
    drawn under the axes. With ``ax`` the curves go on that axes and no figure is made,
    which is how summary() builds its panels. ``labels`` names the curves of ``y`` and
    ``y2`` in order in the legend, in place of the names the file gives. Returns the
    figure with ``return_figure``.
    """
    import matplotlib.pyplot as plt
    st, close = _open(stats)
    try:
        with plt.rc_context(STYLE):
            fig = _plot(plt, st, y, x, y2, xlim, ylim, ylim2, yscale, yscale2,
                        include_layout, include_legend, ax, reduce, labels)
    finally:
        if close:
            st.close()
    return fig if return_figure else None


def _plot(plt, st, y, x, y2, xlim, ylim, ylim2, yscale, yscale2, include_layout,
          include_legend, ax, reduce, labels):
    ys = [y] if isinstance(y, str) else list(y)
    y2s = [y2] if isinstance(y2, str) else list(y2)
    ser = [_series(st, k, x, reduce) for k in ys]
    ser2 = [_series(st, k, x, reduce) for k in y2s]
    if x == "s":
        xfactor, xprefix = _scale("m", [se.s for se in ser + ser2])
        xlab = f"s ({xprefix}m)"
    else:
        xfactor, xlab = 1.0, "record"
    lax = None
    if ax is None:
        if include_layout:
            fig, (ax, lax) = plt.subplots(2, 1, sharex=True, figsize=(8, 4.8), layout="constrained",
                                          gridspec_kw=dict(height_ratios=(5, 1)))
        else:
            fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    else:
        fig = ax.figure
        if include_layout:
            lax = ax.inset_axes((0, -0.3, 1, 0.2), sharex=ax)
    if labels is not None and len(labels) != len(ser) + len(ser2):
        raise PlotError(f"labels names {len(labels)} curves and y and y2 hold {len(ser) + len(ser2)}")
    label, factor, multi = _draw(ax, ser, 0, xfactor, yscale, include_legend, named=bool(ser2),
                                 labels=None if labels is None else labels[:len(ser)])
    ax.set_ylabel(label)
    if ylim is not None:
        ax.set_ylim(*(v / factor if v is not None else None for v in ylim))
    leg_h, leg_l = ax.get_legend_handles_labels()
    if ser2:
        ax2 = ax.twinx()
        # The right axis draws aggregates alone, so the two axes' slices do not overlay.
        label2, factor2, _ = _draw(ax2, ser2, len(ser), xfactor, yscale2, include_legend,
                                   slices=False, named=True,
                                   labels=None if labels is None else labels[len(ser):])
        ax2.grid(False)
        ax2.set_ylabel(label2)
        if ylim2 is not None:
            ax2.set_ylim(*(v / factor2 if v is not None else None for v in ylim2))
        h2, l2 = ax2.get_legend_handles_labels()
        leg_h, leg_l, multi = leg_h + h2, leg_l + l2, True
    if include_legend and multi:
        ax.legend(leg_h, leg_l, frameon=False)
    if xlim is not None:
        ax.set_xlim(*(v / xfactor if v is not None else None for v in xlim))
    if lax is not None:
        layout(st, lax, xfactor=xfactor, xlabel=False)
        lax.set_xlabel(xlab)
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(xlab)
    return fig


def summary(stats, return_figure=False):
    """The standard figure of a run, each panel made with plot(): radiation power and
    window field energy on log and linear scales, the power inside the split angle where
    the file has it, bunching, the energy change and spread, the rms beam and field
    sizes, the beam's normalized and the field's emittances, a panel per polarization
    when the run carries two, a panel per harmonic when it carries harmonics, and the
    layout under them all. Returns the figure with ``return_figure``."""
    import matplotlib.pyplot as plt
    st, close = _open(stats)
    try:
        with plt.rc_context(STYLE):
            fig = _summary(plt, st)
    finally:
        if close:
            st.close()
    return fig if return_figure else None


def _summary(plt, st):
    two_pol = "y" in st.components
    harm = st.harmonics
    rows = [["power", "energy"], ["power_lin", "energy_lin"], ["bunching", "gamma"],
            ["bsize", "fsize"], ["bemit", "femit"]]
    if two_pol:
        rows.append(["pol", "pol"])
    if harm:
        rows.append(["harm", "harm"])
    # A layout under each column, so each lines up with the panels above it.
    rows.append(["layout_left", "layout_right"])
    heights = [1] * (len(rows) - 1) + [0.25]
    fig, axd = plt.subplot_mosaic(rows, figsize=(9.5, 14.5 + 2.6 * (two_pol + bool(harm))),
                                  sharex=True, height_ratios=heights)
    power = ["field/total/power"]
    if "field/total/power_inside_angle" in st:
        power.append("field/total/power_inside_angle")
    panels = {
        "power": dict(y=power, yscale="log"),
        "energy": dict(y="field/total/energy", yscale="log"),
        "power_lin": dict(y=power),
        "energy_lin": dict(y="field/total/energy"),
        "bunching": dict(y="beam/slice/bunching"),
        "gamma": dict(y=["energy_change", "beam/slice/sigma_energy"]),
        "bsize": dict(y=["beam_size_x", "beam_size_y"]),
        "fsize": dict(y=["field_size_x", "field_size_y"]),
        "bemit": dict(y=["norm_emit_x", "norm_emit_y"]),
        "femit": dict(y=["field_emit_x", "field_emit_y"]),
    }
    if two_pol:
        panels["pol"] = dict(y=["field/x/power", "field/y/power"], yscale="log",
                             labels=["x polarization", "y polarization"])
    if harm:
        panels["harm"] = dict(y=["field/total/power"] + [f"field/harm{h}/total/power" for h in harm],
                              yscale="log", labels=["fundamental"] + [f"harmonic {h}" for h in harm])
    for key, kw in panels.items():
        plot(st, ax=axd[key], include_layout=False, **kw)
        axd[key].set_xlabel("")
    xfactor, xprefix = _scale("m", [st.s])
    for key in ("layout_left", "layout_right"):
        layout(st, axd[key], xfactor=xfactor, xlabel=False)
        axd[key].set_xlabel(f"s ({xprefix}m)")
    nslice = st["field/total/power"].shape[1]
    title = st.path.name + (f"   ({nslice} slices; thin lines are individual slices)" if nslice > 1 else "")
    fig.suptitle(title, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def main(argv=None) -> int:
    """python -m lucifer.plot <stats.h5> [-o out.png] [--dpi N] [--show]: write the
    summary figure of a run, by default to <out_root>.png in the current directory."""
    ap = argparse.ArgumentParser(prog="python -m lucifer.plot", description=main.__doc__.split(":")[1].strip())
    ap.add_argument("stats", help="a stats file, <out_root>.stats.h5")
    ap.add_argument("-o", "--output", help="the image to write (default: <out_root>.png)")
    ap.add_argument("--dpi", type=int, default=160, help="output resolution (default: 160)")
    ap.add_argument("--show", action="store_true", help="also open a window")
    args = ap.parse_args(argv)
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    out = args.output or pathlib.Path(args.stats).name.removesuffix(".stats.h5") + ".png"
    fig = summary(args.stats, return_figure=True)
    fig.savefig(out, dpi=args.dpi)
    print(f"wrote {out}")
    if args.show:
        import matplotlib.pyplot as plt
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
