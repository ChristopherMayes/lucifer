"""Plots of a stats file against s, with the lattice layout under them.

    from lucifer.plot import plot, summary, compare
    plot("run.stats.h5", y="field/total/power", yscale="log")
    plot("run.stats.h5", y=["beam_size_x", "beam_size_y"], y2="beam/slice/bunching")
    summary("run.stats.h5")
    compare("run.stats.h5", "Aramis.out.h5")

    python -m lucifer.plot run.stats.h5 -o run.png
    python -m lucifer.plot compare run.stats.h5 Aramis.out.h5 -o run-vs-genesis.png

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

compare() draws a run over a reference run, Genesis4's or another Lucifer run, with the
difference beneath each curve, from the curves lucifer.compare reads: the power on a
linear and a log axis, the bunching, the mean energy change and, in steady state, the
field amplitude across the exit plane or the rms beam sizes, and in a time-dependent run
the pulse energy in place of the power and the power across the window at the exit. The
reference is a thick solid line and the run a thinner dashed one over it. A difference
axis carries its power of ten in its label.

matplotlib and beamphysics.units are imported inside the functions, so importing this
module costs nothing beyond numpy. They come with the extra lucifer[plot].
"""
from __future__ import annotations

import argparse
import contextlib
import pathlib
import re
import sys
import warnings

import numpy as np

__all__ = ["plot", "summary", "compare", "layout", "PlotError", "main"]

M_ELECTRON = 0.51099895069e6   # Bmad's m_electron [eV]

# Series colors from a categorical palette, assigned in order. Per-slice context lines are
# gray, so the color stays with the bold aggregate.
COLORS = ("#2a78d6", "#eb6834", "#2b9f6a", "#8a5cc2", "#c23b4a", "#8c6d3f")
GRAY, INK, INK2 = "#b0afa9", "#0b0b0b", "#52514e"
# The layout's colors, LUME-Genesis's: red for an undulator's aw, blue for a quadrupole's k.
LAYOUT_AW, LAYOUT_K1 = "red", "blue"
# The figures' style, set around each figure by _style(), which restores these keys alone,
# so the user's own settings are never changed.
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


@contextlib.contextmanager
def _style():
    """Set STYLE for the figures drawn inside, and restore those keys alone after.
    matplotlib's rc_context restores every key, interactive among them, and under IPython's
    inline backend that switched interactive mode off after the first figure, so no later
    figure in a notebook was shown."""
    import matplotlib as mpl
    saved = {k: mpl.rcParams[k] for k in STYLE}
    mpl.rcParams.update(STYLE)
    try:
        yield
    finally:
        mpl.rcParams.update(saved)


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


# The symbol of a curve, for labels drawn with tex. A name that openpmd-beamphysics knows is
# rendered by its beamphysics.labels.texlabel, so a beam size reads sigma_x as it does in
# LUME-Genesis; a quantity it has no name for has its own symbol here; any other curve keeps
# the file's own name.
BEAMPHYSICS_KEY = {
    "beam_size_x": "sigma_x", "beam_size_y": "sigma_y",
    "norm_emit_x": "norm_emit_x", "norm_emit_y": "norm_emit_y",
    "beam/slice/energy": "energy", "beam/slice/sigma_energy": "sigma_energy",
}
SYMBOL = {
    "field/total/power": "P", "field/total/power_inside_angle": r"P_\mathrm{inside}",
    "field/total/energy": r"E_\mathrm{window}", "field/total/escaped_energy": r"E_\mathrm{escaped}",
    "field/total/pulse_energy": r"E_\mathrm{pulse}", "field/x/power": "P_x", "field/y/power": "P_y",
    "beam/slice/bunching": "|b|", "beam/slice/current": "I",
    "field_size_x": r"\sigma_{x,\mathrm{field}}", "field_size_y": r"\sigma_{y,\mathrm{field}}",
    "field_emit_x": r"\epsilon_{x,\mathrm{field}}", "field_emit_y": r"\epsilon_{y,\mathrm{field}}",
    "energy_change": r"\Delta\langle E\rangle",
}
# The Twiss quantities of beam/bunch and beam/slice_twiss, by their beamphysics names.
TWISS = {"norm_emit": "norm_emit_{c}", "emit": "twiss_emit_{c}", "beta": "twiss_beta_{c}",
         "alpha": "twiss_alpha_{c}", "gamma": "twiss_gamma_{c}", "eta": "twiss_eta_{c}",
         "etap": "twiss_etap_{c}", "sigma": "sigma_{c}"}


def _symbol(se):
    """The TeX symbol of a curve, without the dollar signs, or None for the file's name."""
    from beamphysics.labels import texlabel
    path = se.name.split("[")[0].strip()
    if path in BEAMPHYSICS_KEY:
        return texlabel(BEAMPHYSICS_KEY[path])
    if path in SYMBOL:
        return SYMBOL[path]
    m = re.fullmatch(r"field/harm(\d+)/total/power", path)
    if m:
        return f"P_{{{m.group(1)}}}"
    m = re.fullmatch(r"field/harm(\d+)/total/(energy|escaped_energy|pulse_energy)", path)
    if m:
        which = {"energy": "window", "escaped_energy": "escaped", "pulse_energy": "pulse"}[m.group(2)]
        return rf"E_{{\mathrm{{{which}}},{m.group(1)}}}"
    m = re.fullmatch(r"beam/(?:bunch|slice_twiss)/(?:twiss|modes)/(\w+)", path)
    if m and m.group(1) in TWISS and se.component:
        return texlabel(TWISS[m.group(1)].format(c=se.component))
    return None


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
          labels=None, tex=True):
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
    symbols = [_symbol(se) if tex else None for se in series]
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
                # Only the first curve of an axis draws its slices: two curves' slices in
                # one gray could not be told apart.
                if slices and k == 0 and not sparse:
                    ax.plot(s, v, color=GRAY, lw=0.6, alpha=0.6, zorder=1)
                agg = np.nansum(v, axis=1) if se.reduce == "sum" else np.nanmean(v, axis=1)
            else:
                agg = v[:, 0]
        agg = np.where(live, agg, np.nan)
        ok = np.isfinite(agg)
        same = len(quantities) == 1 and len(series) > 1 and not any(x.component for x in series)
        if labels is not None:
            label = labels[k]
        elif symbols[k] and (len(series) > 1 or named):
            label = f"${symbols[k]}$"
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
        _log_ticks(ax)
    multi = len(series) > 1 or any(se.values.shape[1] > 1 for se in series)
    # With tex and a symbol for every curve, the axis carries the symbols, as LUME-Genesis's
    # mathlabel does. Else it carries the quantities while their names fit beside the axis,
    # else the first. The legend names each curve.
    if all(symbols):
        name = "$" + ", ".join(dict.fromkeys(symbols)) + "$"
    else:
        joined = ", ".join(quantities)
        name = joined if len(joined) <= 32 else quantities[0]
    return name + _unit_label(prefix, units[0]), factor, multi


def _log_ticks(ax):
    """A log axis labels its decades. One that spans less than a decade would show no
    decade, or crowd its minor ticks with digits, so it takes a few evenly spaced ticks in
    plain numbers instead."""
    from matplotlib.ticker import MaxNLocator, NullFormatter, ScalarFormatter
    ax.yaxis.set_minor_formatter(NullFormatter())
    lo, hi = ax.get_ylim()
    if lo > 0 and hi / lo < 10:
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))


def layout(stats, ax, xfactor=1.0, xlabel=True, ylabels=True):
    """Draw the lattice layout on ``ax`` from the stats file's lattice group, in
    LUME-Genesis's colors: each FEL element's aw as a red box on ``ax``, and each element
    with a k1 as a blue box at its k1, with its sign, on a twin axis at the right. Each box
    has an edge in its color, so a quadrupole a pixel wide stays visible at its true
    length. ``xfactor`` divides s, so the layout matches a plot drawn with a prefix on s,
    and ``ylabels`` labels the two axes, a_w and k_1. Nothing is modeled: the extents and
    strengths are the file's."""
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FormatStrFormatter
    st, close = _open(stats)
    try:
        s0, s1 = st["lattice/s_start"] / xfactor, st["lattice/s_end"] / xfactor
        aw, k1 = st["lattice/aw"], st["lattice/k1"]
        fel = st["lattice/is_fel"].astype(bool)
        quad = (~fel) & (k1 != 0)
        aw_max = float(np.max(np.abs(aw[fel]))) if fel.any() else 1.0
        k_max = float(np.max(np.abs(k1[quad]))) if quad.any() else 1.0
        kax = ax.twinx()
        for a, b, isf, w, k, q in zip(s0, s1, fel, aw, k1, quad):
            if b <= a:
                continue
            if isf:
                ax.add_patch(Rectangle((a, 0), b - a, abs(w), facecolor=LAYOUT_AW, edgecolor=LAYOUT_AW,
                                       alpha=0.5, lw=0.8))
            elif q:
                kax.add_patch(Rectangle((a, 0), b - a, k, facecolor=LAYOUT_K1, edgecolor=LAYOUT_K1,
                                        alpha=0.5, lw=0.8))
        # The two zeros at one height, so the aw boxes and the quadrupoles share a baseline.
        ax.set_ylim(-1.15 * aw_max, 1.15 * aw_max)
        kax.set_ylim(-1.15 * k_max, 1.15 * k_max)
        for a_ in (ax, kax):
            a_.grid(False)
            a_.spines["top"].set_visible(False)
            a_.tick_params(labelsize=7)
        if ylabels:
            ax.set_ylabel("$a_w$", fontsize=8)
            kax.set_ylabel("$k_1$ (1/m$^2$)", fontsize=8)
            ax.set_yticks([0, aw_max])
            kax.set_yticks([-k_max, 0, k_max])
            for a_ in (ax, kax):
                a_.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        else:
            ax.set_yticks([])
            kax.set_yticks([])
        if xlabel:
            ax.set_xlabel("$s$ (m)" if xfactor == 1.0 else "$s$")
    finally:
        if close:
            st.close()


def plot(stats, y="field/total/power", x="s", y2=(), xlim=None, ylim=None, ylim2=None,
         yscale="linear", yscale2="linear", include_layout=True, include_legend=True,
         return_figure=False, ax=None, reduce=None, labels=None, tex=True):
    """Plot curves of a stats file against s, the form of LUME-Genesis's plot().

    ``stats`` is a Stats or a path. ``y`` and ``y2`` are a curve or a list of curves for
    the left and the right axis, each a dataset path, a path with components in
    brackets, or a derived quantity, as the module docstring describes. ``x`` is "s" or
    "record". The limits are in the file's units, before any prefix. ``reduce`` is
    "sum" or "mean" across the slices, and by default the window sum for watts, joules
    and coulombs and the mean otherwise. With ``include_layout`` the lattice layout is
    drawn under the axes. With ``ax`` the curves go on that axes and no figure is made,
    which is how summary() builds its panels. ``labels`` names the curves of ``y`` and
    ``y2`` in order in the legend, in place of the names the file gives. With ``tex`` a
    curve with a symbol is labeled by it, through beamphysics.labels.texlabel where
    openpmd-beamphysics has a name for it, so the beam sizes read sigma_x and sigma_y, and
    without it every curve keeps the file's name. Returns the figure with
    ``return_figure``.
    """
    import matplotlib.pyplot as plt
    st, close = _open(stats)
    try:
        with _style():
            fig = _plot(plt, st, y, x, y2, xlim, ylim, ylim2, yscale, yscale2,
                        include_layout, include_legend, ax, reduce, labels, tex)
    finally:
        if close:
            st.close()
    return fig if return_figure else None


def _plot(plt, st, y, x, y2, xlim, ylim, ylim2, yscale, yscale2, include_layout,
          include_legend, ax, reduce, labels, tex):
    ys = [y] if isinstance(y, str) else list(y)
    y2s = [y2] if isinstance(y2, str) else list(y2)
    ser = [_series(st, k, x, reduce) for k in ys]
    ser2 = [_series(st, k, x, reduce) for k in y2s]
    if x == "s":
        xfactor, xprefix = _scale("m", [se.s for se in ser + ser2])
        xlab = f"$s$ ({xprefix}m)" if tex else f"s ({xprefix}m)"
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
                                 labels=None if labels is None else labels[:len(ser)], tex=tex)
    ax.set_ylabel(label)
    if ylim is not None:
        ax.set_ylim(*(v / factor if v is not None else None for v in ylim))
    leg_h, leg_l = ax.get_legend_handles_labels()
    if ser2:
        ax2 = ax.twinx()
        # The right axis draws aggregates alone, so the two axes' slices do not overlay.
        label2, factor2, _ = _draw(ax2, ser2, len(ser), xfactor, yscale2, include_legend,
                                   slices=False, named=True,
                                   labels=None if labels is None else labels[len(ser):], tex=tex)
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
    elif lax is not None and x == "s":
        # With the layout drawn, s covers the whole line, so curves that start at the
        # first element end still sit over the elements they describe.
        lo = min(float(np.min(st["lattice/s_start"])), *(float(np.min(se.s)) for se in ser + ser2))
        hi = max(float(np.max(st["lattice/s_end"])), *(float(np.max(se.s)) for se in ser + ser2))
        pad = 0.02 * (hi - lo)
        ax.set_xlim((lo - pad) / xfactor, (hi + pad) / xfactor)
    if lax is not None:
        layout(st, lax, xfactor=xfactor, xlabel=False)
        lax.set_xlabel(xlab)
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(xlab)
    return fig


def summary(stats, return_figure=False, tex=True):
    """The standard figure of a run, each panel made with plot(): radiation power and
    field energy on log and linear scales, the energy of the whole pulse drawn over the
    window's where the file records it, the power inside the split angle where
    the file has it, bunching, the energy change and spread, the rms beam and field
    sizes, the beam's normalized and the field's emittances, a panel per polarization
    when the run carries two, a panel per harmonic when it carries harmonics, and the
    layout under them all. ``tex`` is plot()'s. Returns the figure with ``return_figure``."""
    import matplotlib.pyplot as plt
    st, close = _open(stats)
    try:
        with _style():
            fig = _summary(plt, st, tex)
    finally:
        if close:
            st.close()
    return fig if return_figure else None


def _summary(plt, st, tex):
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
    # The whole pulse, the window's energy and the energy slippage carried out of it,
    # drawn over the window's own, so light leaving the window does not read as light lost.
    energy = ["field/total/energy"]
    if "field/total/pulse_energy" in st:
        energy.insert(0, "field/total/pulse_energy")
    panels = {
        "power": dict(y=power, yscale="log"),
        "energy": dict(y=energy, yscale="log"),
        "power_lin": dict(y=power),
        "energy_lin": dict(y=energy),
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
        plot(st, ax=axd[key], include_layout=False, tex=tex, **kw)
        axd[key].set_xlabel("")
    xfactor, xprefix = _scale("m", [st.s])
    for key in ("layout_left", "layout_right"):
        layout(st, axd[key], xfactor=xfactor, xlabel=False, ylabels=False)
        axd[key].set_xlabel(f"$s$ ({xprefix}m)" if tex else f"s ({xprefix}m)")
    nslice = st["field/total/power"].shape[1]
    title = st.path.name + (f"   ({nslice} slices; thin lines are individual slices)" if nslice > 1 else "")
    fig.suptitle(title, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


# The compare figure's two curves: the reference thick and solid, the run thinner and
# dashed over it, so the two read apart in gray as well as in color.
REF_LINE = dict(color=COLORS[0], lw=2.8, ls="-")
RUN_LINE = dict(color=COLORS[1], lw=1.3, ls=(0, (4, 2)))
# The symbol and the plain name of each curve the compare figure draws, and of each
# difference, with the unit the top panel's values are in before any prefix.
COMPARE_LABEL = {
    "power": ("P", "power", "W"), "pulse_energy": (r"E_\mathrm{pulse}", "pulse energy", "J"),
    "bunching": ("|b|", "bunching", "1"),
    "energy": (r"\langle \Delta E \rangle", "mean energy change", "eV"),
    "exit_power": (r"P_\mathrm{exit}", "power at the exit", "W"),
    "size": (r"\sigma", "rms beam size", "m"), "field": ("|E|", "field amplitude", "V/m"),
}
COMPARE_DIFF = {
    "power": (r"\delta P / P", "dP/P"),
    "bunching": (r"\delta |b| \,/\, |b|_\mathrm{max}", "d|b| / max|b|"),
    "energy": (r"\delta \langle \Delta E \rangle \,/\, |\langle \Delta E \rangle|_\mathrm{max}",
               "d<dE> / max|<dE>|"),
    "exit_power": (r"\delta P_\mathrm{exit} \,/\, P_\mathrm{exit,max}", "dP_exit / max P_exit"),
    "size": (r"\delta \sigma / \sigma", "d sigma / sigma"),
    "field": (r"\delta |E| \,/\, |E|_\mathrm{max}", "d|E| / max|E|"),
}


def _power_of_ten(values, label):
    """The values scaled by a power of ten and the label carrying it, as "(10^-6)", so no
    axis shows an offset or a scale factor of its own. A scale between 1e-1 and 1e2 is
    left as it is."""
    m = float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 0.0
    e = int(np.floor(np.log10(m))) if m > 0 else 0
    if -1 <= e <= 2:
        return values, label
    return values / 10.0 ** e, f"{label} ($10^{{{e}}}$)"


def _compare_label(key, prefix, tex):
    sym, name, unit = COMPARE_LABEL[key]
    base = f"${sym}$" if tex else name
    return base + _unit_label(prefix, unit)


def _compare_scale(key, arrays):
    """(factor, prefix) for a compare curve. A quotient unit takes its prefix on its
    numerator, so a field reads GV/m rather than carrying an offset on the axis."""
    unit = COMPARE_LABEL[key][2]
    return _scale(unit.split("/")[0], arrays)


def _diff_label(key, tex):
    return f"${COMPARE_DIFF[key][0]}$" if tex else COMPARE_DIFF[key][1]


def compare(stats, reference, fields=None, labels=None, return_figure=False, tex=True):
    """Draw a run over a reference run with the difference beneath each curve.

    ``stats`` is the run, a Stats, its path or a ``lucifer.compare.Curves``, and
    ``reference`` the run it is drawn over: a Genesis4 ``.out.h5``, a second stats file,
    or Curves. Four blocks, lettered: (a) the power on a linear and a log axis, or the
    pulse energy for a time-dependent run, (b) the bunching, (c) the mean energy change
    from the first record, and (d) in steady state the field amplitude along the row
    through the center of the exit plane when ``fields`` names the run's and the
    reference's field dumps there, else the rms beam sizes, and in a time-dependent run
    the power across the window at the exit. The reference is drawn thick and solid, the
    run thinner and dashed over it. Under each block is Lucifer's value less the
    reference's, divided as ``lucifer.compare.DIVISOR`` says, and the axis label carries
    the power of ten the values are in. A time-dependent run's difference at a record is
    the one over its slices of largest magnitude, with its sign.

    A difference is drawn only where the two runs share their record positions and their
    window, as ``lucifer.compare.aligned`` tests, and the field row only where the two
    fields share a grid. Otherwise the block holds the overlay alone and the reason is
    printed. Nothing interpolates. ``labels`` names the run and the reference in the
    legend, in that order. Without it a run drawn over Genesis4 is Lucifer and a stats
    file is named by its stem. ``tex`` is plot()'s. Returns the figure with
    ``return_figure``.
    """
    import matplotlib.pyplot as plt
    from lucifer import compare as cmp

    run, ref = cmp.read_curves(stats), cmp.read_curves(reference)
    if labels is not None:
        run.name, ref.name = labels
    elif ref.kind == "genesis4" and run.kind == "lucifer":
        run.name = "Lucifer"
    with _style():
        fig = _compare(plt, cmp, run, ref, fields, tex)
    return fig if return_figure else None


def _compare(plt, cmp, run, ref, fields, tex):
    reason = cmp.aligned(run, ref)
    diffs = {}
    if reason:
        print(f"compare: the differences are not drawn: {reason}")
    else:
        diffs = cmp.differences(run, ref)
    td = run.n_slice > 1
    fa = fb = None
    if fields is not None and not td:
        fa, fb = cmp._fields(fields, run.n_slice)
        if not reason:
            try:
                diffs["field"] = cmp.differences(run, ref, fields=(fa, fb))["field"]
            except ValueError as err:
                print(f"compare: the field difference is not drawn: {err}")
    sfactor, sprefix = _scale("m", [run.s, ref.s])
    xs = f"$s$ ({sprefix}m)" if tex else f"s ({sprefix}m)"

    fig = plt.figure(figsize=(11.5, 11), layout="constrained")
    blocks = fig.subfigures(2, 2, hspace=0.04, wspace=0.03)

    def overlay(ax, xr, yr, xv, yv, label, log=False):
        ax.plot(xr, yr, label=ref.name, **REF_LINE)
        ax.plot(xv, yv, label=run.name, **RUN_LINE)
        ax.set_ylabel(label)
        if log:
            ax.set_yscale("log")
            # A dark start's first record is exactly zero. Do not let it stretch the axis
            # below the light. A curve within one decade keeps its own range.
            pos = yr[yr > 0]
            if pos.size and pos.max() / pos.min() > 10:
                ax.set_ylim(bottom=max(pos.min() * 0.3, pos.max() * 1e-14), top=pos.max() * 3)
            _log_ticks(ax)
        ax.tick_params(labelbottom=False)

    def difference(ax, x, d, key, xlabel):
        v, lab = _power_of_ten(d, _diff_label(key, tex))
        ax.axhline(0, color=INK2, lw=0.8)
        ax.plot(x, v, color=INK, lw=1.2)
        ax.set_ylabel(lab)
        ax.set_xlabel(xlabel)

    def letter(ax, tag):
        ax.text(0.0, 1.03, f"({tag})", transform=ax.transAxes, ha="left", va="bottom",
                fontweight="bold", color=INK)

    def rows(block, n_top, with_diff):
        ratios = [1.2] * n_top + ([1] if with_diff else [])
        ax = block.subplots(len(ratios), 1, sharex=True, gridspec_kw={"height_ratios": ratios})
        return list(np.atleast_1d(ax))

    # (a) The power, or the pulse energy, on a linear and a log axis.
    if td:
        key, yr, yv = "pulse_energy", ref.pulse_energy, run.pulse_energy
    else:
        key, yr, yv = "power", ref.power[:, 0], run.power[:, 0]
    factor, prefix = _compare_scale(key, [yr, yv])
    ax = rows(blocks[0, 0], 2, "power" in diffs)
    overlay(ax[0], ref.s / sfactor, yr / factor, run.s / sfactor, yv / factor, _compare_label(key, prefix, tex))
    ax[0].set_ylim(bottom=0)
    ax[0].legend(loc="upper left")
    letter(ax[0], "a")
    overlay(ax[1], ref.s / sfactor, yr / factor, run.s / sfactor, yv / factor,
            _compare_label(key, prefix, tex), log=True)
    if "power" in diffs:
        difference(ax[2], diffs["power"].x / sfactor, diffs["power"].values, "power", xs)
    else:
        ax[1].tick_params(labelbottom=True)
        ax[1].set_xlabel(xs)

    # (b) The bunching and (c) the mean energy change, slice means in a time-dependent run.
    for block, tag, key, yr, yv in (
            (blocks[0, 1], "b", "bunching", ref.bunching.mean(axis=1), run.bunching.mean(axis=1)),
            (blocks[1, 0], "c", "energy", ref.energy_change.mean(axis=1), run.energy_change.mean(axis=1))):
        factor, prefix = _compare_scale(key, [yr, yv])
        ax = rows(block, 1, key in diffs)
        overlay(ax[0], ref.s / sfactor, yr / factor, run.s / sfactor, yv / factor, _compare_label(key, prefix, tex))
        if key == "bunching":
            ax[0].set_ylim(bottom=0)
        letter(ax[0], tag)
        if key in diffs:
            difference(ax[1], diffs[key].x / sfactor, diffs[key].values, key, xs)
        else:
            ax[0].tick_params(labelbottom=True)
            ax[0].set_xlabel(xs)

    # (d) The power across the window at the exit, the field row at the exit, or the sizes.
    block = blocks[1, 1]
    if td:
        xr, xv = np.arange(ref.n_slice) * ref.ds, np.arange(run.n_slice) * run.ds
        xf, xp = _scale("m", [xr, xv])
        yr, yv = ref.power[-1], run.power[-1]
        factor, prefix = _compare_scale("exit_power", [yr, yv])
        ax = rows(block, 1, "exit_power" in diffs)
        ax[0].step(xr / xf, yr / factor, where="mid", label=ref.name, **REF_LINE)
        ax[0].step(xv / xf, yv / factor, where="mid", label=run.name, **RUN_LINE)
        ax[0].set_ylabel(_compare_label("exit_power", prefix, tex))
        ax[0].set_ylim(bottom=0)
        ax[0].tick_params(labelbottom=False)
        letter(ax[0], "d")
        xlabel = (f"$s$ along the window ({xp}m)" if tex else f"s along the window ({xp}m)")
        if "exit_power" in diffs:
            d = diffs["exit_power"]
            v, lab = _power_of_ten(d.values, _diff_label("exit_power", tex))
            ax[1].axhline(0, color=INK2, lw=0.8)
            ax[1].step(d.x / xf, v, where="mid", color=INK, lw=1.2)
            ax[1].set_ylabel(lab)
            ax[1].set_xlabel(xlabel)
        else:
            ax[0].tick_params(labelbottom=True)
            ax[0].set_xlabel(xlabel)
    elif fa is not None:
        ua, ub = fa["u"][0], fb["u"][0]
        ca, cb = ua.shape[0] // 2, ub.shape[0] // 2
        xa = (np.arange(ua.shape[1]) - ua.shape[1] // 2) * fa["dx"]
        xb = (np.arange(ub.shape[1]) - ub.shape[1] // 2) * fb["dx"]
        xf, xp = _scale("m", [xa, xb])
        factor, prefix = _compare_scale("field", [np.abs(ub[cb]), np.abs(ua[ca])])
        ax = rows(block, 1, "field" in diffs)
        ax[0].plot(xb / xf, np.abs(ub[cb]) / factor, label=ref.name, **REF_LINE)
        ax[0].plot(xa / xf, np.abs(ua[ca]) / factor, label=run.name, **RUN_LINE)
        ax[0].set_ylabel(_compare_label("field", prefix, tex))
        ax[0].set_ylim(bottom=0)
        ax[0].tick_params(labelbottom=False)
        letter(ax[0], "d")
        xlabel = f"$x$ ({xp}m)" if tex else f"x ({xp}m)"
        if "field" in diffs:
            difference(ax[1], diffs["field"].x / xf, diffs["field"].values, "field", xlabel)
        else:
            ax[0].tick_params(labelbottom=True)
            ax[0].set_xlabel(xlabel)
    else:
        factor, prefix = _compare_scale("size", [ref.size_x, ref.size_y, run.size_x, run.size_y])
        with_diff = "size_x" in diffs
        ax = rows(block, 1, with_diff)
        for comp, yr, yv, shade in (("x", ref.size_x[:, 0], run.size_x[:, 0], 1.0),
                                     ("y", ref.size_y[:, 0], run.size_y[:, 0], 0.55)):
            ax[0].plot(ref.s / sfactor, yr / factor, label=f"{ref.name} {comp}", alpha=shade, **REF_LINE)
            ax[0].plot(run.s / sfactor, yv / factor, label=f"{run.name} {comp}", alpha=shade, **RUN_LINE)
        ax[0].set_ylabel(_compare_label("size", prefix, tex))
        ax[0].set_ylim(bottom=0)
        ax[0].legend(loc="lower center", ncol=2)
        ax[0].tick_params(labelbottom=False)
        letter(ax[0], "d")
        if with_diff:
            # One power of ten for both planes, so the two lines share an axis.
            both = np.concatenate([diffs["size_x"].values, diffs["size_y"].values])
            v, lab = _power_of_ten(both, _diff_label("size", tex))
            e = len(diffs["size_x"].values)
            ax[1].axhline(0, color=INK2, lw=0.8)
            ax[1].plot(diffs["size_x"].x / sfactor, v[:e], color=INK, lw=1.2, label="x")
            ax[1].plot(diffs["size_y"].x / sfactor, v[e:], color=INK2, lw=1.2, ls=(0, (4, 2)), label="y")
            ax[1].set_ylabel(lab)
            ax[1].set_xlabel(xs)
            ax[1].legend(loc="lower right", ncol=2)
        else:
            ax[0].tick_params(labelbottom=True)
            ax[0].set_xlabel(xs)
    return fig


def main(argv=None) -> int:
    """python -m lucifer.plot <stats.h5> [-o out.png] [--dpi N] [--show]: write the
    summary figure of a run, by default to <out_root>.png in the current directory.
    python -m lucifer.plot compare <stats.h5> <reference> [--fields <run.wf.h5>
    <reference field>] [-o out.png]: write compare()'s figure of the run over the
    reference, by default to <out_root>-compare.png."""
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "compare":
        ap = argparse.ArgumentParser(prog="python -m lucifer.plot compare",
                                     description="write the figure of a run over a reference run, "
                                                 "Genesis4's .out.h5 or a second stats file")
        ap.add_argument("stats", help="the run, <out_root>.stats.h5")
        ap.add_argument("reference", help="the reference, a Genesis4 .out.h5 or a stats file")
        ap.add_argument("--fields", nargs=2, metavar=("RUN_FIELD", "REFERENCE_FIELD"),
                        help="the two runs' field dumps at the exit, for the field row of block (d)")
        ap.add_argument("--labels", nargs=2, metavar=("RUN", "REFERENCE"),
                        help="the legend's names for the run and the reference")
        ap.add_argument("-o", "--output", help="the image to write (default: <out_root>-compare.png)")
        ap.add_argument("--dpi", type=int, default=160, help="output resolution (default: 160)")
        ap.add_argument("--show", action="store_true", help="also open a window")
        args = ap.parse_args(argv[1:])
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        out = args.output or pathlib.Path(args.stats).name.removesuffix(".stats.h5") + "-compare.png"
        fig = compare(args.stats, args.reference, fields=args.fields, labels=args.labels, return_figure=True)
    else:
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
