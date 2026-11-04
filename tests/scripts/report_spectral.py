#!/usr/bin/env python3
"""
Write the tables, the figures and the diagrams of docs/spectral-field.md from the record
spectral_field.py wrote.

  report_spectral.py --json docs/generated/spectral/measured.json
                     --out docs/generated/spectral --page docs/generated/spectral-measured.md
                     [--included-by docs/spectral-field.md]

The record is one JSON file, every number the page states and the curves its figures
draw, with the provenance of the run that measured them. This script turns it into
docs/generated/spectral-measured.md, a page of blocks the spectral-field page includes
between markers, each a table directive or a figure directive with its label and its
caption, and into the figures beside it, docs/generated/spectral/*.svg, through
lucifer.plot.save. The three diagrams of the page are drawn here too, from nothing but
the propagator's constants, so they redraw with the layout.

Each block in the output is

  <!-- begin <label> -->
  ...
  <!-- end <label> -->

and the page includes it with the include directive's start-after and end-before. The
figure paths are written relative to the directory of --included-by, since MyST resolves
an included page's paths against the page that includes it.

The output is deterministic: the same record gives the same bytes, and the keystone
regenerates the page and requires no diff. The figures are SVG written through
lucifer.plot.save, which fixes the layout, salts the element ids and leaves out the
date, so the same arrays give the same bytes, and they stay out of the keystone's diff
as the example figures do, since a plotting library upgrade rewrites them all with no
physics moving. A caption carries the numbers its figure draws, so a moved curve fails
the keystone through its caption.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

# The pad beyond the window's head, as fel_spectral_mod sizes it, for the diagrams.
N_FREE, N_SPONGE = 4, 12
# The slice propagator's rotation threshold, fel_apply_slippage's, for the diagram.
ROTATE_AT = 0.8
# The derived sampling and the cliff, as fractions of the cooperation length
# (fel_spectral_mod), drawn on the sweep's figure.
COOP_FRACTION, COOP_CLIFF = 0.1, 0.25

PROPAGATOR = {"slices": "slice propagator", "spectral": "spectral propagator"}
MACHINE_LABELS = {"aramis": "Aramis, 1280 slices", "sase": "`sase` example", "flash1": "`flash1` example",
                  "lcls": "`lcls` example"}


# --------------------------------------------------------------------------------------
# Numbers as the page writes them.

def sci(x, digits=2):
    """x in scientific notation with digits significant figures and a plain exponent,
    1.8e-5 rather than 1.8e-05, and 0 for zero."""
    if x == 0:
        return "0"
    if not math.isfinite(x):
        return "nan"
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m}e{int(e)}"


def sig(x, n=3):
    """x to n significant figures, written positionally: 0.37, 1.7, 117, 0.047."""
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    r = round(x, n - 1 - e)
    if e >= n - 1:
        return f"{int(round(r)):d}"
    return f"{r:.{n - 1 - e}f}"


PREFIXES = [(1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n"),
            (1e-12, "p"), (1e-15, "f")]


def si(x, unit, n=3):
    """x with an SI prefix so the mantissa sits in [1, 1000): 625 MW, 6.59 kW, 974 nJ."""
    if x == 0:
        return f"0 {unit}"
    for scale, p in PREFIXES:
        if abs(x) >= scale:
            return f"{sig(x / scale, n)} {p}{unit}"
    scale, p = PREFIXES[-1]
    return f"{sig(x / scale, n)} {p}{unit}"


def pct(x, n=2):
    """A relative change as a percentage to n significant figures: 0.5 percent, 31 percent."""
    return f"{sig(100 * x, n)} percent"


# --------------------------------------------------------------------------------------
# The blocks of the page.

def block(label, body):
    return [f"<!-- begin {label} -->", *body, f"<!-- end {label} -->", ""]


def table(label, caption, header, rows):
    body = [f"```{{table}} {caption}", f":name: {label}", "",
            "| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    body += ["| " + " | ".join(r) + " |" for r in rows]
    body.append("```")
    return block(label, body)


def figure(label, path, caption):
    return block(label, [f":::{{figure}} {path}", f":label: {label}", ":class: col-body-outset-right", "", caption, ":::"])


# --------------------------------------------------------------------------------------
# The tables.

def table_provenance(rec):
    p = rec["provenance"]
    rows = [[p["lucifer_version"], p["bmad_version"], p["toolchain"], p["machine"], p["day"],
             str(p["threads"]), ", ".join(str(t) for t in p["cost_threads"]),
             "yes" if p.get("exe_open") else "no"]]
    return table("tab-sp-provenance",
                 "The run that measured every table and figure below, as the all subcommand of "
                 "spectral_field.py recorded it. The thread count is every run's but the cost's, which runs at each count of its own "
                 "column. The sweep's rows above the cliff run on a second tracker built with the cliff "
                 "raised when one was given.",
                 ["Lucifer", "Bmad", "toolchain", "machine", "day", "threads", "cost threads", "cliff-open tracker"],
                 rows)


def table_vacuum(rec):
    v = rec["vacuum"]
    rows = [[sci(r["bandwidth"], 2) + (" (no chirp)" if r["chirp"] == 0 else ""), sci(r["chirp"], 3),
             sci(r["slices"], 3), sci(r["spectral"], 3)] for r in v["rows"]]
    return table("tab-sp-vacuum",
                 f"The vacuum test: a chirped Gaussian pulse of rms length {sig(v['sigma_slices'], 2)} slices and waist "
                 f"{si(v['waist'], 'm')} at {si(v['wavelength'], 'm')}, on {v['n_slice']} slices of one wavelength "
                 f"and a grid of {v['ngrid']} points of {si(v['dx'], 'm')}, drifted {sig(v['length'], 4)} m, two "
                 f"Rayleigh ranges, under each propagator. The error is the largest difference from the exact "
                 f"paraxial solution over the window, over the solution's largest amplitude.",
                 ["Relative rms bandwidth", "Chirp $b$ [1/m²]", "Slice propagator", "Spectral propagator"], rows)


# The differences against the slice run the page states, by lucifer.compare's key: the
# per-record power is left out, since on the SASE tier it is taken at records where the
# slice run holds no light yet.
OWN_KEYS = (("exit_power", "exit power"), ("power", "power"), ("bunching", "bunching"), ("energy", "energy"))


def own_text(own, td):
    keys = OWN_KEYS[0:1] + OWN_KEYS[2:] if td else OWN_KEYS[1:]
    return ", ".join(f"{name} {sci(own[k], 2)}" for k, name in keys if k in own)


def table_null(rec):
    n = rec["null"]
    rows = []
    for t in n["tiers"]:
        own = t["against_slices"]
        rows.append([f"`{t['tier']}`", "slices, the recorded level", sci(t["recorded"], 5), sci(t["tolerance"], 1), ""])
        rows.append([f"`{t['tier']}`", "spectral", sci(t["spectral"], 5), sci(t["tolerance"], 1),
                     own_text(own, t["n_slice"] > 1)])
    f = n["filter_on"]["against_slices"]
    rows.append(["`tdsase`, source filter on", "spectral against slices", "", "", own_text(f, True)])
    return table("tab-sp-null",
                 f"The null test at {n['threads']} threads: the two tiers of the validation suite under the "
                 f"spectral propagator, against Genesis4 at the tier tolerances as the validation suite compares "
                 f"them, and against the validation pass's own slice run through lucifer.compare, each "
                 f"difference the largest over the records as a fraction of the reference. The last row is "
                 f"the SASE tier again under both propagators with the source filter on, the production "
                 f"default, Lucifer against Lucifer.",
                 ["Tier", "Propagator", "Against Genesis4", "Tolerance", "Against the slice run"], rows)


def table_bank(rec):
    b = rec["bank"]
    rows = []
    for key in ("slices", "spectral"):
        r = b["rows"][key]
        esc = r.get("escaped_energy", 0.0)
        rows.append([key, si(r["window_energy"], "J", 4), si(esc, "J", 4), si(r["window_energy"] + esc, "J", 4),
                     f"{r['pulse_planes']} planes, {si(r['pulse_energy'], 'J', 4)}", sci(r["beam_loss"], 2),
                     sci(r["head_share"], 2), sci(r["tail_share"], 2)])
    return table("tab-sp-bank",
                 f"The bank at {b['threads']} threads: the SASE tier's 32-slice window under a line that slips "
                 f"1184 slices, with the escaped field kept and the source filter off, under both propagators. "
                 f"The window's energy at the exit, the escaped energy, their sum, what the pulse file holds, "
                 f"the beam's energy change over the walk, and the head and tail slices' shares of the exit "
                 f"power.",
                 ["Propagator", "Window energy", "Escaped energy", "Sum", "Pulse file", "Beam energy change [J]",
                  "Head slice", "Tail slice"], rows)


def table_cost(rec):
    c = rec["cost"]
    rows = []
    for r in c["rows"]:
        rows.append([str(r["threads"]), r["propagator"], f"{r['solve_seconds']:.1f} s",
                     f"{1e3 * r['solve_seconds'] / r['n_step']:.1f} ms", f"{r['walk_seconds']:.1f} s"])
    f = c["fit"]
    r0 = c["rows"][0]
    return table("tab-sp-cost",
                 f"The cost: the field solve row of the footer's timing table on the SASE tier, {r0['n_field']} "
                 f"planes and a pad of {r0['n_pad'] if r0['propagator'] == 'spectral' else f['n_pad']} on a "
                 f"{r0['ngrid']} by {r0['ngrid']} grid over {r0['n_step']} steps, one run each on an otherwise "
                 f"idle machine. Under the spectral propagator the row holds the deposit, the kernel pass over "
                 f"the window and its pad, and the slip passes. Under the slice propagator it holds the deposit "
                 f"and its transform pair. From the {f['threads']}-thread rows, the excess of {f['excess_seconds']:.1f} s "
                 f"less the pad's transverse share of {f['pad_seconds']:.1f} s gives the rate of the passes along $\\tau$, "
                 f"{sci(f['tau'], 3)} s per point per log2 of the transform's length per pass.",
                 ["Threads", "Propagator", "Field solve", "A step", "Walk"], rows)


def table_seeded(rec):
    s = rec["seeded"]
    names = {"quiet": "seed alone, segment", "segment": "seed and noise, segment", "line": "seed and noise, line"}
    rows = []
    for r in s["rows"]:
        rows.append([f"{100 * r['bandwidth']:.2g} %", names[r["extent"]], f"{r['offset']:.3f}", sci(r["field_raw"], 2),
                     sci(r["field"], 2), sci(r["wavelength"], 2), f"{si(r['power_a'], 'W')} / {si(r['power_b'], 'W')}"])
    return table("tab-sp-seeded",
                 f"A chirped seed through Aramis at {s['threads']} threads: the vacuum test's pulse at the tier's "
                 f"wavelength, rms length {s['sigma_wavelengths']} wavelengths, waist {si(s['waist'], 'm')}, "
                 f"{si(s['peak_power'], 'W')} at its peak, {s['seed_slice']} slices from the tail of a window of "
                 f"{s['n_slice']} slices, with the tier's beam, through the first segment and through the line, "
                 f"with the source filter on, and the segment again on the bunch loaded without its shot noise. "
                 f"The offset is the slippage residual the slice run's field carries, in slices, which the "
                 f"aligned column removes before differencing over the lit slices. The field difference is "
                 f"over the spectral run's peak, the wavelength difference the largest over the seed's core, "
                 f"and the peak powers are the slice run's and the spectral run's.",
                 ["Bandwidth", "Light", "Offset [slices]", "Field, raw", "Field, aligned", "Wavelength",
                  "Peak power, slices / spectral"], rows)


def table_convergence(rec):
    s = rec["seeded"]
    rows = []
    for c in s["convergence"]:
        a, b = c["fine_against_spectral"], c["slices_against_fine"]
        rows.append([c["extent"], sci(a["field"], 2), sci(a["wavelength"], 2),
                     f"{si(a['power_a'], 'W')} / {si(a['power_b'], 'W')}", sci(b["field"], 2)])
    return table("tab-sp-convergence",
                 "The reference's convergence: the 2 percent seed under the spectral propagator with the "
                 "undulator step cut from three periods to one, against the three-period run, through the "
                 "segment and through the line. The last column is the slice run against the fine spectral "
                 "run, aligned.",
                 ["Extent", "Field", "Wavelength", "Peak power, fine / three periods", "Slice run against the fine run"],
                 rows)


def table_chirped(rec):
    c = rec["chirped"]
    rows = []
    for r in c["rows"]:
        rows.append([f"{100 * r['chirp']:.2g} %" if r["chirp"] else "0", sci(r["field"], 2), sci(r["wavelength"], 2),
                     str(r["lit"]), f"{si(r['power_a'], 'W')} / {si(r['power_b'], 'W')}",
                     f"{si(r['pulse_energy_slices'], 'J')} / {si(r['pulse_energy_spectral'], 'J')}",
                     f"{r['slope_slices']:+.4f} / {r['slope_spectral']:+.4f}",
                     f"{r['expected_slope']:+.4f}" if r["chirp"] else "0"])
    return table("tab-sp-chirped",
                 f"An energy-chirped beam from noise at {c['threads']} threads: the same beam with a linear energy "
                 f"chirp along the window, head to tail, the head above the tail, through the full line under "
                 f"both propagators with the source filter on. The slope is that of a line fitted to the exit "
                 f"wavelengths over the lit slices, as a fraction of the mean over the window, against the "
                 f"$-2\\epsilon$ a chirp $\\epsilon$ in energy puts on the resonant wavelength.",
                 ["Chirp", "Field", "Wavelength", "Lit slices", "Peak power, slices / spectral",
                  "Pulse energy, slices / spectral", "Slope, slices / spectral", "$-2\\epsilon$"], rows)


def table_sweep(rec):
    sw = rec["sweep"]
    rows = []
    for key in ("aramis", "sase", "flash1", "lcls"):
        m = sw["machines"].get(key)
        if m is None:
            continue
        first = True
        for r in m["rows"]:
            label = MACHINE_LABELS[key] if first else ""
            first = False
            if "refused" in r:
                rows.append([label, str(r["factor"]), str(r["n_wavelength"]), sig(r["fraction"], 2),
                             "refused: fewer than 8 planes" if "FEWER THAN" in r["refused"] else "refused: above the cliff",
                             "", "", "", "", ""])
                continue
            rows.append([label, str(r["factor"]) + (" (cliff open)" if r["cliff_open"] else ""), str(r["n_wavelength"]),
                         sig(r["fraction"], 2), str(r["planes"]), si(r["window_energy"], "J"), si(r["pulse_energy"], "J"),
                         si(r["peak_power"], "W"), sci(r.get("exit_diff", 0.0), 2) if r["factor"] > 1 else "0",
                         f"{r['walk']:.0f} s"])
    lines = []
    for key in ("aramis", "sase", "flash1", "lcls"):
        m = sw["machines"].get(key)
        if m:
            lines.append(f"{MACHINE_LABELS[key].replace('`', '')}: {sig(m['l_coop'] / m['wavelength'], 3)} wavelengths at "
                         f"$\\rho$ = {sci(m['rho'], 2)}")
    return table("tab-sp-sweep",
                 f"The sweep at {sw['threads']} threads: four lines from noise under the spectral propagator, "
                 f"each at {', '.join(str(f) for f in sw['factors'])} slices a plane, against its own run at "
                 f"one plane a slice. The window's and the pulse's energy at the exit, the peak slice power, "
                 f"the exit power per beam slice against the finest run over its peak, and the walk. The "
                 f"plane's length is also given as the fraction of the cooperation length "
                 f"$\\lambda_s / 4\\pi\\rho$ at the beam's $\\rho$, which is {'; '.join(lines)}. A row marked cliff open "
                 f"ran on the tracker built with the cliff raised.",
                 ["Line", "Slices a plane", "Plane [wavelengths]", "Of $L_c$", "Planes", "Window energy",
                  "Pulse energy", "Peak power", "Exit power against the finest", "Walk"], rows)


# --------------------------------------------------------------------------------------
# The figures.

def style():
    import matplotlib as mpl
    from lucifer.plot import STYLE
    return mpl.rc_context(STYLE)


def colors():
    from lucifer.plot import COLORS
    return COLORS


def fig_vacuum(rec, plt):
    v = rec["vacuum"]
    bw = np.array([r["bandwidth"] for r in v["rows"]])
    sl = np.array([r["slices"] for r in v["rows"]])
    sp = np.array([r["spectral"] for r in v["rows"]])
    fig, ax = plt.subplots(figsize=(6.4, 3.6), layout="constrained")
    c = colors()
    ax.loglog(bw, sl, "o-", color=c[0], label="slice propagator")
    ax.loglog(bw, sp, "s-", color=c[1], label="spectral propagator")
    ax.loglog(bw, 0.77 * bw, ":", color="gray", label="0.77 of the bandwidth")
    from matplotlib.ticker import NullFormatter
    ax.set_xticks(bw)
    ax.set_xticklabels([f"{100 * b:.2g} %" for b in bw])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("relative rms bandwidth of the pulse")
    ax.set_ylabel("largest error over the exact solution's peak")
    ax.legend(frameon=False)
    ratio = sl / bw
    span = (f"{sig(ratio.min(), 2)} of the relative bandwidth at every row" if sig(ratio.min(), 2) == sig(ratio.max(), 2)
            else f"{sig(ratio.min(), 2)} to {sig(ratio.max(), 2)} of the relative bandwidth")
    cap = (f"The vacuum test: the error of each propagator against the exact paraxial solution, from "
           f"[](#tab-sp-vacuum). The slice propagator's error is {span}, the spatio-temporal coupling one "
           f"$k_0$ cannot carry over two Rayleigh ranges, and the spectral propagator's is {sci(sp.min(), 2)} "
           f"to {sci(sp.max(), 2)}, the rounding of its transform pairs.")
    return fig, cap


def fig_null(rec, plt):
    n = rec["null"]
    c = colors()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), layout="constrained")
    for ax, key, title in ((axes[0], "tdsase", "source filter off, as the tier runs"),
                           (axes[1], "tdsase_filter", "source filter on")):
        cv = n["curves"][key]
        a, b = np.array(cv["exit_slices"]), np.array(cv["exit_spectral"])
        i = np.arange(1, a.size + 1)
        scale, prefix = next((sc, pf) for sc, pf in PREFIXES if max(a.max(), b.max()) >= sc)
        ax.plot(i, a / scale, "-", color=c[0], lw=2, label="slice propagator")
        ax.plot(i, b / scale, "--", color=c[1], lw=1.5, label="spectral propagator")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("slice")
        ax.set_ylabel(f"power at the exit [{prefix}W]")
    axes[0].legend(frameon=False)
    sase = next(t for t in n["tiers"] if t["tier"] == "tdsase")
    off = sase["against_slices"]
    on = n["filter_on"]["against_slices"]
    cap = (f"The null test: the SASE tier's power across the window at the exit under both propagators, with "
           f"the source filter off (left) and on (right), from the runs of [](#tab-sp-null). With the filter "
           f"off the two differ by {sci(off['exit_power'], 2)} of the exit power, {sci(off['bunching'], 2)} of "
           f"the bunching and {sci(off['energy'], 2)} of the energy, the wide-angle start-up light the "
           f"frequency-dependent kernel delays. With it on they differ by {sci(on['exit_power'], 2)}, "
           f"{sci(on['bunching'], 2)} and {sci(on['energy'], 2)}.")
    return fig, cap


def fig_bank(rec, plt):
    b = rec["bank"]
    c = colors()
    fig, ax = plt.subplots(figsize=(6.4, 3.6), layout="constrained")
    for key, color, ls, lw in (("slices", c[0], "-", 2), ("spectral", c[1], "--", 1.5)):
        p = np.array(b["rows"][key]["pulse_power"])
        i = np.arange(1, p.size + 1)
        ax.semilogy(i, np.maximum(p, p.max() * 1e-12), ls, color=color, lw=lw, label=PROPAGATOR[key])
    ax.set_xlabel("plane of the pulse file, from its tail")
    ax.set_ylabel("power [W]")
    ax.legend(frameon=False)
    rs, rp = b["rows"]["slices"], b["rows"]["spectral"]
    cap = (f"The bank: the power in each plane of the pulse file the SASE tier's run leaves under both "
           f"propagators, from [](#tab-sp-bank). Both bank {rs['pulse_planes']} planes. The pulse holds the "
           f"window's energy plus the escaped energy to "
           f"{sci(abs(rs['pulse_energy'] / (rs['window_energy'] + rs.get('escaped_energy', 0)) - 1), 1)} under "
           f"the rotation and {sci(abs(rp['pulse_energy'] / (rp['window_energy'] + rp.get('escaped_energy', 0)) - 1), 1)} "
           f"under the ramp, and the spectral run's pulse carries "
           f"{pct(rp['pulse_energy'] / rs['pulse_energy'] - 1)} more energy than the slice run's. The tail "
           f"slice's share of the exit power is {sci(rp['tail_share'], 1)} under the ramp.")
    return fig, cap


def fig_cost(rec, plt):
    cst = rec["cost"]
    c = colors()
    threads = sorted({r["threads"] for r in cst["rows"]})
    fig, ax = plt.subplots(figsize=(6.4, 3.4), layout="constrained")
    width = 0.36
    for k, (prop, color) in enumerate((("slices", c[0]), ("spectral", c[1]))):
        vals = [next(1e3 * r["solve_seconds"] / r["n_step"] for r in cst["rows"]
                     if r["threads"] == t and r["propagator"] == prop) for t in threads]
        x = np.arange(len(threads)) + (k - 0.5) * width
        ax.bar(x, vals, width, color=color, label=PROPAGATOR[prop])
        for xi, v in zip(x, vals):
            ax.text(xi, v, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(np.arange(len(threads)))
    ax.set_xticklabels([f"{t} thread{'s' if t > 1 else ''}" for t in threads])
    ax.set_ylabel("field solve, ms a step")
    ax.legend(frameon=False)
    ratios = []
    for t in threads:
        a = next(r["solve_seconds"] for r in cst["rows"] if r["threads"] == t and r["propagator"] == "slices")
        b = next(r["solve_seconds"] for r in cst["rows"] if r["threads"] == t and r["propagator"] == "spectral")
        ratios.append(f"{sig(b / a, 2)} at {t}")
    cap = (f"The cost: the field phase of a step of the SASE tier under both propagators, from "
           f"[](#tab-sp-cost). The spectral propagator's field phase is {' and '.join(ratios)} "
           f"threads times the slice propagator's on this window of {cst['fit']['n_field']} planes with a pad "
           f"of {cst['fit']['n_pad']}.")
    return fig, cap


def _two_rows(curves, keys, kappa, plt, titles, band=0.03):
    c = colors()
    fig, axes = plt.subplots(2, len(keys), figsize=(2.6 * len(keys) + 0.8, 5.6), layout="constrained")
    for j, key in enumerate(keys):
        cv = curves[key]
        ps, pp = np.array(cv["power_slices"]), np.array(cv["power_spectral"])
        i = np.arange(1, ps.size + 1)
        top, bot = axes[0, j], axes[1, j]
        top.plot(i, ps / 1e9, "-", color=c[0], lw=1.6, label="slice propagator")
        top.plot(i, pp / 1e9, "--", color=c[1], lw=1.2, label="spectral propagator")
        top.set_title(titles[j], fontsize=10)
        top.set_xlabel("slice")
        ss, sp = np.array(cv["spectrum_slices"]), np.array(cv["spectrum_spectral"])
        k = np.array(kappa)
        sel = np.abs(k) <= band
        scale = max(ss[sel].max(), sp[sel].max())
        bot.plot(k[sel], ss[sel] / scale, "-", color=c[0], lw=1.6)
        bot.plot(k[sel], sp[sel] / scale, "--", color=c[1], lw=1.2)
        bot.set_xlabel(r"$\kappa / k_0$")
        if j == 0:
            top.set_ylabel("power at the exit [GW]")
            bot.set_ylabel("exit spectrum over its peak")
            top.legend(frameon=False, fontsize=8)
    return fig


def fig_seeded(rec, plt):
    s = rec["seeded"]
    keys = [k for k in s["curves"] if k.endswith("_line")]
    titles = [f"{100 * s['curves'][k]['bandwidth']:.2g} % bandwidth" for k in keys]
    fig = _two_rows(s["curves"], keys, s["kappa_over_k0"], plt, titles)
    line = [r for r in s["rows"] if r["extent"] == "line"]
    quiet = [r for r in s["rows"] if r["extent"] == "quiet"]
    seg = [r for r in s["rows"] if r["extent"] == "segment"]
    cap = (f"A chirped seed through the Aramis line under both propagators, with the source filter on, at "
           f"each bandwidth of [](#tab-sp-seeded): the power across the window at the exit (top) and the "
           f"exit spectrum about the carrier (bottom). Through the line the propagators differ by "
           f"{sci(min(r['field'] for r in line), 2)} to {sci(max(r['field'] for r in line), 2)} of the peak "
           f"field, {sci(min(r['wavelength'] for r in line), 2)} to {sci(max(r['wavelength'] for r in line), 2)} "
           f"in wavelength and {pct(min(abs(r['power_a'] / r['power_b'] - 1) for r in line), 1)} to "
           f"{pct(max(abs(r['power_a'] / r['power_b'] - 1) for r in line), 1)} in peak power. With the seed "
           f"alone over the first segment they differ by {sci(min(r['field'] for r in quiet), 2)} to "
           f"{sci(max(r['field'] for r in quiet), 2)} of the peak field after alignment, and with the beam's "
           f"shot noise on by {sci(min(r['field'] for r in seg), 2)} to {sci(max(r['field'] for r in seg), 2)}.")
    return fig, cap


def fig_chirped(rec, plt):
    c = rec["chirped"]
    keys = list(c["curves"])
    titles = [(f"{100 * c['curves'][k]['chirp']:.2g} % chirp" if c["curves"][k]["chirp"] else "no chirp") for k in keys]
    fig = _two_rows(c["curves"], keys, c["kappa_over_k0"], plt, titles)
    rows = c["rows"]
    unchirped = rows[0]
    cap = (f"An energy-chirped beam from noise through the Aramis line under both propagators, at each chirp "
           f"of [](#tab-sp-chirped): the power across the window at the exit (top) and the exit spectrum "
           f"about the carrier (bottom). The chirp decides the output: against the unchirped run's "
           f"{si(unchirped['pulse_energy_spectral'], 'J')} the pulse energy is "
           f"{', '.join(si(r['pulse_energy_spectral'], 'J') for r in rows[1:])} at "
           f"{', '.join(f'{100 * r[chr(99) + chr(104) + chr(105) + chr(114) + chr(112)]:.2g}' for r in rows[1:])} "
           f"percent. The two propagators agree on the fitted slope of the exit wavelengths to "
           f"{sci(max(abs(r['slope_slices'] - r['slope_spectral']) for r in rows), 1)}, on the peak power to "
           f"{pct(max(abs(r['power_a'] / r['power_b'] - 1) for r in rows), 1)} and on the slice wavelengths to "
           f"{sci(max(r['wavelength'] for r in rows), 2)}, and their exit fields differ by "
           f"{sci(min(r['field'] for r in rows), 2)} to {sci(max(r['field'] for r in rows), 2)} of the peak.")
    return fig, cap


def fig_sweep(rec, plt):
    sw = rec["sweep"]
    c = colors()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), layout="constrained")
    for k, key in enumerate(("aramis", "sase", "flash1", "lcls")):
        m = sw["machines"].get(key)
        if m is None:
            continue
        rows = [r for r in m["rows"] if "refused" not in r]
        ref = rows[0]
        frac = np.array([r["fraction"] for r in rows])
        energy = np.array([r["pulse_energy"] for r in rows]) / ref["pulse_energy"]
        walk = np.array([r["walk"] for r in rows]) / ref["walk"]
        label = MACHINE_LABELS[key].replace("`", "")
        for ax, y in ((axes[0], energy), (axes[1], walk)):
            ax.plot(frac, y, "-", color=c[k], lw=1.4, label=label)
            closed = np.array([not r["cliff_open"] for r in rows])
            ax.plot(frac[closed], y[closed], "o", color=c[k], ms=5)
            ax.plot(frac[~closed], y[~closed], "o", mfc="white", mec=c[k], ms=5)
    for ax in axes:
        ax.set_xscale("log")
        ax.axvline(COOP_FRACTION, color="gray", ls=":", lw=1)
        ax.axvline(COOP_CLIFF, color="gray", ls="--", lw=1)
        ax.set_xlabel("plane length as a fraction of the cooperation length")
    axes[0].set_ylabel("pulse energy over the finest sampling's")
    axes[1].set_ylabel("walk over the finest sampling's")
    axes[1].set_yscale("log")
    axes[0].legend(frameon=False, fontsize=9)
    a = sw["machines"]["aramis"]
    rows = [r for r in a["rows"] if "refused" not in r]
    ref = rows[0]
    parts = [f"{pct(abs(r['pulse_energy'] / ref['pulse_energy'] - 1), 2)} at {sig(r['fraction'], 2)}" for r in rows[1:]]
    walks = [f"{r['walk']:.0f}" for r in rows]
    cap = (f"The sweep: the pulse energy at the exit (left) and the walk (right) against the field's plane "
           f"as a fraction of the cooperation length, each over the same line's run at one plane a slice, "
           f"from [](#tab-sp-sweep). The dotted line is the derived sampling, {COOP_FRACTION} of the "
           f"cooperation length, and the dashed line the cliff at {COOP_CLIFF}, above which a stated "
           f"sampling is refused. A hollow point ran with the cliff raised. On the Aramis line of "
           f"{ref['planes']} slices the pulse energy moves by {', '.join(parts)} of the cooperation length, "
           f"and the walk is {', '.join(walks)} s over the samplings.")
    return fig, cap


# --------------------------------------------------------------------------------------
# The diagrams.

def _box_row(ax, x0, n, y, h, w, face, edge, text=None, hatch=None):
    import matplotlib.patches as mp
    for i in range(n):
        ax.add_patch(mp.Rectangle((x0 + i * w, y), w, h, facecolor=face, edgecolor=edge, lw=0.8, hatch=hatch))
    if text:
        ax.text(x0 + n * w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)


def diagram_pad(plt):
    """The window, the pad's free slots and the sponge along tau, the plane the bank reads
    and the light the ramp carries."""
    import matplotlib.patches as mp
    c = colors()
    n_win = 10
    w, h = 1.0, 1.0
    fig, ax = plt.subplots(figsize=(9.6, 3.0), layout="constrained")
    _box_row(ax, 0, n_win, 0, h, w, "#dbe7f7", c[0], "the window, slice 1 to $n$")
    _box_row(ax, n_win, N_FREE, 0, h, w, "#f7f0dc", c[5], f"{N_FREE} free")
    _box_row(ax, n_win + N_FREE, N_SPONGE, 0, h, w, "#f3e0dc", c[4], f"{N_SPONGE} sponge slots")
    # The sponge's cos^2 weight, from 1 at its first slot to 0 at the seam.
    x = np.linspace(n_win + N_FREE, n_win + N_FREE + N_SPONGE, 200)
    wgt = np.cos(np.pi / 2 * (x - (n_win + N_FREE)) / N_SPONGE)**2
    ax.plot(x, 1.15 + 0.8 * wgt, color=c[4], lw=1.6)
    ax.text(n_win + N_FREE + N_SPONGE / 2, 2.1, r"$\cos^2$ weight after each slip pass", ha="center", fontsize=9)
    # The ramp carries light toward the head.
    ax.annotate("", xy=(n_win + 0.9, -0.5), xytext=(n_win - 4, -0.5),
                arrowprops=dict(arrowstyle="->", color=c[1], lw=1.6))
    ax.text(n_win - 1.5, -0.95, r"the ramp $e^{-i\kappa\delta}$ moves the light toward the head", ha="center",
            fontsize=9, color=c[1])
    # The plane the bank reads.
    ax.add_patch(mp.Rectangle((n_win, 0), w, h, facecolor="none", edgecolor=c[1], lw=2.2))
    ax.annotate("read into the bank at each whole-slice crossing", xy=(n_win + 0.5, 1.02), xytext=(n_win + 0.5, 2.5),
                ha="center", fontsize=9, arrowprops=dict(arrowstyle="->", color=INK(), lw=1))
    # The periodic seam.
    xs = n_win + N_FREE + N_SPONGE
    ax.plot([xs, xs], [-0.3, 1.3], color="gray", ls="--", lw=1)
    ax.text(xs + 0.15, -0.25, "periodic seam:\nthe transform's end\nmeets slice 1", fontsize=8, va="top")
    ax.annotate("", xy=(0, -2.0), xytext=(xs, -2.0), arrowprops=dict(arrowstyle="->", color="gray", lw=1,
                                                                    connectionstyle="arc3,rad=0.12"))
    ax.text(xs / 2, -1.85, "what passes the seam would reappear at the tail", ha="center", fontsize=8, color="gray")
    ax.annotate("", xy=(xs + 0.3, 0.5), xytext=(-0.6, 0.5))
    ax.text(0, -0.45, r"$\tau \rightarrow$ (toward the head)", fontsize=9, ha="left")
    ax.set_xlim(-0.8, xs + 4.5)
    ax.set_ylim(-2.7, 2.9)
    ax.set_axis_off()
    cap = (f"The window and the pad along $\\tau$. Every transform runs over the window's $n$ slices, "
           f"{N_FREE} free slots and {N_SPONGE} sponge slots together. The ramp carries light toward the "
           f"head and into the free slots, the plane at the first of them is read into the bank each time "
           f"the accumulated slippage crosses a whole slice, and the $\\cos^2$ sponge falls to zero at the "
           f"transform's periodic seam, so nothing comes round to the tail.")
    return fig, cap


def INK():
    from lucifer.plot import INK as ink
    return ink


def diagram_slippage(plt):
    """The rotation's whole-plane slippage against the ramp's continuous one over a few
    steps, on a field of four slices a plane under a step that slips one slice."""
    c = colors()
    n_per = 4
    per_step = 1.0 / n_per
    n = 12
    fig, ax = plt.subplots(figsize=(6.4, 3.4), layout="constrained")
    steps = np.arange(0, n + 1)
    ramp = per_step * steps
    acc, pos, rot = 0.0, 0, [0]
    for k in range(1, n + 1):
        acc += per_step
        if acc >= ROTATE_AT - 1e-12:
            pos += 1
            acc -= 1
        rot.append(pos)
    rot = np.array(rot, dtype=float)
    ax.step(steps, ramp, where="post", color=c[1], lw=2, label="the ramp: the step's own slippage")
    ax.step(steps, rot, where="post", color=c[0], lw=2,
            label=f"the rotation: a whole plane once {ROTATE_AT} of one has accumulated")
    ax.fill_between(steps, rot, ramp, step="post", color=c[0], alpha=0.12)
    ax.set_xlabel("undulator step, each slipping one slice")
    ax.set_ylabel("field position against the beam [planes]")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    cap = (f"The slippage's two forms on a field of {n_per} slices a plane, under a step that slips one "
           f"slice, a quarter of a plane. The ramp moves the field by the step's own slippage. The slice "
           f"propagator's rotation moves the record by whole planes once {ROTATE_AT} of a plane has "
           f"accumulated, so its field would sit up to {ROTATE_AT} of a plane behind the ramp's (shaded), "
           f"which is why a field coarser than the beam's slices needs the ramp.")
    return fig, cap


def diagram_planes(plt):
    """The field's planes over the beam's slices, with the two-plane weights of the
    readback and the deposit."""
    import matplotlib.patches as mp
    c = colors()
    n_per = 4
    n_slice = 9
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.2), layout="constrained")
    for ax, title, down in ((axes[0], "readback: a slice reads its two planes", True),
                            (axes[1], "deposit: a slice's source goes to its two planes", False)):
        for i in range(n_slice):
            ax.add_patch(mp.Rectangle((i, 0), 1, 0.7, facecolor="#dbe7f7", edgecolor=c[0], lw=0.8))
            ax.text(i + 0.5, 0.35, str(i + 1), ha="center", va="center", fontsize=8)
        for j in range(0, n_slice, n_per):
            ax.plot([j, j], [1.4, 2.4], color=c[1], lw=3)
            ax.text(j, 2.55, f"plane {j // n_per + 1}", ha="center", fontsize=8, color=c[1])
        s = 2   # slice 3, between plane 1 at slice 1 and plane 2 at slice 5
        w = (s - 0) / n_per
        ax.add_patch(mp.Rectangle((s, 0), 1, 0.7, facecolor="#f7f0dc", edgecolor=c[5], lw=1.6))
        lo, hi = 0, n_per
        if down:
            ax.annotate("", xy=(s + 0.5, 0.75), xytext=(lo, 1.4), arrowprops=dict(arrowstyle="->", color=INK(), lw=1.2))
            ax.annotate("", xy=(s + 0.5, 0.75), xytext=(hi, 1.4), arrowprops=dict(arrowstyle="->", color=INK(), lw=1.2))
            ax.text(lo - 0.3, 1.5, f"$1 - w$ = {1 - w:.2g}", fontsize=9, ha="right", va="center")
            ax.text(hi + 0.3, 1.5, f"$w$ = {w:.2g}", fontsize=9, ha="left", va="center")
        else:
            ax.annotate("", xy=(lo, 1.4), xytext=(s + 0.5, 0.75), arrowprops=dict(arrowstyle="->", color=INK(), lw=1.2))
            ax.annotate("", xy=(hi, 1.4), xytext=(s + 0.5, 0.75), arrowprops=dict(arrowstyle="->", color=INK(), lw=1.2))
            ax.text(lo - 0.3, 1.5, f"$(1 - w) / m$ = {(1 - w) / n_per:.3g}", fontsize=9, ha="right", va="center")
            ax.text(hi + 0.3, 1.5, f"$w / m$ = {w / n_per:.3g}", fontsize=9, ha="left", va="center")
        ax.text(n_slice / 2, -0.35, f"beam slices, $m$ = {n_per} slices a plane; $w$ = (slice $-$ plane) / $m$",
                ha="center", fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-3.2, n_slice + 1.2)
        ax.set_ylim(-0.7, 2.9)
        ax.set_axis_off()
    cap = (f"A field coarser than the beam's slices, here $m$ = {n_per} slices a plane. A beam slice sits "
           f"between two planes and reads them by the linear weights $1 - w$ and $w$ along $\\tau$ (left), so "
           f"a slice on a plane reads that plane alone. The deposit splits the slice's source between the "
           f"same two planes by the same weights divided by the $m$ slices a plane spans (right), so a "
           f"uniform source over the slices gives every plane one slice's source, and the source filter "
           f"runs once on each plane's sum.")
    return fig, cap


# --------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    ap.add_argument("--json", required=True, help="the record spectral_field.py wrote")
    ap.add_argument("--out", required=True, help="where the figures are written, <out>/<name>.svg")
    ap.add_argument("--page", required=True, help="the Markdown page of blocks to write")
    ap.add_argument("--included-by", default="docs/spectral-field.md",
                    help="the page that includes --page, whose directory the figure paths are relative to")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lucifer.plot import save

    rec = json.loads(Path(args.json).read_text())
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    link_dir = Path(args.included_by).resolve().parent

    lines = ["<!-- Generated by tests/scripts/report_spectral.py from the record spectral_field.py wrote,",
             "     docs/generated/spectral/measured.json. Do not edit. The figures beside it are redrawn by",
             "     the same run, and the captions' numbers come from the record the figures draw. -->", ""]
    tables = [("provenance", table_provenance), ("vacuum", table_vacuum), ("null", table_null), ("bank", table_bank),
              ("cost", table_cost), ("seeded", table_seeded), ("convergence", table_convergence),
              ("chirped", table_chirped), ("sweep", table_sweep)]
    needs = {"provenance": "provenance", "convergence": "seeded"}
    for name, fn in tables:
        if needs.get(name, name) in rec:
            lines += fn(rec)
    figures = [("vacuum", fig_vacuum), ("null", fig_null), ("bank", fig_bank), ("cost", fig_cost),
               ("seeded", fig_seeded), ("chirped", fig_chirped), ("sweep", fig_sweep)]
    diagrams = [("diagram-pad", diagram_pad), ("diagram-slippage", diagram_slippage), ("diagram-planes", diagram_planes)]
    with style():
        for name, fn in diagrams:
            fig, cap = fn(plt)
            image = out / f"{name}.svg"
            save(fig, image)
            plt.close(fig)
            rel = image.resolve().relative_to(link_dir).as_posix()
            lines += figure(f"fig-sp-{name}", rel, cap)
            print(f"{name}: {image.name}")
        for name, fn in figures:
            if name not in rec:
                continue
            fig, cap = fn(rec, plt)
            image = out / f"{name}.svg"
            save(fig, image)
            plt.close(fig)
            rel = image.resolve().relative_to(link_dir).as_posix()
            lines += figure(f"fig-sp-{name}", rel, cap)
            print(f"{name}: {image.name}")
    Path(args.page).write_text("\n".join(lines))
    print(f"wrote {args.page}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
