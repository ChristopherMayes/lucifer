#!/usr/bin/env python3
"""
The saturation demo's automated summary report (PDF). Reads everything from the demo's
work directory -- the Genesis .out.h5, both Bmad diag/stats files, the unaveraged
energy ledger, and the /usr/bin/time files -- recomputes the agreement checks, and
writes a multi-page PDF in which every figure carries the paragraph that explains how
to read it. The energy accounting gets its own page, as labeled decompositions and a
books-close overlay instead of overlapping curves.

Usage:  report_fel_saturation.py <workdir> [-o report.pdf]

Expects in <workdir>: AramisSat.out.h5, sat-avg.diag.txt, sat-unavg.diag.txt,
sat-unavg.ledger.txt, sat-{genesis,avg,unavg}.log.time. Needs numpy, h5py, matplotlib.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import time

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

M_ELECTRON = 0.51099895069e6
C = 299792458.0

BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#b0afa9"
INK, INK2 = "#0b0b0b", "#52514e"
CODES = {"genesis": ("Genesis4 (MPI)", "0.25", 2.6),
         "avg": ("Bmad averaged", ORANGE, 1.5),
         "unavg": ("Bmad unaveraged", BLUE, 1.5)}


def wall(wd, name):
    try:
        for line in (wd / f"{name}.log.time").read_text().splitlines():
            if line.startswith("real"):
                return float(line.split()[1])
    except OSError:
        return None


def load_all(wd):
    d = {}
    with h5py.File(wd / "AramisSat.out.h5") as h5:
        d["nslice"] = h5["Field/power"].shape[1]
        d["z"] = h5["Lattice/zplot"][:]
        d["p_genesis"] = h5["Field/power"][:]
        d["b_genesis"] = h5["Beam/bunching"][:]
        d["e_genesis"] = h5["Beam/energy"][:] * M_ELECTRON        # gamma -> eV
        d["es_genesis"] = h5["Beam/energyspread"][:] * M_ELECTRON
        d["i_genesis"] = h5["Beam/current"][0]
    for tag in ("avg", "unavg"):
        raw = np.loadtxt(wd / f"sat-{tag}.diag.txt").reshape(-1, d["nslice"], 12)
        d[f"p_{tag}"] = raw[:, :, 2]
        d[f"b_{tag}"] = raw[:, :, 4]
        d[f"e_{tag}"] = raw[:, :, 6]
        d[f"es_{tag}"] = raw[:, :, 7]
        d[f"i_{tag}"] = raw[0, :, 10]
    spacing = None
    for line in (wd / "sat-avg.diag.txt").open():
        if line.startswith("# slice_spacing"):
            spacing = float(line.split("=")[1])
        if not line.startswith("#"):
            break
    d["spacing"] = spacing
    d["led"] = np.loadtxt(wd / "sat-unavg.ledger.txt")       # z Eb U dk Uesc Uspont
    d["n"] = min(len(d["z"]), len(d["p_avg"]), len(d["p_unavg"]))
    d["times"] = {k: wall(wd, f"sat-{n}") for k, n in
                  (("genesis", "genesis"), ("avg", "avg"), ("unavg", "unavg"))}
    return d


def beam_given(d, tag):
    """-dE_beam(z) [J]: charge-weighted mean-energy drop summed over slices."""
    q = d[f"i_{tag}"] * d["spacing"] / C
    e = d[f"e_{tag}"][:d["n"]]
    return -((e - e[0]) * q).sum(axis=1)


def style(fig):
    for ax in fig.get_axes():
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, color=INK2, alpha=0.18, lw=0.6)


def caption(fig, text, y=0.02):
    fig.text(0.06, y, text, ha="left", va="bottom", fontsize=9.2, color=INK,
             wrap=True, linespacing=1.45,
             bbox=dict(boxstyle="square,pad=0.6", fc="#f7f6f3", ec="none"))


def page(pdf, fig):
    pdf.savefig(fig)
    plt.close(fig)


def cover(pdf, d, checks):
    fig = plt.figure(figsize=(8.5, 11))
    t = d["times"]
    lines = [
        ("The saturation demo: one SASE case, three trackers, one clock", 15, 1),
        ("", 10, 0),
        (f"Generated {time.strftime('%Y-%m-%d %H:%M')} from lucifer/examples/saturation_demo.", 10, 0),
        ("", 10, 0),
        ("The case.  Genesis's own Benchmark1-SASE configuration run to saturation: the", 10, 0),
        ("full 57 m 6-FODO Aramis line (12 undulator segments, aw = 0.85 rms helical,", 10, 0),
        ("15 mm period, 5.8 GeV, resonant at 1 Angstrom), dark start, growth from shot", 10, 0),
        (f"noise alone. {d['nslice']} slices x 2048 particles, 255^2 field grid. All three", 10, 0),
        ("trackers start from IDENTICAL initial dumps written by an untimed Genesis prep", 10, 0),
        ("run; wall clocks come from one external clock (/usr/bin/time).", 10, 0),
        ("", 10, 0),
        ("Wall clock (12 performance cores each)", 11, 1),
        (f"    Genesis 1.3 v4, 12 MPI ranks              {t['genesis']:8.1f} s", 10, 2),
        (f"    Bmad averaged (bmad_standard default)     {t['avg']:8.1f} s   "
         f"({t['genesis']/t['avg']:.2f}x faster than Genesis)", 10, 2),
        (f"    Bmad unaveraged (fel_tracking selection)  {t['unavg']:8.1f} s   "
         f"({t['unavg']/t['avg']:.0f}x the averaged mode)", 10, 2),
        ("    Each code computes its own in-run diagnostics; the Bmad runs write the", 10, 0),
        ("    full moment-matrix stats file while Genesis writes scalar columns.", 10, 0),
        ("", 10, 0),
        ("The answer checks (recomputed here, must pass before timings mean anything)", 11, 1),
        (f"    exit total power, Genesis:            {checks['p_gen']:.4e} W", 10, 2),
        (f"    Bmad averaged:    {checks['p_avg']:.4e} W    rel. difference {checks['rel']:.2e}   (<= 0.15)", 10, 2),
        (f"    Bmad unaveraged:  {checks['p_unavg']:.4e} W    |ln ratio| {checks['lnr']:.3f}        (<= 1.0)", 10, 2),
        ("", 10, 0),
        ("What to expect in the pages that follow", 11, 1),
        ("    1. Gain curves: the averaged mode lies on Genesis through eight decades of z", 10, 0),
        ("       and three of power (4.9e-4 at the exit). The unaveraged mode -- an", 10, 0),
        ("       independent integrator through the real undulator field, with the", 10, 0),
        ("       coupling factor nowhere in its inputs -- reproduces startup, gain shape", 10, 0),
        ("       and saturation location, riding ~2%/m above the averaged codes while its", 10, 0),
        ("       beam gives up ~14x more energy -- see page 4.", 10, 0),
        ("    2. Pulse structure: per-slice exit power and bunching evolution.", 10, 0),
        ("    3. Beam evolution: energy loss and spread -- where the unaveraged mode's", 10, 0),
        ("       extra physics is visible on the beam itself.", 10, 0),
        ("    4. Energy accounting: the window is an OPEN system (slippage transmits", 10, 0),
        ("       light out of its head), so window energy never equals beam energy given.", 10, 0),
        ("       The decomposition per code, the unaveraged ledger closing exactly, and", 10, 0),
        ("       WHY the two Bmad modes differ: the averaged/KMR model emits spontaneous", 10, 0),
        ("       radiation without debiting the beam; the unaveraged mode conserves and", 10, 0),
        ("       pays, at the analytically correct magnitude.", 10, 0),
        ("", 10, 0),
        ("The measured levels and methodology live in lucifer/doc/validation.md", 10, 0),
        ("(sections: The saturation demo, Diagnostic output) and the physics manual", 10, 0),
        ("(lucifer/doc/fel-physics.md, sec-unaveraged, sec:stats).", 10, 0),
    ]
    y = 0.94
    for text, size, kind in lines:
        weight = "bold" if kind == 1 else "normal"
        family = "monospace" if kind == 2 else "sans-serif"
        fig.text(0.08, y, text, fontsize=size, weight=weight, family=family, color=INK)
        y -= 0.021 if size <= 11 else 0.030
    page(pdf, fig)


def gain_page(pdf, d):
    n = d["n"]
    z = d["z"][:n]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 11), height_ratios=[2, 1])
    fig.subplots_adjust(top=0.94, bottom=0.30, hspace=0.25, left=0.11, right=0.95)
    tot = {k: d[f"p_{k}"][:n].sum(axis=1) for k in CODES}
    for k, (lab, col, lw) in CODES.items():
        a1.semilogy(z, np.maximum(tot[k], 1e-30), color=col, lw=lw, label=lab)
    a1.set_ylim(bottom=max(tot["genesis"][tot["genesis"] > 0].min() * 0.3, 1e-2))
    a1.set_ylabel("total radiation power (W)")
    a1.legend(frameon=False)
    a1.set_title("Gain curves", loc="left", weight="bold")
    for k in ("avg", "unavg"):
        lab, col, lw = CODES[k]
        with np.errstate(divide="ignore"):
            a2.semilogy(z, np.abs(tot[k] / tot["genesis"] - 1), color=col, lw=lw, label=lab)
    a2.set_xlabel("z (m)")
    a2.set_ylabel("|P / P_Genesis − 1|")
    a2.legend(frameon=False)
    style(fig)
    caption(fig,
        "How to read this.  Top: total radiation power against distance, all three\n"
        "trackers from the same initial particles and field. The orange (Bmad averaged)\n"
        "curve lies on the gray (Genesis) curve everywhere -- the bottom panel shows the\n"
        "difference staying at the 1e-4..1e-3 level through saturation. The blue curve is\n"
        "the unaveraged mode: no period averaging, the electrons integrated through the\n"
        "real undulator field with the coupling factor nowhere in its inputs. It\n"
        "reproduces the startup level, the gain shape and the saturation location, and\n"
        "rides ~2%/m above the averaged codes through the exponential regime: that excess\n"
        "is the beam's real shot-noise radiation, which the unaveraged dynamics resolves\n"
        "and the period-averaged model does not track (it injects shot noise once, at\n"
        "load time). The channel is measured independent of particle count and of\n"
        "integrator resolution -- physics, not noise. Its beam-side cost\n"
        "is on the beam-evolution and energy-accounting pages.")
    page(pdf, fig)


def pulse_page(pdf, d):
    n = d["n"]
    z = d["z"][:n]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 11), height_ratios=[1, 1])
    fig.subplots_adjust(top=0.94, bottom=0.28, hspace=0.25, left=0.11, right=0.95)
    s = np.arange(1, d["nslice"] + 1)
    for k, (lab, col, lw) in CODES.items():
        a1.step(s, d[f"p_{k}"][n - 1], where="mid", color=col, lw=lw, label=lab)
    a1.set_xlabel("slice (time-window order; higher index = pulse head)")
    a1.set_ylabel(f"power at z = {z[-1]:.1f} m (W)")
    a1.legend(frameon=False)
    a1.set_title("Pulse structure at the exit", loc="left", weight="bold")
    for k, (lab, col, lw) in CODES.items():
        a2.plot(z, d[f"b_{k}"][:n].mean(axis=1), color=col, lw=lw, label=lab)
    a2.set_xlabel("z (m)")
    a2.set_ylabel("bunching |b| (slice mean)")
    a2.legend(frameon=False)
    style(fig)
    caption(fig,
        "How to read this.  Top: the SASE pulse at the undulator exit, slice by slice\n"
        "across the time window. The first ~12 slices are dark: light continuously slips\n"
        "forward through the bunch (one slice spacing per three undulator periods), so\n"
        "fresh vacuum enters at the tail while radiation leaves at the head -- the\n"
        "leaving light is accounted on the energy page. Genesis and the Bmad averaged\n"
        "mode agree slice by slice; the unaveraged mode carries more power in every lit\n"
        "slice (the same +0.6 ln as the gain page, distributed across the pulse, not an\n"
        "artifact of a few slices). Bottom: the microbunching that drives the FEL, rising\n"
        "from the shot-noise level ~0.006 to a few percent at saturation.")
    page(pdf, fig)


def beam_page(pdf, d):
    n = d["n"]
    z = d["z"][:n]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 11))
    fig.subplots_adjust(top=0.94, bottom=0.30, hspace=0.25, left=0.11, right=0.95)
    for k, (lab, col, lw) in CODES.items():
        e = d[f"e_{k}"][:n].mean(axis=1)
        a1.plot(z, (e - e[0]) / 1e6, color=col, lw=lw, label=lab)
    a1.set_ylabel(r"mean energy change $\Delta\langle E\rangle$ (MeV)")
    a1.legend(frameon=False)
    a1.set_title("Beam evolution", loc="left", weight="bold")
    for k, (lab, col, lw) in CODES.items():
        a2.plot(z, d[f"es_{k}"][:n].mean(axis=1) / 1e6, color=col, lw=lw, label=lab)
    a2.set_xlabel("z (m)")
    a2.set_ylabel(r"rms energy spread $\sigma_E$ (MeV)")
    a2.legend(frameon=False)
    style(fig)
    caption(fig,
        "How to read this.  Top: the energy the beam gives up. Genesis and the Bmad\n"
        "averaged mode show essentially only the coherent FEL exchange. The unaveraged\n"
        "beam loses ~14x more, and the reason is a MODEL difference, not a bug: both\n"
        "models emit spontaneous shot-noise radiation of the same magnitude, but the\n"
        "period-averaged (KMR) model does not charge the beam for it -- its field step\n"
        "adds 2S while the particles are kicked by E, so the 4|S|^2 part of the field\n"
        "energy is created rather than taken from the beam (measured factor 134 on a dark\n"
        "segment; Genesis's optional &sponrad module exists precisely to add the missing\n"
        "loss by hand, and is off by default). The unaveraged mode conserves energy by\n"
        "construction, so its beam pays -- at a rate that agrees to 8% with the analytic\n"
        "spontaneous power (2/3)r_e g^2 ku^2 aw^2 restricted to the grid's angular\n"
        "acceptance. Neither model yet carries the ~90% of spontaneous power radiated\n"
        "outside that acceptance. Bottom: the same physics as heating --\n"
        "the unaveraged spread grows faster (spontaneous diffusion) on top of the\n"
        "FEL-induced spread both models share.")
    page(pdf, fig)


def energy_page(pdf, d, checks):
    n = d["n"]
    fig = plt.figure(figsize=(8.5, 11))
    a1 = fig.add_axes([0.17, 0.72, 0.76, 0.21])
    a2 = fig.add_axes([0.17, 0.52, 0.76, 0.13])

    # -- Panel 1: the decomposition at the exit, one bar pair per code. ------------
    given = {"genesis": None, "avg": beam_given(d, "avg")[-1], "unavg": None}
    # Genesis: charge-weighted from its own arrays.
    qg = d["i_genesis"] * d["spacing"] / C
    eg = d["e_genesis"][:n]
    given["genesis"] = -((eg - eg[0]) * qg).sum(axis=1)[-1]
    led = d["led"]
    given["unavg"] = -(led[-1, 1] - led[0, 1])
    window = {k: d[f"p_{k}"][n - 1].sum() * d["spacing"] / C for k in CODES}
    esc_unavg, spont_unavg = led[-1, 4], led[-1, 5]

    ypos = {"genesis": 2, "avg": 1, "unavg": 0}
    for k, y in ypos.items():
        lab, col, _ = CODES[k]
        g = given[k]
        a1.barh(y + 0.18, g, height=0.30, color=col, alpha=0.9)
        a1.text(g * 1.02, y + 0.18, f"beam gave {g:.2e} J", va="center", fontsize=8)
        w = window[k]
        a1.barh(y - 0.18, w, height=0.30, color=col, alpha=0.45)
        if k == "unavg":
            a1.barh(y - 0.18, esc_unavg, height=0.30, left=w, color=col, alpha=0.2,
                    hatch="//", edgecolor=col)
            a1.text(0.02 * a1.get_xlim()[1] if False else w + esc_unavg + 0.01e-8, y - 0.55,
                    f"window {w:.2e} + escaped {esc_unavg:.2e} J (bookkept exactly, the ledger)",
                    va="center", ha="right", fontsize=8)
        else:
            inf = g - w
            a1.barh(y - 0.18, inf, height=0.30, left=w, color=col, alpha=0.2,
                    hatch="//", edgecolor=col)
            a1.text(g * 1.02, y - 0.18,
                    f"window {w:.2e} + escaped {inf:.2e} J (inferred)", va="center", fontsize=8)
    a1.set_yticks(list(ypos.values()), [CODES[k][0] for k in ypos])
    a1.set_ylim(-0.9, 2.6)
    a1.set_xlabel("energy at the exit (J)")
    a1.set_xlim(0, max(given.values()) * 1.55)
    a1.set_title("Energy accounting at the exit", loc="left", weight="bold")

    # -- Panel 2: the unaveraged books closing along z. ----------------------------
    zl = led[:, 1 - 1]
    given_z = -(led[:, 1] - led[0, 1])
    account_z = (led[:, 2] - led[0, 2]) + led[:, 4] - led[:, 5]
    if led.shape[1] > 6:
        account_z = account_z + led[:, 6]      # E_radiated (bmad_com radiation on).
    a2.plot(zl, given_z, color=BLUE, lw=2.6, label="beam energy given, −ΔE_beam")
    a2.plot(zl, account_z, color=INK, lw=1.0, ls="--",
            label="accounted: U_window + U_escaped − U_spontaneous")
    a2.set_xlabel("z (m)")
    a2.set_ylabel("energy (J)")
    a2.legend(frameon=False, loc="upper left")
    a2.set_title("The unaveraged ledger closes", loc="left", weight="bold")
    resid = np.abs((given_z - account_z)).max()
    turn = np.abs(np.diff(led[:, 2])).sum() + abs(led[-1, 4]) + abs(led[-1, 5])
    a2.text(0.98, 0.08, f"max residual {resid:.2e} J = {resid/turn:.1e} of turnover",
            transform=a2.transAxes, ha="right", fontsize=8.5, color=INK2)
    style(fig)

    caption(fig,
        "How to read this.  The time window travels with the beam, and light travels\n"
        "faster: every three undulator periods the radiation slips one slice forward, so\n"
        "light continuously exits the head of the window and is gone from the simulated\n"
        "volume. The window is an open box: energy in the window never equals the energy\n"
        "the beam gave.\n"
        "\n"
        "Top: per tracker, the upper bar is the total energy the beam handed to radiation\n"
        "by the exit; the lower bar shows where it is now -- still in the window (solid)\n"
        "or slipped out ahead (hatched). Genesis and the averaged mode keep ~72%: their\n"
        "radiation is made late, near saturation, and has not had time to escape. The\n"
        "unaveraged beam gave ~14x more and keeps only ~10%: it is charged for the\n"
        "spontaneous radiation both models emit (the averaged model creates that field\n"
        "energy without debiting the beam, page 4's caption), and\n"
        "that emission is spread uniformly along all 57 m, so anything radiated more than\n"
        "~4 m of undulator ago has fully crossed the 96-slice window and left. For\n"
        "Genesis and the averaged mode the escaped part is INFERRED\n"
        "(given minus window); for the unaveraged mode it is bookkept exactly, slice by\n"
        "slice, as the light leaves (the ledger's U_escaped column).\n"
        "\n"
        "Bottom: the proof the accounting is complete. The beam's cumulative loss (blue)\n"
        "and the accounted sum window + escaped − spontaneous credit (dashed) lie on top\n"
        "of each other over the whole run; U_spontaneous is the one term the kick/deposit\n"
        "energy duality does not charge to the beam (the substep's own |src|² deposit --\n"
        "physically, spontaneous emission; manual eq:ledger).\n"
        f"\n"
        f"At the exit -- Genesis: gave {given['genesis']:.2e} J, keeps "
        f"{100*window['genesis']/given['genesis']:.0f}%.  Averaged: gave {given['avg']:.2e} J, keeps "
        f"{100*window['avg']/given['avg']:.0f}%.  Unaveraged: gave {given['unavg']:.2e} J, keeps\n"
        f"{100*window['unavg']/given['unavg']:.0f}%, escaped {esc_unavg:.2e} J, spontaneous credit "
        f"{spont_unavg:.2e} J.")
    page(pdf, fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    out = args.out or str(wd / "sat-demo-report.pdf")

    d = load_all(wd)
    n = d["n"]
    checks = {"p_gen": d["p_genesis"][n - 1].sum(), "p_avg": d["p_avg"][n - 1].sum(),
              "p_unavg": d["p_unavg"][n - 1].sum()}
    checks["rel"] = abs(checks["p_avg"] - checks["p_gen"]) / checks["p_gen"]
    checks["lnr"] = abs(math.log(checks["p_unavg"] / checks["p_gen"]))

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    with PdfPages(out) as pdf:
        cover(pdf, d, checks)
        gain_page(pdf, d)
        pulse_page(pdf, d)
        beam_page(pdf, d)
        energy_page(pdf, d, checks)
    print(f"wrote {out}")
    if checks["rel"] > 0.15 or checks["lnr"] > 1.0:
        print("WARNING: agreement checks exceed documented levels; see the cover page.")
        return 1
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
