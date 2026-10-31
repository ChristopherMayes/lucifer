#!/usr/bin/env python3
"""
Draw docs/figures/window.svg, the glossary's figure of the beam and the time window.

  python docs/figures/window_figure.py [--out docs/figures/window.svg]

The figure is a cartoon with no data in it: a window of sixteen slices in window order,
the bunch over the low-index slices and the slippage headroom at the head, light slipping
toward the head and leaving there, a slab of four slices, and one slice enlarged to show
its wavelengths, its macroparticles grouped into beamlets, and its one transverse grid.
Every number drawn here is for the picture. The counts a run uses come from the input
file and the derivations on the physics page. The script is run by hand and the SVG it
writes is committed, so the two are kept together.
"""
import argparse
import math
from pathlib import Path

# Tol's bright palette, readable in print and by colorblind readers.
INK = "#333333"
BLUE = "#4477aa"
BLUE_FILL = "#c9d8ea"
ORANGE = "#ee7733"
TEAL = "#009988"
TEAL_FILL = "#cce9e5"
GREY = "#bbbbbb"
FONT = "font-family: Helvetica, Arial, sans-serif"

W, H = 960, 540
NSLICE = 14                # slices drawn
X0, BOX = 130, 50          # left edge and width of a slice box
Y0, Y1 = 180, 240          # top and bottom of the slice row
BUNCH_LAST = 9             # the bunch spans slices 1 to this one
SLAB = (8, 11)             # the slab, first and last slice
MAG = 5                    # the slice enlarged below
NWAVE = 4                  # wavelengths per slice, in the picture
PANEL = (130, 530, 360, 460)   # x0, x1, y0, y1 of the enlarged slice
BEAMLET = 4


def xb(i):
    """Left edge of slice i."""
    return X0 + (i - 1) * BOX


def text(x, y, s, size=13, anchor="middle", fill=INK, style="", weight=""):
    w = f" font-weight: {weight};" if weight else ""
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
            f'fill="{fill}" style="{FONT};{w}{style}">{s}</text>')


def arrow(x0, y0, x1, y1, color=INK, width=1.4, marker="arrow"):
    return (f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" stroke="{color}" '
            f'stroke-width="{width}" marker-end="url(#{marker})"/>')


def bracket(x0, x1, y, depth=8, color=INK):
    """A square bracket below y, opening downward, with its tick at the middle."""
    return (f'<path d="M{x0},{y} v{depth} H{x1} v{-depth}" fill="none" stroke="{color}" '
            f'stroke-width="1.2"/>')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(Path(__file__).with_name("window.svg")))
    args = ap.parse_args()

    out = []
    a = out.append
    a(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
    a('<defs>')
    for name, color in (("arrow", INK), ("arrow-orange", ORANGE), ("arrow-blue", BLUE)):
        a(f'<marker id="{name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
          f'markerHeight="7" orient="auto-start-reverse">'
          f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')
    a('</defs>')
    a(f'<rect width="{W}" height="{H}" fill="white"/>')

    # The slippage arrow, above everything.
    a('<!-- slippage -->')
    a(text(W / 2, 24, "slippage: light moves toward the head, one slice every "
           "n_wavelength undulator periods", fill=ORANGE))
    a(arrow(250, 42, 760, 42, ORANGE, 2, "arrow-orange"))

    # The bunch: a Gaussian current over slices 1 to BUNCH_LAST, sitting low in the window.
    a('<!-- the bunch -->')
    xc = xb(1) + BUNCH_LAST * BOX / 2
    half = BUNCH_LAST * BOX / 2
    sigma = half / 3.6
    peak = 80
    pts = []
    for k in range(0, 121):
        x = xb(1) + k * BUNCH_LAST * BOX / 120
        y = Y0 - peak * math.exp(-0.5 * ((x - xc) / sigma) ** 2)
        pts.append(f"{x:.1f},{y:.1f}")
    a(f'<path d="M{xb(1)},{Y0} L' + " ".join(pts) + f' L{xb(BUNCH_LAST + 1)},{Y0} z" '
      f'fill="{BLUE_FILL}" stroke="{BLUE}" stroke-width="1.5"/>')
    a(text(xc, 66, "the bunch: every particle", fill=BLUE, weight="bold"))
    a(text(xc, 82, "beam.bunch_charge, beam.sig_z", size=12, fill=BLUE,
           style="font-family: Menlo, monospace"))

    # Headroom: the slices ahead of the bunch, where the slippage goes.
    a('<!-- headroom -->')
    hx0, hx1 = xb(BUNCH_LAST + 1), xb(NSLICE + 1)
    a(f'<path d="M{hx0},{Y0 - 14} v-8 H{hx1} v8" fill="none" stroke="{INK}" stroke-width="1.2"/>')
    a(text((hx0 + hx1) / 2, 128, "slippage headroom", size=13))
    a(text((hx0 + hx1) / 2, 144, "what the line slips, in whole slices", size=12))

    # The slice row, in window order.
    a('<!-- the window -->')
    for i in range(1, NSLICE + 1):
        fill = TEAL_FILL if SLAB[0] <= i <= SLAB[1] else "white"
        a(f'<rect x="{xb(i)}" y="{Y0}" width="{BOX}" height="{Y1 - Y0}" fill="{fill}" '
          f'stroke="{INK}" stroke-width="1"/>')
    labels = {1: "1", 2: "2", 3: "3", MAG: "i", 7: "…", NSLICE: "n_slice"}
    for i, s in labels.items():
        size = 11 if s == "n_slice" else 13
        style = "font-family: Menlo, monospace" if s == "n_slice" else ""
        a(text(xb(i) + BOX / 2, (Y0 + Y1) / 2 + 4, s, size=size, style=style))
    a(f'<rect x="{xb(MAG)}" y="{Y0}" width="{BOX}" height="{Y1 - Y0}" fill="none" '
      f'stroke="{BLUE}" stroke-width="2"/>')

    # Tail and head, with what enters and what leaves.
    a('<!-- tail and head -->')
    a(arrow(40, (Y0 + Y1) / 2, X0 - 4, (Y0 + Y1) / 2, GREY, 1.6))
    a(text(X0 / 2 - 2, Y0 - 24, "vacuum enters", size=12))
    a(text(X0 / 2 - 2, Y1 + 20, "tail", size=13, weight="bold"))
    a(text(X0 / 2 - 2, Y1 + 36, "arrives last,", size=12))
    a(text(X0 / 2 - 2, Y1 + 50, "lower z", size=12))
    xe = xb(NSLICE + 1)
    a(arrow(xe + 4, (Y0 + Y1) / 2, xe + 76, (Y0 + Y1) / 2, ORANGE, 2, "arrow-orange"))
    a(text(xe + 60, Y0 - 38, "escaped field,", size=12, fill=ORANGE))
    a(text(xe + 60, Y0 - 24, "banked", size=12, fill=ORANGE))
    a(text(xe + 60, Y1 + 20, "head", size=13, weight="bold"))
    a(text(xe + 60, Y1 + 36, "arrives first,", size=12))
    a(text(xe + 60, Y1 + 50, "higher z", size=12))
    a(text(W / 2, H - 16, "window order: slice index 1 to n_slice, slices Δ apart, "
           "window.n_slice of them", size=13))

    # The slab.
    a('<!-- the slab -->')
    sx0, sx1 = xb(SLAB[0]), xb(SLAB[1] + 1)
    a(bracket(sx0, sx1, Y1 + 4, color=TEAL))
    a(text((sx0 + sx1) / 2, Y1 + 26, "slab: the slices a frame writes", size=12, fill=TEAL))
    a(text((sx0 + sx1) / 2, Y1 + 40, "run.dump_slice_first .. dump_slice_last", size=11,
           fill=TEAL, style="font-family: Menlo, monospace"))

    # One slice enlarged.
    a('<!-- one slice, enlarged -->')
    px0, px1, py0, py1 = PANEL
    a(f'<path d="M{xb(MAG)},{Y1} L{px0},{py0} M{xb(MAG + 1)},{Y1} L{px1},{py0}" '
      f'stroke="{BLUE}" stroke-width="1" stroke-dasharray="4 3" fill="none"/>')
    a(f'<rect x="{px0}" y="{py0}" width="{px1 - px0}" height="{py1 - py0}" fill="white" '
      f'stroke="{BLUE}" stroke-width="1.5"/>')
    lam = (px1 - px0) / NWAVE
    for k in range(1, NWAVE):
        a(f'<line x1="{px0 + k * lam}" y1="{py0}" x2="{px0 + k * lam}" y2="{py1}" '
          f'stroke="{GREY}" stroke-width="1" stroke-dasharray="3 3"/>')
    # The field's phase, one sinusoid per wavelength, faint.
    pts = []
    for k in range(0, 161):
        x = px0 + k * (px1 - px0) / 160
        y = (py0 + py1) / 2 - 10 * math.sin(2 * math.pi * (x - px0) / lam)
        pts.append(f"{x:.1f},{y:.1f}")
    a(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{ORANGE}" stroke-width="1" '
      f'opacity="0.6"/>')
    # Beamlets: BEAMLET macroparticles each, phases spaced evenly over one wavelength, one
    # transverse sample (one row) each.
    rows = [py0 + 18, py0 + 34, py0 + 50, py0 + 66, py0 + 82]
    starts = [0.12, 1.22, 2.10, 3.15, 2.20]
    for j, (y, s) in enumerate(zip(rows, starts)):
        color = TEAL if j == 2 else BLUE
        for m in range(BEAMLET):
            x = px0 + (s + m / BEAMLET) * lam
            a(f'<circle cx="{x:.1f}" cy="{y}" r="3.2" fill="{color}"/>')
    # Dimension lines: one wavelength and the slice.
    a(f'<path d="M{px0},{py0 - 10} v-6 M{px0 + lam},{py0 - 10} v-6" stroke="{INK}" stroke-width="1"/>')
    a(f'<line x1="{px0}" y1="{py0 - 13}" x2="{px0 + lam}" y2="{py0 - 13}" stroke="{INK}" '
      f'stroke-width="1" marker-start="url(#arrow)" marker-end="url(#arrow)"/>')
    a(text(px0 + lam / 2, py0 - 20, "λ_s = field.lambda0", size=12))
    a(f'<path d="M{px0},{py1 + 10} v6 M{px1},{py1 + 10} v6" stroke="{INK}" stroke-width="1"/>')
    a(f'<line x1="{px0}" y1="{py1 + 13}" x2="{px1}" y2="{py1 + 13}" stroke="{INK}" '
      f'stroke-width="1" marker-start="url(#arrow)" marker-end="url(#arrow)"/>')
    a(text((px0 + px1) / 2, py1 + 34, "one slice: Δ = window.n_wavelength · λ_s "
           "(4 here), the slice spacing", size=12))

    # The legend of the enlarged slice, to its right.
    a('<!-- what one slice holds -->')
    lx = px1 + 30
    ly = py0 + 6
    a(f'<circle cx="{lx + 6}" cy="{ly}" r="3.2" fill="{BLUE}"/>')
    a(text(lx + 18, ly + 4, "macroparticle: stands for weight electrons,", size=12, anchor="start"))
    a(text(lx + 18, ly + 18, "beam.n_particle of them a slice", size=12, anchor="start"))
    ly += 38
    for m in range(BEAMLET):
        a(f'<circle cx="{lx + 6 + m * 9}" cy="{ly}" r="3.2" fill="{TEAL}"/>')
    a(text(lx + 46, ly + 4, "beamlet: beam.beamlet_size of them, one transverse",
           size=12, anchor="start"))
    a(text(lx + 46, ly + 18, "sample, one energy, phases evenly over 2π", size=12, anchor="start"))
    ly += 44
    g = 5
    for k in range(g + 1):
        a(f'<line x1="{lx}" y1="{ly - 12 + k * 5}" x2="{lx + g * 5}" y2="{ly - 12 + k * 5}" '
          f'stroke="{ORANGE}" stroke-width="0.8"/>')
        a(f'<line x1="{lx + k * 5}" y1="{ly - 12}" x2="{lx + k * 5}" y2="{ly - 12 + g * 5}" '
          f'stroke="{ORANGE}" stroke-width="0.8"/>')
    a(text(lx + 36, ly + 1, "the field: one transverse grid a slice (a plane),", size=12, anchor="start"))
    a(text(lx + 36, ly + 15, "grid_n_pts² points to ± grid_half_width", size=12, anchor="start"))
    ly += 42
    a(text(lx, ly, "A steady-state run is one slice.", size=12, anchor="start"))

    a('</svg>')
    Path(args.out).write_text("\n".join(out) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
