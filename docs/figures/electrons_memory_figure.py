#!/usr/bin/env python3
"""
Draw docs/figures/electrons-memory.svg, the number of particles and the memory of the
particles in the electrons mode as functions of the charge in the window.

  python docs/figures/electrons_memory_figure.py [--out docs/figures/electrons-memory.svg]

The figure holds no measurement. The number of particles is the charge divided by the
elementary charge. The host holds a particle in 64 bytes of FP64 with a margin of one
sixteenth, 68 bytes per electron. The device holds a particle in 32 bytes of FP32. The
largest number of particles a Bmad bunch indexes is 2**31 - 1. With --table the script
prints the rows of the table on docs/electrons-mode.md. The script is run by hand and the
SVG it writes is committed, so the two are kept together.
"""
import argparse
import math
from pathlib import Path

E_CHARGE = 1.602176634e-19      # C, exact in the SI
HOST_BYTES = 68                 # 64 bytes of FP64 and the margin of one sixteenth
DEVICE_BYTES = 32               # six FP32 values and the 8-byte phase
N_BUNCH_MAX = 2**31 - 1         # the largest default integer
TABLE_PC = [0.1, 1, 10, 100, 1000]


def n_electrons(q_pc):
    return q_pc * 1e-12 / E_CHARGE


def gb(n_bytes):
    """Gigabytes to three significant figures, in plain decimal."""
    x = n_bytes / 1e9
    digits = max(0, 2 - int(math.floor(math.log10(x))))
    return f"{x:.{digits}f}"


def table():
    rows = ["| charge in the window | particles | host, FP64 | device, FP32 |", "|---|---|---|---|"]
    for q in TABLE_PC:
        n = n_electrons(q)
        mant, expo = f"{n:.2e}".split("e")
        rows.append(f"| {q:g} pC | ${mant}\\times10^{{{int(expo)}}}$ | "
                    f"{gb(n * HOST_BYTES)} GB | {gb(n * DEVICE_BYTES)} GB |")
    return "\n".join(rows)


def draw(out):
    import matplotlib
    matplotlib.use("svg")
    matplotlib.rcParams["svg.hashsalt"] = "electrons-memory"
    matplotlib.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt

    q = [0.05, 2000.0]
    q_max = N_BUNCH_MAX * E_CHARGE * 1e12
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(q, [n_electrons(x) * HOST_BYTES / 1e9 for x in q], color="#4477aa",
              label="host, FP64, 68 bytes per electron")
    ax.loglog(q, [n_electrons(x) * DEVICE_BYTES / 1e9 for x in q], color="#ee7733",
              label="device, FP32, 32 bytes per electron")
    ax.axvline(q_max, color="#bbbbbb", linestyle="--")
    ax.text(q_max * 0.88, 0.012, r"$2^{31}-1$ particles", rotation=90, ha="right", va="bottom",
            color="#333333", fontsize=9)
    ax.set_xlim(*q)
    ax.set_xlabel("charge in the window (pC)")
    ax.set_ylabel("memory of the particles (GB)")
    ax.grid(True, which="major", color="#dddddd", linewidth=0.6)
    ax.legend(frameon=False, loc="upper left")
    top = ax.secondary_xaxis("top", functions=(n_electrons, lambda n: n * E_CHARGE * 1e12))
    top.set_xlabel("particles")
    fig.tight_layout()
    fig.savefig(out, metadata={"Date": None, "Creator": None})


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(here / "electrons-memory.svg"))
    ap.add_argument("--table", action="store_true", help="print the table and draw nothing")
    args = ap.parse_args()
    if args.table:
        print(table())
    else:
        draw(args.out)


if __name__ == "__main__":
    main()
