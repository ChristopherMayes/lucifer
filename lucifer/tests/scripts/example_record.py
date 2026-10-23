#!/usr/bin/env python3
"""
Record what an example's run measured, from its stats file and its log.

  example_record.py <deck.in> <log>              # rows of the record, to stdout
  example_record.py --sample <log> --out <md>    # the sample output the user guide shows

The first form reads the deck for its out_root, opens the stats file the run wrote
beside it and the log the runner kept, and prints one row per quantity:

  <directory> <tab> <deck> <tab> <quantity> <tab> <value> <tab> <unit>

The quantities are one fixed set for every deck, so a page can compare its decks column
by column and a number a README once typed from a run comes from the run itself:

  exit_power, exit_pulse_energy          the window totals at the last record
  exit_power_per_slice, exit_power_slice_spread  the exit power over the slices, and the
                                         rms over the mean of the per-slice powers
  peak_power, peak_power_s               the window power's maximum along the line
  exit_power_inside, exit_power_outside  the window power inside and outside the split
                                         angle at the last record that has the split, and
                                         exit_power_inside_s names that record's s when
                                         it is not the last
  wide_angle_ratio, wide_angle_ratio_exit  outside over inside, the run's maximum and its
                                         value at the exit, as the program reports them
  mode_peak_inside, mode_peak_s          the inside power's maximum and where it sits
  split_angle                            the angle the split is taken at
  initial_bunching, exit_bunching, exit_bunching_max
                                         the mean of |b| over slices at the start and the
                                         exit, and the largest slice's at the exit
  bunching_max, bunching_max_s           the maximum along the line of the mean |b| over
                                         slices, and where it sits, the saturation point
  exit_power_slice_max                   the largest slice's power at the exit
  bunching_peak                          the mean |b| at the mode's peak, as the program
                                         reports it
  initial_energy, exit_energy, energy_change
                                         the mean of the slice energies at the start and
                                         the exit, and the difference
  exit_energy_span                       the span of the slice mean energies at the exit,
                                         head to tail
  gain_length                            the shortest e-folding length of the inside
                                         power over any five consecutive records that
                                         have the split, the exponential regime
  power_after_first_element, power_after_first_element_per_slice
                                         the window power at the first element end, and
                                         the same over the slices
  initial_sigma_energy, exit_sigma_energy  the mean over slices of the slice energy spread
  n_slice, n_record, slice_spacing, n_particle, grid_n_pts, grid_cell, fel_steps
  grid_derived, grid_derived_from        whether setup derived the grid, and the count it
                                         rounded up from
  filter_mode_angle, filter_rho_angle, filter_ratio, filter_edge
                                         the source filter's two angles, their ratio and
                                         the edge, as setup printed them

A float prints with five significant digits, which is what a run reproduces, and a
count as itself. A quantity the run does not give, such as a gain length on a line with
no exponential regime, is left out rather than written as zero.

The second form copies the lines of a log that describe the run and its result, the
banner, the beam, the radiation, the work, the load, the exit and the split, into a
fenced block for the user guide to include. Progress rows carry wall times and are left
out, since a wall time belongs to a machine.

run_examples.sh runs both forms and writes the rows to examples/measured.tsv, which the
keystone regenerates and diffs with the generated pages.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
import sys

import h5py
import numpy as np

SAMPLE_PREFIXES = (" Lucifer --", " Beam ", " Radiation ", " Work ", " Load ", " Exit ",
                   " Split ", "             mode peak", "             last record")


def out_root(deck: pathlib.Path) -> str:
    for line in deck.read_text().splitlines():
        m = re.match(r"\s*(?:global%)?out_root\s*=\s*['\"]([^'\"]+)['\"]", line)
        if m:
            return m.group(1)
    raise SystemExit(f"{deck}: no out_root")


def fmt(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    f = float(v)
    if math.isnan(f) or math.isinf(f):
        return repr(f)
    return f"{f:.4e}"


def from_log(log: pathlib.Path):
    """The setup quantities a log states, as (name, value, unit)."""
    out = []
    text = log.read_text() if log.exists() else ""
    m = re.search(r"grid (\d+) points of ([0-9.]+) (um|nm|mm)", text)
    if m:
        scale = {"um": 1e-6, "nm": 1e-9, "mm": 1e-3}[m.group(3)]
        out.append(("grid_cell", float(m.group(2)) * scale, "m"))
    m = re.search(r"Grid (\d+) points, derived:.*?rounded up from (\d+)", text)
    if m:
        out += [("grid_derived", True, ""), ("grid_derived_from", int(m.group(2)), "")]
    elif re.search(r"Grid \d+ points, derived", text):
        out.append(("grid_derived", True, ""))
    else:
        out.append(("grid_derived", False, ""))
    m = re.search(r"Work\s+(\d+) FEL steps", text)
    if m:
        out.append(("fel_steps", int(m.group(1)), ""))
    m = re.search(r"Source filter: mode angle\s+([0-9.E+-]+) rad, rho angle\s+([0-9.E+-]+) rad, ratio\s+([0-9.E+-]+)", text)
    if m:
        out += [("filter_mode_angle", float(m.group(1)), "rad"), ("filter_rho_angle", float(m.group(2)), "rad"),
                ("filter_ratio", float(m.group(3)), "")]
    m = re.search(r"(?:puts the edge at|Edge at)\s+([0-9.E+-]+)\s*rad", text)
    if m:
        out.append(("filter_edge", float(m.group(1)), "rad"))
    return out


def gain_length(s, p_in):
    """The shortest e-folding length of the inside power over any five consecutive records
    that have the split, which is the exponential regime wherever the line has one. A
    window of fixed power fractions would land in the lethargy of one line and past the
    saturation of another."""
    ok = np.isfinite(p_in) & (p_in > 0)
    ss, pp = s[ok], np.log(p_in[ok])
    best = 0.0
    for i in range(len(ss) - 4):
        if ss[i + 4] > ss[i]:
            slope = np.polyfit(ss[i:i + 5], pp[i:i + 5], 1)[0]
            best = max(best, slope)
    # A length past the line's own is no exponential regime, as on a detuned line.
    if best <= 0 or 1.0 / best > float(ss[-1] - ss[0]):
        return None
    return 1.0 / best


def rows_for(deck: pathlib.Path, log: pathlib.Path):
    root = out_root(deck)
    stats = deck.parent / f"{root}.stats.h5"
    if not stats.exists():
        raise SystemExit(f"{deck}: no stats file at {stats}")
    rows = []
    with h5py.File(stats, "r") as h5:
        s = h5["coords/s"][()]
        power = h5["field/total/power"][()]
        energy = h5["field/total/energy"][()]
        inside = h5["field/total/power_inside_angle"][()]
        win = power.sum(axis=1)
        # The split is not computed at every record or every slice, and a slice without it
        # holds NaN. A record's inside power is the sum over the slices that have it, and a
        # record with none is NaN.
        have = np.isfinite(inside)
        win_in = np.where(have.any(axis=1), np.nansum(np.where(have, inside, 0.0), axis=1), np.nan)
        rows += [("exit_power", win[-1], "W"), ("exit_pulse_energy", energy[-1].sum(), "J")]
        n_slice = power.shape[1]
        rows.append(("exit_power_per_slice", win[-1] / n_slice, "W"))
        if n_slice > 1 and win[-1] > 0:
            rows.append(("exit_power_slice_spread", float(np.std(power[-1]) / np.mean(power[-1])), ""))
        k = int(np.argmax(win))
        rows += [("peak_power", win[k], "W"), ("peak_power_s", s[k], "m")]
        last_in = np.nonzero(np.isfinite(win_in))[0]
        if last_in.size:
            j = int(last_in[-1])
            rows += [("exit_power_inside", win_in[j], "W"), ("exit_power_outside", win[j] - win_in[j], "W")]
            if j != len(s) - 1:
                rows.append(("exit_power_inside_s", s[j], "m"))
        run = h5["run"]
        for name, key, unit in (("wide_angle_ratio", "wide_angle_ratio", ""),
                                ("wide_angle_ratio_exit", "wide_angle_ratio_last", ""),
                                ("mode_peak_inside", "power_inside_angle_peak", "W"),
                                ("mode_peak_s", "z_mode_peak", "m"),
                                ("bunching_peak", "bunching_peak", "")):
            if key in run:
                rows.append((name, float(run[key][()]), unit))
        rows.append(("split_angle", float(h5["field/total/split_angle"][()]), "rad"))
        b = h5["beam/slice/bunching"][()]
        bm = np.mean(np.abs(b), axis=1)
        kb = int(np.argmax(bm))
        rows += [("initial_bunching", float(bm[0]), ""), ("exit_bunching", float(bm[-1]), ""),
                 ("exit_bunching_max", float(np.max(np.abs(b[-1]))), ""),
                 ("bunching_max", float(bm[kb]), ""), ("bunching_max_s", s[kb], "m"),
                 ("exit_power_slice_max", float(np.max(power[-1])), "W")]
        e = h5["beam/slice/energy"][()]
        e0, e1 = float(np.nanmean(e[0])), float(np.nanmean(e[-1]))
        rows += [("initial_energy", e0, "eV"), ("exit_energy", e1, "eV"), ("energy_change", e1 - e0, "eV"),
                 ("exit_energy_span", float(np.nanmax(e[-1]) - np.nanmin(e[-1])), "eV")]
        lg = gain_length(s, win_in)
        if lg is not None:
            rows.append(("gain_length", lg, "m"))
        at_end = h5["coords/at_element_end"][()].astype(bool)
        first = np.nonzero(at_end & (s > 0))[0]
        if first.size:
            rows += [("power_after_first_element", win[first[0]], "W"),
                     ("power_after_first_element_per_slice", win[first[0]] / n_slice, "W")]
        se = h5["beam/slice/sigma_energy"][()]
        rows += [("initial_sigma_energy", float(np.mean(se[0])), "eV"), ("exit_sigma_energy", float(np.mean(se[-1])), "eV")]
        rows += [("n_slice", int(run["n_slice"][()]), ""), ("n_record", int(run["n_record"][()]), ""),
                 ("slice_spacing", float(run["slice_spacing"][()]), "m")]
        npart = h5["beam/slice/n_particle_live"][()]
        rows.append(("n_particle", int(npart[0].max()), ""))
        if "params/wavefront_init/grid_n_pts" in h5:
            rows.append(("grid_n_pts", int(h5["params/wavefront_init/grid_n_pts"][()]), ""))
    rows += from_log(log)
    return rows


def sample(log: pathlib.Path, out: pathlib.Path) -> None:
    lines = [l.rstrip() for l in log.read_text().splitlines() if l.startswith(SAMPLE_PREFIXES)]
    if not lines:
        raise SystemExit(f"{log}: no lines to sample")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("<!-- Generated by examples/run_examples.sh from the steady-state example's log. Do not edit. -->\n\n"
                   "```\n" + "\n".join(lines) + "\n```\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("deck", nargs="?", help="The deck the run read")
    ap.add_argument("log", nargs="?", help="The log the runner kept")
    ap.add_argument("--sample", metavar="LOG", help="Write the sample output of this log instead")
    ap.add_argument("--out", help="Where --sample writes")
    a = ap.parse_args()
    if a.sample:
        if not a.out:
            ap.error("--sample needs --out")
        sample(pathlib.Path(a.sample), pathlib.Path(a.out))
        return 0
    if not (a.deck and a.log):
        ap.error("give <deck> <log>, or --sample <log> --out <md>")
    deck, log = pathlib.Path(a.deck), pathlib.Path(a.log)
    for name, value, unit in rows_for(deck, log):
        print(f"{deck.parent.name}\t{deck.name}\t{name}\t{fmt(value)}\t{unit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
