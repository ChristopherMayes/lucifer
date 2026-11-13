#!/usr/bin/env python3
"""
Record what an example's run measured, from its stats file.

  example_record.py <input.toml> <log>              # rows of the record, to stdout
  example_record.py --sample <log> [--sample-ascii <log>] --out <md>   # the sample output the user guide shows
  example_record.py --exit <input.toml> [--since <s>]   # the out_root, and the light at the exit

The first form reads the input file for its out_root, opens the stats file the run wrote
beside it and the log the runner kept, and prints one row per quantity:

  <directory> <tab> <input> <tab> <quantity> <tab> <value> <tab> <unit>

The quantities are one fixed set for every input file, so a page can compare its input files column
by column and a number a README once typed from a run comes from the run itself:

  exit_power, exit_window_energy         the window totals at the last record
  exit_pulse_energy                      the whole pulse's energy at the last record, the
                                         window's and what slippage carried out of it
                                         (field/total/pulse_energy)
  exit_power_<component>                 the exit power of each field component beyond the
                                         total and the x plane, such as y or h3
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
  slice_power_max_at_bunching_max, exit_power_slice_max
                                         the largest slice's power at the saturation point
                                         and at the exit
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
  beamlet_size, beamlet_size_derived, n_particle_derived
                                         the beamlet size, and whether setup derived it and
                                         the macroparticles in a slice
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

# What opens a rule of a box of the printout, in Unicode and in ASCII, and the tag line
# beside the logo, at whose column the logo's other lines carry the version, the date
# and the host.
RULE_STARTS, TAG_LINE = (" \u256d", " \u251c", " \u2570", " +"), "FEL tracking in Bmad"


def out_root(input_file: pathlib.Path) -> str:
    import tomllib
    root = tomllib.loads(input_file.read_text()).get("run", {}).get("out_root")
    if not root:
        raise SystemExit(f"{input_file}: no out_root")
    return root


def si(value: float, unit: str) -> str:
    """A number with an SI prefix, three decimals, as a table row reads best."""
    if value == 0:
        return f"0 {unit}"
    exp = max(-24, min(15, 3 * math.floor(math.log10(abs(value)) / 3)))
    return f"{value / 10.0 ** exp:.3f} {'yzafpnµm kMGTP'[exp // 3 + 8].strip()}{unit}"


def exit_line(input_file: pathlib.Path, since: float) -> str:
    """The out_root, and on a second line the light at the exit as the stats file's run/
    group holds it. Fails where the stats file is absent or older than ``since``."""
    root = out_root(input_file)
    stats = input_file.parent / f"{root}.stats.h5"
    if not stats.exists() or stats.stat().st_mtime < since:
        raise SystemExit(f"{root}\n")
    with h5py.File(stats, "r") as h5:
        run = h5["run"]
        return (f"{root}\npower {si(float(run['exit_power'][()]), 'W')}, window energy "
                f"{si(float(run['exit_window_energy'][()]), 'J')}, pulse energy "
                f"{si(float(run['exit_pulse_energy'][()]), 'J')}")


def fmt(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    f = float(v)
    if math.isnan(f) or math.isinf(f):
        return repr(f)
    return f"{f:.4e}"


def from_setup(h5):
    """The setup quantities, as (name, value, unit), from the stats file's setup/ group
    and the origin each key of params/ carries."""
    out = []
    setup = h5["setup"]

    def origin(key):
        return h5[f"params/{key}"].attrs["origin"].decode() if f"params/{key}" in h5 else ""

    out.append(("grid_cell", float(setup["grid_cell"][()]), "m"))
    out.append(("grid_derived", origin("field/grid_n_pts") == "derived", ""))
    if "grid_n_pts_unrounded" in setup and int(setup["grid_n_pts_unrounded"][()]) != int(setup["grid_n_pts"][()]):
        out.append(("grid_derived_from", int(setup["grid_n_pts_unrounded"][()]), ""))
    out += [("beamlet_size", int(setup["beamlet_size"][()]), ""),
            ("beamlet_size_derived", origin("beam/beamlet_size") == "derived", "")]
    out.append(("n_particle_derived", origin("beam/n_particle") == "derived", ""))
    out.append(("fel_steps", int(setup["n_fel_step"][()]), ""))
    if "filter_mode_angle" in setup:
        mode, rho = float(setup["filter_mode_angle"][()]), float(setup["filter_rho_angle"][()])
        out += [("filter_mode_angle", mode, "rad"), ("filter_rho_angle", rho, "rad"),
                ("filter_ratio", mode / rho if rho > 0 else 0.0, ""), ("filter_edge", float(setup["filter_edge"][()]), "rad")]
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


def rows_for(input_file: pathlib.Path, log: pathlib.Path):
    root = out_root(input_file)
    stats = input_file.parent / f"{root}.stats.h5"
    if not stats.exists():
        raise SystemExit(f"{input_file}: no stats file at {stats}")
    rows = []
    with h5py.File(stats, "r") as h5:
        s = h5["coords/s"][()]
        power = h5["field/total/power"][()]
        energy = h5["field/total/energy"][()]
        pulse = h5["field/total/pulse_energy"][()]
        inside = h5["field/total/power_inside_angle"][()]
        win = power.sum(axis=1)
        # The split is not computed at every record or every slice, and a slice without it
        # holds NaN. A record's inside power is the sum over the slices that have it, and a
        # record with none is NaN.
        have = np.isfinite(inside)
        win_in = np.where(have.any(axis=1), np.nansum(np.where(have, inside, 0.0), axis=1), np.nan)
        rows += [("exit_power", win[-1], "W"), ("exit_window_energy", energy[-1].sum(), "J"),
                 ("exit_pulse_energy", pulse[-1], "J")]
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
        # Each field component beyond the total and the x plane, a second polarization or
        # a harmonic, gets its own exit power, so a page can show what the component did.
        # A polarization plane carries its power directly. A harmonic is a field of its
        # own, with planes under it, and carries its power under its total.
        for comp in sorted(h5["field"]):
            if comp in ("total", "x"):
                continue
            for key in (f"field/{comp}/power", f"field/{comp}/total/power"):
                if key in h5:
                    rows.append((f"exit_power_{comp}", float(h5[key][-1].sum()), "W"))
                    break
        b = h5["beam/slice/bunching"][()]
        bm = np.mean(np.abs(b), axis=1)
        kb = int(np.argmax(bm))
        rows += [("initial_bunching", float(bm[0]), ""), ("exit_bunching", float(bm[-1]), ""),
                 ("exit_bunching_max", float(np.max(np.abs(b[-1]))), ""),
                 ("bunching_max", float(bm[kb]), ""), ("bunching_max_s", s[kb], "m"),
                 ("slice_power_max_at_bunching_max", float(np.max(power[kb])), "W"),
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
        if "params/field/grid_n_pts" in h5:
            rows.append(("grid_n_pts", int(h5["params/field/grid_n_pts"][()]), ""))
        rows += from_setup(h5)
    return rows


def sample_lines(log: pathlib.Path) -> list[str]:
    """The printout of a log as the user guide shows it: the logo, the boxes of the header
    with their legend and warnings, and the result. Left out are what changes from run to
    run or names the machine: the version, date and host beside the logo, the estimate,
    the progress rows, the timing and the sizes of the files."""
    lines = [l.rstrip() for l in log.read_text().splitlines()]
    rules = [i for i, l in enumerate(lines) if l.startswith(RULE_STARTS) and len(l) > 40]
    if not rules or TAG_LINE not in lines[rules[0] - 4]:
        raise SystemExit(f"{log}: no printout to sample")
    first = rules[0]
    column = lines[first - 4].index(TAG_LINE)
    out = [lines[first - 4]] + [l[:column].rstrip() for l in lines[first - 3:first - 1] if l[:column].strip()] + [""]
    i = first
    while i < len(lines) and lines[i].strip():
        if " estimate " not in lines[i][:14]:
            out.append(lines[i])
        i += 1
    while i < len(lines) and (not lines[i].strip() or lines[i].startswith(" ! ")):
        if lines[i].strip():
            out.append(lines[i])
        i += 1
    result = [i for i in rules if " result " in lines[i]]
    if result:
        out.append("")
        i = result[0]
        while " timing " not in lines[i]:
            out.append(lines[i])
            i += 1
        out.append(lines[rules[-1]])
    return out


def sample(log: pathlib.Path, out: pathlib.Path, ascii_log: pathlib.Path | None = None) -> None:
    text = ("<!-- Generated by examples/run_examples.sh from the steady-state example's log. Do not edit. -->\n\n"
            "```\n" + "\n".join(sample_lines(log)) + "\n```\n")
    if ascii_log is not None:
        text += ("\nWith `--ascii`, and `run.load_only=true` on the command line, the same input file "
                 "prints its header and stops:\n\n```\n" + "\n".join(sample_lines(ascii_log)) + "\n```\n")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("input_file", nargs="?", help="The input file the run read")
    ap.add_argument("log", nargs="?", help="The log the runner kept")
    ap.add_argument("--sample", metavar="LOG", help="Write the sample output of this log instead")
    ap.add_argument("--sample-ascii", metavar="LOG", help="With --sample, the log of the same input under --ascii")
    ap.add_argument("--out", help="Where --sample writes")
    ap.add_argument("--exit", metavar="INPUT", help="Print this input file's out_root and its run's light at the exit")
    ap.add_argument("--since", type=float, default=0.0, help="With --exit, the time the run started (s)")
    a = ap.parse_args()
    if a.exit:
        try:
            print(exit_line(pathlib.Path(a.exit), a.since))
        except SystemExit as err:
            print(str(err.code).strip())
            return 1
        return 0
    if a.sample:
        if not a.out:
            ap.error("--sample needs --out")
        sample(pathlib.Path(a.sample), pathlib.Path(a.out),
               pathlib.Path(a.sample_ascii) if a.sample_ascii else None)
        return 0
    if not (a.input_file and a.log):
        ap.error("give <input> <log>, or --sample <log> --out <md>")
    input_file, log = pathlib.Path(a.input_file), pathlib.Path(a.log)
    for name, value, unit in rows_for(input_file, log):
        print(f"{input_file.parent.name}\t{input_file.name}\t{name}\t{fmt(value)}\t{unit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
