#!/usr/bin/env python3
"""
Generate docs/generated/startup-noise-measured.md from the JSON startup_noise.py writes.

The measured tables of docs/startup-noise.md are written here from the record of the
runs that measured them, one JSON per machine under docs/generated/startup-noise/, and
the page includes each table from this file between its markers. A table is a MyST
table directive with its caption and its label, so every link to the page resolves as
it did when the tables were typed by hand.

The output is deterministic: the same JSON gives the same bytes, and the keystone
regenerates this file and requires no diff. The JSON is committed and the driver runs
on demand, so the regeneration costs no run. Each JSON carries the provenance of its
run, the Lucifer commit, the build, the genesis4 version, the backend and the day, and
the file states it in a table of its own.

Usage:

  report_startup_noise.py --json <dir> --out <file.md>

Each table block in the output is

  <!-- begin <label> -->
  ```{table} <caption>
  :name: <label>
  ...
  ```
  <!-- end <label> -->

and the page includes it with the include directive's start-after and end-before.
"""

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import startup_noise as sn  # noqa: E402

# The machines in the order the page presents them, with the name the page uses.
MACHINES = (("aramis", "Aramis"), ("flash", "FLASH probe"), ("flash1", "FLASH1"), ("lcls", "LCLS"))

# The sigmoid's width where the deck sets none, global%source_filter_width's default
# (docs/input-reference.md). The derived default rows of the filter experiment run at it.
DEFAULT_WIDTH = 0.05


def sig(x, n):
    """x to n significant figures, written positionally: 0.37, 1.7, 117, 0.047."""
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    r = round(x, n - 1 - e)
    if e >= n - 1:
        return f"{int(round(r)):d}"
    return f"{r:.{n - 1 - e}f}"


def tex(x, n):
    """x to n significant figures in TeX: $5.1\\times10^{-4}$, or positionally where
    the exponent is small."""
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    if -2 <= e <= 3:
        return sig(x, n)
    m = x / 10 ** e
    return f"${sig(m, n)}\\times10^{{{e}}}$"


PREFIXES = [(1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "u"), (1e-9, "n")]


def si(x, unit, n=3):
    """x with an SI prefix so the mantissa sits in [1, 1000): 625 MW, 6.59 kW, 974 nJ."""
    if x == 0:
        return f"0 {unit}"
    for scale, p in PREFIXES:
        if abs(x) >= scale:
            return f"{g(x / scale, n)} {p}{unit}"
    scale, p = PREFIXES[-1]
    return f"{g(x / scale, n)} {p}{unit}"


def pct(x):
    """A relative change as a signed percentage, whole above ten percent and to a tenth
    below: +102%, -40%, -5.9%, 0%."""
    v = 100 * x
    s = f"{v:+.0f}%" if abs(v) >= 10 else f"{v:+.1f}%"
    return "0%" if s in ("+0%", "-0%", "+0.0%", "-0.0%") else s


def fac(x):
    """A factor: three figures from a hundred up and two below, 1550x, 123x, 12x, 8.7x."""
    return sig(x, 3 if x >= 100 else 2) + "x"


def g(x, n=3):
    """x to n significant figures with trailing zeros dropped: 0.05, 1, 3, 13.7."""
    return f"{x:.{n}g}"


def load_records(directory):
    """{machine: record} for the machines whose JSON is present."""
    out = {}
    for key, _ in MACHINES:
        suffix = "" if key == "aramis" else f"-{key}"
        path = Path(directory) / f"startup-noise{suffix}.json"
        if path.exists():
            out[key] = json.loads(path.read_text())
    return out


def machine(key):
    """The machine's constants, with the driver's module pointed at it."""
    sn.select_machine(key)
    return sn.MACHINES[key]


def angles(key, th):
    """(the mode angle, the rho angle, the derived edge) in radians, as the run derives
    them: four diffraction angles of the mode, the angle where the resonant wavelength
    red-shifts by rho, and the larger of the two."""
    m = sn.MACHINES[key]
    theta_mode = 4 * th["theta_mode"]
    theta_rho = math.sqrt(2 * th["rho"] * m["lambda0"] / m["lambda_u"])
    return theta_mode, theta_rho, max(theta_mode, theta_rho)


def late_index(key):
    """The record index of the second-to-last dump, the one the page reads the bunching
    at, and the undulator it ends."""
    ends = sn.MACHINES[key]["und_end_record"]
    und = sorted(ends)[-2]
    return ends[und], und


def block(label, caption, header, rows):
    L = [f"<!-- begin {label} -->", f"```{{table}} {caption}", f":name: {label}", ""]
    L.append("| " + " | ".join(header) + " |")
    L.append("|" + "---|" * len(header))
    for r in rows:
        L.append("| " + " | ".join(r) + " |")
    L += ["```", f"<!-- end {label} -->", ""]
    return L


def table_provenance(recs):
    rows = []
    for key, name in MACHINES:
        if key not in recs or "provenance" not in recs[key]:
            continue
        p = recs[key]["provenance"]
        walls = recs[key].get("wall_s", {})
        total = sum(v for v in walls.values() if isinstance(v, (int, float)))
        rows.append([name, f"`{p['lucifer_commit']}`", p["build"], p["genesis4"], p["backend"],
                     p["date"], f"{total / 60:.0f} min" if total else ""])
    if not rows:
        return []
    return block("tab-sn-provenance",
                 "The runs that measured the tables below, one JSON per machine, as the driver "
                 "recorded them. A commit marked + was built from a tree with uncommitted changes. "
                 "The wall time is the sum over the experiments of the run that measured each one.",
                 ["machine", "Lucifer commit", "build", "genesis4", "backend", "date", "wall time"],
                 rows)


def table_derived(rec):
    m = machine("aramis")
    th = rec["theory"]
    lam, lam_u = m["lambda0"], m["lambda_u"]
    theta_mode, theta_rho, _ = angles("aramis", th)
    n_und = max(m["und_end_record"])
    periods = round(m["seg_length"] / lam_u)
    n_slice = th["N_lambda"] * sn.LONG[1]
    rows = [
        ["Pierce parameter", r"$\rho$", tex(th["rho"], 2), ""],
        ["cooperation length", r"$\lambda/(4\pi\rho)$", sig(lam / (4 * math.pi * th["rho"]) * 1e9, 2), "nm"],
        ["slippage over the line", "", sig(n_und * periods * lam * 1e9, 2), "nm"],
        ["central cone of one segment", r"$\theta_{con}$", sig(th["theta_cone_one_segment"] * 1e6, 2), "µrad"],
        [r"angle at which the resonant wavelength shifts by $\rho$", "", sig(theta_rho * 1e6, 2), "µrad"],
        ["diffraction angle of the mode", r"$\lambda/(2\pi\sigma_x)$", sig(th["theta_mode"] * 1e6, 2), "µrad"],
        [f"electrons per slice of {sn.LONG[1]} wavelengths", r"$N_\lambda$", tex(n_slice, 2), ""],
        ["rms shot-noise bunching of such a slice", r"$1/\sqrt{N_\lambda}$", tex(1 / math.sqrt(n_slice), 2), ""],
    ]
    return block("tab-sn-derived", "Derived quantities, from the deck's parameters and Ming Xie's fit.",
                 ["quantity", "symbol", "value", "unit"], rows)


def table_floor(rec):
    machine("aramis")
    fl = rec["floor"]
    grids = sorted({r["ngrid"] for r in fl.values() if isinstance(r, dict) and "ngrid" in r})
    loads = sorted({r["npart"] for r in fl.values() if isinstance(r, dict) and "npart" in r})
    z = fl[f"floor_g{grids[0]}_n{loads[0]}"]["dumps"][0]["z"]
    cut = sn.MACHINES["aramis"]["cuts"][1]
    rows = []
    for ng in grids:
        r = [sig(fl[f"floor_g{ng}_n{loads[0]}"]["dx"] * 1e6, 3)]
        for n in loads:
            k = f"floor_g{ng}_n{n}"
            r.append(sig(fl[k]["dumps"][0][sn.cut_key("out")] / 1e6, 2) if k in fl else "")
        rows.append(r)
    header = ["cell size (µm)"] + [f"{loads[0]} particles"] + [str(n) for n in loads[1:]]
    return block("tab-sn-floor",
                 f"Power per slice outside {g(cut * 1e6)} µrad at the end of the first "
                 f"undulator, z = {sig(z, 3)} m, in MW.", header, rows)


def table_beamlets(rec):
    machine("aramis")
    bl = rec["beamlets"]
    cut = sn.MACHINES["aramis"]["cuts"][1]
    keys = [k for k in bl if k.startswith("beamlet_")]
    any_r = bl[keys[0]]
    z37, z56 = any_r["dumps"][3]["z"], any_r["dumps"][4]["z"]
    dx = rec["floor"][f"floor_g{sn.MACHINES['aramis']['filter_grid']}_n1024"]["dx"]
    rows = [[str(bl[k]["npart"]), str(bl[k]["beamlet"]), str(bl[k]["beamlets"]),
             sig(bl[k]["dumps"][3][sn.cut_key("out")] / 1e6, 2),
             sig(bl[k]["dumps"][4][sn.cut_key("out")] / 1e6, 2)] for k in keys]
    return block("tab-sn-beamlets",
                 f"Power per slice outside {g(cut * 1e6)} µrad with the beamlet size and the "
                 f"particle count varied separately, {sig(dx * 1e6, 3)} µm cells, in MW.",
                 ["particles per slice", "beamlet size", "beamlets", f"z = {g(z37)} m", f"z = {g(z56)} m"],
                 rows)


def table_window(rec):
    ln = rec["line"]
    keys = [k for k in ln if k.startswith("line_")]
    grids = sorted({ln[k]["ngrid"] for k in keys})
    loads = sorted({ln[k]["npart"] for k in keys})
    rows = []
    for ng in grids:
        r = [sig(ln[f"line_g{ng}_n{loads[0]}"]["dx"] * 1e6, 3)]
        for n in loads:
            r.append(sig(ln[f"line_g{ng}_n{n}"]["P_exit"] / 1e9, 2))
        rows.append(r)
    return block("tab-sn-window",
                 "Total window power at the exit of the line for a window of 96 slices at three "
                 "wavelengths, in GW.",
                 ["cell size (µm)"] + [f"{n} particles" for n in loads], rows)


def table_genesis(rec):
    m = machine("aramis")
    ge = rec["genesis"]
    grids = sorted({ge[k]["ngrid"] for k in ge})
    rows = []
    for ng in grids:
        fft, adi = ge.get(f"G{ng}fft"), ge.get(f"G{ng}adi")
        dx = 2 * m["half_width"] / (ng - 1)
        rows.append([str(ng), sig(dx * 1e6, 3),
                     sig(fft["genesis_exit"] / 1e6, 3) if fft else "",
                     sig(adi["genesis_exit"] / 1e6, 3) if adi else "",
                     tex(fft["exit_rel"], 2) if fft and "exit_rel" in fft else ""])
    return block("tab-sn-genesis",
                 "Exit window power of Genesis4 4.6.15 for an unseeded simulation of a 32-slice "
                 "window on the line of the table above, with its FFT solver and its "
                 "alternating-direction implicit (ADI) solver, in MW, and the relative difference "
                 "of Lucifer's exit power from the FFT solver's, started from the same dumps.",
                 ["grid points", "cell size (µm)", "FFT solver", "ADI solver", "Lucifer against FFT"],
                 rows)


def filter_rows(key, rec):
    """The filter experiment's cases in the page's order, per load."""
    f = rec["filter"]
    loads = sorted({f[k]["npart"] for k in f if k.startswith("n") and isinstance(f[k], dict)})
    order = ("off", "default", "cone", "cone_sharp", "cut", "cut_sharp")
    return loads, [(n, c, f[f"n{n}_{c}"]) for n in loads for c in order if f"n{n}_{c}" in f]


def edge_label(key, th, f, case):
    _, _, derived = angles(key, th)
    if case == "off":
        return "none", ""
    if case == "default":
        return f"derived, {sig(derived * 1e6, 3)} µrad", g(DEFAULT_WIDTH)
    edge = f["edges"]["cone" if case.startswith("cone") else "cut"]["theta"]
    width = f["sharp_width"] if case.endswith("_sharp") else 1
    return f"{g(edge * 1e6, 2)} µrad", g(width)


def b_rel_late(r):
    return r.get("b_rel_late", r.get("b_rel_z37"))


def table_filter(rec):
    key = "aramis"
    machine(key)
    f, th = rec["filter"], rec["theory"]
    loads, cases = filter_rows(key, rec)
    _, und_late = late_index(key)
    z_late = cases[0][2]["z"][late_index(key)[0]]
    z0 = cases[0][2]["dumps"][0]["z"]
    rows = []
    for n, case, r in cases:
        edge, width = edge_label(key, th, f, case)
        off = case == "off"
        rows.append([str(n), edge, width,
                     "1" if off else fac(r["wide_factor_exit"]),
                     sig(r["dumps"][0][sn.cut_key("in")] / 1e6, 3),
                     "" if off else pct(r["mode_rel_exit"]),
                     "" if off else pct(b_rel_late(r)),
                     sig(r["z_sat"], 3)])
    cut = sn.MACHINES[key]["cuts"][1]
    return block("tab-sn-filter",
                 f"The source filter against the split, {sig(f['dx'] * 1e6, 3)} µm cells, "
                 f"{sn.LONG[0]}-slice window. The wide-angle column is the factor by which the "
                 f"power outside {g(cut * 1e6)} µrad falls at the exit. The in-cone column is "
                 f"the power inside {g(cut * 1e6)} µrad at the end of the first undulator, "
                 f"where the physical spontaneous emission within that angle is "
                 f"{si(th['W_incoh_in_cut_one_segment'], 'W')}. The mode column is the change in "
                 f"the power inside {g(cut * 1e6)} µrad at the exit, against the unfiltered "
                 f"run at the same load, and the bunching column the change in the mean bunching "
                 f"at the end of undulator {und_late}, z = {g(z_late)} m. The derived default's "
                 f"width is the program's default.",
                 ["particles", "edge", "width", "wide-angle at exit", f"in-cone at z = {g(z0)} m (MW)",
                  "mode at exit", f"bunching at {g(z_late)} m", "saturation (m)"], rows)


def table_filter_machine(key, name, rec, label, ordinal):
    machine(key)
    f, th = rec["filter"], rec["theory"]
    _, _, derived = angles(key, th)
    i_late, und_late = late_index(key)
    rows = []
    for case in ("off", "default"):
        r = f[f"n1024_{case}"]
        z_late = r["z"][i_late]
        rows.append(["none" if case == "off" else
                     f"derived default, {sig(derived * 1e6, 3)} µrad, width {g(DEFAULT_WIDTH)}",
                     sig(r["dumps"][0][sn.cut_key("in")] / 1e3, 3),
                     sig(r["dumps"][4][sn.cut_key("in")] / 1e9, 3),
                     si(r["dumps"][4][sn.cut_key("out")], "W"),
                     sig(r["b"][i_late], 3),
                     sig(r["z_sat"], 3)])
    z0 = f["n1024_off"]["dumps"][0]["z"]
    z_late = f["n1024_off"]["z"][i_late]
    return block(label,
                 f"The derived default on the {ordinal} machine. The split is at the derived edge, "
                 f"{sig(derived * 1e6, 3)} µrad, 1024 macroparticles per slice at "
                 f"{sig(f['dx'] * 1e6, 3)} µm cells, and the bunching is read at the end of "
                 f"undulator {und_late}, z = {g(z_late)} m.",
                 ["filter", f"in-cone at z = {g(z0)} m (kW)", "inside at exit (GW)",
                  "outside at exit", f"bunching at {g(z_late)} m", "saturation (m)"], rows)


def table_four(recs):
    rows = []
    for key, name in MACHINES:
        if key not in recs or "filter" not in recs[key]:
            continue
        m = machine(key)
        th, f = recs[key]["theory"], recs[key]["filter"]
        theta_mode, theta_rho, _ = angles(key, th)
        z_r = 4 * math.pi * th["sigma_x"] ** 2 / m["lambda0"]
        loads, _ = filter_rows(key, recs[key])
        for i, n in enumerate(loads):
            off, dflt = f[f"n{n}_off"], f[f"n{n}_default"]
            first = [name, si(m["lambda0"], "m"), sig(th["sigma_x"] * 1e6, 3),
                     sig(z_r / th["Lg"], 3), f"{sig(theta_mode * 1e6, 3)} µrad",
                     f"{sig(theta_rho * 1e6, 3)} µrad", sig(theta_mode / theta_rho, 3)] if i == 0 \
                else [name, "", "", "", "", "", ""]
            rows.append(first + [str(n), fac(dflt["wide_factor_exit"]), pct(dflt["mode_rel_exit"]),
                                 pct(b_rel_late(dflt)), f"{sig(off['z_sat'], 3)} to {sig(dflt['z_sat'], 3)}"])
    return block("tab-sn-filter-four",
                 "The derived default on the four machines. The split is at each machine's derived "
                 "edge. The wide-angle column is the factor by which the power outside that edge falls "
                 "at the exit, the mode column the change in the power inside it at the exit, and the "
                 "bunching column the change in the mean bunching at the end of the second-to-last "
                 "dumped undulator, each against the unfiltered run at the same load. The saturation "
                 "column is the unfiltered run's to the filtered run's.",
                 ["machine", "λ", "σₓ (µm)", "z_R/L_g", "mode angle", "ρ angle", "ratio", "load",
                  "wide-angle", "mode", "bunching", "saturation (m)"], rows)


def table_startup_four(recs):
    rows = []
    order = sorted((k for k, _ in MACHINES if k in recs and "filter" in recs[k]),
                   key=lambda k: sn.MACHINES[k]["seg_length"] / recs[k]["theory"]["Lg"])
    names = dict(MACHINES)
    for key in order:
        m = machine(key)
        th, f = recs[key]["theory"], recs[key]["filter"]
        loads, _ = filter_rows(key, recs[key])
        meas = f[f"n{loads[-1]}_off"]["dumps"][0][sn.cut_key("in")]
        ana = th["W_incoh_in_cut_one_segment"]
        rows.append([names[key], sig(m["seg_length"] / th["Lg"], 3), si(meas, "W"), si(ana, "W"),
                     sig(meas / ana, 3)])
    return block("tab-sn-startup-four",
                 "The power inside the reported cut at the end of the first segment, against eq. (1) "
                 "of Saldin, Schneidmiller and Yurkov for the same cut, at the most converged load of "
                 "each machine.",
                 ["machine", "first segment in power gain lengths", "measured", "analytic", "ratio"], rows)


def table_tolerance(recs):
    rows = []
    grids = []
    for key, name in MACHINES:
        if key not in recs or "tolerance" not in recs[key]:
            continue
        t = recs[key]["tolerance"]
        grids.append(f"{t['ngrid']} on {name}")
        base = t["rows"]["0.75"]
        for tol in t["tolerances"]:
            r = t["rows"][str(tol)]
            if tol == 0.75:
                rows.append([name, g(tol), si(r["accepted_first_segment_W"], "W"),
                             tex(r["wide_over_mode_exit"], 3), f"{sig(r['Lg_fit'], 3)} m",
                             f"{sig(r['z_sat'], 3)} m", si(r["pulse_energy_exit_J"], "J")])
            else:
                rows.append(["", g(tol),
                             pct(r["accepted_first_segment_W"] / base["accepted_first_segment_W"] - 1),
                             "x" + sig(r["wide_over_mode_exit"] / base["wide_over_mode_exit"], 2),
                             pct(r["Lg_fit"] / base["Lg_fit"] - 1), f"{sig(r['z_sat'], 3)} m",
                             pct(r["pulse_energy_exit_J"] / base["pulse_energy_exit_J"] - 1)])
    if not rows:
        return []
    return block("tab-sn-tolerance",
                 "The passband tolerance on the four machines, grid " + ", ".join(grids) +
                 ", 1024 macroparticles per slice, one seed a row, each machine's rows against its "
                 "own 0.75 row. The accepted power is the power inside the protected angle at the "
                 "end of the first segment, and the gain length is the fit through the exponential "
                 "regime.",
                 ["machine", "tolerance", "accepted at the first segment", "wide/mode at the exit",
                  "gain length", "saturation", "pulse energy"], rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", required=True, help="The directory holding startup-noise*.json")
    p.add_argument("--out", required=True, help="The Markdown file to write")
    a = p.parse_args()
    recs = load_records(a.json)
    if "aramis" not in recs:
        raise SystemExit(f"{a.json}: no startup-noise.json")
    ar = recs["aramis"]
    L = ["<!-- Generated by tests/scripts/report_startup_noise.py from the JSON startup_noise.py "
         "writes. Do not edit. -->", ""]
    L += table_provenance(recs)
    L += table_derived(ar)
    missing = []
    for key, fn in (("floor", table_floor), ("beamlets", table_beamlets), ("line", table_window),
                    ("genesis", table_genesis), ("filter", table_filter)):
        if key in ar:
            L += fn(ar)
        else:
            missing.append(key)
    for key, name, label, ordinal in (("flash", "FLASH probe", "tab-sn-filter-flash", "second"),
                                      ("flash1", "FLASH1", "tab-sn-filter-flash1", "third")):
        if key in recs and "filter" in recs[key]:
            L += table_filter_machine(key, name, recs[key], label, ordinal)
        else:
            missing.append(f"filter on {key}")
    L += table_four(recs)
    L += table_startup_four(recs)
    tol = table_tolerance(recs)
    if tol:
        L += tol
    else:
        missing.append("tolerance")
    # Every block the page includes must be present. An include whose markers are missing
    # brings in the whole file, which duplicates every label on the page, so a record
    # lacking an experiment refuses here rather than renders.
    if missing or "provenance" not in ar:
        raise SystemExit(f"{a.json}: the record lacks "
                         + ", ".join(missing + ([] if "provenance" in ar else ["provenance"]))
                         + ", so the page's includes would not resolve. Run the driver.")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(L))
    n = sum(1 for x in L if x.startswith("<!-- begin"))
    print(f"wrote {a.out}: {n} tables from {len(recs)} machines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
