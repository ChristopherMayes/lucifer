#!/usr/bin/env python3
"""
Generate doc/generated/layouts.md, the HDF5 hierarchy of the files this program writes.

The page is an inventory and not a description. The generator runs a few small decks
with the executable it is handed, opens the files they wrote, and prints for every
group, dataset, link and attribute what h5py reports: the hierarchy with its empty
groups, a dataset's shape, on-disk datatype and storage layout, every attribute's name,
datatype, shape and value, and a link's target, with external links left unfollowed.
No dataset values are printed. Nothing is classified or filtered. The only values left
out are the attributes named in OMITTED, which change from run to run or build to build,
and each of those keeps its name, datatype and shape beside a marker. A floating value
prints at six significant figures with no cutoff, so a small nonzero value, a sign or an
imaginary part stays visible.

Deterministic by construction: fixed seeds, one thread, names in sorted order, no
timestamp and no path outside the tree on the page. The keystone regenerates the page and
requires an empty diff, and a moved line is a layout change to review and regenerate
(doc/validation.md).

The page is written only after every fixture ran and every listed file was opened, and
it is written whole, so a failed fixture leaves the committed page as it was. The log
records the executable and the version it reported, and every run's output.

Usage:

  report_layouts.py --exe <lucifer> --out <file.md> [--log <file>]
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

# Attribute names whose values are omitted wherever they occur. Names, not patterns.
OMITTED = ("date", "softwareVersion")

# The lattice head every fixture shares: the benchmark's beam energy and twiss.
LATTICE_HEAD = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 8.53711
beginning[alpha_a] = -0.703306
beginning[beta_b] = 17.3899
beginning[alpha_b] = 1.40348
"""

# Four periods of the benchmark undulator, aw = 0.84853, one comb step per period.
UNDULATOR = ("QU: wiggler, l = 0.06, l_period = 0.015, field_calc = planar_model, "
             "b_max = sqrt(2) * 0.84853 * (twopi / 0.015) * m_electron / c_light, "
             "fel_method = {method}, ds_step = 0.015{extra}")

# The beam and field every fixture shares: 64 particles in beamlets of 8, no seed, no
# shot noise, a 16-point grid, a window of two wavelengths at 200 pm slice spacing.
NAMELIST = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%source_filter = F
  global%ran_seed = 4242
  slicing%n_wavelength = 2
  slicing%window_length = 1.2e-9
{params}/
&fel_beam_init
  beam_init%n_particle = 64
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
{beam}/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 0
  wavefront_init%grid_n_pts = 16
  wavefront_init%grid_half_width = 2e-4
/
"""

GRID_Z = """  beam_init%sig_z = 0
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -2.5e-10
  beam_init%grid(3)%x_max = 2.5e-10
"""

# Each fixture is one run: its root, its lattice body, its deck additions, and the files
# it contributes to the page, each with the case it shows and the file's role.
FIXTURES = (
    dict(root="avg",
         lattice=UNDULATOR.format(method="averaged", extra="") + "\nQL: line = (QU)\nuse, QL\n",
         params="  global%dump_at_comb = T\n  global%dump_orbit = T\n",
         beam=GRID_Z,
         settings="averaged QU, `dump_at_comb = T`, `dump_orbit = T`, the whole window, z on a grid",
         files=(
             ("avg-000002.gc.h5", "interior frame of an averaged undulator, the coordinates the map evolves", "diagnostic"),
             ("avg-000002.beam.h5", "the same record exported as the orbit", "exchange"),
             ("avg-000005.beam.h5", "frame at the end of QU, the whole window, at the same s as the final dump", "exchange, checkpoint"),
             ("avg-final.beam.h5", "final dump at the zero-length END marker Bmad appends to the line, the beam of the checkpoint pair", "exchange, checkpoint"),
             ("avg-final.wf.h5", "final dump at the zero-length END marker Bmad appends to the line, the field of the checkpoint pair", "exchange, checkpoint"),
         )),
    dict(root="unavg",
         lattice=UNDULATOR.format(method="unaveraged", extra="") + "\nQL: line = (QU)\nuse, QL\n",
         params="  global%dump_at_comb = T\n  global%unaveraged_steps_per_period = 40\n  global%unaveraged_ramp_periods = 2\n",
         beam=GRID_Z,
         settings="unaveraged QU, 40 steps per period, a two-period ramp, `dump_at_comb = T`",
         files=(
             ("unavg-000002.beam.h5", "interior frame of an unaveraged undulator", "exchange"),
         )),
    dict(root="crop",
         lattice=UNDULATOR.format(method="averaged", extra=", tilt = 0.4") + "\nQL: line = (QU)\nuse, QL\n",
         params="  global%dump_at_comb = T\n  global%dump_slice_first = 2\n  global%dump_slice_last = 4\n",
         beam=GRID_Z,
         settings="averaged QU with `tilt = 0.4`, `dump_at_comb = T`, `dump_slice_first = 2`, `dump_slice_last = 4`",
         files=(
             ("crop-000002.wf.h5", "field frame of both polarizations over three of six slices", "diagnostic"),
         )),
    dict(root="keep",
         lattice=UNDULATOR.format(method="averaged", extra="") + "\nQL: line = (QU)\nuse, QL\n",
         params="",
         beam="  load_mode = \"keep\"\n  beam_init%sig_z = 1.0e-10\n",
         settings="averaged QU, `load_mode = \"keep\"`, `sig_z = 1.0e-10` in a six-slice window",
         files=(
             ("keep-final.beam.h5", "final dump of several patches, two of them empty", "exchange, checkpoint"),
         )),
)


def anchor(file_name):
    return "layout-" + re.sub(r"[^a-z0-9]+", "-", file_name.lower()).strip("-")


# ---------------------------------------------------------------------------
# Datatypes, layouts and values, as h5py reports them.

def describe_dtype(dt):
    """The on-disk datatype in words: a string's length rule and encoding, a number's
    width and byte order, an enum's members, a compound's fields."""
    info = h5py.check_string_dtype(dt)
    if info is not None:
        length = "variable" if info.length is None else f"fixed({info.length})"
        return f"string {length} {info.encoding}"
    if h5py.check_ref_dtype(dt) is not None:
        return "reference"
    vlen = h5py.check_vlen_dtype(dt)
    if vlen is not None:
        return f"vlen({describe_dtype(np.dtype(vlen))})"
    enum = h5py.check_enum_dtype(dt)
    if enum is not None:
        members = ", ".join(f"{k}={v}" for k, v in sorted(enum.items(), key=lambda kv: kv[1]))
        return f"enum {dt.base.name} {{{members}}}"
    if dt.names:
        fields = ", ".join(f"{n}: {describe_dtype(dt[n])}" for n in dt.names)
        return f"compound{{{fields}}}"
    if dt.kind == "S":
        return f"string fixed({dt.itemsize}) bytes"
    order = {"<": " little-endian", ">": " big-endian"}.get(dt.byteorder, "")
    if dt.byteorder == "=" and dt.itemsize > 1:
        order = f" {sys.byteorder}-endian"
    return f"{dt.name}{order}"


def describe_layout(ds):
    if ds.chunks is None:
        parts = ["contiguous"]
    else:
        parts = [f"chunks {tuple(int(c) for c in ds.chunks)}"]
    if ds.compression is not None:
        opts = "" if ds.compression_opts is None else f"({ds.compression_opts})"
        parts.append(f"{ds.compression}{opts}")
    if ds.shuffle:
        parts.append("shuffle")
    if ds.fletcher32:
        parts.append("fletcher32")
    return ", ".join(parts)


def fmt_scalar(v):
    if isinstance(v, (bytes, np.bytes_)):
        return '"' + v.decode("ascii", "backslashreplace") + '"'
    if isinstance(v, str):
        return '"' + v + '"'
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (complex, np.complexfloating)):
        return f"{v.real:.6g}{v.imag:+.6g}j"
    if isinstance(v, (float, np.floating)):
        return f"{v:.6g}"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return repr(v)


def fmt_value(v):
    if isinstance(v, h5py.Empty):
        return "(empty)"
    if isinstance(v, np.ndarray):
        if v.shape == ():
            return fmt_scalar(v[()])
        return "[" + ", ".join(fmt_scalar(x) for x in v.ravel().tolist()) + "]"
    return fmt_scalar(v)


def shape_text(shape):
    if shape is None:
        return "null"
    return "(" + ", ".join(str(int(s)) for s in shape) + ("," if len(shape) == 1 else "") + ")"


# ---------------------------------------------------------------------------
# The traversal.

def attribute_lines(obj, depth):
    rows = []
    for name in sorted(obj.attrs):
        aid = obj.attrs.get_id(name)
        left = "  " * depth + "@" + name
        kind = f"{describe_dtype(aid.dtype)} {shape_text(aid.shape)}"
        if name in OMITTED:
            rows.append((left, kind, "(omitted)"))
        else:
            rows.append((left, kind, fmt_value(obj.attrs[name])))
    return rows


def walk(group, depth, rows, seen):
    for name in sorted(group.keys()):
        link = group.get(name, getlink=True)
        left = "  " * depth + name
        if isinstance(link, h5py.SoftLink):
            rows.append((left, "soft link", f"-> {link.path}"))
            continue
        if isinstance(link, h5py.ExternalLink):
            rows.append((left, "external link", f"-> {link.filename}:{link.path} (not followed)"))
            continue
        obj = group[name]
        if isinstance(obj, h5py.Dataset):
            rows.append((left, "dataset", f"{shape_text(obj.shape)} {describe_dtype(obj.dtype)}, {describe_layout(obj)}"))
            rows.extend(attribute_lines(obj, depth + 1))
        elif isinstance(obj, h5py.Group):
            addr = h5py.h5o.get_info(obj.id).addr
            if addr in seen:
                rows.append((left + "/", "group", f"hard link to {seen[addr]}"))
                continue
            seen[addr] = obj.name
            rows.append((left + "/", "group", ""))
            rows.extend(attribute_lines(obj, depth + 1))
            walk(obj, depth + 1, rows, seen)
        else:
            rows.append((left, type(obj).__name__.lower(), ""))


def inventory(path):
    """The file's hierarchy as aligned text lines."""
    rows = []
    with h5py.File(path, "r") as f:
        rows.append(("/", "group", ""))
        rows.extend(attribute_lines(f, 1))
        walk(f, 1, rows, {h5py.h5o.get_info(f.id).addr: "/"})
    w0 = max(len(r[0]) for r in rows)
    w1 = max(len(r[1]) for r in rows)
    return [f"{a.ljust(w0)}  {b.ljust(w1)}  {c}".rstrip() for a, b, c in rows]


# ---------------------------------------------------------------------------
# The fixtures.

def run_fixtures(exe, work, log):
    """Run every fixture in work. Returns the version line the program printed, or
    raises RuntimeError with the run that failed."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    version = ""
    for fx in FIXTURES:
        root = fx["root"]
        (work / f"{root}.bmad").write_text(LATTICE_HEAD + fx["lattice"])
        (work / f"{root}.nml").write_text(NAMELIST.format(lat=f"{root}.bmad", root=root,
                                                            params=fx["params"], beam=fx["beam"]))
        r = subprocess.run([str(exe), f"{root}.nml"], cwd=str(work), capture_output=True,
                           text=True, errors="replace", env=env)
        log.append(f"==== fixture {root}: exit {r.returncode}\n{r.stdout}\n{r.stderr}\n")
        if r.returncode != 0:
            raise RuntimeError(f"fixture {root} exited {r.returncode}")
        m = re.search(r"^.*Bmad version.*$", r.stdout, re.M)
        if m and not version:
            version = m.group(0).strip()
        for name, _, _ in fx["files"]:
            if not (work / name).is_file():
                raise RuntimeError(f"fixture {root} did not write {name}")
    return version


def software_version(work):
    """The softwareVersion the first beam file carries, for the log alone."""
    for fx in FIXTURES:
        for name, _, _ in fx["files"]:
            if name.endswith(".beam.h5"):
                with h5py.File(work / name, "r") as f:
                    v = f.attrs.get("softwareVersion")
                return fmt_value(v) if v is not None else "(none)"
    return "(none)"


# ---------------------------------------------------------------------------
# The page.

def page(work):
    L = []
    L.append("---")
    L.append("title: HDF5 file layouts")
    L.append("short_title: HDF5 layouts")
    L.append("---")
    L.append("")
    L.append("<!-- Generated by tests/scripts/report_layouts.py from files the executable wrote. Do not edit. -->")
    L.append("")
    L.append("# HDF5 file layouts")
    L.append("")
    L.append("The hierarchy of each representative file this program writes, printed from the file "
             "itself: every group, dataset, link and attribute in name order, and no dataset values. "
             "A dataset line gives the shape, the on-disk datatype and the storage layout. An attribute "
             "line, marked `@`, gives the datatype, the shape and the value, a floating value at six "
             "significant figures. A string datatype says whether its length is fixed or variable and "
             "names its encoding. A soft or external link shows its target and is not followed. Two "
             "attribute names are left out wherever they occur, `" + "` and `".join(OMITTED) + "`, "
             "since their values change between runs or builds, and each keeps its name, datatype and "
             "shape beside the marker `(omitted)`. Nothing else is left out. What the records mean is "
             "on [Reading the output](../reading-output.md).")
    L.append("")
    L.append("## Contents")
    L.append("")
    for fx in FIXTURES:
        for name, case, _ in fx["files"]:
            L.append(f"- [`{name}`](#{anchor(name)}): {case}")
    L.append("")
    L.append("## Fixtures")
    L.append("")
    L.append("Each file below comes from one of four runs of the executable the regeneration test "
             "hands this generator, in a scratch directory of its own. The runs share the benchmark "
             "beam energy and twiss, four periods of an undulator of `aw = 0.84853` with one comb step "
             "per period, 64 particles in beamlets of 8, no seed and no shot noise, a 16-point grid, "
             "and a window of two wavelengths at 1e-10 m, six slices of 200 pm. The settings column "
             "lists what a run adds to that.")
    L.append("")
    L.append("| Run | Settings | Files and their roles |")
    L.append("|---|---|---|")
    for fx in FIXTURES:
        files = ", ".join(f"[`{name}`](#{anchor(name)}) ({role})" for name, _, role in fx["files"])
        L.append(f"| `{fx['root']}` | {fx['settings']} | {files} |")
    L.append("")
    for fx in FIXTURES:
        for name, case, role in fx["files"]:
            L.append(f"({anchor(name)})=")
            L.append(f"## {name}")
            L.append("")
            L.append(f"Run `{fx['root']}`: {case}. Role: {role}.")
            L.append("")
            L.append("```text")
            L.extend(inventory(work / name))
            L.append("```")
            L.append("")
    return "\n".join(L).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="the lucifer executable to run the fixtures with")
    ap.add_argument("--out", required=True, help="the Markdown page to write")
    ap.add_argument("--log", help="where the runs' output goes (default: stderr on failure)")
    args = ap.parse_args()

    exe = Path(args.exe).resolve()
    out = Path(args.out)
    log = [f"executable: {exe}\n"]
    if not os.access(exe, os.X_OK):
        log.append("==== failed: not an executable\n")
        if args.log:
            Path(args.log).write_text("".join(log))
        print(f"report_layouts: not an executable: {exe}", file=sys.stderr)
        return 2

    status = 0
    with tempfile.TemporaryDirectory(prefix="lucifer-layouts-") as tmp:
        work = Path(tmp)
        try:
            version = run_fixtures(exe, work, log)
            log.insert(1, f"version line: {version}\nsoftwareVersion: {software_version(work)}\n")
            text = page(work)
        except (RuntimeError, OSError) as e:
            log.append(f"==== failed: {e}\n")
            status = 1
        else:
            # Whole or not at all: the committed page is never left partly rewritten.
            tmp_out = out.with_name(out.name + ".tmp")
            tmp_out.write_text(text)
            os.replace(tmp_out, out)
            log.append(f"==== wrote {out}\n")

    if args.log:
        Path(args.log).write_text("".join(log))
    if status:
        print("".join(log)[-4000:], file=sys.stderr)
        print("report_layouts: a fixture failed, the page was not written", file=sys.stderr)
    return status


if __name__ == "__main__":
    sys.exit(main())
