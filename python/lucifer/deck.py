"""Write the TOML deck the program reads.

A deck is one TOML file: the lattice at the top and the tables [run], [window], [beam]
with [beam.resample], [field], [wake], [space_charge] and [bmad]. The check scripts build
their configurations as flat ``key = value`` text, one key a line, in the form the old
namelist decks used, and this module routes each key into its table and renders the
value as TOML. A key may be written bare, ``n_particle = 2048``, and is routed by name,
or dotted, ``beam.n_particle = 2048``, which names its table outright. Values follow the
old form: ``T`` and ``F`` are booleans, a quoted string is a string, a number is a number,
and a comma-separated list is a list. The longitudinal fill, ``distribution_type(3)``
and ``grid(3)%x_min`` and ``grid(3)%x_max``, becomes ``z_distribution``, ``z_min`` and
``z_max``. An old group wrapper, ``&fel_params`` or ``/``, and a table heading, ``[run]``, are
ignored, so an old deck converts as it stands and a flat text may carry headings as
notes.

  from lucifer.deck import to_toml, write
  write("run.toml", "lattice = 'aramis.bmad'\\nn_particle = 2048\\nlambda0 = 1e-10\\n")

The routing table here names every key of the schema. ``lucifer --schema`` is the
source of truth, and tests/scripts/check_input_reference.py holds this table to it.
"""
from __future__ import annotations

import re

TABLES = ("run", "window", "beam", "beam.resample", "field", "wake", "space_charge", "bmad")

# Every key, by table. The order is the deck's.
KEYS = {
    "": ["lattice"],
    "run": [
        "out_root", "interlude_model", "transport_model", "unaveraged_steps_per_period",
        "unaveraged_ramp_periods", "source_model", "source_filter", "source_filter_angle",
        "source_filter_tolerance", "source_filter_xcut", "source_filter_ycut",
        "source_filter_width", "source_filter_mutate", "track_start", "track_end",
        "continuation", "comb_ds_save", "interlude_ds_step", "dump_beam_at", "dump_field_at",
        "dump_at_comb", "dump_slice_first", "dump_slice_last", "dump_reduced", "ran_seed",
        "write_diag", "write_initial", "load_only", "keep_escaped_field", "migrate",
        "migrate_check", "reference_run", "record_environment", "fp32_check", "fp32_mutate",
        "device", "device_timing", "device_dep_mutate",
    ],
    "window": ["window_length", "n_wavelength", "n_slice", "current"],
    "beam": [
        "file", "position_file", "n_particle", "bunch_charge", "sig_z", "sig_pz", "a_norm_emit",
        "b_norm_emit", "z_distribution", "z_min", "z_max", "load_mode", "quiet_start",
        "beamlet_size", "shot_noise", "write_genesis_dist", "write_openpmd_file",
        "split_weights", "swap_beam_xy", "gen_test_weights", "resample_split_weights",
    ],
    "beam.resample": ["slice_width", "n_particle_per_slice", "beamlet_size", "n_slice", "use_beam_init"],
    "field": ["file", "lambda0", "grid_n_pts", "grid_half_width", "seed_power", "seed_waist_size",
              "seed_polarization", "harmonics"],
    "wake": ["on", "model", "loss", "radius", "conductivity", "relaxation", "roundpipe", "material",
             "gap", "lgap", "hrough", "lrough", "write_kernels"],
    "space_charge": ["model", "rmax", "ngrid", "nz", "nphi", "longrange"],
    "bmad": ["radiation_damping_on", "radiation_fluctuations_on", "csr_and_space_charge_on",
             "sr_wakes_on", "lr_wakes_on", "default_ds_step", "space_charge_mesh_size",
             "sc_min_in_bin", "rel_tol_tracking", "particle_sigma_cutoff"],
}

# The old names, as the namelist decks and the check scripts' templates wrote them, to
# the dotted key. A prefix maps a whole structure; a full name maps one key.
OLD_PREFIX = {
    "global%": "run.", "slicing%": "window.", "wavefront_init%": "field.",
    "chamber_wake%": "wake.", "space_charge%": "space_charge.", "resample%": "beam.resample.",
    "bmad_com%": "bmad.", "space_charge_com%": "bmad.", "beam_init%": "beam.",
    "wake_": "wake.", "sc_": "space_charge.",
}
OLD_NAME = {
    "lat_file": "lattice", "beam_file": "beam.file", "field_file": "field.file",
    "beam_init%distribution_type(3)": "beam.z_distribution",
    "beam_init%grid(3)%x_min": "beam.z_min", "beam_init%grid(3)%x_max": "beam.z_max",
    "distribution_type(3)": "beam.z_distribution", "grid(3)%x_min": "beam.z_min", "grid(3)%x_max": "beam.z_max",
    "radiation_damping": "bmad.radiation_damping_on",
    "radiation_fluctuations": "bmad.radiation_fluctuations_on",
}
# A bare name that two tables share is routed to the first here.
BARE_TABLE = {}
for _t in ("", "run", "window", "beam", "field", "wake", "space_charge", "bmad", "beam.resample"):
    for _k in KEYS[_t]:
        BARE_TABLE.setdefault(_k, _t)
# A bare name two tables share routes to the one the old flat decks meant.
BARE_TABLE["beamlet_size"] = "beam"
BARE_TABLE["n_slice"] = "window"
BARE_TABLE["model"] = "wake"

LINE = re.compile(r"^\s*([A-Za-z_][\w%().]*)\s*=\s*(.*?)\s*(?<!\\)$")


class DeckError(ValueError):
    pass


def route(key: str) -> tuple[str, str]:
    """(table, name) for a key as a deck or a template writes it."""
    k = key.strip()
    # A list entry written with an index, field_file(2) = "...", routes by its base name
    # and keeps the index for to_toml to assemble.
    im = re.fullmatch(r"([\w%]+)\((\d+)\)", k)
    if im and im.group(1).lower() not in ("distribution_type", "grid", "beam_init%distribution_type"):
        table, name = route(im.group(1))
        return table, f"{name}({im.group(2)})"
    low = k.lower()
    if low in OLD_NAME:
        k = OLD_NAME[low]
    else:
        for pre, to in OLD_PREFIX.items():
            if low.startswith(pre):
                k = to + k[len(pre):]
                break
    if "." in k:
        table, _, name = k.rpartition(".")
        if table not in TABLES and table != "":
            raise DeckError(f"no table [{table}] for key {key}")
        if name not in KEYS.get(table, []):
            raise DeckError(f"[{table}] has no key {name} (from {key})")
        return table, name
    if k in BARE_TABLE:
        return BARE_TABLE[k], k
    raise DeckError(f"unknown key {key}")


def value_toml(raw: str) -> str:
    """A value in the old form, rendered as TOML."""
    v = raw.strip()
    if "!" in v and not v.startswith(("'", '"')):
        v = v.split("!", 1)[0].strip()
    if v.upper() in ("T", ".TRUE.", "TRUE"):
        return "true"
    if v.upper() in ("F", ".FALSE.", "FALSE"):
        return "false"
    if v.startswith(("'", '"')):
        parts = [p.strip() for p in split_commas(v)]
        strs = []
        for p in parts:
            q = p[0]
            if p[-1] != q:
                raise DeckError(f"unterminated string: {raw}")
            strs.append('"' + p[1:-1].replace("\\", "\\\\").replace('"', '\\"') + '"')
        return strs[0] if len(strs) == 1 else "[" + ", ".join(strs) + "]"
    nums = [p.strip() for p in v.split(",") if p.strip()]
    out = [number_toml(p) for p in nums]
    return out[0] if len(out) == 1 else "[" + ", ".join(out) + "]"


def number_toml(p: str) -> str:
    """One number in Fortran form as TOML. A number TOML already reads, 2048, 2e-4 or
    0.75, is kept as written, so a deck converts into the text its author typed; 1d-3,
    .5 and 5. are rewritten as 1.0e-3, 0.5 and 5.0."""
    s = p.replace("d", "e").replace("D", "e").rstrip("_rp")
    if re.fullmatch(r"[+-]?\d+", s) or re.fullmatch(r"[+-]?\d+(\.\d+)?([eE][+-]?\d+)?", s):
        return s
    m = re.fullmatch(r"([+-]?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?", s)
    if not m or (not m.group(2) and not m.group(3)):
        raise DeckError(f"not a number: {p}")
    sign, whole, frac, exp = m.group(1), m.group(2) or "0", m.group(3) or "", m.group(4)
    out = f"{sign}{whole}.{frac or '0'}"
    if exp is not None:
        out += "e" + exp
    return out


def split_commas(v: str) -> list[str]:
    """Split a comma-separated list of quoted strings, commas inside quotes kept."""
    out, cur, q = [], "", None
    for ch in v:
        if q:
            cur += ch
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            cur += ch
        elif ch == ",":
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def to_toml(flat: str) -> str:
    """Flat key = value text, in the old or the dotted form, as a TOML deck."""
    found: dict[tuple[str, str], str] = {}
    lists: dict[tuple[str, str], dict[int, str]] = {}
    for raw in flat.splitlines():
        line = raw.strip()
        if not line or line in ("/", "&end") or line.startswith(("&", "!", "#", "[")):
            continue
        m = LINE.match(line)
        if not m:
            raise DeckError(f"not a key = value line: {raw}")
        table, name = route(m.group(1))
        val = value_toml(m.group(2))
        im = re.fullmatch(r"(\w+)\((\d+)\)", name)
        if im:
            lists.setdefault((table, im.group(1)), {})[int(im.group(2))] = val
        else:
            found[(table, name)] = val
    for (table, name), entries in lists.items():
        n = max(entries)
        found[(table, name)] = "[" + ", ".join(entries.get(i, '""') for i in range(1, n + 1)) + "]"
    lines = []
    for name in KEYS[""]:
        if ("", name) in found:
            lines.append(f"{name} = {found[('', name)]}")
    for table in TABLES:
        keys = [(k, found[(table, k)]) for k in KEYS[table] if (table, k) in found]
        if not keys:
            continue
        lines.append("")
        lines.append(f"[{table}]")
        lines.extend(f"{k} = {v}" for k, v in keys)
    return "\n".join(lines) + "\n"


def write(path, flat: str) -> None:
    """Write flat text as a TOML deck at path."""
    with open(path, "w") as f:
        f.write(to_toml(flat))


def main(argv=None) -> int:
    """python -m lucifer.deck <flat file or -> <toml file or ->: convert flat text to a deck."""
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: python -m lucifer.deck <flat file | -> <toml file | ->", file=sys.stderr)
        return 2
    text = sys.stdin.read() if argv[0] == "-" else open(argv[0]).read()
    toml = to_toml(text)
    if argv[1] == "-":
        sys.stdout.write(toml)
    else:
        with open(argv[1], "w") as f:
            f.write(toml)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
