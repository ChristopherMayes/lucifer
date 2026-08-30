"""Route flat lucifer namelist text into the three input groups.

The check scripts assemble their run configurations as flat key = value blocks
(many of them templated). lucifer's input is the three groups
&fel_params / &fel_beam_init / &fel_wavefront_init (doc/user-guide.md); this
router rewrites a flat block into them, one key at a time, with the same
old-name -> new-home mapping the program's own retired-group refusal prints.
The configurations themselves are untouched -- same keys, same values.
"""

from __future__ import annotations

# Key roots that keep their name in &fel_params.
PARAMS_LOOSE = {"lat_file", "write_wake_kernels"}
# Key roots that become global%<key> in &fel_params.
PARAMS_GLOBAL = {
    "out_root", "interlude_model", "write_diag", "write_initial",
    "load_only", "keep_escaped_field", "dump_beam_at", "dump_field_at", "ran_seed",
    "migrate", "migrate_check", "reference_run", "comb_ds_save", "track_start", "track_end",
    "source_model",
}
# Retired scalars that became bmad_com switches.
PARAMS_RENAME = {
    "radiation_damping": "bmad_com%radiation_damping_on",
    "radiation_fluctuations": "bmad_com%radiation_fluctuations_on",
}
# wake_<x> -> wake%<x>, sc_<x> -> sc%<x> handled by prefix below.

# &fel_beam_init keys, verbatim.
BEAM = {
    "beam_init", "imp", "beam_file", "dist_file", "write_dist_file", "write_opmd_file",
    "use_beam_init", "nbins", "shotnoise", "split_weights", "swap_beam_xy",
    "gen_test_weights", "imp_split_weights",
}
# Keys that become wavefront_init%<key> in &fel_wavefront_init.
WAVEFRONT = {
    "lambda0", "window_length", "window_sample", "grid_n_pts", "grid_half_width",
    "seed_power", "seed_waist_size", "seed_polarization", "harmonics",
}
WAVEFRONT_LOOSE = {"field_file"}

# Keys already written in new-style form pass through to their group.
NEWSTYLE = {"global": "params", "bmad_com": "params", "space_charge_com": "params",
            "wake": "params", "sc": "params", "wavefront_init": "wavefront"}


def to_groups(flat_text):
    """Old-style (flat) namelist text -> the three-group input text."""
    params, beam, wavefront = [], [], []
    for raw in flat_text.splitlines():
        line = raw.strip()
        if line in ("", "/") or line.startswith("&"):     # wrapper and blanks drop
            continue
        if line.startswith("!"):
            continue
        key = line.split("=", 1)[0].strip()
        root = key.split("%", 1)[0].split("(", 1)[0].strip().lower()
        rest = line[len(line.split("=", 1)[0]):]          # '= value [! comment]'
        subs = key[len(root):]                            # '(2)' etc., kept on the key
        if root in NEWSTYLE:
            {"params": params, "wavefront": wavefront}[NEWSTYLE[root]].append("  " + line)
        elif root in PARAMS_LOOSE:
            params.append("  " + line)
        elif root in PARAMS_GLOBAL:
            params.append(f"  global%{root}{subs} {rest}")
        elif root in PARAMS_RENAME:
            params.append(f"  {PARAMS_RENAME[root]} {rest}")
        elif root.startswith("wake_"):
            params.append(f"  wake%{root[5:]}{subs} {rest}")
        elif root.startswith("sc_"):
            params.append(f"  sc%{root[3:]}{subs} {rest}")
        elif root in BEAM:
            beam.append("  " + line)
        elif root in WAVEFRONT:
            wavefront.append(f"  wavefront_init%{root}{subs} {rest}")
        elif root in WAVEFRONT_LOOSE:
            wavefront.append("  " + line)
        else:
            raise ValueError(f"nml.to_groups: unmapped key: {key!r}")
    out = ["&fel_params"] + params + ["/", "&fel_beam_init"] + beam + ["/",
           "&fel_wavefront_init"] + wavefront + ["/", ""]
    return "\n".join(out)
