#!/usr/bin/env python3
"""Conformance checker for bmad-stats files (lucifer/doc/BMAD-STATS-SPEC.md, 1.0).

    python3 validate_bmad_stats.py run.stats.h5 [more.h5 ...]

Knows NO program and NO extension: everything here is checkable from the spec alone,
which is the property being checked (spec R42). numpy and h5py only.

Findings are reported at three levels. MUST failures mean the file does not conform.
SHOULD findings are recommendations. INFO lines report what could not be checked
generically, so silence is never mistaken for coverage: the validator can verify that
a DECLARED join key resolves (R18), not that every join key is declared; it cannot see
whether R26's strings are really the code's enums, whether R28's NaN means
not-computed, or whether R34's input trees are complete. Those belong in the writing
program's own harness.

Exit status: nonzero if any MUST failure, zero otherwise.
"""

from __future__ import annotations

import argparse
import re
import sys

import h5py
import numpy as np

# The group kinds this core document defines structure for. A kind outside this set is
# legal if declared in @kinds (extensions add their own); its structure is not checked.
CORE_KINDS = {"axis", "input", "run", "table", "provenance", "derived"}

DATASET_ATTRS = ("unit", "long_name", "description", "axes")
LIST_ATTRS = ("components", "derived_from", "extensions", "kinds")
ATTR_CAP = 60000   # bytes; the hard HDF5 limit measured at 65495 (spec R27).


def text(v):
    """An HDF5 string as str, whether bytes, str, or numpy scalar."""
    if isinstance(v, bytes):
        return v.decode("utf8", "replace").strip()
    if isinstance(v, np.ndarray) and v.ndim == 0:
        return text(v[()])
    return str(v).strip()


def as_list(v):
    """A list-valued attribute as a list of str, however it was stored."""
    a = np.atleast_1d(v)
    if a.dtype.kind in "SU":
        return [text(x) for x in a]
    return [text(x) for x in a]


def is_array(v):
    """True when the attribute was stored with a (rank >= 1) array dataspace."""
    return isinstance(v, np.ndarray) and v.ndim >= 1


class Report:
    def __init__(self, path):
        self.path = path
        self.rows = []

    def must(self, where, msg):
        self.rows.append(("MUST", where, msg))

    def should(self, where, msg):
        self.rows.append(("SHOULD", where, msg))

    def info(self, where, msg):
        self.rows.append(("INFO", where, msg))

    @property
    def n_must(self):
        return sum(1 for lv, _, _ in self.rows if lv == "MUST")

    def dump(self):
        order = {"MUST": 0, "SHOULD": 1, "INFO": 2}
        print(f"== {self.path}")
        for lv, where, msg in sorted(self.rows, key=lambda r: order[r[0]]):
            print(f"  {lv:6s} {where:44s} {msg}")
        n_s = sum(1 for lv, _, _ in self.rows if lv == "SHOULD")
        print(f"  -- {self.n_must} MUST failure(s), {n_s} SHOULD finding(s)")


# ----------------------------------------------------------------------------------
# Identity (spec R1-R6).

def check_identity(h5, rep):
    a = h5.attrs
    kinds = []

    if "file_format" not in a:
        rep.must("/", "R1: no @file_format")
    elif text(a["file_format"]) != "bmad-stats":
        rep.must("/", f"R1: @file_format is {text(a['file_format'])!r}, not 'bmad-stats'")

    if "file_format_version" not in a:
        rep.must("/", "R2: no @file_format_version")
    elif not re.fullmatch(r"\d+\.\d+", text(a["file_format_version"])):
        rep.must("/", f"R2: @file_format_version {text(a['file_format_version'])!r} "
                      "is not MAJOR.MINOR")

    if "writer" not in a:
        rep.must("/", "R3: no @writer (program name and version)")

    if "extensions" in a:                      # Absence legally means none (R4).
        if not is_array(a["extensions"]):
            rep.must("/", "R4: @extensions is not an array (R24: shape expresses arity)")
        rep.info("/", f"extensions declared: {as_list(a['extensions'])}; their structure "
                      "is checked by their own documents, not here")

    if "kinds" not in a:
        rep.must("/", "R5: no @kinds vocabulary")
    else:
        if not is_array(a["kinds"]):
            rep.must("/", "R5: @kinds is not an array (R24: shape expresses arity)")
            kinds = [k for k in text(a["kinds"]).split(",") if k]   # Lenient parse.
        else:
            kinds = as_list(a["kinds"])

    if "units_note" not in a:
        rep.should("/", "R6: no @units_note")

    return kinds


# ----------------------------------------------------------------------------------
# One data tree (spec R7-R20 and the per-dataset rules).

def check_tree(tree, prefix, kinds, rep):
    datasets, groups = {}, {}

    def note(name, obj):
        if prefix == "" and name.startswith("branches/"):
            return                          # Sub-trees are checked as their own trees.
        (datasets if isinstance(obj, h5py.Dataset) else groups)[name] = obj

    tree.visititems(note)

    def where(name):
        return f"{prefix}/{name}" if prefix else f"/{name}"

    # The axis table (R14, R15): an axis is a coords/ dataset naming itself.
    axes = {}
    if "coords" not in tree:
        rep.must(where("coords"), "R14: data tree has no coords/ group")
    else:
        for name, dset in tree["coords"].items():
            ax = text(dset.attrs.get("axes", b""))
            if ax == name:
                axes[name] = dset.shape[0]

    # Groups (R10, R5 membership).
    for name, grp in groups.items():
        kind = text(grp.attrs.get("kind", b""))
        if not kind:
            rep.must(where(name), "R10: group carries no @kind")
        elif kinds and kind not in kinds:
            rep.must(where(name), f"R5: @kind {kind!r} is not in the root @kinds")
        if not text(grp.attrs.get("description", b"")):
            rep.must(where(name), "R10: group carries no @description")

        if kind == "derived" and "derived_from" not in grp.attrs:
            rep.must(where(name), "R38: kind 'derived' requires @derived_from")
        if "derived_from" in grp.attrs:                # Any kind may declare it (R38).
            if not is_array(grp.attrs["derived_from"]):
                rep.must(where(name), "R38/R24: @derived_from is not an array")
            else:
                parent = name.rpartition("/")[0]
                sibs = set(tree[parent].keys()) if parent else set(tree.keys())
                for c in as_list(grp.attrs["derived_from"]):
                    if c not in sibs:
                        rep.must(where(name), f"R38: @derived_from names {c!r}, "
                                              "which is not a sibling")

        if "components" in grp.attrs and not is_array(grp.attrs["components"]):
            rep.must(where(name), "R39/R24: @components is not an array")

        if "struct" in grp.attrs:
            rep.info(where(name), f"@struct = {text(grp.attrs['struct'])!r} "
                                  "(verbatim component names are the writer's to check)")

    # Kinds outside the core set are extension-defined: legal when declared in @kinds,
    # their internal structure unchecked here (spec section 13).
    noncore = sorted({text(g.attrs.get("kind", b"")) for g in groups.values()}
                     - CORE_KINDS - {""})
    if noncore:
        rep.info(prefix or "/", "kinds outside the core set: " + ", ".join(noncore) +
                                ". Extension-defined structure, not checked here")

    # Datasets: the four attributes, axes resolution, shape, scalars, booleans, units.
    for name, dset in datasets.items():
        a = dset.attrs
        for req in DATASET_ATTRS:
            if req not in a or not text(a[req]):
                rep.must(where(name), f"R9: no @{req}")
        if "axes" not in a:
            continue

        ax = text(a["axes"])
        if ax == "none":
            if dset.shape != ():
                rep.must(where(name), f"R11: @axes says scalar, shape is {dset.shape}")
        else:
            names = [x for x in ax.split(",") if x]
            if len(names) != len(dset.shape):
                rep.must(where(name), f"R12: {len(dset.shape)} dimensions, "
                                      f"{len(names)} axis names")
            if len(set(names)) != len(names):
                rep.must(where(name), f"R13: an axis is named twice in {ax!r}")
            for k, axn in enumerate(names):
                if axn not in axes:
                    rep.must(where(name), f"R11: @axes names {axn!r}, which is not "
                                          "an axis of this tree")
                elif k < len(dset.shape) and dset.shape[k] != axes[axn]:
                    rep.must(where(name), f"R12: dimension {k} is {dset.shape[k]} "
                                          f"long, axis {axn!r} is {axes[axn]}")
            if dset.shape == ():
                rep.must(where(name), "R11: scalar dataset whose @axes is not 'none'")

        # Booleans (R25), both directions.
        hint = text(a.get("dtype_hint", b""))
        if dset.dtype == np.int8 and hint != "bool":
            rep.must(where(name), "R25: int8 dataset without @dtype_hint = 'bool'")
        if hint == "bool" and dset.dtype != np.int8:
            rep.must(where(name), "R25: @dtype_hint = 'bool' on a non-int8 dataset")

        # Per-entry units (R23).
        if "unit_of_axis" in a:
            fam = text(a["unit_of_axis"])
            if fam not in axes:
                rep.must(where(name), f"R23: @unit_of_axis {fam!r} is not an axis")
            elif "coords" in tree and f"{fam}_unit" not in tree["coords"]:
                rep.must(where(name), f"R23: no coords/{fam}_unit beside axis {fam!r}")
            if "unit_power" not in a:
                rep.must(where(name), "R23: @unit_of_axis without @unit_power")
            elif np.ndim(a["unit_power"]) != 0:
                rep.must(where(name), "R23/R24: @unit_power is not a true scalar")
        elif "unit_power" in a:
            rep.must(where(name), "R23: @unit_power without @unit_of_axis")

        # Declared joins and masks (R18, R19).
        if "indexes" in a:
            tgt = text(a["indexes"])
            if tgt not in axes:
                rep.must(where(name), f"R18: @indexes names {tgt!r}, not an axis")
            elif dset.dtype.kind == "i" and "coords" in tree:
                tv = tree["coords"][tgt][...]
                if tv.dtype.kind == "i" and not np.isin(dset[...], tv).all():
                    rep.must(where(name), f"R18: values fall outside axis {tgt!r}")
        if "selects" in a:
            tgt = text(a["selects"])
            if tgt not in axes:
                rep.must(where(name), f"R19: @selects names {tgt!r}, not an axis")
            elif int(np.count_nonzero(dset[...])) != axes[tgt]:
                rep.must(where(name), f"R19: selects {int(np.count_nonzero(dset[...]))} "
                                      f"entries, axis {tgt!r} has {axes[tgt]}")

        # Named list-valued attributes (R24) wherever they appear.
        for la in LIST_ATTRS:
            if la in a and not is_array(a[la]):
                rep.must(where(name), f"R24: @{la} is not an array")

        # Dataset-level derivation (R38): declared inputs must be siblings.
        if "derived_from" in a and is_array(a["derived_from"]):
            parent = name.rpartition("/")[0]
            sibs = set(tree[parent].keys()) if parent else set(tree.keys())
            for c in as_list(a["derived_from"]):
                if c not in sibs:
                    rep.must(where(name), f"R38: @derived_from names {c!r}, "
                                          "which is not a sibling")
        if "harmonic" in a and np.ndim(a["harmonic"]) != 0:
            rep.must(where(name), "R24: @harmonic is not a true scalar")

    n_joins = sum(1 for d in datasets.values() if "indexes" in d.attrs)
    n_masks = sum(1 for d in datasets.values() if "selects" in d.attrs)
    rep.info(prefix or "/", f"{len(datasets)} datasets, {len(axes)} axes, {n_joins} "
                            f"declared join key(s), {n_masks} declared mask(s). "
                            "UNdeclared joins and masks are invisible here (R42)")

    # Input tree presence (R34), advisory: a program with no inputs is conceivable.
    if not any(text(g.attrs.get("kind", b"")) == "input" for g in groups.values()):
        rep.should(prefix or "/", "R34: no kind='input' group (resolved input structs)")


# ----------------------------------------------------------------------------------
# The attribute cap (R27) and meta/ (R37), file-wide.

def check_caps_and_meta(h5, rep):
    def visit(name, obj):
        for key, val in obj.attrs.items():
            size = val.nbytes if isinstance(val, np.ndarray) else len(np.bytes_(val))
            if size > ATTR_CAP:
                rep.must(f"/{name}@{key}", f"R27: attribute is {size} bytes, near the "
                                           "64 kB cap; provenance text is a dataset")
    h5.visititems(visit)
    for key, val in h5.attrs.items():
        size = val.nbytes if isinstance(val, np.ndarray) else len(np.bytes_(val))
        if size > ATTR_CAP:
            rep.must(f"/@{key}", f"R27: attribute is {size} bytes")

    if "meta" not in h5:
        rep.must("/meta", "R37: no meta/ group")
        return
    for name, obj in h5["meta"].items():
        if not isinstance(obj, h5py.Dataset):
            rep.must(f"/meta/{name}", "R37: meta/ holds a non-dataset child")
    for leak in ("user", "cwd"):
        if leak in h5["meta"]:
            rep.should(f"/meta/{leak}", "R37: machine-local value present; legal only "
                                        "under an explicit opt-in the validator "
                                        "cannot see")


# ----------------------------------------------------------------------------------

def validate(path):
    rep = Report(path)
    with h5py.File(path, "r") as h5:
        kinds = check_identity(h5, rep)
        check_tree(h5, "", kinds, rep)
        if "branches" in h5:
            for bname, grp in h5["branches"].items():
                if isinstance(grp, h5py.Group):
                    check_tree(grp, f"/branches/{bname}", kinds, rep)
        check_caps_and_meta(h5, rep)
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+")
    ap.add_argument("--quiet", action="store_true", help="print only the summary line")
    args = ap.parse_args()

    failed = False
    for path in args.files:
        rep = validate(path)
        if args.quiet:
            print(f"{path}: {rep.n_must} MUST failure(s)")
        else:
            rep.dump()
        failed = failed or rep.n_must > 0
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
