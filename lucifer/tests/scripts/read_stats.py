#!/usr/bin/env python3
"""
Read a Lucifer statistics file. The one reader everything in the tree uses.

    from read_stats import read_stats
    st = read_stats("run.stats.h5")

    st.record                    # (nz,) the record axis, 0..nz-1
    st.s                         # (nz,) path length, a variable on the record axis
    st["beam/slice/current"]     # (nz, ns) A, as the file states it
    st.at_end                    # (nz,) bool, True where a record is an element end
    st["beam/bunch/twiss/beta"]  # (ne, 3) on the element-end axis, plane last
    st.coord("plane")            # ('x', 'y', 'z'), the projected planes
    st.dim_coords("beam/slice/sigma")   # the coordinate of each dimension, in order
    st.ele_name                  # (nz,) element name per record, gathered through ix_ele
    st.run["p0c"]                # what the run produced; st.params["global"]["ran_seed"] is what the user set

The file describes itself (manual sec:stats), so this reader hard-codes almost nothing.
Every dataset carries @unit, @long_name, @description and @axes, and every name in @axes
resolves to a coords/ dataset, so this class only exposes them: there are no unit
conversions here, no reshapes and no table of names, which is why a file that grows a
field/y group or a field/harm5 group needs no change here. Reading a value means reading
a value, and dim_coords labels a dimension without guessing from its length.

Two conveniences on top, both derived rather than stored, and both marked as such:
`ele_name` gathers lattice/name through coords/ix_ele, and `at_end` is the mask as a
bool. The file keeps one source of truth for each.

Units are documentation. The values are already SI and eV, so nothing here scales by
@unit, and neither should a caller.

The version is refused by name. This is internal development, so the format moves
without compatibility machinery: a file this reader does not know is an error rather
than a guess.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np

# The format versions this reader understands. No older ones: see the module docstring.
# 1.0 is bmad-stats with the fel extension, the planned reset from lucifer-stats 2.x.
KNOWN_VERSIONS = ("1.0",)


class Stats:
    """A Lucifer stats file, read lazily and cached."""

    def __init__(self, path):
        self.path = pathlib.Path(path)
        self._h5 = h5py.File(self.path, "r")
        self._cache: dict[str, np.ndarray] = {}

        self.file_format = _text(self._h5.attrs.get("file_format", b"?"))
        self.version = _text(self._h5.attrs.get("file_format_version", b"?"))
        if self.file_format != "bmad-stats":
            raise ValueError(f"{self.path.name} is not a bmad-stats file "
                             f"(@file_format = {self.file_format!r})")
        if self.version not in KNOWN_VERSIONS:
            raise ValueError(f"{self.path.name} is stats format {self.version}, and this "
                             f"reader knows {', '.join(KNOWN_VERSIONS)}")

        # The four documentation attributes of every dataset, and the two of every
        # group, by path. Collected in one walk: this IS the units table, and it is what
        # makes the reader generic.
        self.units: dict[str, str] = {}
        self.long_name: dict[str, str] = {}
        self.description: dict[str, str] = {}
        self.axes: dict[str, tuple[str, ...]] = {}
        self.kind: dict[str, str] = {}
        self.group_description: dict[str, str] = {}

        def note(name, obj):
            if isinstance(obj, h5py.Dataset):
                self.units[name] = _text(obj.attrs.get("unit", b""))
                self.long_name[name] = _text(obj.attrs.get("long_name", b""))
                self.description[name] = _text(obj.attrs.get("description", b""))
                ax = _text(obj.attrs.get("axes", b"none"))
                self.axes[name] = () if ax == "none" else tuple(ax.split(","))
            else:
                self.kind[name] = _text(obj.attrs.get("kind", b""))
                self.group_description[name] = _text(obj.attrs.get("description", b""))

        self._h5.visititems(note)

        # An axis is a coords/ dataset that names itself: coords/record has @axes =
        # 'record'. coords/z and coords/ix_ele name the record axis instead, which is
        # what makes them variables on it rather than axes of their own.
        self.axis_names = tuple(k for k in self._h5["coords"]
                                if self.axes.get(f"coords/{k}") == (k,))

        # run/ holds what the run produced, params/ what the user set: one subgroup
        # per honored input struct, so params is a dict of dicts.
        self.run = _scalars(self._h5["run"]) if "run" in self._h5 else {}
        self.params = {}
        if "params" in self._h5:
            for key, grp in self._h5["params"].items():
                if isinstance(grp, h5py.Group):
                    self.params[key] = _scalars(grp)

    # ------------------------------------------------------------------

    def __getitem__(self, path):
        """The dataset at `path`, as the file states it."""
        if path not in self._cache:
            self._cache[path] = np.asarray(self._h5[path][...])
        return self._cache[path]

    def __contains__(self, path):
        return path in self._h5

    def close(self):
        self._h5.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ------------------------------------------------------------------
    # The axes, reached by name rather than by shape.

    def coord(self, axis):
        """The coordinate of `axis`, as a tuple of str for a label axis and an array
        otherwise. Every name in an @axes attribute resolves here."""
        if axis not in self.axis_names:
            raise KeyError(f"{axis} is not an axis of {self.path.name}")
        val = self[f"coords/{axis}"]
        if val.dtype.kind == "S":
            return tuple(v.decode().strip() for v in val)
        return val

    def dim_coords(self, path):
        """The coordinate of each dimension of `path`, in order. This is the whole point
        of @axes: a (nz, ns, 6, 6) array says which 6 is which without a rule."""
        return tuple(self.coord(ax) for ax in self.axes[path])

    @property
    def record(self):
        """(nz,) the record axis. The axis: coords/s repeats where two records land on
        one plane, so s indexes nothing."""
        return self["coords/record"]

    @property
    def s(self):
        """(nz,) path length along the line [m], Bmad's s. A variable on the record axis.

        Named s, not z: in Bmad s is the position along the lattice and z is the
        phase-space coordinate, which is the fifth entry of coords/bmad and the third of
        coords/plane.
        """
        return self["coords/s"]

    @property
    def slice(self):
        """(ns,) the slice axis, 0..ns-1. The axis, with three positions on it."""
        return self["coords/slice"]

    @property
    def ct_slice(self):
        """(ns,) light-travel distance of each slice ahead of the reference [m].

        exact and free of beta, and the coordinate slippage counts in: one slice is
        window_sample wavelengths of it. See t_slice and z_slice.
        """
        return self["coords/ct_slice"]

    @property
    def t_slice(self):
        """(ns,) arrival time of each slice [s], -ct_slice/c exactly. More negative
        toward the head."""
        return self["coords/t_slice"]

    @property
    def z_slice(self):
        """(ns,) Bmad z of each slice at the reference beta [m].

        The one member of the three that needs a reference: a particle's own offset is
        its own beta times ct_slice, which is why the slice-to-bunch conversion stores
        every entry beta rather than one number.
        """
        return self["coords/z_slice"]

    @property
    def head_direction(self):
        """Which end of the slice index is the window head, as the file states it."""
        return _text(self._h5["coords/ct_slice"].attrs["head_direction"])

    @property
    def at_end(self):
        """(nz,) True where the record is an element end. Selects the element-end axis."""
        return self["coords/at_element_end"].astype(bool)

    @property
    def ix_ele(self):
        """(nz,) lattice element index of each record. Its values index the ele axis."""
        return self["coords/ix_ele"]

    @property
    def ele_name(self):
        """(nz,) element name per record. Derived: lattice/name gathered through ix_ele,
        which is why the file stores neither a name per record nor a second axis."""
        names = np.array([n.decode().strip() for n in self["lattice/name"]])
        return names[self.ix_ele]

    @property
    def s_end(self):
        """(ne,) the element-end positions, which the element-end arrays run over."""
        return self["coords/s_element_end"]

    # ------------------------------------------------------------------
    # What the file holds, as the file says rather than as a name suggests.

    @property
    def components(self):
        """The field polarization components present, from field/@components."""
        return _list(self._h5["field"].attrs.get("components", b""))

    @property
    def harmonics(self):
        """The harmonic numbers present beyond the fundamental, from field/@harmonics."""
        return tuple(int(h) for h in self._h5["field"].attrs.get("harmonics", []))

    def unit_of(self, path):
        """The per-entry units of `path`, or its plain @unit string.

        A dataset whose entries have different units carries @unit_of_axis and
        @unit_power instead of a parseable unit: the units live on the axis, in
        coords/<axis>_unit. Returns (units, power) in that case and (unit, 1) otherwise,
        so a caller never parses a comma list.
        """
        attrs = self._h5[path].attrs
        if "unit_of_axis" not in attrs:
            return (self.units[path],), 1
        axis = _text(attrs["unit_of_axis"])
        power = int(np.atleast_1d(attrs["unit_power"])[0])
        return self.coord(f"{axis}_unit"), power

    @property
    def meta(self):
        """The provenance group as a dict. Datasets, not attributes, so an echoed
        namelist of any size fits: see the manual on HDF5's 64 kB attribute cap. NOT a
        reproducibility record, and meta/lattice_source says so itself."""
        out = {}
        for key, dset in self._h5["meta"].items():
            val = np.atleast_1d(dset[()])[0]
            out[key] = _text(val) if isinstance(val, bytes) else val
        return out

    def derived_from(self, group):
        """What a derived group sums, from its @derived_from. A caller that adds up the
        children of field/ must skip these or it double-counts."""
        return _list(self._h5[group].attrs.get("derived_from", b""))

    def units_table(self):
        """Every dataset with its unit, label, description and axes, in file order."""
        return [(name, self.units[name], self.long_name[name], self.axes[name],
                 self.description[name]) for name in self.units]


def _text(val):
    """An HDF5 string as str, whether it arrived as bytes or as str."""
    if isinstance(val, bytes):
        return val.decode().strip()
    return str(val).strip()


def _scalars(group):
    """A group's datasets as a dict, scalars unwrapped, strings decoded."""
    out = {}
    for key, dset in group.items():
        if not isinstance(dset, h5py.Dataset):
            continue
        val = dset[()]
        if getattr(val, "shape", None) == () or np.isscalar(val):
            val = np.atleast_1d(val)[0]
            out[key] = _text(val) if isinstance(val, bytes) else val
        else:
            out[key] = np.asarray(val)
    return out


def _list(val):
    """A list-valued string attribute as a tuple of str, empty for an absent one.

    The file writes these as string arrays, length one included, so that a
    one-component file parses exactly like a two-component one: shape expresses arity.
    """
    arr = np.atleast_1d(val)
    return tuple(_text(v) for v in arr if _text(v))


def read_stats(path):
    """Open a Lucifer stats file. See the module docstring."""
    return Stats(path)


def same_data(a, b):
    """Whether two arrays hold the same data, counting NaN as equal to NaN.

    The stats file carries NaN wherever a quantity was not computed: the theta moments
    away from element ends, an empty slice's moments, the twiss of a degenerate slice.
    Two files with a NaN in the same place are identical, and plain equality would call
    every such pair different, so every identity check in the harness comes here.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if a.dtype.kind == "f":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


if __name__ == "__main__":
    import sys

    st = read_stats(sys.argv[1])
    print(f"{st.path.name}: {st.file_format} {st.version}")
    print(f"  {len(st.s)} records, {st.at_end.sum()} element ends, "
          f"{st.run['n_slice']} slices, components {st.components}, "
          f"harmonics {st.harmonics}")
    wi = st.params.get("wavefront_init", {})
    print(f"  window: lambda0 {wi.get('lambda0', float('nan')):.4e} m, sample "
          f"{wi.get('window_sample', 0)}, spacing {st.run['slice_spacing']:.4e} m, "
          f"head at {st.head_direction}")
    print(f"  axes: {', '.join(f'{a}[{len(st.coord(a))}]' for a in st.axis_names)}")
    m = st.meta
    print(f"  lattice: {m['lattice_file']}, {m['n_lattice_files']} file(s) parsed, "
          f"source {len(m['lattice_source'])} bytes")
    print(f"  inputs: {', '.join(sorted(st.params))}")
    for name, unit, label, axes, descrip in st.units_table():
        print(f"  {name:44s} [{unit:12s}] ({','.join(axes) or 'none'})  {descrip[:60]}")
