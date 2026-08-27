#!/usr/bin/env python3
"""
Read a Lucifer statistics file. The one reader everything in the tree uses.

    from read_stats import read_stats
    st = read_stats("run.stats.h5")

    st.record                    # (nz,) the record axis, 0..nz-1
    st.z                         # (nz,) path length, a VARIABLE on the record axis
    st["beam/slice/current"]     # (nz, ns) A, as the file states it
    st.at_end                    # (nz,) bool, True where a record is an element end
    st["beam/bunch/twiss/beta"]  # (ne, 6) on the element-end axis, plane last
    st.coord("plane")            # ('x', 'y', 'z', 'a', 'b', 'c')
    st.dim_coords("beam/slice/sigma")   # the coordinate of each dimension, in order
    st.ele_name                  # (nz,) element name per record, gathered through ix_ele
    st.units["beam/slice/t"]     # 's'

THE FILE DESCRIBES ITSELF (manual sec:stats), so this reader hard-codes almost nothing.
Every dataset carries @unit, @long_name, @description and @axes, and every name in @axes
resolves to a coords/ dataset, so this class only exposes them: there are no unit
conversions here, no reshapes and no table of names, which is why a file that grows a
field/y group or a field/harm5 group needs no change here. Reading a value means reading
a value, and dim_coords labels a dimension without guessing from its length.

Two conveniences on top, both derived rather than stored, and both marked as such:
`ele_name` gathers lattice/name through coords/ix_ele, and `at_end` is the mask as a
bool. The file keeps one source of truth for each.

Units are DOCUMENTATION. The values are already SI and eV, so nothing here scales by
@unit, and neither should a caller.

THE VERSION IS REFUSED BY NAME. This is internal development, so the format moves
without compatibility machinery: a file this reader does not know is an error rather
than a guess.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np

# The format versions this reader understands. No older ones: see the module docstring.
KNOWN_VERSIONS = ("2.1",)


class Stats:
    """A Lucifer stats file, read lazily and cached."""

    def __init__(self, path):
        self.path = pathlib.Path(path)
        self._h5 = h5py.File(self.path, "r")
        self._cache: dict[str, np.ndarray] = {}

        self.file_format = _text(self._h5.attrs.get("file_format", b"?"))
        self.version = _text(self._h5.attrs.get("file_format_version", b"?"))
        if self.file_format != "lucifer-stats":
            raise ValueError(f"{self.path.name} is not a Lucifer stats file "
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

        # An axis is a coords/ dataset that names ITSELF: coords/record has @axes =
        # 'record'. coords/z and coords/ix_ele name the record axis instead, which is
        # what makes them variables on it rather than axes of their own.
        self.axis_names = tuple(k for k in self._h5["coords"]
                                if self.axes.get(f"coords/{k}") == (k,))

        self.params = {}
        for key, dset in self._h5["params"].items():
            val = np.atleast_1d(dset[()])[0]
            self.params[key] = _text(val) if isinstance(val, bytes) else val

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
        """(nz,) the record axis. THE axis: coords/z repeats where two records land on
        one plane, so z indexes nothing."""
        return self["coords/record"]

    @property
    def z(self):
        """(nz,) path length along the line [m], a variable on the record axis."""
        return self["coords/z"]

    @property
    def s_slice(self):
        """(ns,) slice position in the time window [m], increasing toward the head."""
        return self["coords/s_slice"]

    @property
    def t_slice(self):
        """(ns,) arrival time of each slice [s], more negative toward the head."""
        return self["coords/t_slice"]

    @property
    def head_direction(self):
        """Which end of the slice index is the window head, as the file states it."""
        return _text(self._h5["coords/s_slice"].attrs["head_direction"])

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
        """(nz,) element name per record. DERIVED: lattice/name gathered through ix_ele,
        which is why the file stores neither a name per record nor a second axis."""
        names = np.array([n.decode().strip() for n in self["lattice/name"]])
        return names[self.ix_ele]

    @property
    def z_end(self):
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


def _list(val):
    """A comma-separated attribute as a tuple of str, empty for an absent one."""
    txt = _text(val)
    return tuple(p for p in txt.split(",") if p) if txt else ()


def read_stats(path):
    """Open a Lucifer stats file. See the module docstring."""
    return Stats(path)


def same_data(a, b):
    """Whether two arrays hold the same data, counting NaN as equal to NaN.

    The stats file carries NaN wherever a quantity was not computed: the theta moments
    away from element ends, an empty slice's moments, the twiss of a degenerate slice.
    Two files with a NaN in the same place ARE identical, and plain equality would call
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
    print(f"  {len(st.z)} records, {st.at_end.sum()} element ends, "
          f"{st.params['n_slice']} slices, components {st.components}, "
          f"harmonics {st.harmonics}")
    print(f"  window: lambda0 {st.params['lambda0']:.4e} m, sample "
          f"{st.params['window_sample']}, spacing {st.params['slice_spacing']:.4e} m, "
          f"head at {st.head_direction}")
    print(f"  axes: {', '.join(f'{a}[{len(st.coord(a))}]' for a in st.axis_names)}")
    for name, unit, label, axes, descrip in st.units_table():
        print(f"  {name:44s} [{unit:12s}] ({','.join(axes) or 'none'})  {descrip[:60]}")
