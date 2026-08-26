#!/usr/bin/env python3
"""
Read a Lucifer statistics file. The one reader everything in the tree uses.

    from read_stats import read_stats
    st = read_stats("run.stats.h5")

    st.z                      # (nz,) path length, the record axis
    st.s_slice                # (ns,) slice position in the window
    st["beam/slice/current"]  # (nz, ns) A, as the file states it
    st.at_end                 # (nz,) bool, True where a record is an element end
    st["beam/bunch/x/beta"]   # (ne,) on the element-end grid, aligned with z[at_end]
    st.ele_name               # (nz,) element name per record, gathered through ix_ele
    st.units["beam/slice/t"]  # 's'

THE FILE DESCRIBES ITSELF (manual sec:stats), so this reader hard-codes almost nothing.
Every dataset carries @unit, @description and @axes, and this class only exposes them:
there are no unit conversions here, no reshapes and no table of names, which is why a
file that grows a field/y group or a field/harm5 group needs no change here. Reading a
value means reading a value.

Two conveniences on top, both derived rather than stored, and both marked as such:
`ele_name` gathers lattice/name through coords/ix_ele, and `at_end` is the mask as a
bool. The file keeps one source of truth for each.

Units are DOCUMENTATION. The values are already SI and eV, so nothing here scales by
@unit, and neither should a caller.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np


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

        # The three documentation attributes of every dataset, by path. Collected in
        # one walk: this IS the units table, and it is what makes the reader generic.
        self.units: dict[str, str] = {}
        self.description: dict[str, str] = {}
        self.axes: dict[str, tuple[str, ...]] = {}

        def note(name, obj):
            if isinstance(obj, h5py.Dataset):
                self.units[name] = _text(obj.attrs.get("unit", b""))
                self.description[name] = _text(obj.attrs.get("description", b""))
                ax = _text(obj.attrs.get("axes", b"none"))
                self.axes[name] = () if ax == "none" else tuple(ax.split(","))

        self._h5.visititems(note)

        self.params = {}
        for key, dset in self._h5["params"].items():
            val = dset[()]
            val = np.atleast_1d(val)[0]
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
    # The axes.

    @property
    def z(self):
        """(nz,) path length along the line, the one record axis [m]."""
        return self["coords/z"]

    @property
    def s_slice(self):
        """(ns,) slice position in the time window [m], increasing toward the head."""
        return self["coords/s_slice"]

    @property
    def t_slice(self):
        """(ns,) arrival-time offset of each slice [s], negative toward the head."""
        return self["coords/t_slice"]

    @property
    def head_direction(self):
        """Which end of the slice index is the window head, as the file states it."""
        return _text(self._h5["coords/s_slice"].attrs["head_direction"])

    @property
    def at_end(self):
        """(nz,) True where the record is an element end. Selects the element-end grid."""
        return self["coords/at_element_end"].astype(bool)

    @property
    def ix_ele(self):
        """(nz,) lattice element index of each record. Indexes the lattice/ table."""
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
        return self.z[self.at_end]

    # ------------------------------------------------------------------
    # What the file holds, for a caller that wants to branch on it rather than guess.

    @property
    def components(self):
        """The field polarization components present, ('x',) or ('x', 'y')."""
        return tuple(c for c in ("x", "y") if f"field/{c}" in self)

    @property
    def harmonics(self):
        """The harmonic numbers present beyond the fundamental, in order."""
        out = [int(k[4:]) for k in self._h5["field"] if k.startswith("harm")]
        return tuple(sorted(out))

    def units_table(self):
        """Every dataset with its unit, description and axes, in file order."""
        return [(name, self.units[name], self.axes[name], self.description[name])
                for name in self.units]


def _text(val):
    """An HDF5 string as str, whether it arrived as bytes or as str."""
    if isinstance(val, bytes):
        return val.decode().strip()
    return str(val).strip()


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
    for name, unit, axes, descrip in st.units_table():
        print(f"  {name:44s} [{unit:12s}] ({','.join(axes) or 'none'})  {descrip[:60]}")
