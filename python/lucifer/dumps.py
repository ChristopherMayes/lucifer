"""
A run's particle and field dumps, reached as files and loaded on first use.

The program names every dump it writes after the run's ``out_root``, and its footer
lists them. ``dump_files`` keys each one by what follows ``<out_root>-``, the form
LUME-Genesis keys its files by:

    final, initial          the state after the last element and before the first
    at<n>-<element>         the dump at the end of element n (dump_beam_at, dump_field_at)
    5, 6, ...               a frame, by the record it was taken at (dump_at_comb)
    pulse, escaped          the whole pulse at the exit (keep_escaped_field), and the light
                            slippage carried out of the window as it left (write_escaped),
                            fields only

A harmonic's field file adds ``-h<h>`` to the key, ``final-h3`` or ``5-h3``. Both
polarizations of a field are components of one file, so they share a key.

``Dumps`` is a mapping from those keys to the loaded objects: an openPMD-beamphysics
``ParticleGroup`` for a beam dump and a ``Wavefront`` for a field dump. A dump is read
on first access and kept, so a run's gigabytes stay on disk until asked for. Its
representation lists every key with its file size and whether it is loaded, and reads
no file. openPMD-beamphysics is imported only when a dump is loaded.

``load_particle_group(path, slice=i)`` reads one slice of a beam dump from its particle
patch alone, the i-th from 0 in the order the file stores them, the order
``lucifer.beam.read_slices`` returns.
"""
from __future__ import annotations

import pathlib
import re
from collections.abc import Mapping

__all__ = ["dump_files", "expand", "Dumps", "load_particle_group", "load_wavefront"]

BEAM, FIELD = ".beam.h5", ".wf.h5"
FRAME = re.compile(r"^(\d{6})(-h\d+)?$")
SERIES = re.compile(r"^(.*)\{(\d+)\.\.(\d+)\}(.*)$")


def expand(name):
    """The file names a footer row stands for: the row itself, or for a frame series
    written as a shell brace range, ``<root>-{000001..000012}.beam.h5``, every name in
    the range with the range's zero padding."""
    m = SERIES.match(name)
    if not m:
        return [name]
    head, first, last, tail = m.groups()
    width = len(first)
    return [f"{head}{i:0{width}d}{tail}" for i in range(int(first), int(last) + 1)]


def _key(stem):
    m = FRAME.match(stem)
    if m:
        return int(m.group(1)) if m.group(2) is None else f"{int(m.group(1))}{m.group(2)}"
    return stem


def dump_files(paths, out_root):
    """The beam and field dumps among ``paths``, as two dicts from key to path, in the
    order given. A path that is not ``<out_root>-<stem>.beam.h5`` or ``.wf.h5`` is passed
    over, the stats file and the text files among them."""
    beams, fields = {}, {}
    prefix = f"{out_root}-"
    for p in map(pathlib.Path, paths):
        name = p.name
        if not name.startswith(prefix):
            continue
        rest = name[len(prefix):]
        if rest.endswith(BEAM):
            beams[_key(rest[:-len(BEAM)])] = p
        elif rest.endswith(FIELD):
            fields[_key(rest[:-len(FIELD)])] = p
    return beams, fields


def load_particle_group(path, slice=None):
    """A beam dump as an openPMD-beamphysics ParticleGroup. With ``slice``, only that
    slice's particles, read from its particle patch alone: the records over the patch's
    span, through beamphysics.readers.particle_array, which applies the units and the
    offsets as a whole load does. So a slice is the same rows of the whole load."""
    from beamphysics import ParticleGroup
    if slice is None:
        return ParticleGroup(h5=str(path))
    import h5py
    import numpy as np
    from beamphysics.readers import particle_array, particle_paths
    with h5py.File(path, "r") as h5:
        found = particle_paths(h5)
        if len(found) != 1:
            raise ValueError(f"{pathlib.Path(path).name} holds {len(found)} particle iterations, not one")
        group = h5[found[0]]
        species = [k for k in group if isinstance(group[k], h5py.Group)]
        if len(species) != 1:
            raise ValueError(f"{pathlib.Path(path).name} holds {len(species)} species, not one")
        g = group[species[0]]
        if "particlePatches" not in g:
            raise ValueError(f"{pathlib.Path(path).name} has no particle patches to read a slice from")
        counts = g["particlePatches/numParticles"][...]
        starts = g["particlePatches/numParticlesOffset"][...]
        if not 0 <= slice < len(counts):
            raise IndexError(f"{pathlib.Path(path).name} holds {len(counts)} slices, 0 to {len(counts) - 1}, "
                             f"and slice {slice} was asked for")
        span = np.s_[int(starts[slice]):int(starts[slice]) + int(counts[slice])]
        n = int(counts[slice])
        data = {key: np.atleast_1d(particle_array(g, key, slice=span)) for key in ("x", "px", "y", "py", "z", "pz")}
        data["t"] = np.atleast_1d(particle_array(g, "t", slice=span))
        data["status"] = (np.atleast_1d(particle_array(g, "particleStatus", slice=span))
                          if "particleStatus" in g else np.ones(n, dtype=int))
        data["weight"] = np.atleast_1d(particle_array(g, "weight", slice=span))
        if "id" in g:
            data["id"] = g["id"][span]
        sp = g.attrs["speciesType"]
        data["species"] = sp.decode() if isinstance(sp, bytes) else str(sp)
    return ParticleGroup(data=data)


def load_wavefront(path):
    """A field dump as an openPMD-beamphysics Wavefront, through Wavefront.from_openpmd."""
    from beamphysics.wavefront import Wavefront
    return Wavefront.from_openpmd(str(path))


def _size(path):
    try:
        n = pathlib.Path(path).stat().st_size
    except OSError:
        return "missing"
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= scale:
            return f"{n / scale:.3g} {unit}"
    return f"{n} B"


class Dumps(Mapping):
    """A run's beam or field dumps, key to loaded object, each read on first access and
    kept. ``files`` maps every key to its path, ``loaded`` names the keys read so far, and
    ``load(key)`` reads one again. The representation lists every key with its file size
    and whether it is loaded, and reads no file."""

    def __init__(self, files, loader, kind):
        self.files = dict(files)
        self._loader = loader
        self._kind = kind
        self._loaded = {}

    def __getitem__(self, key):
        if key not in self._loaded:
            if key not in self.files:
                raise KeyError(f"no {self._kind} dump {key!r}: the run wrote "
                               f"{', '.join(map(repr, self.files)) or 'none'}")
            self._loaded[key] = self._loader(self.files[key])
        return self._loaded[key]

    def __iter__(self):
        return iter(self.files)

    def __len__(self):
        return len(self.files)

    def __contains__(self, key):
        return key in self.files

    @property
    def loaded(self):
        """The keys read so far."""
        return list(self._loaded)

    def load(self, key):
        """Read the dump of ``key`` again, keep it and return it."""
        self._loaded.pop(key, None)
        return self[key]

    def load_all(self):
        """Read every dump not yet read, frames by record after the named dumps, and
        return the keys."""
        keys = sorted(self.files, key=lambda k: (k, "") if isinstance(k, int) else (-1, str(k)))
        for k in keys:
            self[k]
        return keys

    def _rows(self):
        return [(k, self.files[k].name, _size(self.files[k]), "loaded" if k in self._loaded else "on disk")
                for k in self.files]

    def __repr__(self):
        if not self.files:
            return f"Dumps({self._kind}: none)"
        rows = self._rows()
        w = [max(len(str(r[i])) for r in rows) for i in range(4)]
        lines = [f"Dumps({self._kind}: {len(rows)} files, {len(self._loaded)} loaded)"]
        lines += [f"  {str(k):<{w[0]}}  {n:<{w[1]}}  {s:>{w[2]}}  {st}" for k, n, s, st in rows]
        return "\n".join(lines)

    def _repr_html_(self):
        from lucifer._display import table
        return table(f"{self._kind.capitalize()} dumps",
                     [(k, f"{n}, {s}, {st}") for k, n, s, st in self._rows()] or [("none", "")])
