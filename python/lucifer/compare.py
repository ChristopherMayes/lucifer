"""Compare a Lucifer run with a reference run, Genesis 1.3 Version 4's or another Lucifer run.

    from lucifer.compare import read_curves, aligned, differences
    run, ref = read_curves("run.stats.h5"), read_curves("Aramis.out.h5")
    reason = aligned(run, ref)          # None, or why the two cannot be differenced
    d = differences(run, ref)           # each curve's difference, divided as its label says

The curves both codes record along the line are read into one form, ``Curves``: the
record positions s, the power, the bunching, the mean energy and the rms sizes, each per
record and slice, and the slice spacing. A Lucifer run gives them from its stats file and
Genesis4 from its ``.out.h5``, the datasets Lattice/zplot, Field/power, Beam/bunching,
Beam/energy, Beam/xsize, Beam/ysize and Global/s. Genesis4's energy is gamma and is
converted to eV with Bmad's electron mass, so the two codes' energies are in one unit.

A difference is Lucifer's value less the reference's, divided as ``DIVISOR`` says for
each curve: the power elementwise by the reference's power at each record and slice,
floored at 1e-9 of its peak so a slice that holds no light yet cannot divide by zero, the
bunching by the reference's peak bunching, the mean energy change by the reference's
largest change, the exit power across the window by its peak, and a field amplitude by
the reference's peak amplitude. The power and bunching rules are the ones
tests/scripts/compare_fel.py checks the tiers by, so the largest difference here is the
tier's number for that curve. A time-dependent run's difference at a record is the one
over its slices of largest magnitude, with its sign, so the largest over the records is
the largest over everything.

A difference is taken only where the two runs share their record positions to
``S_TOLERANCE`` and their slice count and spacing, which ``aligned`` tests. Nothing here
interpolates: a run recorded at other positions is overlaid and not differenced, since
an interpolated difference would hide the mismatch inside a number that looks small.

numpy and h5py alone. ``lucifer.plot.compare`` draws the figure.
"""
from __future__ import annotations

import dataclasses
import pathlib

import h5py
import numpy as np

__all__ = ["Curves", "Difference", "read_genesis_out", "from_stats", "read_curves",
           "aligned", "differences", "field_level"]

C_LIGHT = 2.99792458e8
M_ELECTRON = 0.51099895069e6   # Bmad's m_electron [eV]
# compare_fel.py's bound on the two codes' record positions [m], and its floor on the
# power denominator as a fraction of the reference's peak.
S_TOLERANCE = 1e-9
POWER_FLOOR = 1e-9

# What each difference is divided by, for a caption.
DIVISOR = {
    "power": "the reference's power at each record and slice",
    "bunching": "the reference's peak bunching",
    "energy": "the reference's largest mean energy change",
    "exit_power": "the reference's peak exit power",
    "size_x": "the reference's rms size at each record",
    "size_y": "the reference's rms size at each record",
    "field": "the reference's peak field amplitude",
}


@dataclasses.dataclass
class Curves:
    """What one run recorded along the line, per record and slice, in one form for both
    codes. The arrays are (nz, ns): nz records at positions ``s`` and ns slices, one in
    steady state."""

    name: str
    """The run's name in a legend: Lucifer, Genesis4, or a stats file's stem."""
    kind: str
    """``lucifer`` or ``genesis4``, the code that wrote the file."""
    s: np.ndarray
    """(nz,) the record positions along the line [m]."""
    power: np.ndarray
    """(nz, ns) the radiation power [W]."""
    bunching: np.ndarray
    """(nz, ns) the bunching factor's magnitude."""
    energy: np.ndarray
    """(nz, ns) the mean electron energy of each slice [eV]."""
    size_x: np.ndarray
    """(nz, ns) the rms beam size in x [m]."""
    size_y: np.ndarray
    """(nz, ns) the rms beam size in y [m]."""
    ds: float | None
    """The slice spacing [m], and None in steady state."""
    source: pathlib.Path | None = None
    """The file the curves were read from."""

    @property
    def n_slice(self):
        """The slice count."""
        return self.power.shape[1]

    @property
    def pulse_energy(self):
        """(nz,) the energy in the window [J], the power summed over the slices times the
        slice spacing over c. None in steady state, where there is no window."""
        if self.ds is None:
            return None
        return self.power.sum(axis=1) * self.ds / C_LIGHT

    @property
    def energy_change(self):
        """(nz, ns) the mean energy less its value at the first record [eV]."""
        return self.energy - self.energy[0]

    def __repr__(self):
        nz, ns = self.power.shape
        return f"Curves({self.name}, {nz} records to s = {self.s[-1]:.4g} m, {ns} slice{'s' if ns > 1 else ''})"


@dataclasses.dataclass
class Difference:
    """One curve's difference between a run and its reference, divided as
    ``DIVISOR[key]`` says, along ``x``."""

    key: str
    """power, bunching, energy, exit_power, size_x, size_y or field."""
    x: np.ndarray
    """(n,) the abscissa: s along the line [m], the position along the window [m], or x
    across the grid [m]."""
    values: np.ndarray
    """(n,) the difference, dimensionless."""

    @property
    def largest(self):
        """The largest magnitude of the difference."""
        return float(np.nanmax(np.abs(self.values))) if self.values.size else 0.0

    def __repr__(self):
        return f"Difference({self.key}, largest {self.largest:.3e})"


def read_genesis_out(path, name="Genesis4"):
    """The curves of a Genesis4 ``.out.h5``. Beam/energy is gamma there and comes back in
    eV. The slice spacing is the step of Global/s, and None for one slice."""
    path = pathlib.Path(path)
    with h5py.File(path) as h5:
        s = np.asarray(h5["Lattice/zplot"][...], dtype=float)
        power = np.asarray(h5["Field/power"][...], dtype=float)
        bunching = np.asarray(h5["Beam/bunching"][...], dtype=float)
        energy = np.asarray(h5["Beam/energy"][...], dtype=float) * M_ELECTRON
        size_x = np.asarray(h5["Beam/xsize"][...], dtype=float)
        size_y = np.asarray(h5["Beam/ysize"][...], dtype=float)
        window = np.asarray(h5["Global/s"][...], dtype=float).ravel()
    ds = float(window[1] - window[0]) if window.size > 1 else None
    return Curves(name, "genesis4", s, power, bunching, energy, size_x, size_y, ds, path)


def from_stats(st, name=None):
    """The curves of a Lucifer run from its ``Stats``: coords/s, field/total/power,
    beam/slice/bunching, beam/slice/energy and the rms sizes from beam/slice/sigma. The
    slice spacing is the step of coords/ct_slice. ``name`` defaults to the file's stem."""
    sigma = np.asarray(st["beam/slice/sigma"], dtype=float)
    ct = np.asarray(st.ct_slice, dtype=float)
    ds = float(abs(ct[1] - ct[0])) if ct.size > 1 else None
    if name is None:
        name = st.path.name.removesuffix(".stats.h5")
    return Curves(name, "lucifer", np.asarray(st.s, dtype=float),
                  np.asarray(st["field/total/power"], dtype=float),
                  np.asarray(st["beam/slice/bunching"], dtype=float),
                  np.asarray(st["beam/slice/energy"], dtype=float),
                  np.sqrt(np.maximum(0.0, sigma[:, :, 0, 0])),
                  np.sqrt(np.maximum(0.0, sigma[:, :, 2, 2])), ds, st.path)


def read_curves(source):
    """The curves of ``source``: a ``Curves`` as it is, a ``Stats`` through
    ``from_stats``, or a path, read by what the file is: a bmad-stats file by its root
    attribute, a Genesis4 output by its Lattice/zplot and Field/power. Any other file is
    refused with what it does hold."""
    from lucifer.stats import Stats, read_stats

    if isinstance(source, Curves):
        return source
    if isinstance(source, Stats):
        return from_stats(source)
    path = pathlib.Path(source)
    with h5py.File(path) as h5:
        fmt = h5.attrs.get("file_format", b"")
        fmt = fmt.decode() if isinstance(fmt, bytes) else str(fmt)
        is_genesis = "Lattice/zplot" in h5 and "Field/power" in h5
        roots = list(h5.keys())
    if fmt == "bmad-stats":
        with read_stats(path) as st:
            return from_stats(st)
    if is_genesis:
        return read_genesis_out(path)
    raise ValueError(f"{path.name} is neither a bmad-stats file nor a Genesis4 output: "
                     f"its root holds {', '.join(roots)}")


def aligned(run, ref):
    """None when the two runs can be differenced record by record and slice by slice,
    else the reason they cannot: a different record count, a record position differing by
    more than ``S_TOLERANCE``, a different slice count, or a different slice spacing."""
    if len(run.s) != len(ref.s):
        return f"{run.name} has {len(run.s)} records and {ref.name} {len(ref.s)}"
    gap = float(np.abs(run.s - ref.s).max()) if len(run.s) else 0.0
    if gap > S_TOLERANCE:
        return f"the record positions of {run.name} and {ref.name} differ by up to {gap:.3g} m"
    if run.n_slice != ref.n_slice:
        return f"{run.name} has {run.n_slice} slices and {ref.name} {ref.n_slice}"
    if run.ds is not None and abs(run.ds - ref.ds) > 1e-9 * abs(ref.ds):
        return f"the slice spacings of {run.name} and {ref.name} are {run.ds:.6g} and {ref.ds:.6g} m"
    return None


def _largest_signed(d):
    """(nz,) the entry of each row of largest magnitude, with its sign."""
    if d.shape[1] == 1:
        return d[:, 0]
    i = np.argmax(np.abs(d), axis=1)
    return d[np.arange(d.shape[0]), i]


def _over(num, den):
    """num / den, and zero where den is zero."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den != 0, num / np.where(den != 0, den, 1.0), 0.0)


def differences(run, ref, fields=None):
    """The differences of ``run`` from ``ref``, by key: power, bunching and energy along
    s; exit_power across the window for a time-dependent run; size_x and size_y along s
    for a steady-state one; and field, the amplitude along the row through the grid's
    center at y = 0, when ``fields`` names the two runs' field dumps at the exit, read
    with ``lucifer.field.read_field``. The runs must be aligned, as ``aligned`` tests,
    and the two fields on one grid, else a ValueError names the reason."""
    reason = aligned(run, ref)
    if reason:
        raise ValueError(reason)
    out = {}
    den = np.maximum(np.abs(ref.power), POWER_FLOOR * np.abs(ref.power).max())
    out["power"] = Difference("power", run.s, _largest_signed(_over(run.power - ref.power, den)))
    peak = np.abs(ref.bunching).max()
    out["bunching"] = Difference("bunching", run.s,
                                 _largest_signed(_over(run.bunching - ref.bunching, np.full_like(ref.bunching, peak))))
    change = ref.energy_change
    scale = np.abs(change).max()
    out["energy"] = Difference("energy", run.s,
                               _largest_signed(_over(run.energy_change - change, np.full_like(change, scale))))
    if run.n_slice > 1:
        x = np.arange(run.n_slice) * run.ds
        exit_peak = np.abs(ref.power[-1]).max()
        out["exit_power"] = Difference("exit_power", x,
                                       _over(run.power[-1] - ref.power[-1], np.full(run.n_slice, exit_peak)))
    else:
        out["size_x"] = Difference("size_x", run.s, _over(run.size_x[:, 0] - ref.size_x[:, 0], ref.size_x[:, 0]))
        out["size_y"] = Difference("size_y", run.s, _over(run.size_y[:, 0] - ref.size_y[:, 0], ref.size_y[:, 0]))
    if fields is not None:
        fa, fb = _fields(fields, run.n_slice)
        ua, ub = fa["u"], fb["u"]
        if ua.shape != ub.shape or abs(fa["dx"] - fb["dx"]) > 1e-9 * fb["dx"]:
            raise ValueError(f"the field grids differ: {ua.shape} at {fa['dx']:.4g} m and "
                             f"{ub.shape} at {fb['dx']:.4g} m")
        i = int(np.argmax([(np.abs(u) ** 2).sum() for u in ub])) if ub.shape[0] > 1 else 0
        c = ub.shape[1] // 2
        x = (np.arange(ub.shape[2]) - ub.shape[2] // 2) * fb["dx"]
        out["field"] = Difference("field", x,
                                  _over(np.abs(ua[i, c]) - np.abs(ub[i, c]), np.full(ub.shape[2], np.abs(ub).max())))
    return out


def _fields(fields, n_slice):
    """The two field dumps ``fields`` names, read, or passed through when already read."""
    from lucifer.field import read_field

    a, b = fields
    fa = a if isinstance(a, dict) else read_field(a, nslice=n_slice)
    fb = b if isinstance(b, dict) else read_field(b, nslice=n_slice)
    return fa, fb


def field_level(fields, n_slice=None):
    """The largest difference of the two complex fields over every slice and grid point,
    over the reference's peak amplitude: compare_fel.py's final-field number. ``fields``
    names the run's and the reference's dumps, or holds them read."""
    fa, fb = _fields(fields, n_slice)
    return float(np.abs(fa["u"] - fb["u"]).max() / np.abs(fb["u"]).max())
