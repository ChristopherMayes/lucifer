#!/usr/bin/env python3
"""
The spontaneous radiation of one undulator against the closed forms of undulator
radiation (docs/undulator-spectrum.md, tests/scripts/undulator_radiation.py). Lucifer
alone: the references are formulas.

Two undulators are run. tests/bmad/und_flash1.bmad is the FLASH1 undulator, planar, 165
periods of 27.3 mm at 668 MeV, radiating at 13.7 nm, and a third of it, 55 periods.
tests/bmad/und_delta.bmad is the LCLS Delta afterburner in its circular mode, helical,
100 periods of 32 mm at 4 GeV, radiating at 707 eV.

The emitter is one loaded slice of one wavelength in a longer window, 62 electrons with
shot noise, far below any gain. The field slips ahead of the slice by one wavelength a
period, so each plane of the pulse file is what the slice radiated in one period, drifted
to the exit. The sum of a plane over x and y is its component at zero angle, which the
drift leaves alone. With S the slice's sum of w exp(-i theta) over its particles, w a
particle's electrons, the script holds:

  the emitter   each plane's zero-angle component is one constant times conj(S) at the
                period that radiated it, so the spectrum on the axis is that of S along
                the undulator and one run gives it with no average over seeds.
  on axis       the line against sinc^2(pi N nu), nu = omega/omega_1 - 1: its width over
                0.8859/N, its level at the zeros nu = m/N, and its largest difference
                from sinc^2 over four lobes. At 165 periods, at 55, and in the helical
                undulator.
  the level     the energy per unit angular frequency and solid angle at the peak, over
                |sum of S|^2, against W1(0) of undulator_radiation.py.
  with angle    a row of the grid's transform is one transverse wavenumber k_perp, and
                x = gamma_z^2 (k_perp/k_0)^2. The ridge of the transform in x, and in y,
                and along the pulse, against nu = -x, which is what the envelope equation
                gives, and its distance at the widest wavenumber from the resonance of
                undulator radiation, nu (1 + nu) = -x. The level along the ridge over the
                spectrum of what the slice deposits in one step, which is 1 where the
                model applies no angular factor, and the formula's W1 at that wavenumber
                beside it. With run.source_filter off and on.
  spectral      the 55-period undulator with field.propagator = "spectral", where each
                frequency diffracts at its own wavenumber: the ridge against the
                resonance.
  the beam      the same lines with FLASH1's emittance, with its energy spread, in a
                window shorter than the slippage and in slices of three wavelengths.
  unaveraged    the 55-period undulator with fel_method = unaveraged.
  a bunch       1024 loaded slices with shot noise at SEEDS seeds: the spectrum on the
                axis of each run against that of its own S, and the energy on the axis
                over what one electron radiates times the electrons.

The curves the page draws are written to undulator_spectrum.json in the work directory,
from which report_undulator_spectrum.py draws them. The pulse files are removed.

Run by the validation suite. Exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from checkrec import record
from lucifer.field import read_pulse
from lucifer.flat import to_toml

import pool
import undulator_radiation as ur

FAILED = False
SEEDS = 8
N_BUNCH = 1024              # Loaded slices of the bunch.

FLASH = dict(lattice="und_flash1.bmad", e_tot=6.684941e8, b_max=0.47, l_period=0.0273, helical=False,
             lambda0=1.37e-8, half_width=1.0e-3)
DELTA = dict(lattice="und_delta.bmad", e_tot=4e9, b_max=0.80015309, l_period=0.032, helical=True,
             lambda0=1.753666e-9, half_width=2.5e-4)

FLAT = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lattice}"
  out_root = "{root}"
  source_filter = {filter}
  lambda0 = {lambda0}
  beam.n_particle = {n_particle}
  beam.bunch_charge = {charge}
  beam.z_distribution = "GRID"
  beam.z_min = {z_min}
  beam.z_max = {z_max}
  beam.sig_pz = {sig_pz}
  beam.a_norm_emit = {emit}
  beam.b_norm_emit = {emit}
  beamlet_size = 8
  grid_n_pts = {grid}
  grid_half_width = {half_width}
  window_length = {window}
  n_wavelength = {n_wavelength}
  shot_noise = T
  ran_seed = {seed}
{extra}&end
"""

WRAP = "call, file = {lattice}\n{lines}\n"

# The bounds. Each is set from what the debug and the production build measured, which
# agree to every printed digit (docs/generated/validation-measured.md). TIGHT holds the
# runs of a beam with no size and no energy spread, where the slice's S does not change
# along the undulator. BEAM holds the runs of FLASH1's emittance or energy spread, where
# each beamlet radiates a line displaced by its own detuning and the sum, with the shot
# noise's phases, is one draw of a speckle: the peak moves by a fraction of a width and
# the zeros fill in, as the seed decides.
TIGHT = dict(emitter=1e-7, peak=1e-3, width=1e-3, zeros=1e-6, sinc=5e-4, level=1e-3)
BEAM = dict(emitter=3e-3, peak=1.0, width=0.2, zeros=0.3, sinc=0.4, level=1e-2)
BOUNDS = dict(
    flash=dict(ridge=2e-3, exact=[4.9, 5.5], deposit=1e-5, factor_x=[0.335, 0.351], factor_y=[1.03, 1.065]),
    delta=dict(ridge=2e-3, exact=[7.6, 8.3], deposit=1e-5, factor_x=[0.375, 0.393], factor_y=[0.375, 0.393]),
    size=dict(BEAM, ridge=4.0, exact=[3.0, 8.0]),
    spread=BEAM,
    window=TIGHT,
    filter=dict(BEAM, ridge=1.5, kept=0.1, cut=1e-3),
    spectral=dict(peak=1e-3, width=4e-2, zeros=2e-3, sinc=4e-2, level=6e-2, ridge=0.1, linear=[-0.55, -0.38]),
    slice3=dict(width=1e-3, level=1e-3),
    unavg=dict(width=[1.03, 1.06], sinc=0.08, level=[0.92, 0.96], ridge=0.06, exact=[1.6, 1.9], deposit=0.03,
               factor_x=[0.335, 0.360], factor_y=[1.03, 1.07]),
    bunch=dict(emitter=1e-4, energy=[0.85, 1.15], scatter=0.3),
)


def check(name, value, bound, kind="max"):
    """Print and record one check of a number against its bound."""
    global FAILED
    if kind == "range":
        ok = bound[0] <= value <= bound[1]
        text = f"range {bound[0]:.3e} .. {bound[1]:.3e}"
    else:
        ok = {"max": value <= bound, "min": value >= bound, "equal": value == bound}[kind]
        text = f"{kind} {bound:.1e}"
    print(f"--- {name}: {value:.4e} ({text}) {'ok' if ok else '** FAIL **'}")
    record(name, value, bound, kind, ok)
    FAILED = FAILED or not ok


class Run:
    """One run of the emitter: what it was started with, and what it radiated."""

    def __init__(self, name, und, n_period, window, grid=128, emit=1e-14, sig_pz=1e-6, filter="F",
                 n_wavelength=1, n_loaded=1, n_particle=512, seed=777, lines="", extra=""):
        self.name, self.und, self.n_period, self.window = name, und, n_period, window
        self.grid, self.emit, self.sig_pz, self.filter = grid, emit, sig_pz, filter
        self.n_wavelength, self.n_loaded, self.n_particle = n_wavelength, n_loaded, n_particle
        self.seed, self.lines, self.extra = seed, lines, extra
        self.u = ur.Undulator(und["e_tot"], und["b_max"], und["l_period"], n_period, und["helical"])

    def start(self, exe, wd, threads="4"):
        """Write the input file and run it."""
        lam = self.und["lambda0"]
        lattice = self.und["lattice"]
        if self.lines:
            lattice = f"{self.name}.bmad"
            (wd / lattice).write_text(WRAP.format(lattice=self.und["lattice"], lines=self.lines))
        # The window is centered on the bunch, and the current of a GRID bunch is
        # evaluated at the slice centers. With as many slices in the window as are loaded,
        # odd or even alike, the loaded slices' centers are half a slice inside the bunch
        # and no center is at an edge of it, where rounding would decide.
        half = 0.5 * self.n_loaded * self.n_wavelength * lam
        if (round(self.window / self.n_wavelength) - self.n_loaded) % 2:
            sys.exit(f"{self.name}: the window's slice count and the loaded slices differ in parity")
        text = FLAT.format(lattice=lattice, root=self.name, filter=self.filter, lambda0=lam,
                           n_particle=self.n_particle, charge=1e-17 * self.n_loaded, z_min=-half, z_max=half,
                           sig_pz=self.sig_pz, emit=self.emit, grid=self.grid,
                           half_width=self.und["half_width"], window=self.window * lam,
                           n_wavelength=self.n_wavelength, seed=self.seed, extra=self.extra)
        (wd / f"{self.name}.toml").write_text(to_toml(text))
        r = subprocess.run([str(exe), f"{self.name}.toml"], cwd=wd, capture_output=True, text=True,
                           env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
        if r.returncode != 0:
            print(f"FAIL: {self.name} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
            sys.exit(1)

    def read(self, wd, keep_planes=True):
        """Read the pulse and the slices' S along the undulator, and remove the pulse."""
        path = wd / f"{self.name}-pulse.wf.h5"
        w = read_pulse(path).wavefront
        field = np.asarray(w.Ex)
        self.dx, self.dz = float(w.dx), float(w.dz)
        self.cell = float(w.dx) * float(w.dy)
        self.a0 = field.sum(axis=(0, 1))                       # The zero-angle component of each plane.
        if keep_planes:
            f2 = np.fft.fft2(field, axes=(0, 1))
            self.along = {"x": f2[:, 0, :], "y": f2[0, :, :]}  # (angle, plane), the other angle zero.
            self.theta = np.fft.fftfreq(field.shape[0], self.dx) * self.und["lambda0"]
        with h5py.File(wd / f"{self.name}.stats.h5") as h:
            b = np.asarray(h["beam/slice/bunching"])
            phase = np.asarray(h["beam/slice/bunching_phase"])
            electrons = np.asarray(h["beam/slice/charge_live"]) / ur.E_CHARGE
        self.s = electrons * b * np.exp(1j * phase)             # (record, slice)
        self.electrons = float(electrons[0].sum())
        self.loaded = np.nonzero(electrons[0])[0]
        path.unlink()
        (wd / f"{self.name}-final.wf.h5").unlink(missing_ok=True)

    # The pulse as the slices' S predicts it. A slice radiates one plane a period, which
    # is then one plane further ahead of the slice at every later period. Record j of the
    # stats file is the beam after period j, and the plane of period j ends n_period - j
    # planes ahead of the plane the last period radiated.

    def predicted(self):
        """sum of conj(S) over the periods and slices that radiated into each plane."""
        n = self.n_period
        out = np.zeros(len(self.a0), dtype=complex)
        first = self.first_plane - self.loaded[0]
        for i in self.loaded:
            for j in range(1, n + 1):
                out[i + first + n - j] += np.conj(self.s[j, i])
        return out

    @property
    def first_plane(self):
        """The plane the first loaded slice's last period radiated."""
        lit = np.nonzero(np.abs(self.a0) > 1e-9 * np.abs(self.a0).max())[0]
        return int(lit.min())

    def spectrum(self, planes, nu):
        """|sum of planes exp(-2 pi i nu p)|^2 at the relative frequencies nu, p in
        wavelengths along the pulse."""
        p = np.arange(len(planes)) * self.dz / self.und["lambda0"]
        return np.abs(np.exp(-2j * np.pi * np.outer(np.atleast_1d(nu), p)) @ planes) ** 2

    def level(self, nu):
        """The energy per unit angular frequency and solid angle on the axis at nu
        [J s / sr], from the field in V/m: (eps0 c / 2) k^2 |E(k = 0, omega)|^2 / (2 pi)^3."""
        k = 2 * np.pi / self.und["lambda0"]
        e = self.spectrum(self.a0, nu) * (self.cell * self.dz / ur.C_LIGHT) ** 2
        return ur.EPS0 * ur.C_LIGHT / 2 * k ** 2 * e / (2 * np.pi) ** 3


def peak_of(f, lo, hi, n=4001):
    """The place and height of the largest value of f on [lo, hi], refined by a parabola."""
    x = np.linspace(lo, hi, n)
    y = f(x)
    k = int(np.clip(np.argmax(y), 1, n - 2))
    a, b, c = y[k - 1], y[k], y[k + 1]
    d = 0.5 * (a - c) / (a - 2 * b + c)
    return x[k] + d * (x[1] - x[0]), b - 0.25 * (a - c) * d


def width_of(f, center, height, span):
    """The full width at half of height of f about center, by bisection."""
    def edge(sign):
        lo, hi = 0.0, span
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if f(np.array([center + sign * mid]))[0] > 0.5 * height:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    return edge(1) + edge(-1)


def on_axis(run, label, bounds, rec):
    """The line on the axis: the emitter, the width, the zeros, sinc^2 and the level."""
    n, u = run.n_period, run.u
    pred = run.predicted()
    if "emitter" in bounds:
        scale = np.vdot(pred, run.a0) / np.vdot(pred, pred)     # The one constant, a least-squares fit.
        check(f"{label}: planes that differ from the constant times conj(S), over the largest",
              float(np.abs(run.a0 - scale * pred).max() / np.abs(run.a0).max()), bounds["emitter"])
    f = lambda nu: run.spectrum(run.a0, nu)
    center, height = peak_of(f, -0.5 / n, 0.5 / n)
    # The resonance of the lattice's undulator, against the carrier the run was given.
    detuning = run.und["lambda0"] / u.wavelength - 1
    check(f"{label}: peak of the line less the undulator's resonance, in line widths",
          abs(center - detuning) * n / ur.LINE_WIDTH, bounds["peak"])
    width = width_of(f, center, height, 1.0 / n)
    check(f"{label}: width of the line over 0.8859/N", width * n / ur.LINE_WIDTH,
          [1 - bounds["width"], 1 + bounds["width"]], "range")
    zeros = f(center + np.arange(1, 4) / n) / height
    check(f"{label}: line at the zeros nu = 1/N, 2/N and 3/N, over the peak", float(zeros.max()), bounds["zeros"])
    nu = center + np.linspace(-4.0 / n, 4.0 / n, 801)
    measured = f(nu) / height
    check(f"{label}: largest difference from sinc^2 over four lobes", float(np.abs(measured - ur.line(n, nu - center)).max()),
          bounds["sinc"])
    ratio = float(run.level(center)[0] / (u.on_axis_level() * run.spectrum(pred, center)[0]))
    check(f"{label}: level on the axis over the closed form", ratio, [1 - bounds["level"], 1 + bounds["level"]], "range")
    rec[run.name] = dict(n_period=n, helical=u.helical, nu=(nu - center).tolist(), line=measured.tolist(),
                         center=center, detuning=detuning, width=width * n / ur.LINE_WIDTH, level=ratio,
                         planes=(np.abs(run.a0) / np.abs(run.a0).max()).tolist(), first_plane=run.first_plane)
    return center


def resonance(x):
    """The relative frequency of the fundamental at the transverse wavenumber k_perp,
    with x = gamma_z^2 (k_perp / k_0)^2. The fundamental at the angle theta has the
    wavenumber k = k_0 / (1 + gamma_z^2 theta^2), and k_perp = k theta, so that
    nu (1 + nu) = -x. It is -x/(1 + x') at x' = gamma_z^2 theta^2, and -x to first order."""
    return 0.5 * (np.sqrt(1 - 4 * x) - 1)


def with_angle(run, label, x_max, bounds, rec, deposit=True, exact=True, spectral=False):
    """The ridge and the level along it, in x and in y. A row of the grid's transform is
    one transverse wavenumber k_perp, and x is gamma_z^2 (k_perp / k_0)^2."""
    n, u = run.n_period, run.u
    lit = slice(run.first_plane, run.first_plane + n)
    out = {}
    for axis in ("x", "y"):
        planes = run.along[axis]
        own = (np.abs(planes[:, lit]) ** 2).mean(axis=1)        # What one step deposits, at each wavenumber.
        rows = []
        on = None
        for i in range(0, len(run.theta) // 2):
            x = u.gamma_z2 * run.theta[i] ** 2
            if x > x_max:
                break
            f = lambda nu, i=i: run.spectrum(planes[i], nu)
            center, height = peak_of(f, resonance(x) - 1.5 / n, -x + 1.5 / n, 2401)
            on = height if i == 0 else on
            # The formula at this wavenumber: at the angle k_perp / k of its resonance,
            # and per unit area of k_perp, which is per solid angle over k^2.
            theta = run.theta[i] / (1 + resonance(x))
            factor = u.angular_factor(theta, 0.0) if axis == "x" else u.angular_factor(0.0, theta)
            rows.append((x, center, height / on, own[i] / own[0], factor / (1 + resonance(x)) ** 2))
        out[axis] = np.array(rows)
    both = np.vstack([out["x"], out["y"]])
    center0 = out["x"][0, 1]
    last = out["x"][-1]
    widths = n / ur.LINE_WIDTH
    if spectral:
        check(f"{label}: ridge less the resonance, in line widths, at its largest over the wavenumbers",
              float(np.abs(both[:, 1] - center0 - resonance(both[:, 0])).max() * widths), bounds["ridge"])
        check(f"{label}: ridge at the widest wavenumber within x = {x_max:g}, less -x, in line widths",
              float((last[1] - center0 + last[0]) * widths), bounds["linear"], "range")
    else:
        check(f"{label}: ridge less -x, in line widths, at its largest over the wavenumbers",
              float(np.abs(both[:, 1] - center0 + both[:, 0]).max() * widths), bounds["ridge"])
        if exact:
            check(f"{label}: ridge at the widest wavenumber within x = {x_max:g}, less the resonance, in line widths",
                  float((last[1] - center0 - resonance(last[0])) * widths), bounds["exact"], "range")
    if deposit:
        check(f"{label}: level along the ridge over the deposit's own spectrum, less 1",
              float(np.abs(both[:, 2] / both[:, 3] - 1).max()), bounds["deposit"])
        for axis in ("x", "y"):
            last = out[axis][-1]
            check(f"{label}: the formula's factor at the widest wavenumber along {axis}, over the model's",
                  float(last[4] / (last[2] / last[3])), bounds[f"factor_{axis}"], "range")
    rec[run.name].update({f"ridge_{a}": out[a].tolist() for a in out})
    return out


def source_filter(on, off, wd, rec):
    """The filter on against the filter off, the same beam: what is kept and what is cut."""
    label = "FLASH1, source filter on"
    with h5py.File(wd / f"{on.name}.stats.h5") as h:
        edge = float(np.asarray(h["setup/filter_edge"]).ravel()[0])
    x_edge = on.u.gamma_z2 * edge ** 2
    on_axis(on, label, BOUNDS["filter"], rec)
    with_angle(on, label, 0.5 * x_edge, BOUNDS["filter"], rec, deposit=False, exact=False)
    lit = slice(on.first_plane, on.first_plane + on.n_period)
    kept = {}
    for axis in ("x", "y"):
        ratio = (np.abs(on.along[axis][:, lit]) ** 2).sum(axis=1) / (np.abs(off.along[axis][:, lit]) ** 2).sum(axis=1)
        half = len(on.theta) // 2
        kept[axis] = np.column_stack([np.abs(on.theta[:half]) / edge, ratio[:half]])
    both = np.vstack([kept["x"], kept["y"]])
    check(f"{label}: energy at an angle over the filter off, less 1, at its largest inside 0.8 of the edge",
          float(np.abs(both[both[:, 0] < 0.8, 1] - 1).max()), BOUNDS["filter"]["kept"])
    check(f"{label}: energy at an angle over the filter off, at its largest outside 1.2 of the edge",
          float(both[both[:, 0] > 1.2, 1].max()), BOUNDS["filter"]["cut"])
    rec[on.name].update(edge=edge, x_edge=x_edge, kept_x=kept["x"].tolist(), kept_y=kept["y"].tolist())


def slice3(run, rec):
    """Slices of three wavelengths: a plane is three periods of the slice's radiation."""
    n, u = run.n_period, run.u
    f = lambda nu: run.spectrum(run.a0, nu)
    center, height = peak_of(f, -0.5 / n, 0.5 / n)
    width = width_of(f, center, height, 1.0 / n)
    label = "FLASH1, slices of 3 wavelengths"
    check(f"{label}: width of the line over 0.8859/N", width * n / ur.LINE_WIDTH,
          [1 - BOUNDS["slice3"]["width"], 1 + BOUNDS["slice3"]["width"]], "range")
    s_sum = abs(run.s[1:n + 1, run.loaded[0]].sum())
    ratio = float(run.level(center)[0] / (u.on_axis_level() * s_sum ** 2))
    check(f"{label}: level on the axis over the closed form", ratio,
          [1 - BOUNDS["slice3"]["level"], 1 + BOUNDS["slice3"]["level"]], "range")
    rec[run.name] = dict(width=width * n / ur.LINE_WIDTH, level=ratio)


def unaveraged(run, rec):
    """The unaveraged mode on the 55-period undulator."""
    n, u = run.n_period, run.u
    b = BOUNDS["unavg"]
    label = "FLASH1, 55 periods, unaveraged"
    f = lambda nu: run.spectrum(run.a0, nu)
    center, height = peak_of(f, -0.5 / n, 0.5 / n)
    width = width_of(f, center, height, 1.0 / n)
    check(f"{label}: width of the line over 0.8859/N", width * n / ur.LINE_WIDTH, b["width"], "range")
    nu = center + np.linspace(-4.0 / n, 4.0 / n, 801)
    measured = f(nu) / height
    check(f"{label}: largest difference from sinc^2 over four lobes",
          float(np.abs(measured - ur.line(n, nu - center)).max()), b["sinc"])
    # S is recorded at every step of the unaveraged mode. Its mean over the undulator
    # times the periods is the sum over periods of the averaged mode.
    s_sum = abs(run.s[1:, run.loaded[0]].mean()) * n
    ratio = float(run.level(center)[0] / (u.on_axis_level() * s_sum ** 2))
    check(f"{label}: level on the axis over the closed form", ratio, b["level"], "range")
    rec[run.name] = dict(n_period=n, helical=u.helical, nu=(nu - center).tolist(), line=measured.tolist(),
                         center=center, width=width * n / ur.LINE_WIDTH, level=ratio,
                         planes=(np.abs(run.a0) / np.abs(run.a0).max()).tolist(), first_plane=run.first_plane)
    with_angle(run, label, 0.14, b, rec)


def a_bunch(runs, wd, rec):
    """Many loaded slices with shot noise: each run against its own S, and the seeds' mean."""
    n, u = runs[0].n_period, runs[0].u
    nu = np.linspace(-4.0 / n, 4.0 / n, 161)
    worst, energy, mean = 0.0, [], np.zeros(len(nu))
    for r in runs:
        r.read(wd, keep_planes=False)
        (wd / f"{r.name}.stats.h5").unlink()                    # 47 MB a seed.
        pred = r.predicted()
        scale = np.vdot(pred, r.a0) / np.vdot(pred, pred)
        worst = max(worst, float(np.abs(r.a0 - scale * pred).max() / np.abs(r.a0).max()))
        # One electron's energy on the axis over the line is W1(0) N. The planes' sum of
        # |a0|^2 is that energy in the units of level(), by Parseval's theorem.
        k = 2 * np.pi / r.und["lambda0"]
        per_plane = ur.EPS0 * ur.C_LIGHT / 2 * k ** 2 * (r.cell * r.dz / ur.C_LIGHT) ** 2 / (2 * np.pi) ** 3
        energy.append(float(per_plane * (np.abs(r.a0) ** 2).sum() / (u.on_axis_level() * n * r.electrons)))
        mean += r.level(nu) / (u.on_axis_level() * n ** 2 * r.electrons) / len(runs)
    b = BOUNDS["bunch"]
    check("a bunch: planes that differ from the constant times the sum of conj(S), over the largest", worst, b["emitter"])
    check(f"a bunch: energy on the axis over one electron's times the electrons, mean of {len(runs)} seeds",
          float(np.mean(energy)), b["energy"], "range")
    check(f"a bunch: the same, the seeds' standard deviation", float(np.std(energy)), b["scatter"])
    rec["bunch"] = dict(n_period=n, slices=N_BUNCH, seeds=len(runs), electrons=runs[0].electrons,
                        nu=nu.tolist(), line=mean.tolist(), energy=energy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    exe = pathlib.Path(args.exe).resolve()
    wd = pathlib.Path(args.workdir).resolve()

    short = "use, SEG55"
    runs = dict(
        flash=Run("us_flash", FLASH, 165, 257),
        filter=Run("us_filter", FLASH, 165, 257, emit=1.5e-6, filter="T"),
        flash55=Run("us_flash55", FLASH, 55, 97, lines=short),
        delta=Run("us_delta", DELTA, 100, 161),
        size=Run("us_size", FLASH, 165, 257, emit=1.5e-6),
        spread=Run("us_spread", FLASH, 165, 257, grid=64, sig_pz=6e-4),
        window=Run("us_window", FLASH, 165, 65, grid=64),
        slice3=Run("us_slice3", FLASH, 165, 261, grid=64, n_wavelength=3),
        spectral=Run("us_spectral", FLASH, 55, 97, lines=short, extra='  propagator = "spectral"\n'),
        unavg=Run("us_unavg", FLASH, 55, 97, lines=short + "\nwiggler::*[FEL_METHOD] = unaveraged"),
    )
    bunch = [Run(f"us_bunch{k}", FLASH, 55, N_BUNCH + 64, grid=32, n_loaded=N_BUNCH, n_particle=64,
                 seed=1000 + k, lines=short) for k in range(SEEDS)]
    pool.run_all([lambda r=r: r.start(exe, wd) for r in list(runs.values()) + bunch])
    rec = {}

    # The line on the axis, the level, the ridge.
    for key, label, x_max in (("flash", "FLASH1, 165 periods", 0.14), ("delta", "Delta, 100 periods", 0.20)):
        r = runs[key]
        r.read(wd)
        on_axis(r, label, TIGHT, rec)
        with_angle(r, label, x_max, BOUNDS[key], rec)
    r = runs["flash55"]
    r.read(wd, keep_planes=False)
    on_axis(r, "FLASH1, 55 periods", TIGHT, rec)
    check("FLASH1: width of the line at 55 periods over that at 165",
          rec["us_flash55"]["width"] * 165 / (rec["us_flash"]["width"] * 55), [2.99, 3.01], "range")

    # What the beam and the window add.
    for key, label in (("size", "FLASH1, emittance 1.5 um"), ("spread", "FLASH1, energy spread 6e-4"),
                       ("window", "FLASH1, window of 65 wavelengths")):
        r = runs[key]
        r.read(wd, keep_planes=(key == "size"))
        on_axis(r, label, BOUNDS[key], rec)
    r = runs["slice3"]
    r.read(wd, keep_planes=False)
    slice3(r, rec)

    # The source filter, on the beam of FLASH1's emittance, whose shot noise radiates
    # into every angle the grid holds.
    with_angle(runs["size"], "FLASH1, emittance 1.5 um", 0.14, BOUNDS["size"], rec, deposit=False)
    r = runs["filter"]
    r.read(wd)
    source_filter(r, runs["size"], wd, rec)

    # The spectral propagator, in which each frequency diffracts at its own wavenumber.
    r = runs["spectral"]
    r.read(wd)
    label = "FLASH1, 55 periods, field.propagator = spectral"
    on_axis(r, label, BOUNDS["spectral"], rec)
    with_angle(r, label, 0.08, BOUNDS["spectral"], rec, deposit=False, spectral=True)

    # The unaveraged mode.
    r = runs["unavg"]
    r.read(wd)
    unaveraged(r, rec)

    # The bunch.
    a_bunch(bunch, wd, rec)

    (wd / "undulator_spectrum.json").write_text(json.dumps(rec))
    if FAILED:
        print("FAIL: undulator-spectrum checks")
        sys.exit(1)
    print("PASS: undulator-spectrum checks")


if __name__ == "__main__":
    main()
