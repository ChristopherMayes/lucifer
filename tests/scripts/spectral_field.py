#!/usr/bin/env python3
"""
Measure the spectral field propagator, field.propagator = "spectral", against the slice
propagator, as docs/spectral-field.md records. The keystone's own check of the mode is
check_spectral.py; nothing here runs in it, and the page is hand written from this
script's report.

    spectral_field.py --exe <lucifer> --work-dir <dir> [--validation-dir <dir>]
                      [--threads N] [--reuse] [vacuum] [null] [bank] [cost]
                      [seeded] [chirped]

The six measurements, each a subcommand and the first four when none is named. With
--reuse a run whose log and final field exist in the work directory is not run again,
so a set that takes hours is re-analyzed from its files.

  vacuum    A chirped Gaussian pulse with no beam, drifted through a pipe under both
            propagators, against the exact paraxial solution with its spatio-temporal
            coupling: each frequency's Gaussian beam diffracts at its own Rayleigh
            range, and the solution is the superposition over the pulse's spectrum. The
            error of each propagator as a function of the pulse's relative bandwidth.
            Needs the tracker alone.
  null      The steady-state tier and the SASE tier under "spectral", compared against
            Genesis4 at the tier tolerances as compare_fel.py compares them and against
            the slice run, with the SASE tier again under both propagators with the
            source filter on. Needs the work directory of a completed validation pass
            (--validation-dir), which holds the tier inputs, the imported dumps and the
            Genesis4 references.
  bank      The SASE tier's window, 32 slices under a line that slips 1184, with the
            escaped field kept under both propagators: the exit energies, the escaped
            energy, the beam's energy loss, and what the pulse file holds. Needs
            --validation-dir.
  cost      The field phase of the SASE tier's walk per step, from the footer's timing
            table, under both propagators, at the thread count given. Needs
            --validation-dir.
  seeded    A chirped Gaussian seed, the vacuum test's pulse at the tier's wavelength on
            the tier's grid, through the first Aramis segment and through the full line
            with the tier's beam continued over a 1280-slice window that holds the line's
            slippage, under both propagators with the source filter on, at four relative
            bandwidths, and the segment again on the bunch loaded without its shot noise,
            the seed the only light. The exit field differenced after the slice run's
            field is moved by the slippage residual it carries, the exit power and the
            per-slice wavelength of each. The spectral run is the reference, and its step
            cut from three periods to one says how far it moves. Needs --validation-dir
            for the lattice.
  chirped   The same beam with a linear energy chirp along the window, from noise
            through the full line under both propagators, unchirped and at three chirps.
            The same three numbers and the slope the exit wavelengths carry against minus
            twice the energy chirp. Needs --validation-dir for the lattice.
  sweep     The field's sampling along tau, field.n_wavelength, at 1, 2, 5, 10, 20 and
            40 slices a plane on four lines from noise under the spectral propagator:
            the SASE tier's beam over 1280 slices, and the sase, flash1 and lcls
            examples. For each run the window's and the pulse's energy, the peak power,
            the exit power per slice against the run at one slice a plane, the exit
            spectrum's rms width and centroid, and the walk. Needs --validation-dir for
            the Aramis lattice.

The tracker is run from the work directory, into which the validation directory's
inputs are linked, so the validation directory is read and never written.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# --------------------------------------------------------------------------------------
# Running the tracker.

def run(exe, toml, overrides, out_root, cwd, threads, log=None):
    """Run the tracker on toml with overrides, out_root set, and return the log's text.
    A failed run returns the text too, and the caller reads the refusal from it."""
    log_path = pathlib.Path(cwd) / (log or f"{out_root}.log")
    if REUSE and log_path.exists() and (pathlib.Path(cwd) / f"{out_root}-final.wf.h5").exists():
        text = log_path.read_text()
        if re.search(r"\n Exit\s+power", text):
            return 0, text
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = [str(exe), str(toml), f"run.out_root={out_root}"] + list(overrides)
    p = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True)
    text = p.stdout + p.stderr
    (pathlib.Path(cwd) / (log or f"{out_root}.log")).write_text(text)
    return p.returncode, text


def footer(text):
    """The footer's exit and escaped energies [J], the field-solve row (regions, seconds)
    and the walk's seconds, from a run's log."""
    out = {}
    m = re.search(r"Exit\s+power\s+(\S+)\s+([A-Za-z]+),\s+window energy\s+(\S+)\s+([A-Za-z]+),\s+pulse energy\s+(\S+)\s+([A-Za-z]+)", text)
    if m:
        out["exit_power"] = si(m.group(1), m.group(2))
        out["window_energy"] = si(m.group(3), m.group(4))
        out["pulse_energy"] = si(m.group(5), m.group(6))
    m = re.search(r"Escaped\s+h1\s+(\S+)\s+J", text)
    if m:
        out["escaped_energy"] = float(m.group(1))
    m = re.search(r"field solve\s+(\d+)\s+([\d.]+)\s+([\d.]+)%", text)
    if m:
        out["solve_regions"], out["solve_seconds"] = int(m.group(1)), float(m.group(2))
    m = re.search(r"\n\s+walk\s+(\d+)\s+([\d.]+)\s+", text)
    if m:
        out["walk_seconds"] = float(m.group(2))
    m = re.search(r"Work\s+(\d+) FEL steps", text)
    if m:
        out["n_step"] = int(m.group(1))
    return out


SI = {"": 1, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15, "a": 1e-18,
      "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}


def si(value, unit):
    """A value printed with an SI prefix on its unit, as the footer prints it."""
    prefix = unit[:-1] if len(unit) > 1 else ""
    return float(value) * SI[prefix]


def link_inputs(vdir, wdir):
    """Link the validation directory's inputs into the work directory: the input files,
    the lattices and the initial dumps both codes start from."""
    for p in pathlib.Path(vdir).iterdir():
        if p.suffix in (".toml", ".bmad", ".lat") or "-initial" in p.name:
            dst = pathlib.Path(wdir) / p.name
            if not dst.exists():
                dst.symlink_to(p.resolve())


# --------------------------------------------------------------------------------------
# The vacuum test.

VACUUM_LATTICE = """! The vacuum test's line (spectral_field.py): one short planar undulator, so the
! lattice has an FEL element, and the pipe the pulse drifts through. The walk tracks
! the pipe alone (run.track_start and track_end). The energy is high enough that the
! pipe's slippage is Genesis4's floor of one wavelength, one slice at n_wavelength = 1,
! which both propagators apply as a whole slice and the exact solution includes.

no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 100e9

beginning[beta_a] = 10
beginning[beta_b] = 10

UND: wiggler, l = 0.03, l_period = 0.015, field_calc = planar_model, &
     b_max = 0.5 * sqrt(2) * (twopi / 0.015) * m_electron / c_light, &
     fel_method = averaged, ds_step = 0.015
DRIFT1: pipe, l = {length}

VAC: line = (UND, DRIFT1)

use, VAC
"""

VACUUM_INPUT = """lattice = "vacuum.bmad"
track_start = "DRIFT1"
track_end = "DRIFT1"
source_filter = F
keep_escaped_field = F
n_wavelength = 1
beam.file = "{beam}"
field.file = "{field}"
lambda0 = {wavelength}
"""


class _Patch:
    """One slice's particles in the form convert_genesis.write_openpmd_beam reads: the
    attributes of a ParticleGroup it uses, and nothing else."""

    def __init__(self, n, k, spacing, pz, weight):
        self.n_particle = n
        self.x = np.zeros(n)
        self.y = np.zeros(n)
        self.px = np.zeros(n)
        self.py = np.zeros(n)
        self.pz = np.full(n, pz)
        self.z = k * spacing + spacing * (np.arange(n) + 0.5) / n
        self.weight = np.full(n, weight)


def vacuum_beam(path, nz, spacing, e_tot):
    """A beam of no consequence, four particles a slice at a weight of 1e-21 C, so the
    walk has a bunch to track while the field drifts: 'no beam' as near as an input file
    can state it, since a field file needs a beam file beside it."""
    from convert_genesis import write_openpmd_beam
    groups = [_Patch(4, k, spacing, e_tot, 1e-21) for k in range(nz)]
    write_openpmd_beam(path, groups, spacing)


def vacuum(args, wdir):
    from beamphysics.wavefront import Wavefront
    from lucifer.field import read_field
    from lucifer.flat import to_toml

    wavelength = 1e-9
    nz, ngrid = 1024, 128
    dtau = wavelength                     # n_wavelength = 1
    sigma = 50 * dtau                     # The pulse's rms length.
    w0 = 20e-6
    half = 200e-6
    k0 = 2 * np.pi / wavelength
    z_r = k0 * w0**2 / 2
    length = 2 * z_r
    dx = 2 * half / (ngrid - 1)
    x = (np.arange(ngrid) - (ngrid - 1) / 2) * dx
    tau = np.arange(nz) * dtau
    tau_c = (nz // 2) * dtau
    r2 = (x[:, None]**2 + x[None, :]**2)

    (wdir / "vacuum.bmad").write_text(VACUUM_LATTICE.format(length=repr(length)))
    vacuum_beam(wdir / "vacuum-initial.beam.h5", nz, dtau, 100e9)
    print("--- vacuum: a chirped Gaussian pulse through a pipe, against the exact solution")
    print(f"  wavelength {wavelength:.1e} m, {nz} slices of one wavelength, pulse sigma {sigma/dtau:.0f} slices,")
    print(f"  waist {w0*1e6:.0f} um, grid {ngrid} points at {dx*1e6:.2f} um, pipe {length:.4f} m = 2 z_R")
    print()
    print("  bandwidth   chirp [1/m^2]    slices        spectral")

    # The relative rms bandwidth B of exp(-xi^2/4 sigma^2 + i b xi^2) is
    # sqrt(1 + 16 b^2 sigma^4) / (2 sigma k0), so a target above the unchirped
    # 1 / (2 sigma k0) fixes b.
    b_unchirped = 1 / (2 * sigma * k0)
    rows = []
    for target in (b_unchirped, 5e-3, 1e-2, 2e-2, 5e-2):
        if target <= b_unchirped * 1.0000001:
            b = 0.0
        else:
            b = np.sqrt((2 * sigma * k0 * target)**2 - 1) / (4 * sigma**2)
        a = 1 / (4 * sigma**2) - 1j * b
        env = np.exp(-a * (tau - tau_c)**2)
        ex = np.exp(-r2 / w0**2)[:, :, None] * env[None, None, :]
        tag = f"chirp{len(rows)}"
        wf = Wavefront(Ex=ex.astype(np.complex128), dx=dx, dy=dx, dz=dtau, wavelength=wavelength)
        wf.write_openpmd(str(wdir / f"{tag}-initial.wf.h5"))
        flat = VACUUM_INPUT.format(beam="vacuum-initial.beam.h5", field=f"{tag}-initial.wf.h5",
                                   wavelength=repr(wavelength))
        (wdir / f"{tag}.toml").write_text(to_toml(flat))

        errors = {}
        for name, over in (("slices", []), ("spectral", ["field.propagator=spectral"])):
            rc, text = run(args.exe, f"{tag}.toml", over, f"{tag}_{name}", wdir, args.threads)
            if rc != 0:
                print(text[-2000:])
                raise SystemExit(f"vacuum: {tag} {name} failed")
            fld = read_field(wdir / f"{tag}_{name}-final.wf.h5")
            u = np.transpose(fld["u"], (2, 1, 0))          # (nx, ny, nz)
            exact = vacuum_exact(a, tau, tau_c, dtau, r2, k0, w0, length, shift=dtau)
            errors[name] = np.abs(u - exact).max() / np.abs(exact).max()
        rows.append((target, b, errors["slices"], errors["spectral"]))
        print(f"  {target:9.1e}   {b:13.4e}    {errors['slices']:.3e}     {errors['spectral']:.3e}")
    print()
    return rows


def vacuum_exact(a, tau, tau_c, dtau, r2, k0, w0, z, shift):
    """The exact paraxial solution at z of exp(-r^2/w0^2) exp(-a (tau - tau_c)^2): the
    pulse's analytic spectrum on an eightfold finer wavenumber grid than the window
    carries, which is the spectrum of a window eight times longer, each component the
    Gaussian beam at its own k = k0 + kappa, superposed by an inverse transform and
    read at the window's slices. The window-frame shift the schedule charges the pipe
    is applied as the slice propagator applies it, one whole slice."""
    nz = len(tau)
    n_fine = 8 * nz
    m = np.fft.fftfreq(n_fine, d=1.0) * n_fine         # -n/2 .. n/2-1 in FFT order
    kappa = 2 * np.pi * m / (n_fine * dtau)
    spectrum = np.sqrt(np.pi / a) * np.exp(-kappa**2 / (4 * a)) * np.exp(-1j * kappa * tau_c)
    spectrum = spectrum * np.exp(-1j * kappa * shift)
    k = k0 + kappa
    z_r = k * w0**2 / 2
    q = 1 + 1j * z / z_r                                   # (n_fine,)
    # E(r, tau) = (1/2 pi) sum_m dkappa spectrum_m g_m(r) exp(i kappa_m tau); on the
    # fine grid the sum is n_fine * ifft of (spectrum * g) times dkappa / 2 pi.
    dkappa = 2 * np.pi / (n_fine * dtau)
    out = np.empty(r2.shape + (nz,), dtype=np.complex128)
    flat_r2 = r2.reshape(-1)
    chunk = 256
    for i0 in range(0, flat_r2.size, chunk):
        rr = flat_r2[i0:i0 + chunk]
        g = np.exp(-rr[:, None] / (w0**2 * q[None, :])) / q[None, :]
        line = np.fft.ifft(spectrum[None, :] * g, axis=1) * n_fine * dkappa / (2 * np.pi)
        out.reshape(-1, nz)[i0:i0 + chunk] = line[:, :nz]
    return out


# --------------------------------------------------------------------------------------
# The null test, the bank and the cost, on the validation directory's tiers.

def tier_inputs(args, wdir):
    vdir = pathlib.Path(args.validation_dir).resolve()
    if not (vdir / "tdsase.toml").exists():
        raise SystemExit(f"{vdir} is not a completed validation pass: no tdsase.toml")
    link_inputs(vdir, wdir)
    return vdir


def null(args, wdir):
    import compare_fel
    from lucifer.compare import differences, read_curves

    vdir = tier_inputs(args, wdir)
    print("--- null: the tiers under the spectral propagator, against Genesis4 and the slice run")
    nslice_sase = compare_fel.load_nslice(f"{vdir}/AramisTDSASE-initial.par.h5")
    tiers = (
        ("tier1", "Aramis1seg", 1, compare_fel.TOLERANCE["tier1"]),
        ("tdsase", "AramisTDSASE", nslice_sase, compare_fel.TOLERANCE["tdsase"]),
    )
    results = {}
    for tier, gen, nslice, tol in tiers:
        slices_log = (vdir / f"fel-{tier}.log").read_text()
        recorded = compare_fel.compare_tier(f"{tier} slices (the validation pass's own run)",
                                            f"{vdir}/{tier}.diag.txt", f"{vdir}/{gen}.out.h5",
                                            f"{vdir}/{tier}-final.wf.h5", f"{vdir}/{gen}-final.fld.h5",
                                            f"{vdir}/{tier}-final.beam.h5", f"{vdir}/{gen}-final.par.h5",
                                            tol, nslice)[0]
        results[(tier, "slices")] = recorded
        root = f"{tier}_spectral"
        rc, text = run(args.exe, f"{tier}.toml", ["field.propagator=spectral"], root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"null: {root} failed")
        worst = compare_fel.compare_tier(f"{tier} spectral",
                                         f"{wdir}/{root}.diag.txt", f"{vdir}/{gen}.out.h5",
                                         f"{wdir}/{root}-final.wf.h5", f"{vdir}/{gen}-final.fld.h5",
                                         f"{wdir}/{root}-final.beam.h5", f"{vdir}/{gen}-final.par.h5",
                                         tol, nslice)[0]
        d = differences(read_curves(f"{wdir}/{root}.stats.h5"), read_curves(f"{vdir}/{tier}.stats.h5"))
        own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
        results[(tier, "spectral")] = (worst, own)
        print(f"  {tier} spectral: against Genesis4 {worst:.6e} (tolerance {tol:.1e}), "
              f"against the slice run: " + ", ".join(f"{k} {v:.2e}" for k, v in own.items()))
        print()
    # The tiers run with the source filter off, since Genesis4 is compared with it off.
    # The filter is the production default, and it removes the wide-angle source, so the
    # SASE tier reruns with it on under both propagators, Lucifer against Lucifer.
    print("  the SASE tier with the source filter on, against the slice run with it on:")
    rc, text = run(args.exe, "tdsase.toml", ["run.source_filter=true"], "tdsase_filter_slices", wdir, args.threads)
    if rc != 0:
        print(text[-2000:])
        raise SystemExit("null: tdsase_filter_slices failed")
    ref = read_curves(f"{wdir}/tdsase_filter_slices.stats.h5")
    root = "tdsase_filter_spectral"
    rc, text = run(args.exe, "tdsase.toml", ["run.source_filter=true", "field.propagator=spectral"],
                   root, wdir, args.threads)
    if rc != 0:
        print(text[-2000:])
        raise SystemExit(f"null: {root} failed")
    d = differences(read_curves(f"{wdir}/{root}.stats.h5"), ref)
    own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
    results[("tdsase_filter", "spectral")] = own
    print(f"    spectral: " + ", ".join(f"{k} {v:.2e}" for k, v in own.items()))
    print()
    return results


def bank(args, wdir):
    from lucifer.stats import read_stats

    tier_inputs(args, wdir)
    print("--- bank: the SASE tier's 32-slice window with the escaped field kept, under both propagators")
    rows = {}
    for key, label, over in (("slices", "the slice propagator", []),
                             ("spectral", "the spectral propagator", ["field.propagator=spectral"])):
        root = f"bank_{key}"
        rc, text = run(args.exe, "tdsase.toml", over + ["run.keep_escaped_field=true"], root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"bank: {root} failed")
        row = dict(label=label, **footer(text))
        st = read_stats(wdir / f"{root}.stats.h5")
        row.update(closure(st))
        pulse = wdir / f"{root}-pulse.wf.h5"
        row["pulse_planes"] = 0
        if pulse.exists():
            import h5py
            with h5py.File(pulse) as h:
                row["pulse_planes"] = h["data/1/meshes/electricField/x"].shape[0]
        rows[key] = row
        print(f"  {label}:")
        print(f"    exit window energy {row['window_energy']:.4e} J, escaped {row.get('escaped_energy', 0):.4e} J, "
              f"their sum {row['window_energy'] + row.get('escaped_energy', 0):.4e} J, pulse {row['pulse_energy']:.4e} J")
        print(f"    beam energy loss {row['beam_loss']:.4e} J, pulse file {row['pulse_planes']} planes")
        print(f"    tail slice's exit power share {row['tail_share']:.2e}, head slice's {row['head_share']:.2e}")
    print()
    return rows


def closure(st):
    """The beam's energy loss over the walk, each slice's live charge times its mean
    energy change, from the stats file, with the tail and head slices' shares of the
    exit power. The averaged equations charge the beam for the field at the particles
    and not for the incoherent emission the deposit adds, so on a run from noise the
    loss is far below the field's energy and the two do not close."""
    charge = np.asarray(st["beam/slice/charge_live"])        # (nz, ns) C
    energy = np.asarray(st["beam/slice/energy"])             # (nz, ns) eV
    beam_loss = float(np.sum(charge[0] * (energy[0] - energy[-1])))   # eV * C = J
    power = np.asarray(st["field/total/power"])
    tot = power[-1].sum()
    return dict(beam_loss=beam_loss,
                tail_share=float(power[-1, 0] / tot) if tot else float("nan"),
                head_share=float(power[-1, -1] / tot) if tot else float("nan"))


def cost(args, wdir):
    tier_inputs(args, wdir)
    print(f"--- cost: the SASE tier's field phase per step, {args.threads} thread(s)")
    rows = {}
    # The full cost under spectral: the kernel pass and the slip pass, the bank off as the
    # tier runs it.
    for key, over in (("slices", []), ("spectral", ["field.propagator=spectral"])):
        root = f"cost_{key}"
        rc, text = run(args.exe, "tdsase.toml", over, root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"cost: {root} failed")
        f = footer(text)
        rows[key] = f
        print(f"  {key:9s} field solve {f['solve_seconds']:.3f} s over {f['n_step']} steps = "
              f"{1e3 * f['solve_seconds'] / f['n_step']:.2f} ms a step, walk {f['walk_seconds']:.3f} s")
    print()
    return rows


# --------------------------------------------------------------------------------------
# The chirped cases: a long flat beam at the SASE tier's parameters, a seed with a chirp
# or a beam with one, through the line under both propagators.

CHIRP_LATTICE = """! The chirped cases' line (spectral_field.py): the Aramis line of the SASE tier behind a
! marker, so a run tracking the marker alone dumps the beam and the dark field the seeded
! runs start from, on the window and grid every run shares.
call, file = aramis.bmad
MARK: marker
CHIRP: line = (MARK, ARAMIS)
use, CHIRP
"""

CHIRP_LATTICE_FINE = """! The chirped cases' line with the undulator step halved to one period, for the
! convergence of the spectral run (spectral_field.py).
call, file = aramis.bmad
wiggler::*[ds_step] = 0.015
MARK: marker
CHIRP: line = (MARK, ARAMIS)
use, CHIRP
"""

# The SASE tier's beam (Aramis-td-sase.in: 3000 A, delgam 1, 4e-7 m emittances), flat
# over a window of 1280 slices of 3 wavelengths, which holds the seed and the line's
# slippage of 1184 slices: the seed sits 30 slices from the tail and reaches slice 1214 at
# the exit. 1024 particles a slice in place of the tier's 2048, since the seeded runs
# stand well above the shot noise.
PREP_INPUT = """lattice = "chirp-line.bmad"
track_start = "MARK"
track_end = "MARK"
keep_escaped_field = F
n_wavelength = 3
window_length = 3.84e-7
window.current = 3000
n_particle = 1024
beamlet_size = 8
sig_pz = 8.8045e-5
a_norm_emit = 4e-7
b_norm_emit = 4e-7
quiet_start = T
shot_noise = T
lambda0 = 1e-10
grid_n_pts = 255
grid_half_width = 2e-4
seed_power = 0
"""

CHIRP_RUN_INPUT = """lattice = "{lattice}"
keep_escaped_field = F
n_wavelength = 3
beam.file = "{beam}"
field.file = "{field}"
lambda0 = 1e-10
"""

SEED_SIGMA_WAVELENGTHS = 15      # The seed's rms length. Its unchirped rms bandwidth is 0.53 percent.
SEED_SLICE = 30                  # Where the seed's center sits, from the tail.
SEED_WAIST = 30e-6               # Its 1/e^2 intensity radius [m], near the mode's.
SEED_POWER = 1e6                 # Its peak power [W], a thousand times the tier's startup noise.
BANDWIDTHS = (0.005, 0.01, 0.02, 0.05)
ENERGY_CHIRPS = (0.0, 0.0025, 0.005, 0.01)   # the first is the unchirped reference
REUSE = False          # --reuse: a run whose log and final field exist is not run again


def chirp_prepare(args, wdir):
    """The line, the beam file and the dark field the chirped cases start from, made once."""
    from lucifer.flat import to_toml
    tier_inputs(args, wdir)
    (wdir / "chirp-line.bmad").write_text(CHIRP_LATTICE)
    (wdir / "chirp-line-fine.bmad").write_text(CHIRP_LATTICE_FINE)
    if not (wdir / "prep-final.beam.h5").exists():
        (wdir / "prep.toml").write_text(to_toml(PREP_INPUT))
        rc, text = run(args.exe, "prep.toml", [], "prep", wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit("chirp: the preparatory run failed")
    return wdir / "prep-final.beam.h5", wdir / "prep-final.wf.h5"


def quiet_beam(args, wdir):
    """The same bunch loaded without its shot noise, so a seed is the only light and the
    propagators' difference on the seed is not under the start-up noise light's."""
    if not (wdir / "prep-quiet-final.beam.h5").exists():
        rc, text = run(args.exe, "prep.toml", ["beam.shot_noise=false"], "prep-quiet", wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit("chirp: the quiet preparatory run failed")
    return wdir / "prep-quiet-final.beam.h5"


def chirp_seed(wdir, dark, tag, bandwidth):
    """The chirped Gaussian seed on the dark field's grid and window, at the relative rms
    bandwidth asked, scaled to SEED_POWER at its peak."""
    from beamphysics.wavefront import Wavefront
    from lucifer.field import read_field
    f = read_field(dark)
    nz, ny, nx = f["u"].shape
    dx, dtau, lam = f["dx"], f["dz"], f["wavelength"]
    k0 = 2 * np.pi / lam
    sigma = SEED_SIGMA_WAVELENGTHS * lam
    x = (np.arange(nx) - (nx - 1) / 2) * dx
    r2 = x[:, None]**2 + x[None, :]**2
    tau = np.arange(nz) * dtau
    tau_c = SEED_SLICE * dtau
    b_unchirped = 1 / (2 * sigma * k0)
    b = 0.0 if bandwidth <= b_unchirped else np.sqrt((2 * sigma * k0 * bandwidth)**2 - 1) / (4 * sigma**2)
    a = 1 / (4 * sigma**2) - 1j * b
    g = np.exp(-r2 / SEED_WAIST**2)
    z0 = 376.730313668
    e0 = np.sqrt(2 * z0 * SEED_POWER / ((g**2).sum() * dx * dx))
    ex = (e0 * g)[:, :, None] * np.exp(-a * (tau - tau_c)**2)[None, None, :]
    path = wdir / f"{tag}-initial.wf.h5"
    Wavefront(Ex=ex.astype(np.complex128), dx=dx, dy=dx, dz=dtau, wavelength=lam).write_openpmd(str(path))
    return path, b, max(bandwidth, b_unchirped)


def chirp_beam(src, dst, chirp):
    """A copy of the beam file with a linear energy chirp along the window, chirp head to
    tail, the head at the high slices and above the tail in energy: every momentum
    component scaled by the same factor, so the angles stay. Bmad's reader takes the
    energy from totalMomentum, the momentum deviation from the reference p0c, when the
    file carries it, so that record is scaled as the whole momentum and not the
    deviation."""
    import h5py
    shutil.copyfile(src, dst)
    with h5py.File(dst, "r+") as h:
        g = h["data/00001/particles/electron"]
        n = np.asarray(g["particlePatches/numParticles"])
        off = np.asarray(g["particlePatches/numParticlesOffset"])
        nz = n.size
        slice_of = np.empty(int(n.sum()), dtype=int)
        for k in range(nz):
            slice_of[off[k]:off[k] + n[k]] = k
        factor = 1 + chirp * (slice_of - (nz - 1) / 2) / (nz - 1)
        for c in ("x", "y", "z"):
            d = g[f"momentum/{c}"]
            d[...] = np.asarray(d) * factor
        p0c = float(np.asarray(g["totalMomentumOffset"].attrs["value"]).ravel()[0])
        d = g["totalMomentum"]
        d[...] = (np.asarray(d) + p0c) * factor - p0c


def slice_wavelengths(u, dtau, lam, floor=1e-3):
    """Each slice's wavelength from the phase advance to its neighbor, 2 pi / (k0 + kappa)
    with kappa from arg(sum u_i conj u_{i+1}), on the slices above floor of the peak
    power; NaN elsewhere."""
    k0 = 2 * np.pi / lam
    power = (np.abs(u)**2).sum(axis=(1, 2))
    lit = power > floor * power.max()
    c = (u[:-1] * np.conj(u[1:])).sum(axis=(1, 2))
    kappa = -np.angle(c) / dtau
    out = np.full(u.shape[0], np.nan)
    ok = lit[:-1] & lit[1:]
    out[:-1][ok] = 2 * np.pi / (k0 + kappa[ok])
    return out, power


def slippage_residual(path):
    """The slippage a run's final field has not yet applied [wavelengths], from the
    file's lucifer/field group; zero under the phase ramp, which leaves none."""
    import h5py
    with h5py.File(path) as h:
        return float(np.asarray(h["lucifer/field"].attrs["slippageResidual"]).ravel()[0])


def shift_slices(u, d):
    """u moved d slices toward the head along tau, by the transform's phase: the slice
    propagator's field carries its residual as a lag, and the comparison removes it."""
    nz = u.shape[0]
    m = np.fft.fftfreq(nz) * nz
    return np.fft.ifft(np.fft.fft(u, axis=0) * np.exp(-2j * np.pi * m * d / nz)[:, None, None], axis=0)


def chirp_compare(wdir, root_a, root_b, floor=0.25):
    """The exit of run a against run b: the field difference over b's peak, over the whole
    window as the runs left it and over the lit slices after a is moved by the slippage
    residual the two runs differ by, each exit power,
    and the largest relative wavelength difference over the slices above floor of the
    peak power. The wavelengths at 1e-3 of the peak are returned for the chirp fit."""
    from lucifer.field import read_field
    fa, fb = read_field(wdir / f"{root_a}-final.wf.h5"), read_field(wdir / f"{root_b}-final.wf.h5")
    ua, ub = fa["u"], fb["u"]
    dtau, lam = fb["dz"], fb["wavelength"]
    n_wavelength = int(round(dtau / lam))
    offset = (slippage_residual(wdir / f"{root_a}-final.wf.h5")
              - slippage_residual(wdir / f"{root_b}-final.wf.h5")) / n_wavelength
    peak = np.abs(ub).max()
    field_raw = float(np.abs(ua - ub).max() / peak)
    if offset != 0:
        ua = shift_slices(ua, offset)
    la, pa = slice_wavelengths(ua, dtau, lam, floor)
    lb, pb = slice_wavelengths(ub, dtau, lam, floor)
    # The aligned difference is taken over the lit slices and four on each side, since the
    # shift's band-limited form rings where the start-up noise light meets the window's ends.
    lit = np.where(pb > floor * pb.max())[0]
    lo, hi = max(lit.min() - 4, 0), min(lit.max() + 5, ub.shape[0])
    z0 = 376.730313668
    ok = np.isfinite(la) & np.isfinite(lb)
    return dict(field=float(np.abs(ua - ub)[lo:hi].max() / peak), field_raw=field_raw, offset=float(offset),
                power_a=float(pa.max() * fb["dx"]**2 / (2 * z0)), power_b=float(pb.max() * fb["dx"]**2 / (2 * z0)),
                wavelength=float(np.abs(la[ok] / lb[ok] - 1).max()) if ok.any() else float("nan"),
                lit=int(ok.sum()), lam_a=slice_wavelengths(ua, dtau, lam, 1e-3)[0],
                lam_b=slice_wavelengths(ub, dtau, lam, 1e-3)[0])


def seeded(args, wdir):
    from lucifer.flat import to_toml
    beam, dark = chirp_prepare(args, wdir)
    quiet = quiet_beam(args, wdir)
    print("--- seeded: a chirped seed through Aramis under both propagators, source filter on")
    print(f"  seed sigma {SEED_SIGMA_WAVELENGTHS} wavelengths at slice {SEED_SLICE}, waist {SEED_WAIST*1e6:.0f} um, "
          f"peak {SEED_POWER:.0e} W; 1280 slices of 3 wavelengths, 1024 particles a slice, {args.threads} threads")
    print()
    print("  bandwidth  extent   offset   field diff raw / aligned   wavelength diff   lit    peak power slices / spectral")
    rows = []
    for i, bw in enumerate(BANDWIDTHS):
        tag = f"seed{i}"
        seed, b, bw_eff = chirp_seed(wdir, dark, tag, bw)
        (wdir / f"{tag}.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad", beam=beam.name, field=seed.name)))
        (wdir / f"{tag}_quiet.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad", beam=quiet.name, field=seed.name)))
        # quiet: the segment with the beam's shot noise off, the seed the only light.
        for extent, over in (("quiet", ["run.track_end=UND##1"]), ("segment", ["run.track_end=UND##1"]), ("line", [])):
            toml = f"{tag}_quiet.toml" if extent == "quiet" else f"{tag}.toml"
            roots = {}
            for prop, more in (("slices", []), ("spectral", ["field.propagator=spectral"])):
                root = f"{tag}_{extent}_{prop}"
                rc, text = run(args.exe, toml, over + more, root, wdir, args.threads)
                if rc != 0:
                    print(text[-2000:])
                    raise SystemExit(f"seeded: {root} failed")
                roots[prop] = root
            r = chirp_compare(wdir, roots["slices"], roots["spectral"])
            rows.append((bw_eff, extent, r))
            print(f"  {bw_eff:8.3f}   {extent:7s}  {r['offset']:6.3f}   {r['field_raw']:.3e} / {r['field']:.3e}    "
                  f"{r['wavelength']:.3e}       {r['lit']:4d}   {r['power_a']:.3e} / {r['power_b']:.3e} W")
    # Convergence of the reference: the 2 percent seed with the undulator step cut from
    # three periods to one, through the segment and through the line.
    tag = "seed2"
    (wdir / f"{tag}_fine.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line-fine.bmad", beam=beam.name,
                                                                           field=f"{tag}-initial.wf.h5")))
    for extent, over in (("segment", ["run.track_end=UND##1"]), ("line", [])):
        rc, text = run(args.exe, f"{tag}_fine.toml", over + ["field.propagator=spectral"],
                       f"{tag}_{extent}_spectral_fine", wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"seeded: the fine-step {extent} run failed")
        r = chirp_compare(wdir, f"{tag}_{extent}_spectral_fine", f"{tag}_{extent}_spectral")
        print(f"  convergence, {extent}: the 2 percent seed, spectral at ds_step 0.015 against 0.045: "
              f"field {r['field']:.3e}, wavelength {r['wavelength']:.3e}, peak power {r['power_a']:.3e} / {r['power_b']:.3e} W")
        r = chirp_compare(wdir, f"{tag}_{extent}_slices", f"{tag}_{extent}_spectral_fine")
        print(f"               and the slice run against that fine spectral run: field {r['field']:.3e} aligned by {r['offset']:.3f} slice")
    print()
    return rows


def chirped(args, wdir):
    from lucifer.flat import to_toml
    beam, dark = chirp_prepare(args, wdir)
    print("--- chirped: an energy-chirped beam from noise through Aramis under both propagators, source filter on")
    print(f"  1280 slices of 3 wavelengths, 1024 particles a slice, {args.threads} threads")
    print()
    print("  chirp    field diff   wavelength diff   lit   exit peak power slices / spectral   fitted slope of the wavelengths slices / spectral")
    rows = []
    for i, chirp in enumerate(ENERGY_CHIRPS):
        tag = f"chirp{i}"
        chirp_beam(beam, wdir / f"{tag}-initial.beam.h5", chirp)
        (wdir / f"{tag}.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad",
                                                                          beam=f"{tag}-initial.beam.h5", field=dark.name)))
        roots = {}
        for prop, more in (("slices", []), ("spectral", ["field.propagator=spectral"])):
            root = f"{tag}_line_{prop}"
            rc, text = run(args.exe, f"{tag}.toml", more, root, wdir, args.threads)
            if rc != 0:
                print(text[-2000:])
                raise SystemExit(f"chirped: {root} failed")
            roots[prop] = root
        r = chirp_compare(wdir, roots["slices"], roots["spectral"])
        fits = []
        for lam_i in (r["lam_a"], r["lam_b"]):
            ok = np.isfinite(lam_i)
            idx = np.arange(lam_i.size)[ok]
            slope = np.polyfit(idx, lam_i[ok] / np.nanmean(lam_i), 1)[0] * (lam_i.size - 1)
            fits.append(slope)
        rows.append((chirp, r, fits))
        print(f"  {chirp:6.4f}   {r['field']:.3e}    {r['wavelength']:.3e}       {r['lit']:4d}  "
              f"{r['power_a']:.3e} / {r['power_b']:.3e} W   {fits[0]:+.4f} / {fits[1]:+.4f} (-2 x chirp {-2*chirp:+.4f})")
    print()
    return rows


# --------------------------------------------------------------------------------------
# The sweep of the field's sampling along tau, field.n_wavelength, on four lines.

SWEEP_FACTORS = (1, 2, 5, 10, 20, 40)   # Slices a plane.


def sweep_measure(wdir, root, ref=None):
    """What one run of the sweep is judged on: the window's and the pulse's energy, the
    peak slice power at the exit, the exit power per beam slice against the reference
    run's over the reference's peak, the rms width of the exit spectrum relative to the
    carrier from the final field file, taken on the field's own planes, and the walk."""
    from lucifer.stats import read_stats
    from lucifer.field import read_field
    text = (wdir / f"{root}.log").read_text()
    f = footer(text)
    with read_stats(wdir / f"{root}.stats.h5") as st:
        power = np.asarray(st["field/total/power"])[-1]
    fld = read_field(wdir / f"{root}-final.wf.h5")
    u, dz, lam = fld["u"], fld["dz"], fld["wavelength"]
    a = u.sum(axis=(1, 2))
    spec = np.abs(np.fft.fft(a))**2
    kappa = 2 * np.pi * np.fft.fftfreq(a.size, d=dz)
    k0 = 2 * np.pi / lam
    tot = spec.sum()
    mean = (spec * kappa).sum() / tot if tot > 0 else float("nan")
    width = np.sqrt((spec * (kappa - mean)**2).sum() / tot) / k0 if tot > 0 else float("nan")
    out = dict(window_energy=f.get("window_energy", float("nan")), pulse_energy=f.get("pulse_energy", float("nan")),
               peak_power=float(power.max()), width=float(width), centroid=float(mean / k0),
               walk=f.get("walk_seconds", float("nan")), planes=int(u.shape[0]), power=power)
    if ref is not None:
        out["exit_diff"] = float(np.abs(power - ref["power"]).max() / np.abs(ref["power"]).max())
    return out


def sweep_machines(args, wdir):
    """The four lines, each as (name, directory, input file, window n_wavelength, extra
    overrides): the SASE tier's beam over 1280 slices from noise, and the sase, flash1 and
    lcls examples as they ship."""
    import shutil
    import tomllib
    here = pathlib.Path(__file__).resolve().parents[2]
    machines = []
    beam, dark = chirp_prepare(args, wdir)
    (wdir / "aramis_sweep.toml").write_text(to_toml_text(CHIRP_RUN_INPUT.format(
        lattice="chirp-line.bmad", beam=beam.name, field=dark.name)))
    machines.append(("aramis", wdir, "aramis_sweep.toml", 3, []))
    for name in ("sase", "flash1", "lcls"):
        d = wdir / name
        d.mkdir(exist_ok=True)
        src = here / "examples" / name
        for f in src.iterdir():
            if f.suffix in (".toml", ".bmad"):
                shutil.copyfile(f, d / f.name)
        shutil.copyfile(here / "examples" / "aramis.bmad", wdir / "aramis.bmad")
        with open(d / "lucifer.toml", "rb") as fh:
            n_wl = int(tomllib.load(fh)["window"]["n_wavelength"])
        machines.append((name, d, "lucifer.toml", n_wl, []))
    return machines


def to_toml_text(flat):
    from lucifer.flat import to_toml
    return to_toml(flat)


def sweep(args, wdir):
    print("--- sweep: the field's sampling along tau on four lines, from noise, under the spectral propagator")
    print(f"  {args.threads} threads; each run against the same line at the window's slice")
    print()
    results = {}
    for name, d, toml, n_wl, extra in sweep_machines(args, wdir):
        print(f"  {name}: window.n_wavelength = {n_wl}")
        print("    slices a plane   field.n_wavelength   planes   window energy   pulse energy   peak power   "
              "exit power vs finest   spectrum rms width   centroid   walk")
        ref = None
        for fac in SWEEP_FACTORS:
            root = f"sweep_{name}_{fac}"
            over = ["field.propagator=spectral", f"field.n_wavelength={fac * n_wl}", "run.keep_escaped_field=true"] + extra
            rc, text = run(args.exe, toml, over, root, d, args.threads)
            if rc != 0:
                m = re.search(r"\n\s*\[ERROR.*?\n(.*?)\n", text)
                reason = text.strip().splitlines()[-1] if text.strip() else "no message"
                for line in text.splitlines():
                    if "FIELD.N_WAVELENGTH" in line:
                        reason = line.strip()
                        break
                print(f"    {fac:6d}   {fac * n_wl:8d}   refused: {reason[:110]}")
                results[(name, fac)] = dict(refused=reason)
                continue
            r = sweep_measure(d, root, ref)
            if ref is None:
                ref = r
            results[(name, fac)] = r
            diff = r.get("exit_diff", 0.0)
            print(f"    {fac:6d}   {fac * n_wl:8d}   {r['planes']:6d}   {r['window_energy']:.4e} J   "
                  f"{r['pulse_energy']:.4e} J   {r['peak_power']:.4e} W   {diff:.3e}   {r['width']:.3e}   "
                  f"{r['centroid']:+.2e}   {r['walk']:.1f} s")
        print()
    return results


# --------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", required=True, help="the lucifer binary")
    p.add_argument("--work-dir", required=True, help="where to run; made if absent")
    p.add_argument("--validation-dir", help="a completed validation pass's work directory")
    p.add_argument("--threads", type=int, default=1, help="OMP_NUM_THREADS for every run")
    p.add_argument("--reuse", action="store_true", help="keep a run whose log and final field exist")
    p.add_argument("what", nargs="*", choices=["vacuum", "null", "bank", "cost", "seeded", "chirped", "sweep"],
                   help="which measurements; the first four by default")
    args = p.parse_args()
    global REUSE
    REUSE = args.reuse
    args.exe = pathlib.Path(args.exe).resolve()
    wdir = pathlib.Path(args.work_dir).resolve()
    wdir.mkdir(parents=True, exist_ok=True)
    what = args.what or ["vacuum", "null", "bank", "cost"]
    if any(w != "vacuum" for w in what) and not args.validation_dir:
        raise SystemExit("null, bank and cost need --validation-dir")
    for w in what:
        globals()[w](args, wdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
