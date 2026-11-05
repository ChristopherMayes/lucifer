#!/usr/bin/env python3
"""
Measure the spectral field propagator, field.propagator = "spectral", against the slice
propagator, and write the record docs/spectral-field.md is generated from. The keystone's
own check of the mode is check_spectral.py; nothing here runs in it.

    spectral_field.py --exe <lucifer> --work-dir <dir> [--validation-dir <dir>]
                      [--threads N] [--cost-threads 1,8] [--exe-open <lucifer>]
                      [--toolchain <text>] [--json <file>] [--reuse]
                      [vacuum] [null] [bank] [cost] [seeded] [chirped] [sweep] [all]

The seven measurements, each a subcommand and the first four when none is named. `all`
runs every one, the cost first on the quiet machine and then the sweep and the two FEL
cases, which take hours at 12 threads. With --json each measurement's result is written
into that file as it finishes, under its own name, and a file already there keeps the
measurements this run does not repeat. The committed record is
docs/generated/spectral/measured.json, which tests/scripts/report_spectral.py turns into the
page's tables and figures. With --reuse a run whose log and final field exist in the work
directory is not run again, so a set that takes hours is re-analyzed from its files.

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
            table, under both propagators, at each thread count of --cost-threads, and
            from the lowest the rate fel_cost_mod's fel_cost_tau$ holds, refit by the
            count that module makes. Needs --validation-dir.
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
  unavg     The unaveraged mode under the spectral propagator, on the sandwich probe of
            check_unaveraged.py over 200 slices of one wavelength from shot noise: the
            mixed line under both propagators and the all-averaged twin under both. For
            each run the ledger's closure where the mode runs, the exit power, the
            window's and the pulse's energy and the walk, and the exit fields of the
            pairs compared as the seeded case compares them. Needs the tracker alone.
  sweep     The field's sampling along tau, field.n_wavelength, at 1, 2, 5, 10, 20 and
            40 slices a plane on four lines from noise under the spectral propagator:
            the SASE tier's beam over 1280 slices, and the sase, flash1 and lcls
            examples. For each run the window's and the pulse's energy, the peak power,
            the exit power per slice against the run at one slice a plane, the exit
            spectrum's rms width and centroid, and the walk. A sampling above the cliff,
            fel_field_coop_cliff$ of the cooperation length, is refused by the tracker,
            and such a row runs on --exe-open when it is given, a tracker built with that
            constant raised, and is marked so in the record. Needs --validation-dir for
            the Aramis lattice.

The tracker is run from the work directory, into which the validation directory's
inputs are linked, so the validation directory is read and never written. Every number
in the record is written to six significant figures.
"""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import pathlib
import platform
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


Z0 = 376.730313668      # The impedance of free space [ohm], so a slice's power is sum |E|^2 dx^2 / 2 Z0.


def spectral(n_wavelength):
    """The overrides that select the spectral propagator at one plane a slice. Left unset,
    field.n_wavelength derives a plane coarser than the slice under this propagator, and
    a comparison against the slice run would measure that sampling and not the
    propagator. The sweep states the sampling itself."""
    return ["field.propagator=spectral", f"field.n_wavelength={n_wavelength}"]


def window_n_wavelength(toml):
    """window.n_wavelength of an input file."""
    import tomllib
    with open(toml, "rb") as fh:
        return int(tomllib.load(fh)["window"]["n_wavelength"])


def rounded(obj):
    """obj with every float written to six significant figures, arrays as lists, so the
    record is short and the same run writes the same bytes."""
    if isinstance(obj, (bool, int, str)) or obj is None:
        return obj
    if isinstance(obj, (float, np.floating)):
        x = float(obj)
        return x if not math.isfinite(x) else float(f"{x:.6g}")
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return rounded(obj.tolist())
    if isinstance(obj, dict):
        return {str(k): rounded(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [rounded(v) for v in obj]
    return obj


def field_curves(path):
    """The exit curves a figure draws from a final field file: the power in each plane
    [W], the power spectrum of the field summed over the grid against kappa / k0, with
    the carrier at zero, and the plane spacing [m]."""
    from lucifer.field import read_field
    fld = read_field(path)
    u, dz, lam = fld["u"], fld["dz"], fld["wavelength"]
    power = (np.abs(u)**2).sum(axis=(1, 2)) * fld["dx"]**2 / (2 * Z0)
    a = u.sum(axis=(1, 2))
    spec = np.fft.fftshift(np.abs(np.fft.fft(a))**2)
    kappa = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(a.size, d=dz)) / (2 * np.pi / lam)
    return dict(power=power, spectrum=spec, kappa_over_k0=kappa, dz=float(dz), wavelength=float(lam))


def provenance(args):
    """What ran: the tracker's version and Bmad's from --version, the build, the
    toolchain stated on the command line, the machine, the day and the thread counts."""
    out = subprocess.run([str(args.exe), "--version"], capture_output=True, text=True).stdout
    versions = dict(re.findall(r"(Lucifer|Bmad) version (\S+)", out))
    machine = platform.machine()
    if platform.system() == "Darwin":
        r = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            machine = r.stdout.strip()
    return dict(lucifer_version=versions.get("Lucifer", "unknown"), bmad_version=versions.get("Bmad", "unknown"),
                toolchain=args.toolchain, machine=machine,
                system=f"{platform.system()} {platform.release()}", day=datetime.date.today().isoformat(),
                threads=args.threads, cost_threads=list(args.cost_threads),
                exe_open=bool(args.exe_open))


def write_record(path, name, value, args):
    """Write value under name into the JSON record at path, keeping what is already
    there, and the provenance of this run."""
    record = {}
    if path.exists():
        record = json.loads(path.read_text())
    record["provenance"] = provenance(args)
    record[name] = rounded(value)
    path.write_text(json.dumps(record, indent=1) + "\n")


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
        for name, over in (("slices", []), ("spectral", spectral(1))):
            rc, text = run(args.exe, f"{tag}.toml", over, f"{tag}_{name}", wdir, args.threads)
            if rc != 0:
                print(text[-2000:])
                raise SystemExit(f"vacuum: {tag} {name} failed")
            fld = read_field(wdir / f"{tag}_{name}-final.wf.h5")
            u = np.transpose(fld["u"], (2, 1, 0))          # (nx, ny, nz)
            exact = vacuum_exact(a, tau, tau_c, dtau, r2, k0, w0, length, shift=dtau)
            errors[name] = np.abs(u - exact).max() / np.abs(exact).max()
        rows.append(dict(bandwidth=float(target), chirp=float(b), slices=float(errors["slices"]),
                         spectral=float(errors["spectral"])))
        print(f"  {target:9.1e}   {b:13.4e}    {errors['slices']:.3e}     {errors['spectral']:.3e}")
    print()
    return dict(wavelength=wavelength, n_slice=nz, sigma_slices=sigma / dtau, waist=w0, ngrid=ngrid, dx=dx,
                half_width=half, length=length, rayleigh=z_r, threads=args.threads, rows=rows)


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
    rows = []
    curves = {}
    for tier, gen, nslice, tol in tiers:
        recorded = compare_fel.compare_tier(f"{tier} slices (the validation pass's own run)",
                                            f"{vdir}/{tier}.diag.txt", f"{vdir}/{gen}.out.h5",
                                            f"{vdir}/{tier}-final.wf.h5", f"{vdir}/{gen}-final.fld.h5",
                                            f"{vdir}/{tier}-final.beam.h5", f"{vdir}/{gen}-final.par.h5",
                                            tol, nslice)[0]
        root = f"{tier}_spectral"
        rc, text = run(args.exe, f"{tier}.toml", spectral(window_n_wavelength(wdir / f"{tier}.toml")), root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"null: {root} failed")
        worst = compare_fel.compare_tier(f"{tier} spectral",
                                         f"{wdir}/{root}.diag.txt", f"{vdir}/{gen}.out.h5",
                                         f"{wdir}/{root}-final.wf.h5", f"{vdir}/{gen}-final.fld.h5",
                                         f"{wdir}/{root}-final.beam.h5", f"{vdir}/{gen}-final.par.h5",
                                         tol, nslice)[0]
        a, b = read_curves(f"{wdir}/{root}.stats.h5"), read_curves(f"{vdir}/{tier}.stats.h5")
        d = differences(a, b)
        own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
        rows.append(dict(tier=tier, recorded=float(recorded), spectral=float(worst), tolerance=float(tol),
                         n_slice=int(nslice), against_slices=own))
        curves[tier] = dict(s=a.s, power_slices=b.power.sum(axis=1), power_spectral=a.power.sum(axis=1),
                            exit_slices=b.power[-1], exit_spectral=a.power[-1])
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
    rc, text = run(args.exe, "tdsase.toml", ["run.source_filter=true"] + spectral(window_n_wavelength(wdir / "tdsase.toml")),
                   root, wdir, args.threads)
    if rc != 0:
        print(text[-2000:])
        raise SystemExit(f"null: {root} failed")
    a = read_curves(f"{wdir}/{root}.stats.h5")
    d = differences(a, ref)
    own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
    curves["tdsase_filter"] = dict(s=a.s, power_slices=ref.power.sum(axis=1), power_spectral=a.power.sum(axis=1),
                                   exit_slices=ref.power[-1], exit_spectral=a.power[-1])
    print(f"    spectral: " + ", ".join(f"{k} {v:.2e}" for k, v in own.items()))
    print()
    return dict(threads=args.threads, tiers=rows, filter_on=dict(against_slices=own), curves=curves)

def bank(args, wdir):
    from lucifer.stats import read_stats

    tier_inputs(args, wdir)
    print("--- bank: the SASE tier's 32-slice window with the escaped field kept, under both propagators")
    rows = {}
    n_wl = window_n_wavelength(wdir / "tdsase.toml")
    for key, label, over in (("slices", "the slice propagator", []),
                             ("spectral", "the spectral propagator", spectral(n_wl))):
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
        if pulse.exists():
            c = field_curves(pulse)
            row["pulse_power"] = c["power"]
        rows[key] = row
        print(f"  {label}:")
        print(f"    exit window energy {row['window_energy']:.4e} J, escaped {row.get('escaped_energy', 0):.4e} J, "
              f"their sum {row['window_energy'] + row.get('escaped_energy', 0):.4e} J, pulse {row['pulse_energy']:.4e} J")
        print(f"    beam energy loss {row['beam_loss']:.4e} J, pulse file {row['pulse_planes']} planes")
        print(f"    tail slice's exit power share {row['tail_share']:.2e}, head slice's {row['head_share']:.2e}")
    print()
    return dict(threads=args.threads, rows=rows)


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
    print(f"--- cost: the SASE tier's field phase per step at {', '.join(str(t) for t in args.cost_threads)} thread(s)")
    rows = []
    # The full cost under spectral: the kernel pass and the slip pass, the bank off as the
    # tier runs it.
    n_wl = window_n_wavelength(wdir / "tdsase.toml")
    for threads in args.cost_threads:
        for key, over in (("slices", []), ("spectral", spectral(n_wl))):
            root = f"cost_{key}_t{threads}"
            rc, text = run(args.exe, "tdsase.toml", over, root, wdir, threads)
            if rc != 0:
                print(text[-2000:])
                raise SystemExit(f"cost: {root} failed")
            f = footer(text)
            f.update(header_counts(text))
            f.update(threads=threads, propagator=key)
            rows.append(f)
            print(f"  {threads:2d} threads {key:9s} field solve {f['solve_seconds']:.3f} s over {f['n_step']} steps = "
                  f"{1e3 * f['solve_seconds'] / f['n_step']:.2f} ms a step, walk {f['walk_seconds']:.3f} s")
    fit = tau_fit(rows)
    print(f"  fel_cost_tau$ from the {fit['threads']}-thread rows: {fit['tau']:.3e} s per point per log2(n_t) per pass, "
          f"the pad's transverse share {fit['pad_seconds']:.2f} s of the excess {fit['excess_seconds']:.2f} s")
    print()
    return dict(rows=rows, fit=fit)


def header_counts(text):
    """The counts of a run the cost model takes, from its header: the slices, the grid,
    the field's planes, the pad's slots and the FEL steps."""
    out = {}
    m = re.search(r"^ Beam\s+(\d+) slices? x (\d+) particles", text, re.M) or \
        re.search(r"^ Beam\s+(\d+) slices?, (\d+) particles", text, re.M)
    if m:
        out["n_slice"], out["n_particle_slice"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"grid (\d+) points", text)
    if m:
        out["ngrid"] = int(m.group(1))
    m = re.search(r"^ Field\s+(\d+) planes of", text, re.M)
    out["n_field"] = int(m.group(1)) if m else out.get("n_slice")
    m = re.search(r"a pad of (\d+) free and (\d+) sponge", text)
    out["n_pad"] = int(m.group(1)) + int(m.group(2)) if m else 0
    return out


def cost_solve_rate():
    """fel_cost_solve$ as src/fel_cost_mod.f90 holds it, the transverse transform's rate
    per grid point per log2(ngrid) per step, read from the module so the fit and the
    estimate share one number."""
    src = pathlib.Path(__file__).resolve().parents[2] / "src" / "fel_cost_mod.f90"
    m = re.search(r"fel_cost_solve\$ = ([0-9.eE+-]+)_rp", src.read_text())
    if not m:
        raise SystemExit(f"{src} holds no fel_cost_solve$")
    return float(m.group(1))


def tau_fit(rows):
    """fel_cost_tau$ from the cost rows at their lowest thread count, by fel_cost_mod's
    count of what the spectral propagator adds to a step: the transverse transform pair
    on the pad's planes at fel_cost_solve$, and two passes along tau over every line of
    the window and pad, each a transform pair of n_t = planes + pad. The excess of the
    spectral row's field phase over the slice row's, less the pad's transverse share, over
    the points the two passes transform, per log2(n_t)."""
    threads = min(r["threads"] for r in rows)
    by = {r["propagator"]: r for r in rows if r["threads"] == threads}
    sp, sl = by["spectral"], by["slices"]
    ngrid, n_field, n_pad, nstep = sp["ngrid"], sp["n_field"], sp["n_pad"], sp["n_step"]
    n_t = n_field + n_pad
    excess = sp["solve_seconds"] - sl["solve_seconds"]
    pad = cost_solve_rate() * n_pad * ngrid**2 * math.log2(ngrid) * nstep
    points = 2 * ngrid**2 * n_t * math.log2(n_t) * nstep
    return dict(threads=threads, excess_seconds=excess, pad_seconds=pad, n_t=n_t, ngrid=ngrid, n_field=n_field,
                n_pad=n_pad, n_step=nstep, tau=(excess - pad) / points)

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

CHIRP_N_WAVELENGTH = 3           # The chirped cases' slice, PREP_INPUT's and CHIRP_RUN_INPUT's n_wavelength.
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
    curves = {}
    kappa = None
    for i, bw in enumerate(BANDWIDTHS):
        tag = f"seed{i}"
        seed, b, bw_eff = chirp_seed(wdir, dark, tag, bw)
        (wdir / f"{tag}.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad", beam=beam.name, field=seed.name)))
        (wdir / f"{tag}_quiet.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad", beam=quiet.name, field=seed.name)))
        # quiet: the segment with the beam's shot noise off, the seed the only light.
        for extent, over in (("quiet", ["run.track_end=UND##1"]), ("segment", ["run.track_end=UND##1"]), ("line", [])):
            toml = f"{tag}_quiet.toml" if extent == "quiet" else f"{tag}.toml"
            roots = {}
            for prop, more in (("slices", []), ("spectral", spectral(CHIRP_N_WAVELENGTH))):
                root = f"{tag}_{extent}_{prop}"
                rc, text = run(args.exe, toml, over + more, root, wdir, args.threads)
                if rc != 0:
                    print(text[-2000:])
                    raise SystemExit(f"seeded: {root} failed")
                roots[prop] = root
            r = chirp_compare(wdir, roots["slices"], roots["spectral"])
            rows.append(dict(bandwidth=float(bw_eff), chirp=float(b), extent=extent, **row_numbers(r)))
            cs, cp = field_curves(wdir / f"{roots['slices']}-final.wf.h5"), field_curves(wdir / f"{roots['spectral']}-final.wf.h5")
            curves[f"{tag}_{extent}"] = dict(bandwidth=float(bw_eff), extent=extent, power_slices=cs["power"],
                                             power_spectral=cp["power"], spectrum_slices=cs["spectrum"],
                                             spectrum_spectral=cp["spectrum"])
            kappa = cp["kappa_over_k0"]
            print(f"  {bw_eff:8.3f}   {extent:7s}  {r['offset']:6.3f}   {r['field_raw']:.3e} / {r['field']:.3e}    "
                  f"{r['wavelength']:.3e}       {r['lit']:4d}   {r['power_a']:.3e} / {r['power_b']:.3e} W")
    # Convergence of the reference: the 2 percent seed with the undulator step cut from
    # three periods to one, through the segment and through the line.
    tag = "seed2"
    convergence = []
    (wdir / f"{tag}_fine.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line-fine.bmad", beam=beam.name,
                                                                           field=f"{tag}-initial.wf.h5")))
    for extent, over in (("segment", ["run.track_end=UND##1"]), ("line", [])):
        rc, text = run(args.exe, f"{tag}_fine.toml", over + spectral(CHIRP_N_WAVELENGTH),
                       f"{tag}_{extent}_spectral_fine", wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"seeded: the fine-step {extent} run failed")
        r = chirp_compare(wdir, f"{tag}_{extent}_spectral_fine", f"{tag}_{extent}_spectral")
        print(f"  convergence, {extent}: the 2 percent seed, spectral at ds_step 0.015 against 0.045: "
              f"field {r['field']:.3e}, wavelength {r['wavelength']:.3e}, peak power {r['power_a']:.3e} / {r['power_b']:.3e} W")
        conv = dict(extent=extent, fine_against_spectral=row_numbers(r))
        r = chirp_compare(wdir, f"{tag}_{extent}_slices", f"{tag}_{extent}_spectral_fine")
        print(f"               and the slice run against that fine spectral run: field {r['field']:.3e} aligned by {r['offset']:.3f} slice")
        conv["slices_against_fine"] = row_numbers(r)
        convergence.append(conv)
    print()
    return dict(threads=args.threads, sigma_wavelengths=SEED_SIGMA_WAVELENGTHS, seed_slice=SEED_SLICE,
                waist=SEED_WAIST, peak_power=SEED_POWER, n_slice=int(kappa.size), kappa_over_k0=kappa,
                rows=rows, convergence=convergence, curves=curves)


def row_numbers(r):
    """The scalar entries of a chirp_compare result, for the record."""
    return {k: (float(v) if not isinstance(v, int) else v) for k, v in r.items() if not isinstance(v, np.ndarray)}


def chirped(args, wdir):
    from lucifer.flat import to_toml
    beam, dark = chirp_prepare(args, wdir)
    print("--- chirped: an energy-chirped beam from noise through Aramis under both propagators, source filter on")
    print(f"  1280 slices of 3 wavelengths, 1024 particles a slice, {args.threads} threads")
    print()
    print("  chirp    field diff   wavelength diff   lit   exit peak power slices / spectral   fitted slope of the wavelengths slices / spectral")
    rows = []
    curves = {}
    kappa = None
    for i, chirp in enumerate(ENERGY_CHIRPS):
        tag = f"chirp{i}"
        chirp_beam(beam, wdir / f"{tag}-initial.beam.h5", chirp)
        (wdir / f"{tag}.toml").write_text(to_toml(CHIRP_RUN_INPUT.format(lattice="chirp-line.bmad",
                                                                          beam=f"{tag}-initial.beam.h5", field=dark.name)))
        roots = {}
        for prop, more in (("slices", []), ("spectral", spectral(CHIRP_N_WAVELENGTH))):
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
        f = footer((wdir / f"{roots['spectral']}.log").read_text())
        rows.append(dict(chirp=float(chirp), slope_slices=float(fits[0]), slope_spectral=float(fits[1]),
                         expected_slope=-2 * chirp, pulse_energy_spectral=f.get("pulse_energy", float("nan")),
                         pulse_energy_slices=footer((wdir / f"{roots['slices']}.log").read_text()).get("pulse_energy", float("nan")),
                         **row_numbers(r)))
        cs, cp = field_curves(wdir / f"{roots['slices']}-final.wf.h5"), field_curves(wdir / f"{roots['spectral']}-final.wf.h5")
        curves[tag] = dict(chirp=float(chirp), power_slices=cs["power"], power_spectral=cp["power"],
                           spectrum_slices=cs["spectrum"], spectrum_spectral=cp["spectrum"],
                           wavelength_slices=r["lam_a"], wavelength_spectral=r["lam_b"])
        kappa = cp["kappa_over_k0"]
        print(f"  {chirp:6.4f}   {r['field']:.3e}    {r['wavelength']:.3e}       {r['lit']:4d}  "
              f"{r['power_a']:.3e} / {r['power_b']:.3e} W   {fits[0]:+.4f} / {fits[1]:+.4f} (-2 x chirp {-2*chirp:+.4f})")
    print()
    return dict(threads=args.threads, n_slice=int(kappa.size), kappa_over_k0=kappa, rows=rows, curves=curves)


# --------------------------------------------------------------------------------------
# The sweep of the field's sampling along tau, field.n_wavelength, on four lines.

SWEEP_FACTORS = (1, 2, 5, 10, 20, 40)   # Slices a plane.

# The unaveraged mode's window on the sandwich probe (tests/bmad/unavg_sandwich.bmad): the
# time-dependent input file of check_unaveraged.py's thread and ledger checks over 200
# slices of one wavelength, so the light the unaveraged segment emits stays in the window
# through the pipes and the trailing segment, the line slipping 139 wavelengths.
UNAVG_INPUT = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lattice}"
  out_root = "{root}"
  source_filter = T
  lambda0 = 1e-10
  beam.n_particle = 512
  beam.bunch_charge = 2.0e-13
  beam.z_distribution = "GRID"
  beam.z_min = -1e-8
  beam.z_max = 1e-8
  beam.sig_pz = 8.8045e-5
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 2e-8
  n_wavelength = 1
  shot_noise = T
  ran_seed = 999
  write_diag = T
&end
"""
UNAVG_SEGMENT = (0.81, 1.41)   # The unaveraged segment's span on the sandwich [m].


def unavg_closure(path):
    """The unaveraged ledger's closure over its rows: the largest excursion of
    E_beam + U_field + U_escaped - U_spont + E_radiated over the turnover."""
    led = np.loadtxt(path)
    etot = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turn = np.abs(np.diff(led[:, 2])).sum() + abs(led[-1, 4]) + abs(led[-1, 5]) + abs(led[-1, 6])
    return float(np.abs(etot - etot[0]).max() / max(turn, 1e-300))


def unavg(args, wdir):
    """The unaveraged mode under the spectral propagator on the sandwich probe: the
    mixed line under both propagators, and the all-averaged twin under both as the level
    the mode's own propagator difference is read against."""
    from lucifer.flat import to_toml
    lat = pathlib.Path(__file__).resolve().parent.parent / "bmad" / "unavg_sandwich.bmad"
    (wdir / "unavg_sandwich.bmad").write_bytes(lat.read_bytes())
    (wdir / "sandwich_avg.bmad").write_text("call, file = unavg_sandwich.bmad\nUNDB[FEL_METHOD] = averaged\n")
    print("--- unavg: the unaveraged mode under the spectral propagator, the sandwich probe over 200 slices")
    print(f"  512 particles a slice, shot noise, a 10 kW seed, {args.threads} threads")
    print()
    print("  run                      closure      exit power    window energy   pulse energy   walk")
    cases = [("unavg_slices", "unavg_sandwich.bmad", []), ("unavg_spectral", "unavg_sandwich.bmad", spectral(1)),
             ("avg_slices", "sandwich_avg.bmad", []), ("avg_spectral", "sandwich_avg.bmad", spectral(1))]
    rows, curves = [], {}
    kappa = None
    for root, lattice, over in cases:
        (wdir / f"{root}.toml").write_text(to_toml(UNAVG_INPUT.format(lattice=lattice, root=root)))
        rc, text = run(args.exe, f"{root}.toml", over, root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"unavg: {root} failed")
        f = footer(text)
        closure = unavg_closure(wdir / f"{root}.ledger.txt") if root.startswith("unavg") else None
        rows.append(dict(run=root, mode="unaveraged" if root.startswith("unavg") else "averaged",
                         propagator="slices" if root.endswith("slices") else "spectral",
                         closure=closure, **f))
        c = field_curves(wdir / f"{root}-final.wf.h5")
        curves[root] = dict(power=c["power"], spectrum=c["spectrum"])
        kappa = c["kappa_over_k0"]
        print(f"  {root:22s}  {closure if closure is not None else float('nan'):.2e}   {f['exit_power']:.4e} W  "
              f"{f['window_energy']:.4e} J  {f['pulse_energy']:.4e} J  {f['walk_seconds']:.1f} s")
    print()
    print("  comparison                                offset   field diff raw / aligned   wavelength diff   lit   peak power a / b")
    compare = []
    pairs = [("unavg_spectral", "unavg_slices"), ("avg_spectral", "avg_slices"), ("unavg_slices", "avg_slices")]
    for a, b in pairs:
        r = chirp_compare(wdir, a, b)
        compare.append(dict(a=a, b=b, **row_numbers(r)))
        print(f"  {a:>14s} against {b:<14s}  {r['offset']:6.3f}   {r['field_raw']:.3e} / {r['field']:.3e}    "
              f"{r['wavelength']:.3e}       {r['lit']:4d}   {r['power_a']:.3e} / {r['power_b']:.3e} W")
    print()
    # This measurement's own provenance, since it was taken after the record's others.
    return dict(threads=args.threads, n_slice=int(kappa.size), segment=list(UNAVG_SEGMENT), kappa_over_k0=kappa,
                measured=provenance(args), rows=rows, compare=compare, curves=curves)


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
               walk=f.get("walk_seconds", float("nan")), planes=int(u.shape[0]), power=power, wavelength=float(lam))
    if ref is not None:
        out["exit_diff"] = float(np.abs(power - ref["power"]).max() / np.abs(ref["power"]).max())
    return out


def sweep_machines(args, wdir):
    """The four lines, each as (name, directory, input file, window n_wavelength, extra
    overrides): the SASE tier's beam over 1280 slices from noise, and the sase, flash1 and
    lcls examples as they ship."""
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


def sweep_derived(args, d, toml, extra):
    """The line's rho, cooperation length [m] and derived sampling [wavelengths], from
    the header of a run that loads the input file and stops (run.load_only)."""
    rc, text = run(args.exe, toml, ["field.propagator=spectral", "run.load_only=true"] + extra, "sweep_derive", d, args.threads)
    if rc != 0:
        print(text[-2000:])
        raise SystemExit("sweep: the load-only run failed")
    m = re.search(r"Field sampling (\d+) wavelengths, (\d+) slices a plane, derived:.*?cooperation length\s+"
                  r"([0-9.Ee+-]+)\s+m at rho =\s+([0-9.Ee+-]+)", text, re.S)
    if not m:
        raise SystemExit("sweep: the header states no derived field sampling")
    return dict(derived_n_wavelength=int(m.group(1)), derived_slices_a_plane=int(m.group(2)),
                l_coop=float(m.group(3)), rho=float(m.group(4)))

def to_toml_text(flat):
    from lucifer.flat import to_toml
    return to_toml(flat)


def sweep(args, wdir):
    print("--- sweep: the field's sampling along tau on four lines, from noise, under the spectral propagator")
    print(f"  {args.threads} threads; each run against the same line at the window's slice")
    print()
    machines = {}
    for name, d, toml, n_wl, extra in sweep_machines(args, wdir):
        derived = sweep_derived(args, d, toml, extra)
        print(f"  {name}: window.n_wavelength = {n_wl}, rho {derived['rho']:.3e}, cooperation length {derived['l_coop']:.3e} m, "
              f"derived sampling {derived['derived_n_wavelength']} wavelengths")
        print("    slices a plane   field.n_wavelength   planes   window energy   pulse energy   peak power   "
              "exit power vs finest   spectrum rms width   centroid   walk")
        ref = None
        rows = []
        wavelength = None
        for fac in SWEEP_FACTORS:
            root = f"sweep_{name}_{fac}"
            over = ["field.propagator=spectral", f"field.n_wavelength={fac * n_wl}", "run.keep_escaped_field=true"] + extra
            exe, cliff_open = args.exe, False
            rc, text = run(exe, toml, over, root, d, args.threads)
            if rc != 0 and args.exe_open and "ABOVE" in text and "COOPERATION" in text:
                exe, cliff_open = args.exe_open, True
                rc, text = run(exe, toml, over, root, d, args.threads)
            row = dict(factor=fac, n_wavelength=fac * n_wl, cliff_open=cliff_open)
            if rc != 0:
                reason = text.strip().splitlines()[-1] if text.strip() else "no message"
                for line in text.splitlines():
                    if "FIELD.N_WAVELENGTH" in line:
                        reason = line.strip()
                        break
                print(f"    {fac:6d}   {fac * n_wl:8d}   refused: {reason[:110]}")
                row["refused"] = reason
                rows.append(row)
                continue
            r = sweep_measure(d, root, ref)
            if ref is None:
                ref = r
            wavelength = r["wavelength"]
            row.update({k: v for k, v in r.items() if k not in ("power", "wavelength")})
            row["fraction"] = fac * n_wl * wavelength / derived["l_coop"]
            rows.append(row)
            diff = r.get("exit_diff", 0.0)
            print(f"    {fac:6d}   {fac * n_wl:8d}   {r['planes']:6d}   {r['window_energy']:.4e} J   "
                  f"{r['pulse_energy']:.4e} J   {r['peak_power']:.4e} W   {diff:.3e}   {r['width']:.3e}   "
                  f"{r['centroid']:+.2e}   {r['walk']:.1f} s{'   (cliff open)' if cliff_open else ''}")
        for row in rows:
            if "fraction" not in row and wavelength:
                row["fraction"] = row["n_wavelength"] * wavelength / derived["l_coop"]
        machines[name] = dict(n_wavelength=n_wl, wavelength=wavelength, rows=rows, **derived)
        print()
    return dict(threads=args.threads, factors=list(SWEEP_FACTORS), machines=machines)

# --------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", required=True, help="the lucifer binary")
    p.add_argument("--exe-open", help="a lucifer built with fel_field_coop_cliff$ raised, for the sweep's rows above the cliff")
    p.add_argument("--work-dir", required=True, help="where to run; made if absent")
    p.add_argument("--validation-dir", help="a completed validation pass's work directory")
    p.add_argument("--threads", type=int, default=1, help="OMP_NUM_THREADS for every run but the cost's")
    p.add_argument("--cost-threads", default="1,8", help="the thread counts the cost rows run at")
    p.add_argument("--toolchain", default="unstated", help="the compiler the tracker was built with, for the record")
    p.add_argument("--json", help="the record to write each measurement into as it finishes")
    p.add_argument("--reuse", action="store_true", help="keep a run whose log and final field exist")
    p.add_argument("what", nargs="*", choices=["vacuum", "null", "bank", "cost", "seeded", "chirped", "sweep", "unavg", "all"],
                   help="which measurements; the first four by default")
    args = p.parse_args()
    global REUSE
    REUSE = args.reuse
    args.exe = pathlib.Path(args.exe).resolve()
    if args.exe_open:
        args.exe_open = pathlib.Path(args.exe_open).resolve()
    args.cost_threads = [int(t) for t in args.cost_threads.split(",")]
    wdir = pathlib.Path(args.work_dir).resolve()
    wdir.mkdir(parents=True, exist_ok=True)
    what = args.what or ["vacuum", "null", "bank", "cost"]
    if "all" in what:
        what = ["vacuum", "cost", "null", "bank", "sweep", "seeded", "chirped", "unavg"]
    if any(w not in ("vacuum", "unavg") for w in what) and not args.validation_dir:
        raise SystemExit("every measurement but vacuum needs --validation-dir")
    for w in what:
        result = globals()[w](args, wdir)
        if args.json:
            write_record(pathlib.Path(args.json), w, result, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
