#!/usr/bin/env python3
"""
Measure the spectral field propagator, field.propagator = "spectral", against the slice
propagator, as docs/spectral-field.md records. Nothing here runs in the keystone: the
prototype is outside it, and the page is hand written from this script's report.

    spectral_field.py --exe <lucifer> --work-dir <dir> [--validation-dir <dir>]
                      [--threads N] [vacuum] [null] [boundary] [cost]

The four measurements, each a subcommand and all four when none is named:

  vacuum    A chirped Gaussian pulse with no beam, drifted through a pipe under both
            propagators, against the exact paraxial solution with its spatio-temporal
            coupling: each frequency's Gaussian beam diffracts at its own Rayleigh
            range, and the solution is the superposition over the pulse's spectrum. The
            error of each propagator as a function of the pulse's relative bandwidth.
            Needs the tracker alone.
  null      The steady-state tier and the SASE tier under "spectral" with whole-slice
            rotation, compared against Genesis4 at the tier tolerances as
            compare_fel.py compares them, and the same under the continuous slippage,
            with the difference from the slice run placed against rho. Needs the work
            directory of a completed validation pass (--validation-dir), which holds
            the tier inputs, the imported dumps and the Genesis4 references.
  boundary  The SASE tier's window, 32 slices under a line that slips about 600, under
            each boundary form: the headroom alone, which refuses the run where the head
            fills; the headroom with the refusal disabled, to see what wraps; the cos^2
            sponge; and the hybrid of whole-slice rotation with zero padding. For each,
            the exit energies, the escaped energy, the beam's energy loss, and what
            the pulse file holds. Needs --validation-dir.
  cost      The field phase of the SASE tier's walk per step, from the footer's timing
            table, under both propagators, at the thread count given. Needs
            --validation-dir.

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
    m = re.search(r"head slice held at most\s+(\S+)\s+of the window energy, at s =\s+(\S+)\s+m", text)
    if m:
        out["head_max"], out["z_head_max"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"HEAD SLICE HOLDS\s+(\S+)\s+OF THE WINDOW'S ENERGY AT s =\s+(\S+)\s+m", text)
    if m:
        out["refused_share"], out["refused_at"] = float(m.group(1)), float(m.group(2))
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
# The null test, the boundary and the cost, on the validation directory's tiers.

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
        # The ramp on the SASE tier's window, 32 slices under a line that slips 1184, needs
        # a boundary at the head, and the padding is the one that discards what crosses as
        # the rotation's zero fill does, with its energy counted.
        for form, boundary in (("slices", "headroom"), ("ramp", "pad")):
            root = f"{tier}_spectral_{form}"
            rc, text = run(args.exe, f"{tier}.toml",
                           ["field.propagator=spectral", f"field.spectral_slippage={form}",
                            f"field.spectral_boundary={boundary}"],
                           root, wdir, args.threads)
            if rc != 0:
                print(text[-2000:])
                raise SystemExit(f"null: {root} failed")
            worst = compare_fel.compare_tier(f"{tier} spectral, slippage {form}",
                                             f"{wdir}/{root}.diag.txt", f"{vdir}/{gen}.out.h5",
                                             f"{wdir}/{root}-final.wf.h5", f"{vdir}/{gen}-final.fld.h5",
                                             f"{wdir}/{root}-final.beam.h5", f"{vdir}/{gen}-final.par.h5",
                                             tol, nslice)[0]
            d = differences(read_curves(f"{wdir}/{root}.stats.h5"), read_curves(f"{vdir}/{tier}.stats.h5"))
            own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
            results[(tier, form)] = (worst, own)
            print(f"  {tier} spectral {form}: against Genesis4 {worst:.6e} (tolerance {tol:.1e}), "
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
    for form, boundary in (("slices", "headroom"), ("ramp", "pad")):
        root = f"tdsase_filter_spectral_{form}"
        rc, text = run(args.exe, "tdsase.toml",
                       ["run.source_filter=true", "field.propagator=spectral",
                        f"field.spectral_slippage={form}", f"field.spectral_boundary={boundary}"],
                       root, wdir, args.threads)
        if rc != 0:
            print(text[-2000:])
            raise SystemExit(f"null: {root} failed")
        d = differences(read_curves(f"{wdir}/{root}.stats.h5"), ref)
        own = {k: float(np.abs(np.asarray(v.values)).max()) for k, v in d.items()}
        results[("tdsase_filter", form)] = own
        print(f"    spectral {form}: " + ", ".join(f"{k} {v:.2e}" for k, v in own.items()))
    print()
    return results


def boundary(args, wdir):
    from lucifer.stats import read_stats

    vdir = tier_inputs(args, wdir)
    print("--- boundary: the SASE tier's 32-slice window under each form at the head")
    forms = (
        ("slices", "the slice propagator", []),
        ("headroom", "spectral, ramp, headroom", ["field.propagator=spectral"]),
        ("headroom-open", "spectral, ramp, headroom, refusal off",
         ["field.propagator=spectral", "field.spectral_head_tolerance=2"]),
        ("sponge", "spectral, ramp, cos^2 sponge",
         ["field.propagator=spectral", "field.spectral_boundary=sponge"]),
        ("hybrid", "spectral, whole-slice rotation, zero padding",
         ["field.propagator=spectral", "field.spectral_slippage=slices", "field.spectral_boundary=pad"]),
    )
    rows = {}
    for key, label, over in forms:
        root = f"boundary_{key}"
        rc, text = run(args.exe, "tdsase.toml", over + ["run.keep_escaped_field=true"], root, wdir, args.threads)
        f = footer(text)
        row = dict(label=label, rc=rc, **f)
        if rc == 0:
            st = read_stats(wdir / f"{root}.stats.h5")
            row.update(closure(st))
            pulse = wdir / f"{root}-pulse.wf.h5"
            row["pulse_file"] = pulse.exists()
            if pulse.exists():
                import h5py
                with h5py.File(pulse) as h:
                    row["pulse_planes"] = h["data/1/meshes/electricField/x"].shape[0]
        rows[key] = row
        print(f"  {label}:")
        if rc != 0:
            print(f"    refused: head slice held {f.get('refused_share', float('nan')):.2e} at s = "
                  f"{f.get('refused_at', float('nan')):.3f} m")
            continue
        print(f"    exit window energy {row['window_energy']:.4e} J, escaped {row.get('escaped_energy', 0):.4e} J, "
              f"their sum {row['window_energy'] + row.get('escaped_energy', 0):.4e} J, pulse {row['pulse_energy']:.4e} J")
        print(f"    beam energy loss {row['beam_loss']:.4e} J")
        if "head_max" in row:
            print(f"    head slice held at most {row['head_max']:.2e} of the window energy, at s = {row['z_head_max']:.3f} m")
        print(f"    pulse file: {'%d planes' % row['pulse_planes'] if row['pulse_file'] else 'none written'}")
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
    # The ramp with the padding, the full cost: the kernel pass and the ramp pass. The
    # SASE tier's window needs the padding under the ramp (see null).
    for key, over in (("slices", []), ("spectral", ["field.propagator=spectral", "field.spectral_boundary=pad"])):
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

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", required=True, help="the lucifer binary")
    p.add_argument("--work-dir", required=True, help="where to run; made if absent")
    p.add_argument("--validation-dir", help="a completed validation pass's work directory")
    p.add_argument("--threads", type=int, default=1, help="OMP_NUM_THREADS for every run")
    p.add_argument("what", nargs="*", choices=["vacuum", "null", "boundary", "cost"],
                   help="which measurements; all four by default")
    args = p.parse_args()
    args.exe = pathlib.Path(args.exe).resolve()
    wdir = pathlib.Path(args.work_dir).resolve()
    wdir.mkdir(parents=True, exist_ok=True)
    what = args.what or ["vacuum", "null", "boundary", "cost"]
    if any(w != "vacuum" for w in what) and not args.validation_dir:
        raise SystemExit("null, boundary and cost need --validation-dir")
    for w in what:
        globals()[w](args, wdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
