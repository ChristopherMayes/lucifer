#!/usr/bin/env python3
"""
Check the spectral field propagator, field.propagator = "spectral", against the slice
propagator (docs/spectral-field.md).

The spectral propagator diffracts each mode of the window at its own wavenumber and
slips the field by a continuous phase ramp, where the slice propagator takes one
wavenumber and whole slices. On a steady-state window the two are the same operator, and
on a SASE run with the source filter on they differ by the wide-angle light the filter
leaves, which docs/spectral-field.md measures. Six checks:

  steady      the steady-state tier under spectral against the validation pass's own
              slice run: the power, the bunching, the energy and the sizes along s. The
              two propagators are one transform pair apart here, so the level is rounding.
  sase        the SASE tier with the source filter on under both propagators, Lucifer
              against Lucifer, the spectral run at one plane a slice: the exit power
              across the window. The tolerance is the measured level with its margin.
  one slice   a Gaussian pulse through a pipe whose slippage is one slice, with no beam,
              under both propagators: the ramp by one slice is the rotation by one slice
              to rounding, since a shift by a whole slot is exact in the transform.
  coarse      the SASE tier under spectral with field.n_wavelength left to derive, so the
              field is sampled coarser than the beam's slices as far as the tier's window
              of 32 slices allows, against the same tier at one plane a slice: the exit
              power per beam slice, which the stats take from the plane each slice reads.
              The bound is the sweep's (docs/spectral-field.md).
  cliff       the SASE tier set up alone, run.load_only, with the window's slice raised to
              192 wavelengths, which puts it at 0.31 of the cooperation length, above the
              cliff at a quarter: a stated field.n_wavelength equal to the slice runs and
              says where the slice sits, the key left unset derives the same value, and a
              stated value of twice the slice is refused with the cliff's message. The loaded
              beam keeps its 32 slices and stretches with the slice, so its current and its
              rho fall as the slice grows.
  bank        the same pulse placed so that it crosses the head over a pipe of 41 slices
              of slippage: the window's energy at the exit plus the escaped energy equals
              the energy that entered, the pulse file holds the 41 planes the header
              states, and the tail slice stays dark, so nothing the ramp carries past the
              head comes round to the tail. The ledger closes to rounding under the
              rotation and to 1e-6 under the ramp: the frequency-dependent kernel moves
              light along tau across the head on its own, by the group delay of the
              wide angles, and the bank samples the ramp's crossings alone.

The steady and SASE checks need the tier inputs and the imported dumps of a validation
pass in the work directory, which the validation suite copies in. The pulse checks need
the tracker alone.

The eleven tracker runs are independent, so they are launched together through
pool.run_all, longest first, and the checks are judged from their files afterwards in the
order above. The three SASE runs take nearly all of the section's time. Each run takes
the threads run_threads gives from the section's allowance, and the tracker writes the
same bits at any thread count, so only the wall time depends on it.

Usage:

  check_spectral.py --exe <lucifer> --workdir <dir> [--threads N] [--sase-tol <rel>]
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pool  # noqa: E402
from checkrec import record  # noqa: E402
from convert_genesis import write_openpmd_beam  # noqa: E402
from lucifer.compare import differences, read_curves  # noqa: E402
from lucifer.field import read_field  # noqa: E402
from lucifer.flat import to_toml  # noqa: E402

STEADY_TOL = 1e-11      # The steady state: one transform pair apart, 1e-12 measured.
COARSE_TOL = 1.5e-1     # The derived coarse field against one plane a slice on the SASE tier, 8.0e-2 measured
                        #   on a window of 32 slices at saturation (docs/spectral-field.md).
ONE_SLICE_TOL = 1e-8    # The ramp by a whole slice against the rotation, over the pulse's peak, 1e-9 measured.
BANK_TOL = {"slices": 1e-9, "spectral": 1e-5}   # The ledger's closure, relative to the energy that
                                                 #   entered: rounding under the rotation, and 1e-6 under
                                                 #   the ramp, the kernel's own transport along tau across
                                                 #   the head, which the crossings do not sample.
TAIL_TOL = 1e-9         # The tail slice's share of the exit energy.
Z0 = 376.730313668
C_LIGHT = 299792458.0

# The pulse's line (the vacuum test of spectral_field.py): one short undulator, so the
# lattice has an FEL element, and the pipe the pulse drifts through. The energy sets the
# pipe's slippage, floor(L / 2 gamma^2 lambda) + 1 wavelengths, one slice each at
# n_wavelength = 1.

PULSE_LATTICE = """! The spectral propagator's pulse line (check_spectral.py).
no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = {e_tot}

beginning[beta_a] = 10
beginning[beta_b] = 10

UND: wiggler, l = 0.03, l_period = 0.015, field_calc = planar_model, &
     b_max = 0.5 * sqrt(2) * (twopi / 0.015) * m_electron / c_light, &
     fel_method = averaged, ds_step = 0.015
DRIFT1: pipe, l = {length}

PULSE: line = (UND, DRIFT1)

use, PULSE
"""

PULSE_INPUT = """lattice = "{lattice}"
track_start = "DRIFT1"
track_end = "DRIFT1"
source_filter = F
keep_escaped_field = T
n_wavelength = 1
beam.file = "{beam}"
field.file = "{field}"
lambda0 = {wavelength}
"""

WAVELENGTH = 1e-9
NZ, NGRID = 1024, 128
SIGMA_SLICES = 50
WAIST = 20e-6
HALF = 200e-6


class _Patch:
    """One slice's particles in the form write_openpmd_beam reads."""

    def __init__(self, n, k, spacing, pz, weight):
        self.n_particle = n
        self.x = np.zeros(n)
        self.y = np.zeros(n)
        self.px = np.zeros(n)
        self.py = np.zeros(n)
        self.pz = np.full(n, pz)
        self.z = k * spacing + spacing * (np.arange(n) + 0.5) / n
        self.weight = np.full(n, weight)


# Every run's exit status and log text, by its root, filled by launch() and read by the
# checks through outcome() and result().
RUNS = {}


def launch(exe, wd, toml, root, overrides, threads):
    """Run the tracker on toml with the overrides into root, keeping its exit status and
    log text under root and writing the log to root.log."""
    cmd = [str(exe), toml, f"run.out_root={root}"] + list(overrides)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    p = subprocess.run(cmd, cwd=str(wd), env=env, capture_output=True, text=True)
    text = p.stdout + p.stderr
    (wd / f"{root}.log").write_text(text)
    RUNS[root] = (p.returncode, text)


def outcome(root):
    """A launched run's exit status and log text."""
    return RUNS[root]


def result(root):
    """A launched run's log text, or None after printing the log's tail when it failed."""
    rc, text = RUNS[root]
    if rc != 0:
        print(text[-3000:], file=sys.stderr)
        print(f"FAIL: {root} returned {rc}", file=sys.stderr)
        return None
    return text


def worst_of(d, keys=None):
    """The largest absolute difference over the keys of a differences() result."""
    out = {}
    for k, v in d.items():
        if keys is None or k in keys:
            out[k] = float(np.abs(np.asarray(v.values)).max())
    return out


def escaped_energy(text):
    """The escaped energy the footer states at full precision [J]."""
    m = re.search(r"Escaped\s+h1\s+(\S+)\s+J", text)
    return float(m.group(1)) if m else 0.0


def field_energy(f):
    """A field's energy [J] as the tracker sums it: |E|^2 over the grid times dx^2 / 2 Z0 is
    a slice's power, times the slice's light time its energy."""
    return float((np.abs(f["u"])**2).sum() * f["dx"]**2 / (2 * Z0) * f["dz"] / C_LIGHT)


def write_pulse(wd, tag, center_slice, e_tot, length):
    """The Gaussian pulse at center_slice, its beam of no consequence, the line and the
    input file, as the vacuum test writes them."""
    from beamphysics.wavefront import Wavefront

    dtau = WAVELENGTH
    dx = 2 * HALF / (NGRID - 1)
    x = (np.arange(NGRID) - (NGRID - 1) / 2) * dx
    r2 = x[:, None]**2 + x[None, :]**2
    tau = np.arange(NZ) * dtau
    sigma = SIGMA_SLICES * dtau
    env = np.exp(-(tau - center_slice * dtau)**2 / (4 * sigma**2))
    ex = np.exp(-r2 / WAIST**2)[:, :, None] * env[None, None, :]
    Wavefront(Ex=ex.astype(np.complex128), dx=dx, dy=dx, dz=dtau,
              wavelength=WAVELENGTH).write_openpmd(str(wd / f"{tag}-initial.wf.h5"))
    write_openpmd_beam(wd / f"{tag}-initial.beam.h5", [_Patch(4, k, dtau, e_tot, 1e-21) for k in range(NZ)], dtau)
    (wd / f"{tag}.bmad").write_text(PULSE_LATTICE.format(e_tot=repr(e_tot), length=repr(length)))
    (wd / f"{tag}.toml").write_text(to_toml(PULSE_INPUT.format(
        lattice=f"{tag}.bmad", beam=f"{tag}-initial.beam.h5", field=f"{tag}-initial.wf.h5",
        wavelength=repr(WAVELENGTH))))
    return float((np.abs(ex)**2).sum())


def steady(wd):
    text = result("spec_tier1")
    if text is None:
        return False
    d = worst_of(differences(read_curves(wd / "spec_tier1.stats.h5"), read_curves(wd / "tier1.stats.h5")))
    worst = max(d.values())
    good = worst <= STEADY_TOL
    print(f"  steady state, spectral against the slice run: " +
          ", ".join(f"{k} {v:.2e}" for k, v in d.items()) + f", worst {worst:.2e} (tolerance {STEADY_TOL:.0e})")
    record("steady state, spectral against the slice run", worst, STEADY_TOL, "max", good)
    if not good:
        print(f"FAIL: the steady state differs by {worst:.2e}", file=sys.stderr)
    return good


def window_n_wavelength(wd, toml):
    """window.n_wavelength of an input file, so a spectral run can be held to one plane a slice."""
    import tomllib
    with open(wd / toml, "rb") as fh:
        return int(tomllib.load(fh)["window"]["n_wavelength"])


def sase(wd, tol):
    roots = {}
    for key in ("slices", "spectral"):
        root = f"spec_tdsase_{key}"
        if result(root) is None:
            return False
        roots[key] = root
    d = worst_of(differences(read_curves(wd / f"{roots['spectral']}.stats.h5"),
                             read_curves(wd / f"{roots['slices']}.stats.h5")), ("exit_power",))
    worst = max(d.values())
    good = worst <= tol
    print(f"  SASE tier, source filter on, spectral against slices: " +
          ", ".join(f"{k} {v:.2e}" for k, v in d.items()) + f", worst {worst:.2e} (tolerance {tol:.0e})")
    record("SASE tier with the source filter on, spectral against slices", worst, tol, "max", good)
    if not good:
        print(f"FAIL: the SASE tier differs by {worst:.2e} between the propagators", file=sys.stderr)
    return good


def coarse(wd):
    root = "spec_tdsase_coarse"
    text = result(root)
    if text is None:
        return False
    if outcome("spec_tdsase_spectral")[0] != 0:
        print("FAIL: the coarse field is compared with the spectral SASE run, which failed", file=sys.stderr)
        return False
    m = re.search(r"Field\s+(\d+) planes of (\d+) wavelengths over the window, (\d+) slices a plane", text)
    planes, n_wl, n_per = (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else (0, 0, 1)
    d = worst_of(differences(read_curves(wd / f"{root}.stats.h5"),
                             read_curves(wd / "spec_tdsase_spectral.stats.h5")), ("exit_power",))
    worst = max(d.values())
    good = worst <= COARSE_TOL and n_per > 1
    print(f"  coarse field, {planes} planes of {n_wl} wavelengths, {n_per} slices a plane, against one a slice: "
          f"exit power {worst:.2e} (tolerance {COARSE_TOL:.1e})")
    record("the derived coarse field against one plane a slice, SASE tier", worst, COARSE_TOL, "max", good)
    record("the derived coarse field's slices a plane, SASE tier", n_per, 1, "min", n_per > 1)
    if not good:
        print(f"FAIL: the coarse field differs by {worst:.2e}, or was not coarser ({n_per} slices a plane)",
              file=sys.stderr)
    return good


CLIFF_N_WAVELENGTH = 192   # The tier's slice raised above the cliff, 0.31 of the cooperation length.
COOP_CLIFF = 0.25          # fel_field_coop_cliff$ (fel_spectral_mod).


def cliff():
    # The setup alone of the SASE tier with its slice raised above the cliff, the key
    # stated at the slice, left unset, and stated at twice the slice.
    out = {key: outcome(f"spec_cliff_{key}") for key in ("stated", "unset", "coarser")}

    rc, text = out["stated"]
    m = re.search(r"Field sampling the window's slice, (\d+) wavelengths, stated: the slice is\s+([0-9.]+) of", text)
    frac = float(m.group(2)) if m else 0.0
    good_s = rc == 0 and m is not None and int(m.group(1)) == CLIFF_N_WAVELENGTH and frac > COOP_CLIFF
    print(f"  cliff, the slice stated: exit {rc}, the slice {frac:.2f} of the cooperation length "
          f"(above {COOP_CLIFF})")
    record("a stated sampling equal to a slice above the cliff runs, the slice over the cooperation length",
           frac, COOP_CLIFF, "min", good_s)

    rc, text = out["unset"]
    m = re.search(r"Field sampling (\d+) wavelengths, (\d+) slices a plane, derived", text)
    n_derived = int(m.group(1)) if m else 0
    good_u = rc == 0 and n_derived == CLIFF_N_WAVELENGTH
    print(f"  cliff, the key unset: exit {rc}, derived {n_derived} wavelengths")
    record("the key unset derives the same slice above the cliff", n_derived, CLIFF_N_WAVELENGTH, "equal", good_u)

    rc, text = out["coarser"]
    good_c = rc != 0 and "IS ABOVE" in text and "OF THE COOPERATION" in text
    print(f"  cliff, twice the slice stated: exit {rc}, refused with the cliff's message {good_c}")
    record("a stated sampling coarser than a slice above the cliff is refused", None, None, "bool", good_c)

    if not (good_s and good_u and good_c):
        print("FAIL: the stated sampling at a slice above the cliff", file=sys.stderr)
    return good_s and good_u and good_c


def one_slice(wd):
    fields = {}
    for key in ("slices", "spectral"):
        if result(f"spec_one_{key}") is None:
            return False
        fields[key] = read_field(wd / f"spec_one_{key}-final.wf.h5")["u"]
    diff = float(np.abs(fields["spectral"] - fields["slices"]).max() / np.abs(fields["slices"]).max())
    good = diff <= ONE_SLICE_TOL
    print(f"  one slice by the ramp against one by the rotation: {diff:.2e} of the peak (tolerance {ONE_SLICE_TOL:.0e})")
    record("one slice by the ramp against one by the rotation", diff, ONE_SLICE_TOL, "max", good)
    if not good:
        print(f"FAIL: the ramp by one slice differs from the rotation by {diff:.2e}", file=sys.stderr)
    return good


# The bank's pulse: 1.14 GeV over 0.4 m, 41 wavelengths of slippage, and the pulse's
# center 100 slices from the head, so most of it crosses.
BANK_E_TOT = 1.14e9
BANK_LENGTH = 0.4


def bank(wd):
    gamma = BANK_E_TOT / 0.51099895e6
    n_slip = int(np.floor(BANK_LENGTH / (2 * gamma**2 * WAVELENGTH))) + 1
    ok = True
    for key in ("slices", "spectral"):
        root = f"spec_bank_{key}"
        text = result(root)
        if text is None:
            return False
        escaped = escaped_energy(text)
        f0 = read_field(wd / "spec_bank-initial.wf.h5")
        f1 = read_field(wd / f"{root}-final.wf.h5")
        e0, window = field_energy(f0), field_energy(f1)
        # The ledger: the window's energy at the exit plus the escaped energy, against
        # the energy that entered, a drift being unitary.
        closure = abs((window + escaped) / e0 - 1)
        tail = float((np.abs(f1["u"][0])**2).sum() / max((np.abs(f1["u"])**2).sum(), 1e-300))
        with h5py.File(wd / f"{root}-pulse.wf.h5") as h:
            planes = int(h["data/1/meshes/electricField/x"].shape[0]) - NZ
        good_c = closure <= BANK_TOL[key]
        good_t = tail <= TAIL_TOL
        good_p = planes == n_slip
        print(f"  bank, {key}: window {window:.6e} J, escaped {escaped:.6e} J, closure {closure:.2e} "
              f"(tolerance {BANK_TOL[key]:.0e}), tail slice {tail:.2e} of the exit energy, pulse file {planes} planes "
              f"of {n_slip} slipped")
        record(f"the ledger closes, {key}", closure, BANK_TOL[key], "max", good_c)
        record(f"the tail slice stays dark, {key}", tail, TAIL_TOL, "max", good_t)
        record(f"the pulse file holds the planes slipped, {key}", planes, n_slip, "equal", good_p)
        if not (good_c and good_t and good_p):
            print(f"FAIL: bank under {key}: closure {closure:.2e}, tail {tail:.2e}, planes {planes} of {n_slip}",
                  file=sys.stderr)
            ok = False
    return ok


def run_threads(budget):
    """The threads each run takes and the threads pool.run_all divides the allowance by,
    from the section's allowance."""
    rule = os.environ.get("LUCIFER_SPECTRAL_THREADS") or "each"
    if rule == "serial":
        return budget, budget
    if rule == "each":
        return budget, max(1, budget // 3)
    t = max(1, budget // 3)
    return t, t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--threads", type=int, default=int(os.environ.get("LUCIFER_CPU_BUDGET", "4") or 4))
    p.add_argument("--sase-tol", type=float, default=2e-2,
                   help="the bound on the SASE tier's exit power between the propagators, 8.7e-3 measured")
    a = p.parse_args()
    wd = pathlib.Path(a.workdir).resolve()
    exe = pathlib.Path(a.exe).resolve()

    for f in ("tier1.toml", "tier1.stats.h5", "tdsase.toml"):
        if not (wd / f).exists():
            print(f"FAIL: the work directory lacks {f} from a validation pass", file=sys.stderr)
            return 1

    # The pulses' inputs. 100 GeV over a micron slips the floor's one wavelength, one slice.
    write_pulse(wd, "spec_one", NZ // 2, 100e9, 1e-6)
    write_pulse(wd, "spec_bank", NZ - 100, BANK_E_TOT, BANK_LENGTH)

    # The runs, longest first. The spectral SASE run is at one plane a slice, since
    # field.n_wavelength would otherwise derive a coarser field, which the coarse run
    # measures on its own.
    n_wl = window_n_wavelength(wd, "tdsase.toml")
    sase_on = ["run.source_filter=true"]
    cliff_base = ["run.load_only=true", "field.propagator=spectral", f"window.n_wavelength={CLIFF_N_WAVELENGTH}"]
    runs = [
        ("tdsase.toml", "spec_tdsase_spectral", sase_on + ["field.propagator=spectral", f"field.n_wavelength={n_wl}"]),
        ("tdsase.toml", "spec_tdsase_slices", sase_on),
        ("tdsase.toml", "spec_tdsase_coarse", sase_on + ["field.propagator=spectral"]),
        ("spec_bank.toml", "spec_bank_spectral", ["field.propagator=spectral"]),
        ("tier1.toml", "spec_tier1", ["field.propagator=spectral"]),
        ("spec_one.toml", "spec_one_spectral", ["field.propagator=spectral"]),
        ("spec_one.toml", "spec_one_slices", []),
        ("spec_bank.toml", "spec_bank_slices", []),
        ("tdsase.toml", "spec_cliff_stated", cliff_base + [f"field.n_wavelength={CLIFF_N_WAVELENGTH}"]),
        ("tdsase.toml", "spec_cliff_unset", cliff_base),
        ("tdsase.toml", "spec_cliff_coarser", cliff_base + [f"field.n_wavelength={2 * CLIFF_N_WAVELENGTH}"]),
    ]
    threads, per_job = run_threads(a.threads)
    pool.run_all([lambda r=r: launch(exe, wd, r[0], r[1], r[2], threads) for r in runs],
                 threads_per_job=per_job)

    ok = steady(wd)
    ok = sase(wd, a.sase_tol) and ok
    ok = coarse(wd) and ok
    ok = cliff() and ok
    ok = one_slice(wd) and ok
    ok = bank(wd) and ok
    print("  spectral checks: PASS" if ok else "  spectral checks: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
