#!/usr/bin/env python3
"""
Harmonic field-set and openPMD wavefront checks (fel-physics.md sec-field-set):

  1. openPMD round trip: a run's field dump must carry every required attribute of
     the standard (verified against the spec text's names and values here,
     independently of the Fortran writer), its complex values must agree with the
     converter's Genesis view of the same field, and the Fortran reader must
     reproduce it exactly: re-importing the .wf.h5 and re-dumping must be
     dataset-identical to the file that was read. The Python Wavefront class reads
     the converted Genesis file through its own path and must agree on the field
     energy, reads the harmonic .wf.h5 to the same complex values, round-trips its
     own writer exactly, and what its writer produces must come back through the
     Fortran reader unchanged.
  2. Harmonic tier vs Genesis4: a planar steady-state segment (the benchmark's
     gamma, wavelength and rms aw -- planar so fc(3) is alive), both codes tracking
     the same Genesis-written starting state with a dark third-harmonic field:
     fundamental and third-harmonic power curves compared per record.
  3. Deposit closed form: from the tier's bunched exit beam, a two-step dark restart
     radiates P_h proportional to (fc(h) |b_h|)^2 -- the ratio P3/P1 against the
     Bessel closed form with b_h measured directly from the dumped particles. No
     gain, no diffraction to speak of: this pins the harmonic deposit normalization
     independently of Genesis.
  4. Thread identity: 1 vs 8 threads byte-identical on a time-dependent harmonic run
     (diag byte-equal, harmonic dumps dataset-equal).
  5. Refusals: harmonics not anchored on the fundamental; harmonic
     fields with an unaveraged element; harmonic fields with two live polarizations;
     an openPMD import declaring the frequency domain; a harmonic import whose
     photonEnergy matches no field of the run; a Genesis-format harmonic import,
     which carries no photon energy to match on.

Run by the benchmark harness. Exits nonzero on failure. Needs the genesis4 binary
(--genesis) for section 2, and a beamphysics checkout carrying
beamphysics/wavefront/openpmd.py (--pyrepo) for the openPMD checks of section 1.
"""

from __future__ import annotations

import os
import argparse
import pathlib
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

import beamio
import convert_genesis
import fieldio
from read_stats import read_stats
from nml import to_groups
from scipy.special import jv

FAILED = False

MU0_C = 1.25663706127e-6 * 2.99792458e8      # Bmad's mu_0_vac (2018 CODATA), not 4pi e-7.
H_PLANCK_EVS = 4.135667696e-15
E_CHARGE = 1.602176634e-19
C_LIGHT = 2.99792458e8

# Tolerances: measured first (values in comments), then set.
TOL_RT_FORMATS = 1e-13       # measured 4e-16 with Bmad's own mu_0: dfl double rounding.
TOL_PY_ENERGY = 1e-9         # measured 1.2e-12: Python class energy via its Genesis path.
TOL_TIER_P1 = 1e-6           # measured 5.3e-8: fundamental power vs Genesis (planar SS).
TOL_TIER_P3 = 1e-3           # measured 1.3e-4: third-harmonic power vs Genesis (dark growth).
TOL_DEPOSIT = 1e-12          # measured 3.3e-16: one-step P3/P1 vs the exact deposit sum.

GENESIS_DECK = """&setup
rootname=H3
lattice=planar.lat
beamline=SEGP
lambda0=1e-10
gamma0=11357.82
delz=0.045000
shotnoise=0
nbins = 8
beam_global_stat = true
field_global_stat = true
&end

&field
power=5e3
dgrid=2.000000e-04
ngrid=255
waist_size=30e-6
&end

&field
harm=3
power=0
dgrid=2.000000e-04
ngrid=255
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
betax=15
betay=15
&end

&write
field = H3-initial
beam = H3-initial
&end

&track
fft_fieldsolver = true
&end

&write
field = H3-final
beam = H3-final
&end
"""

GENESIS_LAT = "UNDP: UNDULATOR = { lambdau=0.015000, nwig=264, aw=0.84853};\nSEGP: LINE={UNDP};\n"

BMAD_LAT = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron

beginning[beta_a] = 15
beginning[beta_b] = 15

UNDP: wiggler, l = {length}, l_period = 0.015, field_calc = planar_model, &
      b_max = sqrt(2) * 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = fel_averaged, ds_step = 0.045

SEGP: line = (UNDP)

use, SEGP
"""

# lambda0 and nbins are the deck's now: an openPMD beam file carries the slice partition
# and not the radiation it was sliced on, nor the beamlet size. Both match the Genesis deck
# above, which is where the imported state comes from.
NML_IMPORT = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  transport_model = "genesis"
  beam_file = "{beam}"
  field_file = "{field}"
  lambda0 = 1e-10
  beamlet_size = 8
  harmonics = 1, 3
  write_diag = T
{extra}&end
"""

NML_TD = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  transport_model = "genesis"
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 63
  grid_half_width = 2e-4
  window_length = 8e-10
  window_sample = 1
  shot_noise = T
  ran_seed = 777
  harmonics = 1, 3
  write_diag = T
{extra}&end
"""


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
        sys.exit(1)


def refuse(exe, wd, name, text, fragment):
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    return r.returncode != 0 and fragment in (r.stdout + r.stderr)


def fc_planar(aw, h):
    """fel_und_coupling's planar closed form: aw * JJ, JJ = J_{(h-1)/2} - J_{(h+1)/2}
    at xi = (h/2) aw^2/(1+aw^2), sign (-1)^((h-1)/2)."""
    xi = 0.5 * aw * aw / (1 + aw * aw) * h
    h0 = (h - 1) // 2
    return aw * (jv(h0, xi) - jv(h0 + 1, xi)) * (-1) ** h0


def fields_identical(fa, fb):
    """Every field component of two openPMD wavefront files, value for value."""
    ua = fieldio.read_field(fa)
    ub = fieldio.read_field(fb)
    if ua["components"] != ub["components"]:
        return False
    return all(np.array_equal(fieldio.read_field(fa, c)["u"], fieldio.read_field(fb, c)["u"])
               for c in ua["components"])


def bessel_p3_over_p1(beamfile, darkfile, aw=0.84853, lambda_u=0.015, lambda0=1e-10):
    """The one-step dark deposit's P3/P1 predicted from the dumped particles alone: the
    bilinear scatter of (sin h*theta + i cos h*theta) sqrt(faw2)/gamma over the dark
    field's grid, weighted by |fc(h)|^2 from the Bessel closed form. Common factors
    (weights, delz, spacing, mu0 c) cancel in the ratio. Shared with check_device.py,
    which holds the device's deposit to the same identity."""
    sl = beamio.read_slices(beamfile, lambda0, lambda0)[0]
    theta, gam, xp, yp = sl["theta"], sl["gamma"], sl["x"], sl["y"]
    dark = fieldio.read_field(darkfile)
    ng, dgrid = dark["u"].shape[1], dark["dx"]
    ku = 2 * np.pi / lambda_u
    ky = ku * ku                      # Planar natural-focusing split: kx = 0, ky = ku^2.
    gmax = (ng - 1) * dgrid / 2

    def predicted_power(h):
        crs = np.zeros((ng, ng), dtype=complex)
        part = np.sqrt(1 + ky * yp * yp) / gam
        cpart = (np.sin(h * theta) + 1j * np.cos(h * theta)) * part
        wx = (xp + gmax) / dgrid
        wy = (yp + gmax) / dgrid
        ix = np.floor(wx).astype(int); iy = np.floor(wy).astype(int)
        fx = 1 + ix - wx; fy = 1 + iy - wy
        on = (np.abs(xp) < gmax) & (np.abs(yp) < gmax)
        for dx_i, dy_i, wgt in ((0, 0, fx * fy), (1, 0, (1 - fx) * fy),
                                (0, 1, fx * (1 - fy)), (1, 1, (1 - fx) * (1 - fy))):
            np.add.at(crs, (ix[on] + dx_i, iy[on] + dy_i), (wgt * cpart)[on])
        return abs(fc_planar(aw, h)) ** 2 * float((abs(crs) ** 2).sum())

    return predicted_power(3) / predicted_power(1)


def dfl_to_vperm(fname, dx_expected=None):
    """A Genesis field dump's slice 1 as complex V/m, plus its grid spacing."""
    with h5py.File(fname) as h5:
        ng = int(h5["gridpoints"][0])
        ds = float(h5["gridsize"][0])
        fr = h5["slice000001/field-real"][:].reshape(ng, ng)
        fi = h5["slice000001/field-imag"][:].reshape(ng, ng)
    return (fr + 1j * fi) / (ds / np.sqrt(2 * MU0_C)), ds


def power_of(E, dx):
    return float((abs(E) ** 2).sum()) * dx * dx / (2 * MU0_C)


def main():
    global FAILED
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--genesis", required=True)
    ap.add_argument("--pyrepo", help="openPMD-beamphysics checkout. "
                    "Default: $OPENPMD_BEAMPHYSICS, or the installed package.")
    args = ap.parse_args()
    wd = pathlib.Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    exe = pathlib.Path(args.exe).resolve()

    (wd / "planar.lat").write_text(GENESIS_LAT)
    (wd / "planar.bmad").write_text(BMAD_LAT.format(length="3.96"))
    (wd / "short.bmad").write_text(BMAD_LAT.format(length="0.045"))
    (wd / "mid.bmad").write_text(BMAD_LAT.format(length="1.98"))

    # ------------------------------------------------------------------
    # 2 first (its dumps feed 1 and 3): the harmonic tier vs Genesis.

    print("== harmonic tier vs Genesis4 (planar SS, harm 1 + dark harm 3) ==")
    (wd / "H3.in").write_text(GENESIS_DECK)
    r = subprocess.run([args.genesis, "H3.in"], cwd=wd, capture_output=True, text=True,
                       env=dict(os.environ, FI_PROVIDER="tcp"))
    if r.returncode != 0:
        print(f"FAIL: genesis exited {r.returncode}:\n{r.stdout[-2000:]}\n{r.stderr[-500:]}")
        sys.exit(1)

    # The tracker reads openPMD only, so Genesis's starting state converts once here and
    # both codes still track the same beam and the same field.
    for src, dst in (("H3-initial.par.h5", "H3-initial.beam.h5"),
                     ("H3-initial.fld.h5", "H3-initial.wf.h5")):
        convert_genesis.to_openpmd(wd / src, wd / dst, args.pyrepo)

    run(exe, wd, "h3bmad", NML_IMPORT.format(lat="planar.bmad", root="h3bmad",
        beam="H3-initial.beam.h5", field="H3-initial.wf.h5", extra=""))

    with h5py.File(wd / "H3.out.h5") as h5:
        zg = h5["Lattice/zplot"][:]
        p1g = h5["Field/power"][:].ravel()
        p3g = h5["Field3/power"][:].ravel()
    with read_stats(wd / "h3bmad.stats.h5") as st:
        zb = st.s
        p1b = st["field/total/power"][:, 0]
        p3b = st["field/harm3/total/power"][:, 0]

    n = min(len(zg), len(zb))
    if not np.allclose(zg[:n], zb[:n], atol=1e-9):
        print("FAIL: z grids disagree between the codes"); sys.exit(1)
    # Relative on the curve maxima, skipping the seed-dominated first steps for P3
    # (a dark field's first records sit at the numerical floor).
    m1 = float(np.max(np.abs(p1b[:n] - p1g[:n])) / np.max(p1g[:n]))
    lo = n // 8
    m3 = float(np.max(np.abs(p3b[lo:n] - p3g[lo:n])) / np.max(p3g[:n]))
    check("tier P1 vs Genesis (planar SS)", m1, TOL_TIER_P1)
    check("tier P3 vs Genesis (dark start, growth)", m3, TOL_TIER_P3,
          note=f"[P3 exit {p3b[n-1]:.3e} W]")

    # ------------------------------------------------------------------
    # 1. openPMD round trip, on the tier's Bmad run re-dumped in both formats.

    print("== openPMD EXT_Wavefront round trip ==")
    run(exe, wd, "h3both", NML_IMPORT.format(lat="planar.bmad", root="h3both",
        beam="H3-initial.beam.h5", field="H3-initial.wf.h5", extra=""))

    # The Genesis view of each dump, through the converter: the second reading of the same
    # field that the complex comparison below needs.
    for fn in ("h3both-final.wf.h5", "h3both-final-h3.wf.h5"):
        convert_genesis.to_openpmd(wd / fn, wd / fn.replace(".wf.h5", ".fld.h5"), args.pyrepo)

    for fn, lam in (("h3both-final.wf.h5", 1e-10), ("h3both-final-h3.wf.h5", 1e-10 / 3)):
        with h5py.File(wd / fn) as f:
            for a, want in (("openPMD", b"2.0.0"), ("openPMDextension", b"Wavefront"),
                            ("basePath", b"/data/%T/"), ("meshesPath", b"meshes/"),
                            ("iterationEncoding", b"groupBased"),
                            ("iterationFormat", b"/data/%T/")):
                if f.attrs[a] != want:
                    print(f"FAIL: {fn} root attribute {a} = {f.attrs[a]!r}, want {want!r}")
                    FAILED = True
            m = f["data/1/meshes/electricField"]
            need = {"geometry": b"cartesian", "temporalDomain": b"time", "spatialDomain": b"r"}
            for a, want in need.items():
                if m.attrs[a] != want:
                    print(f"FAIL: {fn} mesh attribute {a} = {m.attrs[a]!r}, want {want!r}")
                    FAILED = True
            for a in ("axisLabels", "gridSpacing", "gridGlobalOffset", "gridUnitSI",
                      "gridUnitDimension", "unitDimension", "timeOffset", "photonEnergy",
                      "zCoordinate"):
                if a not in m.attrs:
                    print(f"FAIL: {fn} required mesh attribute {a} MISSING")
                    FAILED = True
            if list(m.attrs["axisLabels"]) != [b"z", b"y", b"x"]:
                print(f"FAIL: {fn} axisLabels {m.attrs['axisLabels']}"); FAILED = True
            if list(m.attrs["unitDimension"]) != [1., 1., -3., -1., 0., 0., 0.]:
                print(f"FAIL: {fn} unitDimension is not the V/m exponents"); FAILED = True
            lam_file = H_PLANCK_EVS * C_LIGHT * E_CHARGE / float(m.attrs["photonEnergy"][0])
            if abs(lam_file - lam) > 1e-6 * lam:
                print(f"FAIL: {fn} photonEnergy wavelength {lam_file} vs {lam}"); FAILED = True
            if "unitSI" not in m["x"].attrs or abs(float(m["x"].attrs["unitSI"][0]) - 1) > 0:
                print(f"FAIL: {fn} component unitSI missing or not 1"); FAILED = True
            E_pmd = m["x"][0]                      # native complex128, slice 1 (y, x)

        E_gen, dsp = dfl_to_vperm(wd / fn.replace(".wf.h5", ".fld.h5"))
        d = float(np.max(np.abs(E_gen - E_pmd)) / np.max(np.abs(E_gen)))
        check(f"the converter's Genesis view agrees complex-wise ({fn})", d, TOL_RT_FORMATS)

    # Fortran read-back: import the .wf.h5, write_initial + load_only re-dumps it. The
    # field datasets must be identical to the file that was read.

    run(exe, wd, "h3rt", NML_IMPORT.format(lat="planar.bmad", root="h3rt",
        beam="h3both-final.beam.h5", field="h3both-final.wf.h5",
        extra="  write_initial = T\n  load_only = T\n"))
    check("Fortran openPMD read-back, re-dump dataset-identical",
          0.0 if fields_identical(wd / "h3rt-initial.wf.h5", wd / "h3both-final.wf.h5") else 1.0,
          0.5)

    # The Python Wavefront class, via its own Genesis path, agrees on the energy.

    sys.path.insert(0, args.pyrepo)
    try:
        from beamphysics.wavefront import Wavefront
    except ImportError as exc:
        print(f"FAIL: cannot import the Python Wavefront class from {args.pyrepo}: {exc}")
        print("      The openPMD path needs a beamphysics carrying "
              "beamphysics/wavefront/openpmd.py.")
        sys.exit(1)
    wf_py = Wavefront.from_genesis4(str(wd / "h3both-final.fld.h5"))
    E_gen, dsp = dfl_to_vperm(wd / "h3both-final.fld.h5")
    u_here = power_of(E_gen, dsp) * 1e-10 / C_LIGHT   # SS: dz = wavelength
    d = abs(wf_py.energy - u_here) / u_here
    check("Python Wavefront class energy (its Genesis path)", d, TOL_PY_ENERGY)

    # The Python class's own openPMD path (upstream): its reader must recover the
    # Genesis dump's complex values from the harmonic .wf.h5 through its own axis
    # handling, its writer must round-trip through its reader exactly, and the
    # Fortran reader must accept what it writes.

    wf_h = Wavefront.from_openpmd(wd / "h3both-final-h3.wf.h5")
    E_gen, dsp = dfl_to_vperm(wd / "h3both-final-h3.fld.h5")
    d = float(np.max(np.abs(wf_h.Ex[:, :, 0].T - E_gen)) / np.max(np.abs(E_gen)))
    check("Python openPMD reader (vs Genesis dump)", d, TOL_RT_FORMATS)
    wf_h.write_openpmd(wd / "pyrt.wf.h5")
    wf_rt = Wavefront.from_openpmd(wd / "pyrt.wf.h5")
    ok = (np.array_equal(wf_h.Ex, wf_rt.Ex) and wf_rt.Ey is None
          and wf_h.wavelength == wf_rt.wavelength)
    check("Python openPMD write/read round trip exact", 0.0 if ok else 1.0, 0.5)

    # Python writer to Fortran reader: rewrite both fields of the run through the
    # Python class, import that pair, and re-dump Genesis-format. The fundamental
    # must come back dataset-identical to the run's own Genesis dump.

    for src, dst in (("h3both-final.wf.h5", "pyw-final.wf.h5"),
                     ("h3both-final-h3.wf.h5", "pyw-final-h3.wf.h5")):
        Wavefront.from_openpmd(wd / src).write_openpmd(wd / dst)
    run(exe, wd, "pywrt", NML_IMPORT.format(lat="planar.bmad", root="pywrt",
        beam="h3both-final.beam.h5", field="pyw-final.wf.h5",
        extra="  write_initial = T\n  load_only = T\n"))
    u_py = fieldio.read_field(wd / "pywrt-initial.wf.h5")["u"]
    u_f = fieldio.read_field(wd / "h3both-final.wf.h5")["u"]
    check("Fortran reads the Python writer's file, re-dump identical",
          0.0 if np.array_equal(u_py, u_f) else 1.0, 0.5)

    # ------------------------------------------------------------------
    # 3. Deposit closed form: dark two-step restart from the bunched exit beam.

    print("== harmonic deposit vs the Bessel closed form ==")

    # A hard-seeded Bmad-only buncher: 1 GW over ~2 m drives the pendulum nonlinear,
    # so b3 is real (the tier's 4 m from 5 kW leaves b1 ~ 1e-3 and b3 in the noise).

    buncher = NML_TD.format(lat="mid.bmad", root="h3bunch", extra="")
    buncher = buncher.replace('  beam_init%distribution_type(3) = "GRID"\n', "")
    buncher = buncher.replace("  beam_init%grid(3)%x_min = -4e-10\n", "")
    buncher = buncher.replace("  beam_init%grid(3)%x_max = 4e-10\n", "  beam_init%sig_z = 0\n")
    buncher = buncher.replace("  window_length = 8e-10\n", "")
    buncher = buncher.replace("  window_sample = 1\n", "")
    buncher = buncher.replace("  shot_noise = T\n", "")
    buncher = buncher.replace("beam_init%n_particle = 2048", "beam_init%n_particle = 8192")
    buncher = buncher.replace("bunch_charge = 8.0e-15", "bunch_charge = 1.000692285594e-15")
    buncher = buncher.replace("seed_power = 1e4", "seed_power = 1e9")
    buncher = buncher.replace("grid_n_pts = 63", "grid_n_pts = 151")
    run(exe, wd, "h3bunch", buncher)

    shutil.copy(wd / "h3bunch-final.wf.h5", wd / "dark.wf.h5")
    with h5py.File(wd / "dark.wf.h5", "r+") as h5:
        h5[fieldio.MESH_PATH + "/x"][...] = 0.0

    # one step, dark start: the field at exit is exactly twice the source deposit of
    # the post-advance particles (a zero field diffracts to zero before the add), and
    # a one-step run's final particle dump is exactly the state the deposit read. So
    # the exit powers are predictable from the dump by the deposit sum itself -- the
    # Bessel fc(h) and the harmonic phase h*theta, with no evolution approximation.

    run(exe, wd, "h3dep", NML_IMPORT.format(lat="short.bmad", root="h3dep",
        beam="h3bunch-final.beam.h5", field="dark.wf.h5", extra=""))

    expect = bessel_p3_over_p1(wd / "h3dep-final.beam.h5", wd / "dark.wf.h5")
    with read_stats(wd / "h3dep.stats.h5") as st:
        p1 = float(st["field/total/power"][-1, 0])
        p3 = float(st["field/harm3/total/power"][-1, 0])
    d = abs(p3 / p1 / expect - 1)
    check("one-step deposit P3/P1 vs the Bessel fc + h*theta sum", d, TOL_DEPOSIT,
          note=f"[P3/P1 {p3/p1:.4e}]")

    # ------------------------------------------------------------------
    # 4. Thread identity on a TD harmonic run.

    print("== thread identity (TD, harmonics 1+3) ==")
    run(exe, wd, "h3t1", NML_TD.format(lat="planar.bmad", root="h3t1", extra=""), threads="1")
    run(exe, wd, "h3t8", NML_TD.format(lat="planar.bmad", root="h3t8", extra=""), threads="8")
    same = (wd / "h3t1.diag.txt").read_bytes() == (wd / "h3t8.diag.txt").read_bytes()
    same = same and fields_identical(wd / "h3t1-final-h3.wf.h5", wd / "h3t8-final-h3.wf.h5")
    check("1 vs 8 threads byte/dataset-identical", 0.0 if same else 1.0, 0.5)

    # ------------------------------------------------------------------
    # 5. Refusals.

    print("== refusals ==")
    base = NML_TD.format(lat="planar.bmad", root="rf", extra="{extra}")

    ok = refuse(exe, wd, "rf_anchor",
                base.replace("harmonics = 1, 3", "harmonics = 3").format(extra=""),
                "HARMONICS(1) MUST BE 1")
    print(f"--- refusal harmonics without the fundamental: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    (wd / "planar_uv.bmad").write_text("call, file = planar.bmad\n"
                                       "wiggler::*[TRACKING_METHOD] = fel_unaveraged\n")
    ok = refuse(exe, wd, "rf_unavg",
                NML_TD.format(lat="planar_uv.bmad", root="rf_unavg", extra=""),
                "UNAVERAGED")
    print(f"--- refusal harmonics + unaveraged element: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    ok = refuse(exe, wd, "rf_pol",
                NML_TD.format(lat="planar.bmad", root="rf_pol",
                              extra="  seed_polarization = 'y'\n"),
                "TWO LIVE POLARIZATIONS")
    print(f"--- refusal harmonics + two polarizations: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    shutil.copy(wd / "h3both-final.wf.h5", wd / "freq.wf.h5")
    with h5py.File(wd / "freq.wf.h5", "r+") as h5:
        m = h5["data/1/meshes/electricField"]
        del m.attrs["temporalDomain"]
        m.attrs["temporalDomain"] = np.bytes_("frequency")
    ok = refuse(exe, wd, "rf_freq", NML_IMPORT.format(lat="planar.bmad", root="rf_freq",
                beam="h3both-final.beam.h5", field="freq.wf.h5", extra=""),
                "ONLY THE time DOMAIN")
    print(f"--- refusal frequency-domain openPMD import: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    nomatch = NML_IMPORT.format(lat="planar.bmad", root="rf_match",
                                beam="h3both-final.beam.h5", field="H3-initial.wf.h5", extra="")
    nomatch = nomatch.replace("harmonics = 1, 3", "harmonics = 1, 5")
    nomatch = nomatch.replace('field_file = "H3-initial.wf.h5"',
                              'field_file = "H3-initial.wf.h5", "h3both-final-h3.wf.h5"')
    ok = refuse(exe, wd, "rf_match", nomatch, "MATCHES NO FIELD")
    print(f"--- refusal harmonic import matching no field: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    genharm = NML_IMPORT.format(lat="planar.bmad", root="rf_genh",
                                beam="h3both-final.beam.h5", field="H3-initial.wf.h5", extra="")
    genharm = genharm.replace('field_file = "H3-initial.wf.h5"',
                              'field_file = "H3-initial.wf.h5", "h3both-final-h3.fld.h5"')
    ok = refuse(exe, wd, "rf_genh", genharm, "MUST BE openPMD")
    print(f"--- refusal Genesis-format harmonic import: {'ok' if ok else '** FAIL **'}")
    FAILED = FAILED or not ok

    print("checks: " + ("FAIL" if FAILED else "PASS"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
