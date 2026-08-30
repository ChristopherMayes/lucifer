#!/usr/bin/env python3
"""
Two-polarization checks (manual sec:field vector convention): the radiation carries
(Ex, Ey) when any FEL element is tilted, and an x-planar set followed by a y-planar
set does the right thing. Held by symmetries and physics, not reference files:

  1. Rotation identity: an all-y line (every element tilt = pi/2; symmetric probe --
     round beam, no quads) fed the x/y-SWAPPED beam (the driver's swap_beam_xy check
     instrument: the RNG draws its planes sequentially, so the generated beam itself
     is never swap-symmetric) must reproduce the all-x line fed the original beam:
     total power identical, sizes swapped. Both FEL modes. The only asymmetry left is
     cos(pi/2) = 6e-17, so the level is near machine precision.
  2. Crossed undulator (x-set then y-set): the x field only diffracts through the
     y set -- its polarization-resolved dump power must match the same line with the
     y set replaced by an equal drift (pure-diffraction reference) -- while the y
     field grows from the x-set's microbunching (bunching is longitudinal and
     polarization-blind: the crossed-polarized afterburner), far above its dark
     floor. Both FEL modes, and the two modes agree at a priced level.
  3. Ledger: the unaveraged TD closure holds on the crossed line (kick/deposit are
     exact duals per component).
  4. Thread identity: 1 vs 8 threads byte-identical on the crossed unaveraged run.
  5. Helical re-anchor: a helical run forced onto the two-polarization path (a
     y-polarized zero-power seed makes Ey live) reproduces the scalar-path run's
     diag -- measured level (the vector path evaluates the same physics through
     differently-ordered arithmetic).
  6. Refusals: tilt on a helical element; tilt with the transcribed-Genesis maps.

Run by the benchmark harness; exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

import numpy as np

import fieldio

from nml import to_groups

FAILED = False
TOL_HEL = 1e-12  # measured 7.2e-15: the scalar envelope IS the co-rotating pair

NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 2048
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.804506566858e-05
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 63
  grid_half_width = 2e-4
  window_length = 8e-10
  window_sample = 1
  shotnoise = T
  ran_seed = 777
  write_diag = T
{extra}&end
"""

UV = "fel_unaveraged = 1\nwiggler::*[FEL_TRACKING] = fel_unaveraged\n"


def ss(nml):
    """Steady-state variant: single slice, no slippage, no escape -- the clean frame
    for the rotation and crossed-undulator physics (an 8-slice TD window loses the
    whole field to slippage over this line; measured, which is why these checks are
    SS while the ledger/thread checks stay TD)."""
    out = nml.replace('  beam_init%distribution_type(3) = "GRID"\n', "")
    out = out.replace("  beam_init%grid(3)%x_min = -4e-10\n", "")
    out = out.replace("  beam_init%grid(3)%x_max = 4e-10\n", "  beam_init%sig_z = 0\n")
    out = out.replace("  window_length = 8e-10\n", "")
    out = out.replace("  window_sample = 1\n", "")
    out = out.replace("  shotnoise = T\n", "")
    return out


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


def diag(wd, root):
    fn = wd / f"{root}.diag.txt"
    ns = int(re.search(r"nslice = (\d+)", fn.open().read(400)).group(1))
    return np.loadtxt(fn).reshape(-1, ns, 12)


def dump_power(wd, fname, component="x"):
    """Total power [W] in one polarization of a field dump, summed over every slice.

    One openPMD file carries both transverse polarizations as components of the
    electricField record, so the two are read out of one file rather than out of a file
    each. A component the file does not carry is no power: the tracker allocates the
    second component only for a field that can have one."""
    if component not in fieldio.components(wd / fname):
        return 0.0
    f = fieldio.read_field(wd / fname, component)
    return float(fieldio.field_power(f["u"], f["dx"], f["dy"]).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--latdir", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    latdir = pathlib.Path(args.latdir).resolve()
    wd = pathlib.Path(args.workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    for lat in ("crossed_probe.bmad", "spont_probe.bmad"):
        (wd / lat).write_bytes((latdir / lat).read_bytes())

    (wd / "p2_y.bmad").write_text("call, file = crossed_probe.bmad\nuse, YLINE\n")
    (wd / "p2_c.bmad").write_text("call, file = crossed_probe.bmad\nuse, CROSSED\n")
    (wd / "p2_x_uv.bmad").write_text("call, file = crossed_probe.bmad\n" + UV)
    (wd / "p2_y_uv.bmad").write_text("call, file = crossed_probe.bmad\n" + UV + "use, YLINE\n")
    (wd / "p2_c_uv.bmad").write_text("call, file = crossed_probe.bmad\n" + UV + "use, CROSSED\n")
    (wd / "p2_xd.bmad").write_text(
        "call, file = crossed_probe.bmad\nDEQ: pipe, l = 0.60\n"
        "XDRIFT: line = (UNDX, D1, DEQ)\nuse, XDRIFT\n")
    (wd / "p2_xd_uv.bmad").write_text(
        "call, file = crossed_probe.bmad\n" + UV +
        "DEQ: pipe, l = 0.60\nXDRIFT: line = (UNDX, D1, DEQ)\nuse, XDRIFT\n")

    # 1. Rotation identity.
    print("--- rotation identity (all-x vs all-y line with the swapped beam, SS):")
    for mode, xlat, ylat, tag in (("averaged", "crossed_probe.bmad", "p2_y.bmad", "av"),
                                  ("unaveraged", "p2_x_uv.bmad", "p2_y_uv.bmad", "uv")):
        run(exe, wd, f"p2x_{tag}", ss(NML.format(lat=xlat, root=f"p2x_{tag}", extra="")))
        run(exe, wd, f"p2y_{tag}", ss(NML.format(lat=ylat, root=f"p2y_{tag}",
                                      extra="  swap_beam_xy = T\n  seed_polarization = 'y'\n")))
        dx_, dy_ = diag(wd, f"p2x_{tag}"), diag(wd, f"p2y_{tag}")
        n = min(len(dx_), len(dy_))
        pwr = np.abs(dx_[:n, :, 2] - dy_[:n, :, 2]).max() / dx_[:n, :, 2].max()
        sw = max(np.abs(dx_[-1, :, 8] / dy_[-1, :, 9] - 1).max(),
                 np.abs(dx_[-1, :, 9] / dy_[-1, :, 8] - 1).max())
        # Averaged: machine precision (both runs use the pol-projected coupling).
        # Unaveraged: the x-run takes the scalar path (deposits only the wiggle-
        # direction current, by convention) while the y-run takes the vector path,
        # which also carries the tiny betatron-current radiation -- a real model
        # refinement of the vector path, measured 4.1e-6 here, tolerated not hidden.
        tol_p = 1e-8 if mode == "averaged" else 1e-5
        check(f"rotation, {mode}: power |x-run - y-run| / max, all records", pwr, tol_p)
        check(f"rotation, {mode}: beam-size swap identity at exit", sw, 1e-8)

    # 2. Crossed undulator.
    print("--- crossed undulator (x-set then y-set; x-seeded; SS):")
    lnr_modes = {}
    for mode, clat, dlat, tag in (("averaged", "p2_c.bmad", "p2_xd.bmad", "av"),
                                  ("unaveraged", "p2_c_uv.bmad", "p2_xd_uv.bmad", "uv")):
        run(exe, wd, f"p2c_{tag}", ss(NML.format(lat=clat, root=f"p2c_{tag}", extra="")))
        run(exe, wd, f"p2r_{tag}", ss(NML.format(lat=dlat, root=f"p2r_{tag}", extra="")))
        px = dump_power(wd, f"p2c_{tag}-final.wf.h5", "x")
        py = dump_power(wd, f"p2c_{tag}-final.wf.h5", "y")
        p_ref = dump_power(wd, f"p2r_{tag}-final.wf.h5", "x")
        iso = abs(np.log(px / p_ref))
        check(f"crossed, {mode}: x-field gain ISOLATION through the y set, |ln(Px/P_drift)|",
              iso, 5e-2, note="(the y set must only diffract Ex)")
        check(f"crossed, {mode}: afterburner floor, Px/Py (y must light up)",
              px / max(py, 1e-300), 2e2, note=f"(Py/Px = {py/max(px,1e-300):.3e})")
        lnr_modes[mode] = (px, py)
    lnr = abs(np.log(lnr_modes["averaged"][1] / lnr_modes["unaveraged"][1]))
    check("crossed: averaged vs unaveraged y-power, |ln ratio| (priced)", lnr, 0.5)

    # 3. Ledger on a TD crossed unaveraged run (slippage + escape live).
    run(exe, wd, "p2t_uv", NML.format(lat="p2_c_uv.bmad", root="p2t_uv",
                                      extra="  keep_escaped_field = T\n"))
    led = np.loadtxt(wd / "p2t_uv.ledger.txt")
    etot = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turn = np.abs(np.diff(led[:, 2])).sum() + abs(led[-1, 4]) + abs(led[-1, 5]) + abs(led[-1, 6])
    check("crossed: unaveraged TD ledger closure, both components",
          np.abs(etot - etot[0]).max() / max(turn, 1e-300), 1e-3)

    # 4. Thread identity on the TD crossed unaveraged run.
    run(exe, wd, "p2t_uv1", NML.format(lat="p2_c_uv.bmad", root="p2t_uv1",
                                       extra="  keep_escaped_field = T\n"), threads="1")
    same = all((wd / f"p2t_uv{s}").read_bytes() == (wd / f"p2t_uv1{s}").read_bytes()
               for s in (".diag.txt", ".ledger.txt"))
    check("crossed: 1-thread vs 8-thread byte-identical (1 = yes)", 0.0 if same else 1.0, 0.5)

    # 5. Helical re-anchor: scalar path vs the vector path, on real physics -- a dark
    #    TD run with physical shot noise (a quiet start's dark power is the numerical
    #    floor, meaningless to compare). The helical quiver's spontaneous radiation is
    #    purely co-rotating, which the scalar envelope holds whole and the vector path
    #    splits into (Ex, Ey). The total power must agree. Ey is forced live by a
    #    negligible y seed (1e-30 W).
    (wd / "p2_hel.bmad").write_text("call, file = spont_probe.bmad\n")
    run(exe, wd, "p2h_v0", NML.format(lat="p2_hel.bmad", root="p2h_v0",
                                      extra="  seed_polarization = 'y'\n").replace(
                                      "seed_power = 1e4", "seed_power = 1e-30"))
    run(exe, wd, "p2h_s0", NML.format(lat="p2_hel.bmad", root="p2h_s0", extra="").replace(
                                      "seed_power = 1e4", "seed_power = 1e-30"))
    ds_, dv_ = diag(wd, "p2h_s0"), diag(wd, "p2h_v0")
    n = min(len(ds_), len(dv_))
    hel = np.abs(dv_[:n, :, 2] - ds_[:n, :, 2]).max() / ds_[:n, :, 2].max()
    check("helical re-anchor: vector vs scalar path, shot-noise power", hel, TOL_HEL,
          note="(co-rotating radiation, two representations)")

    # 6. Refusals.
    (wd / "p2_hel_tilt.bmad").write_text("call, file = spont_probe.bmad\nUNDS[tilt] = 0.3\n")
    ok1 = refuse(exe, wd, "p2f1", NML.format(lat="p2_hel_tilt.bmad", root="p2f1", extra=""),
                 "TILT ON A HELICAL FEL ELEMENT")
    (wd / "p2_tr_tilt.bmad").write_text(
        "call, file = crossed_probe.bmad\nfel_transcribed = -1\n"
        "wiggler::*[FEL_TRACKING] = fel_transcribed\nuse, CROSSED\n")
    ok2 = refuse(exe, wd, "p2f2", NML.format(lat="p2_tr_tilt.bmad", root="p2f2", extra=""),
                 "KNOW NO TILT")
    check("refusal: tilt on helical, by name (1 = yes)", 0.0 if ok1 else 1.0, 0.5)
    check("refusal: tilt with transcribed maps, by name (1 = yes)", 0.0 if ok2 else 1.0, 0.5)

    if FAILED:
        print("two-polarization checks: FAIL")
        sys.exit(1)
    print("two-polarization checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
