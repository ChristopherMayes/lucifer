#!/usr/bin/env python3
"""
Two-polarization checks (fel-physics.md, sec-field, the vector convention). The
radiation carries (Ex, Ey) when any FEL element is tilted, and an x-planar set followed
by a y-planar set must behave correctly. Symmetries and physics hold the checks, and no
reference file is involved:

  1. Rotation identity. An all-y line has every element at tilt = pi/2, and the probe is
     symmetric, with a round beam and no quads. The all-y line is fed the x/y-swapped
     beam from the driver's swap_beam_xy check instrument. The RNG draws its planes in
     sequence, so the generated beam itself is never swap-symmetric. The all-y line must
     reproduce the all-x line fed the original beam: total power identical, sizes
     swapped. Both FEL modes run. In the averaged mode the only asymmetry left is
     cos(pi/2) = 6e-17, so the level is near machine precision. The unaveraged mode is
     held at 1e-5, since its vector path carries a betatron-current radiation that the
     scalar path omits.
  2. Crossed undulator, an x set then a y set. The x field only diffracts through the y
     set, so its polarization-resolved dump power must match the same line with the y
     set replaced by an equal drift, a pure-diffraction reference. The y field grows
     from the microbunching of the x set, far above its dark floor. Bunching is
     longitudinal and blind to polarization: this is the crossed-polarized afterburner.
     Both FEL modes run, and the two modes must agree within 0.5 in |ln| of the y power.
  3. Ledger. The unaveraged time-dependent closure holds on the crossed line, since
     kick and deposit are exact duals per component.
  4. Thread identity. 1 and 8 threads are byte-identical on the crossed unaveraged run.
  5. Helical re-anchor. A helical run is forced onto the two-polarization path by a
     y-polarized seed of negligible power, which makes Ey live. The run must reproduce
     the scalar-path run's diag to TOL_HEL. The vector path evaluates the same physics
     through arithmetic in a different order.
     Handedness. A cold beam is driven by the scalar seed and by each circular pair of
     the same intensity. In both modes the pair (1, +i) must reproduce the scalar seed's
     energy modulation, and the pair (1, -i) must leave the beam unmodulated.
  6. Refusals. Tilt on a helical element. Tilt with the transcribed-Genesis maps.
  7. An imported field carrying Ey decides the run's planes, with or without the seed's
     polarization in the input file.

Run by the benchmark. Exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import h5py
import shutil
import pathlib
import re
import subprocess
import sys

import numpy as np

import fieldio
from checkrec import record

from lucifer.flat import to_toml as to_groups

FAILED = False
TOL_HEL = 1e-12  # The scalar envelope is the co-rotating pair. The two paths differ by 7.2e-15.

NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 2048
  beam.bunch_charge = 8.0e-15
  beam.z_distribution = "GRID"
  beam.z_min = -4e-10
  beam.z_max = 4e-10
  beam.sig_pz = 8.804506566858e-05
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 64
  grid_half_width = 2e-4
  window_length = 8e-10
  n_wavelength = 1
  shot_noise = T
  ran_seed = 777
  write_diag = T
{extra}&end
"""

UV = "wiggler::*[FEL_METHOD] = unaveraged\n"


def ss(nml):
    """Steady-state variant: a single slice, no slippage and no escape. This is the clean
    frame for the rotation and crossed-undulator physics. An 8-slice time-dependent window
    loses the whole field to slippage over this line, so these checks run in steady state
    while the ledger and thread checks stay time dependent."""
    out = nml.replace('  beam.z_distribution = "GRID"\n', "")
    out = out.replace("  beam.z_min = -4e-10\n", "")
    out = out.replace("  beam.z_max = 4e-10\n", "  beam.sig_z = 0\n")
    out = out.replace("  window_length = 8e-10\n", "")
    out = out.replace("  n_wavelength = 1\n", "")
    out = out.replace("  shot_noise = T\n", "")
    return out


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    record(name, value, tol, "max", ok)
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".toml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
        sys.exit(1)


def refuse(exe, wd, name, text, fragment):
    (wd / (name + ".toml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    return r.returncode != 0 and fragment in (r.stdout + r.stderr)


def diag(wd, root):
    fn = wd / f"{root}.diag.txt"
    ns = int(re.search(r"nslice = (\d+)", fn.open().read(400)).group(1))
    return np.loadtxt(fn).reshape(-1, ns, 12)


def split_outside(wd, root):
    """The window power outside the split angle at the last record that took the field
    angle moments, the same at the record before it, and the y plane's power at the
    last. On the crossed line the records that take the moments are the element ends,
    so the one before the last is the y set's entrance, where the field is the x plane
    alone and its outside power is the x plane's own."""
    with h5py.File(wd / f"{root}.stats.h5") as h:
        p_in = np.asarray(h["field/total/power_inside_angle"])
        p_tot = np.asarray(h["field/total/power"])
        p_y = np.asarray(h["field/y/power"])
    live = np.where(np.isfinite(p_in).all(axis=1))[0]
    out = p_tot.sum(axis=1) - p_in.sum(axis=1)
    return float(out[live[-1]]), float(out[live[-2]]), float(p_y[live[-1]].sum())


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
        # Averaged: machine precision, since both runs use the polarization-projected
        # coupling. Unaveraged: the x-run takes the scalar path, which by convention
        # deposits only the current in the wiggle direction. The y-run takes the vector
        # path, which also carries the tiny betatron-current radiation. That is a model
        # refinement of the vector path. The difference is 4.1e-6 here, and the tolerance
        # of 1e-5 admits it.
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
        # The split counts both planes. The y plane's share outside the split is the
        # exit's outside power less the x plane's own, read at the y set's entrance, over
        # the y power. A split taken on the x plane alone puts the whole y plane outside
        # and reads 1 here, whatever the y plane's shape.
        out_exit, out_x, p_y = split_outside(wd, f"p2c_{tag}")
        check(f"crossed, {mode}: the y plane's share outside the split, "
              "(outside at exit - the x plane's own) / Py", (out_exit - out_x) / p_y, 0.5,
              note="(a split on the x plane alone reads 1)")
        lnr_modes[mode] = (px, py)
    lnr = abs(np.log(lnr_modes["averaged"][1] / lnr_modes["unaveraged"][1]))
    check("crossed: averaged vs unaveraged y-power, |ln ratio| (priced)", lnr, 0.5)

    # 3. Ledger on a time-dependent crossed unaveraged run, with slippage and escape live.
    run(exe, wd, "p2t_uv", NML.format(lat="p2_c_uv.bmad", root="p2t_uv",
                                      extra="  keep_escaped_field = T\n"))
    led = np.loadtxt(wd / "p2t_uv.ledger.txt")
    etot = led[:, 1] + led[:, 2] + led[:, 4] - led[:, 5] + led[:, 6]
    turn = np.abs(np.diff(led[:, 2])).sum() + abs(led[-1, 4]) + abs(led[-1, 5]) + abs(led[-1, 6])
    check("crossed: unaveraged TD ledger closure, both components",
          np.abs(etot - etot[0]).max() / max(turn, 1e-300), 1e-3)

    # 3b. Both polarizations of the escaped light and of the pulse sit in one file each,
    #     as components x and y of its field record, with the y component's params beside
    #     the x component's.
    for what in ("escaped", "pulse"):
        comps = fieldio.components(wd / f"p2t_uv-{what}.wf.h5")
        check(f"crossed: the {what} file carries both polarizations (0 = yes)",
              0.0 if comps == ("x", "y") else 1.0, 0.5, note=f"[components {comps}]")
    with h5py.File(wd / "p2t_uv-escaped.wf.h5") as h5:
        g = h5["lucifer/escaped"]
        both = "wavefrontParams" in g and "wavefrontParamsY" in g \
            and g["wavefrontParams"].shape == g["wavefrontParamsY"].shape == (len(g["zTransmit"]), 25)
    check("crossed: the escaped file's params come per polarization, one row per plane (0 = yes)",
          0.0 if both else 1.0, 0.5)

    # 4. Thread identity on the time-dependent crossed unaveraged run.
    run(exe, wd, "p2t_uv1", NML.format(lat="p2_c_uv.bmad", root="p2t_uv1",
                                       extra="  keep_escaped_field = T\n"), threads="1")
    same = all((wd / f"p2t_uv{s}").read_bytes() == (wd / f"p2t_uv1{s}").read_bytes()
               for s in (".diag.txt", ".ledger.txt"))
    check("crossed: 1-thread vs 8-thread byte-identical (0 = yes)", 0.0 if same else 1.0, 0.5)

    # 5. Helical re-anchor: the scalar path against the vector path, on real physics.
    #    This is a dark time-dependent run with physical shot noise. A quiet start's
    #    dark power is the numerical floor and means nothing to compare. The helical
    #    quiver's spontaneous radiation is purely co-rotating. The scalar envelope holds
    #    it whole and the vector path splits it into (Ex, Ey). The total power must
    #    agree. A negligible y seed, 1e-30 W, forces Ey live.
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

    # 5b. Handedness. A total power cannot tell the two circular polarizations apart.
    #     The averaged vector once carried the counter-rotating pair while the unaveraged
    #     Lorentz force selected the co-rotating one, so a seed drove one mode and not
    #     the other. A cold, nearly chargeless beam is driven by a uniform seed as the
    #     scalar envelope and as each circular pair of the same total intensity. The rms
    #     energy modulation at the exit says which pair the mode couples to. Both modes
    #     must prefer (1, +i), and each must reproduce its scalar seed with it.
    print("--- helical handedness (seeded, cold, both modes):")
    (wd / "p2_hel_uv.bmad").write_text("call, file = spont_probe.bmad\n" + UV)
    #     The seed is strong and the beam cold enough that the modulation is the seed's:
    #     a uniform 1e9 W seed modulates by about 1 keV where the beam's own spread at
    #     sig_pz = 1e-12 is a few meV, so the counter-rotating pair's residue is read
    #     against the modulation and not against the spread it started with.
    cold = ss(NML).replace("seed_waist_size = 30e-6", "seed_waist_size = 1e-2").replace(
        "seed_power = 1e4", "seed_power = 1e9").replace(
        "beam.sig_pz = 8.804506566858e-05", "beam.sig_pz = 1e-12").replace(
        "beam.bunch_charge = 8.0e-15", "beam.bunch_charge = 8.0e-18")
    run(exe, wd, "p2s_src", cold.format(lat="p2_hel.bmad", root="p2s_src",
                                        extra="  write_initial = T\n  load_only = T\n"))

    def circular(src, dst, sign):
        """The scalar seed as the pair (1, sign i)/sqrt2 of the same total intensity."""
        shutil.copy(src, dst)
        with h5py.File(dst, "r+") as h5:
            m = h5["data"][sorted(h5["data"])[0]]["meshes/electricField"]
            x = m["x"][...] / np.sqrt(2.0)
            m["x"][...] = x
            y = m.create_dataset("y", data=sign * 1j * x)
            for k, v in m["x"].attrs.items():
                y.attrs[k] = v

    circular(wd / "p2s_src-initial.wf.h5", wd / "p2s_plus.wf.h5", +1)
    circular(wd / "p2s_src-initial.wf.h5", wd / "p2s_minus.wf.h5", -1)
    for mode, lat in (("averaged", "p2_hel.bmad"), ("unaveraged", "p2_hel_uv.bmad")):
        sig = {}
        for tag, ff in (("scalar", "p2s_src-initial.wf.h5"), ("plus", "p2s_plus.wf.h5"),
                        ("minus", "p2s_minus.wf.h5")):
            root = f"p2s_{mode[:2]}_{tag}"
            run(exe, wd, root, cold.format(lat=lat, root=root,
                extra=f'  beam.file = "p2s_src-initial.beam.h5"\n  field.file = "{ff}"\n'
                      '  load_mode = "keep"\n'))
            sig[tag] = float(diag(wd, root)[-1, 0, 7])
        check(f"handedness, {mode}: the pair (1, +i) reproduces the scalar seed's energy modulation",
              abs(sig["plus"] - sig["scalar"]) / sig["scalar"], 1e-6,
              note=f"({sig['plus']:.4f} vs {sig['scalar']:.4f} eV)")
        # The counter-rotating pair still drives the electron at the beat of the two
        # rotations. The beat averages to nothing over whole periods and leaves a residue
        # from the ramps: 1e-3 of the modulation in the unaveraged mode here, and 4e-6 in
        # the averaged mode, which carries no beat at all. Both are under the 1e-2 level.
        check(f"handedness, {mode}: the pair (1, -i) leaves the beam unmodulated",
              sig["minus"] / sig["scalar"], 1e-2, note=f"({sig['minus']:.3e} eV)")

    # 6. Refusals.
    (wd / "p2_hel_tilt.bmad").write_text("call, file = spont_probe.bmad\nUNDS[tilt] = 0.3\n")
    ok1 = refuse(exe, wd, "p2f1", NML.format(lat="p2_hel_tilt.bmad", root="p2f1", extra=""),
                 "TILT ON A HELICAL FEL ELEMENT")
    # The transcribed maps are a run switch, transport_model, so the input file asks for them.
    (wd / "p2_tr_tilt.bmad").write_text(
        "call, file = crossed_probe.bmad\nuse, CROSSED\n")
    ok2 = refuse(exe, wd, "p2f2", NML.format(lat="p2_tr_tilt.bmad", root="p2f2",
                 extra='  transport_model = "genesis"\n'), "KNOW NO TILT")
    check("refusal: tilt on helical (0 = yes)", 0.0 if ok1 else 1.0, 0.5)
    check("refusal: tilt with transcribed maps (0 = yes)", 0.0 if ok2 else 1.0, 0.5)

    # 7. An imported vector field decides the run's planes. When the state was set at
    #    lattice setup from the seed's polarization and the tilts alone, a file carrying
    #    Ey under an input file that stated no polarization left the run at one plane. The stats
    #    and the device allocated one plane and the first record wrote past it: an array
    #    bound in the debug build and a bus error in the production one. The source here
    #    writes a y-seeded eight-slice field on the untilted helical line. The restart
    #    from it omits the seed's polarization, against the control that states it.
    print("--- an imported vector field sets the polarization state (TD, helical):")
    run(exe, wd, "p2i_src", NML.format(lat="p2_hel.bmad", root="p2i_src",
                                       extra="  seed_polarization = 'y'\n  write_initial = T\n"
                                             "  load_only = T\n"))
    load = ('  beam.file = "p2i_src-initial.beam.h5"\n  field.file = "p2i_src-initial.wf.h5"\n'
            '  load_mode = "keep"\n')
    run(exe, wd, "p2i_def", NML.format(lat="p2_hel.bmad", root="p2i_def", extra=load))
    run(exe, wd, "p2i_y", NML.format(lat="p2_hel.bmad", root="p2i_y",
                                     extra=load + "  seed_polarization = 'y'\n"))
    same = (wd / "p2i_def.diag.txt").read_bytes() == (wd / "p2i_y.diag.txt").read_bytes()
    check("import: the restart without the seed's polarization runs as the one with it (0 = yes)",
          0.0 if same else 1.0, 0.5)
    py = dump_power(wd, "p2i_def-final.wf.h5", "y")
    px = dump_power(wd, "p2i_def-final.wf.h5", "x")
    check("import: the y plane is carried to the end (Py/Px floor)", 1.0 / max(py / max(px, 1e-300), 1e-300),
          1e3, note=f"(Py/Px = {py/max(px,1e-300):.3e})")

    if FAILED:
        print("two-polarization checks: FAIL")
        sys.exit(1)
    print("two-polarization checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
