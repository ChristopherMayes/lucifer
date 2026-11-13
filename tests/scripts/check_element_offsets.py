#!/usr/bin/env python3
"""
An FEL element's x_offset, y_offset, kx, hkick and vkick, against Bmad's own tracking
(fel-physics.md sec-element, validation.md val-element-offsets). Bmad only.

One short planar wiggler (offset_probe.bmad, 40 periods, Benchmark1-SASE energy and
aw) carries, in turn, nothing, an offset, a kx, a kick, and all three. A beam of 256
particles is tracked through each with the radiation off, a dark run with a quiet
start and no shot noise, so the orbits are the magnetic optics alone, four ways:

  bs   Bmad's own bmad_standard, track_a_wiggler. The averaged FEL mode's transverse
       maps are this routine's, flattened per step, so the two must agree to
       roundoff on what each attribute does.
  rk   Bmad's own runge_kutta through the planar field, the physics reference. The
       exit orbit carries the undulation, which the averaged mode leaves out and the
       hard-edged unaveraged run carries differently, so the comparison is made on
       differences: what the attribute changes at the exit, particle by particle.
  avg  The averaged FEL mode, whose maps take the offset, the kx strengths and the
       half kicks at the faces as track_a_wiggler does.
  uv   The unaveraged FEL mode with the hard-edged field, so that its field is the
       runge_kutta field: the vector potential of Bmad's planar model with kx, about
       the offset axis, with the kicks as a uniform field over the element.

For each attribute and each tracker, d = exit(attribute) - exit(none) per particle,
in x, y, px and py. The checks hold max |d_avg - d_bs| / max |d_bs| near roundoff,
and max |d_uv - d_rk| / max |d_rk| and max |d_avg - d_rk| / max |d_rk| at the model
level. The kick's d_rk itself is checked against the kick: px moves by hkick and x by
L hkick / 2. Two refusals are matched on their messages, kx on a helical element and
a pitched element.

Run by the validation suite. Exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import numpy as np

from checkrec import record
from lucifer.beam import read_slices
from lucifer.flat import to_toml as to_groups

import pool

LAMBDA0 = 1e-10
L_UND = 0.60
L_RAMP = 2 * 0.015                   # the default unaveraged ramp, 2 periods at each end
P0_MC = np.sqrt(11357.82**2 - 1)     # the dumps' px are gamma beta; Bmad's px/p0 divides by this
X_OFF, Y_OFF = 3e-5, -2e-5
KX = 2 * 2 * np.pi / 0.015           # 2 k_u: a strong roll-off, cos(kx x) = 1 - 3e-4 at 30 um
HKICK, VKICK = 2e-6, -3e-6

FAILED = False

NML = """! flat keys; routed into the input file's tables by lucifer.flat.to_toml
  lattice = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam.n_particle = 256
  beam.bunch_charge = 1e-18
  beam.sig_z = 0
  beam.sig_pz = 1e-12
  beam.a_norm_emit = 4e-7
  beam.b_norm_emit = 4e-7
  beamlet_size = 4
  seed_power = 0
  seed_waist_size = 30e-6
  grid_n_pts = 33
  grid_half_width = 2e-4
  shot_noise = F
  ran_seed = 4242
{extra}&end
"""

CASES = {
    "none": "",
    "offset": f"UNDS[x_offset] = {X_OFF}\nUNDS[y_offset] = {Y_OFF}\n",
    "kx": f"UNDS[kx] = {KX}\n",
    "kick": f"UNDS[hkick] = {HKICK}\nUNDS[vkick] = {VKICK}\n",
    "all": (f"UNDS[x_offset] = {X_OFF}\nUNDS[y_offset] = {Y_OFF}\nUNDS[kx] = {KX}\n"
            f"UNDS[hkick] = {HKICK}\nUNDS[vkick] = {VKICK}\n"),
}

# The four trackers: the lattice lines that select each, and the input keys each needs.
TRACKERS = {
    "bs": ("UNDS[fel_method] = off\nUNDS[tracking_method] = bmad_standard\n",
           "  reference_run = T\n"),
    "rk": ("UNDS[fel_method] = off\nUNDS[tracking_method] = runge_kutta\nUNDS[num_steps] = 400\n",
           "  reference_run = T\n"),
    "avg": ("", ""),
    "uv": ("UNDS[fel_method] = unaveraged\n", ""),
}


def check(name, value, bound, note=""):
    """One-sided: value <= bound."""
    global FAILED
    ok = value <= bound
    print(f"--- {name}: {value:.4e} (check <= {bound:.2e}) {note} {'ok' if ok else '** FAIL **'}")
    record(name, value, bound, "max", ok)
    if not ok:
        FAILED = True


def check_message(name, text, needle):
    global FAILED
    ok = needle in text
    print(f"--- {name}: {'matched' if ok else '** FAIL **'}")
    record(name, None, None, "bool", ok)
    if not ok:
        FAILED = True
        print(text[-1500:])


def run(exe, wd, name, text, threads="2", expect_fail=False):
    (wd / (name + ".toml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".toml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if expect_fail:
        return r
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2500:]}\n{r.stderr[-800:]}")
        sys.exit(1)
    return r


def exit_orbit(wd, root):
    """(x, y, px, py) of every particle at the exit, from the final dump, in load order."""
    s = read_slices(wd / f"{root}-final.beam.h5", wavelength=LAMBDA0)
    assert len(s) == 1, f"{root}: {len(s)} slices, one expected"
    return np.stack([s[0]["x"], s[0]["y"], s[0]["px"], s[0]["py"]])


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
    for lat in ("offset_probe.bmad", "spont_probe.bmad"):
        (wd / lat).write_bytes((latdir / lat).read_bytes())

    # One lattice file per tracker and case, each a call of the probe plus its lines.
    jobs = []
    for tr, (lines, extra) in TRACKERS.items():
        for case, attrs in CASES.items():
            lat = f"op_{tr}_{case}.bmad"
            (wd / lat).write_text("call, file = offset_probe.bmad\n" + lines + attrs)
            root = f"op_{tr}_{case}"
            jobs.append(lambda lat=lat, root=root, extra=extra:
                        run(exe, wd, root, NML.format(lat=lat, root=root, extra=extra)))
    pool.run_all(jobs, threads_per_job=2)

    orb = {(tr, case): exit_orbit(wd, f"op_{tr}_{case}") for tr in TRACKERS for case in CASES}
    d = {(tr, case): orb[(tr, case)] - orb[(tr, "none")] for tr in TRACKERS for case in CASES
         if case != "none"}
    for key in d:
        d[key][2:4] /= P0_MC                     # px, py as px/p0 and py/p0, Bmad's vec(2), vec(4)

    def rel(a, b, scale=1.0):
        """max |a - scale b| over max |b|, per coordinate pair (x, y) and (px, py), the larger."""
        out = 0.0
        for rows in ((0, 1), (2, 3)):
            num = np.abs(a[rows, :] - scale * b[rows, :]).max()
            den = np.abs(b[rows, :]).max()
            out = max(out, num / den)
        return out

    # The focusing attributes, offset and kx: the averaged mode against the map it was
    # flattened from, then both modes against Bmad's integrator. The unaveraged field's
    # sin^2 ramps hold the focusing integral at (1 - 5 l_ramp / 4 L) of the hard-edged
    # one (fel-physics.md sec-unaveraged), and the comparison takes that factor; the raw
    # ratio is printed beside it.
    f_ramp = 1 - 1.25 * L_RAMP / L_UND
    for case in ("offset", "kx"):
        print(f"--- {case}: max |d_rk| x,y {np.abs(d[('rk', case)][0:2]).max():.3e} m, "
              f"px,py {np.abs(d[('rk', case)][2:4]).max():.3e} rad")
        check(f"{case}: averaged against bmad_standard, max |d_avg - d_bs| / max |d_bs|",
              rel(d[("avg", case)], d[("bs", case)]), BOUNDS["avg_bs"][case])
        check(f"{case}: averaged against runge_kutta, max |d_avg - d_rk| / max |d_rk|",
              rel(d[("avg", case)], d[("rk", case)]), BOUNDS["avg_rk"][case])
        check(f"{case}: unaveraged against runge_kutta with the ramp factor, "
              f"max |d_uv - f d_rk| / max |d_rk|",
              rel(d[("uv", case)], d[("rk", case)], f_ramp), BOUNDS["uv_rk"][case],
              note=f"(f = {f_ramp:.4f}; without it {rel(d[('uv', case)], d[('rk', case)]):.3e})")

    # The kicks, against the attribute itself: px by hkick and x by L hkick / 2, and the
    # same in y, in both FEL modes. The y plane carries the natural focusing, so its two
    # move by a part in 1e4 from the drift values. Bmad's own two trackers are printed
    # for the record and not held: at release 20261023-0 bmad_standard applies a
    # wiggler's kick twice, once through offset_particle and once through the kick
    # folded into its multipoles, and runge_kutta's planar field carries no kick.
    for tr in ("avg", "uv"):
        dk = d[(tr, "kick")]
        check(f"kick: {tr} |d_px / hkick - 1|", abs(dk[2].mean() / HKICK - 1), BOUNDS["kick"][tr][0])
        check(f"kick: {tr} |d_x / (L hkick / 2) - 1|", abs(dk[0].mean() / (L_UND * HKICK / 2) - 1),
              BOUNDS["kick"][tr][1])
        check(f"kick: {tr} |d_py / vkick - 1|", abs(dk[3].mean() / VKICK - 1), BOUNDS["kick"][tr][2])
        check(f"kick: {tr} |d_y / (L vkick / 2) - 1|", abs(dk[1].mean() / (L_UND * VKICK / 2) - 1),
              BOUNDS["kick"][tr][3])
    for tr in ("bs", "rk"):
        dk = d[(tr, "kick")]
        print(f"--- kick: Bmad's {tr}, for the record: d_px / hkick {dk[2].mean() / HKICK:.4f}, "
              f"d_x / (L hkick / 2) {dk[0].mean() / (L_UND * HKICK / 2):.4f}")

    # All three at once compose: d_all against the sum of the three, in both modes.
    for tr in ("avg", "uv"):
        s3 = d[(tr, "offset")] + d[(tr, "kx")] + d[(tr, "kick")]
        check(f"all: {tr} max |d_all - (d_offset + d_kx + d_kick)| / max |d_all|",
              rel(d[(tr, "all")], s3) * np.abs(s3).max() / np.abs(d[(tr, "all")]).max(),
              BOUNDS["all"][tr])

    # Refusals: kx on a helical element, and a pitched element.
    (wd / "op_helical_kx.bmad").write_text("call, file = spont_probe.bmad\nUNDS[kx] = 100\n")
    r = run(exe, wd, "op_helical_kx", NML.format(lat="op_helical_kx.bmad", root="op_helical_kx", extra=""),
            expect_fail=True)
    check_message("refusal: kx on a helical element", r.stdout + r.stderr,
                  "KX ON A HELICAL FEL ELEMENT HAS NO MEANING")
    (wd / "op_pitch.bmad").write_text("call, file = offset_probe.bmad\nUNDS[x_pitch] = 1e-6\n")
    r = run(exe, wd, "op_pitch", NML.format(lat="op_pitch.bmad", root="op_pitch", extra=""),
            expect_fail=True)
    check_message("refusal: a pitched element", r.stdout + r.stderr,
                  "A PITCHED FEL ELEMENT IS NOT CARRIED")

    if FAILED:
        print("FAIL: element-offsets checks")
        sys.exit(1)
    print("--- element-offsets checks: all ok")


# Measured first with every bound at infinity, then set (validation.md,
# val-element-offsets). The averaged mode against bmad_standard carries the octupole
# term's step discretization alone. The two modes against runge_kutta are the model:
# the flattened periodic map and the ramped integration against Bmad's integrator.
BOUNDS = {
    "avg_bs": {"offset": 1e-4, "kx": 2e-5},
    "avg_rk": {"offset": 3e-3, "kx": 2e-2},
    "uv_rk": {"offset": 3e-3, "kx": 2e-2},
    "kick": {"avg": [1e-12, 1e-10, 2e-4, 2e-4], "uv": [1e-9, 1e-7, 2e-4, 2e-4]},
    "all": {"avg": 5e-2, "uv": 5e-2},
}

if __name__ == "__main__":
    main()
