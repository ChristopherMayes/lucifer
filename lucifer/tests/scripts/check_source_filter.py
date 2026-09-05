#!/usr/bin/env python3
"""
Check the source-term filter against Genesis4's own (fel-physics.md sec-source-filter).

The filter is a transcription, so the check is the one every transcribed path gets: both
codes run the same configuration from the same particles and the same field, and what is
left over is transcription fidelity. Genesis4 4.6.15 carries the filter as source_filter,
so the reference exists rather than having to be argued for. The chain is
Aramis-td-sase-filter.in, the pure-SASE deck with the filter on at Genesis4's own defaults
(xcut = ycut = sigmoid = 1), which puts the sigmoid's edge at half the grid's Nyquist
frequency. The same deck with the filter off is tdsase, recorded at 2.3e-6.

Three runs, each compared against that one reference:

  filtered    the filter on at the reference's own settings, over the whole line. This
              is the level.
  mutate      the sigmoid on the propagated field instead of on the source. Both
              suppress wide angles and both leave a plausible power curve, so only the
              comparison against Genesis4 separates them.
  xcut2       the sigmoid's edge at twice the cut. The shape is right and the position
              is wrong.

The last two, and the unfiltered run beside them, stop at the second undulator.

The two mutations must land far outside the tolerance the first run holds, or the check
reports that it cannot fail and stops.

Usage:

  check_source_filter.py --exe <lucifer> --workdir <dir> [--tol <rel>]
"""

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_fel import load_fortran_diag, load_genesis_out, load_nslice  # noqa: E402

ROOT = "AramisTDSASEF"

# The level is measured over the whole line. The three runs that must fail are stopped at
# the end of the second undulator, since each is wrong from its first step and a full line
# to prove it costs four times the section. The comparison then runs over the records both
# sides have, which is what load_fortran_diag and load_genesis_out already agree on.
SHORT = "UND##2"

NML = """&fel_params
  lat_file = "aramis.bmad"
  global%out_root = "{root}"
  global%interlude_model = "genesis"
  global%transport_model = "genesis"
  global%write_diag = T
  global%source_filter = {on}
  global%source_filter_xcut = {xcut}
  global%source_filter_mutate = {mutate}
{end}/
&fel_beam_init
  beam_file = "{dumps}-initial.beam.h5"
  beamlet_size = 8
/
&fel_wavefront_init
  field_file = "{dumps}-initial.wf.h5"
  wavefront_init%lambda0 = {lam:.12e}
  wavefront_init%window_sample = {sample}
/
"""


def window(par):
    """The wavelength and sample the reference chain sliced on, from its own dump."""
    with h5py.File(par) as f:
        lam = float(f["slicelength"][0])
        spacing = float(f["slicespacing"][0])
    sample = round(spacing / lam)
    assert abs(sample * lam - spacing) < 1e-9 * spacing, f"{par}: spacing is not a whole sample"
    return lam, sample


def run(exe, wd, root, lam, sample, on="T", xcut=1.0, mutate="F", track_end=""):
    """One lucifer run in the work directory. Returns its diag file."""
    nml = wd / f"{root}.nml"
    end = f'  global%track_end = "{track_end}"\n' if track_end else ""
    nml.write_text(NML.format(root=root, dumps=ROOT, lam=lam, sample=sample,
                              on=on, xcut=xcut, mutate=mutate, end=end))
    log = wd / f"fel-{root}.log"

    # Four threads, the same cap check_device.py takes. The two benchmark passes run this
    # section at the same time as three other jobs, and four unbounded runs a pass
    # oversubscribed the machine badly enough to kill an unrelated run in the next section.

    with log.open("w") as fh:
        r = subprocess.run([exe, nml.name], cwd=wd, stdout=fh, stderr=subprocess.STDOUT,
                           env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: lucifer {root} exited {r.returncode}; log tail:", file=sys.stderr)
        print(log.read_text()[-2000:], file=sys.stderr)
        return None
    return wd / f"{root}.diag.txt"


def level(diag, out, nslice):
    """Worst relative difference in power over every record and slice, and at the exit."""
    f = load_fortran_diag(str(diag), nslice)
    g = load_genesis_out(str(out), nslice)
    n = min(len(f["z"]), len(g["z"]))
    fp, gp = np.asarray(f["power"])[:n], np.asarray(g["power"])[:n]
    scale = np.maximum(np.abs(gp), np.abs(gp).max() * 1e-12)
    worst = float(np.max(np.abs(fp - gp) / scale))
    exit_rel = float(abs(fp[-1].sum() - gp[-1].sum()) / abs(gp[-1].sum()))
    return worst, exit_rel, float(gp[-1].sum()), float(fp[-1].sum())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--tol", type=float, default=1.0e-4)
    a = p.parse_args()

    wd = pathlib.Path(a.workdir).resolve()
    exe = str(pathlib.Path(a.exe).resolve())
    out = wd / f"{ROOT}.out.h5"
    par = wd / f"{ROOT}-initial.par.h5"
    for f in (out, par, wd / f"{ROOT}-initial.beam.h5", wd / f"{ROOT}-initial.wf.h5"):
        if not f.exists():
            print(f"FAIL: the filtered reference chain is missing {f.name}", file=sys.stderr)
            return 1

    lam, sample = window(par)
    nslice = load_nslice(str(par))
    ok = True

    diag = run(exe, wd, "sfilter", lam, sample)
    if diag is None:
        return 1
    worst, exit_rel, gexit, fexit = level(diag, out, nslice)
    print(f"  filter on, xcut = 1: exit power {fexit:.6e} W against Genesis4 {gexit:.6e} W, "
          f"relative {exit_rel:.3e}, worst over records {worst:.3e}")
    if worst > a.tol:
        print(f"FAIL: the filtered run differs from Genesis4 by {worst:.3e}, "
              f"above the tolerance {a.tol:.1e}", file=sys.stderr)
        ok = False

    # The filter has to be what moved the answer, or the run above proves nothing about it.

    diag_off = run(exe, wd, "sfilteroff", lam, sample, on="F", track_end=SHORT)
    if diag_off is None:
        return 1
    w_off, rel_off, _, fexit_off = level(diag_off, out, nslice)
    print(f"  filter off against the filtered reference: exit power {fexit_off:.6e} W, "
          f"relative {rel_off:.3e}, worst {w_off:.3e}")
    if w_off <= a.tol:
        print(f"FAIL: the unfiltered run also agrees with the filtered reference "
              f"({w_off:.3e}), so this check cannot see the filter", file=sys.stderr)
        ok = False

    for root, kw, what in (
        ("sfiltermut", dict(mutate="T"), "the sigmoid on the field instead of the source"),
        ("sfilterx2", dict(xcut=2.0), "the sigmoid's edge at twice the cut"),
    ):
        d = run(exe, wd, root, lam, sample, track_end=SHORT, **kw)
        if d is None:
            return 1
        w, rel, _, fe = level(d, out, nslice)
        verdict = "caught" if w > a.tol else "MISSED"
        print(f"  mutation, {what}: exit power {fe:.6e} W, relative {rel:.3e}, "
              f"worst {w:.3e}, {verdict}")
        if w <= a.tol:
            print(f"FAIL: mutation ({what}) is inside the tolerance {a.tol:.1e}, "
                  f"so the check cannot fail", file=sys.stderr)
            ok = False

    print("  source-filter checks: PASS" if ok else "  source-filter checks: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
