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

Then the edge's refusals, where the tolerance puts the edge relative to the angle it
protects, the containment of the seed inside the derived edge, and where the filter
reaches. That last one is measured rather than read off the line the run prints:
the filter is on by default, so an unaveraged element and a coherent source are no-ops
instead of refusals, and a no-op that works and a no-op that quietly filters print the same
line. Each pair runs with the switch on and off. The unaveraged line and the coherent source
must not move, the mixed line must, and the coherent case is the one that caught a fall
through the arming of the per-element state (FINDINGS 7.60).

Usage:

  check_source_filter.py --exe <lucifer> --workdir <dir> [--tol <rel>]
"""

import argparse
import pathlib
import re
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
  global%source_filter_ycut = {xcut}
  global%source_filter_width = 1
  global%source_filter_mutate = {mutate}
  slicing%n_wavelength = {sample}
{end}/
&fel_beam_init
  beam_file = "{dumps}-initial.beam.h5"
  beamlet_size = 8
/
&fel_wavefront_init
  field_file = "{dumps}-initial.wf.h5"
  wavefront_init%lambda0 = {lam:.12e}
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


REFUSE = """&fel_params
  lat_file = "aramis_1seg.bmad"
  global%out_root = "sfref"
  global%source_filter = T
{extra}/
&fel_beam_init
  beam_init%n_particle = 512
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 5e3
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_n_pts = 64
  wavefront_init%grid_half_width = 2e-4
/
"""


def refusals(exe, wd):
    """
    The three settings that must stop a run rather than quietly do something else.

    The soft width is the one that cost half the physical seed with nothing in the output
    saying so. The angle together with the grid-relative cuts has no meaning, since one is
    a property of the physics and the other of the mesh. A negative angle is not an angle.
    """
    ok = True
    for extra, want, what in (
        ("  global%source_filter_width = 1\n", "ON AXIS", "a width that attenuates the axis"),
        ('  global%source_filter_angle = 3e-6\n  global%source_filter_xcut = 1\n',
         "ARE BOTH SET", "the angle and the cut together"),
        ("  global%source_filter_angle = -1e-6\n", "MUST BE POSITIVE", "a negative angle"),
    ):
        (wd / "sfref.nml").write_text(REFUSE.format(extra=extra))
        r = subprocess.run([exe, "sfref.nml"], cwd=wd, capture_output=True, text=True)
        good = r.returncode != 0 and want in r.stdout
        print(f"  refused, {what}: {'ok' if good else 'MISSED'}")
        ok = ok and good
    return ok


def containment(exe, wd):
    """
    The default edge has to contain the mode, and the on-axis transmission cannot say so:
    a sigmoid that is 1 on axis can still clip a gain-guided mode whose divergence exceeds
    the beam's diffraction angle. This measures it instead, on a seeded steady-state run
    with the filter off, read at the first record rather than at the exit. At the exit the
    field is the mode plus the wide-angle emission of the beamlets, which is 10 percent of
    the power even at 2048 beamlets, so a containment measured there reports the artifact
    and not the edge. At the first record the field is the injected Gaussian and nothing
    else. Containing it is the stronger statement in any case, since the gain-guided mode
    that grows out of it is narrower than the seed.
    """
    deck = REFUSE.replace('lat_file = "aramis_1seg.bmad"', 'lat_file = "aramis.bmad"')
    deck = deck.replace("beam_init%n_particle = 512", "beam_init%n_particle = 2048")
    deck = deck.replace('global%out_root = "sfref"', 'global%out_root = "sfcont"')
    deck = deck.replace("  global%source_filter = T\n", "")
    (wd / "sfcont.nml").write_text(deck.format(extra=""))
    r = subprocess.run([exe, "sfcont.nml"], cwd=wd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: containment run exited {r.returncode}:\n{r.stdout[-1500:]}", file=sys.stderr)
        return False
    with h5py.File(wd / "sfcont.stats.h5") as h:
        g = h["field/total"]
        angle = float(np.ravel(g["split_angle"][()])[0])
        p_in = np.asarray(g["power_inside_angle"])
        p_tot = np.asarray(g["power"])
    live = np.isfinite(p_in).all(axis=1) & (p_tot.sum(axis=1) > 0)
    ir = np.where(live)[0]
    frac = float(p_in[ir[0]].sum() / p_tot[ir[0]].sum())
    good = frac > 0.99
    print(f"  seed containment at the default edge {angle:.3e} rad: {frac:.5f} of the "
          f"power inside (check > 0.99)  {'ok' if good else 'FAIL'}")
    return good


EDGE_LINE = re.compile(r"Protected angle\s+([0-9.eE+-]+)\s+rad.*?"
                       r"puts the edge at\s+([0-9.eE+-]+)\s+rad")


def edge_margin(tol, width):
    """The factor the code places the edge outside the protected angle by."""
    t = np.sqrt(1 - tol)
    return 1 - width * (np.log(t) + np.log1p(t) - np.log(tol))


def tolerance_run(exe, wd, extra, root="sftol"):
    """One short run with extra namelist lines. Returns its result and stdout."""
    deck = REFUSE.replace('global%out_root = "sfref"', f'global%out_root = "{root}"')
    (wd / f"{root}.nml").write_text(deck.format(extra=extra))
    return subprocess.run([exe, f"{root}.nml"], cwd=wd, capture_output=True, text=True,
                          env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})


def tolerance(exe, wd):
    """
    Where the edge goes relative to the angle it protects.

    The derived edge used to be the protected angle itself, which puts the sigmoid's
    half-amplitude point on it and passes a quarter of the source intensity there. The
    tolerance names the largest intensity loss allowed inside that angle and the edge
    moves out to hold it.

    Every check here has to fail against a run that reads source_filter_tolerance and
    then derives the old edge anyway, so each one compares the edge the run printed
    against the edge the tolerance asks for. An unchanged tier digit cannot do that
    work: the default reproduces the old edge to the bit, which is the point of it.
    """
    ok = True
    width = 0.05

    # The default and an explicit 0.75 are the edge in use before the tolerance existed,
    # and the two smaller tolerances have to move it by the margin and by nothing else.

    protected = None
    for extra, tol, what in (
        ("", 0.75, "left unset"),
        ("  global%source_filter_tolerance = 0.75\n", 0.75, "stated as 0.75"),
        ("  global%source_filter_tolerance = 0.2\n", 0.2, "stated as 0.2"),
        ("  global%source_filter_tolerance = 0.01\n", 0.01, "stated as 0.01"),
    ):
        r = tolerance_run(exe, wd, extra)
        m = EDGE_LINE.search(r.stdout.replace("\n", " ")) if r.returncode == 0 else None
        if m is None:
            print(f"FAIL: the tolerance {what} did not derive an edge "
                  f"(exit {r.returncode}):\n{r.stdout[-1500:]}", file=sys.stderr)
            return False
        th_p, edge = float(m.group(1)), float(m.group(2))
        if protected is None:
            protected = th_p
        want = th_p / edge_margin(tol, width)

        # The transmission the edge actually delivers at the protected angle, which is
        # what the tolerance was asked for.
        got_loss = 1 - (1 / (1 + np.exp((th_p / edge - 1) / width)))**2
        good = (abs(edge - want) <= 2e-4 * want and abs(th_p - protected) <= 1e-9 * protected
                and abs(got_loss - tol) <= 2e-3)
        print(f"  tolerance {what}: protected {th_p:.5e} rad, edge {edge:.5e} rad against "
              f"{want:.5e}, loss there {got_loss:.4f} against {tol}  "
              f"{'ok' if good else 'FAIL'}")
        ok = ok and good

    # The reported split is the acceptance, so it must not follow the edge out. A scan of
    # the tolerance that moved it would compare each run against a different question.

    angles = []
    for tol in (0.75, 0.2, 0.01):
        tolerance_run(exe, wd, f"  global%source_filter_tolerance = {tol}\n")
        with h5py.File(wd / "sftol.stats.h5") as h:
            angles.append(float(np.ravel(h["field/total/split_angle"][()])[0]))
    spread = max(angles) - min(angles)
    good = spread <= 1e-9 * angles[0]
    print(f"  the reported split across that scan: {angles[0]:.5e} rad, spread "
          f"{spread:.1e}  {'ok' if good else 'FAIL'}")
    ok = ok and good

    # Both ends of the range are singular. At one the margin runs to positive infinity and
    # the edge collapses onto the axis, which a test on the margin's sign would take, so
    # the range is checked before anything takes a logarithm.

    for extra, want, what in (
        ("  global%source_filter_tolerance = 0\n", "MUST LIE STRICTLY", "a tolerance of zero"),
        ("  global%source_filter_tolerance = 1\n", "MUST LIE STRICTLY", "a tolerance of one"),
        ("  global%source_filter_tolerance = 1.5\n", "MUST LIE STRICTLY", "a tolerance above one"),
        ("  global%source_filter_tolerance = -0.1\n", "MUST LIE STRICTLY", "a negative tolerance"),
        ("  global%source_filter_tolerance = 0.01\n  global%source_filter_width = 0.217\n",
         "IS TOO LARGE FOR", "a width too large for the tolerance"),
    ):
        r = tolerance_run(exe, wd, extra)
        good = r.returncode != 0 and want in r.stdout
        print(f"  refused, {what}: {'ok' if good else 'MISSED'}")
        ok = ok and good

    # A tolerance below about 2.2e-16 rounds 1 - tolerance to one, so the amplitude rounds
    # to one and 1 - amplitude to zero. At a width of 0.01 that tolerance still has a
    # margin of 0.60 and is an edge the run can build, so forming 1 - amplitude directly
    # would divide by zero on an input the run has to accept.

    r = tolerance_run(exe, wd, "  global%source_filter_tolerance = 1e-17\n"
                               "  global%source_filter_width = 0.01\n")
    m = EDGE_LINE.search(r.stdout.replace("\n", " ")) if r.returncode == 0 else None
    if m is None:
        print(f"FAIL: a tolerance of 1e-17 at a width of 0.01 is feasible and did not "
              f"derive an edge (exit {r.returncode}):\n{r.stdout[-1500:]}", file=sys.stderr)
        ok = False
    else:
        th_p, edge = float(m.group(1)), float(m.group(2))
        want = th_p / edge_margin(1e-17, 0.01)
        good = abs(edge - want) <= 2e-4 * want
        print(f"  a tolerance of 1e-17 at a width of 0.01: edge {edge:.5e} rad against "
              f"{want:.5e}  {'ok' if good else 'FAIL'}")
        ok = ok and good

    # A stated edge is a stated edge. Neither the angle nor the grid-relative cuts go
    # through the derivation, and a run that applied the tolerance to them would move an
    # edge the deck placed itself.

    for extra, what in (
        ("  global%source_filter_tolerance = 0.01\n  global%source_filter_angle = 5e-6\n",
         "a stated angle"),
        ("  global%source_filter_tolerance = 0.01\n  global%source_filter_xcut = 1\n"
         "  global%source_filter_ycut = 1\n  global%source_filter_width = 1\n",
         "stated cuts"),
    ):
        r = tolerance_run(exe, wd, extra)
        good = (r.returncode == 0 and "is bypassed" in r.stdout
                and EDGE_LINE.search(r.stdout.replace("\n", " ")) is None)
        print(f"  the tolerance is bypassed by {what}: {'ok' if good else 'FAIL'}")
        ok = ok and good

    return ok


SCOPE = """&fel_params
  lat_file = "{lat}"
  global%out_root = "{root}"
  global%source_filter = {on}
{extra}/
&fel_beam_init
  beam_init%n_particle = 512
  beam_init%bunch_charge = 1.000692285594e-15
  beam_init%sig_z = 0
  beam_init%sig_pz = 8.804506566858e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
/
&fel_wavefront_init
  wavefront_init%lambda0 = 1e-10
  wavefront_init%seed_power = 5e3
  wavefront_init%seed_waist_size = 30e-6
  wavefront_init%grid_n_pts = 64
  wavefront_init%grid_half_width = 2e-4
/
"""


def exit_power(wd, root):
    """Total field power at the last record of a run."""
    with h5py.File(wd / f"{root}.stats.h5") as h:
        return float(np.asarray(h["field/total/power"])[-1].sum())


def scope(exe, wd):
    """
    Where the filter reaches, measured rather than read off a printed line.

    The filter is applied to the source term the averaged deposit builds, so an
    unaveraged element has nothing for it to act on and a coherent source is an analytic
    Gaussian that carries no wide-angle content. Both are no-ops rather than refusals,
    since the filter is on by default and a deck that never asked for it must still run.
    A no-op that is announced and a no-op that is silent print the same line, so each
    pair here is run with the filter on and off and the exit powers compared. The mixed
    line is the one that has to move: it holds two averaged segments around one
    unaveraged segment, so a filter that reached neither would leave it unchanged.
    """
    ok = True
    for lat, extra, want_move, what in (
        ("aramis_1seg_unavg.bmad", "", False, "an unaveraged line"),
        ('  global%source_model = "coherent"\n', "", False, "a coherent source"),
        ("unavg_sandwich.bmad", "", True, "a mixed line, two averaged segments of three"),
    ):
        # The coherent case is the one row whose first field is a setting, not a lattice.
        if lat.startswith("  global"):
            lat, extra = "aramis_1seg.bmad", lat
        p = {}
        for on in ("T", "F"):
            root = f"sfscope{on}"
            (wd / f"{root}.nml").write_text(SCOPE.format(lat=lat, root=root, on=on, extra=extra))

            # Four threads, the cap run() takes above and for the same reason: these six
            # runs land in the middle of a pass that shares the machine with four other
            # jobs, and an unbounded one here took an unrelated section down with it.

            r = subprocess.run([exe, f"{root}.nml"], cwd=wd, capture_output=True, text=True,
                               env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
            if r.returncode != 0:
                print(f"FAIL: {what} with global%source_filter = {on} exited "
                      f"{r.returncode}:\n{r.stdout[-1500:]}", file=sys.stderr)
                return False
            p[on] = exit_power(wd, root)
        rel = abs(p["T"] - p["F"]) / abs(p["F"])
        moved = rel > 1e-9
        good = moved == want_move
        verb = "moves" if want_move else "is untouched"
        print(f"  the filter {verb} on {what}: exit power {p['T']:.6e} W on against "
              f"{p['F']:.6e} W off, relative {rel:.3e}  {'ok' if good else 'FAIL'}")
        ok = ok and good
    return ok


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

    ok = refusals(exe, wd) and ok
    ok = tolerance(exe, wd) and ok
    ok = containment(exe, wd) and ok
    ok = scope(exe, wd) and ok

    print("  source-filter checks: PASS" if ok else "  source-filter checks: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
