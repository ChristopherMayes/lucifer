#!/usr/bin/env python3
"""
Compare the Bmad FEL tracker (lucifer) against Genesis 1.3 Version 4, over the
Benchmark1-SASE configuration: steady state and time dependent with slippage. Run through
run_fel_benchmark.sh, which produces all the inputs; this script only reads and compares.

Three tiers, each with its own tolerance sized to what it measures. The comparison
floor is set by fundamental constants: the tracker uses Bmad's (Z0 = mu_0*c =
376.7303...), Genesis carries a truncated impedance (376.73, 8.3e-7 relative), and that
difference enters the coupling and compounds through gain. Transcribing Genesis's own constants
agreed at transcription
level (tier1 2.8e-11, tier2_genesis 5.9e-8, recorded in doc/validation.md); after that
validation was banked the code moved to Bmad constants by decision, and the tiers now
measure against the constants floor.

  tier1         One undulator segment, no interludes. Isolates the FEL core: push,
                deposition, field solve. Observed ~2e-6 (final field), consistent with
                the impedance difference at one segment's growth.

  tier2_genesis The full 6-FODO line with the transcribed Genesis interlude model.
                Observed ~2e-5, the constants difference compounded through full gain.

  tier2_bmad    The full line with the seam: Bmad tracks the interludes, wavefront_drift
                moves the field, theta advances by the exact mapping from Bmad's z. This
                differs from Genesis by a real transport model difference -- Genesis
                samples the path length term px^2+py^2 once at quad mid-element
                (TrackBeam half step, then the theta step, then the second half), Bmad
                integrates it exactly through the quad map -- measured at ~2e-3 rad of
                bunching phase per quadrupole and ~1.3e-2 of power at saturation. The
                tier exists to demonstrate the seam works; its tolerance brackets the
                measured model difference to catch regressions, and tier2_genesis is
                what prices the difference.

Each tier compares the per-record power and bunching curves, and the final field and
particle dumps element by element. Particle ordering is preserved by both codes (no
sorting happens in steady state), so the dumps compare particle by particle.

The two codes write different formats and nothing is converted on disk to compare them.
The tracker writes openPMD: a .beam.h5 particle file carrying one patch per slice, and a
.wf.h5 field file in V/m. Genesis writes its own .par.h5 and .fld.h5. beamio and fieldio
read either format into the same arrays in memory, scaling Genesis's field amplitude to
V/m on the way in, and the comparison happens there. The window the tracker's beam file
does not state, the wavelength and the slice spacing, comes from the Genesis dump it is
being compared against.

A fourth check compares Fortran against itself: tier1 rerun with every particle split
into two coincident copies of weights w/3 and 2w/3 (split_weights = T). Collective
observables are linear in the weights, so the curves and the final field must be
unchanged to round-off. This is the only test of the weighted code paths with nonuniform
weights -- the Genesis dump format carries no weights, so every Genesis comparison sees
the uniform case, where a bug like using one particle's weight for all is invisible.

Three time-dependent tiers mirror the steady-state ones over 32 slices with slippage
active (Aramis-td.in: sample = 3, shot noise on, both codes starting from the same
multi-slice dumps):

  td1           One undulator segment: the FEL core plus slippage -- accumulation,
                threshold, rotation, zero fill, the unguarded end-of-lattice
                autophasing. Constants-floor agreement expected.

  td2_genesis   The full line, transcribed Genesis interludes: adds the drift
                autophasing schedule through twelve undulator/interlude alternations.

  td2_bmad      The full line through the Bmad seam: the transport model difference of
                tier2_bmad, now with slippage interleaved.

  tdsase        The full line, dark start (power = 0), shot noise on: pure SASE, both
                codes tracking the identical noisy beam and identical zero field from
                shared dumps (Aramis-td-sase.in). The one combination the seeded TD
                tiers and the statistical startup check leave uncovered -- the whole
                startup-from-noise path, compared deterministically and elementwise.
                Transcribed interludes, so the check sits at the constants floor.
                plot_fel_compare.py overlays the two codes' curves for any tier.

Per-slice curves compare as full (nrecords, nslice) arrays in time-window order; final
field and particle dumps compare every slice. The power denominator is floored at 1e-9
of the global peak so a freshly slipped-in vacuum slice cannot divide by zero.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

import beamio
import fieldio


def load_genesis_out(fn, nslice=1):
    """
    Genesis .out.h5 per-slice arrays are (nrecords, nslice), slices in time-window
    order. Steady state (nslice=1) returns 1D curves.
    """
    sl = slice(None) if nslice > 1 else 0
    with h5py.File(fn) as h5:
        return {
            "z": h5["Lattice/zplot"][:],
            "power": h5["Field/power"][:, sl],
            "bunching": h5["Beam/bunching"][:, sl],
            "bunchingphase": h5["Beam/bunchingphase"][:, sl],
            "energy": h5["Beam/energy"][:, sl],
            "xsize": h5["Beam/xsize"][:, sl],
            "ysize": h5["Beam/ysize"][:, sl],
        }


def load_nslice(fn):
    """The window's slice count, from either code's particle dump."""
    return beamio.n_slice(fn)


def load_fortran_diag(fn, nslice=1):
    """
    The Fortran diag file: one row per slice per record, columns
    z, slice, power, on_axis, bunching, bunching_phase, mean_energy [eV],
    sigma_energy [eV], sigma_x, sigma_y. Steady state (nslice=1) returns 1D curves;
    time dependent returns (nrecords, nslice) arrays, slices in time-window order,
    matching the Genesis .out.h5 layout. The energy column converts to gamma here
    (Genesis's Beam/energy is gamma; the diag carries Bmad-convention eV).
    """
    M_ELECTRON = 0.51099895069e6
    d = np.loadtxt(fn)
    if nslice > 1:
        assert d.shape[0] % nslice == 0, f"{fn}: {d.shape[0]} rows not divisible by nslice={nslice}"
        d = d.reshape(-1, nslice, d.shape[1])
        return {
            "z": d[:, 0, 0], "power": d[:, :, 2], "bunching": d[:, :, 4],
            "bunchingphase": d[:, :, 5], "energy": d[:, :, 6] / M_ELECTRON,
            "xsize": d[:, :, 8], "ysize": d[:, :, 9],
        }
    return {
        "z": d[:, 0], "power": d[:, 2], "bunching": d[:, 4],
        "bunchingphase": d[:, 5], "energy": d[:, 6] / M_ELECTRON,
        "xsize": d[:, 8], "ysize": d[:, 9],
    }


def load_fld(fn, nslice=1):
    """
    Field dump as (nslice, ny, nx) complex V/m, slices in time-window order (both codes
    unrotate on write). Either format, read through fieldio, so the tracker's openPMD dump
    and Genesis's own dump arrive in the same units and the same axis order.
    """
    return fieldio.read_field(fn, nslice=nslice)["u"]


def load_par(fn, nslice=1, window=None):
    """
    Particle dump with all slices concatenated, slice-major (beam slices never rotate).
    Either format, read through beamio.

    window is the {wavelength, spacing} the run used, needed for an openPMD dump and
    ignored for a Genesis one, which states its own. In a tier comparison it comes from the
    Genesis dump on the other side.
    """
    window = window or {}
    slices = beamio.read_slices(fn, window.get("wavelength"), window.get("spacing"))
    if len(slices) != nslice:
        raise ValueError(f"{fn}: {len(slices)} slices, expected {nslice}")
    return {k: np.concatenate([sl[k] for sl in slices])
            for k in ("gamma", "theta", "x", "y", "px", "py")}


def compare_tier(name, fortran_diag, genesis_out, fortran_fld, genesis_fld,
                 fortran_par, genesis_par, tolerance, nslice=1):
    """
    Compare one tier. Returns (worst_relative_difference, ok). nslice > 1 compares the
    full (nrecords, nslice) per-slice curves, all field slices and all particle slices.
    """
    print(f"--- {name} " + "-" * (74 - len(name)))

    f = load_fortran_diag(fortran_diag, nslice)
    g = load_genesis_out(genesis_out, nslice)

    n = min(len(f["z"]), len(g["z"]))
    if len(f["z"]) != len(g["z"]):
        print(f"  FAIL: record counts differ: fortran {len(f['z'])}, genesis {len(g['z'])}")
        return np.inf, False
    if np.abs(f["z"][:n] - g["z"][:n]).max() > 1e-9:
        print(f"  FAIL: record positions differ, max |dz| = "
              f"{np.abs(f['z'][:n]-g['z'][:n]).max():.3e} m")
        return np.inf, False

    worst = 0.0

    # Curves: elementwise relative for power (it spans decades and every point matters),
    # peak normalized for bunching (near zero at the quiet start, where elementwise
    # relative measures nothing but the quiet loading noise floor). In time-dependent
    # runs the power comparison is per slice per record. The denominator is floored at
    # 1e-9 of the global peak so a freshly slipped-in vacuum slice cannot divide by an
    # exact zero (inactive in steady state).
    den = np.maximum(np.abs(g["power"][:n]), 1e-9 * np.abs(g["power"][:n]).max())
    rel_power = (np.abs(f["power"][:n] - g["power"][:n]) / den).max()
    rel_bunch = np.abs(f["bunching"][:n] - g["bunching"][:n]).max() / np.abs(g["bunching"][:n]).max()
    worst = max(worst, rel_power, rel_bunch)
    slices = f" x {nslice} slices" if nslice > 1 else ""
    print(f"  power curve     ({n} records{slices})   elementwise max rel = {rel_power:.3e}")
    print(f"  bunching curve                  peak normalized     = {rel_bunch:.3e}")

    # Final field dump, peak normalized, all slices.
    uf, ug = load_fld(fortran_fld, nslice), load_fld(genesis_fld, nslice)
    rel_fld = np.abs(uf - ug).max() / np.abs(ug).max()
    worst = max(worst, rel_fld)
    print(f"  final field                     peak normalized     = {rel_fld:.3e}")

    # Final particle dump, particle by particle, all slices. gamma relative to itself,
    # transverse peak normalized. These are checked. The window comes from the Genesis dump,
    # since an openPMD beam file states the slice partition but not the wavelength it was
    # sliced on: that belongs to the run.
    window = beamio.genesis_window(genesis_par)
    pf, pg = load_par(fortran_par, nslice, window), load_par(genesis_par, nslice)
    scales = {"gamma": np.abs(pg["gamma"]).max(),
              "x": np.abs(pg["x"]).max(), "y": np.abs(pg["y"]).max(),
              "px": np.abs(pg["px"]).max(), "py": np.abs(pg["py"]).max()}
    parts = []
    for k in ("gamma", "x", "y", "px", "py"):
        r = np.abs(pf[k] - pg[k]).max() / scales[k]
        worst = max(worst, r)
        parts.append(f"{k} {r:.1e}")
    print(f"  final particles                 per-particle        = " + ", ".join(parts))

    # theta is reported but NOT checked. It is the phase of a particle in its ponderomotive
    # bucket, and at saturation neighboring trajectories separate exponentially, so the
    # worst-particle theta difference measures the Lyapunov amplification of whatever
    # difference exists, not the size of that difference (chaotic growth is not a usable
    # comparison metric). The distribution tells the story: the median is
    # the typical particle, the max is the separatrix tail. theta's collective effect IS
    # checked, through the bunching curve and the final field above.
    #
    # These are absolute phases on both sides. Neither dump format has a place for the
    # run's reference phase, so each code folds it into what it writes: Genesis stores
    # theta itself, and an openPMD beam file stores the lag that theta implies.
    dth = pf["theta"] - pg["theta"]
    print(f"  final theta (not checked; see comment)  max {np.abs(dth).max():.1e}, "
          f"rms {dth.std():.1e}, median {np.median(np.abs(dth)):.1e} rad")

    ok = worst <= tolerance
    print(f"  LARGEST RELATIVE DIFFERENCE: {worst:.6e}   (tolerance {tolerance:.1e})  "
          + ("ok" if ok else "FAIL"))
    print()
    return worst, ok


def compare_split(name, diag_a, diag_b, fld_a, fld_b, tolerance):
    """
    Fortran vs Fortran: split-weight run against the unsplit run.
    """
    print(f"--- {name} " + "-" * (74 - len(name)))
    a, b = load_fortran_diag(diag_a), load_fortran_diag(diag_b)
    worst = 0.0
    n = min(len(a["z"]), len(b["z"]))
    rel_power = (np.abs(a["power"][:n] - b["power"][:n]) / np.abs(b["power"][:n])).max()
    rel_bunch = np.abs(a["bunching"][:n] - b["bunching"][:n]).max() / np.abs(b["bunching"][:n]).max()
    ua, ub = load_fld(fld_a), load_fld(fld_b)
    rel_fld = np.abs(ua - ub).max() / np.abs(ub).max()
    worst = max(rel_power, rel_bunch, rel_fld)
    print(f"  power curve                     elementwise max rel = {rel_power:.3e}")
    print(f"  bunching curve                  peak normalized     = {rel_bunch:.3e}")
    print(f"  final field                     peak normalized     = {rel_fld:.3e}")
    ok = worst <= tolerance
    print(f"  LARGEST RELATIVE DIFFERENCE: {worst:.6e}   (tolerance {tolerance:.1e})  "
          + ("ok" if ok else "FAIL"))
    print()
    return worst, ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("workdir", help="Directory holding all run outputs")
    p.add_argument("--tol-tier1", type=float, default=1.0e-4)
    # The unaveraged tier is a priced model difference, not a transcription check: sin^2
    # end ramps vs Genesis's hard edges, no period averaging, RK4-vs-RK4 of different
    # equations. The self-referenced physics checks (check_unaveraged.py) pin the mode's
    # correctness. This tier pins its distance to Genesis so drift is visible.
    p.add_argument("--tol-tier1-unavg", type=float, default=1.5e-1)  # measured 6.9e-2
    p.add_argument("--tol-tier2-genesis", type=float, default=1.0e-3)
    p.add_argument("--tol-tier2-bmad", type=float, default=1.0e-1)
    p.add_argument("--tol-split", type=float, default=1.0e-10)
    p.add_argument("--tol-td1", type=float, default=1.0e-4)
    p.add_argument("--tol-td2-genesis", type=float, default=1.0e-3)
    p.add_argument("--tol-td2-bmad", type=float, default=1.0e-1)
    # Pure SASE: dark start, so the whole curve is startup-from-noise. Same transcribed
    # interlude model and constants floor as td2_genesis.
    p.add_argument("--tol-tdsase", type=float, default=1.0e-3)
    # Collective tiers: the sc floor is Genesis's truncated epsilon_0 in longRange
    # (8.85e-12, 4.7e-4 relative, measured 2.4e-4). The wake tier measured 8.7e-7.
    p.add_argument("--tol-tdsc", type=float, default=1.0e-3)
    p.add_argument("--tol-tdwk", type=float, default=1.0e-4)
    p.add_argument("--results", help="Append tier levels here for doc generation")
    args = p.parse_args()
    w = args.workdir

    print("=" * 78)
    print("FEL benchmark: Bmad tracker against Genesis 1.3 Version 4")
    print("=" * 78)
    print()

    nslice_td = load_nslice(f"{w}/AramisTD-initial.par.h5")
    nslice_tdsase = load_nslice(f"{w}/AramisTDSASE-initial.par.h5")

    results = []
    for name, diag, out, ffld, gfld, fpar, gpar, tol, nsl in (
        ("tier1: FEL core, one undulator segment",
         f"{w}/tier1.diag.txt", f"{w}/Aramis1seg.out.h5",
         f"{w}/tier1-final.wf.h5", f"{w}/Aramis1seg-final.fld.h5",
         f"{w}/tier1-final.beam.h5", f"{w}/Aramis1seg-final.par.h5",
         args.tol_tier1, 1),
        ("tier1_unavg: one segment, unaveraged dynamics vs Genesis4 (priced model difference)",
         f"{w}/tier1u.diag.txt", f"{w}/Aramis1seg.out.h5",
         f"{w}/tier1u-final.wf.h5", f"{w}/Aramis1seg-final.fld.h5",
         f"{w}/tier1u-final.beam.h5", f"{w}/Aramis1seg-final.par.h5",
         args.tol_tier1_unavg, 1),
        ("tier2_genesis: full line, transcribed interludes",
         f"{w}/tier2g.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2g-final.wf.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2g-final.beam.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_genesis, 1),
        ("tier2_bmad: full line, Bmad seam interludes",
         f"{w}/tier2.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2-final.wf.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2-final.beam.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_bmad, 1),
        ("td1: FEL core + slippage, one undulator segment",
         f"{w}/td1.diag.txt", f"{w}/AramisTD1seg.out.h5",
         f"{w}/td1-final.wf.h5", f"{w}/AramisTD1seg-final.fld.h5",
         f"{w}/td1-final.beam.h5", f"{w}/AramisTD1seg-final.par.h5",
         args.tol_td1, nslice_td),
        ("td2_genesis: full line time dependent, transcribed interludes",
         f"{w}/td2g.diag.txt", f"{w}/AramisTD.out.h5",
         f"{w}/td2g-final.wf.h5", f"{w}/AramisTD-final.fld.h5",
         f"{w}/td2g-final.beam.h5", f"{w}/AramisTD-final.par.h5",
         args.tol_td2_genesis, nslice_td),
        ("td2_bmad: full line time dependent, Bmad seam interludes",
         f"{w}/td2.diag.txt", f"{w}/AramisTD.out.h5",
         f"{w}/td2-final.wf.h5", f"{w}/AramisTD-final.fld.h5",
         f"{w}/td2-final.beam.h5", f"{w}/AramisTD-final.par.h5",
         args.tol_td2_bmad, nslice_td),
        ("tdsase: full line, pure SASE (dark start, growth from shot noise alone)",
         f"{w}/tdsase.diag.txt", f"{w}/AramisTDSASE.out.h5",
         f"{w}/tdsase-final.wf.h5", f"{w}/AramisTDSASE-final.fld.h5",
         f"{w}/tdsase-final.beam.h5", f"{w}/AramisTDSASE-final.par.h5",
         args.tol_tdsase, nslice_tdsase),
        ("tdsc: one segment TD, space charge on (short range + long range)",
         f"{w}/tdsc.diag.txt", f"{w}/AramisTDSC.out.h5",
         f"{w}/tdsc-final.wf.h5", f"{w}/AramisTDSC-final.fld.h5",
         f"{w}/tdsc-final.beam.h5", f"{w}/AramisTDSC-final.par.h5",
         args.tol_tdsc, nslice_td),
        ("tdwk: one segment TD, all wake kernels on",
         f"{w}/tdwk.diag.txt", f"{w}/AramisTDWK.out.h5",
         f"{w}/tdwk-final.wf.h5", f"{w}/AramisTDWK-final.fld.h5",
         f"{w}/tdwk-final.beam.h5", f"{w}/AramisTDWK-final.par.h5",
         args.tol_tdwk, nslice_td),
    ):
        worst, ok = compare_tier(name, diag, out, ffld, gfld, fpar, gpar, tol, nsl)
        results.append((name, worst, ok))

    worst, ok = compare_split(
        "weight_split: nonuniform weights must be invisible",
        f"{w}/tier1s.diag.txt", f"{w}/tier1.diag.txt",
        f"{w}/tier1s-final.wf.h5", f"{w}/tier1-final.wf.h5",
        args.tol_split)
    results.append(("weight_split: nonuniform weights must be invisible", worst, ok))

    print("=" * 78)
    print("Summary")
    print("=" * 78)
    all_ok = True
    for name, worst, ok in results:
        print(f"  {'pass' if ok else 'FAIL'}  {worst:.6e}  {name}")
        all_ok = all_ok and ok
    print()
    print("The tier2_bmad number is a measured transport model difference, not an error:")
    print("Genesis samples the quad path-length term at mid-element, Bmad integrates it")
    print("exactly. tier2_genesis proves it: same code, Genesis's interlude model, and")
    print("the difference collapses by six orders of magnitude. See doc/validation.md.")
    print()
    # The results file feeds doc generation. Levels are written at the same six
    # significant figures the summary prints, so the generated table and the summary
    # cannot disagree.
    if args.results:
        with open(args.results, "a") as fh:
            for name, worst, ok in results:
                key, _, desc = name.partition(":")
                fh.write(f"tier|{key}|{worst:.6e}|{'pass' if ok else 'fail'}|{desc.strip()}\n")

    print("PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
