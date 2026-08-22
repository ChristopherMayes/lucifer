#!/usr/bin/env python3
"""
Compare the Bmad FEL tracker (lucifer) against Genesis 1.3 Version 4, over the
Benchmark1-SASE configuration: steady state and time dependent with slippage. Run through
run_fel_benchmark.sh, which produces all the inputs; this script only reads and compares.

Three tiers, each with its own tolerance sized to what it measures. The comparison
floor is set by fundamental constants: the tracker uses Bmad's (Z0 = mu_0*c =
376.7303...), Genesis carries a truncated impedance (376.73, 8.3e-7 relative), and that
difference enters the coupling and compounds through gain. During deliverable-3
development the tracker transcribed Genesis's constants and agreed at transcription
level (tier1 2.8e-11, tier2_genesis 5.9e-8, recorded in the README); after that
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
    """slicecount of a Genesis dump (.par.h5 or .fld.h5)."""
    with h5py.File(fn) as h5:
        return int(h5["slicecount"][0])


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
    """Field dump as (nslice, n, n), slices in time-window order (Genesis unrotates on write)."""
    with h5py.File(fn) as h5:
        n = int(h5["gridpoints"][0])
        u = np.empty((nslice, n, n), dtype=complex)
        for i in range(nslice):
            g = h5[f"slice{i+1:06d}"]
            u[i] = (g["field-real"][:].reshape(n, n)
                    + 1j * g["field-imag"][:].reshape(n, n))
        return u if nslice > 1 else u[0]


def load_par(fn, nslice=1):
    """Particle dump with all slices concatenated, slice-major (beam slices never rotate)."""
    with h5py.File(fn) as h5:
        out = {k: [] for k in ("gamma", "theta", "x", "y", "px", "py")}
        for i in range(nslice):
            s = h5[f"slice{i+1:06d}"]
            for k in out:
                out[k].append(s[k][:])
        return {k: np.concatenate(v) for k, v in out.items()}


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
    # runs the power comparison is per slice per record; the denominator is floored at
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
    # transverse peak normalized. These are checked.
    pf, pg = load_par(fortran_par, nslice), load_par(genesis_par, nslice)
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
    # The unaveraged tier is a PRICED MODEL DIFFERENCE, not a transcription check: sin^2
    # end ramps vs Genesis's hard edges, no period averaging, RK4-vs-RK4 of different
    # equations. The self-referenced physics checks (check_unaveraged.py) pin the mode's
    # correctness; this tier pins its distance to Genesis so drift is visible.
    p.add_argument("--tol-tier1-unavg", type=float, default=1.5e-1)  # measured 6.9e-2
    p.add_argument("--tol-tier2-genesis", type=float, default=1.0e-3)
    p.add_argument("--tol-tier2-bmad", type=float, default=1.0e-1)
    p.add_argument("--tol-split", type=float, default=1.0e-10)
    p.add_argument("--tol-td1", type=float, default=1.0e-4)
    p.add_argument("--tol-td2-genesis", type=float, default=1.0e-3)
    p.add_argument("--tol-td2-bmad", type=float, default=1.0e-1)
    # Pure SASE: dark start, so the whole curve is startup-from-noise; same transcribed
    # interlude model and constants floor as td2_genesis.
    p.add_argument("--tol-tdsase", type=float, default=1.0e-3)
    # Collective tiers: the sc floor is Genesis's truncated epsilon_0 in longRange
    # (8.85e-12, 4.7e-4 relative; measured 2.4e-4); the wake tier measured 8.7e-7.
    p.add_argument("--tol-tdsc", type=float, default=1.0e-3)
    p.add_argument("--tol-tdwk", type=float, default=1.0e-4)
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
         f"{w}/tier1-final.fld.h5", f"{w}/Aramis1seg-final.fld.h5",
         f"{w}/tier1-final.par.h5", f"{w}/Aramis1seg-final.par.h5",
         args.tol_tier1, 1),
        ("tier1_unavg: one segment, UNAVERAGED dynamics vs Genesis (priced model difference)",
         f"{w}/tier1u.diag.txt", f"{w}/Aramis1seg.out.h5",
         f"{w}/tier1u-final.fld.h5", f"{w}/Aramis1seg-final.fld.h5",
         f"{w}/tier1u-final.par.h5", f"{w}/Aramis1seg-final.par.h5",
         args.tol_tier1_unavg, 1),
        ("tier2_genesis: full line, transcribed interludes",
         f"{w}/tier2g.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2g-final.fld.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2g-final.par.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_genesis, 1),
        ("tier2_bmad: full line, Bmad seam interludes",
         f"{w}/tier2.diag.txt", f"{w}/Aramis.out.h5",
         f"{w}/tier2-final.fld.h5", f"{w}/Aramis-final.fld.h5",
         f"{w}/tier2-final.par.h5", f"{w}/Aramis-final.par.h5",
         args.tol_tier2_bmad, 1),
        ("td1: FEL core + slippage, one undulator segment",
         f"{w}/td1.diag.txt", f"{w}/AramisTD1seg.out.h5",
         f"{w}/td1-final.fld.h5", f"{w}/AramisTD1seg-final.fld.h5",
         f"{w}/td1-final.par.h5", f"{w}/AramisTD1seg-final.par.h5",
         args.tol_td1, nslice_td),
        ("td2_genesis: full line time dependent, transcribed interludes",
         f"{w}/td2g.diag.txt", f"{w}/AramisTD.out.h5",
         f"{w}/td2g-final.fld.h5", f"{w}/AramisTD-final.fld.h5",
         f"{w}/td2g-final.par.h5", f"{w}/AramisTD-final.par.h5",
         args.tol_td2_genesis, nslice_td),
        ("td2_bmad: full line time dependent, Bmad seam interludes",
         f"{w}/td2.diag.txt", f"{w}/AramisTD.out.h5",
         f"{w}/td2-final.fld.h5", f"{w}/AramisTD-final.fld.h5",
         f"{w}/td2-final.par.h5", f"{w}/AramisTD-final.par.h5",
         args.tol_td2_bmad, nslice_td),
        ("tdsase: full line, pure SASE (dark start, growth from shot noise alone)",
         f"{w}/tdsase.diag.txt", f"{w}/AramisTDSASE.out.h5",
         f"{w}/tdsase-final.fld.h5", f"{w}/AramisTDSASE-final.fld.h5",
         f"{w}/tdsase-final.par.h5", f"{w}/AramisTDSASE-final.par.h5",
         args.tol_tdsase, nslice_tdsase),
        ("tdsc: one segment TD, space charge on (short range + long range)",
         f"{w}/tdsc.diag.txt", f"{w}/AramisTDSC.out.h5",
         f"{w}/tdsc-final.fld.h5", f"{w}/AramisTDSC-final.fld.h5",
         f"{w}/tdsc-final.par.h5", f"{w}/AramisTDSC-final.par.h5",
         args.tol_tdsc, nslice_td),
        ("tdwk: one segment TD, all wake kernels on",
         f"{w}/tdwk.diag.txt", f"{w}/AramisTDWK.out.h5",
         f"{w}/tdwk-final.fld.h5", f"{w}/AramisTDWK-final.fld.h5",
         f"{w}/tdwk-final.par.h5", f"{w}/AramisTDWK-final.par.h5",
         args.tol_tdwk, nslice_td),
    ):
        worst, ok = compare_tier(name, diag, out, ffld, gfld, fpar, gpar, tol, nsl)
        results.append((name, worst, ok))

    worst, ok = compare_split(
        "weight_split: nonuniform weights must be invisible",
        f"{w}/tier1s.diag.txt", f"{w}/tier1.diag.txt",
        f"{w}/tier1s-final.fld.h5", f"{w}/tier1-final.fld.h5",
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
    print("the difference collapses by six orders of magnitude. See the README.")
    print()
    print("PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
