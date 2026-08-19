#!/usr/bin/env python3
"""
Reconstruct a Bmad bunch_params_struct-shaped dict from any (record, slice) of the
tracker's stats file <out_root>.stats.h5 -- the per-record beam datasets (centroid,
sigma, charge_live, n_particle_live) are SUFFICIENT statistics, and this module is
the proof by construction: at element ends its output must reproduce the values the
tracker stored from Bmad's own calc_bunch_params (the harness holds that).

The formulas are transcribed from Bmad's calc_emittances_and_twiss_from_sigma_matrix
(beam_utils.f90): projected x/y/z twiss with dispersion removed via sigma(i,6)
columns, and the normal-mode a/b/c emittances as the imaginary parts of the
eigenvalues of sigma.S (Wolski Eq. 32; numpy does the eigensystem). Full normal-mode
BETA functions need Wolski's N matrix and are Bmad's job, not this script's --
projected twiss plus normal-mode emittances is what the reconstruction covers.

Fixed Bmad units throughout, matching the file. norm_emit uses
f_emit = p0c*(1 + <pz>)/m_e c^2 with the file's stored p0c.

Also provides the wavefront-side pooling helper: pulse-level moments from per-slice
wavefront_params (energy-weighted mean of slice sigmas plus the variance of slice
centroids) -- pulse values are DERIVED, never stored; the file stays raw.

Usage as a module:
  from bunch_params_from_stats import bunch_params_at, pool_wavefront
Or from the command line, for a quick look:
  bunch_params_from_stats.py <stats.h5> <record> <slice>
"""

from __future__ import annotations

import sys

import h5py
import numpy as np

M_ELECTRON = 0.51099895069e6   # Bmad's m_electron [eV]


def _projected(sig, i, j, with_dispersion=True):
    """Bmad's projected_twiss_calc for the plane with position index i, momentum j
    (0-based); dispersion removed via the pz column (index 5). The z plane passes
    with_dispersion=False -- Bmad hands it zero dispersion columns by construction."""
    tw = {}
    s66 = sig[5, 5]
    d_i = sig[i, 5] if with_dispersion else 0.0
    d_j = sig[j, 5] if with_dispersion else 0.0
    if s66 != 0:
        tw["eta"] = d_i / s66
        tw["etap"] = d_j / s66
        x2 = sig[i, i] - d_i**2 / s66
        x_px = sig[i, j] - d_i * d_j / s66
        px2 = sig[j, j] - d_j**2 / s66
    else:
        tw["eta"] = 0.0
        tw["etap"] = 0.0
        x2, x_px, px2 = sig[i, i], sig[i, j], sig[j, j]
    tw["sigma"] = np.sqrt(max(0.0, x2))
    tw["sigma_p"] = np.sqrt(max(0.0, px2))
    emit = np.sqrt(max(0.0, x2 * px2 - x_px**2))
    tw["emit"] = emit
    if emit != 0:
        tw["alpha"] = -x_px / emit
        tw["beta"] = x2 / emit
        tw["gamma"] = px2 / emit
    else:
        tw["alpha"] = tw["beta"] = tw["gamma"] = 0.0
    return tw


def bunch_params_from_moments(centroid, sigma36, charge_live, n_live, p0c):
    """The reconstruction proper: moments -> bunch_params dict. sigma36 is the
    flattened 6x6 (symmetric, so the flattening order cannot mislead)."""
    sig = np.asarray(sigma36).reshape(6, 6)
    f_emit = p0c * (1.0 + centroid[5]) / M_ELECTRON

    bp = {"centroid": np.asarray(centroid), "sigma": sig,
          "charge_live": charge_live, "n_particle_live": int(n_live)}
    for name, (i, j) in (("x", (0, 1)), ("y", (2, 3)), ("z", (4, 5))):
        tw = _projected(sig, i, j, with_dispersion = (name != "z"))
        tw["norm_emit"] = f_emit * tw["emit"]
        bp[name] = tw

    # Normal modes: eigenvalues of sigma.S (Wolski Eq. 32). S column-swaps with signs.
    S = np.zeros((6, 6))
    for k in range(3):
        S[2 * k, 2 * k + 1] = 1.0
        S[2 * k + 1, 2 * k] = -1.0
    dim = 6
    cut = 1e-20 * np.abs(sig).max()
    if abs(sig[4, 4]) < cut or abs(sig[5, 5]) < cut:
        dim = 4
    ev = np.linalg.eigvals(sig[:dim, :dim] @ S[:dim, :dim])
    emits = np.sort(np.abs(ev.imag))[::-1]     # Pairs +-i emit; take distinct values.
    emits = emits[::2]                         # One per conjugate pair.
    # NOTE: mode LABELS here are magnitude-descending. Bmad's mat_eigen identifies
    # modes by eigenvector structure, so a/b/c may be permuted relative to Bmad;
    # the harness compares the normal-mode emittances as a sorted set.
    for k, name in enumerate(("a", "b", "c")):
        e = float(emits[k]) if k < len(emits) else 0.0
        bp[name] = {"emit": e, "norm_emit": f_emit * e}
    return bp


def bunch_params_at(h5, record, islice):
    """bunch_params dict from stats file h5 (path or open file) at (record, slice)."""
    close = False
    if not isinstance(h5, h5py.File):
        h5 = h5py.File(h5, "r")
        close = True
    try:
        b = h5["beam"]
        bp = bunch_params_from_moments(
            b["centroid"][record, islice], b["sigma"][record, islice],
            float(b["charge_live"][record, islice]),
            int(b["n_particle_live"][record, islice]), float(h5["p0c"][0]))
    finally:
        if close:
            h5.close()
    return bp


def pool_wavefront(centroids, sigmas, energies):
    """Pulse-level wavefront moments pooled from per-slice params: energy-weighted
    mean of slice sigmas plus the variance of slice centroids. centroids (n,4),
    sigmas (n,16) or (n,4,4), energies (n,). Returns (centroid(4), sigma(4,4),
    energy)."""
    c = np.asarray(centroids, float).reshape(-1, 4)
    s = np.asarray(sigmas, float).reshape(-1, 4, 4)
    w = np.asarray(energies, float)
    wtot = w.sum()
    if wtot <= 0:
        return np.zeros(4), np.zeros((4, 4)), 0.0
    cbar = (w[:, None] * c).sum(axis=0) / wtot
    dc = c - cbar
    pooled = (w[:, None, None] * (s + dc[:, :, None] * dc[:, None, :])).sum(axis=0) / wtot
    return cbar, pooled, wtot


def main():
    fn, rec, isl = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    bp = bunch_params_at(fn, rec, isl)
    print(f"centroid: {bp['centroid']}")
    for m in ("x", "y", "z"):
        tw = bp[m]
        print(f"{m}: beta {tw['beta']:.6g}  alpha {tw['alpha']:.6g}  emit {tw['emit']:.6g}  "
              f"norm_emit {tw['norm_emit']:.6g}  sigma {tw['sigma']:.6g}  eta {tw['eta']:.6g}")
    for m in ("a", "b", "c"):
        print(f"{m}: emit {bp[m]['emit']:.6g}  norm_emit {bp[m]['norm_emit']:.6g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
