#!/usr/bin/env python3
"""
The bunching the shot-noise load imposes on a beamlet, in closed form (docs/shot-noise.md).

A beamlet is M macroparticles at phases evenly spaced over one radiation period, standing
for N real electrons. The load (fel_fawley_noise) kicks every phase by

    theta_j -= sum_k a_k sin(k theta_j + phi_k),   a_k = (2/k) sqrt(x_k / N),   k = 1..H

with phi_k uniform, x_k exponential with unit mean and H = (M-1)/2 in integer division.
The default beamlet is M = 8, so H = 3. Averaged over the
draws, the beamlet's bunching at harmonic h is

    <|b_h|^2> = (1/M) sum_{j=0}^{M-1} cos(2 pi j h / M)
                      exp( -(2 h^2 / N) sum_k (1 - cos(2 pi j k / M)) / k^2 )

exactly. N electrons at random phases have <|b_h|^2> = 1/N, so the ratio of the two,
N <|b_h|^2>, is what the load delivers of the shot-noise power. It is 1 - 2/N to first
order for one harmonic, and 1 as N grows.

Below one electron fel_fawley_noise draws the kick for one electron. The beamlet's bunching is
then the bunching at N = 1, and the ratio is N times the ratio at 1.

This module is arithmetic only. report_shot_noise.py writes the page's table from it,
shot_noise_figures.py draws from it, and check_shot_noise.py holds loaded beams to it.
"""

from __future__ import annotations

import numpy as np

__all__ = ["ratio", "ratio_one_harmonic", "slice_ratio", "bunching_sq"]


def bunching_sq(n, h, m=8, nharm=None):
    """<|b_h|^2> of one beamlet of m phases kicked for n electrons, n at least 1 here.

    nharm is the number of harmonics kicked together, (m-1)/2 as fel_fawley_noise takes it
    when not given. n may be an array."""
    if nharm is None:
        nharm = (m - 1) // 2
    n = np.asarray(n, dtype=float)
    j = np.arange(m)
    k = np.arange(1, nharm + 1)
    spread = ((1 - np.cos(2 * np.pi * np.outer(j, k) / m)) / k**2).sum(axis=1)
    terms = np.cos(2 * np.pi * j * h / m) * np.exp(-2 * h**2 * spread / n[..., None])
    return terms.sum(axis=-1) / m


def ratio(n, h, m=8, nharm=None):
    """What the load delivers of the shot-noise power 1/n at harmonic h, for a beamlet of
    n electrons. Below one electron the kick is the one for a single electron."""
    n = np.asarray(n, dtype=float)
    return np.minimum(n, 1) * np.maximum(n, 1) * bunching_sq(np.maximum(n, 1), h, m, nharm)


def ratio_one_harmonic(n):
    """The ratio with one harmonic kicked and the phases dense, n exp(-2/n) I_1(2/n).

    I_1 is taken from its integral over one period, which the trapezoid rule sums to
    rounding for a periodic integrand."""
    n = np.asarray(n, dtype=float)
    t = 2 * np.pi * np.arange(4096) / 4096
    y = 2 / np.maximum(n, 1)[..., None]
    scaled = (np.exp(y * (np.cos(t) - 1)) * np.cos(t)).mean(axis=-1)
    return np.minimum(n, 1) * np.maximum(n, 1) * scaled


def slice_ratio(electrons, h, m=8, nharm=None):
    """The expected <|b_h|^2> N_lambda of a slice whose beamlets hold the given electron
    counts: the charge-weighted mean of the beamlets' ratios."""
    e = np.asarray(electrons, dtype=float)
    return float(np.sum(e * ratio(e, h, m, nharm)) / np.sum(e))
