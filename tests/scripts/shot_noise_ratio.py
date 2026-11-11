#!/usr/bin/env python3
"""
The bunching the shot-noise load imposes on a beamlet, in closed form (docs/shot-noise.md).

A beamlet is M macroparticles at phases evenly spaced over one radiation period, standing
for N real electrons. The load (fel_fawley_noise) kicks every phase by

    theta_j -= sum_k (z_k / k) sin(k theta_j + phi_k),   k = 1..H

with phi_k uniform and H = (M-1)/2 in integer division. The default beamlet is M = 8, so
H = 3. Each harmonic draws x_k, exponential with unit mean, and asks for the bunching
sqrt(x_k / N). A kick of size z at its own harmonic gives the beamlet the bunching J1(z).
Two kicks are in the program, beam.shot_noise_method:

    "bessel"       J1(z_k) = sqrt(x_k / N), and z_k = Z_MAX where sqrt(x_k / N) is above
                   B_MAX, the maximum of J1
    "first_order"  z_k = 2 sqrt(x_k / N), the small-kick form, with N set to 1 below one
                   electron

Averaged over the phases phi_k, the beamlet's bunching at harmonic h is

    <|b_h|^2> = (1/M) sum_{j=0}^{M-1} cos(2 pi j h / M) prod_k g(2 (h/k) |sin(pi j k / M)|)

exactly, with g(q) the mean of J0(q z) over the draw of x. For the first-order kick
g(q) = exp(-q^2 / N). For the kick from J1 it is one integral over z, taken here by
Gauss-Legendre quadrature as fel_shot_ratio takes it.

N electrons at random phases have <|b_h|^2> = 1/N, so the ratio of the two, N <|b_h|^2>,
is what the load delivers of the shot-noise power.

This module is arithmetic only. report_shot_noise.py writes the page's table from it,
shot_noise_figures.py draws from it, and check_shot_noise.py holds loaded beams to it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import j0, j1

__all__ = ["Z_MAX", "B_MAX", "bessel_kick", "ratio", "ratio_one_harmonic", "slice_ratio", "bunching_sq",
           "ratio_weight_total", "ratio_weight_kept"]

# The kick at the first maximum of J1, and the bunching there: the most a kick gives.
Z_MAX = 1.8411837813406593
B_MAX = 0.5818652242815963

# What is left out of the integral over the draw where the largest kick is not reached.
X_CUT = 40.0
GL_X, GL_W = np.polynomial.legendre.leggauss(16)
PANELS = 4


def _j1_slope(z):
    """J1'(z) = J0(z) - J1(z)/z, 1/2 at zero."""
    z = np.asarray(z, dtype=float)
    safe = np.where(z == 0, 1.0, z)
    return np.where(z == 0, 0.5, j0(z) - j1(z) / safe)


def bessel_kick(bunching):
    """The kick z in [0, Z_MAX] with J1(z) = bunching, Z_MAX where bunching is B_MAX or more.

    Newton's method from the third-order inverse, with a bisection wherever a step leaves
    the bracket, as fel_bessel_kick. Each element stops when its iterate no longer moves."""
    t = np.clip(np.asarray(bunching, dtype=float), 0.0, B_MAX)
    shape = t.shape
    t = t.ravel()
    lo = np.zeros_like(t)
    hi = np.full_like(t, Z_MAX)
    z = np.minimum(t * (2 + t * t), Z_MAX)
    live = np.arange(t.size)
    for _ in range(200):
        zl, tl = z[live], t[live]
        f = j1(zl) - tl
        hi[live] = np.where(f > 0, zl, hi[live])
        lo[live] = np.where(f <= 0, zl, lo[live])
        slope = _j1_slope(zl)
        step = zl - f / np.where(slope > 0, slope, 1.0)
        bad = (slope <= 0) | (step < lo[live]) | (step > hi[live])
        step = np.where(bad, 0.5 * (lo[live] + hi[live]), step)
        z[live] = step
        live = live[step != zl]
        if not live.size:
            break
    return np.where(t >= B_MAX, Z_MAX, z).reshape(shape)


def _bessel_nodes(n):
    """The kick sizes z and their weights for the mean over the draw at n electrons, and
    the probability that the draw is held at Z_MAX."""
    if B_MAX**2 * n <= X_CUT:
        top, tail = Z_MAX, np.exp(-B_MAX**2 * n)
    else:
        top, tail = float(bessel_kick(np.sqrt(X_CUT / n))), 0.0
    edges = np.linspace(0.0, top, PANELS + 1)
    half = 0.5 * (edges[1:] - edges[:-1])[:, None]
    mid = 0.5 * (edges[1:] + edges[:-1])[:, None]
    z = (mid + half * GL_X).ravel()
    w = (half * GL_W).ravel() * 2 * n * j1(z) * _j1_slope(z) * np.exp(-n * j1(z) ** 2)
    return z, w, tail


def bunching_sq(n, h, m=8, nharm=None, kick="bessel"):
    """<|b_h|^2> of one beamlet of m phases kicked for n electrons.

    nharm is the number of harmonics kicked together, (m-1)/2 as fel_fawley_noise takes it
    when not given. n may be an array. Under the first-order kick n is at least 1 here.
    A beamlet of more than 64 phases is taken as 64 under the kick from J1, as
    fel_shot_ratio takes it."""
    if nharm is None:
        nharm = (m - 1) // 2
    n = np.asarray(n, dtype=float)
    k = np.arange(1, nharm + 1)
    if kick == "first_order":
        j = np.arange(m)
        spread = ((1 - np.cos(2 * np.pi * np.outer(j, k) / m)) / k**2).sum(axis=1)
        terms = np.cos(2 * np.pi * j * h / m) * np.exp(-2 * h**2 * spread / n[..., None])
        return terms.sum(axis=-1) / m
    if kick != "bessel":
        raise ValueError(f"kick is bessel or first_order, not {kick!r}")
    m = min(m, 64)
    j = np.arange(m)
    q = 2 * (h / k) * np.abs(np.sin(np.pi * np.outer(j, k) / m))
    out = np.empty(n.shape)
    for idx in np.ndindex(n.shape):
        z, w, tail = _bessel_nodes(float(n[idx]))
        g = (w * j0(q[..., None] * z)).sum(axis=-1) + tail * j0(q * Z_MAX)
        out[idx] = (np.cos(2 * np.pi * j * h / m) * g.prod(axis=1)).sum() / m
    return out


def ratio(n, h, m=8, nharm=None, kick="bessel"):
    """What the load delivers of the shot-noise power 1/n at harmonic h, for a beamlet of
    n electrons. Below one electron the first-order kick is the one for a single electron."""
    n = np.asarray(n, dtype=float)
    if kick == "first_order":
        return np.minimum(n, 1) * np.maximum(n, 1) * bunching_sq(np.maximum(n, 1), h, m, nharm, kick)
    return n * bunching_sq(n, h, m, nharm, kick)


def ratio_one_harmonic(n, kick="bessel"):
    """The ratio with one harmonic kicked and the phases dense.

    The first-order kick gives n exp(-2/n) I_1(2/n). I_1 is taken from its integral over
    one period, which the trapezoid rule sums to rounding for a periodic integrand. The
    kick from J1 gives the mean of min(x, B_MAX^2 n) over the draw, 1 - exp(-B_MAX^2 n)."""
    n = np.asarray(n, dtype=float)
    if kick == "bessel":
        return -np.expm1(-B_MAX**2 * n)
    t = 2 * np.pi * np.arange(4096) / 4096
    y = 2 / np.maximum(n, 1)[..., None]
    scaled = (np.exp(y * (np.cos(t) - 1)) * np.cos(t)).mean(axis=-1)
    return np.minimum(n, 1) * np.maximum(n, 1) * scaled


def slice_ratio(electrons, h, m=8, nharm=None, kick="bessel"):
    """The expected <|b_h|^2> N_lambda of a slice whose beamlets hold the given electron
    counts: the charge-weighted mean of the beamlets' ratios."""
    e = np.asarray(electrons, dtype=float)
    counts, inverse = np.unique(e, return_inverse=True)
    r = ratio(counts, h, m, nharm, kick)[inverse]
    return float(np.sum(e * r) / np.sum(e))


def ratio_weight_total(n_total):
    """What the weight load delivers of the shot-noise power in the sample mode, at every
    harmonic below the beamlet size: the bunch's n_total electrons are counted as the
    nearest whole number, at least one, each of the bunch charge over that number."""
    n_total = np.asarray(n_total, dtype=float)
    return n_total / np.maximum(1.0, np.rint(n_total))


def ratio_weight_kept(n):
    """What the weight load delivers of the shot-noise power in the keep mode, for a particle
    of n electrons: n times the mean of one over its count. The count is floor(n) or one
    more, with the probability of floor(n) that makes that mean 1/n, so the ratio is 1.
    A particle under one electron is counted as one, and the ratio is n."""
    return np.minimum(np.asarray(n, dtype=float), 1.0)
