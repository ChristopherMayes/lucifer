#!/usr/bin/env python3
"""
Monte Carlo of the shot-noise load alone, against the closed forms (docs/shot-noise.md).

No run of the program is taken. A beamlet of M evenly spaced phases stands for N
electrons, and its noise is imposed as the loader imposes it.

By default the noise goes through the weights, as fel_weight_noise does. The mean of
|sum_j c_j exp(-i h theta_j)|^2 over the electrons is printed at every harmonic below M,
in three forms: independent Poisson counts of mean N/M on a beamlet, the counts of a
slice of 16 beamlets drawn with their total held, and the keep mode's count of a particle
of N electrons, floor(N) or one more, spread over its copies.

With --kick the beamlet is kicked by the first-order kick instead, as fel_fawley_noise
kicks it, at a random base phase, and the mean of |b_h|^2 N at each kicked harmonic is
printed beside shot_noise_ratio.py's closed form with the standard error of the mean.

The page quotes this script's numbers with their day. It runs on demand and no test
takes it.

Usage: shot_noise_monte_carlo.py --electrons 150 18 4 3.76 1 0.25 [--samples 1e8] [--seed 1]
       shot_noise_monte_carlo.py --kick --electrons 150 50 18 4 1 [--samples 1e8]
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import shot_noise_ratio as snr  # noqa: E402

M = 8
H = (M - 1) // 2
HARMONICS = np.arange(1, H + 1)
CHUNK = 1_000_000


def kick_monte_carlo(n, samples, rng):
    """The first-order kick on one beamlet: the mean of |b_h|^2 n at the kicked harmonics
    and its standard error."""
    total = np.zeros(H)
    square = np.zeros(H)
    done = 0
    base = 2 * np.pi * np.arange(M) / M
    while done < samples:
        theta = rng.uniform(0, 2 * np.pi, (CHUNK, 1)) + base
        phi = 2 * np.pi * rng.uniform(size=(H, CHUNK))
        # The loader sets n to 1 below one electron and takes a kick above 2 pi modulo 2 pi.
        size = 2 * np.sqrt(-np.log(rng.uniform(size=(H, CHUNK))) / max(n, 1.0)) / HARMONICS[:, None]
        size = np.where(size > 2 * np.pi, np.mod(size, 2 * np.pi), size)
        moved = theta.copy()
        for i, k in enumerate(HARMONICS):
            moved -= size[i][:, None] * np.sin(k * theta + phi[i][:, None])
        for i, h in enumerate(HARMONICS):
            b2 = np.abs(np.exp(-1j * h * moved).mean(axis=1)) ** 2 * n
            total[i] += b2.sum()
            square[i] += (b2 * b2).sum()
        done += CHUNK
    mean = total / done
    return mean, np.sqrt((square / done - mean**2) / done), done


def weight_monte_carlo(n, samples, rng):
    """The three weight forms at n electrons a beamlet. Returns {form: (mean, standard error)}
    of the noise power over the shot-noise power, one entry per harmonic 1 to M - 1."""
    harmonics = np.arange(1, M)
    phase = np.exp(-1j * np.outer(harmonics, 2 * np.pi * np.arange(M) / M))
    n_beamlet = 16
    sums = {k: [np.zeros(M - 1), np.zeros(M - 1), 0] for k in ("poisson", "total held", "kept")}

    def add(key, power, norm):
        sums[key][0] += (power / norm).sum(axis=1)
        sums[key][1] += ((power / norm) ** 2).sum(axis=1)
        sums[key][2] += power.shape[1]

    chunk = 200_000
    done = 0
    while done < samples:
        # Independent Poisson counts on one beamlet, against the mean count n.
        c = rng.poisson(n / M, (chunk, M))
        add("poisson", np.abs(phase @ c.T) ** 2, n)

        # A slice of n_beamlet beamlets at staggered base phases, its total held at the
        # nearest whole number of electrons, against the slice's mean count.
        total = max(1, int(np.rint(n_beamlet * n)))
        c = rng.multinomial(total, np.full(n_beamlet * M, 1 / (n_beamlet * M)), chunk // 8)
        theta = (2 * np.pi * np.arange(M)[None, :] / M + 2 * np.pi * (np.arange(n_beamlet)[:, None] + 0.5)
                 / (M * n_beamlet)).ravel()
        # Each counted electron has the slice's charge over the count, and the power goes as its square.
        add("total held", np.abs(np.exp(-1j * np.outer(harmonics, theta)) @ c.T) ** 2 * (n_beamlet * n / total) ** 2,
            n_beamlet * n)

        # The keep mode: the particle's charge held, n electrons counted as floor(n) or one more.
        low = np.floor(n)
        if low < 1:
            count = np.ones(chunk, dtype=np.int64)
        else:
            count = (low + (rng.random(chunk) >= low * (low + 1 - n) / n)).astype(np.int64)
        c = rng.multinomial(count, np.full(M, 1 / M))
        add("kept", np.abs(phase @ (c / count[:, None]).T) ** 2 * n, 1.0)
        done += chunk
    out = {}
    for key, (total, square, count) in sums.items():
        mean = total / count
        out[key] = (mean, np.sqrt((square / count - mean**2) / count))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--electrons", type=float, nargs="+", required=True, help="Electrons a beamlet holds")
    p.add_argument("--samples", type=float, default=1e7, help="Beamlets at each count")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--kick", action="store_true", help="The first-order kick instead of the weights")
    a = p.parse_args()
    for n in a.electrons:
        rng = np.random.default_rng([a.seed, int(round(1000 * n))])
        if a.kick:
            mean, err, done = kick_monte_carlo(n, int(a.samples), rng)
            print(f"N = {n:g}, first-order kick, {done:.1e} beamlets", flush=True)
            for i, h in enumerate(HARMONICS):
                form = float(snr.ratio(n, int(h), M))
                print(f"  h = {h}: Monte Carlo {mean[i]:.5f} +- {err[i]:.5f}, closed form {form:.5f}, "
                      f"difference {(mean[i] - form) / err[i]:+.2f} standard errors", flush=True)
            continue
        forms = {"poisson": 1.0,
                 "total held": float(snr.ratio_weight_total(16 * n)),
                 "kept": float(snr.ratio_weight_kept(n))}
        for key, (mean, err) in weight_monte_carlo(n, int(a.samples), rng).items():
            # The keep mode under one electron has one copy at the whole charge, and no spread.
            off = np.abs(mean - forms[key]) / np.where(err > 0, err, 1.0)
            worst = np.argmax(off)
            print(f"N = {n:g}, weights, {key}: closed form {forms[key]:.5f}, Monte Carlo at harmonics 1 to {M - 1} "
                  + " ".join(f"{m:.5f}" for m in mean)
                  + f", standard error {err.max():.5f}, worst {off[worst]:.2f} standard errors", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
