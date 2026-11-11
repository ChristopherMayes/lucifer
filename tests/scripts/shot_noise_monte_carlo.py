#!/usr/bin/env python3
"""
Monte Carlo of the shot-noise kick alone, against the closed form (docs/shot-noise.md).

No run of the program is taken. A beamlet of M evenly spaced phases at a random base
phase is kicked as fel_fawley_noise kicks it, with the draws of the kick the input file
names, and its bunching is summed over many beamlets. The mean of |b_h|^2 N at each kicked
harmonic is printed beside shot_noise_ratio.py's closed form, with the standard error of
the mean, and the share of the draws held at the largest bunching a kick gives beside
exp(-B_MAX^2 N).

With --weight the noise is imposed through the weights instead, as fel_weight_noise does,
and the mean of |sum_j c_j exp(-i h theta_j)|^2 over the electrons is printed at every
harmonic below M, in three forms: independent Poisson counts of mean N/M on a beamlet, the
counts of a slice of 16 beamlets drawn with their total held, and the keep mode's count of
a particle of N electrons, floor(N) or one more, spread over its copies.

With --joint the script instead asks whether one kick exists that gives a beamlet the
three bunchings its draws ask for at once, which the kick sized one harmonic at a time
does not give exactly. The three complex kicks are solved by Newton's method on the
beamlet's own phases, continued from zero kick along the straight line to the target.
The share of beamlets on which the solve reaches the target is printed.

The page quotes this script's numbers with their day. It runs on demand and no test
takes it.

Usage: shot_noise_monte_carlo.py --electrons 150 50 18 10 4 1 [--samples 4e8]
                                 [--kick bessel|first_order] [--seed 1]
       shot_noise_monte_carlo.py --weight --electrons 150 18 4 1 0.25 3.76 [--samples 1e8]
       shot_noise_monte_carlo.py --joint --electrons 150 50 18 10 [--samples 2000]
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


def kick_sizes(x, n, kick):
    """The kick's size at its own harmonic for the draws x, and how many were held."""
    if kick == "first_order":
        # The loader takes a kick above 2 pi modulo 2 pi, and sets N to 1 below one electron.
        size = 2 * np.sqrt(x / max(n, 1.0)) / HARMONICS[:, None]
        return np.where(size > 2 * np.pi, np.mod(size, 2 * np.pi), size) * HARMONICS[:, None], 0
    t = np.sqrt(x / n)
    return snr.bessel_kick(t), int(np.count_nonzero(t >= snr.B_MAX))


def monte_carlo(n, samples, kick, rng):
    total = np.zeros(H)
    square = np.zeros(H)
    held = 0
    done = 0
    base = 2 * np.pi * np.arange(M) / M
    while done < samples:
        theta = rng.uniform(0, 2 * np.pi, (CHUNK, 1)) + base
        phi = 2 * np.pi * rng.uniform(size=(H, CHUNK))
        z, n_held = kick_sizes(-np.log(rng.uniform(size=(H, CHUNK))), n, kick)
        held += n_held
        moved = theta.copy()
        for i, k in enumerate(HARMONICS):
            moved -= (z[i] / k)[:, None] * np.sin(k * theta + phi[i][:, None])
        for i, h in enumerate(HARMONICS):
            b2 = np.abs(np.exp(-1j * h * moved).mean(axis=1)) ** 2 * n
            total[i] += b2.sum()
            square[i] += (b2 * b2).sum()
        done += CHUNK
    mean = total / done
    return mean, np.sqrt((square / done - mean**2) / done), held / (H * done), done


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


def bunching(theta, a):
    """The beamlet's bunching at the kicked harmonics under the complex kicks a."""
    moved = theta - np.imag(a[None, :] * np.exp(1j * np.outer(theta, HARMONICS))).sum(axis=1)
    return np.exp(-1j * np.outer(HARMONICS, moved)).mean(axis=1), moved


def newton(theta, target, a):
    """Newton's method for bunching(theta, a) = target from the start a. Returns a, converged."""
    for _ in range(30):
        b, moved = bunching(theta, a)
        f = b - target
        if np.abs(f).max() < 1e-14:
            return a, True
        e = np.exp(-1j * np.outer(HARMONICS, moved)) * (1j * HARMONICS[:, None])
        d_re = (e[:, :, None] * np.sin(np.outer(theta, HARMONICS))[None, :, :]).mean(axis=1)
        d_im = (e[:, :, None] * np.cos(np.outer(theta, HARMONICS))[None, :, :]).mean(axis=1)
        jac = np.block([[d_re.real, d_im.real], [d_re.imag, d_im.imag]])
        try:
            step = np.linalg.solve(jac, -np.concatenate([f.real, f.imag]))
        except np.linalg.LinAlgError:
            return a, False
        if np.abs(step).max() > 1:
            return a, False
        a = a + step[:H] + 1j * step[H:]
    return a, np.abs(bunching(theta, a)[0] - target).max() < 1e-12


def joint_reach(n, samples, rng):
    """The share of beamlets on which the joint solve reaches the three bunchings drawn."""
    reached = 0
    for _ in range(samples):
        theta = rng.uniform(0, 2 * np.pi) + 2 * np.pi * np.arange(M) / M
        target = np.sqrt(-np.log(rng.uniform(size=H)) / n) * np.exp(2j * np.pi * rng.uniform(size=H))
        a = np.zeros(H, complex)
        s, ds = 0.0, 0.25
        while s < 1 and ds > 1e-3:
            s_try = min(1.0, s + ds)
            a_try, good = newton(theta, s_try * target, a)
            if good:
                a, s, ds = a_try, s_try, min(0.25, 1.5 * ds)
            else:
                ds /= 2
        reached += s == 1
    return reached / samples


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--electrons", type=float, nargs="+", required=True, help="Electrons a beamlet holds")
    p.add_argument("--samples", type=float, default=None, help="Beamlets at each count")
    p.add_argument("--kick", default="bessel", choices=("bessel", "first_order"))
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--joint", action="store_true", help="Measure the reach of the joint solve instead")
    p.add_argument("--weight", action="store_true", help="The noise through the weights instead of a kick")
    a = p.parse_args()
    for n in a.electrons:
        rng = np.random.default_rng([a.seed, int(round(1000 * n))])
        if a.weight:
            forms = {"poisson": 1.0,
                     "total held": float(snr.ratio_weight_total(16 * n)),
                     "kept": float(snr.ratio_weight_kept(n))}
            for key, (mean, err) in weight_monte_carlo(n, int(a.samples or 1e7), rng).items():
                # The keep mode under one electron has one copy at the whole charge, and no spread.
                off = np.abs(mean - forms[key]) / np.where(err > 0, err, 1.0)
                worst = np.argmax(off)
                print(f"N = {n:g}, weights, {key}: closed form {forms[key]:.5f}, Monte Carlo at harmonics 1 to {M - 1} "
                      + " ".join(f"{m:.5f}" for m in mean)
                      + f", standard error {err.max():.5f}, worst {off[worst]:.2f} "
                      f"standard errors", flush=True)
            continue
        if a.joint:
            samples = int(a.samples or 2000)
            print(f"N = {n:g}: the joint solve reaches its target on {joint_reach(n, samples, rng):.3f} "
                  f"of {samples} beamlets", flush=True)
            continue
        mean, err, held, done = monte_carlo(n, int(a.samples or 1e7), a.kick, rng)
        want = -np.expm1(-snr.B_MAX**2 * n) if a.kick == "bessel" else 1.0
        print(f"N = {n:g}, {a.kick}, {done:.1e} beamlets, draws held {held:.5f} "
              f"(closed form {1 - want:.5f})", flush=True)
        for i, h in enumerate(HARMONICS):
            form = float(snr.ratio(n, int(h), M, kick=a.kick))
            print(f"  h = {h}: Monte Carlo {mean[i]:.5f} +- {err[i]:.5f}, closed form {form:.5f}, "
                  f"difference {(mean[i] - form) / err[i]:+.2f} standard errors", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
