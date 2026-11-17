"""
The closed forms of undulator radiation that check_undulator_spectrum.py holds a run to
and report_undulator_spectrum.py draws (docs/undulator-spectrum.md).

One electron of Lorentz factor gamma crosses n_period periods of an ideal undulator of
period lambda_u. K is the deflection parameter, e B lambda_u / (2 pi m c): of the peak
field for a planar undulator, and of the field's constant magnitude for a helical one.
Far from the undulator the energy radiated per unit angular frequency and solid angle is
(J. D. Jackson, Classical Electrodynamics, 3rd ed., Eq. 14.67)

    d2W/(d omega d Omega) = e^2 omega^2 / (16 pi^3 eps0 c)
                            |integral of n x (n x beta) exp(i omega (t - n.r/c)) dt|^2

with n the direction of observation. Near the fundamental the integral over n_period
periods is the integral over one period times the sum of n_period phase factors, which
gives the line

    d2W/(d omega d Omega) = W1(theta) n_period^2 sinc^2(pi n_period (omega/omega_1(theta) - 1))

    omega_1(theta) = omega_1(0) / (1 + x),    x = gamma_z^2 theta^2,
    gamma_z^2 = gamma^2 / (1 + K^2/2) planar, gamma^2 / (1 + K^2) helical

and on the axis (K.-J. Kim, AIP Conf. Proc. 184, 565 (1989))

    W1(0) = e^2 gamma^2 / (4 pi eps0 c) K^2 [JJ]^2 / (1 + K^2/2)^2      planar
    W1(0) = e^2 gamma^2 / (4 pi eps0 c) 2 K^2 / (1 + K^2)^2             helical

    [JJ] = J0(xi) - J1(xi),    xi = K^2 / (4 + 2 K^2).

period_level() evaluates W1(theta) by the integral over one period at omega_1(theta), on
the ideal trajectory to second order in 1/gamma, and reduces to on_axis_level() at
theta = 0.
"""

from __future__ import annotations

import numpy as np

C_LIGHT = 299792458.0
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
M_ELECTRON = 0.51099895069e6          # [eV]
H_PLANCK = 6.62607015e-34


class Undulator:
    """An ideal undulator and the electron that crosses it."""

    def __init__(self, e_tot, b_max, l_period, n_period, helical):
        self.gamma = e_tot / M_ELECTRON
        self.l_period = l_period
        self.n_period = n_period
        self.helical = helical
        self.k = b_max * l_period * C_LIGHT / (2 * np.pi * M_ELECTRON)
        self.k2 = self.k ** 2 if helical else self.k ** 2 / 2          # a_w^2
        self.gamma_z2 = self.gamma ** 2 / (1 + self.k2)
        self.wavelength = l_period * (1 + self.k2) / (2 * self.gamma ** 2)

    def on_axis_level(self):
        """W1(0), the energy per unit angular frequency and solid angle at the peak of
        the fundamental on the axis, per electron and per period squared [J s / sr]."""
        a = E_CHARGE ** 2 * self.gamma ** 2 / (4 * np.pi * EPS0 * C_LIGHT)
        if self.helical:
            return a * 2 * self.k ** 2 / (1 + self.k ** 2) ** 2
        xi = self.k ** 2 / (4 + 2 * self.k ** 2)
        return a * self.k ** 2 * (_bessel_j(0, xi) - _bessel_j(1, xi)) ** 2 / (1 + self.k ** 2 / 2) ** 2

    def period_level(self, theta_x, theta_y, n_point=2000):
        """W1 at the angles (theta_x, theta_y) [rad], by the integral over one period
        at the resonant frequency of that direction [J s / sr]."""
        ku = 2 * np.pi / self.l_period
        beta_z = 1 - (1 + self.k2) / (2 * self.gamma ** 2)
        period = self.l_period / (beta_z * C_LIGHT)
        t = (np.arange(n_point) + 0.5) * period / n_point
        phase = ku * beta_z * C_LIGHT * t
        bx = self.k / self.gamma * np.cos(phase)
        x = self.k / (self.gamma * ku * beta_z) * np.sin(phase)
        if self.helical:
            by = self.k / self.gamma * np.sin(phase)
            y = -self.k / (self.gamma * ku * beta_z) * np.cos(phase)
            z = beta_z * C_LIGHT * t
        else:
            by = np.zeros_like(phase)
            y = np.zeros_like(phase)
            z = beta_z * C_LIGHT * t - self.k ** 2 / (8 * self.gamma ** 2 * ku) * np.sin(2 * phase)
        bz = np.sqrt(1 - 1 / self.gamma ** 2 - bx ** 2 - by ** 2)
        nx, ny = theta_x, theta_y
        nz = np.sqrt(1 - nx ** 2 - ny ** 2)
        # The frequency at which one period of the motion advances the phase by 2 pi.
        omega = 2 * np.pi / (period - nz * self.l_period / C_LIGHT)
        kernel = np.exp(1j * omega * (t - (nx * x + ny * y + nz * z) / C_LIGHT)) * period / n_point
        nb = nx * bx + ny * by + nz * bz
        amp = [np.sum((n * nb - b) * kernel) for n, b in ((nx, bx), (ny, by), (nz, bz))]
        return E_CHARGE ** 2 * omega ** 2 / (16 * np.pi ** 3 * EPS0 * C_LIGHT) * sum(abs(a) ** 2 for a in amp)

    def angular_factor(self, theta_x, theta_y):
        """W1 at the angles over W1 on the axis."""
        return self.period_level(theta_x, theta_y) / self.period_level(0.0, 0.0)


def line(n_period, nu):
    """sinc^2(pi n_period nu), the line of n_period periods at the relative frequency
    nu = omega/omega_1 - 1, of unit height."""
    return np.sinc(n_period * np.asarray(nu)) ** 2


def sampled_line(n_period, nu):
    """The line of a pulse of n_period equal samples one wavelength apart, of unit
    height: (sin(pi n_period nu) / (n_period sin(pi nu)))^2. It differs from line() by
    (pi nu)^2 / 3 of itself to lowest order."""
    nu = np.asarray(nu, dtype=float)
    den = n_period * np.sin(np.pi * nu)
    small = np.abs(den) < 1e-12
    return np.where(small, 1.0, np.sin(np.pi * n_period * nu) / np.where(small, 1.0, den)) ** 2


LINE_WIDTH = 0.885892941          # The full width at half maximum of sinc^2(pi u), in u.


def _bessel_j(order, x, n_point=2000):
    """J_order(x) by its integral, (1/pi) integral of cos(order t - x sin t) over 0 to pi."""
    t = (np.arange(n_point) + 0.5) * np.pi / n_point
    return float(np.mean(np.cos(order * t - x * np.sin(t))))
