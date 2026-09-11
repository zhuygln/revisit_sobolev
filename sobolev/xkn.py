"""The xkn-diff secular-ejecta prescription (Ricigliano et al. 2024, MNRAS 529,
647; arXiv:2311.15709 v1 equation numbers) as the boundary condition of the
Paper IV P1 light curve (Phase 10c): the optically thick photospheric
luminosity, the photosphere, the temperature continuation into the thin
layers, and the thin layers' local heating. Every published choice is a
default here; every convention this repository takes differently is named.

THICK LUMINOSITY (their Sec. 2.1). The two-moment comoving-frame transport
of a uniform sphere rho = rho0 (t0/t)^3 (eq. 7-8, rho0 = M / (4/3 pi (v t0)^3))
with Eddington closure, separable E(x, t) = E0 (t0/t)^4 psi(x) phi(t) (eq. 9,
E0 = a T0^4), eigenfunctions psi_n = sqrt(2) sin(n pi x) / x (eq. 18), the
diffusion time tau0 = 3 kappa rho0 (v t0)^2 / c (eq. 12) and the mode ODE
(eq. 22)

    phi_n' + (t / (t0 tau0)) n^2 pi^2 phi_n = (-1)^(n+1) rho0 sqrt2 / (n pi E0) (t / t0) eps(t) f_th(t)

with phi_n(t0) = delta_n1; the surface flux gives (eq. 23)

    L(t) = 4 pi (v t0)^3 sqrt2 E0 / tau0 * sum_n (-1)^(n+1) n pi phi_n(t).

xkn's thick-core thermalisation is f_th = f_th0 (t/t0)^-beta (Sec. 2.1.2).
Here the modes are integrated by an exact integrating-factor recurrence on
a log grid (the kernel exp(-(a(t) - a(t'))) <= 1 never overflows; the
paper's closed form uses the upper incomplete gamma of a negative argument
and takes the real part). The result is rescaled by M_thick / M_ej (eq. 24)
with M_thick the mass inside the photosphere of the (1 - x^2)^3 profile
(eq. 25).

PHOTOSPHERE (eq. 26-27): tau_gamma(x) = kappa R rho0(t) int_x^1 (1 - x'^2)^3 dx'
= 2/3, R = v_max t (the R factor the printed eq. 27 drops and the code
carries), solved exactly here (xkn approximates R_ph by a parabola, 8 %).
It reaches the centre at t2 = sqrt(27 M kappa / (8 pi v_max^2)) (eq. 28).

TEMPERATURES: T_ph = max[(L_thick / (4 pi sigma R_ph^2))^(1/4), T_floor]
(eq. 29); thin layers T_i = T_ph (1 - x_i^2) / (1 - x_ph^2) (eq. 50), floored.

THIN LAYERS (eq. 47, 59): each layer radiates its local deposited heating
f_th,i eps dM_i with Barnes et al. (2016) thermalisation evaluated at
X = t (1 - x^2)^-1 (eq. 59); L = L_thick + L_thin (eq. 48). In the P1 light
curve the thin layers are TRANSPORTED instead (the transport replacement
test): their heating is injected as packets, L_thin here is only the
published quantity our transport replaces.

HEATING: xkn's default is the Perego et al. (2022) SkyNet fit
eps = eps_1d (t / 1 d)^-alpha (eq. 57); its Korobkin option ("K") is
eps = 1.95e10 (t / 1 d)^-1.3 erg s^-1 g^-1. The Korobkin law is the default
here (the SkyNet coefficients are not in this repository); `heating` may
be any callable of t (s).

CONVENTIONS TAKEN DIFFERENTLY (stated in the record): xkn's code uses
v_max = 3 v_rms for the (1 - x^2)^3 profile under a one-dimensional mass
measure and sqrt(5/3) v_rms in the diffusion solution; this repository's
P1 profile uses the spherically consistent v_max = 1.915 v_rms
(`ejecta.xkn_vmax_from_vrms`) and the spherical mass integral, and keeps
sqrt(5/3) v_rms in the diffusion solution as xkn does.

Published secular-component values (their Sec. 5.2, Table 2): M = 2.64e-2
Msun, v_rms = 0.06 c, Y_e = 0.2, s = 10 k_B, tau_exp = 17 ms; the fixed
opacity kappa(Y_e = 0.2) = 22.3 cm^2 g^-1 (Tanaka et al. 2020 parametrisation,
xkn `kappa_2_ye.py`), T_floor = 984 K in the fixed-opacity magnitude fit;
the shipped diffusion constants t0 = 3597 s, T0 = 4.17e4 K, f_th0 = 0.866,
beta = 0.243 (f_th = 0.70 at 0.1 d, 0.40 at 1 d).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from .constants import C
from .source import MSUN, DAY, SIGMA_SB, barnes_params
from .ejecta import xkn_vmax_from_vrms

A_RAD = 4.0 * SIGMA_SB / C


def heating_korobkin_xkn(t, eps_1d=1.95e10, alpha=1.3):
    """xkn's Korobkin option: eps = eps_1d (t / 1 d)^-alpha, erg s^-1 g^-1."""
    return eps_1d * (np.asarray(t, float) / DAY) ** (-alpha)


def _int_profile(x):
    """int_x^1 (1 - u^2)^3 du, closed form."""
    x = np.asarray(x, float)
    F = lambda u: u - u ** 3 + 3.0 * u ** 5 / 5.0 - u ** 7 / 7.0
    return F(1.0) - F(x)


def _int_mass(x):
    """int_x^1 u^2 (1 - u^2)^3 du, closed form (the spherical mass measure)."""
    x = np.asarray(x, float)
    G = lambda u: u ** 3 / 3.0 - 3.0 * u ** 5 / 5.0 + 3.0 * u ** 7 / 7.0 - u ** 9 / 9.0
    return G(1.0) - G(x)


def mass_outside(x):
    """Fraction of the (1 - x^2)^3 profile's mass outside x (spherical measure)."""
    return _int_mass(x) / _int_mass(0.0)


def photosphere_x(t, m, v_max, kappa, tau_ph=2.0 / 3.0):
    """x_ph = v_ph / v_max of the (1 - x^2)^3 profile at epoch t (s): the exact
    root of kappa R rho0 int_x^1 (1 - u^2)^3 du = tau_ph with R = v_max t and
    rho0 = 315 M / (64 pi R^3) (the spherical normalisation). 0 once the
    profile is thin to the centre (t >= t2, eq. 28)."""
    t = float(t); R = v_max * t
    rho0 = 315.0 * m / (64.0 * np.pi * R ** 3)
    tau_c = kappa * R * rho0 * _int_profile(0.0)
    if tau_c <= tau_ph:
        return 0.0
    f = lambda x: kappa * R * rho0 * _int_profile(x) - tau_ph
    return float(brentq(f, 0.0, 1.0, xtol=1e-12))


def t_photosphere_centre(m, v_max, kappa):
    """eq. 28: the epoch at which the photosphere reaches the centre."""
    return float(np.sqrt(27.0 * m * kappa / (8.0 * np.pi * v_max ** 2)))


def thick_luminosity_modes(t, m, v_diff, kappa, T0, t0, heating, f_th0, beta, n_modes=1000, n_grid=4000):
    """L_diff(t) of eq. 23 for the uniform-sphere diffusion solution (before
    the M_thick/M rescaling), on the requested epochs t (s, any order, all
    >= t0). The modes phi_n are integrated on a log grid from t0 to max(t)
    by exponential-Euler steps (exact for a constant rate and source per
    step; the quasi-steady limit phi = S / r for the fast high modes), then
    interpolated to t. Returns L_diff (erg s^-1)."""
    t = np.atleast_1d(np.asarray(t, float))
    rho0 = m / (4.0 / 3.0 * np.pi * (v_diff * t0) ** 3)
    tau0 = 3.0 * kappa * rho0 * (v_diff * t0) ** 2 / C
    E0 = A_RAD * T0 ** 4
    tg = np.geomspace(t0, max(t.max(), t0 * (1.0 + 1e-9)), n_grid)
    n = np.arange(1, n_modes + 1, dtype=float)
    sign = (-1.0) ** (n + 1)
    src = (tg / t0) * np.asarray(heating(tg), float) * f_th0 * (tg / t0) ** (-beta)   # (n_grid,)
    S = (sign * rho0 * np.sqrt(2.0) / (n * np.pi * E0))[:, None] * src[None, :]
    phi = np.zeros((n_modes, n_grid))
    phi[0, 0] = 1.0                                                                  # phi_n(t0) = delta_n1
    # exponential-Euler steps: with the decay rate r_n(t) = n^2 pi^2 t / (t0 tau0)
    # and the source both taken at the interval midpoint, the step is exact for
    # constant r and S, reaches the quasi-steady limit phi = S / r for r dt >> 1
    # (the high modes at late times, where a trapezoid on the kernel fails), and
    # phi + S dt for r dt << 1
    rate = (n[:, None] ** 2) * np.pi ** 2 / (t0 * tau0)                               # per unit t
    for k in range(n_grid - 1):
        dt = tg[k + 1] - tg[k]
        tm = 0.5 * (tg[k] + tg[k + 1])
        r = rate[:, 0] * tm
        Sm = 0.5 * (S[:, k] + S[:, k + 1])
        decay = np.exp(-r * dt)
        phi[:, k + 1] = phi[:, k] * decay + Sm * (-np.expm1(-r * dt)) / r
    L_grid = 4.0 * np.pi * (v_diff * t0) ** 3 * np.sqrt(2.0) * E0 / tau0 * np.sum((sign * n * np.pi)[:, None] * phi, axis=0)
    return np.interp(t, tg, L_grid)


def f_th_thin(t, x, m_msun, v_c):
    """Barnes et al. (2016) eq. 35 thermalisation of a thin layer at x = v/v_max,
    evaluated at X = t (1 - x^2)^-1 (xkn eq. 59); t in s."""
    (a, b, d), _ = barnes_params(m_msun, v_c)
    X = np.asarray(t, float) / DAY / (1.0 - np.asarray(x, float) ** 2)
    y = 2.0 * b * X ** d
    return 0.36 * (np.exp(-a * X) + np.log1p(y) / y)


@dataclass
class XknSecular:
    """The published xkn secular component with the fixed Tanaka opacity, as
    the P1 light curve's boundary condition (defaults = the published values;
    see the module docstring for the conventions)."""
    m_msun: float = 2.64e-2
    v_rms_c: float = 0.06
    kappa: float = 22.3
    T_floor: float = 984.0
    t0: float = 3597.0
    T0: float = 4.17e4
    f_th0: float = 0.866
    beta: float = 0.243
    n_modes: int = 1000                 # 500 in xkn; 1000 is within 0.7 % of 4000 here
    heating: object = field(default=heating_korobkin_xkn)
    heating_name: str = "korobkin_xkn"
    v_max_profile_c: float = None          # default: the spherical rms convention of this repository
    v_max_diff_c: float = None             # default: sqrt(5/3) v_rms as xkn's diffusion solution

    def __post_init__(self):
        if self.v_max_profile_c is None:
            self.v_max_profile_c = float(xkn_vmax_from_vrms(self.v_rms_c * C) / C)
        if self.v_max_diff_c is None:
            self.v_max_diff_c = float(np.sqrt(5.0 / 3.0) * self.v_rms_c)

    @property
    def m(self):
        return self.m_msun * MSUN

    def x_ph(self, t):
        return photosphere_x(t, self.m, self.v_max_profile_c * C, self.kappa)

    def r_ph(self, t):
        return self.x_ph(t) * self.v_max_profile_c * C * float(t)

    def L_diff(self, t):
        return thick_luminosity_modes(t, self.m, self.v_max_diff_c * C, self.kappa, self.T0, self.t0, self.heating,
                                      self.f_th0, self.beta, self.n_modes)

    def L_thick(self, t):
        """eq. 23 rescaled by M_thick / M_ej (eq. 24) at each epoch."""
        t = np.atleast_1d(np.asarray(t, float))
        frac = np.array([1.0 - mass_outside(self.x_ph(ti)) for ti in t])
        return self.L_diff(t) * frac

    def T_ph(self, t):
        """eq. 29 (the floor applied; the photosphere is not re-solved at the floor)."""
        t = np.atleast_1d(np.asarray(t, float))
        L = self.L_thick(t); R = np.array([self.r_ph(ti) for ti in t])
        with np.errstate(divide="ignore", invalid="ignore"):
            T = np.where(R > 0, (L / (4.0 * np.pi * SIGMA_SB * np.maximum(R, 1e-30) ** 2)) ** 0.25, 0.0)
        return np.maximum(T, self.T_floor)

    def thin_temperatures(self, t, x):
        """eq. 50 for layers at x (>= x_ph), floored."""
        T_ph = float(self.T_ph(t)[0]); x_ph = self.x_ph(t)
        T = T_ph * (1.0 - np.asarray(x, float) ** 2) / max(1.0 - x_ph ** 2, 1e-12)
        return np.maximum(T, self.T_floor)

    def thin_heating_rate(self, t, x):
        """Local deposited heating per unit mass of a layer at x: eps(t) f_th(t, x)."""
        return np.asarray(self.heating(t), float) * f_th_thin(t, x, self.m_msun, self.v_max_profile_c)

    def L_thin(self, t, x_edges):
        """eq. 47 on layers bounded by x_edges (those outside x_ph), for the
        record: the quantity the transport replaces."""
        x_edges = np.asarray(x_edges, float); xm = 0.5 * (x_edges[1:] + x_edges[:-1])
        dm = self.m * (mass_outside(x_edges[:-1]) - mass_outside(x_edges[1:]))
        x_ph = self.x_ph(t)
        out = xm >= x_ph
        return float(np.sum(dm[out] * self.thin_heating_rate(t, xm[out])))
