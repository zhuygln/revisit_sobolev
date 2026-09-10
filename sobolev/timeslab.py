"""Time-slab helpers for the time-dependent light curve (Paper IV Phase 10b).

`run_mc(..., relativity="worldline", t_stop=, resume=, launch_energy=,
max_events=)` transports a packet population from its state up to `t_stop`
and returns the paused packets with their state. This module holds the
physics-side pieces around it: the packet population as a checkpoint, the
Fontes et al. (2020) heating law, the trapped radiation field of a state as
an initial population, the escape records as observables (escape time,
observer time, frequency, energy), the radiation temperature of a
population (the feedback rule), and the absolute emergent L_nu of a slab.
Everything here is a pure function of arrays and states; the driver
(`paper4/phase10_fontes/lightcurve.py`) sequences them.
"""
from dataclasses import dataclass

import numpy as np

from .constants import C, H, K_B
from .source import SIGMA_SB, DAY

A_RAD = 4.0 * SIGMA_SB / C                 # radiation constant, erg cm^-3 K^-4


@dataclass
class Population:
    """A packet population at one epoch: the checkpoint between slabs."""
    r: np.ndarray
    mu: np.ndarray
    nu: np.ndarray
    w: np.ndarray
    ctime: np.ndarray
    shell_of: np.ndarray
    n_core_passes: np.ndarray

    @property
    def n(self):
        return int(self.r.size)

    def as_dict(self):
        return dict(r=self.r, mu=self.mu, nu=self.nu, w=self.w, ctime=self.ctime, shell_of=self.shell_of,
                    n_core_passes=self.n_core_passes)

    def energy_lab(self):
        """Sum of the packets' lab-frame energies w h nu."""
        return float(np.sum(self.w * H * self.nu))

    def energy_cm(self):
        """Sum of the comoving energies at the packets' own positions and
        epochs (homologous flow: beta = r / ctime)."""
        beta = self.r / self.ctime
        gam = 1.0 / np.sqrt(1.0 - beta * beta)
        return float(np.sum(self.w * H * self.nu * gam * (1.0 - beta * self.mu)))

    def to_npz(self, path):
        np.savez(path, **self.as_dict())

    @classmethod
    def from_npz(cls, path):
        d = np.load(path)
        return cls(**{k: d[k] for k in ("r", "mu", "nu", "w", "ctime", "shell_of", "n_core_passes")})

    @classmethod
    def empty(cls):
        z = np.zeros(0)
        return cls(z, z.copy(), z.copy(), z.copy(), z.copy(), np.zeros(0, np.int32), np.zeros(0, np.int32))


# ---------------------------------------------------------------------------
# heating
# ---------------------------------------------------------------------------

def fontes_heating_rate(t, eps0=8.2e8, t0=4.0 * DAY, index=-1.3):
    """Specific heating rate of the Fontes et al. (2020) simplified problem,
    erg s^-1 g^-1: eps0 (t/t0)^index with t0 = 4 d, eps0 = 8.2e8."""
    return eps0 * (np.asarray(t, float) / t0) ** index


def heating_energy(m_shell, t_a, t_b, eps0=8.2e8, t0=4.0 * DAY, index=-1.3):
    """Energy released in [t_a, t_b] by a power-law heating rate uniform in
    mass: (E_total, E_per_shell), closed form."""
    m_shell = np.asarray(m_shell, float)
    if index == -1.0:
        integral = eps0 * t0 * np.log(t_b / t_a)
    else:
        integral = eps0 * t0 / (index + 1.0) * ((t_b / t0) ** (index + 1.0) - (t_a / t0) ** (index + 1.0))
    per = m_shell * integral
    return float(per.sum()), per


# ---------------------------------------------------------------------------
# populations from states and from runs
# ---------------------------------------------------------------------------

def _planck_energy_draw(rng, nu_lo, nu_hi, n, T):
    """Frequencies from the energy-weighted Planck law on [nu_lo, nu_hi]
    (the draw `run_mc`'s launch uses)."""
    grid = np.geomspace(nu_lo, nu_hi, 20001)
    x = H * grid / (K_B * T)
    wgt = grid ** 3 / np.expm1(np.minimum(x, 700.0))
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (wgt[1:] + wgt[:-1]) * np.diff(grid))]); cdf /= cdf[-1]
    return np.interp(rng.uniform(0.0, 1.0, n), cdf, grid)


def initial_radiation(state, n, rng, nu_lo, nu_hi, shells, T=None):
    """The trapped radiation field of `state` as a packet population at
    `state.t`: energy a T^4 V per transported shell (T = T_rad unless given),
    packets allocated to shells in proportion to that energy, positions
    uniform in the shell's volume, directions isotropic in the comoving frame
    and aberrated to the lab, frequencies from the energy-weighted Planck law
    at the shell's T in the comoving frame, equal comoving energy per packet
    so that sum(w h nu_cm) == sum(a T^4 V) to roundoff."""
    shells = list(shells)
    T_s = np.asarray(state.T_rad if T is None else T, float)[shells]
    r_e = state.r_edges
    V = 4.0 * np.pi / 3.0 * np.array([r_e[s + 1] ** 3 - r_e[s] ** 3 for s in shells])
    E_s = A_RAD * T_s ** 4 * V
    E_tot = float(E_s.sum())
    if n <= 0 or E_tot <= 0:
        return Population.empty(), E_tot, E_s
    k = rng.choice(len(shells), size=n, p=E_s / E_tot)
    ct = C * float(state.t)
    r3lo = np.array([r_e[s] ** 3 for s in shells]); r3hi = np.array([r_e[s + 1] ** 3 for s in shells])
    r = np.cbrt(r3lo[k] + rng.uniform(0.0, 1.0, n) * (r3hi[k] - r3lo[k]))
    mu_c = rng.uniform(-1.0, 1.0, n)
    nu_cm = np.empty(n)
    for i in range(len(shells)):
        m = k == i
        if m.any():
            nu_cm[m] = _planck_energy_draw(rng, nu_lo, nu_hi, int(m.sum()), float(T_s[i]))
    beta = r / ct; gam = 1.0 / np.sqrt(1.0 - beta * beta); den = 1.0 + beta * mu_c
    mu = (mu_c + beta) / den
    nu = nu_cm * gam * den
    w = (E_tot / n) / (H * nu_cm)
    pop = Population(r, mu, nu, w, np.full(n, ct), k.astype(np.int32), np.zeros(n, np.int32))
    return pop, E_tot, E_s


def carried_population(res, r_edges_next=None):
    """The packets `run_mc` paused (fate 4) as a Population. `shell_of` is
    recomputed on `r_edges_next` (the next slab's grid at the pause epoch) by
    searchsorted, a packet exactly on an edge with mu < 0 belonging to the
    inner shell; r is clipped into the zone."""
    m = np.asarray(res["fate"]) == 4
    r = np.asarray(res["r"], float)[m]; mu = np.asarray(res["mu"], float)[m]
    nu = np.asarray(res["nu_final"], float)[m]; w = np.asarray(res["w"], float)[m]
    ctime = np.asarray(res["ctime"], float)[m]
    ncp = np.asarray(res["n_core_passes"], np.int32)[m]
    if r_edges_next is not None:
        e = np.asarray(r_edges_next, float)
        r = np.clip(r, e[0], e[-1] * (1.0 - 1e-15))
        sh = np.searchsorted(e, r, side="right") - 1
        on_edge = (np.searchsorted(e, r, side="left") == np.searchsorted(e, r, side="right") - 1)   # r == an edge
        sh = np.where(on_edge & (mu < 0.0), sh - 1, sh)
        sh = np.clip(sh, 0, e.size - 2).astype(np.int32)
    else:
        sh = np.asarray(res["shell_of"], np.int32)[m] if res.get("zoned") else np.zeros(int(m.sum()), np.int32)
    return Population(r, mu, nu, w, ctime, sh, ncp)


def escapes(res):
    """The escaped packets of a run as observables: t_esc (s, the light-time
    at the outer boundary / c), t_obs = t_esc - r_esc mu_esc / c (the
    arrival-time correction for a distant observer along the packet's
    direction; Fontes et al. 2020 did not apply it), lab frequency nu (Hz)
    and lab energy e (erg)."""
    m = np.asarray(res["fate"]) == 1
    ct_esc = np.asarray(res["ct_esc"], float)[m]
    mu_esc = np.asarray(res["mu_esc"], float)[m]
    nu = np.asarray(res["nu_out_all"], float)[m]
    e = np.asarray(res["w"], float)[m] * H * nu
    r_esc = res["b_out"] * ct_esc
    return dict(t_esc=ct_esc / C, t_obs=(ct_esc - r_esc * mu_esc) / C, nu=nu, e=e)


def t_rad_from_population(pop, state, shells, floor=300.0, n_min=50):
    """Rule (b) of Phase 10b: the radiation temperature of each transported
    shell from the packets present, T = (sum w h nu_cm / (a V))^(1/4) with
    the comoving energies at the packets' own positions; shells with fewer
    than `n_min` packets get NaN (the caller falls back to the prescribed T
    and counts them). Returns (T per transported shell, counts)."""
    shells = list(shells)
    r_e = np.asarray(state.r_edges, float)
    V = 4.0 * np.pi / 3.0 * np.array([r_e[s + 1] ** 3 - r_e[s] ** 3 for s in shells])
    beta = pop.r / pop.ctime; gam = 1.0 / np.sqrt(1.0 - beta * beta)
    e_cm = pop.w * H * pop.nu * gam * (1.0 - beta * pop.mu)
    E = np.bincount(pop.shell_of, weights=e_cm, minlength=len(shells))[:len(shells)]
    n = np.bincount(pop.shell_of, minlength=len(shells))[:len(shells)]
    T = np.where(n >= n_min, np.maximum((np.maximum(E, 0.0) / (A_RAD * V)) ** 0.25, floor), np.nan)
    return T, n


def lnu_absolute(nu, e, edges, dt):
    """Absolute emergent L_nu (erg s^-1 Hz^-1) of the escapes of one slab of
    duration dt: the energy histogram over `edges` divided by dt dnu."""
    edges = np.asarray(edges, float)
    h, _ = np.histogram(nu, bins=edges, weights=e)
    return h / (dt * np.diff(edges))


def bin_light_curve(t, e, t_edges):
    """Energy and counts per time bin -> (L per bin [erg/s], counts)."""
    t_edges = np.asarray(t_edges, float)
    h, _ = np.histogram(t, bins=t_edges, weights=e)
    n, _ = np.histogram(t, bins=t_edges)
    return h / np.diff(t_edges), n
