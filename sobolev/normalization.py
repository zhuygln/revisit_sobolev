"""Two normalization standards, for two different questions.

Until F35 this project used one recipe everywhere: n_ion such that the strongest
Sobolev depth *inside 3850-3950 A* equals 5. That is fine for a single ion and
wrong for a cross-ion claim, because it is ion-specific by accident -- it works
for La/Ce/Nd only because their strongest lines happen to lie near that window.
Yb II needs n_ion = 1.7e12 cm^-3 to satisfy it, giving tau_max = 1.7e8 and
beta = 1.5e-8, at which point a packet entering that resonance never escapes.

The two questions need different standards and must not be mixed:

CONTROLLED (`global_tau_max`)
    Every ion is placed at the same GLOBAL maximum line strength. This asks:
    *at matched line strength, what does the atomic network's topology change?*
    It is scale-free, applies identically to every ion, and bounds
    beta >= (1 - e^-tau_max)/tau_max by construction. Use for any comparison
    ACROSS ions -- compression, memory, the phase diagram.

ASTROPHYSICAL (`from_conditions`)
    Every ion is placed at the same physical state -- T, rho, epoch, mass
    fractions and ion fraction -- and whatever optical-depth distribution
    physics gives it is what it gets. This asks: *in a real kilonova, what does
    each ion actually do?* Use for anything that claims an observable
    consequence.

A result quoted under one standard does not transfer to the other. §4.32:
Ce II's band-3800 grouped-closure error is +126.7% under the old window recipe
and +12.2% under `global_tau_max`, from a 25% density difference.
"""
import numpy as np

from .atomic_data import load_gsi
from .constants import C
from .optical_depth import SIGMA_CLASSICAL
from .populations import boltzmann_fractions_from_levels, statistical_weight

WINDOW_LEGACY = (3850.0, 3950.0)


def _tau_per_n(levels, lines, temperature, t_exp, window=None):
    """Sobolev depth per unit ion density, per line (optionally windowed)."""
    lam = lines["WV_Transition"].to_numpy()
    m = np.ones(lam.size, bool) if window is None else \
        (lam >= window[0]) & (lam < window[1])
    if not m.any():
        return None
    frac = boltzmann_fractions_from_levels(levels, temperature)
    pop = frac[lines["Lower"].to_numpy()[m]]
    g_l = statistical_weight(lines["J_Lower"].to_numpy()[m])
    f_lu = 10 ** lines["Log(gf)"].to_numpy()[m] / g_l
    return SIGMA_CLASSICAL * f_lu * pop * lam[m] * 1e-8 * t_exp


def global_tau_max(lev_path, tr_path, temperature, t_exp, tau_max=5.0):
    """CONTROLLED standard: n_ion putting the ion's strongest line at `tau_max`.

    Returns (n_ion, diagnostics). `beta_min` is bounded below by
    (1 - e^-tau_max)/tau_max, so the branch chain always terminates.
    """
    levels = load_gsi(lev_path); lines = load_gsi(tr_path)
    tpn = _tau_per_n(levels, lines, temperature, t_exp)
    mx = float(np.nanmax(tpn))
    if not np.isfinite(mx) or mx <= 0:
        return None, {"reason": "no line with positive optical depth"}
    n_ion = tau_max / mx
    return n_ion, {"standard": "global_tau_max", "tau_max_target": tau_max,
                   "beta_min_bound": float(-np.expm1(-tau_max) / tau_max),
                   "n_lines": int(tpn.size)}


def window_tau_max(lev_path, tr_path, temperature, t_exp, tau_max=5.0,
                   window=WINDOW_LEGACY):
    """LEGACY standard: `tau_max` inside a fixed wavelength window.

    Retained so results quoted under it can be reproduced exactly. It is
    ion-specific by accident and will diverge for ions whose window holds only
    weak lines -- `diverges` flags that, and the caller should refuse to run.
    """
    levels = load_gsi(lev_path); lines = load_gsi(tr_path)
    tpn_w = _tau_per_n(levels, lines, temperature, t_exp, window)
    if tpn_w is None or not np.isfinite(np.nanmax(tpn_w)) or np.nanmax(tpn_w) <= 0:
        return None, {"reason": f"no line in {window[0]:.0f}-{window[1]:.0f} A"}
    n_ion = tau_max / float(np.nanmax(tpn_w))
    tpn_all = _tau_per_n(levels, lines, temperature, t_exp)
    tau_max_global = float(np.nanmax(tpn_all)) * n_ion
    beta_min = float(-np.expm1(-tau_max_global) / tau_max_global)
    return n_ion, {"standard": "window_tau_max", "window": list(window),
                   "tau_max_global": tau_max_global, "beta_min": beta_min,
                   # a forest whose strongest line is this thick traps packets
                   "diverges": bool(tau_max_global > 1e3)}


def from_conditions(lev_path, tr_path, temperature, t_exp, rho, mass_frac,
                    ion_frac, atomic_mass):
    """ASTROPHYSICAL standard: n_ion from a physical state.

    n_ion = rho * X_element * f_ion / (A * m_u). No optical depth is targeted;
    whatever tau distribution the state produces is the answer.
    """
    M_U = 1.66053906660e-24
    n_ion = rho * mass_frac * ion_frac / (atomic_mass * M_U)
    return n_ion, {"standard": "from_conditions", "rho": rho,
                   "mass_frac": mass_frac, "ion_frac": ion_frac,
                   "atomic_mass": atomic_mass, "temperature": temperature,
                   "t_exp": t_exp}
