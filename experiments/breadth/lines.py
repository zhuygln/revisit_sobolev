"""Importable line-list construction for the breadth sweep.

`sweep.py` builds its windows, mixes and per-species line tables at import
time and then launches 72 SEDONA runs, so nothing downstream could reuse its
line lists without re-running the sweep. This module carries exactly that
construction -- same data, same Boltzmann populations, same window selection
-- with no side effects, so the recompute scripts and the deterministic
reference can rebuild any condition's analytic legs in seconds.

Two additions over sweep.py: `stim=True` folds the LTE stimulated-emission
factor into each line's population (SEDONA applies it; see
optical_depth.stimulated_emission_factor -- it is 5e-3 at 9100 A), and the
window list is fixed to what the sweep actually used rather than re-scanned.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.atomic_data import load_gsi
from sobolev.constants import C, SIGMA_CLASSICAL
from sobolev.optical_depth import stimulated_emission_factor
from sobolev.populations import boltzmann_fractions_from_levels, statistical_weight

M_P = 1.67262192e-24
T_SHELL = 3000.0
T_CORE = 6000.0
V_CORE, V_MAX = 1.0e8, 3.0e8
V_D = 1.0e7
EPOCHS_DAY = [0.5, 1.0, 3.0]
TAU_REF = 5.0
WINDOWS = [(4300.0, 4400.0), (4900.0, 5000.0), (7000.0, 7100.0), (9100.0, 9200.0)]

SPECIES = {
    "LaII": (57, 139, "57LaII_levels_calib.txt", "57LaII_transitions_calib.txt"),
    "CeII": (58, 140, "58CeII_levels_calib.txt", "58CeII_transitions_calib.txt"),
    "CeIII": (59, 140, "58CeIII_levels_calib.txt", "58CeIII_transitions_calib.txt"),
}
MIXES = {
    "La": ["LaII"],
    "LaCe": ["LaII", "CeII"],
    "LaCeCe3": ["LaII", "CeII", "CeIII"],
}

_DATA = {}


def data(name):
    if name not in _DATA:
        z, a, lev_f, tr_f = SPECIES[name]
        _DATA[name] = dict(
            z=z, a=a,
            levels=load_gsi(ROOT / "data" / lev_f),
            lines=load_gsi(ROOT / "data" / tr_f),
        )
    return _DATA[name]


def window_lines(name, lo, hi, t_exp, stim=False):
    """Window-selected lines with Boltzmann populations and tau per unit rho.
    Identical to sweep.py's, plus the optional stimulated-emission factor."""
    d = data(name)
    lam = d["lines"]["WV_Transition"].to_numpy()
    win = d["lines"][(lam >= lo) & (lam < hi)].reset_index(drop=True)
    if len(win) == 0:
        return None
    frac = boltzmann_fractions_from_levels(d["levels"], T_SHELL)
    pop = frac[win["Lower"].to_numpy()]
    if stim:
        pop = pop * stimulated_emission_factor(C / (win["WV_Transition"].to_numpy() * 1e-8), T_SHELL)
    g_l = statistical_weight(win["J_Lower"].to_numpy())
    f_lu = 10 ** win["Log(gf)"].to_numpy() / g_l
    tau_per_rho = (
        SIGMA_CLASSICAL * f_lu * pop * win["WV_Transition"].to_numpy() * 1e-8
        * t_exp / (d["a"] * M_P)
    )
    return dict(win=win, pop=pop, f_lu=f_lu, tau_per_rho=tau_per_rho, **d)


def condition(mix_name, lo, t_day, stim=False):
    """Everything needed to rebuild one breadth condition's analytic legs:
    (lines_an, rho, r_core, r_out, t_exp, tau_max, n_lines, band, red_margin,
    nu_lo, nu_hi) with the sweep's exact geometry and normalization."""
    hi = lo + 100.0
    t_exp = t_day * 86400.0
    r_core, r_out = V_CORE * t_exp, V_MAX * t_exp
    members = MIXES[mix_name]
    x_frac = 1.0 / len(members)
    specs = [w for w in (window_lines(nm, lo, hi, t_exp, stim) for nm in members) if w is not None]
    ref = window_lines("LaII", 3850.0, 3950.0, 86400.0, stim=False)  # sweep.py's rho definition
    rho = TAU_REF / ref["tau_per_rho"].max() * (1.0 / t_day) ** 3
    tau_max = float(np.concatenate([s["tau_per_rho"] * rho * x_frac for s in specs]).max())
    n_lines = int(sum(len(s["win"]) for s in specs))
    lines_an = [
        (C / (lam_k * 1e-8), f_k, p_k * x_frac / (s["a"] * M_P))
        for s in specs
        for lam_k, f_k, p_k in zip(s["win"]["WV_Transition"], s["f_lu"], s["pop"])
    ]
    lam_blue, lam_red_edge = lo * 0.985, hi * 1.010
    return dict(
        lines=lines_an, rho=rho, r_core=r_core, r_out=r_out, t_exp=t_exp,
        tau_max=tau_max, n_lines=n_lines,
        band=(lo * 0.9885, hi + 1.0),
        # keep the margin inside the grid and clear of the final partial bin
        red_margin=(hi * 1.002, hi * 1.006),
        nu_lo=C / (lam_red_edge * 1e-8), nu_hi=C / (lam_blue * 1e-8),
        tag=f"{mix_name}_w{int(lo)}_t{t_day:g}",
    )


def all_conditions(stim=False):
    for mix in MIXES:
        for lo, _ in WINDOWS:
            for t_day in EPOCHS_DAY:
                yield condition(mix, lo, t_day, stim)
