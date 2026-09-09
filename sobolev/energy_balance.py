"""Normalisation of energy-conserving runs (Paper IV), without touching the
frozen `sobolev/photometry.py`.

Paper III normalised every leg with `photometry._scale(core="conserving")`:
escaped synthetic ergs = the core's window luminosity, a grey rescale that
hid both what returned to the core and what photon-number branching had
left in the gas. With indivisible energy packets (`run_mc(packets="energy")`)
nothing is left in the gas by radiative events, so the only thing to account
for is the energy that returns to the core. For an opaque core that
re-emits what it receives with the launch spectrum, every relaunch is an
independent draw from the launch distribution, and the multi-pass emergent
spectrum is the single-pass one times the geometric series 1/(1 - f_return).
That is exactly `photometry._scale(core="equilibrium")`: E_inj - E_core -
E_abs synthetic ergs are the window luminosity. It is exact, not grey, and
this module names it so.

`run_mc(core="reemit")` performs the relaunches explicitly; it exists to
validate the series (agreement per band < 1 %), not for production.
"""
import numpy as np

from . import photometry
from .energy_packets import energy_accounting


def is_energy_conserving(res, tol=0.0):
    """True when the run carried energy packets and no energy was left in
    the gas: comoving deposit zero (to tol) and no absorbed packets."""
    a = res["accounting"]
    return (res.get("packets") == "energy" and abs(a["E_dep_cm"]) <= tol * a["E_inj"]
            and a["E_abs"] == 0.0)


def scale_exact(res, l_core_window):
    """Synthetic-erg to erg/s for an energy-conserving run: the geometric
    series over core relaunches, `photometry._scale(core="equilibrium")`.
    Raises if the run left energy in the gas (thermal deposits, trapped
    packets), because then the series is not the whole story."""
    if not is_energy_conserving(res):
        a = res["accounting"]
        raise ValueError("scale_exact needs an energy-conserving run: "
                         f"packets={res.get('packets')!r}, E_dep_cm/E_inj={a['E_dep_cm'] / a['E_inj']:.3e}, "
                         f"E_abs/E_inj={a['E_abs'] / a['E_inj']:.3e}")
    return photometry._scale(res, l_core_window, "equilibrium")


def emergent_lnu_exact(res, edges, l_core_window):
    """Absolute emergent L_nu (erg s^-1 Hz^-1) of an energy-conserving run."""
    scale_exact(res, l_core_window)          # the check
    return photometry.emergent_lnu(res, edges, l_core_window, core="equilibrium")


def bolometric_exact(res, l_core_window):
    return res["accounting"]["E_esc"] * scale_exact(res, l_core_window)


def renorm_ratio(res):
    """Paper III's grey "conserving" scale over the exact one, minus 1.

    For a photon-packet run this is the renormalisation the harness needed
    (the level-energy deposits plus the work). For an energy-conserving run
    the deposits are zero and what is left is exactly the Doppler work
    W / E_esc: a physical adiabatic loss of the radiation field, O(v/c),
    which must NOT be re-radiated. It is reported per leg so the two cases
    can be told apart; the Phase 2 acceptance is that it equals W / E_esc."""
    return energy_accounting(res)["renorm_ratio"]


def band_luminosities(res, edges, l_core_window, core):
    """L_nu integrated per bin of `edges` under a photometry core
    convention ("absorbing" | "equilibrium" | "conserving"), erg/s."""
    l_nu = photometry.emergent_lnu(res, edges, l_core_window, core=core)
    return l_nu * np.diff(edges)
