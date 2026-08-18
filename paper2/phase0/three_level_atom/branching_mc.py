"""A small Sobolev Monte Carlo with line branching -- Paper II's instrument.

WHY THIS EXISTS. Paper II asks whether the expansion-opacity bias measured in
Paper I survives radiative redistribution. That requires a code that knows
WHICH line absorbed a photon and can re-emit it in a different line. The
P2-0A audit found that the public SEDONA cannot: `opacity_epsilon` is one
global scalar, and a line interaction is either coherent scattering or a
redraw from the zone's thermal pool. Neither channel carries the identity of
the absorbing transition, so there is no branching to measure
(paper2/phase0/sedona_source_audit/NOTES.md).

TARDIS has the physics but not the other half of the comparison: it has no
expansion-opacity mode, so a TARDIS-vs-SEDONA measurement would be cross-code.
Paper I's standing rule -- earned three separate times, each after a
cross-code artifact was mistaken for physics -- is that same-code differentials
are the robust ones. This module exists so that Phase 1 can vary branching and
opacity treatment INSIDE ONE CODE, with TARDIS held in reserve as an
independent cross-check.

WHAT IT IS. A 1-D spherical, homologously expanding shell around an opaque
lightbulb core, packets propagated under the Sobolev point-interaction
approximation. Deliberately small enough to read line by line, in the same
spirit as the ~200-line formal solver that produced Paper I's most defensible
results.

PHYSICS AND CONVENTIONS, stated because every one of them is a place a
cross-code comparison can silently go wrong:

* Homologous flow, v = r / t_exp. The projection of the local velocity onto a
  ray direction is z / (c t_exp) with z = r . n_hat, which is LINEAR in path
  length. The comoving frequency of a packet is therefore

      nu_cm(s) = nu_lab [1 - (z0 + s) / (c t_exp)],

  strictly decreasing along any ray. This is the same relation, with the same
  first-order Doppler convention, that `sobolev.sobolev_leg.sobolev_attenuation`
  uses for its resonance planes -- deliberately, so the two can be compared
  directly (see `tests/test_branching_mc.py`).

* RELATIVITY (Finding F11). First-order Doppler only. The worldline correction
  is tau -> tau/gamma, i.e. O(beta^2); at the ~0.003-0.03 c relevant here that
  is 5e-6 to 5e-4, far below Monte Carlo noise at any practical packet count.
  Paper I's F11 matters because it distinguishes a 1/gamma law from a spurious
  O(beta) one; there is no O(beta) term to get wrong here.

* SOBOLEV POINT INTERACTION. A line is a delta-function resonance at the plane
  where nu_cm equals its rest frequency. The packet interacts there with
  probability 1 - exp(-tau_S) and is otherwise unaffected. No profile width,
  no overlap -- this is Sobolev proper, the same approximation Paper I
  measured as accurate to <0.5% at physical line widths.

* RESONANCE IS CONSUMED, exactly and not by fiat. Because nu_cm decreases
  monotonically, a packet re-emitted at a line's rest frequency can never
  reach that same resonance again: it sits at s = 0 of that plane and moves
  away from it. The code enforces this only through a strictly-positive
  distance tolerance; nothing special-cases the line just used.

* POPULATIONS ARE FROZEN. Level populations are inputs, not solved for. This
  is Paper I's convention (fixed LTE) carried forward, and it keeps the
  branching measurement free of an NLTE feedback loop. Consequently a line
  whose lower level is unpopulated has no opacity but can still be a
  DOWNWARD channel -- which is exactly what makes a fluorescent escape
  channel escape.

* THE CORE IS OPAQUE. A packet that re-enters r_core is absorbed there. Photon
  number is therefore conserved as
  N_launched = N_escaped + N_core_absorbed, not as N_in = N_out; there is no
  packet destruction anywhere else. `run_mc` returns both counts so the
  identity is checkable rather than assumed.
"""

import numpy as np

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from sobolev.constants import C
from sobolev.optical_depth import tau_sobolev


# --------------------------------------------------------------------------
# Atom
# --------------------------------------------------------------------------


class Atom:
    """Lines plus the branching table that Paper II is about.

    Parameters
    ----------
    lines : sequence of dicts with keys
        nu0     rest frequency (Hz)
        lower   index of the lower level
        upper   index of the upper level
        f_osc   absorption oscillator strength
        n_lower number density of the lower level (cm^-3), frozen
        A       spontaneous rate (s^-1); only RATIOS within an upper level
                matter, since branching probabilities are A_uj / sum_j A_uj

    A line with n_lower = 0 carries no opacity but is still available as a
    downward branch -- the fluorescent escape channel of the three-level atom.
    """

    def __init__(self, lines):
        self.lines = [dict(l) for l in lines]
        self.nu0 = np.array([l["nu0"] for l in self.lines])
        self.upper = np.array([l["upper"] for l in self.lines])

        # Branching table: for each upper level, the lines that depopulate it
        # and their cumulative probabilities.
        self.branches = {}
        for u in sorted(set(self.upper)):
            idx = np.flatnonzero(self.upper == u)
            a = np.array([self.lines[i]["A"] for i in idx], dtype=float)
            if a.sum() <= 0:
                raise ValueError(f"upper level {u} has no downward rate")
            self.branches[u] = (idx, np.cumsum(a / a.sum()))

    def tau(self, t_exp):
        """Sobolev optical depth of every line, from the shared prefactor.

        Uses `sobolev.optical_depth.tau_sobolev` rather than re-deriving it:
        Paper I's very first bug was a hand-typed sigma_classical that was 18%
        wrong, and the fix was to have exactly one definition.
        """
        return np.array([
            tau_sobolev(l["f_osc"], l["n_lower"], C / l["nu0"], t_exp)
            for l in self.lines
        ])

    def branch(self, upper, u_rand):
        """Index of the line the packet leaves `upper` by, given a uniform."""
        idx, cum = self.branches[upper]
        return int(idx[np.searchsorted(cum, u_rand)])


def three_level_atom(nu_13, nu_32, f_osc, n_ground, a31, a32):
    """The P2-0B atom: ground 1, metastable sink 2, upper 3.

    One opacity line (1->3). Two downward channels from level 3: back to
    ground at nu_13 (resonant scattering) and to level 2 at nu_32 < nu_13
    (fluorescence). Level 2 is a sink -- no population, hence no 2->1
    opacity and no way back -- so a fluorescent packet leaves the shell
    without further interaction. That makes the fluorescent yield directly
    countable in the emergent spectrum, which is the whole point of the test.
    """
    if not nu_32 < nu_13:
        raise ValueError("fluorescent line must be redder than the pump")
    return Atom([
        {"nu0": nu_13, "lower": 1, "upper": 3, "f_osc": f_osc,
         "n_lower": n_ground, "A": a31},
        {"nu0": nu_32, "lower": 2, "upper": 3, "f_osc": 0.0,
         "n_lower": 0.0, "A": a32},
    ])


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


def _distance_to_boundary(r, mu, r_core, r_out):
    """(distance, 'out'|'core') to the next shell boundary along the ray."""
    perp2 = r * r * (1.0 - mu * mu)
    s_out = -r * mu + np.sqrt(max(r_out * r_out - perp2, 0.0))
    if mu < 0.0 and perp2 < r_core * r_core:
        s_core = -r * mu - np.sqrt(r_core * r_core - perp2)
        if s_core < s_out:
            return s_core, "core"
    return s_out, "out"


def _advance(r, mu, s):
    """New (r, mu) after flying a distance s. z = r*mu advances as z + s."""
    r_new = np.sqrt(max(r * r + s * s + 2.0 * r * s * mu, 0.0))
    if r_new <= 0.0:
        return 0.0, 1.0
    return r_new, (r * mu + s) / r_new


# --------------------------------------------------------------------------
# Transport
# --------------------------------------------------------------------------

# A resonance must lie strictly ahead. Without this a packet re-emitted at a
# line centre would find that same line at s = 0 forever. The tolerance is in
# cm and is ~1e-8 of the shell size, far below any physical scale here.
_S_MIN = 1.0


def run_mc(
    atom,
    r_core,
    r_out,
    t_exp,
    nu_min,
    nu_max,
    n_packets,
    seed=0,
    interaction="branch",
):
    """Propagate packets from the core through the shell.

    Parameters
    ----------
    atom : Atom
    r_core, r_out : cm. Homologous, so v_core = r_core / t_exp.
    t_exp : s
    nu_min, nu_max : flat launch continuum (Hz). Flat rather than Planck so
        that the emergent spectrum divided by the launch spectrum is directly
        the transmitted fraction, with no shape to divide out.
    n_packets : int
    seed : int, so every reported number is reproducible.
    interaction : 'branch'  -- sample the branching table and re-emit
                  'absorb'  -- destroy the packet on interaction. This is the
                              pure-absorption limit, which is what Paper I's
                              analytic `sobolev_attenuation` computes, so it
                              is how the transport is validated against an
                              independently written calculation.

    Returns a dict with the launched and emergent lab frequencies, the counts
    needed for the conservation identity, and the first-interaction tally used
    to measure branching ratios.
    """
    if interaction not in ("branch", "absorb"):
        raise ValueError(f"unknown interaction mode {interaction!r}")

    rng = np.random.default_rng(seed)
    tau = atom.tau(t_exp)
    ct = C * t_exp

    nu_launch = rng.uniform(nu_min, nu_max, n_packets)
    # Isotropic emission from the core surface gives a direction distribution
    # proportional to mu, hence mu = sqrt(U).
    mu_launch = np.sqrt(rng.uniform(0.0, 1.0, n_packets))

    nu_out = []
    n_core = 0
    n_absorbed = 0
    n_interacted = 0            # packets that interacted at least once
    first_branch = np.zeros(len(atom.lines), dtype=int)
    n_resonance_in_shell = 0    # packets offered at least one resonance

    for i in range(n_packets):
        r, mu, nu = r_core, mu_launch[i], nu_launch[i]
        touched = False
        offered = False

        while True:
            s_bound, kind = _distance_to_boundary(r, mu, r_core, r_out)

            # Next resonance ahead of the packet. The resonance plane z_res is
            # a property of the LAB frequency alone, so this is a direct solve
            # rather than a search along the ray.
            z = r * mu
            s_res, k_res = np.inf, -1
            for k in range(len(atom.nu0)):
                if tau[k] <= 0.0:
                    continue
                s_k = ct * (1.0 - atom.nu0[k] / nu) - z
                if _S_MIN < s_k < s_res:
                    s_res, k_res = s_k, k

            if k_res < 0 or s_res >= s_bound:
                if kind == "core":
                    n_core += 1
                else:
                    nu_out.append(nu)
                break

            offered = True
            r, mu = _advance(r, mu, s_res)
            if rng.uniform() > -np.expm1(-tau[k_res]):
                continue                      # resonance crossed untouched

            if interaction == "absorb":
                n_absorbed += 1
                if not touched:
                    n_interacted += 1
                touched = True
                break

            k_new = atom.branch(int(atom.upper[k_res]), rng.uniform())
            # Tally only each packet's FIRST interaction: subsequent ones are
            # conditioned on it, so pooling them would bias the ratio estimate.
            if not touched:
                n_interacted += 1
                first_branch[k_new] += 1
            touched = True

            # Re-emit isotropically in the comoving frame at the new line's
            # rest frequency, then transform back. To the first-order Doppler
            # convention used throughout, nu_lab = nu_cm / (1 - z/(c t)).
            mu = rng.uniform(-1.0, 1.0)
            z = r * mu
            nu = atom.nu0[k_new] / (1.0 - z / ct)

        if offered:
            n_resonance_in_shell += 1

    return {
        "nu_launch": nu_launch,
        "nu_out": np.array(nu_out),
        "n_packets": n_packets,
        "n_core": n_core,
        "n_absorbed": n_absorbed,
        "n_interacted": n_interacted,
        "n_offered": n_resonance_in_shell,
        "first_branch": first_branch,
        "tau": tau,
    }
