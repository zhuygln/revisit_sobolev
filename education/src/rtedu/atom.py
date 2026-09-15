"""A toy atom built from first principles (chapter 6).

Five levels with energies in cm^-1, every allowed downward transition with an
Einstein coefficient A_ul. The branching probability from level u is
P(u -> l) = A_ul / sum_k A_uk, and a cascade is a sequence of such draws
down to the ground level, emitting one photon per step. Boltzmann
populations at a temperature T turn the atom into a line list with Sobolev
depths (chapter 4's formula), which is the forest every later chapter
transports through.
"""
import numpy as np

from . import C, SIGMA_CLASSICAL
from .sobolev import tau_sobolev_toy

K_B_CM = 0.6950348                    # Boltzmann constant in cm^-1 per kelvin
H_ERG_S = 6.62607015e-27
F_FROM_A = 1.4992e-16                 # f_lu = 1.4992e-16 lambda[A]^2 (g_u/g_l) A_ul


class ToyAtom:
    """levels_cm: level energies (cm^-1, ground first, increasing);
    A: dict {(upper, lower): A_ul [s^-1]} of allowed downward transitions."""

    def __init__(self, levels_cm, A, name="toy"):
        self.E = np.asarray(levels_cm, float)
        assert self.E[0] == 0.0 and np.all(np.diff(self.E) > 0)
        self.n_levels = self.E.size
        self.upper = np.array([u for (u, l) in A]); self.lower = np.array([l for (u, l) in A])
        self.A = np.array([A[(u, l)] for (u, l) in A], float)
        assert np.all(self.upper > self.lower)
        self.nu = C * (self.E[self.upper] - self.E[self.lower])          # Hz
        self.lam_cm = C / self.nu
        self.f_osc = F_FROM_A * (1e8 * self.lam_cm) ** 2 * self.A          # g_u = g_l = 1
        self.n_lines = self.nu.size
        self.name = name

    # ---- level-resolved branching (chapter 6) ----
    def lines_from(self, level):
        return np.flatnonzero(self.upper == level)

    def branching(self, level):
        """(line indices, probabilities) of the downward transitions from
        `level`, P ∝ A_ul; empty for the ground level."""
        k = self.lines_from(level)
        p = self.A[k] / self.A[k].sum() if k.size else np.zeros(0)
        return k, p

    def cascade(self, rng, level):
        """One explicit cascade from `level` to the ground: the list of line
        indices emitted, one photon each."""
        out = []
        while level > 0:
            k, p = self.branching(level)
            j = k[rng.choice(k.size, p=p)]
            out.append(int(j)); level = int(self.lower[j])
        return out

    def photon_energies(self, lines):
        """h nu of each emitted photon (erg)."""
        return H_ERG_S * self.nu[np.asarray(lines, int)]

    # ---- populations and the line list (chapters 5 and later) ----
    def populations(self, T, n_total):
        """Boltzmann populations (g = 1): n_i ∝ exp(-E_i / kT)."""
        w = np.exp(-self.E / (K_B_CM * T))
        return n_total * w / w.sum()

    def line_list(self, T, n_total, t):
        """Sobolev depths of every line at temperature T, total density
        n_total (cm^-3) and epoch t (s): tau = sigma f n_l lambda t."""
        n = self.populations(T, n_total)
        return tau_sobolev_toy(n[self.lower], self.f_osc, self.lam_cm, t)

    def thermal_emissivity(self, T, n_total):
        """Energy emissivity per line, n_u A h nu, normalised: the
        distribution a thermalising interaction re-emits from."""
        n = self.populations(T, n_total)
        e = n[self.upper] * self.A * H_ERG_S * self.nu
        return e / e.sum()

    def describe(self):
        rows = []
        for k in range(self.n_lines):
            rows.append(dict(line=k, upper=int(self.upper[k]), lower=int(self.lower[k]), nm=1e7 * self.lam_cm[k],
                             A=self.A[k], f=self.f_osc[k]))
        return rows


def five_level_atom():
    """The book's atom: five levels spanning the optical to the K band, ten
    lines, one strong route 4 -> 2 -> 0 and weaker side branches."""
    levels = [0.0, 4000.0, 9000.0, 15000.0, 22000.0]
    A = {(4, 0): 2e7, (4, 1): 5e7, (4, 2): 1e8, (4, 3): 3e6,
         (3, 0): 4e7, (3, 1): 2e7, (3, 2): 8e6,
         (2, 0): 6e7, (2, 1): 1e7,
         (1, 0): 5e6}
    return ToyAtom(levels, A, name="five-level")
