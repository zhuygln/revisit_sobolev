"""The macroatom as a Markov process over internal states (chapter 8), and
the two packet bookkeepings it is compared with (chapter 7).

Lucy's rules for a purely radiative, downward-only macroatom (the same
weights as sobolev/macroatom.py::DownwardMacroAtom in the production code):
from level i, the packet either DE-ACTIVATES through line i -> j with weight

    A_ij beta_ij (eps_i - eps_j)          (the line's energy),

emitting one packet of the full energy it carries at the line frequency, or
JUMPS internally to level j with weight

    A_ij beta_ij eps_j                    (the energy that stays in the atom),

and continues from j. eps_i is the level energy above the ground; beta_ij
the Sobolev escape probability of the line, applied once per emission ("beta
once"). With every beta = 1 the expected energy emitted per line equals the
energy the explicit cascade of chapter 6 puts into that line.

Two bookkeepings of the packet's energy after an activation (chapter 7):
  mode="energy": indivisible energy packet, E_out = E_in, w' = E / (h nu').
  mode="photon": photon-number packet, w photons re-emitted at the exit
                 frequency, E_out = w h nu' != E_in.
"""
import numpy as np

from .atom import H_ERG_S
from .sobolev import escape_probability


class ToyMacroAtom:
    def __init__(self, atom, beta=None, mode="energy"):
        self.atom = atom
        self.beta = np.ones(atom.n_lines) if beta is None else np.asarray(beta, float)
        assert mode in ("energy", "photon")
        self.mode = mode
        self.eps = atom.E * 1.0                       # level energies above ground (cm^-1; units cancel)

    @classmethod
    def from_line_list(cls, atom, tau, mode="energy"):
        return cls(atom, escape_probability(tau), mode)

    def probabilities(self, level):
        """(lines, p_deactivate, p_jump) from `level`: the normalised weights
        above; p_deactivate[k] + p_jump[k] summed over k is one."""
        k = self.atom.lines_from(level)
        if k.size == 0:
            return k, np.zeros(0), np.zeros(0)
        ab = self.atom.A[k] * self.beta[k]
        w_de = ab * (self.eps[level] - self.eps[self.atom.lower[k]])
        w_jp = ab * self.eps[self.atom.lower[k]]
        tot = w_de.sum() + w_jp.sum()
        return k, w_de / tot, w_jp / tot

    def activate(self, level):
        return int(level)

    def step(self, rng, level):
        """One Markov step from `level`: ('deactivate', line) or ('jump', level)."""
        k, p_de, p_jp = self.probabilities(level)
        u = rng.random()
        c = np.cumsum(np.concatenate([p_de, p_jp]))
        i = int(np.searchsorted(c, u))
        if i < k.size:
            return "deactivate", int(k[i])
        return "jump", int(self.atom.lower[k[i - k.size]])

    def walk(self, rng, level):
        """From activation at `level` to de-activation: (exit line, number of
        internal jumps)."""
        level = self.activate(level); n_jump = 0
        while True:
            kind, x = self.step(rng, level)
            if kind == "deactivate":
                return x, n_jump
            level = x; n_jump += 1

    def deactivate(self, E_in, w_in, exit_line):
        """The packet after the walk under this bookkeeping:
        (energy out, photon number out, exit frequency)."""
        nu = self.atom.nu[exit_line]
        if self.mode == "energy":
            return E_in, E_in / (H_ERG_S * nu), nu
        return w_in * H_ERG_S * nu, w_in, nu

    def run(self, rng, level, n, E_in=1.0):
        """n activations at `level`, each carrying energy E_in: the exit
        lines, the energies out, and the mean number of internal jumps."""
        exits = np.empty(n, int); E_out = np.empty(n); jumps = np.empty(n, int)
        for i in range(n):
            exits[i], jumps[i] = self.walk(rng, level)
            E_out[i] = E_in if self.mode == "energy" else self._photon_out(E_in, level, exits[i])[0]
        return exits, E_out, float(jumps.mean())

    def _photon_out(self, E_in, level, exit_line):
        nu_act = self.atom.nu[self.atom.upper == level].max()        # activated by the resonance line to the ground-most lower level
        w = E_in / (H_ERG_S * nu_act)
        return self.deactivate(E_in, w, exit_line)

    def energy_per_line(self, rng, level, n, E_in=1.0):
        """Mean energy emitted per line per activation."""
        exits, E_out, _ = self.run(rng, level, n, E_in)
        out = np.zeros(self.atom.n_lines)
        np.add.at(out, exits, E_out)
        return out / n


def cascade_energy_per_line(atom, rng, level, n, E_in=1.0):
    """The explicit cascade of chapter 6, with each emitted photon carrying
    its own h nu: the energy per line per activation, normalised so that the
    activation energy is E_in (the cascade's total is h nu of the level,
    conserved by construction)."""
    out = np.zeros(atom.n_lines)
    E_level = H_ERG_S * (atom.nu[atom.upper == level].max())                           # h nu (level -> ground) = the activation energy
    for _ in range(n):
        for j in atom.cascade(rng, level):
            out[j] += H_ERG_S * atom.nu[j]
    return out / n * (E_in / E_level)
