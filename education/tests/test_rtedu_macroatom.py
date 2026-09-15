"""Chapters 7 and 8: the macroatom's probabilities, its energy bookkeeping,
and its equivalence to the explicit cascade."""
import numpy as np

from rtedu.atom import five_level_atom, H_ERG_S
from rtedu.macroatom import ToyMacroAtom, cascade_energy_per_line


def test_macroatom_probabilities_sum_to_one_and_use_the_lucy_weights():
    a = five_level_atom(); m = ToyMacroAtom(a)
    for level in range(1, a.n_levels):
        k, p_de, p_jp = m.probabilities(level)
        assert abs(p_de.sum() + p_jp.sum() - 1.0) < 1e-12
        w_de = a.A[k] * (a.E[level] - a.E[a.lower[k]]); w_jp = a.A[k] * a.E[a.lower[k]]
        assert np.allclose(p_de, w_de / (w_de.sum() + w_jp.sum())) and np.allclose(p_jp, w_jp / (w_de.sum() + w_jp.sum()))
    # from level 1 only the ground is below: no internal jump possible
    k, p_de, p_jp = m.probabilities(1)
    assert np.allclose(p_de, [1.0]) and np.allclose(p_jp, [0.0])


def test_macroatom_energy_conservation():
    """An indivisible energy packet leaves every activation with the energy
    it brought, whatever line it leaves by; a photon-number packet does not."""
    a = five_level_atom(); rng = np.random.default_rng(3)
    exits, E_out, _ = ToyMacroAtom(a, mode="energy").run(rng, 4, 2000, E_in=1.0)
    assert np.all(E_out == 1.0) and len(set(exits.tolist())) > 1
    exits_p, E_p, _ = ToyMacroAtom(a, mode="photon").run(rng, 4, 2000, E_in=1.0)
    assert np.any(E_p != 1.0) and abs(E_p.mean() - 1.0) > 0.1              # it loses energy on average here


def test_macroatom_reproduces_the_cascade_energy_per_line():
    """With every beta = 1 the expected energy the macroatom emits per line
    equals the energy the explicit cascade puts into that line (Lucy's
    construction), within 4 sigma of the Monte Carlo noise."""
    a = five_level_atom(); n = 40_000
    e_m = ToyMacroAtom(a).energy_per_line(np.random.default_rng(4), 4, n)
    e_c = cascade_energy_per_line(a, np.random.default_rng(5), 4, n)
    assert abs(e_m.sum() - 1.0) < 1e-9 and abs(e_c.sum() - 1.0) < 1e-9
    sig = np.sqrt(np.maximum(e_m, 1e-6) / n) * 2 + 2e-3           # crude per-line noise bound
    assert np.all(np.abs(e_m - e_c) < 4 * sig)


def test_beta_once_changes_the_exit_distribution():
    """A trapped line (small beta) is left less often: the escape
    probability enters the weights once per emission."""
    a = five_level_atom(); rng = np.random.default_rng(6)
    beta = np.ones(a.n_lines); beta[2] = 0.05                       # trap 4 -> 2
    free = ToyMacroAtom(a).energy_per_line(rng, 4, 20_000)
    trapped = ToyMacroAtom(a, beta=beta).energy_per_line(rng, 4, 20_000)
    assert trapped[2] < 0.5 * free[2] and abs(trapped.sum() - 1.0) < 1e-9
