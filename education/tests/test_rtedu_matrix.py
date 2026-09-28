"""Part III: the redistribution matrix built from the macroatom, the models
that use it, and the surrogate."""
import numpy as np

from rtedu import C
from rtedu.atom import five_level_atom
from rtedu.matrix import group_of_line, MacroatomRedistribution, build_R, MatrixRedistribution, low_rank, interpolate_R, row_error
from rtedu.transport import run, emergent_by_line
from rtedu.surrogate import MLP, dataset, split_by_state, predict_R


def _state(T=4000.0, n_tot=30.0):
    a = five_level_atom(); t = 2 * 86400.0
    return a, t, 0.2 * C * t, a.line_list(T, n_tot, t), a.thermal_emissivity(T, n_tot)


def test_groups_are_contiguous_in_frequency_and_capped_at_the_line_count():
    a = five_level_atom()
    g, ng = group_of_line(a.nu, 4)
    order = np.argsort(a.nu)
    assert ng == 4 and np.all(np.diff(g[order]) >= 0) and set(g) == {0, 1, 2, 3}
    g10, ng10 = group_of_line(a.nu, 32)
    assert ng10 == a.n_lines and len(set(g10)) == a.n_lines


def test_R_rows_sum_to_one():
    a, t, r_out, tau, emis = _state()
    macro = MacroatomRedistribution(a, tau)
    run(np.random.default_rng(1), a.nu[0] * 1.001, 400, a.nu, tau, r_out, t, macro)
    for n_g in (1, 2, 4, 10):
        g, ng = group_of_line(a.nu, n_g); R = build_R(macro.events, g, ng)
        assert R.shape == (ng, ng) and np.allclose(R.sum(axis=1), 1.0) and np.all(R >= 0)
    R1 = build_R(macro.events, *group_of_line(a.nu, 1)); assert R1[0, 0] == 1.0


def test_R_reproduces_macroatom_distribution():
    """At N_g = 10 every group is one line, nothing is discarded, and the
    transport with R alone reproduces the macroatom's emergent line
    spectrum within the Monte Carlo noise; at N_g = 2 it does not."""
    a, t, r_out, tau, emis = _state(); n = 3000
    macro = MacroatomRedistribution(a, tau)
    _, last_ref, _ = run(np.random.default_rng(2), a.nu[0] * 1.001, n, a.nu, tau, r_out, t, macro)
    ref = emergent_by_line(last_ref, a.n_lines)
    errs = {}
    for n_g in (10, 2):
        g, ng = group_of_line(a.nu, n_g); R = build_R(macro.events, g, ng)
        _, last, _ = run(np.random.default_rng(3), a.nu[0] * 1.001, n, a.nu, tau, r_out, t, MatrixRedistribution(R, g, emis))
        errs[n_g] = np.abs(emergent_by_line(last, a.n_lines) - ref).sum()
    noise = 4 * np.sqrt(ref * (1 - ref) / n).sum()           # a generous L1 noise bound
    assert errs[10] < noise, (errs, noise)
    assert errs[2] > errs[10]


def test_low_rank_and_interpolation_keep_rows_normalised():
    R = np.array([[0.7, 0.2, 0.1], [0.1, 0.8, 0.1], [0.05, 0.15, 0.8]])
    A = low_rank(R, 1); assert np.allclose(A.sum(axis=1), 1.0) and np.all(A >= 0)
    assert np.allclose(low_rank(R, 3), R, atol=1e-10)
    R2 = np.eye(3)
    Ri = interpolate_R(np.sqrt(3000.0 * 5000.0), [3000.0, 5000.0], [R, R2])
    assert np.allclose(Ri, 0.5 * (R + R2)) and np.allclose(Ri.sum(axis=1), 1.0)
    assert np.allclose(interpolate_R(1000.0, [3000.0, 5000.0], [R, R2]), R)
    assert row_error(R, R) == 0.0 and row_error(R, R2) > 0


def test_surrogate_learns_a_table_and_holds_out_whole_states():
    Ts = [2500.0, 3000.0, 4000.0, 5000.0, 6000.0]
    n_g = 3
    Rs = [np.array([[0.9 - 0.1 * k, 0.1 * k, 0.0], [0.2, 0.6, 0.2], [0.0, 0.1 + 0.05 * k, 0.9 - 0.05 * k]]) for k in range(5)]
    train, held = split_by_state(Ts, held_out=[4000.0])
    assert held == [2] and train == [0, 1, 3, 4]
    X, Y = dataset([Ts[i] for i in train], [Rs[i] for i in train], n_g)
    m = MLP(X.shape[1], 16, n_g, seed=0); hist = m.fit(X, Y, epochs=1500, lr=0.3)
    assert hist[-1] < hist[0]
    P = predict_R(m, 4000.0, n_g)
    assert np.allclose(P.sum(axis=1), 1.0) and row_error(P, Rs[2]) < 0.1
