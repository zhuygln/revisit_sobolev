"""Chapter 13: the three-dimensional toy transport's geometry and accounting."""
import numpy as np

from rtedu import C, DAY
from rtedu.atom import five_level_atom
from rtedu.matrix import MacroatomRedistribution
from rtedu.redistribution import EpsilonRedistribution
from rtedu.transport3d import ToyEjecta, launch, transport, light_curve


def test_a_lineless_sphere_releases_its_packets_over_the_light_crossing_time():
    """With every line removed (tau = 0 through n0 = 0) nothing interacts:
    every packet escapes, none later than the light-crossing time of the
    sphere at its launch, and the frequencies are unchanged."""
    a = five_level_atom(); ej = ToyEjecta(a, n0=0.0)
    rng = np.random.default_rng(1)
    x, d, nu, t = launch(rng, ej, 300, 100, 3 * DAY)
    t_esc, nu_esc, n_int = transport(rng, ej, x, d, nu, t, lambda tt: EpsilonRedistribution(1.0, np.ones(a.n_lines)))
    assert np.all(np.isfinite(t_esc)) and np.all(n_int == 0) and np.allclose(nu_esc, nu)
    assert np.all(t_esc > t) and np.all(t_esc - t <= 2 * ej.r_out(t) / C + 1e-6)


def test_lines_trap_packets_and_the_light_curve_conserves_energy():
    a = five_level_atom(); ej = ToyEjecta(a, n0=30.0)
    rng = np.random.default_rng(2)
    x, d, nu, t = launch(rng, ej, 200, 100, 3 * DAY)
    models = {}
    def at(tt):
        key = round(tt / DAY, 2)
        if key not in models:
            models[key] = MacroatomRedistribution(a, ej.tau(tt))
        return models[key]
    t_esc, nu_esc, n_int = transport(rng, ej, x, d, nu, t, at)
    assert n_int.mean() > 0.5 and np.isfinite(t_esc).all()
    edges = np.linspace(ej.t0, 6 * DAY, 9)
    lc = light_curve(t_esc, edges)
    assert abs((lc * np.diff(edges)).sum() - np.sum(t_esc < 6 * DAY)) < 1e-9     # energy in the bins = packets that escaped within them
    assert np.all(t_esc >= t)
