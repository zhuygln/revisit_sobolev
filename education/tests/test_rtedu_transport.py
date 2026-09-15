"""The radial toy transport with a redistribution hook (chapters 5-12)."""
import numpy as np

from rtedu import C
from rtedu.atom import five_level_atom
from rtedu.redistribution import EpsilonRedistribution
from rtedu.transport import run, sweep, emergent_by_line


def test_coherent_scattering_keeps_the_line_identity_at_every_interaction():
    """With eps = 0 every interaction re-emits in the line that absorbed
    (each history pair is (k, k)); the packet then continues its sweep and
    can still meet lower lines, so the emergent line is not the launch
    line in general: frequency changes only through the Doppler sweep."""
    a = five_level_atom(); t = 2 * 86400.0; r_out = 0.2 * C * t
    tau = a.line_list(4000.0, 30.0, t); emis = a.thermal_emissivity(4000.0, 30.0)
    rng = np.random.default_rng(7); n_hist = 0
    for _ in range(300):
        nu, n_int, hist = sweep(rng, a.nu[0] * 1.001, a.nu, tau, r_out, t, EpsilonRedistribution(0.0, emis))
        assert all(k == j for k, j in hist); n_hist += len(hist)
        assert nu <= a.nu[0] * 1.001 + 1e-6
    assert n_hist > 0
    nu, last, n_int = run(rng, a.nu[0] * 1.001, 200, a.nu, tau, r_out, t, EpsilonRedistribution(0.0, emis))
    assert abs(emergent_by_line(last, a.n_lines).sum() - 1.0) < 1e-12


def test_thermal_redistribution_moves_energy_to_other_lines():
    a = five_level_atom(); t = 2 * 86400.0; r_out = 0.2 * C * t
    tau = a.line_list(4000.0, 30.0, t); emis = a.thermal_emissivity(4000.0, 30.0)
    rng = np.random.default_rng(8)
    nu, last, n_int = run(rng, a.nu[0] * 1.001, 500, a.nu, tau, r_out, t, EpsilonRedistribution(1.0, emis))
    frac = emergent_by_line(last, a.n_lines)
    assert (frac[:a.n_lines] > 0).sum() >= 3 and n_int.max() >= 1
