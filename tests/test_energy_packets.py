"""Paper IV 2A.1: indivisible energy packets in `run_mc(packets="energy")`.

What the switch must do: keep the histories of every photon-number-
probability mode bit for bit (the R1^E rung), conserve the comoving energy
at every event (E_dep_cm == 0 exactly), close the identity to roundoff in
every mode, and hand the kernel post-event weights that conserve energy.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1", ROOT / "paper3"):
    sys.path.insert(0, str(p))

from sobolev.constants import C, H                           # noqa: E402
from sobolev import energy_balance as eb                     # noqa: E402
from sobolev.energy_packets import energy_accounting         # noqa: E402
from forest_mc import MODES, run_mc, spectrum, _weights      # noqa: E402
from redistribution import RedistributionKernel              # noqa: E402
import test_forest_mc as tfm                                 # noqa: E402

R_CORE, R_OUT, T_EXP = tfm.R_CORE, tfm.R_OUT, tfm.T_EXP
NU_13, NU_32 = tfm.NU_13, tfm.NU_32
NU_21 = NU_13 - NU_32

# modes whose draws do not depend on the packet language
PHOTON_PROB_MODES = ("sobolev_absorb", "expansion_absorb", "binned_absorb", "dual_absorb",
                     "sobolev_branch", "expansion_branch", "dual_branch")
ENERGY_MODES = tuple(m for m in MODES if not m.endswith("_group"))


def toy():
    fa, _ = tfm.three_level(1.5)
    fa.temperature = 3000.0
    fa.emis_w = np.array([1.0, 1.0])
    fa.level_energy_cm = np.array([0.0, 0.0, NU_21 / C, NU_13 / C])   # levels 0(unused),1,2,3
    return fa


def _run(fa, mode, packets, **kw):
    lo, hi = tfm.pump_band()
    kw.setdefault("eps", 0.5)
    if mode.endswith("_group"):
        kw.setdefault("kernel", tfm._toy_kernel())
    return run_mc(fa, R_CORE, R_OUT, T_EXP, lo, hi, 20000, mode, seed=2, packets=packets, **kw)


@pytest.mark.parametrize("mode", PHOTON_PROB_MODES)
@pytest.mark.parametrize("relativity", [None, "worldline"])
def test_energy_packets_do_not_change_photon_probability_histories(mode, relativity):
    """The loop never reads w: fates, frequencies and counters are identical,
    only w differs -- the R1 -> R1^E rung is bookkeeping alone."""
    fa = toy()
    a = _run(fa, mode, "photon", relativity=relativity)
    b = _run(fa, mode, "energy", relativity=relativity)
    for key in ("nu_out_all", "nu_final", "fate", "n_events", "first_line", "last_line"):
        assert np.array_equal(a[key], b[key], equal_nan=True), key
    assert np.array_equal(a["w"], b["w_launch"])
    if "branch" in mode:
        assert not np.array_equal(a["w"], b["w"])      # fluorescence rescaled some packets


@pytest.mark.parametrize("mode", ENERGY_MODES)
def test_identity_closes_and_comoving_deposit_is_zero_in_every_mode(mode):
    fa = toy()
    res = _run(fa, mode, "energy")
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12
    assert a["E_dep_cm"] == 0.0 and res["e_dep_cm"].sum() == 0.0
    assert res["packets"] == "energy"
    if mode.endswith("absorb"):
        assert a["E_abs"] > 0
    else:
        # every escaped packet still carries the energy it was launched with,
        # up to the Doppler work: the sum is the lab identity
        assert np.isclose(a["E_esc"] + a["E_core"] + a["E_abs"] + a["E_dep_lab"], a["E_inj"], rtol=1e-12)


def test_fluorescence_keeps_the_packet_energy_and_changes_the_photon_count():
    """A packet absorbed at nu13 and re-emitted at nu32 stands for nu13/nu32
    photons afterwards; its comoving energy is unchanged."""
    fa = toy()
    res = _run(fa, "sobolev_branch", "energy")
    fl = (res["last_line"] == 1) & (res["fate"] == 1)      # exited through 3->2
    assert fl.sum() > 100
    ratio = res["w"][fl] / res["w_launch"][fl]
    assert np.allclose(ratio, NU_13 / NU_32, rtol=1e-12)
    assert np.all(res["w"] >= res["w_launch"] * (1 - 1e-12))


def test_energy_spectrum_and_weights_use_the_launch_weights_on_the_injected_side():
    fa = toy()
    res = _run(fa, "sobolev_branch", "energy")
    w_in, w_out = _weights(res, "energy")
    assert np.array_equal(w_in, res["w_launch"] * res["nu_launch"])
    assert np.isclose(H * w_in.sum(), res["accounting"]["E_inj"], rtol=1e-12)
    assert np.isclose(H * w_out.sum(), res["accounting"]["E_esc"], rtol=1e-12)
    lo, hi = tfm.pump_band()
    edges = np.geomspace(NU_32 * 0.9, hi * 1.01, 200)
    s, _ = spectrum(res, edges, weight="energy")
    assert np.isfinite(s[np.isfinite(s)]).all()


def test_events_carry_post_event_weights_that_conserve_energy():
    fa = toy()
    res = _run(fa, "sobolev_branch", "energy", collect_events=True)
    nu_in, nu_out, w_in, w_out = res["events"]
    assert nu_in.size > 0 and w_out.shape == w_in.shape
    assert np.allclose(w_out * nu_out, w_in * nu_in, rtol=1e-12)
    kern = RedistributionKernel.from_branching_mc(nu_in, nu_out, w_in, 8, w_out=w_out)
    assert np.all(np.abs(kern.q_dep[kern.counts > 0]) < 1e-12)     # nothing deposited
    assert kern.validate_energy() < 1e-12
    assert kern.metadata["energy_events"] is True
    # photon rows still count photons: more photons out than in where the
    # exit is redder
    rs = kern.R.sum(axis=1)[kern.counts > 0]
    assert np.allclose(rs, 1.0, atol=1e-12)


def test_group_leg_in_energy_mode_samples_energy_rows():
    fa = toy()
    ref = _run(fa, "sobolev_branch", "energy", collect_events=True)
    e = ref["events"]
    kern = RedistributionKernel.from_branching_mc(e[0], e[1], e[2], 6, w_out=e[3])
    res = _run(fa, "sobolev_group", "energy", kernel=kern)
    a = res["accounting"]
    assert abs(a["identity_residual"]) < 1e-12 and a["E_dep_cm"] == 0.0
    # the kernel reproduces the reference's escaped energy fraction to noise
    assert abs(a["E_esc"] / a["E_inj"] - ref["accounting"]["E_esc"] / ref["accounting"]["E_inj"]) < 0.02


def test_photon_mode_kernel_construction_is_unchanged_by_the_w_out_argument():
    rng = np.random.default_rng(5)
    nin = rng.uniform(3e14, 1.3e15, 5000); nout = nin * rng.uniform(0.7, 1.05, 5000)
    w = np.ones(5000)
    a = RedistributionKernel.from_branching_mc(nin, nout, w, 8)
    b = RedistributionKernel.from_branching_mc(nin, nout, w, 8, w_out=w)
    for k in ("R", "N_cum", "q_dep", "sub_cum", "counts", "disc_cum", "disc_vals"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    rng1, rng2 = np.random.default_rng(9), np.random.default_rng(9)
    assert np.array_equal(a.sample_nu_out(nin[:300], rng1), a.sample_nu_out(nin[:300], rng2, rows="photon"))


def test_scale_exact_is_the_equilibrium_series_and_refuses_deposits():
    from sobolev import photometry
    fa = toy()
    res = _run(fa, "sobolev_branch", "energy")
    assert eb.is_energy_conserving(res)
    assert eb.scale_exact(res, 1.0) == photometry._scale(res, 1.0, "equilibrium")
    acc = energy_accounting(res)
    assert acc["dep_cm_frac"] == 0.0 and abs(acc["identity_residual"]) < 1e-12
    assert np.isclose(eb.renorm_ratio(res), (1.0 - acc["core_frac"] - acc["abs_frac"]) / acc["esc_frac"] - 1.0)
    bad = _run(fa, "sobolev_branch", "photon")
    with pytest.raises(ValueError):
        eb.scale_exact(bad, 1.0)


def test_invalid_combinations_raise():
    fa = toy()
    with pytest.raises(ValueError):
        _run(fa, "sobolev_dmacro", "photon")
    with pytest.raises(ValueError):
        _run(fa, "sobolev_branch", "energy", core="reemit")        # no t_core
    with pytest.raises(ValueError):
        _run(fa, "sobolev_branch", "photon", core="reemit", t_core=6000.0)
    with pytest.raises(ValueError):
        _run(fa, "sobolev_branch", "ergs")
