"""Paper III plan section 21: the redistribution-kernel test battery."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper3"))
sys.path.insert(0, str(ROOT / "paper2/phase1"))
from sobolev.constants import C
from redistribution import RedistributionKernel
from forest_mc import ForestAtom, band_ratio, run_mc

T_EXP = 86400.0
CT = C * T_EXP
R_CORE, R_OUT = CT / 300.0, CT / 100.0


def _events(n=50000, seed=0, identity=False):
    rng = np.random.default_rng(seed)
    nin = rng.uniform(4e14, 9e14, n)
    nout = nin.copy() if identity else nin * rng.uniform(0.6, 1.1, n)
    return nin, nout, np.ones(n)


def test_energy_conservation_rows():
    k = RedistributionKernel.from_branching_mc(*_events(), 32)
    assert k.validate_energy() < 1e-12
    # q_dep can be negative (net blueward fluorescence); rows still close
    nin, nout, w = _events()
    k2 = RedistributionKernel.from_branching_mc(nin, nin * 1.05, w, 16)
    assert k2.validate_energy() < 1e-12 and (k2.q_dep[~k2.empty_rows] < 0).all()


def test_identity_kernel_reproduces_coherent_scattering():
    """R_ij = delta_ij (trained on nu_out == nu_in events, fine groups, pdf
    sampling) must match the resonant-scattering leg (sobolev_tla eps=0,
    which re-emits at the absorbed line's own rest frequency) within noise."""
    atom = ForestAtom(nu0=np.array([6.0e14, 7.0e14]), f_osc=np.array([0.01, 0.01]),
                      n_lower=np.array([3.0e3, 3.0e3]), n_upper=np.array([0.0, 0.0]),
                      A=np.array([1e8, 1e8]), lower=np.array([0, 0]), upper=np.array([1, 2]),
                      t_exp=T_EXP, stim=False)
    atom.emis_w = np.ones(2); atom.temperature = 3000.0   # tla builds a thermal sampler even at eps=0
    lo, hi = 5.8e14, 7.4e14
    # train the identity kernel on the branch... simplest: synthetic identity
    # events at the two line frequencies (that is where absorption happens)
    nin = np.repeat([6.0e14, 7.0e14], 20000)
    k = RedistributionKernel.from_branching_mc(nin, nin.copy(), np.ones(nin.size), 128,
                                               nu_lo=lo, nu_hi=hi)
    f = {}
    for mode, kw in (("sobolev_tla", {"eps": 0.0}), ("sobolev_group", {"kernel": k})):
        res = run_mc(atom, R_CORE, R_OUT, T_EXP, lo, hi, 200_000, mode, seed=3, **kw)
        f[mode] = (res["fate"] == 1).mean()
    assert abs(f["sobolev_group"] - f["sobolev_tla"]) < 5e-3, f


def test_rebin_invariance():
    """Coarse kernel from events == energy-weighted aggregation of the
    nested fine kernel's flows (exact, same events, nested log edges)."""
    nin, nout, w = _events(seed=4)
    lo, hi = 3.9e14, 9.3e14
    fine = RedistributionKernel.from_branching_mc(nin, nout, w, 64, nu_lo=lo, nu_hi=hi)
    coarse = RedistributionKernel.from_branching_mc(nin, nout, w, 8, nu_lo=lo, nu_hi=hi)
    # aggregate fine: E_in per fine row, flows E_in_row * R_row summed in 8x8 blocks
    gi = fine.group_index(nin)
    E_in_fine = np.bincount(gi, weights=w * nin, minlength=64)
    flow_fine = E_in_fine[:, None] * fine.R
    agg = flow_fine.reshape(8, 8, 8, 8).sum(axis=(1, 3))
    E_in_coarse = E_in_fine.reshape(8, 8).sum(axis=1)
    R_agg = np.where(E_in_coarse[:, None] > 0, agg / np.maximum(E_in_coarse[:, None], 1e-300), 0.0)
    assert np.abs(R_agg - coarse.R).max() < 1e-12


def test_empty_rows_fall_back_to_coherent_and_zero_opacity_limit():
    nin, nout, w = _events(n=1000, seed=5)
    k = RedistributionKernel.from_branching_mc(nin, nout, w, 8, nu_lo=1e14, nu_hi=2e15)
    assert k.empty_rows.any()
    rng = np.random.default_rng(0)
    s = k.sample_nu_out(np.array([1.05e14]), rng)   # far below every event
    assert np.isnan(s[0])                            # -> transport scatters coherently
    # zero opacity: no interactions, kernel never invoked
    atom = ForestAtom(nu0=np.array([6.0e14]), f_osc=np.array([0.01]),
                      n_lower=np.array([0.0]), n_upper=np.array([0.0]),
                      A=np.array([1e8]), lower=np.array([0]), upper=np.array([1]),
                      t_exp=T_EXP, stim=False, tau_min=-1.0)
    res = run_mc(atom, R_CORE, R_OUT, T_EXP, 5.8e14, 6.2e14, 20_000, "sobolev_group",
                 seed=6, kernel=k)
    assert res["n_interactions"] == 0 and res["n_escaped"] + res["n_core"] == 20_000


def test_save_load_roundtrip(tmp_path):
    k = RedistributionKernel.from_branching_mc(*_events(seed=7), 16)
    p = tmp_path / "k.npz"; k.save(p)
    k2 = RedistributionKernel.load(p)
    for a, b in ((k.R, k2.R), (k.N_cum, k2.N_cum), (k.q_dep, k2.q_dep), (k.sub_cum, k2.sub_cum)):
        assert np.array_equal(a, b)
    assert k2.metadata["n_events"] == k.metadata["n_events"]
