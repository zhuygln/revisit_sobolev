"""Interaction count vs transmission along one ray (referee Comment 1).

The question the referee asks is whether expansion opacity reproduces the
statistic it was built for -- Karp's mean free path, the expected number of
line interactions per unit path for a photon that is counted but not removed
-- while failing the deterministic transmission of a resolved realization.
Along a single ray this is not a matter of interpretation; it is two
cumulative sums:

    S(s) = sum_{z_k < s} tau_k                expansion: E(s) = sum (1 - e^{-tau_k})
    survival  e^{-S(s)}  (Bernoulli product)  vs  e^{-E(s)}  (Poisson, mean E)

E(s) IS what alpha_exp integrates to (here drawn as SEDONA's binned ramp at
dnu/nu = 4.17e-5), so the two treatments reach the same expected count at
the far edge of the shell and different survivals -- the whole +38% in one
picture. The first-interaction pdf follows: discrete masses
e^{-S(z_k^-)}(1-e^{-tau_k}) at the resonances against the smooth
alpha_exp e^{-E}. For an absorbed photon the first-interaction probability is
1-e^{-S} vs 1-e^{-E}: expansion opacity UNDER-counts absorptions.

The resolved curves (erf, v_D = 100 and 10 km/s) sit on the Sobolev stairs
away from the edges, which is the F12 statement in the same picture.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.special import erf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from sobolev.constants import C
from sobolev.optical_depth import tau_sobolev
from sobolev.rays import RaySet
from sobolev.sobolev_leg import crossing_depths, expansion_damp, sobolev_attenuation

FOREST = ROOT / "experiments/laII_forest"
T_EXP, R_CORE, R_OUT = 86400.0, 8.64e12, 2.592e13
V_CORE, V_MAX = R_CORE / T_EXP, R_OUT / T_EXP
DNU_BIN = 4.17e-5  # SEDONA transport bin

d = np.load(FOREST / "forest_lines.npz")
N_ION = float(d["n_ion"])
lam0, f_lu, pop = d["lam"], d["f_lu"], d["pop"]
nu0 = C / (lam0 * 1e-8)
tau = np.array([tau_sobolev(f, p * N_ION, C / n0, T_EXP) for n0, f, p in zip(nu0, f_lu, pop)])

# ---- the ray: radial (p = 0). Observer nu chosen so that the ray crosses
# the LARGEST NUMBER of resonances with tau > 0.3 inside the shell (ties
# broken toward the larger summed tau), so the picture is a forest crossing
# and not a single-line one. The rule is stated so the choice is not free.
cands = []
for nu_try in np.linspace(nu0.min() / (1 - V_CORE / C), nu0.max() / (1 - V_MAX / C), 4000):
    zt = C * T_EXP * (1.0 - nu0 / nu_try)
    ins = (zt > R_CORE) & (zt < R_OUT)
    cands.append((int(np.sum(ins & (tau > 0.3))), float(tau[ins].sum()), nu_try))
cands.sort(key=lambda c: (-c[0], -c[1]))
nu_obs = cands[0][2]
z_k = C * T_EXP * (1.0 - nu0 / nu_obs)                      # resonance positions on this ray
inside = (z_k > R_CORE) & (z_k < R_OUT)
order = np.argsort(z_k)
zk, tk = z_k[order][inside[order]], tau[order][inside[order]]
# edge check: no crossed line within one transport bin of either edge
bin_len = C * T_EXP * DNU_BIN
edge_flag = np.any(np.abs(zk - R_CORE) < bin_len) or np.any(np.abs(zk - R_OUT) < bin_len)

s = np.linspace(R_CORE, R_OUT, 4000)
S_stair = np.array([tk[zk < x].sum() for x in s])
E_stair = np.array([(1.0 - np.exp(-tk[zk < x])).sum() for x in s])
# SEDONA-binned expansion ramp: alpha_b = sum_{k in b} (1-e^-tau_k) / (bin length), swept linearly
nu_com = nu_obs * (1.0 - s / (C * T_EXP))          # comoving frequency along the ray
# bin edges in comoving frequency: log-spaced with width dnu/nu
bins_edges = nu_obs * (1.0 - zk / (C * T_EXP))    # each line's resonance comoving freq = nu0
E_ramp = np.zeros_like(s)
for z0, t0 in zip(zk, tk):
    # the bin containing this line spans z in [z_b, z_b + bin_len]; place the
    # bin so the line sits at a random-but-fixed phase (here: bin starts one
    # third of a bin before the resonance), then ramp linearly across it
    zb = z0 - bin_len / 3.0
    E_ramp += (1.0 - np.exp(-t0)) * np.clip((s - zb) / bin_len, 0.0, 1.0)

def S_res(vd):
    delta = vd * T_EXP * (nu0[order][inside[order]] / nu_obs)
    return np.array([np.sum(tk * 0.5 * (erf((x - zk) / delta) - erf((R_CORE - zk) / delta))) for x in s])
S_r100, S_r10 = S_res(1.0e7), S_res(1.0e6)

# first-interaction pdf
mass = np.exp(-np.concatenate([[0.0], np.cumsum(tk)[:-1]])) * (1.0 - np.exp(-tk))   # discrete, at z_k
alpha_exp = np.gradient(E_ramp, s)
pdf_exp = alpha_exp * np.exp(-E_ramp)

# ---- integrated panel: the whole band, all rays
nu_grid = np.geomspace(7.50e14, 7.95e14, 4000)
lam_grid = C / nu_grid * 1e8
band = (lam_grid > 3800) & (lam_grid < 3955)
rays = RaySet.midpoint(R_CORE, R_OUT, 400)
lines = [(n0, f, p) for n0, f, p in zip(nu0, f_lu, pop)]
curve = []
for scale in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0):
    S, E, p, w = crossing_depths(nu_grid, lines, R_CORE, R_OUT, T_EXP, N_ION * scale, rays=rays)
    fs = np.sum(w[:, None] * np.exp(-S), axis=0)[band].mean() / w.sum()
    fe = np.sum(w[:, None] * np.exp(-E), axis=0)[band].mean() / w.sum()
    Eb, Sb = (np.sum(w[:, None] * E, axis=0)[band].mean() / w.sum(), np.sum(w[:, None] * S, axis=0)[band].mean() / w.sum())
    curve.append(dict(tau_max=5.0 * scale, f_sob=fs, f_exp=fe, d_exp=fe / fs - 1, E_over_S=Eb / Sb))
seed = json.loads((ROOT / "experiments/mc_noise/mc_noise_results.json").read_text())
sed_res, sed_exp = seed["tau5_vd100_bb"]["mean"], seed["tau5_vd100_exp"]["mean"]

# ---- figure
fig, axs = plt.subplots(2, 2, figsize=(11, 8))
sv = s / 1e5 / T_EXP  # velocity coordinate, km/s
ax = axs[0, 0]
ax.step(sv, S_stair, where="post", color="C2", lw=1.4, label=r"$S=\sum\tau_k$ (Sobolev)")
ax.plot(sv, S_r100, "C2", ls=":", lw=1, label=r"resolved, $v_D$=100")
ax.plot(sv, S_r10, "C2", ls="--", lw=0.8, alpha=0.6, label=r"resolved, $v_D$=10")
ax.step(sv, E_stair, where="post", color="C0", lw=1.4, label=r"$E=\sum(1-e^{-\tau_k})$ (expected count)")
ax.plot(sv, E_ramp, "C0", ls=":", lw=1, label=r"$\int\alpha_{\rm exp}ds$ (binned)")
ax.set_ylabel("cumulative along the ray"); ax.set_xlabel("v [km/s]"); ax.legend(fontsize=7); ax.grid(alpha=.3)
ax.set_title(f"one radial ray, {len(zk)} resonances; same count, different depth", fontsize=9)
ax = axs[0, 1]
ax.step(sv, np.exp(-S_stair), where="post", color="C2", lw=1.4, label=r"$e^{-S}$ Bernoulli product")
ax.plot(sv, np.exp(-S_r100), "C2", ls=":", lw=1)
ax.step(sv, np.exp(-E_stair), where="post", color="C0", lw=1.4, label=r"$e^{-E}$ Poisson, same mean")
ax.plot(sv, np.exp(-E_ramp), "C0", ls=":", lw=1)
ax.set_ylabel("survival P(no interaction beyond s)"); ax.set_xlabel("v [km/s]"); ax.legend(fontsize=7); ax.grid(alpha=.3)
ax.set_title(f"end of shell: {np.exp(-S_stair[-1]):.3f} vs {np.exp(-E_stair[-1]):.3f}", fontsize=9)
ax = axs[1, 0]
ax.vlines(zk / 1e5 / T_EXP, 0, mass, color="C2", lw=2, label="first-interaction mass at each resonance")
ax.plot(sv, pdf_exp * (s[1] - s[0]) * 20, "C0", lw=1, label=r"$\alpha_{\rm exp}e^{-E}$ (x20 bin)")
ax.set_xlabel("v [km/s]"); ax.set_ylabel("first-interaction probability"); ax.legend(fontsize=7); ax.grid(alpha=.3)
ax.set_title(f"P(absorbed) = {1-np.exp(-S_stair[-1]):.3f} (Bernoulli) vs {1-np.exp(-E_stair[-1]):.3f} (Poisson)", fontsize=9)
ax = axs[1, 1]
tm = [c["tau_max"] for c in curve]
ax.semilogx(tm, [100 * c["d_exp"] for c in curve], "o-", color="C0", label=r"$\Delta_{\rm exp}$ analytic = $E_w[e^{S-E}]-1$")
ax.semilogx(tm, [100 * (1 - c["E_over_S"]) for c in curve], "s--", color="gray", label=r"$1-\langle E\rangle/\langle S\rangle$")
ax.plot([5.0], [100 * (sed_exp / sed_res - 1)], "r*", ms=12, label=f"SEDONA seed mean (+{100*(sed_exp/sed_res-1):.1f}%)")
ax.set_xlabel(r"$\tau_{\max}$"); ax.set_ylabel("[%]"); ax.legend(fontsize=7); ax.grid(alpha=.3)
ax.set_title("band 3800-3955 A, all rays: the gap closes as lines weaken", fontsize=9)
fig.tight_layout()
for out in (ROOT / "outputs/fig_interaction_count.png", ROOT / "docs/figures/fig_interaction_count.png"):
    out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=170)

res = dict(nu_obs=nu_obs, n_crossed=int(len(zk)), edge_flag=bool(edge_flag),
           S_end=float(S_stair[-1]), E_end=float(E_stair[-1]), E_ramp_end=float(E_ramp[-1]),
           survival_bernoulli=float(np.exp(-S_stair[-1])), survival_poisson=float(np.exp(-E_stair[-1])),
           p_absorb_bernoulli=float(1 - np.exp(-S_stair[-1])), p_absorb_poisson=float(1 - np.exp(-E_stair[-1])),
           first_interaction_mass_sum=float(mass.sum()), pdf_exp_integral=float(np.trapezoid(pdf_exp, s)),
           curve=curve, sedona_seed_mean=dict(res=sed_res, exp=sed_exp),
           lines=[dict(v_kms=float(z / 1e5 / T_EXP), tau=float(t), one_minus_etau=float(1 - np.exp(-t))) for z, t in zip(zk, tk)])
(HERE / "ray_diagnostic.json").write_text(json.dumps(res, indent=1))
print(f"ray: lambda_obs = {C/nu_obs*1e8:.2f} A, chosen for most tau>0.3 crossings ({cands[0][0]}); {len(zk)} resonances crossed; edge flag {edge_flag}")
print(f"  S_end {S_stair[-1]:.4f}   E_end {E_stair[-1]:.4f}   binned ramp end {E_ramp[-1]:.4f}   (count preserved: {abs(E_ramp[-1]-E_stair[-1]):.1e})")
print(f"  survival  e^-S {np.exp(-S_stair[-1]):.4f}   e^-E {np.exp(-E_stair[-1]):.4f}")
print(f"  P(absorbed) Bernoulli {1-np.exp(-S_stair[-1]):.4f}  Poisson {1-np.exp(-E_stair[-1]):.4f}")
print(f"  first-interaction masses sum {mass.sum():.4f} (= 1-e^-S);  pdf_exp integral {np.trapezoid(pdf_exp, s):.4f} (= 1-e^-E)")
print("  strongest crossed lines:", ", ".join(f"tau={t:.2f}@{z/1e5/T_EXP:.0f}" for z, t in sorted(zip(zk, tk), key=lambda x: -x[1])[:5]))
print("\n  tau_max   D_exp(analytic)   1-<E>/<S>")
for c in curve: print(f"   {c['tau_max']:6.2f}    {100*c['d_exp']:+7.2f}%        {100*(1-c['E_over_S']):6.2f}%")
print(f"\n  SEDONA seed mean: res {sed_res:.4f} exp {sed_exp:.4f} -> +{100*(sed_exp/sed_res-1):.2f}%  vs analytic at tau_max=5: +{100*[c for c in curve if c['tau_max']==5.0][0]['d_exp']:.2f}%")
