"""Monte Carlo confirmation of the ray diagnostic: Bernoulli vs Poisson.

Same ray, same resonances. Two tallies, each 2e6 photons:
  tally mode -- photons counted, not removed: Bernoulli(1-e^{-tau_k}) at each
                resonance vs Poisson(E_bin) per bin; both means must equal E.
  kill  mode -- photons removed at first interaction: first k with U < p_k
                (survival e^{-S}) vs exponential in the integrated alpha_exp
                (survival e^{-E}).
If the deterministic picture is right, the means agree to N^{-1/2} and the
survivals split exactly as e^{-S} vs e^{-E}.
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
d = json.loads((HERE / "ray_diagnostic.json").read_text())
tau = np.array([l["tau"] for l in d["lines"]])
p = 1.0 - np.exp(-tau)
S, E = tau.sum(), p.sum()
rng = np.random.default_rng(11)
N = 2_000_000

# tally
U = rng.random((N, tau.size))
count_bern = (U < p[None, :]).sum(axis=1).mean()
count_pois = rng.poisson(p[None, :] * np.ones((N, 1))).sum(axis=1).mean()
# kill
hit = U < p[None, :]
surv_bern = 1.0 - hit.any(axis=1).mean()
surv_pois = np.mean(rng.exponential(1.0, N) > E)   # first-interaction depth tau_r ~ Exp(1) beyond E

se = 1 / np.sqrt(N)
print(f"{tau.size} resonances, S = {S:.4f}, E = {E:.4f}")
print(f"tally:  Bernoulli mean count {count_bern:.4f}   Poisson mean count {count_pois:.4f}   target E = {E:.4f}   (se ~ {np.sqrt(E/N):.4f})")
print(f"kill:   Bernoulli survival {surv_bern:.4f} (e^-S = {np.exp(-S):.4f})   Poisson survival {surv_pois:.4f} (e^-E = {np.exp(-E):.4f})   (se ~ {se:.4f})")
assert abs(count_bern - E) < 4 * np.sqrt(E / N) and abs(count_pois - E) < 4 * np.sqrt(E / N)
assert abs(surv_bern - np.exp(-S)) < 4 * se and abs(surv_pois - np.exp(-E)) < 4 * se
print("OK: same expected count, different survival, at MC precision")
(HERE / "mc_check.json").write_text(json.dumps(dict(S=S, E=E, count_bern=count_bern, count_pois=count_pois,
    surv_bern=surv_bern, surv_pois=surv_pois, exp_S=np.exp(-S), exp_E=np.exp(-E), N=N), indent=1))
