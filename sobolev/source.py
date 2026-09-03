"""A physically powered kilonova source: heating, thermalization, one-zone
diffusion, photosphere.

Everything before 2026-09-02 illuminated the line-blanketed shell with an
*imposed* blackbody (T = 5000 K t^-1/2, R = v_core t). On that source the ejecta
mass and the lanthanide fraction enter the transport only through the product
n_ion ~ rho X_lan, so the question the action plan asks -- does the closure
error point along a physical-parameter direction? -- cannot be posed:
dm/dlnM and dm/dlnX are identical by construction.

Here (M_ej, v_ej) set the luminosity and the photosphere, and X_lan sets only
the line opacity. That is the Tier-1 control the user chose: the grey
diffusion opacity `kappa` is FIXED and independent of X_lan, so the lanthanide
fraction is not counted twice (once in the diffusion time, once in the lines).
X_lan is therefore an *effective* or tracer lanthanide abundance in this grid;
a lanthanide-dependent kappa (xkn-like) is the deferred Tier 2.

Ingredients, all cgs (t in s, m in g, v in cm/s):

  heating_rate      eps(t) = 4e18 [1/2 - arctan((t - t0)/sig)/pi]^1.3 erg/g/s
                    (Korobkin+2012, the form Villar+2017 use), composition-
                    independent by design.
  thermalization    Barnes+2016 eq. 35 with (a, b, d) bilinear in
                    (log10 M, v) over their Table 1, as MOSFiT does.
  arnett_luminosity one-zone Arnett (1982) solution,
                    L = (2/tau_d^2) e^{-(t/tau_d)^2} int_0^t t' e^{(t'/tau_d)^2} Qdot dt',
                    tau_d = sqrt(2 kappa M / (beta c v)), beta = 13.7, solved
                    as the ODE dL/dt = (2t/tau_d^2)(Qdot - L) with an exact
                    per-segment exponential integrator on a log-t grid.
  photosphere       R_ph where the GREY optical depth measured inward from
                    the edge of the uniform sphere reaches 2/3 with the same
                    kappa the diffusion uses -- kappa rho (R_out - R_ph) = 2/3
                    -- so the transport launches from the surface the
                    diffusion model itself defines; T_eff = [L/(4 pi sigma
                    R_ph^2)]^{1/4}, T_gas = T_eff. The plan's first
                    convention, v_ph = v_ej/2, is kept as `v_ph_frac=0.5`; the
                    2026-09-02 dry run showed it puts 7/8 of the ejecta mass
                    above the launch surface, a line shell of tau ~ 1e3 that
                    no packet crosses within max_steps. `v_ph_frac="grey"`
                    is the default. The grey photosphere recedes as rho t^-2
                    thins; below `v_ph_min` x v_ej it is floored and flagged.

`SourceModel(m, v).state(t_d)` returns the dict `trajectory.state` returns
(t_exp, rho, T_gas, t_core, r_core, r_out) plus the source diagnostics, so the
transport code is called exactly as before: launch from r_core = R_ph at
t_core = T_eff.

Declared caveats (report them; do not fix them here): the central model has
T_eff = 9700 K at 0.5 d and 7400 K at 1 d, where ION_FRAC = 1 (no Saha)
over-populates the singly ionized stage, and LTE populations at T_gas = T_eff
sample the Boltzmann tail heavily at early epochs.
"""
import numpy as np
from scipy.special import erfcx

from .constants import C

MSUN = 1.98847e33          # g
DAY = 86400.0              # s
SIGMA_SB = 5.670374419e-5  # erg cm^-2 s^-1 K^-4
BETA_ARNETT = 13.7
KAPPA_SRC = 1.0            # cm^2/g, the fixed grey diffusion opacity (Tier 1)
V_PH_FRAC = "grey"         # v_ph / v_ej: a number, or "grey" for tau_grey = 2/3
V_PH_MIN = 0.5             # floor on v_ph / v_ej once the grey photosphere recedes (the plan's v_ej/2)
TAU_PH = 2.0 / 3.0

# Barnes+2016 Table 1: (a, b, d) at M/Msun in rows, v/c in columns.
_BARNES_M = np.array([1e-3, 5e-3, 1e-2, 5e-2])
_BARNES_V = np.array([0.1, 0.2, 0.3])
_BARNES = np.array([
    [[2.01, 0.28, 1.12], [4.52, 0.62, 1.39], [8.16, 1.19, 1.52]],
    [[0.81, 0.19, 0.86], [1.90, 0.28, 1.21], [3.20, 0.45, 1.39]],
    [[0.56, 0.17, 0.74], [1.31, 0.21, 1.13], [2.19, 0.31, 1.32]],
    [[0.27, 0.10, 0.60], [0.55, 0.13, 0.90], [0.95, 0.15, 1.13]],
])


def heating_rate(t, eps0=4.0e18, t0=1.3, sig=0.11, alpha=1.3):
    """Specific r-process heating rate, erg g^-1 s^-1 (Korobkin+2012)."""
    t = np.asarray(t, float)
    return eps0 * (0.5 - np.arctan((t - t0) / sig) / np.pi) ** alpha


def barnes_params(m_msun, v_c):
    """(a, b, d) of Barnes+2016 eq. 35, bilinear in (log10 M, v) over Table 1.

    Returns the parameters and whether the inputs were clamped to the table.
    """
    lm = np.log10(m_msun)
    lms = np.log10(_BARNES_M)
    clamped = not (lms[0] <= lm <= lms[-1] and _BARNES_V[0] <= v_c <= _BARNES_V[-1])
    lm = np.clip(lm, lms[0], lms[-1])
    v = np.clip(v_c, _BARNES_V[0], _BARNES_V[-1])
    i = int(np.clip(np.searchsorted(lms, lm) - 1, 0, len(lms) - 2))
    j = int(np.clip(np.searchsorted(_BARNES_V, v) - 1, 0, len(_BARNES_V) - 2))
    fm = (lm - lms[i]) / (lms[i + 1] - lms[i])
    fv = (v - _BARNES_V[j]) / (_BARNES_V[j + 1] - _BARNES_V[j])
    p = ((1 - fm) * (1 - fv) * _BARNES[i, j] + fm * (1 - fv) * _BARNES[i + 1, j]
         + (1 - fm) * fv * _BARNES[i, j + 1] + fm * fv * _BARNES[i + 1, j + 1])
    return tuple(float(x) for x in p), bool(clamped)


def thermalization(t, m_msun, v_c):
    """Barnes+2016 eq. 35 thermalization efficiency f_th(t); t in s."""
    (a, b, d), _ = barnes_params(m_msun, v_c)
    td = np.asarray(t, float) / DAY
    x = 2.0 * b * td ** d
    return 0.36 * (np.exp(-a * td) + np.log1p(x) / x)


def deposited_power(t, m_msun, v_c):
    """Qdot(t) = M eps(t) f_th(t), erg/s."""
    return m_msun * MSUN * heating_rate(t) * thermalization(t, m_msun, v_c)


def diffusion_time(m_msun, v_c, kappa=KAPPA_SRC):
    """tau_d = sqrt(2 kappa M / (beta c v)), s."""
    return np.sqrt(2.0 * kappa * m_msun * MSUN / (BETA_ARNETT * C * v_c * C))


def arnett_luminosity(t, qdot, tau_d):
    """One-zone Arnett luminosity on the grid `t` (ascending, s) for deposition
    `qdot` (erg/s, same grid).

    Solves dL/dt = (2t/tau_d^2)(Qdot - L) exactly per segment with Qdot held
    at its segment mean: L_{i+1} = Q_i + (L_i - Q_i) exp(-(t_{i+1}^2 - t_i^2)/tau_d^2).
    Unconditionally stable, no exponential overflow, and L -> Qdot as tau_d -> 0.
    """
    t = np.asarray(t, float); qdot = np.asarray(qdot, float)
    L = np.empty_like(t)
    L[0] = -qdot[0] * np.expm1(-(t[0] / tau_d) ** 2)      # constant-Qdot start
    dt2 = np.diff(t ** 2) / tau_d ** 2
    qm = 0.5 * (qdot[1:] + qdot[:-1])
    g = -np.expm1(-dt2)                                    # 1 - e^{-dz}, exact for tiny dz
    for i in range(len(t) - 1):
        L[i + 1] = L[i] + (qm[i] - L[i]) * g[i]
    return L


def radiated_fraction(z):
    """eta(z) = sqrt(pi) z erfcx(z): the fraction of energy deposited at
    t' = z tau_d that is ever radiated (the rest is spent on adiabatic
    expansion). eta(0) = 0, eta(1) = 0.758, eta -> 1 as z -> inf."""
    z = np.asarray(z, float)
    return np.sqrt(np.pi) * z * erfcx(z)


def time_grid(t_min=1e-2, t_max_d=30.0, n=3000):
    return np.geomspace(t_min, t_max_d * DAY, n)


class SourceModel:
    """Heating-powered one-zone kilonova: (M_ej, v_ej) -> L(t), T_eff(t), R_ph(t)."""

    def __init__(self, m_msun, v_ej_c, kappa=KAPPA_SRC, v_ph_frac=V_PH_FRAC,
                 t_min=1e-2, t_max_d=30.0, n=3000, t_scale=1.0, t_scale_gas=False):
        """`t_scale` perturbs the launch spectrum's temperature (t_core) at fixed
        L and R_ph -- the diffusion solution and the photosphere are untouched,
        only what illuminates the shell changes (§4.43's T-direction check);
        `t_scale_gas` scales T_gas with it."""
        self.m_msun, self.v_ej_c, self.kappa, self.v_ph_frac = m_msun, v_ej_c, kappa, v_ph_frac
        self.t_scale, self.t_scale_gas = float(t_scale), bool(t_scale_gas)
        self.tau_d = float(diffusion_time(m_msun, v_ej_c, kappa))
        _, self.fth_clamped = barnes_params(m_msun, v_ej_c)
        self.t = time_grid(t_min, t_max_d, n)
        self.qdot = deposited_power(self.t, m_msun, v_ej_c)
        self.L = arnett_luminosity(self.t, self.qdot, self.tau_d)

    # -- scalar accessors (log-log interpolation on the stored grid) ----------
    def _interp(self, y, t):
        return np.exp(np.interp(np.log(t), np.log(self.t), np.log(np.maximum(y, 1e-300))))

    def luminosity(self, t):
        return self._interp(self.L, t)

    def qdot_at(self, t):
        return self._interp(self.qdot, t)

    def f_th(self, t):
        return thermalization(t, self.m_msun, self.v_ej_c)

    def v_ph(self, t):
        """v_ph / v_ej at time t (array-safe)."""
        t = np.asarray(t, float)
        if self.v_ph_frac != "grey":
            return np.full_like(t, float(self.v_ph_frac))
        r_out = self.v_ej_c * C * t
        f = 1.0 - TAU_PH / (self.kappa * self.rho(t) * r_out)
        return np.clip(f, V_PH_MIN, None)

    def r_ph(self, t):
        return self.v_ph(t) * self.v_ej_c * C * np.asarray(t, float)

    def tau_grey(self, t):
        """Grey optical depth of the whole sphere, kappa rho R_out."""
        return self.kappa * self.rho(t) * self.v_ej_c * C * np.asarray(t, float)

    def t_eff(self, t):
        return (self.luminosity(t) / (4.0 * np.pi * SIGMA_SB * self.r_ph(t) ** 2)) ** 0.25

    @property
    def t_peak(self):
        return float(self.t[np.argmax(self.L)])

    def rho(self, t):
        """Uniform-sphere density inside v_ej at time t."""
        r = self.v_ej_c * C * np.asarray(t, float)
        return self.m_msun * MSUN / ((4.0 / 3.0) * np.pi * r ** 3)

    def state(self, t_d):
        """The `trajectory.state` dict at epoch t_d (days), plus source diagnostics."""
        t = t_d * DAY
        T_grey = float(self.t_eff(t))
        T = T_grey * self.t_scale
        return dict(t_exp=t, rho=float(self.rho(t)),
                    T_gas=T if self.t_scale_gas else T_grey, t_core=T,
                    core_law="source", r_core=float(self.r_ph(t)),
                    r_out=self.v_ej_c * C * t,
                    L=float(self.luminosity(t)), Qdot=float(self.qdot_at(t)),
                    f_th=float(self.f_th(t)), tau_d=self.tau_d, T_eff=T, T_eff_grey=T_grey,
                    t_scale=self.t_scale, t_scale_gas=self.t_scale_gas,
                    R_ph=float(self.r_ph(t)), fth_clamped=self.fth_clamped,
                    v_ph=float(self.v_ph(t)), tau_grey=float(self.tau_grey(t)),
                    v_ph_floored=bool(self.v_ph_frac == "grey"
                                      and 1.0 - TAU_PH / self.tau_grey(t) < V_PH_MIN),
                    kappa=self.kappa, v_ph_frac=self.v_ph_frac,
                    m_ej_msun=self.m_msun, v_ej_c=self.v_ej_c)


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
    for m in (0.003, 0.01, 0.03):
        for v, ls in zip((0.05, 0.1, 0.2), ("-", "--", ":")):
            s = SourceModel(m, v); td = s.t / DAY; ok = td > 0.05
            lab = f"M={m}, v={v}"
            ax[0].loglog(td[ok], s.L[ok], ls, label=lab)
            ax[1].semilogx(td[ok], s.t_eff(s.t[ok]), ls)
            ax[3].semilogx(td[ok], s.v_ph(s.t[ok]), ls)
            ax[2].semilogx(td[ok], s.f_th(s.t[ok]), ls)
    ax[0].set(xlabel="t (d)", ylabel="L (erg/s)", ylim=(1e39, 3e42), title="Arnett L(t)")
    ax[1].set(xlabel="t (d)", ylabel="T_eff (K)", ylim=(0, 15000), title="T_eff = [L/4πσR_ph²]^¼")
    ax[2].set(xlabel="t (d)", ylabel="f_th", title="Barnes+2016 eq. 35")
    ax[3].set(xlabel="t (d)", ylabel="v_ph / v_ej", ylim=(0, 1.05), title="grey photosphere, τ = 2/3")
    ax[0].legend(fontsize=6, ncol=1); fig.tight_layout(); fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None, help="write the L/T_eff/f_th diagnostic figure here")
    a = ap.parse_args()
    for m in (0.003, 0.01, 0.03):
        for v in (0.05, 0.1, 0.2):
            s = SourceModel(m, v)
            st1, st3 = s.state(1.0), s.state(3.0)
            print(f"M={m:5.3f} v={v:4.2f}  tau_d={s.tau_d/DAY:5.2f} d  t_peak={s.t_peak/DAY:5.2f} d  "
                  f"L(1d)={st1['L']:.2e}  T_eff(1d)={st1['T_eff']:6.0f} K  T_eff(3d)={st3['T_eff']:6.0f} K  "
                  f"v_ph/v_ej(1d,7d)={st1['v_ph']:.3f},{s.state(7.0)['v_ph']:.3f}  "
                  f"f_th(1d)={st1['f_th']:.2f}{' (clamped)' if s.fth_clamped else ''}")
    if a.plot:
        _plot(a.plot)
