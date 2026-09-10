"""Radially resolved ejecta states and the single-zone reduction (Paper IV WP1).

Paper III's transport had one zone: a uniform sphere, one density, one
temperature, an equal four-ion blend. Paper IV's benchmarks are coherent
published ejecta models with a radial density profile, a temperature run
and a full lanthanide pattern, so the state is carried per homologous shell:

    {t, v_edges, rho, T_gas, T_rad, X_Z, f_ion, n_e, Y_e}

with `check()` the plan's Gate 1 -- the mass integral, the composition sums
and the ionization-stage sums -- and `plot()` the sanity figure.

The transport is still single-zone (multi-shell run_mc is deferred until
Gate 2), and the PI fixed how a state collapses to one zone: NOT a
mass-weighted average, but the LOCAL state of the shell containing the grey
photosphere tau_grey = 2/3 (the line-forming region), with its neighbours
as a robustness check. `local_zone()` does exactly that and nothing else.
A per-epoch adequacy ratio -- band saturation in that zone against the
shell-resolved value -- is what pre-declares when shells become necessary.

Y_e, s and tau_exp are metadata: the benchmark's pattern table is keyed by
them, the state does not compute from them.
"""
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .constants import C
from .source import MSUN, SIGMA_SB

M_U = 1.66053906660e-24     # g
TAU_PH = 2.0 / 3.0


@dataclass
class EjectaState:
    t: float                              # s since merger
    v_edges: np.ndarray                   # cm/s, n_shell + 1, ascending
    rho: np.ndarray                       # g cm^-3 per shell
    T_gas: np.ndarray                     # K per shell
    T_rad: np.ndarray                     # K per shell (radiation temperature, if known)
    X: dict                               # element symbol -> mass fraction per shell ("bulk" allowed)
    f_ion: dict = field(default_factory=dict)   # "Ce II" -> fraction of Ce per shell
    n_e: np.ndarray | None = None         # cm^-3 per shell
    Y_e: float | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.v_edges = np.asarray(self.v_edges, float)
        self.rho = np.asarray(self.rho, float)
        self.T_gas = np.asarray(self.T_gas, float)
        self.T_rad = np.asarray(self.T_rad, float)
        self.X = {k: np.asarray(v, float) for k, v in self.X.items()}
        self.f_ion = {k: np.asarray(v, float) for k, v in self.f_ion.items()}
        n = self.n_shell
        for name, arr in (("rho", self.rho), ("T_gas", self.T_gas), ("T_rad", self.T_rad)):
            if arr.shape != (n,):
                raise ValueError(f"{name} has shape {arr.shape}, expected ({n},)")
        if np.any(np.diff(self.v_edges) <= 0):
            raise ValueError("v_edges must be strictly increasing")

    # ---- geometry ------------------------------------------------------
    @property
    def n_shell(self):
        return self.v_edges.size - 1

    @property
    def r_edges(self):
        return self.v_edges * self.t

    @property
    def v_mid(self):
        return 0.5 * (self.v_edges[1:] + self.v_edges[:-1])

    @property
    def r_mid(self):
        return self.v_mid * self.t

    def shell_volume(self):
        r = self.r_edges
        return 4.0 * np.pi / 3.0 * (r[1:] ** 3 - r[:-1] ** 3)

    def shell_mass(self):
        """Exact for piecewise-constant rho: 4 pi/3 rho (r_{i+1}^3 - r_i^3)."""
        return self.rho * self.shell_volume()

    def mass(self):
        return float(self.shell_mass().sum())

    # ---- gates ---------------------------------------------------------
    def check(self, m_target=None, tol=1e-6):
        """Gate 1: the mass integral (if a target is given), the composition
        sum and the ionization-stage sums, every shell. Returns a dict of
        the residuals; raises on failure."""
        out = {}
        if m_target is not None:
            out["mass_rel"] = abs(self.mass() - m_target) / m_target
            if out["mass_rel"] > tol:
                raise ValueError(f"mass integral {self.mass():.6e} g vs target {m_target:.6e} g "
                                 f"(rel {out['mass_rel']:.2e} > {tol})")
        xsum = sum(self.X.values())
        out["X_sum_max_dev"] = float(np.max(np.abs(xsum - 1.0)))
        if out["X_sum_max_dev"] > 1e-9:
            raise ValueError(f"mass fractions do not sum to 1 (max dev {out['X_sum_max_dev']:.2e})")
        by_elem = {}
        for ion in self.f_ion:
            el = ion.split()[0]
            by_elem[el] = by_elem.get(el, 0.0) + self.f_ion[ion]
        out["f_ion_max_dev"] = float(max((np.max(np.abs(v - 1.0)) for v in by_elem.values()), default=0.0))
        if out["f_ion_max_dev"] > 1e-9:
            raise ValueError(f"ion fractions do not sum to 1 (max dev {out['f_ion_max_dev']:.2e})")
        return out

    # ---- photosphere and the local zone --------------------------------
    def tau_grey(self, kappa):
        """Grey optical depth at the INNER edge of every shell, integrated
        from the outer boundary inward: tau_i = sum_{j >= i} kappa rho_j dr_j."""
        dr = np.diff(self.r_edges)
        return np.cumsum((kappa * self.rho * dr)[::-1])[::-1]

    def photospheric_shell(self, kappa, tau=TAU_PH):
        """Index of the shell containing the tau_grey = tau surface (the
        first shell, counted from outside, whose inner edge is deeper than
        tau); 0 with a flag if the whole ejecta is thinner than tau."""
        tg = self.tau_grey(kappa)
        deep = np.flatnonzero(tg >= tau)
        if deep.size == 0:
            return 0, False
        return int(deep.max()), True

    def n_ion(self, element, stage, shell, atomic_mass):
        """Number density of one ion in one shell: rho X f_ion / (A m_u)."""
        f = self.f_ion.get(f"{element} {stage}")
        f = 1.0 if f is None else f[shell]
        return float(self.rho[shell] * self.X[element][shell] * f / (atomic_mass * M_U))

    def local_zone(self, shell, t_core=None, v_out=None):
        """The transport contract for ONE shell's own state: `rho`, `T_gas`,
        `X` and `f_ion` of that shell; the core surface at the shell's inner
        edge; the outer boundary at the ejecta's edge (or `v_out`).

        This is the single-zone reduction the PI chose: the local
        line-forming state, never a mass-weighted average."""
        s = int(shell)
        r_out = (self.v_edges[-1] if v_out is None else v_out) * self.t
        return dict(t_exp=float(self.t), rho=float(self.rho[s]), T_gas=float(self.T_gas[s]),
                    t_core=float(self.T_rad[s] if t_core is None else t_core),
                    core_law="local_shell", r_core=float(self.r_edges[s]), r_out=float(r_out),
                    shell=s, v_core=float(self.v_edges[s]), v_out=float(r_out / self.t),
                    X={k: float(v[s]) for k, v in self.X.items()},
                    f_ion={k: float(v[s]) for k, v in self.f_ion.items()},
                    n_e=None if self.n_e is None else float(self.n_e[s]))

    def transport_zone(self, shells, t_core=None):
        """The multi-shell transport contract: the innermost transported
        shell's inner edge is the core, the outermost's outer edge the
        boundary; `shells` is a contiguous list of shell indices."""
        shells = list(shells)
        if shells != list(range(shells[0], shells[-1] + 1)):
            raise ValueError("transported shells must be contiguous")
        s0, s1 = shells[0], shells[-1]
        return dict(t_exp=float(self.t), r_core=float(self.r_edges[s0]), r_out=float(self.r_edges[s1 + 1]),
                    t_core=float(self.T_rad[s0] if t_core is None else t_core), core_law="local_shell",
                    shells=shells, v_edges=[float(v) for v in self.v_edges[s0:s1 + 2]],
                    rho=[float(x) for x in self.rho[s0:s1 + 1]], T_gas=[float(x) for x in self.T_gas[s0:s1 + 1]],
                    v_core=float(self.v_edges[s0]), v_out=float(self.v_edges[s1 + 1]))

    def regrid(self, v_edges_new, profile=None, T_of_v=None):
        """A new state on `v_edges_new`, which must contain every old edge it
        keeps (refine inside old shells, or coarsen across whole old shells).
        Refined sub-shells take rho from `profile` (the model's f(v), exact
        volume averages) when given, else the parent's constant rho; T, X and
        ion fractions are inherited (or T from `T_of_v`). Coarsened shells
        take mass-weighted rho, T and X and ion-mass-weighted f_ion, so the
        mass and the sums stay exact."""
        new = np.asarray(v_edges_new, float)
        old = self.v_edges
        if not (np.isclose(new[0], old[0]) and np.isclose(new[-1], old[-1])):
            raise ValueError("regrid keeps the inner and outer edges")
        n_new = new.size - 1
        parent = np.clip(np.searchsorted(old, 0.5 * (new[1:] + new[:-1]), side="right") - 1, 0, self.n_shell - 1)
        r_new = new * self.t
        vol_new = 4.0 * np.pi / 3.0 * (r_new[1:] ** 3 - r_new[:-1] ** 3)
        rho = np.empty(n_new); T_gas = np.empty(n_new); T_rad = np.empty(n_new)
        X = {k: np.empty(n_new) for k in self.X}; f_ion = {k: np.empty(n_new) for k in self.f_ion}
        n_e = None if self.n_e is None else np.empty(n_new)
        m_old = self.shell_mass()
        for i in range(n_new):
            inside = np.flatnonzero((old[:-1] >= new[i] - 1e-9 * new[i]) & (old[1:] <= new[i + 1] + 1e-9 * new[i + 1]))
            if inside.size >= 1 and np.isclose(old[inside[0]], new[i]) and np.isclose(old[inside[-1] + 1], new[i + 1]):
                w = m_old[inside]; tot = w.sum()
                rho[i] = tot / vol_new[i]
                T_gas[i] = np.sum(w * self.T_gas[inside]) / tot; T_rad[i] = np.sum(w * self.T_rad[inside]) / tot
                for k in X:
                    X[k][i] = np.sum(w * self.X[k][inside]) / tot
                for k in f_ion:
                    el = k.split()[0]
                    wm = w * self.X[el][inside]
                    f_ion[k][i] = np.sum(wm * self.f_ion[k][inside]) / wm.sum() if wm.sum() > 0 else self.f_ion[k][inside[0]]
                if n_e is not None:
                    n_e[i] = np.sum(w * self.n_e[inside]) / tot
            else:
                pj = parent[i]
                if profile is not None:
                    vv = np.linspace(new[i], new[i + 1], 401)
                    f_avg = np.trapezoid(profile(vv) * vv ** 2, vv) / np.trapezoid(vv ** 2, vv)
                    vo = np.linspace(old[pj], old[pj + 1], 401)
                    f_par = np.trapezoid(profile(vo) * vo ** 2, vo) / np.trapezoid(vo ** 2, vo)
                    rho[i] = self.rho[pj] * f_avg / f_par
                else:
                    rho[i] = self.rho[pj]
                T_gas[i] = self.T_gas[pj] if T_of_v is None else T_of_v(0.5 * (new[i] + new[i + 1]))
                T_rad[i] = self.T_rad[pj]
                for k in X:
                    X[k][i] = self.X[k][pj]
                for k in f_ion:
                    f_ion[k][i] = self.f_ion[k][pj]
                if n_e is not None:
                    n_e[i] = self.n_e[pj]
        out = EjectaState(t=self.t, v_edges=new, rho=rho, T_gas=T_gas, T_rad=T_rad, X=X, f_ion=f_ion,
                          n_e=n_e, Y_e=self.Y_e, meta=dict(self.meta, parent_grid=self.v_edges.tolist(),
                                                           parent_n_shell=self.n_shell))
        # profile-refined sub-shells reproduce the parent mass only up to the
        # quadrature; rescale each parent's children so the mass is exact
        if profile is not None:
            for pj in np.unique(parent):
                kids = np.flatnonzero(parent == pj)
                m_kids = (out.rho[kids] * vol_new[kids]).sum()
                if m_kids > 0:
                    out.rho[kids] *= m_old[pj] / m_kids
        return out

    # ---- io ------------------------------------------------------------
    def to_dict(self):
        return dict(t=self.t, v_edges=self.v_edges.tolist(), rho=self.rho.tolist(),
                    T_gas=self.T_gas.tolist(), T_rad=self.T_rad.tolist(),
                    X={k: v.tolist() for k, v in self.X.items()},
                    f_ion={k: v.tolist() for k, v in self.f_ion.items()},
                    n_e=None if self.n_e is None else self.n_e.tolist(),
                    Y_e=self.Y_e, meta=self.meta, mass_g=self.mass(), mass_msun=self.mass() / MSUN)

    def to_json(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=1) + "\n")

    @classmethod
    def from_dict(cls, d):
        return cls(t=d["t"], v_edges=d["v_edges"], rho=d["rho"], T_gas=d["T_gas"], T_rad=d["T_rad"],
                   X=d["X"], f_ion=d.get("f_ion", {}),
                   n_e=None if d.get("n_e") is None else np.asarray(d["n_e"], float),
                   Y_e=d.get("Y_e"), meta=d.get("meta", {}))

    @classmethod
    def from_json(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text()))

    # ---- figure --------------------------------------------------------
    def plot(self, out, kappa=None, atomic_mass=None):
        """rho(v), T(v), X_Z(v) and (if atomic masses are given) n_ion(v):
        the Gate 1 sanity figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        v = self.v_mid / C
        fig, ax = plt.subplots(1, 4, figsize=(17, 3.8))
        ax[0].semilogy(v, self.rho, "k.-"); ax[0].set(xlabel="v/c", ylabel=r"$\rho$ (g cm$^{-3}$)")
        if kappa is not None:
            s, ok = self.photospheric_shell(kappa)
            ax[0].axvline(self.v_edges[s] / C, color="r", ls="--", label=r"$\tau_{\rm grey}=2/3$" + ("" if ok else " (not reached)"))
            ax[0].legend(fontsize=7)
        ax[1].plot(v, self.T_gas, "k.-", label=r"$T_{\rm gas}$"); ax[1].plot(v, self.T_rad, "r.--", label=r"$T_{\rm rad}$")
        ax[1].set(xlabel="v/c", ylabel="T (K)"); ax[1].legend(fontsize=7)
        for el, x in self.X.items():
            if el != "bulk":
                ax[2].semilogy(v, np.maximum(x, 1e-30), ".-", label=el, lw=0.8, ms=3)
        ax[2].set(xlabel="v/c", ylabel="X_Z"); ax[2].legend(fontsize=5, ncol=3)
        if atomic_mass:
            for ion, f in self.f_ion.items():
                el = ion.split()[0]
                if el in atomic_mass and el in self.X:
                    n = self.rho * self.X[el] * f / (atomic_mass[el] * M_U)
                    ax[3].semilogy(v, np.maximum(n, 1e-30), ".-", label=ion, lw=0.8, ms=3)
            ax[3].legend(fontsize=5, ncol=3)
        ax[3].set(xlabel="v/c", ylabel=r"$n_{\rm ion}$ (cm$^{-3}$)")
        fig.suptitle(f"{self.meta.get('name', '')}  t = {self.t / 86400:.2f} d  M = {self.mass() / MSUN:.4e} Msun")
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)
        return out


# --------------------------------------------------------------------------
# profiles
# --------------------------------------------------------------------------

def shells_from_profile(m_ej, t, v_edges, f_of_v):
    """rho per shell for rho(v, t) = rho_0(t) f(v), rho_0 fixed by the mass
    integral over the shells (exact for the piecewise-constant discretisation:
    the shell density is the profile's volume average over the shell, so the
    total mass is exactly m_ej)."""
    v_edges = np.asarray(v_edges, float)
    r = v_edges * t
    vol = 4.0 * np.pi / 3.0 * (r[1:] ** 3 - r[:-1] ** 3)
    # volume-average f over each shell by fine quadrature in v
    f_avg = np.empty(v_edges.size - 1)
    for i in range(f_avg.size):
        vv = np.linspace(v_edges[i], v_edges[i + 1], 401)
        f_avg[i] = np.trapezoid(f_of_v(vv) * vv ** 2, vv) / np.trapezoid(vv ** 2, vv)
    rho0 = m_ej / np.sum(f_avg * vol)
    return rho0 * f_avg


def power_law_profile(m_ej, t, v_min, v_max, gamma, n_shell=32):
    """rho ~ v^gamma between v_min and v_max (Gillanders et al.'s TARDIS
    profiles have gamma = -3), log-spaced shells, normalised to m_ej."""
    v_edges = np.geomspace(v_min, v_max, n_shell + 1)
    return v_edges, shells_from_profile(m_ej, t, v_edges, lambda v: (v / v_max) ** gamma)


def xkn_profile(m_ej, t, v_max, n_shell=32, v_min_frac=0.0):
    """rho ~ (1 - (v/v_max)^2)^3 (Ricigliano et al. 2024, eq. 25), linear
    shells from v_min_frac v_max to v_max, normalised to m_ej."""
    v_edges = np.linspace(v_min_frac * v_max, v_max, n_shell + 1)
    return v_edges, shells_from_profile(m_ej, t, v_edges, lambda v: (1.0 - (v / v_max) ** 2) ** 3)


def xkn_vmax_from_vrms(v_rms):
    """v_max of the (1 - x^2)^3 profile with mass-weighted rms velocity
    v_rms: <v^2> = v_max^2 int x^4 (1-x^2)^3 / int x^2 (1-x^2)^3."""
    x = np.linspace(0.0, 1.0, 20001)
    num = np.trapezoid(x ** 4 * (1 - x ** 2) ** 3, x)
    den = np.trapezoid(x ** 2 * (1 - x ** 2) ** 3, x)
    return v_rms / np.sqrt(num / den)


def grey_temperature(state, kappa, luminosity):
    """The Eddington grey run T^4(v) = 3/4 T_eff^4 (tau_grey(v) + 2/3), with
    T_eff from the luminosity at the photospheric radius -- the same
    photosphere `SourceModel` uses. Returns (T per shell, shell index, T_eff)."""
    s, ok = state.photospheric_shell(kappa)
    r_ph = state.r_edges[s]
    t_eff = (luminosity / (4.0 * np.pi * SIGMA_SB * r_ph ** 2)) ** 0.25
    tg = state.tau_grey(kappa)
    tau_mid = 0.5 * (tg + np.concatenate([tg[1:], [0.0]]))
    return (0.75 * t_eff ** 4 * (tau_mid + TAU_PH)) ** 0.25, s, float(t_eff)
