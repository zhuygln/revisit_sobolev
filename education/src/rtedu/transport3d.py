"""The miniature end-to-end transport of chapter 13: packets in a
homologously expanding sphere, in three dimensions, with time.

Along a straight ray from position x in direction n the comoving frequency
of a packet of lab frequency nu_lab is, to first order,

    nu_com(s) = nu_lab (1 - (x . n + s) / (c t)),

so it falls linearly with the path length s whatever the direction: the
Sobolev sweep of chapter 4 in three dimensions. A packet meets line l at
s_l = c t (1 - nu_l / nu_lab) - x . n if that point is ahead of it and
inside the sphere r_out = v_max t. At an interaction the redistribution
model names the exit line; the packet is re-emitted isotropically (this
is where trapping comes from) with the exit line's comoving frequency at
its new position, and its clock advances by s / c. The state (the line
depths, the temperature, the sphere) is evaluated at the packet's own time,
frozen for each flight: rho ∝ t^-3 so tau ∝ t^-2, and T(t) is prescribed.

Sources: an initial population of packets uniform in the volume at t0
with the thermal emissivity of T(t0) (the stored field), and heating
packets injected uniformly in the volume at times drawn from a t^-1.3
heating law over [t0, t1]. Every packet carries the same energy. The
light curve is the escaped energy per unit escape time; the bands are the
toy bands of rtedu.bands. Nothing here is a kilonova model: the gas has no
energy equation, the sources are prescribed, and the Doppler shift is first
order. It is the same loop as the production light curves, on ten lines.
"""
import numpy as np

from . import C, DAY


class ToyEjecta:
    """rho ∝ t^-3 on a uniform sphere of edge velocity v_max_c; T(t) = T0 (t/t0)^-alpha."""

    def __init__(self, atom, v_max_c=0.2, n0=30.0, t0=2.0 * DAY, T0=4000.0, alpha_T=0.5, heat_index=1.3):
        self.atom = atom; self.v_max = v_max_c * C; self.n0 = n0; self.t0 = t0; self.T0 = T0
        self.alpha_T = alpha_T; self.heat_index = heat_index

    def T(self, t):
        return self.T0 * (t / self.t0) ** (-self.alpha_T)

    def n_total(self, t):
        return self.n0 * (t / self.t0) ** (-3.0)

    def r_out(self, t):
        return self.v_max * t

    def tau(self, t):
        return self.atom.line_list(self.T(t), self.n_total(t), t)

    def emis(self, t):
        return self.atom.thermal_emissivity(self.T(t), self.n_total(t))


def _isotropic(rng, n):
    mu = rng.uniform(-1.0, 1.0, n); phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - mu * mu)
    return np.stack([s * np.cos(phi), s * np.sin(phi), mu], axis=-1)


def launch(rng, ej, n_init, n_heat, t1):
    """The initial population at t0 and the heating packets over [t0, t1]:
    positions uniform in the sphere at their own epoch, isotropic,
    frequencies from the thermal emissivity at their epoch's temperature."""
    t = np.concatenate([np.full(n_init, ej.t0), _heating_times(rng, ej, n_heat, t1)])
    n = t.size
    u = rng.random(n) ** (1.0 / 3.0)
    x = _isotropic(rng, n) * (u * ej.r_out(t))[:, None]
    dirn = _isotropic(rng, n)
    nu = np.empty(n)
    for i in range(n):
        e = ej.emis(t[i]); k = int(np.searchsorted(np.cumsum(e), rng.random()))
        k = min(k, ej.atom.n_lines - 1)
        nu[i] = ej.atom.nu[k] / (1.0 - np.dot(x[i], dirn[i]) / (C * t[i]))     # comoving line frequency seen in the lab
    return x, dirn, nu, t


def _heating_times(rng, ej, n, t1):
    """Times drawn ∝ t^-heat_index on [t0, t1] by inversion."""
    if n == 0:
        return np.zeros(0)
    a = 1.0 - ej.heat_index
    u = rng.random(n)
    return (ej.t0 ** a + u * (t1 ** a - ej.t0 ** a)) ** (1.0 / a)


def transport(rng, ej, x, dirn, nu, t, redistribution_at, max_interactions=400):
    """Every packet to escape. `redistribution_at(t) -> model(rng, k)` gives
    the redistribution model at the packet's time. Returns escape times,
    escape lab frequencies, interaction counts."""
    n = t.size
    t_esc = np.empty(n); nu_esc = np.empty(n); n_int = np.zeros(n, int)
    for i in range(n):
        xi = x[i].copy(); ni = dirn[i].copy(); nui = float(nu[i]); ti = float(t[i])
        for _ in range(max_interactions):
            tau = ej.tau(ti); r_out = ej.r_out(ti)
            xn = float(np.dot(xi, ni))
            # path length to the sphere's edge along the ray
            s_edge = -xn + np.sqrt(max(xn * xn - (np.dot(xi, xi) - r_out * r_out), 0.0))
            s_l = C * ti * (1.0 - ej.atom.nu / nui) - xn
            ok = (s_l > 1e-6) & (s_l < s_edge)
            hit = None
            for k in np.flatnonzero(ok)[np.argsort(s_l[ok])]:
                if rng.random() < 1.0 - np.exp(-tau[k]):
                    hit = int(k); break
            if hit is None:
                t_esc[i] = ti + s_edge / C; nu_esc[i] = nui
                break
            xi = xi + s_l[hit] * ni; ti = ti + s_l[hit] / C; n_int[i] += 1
            j = redistribution_at(ti)(rng, hit)
            ni = _isotropic(rng, 1)[0]
            nui = ej.atom.nu[j] / (1.0 - float(np.dot(xi, ni)) / (C * ti))
        else:
            t_esc[i] = np.nan; nu_esc[i] = nui                 # capped: reported, not hidden
    return t_esc, nu_esc, n_int


def light_curve(t_esc, edges, energy=None):
    """Escaped energy per unit time in the given time bins."""
    energy = np.ones(t_esc.size) if energy is None else energy
    ok = np.isfinite(t_esc)
    h, _ = np.histogram(t_esc[ok], bins=edges, weights=energy[ok])
    return h / np.diff(edges)


def trace_one(rng, ej, x, dirn, nu, t, redistribution_at, max_interactions=400):
    """One packet with its history: the list of (position, time, line index
    or -1 for the launch/escape) at launch, at each interaction and at
    escape (for the chapter-13 animation)."""
    xi = np.array(x, float); ni = np.array(dirn, float); nui = float(nu); ti = float(t)
    trace = [(xi.copy(), ti, -1)]
    for _ in range(max_interactions):
        tau = ej.tau(ti); r_out = ej.r_out(ti); xn = float(np.dot(xi, ni))
        s_edge = -xn + np.sqrt(max(xn * xn - (np.dot(xi, xi) - r_out * r_out), 0.0))
        s_l = C * ti * (1.0 - ej.atom.nu / nui) - xn
        ok = (s_l > 1e-6) & (s_l < s_edge); hit = None
        for k in np.flatnonzero(ok)[np.argsort(s_l[ok])]:
            if rng.random() < 1.0 - np.exp(-tau[k]):
                hit = int(k); break
        if hit is None:
            trace.append((xi + s_edge * ni, ti + s_edge / C, -1)); return trace
        xi = xi + s_l[hit] * ni; ti = ti + s_l[hit] / C
        j = redistribution_at(ti)(rng, hit); trace.append((xi.copy(), ti, int(j)))
        ni = _isotropic(rng, 1)[0]; nui = ej.atom.nu[j] / (1.0 - float(np.dot(xi, ni)) / (C * ti))
    return trace
