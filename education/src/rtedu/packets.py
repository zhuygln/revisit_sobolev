"""Monte Carlo packets in a slab: the optical-depth draw and the random walk.

The one idea of chapter 2: the probability that a packet travels an optical
depth greater than tau without interacting is e^{-tau}, so the interaction
depth of one packet is a draw from the exponential distribution,
tau = -ln(xi) with xi uniform on (0, 1]. Everything else is bookkeeping.
"""
import numpy as np


def sample_tau(rng, n):
    """n interaction optical depths, tau = -ln(xi); rng.random() is in [0, 1),
    so 1 - xi is in (0, 1] and the logarithm is finite."""
    return -np.log(1.0 - rng.random(n))


def propagate_slab(rng, tau_total, n, albedo=0.0, mu0=1.0, max_steps=10_000):
    """n packets through a plane-parallel slab of total optical depth
    tau_total (measured along the normal), launched at the illuminated face
    with direction cosine mu0. At each interaction a packet is scattered
    elastically with probability `albedo` (new direction isotropic, the
    frequency is NOT changed: chapter 3 stays below redistribution) and
    absorbed otherwise. Returns a dict with the transmitted, reflected and
    absorbed fractions, the mean number of interactions of the transmitted
    packets, and their mean path length in units of the mean free path.

    With albedo = 0 and mu0 = 1 the transmitted fraction is the Beer-Lambert
    law e^{-tau_total}: a packet is transmitted iff its single draw exceeds
    the slab depth.
    """
    x = np.zeros(n)                     # depth along the normal, in optical-depth units
    mu = np.full(n, float(mu0))
    alive = np.ones(n, bool)
    n_int = np.zeros(n, int)
    path = np.zeros(n)
    fate = np.zeros(n, int)             # 1 transmitted, 2 reflected, 3 absorbed
    for _ in range(max_steps):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        s = sample_tau(rng, idx.size)               # path to the next interaction
        x_new = x[idx] + mu[idx] * s
        path[idx] += s
        out = x_new >= tau_total
        back = x_new < 0.0
        # transmitted / reflected packets travel only to the face
        path[idx[out]] -= (x_new[out] - tau_total) / mu[idx[out]]
        path[idx[back]] -= x_new[back] / mu[idx[back]]
        fate[idx[out]] = 1; fate[idx[back]] = 2
        alive[idx[out | back]] = False
        stay = idx[~(out | back)]
        x[stay] = x_new[~(out | back)]
        n_int[stay] += 1
        scat = rng.random(stay.size) < albedo
        mu[stay[scat]] = rng.uniform(-1.0, 1.0, int(scat.sum()))
        fate[stay[~scat]] = 3; alive[stay[~scat]] = False
    tr = fate == 1
    return dict(transmitted=float(tr.mean()), reflected=float((fate == 2).mean()),
                absorbed=float((fate == 3).mean()), still_alive=int(alive.sum()),
                mean_interactions_transmitted=float(n_int[tr].mean()) if tr.any() else float("nan"),
                mean_path_transmitted=float(path[tr].mean()) if tr.any() else float("nan"),
                n=int(n))


def walk_2d(rng, radius, max_steps=500):
    """One packet's random walk in a disc of radius `radius` mean free paths:
    exponential step lengths, isotropic directions, from the centre until it
    leaves the disc. Returns the (x, y) vertices, the last one outside."""
    xy = [np.zeros(2)]
    for _ in range(max_steps):
        s = sample_tau(rng, 1)[0]
        phi = rng.uniform(0.0, 2.0 * np.pi)
        xy.append(xy[-1] + s * np.array([np.cos(phi), np.sin(phi)]))
        if np.hypot(*xy[-1]) >= radius:
            break
    return np.array(xy)
