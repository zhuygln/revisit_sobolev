"""Scalar statistics of a line forest, and of the redistribution it produces.

These are the candidate axes for the phase diagram that asks *when* a grouped
opacity closure fails (F30/F31). Five of them existed only as copy-pasted
one-liners at five call sites; two did not exist at all. Collected here so a
forest -- synthetic or GSI or CHIANTI -- is characterized the same way
everywhere.

Ray-resolved S and E already exist as `sobolev.sobolev_leg.crossing_depths`;
this module deliberately does not reimplement them. What it adds is the
forest-level (not ray-level) summary, plus the redistribution-range measure
that F31 identified as the thing separating the two regimes.
"""
import numpy as np

from .constants import C
from .atomic_data import nearest_neighbour_velocity_spacing, overlap_parameter


def saturation_stats(atom):
    """Opacity-line counts and the saturation census.

    Replaces the inline `(t > 1).sum()` idiom repeated across paper2/paper3.
    """
    t = np.asarray(atom.op_tau, float)
    if t.size == 0:
        return dict(n_opacity=0, tau_max=0.0, n_tau_gt1=0, n_tau_gt01=0,
                    frac_saturated=0.0)
    return dict(n_opacity=int(t.size),
                tau_max=float(t.max()),
                n_tau_gt1=int((t > 1).sum()),
                n_tau_gt01=int((t > 0.1).sum()),
                frac_saturated=float((t > 1).mean()))


def SE_sums(atom):
    """Forest-integrated S = sum tau and E = sum (1 - e^-tau), and their ratio.

    F15's two quantities. E is the expected interaction count a photon crossing
    every resonance would register; S is what the Bernoulli survival
    exponentiates. E/S -> 1 for a weak forest and falls as lines saturate, so
    it is a dimensionless saturation coordinate with no free scale.
    """
    t = np.asarray(atom.op_tau, float)
    S = float(t.sum())
    E = float(-np.expm1(-t).sum())
    return dict(S=S, E=E, E_over_S=(E / S if S > 0 else np.nan),
                deficit=S - E)


def crowding(atom, tau_cut=1.0):
    """Saturated lines per unit ln(lambda) -- the density axis, N_sat.

    Counting per decade of wavelength rather than per Angstrom makes the number
    scale-free, so a synthetic forest at 4000 A and one at 12000 A with the same
    crowding land on the same point.
    """
    nu = np.asarray(atom.op_nu, float)
    t = np.asarray(atom.op_tau, float)
    if nu.size < 2:
        return dict(n_sat_per_lnlam=0.0, n_per_lnlam=0.0, span_lnlam=0.0)
    span = float(np.log(nu.max() / nu.min()))
    if span <= 0:
        return dict(n_sat_per_lnlam=np.inf, n_per_lnlam=np.inf, span_lnlam=0.0)
    return dict(n_sat_per_lnlam=float((t > tau_cut).sum()) / span,
                n_per_lnlam=float(t.size) / span,
                span_lnlam=span)


def spacing_stats(atom, v_doppler=None):
    """Nearest-neighbour velocity spacing of the opacity lines, and (optionally)
    the overlap parameter O = v_D / dv against a given Doppler width.

    Note `run_mc` is a point-interaction code and carries no Doppler width, so
    `v_doppler` is only meaningful for the deterministic legs; the spacing
    itself is what the grouped transport actually sees.
    """
    lam = np.sort(C / np.asarray(atom.op_nu, float) * 1e8)
    if lam.size < 2:
        return dict(dv_median=np.nan, dv_min=np.nan, overlap_median=np.nan)
    dv = nearest_neighbour_velocity_spacing(lam)
    out = dict(dv_median=float(np.median(dv)), dv_min=float(dv.min()))
    if v_doppler is not None:
        out["overlap_median"] = float(np.median(overlap_parameter(dv, v_doppler)))
    return out


def redistribution_range(nu_in, nu_out, w=None, edges=None):
    """How far in frequency a redistribution event moves a photon.

    `|d ln lambda| = |ln(nu_in / nu_out)|` per event, from the (absorbed,
    exited) pairs `run_mc(collect_events=True)` returns and
    `paper3/phase0_reference/reference_events_*.npz` stores.

    `same_group_frac` is the quantity F31's two-regime reading points at: if a
    photon re-emerges inside the group it was absorbed in, the group closure has
    no information left to place it, and the ordering of resonances within the
    group is exactly what was thrown away. Requires `edges`.
    """
    nu_in = np.asarray(nu_in, float); nu_out = np.asarray(nu_out, float)
    w = np.ones(nu_in.size) if w is None else np.asarray(w, float)
    good = np.isfinite(nu_in) & np.isfinite(nu_out) & (nu_in > 0) & (nu_out > 0)
    if not good.any():
        return dict(n_events=0)
    a, b, ww = nu_in[good], nu_out[good], w[good]
    d = np.abs(np.log(a / b))
    order = np.argsort(d)
    cw = np.cumsum(ww[order]); cw /= cw[-1]

    def q(f):
        return float(d[order][np.searchsorted(cw, f)])

    out = dict(n_events=int(d.size),
               mean_abs_dlnlam=float(np.average(d, weights=ww)),
               median_abs_dlnlam=q(0.5),
               iqr_abs_dlnlam=q(0.75) - q(0.25),
               # signed: negative = net blueward (fluorescent upscatter)
               mean_dlnlam=float(np.average(np.log(a / b), weights=ww)))
    if edges is not None:
        e = np.asarray(edges, float); ng = e.size - 1
        gi = np.clip(np.searchsorted(e, a, side="right") - 1, 0, ng - 1)
        gj = np.clip(np.searchsorted(e, b, side="right") - 1, 0, ng - 1)
        out["same_group_frac"] = float(np.average(gi == gj, weights=ww))
        out["mean_abs_dgroup"] = float(np.average(np.abs(gj - gi), weights=ww))
    return out


def forest_summary(atom, v_doppler=None, events=None, edges=None):
    """Every scalar for one forest, in one dict -- the phase-diagram row."""
    out = {}
    out.update(saturation_stats(atom))
    out.update(SE_sums(atom))
    out.update(crowding(atom))
    out.update(spacing_stats(atom, v_doppler))
    if events is not None:
        nu_in, nu_out = events[0], events[1]
        w = events[2] if len(events) > 2 else None
        out.update(redistribution_range(nu_in, nu_out, w, edges))
    return out
