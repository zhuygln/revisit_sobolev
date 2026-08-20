"""Analytic limits of the Sobolev attenuation leg.

These pin the leg that carries the "Sobolev approximation proper" claim, so
they check the formula rather than merely exercising the code path.
"""

import numpy as np

from sobolev.constants import C
from sobolev.optical_depth import tau_sobolev
from sobolev.sobolev_leg import expansion_damp, sobolev_attenuation

T_EXP = 864000.0  # 10 days
R_CORE = 1.0e14
R_OUT = 5.0e14
LAMBDA0 = 4000e-8
NU0 = C / LAMBDA0
F_OSC = 0.5
N0 = 2.0


def nu_for_plane(z_res, nu0=NU0):
    """Observer frequency whose resonance plane sits at z_res."""
    return nu0 / (1.0 - z_res / (C * T_EXP))


def test_single_line_fully_shadowing_plane():
    """A plane beyond the whole core (z_res > r_core) shadows every core ray,
    so the transmitted fraction is exactly exp(-tau_S)."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    assert 0.3 < tau < 3.0
    nu = np.array([nu_for_plane(2.0e14)])
    att = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    assert np.isclose(att[0], np.exp(-tau), rtol=1e-12)


def test_off_resonance_is_untouched():
    nu = np.array([NU0 * 0.90, NU0 * 1.30])
    att = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    assert np.allclose(att, 1.0)


def test_two_separated_lines_multiply():
    """Both planes in front of the core: attenuations multiply."""
    nu_b = NU0 * 1.02
    z = 2.0e14
    # Choose a frequency whose plane for line A is at z; line B's plane for
    # the same frequency sits elsewhere, so probe each line separately and
    # then a frequency crossing both.
    lines = [(NU0, F_OSC), (nu_b, F_OSC)]
    tau_a = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    tau_b = tau_sobolev(F_OSC, N0, C / nu_b, T_EXP)

    nu_a_only = np.array([nu_for_plane(z, NU0)])
    att_a = sobolev_attenuation(nu_a_only, [lines[0]], R_CORE, R_OUT, T_EXP, N0)
    att_both = sobolev_attenuation(nu_a_only, lines, R_CORE, R_OUT, T_EXP, N0)
    # Line B's plane for this frequency is far outside the shell, so adding
    # it must not change the result.
    assert np.isclose(att_a[0], np.exp(-tau_a), rtol=1e-12)
    assert np.isclose(att_both[0], att_a[0], rtol=1e-12)

    # A frequency whose plane for B lies in front of the core.
    nu_b_probe = np.array([nu_for_plane(z, nu_b)])
    att_b = sobolev_attenuation(nu_b_probe, [lines[1]], R_CORE, R_OUT, T_EXP, N0)
    assert np.isclose(att_b[0], np.exp(-tau_b), rtol=1e-12)


def test_expansion_damp_reproduces_the_cap():
    """With the expansion damping the single-line result is exp(-(1-e^-tau))."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    nu = np.array([nu_for_plane(2.0e14)])
    att = sobolev_attenuation(
        nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, damp=expansion_damp
    )
    assert np.isclose(att[0], np.exp(-(1.0 - np.exp(-tau))), rtol=1e-12)
    # The cap always transmits more light than true Sobolev for tau > 0.
    assert att[0] > np.exp(-tau)


def test_pop_frac_scales_tau():
    """pop_frac multiplies the density, hence tau, linearly."""
    nu = np.array([nu_for_plane(2.0e14)])
    full = sobolev_attenuation(nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    half = sobolev_attenuation(
        nu, [(NU0, F_OSC, 0.5)], R_CORE, R_OUT, T_EXP, N0
    )
    assert np.isclose(half[0], np.sqrt(full[0]), rtol=1e-12)


def test_partial_shadowing_lies_between():
    """A plane inside the core radius shadows only the outer rays, so the
    transmitted fraction sits strictly between exp(-tau) and 1."""
    tau = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    nu = np.array([nu_for_plane(0.5 * R_CORE)])
    att = sobolev_attenuation(
        nu, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, n_p=400
    )
    assert np.exp(-tau) < att[0] < 1.0


# --------------------------------------------------- relativistic modes
#
# The default (relativity=None) path carries every number in Paper I, so the
# first test here is a bit-for-bit guard on it, not an approximate one.

CT = C * T_EXP
_Z_PROBE = np.array([0.5, 1.5, 2.5, 3.5, 4.5]) * 1e14
NU_PROBE = NU0 / (1.0 - _Z_PROBE / CT)

# Hard-coded outputs of the classical path, captured before the relativistic
# branch was added. array_equal, not allclose: the point is that adding modes
# did not perturb the arithmetic by a single ulp.
_SOB_REF = np.array([
    0.8488409252739633, 0.3996263539825767, 0.3996263539825767,
    0.3996263539825767, 0.3996263539825767,
])
_EXP_REF = np.array([
    0.886350430017603, 0.5486066131172791, 0.5486066131172791,
    0.5486066131172791, 0.5486066131172791,
])


def test_default_path_is_bit_identical():
    sob = sobolev_attenuation(NU_PROBE, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0)
    exp_ = sobolev_attenuation(
        NU_PROBE, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, damp=expansion_damp
    )
    # guard against the probes silently missing the shell, which would make
    # any comparison here pass trivially
    assert (sob < 0.999).all() and (exp_ < 0.999).all()
    assert np.array_equal(sob, _SOB_REF)
    assert np.array_equal(exp_, _EXP_REF)


def test_relativistic_modes_reduce_to_classical_as_beta_goes_to_zero():
    """Both corrections are O(beta), so the deviation must fall like beta."""
    deviations = []
    for b_out in (1.0e-3, 1.0e-4, 1.0e-5):
        r_core, r_out = 1.0e-8 * CT, b_out * CT
        nu = np.array([NU0 / (1.0 - 0.5 * r_out / CT)])
        kw = dict(r_core=r_core, r_out=r_out, t_exp=T_EXP, n_ref=N0)
        base = sobolev_attenuation(nu, [(NU0, F_OSC)], **kw)[0]
        wl = sobolev_attenuation(nu, [(NU0, F_OSC)], relativity="worldline", **kw)[0]
        ex = sobolev_attenuation(nu, [(NU0, F_OSC)], relativity="exact", **kw)[0]
        assert base < 0.999
        deviations.append(max(abs(wl - base), abs(ex - base)) / base)

    assert deviations[-1] < 5.0e-5
    for coarse, fine in zip(deviations, deviations[1:]):
        assert np.isclose(coarse / fine, 10.0, rtol=0.05)  # linear in beta


def test_worldline_leg_gives_one_minus_beta_squared_over_gamma():
    """The physical law, at the leg level rather than the tau level.

    Core made tiny so every ray is radial; then a fully shadowing line must
    transmit exp(-tau_S(t_exp) (1-beta)^2/gamma) at the resonance frequency.
    """
    r_core, r_out = 1.0e-4 * CT, 0.30 * CT
    tau_s = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    for beta in (0.05, 0.10, 0.20, 0.25):
        lorentz = 1.0 / np.sqrt(1.0 - beta**2)
        nu = np.array([NU0 * np.sqrt((1.0 + beta) / (1.0 - beta))])
        att = sobolev_attenuation(
            nu, [(NU0, F_OSC)], r_core, r_out, T_EXP, N0, relativity="worldline"
        )[0]
        assert np.isclose(att, np.exp(-tau_s * (1.0 - beta) ** 2 / lorentz), rtol=2e-4)


def test_exact_leg_gives_one_minus_beta_over_gamma():
    """The frozen snapshot, same geometry -- and a far bigger deficit."""
    r_core, r_out = 1.0e-4 * CT, 0.30 * CT
    tau_s = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP)
    for beta in (0.05, 0.10, 0.20, 0.25):
        lorentz = 1.0 / np.sqrt(1.0 - beta**2)
        nu = np.array([NU0 * np.sqrt((1.0 + beta) / (1.0 - beta))])
        att = sobolev_attenuation(
            nu, [(NU0, F_OSC)], r_core, r_out, T_EXP, N0, relativity="exact"
        )[0]
        assert np.isclose(att, np.exp(-tau_s * (1.0 - beta) / lorentz), rtol=2e-4)


def test_worldline_locus_matches_a_numerical_root_find():
    """z_res is a closed form; check the implementation against brentq.

    D(z) = Z0/sqrt(Z0^2 + 2 Z0 z - p^2), so D = nu0/nu is linear in z. This
    pins the code path, not just the algebra.
    """
    from scipy.optimize import brentq

    for a in (0.0, 0.05, 0.2):
        for y in (1.05, 1.2, 1.4):
            big_z, p_ray = CT, a * CT
            closed = 0.5 * big_z * (y**2 - 1.0) + p_ray**2 / (2.0 * big_z)
            f = lambda z: big_z / np.sqrt(big_z**2 + 2 * big_z * z - p_ray**2) - 1.0 / y
            assert np.isclose(brentq(f, 0.0, 20.0 * CT, xtol=1e-8), closed, rtol=1e-10)


def test_exact_second_root_never_falls_inside_the_shell():
    """Why the plane picture was adequate for Paper I.

    The frozen locus is a two-root quadratic -- Jeffery's CD/CP surfaces. The
    upper root sits above beta_z = 1 - a^2, so for any shell reaching only
    0.35c it is always outside. Carried in the code, never fired in practice.
    """
    b_out = 0.35
    r_core, r_out = 1.0e-3 * CT, b_out * CT
    p_ray = np.linspace(1e-6, r_core, 40)[:, None]
    nu = np.linspace(NU0 * 1.001, NU0 * 1.9, 60)[None, :]

    a = p_ray / CT
    x = NU0 / nu
    disc = x**2 * (1.0 - a**2) - a**2
    u_plus = (1.0 + x * np.sqrt(np.maximum(disc, 0.0))) / (1.0 + x**2)

    real = disc > 0.0
    assert real.any()  # the branch is genuinely exercised
    assert (u_plus[real] > b_out).all()
    z_hi = np.sqrt(np.maximum(r_out**2 - p_ray**2, 0.0))
    assert (u_plus[real] * CT > np.broadcast_to(z_hi, u_plus.shape)[real]).all()


def test_damp_is_applied_after_the_relativistic_tau():
    """The expansion leg must stay a same-code differential at high beta too."""
    r_core, r_out = 1.0e-4 * CT, 0.30 * CT
    beta = 0.2
    lorentz = 1.0 / np.sqrt(1.0 - beta**2)
    tau_rel = tau_sobolev(F_OSC, N0, LAMBDA0, T_EXP) * (1.0 - beta) ** 2 / lorentz
    nu = np.array([NU0 * np.sqrt((1.0 + beta) / (1.0 - beta))])
    att = sobolev_attenuation(
        nu, [(NU0, F_OSC)], r_core, r_out, T_EXP, N0,
        damp=expansion_damp, relativity="worldline",
    )[0]
    assert np.isclose(att, np.exp(-(1.0 - np.exp(-tau_rel))), rtol=2e-4)


def test_n_of_beta_shapes_the_initial_profile():
    """n_of_beta is read as the profile at t_exp, then diluted."""
    r_core, r_out = 1.0e-4 * CT, 0.30 * CT
    beta = 0.2
    nu = np.array([NU0 * np.sqrt((1.0 + beta) / (1.0 - beta))])
    kw = dict(r_core=r_core, r_out=r_out, t_exp=T_EXP, n_ref=N0,
              relativity="worldline")
    flat = sobolev_attenuation(nu, [(NU0, F_OSC)], **kw)[0]
    halved = sobolev_attenuation(
        nu, [(NU0, F_OSC)], n_of_beta=lambda b: 0.5 * np.ones_like(b), **kw
    )[0]
    assert np.isclose(halved, np.sqrt(flat), rtol=1e-9)  # tau halves


def test_relativity_mode_is_validated():
    import pytest

    with pytest.raises(ValueError, match="relativity must be"):
        sobolev_attenuation(
            NU_PROBE, [(NU0, F_OSC)], R_CORE, R_OUT, T_EXP, N0, relativity="worldine"
        )
