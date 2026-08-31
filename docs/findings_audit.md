# Findings audit — what survives, what is conditional, what is superseded

Written 2026-08-30, after F35 showed the project's normalization recipe is
ion-specific by accident. Every cross-ion claim was measured under it, so each
finding has to be re-classified before any of them go into a flagship
manuscript's main text.

## The three classes

**Invariant** — survives a change of normalization standard, because it is a
statement about a mechanism, an identity, or a single fixed physical state.

**Conditional** — true as measured, at the state it was measured at. Quotable
with its conditions attached; not quotable as a general property of the ion or
of line transport.

**Superseded** — a later result replaces it. Kept in the record with a pointer,
never repeated as current.

## Why the reclassification was needed

`setup.py` normalized every atom to τ_max = 5 **inside 3850–3950 Å**. That is
ion-specific by accident: it works for La/Ce/Nd because their strongest lines
happen to lie near that window, and diverges for ions where they do not
(Yb II demands n_ion = 1.7×10¹² cm⁻³, τ_max = 1.7×10⁸, β = 1.5×10⁻⁸, at which
point the branch chain cannot terminate).

The size of the effect on a headline number: **Ce II's band-3800 grouped-closure
error is +126.7% under the window recipe and +12.2% under a global one**, from a
25% density difference (§4.32). Anything comparing ions is suspect until re-run.

`sobolev/normalization.py` now provides two standards which must not be mixed:
`global_tau_max` for controlled cross-ion comparison, `from_conditions` for
astrophysical claims.

---

## Paper I (F1–F18)

Largely **invariant**. These are mechanism and identity results measured against
deterministic references at stated conditions, not cross-ion comparisons.

| # | class | note |
|---|---|---|
| F1, F4, F5, F12 | invariant | structural statements about the approximations |
| F15 | **invariant** | an exact identity, pinned to 1e-13; the S vs E decomposition underlies everything later |
| F2, F8, F16–F18 | conditional | measured at stated conditions, cross-code conventions matched |
| F3 | superseded | by F11 |
| F6, F9, F10, F13 | conditional | quoted with v_D, τ_max and geometry attached, as they already are |
| F11 | invariant | frozen vs worldline is a statement about which problem is being solved |
| F14 | invariant | a statement about what the reference codes contain |

**Not affected by the normalization issue** — Paper I works within single
forests at declared conditions and never claims cross-ion universality.

## Paper II (F19–F24)

| # | class | note |
|---|---|---|
| F19 | **invariant** | escape-probability branching is exact and tested |
| F20 | conditional | La II at one state; the refill number is state-specific |
| F21 | conditional | sign reversal at the stated ε and state |
| F22 | **likely invariant within the tested state** | "no scalar ε reproduces La II fluorescence" is a statement about the TLA's reach, not about density; the ε_best *values* are conditional |
| F23 | conditional | verdict robust to 0.1c, calibration explicitly not transferable — already stated that way |
| F24 | **needs re-run** | the ion-dependence of ε_best and the "+21% La / +113% Ce" density limit are exactly the cross-ion comparison the window recipe distorts |

## Paper III (F25–F35)

| # | class | note |
|---|---|---|
| F25 | **valid** | compression at each ion's own reference state; not a cross-ion claim |
| F26 | conditional | La II only; F28 already narrows it |
| F27 | **promising, needs confirmation** | "compression is generic" is a cross-ion claim and used the window recipe |
| F28 | conditional | the T_src ion-dependence is a cross-ion claim; the τ_scale collapse is structural and likely invariant |
| F29 | conditional | one blend, one ratio, one state |
| F30 | **magnitudes conditional, structure invariant** | "the opacity is the binding constraint, not the redistribution" survives — it is a within-ion comparison at fixed opacity. The *catastrophic Ce magnitude* does not |
| F31 | conditional | superseded in part by F33 (its second clause already retracted) |
| F32 | **invariant** | rank/locality is a property of the operator, measured on stored kernels; a density change moves the kernel but not the argument that coarsening beats truncation at matched parameter count |
| F33 | conditional | the memory-depth null result must be reconfirmed on the new grid |
| F34 | **superseded** | the power law is replaced by F35's sign change |
| F35 | **provisional, potentially central** | the phase boundary; needs the controlled 2-D diagram to become mechanistic |

---

## What this means for the manuscripts

**Main text should contain invariant claims only.** On the current
classification that is: F15's S/E identity, F11's transport-convention
distinction, F19's escape-probability branching, F22's TLA-reach result, F30's
structural claim that the opacity binds rather than the redistribution, F32's
locality-not-low-rank result, and — if it survives the controlled diagram —
F35's sign change.

**Conditional results belong in the results sections with their state
attached**, which is largely how they are already written.

**Two claims need re-running before they can be quoted across ions**: F24's
density limit and F27's Gate 2. `paper3/phase9_audit/audit.py` does exactly
this for La/Ce/Nd/Pr/Yb under the controlled standard.

## Result of the re-run (§4.33, F36)

Done, five ions, controlled standard. Three classifications change:

| finding | provisional class | **after re-run** |
|---|---|---|
| F27 Gate 2 | promising, needs confirmation | **invariant, and stronger** — every ion compresses at four groups to ≤4.3%; "dense ions need 32–64 groups" was the artefact |
| F24 density limit | needs re-run | **superseded** — +113% on Ce II becomes +1.9%; Ce is the *better* case at matched line strength. Ion-dependence survives; "dense forests are where it fails" does not |
| F30 opacity binds | structure invariant, magnitudes conditional | **confirmed on five ions** — redistribution 0.2–1.8% vs grouped opacity −31% to +15% |
| F33 memory null | conditional | **superseded** — memory is the most effective correction found (Ce II +12.2 → +0.2%); the null belonged to an over-dense Ce |
| F35 boundary | provisional | **provisional, and now located three ways** at S ≈ 50 (density scan, 13-ion survey, synthetic sweep) |

The lesson generalizes past this project: a normalization chosen to make one
atom convenient became load-bearing for claims about all atoms, and inverted one
of them. Cross-ion claims need a standard that cannot depend on which ion was
picked first.
