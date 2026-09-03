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

---

## Addendum 2026-09-01 — F36–F41

The audit above was written before the last six findings existed and its table
header still reads "Paper III (F25–F35)". Classifying the rest, on the same
three classes:

| # | class | note |
|---|---|---|
| F36 | **the re-run itself** | not a finding to classify but the instrument that reclassified F24/F27/F30/F33 above |
| F37 | **invariant as a negative** | "the synthetic model has only one side of the boundary" is a statement about that model, and its diagnosis (no net inflow to the band without delocalized exits) was confirmed by the fix that followed |
| F38 | **structure invariant, magnitudes conditional** | "the two approximations are not additive" is a within-ion comparison at fixed opacity and survives; the interaction term's *size* (−6.4% Pr, −4.3% Ce) is one state, one band, one grid. The **sign-flip on Pr II is one ion** and should not be generalized |
| F39 | **existence and orientation invariant; location conditional** | already stated that way in §4.34b — the crossing moves 179 → 688 with exit_tau, so only the sign structure is claimed. **Reproducibility debt**: `boundary.py`'s `main` does not drive the exit_tau scan, so the committed script regenerates only F37 |
| F40 | **conditional** | one density history, one composition, LTE, imposed T(t), and ρ(1 d) chosen so the crossing lands inside 0.5–8 d. The invariant part is structural: a homologous history sweeps S over orders of magnitude and must cross any boundary inside that range. The *epoch* is not a prediction |
| F41 | **invariant in structure, conditional in magnitude** | "the error is chromatic, not bolometric, and the narrow-band residual is a poor proxy for either" is a within-ion statement at fixed opacity, reproduced on two ions, a blend, two core laws and two velocity regimes. The *sizes* (0.74 mag in g−r at 0.5 d, >1 mag at kilonova velocities) belong to their stated states |

### What §10 changes in the earlier classifications

**F30's structural claim is now measured in observer units.** Across every ion,
epoch and velocity regime run, the grouped-*redistribution* leg moves no band by
more than ~0.01 mag. "The opacity binds, not the redistribution" stops being a
statement about band ratios and becomes: *the kernel compression of F25/F27 is
free to an observer.* That strengthens F27 and F32 as well — the compressible
half of the hierarchy is compressible in the quantity people actually measure.

**F40's diagnostic band is demoted.** §4.36's whole argument runs through the
3800–3955 Å residual, and §4.37 shows that residual is a poor proxy for the
photometric error in *both* directions: −59.5% in band with +0.007 mag
bolometric (La II, binned, 0.75 d), and +55% in band with −0.74 mag in g−r
(Ce II, 0.5 d). F40's conclusion survives — near-zero residual still carries no
information about correctness — but the *quantity* in which the boundary was
located should not be read as an observable.

### Standing reproducibility debts

Not classifications, but the same kind of problem the audit exists to catch:

1. **F39's exit_tau scan has no driver and no data file.** `boundary.json` holds
   only the F37 negative result.
2. **The Ce II density scan** — the first of the four independent locations of
   the boundary — exists only as six numbers in `docs/lab_notebook.md`.
3. **The F38 A/B/C table has no generating script**; it is a manual join of
   `phase9_audit/audit.json` with `phase8_survey/survey.json`.
4. **F40's crossing (1.17 d, S = 47.5) was interpolated by hand.** Fixed in §10:
   `observables.py:crossing_epoch` computes and stores it.

### Addendum (2026-09-02): debts 1–3 closed

- Debt 1: `paper3/synthetic/boundary.py --delocalize 1 --n-exit 6 --exit-tau 0.5,2.0`
  → `boundary_exit_tau.json`; crossings S = 179.6 / 689.0 reproduce §4.34b, and
  the dlnlam that was never recorded turns out to be inert at delocalize = 1.
- Debt 2: `paper3/phase8_survey/density_scan.py` → `density_scan_58CeII.json`
  (three seeds; crossing S = 55.7; five of six points within 7 points of the
  single-seed notebook numbers, the steepest one within 9).
- Debt 3: `paper3/phase9_audit/counterfactual_table.py` → `counterfactual_table.json`
  reproduces the §4.35 A/B/C table from `survey.json` alone.
- Debt 4 was closed in §10 (`crossing_epoch`).
- `tests/test_provenance.py` pins all three. Classification of F42 (§4.38):
  **a robustness check, not a new claim** — the F41 chromatic effect is a
  property of the spectra, not of the top-hat filters.

### Addendum (2026-09-02, later): F43 and F44

- **F43 (§4.40)**: a *new measurement* on a new source model — the closure's
  colour error at physically normalized densities. Driver `grid.py` /
  `run_grid.py`, data `grid/model_*.json` (27), tables `grid_table.py`,
  Fig. 2. Provenance complete. Classification: **new claim, bounded** — the
  statement is about colours within the launch window under a conserving
  core normalization; the same JSONs record why no bolometric claim is made
  (`dm_bol_absorbing`, `f_return`, `f_dep`).
- **F44 (§4.40)**: the *Gate 2 verdict*, a pre-registered classification
  (thresholds in the plan before the grid ran; `low_N` added as a flag
  after a partial run, changing no class). Driver `sensitivity.py`, data
  `sensitivity.json`, Fig. 4, tests `tests/test_sensitivity.py`.
  Classification: **decision-grade for the action plan, upper bound as
  physics** — C-B against three one-zone parameters does not bound what a
  richer fit would absorb.
- Deviations from the plan, all recorded in §4.39–§4.40: grey photosphere
  instead of v_ej/2; chain thermalization, wall clock and the conserving
  normalization added to the harness (defaults preserved); the pre-probe
  removed; `low_N`.

### Addendum (2026-09-02, later still): F45 and F46, and an erratum against F44

- **Erratum (F44 / §4.40)**: the floor mask was a no-op (`grid.py` never
  stored `v_ph_floored`); 27 of 153 rows were floored and included. Gate 2
  re-baselined with the mask real: C-B at 20 of 21 analysable points (was
  24 of 24), no class changed at any point analysable under both rules.
  F44's row is restated; the §4.40 numbers stand with the erratum pointer.
  Classification of the error: **a bookkeeping error in the harness, not
  in the physics**; caught by reading the code, not by a test, which is
  why `tests/test_grid_harness.py` and `tests/test_sensitivity.py` now pin
  the flag.
- **F45 (§4.41)**: robustness of F43/F44 to three harness conventions plus
  one layer of source freedom. Drivers `sensitivity.py --floored/--core/--tangent`,
  `robustness.py chain|table`, data `sensitivity_{floored_incl,absorbing,T1,T2,T3}.json`,
  `robustness/chain_*.json`, Fig. 6, tests in `tests/test_sensitivity.py`
  (nuisance absorption, dof accounting, absorbing re-derivation) and
  `tests/test_grid_harness.py`. Classification: **a robustness result with
  an honest converse** — the pre-declared class flips to C-A at most points
  under a free luminosity history, and the report says so next to the
  amplitudes that flip it (1–7 mag per epoch) and the misfit that remains
  (χ²_res/dof ≈ 23). Two diagnostics (`R_nuisance_only`, `a_over_dln`)
  were added *after* seeing the T1 result; they change no class and are
  labelled as post hoc.
- **F46 (§4.42)**: Phase 3A, three pre-declared scenarios and a
  pre-declared Gate 3. Driver `paper3/phase13_observability/observe.py`,
  data `observability.json`, Fig. 5, tests `tests/test_observability.py`.
  Classification: **decision-grade** — the residual is detected everywhere
  and survives (M, v, X) everywhere; under T1 its survival is a measured
  function of NIR coverage (9/18, 3/8, 1/15). The expectation "NIR is the
  leverage" was written before the run and is confirmed in a sharper form
  (detection is optical, distinctness is NIR).
- **F47 (§4.43)**: the T_eff validation, pre-declared (cosine threshold
  0.8 written before the runs). Driver `paper3/phase12_grid/tscale.py`,
  data `tscale.json` + `grid/tscale/*.json`, Fig. 7. Classification:
  **decision-grade for the proxy, and one unplanned physics result** — the
  illumination temperature does not reach the observer at the grid's
  saturation (norm 0.06 of the proxy, within the noise floor), which was
  not predicted and is reported as such; the gas-temperature response has
  the proxy's shape (0.92) and 1.35× its size. Two harness bugs in
  `tscale.py` (row-level `T_gas`, partial JSONs) were found and fixed
  before any number was quoted. One grid point only.
- **F48 (§4.44)**: the grid completed; every downstream number regenerated
  and the §4.40–§4.42 counts superseded with pointers, not overwritten.
  Classification: **bookkeeping with one honest failure** — the chain-cap
  criterion, as pre-declared, fails for the sub-0.3-mag colours at the
  worst-trapped cells and passes for J−K and for the class; the report
  prints the failed half. §4.40's "166 of 170" NIR count could not be
  reproduced under any current rule (it was counted with the dead mask);
  §4.44 states the rule it uses. The chain-cap table is complete (four
  cells × three caps, `robustness/chain_*.json`): the magnitude half of the
  pre-declared criterion fails (4 of 12) and is reported as failed, with the
  non-monotone cap dependence as the reason it does not indicate a
  systematic; the class half passes at all 27 points.
