# Paper IV — the campaign directory

This directory holds the **code, run outputs and record** of Paper IV. The
reusable physics it introduces lives in the package (`sobolev/`), not here:
`paper4/` is drivers, outputs and the frozen record, `sobolev/` is physics,
and `paper2/phase1/forest_mc.py` stays the propagation core with thin hooks.
The write-up half, `docs/paper4/`, does not exist yet.

The program as received is [`plan.md`](plan.md) (2026-09-09, verbatim); the
PI's review of the implementation plan drafted from it, whose corrections
are all adopted below, is [`plan_review.md`](plan_review.md) (same day,
verbatim). Paper III's directory, `paper3/`, is frozen and untouched by this
campaign except for the Phase 0 relabel recorded in its README.

## The question

Paper III measured a 1–3 mag chromatic error of the grouped-opacity closure on
a controlled fixed-atmosphere grid (F43–F49). That harness conserves photon
number, not energy: every fluorescence step books the level-energy
difference into a deposit that is never re-radiated, and each leg is then
rescaled grey to a common luminosity. The near-term question is therefore

> Does grouped line opacity still alter the spectrum when both the resolved
> and grouped calculations use the same energy-conserving downward atomic
> transport?  **(B₂ − R₂)**

and, only if that survives, after radiative equilibrium,

> Does that difference survive when each transport treatment establishes
> its own self-consistent thermal and ionization structure?  **(B_eq − R_eq)**

## Where it stands (2026-09-09)

Phases 0–3 ran the day the plan arrived (`docs/results_report.md` §4.47,
F50–F52; notebook §9ax–§9ba). Gate 2 on P1, the xkn secular component,
is **Green**: B₂ − R₂ = −2.2 to −2.9 mag in z and +0.9 to +1.3 mag in K at
1–5 d, the redistribution control ≤ 0.06 mag, C₂ within 0.08 mag of B₂.
On P2, the published AT2017gfo 3.4-d state, it is **Red** — +0.19 mag in g
with the opposite sign — which is Paper III's saturation boundary
re-found on published states, not a contradiction; the combined verdict
is Gray by the letter of the rule. The single-zone adequacy trigger fired
on P1, so shell-resolved transport is a required check. Photon-number
branching moved the *reference* by 2.5 mag in z (F50). Phases 4–7 proceed.

## Design decisions taken in implementation

(Recorded here because the plan and its review leave them to the
implementation; each names what it omits.)

- **Nomenclature.** R = resolved Sobolev opacity, B = grouped opacity with
  the *same* detailed atomic transport, A = resolved opacity with the
  compressed redistribution kernel, C = grouped opacity with the kernel.
  Subscript 1 = photon-number branching (Paper III), E = energy packets with
  photon-number probabilities (the bookkeeping-only rung), 2 = energy
  packets with the downward macroatom. `expansion_dmacro` is **B₂** from its
  first run; there is no "C₂ (Phase 2)".
- **Energy packets** conserve the *comoving* energy `E_cm = w·hν_cm` at every
  atomic interaction; the lab energy differs by the local Doppler factor and
  is conserved in free flight. Because the transport loop never reads `w`,
  the E rung has bit-identical histories to the photon rung.
- **The Phase 2 atom is a downward macroatom** (`dmacro`): radiative
  deactivation ∝ A β (ε_u − ε_l), internal downward jumps ∝ A β ε_l, no
  internal upward transitions (they need a radiation-field estimator, which
  arrives with WP8–10). Exits use the closed-form A β tables and never
  re-enter the β chain; the populated-level cascade test pins the single
  application of β against the A β² prediction a second pass would give.
  It is energy-conserving downward fluorescence, not a full Lucy macroatom,
  and is named accordingly everywhere.
- **Normalisation.** The geometric-series scale `core="equilibrium"` in
  `sobolev/photometry.py` is exact for i.i.d. relaunches, so it is the
  energy-conserving normalisation; it is wrapped in
  `sobolev/energy_balance.py`, not edited. The re-emitting core exists as a
  validation mode only (it must agree with the scale to < 1 % per band) and
  never runs in production.
- **Single zone = the local shell at τ_grey ≈ 2/3.** Never a global
  mass-weighted average. The per-shell state is built and gated in full; one
  shell's own ρ, T, X and ion fractions feed transport, neighbours are the
  robustness check, and a per-epoch adequacy ratio (band saturation in the
  local zone against the shell-resolved value) pre-declares when shells
  become necessary: outside 0.5–2 in a live band.
- **Benchmarks are exact published models**, pinned numerically as tests:
  P1 the xkn comparison *secular* component (Ricigliano et al. 2024, MNRAS
  529, 647: M = 2.64×10⁻² M☉, v_rms = 0.06c, Y_e = 0.20, s = 10 k_B/baryon,
  τ_exp ≈ 17 ms, its own density prescription); P2 the Gillanders et al.
  2026 3.4-d model (MNRAS 548, stag748: v_in = 0.15c, v_out = 0.35c,
  ρ₀ = 4×10⁻¹⁵ g cm⁻³, ρ ∝ v⁻³, T = 3200 K, X_LN ≈ 2.5×10⁻³). Every value is
  ledgered in `data/README.md` before use; figures are a sanity check, not
  the acceptance criterion. Two inputs are **provisional** and marked so in
  the ledger and the state metadata until the PI confirms them at Gate 1:
  P1's lanthanide fraction (0.10; neither the xkn paper nor Lippuner &
  Roberts 2015 tabulates it) and both lanthanide patterns (the solar
  r-process residuals of Prantzos et al. 2020; the Ye−0.29a element list
  is in the 2022 paper's supplement).
- **The Morag-type treatment (D) is a dual-role closure**: the transport
  quantity (EP93 or Στ per bin) sets where encounters happen, the capped
  quantity `min[κ_l, 1/(ρct)]` limits net energy reprocessing at an
  encounter. It is not a bin cap on the interaction probability and is not
  labelled "Morag" until checked against the paper's equations.
- **The Fontes bridge** compares the resolved, expansion and *binned*
  thermal-redistribution legs before swapping redistribution for `dmacro`.
- **Provenance.** Tags are never moved: `paper3-freeze` is the original
  Paper III freeze, `paper3-freeze-x4ln` the label-only re-freeze of
  Phase 0. Golden SHA-256 histories of every `run_mc` mode
  (`tests/data/golden_run_mc.json`, `tests/test_golden_run_mc.py`) were
  taken before any transport edit and must stay green after each.

## Pre-declared gates

**Gate 1** (Phase 1): `EjectaState.check()` passes at every epoch — mass
integral to rel 10⁻⁶, Σ_Z X = 1 including the bulk, Σ_stage f = 1 — and the
committed P1/P2 states reproduce the pinned published inputs (rel 10⁻⁶ on
mass and v_rms for P1; exact for P2), as tests.

**Gate 2** (Phase 3), primary **B₂ − R₂**, controls A₂ − R₂ and C₂ − R₂,
live bands only:

| outcome | criterion | consequence |
|---|---|---|
| Green | \|Δm_B₂\| ≥ 0.5 mag in ≥ 2 live bands at ≥ 2 epochs; \|Δm_A₂\| ≤ 0.1 mag everywhere; g brighter, K fainter kept. If also C₂ within 0.1 mag of B₂: redistribution compresses, opacity is the culprit | Phases 4–7, then 8–10 |
| Yellow | 0.2 ≤ max\|Δm_B₂\| < 0.5 mag, controls in bounds | mechanism survives at a smaller amplitude; Phases 4–7; 8–10 deferred |
| Red | max\|Δm_B₂\| < 0.2 mag **and** ‖R₁ − R₂‖ > ½‖R₁ − B₁‖ | the Paper III result was the energy treatment; stop the astrophysical interpretation and write that up |
| Gray / unresolved | anything else: a control out of bounds, B₂ − R₂ small with R₁ − R₂ also small, a convergence failure (Phase 6), the adequacy trigger, or a sign pattern contradicting the causal reading | diagnose before interpreting; no narrative is forced |

**Gate 3** (Phase 7): the B-vs-R class and colour signs survive II/III Saha
ionization; the amplitude may move; a sensitivity table over the neutral-
stage policy and the bulk species accompanies the number.

Energy gates in every production row: identity residual < 10⁻¹²;
`E_dep_cm/E_inj = 0` for `dmacro` legs; the conserving/equilibrium scale
ratio equals the Doppler work W/E_esc (a physical O(v/c) loss, reported,
not re-radiated); re-emitting core vs the exact scale < 1 % per band once
per benchmark.

## Layout

Directory numbers are execution order; the last column maps to `plan.md`'s
phases. Directories appear when their phase starts.

```
plan.md                   the program as received, kept unedited         --
plan_review.md            the PI's review of the implementation plan       --
phase1_benchmarks/        build.py -> P1_t{1,2,3,5}.json, P2_t3.4.json,
                          gate1_*.png; Gate 1 in tests/test_benchmarks.py  Phase 1
phase2_energy/            frame tests, toy atoms, the R1/R1E/R2/B1/B1E/B2
                          ladder                                          Phases 2, 2A, 2B
phase3_legs/              verdict.py: the pre-declared Gate 2 on legs_*.json  Phase 3
phase4_thermal_bridge/    the Fontes-like thermal limit, incl. binned     Phase 4
phase5_dual_role/         the dual-role (Morag-type) closure D            Phase 5
phase6_convergence/       tau_min, bin width, packets, table cut          Phase 6
phase7_ionization/        LTE Saha II/III, Gate 3                         Phase 7
phase11_recovery/         global parameter recovery on the paper3 grid   Phase 11
external/tardis/          the external-code record                        Phase 12
```

Phase 0 (this README, the relabel, the golden hashes and the CSR sampler)
has no directory: its record is `paper3/README.md`, `docs/lab_notebook.md`
§9ax and the tag `paper3-freeze-x4ln`. Phases 8–10 (shells, estimators, the
full macroatom, radiative equilibrium, **B_eq − R_eq**) open only on a Green
or Yellow Gate 2.

Tests live in the repo-wide `tests/`.
