# Paper B — the effective redistribution operator: plan

The PI's statement of 2026-09-21 is verbatim in [`plan_review.md`](plan_review.md);
the preregistered gate is [`prl_gate.md`](prl_gate.md). This file is the map.

## The question

> Can detailed energy-conserving lanthanide fluorescence be replaced by a
> compact effective redistribution operator R_ij, while holding the rest of
> radiative transfer — especially resolved Sobolev opacity — fixed?

The decisive comparison, on identical resolved Sobolev opacity and
identical physical states, is

    ε*   versus   R_ij   versus   the full energy-conserving macroatom.

The neural surrogate R_φ(j | i, θ) stays in the plan but only matters if a
fixed or interpolated R_ij cannot carry the state dependence (question 3).

## What is already known, and what is not

- Paper III (F25, F27, F29, F32; `docs/results_report.md` §4.23–4.29) found
  that a group-to-group operator with discrete within-group exit tables
  reproduces explicit branching on the same opacity: La II to 1.6 %, Ce II
  to 7.5 %, Nd II to 0.4 % in the worst band at N_g = 4; that the operator
  is local in frequency, not low-rank; and that an opacity-weighted mixing
  rule predicts a blend at the 5 % level. Those kernel legs use the
  resolved crossing loop and are **not** affected by the bin-inversion bug
  38ebf30 — but they are **pre-F50**: photon-number packets, the
  `sobolev_branch` reference, τ_max = 5 normalised states. F50 moved the
  reference by 1–2.5 mag once the bookkeeping conserved energy. So
  question 1 is a genuine redo, not a confirmation.
- Paper IV's `A2` leg is the energy-packet version at N_g = 32 on the full
  13-ion P1 mixture: within 0.05 mag of the energy-conserving macroatom
  (F61). It is a prior for question 1 at the largest N_g, not an answer at
  small N_g or per ion.
- Paper III's ε findings (F22–F24) are INVALID (38ebf30); ε*'s irreducible
  error is measured fresh in G1.
- **There is no nuclear-heating or composition-realisation ensemble in this
  repository.** The Paper III grid varied (M, v, X_lan) of a source model,
  not nuclear inputs. Question 4's σ_nuclear can only be proxied here by
  the composition patterns that exist (Gillanders Ye-0.21a at X_lan 0.11
  and at its own 0.30, the solar r-process residual at 0.11, and the
  Ye-0.29a pattern); G4 says so and uses that proxy, or waits for an
  external ensemble.

## The four questions, in order, with their kill logic

| gate | question | runs | continue if | stop if |
|---|---|---|---|---|
| **G1** | does the compression survive correct energy physics? | La II, Ce II, Nd II × [ε grid, R_2, R_4, R_8, R_16, R_32, macroatom] on P1 2 d shell 28 | R_4 or R_8 within 0.1 mag of the macroatom for at least two ions and ε* substantially worse (B1, B2 Green) | R needs more than 32 groups, works for one ion only, or ε* is nearly as good (B1 or B2 Red) |
| **G2** | why: frequency locality or parameter count? | the N_g*-group matrix against a rank-k factorisation of the 128-group matrix with the same parameter count, transported | locality wins by more than the noise | rank does as well (a compression result, not a physics one) |
| **G3** | how universal: state, species, mixtures | P1 at 1/3/5 d and a T_gas sweep with fixed vs recomputed kernels; `R_mix = Σ_a α_a R^(a)` against a La+Ce+Nd macroatom run; the 13-ion P1 blend at N_g = 2…32 | the mixing rule predicts the blend without a fit and the state dependence is tabulable | the kernel must be refit per state or per mixture |
| **G4** | negligible against astrophysical uncertainty? | σ_R from R_{N_g*} − macroatom on each composition realisation; σ_nuclear across the realisations | σ_R ≪ σ_nuclear | σ_R comparable to σ_nuclear |

G1 is the first kill gate and is fully preregistered in `prl_gate.md`.
G2–G4 are named there and are detailed, each as its own preregistered
section, only after the previous gate's reading. The surrogate enters only
if G3 finds the state dependence not tabulable.

## The production experiment G1, in one paragraph

The state is `paper4/phase1_benchmarks/P1_t2.json` at its photospheric
shell 28 (T = 3401 K, ρ = 3.30×10⁻¹⁵ g cm⁻³, 2 d, worldline transport),
with one ion carrying the whole opacity at the number density the P1
composition gives it: La II 3.29×10⁴, Ce II 9.54×10⁴, Nd II 1.37×10⁵ cm⁻³.
Every leg uses indivisible energy packets, the same zone, core and launch,
3 seeds × 3×10⁵ packets. The reference is the downward energy-conserving
macroatom (`sobolev_dmacro`, Paper IV's R2). The operator legs are
`sobolev_group` with a kernel built from the reference's own events at
N_g = 2, 4, 8, 16, 32 log-spaced groups with discrete within-group exit
tables and energy rows; a 128-group kernel is built for the event-level
metric but not transported. The scalar legs are `sobolev_tla` at
ε = 0, 0.05, …, 1, and ε* is the grid value minimising the mean |Δm| over
the live bands against the reference. Metrics, readings and gray
conditions are in `prl_gate.md`; the central figure is transport error
against effective-model complexity for the three ions with ε* as a line.

## Provenance rules (inherited from Paper IV)

Every gate has a pre-declared reading with a gray outcome, checked
gray-first. Every number a manuscript quotes is generated from a frozen
record (`paperB/freeze.py` → `FROZEN.json`, when the manuscript starts).
The claim read-through with numbers substituted runs before submission.
Educational results (`education/`) are never cited as evidence.

## Layout

    paperB/plan.md, prl_gate.md, plan_review.md   this plan, the preregistration, the PI's words
    paperB/gate1/run_gate1.py                     the G1 runner (one ion per call)
    paperB/gate1/analyse.py, figure.py            metrics, readings, the central figure
    paperB/gate1/gate1_<ion>.json                 the records; gate1_verdict.json
    tests/test_paperB_gate1.py                    the verdict logic, eps*, m_event, a smoke run
