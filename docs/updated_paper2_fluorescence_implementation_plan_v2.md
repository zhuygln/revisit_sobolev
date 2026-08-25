> **Execution status — updated 2026-08-25.** Every experiment this plan
> defines through Phase 4 has now been run; the gates are all decided.
> Numbers and provenance live in `docs/results_report.md` §4.19–§4.22 and
> the lab notebook §9p–§9s; this block maps the plan to what happened.
>
> | item | status | result |
> |---|---|---|
> | Paper I finalization (§1) | done, pushed | RE claim retracted, regression-guarded; stopping rule enforced (no Phase-1 numbers in Paper I) |
> | E1 energy accounting | done | identity closes to roundoff in every mode; level-energy check 2×10⁻⁵ (GSI wavelength consistency) |
> | E2 thermal width | done | SEDONA resolved RE 0.835 → 0.855 → 0.856 at v_D = 100/10/1 km/s (prediction stated in advance); MC within +0.7%; the MC's Δ_closure is v_D-free by construction |
> | E3 β regression suite | done | incl. the A·β exit-kernel identity |
> | E4/E5 ε sweep + χ² | done | **Outcome B**: ε_best 0.02–0.36 by band on the same opacity; optical 4500–6000 Å unreachable at any ε; χ²/dof ≥ 44 |
> | E6 redistribution matrix | done | branching = blue→optical channel, thermalisation = blue→red; ε interpolates and can never have one without the other |
> | E7 cascade pathways | done | pumps are 3300–4500 Å, not the far-UV (0.00% from 1142–2500 Å); 971 pathways, top 10 = 25% |
> | E13 bulk velocity 0.1c (§8) | done | outcome B survives worldline transport at 0.05–0.15c; ε_best shifts 0.15–0.27 (v_bulk = validity-map axis); frozen-first-order transport overstates blue transmission 3.8× — unusable at high β |
> | E8 line-dependent closure | done | +21% on La II — but **density-limited**: +113% on Ce II (see E9) |
> | E9 Ce II | done | **Outcome C**: ε_best ion-dependent (shifts 0.11–0.24; reachability flips in opposite directions); the Poisson opacity is 3 orders of magnitude too transparent in Ce II's blacked-out band |
> | E10 La+Ce mixture | done | blanketing neither suppresses (−44.6%) nor relocates the redistribution error; ε_best tracks the dominant forest |
> | E11 source-T sweep | open | — |
> | E12 shell-T sweep | open | — |
> | NLTE (§13) | deferred by design | — |
>
> **Gate A** passed; **Gate B** = "no scalar ε works" (outcome B), and E9
> adds outcome C; **Gate C** = qualitative verdict robust at 0.1c, calibration
> not transferable, frozen convention excluded. Paper II's identity per §14:
> *the SN-Ia scalar-closure result does not carry over to lanthanide
> networks — and the fix must carry line identity, within its density limit.*
> Next: E11/E12 if wanted, and the manuscript text
> (`docs/paper2/manuscript.tex` skeleton is structured and waiting).

# Updated Two-Paper Research Plan

## Current scientific split

### Paper I
**Question:** How accurately do resolved finite-profile transfer, per-line Sobolev transfer, and binned expansion opacity agree under controlled kilonova-like conditions?

Hierarchy:

\[
\boxed{
\text{resolved finite-profile transfer}
\rightarrow
\text{per-line Sobolev}
\rightarrow
\text{expansion opacity}
}
\]

Paper I should remain a controlled attenuation/transport-formalism paper.

### Paper II
**Question:** Can the scalar two-level-atom thermalisation closure that works reasonably well for iron-peak SN-Ia fluorescence also reproduce radiative redistribution through lanthanide atomic networks?

The decisive formulation is:

\[
\boxed{
\exists\ \epsilon\in[0,1]
\quad\text{such that}\quad
F_\lambda^{\rm expansion+TLA}(\epsilon)
\approx
F_\lambda^{\rm Sobolev+branching}\ ?
}
\]

The current Phase-1 result establishes only:

\[
\boxed{\epsilon=1\ \text{does not work for the La II reference problem}}
\]

It does **not yet** establish that no scalar \(\epsilon\) can work.

That distinction now defines the next experiments.

---


# Velocity hierarchy — three different velocity scales

The project now needs to keep three physically distinct velocity scales separate. They should **not** be treated as one velocity sweep.

| Scale | Meaning | Representative value | Main physics tested |
|---|---|---:|---|
| \(v_D\) | microscopic Doppler/profile width of an individual line | \(\sim0.5\)–\(1\ {\rm km\,s^{-1}}\) thermally; tests at 1, 10, 100 | Sobolev localization, finite-profile convergence, boundary clipping |
| \(\Delta v_{\rm shell}=v_{\rm out}-v_{\rm in}\) | macroscopic velocity span of the line-forming region | current controlled shell: \(2000\ {\rm km\,s^{-1}}\); realistic cases can be \(10^4\)–\(3\times10^4\ {\rm km\,s^{-1}}\) | frequency sweep through the line forest, number of resonances encountered, finite-region edge effects |
| \(v_{\rm bulk}\) | characteristic bulk ejecta speed relative to the observer | \(\sim0.1c\simeq3\times10^4\ {\rm km\,s^{-1}}\) | relativistic Doppler transport, frozen-snapshot vs worldline evolution, kilonova-scale external validity |

Schematic:

```text
microscopic line profile
        <--- v_D --->

macroscopic line-forming shell in velocity space
v_in  [========================================]  v_out
      <----------- Δv_shell ------------------>

the shell as a whole sits at an observer-frame expansion scale
v_bulk ~ 0.1c
```

A representative kilonova hierarchy is

\[
v_D\sim1\ {\rm km\,s^{-1}},
\qquad
\Delta v_{\rm shell}\sim10^4-3\times10^4\ {\rm km\,s^{-1}},
\qquad
v_{\rm bulk}\sim0.1c\sim3\times10^4\ {\rm km\,s^{-1}}.
\]

The last two scales need not have a strict ordering. They control different pieces of the transport problem.

Two useful dimensionless controls are

\[
\eta_{\rm width}
=
\frac{v_D}{\Delta v_{\rm shell}},
\]

which measures proximity to the zero-width/Sobolev limit and the importance of finite-region profile clipping, and

\[
\beta_{\rm bulk}
=
\frac{v_{\rm bulk}}{c},
\]

which controls relativistic and photon-worldline effects.

For realistic thermal lines,

\[
\eta_{\rm width}\sim10^{-5}-10^{-4},
\]

whereas kilonova ejecta have

\[
\beta_{\rm bulk}\sim0.1.
\]

Therefore

\[
\boxed{v_D=1,10,100\ {\rm km\,s^{-1}}}
\]

is a **microscopic line-width experiment**, whereas

\[
\boxed{v_{\rm bulk}\sim0.1c}
\]

is a **macroscopic transport/realism experiment**.

They answer different questions.

```text
1, 10, 100 km/s
      |
      v
Does the result survive the narrow-line / Sobolev limit?
      |
      v
epsilon sweep
      |
      v
Can one scalar TLA closure reproduce the branching spectrum?
      |
      v
v_bulk ~ 0.1c
      |
      v
Does that conclusion survive realistic kilonova expansion and
worldline-consistent transport?
```

The \(0.1c\) experiment is therefore **not a replacement** for the 1–10–100 km s\(^{-1}\) tests. It is an external-validity test performed after the closure itself is understood in the cleaner low-\(v/c\) harness.

---


# 1. Paper I finalization

## 1.1 Main result

At physically relevant thermal-like widths,

\[
\text{resolved}\approx\text{Sobolev proper},
\]

while expansion opacity introduces a much larger strong-line attenuation error.

The final interpretation should retain these distinctions:

- Sobolev localization and expansion opacity are separate approximations.
- Overlap is inert in the pure-absorption fixed-population harness because optical depths add exactly.
- Broad-\(v_D\) residuals are finite-region boundary effects.
- Refining line width does not remove the strong-line expansion-opacity bias.
- Paper I measures attenuation and controlled transport differences, not a complete kilonova spectrum.

## 1.2 Radiative-equilibrium correction

A previous narrow-grid radiative-equilibrium diagnostic was contaminated by the transport window. The earlier claim that redistribution reduced the expansion differential to roughly \(5\)–\(8\%\) should remain retracted.

With the whole ion's emissivity available, thermal re-emission is redistributed physically across wavelength rather than forced to remain inside the narrow transport window.

Required consistency checks:

- manuscript;
- response letter;
- README;
- findings table;
- `docs/results_report.md`;
- automated structure/grep checks.

Add the obsolete narrow-grid interpretation to the regression guards.

## 1.3 Stopping rule

Do **not** insert the new Phase-1 fluorescence numbers into Paper I.

Paper I ends at:

> attenuation and transport-formalism isolation.

Paper II starts when:

> line identity and wavelength redistribution become active physics.

---

# 2. Paper II revised scientific hypothesis

SN-Ia work showed that complete thermalisation, \(\epsilon=1\), can redistribute too much radiation redward, while an intermediate scalar thermalisation parameter such as \(\epsilon\sim0.3\) may approximate direct iron-peak fluorescence reasonably well.

Paper II asks whether that result transfers to open-\(f\)-shell lanthanide networks.

## Outcome A — SN-Ia result carries over

There exists one global value

\[
\epsilon_{\rm best}
\]

that reproduces the direct branching spectrum across UV, optical, and NIR.

Conclusion:

> Complete thermalisation is inappropriate for La II, but a calibrated scalar TLA closure remains effective.

## Outcome B — scalar closure fails

Different wavelength regions require different values:

\[
\epsilon_{\rm best}^{\rm UV}
\neq
\epsilon_{\rm best}^{\rm optical}
\neq
\epsilon_{\rm best}^{\rm NIR}.
\]

Conclusion:

\[
\boxed{
\text{The SN-Ia scalar-closure result does not carry over to La II.}
}
\]

## Outcome C — composition dependence

For example,

\[
\epsilon_{\rm best}^{\rm La}
\neq
\epsilon_{\rm best}^{\rm Ce}.
\]

Conclusion:

> A universal scalar thermalisation parameter is not transferable across lanthanide composition.

---

# 3. Current Paper-II instrument

## 3.1 Phase-1 implementation

Current vectorized same-code Monte Carlo instrument:

- full La II branching table: 17,743 transitions;
- 949 active opacity lines at 3000 K;
- five treatments in one code;
- 20 acceptance tests;
- pure-absorption legs reproduce Paper-I analytic results to <1%;
- thermal legs reproduce the SEDONA radiative-equilibrium spectrum sub-band by sub-band.

## 3.2 Current Phase-1 measurement

Reference setup:

- La II;
- fixed LTE populations at 3000 K;
- 6000 K Planck incident continuum;
- band: 3800–3955 Å;
- 3 seeds;
- \(2\times10^6\) packets per seed.

| Treatment | Emergent / Incident |
|---|---:|
| Sobolev, pure absorption | 0.183 |
| Expansion opacity, pure absorption | 0.344 |
| Sobolev + thermal re-emission | 0.257 |
| Expansion opacity + thermal | 0.412 |
| Sobolev + A-branching fluorescence | \(0.660\pm0.003\) |

Pure-absorption difference:

\[
\frac{0.344}{0.183}-1\approx+88\%.
\]

Redistributing comparison:

\[
\frac{0.412}{0.660}-1\approx-37.6\%.
\]

Candidate interpretation:

> Expansion opacity is too bright in pure absorption but, when paired with complete thermal redistribution, becomes too faint in the optical relative to direct La II fluorescence.

This sign reversal is **not yet the final headline** because \(\epsilon=1\) is only one TLA choice.

---

# 4. Escape probability: required physics

Radiative branching must include Sobolev escape probability:

\[
\beta_{ul}
=
\frac{1-e^{-\tau_{ul}}}{\tau_{ul}}.
\]

For strong lines,

\[
\tau\gg1
\quad\Rightarrow\quad
\beta\sim1/\tau.
\]

The effective competition among escaping radiative branches is therefore governed by approximately

\[
A_{ul}\beta_{ul},
\]

not just \(A_{ul}\).

This correction should become a named Paper-II methodological finding.

---

# 5. Phase 1.5 — close validation before new physics

## E1 — Packet-energy accounting audit

**Priority: highest.**

Because the instrument uses photon-number packets, verify explicitly that frequency shifts carry the correct energy weighting.

For every run record:

\[
E_{\rm injected},
\qquad
E_{\rm escaped},
\qquad
E_{\rm deposited/internal},
\]

and separately:

\[
N_{\rm injected},
\qquad
N_{\rm escaped}.
\]

Acceptance target:

\[
\left|
\frac{
E_{\rm escaped}+E_{\rm deposited}-E_{\rm injected}
}{
E_{\rm injected}
}
\right|
<10^{-3}
\]

or an appropriately strict Monte Carlo tolerance.

Add tests for:

- UV \(\rightarrow\) optical frequency shifts;
- multi-step cascades;
- packet-weight conversion;
- conservation under resonant reabsorption.

### Documentation updates

Add explicit packet convention and conservation diagnostics to:

- `docs/results_report.md §4.19`;
- lab notebook `§9p`;
- Paper-II Methods.

---

## E2 — Thermal-width reproduction

Repeat the central comparison at

\[
v_D=10\ {\rm km\,s^{-1}}
\]

and, if affordable,

\[
v_D=1\ {\rm km\,s^{-1}}.
\]

Minimum legs:

1. Sobolev pure absorption;
2. expansion pure absorption;
3. Sobolev + direct branching;
4. expansion + thermal.

Define

\[
\Delta_{\rm closure}
=
\frac{
F_{\rm exp+thermal}
-
F_{\rm Sob+branch}
}{
F_{\rm Sob+branch}
}.
\]

Goal:

> show that the fluorescence/thermal discrepancy survives in the regime where Paper I already established resolved \(\approx\) Sobolev.

Suggested figure:

\[
\Delta_{\rm closure}(v_D)
\]

for 1, 10, and 100 km s\(^{-1}\).

---

## E3 — Escape-probability regression suite

Add explicit tests for:

\[
\tau\ll1:\quad\beta\rightarrow1,
\]

\[
\tau\gg1:\quad\beta\rightarrow1/\tau.
\]

Also test:

- two-channel analytic branching;
- one branch becoming highly optically thick;
- repeated reabsorption;
- resonant-scattering invariance;
- effective \(A\beta\) competition.

---

# 6. Phase 2 — decisive literature-driven experiment

## E4 — Global \(\epsilon\) sweep

Run:

\[
\epsilon=
0,\ 0.1,\ 0.2,\ 0.3,\ 0.5,\ 0.7,\ 0.9,\ 1.
\]

For each case compare with direct A-branching.

Measure at minimum:

\[
F_{\rm UV},
\quad
F_{\rm blue},
\quad
F_{\rm optical},
\quad
F_{\rm red},
\quad
F_{\rm NIR}.
\]

For each band \(b\), define

\[
\epsilon_{\rm best}^{(b)}
=
\arg\min_\epsilon
\left|
F_b^{\rm TLA}(\epsilon)
-
F_b^{\rm branch}
\right|.
\]

The central test is:

\[
\boxed{
\epsilon_{\rm best}^{\rm UV}
\stackrel{?}{=}
\epsilon_{\rm best}^{\rm optical}
\stackrel{?}{=}
\epsilon_{\rm best}^{\rm NIR}
}
\]

rather than whether one optical band can be matched.

---

## E5 — Full-spectrum optimization

Compute:

\[
R_\lambda(\epsilon)
=
\frac{
F_\lambda^{\rm TLA}(\epsilon)
-
F_\lambda^{\rm branch}
}{
F_\lambda^{\rm branch}
}.
\]

Also define a global mismatch such as:

\[
\chi^2_{\rm spec}(\epsilon)
=
\sum_j
w_j
\left[
F_j^{\rm TLA}(\epsilon)
-
F_j^{\rm branch}
\right]^2.
\]

A scalar \(\epsilon\) that matches one band only because errors cancel elsewhere should not count as success.

---

# 7. Phase 2.5 — explain the redistribution

## E6 — Redistribution matrix

For escaped packets record:

\[
\lambda_{\rm pump},
\qquad
\lambda_{\rm escape}.
\]

Construct:

\[
P(\lambda_{\rm out}\mid\lambda_{\rm in}).
\]

Compare direct branching with TLA redistribution.

Candidate physical signature:

\[
\text{UV pump}
\rightarrow
\text{lanthanide cascade}
\rightarrow
\text{optical escape}.
\]

The thermal closure is expected to return more energy toward the LTE emissivity peak.

This should become a central explanatory figure.

---

## E7 — Dominant cascade pathways

For optical photons escaping in 3800–3955 Å, record:

- pump wavelength;
- initial upper level;
- branch sequence;
- number of reabsorptions;
- final transition;
- final escape wavelength;
- final-transition optical depth.

Rank pathways by escaped energy:

\[
f_k
=
\frac{
E_{\rm optical,\ pathway\ k}
}{
E_{\rm optical,\ total}
}.
\]

Goal:

> determine whether optical refill is dominated by a few identifiable branching families or by broad network statistics.

---


# 8. Phase 2.75 — realistic bulk-velocity robustness

## E13 — \(v_{\rm bulk}\sim0.1c\) / worldline-consistency test

This experiment is a **robustness test of the Paper-II closure result**, not another Doppler-width point.

The scientific question is

\[
\boxed{
\text{Does the conclusion about the scalar }\epsilon\text{ closure survive at }
v_{\rm bulk}\sim0.1c?
}
\]

### E13.1 Baseline geometry

Start with

\[
v_{\rm in}=0.05c,
\qquad
v_{\rm out}=0.15c,
\]

so

\[
v_{\rm bulk}\approx0.10c,
\qquad
\Delta v_{\rm shell}=0.10c.
\]

Keep the microscopic width independently fixed at

\[
v_D=10\ {\rm km\,s^{-1}},
\]

with an optional \(v_D=1\ {\rm km\,s^{-1}}\) confirmation.

### E13.2 Match line strength

Changing bulk velocity must not accidentally become an optical-depth experiment. Rescale the model as needed so that the high-velocity shell approximately matches the slow-shell Sobolev-depth distribution.

Compare at least

\[
\tau_{\max},
\qquad
N(\tau>1),
\qquad
N(\tau>0.1),
\]

and preferably

\[
P(\log\tau_S).
\]

### E13.3 Single-line worldline validation

Before the full La II run, test one line at

\[
\beta=0.01,\ 0.05,\ 0.10,\ 0.20.
\]

Frozen snapshot:

\[
\frac{\tau_{\rm frozen}}{\tau_S}
=
\frac{1-\beta}{\gamma}.
\]

Worldline-consistent transport:

\[
\frac{\tau_{\rm worldline}}{\tau_S}
=
\frac{1}{\gamma}.
\]

At \(\beta=0.1\), the physical worldline correction is only about \(0.5\%\), while a frozen snapshot can generate an apparent first-order effect of order \(10\%\). No high-velocity fluorescence result should be interpreted until this control passes.

### E13.4 Minimum full-ion comparison

Run:

1. worldline-consistent Sobolev + direct branching;
2. worldline-consistent TLA with \(\epsilon=1\);
3. worldline-consistent TLA with \(\epsilon_{\rm best}^{\rm slow}\);
4. two bracketing \(\epsilon\) values if necessary;
5. frozen versions of key legs as a transport-convention control.

Do not repeat the entire \(\epsilon\)-grid unless the optimum moves.

### E13.5 Main diagnostics

For each band \(b\),

\[
\Delta_b(\epsilon,v)
=
\frac{
F_b^{\rm TLA}(\epsilon,v)
-
F_b^{\rm branch}(v)
}{
F_b^{\rm branch}(v)
}.
\]

Compare

\[
\epsilon_{\rm best}(v_{\rm slow})
\]

with

\[
\epsilon_{\rm best}(0.1c).
\]

Also compare

\[
P(\lambda_{\rm out}\mid\lambda_{\rm in};v_{\rm slow})
\]

and

\[
P(\lambda_{\rm out}\mid\lambda_{\rm in};0.1c).
\]

### E13.6 Decision criteria

**Velocity robust:** if

\[
|\epsilon_{\rm best}(0.1c)-\epsilon_{\rm best}(v_{\rm slow})|
\lesssim0.1
\]

and the UV/optical/NIR residual pattern is unchanged, conclude that the closure result survives realistic kilonova bulk velocity.

**Velocity dependent:** if

\[
|\Delta\epsilon_{\rm best}|\gtrsim0.2
\]

or the redistribution structure changes substantially, promote \(v_{\rm bulk}\) to a validity-map axis.

**Transport convention dominates:** if frozen/worldline differences are comparable to the branching/TLA difference, do not use frozen high-\(v/c\) results to assess the closure.

### E13.7 Documentation

Add a dedicated results-report subsection containing:

- the three-velocity-scale distinction;
- high-velocity shell geometry;
- matched-\(\tau\) construction;
- single-line frozen/worldline validation;
- La II spectra;
- \(\epsilon_{\rm best}\) comparison;
- redistribution-matrix comparison;
- final robustness verdict.

---


# 9. Phase 3 — physics-informed closure

## E8 — Line-dependent fluorescence probability

In the radiative limit, compute:

\[
p_{{\rm fluor},i}
\simeq
\frac{
\sum_{k\neq i}A_{uk}\beta_{uk}
}{
\sum_k A_{uk}\beta_{uk}
}.
\]

Compare:

1. direct A-branching;
2. best global \(\epsilon\);
3. line-dependent \(p_{{\rm fluor},i}\).

High-value outcome:

\[
\boxed{
\text{a branching-aware expansion-opacity closure}
}
\]

that approaches direct-branching accuracy at much lower cost.

---

# 10. Phase 4 — second ion and mixtures

## E9 — Ce II

Repeat the minimum suite for Ce II:

- pure absorption;
- direct branching;
- global \(\epsilon\) sweep;
- spectral best fit;
- line-dependent closure.

Test:

\[
\epsilon_{\rm best}^{\rm La}
\stackrel{?}{=}
\epsilon_{\rm best}^{\rm Ce}.
\]

If not, the scalar closure is composition-dependent.

## E10 — La + Ce

Run the mixture and test the connection with Paper-I F7.

Questions:

1. Does dense blanketing suppress the redistribution error?
2. Does it merely move the discrepancy elsewhere in wavelength?
3. Does \(\epsilon_{\rm best}\) change with forest density?

---

# 11. Phase 5 — incident-spectrum dependence

## E11 — Source-temperature sweep

Run:

\[
T_{\rm src}
=
4000,\ 5000,\ 6000,\ 8000\ {\rm K}.
\]

Track:

\[
f_{\rm UV\rightarrow optical}.
\]

Goal:

> test whether the large optical fluorescence refill is driven primarily by the UV supply of the 6000 K Planck source.

Later, replace the Planck continuum with one realistic kilonova photospheric spectrum.

---

# 12. Phase 6 — LTE-state dependence

## E12 — Shell-temperature sweep

Run:

\[
T_{\rm shell}
=
2500,\ 3000,\ 4000,\ 5000\ {\rm K}.
\]

Recompute LTE populations.

Measure:

\[
\epsilon_{\rm best}(T)
\]

and the full spectral residual.

This probes sensitivity to:

- lower-level populations;
- line optical depths;
- escape probabilities;
- thermal emissivity.

---

# 13. Phase 7 — NLTE only if justified

Do not move to NLTE until the scalar-\(\epsilon\) question is answered.

NLTE introduces:

- statistical equilibrium;
- electron collisions;
- photoionization;
- recombination;
- nonthermal processes;
- much larger atomic-data uncertainty.

NLTE should be a later extension, not part of the minimum viable Paper II.

---

# 14. Decision gates

## Gate A — after E1–E3

Proceed only if:

- energy accounting closes;
- branching implementation is analytically validated;
- the sign reversal survives at \(v_D\le10\) km s\(^{-1}\).

## Gate B — after E4–E6

### One \(\epsilon\) works globally

Paper-II message:

> The SN-Ia TLA result largely carries over to La II, but complete thermalisation \(\epsilon=1\) does not.

Proceed to Ce II to test transferability.

### No scalar \(\epsilon\) works

Paper-II message:

\[
\boxed{
\text{The SN-Ia scalar-closure result does not carry over to La II.}
}
\]

Proceed to:

- redistribution diagnostics;
- line-dependent closure;
- Ce II.

### La and Ce prefer different \(\epsilon\)

Paper-II message:

> The thermalisation closure is ion-dependent and therefore not universal across lanthanide compositions.

---


## Gate C — after E13

The \(0.1c\) experiment is not required to decide whether a scalar \(\epsilon\) works in the controlled La II problem. It is required before making a strong claim that the same closure verdict applies directly to realistic kilonova expansion.

Proceed to broader astrophysical generalization if:

- the worldline implementation passes the single-line analytic controls;
- the high-velocity model is matched in optical-depth distribution to the slow reference;
- the qualitative branching/TLA conclusion survives at \(v_{\rm bulk}\sim0.1c\).

If the conclusion changes materially, treat \(v_{\rm bulk}\) as an explicit validity-map axis rather than treating the slow-shell result as universal.


# 15. Repository/documentation updates

## `docs/results_report.md §4.19`

Expand into:

### §4.19.1 Instrument
- line counts;
- packet number;
- seeds;
- temperatures;
- Doppler width;
- continuum;
- five treatments.

### §4.19.2 Validation
- Paper-I agreement;
- SEDONA thermal comparison;
- escape-probability correction;
- energy conservation.

### §4.19.3 Phase-1 result
Include the five-leg table.

### §4.19.4 Current interpretation
State explicitly:

- \(\epsilon=1\) fails;
- no conclusion yet about all scalar \(\epsilon\);
- sign reversal is provisional pending energy and thermal-width audits.

### §4.19.5 Literature connection
Summarize:

- SN-Ia direct fluorescence versus TLA;
- why \(\epsilon\sim0.3\) matters;
- why the lanthanide test is the new question.

---

## Lab notebook §9p

Record:

- commands;
- seeds;
- commit hashes;
- input hashes;
- line counts;
- branching statistics;
- wall times;
- missing-\(\beta\) diagnosis;
- before/after comparisons;
- energy accounting;
- every \(\epsilon\)-sweep configuration.

---

## README findings

### F19 — Escape-probability branching

Radiative branching must be weighted by

\[
\beta=(1-e^{-\tau})/\tau.
\]

### F20 — Fluorescent optical refill

For the full La II reference model, direct branching raises the 3800–3955 Å emergent flux from 0.183 to \(0.660\pm0.003\).

### F21 — Candidate closure sign reversal

Expansion opacity is too transparent in pure absorption (+88%) but, when paired with complete thermal redistribution (\(\epsilon=1\)), underproduces the same optical band by 37.5% relative to direct branching.

Qualifier:

> Whether this reflects failure of \(\epsilon=1\) specifically or failure of any scalar thermalisation closure is the next experiment.

---

# 16. Paper-II manuscript structure

## 1. Introduction

Narrative:

1. expansion opacity solves the line-count problem;
2. Paper I separated Sobolev from expansion opacity;
3. realistic transport requires wavelength redistribution;
4. SN-Ia work showed scalar TLA closures can approximate iron-peak fluorescence;
5. \(\epsilon=1\) is not generally optimal;
6. lanthanide branching networks are much more complex;
7. central question:
   \[
   \boxed{
   \text{Can one scalar }\epsilon\text{ reproduce them?}
   }
   \]

## 2. Methods

- same-code differential design;
- three velocity scales: $v_D$, $\Delta v_{\rm shell}$, and $v_{\rm bulk}$;
- per-line Sobolev interactions;
- expansion-opacity interactions;
- direct \(A\beta\)-weighted branching;
- TLA thermalisation parameter;
- packet-energy bookkeeping;
- GSI La II model;
- validation against Paper I and SEDONA.

## 3. Validation

- three-level atom;
- branching ratio;
- escape probability;
- energy conservation;
- pure absorption;
- thermal redistribution;
- Monte Carlo noise.

## 4. Results

- Phase-1 sign reversal;
- realistic-$v/c$ robustness and worldline control;
- thermal-width robustness;
- global \(\epsilon\) sweep;
- full-spectrum fit;
- redistribution matrix;
- dominant cascade pathways;
- line-dependent closure;
- Ce II / mixture generalization.

## 5. Discussion

Separate:

\[
\boxed{\text{opacity representation error}}
\]

from

\[
\boxed{\text{redistribution closure error}}.
\]

Position against:

- Kasen et al. direct fluorescence / TLA;
- Fontes et al.;
- Morag;
- TARDIS macroatom work;
- ARTIS line-by-line fluorescence;
- SUMO NLTE calculations.

---

# 17. Recommended execution order

1. **E1 — energy accounting**
2. **E2 — microscopic-width test: \(v_D=10\), then 1 km s\(^{-1}\) if feasible**
3. **E3 — \(\beta\) regression tests**
4. **E4 — global \(\epsilon\) sweep**
5. **E5 — full-spectrum mismatch**
6. **E6 — redistribution matrix**
7. **E7 — cascade pathways**
8. **E13 — macroscopic bulk-velocity test at \(v_{\rm bulk}\sim0.1c\), with frozen/worldline control**
9. **E8 — line-dependent closure**
10. **E9 — Ce II**
11. **E10 — La + Ce**
12. **E11 — incident-spectrum sweep**
13. **E12 — LTE-temperature sweep**
14. NLTE only afterward

---

# Immediate milestone

The next milestone is **not more realism**.

It is:

\[
\boxed{
\text{Determine whether any scalar }\epsilon\text{ can reproduce direct La II fluorescence.}
}
\]

This directly answers:

\[
\boxed{
\text{Does the SN-Ia TLA result carry over to lanthanide atomic networks?}
}
\]

The present work establishes only:

\[
\boxed{\epsilon=1\ \text{does not}.}
\]

The \(\epsilon\)-sweep determines whether the scientific answer becomes:

- **yes, after recalibration**, or
- **no, the scalar closure itself breaks down**.

That decision should determine the final identity of Paper II.
