# Sobolev Approximation in Kilonovae: Baby-Step Implementation Plan and Research Validation Roadmap

## Executive Summary

The best way to approach this project is **not** to begin by building a full kilonova radiative-transfer pipeline.

Because this would be a first hands-on radiative-transfer project, the initial objective should be:

> **Build confidence, falsify the idea early, and validate the technical path before investing heavily.**

The project should therefore proceed through several increasingly realistic stages:

| Phase | Goal | What to Build / Use | Decision |
|---|---|---|---|
| **0 — Weekend** | Prove that the Sobolev limit can be reproduced and deliberately broken | Small Python notebook + one real atomic dataset + SEDONA smoke test | Is the project technically tractable? |
| **1 — Weeks 1–2** | Learn the RT machinery and reproduce known behavior | Python reference solver + SEDONA resolved vs expansion-opacity calculations | Do I understand what the codes are doing? |
| **2 — Weeks 3–4** | Introduce realistic lanthanide line forests | GSI atomic data + narrow wavelength windows | Is there an interesting effect in realistic parameter space? |
| **3 — Months 2–3** | Build a Sobolev validity map | Parameter sweep in \(G\), \(O\), \(\tau\), \(\rho\), \(T\), \(t\) | Is there a journal paper? |
| **4 — Only if justified** | Connect to realistic kilonova observables | realistic ejecta, spectra, NLTE / fluorescence | Can this become a stronger paper? |

A **go/no-go decision** should be made after every phase.

---

# 1. Guiding Principle

Do not make SEDONA, ARTIS, a full kilonova spectrum, or a neutron-star-merger ejecta model the first task.

The first task should be something small enough that every equation and every line of code is understandable.

The preferred progression is:

\[
\boxed{
\text{Understand}
\rightarrow
\text{Reproduce}
\rightarrow
\text{Break deliberately}
\rightarrow
\text{Use real atomic data}
\rightarrow
\text{Scale up}
}
\]

The initial project is therefore not yet:

> “Build a journal-quality kilonova RT model.”

It is:

> **Determine whether this Sobolev-validity question is technically accessible, physically interesting, and worth a larger investment.**

---

# 2. Phase 0: Weekend Validation

The weekend test should avoid the hardest parts:

- no NLTE;
- no full kilonova spectrum;
- no realistic merger ejecta;
- no 3-D geometry;
- no radioactive heating;
- no millions of lines at once;
- no custom Monte Carlo code.

The weekend success criterion is:

> **Can I reproduce the Sobolev limit in a toy problem, make it fail in a controlled problem, inspect real lanthanide line crowding, and run at least one existing RT-code example?**

If yes, the project has a strong technical foothold.

---

# 3. Phase 0A — Single-Line Sobolev Experiment

This should be the **first implementation**.

Use Python only:

```text
numpy
scipy
matplotlib
pandas
```

No Monte Carlo.

No external RT code.

No atomic database initially.

Consider a one-dimensional homologous velocity field,

\[
v(x)=\frac{x}{t}.
\]

For a transition with rest frequency \(\nu_0\), calculate the resolved opacity

\[
\alpha_\nu(x)
=
\frac{\pi e^2}{m_e c}
f n_l(x)
\phi\left[\nu'(x)-\nu_0\right],
\]

where approximately

\[
\nu'(x)
=
\nu
\left(
1-\frac{v(x)}{c}
\right).
\]

Numerically calculate

\[
\tau_{\rm exact}(\nu)
=
\int
\alpha_\nu(x)\,dx.
\]

Compare this with the Sobolev optical depth

\[
\tau_S
=
\frac{\pi e^2}{m_e c}
f
n_l(x_{\rm res})
\lambda_0
t.
\]

---

## 3.1 First Unit Test: Constant Lower-Level Population

Set

\[
n_l(x)=n_0.
\]

The resolved integral should converge toward the Sobolev result as the numerical grid is refined.

Define

\[
E_{\rm Sob}
=
\left|
\frac{
\tau_{\rm exact}
-
\tau_S
}{
\tau_{\rm exact}
}
\right|.
\]

The first numerical result should demonstrate

\[
E_{\rm Sob}\rightarrow0.
\]

This should become a unit test.

If this calculation cannot reproduce the analytic Sobolev result, do not move on.

---

# 4. Phase 0B — Deliberately Break the Locality Assumption

Once the constant-population test works, introduce a gradient.

For example,

\[
n_l(v)
=
n_0
\left[
1
+
A
\tanh
\left(
\frac{v-v_{\rm res}}
{v_{\rm scale}}
\right)
\right].
\]

Define

\[
\epsilon
=
\frac{v_D}{v_{\rm scale}},
\]

where \(v_D\) is the intrinsic Doppler width.

Sweep approximately

\[
\epsilon
=
10^{-4},
10^{-3},
10^{-2},
3\times10^{-2},
0.1,
0.3,
1.
\]

The Sobolev approximation effectively samples the local population

\[
n_l(v_{\rm res}),
\]

whereas the resolved calculation integrates through the actual variation across the resonance region.

Plot

\[
E_{\rm Sob}
\]

against

\[
\frac{v_D}{v_{\rm scale}}.
\]

The expected qualitative behavior is

\[
\epsilon\ll1
\quad\Rightarrow\quad
E_{\rm Sob}\ll1,
\]

and

\[
\epsilon\rightarrow0.1-1
\quad\Rightarrow\quad
E_{\rm Sob}
\text{ grows}.
\]

This becomes the conceptual seed of a later validity map.

---

# 5. Important Caveat About a Two-Line Weekend Test

A pure two-line absorption experiment should **not** be the main Phase 0 line-overlap test.

For two absorption lines,

\[
\alpha_\nu
=
\alpha_{\nu,1}
+
\alpha_{\nu,2},
\]

and the total survival probability can still behave approximately like

\[
e^{-(\tau_1+\tau_2)}.
\]

Therefore line-profile overlap does not automatically create a dramatic difference in every integrated observable.

A naive two-line absorption experiment may produce little difference and create the false impression that line overlap is irrelevant.

For Phase 0, the safer first question is simply:

> **Does physically significant line crowding occur in real lanthanide data on thermal-width velocity scales?**

The transport consequences can be studied later with SEDONA.

---

# 6. Phase 0C — Introduce One Real Lanthanide Dataset

Use modern calibrated lanthanide atomic data.

The preferred starting dataset is the **GSI v2 lanthanide database**, which includes singly and doubly ionized lanthanides La–Yb.

Relevant quantities include:

- lower and upper energy levels;
- angular-momentum information;
- wavelengths;
- \(\log(gf)\);
- Einstein \(A\)-values;
- calibration status.

For the first experiment, choose only one ion, preferably:

\[
\mathrm{La\,II}.
\]

Load quantities such as:

```text
wavelength
log_gf
E_lower
E_upper
method_lower
method_upper
```

Create two subsets.

### High-confidence subset

For example:

- lower level experimentally matched;
- upper level experimentally matched;
- \(\log(gf)>-1\).

### Broader subset

For example:

- experimentally matched or calibrated / shifted transitions;
- \(\log(gf)>-2\).

Sort by wavelength and calculate nearest-neighbor velocity spacing:

\[
\Delta v_i
=
c
\frac{
\lambda_{i+1}
-
\lambda_i
}{
\lambda_i
}.
\]

Plot the distribution

\[
P(\Delta v).
\]

Mark representative intrinsic widths:

\[
v_D
=
1,\ 3,\ 10
\ {\rm km\,s^{-1}}.
\]

Define a simple overlap measure

\[
O_i
=
\frac{v_D}
{\Delta v_i}.
\]

Then

\[
O_i>1
\]

means the neighboring transition lies within approximately one Doppler width.

This does **not** prove Sobolev breakdown.

It answers the more basic Phase 0 question:

> **Does the potentially problematic regime exist in modern calibrated lanthanide data?**

---

# 7. Phase 0D — Get SEDONA Running

SEDONA should be the preferred main research code, but it should not dominate the weekend.

The code provides both:

```text
opacity_line_expansion = 1
```

for Sobolev / expansion-opacity treatment, and

```text
opacity_bound_bound = 1
```

for resolved bound-bound transitions using line profiles.

The first task should be to run existing examples rather than invent a model.

Recommended sequence:

```text
1. spherical_lightbulb
2. Type-Ia expansion-opacity example
3. Type-Ia resolved bound-bound example
```

The goal is simply to learn:

- input structure;
- model structure;
- opacity switches;
- atomic-data format;
- output spectra;
- runtime behavior.

If SEDONA installation starts consuming most of the weekend, stop.

The Python toy model is the core Phase 0 deliverable.

On Windows, use **WSL2 / Linux** rather than making a native Windows build part of the research problem.

---

# 8. Weekend Deliverables

A successful Phase 0 should ideally produce three concrete outputs.

## Figure 1 — Sobolev Convergence

Plot

\[
E_{\rm Sob}
\]

versus

\[
v_D/v_{\rm scale}.
\]

Show the expected asymptotic Sobolev regime.

---

## Figure 2 — Real Lanthanide Line Crowding

For La II, plot

\[
P(\Delta v_{\rm nearest})
\]

with reference widths at

\[
1,\ 3,\ 10
\ {\rm km\,s^{-1}}.
\]

This establishes whether line crowding occurs on intrinsic-width scales.

---

## Result 3 — Existing RT Code Smoke Test

Ideally:

- SEDONA lightbulb test runs;
- expansion-opacity example runs;
- resolved bound-bound example runs.

At minimum:

- understand the directory structure;
- understand the relevant parameters;
- confirm that the code can in principle run both treatments.

---

# 9. Things Not to Do During Phase 0

Avoid:

- realistic merger ejecta;
- full Saha ionization equilibrium;
- NLTE;
- fluorescence;
- a custom Monte Carlo transport code;
- millions of lines simultaneously;
- full UV-to-IR spectra;
- 3-D geometry;
- radioactive heating;
- complete kilonova light curves.

These are distractions before the core approximation is understood.

---

# 10. Component Choices and Alternatives

| Component | Recommended First Choice | Alternative | Reason |
|---|---|---|---|
| Reference RT | **Small Python deterministic solver** | SEDONA resolved BB | Maximum transparency |
| Production comparison | **SEDONA** | Custom solver | Both resolved and expansion-opacity treatments exist in one framework |
| Learning RT code | **SEDONA examples** | TARDIS | TARDIS is easier to inspect but uses Sobolev line interactions |
| Later NLTE / fluorescence validation | **ARTIS** | TARDIS macroatom | Strong NLTE / macroatom capabilities, but line transport remains Sobolev |
| Broad kilonova light curves | **SuperNu** | ARTIS / SEDONA | Useful later, not optimal for testing Sobolev itself |
| Atomic data first | **GSI v2 calibrated** | HFR v2 | Better wavelength reliability |
| Broad atomic coverage | **HFR v2** | GSI | Much wider element / ion coverage |
| Population model first | **Prescribed \(n_l\)** | Boltzmann LTE | Isolate radiative transfer |
| Population model second | **Boltzmann LTE** | Saha + Boltzmann | Introduce thermal populations gradually |
| Later population model | **Saha + Boltzmann LTE** | SEDONA LTE | More realistic but still controlled |
| Line profile first | **Gaussian** | Voigt | Easier for understanding |
| Production line profile | **Voigt** | Gaussian sensitivity test | Includes natural broadening |
| Ejecta first | **Analytic 1-D homologous** | Kasen-style parameterized model | Maximum control |
| Realistic parameter ranges | **Radice ejecta data** | Public kilonova grids | Introduce only after toy models |
| Workflow | **Python + HDF5 / simple scripts** | Snakemake later | Avoid over-engineering early |

---

# 11. Why SEDONA Is the Preferred Main Research Code

For this particular scientific question, SEDONA has a major advantage:

\[
\boxed{
\text{resolved line profiles}
\quad\text{and}\quad
\text{Sobolev / expansion opacity}
}
\]

exist in the same framework.

That means the experiment can hold fixed:

- geometry;
- ejecta structure;
- atomic data;
- populations;
- boundary conditions;
- Monte Carlo treatment;

while changing primarily the line-transfer approximation.

This is ideal for a controlled comparison.

A previous SEDONA application has already demonstrated frequency-dependent line treatment in the comoving frame without relying on the expansion-opacity / Sobolev approximation.

The known challenge is that resolving intrinsic widths of only

\[
1-10\ {\rm km\,s^{-1}}
\]

requires approximately

\[
\delta\nu/\nu
\sim10^{-6},
\]

which makes full-spectrum calculations expensive.

This motivates narrow wavelength windows.

---

# 12. Why ARTIS Is Not the First Code

ARTIS is highly valuable for later stages because it includes:

- multidimensional homologous ejecta;
- macroatoms;
- fluorescence;
- NLTE populations;
- radioactive decay;
- nonthermal physics.

However, its default individual line interactions remain based on Sobolev treatment.

Therefore ARTIS is particularly useful for later questions such as:

> **Does fluorescence alleviate extreme-\(\tau_S\) trapping?**

It is less appropriate as the non-Sobolev reference calculation for the initial paper.

---

# 13. Why SuperNu Is Not the First Code

SuperNu is a strong kilonova transport framework and is widely used for:

- time-dependent transport;
- homologous ejecta;
- light curves;
- spectra;
- multigroup calculations.

However, it naturally works with grouped opacity treatments.

That makes it excellent for kilonova observables but less suitable for the fundamental question:

> **What error occurs when a narrow physical line profile is collapsed into a Sobolev resonance?**

Use SuperNu later if the validity-map results need to be translated into light-curve consequences.

---

# 14. Role of TARDIS

TARDIS is attractive educationally because it provides:

- Python-facing workflows;
- well-organized atomic-data structures;
- accessible Monte Carlo concepts;
- relatively low barrier to experimentation.

However, its line-interaction treatment uses Sobolev optical depths.

Therefore:

\[
\boxed{
\text{TARDIS = good RT learning tool}
}
\]

but not

\[
\boxed{
\text{non-Sobolev reference solver}.
}
\]

It is optional rather than central.

---

# 15. Atomic-Data Strategy

## 15.1 Primary Dataset — GSI v2

Use this first.

Reasons:

- designed for kilonova radiative transfer;
- calibrated wavelength information;
- manageable size;
- covers La–Yb II/III;
- appropriate for line-spacing / overlap studies.

Start with one ion such as

\[
\mathrm{La\,II},
\]

then expand to

\[
\mathrm{Ce\,II/III}
\]

and

\[
\mathrm{Nd\,II/III}.
\]

---

## 15.2 Secondary Dataset — HFR v2

Use later for broader coverage.

Advantages:

\[
Z=20-103,
\]

with neutral through triply ionized species.

This includes:

- light r-process species;
- lanthanides;
- actinides.

However, the dataset is much larger and much more theoretical.

For a study sensitive to separations of only a few km/s, wavelength uncertainty is especially important.

Therefore the sequence should be

\[
\boxed{
\text{GSI}
\rightarrow
\text{establish physics}
\rightarrow
\text{HFR robustness test}.
}
\]

---

# 16. Population Modeling Strategy

Start even simpler than full LTE.

## Stage A — Prescribed Lower-Level Population

Choose one ion and directly set

\[
n_l.
\]

This isolates radiative transfer.

---

## Stage B — Boltzmann Populations

Calculate

\[
\frac{n_l}{n_{\rm ion}}
=
\frac{g_l}{Z(T)}
e^{-E_l/kT}.
\]

This introduces temperature dependence without introducing ionization uncertainty.

---

## Stage C — Saha + Boltzmann LTE

Introduce ionization balance such as

\[
\mathrm{La\,II}
\leftrightarrow
\mathrm{La\,III}.
\]

Cross-check the result with SEDONA's LTE treatment.

---

## Stage D — NLTE

Only after the Sobolev-validity question is already scientifically interesting.

Otherwise discrepancies can become difficult to interpret because they may originate from:

- transfer approximation;
- ionization;
- partition functions;
- collision rates;
- atomic-data mismatch;
- line-selection thresholds.

---

# 17. Ejecta Strategy

Use increasingly realistic models.

## Level 0 — Homogeneous Slab

Use

\[
v=x/t.
\]

Ideal for code validation.

---

## Level 1 — 1-D Homologous Sphere

For example,

\[
\rho(v)
\propto
v^{-n}.
\]

Introduce controlled population or abundance gradients.

---

## Level 2 — Parameterized Kilonova Ejecta

Use a Kasen-like model characterized by quantities such as

\[
M_{\rm ej},
\quad
v_{\rm ej},
\quad
X_{\rm lan}.
\]

This establishes realistic scales without requiring a full merger simulation.

---

## Level 3 — Public Merger Ejecta

Use public merger data containing quantities such as

\[
Y_e,
\quad
v_\infty,
\quad
s,
\quad
\tau_{\rm expansion}.
\]

This becomes relevant once the basic transfer effect is understood.

---

## Level 4 — Multidimensional RT Ejecta

Only much later.

Do not begin with 2-D or 3-D geometry before the one-dimensional validity boundary is known.

---

# 18. Weeks 1–2: Confidence-Building Phase

The next two weeks should be treated as radiative-transfer exercises rather than paper production.

## Week 1

Implement and validate:

```text
1. Gaussian line normalization
2. Voigt line normalization
3. Exact optical-depth integration
4. Sobolev optical depth
5. Boltzmann level populations
```

Add numerical tests.

For example,

\[
\int
\phi_\nu
\,d\nu
=
1.
\]

And verify that constant-population homologous expansion reproduces

\[
\tau_S.
\]

Then run SEDONA examples:

```text
spherical_lightbulb
TypeIa expansion-opacity
TypeIa resolved bound-bound
```

The objective is to become comfortable with:

- model files;
- parameter files;
- atomic-data inputs;
- output spectra;
- runtime and numerical resolution.

---

# 19. Week 2: Build a Minimal Controlled SEDONA Model

Create:

\[
1
\text{ element}
+
1
\text{ ion}
+
1
\text{ line}
+
1\text{D}.
\]

Do not use lanthanides yet.

Compare:

\[
\text{Python deterministic solver}
\]

against

\[
\text{SEDONA resolved bound-bound}.
\]

Then progress through

\[
1
\text{ line}
\rightarrow
2
\text{ lines}
\rightarrow
20
\text{ lines}.
\]

At the end of this stage, the goal is to understand every important term entering the calculation.

This is the largest initial learning barrier.

---

# 20. Weeks 3–4: First Actual Research Calculation

Now introduce realistic GSI data.

Recommended ion sequence:

\[
\mathrm{La\,II}
\rightarrow
\mathrm{Ce\,II/III}
\rightarrow
\mathrm{Nd\,II/III}.
\]

Use narrow wavelength windows, perhaps initially

\[
5-20\ {\rm \AA}.
\]

Sweep quantities such as

\[
T,
\quad
\rho,
\quad
t,
\quad
v_D.
\]

Calculate

\[
F_\lambda^{\rm resolved}
\]

and

\[
F_\lambda^{\rm Sob}.
\]

Start relating the discrepancy to the physical control parameters.

---

# 21. Core Validity Parameters

## 21.1 Gradient Parameter

\[
\boxed{
G
=
v_D
\left|
\frac{d\ln n_l}{dv}
\right|
}
\]

or generalized to quantities such as

\[
\rho,
\quad
T,
\quad
X_i,
\quad
n_{\rm ion}.
\]

---

## 21.2 Line-Overlap Parameter

A simple version is

\[
\boxed{
O
=
\frac{
v_D
}{
\Delta v_{\rm strong}
}
}
\]

or equivalently the number of optically important transitions within a Doppler width.

---

## 21.3 Extreme Optical-Depth Parameter

A diagnostic quantity is

\[
\boxed{
T_{\rm trap}
\sim
\tau_S
\frac{v_D}{c}.
}
\]

This is not yet a complete treatment of fluorescence or redistribution.

It identifies regimes requiring more sophisticated investigation.

---

# 22. First Validity Map

The eventual map should connect

\[
(G,O,T_{\rm trap})
\]

to a numerical Sobolev error.

For example,

\[
\Delta_{\rm Sob}
=
\frac{
Q_{\rm Sobolev}
-
Q_{\rm resolved}
}{
Q_{\rm resolved}
}.
\]

Possible choices for \(Q\) include:

\[
\tau_\lambda,
\quad
\kappa_\lambda,
\quad
J_\nu,
\quad
F_\lambda,
\quad
W_\lambda.
\]

Later this could be extended to:

\[
L_{\rm bol},
\quad
m_{\rm band}.
\]

---

# 23. One-Month Go / No-Go Decision

Define this decision before starting so that sunk cost does not drive the project.

## Strong GO

Physically plausible kilonova conditions produce

\[
|\Delta F_\lambda|
\gtrsim
10\%.
\]

Or the discrepancy shows a strong and reproducible dependence on

\[
G,\quad O,\quad \tau_S.
\]

This is a very strong research direction.

---

## Moderate GO

Differences are only a few percent, but a clean empirical validity boundary emerges.

For example,

\[
G<0.03,
\qquad
O<0.1
\]

might imply

\[
E_{\rm Sob}<2\%.
\]

A robust quantitative validation of a widely used approximation can still support a methodological journal paper.

---

## Reassess

All realistic kilonova conditions give

\[
E_{\rm Sob}<1\%,
\]

and substantial disagreement appears only for clearly unrealistic parameters.

This is scientifically useful, but the publication value should be weighed against the cost of extending the calculation.

---

## Stop or Pivot

If the comparison is dominated by:

- uncertain wavelengths;
- incompatible opacity implementations;
- unisolatable differences between codes;
- numerical-resolution artifacts;

then pivot toward a related but cleaner question such as:

> **line-overlap uncertainty in kilonova opacity calculations**

rather than forcing a Sobolev-validity paper.

---

# 24. Recommended Initial Repository

Keep the codebase intentionally small:

```text
sobolev-kilonova/
│
├── notebooks/
│   ├── 00_single_line.ipynb
│   ├── 01_gradient_breakdown.ipynb
│   └── 02_gsi_line_spacing.ipynb
│
├── src/
│   ├── profiles.py
│   ├── sobolev.py
│   ├── formal_transfer.py
│   └── atomic_data.py
│
├── tests/
│   ├── test_profiles.py
│   └── test_sobolev_limit.py
│
└── README.md
```

Do not begin with a complicated package architecture.

If the project survives the first month, restructure it.

---

# 25. What the First Notebook Should Contain

The first notebook should be:

```text
00_single_line.ipynb
```

It should contain:

1. physical constants;
2. Gaussian profile;
3. profile normalization test;
4. homologous velocity field;
5. exact optical-depth integration;
6. Sobolev optical depth;
7. convergence study;
8. controlled \(n_l(v)\) gradient;
9. Sobolev error versus \(v_D/v_{\rm scale}\);
10. short interpretation.

That notebook should become the foundation for every later calculation.

---

# 26. Recommended Weekend Schedule

A practical weekend sequence is:

## Session 1

Implement:

- Gaussian profile;
- homologous velocity field;
- numerical optical-depth integration;
- Sobolev formula.

Verify constant-population convergence.

---

## Session 2

Introduce the controlled population gradient.

Generate:

\[
E_{\rm Sob}
\quad\text{vs}\quad
v_D/v_{\rm scale}.
\]

This is the most important weekend result.

---

## Session 3

Download / inspect one GSI atomic file.

Calculate:

\[
P(\Delta v_{\rm nearest})
\]

for one ion.

---

## Session 4

Attempt SEDONA installation and run the simplest existing example.

If installation is troublesome, stop rather than sacrificing the whole weekend.

---

# 27. Research-Investment Philosophy

The project should answer increasingly expensive questions.

First:

\[
\boxed{
\text{Can I independently reproduce the Sobolev result?}
}
\]

Then:

\[
\boxed{
\text{Can I make its assumptions fail in a controlled problem?}
}
\]

Then:

\[
\boxed{
\text{Do modern r-process atomic data occupy those failure regimes?}
}
\]

Then:

\[
\boxed{
\text{Does the failure materially affect radiation transport?}
}
\]

Only then:

\[
\boxed{
\text{Does it alter observable kilonova spectra or light curves?}
}
\]

This sequence minimizes research risk.

---

# 28. Overall Recommendation

This is a reasonable direction to investigate, but the first four weeks should be viewed as an **exploration and validation project**, not yet as a committed journal project.

The project has several favorable features:

- modern public calibrated lanthanide data now exist;
- SEDONA already provides both resolved and Sobolev / expansion-opacity treatments;
- the fundamental test can be formulated in simple one-dimensional problems;
- meaningful go/no-go criteria can be established early;
- the initial computation does not require extreme HPC resources;
- a positive or negative result can both be scientifically informative.

The strongest initial strategy is therefore:

\[
\boxed{
\text{small deterministic solver}
\rightarrow
\text{controlled Sobolev failure}
\rightarrow
\text{real line crowding}
\rightarrow
\text{SEDONA comparison}
\rightarrow
\text{validity map}
}
\]

The first implementation target should be:

\[
\boxed{
\texttt{00\_single\_line.ipynb}
}
\]

Once that notebook works and the physical behavior is understood, the project becomes much easier to evaluate rationally.
