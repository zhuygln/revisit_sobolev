# Sobolev Approximation in Kilonovae: Practical Research Requirements, Code, and Data Availability

## Executive Summary

A researcher starting from scratch does **not** need to build a complete non-Sobolev radiative-transfer code to study the validity of the Sobolev approximation in kilonova ejecta.

The most promising strategy is to use existing transport software—especially **SEDONA**, which supports both:

- resolved bound-bound opacity with Voigt profiles, and
- Sobolev/expansion-opacity treatments,

and then build a controlled validation framework around it.

The project therefore becomes:

> **Build a controlled Sobolev-validation framework using existing transport software, public atomic data, a small independent reference solver, and realistic kilonova parameter ranges.**

The main missing-data problem is **not** basic atomic line data. It becomes serious mainly when extending the study to fully self-consistent **NLTE** populations, fluorescence, and extreme optical-depth trapping.

---

# 1. What Code Actually Needs to Be Developed

The project can be divided into six software components.

| Component | Develop Yourself? | Difficulty | Purpose |
|---|---:|---:|---|
| SEDONA Monte Carlo transport | No | — | Existing resolved and Sobolev transport engine |
| Atomic-data converter | Yes | Low–medium | Convert GSI/HFR atomic data into internal/SEDONA formats |
| Controlled ejecta generator | Yes | Low | Generate homologous density, temperature, velocity, and composition models |
| Experiment / parameter-sweep framework | Yes | Medium | Run resolved and Sobolev calculations under identical conditions |
| Validity-map diagnostics | Yes | Medium | Compute gradient, overlap, optical-depth, and error measures |
| Small independent formal-transfer solver | Yes, recommended | Medium | Independently verify controlled non-Sobolev results |

The most important point is that **the radiative-transfer engine itself does not need to be written from scratch**.

---

# 2. SEDONA as the Main Transport Framework

SEDONA is designed for radiative transfer in homologously expanding transient ejecta and accepts model inputs such as

\[
\rho(r), \qquad T(r), \qquad v(r), \qquad X_i(r).
\]

It supports 1-D spherical models directly.

More importantly, its opacity implementation includes two relevant modes.

## 2.1 Resolved bound-bound opacity

SEDONA can represent individual lines using resolved Voigt profiles.

Conceptually,

\[
\alpha_\nu
=
\sum_i
\frac{\pi e^2}{m_e c}
f_i n_{l,i}\phi_i(\nu).
\]

This provides the basis for a non-Sobolev reference calculation.

## 2.2 Sobolev / expansion opacity

SEDONA can also calculate line opacity using Sobolev optical depths and the expansion-opacity construction,

\[
1-e^{-\tau_{\rm Sob}}.
\]

This makes it possible to compare the two approaches in the **same transport framework**.

## 2.3 Why this matters

The desired controlled experiment is therefore approximately

\[
\boxed{
\text{same ejecta}
+
\text{same atomic data}
+
\text{same populations}
}
\]

with only the line-transfer treatment changed:

\[
\text{resolved profiles}
\quad\leftrightarrow\quad
\text{Sobolev}.
\]

A previous SEDONA application has already demonstrated that frequency-dependent bound-bound opacities can be treated in the comoving frame without using the Sobolev or line-expansion approximation.

The main numerical difficulty is frequency resolution rather than development of the transfer algorithm itself.

---

# 3. Proposed Project Repository

A practical project could be organized as follows:

```text
sobolev-kilonova/
    atomic/
        gsi_reader
        hfr_reader
        sedona_converter
        line_selector

    ejecta/
        analytic_homologous
        composition_profiles
        realistic_ejecta_import

    experiments/
        single_line
        two_line_overlap
        line_forest
        gradient_test
        kilonova_grid

    reference_rt/
        formal_solver_1d

    diagnostics/
        sobolev_depth
        gradient_parameter
        overlap_parameter
        trapping_parameter
        error_metrics

    workflows/
        run_resolved
        run_sobolev
        convergence
        validity_map
```

---

# 4. Atomic-Data Infrastructure

A standardized internal atomic-data representation should contain quantities such as

\[
\{
\lambda_i,
E_l,
E_u,
g_l,
g_u,
f_{lu},
A_{ul},
Z,
\text{ion}
\}.
\]

The same exact transition set should then be used for both resolved and Sobolev calculations.

This is important scientifically: otherwise differences between the two calculations could arise from inconsistent line filtering or atomic-data preprocessing rather than the transport approximation itself.

---

# 5. Public Atomic Data Available Now

## 5.1 GSI calibrated lanthanide database

A particularly attractive starting point is the public GSI lanthanide database.

It covers singly and doubly ionized lanthanides:

\[
Z=57-70,
\qquad
\mathrm{II,\ III}.
\]

Available quantities include:

- lower and upper energy levels,
- angular-momentum information,
- wavelengths,
- \(\log(gf)\),
- Einstein \(A\)-values,
- calibration status of levels and transitions.

The underlying calculations contain roughly

\[
146{,}856 \text{ levels}
\]

and

\[
28.69 \text{ million transitions}.
\]

A subset of transitions has experimentally calibrated wavelength information.

This is sufficient to begin a serious LTE Sobolev-vs-resolved study.

## 5.2 Broader HFR database

A broader theoretical database covers approximately

\[
Z=20-103
\]

from Ca through Lr and ionization stages I–IV.

This provides access to:

- light r-process species,
- lanthanides,
- actinides,
- neutral through triply ionized material.

However, it should probably be used **after** the calibrated lanthanide data.

For line-overlap studies, wavelength accuracy matters strongly. If the physical question depends on whether two transitions are separated by only a few km/s, theoretical wavelength uncertainty can change whether they are predicted to overlap.

Therefore the preferred strategy is

\[
\boxed{
\text{calibrated GSI lines first}
\rightarrow
\text{broader theoretical databases later}.
}
\]

---

# 6. Atomic Data That Are Still Missing

For an LTE Sobolev validation study, the available data are already sufficient.

If LTE level populations are imposed,

\[
n_i=n_i(T,\rho,X),
\]

then wavelengths, oscillator strengths, energies, degeneracies, and Einstein coefficients allow calculation of

\[
n_l,
\qquad
\tau_S,
\qquad
\alpha_\nu,
\qquad
\phi_\nu.
\]

Therefore:

\[
\boxed{
\text{atomic-data incompleteness is not a blocker for the first paper.}
}
\]

The major data gap appears when attempting fully self-consistent **NLTE** calculations.

These require quantities such as

\[
C_{ij}(T)
\]

for electron-impact excitation and de-excitation,

\[
\sigma_{\rm PI}(\nu)
\]

for photoionization,

\[
\alpha_{\rm RR}(T)
\]

for radiative recombination,

\[
\alpha_{\rm DR}(T)
\]

for dielectronic recombination,

as well as potentially:

- electron-impact ionization,
- charge exchange,
- detailed collisional coupling among excited states.

These data remain incomplete for many r-process species.

A sensible project boundary is therefore:

\[
\boxed{
\begin{array}{ll}
\text{Paper I:} & \text{fixed LTE populations; test transport approximation} \\[2mm]
\text{Paper II:} & \text{NLTE + fluorescence + extreme optical-depth trapping}
\end{array}
}
\]

---

# 7. Ejecta Models Needed

There are three useful levels of ejecta modeling.

## 7.1 Controlled analytic homologous ejecta

This should be the starting point.

Assume

\[
v=\frac{r}{t}
\]

and specify, for example,

\[
\rho(v)
=
\rho_0
\left(
\frac{v}{v_0}
\right)^{-n}.
\]

Composition can be given a tunable gradient such as

\[
X_{\rm Nd}(v)
=
X_0
\left[
1
+
A\tanh
\left(
\frac{v-v_c}{\Delta v}
\right)
\right].
\]

This lets the researcher directly control the Sobolev-locality parameter

\[
G_X
=
v_D
\left|
\frac{d\ln X}{dv}
\right|.
\]

Controlled models are preferable for establishing causal relationships.

---

# 8. Published Parameterized Kilonova Models

Existing public kilonova grids can be used to establish realistic parameter ranges.

Examples include parameterized models varying:

- ejecta mass,
- characteristic velocity,
- lanthanide fraction,
- geometry,
- viewing angle.

These are useful for answering:

> Are the density, velocity, temperature, and composition ranges used in the controlled experiment representative of kilonova simulations?

However, many public model releases primarily provide synthetic spectra, light curves, and magnitudes rather than complete cell-by-cell radiative-transfer states.

They are therefore useful for context but not necessarily ideal as direct inputs to the first validity experiment.

---

# 9. Public Neutron-Star Merger Ejecta Data

Public merger-simulation datasets do exist.

Useful released quantities can include

\[
v_\infty,
\qquad
Y_e,
\qquad
s,
\qquad
\tau_{\rm expansion},
\]

as well as multidimensional histograms and nucleosynthesis products.

These allow realistic merger conditions to be mapped into the validity study.

However, they are not necessarily provided directly as late-time radiation-transfer grids containing

\[
\rho(v,\theta,t),
\quad
T(v,\theta,t),
\quad
X_Z(v,\theta,t)
\]

at epochs such as

\[
0.1,\ 0.3,\ 1,\ 3~{\rm days}.
\]

Constructing that mapping is a manageable but nontrivial part of the project.

---

# 10. Important Ejecta Data Not Yet Verified as Readily Available

The ideal external dataset would provide, on the same homologous radiative-transfer grid,

\[
\boxed{
\rho(\mathbf v),
\quad
T(\mathbf v),
\quad
Y_e(\mathbf v),
\quad
X_Z(\mathbf v)
}
\]

together with network-derived composition.

Such complete publication-ready late-time datasets may exist for individual studies, but their broad public availability should be verified before relying on them.

This is **not a blocker for the first paper**.

Analytic ejecta plus publicly available merger parameter distributions are sufficient for the initial study.

---

# 11. The Main Computational Challenge: Frequency Resolution

The dominant cost of resolved-line calculations is not the transport algorithm.

It is the frequency resolution required to resolve thermal line profiles.

For a heavy ion with

\[
v_D\sim1~{\rm km\,s^{-1}},
\]

the Doppler width is approximately

\[
\frac{\Delta \nu_D}{\nu}
\sim
\frac{v_D}{c}
\sim
3\times10^{-6}.
\]

Several frequency points are required across a profile, suggesting

\[
\frac{\delta\nu}{\nu}
\lesssim
10^{-6}.
\]

Resolving such widths over the entire UV-to-IR spectrum with millions of transitions becomes extremely expensive.

Therefore the initial study should use

\[
\boxed{
\text{small wavelength windows}
+
\text{realistic line forests}.
}
\]

Examples might include windows such as

\[
4000-4010~\text{\AA}
\]

or

\[
8000-8020~\text{\AA}.
\]

The wavelength-window size can then be increased systematically until results converge.

---

# 12. A Line-Window Engine

A useful new software component is an indexed line-selection engine.

Given a database with tens of millions of transitions, the workflow could be:

1. select ion(s);
2. select wavelength interval;
3. calculate LTE populations;
4. estimate \(\tau_S\);
5. retain transitions above an importance threshold;
6. enlarge the wavelength interval to include velocity Doppler shifts;
7. generate the transport-code atomic input;
8. run both resolved and Sobolev treatments.

The same infrastructure can calculate diagnostics such as

\[
N_{\rm lines}(\Delta v_D),
\]

\[
N_{\tau_S>1}(\Delta v_D),
\]

and the nearest-neighbor separation distribution

\[
P(\Delta v_{\rm nearest}).
\]

These become direct measures of line overlap.

---

# 13. Independent Reference Transfer Solver

Even if SEDONA supplies both transport modes, an independent solver is strongly recommended.

For a prescribed 1-D atmosphere, solve

\[
\frac{dI_\nu}{ds}
=
-\alpha_\nu I_\nu+j_\nu.
\]

Under LTE true absorption,

\[
S_\nu
=
\frac{j_\nu}{\alpha_\nu}
=
B_\nu(T).
\]

For a discretized characteristic,

\[
I_{\nu,k+1}
=
I_{\nu,k}
e^{-\Delta\tau_{\nu,k}}
+
S_{\nu,k}
\left(
1-e^{-\Delta\tau_{\nu,k}}
\right).
\]

The main additional ingredient is transforming frequency into the local comoving frame for

\[
v(r)=\frac{r}{t}.
\]

The validation sequence should be

\[
1~\text{line}
\rightarrow
2~\text{lines}
\rightarrow
10~\text{lines}
\rightarrow
\text{small realistic line forest}.
\]

The solver should recover the Sobolev result when

\[
\frac{v_D}{v_{\rm scale}}
\rightarrow0
\]

and lines are well separated.

Then the assumptions can be broken deliberately.

---

# 14. Validity Diagnostics

The study should calculate several dimensionless parameters.

## 14.1 Gradient/locality parameter

\[
\boxed{
G_Q
=
v_D
\left|
\frac{d\ln Q}{dv}
\right|
}
\]

for quantities such as

\[
Q=
\rho,
T,
n_{\rm ion},
n_l,
X_i.
\]

The Sobolev locality condition is

\[
G_Q\ll1.
\]

## 14.2 Line-overlap parameter

One useful definition is

\[
\boxed{
O
=
N_{\tau>1}(\Delta\nu_D)
}
\]

or alternatively

\[
O
=
\frac{\Delta\nu_D}
{\langle\Delta\nu_{\rm strong}\rangle}.
\]

The isolated-line regime corresponds approximately to

\[
O\ll1.
\]

## 14.3 Optical-depth / trapping diagnostic

A first-order diagnostic is

\[
\boxed{
T_{\rm trap}
\sim
\tau_S
\frac{v_D}{c}.
}
\]

This is not intended as the final treatment of fluorescence or redistribution, but it identifies potentially problematic extreme optical-depth regimes.

---

# 15. Desired Validity Map

The final product should map physically relevant parameters into a measurable Sobolev error.

Conceptually,

\[
(G,O,T_{\rm trap})
\longrightarrow
\Delta_{\rm Sob}.
\]

Possible error measures include

\[
\Delta\kappa_\lambda,
\qquad
\Delta J_\nu,
\qquad
\Delta F_\lambda,
\qquad
\Delta W_\lambda.
\]

A generic relative error can be defined as

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

The goal is to establish quantitative statements such as

\[
G<0.02,
\qquad
O<0.1
\]

implying a small error, while regions with large overlap or sharp gradients produce substantially larger discrepancies.

The numerical boundaries must be determined by the calculations rather than assumed.

---

# 16. Data Availability by Research Phase

| Needed Input | Availability | Blocker? |
|---|---|---|
| \(E_l,E_u,J,\lambda,gf,A\) for La–Yb II/III | Excellent; public calibrated data | No |
| Broad I–IV line data for r-process elements | Public theoretical databases | No |
| Exact experimental wavelengths for millions of lines | Not available | No; quantify uncertainty |
| Analytic homologous ejecta | Generate internally | No |
| Public BNS ejecta \(Y_e,v,s,\tau\) | Available | No |
| Ready-made late-time \(\rho,T,X(v,\theta)\) RT grids | Not broadly verified | No for Paper I |
| LTE populations | Calculate internally or through transport code | No |
| Electron collision strengths for all r-process ions | Incomplete | Yes for ambitious NLTE |
| Photoionization cross sections for all r-process ions | Incomplete | Yes for ambitious NLTE |
| Accurate recombination rates for all heavy ions | Incomplete | Yes for ambitious NLTE |
| Full non-Sobolev transport engine | Existing resolved-line capability in SEDONA | No |
| Independent validation solver | Must be developed | Small–medium effort |

---

# 17. What Should Not Be Built Initially

The first project should **not** require development of:

- a new 3-D Monte Carlo engine;
- a nuclear reaction network;
- a new atomic-structure code;
- a complete NLTE solver;
- a hydrodynamic merger simulation;
- a full opacity-table infrastructure;
- a complete kilonova light-curve framework.

Existing tools already provide substantial parts of this ecosystem.

The research contribution should focus on **testing the approximation itself**, not rebuilding all of kilonova modeling.

---

# 18. Recommended Scope for Paper I

A strong first paper could use:

\[
\boxed{
\begin{array}{l}
\text{1-D homologous expansion} \\
\text{LTE populations} \\
\text{calibrated La–Yb II/III atomic data} \\
\text{selected wavelength windows} \\
\text{resolved Voigt-profile transport} \\
\text{Sobolev transport} \\
\text{small independent formal solver} \\
\text{analytic + realistic ejecta parameter ranges}
\end{array}
}
\]

This is sufficient to answer the central methodological question without introducing the major uncertainties of NLTE modeling.

---

# 19. Topics to Postpone

## 19.1 Full NLTE

Postpone because collision, photoionization, and recombination data remain incomplete for many heavy ions.

## 19.2 Full UV-to-IR profile-resolved transport

Postpone because

\[
\delta\nu/\nu\sim10^{-6}
\]

frequency resolution over the entire spectrum is computationally expensive.

## 19.3 Full 3-D realistic ejecta

Postpone because geometry and composition structure add complexity before the basic approximation has been quantified.

---

# 20. Revised Feasibility Assessment

The project is more feasible than it first appears.

The major software work is approximately

\[
\boxed{
\begin{aligned}
&\text{atomic-data integration} \\
+&\text{controlled model generation} \\
+&\text{resolved-vs-Sobolev experiments} \\
+&\text{small independent transfer solver} \\
+&\text{validity-map analysis}.
\end{aligned}
}
\]

The basic transport engine and large atomic line lists already exist publicly.

The principal missing resources arise mainly for a later NLTE extension.

---

# 21. Recommended Technical Due Diligence Before Major Development

Before committing substantial implementation time, inspect the resolved bound-bound implementation in SEDONA at source-code level.

The main questions are:

1. How are line opacities transformed between comoving and lab frames?
2. How are resolved line absorption and re-emission handled?
3. Can resolved and Sobolev modes use genuinely identical source-function assumptions?
4. How are Voigt profiles truncated?
5. How is frequency-grid resolution controlled?
6. How does memory scale with
   \[
   N_\nu N_{\rm zone}?
   \]
7. Can the same exact transition list be forced into both modes?
8. Can Monte Carlo noise be reduced enough for percent-level comparison?
9. Are overlapping resolved line profiles added directly before transport?
10. Which parts of the resolved implementation were designed for static/stellar applications versus homologously expanding transient ejecta?

The answers determine whether SEDONA can be used essentially unchanged or requires a small patch.

---

# 22. Recommended Research Strategy

The practical sequence should be:

## Stage 1 — Infrastructure

- obtain and parse calibrated lanthanide atomic data;
- reproduce standard Sobolev optical depths;
- build controlled homologous ejecta models;
- inspect and validate SEDONA's two line-opacity modes.

## Stage 2 — Minimal numerical validation

- single resolved line;
- two isolated lines;
- two overlapping lines;
- small synthetic line forest;
- compare with independent deterministic transfer solver.

## Stage 3 — Real atomic line forests

- selected Nd, Ce, Sm, or similar ions;
- selected optical, NIR, and UV windows;
- controlled variation of density, temperature, epoch, and line width.

## Stage 4 — Validity map

Calculate

\[
(G,O,T_{\rm trap})
\rightarrow
\Delta_{\rm Sob}.
\]

## Stage 5 — Kilonova interpretation

Overlay realistic kilonova conditions derived from:

- parameterized ejecta models,
- public merger ejecta distributions,
- published temperature/density histories.

## Stage 6 — Optional extensions

Only if the preceding stages show an important effect:

- wider wavelength coverage;
- fluorescence;
- NLTE populations;
- realistic multidimensional ejecta;
- full spectra or light curves.

---

# 23. Central Research Question

The project should ultimately answer:

\[
\boxed{
\begin{gathered}
\textbf{How large is the actual error introduced by the} \\
\textbf{Sobolev approximation itself in r-process ejecta,} \\
\textbf{and where in kilonova parameter space does that error matter?}
\end{gathered}
}
\]

The most efficient route is not to rebuild kilonova radiative transfer from scratch.

It is to combine:

- existing resolved-line transport,
- existing Sobolev transport,
- modern public r-process atomic data,
- a small independent reference solver,
- carefully controlled experiments,
- realistic kilonova parameter mapping.

That scope is technically manageable and scientifically clean enough to support a serious journal paper.
