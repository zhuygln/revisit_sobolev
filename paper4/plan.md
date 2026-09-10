# Paper IV plan — energy-conserving closure test (as received 2026-09-09)

Kept verbatim as the historical input. The review of the implementation plan
that followed the same day is [`plan_review.md`](plan_review.md); the
implementation plan itself, as approved, is summarised in
[`README.md`](README.md).

---

Yes. I would now reorganize the project around one question:

$$
\boxed{
\textbf{Does the opacity-closure effect survive when the radiation-energy problem is made physically self-consistent?}
}
$$

The current 27-point grid has already done its job as a controlled mechanism experiment. I would **not spend more time extending that grid yet**. The next work should progressively remove the known artificial assumptions, with hard stop/go gates so we do not spend months building a full transport code around an effect that might disappear.

# Phase 0 — Freeze the current experiment as the controlled baseline

Before changing physics, freeze the present state.

Call the existing calculation something like:

> **Controlled fixed-atmosphere closure experiment**

and stop calling its parameter grid "realistic kilonova models."

Preserve the current R/A/B/C outputs, F43–F48, all provenance, and the current equal-four-ion results. Rename the composition parameter in those runs from \(X_{\rm lan}\) conceptually to something like

$$
X_{\rm 4Ln}
=
X_{\rm La}+X_{\rm Ce}+X_{\rm Pr}+X_{\rm Nd}.
$$

That prevents future results from getting mixed with the physically meaningful total lanthanide fraction.

No new MC simulations are needed for this phase.

---

# Phase 1 — Build two physically coherent benchmark ejecta states

Do **not** begin with another Cartesian grid.

First build only two or three benchmark models where all physical quantities come from a coherent source.

I would use two complementary anchors.

| Benchmark       | Purpose                       | Physical character                |
| --------------- | ----------------------------- | --------------------------------- |
| **P1**          | theoretical neutron-rich test | lanthanide-rich merger ejecta     |
| **P2**          | observationally anchored test | AT2017gfo-like photospheric state |
| **P3 optional** | lanthanide-poor control       | higher-\(Y_e\) ejecta             |

For P1, use a published neutron-rich trajectory around \(Y_e\sim0.20\), with consistent mass, velocity/density profile, nucleosynthetic abundances, heating, and lanthanide fraction.

For P2, use a published AT2017gfo-like atmosphere at a specific epoch—roughly 2–4 d is ideal—so that density profile, inner velocity, temperature and composition are already constrained by a real-event model.

The crucial change is that we no longer independently select

$$
M,\quad v,\quad X_{\rm lan},\quad T,\quad \kappa.
$$

Instead, each benchmark is a single coherent physical state.

### Composition

Replace equal La/Ce/Pr/Nd with the full feasible lanthanide pattern:

$$
X_i = w_i X_{\rm lan}.
$$

Use published r-process \(w_i\), ideally for all available La–Yb elements.

For the first comparison, however, keep ionization fixed if necessary so we can distinguish the effect of changing abundances from the effect of changing ionization.

Then add II/III ionization as a separate controlled step.

### Density

Replace the uniform sphere with a radial profile such as

$$
\rho(v,t)=\rho_0(t)f(v),
$$

with \(f(v)\) taken from the benchmark model.

### Output of Phase 1

For every radial shell we should have a machine-readable state:

$$
\{t,r,v,\rho,T_{\rm gas},T_{\rm rad},X_i,n_{i,j},Y_e\}.
$$

This becomes the input to all subsequent RT calculations.

### Gate 1

Before doing transport, verify:

$$
\int 4\pi r^2\rho\,dr=M_{\rm ej}
$$

and that the abundance fractions sum correctly.

Also plot:

$$
\rho(v),\quad T(v),\quad X_i(v),\quad n_{\rm ion}(v).
$$

If those do not resemble the source model, stop.

---

# Phase 2 — Separate energy bookkeeping from opacity closure

This is the most important next technical experiment.

Do **not** immediately implement full radiative equilibrium.

First make the **transport itself energy conserving while keeping the atmosphere fixed**.

We need four calculations:

$$
\begin{array}{ll}
R_1:&\text{resolved Sobolev + current photon-number branching}\\
R_2:&\text{resolved Sobolev + energy-conserving atomic transport}\\
C_1:&\text{grouped opacity + current branching}\\
C_2:&\text{grouped opacity + energy-conserving atomic transport}.
\end{array}
$$

Then:

$$
R_2-R_1
$$

tells us how much of the present result comes from the simplified energy treatment.

And:

$$
C_2-R_2
$$

is the opacity-closure error we actually care about.

That decomposition is crucial.

---

# Phase 2A — Implement indivisible energy packets

Each packet should carry an energy

$$
E_p
$$

rather than simply representing one photon.

A line interaction should never destroy energy.

After absorption, the energy can move between internal atomic states and radiation, but the bookkeeping must satisfy

$$
E_{\rm before}
=
E_{\rm radiative,out}
+
E_{\rm thermal}
+
E_{\rm internal}
$$

to numerical precision.

For a Monte Carlo realization, track globally:

$$
E_{\rm source},
$$

$$
E_{\rm escaped},
$$

$$
E_{\rm deposited},
$$

$$
E_{\rm stored},
$$

and, once included,

$$
E_{\rm adiabatic}.
$$

Initially, with a fixed atmosphere and no thermal feedback, it is enough to enforce

$$
E_{\rm source}
\simeq
E_{\rm escaped}+E_{\rm deposited}.
$$

### Proposed acceptance threshold

I would predeclare:

$$
\frac{|E_{\rm source}-E_{\rm accounted}|}
     {E_{\rm source}}
<10^{-3}
$$

for unit tests and preferably \(<1\%\) in large MC production runs.

No grey renormalization should be necessary for the energy-conserving runs.

---

# Phase 2B — Validate the atomic energy machinery before using lanthanides

Build very small systems first.

### Test A — two-level atom

Known excitation energy:

$$
E_2-E_1=h\nu.
$$

Absorb one energy packet, emit one equal-energy transition.

Verify exact energy conservation.

### Test B — three-level fluorescent atom

For example:

$$
3\rightarrow2\rightarrow1.
$$

Verify:

$$
h\nu_{31}
=
h\nu_{32}+h\nu_{21}.
$$

The MC branching probabilities should reproduce analytic branching ratios.

### Test C — thermal channel

Allow a fraction to thermalize.

Verify:

$$
E_{\rm in}
=
E_{\rm radiation}+E_{\rm thermal}.
$$

### Test D — optically thick repeat-interaction limit

Make sure repeated fluorescence does not create or destroy energy.

This should replace the present ambiguous chain-cap behavior with explicit energy bookkeeping.

Do not move to the lanthanide network until all four pass.

---

# Phase 3 — Repeat the headline differential on one realistic benchmark

Now take **P1 only**.

Do not run a grid.

Run:

$$
R_2,\quad A_2,\quad B_2,\quad C_2.
$$

Where:

* \(R_2\): resolved opacity + detailed energy-conserving branching;
* \(A_2\): resolved opacity + compressed redistribution;
* \(B_2\): grouped opacity + detailed energy-conserving branching;
* \(C_2\): grouped opacity + compressed redistribution.

Then evaluate:

$$
\Delta m_b(t)
=
m_{b,\mathrm{closure}}
-
m_{b,\mathrm{reference}}.
$$

The most important question is **not whether it remains exactly 2.84 mag**.

The important hierarchy is:

$$
|\Delta m_A|
\ll
|\Delta m_B|,
|\Delta m_C|.
$$

And the qualitative signature:

$$
g/r\ {\rm brighter}
$$

versus

$$
J/H/K\ {\rm fainter}
$$

should persist if our interpretation is correct.

### Gate 2

I would use three possible outcomes.

**Green:** The chromatic effect remains large and coherent.

For example, multiple bands remain wrong by several tenths of a magnitude or more and the A control remains much smaller.

Then proceed.

**Yellow:** Effect remains but falls substantially.

Example:

$$
2\ {\rm mag}\rightarrow0.2{-}0.5\ {\rm mag}.
$$

Then the qualitative mechanism survives, but the paper pivots away from the old 1–3 mag headline.

**Red:** \(R_1-R_2\) accounts for most of the old result and \(R_2-C_2\) becomes very small.

Then stop the current astrophysical interpretation.

That would mean the old result was substantially driven by the energy treatment.

---

# Phase 4 — Reproduce the Fontes-like limit

I think this is the second-most-important physics experiment.

We need to demonstrate **why previous line-binned calculations looked successful**.

Use the same benchmark and same atomic data.

Run:

$$
R_{\rm thermal}
=
\text{resolved opacity + complete thermal redistribution}
$$

and

$$
C_{\rm thermal}
=
\text{line-binned/grouped opacity + complete thermal redistribution}.
$$

First ask whether

$$
R_{\rm thermal}\approx C_{\rm thermal}.
$$

If yes, we reproduce the regime in which line-binned opacity appears accurate.

Then change only one thing:

$$
\text{thermal redistribution}
\rightarrow
\text{explicit energy-conserving fluorescence}.
$$

Now compare again.

If:

$$
R_{\rm fluorescence}
\not\approx
C_{\rm fluorescence},
$$

we have a very strong causal result:

$$
\boxed{
\textbf{coarse opacity may appear adequate under thermal redistribution}
\\
\textbf{but fail when the actual fluorescence network is retained.}
}
$$

That may become the conceptual center of the paper.

---

# Phase 5 — Test the Morag closure

Add another opacity treatment rather than merely comparing EP93 and \(\sum\tau\).

Implement the expansion-limited line strength idea schematically as

$$
\kappa_{l,\rm eff}
=
\min
\left[
\kappa_l,\,
\frac{1}{\rho ct}
\right].
$$

Call this treatment \(D\).

Then compare:

$$
R_2,\quad B_2,\quad C_{\rm bin,2},\quad D_2.
$$

The question is:

> Can a physically motivated intermediate closure recover most of the resolved fluorescent spectrum?

If yes, the paper becomes more constructive:

> We identify not only the failure, but the information the closure needs to preserve.

If no, then the conclusion is also interesting:

> Fixing the emissivity-rate problem alone is insufficient because fluorescence carries additional nonlocal wavelength information.

---

# Phase 6 — Numerical convergence on the realistic benchmark

Only after the energy-conserving reference exists.

Do this on P1 and perhaps P2.

Vary:

$$
\tau_{\min}
=
10^{-2},10^{-3},10^{-4}
$$

and opacity bin widths such as

$$
\Delta v
=
1.25,\ 12.5,\ 125\ {\rm km\,s^{-1}}.
$$

Also test packet number:

$$
N,\quad2N,\quad4N.
$$

Measure convergence of actual broadband quantities.

I would ask for approximately

$$
|\Delta m_{\rm convergence}|
<0.03{-}0.05\ {\rm mag}
$$

in the bands used for headline conclusions.

If that is too expensive, quote the measured numerical uncertainty rather than forcing an arbitrary threshold.

---

# Phase 7 — Add realistic ionization

Now upgrade composition from:

$$
X_i
$$

to

$$
X_{i,j}
$$

where \(j\) is the ion stage.

Start with LTE Saha ionization:

$$
\frac{n_{j+1}n_e}{n_j}
=
\frac{2U_{j+1}}{U_j}
\left(
\frac{2\pi m_ekT}{h^2}
\right)^{3/2}
e^{-\chi_j/kT}.
$$

Use at least II/III where atomic data exist.

This is still not NLTE, but it is much more physical than forcing all lanthanides to II.

For every epoch, record:

$$
f_{\rm II},\quad f_{\rm III}.
$$

Then repeat P1/P2.

### Gate 3

The broad R-vs-C conclusion should survive changing the dominant ionization state.

It need not have identical amplitude.

---

# Phase 8 — Now close the radiation–matter feedback loop

Only if the effect survives Phases 2–7 should we build the fully self-consistent thermal calculation.

At this stage the atmosphere is no longer assigned a fixed \(T_{\rm gas}\).

Instead solve:

$$
\frac{du}{dt}
=
\dot q_{\rm nuc}
+
\dot q_{\rm abs}
-
\dot q_{\rm emit}
-
P\nabla\cdot v.
$$

For a quasi-steady epoch, the first implementation can simplify this to radiative equilibrium:

$$
\dot q_{\rm heat}
\simeq
\dot q_{\rm cool}.
$$

The iterative loop becomes:

$$
T^{(0)}
$$

$$
\downarrow
$$

compute ionization/excitation

$$
\downarrow
$$

transport energy packets

$$
\downarrow
$$

estimate heating/cooling

$$
\downarrow
$$

update \(T\)

$$
\downarrow
$$

repeat.

For example:

$$
T^{(n+1)}
=
T^{(n)}
+
\alpha\,
\frac{\dot q_{\rm heat}-\dot q_{\rm cool}}
{\partial\dot q/\partial T}.
$$

Use damping \(0<\alpha<1\) initially for stability.

### Convergence criteria

I would require simultaneously:

$$
\left|
\frac{\dot q_{\rm heat}-\dot q_{\rm cool}}
{\dot q_{\rm heat}}
\right|<1\%
$$

and

$$
\frac{|T^{(n+1)}-T^{(n)}|}
{T^{(n)}}<1\%
$$

in the important line-forming shells.

And broadband magnitudes should stop moving by more than roughly

$$
0.01{-}0.03\ {\rm mag}
$$

between iterations.

---

# Phase 9 — Recompute the opacity itself after each temperature update

This is essential.

When \(T_{\rm gas}\) changes, we cannot keep the original Sobolev optical depths.

Recompute:

$$
T
\rightarrow
\text{ionization}
$$

$$
\rightarrow
\text{level populations}
$$

$$
\rightarrow
\tau_l
$$

$$
\rightarrow
\text{branching / redistribution}
$$

$$
\rightarrow
\text{radiation field}.
$$

So the full iteration is really:

$$
T
\rightarrow
n_i
\rightarrow
\tau_l
\rightarrow
J_\nu
\rightarrow
\dot q
\rightarrow
T.
$$

That is the self-consistent fixed-hydrodynamic radiation problem we actually want.

---

# Phase 10 — Repeat the reference/closure comparison after thermal convergence

This is the definitive experiment.

Do two separate self-consistent solutions:

$$
R_{\rm eq}
$$

using resolved lines, and

$$
C_{\rm eq}
$$

using grouped opacity.

Do **not** force them to share the same \(T_{\rm gas}\).

Each approximation should be allowed to find its own equilibrium.

Then compare not only spectra, but also:

$$
T_R(r)-T_C(r),
$$

$$
f_{{\rm ion},R}-f_{{\rm ion},C},
$$

$$
\dot q_{\rm dep,R}-\dot q_{\rm dep,C},
$$

and

$$
m_{b,R}-m_{b,C}.
$$

This tells us whether the opacity closure changes the **thermal state itself**.

That would be a substantially deeper result than the present fixed-atmosphere work.

---

# Phase 11 — Global inference, not tangent-space inference

There is also a cheap piece we can do earlier using the existing grid, but the final version should use the realistic models.

Instead of only asking whether

$$
d_{\rm RT}
$$

lies in a local derivative space, perform the actual inverse problem:

$$
\theta_{\rm fit}
=
\arg\min_\theta
\chi^2
[
m_{\rm closure}(\theta)
-
m_{\rm reference}(\theta_{\rm true})
].
$$

This determines whether the closure causes:

### A. Good fit, correct parameters

Not astrophysically serious.

### B. Good fit, wrong parameters

$$
\boxed{\text{parameter bias}}
$$

Potentially very important.

### C. Bad fit everywhere

$$
\boxed{\text{distinct model discrepancy}}
$$

which is what the current tangent analysis suggests.

I would run the cheap grid-based version now in parallel, but only make a strong astrophysical claim after the realistic energy-conserving models exist.

---

# Phase 12 — Independent external-code validation

Run this in parallel rather than making it a blocker for every other task.

Use a production energy-conserving macroatom code—ARTIS is the obvious candidate—to reproduce one benchmark state as closely as possible.

The most useful comparison is the resolved reference:

$$
R_{\rm ours}
\quad\mathrm{vs}\quad
R_{\rm external}.
$$

Compare:

* broad SED shape;
* optical/NIR flux ratio;
* dominant fluorescence direction;
* temperature/ionization if available.

We do not need line-by-line agreement at every Å.

We need to know whether our resolved reference is physically in the same regime as an independent implementation.

If possible, later add the grouped-opacity treatment to the external code.

---

# Phase 13 — Only then decide whether another grid is justified

If all the important gates survive, then build a **physical parameter study**.

But I would not return to an arbitrary Cartesian

$$
M\times v\times X_{\rm lan}
$$

grid.

Use trajectories in physically meaningful coordinates, such as

$$
Y_e,\quad s,\quad \tau_{\rm exp},
$$

plus ejecta component parameters.

For example:

### Lanthanide-rich trajectory

$$
Y_e\sim0.15{-}0.20.
$$

### Intermediate

$$
Y_e\sim0.25.
$$

### Lanthanide-poor

$$
Y_e\sim0.30{-}0.35.
$$

Each trajectory supplies its own:

$$
X_i,\quad
X_{\rm lan},\quad
\dot q(t),\quad
\kappa,\quad
\text{ionization environment}.
$$

That is the grid appropriate for an astrophysical paper.

---

# Recommended execution order

If I were running the project, I would sequence it as follows:

1. **Freeze/relabel current grid.**
2. **Construct P1 and P2 physical benchmarks.**
3. **Implement and validate energy packets on toy atoms.**
4. **Run \(R_1/R_2/C_1/C_2\) on P1.**
5. **Fontes thermal-redistribution bridge.**
6. **Morag closure comparison.**
7. **Line-cutoff/bin-width/packet convergence.**
8. **Realistic r-process abundance + II/III ionization.**
9. **Repeat P1 and P2.**
10. **Independent ARTIS reference comparison in parallel.**
11. **Implement radiative-equilibrium temperature iteration.**
12. **Run fully converged \(R_{\rm eq}\) vs \(C_{\rm eq}\).**
13. **Global nonlinear parameter recovery.**
14. Only if all of that survives, **build a broader physically coupled survey and rewrite the flagship paper around it.**

---

## What I would *not* do now

I would specifically avoid four temptations:

* don't rerun the old 27-point grid with only corrected La/Ce/Pr/Nd ratios;
* don't immediately implement full NLTE;
* don't add more tangent-space nuisance directions;
* don't spend weeks polishing the Nature Astronomy manuscript around the current 1–3 mag number.

Those all optimize the wrong level of the problem.

The single most important milestone now is:

$$
\boxed{
R_2-C_2
}
$$

—**resolved versus grouped line transport with exact energy bookkeeping in the same realistic ejecta state**.

If that difference survives, almost every subsequent investment is justified.

If it does not, we will have learned that before rebuilding the entire project around an artifact.

So I would define the next major project gate as:

$$
\boxed{
\textbf{Does the large coherent chromatic closure error survive}
\\
\textbf{energy-conserving fluorescence in a realistic r-process atmosphere?}
}
$$

Everything else should follow from the answer to that.
