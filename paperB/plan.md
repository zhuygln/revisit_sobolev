# Paper B — PRL plan: observable-relevant coarse graining of lanthanide fluorescence

The PI's replacement plan of 2026-10-05, verbatim (the previous plan, the
pre-run experiment map, is in git history; the preregistration
`prl_gate.md` is frozen and is not edited to match this framing).

**Status, 2026-10-05:** G1–G3 are complete, audited, merged and frozen. The
current manuscript is `docs/paperB/manuscript.tex` / `manuscript.pdf`. The
historical discovery path remains in git, `plan_review.md`, and the frozen
preregistration `prl_gate.md`; do not edit the preregistration to match this
new paper framing. This file is now the submission map rather than the
pre-run experiment map.

## The PRL question

The broad problem is not "can a matrix approximate a macroatom?" It is:

> Which microscopic degrees of freedom in a huge atomic fluorescence network
> are actually needed to reproduce the radiation that transport makes
> observable?

The controlled system is lanthanide fluorescence in kilonova ejecta. The
opacity is held fixed at the resolved Sobolev value of every line, so the
experiment isolates the post-absorption redistribution physics.

The central PRL result is:

> The reduced representation that best reproduces microscopic lanthanide
> fluorescence is not the one that best reproduces the radiation it
> transports. A global low-rank representation can fit the microscopic
> redistribution more accurately yet predict the emergent light
> substantially worse than a frequency-local coarse graining. Frequency
> locality, not generic microscopic reconstruction error, is the
> observable-relevant structure.

The second result establishes that this is useful physics rather than an
anchor-state curiosity:

> That local structure remains compact across thermodynamic states and a
> realistic 13-ion mixture, but the effective operator is state dependent:
> its parameters, and for Ce II sometimes its required resolution, vary with
> state; one-dimensional interpolation usually works, while a coupled
> physical trajectory exposes non-separability.

These are the claims the title, abstract, first page, figures and cover
letter must make obvious.

## What is new — and what we must not claim

### The novelty claim

The paper's potentially PRL-level novelty is the controlled inversion
between microscopic fidelity and observable fidelity, together with the
physical identification of frequency locality as the relevant structure.

The decisive experiment keeps fixed:

* the atomic physics,
* the resolved Sobolev opacity,
* the exit tables,
* the physical state,
* the evaluation seeds and observables,

and changes only the reduced representation of the redistribution matrix.
At matched archetype count, the global representation is more expressive
and can fit the fine microscopic redistribution 2–4 times better, yet on the
strongest Nd II comparison it predicts the broadband radiation roughly
twenty times worse than the local representation.

This is stronger than saying merely that the operator is "compressible."

### Claims we explicitly do not make

We do not claim that:

* microscopic fidelity and macroscopic/observable fidelity can differ in
  general — goal-oriented reduction, effective theories and coarse-graining
  already embody that idea;
* redistribution matrices are new — frequency-redistribution formalisms are
  classical;
* low-rank structure in radiative transfer is new — reduced-order and
  low-rank transport methods already exist;
* fluorescence matters in kilonovae — line-by-line macroatom calculations
  already demonstrate this;
* a scalar thermalisation parameter is known a priori to be exact — it is a
  standard heuristic closure, and G1 measures its failure quantitatively;
* millions of transitions "collapse to four numbers" — the conditional exit
  spectrum still contains real line structure, and Ce II can require 32
  groups;
* one universal operator exists — G3 rules that out;
* local coarse-graining is universally superior to low rank in every
  transport problem — the result is established for the physical
  redistribution problem tested here.

The safe novelty statement is:

> In lanthanide fluorescence transport, generic fidelity to the microscopic
> redistribution does not identify the degrees of freedom that control the
> emergent radiation; a frequency-local representation does.

### Literature position

The Introduction should position the result against four neighboring
literatures rather than imply that the whole coarse-graining idea is new.

1. Fontes et al. (2020), MNRAS, line-binned treatment of opacities. This is
   the direct benchmark precedent for coarse opacity transport. Their
   comparison is the starting point, not the result of this paper: it asks
   whether opacity may be frequency-grouped under their redistribution
   assumptions. We instead hold the opacity resolved and ask what can be
   discarded from the fluorescence response itself.
2. Shingles et al. (2023), line-by-line ARTIS kilonova transport. Their
   tens-of-millions-of-lines calculation is the clearest motivation for
   retaining transition identity: line-by-line treatment enables
   fluorescence and element-resolved interactions. Our question begins one
   step later: once absorption is treated line by line, must the full
   post-absorption cascade also be followed to reproduce broadband
   observables?
3. Macroatom / redistribution theory and scalar thermalisation closures.
   Lucy-style macroatoms provide the detailed stochastic reference. Reviews
   of radiative transfer make clear why constant thermalisation parameters
   are heuristic substitutes for wavelength redistribution. G1 is therefore
   a quantitative controlled comparison against a real approximation, but
   it is setup rather than the headline novelty.
4. Reduced-order / low-rank radiative transfer and goal-oriented model
   reduction. These fields establish that high-dimensional transport can
   contain low-dimensional structure and that the "best" reduced model may
   depend on the quantity of interest. Our result is narrower and physical:
   using the same fluorescence events, the globally better microscopic
   reconstruction can be the worse transport model, and the data identify
   the missing organizing principle as frequency locality.

A final ADS/arXiv novelty audit is required immediately before submission,
with searches centered on: "fluorescence redistribution operator radiative
transfer", "low-rank fluorescence macroatom", "observable-aware radiative
transfer reduction", "frequency-local redistribution coarse graining", and
"kilonova macroatom reduced model".

## The evidence: G1 is setup, G2 is the discovery, G3 is the generality proof

### G1 — establish the problem and the compact alternative

Question: with correct energy-conserving fluorescence and resolved opacity,
does a compact frequency-group operator reproduce the macroatom where a
scalar closure does not?

Frozen result: yes.

* The scalar closure misses the dense-ion references by order one
  magnitude.
* The group operator reaches the 0.10 mag band/colour criterion at small
  group count: La II and Nd II at 2 groups, Ce II at 16 at the anchor state.
* A radiation-field-driven macroatom moves the reference substantially (up
  to about 2 mag on Ce II); an operator rebuilt from that reference follows
  it to a few hundredths of a magnitude.
* Therefore the compression is not an accident of the downward macroatom;
  the operator carries the fluorescence physics it was built from.

Role in the PRL: establish that there is a nontrivial reduced model to
explain. Do not spend the paper's conceptual budget here.

### G2 — the PRL centerpiece

Question: why does the compact operator work? Is the microscopic
redistribution globally low rank, or does transport care about local
frequency organization?

Frozen result: frequency locality is the relevant structure.

* Compare local contiguous frequency blocks with global non-negative
  factorisation at matched numbers of archetypal exit distributions.
* The global family is more expressive and uses more free numbers.
* On the decisive ions it needs no fewer archetypes to reach the
  photometric criterion; on Nd II, two local blocks pass while global rank
  32 does not.
* At 32 archetypes, the global representation fits the microscopic
  redistribution 2–4 times better while reproducing the transported light
  2–20 times worse.
* Therefore microscopic reconstruction error is not the correct objective
  for this closure.

Role in the PRL: this is the discovery. Figure 2 is the conceptual center of
the paper and should be the figure an editor can describe from memory.

The section heading should emphasize the broad result first, e.g.

> Observable fidelity is not microscopic fidelity

with frequency locality presented as the physical answer.

### G3 — show that G2 is physical and useful, not a one-state artifact

Question: does the frequency-local effective description remain compact
when state and composition change, and can it be predicted away from the
training state?

Frozen result: compact but state dependent.

* A fixed anchor operator generally fails when gas temperature, density, or
  the physical trajectory changes.
* A fresh compact operator survives across the tested domain; failure of
  transfer is not failure of compression.
* Whole-operator interpolation between bracketing states succeeds at every
  tested single-coordinate interior point.
* Nd II on the coupled physical trajectory is the named boundary: transfer
  and interpolation fail while a fresh compact fit succeeds. The manifold is
  compact but not guaranteed to be separable in independently varied state
  coordinates.
* Ce II's 16-group representation fails at 7/14 tested states but all are
  recovered at 32 groups. Compactness survives; fixed representation size
  does not.
* Simple composition-weighted mixing of separately trained species fails,
  but a directly trained blend operator succeeds.
* On the full 13-ion realistic P1 blend, the directly trained operator
  passes at 4 groups (0.088 mag band, 0.073 mag colour), while the best
  scalar closure misses by 1.69 mag.

Role in the PRL: generality proof. G3 should support G2, not replace it as
the headline.

## The paper argument

The Letter should read as a four-step physical argument, not as a
chronology of gates.

### Step 1 — The microscopic problem is genuinely hard

Lanthanide fluorescence couples an enormous transition network, and
existing high-fidelity transport retains transition identity because
fluorescence depends on which line absorbed. A scalar thermalisation
closure is therefore not guaranteed to preserve the wavelength
redistribution.

Use G1 to show the consequence: the scalar fails while a small structured
operator succeeds.

### Step 2 — The surprising result: the microscopically better model is worse

Introduce the matched-complexity local-versus-global experiment.

Lead with the inversion:

    E_micro^global < E_micro^local,    E_obs^global > E_obs^local.

This is the paper's PRL moment.

Then identify the physical reason: neighboring absorbed frequencies share
transport-relevant response in a way a generic global factorisation does
not respect.

### Step 3 — What transport actually remembers

State the two-level compression carefully:

    P(nu_out | nu_in) = P(j | i) P(nu_out | i, j).

The first factor is extremely coarse; the second still needs a sparse set
of real exit lines. The result is observable-relevant compression, not a
claim that the microscopic atomic network literally contains only a handful
of degrees of freedom.

### Step 4 — The effective object is a state-dependent family

G3 changes the interpretation from one universal operator to

    R = R(theta),

with adaptive resolution where necessary. One-coordinate interpolation
usually works; the coupled trajectory supplies an explicit boundary; the
realistic mixture remains strongly compact.

End on the 13-ion result, because it connects the conceptual finding back
to the physical transport problem.

## Main-paper figures

Keep three figures in the PRL core.

**Figure 1 — existence of the reduced description.** Scalar versus group
operator on La/Ce/Nd, plus the stronger-macroatom robustness panel.
Message: a scalar fails; a compact structured operator follows very
different fluorescence physics.

**Figure 2 — centerpiece: microscopic versus observable fidelity.** Local
frequency coarse-graining versus global NMF: photometric error, microscopic
event error, and the cross-plot. Caption headline: *Better on the
microscopic events, worse on the light.* This is the figure most likely to
carry the PRL decision.

**Figure 3 — generality and boundary.** Transfer versus fresh fit versus
whole-operator interpolation across state axes, with the coupled trajectory
visible. The 13-ion mixture result should either be a compact inset/panel
here or a single sentence plus End Matter table if space is tight.

## What moves out of the PRL core

Current APS guidance limits a PRL Letter to 3750 core words, with figures
and tables carrying word equivalents; up to two pages of End Matter may be
added for specialist detail. The submission should therefore move the
following out of the core:

* full gate tables;
* most per-state G3 numbers;
* implementation-defect history and audit details;
* the full preregistration/gray-condition logic;
* convergence tables;
* exit-table audit details;
* exhaustive per-ion metrics;
* reproducibility machinery beyond one concise paragraph.

Keep enough in the main text for the physics to be independently
understandable. Put the complete audit trail, preregistration and generated
tables in End Matter / Supplemental Material and the public repository.

## Immediate manuscript revisions

Before submission:

1. Make G2 the novelty claim in the title/abstract/introduction. G1 is
   setup; G3 is validation/generality.
2. Change the G2 section heading toward "Observable fidelity is not
   microscopic fidelity"; keep "frequency locality" as the answer.
3. Rewrite the first two Introduction paragraphs so a non-kilonova physicist
   sees the general coarse-graining problem before R_ij.
4. Do not claim the broad principle itself as new. Cite goal-oriented /
   reduced-order context and state that the new result is its controlled
   realization and physical resolution in lanthanide fluorescence.
5. Clarify macroatom language. Do not say every code replaces fluorescence
   by a closure and then list the explicit macroatom as a closure of the
   same kind. Distinguish stochastic detailed-network treatment from
   parameterized/coarse redistribution.
6. Replace internal "Paper IV" dependency. Before submission either post
   Paper IV as a citable preprint and cite it, or make the small amount of
   transport validation Paper B needs self-contained in End
   Matter/Supplement.
7. Fix the F70/report inconsistency about source-spectrum transfer. The
   frozen Ce record and manuscript say Ce fails transfer on every axis,
   while one sentence in `docs/results_report.md` still says J holds for
   every ion. The frozen record is authoritative.
8. Trim the Method aggressively. The main text needs the physical state,
   reference, local/global constructions, independent build/evaluation
   seeds, observables and 0.10-mag criterion. Run-history details go out.
9. End the Letter on the physical statement, not on implementation:
   transport can discard most microscopic fluorescence information, but
   only if the reduced coordinates respect the structure the radiation
   actually samples.

## Candidate title and editorial pitch

Current title remains strong:

> Effective line redistribution in kilonova ejecta: transport-relevant
> structure is frequency-local

A broader alternative, only if the manuscript itself fully supports the
framing:

> Observable-relevant coarse graining of lanthanide fluorescence

The one-sentence editorial pitch should be:

> We find that the reduced model that best reproduces microscopic
> lanthanide fluorescence is not the model that best reproduces the light: a
> globally more accurate low-rank representation gives substantially worse
> transported observables than a frequency-local coarse graining, which
> remains compact across state changes and a realistic 13-ion mixture.

## G4 and new experiments

The original G4 compared operator error with a composition/nuclear-
realisation proxy. It is not a prerequisite for this PRL and should remain
demoted. It does not strengthen the central G2 novelty enough to justify
delaying the Letter, and the repository does not contain a genuine
nuclear-realisation ensemble.

Do not start a neural surrogate before submission. G3 contains one failure
of separable interpolation, not evidence that a learned model is needed.

One follow-up experiment could raise the astrophysical ceiling beyond this
Letter:

> Replace the detailed macroatom by the state-dependent effective operator
> in the existing multi-shell, time-dependent P1/xkn light-curve
> calculation, first with locally rebuilt blend operators and then with a
> sparse predictive table.

That experiment is valuable for a future Nature Astronomy / full transport
paper, but it is not required to establish the present PRL's conceptual
claim. Do not let it delay submission unless the PRL draft itself reveals
that the single-zone limitation prevents the central argument from reading
as complete.

## Submission decision rule

Submit to PRL when all of the following are true:

* Figure 2 is visually and textually the conceptual center.
* The first page explains the coarse-graining question to a broad physicist
  without requiring kilonova-specialist knowledge.
* The novelty statement is literature-safe: no claim that observable-aware
  reduction, redistribution matrices, low rank, or fluorescence itself is
  new.
* G3 is presented as evidence that the G2 mechanism survives beyond one
  anchor state, with the Nd/Ce limitations stated rather than hidden.
* The 13-ion result is visible in the main Letter.
* Paper IV is either citable or Paper B is self-contained.
* The core fits the PRL length budget without moving essential physics into
  the supplement.
* A final ADS/arXiv search finds no paper that already performs the same
  local-versus-global, microscopic-versus-observable inversion experiment.

If PRL declines on breadth rather than correctness, the same manuscript
should transfer cleanly to a specialist physics/astrophysics venue with the
full tables restored. The PRL attempt should therefore optimize the framing
for broad physical significance without changing any frozen scientific
result.
