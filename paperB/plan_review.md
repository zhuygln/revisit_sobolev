# Paper B — PI decisions, verbatim

The rule of this file is the same as `paper4/plan_review.md`: the PI's
decisions are recorded in their own words, in order, before they are executed.

## The research program (2026-09-21)

Yes. Returning to the **research track**, I would now treat the education project purely as preparation and keep the scientific program focused on **Paper B / the PRL attempt**.

The core research question remains: **Can detailed energy-conserving lanthanide fluorescence be replaced by a compact effective redistribution operator R_ij,** while holding the rest of radiative transfer—especially resolved Sobolev opacity—fixed?

The decisive comparison is **ε\* vs R_ij vs full energy-conserving macroatom.**

The optional neural surrogate, R_φ(j | i, θ), stays in the research plan, but only becomes important if a fixed or interpolated R_ij cannot handle state dependence.

At this point I would organize the research into four scientific questions.

1. **Does the compression phenomenon survive correct energy physics?** Redo the old promising Sobolev-only redistribution experiment using indivisible energy packets and the corrected macroatom. Start with La II, Ce II, and Nd II. Sweep roughly N_g = 2, 4, 8, 16, 32. This is the first kill gate. If a small matrix no longer works, the PRL idea largely ends here.

2. **Why does it work?** If a small R succeeds, determine whether the compression is genuinely organized by frequency locality rather than merely by parameter count. Compare contiguous-frequency grouping against an equally constrained global representation. This is where we can turn a numerical observation into a physics statement.

3. **How universal is it?** Test three things separately: state transferability, species composability, realistic mixtures. In particular, test whether R^mix_ij = Σ_a α_{a,i} R^(a)_ij predicts an explicit multi-species macroatom result without fitting the mixture.

4. **Does it remain negligible compared with real astrophysical uncertainty?** Reuse representative nuclear-heating/composition realizations from your earlier work, but run all of them through the same modern RT. Compare σ_R from the reduced operator with σ_nuclear from genuine changes in nuclear inputs. A result like **σ_R ≪ σ_nuclear** would be a strong practical validation.

The PRL decision therefore comes fairly early. I would not start with the NN or the old 32-model ensemble. The immediate production experiment should be:

**La II / Ce II / Nd II × [ε\*, R_2, R_4, R_8, R_16, R_32, macroatom].**

Everything uses identical resolved Sobolev opacity and identical physical states.

The quantities we should preregister before running are event-level redistribution error, SED error, color error, energy closure, and complexity. The most important plot would be **transport error vs effective model complexity.**

What would make me continue toward PRL is a result resembling: ε\* still has substantial irreducible error, while R_4 or R_8 reproduces macroatom transport to ≲ 0.1 mag for multiple structurally different ions. Then we ask *why*.

What would make me stop the PRL attempt is equally clear: if R needs dozens or hundreds of groups, if it only works for one ion, or if optimally tuned ε is nearly as good, then this is probably a specialized methods result rather than a compact-emergent-physics result.

So I think the next research step should be **formalizing this first production gate in `paperB/plan.md` and `paperB/prl_gate.md`, then designing the La/Ce/Nd energy-conserving experiment in enough detail to implement it without changing the criteria after seeing the result.**

## Answers to the design questions (2026-09-21)

- **State:** P1 at 2 d, shell 28, each ion alone at its P1 number density (the published xkn secular state Paper IV used; physical, energy-conserving, comparable to F61's A2 leg).
- **ε\* rule:** the grid value minimising the mean |Δm| over the live bands against the macroatom on the same state; the same metric the R legs are judged on.
- **Budget:** 3 seeds × 3×10⁵ packets per leg, as Paper IV.
- **The radiation-field-driven macroatom:** not in the gate; one robustness run per ion after a pass.
