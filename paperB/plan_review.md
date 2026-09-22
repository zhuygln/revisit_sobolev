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

## The seven pre-run changes to G1 (2026-09-22, verbatim; five blocking)

This is much stronger. I think G1 is now close to being a defensible preregistration rather than just an execution plan. I would **not run production yet**, though. I see **five changes I would treat as blocking**, plus two smaller clarifications. They are easiest to fix now, before any results exist.

The overall experimental logic is excellent: same resolved opacity + {ε\*, R_N, energy-conserving macroatom} on three physically motivated single-ion states, with a genuine kill gate. That is exactly the right first experiment.

**1. "Full macroatom" is too strong for R2.** This is the most important wording issue. Your actual reference is `sobolev_dmacro` — the **downward energy-conserving macroatom**. It is not the full radiation-field-coupled macroatom, and certainly not a self-consistent NLTE macroatom. So I would change everywhere in G1 "full energy-conserving macroatom" to "**energy-conserving downward macroatom reference**" or simply "**detailed energy-conserving reference**". Then your post-pass robustness test is correctly: R2M = radiation-field-driven macroatom with upward transitions under imposed W = ½ B_ν. That gives the scientific hierarchy: G1 reference: downward macroatom → post-pass robustness: radiation-field-driven macroatom. If eventually we want to write "full macroatomic physics" in a PRL, we'll need stronger validation later. But G1 does not need that to answer its narrowly defined question.

**2. The kernel must have a build/validation split.** Right now the R_N matrix and your `m_event` comparison are constructed from **the same R2 events**. That makes `m_event` partly an in-sample reconstruction metric. For a matrix this isn't "ML overfitting" in the usual sense, but a referee could still reasonably say: "Of course the empirical operator resembles the events from which it was constructed. Does it reproduce independent macroatom events?" I would fix this now. The cleanest scheme is: kernel-build macroatom events ≠ kernel-validation macroatom events, while keeping the transport evaluation independent as well. For example: build seeds 101, 102, 103; evaluation/reference seeds 1, 2, 3. Then R_N ← macroatom events from build seeds, but m_event is evaluated against the independent R2 event distribution from seeds 1–3. And transport with R_N also uses seeds 1–3. This makes the interpretation much cleaner: construct effective operator once → predict independent macroatom statistics and transport. If the extra macroatom run is too expensive, three-fold cross-fitting is another option, but dedicated build/evaluation samples are simpler to explain. This also needs to specify whether one kernel is pooled across build seeds or one kernel per seed. I strongly prefer **one pooled kernel per ion/state/N**, followed by independent transport seeds.

**3. `m_sed` needs the spectral-bin widths.** You currently define m_sed = Σ_b |L_ν,b − L_ν,b^R2| / Σ_b L_ν,b^R2. If these are actual spectral densities L_ν, that is only an integral if all bins have equal Δν. Safer definition: m_sed = Σ_b |L_ν,b − L_ν,b^R2| Δν_b / Σ_b L_ν,b^R2 Δν_b, or the corresponding L_λ Δλ version. If the stored quantity is already integrated luminosity per bin, then your original sum is correct — but state that explicitly. This sounds minor but is precisely the kind of definition a methods referee may inspect.

**4. `m_event` should be energy-weighted, not event-count-weighted.** You write "L1 per row, weighted by the row's event count." But these are energy-packet kernels, and `from_branching_mc` receives w_in, w_out. If incoming events carry different energy weights, raw event counts are no longer the correct physical measure. I would define d_i = ½ Σ_j |R^(N)_ij − R^(128)_ij| and then m_event = Σ_i W_i d_i / Σ_i W_i, where W_i = Σ_{e∈i} w_in,e. Then m_event ∈ [0, 1] if the factor ½ is used. That has a clean interpretation: expected total-variation loss in the outgoing redistribution distribution for an energy-weighted incoming interaction. That's much better than "summed row L1."

**5. N_g² is not necessarily the true model complexity.** This may ultimately be the most important scientific issue. Your successful kernel isn't just R_ij. It also contains **discrete within-group exit tables**. Those tables were scientifically important: they fixed the self-absorption artifact in F25. So if an N_g = 4 kernel consists of 4×4 transition probabilities plus thousands of stored microscopic output frequencies, we cannot say: "We compressed millions of transitions into 16 numbers." The operator actually contains more information. You're already recording serialized table size, which is excellent. I would make this explicit now. Call N_g² the **coarse transition-matrix degrees of freedom**, not "model complexity" without qualification. Also record N_exit samples, serialized bytes, and ideally compression ratio relative to the event corpus/reference atomic representation. The central figure could still use N_g², because it is intuitively useful, but label it precisely: "**Error vs coarse matrix size N_g²**" with a companion axis/table giving actual serialized kernel size. This becomes especially important if we're eventually claiming an emergent compact representation. We need to establish that the whole representation is compact — not merely the transition-probability part. Longer term, if the discrete exit tables dominate storage, that may motivate another compression step (P(ν_out | i, j) as quantiles, a small histogram, or another compact conditional representation). But don't add that experiment to G1 yet. Just measure honestly what the current kernel stores.

**6. I would weaken B3's Red condition.** This currently worries me: "Red: m_band rises with N_g for any ion by more than its noise." I would **not infer a transport pathology merely from non-monotonic broadband error**. Even if the underlying operator converges monotonically, a broadband difference can behave like 0.08, 0.04, 0.05, 0.02 because different spectral errors cancel differently when integrated through filters. Transport is nonlinear enough that strict monotonicity of every observable is not guaranteed. I would instead make m_event convergence the strongest structural diagnostic; m_sed convergence strongly expected; m_band non-monotonicity **Yellow**, unless accompanied by an identified pathology such as the old within-bin self-absorption behavior. For example: B3 Green: overall convergence trend in all three; no statistically resolved degradation between N_g = 8 and 32. B3 Yellow: observable metrics show statistically significant local non-monotonicity. B3 Red: higher resolution produces a persistent/worsening discrepancy traceable to a representation pathology — e.g. repeated self-absorption, excessive empty-row fallback, or an event-level metric that does not converge. This prevents an innocent 0.01-mag wiggle from killing the PRL path. Because you have not run production yet, this is exactly the right time to change this preregistration criterion.

**7. Rename G4's σ_nuclear.** Your honesty in the plan is good: there is no nuclear-input ensemble in the repository. Therefore I would **not call the existing quantity** σ_nuclear. The four available composition patterns are not the nuclear-model ensemble from your older work. Call it something like **σ_comp** or σ_composition proxy. Then G4 becomes initially σ_R ≪ σ_comp. And reserve σ_nuclear for the future experiment using the actual external nuclear-realisation ensemble. This avoids a subtle but real overclaim later.

**One smaller concern: the live-band rule.** The current rule combines finite + L_b/L_bol ≥ 1 % + 40-Mpc detection threshold. The first two are natural numerical/physical viability criteria. The 40-Mpc magnitude threshold is an **observational selection criterion** imported from Paper IV. For this single-ion microphysics experiment, I am not sure it belongs in the primary gate. It can make the set of bands — and therefore ε\* — depend on an arbitrary assumed distance/survey depth. I slightly prefer finite + L_b/L_bol ≥ 1 % + acceptable MC precision for G1, with the 40-Mpc detectability cut reported separately. This isn't as serious as points 1–5, but I would decide before running.

**What I think is especially good.** The physical state is excellent: t = 2 d, T = 3401 K, ρ = 3.298×10⁻¹⁵ g cm⁻³, with each ion at its physical P1 abundance-derived density, rather than τ_max-normalized. The scalar baseline is unusually fair: ε\* = argmin over the grid of the mean live-band |Δm_b|. If R_4/R_8 substantially beats that oracle-tuned scalar, the result means something. I also like that B2 requires both e_ε ≥ 0.20 mag and e_ε ≥ 3 e_R. That protects you against claiming victory because R improves an already negligible discrepancy. The Gray logic is also much better than the usual "run until something looks clean." In particular, energy closure, trapped packets, empty rows and seed noise all get checked before scientific interpretation. And this is exactly the right PRL kill logic: B1 = does compact R work? B2 = is scalar ε genuinely insufficient? Only after both succeed do we spend effort asking *why* in G2.

**My recommended pre-run status.** With the changes above, I would freeze G1. I would consider 1–5 blocking before production: rename the reference accurately; separate kernel construction from validation; fix/clarify the integrated SED metric; energy-weight m_event; account honestly for within-cell exit-table complexity. Then I would revise B3 as discussed and rename the G4 proxy before merging the plan. After that, I would be comfortable with: freeze preregistration → run G1. And importantly, I would **not add any more physics to G1**. No matched-τ control, no NN, no extra ions, no full mixture, no extra epochs. Those are later questions. G1 should remain brutally simple: La II, Ce II, Nd II: ε\* vs R_N vs downward energy-conserving macroatom — and either that compression phenomenon survives correct energy physics or it doesn't.

## Freeze and run (2026-09-22)

PI: "run it". The preregistration `paperB/prl_gate.md` is frozen at commit
f9cfbff (the seven changes applied); the production run of G1 starts from
that commit, three ions sequentially, `paperB/gate1/run_<ion>.log`.


## G1 read (2026-09-22)

At 3×10⁵ packets Nd II fired gray condition 2 (g's seed scatter 0.063); the
preregistered remedy (the whole Nd grid at 10⁶ packets per seed) was run,
the 3×10⁵ record kept, and the gate read on the three ions: B1 Green, B2
Green, B3 Yellow → CONTINUE. Nothing in `prl_gate.md` was edited; the
analysis was changed only to accept a raised (never a lowered) packet
count and to read the highest-count record per ion.
