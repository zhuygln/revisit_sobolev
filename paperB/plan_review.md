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

## G1 reading and the path to G2 (2026-09-22)

PI, verbatim (the message was pasted from a rendered page, so each display
formula appears twice — once as rendered text, once as its LaTeX source;
nothing has been edited):

> This is a real G1 pass, not a marginal one. The preregistered decision is satisfied: B1 is Green, B2 is Green on Ce and Nd, and B3 is Yellow rather than Red, so the project should proceed to G2 exactly as planned. The physically important result is that the energy-conserving redo preserves the compression phenomenon: La II and Nd II meet the 0.1-mag criterion already at Ng=2N_g=2, while Ce II reaches it at Ng=16N_g=16. At the same time, the best scalar is essentially useless for Ce and Nd—ϵ⋆=0\epsilon^\star=0 still misses by about 1–1.4 mag.
> The strongest new scientific statement from F67 is therefore narrower and better than the old Paper III claim:
> Under physical P1 conditions and energy-conserving fluorescence,a low-dimensional frequency-group closure reproduces observablesthat no scalar ϵ can reproduce.\boxed{ \text{Under physical P1 conditions and energy-conserving fluorescence,} \quad \text{a low-dimensional frequency-group closure reproduces observables} \quad \text{that no scalar }\epsilon\text{ can reproduce.} }
> That is enough to justify G2.
> There are, however, two things in F67 that materially change how I would frame the next stage.
> First, the observable compression is much stronger than the event-level compression. Even at 32 groups, the event-level total-variation loss remains around 0.16−0.420.16{-}0.42, while the photometric errors are already at the few-hundredths-of-a-magnitude level. The report itself notes that most of the information discarded by the binning never reaches the light.
> That is potentially a deeper result than "the matrix is small." It suggests that the transport problem has an observable-relevant effective dimension much smaller than the microscopic event distribution. I would explicitly carry that into G2. The question becomes not just
> local coarse graining vs low rank\text{local coarse graining vs low rank}
> but
> Which microscopic degrees of freedom actually survive transport into observables?\boxed{ \text{Which microscopic degrees of freedom actually survive transport into observables?} }
> That is a very PRL-like question.
> Second, the current implementation exposes a serious issue with the word compact. The matrix itself is tiny, but the discrete exit tables dominate storage: Ce has about 12,50012{,}500 stored exit samples and Nd about 86,00086{,}000, with Nd's serialized kernel around 2 MB. The report correctly says that the whole representation is not yet compact; only the group-to-group transition matrix is.
> I would not treat that as a failure of the science yet, because it may be an implementation representation issue rather than a fundamental information requirement. If the exit tables literally store repeated Monte Carlo event frequencies, then one obvious next check is whether identical exit frequencies can be aggregated exactly into unique line frequencies plus accumulated weights. That would preserve the empirical discrete distribution while replacing an O(Nevents)O(N_{\rm events}) table with something closer to O(Nunique exit lines)O(N_{\rm unique\,exit\,lines}). If that is true in the implementation, it should be done before we make any final statement about the memory complexity of the operator. This is an inference from the reported table sizes, not something the report itself establishes.
> I would therefore proceed in this order:
>
> 1. Do the preregistered R2M robustness check now. Before investing in G2, rerun the three ions with the radiation-field-driven macroatom/upward transitions and rebuild the Ng=8N_g=8 operator from those events, as already specified. G1 does not depend on this, but it tells us whether the compression result is tied specifically to the downward macroatom.
> 2. Audit the exit-table representation. Determine how many stored samples are duplicates of the same physical exit line and whether exact aggregation by νout\nu_{\rm out} and accumulated energy weight is possible. Do not change the physics or the G1 result; this is a representation audit.
> 3. Then freeze G2. Use the energy-conserving 128-group operator and compare local frequency coarsening with a non-negative constrained global factorization at matched complexity. Score both on transport observables, not only row distance. Given G1, the important comparison is whether a representation with worse event-level fidelity can nevertheless preserve the observable-relevant structure.
> 4. Keep La in G2, but don't let it drive the interpretation. At its physical P1 abundance it has only about 0.7 interactions per packet, so every model is close. Ce and Nd are the discriminating ions. The report explicitly recognizes La as too thin to distinguish the closures strongly.
>
> The central G2 hypothesis I would now preregister as:
> The success of Rij is caused by local smoothness in frequency,̸ by a globally low-rank microscopic redistribution law.\boxed{ \text{The success of }R_{ij}\text{ is caused by local smoothness in frequency,} \not \text{ by a globally low-rank microscopic redistribution law.} }
> But I would add an even more interesting secondary hypothesis suggested directly by G1:
> Transport observables occupy a much smaller effective information space than the underlying fluorescence-event distribution.\boxed{ \text{Transport observables occupy a much smaller effective information space} \text{ than the underlying fluorescence-event distribution.} }
> That second statement is supported by the striking mismatch between the still-large event TV distance and the already-converged magnitudes. If G2 shows that local coarsening preserves observables while equal-parameter global factorization does not, despite both losing large amounts of microscopic event information, then we have something conceptually stronger than a numerical surrogate.
> So my reading is: G1 did exactly what we needed it to do. Continue. But the path to PRL is now less about proving that a matrix can approximate fluorescence—we have that—and more about explaining why so much microscopic information can be discarded without changing the observable radiation field.

Decisions taken from it, in the PI's order:

1. **R2M robustness now**, exactly as the post-pass paragraph of
   `prl_gate.md` specifies (the radiation-field-driven macroatom, W = ½,
   the 8-group operator rebuilt from its events); its run parameters are
   fixed in a new section of `prl_gate.md` before the run. Not a gate.
2. **Exit-table audit**, representation only. Checked before any code was
   written: `RedistributionKernel.from_branching_mc` already aggregates the
   exits by unique frequency (`np.unique(nu_out)`, weights summed), and
   every stored exit is an exact line rest frequency; the stored count is
   the number of *distinct exit lines the sample discovered*, which grows
   with the sample (Nd II 59,270 at 3×10⁵ packets, 86,187 at 10⁶, of 3.34 M
   lines). The audit therefore reports the discovery curve, the energy
   concentration across exit lines and the size decomposition; any further
   shrinking changes the sampled distribution and is an experiment (G2's
   second hypothesis), not a representation fix.
3. **G2 preregistered** with the two hypotheses (H1 locality vs global
   low rank at matched archetype count on shared 128-group exit tables;
   H2 the observable-relevant dimension, by exit-table truncation), scored
   on transport observables against G1's reference; the PI approves the
   criteria before the run.
4. **La II stays in G2 and is reported, never counted** in H1/H2; Ce II
   and Nd II decide.

## The PI's review of PR #10: G2 amended, merge, tag, run (2026-09-22)

PI, verbatim (pasted from a rendered page: each display formula appears
once as rendered text and once as its LaTeX source; nothing edited):

> I would merge PR #10 after one small preregistration amendment, then run G2. The F68 result is exactly the robustness check we needed: the compression follows the redistribution physics of the reference, not a peculiarity of the downward macroatom. The Ce II case is especially decisive: changing the reference by 2.09 mag but recovering the new reference to 0.026 mag with a freshly trained R_8, while the old kernel stays ~1.97 mag away, is a very clean demonstration.
> For G2, I approve the basic design—shared 128-group exit tables, only the transition matrix varying, and comparing local frequency coarse-graining to a global archetype/NMF representation. I would keep the archetype-count comparison rather than reverting to nominal parameter count, because the global model actually gets many more free numbers for the same number of archetypes. If it still needs as many archetypes as the local representation, that is strong evidence for locality.
> I would make one amendment to H1: do not require the global model to "fit the events better" as part of the Green condition. Make event fit a mechanistic diagnostic, not part of the pass/fail logic. The central comparison should be the minimum number of archetypes needed to meet the same transport threshold:
>
> K_{\rm local}^\star
> =
> \min\{N:\text{G1 transport criterion passes}\},
>
> K_{\rm global}^\star
> =
> \min\{k:\text{same transport criterion passes}\}.
>
> Then I would preregister:
>
> * Green: K_{\rm global}^\star \ge K_{\rm local}^\star for both Ce II and Nd II.
> * Red: the global model passes on either Ce or Nd with
>     K_{\rm global}^\star \le K_{\rm local}^\star/4.
> * Yellow: everything between those cases.
>
> Alongside that, report event-level TV/L1. The particularly interesting outcome would be
>
> \boxed{
> \text{global model fits microscopic events better}
> \quad\text{but}\quad
> \text{local model reproduces transport better}.
> }
>
> That would directly support the emerging G1 story: microscopic-distribution fidelity is not the same as observable-relevant fidelity. I would not make that circumstance necessary for Green, though, because it is a stronger and more specific result than the actual locality hypothesis.
> For H2, I like the new exit-table truncation experiment, and the extension to f=0.1,0.2 is legitimate because it followed the audit and preceded any G2 transport. I would keep the headline threshold of 10% of distinct exit lines, but define the reading across the two informative ions:
>
> L^\star_a =
> \text{minimum retained fraction of unique exit lines satisfying G1 transport thresholds}.
>
> Then:
>
> * Green: L^\star\le 0.10 for both Ce II and Nd II.
> * Yellow: at least one is \le0.25, but Green fails.
> * Red: either requires >0.50 of its exit lines.
>
> La can be reported as a control but should not decide H2, because G1 already showed it is too thin to discriminate closures strongly.
> There is one implementation detail I would freeze explicitly before running H2: line selection must be based only on the kernel-building sample, never on transport performance. For example, rank unique exit lines by their accumulated training-event energy contribution within the relevant exit table, retain according to the preregistered f, and renormalize. Do not choose different subsets after seeing which preserve the light best.
> The scientific logic after these changes is very clean:
>
> \boxed{
> \text{H1: Is the transition structure locally smooth rather than globally low-rank?}
> }
>
> \boxed{
> \text{H2: Is the discrete within-group exit spectrum itself highly compressible?}
> }
>
> If both go Green, then the earlier concern about "the matrix is small but the exit table is not" largely disappears. You would then have evidence that both layers of the effective operator are compact.
> F68 also changes my view of the potential paper slightly. The reference-physics shift is large—especially Ce's 2.09 mag—but a small freshly constructed operator follows it. That is stronger than merely saying "a small matrix reproduces one macroatom." It suggests
>
> \boxed{
> R_{ij}\text{ is a compact representation of the chosen microscopic fluorescence physics,}
> }
>
> not a fortuitous calibration to one particular cascade prescription.
> I would therefore make the sequence:
>
> 1. Amend H1 as above and make the H2 selection rule explicit.
> 2. Merge PR #10.
> 3. Tag/freeze the G2 preregistration state.
> 4. Run G2.
> 5. Do not touch G3 until the G2 verdict is written.
>
> On your other open items, I would not let the education PRs #5–#8 block this research run; they are valuable but orthogonal. I also would keep the Paper IV MNRAS submission moving independently—the fact that Paper B is getting stronger is not a reason to hold a completed, frozen Paper IV. The report already treats Paper IV as complete on the corrected transport and tagged paper4-freeze. results_report.md
>
> So my decision is: approve G2 with the H1 amendment above; merge PR #10; run G2.

Applied before any G2 transport existed (`prl_gate.md`, "Amended
2026-09-22"): H1 reads on the archetype counts K*_local and K*_global
alone, with the event-level fit demoted to a named diagnostic
(`events_vs_observables`); H2 reads on the retained fraction L* with the
three-level ladder and Ce II / Nd II decisive, La II a control; the H2
selection rule written out as a function of the build sample alone. One
residual band the PI's three cases leave open — both decisive ions
between 0.25 and 0.50 — is preregistered as Yellow and reported with its
values, the ladder read Red first, then Green, then Yellow.

Sequence as instructed: amend, merge PR #10, tag the G2 preregistration
state, run G2, leave G3 untouched until the G2 verdict is written.

## G2 run, first attempt: two implementation defects, disclosed and repeated (2026-09-22)

Not a PI statement: a record of a decision I took while executing the
approved sequence, because it required touching a tagged preregistration
after part of a run existed.

The G2 run was launched from `paperB-g2-prereg`. La II finished; reading
its record showed `G_k1` and `G_k2` failing the energy identity
(`validate_energy` = 1.0, gray condition 3) because the factorisation sent
one live row to exactly zero, and `G_k16`/`G_k32` failing the convergence
condition (gray 8) because 600 Lee–Seung iterations left the relative
Frobenius error falling by ~3×10⁻³ per 100 iterations. Both are failures
to compute the preregistered object, not properties of it, so the
preregistration's own remedy applies: fix and repeat the run.

I stopped the run with Ce II in progress, deleted every G2 record unread
except La II's (which had already been inspected and which decides
neither reading), fixed both defects, recorded the amendment in
`prl_gate.md` — including exactly what had been seen, and the direction
in which the zeroed-row convention leans — and reran all three ions.

The judgement call worth recording: the zeroed row could be repaired
either by declaring it empty (coherent fallback, the kernel's standing
convention for a row with no information) or by falling back to its own
un-factorised distribution. The first leans toward H1 Green, the second
toward Red and also leaks information the rank-k model does not contain.
I chose the first and wrote the bias into the preregistration rather than
leaving it implicit; on La II the row carried 1.5×10⁻⁶ of the absorbed
energy.
