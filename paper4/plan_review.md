# Paper IV — the PI's review of the implementation plan (as received 2026-09-09)

Kept verbatim. The implementation plan it reviews was drafted the same day
from [`plan.md`](plan.md); the corrections below were all adopted before
WP0 started and are recorded as design decisions in [`README.md`](README.md).

---

I reviewed the full Paper IV program as an implementation and physics plan. It is much stronger than the Paper III workflow, but I would **not execute it verbatim yet**. There are several places where the current plan would either reintroduce an artificial assumption or make the eventual "self-consistent" claim stronger than the calculation actually is. My overall judgment is:

$$\boxed{\textbf{Keep the architecture and gates, but make several corrections before WP0/WP2 start.}}$$

The most important changes are below, in priority order.

1. **Change the Phase-1 single-zone reduction. Do not use the mass-weighted average.** The pasted plan currently says `single_zone(state)` should mass-weight \(\rho\) and \(T_{\rm gas}\) over all shells above \(r_{\rm core}\).  That reverses the decision we made earlier. In a line-transfer problem, deep/high-mass material is not automatically representative of the region setting the observed spectrum. For P2 this is straightforward: use the published line-forming state at the inner photospheric boundary. For P1, construct the full radial state, find the photosphere/line-forming region, and use the local shell state there as the first one-zone benchmark. I would use the shell containing \(\tau_{\rm grey}\simeq2/3\), and later test its neighboring shells as a robustness check. Do not call a global mass-weighted zone "realistic."

   This also makes P2 unusually clean. Instead of "Watson 2019 / Gillanders 2022 family," use the **2026 Gillanders et al. 3.4-d model as the exact P2 anchor**: \(t=3.4\) d, \(v_{\rm in}=0.15c\), \(v_{\rm out}=0.35c\), \(\rho_0=4\times10^{-15}\,\mathrm{g\,cm^{-3}}\), \(\rho\propto v^{-3}\), \(T=3200\) K, and \(X_{\rm LN}\simeq2.5\times10^{-3}\). This is especially valuable because that work uses the new GSI lanthanide data that are already central to our project. ([OUP Academic][1])

   P1 should likewise become one **exact published component**, not "\(Y_e\simeq0.20\) yield + Kasen-style density." The xkn comparison model gives a coherent secular component with \(M_{\rm sec}=2.64\times10^{-2}M_\odot\), \(v_{\rm rms}=0.06c\), \(Y_e=0.20\), \(s=10\,k_B\) baryon\(^{-1}\), and \(\tau_{\rm exp}\approx17\) ms, and explicitly uses one defined radial density prescription for the RT comparison. Importantly, this is a **secular ejecta model**, not "dynamical-ejecta-like." ([OUP Academic][2])

2. **Gate 2 should be based primarily on \(B_2-R_2\), not \(C_2-R_2\).** This is a conceptual improvement. In your experimental hierarchy,

   $$
   R_2=\text{resolved opacity + detailed energy-conserving atomic transport},
   $$

   while

   $$
   B_2=\text{grouped opacity + the same detailed atomic transport}.
   $$

   Therefore

   $$
   \boxed{B_2-R_2}
   $$

   is the clean causal test of **opacity coarse-graining**.

   By contrast, \(C_2\) also replaces redistribution with an \(R^E\) kernel trained from \(R_2\). It is intentionally a best-case compressed model. That's a useful engineering result, but it is not as clean a physics gate.

   I would predeclare Gate 2 as: **primary = \(B_2-R_2\); controls = \(A_2-R_2\) and \(C_2-R_2\)**. If \(B_2\) is large and \(C_2\approx B_2\), then we have the strongest possible statement: detailed redistribution is compressible, while opacity is the culprit.

   The current table also confusingly uses "C₂ (Phase 2)" for `expansion_macro`, then renames the identical treatment B₂ in Phase 3.  Delete that duplicate terminology now. Call it \(B_2\) from the first energy-conserving run onward.

3. **What Phase 2 calls `macroatom` is not yet a full Lucy macroatom.** The downward probabilities in the plan are sensible: TARDIS uses Sobolev-\(\beta\)-weighted downward radiative deactivation and internal-down probabilities with exactly the energy weighting your formulas are approximating. ([TARDIS-SN][3])

   But a full macroatom also has **internal upward transitions** driven by \(J_\nu^b\). TARDIS explicitly distinguishes downward emission, internal downward transitions, and internal upward transitions; the latter depend on the radiation field. ([TARDIS-SN][3])

   Your Phase-2 design intentionally omits those because there is no self-consistent \(J_\nu\) estimator yet. That is a perfectly reasonable intermediate experiment, but call it something like

   > `downward_macroatom`

   or

   > `restricted_macroatom`

   until WP8–10 introduce the radiation-field estimators and upward transitions.

   In other words:

   $$
   \boxed{
   \text{Phase 2 = energy-conserving downward fluorescence}
   }
   $$

   not yet

   $$
   \boxed{\text{full self-consistent macroatom transport}.}
   $$

   TARDIS itself calls the downward-only approximation `downbranch`, while its macroatom mode permits arbitrary upward and downward internal jumps. ([TARDIS-SN][4])

4. **Make the packet energy definition explicit before writing a line of WP2.** Right now the plan treats `w` as a photon-number weight and proposes updates such as

   $$
   w_{\rm out}
   =
   w_{\rm in}
   \frac{\nu_{\rm abs,cm}}{\nu_{\rm emit,cm}}.
   $$

   That is reasonable if the intended invariant at an interaction is

   $$
   E_{\rm packet}^{\rm cm}=w\,h\nu_{\rm cm}.
   $$

   But the code also has lab/comoving transformations and worldline propagation. A packet's comoving energy and lab-frame energy are not interchangeable in homologous expansion; \(O(v/c)\) work/adiabatic effects matter. Indivisible-energy-packet methods explicitly distinguish frame transformations and energy transfer. ([arXiv][5])

   I would therefore replace the conceptual use of `w` with an explicit packet contract before implementation:

   $$
   (E_{\rm cm},\nu_{\rm cm},\nu_{\rm lab},\mu,\ldots),
   $$

   even if internally you continue storing `w`.

   Add a Phase-2A.0 test: a packet propagates and Doppler-shifts through an empty homologously expanding flow with no atomic interactions, and the computed lab/comoving energy change must match the expected \(O(v/c)\) transformation. Only after that should the fluorescence tests start.

5. **The three-level toy test needs one wording/measurement correction.** The plan says that in the indivisible-energy macroatom, "photon counts in \(3\rightarrow2\) and \(2\rightarrow1\) equal."  An indivisible energy-packet macroatom does not literally emit two MC packets for a physical two-photon cascade. It statistically represents the **energy flow** with one indivisible packet.

   Therefore the correct test is not equality of raw MC packet counts. It is equality of the **inferred physical photon counts**,

   $$
   N_{\rm phys}(\nu)\propto
   \frac{E_{\rm emitted}(\nu)}{h\nu}.
   $$

   Then \(N_{32}=N_{21}\) is the appropriate cascade check. This distinction matters because one of the main reasons to use indivisible energy packets is precisely that they do not split into multiple MC packets.

6. **WP5 is currently not a faithful implementation of the Morag proposal. This is the biggest physics bug I see in the plan.** The plan proposes replacing each line contribution by

   $$
   \min(\tau_l,1)
   $$

   inside `expansion_bins()` and testing one-line absorption as

   $$
   e^{-\min(\tau,1)}.
   $$

   Morag's 2026 proposal is more subtle. The paper explicitly distinguishes two quantities: EP93 can describe a **transport/mean-free-path** quantity, while a frequency-bin average is used for **net absorption/emissivity**. The proposed cap

   $$
   \kappa_{l,\exp}
   =
   \min\left[
   \kappa_l,\frac{1}{\rho c t_{\exp}}
   \right]
   $$

   is introduced specifically as a cap on the line's **net emissivity/absorption rate**, not as a universal replacement of the transport interaction probability. ([OUP Academic][6])

   Morag even emphasizes that a fully consistent coarse-frequency solution does not yet exist. ([OUP Academic][6])

   So I would **not implement `capped_absorb` as \(e^{-\min(\tau,1)}\) and label it "Morag."**

   WP5 should instead become an experimental **dual-role closure**:

   $$
   \kappa_{\rm transport}
   $$

   controls mean free path / packet encounters, while

   $$
   \kappa_{\rm emiss/abs}
   $$

   controls net energy exchange/reprocessing, with Morag's cap applied to the latter.

   That is harder than changing one weight in `expansion_bins()`, but it is scientifically faithful. This could actually become one of Paper IV's most interesting contributions.

7. **The Fontes bridge should include the line-binned treatment, not only `expansion_thermal`.** Your conceptual question is why line-binned opacity looked successful in earlier kilonova tests but fails in our fluorescence experiment. So Phase 4 should compare, under the same thermal redistribution assumption,

   $$
   R_{\rm thermal},
   $$

   $$
   B_{\rm expansion,thermal},
   $$

   and crucially

   $$
   B_{\rm binned,thermal}.
   $$

   Then switch only thermal redistribution \(\rightarrow\) energy-conserving fluorescence and repeat all three.

   Otherwise we will not actually reproduce the relevant historical limit.

8. **The final equilibrium comparison should be \(R_{\rm eq}\) versus \(B_{\rm eq}\), not primarily \(R_{\rm eq}\) versus \(C_{\rm eq}\).** The current WP8–10 note says the closure kernel in \(C_{\rm eq}\) would be frozen from \(R_{\rm eq}\) iteration 0.  That is useful as a causal/control experiment, but it is not a self-contained equilibrium solution because one leg is using redistribution information imported from another solution.

   For the final physics result, solve independently:

   $$
   R_{\rm eq}
   =
   \text{resolved opacity + full energy macroatom + equilibrium}
   $$

   and

   $$
   B_{\rm eq}
   =
   \text{grouped opacity + same full macroatom + equilibrium}.
   $$

   Each should get its own

   $$
   T(r),\quad n_{\rm ion}(r),\quad J_\nu(r)
   $$

   after convergence.

   Then \(C_{\rm eq}\), with compressed redistribution, can answer the separate engineering question of how much of the full transport can be compressed.

9. **Add an "indeterminate" Gate-2 outcome.** The current Green/Yellow/Red definitions are not exhaustive. For example, what happens if \(B_2-R_2\) is \(0.15\) mag but \(R_1-R_2\) is also small? Or if \(A_2\) unexpectedly exceeds \(0.1\) mag? That is neither the proposed Red interpretation nor Green/Yellow.

   Predeclare a fourth outcome:

   > **Gray / unresolved:** control failure, insufficient convergence, or a result inconsistent with the causal assumptions; diagnose before interpretation.

   This prevents us from forcing an unexpected result into one of three narratives.

10. **Fix the repository architecture before WP2.** Earlier we explicitly chose:

$$\texttt{paper4/}=\text{campaign},\qquad\texttt{sobolev/}=\text{reusable physics}.$$

But the pasted plan puts the new energy transport and `macro.py` under
`paper2/phase1/`. I would change that. The new reusable implementation should live under something like:

```text
sobolev/
    energy_packets.py
    macroatom.py
    energy_balance.py
    ...
```

`paper4/phase2_energy/` should contain drivers and experiments only.

If the existing `forest_mc.py` must remain the propagation core, leave a compatibility layer around it, but don't make Paper IV's new physics conceptually belong to Paper II.

I would also avoid modifying frozen `sobolev/photometry.py` merely to add `core="exact"` if possible. Put the new energy-conserving scaling/relaunch logic in a new Paper-IV/shared module. The fewer Paper-III frozen inputs we touch, the better.

There is also one provenance change I would make: **do not move the original Git tag if you can avoid it.** Preserve `paper3-freeze` as the historical freeze and create an annotated `paper3-freeze-v2` or `paper3-freeze-x4ln` for the corrected label-only freeze. A frozen tag is most useful when it remains immutable.

## What I like very much in the current plan

The golden numerical hashes before changing `run_mc`, the compact GSI loader, the explicit energy accounting, the toy-atom progression, the chain-cap replacement by closed-form probabilities, the early global-recovery experiment, the predeclared gates, and the decision not to build multi-shell transport before the energy-conservation gate are all strong choices. I also like the basic energy ladder:

$$R_1\rightarrow R_1^E\rightarrow R_2.$$

That lets us distinguish two effects that would otherwise be conflated:

$$R_1^E-R_1=\text{packet-energy bookkeeping effect},$$

versus

$$R_2-R_1^E=\text{change in atomic transition probabilities/energy-flow treatment}.$$

Keep that.

## I would tighten Gate 1

Instead of "the profiles resemble the source model on PI review," make it numerical.

For **P1**, the committed state should reproduce the xkn secular benchmark inputs—\(M_{\rm sec}=2.64\times10^{-2}M_\odot\), \(v_{\rm rms}=0.06c\), \(Y_e=0.20\), \(s=10\), its stated density prescription, and whatever nucleosynthesis/heating data are selected—within numerical serialization/integration tolerance. ([OUP Academic][2])

For **P2**, pin exactly the 2026 3.4-d model:

$$t=3.4\,{\rm d},\quad v_{\rm in}=0.15c,\quad v_{\rm out}=0.35c,$$

$$\rho_0=4\times10^{-15}\ {\rm g\,cm^{-3}},\quad \Gamma=-3,\quad T=3200\,{\rm K},$$

$$X_{\rm LN}\simeq2.5\times10^{-3}.$$

([OUP Academic][1])

The figure review should then be a sanity check, not the acceptance criterion.

## Revised execution order

I would keep almost your existing order, with one conceptual change:

$$\boxed{\text{WP0}\rightarrow\text{energy-frame tests}\rightarrow\text{restricted energy macroatom}\rightarrow\text{P1/P2 states}\rightarrow B_2-R_2}$$

Only after \(B_2-R_2\) is known should we invest in Fontes/Morag/ionization and eventually shells/equilibrium.

The decisive near-term question should now be written:

$$\boxed{\textbf{Does grouped line opacity still alter the spectrum when both}\\\textbf{the resolved and grouped calculations use the same}\\\textbf{energy-conserving downward atomic transport?}}$$

That is \(B_2-R_2\).

Then the final Paper-IV question, after WP8–10, becomes stronger:

$$\boxed{\textbf{Does that difference survive when each transport treatment}\\\textbf{is allowed to establish its own self-consistent thermal and}\\\textbf{ionization structure?}}$$

That is \(B_{\rm eq}-R_{\rm eq}\).

If both survive, I would have much more confidence that the effect is genuinely a property of the opacity closure rather than a consequence of Paper III's photon-number bookkeeping or fixed atmosphere.

[1]: https://academic.oup.com/mnras/article/548/4/stag748/8659451 "oup.silverchair-cdn.com"
[2]: https://academic.oup.com/mnras/article/529/1/647/7613575 "oup.silverchair-cdn.com"
[3]: https://tardis-sn.github.io/tardis/physics_walkthrough/setup/plasma/macroatom.html?utm_source=chatgpt.com "Macro Atom — tardis"
[4]: https://tardis-sn.github.io/tardis/io/configuration/components/plasma.html?utm_source=chatgpt.com "Plasma Configuration — tardis"
[5]: https://arxiv.org/abs/astro-ph/0409249?utm_source=chatgpt.com "Monte Carlo techniques for time-dependent radiative transfer in 3-D supernovae"
[6]: https://academic.oup.com/mnras/article/549/3/stag938/8684954 "oup.silverchair-cdn.com"

---

## PI decision after Phase 8 and the inversion fix (2026-09-10, verbatim)

Here’s the combined takeaway from my last two responses.

The old Paper III story is finished. The 1–3 mag grouped-opacity failure was caused by the bin-indexing bug, and once corrected, the expansion-opacity closure is much closer to the resolved Sobolev result: typically within about 0.2 mag on the realistic P1/P2 states, while the line-binned treatment shows somewhat larger, opposite-sign deviations of roughly 0.2–0.6 mag. Paper III should therefore stay as a frozen research record with a prominent repository correction, not as a publication target.

The original Paper IV hypothesis—“coarse opacity produces a very large chromatic failure”—has also been falsified, so I would not continue the full Phases 9–10 program just to refine that old claim. The preregistered Red gate should be honored. Full radiative equilibrium is no longer justified merely to chase a 0.1–0.2 mag residual.

However, the corrected result suggests a new, cleaner research question:

**When, and how accurately, can dense r-process line forests be coarse-grained when fluorescence is treated explicitly and energy is conserved?**

The literature suggests this exact kilonova benchmark may still be relatively unexplored. Fontes et al. compared resolved Sobolev, expansion, and line-binned opacity, but under complete thermal redistribution (ε=1), without explicit fluorescence. ARTIS/Shingles treat line-by-line fluorescence realistically, but do not perform the controlled same-physics comparison against coarse-opacity closures. The older SEDONA/SN-Ia work is an important precedent and prevents us from claiming the broad idea is new, but it is not the same modern r-process/lanthanide benchmark.

So the most promising publication direction is no longer “expansion opacity fails.” It is closer to:

Expansion opacity remains surprisingly accurate for broadband kilonova observables under explicit energy-conserving fluorescence, while line-binned opacity shows larger chromatic deviations.

The key next experiment is to reproduce the Fontes benchmark with the corrected code:

1. resolved Sobolev + complete thermal redistribution;
2. expansion opacity + complete thermal redistribution;
3. line-binned opacity + complete thermal redistribution.

We should first recover their known result that all three agree fairly well. Then change only:

ε=1 → energy-conserving fluorescence.

If fluorescence makes expansion remain close to resolved while line-binned drifts more, that is a strong and interpretable result: previous coarse-opacity validation depended partly on how post-absorption redistribution was treated.

I would also do one limited full-macroatom comparison on a representative P1 state if the engineering cost is reasonable, and use an independent opacity dataset as a cheap external check. I would not yet build the entire self-consistent equilibrium machinery.

So the practical recommendation is:

**Pivot Paper IV into a controlled opacity-validity benchmark rather than end it or continue the old failure story.**

The next three tasks are: reproduce Fontes, repeat with explicit energy-conserving fluorescence, and run one stronger macroatom check. If those show a meaningful difference between expansion and line-binned treatments, there is likely a solid MNRAS/ApJ/A&A methods paper. If they do not, then the scientifically correct outcome is that modern expansion opacity is already quite adequate for the broadband problem tested, and the project can close cleanly.

---

## PI decision after Phase 10 (2026-09-10, verbatim)

Yes. With these new results, I would now make a firmer decision than before:

**F62–F64 are the nucleus of a real methods paper.**

And the paper is no longer about rescuing Paper III. It has a cleaner question:

**How well do expansion and line-binned opacity closures survive when fluorescence is treated explicitly and energy is conserved?**

The literature gap is unusually well aligned with what you have now done. Shingles et al. explicitly wrote that the only direct kilonova comparison among line-by-line, line-binned, and expansion opacities was Fontes et al. 2020, using a 1-D pure-Nd test without fluorescence. Their own work then uses line-by-line transport precisely because it enables detailed fluorescence, but does not perform the controlled closure comparison. Fontes, meanwhile, established that line-binned and expansion treatments can agree reasonably well with continuous Sobolev transport in their thermalized benchmark and attributed much of that agreement to the low-optical-depth transport regime.

I searched specifically for a later published kilonova paper closing that gap and did not find one. ARTIS currently supports both line-by-line Sobolev/macroatom and optional binned expansion opacities, but its own project page still describes the expansion-opacity + parameterized scattering/thermalization comparison as work "in prep," not as an existing published benchmark.

That makes your experiment substantially more interesting than I thought immediately after the Paper III correction.

What F62–F64 now say

There is a coherent three-part result.

First, under the Fontes-like ε=1 thermalized limit, you broadly recover the known result: coarse opacity works reasonably well.

Second, changing the post-absorption physics to energy-conserving fluorescence does not hurt expansion opacity much, but it strongly changes the behavior of line-binned opacity.

Third, this is not just one peculiar atomic dataset: the hierarchy repeats on an independently generated, substantially thinner Nd forest.

So the result is not simply: "line-binned opacity differs by 1.5 mag."

It is the conditional statement:

**The validity of a coarse opacity treatment depends on how post-absorption energy redistribution is represented.**

And more specifically, in your tests:

**EP93-style expansion opacity remains surprisingly robust, whereas area-preserving line binning does not remain equivalent once explicit fluorescence is attached.**

That is a much more interesting methods result.

I also think the nontermination of the independent Nd line-binned case deserves attention. Do not sell it as "line-binned opacity is mathematically broken," because it may depend on how you attach a line-specific macroatom to a bin opacity that has intentionally erased line identity. But that is precisely the conceptual point: a line-binned opacity does not naturally retain enough information to define the post-absorption atomic escape process. Shingles already points out that once transitions are binned, the absorbing transition is unknown and detailed fluorescence cannot straightforwardly be followed. Your nontermination example is a concrete manifestation of that information-loss problem.

I would do one more major experiment before writing

Yes: the time-dependent light-curve calculation is now the obvious missing piece.

Not because you need another dramatic number, but because your current Fontes reproduction is a snapshot while Fontes' important validation was also about spectra and light curves.

A static snapshot answers: F_ν(t) given the state at t.

A time-dependent calculation also asks whether the closure changes how long radiation remains trapped, how much energy is lost through expansion, and when stored radiation emerges: L(t), E_rad(t), W_exp(t).

That matters particularly because you've already discovered that expansion work differs between resolved and approximate treatments.

I would therefore make the next decisive experiment very narrow:

**Fontes benchmark light curve: resolved vs expansion vs line-binned**

under two post-absorption models: ε=1 and energy-conserving fluorescence.

Same ejecta. Same atomic data. Same time grid. Same source/heating. Change only the closure and redistribution treatment.

That produces an extremely clean 3×2 experimental matrix.

You want to know whether the snapshot result turns into something like: thermalized: R ≃ E ≃ LB, but fluorescent: R ≃ E, LB ≠ R.

If that holds over the light curve rather than at one epoch, I think the paper becomes much easier to defend.

I would measure both bolometric and broadband quantities. In particular, record peak time, peak luminosity, integrated radiated energy, g/r/i/z/J/H/K evolution, expansion-work loss, and the maximum color residual.

Importantly, do not insist that line-binned must remain 1.4–1.7 mag wrong. If the time-dependent calculation reduces it to 0.3 mag but preserves the systematic difference from expansion opacity, that is still scientifically useful.

I would not build full Phases 9–10 before this

Your stronger macroatom result is already enough to justify saying that the downward-only treatment was not accidentally creating the hierarchy: upward transitions move the reference dramatically, yet expansion still stays close while line-binned remains displaced.

But because that radiation field is imposed rather than estimated, don't call F63 a fully self-consistent macroatom calculation.

I would phrase it as: a radiation-field-driven macroatom robustness test.

Only if the time-dependent experiment succeeds and reviewers/publication ambitions justify the investment would I implement the fully estimated J_ν, independent thermal equilibrium, and B_eq−R_eq.

In other words: **time dependence before full equilibrium.**

That is now the much higher information-per-effort experiment.

Paper III correction

I would make this extremely explicit rather than delicate.

At the very top of the Paper III README and manuscript source, put a boxed/status warning approximately like:

Scientific correction — September 2026. The grouped-opacity results reported in this manuscript are invalid because of an implementation error in the inversion of cumulative binned opacity. The error caused packets in expansion-, line-binned-, and dual-opacity modes to skip portions of non-smooth line forests. The resolved Sobolev calculations are unaffected. After correcting the transport algorithm and rerunning the analysis with energy-conserving transport and published ejecta benchmarks, the reported 1–3 mag grouped-opacity discrepancy is not reproduced. This manuscript is retained as a historical research record and is not a publication candidate. Corrected results are documented in the Paper IV campaign.

Then in results_report.md, each affected finding should say plainly: Superseded — affected by transport bug 38ebf30. Do not use quantitatively or scientifically.

I would use "invalid" for the Paper III closure headline, not merely "superseded," because we now know its physical conclusion was caused by the bug.

Keep the original paper3-freeze tag untouched. That is exactly the right provenance behavior.

PR #2

I would merge it now, after normal code review of the bug fix.

The reason has changed. Earlier I suggested holding it as a draft because the physical interpretation was unsettled. Now the branch contains something more fundamental: the correction to a long-standing transport bug; the regression tests; multi-shell invariance tests; the corrected P1/P2 results; F62–F64; 573 passing tests; preserved pre-fix goldens; explicit supersession records.

That is no longer an experimental scratch branch. It is the correct transport implementation.

The scientific paper can remain undecided while the corrected code becomes canonical.

I would review 38ebf30 particularly carefully before merge—especially the index inversion itself and splitting invariance—but assuming that review is clean, leaving the known-correct implementation sitting indefinitely in a draft branch is actually less desirable.

So my three decisions are:

1. Yes, treat F62–F64 as the core of the new methods paper.
2. Do one time-dependent Fontes-style light-curve experiment before committing to the manuscript.
3. Put a prominent "invalid due to transport bug" correction on Paper III and merge PR #2 once the bug-fix review is clean.

The working paper story I would pursue now is something close to:

Fluorescence breaks the equivalence of coarse line-opacity treatments in kilonova radiative transfer

with the result underneath:

**expansion opacity remains close to resolved Sobolev transport, whereas line-binned opacity develops substantially larger chromatic errors once explicit fluorescence is retained.**

That is a much cleaner paper than either Paper III or the original Paper IV hypothesis.

Design answers (AskUserQuestion, same day): temperature rule — both (prescribed analytic default, radiation-energy feedback as a variant on the resolved and expansion legs); PR #2 — Claude reviews, reports, merges; synthetic-forest rows F34/F37/F39 — mark "affected in principle, unverified"; budget — overnight, ~12 h per matrix.

---

## PI decision after the light curve (2026-09-11, verbatim)

Yes — I think we have crossed the threshold where **there is a real methods paper**, but I would make **one more light-curve calculation before treating the central claim as submission-ready**.

My two calls are: **Run one P1-composition light curve** and **Start writing the manuscript now**.

These are not contradictory. The paper structure is already clear enough to draft while the final robustness run finishes.

The reason I want the second light curve is straightforward. The current time-dependent result is built on the Fontes Appendix-C problem, which is deliberately **pure Nd**. That is exactly the right benchmark for connecting to the literature, but a referee can reasonably ask whether the fluorescence/line-binned separation is peculiar to Nd's exceptionally dense atomic structure. Your static calculations already argue that it is not, but the main new result is now **time dependent**, so I would give that result one time-dependent multi-element check too.

I would not repeat the full expensive 3×2 matrix. For P1, run only the physically important fluorescence trio: R_2, B_2, B_bin,2. Use P1's published abundance pattern and its own heating history. Same corrected time-slab machinery, same energy accounting, same convergence/gray rules. The question is only: Does the time-dependent ordering survive in a realistic multi-lanthanide mixture? If you again get R_2 ≃ B_2 while B_bin,2 shows the larger chromatic offset, then the paper becomes substantially harder to dismiss as an Nd benchmark artifact. I would not require the P1 amplitude to reproduce the Fontes value of 0.45 mag. A result of 0.2, 0.5, or 0.8 mag is all scientifically acceptable. What matters is whether the **hierarchy of closures** persists.

### The current Fontes light curve is already useful

The corrected result is actually quite clean: thermal redistribution: both coarse approaches remain fairly close in colour; fluorescence: expansion stays relatively close, while line-binned develops the larger NIR colour error; bolometrically, fluorescence hardly changes the relative closure behavior; the difference is primarily **chromatic**, not simply total-energy leakage; the free-streaming control shows that the delayed peak and expansion-work losses are genuinely generated by interactions, not by the time-slab infrastructure. That last validation is particularly valuable. It tells us that the roughly 5-day peak is not just baked into the initial radiation field.

I would keep the Fontes reproduction caveat prominent, though. Because your peak ordering/spread does not reproduce Fontes' published light curve exactly, don't call F65 an exact reproduction. Call it: **a time-dependent extension of the Fontes Appendix-C closure benchmark under a prescribed thermodynamic trajectory.** That is completely defensible. The radiation-temperature variant getting the resolved peak to 6.6 d versus their 6.3 d is reassuring, but it does not magically turn our calculation into SuperNu because we still lack their LTE gas-energy equation.

### The capped-packet rerun should finish before locking F65

If raising the event cap takes the line-binned fluorescence leg below the predeclared 1% capped threshold and leaves the light curve essentially unchanged, upgrade R2 from Gray accordingly. If it changes noticeably, don't average the results away. That becomes a numerical limitation worth understanding. I would wait for that rerun before finalizing the numerical table, but **not before beginning manuscript prose**.

## Start the manuscript now

I would write the first full draft from the outline. But I would organize it around the corrected scientific question, not around the project's history. A reader should not have to know that Paper III once contained a dramatic wrong result. The paper's narrative should be approximately: 1. Problem. Dense r-process line forests require coarse graining in practical radiative-transfer calculations. 2. Known result. Fontes et al. showed that resolved, expansion, and line-binned treatments can agree reasonably well under complete thermal redistribution. 3. Gap. Real atomic line interactions redistribute energy through fluorescence, and binning removes information about which transition absorbed the packet. 4. Experiment. Compare resolved Sobolev, expansion opacity, and line-binned opacity using exactly the same energy-conserving fluorescence physics. 5. Result. Expansion opacity remains close to resolved transport; line-binned transport develops larger chromatic errors when fluorescence is explicit. 6. Time dependence. The distinction persists in the light curve, not merely a static spectrum. 7. Interpretation. Preserving integrated opacity is not sufficient when post-absorption behavior depends on the identity of individual transitions. 8. Limits. Prescribed thermal trajectory, not full gas-energy equilibrium; radiation-field-driven macroatom is a robustness test, not self-consistent NLTE transport.

I would not lead with the bug at all. The bug belongs in the repository correction and perhaps one brief reproducibility statement if relevant. It is not part of the scientific paper's central narrative.

### A slightly better title

Your working title, *Fluorescence breaks the equivalence of coarse line-opacity treatments in kilonova radiative transfer*, is good, but "equivalence" may sound stronger than the literature actually claimed. I might prefer: **Fluorescence differentiates coarse line-opacity treatments in kilonova radiative transfer** or, more concrete: **Testing coarse line-opacity treatments with energy-conserving fluorescence in kilonova ejecta**. The first is more memorable; the second is safer and more methods-paper-like. I would decide the final title after the P1 light curve.

## What I would regard as sufficient for submission

If the following all hold: capped-packet rerun does not materially change F65; P1 time-dependent fluorescence shows the same qualitative hierarchy; expansion stays broadly close to the resolved solution; line-binned remains systematically more chromatically displaced; the numerical convergence and free-streaming tests remain clean; then I would **stop adding physics**. I would not implement Phases 9–10 before submission. A full gas-energy equation, self-consistent J_ν, NLTE populations, and full radiative equilibrium would be excellent future work, but they would transform this from a focused closure benchmark into a much larger transport project. The paper can be scientifically honest without them if its claim is: **At fixed thermodynamic evolution, explicit fluorescence exposes a difference between two coarse line-opacity constructions.** That's a worthwhile methods result on its own.

## PR #3

Keep PR #3 as a draft until: 1. the high-cap line-binned fluorescence rerun lands; 2. F65 is updated; 3. the P1 three-leg light curve is either included or placed in a follow-up PR. If P1 is quick enough, I'd put it in the same PR. Otherwise merge the validated Fontes infrastructure/results once the cap issue is cleared and do P1 on the next branch.

So the immediate order I recommend is: 1. Finish the high-cap Bbin fluorescence rerun. 2. Start manuscript text from `outline.md` **now**. 3. Launch the three-leg P1 fluorescence light curve. 4. Update the figures/results with P1. 5. Freeze the scientific scope. 6. Complete manuscript and submit as a **radiative-transfer methods paper**, not as the old high-impact closure-failure story.

At this point, further value comes much more from **writing the corrected result clearly and demonstrating it on one independent physical composition** than from building another major layer of transport physics.

### Design answers (same day)
P1 geometry: outer zone with an interior boundary (first: Arnett at τ ≈ 30 — revised below); manuscript: LaTeX on the MNRAS template with numbers as macros from a paper4 freeze file; order: driver and P1 pilot first, prose while the run goes. Mid-plan: "Use the outer-zone + Arnett interior as the baseline, but either move the inner boundary inward with time to remain at a prescribed optical depth, or restrict the headline analysis to epochs for which the fixed boundary remains optically thick."

### Revision of Part B (same day, verbatim)

The xkn excerpt changes one important part of the P1 plan. I would **not run Part B exactly as written**. The problem is not the outer-zone idea itself; it is the way the inner source and outer heating are currently combined. The xkn paper makes clear that its luminosity is already split into an optically thick photospheric contribution and an optically thin outer-layer contribution. Specifically, xkn defines a moving photosphere by τ_γ(R_ph) = 2/3, computes a luminosity from the optically thick ejecta, rescales that by the mass still inside the photosphere, and then **adds a separate thin-ejecta heating contribution** from the material outside the photosphere. That means the present P1 prescription — full-model `SourceModel.luminosity` at an inner boundary **plus** explicit heating in the transported outer zone — can double count the outer material's radioactive energy, especially as the moving boundary recedes inward and more mass joins the transport domain. There is a second issue: the plan calls this an "Arnett interior," but xkn-diff is not just a simple Arnett one-zone source. Its thick-ejecta luminosity comes from a semi-analytic diffusion solution derived from Pinto–Eastman/Wollaeger, and xkn then adds the optically thin layers separately. So I would modify P1 before production.

## Better P1 architecture

Use the **xkn photosphere itself** as the moving inner boundary: **τ_grey = 2/3** rather than an arbitrary τ = 30 surface. Then the calculation becomes very clean: interior: inject only the xkn **thick-ejecta photospheric luminosity** L_thick(t) at R_ph(t), with its photospheric Planck spectrum; exterior: explicitly transport the shells outside R_ph; radioactive source in those shells: inject only the local deposited heating of the **thin/outer mass**; do not add a second full-ejecta luminosity source. This mirrors xkn's own physical decomposition: **L_total = L_thick + L_thin**, except that instead of xkn's approximate thin-shell emission, we let our line transport determine what happens in the outer material. That is actually a very nice experiment. It means: "We take the published xkn thick interior as the boundary condition and replace only its simplified treatment of the outer ejecta by resolved / expansion / line-binned transport." That is much more defensible than inventing a τ = 30 boundary.

## Why this also solves the moving-boundary problem

Your current plan says that when the boundary moves inward, newly exposed shells are initialized with a T⁴ V as added trapped radiation. I would avoid that unless you can derive that energy consistently from the xkn diffusion solution. Otherwise you're injecting an extra energy reservoir by construction. With the xkn photosphere prescription, the interior radiation is represented by L_thick, and the newly transparent material simply moves into the explicitly transported thin region. That is much easier to explain and much harder to double count.

## One more parameter issue: κ = 10 is not the exact published xkn secular model

The paper's Y_e = 0.20 secular benchmark uses the same mass and velocity you've adopted, M_sec = 2.64×10⁻² M⊙, v_rms = 0.06c, but the xkn comparisons do not uniquely correspond to κ = 10. The fixed-opacity secular case uses 22.3 cm² g⁻¹, while fitted xkn-diff values differ depending on whether bolometric luminosity or magnitudes are being fitted. So if you want P1 to be called a **published xkn benchmark**, choose one exact published opacity prescription and document it. I would prefer one of two choices: 1. **Physical fixed-opacity benchmark:** use the published composition-derived value, 22.3 cm² g⁻¹, if the goal is a clean fixed-Y_e state. 2. **Full xkn prescription:** reproduce the actual time/composition-dependent opacity used by the selected xkn-diff model. I would not quietly keep κ = 10 and call it the same xkn state.

## What remains good in your plan

The overall Paper IV manuscript plan is strong. I would keep: Part A exactly as written; manuscript infrastructure with frozen numerical macros; the methods-paper framing; the P1 three-leg fluorescence-only test; II-only at t ≥ 2 d for this robustness calculation; the free-streaming and slab-convergence tests; not implementing Phases 9–10 yet. The key scientific purpose of P1 is still just: **R_2 ≈ B_2 while B_bin,2 shows the larger chromatic deviation** for a multi-lanthanide physical composition. It does **not** need to become a perfect kilonova model.

## My recommended revision to Part B

Replace: outer zone above τ ≈ 30 + full SourceModel luminosity + zone heating with: **outer region above the xkn photosphere + xkn L_thick injected at the photosphere + explicit deposited heating only in the outer region.** Then define the temperature of those outer shells using the published xkn prescription, T_i(t) = T_ph(t) (1 − x_i²)/(1 − x_ph²(t)), rather than a separate ad hoc grey-temperature law. xkn explicitly uses this temperature continuation into the thin layers. That gives you a coherent published thermodynamic trajectory.

## Manuscript wording

If you make that change, I would describe P1 as: "a transport replacement test on the published xkn secular-ejecta structure: the xkn diffusion solution supplies the optically thick photospheric boundary condition, while the outer ejecta are evolved with resolved or coarse-grained line transport." That is much stronger than saying "xkn-like." So my decision is: **Approve Parts A and C; revise Part B before running production.** The current Part B is close, but the energy partition should follow xkn's own thick/thin decomposition. Otherwise the P1 light curve risks becoming internally inconsistent precisely in the way we're trying to avoid in this paper.

### Scope freeze and the review pass (2026-09-11, late; verbatim)

Yes — I think the project is now at the point where you should freeze the scientific scope and move into review/writing mode.

My calls are:

1. Title: use the safer, specific title

I would use:

Testing coarse line-opacity treatments with energy-conserving fluorescence in kilonova ejecta

That is stronger scientifically than the earlier "breaks the equivalence" wording because it does not overstate what the field previously assumed.

A slightly more assertive alternative, if the abstract stays careful, is:

Fluorescence differentiates expansion and line-binned opacity treatments in kilonova radiative transfer

I would probably submit with the first one unless the coauthors strongly prefer a punchier title.

2. F66's Gray bolometric result is not a failure

I would not describe it as "the zone is too thin to order."

What F66 actually says is more interesting: the closures are nearly indistinguishable bolometrically, while they produce different colours, with opposite signs.

That fits the methods-paper story very well.

For P1-xkn: expansion: i-K about 0.32 mag too blue; line-binned: i-K about 0.51 mag too red; peak luminosities differ by less than the MC resolution.

So the correct statement is:

No statistically resolved bolometric ordering is present in the P1 outer-ejecta benchmark, but the two closures produce distinct chromatic biases.

That is actually cleaner than finding another 15% luminosity difference. It tells us the closure choice can alter the SED without strongly changing the total light curve.

I would not rerun deeper shells simply to manufacture a bolometric hierarchy.

3. F65 stays Gray — and that is fine

The high-cap experiment gave exactly the right diagnostic result. The trapped line-binned packets simply consume more interactions when you give them a larger budget: 1786 → 4435 mean interactions per packet, while the trapped fraction barely moves. That shows this is not ordinary insufficient max_events. It's a structural limitation of attaching explicit fluorescence to that particular line-binned representation.

The important observation is that the observable light curve is fairly insensitive to increasing the cap, even though a small trapped population remains.

So I would report:

The line-binned fluorescence light curve is numerically stable at the quoted level, but a small non-terminating packet population prevents interpreting the calculation as fully converged in the microscopic transport sense.

That's more informative than forcing Green/Red. And keep the word Gray. Your preregistered rule did its job.

4. The paper's actual result is now quite clear

I think the manuscript should center around three statements.

First: Under complete thermal redistribution, resolved, expansion, and line-binned treatments remain relatively close. That's the Fontes connection.

Second: With explicit energy-conserving fluorescence, expansion opacity remains close to resolved transport in the tested states.

Third: Line-binned opacity develops larger and oppositely directed chromatic differences because binning has discarded information needed to determine post-absorption line behavior.

The P1 result is especially useful because it shows this isn't just a pure-Nd artifact.

And importantly: bolometric agreement ⇏ spectral agreement. That's now a legitimate conclusion again — but this time from corrected transport.

5. Do not add more major physics before submission

I would now not do: deeper P1 transport; full gas-energy equilibrium; NLTE; Phases 9–10; another abundance grid; another 27-point parameter survey. Those would make the paper harder to finish without materially improving the central claim.

The limitations are already explicit: prescribed thermal trajectory; radiation-field-driven macroatom is not self-consistent; II-only P1; Fontes light curve is a controlled extension, not an exact SuperNu reproduction; line-binned fluorescence has a small non-terminating packet population; some late bands are packet-starved.

That's acceptable for a methods benchmark paper as long as the claims stay inside those boundaries.

6. PR #3 should leave draft now

Yes. You now have: corrected transport merged; the high-cap rerun finished; F65 finalized; F66 finalized; 600 tests passing; make check green; no TODOs; frozen numerical macros; compiled manuscript; all major preregistered conditions reported, including Gray outcomes.

So: Mark PR #3 ready for review.

I would not merge it immediately, though. At this point the purpose of the PR is no longer "does the code run?" It should get one deliberate scientific/code review pass of:

1. the time-slab energy accounting;
2. the moving xkn photosphere bookkeeping;
3. the line-binned fluorescence interpretation;
4. the manuscript claims against F62–F66;
5. the figures/macros/freeze provenance.

Then merge.

7. I would freeze the scope after review

You have had several moments in this project where a more ambitious result disappeared under a better test. That is actually why the current paper is now much more trustworthy.

The final progression is: apparent catastrophic closure failure → transport bug found → corrected expansion opacity works fairly well → but fluorescence reveals a meaningful distinction between expansion and line-binned coarse graining.

That is a solid methods story.

So my final calls are:
* Title: use the safe "Testing coarse line-opacity treatments…" title.
* F66: treat Gray bolometric agreement as a result, not a failed hierarchy.
* F65: retain Gray because of the persistent trapped population.
* PR #3: mark ready for review now.
* Science scope: freeze after that review and finish the manuscript.

### Merge, the snapshot energy disclaimer, and the path to submission (2026-09-12, verbatim)

I'd make three calls.

For the snapshot escape fraction, keep the comoving bracket in the manuscript, because it is actually valuable: it shows the colour conclusion survives three different inner-boundary treatments. But do not let the escaped-energy fraction look like a physical observable of that static experiment.

Add one explicit sentence such as:

In the static snapshot, the absolute escaped-energy fraction is not interpreted physically because a moving thermalising boundary can exchange work with the radiation field without a corresponding time evolution of the ejecta; the boundary brackets are used only to test the robustness of the emergent colours.

Then keep the escape-fraction numbers in the record/SI if useful, but the main-text table/figure should emphasize the colour residuals. I would not drop the comoving bracket. Its agreement in colour is stronger evidence than hiding it because one non-observable diagnostic behaves oddly.

So: **Keep the comoving bracket; explicitly disclaim the snapshot energy fraction.**

The light curve is where energy accounting should be interpreted physically, because there the boundary motion, elapsed time, and expansion work are all represented consistently.

For the review and merge, the review you ran satisfies what I meant before merging. You now have an independent accounting/photosphere/provenance pass, your own claim-by-claim manuscript check, the corrected boundary physics, reruns, 602 passing tests, green structure checks, and a clean tree.

I would therefore: **Merge PR #3 now.**

I would still want a human scientific read before journal submission, but that is a different gate. It should not block merging the now-reviewed implementation and frozen results into `main`.

The human review I would want later is not another code audit. It is a referee-style read asking:

* Is the central claim narrower than the evidence?
* Is the distinction between snapshot and light-curve results obvious?
* Is "line-binned information loss" demonstrated rather than asserted too broadly?
* Are the Gray outcomes presented fairly?
* Would a reader mistake any controlled trajectory for self-consistent ejecta evolution?

That can happen on the merged manuscript.

On the table incident: you handled it correctly by catching it before commit and computing the actual values. It does reinforce one procedural point for submission, though. I would make the final tables and figure annotations machine-generated from the frozen records only, with no manually transcribed numerical rows. You are essentially already there with the macro/freeze infrastructure. Given that you also did the substituted-macro read-through and found several real issues, I would keep that as a permanent pre-submission check.

The comoving-boundary correction actually makes me more comfortable freezing the scope. It was a real relativistic physics issue at v ~ 0.1c, and after fixing it the exact amplitudes move (expansion i-K -0.32 -> -0.33, line-binned +0.51 -> +0.42) but the important qualitative result survives: **expansion and line-binned closures produce opposite chromatic biases**, while the bolometric P1 light curves remain unresolved from one another at the current statistical precision. Likewise, the three Fontes snapshot brackets give essentially the same scientific ordering: |D(r-K)|_bin >> |D(r-K)|_exp. That is exactly the kind of robustness check we wanted.

So I would not reopen the physics program. My final state would be:

1. Add the explicit static-snapshot energy disclaimer.
2. Merge PR #3.
3. Tag/freeze the Paper IV computational result set.
4. Do a human/referee-style manuscript review.
5. Revise wording and presentation only unless that review exposes an actual scientific flaw.
6. Submit.

At this point, another simulation is more likely to expand the scope than materially strengthen the paper. The important thing now is making sure the paper says exactly what the corrected experiments establish, and no more.
