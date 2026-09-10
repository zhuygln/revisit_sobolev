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
