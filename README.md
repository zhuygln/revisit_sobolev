# revisit_sobolev

Kilonova radiative-transfer codes cannot resolve the millions of narrow
lanthanide absorption lines that shape a kilonova's light, so they all use a
shortcut. That shortcut is usually called "the Sobolev approximation." This
project shows it is really **two** approximations with error budgets differing
by roughly an order of magnitude, and measures both against resolved
calculations on calibrated atomic data.

**The result.** Against frequency-resolved transport as truth, on realistic
La II line forests:

| treatment | band-flux error (vs a deterministic finite-profile reference on identical rays) |
|---|---|
| Sobolev approximation proper (`e^−τs` per line) | **≲1%** at v_D ≤ 30 km/s for τ_max ≤ 5, **≲0.3%** at v_D ≤ 10 km/s — a boundary effect ∝ v_D, vanishing at thermal widths |
| expansion opacity (its usual implementation) | **+38%** at thermal widths (the pure Poisson-vs-Bernoulli gap), **+45%** in the Monte Carlo implementation at 100 km/s |

The gap comes from a per-resonance substitution `τ → 1−e^−τ` that belongs to
expansion opacity alone — not to Sobolev's locality or isolation assumptions —
and it does not vanish as line widths shrink. The Sobolev approximation's own
residual, by contrast, switches off entirely at physical line widths: it is a
finite-region boundary effect, not a failure of either Sobolev assumption. Across 36 conditions spanning
four wavelength windows, three epochs and three ion mixtures, the separation
holds, with the realized maximum line optical depth as the controlling
variable at fixed line width and geometry.

## The paper

**[docs/paper/manuscript.pdf](docs/paper/manuscript.pdf)** — 21 pp in MNRAS
format, readable directly on GitHub. Written for a reader with no
radiative-transfer background: a primer builds every concept from scratch and
an appendix derives every formula.

Source is [manuscript.tex](docs/paper/manuscript.tex); rebuild with
`cd docs/paper && make` (needs pdflatex + bibtex — see
[docs/sedona/SETUP.md](docs/sedona/SETUP.md) §4; conda-forge's texlive-core
does not work).

**[docs/paper3/manuscript.pdf](docs/paper3/manuscript.pdf)** — Paper III,
Nature Astronomy Article format: *Coarse-grained line opacity leaves a
detectable chromatic signature in kilonovae that ejecta parameters cannot
mimic.* Every number in it is a macro generated from
[paper3/FROZEN.json](paper3/FROZEN.json) (`docs/paper3/latex_tables.py`),
every figure from the frozen JSONs (`docs/paper3/display_items.py`); `cd
docs/paper3 && make` rebuilds it and runs a structure check that bans literal
result numbers in the prose. The affiliation and repository URL are still
placeholders.

## Findings

| # | Finding | Detail |
|---|---|---|
| F1 | Gradients symmetric about a resonance cancel the Sobolev error at leading order, so gradient-based validity diagnostics are conservative | [§4.2](docs/results_report.md) |
| F2 | Resolved transport cost 378× expansion opacity on a standard example | §4.4 |
| F3 | The O(v_bulk/c) offset is an artifact of frozen-snapshot integration, not relativity — superseded by F11 | §4.5, §4.14 |
| F4 | Expansion opacity attenuates by `e^−(1−e^−τ)` per crossing, not `e^−τ` | §4.5 |
| F5 | That error is per-resonance: it does not average away with line count, only as τ→0 | §4.6 |
| F6 | Δ_expansion = a strength-set floor plus a width-growing term; the floor survives to thermal widths | §4.8 |
| F7 | Δ is non-monotonic in forest density — maximal for strong *sparse* forests, suppressed at full blanketing | §4.10 |
| F8 | Solver-vs-Monte-Carlo offsets traced to the thermal-emission convention, not physics; matched, the codes agree to 0.1% on the La II forest — inside the measured 0.3% Monte Carlo noise | §4.11 |
| F9 | The strength floor is expansion opacity’s alone; Sobolev’s own residual is small and width-driven (see F12) | §4.12 |
| F10 | The separation is universal across windows, epochs and ion mixes; realized τ_max controls it at fixed v_D and geometry | §4.13 |
| F11 | Neglect of light-travel-time evolution is a distinct approximation: frozen τ/τ_S = (1−β)/γ vs worldline 1/γ, which has no O(β) term. SEDONA in steady-iterate mode solves the frozen problem, so transport treatment is a third convention cross-code comparisons must match | §4.14 |
| F12 | Overlap is inert in pure absorption (optical depths add exactly); the Sobolev residual is a finite-region boundary effect ∝ v_D/Δv_shell, negligible at thermal widths | §4.15 |
| F13 | The expansion-opacity error is bin-width invariant from 0.025 to 41 lines per bin — intrinsic to the formalism, not a usage artifact | §4.16 |
| F14 | Neither candidate reference code can answer Paper II alone: SEDONA has no line branching at all, TARDIS has no expansion opacity. A ~250-line branching Sobolev MC, validated to ≤2.1σ against the analytic Sobolev leg, supplies the same-code differential both lack | §4.17 |
| F15 | **The mechanism, made exact (referee revision).** Expansion opacity preserves the expected interaction count per crossing E = Σ(1−e^−τ) identically — Karp's mean free path — and applies Poisson survival e^−E to what is a Bernoulli product e^−S; F_exp/F_Sob = ⟨e^D⟩ weighted by transmission, D = S−E, to 1e-13. F4, F5, F7, F13 are corollaries | §4.18 |
| F16 | Δ_Sobolev against a deterministic reference on identical rays is +3.1% at v_D = 100 km/s in all four transport modes, +0.3% at 10, +0.03% at 1 km/s; the breadth median was a stale normalization (+3.65% → +0.26%); SEDONA resolved validates the reference to ±0.3% | §4.18 |
| F17 | Seed-matched SEDONA pairs correlate at +0.95–0.999; paired Δ_exp scatter 0.02–0.28% against 0.2–0.4% from quadrature. Headline over 10 seeds: +44.90% ± 0.04 | §4.18 |
| F18 | Paper I's radiative-equilibrium check is confined to SEDONA's transport window (re-emission is drawn from the emissivity on the grid), so its +7.7%/+5.0% differentials are bookkeeping within one 230 Å window, not an emergent-band result; retracted as such | §4.18 |
| F19 | **Escape-probability branching (Paper II).** A photon re-emitted in a resonance zone escapes it only with β = (1−e^−τ)/τ and is otherwise re-absorbed and re-drawn, so radiative branches compete as A·β, not A; the chain's exit distribution is A_uj β_uj/Σ A β exactly (tested). Phase 0 lacked it; the SEDONA RE comparison caught its absence in one look | §4.19 |
| F20 | **Fluorescent optical refill.** Full La II, 6000 K continuum: direct A-branching raises the 3800–3955 Å emergent flux from 0.183 to 0.660 ± 0.003; the pumps are 3300–4500 Å (the far-UV contributes nothing), through 5d6p upper levels exiting the strong ground-connected lines, 971 pathways with the top 10 carrying 25% | §4.19 |
| F21 | **Closure sign reversal.** Expansion opacity is too transparent in pure absorption (+88%) but, with complete thermal redistribution (ε = 1, re-emitting from its own κ_exp B_ν), underproduces the same band by 38% relative to direct branching; at ε = 0 it overproduces it by 34%. With line identity carried through the bin (`expansion_branch`, Poisson absorption + exact A·β exit kernel) the opacity-representation error under fluorescence is +21% (density-limited — F24) | §4.19, §4.20 |
| F22 | **No scalar ε reproduces La II fluorescence (outcome B).** Sweeping the two-level-atom thermalisation parameter on the same opacity as the physics leg, ε_best is 0.02–0.06 in the red and blue, 0.36 in the UV, and the optical 4500–6000 Å band is outside the TLA's reach at any ε; the redistribution matrices show why — branching is a blue → optical channel, thermalisation a blue → red one, and ε interpolates between them. The expansion leg's ε_best ≈ 0.3 in the 3800–3955 Å band is a compensation with its +21% opacity error and does not carry to other bands | §4.20 |
| F23 | **The closure verdict survives 0.1c; the frozen shortcut does not.** With worldline-consistent transport (per-packet clock, exact Doppler, epoch-diluted τ, aberrated re-emission) on a 0.05–0.15c shell, no scalar ε reproduces branching either (ε_best 0.20–0.49 by band, red unreachable) — but calibrated ε values shift by 0.15–0.27 vs the slow shell and the redistribution structure changes, so v_bulk is a validity-map axis; a frozen first-order snapshot overstates blue transmission 3.8× at 0.1c and is unusable for closure work | §4.21 |
| F24 | **Outcome C, and the closure's density limit.** ε_best is ion-dependent: La → Ce shifts of 0.11–0.24 per band with reachability flipping in opposite directions, and the La+Ce blend tracks the dominant (Ce) forest — composition sets the calibration. The branching-aware Poisson closure is density-limited: +21% on La II's 949-line forest, +113% on Ce II's 22,960-line one, where per-bin saturation clipping leaves the closure three orders of magnitude too transparent in a band the true opacity blacks out | §4.22 |
| F25 | **The middle method works (Paper III).** A group-to-group redistribution operator with discrete exit-frequency tables reproduces explicit A·β branching on the same opacity: La II to ≤1.6% in every band with a 4×4 matrix, Ce II to 2.9%/0.7% at 32/64 groups, bolometric to ≤0.5%, tables of tens–hundreds of kB. Within-group re-emission must be discrete: a continuous PDF double-counts self-absorption with an error that grows under refinement | §4.23 |
| F26 | **The kernel's state space is (T_gas, τ_scale, ion) — for La II.** Source spectrum transfers freely (≤1.4%, 4000–8000 K); the epoch axis collapses exactly onto τ_scale (a τ-matched 1 d kernel equals each epoch's own kernel, while the fixed kernel fails at 13% where τ_max = 34); T_gas is the one genuine axis (3000→5000 K transfer errs 9.6%, recomputed ≤1.5%). The source-spectrum half does not generalize — see F28 | §4.24, §4.25 |
| F27 | **Gate 2: compression is generic (outcome A).** A discrete-table R_ij reproduces explicit branching for all three ions on the same opacity: La II 1.62% and Ce II 7.48% at N_g = 4, Ce reaching 2.89%/0.74% at 32/64 — and **Nd II 0.44% at N_g = 4**, the easiest of the three despite 3,336,077 transitions, because the fixed τ_max = 5 normalization leaves it only 4,496 opacity lines against Ce's 22,960. The blue→blue block separates the hard ion from the easy ones but does not order them | §4.25 |
| F28 | **The kernel's state space is itself ion-dependent.** F26's structural claims carry to Nd II — the epoch axis still collapses onto τ_scale exactly, T_gas is still the genuine axis (6.4–19.2% fixed, ≤0.51% recomputed). Its third claim does not: Nd's fixed-kernel error is monotone in T_src and flips sign across the training point (+12.44 → −4.88% over 4000–8000 K) where La's is flat at ~1%. Denser groups draw a more evenly weighted absorbing-line mix, so reweighting the continuum reweights the rows | §4.25 |
| F29 | **The composition rule works at the 5% level (P9).** An opacity-weighted mixture of per-ion kernels — weights from the blend's opacity alone, no blend run — reaches 4.27% worst-band on La+Ce at N_g = 64 where a blend-trained kernel gets 1.37%, beating the best single-ion control by 2.4×. The gain is located: `ce_only` fails at −10.3% in the optical, the blue → optical branching channel that La carries despite being 5% of the opacity, and the rule repairs it to −0.7%. Composition therefore leaves the state space at the 5% level — a per-ion library suffices — while explicit blend training is still needed below ~2% | §4.26 |
| F30 | **The opacity is the binding constraint, not the redistribution (P11).** With the identical R_ij, exact line opacity errs 0.92% (La II) / 2.21% (Ce II); grouping the opacity takes that to 14–18% / 91–127% whichever single-scalar rule is used — Στ too opaque on La (−14.5%), Σ(1−e^−τ) too transparent (+17.8%), both far too transparent on Ce. F15 as a design constraint: one scalar per bin cannot carry both the attenuation and the interaction count, and scattering needs both. The κ_grouped + R_ij target architecture fails on dense forests; a bin carrying both is the candidate A bin carrying **both** quantities (survival from S, line draw from p) was then tested: it works on La II (21.32% → **8.66%**, saturated band +21.3% → −0.7%) and fails on Ce II (112.86% → **139.27%**), where more opacity makes the band *brighter* because it is refilled by fluorescence faster than absorbed — redistribution-limited, not attenuation-limited. And `dual_group` is bit-identical to `binned_group`: a pure R_ij closure never draws a line in the bin, so it cannot use the second quantity at all | §4.27 |
| F31 | **One remembered line is not the missing state — and the failure splits by line spacing.** Carrying exactly one extra number per packet (the frequency last emitted at, crediting that line's τ to the next free path) buys a factor 2.2 on La II (14.49% → **6.50%**, saturated band −14.5% → −6.3%), the cheapest grouped-opacity repair found. On Ce II it moves 126.66% → 116.31% and 91.29% → 81.66% — real but no rescue. The events/packet counter shows memory removes a comparable fraction of the excess interactions in both forests (39% La, 34% Ce), so the Ce error is not driven by excess interactions: it is the fluorescent refill of §4.19–4.20, and no local interaction bookkeeping reaches it. Sparse forests need one remembered line; ~~dense forests need the resonance *sequence*~~ — that second clause is **retracted by F33**, which finds no benefit from depth beyond m = 1 on Ce | §4.28, §4.30 |
| F32 | **The redistribution operator is local in frequency, not low-rank — so "few modes" is the wrong explanation for its compressibility.** Effective dimension never saturates: participation ratio grows as N_g^0.64–0.75 across La/Ce/Nd with PR/N_g falling 0.5 → 0.2, and NMF rank-8 of a 25-row operator still misses row-L1 0.47 (23% total variation per row). Rank *anti*-correlates with compressibility — Ce has the lowest energy-operator dimension (PR 1.60) and is the hardest to compress. Decisive at matched parameter count: 16 numbers as a coarse 4×4 matrix give 1.62% where 912 numbers as a rank-16 factorization give 11.30%. Coarsening averages neighbouring groups; truncation projects onto modes; only the first works. This unifies Paper III — redistribution is smooth at the group scale so it coarse-grains (F25/F27), the opacity is a comb whose ordering decides a packet's fate so it does not (F30/F31) | §4.29 |
| F33 | **Resonance-sequence depth is not the missing information, and the density limit is not about density (retracts F31's second clause).** Making memory a depth rather than a switch: La converges by m = 4 and Nd by m = 8, each gaining ≤0.5 points beyond m = 1, and **Ce gains nothing** — 116.31 → 116.79% across a sixteen-fold increase in remembered history. Memory is a between-step correction; the dense-forest failure is within-step. Nd II, run through P11 for the first time, breaks the density reading outright: 4.7× La's opacity lines and **1/12 its error** (expansion + A·β 1.79%, a working closure). What orders the three ions is **band-local saturation** — saturated lines inside the failing band (Nd 1, La 4, Ce 24), or Στ there (11.4, 22.4, 89.6) — not total line count. Three points cannot fix the exponent; that is what the synthetic phase diagram is for | §4.30 |
| F34 | **Band-local saturation controls the grouped-closure failure; redistribution does not.** 96 synthetic forests with independently dialled crowding, saturation, spacing and redistribution range: Spearman ρ = +0.91 for Στ in the band and +0.86 for N_sat there, against **+0.25 and −0.31 for the two redistribution axes** — an independent confirmation of F33 from the opposite direction. The synthetic family collapses as ΔF = 0.162·N_sat^0.58 (scatter ×1.95), and the real atoms sit at ratio 0.71 / 0.40 / 1.25, two of three inside that scatter. Building the forests required first measuring the real τ distribution *inside* the failing band: mostly weak lines with a saturated tail, ln-spread 1.7–2.05 across all three ions. Band-to-forest geometry is eliminated as the residual cause (10× change → 16% effect); the emergent-cascade `ladder` forests match the real interpolation at matched N_sat (59–62% vs 62%) where the dialled ones overshoot (78–94%). A partial collapse: not yet a general law | §4.31 |
| F35 | **The closure error changes sign — there is a phase boundary, not a scaling law (supersedes F34's power law).** Thirteen GSI ions under a uniform normalization do not collapse: La II and Pr II have identical band saturation (S = 13.4 vs 13.8) and differ 5× in error (6.55% vs 31.39%), while Ce II at 5× their saturation errs *less*. A density scan shows why: the binned closure is **too opaque at low density and too transparent at high density**, crossing zero for Ce II between S ≈ 45 and 67 then rising +21% → +125% across a density factor of 1.33. La II (−5.0%) and Pr II (−31.4%) are the same sign at matched saturation, so their 5× difference is within-regime scatter, not a sign flip — that remains open. Also: the project's τ_max = 5 *window* normalization is ion-specific by accident and diverges for most ions (Yb II demands n_ion = 1.7×10¹², β = 1.5×10⁻⁸, and the branch chain cannot terminate) — universality claims need the global normalization used here. On real ions saturation and redistribution range are confounded (ρ = +0.75 vs +0.77), which is why the decorrelated synthetic experiment (F34) is what identifies the cause | §4.32 |
| F36 | **The normalization audit: three cross-ion claims revised.** Re-measuring five ions at matched line strength (`global_tau_max`) instead of the accidental window recipe: **F27 is strengthened** — every ion compresses at *four* groups to ≤4.3%, Ce II included, so "dense ions need 32–64 groups" was an artefact. **F24's density limit inverts** — the branching-aware Poisson closure is +14.7% on La II and **+1.9% on Ce II**, making Ce the better case, not the catastrophic one. **F33's null memory result is superseded** — memory is the most effective correction found (Pr II −31.4→−5.6%, Ce II +12.2→**+0.2%**); its null was a property of an over-dense Ce. **F30's structure survives on five ions**: redistribution 0.2–1.8%, grouped opacity −31% to +15%. Memory's *direction* now has a mechanism — it always adds transparency, which brightens an absorption-limited band (La −5.0→+4.5%, overshooting zero) and dims a refill-limited one (Ce +12.2→+0.2%) | §4.33 |
| F37 | **The synthetic model contains only one side of the boundary — which isolates fluorescent refill as the cause of the other.** Scanning τ across S = 2.7 → 1459 at three redistribution ranges, all 36 controlled conditions are negative and deepen monotonically to −99%; no ΔF = 0 crossing exists in the model. The reference transmission shows why: at matched saturation the synthetic band is 2–3× more opaque than any real ion (S ≈ 55–90: synthetic 0.23–0.25 vs Ce II 0.58), because real forests **feed** the band from outside — Tm II transmits **1.049**, more than the continuum entering it, and Dy III 0.970. `synthetic_forest`'s exits sit at a fixed offset from their own line and carry no opacity, so it redistributes locally and never delivers net inflow; the band only darkens. E4 is therefore blocked on a model deficiency, not a measurement — and the failure is itself the controlled demonstration that fluorescent refill is what produces the too-bright branch. **Fix attempted:** delocalizing the exits (a fraction placed anywhere in the forest) cuts the error 3× (−77.5% → −25.3%) and raises transmission 0.226 → 0.292, and does produce a crossing — but running the **wrong way** (+4.0% at S = 8.1 → −7.9% at S = 54.1, the reverse of every real ion). Necessary, not sufficient; exit lines carrying their own opacity is the next change | §4.34 |
| F38 | **The two approximations are not additive, and their interaction can flip the sign of the total.** Counterfactual legs isolating each — exact opacity + grouped redistribution (A), grouped opacity + exact A·β (B), both (C) — give |A| ≤ 1.4% on every ion, so the redistribution approximation contributes almost nothing (F30 by a third route, measured directly rather than inferred). B dominates, and for the least saturated ions C ≈ B. But at S ≳ 14 an interaction term appears, −6.4% (Pr II) and −4.3% (Ce II), large enough for **Pr II to run +4.3% too bright under the opacity approximation alone and −2.0% too opaque with both**. The zero of §4.32 is therefore not where the dominant error changes sign but where B plus the interaction does — so a closure whose pieces are separately validated can fail, or appear to succeed, for reasons neither piece shows alone | §4.35 |
| F39 | **Recurrent exit opacity recovers the correct boundary orientation — the sign change now appears in three independent settings.** §4.34 predicted the missing physics: exits terminating on unpopulated levels cannot cascade, so the model refilled the band *less* as saturation rose where real forests refill *more*. Giving exit lines their own opacity on shared populated levels flips the boundary the right way: at exit_tau = 0.5 the binned closure runs **−62.3% → +2.5% → +39.1%**, crossing neg→pos at S = 179 (and at S = 688 for exit_tau = 2.0), where terminal exits gave no crossing at all or one running backwards. The boundary is therefore reproduced by a density scan, a 13-ion survey, and a controlled forest that separates the axes real atoms confound. Not claimed: the crossing *location*, which moves with exit_tau, or the expansion leg, which still fails | §4.34b |
| F40 | **A realistic kilonova crosses the cancellation boundary at 1.2 days.** Homologous ejecta (ρ ∝ t⁻³, T ∝ t⁻¹ᐟ², X_lan = 0.1) sweep band saturation across four orders of magnitude in n_ion, and the practical grouped closure runs **+64.6% too bright at 0.5 d → zero at 1.17 d → −28.4% too opaque at 1.5 d** — ninety points across a factor of three in time, straddling the epoch kilonova spectra are taken. The crossing occurs at **S = 47.5**, matching the boundary located independently at S ≈ 50 by a density scan, a 13-ion survey and a controlled synthetic forest: a fourth confirmation from a different construction. Stated carefully, the zero here is the *opacity* error changing sign (B crosses at 1.21 d) rather than two large errors cancelling, with |A| ≤ 2.1% throughout; the cancellation mechanism is separately visible at 2 d, where the binned closure reads −1.1% while its opacity piece alone reads −4.1%. Either way: **near-zero residual at one epoch is not evidence a closure is correct** **La II on the same history sharpens it**: at 0.75 d the expansion closure reads **+0.1%** while the binned closure reads **−55.7%** — same epoch, same ejecta, same atom, differing only in whether a bin carries Σ(1−e^−τ) or Στ. One looks exact, the other is wrong by more than half. | §4.36 |
| F41 | **The closure error is chromatic, not bolometric, and the diagnostic band residual is a proxy for neither.** Converting §4.36's ejecta history into absolute L_ν, escaping luminosity and AB magnitudes: the grouped-**redistribution** approximation is invisible to an observer — worst \|Δm\| of **0.006 / 0.008 / 0.008 mag** on Ce II, La II and a four-ion blend — so F25/F27's kernel compression costs less than a photometric error bar. The grouped-**opacity** approximation is five times larger in colour than in luminosity: Ce II at 0.5 d is **0.14 mag too bright bolometrically and 0.74 mag wrong in g−r**, because the closure moves flux between bands rather than creating or destroying it. And the 3800–3955 Å residual carrying every result from §4.23 to §4.36 tracks neither — **−59.5% in band is +0.007 mag bolometric** (La II binned, 0.75 d), −26.0% is photometrically invisible (Ce II, 2 d), and +55% is 0.46 mag in g (Ce II, 0.5 d). Two defensible groupings give **opposite colour errors 0.15 mag apart at one epoch on one atom**. Every magnitude is a floor: at fixed saturation the photometric error grows **34× (La) to 67× (Ce) from Paper I's 0.01c to a kilonova's 0.3c** while Δm_bol stays under 0.011 mag, and it survives the worldline transport treatment, which changes the band residual by a factor of 4 and the magnitudes not at all. On a **physically normalized** kilonova — M_ej = 0.01 M⊙, v = 0.05–0.2c, X_lan = 10⁻³, ρ derived from the mass rather than tuned, worldline — the practical closure is **0.06 mag bolometric and 0.63–0.77 mag in a band or colour**, and at 3 d its binned variant is bolometrically **exact** while 0.74 mag wrong in r−i. The band diagnostic runs out of signal from 4 d on, exactly where the photometric error is still tenths of a magnitude | §4.37 |
| F42 | **The chromatic closure error survives real filter curves (Gate 1 passes), and the three notebook-only findings now have drivers and data.** DECam g r i z + 2MASS J H Ks from SVO, applied to the stored F41 spectra: every ≥ 0.6 mag top-hat colour error stays at 0.65–0.77 mag; the A_redist floor stays ≤ 0.009 mag. Ce II density scan, F38 table and F39 exit-τ scan reproduced. | [§4.38](docs/results_report.md) |
| F43 | **On a heating-powered kilonova the grouped-opacity closure's colour error is 1–3 mag at every point of a 27-model (M_ej, v_ej, X_lan) grid, and its sign is uniform: too blue.** Worst live colour error per model 0.96–2.84 mag (C_both) against an A_redist floor of 0.02–0.13 mag; the closure's g is 0.2–2.1 mag too bright and its K up to 3.6 mag too faint. Every leg has the same L_bol by construction; the harness makes no bolometric statement at S ≳ 10⁴. | [§4.40](docs/results_report.md) |
| F44 | **Gate 2: the closure error is not degenerate with (M_ej, v_ej, X_lan) — class C-B (distinct residual) at all 24 analysable grid points, every opacity closure, every robustness variant.** χ²_RT/N = 28–549; residual fraction R = 0.46–1.00 after the best three-parameter shift; A_redist is undetectable (C-C) at 23 of 24. A fit measures a residual, not a bias. | [§4.40](docs/results_report.md) |
| F45 | **F43/F44 survive the floor mask and the core convention (chain cap: §4.44); a free luminosity history absorbs the residual as a *class* but not as a *misfit*.** Erratum: the §4.40 floor mask was a no-op (27 floored rows); re-baselined, Gate 2 is C-B at 20 of 21 analysable points. One free grey magnitude per epoch takes median R to 0.28 and 12 of 19 points to C-A — with 1–7 mag luminosity offsets and a leftover χ²_res/dof of 23. | [§4.41](docs/results_report.md) |
| F46 | **Gate 3: every scenario detects the closure error at every eligible point (χ²_RT,obs/N = 41–1005; 30–40σ in single bands at real errors) and it survives the ejecta parameters everywhere (18/18 dense, 11/11 sparse, 16/16 optical); under a free luminosity history its survival is set by the NIR — 9/18 with six-epoch JHK (all four X = 0.1 points), 3/8 with two NIR epochs, 1/15 without NIR.** Counts on the complete grid: §4.44. | [§4.42](docs/results_report.md) |
| F47 | **The Planck temperature proxy is validated as a *gas*-temperature direction (cosine 0.92 with the measured transport response, lever arm 1.35× larger), the illumination temperature alone does not reach the observer (‖d_T^MC‖ = 0.06 of the proxy, every band within the 0.13 mag noise floor), and the central point stays C-B with the measured direction (R 0.37).** The closure residual is not a photospheric-temperature error. | [§4.43](docs/results_report.md) |
| F48 | **With the grid complete (162 of 162 cells), Gate 2 is C-B at 27 of 27 points (median R 0.83, χ²_res/dof 118) and Gate 3 at 26/26, 18/18, 25/25 eligible points; a free luminosity history absorbs the residual at 8 of 9 lanthanide-poor points and 0 of 9 lanthanide-rich ones, and at real errors the residual survives it at 16 of 17 `dense` X ≥ 10⁻² points.** The chain cap moves single bands by 0.14–0.21 mag at the four worst-trapped cells (their own noise floor, non-monotone in the cap); class kept at all 27 points, colour signs kept, the < 25 % magnitude criterion met at 4 of 12 colours there. | [§4.44](docs/results_report.md) |
| F49 | **The closure error is a coherent, signed, one-mode pattern the size of the σ_sys allowance:** C_both exceeds 0.5 mag at 56 % of 524 live observables (18 % beyond 1 mag), (g < 0, K > 0) at 39/39 coepochal pairs, 0.80 of its squared norm in one band-epoch mode against a 0.33 sign-scrambled null and 0.31 for A_redist; median χ²/N 0.56 against a 1 mag allowance, 2.26 against 0.5 mag. A consistency statement about the closure experiment, not a claim about any published fit. Every number the paper quotes is frozen in `paper3/FROZEN.json` (`freeze.py --check`, tag `paper3-freeze`) | [§4.45](docs/results_report.md) |

Full write-up with figures and numbers:
**[docs/results_report.md](docs/results_report.md)**.

## Layout

```
sobolev/       the package: constants, line profiles, optical depths,
               Boltzmann populations, GSI parsing, the deterministic
               formal solver, and the analytic Sobolev leg
experiments/   one directory per experiment; generators and comparison
               scripts committed, SEDONA run outputs gitignored
notebooks/     Phase 0: single-line toy model, GSI line spacing
paper2/        Paper II: Phase 0 (SEDONA audit, TARDIS record, three-level
               branching MC) and Phase 1 (whole-ion Sobolev/expansion MC
               with fluorescence, the La II measurement)
paper3/        Paper III: the reduced redistribution closure (kernel,
               frozen reference, compression sweep, T and epoch transfer),
               the kilonova closure grid (phase12/13) and freeze.py +
               FROZEN.json, the frozen analysis the manuscript quotes
tests/         347 tests pinning the physics of every module
docs/          results report, lab notebook, planning inputs, paper/
data/          raw atomic data (gitignored; provenance in data/README.md)
outputs/       working figures (gitignored; committed copies in docs/figures/)
```

## Setup

Needs Python >= 3.10 and numpy >= 2, so a distro 3.8/3.9 will not do.

```bash
python3 -m venv .venv              # or: conda create -p .venv python=3.12
.venv/bin/python -m pip install -e ".[dev]" h5py
.venv/bin/python -m pytest         # 347 passed (2026-09-03)
```

Atomic data is not committed — see [data/README.md](data/README.md) for the
Zenodo record and re-download instructions. Tests that need it skip cleanly
when `data/` is empty, so a data-less clone reports passes and 5 skips.

SEDONA lives outside this repo and is needed only by `experiments/`
(everything in `paper2/` and `paper3/` is pure Python). Its no-root WSL2
build recipe is in the lab notebook §5, with the working `Makefile.wsl` and
the full machine-setup walkthrough — Python, data, SEDONA, LaTeX — in
**[docs/sedona/SETUP.md](docs/sedona/SETUP.md)**.

The SEDONA experiment drivers find the binary through `SEDONA_HOME`
(default `~/personal/pubsed`) or `SEDONA_EXE`.

Per-experiment reproduction commands are in
[docs/results_report.md](docs/results_report.md) §7.

## How this got built

The project ran as a staged feasibility study — a weekend toy model first,
then cross-code validation, then real atomic data, then parameter maps — with
a go/no-go decision at each stage. That history, including the dead ends and
the bugs that mattered, is in
**[docs/lab_notebook.md](docs/lab_notebook.md)**; the organized results are in
the report §4. The original plans are preserved unmodified as
[docs/babystep_plan.md](docs/babystep_plan.md) and
[docs/research_requirements.md](docs/research_requirements.md).

## Status and open items

Phase 0 through the validity maps is complete, and the manuscript has been
revised for a major-revision report (August 2026): the thesis is now that
the two approximations carried under one name preserve different properties
of a line forest (F15), Δ_Sobolev is formed against a deterministic reference
(F16), every Monte Carlo number is a seed-matched pair (F17), and one
radiative-equilibrium check is included (F18). The response letter is
`docs/paper/response_to_referee.md`. Outstanding:

1. ~~The residual Sobolev error.~~ **Resolved (F12).** Most of it was a
   normalization artifact; the remainder is a finite-region boundary effect
   that vanishes at physical line widths. Overlap was excluded both
   analytically (optical depths add) and numerically.
2. **Realistic velocities.** All results use 1000–3000 km/s shells. The
   machinery now reaches 0.1–0.3c: both the solver and the analytic leg carry
   a co-evolving medium and worldline transport. Two corrections came out of
   building it — the dominant term is light-travel *dilution* ((1−β)² = 51% at
   β = 0.3), not the 4.6% relativistic one, and the CD/CP resonance surfaces
   stay irrelevant under worldline transport, where the locus is linear.
   A pilot finds Δ_Sobolev β-independent but Δ_expansion **growing** with β;
   it is synthetic-forest only and not yet a finding (lab notebook §9n).
3. **Scattering and fluorescence** — Paper II. Phase 0 built and calibrated
   the instrument (F14); Phase 1 (`paper2/phase1/`) carries the whole La II
   ion through eight treatments in one vectorized code: fluorescence, not
   thermal re-emission, refills the optical band, and an expansion-opacity
   code with complete thermalisation lands 38% below the physics once it
   does (F19–F21). Phase 2 answered the decisive question: **no scalar
   thermalisation parameter ε reproduces direct La II fluorescence** across
   the spectrum (F22) — ε_best differs by band and one band is unreachable —
   whereas a closure that keeps line identity through the bin gets within
   21%. The verdict survives realistic bulk velocity: at 0.05–0.15c with
   worldline-consistent transport no scalar ε works either, while the
   calibrated ε values shift by 0.15–0.27 — and a frozen-snapshot code at
   0.1c overstates blue transmission 3.8× (F23). Ce II makes it outcome C:
   ε_best is ion-dependent, the La+Ce blend tracks the dominant forest, and
   the branching-aware closure itself is density-limited (+21% on La II,
   +113% on Ce II, whose 22,960-line forest the Poisson opacity cannot
   black out) (F24). Paper II's identity is that result; manuscript text
   and NLTE come afterwards. Paper III (`paper3/`) opens the constructive
   question — how much redistribution information is actually required? —
   and both its gates are passed. A discrete-table R_ij reproduces
   explicit branching on the same opacity (F25), and it does so for every
   ion tried: La II and Nd II at four groups, Ce II at 32–64, which is
   Gate 2's outcome A (F27). The epoch axis collapses exactly onto
   τ_scale and T_gas is the one genuine thermodynamic axis (F26), both
   ion-independent — but whether the source spectrum can be dropped is
   not: it transfers freely for La and fails at 12% for Nd (F28), so
   T_src must be checked per ion before a kernel is tabulated. Composition,
   by contrast, does leave the state space at the 5% level: an
   opacity-weighted mixture of per-ion kernels reaches 4.27% on the La+Ce
   blend with no blend training, beating the dominant ion alone by 2.4×
   (F29). P11 then put the two halves together and returned a negative
   result worth having: with the identical kernel, exact line opacity errs
   ≤2.2%, but grouping the opacity costs 14–18% on La II and 91–127% on
   Ce II by either single-scalar rule (F30). The redistribution half was
   never the hard half — one scalar per bin cannot carry both the
   attenuation and the interaction count, and scattering needs both. So
   κ_grouped + R_ij is not usable on dense forests. A bin carrying both
   quantities was then tried: it works on La II (21.3% → 8.7%) and fails
   on Ce II (112.9% → 139.3%), where the deep band is redistribution- not
   attenuation-limited — and a pure R_ij closure cannot use the second
   quantity at all, since it never draws a line in the bin. Restoring
   line identity at emission is the remaining lever, and it costs the
   thing grouping was meant to buy. E1 then asked why the redistribution
   half compresses at all, and the answer is not "few modes": the
   operator's effective dimension never saturates (PR ~ N_g^0.65) and
   rank truncation fails badly, while 16 numbers as a coarse matrix beat
   912 as a rank-16 factorization (F32). Redistribution compresses
   because it is *smooth* in frequency, not low-rank — and the opacity
   does not because a comb of resonances is not smooth at any group
   scale. That is one statement covering both halves. E2/E3 then tested the
   resonance-memory reading and it did not survive: memory depth does
   nothing on the dense forest (F33), the redistribution axes are null in
   decorrelated synthetic forests (F34), and thirteen real ions do not
   collapse onto any scaling law because **the closure error changes
   sign** — too opaque at low density, too transparent at high, with a
   boundary between (F35). What the grouped closure lacks is set by
   opacity structure, not by redistribution or by history.
