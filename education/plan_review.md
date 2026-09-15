# The education subproject — PI decisions, verbatim

The rule of this file is the same as `paper4/plan_review.md`: the PI's
decisions are recorded in their own words, in order, before they are executed.

## The specification (2026-09-15)

Yes. I would make this a **separate educational subproject inside `revisit_sobolev`**, designed as a miniature textbook that is executable from first principles and ultimately builds into the PRL calculation.

The important rule is: **educational demos ≠ publication evidence.**

The demos teach and validate concepts. Once those concepts are understood and tested, the PRL production code uses the same ideas with real atomic data and controlled physics.

Proposed repository structure: `education/` with `README.md`, `_quarto.yml`, `index.qmd`, `chapters/00_full_chain.qmd … 12_mini_lightcurve.qmd`, `notebooks/01_beer_lambert.ipynb … 11_nn_surrogate.ipynb`, `src/rtedu/{packets,slab,atom,macroatom,redistribution,sobolev,visualization}.py`, `figures/`, `videos/`, `data/`, `tests/`; and `paperB/` with `plan.md`, `prl_gate.md`. Quarto for the final HTML book (Markdown, LaTeX equations, Jupyter outputs, figures, videos and interactive plots in one navigable HTML document). The notebooks are for exploration; the `.qmd` chapters are the polished teaching material.

**Part I — learn radiative transfer from zero.**
Chapter 0, the entire kilonova chain: nuclear physics → abundances/heating → ejecta → atomic physics → RT → SED → light curves; then zoom into RT: propagate → interact → redistribute → propagate again. One static diagram and one short animation of a photon travelling through ejecta. Educational goal: know exactly where Paper IV, Paper B, and Peng's emulator live.
Chapter 1, the radiative-transfer equation: (1/c) ∂I_ν/∂t + n·∇I_ν = η_ν − χ_ν I_ν; strip it to pure absorption dI/ds = −χ I with I(s) = I_0 e^{−χ s} = I_0 e^{−τ}. Code demo: a slab with fixed opacity; plot I/I_0 versus distance, versus optical depth, and the transmission e^{−τ}. Validation: the numerical calculation reproduces the analytic solution essentially exactly. This establishes the meaning of τ.
Chapter 2, why Monte Carlo gives the same RT: derive P(τ_int > τ) = e^{−τ}, so with ξ ~ U(0,1), τ_MC = −ln ξ. Launch 10², 10³, 10⁴, 10⁵ packets through the slab; compare the Monte Carlo escape fraction with e^{−τ}. Visuals: histogram of sampled interaction optical depths with the analytic exponential overlaid; animated random packet trajectories; convergence error versus packet count. The first conceptual bridge: Monte Carlo packets are a numerical solver for RT.
Chapter 3, add frequency: three frequency groups B, V, IR with κ_B > κ_V > κ_IR, all other physics trivial. Launch equal numbers in each group; show blue photons interact more often, IR photons escape more easily, the emergent spectrum differs from the injected spectrum. Visuals: opacity versus wavelength, input SED, output SED, escape time by frequency. Key idea: changing photon frequency changes its future transport — that is why fluorescence matters.
Chapter 4, expanding ejecta and Sobolev resonance: homologous expansion v(r, t) = r/t; a packet with fixed lab-frame frequency experiences a changing comoving frequency ν_com ≃ ν_lab (1 − n·v/c), so it sweeps through atomic resonances as it travels. Introduce τ_S and P_int = 1 − e^{−τ_S}. Simulation: 10–20 discrete lines in frequency, one packet through an expanding shell. Important video: packet position on top, comoving frequency underneath, line frequencies as horizontal markers, a flash whenever a Sobolev resonance is crossed — you should literally see motion through space ⟹ motion through resonance frequency.

**Part II — what happens after a line interaction** (chapters 5–8: scalar ε; a toy atom built from first principles with Einstein A branching and animated cascades; photon-number versus indivisible energy packets, the chapter that explains why F50 forced the redo of the old Paper II/III results; the macroatom as a stochastic network over internal states with `activate/step/deactivate`, kept independent of the production code).

**Part III — derive the PRL idea yourself** (chapters 9–12: build R_ij = N(i→j)/Σ_k N(i→k) from the macroatom, throw the level network away and run with R alone, always the three-way comparison ε* versus R versus full macroatom; error versus model complexity for N_g = 1, 2, 4, 8, 16, 32 and contiguous grouping versus a constrained low-rank approximation; state dependence R(T_1), R(T_2), R(T_3) and interpolation; an optional neural surrogate (i, T, ρ, J_ν, …) → P(j) compared with ε, R_fixed, R_interpolated, R_NN, with whole physical states withheld — the question is not "neural networks are better" but "what complexity of effective model is actually required by the physics").

**Part IV — put it back into radiative transfer** (chapter 13: source → packet propagation → Sobolev interaction → redistribution → escape, three simulations ε versus R versus macroatom with identical seeds where practical, three animations, emergent SED, toy g, r, i, J, K light curves, colour evolution, error versus complexity; this becomes the miniature Paper B).

**Part V — bridge to the actual code**: a final chapter mapping every educational object to its production equivalent (slab opacity → ejecta opacity; three frequencies → real wavelength groups; discrete toy lines → GSI/atomic line lists; toy Sobolev line → real Sobolev resonance; five-level atom → La/Ce/Nd level network; toy cascade → macroatom; 3×3 R → production R_ij; three-state test → T, ρ, J_ν grid; toy mixture → P1 lanthanide mixture; toy light curve → xkn/Fontes transport benchmark).

**Every chapter has the same structure:** Question / Theory / Limiting cases / Numerical model / Demo / Visualization / Validation / What approximation did we introduce? / How does this map to revisit_sobolev? / What would fail if our assumption were wrong? / Exercises (change parameters and predict the answer before running).

**Visualization plan:** static figures (Beer–Lambert attenuation; sampled optical-depth distribution; frequency-dependent opacity; Sobolev resonance geometry; atomic level diagram; R_ij heat maps; scalar versus matrix redistribution; error-versus-complexity; state-transfer error; final SED/light curves) and short videos (photon random walk; expanding-medium frequency sweep; Sobolev line encounters; macroatom cascade; ε vs R vs macroatom side by side; R(T) changing with state; a full packet from creation to escape), 5–20 s, looping, generated reproducibly from Python (`scripts/make_video_*.py` → `education/videos/`).

**Reproducibility philosophy:** every displayed result comes from code; no hand-calculated result numbers in the HTML. `make education` runs tests, executes notebooks, regenerates figures, regenerates videos if requested, renders HTML; `make education-fast` skips the videos; `quarto render education/` builds the book.

**Tests are part of the teaching:** `test_beer_lambert`, `test_tau_sampling_is_exponential`, `test_branching_probabilities_sum_to_one`, `test_macroatom_energy_conservation`, `test_R_rows_sum_to_one`, `test_R_reproduces_macroatom_distribution`, `test_sobolev_interaction_probability` — each encodes a piece of physics you should understand.

**Two levels of educational code:** Level 1, transparent code inside notebooks (20–50 lines implementing the idea explicitly, optimised for understanding); Level 2, reusable code in `education/src/rtedu/`. Don't immediately reuse the complicated production classes; only at the final bridge chapter compare the toy implementation with the production implementation.

**Milestones:** E1 (chapters 0–4; exit: you can derive and implement how a packet finds its next line interaction); E2 (5–8; exit: you can explain precisely what ε approximates and why a macroatom is more physical); E3 (9–12; exit: you can derive R, identify the information discarded by it, and explain when an NN would actually be justified); E4 (13 + bridge; exit: one complete toy RT simulation produces an SED/light curve using all three: ε, R, macroatom). At that point start the PRL production campaign.

The final product: a small interactive textbook, "From Photon Random Walks to Effective Atomic Fluorescence: A Computational Introduction to Kilonova Radiative Transfer", with equations, derivations, executable examples, animations, exercises, and explicit links from every toy model to the actual `revisit_sobolev` research code.

## Answers to the planning questions (2026-09-15)

Deliver stage by stage, E1 first. Quarto via `pip install quarto-cli` into the venv; animations as GIF + HTML5 (no ffmpeg). Create the `paperB/` placeholder stubs now. (The chapter-12 surrogate implementation was not answered; the plan assumes a numpy MLP and asks before E3.)

## The eight amendments to the E1 plan (2026-09-15, verbatim)

This is a strong E1 specification. I would approve it as the implementation contract for the first PR, with a few amendments before coding so the educational physics stays clean and the build does not become unnecessarily fragile.

1. **Keep Chapter 3 conceptually below redistribution.** With B/V/IR opacities, first demonstrate pure extinction, I_g = I_{g,0} e^{−τ_g}, then optionally add **elastic scattering with frequency held fixed** to demonstrate that larger opacity produces more interactions/longer escape paths. Do not change a packet's frequency yet. The first actual frequency redistribution should happen in Chapter 5 with ε. This gives a very clean conceptual ladder: absorption/extinction → elastic scattering → Sobolev line interaction → redistribution in E2.

2. **Make the Sobolev chapter the real E1 capstone.** The exit criterion is exactly right: by the end of Chapter 4, the reader should be able to derive and implement ν_com(r), find the next resonance satisfying ν_com(r_res) = ν_ℓ, calculate its toy τ_S, and perform the interaction coin P_int = 1 − e^{−τ_S}. I would make the final Chapter-4 exercise implement that entire sequence from scratch in perhaps 30 lines, then compare it against `rtedu.sobolev.resonance_crossings`. Passing that exercise is a much more meaningful PI review than merely having all tests green.

3. **Be careful with `check_no_literals.py`.** The Paper-IV rule cannot be copied literally. An educational text legitimately contains fixed numbers such as 10², 10⁵, three frequency groups, five toy lines, 0.1c, example optical depths, etc. The checker should prohibit **generated numerical results typed into prose**, not pedagogical constants or mathematical examples. I would either have it inspect only specially marked "result" regions or provide an explicit exemption syntax. Otherwise this checker will fight the textbook.

4. **Make results generation transactional.** A single `results.json` is fine for the book, but having every notebook independently merge into it can leave stale keys after a notebook changes. A safer design is `data/generated/ch01.json … ch04.json` → `data/results.json` with a deterministic merge step. At minimum, `make education-notebooks` should delete/reinitialize `results.json` before executing Chapter 1 onward. Then the file genuinely represents the current build.

5. **Use an absolute `PYTHONPATH`.** Since notebooks are executed from subdirectories, make the Makefile export something equivalent to `PYTHONPATH=$(CURDIR)/src` rather than relying on relative `src`. This preserves the deliberate decision not to package/install `rtedu`.

6. **Treat animation determinism separately from physics determinism.** PNG/PDF, JSON and GIF hashes should be byte-stable. Matplotlib's HTML animation machinery can introduce generated IDs or other serialization differences even when the frames are identical. I would require physics/data/GIF reproducibility for E1 and verify HTML semantic/render equivalence, but not make byte-identical HTML animation files a scientific gate unless they happen to be deterministic. Also test the Quarto `include` route early: the generated animation file must be an embeddable HTML fragment rather than a complete nested HTML document.

7. **For Monte Carlo tests, test the statistical prediction rather than a magic tolerance.** For example, if direct counting tests p = 1 − e^{−τ}, use the expected binomial uncertainty σ_p = sqrt(p(1−p)/N) and require agreement within a preregistered multiple of σ_p. Likewise, a fixed-seed KS statistic is preferable to a test whose success depends only on a p-value. This makes the tests themselves educational.

8. **Keep one vocabulary distinction explicit throughout E1:** interaction probability and Sobolev escape probability are different quantities: P_int = 1 − e^{−τ}, β(τ) = (1 − e^{−τ})/τ. The asymptotes β → 1 (τ → 0), β → 1/τ (τ → ∞) apply to β, not to P_int. Given how important "β once" becomes later in the macroatom, establishing this terminology now will pay off in E2.

With those refinements, I would **freeze E1 scope exactly where you have it**. Do not sneak `ToyAtom`, ε, macroatom, or R_ij implementation into this PR. Chapter 0 can show where they will eventually appear in the chain, but E1 should stop precisely at: *Given a packet and a line forest, I understand how the RT solver finds the next Sobolev resonance and decides whether the packet interacts.* That separation will make E2 much more satisfying. E1 answers "How does the photon find an atom?" E2 will answer "Now that it hit the atom, what happens to its energy?"

I also agree with creating the `paperB/` stubs now but leaving them deliberately skeletal. We should not start filling the PRL claims from educational results. The book establishes understanding and validation infrastructure; the later production campaign establishes publication evidence.

So I would consider the proposed **`education: milestone E1` PR ready to implement**, with the eight adjustments above incorporated into the implementation notes.


## Continuation (2026-09-15)

After E1 was delivered as PR #5 the PI wrote: "continue to finish this plan". E2, E3 and E4 are therefore built in sequence as stacked PRs without waiting for the between-stage reviews; each stage keeps its own PR so the exit-criterion reviews can still happen per stage. For the chapter-12 surrogate the plan's default (a numpy MLP, no torch) is used.
