# Response to the referee

We thank the referee for a careful and constructive report. Its central
point — that the manuscript labelled a difference between a deterministic
finite-profile attenuation and a statistical closure as proof that the closure
is intrinsically erroneous — is correct, and acting on it has made the paper
both more defensible and more interesting. We have adopted the reframing the
referee sketches as the paper's thesis, added the targeted calculations asked
for, redefined the Sobolev differential against a deterministic reference,
re-examined the single-line benchmark and the Monte Carlo uncertainties, and
tightened the presentation. Below we answer each point in turn; section and
figure numbers refer to the revised manuscript. Where we disagree in part we
say so.

Two things changed in the numbers, and we flag them up front because the
referee should not have to find them. First, the breadth-sweep Sobolev
column (median +3.6%, max +11.3%) in the submitted draft was computed with a
normalization that the rest of the paper had already abandoned; recomputed
through the one shared routine it is median +0.34%, max +7.8%, and against
the deterministic reference median +0.26%, max +7.6% — the boundary law, not
a floor — while the expansion column is unchanged (median +5.1%, max +62%). Second, the single-line benchmark's deterministic value
(0.1372) was never a discrepancy: it is the frozen-snapshot law the Monte
Carlo code solves, to 0.1%, and the Monte Carlo value (0.1420) was dominated
by a 2×10⁶-packet fluctuation: at the converged grid with 3.2×10⁷ packets
and five seeds it is 0.1383 ± 0.0001, and the +0.9% left is the continuum
offset the zero-opacity control measures (1.008); corrected, 0.1372.

---

## Major Comment 1 — the interpretation of expansion opacity

> The problem is what the manuscript concludes that difference means … I do
> not think the numerical result by itself licenses the current statement that
> the expansion-opacity formalism is "wrong".

We agree, and this is now the paper's thesis rather than an objection to it.

**What we did.** Section 2.5 and Appendix A.4 are rewritten around the
distinction the referee draws. Per ray, the Sobolev leg exponentiates
S = Σ τ_k and the expansion leg exponentiates E = Σ (1 − e^{−τ_k}). E is
exactly what α_exp integrates to — the expected number of line interactions
per crossing for a photon that is counted but not removed, which is the
mean-free-path statistic Karp et al. (1977) built the coefficient from — and
the closure reproduces it identically at every bin width. Transmission is a
different statistic: e^{−S} is a product of independent Bernoulli survivals,
e^{−E} the survival of a Poisson process with the same mean. They agree when
every τ_k ≪ 1 (the many-weak-lines limit in which the formalism was derived
and is exact) and separate, per resonance, by the saturation deficit
D = Σ [τ_k − (1 − e^{−τ_k})] once lines saturate. For a band average the
statement is an identity, F_exp/F_Sob = ⟨e^{D}⟩ weighted by transmission
(Eq. A-identity), pinned in the test suite to 10⁻¹³.

**The diagnostic the referee asked for** is Fig. [count] (new): along one
radial ray through the La II forest (46 resonances), the cumulative expected
count reaches the same value in the Sobolev and expansion treatments at the
far edge of the shell — identical to 10⁻¹⁵ — while the survivals are 0.051 and
0.130; the first-interaction distribution is discrete masses at the resonances
against a smooth α_exp e^{−E}, and integrated, 0.949 against 0.870, so the
closure under-counts absorptions. A 2×10⁶-photon Monte Carlo on the same ray
reproduces every panel to its sampling error. The integrated version is the
same figure's last panel: Δ_exp against τ_max with 1 − ⟨E⟩/⟨S⟩ overlaid, the
gap closing to +0.2% at τ_max = 0.1.

**What it changes.** The words "wrong", "intrinsic error", "fail", "at fault"
are gone from the Abstract, Introduction, Sections 2.5, 5, the Conclusions and
Appendix A.4. The statement is now that the closure is exact in the statistic
it was derived for and departs from deterministic transmission once lines
saturate; we quantify that departure (+38% on the La II forest at thermal
widths, the pure Poisson-versus-Bernoulli number; +45% in the Monte Carlo
implementation at 100 km/s, which additionally carries a ~2% bin systematic
in its expansion leg that we now report rather than fold in) and identify the
ensemble statistic that controls it. Findings F4, F5, F7 and F13 of the
submitted draft are corollaries of this one mechanism and are now collected
under it (R2 in the revised findings table). The connection to Morag (2026),
who reaches the same conclusion from the emissivity side, is made explicit in
Section 2.5 and Appendix B.4.

We did not find, and do not claim, that Karp et al. or Eastman & Pinto
intended the coefficient to reproduce deterministic attenuation; Section 2.5
now states the derivation as what it is.

## Major Comment 2 — novelty relative to Fontes et al. (2020)

> F4 — the existence of the 1 − e^{−τ} cap — cannot reasonably be presented as
> a new discovery.

Conceded. Fontes et al. (2020) write τ_exp = Σ(1 − e^{−τ_i}) against the
line-binned Σ τ_i, state that the former limits a single line's contribution
to one, and call the continuous Monte Carlo Sobolev treatment the more
accurate one. The Introduction (first paragraph of the positioning) and
Appendix B.1 now say so in those terms, F4 is no longer a finding, and the
novelty is stated as the referee suggests: the finite-profile resolved leg
inserted ahead of per-line Sobolev transfer; the controlled map over width,
strength, bin width, line density and calibrated lanthanide lists; and the
identification of the controlling ensemble statistic.

## Major Comment 3 — the scope of the Sobolev claim

> neither of the two Sobolev-intrinsic assumptions is actually stressed …
> There is also an apparent tension between the abstract and the breadth
> sweep.

Both points taken. (i) The Limitations paragraph now states the narrow claim
in the referee's own words: for fixed populations, pure absorption, homologous
expansion and narrow intrinsic profiles, replacing each finite profile by a
localized resonance introduces only the boundary error of Appendix A.6;
locality and isolation are deliberately protected and what that licenses is
said. "Universal" is gone; the breadth result is "robust across the axes
tested". (ii) The headline Sobolev number now carries its width and strength
qualifiers everywhere it appears, including the Abstract: ≲1% at v_D ≤ 30
km/s for τ_max ≤ 5 (2.5% at τ_max = 50), ≲0.3% at v_D ≤ 10 km/s, vanishing at
thermal widths — all against the deterministic reference of Comment 5. The
tension with the breadth sweep was real and had a second cause beyond
wording: the sweep's Sobolev column had been computed with a normalization
the rest of the paper had abandoned (raw luminosity in the margin rather than
the continuum ratio). Recomputed through the shared routine it is median +0.34% and max +7.8%
(against the deterministic reference, +0.26% and +7.6%; +2.6% median where
τ_max > 3 at v_D = 100 km/s), consistent with the boundary law, and the
deterministic reference agrees with SEDONA resolved to −0.04% median, 0.31%
spread, across all 36 conditions.

## Major Comment 4 — emergent-spectrum claims

> I do not see a rigorous basis for the claimed upper bound. Nor do I think
> the sign … is guaranteed.

We agree and have withdrawn both the "upper bound" and the "sign survives"
arguments. Rather than leave the point untested we ran the one
energy-conserving experiment the referee allows for: the La II forest pair
with SEDONA's radiative-equilibrium switch on, resolved against expansion,
seed-matched, in two variants — a single iteration that re-emits at the input
temperature (isolating redistribution from population feedback) and an
iterated run in which the gas temperature converges (to different values in
the two modes, so that differential includes feedback). Three seed-matched
pairs per variant, normalized on a blue margin because re-emission redshifted
by up to 3000 km/s contaminates the red one. In the emergent 3800–3955 Å band
the differential is +7.7 ± 0.02% after one iteration and +5.0 ± 0.7% at
convergence (|δT/T| < 0.7% over the last three iterations; the median gas
temperature rises from 3000 to ~7700 K in both modes), against +44% in pure
absorption on the same normalization. The band fills from 0.34 to 0.84–0.86 in
both treatments; the removed flux reappears immediately redward (3955–3995 Å
at 1.08–1.10 of the continuum); a τ_max = 0.05 control returns both modes to
the continuum within 0.6%. So most of the closure's departure in attenuation
does not survive into the emergent band once absorbed energy is re-emitted —
the differential shrinks by a factor of six to nine — while its sign is
preserved. It is reported in Section 5 as one check on one configuration, explicitly not as the emergent-spectrum error budget, and we
make no claim about the sign or size of the error in a mixed,
multidimensional kilonova spectrum.

## Major Comment 5 — Δ_Sob against the deterministic reference

> The paper already contains the more natural reference calculation.

Done, and it exposed something worth reporting. Δ_Sob is now formed against
the deterministic finite-profile solution on identical rays, with the same
source convention (T_shell → 0), populations and transport treatment; SEDONA's
resolved run, at its seed mean, is the independent validation of that
reference (it agrees to 0.5% on the headline forest, to ±0.3% across the
strength–width grid with a maximum of 0.7%, and to −0.04% median, 0.31%
spread, across the 36 breadth conditions).

Two technical points made this possible. First, the resolved leg for uniform
populations, Gaussian profiles and pure absorption has a closed form — the
erf bracket of Appendix A.6 applied per line and per ray — whose cost is
independent of the Doppler width, so the 1 and 3 km/s frontier points, which
cost the brute-force solver hours, are free; the closed form agrees with the
brute-force solver on identical rays to 8×10⁻⁵. Second, the two legs must
share their impact-parameter quadrature: attenuation is a step function of p,
and with the analytic leg on 200 midpoint rays against the solver's
n_impact/2, Δ_Sob drifted by a percentage point and could change sign with
ray count (we show the table in Section 3). With matched rays it is flat to
0.014 points from 100 to 1600 rays.

The headline becomes +3.1% at v_D = 100 km/s against the reference (+2.6%
against SEDONA's ten-seed mean), rather than +2.1% against a single SEDONA
realization; at v_D ≤ 10 km/s it is ≤ 0.3%, and the whole width axis follows
the boundary law. We consider the larger, noise-free number the right one to
quote.

## Major Comment 6 — the single-line benchmark

> For a paper whose scientific objective is to distinguish effects at the
> <1–2% level, this should not simply be labelled "harness validated".

Agreed. Two things were wrong with how the benchmark was presented. The
target for a code that iterates on a fixed snapshot is not e^{−τ_S} = 0.1353
but the frozen-snapshot law of Appendix A.7, e^{−τ_S(1−β)/γ}, which over the
trough window is 0.1371; the deterministic solver in frozen mode reproduces
it to 0.1% under every variation of ray count, shell emission and profile
truncation, so the 1.4% between it and 0.1353 is the transport treatment, not
an error. SEDONA's production run sat 3.5% above that target, on a transport
grid with under two bins per Doppler width against the eight used in every
forest run and with 2×10⁶ packets.

We ran the convergence ladder the referee asks for — packet number, transport
grid, spectrum grid, zone count, line-to-cell phase, with fixed seeds and a
zero-opacity control (Section 3.2, Table 3). On the production grid alone,
3.2×10⁷ packets give 0.1386 ± 0.0002, so most of the gap was sampling; at the
converged grid (eight bins per Doppler width) five seeds give 0.1383 ± 0.0001,
a residual of +0.9% above the frozen target that does not move with packet
count, grid or zoning. The zero-opacity control identifies it: with the line
removed the same trough window reads 1.0079 ± 0.002 against the red-margin
normalization — SEDONA's continuum is 0.8% high there relative to the
analytic Planck shape, because its core packets are emitted in the comoving
frame of the core surface and boosted to the lab frame — and dividing it out
gives 0.1372, the frozen target to 0.1%. The expansion mode at the converged
grid is 0.4438 ± 0.0003 (0.440 corrected), 4.5% above the Poisson value, the
expansion-leg bin systematic of the bin-width section on a single line. The
table and the Fig. 4 caption now quote the converged values, and the phrase
"harness validated" is gone.

## Major Comment 7 — Monte Carlo uncertainties

> Sharing a deterministic model does not create Monte Carlo covariance.

Correct as the manuscript stood: the production legs were seeded from the
clock, minutes apart, and were independent. With the same fixed seed in both
legs the core-emission stream is identical and the streams diverge only at
the first interaction; measured on seed-matched pairs the resolved and
expansion band fluxes correlate at +0.95 (headline, ten seeds; +0.92 to +0.999
across the grid) and the paired Δ_exp scatters by 0.12% (0.04% on its mean)
against 0.36% from quadrature. Every Δ_exp in the revised paper is a
seed-matched pair — ten seeds for the headline forest, three at every point of
the strength–width grid and the 3 km/s frontier — and every headline band flux
a seed mean with its standard error (resolved 0.3422 ± 0.0003, expansion
0.4958 ± 0.0004); the phrase "partially correlated" is gone and the measured
correlation is quoted. The single-realization seed scatter sits within 1.5×
the Poisson expectation from the packet count at every grid point. The 2σ
offset of the single-realization headline no longer enters anything.

## Major Comment 8 — τ_max as the controlling variable

> A more physically motivated predictor would involve the full saturation
> deficit … appropriately weighted.

Yes — and the appropriate weighting turns out to be exact. The identity of
Comment 1 says F_exp/F_Sob = ⟨e^{D}⟩ with weights ∝ e^{−S}, the transmission-
weighted mean of e^{D}; the unweighted Σ[τ − (1 − e^{−τ})] the referee
suggests is a Jensen lower bound on its logarithm and fails once rays
saturate. Fig. [predictor] (new) plots ln(1 + Δ_exp) for the 36 breadth
conditions and the 12-point grid against the exact predictor (1:1 by
construction for the analytic pair, with SEDONA's points sitting above it by
its expansion-leg bin systematic), against ⟨D⟩_w, and against τ_max. Section
4.5 now says that τ_max is an empirical proxy within a fixed (v_D, geometry)
slice — good because the strongest few lines dominate D on the rays that
still transmit — and not the fundamental variable; the sentence "the error can
be estimated from the strongest line in a band" is removed.

## Major Comment 9 — expansion velocities

> Either one realistic-velocity/time-dependent check should be included, or
> the scope should be stated more narrowly.

We have chosen the narrower scope, stated in the referee's words ("kilonova
line lists and thermodynamic conditions in a controlled homologous test
flow") in Section 2 and the Limitations, with one paragraph added to
Appendix A.7 recording what was learned in extending the machinery to
0.1–0.3c: the controlling term there is not the 1/γ worldline correction
(4.6% at β = 0.3) but light-travel dilution of the medium, (1 − β)² (51%),
confirmed against direct integration; and under worldline transport the
resonance locus is linear in z and τ/τ_S = 1/γ for every impact parameter, so
the Jeffery common-point/common-direction caveat the submitted draft carried
does not arise in the physical treatment. We also corrected the sentence that
SEDONA's steady mode "cannot supply" the time-dependent side: the mode cannot;
the code's time-dependent mode can, and is the route for the follow-up. No
measurement at those velocities is claimed in this paper.

---

## Minor and presentation comments

- **Findings table.** Collapsed from thirteen entries to four results (the
  separation; the mechanism and its controlling statistic; the geometric
  character of the Sobolev error; the three conventions), each listing the
  subsidiary findings it absorbs.
- **Draft-history narrative.** All eleven passages ("previous draft",
  "retracted", "early pass", "we briefly did", "both were wrong", …) are
  removed; the lessons survive as stated validation tests and controls.
- **Primer.** Kept, since the paper is written to be read without a
  radiative-transfer background, but the conversational phrases the referee
  quotes and several others are replaced and the duplication with Appendix A
  reduced.
- **"Transport codes need an opacity per frequency bin."** Now "codes that
  transport on a frequency grid", with line-by-line Monte Carlo schemes named
  as the exception.
- **"Means expansion opacity."** Softened throughout to "a code that works
  from binned opacities … is applying both steps".
- **Stimulated emission.** The factor (1 − n_u g_l/n_l g_u), equal to
  1 − e^{−hν/kT} for the LTE populations used, is now stated in Section 2.4
  and Appendix A.2 with its magnitude (3×10⁻⁶ at 3800 Å and 3000 K,
  5×10⁻³ at 9100 Å); SEDONA applies it, and the analytic legs now fold it into
  the populations so the two agree.
- **Morag.** Cited as the published MNRAS 549(3) 2026 article.

## Where we disagree in part

Nowhere on substance. On two points of emphasis: we kept the primer (the
paper's stated audience includes readers without a radiative-transfer
background, and its length is now reduced rather than removed), and we kept
the cost subsection as a short motivating note rather than a finding.
