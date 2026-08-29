# Paper III plan — reduced redistribution closure (as received 2026-08-29)

Goal: an intermediate treatment between scalar TLA epsilon and explicit
branching/macroatom transport: a group-to-group redistribution operator
R_ij = P(re-emitted in group j | absorbed in group i). Question: can a small
matrix reproduce explicit lanthanide branching well enough for kilonova
light curves? Isolate redistribution first: Sobolev absorption + exact
branching vs Sobolev absorption + reduced R_ij -- identical opacity on both
sides, so any error is redistribution compression alone. Do NOT combine
with expansion opacity until the redistribution approximation is validated
(F24: its absorption placement fails catastrophically in dense forests).

Phases and gates:
- P0  Freeze the reference (La II, T_gas 3000 K, T_src 6000 K, reference
      n_ion, slow shell, Sobolev + A*beta branching, full launch range,
      energy-weighted, 3 seeds x 2e6). Gate 0: the wrapper reproduces the
      current branching results within MC noise, or stop.
- P1  Build the operator from the E6 machinery, made energy conserving:
      R^E_ij = E_{i->j}/E_i,abs with q_dep (and q_core where it applies);
      sum_j R^E_ij + q_i = 1 within noise. Reusable RedistributionKernel
      class (from_branching_mc, sample_output_group, validate_energy,
      save/load). One matrix at one physical state first.
- P2  New transport leg `sobolev_group_redistribution`: same Sobolev opacity
      as sobolev_branch; group of the absorbed packet; sample the row;
      re-emit in the sampled group (uniform or stored within-group PDF),
      same angular convention. NO inspection of atomic levels after
      absorption -- that is what makes it a real compressed closure.
- P3  Compression sweep N_g = 4, 8, 16, 32, 64, 128 vs full branching.
      Metrics: Delta L_bol; Delta F_b for UV/blue/optical/red/NIR/3800-3955;
      E5-style spectral chi^2; matrix row-L1, energy-weighted error, block
      flows (blue->blue/optical/red, optical->red); wall time, memory,
      table size, interactions/packet.
- Gate 1: error vs N_g and vs cost. Strong: N_g <= 32 gives |dF_b| < 5%
      all bands. Excellent: < 2% + good morphology. Weak: needs N_g ~
      100-500. Failure: no grouped matrix works (info lost in lambda_in ->
      lambda_out) -- stop; still a result.
- P4  Ions: Ce II (the strongest scalar/branching-closure failure), then
      Nd II (GSI; path to independent-data/NLTE tests). The matrices need
      not be similar -- both must be compressible.
- Gate 2: is compression generic? A: all three at 16-32 groups. B: La easy,
      Ce/Nd need many more (adaptive representation). C: dense ions cannot
      be compressed -- not a useful universal closure.
- P5  Thermodynamic robustness: T_src = 4000/5000/6000/8000 K with a FIXED
      6000 K kernel (transferability) and with recomputed kernels
      (representation vs state-transfer error); T_gas = 2500-5000 K with
      recomputed LTE state.
- P6  Epoch/density: t = 0.5/1/2/4 d; does R_ij(T,t) collapse onto
      R_ij(T, tau_scale)?
- P7  Smallest state space: which of T, n_e, tau_scale, ion fraction
      actually matter. Ideal: R_ij(T, tau_scale, ion).
- P8  Low-rank: SVD, C_k over k = 1..16, transport with R^(k) -- are there
      only a few macroscopic redistribution modes (stay local, blue->
      optical, optical->red)?
- P9  Mixtures: explicit La+Ce kernel vs opacity-weighted composition
      R_mix = sum_s w_{i,s} R_s. If the rule works, composition is cheap.
- P10 Atomic-data robustness (Nd): R^GSI vs R^independent -- is the coarse
      operator more robust than individual lines?
- P11 Only now, grouped opacity: Sobolev line-by-line vs line-binned vs
      expansion, each with the same R_ij -- separates opacity compression
      from redistribution compression. Target architecture: kappa_grouped
      + R_ij.
- P12 Light-curve toy model (1-D homologous, heating, LTE T iteration,
      16-64 groups, La/Ce/Nd): grey/TLA vs grouped+scalar vs grouped+R vs
      explicit reference; |dL_bol| < 2-5%, |dm_b| < 0.05-0.1 mag.
- P13 External validation vs TARDIS downbranch/macroatom on a controlled
      problem.
- P14 NLTE afterwards, populations decoupled behind the interface
      populations -> opacity + branching -> R_ij.

Required tests: energy conservation of rows; identity kernel == coherent
group scattering; blackbody stationarity (with the paired
emissivity/absorption closure); fine-group convergence to branching; rebin
invariance; mixture limit X_Ce -> 0; zero-opacity limit.

Stop/go: stop after P3 if 32-64 groups cannot reproduce La; stop after P4
if Ce/Nd need line-level resolution; methods paper if 16-64 groups do
several ions at <~5%; light curves if the kernel is accurate across T and
epoch with low-dimensional state.

Immediate experiments R1-R5: build N_g = 8/16/32/64 kernels from the La II
branching events (R1); add the transport leg (R2); compare each against
sobolev_branch (R3); repeat winners on Ce II (R4); add Nd II (R5 --
BLOCKED: no Nd II GSI files in data/ yet).

The expected hierarchy: epsilon -> R^8 -> R^16 -> R^32 -> explicit
branching -> macroatom, with measured accuracy vs cost vs regime.
