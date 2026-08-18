# P2-0A — Is direct line fluorescence usable in the public SEDONA?

**Verdict: NO — absent, not deprecated.** The public `dnkasen/pubsed` release
contains no line-branching machinery at all. Building Paper II's reference
treatment on it would mean writing the physics, not enabling it.

Audited 2026-08-18 against the clone at `~/personal/pubsed`.

## What was searched

```
grep -rn "fluor|branch|downbranch|macroatom|two_level|TLA" --include=*.cpp --include=*.h src/
```

The **only** hit outside `sandbox/` is a comment in `transport/radioactive.h`
about a nuclear decay branch. There is no fluorescence path, no macroatom, no
downbranching, and no per-line re-emission bookkeeping.

## What SEDONA actually does at a line interaction

The chain is short and worth stating precisely, because it determines what
Paper II can and cannot ask of this code.

1. **Opacity split** — `opacity/GasState_opacities.cpp:438`

   ```cpp
   double this_eps = epsilon_;               // opacity_epsilon, a scalar
   atoms[i].line_expansion_opacity(atom_opac, time_);
   opac[j]  += atom_opac[j];                 // total extinction
   aopac[j] += atom_opac[j]*this_eps;        // absorptive part
   ```

   `epsilon` is a **single global number** (optionally zeroed per element via
   `opacity_atom_zero_epsilon`). It is not per line, not per level, and
   carries no information about which transition absorbed the photon.

2. **Interaction** — `transport/scatter.cpp:12`, `do_scatter(p, eps)`

   ```cpp
   if (rangen.uniform() > eps)  isotropic_scatter(p, 0);   // (1-eps)
   else {
     if ((z2 > zone->eps_imc) || radiative_eq) isotropic_scatter(p, 1);
     else fate = absorbed;
   }
   ```

3. **Frequency handling** — `transport/scatter.cpp:303`

   ```cpp
   if (redist) sample_photon_frequency(p);   // redist = the 2nd argument
   ```

So the two channels are:

| channel | probability | frequency outcome |
|---|---|---|
| resonant scattering | `1 − ε` | **coherent in the comoving frame** — direction changes, frequency does not |
| effective scattering | `ε` | redrawn from the **zone's total emissivity**, i.e. the thermal pool |

Neither channel knows which line absorbed the photon. There is no
`λ_in → λ_out` kernel tied to an upper level's radiative branching — exactly
the quantity Paper II exists to measure.

## Why this is a "write it" not a "switch it on"

Implementing direct fluorescence here would require, at minimum:

- upper-level identity carried through the interaction (the expansion-opacity
  path bins lines and discards which line was hit — this is the same
  information loss that produces the F4 per-crossing cap);
- a per-upper-level branching table `A_uj / Σ_j A_uj` built from the atomic
  data;
- a new interaction branch that samples a downward transition and re-emits at
  its wavelength;
- bookkeeping to keep the resolved bound-bound path and the new branching path
  consistent.

That is core packet-interaction work, not a patch. Estimated effort is well
beyond the one-to-two-week Phase 0 budget, and it would put a bespoke
modification of a production code at the centre of Paper II's reference leg,
which is exactly the kind of confounding the project has avoided so far.

## Recommendation

Per the plan's own stop rule ("if it requires rewriting core packet-interaction
machinery: consider TARDIS or a small custom controlled implementation"):

1. **TARDIS as the reference.** It has macroatom, downbranch and scatter modes
   as first-class options, Python-facing, and Carsus already ingests GSI
   lanthanide data. It gives the `λ_in → λ_out` kernel directly.
2. **Keep SEDONA as the expansion-opacity leg.** Its global-`ε` treatment is
   precisely the approximation Paper II is testing, and Paper I has already
   characterised it. This is a feature: the comparison becomes
   *TARDIS macroatom* (reference) vs *SEDONA expansion opacity + global ε*
   (approximation), with the caveat that it is now cross-code — so the three
   conventions Paper I identified (thermal emission F8, spectral normalization,
   transport treatment F11) must all be matched explicitly.
3. **Or a small custom branching MC** for the controlled phase only, where the
   three-level and La II mini-atom tests need no production code at all. Given
   how much of Paper I's value came from a ~200-line solver that could be
   trusted line by line, this is attractive for P2-0B and Phase 1.

## Consequence for P2-0B

The three-level synthetic-atom test cannot be run in SEDONA as planned — there
is no branching to measure. It should be run in TARDIS, or in a small custom
branching Monte Carlo, where the analytic ratios `A₃₁/(A₃₁+A₃₂)` and
`A₃₂/(A₃₁+A₃₂)` are recoverable by construction.
