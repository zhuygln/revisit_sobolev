# Understanding the Thermalisation Parameter \(\epsilon\) in Line Radiative Transfer

## From the classical two-level atom to expansion opacity, fluorescence, macro-atoms, and the planned lanthanide experiments

**Purpose.** This document is a focused study guide for understanding what the line thermalisation parameter \(\epsilon\) means, why values such as \(0\), \(0.01\), \(0.3\), and \(1\) appear in supernova and kilonova radiative-transfer calculations, and how that parameter differs from explicit multilevel fluorescence through a real atomic branching network.

It is written to support the current Sobolev / expansion-opacity / lanthanide-fluorescence research program. The central experimental question is:

\[
\boxed{
\text{Can a scalar two-level-atom thermalisation closure that works reasonably well for}
\atop
\text{iron-peak SN-Ia fluorescence also reproduce lanthanide fluorescence in kilonova ejecta?}
}
\]

The short answer is that \(\epsilon\) has **two related but importantly different meanings** in the literature:

1. in the classical two-level atom, \(\epsilon\) is closely related to a **true destruction or thermalisation probability**, usually set by collisions;
2. in many supernova/kilonova LTE transport schemes, a much larger effective \(\epsilon_0\) is used as a **phenomenological wavelength-redistribution parameter** to mimic fluorescence that the simplified line treatment does not explicitly retain.

The distinction is essential for interpreting an \(\epsilon\)-sweep.

---

# 1. The conceptual map

The progression from textbook line transfer to the current problem can be organized as

\[
\boxed{
\text{two-level atom}
\rightarrow
\text{Sobolev escape}
\rightarrow
\text{TLA redistribution}
\rightarrow
\text{multilevel fluorescence}
\rightarrow
\text{macro-atom / explicit branching}
}
\]

A useful schematic is:

```text
CLASSICAL TWO-LEVEL ATOM

       u
       |
       |  radiative excitation / de-excitation
       |
       l

     collisions can destroy the line photon
     and couple its energy to the thermal pool

             |
             v

SOURCE-FUNCTION CLOSURE

S_nu = (1 - epsilon) J_nu + epsilon B_nu

epsilon ~ true destruction probability
in the simplest textbook interpretation

             |
             v

SUPERNOVA / KILONOVA TLA APPROXIMATION

line interaction
       |
       +-- probability 1-epsilon --> scatter
       |
       +-- probability epsilon --> "absorb" and
                                   thermally redistribute

Here epsilon may be an EFFECTIVE parameter,
not the microscopic collisional destruction probability.

             |
             v

REAL MULTILEVEL ATOM

                   u
              /    |    \
            k1     k2    k3
             \     |     /
       A_uk beta_uk competition
                |
                v
       fluorescence / cascades

The output wavelength depends on the actual level network.
```

The research question is therefore not simply “what is the correct \(\epsilon\)?” It is:

> **Can the enormously structured multilevel redistribution operator of a lanthanide ion be compressed into one scalar \(\epsilon\) without losing important spectral information?**

---

# 2. Classical radiative transfer: where \(\epsilon\) comes from

## 2.1 Transfer equation and source function

For specific intensity \(I_\nu\),

\[
\frac{dI_\nu}{ds}
=
-\chi_\nu I_\nu + \eta_\nu ,
\]

where

- \(\chi_\nu\) is the extinction coefficient;
- \(\eta_\nu\) is the emissivity.

The source function is

\[
S_\nu
\equiv
\frac{\eta_\nu}{\chi_\nu}.
\]

For a simple line that both scatters photons and couples some fraction of their energy to matter, a common two-level-atom form is

\[
\boxed{
S_\nu
=
(1-\epsilon)J_\nu
+
\epsilon B_\nu(T)
}
\]

where

- \(J_\nu\) is the angle-averaged radiation field;
- \(B_\nu(T)\) is the Planck function;
- \(\epsilon\) measures the strength of thermal coupling.

This equation is one of the most useful conceptual equations in line-transfer theory.

---

## 2.2 What \(\epsilon=0\) means

If

\[
\epsilon=0,
\]

then

\[
S_\nu=J_\nu.
\]

The line behaves like a purely scattering line.

A photon interacts with the line, but the interaction does not drive the radiation toward the local thermal distribution.

Schematically:

```text
incoming photon
      |
      v
   line interaction
      |
      v
  scattering-like
      |
      v
same local radiation pool
```

No true thermal destruction occurs in the idealized limit.

---

## 2.3 What \(\epsilon=1\) means

If

\[
\epsilon=1,
\]

then

\[
S_\nu=B_\nu(T).
\]

The line is treated as fully thermalized.

The interaction is effectively interpreted as:

```text
incoming photon
      |
      v
    absorbed
      |
      v
identity forgotten
      |
      v
thermal pool at T
      |
      v
re-emission drawn from thermal emissivity
```

In a simple two-level textbook problem this is the strong-collisional-coupling limit.

In supernova/kilonova expansion-opacity applications, however, \(\epsilon=1\) often plays a more phenomenological role: it can be used to generate strong wavelength redistribution even when the **true microscopic destruction probability is nowhere near unity**.

That is the key subtlety of this entire subject.

---

# 3. Microscopic or “true” thermalisation probability

## 3.1 Radiative versus collisional de-excitation

Consider an excited upper level \(u\).

It can de-excite through:

- spontaneous radiative transitions with Einstein coefficients \(A_{uk}\);
- collisions with thermal electrons, characterized by rates \(n_e C_{uk}\).

A simplified two-level estimate of the destruction probability is often written schematically as

\[
\epsilon_{\rm true}
\sim
\frac{C_{ul}}
{A_{ul}+C_{ul}},
\]

where \(C_{ul}\) here denotes the relevant collisional rate.

For a multilevel atom in an expanding medium, the corresponding expression is more naturally written using all accessible downward channels and Sobolev escape probabilities:

\[
\boxed{
p_{bb}^{\rm abs,true}
=
\frac{
\sum_{\rm down} q_{\rm down}n_e
}{
\sum_{\rm down}
\left(q_{\rm down}n_e + A\beta_s\right)
}
}
\]

in the notation summarized by Jerkstrand (2025).

The important physical point is simple:

\[
\boxed{
\text{true thermalisation competes with radiative escape.}
}
\]

---

## 3.2 Why the true probability is tiny in expanding ejecta

Allowed radiative transitions can have

\[
A\gtrsim10^6\ {\rm s^{-1}},
\]

whereas collisional rates in dilute supernova/kilonova ejecta are often much slower.

Consequently, the modern review by Jerkstrand summarizes typical true line-absorption/thermalisation probabilities as roughly

\[
\boxed{
p_{bb}^{\rm abs,true}
\sim10^{-6}-10^{-4}
}
\]

when allowed radiative de-excitation channels exist.

This is dramatically smaller than

\[
\epsilon=0.3
\]

or

\[
\epsilon=1.
\]

Therefore, when an LTE supernova code adopts a constant

\[
\epsilon_0\sim0.3-1,
\]

one should **not** interpret that as saying that 30–100% of real line interactions are collisionally thermalized.

---

# 4. Why supernova codes nevertheless use a large \(\epsilon_0\)

## 4.1 The missing process: fluorescence

Suppose a photon excites an atom from level \(l\) to an upper level \(u\).

A real multilevel atom may then de-excite not back to \(l\), but to some other lower level \(k\):

\[
l\rightarrow u\rightarrow k,
\qquad
k\neq l.
\]

The output photon then has a different wavelength.

This is fluorescence.

In a complex ion, repeated absorption and radiative de-excitation can produce a cascade:

\[
\lambda_{\rm UV}
\rightarrow
\lambda_1
\rightarrow
\lambda_2
\rightarrow
\cdots
\rightarrow
\lambda_{\rm optical/NIR}.
\]

This process redistributes radiation strongly in wavelength **without requiring collisional thermalisation**.

---

## 4.2 Why a large effective \(\epsilon_0\) can mimic fluorescence

If the transport method discards the identity of the individual atomic line after forming a binned opacity, the detailed branching network is no longer naturally available.

A simple workaround is:

\[
\boxed{
p_{bb}^{\rm abs}=\epsilon_0
}
\]

with

- probability \(1-\epsilon_0\): scatter;
- probability \(\epsilon_0\): absorb and thermally re-emit.

The Jerkstrand review explicitly emphasizes that many LTE codes use this kind of constant effective parameter, often with a value close to unity, because otherwise the calculation lacks enough wavelength redistribution.

Thus:

\[
\boxed{
\epsilon_{\rm true}
\neq
\epsilon_0^{\rm phenomenological}
}
\]

in general.

A useful conceptual distinction is

\[
\begin{array}{lll}
\epsilon_{\rm true}
&
\rightarrow
&
\text{actual matter--radiation thermal coupling},
\\[2mm]
\epsilon_0
&
\rightarrow
&
\text{effective redistribution knob used to mimic omitted fluorescence}.
\end{array}
\]

---

# 5. The two-level-atom expansion-opacity formulation

Kasen, Thomas & Nugent (2006) provide one of the clearest formulations used in supernova transport.

For binned expansion opacity,

\[
\alpha_{\rm exp}(\lambda_c)
=
\frac{1}{ct_{\rm exp}}
\sum_i
\frac{\lambda_i}{\Delta\lambda_c}
\left(1-e^{-\tau_i}\right).
\]

The associated two-level-atom source function is

\[
\boxed{
S_\lambda
=
(1-\epsilon)\bar J_\lambda
+
\epsilon B_\lambda(T).
}
\]

Kasen et al. describe \(\epsilon\) as a probability associated with absorption/redistribution and note that, in principle, it could be unique for every line.

For computational simplicity, one often assigns one common value to every line.

This makes the entire detailed atomic network collapse to one scalar:

\[
\boxed{
\{\text{millions of lines and branching paths}\}
\rightarrow
\epsilon.
}
\]

That is precisely the approximation the current lanthanide study is testing.

---

# 6. Sobolev escape probability and why \(\beta\) matters

## 6.1 A photon emitted in a resonance is not automatically free

For a Sobolev optical depth \(\tau\), the escape probability is

\[
\boxed{
\beta_{\rm esc}
=
\frac{1-e^{-\tau}}{\tau}.
}
\]

Here \(\beta_{\rm esc}\) should not be confused with the relativistic speed parameter \(v/c\).

Two limits are important.

For an optically thin transition,

\[
\tau\ll1
\quad\Rightarrow\quad
\beta_{\rm esc}\approx1.
\]

For a very optically thick transition,

\[
\tau\gg1
\quad\Rightarrow\quad
\beta_{\rm esc}\approx\frac{1}{\tau}.
\]

Thus a photon emitted in a very thick resonance line is likely to be reabsorbed before it escapes the resonance region.

---

## 6.2 Effective radiative branching

Suppose upper level \(u\) has several downward transitions.

The raw spontaneous rates are

\[
A_{u1},A_{u2},A_{u3},\ldots
\]

but photons from some transitions may be strongly trapped.

The competition among **escaping** radiative channels is therefore approximately governed by

\[
\boxed{
A_{uk}\beta_{uk}.
}
\]

This is why explicit branching cannot generally be determined from \(A\)-values alone.

A very optically strong resonance transition may have a huge \(A\), but a small escape probability. Another weaker branch can then become the route through which energy actually leaves the resonance network.

This provides an important physical connection between the two research papers:

```text
Paper I:
optical-depth saturation changes the interaction/attenuation statistics.

Paper II:
the same optical-depth saturation changes which radiative
branch actually allows the energy to escape.
```

---

# 7. Direct fluorescence probability

Kasen et al. derive an approximate line-specific fluorescence probability

\[
\boxed{
p_{\rm fluor}
=
\frac{
\sum_{k\neq l}\beta_{uk}A_{uk}
}{
n_e\sum_k C_{uk}
+
\sum_k\beta_{uk}A_{uk}
}
}
\]

for excitation

\[
l\rightarrow u.
\]

This is the probability that the subsequent successful radiative escape occurs through a transition to a lower level other than the original \(l\).

If collisions are negligible,

\[
n_e C_{uk}
\ll
A_{uk}\beta_{uk},
\]

then approximately

\[
p_{\rm fluor}
\simeq
\frac{
\sum_{k\neq l}\beta_{uk}A_{uk}
}{
\sum_k\beta_{uk}A_{uk}
}.
\]

Kasen et al. explicitly suggested that one could eliminate the free scalar \(\epsilon\) by taking a line-dependent

\[
\epsilon_i=p_{{\rm fluor},i}.
\]

This is directly relevant to the planned “physics-informed closure” experiment.

---

# 8. The classic SN-Ia experiment: Kasen et al. (2006)

This is the most important single paper for interpreting the planned \(\epsilon\)-sweep.

## 8.1 Setup

Kasen et al. compared:

\[
\boxed{
\text{direct treatment of line fluorescence}
}
\]

with

\[
\boxed{
\text{expansion opacity + scalar TLA redistribution}
}
\]

for the W7 Type-Ia supernova model using nearly \(5\times10^5\) lines.

They tested

\[
\epsilon=1.0,\quad0.3,\quad0.01.
\]

---

## 8.2 Result for \(\epsilon=1\)

They found that high redistribution probabilities reproduced the direct-fluorescence spectrum reasonably well overall, but

\[
\epsilon=1
\]

redistributed somewhat too much flux toward red wavelengths.

Conceptually:

```text
direct iron-peak fluorescence
          |
          |   compared with
          v
epsilon = 1 TLA
          |
          v
too much red redistribution
```

---

## 8.3 Result for \(\epsilon=0.3\)

For their SN-Ia model,

\[
\epsilon=0.3
\]

reproduced the colors more accurately than \(\epsilon=1\).

This is the historical reason the current lanthanide experiment cannot stop after demonstrating that \(\epsilon=1\) fails.

A knowledgeable referee can immediately ask:

> “Does an intermediate \(\epsilon\) recover the direct fluorescence, as it did for iron-peak SN-Ia ejecta?”

That question must be answered.

---

## 8.4 Result for \(\epsilon=0.01\)

The nearly scattering case

\[
\epsilon=0.01
\]

failed badly.

Insufficient redistribution left photons trapped at high-opacity UV wavelengths. Their diffusion times became long and adiabatic losses large.

Thus the SN-Ia result was not simply:

> lower \(\epsilon\) is better.

Instead, it showed that substantial wavelength redistribution is necessary, but complete thermal redistribution is not necessarily optimal.

Schematically:

\[
\begin{array}{c|c}
\epsilon & \text{SN-Ia behavior in Kasen et al.}\\
\hline
0.01 & \text{far too little redistribution}\\
0.3 & \text{good approximation to direct fluorescence}\\
1.0 & \text{reasonable, but somewhat too red}
\end{array}
\]

Kasen et al. reported B-band errors of order \(0.1\) mag for the TLA approximation in this test.

---

# 9. Why this does not imply \(\epsilon=0.3\) is universal

The success of a scalar \(\epsilon\) in a Type-Ia calculation can depend on:

- the iron-peak atomic structure;
- line density;
- the frequency distribution of lines;
- repeated line interactions;
- optical depths;
- temperature;
- ionization state;
- incident radiation field;
- ejecta velocity structure.

Kasen et al. argued that the extreme complexity of iron-peak atomic structure and the high interaction rate can drive the radiation toward a quasi-equilibrium, which helps explain why a crude thermal redistribution model can work surprisingly well.

Lanthanides also have extremely complex atomic structures, but the detailed structure is radically different because of their open \(4f\) shells and enormous branching networks.

Therefore the natural question is empirical:

\[
\boxed{
\text{Does the same effective scalar closure work for La II, Ce II, etc.?}
}
\]

There is no reason to assume in advance that

\[
\epsilon_{\rm best}^{\rm Fe}
=
\epsilon_{\rm best}^{\rm La}
=
\epsilon_{\rm best}^{\rm Ce}.
\]

---

# 10. Modern review perspective: Jerkstrand (2025)

The 2025 Living Reviews article by Anders Jerkstrand is the most useful modern review for placing the subject in contemporary supernova and kilonova spectral synthesis.

The review emphasizes several points directly relevant to the present project.

## 10.1 True absorption can be tiny

For allowed line transitions, the true collisional thermalisation probability can be roughly

\[
10^{-6}-10^{-4}.
\]

---

## 10.2 LTE codes often use much larger effective values

Many LTE codes instead adopt

\[
p_{bb}^{\rm abs}
=
\epsilon_0,
\]

with a constant phenomenological value.

The review notes that surprisingly large values, around unity in some applications, can be necessary to obtain realistic spectra and light curves.

---

## 10.3 Why?

Because without explicit fluorescence, the simplified treatment has too little wavelength transformation.

Thus a large \(\epsilon_0\) can act as a heuristic surrogate for omitted fluorescence.

The price is that it can produce an incorrect energy balance because it treats more radiative energy as truly thermalized than the microphysics warrants.

This is particularly relevant for kilonovae, where radioactive material can be distributed throughout the ejecta.

---

# 11. Multilevel atoms: why one scalar may fail

A real ion does not have one line and two states.

For one upper level \(u\),

\[
u
\rightarrow
\{k_1,k_2,k_3,\ldots\}
\]

with many \(A_{uk}\) values and many optical depths.

Moreover, after one de-excitation,

\[
u\rightarrow k_1,
\]

the atom may remain excited, allowing another transition.

Thus the actual redistribution is a network:

\[
\boxed{
P(\lambda_{\rm out}\mid\lambda_{\rm in})
}
\]

rather than one probability.

One scalar \(\epsilon\) can at best approximate some moments of this redistribution operator.

It might reproduce:

- total optical flux;

while failing to reproduce:

- UV flux;
- NIR flux;
- detailed spectral shape;
- color;
- wavelength-dependent redistribution kernel.

This is why the planned experiment must compare **multiple wavelength bands and the full spectrum**, not just one optical band.

---

# 12. Macro-atoms and explicit Monte Carlo branching

Lucy (2002, 2003) developed the macro-atom formalism to handle multilevel atomic redistribution consistently in Monte Carlo radiative transfer.

The conceptual change is profound.

Instead of saying:

```text
line interaction
   |
epsilon
   |
scatter OR thermalize
```

a macro-atom says roughly:

```text
packet activates atomic state
          |
          v
internal transitions through the atomic network
          |
          +--> radiative deactivation
          |
          +--> collisional transition
          |
          +--> another internal state
          |
          v
packet eventually leaves in a channel dictated
by the statistical-equilibrium transition probabilities
```

The transition probabilities are constructed so that, statistically, the Monte Carlo process reproduces the emissivity of a gas in statistical equilibrium.

For understanding modern TARDIS/ARTIS-style treatments, Lucy's work and the Noebauer & Sim review are the natural references.

---

# 13. Energy packets versus photon-number packets

This is especially important for fluorescence.

Kasen et al. use monochromatic **equal-energy packets**.

If a packet has energy \(E_p\) and wavelength \(\lambda\), it represents

\[
N_p
=
\frac{E_p\lambda}{hc}
\]

physical photons.

When fluorescence changes the wavelength, the packet energy remains fixed while the number of represented photons changes.

This guarantees energy conservation through a redistribution event.

Example:

\[
\lambda_{\rm in}=2000~\unicode{x212B},
\qquad
\lambda_{\rm out}=4000~\unicode{x212B}.
\]

A physical photon has

\[
E_\gamma=\frac{hc}{\lambda}.
\]

Thus an individual 4000 Å photon has half the energy of a 2000 Å photon.

A constant-energy Monte Carlo packet therefore represents twice as many photons after the wavelength doubles.

Any photon-number packet implementation must account for this frequency-dependent energy weight explicitly.

Therefore the current branching instrument should distinguish:

\[
N_{\rm packet}
\]

from

\[
E_{\rm packet}.
\]

The correct test is energy conservation, not photon-number conservation alone.

---

# 14. How to interpret different \(\epsilon\) values in the planned experiment

A useful practical table is:

| \(\epsilon\) | Approximate TLA interpretation | What it does **not** necessarily mean |
|---:|---|---|
| 0 | pure scattering | not “zero fluorescence” in a real atom |
| 0.01 | almost pure scattering | not the same as the true \(10^{-4}\) collisional destruction probability in a multilevel fluorescent ion |
| 0.1 | weak effective redistribution | not necessarily 10% physical thermalisation |
| 0.3 | moderate effective redistribution; worked well in Kasen's SN-Ia test | not a universal atomic constant |
| 0.5 | stronger effective redistribution | not automatically closer to real branching |
| 0.7–0.9 | strongly thermalizing closure | still phenomenological if replacing fluorescence |
| 1 | complete thermal redistribution | does not imply real collisions dominate |

The key warning is:

\[
\boxed{
\epsilon_{\rm best}
\text{ is an effective closure parameter, not automatically a microscopic probability.}
}
\]

---

# 15. What the current La II Phase-1 result already says

For the current reference problem:

\[
F_{\rm Sob,abs}=0.183,
\]

\[
F_{\rm exp,abs}=0.344,
\]

\[
F_{\rm Sob+thermal}=0.257,
\]

\[
F_{\rm exp+thermal}=0.412,
\]

\[
F_{\rm Sob+branch}=0.660\pm0.003.
\]

The pure-absorption expansion-opacity leg is brighter than the per-line Sobolev reference:

\[
\frac{0.344}{0.183}-1
\approx+88\%.
\]

But complete thermal redistribution gives much less optical flux than explicit La II branching:

\[
\frac{0.412}{0.660}-1
\approx-37.6\%.
\]

Therefore the current result establishes:

\[
\boxed{
\epsilon=1
\text{ is not a good representation of the direct La II branching result in this band.}
}
\]

It does **not** yet establish:

\[
\boxed{
\text{no scalar }\epsilon\text{ can represent La II.}
}
\]

---

# 16. The decisive experiment

The literature makes the next experiment very clear.

Run

\[
\epsilon=
0,\ 0.1,\ 0.2,\ 0.3,\ 0.5,\ 0.7,\ 0.9,\ 1.
\]

Compare each TLA calculation with explicit branching.

For band \(b\), define

\[
\epsilon_{\rm best}^{(b)}
=
\arg\min_\epsilon
\left|
F_b^{\rm TLA}(\epsilon)-F_b^{\rm branch}
\right|.
\]

Then ask:

\[
\boxed{
\epsilon_{\rm best}^{\rm UV}
\stackrel{?}{=}
\epsilon_{\rm best}^{\rm optical}
\stackrel{?}{=}
\epsilon_{\rm best}^{\rm NIR}.
}
\]

---

# 17. Possible scientific outcomes

## 17.1 One \(\epsilon\) works across the spectrum

Suppose

\[
\epsilon_{\rm best}\approx0.3
\]

reproduces UV, optical, and NIR fluxes reasonably well.

Then:

\[
\boxed{
\text{the SN-Ia scalar-TLA result largely carries over to La II.}
}
\]

The important conclusion would be that \(\epsilon=1\) is inappropriate, not that scalar TLA itself is fundamentally inadequate.

---

## 17.2 Different bands require different \(\epsilon\)

Suppose

\[
\epsilon_{\rm best}^{\rm UV}=0.1,
\]

\[
\epsilon_{\rm best}^{\rm optical}=0.3,
\]

\[
\epsilon_{\rm best}^{\rm NIR}=0.8.
\]

Then:

\[
\boxed{
\text{no scalar }\epsilon\text{ reproduces the La II redistribution spectrum.}
}
\]

This would be a qualitatively stronger result.

It would imply that line identity carries wavelength-dependent information that cannot be represented by one thermalisation probability.

---

## 17.3 One \(\epsilon\) works for La but not Ce

If

\[
\epsilon_{\rm best}^{\rm La}
\neq
\epsilon_{\rm best}^{\rm Ce},
\]

then the closure is composition dependent.

This would matter directly for kilonova ejecta because the composition changes with electron fraction, epoch, and ejecta component.

---

# 18. Why the redistribution matrix is more fundamental than \(\epsilon_{\rm best}\)

The scalar \(\epsilon\) compresses the full redistribution physics.

The real quantity is closer to

\[
\boxed{
P(\lambda_{\rm out}\mid\lambda_{\rm in}).
}
\]

A direct branching calculation can estimate this matrix.

For example, the experiment may reveal strong structure such as

\[
2000\ {\rm \AA}
\rightarrow
3800\text{--}5000\ {\rm \AA}.
\]

Thermal re-emission at 3000 K may instead preferentially distribute energy much farther into the red/NIR.

Two models could accidentally have the same integrated optical flux while having very different matrices.

Thus the hierarchy of diagnostics should be:

1. band flux;
2. full spectral residual;
3. redistribution matrix;
4. cascade pathways.

---

# 19. Connection with expansion opacity

The \(\epsilon\) issue and expansion opacity are related but should not be conflated.

There are at least two separate approximations:

\[
\boxed{
\text{opacity / interaction representation}
}
\]

and

\[
\boxed{
\text{post-interaction redistribution closure}.
}
\]

Paper I studies primarily the first.

Paper II studies primarily the second.

A useful hierarchy is

\[
\text{resolved finite profile}
\rightarrow
\text{per-line Sobolev}
\rightarrow
\text{expansion opacity}
\rightarrow
\text{TLA redistribution}.
\]

Explicit branching instead follows

\[
\text{per-line Sobolev interaction}
\rightarrow
\text{atomic upper level}
\rightarrow
A\beta\text{-weighted branching/cascade}.
\]

Thus a comparison such as

\[
\text{Sobolev+branching}
\]

versus

\[
\text{expansion+thermal}
\]

contains both an opacity-representation difference and a redistribution-closure difference.

The experiment should therefore retain control legs that allow these effects to be separated.

---

# 20. Relation to line-binned opacity work

Fontes et al. (2020) developed a line-binned opacity treatment for kilonova calculations and compared it with traditional expansion opacity and a continuous Monte Carlo Sobolev approach.

Their work is important because it directly addresses how a dense atomic line list can be represented efficiently in a kilonova transport calculation.

For the present project, however, the key complementary issue is that a comparison of opacity representations does not by itself determine whether a simple post-interaction thermalisation closure reproduces the actual multilevel fluorescence network.

This is why explicit branching remains a separate question.

---

# 21. Relation to Monte Carlo radiative-transfer reviews

Noebauer & Sim (2019) is the best broad modern review of Monte Carlo radiative transfer.

It is particularly useful for:

- packet-based radiation transport;
- estimators;
- Monte Carlo noise;
- indivisible energy packets;
- line interactions;
- macro-atoms;
- time dependence;
- expanding media.

For the present project, this review bridges the conceptual gap between textbook source-function language and the algorithms actually used in modern supernova/kilonova codes.

---

# 22. Recommended reading sequence

A full textbook reading is unnecessary. The following sequence is more efficient.

## Stage 1 — understand the two-level atom

### Dimitri Mihalas, *Stellar Atmospheres*, 2nd ed. (1978)

Prioritize:

- Chapter 11: **Non-LTE Line Transfer: The Two-Level Atom**
  - diffusion;
  - destruction;
  - escape;
  - thermalization;
  - source function;
  - thermalization depth.

Then:

- Chapter 12: multilevel atom, especially the equivalent-two-level-atom approach;
- Chapter 13: partial frequency redistribution;
- Chapter 14: moving atmospheres;
- §14-2: Sobolev theory;
- “Escape and Thermalization in an Expanding Medium.”

**Goal:** understand why

\[
S=(1-\epsilon)J+\epsilon B
\]

is physically meaningful.

---

## Stage 2 — modernize the formalism

### Hubeny & Mihalas, *Theory of Stellar Atmospheres* (2014)

Prioritize:

- Chapter 14: NLTE two-level and multilevel atoms;
- Chapter 15: partial redistribution;
- Chapter 19: extended and expanding atmospheres.

**Goal:** understand why a real multilevel atom is not generally reducible to one fixed two-level source function.

---

## Stage 3 — connect to moving media and Sobolev theory

### John I. Castor, *Radiation Hydrodynamics* (2004)

Prioritize:

- Chapter 5: steady-state transfer;
- Chapter 6: comoving-frame picture;
- Chapter 8: radiation–matter interactions;
- Chapter 9: spectral line transport.

**Goal:** connect local line physics to moving, expanding ejecta.

---

## Stage 4 — understand the actual SN-Ia \(\epsilon\) experiment

### Kasen, Thomas & Nugent (2006)

Read:

- §II.2: opacities and emissivities;
- equations for Sobolev optical depth and escape probability;
- TLA source function;
- §II.3: equal-energy Monte Carlo packets;
- §III.6: **Fluorescence Versus Two-Level Atom Redistribution**.

Pay particular attention to:

\[
\epsilon=1,\ 0.3,\ 0.01
\]

and to

\[
p_{\rm fluor}.
\]

**Goal:** understand exactly what our La II experiment is generalizing.

---

## Stage 5 — understand explicit multilevel Monte Carlo transport

### Lucy (2002, 2003)

Read after Kasen.

**Goal:** understand why macro-atoms preserve atomic-state information that a scalar TLA closure discards.

---

## Stage 6 — read the modern overview

### Noebauer & Sim (2019)

Use as a bridge between algorithms.

### Jerkstrand (2025)

Use for the current state of SN/KN spectral synthesis, especially:

- true line thermalisation;
- heuristic \(\epsilon_0\);
- fluorescence;
- LTE versus NLTE;
- expansion opacity;
- energy-equation implications.

**Goal:** understand why the present experiment matters now rather than only historically.

---

# 23. A minimal one-day reading plan

If time is limited:

### 1. Mihalas Chapter 11
Focus on

\[
S=(1-\epsilon)J+\epsilon B.
\]

### 2. Kasen et al. §II.2 and §III.6
Focus on

\[
\epsilon=1,\ 0.3,\ 0.01
\]

and direct fluorescence.

### 3. Jerkstrand 2025, discussion around true line absorption and \(\epsilon_0\)
Focus on the distinction

\[
p_{\rm abs,true}\ll1
\]

versus

\[
\epsilon_0\sim1.
\]

### 4. Lucy / Noebauer & Sim
Focus on energy packets and macro-atoms.

After those readings, the physical interpretation of the planned \(\epsilon\)-sweep should be much clearer.

---

# 24. Key equations to keep on one page

## TLA source function

\[
\boxed{
S_\nu
=
(1-\epsilon)J_\nu+\epsilon B_\nu
}
\]

## Sobolev line interaction probability

\[
\boxed{
P_{\rm int}=1-e^{-\tau_S}
}
\]

## Sobolev escape probability

\[
\boxed{
\beta_{\rm esc}
=
\frac{1-e^{-\tau_S}}{\tau_S}
}
\]

## True line absorption probability, schematic multilevel form

\[
\boxed{
p_{\rm abs,true}
\sim
\frac{\text{collisional de-excitation}}
{\text{collisional de-excitation}+\text{escaping radiative de-excitation}}
}
\]

## Direct fluorescence probability

\[
\boxed{
p_{\rm fluor}
=
\frac{
\sum_{k\ne l}\beta_{uk}A_{uk}
}{
n_e\sum_kC_{uk}+\sum_k\beta_{uk}A_{uk}
}
}
\]

## Expansion opacity

\[
\boxed{
\alpha_{\rm exp}
=
\frac{1}{ct}
\sum_i
\frac{\lambda_i}{\Delta\lambda}
\left(1-e^{-\tau_i}\right)
}
\]

## Full redistribution object

\[
\boxed{
P(\lambda_{\rm out}\mid\lambda_{\rm in})
}
\]

These equations represent increasing levels of physical information.

---

# 25. Terminology recommendations for the papers

To avoid confusion, use explicit terms.

Prefer:

- **true collisional thermalisation probability**
  \[
  p_{\rm abs,true}
  \]
- **TLA redistribution parameter**
  \[
  \epsilon
  \]
- **effective phenomenological thermalisation parameter**
  \[
  \epsilon_0
  \]
- **line-specific fluorescence probability**
  \[
  p_{{\rm fluor},i}
  \]
- **Sobolev escape probability**
  \[
  \beta_{\rm esc}
  \]

Avoid writing simply:

> “\(\epsilon\) is the physical absorption probability”

unless the exact context is the microscopic two-level problem.

For the current Paper-II experiment, wording such as the following is safer:

> “We treat \(\epsilon\) as the scalar redistribution parameter of the two-level-atom closure; it should not be identified with the much smaller microscopic collisional destruction probability.”

---

# 26. What the literature implies for the experimental design

The literature motivates the following sequence:

```text
1. Verify energy-packet bookkeeping.
          |
          v
2. Verify beta-weighted branching.
          |
          v
3. Reach the narrow-line Sobolev regime (v_D <= 10 km/s).
          |
          v
4. Sweep epsilon from 0 to 1.
          |
          v
5. Compare UV / optical / NIR separately.
          |
          v
6. Compare the full spectrum.
          |
          v
7. Build P(lambda_out | lambda_in).
          |
          v
8. Test line-dependent p_fluor.
          |
          v
9. Repeat for Ce II.
```

The decisive question after Step 4–6 is:

\[
\boxed{
\text{Is one scalar }\epsilon\text{ enough?}
}
\]

Only after that should the experiment be generalized to:

- \(v_{\rm bulk}\sim0.1c\);
- multi-ion mixtures;
- realistic source spectra;
- LTE temperature feedback;
- NLTE.

---

# 27. Primary references

## [R1] Mihalas (1978) — classic textbook

**Mihalas, D. (1978). _Stellar Atmospheres_, 2nd ed. W. H. Freeman.**

Key topics:

- two-level atom;
- thermalization;
- multilevel atom;
- partial redistribution;
- moving atmospheres;
- Sobolev theory.

Bibliographic record:

- ISBN: 0-7167-0359-9
- WorldCat / library record: https://search.worldcat.org/title/601343352
- Open Library record: https://openlibrary.org/books/OL17758884M/Stellar_atmospheres

---

## [R2] Hubeny & Mihalas (2014) — modern textbook

**Hubeny, I. & Mihalas, D. (2014). _Theory of Stellar Atmospheres: An Introduction to Astrophysical Non-equilibrium Quantitative Spectroscopic Analysis_. Princeton University Press.**

Key chapters:

- Ch. 14: NLTE two-level and multilevel atoms;
- Ch. 15: partial redistribution;
- Ch. 19: extended and expanding atmospheres.

Publisher/book information:

https://books.google.com/books?id=VA_rAwAAQBAJ

ISBN: 978-0-691-16328-4

---

## [R3] Castor (2004) — moving media / radiation hydrodynamics

**Castor, J. I. (2004). _Radiation Hydrodynamics_. Cambridge University Press.**

DOI:

https://doi.org/10.1017/CBO9780511536182

Publisher page:

https://www.cambridge.org/core/books/radiation-hydrodynamics/A4D7F2A12AE2929A6059D38190234352

Key topics:

- comoving-frame transfer;
- radiation–matter interactions;
- spectral line transport;
- expanding media.

---

## [R4] Kasen, Thomas & Nugent (2006) — essential SN-Ia benchmark

**Kasen, D., Thomas, R. C. & Nugent, P. (2006). “Time-dependent Monte Carlo radiative transfer calculations for three-dimensional supernova spectra, light curves, and polarization.” _The Astrophysical Journal_, 651, 366–380.**

DOI:

https://doi.org/10.1086/506190

Open-access record:

https://escholarship.org/uc/item/4cb8m7zq

arXiv:

https://arxiv.org/abs/astro-ph/0606111

Most important sections:

- §II.2: Sobolev optical depth, escape probability, expansion opacity, TLA source function;
- §II.3: equal-energy packet transport;
- §III.6: direct fluorescence versus TLA with
  \[
  \epsilon=1,\ 0.3,\ 0.01.
  \]

---

## [R5] Lucy (2002) — macro-atom transition probabilities

**Lucy, L. B. (2002). “Monte Carlo transition probabilities.” _Astronomy & Astrophysics_, 384, 725–735.**

DOI:

https://doi.org/10.1051/0004-6361:20011756

arXiv:

https://arxiv.org/abs/astro-ph/0107377

Main relevance:

- indivisible energy packets;
- multilevel transition probabilities;
- statistically correct emissivity;
- foundation of macro-atom-style line redistribution.

---

## [R6] Lucy (2003) — macro-atom continuation

**Lucy, L. B. (2003). “Monte Carlo transition probabilities. II.” _Astronomy & Astrophysics_, 403, 261–275.**

DOI:

https://doi.org/10.1051/0004-6361:20030357

arXiv:

https://arxiv.org/abs/astro-ph/0303202

---

## [R7] Noebauer & Sim (2019) — modern Monte Carlo review

**Noebauer, U. M. & Sim, S. A. (2019). “Monte Carlo radiative transfer.” _Living Reviews in Computational Astrophysics_, 5, 1.**

DOI:

https://doi.org/10.1007/s41115-019-0004-9

Open article:

https://link.springer.com/article/10.1007/s41115-019-0004-9

arXiv:

https://arxiv.org/abs/1907.09840

Main relevance:

- modern MCRT overview;
- packet methods;
- expanding media;
- line transfer;
- macro-atoms;
- estimators and Monte Carlo noise.

---

## [R8] Jerkstrand (2025) — modern SN/KN spectral-synthesis review

**Jerkstrand, A. (2025). “Spectral synthesis techniques for supernovae and kilonovae.” _Living Reviews in Computational Astrophysics_, 11, 1.**

DOI:

https://doi.org/10.1007/s41115-025-00022-2

Open article:

https://link.springer.com/article/10.1007/s41115-025-00022-2

PMC:

https://pmc.ncbi.nlm.nih.gov/articles/PMC12334460/

Main relevance:

- modern comparison of SN/KN spectral-synthesis techniques;
- true line thermalisation probability;
- phenomenological \(\epsilon_0\);
- fluorescence;
- expansion opacity;
- LTE/NLTE energy-equation issues.

Particularly relevant discussion:

\[
p_{bb}^{\rm abs,true}\sim10^{-6}-10^{-4}
\]

versus the much larger constant \(\epsilon_0\) used in approximate LTE line treatments.

---

## [R9] Fontes et al. (2020) — kilonova line-binned opacity

**Fontes, C. J., Fryer, C. L., Hungerford, A. L., Wollaeger, R. T. & Korobkin, O. (2020). “A line-binned treatment of opacities for the spectra and light curves from neutron star mergers.” _Monthly Notices of the Royal Astronomical Society_, 493, 4143.**

arXiv:

https://arxiv.org/abs/1904.08781

Main relevance:

- line-binned opacity;
- comparison with expansion opacity;
- continuous Monte Carlo Sobolev calculations;
- kilonova-specific opacity representation.

---

# 28. Additional historical references

These are useful if a deeper historical treatment is needed.

## Sobolev theory

**Sobolev, V. V.**  
Classical work on radiative transfer in moving media and the approximation now bearing his name.

Mihalas (1978), Castor (2004), and Kasen et al. (2006) provide more accessible derivations for the present purpose.

---

## Castor (1970)

**Castor, J. I. (1970), MNRAS, 149, 111.**

Important early formulation of Sobolev/escape-probability line transfer in expanding atmospheres.

---

## Karp et al. (1977)

**Karp, A. H., Lasher, G., Chan, K. L. & Salpeter, E. E. (1977).**

Classic expansion-opacity work for rapidly expanding supernova envelopes.

This is part of the historical foundation of the binned expansion-opacity formalism later reformulated by Eastman & Pinto.

---

## Eastman & Pinto (1993)

**Eastman, R. G. & Pinto, P. A. (1993). _The Astrophysical Journal_, 412, 731.**

DOI:

https://doi.org/10.1086/172957

A central reference for the expansion-opacity formalism used in supernova modeling.

---

# 29. Suggested citation language for Paper II

A compact literature paragraph could eventually say:

> In a classical two-level atom, the line source function can be written \(S=(1-\epsilon)J+\epsilon B\), with \(\epsilon\) measuring true thermal coupling (e.g. Mihalas 1978; Hubeny & Mihalas 2014). In dilute supernova and kilonova ejecta, however, the microscopic collisional destruction probability of allowed lines is typically very small, while LTE transport calculations often adopt much larger effective values \(\epsilon_0\sim O(1)\) to mimic wavelength redistribution by fluorescence that is absent from the simplified line treatment (Kasen et al. 2006; Jerkstrand 2025). Kasen et al. showed for iron-peak SN-Ia ejecta that a scalar TLA treatment with \(\epsilon\simeq0.3\) could reproduce direct fluorescence substantially better than either nearly pure scattering or complete thermal redistribution. We test whether this scalar closure remains adequate for the much denser and differently structured branching networks of lanthanide ions.

This should be refined once the \(\epsilon\)-sweep results are available.

---

# 30. Bottom line

The most important conceptual distinction is

\[
\boxed{
\epsilon_{\rm true}
\neq
\epsilon_{\rm effective}
}
\]

in the regime relevant to simplified SN/KN line transport.

The classical two-level atom teaches what \(\epsilon\) means as thermal coupling.

The supernova literature then repurposes a large effective \(\epsilon\) as a computational closure for missing fluorescence.

The real multilevel atom instead carries a redistribution network determined by

\[
A_{uk},
\qquad
\tau_{uk},
\qquad
\beta_{uk},
\qquad
\text{level connectivity},
\]

which produces

\[
P(\lambda_{\rm out}\mid\lambda_{\rm in}).
\]

The present experiment asks whether that network can still be represented by one scalar.

The immediate research question is therefore:

\[
\boxed{
\exists\epsilon:
F_\lambda^{\rm TLA}(\epsilon)
\approx
F_\lambda^{\rm direct\ branching}
\quad\text{over the full spectrum?}
}
\]

If yes, the classic SN-Ia result transfers surprisingly well to lanthanides.

If no, then the complexity retained by line identity is not reducible to a universal scalar thermalisation closure, motivating line-dependent or branching-aware treatments.
