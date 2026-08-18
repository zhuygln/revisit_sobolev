# Cover letter — MNRAS submission

Dear Editor,

Please consider the enclosed manuscript, *Two Approximations Under One Name:
Separating Sobolev Line Transfer from Expansion Opacity in Kilonova Ejecta*,
for publication in MNRAS.

Kilonova radiative-transfer calculations cannot resolve the millions of
lanthanide lines that shape their spectra, and so rely on a line-transfer
shortcut that the literature almost universally calls "the Sobolev
approximation". In the codes that work from binned opacities, that phrase
denotes two distinct approximations applied together: Sobolev's own treatment
of an isolated resonance, and the expansion-opacity closure built on top of
it. The manuscript separates them and measures each against frequency-resolved
transport, using a controlled harness in which the same model, atomic data and
level populations are pushed through analytic theory, a purpose-written
deterministic solver, and a Monte Carlo code in two line-treatment modes.

I believe the work suits MNRAS because it is a methods and validity study of a
numerical closure in wide use, of direct interest to the radiative-transfer and
transient-modelling communities the journal serves, and because its central
claim is a controlled measurement rather than an astrophysical inference.

Two points may be useful for selecting a referee and for setting expectations.

First, the configuration is deliberately restricted — a one-dimensional
homologous shell, pure absorption with fixed LTE populations, sub-relativistic
velocities, and calibrated data for two lanthanide ions. That restriction is
not a limitation we tolerated but the design: it is what makes the two errors
separable at all. The manuscript states throughout what the restriction
licenses and what it does not, and in particular does not claim that the
measured band-flux error is the error of a fully mixed, energy-conserving
kilonova spectrum.

Second, the manuscript reports several of its own corrections, including one
result that an earlier draft had wrong and that we retract explicitly in the
text. We have kept those passages rather than quietly removing them, because
in a study whose subject is the reliability of an approximation, the
distinction between a measurement and a convention artifact is part of the
result.

A companion study, examining how much of the effect survives radiative
redistribution, is in preparation; the present paper is self-contained and does
not depend on it.

All data, generating scripts and the full analysis history are publicly
available at `https://github.com/zhuygln/revisit_sobolev`.

The manuscript has not been published elsewhere and is not under consideration
by another journal.

Yours sincerely,

Yonglin Zhu
