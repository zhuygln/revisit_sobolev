# Cover letter — Nature Astronomy submission

Dear Editor,

Please consider the enclosed Article, *Coarse-grained line opacity leaves a
detectable chromatic signature in kilonovae that ejecta parameters cannot
mimic*, for publication in Nature Astronomy.

Kilonova light curves are the main route from a gravitational-wave merger to
the mass, velocity and lanthanide content of what it ejected, and the models
that read them cannot resolve the millions of lanthanide lines that shape the
spectrum. Many production light-curve and inference calculations therefore
coarse-grain the dense line opacity into bins, and inference frameworks often
absorb unresolved model deficiencies into a generic 0.5–1 mag systematic
term, without isolating how much comes from the line-transport closure
itself.

The manuscript tests that assumption directly. It runs one Monte Carlo
transport code, with the same atomic data, level populations and source
model, in a line-by-line reference mode and in three closure modes that
isolate the two approximations usually bundled under the name "expansion
opacity": the coarse-graining of the opacity itself, and the compression of
the fluorescence redistribution. The result is asymmetric. The redistribution
can be compressed to a small matrix without observable consequence. The
opacity coarse-graining cannot: across a complete grid of ejecta mass,
velocity and lanthanide fraction, the closure makes a heating-powered
kilonova too blue by one to three magnitudes in colour, with a coherent
optical-bright, near-infrared-faint pattern at every epoch. A first-order
sensitivity analysis shows that no shift of the three ejecta parameters
reproduces this pattern at any grid point, that a free luminosity history
does not absorb it in lanthanide-rich ejecta, and that it remains detectable
at the distance of GW170817 under realistic optical and near-infrared
sampling.

I believe the work belongs in Nature Astronomy because its message reaches
beyond kilonova codes. It shows, with a controlled measurement rather than an
argument, that a transport approximation used across much of the modelling
community produces a structured, sign-predictable error in the very
observables from which physical parameters are inferred, and that a
substantial part of what those inferences represent as an unstructured
model-error allowance can arise from a single, identifiable cause. The
broader methodological point may apply to other line-rich radiative-transfer
problems: a coarse-grained approximation can produce a coherent error that is
poorly represented by an unstructured nuisance term.

Three points may help in selecting referees and setting expectations.

First, the claim is differential. The manuscript reports the difference
between two treatments of the same model, and does not claim that the
reference treatment is the truth or that any published inference is wrong;
the comparison to the model-error allowance is a consistency statement.

Second, the experiment is deliberately restricted: a one-zone grey
photosphere illuminating a homologous shell, frozen ionization, four
lanthanide ions using calibrated atomic data. The manuscript
states what this licenses and what it does not, and names the extensions
(independent atomic data, a second code, a fit to AT2017gfo) as future work
rather than promising them.

Third, one pre-declared convergence criterion in the numerical control was
not met at a majority of the colours tested; the manuscript reports that
failure and quotes every amplitude with the uncertainty the test measured,
rather than to a precision the control does not support.

Every number in the manuscript is generated from a single frozen analysis
commit, tagged in the public repository, and the build refuses a manuscript
in which a quoted number differs from the frozen value. All data, code and
the complete analysis history are available at
`https://github.com/zhuygln/revisit_sobolev`.

The manuscript has not been published elsewhere and is not under
consideration by another journal.

Yours sincerely,

Yonglin Zhu
