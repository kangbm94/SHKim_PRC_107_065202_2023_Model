# ROOT fit of the E42 polarization

FitE42.cc is the primary fitting macro. It requires ROOT, Eigen, and the
headers in this directory; it does not require Python, NumPy, or SciPy.

From the FitSHKimE42 directory:

    root -l

In ROOT:

    .include /opt/homebrew/include/eigen3
    .x FitE42.cc+(6)

The integer is the number of deterministic multistart fits per hypothesis.
Use /usr/local/include/eigen3 on an Intel Mac when appropriate. If Eigen is
already in ROOT's include path, the .include command is unnecessary.

## Parameters and hypotheses

For every included resonance, Minuit fits three parameters:

1. a signed real amplitudeScale in [-3,3];
2. its mass within the PDG interval;
3. its width within the PDG interval.

The scale is an amplitude multiplier, not a positive-definite intensity.
Allowing either sign is necessary to test interference signs.

| state | mass range (GeV) | width range (GeV) |
|---|---:|---:|
| Sigma(2030) 7/2+ | 2.025--2.040 | 0.150--0.200 |
| Lambda(2100) 7/2- | 2.090--2.110 | 0.100--0.250 |
| Sigma(2250), unknown JP | 2.210--2.280 | 0.060--0.150 |
| Sigma(2230) 3/2+ | 2.213--2.267 | 0.295--0.395 |

The first three rows are PDG Breit-Wigner estimate ranges. Sigma(2230) is a
one-star state with one listed solution, M = 2240 +/- 27 MeV and
Gamma = 345 +/- 50 MeV; those quoted uncertainties define its bounds.

The macro separately fits all twelve hypotheses:

- Sigma(2250): 3/2+, 3/2-, 5/2+, 5/2-, 7/2+, or 7/2-;
- Sigma(2230): included or removed.

Spin and parity are discrete quantum numbers and are not varied continuously.

## Chi-square and angular averaging

Only the 15 E42 polarization points enter chi-square. The polarization in an
angular bin is evaluated as

    integral[P(cos theta) * dsigma/dOmega * dcos theta]
    ---------------------------------------------------
           integral[dsigma/dOmega * dcos theta]

using 12-point Gauss-Legendre quadrature. The old forward differential
cross-section points are drawn beside every polarization fit as an independent
validation and are not included in chi-square.

Because the hypotheses without and with Sigma(2230) contain 9 and 12 fitting
parameters respectively, raw chi-square alone favors the larger model. The
macro therefore also writes

    BIC = chi2 + number_of_parameters * log(15).

## Outputs

- FitE42_hypothesis_results.tsv: chi-square, BIC, status, and all fitted
  parameters for every hypothesis.
- FitE42_best_fit.tsv: the minimum-BIC result and its parameter bounds.
- FitE42_hypothesis_comparison.pdf/png: chi-square and BIC comparison.
- FitE42Plots/FitE42_*.pdf/png: one polarization plus forward-cross-section
  canvas for each hypothesis.
- FitE42Plots/FitE42_spin_observables_*.pdf/png: the spin-non-flip fraction
  and relative f--g phase in the three E42 cos(theta) intervals.
- FitE42Plots/FitE42_amplitudes_*_cth_*.pdf/png: coherent, bin-averaged real
  and imaginary f and g amplitudes for the full model and requested individual
  contributions. Real parts are solid and imaginary parts are dashed.

For a finite angular interval, the macro defines

    R_non-flip = sqrt[ integral |f|^2 dcos(theta)
                       / integral (|f|^2+|g|^2) dcos(theta) ]

and

    relative phase / pi = arg[ integral f* g dcos(theta) ] / pi.

The relative-phase convention is therefore arg(g)-arg(f), matching the model
polarization convention P = 2 Im(f* g)/(|f|^2+|g|^2). The displayed amplitude
itself is the coherent angular average. The black curve is the complete fitted
model sum. Component colors are: Sigma(2250) magenta, Sigma(2030) blue,
Lambda(2100) red, Sigma(2230) green when enabled, Lambda(1890) orange, and the
Reggeized u-channel Lambda cyan.

The supplied `FitE42_reference_results.tsv`,
`FitE42_reference_comparison.pdf/png`, and `FitE42Plots_reference/` files are
a numerical cross-check made with the same objective, hypotheses, and bounds.
The ROOT macro is self-contained and does not require that cross-check.

PDG sources:

- https://pdg.lbl.gov/2025/listings/rpp2025-list-lambda-2100.pdf
- https://pdg.lbl.gov/2025/listings/rpp2025-list-sigma-2030.pdf
- https://pdg.lbl.gov/2025/listings/rpp2025-list-sigma-2250.pdf
