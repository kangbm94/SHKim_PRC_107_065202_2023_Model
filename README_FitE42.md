# Simultaneous ROOT fit of E42 polarization and forward cross section

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

The default call simultaneously fits both datasets. To run the earlier
polarization-only objective for comparison, use

    .x FitE42.cc+(6,false)

The twelve discrete hypotheses are fitted in parallel by default. On macOS and
Linux, `nWorkers=0` uses the available CPU count, capped at the number of
hypotheses. To request four workers explicitly, use

    .x FitE42.cc+(6,true,4)

For sequential execution, use

    .x FitE42.cc+(6,true,1)

Parallelism is implemented with isolated worker processes because legacy
TMinuit keeps global state and is unsafe for concurrent hypothesis fits in
threads. Workers only minimize and return numerical results. The parent process
alone writes TSV files and creates ROOT canvases. Windows falls back to
sequential execution.

### Structural test mode

`gTestMode` near the top of `FitE42.cc` is a smoke-test switch for the
multistart, worker-process, result-writing, and plotting workflow. When it is
nonzero, Minuit uses a small toy objective involving the first three parameters
instead of the polarization-plus-forward-cross-section chi-square. Leave it
enabled while diagnosing program structure, but set it to `0` before producing
physical fit parameters or comparing BIC values. Parameter and chi-square files
written in test mode are diagnostic outputs, not fit results.

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

The default objective contains the 15 E42 polarization points and 25 forward
differential-cross-section points:

    chi2_total = chi2_polarization + chi2_forward.

Each residual is divided by the corresponding quoted experimental uncertainty;
no arbitrary relative dataset weight is applied. The polarization in an angular
bin is evaluated as

    integral[P(cos theta) * dsigma/dOmega * dcos theta]
    ---------------------------------------------------
           integral[dsigma/dOmega * dcos theta]

using 12-point Gauss-Legendre quadrature. Each forward prediction is the
model's laboratory-frame average over theta_lab = 0--20 degrees, matching the
definition of the supplied forward dataset. The angular integral uses 16
subintervals; this was checked against 64 subintervals at representative beam
momenta and agreed to numerical precision.

Because the hypotheses without and with Sigma(2230) contain 9 and 12 fitting
parameters respectively, raw chi-square alone favors the larger model. The
macro therefore also writes

    BIC = chi2 + number_of_parameters * log(40).

In polarization-only mode N=15 is used instead. The result table records total,
polarization, and forward-cross-section chi-square separately.

## Outputs

- FitE42_hypothesis_results.tsv: total and per-dataset chi-square, BIC, status,
  point counts, and all fitted parameters for every hypothesis.
- FitE42_best_fit.tsv: the minimum-BIC result and its parameter bounds.
- FitE42_hypothesis_comparison.pdf/png: chi-square and BIC comparison.
- FitE42Plots/FitE42_*.pdf/png: one polarization plus forward-cross-section
  canvas for each hypothesis. The forward panel contains the black coherent
  total and dashed individual contributions.
- FitE42Plots/FitE42_spin_observables_*.pdf/png: the spin-non-flip fraction
  and relative f--g phase in the three E42 cos(theta) intervals.
- FitE42Plots/FitE42_amplitudes_*_cth_*.pdf/png: coherent, bin-averaged real
  and imaginary f and g amplitudes for the full model and requested individual
  contributions. Real parts are solid and imaginary parts are dashed.
- FitE42Plots/Differential/<hypothesis>/FitE42_differential_*.pdf/png:
  one CM angular-distribution monitor for each experiment and sqrt(s) dataset.
  Each canvas overlays data, the black coherent total, dashed individual
  components, and an inset restricted to 0.5 < cos(theta_CM) < 1.
- FitE42Parameters/FitE42_parameters_<hypothesis>.tsv: the best-fit minimum,
  chi-square decomposition, BIC, status, and parameter values/errors/bounds for
  every spin-parity and Sigma(2230) configuration.

The angular differential-cross-section datasets are monitoring plots only and
do not enter chi-square. The forward integrated dataset still enters the
simultaneous fit. An individual dashed curve is the cross section obtained from
that amplitude alone, |M_i|^2; dashed curves do not sum to the coherent black
curve because the latter also contains interference terms and other enabled
background amplitudes.

The supplied differential-cross-section files contain 26 distinct
experiment/energy datasets. A complete 12-hypothesis scan therefore writes 312
differential canvases, each in both PDF and PNG format. A failure in one
monitoring dataset is reported as a warning and does not prevent later datasets
or hypotheses from being plotted.

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

PDG sources:

- https://pdg.lbl.gov/2025/listings/rpp2025-list-lambda-2100.pdf
- https://pdg.lbl.gov/2025/listings/rpp2025-list-sigma-2030.pdf
- https://pdg.lbl.gov/2025/listings/rpp2025-list-sigma-2250.pdf
