# Kim et al. RPR calculation and E42 fit scaffold

This package reproduces the elementary \(K^-p\to K^+\Xi^-\) calculation of
S.-H. Kim *et al.*, *Phys. Rev. C* **107**, 065202 (2023), through the
Reggeized \(u\)-channel and ground-state/resonant \(s\)-channel amplitudes in
Eqs. (6)--(16).  It also provides the acceptance-folded polarization estimator
needed to fit the model to E42 quasifree data.

## What is implemented

- Reggeized \(\Lambda\), \(\Sigma\), and \(\Sigma(1385)\) exchange.
- Ground-state \(s\)-channel \(\Lambda\) and \(\Sigma\).
- The five resonances in Table II: \(\Lambda(1890)\), \(\Lambda(2100)\),
  \(\Sigma(2030)\), \(\Sigma(2230)\), and \(\Sigma(2250)\).
- Spin-3/2 and momentum-contracted spin-7/2 projectors.
- Differential/total cross sections and recoil polarization.
- Event-by-event acceptance folding and amplitude-scale fitting.

The meson-baryon rescattering calculation in Eqs. (17)--(28) is **not** yet
included.  The paper does not fully specify the numerical/isospin implementation
needed to reproduce that loop unambiguously.  The polarization plots in Fig. 7
use the \(u+s\) calculation, so this is the relevant starting point for E42.

## Conventions fixed by validation

The implementation uses \(\bar u u=1\), the cross-section normalization below
Eq. (8), and the E42 normal
\(\hat y\parallel\mathbf{k}_{K^-}\times\mathbf{k}_{K^+}\).  The normal used in
the paper's Fig. 7 has the opposite sign.  Direct comparison with Figs. 4, 5,
and 7 also selects:

- the squared \(\eta_H\) in Eq. (7);
- relative charged-channel signs \(I_\Lambda=-1\), \(I_\Sigma=+1\);
- \(F_R=\exp[-(s-M_R^2)^2/\Lambda_R^4]\), once per diagram;
- the resonance pole mass in the high-spin projector.

These choices are explicit in `ModelParameters` rather than hidden constants.

## Install and reproduce the paper figures

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce_figures.py --output-dir output
pytest -q
```

The outputs `reproduce_fig4a.png`, `reproduce_fig5a.png`, and
`reproduce_fig7a.png` validate the energy dependence, angular distributions,
and polarized differential cross section.

## Correct E42 fit

E42 bins are quasifree and span \(\sqrt{s}\) and \(\cos\theta\).  The prediction
must therefore be folded over accepted events:

\[
P_b^{\rm pred}=
\frac{\sum_{i\in b}w_i\,(d\sigma_i/d\Omega)\,P_i}
     {\sum_{i\in b}w_i\,(d\sigma_i/d\Omega)}.
\]

Prepare two TSV/CSV files.  The accepted-event file needs `bin`, `sqrt_s`,
`cos_theta`, and optionally `weight`.  The polarization file needs `bin`,
`polarization`, `stat_error`, and optionally `syst_error`.  Templates are in
`examples/`.

```bash
python fit_e42.py \
  --events accepted_e42_events.tsv \
  --polarization e42_polarization.tsv \
  --parameters Lambda2100 Sigma2030 \
  --output-dir fit_output
```

An optional elementary cross-section table may be supplied with
`--cross-section`; its columns are `sqrt_s`, `cos_theta`, `dsigma_domega`, and
`error`.  Add `--fit-normalization` only when such data are present.

The fit parameters multiply complete complex amplitude components.  A scale of
one is the published value.  Polarization alone generally cannot identify many
independent amplitude scales, so start with one or two physically motivated
parameters and inspect the covariance and bound status.

## Bin-center diagnostic

`data/e42_plot_extracted.tsv` contains values digitized from the E42 filled
points.  They are not a substitute for the numerical result table or accepted
event sample.  The following command evaluates a deliberately constrained,
one-parameter diagnostic at bin centers:

```bash
python diagnose_e42_centers.py --output-dir output
```

The constraint
`Lambda2100_scale^2 + Sigma2030_scale^2 = 2` prevents a polarization-only fit
from increasing both amplitudes without limit.  Treat its output as a model-data
tension check, not as a quoted E42 parameter determination.

## Main files

- `rpr_model.py`: amplitudes, kinematics, observables, and parameters.
- `reproduce_figures.py`: validation plots corresponding to Figs. 4(a), 5(a),
  and 7(a).
- `fit_e42.py`: accepted-event folding and least-squares fit.
- `diagnose_e42_centers.py`: non-final comparison to digitized bin centers.
- `tests/`: algebra, kinematics, polarization, and amplitude-sum tests.
