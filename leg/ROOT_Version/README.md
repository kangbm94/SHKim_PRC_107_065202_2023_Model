# ROOT/C++ implementation of the S.-H. Kim RPR model

This is a header-only C++17 implementation of the elementary
`K- p -> K+ Xi-` amplitude in S.-H. Kim et al., PRC 107, 065202 (2023),
Eqs. (6)-(16).  It includes the Reggeized u-channel, ground-state s-channel,
and the five selected spin-3/2 and spin-7/2 resonances.  It does not include
the paper's meson-baryon rescattering loop.

## Run in ROOT

Eigen3 is the matrix package.  It is header-only and must be available as
either `Eigen/Dense` or `eigen3/Eigen/Dense`.

```bash
root -l
root [0] .L RunRPR.C+
root [1] RunRPR(2.15, 0.95)
```

If Eigen is installed in a nonstandard place, add it before compilation:

```cpp
gSystem->AddIncludePath("-I/path/to/eigen3");
.L RunRPR.C+
```

To export an angular table:

```bash
root -l -q 'ExportRPR.C(2.15,"rpr_amplitudes.tsv",201)'
```

## Select and adjust contributions

The default names are:

```text
u_lambda, u_sigma, u_sigma1385, s_lambda, s_sigma,
Lambda1890, Lambda2100, Sigma2030, Sigma2230, Sigma2250
```

Typical controls are:

```cpp
RPRModel model;

model.SetEnabled("Sigma2250", false);
model.EnableOnly({"u_lambda", "u_sigma", "Lambda2100", "Sigma2030"});

// Phenomenological multiplier of one complete diagram:
model.SetScale("Lambda2100", 0.80);

// A complex scale also permits an extra fitted phase:
model.SetScale("Sigma2030", Complex(0.90, 0.10));

// Change the two vertex couplings themselves:
model.SetCouplings("Lambda2100", 2.41, 2.95);
```

Only the product `gKN*gKXi` occurs in this reaction, so these two couplings
cannot be determined separately from `K- p -> K+ Xi-` data.  For a fit, the
single `amplitudeScale` is usually the clearer parameter.  A polarization-only
fit also cannot determine a common scale applied to every amplitude.

A hypothetical spin-1/2 s-channel resonance can be added with

```cpp
model.AddResonance(MakeSpinHalfResonance(
    "Lambda2000", "Lambda", -1, 2.000, 0.150));
model.SetScale("Lambda2000", 0.5);
```

The published selected set itself contains spin-3/2 and spin-7/2 excited
states; the spin-1/2 helper is provided for model tests or extensions.

## Extract f and g

With the reaction plane chosen as x-z and

```text
+y = k(K-) cross k(K+),
```

the code defines

```text
M = f I - i g sigma_y = [ f  -g ]
                         [ g   f ] .
```

Consequently,

```cpp
ModelResult r = model.Evaluate(sqrtS, cosTheta);
Complex f = r.amplitudes.f;
Complex g = r.amplitudes.g;

double ReF = f.real();
double ImF = f.imag();
double ReG = g.real();
double ImG = g.imag();
```

The reconstruction formulas are

```text
dSigma/dOmega = C(s) (|f|^2 + |g|^2),
P_Xi(+y)      = 2 Im(f* g) / (|f|^2 + |g|^2),
```

where `f*` denotes the complex conjugate and `C(s)` is the phase-space and
unit-conversion factor.  `ModelResult` also evaluates both observables directly
from the full 2x2 density matrix.  `parityResidual` and the two `...Matrix`
members allow the two calculations to be checked against each other.

The normal used for Fig. 7 of the Kim paper is opposite to the E42 convention,
so change the sign of the plotted polarization when comparing directly with
that figure.

## File roles

- `ComplexLib.hh`: complex Eigen matrices, gamma matrices, spinors, and
  contracted high-spin projectors.
- `MathLib.hh`: CM kinematics, form factors, Regge factor, and f/g extraction.
- `ResonanceLib.hh`: default contribution table and spin-1/2 helper.
- `RPR_Model.hh`: coherent amplitude, setters, and observables.
- `RunRPR.C`: single-point output and angular plots.
- `ExportRPR.C`: TSV export.
- `ValidateRPR.cpp`: comparison with the independent Python implementation.

## Standalone validation

Outside ROOT, compile `ValidateRPR.cpp` using any C++17 compiler:

```bash
g++ -std=c++17 -O2 -I/path/to/eigen3 ValidateRPR.cpp -o ValidateRPR
./ValidateRPR
```
