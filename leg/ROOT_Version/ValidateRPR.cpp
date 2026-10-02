#include "RPR_Model.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

int main()
{
    using namespace rpr;
    RPRModel model;
    const ModelResult r = model.Evaluate(2.15, 0.95);

    // Reference numbers from the independently implemented Python version.
    const Complex expectedF(-0.11488087566997916, -2.355210174219978);
    const Complex expectedG( 2.4345139868314307,  -3.1652338109773357);
    const double expectedCrossSection = 9.638284881115913;
    const double expectedPolarization = 0.5670492373744597;

    double error = 0.0;
    error = std::max(error, std::abs(r.amplitudes.f-expectedF));
    error = std::max(error, std::abs(r.amplitudes.g-expectedG));
    error = std::max(error, std::abs(r.differentialCrossSection-expectedCrossSection));
    error = std::max(error, std::abs(r.polarization-expectedPolarization));
    error = std::max(error, r.amplitudes.parityResidual);
    error = std::max(error, std::abs(r.differentialCrossSection
                                    -r.differentialCrossSectionMatrix));
    error = std::max(error, std::abs(r.polarization-r.polarizationMatrix));

    // Check that an individual amplitude scale is applied linearly.
    RPRModel scaled;
    scaled.EnableOnly({"Lambda2100"});
    const CMatrix2 original = scaled.Amplitude(2.15, 0.95);
    scaled.SetScale("Lambda2100", Complex(0.8, 0.1));
    const CMatrix2 changed = scaled.Amplitude(2.15, 0.95);
    error = std::max(error,
        (changed-Complex(0.8,0.1)*original).cwiseAbs().maxCoeff());

    // Check the optional finite-width spin-1/2 resonance path.
    RPRModel spinHalf;
    spinHalf.EnableOnly({});
    spinHalf.AddResonance(MakeSpinHalfResonance(
        "Lambda2000", "Lambda", -1, 2.000, 0.150));
    const ModelResult spinHalfResult = spinHalf.Evaluate(2.15, 0.95);
    if (!std::isfinite(spinHalfResult.differentialCrossSection)
        || !std::isfinite(spinHalfResult.polarization)) return 2;

    std::cout << "f = " << r.amplitudes.f << '\n'
              << "g = " << r.amplitudes.g << '\n'
              << "dSigma/dOmega = " << r.differentialCrossSection << '\n'
              << "P_Xi = " << r.polarization << '\n'
              << "maximum check error = " << error << '\n';
    return error < 1.0e-11 ? 0 : 1;
}
