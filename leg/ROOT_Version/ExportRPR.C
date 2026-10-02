#include "RPR_Model.hh"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace rpr;

// Export f, g, differential cross section, and polarization at fixed sqrt(s).
// Example:
//   root -l -q 'ExportRPR.C(2.15,"rpr_amplitudes.tsv",201)'
void ExportRPR(double sqrtS=2.15,
               const char* outputName="rpr_amplitudes.tsv",
               int nPoints=201)
{
    if (nPoints < 2) nPoints = 2;
    RPRModel model;

    // Apply the same choices here as in RunRPR.C, for example:
    // model.SetEnabled("Sigma2250", false);
    // model.SetScale("Lambda2100", 0.8);

    std::ofstream output(outputName);
    if (!output) {
        std::cerr << "Cannot open " << outputName << '\n';
        return;
    }

    output << "#sqrt_s\tcos_theta\tRe_f\tIm_f\tRe_g\tIm_g"
              "\tdsigma_domega_ub_sr\tP_Xi\n";
    output << std::setprecision(12);
    for (int i = 0; i < nPoints; ++i) {
        const double cosTheta = -1.0+2.0*i/(nPoints-1.0);
        const ModelResult r = model.Evaluate(sqrtS, cosTheta);
        output << sqrtS << '\t' << cosTheta << '\t'
               << r.amplitudes.f.real() << '\t' << r.amplitudes.f.imag() << '\t'
               << r.amplitudes.g.real() << '\t' << r.amplitudes.g.imag() << '\t'
               << r.differentialCrossSection << '\t' << r.polarization << '\n';
    }
    std::cout << "Wrote " << outputName << '\n';
}
