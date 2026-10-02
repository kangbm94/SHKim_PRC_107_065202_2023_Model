#include "RPR_Model.hh"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>

#include <iomanip>
#include <iostream>
#include <vector>

using namespace rpr;

void RunRPR(double sqrtS=2.15, double cosTheta=0.95,
            bool draw=true, bool includeRescattering=false)
{
    RPRModel model;
    //model.EnableRescattering(includeRescattering);
    //model.EnableRescattering();

    // Examples of adjustable model content.  Uncomment what you need.
    // model.SetEnabled("Sigma2250", false);
    // model.EnableOnly({"u_lambda", "u_sigma", "Lambda2100", "Sigma2030"});
    // model.SetScale("Lambda2100", 0.80);              // real amplitude scale
    // model.SetScale("Sigma2030", Complex(0.90,0.10)); // scale plus extra phase
    // model.SetCouplings("Lambda2100", 2.41, 2.95);
    // model.AddResonance(MakeSpinHalfResonance(
    //     "Lambda2000", "Lambda", -1, 2.000, 0.150));
    // model.EnableOnlyRescatterings({"resc_phi_lambda"});
    // model.SetRescatteringScale("resc_phi_lambda", Complex(0.8,0.1));

    const ModelResult value = model.Evaluate(sqrtS, cosTheta);
    const Complex f = value.amplitudes.f;
    const Complex g = value.amplitudes.g;

    std::cout << std::setprecision(10);
    std::cout << "sqrt(s) = " << sqrtS << " GeV, cos(theta) = "
              << cosTheta << "\n\n";
    std::cout << "M(up,up)     = " << value.spinMatrix(0,0) << '\n';
    std::cout << "M(up,down)   = " << value.spinMatrix(0,1) << '\n';
    std::cout << "M(down,up)   = " << value.spinMatrix(1,0) << '\n';
    std::cout << "M(down,down) = " << value.spinMatrix(1,1) << "\n\n";

    std::cout << "f = " << f.real() << " + i " << f.imag() << '\n';
    std::cout << "g = " << g.real() << " + i " << g.imag() << '\n';
    std::cout << "parity reconstruction residual = "
              << value.amplitudes.parityResidual << "\n\n";

    std::cout << "dSigma/dOmega = " << value.differentialCrossSection
              << " microbarn/sr\n";
    std::cout << "P_Xi(+y)      = " << value.polarization << '\n';
    std::cout << "matrix check  = (" << value.differentialCrossSectionMatrix
              << ", " << value.polarizationMatrix << ")\n";

    std::cout << "\nIndividual coherent components:\n";
    for (const auto& item : model.AmplitudeComponents(sqrtS, cosTheta)) {
        const FGAmplitudes fg = ExtractFG(item.second);
        std::cout << std::setw(14) << item.first
                  << "  f=" << std::setw(24) << fg.f
                  << "  g=" << std::setw(24) << fg.g << '\n';
    }

    if (!draw) return;

    constexpr int nPoints = 201;
    std::vector<double> c(nPoints), ds(nPoints), polarization(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        c[i] = -1.0 + 2.0*i/(nPoints-1.0);
        cout<<Form("Evaluating sqrtS = %g, cth = %g...", sqrtS, c[i])<<endl;
        const ModelResult point = model.Evaluate(sqrtS, c[i]);
        ds[i] = point.differentialCrossSection;
        polarization[i] = point.polarization;
        cout<<Form("Model: sqrtS = %g, cth = %g, dsig = %g, pol = %g", sqrtS, c[i], ds[i], polarization[i])<<endl;
    }

    auto canvas = new TCanvas("cRPR", "RPR model", 800, 800);
    canvas->Divide(1,2);
    canvas->cd(1);
    gPad->SetGrid();
    auto graphCross = new TGraph(nPoints, c.data(), ds.data());
    graphCross->SetTitle(Form("K^{-}p #rightarrow K^{+}#Xi^{-}, #sqrt{s}=%.3f GeV;cos#theta^{*}_{K^{+}};d#sigma/d#Omega [#mub/sr]", sqrtS));
    graphCross->SetLineWidth(3);
    graphCross->SetLineColor(kBlue+1);
    graphCross->Draw("AL");

    canvas->cd(2);
    gPad->SetGrid();
    auto graphPolarization = new TGraph(nPoints, c.data(), polarization.data());
    graphPolarization->SetTitle(";cos#theta^{*}_{K^{+}};P_{#Xi}(+y)");
    graphPolarization->SetLineWidth(3);
    graphPolarization->SetLineColor(kRed+1);
    graphPolarization->SetMinimum(-1.05);
    graphPolarization->SetMaximum(1.05);
    graphPolarization->Draw("AL");
    canvas->Update();
}
