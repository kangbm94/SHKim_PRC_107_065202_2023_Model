#include "FitE42.hh"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMinuit.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using rpr::Complex;
using rpr::ModelResult;
using rpr::RPRModel;

namespace {

constexpr double kScaleLower = -3.0;
constexpr double kScaleUpper =  3.0;
constexpr int kForwardFitIntervals = 16;

struct ResonanceBounds {
    const char* name;
    double scaleStart;
    double massStart;
    double massLow;
    double massHigh;
    double widthStart;
    double widthLow;
    double widthHigh;
};

// PDG 2025 Breit-Wigner ranges in GeV.  Sigma(2230) has only one
// listed solution, so its quoted value +/- uncertainty defines the bounds.
constexpr ResonanceBounds kSigma2030 = {
    "Sigma2030", 1.0, 2.030, 2.025, 2.040, 0.180, 0.150, 0.200
};
constexpr ResonanceBounds kLambda2100 = {
    "Lambda2100", 1.0, 2.100, 2.090, 2.110, 0.200, 0.100, 0.250
};
constexpr ResonanceBounds kSigma2250 = {
    "Sigma2250", 1.0, 2.250, 2.210, 2.280, 0.100, 0.060, 0.150
};
constexpr ResonanceBounds kSigma2230 = {
    "Sigma2230", 1.0, 2.240, 2.213, 2.267, 0.345, 0.295, 0.395
};

constexpr double kNodes[12] = {
    -0.9815606342467192, -0.9041172563704749,
    -0.7699026741943047, -0.5873179542866175,
    -0.3678314989981802, -0.1252334085114689,
     0.1252334085114689,  0.3678314989981802,
     0.5873179542866175,  0.7699026741943047,
     0.9041172563704749,  0.9815606342467192
};
constexpr double kWeights[12] = {
    0.0471753363865118, 0.1069393259953184,
    0.1600783285433462, 0.2031674267230659,
    0.2334925365383548, 0.2491470458134029,
    0.2491470458134029, 0.2334925365383548,
    0.2031674267230659, 0.1600783285433462,
    0.1069393259953184, 0.0471753363865118
};

std::vector<E42Point> gE42Data;
std::vector<ForwardPoint> gForwardData;
FitHypothesis gCurrentHypothesis;
bool gIncludeForwardCrossSection = true;

constexpr std::array<std::array<double,2>,3> kAngularBins = {{
    {{0.85, 0.90}}, {{0.90, 0.97}}, {{0.97, 1.00}}
}};

std::vector<std::string> SplitTabs(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

std::vector<E42Point> LoadE42Data(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open " + filename);
    std::vector<E42Point> result;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        E42Point point;
        double left = 0.0, right = 0.0;
        if (!(row >> point.sqrtS >> point.cosTheta >> point.polarization
                  >> point.error >> left >> right)) continue;
        point.cosThetaLow = point.cosTheta - left;
        point.cosThetaHigh = point.cosTheta + right;
        result.push_back(point);
    }
    if (result.empty()) throw std::runtime_error("No E42 points in " + filename);
    return result;
}

std::vector<ForwardPoint> LoadForwardData(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open " + filename);
    std::vector<ForwardPoint> result;
    std::string line;
    std::getline(input, line);
    constexpr double mK = 0.493677;
    constexpr double mP = 0.9382720813;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = SplitTabs(line);
        if (fields.size() < 4) continue;
        ForwardPoint point;
        point.experiment = fields[0];
        point.pLab = std::stod(fields[1]);
        point.crossSection = std::stod(fields[2]);
        point.error = std::stod(fields[3]);
        const double eK = std::sqrt(point.pLab*point.pLab + mK*mK);
        point.sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
        result.push_back(point);
    }
    return result;
}

std::vector<ResonanceBounds> ActiveResonances(const FitHypothesis& hypothesis)
{
    std::vector<ResonanceBounds> result = {
        kSigma2030, kLambda2100, kSigma2250
    };
    if (hypothesis.includeSigma2230) result.push_back(kSigma2230);
    return result;
}

void ConfigureModel(RPRModel& model, const FitHypothesis& hypothesis,
                    const double* parameters)
{
    const auto active = ActiveResonances(hypothesis);
    for (std::size_t i = 0; i < active.size(); ++i) {
        auto& resonance = model.Get(active[i].name);
        resonance.amplitudeScale = Complex(parameters[3*i], 0.0);
        resonance.mass = parameters[3*i + 1];
        resonance.width = parameters[3*i + 2];
    }
    auto& sigma2250 = model.Get("Sigma2250");
    sigma2250.spin = 0.5*hypothesis.sigma2250TwiceSpin;
    sigma2250.parity = hypothesis.sigma2250Parity;

    auto& sigma2230 = model.Get("Sigma2230");
    sigma2230.enabled = hypothesis.includeSigma2230;
    if (!hypothesis.includeSigma2230)
        sigma2230.amplitudeScale = Complex(0.0, 0.0);
}

void ConfigureModel(RPRModel& model, const HypothesisResult& result)
{
    std::vector<double> values;
    values.reserve(result.parameters.size());
    for (const auto& parameter : result.parameters)
        values.push_back(parameter.value);
    ConfigureModel(model, result.hypothesis, values.data());
}

double AveragePolarization(RPRModel& model, double sqrtS,
                           double cosThetaLow, double cosThetaHigh)
{
    const double middle = 0.5*(cosThetaLow + cosThetaHigh);
    const double halfWidth = 0.5*(cosThetaHigh - cosThetaLow);
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < 12; ++i) {
        const double cosTheta = middle + halfWidth*kNodes[i];
        const ModelResult value = model.Evaluate(sqrtS, cosTheta);
        numerator += kWeights[i]*value.polarization
                   * value.differentialCrossSection;
        denominator += kWeights[i]*value.differentialCrossSection;
    }
    if (!(denominator > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return numerator/denominator;
}

struct Chi2Contributions {
    double polarization = 0.0;
    double forwardCrossSection = 0.0;

    double Total() const
    {
        return polarization
             + (gIncludeForwardCrossSection ? forwardCrossSection : 0.0);
    }
};

Chi2Contributions CalculateChi2(const FitHypothesis& hypothesis,
                                const double* parameters)
{
    RPRModel model;
    ConfigureModel(model, hypothesis, parameters);
    Chi2Contributions chi2;
    for (const auto& point : gE42Data) {
        const double prediction = AveragePolarization(
            model, point.sqrtS, point.cosThetaLow, point.cosThetaHigh);
        if (!std::isfinite(prediction)) {
            chi2.polarization = 1.0e100;
            return chi2;
        }
        const double pull = (prediction - point.polarization)/point.error;
        chi2.polarization += pull*pull;
    }
    if (gIncludeForwardCrossSection) {
        for (const auto& point : gForwardData) {
            const double prediction = model.AverageDifferentialCrossSectionLabTheta(
                point.sqrtS, 0.0, 20.0, kForwardFitIntervals);
            if (!std::isfinite(prediction)) {
                chi2.forwardCrossSection = 1.0e100;
                return chi2;
            }
            const double pull = (prediction-point.crossSection)/point.error;
            chi2.forwardCrossSection += pull*pull;
        }
    }
    return chi2;
}
int nstep = 0;
void MinuitFCN(int&, double*, double& value, double* parameters, int)
{
    const Chi2Contributions chi2 =
        CalculateChi2(gCurrentHypothesis, parameters);

    nstep++;
    if(nstep %10 == 0){
        std::cout << "Step " << nstep << ": chi2 = " << chi2.Total() << std::endl;
    }
    // This is the objective seen by TMinuit.  Keep the two contributions
    // explicit here so a simultaneous fit cannot be confused with the older
    // polarization-only FCN.
    value = chi2.polarization;
    if (gIncludeForwardCrossSection)
        value += chi2.forwardCrossSection;
}

std::vector<double> MakeStart(const FitHypothesis& hypothesis, int iStart,
                              std::mt19937_64& random)
{
    const auto active = ActiveResonances(hypothesis);
    std::vector<double> result(3*active.size());
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (std::size_t i = 0; i < active.size(); ++i) {
        if (iStart == 0) {
            result[3*i] = active[i].scaleStart;
            result[3*i + 1] = active[i].massStart;
            result[3*i + 2] = active[i].widthStart;
        }
        else {
            result[3*i] = kScaleLower
                        + (kScaleUpper-kScaleLower)*unit(random);
            result[3*i + 1] = active[i].massLow
                            + (active[i].massHigh-active[i].massLow)*unit(random);
            result[3*i + 2] = active[i].widthLow
                            + (active[i].widthHigh-active[i].widthLow)*unit(random);
        }
    }
    return result;
}

HypothesisResult FitOneHypothesis(const FitHypothesis& hypothesis, int nStarts)
{
    gCurrentHypothesis = hypothesis;
    const auto active = ActiveResonances(hypothesis);
    const int nParameters = 3*static_cast<int>(active.size());
    HypothesisResult best;
    best.hypothesis = hypothesis;
    best.chi2 = std::numeric_limits<double>::infinity();
    best.numberOfPolarizationPoints = static_cast<int>(gE42Data.size());
    best.numberOfForwardCrossSectionPoints = gIncludeForwardCrossSection
        ? static_cast<int>(gForwardData.size()) : 0;
    best.numberOfPoints = best.numberOfPolarizationPoints
                        + best.numberOfForwardCrossSectionPoints;
    best.numberOfParameters = nParameters;
    best.includesForwardCrossSection = gIncludeForwardCrossSection;
    std::mt19937_64 random(20250915ULL
        + 100*hypothesis.sigma2250TwiceSpin
        + 10*(hypothesis.sigma2250Parity > 0)
        + hypothesis.includeSigma2230);

    for (int iStart = 0; iStart < std::max(1, nStarts); ++iStart) {
        nstep = 0;
        TMinuit minuit(nParameters);
        minuit.SetFCN(MinuitFCN);
        minuit.SetPrintLevel(-1);
        minuit.SetErrorDef(1.0);
        const auto start = MakeStart(hypothesis, iStart, random);
        int errorFlag = 0;
        for (std::size_t i = 0; i < active.size(); ++i) {
            const std::string prefix = active[i].name;
            minuit.mnparm(3*i, (prefix + ".scale").c_str(), start[3*i],
                          0.02, kScaleLower, kScaleUpper, errorFlag);
            minuit.mnparm(3*i + 1, (prefix + ".mass").c_str(), start[3*i + 1],
                          0.001, active[i].massLow, active[i].massHigh, errorFlag);
            minuit.mnparm(3*i + 2, (prefix + ".width").c_str(), start[3*i + 2],
                          0.002, active[i].widthLow, active[i].widthHigh, errorFlag);
        }
        double arguments[2] = {3000.0, 0.01};
        minuit.mnexcm("MIGRAD", arguments, 2, errorFlag);

        double amin = 0.0, edm = 0.0, errdef = 0.0;
        int nvpar = 0, nparx = 0, status = 0;
        minuit.mnstat(amin, edm, errdef, nvpar, nparx, status);
        if (!std::isfinite(amin) || amin >= best.chi2) continue;

        best.chi2 = amin;
        best.minimizerStatus = status;
        best.parameters.clear();
        for (std::size_t i = 0; i < active.size(); ++i) {
            for (int offset = 0; offset < 3; ++offset) {
                double value = 0.0, error = 0.0;
                minuit.GetParameter(3*i + offset, value, error);
                FittedParameter parameter;
                parameter.name = std::string(active[i].name)
                    + (offset == 0 ? ".scale"
                       : offset == 1 ? ".mass_GeV" : ".width_GeV");
                parameter.value = value;
                parameter.error = error;
                parameter.lower = offset == 0 ? kScaleLower
                    : offset == 1 ? active[i].massLow : active[i].widthLow;
                parameter.upper = offset == 0 ? kScaleUpper
                    : offset == 1 ? active[i].massHigh : active[i].widthHigh;
                best.parameters.push_back(parameter);
            }
        }

        std::vector<double> bestValues;
        bestValues.reserve(best.parameters.size());
        for (const auto& parameter : best.parameters)
            bestValues.push_back(parameter.value);
        const auto contributions = CalculateChi2(
            hypothesis, bestValues.data());
        best.polarizationChi2 = contributions.polarization;
        best.forwardCrossSectionChi2 = contributions.forwardCrossSection;
        best.chi2 = contributions.Total();
    }

    best.bic = best.chi2 + nParameters*std::log(best.numberOfPoints);
    for (std::size_t i = 0; i < best.parameters.size(); i += 3) {
        const auto& scale = best.parameters[i];
        if (std::abs(scale.value-scale.lower) < 0.01*(scale.upper-scale.lower)
            || std::abs(scale.upper-scale.value) < 0.01*(scale.upper-scale.lower))
            best.scaleAtBoundary = true;
    }
    return best;
}

int ColorForBin(double cosTheta)
{
    if (cosTheta < 0.90) return kBlue + 1;
    if (cosTheta < 0.97) return kRed + 1;
    return kGreen + 2;
}

struct AveragedAmplitudes {
    Complex f = Complex(0.0, 0.0);
    Complex g = Complex(0.0, 0.0);
    Complex fConjugateG = Complex(0.0, 0.0);
    double fIntensity = 0.0;
    double gIntensity = 0.0;
};

// Coherent bin averages are used for the amplitudes.  The intensity ratio and
// relative phase are formed from separately integrated bilinears, which avoids
// phase-wrap artefacts and is the direct angular-bin analogue of the
// polarization average.
AveragedAmplitudes AverageAmplitudes(RPRModel& model, double sqrtS,
                                     double cosThetaLow,
                                     double cosThetaHigh,
                                     const std::string* component = nullptr)
{
    const double middle = 0.5*(cosThetaLow + cosThetaHigh);
    const double halfWidth = 0.5*(cosThetaHigh - cosThetaLow);
    AveragedAmplitudes result;
    for (int i = 0; i < 12; ++i) {
        const double cosTheta = middle + halfWidth*kNodes[i];
        const ModelResult value = component
            ? model.Evaluate(sqrtS, cosTheta, *component)
            : model.Evaluate(sqrtS, cosTheta);
        const double normalizedWeight = 0.5*kWeights[i];
        result.f += normalizedWeight*value.amplitudes.f;
        result.g += normalizedWeight*value.amplitudes.g;
        result.fIntensity += normalizedWeight*std::norm(value.amplitudes.f);
        result.gIntensity += normalizedWeight*std::norm(value.amplitudes.g);
        result.fConjugateG += normalizedWeight
                           * std::conj(value.amplitudes.f)
                           * value.amplitudes.g;
    }
    return result;
}

double NonFlipRatio(const AveragedAmplitudes& amplitudes)
{
    const double denominator = amplitudes.fIntensity + amplitudes.gIntensity;
    return denominator > 0.0
        ? std::sqrt(amplitudes.fIntensity/denominator) : 0.0;
}

double RelativePhaseOverPi(const AveragedAmplitudes& amplitudes)
{
    return std::abs(amplitudes.fConjugateG) > 0.0
        ? std::arg(amplitudes.fConjugateG)/rpr::Pi : 0.0;
}

void DrawSpinObservables(RPRModel& model, const HypothesisResult& result)
{
    const std::string tag = result.hypothesis.FileLabel();
    auto* canvas = new TCanvas(("c_spin_" + tag).c_str(),
        ("Spin observables: " + result.hypothesis.Label()).c_str(), 1250, 520);
    canvas->Divide(2, 1);

    canvas->cd(1);
    gPad->SetLeftMargin(0.13);
    auto* ratioFrame = gPad->DrawFrame(1.98, 0.0, 2.29, 1.02);
    ratioFrame->SetTitle((result.hypothesis.Label()
        + ";#sqrt{s} [GeV];"
          "R_{non-flip} = #sqrt{#int|f|^{2}/#int(|f|^{2}+|g|^{2})}").c_str());
    auto* ratioLegend = new TLegend(0.14, 0.16, 0.55, 0.38);
    ratioLegend->SetBorderSize(0);
    ratioLegend->SetFillStyle(0);

    canvas->cd(2);
    gPad->SetLeftMargin(0.13);
    auto* phaseFrame = gPad->DrawFrame(1.98, -1.05, 2.29, 1.05);
    phaseFrame->SetTitle((result.hypothesis.Label()
        + ";#sqrt{s} [GeV];arg[#int f^{*}g dcos#theta]/#pi").c_str());
    auto* phaseLegend = new TLegend(0.14, 0.16, 0.55, 0.38);
    phaseLegend->SetBorderSize(0);
    phaseLegend->SetFillStyle(0);

    for (const auto& bin : kAngularBins) {
        const double centre = 0.5*(bin[0] + bin[1]);
        const int color = ColorForBin(centre);
        auto* ratio = new TGraph();
        auto* phase = new TGraph();
        for (double sqrtS = 1.99; sqrtS <= 2.281; sqrtS += 0.004) {
            const auto amplitudes = AverageAmplitudes(
                model, sqrtS, bin[0], bin[1]);
            ratio->SetPoint(ratio->GetN(), sqrtS, NonFlipRatio(amplitudes));
            phase->SetPoint(phase->GetN(), sqrtS,
                            RelativePhaseOverPi(amplitudes));
        }
        ratio->SetLineColor(color);
        ratio->SetLineWidth(2);
        phase->SetLineColor(color);
        phase->SetLineWidth(2);
        canvas->cd(1);
        ratio->Draw("L SAME");
        canvas->cd(2);
        phase->Draw("L SAME");
        const std::string entry = Form("%.2f < cos#theta < %.2f",
                                       bin[0], bin[1]);
        ratioLegend->AddEntry(ratio, entry.c_str(), "l");
        phaseLegend->AddEntry(phase, entry.c_str(), "l");
    }
    canvas->cd(1);
    ratioLegend->Draw();
    canvas->cd(2);
    auto* phaseZero = new TLine(1.98, 0.0, 2.29, 0.0);
    phaseZero->SetLineColor(kGray + 1);
    phaseZero->SetLineStyle(3);
    phaseZero->Draw();
    phaseLegend->Draw();

    canvas->SaveAs(("FitE42_spin_observables_" + tag + ".pdf").c_str());
    canvas->SaveAs(("FitE42_spin_observables_" + tag + ".png").c_str());
}

struct AmplitudeComponentStyle {
    std::string name;
    std::string label;
    int color = kBlack;
    bool fullSum = false;
};

std::vector<AmplitudeComponentStyle>
AmplitudeComponents(const FitHypothesis& hypothesis)
{
    std::vector<AmplitudeComponentStyle> result = {
        {"", "full coherent sum", kBlack, true},
        {"Sigma2250", "#Sigma(2250)", kMagenta + 1, false},
        {"Sigma2030", "#Sigma(2030)", kBlue + 1, false},
        {"Lambda2100", "#Lambda(2100)", kRed + 1, false}
    };
    if (hypothesis.includeSigma2230)
        result.push_back({"Sigma2230", "#Sigma(2230)", kGreen + 2, false});
    result.push_back({"Lambda1890", "#Lambda(1890)", kOrange + 7, false});
    result.push_back({"u_lambda", "Reggeized u-channel #Lambda",
                      kCyan + 2, false});
    return result;
}

struct AmplitudeGraphs {
    AmplitudeComponentStyle style;
    TGraph* fReal = nullptr;
    TGraph* fImaginary = nullptr;
    TGraph* gReal = nullptr;
    TGraph* gImaginary = nullptr;
};

void DrawAmplitudeComponents(RPRModel& model, const HypothesisResult& result,
                             double cosThetaLow, double cosThetaHigh)
{
    std::vector<AmplitudeGraphs> curves;
    double fMaximum = 0.0;
    double gMaximum = 0.0;
    for (const auto& style : AmplitudeComponents(result.hypothesis)) {
        AmplitudeGraphs graphs;
        graphs.style = style;
        graphs.fReal = new TGraph();
        graphs.fImaginary = new TGraph();
        graphs.gReal = new TGraph();
        graphs.gImaginary = new TGraph();
        for (double sqrtS = 1.99; sqrtS <= 2.281; sqrtS += 0.004) {
            const std::string* component = style.fullSum ? nullptr : &style.name;
            const auto amplitudes = AverageAmplitudes(
                model, sqrtS, cosThetaLow, cosThetaHigh, component);
            graphs.fReal->SetPoint(graphs.fReal->GetN(), sqrtS,
                                   amplitudes.f.real());
            graphs.fImaginary->SetPoint(graphs.fImaginary->GetN(), sqrtS,
                                        amplitudes.f.imag());
            graphs.gReal->SetPoint(graphs.gReal->GetN(), sqrtS,
                                   amplitudes.g.real());
            graphs.gImaginary->SetPoint(graphs.gImaginary->GetN(), sqrtS,
                                        amplitudes.g.imag());
            fMaximum = std::max(fMaximum, std::abs(amplitudes.f.real()));
            fMaximum = std::max(fMaximum, std::abs(amplitudes.f.imag()));
            gMaximum = std::max(gMaximum, std::abs(amplitudes.g.real()));
            gMaximum = std::max(gMaximum, std::abs(amplitudes.g.imag()));
        }
        for (auto* graph : {graphs.fReal, graphs.fImaginary,
                            graphs.gReal, graphs.gImaginary}) {
            graph->SetLineColor(style.color);
            graph->SetLineWidth(style.fullSum ? 3 : 2);
        }
        graphs.fImaginary->SetLineStyle(2);
        graphs.gImaginary->SetLineStyle(2);
        curves.push_back(graphs);
    }
    fMaximum = std::max(1.0e-6, 1.12*fMaximum);
    gMaximum = std::max(1.0e-6, 1.12*gMaximum);

    const int lowTag = static_cast<int>(std::lround(100.0*cosThetaLow));
    const int highTag = static_cast<int>(std::lround(100.0*cosThetaHigh));
    const std::string tag = Form("%s_cth_%02d_%02d",
        result.hypothesis.FileLabel().c_str(), lowTag, highTag);
    auto* canvas = new TCanvas(("c_amplitude_" + tag).c_str(),
        ("Amplitude components: " + result.hypothesis.Label()).c_str(),
        1350, 560);
    canvas->Divide(2, 1);

    canvas->cd(1);
    gPad->SetLeftMargin(0.13);
    auto* fFrame = gPad->DrawFrame(1.98, -fMaximum, 2.29, fMaximum);
    fFrame->SetTitle(Form("%.2f < cos#theta < %.2f;#sqrt{s} [GeV];"
                          "spin-non-flip amplitude f [model units]",
                          cosThetaLow, cosThetaHigh));
    for (auto& curve : curves) {
        curve.fReal->Draw("L SAME");
        curve.fImaginary->Draw("L SAME");
    }
    auto* fZero = new TLine(1.98, 0.0, 2.29, 0.0);
    fZero->SetLineColor(kGray + 1);
    fZero->SetLineStyle(3);
    fZero->Draw();
    auto* componentLegend = new TLegend(0.14, 0.58, 0.58, 0.89);
    componentLegend->SetBorderSize(0);
    componentLegend->SetFillStyle(0);
    for (auto& curve : curves)
        componentLegend->AddEntry(curve.fReal, curve.style.label.c_str(), "l");
    componentLegend->Draw();

    canvas->cd(2);
    gPad->SetLeftMargin(0.13);
    auto* gFrame = gPad->DrawFrame(1.98, -gMaximum, 2.29, gMaximum);
    gFrame->SetTitle(Form("%.2f < cos#theta < %.2f;#sqrt{s} [GeV];"
                          "spin-flip amplitude g [model units]",
                          cosThetaLow, cosThetaHigh));
    for (auto& curve : curves) {
        curve.gReal->Draw("L SAME");
        curve.gImaginary->Draw("L SAME");
    }
    auto* gZero = new TLine(1.98, 0.0, 2.29, 0.0);
    gZero->SetLineColor(kGray + 1);
    gZero->SetLineStyle(3);
    gZero->Draw();
    auto* styleLegend = new TLegend(0.68, 0.77, 0.88, 0.89);
    styleLegend->SetBorderSize(0);
    styleLegend->SetFillStyle(0);
    styleLegend->AddEntry(curves.front().gReal, "Re", "l");
    styleLegend->AddEntry(curves.front().gImaginary, "Im", "l");
    styleLegend->Draw();

    canvas->SaveAs(("FitE42_amplitudes_" + tag + ".pdf").c_str());
    canvas->SaveAs(("FitE42_amplitudes_" + tag + ".png").c_str());
}

void DrawHypothesis(const HypothesisResult& result,
                    const std::vector<ForwardPoint>& forwardData)
{
    RPRModel model;
    ConfigureModel(model, result);
    const std::string canvasName = "c_" + result.hypothesis.FileLabel();
    auto* canvas = new TCanvas(canvasName.c_str(),
        result.hypothesis.Label().c_str(), 1250, 520);
    canvas->Divide(2, 1);

    canvas->cd(1);
    gPad->SetLeftMargin(0.13);
    auto* frameP = gPad->DrawFrame(1.98, -1.05, 2.29, 0.45);
    frameP->SetTitle((result.hypothesis.Label()
        + ";#sqrt{s} [GeV];P_{#Xi}").c_str());
    auto* legendP = new TLegend(0.14, 0.65, 0.55, 0.89);
    legendP->SetBorderSize(0);
    legendP->SetFillStyle(0);

    for (const auto& bin : kAngularBins) {
        const double centre = 0.5*(bin[0] + bin[1]);
        const int color = ColorForBin(centre);
        auto* data = new TGraphErrors();
        auto* curve = new TGraph();
        for (const auto& point : gE42Data) {
            if (std::abs(point.cosTheta-centre) > 1.0e-6) continue;
            const int n = data->GetN();
            data->SetPoint(n, point.sqrtS, point.polarization);
            data->SetPointError(n, 0.0, point.error);
        }
        for (double sqrtS = 1.99; sqrtS <= 2.281; sqrtS += 0.004) {
            const double p = AveragePolarization(model, sqrtS,
                                                  bin[0], bin[1]);
            curve->SetPoint(curve->GetN(), sqrtS, p);
        }
        data->SetMarkerStyle(20);
        data->SetMarkerColor(color);
        data->SetLineColor(color);
        curve->SetLineColor(color);
        curve->SetLineWidth(2);
        curve->Draw("L SAME");
        data->Draw("P SAME");
        const std::string entry = Form("%.2f < cos#theta < %.2f",
                                       bin[0], bin[1]);
        legendP->AddEntry(data, entry.c_str(), "p");
    }
    auto* zero = new TLine(1.98, 0.0, 2.29, 0.0);
    zero->SetLineColor(kGray + 1);
    zero->SetLineStyle(3);
    zero->Draw();
    legendP->Draw();

    canvas->cd(2);
    gPad->SetLeftMargin(0.13);
    auto* frameX = gPad->DrawFrame(1.15, 0.0, 2.72, 80.0);
    const std::string forwardTitle =
        (result.includesForwardCrossSection
         ? "Forward cross-section simultaneous fit;"
         : "Forward cross-section validation;")
        + std::string("p_{K^{-}}^{lab} [GeV/c];"
                      "d#sigma/d#Omega_{lab} [#mub/sr]");
    frameX->SetTitle(forwardTitle.c_str());
    auto* dataX = new TGraphErrors();
    for (const auto& point : forwardData) {
        const int n = dataX->GetN();
        dataX->SetPoint(n, point.pLab, point.crossSection);
        dataX->SetPointError(n, 0.0, point.error);
    }
    dataX->SetMarkerStyle(20);
    dataX->SetMarkerSize(0.8);
    dataX->Draw("P SAME");
    auto* curveX = new TGraph();
    constexpr double mK = 0.493677;
    constexpr double mP = 0.9382720813;
    for (double pLab = 1.15; pLab <= 2.72; pLab += 0.015) {
        const double eK = std::sqrt(pLab*pLab + mK*mK);
        const double sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
        const double prediction = model.AverageDifferentialCrossSectionLabTheta(
            sqrtS, 0.0, 20.0);
        curveX->SetPoint(curveX->GetN(), pLab, prediction);
    }
    curveX->SetLineColor(kRed + 1);
    curveX->SetLineWidth(2);
    curveX->Draw("L SAME");
    auto* legendX = new TLegend(0.58, 0.75, 0.88, 0.88);
    legendX->SetBorderSize(0);
    legendX->SetFillStyle(0);
    legendX->AddEntry(dataX, "forward data", "p");
    legendX->AddEntry(curveX,
        result.includesForwardCrossSection
            ? "simultaneously fitted model"
            : "validation only (not in #chi^{2})", "l");
    legendX->Draw();

    canvas->SaveAs(("FitE42_" + result.hypothesis.FileLabel() + ".pdf").c_str());
    canvas->SaveAs(("FitE42_" + result.hypothesis.FileLabel() + ".png").c_str());

    DrawSpinObservables(model, result);
    for (const auto& bin : kAngularBins)
        DrawAmplitudeComponents(model, result, bin[0], bin[1]);
}

void WriteResults(const std::vector<HypothesisResult>& results)
{
    std::ofstream table("FitE42_hypothesis_results.tsv");
    table << "hypothesis\tchi2_total\tchi2_polarization"
             "\tchi2_forward_cross_section\tN_polarization"
             "\tN_forward_cross_section\tndof\tchi2_per_ndof\tnpar\tBIC"
             "\tminuit_status\tscale_at_boundary";
    const std::vector<std::string> parameterNames = {
        "Sigma2030.scale", "Sigma2030.mass_GeV", "Sigma2030.width_GeV",
        "Lambda2100.scale", "Lambda2100.mass_GeV", "Lambda2100.width_GeV",
        "Sigma2250.scale", "Sigma2250.mass_GeV", "Sigma2250.width_GeV",
        "Sigma2230.scale", "Sigma2230.mass_GeV", "Sigma2230.width_GeV"
    };
    for (const auto& name : parameterNames)
        table << '\t' << name << '\t' << name << "_error";
    table << '\n';
    table << std::setprecision(10);
    for (const auto& result : results) {
        const int ndof = result.numberOfPoints-result.numberOfParameters;
        table << result.hypothesis.Label() << '\t' << result.chi2 << '\t'
              << result.polarizationChi2 << '\t'
              << result.forwardCrossSectionChi2 << '\t'
              << result.numberOfPolarizationPoints << '\t'
              << result.numberOfForwardCrossSectionPoints << '\t'
              << ndof << '\t' << (ndof > 0 ? result.chi2/ndof : 0.0) << '\t'
              << result.numberOfParameters << '\t' << result.bic << '\t'
              << result.minimizerStatus << '\t' << result.scaleAtBoundary;
        std::map<std::string, FittedParameter> byName;
        for (const auto& parameter : result.parameters)
            byName[parameter.name] = parameter;
        for (const auto& name : parameterNames) {
            const auto found = byName.find(name);
            if (found == byName.end()) table << "\t\t";
            else table << '\t' << found->second.value
                       << '\t' << found->second.error;
        }
        table << '\n';
    }

    const auto best = std::min_element(results.begin(), results.end(),
        [](const HypothesisResult& a, const HypothesisResult& b) {
            return a.bic < b.bic;
        });
    std::ofstream summary("FitE42_best_fit.tsv");
    summary << "parameter\tvalue\terror\tlower\tupper\n";
    summary << "hypothesis\t" << best->hypothesis.Label() << "\t\t\t\n";
    summary << "chi2_total\t" << best->chi2 << "\t\t\t\n";
    summary << "chi2_polarization\t" << best->polarizationChi2 << "\t\t\t\n";
    summary << "chi2_forward_cross_section\t"
            << best->forwardCrossSectionChi2 << "\t\t\t\n";
    summary << "N_total\t" << best->numberOfPoints << "\t\t\t\n";
    summary << "BIC\t" << best->bic << "\t\t\t\n";
    summary << std::setprecision(10);
    for (const auto& parameter : best->parameters)
        summary << parameter.name << '\t' << parameter.value << '\t'
                << parameter.error << '\t' << parameter.lower << '\t'
                << parameter.upper << '\n';
}

void DrawComparison(const std::vector<HypothesisResult>& results)
{
    const int n = static_cast<int>(results.size());
    auto* canvas = new TCanvas("c_hypothesis_comparison",
                               "Hypothesis comparison", 1300, 650);
    canvas->SetBottomMargin(0.30);
    auto* chi2 = new TH1D("h_hypothesis_chi2",
        "Simultaneous-fit hypothesis comparison;hypothesis;#chi^{2} or BIC",
        n, 0, n);
    auto* chi2Polarization = new TH1D(
        "h_hypothesis_chi2_polarization", "", n, 0, n);
    auto* chi2Forward = new TH1D(
        "h_hypothesis_chi2_forward", "", n, 0, n);
    auto* bic = new TH1D("h_hypothesis_bic", "", n, 0, n);
    for (int i = 0; i < n; ++i) {
        chi2->SetBinContent(i+1, results[i].chi2);
        chi2Polarization->SetBinContent(i+1, results[i].polarizationChi2);
        chi2Forward->SetBinContent(i+1, results[i].forwardCrossSectionChi2);
        bic->SetBinContent(i+1, results[i].bic);
        chi2->GetXaxis()->SetBinLabel(i+1,
            results[i].hypothesis.Label().c_str());
    }
    chi2->SetFillColor(kAzure - 9);
    chi2->SetLineColor(kBlue + 1);
    chi2Polarization->SetMarkerStyle(22);
    chi2Polarization->SetMarkerColor(kGreen + 2);
    chi2Polarization->SetLineColor(kGreen + 2);
    chi2Forward->SetMarkerStyle(23);
    chi2Forward->SetMarkerColor(kMagenta + 1);
    chi2Forward->SetLineColor(kMagenta + 1);
    bic->SetMarkerStyle(20);
    bic->SetMarkerColor(kRed + 1);
    bic->SetLineColor(kRed + 1);
    chi2->GetXaxis()->LabelsOption("v");
    chi2->SetMaximum(1.20*std::max(chi2->GetMaximum(), bic->GetMaximum()));
    chi2->Draw("HIST");
    chi2Polarization->Draw("P SAME");
    if (results.front().includesForwardCrossSection)
        chi2Forward->Draw("P SAME");
    bic->Draw("P SAME");
    auto* legend = new TLegend(0.68, 0.70, 0.90, 0.90);
    legend->SetBorderSize(0);
    legend->AddEntry(chi2, "total #chi^{2}", "f");
    legend->AddEntry(chi2Polarization, "polarization #chi^{2}", "p");
    if (results.front().includesForwardCrossSection)
        legend->AddEntry(chi2Forward, "forward #chi^{2}", "p");
    legend->AddEntry(bic, "BIC", "p");
    legend->Draw();
    canvas->SaveAs("FitE42_hypothesis_comparison.pdf");
    canvas->SaveAs("FitE42_hypothesis_comparison.png");
}

} // namespace

std::string FitHypothesis::Label() const
{
    return Form("S2250 %d/2%c, S2230 %s", sigma2250TwiceSpin,
                sigma2250Parity > 0 ? '+' : '-',
                includeSigma2230 ? "on" : "off");
}

std::string FitHypothesis::FileLabel() const
{
    return Form("S2250_%dhalf_%s_S2230_%s", sigma2250TwiceSpin,
                sigma2250Parity > 0 ? "plus" : "minus",
                includeSigma2230 ? "on" : "off");
}

void FitE42(int nStarts, bool includeForwardCrossSection)
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("FitE42Plots", true);
    const std::string oldDirectory = gSystem->WorkingDirectory();

    gE42Data = LoadE42Data(
        "Datapoints/Polarization/E42_polarization_CM.tsv");
    gForwardData = LoadForwardData(
        "Datapoints/ForwardCrossSection/Forward_differential_cross_sections.tsv");
    gIncludeForwardCrossSection = includeForwardCrossSection;
    if (gE42Data.empty())
        throw std::runtime_error("The polarization dataset is empty");
    if (gIncludeForwardCrossSection && gForwardData.empty())
        throw std::runtime_error(
            "The simultaneous fit requires non-empty forward data");
    std::cout << "Fit mode: "
              << (gIncludeForwardCrossSection
                  ? "E42 polarization + forward differential cross section"
                  : "E42 polarization only") << std::endl;
    std::cout << "  N(polarization) = " << gE42Data.size()
              << ", N(forward) used in FCN = "
              << (gIncludeForwardCrossSection ? gForwardData.size() : 0)
              << std::endl;

    std::vector<HypothesisResult> results;
    for (int twiceSpin : {3, 5, 7}) {
        for (int parity : {-1, +1}) {
            for (bool include2230 : {false, true}) {
                FitHypothesis hypothesis;
                hypothesis.sigma2250TwiceSpin = twiceSpin;
                hypothesis.sigma2250Parity = parity;
                hypothesis.includeSigma2230 = include2230;
                std::cout << "Fitting " << hypothesis.Label() << std::endl;
                auto result = FitOneHypothesis(hypothesis, nStarts);
                std::cout << "  chi2(total) = " << result.chi2
                          << ", chi2(P) = " << result.polarizationChi2
                          << ", chi2(forward) = "
                          << result.forwardCrossSectionChi2
                          << ", BIC = " << result.bic
                          << ", status = " << result.minimizerStatus
                          << std::endl;
                results.push_back(result);
            }
        }
    }

    std::sort(results.begin(), results.end(),
        [](const HypothesisResult& a, const HypothesisResult& b) {
            return a.bic < b.bic;
        });
    WriteResults(results);
    DrawComparison(results);

    gSystem->ChangeDirectory("FitE42Plots");
    for (const auto& result : results) DrawHypothesis(result, gForwardData);
    gSystem->ChangeDirectory(oldDirectory.c_str());

    std::cout << "\nBest by BIC: " << results.front().hypothesis.Label()
              << "  chi2(total)=" << results.front().chi2
              << "  chi2(P)=" << results.front().polarizationChi2
              << "  chi2(forward)="
              << results.front().forwardCrossSectionChi2
              << "  BIC=" << results.front().bic << std::endl;
    for (const auto& parameter : results.front().parameters)
        std::cout << "  " << parameter.name << " = " << parameter.value
                  << " +/- " << parameter.error << std::endl;
    std::cout << "Results: FitE42_hypothesis_results.tsv, "
                 "FitE42_best_fit.tsv, FitE42Plots/" << std::endl;
}
