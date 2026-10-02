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
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
// Structural smoke-test switch.  Keep this enabled only when checking the
// fit/IPC/output workflow; its toy objective does not fit the physics data.
bool gTestMode = 1;
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

struct DifferentialDataset {
    std::string experiment;
    double sqrtS = 0.0;
    std::vector<double> cosTheta;
    std::vector<double> crossSection;
    std::vector<double> error;
};

std::string SanitizedFilename(std::string value)
{
    for (char& character : value) {
        const unsigned char c = static_cast<unsigned char>(character);
        if (!std::isalnum(c) && character != '-' && character != '_')
            character = '_';
    }
    return value;
}

std::vector<DifferentialDataset> LoadDifferentialDatasets()
{
    const std::vector<std::pair<std::string,std::string>> sources = {
        {"Datapoints/DiffCrossSection/Berge_differential_cross_section.tsv",
         "Berge"},
        {"Datapoints/DiffCrossSection/Burgun_differential_cross_section.tsv",
         "Burgun"},
        {"Datapoints/DiffCrossSection/Dauber_differential_cross_section.tsv",
         "Dauber"},
        {"Datapoints/DiffCrossSection/Trippe_differential_cross_section.tsv",
         "Trippe"},
        {"Datapoints/DiffCrossSection/Baltay_differential_cross_section.tsv",
         "Baltay"},
        {"Datapoints/DiffCrossSection/E05_dif_CM.tsv", "J-PARC_E05"}
    };

    // The integer energy key prevents floating-point map-key ambiguities while
    // retaining sub-MeV separation if it exists in an input file.
    using Key = std::pair<std::string,long long>;
    std::map<Key,DifferentialDataset> grouped;
    for (const auto& source : sources) {
        std::ifstream input(source.first);
        if (!input)
            throw std::runtime_error("Cannot open " + source.first);
        std::string line;
        while (std::getline(input, line)) {
            const std::size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#') continue;
            std::istringstream row(line);
            double sqrtS = 0.0, cosTheta = 0.0;
            double crossSection = 0.0, error = 0.0;
            if (!(row >> sqrtS >> cosTheta >> crossSection >> error))
                throw std::runtime_error(
                    "Cannot parse differential-cross-section row in "
                    + source.first);
            const Key key{source.second,
                          static_cast<long long>(std::llround(1.0e6*sqrtS))};
            auto& dataset = grouped[key];
            dataset.experiment = source.second;
            dataset.sqrtS = sqrtS;
            dataset.cosTheta.push_back(cosTheta);
            dataset.crossSection.push_back(crossSection);
            dataset.error.push_back(error);
        }
    }

    std::vector<DifferentialDataset> result;
    result.reserve(grouped.size());
    for (auto& item : grouped) {
        auto& dataset = item.second;
        std::vector<std::size_t> order(dataset.cosTheta.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&dataset](std::size_t a,
                                                         std::size_t b) {
            return dataset.cosTheta[a] < dataset.cosTheta[b];
        });
        DifferentialDataset sorted;
        sorted.experiment = dataset.experiment;
        sorted.sqrtS = dataset.sqrtS;
        for (std::size_t i : order) {
            sorted.cosTheta.push_back(dataset.cosTheta[i]);
            sorted.crossSection.push_back(dataset.crossSection[i]);
            sorted.error.push_back(dataset.error[i]);
        }
        result.push_back(std::move(sorted));
    }
    std::sort(result.begin(), result.end(),
        [](const DifferentialDataset& a, const DifferentialDataset& b) {
            if (a.sqrtS != b.sqrtS) return a.sqrtS < b.sqrtS;
            return a.experiment < b.experiment;
        });
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
    int numberOfPoints = 0;

    double Total() const
    {
        return polarization
             + (gIncludeForwardCrossSection ? forwardCrossSection : 0.0);
    }
    double NumberOfPoints() const
    {
        return numberOfPoints;
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
        chi2.numberOfPoints++;
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
            chi2.numberOfPoints++;
        }
    }
    return chi2;
}



int ColorForBin(double cosTheta)
{
    if (cosTheta < 0.90) return kBlue + 1;
    if (cosTheta < 0.97) return kRed + 1;
    return kGreen + 2;
}


// Coherent bin averages are used for the amplitudes.  The intensity ratio and
// relative phase are formed from separately integrated bilinears, which avoids
// phase-wrap artefacts and is the direct angular-bin analogue of the
// polarization average.



struct AmplitudeComponentStyle {
    std::string name;
    std::string label;
    int color = kBlack;
    bool fullSum = false;
};
struct CrossSectionCurve {
    AmplitudeComponentStyle style;
    TGraph* graph = nullptr;
};
void StyleCrossSectionCurve(CrossSectionCurve& curve)
{
    curve.graph->SetLineColor(curve.style.color);
    curve.graph->SetLineWidth(curve.style.fullSum ? 3 : 2);
    curve.graph->SetLineStyle(curve.style.fullSum ? 1 : 2);
}



struct AmplitudeGraphs {
    AmplitudeComponentStyle style;
    TGraph* fReal = nullptr;
    TGraph* fImaginary = nullptr;
    TGraph* gReal = nullptr;
    TGraph* gImaginary = nullptr;
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

void DrawHypothesis(const FitHypothesis& hypothesis, double* parameters,
                    const std::vector<ForwardPoint>& forwardData,
                    const std::vector<DifferentialDataset>& differentialData)
{
    RPRModel model;
    ConfigureModel(model, hypothesis, parameters);
    const std::string canvasName = "c_" + hypothesis.FileLabel();
    auto* canvas = new TCanvas(canvasName.c_str(),
        hypothesis.Label().c_str(), 1250, 520);
    canvas->Divide(2, 1);

    canvas->cd(1);
    gPad->SetLeftMargin(0.13);
    auto* frameP = gPad->DrawFrame(1.9, -1.05, 2.4, 1.05);
    frameP->SetTitle((hypothesis.Label()
        + ";#sqrt{s} [GeV];P_{#Xi}").c_str());
    auto* legendP = new TLegend(0.14, 0.65, 0.55, 0.89);
    legendP->SetBorderSize(0);
    legendP->SetFillStyle(0);

    std::vector<TGraphErrors*> polarizationData;
    std::vector<TGraph*> polarizationCurves;
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
        polarizationData.push_back(data);
        polarizationCurves.push_back(curve);
    }
    auto* zero = new TLine(1.98, 0.0, 2.29, 0.0);
    zero->SetLineColor(kGray + 1);
    zero->SetLineStyle(3);
    zero->Draw();
    legendP->Draw();

    canvas->cd(2);
    gPad->SetLeftMargin(0.13);
    const std::string forwardTitle =
          ("Forward cross-section validation;")
        + std::string("p_{K^{-}}^{lab} [GeV/c];"
                      "d#sigma/d#Omega_{lab} [#mub/sr]");
    auto* dataX = new TGraphErrors();
    double forwardMaximum = 0.0;
    for (const auto& point : forwardData) {
        const int n = dataX->GetN();
        dataX->SetPoint(n, point.pLab, point.crossSection);
        dataX->SetPointError(n, 0.0, point.error);
        forwardMaximum = std::max(
            forwardMaximum, point.crossSection+point.error);
    }
    dataX->SetMarkerStyle(20);
    dataX->SetMarkerSize(0.8);
    std::vector<CrossSectionCurve> forwardCurves;
    for (const auto& style : AmplitudeComponents(hypothesis)) {
        CrossSectionCurve curve;
        curve.style = style;
        curve.graph = new TGraph();
        forwardCurves.push_back(curve);
    }
    constexpr double mK = 0.493677;
    constexpr double mP = 0.9382720813;
    for (double pLab = 1.15; pLab <= 2.72; pLab += 0.015) {
        const double eK = std::sqrt(pLab*pLab + mK*mK);
        const double sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
        for (auto& curve : forwardCurves) {
            const double prediction = curve.style.fullSum
                ? model.AverageDifferentialCrossSectionLabTheta(
                    sqrtS, 0.0, 20.0, kForwardFitIntervals)
                : model.AverageDifferentialCrossSectionLabTheta(
                    sqrtS, 0.0, 20.0, curve.style.name,
                    kForwardFitIntervals);
            curve.graph->SetPoint(curve.graph->GetN(), pLab, prediction);
            forwardMaximum = std::max(forwardMaximum, prediction);
        }
    }
    auto* frameX = gPad->DrawFrame(
        //1.15, 0.0, 2.72, std::max(1.0, 1.15*forwardMaximum));
        1., 0.0, 3.0, std::max(1.0, 1.15*forwardMaximum));
    frameX->SetTitle(forwardTitle.c_str());
    for (auto& curve : forwardCurves) {
        StyleCrossSectionCurve(curve);
        curve.graph->Draw("L SAME");
    }
    dataX->Draw("P SAME");
    auto* legendX = new TLegend(0.52, 0.48, 0.89, 0.89);
    legendX->SetBorderSize(0);
    legendX->SetFillStyle(0);
    legendX->AddEntry(dataX, "forward data", "p");
    legendX->AddEntry(forwardCurves.front().graph,
             "full model (not in #chi^{2})", "l");
    for (std::size_t i = 1; i < forwardCurves.size(); ++i)
        legendX->AddEntry(forwardCurves[i].graph,
                          forwardCurves[i].style.label.c_str(), "l");
    legendX->Draw();

    canvas->SaveAs(("FitE42_" + hypothesis.FileLabel() + ".pdf").c_str());

    delete legendP;
    delete legendX;
    delete zero;
    for (auto* graph : polarizationData) delete graph;
    for (auto* graph : polarizationCurves) delete graph;
    for (auto& curve : forwardCurves) delete curve.graph;
    delete dataX;
    delete canvas;

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

void ExtractChi2()
{
    gStyle->SetOptStat(0);
    TString figdir = "figs/Chi2Calc/";
    gSystem->mkdir(figdir, true);

    gE42Data = LoadE42Data(
        "Datapoints/Polarization/E42_polarization_CM.tsv");
    gForwardData = LoadForwardData(
        "Datapoints/ForwardCrossSection/Forward_differential_cross_sections.tsv");
    gIncludeForwardCrossSection = 1;
    if (gE42Data.empty())
        throw std::runtime_error("The polarization dataset is empty");

    //const auto differentialData = LoadDifferentialDatasets();
    vector<DifferentialDataset> differentialData;
    std::vector<FitHypothesis> hypotheses;
    RPRModel model;
    double parameters[12] = {
        //1., 2.04, 0.2,// Sigma(2030) scale, mass, width
        1., 2.03, 0.18,// Sigma(2030) scale, mass, width
        //1., 2.09, 0.15,// Lambda(2100) scale, mass, width
        1., 2.10, 0.20,// Lambda(2100) scale, mass, width
        1.,2.29, 0.1, // Sigma(2250) scale, mass, width
        1., 2.23, 0.345// Sigma(2230) scale, mass, width
    };
    double n_param = 5;// 5 were used in the fit: Sigma(2030) parameters, and Mass/width of Lambda(2100).
    for (int twiceSpin : { 7}) {
        for (int parity : {-1}) {
            for (bool include2230 : {true}) {
                FitHypothesis hypothesis;
                hypothesis.sigma2250TwiceSpin = twiceSpin;
                hypothesis.sigma2250Parity = parity;
                hypothesis.includeSigma2230 = include2230;
                auto Chi2 = CalculateChi2(hypothesis, parameters);
                cout<< "Hypothesis: " << hypothesis.Label()
                    << ", chi2 / ndf = " << Chi2.Total() << " / "
                    << Chi2.NumberOfPoints() - n_param
                    << " = " << Chi2.Total() / (Chi2.NumberOfPoints() - n_param) << endl;
                
                DrawHypothesis(hypothesis, parameters, gForwardData, differentialData);
            }
        }
    }

}
