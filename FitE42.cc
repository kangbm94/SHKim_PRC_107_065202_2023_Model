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

#if !defined(_WIN32)
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif


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

void MinuitFCN(int&, double*, double& value, double* parameters, int)
{
    if (gTestMode) {
        value = 0.0;
        for (int i = 0; i < 3; ++i)
            value += parameters[i]*parameters[i];
        return;
    }

    const Chi2Contributions chi2 =
        CalculateChi2(gCurrentHypothesis, parameters);

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

    if (best.parameters.size() != static_cast<std::size_t>(nParameters))
        throw std::runtime_error(
            "all Minuit starts failed for " + hypothesis.Label());

    best.bic = best.chi2 + nParameters*std::log(best.numberOfPoints);
    for (std::size_t i = 0; i < best.parameters.size(); i += 3) {
        const auto& scale = best.parameters[i];
        if (std::abs(scale.value-scale.lower) < 0.01*(scale.upper-scale.lower)
            || std::abs(scale.upper-scale.value) < 0.01*(scale.upper-scale.lower))
            best.scaleAtBoundary = true;
    }
    return best;
}

std::vector<HypothesisResult> FitHypothesesSequentially(
    const std::vector<FitHypothesis>& hypotheses, int nStarts)
{
    std::vector<HypothesisResult> results;
    results.reserve(hypotheses.size());
    for (const auto& hypothesis : hypotheses) {
        std::cout << "Fitting " << hypothesis.Label() << std::endl;
        auto result = FitOneHypothesis(hypothesis, nStarts);
        std::cout << "  chi2(total) = " << result.chi2
                  << ", chi2(P) = " << result.polarizationChi2
                  << ", chi2(forward) = "
                  << result.forwardCrossSectionChi2
                  << ", BIC = " << result.bic
                  << ", status = " << result.minimizerStatus
                  << std::endl;
        results.push_back(std::move(result));
    }
    return results;
}

#if !defined(_WIN32)

constexpr int kMaximumFitParameters = 12;

// Fixed-size packet used only for parent/worker communication. ROOT objects,
// STL strings, plotting, and file output remain in the parent process.
struct HypothesisPacket {
    std::int32_t success = 0;
    std::int32_t minimizerStatus = 0;
    std::int32_t parameterCount = 0;
    std::int32_t numberOfPoints = 0;
    std::int32_t numberOfPolarizationPoints = 0;
    std::int32_t numberOfForwardPoints = 0;
    std::int32_t includesForward = 0;
    std::int32_t scaleAtBoundary = 0;
    double chi2 = 0.0;
    double polarizationChi2 = 0.0;
    double forwardChi2 = 0.0;
    double bic = 0.0;
    double values[kMaximumFitParameters] = {};
    double errors[kMaximumFitParameters] = {};
};

bool WriteAll(int descriptor, const void* buffer, std::size_t bytes)
{
    const auto* current = static_cast<const unsigned char*>(buffer);
    while (bytes > 0) {
        const ssize_t written = ::write(descriptor, current, bytes);
        if (written < 0 && errno == EINTR) continue;
        if (written <= 0) return false;
        current += written;
        bytes -= static_cast<std::size_t>(written);
    }
    return true;
}

bool ReadAll(int descriptor, void* buffer, std::size_t bytes)
{
    auto* current = static_cast<unsigned char*>(buffer);
    while (bytes > 0) {
        const ssize_t count = ::read(descriptor, current, bytes);
        if (count < 0 && errno == EINTR) continue;
        if (count <= 0) return false;
        current += count;
        bytes -= static_cast<std::size_t>(count);
    }
    return true;
}

HypothesisPacket MakePacket(const HypothesisResult& result)
{
    HypothesisPacket packet;
    if (result.parameters.size() > kMaximumFitParameters) return packet;
    packet.success = 1;
    packet.minimizerStatus = result.minimizerStatus;
    packet.parameterCount = static_cast<std::int32_t>(result.parameters.size());
    packet.numberOfPoints = result.numberOfPoints;
    packet.numberOfPolarizationPoints = result.numberOfPolarizationPoints;
    packet.numberOfForwardPoints = result.numberOfForwardCrossSectionPoints;
    packet.includesForward = result.includesForwardCrossSection;
    packet.scaleAtBoundary = result.scaleAtBoundary;
    packet.chi2 = result.chi2;
    packet.polarizationChi2 = result.polarizationChi2;
    packet.forwardChi2 = result.forwardCrossSectionChi2;
    packet.bic = result.bic;
    for (std::size_t i = 0; i < result.parameters.size(); ++i) {
        packet.values[i] = result.parameters[i].value;
        packet.errors[i] = result.parameters[i].error;
    }
    return packet;
}

HypothesisResult ResultFromPacket(const FitHypothesis& hypothesis,
                                  const HypothesisPacket& packet)
{
    if (!packet.success)
        throw std::runtime_error("parallel fit failed for " + hypothesis.Label());
    const auto active = ActiveResonances(hypothesis);
    if (packet.parameterCount != 3*static_cast<int>(active.size()))
        throw std::runtime_error(
            "invalid worker result for " + hypothesis.Label());

    HypothesisResult result;
    result.hypothesis = hypothesis;
    result.minimizerStatus = packet.minimizerStatus;
    result.numberOfPoints = packet.numberOfPoints;
    result.numberOfPolarizationPoints = packet.numberOfPolarizationPoints;
    result.numberOfForwardCrossSectionPoints = packet.numberOfForwardPoints;
    result.numberOfParameters = packet.parameterCount;
    result.includesForwardCrossSection = packet.includesForward;
    result.scaleAtBoundary = packet.scaleAtBoundary;
    result.chi2 = packet.chi2;
    result.polarizationChi2 = packet.polarizationChi2;
    result.forwardCrossSectionChi2 = packet.forwardChi2;
    result.bic = packet.bic;

    int index = 0;
    for (const auto& resonance : active) {
        for (int offset = 0; offset < 3; ++offset, ++index) {
            FittedParameter parameter;
            parameter.name = std::string(resonance.name)
                + (offset == 0 ? ".scale"
                   : offset == 1 ? ".mass_GeV" : ".width_GeV");
            parameter.value = packet.values[index];
            parameter.error = packet.errors[index];
            parameter.lower = offset == 0 ? kScaleLower
                : offset == 1 ? resonance.massLow : resonance.widthLow;
            parameter.upper = offset == 0 ? kScaleUpper
                : offset == 1 ? resonance.massHigh : resonance.widthHigh;
            result.parameters.push_back(parameter);
        }
    }
    return result;
}

struct ActiveWorker {
    pid_t process = -1;
    int readDescriptor = -1;
    std::size_t hypothesisIndex = 0;
};

int AvailableWorkerCount(int requested, int numberOfHypotheses)
{
    if (requested == 1) return 1;
    long available = requested > 0 ? requested : ::sysconf(_SC_NPROCESSORS_ONLN);
    if (available < 1) available = 1;
    return std::min(numberOfHypotheses, static_cast<int>(available));
}

std::vector<HypothesisResult> FitHypothesesInParallel(
    const std::vector<FitHypothesis>& hypotheses, int nStarts, int nWorkers)
{
    const int workers = AvailableWorkerCount(
        nWorkers, static_cast<int>(hypotheses.size()));
    if (workers <= 1)
        return FitHypothesesSequentially(hypotheses, nStarts);

    std::cout << "Running " << hypotheses.size() << " hypotheses with "
              << workers << " isolated worker processes" << std::endl;
    std::vector<HypothesisResult> results(hypotheses.size());
    std::vector<bool> completed(hypotheses.size(), false);
    std::vector<ActiveWorker> active;
    std::size_t next = 0;

    while (next < hypotheses.size() || !active.empty()) {
        while (next < hypotheses.size()
               && static_cast<int>(active.size()) < workers) {
            int descriptors[2] = {-1, -1};
            if (::pipe(descriptors) != 0)
                throw std::runtime_error("cannot create parallel-fit pipe");
            const pid_t process = ::fork();
            if (process < 0) {
                ::close(descriptors[0]);
                ::close(descriptors[1]);
                throw std::runtime_error("cannot fork parallel-fit worker");
            }
            if (process == 0) {
                ::close(descriptors[0]);
                HypothesisPacket packet;
                try {
                    packet = MakePacket(
                        FitOneHypothesis(hypotheses[next], nStarts));
                }
                catch (...) {
                    packet.success = 0;
                }
                const bool written = WriteAll(
                    descriptors[1], &packet, sizeof(packet));
                ::close(descriptors[1]);
                ::_exit(written && packet.success ? 0 : 2);
            }

            ::close(descriptors[1]);
            active.push_back({process, descriptors[0], next});
            std::cout << "  launched " << hypotheses[next].Label()
                      << " (pid " << process << ")" << std::endl;
            ++next;
        }

        int status = 0;
        const pid_t finished = ::waitpid(-1, &status, 0);
        if (finished < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("waitpid failed during parallel fit");
        }
        const auto worker = std::find_if(active.begin(), active.end(),
            [finished](const ActiveWorker& value) {
                return value.process == finished;
            });
        if (worker == active.end())
            throw std::runtime_error("unknown parallel-fit worker completed");

        HypothesisPacket packet;
        const bool read = ReadAll(
            worker->readDescriptor, &packet, sizeof(packet));
        ::close(worker->readDescriptor);
        const std::size_t index = worker->hypothesisIndex;
        active.erase(worker);
        if (!read || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
            throw std::runtime_error(
                "worker process failed for " + hypotheses[index].Label());
        results[index] = ResultFromPacket(hypotheses[index], packet);
        completed[index] = true;
        const auto& result = results[index];
        std::cout << "  completed " << result.hypothesis.Label()
                  << ": chi2(total)=" << result.chi2
                  << ", chi2(P)=" << result.polarizationChi2
                  << ", chi2(forward)=" << result.forwardCrossSectionChi2
                  << ", BIC=" << result.bic
                  << ", status=" << result.minimizerStatus << std::endl;
    }

    if (std::find(completed.begin(), completed.end(), false) != completed.end())
        throw std::runtime_error("not all parallel hypotheses completed");
    return results;
}

#endif

std::vector<HypothesisResult> FitAllHypotheses(
    const std::vector<FitHypothesis>& hypotheses, int nStarts, int nWorkers)
{
#if defined(_WIN32)
    if (nWorkers != 1)
        std::cout << "Parallel hypothesis fitting is unavailable on Windows; "
                     "running sequentially." << std::endl;
    return FitHypothesesSequentially(hypotheses, nStarts);
#else
    return FitHypothesesInParallel(hypotheses, nStarts, nWorkers);
#endif
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

    std::vector<TGraph*> ratioGraphs;
    std::vector<TGraph*> phaseGraphs;
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
        ratioGraphs.push_back(ratio);
        phaseGraphs.push_back(phase);
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

    delete ratioLegend;
    delete phaseLegend;
    delete phaseZero;
    for (auto* graph : ratioGraphs) delete graph;
    for (auto* graph : phaseGraphs) delete graph;
    delete canvas;
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

void DrawDifferentialDataset(RPRModel& model,
                             const HypothesisResult& result,
                             const DifferentialDataset& dataset,
                             const std::string& outputDirectory)
{
    auto* data = new TGraphErrors();
    double maximum = 0.0;
    double insetMaximum = 0.0;
    for (std::size_t i = 0; i < dataset.cosTheta.size(); ++i) {
        data->SetPoint(data->GetN(), dataset.cosTheta[i],
                       dataset.crossSection[i]);
        data->SetPointError(data->GetN()-1, 0.0, dataset.error[i]);
        maximum = std::max(maximum,
                           dataset.crossSection[i]+dataset.error[i]);
        if (dataset.cosTheta[i] >= 0.5)
            insetMaximum = std::max(insetMaximum,
                dataset.crossSection[i]+dataset.error[i]);
    }
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
    data->SetMarkerColor(kBlack);
    data->SetLineColor(kBlack);

    std::vector<CrossSectionCurve> curves;
    for (const auto& style : AmplitudeComponents(result.hypothesis)) {
        CrossSectionCurve curve;
        curve.style = style;
        curve.graph = new TGraph();
        // Use an integer grid so floating-point accumulation can never send
        // cos(theta) slightly outside the model's required [-1,1] domain.
        for (int iCosTheta = 0; iCosTheta <= 200; ++iCosTheta) {
            const double cosTheta = -1.0 + 0.01*iCosTheta;
            const double prediction = style.fullSum
                ? model.DifferentialCrossSectionCM(dataset.sqrtS, cosTheta)
                : model.DifferentialCrossSectionCM(
                    dataset.sqrtS, cosTheta, style.name);
            curve.graph->SetPoint(curve.graph->GetN(), cosTheta, prediction);
            maximum = std::max(maximum, prediction);
            if (cosTheta >= 0.5)
                insetMaximum = std::max(insetMaximum, prediction);
        }
        StyleCrossSectionCurve(curve);
        curves.push_back(curve);
    }
    maximum = std::max(1.0, 1.15*maximum);
    insetMaximum = std::max(1.0, 1.15*insetMaximum);

    const int energyMeV = static_cast<int>(std::lround(1000.0*dataset.sqrtS));
    const std::string dataTag = SanitizedFilename(dataset.experiment)
        + Form("_sqrtS_%04dMeV", energyMeV);
    auto* canvas = new TCanvas(("c_differential_" + dataTag).c_str(),
        ("Differential cross section: " + dataTag).c_str(), 900, 700);
    canvas->SetLeftMargin(0.12);
    auto* frame = canvas->DrawFrame(-1.0, 0.0, 1.0, maximum);
    frame->SetTitle(Form("%s, #sqrt{s} = %.3f GeV;cos#theta_{CM};"
                         "d#sigma/d#Omega_{CM} [#mub/sr]",
                         dataset.experiment.c_str(), dataset.sqrtS));
    for (auto& curve : curves) curve.graph->Draw("L SAME");
    data->Draw("P SAME");

    auto* legend = new TLegend(0.13, 0.55, 0.48, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(data, dataset.experiment.c_str(), "p");
    for (auto& curve : curves)
        legend->AddEntry(curve.graph, curve.style.label.c_str(), "l");
    legend->Draw();

    auto* inset = new TPad(("inset_" + dataTag).c_str(), "", 0.53, 0.48,
                           0.94, 0.91);
    inset->SetFillStyle(0);
    inset->SetLeftMargin(0.15);
    inset->SetBottomMargin(0.16);
    inset->Draw();
    inset->cd();
    auto* insetFrame = inset->DrawFrame(0.5, 0.0, 1.0, insetMaximum);
    insetFrame->SetTitle("forward zoom;cos#theta_{CM};d#sigma/d#Omega_{CM}");
    for (auto& curve : curves) curve.graph->Draw("L SAME");
    data->Draw("P SAME");

    canvas->cd(0);
    const std::string outputStem = outputDirectory
        + "/FitE42_differential_" + dataTag;
    canvas->SaveAs((outputStem + ".pdf").c_str());
    canvas->SaveAs((outputStem + ".png").c_str());

    // Hundreds of monitoring canvases are produced in a full hypothesis scan;
    // release them immediately instead of retaining every canvas in ROOT.
    delete legend;
    for (auto& curve : curves) delete curve.graph;
    delete data;
    delete canvas;
}

void DrawDifferentialMonitoring(
    RPRModel& model, const HypothesisResult& result,
    const std::vector<DifferentialDataset>& differentialData)
{
    const std::string directory =
        "Differential/" + result.hypothesis.FileLabel();
    gSystem->mkdir(directory.c_str(), true);
    if (gSystem->AccessPathName(directory.c_str()))
        throw std::runtime_error(
            "Cannot create differential-plot directory " + directory);

    int failures = 0;
    for (const auto& dataset : differentialData) {
        try {
            DrawDifferentialDataset(
                model, result, dataset, directory);
        }
        catch (const std::exception& error) {
            ++failures;
            std::cerr << "WARNING: differential plot failed for "
                      << result.hypothesis.Label() << ", "
                      << dataset.experiment << " at sqrt(s)="
                      << dataset.sqrtS << " GeV: " << error.what()
                      << std::endl;
        }
    }
    if (failures != 0)
        std::cerr << "WARNING: " << failures << " of "
                  << differentialData.size()
                  << " differential plots failed for "
                  << result.hypothesis.Label() << std::endl;
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

    delete componentLegend;
    delete styleLegend;
    delete fZero;
    delete gZero;
    for (auto& curve : curves) {
        delete curve.fReal;
        delete curve.fImaginary;
        delete curve.gReal;
        delete curve.gImaginary;
    }
    delete canvas;
}

void DrawHypothesis(const HypothesisResult& result,
                    const std::vector<ForwardPoint>& forwardData,
                    const std::vector<DifferentialDataset>& differentialData)
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
        (result.includesForwardCrossSection
         ? "Forward cross-section simultaneous fit;"
         : "Forward cross-section validation;")
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
    for (const auto& style : AmplitudeComponents(result.hypothesis)) {
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
        1.15, 0.0, 2.72, std::max(1.0, 1.15*forwardMaximum));
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
        result.includesForwardCrossSection
            ? "full simultaneous fit"
            : "full model (not in #chi^{2})", "l");
    for (std::size_t i = 1; i < forwardCurves.size(); ++i)
        legendX->AddEntry(forwardCurves[i].graph,
                          forwardCurves[i].style.label.c_str(), "l");
    legendX->Draw();

    canvas->SaveAs(("FitE42_" + result.hypothesis.FileLabel() + ".pdf").c_str());
    canvas->SaveAs(("FitE42_" + result.hypothesis.FileLabel() + ".png").c_str());

    delete legendP;
    delete legendX;
    delete zero;
    for (auto* graph : polarizationData) delete graph;
    for (auto* graph : polarizationCurves) delete graph;
    for (auto& curve : forwardCurves) delete curve.graph;
    delete dataX;
    delete canvas;

    DrawSpinObservables(model, result);
    for (const auto& bin : kAngularBins)
        DrawAmplitudeComponents(model, result, bin[0], bin[1]);
    DrawDifferentialMonitoring(model, result, differentialData);
}

void WriteResults(const std::vector<HypothesisResult>& results)
{
    gSystem->mkdir("FitE42Parameters", true);
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

        const std::string parameterFilename =
            "FitE42Parameters/FitE42_parameters_"
            + result.hypothesis.FileLabel() + ".tsv";
        std::ofstream parameters(parameterFilename);
        if (!parameters)
            throw std::runtime_error("Cannot create " + parameterFilename);
        parameters << "quantity\tvalue\terror\tlower\tupper\n";
        parameters << "hypothesis\t" << result.hypothesis.Label()
                   << "\t\t\t\n";
        parameters << std::setprecision(10);
        parameters << "chi2_total\t" << result.chi2 << "\t\t\t\n";
        parameters << "chi2_polarization\t" << result.polarizationChi2
                   << "\t\t\t\n";
        parameters << "chi2_forward_cross_section\t"
                   << result.forwardCrossSectionChi2 << "\t\t\t\n";
        parameters << "BIC\t" << result.bic << "\t\t\t\n";
        parameters << "N_total\t" << result.numberOfPoints << "\t\t\t\n";
        parameters << "N_parameters\t" << result.numberOfParameters
                   << "\t\t\t\n";
        parameters << "minuit_status\t" << result.minimizerStatus
                   << "\t\t\t\n";
        for (const auto& parameter : result.parameters)
            parameters << parameter.name << '\t' << parameter.value << '\t'
                       << parameter.error << '\t' << parameter.lower << '\t'
                       << parameter.upper << '\n';
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

    delete legend;
    delete chi2;
    delete chi2Polarization;
    delete chi2Forward;
    delete bic;
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

void FitE42(int nStarts, bool includeForwardCrossSection, int nWorkers)
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

    std::vector<FitHypothesis> hypotheses;
    for (int twiceSpin : {3, 5, 7}) {
        for (int parity : {-1, +1}) {
            for (bool include2230 : {false, true}) {
                FitHypothesis hypothesis;
                hypothesis.sigma2250TwiceSpin = twiceSpin;
                hypothesis.sigma2250Parity = parity;
                hypothesis.includeSigma2230 = include2230;
                hypotheses.push_back(hypothesis);
            }
        }
    }
    std::vector<HypothesisResult> results =
        FitAllHypotheses(hypotheses, nStarts, nWorkers);

    std::sort(results.begin(), results.end(),
        [](const HypothesisResult& a, const HypothesisResult& b) {
            return a.bic < b.bic;
        });
    WriteResults(results);
    DrawComparison(results);

    const auto differentialData = LoadDifferentialDatasets();
    std::cout << "Loaded " << differentialData.size()
              << " differential-cross-section energy datasets for monitoring"
              << std::endl;

    gSystem->ChangeDirectory("FitE42Plots");
    int hypothesisPlotFailures = 0;
    for (const auto& result : results) {
        std::cout << "Plotting " << result.hypothesis.Label() << std::endl;
        try {
            DrawHypothesis(result, gForwardData, differentialData);
        }
        catch (const std::exception& error) {
            ++hypothesisPlotFailures;
            std::cerr << "WARNING: plotting failed for "
                      << result.hypothesis.Label() << ": "
                      << error.what() << std::endl;
        }
    }
    gSystem->ChangeDirectory(oldDirectory.c_str());
    if (hypothesisPlotFailures != 0)
        std::cerr << "WARNING: plotting failed for "
                  << hypothesisPlotFailures << " of " << results.size()
                  << " hypotheses; completed plots were retained."
                  << std::endl;

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
