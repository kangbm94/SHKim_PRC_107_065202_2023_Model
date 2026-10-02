#pragma once

#include "RPRUtils.hh"

#include <string>
#include <vector>
std::vector<std::string> SplitTabs(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

struct E42Point {
    double sqrtS = 0.0;
    double cosTheta = 0.0;
    double polarization = 0.0;
    double error = 0.0;
    double cosThetaLow = 0.0;
    double cosThetaHigh = 0.0;
};

struct ForwardPoint {
    std::string experiment;
    double pLab = 0.0;
    double sqrtS = 0.0;
    double crossSection = 0.0;
    double error = 0.0;
};
struct DifferentialDataset {
    std::string experiment;
    double sqrtS = 0.0;
    std::vector<double> cosTheta;
    std::vector<double> crossSection;
    std::vector<double> error;
};

struct FitHypothesis {
    int sigma2250TwiceSpin = 7;
    int sigma2250Parity = -1;
    bool includeSigma2230 = false;

    std::string Label() const;
    std::string FileLabel() const;
};

struct FittedParameter {
    std::string name;
    double value = 0.0;
    double error = 0.0;
    double lower = 0.0;
    double upper = 0.0;
};

struct HypothesisResult {
    FitHypothesis hypothesis;
    std::vector<FittedParameter> parameters;
    double chi2 = 0.0;
    double polarizationChi2 = 0.0;
    double forwardCrossSectionChi2 = 0.0;
    double bic = 0.0;
    int numberOfPoints = 0;
    int numberOfPolarizationPoints = 0;
    int numberOfForwardCrossSectionPoints = 0;
    int numberOfParameters = 0;
    int minimizerStatus = 0;
    bool includesForwardCrossSection = true;
    bool scaleAtBoundary = false;
};

// Scan all supported Sigma(2250) J^P and Sigma(2230) on/off hypotheses.
// By default, simultaneously fit E42 polarization and the forward laboratory
// differential cross section. Pass false as the second argument to reproduce
// the earlier polarization-only fit. nWorkers=0 uses the available CPU count;
// nWorkers=1 runs sequentially.
