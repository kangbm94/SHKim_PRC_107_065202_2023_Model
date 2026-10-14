#pragma once

#include "RPR_Model.hh"
#include "MathLib.hh"
inline Amplitudes AverageAmplitudes(
    RPRModel& model, double sqrtS,
    double cosThetaLow, double cosThetaHigh,
    const std::string* component = nullptr)
{
    Amplitudes result;
    double width = cosThetaHigh - cosThetaLow;
    double dcth = 0.001;
    for (double cth_0 = cosThetaLow; cth_0 < cosThetaHigh; cth_0 += dcth) {
        //double cth = cth_0 + dcth < cosThetaHigh ? cth_0 + 0.5 * dcth : 0.5 * (cth_0 + cosThetaHigh);
        //double dcth_width = 2 *(cth - cth_0);
        double upper = std::min(cth_0 + dcth, cosThetaHigh);
        double dcth_width = upper - cth_0;
        double cth = cth_0 + 0.5 * dcth_width;
        double weight = dcth_width / width;
        const ModelResult value = component
            ? model.Evaluate(sqrtS, cth, *component)
            : model.Evaluate(sqrtS, cth);
        result.f += weight*value.amplitudes.f;
        result.g += weight*value.amplitudes.g;
        result.MatrixElement += weight*value.spinMatrix;
        result.fIntensity += weight*std::norm(value.amplitudes.f);
        result.gIntensity += weight*std::norm(value.amplitudes.g);
        result.fConjugateG += weight
                   * std::conj(value.amplitudes.f)
                   * value.amplitudes.g;
        //result.fIntensity += weight*std::norm(value.amplitudes.f);
        //result.gIntensity += weight*std::norm(value.amplitudes.g);
        //result.fConjugateG += weight
        //                   * std::conj(value.amplitudes.f)
        //                   * value.amplitudes.g;
    }
    return result;
}
Amplitudes PointAmplitude(RPRModel& model,
double sqrtS, double cosTheta, const std::string* component = nullptr){
    Amplitudes result;
    const ModelResult value = component?
        model.Evaluate(sqrtS, cosTheta, *component)
        : model.Evaluate(sqrtS, cosTheta);
    result.f = value.amplitudes.f;
    result.g = value.amplitudes.g;
    result.MatrixElement = value.amplitudes.MatrixElement;
    result.fIntensity = value.amplitudes.fIntensity;
    result.gIntensity = value.amplitudes.gIntensity;
    result.fConjugateG = value.amplitudes.fConjugateG;
    return result;
}
Amplitudes GetAmplitudes(RPRModel& model,
     double sqrtS, double cosThetaLow, double cosThetaHigh,
     const std::string* component = nullptr){
    Amplitudes result;
    vector<std::string> components;
    if(!component or *component == "total"){
        return AverageAmplitudes(model, sqrtS, cosThetaLow, cosThetaHigh);
    }
    if(*component == "uch_regge"){
        components.push_back("u_lambda");
        components.push_back("u_sigma");
        components.push_back("u_sigma1385");
    }
    else if(*component == "born_term"){
        components.push_back("s_lambda");
        components.push_back("s_sigma");
    }
    else{
        return AverageAmplitudes(
            model, sqrtS, cosThetaLow, cosThetaHigh, component);
    }

    RPRModel selectedModel = model;
    selectedModel.EnableOnly(components);
    selectedModel.EnableRescattering(false);

    return AverageAmplitudes(
        selectedModel, sqrtS, cosThetaLow, cosThetaHigh);
}
Amplitudes GetAmplitudes(RPRModel& model,
     double sqrtS, double cosThetaLow, double cosThetaHigh,
     TString component){
    std::string component_str = component.Data();
    return GetAmplitudes(model, sqrtS, cosThetaLow, cosThetaHigh, &component_str);
}
Amplitudes GetAmplitudes(RPRModel& model,
     double sqrtS, double cosTheta,
    const std::string* component = nullptr){
    Amplitudes result;
    vector<std::string> components;
    if(!component or *component == "total"){
        return PointAmplitude(model, sqrtS, cosTheta);
    }
    if(*component == "uch_regge"){
        components.push_back("u_lambda");
        components.push_back("u_sigma");
        components.push_back("u_sigma1385");
    }
    else if(*component == "born_term"){
        components.push_back("s_lambda");
        components.push_back("s_sigma");
    }
    else{
        return PointAmplitude(
            model, sqrtS, cosTheta, component);
    }
    RPRModel selectedModel = model;
    selectedModel.EnableOnly(components);
    selectedModel.EnableRescattering(false);
    return PointAmplitude(
        selectedModel, sqrtS, cosTheta);
}
Amplitudes GetAmplitudes(RPRModel& model,
     double sqrtS, double cosTheta,
    TString component){
    std::string component_str = component.Data();
    return GetAmplitudes(model, sqrtS, cosTheta, &component_str);
}

inline double NonFlipRatio(const Amplitudes& amplitudes)
{
    const double denominator = amplitudes.fIntensity + amplitudes.gIntensity;
    return denominator > 0.0
        ? std::sqrt(amplitudes.fIntensity/denominator) : 0.0;
}

inline double RelativePhaseOverPi(const Amplitudes& amplitudes)
{
    const Complex product = std::conj(amplitudes.f)*amplitudes.g;
    return std::abs(product) > 0.0
        ? std::arg(product)/Pi : 0.0;
}
