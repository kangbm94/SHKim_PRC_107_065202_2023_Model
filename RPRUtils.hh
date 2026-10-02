#pragma once

#include "RPR_Model.hh"
#include "MathLib.hh"
using rpr::Complex;
using rpr::ModelResult;
using rpr::RPRModel;
using namespace rpr;
inline Amplitudes AverageAmplitudes(
    RPRModel& model, double sqrtS,
    double cosThetaLow, double cosThetaHigh,
    const std::string* component = nullptr)
{
    Amplitudes result;
    double width = cosThetaHigh - cosThetaLow;
    double dcth = 0.001;
    for (double cth_0 = cosThetaLow; cth_0 < cosThetaHigh; cth_0 += dcth) {
        double cth = cth_0 + dcth < cosThetaHigh ? cth_0 + 0.5 * dcth : 0.5 * (cth_0 + cosThetaHigh);
        double dcth_width = 2 *(cth - cth_0);
        double weight = dcth_width / width;
        const ModelResult value = component
            ? model.Evaluate(sqrtS, cth, *component)
            : model.Evaluate(sqrtS, cth);
        result.f += weight*value.amplitudes.f;
        result.g += weight*value.amplitudes.g;
        //result.fIntensity += weight*std::norm(value.amplitudes.f);
        //result.gIntensity += weight*std::norm(value.amplitudes.g);
        //result.fConjugateG += weight
        //                   * std::conj(value.amplitudes.f)
        //                   * value.amplitudes.g;
    }
    result.fIntensity = std::norm(result.f);
    result.gIntensity = std::norm(result.g);
    result.fConjugateG = std::conj(result.f) * result.g;
    return result;
}
Amplitudes GetAmplitudes(RPRModel& model,
     double sqrtS, double cosThetaLow, double cosThetaHigh,
     const std::string* component = nullptr){
    Amplitudes result;
    vector<std::string> components;
    if(!component || *component == "total"){
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
        components.push_back(*component);
    }
    for(const auto& comp : components){
        result = result + AverageAmplitudes(model, sqrtS, cosThetaLow, cosThetaHigh, &comp);
    }
    return result;
}
inline double NonFlipRatio(const Amplitudes& amplitudes)
{
    const double denominator = amplitudes.fIntensity + amplitudes.gIntensity;
    return denominator > 0.0
        ? std::sqrt(amplitudes.fIntensity/denominator) : 0.0;
}

inline double RelativePhaseOverPi(const Amplitudes& amplitudes)
{
    return std::abs(amplitudes.fConjugateG) > 0.0
        ? std::arg(amplitudes.fConjugateG)/rpr::Pi : 0.0;
}
