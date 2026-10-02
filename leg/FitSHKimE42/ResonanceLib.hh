#ifndef SHKIM_RPR_RESONANCE_LIB_HH
#define SHKIM_RPR_RESONANCE_LIB_HH

#include "ComplexLib.hh"

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace rpr {

enum class Channel { S, U };
enum class VertexType {
    Pseudovector,       // ground-state 1/2+ Lambda/Sigma
    Pseudoscalar,       // optional spin-1/2 excited state
    HighSpinDerivative  // spin-3/2, spin-5/2, or spin-7/2 state
};

struct Resonance {
    std::string name;
    std::string family;
    Channel channel = Channel::S;
    VertexType vertex = VertexType::HighSpinDerivative;
    double spin = 0.5;
    int parity = +1;
    double mass = 0.0;
    double width = 0.0;
    double gKN = 0.0;
    double gKXi = 0.0;
    double isospin = 1.0;
    bool reggeized = false;
    double reggeIntercept = 0.0;
    double reggeSlope = 0.0;
    double reggeEta = 1.0;
    bool enabled = true;

    // A convenient fit parameter multiplying the complete complex diagram.
    // It is deliberately separate from gKN and gKXi.
    Complex amplitudeScale = Complex(1.0, 0.0);
};

using ResonanceMap = std::map<std::string, Resonance>;

enum class IntermediateMeson { Vector, Pseudoscalar };

// One M B intermediate state in Eq. (17) of Kim et al.
//
// For Vector channels, gMesonKK and gMesonKStarK generate the mixed
// K/K* exchange products retained in the paper.  For Pseudoscalar channels,
// gMesonKStarK is the K*K P coupling and both sides use K* exchange.
// isospinScale contains the charge-channel Clebsch factor.  amplitudeScale is
// intentionally separate so the complete loop can be fitted without giving
// the fitted number the interpretation of a fundamental coupling.
struct RescatteringChannel {
    std::string name;
    std::string mesonName;
    std::string baryonName;
    IntermediateMeson mesonType = IntermediateMeson::Vector;
    double mesonMass = 0.0;
    double baryonMass = 0.0;

    double gMesonKK = 0.0;
    double gMesonKStarK = 0.0;
    double gKNB = 0.0;
    double gKXiB = 0.0;
    double gKStarNB = 0.0;
    double kappaKStarNB = 0.0;
    double gKStarXiB = 0.0;
    double kappaKStarXiB = 0.0;

    double isospinScale = 1.0;
    double cutoff = 0.85; // GeV, Lambda_MB in Eq. (28)
    bool enabled = true;
    Complex amplitudeScale = Complex(1.0, 0.0);
};

using RescatteringMap = std::map<std::string, RescatteringChannel>;

inline RescatteringMap DefaultRescatteringChannels()
{
    // Couplings are from Eqs. (20)-(25).  The Lambda/Sigma pseudoscalar
    // couplings are the same SU(3) values used by the tree amplitudes.
    constexpr double gRhoPiPi = 5.94;
    constexpr double gRhoOmegaPi = 14.4; // GeV^-1
    constexpr double gKStarKPi = 6.56;
    constexpr double sqrt2 = 1.4142135623730950488;
    constexpr double sqrt3 = 1.7320508075688772935;

    RescatteringMap result;
    auto add = [&result](const RescatteringChannel& channel) {
        result[channel.name] = channel;
    };

    auto addVector = [&add](const std::string& meson, double mesonMass,
                            double gVKK, double gVKStarK,
                            const std::string& baryon, double baryonMass,
                            double gKNB, double gKXiB,
                            double gKStarNB, double kappaNB,
                            double gKStarXiB, double kappaXiB) {
        RescatteringChannel c;
        c.name = "resc_" + meson + "_" + baryon;
        c.mesonName = meson;
        c.baryonName = baryon;
        c.mesonType = IntermediateMeson::Vector;
        c.mesonMass = mesonMass;
        c.baryonMass = baryonMass;
        c.gMesonKK = gVKK;
        c.gMesonKStarK = gVKStarK;
        c.gKNB = gKNB;
        c.gKXiB = gKXiB;
        c.gKStarNB = gKStarNB;
        c.kappaKStarNB = kappaNB;
        c.gKStarXiB = gKStarXiB;
        c.kappaKStarXiB = kappaXiB;
        c.cutoff = 0.85;
        add(c);
    };

    auto addPseudoscalar = [&add](const std::string& meson, double mesonMass,
                                  double gKStarKP,
                                  const std::string& baryon,
                                  double baryonMass,
                                  double gKStarNB, double kappaNB,
                                  double gKStarXiB, double kappaXiB) {
        RescatteringChannel c;
        c.name = "resc_" + meson + "_" + baryon;
        c.mesonName = meson;
        c.baryonName = baryon;
        c.mesonType = IntermediateMeson::Pseudoscalar;
        c.mesonMass = mesonMass;
        c.baryonMass = baryonMass;
        c.gMesonKStarK = gKStarKP;
        c.gKStarNB = gKStarNB;
        c.kappaKStarNB = kappaNB;
        c.gKStarXiB = gKStarXiB;
        c.kappaKStarXiB = kappaXiB;
        c.cutoff = 0.50;
        // The paper finds PB loops nearly negligible and does not publish
        // enough charge-channel phase information to reproduce their small
        // coherent remainder uniquely.  Keep them available but opt-in.
        c.enabled = false;
        add(c);
    };

    // Vector-meson couplings: g_rhoKK=g_omegaKK=g_rhopipi/2,
    // g_phiK*K=g_rhoomega pi/sqrt(2), and
    // g_rhoK*K=g_omegaK*K=g_rhoomega pi/2.
    const double gRhoKK = gRhoPiPi/2.0;
    const double gPhiKStarK = gRhoOmegaPi/sqrt2;
    const double gRhoKStarK = gRhoOmegaPi/2.0;

    // B = Lambda.
    for (const auto& v : {
             std::array<double,3>{0.77526,  gRhoKK, gRhoKStarK},
             std::array<double,3>{0.78265,  gRhoKK, gRhoKStarK},
             std::array<double,3>{1.019461, 4.48,   gPhiKStarK}}) {
        const std::string name = (v[0] < 0.78) ? "rho" :
                                 (v[0] < 0.90) ? "omega" : "phi";
        addVector(name, v[0], v[1], v[2], "lambda", 1.115683,
                  -13.24, 3.52, 4.26, 2.66, 4.26, 1.10);
        addVector(name, v[0], v[1], v[2], "sigma", 1.192642,
                  3.58, -13.26, -2.46, -0.467, -2.46, 4.22);
    }

    addPseudoscalar("pi", 0.13804, gKStarKPi, "lambda", 1.115683,
                    4.26, 2.66, 4.26, 1.10);
    addPseudoscalar("pi", 0.13804, gKStarKPi, "sigma", 1.192642,
                    -2.46, -0.467, -2.46, 4.22);
    addPseudoscalar("eta", 0.547862, sqrt3*gKStarKPi,
                    "lambda", 1.115683, 4.26, 2.66, 4.26, 1.10);
    addPseudoscalar("eta", 0.547862, sqrt3*gKStarKPi,
                    "sigma", 1.192642, -2.46, -0.467, -2.46, 4.22);
    return result;
}

inline ResonanceMap DefaultResonances()
{
    ResonanceMap result;
    auto add = [&result](const Resonance& r) { result[r.name] = r; };
     // name, family, channel, vertextype,
     // spin, parity, mass, width, gKN, gKXi, isospin,
     // reggeized, reggeIntercept, reggeSlope, reggeEta, enabled
    add({"u_lambda", "Lambda", Channel::U, VertexType::Pseudovector,
         0.5, +1, 1.115683, 0.0, -13.24, 3.52, -1.0,
         true, -0.65, 0.94, 2.60, true, Complex(1.0,0.0)});
    add({"u_sigma", "Sigma", Channel::U, VertexType::Pseudovector,
         0.5, +1, 1.192642, 0.0, 3.58, -13.26, +1.0,
         true, -0.79, 0.87, 0.66, true, Complex(1.0,0.0)});
    add({"u_sigma1385", "Sigma", Channel::U, VertexType::HighSpinDerivative,
         1.5, +1, 1.3837, 0.0, -3.22, -3.22, +1.0,
         true, -0.27, 0.90, 0.66, true, Complex(1.0,0.0)});

    add({"s_lambda", "Lambda", Channel::S, VertexType::Pseudovector,
         0.5, +1, 1.115683, 0.0, -13.24, 3.52, -1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    add({"s_sigma", "Sigma", Channel::S, VertexType::Pseudovector,
         0.5, +1, 1.192642, 0.0, 3.58, -13.26, +1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});

    add({"Lambda1890", "Lambda", Channel::S, VertexType::HighSpinDerivative,
         1.5, +1, 1.890, 0.120, 0.84, -0.26, -1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    add({"Lambda2100", "Lambda", Channel::S, VertexType::HighSpinDerivative,
         3.5, -1, 2.100, 0.200, 2.41, 2.95, -1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    add({"Sigma2030", "Sigma", Channel::S, VertexType::HighSpinDerivative,
         3.5, +1, 2.030, 0.180, 0.82, -0.93, +1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    add({"Sigma2230", "Sigma", Channel::S, VertexType::HighSpinDerivative,
         1.5, +1, 2.230, 0.345, 0.41, 0.34, +1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    add({"Sigma2250", "Sigma", Channel::S, VertexType::HighSpinDerivative,
         3.5, -1, 2.290, 0.100, 0.59, 0.88, +1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)});
    return result;
}

// Helper for adding a hypothetical spin-1/2 s-channel resonance.  The
// nonderivative parity convention is Gamma=gamma5 for J^P=1/2+ and Gamma=1
// for J^P=1/2-.  Its amplitudeScale can be used if only a phenomenological
// strength, rather than physical couplings, is wanted.
inline Resonance MakeSpinHalfResonance(const std::string& name,
                                       const std::string& family,
                                       int parity, double mass, double width,
                                       double gKN=1.0, double gKXi=1.0,
                                       double isospin=1.0)
{
    return {name, family, Channel::S, VertexType::Pseudoscalar,
            0.5, parity, mass, width, gKN, gKXi, isospin,
            false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)};
}

// Helper for adding a derivative-coupled s-channel resonance.  The amplitude
// implementation supports every high-spin case retained in Eq. (14).
inline Resonance MakeHighSpinResonance(const std::string& name,
                                       const std::string& family,
                                       double spin, int parity,
                                       double mass, double width,
                                       double gKN=1.0, double gKXi=1.0,
                                       double isospinFactor=1.0)
{
    if (std::abs(spin-1.5) > 1.0e-12
        && std::abs(spin-2.5) > 1.0e-12
        && std::abs(spin-3.5) > 1.0e-12)
        throw std::invalid_argument(
            "high-spin resonance must have spin 3/2, 5/2, or 7/2");
    if (parity != +1 && parity != -1)
        throw std::invalid_argument("resonance parity must be +1 or -1");
    return {name, family, Channel::S, VertexType::HighSpinDerivative,
            spin, parity, mass, width, gKN, gKXi, isospinFactor,
            false, 0.0, 0.0, 1.0, true, Complex(1.0,0.0)};
}

} // namespace rpr

#endif
