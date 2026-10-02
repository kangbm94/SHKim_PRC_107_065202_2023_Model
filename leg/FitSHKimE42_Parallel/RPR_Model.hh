#ifndef SHKIM_RPR_MODEL_HH
#define SHKIM_RPR_MODEL_HH

#include "MathLib.hh"
#include "ResonanceLib.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace rpr {

inline constexpr double GeV2ToMicrobarn = 389.3793721;
inline constexpr double Pi = 3.14159265358979323846;

struct ModelParameters {
    // Particle masses in GeV.
    double mK = 0.493677;
    double mP = 0.9382720813;
    double mXi = 1.32171;

    // Form-factor and Regge parameters.
    double reggeS0 = 1.0;
    double reggeCutoff = 1.0;
    double bornCutoff = 0.85;
    double bornN = 2.0;
    double resonanceCutoff = 0.73;
    int formFactorPower = 1;
    ResonanceForm resonanceForm = ResonanceForm::GaussianSquared;
    //ResonanceForm resonanceForm = ResonanceForm::PaperLiteral;

    // false reproduces the validated implementation: use the pole mass in
    // the spin-3/2 and spin-7/2 projectors.
    bool useSqrtSAsProjectorMass = false;

    // Extra phenomenological normalization.  Keep this at one when fitting
    // polarization alone; it affects cross sections but not polarization.
    double crossSectionNormalization = 1.0;

    // Meson-baryon rescattering loop of Eqs. (17)-(28).  It is off by
    // default so constructing RPRModel preserves the original tree-level
    // result.  EnableRescattering() turns it on without changing any of the
    // existing amplitude/evaluation calls.
    bool includeRescattering = false;
    bool includeRescatteringPrincipalPart = true;
    bool includeRescatteringSingularPart = true;
    // Remove numerical components forbidden by parity after loop integration.
    // Set false only when inspecting the raw quadrature residual.
    bool projectRescatteringToFG = true;
    double loopMomentumMax = 3.0; // GeV
    int loopRadialPoints = 12;
    int loopCosThetaPoints = 10;
    int loopPhiPoints = 12;
    double mKStar = 0.89166;
    Complex rescatteringScale = Complex(1.0,0.0);
    // Finite i*epsilon used for the secondary K/K* propagator poles.  The
    // paper does not publish its numerical pole prescription; 0.02 GeV^2
    // gives stable low-order quadrature and can be reduced in convergence
    // studies.
    double exchangeEpsilon = 2.0e-2; // GeV^2
};

struct ModelResult {
    Kinematics kinematics;
    CMatrix2 spinMatrix = CMatrix2::Zero();
    FGAmplitudes amplitudes;

    // Reconstructed from f and g.
    double differentialCrossSection = 0.0; // microbarn/sr
    double polarization = 0.0;             // along k(K-) x k(K+)
    double polarizedCrossSection = 0.0;

    // Independent density-matrix calculation, useful as a consistency test.
    double differentialCrossSectionMatrix = 0.0;
    double polarizationMatrix = 0.0;
};

class RPRModel {
public:
    explicit RPRModel(const ModelParameters& parameters = ModelParameters())
        : parameters_(parameters), resonances_(DefaultResonances()),
          rescatterings_(DefaultRescatteringChannels()) {}

    ModelParameters& Parameters() { return parameters_; }
    const ModelParameters& Parameters() const { return parameters_; }

    ResonanceMap& Resonances() { return resonances_; }
    const ResonanceMap& Resonances() const { return resonances_; }

    RescatteringMap& Rescatterings() { return rescatterings_; }
    const RescatteringMap& Rescatterings() const { return rescatterings_; }

    bool HasRescattering(const std::string& name) const
    {
        return rescatterings_.find(name) != rescatterings_.end();
    }

    RescatteringChannel& GetRescattering(const std::string& name)
    {
        auto it = rescatterings_.find(name);
        if (it == rescatterings_.end())
            throw std::invalid_argument("unknown rescattering channel: " + name);
        return it->second;
    }

    const RescatteringChannel& GetRescattering(const std::string& name) const
    {
        auto it = rescatterings_.find(name);
        if (it == rescatterings_.end())
            throw std::invalid_argument("unknown rescattering channel: " + name);
        return it->second;
    }

    void EnableRescattering(bool enabled=true)
    {
        parameters_.includeRescattering = enabled;
    }

    void SetOverallRescatteringScale(Complex scale)
    {
        parameters_.rescatteringScale = scale;
    }

    void SetOverallRescatteringScale(double scale)
    {
        parameters_.rescatteringScale = Complex(scale,0.0);
    }

    void SetRescatteringEnabled(const std::string& name, bool enabled=true)
    {
        GetRescattering(name).enabled = enabled;
    }

    void SetRescatteringScale(const std::string& name, double scale)
    {
        GetRescattering(name).amplitudeScale = Complex(scale, 0.0);
    }

    void SetRescatteringScale(const std::string& name, Complex scale)
    {
        GetRescattering(name).amplitudeScale = scale;
    }

    void EnableOnlyRescatterings(const std::vector<std::string>& names)
    {
        const std::set<std::string> selected(names.begin(), names.end());
        for (auto& item : rescatterings_)
            item.second.enabled = selected.count(item.first) != 0;
        for (const auto& name : names)
            if (!HasRescattering(name))
                throw std::invalid_argument("unknown rescattering channel: " + name);
    }

    void EnableAllRescatterings()
    {
        for (auto& item : rescatterings_) item.second.enabled = true;
    }

    void AddRescattering(const RescatteringChannel& channel)
    {
        if (channel.name.empty())
            throw std::invalid_argument(
                "a rescattering channel must have a non-empty name");
        rescatterings_[channel.name] = channel;
    }

    bool Has(const std::string& name) const
    {
        return resonances_.find(name) != resonances_.end();
    }

    Resonance& Get(const std::string& name)
    {
        auto it = resonances_.find(name);
        if (it == resonances_.end())
            throw std::invalid_argument("unknown RPR contribution: " + name);
        return it->second;
    }

    const Resonance& Get(const std::string& name) const
    {
        auto it = resonances_.find(name);
        if (it == resonances_.end())
            throw std::invalid_argument("unknown RPR contribution: " + name);
        return it->second;
    }

    void SetEnabled(const std::string& name, bool enabled=true)
    {
        Get(name).enabled = enabled;
    }

    void EnableOnly(const std::vector<std::string>& names)
    {
        const std::set<std::string> selected(names.begin(), names.end());
        for (auto& item : resonances_)
            item.second.enabled = selected.count(item.first) != 0;
        for (const auto& name : names)
            if (!Has(name)) throw std::invalid_argument("unknown RPR contribution: " + name);
    }

    void EnableAll()
    {
        for (auto& item : resonances_) item.second.enabled = true;
    }

    void SetScale(const std::string& name, double scale)
    {
        Get(name).amplitudeScale = Complex(scale, 0.0);
    }

    void SetScale(const std::string& name, Complex scale)
    {
        Get(name).amplitudeScale = scale;
    }

    void SetCouplings(const std::string& name, double gKN, double gKXi)
    {
        Resonance& r = Get(name);
        r.gKN = gKN;
        r.gKXi = gKXi;
    }

    void AddResonance(const Resonance& resonance)
    {
        if (resonance.name.empty())
            throw std::invalid_argument("a resonance must have a non-empty name");
        resonances_[resonance.name] = resonance;
    }

    Kinematics GetKinematics(double sqrtS, double cosTheta) const
    {
        return MakeKinematics(sqrtS, cosTheta,
                              parameters_.mK, parameters_.mP, parameters_.mXi);
    }

    // Return each enabled diagram before coherent summation.
    std::map<std::string, CMatrix2>
    AmplitudeComponents(double sqrtS, double cosTheta) const
    {
        const Kinematics k = GetKinematics(sqrtS, cosTheta);
        std::map<std::string, CMatrix2> result;
        for (const auto& item : resonances_) {
            const Resonance& resonance = item.second;
            if (resonance.enabled)
                result[item.first] = ContributionAmplitude(resonance, k);
        }
        if (parameters_.includeRescattering) {
            for (const auto& item : RescatteringAmplitudeComponents(k))
                result[item.first] = item.second;
        }
        return result;
    }

    std::map<std::string, CMatrix2>
    RescatteringAmplitudeComponents(double sqrtS, double cosTheta) const
    {
        return RescatteringAmplitudeComponents(GetKinematics(sqrtS, cosTheta));
    }

    CMatrix2 RescatteringAmplitude(double sqrtS, double cosTheta) const
    {
        CMatrix2 result = CMatrix2::Zero();
        for (const auto& item :
             RescatteringAmplitudeComponents(sqrtS,cosTheta))
            result += item.second;
        return result;
    }

    CMatrix2 Amplitude(double sqrtS, double cosTheta) const
    {
        CMatrix2 result = CMatrix2::Zero();
        for (const auto& item : AmplitudeComponents(sqrtS, cosTheta))
            result += item.second;
        return result;
    }

    ModelResult Evaluate(double sqrtS, double cosTheta) const
    {
        const Kinematics k = GetKinematics(sqrtS,cosTheta);
        return MakeResult(k,Amplitude(sqrtS,cosTheta));
    }

    // Evaluate one named tree or rescattering contribution in isolation.
    // An explicit component request is independent of its enabled flag and,
    // for a rescattering channel, independent of the global loop switch.
    ModelResult Evaluate(double sqrtS, double cosTheta,
                         const std::string& componentName) const
    {
        const Kinematics k = GetKinematics(sqrtS,cosTheta);
        if (Has(componentName))
        {
            return MakeResult(k,ContributionAmplitude(Get(componentName),k));
        }
        if (HasRescattering(componentName)) {
            const CMatrix2 amplitude = parameters_.rescatteringScale
                *RescatteringContribution(GetRescattering(componentName),k);
            return MakeResult(k,amplitude);
        }
        throw std::invalid_argument(
            "unknown RPR contribution: " + componentName);
    }

    ModelResult Evaluate(double sqrtS, double cosTheta,
                         const char* componentName) const
    {
        if (componentName == nullptr)
            throw std::invalid_argument("component name must not be null");
        return Evaluate(sqrtS,cosTheta,std::string(componentName));
    }

    // ROOT TString adapter without introducing a ROOT dependency into this
    // header.  SFINAE enables it for string-like classes exposing Data().
    template<typename RootString,
             typename = decltype(std::declval<const RootString&>().Data())>
    ModelResult Evaluate(double sqrtS, double cosTheta,
                         const RootString& componentName) const
    {
        return Evaluate(sqrtS,cosTheta,componentName.Data());
    }

    // Differential cross section in the center-of-mass frame [microbarn/sr].
    double DifferentialCrossSectionCM(double sqrtS,
                                      double cosThetaCM) const
    {
        return Evaluate(sqrtS,cosThetaCM).differentialCrossSection;
    }

    double DifferentialCrossSectionCM(double sqrtS,
                                      double cosThetaCM,
                                      const std::string& componentName) const
    {
        return Evaluate(sqrtS,cosThetaCM,componentName)
            .differentialCrossSection;
    }

    double DifferentialCrossSectionCM(double sqrtS,
                                      double cosThetaCM,
                                      const char* componentName) const
    {
        if (componentName == nullptr)
            throw std::invalid_argument("component name must not be null");
        return DifferentialCrossSectionCM(
            sqrtS,cosThetaCM,std::string(componentName));
    }

    // Convert the outgoing K+ angle from the CM frame to the proton-rest
    // laboratory frame. Both arguments are cosines, not angles in degrees.
    double CMToLabCosTheta(double sqrtS, double cosThetaCM) const
    {
        ValidateCosTheta(cosThetaCM,"CM cos(theta)");
        const LabFrameKinematics lab = GetLabFrameKinematics(sqrtS);
        const double transverse = lab.pFinalCM
            *std::sqrt(std::max(0.0,1.0-cosThetaCM*cosThetaCM));
        const double longitudinal = lab.gamma
            *(lab.pFinalCM*cosThetaCM+lab.beta*lab.eKFinalCM);
        const double momentum = std::hypot(longitudinal,transverse);
        if (momentum == 0.0)
            throw std::runtime_error(
                "outgoing K+ is at rest in the laboratory frame");
        return std::clamp(longitudinal/momentum,-1.0,1.0);
    }

    // A laboratory angle can correspond to two CM angles near threshold.
    // Return every physical solution, ordered from backward to forward CM.
    std::vector<double> LabToCMCosTheta(double sqrtS,
                                       double cosThetaLab) const
    {
        ValidateCosTheta(cosThetaLab,"laboratory cos(theta)");
        const LabFrameKinematics lab = GetLabFrameKinematics(sqrtS);
        const double p = lab.pFinalCM;
        const double b = lab.beta*lab.eKFinalCM;
        const double oneMinusCos2 = std::max(
            0.0,1.0-cosThetaLab*cosThetaLab);
        const double gamma2 = lab.gamma*lab.gamma;

        const double quadraticA = p*p
            *(oneMinusCos2*gamma2+cosThetaLab*cosThetaLab);
        const double quadraticB = 2.0*oneMinusCos2*gamma2*p*b;
        const double quadraticC = oneMinusCos2*gamma2*b*b
            -cosThetaLab*cosThetaLab*p*p;
        double discriminant = quadraticB*quadraticB
            -4.0*quadraticA*quadraticC;
        const double scale = quadraticB*quadraticB
            +std::abs(4.0*quadraticA*quadraticC)+1.0;
        if (discriminant < -1.0e-13*scale) return {};
        discriminant = std::max(0.0,discriminant);

        std::vector<double> result;
        const double root = std::sqrt(discriminant);
        for (double candidate : {
                 (-quadraticB-root)/(2.0*quadraticA),
                 (-quadraticB+root)/(2.0*quadraticA)}) {
            if (candidate < -1.0-1.0e-11
                || candidate > 1.0+1.0e-11) continue;
            candidate = std::clamp(candidate,-1.0,1.0);
            // Squaring the angular equation can introduce a solution with
            // the wrong sign of the longitudinal laboratory momentum.
            if (std::abs(CMToLabCosTheta(sqrtS,candidate)-cosThetaLab)
                > 2.0e-8) continue;
            if (result.empty()
                || std::abs(candidate-result.back()) > 1.0e-10)
                result.push_back(candidate);
        }
        std::sort(result.begin(),result.end());
        return result;
    }

    // Laboratory differential cross section [microbarn/sr]. Contributions
    // from all CM branches corresponding to the requested lab angle are
    // summed with the solid-angle Jacobian.
    double DifferentialCrossSectionLab(double sqrtS,
                                       double cosThetaLab) const
    {
        return DifferentialCrossSectionLabImpl(
            sqrtS,cosThetaLab,nullptr);
    }

    double DifferentialCrossSectionLab(double sqrtS,
                                       double cosThetaLab,
                                       const std::string& componentName) const
    {
        return DifferentialCrossSectionLabImpl(
            sqrtS,cosThetaLab,&componentName);
    }

    double DifferentialCrossSectionLab(double sqrtS,
                                       double cosThetaLab,
                                       const char* componentName) const
    {
        if (componentName == nullptr)
            throw std::invalid_argument("component name must not be null");
        const std::string name(componentName);
        return DifferentialCrossSectionLabImpl(sqrtS,cosThetaLab,&name);
    }

    // Angular averages <dSigma/dOmega> over a cos(theta) interval. Bounds
    // may be supplied in either order. The azimuthal factor 2*pi cancels
    // between the integrated cross section and the interval solid angle.
    double AverageDifferentialCrossSectionCM(
        double sqrtS, double cosThetaMin, double cosThetaMax,
        int intervals=64) const
    {
        return AverageDifferentialCrossSectionCMImpl(
            sqrtS,cosThetaMin,cosThetaMax,intervals,nullptr);
    }

    double AverageDifferentialCrossSectionCM(
        double sqrtS, double cosThetaMin, double cosThetaMax,
        const std::string& componentName, int intervals=64) const
    {
        return AverageDifferentialCrossSectionCMImpl(
            sqrtS,cosThetaMin,cosThetaMax,intervals,&componentName);
    }

    double AverageDifferentialCrossSectionCM(
        double sqrtS, double cosThetaMin, double cosThetaMax,
        const char* componentName, int intervals=64) const
    {
        if (componentName == nullptr)
            throw std::invalid_argument("component name must not be null");
        const std::string name(componentName);
        return AverageDifferentialCrossSectionCMImpl(
            sqrtS,cosThetaMin,cosThetaMax,intervals,&name);
    }

    double AverageDifferentialCrossSectionLab(
        double sqrtS, double cosThetaLabMin, double cosThetaLabMax,
        int intervals=64) const
    {
        return AverageDifferentialCrossSectionLabImpl(
            sqrtS,cosThetaLabMin,cosThetaLabMax,intervals,nullptr);
    }

    double AverageDifferentialCrossSectionLab(
        double sqrtS, double cosThetaLabMin, double cosThetaLabMax,
        const std::string& componentName, int intervals=64) const
    {
        return AverageDifferentialCrossSectionLabImpl(
            sqrtS,cosThetaLabMin,cosThetaLabMax,intervals,&componentName);
    }

    double AverageDifferentialCrossSectionLab(
        double sqrtS, double cosThetaLabMin, double cosThetaLabMax,
        const char* componentName, int intervals=64) const
    {
        if (componentName == nullptr)
            throw std::invalid_argument("component name must not be null");
        const std::string name(componentName);
        return AverageDifferentialCrossSectionLabImpl(
            sqrtS,cosThetaLabMin,cosThetaLabMax,intervals,&name);
    }

    // Convenience wrappers accepting theta limits in degrees.
    double AverageDifferentialCrossSectionCMTheta(
        double sqrtS, double thetaMinDeg, double thetaMaxDeg,
        int intervals=64) const
    {
        return AverageDifferentialCrossSectionCM(
            sqrtS,CosThetaLowerBound(thetaMinDeg,thetaMaxDeg),
            CosThetaUpperBound(thetaMinDeg,thetaMaxDeg),intervals);
    }

    double AverageDifferentialCrossSectionCMTheta(
        double sqrtS, double thetaMinDeg, double thetaMaxDeg,
        const std::string& componentName, int intervals=64) const
    {
        return AverageDifferentialCrossSectionCM(
            sqrtS,CosThetaLowerBound(thetaMinDeg,thetaMaxDeg),
            CosThetaUpperBound(thetaMinDeg,thetaMaxDeg),
            componentName,intervals);
    }

    double AverageDifferentialCrossSectionLabTheta(
        double sqrtS, double thetaMinDeg, double thetaMaxDeg,
        int intervals=64) const
    {
        return AverageDifferentialCrossSectionLab(
            sqrtS,CosThetaLowerBound(thetaMinDeg,thetaMaxDeg),
            CosThetaUpperBound(thetaMinDeg,thetaMaxDeg),intervals);
    }

    double AverageDifferentialCrossSectionLabTheta(
        double sqrtS, double thetaMinDeg, double thetaMaxDeg,
        const std::string& componentName, int intervals=64) const
    {
        return AverageDifferentialCrossSectionLab(
            sqrtS,CosThetaLowerBound(thetaMinDeg,thetaMaxDeg),
            CosThetaUpperBound(thetaMinDeg,thetaMaxDeg),
            componentName,intervals);
    }

    double TotalCrossSection(double sqrtS, int intervals=400) const
    {
        // Simpson integration over cos(theta), followed by the azimuthal 2*pi.
        if (intervals < 2) intervals = 2;
        if (intervals % 2 != 0) ++intervals;
        const double h = 2.0/intervals;
        double sum = Evaluate(sqrtS, -1.0).differentialCrossSection
                     + Evaluate(sqrtS, 1.0).differentialCrossSection;
        for (int i = 1; i < intervals; ++i) {
            const double c = -1.0+i*h;
            sum += (i % 2 ? 4.0 : 2.0)
                   * Evaluate(sqrtS, c).differentialCrossSection;
        }
        return 2.0*Pi*h*sum/3.0;
    }

private:
    struct LabFrameKinematics {
        double beta = 0.0;
        double gamma = 1.0;
        double pFinalCM = 0.0;
        double eKFinalCM = 0.0;
    };

    ModelParameters parameters_;
    ResonanceMap resonances_;
    RescatteringMap rescatterings_;

    ModelResult MakeResult(const Kinematics& k,
                           const CMatrix2& spinMatrix) const
    {
        ModelResult result;
        result.kinematics = k;
        result.spinMatrix = spinMatrix;
        result.amplitudes = ExtractFG(result.spinMatrix);

        const double f2g2 = std::norm(result.amplitudes.f)
                            +std::norm(result.amplitudes.g);
        const double phaseSpace = CrossSectionFactor(result.kinematics);// 1./(64*pi^2 s)*p_F/p_i -> m1m2/(16*pi^2 s)*p_F/p_i in this work convention.
        result.differentialCrossSection = phaseSpace*f2g2;
        result.polarization = (f2g2 == 0.0)
            ? 0.0
            : 2.0*std::imag(std::conj(result.amplitudes.f)
                            *result.amplitudes.g)/f2g2;
        result.polarizedCrossSection =
            result.differentialCrossSection*result.polarization;

        const CMatrix2 density = result.spinMatrix*result.spinMatrix.adjoint();
        const double trace = density.trace().real();
        result.differentialCrossSectionMatrix = phaseSpace*0.5*trace;
        result.polarizationMatrix = (trace == 0.0)
            ? 0.0 : (density*SigmaY).trace().real()/trace;
        return result;
    }

    double CrossSectionFactor(const Kinematics& k) const
    {
        // In this code, the spinnor normalization ubar u = 1 is used, while
        // conventional normalization is ubar_{conv} u_{conv} = 2*m. 
        // The Matrix element becomes
        // M = ubar_{conv} Oper. u_{conv} = sqrt(2*m1*2*m2) ubar Oper. u. 
        // So, sigma =  1./(64pi^2 s)*p_F/p_i Int |M|^2 dOmega =  (4*m1*m2)/(64pi^2 s)*p_F/p_i Int |M|^2 dOmega
        // Hence, the factor  should be (m1*m2)/(16pi^2 s)).
        if (k.pInitial == 0.0) return 0.0;
        return parameters_.crossSectionNormalization
               * parameters_.mP * parameters_.mXi
               / (16.0*Pi*Pi*k.s)
               * k.pFinal/k.pInitial
               * GeV2ToMicrobarn;
    }

    static void ValidateCosTheta(double value, const char* description)
    {
        if (value < -1.0 || value > 1.0)
            throw std::invalid_argument(
                std::string(description)+" must be in [-1,1]");
    }

    static std::pair<double,double> OrderedCosInterval(
        double first, double second, const char* description)
    {
        ValidateCosTheta(first,description);
        ValidateCosTheta(second,description);
        if (first == second)
            throw std::invalid_argument(
                std::string(description)+" interval has zero width");
        return {std::min(first,second),std::max(first,second)};
    }

    static double CosThetaLowerBound(double thetaFirstDeg,
                                     double thetaSecondDeg)
    {
        if (thetaFirstDeg < 0.0 || thetaFirstDeg > 180.0
            || thetaSecondDeg < 0.0 || thetaSecondDeg > 180.0)
            throw std::invalid_argument(
                "theta in degrees must be in [0,180]");
        return std::min(std::cos(thetaFirstDeg*Pi/180.0),
                        std::cos(thetaSecondDeg*Pi/180.0));
    }

    static double CosThetaUpperBound(double thetaFirstDeg,
                                     double thetaSecondDeg)
    {
        if (thetaFirstDeg < 0.0 || thetaFirstDeg > 180.0
            || thetaSecondDeg < 0.0 || thetaSecondDeg > 180.0)
            throw std::invalid_argument(
                "theta in degrees must be in [0,180]");
        return std::max(std::cos(thetaFirstDeg*Pi/180.0),
                        std::cos(thetaSecondDeg*Pi/180.0));
    }

    LabFrameKinematics GetLabFrameKinematics(double sqrtS) const
    {
        if (sqrtS <= parameters_.mK+parameters_.mXi)
            throw std::invalid_argument(
                "laboratory angular conversion requires sqrt(s) above threshold");
        const double s = sqrtS*sqrtS;
        const double eKBeamLab =
            (s-parameters_.mP*parameters_.mP
             -parameters_.mK*parameters_.mK)/(2.0*parameters_.mP);
        const double pKBeamLab = std::sqrt(std::max(
            0.0,eKBeamLab*eKBeamLab-parameters_.mK*parameters_.mK));

        LabFrameKinematics result;
        result.beta = pKBeamLab/(eKBeamLab+parameters_.mP);
        result.gamma = (eKBeamLab+parameters_.mP)/sqrtS;
        result.pFinalCM = TwoBodyMomentum(
            sqrtS,parameters_.mK,parameters_.mXi);
        result.eKFinalCM = std::hypot(result.pFinalCM,parameters_.mK);
        return result;
    }

    double CMToLabCosThetaJacobian(double sqrtS,
                                   double cosThetaCM) const
    {
        const LabFrameKinematics lab = GetLabFrameKinematics(sqrtS);
        const double p = lab.pFinalCM;
        const double b = lab.beta*lab.eKFinalCM;
        const double transverse2 = p*p
            *std::max(0.0,1.0-cosThetaCM*cosThetaCM);
        const double longitudinal = lab.gamma*(p*cosThetaCM+b);
        const double momentum = std::sqrt(
            longitudinal*longitudinal+transverse2);
        if (momentum == 0.0) return 0.0;
        return lab.gamma*p*p*(p+b*cosThetaCM)
               /(momentum*momentum*momentum);
    }

    double DifferentialCrossSectionCMImpl(
        double sqrtS, double cosThetaCM,
        const std::string* componentName) const
    {
        return componentName
            ? Evaluate(sqrtS,cosThetaCM,*componentName)
                .differentialCrossSection
            : Evaluate(sqrtS,cosThetaCM).differentialCrossSection;
    }

    double DifferentialCrossSectionLabImpl(
        double sqrtS, double cosThetaLab,
        const std::string* componentName) const
    {
        double result = 0.0;
        for (double cosThetaCM :
             LabToCMCosTheta(sqrtS,cosThetaLab)) {
            const double derivative = std::abs(
                CMToLabCosThetaJacobian(sqrtS,cosThetaCM));
            if (derivative <= 1.0e-14)
                return std::numeric_limits<double>::infinity();
            result += DifferentialCrossSectionCMImpl(
                sqrtS,cosThetaCM,componentName)/derivative;
        }
        return result;
    }

    double AverageDifferentialCrossSectionCMImpl(
        double sqrtS, double first, double second, int intervals,
        const std::string* componentName) const
    {
        const auto bounds = OrderedCosInterval(
            first,second,"CM cos(theta)");
        double integral = 0.0;
        for (const auto& point :
             GaussLegendre(intervals,bounds.first,bounds.second))
            integral += point.second*DifferentialCrossSectionCMImpl(
                sqrtS,point.first,componentName);
        return integral/(bounds.second-bounds.first);
    }

    double AverageDifferentialCrossSectionLabImpl(
        double sqrtS, double first, double second, int intervals,
        const std::string* componentName) const
    {
        const auto bounds = OrderedCosInterval(
            first,second,"laboratory cos(theta)");
        const LabFrameKinematics lab = GetLabFrameKinematics(sqrtS);

        // Partition [-1,1] at every preimage of the lab interval boundaries
        // and at the possible CM-to-lab turning point. On each resulting
        // segment, membership in the requested lab interval is constant.
        std::vector<double> cuts = {-1.0,1.0};
        for (double boundary : {bounds.first,bounds.second}) {
            const auto roots = LabToCMCosTheta(sqrtS,boundary);
            cuts.insert(cuts.end(),roots.begin(),roots.end());
        }
        const double b = lab.beta*lab.eKFinalCM;
        if (b != 0.0) {
            const double turningPoint = -lab.pFinalCM/b;
            if (turningPoint > -1.0 && turningPoint < 1.0)
                cuts.push_back(turningPoint);
        }
        std::sort(cuts.begin(),cuts.end());
        cuts.erase(std::unique(cuts.begin(),cuts.end(),
            [](double x, double y) { return std::abs(x-y) < 1.0e-11; }),
            cuts.end());

        double integral = 0.0;
        for (std::size_t i = 0; i+1 < cuts.size(); ++i) {
            const double lower = cuts[i];
            const double upper = cuts[i+1];
            if (upper-lower <= 1.0e-12) continue;
            const double midpoint = 0.5*(lower+upper);
            const double cosThetaLab = CMToLabCosTheta(sqrtS,midpoint);
            if (cosThetaLab < bounds.first-1.0e-11
                || cosThetaLab > bounds.second+1.0e-11) continue;
            for (const auto& point :
                 GaussLegendre(intervals,lower,upper))
                integral += point.second*DifferentialCrossSectionCMImpl(
                    sqrtS,point.first,componentName);
        }
        return integral/(bounds.second-bounds.first);
    }

    static std::vector<std::pair<double,double>>
    GaussLegendre(int n, double lower, double upper)
    {
        if (n < 2) n = 2;
        std::vector<std::pair<double,double>> result(n);
        const int half = (n+1)/2;
        const double midpoint = 0.5*(upper+lower);
        const double halfWidth = 0.5*(upper-lower);
        for (int i = 0; i < half; ++i) {
            double z = std::cos(Pi*(i+0.75)/(n+0.5));
            double previous = 0.0;
            double derivative = 0.0;
            do {
                double p0 = 1.0;
                double p1 = z;
                for (int j = 2; j <= n; ++j) {
                    const double p = ((2.0*j-1.0)*z*p1-(j-1.0)*p0)/j;
                    p0 = p1;
                    p1 = p;
                }
                derivative = n*(z*p1-p0)/(z*z-1.0);
                previous = z;
                z = previous-p1/derivative;
            } while (std::abs(z-previous) > 2.0e-15);

            const double weight = 2.0/((1.0-z*z)*derivative*derivative);
            result[i] = {midpoint-halfWidth*z, halfWidth*weight};
            result[n-1-i] = {midpoint+halfWidth*z, halfWidth*weight};
        }
        return result;
    }

    static double TwoBodyMomentum(double totalEnergy,
                                  double m1, double m2)
    {
        if (totalEnergy <= m1+m2) return 0.0;
        const double e2 = totalEnergy*totalEnergy;
        return std::sqrt(std::max(0.0,
            KallenFunction(e2, m1*m1, m2*m2)/(4.0*e2)));
    }

    static double LoopFormFactor(const CVector4& transfer, double cutoff)
    {
        const double cutoff2 = cutoff*cutoff;
        const double ratio = cutoff2/(cutoff2+SpatialNorm2(transfer));
        return ratio*ratio;
    }

    std::array<CMatrix4,4> KStarCurrent(const CVector4& transfer,
                                       double kappa,
                                       double referenceMass) const
    {
        std::array<CMatrix4,4> current;
        for (int mu = 0; mu < 4; ++mu) {
            current[mu] = GammaLower[mu];
            for (int nu = 0; nu < 4; ++nu)
                current[mu] += I*kappa/(2.0*referenceMass)
                               *SigmaLowerMuNu(mu,nu)*transfer(nu);
        }
        return current;
    }

    static std::array<CVector4,3>
    VectorPolarizations(const CVector4& momentum, double mass)
    {
        const double qx = momentum(1).real();
        const double qy = momentum(2).real();
        const double qz = momentum(3).real();
        const double q = std::sqrt(qx*qx+qy*qy+qz*qz);
        std::array<double,3> n = {0.0, 0.0, 1.0};
        if (q > 1.0e-14) n = {qx/q, qy/q, qz/q};

        const std::array<double,3> reference =
            (std::abs(n[2]) < 0.9)
                ? std::array<double,3>{0.0,0.0,1.0}
                : std::array<double,3>{0.0,1.0,0.0};
        std::array<double,3> e1 = {
            reference[1]*n[2]-reference[2]*n[1],
            reference[2]*n[0]-reference[0]*n[2],
            reference[0]*n[1]-reference[1]*n[0]};
        const double e1norm = std::sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
        for (double& x : e1) x /= e1norm;
        const std::array<double,3> e2 = {
            n[1]*e1[2]-n[2]*e1[1],
            n[2]*e1[0]-n[0]*e1[2],
            n[0]*e1[1]-n[1]*e1[0]};

        std::array<CVector4,3> result;
        result[0] << 0.0, e1[0], e1[1], e1[2];
        result[1] << 0.0, e2[0], e2[1], e2[2];
        const double energy = momentum(0).real();
        result[2] << q/mass,
                     energy*n[0]/mass, energy*n[1]/mass, energy*n[2]/mass;
        return result;
    }

    CMatrix2 VectorProductionK(const RescatteringChannel& c,
                               const Kinematics& k,
                               const CVector4& meson,
                               const CVector4& baryon,
                               const CVector4& epsilon) const
    {
        const CVector4 transfer = k.k1-meson;
        const Complex denominator(Dot(transfer,transfer)
                                  -parameters_.mK*parameters_.mK,
                                  parameters_.exchangeEpsilon);
        const Complex scalar = I*c.gKNB*c.gMesonKK
            *Dot(k.k1+transfer,epsilon)
            /((parameters_.mP+c.baryonMass)*denominator)
            *LoopFormFactor(transfer,c.cutoff);
        const CMatrix4 op = Slash(transfer)*Gamma5;
        return scalar*SpinMatrixBetween(op,k.p1,parameters_.mP,
                                        baryon,c.baryonMass);
    }

    CMatrix2 VectorAbsorptionK(const RescatteringChannel& c,
                               const Kinematics& k,
                               const CVector4& meson,
                               const CVector4& baryon,
                               const CVector4& epsilon) const
    {
        const CVector4 transfer = meson-k.k2;
        const Complex denominator(Dot(transfer,transfer)
                                  -parameters_.mK*parameters_.mK,
                                  parameters_.exchangeEpsilon);
        const Complex scalar = I*c.gKXiB*c.gMesonKK
            *Dot(k.k2-transfer,epsilon)
            /((parameters_.mXi+c.baryonMass)*denominator)
            *LoopFormFactor(transfer,c.cutoff);
        const CMatrix4 op = Slash(transfer)*Gamma5;
        return scalar*SpinMatrixBetween(op,baryon,c.baryonMass,
                                        k.p2,parameters_.mXi);
    }

    CMatrix2 VectorProductionKStar(const RescatteringChannel& c,
                                   const Kinematics& k,
                                   const CVector4& meson,
                                   const CVector4& baryon,
                                   const CVector4& epsilon) const
    {
        const CVector4 transfer = k.k1-meson;
        const Complex denominator(Dot(transfer,transfer)
                                  -parameters_.mKStar*parameters_.mKStar,
                                  parameters_.exchangeEpsilon);
        const CVector4 transferLower = Lower(transfer);
        const CVector4 mesonLower = Lower(meson);
        const CVector4 epsilonLower = Lower(epsilon);
        const auto current = KStarCurrent(transfer,c.kappaKStarNB,
                                          parameters_.mP);
        CMatrix4 op = CMatrix4::Zero();
        for (int mu = 0; mu < 4; ++mu) {
            Complex coefficient = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                for (int alpha = 0; alpha < 4; ++alpha)
                    for (int beta = 0; beta < 4; ++beta)
                        coefficient += static_cast<double>(LeviCivita(mu,nu,alpha,beta))
                            *transferLower(nu)*mesonLower(alpha)*epsilonLower(beta);
            op += coefficient*current[mu];
        }
        const Complex scalar = 2.0*c.gMesonKStarK*c.gKStarNB/denominator
                               *LoopFormFactor(transfer,c.cutoff);
        return scalar*SpinMatrixBetween(op,k.p1,parameters_.mP,
                                        baryon,c.baryonMass);
    }

    CMatrix2 VectorAbsorptionKStar(const RescatteringChannel& c,
                                   const Kinematics& k,
                                   const CVector4& meson,
                                   const CVector4& baryon,
                                   const CVector4& epsilon) const
    {
        const CVector4 transfer = meson-k.k2;
        const Complex denominator(Dot(transfer,transfer)
                                  -parameters_.mKStar*parameters_.mKStar,
                                  parameters_.exchangeEpsilon);
        const CVector4 transferLower = Lower(transfer);
        const CVector4 mesonLower = Lower(meson);
        const CVector4 epsilonLower = Lower(epsilon);
        const auto current = KStarCurrent(transfer,c.kappaKStarXiB,
                                          parameters_.mXi);
        CMatrix4 op = CMatrix4::Zero();
        for (int mu = 0; mu < 4; ++mu) {
            Complex coefficient = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                for (int alpha = 0; alpha < 4; ++alpha)
                    for (int beta = 0; beta < 4; ++beta)
                        coefficient += static_cast<double>(LeviCivita(mu,nu,alpha,beta))
                            *transferLower(nu)*mesonLower(alpha)*epsilonLower(beta);
            op += coefficient*current[mu];
        }
        const Complex scalar = 2.0*c.gMesonKStarK*c.gKStarXiB/denominator
                               *LoopFormFactor(transfer,c.cutoff);
        return scalar*SpinMatrixBetween(op,baryon,c.baryonMass,
                                        k.p2,parameters_.mXi);
    }

    CMatrix2 PseudoscalarKStar(const RescatteringChannel& c,
                               const CVector4& kaon,
                               const CVector4& meson,
                               const CVector4& pInitial,
                               double mInitial,
                               const CVector4& pFinal,
                               double mFinal,
                               double gKStarB,
                               double kappa,
                               bool production) const
    {
        const CVector4 transfer = production ? kaon-meson : meson-kaon;
        const double mass2 = parameters_.mKStar*parameters_.mKStar;
        const Complex denominator(Dot(transfer,transfer)-mass2,
                                  parameters_.exchangeEpsilon);
        const auto current = KStarCurrent(transfer,kappa,
                                          production ? parameters_.mP
                                                     : parameters_.mXi);
        const Complex transferDotMeson = MinkowskiDotComplex(transfer,meson);
        CMatrix4 op = CMatrix4::Zero();
        for (int mu = 0; mu < 4; ++mu) {
            const Complex coefficient = -meson(mu)
                +transfer(mu)*transferDotMeson/mass2;
            op += coefficient*current[mu];
        }
        const Complex scalar = c.gMesonKStarK*gKStarB/denominator
                               *LoopFormFactor(transfer,c.cutoff);
        return scalar*SpinMatrixBetween(op,pInitial,mInitial,pFinal,mFinal);
    }

    CMatrix2 LoopNumerator(const RescatteringChannel& c,
                           const Kinematics& k,
                           const CVector4& meson,
                           const CVector4& baryon) const
    {
        CMatrix2 result = CMatrix2::Zero();
        if (c.mesonType == IntermediateMeson::Vector) {
            for (const CVector4& epsilon :
                 VectorPolarizations(meson,c.mesonMass)) {
                const CMatrix2 productionK =
                    VectorProductionK(c,k,meson,baryon,epsilon);
                const CMatrix2 productionKStar =
                    VectorProductionKStar(c,k,meson,baryon,epsilon);
                const CMatrix2 absorptionK =
                    VectorAbsorptionK(c,k,meson,baryon,epsilon);
                const CMatrix2 absorptionKStar =
                    VectorAbsorptionKStar(c,k,meson,baryon,epsilon);

                // These are the two mixed products retained in Sec. III B.
                result += absorptionKStar*productionK
                          +absorptionK*productionKStar;
            }
        }
        else {
            const CMatrix2 production = PseudoscalarKStar(
                c,k.k1,meson,k.p1,parameters_.mP,baryon,c.baryonMass,
                c.gKStarNB,c.kappaKStarNB,true);
            const CMatrix2 absorption = PseudoscalarKStar(
                c,k.k2,meson,baryon,c.baryonMass,k.p2,parameters_.mXi,
                c.gKStarXiB,c.kappaKStarXiB,false);
            result = absorption*production;
        }
        return c.amplitudeScale*c.isospinScale*result;
    }

    CMatrix2 RescatteringContribution(const RescatteringChannel& c,
                                      const Kinematics& k) const
    {
        if (!parameters_.includeRescatteringPrincipalPart
            && !parameters_.includeRescatteringSingularPart)
            return CMatrix2::Zero();
        if (parameters_.loopMomentumMax <= 0.0)
            throw std::invalid_argument("loopMomentumMax must be positive");

        const double threshold = c.mesonMass+c.baryonMass;
        const double qMax = parameters_.loopMomentumMax;
        const double energyMax = std::hypot(qMax,c.mesonMass)
                                 +std::hypot(qMax,c.baryonMass);
        const auto radial = GaussLegendre(parameters_.loopRadialPoints,
                                           threshold,energyMax);
        const auto polar = GaussLegendre(parameters_.loopCosThetaPoints,
                                         -1.0,1.0);
        int nPhi = std::max(4,parameters_.loopPhiPoints);
        // Pair every azimuth with its pi-rotated partner.  This prevents an
        // odd user-supplied grid from breaking the exact spin symmetry.
        if (nPhi % 2 != 0) ++nPhi;
        const double phiWeight = 2.0*Pi/nPhi;
        const double phaseSpaceDenominator = std::pow(2.0*Pi,3);

        auto weightedNumerator = [&](double energy, double cosQ,
                                     double phi) -> CMatrix2 {
            const double q = TwoBodyMomentum(energy,c.mesonMass,c.baryonMass);
            const double sinQ = std::sqrt(std::max(0.0,1.0-cosQ*cosQ));
            const double qx = q*sinQ*std::cos(phi);
            const double qy = q*sinQ*std::sin(phi);
            const double qz = q*cosQ;
            const double eM = std::hypot(q,c.mesonMass);
            const double eB = std::hypot(q,c.baryonMass);
            CVector4 meson;
            CVector4 baryon;
            meson << eM,qx,qy,qz;
            baryon << eB,-qx,-qy,-qz;
            const double jacobian = q*c.baryonMass/(2.0*energy)
                                    /phaseSpaceDenominator;
            // The explicit CMatrix2 return type is essential: with an auto
            // return type Eigen would return a lazy product expression that
            // refers to the temporary produced by LoopNumerator().
            return (jacobian*LoopNumerator(c,k,meson,baryon)).eval();
        };

        CMatrix2 result = CMatrix2::Zero();
        const bool poleInside = k.sqrtS > threshold+1.0e-12
                                && k.sqrtS < energyMax-1.0e-12;
        for (const auto& polarPoint : polar) {
            const double cosQ = polarPoint.first;
            for (int iPhi = 0; iPhi < nPhi; ++iPhi) {
                const double phi = (iPhi+0.5)*phiWeight;
                CMatrix2 angular = CMatrix2::Zero();

                if (poleInside) {
                    const CMatrix2 atPole = weightedNumerator(k.sqrtS,cosQ,phi);
                    if (parameters_.includeRescatteringPrincipalPart) {
                        for (const auto& radialPoint : radial) {
                            const double energy = radialPoint.first;
                            const double denominator = k.sqrtS-energy;
                            CMatrix2 quotient;
                            if (std::abs(denominator) > 1.0e-9) {
                                quotient = (weightedNumerator(energy,cosQ,phi)
                                            -atPole)/denominator;
                            }
                            else {
                                const double delta = 1.0e-6;
                                quotient = -(weightedNumerator(k.sqrtS+delta,cosQ,phi)
                                             -weightedNumerator(k.sqrtS-delta,cosQ,phi))
                                            /(2.0*delta);
                            }
                            angular += radialPoint.second*quotient;
                        }
                        angular += atPole*std::log((k.sqrtS-threshold)
                                                  /(energyMax-k.sqrtS));
                    }
                    if (parameters_.includeRescatteringSingularPart)
                        angular += Complex(0.0,-Pi)*atPole;
                }
                else if (parameters_.includeRescatteringPrincipalPart) {
                    for (const auto& radialPoint : radial) {
                        const double energy = radialPoint.first;
                        angular += radialPoint.second
                            *weightedNumerator(energy,cosQ,phi)/(k.sqrtS-energy);
                    }
                }
                result += polarPoint.second*phiWeight*angular;
            }
        }
        return parameters_.projectRescatteringToFG ? ProjectToFG(result)
                                                   : result;
    }

    std::map<std::string,CMatrix2>
    RescatteringAmplitudeComponents(const Kinematics& k) const
    {
        std::map<std::string,CMatrix2> result;
        for (const auto& item : rescatterings_) {
            if (item.second.enabled)
                result[item.first] = parameters_.rescatteringScale
                                     *RescatteringContribution(item.second,k);
        }
        return result;
    }

    CMatrix2 ContributionAmplitude(const Resonance& r,
                                   const Kinematics& k) const
    {
        const double coupling = r.isospin*r.gKN*r.gKXi;//isospin is the isospin factor I^mu_Y in Eq. (8)
        CMatrix4 operator4 = CMatrix4::Zero();
        Complex scalar = r.amplitudeScale*coupling;

        if (r.channel == Channel::U && r.vertex == VertexType::Pseudovector) {
            operator4 = Slash(k.k1)*Gamma5
                        *(Slash(k.qU)+r.mass*Identity4)
                        *Slash(k.k2)*Gamma5;
            scalar /= (parameters_.mP+r.mass)*(parameters_.mXi+r.mass);
        }
        else if (r.channel == Channel::U
                 && r.vertex == VertexType::HighSpinDerivative) {
            if (std::abs(r.spin-1.5) > 1.0e-12)
                throw std::invalid_argument("only spin-3/2 high-spin u-channel is implemented");
            const CMatrix4 delta = ContractDelta3Half(k.k1, k.k2, k.qU, r.mass);
            operator4 = (Slash(k.qU)+r.mass*Identity4)*delta
                        / std::pow(parameters_.mK, 2);
        }
        else if (r.channel == Channel::S
                 && r.vertex == VertexType::Pseudovector) {
            operator4 = Slash(k.k2)*Gamma5
                        *(Slash(k.qS)+r.mass*Identity4)
                        *Slash(k.k1)*Gamma5;
            scalar /= Complex(k.s-r.mass*r.mass, r.mass*r.width);
            scalar /= (parameters_.mP+r.mass)*(parameters_.mXi+r.mass);
            scalar *= std::pow(BornFormFactor(k.s, r.mass,
                                              parameters_.bornN,
                                              parameters_.bornCutoff),
                               parameters_.formFactorPower);
        }
        else if (r.channel == Channel::S
                 && r.vertex == VertexType::Pseudoscalar) {
            // Optional J=1/2 state.  Parity changes the vertex, not the
            // Dirac propagator. Overall i factors common to all diagrams are omitted.
            const CMatrix4 vertexGamma = (r.parity > 0) ? Gamma5 : Identity4;
            operator4 = vertexGamma*(Slash(k.qS)+r.mass*Identity4)*vertexGamma;
            scalar /= Complex(k.s-r.mass*r.mass, r.mass*r.width);
            scalar *= std::pow(ResonanceFormFactor(k.s, r.mass,
                                                   parameters_.resonanceCutoff,
                                                   parameters_.resonanceForm),
                               parameters_.formFactorPower);
        }
        else if (r.channel == Channel::S
                 && r.vertex == VertexType::HighSpinDerivative) {
            const double projectorMass = parameters_.useSqrtSAsProjectorMass
                                         ? k.sqrtS : r.mass;
            if (std::abs(r.spin-1.5) < 1.0e-12) {
                // Eq. (14): Gamma^(-P) for J=3/2.
                const CMatrix4 vertexGamma =
                    (r.parity > 0) ? Identity4 : Gamma5;
                const CMatrix4 delta =
                    ContractDelta3Half(k.k2, k.k1, k.qS, projectorMass);
                operator4 = vertexGamma*(Slash(k.qS)+r.mass*Identity4)
                            *delta*vertexGamma/std::pow(parameters_.mK, 2);
            }
            else if (std::abs(r.spin-2.5) < 1.0e-12) {
                // Eq. (14): -Gamma^(P) for J=5/2.  The parity vertex is
                // opposite to the J=3/2 and J=7/2 cases.
                const CMatrix4 vertexGamma =
                    (r.parity > 0) ? Gamma5 : Identity4;
                const CMatrix4 delta =
                    ContractDelta5Half(k.k2, k.k1, k.qS, projectorMass);
                operator4 = -vertexGamma*(Slash(k.qS)+r.mass*Identity4)
                            *delta*vertexGamma/std::pow(parameters_.mK, 4);
            }
            else if (std::abs(r.spin-3.5) < 1.0e-12) {
                // Eq. (14): Gamma^(-P) for J=7/2.
                const CMatrix4 vertexGamma =
                    (r.parity > 0) ? Identity4 : Gamma5;
                const CMatrix4 delta =
                    ContractDelta7Half(k.k2, k.k1, k.qS, projectorMass);
                operator4 = vertexGamma*(Slash(k.qS)+r.mass*Identity4)
                            *delta*vertexGamma/std::pow(parameters_.mK, 6);
            }
            else {
                throw std::invalid_argument(
                    "high-spin s-channel must be spin 3/2, 5/2, or 7/2");
            }
            scalar /= Complex(k.s-r.mass*r.mass, r.mass*r.width);
            scalar *= std::pow(ResonanceFormFactor(k.s, r.mass,
                                                   parameters_.resonanceCutoff,
                                                   parameters_.resonanceForm),
                               parameters_.formFactorPower);
        }
        else {
            throw std::invalid_argument("unsupported channel/vertex combination for " + r.name);
        }

        if (r.reggeized) {
            const double spinOffset = r.spin;
            scalar *= ReggeFactor(k.s, k.u, r.reggeIntercept, r.reggeSlope,
                                  spinOffset, r.reggeEta,
                                  parameters_.reggeS0,
                                  parameters_.reggeCutoff);
        }
        return scalar*SpinMatrix(operator4, k, parameters_.mP, parameters_.mXi);
    }
};

inline double BeamMomentumToSqrtS(double pLab,
                                  double mK=0.493677,
                                  double mP=0.9382720813)
{
    const double eK = std::hypot(pLab, mK);
    return std::sqrt(mK*mK+mP*mP+2.0*mP*eK);
}

} // namespace rpr

#endif
