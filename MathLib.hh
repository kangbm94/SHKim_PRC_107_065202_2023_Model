#ifndef SHKIM_RPR_MATH_LIB_HH
#define SHKIM_RPR_MATH_LIB_HH

#include "ComplexLib.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rpr {
inline constexpr double GeV2ToMicrobarn = 389.3793721;
inline constexpr double Pi = 3.14159265358979323846;
static std::vector<std::pair<double,double>>
GaussLegendre(int n, double lower, double upper){
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

inline double KallenFunction(double x, double y, double z)
{
    return x*x + y*y + z*z - 2.0*(x*y + y*z + z*x);
}
inline double BeamMomentumToSqrtS(double pLab,
                                  double mK=0.493677,
                                  double mP=0.9382720813)
{
    const double eK = std::hypot(pLab, mK);
    return std::sqrt(mK*mK+mP*mP+2.0*mP*eK);
}
static double TwoBodyMomentum(double totalEnergy,
                                double m1, double m2)
{
    if (totalEnergy <= m1+m2) return 0.0;
    const double e2 = totalEnergy*totalEnergy;
    return std::sqrt(std::max(0.0,
        KallenFunction(e2, m1*m1, m2*m2)/(4.0*e2)));
}

struct Kinematics {
    double sqrtS = 0.0;
    double cosTheta = 0.0;
    CVector4 k1 = CVector4::Zero(); // incoming K-
    CVector4 k2 = CVector4::Zero(); // outgoing K+
    CVector4 p1 = CVector4::Zero(); // incoming proton
    CVector4 p2 = CVector4::Zero(); // outgoing Xi-
    CVector4 qS = CVector4::Zero();
    CVector4 qU = CVector4::Zero();
    double s = 0.0;
    double u = 0.0;
    double pInitial = 0.0;
    double pFinal = 0.0;
};

inline Kinematics MakeKinematics(double sqrtS, double cosTheta,
                                 double mK, double mP, double mXi)
{
    if (sqrtS < mK + mXi)
        throw std::invalid_argument("sqrt(s) is below the K+ Xi- threshold");
    if (cosTheta < -1.0 || cosTheta > 1.0)
        throw std::invalid_argument("cos(theta) must be in [-1,1]");

    Kinematics k;
    k.sqrtS = sqrtS;
    k.cosTheta = cosTheta;
    k.s = sqrtS * sqrtS;

    const double pInitial2 = KallenFunction(k.s, mK*mK, mP*mP) / (4.0*k.s);
    const double pFinal2   = KallenFunction(k.s, mK*mK, mXi*mXi) / (4.0*k.s);
    k.pInitial = std::sqrt(std::max(0.0, pInitial2));
    k.pFinal   = std::sqrt(std::max(0.0, pFinal2));

    const double sinTheta = std::sqrt(std::max(0.0, 1.0-cosTheta*cosTheta));
    const double eKInitial = std::hypot(k.pInitial, mK);
    const double ePInitial = std::hypot(k.pInitial, mP);
    const double eKFinal   = std::hypot(k.pFinal, mK);
    const double eXiFinal  = std::hypot(k.pFinal, mXi);

    // The reaction plane is x-z.  k1 x k2 points along +y.
    k.k1 << eKInitial, 0.0, 0.0, k.pInitial;
    k.p1 << ePInitial, 0.0, 0.0,-k.pInitial;
    k.k2 << eKFinal, k.pFinal*sinTheta, 0.0, k.pFinal*cosTheta;
    k.p2 << eXiFinal,-k.pFinal*sinTheta, 0.0,-k.pFinal*cosTheta;
    k.qS = k.k1 + k.p1;
    k.qU = k.p2 - k.k1;
    k.u = Dot(k.qU, k.qU);
    return k;
}
struct LabFrameKinematics {
    double beta = 0.0;
    double gamma = 1.0;
    double pFinalCM = 0.0;
    double eKFinalCM = 0.0;
};

inline CMatrix2 SpinMatrix(const CMatrix4& operator4, const Kinematics& k,
                           double mInitial, double mFinal)
{
    CMatrix2 result = CMatrix2::Zero();
    for (int sf = 0; sf < 2; ++sf) {
        const CVector4 uFinal = DiracSpinor(k.p2, mFinal, sf);
        const Eigen::RowVector4cd uFinalBar = DiracAdjoint(uFinal);
        for (int si = 0; si < 2; ++si) {
            const CVector4 uInitial = DiracSpinor(k.p1, mInitial, si);
            result(sf, si) = (uFinalBar * operator4 * uInitial)(0, 0);
        }
    }
    return result;
}

inline CMatrix2 SpinMatrixBetween(const CMatrix4& operator4,
                                  const CVector4& pInitial,
                                  double mInitial,
                                  const CVector4& pFinal,
                                  double mFinal)
{
    CMatrix2 result = CMatrix2::Zero();
    for (int sf = 0; sf < 2; ++sf) {
        const CVector4 uFinal = DiracSpinor(pFinal, mFinal, sf);
        const Eigen::RowVector4cd uFinalBar = DiracAdjoint(uFinal);
        for (int si = 0; si < 2; ++si) {
            const CVector4 uInitial = DiracSpinor(pInitial, mInitial, si);
            result(sf, si) = (uFinalBar*operator4*uInitial)(0,0);
        }
    }
    return result;
}

inline double BornFormFactor(double s, double mass,
                             double n, double cutoff)
{
    const double numerator = n * std::pow(cutoff, 4);
    const double base = numerator / (numerator + std::pow(s-mass*mass, 2));
    return std::pow(base, n);
}

enum class ResonanceForm { PaperLiteral, GaussianSquared };

inline double ResonanceFormFactor(double s, double mass, double cutoff,
                                  ResonanceForm form)
{
    const double x = (s-mass*mass)/(cutoff*cutoff);
    if (form == ResonanceForm::PaperLiteral) return std::exp(-x);
    return std::exp(-x*x);
}

inline double ReggeFactor(double s, double u, double intercept, double slope,
                          double spinOffset, double eta,
                          double s0, double cutoff)
{
    const double alpha = intercept + slope*u;
    const double cutoff2 = cutoff*cutoff;
    const double cH = std::pow(eta*cutoff2/(cutoff2-u), 2);
    return cH * slope * std::pow(s/s0, alpha-spinOffset)
           * std::tgamma(spinOffset-alpha);
}
struct Amplitudes {
    Complex f = Complex(0.0, 0.0);
    Complex g = Complex(0.0, 0.0);
    Complex fConjugateG = Complex(0.0, 0.0);
    CMatrix2 MatrixElement = CMatrix2::Zero();
    double fIntensity = 0.0;
    double gIntensity = 0.0;
};
inline Amplitudes operator +(const Amplitudes& a,
                                     const Amplitudes& b)
{
    Amplitudes result;
    result.f = a.f + b.f;
    result.g = a.g + b.g;
    result.MatrixElement = a.MatrixElement + b.MatrixElement;
    result.fIntensity = std::norm(result.f);
    result.gIntensity = std::norm(result.g);
    result.fConjugateG = std::conj(result.f) * result.g;
    return result;
}
inline Amplitudes ExtractAmplitudes(const CMatrix2& matrix)
{
    Amplitudes result;
    result.f = 0.5*(matrix(0,0)+matrix(1,1));
    result.g = 0.5*(matrix(1,0)-matrix(0,1));
    result.MatrixElement = matrix;
    result.fIntensity = std::norm(result.f);
    result.gIntensity = std::norm(result.g);
    result.fConjugateG = std::conj(result.f) * result.g;
    return result;
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
} // namespace rpr

#endif
