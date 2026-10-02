#ifndef SHKIM_RPR_MATH_LIB_HH
#define SHKIM_RPR_MATH_LIB_HH

#include "ComplexLib.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rpr {

inline double KallenFunction(double x, double y, double z)
{
    return x*x + y*y + z*z - 2.0*(x*y + y*z + z*x);
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

struct FGAmplitudes {
    Complex f = 0.0; // spin non-flip
    Complex g = 0.0; // spin flip in M = f I - i g sigma_y
    CMatrix2 reconstructed = CMatrix2::Zero();
    double parityResidual = 0.0;
};

inline FGAmplitudes ExtractFG(const CMatrix2& amplitude)
{
    FGAmplitudes result;
    result.f = 0.5*(amplitude(0,0)+amplitude(1,1));
    result.g = 0.5*(amplitude(1,0)-amplitude(0,1));
    result.reconstructed << result.f, -result.g,
                            result.g,  result.f;
    result.parityResidual =
        (amplitude-result.reconstructed).cwiseAbs().maxCoeff();
    return result;
}

// Orthogonal projection onto the parity-allowed spin structure for
// pseudoscalar-meson + spin-1/2-baryon scattering:
// M = f I - i g sigma_y = [[f,-g],[g,f]].
inline CMatrix2 ProjectToFG(const CMatrix2& amplitude)
{
    return ExtractFG(amplitude).reconstructed;
}

} // namespace rpr

#endif
