#ifndef SHKIM_RPR_COMPLEX_LIB_HH
#define SHKIM_RPR_COMPLEX_LIB_HH

// Complex 2x2/4x4 algebra used by the RPR model.
// Eigen is header-only.  The second include form works on systems where
// Eigen is installed below /usr/include/eigen3.
#if __has_include(<Eigen/Dense>)
#  include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#  include <eigen3/Eigen/Dense>
#else
#  error "Eigen/Dense was not found. Install Eigen3 or add its include path."
#endif

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace rpr {

using Complex  = std::complex<double>;
using CMatrix2 = Eigen::Matrix2cd;
using CMatrix4 = Eigen::Matrix4cd;
using CVector4 = Eigen::Vector4cd;

inline const Complex I(0.0, 1.0);

inline const CMatrix2 SigmaX =
    (CMatrix2() << 0.0, 1.0,
                   1.0, 0.0).finished();
inline const CMatrix2 SigmaY =
    (CMatrix2() << 0.0, -I,
                   I,    0.0).finished();
inline const CMatrix2 SigmaZ =
    (CMatrix2() << 1.0,  0.0,
                   0.0, -1.0).finished();
inline const std::array<CMatrix2, 3> Sigma = {SigmaX, SigmaY, SigmaZ};

inline const CMatrix4 Identity4 = CMatrix4::Identity();
inline const CMatrix4 Gamma0 =
    (CMatrix4() << 1.0, 0.0, 0.0, 0.0,
                   0.0, 1.0, 0.0, 0.0,
                   0.0, 0.0,-1.0, 0.0,
                   0.0, 0.0, 0.0,-1.0).finished();
inline const CMatrix4 Gamma1 =
    (CMatrix4() << 0.0, 0.0, 0.0, 1.0,
                   0.0, 0.0, 1.0, 0.0,
                   0.0,-1.0, 0.0, 0.0,
                  -1.0, 0.0, 0.0, 0.0).finished();
inline const CMatrix4 Gamma2 =
    (CMatrix4() << 0.0, 0.0, 0.0, -I,
                   0.0, 0.0, I,    0.0,
                   0.0, I,   0.0,  0.0,
                  -I,   0.0, 0.0,  0.0).finished();
inline const CMatrix4 Gamma3 =
    (CMatrix4() << 0.0, 0.0, 1.0, 0.0,
                   0.0, 0.0, 0.0,-1.0,
                  -1.0, 0.0, 0.0, 0.0,
                   0.0, 1.0, 0.0, 0.0).finished();
inline const std::array<CMatrix4, 4> Gamma = {Gamma0, Gamma1, Gamma2, Gamma3};
inline const CMatrix4 Gamma5 = I * Gamma0 * Gamma1 * Gamma2 * Gamma3;

// g_{mu nu} = diag(+1,-1,-1,-1).
inline const std::array<double, 4> MetricSign = {1.0, -1.0, -1.0, -1.0};
inline const std::array<CMatrix4, 4> GammaLower = {
    Gamma0, -Gamma1, -Gamma2, -Gamma3
};

inline Complex MinkowskiDotComplex(const CVector4& a, const CVector4& b)
{
    // g^nu_mu a^mu b_nu
    Complex result = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        result += MetricSign[mu] * a(mu) * b(mu);
    return result;
}

inline double Dot(const CVector4& a, const CVector4& b)
{
    return MinkowskiDotComplex(a, b).real();
}

inline CVector4 Lower(const CVector4& a)
{
    CVector4 result;
    for (int mu = 0; mu < 4; ++mu) result(mu) = MetricSign[mu] * a(mu);
    return result;
}

inline CMatrix4 Slash(const CVector4& a)
{
    const CVector4 aLower = Lower(a);
    CMatrix4 result = CMatrix4::Zero();
    for (int mu = 0; mu < 4; ++mu) result += Gamma[mu] * aLower(mu);
    return result;
}

inline CMatrix4 SigmaMuNu(int mu, int nu)
{
    if (mu < 0 || mu > 3 || nu < 0 || nu > 3)
        throw std::out_of_range("Lorentz index outside [0,3]");
    return 0.5*I*(Gamma[mu]*Gamma[nu]-Gamma[nu]*Gamma[mu]);
}

inline CMatrix4 SigmaLowerMuNu(int mu, int nu)
{
    return MetricSign[mu]*MetricSign[nu]*SigmaMuNu(mu, nu);
}

inline int LeviCivita(int a, int b, int c, int d)
{
    const int index[4] = {a,b,c,d};
    for (int i = 0; i < 4; ++i)
        for (int j = i+1; j < 4; ++j)
            if (index[i] == index[j]) return 0;
    int sign = 1;
    for (int i = 0; i < 4; ++i)
        for (int j = i+1; j < 4; ++j)
            if (index[i] > index[j]) sign = -sign;
    return sign; // epsilon^{0123}=+1
}

inline double SpatialNorm2(const CVector4& a)
{
    return std::norm(a(1))+std::norm(a(2))+std::norm(a(3));
}

// Positive-energy spinor with spin quantized along the fixed z axis.
// spin=0 means up; spin=1 means down.  The normalization is ubar*u=1.
inline CVector4 DiracSpinor(const CVector4& p, double mass, int spin)
{
    if (spin != 0 && spin != 1)
        throw std::invalid_argument("spin must be 0 (up) or 1 (down)");

    const double energy = p(0).real();
    if (mass <= 0.0 || energy + mass <= 0.0)
        throw std::invalid_argument("invalid mass or energy in DiracSpinor");

    Eigen::Vector2cd chi;
    chi << (spin == 0 ? 1.0 : 0.0), (spin == 1 ? 1.0 : 0.0);

    const CMatrix2 sigmaP = SigmaX * p(1) + SigmaY * p(2) + SigmaZ * p(3);
    const double normalization = std::sqrt((energy + mass) / (2.0 * mass));

    CVector4 result;
    result.head<2>() = normalization * chi;
    result.tail<2>() = normalization * sigmaP * chi / (energy + mass);
    return result;
}

inline Eigen::RowVector4cd DiracAdjoint(const CVector4& u)
{
    return u.adjoint() * Gamma0;
}

// Direct decomposition of the spin-3/2 projector in Eq. (9) of
// Kim et al.:
//   P_{mu nu} = g_{mu nu},
//   R_{mu nu} = gamma_mu gamma_nu
//               +(gamma_mu q_nu-q_mu gamma_nu)/M
//               +2 q_mu q_nu/M^2,
//   Delta_{mu nu} = -P_{mu nu}+R_{mu nu}/3.
//
// This grouping deliberately keeps every term of Eq. (9) explicit.  It is
// algebraically identical to the earlier transverse-P decomposition, in
// which q_mu q_nu/M^2 appeared once in both P and R and combined to 2/3.
inline double ContractP(const CVector4& a, const CVector4& b,
                        const CVector4& q, double mass)
{
    (void)q;
    (void)mass;
    Complex result = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double metric = (mu == nu) ? MetricSign[mu] : 0.0;
            result += a(mu)*metric*b(nu);
        }
    }
    return result.real();
}

inline CMatrix4 ContractR(const CVector4& a, const CVector4& b,
                          const CVector4& q, double mass)
{
    const CVector4 qLower = Lower(q);
    CMatrix4 result = CMatrix4::Zero();
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const CMatrix4 rMuNu =
                GammaLower[mu]*GammaLower[nu]
                +(GammaLower[mu]*qLower(nu)
                  -GammaLower[nu]*qLower(mu))/mass
                +2.0*Identity4*qLower(mu)*qLower(nu)/(mass*mass);
            result += a(mu) * b(nu) * rMuNu;
        }
    }
    return result;
}

inline CMatrix4 ContractDelta3Half(const CVector4& a, const CVector4& b,
                                   const CVector4& q, double mass)
{
    return -ContractP(a,b,q,mass)*Identity4
           +ContractR(a,b,q,mass)/3.0;
}

// Building blocks of Eq. (A4) in Man, Oh, and Nakayama,
// Phys. Rev. C 83, 055201 (2011):
//   theta_{mu nu} = g_{mu nu} - q_mu q_nu/M^2,
//   barGamma_mu   = gamma_mu - q_mu slash(q)/M^2.
// No complex conjugation is taken in these Lorentz contractions.
inline Complex ContractTheta(const CVector4& a, const CVector4& b,
                             const CVector4& q, double mass)
{
    const double mass2 = mass*mass;
    // g^nu_mu a^mu b_nu - p_mu p^nu a^mu b_nu / M^2
    return MinkowskiDotComplex(a,b)
           -MinkowskiDotComplex(a,q)*MinkowskiDotComplex(q,b)/mass2;
}

inline CMatrix4 ContractBarGammaPair(const CVector4& a, const CVector4& b,
                                     const CVector4& q, double mass)
{
    const double mass2 = mass*mass;
    const CMatrix4 slashQ = Slash(q);
    const CMatrix4 aBarGamma =
        Slash(a)-MinkowskiDotComplex(a,q)*slashQ/mass2;
        // g_mu a^mu - q_mu slash(q) a^mu / M^2
    const CMatrix4 barGammaB =
        Slash(b)-MinkowskiDotComplex(q,b)*slashQ/mass2;
        // g^nu b_nu - q^nu slash(q) b_nu / M^2
    return aBarGamma*barGammaB;
}

inline CMatrix4 ContractDelta5Half(const CVector4& a, const CVector4& b,
                                   const CVector4& q, double mass)
{
    // Exact contraction of the spin-5/2 projector in Eq. (A3):
    // a^{alpha1} a^{alpha2} Delta_{alpha1 alpha2}^{beta1 beta2}
    // b_{beta1} b_{beta2}.  Because both vertices contain repeated copies
    // of a and b, the explicitly symmetrized terms reduce to this form.
    const Complex aThetaB = ContractTheta(a,b,q,mass);//th_a^b
    const Complex aThetaA = ContractTheta(a,a,q,mass);//th_a^a
    const Complex bThetaB = ContractTheta(b,b,q,mass);
    const CMatrix4 aBarGammaBarGammaB =
        ContractBarGammaPair(a,b,q,mass);
        // 1/2 (th_a1^b1 th_a2^b2 + th_a1^b2 th_a2^b1) = th_a^b* th_a^b
        // - 1/5 (th_a1a2 th^b1b2) = -1/5 th_a^a* th_b^b
        // - 1/10 ( g_a1 g^b1 th_a2^b2 + g_a1 g^b2 th_a_2b^1 + g_a2 g^b1 th_a1^b2 + g_a2 g^b2 th_a1^b1 )
        // = -1/10 (4 * g_a g^b th_a^b ) = -2/5 th_a^b * g_a g^b
    return ((aThetaB*aThetaB-aThetaA*bThetaB/5.0)*Identity4
           -(2.0/5.0)*aThetaB*aBarGammaBarGammaB);
}

inline CMatrix4 ContractDelta7Half(const CVector4& a, const CVector4& b,
                                   const CVector4& q, double mass)
{
    // Exact contraction of the fully symmetrized spin-7/2 projector in
    // Eq. (A3).  The 1/36 permutation sum collapses because the three
    // momenta at each derivative vertex are identical.
    //  1/36 sigma_{a1a2a3}^{b1b2b3} -> no symmetrization needed because a1=a2=a3 and b1=b2=b3 
    // D = th_a^b*th_a^b*th_a^b - 3/7 th_a^b th_a^a th_b^b - 3/7 g_a g^b th_a^b th_a^b + 3/35 g_a g^b th_a^a th_b^b
    const Complex aThetaB = ContractTheta(a,b,q,mass);
    const Complex aThetaA = ContractTheta(a,a,q,mass);
    const Complex bThetaB = ContractTheta(b,b,q,mass);
    const CMatrix4 aBarGammaBarGammaB =
        ContractBarGammaPair(a,b,q,mass);

    const Complex scalar =
        aThetaB*aThetaB*aThetaB
        -(3.0/7.0)*aThetaA*bThetaB*aThetaB;
    const Complex gammaCoefficient =
        -(3.0/7.0)*aThetaB*aThetaB
        +(3.0/35.0)*aThetaA*bThetaB;
    return -scalar*Identity4
           -gammaCoefficient*aBarGammaBarGammaB;
}

} // namespace rpr

#endif
