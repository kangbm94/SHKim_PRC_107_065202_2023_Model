
#include <Eigen/Dense>
#include <complex>
#include <iostream>
using Complex = std::complex<double>;
using CMatrix2 = Eigen::Matrix2cd;
using CMatrix4 = Eigen::Matrix4cd;
using CVector4 = Eigen::Vector4cd;
const Complex I(0.0, 1.0);
const CMatrix2 SigmaX =
    (CMatrix2() << 0, 1,
                   1, 0).finished();
const CMatrix2 SigmaY =
    (CMatrix2() << 0, -I,
                   I,  0).finished();
const CMatrix2 SigmaZ =
    (CMatrix2() << 1, 0,
                   0, -1).finished();
std::array<CMatrix2, 3> Sigma;
Sigma[0] = SigmaX;
Sigma[1] = SigmaY;
Sigma[2] = SigmaZ;

const CMatrix4 Identity4 = CMatrix4::Identity();

const CMatrix4 Gamma0 =
    (CMatrix4() << 1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0,-1, 0,
                   0, 0, 0,-1).finished();
const CMatrix4 Gamma1 =
    (CMatrix4() << 0, 0, 0, 1,
                   0, 0, 1, 0,
                   0,-1, 0, 0,
                  -1, 0, 0, 0).finished();
const CMatrix4 Gamma2 =
    (CMatrix4() << 0, 0, 0,-I,
                   0, 0, I, 0,
                   0, I, 0, 0,
                  -I, 0, 0, 0).finished();
const CMatrix4 Gamma3 =
    (CMatrix4() << 0, 0, 1, 0,
                   0, 0, 0,-1,
                  -1, 0, 0, 0,
                   0, 1, 0, 0).finished();  
std::array<CMatrix4, 4> Gamma;
Gamma[0] = Gamma0;
Gamma[1] = Gamma1;
Gamma[2] = Gamma2;
Gamma[3] = Gamma3;
const CMatrix4 Gamma5 = I* Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];
const CMatrix4 Metric = (CMatrix4() << 1, 0, 0, 0,
                               0,-1, 0, 0,
                               0, 0,-1, 0,
                               0, 0, 0,-1).finished();
std::array<CMatrix4, 4> GammaBar;
GammaBar[0] = Metric * Gamma[0];
GammaBar[1] = Metric * Gamma[1];
GammaBar[2] = Metric * Gamma[2];
GammaBar[3] = Metric * Gamma[3];

Complex ComplexDot(const CVector4& a, const CVector4& b){
    return a.adjoint() * Metric * b;
}
double Dot(const CVector4& a, const CVector4& b){
    return ComplexDot(a, b).real();
}
CVector4 Lower(const CVector4& a){
    /*
    a_mu = g_{mu nu} a^nu
    */
    return Metric * a;
}
CMatrix4 Slash(const CVector4& a){
    /*
    a_mu = gamma^mu_{ab} a^a
    slash(a) = gamma^mu a_mu
    = gamma^0 a^0 - gamma^1 a^1 - gamma^2 a^2 - gamma^3 a^3
    */
    CVector4 a_low = Lower(a);
    CMatrix4 ret = CMatrix4::Zero();
    for(int mu = 0; mu < 4; ++mu){
        ret +=  Gamma[mu] * a_low(mu);
    } 
    return ret;
}
CVector4 DiracSpinor(const CVector4& p, double m, int s){
    // s = 0: spin up along z, s = 1: spin down along z
    /*
    u(p,s) = sqrt((E+m)/2m) * ( chi_s, (sigma . p)/(E+m) chi_s )^T
    where chi_s is the two-component spinor for spin s 
    */
    if(s != 0 && s != 1){
        throw std::invalid_argument("spin must be 0 (up) or 1 (down)");
    }

    const double E = p(0).real();

    Eigen::Vector2cd chi;
    if(s == 0){
        chi << 1.0, 0.0;
    } else {
        chi << 0.0, 1.0;
    }

    // sigma . p = sigma_x p_x + sigma_y p_y + sigma_z p_z
    const CMatrix2 sigma_p =
        SigmaX * p(1) +
        SigmaY * p(2) +
        SigmaZ * p(3);

    const double norm = std::sqrt((E + m) / (2.0 * m));

    CVector4 u;
    u.head<2>() = norm * chi;
    u.tail<2>() = norm * (sigma_p * chi) / (E + m);

    return u;
}

double Contract_1(const CVector4& a, const CVector4& b,
     const CVector4& q, const double& mass){
    /*
        Spin-1 propagator contraction with two vectors a and b:
        The propagator is 
        P_{mu nu} = -g_{mu nu} + q_mu q_nu / m^2
        and this returns
        a^mu P_{mu nu} b^nu
    */
    CVector4 q_low = Lower(q);
    double val = 0.;
    Complex p_munu = 0.;
    for(int mu = 0; mu < 4; ++mu){
        for(int nu = 0; nu < 4; ++nu){
            p_munu = -Metric(mu,nu) + q_low(mu) * q_low(nu) / (mass*mass);
            val += (a(mu) * p_munu * b(nu)).real();
        }
    }
    return val;
}

CMatrix4 Contract_3_2(const CVector4& a, const CVector4& b,
     const CVector4& q, const double& mass){
    /*
        Spin-3/2 propagator contraction with two vectors a and b:
        r_{mu nu} = -gamma_bar_mu gamma_nu - (gamma_bar_mu q_nu - gamma_bar_nu q_mu)/m + q_mu q_nu / m^2
        and this returns
    */
    CVector4 q_low = Lower(q);
    CMatrix4 val = CMatrix4::Zero();
    CMatrix4 r_munu = CMatrix4::Zero();
    for(int mu = 0; mu < 4; ++mu){
        for(int nu = 0; nu < 4; ++nu){
            r_munu =  (
                -GammaBar[mu] * Gamma[nu] + 
                -(GammaBar[mu] * q_low(nu) - GammaBar[nu] * q_low(mu)) / mass 
                + Identity4 * q_low(mu) * q_low(nu) / (mass*mass)
            );
            val += a(mu) * r_munu * b(nu);
        }
    }
    return val;
}
CMatrix4 Contract_Delta3_2(const CVector4& a, const CVector4& b,
     const CVector4& q, const double& mass){
    
    double p = Contract_1(a, b, q, mass);
    CMatrix4 val = Contract_3_2(a, b, q, mass);
    return val - p * Identity4 / 3.;
}
CMatrix4 Contract_Delta7_2(const CVector4& a, const CVector4& b,
     const CVector4& q, const double& mass){
    double apb = Contract_1(a, b, q, mass);
    double apa = Contract_1(a, a, q, mass);
    double bpb = Contract_1(b, b, q, mass);
    CMatrix4 arb = Contract_3_2(a, b, q, mass);
    double scaler = pow(apb,3) - 3./7. * apa*bpb*apb;
    double matrix_coeff = -3./7.*pow(apb,2)  + 3./35 * apa*bpb;
    return scaler*Identity4 + matrix_coeff*arb;
    
}
    