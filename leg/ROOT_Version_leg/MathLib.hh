//#include "../../Include/ComplexLib.hh"
#include "ComplexLib.hh"
double KallenFunction(double x, double y, double z){
    return x*x + y*y + z*z - 2*x*y - 2*y*z - 2*z*x;
}

struct Kinematics{
    double sqrt_s,cos_th;
    CVector4 k1,k2,p1,p2;
    CVector4 q_s,q_u;
    double s, u;
    double p_ini, p_fin;
};
CMatrix2 SpinMatrix(const CMatrix4& Operator, const Kinematics& kin, double m_i, double m_f){
    CMatrix2 val = CMatrix2::Zero();
    for(int sf = 0; sf < 2; ++sf){
        CVector4 uf = DiracSpinor(kin.p2, m_f, sf);
        Eigen::RowVector4cd uf_bar = uf.adjoint() * Gamma0;
        for(int si = 0; si < 2; ++si){
            CVector4 ui = DiracSpinor(kin.p1, m_i, si);
            val(sf,si) = (uf_bar * Operator * ui)(0,0);
        } 
    }
    return val;
}
double BornFormFactor(double s, double m,
    double born_n, double born_cutoff){
        double numerator = born_n * pow(born_cutoff, 4);
        double base = numerator / (numerator + pow(s - m*m, 2));
        return pow(base, born_n);
}
double ResonanceFormFactorExpo(double s, double m, double res_cutoff, int mode){
    double x = (s - m*m) / (res_cutoff * res_cutoff);
    // mode = 0: exponential, mode = 1: gaussian
    if(mode == 0) return exp(-x);
    else return exp(-x*x);
}
double ReggeFactor(double s, double u, double intercept, double slope,
     double spin_offset, double eta, double reg_s0, double reg_cutoff){
    double alpha = intercept + slope * u;
    double reg2 = reg_cutoff * reg_cutoff;
    double c_h = pow(eta * reg2/(reg2-u),2);
    return c_h * slope * pow(s/reg_s0, alpha - spin_offset) * tgamma(spin_offset - alpha);
}
