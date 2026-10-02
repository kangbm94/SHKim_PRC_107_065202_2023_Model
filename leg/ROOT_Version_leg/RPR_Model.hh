#define GeV2_to_mub 389.379721
#include "ResonanceLib.hh"
#include "MathLib.hh"
//#include "../../Include/Dir.hh"

const double regge_s0 = 1.0;
const double regge_cutoff = 1.0;
const double born_cutoff = 0.85;
const double born_n = 2.;
const double resonance_cutoff = 0.73;

class RPRModel{
    /*
    Free proton K-p -> K+ Xi- RPR calculation without rescattering.
    */
    public:
        RPRModel(){};
        ~RPRModel(){};
    private:
        map<TString, Resonance> resonances;
        Kinematics kinematics;
    public:
        void SetKinematics(double sqrt_s, double cos_th);
        void SetResonances(map<TString, Resonance> resonances);
        CVector4 GetAmplitude();
        CMatrix2 AmplitudeComponents(TString channel, double m_i, double m_f, double scale){

            double mass, g1, g2, iso, intercept, slope, eta;
            if(channel == "u_lambda"){
                mass = ResonanceLists["u_Lambda"].mass;
                g1 = ResonanceLists["u_Lambda"].g_kn;
                g2 = ResonanceLists["u_Lambda"].g_kxi;
                iso = 1.;
                intercept = 0.5;
                slope = 0.9;
                eta = 1.;
            }
            else if(channel == "u_sigma"){
                mass = ResonanceLists["u_Sigma"].mass;
                g1 = ResonanceLists["u_Sigma"].g_kn;
                g2 = ResonanceLists["u_Sigma"].g_kxi;
                iso = 1.;
                intercept = 0.5;
                slope = 0.9;
                eta = -1.;
            }
            double prefactor = iso * g1 * g2 / (m_i + mass)*(m_f + mass);
            CMatrix4 Operator = Slash(kinematics.k1) *
                Gamma5 *(Slash(kinematics.q_u) + mass * Identity4)
                * Slash(kinematics.k2) * Gamma5;
            double reg_factor = ReggeFactor(kinematics.s, kinematics.u,
            intercept, slope, 0.5, eta, regge_s0, regge_cutoff);
            CMatrix2 amplitude = scale * prefactor * reg_factor * SpinMatrix(Operator, kinematics, m_i, m_f);
            
            return amplitude;

        }

};
void
RPRModel::SetKinematics(double sqrt_s, double cos_th){
    kinematics.sqrt_s = sqrt_s;
    kinematics.cos_th = cos_th;
    kinematics.s = sqrt_s*sqrt_s;
    kinematics.p_ini = sqrt(KallenFunction(kinematics.s, mk*mk, mp*mp))/(2.*sqrt_s);
    kinematics.p_fin = sqrt(KallenFunction(kinematics.s, mk*mk, mXi*mXi))/(2.*sqrt_s);
    double sin_th = sqrt(1. - cos_th*cos_th);
    double e_kini = hypot(kinematics.p_ini, mk);
    double e_pini = hypot(kinematics.p_ini, mp);
    double e_kfin = hypot(kinematics.p_fin, mk);
    double e_pfin = hypot(kinematics.p_fin, mXi);
    kinematics.k1 = CVector4(e_kini, 0., 0., kinematics.p_ini);
    kinematics.p1 = CVector4(e_pini, 0., 0., -kinematics.p_ini);
    kinematics.k2 = CVector4(e_kfin, kinematics.p_fin*sin_th, 0., kinematics.p_fin*cos_th);
    kinematics.p2 = CVector4(e_pfin, -kinematics.p_fin*sin_th, 0., -kinematics.p_fin*cos_th);
    kinematics.q_s = kinematics.k1 + kinematics.p1;
    kinematics.q_u = kinematics.p2 - kinematics.k1;
    kinematics.u = Dot(kinematics.q_u, kinematics.q_u);
}



