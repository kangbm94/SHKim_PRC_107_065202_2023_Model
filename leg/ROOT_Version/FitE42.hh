#include "RPR_Model.hh"
#include "DiffCrossSectionReader.hh"
using namespace rpr;
//TString base_dir = "/Users/MIN/Desktop/JPARC/E42/Analysis/";
TString base_dir = "./";
TString input_dir = base_dir + "Datapoints/Polarization/";
TString E42pol_name(double cth0, double cth1){
    return input_dir + Form("E42_polarization_CM_%g_%g.tsv", cth0*100, cth1*100);
}
vector<TString> diff_files;
vector<TString> diff_names;

static vector<double> sqrtSs;
static vector<double> cth0s;
static vector<double> cth1s;
static vector<vector<double>> E42pol;// [cth][sqrtS]
static vector<vector<double>> E42polErr;
int iss = 0;
double gS2030 = 1.0;
double gS2230 = 1.0;
double gL2100 = 1.0;
double gS2250 = 1.0;
double guLambda = 1.0;
void SetModelParameters(RPRModel& model){
    model.SetScale("u_lambda", Complex(guLambda, 0.0));
    model.SetScale("Sigma2030", Complex(gS2030, 0.0));
    model.SetScale("Lambda2100", Complex(gL2100, 0.0));
    model.SetScale("Sigma2250", Complex(gS2250, 0.0));
    model.SetScale("Sigma2230", Complex(gS2230, 0.0));
    //model.SetScale("Sigma2230", Complex(0.0, 0.0));
}

static void Chi2E42Pol(int& npar, double* gin, double& f, double* par, int iflag){
    gS2030 = par[0];
    gL2100 = par[1];
    gS2250 = par[2];
    f = 0.0;
    double d_cth = 0.005;
    rpr::RPRModel model;
    SetModelParameters(model);
    //model.EnableRescattering();
    iss++;
    if(iss % 10 == 0) std::cout<<"Chi2E42Step: "<<iss<<std::endl;
    for(size_t i=0; i<cth0s.size(); i++){
        for(size_t j=0; j<sqrtSs.size(); j++){
            double sqrtS = sqrtSs[j];
            double cth0 = cth0s[i];
            double cth1 = cth1s[i];
            double pol = 0;
            double dsig = 0;
            for(double cth=cth0; cth<cth1; cth+=d_cth){
                ModelResult value = model.Evaluate(sqrtS, cth);
                pol += value.polarization*value.differentialCrossSection*d_cth;
                dsig += value.differentialCrossSection*d_cth;
            }
            pol /= dsig;
            if(dsig == 0) continue;
            const double err = E42polErr[i][j];
            f += std::pow((pol - E42pol[i][j])/err, 2);
        }
    }

}
static void Chi2E42DifPol(int& npar, double* gin, double& f, double* par, int iflag){
    gS2030 = par[0];
    gL2100 = par[1];
    gS2250 = par[2];
    //gS2230 = par[3];
    f = 0.0;
    RPRModel model;
    SetModelParameters(model);
    //model.EnableRescattering();
    iss++;
    if(iss % 10 == 0) std::cout<<"Chi2E42Step: "<<iss<<std::endl;
    for(auto d:DiffDatasets){
        double sqrtS = d.second.sqrtS;
        if(sqrtS < 1.9 or sqrtS > 2.4) continue;
        int npoints = d.second.cth.size();
        for(int ipt = 0; ipt<npoints; ipt++){
            double cth = d.second.cth[ipt];
            if(cth < 0.8) continue;
            double sig = d.second.dsig[ipt];
            double sig_err = d.second.dsigErr[ipt];
            ModelResult value = model.Evaluate(sqrtS, cth);
            double model_sig = value.differentialCrossSection;
            f += std::pow((model_sig - sig)/sig_err, 2);
        }
    }
    
    double d_cth = 0.005;
    for(size_t i=0; i<cth0s.size(); i++){
        for(size_t j=0; j<sqrtSs.size(); j++){
            double sqrtS = sqrtSs[j];
            double cth0 = cth0s[i];
            double cth1 = cth1s[i];
            double pol = 0;
            double dsig = 0;
            for(double cth=cth0; cth<cth1; cth+=d_cth){
                ModelResult value = model.Evaluate(sqrtS, cth);
                pol += value.polarization*value.differentialCrossSection*d_cth;
                dsig += value.differentialCrossSection*d_cth;
            }
            pol /= dsig;
            if(dsig == 0) continue;
            const double err = E42polErr[i][j];
            f += std::pow((pol - E42pol[i][j])/err, 2);
        }
    }

}

static void Chi2E42FwdDifPol(int& npar, double* gin, double& f, double* par, int iflag){
    gS2030 = par[0];
    gL2100 = par[1];
    gS2250 = par[2];
    gS2230 = par[3];
    f = 0.0;
    RPRModel model;
    SetModelParameters(model);
    //model.EnableRescattering();
    iss++;
    if(iss % 10 == 0) std::cout<<"Chi2E42Step: "<<iss<<std::endl;
    for(int idat = 0; idat < OldFwdData.name.size(); idat++){
        double sqrtS = OldFwdData.sqrts[idat];
        if(sqrtS < 1.9 or sqrtS > 2.4) continue;
        double sig = OldFwdData.dsig[idat];
        double sig_err = OldFwdData.dsigerr[idat];
        double model_sig = model.AverageDifferentialCrossSectionLabTheta(sqrtS, 0.0,20.0);
        f += std::pow((model_sig - sig)/sig_err, 2);
    }

    double d_cth = 0.005;
    for(size_t i=0; i<cth0s.size(); i++){
        for(size_t j=0; j<sqrtSs.size(); j++){
            double sqrtS = sqrtSs[j];
            double cth0 = cth0s[i];
            double cth1 = cth1s[i];
            double pol = 0;
            double dsig = 0;
            for(double cth=cth0; cth<cth1; cth+=d_cth){
                ModelResult value = model.Evaluate(sqrtS, cth);
                pol += value.polarization*value.differentialCrossSection*d_cth;
                dsig += value.differentialCrossSection*d_cth;
            }
            pol /= dsig;
            if(dsig == 0) continue;
            const double err = E42polErr[i][j];
            f += std::pow((pol - E42pol[i][j])/err, 2);
        }
    }

}
