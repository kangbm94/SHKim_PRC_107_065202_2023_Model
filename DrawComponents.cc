#include "../Include/Dir.hh"
#include "FitE42.hh"
#define ModelConfig 1 // 0 -> Default, 1 -> PolaFit, 2 -> L2110
std::vector<E42Point> gE42Data;
std::vector<ForwardPoint> gForwardData;
double gChi2Pol = 0.0;
double gChi2Forward = 0.0;
double gChi2 = 0.0;
int gnpointsPol = 0;
int gnpointsForward = 0;
int gnpoints = 0;
bool HelicityBasis = 0;
double crosssection_min = 10;
figdir = "figs_norm/";
constexpr int kForwardFitIntervals = 16;
std::vector<E42Point> LoadE42Data(const std::string& filename){
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open " + filename);
    std::vector<E42Point> result;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        E42Point point;
        double left = 0.0, right = 0.0;
        if (!(row >> point.sqrtS >> point.cosTheta >> point.polarization
                  >> point.error >> left >> right)) continue;
        point.cosThetaLow = point.cosTheta - left;
        point.cosThetaHigh = point.cosTheta + right;
        result.push_back(point);
    }
    if (result.empty()) throw std::runtime_error("No E42 points in " + filename);
    return result;
}

std::vector<ForwardPoint> LoadForwardData(const std::string& filename){
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open " + filename);
    std::vector<ForwardPoint> result;
    std::string line;
    std::getline(input, line);
    constexpr double mK = 0.493677;
    constexpr double mP = 0.9382720813;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = SplitTabs(line);
        if (fields.size() < 4) continue;
        ForwardPoint point;
        point.experiment = fields[0];
        point.pLab = std::stod(fields[1]);
        point.crossSection = std::stod(fields[2]);
        point.error = std::stod(fields[3]);
        const double eK = std::sqrt(point.pLab*point.pLab + mK*mK);
        point.sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
        result.push_back(point);
    }
    return result;
}
constexpr std::array<std::array<double,2>,4> kAngularBins = {{
    {{0.85, 0.90}}, {{0.90, 0.97}}, {{0.97, 1.00}}, {{0.85, 1.00}}
}};
void ConfigureModel(RPRModel& model){
    cout<<"Configuring model..."<<endl;
    model.SetEnabled("Sigma2230", 0);
    auto& sigma2030 = model.Get("Sigma2030");
    sigma2030.amplitudeScale = Complex(0.55, 0.0);
    sigma2030.mass = 2.040;
    sigma2030.width = 0.200;
    auto& lambda2100 = model.Get("Lambda2100");
    lambda2100.mass = 2.090;
    lambda2100.width = 0.150;
}
void ConfigureModelLambda2110(RPRModel& model){
    model.SetEnabled("Sigma2230", 0);
    cout<<"Configuring L2110 model..."<<endl;
    model.SetEnabled("Sigma2230", 0);
    #if 1
    auto& sigma2030 = model.Get("Sigma2030");
//    sigma2030.amplitudeScale = Complex(0.641964, 0.0);
    sigma2030.mass = 2.040;
    sigma2030.width = 0.200;
    auto& lambda2100 = model.Get("Lambda2100");
//    lambda2100.amplitudeScale = Complex(0.90001, 0.0);
    lambda2100.mass = 2.090;
    lambda2100.width = 0.150;
    #endif
    Resonance lambda2110 = {"Lambda2110", "Lambda", Channel::S, VertexType::HighSpinDerivative,
         2.5, 1, 2.110, 0.200, 1. , 1. , -1.0,
         false, 0.0, 0.0, 1.0, true, Complex(1,0.0)};
    lambda2110.mass = 2.050;
    lambda2110.width = 0.200;
    //lambda2110.amplitudeScale = Complex(-2.38519, 0.0);
    lambda2110.amplitudeScale = Complex(2.38519, 0.0);
    model.AddResonance(lambda2110);
}
string GetKey(const std::string& component,int flip_nonflip, int real_imag,   double cosThetaLow, double cosThetaHigh){
    string key = "";
    if(flip_nonflip == 0) key += "f_";
    else key += "g_";
    key += component;
    if(real_imag == 0) key += "_real_";
    else key += "_imag_";
    key += to_string(cosThetaLow) + "_" + to_string(cosThetaHigh);
    return key;
}
string GetPairedKey(const std::string& component1, const std::string& component2, int config,   double cosThetaLow, double cosThetaHigh){
    string key = "";
    if(config == 0) key += "fg_";
    else if(config == 1) key += "intensity_";
    else if(config == 2) key += "Polarization_";
    else if(config == 3) key += "nonflip2_";
    else if(config == 4) key += "flip2_";
    key += component1 + "_" + component2 + "_";
    key += to_string(cosThetaLow) + "_" + to_string(cosThetaHigh);
    return key;
}
int ColorForBin(double cosTheta){
    if (abs(cosTheta - (1+0.85)/2) < 1.0e-6) return kBlack;
    if (0.85 < cosTheta and cosTheta < 0.90 ) return kBlue + 1;
    if (0.90 < cosTheta and cosTheta < 0.97) return kRed + 1;
    if (0.97 < cosTheta and cosTheta < 1.00) return kGreen + 2;
    else return kGreen + 2;
}
struct AmplitudeComponentStyle {
    std::string name;
    std::string label;
    int color = kBlack;
    bool fullSum = false;
};
struct CrossSectionCurve {
    AmplitudeComponentStyle style;
    TGraph* graph = nullptr;
};
void StyleCrossSectionCurve(CrossSectionCurve& curve)
{
    curve.graph->SetLineColor(curve.style.color);
    curve.graph->SetLineWidth(curve.style.fullSum ? 3 : 2);
    curve.graph->SetLineStyle(curve.style.fullSum ? 1 : 2);
}

struct AmplitudeGraphs {
    AmplitudeComponentStyle style;
    TGraph* fReal = nullptr;
    TGraph* fImaginary = nullptr;
    TGraph* gReal = nullptr;
    TGraph* gImaginary = nullptr;
};
std::vector<AmplitudeComponentStyle>
AmplitudeComponents()
{
    std::vector<AmplitudeComponentStyle> result = {
        {"", "full coherent sum", kBlack, true},
        {"Sigma2250", "#Sigma(2250)", kGreen + 2, false},
        {"Sigma2030", "#Sigma(2030)", kRed, false},
        {"Lambda2100", "#Lambda(2100)", kBlue, false},
    };
    result.push_back({"Lambda1890", "#Lambda(1890)", kOrange + 7, false});
    result.push_back({"u_lambda", "Reggeized u-channel #Lambda",
                      kCyan + 2, false});
    return result;
}
double AveragePolarization(RPRModel& model, double sqrtS,
                           double cosThetaLow, double cosThetaHigh)
{
    const double middle = 0.5*(cosThetaLow + cosThetaHigh);
    const double halfWidth = 0.5*(cosThetaHigh - cosThetaLow);
    double numerator = 0.0;
    double denominator = 0.0;
    double dcth = 0.001;
    for (double cth_0 = cosThetaLow; cth_0 + dcth <= cosThetaHigh; cth_0 += dcth) {
        double cth = cth_0 + 0.5*dcth;
        const ModelResult value = model.Evaluate(sqrtS, cth);
        numerator += value.polarization * value.differentialCrossSection*dcth;
        denominator += value.differentialCrossSection*dcth;
    }
    if (!(denominator > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return numerator/denominator;
}
void DrawWithData(RPRModel& model,
                    const std::vector<ForwardPoint>& forwardData){
    auto* canvasPol = new TCanvas("canvasPol", "canvasPol", 2000, 1000);
    auto* frameP = gPad->DrawFrame(1.9, -1.05, 2.4, 1.05);
    frameP->SetTitle(";#sqrt{s} [GeV];P_{#Xi}");
    auto* legendP = new TLegend(0.25, 0.65, 0.65, 0.89);
    legendP->SetBorderSize(0);
    legendP->SetFillStyle(0);
    std::vector<TGraphErrors*> polarizationData;
    std::vector<TGraph*> polarizationCurves;
    gnpointsPol = 0;
    gnpointsForward = 0;
    gnpoints = 0;
    gChi2 = 0.0;
    gChi2Pol = 0.0;
    gChi2Forward = 0.0;
    for (const auto& bin : kAngularBins) {
        const double centre = 0.5*(bin[0] + bin[1]);
        const int color = ColorForBin(centre);
        auto* data = new TGraphErrors();
        auto* curve = new TGraph();
        for (double sqrtS = 1.9; sqrtS <= 2.35; sqrtS += 0.005) {
            const double p = AveragePolarization(model, sqrtS,
                                                  bin[0], bin[1]);
            curve->SetPoint(curve->GetN(), sqrtS, p);
        }
        for (const auto& point : gE42Data) {
            if (std::abs(point.cosTheta-centre) > 1.0e-6) continue;
            const int n = data->GetN();
            gnpointsPol++;
            data->SetPoint(n, point.sqrtS, point.polarization);
            data->SetPointError(n, 0.0, point.error);
            gChi2Pol += pow((point.polarization - curve->Eval(point.sqrtS))/point.error,2);
        }
        data->SetMarkerStyle(20);
        data->SetMarkerSize(3.0);
        data->SetMarkerColor(color);
        data->SetLineColor(color);
        curve->SetLineColor(color);
        curve->SetLineWidth(2);
        curve->Draw("L SAME");
        data->Draw("P SAME");
        const std::string entry = Form("%.2f < cos#theta < %.2f",
                                       bin[0], bin[1]);
        legendP->AddEntry(data, entry.c_str(), "p");
        polarizationData.push_back(data);
        polarizationCurves.push_back(curve);
    }
    legendP->Draw();
    canvasPol->SaveAs(figdir + "Polarization_wData.pdf");
    TCanvas* canvasCrossSection = new TCanvas("canvasCrossSection", "canvasCrossSection", 2000, 1000);
    canvasCrossSection->cd();
    TString forwardTitle = "p_{K^{-}}^{lab} [GeV/c];d#sigma/d#Omega_{lab} [#mub/sr]";
    double forwardMaximum = 0.0;
    auto* dataX = new TGraphErrors();
    std::vector<CrossSectionCurve> forwardCurves;
    for (const auto& style : AmplitudeComponents()) {
        CrossSectionCurve curve;
        curve.style = style;
        curve.graph = new TGraph();
        forwardCurves.push_back(curve);
    }
    constexpr double mK = 0.493677;
    constexpr double mP = 0.9382720813;
    ofstream ForwardCrossSection(figdir + "forward_curve.dat");
    ForwardCrossSection<<"#sqrtS [GeV]  dSigma/dOmega [microbarn/sr]"<<endl;
    for (double pLab = 1.15; pLab <= 2.72; pLab += 0.015) {
        const double eK = std::sqrt(pLab*pLab + mK*mK);
        const double sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
        int ic = 0;
        for (auto& curve : forwardCurves) {
            const double prediction = curve.style.fullSum
                ? model.AverageDifferentialCrossSectionLabTheta(
                    sqrtS, 0.0, 20.0, kForwardFitIntervals)
                : model.AverageDifferentialCrossSectionLabTheta(
                    sqrtS, 0.0, 20.0, curve.style.name,
                    kForwardFitIntervals);
            //curve.graph->SetPoint(curve.graph->GetN(), pLab, prediction);
            curve.graph->SetPoint(curve.graph->GetN(), sqrtS, prediction);
            forwardMaximum = std::max(forwardMaximum, prediction);
            if(ic == 0) ForwardCrossSection<<sqrtS<<"  "<<prediction<<endl;
            ic++;
        }
    }
    for (const auto& point : forwardData) {
        const int n = dataX->GetN();
        //dataX->SetPoint(n, point.pLab, point.crossSection);
        dataX->SetPoint(n, point.sqrtS, point.crossSection);
        bool exclude = 0;
        if(point.crossSection < crosssection_min and point.sqrtS < 2.0) exclude = 1;
        if(exclude) continue;
        gnpointsForward++;
        dataX->SetPointError(n, 0.0, point.error);
        forwardMaximum = std::max(
            forwardMaximum, point.crossSection+point.error);
        gChi2Forward += pow((point.crossSection - forwardCurves[0].graph->Eval(point.sqrtS))/point.error,2);
    }
    dataX->SetMarkerStyle(20);
    dataX->SetMarkerSize(3.0);
    auto* frameX = gPad->DrawFrame(
        //1.15, 0.0, 2.72, std::max(1.0, 1.15*forwardMaximum));
        //1., 0.0, 3.0, std::max(1.0, 1.15*forwardMaximum));
        1.8, 0.0, 2.6, std::max(1.0, 1.15*forwardMaximum));
    frameX->SetTitle(forwardTitle);
    for (auto& curve : forwardCurves) {
        StyleCrossSectionCurve(curve);
        curve.graph->Draw("L SAME");
    }
    dataX->Draw("P SAME");
    auto* legendX = new TLegend(0.55, 0.48, 0.8, 0.89);
    legendX->SetBorderSize(0);
    legendX->SetFillStyle(0);
    legendX->AddEntry(dataX, "forward data", "p");
    legendX->AddEntry(forwardCurves.front().graph,
             "full model (not in #chi^{2})", "l");
    for (std::size_t i = 1; i < forwardCurves.size(); ++i)
        legendX->AddEntry(forwardCurves[i].graph,
                          forwardCurves[i].style.label.c_str(), "l");
    legendX->Draw();
    gnpoints = gnpointsPol + gnpointsForward;
    gChi2 = gChi2Pol + gChi2Forward;
    canvasCrossSection->SaveAs(figdir + "CrossSection_wData.pdf");

}
void DrawWithTheory(RPRModel& model){
    vector<TString> components = 
    {
        "s_lambda", "s_sigma", "u_lambda", "u_sigma", "u_sigma1385",
        "Lambda2100", "Sigma2030", "Sigma2250", "Lambda1890", "total"
    };
    vector<double> angle_bin = {0.97, 1.00};
    vector<TString> NF_F = {"nonflip", "flip"};
    vector<TString> Re_Im = {"real", "imag"};
    map<TString, TGraph*> Graphs_Local;
    auto key = [](TString component, TString f_nf, TString r_i){
        TString val = f_nf + "_" + component + "_" + r_i;
        return val;
    };
    double scale_factor = 1.;
    double s_min = 1.83, s_max = 2.61, ds = 0.01;
    for(const auto& component : components){
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                TString graph_key = key(component, flip_nonflip, real_imag);
                Graphs_Local[graph_key] = new TGraph();
                Graphs_Local[graph_key] -> SetLineWidth(3);
                int color = flip_nonflip == "nonflip" ? kBlue : kRed;
                int style = real_imag == "real" ? kSolid : kDashed;
                Graphs_Local[graph_key]->SetLineColor(color);
                Graphs_Local[graph_key]->SetLineStyle(style);
                for(double sqrtS = s_min; sqrtS <= s_max; sqrtS += ds){
                    double cth_low = angle_bin[0], cth_high = angle_bin[1];
                    auto amplitude = GetAmplitudes(model, sqrtS, cth_low, cth_high, component);
                    double value = 0.0;
                    if(flip_nonflip == "nonflip" and real_imag == "real") value = amplitude.f.real();
                    else if(flip_nonflip == "nonflip" and real_imag == "imag") value = amplitude.f.imag();
                    else if(flip_nonflip == "flip" and real_imag == "real") value = amplitude.g.real();
                    else if(flip_nonflip == "flip" and real_imag == "imag") value = amplitude.g.imag();
                    Graphs_Local[graph_key]->SetPoint(Graphs_Local[graph_key]->GetN(), sqrtS, value*scale_factor);
                }
            }
        }
    }

    TString dat_dir = "KN-KXi_data/";
    vector<TString> components_SHKim = 
    {
        "schL", "schS", "uchL", "uchS", "uchSs",
        "L2100", "S2030", "S2250", "L1890", "full"
    };
    map<TString, TString> components_map = {
        {"schL", "s_lambda"}, {"schS", "s_sigma"}, {"uchL", "u_lambda"}, {"uchS", "u_sigma"},
        {"uchSs", "u_sigma1385"}, {"L2100", "Lambda2100"}, {"S2030", "Sigma2030"},
        {"S2250", "Sigma2250"}, {"L1890", "Lambda1890"}, {"full", "total"}
    };
    auto key_shkim = [] (TString f_nf, TString r_i, TString component){
        TString val = f_nf + "_" + r_i + "/" + "PLab-Py_"+component + ".dat";
        return val;
    };
    map<TString, TGraph*> Graphs_SHKim;
    for(const auto& flip_nonflip : NF_F){
        for(const auto& real_imag : Re_Im){
            for(const auto& component : components_SHKim){
                TString graph_key = key(component, flip_nonflip, real_imag);
                Graphs_SHKim[graph_key] = new TGraph();
                Graphs_SHKim[graph_key]->SetLineWidth(2);
                int color = flip_nonflip == "nonflip" ? kBlue + 3 : kRed + 3;
                int style = real_imag == "real" ? kSolid : kDashed;
                Graphs_SHKim[graph_key]->SetLineColor(color);
                Graphs_SHKim[graph_key]->SetLineStyle(style);
                if(component == "L1890") continue;
                TString file_name = key_shkim(flip_nonflip, real_imag, component);
                ifstream input(dat_dir + file_name);
                if(!input) throw runtime_error("Cannot open "+dat_dir+file_name);
                string line;
                while(getline(input, line)){
                    if(line.empty() or line[0] == '#') continue;
                    istringstream row(line);
                    double sqrtS, value;
                    if(!(row >> sqrtS >> value)) continue;
                    Graphs_SHKim[graph_key]->SetPoint(Graphs_SHKim[graph_key]->GetN(), sqrtS, value);
                }
            }
            TString graph_total_key = key("full", flip_nonflip, real_imag);
            TString graph_L1890_key = key("L1890", flip_nonflip, real_imag);
            for(int ip = 0; ip < Graphs_SHKim[graph_total_key]->GetN(); ++ip){
                double sqrtS, value;
                Graphs_SHKim[graph_total_key]->GetPoint(ip, sqrtS, value);
                for(const auto& component : components_SHKim){
                    if(component == "full") continue;
                    if(component == "L1890") continue;
                    TString graph_key = key(component, flip_nonflip, real_imag);
                    double value_component = Graphs_SHKim[graph_key]->Eval(sqrtS);
                    value -= value_component;
                }
                Graphs_SHKim[graph_L1890_key]->SetPoint(ip, sqrtS, value);
            }
        }
    }
    gSystem->mkdir(figdir + "Components_Comparison/", true);
    gSystem->mkdir(figdir + "Components_Comparison/DatFiles/", true);
    for(auto comp_map: components_map){
        TString component_SHKim = comp_map.first;
        TString component_local = comp_map.second;
        TString ct = Form("Amplitude_%s", component_local.Data());
        TCanvas* canvas = new TCanvas(ct,ct, 2000, 1000);
        TLegend* legend = new TLegend(0.55, 0.75, 0.85, 0.89);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        double y_min = 99990.0, y_max = -99990.0, y_range = 0.0;
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                if(real_imag == "real"){
                    legend->AddEntry(Graphs_Local[key(component_local, flip_nonflip, real_imag)],
                     Form("KBM_%s",flip_nonflip.Data()), "l");
                    legend->AddEntry(Graphs_SHKim[key(component_SHKim, flip_nonflip, real_imag)],
                        Form("SHKim_%s",flip_nonflip.Data()), "l");
                }
                TString graph_key_local = key(component_local, flip_nonflip, real_imag);
                TString graph_key_SHKim = key(component_SHKim, flip_nonflip, real_imag);
                TGraph* graph_local = Graphs_Local[graph_key_local];
                TGraph* graph_SHKim = Graphs_SHKim[graph_key_SHKim];
                for(int ip = 0; ip < graph_local->GetN(); ++ip){
                    double sqrtS, value;
                    graph_local->GetPoint(ip, sqrtS, value);
                    y_min = std::min(y_min, value);
                    y_max = std::max(y_max, value);
                }
                for(int ip = 0; ip < graph_SHKim->GetN(); ++ip){
                    double sqrtS, value;
                    graph_SHKim->GetPoint(ip, sqrtS, value);
                    y_min = std::min(y_min, value);
                    y_max = std::max(y_max, value);
                }
                y_range = y_max - y_min;
            }
        }
        gPad->DrawFrame(1.83, y_min - 0.1*y_range, 2.61, y_max + 0.1*y_range);
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                TString graph_key_local = key(component_local, flip_nonflip, real_imag);
                TString graph_key_SHKim = key(component_SHKim, flip_nonflip, real_imag);
                TGraph* graph_local = Graphs_Local[graph_key_local];
                TGraph* graph_SHKim = Graphs_SHKim[graph_key_SHKim];
                graph_local->Draw("L SAME");
                graph_SHKim->Draw("L SAME");
            }
        }
        legend->Draw();
        canvas->SaveAs(figdir + "Components_Comparison/" + Form("/Amplitude_%s.pdf", component_local.Data()));
        if(component_local == "s_lambda" or component_local == "s_sigma"){
            TCanvas* canvas_local = new TCanvas(ct+"_local",ct+"_local", 2000, 1000);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            double y_min = 9999, y_max = -9999, y_range = 0.0;
            for(const auto& flip_nonflip : NF_F){
                for(const auto& real_imag : Re_Im){
                    TString graph_key_local = key(component_local, flip_nonflip, real_imag);
                    TGraph* graph_local = Graphs_Local[graph_key_local];
                    for(int ip = 0; ip < graph_local->GetN(); ++ip){
                        double sqrtS, value;
                        graph_local->GetPoint(ip, sqrtS, value);
                        y_min = std::min(y_min, value);
                        y_max = std::max(y_max, value);   
                    }
                }
            }
            gPad->DrawFrame(1.83, y_min - 0.1*y_range, 2.61, y_max + 0.1*y_range);
            for(const auto& flip_nonflip : NF_F){
                for(const auto& real_imag : Re_Im){
                    TString graph_key_local = key(component_local, flip_nonflip, real_imag);
                    TGraph* graph_local = Graphs_Local[graph_key_local];
                    graph_local->Draw("L SAME");  
                }
            }
            canvas_local->SaveAs(figdir + "Components_Comparison/" + Form("/Amplitude_%s_local.pdf", component_local.Data()));
        }
        ofstream outfile_local(figdir + "Components_Comparison/DatFiles/KBM_" + component_local + ".dat");
        outfile_local<<"# sqrtS [GeV]  f_real  f_imag  g_real  g_imag"<<endl;
        ofstream outfile_shkim(figdir + "Components_Comparison/DatFiles/SHKim_" + component_local + ".dat");
        outfile_shkim<<"# sqrtS [GeV]  f_real  f_imag  g_real  g_imag"<<endl;
        for(double sqrtS = s_min; sqrtS <= s_max; sqrtS += ds){
            double f_real_local = Graphs_Local[key(component_local, "nonflip", "real")]->Eval(sqrtS);
            double f_imag_local = Graphs_Local[key(component_local, "nonflip", "imag")]->Eval(sqrtS);
            double g_real_local = Graphs_Local[key(component_local, "flip", "real")]->Eval(sqrtS);
            double g_imag_local = Graphs_Local[key(component_local, "flip", "imag")]->Eval(sqrtS);
            outfile_local<<sqrtS<<" "<<f_real_local<<" "<<f_imag_local<<" "<<g_real_local<<" "<<g_imag_local<<endl;
            double f_real_shkim = Graphs_SHKim[key(component_SHKim, "nonflip", "real")]->Eval(sqrtS);
            double f_imag_shkim = Graphs_SHKim[key(component_SHKim, "nonflip", "imag")]->Eval(sqrtS);
            double g_real_shkim = Graphs_SHKim[key(component_SHKim, "flip", "real")]->Eval(sqrtS);
            double g_imag_shkim = Graphs_SHKim[key(component_SHKim, "flip", "imag")]->Eval(sqrtS);
            outfile_shkim<<sqrtS<<" "<<f_real_shkim<<" "<<f_imag_shkim<<" "<<g_real_shkim<<" "<<g_imag_shkim<<endl;
        }
    }
}
void DrawAngularComponents(RPRModel& model, double sqrtS){
    vector<TString> components = 
    {
        "s_lambda", "s_sigma", "u_lambda", "u_sigma", "u_sigma1385",
        "Lambda2100", "Sigma2030", "Sigma2250", "Lambda1890"
    };
    vector<TString> NF_F = {"nonflip", "flip"};
    vector<TString> Re_Im = {"real", "imag"};
    map<TString, TGraph*> Graphs_Local;
    auto key = [](TString component, TString f_nf, TString r_i){
        TString val = f_nf + "_" + component + "_" + r_i;
        return val;
    };
    double cth_min = -1.0, cth_max = 1.0, dcth = 0.01;
    for(const auto& component : components){
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                TString graph_key = key(component, flip_nonflip, real_imag);
                Graphs_Local[graph_key] = new TGraph();
                Graphs_Local[graph_key] -> SetLineWidth(3);
                int color = flip_nonflip == "nonflip" ? kBlue : kRed;
                int style = real_imag == "real" ? kSolid : kDashed;
                Graphs_Local[graph_key]->SetLineColor(color);
                Graphs_Local[graph_key]->SetLineStyle(style);
                for(double cth = cth_min; cth <= cth_max; cth += dcth){
                    auto amplitude = GetAmplitudes(model, sqrtS, cth,component);
                    double value = 0.0;
                    if(flip_nonflip == "nonflip" and real_imag == "real") value = amplitude.f.real();
                    else if(flip_nonflip == "nonflip" and real_imag == "imag") value = amplitude.f.imag();
                    else if(flip_nonflip == "flip" and real_imag == "real") value = amplitude.g.real();
                    else if(flip_nonflip == "flip" and real_imag == "imag") value = amplitude.g.imag();
                    Graphs_Local[graph_key]->SetPoint(Graphs_Local[graph_key]->GetN(), cth, value);
                }
            }
        }
    }
    gSystem->mkdir(figdir + "Angular/", true);
    gSystem->mkdir(figdir + "Angular/DatFiles/", true);
    for(const auto& component : components){
        TString ct = Form("Angular_%s", component.Data());
        TCanvas* canvas = new TCanvas(ct,ct, 2000, 1000);
        TLegend* legend = new TLegend(0.55, 0.75, 0.85, 0.89);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        double y_min = 99990.0, y_max = -99990.0, y_range = 0.0;
        ofstream outfile_angular(figdir + "Angular/DatFiles/Angular_" + component + ".dat");
        outfile_angular<<"# cosTheta  f_real  f_imag  g_real  g_imag"<<endl;
        for(double cth = cth_min; cth <= cth_max; cth += dcth){
            double f_real = Graphs_Local[key(component, "nonflip", "real")]->Eval(cth);
            double f_imag = Graphs_Local[key(component, "nonflip", "imag")]->Eval(cth);
            double g_real = Graphs_Local[key(component, "flip", "real")]->Eval(cth);
            double g_imag = Graphs_Local[key(component, "flip", "imag")]->Eval(cth);
            outfile_angular<<cth<<" "<<f_real<<" "<<f_imag<<" "<<g_real<<" "<<g_imag<<endl;
        }
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                TString graph_key_local = key(component, flip_nonflip, real_imag);
                TGraph* graph_local = Graphs_Local[graph_key_local];
                for(int ip = 0; ip < graph_local->GetN(); ++ip){
                    double cth, value;
                    graph_local->GetPoint(ip, cth, value);
                    y_min = std::min(y_min, value);
                    y_max = std::max(y_max, value);   
                }
                y_range = y_max - y_min;
            }
        }
        gPad->DrawFrame(-1.0, y_min - 0.1*y_range, 1.0, y_max + 0.1*y_range);
        for(const auto& flip_nonflip : NF_F){
            for(const auto& real_imag : Re_Im){
                TString graph_key_local = key(component, flip_nonflip, real_imag);
                TGraph* graph_local = Graphs_Local[graph_key_local];
                graph_local->Draw("L SAME"); 
                legend->AddEntry(graph_local, Form("%s %s", flip_nonflip.Data(), real_imag.Data()), "l"); 
            }
        }
        legend->Draw();
        canvas->SaveAs(figdir + "Angular/" + Form("/Angular_%s.pdf", component.Data()));
    }
}


void DrawComponents(){
    SetPaperStyle();
    cout<<"Drawing components..."<<endl;
#if ModelConfig == 0
    figdir ="figs/Default/";
#elif ModelConfig == 1
    figdir ="figs/PolaFit/";
#elif ModelConfig == 2:
    figdir ="figs/Lambda2110_Inv/";
#endif
    gSystem->mkdir(figdir, true);
    RPRModel model;
    if(HelicityBasis){
        model.SetHelicityBasis(1);
        figdir ="figs/PolaFit_HelicityBasis/";
    }
#if ModelConfig == 0
#elif ModelConfig == 1
    ConfigureModel(model);
#elif ModelConfig == 2
    ConfigureModelLambda2110(model);
#endif
    cout<<"Model configured."<<endl;
    cout<<"Drawing with theory..."<<endl;
    DrawWithTheory(model);
    DrawAngularComponents(model, 2.3);
    cout<<"Done drawing components."<<endl;
    vector<std::string> components = {
        "total", "uch_regge", "Sigma2030", "Lambda2100", "Sigma2250","Lambda1890", "s_lambda", "s_sigma"
        //"total", "uch_regge", "Sigma2030", "Lambda2100", "Sigma2250"
        //"total", "Sigma2030", "Lambda2100", "Sigma2250"
    };
    #if ModelConfig == 2
    components.push_back("Lambda2110");
    #endif
    vector<array<string,2>> fg_combi = {
        {"total", "total"},
        {"uch_regge", "Sigma2030"},
        {"uch_regge", "Lambda2100"},
        {"uch_regge", "Sigma2250"},
        {"Sigma2030", "Lambda2100"}
    };
    vector<int> colors = 
    {
        kBlack, kCyan + 2,kMagenta, kBlue , kRed ,  kGreen + 2, kOrange + 7, kGray+2
    };
    vector< std::array<double,2> > angularBins = {
        {{0.85, 0.90}}, {{0.90, 0.97}}, {{0.97, 1.00}},{{0.85, 1.00}}
    };
    map<string, TGraph*> Graphs;
    //double sqrtS_low = 1.85, sqrtS_high = 2.4, dS = 0.01;
    double sqrtS_low = 1.83, sqrtS_high = 2.61, dS = 0.01;
    #if 1
    for(const auto& bin : angularBins){
        cout<<"Processing angular bin: "<<bin[0]<<" < cosTheta < "<<bin[1]<<endl;
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            gSystem->mkdir(figdir +Form("Amplitudes/%s/",components[ic].c_str()), true);
            ofstream FComp_Real(figdir +Form("Amplitudes/%s/",components[ic].c_str())+ Form("nonflip_real_%s_cth_%g_%g.dat", components[ic].c_str(),bin[0]*100, bin[1]*100));
            ofstream FComp_Imag(figdir +Form("Amplitudes/%s/",components[ic].c_str())+ Form("nonflip_imag_%s_cth_%g_%g.dat", components[ic].c_str(),bin[0]*100, bin[1]*100));
            ofstream GComp_Real(figdir +Form("Amplitudes/%s/",components[ic].c_str())+ Form("flip_real_%s_cth_%g_%g.dat", components[ic].c_str(),bin[0]*100, bin[1]*100));
            ofstream GComp_Imag(figdir +Form("Amplitudes/%s/",components[ic].c_str())+ Form("flip_imag_%s_cth_%g_%g.dat", components[ic].c_str(),bin[0]*100, bin[1]*100));
            cout<<"Processing component: "<<components[ic]<<endl;
            const auto& component = components[ic];
            double cosThetaLow = bin[0], cosThetaHigh = bin[1];
            string key = GetKey(component, 0, 0, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kSolid);
            key = GetKey(component, 0, 1, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kDashed);
            key = GetKey(component, 1, 0, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kSolid);
            key = GetKey(component, 1, 1, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kDashed);


            int is = 0;
            for(double sqrtS_0 = sqrtS_low; sqrtS_0 <= sqrtS_high; sqrtS_0 += dS){
                double sqrtS = sqrtS_0 + dS < sqrtS_high ?
                sqrtS_0 + 0.5*dS : 0.5*(sqrtS_0 + sqrtS_high);
                if(is%10 == 0){
                    cout << "Processing sqrtS = " << sqrtS << " GeV for component " << component << endl;
                }
                is++;
                const auto amplitudes = GetAmplitudes(
                    model, sqrtS, cosThetaLow, cosThetaHigh, &component);
                string key = GetKey(component, 0, 0, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitudes.f.real());
                key = GetKey(component, 0, 1, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitudes.f.imag());
                key = GetKey(component, 1, 0, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitudes.g.real());
                key = GetKey(component, 1, 1, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitudes.g.imag());
                double ampl = amplitudes.fIntensity + amplitudes.gIntensity;
                ampl = sqrt(ampl);
                FComp_Real << sqrtS << " " << amplitudes.f.real() << endl;
                FComp_Imag << sqrtS << " " << amplitudes.f.imag() << endl;
                GComp_Real << sqrtS << " " << amplitudes.g.real() << endl;
                GComp_Imag << sqrtS << " " << amplitudes.g.imag() << endl;
            }
        }
        TString ct_nonflip = Form("cth_%g_%g_nonflip", bin[0]*100, bin[1]*100);
        TString ct_flip = Form("cth_%g_%g_flip", bin[0]*100, bin[1]*100);
        double f_max = 0.0, g_max = 0.0;
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            const auto& component = components[ic];
            string key = GetKey(component, 0, 0, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                f_max = std::max(f_max, std::abs(Graphs[key]->GetY()[i]));
            }
            key = GetKey(component, 0, 1, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                f_max = std::max(f_max, std::abs(Graphs[key]->GetY()[i]));
            }
            key = GetKey(component, 1, 0, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                g_max = std::max(g_max, std::abs(Graphs[key]->GetY()[i]));
            }
            key = GetKey(component, 1, 1, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                g_max = std::max(g_max, std::abs(Graphs[key]->GetY()[i]));
            }
        }
        f_max *= 1.2;
        g_max *= 1.2;
        TCanvas* c_nonflip = new TCanvas(ct_nonflip, ct_nonflip, 2000,1000);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            const auto& component = components[ic];
            string key = GetKey(component, 0, 0, bin[0], bin[1]);
            Graphs[key]->SetTitle(Form("Non-flip amplitudes: %.2f < cos#theta < %.2f", bin[0], bin[1]));
            Graphs[key]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
            Graphs[key]->GetYaxis()->SetRangeUser(-f_max, f_max);
            Graphs[key]->Draw(ic == 0 ? "AL" : "L SAME");
            key = GetKey(component, 0, 1, bin[0], bin[1]);
            Graphs[key]->Draw("L SAME");
        }
        TLegend* legend = new TLegend(0.25, 0.7, 0.5, 0.9);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            const auto& component = components[ic];
            string key = GetKey(component, 0, 0, bin[0], bin[1]);
            legend->AddEntry(Graphs[key], component.c_str(), "l"); 
        }
        legend->Draw();
        c_nonflip->SaveAs(figdir + Form("nonflip_amplitudes_cth_%g_%g.pdf", bin[0]*100, bin[1]*100));
        TCanvas* c_flip = new TCanvas(ct_flip, ct_flip, 2000,1000);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            const auto& component = components[ic];
            string key = GetKey(component, 1, 0, bin[0], bin[1]);
            Graphs[key]->SetTitle(Form("Flip amplitudes: %.2f < cos#theta < %.2f", bin[0], bin[1]));
            Graphs[key]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
            Graphs[key]->GetYaxis()->SetRangeUser(-g_max, g_max);
            Graphs[key]->Draw(ic == 0 ? "AL" : "L SAME");
            key = GetKey(component, 1, 1, bin[0], bin[1]);
            Graphs[key]->Draw("L SAME");
        }
        TLegend* legend_flip = new TLegend(0.25, 0.7, 0.5, 0.9);
        legend_flip->SetBorderSize(0);
        legend_flip->SetFillStyle(0);
        for(std::size_t ic = 0; ic < components.size(); ++ic){
            const auto& component = components[ic];
            string key = GetKey(component, 1, 0, bin[0], bin[1]);
            legend_flip->AddEntry(Graphs[key], component.c_str(), "l");
        }
        legend_flip->Draw();
        c_flip->SaveAs(figdir + Form("flip_amplitudes_cth_%g_%g.pdf", bin[0]*100, bin[1]*100));
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            double cosThetaLow = bin[0], cosThetaHigh = bin[1];
            string key = GetPairedKey(combi[0], combi[1], 0, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kDashed);
            key = GetPairedKey(combi[0], combi[1], 1, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            Graphs[key]->SetLineStyle(kSolid);
            key = GetPairedKey(combi[0], combi[1], 2, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            key = GetPairedKey(combi[0], combi[1], 3, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            key = GetPairedKey(combi[0], combi[1], 4, cosThetaLow, cosThetaHigh);
            Graphs[key] = new TGraph();
            Graphs[key]->SetLineColor(colors[ic]);
            int is = 0;
            for(double sqrtS_0 = sqrtS_low; sqrtS_0 <= sqrtS_high; sqrtS_0 += dS){
                double sqrtS = sqrtS_0 + dS < sqrtS_high ?
                sqrtS_0 + 0.5*dS : 0.5*(sqrtS_0 + sqrtS_high);
                is++;
                if(is%10 == 0){
                    cout << "Processing sqrtS = " << sqrtS << " GeV for component pair " << combi[0] << " and " << combi[1] << endl;
                }
                const auto amplitudes1 = GetAmplitudes(
                    model, sqrtS, cosThetaLow, cosThetaHigh, &combi[0]);
                const auto amplitudes2 = GetAmplitudes(
                    model, sqrtS, cosThetaLow, cosThetaHigh, &combi[1]);
                auto amplitude_sum = amplitudes1 + amplitudes2;
                if(combi[0] == "total" and combi[1] == "total"){
                    amplitude_sum = amplitudes1;
                }
                key = GetPairedKey(combi[0], combi[1], 0, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, 2 * amplitude_sum.fConjugateG.imag());
                key = GetPairedKey(combi[0], combi[1], 1, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitude_sum.fIntensity + amplitude_sum.gIntensity);
                key = GetPairedKey(combi[0], combi[1], 2, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, 2 * amplitude_sum.fConjugateG.imag() / (amplitude_sum.fIntensity + amplitude_sum.gIntensity));
                key = GetPairedKey(combi[0], combi[1], 3, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitude_sum.fIntensity);
                key = GetPairedKey(combi[0], combi[1], 4, cosThetaLow, cosThetaHigh);
                Graphs[key]->SetPoint(Graphs[key]->GetN(), sqrtS, amplitude_sum.gIntensity);
            }
        }
        TString ct_ampl = Form("cth_%g_%g_ampl", bin[0]*100, bin[1]*100);
        TCanvas* c_ampl = new TCanvas(ct_ampl, ct_ampl, 2000,1000);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        double ampl_max = 0.0;
        double inter_max = 0, inter_min = 0;
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 1, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                ampl_max = std::max(ampl_max, std::abs(Graphs[key]->GetY()[i]));
            }
            key = GetPairedKey(combi[0], combi[1], 0, bin[0], bin[1]);
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                inter_max = std::max(inter_max, Graphs[key]->GetY()[i]);
                inter_min = std::min(inter_min, Graphs[key]->GetY()[i]);
            }
        }
        ampl_max *= 1.2;
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 1, bin[0], bin[1]);
            Graphs[key]->SetTitle(Form("Interference and Polarization: %.2f < cos#theta < %.2f", bin[0], bin[1]));
            Graphs[key]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
            Graphs[key]->GetYaxis()->SetRangeUser(-0.1*ampl_max, ampl_max);
            Graphs[key]->Draw(ic == 0 ? "AL" : "L SAME");
            key = GetPairedKey(combi[0], combi[1], 3, bin[0], bin[1]);
            Graphs[key]->SetLineStyle(kDashed);
            Graphs[key]->Draw("L SAME");
            key = GetPairedKey(combi[0], combi[1], 4, bin[0], bin[1]);
            Graphs[key]->SetLineStyle(5);
            Graphs[key]->Draw("L SAME");
        }
        TLegend* legend_ampl = new TLegend(0.15, 0.7, 0.5, 0.9);
        legend_ampl->SetBorderSize(0);
        legend_ampl->SetFillStyle(0);
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 1, bin[0], bin[1]);
            legend_ampl->AddEntry(Graphs[key], Form("%s + %s", combi[0].c_str(), combi[1].c_str()), "l");
            if(ic == 0){
                key = GetPairedKey(combi[0], combi[1], 3, bin[0], bin[1]);
                legend_ampl->AddEntry(Graphs[key], "nonflip^{2}", "l");
                key = GetPairedKey(combi[0], combi[1], 4, bin[0], bin[1]);
                legend_ampl->AddEntry(Graphs[key], "flip^{2}", "l");
            }
        }
        legend_ampl->Draw();
        c_ampl->SaveAs(figdir + Form("amplitude2_cth_%g_%g.pdf", bin[0]*100, bin[1]*100));
        
        TString ct_inter = Form("cth_%g_%g_inter", bin[0]*100, bin[1]*100);
        TCanvas* c_inter = new TCanvas(ct_inter, ct_inter, 2000,1000);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        TLegend* legend_inter = new TLegend(0.15, 0.7, 0.5, 0.9);
        legend_inter->SetBorderSize(0);
        legend_inter->SetFillStyle(0);
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 0, bin[0], bin[1]);
            TGraph* gr = (TGraph*)Graphs[key]->Clone(Form("gr_%s_%s", combi[0].c_str(), combi[1].c_str()));
            gr->SetTitle(Form("Interference and Polarization: %.2f < cos#theta < %.2f", bin[0], bin[1]));
            gr->SetLineStyle(kSolid);
            gr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
            gr->GetYaxis()->SetRangeUser(1.2*inter_min, 1.2*inter_max);
            gr->Draw(ic == 0 ? "AL" : "L SAME");
            legend_inter->AddEntry(gr, Form("%s + %s", combi[0].c_str(), combi[1].c_str()), "l");
        }
        legend_inter->Draw();
        c_inter->SaveAs(figdir + Form("interference_cth_%g_%g.pdf", bin[0]*100, bin[1]*100));
        TString ct_pol = Form("cth_%g_%g_pol", bin[0]*100, bin[1]*100);
        TCanvas* c_pol = new TCanvas(ct_pol, ct_pol, 2000 ,1000);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.05);
        map<string, int> res_color = {
            {"Sigma2030", kRed},
            {"Lambda2100", kBlue},
            {"Sigma2250", kGreen + 2},
            {"total", kBlack}
        };
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 2, bin[0], bin[1]);
            if(combi[0] == "uch_regge" or combi[1] == "uch_regge") Graphs[key]->SetLineStyle(kDashed);
            else Graphs[key]->SetLineStyle(kSolid);
            Graphs[key]->SetLineColor(res_color[combi[1]]);
            Graphs[key]->SetTitle(Form("Polarization: %.2f < cos#theta < %.2f", bin[0], bin[1]));
            Graphs[key]->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
            Graphs[key]->GetYaxis()->SetRangeUser(-1.2, 1.2);
            Graphs[key]->Draw(ic == 0 ? "AL" : "L SAME");
            ofstream outfile(figdir + Form("PXi_%s_%s_cth_%g_%g.txt"
            , combi[0].c_str(), combi[1].c_str(), bin[0]*100, bin[1]*100));
            for(int i = 0; i < Graphs[key]->GetN(); i++){
                outfile << Graphs[key]->GetX()[i] << " " << Graphs[key]->GetY()[i] << endl;
            }
        }
        TLegend* legend_pol = new TLegend(0.15, 0.7, 0.5, 0.9);
        legend_pol->SetBorderSize(0);
        legend_pol->SetFillStyle(0);
        for(std::size_t ic = 0; ic < fg_combi.size(); ++ic){
            const auto& combi = fg_combi[ic];
            string key = GetPairedKey(combi[0], combi[1], 2, bin[0], bin[1]);
            legend_pol->AddEntry(Graphs[key], Form("%s + %s", combi[0].c_str(), combi[1].c_str()), "l");
        }
        legend_pol->Draw();
        c_pol->SaveAs(figdir + Form("polarization_cth_%g_%g.pdf", bin[0]*100, bin[1]*100));
    }
    #endif
    double sqrtS_ref = 2.0;
    double cth_low = -0.9, cth_high = 0.9, dcth = 0.1;
    TGraph* g_regge_ampl = new TGraph();
    string channel = "uch_regge";
    for(double cth_0 = cth_low; cth_0 <= cth_high; cth_0 += dcth){
        double cth = cth_0 + dcth < cth_high ?
        cth_0 + 0.5*dcth : 0.5*(cth_0 + cth_high);
        cout << "cth = " << cth;
        const auto amplitudes = GetAmplitudes(model, sqrtS_ref, cth, cth + dcth, &channel);
        cout<< ", fIntensity = " << amplitudes.fIntensity << ", gIntensity = " << amplitudes.gIntensity << endl;
        g_regge_ampl->SetPoint(g_regge_ampl->GetN(), cth, amplitudes.fIntensity + amplitudes.gIntensity);
    }
    TCanvas* c_regge_ampl = new TCanvas("regge_ampl", "regge_ampl", 2000,1000);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    g_regge_ampl->SetTitle(Form("Regge amplitude: #sqrt{s} = %.2f GeV", sqrtS_ref));
    g_regge_ampl->GetXaxis()->SetTitle("cos#theta");
    g_regge_ampl->GetYaxis()->SetTitle("Intensity");
    g_regge_ampl->Draw("AL");
    c_regge_ampl->SaveAs(figdir + Form("regge_amplitude_sqrtS_%g.pdf", sqrtS_ref));
    
    gE42Data = LoadE42Data(
        "Datapoints/Polarization/E42_polarization_CM.tsv");
    gForwardData = LoadForwardData(
        "Datapoints/ForwardCrossSection/Forward_differential_cross_sections.tsv");

    DrawWithData(model, gForwardData);
    cout << "Finished drawing components:" << endl;
    cout << Form(" Chi2 / Ndf = %.2f / %d, %.2f /%d, %.2f /%d",
         gChi2Pol, gnpointsPol, gChi2Forward, gnpointsForward, gChi2, gnpoints) << endl;
}
