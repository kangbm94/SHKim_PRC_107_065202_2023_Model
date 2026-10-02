#include "../../Include/Dir.hh"
#include "FitE42.hh"
std::vector<E42Point> gE42Data;
std::vector<ForwardPoint> gForwardData;
double gChi2Pol = 0.0;
double gChi2Forward = 0.0;
double gChi2 = 0.0;
int gnpointsPol = 0;
int gnpointsForward = 0;
int gnpoints = 0;
double crosssection_min = 10;
constexpr int kForwardFitIntervals = 16;
constexpr double kNodes[12] = {
    -0.9815606342467192, -0.9041172563704749,
    -0.7699026741943047, -0.5873179542866175,
    -0.3678314989981802, -0.1252334085114689,
     0.1252334085114689,  0.3678314989981802,
     0.5873179542866175,  0.7699026741943047,
     0.9041172563704749,  0.9815606342467192
};
constexpr double kWeights[12] = {
    0.0471753363865118, 0.1069393259953184,
    0.1600783285433462, 0.2031674267230659,
    0.2334925365383548, 0.2491470458134029,
    0.2491470458134029, 0.2334925365383548,
    0.2031674267230659, 0.1600783285433462,
    0.1069393259953184, 0.0471753363865118
};
std::vector<E42Point> LoadE42Data(const std::string& filename)
{
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

std::vector<ForwardPoint> LoadForwardData(const std::string& filename)
{
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
    model.SetEnabled("Sigma2230", 0);
    auto& sigma2030 = model.Get("Sigma2030");
    sigma2030.amplitudeScale = Complex(0.55, 0.0);
    sigma2030.mass = 2.040;
    sigma2030.width = 0.200;
    auto& lambda2100 = model.Get("Lambda2100");
    lambda2100.mass = 2.090;
    lambda2100.width = 0.150;
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
                    //const std::vector<DifferentialDataset>& differentialData){
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
    canvasPol->SaveAs("figs/Polarization_wData.pdf");
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
    for (double pLab = 1.15; pLab <= 2.72; pLab += 0.015) {
        const double eK = std::sqrt(pLab*pLab + mK*mK);
        const double sqrtS = std::sqrt(mK*mK + mP*mP + 2.0*mP*eK);
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
    canvasCrossSection->SaveAs("figs/CrossSection_wData.pdf");

}
void DrawComponents(){
    SetPaperStyle();
    TString figdir = "figs/";
    gSystem->mkdir(figdir, true);
    RPRModel model;
    ConfigureModel(model);
    vector<std::string> components = {
        "total", "born_term", "uch_regge", "Sigma2030", "Lambda2100", "Sigma2250","Lambda1890"
    };
    vector<array<string,2>> fg_combi = {
 //       {"born_term", "Sigma2030"},
 //       {"born_term", "Lambda2100"},
 //       {"born_term", "Sigma2250"},
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
        ofstream F_Real(figdir + Form("nonflip_real_total_cth_%g_%g.dat", bin[0]*100, bin[1]*100));
        ofstream F_Imag(figdir + Form("nonflip_imag_total_cth_%g_%g.dat", bin[0]*100, bin[1]*100));
        ofstream G_Real(figdir + Form("flip_real_total_cth_%g_%g.dat", bin[0]*100, bin[1]*100));
        ofstream G_Imag(figdir + Form("flip_imag_total_cth_%g_%g.dat", bin[0]*100, bin[1]*100));
        for(std::size_t ic = 0; ic < components.size(); ++ic){
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
                is++;
                if(is%10 == 0){
                    cout << "Processing sqrtS = " << sqrtS << " GeV for component " << component << endl;
                }
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
                if(component == "total"){
                    double ampl = amplitudes.fIntensity + amplitudes.gIntensity;
                    ampl = sqrt(ampl);
                    F_Real << sqrtS << " " << amplitudes.f.real() << endl;
                    F_Imag << sqrtS << " " << amplitudes.f.imag() << endl;
                    G_Real << sqrtS << " " << amplitudes.g.real() << endl;
                    G_Imag << sqrtS << " " << amplitudes.g.imag() << endl;
                }
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
