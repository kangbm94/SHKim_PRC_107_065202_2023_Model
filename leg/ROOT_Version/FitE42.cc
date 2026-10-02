#include "FitE42.hh"
#include <TMinuit.h>
#include "../../Include/GraphReader.hh"
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>

#include <iomanip>
#include <iostream>
#include <vector>

using namespace rpr;
diff_files = {
    base_dir + "Datapoints/DiffCrossSection/Berge_differential_cross_section.tsv",
    base_dir + "Datapoints/DiffCrossSection/Burgun_differential_cross_section.tsv",
    base_dir + "Datapoints/DiffCrossSection/Dauber_differential_cross_section.tsv",
    //base_dir + "Datapoints/DiffCrossSection/Baltay_differential_cross_section.tsv",
    base_dir + "Datapoints/DiffCrossSection/Trippe_differential_cross_section.tsv",
    base_dir + "Datapoints/DiffCrossSection/E05_dif_CM.tsv",
};
diff_names = {
    "Berge",
    "Burgun",
    "Dauber",
    //"Baltay",
    "Trippe",
    "E05"
};

void FitE42(){
    gS2030 = 1.;
    gL2100 = 1.;
    gS2250 = 1.;
    gS2230 = 1.;
    cth0s = {0.85,0.9,0.97};
    cth1s = {0.9,0.97,1.0};
    vector<int> colors = {kBlue, kRed, kBlack};
    TGraphErrors* g_e42[cth0s.size()];
    for(int ith = 0; ith<cth0s.size(); ith++){
        TString filename = E42pol_name(cth0s[ith], cth1s[ith]);
        vector<double> x_e42, y_e42, yerr_e42;
        GraphReader::LoadE42Polarization(filename, x_e42, y_e42, yerr_e42);
        if(ith == 0) sqrtSs = x_e42;
        E42pol.push_back(y_e42);
        E42polErr.push_back(yerr_e42);
        g_e42[ith] = new TGraphErrors(x_e42.size(), &x_e42[0], &y_e42[0], nullptr, &yerr_e42[0]);
        g_e42[ith]->SetMarkerStyle(20);
        g_e42[ith]->SetMarkerColor(colors[ith]);
        g_e42[ith]->SetLineColor(colors[ith]);
        for(int ipt=0; ipt<x_e42.size(); ipt++){
            cout<<Form("E42 point: cth [%g,%g],  sqrtS = %g, pol = %g +/- %g",
                 cth0s[ith], cth1s[ith], x_e42[ipt], y_e42[ipt], yerr_e42[ipt])<<endl;
        }
    }
    DiffCrossSection diffs;
    LoadDiffCrossSections(diff_files, diff_names, diffs);
    SortDiffCrossSection(diffs);
    TString fwdData = base_dir + "Datapoints/ForwardCrossSection/Forward_differential_cross_sections.tsv";
    LoadFwdData(fwdData);
    TMinuit* minuit = new TMinuit(3);
    //minuit->SetFCN(Chi2E42Pol);
    //minuit->SetFCN(Chi2E42DifPol);
    minuit->SetFCN(Chi2E42FwdDifPol);
    minuit->SetPrintLevel(-1);
    int ierflg = 0;
    minuit->mnparm(0, "S2030", 1.0, 0.1, 0.4, 1.2, ierflg);
    minuit->mnparm(1, "L2100", 1.0, 0.1, 0.4, 1.2, ierflg);
    minuit->mnparm(2, "S2250", 1.0, 0.1, 0.4, 1.5, ierflg);
    minuit->mnparm(3, "S2230", 1.0, 0.1, 0.4, 1.2, ierflg);
    //minuit->Migrad();
    double S2030_err, L2100_err, S2250_err, S2230_err;
    minuit->GetParameter(0, gS2030, S2030_err);
    minuit->GetParameter(1, gL2100, L2100_err);
    minuit->GetParameter(2, gS2250, S2250_err);
    minuit->GetParameter(3, gS2230, S2230_err);
    std::cout<<"Fit Result: S2030 = "<<gS2030<<"+/-"<<S2030_err<<", L2100 = "<<gL2100<<"+/-"<<L2100_err<<", S2250 = "<<gS2250<<"+/-"<<S2250_err<<", S2230 = "<<gS2230<<"+/-"<<S2230_err<<std::endl;
    ofstream result_file("FitE42_result.txt");
    result_file<<"Fit Result: S2030 = "<<gS2030<<"+/-"<<S2030_err<<", L2100 = "<<gL2100<<"+/-"<<L2100_err<<", S2250 = "<<gS2250<<"+/-"<<S2250_err<<", S2230 = "<<gS2230<<"+/-"<<S2230_err<<std::endl;
    result_file.close();
    RPRModel model;
    gS2030 = 1.;
    gL2100 = 0.75;
    gS2250 = 1.;
    gS2230 = 0.;
    guLambda = 1.;
    SetModelParameters(model);

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    double sqrtS_0 = 1.85;
    double sqrtS_1 = 2.5;
    double d_sqrtS = 0.01;
    
    for(int ith = 0; ith<cth0s.size(); ith++){
        TGraph* g_model = new TGraph();
        g_model->SetLineColor(colors[ith]);
        g_model->SetLineWidth(2);
        double d_cth = 0.001;
        for(double sqrtS=sqrtS_0; sqrtS<sqrtS_1; sqrtS+=d_sqrtS){
            double pol = 0;
            double dsig = 0;
            for(double cth=cth0s[ith]; cth<cth1s[ith]; cth+=d_cth){
                //cout<<Form("Model: cth [%g,%g],  sqrtS = %g, cth = %g", cth0s[ith], cth1s[ith], sqrtS, cth)<<endl;
                ModelResult value = model.Evaluate(sqrtS, cth);
                pol += value.polarization*value.differentialCrossSection*d_cth;
                dsig += value.differentialCrossSection*d_cth;
            }
            pol /= dsig;
            g_model->SetPoint(g_model->GetN(), sqrtS, pol);
        }
        g_e42[ith]->GetXaxis()->SetTitle("sqrt(s) [GeV]");
        g_e42[ith]->GetXaxis()->SetLimits(sqrtS_0, sqrtS_1);
        g_e42[ith]->GetYaxis()->SetTitle("Polarization");   
        g_e42[ith]->GetYaxis()->SetLimits(-1, 1);
        g_e42[ith]->GetYaxis()->SetRangeUser(-1, 1);
        ith == 0? g_e42[ith]->Draw("APE"): g_e42[ith]->Draw("PEsame");
        g_model->Draw("Lsame");
    }
    TLegend* legend = new TLegend(0.2, 0.7, 0.5, 0.9);
    legend->AddEntry(g_e42[0], "0.85 < cos#theta < 0.9", "p");
    legend->AddEntry(g_e42[1], "0.9 < cos#theta < 0.97", "p");
    legend->AddEntry(g_e42[2], "0.97 < cos#theta < 1.0", "p");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("E42_fit_result.pdf");   
    vector<TString> components = {"Sigma2030", "Lambda2100", "Sigma2250", "Lambda1890"};
    vector<int> comp_colors = {kGreen+4, kRed, kBlue, kOrange+7};



    for(auto d:DiffDatasets){
        if(d.second.sqrtS < 1.9 or d.second.sqrtS > 2.5) continue;
        TCanvas* c_diff = new TCanvas(Form("c_diff_%s", d.first.Data()), Form("c_diff_%s", d.first.Data()), 800, 600);
        TGraphErrors* g_diff = new TGraphErrors(d.second.cth.size(), &d.second.cth[0], &d.second.dsig[0], nullptr, &d.second.dsigErr[0]);
        double maxi = 0;
        for(int ipt=0; ipt<d.second.cth.size(); ipt++){
            if(d.second.dsig[ipt] > maxi) maxi = d.second.dsig[ipt];
        }
        g_diff->SetMarkerStyle(20);
        g_diff->SetMarkerColor(kBlue);
        g_diff->SetLineColor(kBlue);
        g_diff->GetXaxis()->SetTitle("cos#theta");
        g_diff->GetYaxis()->SetTitle("d#sigma/d#Omega [#mu b/sr]");
        g_diff->GetYaxis()->SetRangeUser(0, maxi*2);
        g_diff->GetXaxis()->SetLimits(0.5, 1);
        g_diff->Draw("APE");
        TGraph* g_model_diff = new TGraph();
        g_model_diff->SetLineColor(kBlack);
        g_model_diff->SetLineWidth(2);
        map<TString, TGraph*> model_diff_components;
        double sqrtS = d.second.sqrtS;
        TLegend* legend_diff = new TLegend(0.2, 0.7, 0.5, 0.9);
        legend_diff->AddEntry(g_diff, "Data", "p");
        legend_diff->SetBorderSize(0);
        legend_diff->SetFillStyle(0);
        for(double cth=-1; cth<1; cth+=0.01){
            ModelResult value = model.Evaluate(sqrtS, cth);
            g_model_diff->SetPoint(g_model_diff->GetN(), cth, value.differentialCrossSection);
        }
        for(int icomp = 0; icomp<components.size(); icomp++){
            TString comp = components[icomp];
            model_diff_components[comp] = new TGraph();
            model_diff_components[comp]->SetLineColor(comp_colors[icomp]);
            model_diff_components[comp]->SetLineStyle(kDashed);
            model_diff_components[comp]->SetLineWidth(2);
            legend_diff->AddEntry(model_diff_components[comp], comp.Data(), "l");
            for(double cth=-1; cth<1; cth+=0.01){
                ModelResult value_comp = model.Evaluate(sqrtS, cth, comp.Data());
                model_diff_components[comp]->SetPoint(model_diff_components[comp]->GetN(), cth, value_comp.differentialCrossSection);
            }
        }
        legend_diff->Draw();
        g_model_diff->Draw("Lsame");
        for(int icomp = 0; icomp<components.size(); icomp++){
            TString comp = components[icomp];
            model_diff_components[comp]->Draw("Lsame");
        }
        TLatex* latex = new TLatex();
        latex->SetTextSize(0.04);
        latex->DrawLatexNDC(0.55, 0.85, Form("sqrt(s) = %.2f GeV", sqrtS));
        c_diff->Update();
        c_diff->Modified();
        c_diff->SaveAs(Form("DiffCrossSection_fit_result_%s.pdf", d.first.Data()));
    }
    TCanvas* c_FwdDiff = new TCanvas("c_FwdDiff", "c_FwdDiff", 800, 600);
    //TGraphErrors* g_FwdDiff = new TGraphErrors(OldFwdData.name.size(), &OldFwdData.sqrts[0], &OldFwdData.dsig[0], nullptr, &OldFwdData.dsigerr[0]);
    TGraphErrors* g_FwdDiff = new TGraphErrors(OldFwdData.name.size(), &OldFwdData.beam_mom[0], &OldFwdData.dsig[0], nullptr, &OldFwdData.dsigerr[0]);
    TLegend* legend_FwdDiff = new TLegend(0.7, 0.4, 0.9, 0.9);
    g_FwdDiff->SetMarkerStyle(20);
    g_FwdDiff->SetMarkerColor(kBlue);
    g_FwdDiff->SetLineColor(kBlue);
    //g_FwdDiff->GetXaxis()->SetTitle("sqrt(s) [GeV]");
    g_FwdDiff->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
    g_FwdDiff->GetYaxis()->SetTitle("d#sigma/d#Omega [#mu b/sr]");
    legend_FwdDiff->AddEntry(g_FwdDiff, "Data", "p");
    legend_FwdDiff->SetBorderSize(0);
    legend_FwdDiff->SetFillStyle(0);
    g_FwdDiff->Draw("APE");
    c_FwdDiff->Update();
    TGraph* g_model_FwdDiff = new TGraph();
    g_model_FwdDiff->SetLineColor(kBlack);
    g_model_FwdDiff->SetLineWidth(2);
    map<TString, TGraph*> model_FwdDiff_components;
    for(int icomp = 0; icomp<components.size(); icomp++){
        TString comp = components[icomp];
        model_FwdDiff_components[comp] = new TGraph();
        model_FwdDiff_components[comp]->SetLineColor(comp_colors[icomp]);
        model_FwdDiff_components[comp]->SetLineStyle(kDashed);
        model_FwdDiff_components[comp]->SetLineWidth(2);
        legend_FwdDiff->AddEntry(model_FwdDiff_components[comp], comp.Data(), "l");
    }
    for(double sqrtS=sqrtS_0; sqrtS<sqrtS_1; sqrtS+=d_sqrtS){
        double model_sig = model.AverageDifferentialCrossSectionLabTheta(sqrtS, 0.0,20.0);
        double pK = CalculateBeamMomentum(sqrtS);
        g_model_FwdDiff->SetPoint(g_model_FwdDiff->GetN(), pK, model_sig);
        for(auto& comp:components){
            double comp_sig = model.AverageDifferentialCrossSectionLabTheta(sqrtS, 0.0,20.0, comp.Data());
            model_FwdDiff_components[comp]->SetPoint(model_FwdDiff_components[comp]->GetN(), pK, comp_sig);
        }
        cout<<Form("Model FwdDiff: sqrtS = %g, sig = %g", sqrtS, model_sig)<<endl;
    }
    g_model_FwdDiff->Draw("Lsame");
    for(int icomp = 0; icomp<components.size(); icomp++){
        TString comp = components[icomp];
        model_FwdDiff_components[comp]->Draw("Lsame");
    }
    legend_FwdDiff->Draw();
    c_FwdDiff->Update();
    c_FwdDiff->Modified();
    c_FwdDiff->SaveAs("FwdDiff_fit_result.pdf");
}
