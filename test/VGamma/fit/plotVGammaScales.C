#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooNumIntConfig.h"
#include "TROOT.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/tdrstyle.C"

using namespace std;

struct GaussResult {
    double mu;
    double err;
    double chi2ndf;
};

GaussResult extractGaussFromTxt(const char* filename, bool isQCD, int RunYear,
                              int preVFP,
                              int ichannel,
                              double lowpt,
                              double highpt,
                              double lowmet,
                              double highmet)
{

ifstream file(filename);
float leplow, lephigh;
float fakescale, fakescaleerror;
float vgammascale, vgammascaleerror;
vector<double> values;
double fitError = 0;
    bool firstLine=true;

    while(firstLine == true){

        file >> leplow >> lephigh
             >> fakescale >> fakescaleerror
             >> vgammascale >> vgammascaleerror;

        if(!file.good()) break;

        float val  = isQCD ? fakescale      : vgammascale;
        float err  = isQCD ? fakescaleerror : vgammascaleerror;

        if(!std::isfinite(val)) continue;
        if(val<=0) continue;

        values.push_back(val);

        if(firstLine){
            fitError = err;
            firstLine=false;
        }
    }
for(int i=0;i<1000;i++){

    file >> leplow >> lephigh
         >> fakescale >> fakescaleerror
         >> vgammascale >> vgammascaleerror;

    if(!file.good()) break;

    float val = isQCD ? fakescale : vgammascale;
    values.push_back(val);
}

sort(values.begin(), values.end());
double lowbound  = values.front();
double highbound = values.back();

int n = values.size();

double q05 = values[int(0.05*n)];
double q95 = values[int(0.95*n)];

double width  = q95 - q05;
double margin = 0.30 * width;

double lowRange  = q05 - margin;
double highRange = q95 + margin;

if(highRange<=lowRange){
    lowRange  = values.front() - 0.1*fabs(values.front());
    highRange = values.back()  + 0.1*fabs(values.back());
}

int nbins = max(30, int(2*sqrt(n)));

TH1F *h = new TH1F(
    Form("h_%.0f_%.0f_%.0f_%.0f",lowpt,highpt,lowmet,highmet),
    "",
    nbins,
    lowRange,
    highRange
);

for(auto v : values)
    h->Fill(v);

while(true){
    int zeroBins=0;
    for(int i=1;i<=h->GetNbinsX();i++)
        if(h->GetBinContent(i)==0) zeroBins++;

    double frac = (double)zeroBins/h->GetNbinsX();

    if(frac < 0.15 || h->GetNbinsX()<15) break;

    h->Rebin(2);}

/*int zeroBins=0;
for(int i=1;i<=h->GetNbinsX();i++)
    if(h->GetBinContent(i)==0) zeroBins++;

if(zeroBins > 0.3*h->GetNbinsX())
    h->Rebin(2);*/

    double mean = h->GetMean();
    double rms = h->GetRMS();
    double fitLow  = mean - 1.2*rms;
    double fitHigh = mean + 1.2*rms;

    h->Fit("gaus","Q","",fitLow,fitHigh);

    TF1 *f = h->GetFunction("gaus");
    double mu    = f->GetParameter(1);
    double sigma = f->GetParameter(2);

    double chi2ndf = (f->GetNDF() > 0) ? f->GetChisquare()/f->GetNDF() : 0;

    double syst = max(highbound-mu, mu-lowbound);
    double totalErr = sqrt(syst*syst + fitError*fitError);

    float peakHeight = f->GetParameter(0);
    float histMax    = h->GetMaximum();
    float yMax       = std::max(peakHeight,histMax)*1.2;
    h->SetMaximum(yMax);

    TCanvas *canscale = new TCanvas("canscale","",600,420);
    h->Draw();    
    TLatex title;
    title.SetTextSize(0.1);
    title.SetTextAlign(22);
    title.SetTextFont(42);

    if(isQCD) title.DrawLatexNDC(0.5,0.95,"a_{fake}");
    else      title.DrawLatexNDC(0.5,0.95,"a_{V#gamma}");

    TLatex chi2Text;
    chi2Text.SetTextSize(0.06);
    chi2Text.SetNDC();
    chi2Text.DrawLatex(0.6,0.8,Form("#chi^{2}/ndf = %.2f",chi2ndf));

    if(RunYear==2016 && preVFP==1) CMS_lumi(canscale,1,ichannel,11);
    else if(RunYear==2016)         CMS_lumi(canscale,2,ichannel,11);
    else if(RunYear==2017)         CMS_lumi(canscale,3,ichannel,11);
    else if(RunYear==2018)         CMS_lumi(canscale,4,ichannel,11);

    string whichVFP="";
    if(RunYear==2016 && preVFP==1) whichVFP="preVFP";
    if(RunYear==2016 && preVFP==0) whichVFP="postVFP";

    string ch=(ichannel==1)?"eg":"mg";
    string bkg=(isQCD)?"QCD":"VGamma";
    canscale->Update();
    canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/%sScale_%s_%d%s_MET%.0f-%.0f_lep%.0f-%.0f.png", RunYear,whichVFP.c_str(), bkg.c_str(), ch.c_str(), RunYear,whichVFP.c_str(), lowmet,highmet, lowpt,highpt));

    delete canscale;
    delete h;
    return {mu, totalErr, chi2ndf};
}

/* ============================================================
   Generic dependence plotter (pT or MET)
   ============================================================ */
void makeDependencePlot(int RunYear, int preVFP, int ichannel,
    const vector<string>& files, const vector<double>& x,
    const vector<double>& xerr, const vector<double>& ptLow,
    const vector<double>& ptHigh,
    const vector<double>& metLow,
    const vector<double>& metHigh,
    const char* fullFile,
    const char* xTitle,
    double fullXlow, double fullXhigh,
    const char* outName,
    bool isQCD,
    const char* legendLabel, TH1F *hSF){

    TGraphErrors *g_bins = new TGraphErrors(files.size());

    for(size_t i=0;i<files.size();i++){
	auto res = extractGaussFromTxt(files[i].c_str(), isQCD, RunYear, preVFP, ichannel, ptLow[i],ptHigh[i], metLow[i],metHigh[i]);
        double sf    = res.mu;
    	double sferr = res.err;

    	g_bins->SetPoint(i, x[i], sf);
    	g_bins->SetPointError(i, xerr[i], sferr);

        cout << x[i]-xerr[i] << " " << x[i]+xerr[i]
             << "  SF = " << res.mu << " ± " << res.err
             << "  chi2/ndf = " << res.chi2ndf << endl;
    	
	int bin = hSF->FindBin(x[i]);
    	hSF->SetBinContent(bin, sf);
    	hSF->SetBinError(bin, sferr);
    }

    // ---- full range
    auto full = extractGaussFromTxt(fullFile, isQCD, RunYear, preVFP, ichannel, 0, 1000, 40, 70);
    double fullXcenter = 0.5*(fullXlow+fullXhigh);
    double fullXerr    = 0.5*(fullXhigh-fullXlow);

    TGraphErrors *g_full = new TGraphErrors(1);
    g_full->SetPoint(0, fullXcenter, full.mu);
    g_full->SetPointError(0, fullXerr, full.err);

    // ---- systematic band
    TGraphErrors *g_band = new TGraphErrors(1);
    g_band->SetPoint(0, fullXcenter, full.mu);
    g_band->SetPointError(0, fullXerr, full.err);
    g_band->SetFillStyle(3344);
    g_band->SetFillColor(kGray+1);

    // ---- draw
    TCanvas *c = new TCanvas(outName,"",600,600);
    TH1D *frame = new TH1D(Form("frame_%s", outName), Form(";%s;scalefactor",xTitle), 1, fullXlow, fullXhigh);
    frame->SetMinimum(0.0);
    frame->SetMaximum(2.5);
    frame->Draw();

    g_band->Draw("E2 same");

    g_bins->SetMarkerStyle(20);
    g_bins->SetMarkerColor(kBlue);
    g_bins->SetLineColor(kBlue);
    g_bins->Draw("EP same");

    g_full->SetMarkerStyle(20);
    g_full->SetMarkerColor(kBlack);
    g_full->Draw("P same");

    TLegend *leg = new TLegend(0.55,0.72,0.88,0.88);
    leg->AddEntry(g_bins, Form("%s (binned)",legendLabel),"lp");
    leg->AddEntry(g_full, Form("%s (full)",legendLabel),"p");
    leg->AddEntry(g_band,"Syst. unc.","f");
    leg->Draw();

    if(RunYear==2016 && preVFP)      CMS_lumi(c,1,ichannel,11);
    else if(RunYear==2016)           CMS_lumi(c,2,ichannel,11);
    else if(RunYear==2017)           CMS_lumi(c,3,ichannel,11);
    else if(RunYear==2018)           CMS_lumi(c,4,ichannel,11);

    c->SaveAs(outName);

}

/* ============================================================
   Main driver
   ============================================================ */
void plotVGammaScales(int RunYear, int preVFP, int ichannel)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    string base = "/eos/uscms/store/user/tmishra/VGamma/";

    string yearTag;
    if (RunYear==2016 && preVFP)       yearTag="2016preVFP";
    else if (RunYear==2016)            yearTag="2016postVFP";
    else                              yearTag=Form("%d",RunYear);

    const int nPtBins = 2;
    double ptBins[nPtBins+1] = {0, 70, 1000};
    const int nMetBins = 3;
    double metBins[nMetBins+1] = {40, 50, 60, 70};

    string chTag = (ichannel==1 ? "_eg_" : "_mg_");
    string fullFile = base + "VGamma_scalefactor" + chTag + yearTag + ".txt";

    TFile *fout = new TFile(Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGammaScaleHist%s%s.root",
                           yearTag.c_str(), chTag.c_str(), yearTag.c_str()), "RECREATE");

    TH1F *hVG_pt  = new TH1F("VGamma_pt","VGamma SF vs p_{T};p_{T} [GeV];SF",nPtBins, ptBins);
    TH1F *hVG_met = new TH1F("VGamma_met","VGamma SF vs MET;MET [GeV];SF",nMetBins, metBins);

    TH1F *hQCD_pt  = new TH1F("QCDfake_pt","Fake lepton SF vs p_{T};p_{T} [GeV];SF",nPtBins, ptBins);
    TH1F *hQCD_met = new TH1F("QCDfake_met","Fake lepton SF vs MET;MET [GeV];SF",nMetBins, metBins);

    hVG_pt->Sumw2();   hVG_met->Sumw2();
    hQCD_pt->Sumw2();  hQCD_met->Sumw2();

    vector<string> ptFiles = {
        base + "VGamma_scalefactor" + chTag + yearTag + "_MET40-70_lep0-70.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_MET40-70_lep70-1000.txt"
    };

	
	vector<double> ptX    = {35, 235};
	vector<double> ptXerr = {35, 165};
	vector<double> ptLow  = {0,70};
	vector<double> ptHigh = {70,1000};

    vector<double> metLow_pt  = {40,40,40,40};
    vector<double> metHigh_pt = {70,70,70,70};

    vector<string> metFiles = {
        base + "VGamma_scalefactor" + chTag + yearTag + "_MET40-50_lep0-1000.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_MET50-60_lep0-1000.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_MET60-70_lep0-1000.txt"
    };

    vector<double> metX    = {45, 55, 65};
    vector<double> metXerr = {5, 5, 5};

    vector<double> metLow  = {40,50,60};
    vector<double> metHigh = {50,60,70};

    vector<double> ptLow_met  = {0,0,0};
    vector<double> ptHigh_met = {1000,1000,1000};
    
    makeDependencePlot(RunYear,preVFP,ichannel,
        ptFiles,ptX,ptXerr,ptLow,ptHigh,metLow_pt,metHigh_pt, fullFile.c_str(),
        "Lepton p_{T} (GeV)",0,400,
        Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGamma_ptDependence%s%s.png", yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
        false,"V#gamma scale",hVG_pt);

//    makeDependencePlot(RunYear,preVFP,ichannel,
//        metFiles,metX,metXerr,
//        ptLow_met,ptHigh_met,metLow,metHigh,
//        fullFile.c_str(),
//        "p_{T}^{miss} (GeV)",40,70,
//        Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGamma_METDependence%s%s.png",
//             yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
//        false,"V#gamma scale",hVG_met);

    makeDependencePlot(RunYear,preVFP,ichannel,
        ptFiles,ptX,ptXerr,
        ptLow,ptHigh,metLow_pt,metHigh_pt,
        fullFile.c_str(),
        "Lepton p_{T} (GeV)",0,400,
        Form("/eos/uscms/store/user/tmishra/VGamma/%s/QCD_ptDependence%s%s.png",
             yearTag.c_str(),chTag.c_str(),yearTag.c_str()), true,"fake-lepton scale",hQCD_pt);

//    makeDependencePlot(RunYear,preVFP,ichannel,
//        metFiles,metX,metXerr,
//        ptLow_met,ptHigh_met,metLow,metHigh,
//        fullFile.c_str(),
//        "p_{T}^{miss} (GeV)",40,70,
//        Form("/eos/uscms/store/user/tmishra/VGamma/%s/QCD_METDependence%s%s.png",
//             yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
//        true,"fake-lepton scale",hQCD_met);

    fout->cd();
    hVG_pt->Write();
    hVG_met->Write();
    hQCD_pt->Write();
    hQCD_met->Write();
    fout->Close();
}
int main(int argc, char** argv)
{

    	plotVGammaScales(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
    	return 0;
}
