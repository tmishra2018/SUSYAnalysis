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
#include "TGraphAsymmErrors.h"

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

// mu: nominal fit (toy 0); errLo/errHi: 1 sigma, toys (16%/84% quantiles about their median) in quadrature with the
// MINOS error of the nominal fit (<file>_minos.txt from Fitfraction*.C, HESSE error if that file is missing)
struct GaussResult {
    double mu;
    double errLo;
    double errHi;
    double chi2ndf;
};

// name of the split variable in the per-bin png names ("HT", or "pho" for the photon-pT bins)
std::string g_splitName = "HT";

GaussResult extractGaussFromTxt(const char* filename, bool isQCD, int RunYear,
                              int preVFP,
                              int ichannel,
                              double lowpt,
                              double highpt,
                              double lowmet,
                              double highmet,
                              double lowht = -1,
                              double highht = -1)
{

ifstream file(filename);
float leplow, lephigh;
float fakescale, fakescaleerror;
float vgammascale, vgammascaleerror;
vector<double> values;
double fitError = 0, nominal = 0;
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
            nominal = val;
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

// toys only (values[0] is the nominal fit), for the 1 sigma toy spread
vector<double> toys(values.begin() + 1, values.end());
sort(toys.begin(), toys.end());
sort(values.begin(), values.end());

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
    Form("h_%.0f_%.0f_%.0f_%.0f_HT%.0f_%.0f",lowpt,highpt,lowmet,highmet,lowht,highht),
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
    // A Gaussian fit to a histogram of only a few toys (10 per bin) fails at random (gave e.g. -5.5 or 2.5
    // for toys around 1.2-1.4): use the toy mean then. syst below stays the max deviation of the toys.
    if(n < 50){
        double sum = 0;
        for(auto v : values) sum += v;
        mu = sum/n;
    }

    double chi2ndf = (f->GetNDF() > 0) ? f->GetChisquare()/f->GetNDF() : 0;

    // 1 sigma spread of the toys (template variations): 16% / 84% quantiles about the toy median
    int nt = toys.size();
    double toyMed = toys[nt/2], toyLo = toyMed - toys[int(0.16*nt)], toyHi = toys[int(0.84*nt)] - toyMed;
    // 1 sigma fit uncertainty of the nominal fit: MINOS (asymmetric), else the symmetric HESSE error
    double fitLo = fitError, fitHi = fitError;
    string fname(filename);
    ifstream minosfile(fname.substr(0, fname.size() - 4) + "_minos.txt");
    double mfake, mfakeLo, mfakeHi, mvg, mvgLo, mvgHi;
    if(minosfile >> mfake >> mfakeLo >> mfakeHi >> mvg >> mvgLo >> mvgHi){
        fitLo = isQCD ? mfakeLo : mvgLo;
        fitHi = isQCD ? mfakeHi : mvgHi;
    } else cout << "WARNING: no MINOS file for " << filename << ", using the HESSE fit error" << endl;
    double errLo = sqrt(toyLo*toyLo + fitLo*fitLo), errHi = sqrt(toyHi*toyHi + fitHi*fitHi);

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

    string htSuffix = "";
    if(lowht>=0 && highht>=0)
        htSuffix = Form("_%s%.0f-%.0f", g_splitName.c_str(), lowht, highht);

    canscale->Update();
    canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/%sScale_%s_%d%s_MET%.0f-%.0f_lep%.0f-%.0f%s.png", RunYear,whichVFP.c_str(), bkg.c_str(), ch.c_str(), RunYear,whichVFP.c_str(), lowmet,highmet, lowpt,highpt, htSuffix.c_str()));

    delete canscale;
    delete h;
    return {nominal, errLo, errHi, chi2ndf};
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
    const char* legendLabel, TH1F *hSF,
    const vector<double>& htLow = {},
    const vector<double>& htHigh = {},
    double fullHTlow = -1,
    double fullHThigh = -1,
    TH1F *hErrDown = nullptr, TH1F *hErrUp = nullptr){

    TGraphAsymmErrors *g_bins = new TGraphAsymmErrors(files.size());

    for(size_t i=0;i<files.size();i++){
	double lht = (i<htLow.size())  ? htLow[i]  : -1;
	double hht = (i<htHigh.size()) ? htHigh[i] : -1;
	auto res = extractGaussFromTxt(files[i].c_str(), isQCD, RunYear, preVFP, ichannel, ptLow[i],ptHigh[i], metLow[i],metHigh[i], lht, hht);
        double sf    = res.mu;

    	g_bins->SetPoint(i, x[i], sf);
    	g_bins->SetPointError(i, xerr[i], xerr[i], res.errLo, res.errHi);

        cout << x[i]-xerr[i] << " " << x[i]+xerr[i]
             << "  SF = " << res.mu << " -" << res.errLo << " +" << res.errHi
             << "  chi2/ndf = " << res.chi2ndf << endl;
    	
	// bin error = +1 sigma (the "up" variation in pred_*Bkg.C, Background/analysis_*Bkg.C); -/+ 1 sigma in hErrDown/hErrUp
	int bin = hSF->FindBin(x[i]);
    	hSF->SetBinContent(bin, sf);
    	hSF->SetBinError(bin, res.errHi);
	if(hErrDown) hErrDown->SetBinContent(bin, res.errLo);
	if(hErrUp)   hErrUp->SetBinContent(bin, res.errHi);
    }

    // ---- full range
    auto full = extractGaussFromTxt(fullFile, isQCD, RunYear, preVFP, ichannel, 0, 1000, 40, 70, fullHTlow, fullHThigh);
    double fullXcenter = 0.5*(fullXlow+fullXhigh);
    double fullXerr    = 0.5*(fullXhigh-fullXlow);

    TGraphAsymmErrors *g_full = new TGraphAsymmErrors(1);
    g_full->SetPoint(0, fullXcenter, full.mu);
    g_full->SetPointError(0, fullXerr, fullXerr, full.errLo, full.errHi);

    // ---- 1 sigma band
    TGraphAsymmErrors *g_band = new TGraphAsymmErrors(1);
    g_band->SetPoint(0, fullXcenter, full.mu);
    g_band->SetPointError(0, fullXerr, fullXerr, full.errLo, full.errHi);
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
    leg->AddEntry(g_band,"Unc. (1#sigma)","f");
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
    const int nHtBins = 3;
    double htBins[nHtBins+1] = {0, 100, 400, 3000};

    string chTag = (ichannel==1 ? "_eg_" : "_mg_");
    string fullFile = base + "VGamma_scalefactor" + chTag + yearTag + ".txt";

    TFile *fout = new TFile(Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGammaScaleHist%s%s.root",
                           yearTag.c_str(), chTag.c_str(), yearTag.c_str()), "RECREATE");

    TH1F *hVG_pt  = new TH1F("VGamma_pt","VGamma SF vs p_{T};p_{T} [GeV];SF",nPtBins, ptBins);
    TH1F *hVG_met = new TH1F("VGamma_met","VGamma SF vs MET;MET [GeV];SF",nMetBins, metBins);

    TH1F *hQCD_pt  = new TH1F("QCDfake_pt","Fake lepton SF vs p_{T};p_{T} [GeV];SF",nPtBins, ptBins);
    TH1F *hQCD_met = new TH1F("QCDfake_met","Fake lepton SF vs MET;MET [GeV];SF",nMetBins, metBins);

    TH1F *hVG_ht   = new TH1F("VGamma_ht","VGamma SF vs H_{T};H_{T} [GeV];SF",nHtBins, htBins);
    TH1F *hQCD_ht  = new TH1F("QCDfake_ht","Fake lepton SF vs H_{T};H_{T} [GeV];SF",nHtBins, htBins);

    // photon-pT binned SFs (default binning of the V+gamma and fake-lepton factors)
    double phoBins[4] = {0, 100, 200, 5000};
    TH1F *hVG_pho  = new TH1F("VGamma_pho","VGamma SF vs photon p_{T};p_{T}^{#gamma} [GeV];SF",3, phoBins);
    TH1F *hQCD_pho = new TH1F("QCDfake_pho","Fake lepton SF vs photon p_{T};p_{T}^{#gamma} [GeV];SF",3, phoBins);
    hVG_pho->Sumw2(); hQCD_pho->Sumw2();
    // -1 sigma / +1 sigma of the photon-pT binned factors (bin error of VGamma_pho / QCDfake_pho is the +1 sigma)
    TH1F *hVG_pho_dn  = new TH1F("VGamma_pho_errDown","",3, phoBins), *hVG_pho_up  = new TH1F("VGamma_pho_errUp","",3, phoBins);
    TH1F *hQCD_pho_dn = new TH1F("QCDfake_pho_errDown","",3, phoBins), *hQCD_pho_up = new TH1F("QCDfake_pho_errUp","",3, phoBins);

    hVG_pt->Sumw2();   hVG_met->Sumw2();  hVG_ht->Sumw2();
    hQCD_pt->Sumw2();  hQCD_met->Sumw2(); hQCD_ht->Sumw2();

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

    vector<string> htFiles = {
        base + "VGamma_scalefactor" + chTag + yearTag + "_HT0_100.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_HT100_400.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_HT400_3000.txt"
    };

    vector<double> htLow   = {0, 100, 400};
    vector<double> htHigh  = {100, 400, 3000};
    vector<double> htX     = {50, 250, 1700};
    vector<double> htXerr  = {50, 150, 1300};

    // HT bins are integrated over the full lepton pT and MET ranges
    vector<double> ptLow_ht  = {0,0,0};
    vector<double> ptHigh_ht = {1000,1000,1000};
    vector<double> metLow_ht  = {40,40,40};
    vector<double> metHigh_ht = {70,70,70};

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

    makeDependencePlot(RunYear,preVFP,ichannel,
        htFiles,htX,htXerr,
        ptLow_ht,ptHigh_ht,metLow_ht,metHigh_ht,
        fullFile.c_str(),
        "H_{T} (GeV)",0,3000,
        Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGamma_HTDependence%s%s.png",
             yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
        false,"V#gamma scale",hVG_ht, htLow, htHigh);

    makeDependencePlot(RunYear,preVFP,ichannel,
        htFiles,htX,htXerr,
        ptLow_ht,ptHigh_ht,metLow_ht,metHigh_ht,
        fullFile.c_str(),
        "H_{T} (GeV)",0,3000,
        Form("/eos/uscms/store/user/tmishra/VGamma/%s/QCD_HTDependence%s%s.png",
             yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
        true,"fake-lepton scale",hQCD_ht, htLow, htHigh);

    // ---- photon-pT split (only if all its fit outputs exist)
    auto allExist = [](const vector<string>& fs){ for(auto &f : fs){ ifstream t(f); if(!t.good()) return false; } return true; };
    vector<string> phoFiles = {
        base + "VGamma_scalefactor" + chTag + yearTag + "_pho0_100.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_pho100_200.txt",
        base + "VGamma_scalefactor" + chTag + yearTag + "_pho200_5000.txt"
    };
    vector<double> phoX = {50, 150, 600}, phoXerr = {50, 50, 400};
    vector<double> phoLow = {0, 100, 200}, phoHigh = {100, 200, 5000};
    vector<double> lep3Low = {0,0,0}, lep3High = {1000,1000,1000}, met3Low = {40,40,40}, met3High = {70,70,70};
    if(allExist(phoFiles)){
        g_splitName = "pho";   // per-bin png names _pho<lo>-<hi> (photon window passed through the HT slots, used for naming only)
        makeDependencePlot(RunYear,preVFP,ichannel, phoFiles,phoX,phoXerr, lep3Low,lep3High,met3Low,met3High, fullFile.c_str(),
            "Photon p_{T} (GeV)",0,1000,
            Form("/eos/uscms/store/user/tmishra/VGamma/%s/VGamma_phoPtDependence%s%s.png", yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
            false,"V#gamma scale",hVG_pho, phoLow, phoHigh, -1, -1, hVG_pho_dn, hVG_pho_up);
        makeDependencePlot(RunYear,preVFP,ichannel, phoFiles,phoX,phoXerr, lep3Low,lep3High,met3Low,met3High, fullFile.c_str(),
            "Photon p_{T} (GeV)",0,1000,
            Form("/eos/uscms/store/user/tmishra/VGamma/%s/QCD_phoPtDependence%s%s.png", yearTag.c_str(),chTag.c_str(),yearTag.c_str()),
            true,"fake-lepton scale",hQCD_pho, phoLow, phoHigh, -1, -1, hQCD_pho_dn, hQCD_pho_up);
        g_splitName = "HT";
    } else cout << "photon-pT split fits not found, VGamma_pho/QCDfake_pho left empty" << endl;

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
    hVG_ht->Write();
    hQCD_pt->Write();
    hQCD_met->Write();
    hQCD_ht->Write();
    hVG_pho->Write();  hQCD_pho->Write();
    hVG_pho_dn->Write(); hVG_pho_up->Write(); hQCD_pho_dn->Write(); hQCD_pho_up->Write();
    fout->Close();
}
int main(int argc, char** argv)
{

    	plotVGammaScales(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
    	return 0;
}
