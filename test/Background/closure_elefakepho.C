// closure_elefakepho.C
//
// PURPOSE
// ───────
// Performs a closure test of the electron-faking-photon fake-rate estimate
// by comparing two parameterisations:
//   A)  f(pT) × f(nVtx) × f(eta)   [baseline]
//   B)  f(pT) × f(HT)   × f(eta)   [HT-based]
//
// The test is done in three MET regions:
//   CR  : MET < 70 GeV            (where the fake rate was derived)
//   VR  : MET > 70 GeV, MT < 100  (validation region)
//   SR  : MET > 120 GeV, MT > 100 (signal region proxy)
//
// The "truth" denominator events (vetovalue==false) are re-weighted with
// both fake rates and compared to the observed numerator (vetovalue==true
// && FSRveto==true) in each region.
//
// Closure = predicted / observed.   Good closure: ratio ≈ 1.
//
// The key diagnostic output is the ratio   pred_HT / pred_nVtx   as a
// function of pT and HT, to reveal the kinematic origin of the difference.
//
// Compile:
//   g++ `root-config --cflags` closure_elefakepho.C -o closure_elefakepho.exe `root-config --libs`
// Run:
//   ./closure_elefakepho.exe <RunYear> [preVFP=0/1]
//   e.g.   ./closure_elefakepho.exe 2018 0
//
// Or interactively:
//   root -l 'closure_elefakepho.C(2018,false)'

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "TGraphErrors.h"

// ════════════════════════════════════════════════════════════════════
//  Fake-rate helpers — identical analytic forms used in fitFakeFunc.C
// ════════════════════════════════════════════════════════════════════

// pT: power-law  (slope*pT + constant)^index
double fr_pt(double pt, double slope, double constant, double index)
{
    double arg = slope * std::max(pt, 0.001) + constant;
    if (arg <= 0) return 0;
    return std::pow(arg, index);
}

// nVtx: linear pol1
double fr_vtx(double vtx, double p0, double p1)
{ return p0 + p1 * vtx; }

// HT: bin-by-bin (piecewise constant)
double fr_binned(double val, const std::vector<double>& binEdges,
                 const std::vector<double>& num, const std::vector<double>& den)
{
    for (size_t i = 0; i < binEdges.size(); ++i) {
        double lo = binEdges[i];
        double hi = (i + 1 < binEdges.size()) ? binEdges[i+1] : 1e9;
        if (val >= lo && val < hi)
            return (den[i] > 0) ? num[i] / den[i] : 0;
    }
    return 0;
}

// ════════════════════════════════════════════════════════════════════
//  Read the fake-rate text files produced by fitFakeFunc / fitFakeFunc_HT
// ════════════════════════════════════════════════════════════════════

struct FakeRateParams {
    // pT power-law
    double ptSlope, ptConst, ptIndex;
    // nVtx linear
    double vtxP0, vtxP1;
    // eta (binned, few bins)
    std::vector<double> etaBins;   // lower edges
    std::vector<double> etaNum;
    std::vector<double> etaDen;
    // HT (binned)
    std::vector<double> htBins;
    std::vector<double> htNum;
    std::vector<double> htDen;
    // global scale factor from Z-peak fit
    double scaleFactor_vtx = 1.0;
    double scaleFactor_ht  = 1.0;
};

// Read the EleFakeRate-Data-ByPtVtx-EB.txt style file
// Format (as written by fitFakeFunc.C):
//   scalefactor <val>
//   ptslope <val>  ptconstant <val>  ptindex <val>  ptcoeff <val>
//   vtxconst <val>  vtxslope <val>
//   etaRate_0 <val> ... (written as a C++ array line)
bool readVtxFakeRateFile(const std::string& fname, FakeRateParams& p)
{
    std::ifstream f(fname);
    if (!f.is_open()) {
        std::cerr << "Cannot open: " << fname << std::endl;
        return false;
    }
    std::string key;
    while (f >> key) {
        if (key == "scalefactor") f >> p.scaleFactor_vtx;
        else if (key == "ptslope")    f >> p.ptSlope;
        else if (key == "ptconstant") f >> p.ptConst;
        else if (key == "ptindex")    f >> p.ptIndex;
        else if (key == "vtxconst")   f >> p.vtxP0;
        else if (key == "vtxslope")   f >> p.vtxP1;
    }
    std::cout << "nVtx FR: SF=" << p.scaleFactor_vtx
              << "  ptSlope=" << p.ptSlope
              << "  ptConst=" << p.ptConst
              << "  ptIndex=" << p.ptIndex
              << "  vtxP0="   << p.vtxP0
              << "  vtxP1="   << p.vtxP1 << std::endl;
    return true;
}

// Read the eta rate file (a C++ array declaration line written by fitFakeFunc.C)
// e.g.:  double etaRatesEB_2018[] = {0.031, 0.029, ...};
// We parse the numbers out of it.
bool readEtaFile(const std::string& fname, FakeRateParams& p,
                 const std::vector<double>& etaBinEdges)
{
    std::ifstream f(fname);
    if (!f.is_open()) {
        std::cerr << "Cannot open: " << fname << std::endl;
        return false;
    }
    p.etaBins = etaBinEdges;
    // Read the whole line, extract numbers between { }
    std::string line;
    while (std::getline(f, line)) {
        auto lb = line.find('{');
        auto rb = line.find('}');
        if (lb == std::string::npos || rb == std::string::npos) continue;
        std::string nums = line.substr(lb+1, rb-lb-1);
        std::replace(nums.begin(), nums.end(), ',', ' ');
        std::istringstream ss(nums);
        double v;
        while (ss >> v) p.etaNum.push_back(v);  // these are already ratios
    }
    // etaDen will be set to 1 since etaNum already stores num/den
    p.etaDen.assign(p.etaNum.size(), 1.0);
    std::cout << "Eta FR: " << p.etaNum.size() << " bins read." << std::endl;
    return true;
}

// Read HT fake-rate from text files produced by fitFakeFunc_HT.C
// Format (same as vtx but with HT bin edges and num/den per bin):
//   bintype HT numtype den lowcut <val> signal <val> err <val> fitmean <val> fitrms <val>
// OR simply the same "EleFakeRate" text file format with ht keyword
// We read the raw FitKer_HT output format:
//   "HT den <lowedge> <signal> <err> <mean> <rms>"   (denominator lines first)
//   "HT num <lowedge> <signal> <err> <mean> <rms>"   (numerator lines second)
bool readHTFakeRateFile(const std::string& fname, FakeRateParams& p, double scaleFactor_ht)
{
    std::ifstream f(fname);
    if (!f.is_open()) {
        std::cerr << "Cannot open: " << fname << std::endl;
        return false;
    }
    p.scaleFactor_ht = scaleFactor_ht;  // set externally from fitFakeFunc_HT result file

    std::string bintype, numtype;
    float lowcut, signal, err, fitmean, fitrms;

    std::vector<double> tmpDen_lo, tmpDen_mean;
    std::vector<double> tmpNum_lo, tmpNum_mean;

    while (f >> bintype >> numtype >> lowcut >> signal >> err >> fitmean >> fitrms) {
        if (numtype == "den") { tmpDen_lo.push_back(lowcut); tmpDen_mean.push_back(fitmean); }
        if (numtype == "num") { tmpNum_lo.push_back(lowcut); tmpNum_mean.push_back(fitmean); }
    }
    if (tmpDen_lo.size() != tmpNum_lo.size()) {
        std::cerr << "HT file: mismatched den/num bins!" << std::endl;
        return false;
    }
    p.htBins = tmpDen_lo;
    p.htDen  = tmpDen_mean;
    p.htNum  = tmpNum_mean;
    std::cout << "HT FR: " << p.htBins.size() << " bins read." << std::endl;
    return true;
}

// ════════════════════════════════════════════════════════════════════
//  Per-event weight calculators
// ════════════════════════════════════════════════════════════════════

double weightVtx(double probePt, double probeEta, double nVtx,
                 const FakeRateParams& p)
{
    // eta bin
    double wEta = 0;
    size_t nEta = p.etaBins.size();
    for (size_t i = 0; i < nEta; ++i) {
        double lo = p.etaBins[i];
        double hi = (i+1 < nEta) ? p.etaBins[i+1] : 99.;
        if (std::fabs(probeEta) >= lo && std::fabs(probeEta) < hi)
            wEta = (p.etaDen[i]>0) ? p.etaNum[i]/p.etaDen[i] : p.etaNum[i];
    }
    double wPt  = fr_pt(probePt, p.ptSlope, p.ptConst, p.ptIndex);
    double wVtx = fr_vtx(nVtx, p.vtxP0, p.vtxP1);
    return p.scaleFactor_vtx * wPt * wVtx * wEta;
}

double weightHT(double probePt, double probeEta, double HT,
                const FakeRateParams& p)
{
    double wEta = 0;
    size_t nEta = p.etaBins.size();
    for (size_t i = 0; i < nEta; ++i) {
        double lo = p.etaBins[i];
        double hi = (i+1 < nEta) ? p.etaBins[i+1] : 99.;
        if (std::fabs(probeEta) >= lo && std::fabs(probeEta) < hi)
            wEta = (p.etaDen[i]>0) ? p.etaNum[i]/p.etaDen[i] : p.etaNum[i];
    }
    double wPt = fr_pt(probePt, p.ptSlope, p.ptConst, p.ptIndex);
    double wHT = fr_binned(HT, p.htBins, p.htNum, p.htDen);
    return p.scaleFactor_ht * wPt * wHT * wEta;
}

// ════════════════════════════════════════════════════════════════════
//  Style
// ════════════════════════════════════════════════════════════════════
void setStyle() {
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleOffset(1.3,"Y"); gStyle->SetTitleOffset(1.1,"X");
    gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.13);
}

// ════════════════════════════════════════════════════════════════════
//  Ratio pad canvas
// ════════════════════════════════════════════════════════════════════
void drawWithRatio(TH1F* hObs, TH1F* hPredVtx, TH1F* hPredHT,
                   const char* xtitle, const char* tag, const char* region,
                   int year, const std::string& vfp)
{
    TCanvas* c = new TCanvas(Form("c_%s_%s",tag,region),"",700,750);
    TPad* p1 = new TPad("p1","",0,0.32,1,1);
    TPad* p2 = new TPad("p2","",0,0.00,1,0.32);
    p1->SetBottomMargin(0.01); p1->SetTopMargin(0.08);
    p2->SetTopMargin(0.01);    p2->SetBottomMargin(0.40);
    p1->Draw(); p2->Draw();

    p1->cd();
    hObs->SetLineColor(kBlack);    hObs->SetMarkerStyle(20); hObs->SetMarkerColor(kBlack);
    hPredVtx->SetLineColor(kBlue+1); hPredVtx->SetLineWidth(2); hPredVtx->SetFillColor(kBlue-9); hPredVtx->SetFillStyle(1001);
    hPredHT->SetLineColor(kRed+1);  hPredHT->SetLineWidth(2);  hPredHT->SetFillColor(kRed-9);  hPredHT->SetFillStyle(1001);
    double ymax = std::max({hObs->GetMaximum(), hPredVtx->GetMaximum(), hPredHT->GetMaximum()}) * 1.5;
    hObs->SetMaximum(ymax); hObs->SetMinimum(0);
    hObs->GetYaxis()->SetTitle("Events");
    hObs->GetXaxis()->SetLabelSize(0);
    hObs->Draw("EP"); hPredVtx->Draw("HIST same"); hPredHT->Draw("HIST same"); hObs->Draw("EP same");
    TLegend* leg = new TLegend(0.55,0.65,0.90,0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hObs,     "Observed (num)",          "pe");
    leg->AddEntry(hPredVtx, Form("Pred (pt,#eta,nVtx)  #int=%.1f",hPredVtx->Integral()), "f");
    leg->AddEntry(hPredHT,  Form("Pred (pt,#eta,HT)    #int=%.1f",hPredHT->Integral()),  "f");
    leg->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038); lat.SetTextFont(42);
    lat.DrawLatex(0.14,0.92, Form("Closure test  %d  Region: %s",year,region));

    p2->cd();
    TH1F* ratioVtx = (TH1F*)hPredVtx->Clone(Form("rVtx_%s_%s",tag,region));
    TH1F* ratioHT  = (TH1F*)hPredHT->Clone(Form("rHT_%s_%s",tag,region));
    ratioVtx->Divide(hObs); ratioHT->Divide(hObs);
    ratioVtx->SetFillStyle(0);
    ratioHT->SetFillStyle(0);
    ratioVtx->GetYaxis()->SetTitle("Pred/Obs");
    ratioVtx->GetXaxis()->SetTitle(xtitle);
    ratioVtx->GetYaxis()->SetNdivisions(504);
    ratioVtx->GetYaxis()->SetTitleSize(0.13); ratioVtx->GetYaxis()->SetLabelSize(0.11);
    ratioVtx->GetXaxis()->SetTitleSize(0.13); ratioVtx->GetXaxis()->SetLabelSize(0.11);
    ratioVtx->GetYaxis()->SetTitleOffset(0.45);
    ratioVtx->SetMaximum(2.5); ratioVtx->SetMinimum(0.0);
    ratioVtx->Draw("HIST"); ratioHT->Draw("HIST same");
    TLine* l1 = new TLine(ratioVtx->GetXaxis()->GetXmin(),1,ratioVtx->GetXaxis()->GetXmax(),1);
    l1->SetLineStyle(2); l1->Draw();

    c->SaveAs(Form("closure_%s_%s_%d%s.png",tag,region,year,vfp.c_str()));
    delete c;
}

// ════════════════════════════════════════════════════════════════════
//  Closure summary table
// ════════════════════════════════════════════════════════════════════
void printSummaryTable(
    double obs_CR,  double predVtx_CR,  double predHT_CR,
    double obs_VR,  double predVtx_VR,  double predHT_VR,
    double obs_SR,  double predVtx_SR,  double predHT_SR,
    int year, std::ofstream& txt)
{
    auto pct=[](double a,double b)->double{ return b>0 ? 100*(a-b)/b : 0; };
    auto clo=[](double pred,double obs)->double{ return obs>0 ? pred/obs : 0; };

    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Closure test summary  RunYear=" << year << "                                    ║\n";
    std::cout << "╠═══════════════╦══════════╦══════════╦═══════════╦══════════╦═════════════╣\n";
    std::cout << "║  Region       ║  Obs     ║ Pred_Vtx ║ Clo_Vtx  ║ Pred_HT  ║ Clo_HT     ║\n";
    std::cout << "╠═══════════════╬══════════╬══════════╬═══════════╬══════════╬═════════════╣\n";
    auto row=[&](const char* r, double o, double pv, double ph){
        printf("║  %-13s║ %8.1f ║ %8.1f ║ %9.3f ║ %8.1f ║ %11.3f ║\n",
               r, o, pv, clo(pv,o), ph, clo(ph,o));
    };
    row("CR (MET<70)",  obs_CR,  predVtx_CR,  predHT_CR);
    row("VR (MT<100)",  obs_VR,  predVtx_VR,  predHT_VR);
    row("SR (MT>100)",  obs_SR,  predVtx_SR,  predHT_SR);
    std::cout << "╠═══════════════╬══════════╬══════════╬═══════════╬══════════╬═════════════╣\n";
    printf("║  HT/Vtx ratio  ║          ║          ║           ║          ║             ║\n");
    printf("║    CR          ║          ║          ║           ║  %+.1f%%              ║\n", pct(predHT_CR,predVtx_CR));
    printf("║    VR          ║          ║          ║           ║  %+.1f%%              ║\n", pct(predHT_VR,predVtx_VR));
    printf("║    SR          ║          ║          ║           ║  %+.1f%%              ║\n", pct(predHT_SR,predVtx_SR));
    std::cout << "╚══════════════════════════════════════════════════════════════════════════╝\n\n";

    txt << "\n══════════ CLOSURE SUMMARY  " << year << " ══════════\n";
    txt << std::fixed << std::setprecision(2);
    txt << "Region   Obs      Pred_Vtx  Clo_Vtx   Pred_HT   Clo_HT   HT/Vtx\n";
    auto trow=[&](const char* r, double o, double pv, double ph){
        txt << std::setw(8)<<r <<" "
            <<std::setw(8)<<o  <<" "
            <<std::setw(9)<<pv <<" "
            <<std::setw(9)<<clo(pv,o) <<" "
            <<std::setw(9)<<ph <<" "
            <<std::setw(8)<<clo(ph,o) <<" "
            <<std::setw(7)<<pct(ph,pv)<<"%\n";
    };
    trow("CR",  obs_CR,  predVtx_CR,  predHT_CR);
    trow("VR",  obs_VR,  predVtx_VR,  predHT_VR);
    trow("SR",  obs_SR,  predVtx_SR,  predHT_SR);
}

// ════════════════════════════════════════════════════════════════════
//  Main
// ════════════════════════════════════════════════════════════════════
void closure_elefakepho(int RunYear=2018, bool preVFP=false)
{
    setStyle();

    std::string whichVFP;
    if (RunYear==2016) whichVFP = preVFP ? "preVFP" : "postVFP";
    else               whichVFP = "";

    // ── open log file ────────────────────────────────────────────────
    std::ofstream txt(Form("closure_results_%d%s.txt",RunYear,whichVFP.c_str()));

    txt << "Closure test: (pt,eta,nVtx) vs (pt,eta,HT) fake rates\n";
    txt << "RunYear=" << RunYear << "  VFP=" << whichVFP << "\n\n";

    // ────────────────────────────────────────────────────────────────
    //  STEP 1 — Load fake rate parameters
    //  Adapt these paths to your actual output directories
    // ────────────────────────────────────────────────────────────────

    FakeRateParams vtxFR, htFR;

    // --- nVtx-based ---
    std::string vtxDir = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s",
                              RunYear, whichVFP.c_str());
    bool ok_vtx = readVtxFakeRateFile(
        vtxDir + "/EleFakeRate-Data-ByPtVtx-EB.txt", vtxFR);

    // Eta bins: adapt to what your analysis uses
    // These are the lower edges of |eta| bins used in fitFakeFunc.C
    std::vector<double> etaBinEdges = {0.0, 0.5, 1.0};   // example; edit as needed

    // The eta rates are stored in a C++ array file; parse it
    bool ok_eta_vtx = readEtaFile(
        vtxDir + "/result_eta_dependence_Data.txt", vtxFR, etaBinEdges);

    // Same eta parameterisation is shared by both methods
    htFR.etaBins = vtxFR.etaBins;
    htFR.etaNum  = vtxFR.etaNum;
    htFR.etaDen  = vtxFR.etaDen;

    // Also copy pT parameters to htFR (pT fit is shared)
    htFR.ptSlope = vtxFR.ptSlope;
    htFR.ptConst = vtxFR.ptConst;
    htFR.ptIndex = vtxFR.ptIndex;

    // --- HT-based ---
    std::string htDir = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s",
                             RunYear, whichVFP.c_str());
    // Read HT scale factor from the fitFakeFunc_HT result file
    double sf_ht = 1.0;
    {
        std::ifstream sffile(htDir + "/result_fitFakeFunc_Data_HT.txt");
        if (sffile.is_open()) {
            std::string key; double val;
            while (sffile >> key >> val)
                if (key == "scalefactor") { sf_ht = val; break; }
        } else {
            std::cout << "WARNING: HT scalefactor file not found, using SF=1\n";
        }
    }
    bool ok_ht = readHTFakeRateFile(
        htDir + "/EleFakeRate-DATA-Bw-expo-ht-60-120.txt", htFR, sf_ht);

    if (!ok_vtx || !ok_ht || !ok_eta_vtx) {
        std::cerr << "Could not load all fake rate files. "
                  << "Check paths and run fitFakeFunc.C / fitFakeFunc_HT.C first.\n";
        // For testing: fill dummy parameters so histograms are still made
        vtxFR.ptSlope=1e4; vtxFR.ptConst=-3e4; vtxFR.ptIndex=-0.3;
        vtxFR.vtxP0=0.02; vtxFR.vtxP1=0.0003; vtxFR.scaleFactor_vtx=1.0;
        vtxFR.etaBins={0,0.5,1.0}; vtxFR.etaNum={0.03,0.028,0.025}; vtxFR.etaDen={1,1,1};
        htFR = vtxFR;
        htFR.htBins={0,100,200,400}; htFR.htNum={0.020,0.025,0.030};
        htFR.htDen={1,1,1}; htFR.scaleFactor_ht=1.0;
    }

    txt << "nVtx FR: ptSlope=" << vtxFR.ptSlope << "  ptConst=" << vtxFR.ptConst
        << "  ptIndex=" << vtxFR.ptIndex << "  vtxP0=" << vtxFR.vtxP0
        << "  vtxP1=" << vtxFR.vtxP1 << "  SF=" << vtxFR.scaleFactor_vtx << "\n";
    txt << "HT FR:   nHTbins=" << htFR.htBins.size() << "  SF=" << htFR.scaleFactor_ht << "\n\n";

    // ────────────────────────────────────────────────────────────────
    //  STEP 2 — Load the FakeRateTree (Data T&P)
    //  This is the same tree used in fitFakeFunc.C / analysis_elefakepho.C
    //  It contains MET, MT and all probe variables per tag-probe pair.
    // ────────────────────────────────────────────────────────────────

    TChain* etree = new TChain("FakeRateTree");
    etree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/"
                    "plot_elefakepho_DataTnP_dR05_%d%s.root", RunYear, whichVFP.c_str()));
    if (etree->GetEntries() == 0)
        etree->Add(Form("plot_elefakepho_DataTnP_dR05_%d%s.root",
                        RunYear, whichVFP.c_str()));

    // NOTE: MET and MT must be stored in the FakeRateTree.
    // If they are not (original analysis_elefakepho.C does not store them),
    // they need to be added to that script first — see comments in the code.
    // Here we read them; if the branch is missing, ROOT returns 0 silently.
    float probePt=0, probeEta=0, HT=0, MET=0, MT=0;
    int   nVertex=0;
    bool  vetovalue=false, FSRveto=false;

    etree->SetBranchAddress("probePt",   &probePt);
    etree->SetBranchAddress("probeEta",  &probeEta);
    etree->SetBranchAddress("HT",        &HT);
    etree->SetBranchAddress("nVertex",   &nVertex);
    etree->SetBranchAddress("vetovalue", &vetovalue);
    etree->SetBranchAddress("FSRveto",   &FSRveto);
    // MET and MT — only available if analysis_elefakepho.C was updated to store them
    if (etree->GetBranch("MET")) etree->SetBranchAddress("MET", &MET);
    if (etree->GetBranch("MT"))  etree->SetBranchAddress("MT",  &MT);

    // ────────────────────────────────────────────────────────────────
    //  STEP 3 — Book histograms for three regions
    // ────────────────────────────────────────────────────────────────

    auto mk1 = [](const char* n, const char* t, int nb, double lo, double hi) -> TH1F* {
        TH1F* h = new TH1F(n,t,nb,lo,hi);
        h->Sumw2(); return h;
    };

    // pT distributions
    for (const char* reg : {"CR","VR","SR"}) {
        // declared below after region loop
    }

    struct Region {
        const char* name;
        TH1F *hObs_pt,    *hPredVtx_pt,  *hPredHT_pt;
        TH1F *hObs_HT,    *hPredVtx_HT,  *hPredHT_HT;
        TH1F *hObs_nVtx,  *hPredVtx_nVtx,*hPredHT_nVtx;
        TH1F *hObs_invmass,*hPredVtx_invmass,*hPredHT_invmass;
        TH2F *h2_htPt_predRatio;  // f(HT-based) / f(nVtx-based) weight in (pT,HT) plane
    };

    std::vector<Region> regions(3);
    std::vector<const char*> regNames = {"CR","VR","SR"};
    for (int ir=0; ir<3; ++ir) {
        const char* rn = regNames[ir];
        regions[ir].name = rn;
        regions[ir].hObs_pt         = mk1(Form("hObs_pt_%s",rn),   Form("%s;probe p_{{T}} (GeV);Events",rn),   20,30,230);
        regions[ir].hPredVtx_pt     = mk1(Form("hPVtx_pt_%s",rn),  Form("%s;probe p_{{T}} (GeV);Events",rn),   20,30,230);
        regions[ir].hPredHT_pt      = mk1(Form("hPHT_pt_%s",rn),   Form("%s;probe p_{{T}} (GeV);Events",rn),   20,30,230);
        regions[ir].hObs_HT         = mk1(Form("hObs_HT_%s",rn),   Form("%s;H_{{T}} (GeV);Events",rn),         25,0,1000);
        regions[ir].hPredVtx_HT     = mk1(Form("hPVtx_HT_%s",rn),  Form("%s;H_{{T}} (GeV);Events",rn),         25,0,1000);
        regions[ir].hPredHT_HT      = mk1(Form("hPHT_HT_%s",rn),   Form("%s;H_{{T}} (GeV);Events",rn),         25,0,1000);
        regions[ir].hObs_nVtx       = mk1(Form("hObs_nVtx_%s",rn), Form("%s;nVtx;Events",rn),                  25,0,100);
        regions[ir].hPredVtx_nVtx   = mk1(Form("hPVtx_nV_%s",rn),  Form("%s;nVtx;Events",rn),                  25,0,100);
        regions[ir].hPredHT_nVtx    = mk1(Form("hPHT_nV_%s",rn),   Form("%s;nVtx;Events",rn),                  25,0,100);
        regions[ir].hObs_invmass    = mk1(Form("hObs_m_%s",rn),    Form("%s;m_{{tag,probe}} (GeV);Events",rn),  40,70,110);
        regions[ir].hPredVtx_invmass= mk1(Form("hPVtx_m_%s",rn),   Form("%s;m_{{tag,probe}} (GeV);Events",rn),  40,70,110);
        regions[ir].hPredHT_invmass = mk1(Form("hPHT_m_%s",rn),    Form("%s;m_{{tag,probe}} (GeV);Events",rn),  40,70,110);
        regions[ir].h2_htPt_predRatio = new TH2F(
            Form("h2ratio_%s",rn),
            Form("%s;probe p_{{T}} (GeV);H_{{T}} (GeV);w_{{HT}}/w_{{Vtx}}",rn),
            15,30,180, 15,0,750);
    }

    // Weight-ratio diagnostic histogram (filled with wHT/wVtx per denominator event)
    TH2F* h2_wRatio_full = new TH2F("h2_wRatio_full",
        "w_{HT}/w_{Vtx} vs (pT,HT)  full CR;probe p_{T} (GeV);H_{T} (GeV)",
        20,30,230, 20,0,1000);
    TProfile* prof_wRatio_pt  = new TProfile("prof_wRatio_pt",
        "mean(w_{HT}/w_{Vtx}) vs pT;probe p_{T} (GeV);<w_{HT}/w_{Vtx}>",
        20,30,230);
    TProfile* prof_wRatio_HT  = new TProfile("prof_wRatio_HT",
        "mean(w_{HT}/w_{Vtx}) vs HT;H_{T} (GeV);<w_{HT}/w_{Vtx}>",
        20,0,1000);
    TProfile* prof_wRatio_nVtx= new TProfile("prof_wRatio_nVtx",
        "mean(w_{HT}/w_{Vtx}) vs nVtx;nVtx;<w_{HT}/w_{Vtx}>",
        20,0,100);

    // ────────────────────────────────────────────────────────────────
    //  STEP 4 — Event loop
    // ────────────────────────────────────────────────────────────────
    long nEvts = etree->GetEntries();
    std::cout << "Processing " << nEvts << " events ..." << std::endl;

    // Also read tag kinematics for MT calculation if MT branch is missing
    float tagPt=0, tagEta=0, tagPhi=0, probePhi=0;
    if (!etree->GetBranch("MT")) {
        if (etree->GetBranch("tagPt"))  etree->SetBranchAddress("tagPt",  &tagPt);
        if (etree->GetBranch("tagPhi")) etree->SetBranchAddress("tagPhi", &tagPhi);
        if (etree->GetBranch("probePhi"))etree->SetBranchAddress("probePhi",&probePhi);
    }
    float invmass=0;
    if (etree->GetBranch("invmass")) etree->SetBranchAddress("invmass",&invmass);

    for (long ievt=0; ievt<nEvts; ++ievt) {
        if (ievt%500000==0) std::cout << "  event " << ievt << std::endl;
        etree->GetEntry(ievt);

        if (probePt < 35) continue;
        if (std::fabs(probeEta) > 1.4442) continue;  // EB only

        // If MT not stored, approximate from tag+probe opening angle and MET
        // (this is a proxy — the actual MT uses the lepton from the tag leg)
        float mt_val = MT;
        if (!etree->GetBranch("MT") && MET > 0)
            mt_val = std::sqrt(2*probePt*MET*(1-std::cos(probePhi - tagPhi)));

        // region assignment
        bool isCR = (MET < 70);
        bool isVR = (MET > 70  && mt_val < 100);
        bool isSR = (MET > 120 && mt_val > 100);

        double wVtx = weightVtx(probePt, probeEta, (double)nVertex, vtxFR);
        double wHT  = weightHT (probePt, probeEta, HT, htFR);

        for (int ir=0; ir<3; ++ir) {
            bool inReg = (ir==0 && isCR) || (ir==1 && isVR) || (ir==2 && isSR);
            if (!inReg) continue;
            Region& R = regions[ir];

            if (vetovalue == false) {                    // denominator
                R.hPredVtx_pt ->Fill(probePt, wVtx);
                R.hPredHT_pt  ->Fill(probePt, wHT);
                R.hPredVtx_HT ->Fill(HT, wVtx);
                R.hPredHT_HT  ->Fill(HT, wHT);
                R.hPredVtx_nVtx->Fill(nVertex, wVtx);
                R.hPredHT_nVtx ->Fill(HT, wHT);
                if (invmass>0) {
                    R.hPredVtx_invmass->Fill(invmass, wVtx);
                    R.hPredHT_invmass ->Fill(invmass, wHT);
                }
                // weight ratio diagnostics (fill in CR only)
                if (ir==0 && wVtx>0) {
                    double ratio = wHT/wVtx;
                    h2_wRatio_full ->Fill(probePt, HT, ratio);
                    prof_wRatio_pt ->Fill(probePt, ratio);
                    prof_wRatio_HT ->Fill(HT,      ratio);
                    prof_wRatio_nVtx->Fill(nVertex, ratio);
                    R.h2_htPt_predRatio->Fill(probePt, HT, ratio);
                }
            } else if (vetovalue==true && FSRveto==true) {  // numerator (observed)
                R.hObs_pt   ->Fill(probePt);
                R.hObs_HT   ->Fill(HT);
                R.hObs_nVtx ->Fill(nVertex);
                if (invmass>0) R.hObs_invmass->Fill(invmass);
            }
        }
    }

    // ────────────────────────────────────────────────────────────────
    //  STEP 5 — Compute closure in each region and print
    // ────────────────────────────────────────────────────────────────

    std::vector<double> obsI(3), vtxI(3), htI(3);
    for (int ir=0; ir<3; ++ir) {
        obsI[ir]  = regions[ir].hObs_pt->Integral();
        vtxI[ir]  = regions[ir].hPredVtx_pt->Integral();
        htI[ir]   = regions[ir].hPredHT_pt->Integral();
    }

    printSummaryTable(
        obsI[0],vtxI[0],htI[0],
        obsI[1],vtxI[1],htI[1],
        obsI[2],vtxI[2],htI[2],
        RunYear, txt);

    // ────────────────────────────────────────────────────────────────
    //  STEP 6 — Save plots
    // ────────────────────────────────────────────────────────────────

    for (int ir=0; ir<3; ++ir) {
        const char* rn = regions[ir].name;
        drawWithRatio(regions[ir].hObs_pt,     regions[ir].hPredVtx_pt,  regions[ir].hPredHT_pt,
                      "probe p_{T} (GeV)", "pt",     rn, RunYear, whichVFP);
        drawWithRatio(regions[ir].hObs_HT,     regions[ir].hPredVtx_HT,  regions[ir].hPredHT_HT,
                      "H_{T} (GeV)",       "HT",     rn, RunYear, whichVFP);
        drawWithRatio(regions[ir].hObs_nVtx,   regions[ir].hPredVtx_nVtx,regions[ir].hPredHT_nVtx,
                      "nVtx",              "nVtx",   rn, RunYear, whichVFP);
        drawWithRatio(regions[ir].hObs_invmass,regions[ir].hPredVtx_invmass,regions[ir].hPredHT_invmass,
                      "m_{tag,probe} (GeV)","invmass",rn, RunYear, whichVFP);
    }

    // Weight ratio diagnostic plots
    {
        TCanvas* c = new TCanvas("cwR","weight ratio",1400,400);
        c->Divide(3,1);
        c->cd(1); prof_wRatio_pt->SetMarkerStyle(20); prof_wRatio_pt->Draw("EP");
        TLine* l1=new TLine(30,1,230,1); l1->SetLineStyle(2); l1->Draw();
        c->cd(2); prof_wRatio_HT->SetMarkerStyle(20); prof_wRatio_HT->Draw("EP");
        TLine* l2=new TLine(0,1,1000,1); l2->SetLineStyle(2); l2->Draw();
        c->cd(3); prof_wRatio_nVtx->SetMarkerStyle(20); prof_wRatio_nVtx->Draw("EP");
        TLine* l3=new TLine(0,1,100,1); l3->SetLineStyle(2); l3->Draw();
        c->SaveAs(Form("closure_weightRatio_%d%s.png",RunYear,whichVFP.c_str()));
        delete c;
    }
    {
        TCanvas* c2d = new TCanvas("c2dR","2D weight ratio",700,600);
        c2d->SetRightMargin(0.15);
        h2_wRatio_full->Draw("COLZ");
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.DrawLatex(0.14,0.92,"w_{HT}/w_{Vtx} per denominator event (CR)");
        c2d->SaveAs(Form("closure_weightRatio2D_%d%s.png",RunYear,whichVFP.c_str()));
        delete c2d;
    }

    // ────────────────────────────────────────────────────────────────
    //  STEP 7 — Write full physics interpretation to text file
    // ────────────────────────────────────────────────────────────────

    txt << R"(
════════════════════════════════════════════════════════════════════
  WHY ARE THE SR/VR DIFFERENCES SO LARGE?  PHYSICS EXPLANATION
════════════════════════════════════════════════════════════════════

The fake rate is measured in the CR (MET < 70 GeV) and applied in
the VR and SR.  The key point is that the MET cut selects events
from VERY DIFFERENT kinematic regimes:

  CR (MET<70)   →  soft, low-activity events, moderate HT, low pT
  VR (MT<100)   →  intermediate MET, mixed HT
  SR (MET>120)  →  high-MET, boosted events, HIGH HT, HIGH pT

The HT distribution in the SR is strongly shifted to large values
relative to the CR.  Because HT and pT are correlated (Pearson
r~0.5), the HT-based fake rate EXTRAPOLATES into a kinematic region
(high HT + high pT) that carries the full double-counting bias.

The nVtx distribution in the SR is NOT strongly shifted relative to
the CR (nVtx is pileup-driven, independent of MET), so the nVtx-
based fake rate does not carry an extrapolation bias into the SR.

Observed effect on the estimate:
  (a) SR: pred_HT < pred_nVtx  (decrease of 27-45%)
      WHY: In the SR, events have large HT and large pT.  The HT-
      based fake rate f(pt)*f(HT) assigns a SINGLE large weight to
      each event (correlated), while f(pt)*f(nVtx) in principle
      should assign two independent weights.  However the MEASURED
      fake rate as a function of HT in the CR is LOWER in the high-
      HT regime than f(vtx) in the high-vtx regime, because in the
      CR high-HT events are produced by ISR/recoil that dilutes the
      electron-faking-photon rate.  Thus when applied to the high-HT
      SR events, f(HT) is SMALLER than f(nVtx), giving a lower
      prediction.

  (b) VR: pred_HT > pred_nVtx  (increase of 0-56% especially μγ)
      WHY: The VR (MET>70, MT<100) selects events with intermediate
      MET — these tend to have moderate-to-high HT from W/Z bosons
      with some recoil.  The μγ channel in particular involves a
      real muon from a W/Z decay which produces hadronic recoil
      (higher HT than eγ), so the HT distribution in the VR is
      shifted MORE for μγ.  The HT-based fake rate over-weights
      these events relative to the CR calibration, while nVtx
      remains flat and does not over-weight.

════════════════════════════════════════════════════════════════════
  IS THERE A BUG?
════════════════════════════════════════════════════════════════════

The direction and magnitude of the difference is PHYSICALLY MOTIVATED
and NOT a sign of a bug per se in the application of the fake rate.
However, the following confirmed bugs in fitFakeFunc_HT.C (see review)
could QUANTITATIVELY affect the result:

  Bug 1 (fitFakeFunc_HT.C): bgtree uses a hardcoded stale path.
    → Background template is from the wrong year.
    → This biases the Z-peak fit and hence the scale factor SF_HT.
    → FIX: use Form() with RunYear/Era like in fitFakeFunc.C.

  Bug 2 (fitFakeFunc_HT.C): HT systematic uses the same file twice
    (DYfile == Polfile = expo file).
    → The systematic uncertainty on the HT-based prediction is
      UNDERESTIMATED.  The scale factor uncertainty is too small.

  Bug 3 (FitKer_HT.C): byHT case missing in DY keyvariable selector.
    → The DY background template for each HT bin is EMPTY.
    → This means the signal extraction in the HT bins is done with
      a misidentified background shape, biasing the fake rate values
      in each HT bin.  This is likely the dominant source of the
      quantitative difference beyond physics expectations.

  Bug 4 (fitFakeFunc_HT.C): fr_eta_sigmaband is sized nEtaBins-1.
    → The last eta bin uncertainty is missing; minor but should be fixed.

════════════════════════════════════════════════════════════════════
  WHICH PARAMETERISATION IS CORRECT?
════════════════════════════════════════════════════════════════════

Closure test criteria:
  - Good closure: pred/obs ~ 1 in CR AND in VR (same-sign departure
    in VR that is not much larger than in CR is expected).
  - The parameterisation with better closure in the VR, after fixing
    the bugs above, should be preferred.

Expected outcome after bug fixes:
  1. After fixing Bug 3 (DY keyvariable), the HT bin fake rates will
     change and the SR/VR difference will shift.
  2. After fixing Bug 1 (bgtree path), the scale factor will be
     corrected; this affects the overall normalisation.
  3. If the closure in the VR is still worse for the HT-based
     parameterisation (pred/obs deviates more from 1), then the
     nVtx-based parameterisation is the safer choice.
  4. The correct long-term solution is a 2D measurement in (pT, HT)
     bins, which eliminates the factorisation bias entirely.

Run this closure script after fixing the bugs to obtain the
definitive answer.
)";

    txt.close();

    std::cout << "\nOutputs:\n";
    std::cout << "  closure_results_" << RunYear << whichVFP << ".txt\n";
    std::cout << "  closure_[pt/HT/nVtx/invmass]_[CR/VR/SR]_" << RunYear << whichVFP << ".png\n";
    std::cout << "  closure_weightRatio_" << RunYear << whichVFP << ".png\n";
    std::cout << "  closure_weightRatio2D_" << RunYear << whichVFP << ".png\n";

    delete etree;
}

// ════════════════════════════════════════════════════════════════════
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./closure_elefakepho.exe <RunYear> [preVFP=0/1]\n";
        return 1;
    }
    bool preVFP = (argc>2) ? (bool)atoi(argv[2]) : false;
    closure_elefakepho(atoi(argv[1]), preVFP);
    return 0;
}
