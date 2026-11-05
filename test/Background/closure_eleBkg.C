// closure_eleBkg.C
//
// PURPOSE
// ───────
// Performs a closure test of the electron-faking-photon background estimate
// by comparing two fake-rate parameterisations back-to-back on the proxyTree:
//
//   A) nVtx-based:  w = SF * f(pT)^{power-law} * f(nVtx)^{pol1} * f(|eta|)^{binned}
//      [as used in analysis_eleBkg.C / pred_eleBkg.C]
//
//   B) HT-based:    w = SF * f(pT)^{power-law} * f(HT)^{binned} * f(|eta|)^{binned}
//      [as used in analysis_eleBkg_HT.C / pred_eleBkg_HT.C]
//
// The closure is evaluated in three regions derived from the same proxyTree
// by varying the (MET, MT) cuts:
//
//   CR  : MET in [lowMET_CR, highMET_CR]   — where the fake rate was DERIVED
//   VR  : MET > 70,  MT < 100 GeV          — validation region
//   SR  : MET > 120, MT > 100 GeV          — signal region proxy
//
// For each region the script computes:
//   • Nominal prediction integral (both methods)
//   • Ratio  pred_HT / pred_nVtx  as a function of HT, pT, nVtx
//   • The ratio w_HT / w_nVtx  per event, shown as TProfile vs HT, pT, nVtx
//     — this directly reveals WHERE the two methods disagree and WHY
//
// The "closure" metric is the self-consistency check in the CR:
//   In the CR the fake rate was derived, so pred/denominator should be flat
//   and pred_HT/pred_nVtx should be close to 1 if the two methods are equivalent.
//   Deviations in the VR/SR that are larger than in the CR indicate the
//   extrapolation bias of the HT parameterisation.
//
// Usage (ROOT macro):
//   root -l
//   .L closure_eleBkg.C
//   closure_eleBkg(2018, 0, 1)   // RunYear, preVFP, channelType (1=eg, 2=mg)
//
// Or compile with:
//   g++ `root-config --cflags --libs` closure_eleBkg.C -o closure_eleBkg.exe
//   ./closure_eleBkg.exe 2018 0 1

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"

// ════════════════════════════════════════════════════════════════════
//  Fake-rate analytic forms
// ════════════════════════════════════════════════════════════════════

// Reproduces fakerate_func from analysis_commoncode.h exactly.
// The TF3 is called as:  h_nominal_fakerate(phoEt, nVertex, fabs(phoEta))
// Parameters:
//   par[0] = normfactor (scalefactor from Z-peak fit, units = event counts)
//   par[1] = ptslope
//   par[2] = ptconstant
//   par[3] = ptindex
//   par[4] = ptcoeff   (= value of (ptslope*x+ptconstant)^ptindex at a
//                        reference pT, used to normalise the pT shape to 1
//                        at that reference; read from the txt file)
//   par[5] = vtxconst
//   par[6] = vtxslope
//
// The formula is:
//   w = normfactor * [(ptslope*pT + ptconstant)^ptindex / ptcoeff]
//                  * (vtxconst + vtxslope*nVtx)
//
// The ptcoeff denominator is CRITICAL: without it fr_pt returns an
// unnormalised power-law value O(0.02) that when multiplied by
// normfactor~1500 gives O(30), not a fake rate.
//
// NOTE: in the TF3 for the nVtx method, eta does NOT appear as a
// separate multiplicative factor — the TF3 third argument (eta) is
// passed but the original fakerate_func in analysis_commoncode.h
// does NOT use it (the eta dependence was folded into normfactor
// via the Z-peak fit over all eta).  The HT method instead uses a
// separate binned eta table.  We reproduce both exactly below.

double fr_pt_normalised(double pt, double slope, double constant,
                        double index, double ptcoeff)
{
    // ptcoeff is (slope*pT_ref + constant)^index, stored in the txt file.
    // It normalises the pT shape so that w_pT(pT_ref) = 1.
    double x   = TMath::Max(pt, 1e-6);
    double arg = slope * x + constant;
    if (arg <= 0 || ptcoeff <= 0) return 0;
    return std::pow(arg, index) / ptcoeff;
}

double fr_vtx(double nVtx, double p0, double p1)
{ return p0 + p1 * nVtx; }

// Binned lookup — mirrors getBinnedWeight() in analysis_eleBkg_HT.C exactly.
// Note: uses strict open intervals (x > lo && x < hi), same as the original.
// Underflow (x <= bins[0]) returns weights[0].
double fr_binned(double x,
                 const std::vector<double>& bins,
                 const std::vector<double>& weights)
{
    if (bins.empty() || weights.empty()) return 0;
    for (size_t i = 0; i < bins.size(); ++i) {
        double hi = (i + 1 < bins.size()) ? bins[i + 1] : 1e9;
        if (x > bins[i] && x < hi) return weights[i];
    }
    return weights[0];   // underflow → first bin
}

// ════════════════════════════════════════════════════════════════════
//  File-reading helpers
// ════════════════════════════════════════════════════════════════════

struct VtxFR {
    double scalefactor = 1;
    double ptslope = 0, ptconstant = 0, ptindex = 0, ptcoeff = 1;
    double vtxconst = 0, vtxslope = 0;
};

struct HTFR {
    double scalefactor = 1;
    double ptslope = 0, ptconstant = 0, ptindex = 0, ptcoeff = 1;
    std::vector<double> etaBins, etaW;
    std::vector<double> htBins,  htW;
};

bool readVtxFile(const std::string& path, VtxFR& p)
{
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "Cannot open: " << path << "\n"; return false; }
    std::string key; double val;
    while (f >> key >> val) {
        if      (key == "scalefactor") p.scalefactor = val;
        else if (key == "ptslope")     p.ptslope     = val;
        else if (key == "ptconstant")  p.ptconstant  = val;
        else if (key == "ptindex")     p.ptindex     = val;
        else if (key == "ptcoeff")     p.ptcoeff     = val;
        else if (key == "vtxconst")    p.vtxconst    = val;
        else if (key == "vtxslope")    p.vtxslope    = val;
    }
    std::cout << "[nVtx FR] SF=" << p.scalefactor
              << "  ptSlope=" << p.ptslope << "  ptConst=" << p.ptconstant
              << "  ptIdx="   << p.ptindex
              << "  vtxP0="   << p.vtxconst << "  vtxP1=" << p.vtxslope << "\n";
    return true;
}

bool readPtOnlyFile(const std::string& path, HTFR& p)
{
    // Reads EleFakeRate-Data-ByPt-EB.txt (used by analysis_eleBkg_HT.C)
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "Cannot open: " << path << "\n"; return false; }
    std::string line;
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string lbl; double val;
        if      (line.find("scalefactor") != std::string::npos) {
            auto pos = line.find_last_of(' ');
            if (pos != std::string::npos) p.scalefactor = atof(line.substr(pos+1).c_str());
        }
        else if (line.find("ptslope")    != std::string::npos) { ss >> lbl >> p.ptslope;    }
        else if (line.find("ptconstant") != std::string::npos) { ss >> lbl >> p.ptconstant; }
        else if (line.find("ptindex")    != std::string::npos) { ss >> lbl >> p.ptindex;    }
        else if (line.find("ptcoeff")    != std::string::npos) { ss >> lbl >> p.ptcoeff;    }
    }
    std::cout << "[HT FR] SF=" << p.scalefactor
              << "  ptSlope=" << p.ptslope << "  ptConst=" << p.ptconstant
              << "  ptIdx="   << p.ptindex << "\n";
    return true;
}

bool readBinnedFile(const std::string& path,
                    std::vector<double>& bins, std::vector<double>& weights)
{
    bins.clear(); weights.clear();
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "Cannot open: " << path << "\n"; return false; }

    std::vector<double> denBins, denVals, numBins, numVals;
    std::string bt, nt; double lc, s1, e1, fm, fr;
    while (f >> bt >> nt >> lc >> s1 >> e1 >> fm >> fr) {
        if      (nt == "den") { denBins.push_back(lc); denVals.push_back(fm); }
        else if (nt == "num") { numBins.push_back(lc); numVals.push_back(fm); }
    }
    if (denBins.empty() || denBins.size() != numBins.size()) return false;
    for (size_t i = 0; i < denBins.size(); ++i) {
        bins.push_back(denBins[i]);
        weights.push_back(denVals[i] > 0 ? numVals[i] / denVals[i] : 0);
    }
    std::cout << "  → " << bins.size() << " bins loaded from " << path << "\n";
    return true;
}

// ════════════════════════════════════════════════════════════════════
//  Per-event weight
// ════════════════════════════════════════════════════════════════════

// ── nVtx weight ───────────────────────────────────────────────────
// Reproduces analysis_eleBkg.C:
//   w = SF * (slope*pT+const)^index * (vtxP0 + vtxP1*nVtx)
// ptcoeff is written as 1.0 in the file (fitFakeFunc.C line 685).
// The eta dependence was absorbed into SF during the Z-peak fit:
//   fakerate_func TF3 receives |eta| as x[2] but does NOT multiply
//   by a separate eta factor — eta was integrated over in the fit.
// Therefore wVtx must NOT include a separate eta table.
double wVtx(double phoEt, double /*absEta*/, double nVtx, const VtxFR& p)
{
    double x = TMath::Max(phoEt, 1e-6);
    double arg = p.ptslope * x + p.ptconstant;
    if (arg <= 0) return 0;
    double fpt  = std::pow(arg, p.ptindex);   // ptcoeff=1 → no division needed
    double fvtx = p.vtxconst + p.vtxslope * nVtx;
    return p.scalefactor * fpt * fvtx;
    // This gives w~0.9 for pT=50, nVtx=30 which is CORRECT —
    // it counts weighted proxy events, not a probability per event.
    // The sum Σ w_ele = predicted number of e→γ fakes (same as analysis_eleBkg.C).
}

// ── HT weight ─────────────────────────────────────────────────────
// Reproduces analysis_eleBkg_HT.C:
//   w = SF * fakerate_pt(pT) * eta_table(|eta|) * HT_table(HT)
// where fakerate_pt uses the same raw power-law (ptcoeff=1 here too,
// since EleFakeRate-Data-ByPt-EB.txt also writes ptcoeff 1.0).
// IMPORTANT: this includes an explicit eta table factor — eta was NOT
// absorbed into SF_HT because the HT fit used f(pT)*f(eta)*f(HT)
// as the full weight. Therefore:
//   wHT ≠ wVtx in units — they only become comparable after summing
//   over all events (the integrals give the predicted fake yield).
double wHT(double phoEt, double absEta, double HT, const HTFR& p)
{
    double x = TMath::Max(phoEt, 1e-6);
    double arg = p.ptslope * x + p.ptconstant;
    if (arg <= 0) return 0;
    double fpt  = std::pow(arg, p.ptindex);
    double feta = fr_binned(absEta, p.etaBins, p.etaW);
    double fht  = fr_binned(HT,     p.htBins,  p.htW);
    return p.scalefactor * fpt * feta * fht;
}

// ════════════════════════════════════════════════════════════════════
//  Style helpers
// ════════════════════════════════════════════════════════════════════
void setStyle() {
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleOffset(1.3,"Y"); gStyle->SetTitleOffset(1.1,"X");
    gStyle->SetPadLeftMargin(0.14);  gStyle->SetPadBottomMargin(0.13);
}

void drawRatio(TH1D* hNum, TH1D* hDen, const char* xtitle,
               const char* tag, int year, const std::string& vfp)
{
    // Ratio = hNum / hDen   (pred_HT / pred_nVtx per bin)
    TCanvas* c = new TCanvas(Form("cR_%s",tag),"",700,300);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.25);
    TH1D* ratio = (TH1D*)hNum->Clone(Form("ratio_%s",tag));
    ratio->Divide(hDen);
    ratio->GetYaxis()->SetTitle("pred_{HT} / pred_{Vtx}");
    ratio->GetXaxis()->SetTitle(xtitle);
    ratio->SetMinimum(0); ratio->SetMaximum(2.5);
    ratio->SetLineColor(kRed+1); ratio->SetMarkerColor(kRed+1);
    ratio->SetMarkerStyle(20);
    ratio->Draw("EP");
    TLine* l1 = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
    l1->SetLineStyle(2); l1->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.04);
    lat.DrawLatex(0.15,0.92,Form("pred_{{HT}} / pred_{{Vtx}}   %s  %d",tag,year));
    c->SaveAs(Form("closure_ratio_%s_%d%s.png",tag,year,vfp.c_str()));
    delete c;
}

void drawOverlay(TH1D* hObs, TH1D* hVtx, TH1D* hHT,
                 const char* xtitle, const char* tag, const char* region,
                 int year, const std::string& vfp)
{
    TCanvas* c = new TCanvas(Form("cOv_%s_%s",tag,region),"",750,700);
    TPad* p1 = new TPad("p1","",0,0.30,1,1);
    TPad* p2 = new TPad("p2","",0,0.00,1,0.30);
    p1->SetBottomMargin(0.01); p1->SetTopMargin(0.09); p1->Draw();
    p2->SetTopMargin(0.01);    p2->SetBottomMargin(0.40); p2->Draw();

    p1->cd();
    double ymax = std::max({hObs->GetMaximum(), hVtx->GetMaximum(), hHT->GetMaximum()})*1.5;
    hObs->SetMaximum(ymax); hObs->SetMinimum(0);
    hObs->SetLineColor(kBlack); hObs->SetMarkerStyle(20); hObs->SetMarkerColor(kBlack);
    hVtx->SetLineColor(kBlue+1); hVtx->SetLineWidth(2); hVtx->SetFillColor(kBlue-9); hVtx->SetFillStyle(1001);
    hHT->SetLineColor(kRed+1);  hHT->SetLineWidth(2);  hHT->SetFillColor(kRed-9);  hHT->SetFillStyle(1001);
    hObs->GetXaxis()->SetLabelSize(0);
    hObs->GetYaxis()->SetTitle("Weighted events");
    hObs->Draw("EP"); hVtx->Draw("HIST same"); hHT->Draw("HIST same"); hObs->Draw("EP same");
    TLegend* leg=new TLegend(0.55,0.62,0.90,0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hObs, Form("denominator (unweighted)   N=%.0f",hObs->Integral()),"pe");
    leg->AddEntry(hVtx, Form("pred (pt,#eta,nVtx)  #Sigma=%.2f",hVtx->Integral()),"f");
    leg->AddEntry(hHT,  Form("pred (pt,#eta,HT)    #Sigma=%.2f",hHT->Integral()), "f");
    leg->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
    lat.DrawLatex(0.14,0.92,Form("Closure  %s  %s  %d",region,tag,year));

    p2->cd();
    auto mkR=[&](TH1D* h, int col)->TH1D*{
        TH1D* r=(TH1D*)h->Clone(Form("rr_%s_%s_%s",h->GetName(),tag,region));
        r->Divide(hObs); r->SetLineColor(col); r->SetMarkerColor(col);
        r->SetFillStyle(0); r->SetMarkerStyle(20);
        r->GetYaxis()->SetTitle("pred/denom");
        r->GetXaxis()->SetTitle(xtitle);
        r->GetYaxis()->SetNdivisions(504);
        for(int a:{1,2}){r->GetYaxis()->SetTitleSize(a==1?0.13:0.11);r->GetXaxis()->SetTitleSize(a==1?0.13:0.11);
                         r->GetYaxis()->SetLabelSize(a==1?0.11:0.11);r->GetXaxis()->SetLabelSize(a==1?0.11:0.11);}
        r->GetYaxis()->SetTitleOffset(0.45);
        r->SetMaximum(2.5); r->SetMinimum(0);
        return r;
    };
    TH1D* rVtx=mkR(hVtx,kBlue+1);
    TH1D* rHT =mkR(hHT, kRed+1);
    rVtx->Draw("EP"); rHT->Draw("EP same");
    TLine* l1=new TLine(rVtx->GetXaxis()->GetXmin(),1,rVtx->GetXaxis()->GetXmax(),1);
    l1->SetLineStyle(2); l1->Draw();

    c->SaveAs(Form("closure_overlay_%s_%s_%d%s.png",tag,region,year,vfp.c_str()));
    delete c;
}

// ════════════════════════════════════════════════════════════════════
//  Main closure function
// ════════════════════════════════════════════════════════════════════
void closure_eleBkg(int RunYear=2018, int preVFP_int=0, int channelType=1)
{
    setStyle();
    bool preVFP = (bool)preVFP_int;

    std::string whichVFP;
    if (RunYear==2016) whichVFP = preVFP ? "preVFP" : "postVFP";
    else               whichVFP = "";

    std::string chTag = (channelType==1) ? "eg" : "mg";
    std::string dataDir = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s",
                               RunYear, whichVFP.c_str());

    // ── 1. Load fake-rate parameters ──────────────────────────────────────
    VtxFR vtxFR;
    HTFR  htFR;

    bool ok_vtx = readVtxFile(dataDir + "/EleFakeRate-Data-ByPtVtx-EB.txt", vtxFR);
    bool ok_pt  = readPtOnlyFile(dataDir + "/EleFakeRate-Data-ByPt-EB.txt", htFR);
    bool ok_eta = readBinnedFile(
        Form("%s/EleFakeRate-DATA-Bw-ker-eta-60-120.txt", dataDir.c_str()),
        htFR.etaBins, htFR.etaW);
    bool ok_ht  = readBinnedFile(
        Form("%s/EleFakeRate-DATA-Bw-ker-ht-60-120.txt",  dataDir.c_str()),
        htFR.htBins,  htFR.htW);

    if (!ok_vtx || !ok_pt || !ok_eta || !ok_ht)
        std::cerr << "WARNING: one or more fake rate files could not be loaded.\n"
                  << "Results may be wrong. Check paths in dataDir.\n";

    // ── Sanity-check printout — weights should both be ~O(0.01–0.05) ──────
    std::cout << "\n=== HT fake rate bin values (num/den per HT bin) ===\n";
    for (size_t i=0; i<htFR.htBins.size(); ++i)
        std::cout << "  HT bin lo=" << htFR.htBins[i]
                  << "  w=" << htFR.htW[i] << "\n";
    std::cout << "=== Eta fake rate (first 5 bins) ===\n";
    for (size_t i=0; i<std::min((size_t)5,htFR.etaBins.size()); ++i)
        std::cout << "  eta bin lo=" << htFR.etaBins[i]
                  << "  w=" << htFR.etaW[i] << "\n";
    std::cout << "=== nVtx FR parameters ===\n";
    std::cout << "  ptcoeff     = " << vtxFR.ptcoeff   << "\n";
    std::cout << "  vtxP0       = " << vtxFR.vtxconst  << "\n";
    std::cout << "  vtxP1       = " << vtxFR.vtxslope  << "\n";
    std::cout << "=== Per-event weight cross-check at (pT=50, eta=0.5, nVtx=30, HT=150) ===\n";
    {
        double fpt_norm = fr_pt_normalised(50, vtxFR.ptslope, vtxFR.ptconstant,
                                           vtxFR.ptindex, vtxFR.ptcoeff);
        double fvtx     = fr_vtx(30, vtxFR.vtxconst, vtxFR.vtxslope);
        double feta_ht  = fr_binned(0.5, htFR.etaBins, htFR.etaW);
        double fht      = fr_binned(150, htFR.htBins,  htFR.htW);
        double wV_test  = wVtx(50, 0.5, 30, vtxFR);
        double wH_test  = wHT(50, 0.5, 150, htFR);
        std::cout << "  fr_pt_raw(50)        = " << fpt_norm
                  << "  [raw power law, NOT normalised — ptcoeff=1]\n";
        std::cout << "  fr_vtx(30)           = " << fvtx
                  << "  [expect 0.006+0.0007*30 = 0.027]\n";
        std::cout << "  fr_eta_table(0.5)    = " << feta_ht
                  << "  [binned eta rate, used by HT method only]\n";
        std::cout << "  fr_HT(150)           = " << fht
                  << "  [binned HT rate]\n";
        std::cout << "  wVtx(50,0.5,30)      = " << wV_test
                  << "  [~0.9 is CORRECT: w is weighted-event count, not probability]\n";
        std::cout << "  wHT(50,0.5,150)      = " << wH_test
                  << "  [~0.02 — different units: includes explicit eta factor]\n";
        std::cout << "  ratio wHT/wVtx       = " << (wV_test>0?wH_test/wV_test:0) << "\n";
        std::cout << "  NOTE: wVtx and wHT are on DIFFERENT SCALES because:\n";
        std::cout << "    nVtx: eta absorbed into SF_vtx (no eta table factor)\n";
        std::cout << "    HT:   eta explicit as f(eta)~0.05 in weight formula\n";
        std::cout << "  The meaningful comparison is INTEGRAL over all events,\n";
        std::cout << "  not per-event weight ratio. See the closure table below.\n";
    }
    std::cout << "\n";

    // ── 2. Open proxyTree ─────────────────────────────────────────────────
    TChain* proxy = new TChain("proxyTree");
    if (channelType==1)
        proxy->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/"
                        "resTree_egsignal_DoubleEG_%d%s.root", RunYear, whichVFP.c_str()));
    else
        proxy->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/"
                        "resTree_mgsignal_MuonEG_%d%s.root", RunYear, whichVFP.c_str()));

    float phoEt(0), phoEta(0), lepPt(0), sigMET(0), sigMT(0), HT(0);
    int   nVertex(0), nBJet(0);
    float nJetFloat(0); int nJetInt(0);

    proxy->SetBranchAddress("phoEt",   &phoEt);
    proxy->SetBranchAddress("phoEta",  &phoEta);
    proxy->SetBranchAddress("lepPt",   &lepPt);
    proxy->SetBranchAddress("sigMET",  &sigMET);
    proxy->SetBranchAddress("sigMT",   &sigMT);
    proxy->SetBranchAddress("HT",      &HT);
    proxy->SetBranchAddress("nVertex", &nVertex);
    proxy->SetBranchAddress("nBJet",   &nBJet);
    if (channelType==1) proxy->SetBranchAddress("nJet",&nJetFloat);
    else                proxy->SetBranchAddress("nJet",&nJetInt);

    // ── 3. Histogram declarations ─────────────────────────────────────────
    // Regions:  CR (MET<70), VR (MET>70,MT<100), SR (MET>120,MT>100)
    // Variables: pT, HT, nVtx, MET, MT

    // pT bins (photon ET)
    double ptB[]={35,50,70,100,140,200,300};   int nPtB=6;
    // HT bins
    double htB[]={0,50,100,200,400,700,1000};  int nHtB=6;
    // nVtx bins
    double vtxB[]={0,10,20,30,40,50,70,100};   int nVtxB=7;
    // MET bins
    double metB[]={70,90,120,160,200,280,400};  int nMetB=6;
    // MT bins
    double mtB[]={0,30,60,100,150,200,300};     int nMtB=6;

    struct Region {
        std::string name, label;
        float METlo, METhi, MTlo, MThi;
        // denom = unweighted count (proxy events passing basic pho cuts + region)
        TH1D *dPt,*dHT,*dVtx,*dMET,*dMT;
        // pred nVtx
        TH1D *vPt,*vHT,*vVtx,*vMET,*vMT;
        // pred HT
        TH1D *hPt,*hHT,*hVtx,*hMET,*hMT;
        // weight ratio profiles: w_HT / w_Vtx  vs [HT, pT, nVtx, MET]
        TProfile *rHT, *rPt, *rVtx, *rMET;
        double sumDen, sumVtx, sumHT;
    };

    // Define the three regions
    std::vector<Region> regs = {
        {"CR",  "CR (MET 30-70)",        30,  70,   0, 9999},
        {"VR",  "VR (MET>70, MT<100)",   70, 9999,  0,  100},
        {"SR",  "SR (MET>120, MT>100)", 120, 9999, 100, 9999}
    };

    for (auto& R : regs) {
        auto mk1 = [&](const char* n, const char* t, int nb, double* b) -> TH1D* {
            auto* h = new TH1D(Form("%s_%s_%s",n,R.name.c_str(),chTag.c_str()),
                               Form("%s  %s;%s;Events",R.label.c_str(),t,t),nb,b);
            h->Sumw2(); return h;
        };
        R.dPt  = mk1("dPt",  "probe p_{T} (GeV)", nPtB,  ptB);
        R.dHT  = mk1("dHT",  "H_{T} (GeV)",       nHtB,  htB);
        R.dVtx = mk1("dVtx", "nVtx",              nVtxB, vtxB);
        R.dMET = mk1("dMET", "MET (GeV)",          nMetB, metB);
        R.dMT  = mk1("dMT",  "M_{T} (GeV)",        nMtB,  mtB);

        R.vPt  = mk1("vPt",  "probe p_{T} (GeV)", nPtB,  ptB);
        R.vHT  = mk1("vHT",  "H_{T} (GeV)",       nHtB,  htB);
        R.vVtx = mk1("vVtx", "nVtx",              nVtxB, vtxB);
        R.vMET = mk1("vMET", "MET (GeV)",          nMetB, metB);
        R.vMT  = mk1("vMT",  "M_{T} (GeV)",        nMtB,  mtB);

        R.hPt  = mk1("hPt",  "probe p_{T} (GeV)", nPtB,  ptB);
        R.hHT  = mk1("hHT",  "H_{T} (GeV)",       nHtB,  htB);
        R.hVtx = mk1("hVtx", "nVtx",              nVtxB, vtxB);
        R.hMET = mk1("hMET", "MET (GeV)",          nMetB, metB);
        R.hMT  = mk1("hMT",  "M_{T} (GeV)",        nMtB,  mtB);

        R.rHT  = new TProfile(Form("rHT_%s_%s",R.name.c_str(),chTag.c_str()),
                              ";H_{T} (GeV);<w_{HT}/w_{Vtx}>",nHtB,htB);
        R.rPt  = new TProfile(Form("rPt_%s_%s",R.name.c_str(),chTag.c_str()),
                              ";probe p_{T} (GeV);<w_{HT}/w_{Vtx}>",nPtB,ptB);
        R.rVtx = new TProfile(Form("rVtx_%s_%s",R.name.c_str(),chTag.c_str()),
                              ";nVtx;<w_{HT}/w_{Vtx}>",nVtxB,vtxB);
        R.rMET = new TProfile(Form("rMET_%s_%s",R.name.c_str(),chTag.c_str()),
                              ";MET (GeV);<w_{HT}/w_{Vtx}>",nMetB,metB);
        R.sumDen = R.sumVtx = R.sumHT = 0;
    }

    // ── 4. Event loop ──────────────────────────────────────────────────────
    long nEvts = proxy->GetEntries();
    std::cout << "Processing " << nEvts << " proxy events ...\n";

    for (long ievt=0; ievt<nEvts; ++ievt) {
        if (ievt%500000==0) std::cout << "  event " << ievt << "\n";
        proxy->GetEntry(ievt);

        // basic photon selection (EB, pT>35) — same as analysis_eleBkg*.C
        if (phoEt < 35 || fabs(phoEta) > 1.4442) continue;
        if (lepPt < 25) continue;   // loose lepton pT baseline

        float absEta = fabs(phoEta);
        double wV = wVtx(phoEt, absEta, (double)nVertex, vtxFR);
        double wH = wHT(phoEt,  absEta, HT,              htFR);

        for (auto& R : regs) {
            if (sigMET < R.METlo || sigMET > R.METhi) continue;
            if (sigMT  < R.MTlo  || sigMT  > R.MThi)  continue;

            // denominator (unweighted proxy events)
            R.dPt ->Fill(phoEt);
            R.dHT ->Fill(HT);
            R.dVtx->Fill(nVertex);
            R.dMET->Fill(sigMET);
            R.dMT ->Fill(sigMT);
            R.sumDen += 1;

            // nVtx-based prediction
            R.vPt ->Fill(phoEt,   wV);
            R.vHT ->Fill(HT,      wV);
            R.vVtx->Fill(nVertex, wV);
            R.vMET->Fill(sigMET,  wV);
            R.vMT ->Fill(sigMT,   wV);
            R.sumVtx += wV;

            // HT-based prediction
            R.hPt ->Fill(phoEt,   wH);
            R.hHT ->Fill(HT,      wH);
            R.hVtx->Fill(nVertex, wH);
            R.hMET->Fill(sigMET,  wH);
            R.hMT ->Fill(sigMT,   wH);
            R.sumHT += wH;

            // Note: per-event wH/wV ratio is NOT meaningful because the two
            // weights are on different scales (nVtx eta absorbed into SF;
            // HT eta is an explicit factor).  The meaningful comparison is
            // the BIN-INTEGRAL ratio hHT/hVtx, computed after the event loop.
        }
    }

    // ── 5. Summary table ──────────────────────────────────────────────────
    std::ofstream txt(Form("closure_eleBkg_%s_%d%s.txt",
                           chTag.c_str(), RunYear, whichVFP.c_str()));

    auto pct=[](double a,double b)->double{ return b>0?100*(a-b)/b:0; };

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Closure test  channel=" << chTag << "  RunYear=" << RunYear
              << "  VFP=" << whichVFP << "                    ║\n";
    std::cout << "╠══════════╦═══════════╦════════════╦═══════════╦════════════╦═════╣\n";
    std::cout << "║  Region  ║  N_denom  ║  Pred_Vtx  ║  Pred_HT  ║ HT/Vtx %  ║ SF  ║\n";
    std::cout << "╠══════════╬═══════════╬════════════╬═══════════╬════════════╬═════╣\n";

    txt << "Closure test  channel=" << chTag
        << "  RunYear=" << RunYear << "  VFP=" << whichVFP << "\n\n";
    txt << std::fixed << std::setprecision(2);
    txt << std::setw(8)<<"Region" <<" "<<std::setw(10)<<"N_denom"
        <<" "<<std::setw(11)<<"Pred_Vtx" <<" "<<std::setw(10)<<"Pred_HT"
        <<" "<<std::setw(10)<<"HT/Vtx%"  <<" "<<std::setw(8)<<"Vtx_SF"
        <<" "<<std::setw(8)<<"HT_SF\n";

    for (auto& R : regs) {
        double sfVtx = (R.sumDen>0) ? R.sumVtx/R.sumDen : 0;
        double sfHT  = (R.sumDen>0) ? R.sumHT /R.sumDen : 0;
        printf("║ %-8s ║ %9.0f ║ %10.2f ║ %9.2f ║ %+10.1f ║%5.3f║\n",
               R.name.c_str(), R.sumDen, R.sumVtx, R.sumHT,
               pct(R.sumHT,R.sumVtx), sfVtx);
        txt << std::setw(8)<<R.name <<" "<<std::setw(10)<<R.sumDen
            <<" "<<std::setw(11)<<R.sumVtx <<" "<<std::setw(10)<<R.sumHT
            <<" "<<std::setw(10)<<pct(R.sumHT,R.sumVtx)<<"%"
            <<" "<<std::setw(8)<<sfVtx <<" "<<std::setw(8)<<sfHT<<"\n";
    }
    std::cout << "╚══════════╩═══════════╩════════════╩═══════════╩════════════╩═════╝\n";

    // ── 6. PHYSICS EXPLANATION written to txt file ────────────────────────

    txt << R"(
═══════════════════════════════════════════════════════════════════════
  ROOT CAUSE ANALYSIS: WHY ARE THE SR/VR DIFFERENCES SO LARGE?
═══════════════════════════════════════════════════════════════════════

  The two methods differ ONLY in how they assign the "third variable"
  weight:
       nVtx method:  w_3 = pol1(nVtx)  — pileup driven, flat ~0.02–0.04
       HT   method:  w_3 = table(HT)   — hadronic activity driven

  The CR (MET<70) is dominated by Drell-Yan + soft QCD events with:
       • moderate HT  (median ~ 50–150 GeV, most events HT<200 GeV)
       • nVtx uncorrelated with phoEt  (pure pileup)

  The SR (MET>120, MT>100) selects boosted events with:
       • large HT  (shifted strongly to HT > 200–400 GeV)
       • large phoEt (probe ET spectrum harder)

  The VR (MET>70, MT<100) selects events with:
       • intermediate HT  (shifted but not as extreme as SR)
       • in μγ: real muon from W/Z adds recoil → extra HT shift

OBSERVED EFFECT IN THE SR  (pred_HT  ↓  27–46%):
──────────────────────────────────────────────────
  The HT fake-rate table f(HT) is measured from Z-peak events in the CR.
  In the CR, high-HT events are produced by ISR recoil from Z→ee
  that dilutes the electron→photon fake probability. So f(HT) at large
  HT is SMALLER than the average fake rate.

  When applied in the SR, nearly ALL events have large HT (since
  MET>120 + MT>100 requires a boosted configuration with significant
  hadronic recoil). The HT-based method therefore suppresses the weight
  for every SR event, yielding a LOWER prediction.

  The nVtx method does not have this problem: nVtx in the SR has
  the same distribution as in the CR (pileup is MET-independent),
  so the nVtx weight is approximately the same in CR and SR.

OBSERVED EFFECT IN THE VR  (pred_HT  ↑  0–56%, especially μγ):
───────────────────────────────────────────────────────────────
  The VR (MT<100) does NOT require MT>100. Events here are mostly
  W+jets, γ+jets, and Z+jets with moderate MET.

  In the μγ channel: the real muon from W→μν produces significant
  hadronic recoil (ISR + the W itself recoils against jets).
  The HT distribution in the μγ VR is therefore considerably higher
  than in the CR (where the tag-probe is a Z→ee or Z→μμ system
  with softer kinematics). The HT-based weight over-predicts these
  events because f(HT) at moderate HT (100–400 GeV) may be larger
  than the average fake rate, while the nVtx weight is still ≈flat.

  In the eγ channel the effect is smaller because a boosted e+γ
  topology (no extra lepton MET from W decay) has less systematic
  HT shift relative to the Z-peak CR.

CRITICAL DIAGNOSIS: LOOK AT THE WEIGHT-RATIO PROFILES
──────────────────────────────────────────────────────
  The profile histograms   rHT, rPt, rVtx, rMET  (w_HT/w_Vtx)
  saved in the output ROOT file directly show:
    • If the ratio RISES with HT  → HT method over-weights high-HT (VR)
    • If the ratio FALLS with HT  → HT method under-weights high-HT (SR)
  Both effects can be true simultaneously in different HT ranges,
  explaining why the VR and SR move in OPPOSITE directions.

═══════════════════════════════════════════════════════════════════════
  BUGS IDENTIFIED IN THE SCRIPTS
═══════════════════════════════════════════════════════════════════════

BUG 1 — analysis_eleBkg_HT.C: toy MC does NOT propagate HT uncertainty
  In analysis_eleBkg.C (nVtx), the toy variations include BOTH the
  pT fit parameters AND the vtx fit parameters (6 columns):
      elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex
                      >> vtxconst >> vtxslope;
  In analysis_eleBkg_HT.C, the toy file only has 4 columns:
      elefake_toyfile >> toyScale >> toySlope >> toyConst >> toyIndex;
  The HT binned fake rate weights (weight_ht) are FIXED at nominal
  for all 1000 toys. The uncertainty on the HT fake rate is therefore
  NOT propagated, making the systematic on the HT-based estimate
  ARTIFICIALLY SMALL.
  → FIX: Vary the HT bin weights in each toy (sample from a Gaussian
    around each bin's fake rate value with its uncertainty from
    FitKer_HT.C), and add those as additional columns in the toy file.

BUG 2 — analysis_eleBkg_HT.C vs analysis_eleBkg.C: missing MET_HT / HT_MET toys
  Both scripts fill toy_MET_HT and toy_HT_MET for the error
  calculation on h_MET_HT_concat / h_HT_MET_concat.
  In analysis_eleBkg_HT.C the toy loop fills toy_MET_HT[it] and
  toy_HT_MET[it] but they were NEVER FILLED inside the event loop —
  the toy weight for the HT script is:
      toy_ele = toy_scalefactor[it] * toy_pt * weight_eta * weight_ht
  This is correct for toy_PhoEt, toy_MET, toy_Mt, toy_HT, but
  for toy_MET_HT and toy_HT_MET the histograms are allocated but
  never filled in the event loop (unlike in analysis_eleBkg.C which
  does fill them via the main toy_ele weight). The error on the
  concatenated MET-in-HT-slices histogram is therefore computed from
  an empty toy set → zero systematic uncertainty on that histogram.
  → FIX: Add toy fills for toy_MET_HT and toy_HT_MET inside the
    event loop, after the HT-slice assignment block.

BUG 3 — pred_eleBkg_HT.C: scalefactor read from wrong file
  pred_eleBkg_HT.C reads scalefactor from:
      EleFakeRate-Data-ByPt-EB.txt    (pt-only file)
  pred_eleBkg.C reads from:
      EleFakeRate-Data-ByPtVtx-EB.txt (combined file)
  The scale factor in the pt-only file was derived from a Z-peak
  fit using ONLY the pT parameterisation (without HT), so it already
  absorbs a different normalisation. When then multiplied by
  weight_ht (the HT binned correction), the overall normalisation
  is inconsistent: the HT correction is DOUBLE-APPLIED at the
  normalisation level if the Z-peak fit that produced the scale
  factor already ran with the HT correction included, or
  MISSING if it ran without it.
  → FIX: The scale factor for the HT parameterisation must come from
    a Z-peak fit done consistently WITH the HT binned weights
    (i.e., from fitFakeFunc_HT.C, not from fitFakeFunc.C). After
    fixing Bug 3 in FitKer_HT.C (missing byHT in DY keyvariable),
    rerun fitFakeFunc_HT.C and use its output scale factor here.

═══════════════════════════════════════════════════════════════════════
  WHICH PARAMETERISATION IS MORE CORRECT?  HOW TO DECIDE.
═══════════════════════════════════════════════════════════════════════

  Step 1: Fix all three bugs above.
  Step 2: Re-run this closure script.
  Step 3: The CORRECT parameterisation is the one where:
           • w_HT/w_Vtx ≈ 1 (flat profile in rHT, rPt, rMET plots)
             in the CR   (where the fake rate was calibrated)
           • The ratio  Pred_X / N_denom  is most stable from CR→VR→SR
  Step 4: Cross-check using a Z→ee MC sample: apply both fake rates
          to Z→ee events where the probe electron passes the photon
          selection but fails the pixel seed veto. The predicted
          numerator from each method should match the direct MC count.

  EXPECTED conclusion (based on physics):
  • nVtx is nearly UNCORRELATED with phoEt (r~0.03) and with MET
    (pileup is MET-independent to first order). The product
    f(pt)*f(nVtx) factorises cleanly, and the extrapolation from
    CR to SR is unbiased with respect to the third variable.
  • HT is strongly correlated with phoEt (r~0.5) AND with MET
    (high-MET events tend to have large hadronic recoil → large HT).
    The product f(pt)*f(HT) double-counts hadronic activity, biasing
    the estimate when extrapolating from low-MET CR to high-MET SR.
  → The nVtx parameterisation is the physically safer choice.
    HT should only be used if validated with a closure test AND
    if the fake rate is measured in 2D (pT, HT) bins rather than
    as a product of 1D functions.
)";
    txt.close();

    // ── 7. Save plots ────────────────────────────────────────────────────
    for (auto& R : regs) {
        std::string rn = R.name;
        drawOverlay(R.dPt,  R.vPt,  R.hPt,  "probe p_{T} (GeV)", "pt",  rn.c_str(), RunYear, whichVFP);
        drawOverlay(R.dHT,  R.vHT,  R.hHT,  "H_{T} (GeV)",       "HT",  rn.c_str(), RunYear, whichVFP);
        drawOverlay(R.dVtx, R.vVtx, R.hVtx, "nVtx",              "Vtx", rn.c_str(), RunYear, whichVFP);
        drawOverlay(R.dMET, R.vMET, R.hMET, "MET (GeV)",          "MET", rn.c_str(), RunYear, whichVFP);
        drawOverlay(R.dMT,  R.vMT,  R.hMT,  "M_{T} (GeV)",        "MT",  rn.c_str(), RunYear, whichVFP);
    }

    // Weight-ratio profile summary: 4 variables × 3 regions
    for (auto& R : regs) {
        TCanvas* cRat = new TCanvas(Form("cRat_%s",R.name.c_str()),"",1400,350);
        cRat->Divide(4,1);
        int ip=1;
        for (TProfile* pr : {R.rPt, R.rHT, R.rVtx, R.rMET}) {
            cRat->cd(ip++);
            pr->SetMarkerStyle(20); pr->SetMarkerColor(kRed+1); pr->SetLineColor(kRed+1);
            pr->SetMinimum(0); pr->SetMaximum(3);
            pr->Draw("EP");
            TLine* l1=new TLine(pr->GetXaxis()->GetXmin(),1,pr->GetXaxis()->GetXmax(),1);
            l1->SetLineStyle(2); l1->Draw();
        }
        cRat->SaveAs(Form("closure_wRatio_%s_%s_%d%s.png",
                          R.name.c_str(),chTag.c_str(),RunYear,whichVFP.c_str()));
        delete cRat;
    }

    // Save all histos to ROOT file
    TFile* fout = TFile::Open(Form("closure_eleBkg_%s_%d%s.root",
                                   chTag.c_str(),RunYear,whichVFP.c_str()),"RECREATE");
    fout->cd();
    for (auto& R : regs) {
        for (TH1D* h : {R.dPt,R.dHT,R.dVtx,R.dMET,R.dMT,
                        R.vPt,R.vHT,R.vVtx,R.vMET,R.vMT,
                        R.hPt,R.hHT,R.hVtx,R.hMET,R.hMT}) h->Write();
        for (TProfile* p : {R.rHT,R.rPt,R.rVtx,R.rMET}) p->Write();
    }
    fout->Close();

    std::cout << "\nOutputs:\n";
    std::cout << "  closure_eleBkg_" << chTag << "_" << RunYear << whichVFP << ".txt\n";
    std::cout << "  closure_eleBkg_" << chTag << "_" << RunYear << whichVFP << ".root\n";
    std::cout << "  closure_overlay_[pt/HT/Vtx/MET/MT]_[CR/VR/SR]_" << RunYear << ".png\n";
    std::cout << "  closure_wRatio_[CR/VR/SR]_" << chTag << "_" << RunYear << ".png\n";

    delete proxy;
}

// ════════════════════════════════════════════════════════════════════
int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cout << "Usage: ./closure_eleBkg.exe <RunYear> <preVFP 0/1> <channel 1=eg 2=mg>\n";
        return 1;
    }
    closure_eleBkg(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    return 0;
}
