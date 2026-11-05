// =============================================================================
// compare_elefakepho_methods.C
//
// Overlays the e->gamma fake photon background estimate (h_elefakepho_norm,
// combined into NBIN*2 = 36 search bins: mg in [1,NBIN], eg in
// [NBIN+1,2*NBIN]) from TWO fake-rate methods, read from two sets of
// pred_eleBkg output files that differ only by a filename suffix:
//
//   default method : /uscms_data/d3/tmishra/Output/signalTree_<channel>_eleBkg_<year><VFP>.root
//   new method     : /uscms_data/d3/tmishra/Output/signalTree_<channel>_eleBkg_<year><VFP>_HT_pT_nVtx.root
//
// where <channel> is "mg" (muon+gamma) or "egamma" (electron+gamma), exactly
// the file-naming/bin-combination convention used in plot_eventct.C (see
// its h_elefakepho_norm construction, ~line 601 & ~684). Only the
// e->gamma-fake histograms are read here (h_elefakepho_norm +
// h_elefakepho_norm's per-bin stat error) -- this script does not touch any
// of the other backgrounds/systematics plot_eventct.C also handles.
//
// Produces a two-pad canvas: top = both estimates overlaid (log-y, stat
// error bars), bottom = ratio (new method / default method) per bin, plus a
// per-bin printout to stdout.
// =============================================================================

// ----- Standard library includes -----
#include <sstream>
#include <string>
#include <iostream>

// ----- ROOT includes -----
#include <TROOT.h>
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"

// ----- Analysis framework includes -----
// Provides: setTDRStyle(), SetSignalConfig(), CMS_lumi(), setTopPad(),
//           and global variables RunYear, preVFP
#include "../../include/analysis_commoncode.h"

// total16 = true  -> treat 2016preVFP + 2016postVFP as a single 2016 dataset
// total16 = false -> use preVFP / postVFP split as set by analysis_commoncode
// (same flag/convention as plot_eventct.C)
bool total16 = false;

// -----------------------------------------------------------------------------
// Reads one method's e->gamma fake estimate and combines mg+eg into a single
// NBIN*2-bin histogram, exactly the way plot_eventct.C builds
// h_elefakepho_norm: bins [1,NBIN] <- mg channel, bins [NBIN+1,2*NBIN] <- eg
// channel. fileTag is appended right before ".root" (empty string for the
// default method, "_HT_pT_nVtx" for the new method).
// -----------------------------------------------------------------------------
TH1D* ReadCombinedEleFakeEstimate(int NBIN, int RunYear, const std::string &whichVFP,
                                   const std::string &fileTag, const std::string &histName)
{
    TFile *mgfile = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_eleBkg_%d%s%s.root",
                                      RunYear, whichVFP.c_str(), fileTag.c_str()));
    TFile *egfile = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_eleBkg_%d%s%s.root",
                                      RunYear, whichVFP.c_str(), fileTag.c_str()));
    if(!mgfile || mgfile->IsZombie()){
        std::cerr << "ERROR: could not open mg file for tag \"" << fileTag << "\"" << std::endl;
        return nullptr;
    }
    if(!egfile || egfile->IsZombie()){
        std::cerr << "ERROR: could not open eg file for tag \"" << fileTag << "\"" << std::endl;
        return nullptr;
    }

    TH1D *mg_norm = (TH1D*)mgfile->Get("mg_elefakepho_norm");
    TH1D *eg_norm = (TH1D*)egfile->Get("eg_elefakepho_norm");
    if(!mg_norm || !eg_norm){
        std::cerr << "ERROR: could not read mg_elefakepho_norm / eg_elefakepho_norm for tag \"" << fileTag << "\"" << std::endl;
        return nullptr;
    }

    TH1D *h = new TH1D(histName.c_str(), "", NBIN*2, 1, NBIN*2+1);
    for(int ibin=1; ibin<=NBIN; ibin++){
        h->SetBinContent(ibin,      mg_norm->GetBinContent(ibin));
        h->SetBinError(ibin,        mg_norm->GetBinError(ibin));
        h->SetBinContent(ibin+NBIN, eg_norm->GetBinContent(ibin));
        h->SetBinError(ibin+NBIN,   eg_norm->GetBinError(ibin));
    }
    h->SetDirectory(0); // detach from file so it survives mgfile/egfile being closed
    mgfile->Close();
    egfile->Close();
    return h;
}

// =============================================================================
// Main function
//   NBIN : number of signal bins per lepton channel (same convention as
//          plot_eventct.C; typically 18, giving 2*NBIN = 36 total bins)
// =============================================================================
void compare_elefakepho_methods(int NBIN = 18)
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    setTDRStyle();
    SetSignalConfig();
    gStyle->SetErrorX(0);
    gSystem->Load("../../lib/libAnaClasses.so");

    std::string whichVFP;
    if      (RunYear == 2016 && preVFP == 1) whichVFP = "preVFP";
    else if (RunYear == 2016 && preVFP == 0) whichVFP = "postVFP";
    else                                     whichVFP = "";
    if (RunYear == 2016 && total16 == true) whichVFP = "";

    std::cout << "Year: " << RunYear << "\t" << whichVFP << std::endl;

    // ----- Read both methods' combined (mg+eg) e->gamma fake estimates -----
    TH1D *h_default = ReadCombinedEleFakeEstimate(NBIN, RunYear, whichVFP, "",             "h_elefakepho_default");
    TH1D *h_newmethod = ReadCombinedEleFakeEstimate(NBIN, RunYear, whichVFP, "_HT_pT_nVtx", "h_elefakepho_HT_pT_nVtx");
    if(!h_default || !h_newmethod){
        std::cerr << "ERROR: failed to build one or both combined histograms -- aborting." << std::endl;
        return;
    }

    // ----- Per-bin printout -----
    std::cout << "\nbin\tdefault\t\tHT_pT_nVtx\tratio(new/default)\n";
    for(int ibin=1; ibin<=NBIN*2; ibin++){
        double d = h_default->GetBinContent(ibin);
        double n = h_newmethod->GetBinContent(ibin);
        std::cout << ibin << "\t" << d << "\t" << n << "\t" << (d>0 ? n/d : 0) << "\n";
    }

    // ----- Ratio graph: new method / default method, per bin -----
    TGraphErrors *g_ratio = new TGraphErrors();
    int nRatioPts = 0;
    for(int ibin=1; ibin<=NBIN*2; ibin++){
        double d  = h_default->GetBinContent(ibin);
        double n  = h_newmethod->GetBinContent(ibin);
        double de = h_default->GetBinError(ibin);
        double ne = h_newmethod->GetBinError(ibin);
        if(d <= 0) continue;
        double r  = n/d;
        // standard error propagation for a ratio of two (uncorrelated) estimates
        double re = r * sqrt( (de>0? (de/d)*(de/d) : 0) + (ne>0? (ne/n)*(ne/n) : 0) );
        g_ratio->SetPoint(nRatioPts, ibin+0.5, r);
        g_ratio->SetPointError(nRatioPts, 0.5, re);
        nRatioPts++;
    }

    // =========================================================================
    // Canvas: top pad overlays both estimates, bottom pad shows the ratio.
    // =========================================================================
    TCanvas *can = new TCanvas("can", "", 1200, 800);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.04);
    can->SetTopMargin(0.08);
    can->SetBottomMargin(0.1);
    can->cd();

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    setTopPad(pad1);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();

    h_default->SetLineColor(kBlack);
    h_default->SetMarkerColor(kBlack);
    h_default->SetMarkerStyle(20);
    h_default->SetLineWidth(2);
    h_default->GetXaxis()->SetRangeUser(1, NBIN*2+1);
    // Auto y-range with headroom, floor at 0.01 so log-scale doesn't choke on zero-content bins
    double ymax = std::max(h_default->GetMaximum(), h_newmethod->GetMaximum());
    h_default->GetYaxis()->SetRangeUser(0.01, ymax*10);
    h_default->GetYaxis()->SetTitle("e#rightarrow#gamma fake estimate (events / bin)");
    h_default->GetYaxis()->SetTitleOffset(1.4);
    h_default->GetXaxis()->SetTitle("");
    h_default->Draw("E1P");

    h_newmethod->SetLineColor(kRed+1);
    h_newmethod->SetMarkerColor(kRed+1);
    h_newmethod->SetMarkerStyle(21);
    h_newmethod->SetLineWidth(2);
    h_newmethod->Draw("E1P same");

    TLegend *leg = new TLegend(0.55, 0.80, 0.94, 0.90);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_default,   "Default: R = N #upoint f(p_{T}) #upoint f(N_{vtx}) #upoint f(#eta)", "ep");
    leg->AddEntry(h_newmethod, "New: R = f_{p_{T}}(p_{T} | N_{vtx} slice) #upoint [f(#eta)/#LTf(#eta)#GT]", "ep");
    leg->Draw("same");

    // ----- Dividing lines / labels between kinematic regions -----
    // Generalizes plot_eventct.C's hardcoded NBIN=18 dividers to arbitrary
    // (even) NBIN: each channel's NBIN bins split into two equal pT(gamma)
    // blocks (< 200 GeV, > 200 GeV).
    int half = NBIN/2;
    double lineTop = ymax*10;
    TLine *line_mg_pt  = new TLine(half+1,      0, half+1,      lineTop); line_mg_pt->SetLineStyle(2);
    TLine *line_mgeg   = new TLine(NBIN+1,      0, NBIN+1,      lineTop); line_mgeg->SetLineStyle(2);
    TLine *line_eg_pt  = new TLine(NBIN+half+1, 0, NBIN+half+1, lineTop); line_eg_pt->SetLineStyle(2);
    line_mg_pt->Draw("same");
    line_mgeg->Draw("same");
    line_eg_pt->Draw("same");

    TLatex *latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->DrawLatex(1.5,           ymax*3, "#mu#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(half+1.5,      ymax*3, "#mu#gamma, p_{T}^{#gamma} > 200 GeV");
    latex->DrawLatex(NBIN+1.5,      ymax*3, "e#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(NBIN+half+1.5, ymax*3, "e#gamma, p_{T}^{#gamma} > 200 GeV");

    if      (RunYear == 2016 && preVFP == 1)  CMS_lumi(pad1, 1, 1, 11);
    else if (RunYear == 2016 && preVFP == 0)  CMS_lumi(pad1, 2, 1, 11);
    else if (RunYear == 2017)                  CMS_lumi(pad1, 3, 1, 11);
    else if (RunYear == 2018)                  CMS_lumi(pad1, 4, 1, 11);
    else if (RunYear == 678)                   CMS_lumi(pad1, 5, 1, 11);
    pad1->RedrawAxis();

    // ----- Ratio pad -----
    can->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    TH1F *ratioframe = new TH1F("frame", ";Search region;#frac{HT_pT_nVtx}{Default}", NBIN*2, 1, NBIN*2+1);
    ratioframe->GetXaxis()->SetRangeUser(1, NBIN*2+1);
    ratioframe->GetYaxis()->SetTitleOffset(1.4);
    ratioframe->GetYaxis()->SetRangeUser(0, 2.0);
    std::ostringstream label;
    for(int i=1; i<=NBIN*2; i++){
        label.str("");
        label << i;
        ratioframe->GetXaxis()->SetBinLabel(i, label.str().c_str());
    }
    ratioframe->Draw();

    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerColor(kRed+1);
    g_ratio->SetLineColor(kRed+1);
    g_ratio->Draw("P same");

    TLine *flatratio = new TLine(1, 1, NBIN*2+1, 1);
    flatratio->SetLineStyle(2);
    flatratio->Draw("same");

    can->Update();
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/elefakepho_compare_%d%s.png", RunYear, whichVFP.c_str()));
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/elefakepho_compare_%d%s.pdf", RunYear, whichVFP.c_str()));
}
