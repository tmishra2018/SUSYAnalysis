#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "../../include/analysis_commoncode.h"   


struct SampleDef {
    TString procName;   // used to build the file name
    TString label;      // legend label
    TString colorHex;   // fill color
};

struct SignalDef {
    TString label;       // legend label
    TString filepath;    // full path to the signal histogram file
    TString histname;    // histogram key, WITHOUT the trailing _eg/_mg
    int     lineColor;
    int     lineStyle;
};

struct HistBundle {
    TH1D *hMET = nullptr;
    TH1D *hMT  = nullptr;
    TH2D *h2D  = nullptr;
};


HistBundle MakeBkgHistBundle(const char* tag, const TString& file, const char* treeName,
                              int channelType, int RunYear, bool preVFP, esfScaleFactor &objectESF)
{
    HistBundle b;
    b.hMET = new TH1D(Form("p_MET_%s", tag), "; p_{T}^{miss} (GeV); Events", nBkgMETBins, bkgMETBins);
    b.hMT  = new TH1D(Form("p_MT_%s",  tag), "; M_{T} (GeV); Events", nBkgMtBins, bkgMtBins);
    b.h2D  = new TH2D(Form("p_MET_MT_%s", tag), "; p_{T}^{miss} (GeV); M_{T} (GeV)",
                       nBkgMETBins, bkgMETBins, nBkgMtBins, bkgMtBins);
    b.hMET->Sumw2(); b.hMT->Sumw2(); b.h2D->Sumw2();

    TChain *tree = new TChain(treeName, treeName);
    tree->Add(file);
    if (tree->GetEntries() == 0) {
        std::cout << "[WARNING] no entries found for " << file << " -- skipping." << std::endl;
        return b;
    }

    float crosssection(0), ntotalevent(0), PUweight(1);
    float phoEt(0), phoEta(0), lepPt(0), lepEta(0), sigMET(0), sigMT(0), phoPhi(0);
    std::vector<int>   *mcPID = 0, *mcMomPID = 0;
    std::vector<float> *mcEta = 0, *mcPhi = 0, *mcPt = 0;

    tree->SetBranchAddress("crosssection", &crosssection);
    tree->SetBranchAddress("ntotalevent",  &ntotalevent);
    tree->SetBranchAddress("PUweight",     &PUweight);
    tree->SetBranchAddress("phoEt",        &phoEt);
    tree->SetBranchAddress("phoEta",       &phoEta);
    tree->SetBranchAddress("phoPhi",       &phoPhi);
    tree->SetBranchAddress("lepPt",        &lepPt);
    tree->SetBranchAddress("lepEta",       &lepEta);
    tree->SetBranchAddress("sigMET",       &sigMET);
    tree->SetBranchAddress("sigMT",        &sigMT);
    tree->SetBranchAddress("mcPID",    &mcPID);
    tree->SetBranchAddress("mcMomPID", &mcMomPID);
    tree->SetBranchAddress("mcEta",    &mcEta);
    tree->SetBranchAddress("mcPhi",    &mcPhi);
    tree->SetBranchAddress("mcPt",     &mcPt);

    for (Long64_t ievt = 0; ievt < tree->GetEntries(); ievt++) {
        tree->GetEntry(ievt);

        if (phoEt < 35 || fabs(phoEta) > 1.4442) continue;

        double scalefactor = 1.0;
        if (channelType == 1) {
            scalefactor = objectESF.getElectronESF(lepPt, lepEta) * objectESF.getPhotonESF(phoEt, phoEta)
                        * objectESF.getegPhotonTRGESF(phoEt, phoEta) * objectESF.getElectronTRGESF(lepPt, lepEta);
        } else if (channelType == 2) {
            scalefactor = objectESF.getMuonESF(lepPt, lepEta) * objectESF.getPhotonESF(phoEt, phoEta)
                        * objectESF.getMuonEGTRGESF(phoEt, lepPt);
        }

        float XS_weight = 1;
        if (channelType == 1) {
            if      (RunYear == 2016 && preVFP)  XS_weight = lumi_2016preVFP_DoubleEG  * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2016 && !preVFP) XS_weight = lumi_2016postVFP_DoubleEG * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2017)            XS_weight = lumi_2017_DoubleEG        * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2018)            XS_weight = lumi_2018_DoubleEG        * 1000 * crosssection / ntotalevent;
        } else if (channelType == 2) {
            if      (RunYear == 2016 && preVFP)  XS_weight = lumi_2016preVFP_MuonEG  * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2016 && !preVFP) XS_weight = lumi_2016postVFP_MuonEG * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2017)            XS_weight = lumi_2017_MuonEG        * 1000 * crosssection / ntotalevent;
            else if (RunYear == 2018)            XS_weight = lumi_2018_MuonEG        * 1000 * crosssection / ntotalevent;
        }

        float weight = PUweight * XS_weight * scalefactor;

        double mindRpho(0.3);
        int phoIndex(0), anyphoIndex(-1);
        for (unsigned iMC(0); iMC < mcPID->size(); iMC++) {
            double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
            if (dR < mindRpho) { mindRpho = dR; phoIndex = iMC; }
            if (dR < 0.3 && fabs((*mcPID)[iMC]) == 22) anyphoIndex = iMC;
        }
        bool isTruePho(false);
        auto passMom = [](int momPID){
            return fabs(momPID) <= 6 || fabs(momPID) == 21 || fabs(momPID) == 999
                || fabs(momPID) == 11 || fabs(momPID) == 13 || fabs(momPID) == 15
                || fabs(momPID) == 23 || fabs(momPID) == 24;
        };
        if (mindRpho < 0.1) {
            if ((*mcPID)[phoIndex] == 22 && passMom((*mcMomPID)[phoIndex])) isTruePho = true;
        } else if (anyphoIndex >= 0) {
            if ((*mcPID)[anyphoIndex] == 22 && passMom((*mcMomPID)[anyphoIndex])) isTruePho = true;
        }
        if (!isTruePho) continue;

        b.hMET->Fill(sigMET, weight);
        b.hMT->Fill(sigMT, weight);
        b.h2D->Fill(sigMET, sigMT, weight);
    }

    delete tree;
    return b;
}

void plot_SignalContam_vs_MET(int channelType = 1, int RunYear = 2018, bool preVFP = false)
{
    gROOT->SetBatch(kTRUE);
    gSystem->Load("../../lib/libAnaClasses.so");
    setTDRStyle();
    gStyle->SetOptStat(0);
    esfScaleFactor objectESF;

    std::string whichVFP;
    if (RunYear == 2016 && preVFP)  whichVFP = "preVFP";
    if (RunYear == 2016 && !preVFP) whichVFP = "postVFP";
    if (RunYear == 2017 || RunYear == 2018) whichVFP = "";

    const char *treeName = (channelType == 1) ? "egTree" : "mgTree";

    const TString baseDir = "/eos/uscms/store/user/tmishra/VGamma/";
    std::vector<SampleDef> samples = {
        {"TTJets",     "t#bar{t}+jets",        "#3f90da"},
        {"TTGJets",    "t#bar{t}#gamma+jets",  "#ffa90e"},
        {"WWG",        "WW#gamma",             "#bd1f01"},
        {"WZG",        "WZ#gamma",             "#94a4a2"},
        {"WW",         "WW",                   "#832db6"},
        {"WZ",         "WZ",                   "#a96b59"},
        {"WGToLNuG",   "W#gamma (l#nu)",       "#e76300"},
        {"WGJet40",    "W+jets (#gamma p_{T}>40)",  "#b9ac70"},
        {"WGJet130",   "W+jets (#gamma p_{T}>130)", "#717581"},
        {"DYJetsToLL", "DY+jets#rightarrowll",  "#92dadd"},
        {"ZGToLLG",    "Z#gamma#rightarrowll#gamma", "#9c9ca1"},
	{"QCD",        "QCD",                       "#e377c2"},
	{"WJetsToLNu", "W+Jets (l#nu)",             "#3288bd"}
    };

    const TString signalBaseDir = "/uscms/home/tmishra/nobackup/signal_trees/";
    // MET comes from the premade 1D histograms (confirmed working source).
    std::vector<SignalDef> signals = {
        {"TChiWg(800)",      Form("%sTree_TChiWG_%d%s.root", signalBaseDir.Data(), RunYear, whichVFP.c_str()), "p_tchiwg_MET_signal",         kRed+1,   1},
        {"T5Wg(1700,1000)",  Form("%sTree_T5WG_%d%s.root",   signalBaseDir.Data(), RunYear, whichVFP.c_str()), "p_t5wg_MET_signal_1700_1000", kBlue+1,  2},
        {"T5Wg(1800,600)",   Form("%sTree_T5WG_%d%s.root",   signalBaseDir.Data(), RunYear, whichVFP.c_str()), "p_t5wg_MET_signal_1800_600",  kGreen+2, 4},
        {"T5Wg(1900,100)",   Form("%sTree_T5WG_%d%s.root",   signalBaseDir.Data(), RunYear, whichVFP.c_str()), "p_t5wg_MET_signal_1900_100",  kMagenta+1, 5}
    };
    // EDIT ME: raw signal ntuple (egTree/mgTree with sigMET & sigMT branches),
    // same file layout as the background samples above. Needed for the MT
    // and 2D (MET,MT) plots below -- the premade files above only hold a
    // 1D MET histogram, so they can't supply MT for signal. Leave this
    // pointing at a real file (or set hasSignalTree=false to skip signal
    // on those two plots) once you tell me where it lives.
    bool hasSignalTree = false;
    const TString signalTreeFile = "";  

    // ------------------------------------------------------------------
    // Build background histogram bundles (MET, MT, 2D) and the MET stack
    // (smallest integral at bottom)
    // ------------------------------------------------------------------
    std::vector<std::pair<HistBundle, SampleDef>> bkgHists;
    for (auto &s : samples) {
        TString file = Form("%sresTree_VGamma_%s_%d%s.root", baseDir.Data(), s.procName.Data(), RunYear, whichVFP.c_str());
        HistBundle b = MakeBkgHistBundle(s.procName.Data(), file, treeName, channelType, RunYear, preVFP, objectESF);
        b.hMET->SetFillColor(TColor::GetColor(s.colorHex.Data())); b.hMET->SetLineColor(TColor::GetColor(s.colorHex.Data())); b.hMET->SetFillStyle(1001);
        b.hMT->SetFillColor(TColor::GetColor(s.colorHex.Data()));  b.hMT->SetLineColor(TColor::GetColor(s.colorHex.Data()));  b.hMT->SetFillStyle(1001);
        bkgHists.push_back({b, s});
    }
    std::sort(bkgHists.begin(), bkgHists.end(), [](const auto&a, const auto&b){
        return a.first.hMET->Integral() < b.first.hMET->Integral();
    });

    THStack *stack = new THStack("metStack", "");
    for (auto &hp : bkgHists) stack->Add(hp.first.hMET);

    // Summed total-background MT and 2D histograms (used for the new plots)
    TH1D *totalBkgMT = (TH1D*)bkgHists.front().first.hMT->Clone("totalBkgMT");
    TH2D *totalBkg2D = (TH2D*)bkgHists.front().first.h2D->Clone("totalBkg2D");
    for (size_t i = 1; i < bkgHists.size(); i++) {
        totalBkgMT->Add(bkgHists[i].first.hMT);
        totalBkg2D->Add(bkgHists[i].first.h2D);
    }

    // ------------------------------------------------------------------
    // Build signal MET histograms (unstacked line overlays) -- these are
    // pre-made TH1D's already sitting in the signal files, not trees,
    // so read them with TFile::Get() rather than a tree loop.
    // ------------------------------------------------------------------
    const char *chanSuffix = (channelType == 1) ? "eg" : "mg";
    std::vector<std::pair<TH1D*, SignalDef>> sigHists;
    for (auto &sig : signals) {
        TFile *fsig = TFile::Open(sig.filepath);
        if (!fsig || fsig->IsZombie()) {
            std::cout << "[WARNING] could not open " << sig.filepath << " -- skipping " << sig.label << std::endl;
            continue;
        }
        TString keyname = Form("%s_%s", sig.histname.Data(), chanSuffix);
        TH1D *hsrc = (TH1D*)fsig->Get(keyname);
        if (!hsrc) {
            std::cout << "[WARNING] histogram " << keyname << " not found in " << sig.filepath
                       << " -- skipping " << sig.label << std::endl;
            fsig->Close();
            continue;
        }
        TH1D *h = (TH1D*)hsrc->Clone(Form("p_MET_sig_%s", sig.label.Data()));
        h->SetDirectory(0);
        fsig->Close();

        h->SetLineColor(sig.lineColor);
        h->SetLineStyle(sig.lineStyle);
        h->SetLineWidth(3);
        h->SetFillStyle(0);
        sigHists.push_back({h, sig});
    }

    // Signal MT histograms -- premade HT histograms, same file/Get() pattern
    // as the MET ones above (key name: MET -> HT substitution).
    std::vector<std::pair<TH1D*, SignalDef>> sigMTHists;
    for (auto &sig : signals) {
        TFile *fsig = TFile::Open(sig.filepath);
        if (!fsig || fsig->IsZombie()) { std::cout << "[WARNING] could not open " << sig.filepath << " for MT -- skipping " << sig.label << std::endl; continue; }
        TString histnameMT = sig.histname; histnameMT.ReplaceAll("MET", "Mt");   // was "HT" (bug: signal HT drawn as MT), fixed 2026-11-02
        TString keynameMT = Form("%s_%s", histnameMT.Data(), chanSuffix);
        TH1D *hsrc = (TH1D*)fsig->Get(keynameMT);
        if (!hsrc) {
            std::cout << "[WARNING] histogram " << keynameMT << " not found in " << sig.filepath << " -- skipping " << sig.label << " for MT" << std::endl;
            fsig->Close();
            continue;
        }
        TH1D *h = (TH1D*)hsrc->Clone(Form("p_MT_sig_%s", sig.label.Data()));
        h->SetDirectory(0);
        fsig->Close();
        h->SetLineColor(sig.lineColor); h->SetLineStyle(sig.lineStyle); h->SetLineWidth(3); h->SetFillStyle(0);
        sigMTHists.push_back({h, sig});
    }

    // Signal 2D (MET,MT) histograms -- still needs the raw signal tree
    // (a 1D MET histogram + a 1D MT/HT histogram cannot supply their
    // correlation), only built if a raw signal tree is configured.
    std::vector<std::pair<TH2D*, SignalDef>> sig2DHists;
    if (hasSignalTree) {
        for (auto &sig : signals) {
            HistBundle sb = MakeBkgHistBundle(("sig_"+sig.label).Data(), signalTreeFile, treeName,
                                               channelType, RunYear, preVFP, objectESF);
            // NOTE: reuses MakeBkgHistBundle purely for its tree-reading loop;
            // signal doesn't need the isTruePho requirement, but leaving it in
            // is harmless since signal photons are genuine prompt photons.
            sig2DHists.push_back({sb.h2D, sig});
        }
    } else {
        std::cout << "[INFO] hasSignalTree=false -- the 2D plot below will show BACKGROUND ONLY for signal. "
                      "Set signalTreeFile to a raw egTree/mgTree signal ntuple to add signal there too." << std::endl;
    }

    // ------------------------------------------------------------------
    // Draw: single log-y pad, stacked bkg + signal lines, CR/SR lines
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("SignalContam_MET", "Signal contamination vs MET", 700, 650);
    c->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.13);

    stack->Draw("hist");   // must draw once before axes/max are available
    stack->SetMaximum(std::max(stack->GetMaximum(), 1.0) * 500);
    stack->SetMinimum(1e-2);
    stack->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    stack->GetYaxis()->SetTitle("Events / bin");
    stack->GetXaxis()->SetRangeUser(0, 600);
    stack->Draw("hist");

    for (auto &hp : sigHists) hp.first->Draw("hist same");

    // CR / SR boundary lines: MET = 70 (CR upper edge) and MET = 120 (SR lower edge)
    double ymin = stack->GetMinimum();
    double ymax = stack->GetMaximum();
    TLine *lineCR = new TLine(70, 0, 70, 1e7);
    lineCR->SetLineStyle(2); lineCR->SetLineWidth(2); lineCR->SetLineColor(kGray+2);
    lineCR->Draw("same");
    TLine *lineSR = new TLine(120, 0, 120, 1e7);
    lineSR->SetLineStyle(2); lineSR->SetLineWidth(2); lineSR->SetLineColor(kGray+2);
    lineSR->Draw("same");

    TLatex tex;
    tex.SetNDC(); tex.SetTextSize(0.025);
    tex.DrawLatex(0.16, 0.87, "MET<70");
    tex.DrawLatex(0.27, 0.87, "MET: 70-120");
    tex.DrawLatex(0.50, 0.87, "SR: MET>120");

    TLegend *leg = new TLegend(0.39, 0.6, 0.93, 0.84);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);

    for (auto &hp : bkgHists) {
    if (hp.first.hMET) {
        hp.first.hMET->SetMarkerSize(0);
        hp.first.hMET->SetLineWidth(0);
        leg->AddEntry(hp.first.hMET, hp.second.label, "f");
    }
    }

    for (auto &hp : sigHists) {
    if (hp.first) {
        hp.first->SetMarkerSize(0);
        leg->AddEntry(hp.first, hp.second.label, "l");
    }
    }
    leg->Draw("same");

    //if      (RunYear == 2016 && preVFP)  CMS_lumi(gPad, 1, channelType, 11);
    //else if (RunYear == 2016 && !preVFP) CMS_lumi(gPad, 2, channelType, 11);
    //else if (RunYear == 2017)            CMS_lumi(gPad, 3, channelType, 11);
    //else if (RunYear == 2018)            CMS_lumi(gPad, 4, channelType, 11);
    {   // "CMS Simulation" + period luminosity (simulation only, no data)
        double L = 0;
        if      (RunYear == 2016 && preVFP)  L = (channelType == 1) ? lumi_2016preVFP_DoubleEG  : lumi_2016preVFP_MuonEG;
        else if (RunYear == 2016 && !preVFP) L = (channelType == 1) ? lumi_2016postVFP_DoubleEG : lumi_2016postVFP_MuonEG;
        else if (RunYear == 2017)            L = (channelType == 1) ? lumi_2017_DoubleEG        : lumi_2017_MuonEG;
        else if (RunYear == 2018)            L = (channelType == 1) ? lumi_2018_DoubleEG        : lumi_2018_MuonEG;
        TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.04);
        lt.DrawLatex(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin() + 0.01, "#bf{CMS} #it{Simulation}");
        lt.SetTextAlign(31);
        lt.DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("%.1f fb^{-1} (13 TeV)", L));
    }
    gPad->RedrawAxis();

    TString outname = Form("SignalContam_MET_%s_%d%s.pdf",
                            (channelType == 1) ? "egamma" : "mgamma", RunYear, whichVFP.c_str());
    c->SaveAs(outname);
    c->SaveAs(TString(outname).ReplaceAll(".pdf", ".png"));

    // ------------------------------------------------------------------
    // Quantitative contamination table -- the number the reviewer wants:
    // what fraction of each signal's yield falls in MET<70 / 70-120 / >120
    // ------------------------------------------------------------------
    std::cout << "\n=== Signal fraction by MET region (this is the number to quote to the reviewer) ===\n";
    std::cout << Form("%-35s %12s %12s %12s\n", "Signal point", "MET<70(%)", "70-120(%)", "MET>120(%)");
    for (auto &hp : sigHists) {
        TH1D *h = hp.first;
        double total = h->Integral(0, h->GetNbinsX() + 1);
        if (total <= 0) { std::cout << "  [skip: zero yield] " << hp.second.label << std::endl; continue; }
        double lowCR  = h->Integral(h->FindBin(0),   h->FindBin(70 - 1e-3));
        double midBand= h->Integral(h->FindBin(70),  h->FindBin(120 - 1e-3));
        double highSR = h->Integral(h->FindBin(120), h->GetNbinsX() + 1);
        std::cout << Form("%-35s %12.2f %12.2f %12.2f\n", hp.second.label.Data(),
                           100 * lowCR / total, 100 * midBand / total, 100 * highSR / total);
    }
    std::cout << "=====================================================================\n" << std::endl;

    // ==================================================================
    // 1. MT plot
    //    (VR: MT<100  vs  MT>100)
    // ==================================================================
    THStack *stackMT = new THStack("mtStack", "");
    for (auto &hp : bkgHists) stackMT->Add(hp.first.hMT);

    TCanvas *cMT = new TCanvas("SignalContam_MT", "Signal contamination vs MT", 700, 650);
    cMT->cd();
    gPad->SetLogy();
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.13);

    stackMT->Draw("hist");
    stackMT->SetMaximum(std::max(stackMT->GetMaximum(), 1.0) * 500);
    stackMT->SetMinimum(1e-2);
    stackMT->GetXaxis()->SetTitle("M_{T} (GeV)");
    stackMT->GetYaxis()->SetTitle("Events / bin");
    stackMT->GetXaxis()->SetRangeUser(0, 1000);
    stackMT->Draw("hist");

    for (auto &hp : sigMTHists) hp.first->Draw("hist same");

    TLine *lineVR = new TLine(100, 0, 100, 1e7);
    lineVR->SetLineStyle(2); lineVR->SetLineWidth(2); lineVR->SetLineColor(kGray+2);
    lineVR->Draw("same");

    TLatex texMT;
    texMT.SetNDC(); texMT.SetTextSize(0.025);
    texMT.DrawLatex(0.16, 0.87, "VR: M_{T}<100");
    texMT.DrawLatex(0.55, 0.87, "M_{T}>100");

    TLegend *legMT = new TLegend(0.39, 0.6, 0.93, 0.84);
    legMT->SetNColumns(2);
    legMT->SetBorderSize(0);
    legMT->SetFillStyle(0);
    legMT->SetFillColor(0);
    for (auto &hp : bkgHists) {
        if (hp.first.hMT) legMT->AddEntry(hp.first.hMT, hp.second.label, "f");
    }
    for (auto &hp : sigMTHists) {
        if (hp.first) legMT->AddEntry(hp.first, hp.second.label, "l");
    }
    legMT->Draw("same");
    {   // "CMS Simulation" + period luminosity (simulation only, no data)
        double L = 0;
        if      (RunYear == 2016 && preVFP)  L = (channelType == 1) ? lumi_2016preVFP_DoubleEG  : lumi_2016preVFP_MuonEG;
        else if (RunYear == 2016 && !preVFP) L = (channelType == 1) ? lumi_2016postVFP_DoubleEG : lumi_2016postVFP_MuonEG;
        else if (RunYear == 2017)            L = (channelType == 1) ? lumi_2017_DoubleEG        : lumi_2017_MuonEG;
        else if (RunYear == 2018)            L = (channelType == 1) ? lumi_2018_DoubleEG        : lumi_2018_MuonEG;
        TLatex lt; lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(0.04);
        lt.DrawLatex(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin() + 0.01, "#bf{CMS} #it{Simulation}");
        lt.SetTextAlign(31);
        lt.DrawLatex(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin() + 0.01, Form("%.1f fb^{-1} (13 TeV)", L));
    }
    gPad->RedrawAxis();

    TString outnameMT = Form("SignalContam_MT_%s_%d%s.pdf",
                              (channelType == 1) ? "egamma" : "mgamma", RunYear, whichVFP.c_str());
    cMT->SaveAs(outnameMT);
    cMT->SaveAs(TString(outnameMT).ReplaceAll(".pdf", ".png"));
    // signal / background fractions passing the SR cuts (for the ARC reply)
    for (auto &hp : sigHists)   { TH1D *h = hp.first; std::cout << "[eff] " << hp.second.label << "  MET>120: " << h->Integral(h->FindBin(120.01), h->GetNbinsX()+1) / h->Integral(0, h->GetNbinsX()+1) << std::endl; }
    for (auto &hp : sigMTHists) { TH1D *h = hp.first; std::cout << "[eff] " << hp.second.label << "  MT>100: "  << h->Integral(h->FindBin(100.01), h->GetNbinsX()+1) / h->Integral(0, h->GetNbinsX()+1) << std::endl; }
    { TH1D *h = (TH1D*)stack->GetStack()->Last();   std::cout << "[eff] SM bkg (true-photon MC)  MET>120: " << h->Integral(h->FindBin(120.01), h->GetNbinsX()+1) / h->Integral(0, h->GetNbinsX()+1) << std::endl; }
    { TH1D *h = totalBkgMT;                          std::cout << "[eff] SM bkg (true-photon MC)  MT>100: "  << h->Integral(h->FindBin(100.01), h->GetNbinsX()+1) / h->Integral(0, h->GetNbinsX()+1) << std::endl; }

    // MT fraction table: background (per process + total) and signal (if available)
    std::cout << "\n=== MT<100 (VR) vs MT>100 fraction ===\n";
    std::cout << Form("%-35s %12s %12s\n", "Sample", "MT<100(%)", "MT>100(%)");
    for (auto &hp : bkgHists) {
        TH1D *h = hp.first.hMT;
        double total = h->Integral(0, h->GetNbinsX() + 1);
        if (total <= 0) continue;
        double lowVR  = h->Integral(h->FindBin(0),   h->FindBin(100 - 1e-3));
        double highSR = h->Integral(h->FindBin(100), h->GetNbinsX() + 1);
        std::cout << Form("%-35s %12.2f %12.2f\n", ("[bkg] "+hp.second.label).Data(),
                           100 * lowVR / total, 100 * highSR / total);
    }
    {
        double total = totalBkgMT->Integral(0, totalBkgMT->GetNbinsX() + 1);
        double lowVR  = totalBkgMT->Integral(totalBkgMT->FindBin(0),   totalBkgMT->FindBin(100 - 1e-3));
        double highSR = totalBkgMT->Integral(totalBkgMT->FindBin(100), totalBkgMT->GetNbinsX() + 1);
        std::cout << Form("%-35s %12.2f %12.2f\n", "[bkg] TOTAL", 100 * lowVR / total, 100 * highSR / total);
    }
    for (auto &hp : sigMTHists) {
        TH1D *h = hp.first;
        double total = h->Integral(0, h->GetNbinsX() + 1);
        if (total <= 0) { std::cout << "  [skip: zero yield] " << hp.second.label << std::endl; continue; }
        double lowVR  = h->Integral(h->FindBin(0),   h->FindBin(100 - 1e-3));
        double highSR = h->Integral(h->FindBin(100), h->GetNbinsX() + 1);
        std::cout << Form("%-35s %12.2f %12.2f\n", ("[sig] "+hp.second.label).Data(),
                           100 * lowVR / total, 100 * highSR / total);
    }
    if (sigMTHists.empty()) std::cout << "  (no signal MT rows -- check that the p_t5wg_HT_signal_*/p_tchiwg_HT_signal_* keys exist)\n";
    std::cout << "=====================================================================\n" << std::endl;

    // ==================================================================
    // 2. 2D plot: MET (x) vs MT (y), total background as colz map,
    //    signal (if available) overlaid as contour lines, CR/VR/SR
    //    regions marked exactly as defined at the top of the analysis:
    //      CR: MET<70 (any MT)   VR: MT<100 (any MET)   SR: MT>100 & MET>120
    // ==================================================================
    TCanvas *c2D = new TCanvas("SignalContam_MET_vs_MT", "MET vs MT", 700, 650);
    c2D->cd();
    gPad->SetLogz();
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.13);

    totalBkg2D->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    totalBkg2D->GetYaxis()->SetTitle("M_{T} (GeV)");
    totalBkg2D->GetXaxis()->SetRangeUser(0, 400);
    totalBkg2D->GetYaxis()->SetRangeUser(0, 400);
    totalBkg2D->SetTitle("Total background");
    totalBkg2D->Draw("colz");

    for (auto &hp : sig2DHists) {
        hp.first->SetLineColor(hp.second.lineColor);
        hp.first->SetLineWidth(2);
        hp.first->Draw("cont3 same");
    }

    double xmaxPlot = 400, ymaxPlot = 400;
    TLine *l2_CR = new TLine(70, 0, 70, ymaxPlot);
    l2_CR->SetLineStyle(2); l2_CR->SetLineWidth(2); l2_CR->SetLineColor(kBlack); l2_CR->Draw("same");
    TLine *l2_VR = new TLine(0, 100, xmaxPlot, 100);
    l2_VR->SetLineStyle(2); l2_VR->SetLineWidth(2); l2_VR->SetLineColor(kBlack); l2_VR->Draw("same");
    TBox *box_SR = new TBox(120, 100, xmaxPlot, ymaxPlot);
    box_SR->SetFillStyle(0); box_SR->SetLineStyle(1); box_SR->SetLineWidth(3); box_SR->SetLineColor(kRed+1);
    box_SR->Draw("same");

    TLatex tex2D;
    tex2D.SetNDC(); tex2D.SetTextSize(0.025); tex2D.SetTextColor(kBlack);
    tex2D.DrawLatex(0.16, 0.93, "CR: MET<70");
    tex2D.DrawLatex(0.55, 0.30, "VR: M_{T}<100");
    tex2D.SetTextColor(kRed+1);
    tex2D.DrawLatex(0.65, 0.80, "SR");

    TLegend *leg2D = new TLegend(0.16, 0.65, 0.5, 0.85);
    leg2D->SetBorderSize(0);
    leg2D->SetFillStyle(0);
    for (auto &hp : sig2DHists) leg2D->AddEntry(hp.first, hp.second.label, "l");
    if (!sig2DHists.empty()) leg2D->Draw("same");
    gPad->RedrawAxis();

    TString outname2D = Form("SignalContam_MET_vs_MT_%s_%d%s.pdf",
                              (channelType == 1) ? "egamma" : "mgamma", RunYear, whichVFP.c_str());
    c2D->SaveAs(outname2D);

    // 2D CR/VR/SR fraction table -- background (per process + total) and signal
    auto region2DFractions = [](TH2D *h, double &crPct, double &vrPct, double &srPct) {
        double total = h->Integral(0, h->GetNbinsX() + 1, 0, h->GetNbinsY() + 1);
        if (total <= 0) { crPct = vrPct = srPct = 0; return false; }
        double cr = h->Integral(h->GetXaxis()->FindBin(0.0),   h->GetXaxis()->FindBin(70 - 1e-3),
                                 0, h->GetNbinsY() + 1);
        double vr = h->Integral(0, h->GetNbinsX() + 1,
                                 h->GetYaxis()->FindBin(0.0),   h->GetYaxis()->FindBin(100 - 1e-3));
        double sr = h->Integral(h->GetXaxis()->FindBin(120), h->GetNbinsX() + 1,
                                 h->GetYaxis()->FindBin(100), h->GetNbinsY() + 1);
        crPct = 100 * cr / total; vrPct = 100 * vr / total; srPct = 100 * sr / total;
        return true;
    };

    std::cout << "\n=== 2D region fraction: CR(MET<70) / VR(MT<100) / SR(MT>100 & MET>120) ===\n";
    std::cout << Form("%-35s %10s %10s %10s\n", "Sample", "CR(%)", "VR(%)", "SR(%)");
    for (auto &hp : bkgHists) {
        double cr, vr, sr;
        if (region2DFractions(hp.first.h2D, cr, vr, sr))
            std::cout << Form("%-35s %10.2f %10.2f %10.2f\n", ("[bkg] "+hp.second.label).Data(), cr, vr, sr);
    }
    {
        double cr, vr, sr;
        region2DFractions(totalBkg2D, cr, vr, sr);
        std::cout << Form("%-35s %10.2f %10.2f %10.2f\n", "[bkg] TOTAL", cr, vr, sr);
    }
    for (auto &hp : sig2DHists) {
        double cr, vr, sr;
        if (region2DFractions(hp.first, cr, vr, sr))
            std::cout << Form("%-35s %10.2f %10.2f %10.2f\n", ("[sig] "+hp.second.label).Data(), cr, vr, sr);
        else
            std::cout << "  [skip: zero yield] " << hp.second.label << std::endl;
    }
    if (!hasSignalTree) std::cout << "  (signal rows above are empty -- set signalTreeFile/hasSignalTree to include them)\n";
    std::cout << "=====================================================================\n" << std::endl;
}
