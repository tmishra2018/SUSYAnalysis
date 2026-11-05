// plotHEMEffect.C
//
// Check the effect of the HEM15/16 veto (passHEMVeto branch) on 2018
// DoubleEG and MuonEG data, across all four trees and four kinematic
// variables. For each of the 2 x 4 x 4 = 32 combinations, draws MET/HT/
// Photon pT/MT with no requirement on passHEMVeto ("without HEM veto") vs
// with passHEMVeto==1 required ("with HEM veto"), in the two-pad
// log-scale/ratio style from the reference plot.
//
// Run with:  root -l -b -q plotHEMEffect.C
//
// NOTE: the channel -> input file mapping below (DoubleEG -> *_egsignal_*,
// MuonEG -> *_mgsignal_*) is inferred from the file listing you gave.
// Double check it matches your actual egsignal/mgsignal <-> DoubleEG/MuonEG
// convention before trusting the output.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <utility>

// --- Binning, exactly as given ---
Double_t bkgEtBins[]  = {35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
Double_t bkgMETBins[] = {0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
Double_t bkgMtBins[]  = {0,20,40,60,80,100,120,140,160,180,200,300,400,500,600,800,1000};
Double_t bkgHTBins[]  = {0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};

struct VarDef {
    TString branch;   // TTree branch name
    TString label;    // x-axis label
    Double_t* bins;   // variable bin edges
    Int_t nBins;       // number of bins = len(bins)-1
};

void MakeHEMCheckPlot(TTree* tree, const VarDef& v,
                       const TString& channel, const TString& treeName,
                       const TString& outdir){

    if(!tree){
        std::cerr << "WARNING: null tree for " << channel << " " << treeName
                   << " " << v.branch << " -- skipping." << std::endl;
        return;
    }

    TString tag         = channel + "_" + treeName + "_" + v.branch;
    TString hNoVetoName = "h_noveto_" + tag;
    TString hVetoName   = "h_veto_"   + tag;

    TH1D* hNoVeto = new TH1D(hNoVetoName, "", v.nBins, v.bins);
    TH1D* hVeto   = new TH1D(hVetoName,   "", v.nBins, v.bins);
    hNoVeto->Sumw2();
    hVeto->Sumw2();

    // "without HEM veto": all events, no requirement on passHEMVeto.
    // "with HEM veto": passHEMVeto == 1 required.
    // Both unweighted (raw event counts), as requested.
    tree->Draw((v.branch + ">>" + hNoVetoName).Data(), "", "goff");
    tree->Draw((v.branch + ">>" + hVetoName).Data(), "passHEMVeto==1", "goff");

    if(hNoVeto->GetEntries() == 0){
        std::cerr << "WARNING: zero entries (without-veto) for " << tag
                   << " -- check branch/tree/selection." << std::endl;
        delete hNoVeto; delete hVeto;
        return;
    }

    Double_t nNoVeto = hNoVeto->Integral();
    Double_t nVeto   = hVeto->Integral();

    TCanvas* c = new TCanvas("c_" + tag, "", 580, 620);
    c->cd();

    TPad* padTop = new TPad("padTop_" + tag, "", 0, 0.30, 1, 1.0);
    TPad* padBot = new TPad("padBot_" + tag, "", 0, 0.0,  1, 0.30);
    padTop->SetBottomMargin(0.02);
    padBot->SetTopMargin(0.03);
    padBot->SetBottomMargin(0.35);
    padTop->Draw();
    padBot->Draw();

    // --- Top pad: log-scale overlay ---
    padTop->cd();
    padTop->SetLogy();

    hNoVeto->SetMarkerStyle(20);
    hNoVeto->SetMarkerColor(kBlack);
    hNoVeto->SetLineColor(kBlack);
    hNoVeto->SetTitle("");
    hNoVeto->GetYaxis()->SetTitle("Events");
    hNoVeto->GetYaxis()->SetTitleOffset(0.75);
    hNoVeto->GetYaxis()->SetLabelSize(0.05);
    hNoVeto->GetYaxis()->SetTitleSize(0.06);
    hNoVeto->GetXaxis()->SetLabelSize(0);
    Double_t ymax = hNoVeto->GetMaximum();
    hNoVeto->GetYaxis()->SetRangeUser(std::max(0.5, hNoVeto->GetMinimum(0.0)*0.3), ymax*20);
    hNoVeto->Draw("PE");

    hVeto->SetLineColor(kBlue+1);
    hVeto->SetLineWidth(2);
    hVeto->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.40, 0.70, 0.88, 0.86);
    leg->SetFillStyle(0);
    leg->SetLineColor(kWhite);
    leg->SetTextSize(0.042);
    leg->AddEntry(hNoVeto, Form("No veto [N=%.0f]", nNoVeto), "pe");
    leg->AddEntry(hVeto,   Form("HEM veto [N=%.0f]",   nVeto),   "l");
    leg->Draw();

    TLatex la;
    la.SetNDC();
    la.SetTextFont(62);
    la.SetTextSize(0.05);
    la.DrawLatex(0.12, 0.91, "CMS");
    la.SetTextFont(42);
    la.DrawLatex(0.22, 0.91, "Preliminary");
    la.SetTextAlign(31);
    la.DrawLatex(0.90, 0.91, channel);

    // --- Bottom pad: ratio ---
    // Use binomial errors ("B"): the veto sample is a SUBSET of the
    // no-veto sample (same underlying events, veto just removes some), so
    // treating the two histograms' errors as independent (default Divide
    // behavior) would overestimate the ratio uncertainty. Binomial
    // (efficiency-style) errors are the statistically correct choice here.
    padBot->cd();
    TH1D* hRatio = (TH1D*)hNoVeto->Clone("hRatio_" + tag);
    hRatio->Reset();
    hRatio->Divide(hVeto, hNoVeto, 1.0, 1.0, "B");
    hRatio->SetTitle("");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.1);
    hRatio->GetYaxis()->SetTitle("veto / no veto");
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetTitleOffset(0.35);
    hRatio->GetYaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetXaxis()->SetTitle(v.label);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetXaxis()->SetTitleOffset(1.05);
    hRatio->Draw("PE");

    TLine* line1 = new TLine(v.bins[0], 1.0, v.bins[v.nBins], 1.0);
    line1->SetLineColor(kRed);
    line1->SetLineStyle(2);
    line1->Draw();

    gSystem->mkdir(outdir, kTRUE);
    TString outname = outdir + "/" + tag + "_HEMcheck";
    c->SaveAs(outname + ".pdf");
    c->SaveAs(outname + ".png");

    std::cout << "Wrote " << outname << ".pdf  (N_without=" << nNoVeto
               << ", N_with=" << nVeto << ")" << std::endl;

    delete c;
}

void plotHEMEffect(){

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // channel -> input file. See NOTE at top of file about the
    // egsignal/mgsignal <-> DoubleEG/MuonEG mapping.
    std::vector<std::pair<TString,TString>> channels = {
        {"Ele-Gamma", "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"},
        {"Mu-Gamma",   "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root"}
    };

    std::vector<TString> treeNames = {"signalTree", "proxyTree", "jetTree", "fakeLepTree"};

    // Compute nBins from the arrays themselves instead of hardcoding counts.
    const Int_t nEtBins  = sizeof(bkgEtBins)  / sizeof(bkgEtBins[0])  - 1;
    const Int_t nMETBins = sizeof(bkgMETBins) / sizeof(bkgMETBins[0]) - 1;
    const Int_t nMtBins  = sizeof(bkgMtBins)  / sizeof(bkgMtBins[0])  - 1;
    const Int_t nHTBins  = sizeof(bkgHTBins)  / sizeof(bkgHTBins[0])  - 1;

    std::vector<VarDef> vars = {
        {"phoEt",  "Photon p_{T} (GeV)", bkgEtBins,  nEtBins},
        {"sigMET", "MET (GeV)",          bkgMETBins, nMETBins},
        {"sigMT",  "M_{T} (GeV)",        bkgMtBins,  nMtBins},
        {"HT",     "H_{T} (GeV)",        bkgHTBins,  nHTBins}
    };

    TString outdir = "HEMCheck_plots";
    Int_t nMade = 0, nSkipped = 0;

    for(auto& ch : channels){
        TString channel  = ch.first;
        TString filepath = ch.second;

        TFile* f = TFile::Open(filepath);
        if(!f || f->IsZombie()){
            std::cerr << "WARNING: could not open " << filepath
                       << " -- skipping channel " << channel << std::endl;
            nSkipped += 4 * (Int_t)vars.size();
            continue;
        }

        for(auto& treeName : treeNames){
            TTree* tree = (TTree*)f->Get(treeName);
            if(!tree){
                std::cerr << "WARNING: " << treeName << " not found in " << filepath
                           << " -- skipping." << std::endl;
                nSkipped += (Int_t)vars.size();
                continue;
            }

            if(!tree->GetBranch("passHEMVeto")){
                std::cerr << "WARNING: passHEMVeto branch not found in " << channel
                           << " " << treeName
                           << " -- skipping (cannot make with/without-veto comparison)."
                           << std::endl;
                nSkipped += (Int_t)vars.size();
                continue;
            }

            for(auto& v : vars){
                if(!tree->GetBranch(v.branch)){
                    std::cerr << "WARNING: branch " << v.branch << " not found in "
                               << channel << " " << treeName << " -- skipping." << std::endl;
                    nSkipped++;
                    continue;
                }
                MakeHEMCheckPlot(tree, v, channel, treeName, outdir);
                nMade++;
            }
        }
        f->Close();
    }

    std::cout << "Done: " << nMade << " plot(s) written to " << outdir << "/, "
               << nSkipped << " skipped (see WARNINGs above)." << std::endl;
}
