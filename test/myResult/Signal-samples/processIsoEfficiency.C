#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <vector>

void processIsoEfficiency(const char* inputFilePath,
                           const std::string& leptonType,
                           const std::string& outputPrefix,
                           const std::string& isoBranchName = "lepMiniIso",
                           Float_t isoThreshold = 0.1) {

    gROOT->SetBatch(kTRUE); // Disable GUI windows
    TFile* inputFile = TFile::Open(inputFilePath);

    const std::string treeName = (leptonType == "muon") ? "mgTree" : "egTree";
    TTree* tree = (TTree*)inputFile->Get(treeName.c_str());

    Float_t lepPt, lepEta, isoVal;
    tree->SetBranchAddress("lepPt", &lepPt);
    tree->SetBranchAddress("lepEta", &lepEta);
    tree->SetBranchAddress(isoBranchName.c_str(), &isoVal);

    const int nBinsPt = 6;
    double ptBins[] = {10, 20, 35, 50, 90, 150, 500};
    Double_t etaBins[] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5};

    TH1F* hDenominator1D = new TH1F("hDenominator1D", (";" + leptonType + " p_{T} [GeV];Events").c_str(), nBinsPt, ptBins);
    TH1F* hNumerator1D   = new TH1F("hNumerator1D",   (";" + leptonType + " p_{T} [GeV];Events").c_str(), nBinsPt, ptBins);
    TH2F* hNumerator2D   = new TH2F("hNumerator2D",   ";p_{T} (GeV);#eta", 6, ptBins, 10, etaBins);
    TH2F* hDenominator2D = new TH2F("hDenominator2D", ";p_{T} (GeV);#eta", 6, ptBins, 10, etaBins);

    Long64_t nEntries = tree->GetEntries();
    Long64_t nPassIso = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Denominator: all reconstructed leptons (before isolation cut)
        hDenominator1D->Fill(lepPt);
        hDenominator2D->Fill(lepPt, lepEta);

        // Numerator: leptons passing the isolation requirement
        if (isoVal < isoThreshold) {
            hNumerator1D->Fill(lepPt);
            hNumerator2D->Fill(lepPt, lepEta);
            ++nPassIso;
        }
    }

    Double_t overallEfficiency = (nEntries > 0)
        ? 100. * static_cast<Double_t>(nPassIso) / static_cast<Double_t>(nEntries)
        : 0.;

    std::cout << "[" << outputPrefix << ", " << leptonType << "] "
              << "isolation efficiency: " << overallEfficiency << "% "
              << "(" << nPassIso << " / " << nEntries << ")\n";

    // Flag whenever isolation efficiency exceeds 99%
    if (overallEfficiency > 99.) {
        std::cout << "  >>> PASSES 99% threshold: " << outputPrefix
                   << " (" << leptonType << ") isolation efficiency = "
                   << overallEfficiency << "%\n";
    }

    // Create efficiency object (differential vs pT)
    TEfficiency* efficiency1D = nullptr;
    if (TEfficiency::CheckConsistency(*hNumerator1D, *hDenominator1D)) {
        efficiency1D = new TEfficiency(*hNumerator1D, *hDenominator1D);
        efficiency1D->SetTitle((leptonType + " Isolation Efficiency vs p_{T};" +
                                 leptonType + " p_{T} [GeV];Efficiency").c_str());
    } else {
        std::cerr << "ERROR: TEfficiency histograms are inconsistent. Skipping efficiency calculation.\n";
    }

    // Draw efficiency plots
    TCanvas* c1 = new TCanvas("c1", "1D Isolation Efficiency Plot", 800, 600);
    c1->SetGrid();
    if (efficiency1D) efficiency1D->Draw("AP");
    c1->SaveAs(("isoEfficiency_1D_" + leptonType + "Pt_" + outputPrefix + ".png").c_str());

    TCanvas* c2 = new TCanvas("c2", "2D Isolation Efficiency Plot", 800, 600);
    gStyle->SetOptStat(0);
    hNumerator2D->Divide(hDenominator2D);
    hNumerator2D->GetZaxis()->SetLabelSize(0.02);
    hNumerator2D->SetMarkerSize(0.8);
    gStyle->SetPaintTextFormat(".2f");
    hNumerator2D->Draw("colz text");
    c2->SaveAs(("isoEfficiency_2D_" + leptonType + "PtEta_" + outputPrefix + ".png").c_str());

    // Clean up
    delete hDenominator1D;
    delete hNumerator1D;
    delete hNumerator2D;
    delete hDenominator2D;
    delete efficiency1D;
    delete c1;
    delete c2;
    inputFile->Close();
    delete inputFile;
}

void processIsoEfficiency() {
    // Signal trees per era/model. Note: these are the baseline (uncut) resTree
    // files; the isolation cut is applied in-macro rather than via a pre-skim.
    const std::vector<std::pair<std::string, std::string>> inputFiles = {
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016postVFP.root", "TChiWG_2016postVFP"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016preVFP.root",  "TChiWG_2016preVFP"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2017.root",        "TChiWG_2017"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2018.root",        "TChiWG_2018"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2016postVFP.root",   "T5WG_2016postVFP"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2016preVFP.root",    "T5WG_2016preVFP"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2017.root",          "T5WG_2017"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2018.root",          "T5WG_2018"}
    };

    const std::string isoBranchNameMuon     = "MuMiniIso";
    const std::string isoBranchNameElectron = "EleMiniIso";
    const Float_t isoThreshold = 0.1;

    for (const auto& entry : inputFiles) {
        processIsoEfficiency(entry.first.c_str(), "muon",     entry.second, isoBranchNameMuon,     0.2);
        processIsoEfficiency(entry.first.c_str(), "electron", entry.second, isoBranchNameElectron, 0.1);
    }
}
