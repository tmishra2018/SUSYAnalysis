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


void processEfficiency(const char* numeratorFilePath, const char* denominatorFilePath, const std::string& outputPrefix) {
    // Run ROOT in batch mode
    gROOT->SetBatch(kTRUE); // Disable GUI windows

    TFile* numeratorFile = TFile::Open(numeratorFilePath);
    TTree* numeratorTree = (TTree*)numeratorFile->Get("mgTree");

    TFile* denominatorFile = TFile::Open(denominatorFilePath);
    TTree* denominatorTree = (TTree*)denominatorFile->Get("mgTree");

    Float_t lepPt_numerator, lepEta_numerator;
    Float_t lepPt_denominator, lepEta_denominator;

    numeratorTree->SetBranchAddress("lepPt", &lepPt_numerator);
    numeratorTree->SetBranchAddress("lepEta", &lepEta_numerator);
    denominatorTree->SetBranchAddress("lepPt", &lepPt_denominator);
    denominatorTree->SetBranchAddress("lepEta", &lepEta_denominator);

    const int nBinsPt = 6;
    double ptBins[] = {10, 20, 35, 50, 90, 150, 500};
    Double_t etaBins[] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5};

    TH1F* hDenominator1D = new TH1F("hDenominator1D", ";Muon p_{T} [GeV];Events", nBinsPt, ptBins);
    TH1F* hNumerator1D = new TH1F("hNumerator1D", ";Muon p_{T} [GeV];Events", nBinsPt, ptBins);
    TH2F* hNumerator2D = new TH2F("hNumerator2D", ";p_{T} (GeV);#eta", 6, ptBins, 10, etaBins);
    TH2F* hDenominator2D = new TH2F("hDenominator2D", ";p_{T} (GeV);#eta", 6, ptBins, 10, etaBins);

    // Fill histograms for numerator
    Long64_t nEntriesNumerator = numeratorTree->GetEntries();
    for (Long64_t i = 0; i < nEntriesNumerator; ++i) {
        numeratorTree->GetEntry(i);
        hNumerator1D->Fill(lepPt_numerator);
        hNumerator2D->Fill(lepPt_numerator, lepEta_numerator);
    }

    // Fill histograms for denominator
    Long64_t nEntriesDenominator = denominatorTree->GetEntries();
    for (Long64_t i = 0; i < nEntriesDenominator; ++i) {
        denominatorTree->GetEntry(i);
        hDenominator1D->Fill(lepPt_denominator);
        hDenominator2D->Fill(lepPt_denominator, lepEta_denominator);
    }
	
    std::cout<< " ratio: "<< (1 - (static_cast<Double_t>(nEntriesNumerator) / static_cast<Double_t>(nEntriesDenominator)))*100.<<"\n";
    // Check bin contents before efficiency calculation
    for (int bin = 1; bin <= hNumerator1D->GetNbinsX(); ++bin) {
        double num = hNumerator1D->GetBinContent(bin);
        double den = hDenominator1D->GetBinContent(bin);
    }

    // Create efficiency object
    TEfficiency* efficiency1D = nullptr;
    if (TEfficiency::CheckConsistency(*hNumerator1D, *hDenominator1D)) {
        efficiency1D = new TEfficiency(*hNumerator1D, *hDenominator1D);
        efficiency1D->SetTitle("Muon Efficiency vs p_{T};Muon p_{T} [GeV];Efficiency");
    } else {
        std::cerr << "ERROR: TEfficiency histograms are inconsistent. Skipping efficiency calculation.\n";
    }

    // Draw efficiency plots
    TCanvas* c1 = new TCanvas("c1", "1D Efficiency Plot", 800, 600);
    c1->SetGrid();
    if (efficiency1D) efficiency1D->Draw("AP");
    c1->SaveAs(("efficiency_1D_muonPt_" + outputPrefix + ".png").c_str());

    TCanvas* c2 = new TCanvas("c2", "2D Efficiency Plot", 800, 600);
    gStyle->SetOptStat(0);
    hNumerator2D->Divide(hDenominator2D);
    hNumerator2D->GetZaxis()->SetLabelSize(0.02); 
    hNumerator2D->SetMarkerSize(0.8); 
    gStyle->SetPaintTextFormat(".2f"); 
    hNumerator2D->Draw("colz text");
    c2->SaveAs(("efficiency_2D_muonPtEta_" + outputPrefix + ".png").c_str());

    // Clean up
    delete hDenominator1D;
    delete hNumerator1D;
    delete hNumerator2D;
    delete hDenominator2D;
    delete efficiency1D;
    delete c1;
    delete c2;
    numeratorFile->Close();
    denominatorFile->Close();
    delete numeratorFile;
    delete denominatorFile;
}

void processEfficiency() {
    // Define input file paths and output prefixes
    const std::vector<std::pair<std::string, std::string>> inputFiles = {
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016postVFP_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016postVFP.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016preVFP_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016preVFP.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2017_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2017.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2018_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2018.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2016postVFP_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2016postVFP.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2016preVFP_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2016preVFP.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2017_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2017.root"},
        {"/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2018_muon_dxy_0p02_cut.root",
         "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2018.root"}
    };

    const std::vector<std::string> outputPrefixes = {
        "TChiWG_2016postVFP", "TChiWG_2016preVFP", "TChiWG_2017", "TChiWG_2018",
        "T5WG_2016postVFP", "T5WG_2016preVFP", "T5WG_2017", "T5WG_2018"
    };

    // Loop over all input files and process efficiency
    for (size_t i = 0; i < inputFiles.size(); ++i) {
        processEfficiency(inputFiles[i].first.c_str(), inputFiles[i].second.c_str(), outputPrefixes[i]);
    }
}
