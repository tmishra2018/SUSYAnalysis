#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>

void plotSigMETSigMT() {
    // Path to the ROOT file
    const char* filePath = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP.root";

    // Open the ROOT file
    TFile* file = TFile::Open(filePath, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return;
    }

    // Get the TTree named "signalTree"
    TTree* signalTree = dynamic_cast<TTree*>(file->Get("signalTree"));
    if (!signalTree) {
        std::cerr << "Error: Could not find TTree 'signalTree' in the file." << std::endl;
        file->Close();
        return;
    }

    // Variables to hold the values of sigMET and sigMT
    float sigMET = 0;
    float sigMT = 0;

    // Set branch addresses
    signalTree->SetBranchAddress("sigMET", &sigMET);
    signalTree->SetBranchAddress("sigMT", &sigMT);

    // Create a 2D histogram
    TH2F* h2D = new TH2F("h2D", "MET vs MT, 2016 preVFP e#gamma events;sigMET;sigMT", 70, 0, 70, 100, 0, 100);

    // Loop over the entries in the TTree
    Long64_t nEntries = signalTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        signalTree->GetEntry(i);
        h2D->Fill(sigMET, sigMT);
    }

    // Draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "2D Plot of sigMET vs sigMT", 800, 600);
    h2D->Draw("COLZ");

    // Move stat box to avoid overlap with z-axis
    gStyle->SetOptStat(1111);
    h2D->SetStats(true);
    TPaveStats* stats = (TPaveStats*)h2D->GetListOfFunctions()->FindObject("stats");
    if (stats) {
        stats->SetX1NDC(0.7); // Move stat box to the right
        stats->SetX2NDC(0.9);
        stats->SetY1NDC(0.7);
        stats->SetY2NDC(0.9);
    }

    // Save the canvas as a PNG image
    canvas->SaveAs("sigMET_vs_sigMT_egamma.pdf");

}
