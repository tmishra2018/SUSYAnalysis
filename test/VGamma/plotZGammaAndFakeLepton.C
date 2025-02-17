#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <algorithm> // For std::max

void plotZGammaAndFakeLepton() {
        gROOT->SetBatch(1);
    // File paths
    const char *file1 = "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_VGBkg_met0_70_pt0_1000_2018.root";
    const char *file2 = "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_qcd_met0_70_pt0_1000_2018.root";

    // Variables to plot
    const char *hist1Name = "p_dPhiEleMET_ZG"; // From file1
    const char *hist2Name = "p_dPhiEleMET";    // From file2

    // Open files and retrieve histograms
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);

    if (!f1 || !f2 || f1->IsZombie() || f2->IsZombie()) {
        std::cerr << "Error: One or both files could not be opened!" << std::endl;
        return;
    }

    TH1D *hZG = (TH1D *)f1->Get(hist1Name); // Zgamma
    TH1D *hFake = (TH1D *)f2->Get(hist2Name); // Fake lepton

    if (!hZG || !hFake) {
        std::cerr << "Error: Could not retrieve one or both histograms!" << std::endl;
        f1->Close();
        f2->Close();
        return;
    }

    // Set histogram styles
    hZG->SetLineColor(kRed);      // Zgamma
    hZG->SetLineWidth(2);
    hZG->SetTitle("Overlay of Zgamma and Fake Lepton;#Delta#phi(Ele, MET);Events");

    hFake->SetLineColor(kBlue);   // Fake lepton
    hFake->SetLineWidth(2);

    // Calculate global maximum and set axis limits
    double maxZG = hZG->GetMaximum();
    double maxFake = hFake->GetMaximum();
    double globalMax = std::max(maxZG, maxFake);
    double globalMin = 0.1; // Minimum y-value for better visibility

    hZG->SetMaximum(1.2 * globalMax); // Add margin above the maximum
    hZG->SetMinimum(globalMin);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay of Distributions", 800, 600);
    c1->cd();

    // Draw histograms
    hZG->Draw("HIST");         // Draw first histogram
    hFake->Draw("HIST SAME");  // Overlay second histogram

    // Add legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.4, 0.85); // Top-left corner
    legend->AddEntry(hZG, "Zgamma", "l");
    legend->AddEntry(hFake, "Fake lepton", "l");
    legend->Draw();

    // Save the canvas to a file
    c1->SaveAs("zgamma_fakelepton_overlay.png");

    // Clean up
    f1->Close();
    f2->Close();
    delete c1;
}
