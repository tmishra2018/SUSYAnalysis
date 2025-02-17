#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <algorithm> // For std::max

void plotDistributions() {
	    gROOT->SetBatch(1);
    // File paths
    const char *file1 = "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_VGBkg_met0_70_pt0_1000_2018.root";
    const char *file2 = "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_eleBkg_met0_70_pt0_1000_2018.root";
    const char *file3 = "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_jetbkg_met0_70_pt0_1000_2018.root";

    // Histogram names
    const char *histName = "p_dPhiEleMET";

    // Open files and retrieve histograms
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);

    if (!f1 || !f2 || !f3 || f1->IsZombie() || f2->IsZombie() || f3->IsZombie()) {
        std::cerr << "Error: One or more files could not be opened!" << std::endl;
        return;
    }

    TH1D *h1 = (TH1D *)f1->Get(histName);
    TH1D *h2 = (TH1D *)f2->Get(histName);
    TH1D *h3 = (TH1D *)f3->Get(histName);

    if (!h1 || !h2 || !h3) {
        std::cerr << "Error: Could not retrieve one or more histograms!" << std::endl;
        f1->Close();
        f2->Close();
        f3->Close();
        return;
    }

    // Set histogram styles
    h1->SetLineColor(kRed);    // WG/ZG (file1)
    h1->SetLineWidth(2);

    h2->SetLineColor(kBlue);   // Electron fake photon (file2)
    h2->SetLineWidth(2);

    h3->SetLineColor(kGreen);  // Jet fake photon (file3)
    h3->SetLineWidth(2);

    // Calculate global maximum and set axis limits
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double max3 = h3->GetMaximum();

    double globalMax = std::max({max1, max2, max3});
    double globalMin = 0.1; // Minimum y-value for better visibility of histograms

    h1->SetMaximum(1.2 * globalMax); // Add margin above the maximum
    h1->SetMinimum(globalMin);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay of Distributions", 800, 600);
    c1->cd();

    // Draw histograms
    h1->Draw("HIST");         // Draw first histogram
    h2->Draw("HIST SAME");    // Overlay second histogram
    h3->Draw("HIST SAME");    // Overlay third histogram

    // Add legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.4, 0.85); // Top-left corner
    legend->AddEntry(h1, "WG/ZG", "l");
    legend->AddEntry(h2, "Electron fake photon", "l");
    legend->AddEntry(h3, "Jet fake photon", "l");
    legend->Draw();

    // Save the canvas to a file
    c1->SaveAs("overlay_distributions.png");

    // Clean up
    f1->Close();
    f2->Close();
    f3->Close();
    delete c1;
}
