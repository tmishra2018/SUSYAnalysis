#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void plotHistograms() {
	gROOT->SetBatch(1);
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_VGBkg_met40_70_pt0_1000_2018.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }

    // Retrieve the histograms
    TH1D *hWG = (TH1D *)file->Get("p_dPhiEleMET_WG");
    TH1D *hZG = (TH1D *)file->Get("p_dPhiEleMET_ZG");

    if (!hWG) {
        std::cerr << "Error: Histogram 'p_dPhiEleMET_WG' not found in the file!" << std::endl;
        file->Close();
        return;
    }

    if (!hZG) {
        std::cerr << "Error: Histogram 'p_dPhiEleMET_ZG' not found in the file!" << std::endl;
        file->Close();
        return;
    }

    // Check if histograms have data
    std::cout << "Entries in p_dPhiEleMET_WG: " << hWG->GetEntries() << std::endl;
    std::cout << "Entries in p_dPhiEleMET_ZG: " << hZG->GetEntries() << std::endl;

    if (hWG->GetEntries() == 0) {
        std::cerr << "Warning: Histogram 'p_dPhiEleMET_WG' is empty!" << std::endl;
    }

    if (hZG->GetEntries() == 0) {
        std::cerr << "Warning: Histogram 'p_dPhiEleMET_ZG' is empty!" << std::endl;
    }

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay Histograms", 800, 600);
    c1->cd();

    // Set histogram styles
    hWG->SetLineColor(kRed);
    hWG->SetLineWidth(2);
    hWG->SetTitle("Overlay of WGamma and ZGamma Histograms;#Delta#phi(Ele, MET);Events");

    hZG->SetLineColor(kBlue);
    hZG->SetLineWidth(2);

    // Adjust the maximum
    double maxWG = hWG->GetMaximum();
    double maxZG = hZG->GetMaximum();
    double globalMax = std::max(maxWG, maxZG);

    hWG->SetMaximum(1.2 * globalMax); // Add some margin above the maximum
    hWG->SetMinimum(0.1);

    // Draw histograms
    hWG->Draw("HIST");
    hZG->Draw("HIST SAME");

    // Add a legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.85);
    legend->AddEntry(hWG, "WGamma", "l");
    legend->AddEntry(hZG, "ZGamma", "l");
    legend->Draw();

    // Save the canvas to a file
    c1->SaveAs("overlay_histograms_egamma.png");

    // Clean up
    file->Close();
    delete c1;
}
