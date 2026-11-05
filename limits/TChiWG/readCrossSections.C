#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <iostream>

void readCrossSections() {
    // Load the required library
    gSystem->Load("../../lib/libAnaClasses.so");

    // Open the ROOT file
    TFile xSecFile("../../test/cross/susyCrossSection.root");
    
    if (xSecFile.IsZombie()) {
        std::cerr << "Error: Failed to open the ROOT file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1D *p_crosssection_tchiwg = (TH1D*)xSecFile.Get("p_charginoSec");
    TH1D *p_crosssection_t5wg   = (TH1D*)xSecFile.Get("p_gluinoxSec");
    TH1D *p_crosssection_t6wg   = (TH1D*)xSecFile.Get("p_squarkxSec");

    // Function to print histogram content
    auto printHistogram = [](TH1D* hist, const std::string& name) {
        if (!hist) {
            std::cerr << "Error: Histogram '" << name << "' not found!" << std::endl;
            return;
        }

        std::cout << "\nBin-wise content of " << name << ":\n";
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
            double binCenter = hist->GetXaxis()->GetBinCenter(bin);
            double binContent = hist->GetBinContent(bin);
            double binError = hist->GetBinError(bin);

            std::cout << "Bin " << bin << " (M = " << binCenter << " GeV): "
                      << "Content = " << binContent << ", Error = " << binError << std::endl;
        }
    };

    // Print histograms
    printHistogram(p_crosssection_tchiwg, "p_charginoSec");
    printHistogram(p_crosssection_t5wg, "p_gluinoxSec");
    printHistogram(p_crosssection_t6wg, "p_squarkxSec");

    // Close the file
    xSecFile.Close();
}
