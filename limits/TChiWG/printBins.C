#include <TFile.h>
#include <TH1F.h>
#include <iostream>

void printBins() {  // Define the function only once
    TFile xSecFile("susyxSec.root");

    TH1F *p_crosssection = (TH1F*)xSecFile.Get("p_charginoSec");
    
    if (!p_crosssection) {
        std::cerr << "Error: Histogram 'p_charginoSec' not found!" << std::endl;
        return;
    }

    std::cout << "Bin-wise content of p_charginoSec:\n";
    for (int bin = 1; bin <= p_crosssection->GetNbinsX(); ++bin) {
        double binContent = p_crosssection->GetBinContent(bin);
        double binError = p_crosssection->GetBinError(bin);
        double binCenter = p_crosssection->GetXaxis()->GetBinCenter(bin);
        std::cout <<  " (M = " << binCenter << " GeV): Content = " 
                  << binContent  << std::endl;
    }

    xSecFile.Close();
}

// Automatically execute when running `root -l printBins.C`
void printBins();

