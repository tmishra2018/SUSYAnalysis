// Run against the OUTPUT file you already have, e.g.:
//   root -l -b -q 'check_susymass_hist.C("/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_2018.root")'
//
// This tells us directly whether the fixed (extended, 1975-split) binning is
// actually the one embedded in this file, or whether it's still the old
// 13-bin array -- without needing to touch the raw ntuple at all.

#include "TFile.h"
#include "TH2D.h"
#include <iostream>
#include <iomanip>

void check_susymass_hist(const char* infile) {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { std::cout << "Could not open " << infile << std::endl; return; }
    TH2D *h = (TH2D*)f->Get("SUSYMass");
    if (!h) { std::cout << "SUSYMass histogram not found in file." << std::endl; return; }

    int nx = h->GetXaxis()->GetNbins();
    std::cout << "\n=== SUSYMass x-axis: " << nx << " bins ===\n";
    std::cout << "(old/broken binning has 13 bins ending at 2125; "
                 "fixed/extended binning has 28 bins ending at 2825)\n\n";

    std::cout << "Bin edges:\n";
    for (int i = 1; i <= nx + 1; i++) {
        std::cout << h->GetXaxis()->GetBinLowEdge(i) << " ";
    }
    std::cout << "\n\nBin-by-bin: low edge, high edge, center, content (for bins with any entries):\n";
    for (int i = 1; i <= nx; i++) {
        double lo = h->GetXaxis()->GetBinLowEdge(i);
        double hi = h->GetXaxis()->GetBinUpEdge(i);
        double c  = h->GetXaxis()->GetBinCenter(i);
        double sumInY = 0;
        for (int j = 1; j <= h->GetYaxis()->GetNbins(); j++) sumInY += h->GetBinContent(i,j);
        std::cout << "bin " << i << ":  [" << lo << ", " << hi << ")  center=" << c
                   << "  totalEntriesThisGluinoMass=" << sumInY << "\n";
    }
}
