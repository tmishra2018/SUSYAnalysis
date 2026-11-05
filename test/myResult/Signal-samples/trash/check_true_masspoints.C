// Diagnostic: print the ACTUAL discrete (MsGsQ, Mneutralino/Mchargino) mass
// values present in the T5Wg ntuple, independent of any TH2D binning.
// This tells us the real underlying mass grid, so we're not guessing based
// on histogram bin centers (which is what got us into trouble last time).
//
// Usage -- point this at the RAW input tree analysis_TChiWG.C reads from
// (NOT signalTree_T5WG_*.root, which is the output). For T5WG that's:
//   /uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_<year><VFP>.root
// tree name "SUSYtree":
//   root -l -b -q 'check_true_masspoints.C("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2018.root", "SUSYtree")'

#include <map>
#include <set>
#include <utility>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"

void check_true_masspoints(const char* infile, const char* treename = "SUSYtree") {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { std::cout << "Could not open " << infile << std::endl; return; }
    TTree *t = (TTree*)f->Get(treename);
    if (!t) { std::cout << "Tree '" << treename << "' not found -- check the tree name." << std::endl; return; }

    float MsGsQ(0), Mneutralino(0), Mchargino(0);
    t->SetBranchAddress("MsGsQ", &MsGsQ);
    t->SetBranchAddress("Mneutralino", &Mneutralino);
    if (t->GetBranch("Mchargino")) t->SetBranchAddress("Mchargino", &Mchargino);

    // key = (gluino mass rounded to nearest integer, LSP mass rounded to nearest integer)
    std::map<std::pair<int,int>, long long> counts;
    std::set<int> gluinoRawSet;
    std::set<int> neutralinoRawSet;
    std::set<int> charginoRawSet;

    Long64_t nent = t->GetEntries();
    for (Long64_t i = 0; i < nent; i++) {
        t->GetEntry(i);
        // OLD (matches current analysis_TChiWG.C logic): prefer chargino if >0
        float lsp = (Mchargino > 0) ? Mchargino : Mneutralino;
        if (MsGsQ <= 0 || lsp <= 0) continue;
        int gRound = (int) std::lround(MsGsQ);   // round, don't truncate
        int lRound = (int) std::lround(lsp);
        counts[{gRound, lRound}]++;
        gluinoRawSet.insert(gRound);
        if (Mneutralino > 0) neutralinoRawSet.insert((int)std::lround(Mneutralino));
        if (Mchargino   > 0) charginoRawSet.insert((int)std::lround(Mchargino));
    }

    std::cout << "\n=== Distinct Mneutralino values (rounded), independent of Mchargino ===\n";
    std::cout << "count=" << neutralinoRawSet.size() << ": ";
    for (int v : neutralinoRawSet) std::cout << v << " ";
    std::cout << "\n\n=== Distinct Mchargino values (rounded), independent of Mneutralino ===\n";
    std::cout << "count=" << charginoRawSet.size() << ": ";
    for (int v : charginoRawSet) std::cout << v << " ";
    std::cout << "\n";

    std::cout << "\n=== Distinct gluino mass values seen (rounded to nearest GeV) ===\n";
    for (int g : gluinoRawSet) std::cout << g << " ";
    std::cout << "\n\n=== Sample of raw (un-rounded) MsGsQ near 1900-2100 GeV ===\n";
    std::cout << "(rerun of the loop, print first 20 raw values in that window)\n";
    t->ResetBranchAddresses();
    t->SetBranchAddress("MsGsQ", &MsGsQ);
    t->SetBranchAddress("Mneutralino", &Mneutralino);
    if (t->GetBranch("Mchargino")) t->SetBranchAddress("Mchargino", &Mchargino);
    int printed = 0;
    for (Long64_t i = 0; i < nent && printed < 20; i++) {
        t->GetEntry(i);
        if (MsGsQ > 1900 && MsGsQ < 2100) {
            std::cout << std::fixed << std::setprecision(4) << MsGsQ << "  ";
            printed++;
        }
    }
    std::cout << "\n\n=== Full (gluino, LSP) point counts, sorted ===\n";
    for (auto &kv : counts) {
        std::cout << "mGluino=" << kv.first.first
                   << "  mLSP=" << kv.first.second
                   << "  DM=" << (kv.first.first - kv.first.second)
                   << "  nEntries=" << kv.second << "\n";
    }
}
