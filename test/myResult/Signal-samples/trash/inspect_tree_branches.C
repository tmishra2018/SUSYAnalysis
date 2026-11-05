// Usage:
//   root -l -b -q 'inspect_tree_branches.C("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2018.root", "SUSYtree")'

#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TBranch.h"
#include <iostream>
#include <iomanip>

void inspect_tree_branches(const char* infile, const char* treename = "SUSYtree") {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { std::cout << "Could not open " << infile << std::endl; return; }
    TTree *t = (TTree*)f->Get(treename);
    if (!t) { std::cout << "Tree '" << treename << "' not found." << std::endl; return; }

    std::cout << "\n=== ALL branches in tree '" << treename << "' (name : leaf type) ===\n";
    TObjArray *branches = t->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *b = (TBranch*)branches->At(i);
        std::string btitle = b->GetTitle();
        std::cout << std::setw(28) << std::left << b->GetName()
                  << "  title=" << btitle << "\n";
    }

    // Now print raw per-event values for the mass-related branches, side by side,
    // for the first N entries -- to visually check for garbage/misalignment vs.
    // genuinely continuous physics values.
    float MsGsQ(0), Mneutralino(0), Mchargino(0);
    t->SetBranchAddress("MsGsQ", &MsGsQ);
    if (t->GetBranch("Mneutralino")) t->SetBranchAddress("Mneutralino", &Mneutralino);
    if (t->GetBranch("Mchargino"))   t->SetBranchAddress("Mchargino", &Mchargino);

    std::cout << "\n=== Raw per-event values, first 40 entries, full precision ===\n";
    std::cout << std::setw(6) << "evt" << std::setw(16) << "MsGsQ"
               << std::setw(16) << "Mneutralino" << std::setw(16) << "Mchargino" << "\n";
    int n = std::min((Long64_t)40, t->GetEntries());
    for (int i = 0; i < n; i++) {
        t->GetEntry(i);
        std::cout << std::setw(6) << i
                   << std::setw(16) << std::fixed << std::setprecision(5) << MsGsQ
                   << std::setw(16) << Mneutralino
                   << std::setw(16) << Mchargino << "\n";
    }

    // Also: for a FIXED gluino mass (say the first one seen), print the spread of
    // Mneutralino/Mchargino for just that slice, to see if it's discrete or continuous
    // within a single, controlled gluino-mass bucket.
    t->GetEntry(0);
    int fixedG = (int) std::lround(MsGsQ);
    std::cout << "\n=== For events with MsGsQ rounding to " << fixedG
              << ", print Mneutralino/Mchargino values (first 40 matches) ===\n";
    int shown = 0;
    for (Long64_t i = 0; i < t->GetEntries() && shown < 40; i++) {
        t->GetEntry(i);
        if ((int)std::lround(MsGsQ) != fixedG) continue;
        std::cout << "evt=" << i << "  Mneutralino=" << Mneutralino
                  << "  Mchargino=" << Mchargino << "\n";
        shown++;
    }
}
