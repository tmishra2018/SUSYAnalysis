#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <limits>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <limits>

void findMinLepPt2() {
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file." << std::endl;
        return;
    }

    // Names of the TTrees to process
    const char* treeNames[] = {"signalTree", "proxyTree", "jetTree", "fakeLepTree"};
    
    // Loop over each TTree
    for (const char* treeName : treeNames) {
        // Get the TTree
        TTree *tree = (TTree*)file->Get(treeName);
        if (!tree) {
            std::cerr << "Error: TTree " << treeName << " not found in the file." << std::endl;
            continue;
        }

        // Set up the branch
        float lepPt;
        tree->SetBranchAddress("lepPt", &lepPt);

        // Find the minimum value of lepPt
        float minLepPt = std::numeric_limits<float>::max();
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            if (lepPt < minLepPt) {
                minLepPt = lepPt;
            }
        }

        // Print the minimum value for this TTree
        std::cout << "Minimum lepPt in " << treeName << ": " << minLepPt << std::endl;
    }

    // Close the file
    file->Close();
    delete file;
}

