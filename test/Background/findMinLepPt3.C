#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <limits>

void findMinLepPt3() {
    std::vector<std::string> files = {
	"/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_DY_2016postVFP.root",
      	"/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_DY_2017.root",
	"/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_DY_2016preVFP.root",
     	"/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_DY_2018.root"

    };

    // List of TTrees to process
    const char* treeNames[] = {"ZTree", "mgTree"};

    // Loop over each file
    for (const std::string& fileName : files) {
        // Open the ROOT file
        TFile* file = TFile::Open(fileName.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Unable to open file " << fileName << std::endl;
            continue;
        }

        std::cout << "File: " << fileName << std::endl;

        // Loop over each TTree
        for (const char* treeName : treeNames) {
            // Get the TTree
            TTree* tree = (TTree*)file->Get(treeName);
            if (!tree) {
                std::cerr << "  Error: TTree " << treeName << " not found in the file." << std::endl;
                continue;
            }

            // Set up the branch
            float muPt;
            tree->SetBranchAddress("muPt", &muPt);

            // Find the minimum value of lepPt
            float minLepPt = std::numeric_limits<float>::max();
            Long64_t nEntries = tree->GetEntries();
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                if (muPt < minLepPt) {
                    minLepPt = muPt;
                }
            }

            // Print the minimum value for this TTree
            std::cout << "  TTree: " << treeName << ", Minimum lepPt: " << minLepPt << std::endl;
        }

        // Close the file
        file->Close();
        delete file;
    }
}
