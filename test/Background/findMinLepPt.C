#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <limits>

void findMinLepPt() {
    // List of ROOT files
    std::vector<std::string> files = {
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGJet130_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGJet40_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_DYJetsToLL_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTGJets_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTJets_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGToLNuG_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WWG_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WW_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZG_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZ_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_ZGToLLG_2016postVFP.root"

    };

    // List of TTrees to process
    const char* treeNames[] = {"mgTree"};
    //const char* treeNames[] = {"mgTree", "egTree"};

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
            std::cout << "  TTree: " << treeName << ", Minimum lepPt: " << minLepPt << std::endl;
        }

        // Close the file
        file->Close();
        delete file;
    }
}
