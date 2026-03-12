void printBranchDataType(const char* fileName, const char* treeName, const char* branchName) {
    // Open the ROOT file
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
    }

    // Get the TTree
    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Error: Could not find tree " << treeName << " in file " << fileName << std::endl;
        file->Close();
        return;
    }

    // Get the branch
    TBranch* branch = tree->GetBranch(branchName);
    if (!branch) {
        std::cerr << "Error: Could not find branch " << branchName << " in tree " << treeName << std::endl;
        file->Close();
        return;
    }

    // Get the leaf associated with the branch
    TLeaf* leaf = branch->GetLeaf(branchName);
    if (!leaf) {
        std::cerr << "Error: Could not find leaf for branch " << branchName << std::endl;
        file->Close();
        return;
    }

    // Print the data type of the branch using the leaf's type
    std::cout << "Data type of branch '" << branchName << "' in tree '" << treeName << "' from file '" << fileName << "': "
              << leaf->GetTypeName() << std::endl;

    // Close the file
    file->Close();
}

void checkNJetBranch() {
    // List of files to check
    const char* files[] = {
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_DYJetsToLL_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTGJets_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTJets_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGJet130_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGJet40_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WGToLNuG_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WWG_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WW_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZG_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZ_2018.root",
        "/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_ZGToLLG_2018.root"
    };

    // Check 'nJet' branch in egTree and mgTree separately for all files
    for (int i = 0; i < sizeof(files)/sizeof(files[0]); ++i) {
        std::cout << "Checking file: " << files[i] << std::endl;

        // For 'egTree'
        printBranchDataType(files[i], "egTree", "nISRJet");

        // For 'mgTree'
        printBranchDataType(files[i], "mgTree", "nISRJet");

        std::cout << std::endl; // To separate outputs
    }
}
