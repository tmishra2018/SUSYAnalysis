void printPhoEtMinMax() {
    // List of tree names to process
    std::vector<std::string> treeNames = {"signalTree", "proxyTree", "jetTree", "fakeLepTree"};

    // Open the two ROOT files
    TFile *file1 = TFile::Open("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root");
    TFile *file2 = TFile::Open("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root");

    if (!file1 || !file2) {
        std::cout << "Error: One or both files could not be opened!" << std::endl;
        return;
    }

    // Loop over the tree names and process each tree
    for (const auto& treeName : treeNames) {
        std::cout << "Processing tree: " << treeName << std::endl;

        // Get the tree from both files
        TTree *tree1 = (TTree*)file1->Get(treeName.c_str());
        TTree *tree2 = (TTree*)file2->Get(treeName.c_str());

        if (!tree1 || !tree2) {
            std::cout << "Error: Tree '" << treeName << "' not found in one or both files!" << std::endl;
            continue;
        }

        // Declare the phoEt variable (type float)
        float phoEt1, phoEt2;

        // Set branch addresses for phoEt in both trees
        tree1->SetBranchAddress("phoEt", &phoEt1);
        tree2->SetBranchAddress("phoEt", &phoEt2);

        // Variables to track min and max values for both files
        float minPhoEt1 = FLT_MAX, maxPhoEt1 = -FLT_MAX;
        float minPhoEt2 = FLT_MAX, maxPhoEt2 = -FLT_MAX;

        // Loop through tree1
        int nEntries1 = tree1->GetEntries();
        for (int i = 0; i < nEntries1; ++i) {
            tree1->GetEntry(i);
            if (phoEt1 < minPhoEt1) minPhoEt1 = phoEt1;
            if (phoEt1 > maxPhoEt1) maxPhoEt1 = phoEt1;
        }

        // Loop through tree2
        int nEntries2 = tree2->GetEntries();
        for (int i = 0; i < nEntries2; ++i) {
            tree2->GetEntry(i);
            if (phoEt2 < minPhoEt2) minPhoEt2 = phoEt2;
            if (phoEt2 > maxPhoEt2) maxPhoEt2 = phoEt2;
        }

        // Print the results for both files
        std::cout << "File 1: /eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root" << std::endl;
        std::cout << "  Minimum phoEt: " << minPhoEt1 << std::endl;
        std::cout << "  Maximum phoEt: " << maxPhoEt1 << std::endl;

        std::cout << "File 2: /eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root" << std::endl;
        std::cout << "  Minimum phoEt: " << minPhoEt2 << std::endl;
        std::cout << "  Maximum phoEt: " << maxPhoEt2 << std::endl;

        std::cout << "----------------------------------------------------" << std::endl;
    }

    // Close the files
    file1->Close();
    file2->Close();
}
