#include <iostream>
#include <TFile.h>
#include <TTree.h>

void printEventContent(const char* filename, const char* treeName) {
    TFile file(filename, "READ");
    if (file.IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)file.Get(treeName);
    if (!tree) {
        std::cerr << "Error: TTree " << treeName << " not found in file " << filename << std::endl;
        return;
    }

    std::cout << "File: " << filename << ", TTree: " << treeName << ", Entries: " << tree->GetEntries() << std::endl;
    file.Close();
}

void compareEventContent() {
    const char* file1 = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_Muon20.root";
    const char* file2 = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root";
    const char* treeNames[] = {"signalTree", "proxyTree", "jetTree", "fakeLepTree", "hadronTree"};

    for (const char* treeName : treeNames) {
        TFile f1(file1, "READ");
        TFile f2(file2, "READ");

        TTree* tree1 = (TTree*)f1.Get(treeName);
        TTree* tree2 = (TTree*)f2.Get(treeName);

        if (!tree1 || !tree2) {
            std::cerr << "Error: TTree " << treeName << " not found in one of the files." << std::endl;
            continue;
        }

        Long64_t entries1 = tree1->GetEntries();
        Long64_t entries2 = tree2->GetEntries();

        std::cout << "TTree: " << treeName << std::endl;
        std::cout << "File 1 Entries: " << entries1 << std::endl;
        std::cout << "File 2 Entries: " << entries2 << std::endl;
        if (entries2 != 0) {
            std::cout << "Ratio (File 1 / File 2): " << static_cast<double>(entries1) / entries2 << std::endl;
        } else {
            std::cout << "Ratio (File 1 / File 2): Undefined (File 2 has zero entries)" << std::endl;
        }

        f1.Close();
        f2.Close();
    }
}

void run() {
    compareEventContent();
}
