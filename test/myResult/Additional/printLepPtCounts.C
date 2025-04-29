#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <vector>

void printLepPtCounts() {
    // Open the ROOT file
    TFile *file = TFile::Open("/eos/uscms/store/user/tmishra/InputFilesMC/WZG/WZG_2018.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return;
    }
    
    // Get the TTree
    TTree *tree = dynamic_cast<TTree*>(file->Get("ggNtuplizer/EventTree"));
    if (!tree) {
        std::cerr << "Error: TTree 'ggNtuplizer/EventTree' not found!" << std::endl;
        file->Close();
        return;
    }
    
    // Define a vector to hold muPt values
    std::vector<float> *muPt = nullptr;
    tree->SetBranchAddress("muPt", &muPt);
    
    // Loop over tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (muPt && muPt->size() > 2) {
            std::cout << "Event " << i << ": muPt count = " << muPt->size() << std::endl;
        }
    }
    
    // Cleanup
    file->Close();
    delete file;
}

