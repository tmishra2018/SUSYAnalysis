#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include "TROOT.h"

using namespace std;

void PUReweight(int RunYear, const char *whichVFP)
{
    gROOT->SetBatch(kTRUE);
    TFile *fMC = TFile::Open(Form("/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/DYJetsToLL/DYJetsToLL_%d%s.root",RunYear, whichVFP));
    TTree *tree = (TTree *)fMC->Get("ggNtuplizer/EventTree");
    
    vector<int>   *puBX   = nullptr;
    vector<float> *puTrue = nullptr;

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("puBX",1);
    tree->SetBranchStatus("puTrue",1);

    tree->SetBranchAddress("puBX",   &puBX);
    tree->SetBranchAddress("puTrue", &puTrue);

    const int   nBins = 200;
    const float puMax = 200.0;

    TH1F *hPUTrue_MC = new TH1F("hPUTrue_MC", "MC true PU", nBins, 0, puMax);
    hPUTrue_MC->Sumw2();

    Long64_t nEntries = tree->GetEntries();
    cout << "Processing " << nEntries << " events" << endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        if (!puBX || !puTrue) continue;
        if (puBX->size() != puTrue->size()) continue;

        for (size_t j = 0; j < puBX->size(); j++) {
            if (puBX->at(j) == 0) {   // in-time BX
                hPUTrue_MC->Fill(puTrue->at(j));
                break;
            }
        }
    }

    if (hPUTrue_MC->Integral() > 0)
        hPUTrue_MC->Scale(1.0 / hPUTrue_MC->Integral());

TFile *fDataNom = TFile::Open(Form("PUfiles/dataPU_%d%s.root",RunYear, whichVFP));
TH1F *hPUTrue_Data   = (TH1F*) fDataNom->Get("pileup");
hPUTrue_Data->Scale(1.0 / hPUTrue_Data->Integral());
TH1F *hPUWeight     = (TH1F*)hPUTrue_Data->Clone("PUweight");

    hPUWeight->Reset();

    for (int i = 1; i <= nBins; i++) {
        double mc = hPUTrue_MC->GetBinContent(i);
        if (mc <= 0) continue;
        hPUWeight->SetBinContent(i, hPUTrue_Data->GetBinContent(i) / mc);
    }

    TFile *fOut = new TFile(Form("PUfiles/PUWeights_%d%s.root", RunYear, whichVFP), "RECREATE");
    hPUTrue_MC->Write();
    hPUTrue_Data->Write();
    hPUWeight->Write();
    fOut->Close();
    cout << "PU reweighting done successfully ✔" << endl;
}
int main(int argc, char** argv)
{
	    PUReweight(atoi(argv[1]),argv[2]);
		    return 0;
}
