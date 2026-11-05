#include "TFile.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TArrayI.h"
#include <iostream>
#include <fstream>  
using namespace std;

void count() {
  // Open output file
  ofstream out("pMSSM_yields_2018.txt");
  if (!out.is_open()) {
    cout << "Failed to open output file!" << endl;
    return;
  }

  // Open the ROOT file
  TFile *file = TFile::Open("/uscms/home/tmishra/nobackup/pMSSM_output/PMSSM_set_1_prompt_3_UL18_yields.root", "READ");
  if (!file || file->IsZombie()) {
    out << "Error opening ROOT file!" << endl;
    out.close();
    return;
  }
  
  // Get the THnSparse object
  THnSparse *sparse = (THnSparse*)file->Get("Nominal");
  if (!sparse) {
    out << "Could not find THnSparse 'Nominal' in the file!" << endl;
    file->Close();
    out.close();
    return;
  }

  // Get the number of bins for each axis
  Int_t nbinsX = sparse->GetAxis(0)->GetNbins();
  Int_t nbinsY = sparse->GetAxis(1)->GetNbins();
  Int_t nbinsZ = sparse->GetAxis(2)->GetNbins();
  
  out << "Dimensions: " << nbinsX << " x " << nbinsY << " x " << nbinsZ << " bins" << 
endl;
  out << "Non-zero bin contents:" << endl;
  out << "Axis0\tAxis1\tAxis2\tContent" << endl;
  
  // Prepare for bin iteration
  TH1::SetDefaultSumw2();
  TArrayI bins(3);
  Long64_t nbins = sparse->GetNbins();
  
  for (Long64_t i = 0; i < nbins; i++) {
    Double_t content = sparse->GetBinContent(i, bins.GetArray());
    if (content == 0) continue;
    out << bins[0] << "\t" << bins[1] << "\t" << bins[2] << "\t" << content << endl;
  }

  // Clean up
  file->Close();
  out.close();
}
