#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace std;

// ============================================
// Helper: DeltaR
// ============================================
double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = fabs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  double deta = eta1 - eta2;
  return sqrt(deta*deta + dphi*dphi);
}

// ============================================
// Helper: prompt lepton origin
// ============================================
bool isFromHardProcess(int momPID, int gmomPID) {
  int amom  = abs(momPID);
  int agmom = abs(gmomPID);

  // W / Z
  if (amom == 23 || amom == 24) return true;
  if (agmom == 23 || agmom == 24) return true;

  // top
  if (amom == 6 || agmom == 6) return true;

  return false;
}

// ============================================
// Helper: prompt photon definition
// ============================================
bool isPromptPhoton(int momPID, int gmomPID) {
  int amom  = abs(momPID);
  int agmom = abs(gmomPID);

  // reject hadron decays
  if (amom == 111 || amom == 221) return false;   // pi0, eta
  if (agmom == 111 || agmom == 221) return false;

  // accept photons from hard partons or leptons
  if (amom <= 6 || amom == 11 || amom == 13) return true;
  if (agmom <= 6 || agmom == 11 || agmom == 13) return true;

  return false;
}

// ============================================
// Main analysis
// ============================================
void analysis_doubleFake(const char *infile = "input.root") {

  TFile *fin = TFile::Open(infile);
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input file" << endl;
    return;
  }

  TTree *tree = (TTree*)fin->Get("sigtree");
  if (!tree) {
    cout << "ERROR: sigtree not found" << endl;
    return;
  }

  // ============================================
  // Branches (following analysis_fake.C style)
  // ============================================
  float lepEta, lepPhi;
  float phoEta, phoPhi;

  vector<int>   *mcPID     = 0;
  vector<int>   *mcMomPID  = 0;
  vector<int>   *mcGMomPID = 0;
  vector<float> *mcEta     = 0;
  vector<float> *mcPhi     = 0;

  tree->SetBranchAddress("lepEta", &lepEta);
  tree->SetBranchAddress("lepPhi", &lepPhi);
  tree->SetBranchAddress("phoEta", &phoEta);
  tree->SetBranchAddress("phoPhi", &phoPhi);

  tree->SetBranchAddress("mcPID", &mcPID);
  tree->SetBranchAddress("mcMomPID", &mcMomPID);
  tree->SetBranchAddress("mcGMomPID", &mcGMomPID);
  tree->SetBranchAddress("mcEta", &mcEta);
  tree->SetBranchAddress("mcPhi", &mcPhi);

  // ============================================
  // Counters (item 4)
  // ============================================
  Long64_t nTotal        = 0;
  Long64_t nFakeLep      = 0;
  Long64_t nFakePho      = 0;
  Long64_t nDoubleFake   = 0;

  // ============================================
  // Event loop (item 5)
  // ============================================
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ievt++) {
    tree->GetEntry(ievt);
    nTotal++;

    bool lepPrompt = false;
    bool phoPrompt = false;

    // ---- electron truth matching ----
    for (unsigned i = 0; i < mcPID->size(); i++) {
      if (abs(mcPID->at(i)) != 11) continue;

      double dR = DeltaR(mcEta->at(i), mcPhi->at(i), lepEta, lepPhi);
      if (dR > 0.3) continue;

      if (isFromHardProcess(mcMomPID->at(i), mcGMomPID->at(i))) {
        lepPrompt = true;
        break;
      }
    }

    // ---- photon truth matching ----
    for (unsigned i = 0; i < mcPID->size(); i++) {
      if (mcPID->at(i) != 22) continue;

      double dR = DeltaR(mcEta->at(i), mcPhi->at(i), phoEta, phoPhi);
      if (dR > 0.3) continue;

      if (isPromptPhoton(mcMomPID->at(i), mcGMomPID->at(i))) {
        phoPrompt = true;
        break;
      }
    }

    bool lepFake = !lepPrompt;
    bool phoFake = !phoPrompt;

    if (lepFake) nFakeLep++;
    if (phoFake) nFakePho++;
    if (lepFake && phoFake) nDoubleFake++;
  }

  // ============================================
  // Final quantification (item 6)
  // ============================================
  cout << "======================================" << endl;
  cout << "Total events           : " << nTotal << endl;
  cout << "Fake lepton events     : " << nFakeLep << endl;
  cout << "Fake photon events     : " << nFakePho << endl;
  cout << "Double fake events     : " << nDoubleFake << endl;

  double fracDoubleAll = (double)nDoubleFake / nTotal;
  double fracDoubleInFakePho =
      (nFakePho > 0 ? (double)nDoubleFake / nFakePho : 0.0);

  cout << "Double fake fraction (all events)      : "
       << fracDoubleAll << endl;
  cout << "Double fake fraction (within fake pho) : "
       << fracDoubleInFakePho << endl;

  cout << "======================================" << endl;

  fin->Close();
}
