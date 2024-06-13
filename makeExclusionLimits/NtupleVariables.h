//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  4 01:48:45 2016 by ROOT version 6.06/01
// from TTree PreSelection/PreSelection
// found on file: root://cmseos.fnal.gov//store/user/vhegde/myProduction_V11/Spring16.WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_96_RA2AnalysisTree.root
//////////////////////////////////////////////////////////
#ifndef NtupleVariables_h
#define NtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;
class NtupleVariables : public TSelector {
 public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        limit;
   Double_t        limitErr;
   Double_t        mh;
   Int_t           syst;
   Int_t           iToy;
   Int_t           iSeed;
   Int_t           iChannel;
   Float_t         t_cpu;
   Float_t         t_real;
   Float_t         quantileExpected;

   /* double          GluinoMass; */
   /* double          NLSPMass; */
   // List of branches
   TBranch        *b_limit;   //!
   TBranch        *b_limitErr;   //!
   TBranch        *b_mh;   //!
   TBranch        *b_syst;   //!
   TBranch        *b_iToy;   //!
   TBranch        *b_iSeed;   //!
   TBranch        *b_iChannel;   //!
   TBranch        *b_t_cpu;   //!
   TBranch        *b_t_real;   //!
   TBranch        *b_quantileExpected;   //!
   
   /* TBranch        *b_GluinoMass;   //! */
   /* TBranch        *b_NLSPMass;   //! */
 NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   void    Init(TTree *tree, string);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   double  DeltaPhi(double, double);
   double  DeltaR(double eta1, double phi1, double eta2, double phi2);
};

#endif
#ifdef NtupleVariables_cxx
void NtupleVariables::Init(TTree *tree, string nameData)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

   // Set object pointer
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("limit", &limit, &b_limit);
   fChain->SetBranchAddress("limitErr", &limitErr, &b_limitErr);
   fChain->SetBranchAddress("mh", &mh, &b_mh);
   fChain->SetBranchAddress("syst", &syst, &b_syst);
   fChain->SetBranchAddress("iToy", &iToy, &b_iToy);
   fChain->SetBranchAddress("iSeed", &iSeed, &b_iSeed);
   fChain->SetBranchAddress("iChannel", &iChannel, &b_iChannel);
   fChain->SetBranchAddress("t_cpu", &t_cpu, &b_t_cpu);
   fChain->SetBranchAddress("t_real", &t_real, &b_t_real);
   fChain->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);

   Notify();
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef temp_cxx
