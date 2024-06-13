//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  2 04:11:41 2023 by ROOT version 6.12/07
// from TTree egTree/egTree
// found on file: /uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root
//////////////////////////////////////////////////////////

#ifndef egTree_h
#define egTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class egTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         crosssection;
   Float_t         ntotalevent;
   Float_t         ISRWeight;
   Float_t         pdfWeight;
   vector<float>   *pdfSystWeight;
   vector<float>   *ScaleSystWeight;
   Float_t         phoEt;
   Float_t         phoEta;
   Float_t         phoPhi;
   Float_t         lepPt;
   Float_t         lepEta;
   Float_t         lepPhi;
   Float_t         PUweight;
   Float_t         sigMT;
   Float_t         sigMET;
   Float_t         sigMETPhi;
   Float_t         dPhiLepMET;
   Int_t           nVertex;
   Float_t         dRPhoLep;
   Float_t         HT;
   Float_t         nJet;
   Int_t           nBJet;
   Int_t           nISRJet;
   Float_t         ISRJetPt;
   Float_t         invmass;
   Float_t         bosonPt;
   Float_t         sigMETJESup;
   Float_t         sigMETJESdo;
   Float_t         sigMETJERup;
   Float_t         sigMETJERdo;
   Float_t         sigMTJESup;
   Float_t         sigMTJESdo;
   Float_t         sigMTJERup;
   Float_t         sigMTJERdo;
   Float_t         dPhiLepMETJESup;
   Float_t         dPhiLepMETJESdo;
   Float_t         dPhiLepMETJERup;
   Float_t         dPhiLepMETJERdo;
   Float_t         HTJESup;
   Float_t         HTJESdo;
   Int_t           mcType;
   vector<int>     *mcPID;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcPt;
   vector<int>     *mcMomPID;
   vector<int>     *mcGMomPID;
   vector<int>     *mcStatus;
   Float_t         llmass;

   // List of branches
   TBranch        *b_crosssection;   //!
   TBranch        *b_ntotalevent;   //!
   TBranch        *b_ISRWeight;   //!
   TBranch        *b_pdfWeight;   //!
   TBranch        *b_pdfSystWeight;   //!
   TBranch        *b_ScaleSystWeight;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_lepPt;   //!
   TBranch        *b_lepEta;   //!
   TBranch        *b_lepPhi;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_sigMT;   //!
   TBranch        *b_sigMET;   //!
   TBranch        *b_sigMETPhi;   //!
   TBranch        *b_dPhiLepMET;   //!
   TBranch        *b_nVertex;   //!
   TBranch        *b_dRPhoLep;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_nBJet;   //!
   TBranch        *b_nISRJet;   //!
   TBranch        *b_ISRJetPt;   //!
   TBranch        *b_invmass;   //!
   TBranch        *b_bosonPt;   //!
   TBranch        *b_sigMETJESup;   //!
   TBranch        *b_sigMETJESdo;   //!
   TBranch        *b_sigMETJERup;   //!
   TBranch        *b_sigMETJERdo;   //!
   TBranch        *b_sigMTJESup;   //!
   TBranch        *b_sigMTJESdo;   //!
   TBranch        *b_sigMTJERup;   //!
   TBranch        *b_sigMTJERdo;   //!
   TBranch        *b_dPhiLepMETJESup;   //!
   TBranch        *b_dPhiLepMETJESdo;   //!
   TBranch        *b_dPhiLepMETJERup;   //!
   TBranch        *b_dPhiLepMETJERdo;   //!
   TBranch        *b_HTJESup;   //!
   TBranch        *b_HTJESdo;   //!
   TBranch        *b_mcType;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_llmass;   //!

   egTree(TTree *tree=0);
   virtual ~egTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef egTree_cxx
egTree::egTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root");
      }
      f->GetObject("egTree",tree);

   }
   Init(tree);
}

egTree::~egTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t egTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t egTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void egTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdfSystWeight = 0;
   ScaleSystWeight = 0;
   mcPID = 0;
   mcEta = 0;
   mcPhi = 0;
   mcPt = 0;
   mcMomPID = 0;
   mcGMomPID = 0;
   mcStatus = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("crosssection", &crosssection, &b_crosssection);
   fChain->SetBranchAddress("ntotalevent", &ntotalevent, &b_ntotalevent);
   fChain->SetBranchAddress("ISRWeight", &ISRWeight, &b_ISRWeight);
   fChain->SetBranchAddress("pdfWeight", &pdfWeight, &b_pdfWeight);
   fChain->SetBranchAddress("pdfSystWeight", &pdfSystWeight, &b_pdfSystWeight);
   fChain->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight, &b_ScaleSystWeight);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("lepPt", &lepPt, &b_lepPt);
   fChain->SetBranchAddress("lepEta", &lepEta, &b_lepEta);
   fChain->SetBranchAddress("lepPhi", &lepPhi, &b_lepPhi);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("sigMT", &sigMT, &b_sigMT);
   fChain->SetBranchAddress("sigMET", &sigMET, &b_sigMET);
   fChain->SetBranchAddress("sigMETPhi", &sigMETPhi, &b_sigMETPhi);
   fChain->SetBranchAddress("dPhiLepMET", &dPhiLepMET, &b_dPhiLepMET);
   fChain->SetBranchAddress("nVertex", &nVertex, &b_nVertex);
   fChain->SetBranchAddress("dRPhoLep", &dRPhoLep, &b_dRPhoLep);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("nBJet", &nBJet, &b_nBJet);
   fChain->SetBranchAddress("nISRJet", &nISRJet, &b_nISRJet);
   fChain->SetBranchAddress("ISRJetPt", &ISRJetPt, &b_ISRJetPt);
   fChain->SetBranchAddress("invmass", &invmass, &b_invmass);
   fChain->SetBranchAddress("bosonPt", &bosonPt, &b_bosonPt);
   fChain->SetBranchAddress("sigMETJESup", &sigMETJESup, &b_sigMETJESup);
   fChain->SetBranchAddress("sigMETJESdo", &sigMETJESdo, &b_sigMETJESdo);
   fChain->SetBranchAddress("sigMETJERup", &sigMETJERup, &b_sigMETJERup);
   fChain->SetBranchAddress("sigMETJERdo", &sigMETJERdo, &b_sigMETJERdo);
   fChain->SetBranchAddress("sigMTJESup", &sigMTJESup, &b_sigMTJESup);
   fChain->SetBranchAddress("sigMTJESdo", &sigMTJESdo, &b_sigMTJESdo);
   fChain->SetBranchAddress("sigMTJERup", &sigMTJERup, &b_sigMTJERup);
   fChain->SetBranchAddress("sigMTJERdo", &sigMTJERdo, &b_sigMTJERdo);
   fChain->SetBranchAddress("dPhiLepMETJESup", &dPhiLepMETJESup, &b_dPhiLepMETJESup);
   fChain->SetBranchAddress("dPhiLepMETJESdo", &dPhiLepMETJESdo, &b_dPhiLepMETJESdo);
   fChain->SetBranchAddress("dPhiLepMETJERup", &dPhiLepMETJERup, &b_dPhiLepMETJERup);
   fChain->SetBranchAddress("dPhiLepMETJERdo", &dPhiLepMETJERdo, &b_dPhiLepMETJERdo);
   fChain->SetBranchAddress("HTJESup", &HTJESup, &b_HTJESup);
   fChain->SetBranchAddress("HTJESdo", &HTJESdo, &b_HTJESdo);
   fChain->SetBranchAddress("mcType", &mcType, &b_mcType);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("llmass", &llmass, &b_llmass);
   Notify();
}

Bool_t egTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void egTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t egTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef egTree_cxx
