#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"


void analysis_SMMC(){//main  
  gSystem->Load("../lib/libAnaClasses.so");
  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/WGToLNuG/WGToLNuG_2017.root");

  std::ostringstream outputname;
  outputname << "/eos/uscms/store/user/tmishra/MC/SM_MC_WGToLNuG_2017.root";
  TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
  outputfile->cd();
  TTree *tree = new TTree("MTtree","MTtree");
  double MT_(0), ThreeBodyMass_(0);
  tree->Branch("MT",&MT_);
  tree->Branch("ThreeBodyMass",&ThreeBodyMass_);

 /// histo list
  TH1F *p_unprescaleMCHLT = new TH1F("p_unprescaleMCHLT","passed HLT event; HLT; entries",38,0,38);
  TH1F *p_prescaleMCHLT = new TH1F("p_prescaleMCHLT","passed HLT event; HLT; entries",38,0,38);

  TH1F *p_mcphotonET = new TH1F("photonET","signal #gamma (#tilde{#chi}^{0}#rightarrow#gamma) E_{T}; E_{T} (GeV)",500,0,1500);
  TH1F *p_mcphotonEta = new TH1F("p_mcphotonEta","signal #gamma (#tilde{#chi}^{0}#rightarrow#gamma) #eta; #eta;",60,-3,3);

  TProfile *p_phoMatchEff= new TProfile("p_phoMatchEff","p_phoMatchEff",150,20,320);

  TH1F *p_mceleET = new TH1F("eleET","signal e (#tilde{#chi}^{#pm}#rightarrowW^{#pm}#tilde{#chi}^{0}, W^{#pm}#rightarrowe) E_{T}; E_{T} (GeV)",500,0,1000);
  TH1F *p_mceleEta = new TH1F("p_mceleEta","signal e (#tilde{#chi}^{#pm}#rightarrowW^{#pm}#tilde{#chi}^{0}, W^{#pm}#rightarrowe) #eta; #eta;",60,-3,3);

  TH1F *p_phoEBR9_true = new TH1F("p_phoEBR9","gen #gamma (Barrel) R9;R9;",100,0,1);
  TH1F *p_phoEER9_true = new TH1F("p_phoEER9","gen #gamma (Endcap) R9;R9;",100,0,1);
  TH1F *p_phoEBHoverE_true = new TH1F("p_phoEBHoverE","gen #gamma (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_phoEEHoverE_true = new TH1F("p_phoEEHoverE","gen #gamma (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_phoEBsigma_true = new TH1F("p_phoEBsigma","gen #gamma (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_phoEEsigma_true = new TH1F("p_phoEEsigma","gen #gamma (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_phoEBChIso_true = new TH1F("p_phoEBChIso","gen #gamma (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEEChIso_true = new TH1F("p_phoEEChIso","gen #gamma (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEBNeuIso_true = new TH1F("p_phoEBNeuIso","gen #gamma (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEENeuIso_true = new TH1F("p_phoEENeuIso","gen #gamma (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEBPhoIso_true = new TH1F("p_phoEBPhoIso","gen #gamma (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_phoEEPhoIso_true = new TH1F("p_phoEEPhoIso","gen #gamma (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  
  TH1F *p_eleEBR9_true = new TH1F("p_eleEBR9","gen e (Barrel) R9;R9;",100,0,1);
  TH1F *p_eleEER9_true = new TH1F("p_eleEER9","gen e (Endcap) R9;R9;",100,0,1);
  TH1F *p_eleEBHoverE_true = new TH1F("p_eleEBHoverE","gen e (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_eleEEHoverE_true = new TH1F("p_eleEEHoverE","gen e (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_eleEBsigma_true = new TH1F("p_eleEBsigma","gen e (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_eleEEsigma_true = new TH1F("p_eleEEsigma","gen e (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_eleEBChIso_true = new TH1F("p_eleEBChIso","gen e (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEEChIso_true = new TH1F("p_eleEEChIso","gen e (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEBNeuIso_true = new TH1F("p_eleEBNeuIso","gen e (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEENeuIso_true = new TH1F("p_eleEENeuIso","gen e (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEBPhoIso_true = new TH1F("p_eleEBPhoIso","gen e (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_eleEEPhoIso_true = new TH1F("p_eleEEPhoIso","gen e (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_eleEBRelIso_true = new TH1F("p_eleEBRelIso_true","gen e (Barrel) Relative Iso; Iso_{Rel}(GeV);",60,0,0.3);
  TH1F *p_eleEERelIso_true = new TH1F("p_eleEERelIso_true","gen e (Endcap) Relative Iso; Iso_{Rel}(GeV);",60,0,0.3);
  TH1F *p_eleEBdEtaIn_true = new TH1F("p_eleEBdEtaIn_true","gen e (Barrel) dEtaIn; dEtaIn;",50,0,0.02);
  TH1F *p_eleEEdEtaIn_true = new TH1F("p_eleEEdEtaIn_true","gen e (Endcap) dEtaIn; dEtaIn;",50,0,0.02);
  TH1F *p_eleEBdPhiIn_true = new TH1F("p_eleEBdPhiIn_true","gen e (Barrel) dPhiIn; dPhiIn;",60,0,0.3);
  TH1F *p_eleEEdPhiIn_true = new TH1F("p_eleEEdPhiIn_true","gen e (Endcap) dPhiIn; dPhiIn;",60,0,0.3);
  TH1F *p_eleEBD0_true = new TH1F("p_eleEBD0_true","gen e (Barrel) d0; d0;",60,0,0.06);
  TH1F *p_eleEBDz_true = new TH1F("p_eleEBDz_true","gen e (Barrel) dz; dz;",50,0,0.5);
  TH1F *p_eleEBooEmooP_true = new TH1F("p_eleEBooEmooP_true","gen e (Barrel) 1/E-1/P; 1/E-1/P;",60,0,0.3);
  TH1F *p_eleEED0_true = new TH1F("p_eleEED0_true","gen e (Endcap) d0; d0;",60,0,0.06);
  TH1F *p_eleEEDz_true = new TH1F("p_eleEEDz_true","gen e (Endcap) dz; dz;",50,0,0.5);
  TH1F *p_eleEEooEmooP_true = new TH1F("p_eleEEooEmooP_true","gen e (Endcap) 1/E-1/P; 1/E-1/P;",60,0,0.3);

  TH1F *p_phoEBR9_ptselect = new TH1F("p_phoEBR9_ptselect","gen #gamma (Barrel) R9;R9;",100,0,1);
  TH1F *p_phoEER9_ptselect = new TH1F("p_phoEER9_ptselect","gen #gamma (Endcap) R9;R9;",100,0,1);
  TH1F *p_phoEBHoverE_ptselect = new TH1F("p_phoEBHoverE_ptselect","gen #gamma (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_phoEEHoverE_ptselect = new TH1F("p_phoEEHoverE_ptselect","gen #gamma (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_phoEBsigma_ptselect = new TH1F("p_phoEBsigma_ptselect","gen #gamma (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_phoEEsigma_ptselect = new TH1F("p_phoEEsigma_ptselect","gen #gamma (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_phoEBChIso_ptselect = new TH1F("p_phoEBChIso_ptselect","gen #gamma (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEEChIso_ptselect = new TH1F("p_phoEEChIso_ptselect","gen #gamma (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEBNeuIso_ptselect = new TH1F("p_phoEBNeuIso_ptselect","gen #gamma (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEENeuIso_ptselect = new TH1F("p_phoEENeuIso_ptselect","gen #gamma (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEBPhoIso_ptselect = new TH1F("p_phoEBPhoIso_ptselect","gen #gamma (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_phoEEPhoIso_ptselect = new TH1F("p_phoEEPhoIso_ptselect","gen #gamma (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  
  TH1F *p_eleEBR9_ptselect = new TH1F("p_eleEBR9_ptselect","gen e (Barrel) R9;R9;",100,0,1);
  TH1F *p_eleEER9_ptselect = new TH1F("p_eleEER9_ptselect","gen e (Endcap) R9;R9;",100,0,1);
  TH1F *p_eleEBHoverE_ptselect = new TH1F("p_eleEBHoverE_ptselect","gen e (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_eleEEHoverE_ptselect = new TH1F("p_eleEEHoverE_ptselect","gen e (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_eleEBsigma_ptselect = new TH1F("p_eleEBsigma_ptselect","gen e (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_eleEEsigma_ptselect = new TH1F("p_eleEEsigma_ptselect","gen e (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_eleEBChIso_ptselect = new TH1F("p_eleEBChIso_ptselect","gen e (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEEChIso_ptselect = new TH1F("p_eleEEChIso_ptselect","gen e (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEBNeuIso_ptselect = new TH1F("p_eleEBNeuIso_ptselect","gen e (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEENeuIso_ptselect = new TH1F("p_eleEENeuIso_ptselect","gen e (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEBPhoIso_ptselect = new TH1F("p_eleEBPhoIso_ptselect","gen e (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_eleEEPhoIso_ptselect = new TH1F("p_eleEEPhoIso_ptselect","gen e (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_eleEBRelIso_ptselect = new TH1F("p_eleEBRelIso_ptselect","gen e (Barrel) Relative Iso; Iso_{Rel}(GeV);",60,0,0.3);
  TH1F *p_eleEERelIso_ptselect = new TH1F("p_eleEERelIso_ptselect","gen e (Endcap) Relative Iso; Iso_{Rel}(GeV);",60,0,0.3);
  TH1F *p_eleEBdEtaIn_ptselect = new TH1F("p_eleEBdEtaIn_ptselect","gen e (Barrel) dEtaIn; dEtaIn;",50,0,0.02);
  TH1F *p_eleEEdEtaIn_ptselect = new TH1F("p_eleEEdEtaIn_ptselect","gen e (Endcap) dEtaIn; dEtaIn;",50,0,0.02);
  TH1F *p_eleEBdPhiIn_ptselect = new TH1F("p_eleEBdPhiIn_ptselect","gen e (Barrel) dPhiIn; dPhiIn;",60,0,0.3);
  TH1F *p_eleEEdPhiIn_ptselect = new TH1F("p_eleEEdPhiIn_ptselect","gen e (Endcap) dPhiIn; dPhiIn;",60,0,0.3);
  TH1F *p_eleEBD0_ptselect = new TH1F("p_eleEBD0_ptselect","gen e (Barrel) d0; d0;",60,0,0.06);
  TH1F *p_eleEBDz_ptselect = new TH1F("p_eleEBDz_ptselect","gen e (Barrel) dz; dz;",50,0,0.5);
  TH1F *p_eleEBooEmooP_ptselect = new TH1F("p_eleEBooEmooP_ptselect","gen e (Barrel) 1/E-1/P; 1/E-1/P;",60,0,0.3);
  TH1F *p_eleEED0_ptselect = new TH1F("p_eleEED0_ptselect","gen e (Endcap) d0; d0;",60,0,0.06);
  TH1F *p_eleEEDz_ptselect = new TH1F("p_eleEEDz_ptselect","gen e (Endcap) dz; dz;",50,0,0.5);
  TH1F *p_eleEEooEmooP_ptselect = new TH1F("p_eleEEooEmooP_ptselect","gen e (Endcap) 1/E-1/P; 1/E-1/P;",60,0,0.3);
  TH1F *p_eleConvVeto_ptselect = new TH1F("p_eleConvVeto_ptselect","gen e Conversion veto",5,0,5);
  TH2F *p_eleEBR9vsPt = new TH2F("p_eleEBR9vsPt","gen e (Barrel) R9 vs P_{T}; P_{T}(GeV); R9",250,0,500,20,0,1);
  TH2F *p_eleEER9vsPt = new TH2F("p_eleEER9vsPt","gen e (Endcap) R9 vs P_{T}; P_{T}(GeV); R9",250,0,500,20,0,1);

  TH1F *p_phoEBR9_idselect = new TH1F("p_phoEBR9_idselect","gen #gamma (Barrel) R9;R9;",100,0,1);
  TH1F *p_phoEER9_idselect = new TH1F("p_phoEER9_idselect","gen #gamma (Endcap) R9;R9;",100,0,1);
  TH1F *p_phoEBHoverE_idselect = new TH1F("p_phoEBHoverE_idselect","gen #gamma (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_phoEEHoverE_idselect = new TH1F("p_phoEEHoverE_idselect","gen #gamma (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_phoEBsigma_idselect = new TH1F("p_phoEBsigma_idselect","gen #gamma (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_phoEEsigma_idselect = new TH1F("p_phoEEsigma_idselect","gen #gamma (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_phoEBChIso_idselect = new TH1F("p_phoEBChIso_idselect","gen #gamma (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEEChIso_idselect = new TH1F("p_phoEEChIso_idselect","gen #gamma (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_phoEBNeuIso_idselect = new TH1F("p_phoEBNeuIso_idselect","gen #gamma (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEENeuIso_idselect = new TH1F("p_phoEENeuIso_idselect","gen #gamma (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_phoEBPhoIso_idselect = new TH1F("p_phoEBPhoIso_idselect","gen #gamma (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_phoEEPhoIso_idselect = new TH1F("p_phoEEPhoIso_idselect","gen #gamma (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  
  TH1F *p_eleEBR9_idselect = new TH1F("p_eleEBR9_idselect","gen e (Barrel) R9;R9;",100,0,1);
  TH1F *p_eleEER9_idselect = new TH1F("p_eleEER9_idselect","gen e (Endcap) R9;R9;",100,0,1);
  TH1F *p_eleEBHoverE_idselect = new TH1F("p_eleEBHoverE_idselect","gen e (Barrel) H/E;H/E;",51,-0.001,0.05);
  TH1F *p_eleEEHoverE_idselect = new TH1F("p_eleEEHoverE_idselect","gen e (Endcap) H/E;H/E;",51,-0.01,0.05);
  TH1F *p_eleEBsigma_idselect = new TH1F("p_eleEBsigma_idselect","gen e (Barrel) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",40,0,0.02);
  TH1F *p_eleEEsigma_idselect = new TH1F("p_eleEEsigma_idselect","gen e (Endcap) #sigma_{i#etai#eta};#sigma_{i#etai#eta};",50,0,0.05);
  TH1F *p_eleEBChIso_idselect = new TH1F("p_eleEBChIso_idselect","gen e (Barrel) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEEChIso_idselect = new TH1F("p_eleEEChIso_idselect","gen e (Endcap) Charged Hadron Iso; Iso_{Ch}(GeV);",62,-0.1,3);
  TH1F *p_eleEBNeuIso_idselect = new TH1F("p_eleEBNeuIso_idselect","gen e (Barrel) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEENeuIso_idselect = new TH1F("p_eleEENeuIso_idselect","gen e (Endcap) Neutral Hadron Iso;Iso_{Neu}(GeV);",61,-0.1,6);
  TH1F *p_eleEBPhoIso_idselect = new TH1F("p_eleEBPhoIso_idselect","gen e (Barrel) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);
  TH1F *p_eleEEPhoIso_idselect = new TH1F("p_eleEEPhoIso_idselect","gen e (Endcap) Photon Iso;Iso_{Pho}(GeV);",62,-0.1,3);

  TH2F *p_SUSYMass = new TH2F("Mass","Mass; M_{#tilde{#chi}^{0}/#tilde{#chi}^{#pm}} (GeV); M_{#tilde{g}} (GeV)",32,0,1600,16,800,1600);
  TH1F *p_Mt = new TH1F("Mt","M_{T}; M_{T} (GeV);",200,0,400);
  TH1F *p_PhoELeDeltaR = new TH1F("p_PhoELeDeltaR","#DeltaR(e, #gamma); #DeltaR(e, #gamma);",100,0,10);
  TH1F *p_PhoEleMass = new TH1F("p_PhoEleMass","M_{e#gamma}; M_{e#gamma}(GeV);", 200,0,400);
  TH2F *p_PhoEleDRvsM = new TH2F("p_PhoEleDRvsM",";M_{e#gamma}(GeV);#DeltaR(e, #gamma)",200,0,400,100,0,5);
  TH1F *p_selectPhoEleMass = new TH1F("p_selectPhoEleMass","M_{e#gamma}; M_{e#gamma}(GeV);", 200,0,400);
  TProfile *p_eleIDeff = new TProfile("p_eleIDeff","e#pm ID efficiency; p_{T}(GeV);",500,0,500);  
  TProfile *p_phoIDeff = new TProfile("p_phoIDeff","#gamma ID efficiency; E_{T}(GeV);",500,0,500);  
  
  TH1F *p_checkNeuEleMass = new TH1F("p_checkNeuEleMass","p_checkNeuEleMass",200,0,200);
  TH1F *p_checkPhoMatch = new TH1F("p_checkPhoMatch","p_checkPhoMatch",10,0,10);
  TH1F *p_checkEleMatch = new TH1F("p_checkEleMatch","p_checkEleMatch",10,0,10);

  TH1F *p_checkMETdiff = new TH1F("p_checkMETdiff","p_checkMETdiff",500,-250,250);

  TH1F *p_threebodymass = new TH1F("p_threebodymass","p_threebodymass",1000,0,1000);
  TH1F *p_threebodymass_W = new TH1F("p_threebodymass_W","p_threebodymass_W",1000,0,1000);
  TH2F *p_threebodymass_MT = new TH2F("p_threebodymass_MT","p_threebodymass_MT",1000,0,1000,1000,0,1000);
  RunType datatype(MC);
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoEle> Ele;
  float MET(0);
  float METPhi(0);

    const unsigned nEvts = es->GetEntries(); 
    //const unsigned nEvts = 10000; 
    std::cout << "total event : " << nEvts << std::endl;

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

        // Get data         
        raw.GetData(es, ievt);
        MCData.clear();
        Photon.clear();
        Ele.clear();
        for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
        for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
        MET = raw.pfMET;
        METPhi = raw.pfMETPhi;      

        //check HLT bit
        for(int iHLT(0); iHLT<38; iHLT++){
          if(((raw.HLTPho >> iHLT) &1) ==1){
            if(((raw.HLTPhoIsPrescaled >> iHLT) &1) ==0)p_unprescaleMCHLT->Fill(iHLT);
            else if(((raw.HLTPhoIsPrescaled >> iHLT) &1) ==1)p_prescaleMCHLT->Fill(iHLT);
          }
        }

	  std::vector<mcData>::iterator PhoLeg;
	  std::vector<mcData>::iterator EleLeg1;
          std::vector<mcData>::iterator NvLeg;
	  int nPhoLeg(0), nEleLeg1(0), nNvLeg(0);
 
          std::vector< std::vector<mcData>::iterator> PhoCollection;
          PhoCollection.clear();
          float sum_ex(0),sum_ey(0);
          bool initialPho(true);
	  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
	    if(itMC->getPt() < 0.001)continue;
            if(itMC->getEta() > 5.0)continue; 
            if(fabs(itMC->getPID()) != 12 && fabs(itMC->getPID()) != 14 && fabs(itMC->getPID()) != 16){
		// 12, 14, 16 for neutrinos
              float pt = itMC->getPt();
              float phi = itMC->getPhi();
              sum_ex += pt*cos(phi);
              sum_ey += pt*sin(phi);
            }
		// 11 for electron, 12 for neutrino, 24 for W+, 22 for gamma
	    if(fabs(itMC->getPID()) == 11 && fabs(itMC->getMomPID())== 24){EleLeg1=itMC;nEleLeg1+=1;}
	    if(fabs(itMC->getPID()) == 12 && fabs(itMC->getMomPID())== 24){NvLeg=itMC; nNvLeg+=1;}
	    //else if(itMC->getPID() == -11 && itMC->getMomPID()== 23){EleLeg2=itMC;nEleLeg2+=1;}
	    //else if(itMC->getPID() == 22 && itMC->getEt()>10 && fabs(itMC->getMomPID())== 11 && itMC->getGMomPID()== 23){PhoLeg = itMC; nPhoLeg+=1;}
	    else if(itMC->getPID() == 22 && itMC->getEt()>10){
              switch((int)fabs(itMC->getMomPID())){
		// 1-6 for quarks, 21:gluon, 24:W+, 11:electron
                case 1: PhoLeg = itMC; nPhoLeg+=1; break;
                case 2: PhoLeg = itMC; nPhoLeg+=1; break;
                case 3: PhoLeg = itMC; nPhoLeg+=1; break;
                case 4: PhoLeg = itMC; nPhoLeg+=1; break;
                case 5: PhoLeg = itMC; nPhoLeg+=1; break;
                case 6: PhoLeg = itMC; nPhoLeg+=1; break;
                case 21: PhoLeg = itMC; nPhoLeg+=1; break;
                case 24: PhoLeg = itMC; nPhoLeg+=1; break;
                case 11: PhoLeg = itMC; nPhoLeg+=1; break;
                default: break;
              }
            }
	  }
        
          double trueMET = sqrt(sum_ex*sum_ex + sum_ey*sum_ey);

          p_checkMETdiff->Fill(MET - trueMET); 
        
          //(((((((((((
          //Look for neutralino-decayed photon
          //if(nEleLeg1==1 && nEleLeg2==1 && nPhoLeg==1){
          if(nEleLeg1==1 && nPhoLeg==1){
             p_mcphotonET->Fill(PhoLeg->getEt());
             p_mcphotonEta->Fill(PhoLeg->getEta());                      
             // Match reco photon to signal photon
             float mindR(1);
             int nphomatch(0);
             std::vector<recoPhoton>::iterator recopho;
	     for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
               float tmpdE = fabs(itpho->getEt() - PhoLeg->getEt())/PhoLeg->getEt();
               if(tmpdE < 0.1){
                  float tmpdR = DeltaR(itpho->getEta(), itpho->getPhi(), PhoLeg->getEta(), PhoLeg->getPhi());
                  if(tmpdR < mindR){mindR = tmpdR;recopho = itpho;}
                  if(mindR < 0.15)nphomatch+=1;
               }
             }
             p_checkPhoMatch->Fill(nphomatch);

             if(mindR < 0.15){
	       p_phoMatchEff->Fill(PhoLeg->getEt(),1);

	       if(recopho->isEB()){
	         if(recopho->isLoose())p_phoIDeff->Fill(recopho->getEt(),1);
	         else p_phoIDeff->Fill(recopho->getEt(),0);
		 p_phoEBR9_true->Fill(recopho->getR9());
		 p_phoEBHoverE_true->Fill(recopho->getHoverE());
		 p_phoEBsigma_true->Fill(recopho->getSigma());
		 p_phoEBChIso_true->Fill(recopho->getChIso());
		 p_phoEBNeuIso_true->Fill(recopho->getNeuIso());
		 p_phoEBPhoIso_true->Fill(recopho->getPhoIso());
		 if(recopho->getEt()>40){
		   p_phoEBR9_ptselect->Fill(recopho->getR9());
		   p_phoEBHoverE_ptselect->Fill(recopho->getHoverE());
		   p_phoEBsigma_ptselect->Fill(recopho->getSigma());
		   p_phoEBChIso_ptselect->Fill(recopho->getChIso());
		   p_phoEBNeuIso_ptselect->Fill(recopho->getNeuIso());
		   p_phoEBPhoIso_ptselect->Fill(recopho->getPhoIso());
		   if(recopho->isLoose()){
		     p_phoEBR9_idselect->Fill(recopho->getR9());
		     p_phoEBHoverE_idselect->Fill(recopho->getHoverE());
		     p_phoEBsigma_idselect->Fill(recopho->getSigma());
		     p_phoEBChIso_idselect->Fill(recopho->getChIso());
		     p_phoEBNeuIso_idselect->Fill(recopho->getNeuIso());
		     p_phoEBPhoIso_idselect->Fill(recopho->getPhoIso());
		   } 
		 }
	       }
	       if(recopho->isEE()){
	         if(recopho->isLoose())p_phoIDeff->Fill(recopho->getEt(),1);
	         else p_phoIDeff->Fill(recopho->getEt(),0);
		 p_phoEER9_true->Fill(recopho->getR9());
		 p_phoEEHoverE_true->Fill(recopho->getHoverE());
		 p_phoEEsigma_true->Fill(recopho->getSigma());
		 p_phoEEChIso_true->Fill(recopho->getChIso());
		 p_phoEENeuIso_true->Fill(recopho->getNeuIso());
		 p_phoEEPhoIso_true->Fill(recopho->getPhoIso());
		 if(recopho->getEt()>40){
		   p_phoEER9_ptselect->Fill(recopho->getR9());
		   p_phoEEHoverE_ptselect->Fill(recopho->getHoverE());
		   p_phoEEsigma_ptselect->Fill(recopho->getSigma());
		   p_phoEEChIso_ptselect->Fill(recopho->getChIso());
		   p_phoEENeuIso_ptselect->Fill(recopho->getNeuIso());
		   p_phoEEPhoIso_ptselect->Fill(recopho->getPhoIso());
		   if(recopho->isLoose()){
		     p_phoEER9_idselect->Fill(recopho->getR9());
		     p_phoEEHoverE_idselect->Fill(recopho->getHoverE());
		     p_phoEEsigma_idselect->Fill(recopho->getSigma());
		     p_phoEEChIso_idselect->Fill(recopho->getChIso());
		     p_phoEENeuIso_idselect->Fill(recopho->getNeuIso());
		     p_phoEEPhoIso_idselect->Fill(recopho->getPhoIso());
		   } 
		 }
	       }
	      }//endif: dR<0.3 dE<0.1
	      else p_phoMatchEff->Fill(PhoLeg->getEt(),0);
              //End of signal photon
            //))))))))))))))))))))))
           

	     std::vector<mcData>::iterator signalEle;
             //if(EleLeg1->getEt() > EleLeg2->getEt())signalEle=EleLeg1;
             //else signalEle=EleLeg2;
             signalEle=EleLeg1;
             p_mceleET->Fill(signalEle->getEt());
             p_mceleEta->Fill(signalEle->getEta());                      
             // Match reco photon to signal photon
             float minEledR(1);
             int nelematch(0);
             std::vector<recoEle>::iterator recoele;
	     for(std::vector<recoEle>::iterator itele = Ele.begin() ; itele != Ele.end(); ++itele){
               float tmpdE = fabs(itele->getPt() - signalEle->getPt())/signalEle->getPt();
               if(tmpdE < 0.1){
                  float tmpdR = DeltaR(itele->getEta(), itele->getPhi(), signalEle->getEta(), signalEle->getPhi());
                  if(tmpdR < minEledR){minEledR = tmpdR;recoele = itele;}
                  if(minEledR < 0.15)nelematch+=1;
               }
             }
             p_checkEleMatch->Fill(nelematch);
	     if(minEledR < 0.15){

	       if(recoele->isEB()){
	         if(recoele->isMedium())p_eleIDeff->Fill(recoele->getEt(),1);
	         else p_eleIDeff->Fill(recoele->getEt(),0);
		 p_eleEBR9_true->Fill(recoele->getR9());
		 p_eleEBHoverE_true->Fill(recoele->getHoverE());
		 p_eleEBsigma_true->Fill(recoele->getSigma());
		 p_eleEBChIso_true->Fill(recoele->getChIso());
		 p_eleEBNeuIso_true->Fill(recoele->getNeuIso());
		 p_eleEBPhoIso_true->Fill(recoele->getPhoIso());
		 p_eleEBRelIso_true->Fill(recoele->getRelIso()); 
		 p_eleEBdEtaIn_true->Fill(recoele->getdEtaIn());
		 p_eleEBdPhiIn_true->Fill(recoele->getdPhiIn());
		 p_eleEBD0_true->Fill(recoele->getD0());
		 p_eleEBDz_true->Fill(recoele->getDz());
		 p_eleEBooEmooP_true->Fill(1/(recoele->getE())*(1-recoele->getEoverPInv()));
		 p_eleEBR9vsPt->Fill(recoele->getPt(), recoele->getR9());

		 if(recoele->getPt()>25){
		   p_eleConvVeto_ptselect->Fill(recoele->getConvVeto());
		   p_eleEBR9_ptselect->Fill(recoele->getR9());
		   p_eleEBHoverE_ptselect->Fill(recoele->getHoverE());
		   p_eleEBsigma_ptselect->Fill(recoele->getSigma());
		   p_eleEBChIso_ptselect->Fill(recoele->getChIso());
		   p_eleEBNeuIso_ptselect->Fill(recoele->getNeuIso());
		   p_eleEBPhoIso_ptselect->Fill(recoele->getPhoIso());
		   p_eleEBRelIso_ptselect->Fill(recoele->getRelIso()); 
		   p_eleEBdEtaIn_ptselect->Fill(recoele->getdEtaIn());
		   p_eleEBdPhiIn_ptselect->Fill(recoele->getdPhiIn());
		   p_eleEBD0_ptselect->Fill(recoele->getD0());
		   p_eleEBDz_ptselect->Fill(recoele->getDz());
		   p_eleEBooEmooP_ptselect->Fill(1/(recoele->getE())*(1-recoele->getEoverPInv()));
		   if(recoele->isMedium()){
		     p_eleEBR9_idselect->Fill(recoele->getR9());
		     p_eleEBHoverE_idselect->Fill(recoele->getHoverE());
		     p_eleEBsigma_idselect->Fill(recoele->getSigma());
		     p_eleEBChIso_idselect->Fill(recoele->getChIso());
		     p_eleEBNeuIso_idselect->Fill(recoele->getNeuIso());
		     p_eleEBPhoIso_idselect->Fill(recoele->getPhoIso());
		   }
		 }
	       }
	       if(recoele->isEE()){
	         if(recoele->isMedium())p_eleIDeff->Fill(recoele->getEt(),1);
	         else p_eleIDeff->Fill(recoele->getEt(),0);
		 p_eleEER9_true->Fill(recoele->getR9());
		 p_eleEEHoverE_true->Fill(recoele->getHoverE());
		 p_eleEEsigma_true->Fill(recoele->getSigma());
		 p_eleEEChIso_true->Fill(recoele->getChIso());
		 p_eleEENeuIso_true->Fill(recoele->getNeuIso());
		 p_eleEEPhoIso_true->Fill(recoele->getPhoIso());
		 p_eleEERelIso_true->Fill(recoele->getRelIso()); 
		 p_eleEEdEtaIn_true->Fill(recoele->getdEtaIn());
		 p_eleEEdPhiIn_true->Fill(recoele->getdPhiIn());
		 p_eleEED0_true->Fill(recoele->getD0());
		 p_eleEEDz_true->Fill(recoele->getDz());
		 p_eleEEooEmooP_true->Fill(1/(recoele->getE())*(1-recoele->getEoverPInv()));
		 p_eleEER9vsPt->Fill(recoele->getPt(), recoele->getR9());
		 if(recoele->getPt()>25){
		   p_eleConvVeto_ptselect->Fill(recoele->getConvVeto());
		   p_eleEER9_ptselect->Fill(recoele->getR9());
		   p_eleEEHoverE_ptselect->Fill(recoele->getHoverE());
		   p_eleEEsigma_ptselect->Fill(recoele->getSigma());
		   p_eleEEChIso_ptselect->Fill(recoele->getChIso());
		   p_eleEENeuIso_ptselect->Fill(recoele->getNeuIso());
		   p_eleEEPhoIso_ptselect->Fill(recoele->getPhoIso());
		   p_eleEERelIso_ptselect->Fill(recoele->getRelIso()); 
		   p_eleEEdEtaIn_ptselect->Fill(recoele->getdEtaIn());
		   p_eleEEdPhiIn_ptselect->Fill(recoele->getdPhiIn());
		   p_eleEED0_ptselect->Fill(recoele->getD0());
		   p_eleEEDz_ptselect->Fill(recoele->getDz());
		   p_eleEEooEmooP_ptselect->Fill(1/(recoele->getE())*(1-recoele->getEoverPInv()));
		   if(recoele->isMedium()){
		     p_eleEER9_idselect->Fill(recoele->getR9());
		     p_eleEEHoverE_idselect->Fill(recoele->getHoverE());
		     p_eleEEsigma_idselect->Fill(recoele->getSigma());
		     p_eleEEChIso_idselect->Fill(recoele->getChIso());
		     p_eleEENeuIso_idselect->Fill(recoele->getNeuIso());
		     p_eleEEPhoIso_idselect->Fill(recoele->getPhoIso());
		   }
		 }
	       };
             }//endif: dR<0.15 dE<0.1

             double dRPhoEle = DeltaR(PhoLeg->getEta(), PhoLeg->getPhi(), signalEle->getEta(), signalEle->getPhi());
             double massPhoEle = (PhoLeg->getP4()+signalEle->getP4()).M();
             double MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(DeltaR(0, signalEle->getPhi(), 0, METPhi))));
             MT_ = MT;
             p_PhoELeDeltaR->Fill(dRPhoEle);
             p_PhoEleMass->Fill(massPhoEle);
             p_PhoEleDRvsM->Fill(massPhoEle, dRPhoEle);
             p_Mt->Fill(MT);
             double trimass = sqrt(2*MET*(PhoLeg->getP4()+signalEle->getP4()).Pt()*(1-std::cos(DeltaR(0, (PhoLeg->getP4()+signalEle->getP4()).Phi(), 0, METPhi))));
             ThreeBodyMass_ = trimass; 
             p_threebodymass->Fill(trimass);
             if(fabs(PhoLeg->getMomPID())==24)p_threebodymass_W->Fill(trimass);
             p_threebodymass_MT->Fill(trimass, MT);
             if(PhoLeg->getEt() > 40 && signalEle->getPt() > 25 && dRPhoEle > 0.8)p_selectPhoEleMass->Fill(massPhoEle);
             tree->Fill(); 
            }//endif:: found signal ELectron
            //)))))))))))))))))))))))


    }//loop on entries


    outputfile->Write();
}
