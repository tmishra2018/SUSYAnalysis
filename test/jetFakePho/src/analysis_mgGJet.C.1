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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"

void analysis_mgGJet(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/Sep1/plot_hadron_mgGJet.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/Sep1/plot_hadron_mgGJet.log"); 

  logfile << "analysis_hadron()" << std::endl;

  RunType datatype(MC); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("root://cmseos.fnal.gov///store/group/lpcsusystealth/ggNtuple_leppho/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root");
  logfile << "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root" << std::endl;

  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

  int mcType = MCType::generalMC;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;
//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float mg_phoEt(0);
  float mg_phoEta(0);
  float mg_phoPhi(0);
  float mg_phoSigma(0);
  float mg_phoChIso(0);
  float mg_sigMT(0);
  float mg_sigMET(0);
  float mg_sigMETPhi(0);
  float mg_dPhiLepMET(0);
  int   mg_nVertex(0);
  float mg_HT(0);
  float mg_nJet(0);
  std::vector<int> mg_mcPID;
  std::vector<float> mg_mcEta;
  std::vector<float> mg_mcPhi;
  std::vector<float> mg_mcPt;
  std::vector<int> mg_mcMomPID;
  
  egtree->Branch("phoEt",     &mg_phoEt);
  egtree->Branch("phoEta",    &mg_phoEta);
  egtree->Branch("phoPhi",    &mg_phoPhi);
  egtree->Branch("phoSigma",  &mg_phoSigma);
  egtree->Branch("phoChIso",  &mg_phoChIso);
  egtree->Branch("sigMT",     &mg_sigMT);
  egtree->Branch("sigMET",    &mg_sigMET);
  egtree->Branch("sigMETPhi", &mg_sigMETPhi);
  egtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
  egtree->Branch("nVertex",   &mg_nVertex);
  egtree->Branch("HT",        &mg_HT);
  egtree->Branch("nJet",      &mg_nJet);
  egtree->Branch("mcPID",     &mg_mcPID);
  egtree->Branch("mcEta",     &mg_mcEta);
  egtree->Branch("mcPhi",     &mg_mcPhi);
  egtree->Branch("mcPt",      &mg_mcPt);
  egtree->Branch("mcMomPID",  &mg_mcMomPID);
  egtree->Branch("mcType",    &mcType);

  TTree *hadtree = new TTree("hadTree","hadTree");
  float had_phoEt(0);
  float had_phoEta(0);
  float had_phoPhi(0);
  float had_phoSigma(0);
  float had_phoChIso(0);
  float had_sigMT(0);
  float had_sigMET(0);
  float had_sigMETPhi(0);
  float had_dPhiLepMET(0);
  int   had_nVertex(0);
  float had_HT(0);
  float had_nJet(0);
  
  hadtree->Branch("phoEt",     &had_phoEt);
  hadtree->Branch("phoEta",    &had_phoEta);
  hadtree->Branch("phoPhi",    &had_phoPhi);
  hadtree->Branch("phoSigma",  &had_phoSigma);
  hadtree->Branch("phoChIso",  &had_phoChIso);
  hadtree->Branch("sigMT",     &had_sigMT);
  hadtree->Branch("sigMET",    &had_sigMET);
  hadtree->Branch("sigMETPhi", &had_sigMETPhi);
  hadtree->Branch("dPhiLepMET",&had_dPhiLepMET);
  hadtree->Branch("nVertex",   &had_nVertex);
  hadtree->Branch("HT",        &had_HT);
  hadtree->Branch("nJet",      &had_nJet);

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
  int nVtx(0);
  int jetNumber(0);
  int METFilter(0);
  logfile << "RunType: " << datatype << std::endl;

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "total " << nEvts << std::endl;
  logfile <<   "total " << nEvts << std::endl;

  for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
    if(ievt%10000==0)std::cout << " -- Processing event " << ievt << std::endl;
    if(ievt%10000==0)logfile << " -- Processing event " << ievt << std::endl;

	raw.GetData(es, ievt);
	MCData.clear();
	Photon.clear();
	Muon.clear();
	Ele.clear();
	if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
	for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
	for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
	for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
	MET = raw.pfMET;
	METPhi = raw.pfMETPhi;
	METFilter = raw.metFilters;
	nVtx = raw.nVtx;
	jetNumber = raw.nJet;

	if(raw.nPho <1)continue;
   
	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
    bool ismgCandidate(false), ishadCandidate(false);
    bool hasmgCandidate(false),  hashadCandidate(false);
		if(itpho->getR9() < 0.5)continue; 
	  if(itpho->getCalibEt() < 25 || !itpho->passHoverE(1)|| !itpho->passNeuIso(1) || !itpho->passPhoIso(1))continue;
    if(itpho->getChIso() > 20.0)continue;
	  bool PixelVeto = itpho->PixelSeed()==0? true: false;
	  bool GSFveto(true);
	  bool FSRVeto(true);
	  for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
		 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
	  }
	  for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
		 if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
	  if(GSFveto && PixelVeto && FSRVeto){
      ismgCandidate = true;
      ishadCandidate = true;
	  }
    if(!itpho->fireL1Trg(12) && !itpho->fireL1Trg(17))ismgCandidate=false;

    if(ismgCandidate && !hasmgCandidate){
			mg_phoEt = itpho->getCalibEt();
			mg_phoEta = itpho->getEta();
			mg_phoPhi = itpho->getPhi();
			mg_phoSigma = itpho->getSigma();
			mg_phoChIso = itpho->getChIso();
			//mg_sigMT = MT;
			mg_sigMET = MET;
			mg_sigMETPhi = METPhi;
			//mg_dPhiLepMET = deltaPhi;
			mg_nVertex = nVtx;
			mg_HT = 0;
			mg_nJet = jetNumber;
			mg_mcPID.clear();
			mg_mcEta.clear();
			mg_mcPhi.clear();
			mg_mcPt.clear();
			mg_mcMomPID.clear();
			if(datatype == MC){
				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					if(itMC->getEt() < 1.0)continue;
					float mcdR = DeltaR(itpho->getEta(), itpho->getPhi(), itMC->getEta(), itMC->getPhi());
					if(mcdR < 0.3){
						mg_mcPID.push_back(itMC->getPID());
						mg_mcMomPID.push_back(itMC->getMomPID());
						mg_mcEta.push_back(itMC->getEta());      
						mg_mcPhi.push_back(itMC->getPhi());
						mg_mcPt.push_back(itMC->getEt());
					}
		  	}
			}
			egtree->Fill();
			hasmgCandidate = true;
    }
    if(ishadCandidate && !hashadCandidate){
			had_phoEt = itpho->getCalibEt();
			had_phoEta = itpho->getEta();
			had_phoPhi = itpho->getPhi();
			had_phoSigma = itpho->getSigma();
			had_phoChIso = itpho->getChIso();
	//	had_sigMT = MT;
			had_sigMET = MET;
			had_sigMETPhi = METPhi;
//		had_dPhiLepMET = deltaPhi;
			had_nVertex = nVtx;
			had_HT = 0;
			had_nJet = jetNumber;
      hadtree->Fill();
      hashadCandidate = true;
    }
	}

  } 
outputfile->Write();
}


