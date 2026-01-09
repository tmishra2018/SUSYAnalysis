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
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"


void analysis_HadronMC(int RunYear, bool ISpreVFP){//main

        std::string whichVFP;
        if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	gSystem->Load("../../../lib/libAnaClasses.so");
	ofstream logfile;
	logfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/logs/plot_hadron_GJets_%d%s.log",RunYear,whichVFP.c_str()),ios::trunc);

	logfile << "analysis_hadron()" << std::endl;
	

	RunType datatype(MC); 

	TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add(Form("/eos/uscms/store/user/tmishra/InputFilesMC/GJets/GJets_DoubleEM_%d%s.root",RunYear,whichVFP.c_str()));

	TFile *outputfile;
	outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_GJet_%d%s.root",RunYear,whichVFP.c_str()),"RECREATE");
	outputfile->cd();

  //************ Signal Tree ********************** /
  TTree *egtree = new TTree("egTree","egTree");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_phoSigma(0);
  float eg_phoChIso(0);
  std::vector<int> eg_mcPID;
  std::vector<float> eg_mcEta;
  std::vector<float> eg_mcPhi;
  std::vector<float> eg_mcPt;
  std::vector<int> eg_mcMomPID;
  
  egtree->Branch("phoEt",     &eg_phoEt);
  egtree->Branch("phoEta",    &eg_phoEta);
  egtree->Branch("phoPhi",    &eg_phoPhi);
  egtree->Branch("phoSigma",  &eg_phoSigma);
  egtree->Branch("phoChIso",  &eg_phoChIso);
  egtree->Branch("mcPID",     &eg_mcPID);
  egtree->Branch("mcEta",     &eg_mcEta);
  egtree->Branch("mcPhi",     &eg_mcPhi);
  egtree->Branch("mcPt",      &eg_mcPt);
  egtree->Branch("mcMomPID",  &eg_mcMomPID);


  TTree *mgtree = new TTree("mgTree","mgTree");
  float mg_phoEt(0);
  float mg_phoEta(0);
  float mg_phoPhi(0);
  float mg_phoSigma(0);
  float mg_phoChIso(0);
  std::vector<int> mg_mcPID;
  std::vector<float> mg_mcEta;
  std::vector<float> mg_mcPhi;
  std::vector<float> mg_mcPt;
  std::vector<int> mg_mcMomPID;
  
  mgtree->Branch("phoEt",     &mg_phoEt);
  mgtree->Branch("phoEta",    &mg_phoEta);
  mgtree->Branch("phoPhi",    &mg_phoPhi);
  mgtree->Branch("phoSigma",  &mg_phoSigma);
  mgtree->Branch("phoChIso",  &mg_phoChIso);
  mgtree->Branch("mcPID",     &mg_mcPID);
  mgtree->Branch("mcEta",     &mg_mcEta);
  mgtree->Branch("mcPhi",     &mg_mcPhi);
  mgtree->Branch("mcPt",      &mg_mcPt);
  mgtree->Branch("mcMomPID",  &mg_mcMomPID);


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
  
		if(ievt%100000==0)std::cout << " -- Processing event " << ievt << std::endl;

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
		
		if(!raw.passHLT())continue;
		if(raw.nPho <1)continue;
		
		bool hasDiEM(false);
		if(Photon.size() >= 2){
			if(Photon[0].getCalibEt() > 25 && Photon[1].getCalibEt() > 25)hasDiEM = true;
		}

		//Get reco-level photons
		for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
			if(itpho->getR9() < 0.5)continue;
			if(itpho->getCalibEt() < 25)continue; 
			if(!itpho->passBasicSelection())continue;
			
			bool isegCandidate(false), ismgCandidate(false);
			bool hasegCandidate(false), hasmgCandidate(false);
			if(itpho->getChIso() > 20.0)continue;
			bool PixelVeto = itpho->PixelSeed()==0? true: false;
			bool GSFveto(true);
			bool FSRVeto(true);
			for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
			}
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
			if(GSFveto && PixelVeto && FSRVeto){
				isegCandidate = true;
				ismgCandidate = true;
			}
			if(!itpho->fireDoubleTrg(5) && !itpho->fireDoubleTrg(6))isegCandidate=false;
			if (RunYear == 2016) {
				    if (!itpho->fireL1Trg(12)) {
        				ismgCandidate = false;
    					}
			}else {
    				if (!itpho->fireL1Trg(32)) {
        				ismgCandidate = false;
    				}
			}

			if(isegCandidate && !hasegCandidate && hasDiEM){
				eg_phoEt = itpho->getCalibEt();
				eg_phoEta = itpho->getEta();
				eg_phoPhi = itpho->getPhi();
				eg_phoSigma = itpho->getSigma();
				eg_phoChIso = itpho->getChIso();
				eg_mcPID.clear();
				eg_mcEta.clear();
				eg_mcPhi.clear();
				eg_mcPt.clear();
				eg_mcMomPID.clear();
				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					if(itMC->getEt() < 1.0)continue;
					float mcdR = DeltaR(itpho->getEta(), itpho->getPhi(), itMC->getEta(), itMC->getPhi());
					if(mcdR < 0.3){
						eg_mcPID.push_back(itMC->getPID());
						eg_mcMomPID.push_back(itMC->getMomPID());
						eg_mcEta.push_back(itMC->getEta());      
						eg_mcPhi.push_back(itMC->getPhi());
						eg_mcPt.push_back(itMC->getEt());
					}
				}
				egtree->Fill();
				hasegCandidate = true;
			}
      			if(ismgCandidate && !hasmgCandidate && hasDiEM){
				mg_phoEt = itpho->getCalibEt();
				mg_phoEta = itpho->getEta();
				mg_phoPhi = itpho->getPhi();
				mg_phoSigma = itpho->getSigma();
				mg_phoChIso = itpho->getChIso();
				mg_mcPID.clear();
				mg_mcEta.clear();
				mg_mcPhi.clear();
				mg_mcPt.clear();
				mg_mcMomPID.clear();
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
				mgtree->Fill();
				hasmgCandidate = true;
			}
		}

	}// loop over events
	cout<<egtree->GetEntries()<<endl;
	cout<<mgtree->GetEntries()<<endl;
        logfile << "egTree: " << egtree->GetEntries() << "  percentage  " <<egtree->GetEntries()/nEvts<<std::endl;
        logfile << "mgTree: " << mgtree->GetEntries() << "  percentage  " <<mgtree->GetEntries()/nEvts<<std::endl;
	outputfile->Write();
}
int main(int argc, char** argv)
{
    bool ISpreVFP = (atoi(argv[2]) == 1);
    analysis_HadronMC(atoi(argv[1]), ISpreVFP);
    return 0;
}
