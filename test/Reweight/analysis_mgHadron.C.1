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
#include "TFileCollection.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_mcData.h"

void analysis_mgHadron(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/plot_hadron_ISR.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/plot_hadron_ISR.log"); 

  logfile << "analysis_hadron()" << std::endl;

  RunType datatype(MuonEG2016); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_BCD.root");
  es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_EFG.root");
  es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_H.root");

  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

  int mcType = MCType::NOMC;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;
//************ Signal Tree **********************//
	TTree *signaltree = new TTree("signalTree","signalTree");
	float signal_phoEt(0);
	float signal_phoEta(0);
	float signal_phoPhi(0);
	float signal_phoSigma(0);
	float signal_phoChIso(0);
	float signal_sigMT(0);
	float signal_sigMET(0);
	float signal_sigMETPhi(0);
	float signal_dPhiLepMET(0);
	int   signal_nVertex(0);
	float signal_HT(0);
	float signal_nJet(0);
	std::vector<int> signal_mcPID;
	std::vector<float> signal_mcEta;
	std::vector<float> signal_mcPhi;
	std::vector<float> signal_mcPt;
	std::vector<int> signal_mcMomPID;
	float  signal_dilepMass(0);

	signaltree->Branch("phoEt",     &signal_phoEt);
	signaltree->Branch("phoEta",    &signal_phoEta);
	signaltree->Branch("phoPhi",    &signal_phoPhi);
	signaltree->Branch("phoSigma",  &signal_phoSigma);
	signaltree->Branch("phoChIso",  &signal_phoChIso);
	signaltree->Branch("sigMT",     &signal_sigMT);
	signaltree->Branch("sigMET",    &signal_sigMET);
	signaltree->Branch("sigMETPhi", &signal_sigMETPhi);
	signaltree->Branch("dPhiLepMET",&signal_dPhiLepMET);
	signaltree->Branch("nVertex",   &signal_nVertex);
	signaltree->Branch("HT",        &signal_HT);
	signaltree->Branch("nJet",      &signal_nJet);
	signaltree->Branch("mcPID",     &signal_mcPID);
	signaltree->Branch("mcEta",     &signal_mcEta);
	signaltree->Branch("mcPhi",     &signal_mcPhi);
	signaltree->Branch("mcPt",      &signal_mcPt);
	signaltree->Branch("mcMomPID",  &signal_mcMomPID);
	signaltree->Branch("mcType",    &mcType);
	signaltree->Branch("dilepMass", &signal_dilepMass);

	TTree *eleproxytree = new TTree("eleproxyTree","eleproxyTree");
	float eleproxy_phoEt(0);
	float eleproxy_phoEta(0);
	float eleproxy_phoPhi(0);
	float eleproxy_phoSigma(0);
	float eleproxy_phoChIso(0);
	float eleproxy_sigMT(0);
	float eleproxy_sigMET(0);
	float eleproxy_sigMETPhi(0);
	float eleproxy_dPhiLepMET(0);
	int   eleproxy_nVertex(0);
	float eleproxy_HT(0);
	float eleproxy_nJet(0);
	float eleproxy_dilepMass(0);

	eleproxytree->Branch("phoEt",     &eleproxy_phoEt);
	eleproxytree->Branch("phoEta",    &eleproxy_phoEta);
	eleproxytree->Branch("phoPhi",    &eleproxy_phoPhi);
	eleproxytree->Branch("phoSigma",  &eleproxy_phoSigma);
	eleproxytree->Branch("phoChIso",  &eleproxy_phoChIso);
	eleproxytree->Branch("sigMT",     &eleproxy_sigMT);
	eleproxytree->Branch("sigMET",    &eleproxy_sigMET);
	eleproxytree->Branch("sigMETPhi", &eleproxy_sigMETPhi);
	eleproxytree->Branch("dPhiLepMET",&eleproxy_dPhiLepMET);
	eleproxytree->Branch("nVertex",   &eleproxy_nVertex);
	eleproxytree->Branch("HT",        &eleproxy_HT);
	eleproxytree->Branch("nJet",      &eleproxy_nJet);
	eleproxytree->Branch("dilepMass", &eleproxy_dilepMass);

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
	float had_dilepMass(0);

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
	hadtree->Branch("dilepMass", &had_dilepMass);

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
		if(raw.nGoodVtx < 1)continue;

		int nGoodMu(0);
		for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
			if(im->isMedium() && im->getMiniIso() < 0.2 &&  im->getPt() > 15){nGoodMu+=1;}
		}

		if(nGoodMu != 2 || raw.nPho <1)continue;

		std::vector<recoPhoton>::iterator sigPho;
		bool hasSigPho(false);
		std::vector< std::vector<recoPhoton>::iterator> eleproxyVec;
		eleproxyVec.clear();
		std::vector< std::vector<recoPhoton>::iterator> hadPhoVec;
		hadPhoVec.clear();
 
    for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
      if(itpho->getR9() < 0.5)continue;
			if(!itpho->passHLTSelection())continue;
			if(!itpho->passBasicSelection())continue; 
      if(itpho->getCalibEt() < 25)continue;
      if(itpho->getChIso() > 20.0)continue;
      bool PixelVeto = itpho->PixelSeed()==0? true: false;
      bool GSFveto(true);
      bool FSRVeto(true);
      for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
        if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
        if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
			}
      for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
        if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
      if(GSFveto && PixelVeto && FSRVeto){
				if(!hasSigPho){
        	sigPho = itpho;
        	hasSigPho = true;
				}
				hadPhoVec.push_back(itpho);
      }
      if((!PixelVeto || !GSFveto))eleproxyVec.push_back(itpho);
    }

		bool hasMu(false);
		std::vector<recoMuon>::iterator signalMu = Muon.begin();
		for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
			if(itMu->getPt() < 25)continue;
			if(!itMu->passHLTSelection())continue;
			if(itMu->passSignalSelection()){
				if(!hasMu){
					hasMu=true; 
					signalMu = itMu;
				}
			}
		}
     
    if(hasSigPho && hasMu){
      double DeltaPhoLep = DeltaR(sigPho->getEta(), sigPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
      if(DeltaPhoLep > 0.8){ 
					bool foundZG(false); 
					std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == signalMu)continue;
						if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
						double llmass = (im->getP4()+signalMu->getP4()).M();
						double llgmass = (im->getP4()+signalMu->getP4()+sigPho->getP4()).M();
						if( (80 < llmass && llmass < 100) || ( 80 < llgmass && llgmass < 100)){
							foundZG = true;
							trailLep = im;
						}
					}	
					
					if(foundZG){	
						float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));
						signal_dilepMass = (trailLep->getP4()+signalMu->getP4()).M();
						signal_phoEt = sigPho->getCalibEt();
						signal_phoEta = sigPho->getEta();
						signal_phoPhi = sigPho->getPhi();
						signal_phoSigma = sigPho->getSigma();
						signal_phoChIso = sigPho->getChIso();
						signal_sigMT = MT;
						signal_sigMET = MET;
						signal_sigMETPhi = METPhi;
						signal_dPhiLepMET = deltaPhi;
						signal_nVertex = nVtx;
						signal_HT = 0;
						signal_nJet = jetNumber;
						signal_mcPID.clear();
						signal_mcEta.clear();
						signal_mcPhi.clear();
						signal_mcPt.clear();
						signal_mcMomPID.clear();
						if(datatype == MC){
							for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
								if(itMC->getEt() < 1.0)continue;
								float mcdR = DeltaR(sigPho->getEta(), sigPho->getPhi(), itMC->getEta(), itMC->getPhi());
								if(mcdR < 0.3){
									signal_mcPID.push_back(itMC->getPID());
									signal_mcMomPID.push_back(itMC->getMomPID());
									signal_mcEta.push_back(itMC->getEta());
									signal_mcPhi.push_back(itMC->getPhi());
									signal_mcPt.push_back(itMC->getEt());
								}
							}
						}
						signaltree->Fill();
					}
      }
    }


		if(hasMu){
    	for(unsigned ip(0); ip < eleproxyVec.size(); ip++){
        std::vector<recoPhoton>::iterator eleproxyPho= eleproxyVec[ip];
        double DeltaPhoLep = DeltaR(eleproxyPho->getEta(), eleproxyPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
        if(DeltaPhoLep > 0.8){
					bool foundZG(false); 
					std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == signalMu)continue;
						if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
						double llmass = (im->getP4()+signalMu->getP4()).M();
						double llgmass = (im->getP4()+signalMu->getP4()+eleproxyPho->getP4()).M();
						if( (80 < llmass && llmass < 100) || ( 80 < llgmass && llgmass < 100)){
							foundZG = true;
							trailLep = im;
						}
					}	
					
					if(foundZG){	
					eleproxy_dilepMass = (trailLep->getP4()+signalMu->getP4()).M();
          eleproxy_phoEt = eleproxyPho->getCalibEt();
          eleproxy_phoEta = eleproxyPho->getEta();
          eleproxy_phoPhi = eleproxyPho->getPhi();
          eleproxy_phoSigma = eleproxyPho->getSigma();
          eleproxy_phoChIso = eleproxyPho->getChIso();
          eleproxy_sigMT = 0;
          eleproxy_sigMET = MET;
          eleproxy_sigMETPhi = METPhi;
          eleproxy_dPhiLepMET = 0;
          eleproxy_nVertex = nVtx;
          eleproxy_HT = 0;
          eleproxy_nJet = jetNumber;
          eleproxytree->Fill();
					}
        }
    	}
		}

    if(hasMu){
      for(unsigned ip(0); ip < hadPhoVec.size(); ip++){
        std::vector<recoPhoton>::iterator hadPho= hadPhoVec[ip];
        double DeltaPhoLep = DeltaR(hadPho->getEta(), hadPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
        if(DeltaPhoLep > 0.8){
					bool foundZG(false); 
					std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == signalMu)continue;
						if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
						double llmass = (im->getP4()+signalMu->getP4()).M();
						double llgmass = (im->getP4()+signalMu->getP4()+hadPho->getP4()).M();
						if( (80 < llmass && llmass < 100) || ( 80 < llgmass && llgmass < 100)){
							foundZG = true;
							trailLep = im;
						}
					}	
					
					if(foundZG){	
					had_dilepMass = (trailLep->getP4()+signalMu->getP4()).M();
          had_phoEt = hadPho->getCalibEt();
          had_phoEta = hadPho->getEta();
          had_phoPhi = hadPho->getPhi();
          had_phoSigma = hadPho->getSigma();
          had_phoChIso = hadPho->getChIso();
          had_sigMT = 0;
          had_sigMET = MET;
          had_sigMETPhi = METPhi;
          had_dPhiLepMET = 0;
          had_nVertex = nVtx;
          had_HT = 0;
          had_nJet = jetNumber;
          hadtree->Fill();
					}
        }
      }
    }

 

  } 
outputfile->Write();
}


