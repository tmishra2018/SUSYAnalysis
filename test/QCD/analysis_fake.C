#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

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
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_tools.h"

int RunYear = 2016;
bool passEleSignalSelection(float eta, float sigma, float dEta, float dPhi, float HoverE, float EoverP, int MissHit, int ConvVeto, float miniIso ){
  bool passMiniMedium(true);
 
  if( fabs(eta) > 1.4442 && fabs(eta) < 1.56) {passMiniMedium=false; return passMiniMedium; }
  if(fabs(eta) < 1.4442){
		if(sigma > 0.00998) {passMiniMedium=false; return passMiniMedium; }
		if(fabs(dEta) > 0.00311){passMiniMedium=false; return passMiniMedium;}
		if(fabs(dPhi) > 0.103){passMiniMedium=false; return passMiniMedium;}
		if(HoverE > 0.253){passMiniMedium=false; return passMiniMedium;}
		if(fabs(EoverP) > 0.134){passMiniMedium=false; return passMiniMedium;}
		if(MissHit > 1){passMiniMedium=false; return passMiniMedium;}
		if(ConvVeto == 0){passMiniMedium=false; return passMiniMedium;}
		//if(miniIso > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  else if(fabs(eta) > 1.56){
		if(sigma > 0.0298) {passMiniMedium=false; return passMiniMedium; }
		if(fabs(dEta) > 0.00609){passMiniMedium=false; return passMiniMedium;}
		if(fabs(dPhi) > 0.045){passMiniMedium=false; return passMiniMedium;}
		if(HoverE > 0.0878){passMiniMedium=false; return passMiniMedium;}
		if(fabs(EoverP) > 0.13){passMiniMedium=false; return passMiniMedium;}
		if(MissHit > 1){passMiniMedium=false; return passMiniMedium;}
		if(ConvVeto == 0){passMiniMedium=false; return passMiniMedium;}
		//if(miniIso  > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  return passMiniMedium;
}

void analysis_fake(){//main 

  gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/eos/uscms/store/user/tmishra/fakeLep/test_egsignal_QCDfake_2.root";
  ofstream logfile;
  logfile.open("/eos/uscms/store/user/tmishra/fakeLep/test_egsignal_QCDfake.log"); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;
  //logfile << "Loose the proxy definition: no upper bounds for photon; LooseFakeProxy for electron" << std::endl;

  RunType datatype(MCDoubleEG2016); 
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add(Form("/eos/uscms/store/user/tmishra/InputFilesMC/QCD_DoubleEM/QCD_DoubleEM_%d.root",RunYear));
	//es->Add("/eos/uscms/store/user/mengleis/copied/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TrancheIV_v6-v1.root");

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

	TProfile *eff_pass = new TProfile("eff_pass","",40,0,400);
//************ Signal Tree **********************//
  TTree *sigtree = new TTree("signalTree","signalTree");
  int   run(0);
  Long64_t  event(0);
  int   lumis(0);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float lepMiniIso(0);
	float lepRelIso(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
  std::vector<int>   mcPID;
  std::vector<float> mcEta;
  std::vector<float> mcPhi;
  std::vector<float> mcPt;
  std::vector<int>   mcMomPID;
  std::vector<int>   mcGMomPID;

  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("lumis",     &lumis);
  sigtree->Branch("phoEt",     &phoEt);
  sigtree->Branch("phoEta",    &phoEta);
  sigtree->Branch("phoPhi",    &phoPhi);
  sigtree->Branch("lepPt",     &lepPt);
  sigtree->Branch("lepEta",    &lepEta);
  sigtree->Branch("lepPhi",    &lepPhi);
  sigtree->Branch("lepMiniIso",&lepMiniIso);
  sigtree->Branch("lepRelIso", &lepRelIso);
  sigtree->Branch("sigMT",     &sigMT);
  sigtree->Branch("sigMET",    &sigMET);
  sigtree->Branch("sigMETPhi", &sigMETPhi);
  sigtree->Branch("dPhiLepMET",&dPhiLepMET);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
	if(isMC){
  	sigtree->Branch("mcPID",     &mcPID);
  	sigtree->Branch("mcEta",     &mcEta);
  	sigtree->Branch("mcPhi",     &mcPhi);
  	sigtree->Branch("mcPt",      &mcPt);
  	sigtree->Branch("mcMomPID",  &mcMomPID);
  	sigtree->Branch("mcGMomPID", &mcGMomPID);
	}

//*********** fake lepton *********************//
  TTree *fakeLeptree = new TTree("fakeLepTree","fakeLepTree");
  float fakeLepphoEt(0);
  float fakeLepphoEta(0);
  float fakeLepphoPhi(0);
  float fakeLepPt(0);
  float fakeLepEta(0);
  float fakeLepPhi(0);
	float fakeLepMiniIso(0);
	float fakeLepSigma(0);
	float fakeLepdEta(0);
	float fakeLepdPhi(0);
	float	fakeLepRelIso(0);
	float	fakeLepD0(0);
	float	fakeLepDz(0);
	float	fakeLepHoverE(0);
	float	fakeLepEoverP(0);
	int		fakeLepVeto(0);
	int		fakeLepMissHits(0);
	float fakeLepR9(0);
	unsigned fakeLepIndex(0);
  float fakeLepsigMT(0);
  float fakeLepsigMET(0);
  float fakeLepsigMETPhi(0);
  float fakeLepdPhiLepMET(0);
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  float fakeLepnJet(0);
  std::vector<int>   fakeLep_mcPID;
  std::vector<float> fakeLep_mcEta;
  std::vector<float> fakeLep_mcPhi;
  std::vector<float> fakeLep_mcPt;
  std::vector<int>   fakeLep_mcMomPID;
  std::vector<int>   fakeLep_mcGMomPID;
  
  
  fakeLeptree->Branch("phoEt",     &fakeLepphoEt);
  fakeLeptree->Branch("phoEta",    &fakeLepphoEta);
  fakeLeptree->Branch("phoPhi",    &fakeLepphoPhi);
  fakeLeptree->Branch("lepPt",     &fakeLepPt);
  fakeLeptree->Branch("lepEta",    &fakeLepEta);
  fakeLeptree->Branch("lepPhi",    &fakeLepPhi);
  fakeLeptree->Branch("fakeLepMiniIso",&fakeLepMiniIso);
	fakeLeptree->Branch("fakeLepSigma",&fakeLepSigma);
	fakeLeptree->Branch("fakeLepdEta", &fakeLepdEta);
	fakeLeptree->Branch("fakeLepdPhi", &fakeLepdPhi);
	fakeLeptree->Branch("fakeLepRelIso",&fakeLepRelIso);
	fakeLeptree->Branch("fakeLepD0",  &fakeLepD0);
	fakeLeptree->Branch("fakeLepDz",  &fakeLepDz);
	fakeLeptree->Branch("fakeLepHoverE", &fakeLepHoverE);
	fakeLeptree->Branch("fakeLepEoverP", &fakeLepEoverP);
	fakeLeptree->Branch("fakeLepVeto",   &fakeLepVeto);
	fakeLeptree->Branch("fakeLepMissHits", &fakeLepMissHits);
	fakeLeptree->Branch("fakeLepR9",  &fakeLepR9);
	fakeLeptree->Branch("fakeLepIndex",&fakeLepIndex);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
	if(isMC){
  	fakeLeptree->Branch("mcPID",     &fakeLep_mcPID);
  	fakeLeptree->Branch("mcEta",     &fakeLep_mcEta);
  	fakeLeptree->Branch("mcPhi",     &fakeLep_mcPhi);
  	fakeLeptree->Branch("mcPt",      &fakeLep_mcPt);
  	fakeLeptree->Branch("mcMomPID",  &fakeLep_mcMomPID);
  	fakeLeptree->Branch("mcGMomPID", &fakeLep_mcGMomPID);
	}

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  float MET(0);
  float METPhi(0);
  int nVtx(0);
  int jetNumber(0);
  int METFilter(0);
  logfile << "RunType: " << datatype << std::endl;

  std::cout << "Total evetns : " << nEvts << std::endl;
  logfile << "Total evetns : " << nEvts << std::endl;
	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

		if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if (ievt%100000==0) logfile  << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			MCData.clear();
			Photon.clear();
			Muon.clear();
			Ele.clear();
			JetCollection.clear();
			if(isMC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			nVtx = raw.nVtx;
			run=raw.run;
			event=raw.event;
			lumis=raw.lumis;

			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Lepton                                              **********************************************/
			int eleOrder(0);
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				eleOrder += 1;
				if(itEle->getCalibPt() < 25)continue;
				std::vector<recoPhoton>::iterator matchPho = Photon.begin();
				bool hasmatch(false);
				for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), itEle->getEta(), itEle->getPhi()) < 0.05){
						matchPho = itpho;
						hasmatch = true; // ele match pho
					}
				}
				bool isSig = passEleSignalSelection( itEle->getEta(), itEle->getSigma(),fabs(itEle->getdEtaIn()) , fabs(itEle->getdPhiIn()), itEle->getHoverE(), itEle->getEoverPInv(), itEle->getMissHits(), itEle->getConvVeto(), itEle->getMiniIso());
				//if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))isSig=false;
				//if(itEle->passSignalSelection() && isSig){

				
				// electron signal selection
				if(isSig){
					if(eleOrder == 1)eff_pass->Fill(itEle->getCalibPt(), 1);
					phoEt = 100; 
					phoEta= 0.5; 
					phoPhi= 0.5; 
					lepPt = itEle->getCalibPt();
					lepEta= itEle->getEta();
					lepPhi= itEle->getPhi();
					lepMiniIso = itEle->getMiniIso();
					lepRelIso= itEle->getRelIso();
					sigMT = 100;
					sigMET= 100;
					sigMETPhi = 0.5;
					dPhiLepMET = 0.5; 
					nVertex = 1; 
					dRPhoLep= 0.8;
					if(hasmatch){
						if(matchPho->PixelSeed()==0)nJet = 0;
						else nJet = 1;
					}
					else nJet = -1;
					HT = 0;
          mcPID.clear();
          mcEta.clear();
          mcPhi.clear();
          mcPt.clear();
          mcMomPID.clear();
          mcGMomPID.clear();
					if(isMC){
           for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
             if(itMC->getEt() < 1.0)continue;
						 double dR = DeltaR( itMC->getEta(), itMC->getPhi(), itEle->getEta(), itEle->getPhi());
						 if(dR > 0.3)continue;
             mcPID.push_back(itMC->getPID());
             mcMomPID.push_back(itMC->getMomPID());
             mcGMomPID.push_back(itMC->getGMomPID());
             mcEta.push_back(itMC->getEta());
             mcPhi.push_back(itMC->getPhi());
             mcPt.push_back(itMC->getEt());
           }
					}
					sigtree->Fill();
				}
				//else if(!itEle->passSignalSelection()){

				// electron does n't pass signal selection
				else{
					if(eleOrder == 1)eff_pass->Fill(itEle->getCalibPt(), 0);
					fakeLepIndex = eleOrder;
					fakeLepphoEt = 100; 
 					fakeLepphoEta= 0.5; 
					fakeLepphoPhi= 0.5;
					fakeLepPt = itEle->getCalibPt();
					fakeLepEta= itEle->getEta();
					fakeLepPhi= itEle->getPhi();
					fakeLepMiniIso = itEle->getMiniIso();
					fakeLepSigma = itEle->getSigma();
					fakeLepdEta  = fabs(itEle->getdEtaIn());
					fakeLepdPhi  = fabs(itEle->getdPhiIn());
					fakeLepRelIso= itEle->getRelIso();
					fakeLepD0    = itEle->getD0();
					fakeLepDz    = itEle->getDz();
					fakeLepHoverE= itEle->getHoverE();
					fakeLepEoverP= itEle->getEoverPInv();   
					fakeLepVeto  = itEle->getConvVeto();
					fakeLepMissHits=itEle->getMissHits();
					fakeLepR9  = itEle->getR9();
					if(hasmatch){
						if(matchPho->PixelSeed()==0)fakeLepnJet = 0;
						else fakeLepnJet = 1;
					}
					else fakeLepnJet = -1;
					fakeLepHT = 0;
          fakeLep_mcPID.clear();
          fakeLep_mcEta.clear();
          fakeLep_mcPhi.clear();
          fakeLep_mcPt.clear();
          fakeLep_mcMomPID.clear();
          fakeLep_mcGMomPID.clear();
					if(isMC){
           for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
             if(itMC->getEt() < 1.0)continue;
						 double dR = DeltaR( itMC->getEta(), itMC->getPhi(), itEle->getEta(), itEle->getPhi());
						 if(dR > 0.3)continue;
             fakeLep_mcPID.push_back(itMC->getPID());
             fakeLep_mcMomPID.push_back(itMC->getMomPID());
             fakeLep_mcGMomPID.push_back(itMC->getGMomPID());
             fakeLep_mcEta.push_back(itMC->getEta());
             fakeLep_mcPhi.push_back(itMC->getPhi());
             fakeLep_mcPt.push_back(itMC->getEt());
           }
					}
					fakeLeptree->Fill();
				}
			}


	}

	outputfile->Write();
	outputfile->Close();
	logfile.close();
}


