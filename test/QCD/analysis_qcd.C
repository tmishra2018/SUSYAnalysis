// Here jet is taken as lepton only signalTree , signal events from QCD sample
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
#include "../../include/analysis_jet.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"


void analysis_qcd(){//main 

  gSystem->Load("../../lib/libAnaClasses.so");

  char outputname[100] = "/eos/uscms/store/user/tmishra/fakeLep/fakelep_QCD.root";
  ofstream logfile;
  logfile.open("/eos/uscms/store/user/tmishra/fakeLep/fakelep_egsignal_QCD.log"); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;
  //logfile << "Loose the proxy definition: no upper bounds for photon; LooseFakeProxy for electron" << std::endl;

  RunType datatype(MCDoubleEG2016); 
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("/eos/uscms/store/user/tmishra/InputFilesMC/QCD_DoubleEM/QCD_DoubleEM_2016.root");

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

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
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	float trailPt(0);
	float trailEta(0);
	float trailPhi(0);
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
  sigtree->Branch("sigMT",     &sigMT);
  sigtree->Branch("sigMET",    &sigMET);
  sigtree->Branch("sigMETPhi", &sigMETPhi);
  sigtree->Branch("dPhiLepMET",&dPhiLepMET);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
	sigtree->Branch("trailPt",   &trailPt);
	sigtree->Branch("trailEta",  &trailEta);
	sigtree->Branch("trailPhi",  &trailPhi);
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
	int   fakeLepIsStandardProxy(0);
	float fakeLepMiniIso(0);
	float fakeLepSigma(0);
	float fakeLepdEta(0);
	float fakeLepdPhi(0);
	unsigned fakeLepIndex(0);
  float fakeLepsigMT(0);
  float fakeLepsigMET(0);
  float fakeLepsigMETPhi(0);
  float fakeLepdPhiLepMET(0);
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  float fakeLepnJet(0);
	std::vector<float> fakeLep_JetPt;
	std::vector<float> fakeLep_JetEta;
	std::vector<float> fakeLep_JetPhi;
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
  fakeLeptree->Branch("fakeLepIsStandardProxy", &fakeLepIsStandardProxy);
  fakeLeptree->Branch("fakeLepMiniIso",&fakeLepMiniIso);
	fakeLeptree->Branch("fakeLepSigma",&fakeLepSigma);
	fakeLeptree->Branch("fakeLepdEta", &fakeLepdEta);
	fakeLeptree->Branch("fakeLepdPhi", &fakeLepdPhi);
	fakeLeptree->Branch("fakeLepIndex",&fakeLepIndex);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
	fakeLeptree->Branch("JetPt",     &fakeLep_JetPt);
	fakeLeptree->Branch("JetEta",    &fakeLep_JetEta);
	fakeLeptree->Branch("JetPhi",    &fakeLep_JetPhi);
	if(isMC){
  	fakeLeptree->Branch("mcPID",     &fakeLep_mcPID);
  	fakeLeptree->Branch("mcEta",     &fakeLep_mcEta);
  	fakeLeptree->Branch("mcPhi",     &fakeLep_mcPhi);
  	fakeLeptree->Branch("mcPt",      &fakeLep_mcPt);
  	fakeLeptree->Branch("mcMomPID",  &fakeLep_mcMomPID);
  	fakeLeptree->Branch("mcGMomPID", &fakeLep_mcGMomPID);
	}

//*********** histo list **********************//
  TH1F *p_eventcount = new TH1F("p_eventcount","p_eventcount",7,0,7);

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

			nTotal+=1;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;


		//	bool hasPho(false);
		//	std::vector<recoJet>::iterator signalPho = JetCollection.begin();
		//	for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
		//		if(!itJet->passSignalSelection())continue;
		//		if(itJet->getPt() < 35 || fabs(itJet->getEta())>1.442)continue;
		//		if(!hasPho){
		//			signalPho = itJet;
		//			hasPho = true;
		//		}
		//	}	

			bool hasLep(false);
			std::vector<recoJet>::iterator signalLep = JetCollection.begin();
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(!itJet->passSignalSelection())continue;
				if(itJet->getPt() < 25)continue;
			
				if(!hasLep){
					hasLep=true; 
					npassLep +=1;
					signalLep = itJet;
				}
			}

			if(hasLep){
				if(raw.passMETFilter(METFilter)){
					npassMETFilter +=1;

					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));

					lepPt = signalLep->getPt();
					lepEta= signalLep->getEta();
					lepPhi= signalLep->getPhi();
					sigMT = MT;
					sigMET= MET;
					sigMETPhi = METPhi;
					dPhiLepMET = deltaPhi; 
					nVertex = nVtx; 

					nJet = 0;
					HT = 0;
					for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
						if(!itJet->passSignalSelection())continue;
						if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
						nJet += 1;
						HT += itJet->getPt();
					}	

          mcPID.clear();
          mcEta.clear();
          mcPhi.clear();
          mcPt.clear();
          mcMomPID.clear();
          mcGMomPID.clear();
					if(isMC){
         	 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
         	   if(itMC->getEt() < 1.0)continue;
         	   mcPID.push_back(itMC->getPID());
         	   mcMomPID.push_back(itMC->getMomPID());
         	   mcGMomPID.push_back(itMC->getGMomPID());
         	   mcEta.push_back(itMC->getEta());
         	   mcPhi.push_back(itMC->getPhi());
         	   mcPt.push_back(itMC->getEt());
         	 }
					}
					sigtree->Fill();

				}//MET Filter
			}// ele + pho candidate
	
	}//loop on  events

  p_eventcount->GetXaxis()->SetBinLabel(1,"nTotal");
  p_eventcount->GetXaxis()->SetBinLabel(2,"npassHLT");
  p_eventcount->GetXaxis()->SetBinLabel(3,"npassPho");
  p_eventcount->GetXaxis()->SetBinLabel(4,"npassLep");
  p_eventcount->GetXaxis()->SetBinLabel(5,"npassdR");
  p_eventcount->GetXaxis()->SetBinLabel(6,"npassZ");
  p_eventcount->GetXaxis()->SetBinLabel(7,"npassMETFilter");
  p_eventcount->Fill(0.5, nTotal);
  p_eventcount->Fill(1.5, npassHLT);
  p_eventcount->Fill(2.5, npassPho);
  p_eventcount->Fill(3.5, npassLep);
  p_eventcount->Fill(4.5, npassdR);
  p_eventcount->Fill(5.5, npassZ);
  p_eventcount->Fill(6.5, npassMETFilter);

	outputfile->Write();
	outputfile->Close();
	logfile.close();
}


