// g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_egMC.C -o analysis_egMC.exe `root-config --libs`

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
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"

int ichannel = 1;
void analysis_egMC(int RunYear, bool preVFP, const char *Sample){//main 

  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_%s_%d%s_.log",Sample,RunYear,whichVFP.c_str())); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;

  RunType datatype;
  if(RunYear==2016) datatype = MCDoubleEG2016;
  if(RunYear==2017) datatype = MCDoubleEG2017;
  if(RunYear==2018) datatype = MCDoubleEG2018;
  bool  isMC(true);

  TChain* es = new TChain("ggNtuplizer/EventTree");
	
  char* inputfile = new char[300];

  if (strstr(Sample, "DYJetsToLL") != NULL or strstr(Sample, "TTJets") != NULL)
        sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  else
        sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  es->Add(inputfile);

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;

  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile* outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_%s_%d%s.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
  outputfile->cd();
  
  int mcType;

  if(strstr(inputfile, "WGToLNuG") != NULL){
                std::cout << "WGToLNuG sample !" << std::endl;
                mcType = MCType::WGJetInclusive;
  }
  else if(strstr(inputfile, "WGJet40") != NULL){
                std::cout << "WGJet40 sample !" << std::endl;
                mcType = MCType::WGJet40;
  }
  else if(strstr(inputfile, "WGJet130") != NULL){
                std::cout << "WGJet130 sample !" << std::endl;
                mcType = MCType::WGJet130;
  }
  else if(strstr(inputfile, "ZGToLLG") != NULL){
                std::cout << "ZGInclusive sample !" << std::endl;
                mcType = MCType::ZGInclusive;
  }
  else if(strstr(inputfile, "DYJetsToLL") != NULL){
                std::cout << "DYJetsToLL sample !" << std::endl;
                mcType = MCType::DYLL50;
  }
  else  if(strstr(inputfile, "TTGJets") != NULL){
                std::cout << "TTGJets sample !" << std::endl;
                mcType = MCType::TTG;
  }
  else if(strstr(inputfile, "WWG") != NULL){
                std::cout << "WWG sample !" << std::endl;
                mcType = MCType::WWG;
  }
  else if(strstr(inputfile, "WZG") != NULL){
                std::cout << "WZG sample !" << std::endl;
                mcType = MCType::WZG;
  }
  else if(strstr(inputfile, "TTJets") != NULL){
                std::cout << "TTJets sample !" << std::endl;
                mcType = MCType::TT;
  }
  else if(strstr(inputfile, "WW") != NULL){
                std::cout << "WW sample !" << std::endl;
                mcType = MCType::WW;
  }
  else if(strstr(inputfile, "WZ") != NULL){
                std::cout << "WZ sample !" << std::endl;
                mcType = MCType::WZ;
  }


  else if(strstr(inputfile, "WJetsToLNu") != NULL){
                std::cout << "WJetsToLNu sample !" << std::endl;
		mcType = MCType::W;
  }
  else if(strstr(inputfile, "GJet") != NULL){
                std::cout << "GJet sample !" << std::endl;
                mcType = MCType::GJet;
  }
  else if(strstr(inputfile, "QCD_DoubleEM") != NULL){
                std::cout << "QCD_DoubleEM sample !" << std::endl;
                mcType = MCType::QCDEM40;
  }
  else {
                std::cout << "not specific MC !" << std::endl;
                mcType = MCType::NOMC;
  }

  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();
 
  Double_t  L1ECALPrefire;
  float lumiWeight = getEvtWeight(ichannel, RunYear, preVFP, crosssection, ntotalevent);

  cout<<"crosssection = "<<crosssection<<endl;
  cout<<"ntotalevent = "<<ntotalevent<<endl;
  cout<<"lumiWeight = "<<lumiWeight<<endl;

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

  sigtree->Branch("crosssection",&crosssection);
  sigtree->Branch("ntotalevent", &ntotalevent);
  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("L1ECALPrefire",     &L1ECALPrefire);
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
                                                                                                
//************ Signal Tree **********************//
  TTree *proxytree = new TTree("proxyTree","proxyTree");
  float proxyphoEt(0);
  float proxyphoEta(0);
  float proxyphoPhi(0);
  float proxylepPt(0);
  float proxylepEta(0);
  float proxylepPhi(0);
  float proxysigMT(0);
  float proxysigMET(0);
  float proxysigMETPhi(0);
  float proxydPhiLepMET(0);
	float proxythreeMass(0);
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  float proxynJet(0);
	int   proxyFSRVeto(0);
	float proxytrailPt(0);
	float proxytrailEta(0);
	float proxytrailPhi(0);
  
  proxytree->Branch("crosssection",&crosssection);
  proxytree->Branch("L1ECALPrefire",     &L1ECALPrefire);
  proxytree->Branch("ntotalevent", &ntotalevent);
  proxytree->Branch("phoEt",     &proxyphoEt);
  proxytree->Branch("phoEta",    &proxyphoEta);
  proxytree->Branch("phoPhi",    &proxyphoPhi);
  proxytree->Branch("lepPt",     &proxylepPt);
  proxytree->Branch("lepEta",    &proxylepEta);
  proxytree->Branch("lepPhi",    &proxylepPhi);
  proxytree->Branch("sigMT",     &proxysigMT);
  proxytree->Branch("sigMET",    &proxysigMET);
  proxytree->Branch("sigMETPhi", &proxysigMETPhi);
  proxytree->Branch("dPhiLepMET",&proxydPhiLepMET);
	proxytree->Branch("threeMass", &proxythreeMass);
  proxytree->Branch("nVertex",   &proxynVertex);
  proxytree->Branch("dRPhoLep",  &proxydRPhoLep);
  proxytree->Branch("HT",        &proxyHT);
  proxytree->Branch("nJet",      &proxynJet);
	proxytree->Branch("FSRVeto",   &proxyFSRVeto);
	proxytree->Branch("trailPt",   &proxytrailPt);
	proxytree->Branch("trailEta",  &proxytrailEta);
	proxytree->Branch("trailPhi",  &proxytrailPhi);

//************ Signal Tree **********************//
  TTree *jettree = new TTree("jetTree","jetTree");
  float jetphoEt(0);
  float jetphoEta(0);
  float jetphoPhi(0);
  float jetlepPt(0);
  float jetlepEta(0);
  float jetlepPhi(0);
  float jetsigMT(0);
  float jetsigMET(0);
  float jetsigMETPhi(0);
  float jetdPhiLepMET(0);
  int   jetnVertex(0);
  float jetdRPhoLep(0);
  float jetHT(0);
  float jetnJet(0);
	float jettrailPt(0);
	float jettrailEta(0);
	float jettrailPhi(0);
  
  jettree->Branch("crosssection",&crosssection);
  jettree->Branch("ntotalevent", &ntotalevent);
  jettree->Branch("L1ECALPrefire",     &L1ECALPrefire);
  jettree->Branch("phoEt",     &jetphoEt);
  jettree->Branch("phoEta",    &jetphoEta);
  jettree->Branch("phoPhi",    &jetphoPhi);
  jettree->Branch("lepPt",     &jetlepPt);
  jettree->Branch("lepEta",    &jetlepEta);
  jettree->Branch("lepPhi",    &jetlepPhi);
  jettree->Branch("sigMT",     &jetsigMT);
  jettree->Branch("sigMET",    &jetsigMET);
  jettree->Branch("sigMETPhi", &jetsigMETPhi);
  jettree->Branch("dPhiLepMET",&jetdPhiLepMET);
  jettree->Branch("nVertex",   &jetnVertex);
  jettree->Branch("dRPhoLep",  &jetdRPhoLep);
  jettree->Branch("HT",        &jetHT);
  jettree->Branch("nJet",      &jetnJet);
	jettree->Branch("trailPt",   &jettrailPt);
	jettree->Branch("trailEta",  &jettrailEta);
	jettree->Branch("trailPhi",  &jettrailPhi);
  
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
	unsigned fakeLepIndex(0);
  float fakeLepsigMT(0);
  float fakeLepsigMET(0);
  float fakeLepsigMETPhi(0);
  float fakeLepdPhiLepMET(0);
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  float fakeLepnJet(0);
  
  
  fakeLeptree->Branch("crosssection",&crosssection);
  fakeLeptree->Branch("ntotalevent", &ntotalevent);
  fakeLeptree->Branch("L1ECALPrefire",     &L1ECALPrefire);
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
	fakeLeptree->Branch("fakeLepIndex",&fakeLepIndex);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);

//*************** for jet-photon fake rate ***********************//
	TTree *hadrontree = new TTree("hadronTree","hadronTree");
	float hadron_phoEt(0);
	float hadron_phoEta(0);
	float hadron_phoPhi(0);
	float hadron_phoSigma(0);
	float hadron_phoChIso(0);
  std::vector<float> hadron_eleproxyEt;
	std::vector<float> hadron_eleproxyEta;
	std::vector<float> hadron_eleproxyPhi;
	std::vector<float> hadron_eleproxySigma;
	std::vector<float> hadron_eleproxyChIso;
	std::vector<int> hadron_eleproxynVertex;
	float hadron_sigMT(0);
	float hadron_sigMET(0);
	float hadron_sigMETPhi(0);
	float hadron_dPhiLepMET(0);
	int   hadron_nVertex(0);
	float hadron_HT(0);
	float hadron_nJet(0);
	std::vector<int> hadron_mcPID;
	std::vector<float> hadron_mcEta;
	std::vector<float> hadron_mcPhi;
	std::vector<float> hadron_mcPt;
	std::vector<int> hadron_mcMomPID;

	hadrontree->Branch("crosssection",&crosssection);
	hadrontree->Branch("ntotalevent", &ntotalevent);
  	hadrontree->Branch("L1ECALPrefire",     &L1ECALPrefire);
	hadrontree->Branch("phoEt",     &hadron_phoEt);
	hadrontree->Branch("phoEta",    &hadron_phoEta);
	hadrontree->Branch("phoPhi",    &hadron_phoPhi);
	hadrontree->Branch("phoSigma",  &hadron_phoSigma);
	hadrontree->Branch("phoChIso",  &hadron_phoChIso);
	hadrontree->Branch("eleproxyEt",&hadron_eleproxyEt);
	hadrontree->Branch("eleproxyEta",&hadron_eleproxyEta);
	hadrontree->Branch("eleproxyPhi",&hadron_eleproxyPhi);
	hadrontree->Branch("eleproxySigma", &hadron_eleproxySigma);
	hadrontree->Branch("eleproxyChIso", &hadron_eleproxyChIso);
	hadrontree->Branch("eleproxynVertex",&hadron_eleproxynVertex);
	hadrontree->Branch("sigMT",     &hadron_sigMT);
	hadrontree->Branch("sigMET",    &hadron_sigMET);
	hadrontree->Branch("sigMETPhi", &hadron_sigMETPhi);
	hadrontree->Branch("dPhiLepMET",&hadron_dPhiLepMET);
	hadrontree->Branch("nVertex",   &hadron_nVertex);
	hadrontree->Branch("HT",        &hadron_HT);
	hadrontree->Branch("nJet",      &hadron_nJet);
	if(isMC){
		hadrontree->Branch("mcPID",     &hadron_mcPID);
		hadrontree->Branch("mcEta",     &hadron_mcEta);
		hadrontree->Branch("mcPhi",     &hadron_mcPhi);
		hadrontree->Branch("mcPt",      &hadron_mcPt);
		hadrontree->Branch("mcMomPID",  &hadron_mcMomPID);
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

  std::cout << "Total events : " << nEvts << std::endl;
  logfile << "Total events : " << nEvts << std::endl;
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
			L1ECALPrefire=raw.L1ECALPrefire;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nEle < 1 || raw.nPho <1)continue;


			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Photon                                              **********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			std::vector< int > proxyPhoFSRVeto;
			proxyPhoCollection.clear();
			proxyPhoFSRVeto.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			// ******* for jet-photon fake rate***************//
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  hadeleproxyPhoCollection;
			hadeleproxyPhoCollection.clear();

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool photonFSRVeto(true);
				bool eleFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getEt()>2.0)photonFSRVeto=false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.02)eleFSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0){
						photonFSRVeto=false;
						eleFSRVeto=false;
					}
				}
				// ******** very loose, before sigma and isolation cut.  For jet-photon fake rate, and hadron proxy ************//	
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasHadronPho){
						hasHadronPho = true;
						hadronPho = itpho;
					}

					if(!passSigma || !passChIso){
						if( (itpho->getSigma()< 0.02 && itpho->isEB()) || (itpho->getSigma()< 0.04 && itpho->isEE()) )jetPhoCollection.push_back(itpho);
					}
				}
				else if( (!GSFveto || !PixelVeto) )hadeleproxyPhoCollection.push_back(itpho);

        // ****************  standard ID ************************************//
				if(!itpho->passSignalSelection())continue;
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

				if((!PixelVeto || !GSFveto)){
						proxyPhoCollection.push_back(itpho);
						proxyPhoFSRVeto.push_back( eleFSRVeto? 1:0);
				}
			}





			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Lepton                                              **********************************************/
			bool hasLep(false);
			std::vector<recoEle>::iterator signalLep = Ele.begin();
			bool hasTrail(false);
			std::vector<recoEle>::iterator trailLep = Ele.begin();
			std::vector< std::vector<recoEle>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoEle>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			std::vector< std::vector<recoEle>::iterator > miniisoLep;
			miniisoLep.clear();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if(itEle->isMedium() && itEle->getPt() > 15 && itEle->getMiniIso() < 0.2)miniisoLep.push_back(itEle);
				if(itEle->getCalibPt() < 25)continue;
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;

				if(!itEle->passHLTSelection())continue;
				if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
				//if(itEle->isLooseFakeProxy())fakeLepCollection.push_back(itEle);//Loose the proxy definition	
				if(itEle->passSignalSelection()){
					proxyLepCollection.push_back(itEle);
					if(hasLep && !hasTrail){
						hasTrail = true;
						trailLep = itEle;
					}
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itEle;
					}
				}
			}

			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi()); 
				if(dRlepphoton > 0.8){
					npassdR+=1;
					if(fabs((signalPho->getCalibP4()+signalLep->getCalibP4()).M() - 91.188) > 10.0){

						npassZ+=1;
						if(raw.passMETFilter(RunYear,METFilter)){
							npassMETFilter +=1;

							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getCalibPt()*(1-std::cos(deltaPhi)));

							phoEt = signalPho->getCalibEt();
							phoEta= signalPho->getEta();
							phoPhi= signalPho->getPhi();
							lepPt = signalLep->getCalibPt();
							lepEta= signalLep->getEta();
							lepPhi= signalLep->getPhi();
							sigMT = MT;
							sigMET= MET;
							sigMETPhi = METPhi;
							dPhiLepMET = deltaPhi; 
							nVertex = nVtx; 
							dRPhoLep= dRlepphoton;
							if(hasTrail){
								trailPt = trailLep->getPt();
								trailEta = trailLep->getEta();
								trailPhi = trailLep->getPhi();
							}
							else{
								trailPt  = 0; 
								trailEta = 0;
								trailPhi = 0;
							}

							nJet = 0;
							HT = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
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
          	 	   float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
          	 	   if(mcdR < 0.3){
          	 	     mcPID.push_back(itMC->getPID());
          	 	     mcMomPID.push_back(itMC->getMomPID());
          	 	     mcGMomPID.push_back(itMC->getGMomPID());
          	 	     mcEta.push_back(itMC->getEta());
          	 	     mcPhi.push_back(itMC->getPhi());
          	 	     mcPt.push_back(itMC->getEt());
          	 	   }
          	 	 }
							}
							sigtree->Fill();

						}//MET Filter
					}// Z mass Filter
				}//dR filter
			}// ele + pho candidate
	 
			if(!hasPho){
			for(unsigned ip(0); ip < proxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = proxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyLepCollection[ie];
					double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dRlepphoton>0.8){
						if(fabs((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(RunYear,METFilter)){

								float proxy_deltaPhi = DeltaPhi(proxyEle->getPhi(), METPhi);
								float proxy_MT = sqrt(2*MET*proxyEle->getCalibPt()*(1-std::cos(proxy_deltaPhi)));
								proxyphoEt = proxyPho->getCalibEt();
								proxyphoEta= proxyPho->getEta();
								proxyphoPhi= proxyPho->getPhi();
								proxylepPt = proxyEle->getCalibPt();
								proxylepEta= proxyEle->getEta();
								proxylepPhi= proxyEle->getPhi();
								proxysigMT = proxy_MT;
								proxysigMET= MET;
								proxysigMETPhi = METPhi;
								proxydPhiLepMET = proxy_deltaPhi; 
								proxynVertex = nVtx; 
								proxydRPhoLep= dRlepphoton;
								proxyFSRVeto = proxyPhoFSRVeto[ip];
								if(hasTrail){
									proxytrailPt = trailLep->getPt();
									proxytrailEta = trailLep->getEta();
									proxytrailPhi = trailLep->getPhi();
								}
								else{
									proxytrailPt  = 0; 
									proxytrailEta = 0;
									proxytrailPhi = 0;
								}

								proxynJet = 0;
								proxyHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyEle->getEta(),proxyEle->getPhi()) <= 0.4)continue;
									proxynJet += 1;
									proxyHT += itJet->getPt();
								}
								proxytree->Fill();

	
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection
			}
			
			if(hasLep){
				std::vector<recoEle>::iterator jetEle = signalLep; 
				for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
					std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
					double dRlepphoton = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetEle->getEta(), jetEle->getPhi());
					if(dRlepphoton>0.8){
						if(fabs((jetPho->getCalibP4()+jetEle->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(RunYear,METFilter)){


								float jet_deltaPhi = DeltaPhi(jetEle->getPhi(), METPhi);
								float jet_MT = sqrt(2*MET*jetEle->getCalibPt()*(1-std::cos(jet_deltaPhi)));
								jetphoEt = jetPho->getCalibEt();
								jetphoEta= jetPho->getEta();
								jetphoPhi= jetPho->getPhi();
								jetlepPt = jetEle->getCalibPt();
								jetlepEta= jetEle->getEta();
								jetlepPhi= jetEle->getPhi();
								jetsigMT = jet_MT;
								jetsigMET= MET;
								jetsigMETPhi = METPhi;
								jetdPhiLepMET = jet_deltaPhi; 
								jetnVertex = nVtx; 
								jetdRPhoLep= dRlepphoton;
								if(hasTrail){
									jettrailPt = trailLep->getPt();
									jettrailEta = trailLep->getEta();
									jettrailPhi = trailLep->getPhi();
								}
								else{
									jettrailPt  = 0; 
									jettrailEta = 0;
									jettrailPhi = 0;
								}

								jetnJet = 0;
								jetHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetEle->getEta(),jetEle->getPhi()) <= 0.4)continue;
									jetnJet += 1;
									jetHT += itJet->getPt();	
								}
								jettree->Fill();// jet fake photon tree

							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}

			if(hasPho && !hasLep){
				std::vector<recoPhoton>::iterator fakeLepPho = signalPho;
				for(unsigned ip(0); ip < fakeLepCollection.size(); ip++){
					std::vector<recoEle>::iterator fakeLep = fakeLepCollection[ip];
					double dRlepphoton = DeltaR(fakeLepPho->getEta(), fakeLepPho->getPhi(), fakeLep->getEta(), fakeLep->getPhi());
					if(dRlepphoton>0.8){
						if(fabs((fakeLepPho->getCalibP4()+fakeLep->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(RunYear,METFilter)){

								fakeLepdRPhoLep = 3;
								for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
									double dR = DeltaR(itpho->getEta(), itpho->getPhi(), fakeLep->getEta(), fakeLep->getPhi());
									if(dR < 0.05 && fabs(itpho->getEt() - fakeLep->getPt())/fakeLep->getPt() < 0.1)continue;
									if(dR < fakeLepdRPhoLep)fakeLepdRPhoLep = dR;
								}	
		
								float fakeLep_deltaPhi = DeltaPhi(fakeLep->getPhi(), METPhi);
								float fakeLep_MT = sqrt(2*MET*fakeLep->getCalibPt()*(1-std::cos(fakeLep_deltaPhi)));
								fakeLepphoEt = fakeLepPho->getCalibEt();
								fakeLepphoEta= fakeLepPho->getEta();
								fakeLepphoPhi= fakeLepPho->getPhi();
								fakeLepPt = fakeLep->getCalibPt();
								fakeLepEta= fakeLep->getEta();
								fakeLepPhi= fakeLep->getPhi();
								fakeLepMiniIso = fakeLep->getMiniIso();
								fakeLepsigMT = fakeLep_MT;
								fakeLepsigMET= MET;
								fakeLepsigMETPhi = METPhi;
								fakeLepdPhiLepMET = fakeLep_deltaPhi; 
								fakeLepnVertex = nVtx; 
								//fakeLepdRPhoLep= dRlepphoton;
								fakeLepSigma = fakeLep->getSigma();
								fakeLepdEta  = fabs(fakeLep->getdEtaIn());
								fakeLepdPhi  = fabs(fakeLep->getdPhiIn());
								fakeLepIndex = ip;

								fakeLepnJet = 0;
								fakeLepHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLep->getEta(),  fakeLep->getPhi()) <= 0.4)continue;
									fakeLepnJet += 1;
									fakeLepHT += itJet->getPt();
								}	
								fakeLeptree->Fill();// fake lepton tree

							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}


		if(hasHadronPho || hadeleproxyPhoCollection.size() > 0){
			hadron_phoEt = 0;
			hadron_phoEta = 0;
			hadron_phoPhi = 0;
			hadron_phoSigma = 0;
			hadron_phoChIso = 0;
			hadron_sigMT = 0;
			hadron_sigMET = 0;
			hadron_sigMETPhi = 0;
			hadron_dPhiLepMET = 0;
			hadron_nVertex = 0;
			hadron_HT = 0;
			hadron_nJet = 0;
  		hadron_eleproxyEt.clear();
			hadron_eleproxyEta.clear();
			hadron_eleproxyPhi.clear();
			hadron_eleproxySigma.clear();
			hadron_eleproxyChIso.clear();
			hadron_eleproxynVertex.clear();

			if(hasLep){
				double DeltaPhoLep = DeltaR(hadronPho->getEta(), hadronPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				double DoubleMass  = (hadronPho->getP4()+signalLep->getP4()).M();
				if(DeltaPhoLep > 0.8 && DoubleMass > 90 && fabs(DoubleMass - 91.188) > 10.0){
				if(raw.passMETFilter(RunYear,METFilter)){
					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
					hadron_phoEt = hadronPho->getCalibEt();
					hadron_phoEta = hadronPho->getEta();
					hadron_phoPhi = hadronPho->getPhi();
					hadron_phoSigma = hadronPho->getSigma();
					hadron_phoChIso = hadronPho->getChIso();
					hadron_sigMT = MT;
					hadron_sigMET = MET;
					hadron_sigMETPhi = METPhi;
					hadron_dPhiLepMET = deltaPhi;
					hadron_nVertex = nVtx;
					hadron_HT = 0;
					hadron_nJet = jetNumber;
					// photon passing all selection, except signal selection; hadron photon
				}}
			}
			for(unsigned ip(0); ip < hadeleproxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = hadeleproxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyLepCollection[ie];
					double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dRlepphoton>0.8 && fabs((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0 && raw.passMETFilter(RunYear,METFilter)){
						hadron_eleproxyEt.push_back(proxyPho->getCalibEt());
						hadron_eleproxyEta.push_back(proxyPho->getEta());
						hadron_eleproxyPhi.push_back(proxyPho->getPhi());
						hadron_eleproxySigma.push_back(proxyPho->getSigma());
						hadron_eleproxyChIso.push_back(proxyPho->getChIso());
						hadron_eleproxynVertex.push_back(nVtx); 
						// electron-fake-photon
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection

      hadron_mcPID.clear();
      hadron_mcEta.clear();
      hadron_mcPhi.clear();
      hadron_mcPt.clear();
      hadron_mcMomPID.clear();
			if(isMC){
       for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
         if(itMC->getPt() < 10.0)continue;
         mcPID.push_back(itMC->getPID());
         mcMomPID.push_back(itMC->getMomPID());
         mcEta.push_back(itMC->getEta());
         mcPhi.push_back(itMC->getPhi());
         mcPt.push_back(itMC->getEt());
       }
			}

			hadrontree->Fill();
		}
	
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
cout<< "signal tree entries "<< sigtree->GetEntries()<<endl;
cout<< "proxy tree entries "<< proxytree->GetEntries()<<endl;


	outputfile->Write();
	outputfile->Close();
	logfile.close();
}

int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    bool preVFP = (atoi(argv[2]) == 1);
    analysis_egMC(atoi(argv[1]), preVFP, argv[3]);
    return 0;
}
