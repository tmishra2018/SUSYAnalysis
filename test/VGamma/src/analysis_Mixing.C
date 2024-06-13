// look entries with photon(all, no pixel, id, eg and mg) cases and store photon, eg and mg trees
// g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_Mixing.C -o analysis_Mixing.exe `root-config --libs`
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
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
#include "TProfile2D.h"
#include "TCut.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"

bool nisrMatch(float jetEta, float jetPhi, std::vector<mcData>& genParticles){

    bool matched=false;
    for(std::vector<mcData>::iterator itMC = genParticles.begin(); itMC!= genParticles.end(); itMC++){
      	if(matched)break;
      	bool isPromptFinal(false);
      	int momid = abs(itMC->getMomPID());
      	if(abs(itMC->getPID())<=5 && (momid==6 || momid==23 || momid==24 || momid==25))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 11 && (momid<=6 || momid==23 || momid==24 || momid== 22))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 13 && (momid<=6 || momid==23 || momid==24 || momid== 22))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 15 && (momid<=6 || momid==23 || momid==24 || momid== 22))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 22 && (momid == 11 || momid == 13 || momid == 15 || momid==23 || momid==24))isPromptFinal = true;
			if(!isPromptFinal)continue;
			if(DeltaR(jetEta, jetPhi, itMC->getEta(), itMC->getPhi()) < 0.3){
				//std::cout << "match " << itMC->getPID() << " " << itMC->getMomPID() << " status:" << itMC->getStatus() << std::endl; 
      				matched = true;
			}
    } // Loop over MC particles
    return matched;
} 
		
bool preVFP = true;
void analysis_Mixing(int RunYear, const char *Sample){//main 
	
	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

  	RunType datatype(MC); 
  	TChain* es = new TChain("ggNtuplizer/EventTree");
	char* inputfile = new char[300];
        sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
        es->Add(inputfile);

  	const unsigned nEvts = es->GetEntries(); 
  	TFile* outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/egMC/mixing_%s_%d%s_TH1D.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
  	outputfile->cd();

	TH1D *all_phoEt = new TH1D("all_phoEt","",400,0,400);
	TH1D *nopixel_phoEt = new TH1D("nopixel_phoEt","",400,0,400);
	TH1D *id_phoEt = new TH1D("id_phoEt","",400,0,400);
	TH1D *egamma_phoEt = new TH1D("egamma_phoEt","",400,0,400);
  	TH1D *mugamma_phoEt = new TH1D("mugamma_phoEt","",400,0,400);

	double all_totalEvent(0), all_reweight(0);
	double nopixel_totalEvent(0), nopixel_reweight(0);
	double id_totalEvent(0), id_reweight(0);
	double egamma_totalEvent(0), egamma_reweight(0);
	double mugamma_totalEvent(0), mugamma_reweight(0);

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

  	if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 

  	float crosssection = MC_XS[mcType];
  	float ntotalevent = es->GetEntries();
	float PUweight(1);

	//************ Signal Tree **********************//
  	TTree *photontree = new TTree("phoTree","phoTree");
  	std::vector<float> pho_phoEt(0);
  	std::vector<float> pho_phoEta(0);
  	std::vector<float> pho_phoPhi(0);
  	std::vector<int>   pho_phoID(0);
  	std::vector<int>   pho_phoveto(0);
  	std::vector<int> pho_mcPID;
  	std::vector<float> pho_mcEta;
  	std::vector<float> pho_mcPhi;
  	std::vector<float> pho_mcPt;
  	std::vector<int> pho_mcMomPID;
  	std::vector<int> pho_mcGMomPID;
  	std::vector<int> pho_mcStatus;

	photontree->Branch("crosssection",&crosssection);
	photontree->Branch("ntotalevent", &ntotalevent);
  	photontree->Branch("phoEt",     &pho_phoEt);
  	photontree->Branch("phoEta",    &pho_phoEta);
  	photontree->Branch("phoPhi",    &pho_phoPhi);
  	photontree->Branch("phoID",    &pho_phoID);
  	photontree->Branch("phoveto",  &pho_phoveto);
  	photontree->Branch("mcType",   &mcType);
  	photontree->Branch("mcPID",    &pho_mcPID);
  	photontree->Branch("mcEta",    &pho_mcEta);
  	photontree->Branch("mcPhi",    &pho_mcPhi);
  	photontree->Branch("mcPt",     &pho_mcPt);
  	photontree->Branch("mcMomPID", &pho_mcMomPID);
  	photontree->Branch("mcGMomPID",&pho_mcGMomPID);
  	photontree->Branch("mcStatus", &pho_mcStatus);

  	TTree *egtree = new TTree("egTree","egTree");
  	float eg_phoEt(0);
  	float eg_phoEta(0);
  	float eg_phoPhi(0);
 	float eg_lepPt(0);
  	float eg_lepEta(0);
  	float eg_lepPhi(0);
  	float eg_sigMT(0);
  	float eg_sigMET(0);
  	float eg_sigMETPhi(0);
  	float eg_dPhiLepMET(0);
  	int   eg_nVertex(0);
  	float eg_dRPhoLep(0);
  	float eg_HT(0);
  	float eg_nJet(0);
  	float eg_invmass(0);
	float eg_sigMETJESup(0);
	float eg_sigMETJESdo(0);
	float eg_sigMETJERup(0);
	float eg_sigMETJERdo(0);
	float eg_sigMTJESup(0);
	float eg_sigMTJESdo(0);
	float eg_sigMTJERup(0);
	float eg_sigMTJERdo(0);
	float eg_HTJESup(0);
	float eg_HTJESdo(0);
	float eg_dPhiLepMETJESup(0);
	float eg_dPhiLepMETJESdo(0);
	float eg_dPhiLepMETJERup(0);
	float eg_dPhiLepMETJERdo(0);
  	std::vector<int> eg_mcPID;
  	std::vector<float> eg_mcEta;
  	std::vector<float> eg_mcPhi;
  	std::vector<float> eg_mcPt;
  	std::vector<int> eg_mcMomPID;
  	std::vector<int> eg_mcGMomPID;

	egtree->Branch("crosssection",&crosssection);
	egtree->Branch("ntotalevent", &ntotalevent);
  	egtree->Branch("phoEt",     &eg_phoEt);
  	egtree->Branch("phoEta",    &eg_phoEta);
  	egtree->Branch("phoPhi",    &eg_phoPhi);
  	egtree->Branch("lepPt",     &eg_lepPt);
  	egtree->Branch("lepEta",    &eg_lepEta);
  	egtree->Branch("lepPhi",    &eg_lepPhi);
	egtree->Branch("PUweight",      &PUweight);
  	egtree->Branch("sigMT",     &eg_sigMT);
  	egtree->Branch("sigMET",    &eg_sigMET);
  	egtree->Branch("sigMETPhi", &eg_sigMETPhi);
  	egtree->Branch("dPhiLepMET",&eg_dPhiLepMET);
  	egtree->Branch("nVertex",   &eg_nVertex);
  	egtree->Branch("dRPhoLep",  &eg_dRPhoLep);
  	egtree->Branch("HT",        &eg_HT);
  	egtree->Branch("nJet",      &eg_nJet);
  	egtree->Branch("invmass",   &eg_invmass);
	egtree->Branch("sigMETJESup",     &eg_sigMETJESup);
	egtree->Branch("sigMETJESdo",     &eg_sigMETJESdo);
	egtree->Branch("sigMETJERup",     &eg_sigMETJERup);
	egtree->Branch("sigMETJERdo",     &eg_sigMETJERdo);
	egtree->Branch("sigMTJESup",      &eg_sigMTJESup);
	egtree->Branch("sigMTJESdo",      &eg_sigMTJESdo);
	egtree->Branch("sigMTJERup",      &eg_sigMTJERup);
	egtree->Branch("sigMTJERdo",      &eg_sigMTJERdo);
	egtree->Branch("dPhiLepMETJESup", &eg_dPhiLepMETJESup);
	egtree->Branch("dPhiLepMETJESdo", &eg_dPhiLepMETJESdo);
	egtree->Branch("dPhiLepMETJERup", &eg_dPhiLepMETJERup);
	egtree->Branch("dPhiLepMETJERdo", &eg_dPhiLepMETJERdo);
	egtree->Branch("HTJESup",     &eg_HTJESup);
	egtree->Branch("HTJESdo",     &eg_HTJESdo);
  	egtree->Branch("mcType",   &mcType);
  	egtree->Branch("mcPID",    &eg_mcPID);
  	egtree->Branch("mcEta",    &eg_mcEta);
  	egtree->Branch("mcPhi",    &eg_mcPhi);
  	egtree->Branch("mcPt",     &eg_mcPt);
  	egtree->Branch("mcMomPID", &eg_mcMomPID);
  	egtree->Branch("mcGMomPID", &eg_mcGMomPID);

  	TTree *mgtree = new TTree("mgTree","mgTree");
  	float mg_phoEt(0);
  	float mg_phoEta(0);
  	float mg_phoPhi(0);
  	float mg_lepPt(0);
  	float mg_lepEta(0);
  	float mg_lepPhi(0);
  	float mg_sigMT(0);
  	float mg_sigMET(0);
  	float mg_sigMETPhi(0);
  	float mg_dPhiLepMET(0);
  	int   mg_nVertex(0);
  	float mg_dRPhoLep(0);
  	float mg_HT(0);
  	float mg_nJet(0);
	float mg_sigMETJESup(0);
	float mg_sigMETJESdo(0);
	float mg_sigMETJERup(0);
	float mg_sigMETJERdo(0);
	float mg_sigMTJESup(0);
	float mg_sigMTJESdo(0);
	float mg_sigMTJERup(0);
	float mg_sigMTJERdo(0);
	float mg_HTJESup(0);
	float mg_HTJESdo(0);
	float mg_dPhiLepMETJESup(0);
	float mg_dPhiLepMETJESdo(0);
	float mg_dPhiLepMETJERup(0);
	float mg_dPhiLepMETJERdo(0);
  	std::vector<int> mg_mcPID;
  	std::vector<float> mg_mcEta;
  	std::vector<float> mg_mcPhi;
  	std::vector<float> mg_mcPt;
  	std::vector<int> mg_mcMomPID;
  
	mgtree->Branch("crosssection",&crosssection);
	mgtree->Branch("ntotalevent", &ntotalevent);
  	mgtree->Branch("phoEt",     &mg_phoEt);
  	mgtree->Branch("phoEta",    &mg_phoEta);
  	mgtree->Branch("phoPhi",    &mg_phoPhi);
  	mgtree->Branch("lepPt",     &mg_lepPt);
  	mgtree->Branch("lepEta",    &mg_lepEta);
  	mgtree->Branch("lepPhi",    &mg_lepPhi);
	mgtree->Branch("PUweight",      &PUweight);
  	mgtree->Branch("sigMT",     &mg_sigMT);
  	mgtree->Branch("sigMET",    &mg_sigMET);
  	mgtree->Branch("sigMETPhi", &mg_sigMETPhi);
  	mgtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
  	mgtree->Branch("nVertex",   &mg_nVertex);
  	mgtree->Branch("dRPhoLep",  &mg_dRPhoLep);
  	mgtree->Branch("HT",        &mg_HT);
  	mgtree->Branch("nJet",      &mg_nJet);
	mgtree->Branch("sigMETJESup",     &mg_sigMETJESup);
	mgtree->Branch("sigMETJESdo",     &mg_sigMETJESdo);
	mgtree->Branch("sigMETJERup",     &mg_sigMETJERup);
	mgtree->Branch("sigMETJERdo",     &mg_sigMETJERdo);
	mgtree->Branch("sigMTJESup",      &mg_sigMTJESup);
	mgtree->Branch("sigMTJESdo",      &mg_sigMTJESdo);
	mgtree->Branch("sigMTJERup",      &mg_sigMTJERup);
	mgtree->Branch("sigMTJERdo",      &mg_sigMTJERdo);
	mgtree->Branch("dPhiLepMETJESup", &mg_dPhiLepMETJESup);
	mgtree->Branch("dPhiLepMETJESdo", &mg_dPhiLepMETJESdo);
	mgtree->Branch("dPhiLepMETJERup", &mg_dPhiLepMETJERup);
	mgtree->Branch("dPhiLepMETJERdo", &mg_dPhiLepMETJERdo);
	mgtree->Branch("HTJESup",     &mg_HTJESup);
	mgtree->Branch("HTJESdo",     &mg_HTJESdo);
  	mgtree->Branch("mcType",   &mcType);
  	mgtree->Branch("mcPID",    &mg_mcPID);
  	mgtree->Branch("mcEta",    &mg_mcEta);
  	mgtree->Branch("mcPhi",    &mg_mcPhi);
  	mgtree->Branch("mcPt",     &mg_mcPt);
  	mgtree->Branch("mcMomPID", &mg_mcMomPID);
	//*********** histo list **********************//
  	rawData raw(es, datatype);
  	std::vector<mcData>  MCData;
  	std::vector<recoPhoton> Photon;
  	std::vector<recoMuon>   Muon;
  	std::vector<recoEle>   Ele;
  	std::vector<recoJet>   JetCollection;
	/*********************************************/
  	float MET(0);
  	float METPhi(0);
	float MET_T1JERUp(0);
	float MET_T1JERDo(0);
	float MET_T1JESUp(0);
	float MET_T1JESDo(0);	
	float	METPhi_T1JESUp(0);
	float	METPhi_T1JESDo(0);
	float	METPhi_T1UESUp(0);
	float	METPhi_T1UESDo(0);
  	int nVtx(0);
  	int jetNumber(0);
  	int METFilter(0);

	int theoryCut(0);
	int highEta(0);
  	std::cout << "Total evetns : " << nEvts << std::endl;
  
    	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      		if (ievt%100000==0){ std::cout << " -- Processing event " << ievt << std::endl;  std::cout << "highEta " << highEta << " theory " << theoryCut << std::endl;}

        	raw.GetData(es, ievt);
        	MCData.clear();
        	Photon.clear();
        	Muon.clear();
        	Ele.clear();
		JetCollection.clear();
        	if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
        	for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        	for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
        	for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
		for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
        	MET = raw.pfMET;
        	METPhi = raw.pfMETPhi;
		MET_T1JERUp = raw.pfMET_T1JERUp;
		MET_T1JERDo = raw.pfMET_T1JERDo;
		MET_T1JESUp = raw.pfMET_T1JESUp;
		MET_T1JESDo = raw.pfMET_T1JESDo;
		METPhi_T1JESUp = raw.pfMETPhi_T1JESUp;
		METPhi_T1JESDo = raw.pfMETPhi_T1JESDo;
		METPhi_T1UESUp = raw.pfMETPhi_T1UESUp;
		METPhi_T1UESDo = raw.pfMETPhi_T1UESDo;
        	METFilter = raw.metFilters;
        	nVtx = raw.nVtx;
		if(RunYear==2016 and preVFP==1)		PUweight = getPUESF16preVFP(nVtx);
		else if(RunYear==2016 and preVFP==0)	PUweight = getPUESF16(nVtx);
               	else if(RunYear==2017)			PUweight = getPUESF17(nVtx);
                else if(RunYear==2018)			PUweight = getPUESF18(nVtx);
	
		double ISRJetPt = 0;
		TLorentzVector JetVec(0,0,0,0);	
		for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
			if(!itJet->passSignalSelection())continue;
			if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))JetVec = JetVec + itJet->getP4();}

		ISRJetPt = JetVec.Pt();
		double reweightF(1);	
		/*if(ISRJetPt < 50)reweightF = 1.015;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF  = 1.110;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.845;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.715;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF =	0.730;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF =	0.732;
		else if(ISRJetPt >= 300)reweightF =  0.642; 
		*/
		bool hasHighPho(false);
		std::vector<mcData>::iterator sigMCNeu = MCData.begin();
		std::vector<mcData>::iterator sigMCLep = MCData.begin();
		std::vector<mcData>::iterator sigMCPho = MCData.begin();
		bool hassigMCNeu(false), hassigMCLep(false), hassigMCPho(false);	
		for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 1.0)continue;
				if(!hassigMCLep){
					if( abs(itMC->getPID()) == 15){sigMCLep = itMC; hassigMCLep = true;}
					else if(abs(itMC->getPID()) == 11 || abs(itMC->getPID()) == 13){sigMCLep = itMC; hassigMCLep = true;}
				}
				if(!hassigMCNeu){
					if( abs(itMC->getPID()) == 16){sigMCNeu = itMC; hassigMCNeu = true;}
					else if( abs(itMC->getPID()) == 12 || abs(itMC->getPID()) == 14 ){sigMCNeu = itMC; hassigMCNeu = true;}
				}
				if( abs(itMC->getPID()) == 22 && itMC->getPt() > 15 && fabs(itMC->getEta()) < 2.5){ sigMCPho = itMC; hassigMCPho = true;}
				if( abs(itMC->getPID()) == 22 && itMC->getPt() > 175 && fabs(itMC->getEta()) < 1.4442)hasHighPho = true;
		}

		if(hasHighPho)highEta+= 1;
		if( hassigMCNeu && hassigMCLep && hassigMCPho){
				if( sigMCLep->getPt() > 25 && sigMCNeu->getPt() > 25 && fabs(sigMCLep->getEta()) < 2.5 && DeltaR( sigMCLep->getEta(), sigMCLep->getPhi(), sigMCPho->getEta(), sigMCPho->getPhi()) > 0.7 &&  (sigMCLep->getP4()+sigMCNeu->getP4()).M() > 40){
						theoryCut += 1;
					 // if( sigMCPho->getPt() > 130)highEta+= 1;
				}
		}	


        	if(raw.nPho <1)continue;

  		pho_phoEt.clear();
  		pho_phoEta.clear();
  		pho_phoPhi.clear();
  		pho_phoID.clear();
  		pho_phoveto.clear();
		pho_mcPID.clear();
		pho_mcEta.clear();
		pho_mcPhi.clear();
		pho_mcPt.clear();
		pho_mcMomPID.clear();
		pho_mcGMomPID.clear();
		pho_mcStatus.clear();

		{all_totalEvent += 1; all_reweight += reweightF; }
		bool hasnopixel(false);
		bool hasid(false);
		bool hasleading(false);
        	bool hasegPho(false);
        	bool hasmgPho(false);
        	std::vector<recoPhoton>::iterator egsignalPho = Photon.begin();
        	std::vector<recoPhoton>::iterator mgsignalPho = Photon.begin();
		std::vector<recoPhoton>::iterator leadingphoton = Photon.begin();
       		for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){

					all_phoEt->Fill(itpho->getCalibEt());
					pho_phoEt.push_back(itpho->getEt());
					pho_phoEta.push_back(itpho->getEta());
					pho_phoPhi.push_back(itpho->getPhi());
					if(itpho->PixelSeed()==0){
						nopixel_phoEt->Fill(itpho->getCalibEt());
						pho_phoveto.push_back(1);
						if(!hasnopixel){nopixel_totalEvent += 1; nopixel_reweight += reweightF; hasnopixel=true;}
					}
					else pho_phoveto.push_back(0);
					if(itpho->isLoose())pho_phoID.push_back(1);
					else pho_phoID.push_back(0);
					if(itpho->passSignalSelection() && itpho->PixelSeed()==0){
						id_phoEt->Fill(itpho->getCalibEt());
						if(!hasid){id_totalEvent += 1; id_reweight += reweightF; hasid = true;}
					}

					if(itpho->getR9() < 0.5)continue;
          				if(!itpho->passSignalSelection())continue;
					bool PixelVeto = itpho->PixelSeed()==0? true: false;
					bool GSFveto(true);
					bool FSRVeto(true);
					for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					 	if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
					 	if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
					}
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					 	if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;
					if(PixelVeto && !hasleading){
						hasleading = true;
						leadingphoton = itpho;
					}
					if(GSFveto && PixelVeto && FSRVeto){
						if((itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6))){
					 		if(!hasegPho){
								hasegPho=true;
					 			egsignalPho = itpho;
								egamma_phoEt->Fill(itpho->getCalibEt());
							}
						}
						if((RunYear == 2016 && (itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30))) || ((RunYear == 2017 || RunYear == 2018) && (itpho->fireDoubleTrg(27) || itpho->fireDoubleTrg(36)))) {
					 		if(!hasmgPho){
								hasmgPho=true;
					 			mgsignalPho = itpho;
								mugamma_phoEt->Fill(itpho->getCalibEt());
							}
						}
					}
        	}
	        for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 1.0)continue;
				pho_mcPID.push_back(itMC->getPID());
				pho_mcMomPID.push_back(itMC->getMomPID());
				pho_mcEta.push_back(itMC->getEta());
				pho_mcPhi.push_back(itMC->getPhi());
				pho_mcPt.push_back(itMC->getEt());
				pho_mcStatus.push_back(itMC->getStatus());
		}
		photontree->Fill();
       
				//if(leadingphoton->getCalibEt() > 50)continue;
 
        	bool hasEle(false);
        	std::vector<recoEle>::iterator signalEle = Ele.begin();
        	if(hasegPho){
          		for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
            			if(hasEle)break;
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
				if((RunYear==2016 && (itEle->fireTrgs(21) || itEle->fireTrgs(22))) || ((RunYear == 2017 || RunYear == 2018) && (itEle->fireTrgs(43) || itEle->fireTrgs(44)))) {
					if(itEle->passSignalSelection()){
						hasEle=true; 
						signalEle = itEle;
					}
				}
          		}
        	}

        	bool hasMu(false);
        	std::vector<recoMuon>::iterator signalMu = Muon.begin();
        	if(hasmgPho){
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
			if(hasMu)break;
			if(itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || ((RunYear == 2016 && itMu->fireSingleTrg(22)) || (RunYear == 2017 && itMu->fireSingleTrg(32)) || (RunYear == 2018 && itMu->fireSingleTrg(32)))) {
				if(itMu->passSignalSelection()){
					hasMu=true; 
					signalMu = itMu;
				}
			}
			}
        	}

///*************************   eg filters *****************************//
 //       if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){
 //         double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
 //         if(dReg>0.8){
 //           if(fabs((egsignalPho->getP4()+signalEle->getP4()).M() - 91.188) > 10.0){
 //             if(METFilter == 0){

        
		if(hasegPho && hasEle){
			{egamma_totalEvent += 1; egamma_reweight += reweightF;}
          		double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
          		if(dReg>0.8){
            		{
              		if(METFilter == 0){
				float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
				float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
	        	    	eg_phoEt = egsignalPho->getCalibEt();
	    			eg_phoEta= egsignalPho->getEta();
        		        eg_phoPhi= egsignalPho->getPhi();
			 	eg_lepPt = signalEle->getCalibPt();
 			    	eg_lepEta= signalEle->getEta();
                		eg_lepPhi= signalEle->getPhi();
 			    	eg_sigMT = MT;
			    	eg_sigMET= MET;
                		eg_sigMETPhi = METPhi;
			    	eg_dPhiLepMET = deltaPhi; 
			    	eg_nVertex = nVtx; 
			    	eg_dRPhoLep= dReg;
                		eg_invmass = (egsignalPho->getP4()+signalEle->getP4()).M();
				eg_sigMETJESup = MET_T1JESUp;
				eg_sigMETJESdo = MET_T1JESDo;
				eg_sigMETJERup = MET_T1JERUp;
				eg_sigMETJERdo = MET_T1JERDo;
				eg_dPhiLepMETJESup = DeltaPhi(signalEle->getPhi(), METPhi_T1JESUp);
				eg_dPhiLepMETJESdo = DeltaPhi(signalEle->getPhi(), METPhi_T1JESDo);
				eg_dPhiLepMETJERup = deltaPhi;
				eg_dPhiLepMETJERdo = deltaPhi;
				eg_sigMTJESup = sqrt(2*MET_T1JESUp*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESup)));
				eg_sigMTJESdo = sqrt(2*MET_T1JESDo*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESdo)));
				eg_sigMTJERup = sqrt(2*MET_T1JERUp*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERup)));
				eg_sigMTJERdo = sqrt(2*MET_T1JERDo*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERdo)));
							
				eg_nJet = 0;
				eg_HT = 0;
				eg_HTJESup = 0;
				eg_HTJESdo = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalPho->getEta(),egsignalPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(),signalEle->getPhi()) <= 0.4)continue;
					eg_nJet += 1;
					eg_HT += itJet->getPt();
					eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
					eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
				}	
 				if(datatype == MC){
					eg_mcPID.clear();
					eg_mcEta.clear();
					eg_mcPhi.clear();
					eg_mcPt.clear();
					eg_mcMomPID.clear();
					eg_mcGMomPID.clear();
					for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						if(itMC->getEt() < 5.0)continue;
						float mcdR = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						float mcdRe = DeltaR(signalEle->getEta(),signalEle->getPhi(), itMC->getEta(), itMC->getPhi());
						//if(mcdR < 0.3 || mcdRe < 0.3){
						eg_mcPID.push_back(itMC->getPID());
						eg_mcMomPID.push_back(itMC->getMomPID());
						eg_mcEta.push_back(itMC->getEta());      
						eg_mcPhi.push_back(itMC->getPhi());
						eg_mcPt.push_back(itMC->getEt());
						eg_mcGMomPID.push_back(itMC->getGMomPID());
										//}
					}
				}

                		egtree->Fill();
 
              }//MET Filter
            }// Z mass Filter
          }//dR filter
        }// ele + pho candidate
       
//**********************  mg filter **************************************//         
//        if(hasmgPho && hasMu && (((raw.HLTEleMuX >> 8)&1)!=0 || ((raw.HLTEleMuX >> 41)&1)!=0)){
//					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
//					if(dRmg>0.8){
//						if(METFilter == 0){ 
 
       		if(hasmgPho && hasMu){
			{mugamma_totalEvent += 1; mugamma_reweight += reweightF;}
			double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
			if(dRmg>0.8){
				if(METFilter == 0){ 

					float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));

					mg_phoEt = mgsignalPho->getCalibEt();
					mg_phoEta= mgsignalPho->getEta();
					mg_phoPhi= mgsignalPho->getPhi();
					mg_lepPt = signalMu->getPt();
					mg_lepEta= signalMu->getEta();
					mg_lepPhi= signalMu->getPhi();
					mg_sigMT = MT;
					mg_sigMET= MET;
					mg_sigMETPhi = METPhi;
					mg_dPhiLepMET = deltaPhi;
					mg_nVertex = nVtx;
					mg_dRPhoLep= dRmg;
					mg_sigMETJESup = MET_T1JESUp;
					mg_sigMETJESdo = MET_T1JESDo;
					mg_sigMETJERup = MET_T1JERUp;
					mg_sigMETJERdo = MET_T1JERDo;
					mg_dPhiLepMETJESup = DeltaPhi(signalMu->getPhi(), METPhi_T1JESUp);
					mg_dPhiLepMETJESdo = DeltaPhi(signalMu->getPhi(), METPhi_T1JESDo);
					mg_dPhiLepMETJERup = deltaPhi;
					mg_dPhiLepMETJERdo = deltaPhi;
					mg_sigMTJESup = sqrt(2*MET_T1JESUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESup)));
					mg_sigMTJESdo = sqrt(2*MET_T1JESDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESdo)));
					mg_sigMTJERup = sqrt(2*MET_T1JERUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERup)));
					mg_sigMTJERdo = sqrt(2*MET_T1JERDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERdo)));

					mg_nJet = 0;
					mg_HT = 0;
					mg_HTJESup = 0;
					mg_HTJESdo = 0;
					for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), mgsignalPho->getEta(),mgsignalPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= 0.4)continue;
							mg_nJet += 1;
							mg_HT += itJet->getPt();
							mg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
							mg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
					}	

					if(datatype == MC){
						 mg_mcPID.clear();
						 mg_mcEta.clear();
						 mg_mcPhi.clear();
						 mg_mcPt.clear();
						 mg_mcMomPID.clear();
						 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	if(itMC->getEt() < 1.0)continue;
							float mcdRmu = DeltaR(signalMu->getEta(), signalMu->getPhi(), itMC->getEta(), itMC->getPhi());
						 	float mcdR = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						 	if(mcdR < 0.3){
						 		mg_mcPID.push_back(itMC->getPID());
						 		mg_mcMomPID.push_back(itMC->getMomPID());
						 		mg_mcEta.push_back(itMC->getEta());
						 		mg_mcPhi.push_back(itMC->getPhi());
						 		mg_mcPt.push_back(itMC->getEt());
						 	}
						}
					}

				 	mgtree->Fill();
           }//MET Filter
         }//dR Filter
       	}//Candidate Filter
 
	}//loop on  events

	outputfile->Write();
	outputfile->Close();

	std::cout << "all_totalEvent" << all_totalEvent << " all_reweight " << all_reweight << " ratio " << all_totalEvent/all_reweight << std::endl;
	std::cout << "nopixel_totalEvent" << nopixel_totalEvent << " nopixel_reweight " << nopixel_reweight << " ratio " << nopixel_totalEvent/nopixel_reweight << std::endl;
	std::cout << "id_totalEvent" << id_totalEvent << " id_reweight " << id_reweight << " ratio " << id_totalEvent/id_reweight << std::endl;
	std::cout << "egamma_totalEvent" << egamma_totalEvent << " egamma_reweight " << egamma_reweight << " ratio " << egamma_totalEvent/egamma_reweight << std::endl;
	std::cout << "mugamma_totalEvent" << mugamma_totalEvent << " mugamma_reweight " << mugamma_reweight << " ratio " << mugamma_totalEvent/mugamma_reweight << std::endl;
}

int main(int argc, char** argv)
{
    if(argc < 2)
      cout << "You have to provide two arguments!!\n";
    analysis_Mixing(atoi(argv[1]),argv[2]);
    return 0;
}
