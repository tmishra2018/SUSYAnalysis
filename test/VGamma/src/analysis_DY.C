// similar to analysis_VGamma.C script
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
#include "TProfile2D.h"
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_jet.h"

bool isHardLepton(int momID){
		if((fabs(momID) >= 0 && fabs(momID) <= 6) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
		else return false;
}
		
bool isHardPhoton(int momID){
		if((fabs(momID) >= 0 && fabs(momID) <= 16) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
		else return false;
}

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
		

void analysis_DY(){//main 
	
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG130Calib.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG130Calib.log"); 

  logfile << "analysis_VGamma()" << std::endl;

  RunType datatype(MC);
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WGJets_MonoPhoton_PtG-130_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");
  logfile << "root://cmseos.fnal.gov//store/user/msun/MCSummer16/WGJets_MonoPhoton_PtG-130_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root" << std::endl;
	logfile << "muon MiniIso" << std::endl;
 
  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

  //int mcType = MCType::WGJet40;
  int mcType = MCType::WGJet130;
  //int mcType = MCType::ZGInclusive;
  //int mcType = MCType::WGJetInclusive;
  //int mcType = MCType::TT;
  //int mcType = MCType::DYLL50;
  //int mcType = MCType::TTG;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;

  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();
	float PUweight(1);
	float llmass(0);
	int   nBJet(0); 
//************ Signal Tree **********************//
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
	int   eg_nISRJet(0);
	float eg_ISRJetPt(0);
  float eg_invmass(0);
	float eg_bosonPt(0);
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
  std::vector<int> eg_mcStatus;

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
  egtree->Branch("nBJet",     &nBJet);
  egtree->Branch("nISRJet",   &eg_nISRJet);
	egtree->Branch("ISRJetPt",  &eg_ISRJetPt);
  egtree->Branch("invmass",   &eg_invmass);
	egtree->Branch("bosonPt",   &eg_bosonPt);
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
  egtree->Branch("mcStatus", &eg_mcStatus);
	egtree->Branch("llmass",   &llmass);

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
	float mg_threeMass(0);
  int   mg_nVertex(0);
  float mg_dRPhoLep(0);
  float mg_HT(0);
  float mg_nJet(0);
  float mg_nISRJet(0);
	float mg_ISRJetPt(0);
	float mg_bosonPt(0);
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
  std::vector<int> mg_mcGMomPID;
	std::vector<int> mg_mcStatus;
 
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
	mgtree->Branch("threeMass", &mg_threeMass);
  mgtree->Branch("nVertex",   &mg_nVertex);
  mgtree->Branch("dRPhoLep",  &mg_dRPhoLep);
  mgtree->Branch("HT",        &mg_HT);
  mgtree->Branch("nJet",      &mg_nJet);
  mgtree->Branch("nBJet",     &nBJet);
  mgtree->Branch("nISRJet",   &mg_nISRJet);
	mgtree->Branch("ISRJetPt",  &mg_ISRJetPt);
	mgtree->Branch("bosonPt",   &mg_bosonPt);
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
  mgtree->Branch("mcGMomPID", &mg_mcGMomPID);
  mgtree->Branch("mcStatus", &mg_mcStatus);
	mgtree->Branch("llmass",   &llmass);
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
	float ISRJetPt(0);
  logfile << "RunType: " << datatype << std::endl;

  std::cout << "Total evetns : " << nEvts << std::endl;

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;

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
				
				if(RunYear==2016 and preVFP==1)PUweight = getPUESF16preVFP(nVtx);
				else if(RunYear==2016 and preVFP==0)PUweight = getPUESF16(nVtx);
                        	else if(RunYear==2017)PUweight = getPUESF17(nVtx);
                        	else if(RunYear==2018)PUweight = getPUESF18(nVtx);
	
        			if(raw.nPho <1)continue;

				nBJet = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(itJet->getPt() < 20)continue;
					if(itJet->isBJet())nBJet+=1;
				}
				//process MC first  ****************************/
				llmass = -1;
				double WGPt = -1;
				if(datatype == MC){
				 std::vector<mcData>::iterator posMcMu; 
				 bool hasposMcMu(false);
				 std::vector<mcData>::iterator negMcMu;
				 bool hasnegMcMu(false);
				 std::vector<mcData>::iterator posMcTau; 
				 bool hasposMcTau(false);
				 std::vector<mcData>::iterator negMcTau;
				 bool hasnegMcTau(false);
				 std::vector<mcData>::iterator posMcEle; 
				 bool hasposMcEle(false);
				 std::vector<mcData>::iterator negMcEle;
				 bool hasnegMcEle(false);
				 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					switch(itMC->getPID()){
						case 13: if(isHardLepton(itMC->getMomPID())){posMcMu = itMC; hasposMcMu=true;};break;
						case -13:if(isHardLepton(itMC->getMomPID())){negMcMu = itMC; hasnegMcMu=true;};break;
						case 15: if(isHardLepton(itMC->getMomPID())){posMcTau = itMC; hasposMcTau=true;};break;
						case -15:if(isHardLepton(itMC->getMomPID())){negMcTau = itMC; hasnegMcTau=true;};break;
						case 11: if(isHardLepton(itMC->getMomPID())){posMcEle = itMC; hasposMcEle=true;};break;
						case -11:if(isHardLepton(itMC->getMomPID())){negMcEle = itMC; hasnegMcEle=true;};break;
						case 23: WGPt = itMC->getPt(); break;
						case 24: WGPt = itMC->getPt(); break;
						case -24: WGPt = itMC->getPt(); break;
						default: break;
					}
				 }
				 if(hasposMcMu && hasnegMcMu)llmass = (posMcMu->getP4() + negMcMu->getP4()).M();
				 else if(hasposMcTau && hasnegMcTau)llmass = (posMcTau->getP4() + negMcTau->getP4()).M();
				 else if(hasposMcEle && hasnegMcEle)llmass = (posMcEle->getP4() + negMcEle->getP4()).M();
				}

        bool hasegPho(false);
        bool hasmgPho(false);
        std::vector<recoPhoton>::iterator egsignalPho = Photon.begin();
        std::vector<recoPhoton>::iterator mgsignalPho = Photon.begin();
       	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
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
					if(GSFveto && PixelVeto && FSRVeto){
						if((itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6))){
					 		if(!hasegPho){
								hasegPho=true;
					 			egsignalPho = itpho;
							}
						}
						if(itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30)){
					 		if(!hasmgPho){
								hasmgPho=true;
					 			mgsignalPho = itpho;
							}
						}
					}
        }
     
				if(hasegPho){ 
					if(egsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasegPho=false;
					if(egsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasegPho=false;
				}
				if(hasmgPho){ 
					if(mgsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasmgPho=false;
					if(mgsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasmgPho=false;
				}
 
        bool hasEle(false);
        std::vector<recoEle>::iterator signalEle = Ele.begin();
        if(hasegPho){
          for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
            if(hasEle)break;
						if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
						if(itEle->fireTrgs(21) || itEle->fireTrgs(22)){
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
						if(itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22)){
							if(itMu->passSignalSelection()){
								hasMu=true; 
								signalMu = itMu;
							}
						}
					}
        }
	
				ISRJetPt = 0;
				TLorentzVector JetVec(0,0,0,0);	
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))JetVec = JetVec + itJet->getP4();
				}
				ISRJetPt = JetVec.Pt();	
///*************************   eg filters *****************************//
 //       if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){
 //         double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
 //         if(dReg>0.8){
 //           if(fabs((egsignalPho->getP4()+signalEle->getP4()).M() - 91.188) > 10.0){
 //             if(METFilter == 0){

        if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){
          double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
          if(dReg>0.8){
            	if(fabs((egsignalPho->getP4()+signalEle->getP4()).M() - 91.188) > 10.0){
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
								eg_bosonPt = WGPt;
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
								eg_nISRJet = 0;
								eg_HT = 0;
								eg_HTJESup = 0;
								eg_HTJESdo = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
								//	bool isLep(false);
          			//	for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
								//		if(!itEle->isMedium() || itEle->getMiniIso() > 0.1 || itEle->getPt() < 25)continue;
								//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itEle->getEta(), itEle->getPhi()) <= 0.4 && fabs(itJet->getPt() - itEle->getPt())/itEle->getPt() < 1)isLep=true;
								//	}
								//	for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
								//		if(!itMu->isMedium() || itMu->getMiniIso() > 0.2 || itMu->getPt() < 25)continue;
								//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4 && fabs(itJet->getPt() - itMu->getPt())/itMu->getPt() < 1)isLep=true;
								//	}
								//	if(isLep)continue;
									
									if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))eg_nISRJet += 1;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(), signalEle->getPhi()) <= 0.4)continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalPho->getEta(),egsignalPho->getPhi()) <= 0.4)continue;	
									eg_nJet += 1;
									eg_HT += itJet->getPt();
									eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
									eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
								}
								eg_ISRJetPt = ISRJetPt;	
 
								if(datatype == MC){
									eg_mcPID.clear();
									eg_mcEta.clear();
									eg_mcPhi.clear();
									eg_mcPt.clear();
									eg_mcMomPID.clear();
									eg_mcGMomPID.clear();
									eg_mcStatus.clear();
									for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
										if(itMC->getEt() < 5.0)continue;
										eg_mcPID.push_back(itMC->getPID());
										eg_mcMomPID.push_back(itMC->getMomPID());
										eg_mcEta.push_back(itMC->getEta());      
										eg_mcPhi.push_back(itMC->getPhi());
										eg_mcPt.push_back(itMC->getEt());
										eg_mcGMomPID.push_back(itMC->getGMomPID());
										eg_mcStatus.push_back(itMC->getStatus());
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
 
       if(hasmgPho && hasMu && (((raw.HLTEleMuX >> 8)&1)!=0 || ((raw.HLTEleMuX >> 41)&1)!=0)){
					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
					if(dRmg>0.8){
						if(METFilter == 0){ 

						double deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
						double MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));
						float ThreeBodyMass = sqrt(2*MET*(mgsignalPho->getP4()+ signalMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (mgsignalPho->getP4()+signalMu->getP4()).Phi(), 0, METPhi))));
						std::cout << MET<< "  "<< signalMu->getPt() << " " <<  deltaPhi <<" " << std::cos(deltaPhi) << " " << MT << std::endl;

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
						mg_threeMass = ThreeBodyMass;
						mg_nVertex = nVtx;
						mg_dRPhoLep= dRmg;
						mg_bosonPt = WGPt;
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
						mg_nISRJet = 0;
						mg_HT = 0;
						mg_HTJESup = 0;
						mg_HTJESdo = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
						//	bool isLep(false);
						//	for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
						//		if(!itEle->isMedium() || itEle->getMiniIso() > 0.1 || itEle->getPt() < 25)continue;
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itEle->getEta(), itEle->getPhi()) <= 0.4 && fabs(itJet->getPt() - itEle->getPt())/itEle->getPt() < 1)isLep=true;
						//	}
						//	for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
						//		if(!itMu->isMedium() || itMu->getMiniIso() > 0.2 || itMu->getPt() < 25)continue;
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4 && fabs(itJet->getPt() - itMu->getPt())/itMu->getPt() < 1)isLep=true;
						//	}
						//	if(isLep)continue;

							if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))mg_nISRJet += 1;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= 0.4)continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), mgsignalPho->getEta(),mgsignalPho->getPhi()) <= 0.4)continue;	
							mg_nJet += 1;
							mg_HT += itJet->getPt();
							mg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
							mg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
						}	
						mg_ISRJetPt = ISRJetPt;

						if(datatype == MC){
						 mg_mcPID.clear();
						 mg_mcEta.clear();
						 mg_mcPhi.clear();
						 mg_mcPt.clear();
						 mg_mcMomPID.clear();
						 mg_mcGMomPID.clear();
						 mg_mcStatus.clear();
						 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	if(itMC->getEt() < 1.0)continue;
							float mcdRmu = DeltaR(signalMu->getEta(), signalMu->getPhi(), itMC->getEta(), itMC->getPhi());
						 	float mcdR = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						 	if(mcdR < 0.3 || mcdRmu < 0.3){
						 		mg_mcPID.push_back(itMC->getPID());
						 		mg_mcMomPID.push_back(itMC->getMomPID());
						 		mg_mcGMomPID.push_back(itMC->getGMomPID());
						 		mg_mcEta.push_back(itMC->getEta());
						 		mg_mcPhi.push_back(itMC->getPhi());
						 		mg_mcPt.push_back(itMC->getEt());
								mg_mcStatus.push_back(itMC->getStatus());
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
logfile.close();

}


