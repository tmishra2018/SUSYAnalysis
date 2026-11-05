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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_tools.h"

int RunYear = 2016;
void analysis_TChiWG(){//main  

  gSystem->Load("../../lib/libAnaClasses.so");

  TChain* es = new TChain("ggNtuplizer/EventTree");
  //es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/SMS-TChiWg_mChi-1000_mLSP-1/SMS-TChiWg_mChi-1000_mLSP-1_2016.root"); // working
  //es->Add(Form("/eos/uscms/store/group/lpcsusyhad/Tribeni/TChiWG_FastSim/TChiWG_FastSim_%d.root",RunYear)); // not working
  //es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/TChiWG_FastSim/TChiWG_FastSim_2018.root"); // not working
  es->Add("/eos/uscms/store/user/msun/Signal/SMS-TChiWG_TuneCUETP8M1_RunIISummer16MiniAODv2.root"); // working
  //es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/ggtree_mc_t5wg_fullTree.root"); // working


  RunType datatype(MC); 
  TFile *outputfile = TFile::Open(Form("/eos/uscms/store/group/lpcsusyhad/Tribeni/signal_trees/resTree_TChiWG_%d.root",RunYear),"RECREATE");
  outputfile->cd();
  TTree *tree = new TTree("SUSYtree","SUSYtree");
  float mcPhotonEt(0);
  float mcPhotonEta(0);
  float mcPhotonPhi(0);
  float recoPhotonEt(0);
  float recoPhotonEta(0);
  float recoPhotonPhi(0);
  float PhoR9(0);
  float PhoHoverE(0);
  float PhoSigma(0);
  float PhoChIso(0);
  float PhoNeuIso(0);
  float PhoPhoIso(0);
  bool  PhoPassID;
  float mcElePt(0);
  float mcEleEta(0);
  float mcElePhi(0);
  float recoEleEt(0);
  float recoEleEta(0);
  float recoElePhi(0);
  float EleR9(0);
  float EleHoverE(0);
  float EleSigma(0);
  float EleChIso(0);
  float EleNeuIso(0);
  float ElePhoIso(0);
  float EleMiniIso(0);
  float EledEtaIn(0);
  float EledPhiIn(0);
  float EleD0(0);
  float EleDz(0);
  float EleooEmooP(0);
  bool  ElePassID;
  float Mchagino(0);
  float Mneutralino(0);
  float dRPhoEle(0);
  float Invmass(0);
  int   nJet(0);
  float HT(0);
	float sigMET(0);
  double MT_(0), ThreeBodyMass_(0);
  tree->Branch("MT",&MT_);
  tree->Branch("ThreeBodyMass",&ThreeBodyMass_);
  tree->Branch("mcPhotonEt"     ,&mcPhotonEt); 
  tree->Branch("mcPhotonEta"    ,&mcPhotonEta);
  tree->Branch("mcPhotonPhi"    ,&mcPhotonPhi);
  tree->Branch("recoPhotonEt",   &recoPhotonEt);
  tree->Branch("recoPhotonEta",  &recoPhotonEta);
  tree->Branch("recoPhotonPhi",  &recoPhotonPhi);
  tree->Branch("PhoR9"          ,&PhoR9);
  tree->Branch("PhoHoverE"      ,&PhoHoverE);
  tree->Branch("PhoSigma"       ,&PhoSigma);
  tree->Branch("PhoChIso"       ,&PhoChIso);
  tree->Branch("PhoNeuIso"      ,&PhoNeuIso);
  tree->Branch("PhoPhoIso"      ,&PhoPhoIso);
  tree->Branch("PhoPassID"      ,&PhoPassID);
  tree->Branch("mcElePt"        ,&mcElePt);
  tree->Branch("mcEleEta"       ,&mcEleEta);
  tree->Branch("mcElePhi"       ,&mcElePhi);
  tree->Branch("recoEleEt",     &recoEleEt);
  tree->Branch("recoEleEta",    &recoEleEta);
  tree->Branch("recoElePhi",    &recoElePhi);
  tree->Branch("EleR9"          ,&EleR9);
  tree->Branch("EleHoverE"      ,&EleHoverE);
  tree->Branch("EleSigma"       ,&EleSigma);
  tree->Branch("EleChIso"       ,&EleChIso);
  tree->Branch("EleNeuIso"      ,&EleNeuIso);
  tree->Branch("ElePhoIso"      ,&ElePhoIso);
  tree->Branch("EleMiniIso"     ,&EleMiniIso);
  tree->Branch("EledEtaIn"      ,&EledEtaIn);
  tree->Branch("EledPhiIn"      ,&EledPhiIn);
  tree->Branch("EleD0"          ,&EleD0);
  tree->Branch("EleDz"          ,&EleDz);
  tree->Branch("EleooEmooP"     ,&EleooEmooP);
  tree->Branch("ElePassID"      ,&ElePassID);
  tree->Branch("Mchagino"       ,&Mchagino);
  tree->Branch("Mneutralino"        ,&Mneutralino);
  tree->Branch("dRPhoEle"       ,&dRPhoEle);
  tree->Branch("Invmass"        ,&Invmass);    
  tree->Branch("nJet"           ,&nJet);
  tree->Branch("HT"             ,&HT);
	tree->Branch("sigMET",         &sigMET);

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

  egtree->Branch("Mchagino",  &Mchagino);
  egtree->Branch("Mneutralino",&Mneutralino);
  egtree->Branch("phoEt",     &eg_phoEt);
  egtree->Branch("phoEta",    &eg_phoEta);
  egtree->Branch("phoPhi",    &eg_phoPhi);
  egtree->Branch("lepPt",     &eg_lepPt);
  egtree->Branch("lepEta",    &eg_lepEta);
  egtree->Branch("lepPhi",    &eg_lepPhi);
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

  mgtree->Branch("Mchagino",  &Mchagino);
  mgtree->Branch("Mneutralino",&Mneutralino);
  mgtree->Branch("phoEt",     &mg_phoEt);
  mgtree->Branch("phoEta",    &mg_phoEta);
  mgtree->Branch("phoPhi",    &mg_phoPhi);
  mgtree->Branch("lepPt",     &mg_lepPt);
  mgtree->Branch("lepEta",    &mg_lepEta);
  mgtree->Branch("lepPhi",    &mg_lepPhi);
  mgtree->Branch("sigMT",     &mg_sigMT);
  mgtree->Branch("sigMET",    &mg_sigMET);
  mgtree->Branch("sigMETPhi", &mg_sigMETPhi);
  mgtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
	mgtree->Branch("threeMass", &mg_threeMass);
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

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoEle> Ele;
  std::vector<recoMuon>   Muon;
  std::vector<recoJet>   JetCollection;
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

    const unsigned nEvts = es->GetEntries(); 
    std::cout << "total event : " << nEvts << std::endl;


    int  n750(0);
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

			mcPhotonEt=0;
			mcPhotonEta=0;
			mcPhotonPhi=0;
			recoPhotonEt=0;
			recoPhotonEta=0;
			recoPhotonPhi=0;
			PhoR9=0;
			PhoHoverE=0;
			PhoSigma=0;
			PhoChIso=0;
			PhoNeuIso=0;
			PhoPhoIso=0;
			PhoPassID=0;
			mcElePt=0;
			mcEleEta=0;
			mcElePhi=0;
			recoEleEt=0;
			recoEleEta=0;
			recoElePhi=0;
			EleR9=0;
			EleHoverE=0;
			EleSigma=0;
			EleChIso=0;
			EleNeuIso=0;
			ElePhoIso=0;
			EleMiniIso=0;
			EledEtaIn=0;
			EledPhiIn=0;
			EleD0=0;
			EleDz=0;
			EleooEmooP=0;
			ElePassID=0;
			Mchagino=0;
			Mneutralino=0;
			dRPhoEle=0;
			Invmass=0;
			MT_=0;
			ThreeBodyMass_=0;
			nJet = 0;
			HT = 0;

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
      nVtx = raw.nVtx;

      bool hasSignalPho(false), hasSignalEle(false),hasSignalNeu(false);
      bool matchRecoPho(false), matchRecoEle(false);
      float charginoMass(0), neutralinoMass(0);
      std::vector<mcData>::iterator signalPho;
      std::vector<mcData>::iterator signalEle;
      std::vector<mcData>::iterator signalNeu;
      std::vector<recoPhoton>::iterator recopho;
      std::vector<recoEle>::iterator recoele;
      for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){ 

          if(hasSignalPho && hasSignalEle && hasSignalNeu && charginoMass>0 && neutralinoMass>0)break; //if photon, electron and SUSY particle have been found, skip the loop
          //Look for neutralino
          if(itMC->getPID() == 1000022 && itMC->getMomPID()== 1000023)neutralinoMass = itMC->getmomMass();
          //Look for neutralino
          if(itMC->getPID()== 1000022 && fabs(itMC->getMomPID())== 1000024)charginoMass = itMC->getmomMass();

          //(((((((((((
          //Look for neutralino-decayed photon
          if(!hasSignalPho && itMC->isPhoton() && itMC->decayFromNeu()){
             hasSignalPho=true;
             signalPho = itMC;
             // Match reco photon to signal photon
             float mindR(1);
             int nphomatch(0);
	     	     for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
               float tmpdE = fabs(itpho->getEt() - signalPho->getEt())/signalPho->getEt();
               if(tmpdE < 0.1){
                  float tmpdR = DeltaR(itpho->getEta(), itpho->getPhi(), signalPho->getEta(), signalPho->getPhi());
                  if(tmpdR < mindR){mindR = tmpdR;recopho = itpho;}
                  if(mindR < 0.15)nphomatch+=1;
               }
             }
             if(mindR < 0.15)matchRecoPho=true;
             else matchRecoPho=false;              
           }//endif: found signal photon 
           //))))))))))))))))))))))
           
           if(abs(itMC->getPID()) == 12 && abs(itMC->getGMomPID())== 1000024){
              signalNeu = itMC; hasSignalNeu = true;
           }
           //(((((((((((((((((((((((        
           if(!hasSignalEle && itMC->isElectron() && abs(itMC->getMomPID())== 24){

             hasSignalEle=true;
             signalEle = itMC;
             // Match reco photon to signal photon
             float mindR(1);
             int nelematch(0);
	     	     for(std::vector<recoEle>::iterator itele = Ele.begin() ; itele != Ele.end(); ++itele){
               float tmpdE = fabs(itele->getPt() - signalEle->getPt())/signalEle->getPt();
               if(tmpdE < 0.1){
                  float tmpdR = DeltaR(itele->getEta(), itele->getPhi(), signalEle->getEta(), signalEle->getPhi());
                  if(tmpdR < mindR){mindR = tmpdR;recoele = itele;}
                  if(mindR < 0.15)nelematch+=1;
               }
             }
	     	     if(mindR < 0.15)matchRecoEle=true;
             else matchRecoEle=false;

           }//endif:: found signal ELectron
           //)))))))))))))))))))))))
      }//loop on MC particles

      sigMET = MET;     
      nJet = 0;
      HT = 0;
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(!itJet->passSignalSelection())continue;
				bool matchEG(false);
				if(hasSignalPho)
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)matchEG=true;	
				if(hasSignalEle)
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(),signalEle->getPhi()) <= 0.4)matchEG=true;
				if(matchEG)continue;
				nJet += 1;
				HT += itJet->getPt();
			}	
	  	Mchagino=neutralinoMass;
	  	Mneutralino=neutralinoMass;
      if(hasSignalPho){
				mcPhotonEt=signalPho->getEt();
				mcPhotonEta=signalPho->getEta();
				mcPhotonPhi=signalPho->getPhi();
      }
      if(matchRecoPho){ 
				recoPhotonEt=recopho->getEt();
				recoPhotonEta=recopho->getEta();
				recoPhotonPhi=recopho->getPhi();
				PhoR9=recopho->getR9();
				PhoHoverE=recopho->getHoverE();
				PhoSigma=recopho->getSigma();
				PhoChIso=recopho->getChIso();
				PhoNeuIso=recopho->getNeuIso();
				PhoPhoIso=recopho->getPhoIso();
				PhoPassID=recopho->isLoose();
      }
      if(hasSignalEle){
				mcElePt=signalEle->getPt();
				mcEleEta=signalEle->getEta();
				mcElePhi=signalEle->getPhi();
      }
      if(matchRecoEle){ 
				recoEleEt=recoele->getPt();
				recoEleEta=recoele->getEta();
				recoElePhi=recoele->getPhi();
				EleR9=recoele->getR9();
				EleHoverE=recoele->getHoverE();
				EleSigma=recoele->getSigma();
				EleChIso=recoele->getChIso();
				EleNeuIso=recoele->getNeuIso();
				ElePhoIso=recoele->getPhoIso();
				EleMiniIso=recoele->getMiniIso();
				EledEtaIn=recoele->getdEtaIn();
				EledPhiIn=recoele->getdPhiIn();
				EleD0=recoele->getD0();
				EleDz=recoele->getDz();
				EleooEmooP=recoele->getEoverPInv();
				ElePassID=recoele->isMiniMedium(); 
      }
      if(hasSignalPho && hasSignalEle){
	    	dRPhoEle=DeltaR(signalPho->getEta(), signalPho->getPhi(), signalEle->getEta(), signalEle->getPhi());
	    	Invmass=(signalPho->getP4() + signalEle->getP4()).M();
	    	MT_=sqrt(2*MET*signalEle->getPt()*(1-std::cos(DeltaR(0, signalEle->getPhi(), 0, METPhi))));
        ThreeBodyMass_=sqrt(2*MET*(signalPho->getP4()+signalEle->getP4()).Pt()*(1-std::cos(DeltaR(0, (signalPho->getP4()+signalEle->getP4()).Phi(), 0, METPhi))));
      }
      tree->Fill();


			if( fabs(Mchagino-750) < 5){
        n750 += 1;
        bool hasMuFromCh(false);
        bool hasEleFromCh(false);
        bool hasQFromCh(false);
        std::cout << "evt " << n750 << std::endl;     
				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
             std::cout << itMC->getPID() << "  " << itMC->getMomPID() << "  pt = " << itMC->getPt() << std::endl;
             if(fabs(itMC->getGMomPID())== 1000024){
                if(fabs(itMC->getPID()) == 11 && fabs(itMC->getMomPID()) == 24)hasEleFromCh=true;
                else if(fabs(itMC->getPID()) == 13 && fabs(itMC->getMomPID()) == 24)hasMuFromCh=true;
                else if(fabs(itMC->getPID()) <=6 && fabs(itMC->getMomPID()) == 24)hasQFromCh=true;
             }
				}//loop on MC particles
        if(hasMuFromCh)std::cout << "hasMuFromCh ";
        if(hasEleFromCh)std::cout << "hasEleFromCh ";
        if(hasQFromCh)std::cout << "hasQFromCh ";
        if(!hasMuFromCh && !hasEleFromCh && !hasQFromCh)std::cout << "NoDecayFromCh";
        std::cout << std::endl;
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
						if(!hasegPho){
							hasegPho=true;
							egsignalPho = itpho;
						}
						if(!hasmgPho){
							hasmgPho=true;
							mgsignalPho = itpho;
						}
					}
        }
        bool hasEle(false);
        std::vector<recoEle>::iterator egsignalEle = Ele.begin();
        if(hasegPho){
          for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
            if(hasEle)continue;
						if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
						if(itEle->passSignalSelection()){
							hasEle=true; 
							egsignalEle = itEle;
						}
          }
        }

        bool hasMu(false);
        std::vector<recoMuon>::iterator signalMu = Muon.begin();
        if(hasmgPho){
					for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
					if(hasMu)continue;
						if(itMu->passSignalSelection()){
							hasMu=true; 
							signalMu = itMu;
						}
					}
        }

			if( fabs(Mchagino-750) < 5){
         if(hasMu && hasmgPho){ 
					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
					float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));
          std::cout << "has mu+gamma :   dR = " << dRmg << "  MT= " << MT << "  MET= " << MET << std::endl;
         }
         else std::cout << "No mu+gamma" << std::endl;
      }




        if(hasegPho && hasEle){
          double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), egsignalEle->getEta(), egsignalEle->getPhi()); 
          if(dReg>0.8){
            	if(fabs((egsignalPho->getP4()+egsignalEle->getP4()).M() - 91.188) > 10.0){

								float deltaPhi = DeltaPhi(egsignalEle->getPhi(), METPhi);
								float MT = sqrt(2*MET*egsignalEle->getPt()*(1-std::cos(deltaPhi)));
        	    	eg_phoEt = egsignalPho->getCalibEt();
			    			eg_phoEta= egsignalPho->getEta();
                eg_phoPhi= egsignalPho->getPhi();
			    			eg_lepPt = egsignalEle->getCalibPt();
 			    			eg_lepEta= egsignalEle->getEta();
                eg_lepPhi= egsignalEle->getPhi();
 			    			eg_sigMT = MT;
			    			eg_sigMET= MET;
                eg_sigMETPhi = METPhi;
			    			eg_dPhiLepMET = deltaPhi; 
			    			eg_nVertex = nVtx; 
			    			eg_dRPhoLep= dReg;
                eg_invmass = (egsignalPho->getP4()+egsignalEle->getP4()).M();
								eg_sigMETJESup = MET_T1JESUp;
								eg_sigMETJESdo = MET_T1JESDo;
								eg_sigMETJERup = MET_T1JERUp;
								eg_sigMETJERdo = MET_T1JERDo;
								eg_dPhiLepMETJESup = DeltaPhi(egsignalEle->getPhi(), METPhi_T1JESUp);
								eg_dPhiLepMETJESdo = DeltaPhi(egsignalEle->getPhi(), METPhi_T1JESDo);
								eg_dPhiLepMETJERup = deltaPhi;
								eg_dPhiLepMETJERdo = deltaPhi;
								eg_sigMTJESup = sqrt(2*MET_T1JESUp*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESup)));
								eg_sigMTJESdo = sqrt(2*MET_T1JESDo*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESdo)));
								eg_sigMTJERup = sqrt(2*MET_T1JERUp*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERup)));
								eg_sigMTJERdo = sqrt(2*MET_T1JERDo*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERdo)));
							
								eg_nJet = 0;
								eg_HT = 0;
								eg_HTJESup = 0;
								eg_HTJESdo = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalPho->getEta(),egsignalPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalEle->getEta(),egsignalEle->getPhi()) <= 0.4)continue;
									eg_nJet += 1;
									eg_HT += itJet->getPt();
									eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
									eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
								}	
 

                egtree->Fill();
 
            }// Z mass Filter
          }//dR filter
        }// ele + pho candidate
       
 
       if(hasmgPho && hasMu){
					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
					if(dRmg>0.8){
						float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));
						float ThreeBodyMass = sqrt(2*MET*(mgsignalPho->getP4()+ signalMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (mgsignalPho->getP4()+signalMu->getP4()).Phi(), 0, METPhi))));

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
				 		mgtree->Fill();
         }//dR Filter
       }//Candidate Filter

    }//loop on entries


    outputfile->Write();
}
