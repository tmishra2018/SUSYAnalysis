// added some branches scalefactor, else rest is same
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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_scalefactor.h"

void analysis_VGamma(){//main 
	float crosssection = 3.772;
	float ntotalevent  = 9885348;
	float PUweight(1);
	
  gSystem->Load("../../../lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/resTree_VGamma_TTGJets.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/resTree_VGamma_TTGJets.log"); 

  logfile << "analysis_VGamma()" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTGJets_RunIISummer16MiniAODv2-TrancheIV_v6_ext1.root");
  logfile << "/store/user/msun/MCSummer16/TTGJets_RunIISummer16MiniAODv2-TrancheIV_v6_ext1.root" << std::endl;
    
  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

  int mcType = MCType::ZG;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;

//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_lepPt(0);
  float eg_lepEta(0);
  float eg_lepPhi(0);
	float eg_scalefactor(0);
	float eg_scalefactorup(0);
	float eg_scalefactordo(0);
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
	egtree->Branch("scalefactor", &eg_scalefactor);
	egtree->Branch("scalefactorup", &eg_scalefactorup);
	egtree->Branch("scalefactordo", &eg_scalefactordo);
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
	float mg_scalefactor(0);
	float mg_scalefactorup(0);
	float mg_scalefactordo(0);
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
	mgtree->Branch("scalefactor",   &mg_scalefactor);
	mgtree->Branch("scalefactorup", &mg_scalefactorup);
	mgtree->Branch("scalefactordo", &mg_scalefactordo);
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
	esfScaleFactor  objectESF;
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
  logfile << "RunType: " << datatype << std::endl;

  std::cout << "Total evetns : " << nEvts << std::endl;

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;

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
        jetNumber = raw.nJet;
				 
				if(RunYear==2016 and preVFP==1)PUweight = getPUESF16preVFP(nVtx);
				else if(RunYear==2016 and preVFP==0)PUweight = getPUESF16(nVtx);
                        	else if(RunYear==2017)PUweight = getPUESF17(nVtx);
                        	else if(RunYear==2018)PUweight = getPUESF18(nVtx);

        if(raw.nPho <1)continue;

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
					 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getEt()>2.0)FSRVeto=false;
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

        if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){
          double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
          if(dReg>0.8){
            if(fabs((egsignalPho->getP4()+signalEle->getP4()).M() - 91.188) > 10.0){
              if(METFilter == 0){

								float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
								float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
        	    	eg_phoEt = egsignalPho->getEt();
			    			eg_phoEta= egsignalPho->getEta();
                eg_phoPhi= egsignalPho->getPhi();
			    			eg_lepPt = signalEle->getPt();
 			    			eg_lepEta= signalEle->getEta();
                eg_lepPhi= signalEle->getPhi();
								eg_scalefactor = objectESF.getElectronESF(signalEle->getEt(), signalEle->getEta())*objectESF.getPhotonESF(egsignalPho->getEt(), egsignalPho->getEta()); 
								eg_scalefactorup = eg_scalefactor + objectESF.getElectronESFError(signalEle->getEt(), signalEle->getEta())+ objectESF.getPhotonESFError(egsignalPho->getEt(), egsignalPho->getEta());
								eg_scalefactordo = eg_scalefactor - objectESF.getElectronESFError(signalEle->getEt(), signalEle->getEta())- objectESF.getPhotonESFError(egsignalPho->getEt(), egsignalPho->getEta());
 			    			eg_sigMT = MT;
			    			eg_sigMET= MET;
                eg_sigMETPhi = METPhi;
			    			eg_dPhiLepMET = deltaPhi; 
			    			eg_nVertex = nVtx; 
			    			eg_dRPhoLep= dReg;
			    			eg_nJet = jetNumber; 
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
 
								if(datatype == MC){
									eg_mcPID.clear();
									eg_mcEta.clear();
									eg_mcPhi.clear();
									eg_mcPt.clear();
									eg_mcMomPID.clear();
									eg_mcGMomPID.clear();
									for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
										if(itMC->getEt() < 1.0)continue;
										float mcdR = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
										float mcdRe = DeltaR(signalEle->getEta(),signalEle->getPhi(), itMC->getEta(), itMC->getPhi());
										if(mcdR < 0.3 || mcdRe < 0.3){
											eg_mcPID.push_back(itMC->getPID());
											eg_mcMomPID.push_back(itMC->getMomPID());
											eg_mcEta.push_back(itMC->getEta());      
											eg_mcPhi.push_back(itMC->getPhi());
											eg_mcPt.push_back(itMC->getEt());
											eg_mcGMomPID.push_back(itMC->getGMomPID());
										}
									}
								}

                egtree->Fill();
 
              }//MET Filter
            }// Z mass Filter
          }//dR filter
        }// ele + pho candidate
       
         
        if(hasmgPho && hasMu && (((raw.HLTEleMuX >> 8)&1)!=0 || ((raw.HLTEleMuX >> 41)&1)!=0)){
					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
					if(dRmg>0.8){
						if(METFilter == 0){ 

						float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));

						mg_phoEt = mgsignalPho->getEt();
						mg_phoEta= mgsignalPho->getEta();
						mg_phoPhi= mgsignalPho->getPhi();
						mg_lepPt = signalMu->getPt();
						mg_lepEta= signalMu->getEta();
						mg_lepPhi= signalMu->getPhi();
						mg_scalefactor = objectESF.getMuonESF(signalMu->getEt(), signalMu->getEta())*objectESF.getPhotonESF(mgsignalPho->getEt(), mgsignalPho->getEta()); 
						mg_scalefactorup = mg_scalefactor + objectESF.getMuonESFError(signalMu->getEt(), signalMu->getEta())+ objectESF.getPhotonESFError(mgsignalPho->getEt(), mgsignalPho->getEta());
						mg_scalefactordo = mg_scalefactor - objectESF.getMuonESFError(signalMu->getEt(), signalMu->getEta())- objectESF.getPhotonESFError(mgsignalPho->getEt(), mgsignalPho->getEta());
						mg_sigMT = MT;
						mg_sigMET= MET;
						mg_sigMETPhi = METPhi;
						mg_dPhiLepMET = deltaPhi;
						mg_nVertex = nVtx;
						mg_dRPhoLep= dRmg;
						mg_nJet = jetNumber;
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

						if(datatype == MC){
						 mg_mcPID.clear();
						 mg_mcEta.clear();
						 mg_mcPhi.clear();
						 mg_mcPt.clear();
						 mg_mcMomPID.clear();
						 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	if(itMC->getEt() < 1.0)continue;
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
logfile.close();

}


