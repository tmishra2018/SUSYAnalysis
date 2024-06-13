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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_tools.h"

void analysis_AbsEff(){//main 
	
  gSystem->Load("../../lib/libAnaClasses.so");

  char outputname[100] = "/eos/uscms/store/user/tmishra/Trigger/Efficiency_TT.root";
  RunType datatype(MC);
  TChain* es = new TChain("ggNtuplizer/EventTree");
        es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/TTJets/TTJets_2016.root");
	//es->Add("/eos/uscms/store/user/msun/copied/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root");
 
  //const unsigned nEvts = es->GetEntries(); 
  const unsigned nEvts = 100000; 

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

  int mcType = MCType::TT;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 

  float ntotalevent = es->GetEntries();
//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float eg_lepPt(0);
  float eg_lepEta(0);
  float eg_lepPhi(0);
  int   eg_nVertex(0);
  int   eg_mcPID;
  float eg_mcEta;
  float eg_mcPhi;
  float eg_mcPt;
  int   eg_mcMomPID;
  int   eg_mcGMomPID;
  int   eg_mcStatus;

  egtree->Branch("lepPt",     &eg_lepPt);
  egtree->Branch("lepEta",    &eg_lepEta);
  egtree->Branch("lepPhi",    &eg_lepPhi);
  egtree->Branch("nVertex",   &eg_nVertex);
  egtree->Branch("mcType",   &mcType);
  egtree->Branch("mcPID",    &eg_mcPID);
  egtree->Branch("mcEta",    &eg_mcEta);
  egtree->Branch("mcPhi",    &eg_mcPhi);
  egtree->Branch("mcPt",     &eg_mcPt);
  egtree->Branch("mcMomPID", &eg_mcMomPID);
  egtree->Branch("mcGMomPID", &eg_mcGMomPID);
  egtree->Branch("mcStatus", &eg_mcStatus);

  TTree *mgtree = new TTree("mgTree","mgTree");
  float mg_lepPt(0);
  float mg_lepEta(0);
  float mg_lepPhi(0);
  int   mg_nVertex(0);
  int   mg_mcPID;
  float mg_mcEta;
  float mg_mcPhi;
  float mg_mcPt;
  int   mg_mcMomPID;
  int   mg_mcGMomPID;
	int   mg_mcStatus;
 
  mgtree->Branch("lepPt",     &mg_lepPt);
  mgtree->Branch("lepEta",    &mg_lepEta);
  mgtree->Branch("lepPhi",    &mg_lepPhi);
  mgtree->Branch("nVertex",   &mg_nVertex);
  mgtree->Branch("mcType",   &mcType);
  mgtree->Branch("mcPID",    &mg_mcPID);
  mgtree->Branch("mcEta",    &mg_mcEta);
  mgtree->Branch("mcPhi",    &mg_mcPhi);
  mgtree->Branch("mcPt",     &mg_mcPt);
  mgtree->Branch("mcMomPID", &mg_mcMomPID);
  mgtree->Branch("mcGMomPID", &mg_mcGMomPID);
  mgtree->Branch("mcStatus", &mg_mcStatus);

//*********** histo list **********************//
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
/*********************************************/
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
        int nVtx = raw.nVtx;

				bool hasMC_ele(false);
				bool hasMC_mu(false);
				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					if(itMC->getEt() < 5.0)continue;
					if(fabs(itMC->getPID()) == 11 && fabs(itMC->getMomPID()) == 24){
 						eg_mcPID = itMC->getPID();
 						eg_mcEta = itMC->getEta();
 						eg_mcPhi = itMC->getPhi();
 						eg_mcPt  = itMC->getPt();
 						eg_mcMomPID = itMC->getMomPID();
 						eg_mcGMomPID = itMC->getGMomPID();

       			bool hasEle(false);
       			std::vector<recoEle>::iterator signalEle = Ele.begin();
       			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
       			  if(hasEle)break;
			 				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
			 				if(DeltaR(eg_mcEta, eg_mcPhi, itEle->getEta(), itEle->getPhi()) < 0.5 &&  itEle->isMiniMedium()){
			 					hasEle=true; 
			 					signalEle = itEle;
			 				}
       			}
 	    		  if(hasEle){ 
 	 				    eg_lepPt = signalEle->getPt();
 							eg_lepEta= signalEle->getEta();
 	    		    eg_lepPhi= signalEle->getPhi();
 	 				    eg_nVertex = nVtx; 
						}
						else{
 	 				    eg_lepPt = 0; 
 							eg_lepEta= 0;
 	    		    eg_lepPhi= 0;
 	 				    eg_nVertex = 0;
						}
 	    		  egtree->Fill();
					}
					if(fabs(itMC->getPID()) == 13 && fabs(itMC->getMomPID()) == 24){
 						mg_mcPID = itMC->getPID();
 						mg_mcEta = itMC->getEta();
 						mg_mcPhi = itMC->getPhi();
 						mg_mcPt  = itMC->getPt();
 						mg_mcMomPID = itMC->getMomPID();
 						mg_mcGMomPID = itMC->getGMomPID();

  			    bool hasMu(false);
  			    std::vector<recoMuon>::iterator signalMu = Muon.begin();
							for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
							if(hasMu)break;
							if(DeltaR(mg_mcEta, mg_mcPhi, itMu->getEta(), itMu->getPhi()) < 0.5 && itMu->passSignalSelection()){
								hasMu=true; 
								signalMu = itMu;
							}
						}
						if(hasMu){
	  					mg_lepPt = signalMu->getPt();
	  					mg_lepEta= signalMu->getEta();
	  					mg_lepPhi= signalMu->getPhi();
						}
						else{
							mg_lepPt = 0;
							mg_lepEta= 0;
							mg_lepPhi= 0;
						}
	  		 		mgtree->Fill();
					}	
				}
 
///*************************   eg filters *****************************//
       
//**********************  mg filter **************************************//         
 
	}//loop on  events


outputfile->Write();
outputfile->Close();

}
