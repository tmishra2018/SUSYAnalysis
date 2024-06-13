// g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_egMC.C -o analysis_egMC.exe `root-config --libs`

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
#include "TRandom3.h"
#include "TFileCollection.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_jet.h"
#include "../include/analysis_tools.h"
#include "../src/analysis_rawData.cc"
#include "../src/analysis_ele.cc"
#include "../src/analysis_photon.cc"


void analysis_egMC(int RunYear, const char *Era, const char *Sample){//main  
 
  RunType datatype;
  if(RunYear==2016) datatype = DoubleEG2016;
  if(RunYear==2017) datatype = DoubleEG2017;
  if(RunYear==2018) datatype = DoubleEG2018;

  TChain* es = new TChain("ggNtuplizer/EventTree");
  char* inputfile = new char[300];
  sprintf(inputfile,"/eos/uscms/store/group/lpcsusyhad/Tribeni/%s/%s_%d%s.root",Sample,Sample,RunYear,Era); 
  es->Add(inputfile);
  TFile* outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/egMC/analysis_egMC_%s_%d%s.root",Sample,RunYear,Era),"RECREATE");
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
  else if(strstr(inputfile, "TTGJets") != NULL){
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
  else if(strstr(inputfile, "T5Wg_mG-1800_mLSP-800") != NULL){
		std::cout << "T5Wg_mG-1800_mLSP-800 sample !" << std::endl;
  		mcType = MCType::T5Wg_mG-1800_mLSP800; 
  }
  else if(strstr(inputfile, "T5Wg_mG-1800_mLSP-1600") != NULL){
		std::cout << "T5Wg_mG-1800_mLSP-1600 sample !" << std::endl;
  		mcType = MCType::T5Wg_mG-1800_mLSP1600; 
  }
  else if(strstr(inputfile, "T5Wg_mG-1500_mLSP-1") != NULL){
		std::cout << "T5Wg_mG-1500_mLSP-1 sample !" << std::endl;
  		mcType = MCType::T5Wg_mG-1500_mLSP1; 
  }
  else if(strstr(inputfile, "SMS-TChiWg_mChi-1000_mLSP-1") != NULL){
		std::cout << "SMS-TChiWg_mChi-1000_mLSP-1 sample !" << std::endl;
  		mcType = MCType::TChiWg_mChi-1000_mLSP1; 
  }
  else {
		std::cout << "not specific MC !" << std::endl;
		mcType = MCType::NOMC;
  }
  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();
  float lumiWeight = getEvtWeight(mcType, RunYear,ntotalevent);
  
  cout<<"crosssection = "<<crosssection<<endl;
  cout<<"ntotalevent = "<<ntotalevent<<endl;
  cout<<"lumiWeight = "<<lumiWeight<<endl;

  TH1F *p_eventcount = new TH1F("p_eventcount","p_eventcount",7,0,7);
  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);
  TTree *newtree = (TTree*)es->Clone();
  rawData raw(es, datatype); 
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  float MET(0);
  float METPhi(0);
  int   nVtx(0);
  int METFilter(0);
  const unsigned nEvts = newtree->GetEntries(); 
  cout << "Total event: " << nEvts << std::endl;

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
  int nBJet(0);
  float eleR9(0);
  float phoR9(0);
  float MiniIso(0);

  sigtree->Branch("crosssection",       &crosssection);
  sigtree->Branch("ntotalevent",       &ntotalevent);
  sigtree->Branch("lumiWeight",       &lumiWeight);
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
  sigtree->Branch("nBJet",     &nBJet);
  sigtree->Branch("eleR9",     &eleR9);
  sigtree->Branch("phoR9",     &phoR9);
  sigtree->Branch("MiniIso",     &MiniIso);
   
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      	if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
        raw.GetData(es, ievt);
	Photon.clear();
        Muon.clear();
        Ele.clear();
        JetCollection.clear();
	for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
        for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
	for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}

        MET = raw.pfMET;
        METPhi = raw.pfMETPhi;
        nVtx = raw.nVtx;
	METFilter = raw.metFilters;
	nVtx = raw.nVtx;
	run=raw.run;
	event=raw.event;
	lumis=raw.lumis;
	lumiWeight = getEvtWeight(mcType, RunYear,ntotalevent);
	nTotal+=1;
	if(!raw.passHLT())continue;
	if(raw.nGoodVtx < 1)continue;
	npassHLT+=1;
	if(raw.nEle < 1 || raw.nPho < 1)continue;      


	nBJet = 0;
	for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
			if(itJet->getPt() < 20)continue;
			if(itJet->isBJet())nBJet+=1;
	}

	bool hasPho(false);
	std::vector<recoPhoton>::iterator signalPho = Photon.begin();
	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
		if(!itpho->isEB()) continue;
		if(itpho->getR9() < 0.5)continue;
		if(!itpho->passHLTSelection())continue;
		if(!itpho->passSignalSelection())continue;
		bool PixelVeto = itpho->PixelSeed()==0? true: false;
		bool GSFveto(true);
		bool FSRVeto(true);
		for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
			if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
			if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
		}
		for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
			if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;
		
		if(GSFveto && PixelVeto && FSRVeto){
					// if not found a photon
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
		}
	}
	
	bool hasLep(false);
	std::vector<recoEle>::iterator signalLep = Ele.begin();
	for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
			if(!itEle->passHLTSelection())continue;
			if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
			if(!itEle->passSignalSelection())continue;
			// if not found a lep
			if(!hasLep){
				hasLep = true;
				npassLep +=1;
				signalLep = itEle;
			}
	}	



	if(hasPho && hasLep){
		double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
		if(dRlepphoton > 0.8){
			npassdR+=1;
			if(fabs((signalPho->getCalibP4()+signalLep->getP4()).M() - 91.188) > 10.0){
				npassZ+=1;
				if(raw.passMETFilter(METFilter)){
					npassMETFilter +=1;

					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getCalibPt()*(1-std::cos(deltaPhi)));
					phoEt= signalPho->getCalibEt();
					phoEta= signalPho->getEta();
					phoPhi= signalPho->getPhi();
					phoR9= signalPho->getR9();

					lepPt= signalLep->getCalibPt();
					lepEta= signalLep->getEta();
					lepPhi= signalLep->getPhi();
					eleR9= signalLep->getR9();
					MiniIso= signalLep->getMiniIso();
					sigMT= MT;
					sigMET= MET;
					sigMETPhi = METPhi;
					dPhiLepMET = deltaPhi; 
					nVertex = nVtx; 
					dRPhoLep= dRlepphoton;
					nJet = 0;
					HT = 0;
					for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
						if(!itJet->passSignalSelection())continue;
						if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
						if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
						nJet += 1;
						HT += itJet->getPt();
					}
					sigtree->Fill();
				} // MET Filter
			} //  Z mass Filter
		}  //dR filter
	}// ele + pho candidate


  }//loop on  events
  cout<<" nTotal : "<<nTotal<<endl;	
  cout<<" npassHLT : "<<npassHLT<<endl;	
  cout<<" npassPho : "<<npassPho<<endl;	
  cout<<" npassLep : "<<npassLep<<endl;	
  cout<<" npassdR : "<<npassdR<<endl;	
  cout<<" npassZ : "<<npassZ<<endl;	
  cout<<" npassMETFilter : "<<npassMETFilter<<endl;	

  p_eventcount->Fill(0.5, nTotal);
  p_eventcount->Fill(1.5, npassHLT);
  p_eventcount->Fill(2.5, npassPho);
  p_eventcount->Fill(3.5, npassLep);
  p_eventcount->Fill(4.5, npassdR);
  p_eventcount->Fill(5.5, npassZ);
  p_eventcount->Fill(6.5, npassMETFilter);
  outputfile->Write();
  outputfile->Close();
}

int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_egMC(atoi(argv[1]),argv[2],argv[3]);
    return 0;
}
