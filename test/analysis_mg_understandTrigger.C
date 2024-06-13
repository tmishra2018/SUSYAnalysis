// g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_mg.C -o analysis_mg.exe `oot-config --libs`
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

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_jet.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"
bool apply_HEMveto=false;

void analysis_mg_understandTrigger(int RunYear, const char *Era){//main

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_mg()" << std::endl;
  logfile << "miniIso; one lepton for fakephoton background" << std::endl;
  RunType datatype;
  if(RunYear==2016) datatype = MuonEG2016;
  if(RunYear==2017) datatype = MuonEG2017;
  if(RunYear==2018) datatype = MuonEG2018;
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add(Form("/eos/uscms/store/user/tmishra/InputFilesDATA/MuonEG/MuonEG_%d%s.root",RunYear,Era));

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);


	TH1D *p_METFilter = new TH1D("p_METFilter","",12,-2,10);	

	int nBJet(0);
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
	float threeMass(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  int   nJet(0);
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
	sigtree->Branch("threeMass", &threeMass);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
  sigtree->Branch("nBJet",     &nBJet);
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

  int passHEM(0);
  std::cout << "Total evetns : " << nEvts << std::endl;
  logfile << "Total evetns : " << nEvts << std::endl;
	for (unsigned ievt(0); ievt<2000000; ++ievt){//loop on entries

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
			
			if(RunYear==2018 && !passHEMVeto(0,raw)) continue;
                        passHEM++;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nMu < 1 || raw.nPho <1)continue;

  //    int Nmedpho(0);
  //    for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
  //      if(!itpho->isMedium())continue;
  //      if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
  //      if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
  //      if(itpho->getSigma() < 0.005)continue;
  //      if(itpho->PixelSeed() != 0)continue;

  //      Nmedpho+=1;
  //    }
  //    if(Nmedpho >= 2)continue;

			nBJet = 0;
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(itJet->getPt() < 20)continue;
				if(itJet->isBJet())nBJet+=1;
			}
			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Photon                                              **********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			// ******* for jet-photon fake rate***************//
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  hadeleproxyPhoCollection;
			hadeleproxyPhoCollection.clear();

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				//if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

				// ******** very loose, before sigma and isolation cut.  For jet-photon fake rate, and hadron proxy ************//	
				if(GSFveto && PixelVeto && FSRVeto){
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
				if(GSFveto && PixelVeto && FSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

				if((!PixelVeto || !GSFveto)){
						proxyPhoCollection.push_back(itpho);
				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			bool hasTrail(false);
			std::vector<recoMuon>::iterator trailLep = Muon.begin();
			std::vector< std::vector<recoMuon>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoMuon>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			std::vector< std::vector<recoMuon>::iterator > miniisoLep;
			miniisoLep.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->isMedium() && itMu->getPt() > 15 && itMu->getMiniIso() < 0.2)miniisoLep.push_back(itMu);
				if(itMu->getPt() < 25)continue;
				if(itMu->isFakeProxy())fakeLepCollection.push_back(itMu);
				if(!itMu->passHLTSelection())continue;
				//if(itMu->isFakeProxy())fakeLepCollection.push_back(itMu);
				if(itMu->passSignalSelection()){
					if(proxyLepCollection.size() == 0)proxyLepCollection.push_back(itMu);
					//proxyLepCollection.push_back(itMu);
					if(hasLep && !hasTrail){
						hasTrail = true;
						trailLep = itMu;
					}
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itMu;
					}
				}
			}

			bool hasSigEle(false);
			bool hasDoubleEG(false);
			std::vector<recoEle>::iterator signalEle = Ele.begin();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if(itEle->getCalibPt() < 25)continue;
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
				if(!itEle->passHLTSelection())continue;
				if(!hasSigEle && itEle->passSignalSelection()){
					hasSigEle = true;
					signalEle = itEle;
					if(hasPho){
						double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
						if(dRlepphoton > 0.8 && ((signalPho->getCalibP4()+signalEle->getCalibP4()).M() - 91.188) > 10.0)hasDoubleEG = true;
					}
				}
			}
	
			if(hasDoubleEG)continue;

			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.8){
					npassdR+=1;
					p_METFilter->Fill(-2);
					p_METFilter->Fill(raw.failFilterStep(METFilter));	
					if(raw.passMETFilter(METFilter)){ 
						npassMETFilter +=1;
						if(fabs((signalPho->getCalibP4()+signalLep->getP4()).M() - 91.188) > 10.0)npassZ+=1;

						float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
						float ThreeBodyMass = sqrt(2*MET*(signalPho->getP4()+ signalLep->getP4()).Pt()*(1-std::cos(DeltaR(0, (signalPho->getP4()+signalLep->getP4()).Phi(), 0, METPhi))));

						phoEt = signalPho->getCalibEt();
						phoEta= signalPho->getEta();
						phoPhi= signalPho->getPhi();
						lepPt = signalLep->getPt();
						lepEta= signalLep->getEta();
						lepPhi= signalLep->getPhi();
						sigMT = MT;
						sigMET= MET;
						sigMETPhi = METPhi;
						dPhiLepMET = deltaPhi;
						threeMass = ThreeBodyMass;
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
				}//dR Filter
			}//Candidate Filter
		 
 
	}//loop on  events
cout<< sigtree->GetEntries()<<endl;
logfile << "sigTree events: " << sigtree->GetEntries() <<"; "<<100*sigtree->GetEntries()/nEvts<<"\%"<<std::endl;
cout<<"HLT passed : "<<npassHLT<<endl;

p_eventcount->Fill("Total",nTotal);
p_eventcount->Fill("passHLT",npassHLT);
p_eventcount->Fill("passPho",npassPho);
p_eventcount->Fill("passMuon",npassLep);
p_eventcount->Fill("passdR",npassdR);
p_eventcount->Fill("passMETFilter",npassMETFilter);
p_eventcount->Fill("passZ",npassZ);
	if(RunYear==2018) logfile << "pass HEM cut:  " << passHEM*100/nEvts<<endl;
}

int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_mg_understandTrigger(atoi(argv[1]),argv[2]);
    return 0;
}
