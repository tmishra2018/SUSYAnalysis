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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

void validateProxyLep(){//main 
	float crosssection(1);
	float totalEvent(1);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/resTree_egsignal_GJets_EMall.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/resTree_egsignal_GJets_Pt-EM.log"); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
  //es->Add("root://cmseos.fnal.gov//store/user/msun/MC/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root");

  //const unsigned nEvts = es->GetEntries(); 
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
	int   lepPixel(0);
	float lepRelIso(0);
  float miniIso(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	int   iPho(0);
  std::vector<int> mcPID;
  std::vector<float> mcEta;
  std::vector<float> mcPhi;
  std::vector<float> mcPt;
  std::vector<int> mcMomPID;
 
  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("lumis",     &lumis);
  sigtree->Branch("phoEt",     &phoEt);
  sigtree->Branch("phoEta",    &phoEta);
  sigtree->Branch("phoPhi",    &phoPhi);
  sigtree->Branch("lepPt",     &lepPt);
  sigtree->Branch("lepEta",    &lepEta);
  sigtree->Branch("lepPhi",    &lepPhi);
	sigtree->Branch("lepPixel",  &lepPixel);
	sigtree->Branch("lepRelIso", &lepRelIso);
	sigtree->Branch("miniIso",   &miniIso);
  sigtree->Branch("sigMT",     &sigMT);
  sigtree->Branch("sigMET",    &sigMET);
  sigtree->Branch("sigMETPhi", &sigMETPhi);
  sigtree->Branch("dPhiLepMET",&dPhiLepMET);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
  sigtree->Branch("mcPID",           &mcPID);
  sigtree->Branch("mcEta",           &mcEta);
  sigtree->Branch("mcPhi",           &mcPhi);
  sigtree->Branch("mcPt",            &mcPt);
  sigtree->Branch("mcMomPID",        &mcMomPID);
	sigtree->Branch("crosssection",    &crosssection);
	sigtree->Branch("totalEvent",      &totalEvent);
	sigtree->Branch("iPho",            &iPho);

//*********** fake lepton *********************//
  TTree *fakeEtree = new TTree("fakeETree","fakeETree");
  float fakeEphoEt(0);
  float fakeEphoEta(0);
  float fakeEphoPhi(0);
  float fakeElepPt(0);
  float fakeElepEta(0);
  float fakeElepPhi(0);
  float fakeEsigMT(0);
  float fakeEsigMET(0);
  float fakeEsigMETPhi(0);
  float fakeEdPhiLepMET(0);
  int   fakeEnVertex(0);
  float fakeEdRPhoLep(0);
  float fakeEHT(0);
  float fakeEnJet(0);
  std::vector<int> fakeEmcPID;
  std::vector<float> fakeEmcEta;
  std::vector<float> fakeEmcPhi;
  std::vector<float> fakeEmcPt;
  std::vector<int> fakeEmcMomPID;
  
  fakeEtree->Branch("phoEt",     &fakeEphoEt);
  fakeEtree->Branch("phoEta",    &fakeEphoEta);
  fakeEtree->Branch("phoPhi",    &fakeEphoPhi);
  fakeEtree->Branch("lepPt",     &fakeElepPt);
  fakeEtree->Branch("lepEta",    &fakeElepEta);
  fakeEtree->Branch("lepPhi",    &fakeElepPhi);
  fakeEtree->Branch("sigMT",     &fakeEsigMT);
  fakeEtree->Branch("sigMET",    &fakeEsigMET);
  fakeEtree->Branch("sigMETPhi", &fakeEsigMETPhi);
  fakeEtree->Branch("dPhiLepMET",&fakeEdPhiLepMET);
  fakeEtree->Branch("nVertex",   &fakeEnVertex);
  fakeEtree->Branch("dRPhoLep",  &fakeEdRPhoLep);
  fakeEtree->Branch("HT",        &fakeEHT);
  fakeEtree->Branch("nJet",      &fakeEnJet);
  fakeEtree->Branch("mcPID",           &fakeEmcPID);
  fakeEtree->Branch("mcEta",           &fakeEmcEta);
  fakeEtree->Branch("mcPhi",           &fakeEmcPhi);
  fakeEtree->Branch("mcPt",            &fakeEmcPt);
  fakeEtree->Branch("mcMomPID",        &fakeEmcMomPID);
	fakeEtree->Branch("crosssection",    &crosssection);
	fakeEtree->Branch("totalEvent",      &totalEvent);
	fakeEtree->Branch("iPho",            &iPho);

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
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  float proxynJet(0);
  
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
  proxytree->Branch("nVertex",   &proxynVertex);
  proxytree->Branch("dRPhoLep",  &proxydRPhoLep);
  proxytree->Branch("HT",        &proxyHT);
  proxytree->Branch("nJet",      &proxynJet);
	proxytree->Branch("crosssection",    &crosssection);
	proxytree->Branch("totalEvent",      &totalEvent);

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
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
			METFilter = raw.metFilters;
			nVtx = raw.nVtx;
			jetNumber = raw.nJet;
			run=raw.run;
			event=raw.event;
			lumis=raw.lumis;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(((raw.HLTPho >> 14) &1) ==0)continue;
			npassHLT+=1;

			if(raw.nEle < 1 || raw.nPho <1)continue;

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			iPho = 0;
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(!itpho->isEB())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.03)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.03 && ie->getEt()>2.0)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

				if(itpho->getChIso()<20 && itpho->getSigma()< 0.02){
					if(!passSigma || !passChIso){
						if(GSFveto && PixelVeto && FSRVeto)jetPhoCollection.push_back(itpho);
					}
				}

				if(!itpho->passSignalSelection())continue;
				if(!itpho->fireDoubleTrg(5) && !itpho->fireDoubleTrg(6))continue;
				if(GSFveto && PixelVeto && FSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
						iPho = std::distance( Photon.begin(), itpho);
					}
				}

				if(FSRVeto && (!PixelVeto || !GSFveto))proxyPhoCollection.push_back(itpho);
			}

			bool hasEle(false);
			std::vector<recoEle>::iterator signalEle = Ele.begin();
			std::vector< std::vector<recoEle>::iterator > proxyEleCollection;
			proxyEleCollection.clear();
			std::vector< std::vector<recoEle>::iterator > fakeEleCollection;
			fakeEleCollection.clear();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if(!itEle->fireTrgs(21) && !itEle->fireTrgs(22))continue;
				if(itEle->getCalibPt() < 25)continue;
				if(itEle->isFakeProxy())fakeEleCollection.push_back(itEle);	
				if(itEle->passSignalSelection()){
				//if(itEle->passBasicID() && itEle->getMiniIso() < 0.4){
					proxyEleCollection.push_back(itEle);
					if(!hasEle){
						for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
							if(DeltaR(itpho->getEta(), itpho->getPhi(), itEle->getEta(), itEle->getPhi()) < 0.05){
								if(itpho->PixelSeed()==0)lepPixel = 2;
								else lepPixel = 1;
							}
						}
						hasEle=true; 
						npassLep +=1;
						signalEle = itEle;
					}
				}
			}

			bool saveMC(false);
			 
			if(hasPho && hasEle){
				double dReg = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
				if(dReg > 0.8){
					npassdR+=1;
					if(fabs((signalPho->getCalibP4()+signalEle->getCalibP4()).M() - 91.188) > 10.0){

						npassZ+=1;
						if(METFilter == 0){
							npassMETFilter +=1;

							float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalEle->getCalibPt()*(1-std::cos(deltaPhi)));

							phoEt = signalPho->getCalibEt();
							phoEta= signalPho->getEta();
							phoPhi= signalPho->getPhi();
							lepPt = signalEle->getCalibPt();
							lepEta= signalEle->getEta();
							lepPhi= signalEle->getPhi();
							lepRelIso = signalEle->getRelIso();
							miniIso = signalEle->getMiniIso();
							sigMT = MT;
							sigMET= MET;
							sigMETPhi = METPhi;
							dPhiLepMET = deltaPhi; 
							nVertex = nVtx; 
							dRPhoLep= dReg;
							nJet = jetNumber; 

							mcPID.clear();
							mcEta.clear();
							mcPhi.clear();
							mcPt.clear();
							mcMomPID.clear();
							for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
								if(itMC->getEt() < 1.0)continue;
								float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
								float mcdR2 = DeltaR(signalEle->getEta(),signalEle->getPhi(), itMC->getEta(), itMC->getPhi());
								if(mcdR < 0.3 || mcdR2<0.3){
									mcPID.push_back(itMC->getPID());
									mcMomPID.push_back(itMC->getMomPID());
									mcEta.push_back(itMC->getEta());      
									mcPhi.push_back(itMC->getPhi());
									mcPt.push_back(itMC->getEt());
								}
							}
							sigtree->Fill();

							if(datatype == MC)saveMC=true;
						}//MET Filter
					}// Z mass Filter
				}//dR filter
			}// ele + pho candidate
	 

			if(hasPho && !hasEle){
				std::vector<recoPhoton>::iterator fakeEPho = signalPho; 
				for(unsigned ip(0); ip < fakeEleCollection.size(); ip++){
					std::vector<recoEle>::iterator fakeEle = fakeEleCollection[ip];
					double dReg = DeltaR(fakeEPho->getEta(), fakeEPho->getPhi(), fakeEle->getEta(), fakeEle->getPhi());
					if(dReg>0.8){
						if(fabs((fakeEPho->getCalibP4()+fakeEle->getCalibP4()).M() - 91.188) > 10.0){
							if(METFilter == 0){

								float fakeE_deltaPhi = DeltaPhi(fakeEle->getPhi(), METPhi);
								float fakeE_MT = sqrt(2*MET*fakeEle->getCalibPt()*(1-std::cos(fakeE_deltaPhi)));
								fakeEphoEt = fakeEPho->getCalibEt();
								fakeEphoEta= fakeEPho->getEta();
								fakeEphoPhi= fakeEPho->getPhi();
								fakeElepPt = fakeEle->getCalibPt();
								fakeElepEta= fakeEle->getEta();
								fakeElepPhi= fakeEle->getPhi();
								fakeEsigMT = fakeE_MT;
								fakeEsigMET= MET;
								fakeEsigMETPhi = METPhi;
								fakeEdPhiLepMET = fakeE_deltaPhi; 
								fakeEnVertex = nVtx; 
								fakeEdRPhoLep= dReg;
								fakeEnJet = jetNumber; 
								fakeEmcPID.clear();
								fakeEmcEta.clear();
								fakeEmcPhi.clear();
								fakeEmcPt.clear();
								fakeEmcMomPID.clear();
								for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
									if(itMC->getEt() < 1.0)continue;
									float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
									float mcdR2 = DeltaR(signalEle->getEta(),signalEle->getPhi(), itMC->getEta(), itMC->getPhi());
									if(mcdR < 0.3 || mcdR2<0.3){
										fakeEmcPID.push_back(itMC->getPID());
										fakeEmcMomPID.push_back(itMC->getMomPID());
										fakeEmcEta.push_back(itMC->getEta());      
										fakeEmcPhi.push_back(itMC->getPhi());
										fakeEmcPt.push_back(itMC->getEt());
									}
								}
								fakeEtree->Fill();

								if(datatype == MC)saveMC=true;
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}

			for(unsigned ip(0); ip < proxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyEleCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = proxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyEleCollection[ie];
					double dReg = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dReg>0.8){
						if(fabs((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0){
							if(METFilter == 0){
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
								proxydRPhoLep= dReg;
								proxynJet = jetNumber; 
								proxytree->Fill();

								if(datatype == MC)saveMC=true;
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection

	}//loop on  events

	outputfile->Write();
	outputfile->Close();
	logfile.close();
}


