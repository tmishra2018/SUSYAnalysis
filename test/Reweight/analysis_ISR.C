// g++ `root-config --cflags` analysis_ISR.C -o analysis_ISR.exe `root-config --libs`
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
#include "../../include/analysis_tools.h"
#include "../../include/analysis_mcData.h"
#include "../../src/analysis_rawData.cc"
#include "../../src/analysis_ele.cc"
#include "../../src/analysis_muon.cc"
#include "../../src/analysis_photon.cc"


bool passMETFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 9; im++){
	if(im != 8){ 
    		if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}
	}
  }
  return passfilter;
}

void analysis_ISR(int RunYear, const char *Era){//main  

  gSystem->Load("../../lib/libAnaClasses.so");

  RunType datatype;
  if(RunYear==2016)    datatype = MuonEG2016;
  if(RunYear==2017)    datatype = MuonEG2017;
  if(RunYear==2018)    datatype = MuonEG2018;

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_%d%s.log",RunYear,Era)); 
  logfile << "analysis_mg()" << std::endl;

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add(Form("/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_%d%s.root",RunYear,Era));

  const unsigned nEvts = es->GetEntries(); 
  float MCweight = 1;	
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile;
  outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_%d%s.root",RunYear,Era),"RECREATE");
  outputfile->cd();

  //************ ZGamma Tree **********************//
  TTree *Ztree = new TTree("ZTree","ZTree");
  float ZphoEt(0);
  float ZphoEta(0);
  float ZphoPhi(0);
  float ZlepPt(0);
  float ZlepEta(0);
  float ZlepPhi(0);
  float ZtrailPt(0);
  float ZtrailEta(0);
  float ZtrailPhi(0);
  float ZsigMT(0);
  float ZsigMET(0);
  float ZsigMETPhi(0);
  float ZdPhiLepMET(0);
  float ZthreeMass(0);
  float ZdilepMass(0);
  int   ZnVertex(0);
  float ZdRPhoLep(0);
  float ZHT(0);
  float ZnJet(0);
  float ZISRJetPt(0);
 
  Ztree->Branch("MCweight",  &MCweight); 
  Ztree->Branch("phoEt",     &ZphoEt);
  Ztree->Branch("phoEta",    &ZphoEta);
  Ztree->Branch("phoPhi",    &ZphoPhi);
  Ztree->Branch("lepPt",     &ZlepPt);
  Ztree->Branch("lepEta",    &ZlepEta);
  Ztree->Branch("lepPhi",    &ZlepPhi);
  Ztree->Branch("trailPt",   &ZtrailPt);
  Ztree->Branch("trailEta",  &ZtrailEta);
  Ztree->Branch("trailPhi",  &ZtrailPhi);
  Ztree->Branch("sigMT",     &ZsigMT);
  Ztree->Branch("sigMET",    &ZsigMET);
  Ztree->Branch("sigMETPhi", &ZsigMETPhi);
  Ztree->Branch("dPhiLepMET",&ZdPhiLepMET);
  Ztree->Branch("threeMass", &ZthreeMass);
  Ztree->Branch("dilepMass", &ZdilepMass);
  Ztree->Branch("nVertex",   &ZnVertex);
  Ztree->Branch("dRPhoLep",  &ZdRPhoLep);
  Ztree->Branch("HT",        &ZHT);
  Ztree->Branch("nJet",      &ZnJet);
  Ztree->Branch("ISRJetPt",     &ZISRJetPt);

//*************** for jet-photon fake rate ***********************//
  TTree *hadrontree = new TTree("hadronTree","hadronTree");
  float hadron_phoEt(0);
  float hadron_phoEta(0);
  float hadron_phoPhi(0);
  float hadron_phoSigma(0);
  float hadron_phoChIso(0);

  hadrontree->Branch("phoEt",     &hadron_phoEt);
  hadrontree->Branch("phoEta",    &hadron_phoEta);
  hadrontree->Branch("phoPhi",    &hadron_phoPhi);
  hadrontree->Branch("phoSigma",  &hadron_phoSigma);
  hadrontree->Branch("phoChIso",  &hadron_phoChIso);

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
			if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			MET_T1JERUp = raw.pfMET_T1JERUp;
			MET_T1JERDo = raw.pfMET_T1JERDo;
			MET_T1JESUp = raw.pfMET_T1JESUp;
			MET_T1JESDo = raw.pfMET_T1JESDo;
			METPhi_T1JESUp = raw.pfMETPhi_T1JESUp;
			METPhi_T1JESDo = raw.pfMETPhi_T1JESDo;
			METPhi_T1UESUp = raw.pfMETPhi_T1UESUp;
			METPhi_T1UESDo = raw.pfMETPhi_T1UESDo;
			nVtx = raw.nVtx;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			int nGoodMu(0);
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
				if(im->isMedium() && im->getMiniIso() < 0.2 &&  im->getPt() > 15){nGoodMu+=1;}
			}

			if(nGoodMu != 2 || raw.nPho <1)continue;

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			// ******* for jet-photon fake rate***************//
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
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
				}

				if(!itpho->passSignalSelection())continue;
				if(GSFveto && PixelVeto && FSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 20)continue;
				if(!itMu->passHLTSelection())continue;
				if(itMu->passSignalSelection()){
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itMu;
					}
				}
			}

			/*********  ZG tree************/ 
			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.8){
					if(passMETFilter(METFilter)){ 
						
						bool foundZG(false); 
						std::vector<recoMuon>::iterator trailLep=Muon.begin();	
						for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
							if(im == signalLep)continue;
							if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
							double llmass = (im->getP4()+signalLep->getP4()).M();
							double llgmass = (im->getP4()+signalLep->getP4()+signalPho->getP4()).M();
							if( 80 < llmass && llmass < 100){
								foundZG = true;
								trailLep = im;
							}
						}	
					
						if(foundZG){	
							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
							float ThreeBodyMass = (trailLep->getP4()+signalLep->getP4()+signalPho->getCalibP4()).M(); 

							ZphoEt = signalPho->getCalibEt();
							ZphoEta= signalPho->getEta();
							ZphoPhi= signalPho->getPhi();
							ZlepPt = signalLep->getPt();
							ZlepEta= signalLep->getEta();
							ZlepPhi= signalLep->getPhi();
							ZtrailPt = trailLep->getPt();
							ZtrailEta = trailLep->getEta();
							ZtrailPhi = trailLep->getPhi();
							ZsigMT = MT;
							ZsigMET= MET;
							ZsigMETPhi = METPhi;
							ZdPhiLepMET = deltaPhi;
							ZthreeMass = ThreeBodyMass;
							ZdilepMass = (trailLep->getP4()+signalLep->getP4()).M();
							ZnVertex = nVtx;
							ZdRPhoLep= dRlepphoton;

							ZnJet = 0;
							ZHT = 0;
							TLorentzVector JetVec(0,0,0,0);	
							ZISRJetPt = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), trailLep->getEta(),trailLep->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
								ZnJet += 1;
								ZHT += itJet->getPt();
								JetVec = JetVec + itJet->getP4();
							}
							ZISRJetPt = JetVec.Pt();

							
							Ztree->Fill();
						}	

					}//MET Filter
				}//dR Filter
			}//Candidate Filter

		if(hasHadronPho && hasLep){ 
			double DeltaPhoLep = DeltaR(hadronPho->getEta(), hadronPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
			if(DeltaPhoLep > 0.8){
				if(raw.passMETFilter(RunYear, METFilter)){

					bool foundZG(false); 
					std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == signalLep)continue;
						if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
						double llmass = (im->getP4()+signalLep->getP4()).M();
						if( 80 < llmass && llmass < 100){
							foundZG = true;
							trailLep = im;
						}
					}

					if(foundZG){	
						hadron_phoEt = hadronPho->getCalibEt();
						hadron_phoEta = hadronPho->getEta();
						hadron_phoPhi = hadronPho->getPhi();
						hadron_phoSigma = hadronPho->getSigma();
						hadron_phoChIso = hadronPho->getChIso();
			
						hadrontree->Fill();
					}
				}
			}
		}

	}//loop on  events

	p_eventcount->Fill("Total",nTotal);
	p_eventcount->Fill("passHLT",npassHLT);
	p_eventcount->Fill("passPho",npassPho);
	p_eventcount->Fill("passMuon",npassLep);
	p_eventcount->Fill("passdR",npassdR);
	p_eventcount->Fill("passMETFilter",npassMETFilter);
	p_eventcount->Fill("passZ",npassZ);

	outputfile->Write();
	logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_ISR(atoi(argv[1]), argv[2]);
    return 0;
}

