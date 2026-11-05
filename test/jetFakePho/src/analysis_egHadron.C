//  g++ `root-config --cflags` analysis_egHadron.C  -o analysis_egHadron.exe `root-config --libs`
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
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
#include "../../../src/analysis_rawData.cc"
#include "../../../src/analysis_ele.cc"
#include "../../../src/analysis_photon.cc"
#include "../../../src/analysis_muon.cc"

void analysis_egHadron(int RunYear, const char *Era){//main
  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/logs/plot_hadron_DoubleEG_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_hadron()" << std::endl;
  RunType datatype;
  if(RunYear==2016) datatype = DoubleEG2016;
  if(RunYear==2017) datatype = DoubleEG2017;
  if(RunYear==2018) datatype = DoubleEG2018;

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add(Form("/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_%d%s.root",RunYear,Era));

  TFile *output = TFile::Open(Form("/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_DoubleEG_%d%s.root",RunYear,Era),"RECREATE");
  output->cd();

  TH1F *size_eleproxy = new TH1F("size_eleproxy","",10,0,10);
	TH1F *ele_et = new TH1F("ele_et","",100,0,100);
  int mcType = MCType::NOMC;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;
//************ Signal Tree **********************//
	TTree *signaltree = new TTree("signalTree","signalTree");
	float signal_phoEt(0);
	float signal_phoEta(0);
	float signal_phoPhi(0);
	float signal_phoSigma(0);
	float signal_phoChIso(0);
	float signal_sigMT(0);
	float signal_sigMET(0);
	float signal_sigMETPhi(0);
	float signal_dPhiLepMET(0);
	int   signal_nVertex(0);
	float signal_HT(0);
	float signal_nJet(0);
	std::vector<int> signal_mcPID;
	std::vector<float> signal_mcEta;
	std::vector<float> signal_mcPhi;
	std::vector<float> signal_mcPt;
	std::vector<int> signal_mcMomPID;

	signaltree->Branch("phoEt",     &signal_phoEt);
	signaltree->Branch("phoEta",    &signal_phoEta);
	signaltree->Branch("phoPhi",    &signal_phoPhi);
	signaltree->Branch("phoSigma",  &signal_phoSigma);
	signaltree->Branch("phoChIso",  &signal_phoChIso);
	signaltree->Branch("sigMT",     &signal_sigMT);
	signaltree->Branch("sigMET",    &signal_sigMET);
	signaltree->Branch("sigMETPhi", &signal_sigMETPhi);
	signaltree->Branch("dPhiLepMET",&signal_dPhiLepMET);
	signaltree->Branch("nVertex",   &signal_nVertex);
	signaltree->Branch("HT",        &signal_HT);
	signaltree->Branch("nJet",      &signal_nJet);
	signaltree->Branch("mcPID",     &signal_mcPID);
	signaltree->Branch("mcEta",     &signal_mcEta);
	signaltree->Branch("mcPhi",     &signal_mcPhi);
	signaltree->Branch("mcPt",      &signal_mcPt);
	signaltree->Branch("mcMomPID",  &signal_mcMomPID);
	signaltree->Branch("mcType",    &mcType);

	TTree *eleproxytree = new TTree("eleproxyTree","eleproxyTree");
	float eleproxy_phoEt(0);
	float eleproxy_phoEta(0);
	float eleproxy_phoPhi(0);
	float eleproxy_phoSigma(0);
	float eleproxy_phoChIso(0);
	float eleproxy_sigMT(0);
	float eleproxy_sigMET(0);
	float eleproxy_sigMETPhi(0);
	float eleproxy_dPhiLepMET(0);
	int   eleproxy_nVertex(0);
	float eleproxy_HT(0);
	float eleproxy_nJet(0);

	eleproxytree->Branch("phoEt",     &eleproxy_phoEt);
	eleproxytree->Branch("phoEta",    &eleproxy_phoEta);
	eleproxytree->Branch("phoPhi",    &eleproxy_phoPhi);
	eleproxytree->Branch("phoSigma",  &eleproxy_phoSigma);
	eleproxytree->Branch("phoChIso",  &eleproxy_phoChIso);
	eleproxytree->Branch("sigMT",     &eleproxy_sigMT);
	eleproxytree->Branch("sigMET",    &eleproxy_sigMET);
	eleproxytree->Branch("sigMETPhi", &eleproxy_sigMETPhi);
	eleproxytree->Branch("dPhiLepMET",&eleproxy_dPhiLepMET);
	eleproxytree->Branch("nVertex",   &eleproxy_nVertex);
	eleproxytree->Branch("HT",        &eleproxy_HT);
	eleproxytree->Branch("nJet",      &eleproxy_nJet);

	TTree *hadtree = new TTree("hadTree","hadTree");
	float had_phoEt(0);
	float had_phoEta(0);
	float had_phoPhi(0);
	float had_phoSigma(0);
	float had_phoChIso(0);
	float had_sigMT(0);
	float had_sigMET(0);
	float had_sigMETPhi(0);
	float had_dPhiLepMET(0);
	int   had_nVertex(0);
	float had_HT(0);
	float had_nJet(0);

	hadtree->Branch("phoEt",     &had_phoEt);
	hadtree->Branch("phoEta",    &had_phoEta);
	hadtree->Branch("phoPhi",    &had_phoPhi);
	hadtree->Branch("phoSigma",  &had_phoSigma);
	hadtree->Branch("phoChIso",  &had_phoChIso);
	hadtree->Branch("sigMT",     &had_sigMT);
	hadtree->Branch("sigMET",    &had_sigMET);
	hadtree->Branch("sigMETPhi", &had_sigMETPhi);
	hadtree->Branch("dPhiLepMET",&had_dPhiLepMET);
	hadtree->Branch("nVertex",   &had_nVertex);
	hadtree->Branch("HT",        &had_HT);
	hadtree->Branch("nJet",      &had_nJet);

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

	const unsigned nEvts = es->GetEntries();
	//const unsigned nEvts = 1000000;
	std::cout << "total " << nEvts << std::endl;
	logfile <<   "total " << nEvts << std::endl;

	const bool applyGolden = (datatype != MC);  // golden JSON lumi mask for data
	if(applyGolden && !loadGolden(RunYear)) return;
	for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if(ievt%10000==0)std::cout << " -- Processing event " << ievt << std::endl;
		if(ievt%10000==0)logfile << " -- Processing event " << ievt << std::endl;

		raw.GetData(es, ievt);
		if(applyGolden && !passGolden(raw.run, raw.lumis)) continue;
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

		if(!raw.passHLT())continue;
		if(raw.nGoodVtx < 1)continue;
		if(raw.nPho <1)continue;
		if(MET>70)continue;
  
		std::vector<recoPhoton>::iterator sigPho;
		bool hasSigPho(false);
		std::vector< std::vector<recoPhoton>::iterator> eleproxyVec;
		eleproxyVec.clear(); 
		std::vector< std::vector<recoPhoton>::iterator> hadPhoVec;
		hadPhoVec.clear(); 
		for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
			if(itpho->getR9() < 0.5)continue;
			if(!itpho->passHLTSelection())continue;
			if(!itpho->passBasicSelection())continue; 
      if(itpho->getCalibEt() < 25)continue;
			if(itpho->getChIso() > 20.0)continue;
			bool PixelVeto = itpho->PixelSeed()==0? true: false;
			bool GSFveto(true);
			bool FSRVeto(true);
			for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
			}
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
      if(GSFveto && PixelVeto && FSRVeto){
				if(!hasSigPho){
        	sigPho = itpho;
        	hasSigPho = true;
				}
				hadPhoVec.push_back(itpho);
      }
      if((!PixelVeto || !GSFveto))eleproxyVec.push_back(itpho);
		}
 
		std::vector<std::vector<recoEle>::iterator > signalEleVec;
		signalEleVec.clear();
		for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
			if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
			if(!itEle->passHLTSelection())continue;
			if(itEle->passSignalSelection()){
				signalEleVec.push_back(itEle);
			}
		}
  
		if(hasSigPho && signalEleVec.size() > 0){
			std::vector<recoEle>::iterator signalEle = signalEleVec[0];
			double DeltaPhoLep = DeltaR(sigPho->getEta(), sigPho->getPhi(), signalEle->getEta(), signalEle->getPhi());
			double DoubleMass  = (sigPho->getP4()+signalEle->getP4()).M();
			if(DeltaPhoLep > 0.8 && DoubleMass > 90 && fabs(DoubleMass - 91.188) > 10.0 && sigPho->isEB()){
				float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
				float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
				signal_phoEt = sigPho->getCalibEt();
				signal_phoEta = sigPho->getEta();
				signal_phoPhi = sigPho->getPhi();
				signal_phoSigma = sigPho->getSigma();
				signal_phoChIso = sigPho->getChIso();
				signal_sigMT = MT;
				signal_sigMET = MET;
				signal_sigMETPhi = METPhi;
				signal_dPhiLepMET = deltaPhi;
				signal_nVertex = nVtx;
				signal_HT = 0;
				signal_nJet = jetNumber;
				signal_mcPID.clear();
				signal_mcEta.clear();
				signal_mcPhi.clear();
				signal_mcPt.clear();
				signal_mcMomPID.clear();
				if(datatype == MC){
					for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						if(itMC->getEt() < 1.0)continue;
						float mcdR = DeltaR(sigPho->getEta(), sigPho->getPhi(), itMC->getEta(), itMC->getPhi());
						if(mcdR < 0.3){
							signal_mcPID.push_back(itMC->getPID());
							signal_mcMomPID.push_back(itMC->getMomPID());
							signal_mcEta.push_back(itMC->getEta());      
							signal_mcPhi.push_back(itMC->getPhi());
							signal_mcPt.push_back(itMC->getEt());
						}
					}
				}
				signaltree->Fill();
			}
		}
	
    int nsizeele(0);	
		for(unsigned ip(0); ip < eleproxyVec.size(); ip++){
			for(unsigned ie(0); ie < signalEleVec.size(); ie++){
				std::vector<recoPhoton>::iterator eleproxyPho= eleproxyVec[ip];  
				std::vector<recoEle>::iterator    eleproxyEle= signalEleVec[ie];
				if(!eleproxyPho->isEB())continue;
				double DeltaPhoLep = DeltaR(eleproxyPho->getEta(), eleproxyPho->getPhi(), eleproxyEle->getEta(), eleproxyEle->getPhi());
				double DoubleMass  = (eleproxyPho->getP4()+eleproxyEle->getP4()).M();
				if(DeltaPhoLep > 0.8 && DoubleMass > 101.188){ 
					eleproxy_phoEt = eleproxyPho->getCalibEt();
					eleproxy_phoEta = eleproxyPho->getEta();
					eleproxy_phoPhi = eleproxyPho->getPhi();
					eleproxy_phoSigma = eleproxyPho->getSigma();
					eleproxy_phoChIso = eleproxyPho->getChIso();
					eleproxy_sigMT = 0; 
					eleproxy_sigMET = MET;
					eleproxy_sigMETPhi = METPhi;
					eleproxy_dPhiLepMET = 0;
					eleproxy_nVertex = nVtx;
					eleproxy_HT = 0;
					eleproxy_nJet = jetNumber;
					nsizeele += 1;
					eleproxytree->Fill();
				}
			}
		}
		size_eleproxy->Fill(nsizeele + 0.5);     
		if(nsizeele == 1){
			for(unsigned ie(0); ie < signalEleVec.size(); ie++)ele_et->Fill(signalEleVec[ie]->getPt());
		}
		if(signalEleVec.size() > 0 && !hasSigPho){
			std::vector<recoEle>::iterator hadEle = signalEleVec[0];
			for(unsigned ip(0); ip < hadPhoVec.size(); ip++){
				std::vector<recoPhoton>::iterator hadPho= hadPhoVec[ip];  
				if(!hadPho->isEB())continue;
				double DeltaPhoLep = DeltaR(hadPho->getEta(), hadPho->getPhi(), hadEle->getEta(), hadEle->getPhi());
				double DoubleMass  = (hadPho->getP4()+hadEle->getP4()).M();
				if(DeltaPhoLep > 0.8 && DoubleMass > 101.188){ 
					had_phoEt = hadPho->getCalibEt();
					had_phoEta = hadPho->getEta();
					had_phoPhi = hadPho->getPhi();
					had_phoSigma = hadPho->getSigma();
					had_phoChIso = hadPho->getChIso();
					had_sigMT = 0; 
					had_sigMET = MET;
					had_sigMETPhi = METPhi;
					had_dPhiLepMET = 0;
					had_nVertex = nVtx;
					had_HT = 0;
					had_nJet = jetNumber;
					hadtree->Fill();
				}
			}
		}

	} 
cout<<signaltree->GetEntries()<<endl;
cout<<eleproxytree->GetEntries()<<endl;
cout<<hadtree->GetEntries()<<endl;
logfile << "signalTree: " << signaltree->GetEntries() <<std::endl;
logfile << "eleproxyTree: " << eleproxytree->GetEntries() <<std::endl;
logfile << "hadTree: " << hadtree->GetEntries() <<std::endl;
output->Write();
output->Close();
logfile.close();
}

int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_egHadron(atoi(argv[1]),argv[2]);
    return 0;
}
