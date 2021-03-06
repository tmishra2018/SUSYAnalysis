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
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TFileCollection.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"

void analysis_mgTrigger(){//main  

	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	char outputname[100] = "/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger_ReMiniAOD.root";
	ofstream logfile;
	logfile.open("/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger_ReMiniAOD.log"); 

	logfile << "analysis_mgTrigger()" << std::endl;

	RunType datatype(SingleMuon2016); 

	TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_C.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_D.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_E.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_F.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_G.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_H.root");


	TFile *outputfile = TFile::Open(outputname,"RECREATE");
	outputfile->cd();

	TTree *Ztree = new TTree("ZTree","ZTree");
	int   Z_run(0);
	float Z_phoEt(0);
	float Z_phoEta(0);
	float Z_phoPhi(0);
	float Z_phoChIso(0);
	float Z_phoNeuIso(0);
	float Z_phoSigma(0);
	float Z_phoR9(0);
	float Z_muPt(0);
	float Z_muMiniIso(0);
	float Z_dR(0); 
	int   Z_phofireL1;
	int   Z_mufireL1;
	int   Z_phofireL1_2;
	int   Z_mufireL1_2;
	int   Z_phofireHLT;
	int   Z_mufireHLT;
	int   Z_mufireHLT2;
	int   Z_phofireHLT2;
	int   Z_nVtx(0);
	float Z_rho(0); 

	Ztree->Branch("runN",      &Z_run);
	Ztree->Branch("phoEt",     &Z_phoEt);
	Ztree->Branch("phoEta",    &Z_phoEta);
	Ztree->Branch("phoPhi",    &Z_phoPhi);
	Ztree->Branch("phoSigma",  &Z_phoSigma);
	Ztree->Branch("phoR9",     &Z_phoR9);
	Ztree->Branch("phoChIso",  &Z_phoChIso);
	Ztree->Branch("phoNeuIso", &Z_phoNeuIso);
	Ztree->Branch("muPt",      &Z_muPt);
	Ztree->Branch("muMiniIso", &Z_muMiniIso);
	Ztree->Branch("dR",        &Z_dR);
	Ztree->Branch("phofireL1", &Z_phofireL1);
	Ztree->Branch("mufireL1",  &Z_mufireL1); 
	Ztree->Branch("phofireL1_2", &Z_phofireL1_2);
	Ztree->Branch("mufireL1_2",  &Z_mufireL1_2); 
	Ztree->Branch("phofireHLT",&Z_phofireHLT);
	Ztree->Branch("mufireHLT", &Z_mufireHLT);
	Ztree->Branch("phofireHLT2",&Z_phofireHLT2);
	Ztree->Branch("mufireHLT2", &Z_mufireHLT2);
	Ztree->Branch("nVtx",      &Z_nVtx);
	Ztree->Branch("rho",       &Z_rho);

	TTree *mgtree = new TTree("mgTree","mgTree");
	float mg_phoEt(0);
	float mg_phoEta(0);
	float mg_phoPhi(0);
	float mg_muPt(0);
	float mg_muMiniIso(0);
	float mg_dR(0); 
	int   mg_phofireL1;
	int   mg_mufireL1;
	int   mg_phofireL1_2;
	int   mg_mufireL1_2;
	int   mg_phofireHLT;
	int   mg_mufireHLT;
	int   mg_phofireHLT2;
	int   mg_mufireHLT2;

	mgtree->Branch("phoEt",     &mg_phoEt);
	mgtree->Branch("phoEta",    &mg_phoEta);
	mgtree->Branch("phoPhi",    &mg_phoPhi);
	mgtree->Branch("muPt",      &mg_muPt);
	mgtree->Branch("muMiniIso", &mg_muMiniIso);
	mgtree->Branch("dR",        &mg_dR);
	mgtree->Branch("phofireL1", &mg_phofireL1);
	mgtree->Branch("mufireL1",  &mg_mufireL1); 
	mgtree->Branch("phofireL1_2",&mg_phofireL1_2);
	mgtree->Branch("mufireL1_2", &mg_mufireL1_2); 
	mgtree->Branch("phofireHLT", &mg_phofireHLT);
	mgtree->Branch("mufireHLT",  &mg_mufireHLT);
	mgtree->Branch("phofireHLT2",&mg_phofireHLT2);
	mgtree->Branch("mufireHLT2", &mg_mufireHLT2);

	TH1F *p_dimuon = new TH1F("p_dimuon","di-muon invmass; #mu#mu mass(GeV);",200,0,200);
	TH1F *p_invmass = new TH1F("p_invmass","mu#mu#gamma invmass; #mu#mu#gamma mass(GeV);",60,60,120);

	const unsigned nEvts = es->GetEntries(); 
	std::cout << "total=" << es->GetEntries() << std::endl;

	rawData raw(es, datatype);
	std::vector<recoPhoton> Photon;
	std::vector<recoMuon>   Muon;
	std::vector<recoEle>   Ele;
	float MET(0);

	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if(ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			Photon.clear();
			Muon.clear();
			Ele.clear();
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			MET = raw.pfMET;

			if(!raw.passHLT())continue;
			Z_nVtx= raw.nVtx;
			Z_rho = raw.rho;
			Z_run = raw.run; 

			std::vector<std::vector<recoMuon>::iterator> tagMuVec;
			std::vector<std::vector<recoPhoton>::iterator> probePhoVec; 
			std::vector<std::pair<std::vector<recoMuon>::iterator, std::vector<recoMuon>::iterator> > DiMuVec;
			std::vector<recoMuon>::iterator  itIsoMu;
			bool hasMu(false);

			tagMuVec.clear();
			DiMuVec.clear();
			probePhoVec.clear();
			for(std::vector<recoMuon>::iterator itLeadMu = Muon.begin(); itLeadMu!= Muon.end(); ++itLeadMu){
				if(!itLeadMu->passHLTSelection())continue;
				if(itLeadMu->passSignalSelection()){
					if(!hasMu){ itIsoMu = itLeadMu; hasMu = true; }
				}
				if(!itLeadMu->isLoose() || itLeadMu->getPt() < 25.0 )continue;
				tagMuVec.push_back(itLeadMu);
			}

			for(int iT(0); iT < tagMuVec.size(); iT++){
				for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu!= Muon.end(); itMu++){

					if( itMu == tagMuVec[iT] )continue;
					bool sameCharge = ((tagMuVec[iT]->isPosi() && itMu->isPosi()) || (!tagMuVec[iT]->isPosi() && !itMu->isPosi()));        
					if( sameCharge || !itMu->isMedium() || itMu->getMiniIso() > 0.2 )continue; 
					float dimuonmass = (tagMuVec[iT]->getP4()+itMu->getP4()).M(); 
					p_dimuon->Fill(dimuonmass);
					if(dimuonmass > 35)DiMuVec.push_back(std::make_pair(tagMuVec[iT],itMu));
				}
			}

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(!itpho->isLoose())continue;
				//if(!itpho->isEB())continue;
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.05)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
				if(GSFveto && PixelVeto && FSRVeto){
					probePhoVec.push_back(itpho);
				}
			}



			for(int iDiMu(0); iDiMu < DiMuVec.size(); iDiMu++){
				for(int ipho(0); ipho < probePhoVec.size(); ipho++){

					std::vector<recoMuon>::iterator tagMu = DiMuVec[iDiMu].first;
					std::vector<recoMuon>::iterator probeMu = DiMuVec[iDiMu].second;
					std::vector<recoPhoton>::iterator probePho = probePhoVec[ipho];
					float InvMass(0);
					float dimuonmass(0);
					float dR1(0),dR2(0);
					dR1 = DeltaR(probePho->getEta(), probePho->getPhi(),tagMu->getEta(), tagMu->getPhi());
					dR2 = DeltaR(probePho->getEta(), probePho->getPhi(),probeMu->getEta(),probeMu->getPhi());
			//		if(dR1 < 0.8 || dR2 < 0.8){
						InvMass = (tagMu->getP4()+probeMu->getP4()+ probePho->getP4()).M();
						dimuonmass= (tagMu->getP4()+probeMu->getP4()).M();
						if( dimuonmass > 35 && dimuonmass  < 80)p_invmass->Fill(InvMass);
			//		}

					if(dimuonmass < 80 && InvMass > 60 && InvMass < 120){
      
						Z_dR = dR1<dR2? dR1: dR2; 
						Z_phoEt = probePho->getEt();
						Z_phoEta=probePho->getEta();
						Z_phoPhi=probePho->getPhi();
						Z_phoChIso=probePho->getChIso();
						Z_phoNeuIso=probePho->getNeuIso();
						Z_phoSigma=probePho->getSigma();
						Z_phoR9 = probePho->getR9();
						Z_muPt=probeMu->getPt();
						Z_muMiniIso=probeMu->getMiniIso();

						if(probePho->fireL1Trg(12))Z_phofireL1=12;
						else if(probePho->fireL1Trg(17))Z_phofireL1=17;
						else Z_phofireL1=0;

						if(probeMu->fireL1Trg(12))Z_mufireL1=12;
						else if(probeMu->fireL1Trg(17))Z_mufireL1=17;
						else Z_mufireL1=0;

						if(probePho->fireL1Trg(29))Z_phofireL1_2=29;
						else Z_phofireL1_2=0;

						if(probeMu->fireL1Trg(29))Z_mufireL1_2=29;
						else Z_mufireL1_2=0;

            if(probePho->fireDoubleTrg(28) || probePho->fireDoubleTrg(29))Z_phofireHLT=1;
						else Z_phofireHLT=0;

						if(probeMu->fireSingleTrg(2) || probeMu->fireSingleTrg(21))Z_mufireHLT=1;
						else Z_mufireHLT=0;

            if(probePho->fireDoubleTrg(30))Z_phofireHLT2=1;
						else Z_phofireHLT2=0;

						if(probeMu->fireSingleTrg(22))Z_mufireHLT2=1;
						else Z_mufireHLT2=0;

            Ztree->Fill();
					}
				}
			}

			if(hasMu){
				for(int ipho(0); ipho < probePhoVec.size(); ipho++){
					std::vector<recoPhoton>::iterator itIsoPho = probePhoVec[ipho];
					float dR = DeltaR(itIsoPho->getEta(), itIsoPho->getPhi(),itIsoMu->getEta(), itIsoMu->getPhi());
					if(dR > 0.8){ 
						mg_dR = dR;
						mg_phoEt =itIsoPho->getEt();
						mg_phoEta=itIsoPho->getEta();
						mg_phoPhi=itIsoPho->getPhi();
						mg_muPt=itIsoMu->getPt();
						mg_muMiniIso=itIsoMu->getMiniIso();

						if(itIsoPho->fireL1Trg(12))mg_phofireL1=12;
						else if(itIsoPho->fireL1Trg(17))mg_phofireL1=17;
						else mg_phofireL1=0;

						if(itIsoMu->fireL1Trg(12))mg_mufireL1=12;
						else if(itIsoMu->fireL1Trg(17))mg_mufireL1=17;
						else mg_mufireL1=0;
						
						if(itIsoPho->fireL1Trg(29))mg_phofireL1_2=29;
						else mg_phofireL1_2=0;

						if(itIsoMu->fireL1Trg(29))mg_mufireL1_2=29;
						else mg_mufireL1_2=0;

						if(itIsoPho->fireDoubleTrg(28) || itIsoPho->fireDoubleTrg(29))mg_phofireHLT=1;
						else mg_phofireHLT=0;

						if(itIsoMu->fireSingleTrg(2) || itIsoMu->fireSingleTrg(21))mg_mufireHLT=1;
						else mg_mufireHLT=0;

						if(itIsoPho->fireDoubleTrg(30))mg_phofireHLT2=1;
						else mg_phofireHLT2=0;

						if(itIsoMu->fireSingleTrg(22))mg_mufireHLT2=1;
						else mg_mufireHLT2=0;

						mgtree->Fill();
					}
				}
			}


    }

outputfile->Write();
}


