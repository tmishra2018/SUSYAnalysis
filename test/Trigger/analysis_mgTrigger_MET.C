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
#include "../../include/analysis_jet.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_mcData.h"

#include "../../src/analysis_rawData.cc"
#include "../../src/analysis_ele.cc"
#include "../../src/analysis_muon.cc"
#include "../../src/analysis_photon.cc"


void analysis_mgTrigger_MET(int RunYear, bool preVFP){

	gSystem->Load("../../lib/libAnaClasses.so");
	std::string whichVFP;
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";


	ofstream logfile;
        logfile.open(Form("/eos/uscms/store/user/tmishra/Trigger/plot_MuonTrigger_MET_%d%s.log",RunYear,whichVFP.c_str()));

        logfile << "analysis_mgTrigger()" << std::endl;

        RunType datatype(MC);

        TChain* es = new TChain("ggNtuplizer/EventTree");
        es->Add(Form("/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/MET/MET_%d%s.root",RunYear,whichVFP.c_str()));

        TFile *outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/Trigger/plot_MuonTrigger_MET_%d%s.root",RunYear,whichVFP.c_str()),"RECREATE");
        outputfile->cd();

	TTree *mgtree = new TTree("mgTree","mgTree");
	float mg_phoEt(0);
	float mg_phoEta(0);
	float mg_phoPhi(0);
	float mg_muPt(0);
	float mg_muMiniIso(0);
	float mg_dR(0); 
	std::vector<int>   mg_phofireL1;
	std::vector<int>   mg_mufireL1;
	int   mg_phofireHLT;
	int   mg_mufireHLT;
	int   mg_phofireHLT2;
	int   mg_mufireHLT2;
	int   mg_passHLT;
	int   mg_passHLT2;

	mgtree->Branch("phoEt",     &mg_phoEt);
	mgtree->Branch("phoEta",    &mg_phoEta);
	mgtree->Branch("phoPhi",    &mg_phoPhi);
	mgtree->Branch("muPt",      &mg_muPt);
	mgtree->Branch("muMiniIso", &mg_muMiniIso);
	mgtree->Branch("dR",        &mg_dR);
	mgtree->Branch("phofireL1", &mg_phofireL1);
	mgtree->Branch("mufireL1",  &mg_mufireL1); 
	mgtree->Branch("phofireHLT", &mg_phofireHLT);
	mgtree->Branch("mufireHLT",  &mg_mufireHLT);
	mgtree->Branch("phofireHLT2",&mg_phofireHLT2);
	mgtree->Branch("mufireHLT2", &mg_mufireHLT2);
  	mgtree->Branch("passHLT",    &mg_passHLT);
  	mgtree->Branch("passHLT2",   &mg_passHLT2);

	const unsigned nEvts = es->GetEntries(); 
	std::cout << "total=" << es->GetEntries() << std::endl;

	rawData raw(es, datatype);
	std::vector<recoPhoton> Photon;
	std::vector<recoMuon>   Muon;
	std::vector<recoEle>   Ele;
	float MET(0);

	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if(ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			Photon.clear();
			Muon.clear();
			Ele.clear();
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			MET = raw.pfMET;

			if(!raw.passHLT())continue;
			bool hasMu(false);
			std::vector<recoMuon>::iterator tagMu;
			std::vector<std::vector<recoPhoton>::iterator> probePhoVec;
			probePhoVec.clear(); 

			for(std::vector<recoMuon>::iterator itLeadMu = Muon.begin(); itLeadMu!= Muon.end(); ++itLeadMu){
				if(!itLeadMu->isMedium() || itLeadMu->getMiniIso() > 0.2 || fabs(itLeadMu->getEta()) > 2.4 )continue;
				if(itLeadMu->getD0() > 0.05 || itLeadMu->getDz() > 0.1)continue; 
				if(!hasMu){
					tagMu = itLeadMu;
					hasMu = true;
				}
			}

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(!itpho->isLoose())continue;
				if(!itpho->isEB())continue;
				if(itpho->getR9() < 0.5)continue;
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.05)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
				//if(GSFveto && PixelVeto && FSRVeto){
				if (FSRVeto){
					probePhoVec.push_back(itpho);
				}
			}

			mg_passHLT = 0;
			mg_passHLT2 = 0;
		  mg_phofireL1.clear();
			mg_mufireL1.clear();	
			if(hasMu){
				for(std::size_t ipho(0); ipho < probePhoVec.size(); ipho++){
					std::vector<recoPhoton>::iterator probePho = probePhoVec[ipho];
					float dR = DeltaR(probePho->getEta(), probePho->getPhi(),tagMu->getEta(), tagMu->getPhi());
							mg_dR = dR;
							mg_phoEt =probePho->getCalibEt();
							mg_phoEta=probePho->getEta();
							mg_phoPhi=probePho->getPhi();
							mg_muPt=tagMu->getPt();
							mg_muMiniIso=tagMu->getMiniIso();

    							if(((raw.HLTEleMuX >> 8) &1) !=0)mg_passHLT = 8; //  HLT_Mu17_Photon30_CaloIdL_L1ISO_v
							// or HLT_Mu17_Photon30_IsoCaloId_v for 2017 and 2018
							if(RunYear == 2016 and ((raw.HLTEleMuX >> 51) &1) !=0) mg_passHLT2 = 51; // HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_v
							if(RunYear == 2017 and ((raw.HLTEleMuX >> 57) &1) !=0) mg_passHLT2 = 57; // HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v
							if(RunYear == 2018 and ((raw.HLTEleMuX >> 57) &1) !=0) mg_passHLT2 = 57; // HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v
						
							for(int iT(0); iT < 30; iT++){
								if(probePho->fireL1Trg(iT))mg_phofireL1.push_back(iT);
								if(tagMu->fireL1Trg(iT))mg_mufireL1.push_back(iT);
							}	

							if(probePho->fireDoubleTrg(28) || probePho->fireDoubleTrg(29))mg_phofireHLT=28;
							else mg_phofireHLT=0;

							if(tagMu->fireSingleTrg(2) || tagMu->fireSingleTrg(21))mg_mufireHLT=21;
							else mg_mufireHLT=0;

							if(probePho->fireDoubleTrg(30))mg_phofireHLT2=30;
							else mg_phofireHLT2=0;

							if(tagMu->fireSingleTrg(22))mg_mufireHLT2=22;
							else mg_mufireHLT2=0;

							mgtree->Fill();
				}
			}

    }

outputfile->Write();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    bool preVFP = (atoi(argv[2]) == 1);
    analysis_mgTrigger_MET(atoi(argv[1]), preVFP);
    return 0;
}
