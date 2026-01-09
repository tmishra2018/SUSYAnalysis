
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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

double Normalization = 0.98;

bool isHardLepton(int momID){
		if((fabs(momID) >= 0 && fabs(momID) <= 6) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
		else return false;
}
		
bool isHardPhoton(int momID){
		if((fabs(momID) >= 0 && fabs(momID) <= 16) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
		else return false;
}

bool nisrMatch(float jetEta, float jetPhi, std::vector<mcData>& genParticles){
    bool matched=false;
    for(std::vector<mcData>::iterator itMC = genParticles.begin(); itMC!= genParticles.end(); itMC++){
      			if(matched)break;
			bool isPromptFinal(false);
      			int momid = abs(itMC->getMomPID());
      			if(abs(itMC->getPID())<=5 && (momid==6 || momid==23 || momid==24 || momid==25))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 11 && (momid<=6 || momid==23 || momid==24 || momid== 22 || momid== 15))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 13 && (momid<=6 || momid==23 || momid==24 || momid== 22 || momid== 15))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 15 && (momid<=6 || momid==23 || momid==24 || momid== 22))isPromptFinal = true;
			else if(abs(itMC->getPID()) == 22 && (momid == 11 || momid == 13 || momid == 15 || momid==23 || momid==24))isPromptFinal = true;
			if(!isPromptFinal)continue;
			if(DeltaR(jetEta, jetPhi, itMC->getEta(), itMC->getPhi()) < 0.3){
				//std::cout << "match " << itMC->getPID() << " " << itMC->getMomPID() << " status:" << itMC->getStatus() << std::endl; 
      				matched = true;
			}
    } // Loop over MC particles
    return matched;
} 
		
bool passFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 8; im++){
		if(im != 5){ // no ee Bad Sc filter
    	if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}
		}
	}

  return passfilter;
}


void analysis_VGamma(int RunYear, bool preVFP, const char *Sample){//main 
  
  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_%s_%d%s.log",Sample,RunYear,whichVFP.c_str()));
  logfile << "analysis_VGamma()" << std::endl;

  RunType datatype(MC);
  TChain* es = new TChain("ggNtuplizer/EventTree");
  char* inputfile = new char[300];

  if (strstr(Sample, "DYJetsToLL") != NULL or strstr(Sample, "TTJets") != NULL or strstr(Sample, "WJetsToLNu"))
        sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  else
        sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  es->Add(inputfile);

  logfile << "VG";
  logfile << "muon MiniIso" << std::endl;

  const unsigned nEvts = es->GetEntries();
  logfile << "Total events : " << nEvts << std::endl;
  
  int npassEGselection(0), npassdR(0), npassZ(0), npassMETFilter(0), npassPho(0), npassEle(0), npassHLTPho(0);
  
  TFile* outputfile;
  outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_%s_%d%s.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
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
  else  if(strstr(inputfile, "TTGJets") != NULL){
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
  else if(strstr(inputfile, "ZZ") != NULL){
                std::cout << "ZZ sample !" << std::endl;
                mcType = MCType::ZZ;
  }


  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;

  TH1F *h_dReg = new TH1F("h_dReg", "DeltaR(electron, photon);#DeltaR(e, #gamma);Events", 50, 0, 5.0);
  TH1F *h_dRmg = new TH1F("h_dRmg", "DeltaR(muon, photon);#DeltaR(#mu, #gamma);Events", 50, 0, 5.0);

  float eg_Z_mass(0);
  float eg_Z_lep1_pt(0); 
  float eg_Z_lep2_pt(0);
  
  float mg_Z_mass(0);
  float mg_Z_lep1_pt(0); 
  float mg_Z_lep2_pt(0);
  
  bool eg_hasZCandidate = false;
  bool mg_hasZCandidate = false;
  
  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();
  float PUweight(1);
  float llmass(0);
  int   nBJet(0);
  float ISRWeight(1);
  float pdfWeight(0);
  std::vector<float> pdfSystWeight;
  std::vector<float> ScaleSystWeight; 
//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_lepPt(0);
  float eg_lepEta(0);
  float eg_lepPhi(0);
  float trailPt(0);
  float trailEta(0);
  float trailPhi(0);
  float eg_sigMT(0);
  float eg_sigMET(0);
  float eg_sigMETPhi(0);
  float eg_dPhiLepMET(0);
  int   eg_nVertex(0);
  float eg_dRPhoLep(0);
  float eg_HT(0);
  float eg_nJet(0);
  int   eg_nISRJet(0);
  float eg_ISRJetPt(0);
  float eg_invmass(0);
  float eg_bosonPt(0);
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
  std::vector<int> eg_mcPID;
  std::vector<float> eg_mcEta;
  std::vector<float> eg_mcPhi;
  std::vector<float> eg_mcPt;
  std::vector<int> eg_mcMomPID;
  std::vector<int> eg_mcGMomPID;
  std::vector<int> eg_mcStatus;

  egtree->Branch("crosssection",&crosssection);
  egtree->Branch("ntotalevent", &ntotalevent);
  egtree->Branch("ISRWeight", &ISRWeight);
  egtree->Branch("pdfWeight", &pdfWeight);
  egtree->Branch("pdfSystWeight", &pdfSystWeight);
  egtree->Branch("ScaleSystWeight", &ScaleSystWeight);
  egtree->Branch("phoEt",     &eg_phoEt);
  egtree->Branch("phoEta",    &eg_phoEta);
  egtree->Branch("phoPhi",    &eg_phoPhi);
  egtree->Branch("lepPt",     &eg_lepPt);
  egtree->Branch("lepEta",    &eg_lepEta);
  egtree->Branch("lepPhi",    &eg_lepPhi);
  egtree->Branch("trailPt",   &trailPt);
  egtree->Branch("trailEta",  &trailEta);
  egtree->Branch("trailPhi",  &trailPhi);
  egtree->Branch("PUweight",  &PUweight);
  egtree->Branch("sigMT",     &eg_sigMT);
  egtree->Branch("sigMET",    &eg_sigMET);
  egtree->Branch("sigMETPhi", &eg_sigMETPhi);
  egtree->Branch("dPhiLepMET",&eg_dPhiLepMET);
  egtree->Branch("nVertex",   &eg_nVertex);
  egtree->Branch("dRPhoLep",  &eg_dRPhoLep);
  egtree->Branch("HT",        &eg_HT);
  egtree->Branch("nJet",      &eg_nJet);
  egtree->Branch("nBJet",     &nBJet);
  egtree->Branch("nISRJet",   &eg_nISRJet);
  egtree->Branch("ISRJetPt",  &eg_ISRJetPt);
  egtree->Branch("invmass",   &eg_invmass);
  egtree->Branch("bosonPt",   &eg_bosonPt);
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
  egtree->Branch("mcType",   &mcType);
  egtree->Branch("mcPID",    &eg_mcPID);
  egtree->Branch("mcEta",    &eg_mcEta);
  egtree->Branch("mcPhi",    &eg_mcPhi);
  egtree->Branch("mcPt",     &eg_mcPt);
  egtree->Branch("mcMomPID", &eg_mcMomPID);
  egtree->Branch("mcGMomPID", &eg_mcGMomPID);
  egtree->Branch("mcStatus", &eg_mcStatus);
  egtree->Branch("llmass",   &llmass);

  egtree->Branch("Z_mass", &eg_Z_mass);
  egtree->Branch("Z_lep1_pt", &eg_Z_lep1_pt);
  egtree->Branch("Z_lep2_pt", &eg_Z_lep2_pt);
  egtree->Branch("hasZCandidate", &eg_hasZCandidate);

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
  float mg_nISRJet(0);
  float mg_ISRJetPt(0);
  float mg_bosonPt(0);
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
  std::vector<int> mg_mcPID;
  std::vector<float> mg_mcEta;
  std::vector<float> mg_mcPhi;
  std::vector<float> mg_mcPt;
  std::vector<int> mg_mcMomPID;
  std::vector<int> mg_mcGMomPID;
  std::vector<int> mg_mcStatus;
 
  mgtree->Branch("crosssection",&crosssection);
  mgtree->Branch("ntotalevent", &ntotalevent);
  mgtree->Branch("ISRWeight", &ISRWeight);
  mgtree->Branch("pdfWeight", &pdfWeight);
  mgtree->Branch("pdfSystWeight", &pdfSystWeight);
  mgtree->Branch("ScaleSystWeight", &ScaleSystWeight);
  mgtree->Branch("phoEt",     &mg_phoEt);
  mgtree->Branch("phoEta",    &mg_phoEta);
  mgtree->Branch("phoPhi",    &mg_phoPhi);
  mgtree->Branch("lepPt",     &mg_lepPt);
  mgtree->Branch("lepEta",    &mg_lepEta);
  mgtree->Branch("lepPhi",    &mg_lepPhi);
  mgtree->Branch("PUweight",      &PUweight);
  mgtree->Branch("sigMT",     &mg_sigMT);
  mgtree->Branch("sigMET",    &mg_sigMET);
  mgtree->Branch("sigMETPhi", &mg_sigMETPhi);
  mgtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
  mgtree->Branch("threeMass", &mg_threeMass);
  mgtree->Branch("nVertex",   &mg_nVertex);
  mgtree->Branch("dRPhoLep",  &mg_dRPhoLep);
  mgtree->Branch("HT",        &mg_HT);
  mgtree->Branch("nJet",      &mg_nJet);
  mgtree->Branch("nBJet",     &nBJet);
  mgtree->Branch("nISRJet",   &mg_nISRJet);
  mgtree->Branch("ISRJetPt",  &mg_ISRJetPt);
  mgtree->Branch("bosonPt",   &mg_bosonPt);
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
  mgtree->Branch("mcType",   &mcType);
  mgtree->Branch("mcPID",    &mg_mcPID);
  mgtree->Branch("mcEta",    &mg_mcEta);
  mgtree->Branch("mcPhi",    &mg_mcPhi);
  mgtree->Branch("mcPt",     &mg_mcPt);
  mgtree->Branch("mcMomPID", &mg_mcMomPID);
  mgtree->Branch("mcGMomPID", &mg_mcGMomPID);
  mgtree->Branch("mcStatus", &mg_mcStatus);
  mgtree->Branch("llmass",   &llmass);
  
  mgtree->Branch("Z_mass", &mg_Z_mass);
  mgtree->Branch("Z_lep1_pt", &mg_Z_lep1_pt);
  mgtree->Branch("Z_lep2_pt", &mg_Z_lep2_pt);
  mgtree->Branch("hasZCandidate", &mg_hasZCandidate);

  //*********** histo list **********************//
  TH1F *p_eventcount = new TH1F("p_eventcount","p_eventcount",7,0,7);
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  /*********************************************/

  std::vector<recoEle>::iterator ZLep1_e = Ele.begin();
  std::vector<recoEle>::iterator ZLep2_e = Ele.begin();
  std::vector<recoMuon>::iterator ZLep1_mu = Muon.begin();
  std::vector<recoMuon>::iterator ZLep2_mu = Muon.begin();
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
  float ISRJetPt(0);
  logfile << "RunType: " << datatype << std::endl;

  std::cout << "Total events : " << nEvts << std::endl;

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
				
				
        if(preVFP==true and RunYear==2016) PUweight = getPUESF16preVFP(nVtx);
        if(preVFP==false and RunYear==2016) PUweight = getPUESF16(nVtx);
        if(RunYear==2017) PUweight = getPUESF17(nVtx);
        if(RunYear==2018) PUweight = getPUESF18(nVtx);
	
        if(raw.nPho <1)continue;

	int Nmedpho(0);
	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
	        if(!itpho->isMedium())continue;
	        if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
	        if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
	        if(itpho->getSigma() < 0.005)continue;
	        if(itpho->PixelSeed() != 0)continue;
	
	        Nmedpho+=1;
	}
	if(Nmedpho >= 2)continue;

	nBJet = 0;
	for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
		if(itJet->getPt() < 20)continue;
		if(itJet->isBJet())nBJet+=1;
	}
	//pdfWeight = raw.pdfWeight;
	pdfWeight = 0; 
	pdfSystWeight.clear();
	ScaleSystWeight.clear(); 
	//	for(unsigned i(0); i < raw.pdfSystWeight->size(); i++){
	//		pdfSystWeight.push_back( (*raw.pdfSystWeight)[i] );
	//	}
	//	for(unsigned i(0); i < raw.genScaleSystWeights->size(); i++){
	//		ScaleSystWeight.push_back( (*raw.genScaleSystWeights)[i]);
	//	}
	//process MC first  ****************************/
	llmass = -1;
	double WGPt = -1;
	if(datatype == MC){
		 std::vector<mcData>::iterator posMcMu; 
		 bool hasposMcMu(false);
		 std::vector<mcData>::iterator negMcMu;
		 bool hasnegMcMu(false);
		 std::vector<mcData>::iterator posMcTau; 
		 bool hasposMcTau(false);
		 std::vector<mcData>::iterator negMcTau;
		 bool hasnegMcTau(false);
		 std::vector<mcData>::iterator posMcEle; 
		 bool hasposMcEle(false);
		 std::vector<mcData>::iterator negMcEle;
		 bool hasnegMcEle(false);
		 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					switch(itMC->getPID()){
						case 13: if(isHardLepton(itMC->getMomPID())){posMcMu = itMC; hasposMcMu=true;};break;
						case -13:if(isHardLepton(itMC->getMomPID())){negMcMu = itMC; hasnegMcMu=true;};break;
						case 15: if(isHardLepton(itMC->getMomPID())){posMcTau = itMC; hasposMcTau=true;};break;
						case -15:if(isHardLepton(itMC->getMomPID())){negMcTau = itMC; hasnegMcTau=true;};break;
						case 11: if(isHardLepton(itMC->getMomPID())){posMcEle = itMC; hasposMcEle=true;};break;
						case -11:if(isHardLepton(itMC->getMomPID())){negMcEle = itMC; hasnegMcEle=true;};break;
						case 23: WGPt = itMC->getPt(); break;
						case 24: WGPt = itMC->getPt(); break;
						case -24: WGPt = itMC->getPt(); break;
						default: break;
					}
	         }
		 if(hasposMcMu && hasnegMcMu)llmass = (posMcMu->getP4() + negMcMu->getP4()).M();
		 else if(hasposMcTau && hasnegMcTau)llmass = (posMcTau->getP4() + negMcTau->getP4()).M();
		 else if(hasposMcEle && hasnegMcEle)llmass = (posMcEle->getP4() + negMcEle->getP4()).M();
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
		// No leptonic activities surrounding photons
		for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
			if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
			if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
		}
	
		for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
			if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
		                
		if(GSFveto && PixelVeto && FSRVeto){
	//		if((itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6))){ // Year = 2016 Year = 2017
		 		if(!hasegPho){
					hasegPho=true;
		 			egsignalPho = itpho;
				}
	//		}
			if(itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30)){ // Year = 2016
		 		if(!hasmgPho){
					hasmgPho=true;
		 			mgsignalPho = itpho;
				}
			}
		}
        }
     
	
	if(hasegPho){ 
		if(egsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasegPho=false;
		if(egsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasegPho=false;
	}
	if(hasmgPho){ 
		if(mgsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasmgPho=false;
		if(mgsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasmgPho=false;
	}
 
        bool hasEle(false);
        std::vector<recoEle>::iterator signalEle = Ele.begin();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
        	if(hasEle)break;
		if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
		//if((itEle->isEB() && (itEle->getD0() > 0.05 || itEle->getDz() > 0.10)) || (itEle->isEE() && ( itEle->getD0() > 0.10 || itEle->getDz() > 0.20)) )continue;
		if(itEle->fireTrgs(21) || itEle->fireTrgs(22)){   // Year = 2016
			if(itEle->passSignalSelection()){
				hasEle=true; 
				signalEle = itEle;
			}
		}
        }

        bool hasMu(false);
        std::vector<recoMuon>::iterator signalMu = Muon.begin();
        
	if(hasmgPho){
		for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
			if(hasMu)break;
			if(itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22)){  // Year = 2016
				if(itMu->passSignalSelection()){
					hasMu=true; 
					signalMu = itMu;
				}
			}
		}
        }

	bool hasTrail(false);
	std::vector<recoMuon>::iterator trailLep = Muon.begin();
	for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
		if(itMu->getPt() < 25)continue;
		if(itMu->passSignalSelection()){
			if(!hasTrail){
				hasTrail = true;
				trailLep = itMu;
			}
		}	
	}

	bool hasDoubleEG(false);
	
	ISRJetPt = 0;
	TLorentzVector JetVec(0,0,0,0);	
	for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
		if(!itJet->passSignalSelection())continue;
		if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))JetVec = JetVec + itJet->getP4();
	}
	ISRJetPt = JetVec.Pt();
	double reweightF(1);	
				
	if(RunYear==2016 && preVFP==1){
	        if(ISRJetPt < 50)reweightF = 0.967889;
             	else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF  = 1.11724;
                else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.979645;
                else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.833069;
                else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.881191;
                else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.895074;
                else if(ISRJetPt >= 300)reweightF = 0.816445;
                Normalization = 0.988597;
        }
        else if(RunYear==2016 && preVFP==0){
              	if(ISRJetPt < 50)reweightF = 1.0078;
                else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF  = 1.18628;
                else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.940295;
                else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.926074;
                else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.856512;
                else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.752113;
                else if(ISRJetPt >= 300)reweightF = 0.671831;
                 Normalization = 0.955131;
        }

       	ISRWeight = reweightF*Normalization;

std::vector<std::pair<recoEle, recoEle>> ZeeCandidates;
std::vector<std::pair<recoMuon, recoMuon>> ZmmCandidates;

for (std::vector<recoEle>::iterator itEle1 = Ele.begin(); itEle1 != Ele.end(); ++itEle1) {
    for (std::vector<recoEle>::iterator itEle2 = itEle1 + 1; itEle2 != Ele.end(); ++itEle2) { 
        if (itEle1->passSignalSelection() && itEle2->passSignalSelection()) {
            if ((itEle1->isPosi() && itEle2->isPosi()) || (!itEle1->isPosi() && !itEle2->isPosi())) continue;
            float dimass = (itEle1->getP4() + itEle2->getP4()).M();
            if (fabs(dimass - 91.1876) < 10) {
                ZeeCandidates.emplace_back(*itEle1, *itEle2);
		goto EndZSearch;
            }
        }
    }
}

for (std::vector<recoMuon>::iterator itMu1 = Muon.begin(); itMu1 != Muon.end(); ++itMu1) {
    for (std::vector<recoMuon>::iterator itMu2 = itMu1 + 1; itMu2 != Muon.end(); ++itMu2) { 
        if (itMu1->passSignalSelection() && itMu2->passSignalSelection()) {
            if ((itMu1->isPosi() && itMu2->isPosi()) || (!itMu1->isPosi() && !itMu2->isPosi())) continue;
            float dimuonmass = (itMu1->getP4() + itMu2->getP4()).M();
            if (fabs(dimuonmass - 91.1876) < 10) {
                ZmmCandidates.emplace_back(*itMu1, *itMu2);
		goto EndZSearch;
            }
        }
    }
}
EndZSearch:

bool hasZCandidate = (!ZeeCandidates.empty() || !ZmmCandidates.empty());
int totalZCandidates = ZeeCandidates.size() + ZmmCandidates.size();

if (totalZCandidates > 1) 
    std::cout << "More than one Z candidate found in this event: " << totalZCandidates << std::endl;



	///*************************   eg filters *****************************
	if(hasegPho) npassPho+=1;
	if(hasEle)   npassEle+=1;
	if(((raw.HLTPho >> 14)&1)==1) npassHLTPho+=1;


        if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){ // Year = 2016 Year = 2017
	  	npassEGselection+=1;
          	double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
		h_dReg->Fill(dReg);
          	if(dReg>0.8){
	  		npassdR+=1;
            		if(((egsignalPho->getP4()+signalEle->getP4()).M() - 91.188) > 10.0){
	  		npassZ+=1;
                	if(passFilter(METFilter)){
				npassMETFilter+=1;
				hasDoubleEG = true;
				float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
				float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
        	    		eg_phoEt = egsignalPho->getEt();
				eg_phoEta= egsignalPho->getEta();
                		eg_phoPhi= egsignalPho->getPhi();
				eg_lepPt = signalEle->getPt();
 				eg_lepEta= signalEle->getEta();
                		eg_lepPhi= signalEle->getPhi();
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
 				eg_sigMT = MT;
				eg_sigMET= MET;
                		eg_sigMETPhi = METPhi;
				eg_dPhiLepMET = deltaPhi; 
				eg_nVertex = nVtx; 
				eg_dRPhoLep= dReg;
                		eg_invmass = (egsignalPho->getP4()+signalEle->getP4()).M();
				eg_bosonPt = WGPt;
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
							
				eg_nJet = 0;
				eg_nISRJet = 0;
				eg_HT = 0;
				eg_HTJESup = 0;
				eg_HTJESdo = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
									
					if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))eg_nISRJet += 1;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(), signalEle->getPhi()) <= 0.4)continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalPho->getEta(),egsignalPho->getPhi()) <= 0.4)continue;	
					eg_nJet += 1;
					eg_HT += itJet->getPt();
					eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
					eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
				}
				eg_ISRJetPt = ISRJetPt;	
 
				if(datatype == MC){
					eg_mcPID.clear();
					eg_mcEta.clear();
					eg_mcPhi.clear();
					eg_mcPt.clear();
					eg_mcMomPID.clear();
					eg_mcGMomPID.clear();
					eg_mcStatus.clear();
					for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						if(itMC->getEt() < 5.0)continue;
						eg_mcPID.push_back(itMC->getPID());
						eg_mcMomPID.push_back(itMC->getMomPID());
						eg_mcEta.push_back(itMC->getEta());      
						eg_mcPhi.push_back(itMC->getPhi());
						eg_mcPt.push_back(itMC->getEt());
						eg_mcGMomPID.push_back(itMC->getGMomPID());
						eg_mcStatus.push_back(itMC->getStatus());
					}
				}

eg_hasZCandidate = false;
eg_Z_mass = -1;
eg_Z_lep1_pt = -1;
eg_Z_lep2_pt = -1;
if (hasZCandidate) {
    for (auto& pair : ZeeCandidates) {
        if (&pair.first != &(*signalEle) && &pair.second != &(*signalEle)) {
            eg_hasZCandidate = true;
            eg_Z_mass = (pair.first.getP4() + pair.second.getP4()).M();
            eg_Z_lep1_pt = pair.first.getPt();
            eg_Z_lep2_pt = pair.second.getPt();
            break;
        }
    }
    if (!eg_hasZCandidate) {
        for (auto& pair : ZmmCandidates) {
            eg_hasZCandidate = true;
            eg_Z_mass = (pair.first.getP4() + pair.second.getP4()).M();
            eg_Z_lep1_pt = pair.first.getPt();
            eg_Z_lep2_pt = pair.second.getPt();
            break;
        }
    }
}

                		egtree->Fill();
 
              }//MET Filter
            }// Z mass Filter
          }//dR filter
        }// ele + pho candidate
       
//**********************  mg filter **************************************//         
 
       if(hasmgPho && hasMu && (((raw.HLTEleMuX >> 8)&1)!=0 || ((raw.HLTEleMuX >> 51)&1)!=0) && !hasDoubleEG){  // Year = 2016
					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
					h_dRmg->Fill(dRmg);
					if(dRmg>0.8){
						if(passFilter(METFilter)){ 

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
						mg_bosonPt = WGPt;
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
						mg_nISRJet = 0;
						mg_HT = 0;
						mg_HTJESup = 0;
						mg_HTJESdo = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
						//	bool isLep(false);
						//	for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
						//		if(!itEle->isMedium() || itEle->getMiniIso() > 0.1 || itEle->getPt() < 25)continue;
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itEle->getEta(), itEle->getPhi()) <= 0.4 && fabs(itJet->getPt() - itEle->getPt())/itEle->getPt() < 1)isLep=true;
						//	}
						//	for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
						//		if(!itMu->isMedium() || itMu->getMiniIso() > 0.2 || itMu->getPt() < 25)continue;
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4 && fabs(itJet->getPt() - itMu->getPt())/itMu->getPt() < 1)isLep=true;
						//	}
						//	if(isLep)continue;

							if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))mg_nISRJet += 1;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= 0.4)continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), mgsignalPho->getEta(),mgsignalPho->getPhi()) <= 0.4)continue;	
							mg_nJet += 1;
							mg_HT += itJet->getPt();
							mg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
							mg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
						}	
						mg_ISRJetPt = ISRJetPt;

						if(datatype == MC){
						 mg_mcPID.clear();
						 mg_mcEta.clear();
						 mg_mcPhi.clear();
						 mg_mcPt.clear();
						 mg_mcMomPID.clear();
						 mg_mcGMomPID.clear();
						 mg_mcStatus.clear();
						 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	if(itMC->getEt() < 1.0)continue;
							float mcdRmu = DeltaR(signalMu->getEta(), signalMu->getPhi(), itMC->getEta(), itMC->getPhi());
						 	float mcdR = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						 	if(mcdR < 0.3 || mcdRmu < 0.3){
						 		mg_mcPID.push_back(itMC->getPID());
						 		mg_mcMomPID.push_back(itMC->getMomPID());
						 		mg_mcGMomPID.push_back(itMC->getGMomPID());
						 		mg_mcEta.push_back(itMC->getEta());
						 		mg_mcPhi.push_back(itMC->getPhi());
						 		mg_mcPt.push_back(itMC->getEt());
								mg_mcStatus.push_back(itMC->getStatus());
						 	}
						  }
						}
mg_hasZCandidate = false;
mg_Z_mass = -1;
mg_Z_lep1_pt = -1;
mg_Z_lep2_pt = -1;
if (hasZCandidate) {
    for (auto& pair : ZmmCandidates) {
        if (&pair.first != &(*signalMu) && &pair.second != &(*signalMu)) {
            mg_hasZCandidate = true;
            mg_Z_mass = (pair.first.getP4() + pair.second.getP4()).M();
            mg_Z_lep1_pt = pair.first.getPt();
            mg_Z_lep2_pt = pair.second.getPt();
            break;
        }
    }
    if (!mg_hasZCandidate) {
        for (auto& pair : ZeeCandidates) {
            mg_hasZCandidate = true;
            mg_Z_mass = (pair.first.getP4() + pair.second.getP4()).M();
            mg_Z_lep1_pt = pair.first.getPt();
            mg_Z_lep2_pt = pair.second.getPt();
            break;
        }
    }
}

				 		mgtree->Fill();
           }//MET Filter
         }//dR Filter
       }//Candidate Filter
 
	}//loop on  events

  p_eventcount->GetXaxis()->SetBinLabel(1,"npassPho");
  p_eventcount->GetXaxis()->SetBinLabel(2,"npassEle");
  p_eventcount->GetXaxis()->SetBinLabel(3,"npassHLTPho");
  p_eventcount->GetXaxis()->SetBinLabel(4,"npassEGselection");
  p_eventcount->GetXaxis()->SetBinLabel(5,"npassdR");
  p_eventcount->GetXaxis()->SetBinLabel(6,"npassZ");
  p_eventcount->GetXaxis()->SetBinLabel(7,"npassMETFilter");
  
  p_eventcount->Fill(0.5, npassPho);
  p_eventcount->Fill(1.5, npassEle);
  p_eventcount->Fill(2.5, npassHLTPho);
  p_eventcount->Fill(3.5, npassEGselection);
  p_eventcount->Fill(4.5, npassdR);
  p_eventcount->Fill(5.5, npassZ);
  p_eventcount->Fill(6.5, npassMETFilter);
  cout<<egtree->GetEntries()<<"\t"<<npassPho<<"\t"<<npassEle<<"\t"<<npassHLTPho<<"\t"<<npassEGselection<<"\t"<<npassdR<<"\t"<<npassZ<<"\t"<<npassMETFilter<<endl;
  outputfile->Write();
  h_dReg->Write();
  h_dRmg->Write();
  outputfile->Close();
  logfile.close();

}

int main(int argc, char** argv)
{

    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    bool preVFP = (atoi(argv[2]) == 1);
    analysis_VGamma(atoi(argv[1]), preVFP, argv[3]);
    return 0;

}
