#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"
#include "../../../include/analysis_scalefactor.h"

void analysis_rareBkg(){

	std::ifstream configfile("BkgPredConfig.txt");
	int ichannel(1);
	int anatype(0);
	int lowMt(0);
	int highMt(-1);
	int lowMET(0);
	int highMET(-1);
	int lowPt(25);
	int highPt(-1);
	int lepIso(4);
	std::string conftype;
	double confvalue;
	if(configfile.is_open()){
  	for(int i(0); i<9; i++){ 
			configfile >> conftype >> confvalue; 
			if(conftype.find("ichannel")!=std::string::npos)ichannel = confvalue;
			if(conftype.find("anatype")!=std::string::npos)anatype = confvalue;
			if(conftype.find("lowMt")!=std::string::npos)lowMt = confvalue;
			if(conftype.find("highMt")!=std::string::npos)highMt = confvalue;
			if(conftype.find("lowMET")!=std::string::npos)lowMET = confvalue;
			if(conftype.find("highMET")!=std::string::npos)highMET = confvalue;
			if(conftype.find("lowPt")!=std::string::npos)lowPt = confvalue;
			if(conftype.find("highPt")!=std::string::npos)highPt = confvalue;
			if(conftype.find("lepIso")!=std::string::npos)lepIso = confvalue;
	  }
	}
	configfile.close();

	std::ifstream binfile("binConfig.txt");
	float METbin1(200), METbin2(300);
	if(binfile.is_open()){
		for(int i(0); i<2; i++){
			binfile >> conftype >> confvalue;
			if(conftype.find("METbin1")!=std::string::npos)METbin1= confvalue;
			if(conftype.find("METbin2")!=std::string::npos)METbin2= confvalue;
		}
	}

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	esfScaleFactor  objectESF;

  int channelType = ichannel; // eg = 1; mg =2;
	//*********** histo list **********************//
	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_rareBkg";
	else if(channelType==2)outputname << "mg_rareBkg";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	std::ostringstream histname;

	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);

	TH1D *p_reweight_PhoEt = new TH1D("p_reweight_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_reweight_PhoEta = new TH1D("p_reweight_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_reweight_LepPt = new TH1D("p_reweight_LepPt","p_reweight_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_reweight_LepEta = new TH1D("p_reweight_LepEta","p_reweight_LepEta",60,-3,3);
	TH1D *p_reweight_HT = new TH1D("p_reweight_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_reweight_Mt = new TH1D("p_reweight_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *p_reweight_dPhiEleMET = new TH1D("p_reweight_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_reweight_MET = new TH1D("p_reweight_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);

	TH1D *h_rare_norm            = new TH1D("h_rare_norm","eventcount",9,0,9);
	TH1D *h_rare_jesUp           = new TH1D("h_rare_jesUp","eventcount",9,0,9);
	TH1D *h_rare_jesDown         = new TH1D("h_rare_jesDown","eventcount",9,0,9);
	TH1D *h_rare_jerUp           = new TH1D("h_rare_jerUp","eventcount",9,0,9);
	TH1D *h_rare_jerDown         = new TH1D("h_rare_jerDown","eventcount",9,0,9);
	TH1D *h_rare_esfUp           = new TH1D("h_rare_esfUp","eventcount",9,0,9);
	TH1D *h_rare_esfDown         = new TH1D("h_rare_esfDown","eventcount",9,0,9);
	TH1D *h_rare_syserr_jes      = new TH1D("h_rare_syserr_jes","",9,0,9);	
	TH1D *h_rare_syserr_jer      = new TH1D("h_rare_syserr_jer","",9,0,9);	
	TH1D *h_rare_syserr_esf      = new TH1D("h_rare_syserr_esf","",9,0,9);	
	TH1D *h_rare_syserr_scale    = new TH1D("h_rare_syserr_scale","",9,0,9);	
	TH1D *h_rare_syserr_eleshape = new TH1D("h_rare_syserr_eleshape","",9,0,9);	
	TH1D *h_rare_syserr_jetshape = new TH1D("h_rare_syserr_jetshape","",9,0,9);	
	TH1D *h_rare_syserr_xs       = new TH1D("h_rare_syserr_xs","",9,0,9);	
	TH1D *h_rare_syserr_lumi     = new TH1D("h_rare_syserr_lumi","",9,0,9);

	TH1D *jesup_MET = new TH1D("jesup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesup_Mt = new TH1D("jesup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *jesup_HT = new TH1D("jesup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesup_dPhiEleMET = new TH1D("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jesdo_MET = new TH1D("jesdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesdo_Mt = new TH1D("jesdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jesdo_HT = new TH1D("jesdo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesdo_dPhiEleMET = new TH1D("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerup_MET = new TH1D("jerup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerup_Mt = new TH1D("jerup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerup_dPhiEleMET = new TH1D("jerup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerdo_MET = new TH1D("jerdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerdo_Mt = new TH1D("jerdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerdo_dPhiEleMET = new TH1D("jerdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *scaleup_PhoEt = new TH1D("scaleup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *scaleup_PhoEta = new TH1D("scaleup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *scaleup_LepPt = new TH1D("scaleup_LepPt","scaleup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *scaleup_LepEta = new TH1D("scaleup_LepEta","scaleup_LepEta",60,-3,3);
	TH1D *scaleup_MET = new TH1D("scaleup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *scaleup_Mt = new TH1D("scaleup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *scaleup_HT = new TH1D("scaleup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *scaleup_dPhiEleMET = new TH1D("scaleup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *scaledo_PhoEt = new TH1D("scaledo_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *scaledo_PhoEta = new TH1D("scaledo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *scaledo_LepPt = new TH1D("scaledo_LepPt","scaledo_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *scaledo_LepEta = new TH1D("scaledo_LepEta","scaledo_LepEta",60,-3,3);
	TH1D *scaledo_MET = new TH1D("scaledo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *scaledo_Mt = new TH1D("scaledo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *scaledo_HT = new TH1D("scaledo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *scaledo_dPhiEleMET = new TH1D("scaledo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

// ********  MC *************************//
  TChain *mctree;
	if(channelType == 1)mctree = new TChain("egTree","egTree");
  else if(channelType == 2)mctree = new TChain("mgTree","mgTree");
  mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WWG.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WZG.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep1/resTree_VGamma_WW.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep1/resTree_VGamma_WZ.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_TT.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_TTG.root");
	float crosssection(0);
	float ntotalevent(0);
	float PUweight(1);
	int   mcType(0);
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
	float threeMass(0);
  float HT(0);
  float nJet(0);
  int   nBJet(0);
	int   nISRJet(0);
	float bosonPt(0);
  //float invmass(0);  
	float sigMETJESup(0);
	float sigMETJESdo(0);
	float sigMETJERup(0);
	float sigMETJERdo(0);
	float sigMTJESup(0);
	float sigMTJESdo(0);
	float sigMTJERup(0);
	float sigMTJERdo(0);
	float HTJESup(0);
	float HTJESdo(0);
	float dPhiLepMETJESup(0);
	float dPhiLepMETJESdo(0);
	float dPhiLepMETJERup(0);
	float dPhiLepMETJERdo(0);
  std::vector<int> *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
  std::vector<int> *mcStatus=0;

	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	mctree->SetBranchAddress("mcType",    &mcType);
  mctree->SetBranchAddress("phoEt",     &phoEt);
  mctree->SetBranchAddress("phoEta",    &phoEta);
  mctree->SetBranchAddress("phoPhi",    &phoPhi);
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("lepPhi",    &lepPhi);
  mctree->SetBranchAddress("sigMT",     &sigMT);
  mctree->SetBranchAddress("sigMET",    &sigMET);
  mctree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  mctree->SetBranchAddress("nVertex",   &nVertex);
  mctree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	mctree->SetBranchAddress("threeMass", &threeMass);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("nBJet",     &nBJet);
  mctree->SetBranchAddress("nISRJet",   &nISRJet);
  mctree->SetBranchAddress("ISRJetPt",     &bosonPt);
  //mctree->SetBranchAddress("invmass",   &invmass);
	mctree->SetBranchAddress("sigMETJESup",     &sigMETJESup);
	mctree->SetBranchAddress("sigMETJESdo",     &sigMETJESdo);
	mctree->SetBranchAddress("sigMETJERup",     &sigMETJERup);
	mctree->SetBranchAddress("sigMETJERdo",     &sigMETJERdo);
	mctree->SetBranchAddress("sigMTJESup",      &sigMTJESup);
	mctree->SetBranchAddress("sigMTJESdo",      &sigMTJESdo);
	mctree->SetBranchAddress("sigMTJERup",      &sigMTJERup);
	mctree->SetBranchAddress("sigMTJERdo",      &sigMTJERdo);
	mctree->SetBranchAddress("HTJESup",         &HTJESup);
	mctree->SetBranchAddress("HTJESdo",         &HTJESdo);
	mctree->SetBranchAddress("dPhiLepMETJESup", &dPhiLepMETJESup);
	mctree->SetBranchAddress("dPhiLepMETJESdo", &dPhiLepMETJESdo);
	mctree->SetBranchAddress("dPhiLepMETJERup", &dPhiLepMETJERup);
	mctree->SetBranchAddress("dPhiLepMETJERdo", &dPhiLepMETJERdo);
  mctree->SetBranchAddress("mcPID",     &mcPID);
  mctree->SetBranchAddress("mcEta",     &mcEta);
  mctree->SetBranchAddress("mcPhi",     &mcPhi);
  mctree->SetBranchAddress("mcPt",      &mcPt);
  mctree->SetBranchAddress("mcMomPID",  &mcMomPID);
  mctree->SetBranchAddress("mcStatus",  &mcStatus);

	double totalevent(0);
	double totalreweight(0);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(0);
		double scalefactorup(0);
		double scalefactordo(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
			scalefactorup = scalefactor + s_error; 
			scalefactordo = scalefactor - s_error;
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
			scalefactorup = scalefactor + s_error; 
			scalefactordo = scalefactor - s_error;
		}
		float XS_weight = getEvtWeight(RunYear,crosssection,ntotalevent);
		float weight = PUweight*XS_weight*scalefactor;
		float weight_scaleup = PUweight*XS_weight*scalefactorup;
		float weight_scaledo = PUweight*XS_weight*scalefactordo;
		/** cut flow *****/
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) < 1.56 || fabs(lepEta) > 2.5)continue;
		if(fabs(phoEta) > 2.4)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		totalevent += 1;
		double reweightF(1);
		if(bosonPt < 50){reweightF = 1.08893;} 
		else if(bosonPt >= 50 && bosonPt < 80){reweightF  = 1.1675;    }
		else if(bosonPt >= 80 && bosonPt < 100){reweightF  = 0.908314; }
		else if(bosonPt >= 100 && bosonPt < 125){reweightF = 0.848379; }
		else if(bosonPt >= 125 && bosonPt < 150){reweightF = 0.747166; }
		else if(bosonPt >= 150 && bosonPt < 200){reweightF = 0.714054; }
		else if(bosonPt >= 200 && bosonPt < 250){reweightF =	0.729144;}
		else if(bosonPt >= 250 && bosonPt < 300){reweightF =	0.709433;}
		else if(bosonPt >= 300 && bosonPt < 400){reweightF = 0.768222; }
		else if(bosonPt >= 400 && bosonPt < 600){reweightF =  0.452023;}
		else if(bosonPt >= 600){reweightF =  0.252564;}

		totalreweight += reweightF;
		reweightF = reweightF*weight;	

		bool isFake(true);
		double mindRele(0.3), mindRpho(0.3);
		unsigned eleIndex(0), phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			//std::cout << (*mcPID)[iMC] << " " << (*mcMomPID)[iMC] << " " << (*mcStatus)[iMC] << std::endl;
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], lepEta, lepPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE1 = fabs((*mcPt)[iMC] - lepPt)/lepPt;
			double dE2 = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR1 < mindRele){mindRele=dR1; eleIndex=iMC;}
			if(dR2 < mindRpho){mindRpho=dR2; phoIndex=iMC;}
		}
		bool isTrueEle(false);
		if(mindRele < 0.1){
			//std::cout << (*mcPID)[eleIndex] << " " << fabs((*mcMomPID)[eleIndex]) << std::endl;
			if(fabs((*mcPID)[eleIndex]) == 11 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23 || fabs((*mcMomPID)[eleIndex]) == 15))isTrueEle = true;
			else if(fabs((*mcPID)[eleIndex]) == 13 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23 || fabs((*mcMomPID)[eleIndex]) == 15))isTrueEle = true;
			else if(fabs((*mcPID)[eleIndex]) == 15 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23))isTrueEle = true;
		}
		bool isTruePho(false);
		if(mindRpho < 0.1){
			//std::cout << (*mcPID)[phoIndex] << " " << fabs((*mcMomPID)[phoIndex]) << std::endl;
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex])==21 || fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex])== 999 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15 || fabs((*mcMomPID)[phoIndex]) == 24))isTruePho = true;
		}
		bool isFSRPho(false);
		if(mindRpho < 0.1){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15 || fabs((*mcMomPID)[phoIndex])==24 || fabs((*mcMomPID)[phoIndex])==6 ))isFSRPho = true;
			//if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15  ))isFSRPho = true;
		}
		if(isTrueEle && isTruePho)isFake = false;
		//if(isTruePho)isFake = false;
		if(isFake)continue;

		if(mcType >=11 && !isFSRPho)continue; 

		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta,weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta,weight);
		p_MET->Fill(sigMET, weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
		p_nJet->Fill( nJet, weight);

		//if(phoEt > 200 && sigMET > 200 && HT > 400)p_nBJet->Fill(nBJet, weight);
		p_nBJet->Fill(nBJet, weight);

		p_reweight_MET->Fill(sigMET, reweightF);
		p_reweight_PhoEt->Fill(phoEt, reweightF);
		p_reweight_PhoEta->Fill(phoEta,reweightF);
		p_reweight_LepPt->Fill(lepPt, reweightF);
		p_reweight_LepEta->Fill(lepEta,reweightF);
		p_reweight_Mt->Fill(sigMT, reweightF);
		p_reweight_HT->Fill(HT, reweightF);
		p_reweight_dPhiEleMET->Fill(fabs(dPhiLepMET), reweightF);

		jesup_MET->Fill(sigMETJESup, weight);
		jesup_Mt->Fill(sigMTJESup, weight);
		jesup_HT->Fill(HTJESup, weight);
		jesup_dPhiEleMET->Fill(fabs(dPhiLepMETJESup), weight);

		jesdo_MET->Fill(sigMETJESdo, weight);
		jesdo_Mt->Fill(sigMTJESdo, weight);
		jesdo_HT->Fill(HTJESdo, weight);
		jesdo_dPhiEleMET->Fill(fabs(dPhiLepMETJESdo), weight);

		jerup_MET->Fill(sigMETJERup, weight);
		jerup_Mt->Fill(sigMTJERup, weight);
		jerup_dPhiEleMET->Fill(fabs(dPhiLepMETJERup), weight);

		jerdo_MET->Fill(sigMETJERdo, weight);
		jerdo_Mt->Fill(sigMTJERdo, weight);
		jerdo_dPhiEleMET->Fill(fabs(dPhiLepMETJERdo), weight);

		scaleup_PhoEt->Fill(phoEt, weight_scaleup);
		scaleup_PhoEta->Fill(phoEta,weight_scaleup);
		scaleup_LepPt->Fill(lepPt, weight_scaleup);
		scaleup_LepEta->Fill(lepEta,weight_scaleup);
		scaleup_MET->Fill(sigMET, weight_scaleup);
		scaleup_Mt->Fill(sigMT, weight_scaleup);
		scaleup_HT->Fill(HT, weight_scaleup);
		scaleup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaleup);

		scaledo_PhoEt->Fill(phoEt, weight_scaledo);
		scaledo_PhoEta->Fill(phoEta,weight_scaledo);
		scaledo_LepPt->Fill(lepPt, weight_scaledo);
		scaledo_LepEta->Fill(lepEta,weight_scaledo);
		scaledo_MET->Fill(sigMET, weight_scaledo);
		scaledo_Mt->Fill(sigMT, weight_scaledo);
		scaledo_HT->Fill(HT, weight_scaledo);
		scaledo_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaledo);
	}
	p_reweight_MET->Scale(totalevent/totalreweight);
	p_reweight_PhoEt->Scale(totalevent/totalreweight);
	p_reweight_PhoEta->Scale(totalevent/totalreweight);
	p_reweight_LepPt->Scale(totalevent/totalreweight);
	p_reweight_LepEta->Scale(totalevent/totalreweight);
	p_reweight_Mt->Scale(totalevent/totalreweight);
	p_reweight_HT->Scale(totalevent/totalreweight);
	p_reweight_dPhiEleMET->Scale(totalevent/totalreweight);

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((scaledo_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((p_PhoEt->GetBinContent(ibin)*0.3),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
		p_reweight_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((scaledo_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((p_LepPt->GetBinContent(ibin)*0.3),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
		p_reweight_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((scaledo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((jerup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((jerdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((p_MET->GetBinContent(ibin)*0.3),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
		p_reweight_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((scaledo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jerup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jerdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((p_Mt->GetBinContent(ibin)*0.3),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
		p_reweight_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((scaledo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((p_HT->GetBinContent(ibin)*0.3),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
		p_reweight_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow(p_dPhiEleMET->GetBinError(ibin),2);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((scaledo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jesup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jesdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jerup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jerdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	

	p_PhoEt->Sumw2();
	outputfile->Write();
	outputfile->Close();

}


