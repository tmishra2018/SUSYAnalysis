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

void analysis_VGBkg(){
	
	bool toDeriveScale(false);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	esfScaleFactor  objectESF;
	double WG40Factor(0.86);
	double WG130Factor(0.66);
	double factor_egVGamma(1.29);
	double factorerror_egVGamma(0.25);
	double factor_mgVGamma(1.38);
	double factorerror_mgVGamma(0.25);

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
  	for(int i(0); i<8; i++){ 
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

	if(anatype == 0){
		toDeriveScale = true;
		factor_egVGamma = 1;
	  factorerror_egVGamma = 0;
	  factor_mgVGamma = 1;
	  factorerror_mgVGamma = 0;
	}


  int channelType = ichannel; // eg = 1; mg =2;
	double factorMC(1);
	double factorMCUP = factorMC*(1+0);
	double factorMCDO = factorMC*(1-0);
	if(toDeriveScale){
		factorMC = 1;
		factorMCUP = 1;
		factorMCDO = 1;
	}
	else{
		if(channelType == 1){
			factorMC = factor_egVGamma;
			factorMCUP = factor_egVGamma+factorerror_egVGamma;
			factorMCDO = factor_egVGamma-factorerror_egVGamma;
		}
		else if(channelType == 2){
			factorMC = factor_mgVGamma;
			factorMCUP = factor_mgVGamma + factorerror_mgVGamma;
			factorMCDO = factor_mgVGamma - factorerror_mgVGamma;
		}
	}

	//*********** histo list **********************//
	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_VGBkg";
	else if(channelType==2)outputname << "mg_VGBkg";
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

	TH1D *p_dPhiEleMET_WG = new TH1D("p_dPhiEleMET_WG","dPhiEleMET",32,0,3.2); 
	TH1D *p_dPhiEleMET_ZG = new TH1D("p_dPhiEleMET_ZG","dPhiEleMET",32,0,3.2); 

//	TH1D *p_reweight_PhoEt = new TH1D("p_reweight_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
//	TH1D *p_reweight_PhoEta = new TH1D("p_reweight_PhoEta","#gamma #eta; #eta;",60,-3,3);
//	TH1D *p_reweight_LepPt = new TH1D("p_reweight_LepPt","p_reweight_LepPt",nBkgPtBins,bkgPtBins);
//	TH1D *p_reweight_LepEta = new TH1D("p_reweight_LepEta","p_reweight_LepEta",60,-3,3);
//	TH1D *p_reweight_HT = new TH1D("p_reweight_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
//	TH1D *p_reweight_Mt = new TH1D("p_reweight_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
//	TH1D *p_reweight_dPhiEleMET = new TH1D("p_reweight_dPhiEleMET","dPhiEleMET",32,0,3.2); 
//	TH1D *p_reweight_MET = new TH1D("p_reweight_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);

	TH1D *p_PhoEt_WG = new TH1D("p_PhoEt_WG ","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt_WG = new TH1D("p_LepPt_WG ","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET_WG = new TH1D("p_MET_WG ","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_WG = new TH1D("p_Mt_WG ","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *p_PhoEt_ZG = new TH1D("p_PhoEt_ZG ","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt_ZG = new TH1D("p_LepPt_ZG ","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET_ZG = new TH1D("p_MET_ZG ","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_ZG = new TH1D("p_Mt_ZG ","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 

	TH1D *h_VGamma_norm            = new TH1D("h_VGamma_norm","eventcount",9,0,9);
	TH1D *h_VGamma_jesUp           = new TH1D("h_VGamma_jesUp","eventcount",9,0,9);
	TH1D *h_VGamma_jesDown         = new TH1D("h_VGamma_jesDown","eventcount",9,0,9);
	TH1D *h_VGamma_jerUp           = new TH1D("h_VGamma_jerUp","eventcount",9,0,9);
	TH1D *h_VGamma_jerDown         = new TH1D("h_VGamma_jerDown","eventcount",9,0,9);
	TH1D *h_VGamma_esfUp           = new TH1D("h_VGamma_esfUp","eventcount",9,0,9);
	TH1D *h_VGamma_esfDown         = new TH1D("h_VGamma_esfDown","eventcount",9,0,9);
	TH1D *h_VGamma_normUp          = new TH1D("h_VGamma_normUp","eventcount",9,0,9);
	TH1D *h_VGamma_normDown        = new TH1D("h_VGamma_normDown","eventcount",9,0,9);
	TH1D *h_VGamma_isrUp          = new TH1D("h_VGamma_isrUp","eventcount",9,0,9);
	TH1D *h_VGamma_isrDown        = new TH1D("h_VGamma_isrDown","eventcount",9,0,9);
	TH1D *h_VGamma_syserr_jes      = new TH1D("h_VGamma_syserr_jes","",9,0,9);	
	TH1D *h_VGamma_syserr_jer      = new TH1D("h_VGamma_syserr_jer","",9,0,9);	
	TH1D *h_VGamma_syserr_esf      = new TH1D("h_VGamma_syserr_esf","",9,0,9);	
	TH1D *h_VGamma_syserr_scale    = new TH1D("h_VGamma_syserr_scale","",9,0,9);	
	TH1D *h_VGamma_syserr_eleshape = new TH1D("h_VGamma_syserr_eleshape","",9,0,9);	
	TH1D *h_VGamma_syserr_jetshape = new TH1D("h_VGamma_syserr_jetshape","",9,0,9);	
	TH1D *h_VGamma_syserr_xs       = new TH1D("h_VGamma_syserr_xs","",9,0,9);	
	TH1D *h_VGamma_syserr_lumi     = new TH1D("h_VGamma_syserr_lumi","",9,0,9);
	TH1D *h_VGamma_syserr_isr     = new TH1D("h_VGamma_syserr_isr","",9,0,9);

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

	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *normdo_PhoEt = new TH1D("normdo_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normdo_PhoEta = new TH1D("normdo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normdo_LepPt = new TH1D("normdo_LepPt","normdo_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normdo_LepEta = new TH1D("normdo_LepEta","normdo_LepEta",60,-3,3);
	TH1D *normdo_MET = new TH1D("normdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normdo_Mt = new TH1D("normdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normdo_HT = new TH1D("normdo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normdo_dPhiEleMET = new TH1D("normdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *isrup_PhoEt = new TH1D("isrup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *isrup_PhoEta = new TH1D("isrup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *isrup_LepPt = new TH1D("isrup_LepPt","isrup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *isrup_LepEta = new TH1D("isrup_LepEta","isrup_LepEta",60,-3,3);
	TH1D *isrup_MET = new TH1D("isrup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *isrup_Mt = new TH1D("isrup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *isrup_HT = new TH1D("isrup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *isrup_dPhiEleMET = new TH1D("isrup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *isrdo_PhoEt = new TH1D("isrdo_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *isrdo_PhoEta = new TH1D("isrdo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *isrdo_LepPt = new TH1D("isrdo_LepPt","isrdo_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *isrdo_LepEta = new TH1D("isrdo_LepEta","isrdo_LepEta",60,-3,3);
	TH1D *isrdo_MET = new TH1D("isrdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *isrdo_Mt = new TH1D("isrdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *isrdo_HT = new TH1D("isrdo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *isrdo_dPhiEleMET = new TH1D("isrdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 
// ********  MC *************************//
  TChain *mctree;
 	if(channelType == 1)mctree = new TChain("egTree","egTree");
  else if(channelType == 2)mctree = new TChain("mgTree","mgTree");
	mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG_Pt50.root");
  mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG_Pt35.root");
	mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG_Pt130.root");
	mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_ZG.root");
	mctree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_DY.root");
	//mctree->Add("/uscms_data/d3/mengleis/test/resTree_VGamma_DY10LO.root");
	float crosssection(0);
	float ntotalevent(0);
	int   mcType(0);
	float PUweight(1);
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
  float invmass(0);
	float llmass(0);
	float bosonPt(0); 
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
//  std::vector<int> *mcGMomPID=0;
//
	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("mcType",    &mcType);
	mctree->SetBranchAddress("PUweight",  &PUweight);
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
	mctree->SetBranchAddress("llmass",    &llmass);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  //mctree->SetBranchAddress("invmass",   &invmass);
  mctree->SetBranchAddress("ISRJetPt",     &bosonPt);
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
//  mctree->SetBranchAddress("mcGMomPID", &mcGMomPID);
//
	double totalevent[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double totalreweight[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double totalreweight_up[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double totalreweight_do[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);

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

		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){mindRpho=dR; phoIndex=iMC;}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 23)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 24)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 1)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 2)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 3)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 4)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 5)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 6)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 11)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 13)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 15)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 21)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 999)istruepho=true;
		}
		if(!istruepho)continue;

		totalevent[mcType] += 1;
		double reweightF = 1;
		double reweightF_up = 1;
		double reweightF_do = 1;
		if(bosonPt < 50){reweightF = 1.08893; reweightF_up = reweightF+0.024; reweightF_do = reweightF-0.024; } 
		else if(bosonPt >= 50 && bosonPt < 80){reweightF  = 1.1675; reweightF_up = reweightF+0.020; reweightF_do = reweightF-0.020; }
		else if(bosonPt >= 80 && bosonPt < 100){reweightF  = 0.908314; reweightF_up = reweightF+0.026; reweightF_do = reweightF-0.026; }
		else if(bosonPt >= 100 && bosonPt < 125){reweightF = 0.848379; reweightF_up = reweightF+0.029;reweightF_do = reweightF-0.029; }
		else if(bosonPt >= 125 && bosonPt < 150){reweightF = 0.747166; reweightF_up = reweightF+0.034;reweightF_do = reweightF-0.034; }
		else if(bosonPt >= 150 && bosonPt < 200){reweightF = 0.714054; reweightF_up = reweightF+0.034; reweightF_do = reweightF-0.034;  }
		else if(bosonPt >= 200 && bosonPt < 250){reweightF =	0.729144;reweightF_up = reweightF+0.055; reweightF_do = reweightF-0.055;  }
		else if(bosonPt >= 250 && bosonPt < 300){reweightF =	0.709433;reweightF_up = reweightF+0.12; reweightF_do = reweightF+0.12;  }
		else if(bosonPt >= 300 && bosonPt < 400){reweightF = 0.768222;reweightF_up = reweightF+0.15;  reweightF_do = reweightF+0.15;  }
		else if(bosonPt >= 400 && bosonPt < 600){reweightF =  0.452023; reweightF_up = reweightF+0.22; reweightF_do = reweightF+0.22; }
		else if(bosonPt >= 600){reweightF =  0.252564; reweightF_up = reweightF+0.2; reweightF_do = reweightF-0.2;}
		totalreweight[mcType] += reweightF;
		totalreweight_up[mcType] += reweightF_up;
		totalreweight_do[mcType] += reweightF_do;
	}



	//   start filling ///
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(0);
		double scalefactorup(0);
		double scalefactordo(0);
		if(channelType == 1){
				scalefactor = 1;
				scalefactorup = scalefactor;
				scalefactordo = scalefactor;
	//		scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
	//		double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
	//		double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
	//		double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
	//		double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
	//		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
	//		scalefactorup = scalefactor + s_error; 
	//		scalefactordo = scalefactor - s_error;
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
    if(mcType == 7)crosssection = 18760;
		if((mcType == 7) && llmass > 30)continue;

		float XS_weight = getEvtWeight(RunYear,crosssection,ntotalevent);
	
	
		if(mcType == 2)XS_weight = XS_weight*WG40Factor;
		else if(mcType == 3)XS_weight = XS_weight*WG130Factor;
		
		if(mcType <= 3)XS_weight = XS_weight/WG130Factor;
		if(mcType == 4)XS_weight = XS_weight*(0.14/0.26);
		if(mcType == 5)XS_weight = XS_weight*(4895.0/5670.0);

		if(mcType == 4 && llmass < 30)continue;
		if((mcType == 5 || mcType == 9) && llmass > 30)continue;


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

		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){mindRpho=dR; phoIndex=iMC;}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 23)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 24)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 1)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 2)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 3)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 4)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 5)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 6)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 11)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 13)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 15)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 21)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 999)istruepho=true;
			else{ if(sigMET> 300)std::cout <<mcType << "  " << (*mcPID)[phoIndex] << " " << fabs((*mcMomPID)[phoIndex]) << std::endl;}
		}
		else{ if(sigMET> 300)std::cout << "no match" << std::endl; }
		if(!istruepho)continue;

		double reweightF = 1;
		double reweightF_up = 1;
		double reweightF_do = 1;
		if(bosonPt < 50){reweightF = 1.08893; reweightF_up = reweightF+0.024; reweightF_do = reweightF-0.024; } 
		else if(bosonPt >= 50 && bosonPt < 80){reweightF  = 1.1675; reweightF_up = reweightF+0.020; reweightF_do = reweightF-0.020; }
		else if(bosonPt >= 80 && bosonPt < 100){reweightF  = 0.908314; reweightF_up = reweightF+0.026; reweightF_do = reweightF-0.026; }
		else if(bosonPt >= 100 && bosonPt < 125){reweightF = 0.848379; reweightF_up = reweightF+0.029;reweightF_do = reweightF-0.029; }
		else if(bosonPt >= 125 && bosonPt < 150){reweightF = 0.747166; reweightF_up = reweightF+0.034;reweightF_do = reweightF-0.034; }
		else if(bosonPt >= 150 && bosonPt < 200){reweightF = 0.714054; reweightF_up = reweightF+0.034; reweightF_do = reweightF-0.034;  }
		else if(bosonPt >= 200 && bosonPt < 250){reweightF =	0.729144;reweightF_up = reweightF+0.055; reweightF_do = reweightF-0.055;  }
		else if(bosonPt >= 250 && bosonPt < 300){reweightF =	0.709433;reweightF_up = reweightF+0.12; reweightF_do = reweightF+0.12;  }
		else if(bosonPt >= 300 && bosonPt < 400){reweightF = 0.768222;reweightF_up = reweightF+0.15;  reweightF_do = reweightF+0.15;  }
		else if(bosonPt >= 400 && bosonPt < 600){reweightF =  0.452023; reweightF_up = reweightF+0.22; reweightF_do = reweightF+0.22; }
		else if(bosonPt >= 600){reweightF =  0.252564; reweightF_up = reweightF+0.2; reweightF_do = reweightF-0.2;}
		reweightF = reweightF*weight*totalevent[mcType]/totalreweight[mcType];	
		reweightF_up = reweightF_up*weight*totalevent[mcType]/totalreweight_up[mcType];	
		reweightF_do = reweightF_do*weight*totalevent[mcType]/totalreweight_do[mcType];	
		weight = reweightF;

//		p_reweight_MET->Fill(sigMET, reweightF*factorMC);
//		p_reweight_PhoEt->Fill(phoEt, reweightF*factorMC);
//		p_reweight_PhoEta->Fill(phoEta,reweightF*factorMC);
//		p_reweight_LepPt->Fill(lepPt, reweightF*factorMC);
//		p_reweight_LepEta->Fill(lepEta,reweightF*factorMC);
//		p_reweight_Mt->Fill(sigMT, reweightF*factorMC);
//		p_reweight_HT->Fill(HT, reweightF*factorMC);
//		p_reweight_dPhiEleMET->Fill(fabs(dPhiLepMET), reweightF*factorMC);
	
		p_PhoEt->Fill(phoEt, weight*factorMC);
		p_PhoEta->Fill(phoEta,weight*factorMC);
		p_LepPt->Fill(lepPt, weight*factorMC);
		p_LepEta->Fill(lepEta,weight*factorMC);
		p_MET->Fill(sigMET, weight*factorMC);
		p_Mt->Fill(sigMT, weight*factorMC);
		p_HT->Fill(HT, weight*factorMC);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMC);
		p_nJet->Fill(nJet, weight*factorMC);


		if(mcType <= 3){
			p_PhoEt_WG->Fill(phoEt, weight*factorMC);
			p_LepPt_WG->Fill(lepPt, weight*factorMC);
			p_MET_WG->Fill(sigMET, weight*factorMC);
			p_Mt_WG->Fill(sigMT, weight*factorMC);
			p_dPhiEleMET_WG->Fill(fabs(dPhiLepMET), weight*factorMC);
    }
		else{
			p_PhoEt_ZG->Fill(phoEt, weight*factorMC);
			p_LepPt_ZG->Fill(lepPt, weight*factorMC);
			p_MET_ZG->Fill(sigMET, weight*factorMC);
			p_Mt_ZG->Fill(sigMT, weight*factorMC);
			p_dPhiEleMET_ZG->Fill(fabs(dPhiLepMET), weight*factorMC);
    }

		jesup_MET->Fill(sigMETJESup, weight*factorMC);
		jesup_Mt->Fill(sigMTJESup, weight*factorMC);
		jesup_HT->Fill(HTJESup, weight*factorMC);
		jesup_dPhiEleMET->Fill(fabs(dPhiLepMETJESup), weight*factorMC);

		jesdo_MET->Fill(sigMETJESdo, weight*factorMC);
		jesdo_Mt->Fill(sigMTJESdo, weight*factorMC);
		jesdo_HT->Fill(HTJESdo, weight*factorMC);
		jesdo_dPhiEleMET->Fill(fabs(dPhiLepMETJESdo), weight*factorMC);

		jerup_MET->Fill(sigMETJERup, weight*factorMC);
		jerup_Mt->Fill(sigMTJERup, weight*factorMC);
		jerup_dPhiEleMET->Fill(fabs(dPhiLepMETJERup), weight*factorMC);

		jerdo_MET->Fill(sigMETJERdo, weight*factorMC);
		jerdo_Mt->Fill(sigMTJERdo, weight*factorMC);
		jerdo_dPhiEleMET->Fill(fabs(dPhiLepMETJERdo), weight*factorMC);

		scaleup_PhoEt->Fill(phoEt, weight_scaleup*factorMC);
		scaleup_PhoEta->Fill(phoEta,weight_scaleup*factorMC);
		scaleup_LepPt->Fill(lepPt, weight_scaleup*factorMC);
		scaleup_LepEta->Fill(lepEta,weight_scaleup*factorMC);
		scaleup_MET->Fill(sigMET, weight_scaleup*factorMC);
		scaleup_Mt->Fill(sigMT, weight_scaleup*factorMC);
		scaleup_HT->Fill(HT, weight_scaleup*factorMC);
		scaleup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaleup*factorMC);

		scaledo_PhoEt->Fill(phoEt, weight_scaledo*factorMC);
		scaledo_PhoEta->Fill(phoEta,weight_scaledo*factorMC);
		scaledo_LepPt->Fill(lepPt, weight_scaledo*factorMC);
		scaledo_LepEta->Fill(lepEta,weight_scaledo*factorMC);
		scaledo_MET->Fill(sigMET, weight_scaledo*factorMC);
		scaledo_Mt->Fill(sigMT, weight_scaledo*factorMC);
		scaledo_HT->Fill(HT, weight_scaledo*factorMC);
		scaledo_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaledo*factorMC);

		normup_PhoEt->Fill(phoEt, weight*factorMCUP);
		normup_PhoEta->Fill(phoEta,weight*factorMCUP);
		normup_LepPt->Fill(lepPt, weight*factorMCUP);
		normup_LepEta->Fill(lepEta,weight*factorMCUP);
		normup_MET->Fill(sigMET, weight*factorMCUP);
		normup_Mt->Fill(sigMT, weight*factorMCUP);
		normup_HT->Fill(HT, weight*factorMCUP);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMCUP);

		normdo_PhoEt->Fill(phoEt, weight*factorMCDO);
		normdo_PhoEta->Fill(phoEta,weight*factorMCDO);
		normdo_LepPt->Fill(lepPt, weight*factorMCDO);
		normdo_LepEta->Fill(lepEta,weight*factorMCDO);
		normdo_MET->Fill(sigMET, weight*factorMCDO);
		normdo_Mt->Fill(sigMT, weight*factorMCDO);
		normdo_HT->Fill(HT, weight*factorMCDO);
		normdo_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMCDO);

		isrup_PhoEt->Fill(phoEt, reweightF_up*factorMC);
		isrup_PhoEta->Fill(phoEta,reweightF_up*factorMC);
		isrup_LepPt->Fill(lepPt, reweightF_up*factorMC);
		isrup_LepEta->Fill(lepEta,reweightF_up*factorMC);
		isrup_MET->Fill(sigMET, 0.5*reweightF_up*factorMC);
		isrup_Mt->Fill(sigMT, reweightF_up*factorMC);
		isrup_HT->Fill(HT, reweightF_up*factorMC);
		isrup_dPhiEleMET->Fill(fabs(dPhiLepMET), reweightF_up*factorMC);

		isrdo_PhoEt->Fill(phoEt, reweightF_do*factorMC);
		isrdo_PhoEta->Fill(phoEta,reweightF_do*factorMC);
		isrdo_LepPt->Fill(lepPt, reweightF_do*factorMC);
		isrdo_LepEta->Fill(lepEta,reweightF_do*factorMC);
		isrdo_MET->Fill(sigMET, reweightF_do*factorMC);
		isrdo_Mt->Fill(sigMT, reweightF_do*factorMC);
		isrdo_HT->Fill(HT, reweightF_do*factorMC);
		isrdo_dPhiEleMET->Fill(fabs(dPhiLepMET), reweightF_do*factorMC);
	}
//	p_reweight_MET->Scale(totalevent/totalreweight);
//	p_reweight_PhoEt->Scale(totalevent/totalreweight);
//	p_reweight_PhoEta->Scale(totalevent/totalreweight);
//	p_reweight_LepPt->Scale(totalevent/totalreweight);
//	p_reweight_LepEta->Scale(totalevent/totalreweight);
//	p_reweight_Mt->Scale(totalevent/totalreweight);
//	p_reweight_HT->Scale(totalevent/totalreweight);
//	p_reweight_dPhiEleMET->Scale(totalevent/totalreweight);
//
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((isrup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((isrup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		//syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((isrup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		//syserror += pow((jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		//syserror += pow((jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		//syserror += pow((jerup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		//syserror += pow((jerdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((isrup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jerup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jerdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((isrup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_dPhiEleMET->GetBinError(ibin)* p_dPhiEleMET->GetBinError(ibin);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jesup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jesdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jerup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((jerdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	

	outputfile->Write();
	outputfile->Close();

}


