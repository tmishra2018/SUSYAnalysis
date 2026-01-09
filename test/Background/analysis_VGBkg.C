#include<string>
#include "../../include/analysis_commoncode.h"



void analysis_VGBkg(){
	
	SetRunConfig();
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
	bool noISRweight = true;
	
	esfScaleFactor  objectESF;
	bool toDeriveScale(false);
	if(anatype == 0)toDeriveScale = true;

  	int channelType = ichannel; // eg = 1; mg =2;
	double factorMC(1);
	double factorMCUP = factorMC*(1+0);
	
	 // make this ON when you have to derive the QCD and VGamma scales
	if(toDeriveScale){
		factorMC = 1;
		factorMCUP = 1;
	}
	else{
		if(channelType == 1){
			if(RunYear==2016 and preVFP == 1){
				factorMC = factor_egVGamma_2016preVFP;
				factorMCUP = factor_egVGamma_2016preVFP+factorerror_egVGamma_2016preVFP;}
			if(RunYear==2016 and preVFP == 0){
				factorMC = factor_egVGamma_2016postVFP;
				factorMCUP = factor_egVGamma_2016postVFP+factorerror_egVGamma_2016postVFP;}
			if(RunYear==2017){
				factorMC = factor_egVGamma_2017;
				factorMCUP = factor_egVGamma_2017+factorerror_egVGamma_2017;}
			if(RunYear==2018){
				factorMC = factor_egVGamma_2018;
				factorMCUP = factor_egVGamma_2018+factorerror_egVGamma_2018;}
		}
		else if(channelType == 2){
			if(RunYear==2016 and preVFP == 1){
				factorMC = factor_mgVGamma_2016preVFP;
				factorMCUP = factor_mgVGamma_2016preVFP+factorerror_mgVGamma_2016preVFP;}
			if(RunYear==2016 and preVFP == 0){
				factorMC = factor_mgVGamma_2016postVFP;
				factorMCUP = factor_mgVGamma_2016postVFP+factorerror_mgVGamma_2016postVFP;}
			if(RunYear==2017){
				factorMC = factor_mgVGamma_2017;
				factorMCUP = factor_mgVGamma_2017+factorerror_mgVGamma_2017;}
			if(RunYear==2018){
				factorMC = factor_mgVGamma_2018;
				factorMCUP = factor_mgVGamma_2018+factorerror_mgVGamma_2018;}
		}
	}

	//*********** histo list ********************* 
	std::ostringstream outputname;
	outputname << "/eos/uscms/store/user/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_VGBkg";
	else if(channelType==2)outputname << "mg_VGBkg";
	// MET and lepton pT range
	if(anatype ==0 or anatype ==1)	outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname <<"_" << RunYear<<whichVFP<<".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	std::ostringstream histname;

	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","Lepton p_{T}",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2);
	TH1D *p_dPhiEleMET_WG = new TH1D("p_dPhiEleMET_WG","dPhiEleMET",32,0,3.2); 
	TH1D *p_dPhiEleMET_ZG = new TH1D("p_dPhiEleMET_ZG","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);

	TH1D *p_LepPt_TT = new TH1D("p_LepPt_TT","p_LepPt",nBkgPtBins,bkgPtBins);
        TH1D *p_nJet_TT = new TH1D("p_nJet_TT","p_nJet",10,0,10);
        TH1D *p_nBJet_TT = new TH1D("p_nBJet_TT","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET_TT = new TH1D("p_dPhiEleMET_TT","dPhiEleMET",32,0,3.2);

	TH1D *jesup_MET = new TH1D("jesup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesup_Mt = new TH1D("jesup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *jesup_HT = new TH1D("jesup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesup_dPhiEleMET = new TH1D("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	TH1D *jesup_MET_TT = new TH1D("jesup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesup_Mt_TT = new TH1D("jesup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *jesup_HT_TT = new TH1D("jesup_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesup_dPhiEleMET_TT = new TH1D("jesup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *jesdo_MET = new TH1D("jesdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesdo_Mt = new TH1D("jesdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jesdo_HT = new TH1D("jesdo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesdo_dPhiEleMET = new TH1D("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	TH1D *jesdo_MET_TT = new TH1D("jesdo_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jesdo_Mt_TT = new TH1D("jesdo_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jesdo_HT_TT = new TH1D("jesdo_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *jesdo_dPhiEleMET_TT = new TH1D("jesdo_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *jerup_MET = new TH1D("jerup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerup_Mt = new TH1D("jerup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerup_dPhiEleMET = new TH1D("jerup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerup_MET_TT = new TH1D("jerup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerup_Mt_TT = new TH1D("jerup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerup_dPhiEleMET_TT = new TH1D("jerup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *jerdo_MET = new TH1D("jerdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerdo_Mt = new TH1D("jerdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerdo_dPhiEleMET = new TH1D("jerdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerdo_MET_TT = new TH1D("jerdo_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *jerdo_Mt_TT = new TH1D("jerdo_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *jerdo_dPhiEleMET_TT = new TH1D("jerdo_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *scaleup_PhoEt = new TH1D("scaleup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *scaleup_PhoEta = new TH1D("scaleup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *scaleup_LepPt = new TH1D("scaleup_LepPt","scaleup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *scaleup_LepEta = new TH1D("scaleup_LepEta","scaleup_LepEta",60,-3,3);
	TH1D *scaleup_MET = new TH1D("scaleup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *scaleup_Mt = new TH1D("scaleup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *scaleup_HT = new TH1D("scaleup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *scaleup_dPhiEleMET = new TH1D("scaleup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *scaleup_PhoEt_TT = new TH1D("scaleup_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *scaleup_LepPt_TT = new TH1D("scaleup_LepPt_TT","scaleup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *scaleup_MET_TT = new TH1D("scaleup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *scaleup_Mt_TT = new TH1D("scaleup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *scaleup_HT_TT = new TH1D("scaleup_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *scaleup_dPhiEleMET_TT = new TH1D("scaleup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	TH1D *normup_PhoEt_TT  = new TH1D("normup_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_LepPt_TT  = new TH1D("normup_LepPt_TT","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_MET_TT  = new TH1D("normup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt_TT  = new TH1D("normup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT_TT  = new TH1D("normup_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET_TT = new TH1D("normup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *isrup_PhoEt = new TH1D("isrup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *isrup_PhoEta = new TH1D("isrup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *isrup_LepPt = new TH1D("isrup_LepPt","isrup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *isrup_LepEta = new TH1D("isrup_LepEta","isrup_LepEta",60,-3,3);
	TH1D *isrup_MET = new TH1D("isrup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *isrup_Mt = new TH1D("isrup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *isrup_HT = new TH1D("isrup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *isrup_dPhiEleMET = new TH1D("isrup_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	TH1D *isrup_PhoEt_TT  = new TH1D("isrup_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *isrup_LepPt_TT  = new TH1D("isrup_LepPt_TT","isrup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *isrup_MET_TT  = new TH1D("isrup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *isrup_Mt_TT  = new TH1D("isrup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *isrup_HT_TT  = new TH1D("isrup_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *isrup_dPhiEleMET_TT = new TH1D("isrup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	TH1D *toy_dPhiEleMET[500];
	for(unsigned ih(0); ih < 500; ih++){
		histname.str("");
		histname << "toy_VGdPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}
	TH1D *toy_dPhiEleMET_TT[500];
	for(unsigned ih(0); ih < 500; ih++){
		histname.str("");
		histname << "toy_VGdPhiEleMET_TT_" << ih;
		toy_dPhiEleMET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}
	// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
	// ZG bkg is directly from simulation, mctree
  	TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
 	
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
        mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
        mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
        mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_%d%s.root",RunYear,whichVFP.c_str()));
        mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
	

	float crosssection(0);
	float ntotalevent(0);
	float ISRWeight(0);
	float pdfWeight(0);
	std::vector<float> *pdfSystWeight = 0;
	std::vector<float> *ScaleSystWeight = 0; 
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
  	float HT(0);
  	float nJet(0);
	int   nBJet(0);
	float llmass(0);
	float ISRJetPt(0); 
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
	
	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("ISRWeight", &ISRWeight);
	mctree->SetBranchAddress("pdfWeight", &pdfWeight);
	mctree->SetBranchAddress("pdfSystWeight", &pdfSystWeight);
	mctree->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight);
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
	mctree->SetBranchAddress("llmass",    &llmass);
  	mctree->SetBranchAddress("HT",        &HT);
  	mctree->SetBranchAddress("nJet",      &nJet);
  	mctree->SetBranchAddress("nBJet",     &nBJet);
  	mctree->SetBranchAddress("ISRJetPt",        &ISRJetPt);
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

	//   start filling ///
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(1);
		double scalefactorup(1);
		if(channelType == 1){
			// from include/analysis_scalefactor.h, src/analysis_scalefactor.cc
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
			scalefactorup = scalefactor + s_error; 
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      			double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
			scalefactorup = scalefactor + s_error; 
		}
		// ZG sample has llmass > 30, we use DY for llmass < 30
		if(mcType == 4 && llmass < 30)continue;
		if(mcType == 5 && llmass > 30)continue;

		float XS_weight = 1;
		
                if(channelType == 1){
                        if(RunYear == 2016 and preVFP == 1)             XS_weight = lumi_2016preVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        XS_weight = lumi_2016postVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2017)                        XS_weight = lumi_2017_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2018)                        XS_weight = lumi_2018_DoubleEG*1000*crosssection/ntotalevent;}

                else if(channelType == 2){
                        if(RunYear == 2016 and preVFP == 1)             XS_weight = lumi_2016preVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        XS_weight = lumi_2016postVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2017)                        XS_weight = lumi_2017_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2018)                        XS_weight = lumi_2018_MuonEG*1000*crosssection/ntotalevent;}


		double reweightPt=1.0;
		double NormPt=1.0;
		
		if(RunYear==2016 && preVFP==1){
			if(ichannel == 1){     
				NormPt = 0.902349;
                 		if(phoEt < 50)reweightPt = 1.28212;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 1.03049;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.9531;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 0.953816;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 1.03967;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 1.07921;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 1.10223;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 1.42878;}
			else if(ichannel == 2){
				NormPt = 1.01353;
                 		if(phoEt < 50)reweightPt = 1.09103;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 0.890823;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.689855;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 0.617876;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 0.575871;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 0.634576;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 0.632047;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 0.873833;}
		}
                else if(RunYear==2016 && preVFP==0){
			if(ichannel == 1){
				NormPt = 1.11914;
                 		if(phoEt < 50)reweightPt = 1.06132;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 0.748704;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.870131;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 1.18874;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 1.18161;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 1.08044;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 1.30831;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 0.481584;}

			else if(ichannel == 2){
				NormPt = 1.04551;
                 		if(phoEt < 50)reweightPt = 1.08342;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 0.820422;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.723338;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 0.56966;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 0.460661;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 0.440407;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 0.583253;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 0.432163;}
		}
                else if(RunYear==2017){
			if(ichannel == 1){
				NormPt = 0.942539;
                 		if(phoEt < 50)reweightPt = 1.14592;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 0.950863;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 1.12263;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 1.3817;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 1.39352;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 1.4926;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 1.27273;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 1.04809;}
			else if(ichannel == 2){
				NormPt = 0.955486;
                 		if(phoEt < 50)reweightPt = 1.05087;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 1.07524;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.92892;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 0.825291;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 0.899771;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 0.837221;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 0.732376;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 0.685758;}
		}
                else if(RunYear==2018){
			if(ichannel == 1){
				NormPt = 0.907516;
                 		if(phoEt < 50)reweightPt = 1.21675;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 1.03955;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 1.00025;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 1.08549;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 1.14779;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 1.31126;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 1.08117;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 1.3491;}
			else if(ichannel == 2){
				NormPt = 0.989044;
                 		if(phoEt < 50)reweightPt = 1.09491;
                 		else if(phoEt >= 50 && phoEt < 100)reweightPt  = 0.939154;
                 		else if(phoEt >= 100 && phoEt < 150)reweightPt = 0.732848;
                 		else if(phoEt >= 150 && phoEt < 200)reweightPt = 0.691871;
                 		else if(phoEt >= 200 && phoEt < 250)reweightPt = 0.676227;
                 		else if(phoEt >= 250 && phoEt < 300)reweightPt = 0.523843;
                 		else if(phoEt >= 300 && phoEt < 500)reweightPt = 0.569649;
                 		else if(phoEt >= 500 && phoEt < 800)reweightPt = 0.347228;}
		}
		// check all weights used
		
		float weight = 1.0;
		float weight_scaleup = 1.0;
		float weight_normup = 1.0;
		float weight_noisr = 1.0;
		
                weight = PUweight*XS_weight*factorMC*scalefactor;
                weight_scaleup = PUweight*XS_weight*factorMC*scalefactorup;
                weight_normup = PUweight*XS_weight*factorMCUP*scalefactor;
                weight_noisr = PUweight*XS_weight*factorMC*scalefactor;

	        if (!toDeriveScale){
			weight = weight*reweightPt*NormPt;
			weight_scaleup = weight_scaleup*reweightPt*NormPt;
			weight_normup = weight_normup*reweightPt*NormPt;
			weight_noisr = weight_noisr*reweightPt*NormPt;
		}

		double weight_toy[500];
		for(unsigned ii(0); ii < 500; ii++){
			weight_toy[ii] = 1;
			//if(pdfSystWeight->size() < ii)weight_toy[ii] = weight;
			//else weight_toy[ii] = weight*(*pdfSystWeight)[ii]/pdfWeight;
		}
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		// MET, MT and lepton pT range
		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){mindRpho=dR; phoIndex=iMC;}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )istruepho=true;
		}
		if(!istruepho)continue;
		
		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta,weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta,weight);
		p_MET->Fill(sigMET, weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
		// WG samples
		if(mcType <= 3)p_dPhiEleMET_WG->Fill(fabs(dPhiLepMET), weight);
		// ZG samples
		else if( mcType == 4 || mcType == 5)p_dPhiEleMET_ZG->Fill(fabs(dPhiLepMET), weight);
		p_nJet->Fill(nJet, weight);
		p_nBJet->Fill(nBJet, weight);


		// JES and JER are stored in ntuple
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

		normup_PhoEt->Fill(phoEt,  weight_normup);
		normup_PhoEta->Fill(phoEta,weight_normup);
		normup_LepPt->Fill(lepPt,  weight_normup);
		normup_LepEta->Fill(lepEta,weight_normup);
		normup_MET->Fill(sigMET,   weight_normup);
		normup_Mt->Fill(sigMT,     weight_normup);
		normup_HT->Fill(HT,        weight_normup);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_normup);

		isrup_PhoEt->Fill(phoEt, weight_noisr);
		isrup_PhoEta->Fill(phoEta,weight_noisr);
		isrup_LepPt->Fill(lepPt, weight_noisr);
		isrup_LepEta->Fill(lepEta,weight_noisr);
		isrup_MET->Fill(sigMET, weight_noisr);
		isrup_Mt->Fill(sigMT, weight_noisr);
		isrup_HT->Fill(HT, weight_noisr);
		isrup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_noisr);
		
		
		
		// when >= 2 bJets, tt events
		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  weight);
			p_MET_TT->Fill(sigMET,  weight);
			p_Mt_TT->Fill(sigMT,  weight);
			p_HT_TT->Fill(HT,  weight);
			p_LepPt_TT->Fill(lepPt,  weight);
                        p_nJet_TT->Fill(nJet,  weight);
                        p_nBJet_TT->Fill(nBJet, weight);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), weight);

			jesup_MET_TT->Fill(sigMETJESup, weight);
			jesup_Mt_TT->Fill(sigMTJESup, weight);
			jesup_HT_TT->Fill(HTJESup, weight);
			jesup_dPhiEleMET_TT->Fill(fabs(dPhiLepMETJESup), weight);

			jesdo_MET_TT->Fill(sigMETJESdo, weight);
			jesdo_Mt_TT->Fill(sigMTJESdo, weight);
			jesdo_HT_TT->Fill(HTJESdo, weight);
			jesdo_dPhiEleMET_TT->Fill(fabs(dPhiLepMETJESdo), weight);

			jerup_MET_TT->Fill(sigMETJERup, weight);
			jerup_Mt_TT->Fill(sigMTJERup, weight);
			jerup_dPhiEleMET_TT->Fill(fabs(dPhiLepMETJERup), weight);

			jerdo_MET_TT->Fill(sigMETJERdo, weight);
			jerdo_Mt_TT->Fill(sigMTJERdo, weight);
			jerdo_dPhiEleMET_TT->Fill(fabs(dPhiLepMETJERdo), weight);

			scaleup_PhoEt_TT->Fill(phoEt, weight_scaleup);
			scaleup_LepPt_TT->Fill(lepPt, weight_scaleup);
			scaleup_MET_TT->Fill(sigMET, weight_scaleup);
			scaleup_Mt_TT->Fill(sigMT, weight_scaleup);
			scaleup_HT_TT->Fill(HT, weight_scaleup);
			scaleup_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), weight_scaleup);

			normup_PhoEt_TT->Fill(phoEt,  weight_normup);
			normup_LepPt_TT->Fill(lepPt,  weight_normup);
			normup_MET_TT->Fill(sigMET,   weight_normup);
			normup_Mt_TT->Fill(sigMT,     weight_normup);
			normup_HT_TT->Fill(HT,        weight_normup);
			normup_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), weight_normup);

			isrup_PhoEt_TT->Fill(phoEt, weight_noisr);
			isrup_LepPt_TT->Fill(lepPt, weight_noisr);
			isrup_MET_TT->Fill(sigMET, weight_noisr);
			isrup_Mt_TT->Fill(sigMT, weight_noisr);
			isrup_HT_TT->Fill(HT, weight_noisr);
			isrup_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), weight_noisr);
		}
			for(unsigned ii(0);  ii < 500; ii++)toy_dPhiEleMET[ii]->Fill(fabs(dPhiLepMET), weight_toy[ii]); 
			for(unsigned ii(0);  ii < 500; ii++)toy_dPhiEleMET_TT[ii]->Fill(fabs(dPhiLepMET), weight_toy[ii]); 

	}


	// Stat and Syst error
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);                    // stat
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);  // ID and trigger ESF
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);   // normalisation scale
		syserror += pow((isrup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);    // ISR corrections
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt_TT->GetBinError(ibin)* p_PhoEt_TT->GetBinError(ibin);                    // stat
		syserror += pow((scaleup_PhoEt_TT->GetBinContent(ibin)-p_PhoEt_TT->GetBinContent(ibin)),2);  // ID and trigger ESF
		syserror += pow((normup_PhoEt_TT->GetBinContent(ibin)-p_PhoEt_TT->GetBinContent(ibin)),2);   // normalisation scale
		syserror += pow((isrup_PhoEt_TT->GetBinContent(ibin)-p_PhoEt_TT->GetBinContent(ibin)),2);    // ISR corrections
		p_PhoEt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((isrup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt_TT->GetBinError(ibin)* p_LepPt_TT->GetBinError(ibin);
		syserror += pow((scaleup_LepPt_TT->GetBinContent(ibin)-p_LepPt_TT->GetBinContent(ibin)),2);
		syserror += pow((normup_LepPt_TT->GetBinContent(ibin)-p_LepPt_TT->GetBinContent(ibin)),2);
		syserror += pow((isrup_LepPt_TT->GetBinContent(ibin)-p_LepPt_TT->GetBinContent(ibin)),2);
		p_LepPt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((isrup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin))); // JEC
		double jererror = max( fabs(jerup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jerdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin))); // JER
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET_TT->GetBinError(ibin)* p_MET_TT->GetBinError(ibin);
		syserror += pow((scaleup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)),2);
		syserror += pow((normup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)),2);
		syserror += pow((isrup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)), fabs(jesdo_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin))); // JEC
		double jererror = max( fabs(jerup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)), fabs(jerdo_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin))); // JER
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_MET_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((isrup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jerdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt_TT->GetBinError(ibin)* p_Mt_TT->GetBinError(ibin);
		syserror += pow((scaleup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)),2);
		syserror += pow((normup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)),2);
		syserror += pow((isrup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)), fabs(jesdo_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)), fabs(jerdo_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_Mt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((isrup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)), fabs(jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT_TT->GetBinError(ibin)* p_HT_TT->GetBinError(ibin);
		syserror += pow((normup_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)),2);
		syserror += pow((isrup_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)),2);
		syserror += pow((scaleup_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)), fabs(jesdo_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		p_HT_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_dPhiEleMET->GetBinError(ibin)* p_dPhiEleMET->GetBinError(ibin);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((normup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((isrup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jesdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jerdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_dPhiEleMET_TT->GetBinError(ibin)* p_dPhiEleMET_TT->GetBinError(ibin);
		syserror += pow((scaleup_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)),2);
		syserror += pow((normup_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)),2);
		syserror += pow((isrup_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)), fabs(jesdo_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)), fabs(jerdo_dPhiEleMET_TT->GetBinContent(ibin)-p_dPhiEleMET_TT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_dPhiEleMET_TT->SetBinError(ibin,sqrt(syserror));
	}	
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	p_dPhiEleMET_TT->Write();
	outputfile->Write();
	outputfile->Close();

}
