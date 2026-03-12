#include<string>
#include "../../include/analysis_commoncode.h"

void analysis_rareBkg(){

	SetRunConfig();
	setTDRStyle();

  gSystem->Load("../../lib/libAnaClasses.so");
	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	esfScaleFactor  objectESF;

  int channelType = ichannel; // eg = 1; mg =2;
	//*********** histo list **********************//
	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_rareBkg";
	else if(channelType==2)outputname << "mg_rareBkg";
	if(anatype ==0 or anatype ==1)	outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname <<"_" << RunYear<<whichVFP <<".root";

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

	const int nHTSlices = 3;
	double HTsliceMin[nHTSlices] = {0.,   100., 400.};
	double HTsliceMax[nHTSlices] = {100., 400., 1e6};
	TH2D *p_MET_vs_HTslice = new TH2D("p_MET_vs_HTslice", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
	p_MET_vs_HTslice->GetYaxis()->SetBinLabel(1,"HT: 0-100");
	p_MET_vs_HTslice->GetYaxis()->SetBinLabel(2,"HT: 100-400");
	p_MET_vs_HTslice->GetYaxis()->SetBinLabel(3,"HT: >400");
	        
	TH2D *scaleup_MET_HT = new TH2D("scaleup_MET_HT", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        TH2D *jesup_MET_HT = new TH2D("jesup_MET_HT", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        TH2D *jerup_MET_HT = new TH2D("jerup_MET_HT", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        TH2D *jesdo_MET_HT = new TH2D("jesdo_MET_HT", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        TH2D *jerdo_MET_HT = new TH2D("jerdo_MET_HT", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);

	const int nMETSlices = 3;
        double METsliceMin[nMETSlices] = {120., 200., 400.};
        double METsliceMax[nMETSlices] = {200., 400., 1e6};
        TH2D *p_HT_vs_METslice = new TH2D("p_HT_vs_METslice", "HT in MET slices; H_{T} (GeV); MET slice", nBkgHTBins, bkgHTBins, nMETSlices, 0, nMETSlices);
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(1,"MET: 120-200");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(2,"MET: 200-400");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(3,"MET: >400");
        TH2D *scaleup_HT_MET = (TH2D*)p_HT_vs_METslice->Clone("scaleup_HT_MET");
        TH2D *normup_HT_MET  = (TH2D*)p_HT_vs_METslice->Clone("normup_HT_MET");
        TH2D *jesup_HT_MET   = (TH2D*)p_HT_vs_METslice->Clone("jesup_HT_MET");
        TH2D *jesdo_HT_MET   = (TH2D*)p_HT_vs_METslice->Clone("jesdo_HT_MET");

// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
	// background from directly simulations, mctree
  	TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());

		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTGJets_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WWG_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZG_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WW_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WZ_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_ZZ_%d%s.root",RunYear,whichVFP.c_str()));
  		mctree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTJets_%d%s.root",RunYear,whichVFP.c_str()));
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
  	float HT(0);
	float nJet(0);
  	int   nBJet(0);
	
	float nISRJetFloat(0);
  	int nISRJetInt(0);

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
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("nBJet",     &nBJet);
  if (channelType == 1) mctree->SetBranchAddress("nISRJet", &nISRJetInt);
  else mctree->SetBranchAddress("nISRJet", &nISRJetFloat);
  
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
  mctree->SetBranchAddress("mcStatus",  &mcStatus);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(0);
		double scalefactorup(0);
		// scale factors
		if(channelType == 1){
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

		float weight = PUweight*XS_weight*scalefactor;
		float weight_scaleup = PUweight*XS_weight*scalefactorup;
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		// different MET, MT, lepton pT range
		double mindRpho(0.3);
		int phoIndex(0), anyphoIndex(-1);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho){mindRpho=dR; phoIndex=iMC;}
			if(dR < 0.3 && fabs((*mcPID)[iMC]) == 22)anyphoIndex = iMC;
		}
		bool isTruePho(false);
		// checking dR and momID for true photon identification
		if(mindRpho < 0.1){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )isTruePho = true;
		}
		else if(anyphoIndex >= 0){
			if((*mcPID)[anyphoIndex] == 22 && (fabs((*mcMomPID)[anyphoIndex]) <= 6 || fabs((*mcMomPID)[anyphoIndex]) == 21 || fabs((*mcMomPID)[anyphoIndex]) == 999 || fabs((*mcMomPID)[anyphoIndex])== 11 || fabs((*mcMomPID)[anyphoIndex])== 13 || fabs((*mcMomPID)[anyphoIndex])== 15 || fabs((*mcMomPID)[anyphoIndex])== 23 || fabs((*mcMomPID)[anyphoIndex])== 24)  )isTruePho = true;
		}

		bool isFSRPho(false);
		// FSR photon if particles found within dR < 0.3 around photon, and momID is e, mu, tau, W, top
		if(mindRpho < 0.3){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15 || fabs((*mcMomPID)[phoIndex])==24 || fabs((*mcMomPID)[phoIndex])==6 ))isFSRPho = true;
		}

	//	std::cout << "Type " << mcType;
	//	if(mindRpho < 0.1)std::cout << " " << (*mcPID)[phoIndex] << " " << fabs((*mcMomPID)[phoIndex]);
	//	else if(anyphoIndex >= 0)std::cout << " " << (*mcPID)[anyphoIndex] << " " << fabs((*mcMomPID)[phoIndex]);
  //  else std::cout << " no match";
	//	std::cout << std::endl;
	//	if(isTruePho)isFake = false;

		if(!isTruePho)continue;
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
		p_nBJet->Fill(nBJet, weight);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt, weight);
			p_MET_TT->Fill(sigMET,  weight);
			p_Mt_TT->Fill(sigMT,  weight);
			p_HT_TT->Fill(HT,  weight);
			p_LepPt_TT->Fill(lepPt,  weight);
                        p_nJet_TT->Fill(nJet,  weight);
                        p_nBJet_TT->Fill(nBJet, weight);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), weight);}
		
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
		if(anatype != 0){
			int htslice = -1;
  			if(HT < 100)        htslice = 0;
  			else if(HT < 400)   htslice = 1;
  			else                htslice = 2;
			p_MET_vs_HTslice->Fill(sigMET, htslice + 0.5, weight);
			scaleup_MET_HT->Fill(sigMET, htslice + 0.5, weight_scaleup);
                        jesup_MET_HT->Fill(sigMETJESup, htslice + 0.5, weight);
                        jesdo_MET_HT->Fill(sigMETJESdo, htslice + 0.5, weight);
                        jerup_MET_HT->Fill(sigMETJERup, htslice + 0.5, weight);
                        jerdo_MET_HT->Fill(sigMETJERdo, htslice + 0.5, weight);}

                if (anatype != 0) {
                        int metslice = -1;
                        if      (sigMET >= 120 && sigMET < 200) metslice = 0;
                        else if (sigMET >= 200 && sigMET < 400) metslice = 1;
                        else if (sigMET >= 400)                 metslice = 2;
                        p_HT_vs_METslice->Fill(HT, metslice + 0.5, weight);
                        scaleup_HT_MET->Fill(HT, metslice + 0.5, weight_scaleup);
                        jesup_HT_MET->Fill(HTJESup, metslice + 0.5, weight);
                        jesdo_HT_MET->Fill(HTJESdo, metslice + 0.5, weight);}

	}
	
	TH1D *h_MET_HT0 = p_MET_vs_HTslice->ProjectionX("h_MET_HT0",1,1);
	TH1D *h_MET_HT1 = p_MET_vs_HTslice->ProjectionX("h_MET_HT1",2,2);
	TH1D *h_MET_HT2 = p_MET_vs_HTslice->ProjectionX("h_MET_HT2",3,3);

        TH1D *scaleup_MET_HT_proj[3];
        scaleup_MET_HT_proj[0] = scaleup_MET_HT->ProjectionX("scaleup_HT0",1,1);
        scaleup_MET_HT_proj[1] = scaleup_MET_HT->ProjectionX("scaleup_HT1",2,2);
        scaleup_MET_HT_proj[2] = scaleup_MET_HT->ProjectionX("scaleup_HT2",3,3);

        TH1D *jesup_MET_HT_proj[3];
        jesup_MET_HT_proj[0] = jesup_MET_HT->ProjectionX("jesup_HT0",1,1);
        jesup_MET_HT_proj[1] = jesup_MET_HT->ProjectionX("jesup_HT1",2,2);
        jesup_MET_HT_proj[2] = jesup_MET_HT->ProjectionX("jesup_HT2",3,3);

        TH1D *jerup_MET_HT_proj[3];
        jerup_MET_HT_proj[0] = jerup_MET_HT->ProjectionX("jerup_HT0",1,1);
        jerup_MET_HT_proj[1] = jerup_MET_HT->ProjectionX("jerup_HT1",2,2);
        jerup_MET_HT_proj[2] = jerup_MET_HT->ProjectionX("jerup_HT2",3,3);

        TH1D *jesdo_MET_HT_proj[3];
        jesdo_MET_HT_proj[0] = jesdo_MET_HT->ProjectionX("jesdo_HT0",1,1);
        jesdo_MET_HT_proj[1] = jesdo_MET_HT->ProjectionX("jesdo_HT1",2,2);
        jesdo_MET_HT_proj[2] = jesdo_MET_HT->ProjectionX("jesdo_HT2",3,3);

        TH1D *jerdo_MET_HT_proj[3];
        jerdo_MET_HT_proj[0] = jerdo_MET_HT->ProjectionX("jerdo_HT0",1,1);
        jerdo_MET_HT_proj[1] = jerdo_MET_HT->ProjectionX("jerdo_HT1",2,2);
        jerdo_MET_HT_proj[2] = jerdo_MET_HT->ProjectionX("jerdo_HT2",3,3);

	int nTotBins = nHTSlices * nBkgMETBins;
	TH1D *h_MET_HT_concat = new TH1D("h_MET_HT_concat", "MET in HT slices; MET / HT slice index;Events", nTotBins, 0, nTotBins);
	TH1D *hists[3] = {h_MET_HT0, h_MET_HT1, h_MET_HT2};
	for(int i=0;i<3;i++){
  		for(int b=1;b<=nBkgMETBins;b++){
    			int newBin = i*nBkgMETBins + b;
    			h_MET_HT_concat->SetBinContent(newBin, hists[i]->GetBinContent(b));
    			h_MET_HT_concat->SetBinError  (newBin, hists[i]->GetBinError(b));
  		}
	}

	TH1D *p_MET_HT[3] = {h_MET_HT0, h_MET_HT1, h_MET_HT2};
        for (int htSlice = 0; htSlice < 3; htSlice++) {
                for (int ibin = 1; ibin <= nBkgMETBins; ibin++) {
                        int concatBin = htSlice * nBkgMETBins + ibin;
                        double nominal = p_MET_HT[htSlice]->GetBinContent(ibin);
                        double jeserr = std::max(fabs(jesup_MET_HT_proj[htSlice]->GetBinContent(ibin) - nominal), fabs(jesdo_MET_HT_proj[htSlice]->GetBinContent(ibin) - nominal));
                        double jererr = std::max(fabs(jerup_MET_HT_proj[htSlice]->GetBinContent(ibin) - nominal), fabs(jerdo_MET_HT_proj[htSlice]->GetBinContent(ibin) - nominal));
                        double syserror2 = 0.0;
                        syserror2 += pow(scaleup_MET_HT_proj[htSlice]->GetBinContent(ibin) - nominal, 2);
                        syserror2 += jeserr * jeserr;
                        syserror2 += jererr * jererr;
			syserror2 += pow((p_MET_HT[htSlice]->GetBinContent(ibin)*0.5), 2);
                        h_MET_HT_concat->SetBinError(concatBin, std::sqrt(syserror2));
                }
        }

	delete h_MET_HT0;
	delete h_MET_HT1;
	delete h_MET_HT2;
	delete p_MET_vs_HTslice;


	TH1D *h_HT_MET0 = p_HT_vs_METslice->ProjectionX("h_HT_MET0",1,1);
        TH1D *h_HT_MET1 = p_HT_vs_METslice->ProjectionX("h_HT_MET1",2,2);
        TH1D *h_HT_MET2 = p_HT_vs_METslice->ProjectionX("h_HT_MET2",3,3);
	
	TH1D *scaleup_HT_MET_proj[3] = {scaleup_HT_MET->ProjectionX("scaleup_MET0",1,1), scaleup_HT_MET->ProjectionX("scaleup_MET1",2,2), scaleup_HT_MET->ProjectionX("scaleup_MET2",3,3)};
        TH1D *normup_HT_MET_proj[3] = {normup_HT_MET->ProjectionX("normup_MET0",1,1), normup_HT_MET->ProjectionX("normup_MET1",2,2), normup_HT_MET->ProjectionX("normup_MET2",3,3)};
        TH1D *jesup_HT_MET_proj[3] = {jesup_HT_MET->ProjectionX("jesup_MET0",1,1), jesup_HT_MET->ProjectionX("jesup_MET1",2,2), jesup_HT_MET->ProjectionX("jesup_MET2",3,3)};
        TH1D *jesdo_HT_MET_proj[3] = {jesdo_HT_MET->ProjectionX("jesdo_MET0",1,1), jesdo_HT_MET->ProjectionX("jesdo_MET1",2,2), jesdo_HT_MET->ProjectionX("jesdo_MET2",3,3)};

	nTotBins = nMETSlices * nBkgHTBins;
        TH1D *h_HT_MET_concat = new TH1D("h_HT_MET_concat","HT in MET slices; HT / MET slice index; Events",nTotBins, 0, nTotBins);
        TH1D *p_HT_MET[3] = {h_HT_MET0, h_HT_MET1, h_HT_MET2};
        for (int i = 0; i < 3; i++) {
                for (int b = 1; b <= nBkgHTBins; b++) {
                        int newBin = i*nBkgHTBins + b;
                        h_HT_MET_concat->SetBinContent(newBin, p_HT_MET[i]->GetBinContent(b));
                        h_HT_MET_concat->SetBinError  (newBin, p_HT_MET[i]->GetBinError(b));
                }
        }

        for (int metSlice = 0; metSlice < 3; metSlice++) {
                for (int ibin = 1; ibin <= nBkgHTBins; ibin++) {
                        int concatBin = metSlice * nBkgHTBins + ibin;
                        double nominal = p_HT_MET[metSlice]->GetBinContent(ibin);
                        double jeserr = std::max(fabs(jesup_HT_MET_proj[metSlice]->GetBinContent(ibin) - nominal), fabs(jesdo_HT_MET_proj[metSlice]->GetBinContent(ibin) - nominal));
                        double err2 = pow((p_HT_MET[metSlice]->GetBinContent(ibin)*0.5), 2);
                        err2 += pow(scaleup_HT_MET_proj[metSlice]->GetBinContent(ibin) - nominal, 2);
                        err2 += jeserr * jeserr;
                        h_HT_MET_concat->SetBinError(concatBin, std::sqrt(err2));
                }
        }
        delete h_HT_MET0;
        delete h_HT_MET1;
        delete h_HT_MET2;
        delete p_HT_vs_METslice;

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((p_PhoEt->GetBinContent(ibin)*0.3),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((p_LepPt->GetBinContent(ibin)*0.3),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jerdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		syserror += pow((p_MET->GetBinContent(ibin)*0.5),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((p_Mt->GetBinContent(ibin)*0.5),2);
		// 50%  uncertainty  on cross section 
		double jeserror = max( fabs(jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jerdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)), fabs(jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow((p_HT->GetBinContent(ibin)*0.5),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow(p_dPhiEleMET->GetBinError(ibin),2);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jesdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jerdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		syserror += pow((p_dPhiEleMET->GetBinContent(ibin)*0.5),2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	

	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((p_PhoEt_TT->GetBinContent(ibin)*0.5),2);
		p_PhoEt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((p_MET_TT->GetBinContent(ibin)*0.5),2);
		p_MET_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((p_Mt_TT->GetBinContent(ibin)*0.5),2);
		p_Mt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((p_HT_TT->GetBinContent(ibin)*0.5),2);
		p_HT_TT->SetBinError(ibin,sqrt(syserror));
	}	
	//p_PhoEt->Sumw2();
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	p_dPhiEleMET_TT->Write();
	outputfile->Write();
	outputfile->Close();

}


