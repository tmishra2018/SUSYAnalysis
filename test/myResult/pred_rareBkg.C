#include "../../include/analysis_commoncode.h"
void pred_rareBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
	
	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
	else whichVFP = "";


  	gSystem->Load("../../lib/libAnaClasses.so");
	esfScaleFactor  objectESF;

  	int channelType = ichannel; // eg = 1; mg =2;
	//*********** histo list **********************//
	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_rareBkg";
	else if(channelType==2)outputname << "mg_rareBkg";
	if(anatype ==0)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << "_" <<RunYear<<whichVFP<<".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	std::ostringstream histname;

	TH1D *p_PhoEt = new TH1D("p_PhoEt","; p_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nSigPtBins,sigPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	TH1D *h_rare_norm = 0;
	TH1D *h_rare_control = 0;
	TH1D *h_rare_jesUp = 0;
	TH1D *h_rare_jesDown = 0;
	TH1D *h_rare_jerUp = 0;
	TH1D *h_rare_jerDown = 0;
	TH1D *h_rare_esfUp = 0;
	TH1D *h_rare_syserr_jes = 0;
	TH1D *h_rare_syserr_jer = 0;
	TH1D *h_rare_syserr_esf = 0;
	TH1D *h_rare_syserr_scale = 0;
	TH1D *h_rare_syserr_eleshape = 0;
	TH1D *h_rare_syserr_jetshape = 0;
	TH1D *h_rare_syserr_qcdshape = 0;
	TH1D *h_rare_syserr_xs = 0;
	TH1D *h_rare_syserr_lumi = 0;
	TH1D *h_rare_syserr_isr = 0;
	if(channelType==1){
		h_rare_norm            = new TH1D("eg_rare_norm","eventcount",NBIN,0,NBIN);
		h_rare_control         = new TH1D("eg_rare_control","eventcount",NBIN,0,NBIN);
		h_rare_jesUp           = new TH1D("eg_rare_jesUp","eventcount",NBIN,0,NBIN);
		h_rare_jesDown         = new TH1D("eg_rare_jesDown","eventcount",NBIN,0,NBIN);
		h_rare_jerUp           = new TH1D("eg_rare_jerUp","eventcount",NBIN,0,NBIN);
		h_rare_jerDown         = new TH1D("eg_rare_jerDown","eventcount",NBIN,0,NBIN);
		h_rare_esfUp           = new TH1D("eg_rare_esfUp","eventcount",NBIN,0,NBIN);
		h_rare_syserr_jes      = new TH1D("eg_rare_syserr_jes","",NBIN,0,NBIN);	
		h_rare_syserr_jer      = new TH1D("eg_rare_syserr_jer","",NBIN,0,NBIN);	
		h_rare_syserr_esf      = new TH1D("eg_rare_syserr_esf","",NBIN,0,NBIN);	
		h_rare_syserr_scale    = new TH1D("eg_rare_syserr_scale","",NBIN,0,NBIN);	
		h_rare_syserr_eleshape = new TH1D("eg_rare_syserr_eleshape","",NBIN,0,NBIN);	
		h_rare_syserr_jetshape = new TH1D("eg_rare_syserr_jetshape","",NBIN,0,NBIN);	
		h_rare_syserr_qcdshape = new TH1D("eg_rare_syserr_qcdshape","",NBIN,0,NBIN);	
		h_rare_syserr_xs       = new TH1D("eg_rare_syserr_xs","",NBIN,0,NBIN);	
		h_rare_syserr_lumi     = new TH1D("eg_rare_syserr_lumi","",NBIN,0,NBIN);
		h_rare_syserr_isr      = new TH1D("eg_rare_syserr_isr","",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_rare_norm            = new TH1D("mg_rare_norm","eventcount",NBIN,0,NBIN);
		h_rare_control         = new TH1D("mg_rare_control","eventcount",NBIN,0,NBIN);
		h_rare_jesUp           = new TH1D("mg_rare_jesUp","eventcount",NBIN,0,NBIN);
		h_rare_jesDown         = new TH1D("mg_rare_jesDown","eventcount",NBIN,0,NBIN);
		h_rare_jerUp           = new TH1D("mg_rare_jerUp","eventcount",NBIN,0,NBIN);
		h_rare_jerDown         = new TH1D("mg_rare_jerDown","eventcount",NBIN,0,NBIN);
		h_rare_esfUp           = new TH1D("mg_rare_esfUp","eventcount",NBIN,0,NBIN);
		h_rare_syserr_jes      = new TH1D("mg_rare_syserr_jes","",NBIN,0,NBIN);	
		h_rare_syserr_jer      = new TH1D("mg_rare_syserr_jer","",NBIN,0,NBIN);	
		h_rare_syserr_esf      = new TH1D("mg_rare_syserr_esf","",NBIN,0,NBIN);	
		h_rare_syserr_scale    = new TH1D("mg_rare_syserr_scale","",NBIN,0,NBIN);	
		h_rare_syserr_eleshape = new TH1D("mg_rare_syserr_eleshape","",NBIN,0,NBIN);	
		h_rare_syserr_jetshape = new TH1D("mg_rare_syserr_jetshape","",NBIN,0,NBIN);	
		h_rare_syserr_qcdshape = new TH1D("mg_rare_syserr_qcdshape","",NBIN,0,NBIN);	
		h_rare_syserr_xs       = new TH1D("mg_rare_syserr_xs","",NBIN,0,NBIN);	
		h_rare_syserr_lumi     = new TH1D("mg_rare_syserr_lumi","",NBIN,0,NBIN);
		h_rare_syserr_isr      = new TH1D("mg_rare_syserr_isr","",NBIN,0,NBIN);
	} 

	TH1D *jesup_MET = new TH1D("jesup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jesup_Mt = new TH1D("jesup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *jesup_HT = new TH1D("jesup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *jesup_dPhiEleMET = new TH1D("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jesdo_MET = new TH1D("jesdo_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jesdo_Mt = new TH1D("jesdo_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jesdo_HT = new TH1D("jesdo_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *jesdo_dPhiEleMET = new TH1D("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerup_MET = new TH1D("jerup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jerup_Mt = new TH1D("jerup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jerup_dPhiEleMET = new TH1D("jerup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerdo_MET = new TH1D("jerdo_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jerdo_Mt = new TH1D("jerdo_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jerdo_dPhiEleMET = new TH1D("jerdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *scaleup_PhoEt = new TH1D("scaleup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *scaleup_PhoEta = new TH1D("scaleup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *scaleup_LepPt = new TH1D("scaleup_LepPt","scaleup_LepPt",nSigPtBins,sigPtBins);
	TH1D *scaleup_LepEta = new TH1D("scaleup_LepEta","scaleup_LepEta",60,-3,3);
	TH1D *scaleup_MET = new TH1D("scaleup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *scaleup_Mt = new TH1D("scaleup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *scaleup_HT = new TH1D("scaleup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *scaleup_dPhiEleMET = new TH1D("scaleup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  	TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
	
/*	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str())); //
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_%d%s.root",RunYear,whichVFP.c_str())); //
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str())); //
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str())); //
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str())); //
*/
//      mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root");
//      mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG50_VetoEle.root");
//      mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG130_VetoEle.root");
//      mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_ZG_VetoEle.root");
//      mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_DY.root");

	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTGJets_%d%s.root",RunYear,whichVFP.c_str())); // looks correct
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTJets_%d%s.root",RunYear,whichVFP.c_str())); // 25% more
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WWG_%d%s.root",RunYear,whichVFP.c_str())); // 10 times larger
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZG_%d%s.root",RunYear,whichVFP.c_str()));   //  6 times large
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WW_%d%s.root",RunYear,whichVFP.c_str()));     // ~ 0 entries
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_%d%s.root",RunYear,whichVFP.c_str()));     // 0 entries

/*
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_TTG_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WWG_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WZG_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WW_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WZ_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_TT_VetoEle.root");
*/
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
	int   nISRJet(0);
	float ISRPt(0);
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
  	mctree->SetBranchAddress("nISRJet",   &nISRJet);
  	mctree->SetBranchAddress("ISRJetPt",     &ISRPt);
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
                        else if(RunYear == 2018)                        XS_weight = lumi_2018_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 678)                         XS_weight = lumi_678_DoubleEG*1000*crosssection/ntotalevent;}


                else if(channelType == 2){
                        if(RunYear == 2016 and preVFP == 1)             XS_weight = lumi_2016preVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        XS_weight = lumi_2016postVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2017)                        XS_weight = lumi_2017_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2018)                        XS_weight = lumi_2018_MuonEG*1000*crosssection/ntotalevent;
		        else if(RunYear == 678)                         XS_weight = lumi_678_MuonEG*1000*crosssection/ntotalevent;}

		
		//XS_weight = (lumi_2016preVFP_DoubleEG+lumi_2016postVFP_DoubleEG)*1000*crosssection/ntotalevent;
		//crosssection = 0.8099;
		
		double weight = PUweight*XS_weight*scalefactor;
		double weight_scaleup = PUweight*XS_weight*scalefactorup;
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		double mindRpho(0.3);
		int phoIndex(0), anyphoIndex(-1);
		int posW(0), negW(0), posZ(0);
		int decaytype(-1);
		bool phoCloseT(false);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho){mindRpho=dR; phoIndex=iMC;}
			if(dR < 0.3 && fabs((*mcPID)[iMC]) == 22)anyphoIndex = iMC;
			if(dR < 0.3 && fabs((*mcPID)[iMC]) == 5 && fabs((*mcMomPID)[iMC]) == 6)phoCloseT=true;

			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == 24 )posW = 1;
			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == -24)negW = 1;
			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == 23 )posZ = 1;
		}
		if( (posW == 1 && negW ==1) || ( posW == 1 && posZ == 1) || (negW == 1 && posZ == 1))decaytype = 1; 

		bool isTruePho(false);
		if(mindRpho < 0.1){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )isTruePho = true;
		}
		else if(anyphoIndex >= 0){
			if((*mcPID)[anyphoIndex] == 22 && (fabs((*mcMomPID)[anyphoIndex]) <= 6 || fabs((*mcMomPID)[anyphoIndex]) == 21 || fabs((*mcMomPID)[anyphoIndex]) == 999 || fabs((*mcMomPID)[anyphoIndex])== 11 || fabs((*mcMomPID)[anyphoIndex])== 13 || fabs((*mcMomPID)[anyphoIndex])== 15 || fabs((*mcMomPID)[anyphoIndex])== 23 || fabs((*mcMomPID)[anyphoIndex])== 24)  )isTruePho = true;
		}

		bool isFSRPho(false);
		if(mindRpho < 0.3){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15 || fabs((*mcMomPID)[phoIndex])==24 || fabs((*mcMomPID)[phoIndex])==6 ))isFSRPho = true;
			if((*mcPID)[phoIndex] == 22 &&  fabs((*mcMomPID)[phoIndex])==5 && phoCloseT )isFSRPho = true;
		}

		if(!isTruePho)continue;
		if(mcType ==11 && !isFSRPho)continue;
		if((mcType == 12 || mcType == 13) &&  decaytype < 1)continue;

		if(sigMT >= lowMt && (highMt < 0 || (highMt > 0 && sigMT < highMt)) ){

			p_MET->Fill(sigMET, weight);
			scaleup_MET->Fill(sigMET, weight_scaleup);

			/** cut flow *****/
			if(sigMET >= lowMET && (highMET < 0|| (highMET > 0 && sigMET < highMET))){

				p_PhoEt->Fill(phoEt, weight);
				p_PhoEta->Fill(phoEta,weight);
				p_LepPt->Fill(lepPt, weight);
				p_LepEta->Fill(lepEta,weight);
				p_Mt->Fill(sigMT, weight);
				p_HT->Fill(HT, weight);
				p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
				p_nJet->Fill( nJet, weight);

				int SigBinIndex(-1);
				SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt);
				if( SigBinIndex >= 0){ 
					  h_rare_norm->Fill( SigBinIndex, weight);
					  h_rare_control->Fill( SigBinIndex);
					  h_rare_esfUp->Fill( SigBinIndex, weight_scaleup);
				}

				scaleup_PhoEt->Fill(phoEt, weight_scaleup);
				scaleup_PhoEta->Fill(phoEta,weight_scaleup);
				scaleup_LepPt->Fill(lepPt, weight_scaleup);
				scaleup_LepEta->Fill(lepEta,weight_scaleup);
				scaleup_Mt->Fill(sigMT, weight_scaleup);
				scaleup_HT->Fill(HT, weight_scaleup);
				scaleup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaleup);
			}
		}
		if(sigMTJESup >= lowMt && (highMt < 0 || (highMt > 0 && sigMTJESup < highMt)) ){
			if(sigMETJESup >= lowMET && (highMET < 0|| (highMET > 0 && sigMETJESup < highMET))){
				if(Bin.findSignalBin(sigMETJESup, HTJESup, phoEt)>=0)h_rare_jesUp->Fill( Bin.findSignalBin(sigMETJESup, HTJESup, phoEt),  weight);
				jesup_MET->Fill(sigMETJESup, weight);
				jesup_Mt->Fill(sigMTJESup, weight);
				jesup_HT->Fill(HTJESup, weight);
				jesup_dPhiEleMET->Fill(fabs(dPhiLepMETJESup), weight);
			}
		}
		if(sigMTJESdo >= lowMt && (highMt < 0 || (highMt > 0 && sigMTJESdo < highMt)) ){
			if(sigMETJESdo >= lowMET && (highMET < 0|| (highMET > 0 && sigMETJESdo < highMET))){
				if(Bin.findSignalBin(sigMETJESdo, HTJESdo, phoEt)>=0)h_rare_jesDown->Fill( Bin.findSignalBin(sigMETJESdo, HTJESdo, phoEt),  weight);
				jesdo_MET->Fill(sigMETJESdo, weight);
				jesdo_Mt->Fill(sigMTJESdo, weight);
				jesdo_HT->Fill(HTJESdo, weight);
				jesdo_dPhiEleMET->Fill(fabs(dPhiLepMETJESdo), weight);
			}
		}
		if(sigMTJERup >= lowMt && (highMt < 0 || (highMt > 0 && sigMTJERup < highMt)) ){
			if(sigMETJERup >= lowMET && (highMET < 0|| (highMET > 0 && sigMETJERup < highMET))){
				if(Bin.findSignalBin(sigMETJERup, HT, phoEt)>=0)h_rare_jerUp->Fill( Bin.findSignalBin(sigMETJERup, HT, phoEt),       weight);
				jerup_MET->Fill(sigMETJERup, weight);
				jerup_Mt->Fill(sigMTJERup, weight);
				jerup_dPhiEleMET->Fill(fabs(dPhiLepMETJERup), weight);
			}
		}
		if(sigMTJERdo >= lowMt && (highMt < 0 || (highMt > 0 && sigMTJERdo < highMt)) ){
			if(sigMETJERdo >= lowMET && (highMET < 0|| (highMET > 0 && sigMETJERdo < highMET))){
				if(Bin.findSignalBin(sigMETJERdo, HT, phoEt)>=0)h_rare_jerDown->Fill( Bin.findSignalBin(sigMETJERdo, HT, phoEt),       weight); 
					jerdo_MET->Fill(sigMETJERdo, weight);
					jerdo_Mt->Fill(sigMTJERdo, weight);
					jerdo_dPhiEleMET->Fill(fabs(dPhiLepMETJERdo), weight);
			}
		}

	}

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
	
	//h_rare_norm->Sumw2();
	float binSum=0;
	for(int contbin(1); contbin <=NBIN; contbin++){
		binSum= binSum+h_rare_norm->GetBinContent(contbin);
		float nominalsig = h_rare_norm->GetBinContent(contbin); 
	//	std::cout << " bin " << contbin << " rare " << nominalsig << " error " << h_rare_norm->GetBinError(contbin) << " rel error " << h_rare_norm->GetBinError(contbin)/h_rare_norm->GetBinContent(contbin) << " rel error in control region " << sqrt(h_rare_control->GetBinContent(contbin))/h_rare_control->GetBinContent(contbin) << std::endl;
		float jesuperror = fabs(h_rare_jesUp->GetBinContent(contbin)- nominalsig);
		float jesdoerror = fabs(h_rare_jesDown->GetBinContent(contbin)- nominalsig);
		float jeruperror = fabs(h_rare_jerUp->GetBinContent(contbin)- nominalsig);
		float jerdoerror = fabs(h_rare_jerDown->GetBinContent(contbin)- nominalsig);
		float esfuperror = fabs(h_rare_esfUp->GetBinContent(contbin)- nominalsig); 	
	//	std::cout << " jesup " << jesuperror << " jesdo " << jesdoerror << " jerup " << jeruperror << " jerdo " << jerdoerror << std::endl;

		h_rare_syserr_jes->SetBinContent(contbin, max(jesuperror,jesdoerror));
		h_rare_syserr_jer->SetBinContent(contbin, max(jeruperror,jerdoerror));
		h_rare_syserr_esf->SetBinContent(contbin, esfuperror);
		h_rare_syserr_scale->SetBinContent(contbin, -1);
		h_rare_syserr_eleshape->SetBinContent(contbin, -1);
		h_rare_syserr_jetshape->SetBinContent(contbin, -1);
		h_rare_syserr_qcdshape->SetBinContent(contbin, -1);
		h_rare_syserr_xs->SetBinContent(contbin, 0.5*h_rare_norm->GetBinContent(contbin));     
		h_rare_syserr_lumi->SetBinContent(contbin, 0.026*h_rare_norm->GetBinContent(contbin));      
		h_rare_syserr_isr->SetBinContent(contbin, -1);      
	}
	
	cout<<"total entries : "<<binSum<<endl;	
	outputfile->Write();
	outputfile->Close();

}
