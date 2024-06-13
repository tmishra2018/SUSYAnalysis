#include "../../include/analysis_commoncode.h"

struct runinfo{
	int runN;
	int lumiN;
	Long_t eventN;
	int binN;
};

bool compareByRun(const runinfo &a, const runinfo &b)
{
		if(a.runN < b.runN)
			return true;
		else if(a.runN == b.runN)
			return a.lumiN < b.lumiN;
		else return false;
}


void pred_qcdBkg(){

	bool toDeriveScale(false);

	SetSignalConfig();
	//binning Bin(NBIN, METbin1, METbin2, METbin3, HTbin1, HTbin2, HTbin3, PHOETbin, PHOETBin2);
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  gSystem->Load("../../lib/libAnaClasses.so");

  int channelType = ichannel; // eg = 1; mg =2;
	double factorQCD(1);
	double factorQCDUP = 1; 

	if(channelType == 1){
		factorQCD = factor_egQCD;
		factorQCDUP = factor_egQCD + factorerror_egQCD;
	}
	else if(channelType == 2){
		factorQCD = factor_mgQCD;
		factorQCDUP = factor_mgQCD + factorerror_mgQCD;
	}
	if(anatype == 0){
		toDeriveScale = true;
		factorQCD = 1;
		factorQCDUP = 1;
	}

	// Corrections on electron proxy sample as described in Section 5.3.3 in AN
	TFile *scaleFile;
	if(channelType == 1)scaleFile = TFile::Open("qcd_eg_scale.root");
	else if(channelType == 2)scaleFile = TFile::Open("qcd_mg_scale.root");
	TH1D *p_scale = 0;
	if(channelType == 1)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	else if(channelType == 2)p_scale = (TH1D*)scaleFile->Get("transfer_factor");

	//*********** histo list **********************//
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

	TH1D *h_qcdfakelep_norm =0;
	TH1D *h_qcdfakelep_controlsample =0;
	TH1D *h_qcdfakelep_transferfactor =0;
	TH1D *h_qcdfakelep_normup =0;
	TH1D *h_qcdfakelep_unweight =0;
	TH1D *h_qcdfakelep_syserr_jes =0;
	TH1D *h_qcdfakelep_syserr_jer =0;
	TH1D *h_qcdfakelep_syserr_esf =0;
	TH1D *h_qcdfakelep_syserr_scale =0;
	TH1D *h_qcdfakelep_syserr_e_to_pho =0;
	TH1D *h_qcdfakelep_syserr_j_to_pho =0;
	TH1D *h_qcdfakelep_syserr_j_to_lep =0;
	TH1D *h_qcdfakelep_syserr_xs =0;
	TH1D *h_qcdfakelep_syserr_lumi =0;
	TH1D *h_qcdfakelep_syserr_isr =0;
	if(channelType==1){
		h_qcdfakelep_norm            = new TH1D("eg_qcdfakelep_norm","eventcount",NBIN,0,NBIN);
		h_qcdfakelep_controlsample   = new TH1D("eg_qcdfakelep_controlsample","",NBIN,0,NBIN);
		h_qcdfakelep_transferfactor  = new TH1D("eg_qcdfakelep_transferfactor","",NBIN,0,NBIN);
		h_qcdfakelep_normup          = new TH1D("eg_qcdfakelep_normup","eventcount; eventcount (GeV);",NBIN,0,NBIN);
		h_qcdfakelep_unweight        = new TH1D("eg_qcdfakelep_unweight","eventcount; eventcount (GeV);",NBIN,0,NBIN);
		h_qcdfakelep_syserr_jes      = new TH1D("eg_qcdfakelep_syserr_jes","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_jer      = new TH1D("eg_qcdfakelep_syserr_jer","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_esf      = new TH1D("eg_qcdfakelep_syserr_esf","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_scale    = new TH1D("eg_qcdfakelep_syserr_scale","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_e_to_pho = new TH1D("eg_qcdfakelep_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_j_to_pho = new TH1D("eg_qcdfakelep_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_j_to_lep = new TH1D("eg_qcdfakelep_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_xs       = new TH1D("eg_qcdfakelep_syserr_xs","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_lumi     = new TH1D("eg_qcdfakelep_syserr_lumi","",NBIN,0,NBIN);
		h_qcdfakelep_syserr_isr      = new TH1D("eg_qcdfakelep_syserr_isr","",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_qcdfakelep_norm            = new TH1D("mg_qcdfakelep_norm","eventcount",NBIN,0,NBIN);
		h_qcdfakelep_controlsample   = new TH1D("mg_qcdfakelep_controlsample","",NBIN,0,NBIN);
		h_qcdfakelep_transferfactor  = new TH1D("mg_qcdfakelep_transferfactor","",NBIN,0,NBIN);
		h_qcdfakelep_normup          = new TH1D("mg_qcdfakelep_normup","eventcount; eventcount (GeV);",NBIN,0,NBIN);
		h_qcdfakelep_unweight        = new TH1D("mg_qcdfakelep_unweight","eventcount; eventcount (GeV);",NBIN,0,NBIN);
		h_qcdfakelep_syserr_jes      = new TH1D("mg_qcdfakelep_syserr_jes","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_jer      = new TH1D("mg_qcdfakelep_syserr_jer","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_esf      = new TH1D("mg_qcdfakelep_syserr_esf","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_scale    = new TH1D("mg_qcdfakelep_syserr_scale","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_e_to_pho = new TH1D("mg_qcdfakelep_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_j_to_pho = new TH1D("mg_qcdfakelep_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_j_to_lep = new TH1D("mg_qcdfakelep_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_xs       = new TH1D("mg_qcdfakelep_syserr_xs","",NBIN,0,NBIN);	
		h_qcdfakelep_syserr_lumi     = new TH1D("mg_qcdfakelep_syserr_lumi","",NBIN,0,NBIN);
		h_qcdfakelep_syserr_isr      = new TH1D("mg_qcdfakelep_syserr_isr","",NBIN,0,NBIN);
	} 
	
	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nSigPtBins,sigPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *unweight_PhoEt = new TH1D("unweight_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *unweight_PhoEta = new TH1D("unweight_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *unweight_LepPt = new TH1D("unweight_LepPt","unweight_LepPt",nSigPtBins,sigPtBins);
	TH1D *unweight_LepEta = new TH1D("unweight_LepEta","unweight_LepEta",60,-3,3);
	TH1D *unweight_MET = new TH1D("unweight_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *unweight_Mt = new TH1D("unweight_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *unweight_HT = new TH1D("unweight_HT","HT; HT (GeV);",nSigHTBins, sigHTBins);
	TH1D *unweight_dPhiEleMET = new TH1D("unweight_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	std::vector<runinfo> sig_runV;
	sig_runV.clear();
	
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal_newEta.root");
        if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");

  int   run(0);
  Long_t event(0);
  int   lumis(0);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
//	float threeMass(0);
  float HT(0);
  float nJet(0);
  
  fakeEtree->SetBranchAddress("run",       &run);	
  fakeEtree->SetBranchAddress("event",     &event);
  fakeEtree->SetBranchAddress("lumis",     &lumis);
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
//	fakeEtree->SetBranchAddress("threeMass", &threeMass);
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nJet",      &nJet);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		//if(ievt%1000 ==0)std::cout <<"event " << ievt << std::endl;

		double w_qcd = 1.; 
		double w_qcd_up = 1.; 
		double w_qcd_unweight = 1.;

		// Corrections on electron proxy sample, not on muon proxy sample as described in Section 5.3.3 in AN
		//if(channelType == 1){ // eg channel
		//	w_qcd = factorQCD*p_scale->GetBinContent(p_scale->FindBin(lepPt));
		//	w_qcd_up = factorQCDUP*p_scale->GetBinContent(p_scale->FindBin(lepPt));
		//	w_qcd_unweight = factorQCD;
		//}
		//else{
		//	w_qcd = factorQCD;
		//	w_qcd_up = factorQCDUP;
		//	w_qcd_unweight = factorQCD;
		//}
	
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		
		bool isProxy(false);
		if(channelType==1){if(fakeLepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;

		p_MET->Fill(sigMET, w_qcd);
		normup_MET->Fill(sigMET, w_qcd_up);
		unweight_MET->Fill(sigMET, w_qcd_unweight);

		/** cut flow *****/
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		p_PhoEt->Fill(phoEt, w_qcd);
		p_PhoEta->Fill(phoEta, w_qcd);
		p_LepPt->Fill(lepPt, w_qcd);
		p_LepEta->Fill(lepEta, w_qcd);
		p_Mt->Fill(sigMT, w_qcd);
		p_HT->Fill(HT, w_qcd);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd);
		p_nJet->Fill(nJet, w_qcd);	
	
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			runinfo sig_runinfo;
			sig_runinfo.runN = run;
			sig_runinfo.eventN = event;
			sig_runinfo.lumiN = lumis;
			sig_runinfo.binN = SigBinIndex;
			sig_runV.push_back(sig_runinfo);
			h_qcdfakelep_norm->Fill( SigBinIndex, w_qcd);
			h_qcdfakelep_controlsample->Fill( SigBinIndex );
			h_qcdfakelep_normup->Fill( SigBinIndex, w_qcd_up); 
			h_qcdfakelep_unweight->Fill( SigBinIndex, w_qcd_unweight);
		}

		normup_PhoEt->Fill(phoEt, w_qcd_up);
		normup_PhoEta->Fill(phoEta,w_qcd_up);
		normup_LepPt->Fill(lepPt, w_qcd_up);
		normup_LepEta->Fill(lepEta,w_qcd_up);
		normup_Mt->Fill(sigMT, w_qcd_up);
		normup_HT->Fill(HT, w_qcd_up);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd_up);

		unweight_PhoEt->Fill(phoEt, w_qcd_unweight);
		unweight_PhoEta->Fill(phoEta,w_qcd_unweight);
		unweight_LepPt->Fill(lepPt, w_qcd_unweight);
		unweight_LepEta->Fill(lepEta,w_qcd_unweight);
		unweight_Mt->Fill(sigMT, w_qcd_unweight);
		unweight_HT->Fill(HT, w_qcd_unweight);
		unweight_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd_unweight);
	}

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((unweight_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((unweight_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((unweight_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((unweight_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((unweight_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < h_qcdfakelep_norm->GetSize(); ibin++){
		if(h_qcdfakelep_controlsample->GetBinContent(ibin) > 0){
			float normerror = fabs(h_qcdfakelep_normup->GetBinContent(ibin) - h_qcdfakelep_norm->GetBinContent(ibin)); 
			float j_to_leperror = fabs(h_qcdfakelep_unweight->GetBinContent(ibin) - h_qcdfakelep_norm->GetBinContent(ibin));

			h_qcdfakelep_transferfactor->SetBinContent(ibin,h_qcdfakelep_norm->GetBinContent(ibin)/h_qcdfakelep_controlsample->GetBinContent(ibin)); 
			h_qcdfakelep_transferfactor->SetBinError(ibin,normerror/h_qcdfakelep_norm->GetBinContent(ibin));
			h_qcdfakelep_syserr_jes->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_jer->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_esf->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_scale->SetBinContent(ibin, normerror);
			h_qcdfakelep_syserr_e_to_pho->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_j_to_pho->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_j_to_lep->SetBinContent(ibin, j_to_leperror);
			h_qcdfakelep_syserr_xs->SetBinContent(ibin, -1);     
			h_qcdfakelep_syserr_lumi->SetBinContent(ibin, -1);      
			h_qcdfakelep_syserr_isr->SetBinContent(ibin, -1);     
		} 
		else{	
			float normerror = fabs(h_qcdfakelep_normup->GetBinContent(ibin) - h_qcdfakelep_norm->GetBinContent(ibin)); 
			float j_to_leperror = fabs(h_qcdfakelep_unweight->GetBinContent(ibin) - h_qcdfakelep_norm->GetBinContent(ibin));

			h_qcdfakelep_transferfactor->SetBinContent(ibin,factorQCD);
			h_qcdfakelep_transferfactor->SetBinError(ibin, (factorQCDUP-factorQCD)/factorQCD);
			h_qcdfakelep_syserr_jes->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_jer->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_esf->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_scale->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_e_to_pho->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_j_to_pho->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_j_to_lep->SetBinContent(ibin, -1);
			h_qcdfakelep_syserr_xs->SetBinContent(ibin, -1);     
			h_qcdfakelep_syserr_lumi->SetBinContent(ibin, -1);      
			h_qcdfakelep_syserr_isr->SetBinContent(ibin, -1);     
		} 
	}	

	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_qcd";
	else if(channelType==2)outputname << "mg_qcd";
	if(anatype ==0)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << "_iso" << lepIso;
	outputname << ".root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	h_qcdfakelep_norm->Write();             
	h_qcdfakelep_controlsample->Write();
	h_qcdfakelep_transferfactor->Write();
	h_qcdfakelep_syserr_jes->Write();       
	h_qcdfakelep_syserr_jer->Write();       
	h_qcdfakelep_syserr_esf->Write();       
	h_qcdfakelep_syserr_scale->Write();     
	h_qcdfakelep_syserr_e_to_pho->Write();  
	h_qcdfakelep_syserr_j_to_pho->Write();  
	h_qcdfakelep_syserr_j_to_lep->Write();  
	h_qcdfakelep_syserr_xs->Write();        
	h_qcdfakelep_syserr_lumi->Write();      
	h_qcdfakelep_syserr_isr->Write();     

	float binSum=0;
        for(int ibin(1); ibin < h_qcdfakelep_controlsample->GetSize(); ibin++) binSum= binSum+h_qcdfakelep_norm->GetBinContent(ibin);

	//for(int ibin(1); ibin < h_qcdfakelep_controlsample->GetSize(); ibin++)std::cout <<"bin " << ibin << " " <<  h_qcdfakelep_controlsample->GetBinContent(ibin) << std::endl; 
	std::sort(sig_runV.begin(), sig_runV.end(), compareByRun);
//	for(unsigned i=0; i < sig_runV.size(); i++)std::cout <<  sig_runV[i].binN+NBIN << " " <<  sig_runV[i].runN << " " << sig_runV[i].lumiN << " " << sig_runV[i].eventN << std::endl;	
	cout<<"total entries : "<<binSum<<endl;
	outputfile->Write();
	outputfile->Close();
}


