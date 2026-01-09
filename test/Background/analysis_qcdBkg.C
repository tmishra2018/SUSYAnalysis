#include<string>
#include "TString.h"
#include "../../include/analysis_commoncode.h"

void analysis_qcdBkg(){

	SetRunConfig();
	setTDRStyle();
	
  	gSystem->Load("../../lib/libAnaClasses.so");

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
	bool toDeriveScale(false);
	if(anatype == 0)toDeriveScale = true;

  	int channelType = ichannel; // eg = 1; mg =2;
	double factorQCD(1);
	double factorQCDUP = 1; 

	if(toDeriveScale){
		factorQCD = 1;
		factorQCDUP = 1;
	}
	else{
	        if(channelType == 1){
                        if(RunYear==2016 and preVFP == 1){
                                factorQCD = factor_egQCD_2016preVFP;
                                factorQCDUP = factor_egQCD_2016preVFP+factorerror_egQCD_2016preVFP;}
                        if(RunYear==2016 and preVFP == 0){
                                factorQCD = factor_egQCD_2016postVFP;
                                factorQCDUP = factor_egQCD_2016postVFP+factorerror_egQCD_2016postVFP;}
                        if(RunYear==2017 and preVFP == 0){
                                factorQCD = factor_egQCD_2017;
                                factorQCDUP = factor_egQCD_2017+factorerror_egQCD_2017;}
                        if(RunYear==2018 and preVFP == 0){
                                factorQCD = factor_egQCD_2018;
                                factorQCDUP = factor_egQCD_2018+factorerror_egQCD_2018;}
                }
                else if(channelType == 2){
                        if(RunYear==2016 and preVFP == 1){
                                factorQCD = factor_mgQCD_2016preVFP;
                                factorQCDUP = factor_mgQCD_2016preVFP+factorerror_mgQCD_2016preVFP;}
                        if(RunYear==2016 and preVFP == 0){
                                factorQCD = factor_mgQCD_2016postVFP;
                                factorQCDUP = factor_mgQCD_2016postVFP+factorerror_mgQCD_2016postVFP;}
                        if(RunYear==2017 and preVFP == 0){
                                factorQCD = factor_mgQCD_2017;
                                factorQCDUP = factor_mgQCD_2017+factorerror_mgQCD_2017;}
                        if(RunYear==2018 and preVFP == 0){
                                factorQCD = factor_mgQCD_2018;
                                factorQCDUP = factor_mgQCD_2018+factorerror_mgQCD_2018;}
                }
	}
	// make this ON when you have to derive the QCD and VGamma scales
	// reading scale factor from sf files
	TFile *scaleFile;
		if(channelType == 1)scaleFile = TFile::Open(Form("/eos/uscms/store/user/tmishra/fakeLep/qcd_eg_scale_%d%s.root",RunYear,whichVFP.c_str()));
		else if(channelType == 2)scaleFile = TFile::Open(Form("/eos/uscms/store/user/tmishra/fakeLep/qcd_mg_scale_%d%s.root",RunYear,whichVFP.c_str()));
	TH1D *p_scale = 0;
	if(channelType == 1)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	else if(channelType == 2)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	
	//*********** histo list ********************** 
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

	TH1D *p_dPhiEleMET_TT = new TH1D("p_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 
	TH1D *p_LepPt_TT = new TH1D("p_LepPt_TT","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_nJet_TT = new TH1D("p_nJet_TT","p_nJet",10,0,10);
	TH1D *p_nBJet_TT = new TH1D("p_nBJet_TT","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	
	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *normup_PhoEt_TT = new TH1D("normup_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta_TT = new TH1D("normup_PhoEta_TT","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt_TT = new TH1D("normup_LepPt_TT","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta_TT = new TH1D("normup_LepEta_TT","normup_LepEta",60,-3,3);
	TH1D *normup_MET_TT = new TH1D("normup_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt_TT = new TH1D("normup_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT_TT = new TH1D("normup_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET_TT = new TH1D("normup_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 
	
	TH1D *unweight_PhoEt = new TH1D("unweight_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *unweight_PhoEta = new TH1D("unweight_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *unweight_LepPt = new TH1D("unweight_LepPt","unweight_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *unweight_LepEta = new TH1D("unweight_LepEta","unweight_LepEta",60,-3,3);
	TH1D *unweight_MET = new TH1D("unweight_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *unweight_Mt = new TH1D("unweight_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *unweight_HT = new TH1D("unweight_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *unweight_dPhiEleMET = new TH1D("unweight_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	
	TH1D *unweight_PhoEt_TT = new TH1D("unweight_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *unweight_PhoEta_TT = new TH1D("unweight_PhoEta_TT","#gamma #eta; #eta;",60,-3,3);
	TH1D *unweight_LepPt_TT = new TH1D("unweight_LepPt_TT","unweight_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *unweight_LepEta_TT = new TH1D("unweight_LepEta_TT","unweight_LepEta",60,-3,3);
	TH1D *unweight_MET_TT = new TH1D("unweight_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *unweight_Mt_TT = new TH1D("unweight_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *unweight_HT_TT = new TH1D("unweight_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *unweight_dPhiEleMET_TT = new TH1D("unweight_dPhiEleMET_TT","dPhiEleMET",32,0,3.2); 

	// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
  // fake lepton is predicted from data, fakeLeptree

  if(channelType==1)fakeEtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
  if(channelType==2)fakeEtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
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
  float HT(0);
  float nJetFloat(0);
  int nJetInt(0);
  int   nBJet(0); 


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
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nBJet",     &nBJet);
  if (channelType == 1) fakeEtree->SetBranchAddress("nJet", &nJetFloat);
  else fakeEtree->SetBranchAddress("nJet", &nJetInt);


  for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		double w_qcd = 0; 
		double w_qcd_up = 0; 
		double w_qcd_unweight = 0;
		w_qcd = factorQCD;
		w_qcd_up = factorQCDUP;
		w_qcd_unweight = factorQCD;
		
	//	w_qcd = factorQCD*p_scale->GetBinContent(p_scale->FindBin(lepPt));
	//	w_qcd_up = factorQCDUP*p_scale->GetBinContent(p_scale->FindBin(lepPt));
	//	w_qcd_unweight = factorQCD;


		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		// different MET, MT, lepton pT ranges
		bool isProxy(false);
		// lepIso = 4
		if(channelType==1){if(fakeLepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;

		p_PhoEt->Fill(phoEt, w_qcd);
		p_PhoEta->Fill(phoEta, w_qcd);
		p_LepPt->Fill(lepPt, w_qcd);
		p_LepEta->Fill(lepEta, w_qcd);
		p_MET->Fill(sigMET, w_qcd);
		p_Mt->Fill(sigMT, w_qcd);
		p_HT->Fill(HT, w_qcd);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd);
		//p_nJet->Fill(nJet, w_qcd);	
		if (channelType == 1) p_nJet->Fill(nJetFloat, w_qcd);
                if (channelType == 2) p_nJet->Fill(nJetInt, w_qcd);
		p_nBJet->Fill(nBJet, w_qcd);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  w_qcd);
			p_MET_TT->Fill(sigMET,  w_qcd);
			p_Mt_TT->Fill(sigMT,  w_qcd);
			p_HT_TT->Fill(HT,  w_qcd);
			p_LepPt_TT->Fill(lepPt,  w_qcd);
                        if (channelType == 1) p_nJet_TT->Fill(nJetFloat,  w_qcd);
                        if (channelType == 2) p_nJet_TT->Fill(nJetInt,  w_qcd);
                        p_nBJet_TT->Fill(nBJet, w_qcd);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), w_qcd);
		
			normup_PhoEt_TT->Fill(phoEt, w_qcd_up);
			normup_LepPt_TT->Fill(lepPt, w_qcd_up);
			normup_MET_TT->Fill(sigMET, w_qcd_up);
			normup_Mt_TT->Fill(sigMT, w_qcd_up);
			normup_HT_TT->Fill(HT, w_qcd_up);
			normup_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), w_qcd_up);

			unweight_PhoEt_TT->Fill(phoEt, w_qcd_unweight);
			unweight_LepPt_TT->Fill(lepPt, w_qcd_unweight);
			unweight_MET_TT->Fill(sigMET, w_qcd_unweight);
			unweight_Mt_TT->Fill(sigMT, w_qcd_unweight);
			unweight_HT_TT->Fill(HT, w_qcd_unweight);
			unweight_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), w_qcd_unweight);
		}

		normup_PhoEt->Fill(phoEt, w_qcd_up);
		normup_PhoEta->Fill(phoEta,w_qcd_up);
		normup_LepPt->Fill(lepPt, w_qcd_up);
		normup_LepEta->Fill(lepEta,w_qcd_up);
		normup_MET->Fill(sigMET, w_qcd_up);
		normup_Mt->Fill(sigMT, w_qcd_up);
		normup_HT->Fill(HT, w_qcd_up);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd_up);

		unweight_PhoEt->Fill(phoEt, w_qcd_unweight);
		unweight_PhoEta->Fill(phoEta,w_qcd_unweight);
		unweight_LepPt->Fill(lepPt, w_qcd_unweight);
		unweight_LepEta->Fill(lepEta,w_qcd_unweight);
		unweight_MET->Fill(sigMET, w_qcd_unweight);
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
	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt_TT->GetBinError(ibin)* p_PhoEt_TT->GetBinError(ibin);
		syserror += pow((normup_PhoEt_TT->GetBinContent(ibin)-p_PhoEt_TT->GetBinContent(ibin)),2);
		syserror += pow((unweight_PhoEt_TT->GetBinContent(ibin)-p_PhoEt_TT->GetBinContent(ibin)),2);
		p_PhoEt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((unweight_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt_TT->GetBinError(ibin)* p_LepPt_TT->GetBinError(ibin);
		syserror += pow((normup_LepPt_TT->GetBinContent(ibin)-p_LepPt_TT->GetBinContent(ibin)),2);
		syserror += pow((unweight_LepPt_TT->GetBinContent(ibin)-p_LepPt_TT->GetBinContent(ibin)),2);
		p_LepPt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);// stat
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2); // error of normalisation scale and lepton pt scaling
		syserror += pow((unweight_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2); // central value w/o error of normalisation and no pt scaling
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET_TT->GetBinError(ibin)* p_MET_TT->GetBinError(ibin);// stat
		syserror += pow((normup_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)),2); // error of normalisation scale and lepton pt scaling
		syserror += pow((unweight_MET_TT->GetBinContent(ibin)-p_MET_TT->GetBinContent(ibin)),2); // central value w/o error of normalisation and no pt scaling
		p_MET_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((unweight_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt_TT->GetBinError(ibin)* p_Mt_TT->GetBinError(ibin);
		syserror += pow((normup_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)),2);
		syserror += pow((unweight_Mt_TT->GetBinContent(ibin)-p_Mt_TT->GetBinContent(ibin)),2);
		p_Mt_TT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((unweight_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT_TT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT_TT->GetBinError(ibin)* p_HT_TT->GetBinError(ibin);
		syserror += pow((normup_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)),2);
		syserror += pow((unweight_HT_TT->GetBinContent(ibin)-p_HT_TT->GetBinContent(ibin)),2);
		p_HT_TT->SetBinError(ibin,sqrt(syserror));
	}	

	std::ostringstream outputname;
	outputname << "/eos/uscms/store/user/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_qcd";
	else if(channelType==2)outputname << "mg_qcd";
	if(anatype ==0 or anatype ==1) outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
        outputname <<"_" << RunYear<<whichVFP <<".root";
	
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
	p_nBJet->Write();
	p_PhoEt_TT->Write();
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	p_dPhiEleMET_TT->Write();
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	outputfile->Write();
	outputfile->Close();
}
