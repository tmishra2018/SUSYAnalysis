// root -l -q "closure_jetfakepho.C+(1, 2016, 1)"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TROOT.h"
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
#include "TGraphErrors.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_commoncode.h"

//#define NTOY 1000
#define NTOY 1

#define MAXET 199
#define MAXMT 399
#define MAXMET 399
#define MAXHT 399

void closure_jetfakepho(int ichannel, int RunYear, bool preVFP){

  gSystem->Load("../../lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;

  int lepPtCut = 25;
  if (ichannel == 2) lepPtCut = 20;

  std::string whichVFP;
  if(RunYear==2016 and preVFP == true) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == false) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
  
  gStyle->SetTitleXOffset(2.5);

  gROOT->SetBatch(kTRUE);
  gRandom = new TRandom3(0);
  gRandom->SetSeed(0);
	double randomweight_jet[1000];
	randomweight_jet[0] = 0;
	for(unsigned ir(1); ir<1000; ir++)	
		randomweight_jet[ir] = -1+ gRandom->Rndm()*2.0;

	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,1000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,1000,4);
//	double jetfake_numerror[264];
//	double jetfake_denerror[264];
	
	std::ifstream jetfakefile;
	if(channelType == 1)
		jetfakefile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EB.txt",RunYear,whichVFP.c_str()));
	if(channelType == 2)
		jetfakefile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EB.txt",RunYear,whichVFP.c_str()));
	
	std::string paratype;
	float paravalue;	
	if(jetfakefile.is_open()){
		for(int i(0); i < 4; i++){
			jetfakefile >> paratype >> paravalue;
			fitfunc_den->SetParameter(i, paravalue);
		}
		for(int i(0); i < 4; i++){
			jetfakefile >> paratype >> paravalue;
			fitfunc_num->SetParameter(i, paravalue);
		}
	}
//	int binnumber;
//	for(int i(0); i < 264; i++){
//		jetfakefile >> paratype >> binnumber >> paravalue;
//		jetfake_numerror[i] = paravalue/2;
//	}
//	for(int i(0); i < 264; i++){
//		jetfakefile >> paratype >> binnumber >>  paravalue;
//		jetfake_denerror[i] = paravalue/2;
//	}

	// Signal Tree //
	//*********** histo list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",11,35,200);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",11,35,200);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",10,0,200);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",10,0,200);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",10,0,400);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",9,0,9);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	// signal events directly from simulation
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));

	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));

	float crosssection(0);
	float ntotalevent(0);
	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float sigMETPhi(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float nJet(0);
  std::vector<int>   *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int>   *mcMomPID=0;
  std::vector<int>   *mcGMomPID=0;

	sigtree->SetBranchAddress("crosssection",&crosssection);
	sigtree->SetBranchAddress("ntotalevent", &ntotalevent);
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("HT",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	sigtree->SetBranchAddress("nJet",      &nJet);
  sigtree->SetBranchAddress("mcPID",     &mcPID);
  sigtree->SetBranchAddress("mcEta",     &mcEta);
  sigtree->SetBranchAddress("mcPhi",     &mcPhi);
  sigtree->SetBranchAddress("mcPt",      &mcPt);
  sigtree->SetBranchAddress("mcMomPID",  &mcMomPID);
  sigtree->SetBranchAddress("mcGMomPID", &mcGMomPID);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		double weight = 1;
		if(channelType == 1){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_DoubleEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_DoubleEG*1000*crosssection/ntotalevent;}
		else if(channelType == 2){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_MuonEG*1000*crosssection/ntotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_MuonEG*1000*crosssection/ntotalevent;}

		/** cut flow *****/
		if(phoEt < 35 || lepPt < lepPtCut)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;


		bool isFakePho(true);
		double mindR(0.3),deltaE(1);
		unsigned matchIndex(0);
		for(unsigned iMC(0); iMC < mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindR && dE < 0.5){mindR=dR; matchIndex=iMC;deltaE = dE;}
		}
		if(mindR < 0.1){
			if( (fabs((*mcPID)[matchIndex]) == 11 || fabs((*mcPID)[matchIndex]) == 22)  && (fabs((*mcMomPID)[matchIndex]) <= 6 || fabs((*mcMomPID)[matchIndex])==21 || fabs((*mcMomPID)[matchIndex])==11 || fabs((*mcMomPID)[matchIndex])== 999 || fabs((*mcMomPID)[matchIndex]) == 13 ||  fabs((*mcMomPID)[matchIndex]) == 15 || fabs((*mcMomPID)[matchIndex]) == 24 || fabs((*mcMomPID)[matchIndex]) == 23))isFakePho = false;
		}
		if(!isFakePho)continue;
		// only fake photon considered
		if(phoEt > MAXET)phoEt = MAXET;
		if(sigMET > MAXMET)sigMET = MAXMET;
		if(sigMT > MAXMT)sigMT = MAXMT;
		if(HT > MAXHT)HT = MAXHT;	
		// XSectional weight
		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta, weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta, weight);
		p_MET->Fill(sigMET, weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
		p_nJet->Fill(nJet, weight);
		
	}        


	// Fake Tree //
	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *pred_PhoEt = new TH1D("pred_PhoEt","#gamma E_{T}; E_{T} (GeV)",11,35,200);
	TH1D *pred_LepPt = new TH1D("pred_LepPt","pred_LepPt",11,35,200);
	TH1D *pred_MET = new TH1D("pred_MET","MET; MET (GeV);",10,0,200);
	TH1D *pred_Mt = new TH1D("pred_Mt","M_{T}; M_{T} (GeV);",10,0,200);
	TH1D *pred_HT = new TH1D("pred_HT","HT; HT (GeV);",10,0,400);
	TH1D *pred_PhoEta = new TH1D("pred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *pred_LepEta = new TH1D("pred_LepEta","pred_LepEta",60,-3,3);
	TH1D *pred_dPhiEleMET = new TH1D("pred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *pred_nJet = new TH1D("pred_nJet","pred_nJet",10,0,10);

	TH1D *DY_PhoEt = new TH1D("DY_PhoEt","#gamma E_{T}; E_{T} (GeV)",11,35,200);
	TH1D *DY_LepPt = new TH1D("DY_LepPt","DY_LepPt",11,35,200);
	TH1D *DY_MET = new TH1D("DY_MET","MET; MET (GeV);",10,0,200);
	TH1D *DY_Mt = new TH1D("DY_Mt","M_{T}; M_{T} (GeV);",10,0,200);
	TH1D *DY_HT = new TH1D("DY_HT","HT; HT (GeV);",10,0,400);
	TH1D *DY_PhoEta = new TH1D("DY_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *DY_LepEta = new TH1D("DY_LepEta","DY_LepEta",60,-3,3);
	TH1D *DY_dPhiEleMET = new TH1D("DY_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *DY_nJet = new TH1D("DY_nJet","DY_nJet",10,0,10);

//	TH1D *toy_PhoEt[NTOY];
//	TH1D *toy_LepPt[NTOY];
//	TH1D *toy_HT[NTOY];
//	TH1D *toy_MET[NTOY];
//	TH1D *toy_Mt[NTOY];
//	TH1D *toy_dPhiEleMET[NTOY];
//	for(unsigned ih(0); ih < NTOY; ih++){
//		histname.str("");
//		histname << "toy_PhoEt_ " << ih;
//		toy_PhoEt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
//		histname.str("");
//		histname << "toy_LepPt_" << ih;
//		toy_LepPt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgPtBins,bkgPtBins);
//		histname.str("");
//		histname << "toy_MET_ " << ih;
//		toy_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
//		histname.str("");
//		histname << "toy_Mt_ " << ih;
//		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
//		histname.str("");
//		histname << "toy_HT_ " << ih;
//		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
//		histname.str("");
//		histname << "toy_eledPhiEleMET_" << ih;
//		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
//	}
	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("jetTree");
	// jetTree for proxy events
	// proxy events weighted by fake rate, WJet has major contribution
	if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));

	float proxycrosssection(0);
	float proxyntotalevent(0);
	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxyphoPhi(0);
	float proxylepPt(0);
	float proxylepEta(0);
	float proxylepPhi(0);
	float proxysigMT(0);
	float proxysigMET(0);
	float proxysigMETPhi(0);
	float proxydPhiLepMET(0);
	int   proxynVertex(0);
	float proxydRPhoLep(0);
	float proxyHT(0);
	float proxynJet(0);
	
	proxytree->SetBranchAddress("crosssection",&proxycrosssection);
	proxytree->SetBranchAddress("ntotalevent", &proxyntotalevent);
	proxytree->SetBranchAddress("phoEt",     	 &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    	 &proxyphoEta);
	proxytree->SetBranchAddress("phoPhi",    	 &proxyphoPhi);
	proxytree->SetBranchAddress("lepPt",     	 &proxylepPt);
	proxytree->SetBranchAddress("lepEta",    	 &proxylepEta);
	proxytree->SetBranchAddress("lepPhi",    	 &proxylepPhi);
	proxytree->SetBranchAddress("sigMT",     	 &proxysigMT);
	proxytree->SetBranchAddress("sigMET",    	 &proxysigMET);
	proxytree->SetBranchAddress("sigMETPhi", 	 &proxysigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",	 &proxydPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   	 &proxynVertex);
	proxytree->SetBranchAddress("dRPhoLep",  	 &proxydRPhoLep);
	proxytree->SetBranchAddress("HT",        	 &proxyHT);
	proxytree->SetBranchAddress("nJet",      	 &proxynJet);

	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		double weight = 1;
		if(channelType == 1){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_DoubleEG*1000*proxycrosssection/proxyntotalevent;}
		else if(channelType == 2){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_MuonEG*1000*proxycrosssection/proxyntotalevent;}

		if(proxyphoEt > MAXET)proxyphoEt = MAXET;
		if(proxysigMET > MAXMET)proxysigMET = MAXMET;
		if(proxysigMT > MAXMT)proxysigMT = MAXMT;
		if(proxyHT > MAXHT)proxyHT = MAXHT;	
		
		/** cut flow *****/
		if(proxyphoEt < 35 || proxylepPt < lepPtCut)continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;
		double w_jet(1);
		w_jet = fitfunc_num->Eval(proxyphoEt)/fitfunc_den->Eval(proxyphoEt);
//		double jetfakeerror(0);
//		for(int ipt(0); ipt < 264; ipt++){
//			if(proxyphoEt >= ipt+35 && proxyphoEt < ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt]/fitfunc_den->Eval(proxyphoEt)/fitfunc_den->Eval(proxyphoEt) + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(proxyphoEt);
//		}
//		if(proxyphoEt >= 264)jetfakeerror = sqrt(jetfake_numerror[263]*jetfake_numerror[263]/fitfunc_den->Eval(300)/fitfunc_den->Eval(300) + jetfake_denerror[263]*jetfake_denerror[263]*w_jet*w_jet)/fitfunc_den->Eval(300); 
//		double sysJetFakePho = jetfakeerror/w_jet;
		w_jet = w_jet*weight;
		// XSec weight * fake rate weight

		pred_PhoEt->Fill(proxyphoEt,w_jet);
		pred_PhoEta->Fill(proxyphoEta, w_jet);
		pred_MET->Fill(proxysigMET, w_jet);
		pred_Mt->Fill(proxysigMT, w_jet);
		pred_HT->Fill(proxyHT, w_jet);
		pred_LepPt->Fill(proxylepPt, w_jet);
		pred_LepEta->Fill(proxylepEta, w_jet);
		pred_dPhiEleMET->Fill(fabs(proxydPhiLepMET), w_jet);
		pred_nJet->Fill(proxynJet, w_jet);

		DY_PhoEt->Fill(proxyphoEt,w_jet);
		DY_PhoEta->Fill(proxyphoEta, w_jet);
		DY_MET->Fill(proxysigMET, w_jet);
		DY_Mt->Fill(proxysigMT, w_jet);
		DY_HT->Fill(proxyHT, w_jet);
		DY_LepPt->Fill(proxylepPt, w_jet);
		DY_LepEta->Fill(proxylepEta, w_jet);
		DY_dPhiEleMET->Fill(fabs(proxydPhiLepMET), w_jet);
		DY_nJet->Fill(proxynJet, w_jet);

//		for(unsigned it(0); it < NTOY; it++){
//			double toy_jet = w_jet*(1+sysJetFakePho*randomweight_jet[ih]);
//			toy_PhoEt[it]->Fill(proxyphoEt,toy_jet);
//			toy_MET[it]->Fill(proxysigMET, toy_jet);
//			toy_Mt[it]->Fill(proxysigMT, toy_jet);
//			toy_HT[it]->Fill(proxyHT, toy_jet);
//			toy_LepPt[it]->Fill(proxylepPt, toy_jet);
//			toy_dPhiEleMET[it]->Fill(fabs(proxydPhiLepMET), toy_jet);
//		}
	}

	//************ Proxy Tree **********************//
	TChain *raretree = new TChain("jetTree");
	// rare contribution from DY
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));

	//if(channelType==2)raretree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_DY.root");

	float rarecrosssection(0);
	float rarentotalevent(0);
	float rarephoEt(0);
	float rarephoEta(0);
	float rarephoPhi(0);
	float rarelepPt(0);
	float rarelepEta(0);
	float rarelepPhi(0);
	float raresigMT(0);
	float raresigMET(0);
	float raresigMETPhi(0);
	float raredPhiLepMET(0);
	int   rarenVertex(0);
	float raredRPhoLep(0);
	float rareHT(0);
	float rarenJet(0);
	
	raretree->SetBranchAddress("crosssection",&rarecrosssection);
	raretree->SetBranchAddress("ntotalevent", &rarentotalevent);
	raretree->SetBranchAddress("phoEt",     	 &rarephoEt);
	raretree->SetBranchAddress("phoEta",    	 &rarephoEta);
	raretree->SetBranchAddress("phoPhi",    	 &rarephoPhi);
	raretree->SetBranchAddress("lepPt",     	 &rarelepPt);
	raretree->SetBranchAddress("lepEta",    	 &rarelepEta);
	raretree->SetBranchAddress("lepPhi",    	 &rarelepPhi);
	raretree->SetBranchAddress("sigMT",     	 &raresigMT);
	raretree->SetBranchAddress("sigMET",    	 &raresigMET);
	raretree->SetBranchAddress("sigMETPhi", 	 &raresigMETPhi);
	raretree->SetBranchAddress("dPhiLepMET",	 &raredPhiLepMET);
	raretree->SetBranchAddress("nVertex",   	 &rarenVertex);
	raretree->SetBranchAddress("dRPhoLep",  	 &raredRPhoLep);
	raretree->SetBranchAddress("HT",        	 &rareHT);
	raretree->SetBranchAddress("nJet",      	 &rarenJet);

	for (unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		double weight=1;
		if(channelType == 1){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_DoubleEG*1000*rarecrosssection/rarentotalevent;}
		else if(channelType == 2){
                        if(RunYear == 2016 and preVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2016 and preVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2017)                        weight = lumi_2017_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(RunYear == 2018)                        weight = lumi_2018_MuonEG*1000*rarecrosssection/rarentotalevent;}

		if(rarephoEt > MAXET)rarephoEt = MAXET;
		if(raresigMET > MAXMET)raresigMET = MAXMET;
		if(raresigMT > MAXMT)raresigMT = MAXMT;
		if(rareHT > MAXHT)rareHT = MAXHT;	
		

		/** cut flow *****/
		if(rarephoEt < 35 || rarelepPt < lepPtCut)continue;
		if(fabs(rarephoEta) > 1.4442 || fabs(rarelepEta) > 2.5)continue;
		double w_jet(1);
		w_jet = fitfunc_num->Eval(rarephoEt)/fitfunc_den->Eval(rarephoEt);

	//	double jetfakeerror(0);
	//	for(int ipt(0); ipt < 264; ipt++){
	//		if(rarephoEt >= ipt+35 && rarephoEt < ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt]/fitfunc_den->Eval(rarephoEt)/fitfunc_den->Eval(rarephoEt) + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(rarephoEt);
	//	}
	//	if(rarephoEt >= 264)jetfakeerror = sqrt(jetfake_numerror[263]*jetfake_numerror[263]/fitfunc_den->Eval(300)/fitfunc_den->Eval(300) + jetfake_denerror[263]*jetfake_denerror[263]*w_jet*w_jet)/fitfunc_den->Eval(300); 
	//	double sysJetFakePho = jetfakeerror/w_jet;
		w_jet = w_jet*weight;
		// XSec weight * fake rate weight


		pred_PhoEt->Fill(rarephoEt,w_jet);
		pred_PhoEta->Fill(rarephoEta, w_jet);
		pred_MET->Fill(raresigMET, w_jet);
		pred_Mt->Fill(raresigMT, w_jet);
		pred_HT->Fill(rareHT, w_jet);
		pred_LepPt->Fill(rarelepPt, w_jet);
		pred_LepEta->Fill(rarelepEta, w_jet);
		pred_dPhiEleMET->Fill(fabs(raredPhiLepMET), w_jet);
		pred_nJet->Fill(rarenJet, w_jet);

	//	for(unsigned it(0); it < NTOY; it++){
	//		double toy_jet = w_jet*(1+sysJetFakePho*randomweight_jet[ih]);
	//		toy_PhoEt[it]->Fill(rarephoEt,toy_jet);
	//		toy_MET[it]->Fill(raresigMET, toy_jet);
	//		toy_Mt[it]->Fill(raresigMT, toy_jet);
	//		toy_HT[it]->Fill(rareHT, toy_jet);
	//		toy_LepPt[it]->Fill(rarelepPt, toy_jet);
	//		toy_dPhiEleMET[it]->Fill(fabs(raredPhiLepMET), toy_jet);
	//	}
	}


//	for(int ibin(1); ibin < pred_PhoEt->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_PhoEt->GetBinContent(ibin),1.5*pred_PhoEt->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_PhoEt[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_PhoEt->GetBinError(ibin)*pred_PhoEt->GetBinError(ibin));
//		pred_PhoEt->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
//	for(int ibin(1); ibin < pred_LepPt->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_LepPt->GetBinContent(ibin),1.5*pred_LepPt->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_LepPt[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_LepPt->GetBinError(ibin)*pred_LepPt->GetBinError(ibin));
//		pred_LepPt->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
//	for(int ibin(1); ibin < pred_MET->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_MET->GetBinContent(ibin),1.5*pred_MET->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_MET[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_MET->GetBinError(ibin)*pred_MET->GetBinError(ibin));
//		pred_MET->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
//	for(int ibin(1); ibin < pred_Mt->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_Mt->GetBinContent(ibin),1.5*pred_Mt->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_Mt[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_Mt->GetBinError(ibin)*pred_Mt->GetBinError(ibin));
//		pred_Mt->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
//	for(int ibin(1); ibin < pred_dPhiEleMET->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_dPhiEleMET->GetBinContent(ibin),1.5*pred_dPhiEleMET->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_dPhiEleMET[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_dPhiEleMET->GetBinError(ibin)*pred_dPhiEleMET->GetBinError(ibin));
//		pred_dPhiEleMET->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
//	for(int ibin(1); ibin < pred_HT->GetSize(); ibin++){
//		TH1D *temphist = new TH1D("temphist","",500,0.5*pred_HT->GetBinContent(ibin),1.5*pred_HT->GetBinContent(ibin));
//		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_HT[it]->GetBinContent(ibin));
//		temphist->Fit("gaus");
//		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
//		double totalerror = sqrt(syserr*syserr + pred_HT->GetBinError(ibin)*pred_HT->GetBinError(ibin));
//		pred_HT->SetBinError(ibin, totalerror);
//		delete temphist;
//	}
		
	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);

	gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "",600,600);
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.35, 1, 1.0);
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_PhoEt->SetTitle("");
	p_PhoEt->GetXaxis()->SetRangeUser(35,200);
	p_PhoEt->Draw();
	p_PhoEt->SetLineColor(kBlack);
	p_PhoEt->SetMarkerStyle(20);
	pred_PhoEt->SetFillStyle(1001);
	pred_PhoEt->SetLineColor(kRed);
	pred_PhoEt->SetFillColor(kRed);
	DY_PhoEt->SetFillStyle(1001);
	DY_PhoEt->SetLineColor(kYellow);
	DY_PhoEt->SetFillColor(kYellow);
	pred_PhoEt->Draw("hist same");
	DY_PhoEt->Draw("hist same");
	  TLegend *leg =  new TLegend(0.5,0.55,0.9,0.8);
        leg->SetFillStyle(0);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetLegendFillColor(0);
	
	leg->AddEntry(p_PhoEt,"observed");
        leg->AddEntry(pred_PhoEt,"W+jets");
        leg->AddEntry(DY_PhoEt,"DY");
	leg->AddEntry(error_PhoEt, "Syst. Unc.");
	leg->Draw("same");
	
	p_PhoEt->Draw("same");
	for(int ibin(1); ibin < pred_PhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), pred_PhoEt->GetBinContent(ibin));
		float prederror = pred_PhoEt->GetBinError(ibin);
		error_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2, prederror/pred_PhoEt->GetBinContent(ibin)); 
	}
	error_PhoEt->SetFillColor(15);
  error_PhoEt->SetFillStyle(3345);
	error_PhoEt->Draw("E2 same");

        TLatex chantex;
        chantex.SetNDC();
        chantex.SetTextFont(42);
        chantex.SetTextSize(0.07);
        if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
        if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
        gPad->RedrawAxis();
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( pt_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( pt_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                  CMS_lumi( pt_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                  CMS_lumi( pt_pad1, 4, ichannel, 11 );

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.35);
	pt_pad2->SetTopMargin(0);
        pt_pad2->SetBottomMargin(0.3);

	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,200,1);
	TH1F *ratio=(TH1F*)p_PhoEt->Clone("transfer factor");
	ratio->SetMinimum(0);
	ratio->SetMaximum(2);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(pred_PhoEt);
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitleOffset(1.0);
	ratio->GetYaxis()->SetTitleOffset(0.5);

	ratio->GetYaxis()->SetTitleSize(0.08);
	ratio->GetXaxis()->SetTitleSize(0.08);
	ratio->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
        ratio->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(11);
	ratio->Draw();
	ratioerror_PhoEt->SetFillColor(15);
	ratioerror_PhoEt->SetFillStyle(3345);
	ratioerror_PhoEt->Draw("E2 same");
	ratio->Draw("same");
	flatratio->Draw("same");
	if(channelType==1) 	c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_PhotonEt_eg_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) 	c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_PhotonEt_mg_%d%s.pdf",RunYear,whichVFP.c_str()));



// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.35, 1, 1.0);
	met_pad1->SetBottomMargin(0);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_MET->GetYaxis()->SetRangeUser(1,1000000);
	p_MET->SetMinimum(0.001);
	p_MET->GetXaxis()->SetRangeUser(0,400);
	p_MET->SetLineColor(1);
	p_MET->SetMarkerStyle(20);
	p_MET->SetTitle("");

	p_MET->Draw("P");
	DY_MET->SetFillStyle(1001);
	DY_MET->SetLineColor(kYellow-4);
	DY_MET->SetFillColor(kYellow-4);
	pred_MET->SetFillStyle(1001);
	pred_MET->SetLineColor(kRed);
	pred_MET->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_MET->GetSize(); ibin++){
		float prederror = pred_MET->GetBinError(ibin);
		error_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), pred_MET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), 1); 
		ratioerror_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2, prederror/pred_MET->GetBinContent(ibin)); 
	}
	cout<<"direct simulation "<< p_MET->Integral()<<endl;
//      cout<<"DY " <<DY_MET->Integral()<<endl;
        cout<<"DY and W + jets "<<pred_MET->Integral()<<endl;
	cout<<"Ratio "<<p_MET->Integral()/pred_MET->Integral()<<endl;

	pred_MET->Draw("hist same");
	DY_MET->Draw("hist same");
  	error_MET->SetFillColor(15);
  	error_MET->SetFillStyle(3345);
	error_MET->Draw("E2 same");
	leg->Draw("same");
	p_MET->Draw("E same");

	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
        if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
        //TLine *line_met = new TLine(70,0,70,10000);
        TLine *line_met = new TLine(70,0.05,70,1000*p_MET->GetBinContent(p_MET->GetMaximumBin()));
        line_met->SetLineStyle(2);
        line_met->Draw("same");
        TLatex* latex = new TLatex();
        latex->SetTextSize(0.05);
        latex->DrawLatex(20, 60000,"control");
        latex->DrawLatex(20, 30000,"region");
        gPad->RedrawAxis();
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( met_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( met_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                  CMS_lumi( met_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                  CMS_lumi( met_pad1, 4, ichannel, 11 );

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.35);
	met_pad2->SetTopMargin(0);
        met_pad2->SetBottomMargin(0.3);

	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,400,1);
	TH1F *ratio_met=(TH1F*)p_MET->Clone("transfer factor");
	ratio_met->GetXaxis()->SetRangeUser(0,400);

        ratio_met->GetYaxis()->SetNdivisions(504);
        ratio_met->SetLineColor(kBlack);
        ratio_met->SetMarkerStyle(20);
        ratio_met->Divide(pred_MET);
        ratio_met->SetTitle("");
	ratio_met->GetXaxis()->SetTitleOffset(1.0);
        ratio_met->GetYaxis()->SetTitleOffset(0.5);

        ratio_met->GetYaxis()->SetTitleSize(0.08);
        ratio_met->GetXaxis()->SetTitleSize(0.08);

        ratio_met->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
        ratio_met->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_met->GetYaxis()->SetRangeUser(0,2.1);
	ratio_met->GetXaxis()->SetLabelFont(63);
	ratio_met->GetXaxis()->SetLabelSize(14);
	ratio_met->GetYaxis()->SetLabelFont(63);
	ratio_met->GetYaxis()->SetLabelSize(11);
	ratio_met->Draw();
	ratioerror_MET->SetFillColor(15);
	ratioerror_MET->SetFillStyle(3345);
	ratioerror_MET->Draw("E2 same");
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	if(channelType==1) 	c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_MET_eg_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) 	c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_MET_mg_%d%s.pdf",RunYear,whichVFP.c_str()));

// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.35, 1, 1.0);
	mt_pad1->SetBottomMargin(0);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_Mt->SetMinimum(0.001);
  p_Mt->SetMaximum(1000000);
	p_Mt->GetXaxis()->SetRangeUser(0,400);
	p_Mt->SetLineColor(1);
	p_Mt->SetMarkerStyle(20);
	p_Mt->SetTitle("");
	p_Mt->Draw("P");
	DY_Mt->SetFillStyle(1001);
	DY_Mt->SetLineColor(kYellow-4);
	DY_Mt->SetFillColor(kYellow-4);
	pred_Mt->SetFillStyle(1001);
	pred_Mt->SetLineColor(kRed);
	pred_Mt->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_Mt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), pred_Mt->GetBinContent(ibin));
		float prederror = pred_Mt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2, prederror/pred_Mt->GetBinContent(ibin)); 
	}
	pred_Mt->Draw("hist same");
	DY_Mt->Draw("hist same");
  error_Mt->SetFillColor(15);
  error_Mt->SetFillStyle(3345);
	error_Mt->Draw("E2 same");
	leg->Draw("same");
	p_Mt->Draw("E same");

	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
        if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
        gPad->RedrawAxis();
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( mt_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( mt_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                  CMS_lumi( mt_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                  CMS_lumi( mt_pad1, 4, ichannel, 11 );

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.35);
	mt_pad2->SetTopMargin(0);
        mt_pad2->SetBottomMargin(0.3);

	mt_pad2->Draw();
	mt_pad2->cd();
  TLine *flatratio_mt = new TLine(0,1,400,1);
	TH1F *ratio_mt=(TH1F*)p_Mt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,400);
	ratio_mt->GetYaxis()->SetRangeUser(0,2.1);
	ratio_mt->SetMinimum(0);
	ratio_mt->SetMaximum(2);
	ratio_mt->Divide(pred_Mt);
	ratio_mt->SetTitle("");
	ratio_mt->GetXaxis()->SetTitleOffset(1.0);
        ratio_mt->GetYaxis()->SetTitleOffset(0.5);

        ratio_mt->GetYaxis()->SetTitleSize(0.08);
        ratio_mt->GetXaxis()->SetTitleSize(0.08);

	ratio_mt->GetXaxis()->SetTitle("M_{T} (GeV)");
        ratio_mt->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(11);
	ratio_mt->Draw();
	ratioerror_Mt->SetFillColor(15);
	ratioerror_Mt->SetFillStyle(3345);
	ratioerror_Mt->Draw("E2 same");
	ratio_mt->Draw("same");
	flatratio_mt->Draw("same");
	if(channelType==1)	c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_MT_eg_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2)	c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_MT_mg_%d%s.pdf",RunYear,whichVFP.c_str()));

// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.35, 1, 1.0);
	HT_pad1->SetBottomMargin(0);
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_HT->GetXaxis()->SetRangeUser(0,400);
	p_HT->SetMinimum(1);
	p_HT->SetLineColor(1);
	p_HT->SetMarkerStyle(20);
	p_HT->SetTitle("");
	p_HT->Draw("P");
	DY_HT->SetFillStyle(1001);
	DY_HT->SetLineColor(kYellow-4);
	DY_HT->SetFillColor(kYellow-4);
	pred_HT->SetFillStyle(1001);
	pred_HT->SetLineColor(kRed);
	pred_HT->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_HT->GetSize(); ibin++){
		error_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), pred_HT->GetBinContent(ibin));
		float prederror = pred_HT->GetBinError(ibin);
		error_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), 1); 
		ratioerror_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2, prederror/pred_HT->GetBinContent(ibin)); 
	}
	pred_HT->Draw("hist same");
	DY_HT->Draw("hist same");
	error_HT->SetFillColor(15);
  error_HT->SetFillStyle(3345);
	error_HT->Draw("E2 same");
	leg->Draw("same");
	p_HT->Draw("E same");
	
	    if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
          if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
        gPad->RedrawAxis();
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( HT_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( HT_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                  CMS_lumi( HT_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                  CMS_lumi( HT_pad1, 4, ichannel, 11 );

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.35);
	HT_pad2->SetTopMargin(0);
        HT_pad2->SetBottomMargin(0.3);

	HT_pad2->Draw();
	HT_pad2->cd();
  TLine *flatratio_HT = new TLine(0,1,400,1);
	TH1F *ratio_HT=(TH1F*)p_HT->Clone("transfer factor");
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,400);
	ratio_HT->GetYaxis()->SetRangeUser(0,2.1);
	ratio_HT->SetMinimum(0);
	ratio_HT->SetMaximum(2);
	ratio_HT->Divide(pred_HT);
	ratio_HT->GetXaxis()->SetTitleOffset(1.0);
        ratio_HT->GetYaxis()->SetTitleOffset(0.5);

        ratio_HT->GetYaxis()->SetTitleSize(0.08);
        ratio_HT->GetXaxis()->SetTitleSize(0.08);

	ratio_HT->SetTitle("");
	ratio_HT->GetXaxis()->SetTitle("H_{T} (GeV)");

        ratio_HT->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_HT->GetXaxis()->SetLabelFont(63);
	ratio_HT->GetXaxis()->SetLabelSize(14);
	ratio_HT->GetYaxis()->SetLabelFont(63);
	ratio_HT->GetYaxis()->SetLabelSize(11);
	ratio_HT->Draw();
	ratioerror_HT->SetFillColor(15);
	ratioerror_HT->SetFillStyle(3345);
	ratioerror_HT->Draw("E2 same");
	flatratio_HT->Draw("same");
	if(channelType==1)	c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_HT_eg_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2)	c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Closure/closure_jetfakepho_HT_mg_%d%s.pdf",RunYear,whichVFP.c_str()));

}


