// run using root -l "closure_jetfakelep.C(1)"
// 	Scale factor due to pT dependence of miniIsolation leads to disagreement in fake electron pT shape. 
	// this scale factor is to correct the electron proxy sample, as described in section 5.3.3 AN
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
#include "TGraphErrors.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/tdrstyle.C"

#define NTOY 1000
int RunYear = 2018;
bool ISpreVFP = false;
void closure_jetfakelep(int ichannel){

	setTDRStyle();   
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
  	gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;
	gROOT->SetBatch(1);
	// Signal Tree //
	//*********** hist o list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt","; E_{T} (GeV);",20,0,200);
	TH1D *p_LepPt = new TH1D("p_LepPt","LepPt; p_{T} (GeV);",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET",";p_{T}^{miss} (GeV);",20,0,100);
	TH1D *p_Mt = new TH1D("p_Mt",";M_{T}(l,p_{T}^{miss}) (GeV);",40,0,200);
	TH1D *p_HT = new TH1D("p_HT",";HT (GeV);",40,0,400); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET",";#Delta#phi(l,p_{T}^{miss});",32,0,3.2); 

	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",9,0,9);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	//************ Signal Tree **********************//
	
	std::string whichVFP;
	if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	// from QCD simulation
	// Inputs from ./QCD/analysis_egMC.C
	TChain *sigtree = new TChain("signalTree");
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/fakeLep/fakelep_egsignal_QCD_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/fakeLep/fakelep_mgsignal_QCD_%d%s.root",RunYear,whichVFP.c_str())); // looks correct from stat
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
		/** cut flow *****/
		if(phoEt < 40 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

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
	//*********** hist o list **********************//
	std::ostringstream histname;
	TH1D *pred_PhoEt = new TH1D("pred_PhoEt","#gamma E_{T}; E_{T} (GeV)",20,0,200);
	TH1D *pred_LepPt = new TH1D("pred_LepPt","pred_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *pred_MET = new TH1D("pred_MET","MET; MET (GeV);",20,0,100);
	TH1D *pred_Mt = new TH1D("pred_Mt","M_{T}; M_{T} (GeV);",40,0,200);
	TH1D *pred_HT = new TH1D("pred_HT","HT; HT (GeV);",40,0,400); 
	TH1D *pred_PhoEta = new TH1D("pred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *pred_LepEta = new TH1D("pred_LepEta","pred_LepEta",60,-3,3);
	TH1D *pred_dPhiEleMET = new TH1D("pred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *pred_nJet = new TH1D("pred_nJet","pred_nJet",10,0,10);

	//************ Proxy Tree **********************//
	// proxy events enriched in fake leptons from data
	TChain *proxytree = new TChain("fakeLepTree");
        if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
        if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

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
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
	
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
  	proxytree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  	proxytree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);
	proxytree->SetBranchAddress("nVertex",   	 &proxynVertex);
	proxytree->SetBranchAddress("dRPhoLep",  	 &proxydRPhoLep);
	proxytree->SetBranchAddress("HT",        	 &proxyHT);
	proxytree->SetBranchAddress("nJet",      	 &proxynJet);

	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);

		double weight = 1; 
		
		/** cut flow *****/
		if(proxyphoEt < 40 || proxylepPt < 25)continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;

		bool isProxy(false);
		// for ele proxy, miniIso<0.4 and for muon proxy 0.2<miniIso<0.4
		if(channelType==1){if(fakeLepMiniIso < 0.4)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < 0.4))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		// only take lepton proxies
		if(!isProxy)continue;

		pred_PhoEt->Fill(proxyphoEt,weight);
		pred_PhoEta->Fill(proxyphoEta, weight);
		pred_MET->Fill(proxysigMET, weight);
		pred_Mt->Fill(proxysigMT, weight);
		pred_HT->Fill(proxyHT, weight);
		pred_LepPt->Fill(proxylepPt, weight);
		pred_LepEta->Fill(proxylepEta, weight);
		pred_dPhiEleMET->Fill(fabs(proxydPhiLepMET), weight);
		pred_nJet->Fill(proxynJet, weight);
	}
		
	//MC Fake Tree //
	//*********** hist o list **********************//
	TH1D *mcpred_PhoEt = new TH1D("mcpred_PhoEt","#gamma E_{T}; E_{T} (GeV)",20,0,200);
	TH1D *mcpred_LepPt = new TH1D("mcpred_LepPt","mcpred_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *mcpred_MET = new TH1D("mcpred_MET","MET; MET (GeV);",20,0,100);
	TH1D *mcpred_Mt = new TH1D("mcpred_Mt","M_{T}; M_{T} (GeV);",40,0,200);
	TH1D *mcpred_HT = new TH1D("mcpred_HT","HT; HT (GeV);",40,0,400); 
	TH1D *mcpred_PhoEta = new TH1D("mcpred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *mcpred_LepEta = new TH1D("mcpred_LepEta","mcpred_LepEta",60,-3,3);
	TH1D *mcpred_dPhiEleMET = new TH1D("mcpred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *mcpred_nJet = new TH1D("mcpred_nJet","mcpred_nJet",10,0,10);

	//************ Proxy Tree **********************//
	TChain *mcproxytree = new TChain("fakeLepTree");
	// proxy events enriched in fake leptons from MC
	if(channelType==1)	mcproxytree->Add(Form("/eos/uscms/store/user/tmishra/fakeLep/fakelep_egsignal_QCD_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)	mcproxytree->Add(Form("/eos/uscms/store/user/tmishra/fakeLep/fakelep_mgsignal_QCD_%d%s.root",RunYear,whichVFP.c_str()));

	float mcproxyphoEt(0);
	float mcproxyphoEta(0);
	float mcproxyphoPhi(0);
	float mcproxylepPt(0);
	float mcproxylepEta(0);
	float mcproxylepPhi(0);
	float mcproxysigMT(0);
	float mcproxysigMET(0);
	float mcproxysigMETPhi(0);
	float mcproxydPhiLepMET(0);
	int   mcproxynVertex(0);
	float mcproxydRPhoLep(0);
	float mcproxyHT(0);
	float mcproxynJet(0);
	float mcfakeLepMiniIso(0);
	
	mcproxytree->SetBranchAddress("phoEt",     	 &mcproxyphoEt);
	mcproxytree->SetBranchAddress("phoEta",    	 &mcproxyphoEta);
	mcproxytree->SetBranchAddress("phoPhi",    	 &mcproxyphoPhi);
	mcproxytree->SetBranchAddress("lepPt",     	 &mcproxylepPt);
	mcproxytree->SetBranchAddress("lepEta",    	 &mcproxylepEta);
	mcproxytree->SetBranchAddress("lepPhi",    	 &mcproxylepPhi);
	mcproxytree->SetBranchAddress("sigMT",     	 &mcproxysigMT);
	mcproxytree->SetBranchAddress("sigMET",    	 &mcproxysigMET);
	mcproxytree->SetBranchAddress("sigMETPhi", 	 &mcproxysigMETPhi);
	mcproxytree->SetBranchAddress("dPhiLepMET",	 &mcproxydPhiLepMET);
  	mcproxytree->SetBranchAddress("fakeLepMiniIso", &mcfakeLepMiniIso);
	mcproxytree->SetBranchAddress("nVertex",   	 &mcproxynVertex);
	mcproxytree->SetBranchAddress("dRPhoLep",  	 &mcproxydRPhoLep);
	mcproxytree->SetBranchAddress("HT",        	 &mcproxyHT);
	mcproxytree->SetBranchAddress("nJet",      	 &mcproxynJet);

	for (unsigned ievt(0); ievt<mcproxytree->GetEntries(); ++ievt){//loop on entries
		mcproxytree->GetEntry(ievt);

		double weight = 1; 
		
		/** cut flow *****/
		if(mcproxyphoEt < 40 || mcproxylepPt < 25)continue;
		if(fabs(mcproxyphoEta) > 1.4442 || fabs(mcproxylepEta) > 2.5)continue;

		bool isProxy(false);
		// for ele proxy, miniIso<0.4 and for muon proxy 0.2<miniIso<0.4
		if(channelType==1){if(mcfakeLepMiniIso < 0.4)isProxy=true;}
		else if(channelType==2){if((mcfakeLepMiniIso > 0.2 && mcfakeLepMiniIso < 0.4))isProxy=true;}
		// only take lepton proxies
		if(!isProxy)continue;

		mcpred_PhoEt->Fill(mcproxyphoEt,weight);
		mcpred_PhoEta->Fill(mcproxyphoEta, weight);
		mcpred_MET->Fill(mcproxysigMET, weight);
		mcpred_Mt->Fill(mcproxysigMT, weight);
		mcpred_HT->Fill(mcproxyHT, weight);
		mcpred_LepPt->Fill(mcproxylepPt, weight);
		mcpred_LepEta->Fill(mcproxylepEta, weight);
		mcpred_dPhiEleMET->Fill(fabs(mcproxydPhiLepMET), weight);
		mcpred_nJet->Fill(mcproxynJet, weight);
	}


	p_PhoEt->Sumw2();
	p_MET->Sumw2();
	p_Mt->Sumw2();
	p_HT->Sumw2();
	p_dPhiEleMET->Sumw2();
	// unit normalised
	p_PhoEt->Scale(1.0/p_PhoEt->GetEntries());
	p_MET->Scale(1.0/p_MET->GetEntries());
	p_Mt->Scale(1.0/p_Mt->GetEntries());
	p_HT->Scale(1.0/p_HT->GetEntries());
	p_dPhiEleMET->Scale(1.0/p_dPhiEleMET->GetEntries());

	pred_PhoEt->Scale(1.0/pred_PhoEt->GetEntries());
	pred_MET->Scale(1.0/pred_MET->GetEntries());
	pred_Mt->Scale(1.0/pred_Mt->GetEntries());
	pred_HT->Scale(1.0/pred_HT->GetEntries());
	pred_dPhiEleMET->Scale(1.0/pred_dPhiEleMET->GetEntries());

	mcpred_PhoEt->Scale(1.0/mcpred_PhoEt->GetEntries());
	mcpred_MET->Scale(1.0/mcpred_MET->GetEntries());
	mcpred_Mt->Scale(1.0/mcpred_Mt->GetEntries());
	mcpred_HT->Scale(1.0/mcpred_HT->GetEntries());
	mcpred_dPhiEleMET->Scale(1.0/mcpred_dPhiEleMET->GetEntries());

	gStyle->SetOptStat(0);
	TCanvas *c_dphi = new TCanvas("dPhi","dPhi",600,600);
	c_dphi->cd();
	p_dPhiEleMET->GetYaxis()->SetRangeUser(0, 1.5*p_dPhiEleMET->GetBinContent(1));
	// directly from simulation
	p_dPhiEleMET->Draw();
	p_dPhiEleMET->SetLineColor(kBlack);
	p_dPhiEleMET->SetMarkerStyle(20);
	// proxy template from data
	pred_dPhiEleMET->SetLineColor(kRed);
	pred_dPhiEleMET->Draw("hist  same");
	// proxy template from MC
	//mcpred_dPhiEleMET->SetLineColor(kCyan);
	//mcpred_dPhiEleMET->Draw("hist  same");
	TLegend *leg =  new TLegend(0.5,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	pred_dPhiEleMET->SetMarkerColor(kRed);
	leg->AddEntry(p_dPhiEleMET,"Simulation");
	leg->AddEntry(pred_dPhiEleMET,"data proxy template");
	leg->Draw("same");
	p_dPhiEleMET->Draw("P same");
	if(channelType==1) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_electron_dPhi_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_muon_dPhi_%d%s.pdf",RunYear,whichVFP.c_str()));

	// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->SetBottomMargin(0.12);
	c_met->SetLeftMargin(0.12);
	c_met->cd();
	gPad->SetLogy();
	p_MET->SetMinimum(0.0001);
	p_MET->SetMaximum(1);
	p_MET->GetXaxis()->SetTitleOffset(0.9);
	p_MET->GetXaxis()->SetRangeUser(0,100);
	p_MET->SetLineColor(1);
	p_MET->SetMarkerStyle(20);
	p_MET->Draw();
	pred_MET->SetLineColor(kRed);
	pred_MET->Draw("hist  same");
	//	mcpred_MET->SetLineColor(kCyan);
	//	mcpred_MET->Draw("hist  same");
	leg->Draw("same");
	p_MET->Draw("E same");

	if(channelType==1) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_electron_met_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_muon_met_%d%s.pdf",RunYear,whichVFP.c_str()));

	// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	c_mt->SetBottomMargin(0.12);
	c_mt->SetLeftMargin(0.12);
	gPad->SetLogy();
	p_Mt->SetMinimum(0.00001);
  	//p_Mt->SetMaximum(1000000);
	p_Mt->GetXaxis()->SetRangeUser(0,200);
	p_Mt->GetXaxis()->SetTitleOffset(0.9);
	p_Mt->SetLineColor(1);
	p_Mt->SetMarkerStyle(20);
	p_Mt->Draw("P");
	pred_Mt->SetLineColor(kRed);
	pred_Mt->Draw("hist  same");
	leg->Draw("same");
	p_Mt->Draw("E same");

	if(channelType==1) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_electron_mt_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_muon_mt_%d%s.pdf",RunYear,whichVFP.c_str()));

	// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	c_HT->cd();
	c_HT->SetBottomMargin(0.12);
        c_HT->SetLeftMargin(0.12);
	gPad->SetLogy();
	//p_HT->GetYaxis()->SetRangeUser(1,10*p_HT->GetBinContent(1));
	p_HT->GetXaxis()->SetRangeUser(0,400);
	p_HT->SetMinimum(1);
	p_HT->SetLineColor(1);
	p_HT->SetMarkerStyle(20);
	p_HT->Draw("P");
	pred_HT->SetLineColor(kRed);
	pred_HT->Draw("hist  same");
	leg->Draw("same");
	p_HT->Draw("E same");

	if(channelType==1) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_electron_ht_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_muon_ht_%d%s.pdf",RunYear,whichVFP.c_str()));

	// ******** LepPt ************************//
	p_LepPt->Sumw2();
	pred_LepPt->Sumw2();
	mcpred_LepPt->Sumw2();
	p_LepPt->Scale(1.0/p_LepPt->GetEntries());
	pred_LepPt->Scale(1.0/pred_LepPt->GetEntries());
	mcpred_LepPt->Scale(1.0/mcpred_LepPt->GetEntries());
	mcpred_LepPt->Scale(pred_LepPt->GetBinContent(1)/mcpred_LepPt->GetBinContent(1));
	gStyle->SetOptStat(0);
	TCanvas *c_LepPt = new TCanvas("LepPt", "LepPt",600,600);
	c_LepPt->cd();
	c_LepPt->SetBottomMargin(0.12);
        c_LepPt->SetLeftMargin(0.12);
	TPad *LepPt_pad1 = new TPad("LepPt_pad1", "LepPt_pad1", 0, 0.3, 1, 1.0);
	LepPt_pad1->SetBottomMargin(0);
	LepPt_pad1->Draw();  
	LepPt_pad1->cd();  
	gPad->SetLogy();
	//p_LepPt->GetYaxis()->SetRangeUser(1,10*p_LepPt->GetBinContent(1));
	p_LepPt->GetXaxis()->SetRangeUser(25,800);
	p_LepPt->SetMinimum(1);
	p_LepPt->SetLineColor(1);
	p_LepPt->SetMarkerStyle(20);
	p_LepPt->Draw("P");
	pred_LepPt->SetLineColor(kRed);
	pred_LepPt->SetLineWidth(2);
	pred_LepPt->Draw("hist same");
	mcpred_LepPt->SetLineColor(kCyan);
	mcpred_LepPt->SetLineWidth(2);
	//mcpred_LepPt->Draw("hist same");
	TLegend *leg_pt =  new TLegend(0.5,0.7,0.85,0.85);
	leg_pt->SetFillStyle(0);
	leg_pt->AddEntry(p_LepPt,"Simulation");
	leg_pt->AddEntry(pred_LepPt,"data proxy template");
	//leg_pt->AddEntry(mcpred_LepPt,"MC proxy template");
	leg_pt->Draw("same");
	p_LepPt->Draw("E same");

	c_LepPt->cd();
	TPad *LepPt_pad2 = new TPad("LepPt_pad2", "LepPt_pad2", 0, 0.05, 1, 0.3);
	LepPt_pad2->SetBottomMargin(0.4);
	LepPt_pad2->SetTopMargin(0);
	LepPt_pad2->Draw();
	LepPt_pad2->cd();
  	TLine *flatratio_LepPt = new TLine(25,1,800,1);
	// signal trees from QCD simulation divided by fakeLepTree from data
	TH1D *ratio_LepPt=(TH1D*)p_LepPt->Clone("transfer_factor");
	ratio_LepPt->SetMarkerStyle(20);
	ratio_LepPt->SetLineColor(kBlack);
	ratio_LepPt->GetXaxis()->SetRangeUser(25,800);
	ratio_LepPt->GetXaxis()->SetTitle("p_{T} (GeV)");
	ratio_LepPt->GetYaxis()->SetRangeUser(0,4);
	ratio_LepPt->SetMinimum(0);
	ratio_LepPt->SetMaximum(4.2);
	ratio_LepPt->GetYaxis()->SetLabelSize(11);
	ratio_LepPt->Divide(pred_LepPt);
	ratio_LepPt->SetTitle("");
	ratio_LepPt->GetYaxis()->SetTitleOffset(2);
	ratio_LepPt->GetYaxis()->SetTitleSize(12);
	ratio_LepPt->GetYaxis()->SetTitle("#frac{Simulation}{data proxy template}   ");
	ratio_LepPt->Draw();
	flatratio_LepPt->Draw("same");
	if(channelType==1) c_LepPt->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_electron_LepPt_%d%s.pdf",RunYear,whichVFP.c_str()));
	if(channelType==2) c_LepPt->SaveAs(Form("/eos/uscms/store/user/tmishra/fakeLep/faketemp_muon_LepPt_%d%s.pdf",RunYear,whichVFP.c_str()));
	TFile *outputfile;
	if(channelType==1)  outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/fakeLep/qcd_eg_scale_%d%s.root",RunYear,whichVFP.c_str()),"RECREATE");
	if(channelType==2)  outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/fakeLep/qcd_mg_scale_%d%s.root",RunYear,whichVFP.c_str()),"RECREATE");
	// this scale factor is to correct the electron proxy sample, as described in section 4.3.3 AN
	outputfile->cd();
	ratio_LepPt->Write();
	outputfile->Write();
	outputfile->Close(); 
}
