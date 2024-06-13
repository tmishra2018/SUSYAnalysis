// g++ `root-config --cflags` plot_dPhiTemplate.C -o plot_dPhiTemplate.exe `root-config --libs`
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
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/tdrstyle.C"

#define NTOY 1000
//int channelType = 1; // eg = 1; mg =2;
//int RunYear = 2016;
void plot_dPhiTemplate(int RunYear, int channelType){

	//setTDRStyle();   
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
  	gSystem->Load("../../lib/libAnaClasses.so");
        gROOT->SetBatch(kTRUE);
	// Signal Tree //
	//*********** hist o list **********************//
	TH1D *p_MET_WG = new TH1D("p_MET_WG","; p_{T}^{miss} (GeV);",40,0,200);
	TH1D *p_MET_ZG = new TH1D("p_MET_ZG","; p_{T}^{miss} (GeV);",40,0,200);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET",";|#Delta#phi(l,p_{T}^{miss})|;",32,0,3.2); 
	//************ Signal Tree **********************//
	TChain *mctree;
 	if(channelType == 1)mctree = new TChain("egTree","egTree");
  	else if(channelType == 2)mctree = new TChain("mgTree","mgTree");
	
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d.root",RunYear));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d.root",RunYear));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d.root",RunYear));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_%d.root",RunYear));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d.root",RunYear));
	
	float crosssection(0);
	float ntotalevent(0);
	int   mcType(0);
	float phoEt(0);
	float phoEta(0);
	float lepPt(0);
	float lepEta(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float llmass(0);
	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("mcType",    &mcType);
	mctree->SetBranchAddress("phoEt",     &phoEt);
	mctree->SetBranchAddress("phoEta",    &phoEta);
	mctree->SetBranchAddress("lepPt",     &lepPt);
	mctree->SetBranchAddress("lepEta",    &lepEta);
	mctree->SetBranchAddress("sigMET",    &sigMET);
	mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	mctree->SetBranchAddress("llmass",    &llmass);

	for (unsigned ievt(0); ievt<mctree->GetEntries(); ++ievt){//loop on entries
		mctree->GetEntry(ievt);

		double weight = crosssection/ntotalevent; 
		/** cut flow *****/
		if(phoEt < 40 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		if((mcType == 7) && llmass > 30)continue;
		if(mcType == 4 && llmass < 30)continue;
		if((mcType == 5 || mcType == 9) && llmass > 30)continue;

		if(mcType <= 3)p_MET_WG->Fill(sigMET, weight);
		else p_MET_ZG->Fill(sigMET, weight);

		if(sigMET < 40 || sigMET > 70)continue;
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
	}        


	// Fake Tree //
	//*********** hist o list **********************//
	std::ostringstream histname;
	TH1D *pred_dPhiEleMET = new TH1D("pred_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("fakeLepTree");
	if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d_April.root",RunYear));
	if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d_April.root",RunYear));

	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxylepPt(0);
	float proxylepEta(0);
	float proxysigMET(0);
	float proxydPhiLepMET(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
	
	proxytree->SetBranchAddress("phoEt",     	 &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    	 &proxyphoEta);
	proxytree->SetBranchAddress("lepPt",     	 &proxylepPt);
	proxytree->SetBranchAddress("lepEta",    	 &proxylepEta);
	proxytree->SetBranchAddress("sigMET",    	 &proxysigMET);
	proxytree->SetBranchAddress("dPhiLepMET",	 &proxydPhiLepMET);
  proxytree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  proxytree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);

	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);

		double weight = 1; 
		
		/** cut flow *****/
		if(proxyphoEt < 40 || proxylepPt < 25)continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;

		bool isProxy(false);
		if(channelType==1){if(fakeLepMiniIso < 0.4)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < 0.4))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;
	
		if(proxysigMET < 40 || proxysigMET > 70)continue;
		pred_dPhiEleMET->Fill(fabs(proxydPhiLepMET), weight);
	}
		
	p_MET_WG->Sumw2();
	p_MET_ZG->Sumw2();
	p_dPhiEleMET->Sumw2();
	pred_dPhiEleMET->Sumw2();
	p_MET_WG->Scale(1.0/p_MET_WG->Integral(1,100));
	p_MET_ZG->Scale(1.0/p_MET_ZG->Integral(1,100));
	p_dPhiEleMET->Scale(1.0/p_dPhiEleMET->Integral(1,100));
	pred_dPhiEleMET->Scale(1.0/pred_dPhiEleMET->Integral(1,100));

	gStyle->SetOptStat(0);
	TCanvas *c_dphi = new TCanvas("dPhi","dPhi",600,600);
	c_dphi->cd();
	p_dPhiEleMET->GetYaxis()->SetRangeUser(0, 4*p_dPhiEleMET->GetBinContent(1));
	p_dPhiEleMET->GetXaxis()->SetLabelFont(63);
	p_dPhiEleMET->GetXaxis()->SetLabelSize(14);
	p_dPhiEleMET->GetYaxis()->SetLabelFont(63);
	p_dPhiEleMET->GetYaxis()->SetLabelSize(14);
	p_dPhiEleMET->Draw();
	p_dPhiEleMET->SetLineColor(kBlue);
	p_dPhiEleMET->SetMarkerColor(kBlue);
	p_dPhiEleMET->SetMarkerStyle(20);
	pred_dPhiEleMET->SetLineColor(kRed);
	pred_dPhiEleMET->SetMarkerColor(kRed);
	pred_dPhiEleMET->Draw("P same");
	TLegend *leg =  new TLegend(0.5,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->AddEntry(p_dPhiEleMET,"W#gamma+Z#gamma template");
	if(channelType==1) leg->AddEntry(pred_dPhiEleMET,"jet-fake-electron template");
	if(channelType==2) leg->AddEntry(pred_dPhiEleMET,"jet-fake-muon template");
	leg->Draw("same");
	p_dPhiEleMET->Draw("P same");
	if(channelType==1) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_egChannel_%d.pdf",RunYear,RunYear));
	if(channelType==2) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_mgChannel_%d.pdf",RunYear,RunYear));
	if(channelType==1) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_egChannel_%d.png",RunYear,RunYear));
	if(channelType==2) c_dphi->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_mgChannel_%d.png",RunYear,RunYear));

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
	p_MET_ZG->GetXaxis()->SetLabelFont(63);
	p_MET_ZG->GetXaxis()->SetLabelSize(14);
	p_MET_ZG->GetYaxis()->SetLabelFont(63);
	p_MET_ZG->GetYaxis()->SetLabelSize(14);
	p_MET_ZG->SetLineColor(kBlue);
	p_MET_ZG->Draw("hist");
	p_MET_WG->SetLineColor(kRed);
	p_MET_WG->Draw("hist  same");
	TLegend *leg_met =  new TLegend(0.7,0.75,0.9,0.85);
	leg_met->SetFillStyle(0);
	leg_met->AddEntry(p_MET_ZG,"Z#gamma");
	leg_met->AddEntry(p_MET_WG,"W#gamma");
	leg_met->Draw("same");
	if(channelType==1) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_MET_egChannel_%d.pdf",RunYear,RunYear));
	if(channelType==2) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_MET_mgChannel_%d.pdf",RunYear,RunYear));
	if(channelType==1) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_MET_egChannel_%d.png",RunYear,RunYear));
	if(channelType==2) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d/dphiTemplate_MET_mgChannel_%d.png",RunYear,RunYear));
}

int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    plot_dPhiTemplate(atoi(argv[1]),atoi(argv[2]));
    return 0;
}
