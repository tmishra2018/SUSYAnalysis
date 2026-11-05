#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TF3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFitResult.h"

#include "../../../include/analysis_commoncode.h"

int
plotChIsoMCFSR(){

	setTDRStyle();   
	TChain *datatree = new TChain("ZeeTree");
	//datatree->Add("/uscms_data/d3/mengleis/plotZ_DoubleMu_2016Rereco_eg.root");
	datatree->Add("/uscms_data/d3/mengleis/plotZ_DoubleMu_eg.root");
	TChain *mctree = new TChain("egTree");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_hadron_GJet.root");
	

	TH1D *mc_sig = new TH1D("mc_sig","photon I_{Ch}; charged isolation (GeV); ",20,0,20);
	TH1D *data_sig = new TH1D("data_sig","photon I_{Ch}; charged isolation (GeV); ",20,0,20);
	TH1D *mc_sigma = new TH1D("mc_sigma","Sieie; Sieie; ",20,0.005,0.015);
	TH1D *mc_fullsigma = new TH1D("mc_fullsigma","Sieie; Sieie; ",30,0.005,0.02);
	TH1D *data_sigma = new TH1D("data_sigma","Sieie; Sieie; ",20,0.005,0.015);
	
//************ Signal Tree **********************//
	float mc_phoEt(0);
	float mc_phoEta(0); 
	float mc_phoPhi(0); 
	float mc_phoSigma(0);
	float mc_phoChIso(0);
	std::vector<int> *mc_mcPID=0;
	std::vector<float> *mc_mcEta=0;
	std::vector<float> *mc_mcPhi=0;
	std::vector<float> *mc_mcPt=0;
	std::vector<int> *mc_mcMomPID=0; 

	mctree->SetBranchAddress("phoEt",     &mc_phoEt);
	mctree->SetBranchAddress("phoEta",    &mc_phoEta);
	mctree->SetBranchAddress("phoPhi",    &mc_phoPhi);
	mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);
	mctree->SetBranchAddress("mcPID",   &mc_mcPID);
	mctree->SetBranchAddress("mcEta",   &mc_mcEta);
	mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
	mctree->SetBranchAddress("mcPt",    &mc_mcPt);
	mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);

	float phoEt(0);
	float phoEta(0);
	float sigMET(0);
	float phoSigma(0);
	float phoChIso(0);
	float invmass(0);
	datatree->SetBranchAddress("phoEt",     &phoEt);
	datatree->SetBranchAddress("phoEta",    &phoEta);
	datatree->SetBranchAddress("sigMET",    &sigMET);
	datatree->SetBranchAddress("phoSigma",  &phoSigma);
	datatree->SetBranchAddress("phoChIso",  &phoChIso);
	datatree->SetBranchAddress("invmass",   &invmass);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);

		if( mc_phoEt < 35)continue;
		if( fabs(mc_phoEta) > 1.4442)continue;

		bool isFake(true);
		unsigned mcIndex(0);
		unsigned mcPhoIndex(0);
		float mindR(0.7);
		float phodR(0.7);
		bool hasMatch(false);
		for(unsigned ii(0); ii < mc_mcPID->size(); ii++){
			float dR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
			float dE = fabs(mc_phoEt - (*mc_mcPt)[ii])/mc_phoEt;
			if((*mc_mcPID)[ii] == 22){
				phodR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
				mcPhoIndex = ii;
			}
			if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;}
		}
		if(phodR < 0.1)mcIndex = mcPhoIndex;
		if(hasMatch)isFake = isHad(fabs((*mc_mcPID)[mcIndex]), fabs((*mc_mcMomPID)[mcIndex]));

		bool isTrueTemplate = (mc_phoSigma <= 0.0106 && !isFake);
		if(isTrueTemplate){
			mc_sig->Fill(mc_phoChIso);
		}
		if( mc_phoChIso < 1.694)mc_sigma->Fill(mc_phoSigma);
	}


	for(unsigned ievt(0); ievt < datatree->GetEntries(); ievt++){
		datatree->GetEntry(ievt);

		if( phoEt < 35)continue;
		if( fabs(phoEta) > 1.4442)continue;
		if(invmass < 80 || invmass > 100)continue;

		if(phoSigma < 0.0106){
			data_sig->Fill(phoChIso);
		}
		if( phoChIso < 1.694)data_sigma->Fill(phoSigma);
	}

	TCanvas *can = new TCanvas("HT", "HT",600,600);
	can->cd();
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
	setTopPad(pad1); 
	pad1->SetBottomMargin(0);
	pad1->Draw();  
	pad1->cd();  
	gPad->SetLogy();
	data_sig->SetLineColor(kRed);
	data_sig->SetMarkerColor(kRed);
  data_sig->Scale(mc_sig->GetEntries()/data_sig->GetEntries());
  data_sig->Draw();
  mc_sig->Draw("same");
	TLegend *leg= new TLegend(0.5,0.6,0.9,0.9);
	leg->AddEntry(data_sig, "Z->mu+mu+g");
	leg->AddEntry(mc_sig, "GJets MC");
	leg->Draw("same");
 	gPad->RedrawAxis();
	
	can->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	TH1D *ratio = (TH1D*)data_sig->Clone("ratio");
	ratio->SetMarkerColor(kBlack);
	ratio->GetYaxis()->SetTitle("data/MC");
	ratio->Divide(mc_sig);
  TLine *flatratio = new TLine(0,1,400,1);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->GetYaxis()->SetRangeUser(0,5);
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->Draw("EP");
	flatratio->Draw("same");
	can->SaveAs("MC_Data_ChIso.pdf");

	TCanvas *cansigma = new TCanvas("Sigma", "Sigma",600,600);
	cansigma->cd();
	TPad *padsigma1 = new TPad("padsigma1", "padsigma1", 0, 0.35, 1, 1.0);
	setTopPad(padsigma1); 
	padsigma1->SetBottomMargin(0);
	padsigma1->Draw();  
	padsigma1->cd();  
	gPad->SetLogy();
	data_sigma->SetLineColor(kRed);
	data_sigma->SetMarkerColor(kRed);
	data_sigma->Sumw2();
  data_sigma->Scale(mc_sigma->GetEntries()/data_sigma->GetEntries());
  data_sigma->Draw("EP same");
  mc_sigma->Draw("same");
	//leg->Draw("same");
 	gPad->RedrawAxis();
	
	cansigma->cd();
	TPad *padsigma2 = new TPad("padsigma2", "padsigma2", 0, 0, 1, 0.35);
	padsigma2->SetBottomMargin(0.3);
	padsigma2->Draw();
	padsigma2->cd();
	TH1D *ratiosigma = (TH1D*)data_sigma->Clone("ratiosigma");
	ratiosigma->SetMarkerColor(kBlack);
	ratiosigma->GetYaxis()->SetTitle("data/MC");
	ratiosigma->Divide(mc_sigma);
  TLine *flatratiosigma = new TLine(0,1,400,1);
	ratiosigma->SetMarkerStyle(20);
	ratiosigma->SetLineColor(kBlack);
	ratiosigma->GetYaxis()->SetRangeUser(0,5);
	ratiosigma->GetYaxis()->SetNdivisions(504);
	ratiosigma->Draw("EP");
	flatratiosigma->Draw("same");
	cansigma->SaveAs("MC_Data_Sigma.pdf");

	TCanvas *can_mc = new TCanvas("can_mc","",600,600);
	can_mc->cd();
	gPad->SetLogy();
	mc_sigma->Draw();

	return 1;
}
