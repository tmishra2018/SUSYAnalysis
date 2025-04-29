#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1F.h"
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

#include "../../include/tdrstyle.C"
#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

int
plot_JetPhoTemplate(){

	setTDRStyle();   
//************ Signal Tree **********************//
	TChain *mctree = new TChain("egTree");
	mctree->Add("/eos/uscms/store/user/tmishra/InputFilesDATA/2017/plot_hadron_GJet_2017.root");
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

	mctree->SetBranchAddress("phoEt",   &mc_phoEt);
	mctree->SetBranchAddress("phoEta",  &mc_phoEta);
	mctree->SetBranchAddress("phoPhi",  &mc_phoPhi);
	mctree->SetBranchAddress("phoSigma",&mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",&mc_phoChIso);
	mctree->SetBranchAddress("mcPID",   &mc_mcPID);
	mctree->SetBranchAddress("mcEta",   &mc_mcEta);
	mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
	mctree->SetBranchAddress("mcPt",    &mc_mcPt);
	mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);

	TH1F *p_true = new TH1F("p_true",";Iso_{h^{#pm}};Events",20,0,20);
	TH1F *p_fake = new TH1F("p_fake",";Iso_{h^{#pm}};Events",20,0,20);
	TH1F *p_truesigma = new TH1F("p_truesigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	TH1F *p_sigsigma = new TH1F("p_sigsigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	TH1F *p_sbsigma = new TH1F("p_sbsigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
	
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

		if(!isFake){
			p_truesigma->Fill(mc_phoSigma);
			if(mc_phoSigma <= 0.0106)p_sigsigma->Fill(mc_phoSigma);
			else if(mc_phoSigma > 0.0106 && mc_phoSigma < 0.016)p_sbsigma->Fill(mc_phoSigma);
		}

		if(mc_phoSigma > 0.0106)continue;
		if(!isFake)p_true->Fill(mc_phoChIso);
		if(isFake)p_fake->Fill(mc_phoChIso);
	}

	TCanvas *can_iso = new TCanvas("can_iso","",600,600);
	can_iso->cd();
	gPad->SetLogy();
	p_true->Scale(1.0/p_true->GetEntries());
	p_fake->Scale(1.0/p_true->GetEntries());	
	p_true->GetXaxis()->SetTitle("Iso_{h^{#pm}} (GeV)");
	p_true->GetXaxis()->SetTitleOffset(1);
	p_true->GetXaxis()->SetTitleSize(20);
	p_true->SetLineWidth(2);
	p_fake->SetLineWidth(2);
	p_true->Draw("hist");
	p_fake->SetLineColor(kRed);
	p_fake->Draw("hist same");
	TLegend *leg= new TLegend(0.5,0.7,0.85,0.85);
	leg->AddEntry(p_true, "true photons");
	leg->AddEntry(p_fake, "jet->#gamma fakes");
	leg->Draw("hist same");
	can_iso->SaveAs("/eos/uscms/store/user/tmishra/jetfakepho/Plots2017/PLOT_JetPhoTemplate_2017.pdf");

	TLatex* latex = new TLatex();
	TCanvas *can_sigma = new TCanvas("can_sigma","",600,600);
	can_sigma->cd();
	gPad->SetLogy();	
	p_truesigma->GetXaxis()->SetTitleOffset(1);
	p_truesigma->GetXaxis()->SetTitleSize(20);
	p_truesigma->Scale(1/p_truesigma->GetEntries());
	p_sigsigma->Scale(1/p_truesigma->GetEntries());
	p_sbsigma->Scale(1/p_truesigma->GetEntries());
	p_truesigma->Draw("hist");
	p_sigsigma->SetFillColor(kBlue);
	p_sbsigma->SetFillColor(kOrange);
	p_sigsigma->Draw("hist same");
	p_sbsigma->Draw("hist same");
	latex->DrawLatex(0.006,0.000005,"p_{signal}");
	latex->DrawLatex(0.011,0.000005,"p_{sideband}");
	can_sigma->SaveAs("/eos/uscms/store/user/tmishra/jetfakepho/Plots2017/PLOT_Iteration_2017.pdf");

	return 1;
}
