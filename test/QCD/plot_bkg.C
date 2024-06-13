// inputs from analysis_qcdBkg.C  and analysis_sig.C
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1F.h"
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
void plot_bkg(){//main  

  gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");

	TFile *file_sig = TFile::Open("signalTree_sig.root");
	TFile *file_qcd = TFile::Open("signalTree_qcd.root");

	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_dPhiEleMET");
	//TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	//TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepEta");
	p_allLepPt->SetMarkerColor(1);

	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_dPhiEleMET");
	//TH1F *p_qcdLepPt_2 = (TH1F*)file_qcd->Get("p_dPhiEleMET_loose");
	//TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	//TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepEta");

// ******** LepPt ************************//
	gStyle->SetOptStat(0);
	gROOT->SetBatch(1);
	TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
	c_leppt->cd();
	//gPad->SetLogy();
	p_allLepPt->SetMinimum(0.001);
	p_allLepPt->SetLineColor(1);
	p_allLepPt->SetMarkerStyle(20);
	p_allLepPt->Draw("P");

	p_qcdLepPt->Scale(p_allLepPt->GetEntries()/p_qcdLepPt->GetEntries());
	p_qcdLepPt->SetFillStyle(1001);
	p_qcdLepPt->SetLineColor(kBlue);
	p_qcdLepPt->SetFillColor(kBlue);
	p_qcdLepPt->Draw("same");

//	p_qcdLepPt_2->Scale(p_allLepPt->GetEntries()/p_qcdLepPt_2->GetEntries());
//	p_qcdLepPt_2->SetLineColor(kRed);	
//	p_qcdLepPt_2->Draw("hist same");
	p_allLepPt->Draw("P same");
	c_leppt->SaveAs("dPhiEleMET_2016.pdf");

	TFile *outputfile = TFile::Open("qcd_eg_scale.root","RECREATE");
	outputfile->cd();
	TH1F *p_scale = (TH1F*)p_allLepPt->Clone("p_scale");
	p_scale->Divide(p_qcdLepPt);
	p_scale->Write();
	outputfile->Write();
	outputfile->Close(); 
}
