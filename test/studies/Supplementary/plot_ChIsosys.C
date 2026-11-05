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
#include "../../../include/tdrstyle.C"
void plot_sys(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	setTDRStyle();

	TFile *file = TFile::Open("signalTree_mg_jetbkg.root");
	
	TH1D *p_MET_nom  = (TH1D*)file->Get("p_MET_nom"); 
	TH1D *p_MET_stat = (TH1D*)file->Get("p_MET_stat"); 
	TH1D *p_MET_alt  = (TH1D*)file->Get("p_MET_alt");  
	TH1D *p_HT_nom   = (TH1D*)file->Get("p_HT_nom");   
	TH1D *p_HT_stat  = (TH1D*)file->Get("p_HT_stat");  
	TH1D *p_HT_alt   = (TH1D*)file->Get("p_HT_alt");   
	TH1D *p_ET_nom   = (TH1D*)file->Get("p_ET_nom");   
	TH1D *p_ET_stat  = (TH1D*)file->Get("p_ET_stat");  
	TH1D *p_ET_alt   = (TH1D*)file->Get("p_ET_alt");  

	p_MET_nom->GetXaxis()->SetTitle("MET");
	p_HT_nom->GetXaxis()->SetTitle("HT");
	p_ET_nom->GetXaxis()->SetTitle("photon ET");
	p_MET_nom->SetTitle("MET");
	p_HT_nom->SetTitle("HT");
	p_ET_nom->SetTitle("photon ET");
  p_MET_nom->GetXaxis()->SetTitleSize(0.04);
	p_HT_nom->GetXaxis()->SetTitleSize(0.04);
  p_ET_nom->GetXaxis()->SetTitleSize(0.04);
	setTDRStyle();
	
	TH1D *p_flat_sys = new TH1D("p_flat_sys","",1,0,1000);
	p_flat_sys->SetBinContent(1,1);
	p_flat_sys->SetBinError(1,0.3);
  p_flat_sys->SetFillColor(12);
  p_flat_sys->SetFillStyle(3345);
	p_flat_sys->SetMarkerSize(0);
	
	TCanvas *can_met = new TCanvas("can_met","",600,600);
	setCanvas(can_met); 
	can_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_MET_nom->Draw();
	for(unsigned i(1); i < p_MET_alt->GetSize(); i++)p_MET_alt->SetBinError(i, sqrt(p_MET_stat->GetBinContent(i))/p_MET_stat->GetBinContent(i)*p_MET_alt->GetBinContent(i));
	p_MET_alt->SetLineColor(kRed);
	p_MET_alt->Draw("same");
	TLegend *leg = new TLegend(0.5,0.6,0.87,0.87);
  leg->AddEntry(p_MET_nom, "ChIso < 15 GeV");
  leg->AddEntry(p_MET_alt, "ChIso < 5 GeV");
  leg->Draw("same");
	
	can_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
	met_pad2->SetBottomMargin(0.3);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,1000,1);
	TH1D *unc_met = (TH1D*)p_MET_alt->Clone("unc_met");
	unc_met->GetXaxis()->SetTitle("MET");
  unc_met->GetYaxis()->SetRangeUser(0,2);
	unc_met->Divide(p_MET_nom);
	unc_met->Draw();
	flatratio_met->Draw("same");
	//p_flat_sys->Draw("E2 same");	



	TCanvas *can_ht = new TCanvas("can_ht","",600,600);
	setCanvas(can_ht); 
	can_ht->cd();
	TPad *ht_pad1 = new TPad("ht_pad1", "ht_pad1", 0, 0.3, 1, 1.0);
	setTopPad(ht_pad1); 
	ht_pad1->Draw();  
	ht_pad1->cd();  
	gPad->SetLogy();
	p_HT_nom->Draw();
	for(unsigned i(1); i < p_HT_alt->GetSize(); i++)p_HT_alt->SetBinError(i, sqrt(p_HT_stat->GetBinContent(i))/p_HT_stat->GetBinContent(i)*p_HT_alt->GetBinContent(i));
	p_HT_alt->SetLineColor(kRed);
	p_HT_alt->Draw("same");
  leg->Draw("same");
	
	can_ht->cd();
	TPad *ht_pad2 = new TPad("ht_pad2", "ht_pad2", 0, 0, 1, 0.3);
	ht_pad2->SetBottomMargin(0.3);
	ht_pad2->Draw();
	ht_pad2->cd();
  TLine *flatratio_ht = new TLine(0,1,1000,1);
	TH1D *unc_ht = (TH1D*)p_HT_alt->Clone("unc_ht");
	unc_ht->GetXaxis()->SetTitle("H_{T}");
  unc_ht->GetYaxis()->SetRangeUser(0,2);
	unc_ht->Divide(p_HT_nom);
	unc_ht->Draw();
	flatratio_ht->Draw("same");
	//p_flat_sys->Draw("E2 same");	


	TCanvas *can_et = new TCanvas("can_et","",600,600);
	setCanvas(can_et); 
	can_et->cd();
	TPad *et_pad1 = new TPad("et_pad1", "et_pad1", 0, 0.3, 1, 1.0);
	setTopPad(et_pad1); 
	et_pad1->Draw();  
	et_pad1->cd();  
	gPad->SetLogy();
	p_ET_nom->Draw();
	for(unsigned i(1); i < p_ET_alt->GetSize(); i++)p_ET_alt->SetBinError(i, sqrt(p_ET_stat->GetBinContent(i))/p_ET_stat->GetBinContent(i)*p_ET_alt->GetBinContent(i));
	p_ET_alt->SetLineColor(kRed);
	p_ET_alt->Draw("same");
  leg->Draw("same");
	
	can_et->cd();
	TPad *et_pad2 = new TPad("et_pad2", "et_pad2", 0, 0, 1, 0.3);
	et_pad2->SetBottomMargin(0.3);
	et_pad2->Draw();
	et_pad2->cd();
  TLine *flatratio_et = new TLine(0,1,1000,1);
	TH1D *unc_et = (TH1D*)p_ET_alt->Clone("unc_et");
	unc_et->GetXaxis()->SetTitle("photon E_{T}");
  unc_et->GetYaxis()->SetRangeUser(0,2);
	unc_et->Divide(p_ET_nom);
	unc_et->Draw();
	flatratio_et->Draw("same");
	//p_flat_sys->Draw("E2 same");	

	can_met->SaveAs("ChIso_mg_MET.pdf");
	can_ht->SaveAs("ChIso_mg_HT.pdf");
	can_et->SaveAs("ChIso_mg_ET.pdf");


}


