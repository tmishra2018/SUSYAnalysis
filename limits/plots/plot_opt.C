#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"

void plot_opt(){
	gROOT->SetBatch();
	TFile *file1 = TFile::Open("T5WG_ANbins_2016.root", "read");
	TGraph *g_expectedLimit_1 = (TGraph*)file1->Get("hExpLimMid");
	TH2D   *h_xSecLimit_1 = (TH2D*)file1->Get("obs");
	g_expectedLimit_1->SetLineColor(kBlack);
	g_expectedLimit_1->SetLineWidth(1);

	TFile *file2 = TFile::Open("T5WG_2016_preVFP.root", "read");
	TGraph *g_expectedLimit_2 = (TGraph*)file2->Get("hExpLimMid");
	TH2D   *h_xSecLimit_2 = (TH2D*)file2->Get("obs");
	g_expectedLimit_2->SetLineColor(kRed);
	g_expectedLimit_2->SetLineWidth(1);

	TFile *file3 = TFile::Open("T5WG_2017.root", "read");
	TGraph *g_expectedLimit_3 = (TGraph*)file3->Get("hExpLimMid");
	TH2D   *h_xSecLimit_3 = (TH2D*)file3->Get("obs");
	g_expectedLimit_3->SetLineColor(kMagenta);
	g_expectedLimit_3->SetLineWidth(1);

	TFile *file4 = TFile::Open("T5WG_2018.root", "read");
	TGraph *g_expectedLimit_4 = (TGraph*)file4->Get("hExpLimMid");
	TH2D   *h_xSecLimit_4 = (TH2D*)file4->Get("obs");
	g_expectedLimit_4->SetLineColor(kGreen);
	g_expectedLimit_4->SetLineWidth(1);

	TFile *file5 = TFile::Open("T5WG_Run2.root", "read");
	TGraph *g_expectedLimit_5 = (TGraph*)file5->Get("hExpLimMid");
	TH2D   *h_xSecLimit_5 = (TH2D*)file5->Get("obs");
	g_expectedLimit_5->SetLineColor(kCyan);
	g_expectedLimit_5->SetLineWidth(1);

	Double_t red[] = {0.5, 0.5, 1.0, 1.0, 1.0};
	Double_t green[] = {0.5, 1.0, 1.0, 0.6, 0.5};
	Double_t blue[] = {1.0, 1.0, 0.5, 0.4, 0.5};
	Double_t length[] = {0.0, 0.34, 0.61, 0.84, 1.0};
	TColor::CreateGradientColorTable(5, length, red, green, blue, 255);

	TCanvas *can_exclude = new TCanvas("can","",700,600);	
	can_exclude->cd();
	gPad->SetLogz();
	gPad->SetLeftMargin(0.15);
	h_xSecLimit_1->GetZaxis()->SetRangeUser(0.5,40);
	h_xSecLimit_1->Draw("colz");
	g_expectedLimit_1->Draw("same");
	g_expectedLimit_2->Draw("same");
	g_expectedLimit_3->Draw("same");
	g_expectedLimit_4->Draw("same");
	g_expectedLimit_5->Draw("same");

	g_expectedLimit_1->SetFillStyle(0);
	g_expectedLimit_2->SetFillStyle(0);
	g_expectedLimit_3->SetFillStyle(0);
	g_expectedLimit_4->SetFillStyle(0);
	g_expectedLimit_5->SetFillStyle(0);

	TLegend *leg=new TLegend(0.15,0.7,0.5,0.9);
	leg->AddEntry(g_expectedLimit_1, "Exp. limit (2016)");
	leg->AddEntry(g_expectedLimit_2, "2016 preVFP");
	leg->AddEntry(g_expectedLimit_3, "2017 ");
	leg->AddEntry(g_expectedLimit_4, "2018 ");
	leg->AddEntry(g_expectedLimit_5, "Run2 ");
	leg->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg->Draw("same");
	can_exclude->SaveAs("t5wg_checkLimits.pdf");	
/*
	TCanvas *can2 = new TCanvas("can2","",600,600);	
	can2->cd();
	gStyle->SetOptStat(0);
	TH2D *ratio_2 = new TH2D("ratio_2","plan2/plan1 ratio",23, 975.0, 2125.0, 80, 12.5, 2012.5);	
	for(unsigned ii(1); ii < h_xSecLimit_2->GetXaxis()->GetNbins() + 1; ii++){
		for(unsigned jj(0);  jj < h_xSecLimit_2->GetYaxis()->GetNbins() + 1; jj++){
			if(h_xSecLimit_1->GetBinContent(ii, jj) > 0)ratio_2->SetBinContent(ii, jj, h_xSecLimit_2->GetBinContent(ii,jj)/h_xSecLimit_1->GetBinContent(ii,jj));
		}
	}
	ratio_2->GetXaxis()->SetTitle("M_{#tilde{g}} (GeV)");
	ratio_2->GetYaxis()->SetTitle("M_{#tilde{#chi}_{1}^{0}} (GeV) ");
	ratio_2->GetZaxis()->SetRangeUser(0.1,1);
	ratio_2->GetZaxis()->SetLabelOffset(-0.003);
	ratio_2->GetZaxis()->SetLabelSize(0.03);
	ratio_2->GetZaxis()->SetTitleOffset(0.9);
	ratio_2->GetZaxis()->SetLabelSize(0.04);
	ratio_2->SetContour(255);
	ratio_2->Draw("colz");
	can2->SaveAs("opt2.pdf");	


	TCanvas *can3 = new TCanvas("can3","",600,600);	
	can3->cd();
	gStyle->SetOptStat(0);
	TH2D *ratio_3 = new TH2D("ratio_3","plan3/plan1 ratio",23, 975.0, 2125.0, 80, 12.5, 2012.5);	
	for(unsigned ii(1); ii < h_xSecLimit_3->GetXaxis()->GetNbins() + 1; ii++){
		for(unsigned jj(0);  jj < h_xSecLimit_3->GetYaxis()->GetNbins() + 1; jj++){
			if(h_xSecLimit_1->GetBinContent(ii, jj) > 0)ratio_3->SetBinContent(ii, jj, h_xSecLimit_3->GetBinContent(ii,jj)/h_xSecLimit_1->GetBinContent(ii,jj));
		}
	}
	ratio_3->GetXaxis()->SetTitle("M_{#tilde{g}} (GeV)");
	ratio_3->GetYaxis()->SetTitle("M_{#tilde{#chi}_{1}^{0}} (GeV) ");
	ratio_3->GetZaxis()->SetRangeUser(0.1,1);
	ratio_3->GetZaxis()->SetLabelOffset(-0.003);
	ratio_3->GetZaxis()->SetLabelSize(0.03);
	ratio_3->GetZaxis()->SetTitleOffset(0.9);
	ratio_3->GetZaxis()->SetLabelSize(0.04);
	ratio_3->SetContour(255);
	ratio_3->Draw("colz");
	can3->SaveAs("opt3.pdf");	


	TCanvas *can4 = new TCanvas("can4","",600,600);	
	can4->cd();
	gStyle->SetOptStat(0);
	TH2D *ratio_4 = new TH2D("ratio_4","plan4/plan1 ratio",23, 975.0, 2125.0, 80, 12.5, 2012.5);	
	for(unsigned ii(1); ii < h_xSecLimit_4->GetXaxis()->GetNbins() + 1; ii++){
		for(unsigned jj(0);  jj < h_xSecLimit_4->GetYaxis()->GetNbins() + 1; jj++){
			if(h_xSecLimit_1->GetBinContent(ii, jj) > 0)ratio_4->SetBinContent(ii, jj, h_xSecLimit_4->GetBinContent(ii,jj)/h_xSecLimit_1->GetBinContent(ii,jj));
		}
	}
	ratio_4->GetXaxis()->SetTitle("M_{#tilde{g}} (GeV)");
	ratio_4->GetYaxis()->SetTitle("M_{#tilde{#chi}_{1}^{0}} (GeV) ");
	ratio_4->GetZaxis()->SetRangeUser(0.1,1);
	ratio_4->GetZaxis()->SetLabelOffset(-0.003);
	ratio_4->GetZaxis()->SetLabelSize(0.03);
	ratio_4->GetZaxis()->SetTitleOffset(0.9);
	ratio_4->GetZaxis()->SetLabelSize(0.04);
	ratio_4->SetContour(255);
	ratio_4->Draw("colz");
	can4->SaveAs("opt4.pdf");	


	TCanvas *can5 = new TCanvas("can5","",600,600);	
	can5->cd();
	gStyle->SetOptStat(0);
	TH2D *ratio_5 = new TH2D("ratio_5","plan5/plan1 ratio",23, 975.0, 2125.0, 80, 12.5, 2012.5);	
	for(unsigned ii(1); ii < h_xSecLimit_5->GetXaxis()->GetNbins() + 1; ii++){
		for(unsigned jj(0);  jj < h_xSecLimit_5->GetYaxis()->GetNbins() + 1; jj++){
			if(h_xSecLimit_1->GetBinContent(ii, jj) > 0)ratio_5->SetBinContent(ii, jj, h_xSecLimit_5->GetBinContent(ii,jj)/h_xSecLimit_1->GetBinContent(ii,jj));
		}
	}
	ratio_5->GetXaxis()->SetTitle("M_{#tilde{g}} (GeV)");
	ratio_5->GetYaxis()->SetTitle("M_{#tilde{#chi}_{1}^{0}} (GeV) ");
	ratio_5->GetZaxis()->SetRangeUser(0.1,1);
	ratio_5->GetZaxis()->SetLabelOffset(-0.003);
	ratio_5->GetZaxis()->SetLabelSize(0.03);
	ratio_5->GetZaxis()->SetTitleOffset(0.9);
	ratio_5->GetZaxis()->SetLabelSize(0.04);
	ratio_5->SetContour(255);
	ratio_5->Draw("colz");
	can5->SaveAs("opt5.pdf");	*/
}
