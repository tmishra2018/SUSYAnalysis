#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooNumIntConfig.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

void
plotFit_menglei(){
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	TGraphErrors *p_frac = new TGraphErrors(4);
	TGraphErrors *p_frac_total = new TGraphErrors(1);
	TGraphErrors *p_error_total = new TGraphErrors(1);
	TH1F *p_frac_0 = new TH1F("p_frac_0","a_{V#gamma}",50,1,1.5);

	std::ifstream vgammascalefile("/uscms/homes/m/mengleis/work/SUSY2016/SUSYAnalysis/test/VGamma/fit/VGamma_scalefactor_mg.txt");
	float fakescale, vgammascale, fakescaleerror, vgammascaleerror;
	float leplow, lephigh;

	float highest(0), lowest(1);
	float fittingerror(0), systematicerror(0), totalerror(0);

	TCanvas *cantemp = new TCanvas("cantemp","",1200,1200);
	TH1D *temphist[4];
	temphist[0] = new TH1D("temphist0","",100,0.5,2);
	temphist[1] = new TH1D("temphist1","",100,0.5,2);
	temphist[2] = new TH1D("temphist2","",100,0.5,2);
	temphist[3] = new TH1D("temphist3","",100,0.5,2);
	cantemp->Divide(2,2);
	for(unsigned i(0);  i < 4; i++){
		double xvalue, xerror, yvalue, yerror;
		for(unsigned j(0);  j < 1000; j++){
			vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
			if(j==0){
				xvalue = (leplow + lephigh) < 200? (leplow + lephigh)/2 : (leplow + 200)/2;
				xerror = (lephigh - leplow) < 100? (lephigh - leplow)/2 : (200 - leplow)/2;
				yvalue = vgammascale;
			}
			temphist[i]->Fill(vgammascale);
		}
		cantemp->cd(i+1);
		temphist[i]->Draw();
		temphist[i]->Fit("gaus","","",0.5,2);
		yerror = temphist[i]->GetFunction("gaus")->GetParameter(2);
		if(yerror > 0.5)yerror= temphist[i]->GetRMS()/2;
		p_frac->SetPoint(i, xvalue, yvalue);
		p_frac->SetPointError(i, xerror, yerror);
	}
	double lowbound(2), highbound(0), norm(0);
	for(unsigned i(0);  i < 1000; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
		if(i == 0)p_frac_total->SetPoint(0, 100, vgammascale);
		if(i == 0)p_error_total->SetPoint(0, 100, vgammascale);
		if(i == 0)norm = vgammascale;
		if(i == 0)fittingerror = vgammascaleerror;
		p_frac_0->Fill(vgammascale);
		if(vgammascale < lowbound)lowbound = vgammascale;
		if(vgammascale > highbound)highbound = vgammascale;
	}
	p_frac_0->Fit("gaus");
	systematicerror = 3*p_frac_0->GetFunction("gaus")->GetParameter(2);
	//systematicerror = (highbound-norm) > (norm - lowbound)? (highbound-norm):(norm - lowbound); 
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror;
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror;
	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);
	p_frac_total->SetPointError(0,100,totalerror);
	p_error_total->SetPointError(0,100,totalerror);

	std::cout << "highbound = " << highbound << std::endl;
	std::cout << "lowbound = " << lowbound << std::endl;
	std::cout << "scale factor = " << p_frac_0->GetFunction("gaus")->GetParameter(1)  << std::endl;
	std::cout << "error = " << totalerror << std::endl;

	TCanvas *can = new TCanvas("can","",600,600);
	can->cd();
	TH1D *dummy = new TH1D("dummy","VGamma scale; p_{T} (GeV); scalefactor",1,0,200);
	dummy->SetMaximum(2);
	dummy->Draw();
	p_frac->SetLineColor(kBlue);
	p_frac->SetMarkerColor(kBlue);
	p_frac->SetMarkerStyle(20);
	p_frac->SetFillStyle(0);
	p_frac->SetFillColor(0);
	p_frac->Draw("EP same");
	p_frac_total->SetMarkerStyle(20);
	p_frac_total->SetFillStyle(0);
	p_frac_total->SetFillColor(0);
	p_frac_total->Draw("P same");
	p_error_total->SetLineColor(0);
	p_error_total->SetFillColor(15);
	p_error_total->SetFillStyle(3344);
	p_error_total->Draw("E2 same");
	TLegend *leg = new TLegend(0.5,0.75,0.85,0.9);
	leg->AddEntry(p_frac, "pt-binned VG scales");
	leg->AddEntry(p_frac_total, "VG scale in full range");
	leg->AddEntry(p_error_total, "Syst. Unc.");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	leg->Draw("same");
	can->SaveAs("scale_ptDependence_mg.pdf");	

	TCanvas *canscale = new TCanvas("canscale","",600,600);
	canscale->cd();
	p_frac_0->Draw();	
	canscale->SaveAs("VGammaScale_mg.pdf");	
}
