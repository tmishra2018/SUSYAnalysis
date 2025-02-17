#include "TChain.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TSystem.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TPad.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_tools.h"

#define NBIN 13
int ichannel = 2;

void plotJetFake(){//main 
	gROOT->SetBatch(kTRUE);
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	bool useMC(false);
	bool doCompare(true);
	 
	gSystem->Load("../../../lib/libAnaClasses.so");
	gStyle->SetOptStat(0);
	TGraphErrors *jetfrac_16preVFP = new TGraphErrors(NBIN);
	TGraphErrors *jetfracsys_16preVFP = new TGraphErrors(NBIN);
	TGraphErrors *jetfracTrue_16preVFP = new TGraphErrors(NBIN);
	TGraphErrors *ratio_16preVFP = new TGraphErrors(NBIN);
	
	TGraphErrors *jetfrac_16postVFP = new TGraphErrors(NBIN);
	TGraphErrors *jetfracsys_16postVFP = new TGraphErrors(NBIN);
	TGraphErrors *jetfracTrue_16postVFP = new TGraphErrors(NBIN);
	TGraphErrors *ratio_16postVFP = new TGraphErrors(NBIN);
	
	TGraphErrors *jetfrac_menglei = new TGraphErrors(NBIN);
	TGraphErrors *ratio_menglei= new TGraphErrors(NBIN);	

	TGraphErrors *jetfrac_17 = new TGraphErrors(NBIN);
	TGraphErrors *jetfracsys_17 = new TGraphErrors(NBIN);
	TGraphErrors *jetfracTrue_17 = new TGraphErrors(NBIN);
	TGraphErrors *ratio_17 = new TGraphErrors(NBIN);

	TGraphErrors *jetfrac_18 = new TGraphErrors(NBIN);
	TGraphErrors *jetfracsys_18 = new TGraphErrors(NBIN);
	TGraphErrors *jetfracTrue_18 = new TGraphErrors(NBIN);
	TGraphErrors *ratio_18 = new TGraphErrors(NBIN);
	int i(0);
	float pt_lower(0), pt_upper(0);
	float fakerate(0), error(0), systematic(0); 
	float truerate(0);

	cout<<endl<<endl<<"2016 pre VFP "<<endl<<endl;
	std::ifstream jetfake_file16preVFP;
	if(ichannel == 1)  jetfake_file16preVFP.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2016preVFP/JetFakeRate-DoubleEG-EB.txt");
	if(ichannel == 2)  jetfake_file16preVFP.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2016preVFP/JetFakeRate-MuonEG-EB.txt");
	if(jetfake_file16preVFP.is_open()){
		for(int i(0); i < NBIN; i++){ 
			if(!useMC)jetfake_file16preVFP >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			if(useMC)jetfake_file16preVFP >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			jetfrac_16preVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        cout<<fakerate<<"\t";
			jetfrac_16preVFP ->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
			jetfracsys_16preVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
			jetfracsys_16preVFP ->SetPointError(i, 0, systematic);
			if(useMC)jetfracTrue_16preVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, truerate);
			if(useMC)ratio_16preVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
			if(useMC)ratio_16preVFP ->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
		}
		jetfake_file16preVFP.close(); 
	}

	cout<<endl<<endl<<"2016 post VFP "<<endl<<endl;
	std::ifstream jetfake_file16postVFP;
	if(ichannel == 1)  jetfake_file16postVFP.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2016postVFP/JetFakeRate-DoubleEG-EB.txt");
	if(ichannel == 2)  jetfake_file16postVFP.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2016postVFP/JetFakeRate-MuonEG-EB.txt");
	if(jetfake_file16postVFP.is_open()){
		for(int i(0); i < NBIN; i++){ 
			if(!useMC)jetfake_file16postVFP >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			if(useMC)jetfake_file16postVFP >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			jetfrac_16postVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        cout<<fakerate<<"\t";
			jetfrac_16postVFP ->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
			jetfracsys_16postVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
			jetfracsys_16postVFP ->SetPointError(i, 0, systematic);
			if(useMC)jetfracTrue_16postVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, truerate);
			if(useMC)ratio_16postVFP ->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
			if(useMC)ratio_16postVFP ->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
		}
		jetfake_file16postVFP.close(); 
	}
	
	cout<<endl<<endl<<"2017 "<<endl<<endl;
	std::ifstream jetfake_file17;
	if(ichannel == 1)  jetfake_file17.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2017/JetFakeRate-DoubleEG-EB.txt");
	if(ichannel == 2)  jetfake_file17.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2017/JetFakeRate-MuonEG-EB.txt");
	if(jetfake_file17.is_open()){
		for(int i(0); i < NBIN; i++){ 
			if(!useMC)jetfake_file17 >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			if(useMC)jetfake_file17 >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			jetfrac_17->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        cout<<fakerate<<"\t";
			jetfrac_17->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
			jetfracsys_17->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
			jetfracsys_17->SetPointError(i, 0, systematic);
			if(useMC)jetfracTrue_17->SetPoint(i, (pt_lower+pt_upper)/2.0, truerate);
			if(useMC)ratio_17->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
			if(useMC)ratio_17->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
		}
		jetfake_file17.close(); 
	}
	
	cout<<endl<<endl<<"2018 "<<endl<<endl;
	std::ifstream jetfake_file18;
	if(ichannel == 1)  jetfake_file18.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2018/JetFakeRate-DoubleEG-EB.txt");
	if(ichannel == 2)  jetfake_file18.open("/eos/uscms/store/user/tmishra/jetfakepho/txt2018/JetFakeRate-MuonEG-EB.txt");
	if(jetfake_file18.is_open()){
		for(int i(0); i < NBIN; i++){ 
			if(!useMC)jetfake_file18 >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			if(useMC)jetfake_file18 >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			jetfrac_18->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        cout<<fakerate<<"\t";
			jetfrac_18->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
			jetfracsys_18->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
			jetfracsys_18->SetPointError(i, 0, systematic);
			if(useMC)jetfracTrue_18->SetPoint(i, (pt_lower+pt_upper)/2.0, truerate);
			if(useMC)ratio_18->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
			if(useMC)ratio_18->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
		}
		jetfake_file18.close(); 
	}
	
	cout<<endl<<endl<<"menglei "<<endl<<endl;
	if(doCompare){
		std::ifstream jetfake_menglei;
		if(ichannel == 1)  jetfake_menglei.open("/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-DoubleEG-EB.txt");
		if(ichannel == 2)  jetfake_menglei.open("/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-MuonEG-EB.txt");
		if(jetfake_menglei.is_open()){
			for(int i(0); i < NBIN; i++){
				if(!useMC)jetfake_menglei >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
				if(useMC)jetfake_menglei >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
				jetfrac_menglei->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        	cout<<fakerate<<"\t";
				jetfrac_menglei->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
				if(useMC)ratio_menglei->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
				if(useMC)ratio_menglei->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
			}
			jetfake_menglei.close();
		}
	}
	cout<<endl;

	jetfrac_16preVFP->SetMarkerStyle(20);
	jetfrac_16preVFP->SetMarkerColor(kRed);
	jetfrac_16preVFP->SetLineColor(kRed);
	jetfracsys_16preVFP->SetMarkerStyle(20);
	jetfracsys_16preVFP->SetMarkerColor(kRed);
	jetfracsys_16preVFP->SetLineColor(kRed);
	jetfrac_16preVFP->SetLineWidth(2);
	jetfracsys_16preVFP->SetLineWidth(2);
	jetfrac_16preVFP->SetFillColor(0);
	jetfracTrue_16preVFP->SetFillColor(0);
	jetfracTrue_16preVFP->SetLineColor(kRed);

	jetfrac_menglei->SetFillColor(0);
	jetfrac_menglei->SetMarkerStyle(20);
	jetfrac_menglei->SetMarkerColor(kBlack);
	jetfrac_menglei->SetLineColor(kBlack);

	jetfrac_16postVFP->SetMarkerStyle(20);
	jetfrac_16postVFP->SetMarkerColor(kViolet);
	jetfrac_16postVFP->SetLineColor(kViolet);
	
	jetfrac_17->SetMarkerStyle(20);
	jetfrac_17->SetMarkerColor(kBlue);
	jetfrac_17->SetLineColor(kBlue);
	
	jetfrac_18->SetMarkerStyle(20);
	jetfrac_18->SetMarkerColor(kGreen);
	jetfrac_18->SetLineColor(kGreen);
	TCanvas *canvas = new TCanvas("Hadron Fraction"," fake rate",600,600);
	canvas->cd();
//	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
//	pad1->SetBottomMargin(0);
//	pad1->Draw();  
//	pad1->cd();  
	TH1F *dummy = new TH1F("Hadron Fraction","e#gamma channel;p_{T}(GeV);hadron fraction",17,30,200);
	dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	dummy->GetXaxis()->SetTitleOffset(1);
	dummy->SetTitleOffset(1);
	dummy->SetTitle("e#gamma channel");
	if(ichannel == 1) dummy->SetMaximum(0.6);
	if(ichannel == 2) dummy->SetMaximum(0.3);
	dummy->Draw();
	jetfrac_16preVFP->Draw("P same");
	jetfrac_16postVFP->Draw("P same");
//	jetfracsys->Draw("P same");
	if(doCompare)jetfrac_menglei->SetFillStyle(0);
//	if(doCompare)jetfrac_menglei->Draw("P same");
	if(doCompare)jetfrac_17->Draw("P same");
	if(doCompare)jetfrac_18->Draw("P same");
	//if(useMC)jetfracTrue->Draw("P same");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	TLegend *leg =  new TLegend(0.6,0.75,0.85,0.9);
	jetfrac_16preVFP->SetFillStyle(0);
	leg->SetFillStyle(0);
//	leg->AddEntry(jetfrac_menglei, "2016 : SUS-17-012");
	leg->AddEntry(jetfrac_16preVFP, "2016 preVFP");
	leg->AddEntry(jetfrac_16postVFP, "2016 postVFP");
	leg->AddEntry(jetfrac_17, "2017");
	leg->AddEntry(jetfrac_18, "2018");
	//if(useMC)leg->AddEntry(jetfracTrue,"MC truth");
	leg->Draw("same");
	CMS_lumi( canvas, 7,ichannel, 11 );
//	canvas->cd();
//	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
//	pad2->Draw();
//	pad2->cd();
//	pad2->SetTopMargin(0);
//	pad2->SetBottomMargin(0.3);
//	TH1F *dummy2 = new TH1F("dummy2",";p_{T}(GeV);estimated/true",17,30,200);
//	dummy2->SetMinimum(0);
//	dummy2->SetMaximum(5);
//	dummy2->SetMarkerStyle(20);
//	dummy2->Draw();
//  TLine *flatratio = new TLine(30,1,200,1);
//	ratio->Draw("P same");
//	if(doCompare)ratio_menglei->SetMarkerStyle(20);
//	if(doCompare)ratio_menglei->SetMarkerColor(kBlue);
//	if(doCompare)ratio_menglei->Draw("P same");
//	flatratio->Draw("same");
	if(ichannel == 1)  	canvas->SaveAs("/eos/uscms/store/user/tmishra/jetfakepho/Plots2016/JetFakePho_DoubleEG_comparison.pdf");
	if(ichannel == 2)  	canvas->SaveAs("/eos/uscms/store/user/tmishra/jetfakepho/Plots2016/JetFakePho_MuonEG_comparison.pdf");
}
