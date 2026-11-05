#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLorentzVector.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

void plotGJet(){//main 

	TFile *file1 = TFile::Open("plot_eg_GJet-40to100.root");
	TH1F *p_dphi1 = (TH1F*)file1->Get("p_dphi");
	TH1F *p_proxy1= (TH1F*)file1->Get("p_proxy");
	TH1F *p_proxy_data1=(TH1F*)file1->Get("p_proxy_data");
	TH1F *p_phoEt1=(TH1F*)file1->Get("p_phoEt");
	TH1F *p_proxyEt1=(TH1F*)file1->Get("p_proxyEt");
	TH1F *p_proxyEt_data1=(TH1F*)file1->Get("p_proxyEt_data");
	
	TFile *file2 = TFile::Open("plot_eg_GJet-100To200.root");
	TH1F *p_dphi2 = (TH1F*)file2->Get("p_dphi");
	TH1F *p_proxy2= (TH1F*)file2->Get("p_proxy");
	TH1F *p_proxy_data2=(TH1F*)file2->Get("p_proxy_data");
	TH1F *p_phoEt2=(TH1F*)file2->Get("p_phoEt");
	TH1F *p_proxyEt2=(TH1F*)file2->Get("p_proxyEt");
	TH1F *p_proxyEt_data2=(TH1F*)file2->Get("p_proxyEt_data");


	TFile *file3 = TFile::Open("plot_eg_GJet-200To400.root");
	TH1F *p_dphi3 = (TH1F*)file3->Get("p_dphi");
	TH1F *p_proxy3= (TH1F*)file3->Get("p_proxy");
	TH1F *p_proxy_data3=(TH1F*)file3->Get("p_proxy_data");
	TH1F *p_phoEt3=(TH1F*)file3->Get("p_phoEt");
	TH1F *p_proxyEt3=(TH1F*)file3->Get("p_proxyEt");
	TH1F *p_proxyEt_data3=(TH1F*)file3->Get("p_proxyEt_data");

	TFile *file4 = TFile::Open("plot_eg_GJet-400To600.root");
	TH1F *p_dphi4 = (TH1F*)file4->Get("p_dphi");
	TH1F *p_proxy4= (TH1F*)file4->Get("p_proxy");
	TH1F *p_proxy_data4=(TH1F*)file4->Get("p_proxy_data");
	TH1F *p_phoEt4=(TH1F*)file4->Get("p_phoEt");
	TH1F *p_proxyEt4=(TH1F*)file4->Get("p_proxyEt");
	TH1F *p_proxyEt_data4=(TH1F*)file4->Get("p_proxyEt_data");


	TFile *file5 = TFile::Open("plot_eg_GJet-600ToInf.root");
	TH1F *p_dphi5 = (TH1F*)file5->Get("p_dphi");
	TH1F *p_proxy5= (TH1F*)file5->Get("p_proxy");
	TH1F *p_proxy_data5=(TH1F*)file5->Get("p_proxy_data");
	TH1F *p_phoEt5=(TH1F*)file5->Get("p_phoEt");
	TH1F *p_proxyEt5=(TH1F*)file5->Get("p_proxyEt");
	TH1F *p_proxyEt_data5=(TH1F*)file5->Get("p_proxyEt_data");

	TH1F *p_dphi = new TH1F("p_dphi","",32,0,3.2);
	TH1F *p_proxy= new TH1F("p_proxy","",32,0,3.2);
	TH1F *p_proxy_data= new TH1F("p_proxy_data","",32,0,3.2);
  TH1F *p_phoEt = new TH1F("p_phoEt","",100,0,300);
	TH1F *p_proxyEt=new TH1F("p_proxyEt","",100,0,300);
	TH1F *p_proxyEt_data=new TH1F("p_proxyEt_data","",100,0,300);	
	
	p_dphi->Add(p_dphi1);
	p_dphi->Add(p_dphi2);
	p_dphi->Add(p_dphi3);
	p_dphi->Add(p_dphi4);
	p_dphi->Add(p_dphi5);
	
	p_proxy->Add(p_proxy1);
	p_proxy->Add(p_proxy2);
	p_proxy->Add(p_proxy3);
	p_proxy->Add(p_proxy4);
	p_proxy->Add(p_proxy5);

	p_proxy_data->Add(p_proxy_data1);
	p_proxy_data->Add(p_proxy_data2);
	p_proxy_data->Add(p_proxy_data3);
	p_proxy_data->Add(p_proxy_data4);
	p_proxy_data->Add(p_proxy_data5);

	p_phoEt->Add(p_phoEt1);
	p_phoEt->Add(p_phoEt2);
	p_phoEt->Add(p_phoEt3);
	p_phoEt->Add(p_phoEt4);
	p_phoEt->Add(p_phoEt5);
	
	p_proxyEt->Add(p_proxyEt1);
	p_proxyEt->Add(p_proxyEt2);
	p_proxyEt->Add(p_proxyEt3);
	p_proxyEt->Add(p_proxyEt4);
	p_proxyEt->Add(p_proxyEt5);
	
	p_proxyEt_data->Add(p_proxyEt_data1);
	p_proxyEt_data->Add(p_proxyEt_data2);
	p_proxyEt_data->Add(p_proxyEt_data3);
	p_proxyEt_data->Add(p_proxyEt_data4);
	p_proxyEt_data->Add(p_proxyEt_data5);

	TCanvas *can1=new TCanvas("can1","",600,600);
	can1->cd();
  p_dphi->Draw();
	p_proxy->SetLineColor(kRed);
	p_proxy->Scale(p_dphi->Integral(1,32)/p_proxy->Integral(1,32));
	p_proxy->Draw("same");
	p_proxy_data->SetLineColor(kGreen);
	p_proxy_data->Scale(p_dphi->Integral(1,32)/p_proxy_data->Integral(1,32));
	p_proxy_data->Draw("same");

	TCanvas *can2=new TCanvas("can2","",600,600);
	can2->cd();
  p_phoEt->Draw();
	p_proxyEt->SetLineColor(kRed);
	p_proxyEt->Scale(p_phoEt->Integral(1,100)/p_proxyEt->Integral(1,100));
	p_proxyEt->Draw("same");
	p_proxyEt_data->SetLineColor(kGreen);
	p_proxyEt_data->Scale(p_phoEt->Integral(1,100)/p_proxyEt_data->Integral(1,100));
	p_proxyEt_data->Draw("same");

}


