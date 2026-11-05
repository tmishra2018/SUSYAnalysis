#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"
#include "../../../include/analysis_scalefactor.h"
#include "../../../include/tdrstyle.C"

bool doEB=false;

void plot_ISRweight(){//main 

	setTDRStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetErrorX(0.5);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={0,50,100,150,200,250,300,800};
	TH1F *p_JetPt_data     = new TH1F("p_JetPt_data","",7,plotEtBins);
	TH1F *p_JetPt_ZG     = new TH1F("p_JetPt_ZG","",7,plotEtBins);
	TH1F *p_JetPt_rare     = new TH1F("p_JetPt_rare","",7,plotEtBins);
//************ Signal Tree **********************//
  TChain *tree = new TChain("ZTree");
  tree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_data.root");
  float phoEt(0);
  float phoEta(0);
  float dRPhoLep(0);
	float JetPt=0;

  tree->SetBranchAddress("phoEt",     &phoEt);
  tree->SetBranchAddress("phoEta",    &phoEta);
  tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	tree->SetBranchAddress("ISRJetPt",  &JetPt);

  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
		
		if(doEB && fabs(phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4))continue;

		if(JetPt > 799)JetPt = 799;
		if(dRPhoLep < 0.8)continue;
		double weight = 1;

		if( fabs(phoEta) < 1.4442){
			if(phoEt > 35 && phoEt < 40)weight = 1-0.225;            
			else if(phoEt > 40 && phoEt < 50)weight = 1- 0.188;      
			else if(phoEt > 50 && phoEt < 60)weight = 1- 0.181;      
			else if(phoEt > 60 && phoEt < 70)weight = 1- 0.135;      
			else if(phoEt > 70 && phoEt < 80)weight = 1- 0.117;      
			else if(phoEt > 80 && phoEt < 100)weight = 1- 0.054;     
			else if(phoEt > 100 && phoEt < 150)weight = 1- 0.065;    
			else if(phoEt > 150 && phoEt < 200)weight = 1- 0.071;    
			else if(phoEt > 200 && phoEt < 250)weight = 1- 0.011;    
			else weight=1;	                                         
		}	
		else{
			if(phoEt > 35 && phoEt < 40)weight = 1 - 0.275 ;
			else if(phoEt > 40 && phoEt < 50)weight = 1- 0.217 ;
			else if(phoEt > 50 && phoEt < 60)weight = 1- 0.126;
			else if(phoEt > 60 && phoEt < 70)weight = 1- 0.207;
			else if(phoEt > 70 && phoEt < 80)weight = 1- 0.175;
			else if(phoEt > 80 && phoEt < 100)weight = 1- 0.180;
			else if(phoEt > 100 && phoEt < 150)weight = 1- 0.033;
			else if(phoEt > 150 && phoEt < 200)weight = 1-0.03;
			else if(phoEt > 200 && phoEt < 250)weight = 1-0.03;
			else weight=1;	
		}

	 	p_JetPt_data->Fill(JetPt, weight);
	}//loop on  events

//************ Signal Tree **********************//
  TChain *ZGtree = new TChain("ZTree");
  ZGtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_ZG.root");
	float ZG_MCweight(0);
	float ZG_PUweight(0);
  float ZG_phoEt(0);
  float ZG_phoEta(0);
  float ZG_phoPhi(0);
  float ZG_lepPt(0);
  float ZG_lepEta(0);
	float ZG_bosonPt(0);
  float ZG_dRPhoLep(0);
	float ZG_JetPt(0);
  std::vector<int>   *ZG_mcPID=0;
  std::vector<float> *ZG_mcEta=0;
  std::vector<float> *ZG_mcPhi=0;
  std::vector<float> *ZG_mcPt=0;
  std::vector<int>   *ZG_mcMomPID=0;
  std::vector<int>   *ZG_mcGMomPID=0;

	ZGtree->SetBranchAddress("MCweight",  &ZG_MCweight); 
	ZGtree->SetBranchAddress("PUweight",  &ZG_PUweight);
  ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
  ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
  ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
  ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
	ZGtree->SetBranchAddress("bosonPt",   &ZG_bosonPt);
  ZGtree->SetBranchAddress("dRPhoLep",  &ZG_dRPhoLep);
	ZGtree->SetBranchAddress("JetPt",     &ZG_JetPt);
  ZGtree->SetBranchAddress("mcPID",     &ZG_mcPID);
  ZGtree->SetBranchAddress("mcEta",     &ZG_mcEta);
  ZGtree->SetBranchAddress("mcPhi",     &ZG_mcPhi);
  ZGtree->SetBranchAddress("mcPt",      &ZG_mcPt);
  ZGtree->SetBranchAddress("mcMomPID",  &ZG_mcMomPID);
  ZGtree->SetBranchAddress("mcGMomPID", &ZG_mcGMomPID);

  for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		if(doEB && fabs(ZG_phoEta) > 1.4442)continue;
		//if(!doEB && (fabs(ZG_phoEta)< 1.56 || fabs(ZG_phoEta) > 2.4))continue;

		double scalefactor = objectESF.getMuonESF(ZG_lepPt,ZG_lepEta)*objectESF.getPhotonESF(ZG_phoEt,ZG_phoEta)*objectESF.getMuonEGTRGESF(ZG_phoEt, ZG_lepPt);

		if(ZG_JetPt > 799)ZG_JetPt = 799;
		if(ZG_dRPhoLep < 0.8)continue;

		double weight = ZG_MCweight*ZG_PUweight*scalefactor;
		p_JetPt_ZG->Fill(ZG_JetPt, weight);	
	}//loop on  events


//************ Signal Tree **********************//
  TChain *raretree = new TChain("ZTree");
  raretree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_TTG.root");
  raretree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_TT.root");
  raretree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_WWG.root");
  raretree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_WZG.root");
	float rare_MCweight(0);
	float rare_PUweight(0);
  float rare_phoEt(0);
  float rare_phoEta(0);
  float rare_phoPhi(0);
  float rare_lepPt(0);
  float rare_lepEta(0);
  float rare_dRPhoLep(0);
	float rare_JetPt(0);
  std::vector<int>   *rare_mcPID=0;
  std::vector<float> *rare_mcEta=0;
  std::vector<float> *rare_mcPhi=0;
  std::vector<float> *rare_mcPt=0;
  std::vector<int>   *rare_mcMomPID=0;
  std::vector<int>   *rare_mcGMomPID=0;

	raretree->SetBranchAddress("MCweight",  &rare_MCweight); 
	raretree->SetBranchAddress("PUweight",  &rare_PUweight);
  raretree->SetBranchAddress("phoEt",     &rare_phoEt);
  raretree->SetBranchAddress("phoEta",    &rare_phoEta);
  raretree->SetBranchAddress("phoPhi",    &rare_phoPhi);
  raretree->SetBranchAddress("lepPt",     &rare_lepPt);
  raretree->SetBranchAddress("lepEta",    &rare_lepEta);
  raretree->SetBranchAddress("dRPhoLep",  &rare_dRPhoLep);
	raretree->SetBranchAddress("JetPt",     &rare_JetPt);
  raretree->SetBranchAddress("mcPID",    &rare_mcPID);
  raretree->SetBranchAddress("mcEta",    &rare_mcEta);
  raretree->SetBranchAddress("mcPhi",    &rare_mcPhi);
  raretree->SetBranchAddress("mcPt",     &rare_mcPt);
  raretree->SetBranchAddress("mcMomPID", &rare_mcMomPID);
  raretree->SetBranchAddress("mcGMomPID",&rare_mcGMomPID);

  for(unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		if(doEB && fabs(rare_phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(rare_phoEta) < 1.56|| fabs(rare_phoEta) > 2.4))continue;

		double scalefactor = scalefactor = objectESF.getMuonESF(rare_lepPt,rare_lepEta)*objectESF.getPhotonESF(rare_phoEt,rare_phoEta)*objectESF.getMuonEGTRGESF(rare_phoEt, rare_lepPt);

		if(rare_JetPt > 799)rare_JetPt = 799;
		if(rare_dRPhoLep < 0.8)continue;

		double weight = rare_MCweight*rare_PUweight*scalefactor;
		p_JetPt_rare->Fill(rare_JetPt, weight);
	}//loop on  events

	TCanvas *can_JetPt        = new TCanvas("can_JetPt",          "can_JetPt",600,600); 
	float scalefactor = 1.0; 
	double ISRwgt_norm[11];
	double ISRwgt_alter[11];
	for(int i(1); i<=11; i++){
		std::cout << "norm ratio " << i << " " << p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) 
							<< " stat " <<   p_JetPt_data->GetBinError(i)/p_JetPt_ZG->GetBinContent(i)
              << std::endl;
		ISRwgt_norm[i-1] = p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
	}

	gStyle->SetOptStat(0);
	can_JetPt->cd();
	TPad *JetPt_pad1 = new TPad("JetPt_pad1", "JetPt_pad1", 0, 0.3, 1, 1.0);
	JetPt_pad1->SetBottomMargin(0);
	JetPt_pad1->Draw();  
	JetPt_pad1->cd();  
	gStyle->SetOptStat(0);
	JetPt_pad1->SetLogy();
	p_JetPt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
	p_JetPt_data->SetMinimum(0.1);
	p_JetPt_data->SetLineColor(kBlack);
	p_JetPt_data->SetMarkerStyle(20);
	p_JetPt_data->SetMarkerColor(kBlack);
	p_JetPt_data->Draw("P");
	p_JetPt_ZG->Add(p_JetPt_rare);
	p_JetPt_ZG->SetLineColor(6);
	p_JetPt_ZG->SetFillStyle(1001);
	p_JetPt_ZG->SetFillColor(6);
	p_JetPt_ZG->Scale(scalefactor);
	p_JetPt_ZG->Draw("hist same");
	p_JetPt_rare->SetLineColor(8);
	p_JetPt_rare->SetFillStyle(1001);
	p_JetPt_rare->SetFillColor(8);
	p_JetPt_rare->Scale(scalefactor);
	p_JetPt_rare->Draw("hist same");
	p_JetPt_data->Draw("EP same");
	TLegend *leg_JetPt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_JetPt->SetFillStyle(0);
	leg_JetPt->AddEntry(p_JetPt_data, "Data");
	leg_JetPt->AddEntry(p_JetPt_ZG,   "Z#gamma");
	leg_JetPt->AddEntry(p_JetPt_rare, "t#bar{t}#gamma,WW#gamma,WZ#gamma");   
	leg_JetPt->Draw("same");

	can_JetPt->cd();
	TPad *JetPt_pad2 = new TPad("JetPt_pad2", "JetPt_pad2", 0, 0.05, 1, 0.3);
	JetPt_pad2->SetTopMargin(0);
	JetPt_pad2->SetBottomMargin(0.2);
	JetPt_pad2->Draw();
	JetPt_pad2->cd();
  TLine *flatratio_JetPt = new TLine(0,1,800,1);
	TH1F *ratio_JetPt=(TH1F*)p_JetPt_data->Clone("transfer factor");
	ratio_JetPt->SetMarkerStyle(20);
	ratio_JetPt->SetLineColor(kBlack);
	ratio_JetPt->GetXaxis()->SetRangeUser(0,800);
	ratio_JetPt->GetYaxis()->SetRangeUser(0,2);
	ratio_JetPt->SetMinimum(0);
	ratio_JetPt->SetMaximum(2);
	ratio_JetPt->Divide(p_JetPt_ZG);
	ratio_JetPt->SetTitle("");
	ratio_JetPt->GetYaxis()->SetTitle("Data/MC");
	ratio_JetPt->GetXaxis()->SetLabelFont(63);
	ratio_JetPt->GetXaxis()->SetLabelSize(14);
	ratio_JetPt->GetYaxis()->SetLabelFont(63);
	ratio_JetPt->GetYaxis()->SetLabelSize(14);
	ratio_JetPt->Draw();
	flatratio_JetPt->Draw("same");
	can_JetPt->SaveAs("PLOT_ISRweight.pdf");
}


