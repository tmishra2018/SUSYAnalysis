// input is VGamma/src/analysis_Mixing.C
// g++ `root-config --cflags` plot_Mixing.C -o plot_Mixing.exe `root-config --libs`
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile2D.h"
#include "../../include/tdrstyle.C"
bool preVFP = true;

float lumi_2016preVFP_DoubleEG = 19.499256;
float lumi_2016postVFP_DoubleEG = 16.810813;
float lumi_2017_DoubleEG = 41.48;
float lumi_2018_DoubleEG = 59.816229;

float lumi_2016preVFP_MuonEG = 19.495077;
float lumi_2016postVFP_MuonEG = 16.810813;
float lumi_2017_MuonEG = 41.48;
float lumi_2018_MuonEG = 59.810698;


float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

void plot_Mixing(int RunYear){//main 
	
	std::string whichVFP;
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	// for pt<50, pt>50 pt > 130 samples
  	double scalefactor1,scalefactor2,scalefactor3;
  	double lumi_DoubleEG, lumi_MuonEG;
  	if(RunYear==2016 and preVFP == 1){
		lumi_DoubleEG = lumi_2016preVFP_DoubleEG;
		lumi_MuonEG = lumi_2016preVFP_MuonEG;}
  	else if(RunYear==2016 and preVFP == 0){
		lumi_DoubleEG = lumi_2016postVFP_DoubleEG;
		lumi_MuonEG = lumi_2016postVFP_MuonEG;}
  	else if(RunYear==2017){
		lumi_DoubleEG = lumi_2017_DoubleEG;
		lumi_MuonEG = lumi_2017_MuonEG;}
   	else if(RunYear==2018){
		lumi_DoubleEG = lumi_2018_DoubleEG;
		lumi_MuonEG = lumi_2018_MuonEG;}

    	//gROOT->SetBatch(kTRUE);
  	TH1D *mugamma_phoEt_1 = new TH1D("mugamma_phoEt_1","",165,35,200); 
  	TH1D *mugamma_phoEt_2 = new TH1D("mugamma_phoEt_2","",165,35,200); 
  	TH1D *mugamma_phoEt_3 = new TH1D("mugamma_phoEt_3","",165,35,200); 
	TH1D *egamma_phoEt_1 =  new TH1D("egamma_phoEt_1","", 165,35,200); 
	TH1D *egamma_phoEt_2 =  new TH1D("egamma_phoEt_2","", 165,35,200); 
	TH1D *egamma_phoEt_3 =  new TH1D("egamma_phoEt_3","", 165,35,200); 
  	TH1D *mugamma_phoEt_total = new TH1D("mugamma_phoEt_total",";p_{T} (GeV);Events",165,35,200); 
	TH1D *egamma_phoEt_total  = new TH1D("egamma_phoEt_total",";p_{T} (GeV);Events", 165,35,200); 

  	mugamma_phoEt_1->Sumw2();
  	mugamma_phoEt_2->Sumw2();
  	mugamma_phoEt_3->Sumw2();
	egamma_phoEt_1->Sumw2();
	egamma_phoEt_2->Sumw2();
	egamma_phoEt_3->Sumw2();
  	mugamma_phoEt_total->Sumw2();
	egamma_phoEt_total->Sumw2();
  	// WGToLNuG tree
  	TChain *mgtree = new TChain("mgTree");
	mgtree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg_phoEt=0;
  	float mg_phoEta=0;
  	float mg_phoPhi=0;
  	float mg_ntotalevent=0;
  	float mg_crosssection=0;

  	mgtree->SetBranchAddress("phoEt",    &mg_phoEt);
  	mgtree->SetBranchAddress("phoEta",   &mg_phoEta);
  	mgtree->SetBranchAddress("phoPhi",   &mg_phoPhi);
  	mgtree->SetBranchAddress("ntotalevent",   &mg_ntotalevent);
  	mgtree->SetBranchAddress("crosssection",   &mg_crosssection);

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);
		if(mg_phoEt > 50)continue;
		//if(fabs(mg_phoEta) > 1.4442)continue;
		mugamma_phoEt_1->Fill(mg_phoEt, lumi_MuonEG*1000*mg_crosssection/mg_ntotalevent);	
		mugamma_phoEt_total->Fill(mg_phoEt, lumi_MuonEG*1000*mg_crosssection/mg_ntotalevent);	
	}

  	// WGJets_PtG-40-130 tree
  	TChain *mg40tree = new TChain("mgTree");
	mg40tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg40_phoEt=0;
  	float mg40_phoEta=0;
  	float mg40_phoPhi=0;
  	float mg40_ntotalevent=0;
  	float mg40_crosssection=0;

  	mg40tree->SetBranchAddress("phoEt",    &mg40_phoEt);
  	mg40tree->SetBranchAddress("phoEta",   &mg40_phoEta);
  	mg40tree->SetBranchAddress("phoPhi",   &mg40_phoPhi);
  	mg40tree->SetBranchAddress("ntotalevent",   &mg40_ntotalevent);
  	mg40tree->SetBranchAddress("crosssection",   &mg40_crosssection);

	for(unsigned ievt(0); ievt < mg40tree->GetEntries(); ievt++){
		mg40tree->GetEntry(ievt);
		if(mg40_phoEt <= 50)continue;
		//if(fabs(mg40_phoEta) > 1.4442)continue;
		mugamma_phoEt_2->Fill(mg40_phoEt, lumi_MuonEG*1000*mg40_crosssection/mg40_ntotalevent);			
		mugamma_phoEt_total->Fill(mg40_phoEt, lumi_MuonEG*1000*mg40_crosssection/mg40_ntotalevent);
	}

  	// WGJets_PtG-130 tree
  	TChain *mg130tree = new TChain("mgTree");
	mg130tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg130_phoEt=0;
  	float mg130_phoEta=0;
  	float mg130_phoPhi=0;
  	float mg130_ntotalevent=0;
  	float mg130_crosssection=0;

  	mg130tree->SetBranchAddress("phoEt",    &mg130_phoEt);
  	mg130tree->SetBranchAddress("phoEta",   &mg130_phoEta);
  	mg130tree->SetBranchAddress("phoPhi",   &mg130_phoPhi);
  	mg130tree->SetBranchAddress("ntotalevent",   &mg130_ntotalevent);
  	mg130tree->SetBranchAddress("crosssection",   &mg130_crosssection);

	for(unsigned ievt(0); ievt < mg130tree->GetEntries(); ievt++){
		mg130tree->GetEntry(ievt);
		//if(fabs(mg130_phoEta) > 1.4442)continue;
		mugamma_phoEt_3->Fill(mg130_phoEt, lumi_MuonEG*1000*mg130_crosssection/mg130_ntotalevent);			
		mugamma_phoEt_total->Fill(mg130_phoEt, lumi_MuonEG*1000*mg130_crosssection/mg130_ntotalevent);
	}

  	TChain *egtree = new TChain("egTree");
	egtree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg_phoEt=0;
  	float eg_phoEta=0;
  	float eg_phoPhi=0;
  	float eg_ntotalevent=0;
  	float eg_crosssection=0;

  	egtree->SetBranchAddress("phoEt",    &eg_phoEt);
  	egtree->SetBranchAddress("phoEta",   &eg_phoEta);
  	egtree->SetBranchAddress("phoPhi",   &eg_phoPhi);
  	egtree->SetBranchAddress("ntotalevent",   &eg_ntotalevent);
  	egtree->SetBranchAddress("crosssection",   &eg_crosssection);
  	//cout<<egtree->GetEntries()<<endl;

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);
		//if(eg_phoEt > 50)continue;
		//if(fabs(eg_phoEta) > 1.4442)continue;
		egamma_phoEt_1->Fill(eg_phoEt, lumi_DoubleEG*1000*eg_crosssection/eg_ntotalevent);		
		egamma_phoEt_total->Fill(eg_phoEt, lumi_DoubleEG*1000*eg_crosssection/eg_ntotalevent);	
	}


  	TChain *eg40tree = new TChain("egTree");
	eg40tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg40_phoEt=0;
  	float eg40_phoEta=0;
  	float eg40_phoPhi=0;
  	float eg40_ntotalevent=0;
  	float eg40_crosssection=0; 

  	eg40tree->SetBranchAddress("phoEt",    &eg40_phoEt);
  	eg40tree->SetBranchAddress("phoEta",   &eg40_phoEta);
  	eg40tree->SetBranchAddress("phoPhi",   &eg40_phoPhi);
  	eg40tree->SetBranchAddress("ntotalevent",   &eg40_ntotalevent);
  	eg40tree->SetBranchAddress("crosssection",   &eg40_crosssection);
  	//cout<<eg40tree->GetEntries()<<endl;
	for(unsigned ievt(0); ievt < eg40tree->GetEntries(); ievt++){
		eg40tree->GetEntry(ievt);
	//	if(eg40_phoEt <= 50)continue;
		//if(fabs(eg40_phoEta) > 1.4442)continue;
		egamma_phoEt_2->Fill(eg40_phoEt, lumi_DoubleEG*1000*eg40_crosssection/eg40_ntotalevent);			
		egamma_phoEt_total->Fill(eg40_phoEt, lumi_DoubleEG*1000*eg40_crosssection/eg40_ntotalevent);
	}

  	TChain *eg130tree = new TChain("egTree");
	eg130tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg130_phoEt=0;
  	float eg130_phoEta=0;
  	float eg130_phoPhi=0;
  	float eg130_ntotalevent=0;
  	float eg130_crosssection=0;

  	eg130tree->SetBranchAddress("phoEt",    &eg130_phoEt);
  	eg130tree->SetBranchAddress("phoEta",   &eg130_phoEta);
  	eg130tree->SetBranchAddress("phoPhi",   &eg130_phoPhi);
  	eg130tree->SetBranchAddress("ntotalevent",   &eg130_ntotalevent);
  	eg130tree->SetBranchAddress("crosssection",   &eg130_crosssection);
  	//cout<<eg130tree->GetEntries()<<endl;

	for(unsigned ievt(0); ievt < eg130tree->GetEntries(); ievt++){
		eg130tree->GetEntry(ievt);
		//if(fabs(eg130_phoEta) > 1.4442)continue;
		egamma_phoEt_3->Fill(eg130_phoEt, lumi_DoubleEG*1000*eg130_crosssection/eg130_ntotalevent);			
		egamma_phoEt_total->Fill(eg130_phoEt, lumi_DoubleEG*1000*eg130_crosssection/eg130_ntotalevent);
	}

  	double ratio1 = egamma_phoEt_1->Integral(6,8)/egamma_phoEt_2->Integral(6,8);
  	double ratio2 = mugamma_phoEt_1->Integral(6,8)/mugamma_phoEt_2->Integral(6,8);

  	double ratio1_130 = egamma_phoEt_1->Integral(15,30)/egamma_phoEt_3->Integral(15,30);
  	double ratio2_130 = mugamma_phoEt_1->Integral(15,30)/mugamma_phoEt_3->Integral(15,30);

	gStyle->SetOptStat(0);
	TCanvas *can1 = new TCanvas("can1","",600,600);
  	can1->cd();
	gPad->SetLogy();
	egamma_phoEt_total->SetFillColor(17);
	egamma_phoEt_total->SetLineColor(17);
	egamma_phoEt_total->SetFillStyle(1001);
	egamma_phoEt_total->Draw("hist");
	egamma_phoEt_1->SetMarkerStyle(7);
	egamma_phoEt_1->SetMarkerColor(kBlue);
	egamma_phoEt_1->SetLineColor(kBlue);
	egamma_phoEt_1->Draw("EP same");
	egamma_phoEt_2->SetMarkerStyle(7);
	egamma_phoEt_2->SetMarkerColor(kRed);
	egamma_phoEt_2->SetLineColor(kRed);
	egamma_phoEt_2->Draw("EP same");
	egamma_phoEt_3->SetMarkerStyle(7);
	egamma_phoEt_3->SetMarkerColor(kGreen);
	egamma_phoEt_3->SetLineColor(kGreen);
	egamma_phoEt_3->Draw("EP same");
	TLegend *leg =  new TLegend(0.4,0.6,0.87,0.9);
	leg->SetTextSize(0.025);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	std::cout<<egamma_phoEt_total->Integral()<<std::endl;
	leg->AddEntry(egamma_phoEt_total, "mixed W#gamma");
	leg->AddEntry(egamma_phoEt_1,"WGToLNuG");
	leg->AddEntry(egamma_phoEt_2,"WGJets_MonoPhoton_PtG-40to130");
	leg->AddEntry(egamma_phoEt_3,"WGJets_MonoPhoton_PtG-130");
	leg->Draw("same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( can1,1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( can1,2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( can1,3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( can1,4,1, 11 );

	can1->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_%d%s.pdf",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	can1->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
}
int main(int argc, char** argv)
{
    plot_Mixing(atoi(argv[1]));
    return 0;
}
