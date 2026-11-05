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
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/tdrstyle.C"

void plot_Mixing(int RunYear, bool preVFP){//main 
	
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

	// -------------------------------------------------------
	// Histogram declarations — egamma channel
	// -------------------------------------------------------
	TH1D *egamma_phoEt_1 =  new TH1D("egamma_phoEt_1","", 165,35,200); 
	TH1D *egamma_phoEt_2 =  new TH1D("egamma_phoEt_2","", 165,35,200); 
	TH1D *egamma_phoEt_3 =  new TH1D("egamma_phoEt_3","", 165,35,200); 
	TH1D *egamma_phoEt_total = new TH1D("egamma_phoEt_total",";photon p_{T} (GeV);Events", 165,35,200); 

	egamma_phoEt_1->Sumw2();
	egamma_phoEt_2->Sumw2();
	egamma_phoEt_3->Sumw2();
	egamma_phoEt_total->Sumw2();

	// -------------------------------------------------------
	// Histogram declarations — mugamma channel
	// -------------------------------------------------------
  	TH1D *mugamma_phoEt_1 = new TH1D("mugamma_phoEt_1","",165,35,200); 
  	TH1D *mugamma_phoEt_2 = new TH1D("mugamma_phoEt_2","",165,35,200); 
  	TH1D *mugamma_phoEt_3 = new TH1D("mugamma_phoEt_3","",165,35,200); 
  	TH1D *mugamma_phoEt_total = new TH1D("mugamma_phoEt_total",";photon p_{T} (GeV);Events",165,35,200); 

  	mugamma_phoEt_1->Sumw2();
  	mugamma_phoEt_2->Sumw2();
  	mugamma_phoEt_3->Sumw2();
  	mugamma_phoEt_total->Sumw2();

	// =======================================================
	// WGToLNuG tree  (pt < 50 cut for mugamma)
	// =======================================================
  	TChain *mgtree = new TChain("mgTree");
	mgtree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg_phoEt=0, mg_phoEta=0, mg_phoPhi=0, mg_ntotalevent=0, mg_crosssection=0;

  	mgtree->SetBranchAddress("phoEt",        &mg_phoEt);
  	mgtree->SetBranchAddress("phoEta",        &mg_phoEta);
  	mgtree->SetBranchAddress("phoPhi",        &mg_phoPhi);
  	mgtree->SetBranchAddress("ntotalevent",   &mg_ntotalevent);
  	mgtree->SetBranchAddress("crosssection",  &mg_crosssection);

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);
		if(mg_phoEt > 50) continue;
		mugamma_phoEt_1->Fill(mg_phoEt, lumi_MuonEG*1000*mg_crosssection/mg_ntotalevent);	
		mugamma_phoEt_total->Fill(mg_phoEt, lumi_MuonEG*1000*mg_crosssection/mg_ntotalevent);	
	}

	// =======================================================
	// WGJets_PtG-40-130 tree  (pt > 50 cut for mugamma)
	// =======================================================
  	TChain *mg40tree = new TChain("mgTree");
	mg40tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg40_phoEt=0, mg40_phoEta=0, mg40_phoPhi=0, mg40_ntotalevent=0, mg40_crosssection=0;

  	mg40tree->SetBranchAddress("phoEt",        &mg40_phoEt);
  	mg40tree->SetBranchAddress("phoEta",        &mg40_phoEta);
  	mg40tree->SetBranchAddress("phoPhi",        &mg40_phoPhi);
  	mg40tree->SetBranchAddress("ntotalevent",   &mg40_ntotalevent);
  	mg40tree->SetBranchAddress("crosssection",  &mg40_crosssection);

	for(unsigned ievt(0); ievt < mg40tree->GetEntries(); ievt++){
		mg40tree->GetEntry(ievt);
		if(mg40_phoEt <= 50) continue;
		mugamma_phoEt_2->Fill(mg40_phoEt, lumi_MuonEG*1000*mg40_crosssection/mg40_ntotalevent);			
		mugamma_phoEt_total->Fill(mg40_phoEt, lumi_MuonEG*1000*mg40_crosssection/mg40_ntotalevent);
	}

	// =======================================================
	// WGJets_PtG-130 tree  (mugamma, no pt cut)
	// =======================================================
  	TChain *mg130tree = new TChain("mgTree");
	mg130tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
  	float mg130_phoEt=0, mg130_phoEta=0, mg130_phoPhi=0, mg130_ntotalevent=0, mg130_crosssection=0;

  	mg130tree->SetBranchAddress("phoEt",        &mg130_phoEt);
  	mg130tree->SetBranchAddress("phoEta",        &mg130_phoEta);
  	mg130tree->SetBranchAddress("phoPhi",        &mg130_phoPhi);
  	mg130tree->SetBranchAddress("ntotalevent",   &mg130_ntotalevent);
  	mg130tree->SetBranchAddress("crosssection",  &mg130_crosssection);

	for(unsigned ievt(0); ievt < mg130tree->GetEntries(); ievt++){
		mg130tree->GetEntry(ievt);
		mugamma_phoEt_3->Fill(mg130_phoEt, lumi_MuonEG*1000*mg130_crosssection/mg130_ntotalevent);			
		mugamma_phoEt_total->Fill(mg130_phoEt, lumi_MuonEG*1000*mg130_crosssection/mg130_ntotalevent);
	}

	// =======================================================
	// egTree — WGToLNuG  (egamma channel, no pt cut)
	// =======================================================
  	TChain *egtree = new TChain("egTree");
	egtree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg_phoEt=0, eg_phoEta=0, eg_phoPhi=0, eg_ntotalevent=0, eg_crosssection=0;

  	egtree->SetBranchAddress("phoEt",        &eg_phoEt);
  	egtree->SetBranchAddress("phoEta",        &eg_phoEta);
  	egtree->SetBranchAddress("phoPhi",        &eg_phoPhi);
  	egtree->SetBranchAddress("ntotalevent",   &eg_ntotalevent);
  	egtree->SetBranchAddress("crosssection",  &eg_crosssection);

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);
		egamma_phoEt_1->Fill(eg_phoEt, lumi_DoubleEG*1000*eg_crosssection/eg_ntotalevent);		
		egamma_phoEt_total->Fill(eg_phoEt, lumi_DoubleEG*1000*eg_crosssection/eg_ntotalevent);	
	}

	// =======================================================
	// egTree — WGJets_PtG-40-130  (egamma channel, no pt cut)
	// =======================================================
  	TChain *eg40tree = new TChain("egTree");
	eg40tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg40_phoEt=0, eg40_phoEta=0, eg40_phoPhi=0, eg40_ntotalevent=0, eg40_crosssection=0; 

  	eg40tree->SetBranchAddress("phoEt",        &eg40_phoEt);
  	eg40tree->SetBranchAddress("phoEta",        &eg40_phoEta);
  	eg40tree->SetBranchAddress("phoPhi",        &eg40_phoPhi);
  	eg40tree->SetBranchAddress("ntotalevent",   &eg40_ntotalevent);
  	eg40tree->SetBranchAddress("crosssection",  &eg40_crosssection);

	for(unsigned ievt(0); ievt < eg40tree->GetEntries(); ievt++){
		eg40tree->GetEntry(ievt);
		egamma_phoEt_2->Fill(eg40_phoEt, lumi_DoubleEG*1000*eg40_crosssection/eg40_ntotalevent);			
		egamma_phoEt_total->Fill(eg40_phoEt, lumi_DoubleEG*1000*eg40_crosssection/eg40_ntotalevent);
	}

	// =======================================================
	// egTree — WGJets_PtG-130  (egamma channel, no pt cut)
	// =======================================================
  	TChain *eg130tree = new TChain("egTree");
	eg130tree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
  	float eg130_phoEt=0, eg130_phoEta=0, eg130_phoPhi=0, eg130_ntotalevent=0, eg130_crosssection=0;

  	eg130tree->SetBranchAddress("phoEt",        &eg130_phoEt);
  	eg130tree->SetBranchAddress("phoEta",        &eg130_phoEta);
  	eg130tree->SetBranchAddress("phoPhi",        &eg130_phoPhi);
  	eg130tree->SetBranchAddress("ntotalevent",   &eg130_ntotalevent);
  	eg130tree->SetBranchAddress("crosssection",  &eg130_crosssection);

	for(unsigned ievt(0); ievt < eg130tree->GetEntries(); ievt++){
		eg130tree->GetEntry(ievt);
		egamma_phoEt_3->Fill(eg130_phoEt, lumi_DoubleEG*1000*eg130_crosssection/eg130_ntotalevent);			
		egamma_phoEt_total->Fill(eg130_phoEt, lumi_DoubleEG*1000*eg130_crosssection/eg130_ntotalevent);
	}

	// =======================================================
	// Overlap ratios (kept from original)
	// =======================================================
  	double ratio1     = egamma_phoEt_1->Integral(6,8)  / egamma_phoEt_2->Integral(6,8);
  	double ratio2     = mugamma_phoEt_1->Integral(6,8) / mugamma_phoEt_2->Integral(6,8);
  	double ratio1_130 = egamma_phoEt_1->Integral(15,30) / egamma_phoEt_3->Integral(15,30);
  	double ratio2_130 = mugamma_phoEt_1->Integral(15,30)/ mugamma_phoEt_3->Integral(15,30);

	gStyle->SetOptStat(0);

	// =======================================================
	// Canvas 1: egamma channel
	// =======================================================
	TCanvas *can1 = new TCanvas("can1","egamma channel",600,600);
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
	TLegend *leg1 = new TLegend(0.4,0.6,0.87,0.9);
	leg1->SetTextSize(0.025);
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	std::cout << "egamma total integral: " << egamma_phoEt_total->Integral() << std::endl;
	leg1->AddEntry(egamma_phoEt_total, "mixed W#gamma (e#gamma)");
	leg1->AddEntry(egamma_phoEt_1,"WGToLNuG");
	leg1->AddEntry(egamma_phoEt_2,"WGJets_MonoPhoton_PtG-40to130");
	leg1->AddEntry(egamma_phoEt_3,"WGJets_MonoPhoton_PtG-130");
	leg1->Draw("same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( can1,1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( can1,2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( can1,3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( can1,4,1, 11 );

	can1->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_egamma_%d%s.pdf",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	can1->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_egamma_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));

	// =======================================================
	// Canvas 2: mugamma channel
	// =======================================================
	TCanvas *can2 = new TCanvas("can2","mugamma channel",600,600);
	can2->cd();
	gPad->SetLogy();
	mugamma_phoEt_total->SetFillColor(17);
	mugamma_phoEt_total->SetLineColor(17);
	mugamma_phoEt_total->SetFillStyle(1001);
	mugamma_phoEt_total->Draw("hist");
	mugamma_phoEt_1->SetMarkerStyle(7);
	mugamma_phoEt_1->SetMarkerColor(kBlue);
	mugamma_phoEt_1->SetLineColor(kBlue);
	mugamma_phoEt_1->Draw("EP same");
	mugamma_phoEt_2->SetMarkerStyle(7);
	mugamma_phoEt_2->SetMarkerColor(kRed);
	mugamma_phoEt_2->SetLineColor(kRed);
	mugamma_phoEt_2->Draw("EP same");
	mugamma_phoEt_3->SetMarkerStyle(7);
	mugamma_phoEt_3->SetMarkerColor(kGreen);
	mugamma_phoEt_3->SetLineColor(kGreen);
	mugamma_phoEt_3->Draw("EP same");
	TLegend *leg2 = new TLegend(0.4,0.6,0.87,0.9);
	leg2->SetTextSize(0.025);
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	leg2->SetFillColor(0);
	std::cout << "mugamma total integral: " << mugamma_phoEt_total->Integral() << std::endl;
	leg2->AddEntry(mugamma_phoEt_total, "mixed W#gamma (#mu#gamma)");
	leg2->AddEntry(mugamma_phoEt_1,"WGToLNuG");
	leg2->AddEntry(mugamma_phoEt_2,"WGJets_MonoPhoton_PtG-40to130");
	leg2->AddEntry(mugamma_phoEt_3,"WGJets_MonoPhoton_PtG-130");
	leg2->Draw("same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( can2,1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( can2,2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( can2,3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( can2,4,1, 11 );

	can2->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_mugamma_%d%s.pdf",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	can2->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/WGMixing_mugamma_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
}

int main(int argc, char** argv)
{
    bool preVFP = (atoi(argv[2]) == 1);
    plot_Mixing(atoi(argv[1]), preVFP);
    return 0;
}
