#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "../../include/tdrstyle.C"
#include "../../include/analysis_rawData.h"
#include "TROOT.h"

void plot_PUReweight(int RunYear, const char *whichVFP)
{
    setTDRStyle();
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gROOT->SetBatch(kTRUE);

    const int NBINS = 120;
    const double XMIN = 0;
    const double XMAX = 120;
    
    bool preVFP =1;
    if(whichVFP == "preVFP")	preVFP == 1;
    else if(whichVFP == "postVFP")     preVFP == 0;

    TH1D *p_PU_data   = new TH1D("p_PU_data",";N_{vtx};",NBINS,XMIN,XMAX);
    TH1D *p_PU_MC     = new TH1D("p_PU_MC","",NBINS,XMIN,XMAX);
    TH1D *p_PU_MC_raw = new TH1D("p_PU_MC_raw","",NBINS,XMIN,XMAX);

    TChain *sigtree = new TChain("signalTree");
    sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP));
    sigtree->Draw("nVertex >> p_PU_data");

    	TFile *fPU = nullptr;
        fPU = TFile::Open(Form("PUfiles/PUWeights_%d%s.root",RunYear, whichVFP));
    	TH1D *hPUweight = (TH1D*)fPU->Get("PUweight");

        TChain *mctree;
        mctree = new TChain("mgTree","mgTree");
        mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d%s.root",RunYear,whichVFP));

        int   nVertex(0);
        float PUweight(1);
        
	mctree->SetBranchStatus("*", 0);
        mctree->SetBranchStatus("nVertex", 1);
        mctree->SetBranchAddress("nVertex", &nVertex);
    	
	cout << "Total events : " << mctree->GetEntries() << std::endl;
    	for (Long64_t ievt = 0; ievt < mctree->GetEntries(); ++ievt) {
	        if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
                mctree->GetEntry(ievt);
		if(nVertex == 0) continue;
        	int bin = hPUweight->FindBin(nVertex);
        	PUweight = hPUweight->GetBinContent(bin);
        	p_PU_MC->Fill(nVertex, PUweight);
        	p_PU_MC_raw->Fill(nVertex);
		//if (nVertex > 0)  cout << nVertex << "\t" << PUweight << "\n"; 
    	}

    p_PU_data->Sumw2();
    p_PU_MC->Sumw2();
    p_PU_MC_raw->Sumw2();

    p_PU_data->Scale(1.0 / p_PU_data->Integral(1,120));
    p_PU_MC->Scale(1.0 / p_PU_MC->Integral(1,120));
    p_PU_MC_raw->Scale(1.0 / p_PU_MC_raw->Integral(1,120));

    std::cout << "p_PU_data entries = " << p_PU_data->GetEntries() << std::endl;
    std::cout << "p_PU_MC entries   = " << p_PU_MC->GetEntries() << std::endl;
    std::cout << "p_PU_MC_raw entries = " << p_PU_MC_raw->GetEntries() << std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c_PU = new TCanvas("PU", "PU", 600, 600);

    TPad *PU_pad1 = new TPad("PU_pad1","PU_pad1",0,0.3,1,0.98);
    PU_pad1->SetBottomMargin(0);
    PU_pad1->Draw();
    PU_pad1->cd();

    p_PU_data->SetMaximum(1.4 * p_PU_data->GetMaximum());
    p_PU_data->SetMarkerStyle(20);
    p_PU_data->SetLineColor(kBlack);
    p_PU_data->Draw("P");

    p_PU_MC_raw->SetLineColor(kBlue);
    p_PU_MC_raw->Draw("hist same");

    p_PU_MC->SetLineColor(kMagenta);
    p_PU_MC->Draw("hist same");

    TLegend *leg_PU = new TLegend(0.45,0.7,0.8,0.9);
    leg_PU->AddEntry(p_PU_data,"Data","pl");
    leg_PU->AddEntry(p_PU_MC_raw,"Simulation (no correction)","l");
    leg_PU->AddEntry(p_PU_MC,"Simulation (PU corrected)","l");
    leg_PU->Draw();

    if (RunYear == 2016 && preVFP == 1)       CMS_lumi(PU_pad1,1,1,11);
    else if (RunYear == 2016 && preVFP == 0)  CMS_lumi(PU_pad1,2,1,11);
    else if (RunYear == 2017)                 CMS_lumi(PU_pad1,3,1,11);
    else if (RunYear == 2018)                 CMS_lumi(PU_pad1,4,1,11);

    c_PU->cd();
    TPad *PU_pad2 = new TPad("PU_pad2","PU_pad2",0,0.0,1,0.3);
    PU_pad2->SetTopMargin(0);
    PU_pad2->SetBottomMargin(0.4);
    PU_pad2->Draw();
    PU_pad2->cd();

    TH1D *ratio_PU = (TH1D*)p_PU_data->Clone("ratio_PU");
    ratio_PU->Divide(p_PU_MC);
    ratio_PU->SetMarkerStyle(20);
    ratio_PU->SetLineColor(kBlack);
    ratio_PU->GetYaxis()->SetTitle("Data / MC");
    ratio_PU->GetYaxis()->SetRangeUser(0,2.7);
    ratio_PU->Draw("P");

    TLine *line = new TLine(0,1,120,1);
    line->Draw("same");

    c_PU->cd();
    c_PU->Update();
    c_PU->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_nVtx_reweighted_%d%s.png", RunYear, whichVFP));

    /***********************************************************/

    	TCanvas *c_PUw = new TCanvas("c_PUw","PU weight",600,600);
	c_PUw->cd();
	TH1D *hFrame = new TH1D("hFrame", "", 60, 0, 60);
	hFrame->SetMinimum(0.0);
	hFrame->SetMaximum(3.5);
	hFrame->GetXaxis()->SetTitle("N_{vtx}");
	hFrame->GetYaxis()->SetTitle("PU weight");
	hFrame->Draw("axis");

    hPUweight->SetLineWidth(2);
    hPUweight->SetLineColor(kRed+1);
    hPUweight->SetTitle(";N_{vtx};PU weight");
    hPUweight->Draw("hist same");

    if (RunYear == 2016 && preVFP)       CMS_lumi(c_PUw,1,1,11);
    else if (RunYear == 2016 && !preVFP) CMS_lumi(c_PUw,2,1,11);
    else if (RunYear == 2017)            CMS_lumi(c_PUw,3,1,11);
    else if (RunYear == 2018)            CMS_lumi(c_PUw,4,1,11);

    c_PUw->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_PUweight_%d%s.png",RunYear,whichVFP));
    cout << "PU reweighting plot produced successfully" << endl;
}
int main(int argc, char** argv)
{
	plot_PUReweight(atoi(argv[1]),argv[2]);
	return 0;
}
