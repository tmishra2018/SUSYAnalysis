// g++ `root-config --cflags` plotMuonTrigger.C -o plotMuonTrigger.exe `root-config --libs`
#include <TROOT.h>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "../include/analysis_commoncode.h"

using namespace std;
bool plotNum = 1;
bool doData = 1;

void MuonTrigger_ExtraPlots (int RunYear, bool preVFP){
	gROOT->SetBatch(kTRUE);
  	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("4.4f");
	setTDRStyle();   
  	Double_t xaxis2d[] = {35,40,45,60,200};
  	Double_t yaxis2d[] = {20,38,45,60,200};
	TH2F 			 *p_crosseff = new TH2F("p_mgEff","MuonEG trigger efficiency; photon p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,yaxis2d); 
	TH1F *p_dimuon = new TH1F("p_dimuon","di-muon invmass; #mu#mu mass(GeV);",200,0,200);
        TH1F *p_invmass = new TH1F("p_invmass","mu#mu#gamma invmass; #mu#mu#gamma mass(GeV);",60,60,120);
	TH1D *p_PhoEt = new TH1D("p_PhoEt"," ;#gamma p_{T} (GeV)",nBkgEtBins,bkgEtBins);
        TH1D *p_LepPt = new TH1D("p_LepPt"," ;#mu p_{T} (GeV)",nBkgPtBins,bkgPtBins);

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	std::string isData;
	if(doData==1) isData = "data";
	else isData = "DY";
	// Ztree from data
	TChain *Ztree = new TChain("ZTree","ZTree");
	Ztree->Add(Form("/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_Data_%d%s.root",RunYear,whichVFP.c_str()));
	float Z_phoEt(0);
	float Z_muPt(0);
	float Z_phoEta(0);
	float Z_dR(0);
	int   Z_phofireHLT;
	int   Z_mufireHLT;
	int   Z_phofireHLT2;
	int   Z_mufireHLT2;
	float Z_dimuonmass;
	float Z_InvMass;

	Ztree->SetBranchAddress("Z_dimuonmass", &Z_dimuonmass);
        Ztree->SetBranchAddress("Z_InvMass", &Z_InvMass);


	Ztree->SetBranchAddress("phoEt",     &Z_phoEt);
	Ztree->SetBranchAddress("phoEta",    &Z_phoEta);
	Ztree->SetBranchAddress("muPt",      &Z_muPt);
	Ztree->SetBranchAddress("dR",        &Z_dR);
	Ztree->SetBranchAddress("phofireHLT",&Z_phofireHLT);
	Ztree->SetBranchAddress("mufireHLT", &Z_mufireHLT);
	Ztree->SetBranchAddress("phofireHLT2",&Z_phofireHLT2);
	Ztree->SetBranchAddress("mufireHLT2", &Z_mufireHLT2);



	// Data tree
if(doData){
	for(unsigned ievt(0); ievt < Ztree->GetEntries(); ievt++){
		Ztree->GetEntry(ievt);
		if(fabs(Z_phoEta) > 1.4442)continue;
		if(Z_phoEt > 200)Z_phoEt = 199;
                if(Z_muPt > 200)Z_muPt = 199;

		if (plotNum==true){
			if((Z_phofireHLT > 0 && Z_mufireHLT > 0) || (Z_phofireHLT2 > 0 && Z_mufireHLT2 > 0)){
				p_crosseff->Fill(Z_phoEt, Z_muPt);
				p_dimuon->Fill(Z_dimuonmass);
				p_invmass->Fill(Z_InvMass);
				p_PhoEt->Fill(Z_phoEt);
				p_LepPt->Fill(Z_muPt);
			}
		}
		else {
			p_crosseff->Fill(Z_phoEt, Z_muPt);
			p_dimuon->Fill(Z_dimuonmass);
			p_invmass->Fill(Z_InvMass);
			p_PhoEt->Fill(Z_phoEt);
			p_LepPt->Fill(Z_muPt);

		}
	}
}



	TChain *DYtree = new TChain("mgTree","mgTree");
        DYtree->Add(Form("/eos/uscms/store/user/tmishra/Trigger/files/plot_MuonTrigger_DY_%d%s.root",RunYear,whichVFP.c_str()));
        float DY_phoEt(0);
        float DY_phoEta(0);
        float DY_muPt(0);
        float DY_muMiniIso(0);
        int   DY_phofireHLT;
        int   DY_mufireHLT;
        int   DY_phofireHLT2;
        int   DY_mufireHLT2;

        DYtree->SetBranchAddress("phoEt",     &DY_phoEt);
        DYtree->SetBranchAddress("phoEta",    &DY_phoEta);
        DYtree->SetBranchAddress("muPt",      &DY_muPt);
        DYtree->SetBranchAddress("muMiniIso", &DY_muMiniIso);
        DYtree->SetBranchAddress("phofireHLT",&DY_phofireHLT);
        DYtree->SetBranchAddress("mufireHLT", &DY_mufireHLT);
        DYtree->SetBranchAddress("phofireHLT2",&DY_phofireHLT2);
        DYtree->SetBranchAddress("mufireHLT2", &DY_mufireHLT2);

        // MC tree
	
if(!doData){
        for(unsigned ievt(0); ievt < DYtree->GetEntries(); ievt++){
                DYtree->GetEntry(ievt);
                if(DY_phoEt > 200)DY_phoEt = 199;
                if(DY_muPt > 200)DY_muPt = 199;
			
		if (plotNum==true){
                        if((DY_phofireHLT > 0 && DY_mufireHLT > 0) || (DY_phofireHLT2 > 0 && DY_mufireHLT2 > 0)){
                                p_crosseff->Fill(DY_phoEt, DY_muPt);}
                }
                else p_crosseff->Fill(DY_phoEt, DY_muPt);
                
        }
}


if(doData){


	TCanvas *can_phoEt = new TCanvas("can_phoEt","",600,600);
	can_phoEt->cd();
	can_phoEt->SetLogy();
	p_PhoEt->Draw();
	 if(plotNum==true) can_phoEt->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_phoEt_Numerator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));
	 else can_phoEt->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_phoEt_Denominator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));

	TCanvas *can_lepPt = new TCanvas("can_lepPt","",600,600);
	can_lepPt->cd();
	can_lepPt->SetLogy();
	p_LepPt->Draw();
	 if(plotNum==true) can_lepPt->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_lepPt_Numerator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));
	 else can_lepPt->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_lepPt_Denominator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));

	TCanvas *can_dimuon = new TCanvas("can_dimuon","",600,600);
	can_dimuon->cd();
	p_dimuon->Draw();
	 if(plotNum==true) can_dimuon->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_dimuon_Numerator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));
	 else can_dimuon->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_dimuon_Denominator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));

	TCanvas *can_invmass = new TCanvas("can_invmass","",600,600);
	can_invmass->cd();
	p_invmass->Draw();
	 if(plotNum==true) can_invmass->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_invmass_Numerator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));
	 else can_invmass->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_invmass_Denominator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));

}

	TCanvas *canZ = new TCanvas("canZ","",600,600);
	canZ->cd();
	canZ->SetRightMargin(0.15);
	gPad->SetLogx();
	gPad->SetLogy();
	gStyle->SetPaintTextFormat("4.2f");
	gStyle->SetTextSize(0.002);
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	p_crosseff->SetMarkerSize(0.7); // Reduce text size for text painting

	//p_crosseff->GetXaxis()->SetTitle("#gamma p_{T} (GeV)");
	// p_crosseff is from data
	p_crosseff->Draw("E colz text");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi(canZ, 1, 2, 11);
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(canZ, 2, 2, 11);
        else if(RunYear==2017)                  CMS_lumi(canZ, 3, 2, 11);
        else if(RunYear==2018)                  CMS_lumi(canZ, 4, 2, 11);

	if(plotNum==true) canZ->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_Statistics_Numerator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));
	else canZ->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_Statistics_Denominator_%s_%d%s.pdf",isData.c_str(),RunYear,whichVFP.c_str()));

}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    bool preVFP = (atoi(argv[2]) == 1);
    MuonTrigger_ExtraPlots (atoi(argv[1]), preVFP);
    return 0;
}
