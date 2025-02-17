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
#include "../include/tdrstyle.C"

int RunYear = 2016;
bool preVFP = true;

void plotMuonTrigger(){//main  
	gROOT->SetBatch(kTRUE);
  	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("4.4f");
	setTDRStyle();   
  	Double_t xaxis2d[] = {35,40,45,60,200};
  	Double_t yaxis2d[] = {20,38,45,60,200};
	TProfile2D *p_HLTeff_Z  = new TProfile2D("p_HLTeff_Z", "#mu#gamma trigger efficiency;photon p_{T} (GeV); #mu p_{T} (GeV)",4,xaxis2d,4,yaxis2d);
	TProfile2D *p_HLTeff_mg = new TProfile2D("p_HLTeff_mg","MET dataset",4,xaxis2d,4,yaxis2d);
	TProfile2D *p_HLTeff_DY = new TProfile2D("p_HLTeff_DY","DY efficiency",4,xaxis2d,4,yaxis2d);
	TH2F 			 *p_crosseff = new TH2F("p_mgEff","MuonEG trigger efficiency; photon p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,yaxis2d); 
	TH2F 			 *p_mgESF = new TH2F("p_mgESF","MuonEG trigger ESF; photon p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,yaxis2d); 

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	TChain *mgtree = new TChain("mgTree","mgTree");
	mgtree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger_MET.root");
	float mg_phoEt(0);
	float mg_phoEta(0);
	float mg_muPt(0);
	float mg_muMiniIso(0);
	float mg_dR(0);
	int   mg_phofireHLT;
	int   mg_mufireHLT;
	int   mg_phofireHLT2;
	int   mg_mufireHLT2;
	int   mg_passHLT;
	int   mg_passHLT2; 

	mgtree->SetBranchAddress("phoEt",     &mg_phoEt);
	mgtree->SetBranchAddress("phoEta",    &mg_phoEta);
	mgtree->SetBranchAddress("muPt",      &mg_muPt);
	mgtree->SetBranchAddress("muMiniIso", &mg_muMiniIso);
	mgtree->SetBranchAddress("dR",        &mg_dR);
	mgtree->SetBranchAddress("phofireHLT",&mg_phofireHLT);
	mgtree->SetBranchAddress("mufireHLT", &mg_mufireHLT);
	mgtree->SetBranchAddress("phofireHLT2",&mg_phofireHLT2);
	mgtree->SetBranchAddress("mufireHLT2", &mg_mufireHLT2);
  	mgtree->SetBranchAddress("passHLT",    &mg_passHLT);
  	mgtree->SetBranchAddress("passHLT2",   &mg_passHLT2);

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

	Ztree->SetBranchAddress("phoEt",     &Z_phoEt);
	Ztree->SetBranchAddress("phoEta",    &Z_phoEta);
	Ztree->SetBranchAddress("muPt",      &Z_muPt);
	Ztree->SetBranchAddress("dR",        &Z_dR);
	Ztree->SetBranchAddress("phofireHLT",&Z_phofireHLT);
	Ztree->SetBranchAddress("mufireHLT", &Z_mufireHLT);
	Ztree->SetBranchAddress("phofireHLT2",&Z_phofireHLT2);
	Ztree->SetBranchAddress("mufireHLT2", &Z_mufireHLT2);

	for(unsigned ievt(0); ievt < Ztree->GetEntries(); ievt++){
		Ztree->GetEntry(ievt);
		if(fabs(Z_phoEta) > 1.4442)continue;
		if(Z_phoEt > 200)Z_phoEt = 199;
		if(Z_muPt > 200)Z_muPt = 199;
		// [(probePho->fireDoubleTrg(28) || probePho->fireDoubleTrg(29)) && (probeMu->fireSingleTrg(2)] or [probeMu->fireSingleTrg(21) || (probePho->fireDoubleTrg(30)&&(probeMu->fireSingleTrg(22))]
		// hltMu17Photon30CaloIdLL1ISOHEFilter || hltMu17Photon30CaloIdLL1ISOORHEFilter && hltL3fL1sL1Mu5IsoEG18L1f5L2f7L3Filtered17 || hltL3fL1sL1Mu5IsoEG18ORL1Mu5IsoEG20L1f5L2f7L3Filtered17
		// hltMu38NoFiltersNoVtxPhoton38CaloIdLHEFilter && hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered38
		//     sometthing wrong in Z_mufireHLT or Z_mufireHLT2
//		if((Z_phofireHLT > 0 && Z_mufireHLT > 0) || (Z_phofireHLT2 > 0 && Z_mufireHLT2 > 0))
		if((Z_phofireHLT > 0 || Z_phofireHLT2 > 0) && (Z_mufireHLT > 0 || Z_mufireHLT2 > 0))
			p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 1);
		else p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 0);
	}

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);
		//if(mg_phoEt > 200)mg_phoEt = 199;
    		//if(mg_muPt > 200)mg_muPt = 199;
		if(mg_passHLT > 0 || mg_passHLT2 > 0)p_HLTeff_mg->Fill(mg_phoEt, mg_muPt, 1);
		else p_HLTeff_mg->Fill(mg_phoEt, mg_muPt, 0);
	}

	TCanvas *canmg = new TCanvas("canmg","",600,600);
	canmg->cd();
	p_HLTeff_mg->Draw("colz text");
  	for(int binx(1); binx <= 4; binx++){
		for(int biny(1); biny <= 4; biny++){
			p_crosseff->SetBinContent(binx, biny, p_HLTeff_Z->GetBinContent(binx, biny));
			p_crosseff->SetBinError(binx, biny, fabs(p_HLTeff_mg->GetBinContent(binx, biny) - p_HLTeff_Z->GetBinContent(binx, biny)));
			// error is using MET dataset
		}
	}	

	TCanvas *canZ = new TCanvas("canZ","",600,600);
	canZ->cd();
	canZ->SetRightMargin(0.12);
	gPad->SetLogx();
	gPad->SetLogy();
	gStyle->SetPaintTextFormat("4.2f");
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	//p_crosseff->GetXaxis()->SetTitle("#gamma p_{T} (GeV)");
	// p_crosseff is from data
	p_crosseff->Draw("E colz text");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi(canZ, 1, 2, 11);
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(canZ, 2, 2, 11);
        else if(RunYear==2017)                  CMS_lumi(canZ, 3, 2, 11);
        else if(RunYear==2018)                  CMS_lumi(canZ, 4, 2, 11);
	canZ->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_efficiency_%d%s.pdf",RunYear,whichVFP.c_str()));

	// DYtree from DY mc
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

	for(unsigned ievt(0); ievt < DYtree->GetEntries(); ievt++){
		DYtree->GetEntry(ievt);
		if(DY_phoEt > 200)DY_phoEt = 199;
		if(DY_muPt > 200)DY_muPt = 199;
		if((DY_phofireHLT > 0 || DY_phofireHLT2 > 0) && (DY_mufireHLT > 0 || DY_mufireHLT2 > 0))
	//	if((DY_phofireHLT > 0 && DY_mufireHLT > 0) || (DY_phofireHLT2 > 0 && DY_mufireHLT2 > 0))
			p_HLTeff_DY->Fill(DY_phoEt, DY_muPt, 1);
		else p_HLTeff_DY->Fill(DY_phoEt, DY_muPt, 0);
	}
  	for(int binx(1); binx <= 4; binx++){
		for(int biny(1); biny <= 4; biny++){
			// Data/mc scalefactor
			p_mgESF->SetBinContent(binx, biny, p_crosseff->GetBinContent(binx, biny)/p_HLTeff_DY->GetBinContent(binx, biny));
			p_mgESF->SetBinError(binx, biny, p_crosseff->GetBinError(binx, biny)/p_HLTeff_DY->GetBinContent(binx, biny));
		}
	}	
	
	TCanvas *canDY = new TCanvas("canDY","",600,600);
	canDY->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	gStyle->SetPaintTextFormat("4.4f");
	p_HLTeff_DY->Draw("E colz text");

	TCanvas *canESF = new TCanvas("canESF","",600,600);
	canESF->SetRightMargin(0.14);
	gPad->SetLogx();
	gPad->SetLogy();
	gStyle->SetPaintTextFormat("4.4f");
	p_mgESF->Draw("E colz text");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi(canESF, 1, 11);
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(canESF, 2, 11);
        else if(RunYear==2017)                  CMS_lumi(canESF, 5, 11);
        else if(RunYear==2018)                  CMS_lumi(canESF, 4, 11);

	canESF->SaveAs(Form("/eos/uscms/store/user/tmishra/Trigger/mgTrigger_ESF_%d%s.pdf",RunYear,whichVFP.c_str()));
	
	TFile *outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_%d%s.root",RunYear,whichVFP.c_str()),"RECREATE");
	outputfile->cd();
	p_crosseff->Write();
	p_mgESF->Write();
	outputfile->Write();
	outputfile->Close();
}
