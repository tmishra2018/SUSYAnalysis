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
#include "TPad.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"


void validateDPhi(){//main 

  gSystem->Load("../../../lib/libAnaClasses.so");
  std::ostringstream histname;
 
  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/WGToLNuG_corr.root");
  float MET;
  float METphi;
  float genMET;
  float genMETphi;
  float pfMET_T1TxyPhi;
  float pfMET_T1TxyPt;
  std::vector<float> *elePhi = 0;

  es->SetBranchAddress("pfMET", &MET);
  es->SetBranchAddress("pfMETPhi", &METphi);
  es->SetBranchAddress("genMET",&genMET);
  es->SetBranchAddress("genMETPhi",&genMETphi);
  es->SetBranchAddress("elePhi",  &elePhi);
  es->SetBranchAddress("pfMET_T1TxyPhi",&pfMET_T1TxyPhi);
  es->SetBranchAddress("pfMET_T1TxyPt", &pfMET_T1TxyPt);

//  TChain* es = new TChain("egTree","egTree");
//  es->Add("../../VGamma/src/resTree_VGamma_WG.root");
//  float METPhi(0);
//  float genMETPhi(0);
//  float lepPhi(0);
//  float MET(0);
//  es->SetBranchAddress("lepPhi",    &lepPhi);
//  es->SetBranchAddress("sigMET",    &MET);
//  es->SetBranchAddress("sigMETPhi", &METPhi);
//  es->SetBranchAddress("genMETPhi", &genMETPhi);
  
	TH1F *p_dphi_gen = new TH1F("p_dphi_gen","",64,-3.2,3.2);
	TH1F *p_dphi_pf  = new TH1F("p_dphi_pf","", 64,-3.2,3.2);
	TH1F *p_dphi_corr= new TH1F("p_dphi_corr","",64,-3.2,3.2);
	for (unsigned ievt(0); ievt< 500000; ++ievt){//loop on entries

		if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
		es->GetEntry(ievt);
 //   if(MET > 40){
 // 	  p_dphi_gen->Fill(DeltaPhi(lepPhi, METPhi));
 // 	  p_dphi_pf->Fill(DeltaPhi(lepPhi,  genMETPhi));
 //   }
 //   if(MET < 40){
 // 		p_diffdphi_lowmet->Fill(DeltaPhi(lepPhi, METPhi) - DeltaPhi(lepPhi,  genMETPhi));
 // 	}
 //   else if(MET > 70){
 // 		p_diffdphi_highmet->Fill(DeltaPhi(lepPhi, METPhi) - DeltaPhi(lepPhi,  genMETPhi));
 // 	}
 // 	if(fabs(genMETPhi)<3.1415/2)p_dphi_pf_lowmet->Fill(METPhi);
 // 	else p_dphi_pf_highmet->Fill(METPhi);
  		if(elePhi->size() > 0){
	 			if(pfMET_T1TxyPt > 40)p_dphi_corr->Fill(DeltaPhi((*elePhi)[0], pfMET_T1TxyPhi));			
				if(MET > 40)p_dphi_pf->Fill(DeltaPhi((*elePhi)[0], METphi));
				if(genMET > 40)p_dphi_gen->Fill(DeltaPhi((*elePhi)[0], genMETphi));
			}
	}
 // p_dphi_gen->Draw();
 //p_dphi_pf->Draw("EP same");	
  // p_diffdphi_lowmet->Draw();
  // p_diffdphi_highmet->SetLineColor(kRed);
  // p_diffdphi_highmet->Draw("same");
//   p_dphi_pf_low->Draw();
//	p_dphi_pf_high->Scale(p_dphi_pf_high->GetEntries()/p_dphi_pf_low->GetEntries());
//   p_dphi_pf_high->SetLineColor(kRed);
//   p_dphi_pf_high->Draw("same");
	p_dphi_gen->Scale(p_dphi_pf->GetEntries()/p_dphi_gen->GetEntries());
	p_dphi_gen->Draw();
	p_dphi_pf->SetLineColor(kRed);
	p_dphi_pf->Draw("same");
	p_dphi_corr->SetLineColor(kGreen);
	p_dphi_corr->Draw("same");

	gStyle->SetOptStat(0);
	TLegend *leg= new TLegend(0.5,0.7,0.9,0.9);
	leg->AddEntry(p_dphi_gen, "corrected #Delta#phi, MET > 40");
	leg->AddEntry(p_dphi_gen, "#Delta#phi, MET > 40");
	leg->Draw("same");
}
