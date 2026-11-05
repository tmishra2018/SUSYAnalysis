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
#include "TLatex.h"
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


void validateJetFake(){//main 

  double Etreweight[]={9.8084,9.3494,8.69902,7.79226,6.49191,5.19969,4.74666,3.51645,3.07053,2.61752,1.93152,1.5401,1.409,1.09486,0.909713,0.936577,0.864131,0.564228,0.602222,0.499165,0.494812,0.481037,0.413251,0.352619,0.298637,0.257747,0.247043,0.292624,0.272811,0.241988,0.231593,0.186698,0.239414,0.283791,0.179729,0.212801,0.219288,0.172258,0.114328,0.170062,0.154112,0.0960584,0.15714,0.172087,0.134208,0.120842,0.0870008,0.16059,0.104379,0.164007,0.146831,0.129771,0.146885,0.153311,0.105834,0.096875,0.101025,0.089973,0.0933016,0.130326,0.134419,0.122936,0.127599,0.15147,0.0790691,0.0614831,0.10678,0.0878018,0.091917,0.0959308,0.0990909};

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  TH1F *p_ChIso_Z[10];
  TH1F *p_ChIso_mc[10];
  TH1F *p_Sigma_Z[10];
  TH1F *p_Sigma_mc[10];
	TH1F *p_iterative_correction = new TH1F("p_iterative_correction","True Photon #sigma_{i#etai#eta};#sigma_{i#etai#eta};Event",140,0.006,0.02);
	TH1F *p_iterative_signal = new TH1F("p_iterative_signal","Sigmaietaieta;#sigma_{i#etai#eta};Event",140,0.006,0.02);
	TH1F *p_iterative_sideband = new TH1F("p_iterative_sideband","Sigmaietaieta;#sigma_{i#etai#eta};Event",140,0.006,0.02);
  std::ostringstream histname;
  for(unsigned i(0); i<10; i++){
    histname.str("");
    histname << "p_ChIso_Z_" << i;
    p_ChIso_Z[i]= new TH1F(histname.str().c_str(), histname.str().c_str(), 10,0,10);
    p_ChIso_Z[i]->SetLineColor(kBlue);
    histname.str("");
    histname << "p_ChIso_mc_" << i;
    p_ChIso_mc[i]= new TH1F(histname.str().c_str(), histname.str().c_str(), 10,0,10);
    p_ChIso_mc[i]->SetLineColor(i+1);
    histname.str("");
    histname << "p_Sigma_Z_" << i;
    p_Sigma_Z[i]= new TH1F(histname.str().c_str(), histname.str().c_str(), 10,0,0.02);
    p_Sigma_Z[i]->SetLineColor(kBlue);
    histname.str("");
    histname << "p_Sigma_mc_" << i;
    p_Sigma_mc[i]= new TH1F(histname.str().c_str(), histname.str().c_str(), 10,0,0.02);
    p_Sigma_mc[i]->SetLineColor(i+1);
  }  
 
  TChain *egtree = new TChain("egTree","egTree");
  egtree->Add("/uscms_data/d3/mengleis/data/plot_hadron_GJet.root");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_phoSigma(0);
  float eg_phoChIso(0);
  float eg_sigMT(0);
  float eg_sigMET(0);
  float eg_sigMETPhi(0);
  float eg_dPhiLepMET(0);
  int   eg_nVertex(0);
  float eg_HT(0);
  float eg_nJet(0);
  std::vector<int> *eg_mcPID=0;
  std::vector<float> *eg_mcEta=0;
  std::vector<float> *eg_mcPhi=0;
  std::vector<float> *eg_mcPt=0;
  std::vector<int> *eg_mcMomPID=0;
  
  egtree->SetBranchAddress("phoEt",     &eg_phoEt);
  egtree->SetBranchAddress("phoEta",    &eg_phoEta);
  egtree->SetBranchAddress("phoPhi",    &eg_phoPhi);
  egtree->SetBranchAddress("phoSigma",  &eg_phoSigma);
  egtree->SetBranchAddress("phoChIso",  &eg_phoChIso);
  egtree->SetBranchAddress("sigMT",     &eg_sigMT);
  egtree->SetBranchAddress("sigMET",    &eg_sigMET);
  egtree->SetBranchAddress("sigMETPhi", &eg_sigMETPhi);
  egtree->SetBranchAddress("dPhiLepMET",&eg_dPhiLepMET);
  egtree->SetBranchAddress("nVertex",   &eg_nVertex);
  egtree->SetBranchAddress("HT",        &eg_HT);
  egtree->SetBranchAddress("nJet",      &eg_nJet);
  egtree->SetBranchAddress("mcPID",     &eg_mcPID);
  egtree->SetBranchAddress("mcEta",     &eg_mcEta);
  egtree->SetBranchAddress("mcPhi",     &eg_mcPhi);
  egtree->SetBranchAddress("mcPt",      &eg_mcPt);
  egtree->SetBranchAddress("mcMomPID",  &eg_mcMomPID);

  TChain *Zeetree = new TChain("ZeeTree");
  Zeetree->Add("/uscms_data/d3/mengleis/Rereco/plotZ_DoubleMu_2016Rereco_eg.root");
  float Z_phoEt(0);
  float Z_phoEta(0);
  float Z_phoPhi(0);
  float Z_phoSigma(0);
  float Z_invmass(0);
  float Z_phoChIso(0);
  bool  Z_ismg(0); 
  int   Z_nVtx(0); 

  Zeetree->SetBranchAddress("phoEt",     &Z_phoEt);
  Zeetree->SetBranchAddress("phoEta",    &Z_phoEta);
  Zeetree->SetBranchAddress("phoPhi",    &Z_phoPhi);
  Zeetree->SetBranchAddress("phoSigma",  &Z_phoSigma);
  Zeetree->SetBranchAddress("invmass",   &Z_invmass);
  Zeetree->SetBranchAddress("phoChIso",  &Z_phoChIso);
  Zeetree->SetBranchAddress("ismg",      &Z_ismg);
  Zeetree->SetBranchAddress("nVtx",      &Z_nVtx);

  for(unsigned ievt(0); ievt < Zeetree->GetEntries(); ievt++){
	Zeetree->GetEntry(ievt);
	if(Z_invmass > 70 && Z_invmass < 110){
      for(int iE(0); iE < 10; iE++){
        if(Z_phoEt > 40 + iE*10 && Z_phoEt < 70 + iE*10 && Z_phoSigma < 0.0106)p_ChIso_Z[iE]->Fill(Z_phoChIso);
        if(Z_phoEt > 40 + iE*10 && Z_phoEt < 50 + iE*10 && Z_phoChIso < 1.694)p_Sigma_Z[iE]->Fill(Z_phoSigma);
	  }
    }
  }

  int nEvts = egtree->GetEntries();
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    egtree->GetEntry(ievt);
	if(ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;

    double wgt = 1;
    //if(eg_phoEt > 30 && eg_phoEt < 100){
    //  for(int iE(0); iE<70; iE++)
    //    if(eg_phoEt > 30 + iE && eg_phoEt < 31+ iE)wgt = Etreweight[iE];
    //}
    //else wgt = 1;

	bool isFake(true);
	unsigned mcIndex(0);
    unsigned mcPhoIndex(0);
	float mindR(0.7);
    float phodR(0.7);
	bool hasMatch(false);
	for(unsigned ii(0); ii < eg_mcPID->size(); ii++){
	  float dR = DeltaR(eg_phoEta, eg_phoPhi, (*eg_mcEta)[ii], (*eg_mcPhi)[ii]);
	  float dE = fabs(eg_phoEt - (*eg_mcPt)[ii])/eg_phoEt;
      if((*eg_mcPID)[ii] ==22){
        phodR = DeltaR(eg_phoEta, eg_phoPhi, (*eg_mcEta)[ii], (*eg_mcPhi)[ii]);
        mcPhoIndex = ii;
      }
	  if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;} 
	}
    if(phodR < 0.1)mcIndex = mcPhoIndex;
    if(hasMatch)isFake = isHad(fabs((*eg_mcPID)[mcIndex]), fabs((*eg_mcMomPID)[mcIndex]));

    for(int iE(0); iE < 10; iE++){
      if(eg_phoEt > 40 + iE*10 && eg_phoEt < 70 + iE*10){
    	if(eg_phoSigma < 0.0106 && !isFake)p_ChIso_mc[iE]->Fill(eg_phoChIso);
    	if(eg_phoChIso < 1.694 && !isFake)p_Sigma_mc[iE]->Fill(eg_phoSigma);
      }
    }
		if(!isFake){
			p_iterative_correction->Fill(eg_phoSigma);
			if(eg_phoSigma < 0.0106)p_iterative_signal->Fill(eg_phoSigma);
			else if(eg_phoSigma > 0.0106 && eg_phoSigma < 0.016)p_iterative_sideband->Fill(eg_phoSigma);
		}
  }//loop on  events
	TCanvas *cansb= new TCanvas("sbcan","",600,600);
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	p_iterative_correction->GetXaxis()->SetTitleSize(0.1);
	p_iterative_correction->Draw();
	p_iterative_signal->SetFillColor(kBlue);
	p_iterative_sideband->SetFillColor(kYellow);
	p_iterative_signal->Draw("same");
	p_iterative_sideband->Draw("same");	
	TLatex* latex = new TLatex();
	latex->DrawLatex(0.007,100,"n_{Signal}"); 
  latex->DrawLatex(0.011,100,"n_{Sideband}");
	cansb->SaveAs("IterateCorrection_JetFakePho.pdf");

  TCanvas *can1=new TCanvas("can1","",800,800);
  TH1F *ratio[10];
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.1);
  pad1->Draw();  
  pad1->cd();  
  gPad->SetLogy();
  p_ChIso_Z[0]->Draw();
  p_ChIso_mc[0]->Scale(p_ChIso_Z[0]->GetBinContent(1)/p_ChIso_mc[0]->GetBinContent(1));
  p_ChIso_mc[0]->Draw("same");
  can1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
  pad2->Draw();
  pad2->cd();
	ratio[0]=(TH1F*)p_ChIso_mc[0]->Clone("transfer factor");
	ratio[0]->SetLineColor(2);
	ratio[0]->Divide(p_ChIso_Z[0]);
	ratio[0]->Draw("same");

  TCanvas *can2=new TCanvas("can2","",800,800);
  can2->cd();
  TH1F *ratio2[10];
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
  pad3->SetBottomMargin(0.1);
  pad3->Draw();  
  pad3->cd();  
  gPad->SetLogy();
  p_Sigma_mc[0]->Draw();
  for(unsigned iS(1); iS < 9; iS++){
    p_Sigma_mc[iS]->Scale(p_Sigma_mc[0]->GetBinContent(5)/p_Sigma_mc[iS]->GetBinContent(5));
    p_Sigma_mc[iS]->Draw("same");
  }
  can2->cd();
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.25);
  pad4->Draw();
  pad4->cd();
  for(unsigned iS(0); iS < 9; iS++){
    ratio2[iS]=(TH1F*)p_Sigma_mc[iS]->Clone("transfer factor");
    if(iS == 0)ratio2[0]->Draw();
    ratio2[iS]->SetLineColor(iS+1);
    ratio2[iS]->Divide(p_Sigma_mc[0]);
    ratio2[iS]->Draw("same");
  }
}


