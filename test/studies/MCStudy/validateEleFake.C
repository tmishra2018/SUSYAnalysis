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


void validateEleFake(){//main 

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  std::ostringstream histname;
 
  TChain *egtree = new TChain("signalTree","signalTree");
  egtree->Add("/uscms_data/d3/mengleis/plot_hadron_DY2.root");
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

	TChain *eleproxytree = new TChain("eleproxyTree","eleproxyTree");
  eleproxytree->Add("/uscms_data/d3/mengleis/plot_hadron_DY2.root");
	float eleproxy_phoEt(0);
	float eleproxy_phoEta(0);
	float eleproxy_phoSigma(0);
	eleproxytree->SetBranchAddress("phoEt",     &eleproxy_phoEt);
	eleproxytree->SetBranchAddress("phoEta",    &eleproxy_phoEta);
	eleproxytree->SetBranchAddress("phoSigma",  &eleproxy_phoSigma);

	TH1F *p_target = new TH1F("p_target","",100,0,100);
	TH1F *p_fake = new TH1F("p_fake","",100,0,100);
	TH1F *p_proxy = new TH1F("p_proxy","",100,0,100);
	
  int nEvts = egtree->GetEntries();
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    egtree->GetEntry(ievt);
		if(ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(eg_phoSigma > 0.0106)continue;

		p_target->Fill(eg_phoEt);
	bool isFake(true);
	unsigned mcIndex(0);
    unsigned mcPhoIndex(0);
	float mindR(0.3);
    float phodR(0.3);
	bool hasMatch(false);
	for(unsigned ii(0); ii < eg_mcPID->size(); ii++){
	  float dR = DeltaR(eg_phoEta, eg_phoPhi, (*eg_mcEta)[ii], (*eg_mcPhi)[ii]);
	  float dE = fabs(eg_phoEt - (*eg_mcPt)[ii])/eg_phoEt;
      if((*eg_mcPID)[ii] ==22){
        phodR = DeltaR(eg_phoEta, eg_phoPhi, (*eg_mcEta)[ii], (*eg_mcPhi)[ii]);
        mcPhoIndex = ii;
      }
	  if(dR < mindR && dE < 0.1){mindR = dR; mcIndex = ii; hasMatch = true;} 
	}
    if(hasMatch)isFake = isHad(fabs((*eg_mcPID)[mcIndex]), fabs((*eg_mcMomPID)[mcIndex]));

	if(hasMatch && fabs((*eg_mcPID)[mcIndex]) == 11)p_fake->Fill(eg_phoEt);
	//	if(hasMatch)std::cout << " ID=" << fabs((*eg_mcPID)[mcIndex]) << "  Mom = " << fabs((*eg_mcMomPID)[mcIndex]) << " " << std::endl;
	//	else std::cout << "no match" << std::endl;

  }//loop on  events

	int nelefakes(0);
  for (unsigned ievt(0); ievt< eleproxytree->GetEntries(); ++ievt){//loop on entries

    eleproxytree->GetEntry(ievt);
		if(eleproxy_phoSigma > 0.0106)continue;
		float w_ele;
		if(eleproxy_phoEt > 30 && eleproxy_phoEt < 40)w_ele=0.018;
		else if(eleproxy_phoEt > 40 && eleproxy_phoEt < 60)w_ele = 0.015;
		else if(eleproxy_phoEt > 60)w_ele = 0.01;
		p_proxy->Fill(eleproxy_phoEt, w_ele);
		nelefakes += 1;
  }
	
  std::cout << " ele fakes = " << nelefakes*0.015 << std::endl;
  p_target->Draw();
	p_fake->SetLineColor(kBlue);
	p_fake->Draw("same");
	p_proxy->Sumw2();
	p_proxy->SetLineColor(kRed);
	p_proxy->Draw("same");
}


