#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
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
#include "TRandom3.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
void analysis_qcdBkg(){

  gROOT->SetBatch(kTRUE);
  gSystem->Load("../../lib/libAnaClasses.so");

	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);

	TH1D *p_LepPt_sigma = new TH1D("p_LepPt_sigma","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_LepPt_dEta = new TH1D("p_LepPt_dEta","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_LepPt_dPhi = new TH1D("p_LepPt_dPhi","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_LepPt_D0   = new TH1D("p_LepPt_D0","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_LepPt_Dz   = new TH1D("p_LepPt_Dz","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_LepPt_R9   = new TH1D("p_LepPt_R9","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_dPhiEleMET_loose = new TH1D("p_dPhiEleMET_loose","dPhiEleMET",32,0,3.2); 


	TH1D *p_dPhiEleMET_1 = new TH1D("p_dPhiEleMET_1","dPhiEleMET",32,0,3.2); 
	TH1D *p_dPhiEleMET_2 = new TH1D("p_dPhiEleMET_2","dPhiEleMET",32,0,3.2); 
	TH1D *p_dPhiEleMET_3 = new TH1D("p_dPhiEleMET_3","dPhiEleMET",32,0,3.2); 
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	//fakeEtree->Add("/uscms_data/d3/mengleis/Sep1/test_egsignal_GJet.root");
	fakeEtree->Add("/eos/uscms/store/user/tmishra/fakeLep/fakelep_egsignal_GJet.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsMedium(0);
	float fakeLepSigma(0);
	float fakeLepdEta(0);
	float fakeLepdPhi(0);
	float	fakeLepRelIso(0);
	float	fakeLepD0(0);
	float	fakeLepDz(0);
	float	fakeLepHoverE(0);
	float	fakeLepEoverP(0);
	int		fakeLepVeto(0);
	int		fakeLepMissHits(0);
	float fakeLepR9(0);
	unsigned fakeLepIndex(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
	float threeMass(0);
  float HT(0);
  float nJet(0);
	int   nBJet(0); 
 
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsMedium",&fakeLepIsMedium);
	fakeEtree->SetBranchAddress("fakeLepSigma",&fakeLepSigma);
	fakeEtree->SetBranchAddress("fakeLepdEta", &fakeLepdEta);
	fakeEtree->SetBranchAddress("fakeLepdPhi", &fakeLepdPhi);
	fakeEtree->SetBranchAddress("fakeLepRelIso",&fakeLepRelIso);
	fakeEtree->SetBranchAddress("fakeLepD0",  &fakeLepD0);
	fakeEtree->SetBranchAddress("fakeLepDz",  &fakeLepDz);
	fakeEtree->SetBranchAddress("fakeLepHoverE", &fakeLepHoverE);
	fakeEtree->SetBranchAddress("fakeLepEoverP", &fakeLepEoverP);
	fakeEtree->SetBranchAddress("fakeLepVeto",   &fakeLepVeto);
	fakeEtree->SetBranchAddress("fakeLepMissHits", &fakeLepMissHits);
	fakeEtree->SetBranchAddress("fakeLepR9",  &fakeLepR9);
	fakeEtree->SetBranchAddress("fakeLepIndex",&fakeLepIndex);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	fakeEtree->SetBranchAddress("threeMass", &threeMass);
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nJet",      &nJet);
  fakeEtree->SetBranchAddress("nBJet",     &nBJet);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		if(ievt%1000 ==0)std::cout <<"event " << ievt << std::endl;
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 1.5)continue;

//		if(fakeLepVeto != 1)continue;
//		if(fakeLepMissHits > 1)continue;
//		if(fakeLepHoverE > 0.253)continue;
//		if(fakeLepEoverP > 0.134)continue;
//		if(fakeLepdEta > 0.00311 )continue;
//		if(fakeLepdPhi > 0.103   )continue;
//		if(fakeLepD0 > 0.05)continue;
//		if(fakeLepDz > 0.1)continue;
//		if(fakeLepR9 > 0.5)continue;
//		if(fakeLepIndex != 1)continue;
//		if(fakeLepSigma < 0.00998 && fakeLepMiniIso < 0.1)continue;
		if(fakeLepMiniIso > 0.4)continue;
		// mini-isolation < 0.4 for fake electron proxies
		p_PhoEt->Fill(phoEt, 1);
		p_PhoEta->Fill(phoEta, 1);
		p_LepPt->Fill(lepPt, 1);
		p_LepEta->Fill(lepEta, 1);
		p_MET->Fill(sigMET, 1);
		p_Mt->Fill(sigMT, 1);
		p_HT->Fill(HT, 1);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), 1);
		if(lepPt < 50)p_dPhiEleMET_1->Fill(fabs(dPhiLepMET), 1);
		else if(lepPt >= 50 && lepPt < 100)p_dPhiEleMET_2->Fill(fabs(dPhiLepMET), 1);
		else if(lepPt >= 100)p_dPhiEleMET_3->Fill(fabs(dPhiLepMET), 1);
		p_nJet->Fill(nJet, 1);	

		p_nBJet->Fill(nBJet, 1);
	

	}
	p_LepPt->Draw();
	p_LepPt_sigma->Scale(p_LepPt->GetEntries()/p_LepPt_sigma->GetEntries());
	p_LepPt_dEta->Scale(p_LepPt->GetEntries()/p_LepPt_dEta->GetEntries());
	p_LepPt_dPhi->Scale(p_LepPt->GetEntries()/p_LepPt_dPhi->GetEntries());
	p_LepPt_D0->Scale(p_LepPt->GetEntries()/p_LepPt_D0->GetEntries());
	p_LepPt_Dz->Scale(p_LepPt->GetEntries()/p_LepPt_Dz->GetEntries());
	p_LepPt_R9->Scale(p_LepPt->GetEntries()/p_LepPt_R9->GetEntries());
	p_LepPt_sigma->SetLineColor(kRed);
	p_LepPt_sigma->Draw("same");
	p_LepPt_dEta->SetLineColor(kMagenta);
	p_LepPt_dEta->Draw("same");
	p_LepPt_dPhi->SetLineColor(kGreen);
	p_LepPt_dPhi->Draw("same");
	p_LepPt_D0->SetLineColor(kYellow);
  p_LepPt_Dz->SetLineColor(kCyan);
  p_LepPt_R9->SetLineColor(28);
	p_LepPt_D0->Draw("same");
  p_LepPt_Dz->Draw("same");
  p_LepPt_R9->Draw("same");
	

	TCanvas *can = new TCanvas("can","",600,600);
	p_dPhiEleMET_1->SetLineColor(1);
	p_dPhiEleMET_2->SetLineColor(2);
	p_dPhiEleMET_3->SetLineColor(3);
	p_dPhiEleMET_1->SetMarkerColor(1);
	p_dPhiEleMET_2->SetMarkerColor(2);
	p_dPhiEleMET_3->SetMarkerColor(3);
	p_dPhiEleMET_2->Scale(p_dPhiEleMET_1->GetEntries()/p_dPhiEleMET_2->GetEntries());	
	p_dPhiEleMET_3->Scale(p_dPhiEleMET_1->GetEntries()/p_dPhiEleMET_3->GetEntries());	
	p_dPhiEleMET_1->SetMarkerStyle(20);
	p_dPhiEleMET_2->SetMarkerStyle(20);
	p_dPhiEleMET_3->SetMarkerStyle(20);
	p_dPhiEleMET_1->Draw();
	p_dPhiEleMET_2->Draw("EP same");
	p_dPhiEleMET_3->Draw("EP same");

	TFile *outputfile = TFile::Open("/eos/uscms/store/user/tmishra/fakeLep/signalTree_qcd.root","RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_dPhiEleMET_loose->Write();
	p_PU->Write();
	p_nJet->Write();
	p_nBJet->Write();
	outputfile->Write();
	outputfile->Close();
}


