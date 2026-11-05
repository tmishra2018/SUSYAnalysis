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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

void plotProxyLep(){//main 

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
//************ Signal Tree **********************//
  TChain *sigtree = new TChain("signalTree","signalTree");
	sigtree->Add("/uscms_data/d3/mengleis/resTree_egsignal_GJets_EM.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	//float lepRelIso;
	//float miniIso;
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
	int   lepPixel(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float nJet(0);
  std::vector<int> *mcPID =0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
	float crosssection;
	float totalEvent;
	int  iPho;
 
  sigtree->SetBranchAddress("phoEt",     &phoEt);
  sigtree->SetBranchAddress("phoEta",    &phoEta);
  sigtree->SetBranchAddress("phoPhi",    &phoPhi);
  sigtree->SetBranchAddress("lepPt",     &lepPt);
  sigtree->SetBranchAddress("lepEta",    &lepEta);
  sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("lepPixel",  &lepPixel);
	//sigtree->SetBranchAddress("lepRelIso", &lepRelIso);
	//sigtree->SetBranchAddress("miniIso",   &miniIso);
  sigtree->SetBranchAddress("sigMT",     &sigMT);
  sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  sigtree->SetBranchAddress("nVertex",   &nVertex);
  sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  sigtree->SetBranchAddress("nJet",      &nJet);
  sigtree->SetBranchAddress("mcPID",           &mcPID);
  sigtree->SetBranchAddress("mcEta",           &mcEta);
  sigtree->SetBranchAddress("mcPhi",           &mcPhi);
  sigtree->SetBranchAddress("mcPt",            &mcPt);
  sigtree->SetBranchAddress("mcMomPID",        &mcMomPID);
	sigtree->SetBranchAddress("crosssection",    &crosssection);
	sigtree->SetBranchAddress("totalEvent",      &totalEvent);
	sigtree->SetBranchAddress("iPho",            &iPho);

//*********** fake lepton *********************//
  TChain *fakeEtree = new TChain("fakeETree","fakeETree");
	fakeEtree->Add("/uscms_data/d3/mengleis/resTree_egsignal_GJets_EM.root");
  float fakeEphoEt(0);
  float fakeEphoEta(0);
  float fakeEphoPhi(0);
  float fakeElepPt(0);
  float fakeElepEta(0);
  float fakeElepPhi(0);
  float fakeEsigMT(0);
  float fakeEsigMET(0);
  float fakeEsigMETPhi(0);
  float fakeEdPhiLepMET(0);
  int   fakeEnVertex(0);
  float fakeEdRPhoLep(0);
  float fakeEEM(0);
  float fakeEnJet(0);
  std::vector<int> *fakeEmcPID=0;
  std::vector<float> *fakeEmcEta=0;
  std::vector<float> *fakeEmcPhi=0;
  std::vector<float> *fakeEmcPt=0;
  std::vector<int> *fakeEmcMomPID=0;
	float fakeEcrosssection;
	float fakeEtotalEvent;

  
  fakeEtree->SetBranchAddress("phoEt",     &fakeEphoEt);
  fakeEtree->SetBranchAddress("phoEta",    &fakeEphoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &fakeEphoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &fakeElepPt);
  fakeEtree->SetBranchAddress("lepEta",    &fakeElepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &fakeElepPhi);
  fakeEtree->SetBranchAddress("sigMT",     &fakeEsigMT);
  fakeEtree->SetBranchAddress("sigMET",    &fakeEsigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &fakeEsigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&fakeEdPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &fakeEnVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &fakeEdRPhoLep);
  fakeEtree->SetBranchAddress("nJet",      &fakeEnJet);
  fakeEtree->SetBranchAddress("mcPID",           &fakeEmcPID);
  fakeEtree->SetBranchAddress("mcEta",           &fakeEmcEta);
  fakeEtree->SetBranchAddress("mcPhi",           &fakeEmcPhi);
  fakeEtree->SetBranchAddress("mcPt",            &fakeEmcPt);
  fakeEtree->SetBranchAddress("mcMomPID",        &fakeEmcMomPID);
	fakeEtree->SetBranchAddress("crosssection",    &fakeEcrosssection);
	fakeEtree->SetBranchAddress("totalEvent",      &fakeEtotalEvent);

	TH1F *p_dphi = new TH1F("p_dphi","",32,0,3.2);
	TH1F *p_proxy= new TH1F("p_proxy","",32,0,3.2);
	TH1F *p_proxy_data= new TH1F("p_proxy_data","",32,0,3.2);
	
  TH1F *p_phoEt = new TH1F("p_phoEt","",100,0,300);
	TH1F *p_proxyEt=new TH1F("p_proxyEt","",100,0,300);
	TH1F *p_proxyEt_data=new TH1F("p_proxyEt_data","",100,0,300);	

	TH1F *p_met = new TH1F("p_met","",100,0,300);
	TH1F *p_proxymet = new TH1F("p_proxymet","",100,0,300);
	TH1F *p_proxymet_data=new TH1F("p_proxymet_data","",100,0,300);
	//TProfile *p_fakerate = new TProfile("p_fakerate","",100,0,500);

	//TProfile *p_bkgeff_RelIso = new TProfile("p_bkgeff_RelIso","",100,0,0.5);
	//TProfile *p_bkgeff_miniIso = new TProfile("p_bkgeff_miniIso","",100,0,0.5);
	//TProfile *p_bkgeff_RelIsoMedium = new TProfile("p_bkgeff_RelIsoMedium","",100,0,100);
	//TProfile *p_bkgeff_miniIsoMedium = new TProfile("p_bkgeff_miniIsoMedium","",100,0,100);
	//TProfile *p_bkgeff_miniIsoTight = new TProfile("p_bkgeff_miniIsoTight","",100,0,100);
	//TProfile *p_profile_RelIso = new TProfile("p_profile_RelIso","",100,0,100);
	//TProfile *p_profile_miniIso = new TProfile("p_profile_miniIso","",100,0,100);

	sigtree->Draw("phoEt >> p_phoEt");
	fakeEtree->Draw("phoEt >> p_proxyEt");
	
	float weight(0);
	for (unsigned ievt(0); ievt< sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		weight = crosssection*1000*36.8/totalEvent;
			p_dphi->Fill(fabs(dPhiLepMET),weight);
			p_met->Fill(sigMET, weight);
		//	p_phoEt->Fill(phoEt,weight);
	  //if(sigMET < 40 || sigMET > 70)continue;	
	  if(iPho != 0){
	  
		std::cout << std::endl;
		std::cout << "e pt " << lepPt << "  eta " << lepEta << " phi " << lepPhi << " pixel " << lepPixel << std::endl;
		bool isFake(true);
		unsigned mcIndex(0);
    unsigned mcPhoIndex(0);
		float mindR(0.7);
    float phodR(0.7);
		bool hasMatch(false);
		for(unsigned ii(0); ii < mcPID->size(); ii++){
	  	float dR = DeltaR(lepEta, lepPhi, (*mcEta)[ii], (*mcPhi)[ii]);
	  	float dE = fabs(lepPt - (*mcPt)[ii])/lepPt;
      if(fabs((*mcPID)[ii]) ==11){
        phodR = DeltaR(lepEta, lepPhi, (*mcEta)[ii], (*mcPhi)[ii]);
        mcPhoIndex = ii;
      }
	 		if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;}
			if(dR < 0.3)std::cout << "	around e  ID " << fabs((*mcPID)[ii]) << " Mom " << fabs((*mcMomPID)[ii]) << " pt " << (*mcPt)[ii] << " eta " << (*mcEta)[ii] << " phi " <<  (*mcPhi)[ii]  << std::endl; 
		}
		if(hasMatch)std::cout << "ele ID " << fabs((*mcPID)[mcIndex]) << "  Mom " << fabs((*mcMomPID)[mcIndex]) << std::endl;
		else std::cout << "no match" << std::endl;
	//		p_fakerate->Fill(lepPt, 1);
	//		p_fakerate->Fill(lepPt, 0);

  //    for(int iCut(0); iCut < 100; iCut++){
	//			if(lepRelIso > iCut*0.5/100)p_bkgeff_RelIso->Fill(iCut*0.5/100,1);
	//			else p_bkgeff_RelIso->Fill(iCut*0.5/100,0);
	//			if(miniIso > iCut*0.5/100)p_bkgeff_miniIso->Fill(iCut*0.5/100,1);
	//			else p_bkgeff_miniIso->Fill(iCut*0.5/100,0);
	//		}
	//		if(lepRelIso > 0.0695)p_bkgeff_RelIsoMedium->Fill(lepPt, 1);
	//		else p_bkgeff_RelIsoMedium->Fill(lepPt, 0);
	//		if(miniIso > 0.2)p_bkgeff_miniIsoMedium->Fill(lepPt, 1);
	//		else p_bkgeff_miniIsoMedium->Fill(lepPt, 0);
	//		if(miniIso > 0.1)p_bkgeff_miniIsoTight->Fill(lepPt, 1);
	//		else p_bkgeff_miniIsoTight->Fill(lepPt, 0);

	//		p_profile_RelIso->Fill(lepPt, lepRelIso);
	//		p_profile_miniIso->Fill(lepPt, miniIso);

		std::cout << "pho pt " << phoEt << " eta " << phoEta << " phi " << phoPhi << std::endl;
		isFake = true;
		mcIndex = 0;
    mcPhoIndex = 0;
		mindR = 0.7;
    phodR = 0.7;
		hasMatch = false;
		for(unsigned ii(0); ii < mcPID->size(); ii++){
	  	float dR = DeltaR(phoEta, phoPhi, (*mcEta)[ii], (*mcPhi)[ii]);
	  	float dE = fabs(phoEt - (*mcPt)[ii])/phoEt;
      if(fabs((*mcPID)[ii]) ==22){
        phodR = DeltaR(phoEta, phoPhi, (*mcEta)[ii], (*mcPhi)[ii]);
        mcPhoIndex = ii;
      }
	 		if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;} 
			if(dR < 0.3)std::cout << "	around pho  ID " << fabs((*mcPID)[ii]) << " Mom " << fabs((*mcMomPID)[ii]) << " pt " << (*mcPt)[ii] << " eta " << (*mcEta)[ii] << " phi " <<  (*mcPhi)[ii]  << std::endl; 
		}
		isFake = isHad(fabs((*mcPID)[mcIndex]), fabs((*mcMomPID)[mcIndex]));

		if(hasMatch)std::cout << "pho ID " << fabs((*mcPID)[mcIndex]) << "  Mom " << fabs((*mcMomPID)[mcIndex]) << std::endl;
		}
	}
	p_dphi->Sumw2();
	p_phoEt->Sumw2();
	p_met->Sumw2();

	for (unsigned ievt(0); ievt< fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		weight = fakeEcrosssection*1000*36.8/fakeEtotalEvent;
		//if(fakeEsigMET < 40 || fakeEsigMET > 70)continue;
/////		p_proxy->Fill(fakeEdPhiLepMET,weight);
/////		p_proxymet->Fill(fakeEsigMET,weight);
/////		p_proxyEt->Fill(fakeElepPt, weight);
//		std::cout << "fake  et " << fakeEphoEt << "  eta " << fakeEphoEta << std::endl;
//		p_fakerate->Fill(fakeElepPt, 0);
		bool isFake(true);
		unsigned fakeEmcIndex(0);
    unsigned fakeEmcPhoIndex(0);
		float mindR(0.7);
    float phodR(0.7);
		bool hasMatch(false);
		for(unsigned ii(0); ii < fakeEmcPID->size(); ii++){
	  	float dR = DeltaR(fakeEphoEta, fakeEphoPhi, (*fakeEmcEta)[ii], (*fakeEmcPhi)[ii]);
	  	float dE = fabs(fakeEphoEt - (*fakeEmcPt)[ii])/fakeEphoEt;
      if(fabs((*fakeEmcPID)[ii]) ==22){
        phodR = DeltaR(fakeEphoEta, fakeEphoPhi, (*fakeEmcEta)[ii], (*fakeEmcPhi)[ii]);
        fakeEmcPhoIndex = ii;
      }
	 		if(dR < mindR && dE < 0.7){mindR = dR; fakeEmcIndex = ii; hasMatch = true;} 
		}
		std::cout << std::endl;
		if(hasMatch)std::cout << "fake ID " << fabs((*fakeEmcPID)[fakeEmcIndex]) << "  Mom " << fabs((*fakeEmcMomPID)[fakeEmcIndex]) << std::endl;
		else std::cout << "fake no match" << std::endl;


		isFake = true;
		fakeEmcIndex = 0;
    fakeEmcPhoIndex = 0;
		mindR = 0.7;
    phodR = 0.7;
		hasMatch = false;
		for(unsigned ii(0); ii < fakeEmcPID->size(); ii++){
	  	float dR = DeltaR(fakeElepEta, fakeElepPhi, (*fakeEmcEta)[ii], (*fakeEmcPhi)[ii]);
	  	float dE = fabs(fakeElepPt - (*fakeEmcPt)[ii])/fakeElepPt;
      if(fabs((*fakeEmcPID)[ii]) ==11){
        phodR = DeltaR(fakeElepEta, fakeElepPhi, (*fakeEmcEta)[ii], (*fakeEmcPhi)[ii]);
        fakeEmcPhoIndex = ii;
      }
	 		if(dR < mindR && dE < 0.7){mindR = dR; fakeEmcIndex = ii; hasMatch = true;} 
		}
		if(hasMatch)std::cout << "fakeE ID " << fabs((*fakeEmcPID)[fakeEmcIndex]) << "  Mom " << fabs((*fakeEmcMomPID)[fakeEmcIndex]) << std::endl;
		else std::cout << "no match" << std::endl;

		p_proxy->Fill(fakeEdPhiLepMET,weight);
		p_proxymet->Fill(fakeEsigMET,weight);
//		p_proxyEt->Fill(fakeEphoEt, weight);
	}
  p_proxy->Sumw2();
	p_proxyEt->Sumw2();
	p_proxymet->Sumw2();

//*********** fake lepton *********************//
	TFile *datafile = TFile::Open("root://cmseos.fnal.gov//store/user/msun/ReMiniAOD/resTree_egsignal_ReMiniAOD.root");
  TTree *datatree = (TTree*)datafile->Get("fakeETree");
	float datadPhiLepMET(0);
	float dataphoEt(0);
	float dataMET(0);
	datatree->SetBranchAddress("phoEt",     &dataphoEt);
  datatree->SetBranchAddress("sigMET",    &dataMET);
  datatree->SetBranchAddress("dPhiLepMET",&datadPhiLepMET);
	for (unsigned ievt(0); ievt< datatree->GetEntries(); ievt++){		
		datatree->GetEntry(ievt);
		//if(dataMET < 40 || dataMET > 70)continue;
		p_proxy_data->Fill(datadPhiLepMET);
		p_proxyEt_data->Fill(dataphoEt);
		p_proxymet_data->Fill(dataMET);
	}		


	gStyle->SetOptStat(0);

//	TFile *outputfile = TFile::Open("plot_eg_GJet.root","RECREATE");
//	outputfile->cd();
//	p_dphi->Write();
//	p_bkgeff_RelIsoMedium->Write();
//	p_bkgeff_miniIsoMedium->Write();
//	p_bkgeff_miniIsoTight->Write();
//	outputfile->Write();
//	outputfile->Close();
	TCanvas *can1=new TCanvas("can1","",600,600);
	can1->cd();
  p_dphi->Draw();
	p_proxy->SetLineColor(kRed);
	p_proxy->Scale(p_dphi->Integral(1,32)/p_proxy->Integral(1,32));
	p_proxy->Draw("same");
	p_proxy_data->SetLineColor(kGreen);
	p_proxy_data->Scale(p_dphi->Integral(1,32)/p_proxy_data->Integral(1,32));
	p_proxy_data->Draw("same");

	TCanvas *can2=new TCanvas("can2","",600,600);
	can2->cd();
  p_phoEt->Draw();
	p_proxyEt->SetLineColor(kRed);
	p_proxyEt->Scale(p_phoEt->GetEntries()/p_proxyEt->GetEntries());
	p_proxyEt->Draw("same");
	p_proxyEt_data->SetLineColor(kGreen);
	p_proxyEt_data->Scale(p_phoEt->Integral(1,100)/p_proxyEt_data->Integral(1,100));
	p_proxyEt_data->Draw("same");

	TCanvas *can3=new TCanvas("can3","",600,600);
	can3->cd();
	p_met->Draw();
	p_proxymet->SetLineColor(kRed);
	p_proxymet->Scale(p_met->Integral(1,100)/p_proxymet->Integral(1,100));
	p_proxymet->Draw("same");
	p_proxymet_data->SetLineColor(kGreen);
	p_proxymet_data->Scale(p_met->Integral(1,100)/p_proxymet_data->Integral(1,100));
	p_proxymet_data->Draw("same");
	
//	TFile *outputfile = TFile::Open("plot_eg_GJet-600ToInf.root","RECREATE");
//	outputfile->cd();
//	p_dphi->Write();
//	p_proxy->Write();
//	p_proxy_data->Write();
//	p_phoEt->Write();
//	p_proxyEt->Write();	
//	p_proxyEt_data->Write();
//	outputfile->Write();
//	outputfile->Close();
}


