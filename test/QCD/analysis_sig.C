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

void analysis_sig(){

  gSystem->Load("../../lib/libAnaClasses.so");
	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",9,0,9);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);

	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	sigtree->Add("/eos/uscms/store/user/tmishra/fakeLep/fakelep_egsignal_GJet.root");
	//sigtree->Add("/uscms_data/d3/mengleis/Sep1/test_egsignal_QCDfake_2.root");

	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float lepMiniIso(0);
	float lepRelIso(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float sigMETPhi(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float threeMass(0);
	float nJet(0);
	int   nBJet(0);	
  std::vector<int> *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
  std::vector<int> *mcGMomPID=0;
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("lepMiniIso", &lepMiniIso);
	sigtree->SetBranchAddress("lepRelIso", &lepRelIso);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("HT",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	sigtree->SetBranchAddress("threeMass", &threeMass);
	sigtree->SetBranchAddress("nJet",      &nJet);
	sigtree->SetBranchAddress("nBJet",     &nBJet);
  sigtree->SetBranchAddress("mcPID",     &mcPID);
  sigtree->SetBranchAddress("mcEta",     &mcEta);
  sigtree->SetBranchAddress("mcPhi",     &mcPhi);
  sigtree->SetBranchAddress("mcPt",      &mcPt);
  sigtree->SetBranchAddress("mcMomPID",  &mcMomPID);
  sigtree->SetBranchAddress("mcGMomPID", &mcGMomPID);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 1.5)continue;

		//if(lepMiniIso > 0.1)continue;
		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		p_nJet->Fill(nJet);
		
		p_nBJet->Fill(nBJet);
	}        

	std::ostringstream outputname;

	outputname << "/eos/uscms/store/user/tmishra/fakeLep/signalTree_sig.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_eventcount->Write();
	p_nJet->Write();
	p_nBJet->Write();
	outputfile->Write();
	outputfile->Close();
}


//		bool isFake(true);
//		double mindRele(0.3);
//		unsigned eleIndex(0);
//		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
//			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], lepEta, lepPhi);
//			double dE1 = fabs((*mcPt)[iMC] - lepPt)/lepPt;
//			if(dR1 < mindRele){mindRele=dR1; eleIndex=iMC;}
//		}
//		bool isTrueEle(false);
//		if(mindRele < 0.1){
//			//std::cout << (*mcPID)[eleIndex] << " " << fabs((*mcMomPID)[eleIndex]) << std::endl;
//			if(fabs((*mcPID)[eleIndex]) == 11 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23 || fabs((*mcMomPID)[eleIndex]) == 15 || fabs((*mcMomPID)[eleIndex]) == 22))isTrueEle = true;
//			else if(fabs((*mcPID)[eleIndex]) == 13 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23 || fabs((*mcMomPID)[eleIndex]) == 15))isTrueEle = true;
//			else if(fabs((*mcPID)[eleIndex]) == 15 && (fabs((*mcMomPID)[eleIndex]) == 24 || fabs((*mcMomPID)[eleIndex]) == 23))isTrueEle = true;
//		}
//		//if(lepPt > 200){
//		if(!isTrueEle){
//			if(mindRele < 0.1){
//				if((*mcPID)[eleIndex] ==22){
//					std::cout << "pho " <<  (*mcPID)[eleIndex] << " " << fabs((*mcMomPID)[eleIndex]) << " " << fabs((*mcGMomPID)[eleIndex]) <<  std::endl;
//					std::cout << "    " <<  "pixel " << nJet << std::endl;
//					for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
//						std::cout << "    " <<  (*mcPID)[iMC] << " " << fabs((*mcMomPID)[iMC]) << " " << fabs((*mcGMomPID)[iMC]) <<  std::endl;	
//					}
//				}
//				else std::cout << "ele " <<  (*mcPID)[eleIndex] << " " << fabs((*mcMomPID)[eleIndex]) << " " << fabs((*mcGMomPID)[eleIndex]) <<  std::endl;
//			}
//			else std::cout << "no match" << std::endl;
//		}
//		else std::cout << "true electron" << std::endl;
//		//}
