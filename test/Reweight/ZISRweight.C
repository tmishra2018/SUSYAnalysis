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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"

void ZISRweight(){//main 

	int channelType = 2;
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");


	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={35,50,100,150,200,250,300,400,600,800};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	//TH1D *p_lepPt_data     = new TH1D("p_lepPt_data","",9,plotEtBins);
	TH1D *p_lepPt_data     = new TH1D("p_lepPt_data","",15,25,400);
	TH1D *p_trailPt_data     = new TH1D("p_trail_data","",9,plotEtBins);
	TH1D *p_lepEta_data    = new TH1D("p_lepEta_data","",60,-3,3);
	TH1D *p_lepPhi_data    = new TH1D("p_lepPhi_data","",64,-3.2, 3.2);
	TH1D *p_sigMET_data    = new TH1D("p_sigMET_data","",100,0,400);
	TH1D *p_sigMT_data     = new TH1D("p_sigMT_data","",100,0,400);
	TH1D *p_sigMETPhi_data = new TH1D("p_sigMETPhi_data","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_data= new TH1D("p_dPhiLepMET_data","",64,-3.2, 3.2);
	TH1D *p_nVertex_data   = new TH1D("p_nVertex_data","",100,0,100);
	TH1D *p_HT_data        = new TH1D("p_HT_data","",100,0,1000);
	TH1D *p_nJet_data      = new TH1D("p_nJet_data","",10,0,10);
	TH1D *p_llmass_data   = new TH1D("p_llmass_data","",160,70,110);
	TH1D *p_JetPt_data     = new TH1D("p_JetPt_data","",9,plotEtBins);

	//TH1D *p_lepPt_ZG     = new TH1D("p_lepPt_ZG","",9,plotEtBins);
	TH1D *p_lepPt_ZG     = new TH1D("p_lepPt_ZG","",15,25,400);
	TH1D *p_trailPt_ZG     = new TH1D("p_trail_ZG","",9,plotEtBins);
	TH1D *p_lepEta_ZG    = new TH1D("p_lepEta_ZG","",60,-3,3);
	TH1D *p_lepPhi_ZG    = new TH1D("p_lepPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_sigMET_ZG    = new TH1D("p_sigMET_ZG","",100,0,400);
	TH1D *p_sigMT_ZG     = new TH1D("p_sigMT_ZG","",100,0,400);
	TH1D *p_sigMETPhi_ZG = new TH1D("p_sigMETPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_ZG= new TH1D("p_dPhiLepMET_ZG","",64,-3.2, 3.2);
	TH1D *p_nVertex_ZG   = new TH1D("p_nVertex_ZG","",100,0,100);
	TH1D *p_HT_ZG        = new TH1D("p_HT_ZG","",100,0,1000);
	TH1D *p_nJet_ZG      = new TH1D("p_nJet_ZG","",10,0,10);
	TH1D *p_llmass_ZG   = new TH1D("p_llmass_ZG","",160,70,110);
	TH1D *p_JetPt_ZG     = new TH1D("p_JetPt_ZG","",9,plotEtBins);

	//TH1D *p_lepPt_rare     = new TH1D("p_lepPt_rare","",9,plotEtBins);
	TH1D *p_lepPt_rare     = new TH1D("p_lepPt_rare","",15,25,400);
	TH1D *p_lepEta_rare    = new TH1D("p_lepEta_rare","",60,-3,3);
	TH1D *p_lepPhi_rare    = new TH1D("p_lepPhi_rare","",64,-3.2, 3.2);
	TH1D *p_sigMET_rare    = new TH1D("p_sigMET_rare","",100,0,400);
	TH1D *p_sigMT_rare     = new TH1D("p_sigMT_rare","",100,0,400);
	TH1D *p_sigMETPhi_rare = new TH1D("p_sigMETPhi_rare","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_rare= new TH1D("p_dPhiLepMET_rare","",64,-3.2, 3.2);
	TH1D *p_nVertex_rare   = new TH1D("p_nVertex_rare","",100,0,100);
	TH1D *p_HT_rare        = new TH1D("p_HT_rare","",100,0,1000);
	TH1D *p_nJet_rare      = new TH1D("p_nJet_rare","",10,0,10);
	TH1D *p_llmass_rare   = new TH1D("p_llmass_rare","",160,70,110);

//************ Signal Tree **********************//
  TChain *tree = new TChain("ZTree");
  tree->Add("/uscms_data/d3/mengleis/test/resTree_ZISR_data.root");
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float trailPt(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
	int   nGoodMu(0);
  float HT(0);
  float nJet(0);
	float dilepMass(0);
	std::vector<float> *JetPt=0;

  tree->SetBranchAddress("lepPt",     &lepPt);
  tree->SetBranchAddress("lepEta",    &lepEta);
  tree->SetBranchAddress("lepPhi",    &lepPhi);
  tree->SetBranchAddress("trailPt",   &trailPt);
  tree->SetBranchAddress("sigMT",     &sigMT);
  tree->SetBranchAddress("sigMET",    &sigMET);
  tree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  tree->SetBranchAddress("nVertex",   &nVertex);
  tree->SetBranchAddress("HT",        &HT);
  tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("nGoodMu",   &nGoodMu);
	tree->SetBranchAddress("dilepMass", &dilepMass);
	tree->SetBranchAddress("JetPt",     &JetPt);

  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);

		if(nGoodMu !=2)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		if(trailPt < 10)continue;
		double weight=1;

		p_lepPt_data->Fill(lepPt, weight);
		p_lepEta_data->Fill(lepEta, weight);
		p_lepPhi_data->Fill(lepPhi, weight);
		p_trailPt_data->Fill(trailPt, weight);
		p_sigMET_data->Fill(sigMET, weight);
		p_sigMT_data->Fill(sigMT, weight);
		p_sigMETPhi_data->Fill(sigMETPhi, weight);
		p_dPhiLepMET_data->Fill(dPhiLepMET, weight);
		p_nVertex_data->Fill(nVertex, weight);
		p_HT_data->Fill(HT, weight);
		p_nJet_data->Fill(nJet, weight);
		p_llmass_data->Fill(dilepMass, weight);
	}//loop on  events




//************ Signal Tree **********************//
  TChain *ZGtree = new TChain("ZTree");
  ZGtree->Add("/uscms_data/d3/mengleis/test/resTree_ZISR_DY.root");
  float ZG_lepPt(0);
  float ZG_lepEta(0);
  float ZG_lepPhi(0);
	float ZG_trailPt(0);
	float ZG_trailEta(0);
	float ZG_trailPhi(0);
  float ZG_sigMT(0);
  float ZG_sigMET(0);
  float ZG_sigMETPhi(0);
  float ZG_dPhiLepMET(0);
  int   ZG_nVertex(0);
	int   ZG_nGoodMu(0);
  float ZG_HT(0);
  float ZG_nJet(0);
	float ZG_dilepMass(0);
	std::vector<float> *ZG_JetPt=0;
  std::vector<int>   *ZG_mcPID=0;
  std::vector<float> *ZG_mcEta=0;
  std::vector<float> *ZG_mcPhi=0;
  std::vector<float> *ZG_mcPt=0;
  std::vector<int>   *ZG_mcMomPID=0;
  std::vector<int>   *ZG_mcGMomPID=0;
	std::vector<int>   *ZG_mcStatus=0;

  ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
	ZGtree->SetBranchAddress("trailPt",   &ZG_trailPt);
	ZGtree->SetBranchAddress("trailEta",  &ZG_trailEta);
	//ZGtree->SetBranchAddress("trailPhi",  &ZG_trailPhi);
  ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  ZGtree->SetBranchAddress("sigMETPhi", &ZG_sigMETPhi);
  ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  ZGtree->SetBranchAddress("nVertex",   &ZG_nVertex);
  ZGtree->SetBranchAddress("HT",        &ZG_HT);
  ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
	ZGtree->SetBranchAddress("nGoodMu",   &ZG_nGoodMu);
	ZGtree->SetBranchAddress("dilepMass", &ZG_dilepMass);
	ZGtree->SetBranchAddress("JetPt",     &ZG_JetPt);
  ZGtree->SetBranchAddress("mcPID",    &ZG_mcPID);
  ZGtree->SetBranchAddress("mcEta",    &ZG_mcEta);
  ZGtree->SetBranchAddress("mcPhi",    &ZG_mcPhi);
  ZGtree->SetBranchAddress("mcPt",     &ZG_mcPt);
  ZGtree->SetBranchAddress("mcMomPID", &ZG_mcMomPID);
  ZGtree->SetBranchAddress("mcGMomPID",&ZG_mcGMomPID);
	ZGtree->SetBranchAddress("mcStatus", &ZG_mcStatus);

  for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		if(ZG_nGoodMu != 2)continue;
		double weight = 35.8*1000*1921.8*3.0/122053259; 
		double scalefactor(0);
		scalefactor = objectESF.getMuonESF(ZG_lepPt,ZG_lepEta)*objectESF.getMuonESF(ZG_trailPt,ZG_trailEta);

		if(ZG_dilepMass < 80 || ZG_dilepMass >100)continue;
		if(ZG_trailPt < 10)continue;
		weight = weight*scalefactor;

		float PUweight = getPUESF(ZG_nVertex);
		weight = weight*PUweight;
 
		p_lepPt_ZG->Fill(ZG_lepPt, weight);
		p_lepEta_ZG->Fill(ZG_lepEta, weight);
		p_lepPhi_ZG->Fill(ZG_lepPhi, weight);
		p_trailPt_ZG->Fill(ZG_trailPt, weight);
		p_sigMET_ZG->Fill(ZG_sigMET, weight);
		p_sigMT_ZG->Fill(ZG_sigMT, weight);
		p_sigMETPhi_ZG->Fill(ZG_sigMETPhi, weight);
		p_dPhiLepMET_ZG->Fill(ZG_dPhiLepMET, weight);
		p_nVertex_ZG->Fill(ZG_nVertex, weight);
		p_HT_ZG->Fill(ZG_HT, weight);
		p_nJet_ZG->Fill(ZG_nJet, weight);
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events

//************ Signal Tree **********************//
  TChain *raretree = new TChain("ZTree");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ZISR_TT.root");
  //raretree->Add("/uscms_data/d3/mengleis/test/resTree_ZGamma_WWG.root");
  //raretree->Add("/uscms_data/d3/mengleis/test/resTree_ZGamma_WZG.root");
	float MCweight(0);
  float rare_lepPt(0);
  float rare_lepEta(0);
  float rare_lepPhi(0);
  float rare_sigMT(0);
  float rare_sigMET(0);
  float rare_sigMETPhi(0);
  float rare_dPhiLepMET(0);
  int   rare_nVertex(0);
  float rare_HT(0);
  float rare_nJet(0);
	float rare_dilepMass(0);
  std::vector<int>   *rare_mcPID=0;
  std::vector<float> *rare_mcEta=0;
  std::vector<float> *rare_mcPhi=0;
  std::vector<float> *rare_mcPt=0;
  std::vector<int>   *rare_mcMomPID=0;
  std::vector<int>   *rare_mcGMomPID=0;

  raretree->SetBranchAddress("MCweight",  &MCweight);
  raretree->SetBranchAddress("lepPt",     &rare_lepPt);
  raretree->SetBranchAddress("lepEta",    &rare_lepEta);
  raretree->SetBranchAddress("lepPhi",    &rare_lepPhi);
  raretree->SetBranchAddress("sigMT",     &rare_sigMT);
  raretree->SetBranchAddress("sigMET",    &rare_sigMET);
  raretree->SetBranchAddress("sigMETPhi", &rare_sigMETPhi);
  raretree->SetBranchAddress("dPhiLepMET",&rare_dPhiLepMET);
  raretree->SetBranchAddress("nVertex",   &rare_nVertex);
  raretree->SetBranchAddress("HT",        &rare_HT);
  raretree->SetBranchAddress("nJet",      &rare_nJet);
	raretree->SetBranchAddress("dilepMass", &rare_dilepMass);
  raretree->SetBranchAddress("mcPID",    &rare_mcPID);
  raretree->SetBranchAddress("mcEta",    &rare_mcEta);
  raretree->SetBranchAddress("mcPhi",    &rare_mcPhi);
  raretree->SetBranchAddress("mcPt",     &rare_mcPt);
  raretree->SetBranchAddress("mcMomPID", &rare_mcMomPID);
  raretree->SetBranchAddress("mcGMomPID",&rare_mcGMomPID);

  for(unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);

		double weight = MCweight; 
		double scalefactor(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(rare_lepPt,rare_lepEta);
		}

		if(rare_dilepMass < 80 || rare_dilepMass > 100)continue;
		weight = weight*scalefactor;
		float PUweight = getPUESF(rare_nVertex);
		weight = weight*PUweight;

		p_lepPt_rare->Fill(rare_lepPt, weight);
		p_lepEta_rare->Fill(rare_lepEta, weight);
		p_lepPhi_rare->Fill(rare_lepPhi, weight);
		p_sigMET_rare->Fill(rare_sigMET, weight);
		p_sigMT_rare->Fill(rare_sigMT, weight);
		p_sigMETPhi_rare->Fill(rare_sigMETPhi, weight);
		p_dPhiLepMET_rare->Fill(rare_dPhiLepMET, weight);
		p_nVertex_rare->Fill(rare_nVertex, weight);
		p_HT_rare->Fill(rare_HT, weight);
		p_nJet_rare->Fill(rare_nJet, weight);
		p_llmass_rare->Fill(rare_dilepMass, weight);
	}//loop on  events
	
	TCanvas *can_lepPt     = new TCanvas("can_lepPt",       "can_lepPt",600,600); 
	TCanvas *can_lepEta    = new TCanvas("can_lepEta",      "can_lepEta",600,600); 
	TCanvas *can_lepPhi    = new TCanvas("can_lepPhi",      "can_lepPhi",600,600); 
	TCanvas *can_sigMET    = new TCanvas("can_sigMET",      "can_sigMET",600,600); 
	TCanvas *can_sigMT     = new TCanvas("can_sigMT",       "can_sigMT",600,600); 
	TCanvas *can_sigMETPhi = new TCanvas("can_sigMETPhi",   "can_sigMETPhi",600,600); 
	TCanvas *can_dPhiLepMET= new TCanvas("can_dPhiLepMET",  "can_dPhiLepMET",600,600); 
	TCanvas *can_nVertex   = new TCanvas("can_nVertex",     "can_nVertex",600,600); 
	TCanvas *can_HT        = new TCanvas("can_HT",          "can_HT",600,600); 
	TCanvas *can_nJet      = new TCanvas("can_nJet",        "can_nJet",600,600); 
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   

	//float scalefactor = p_llmass_data->Integral(20,50)/p_llmass_ZG->Integral(20,50); 
	//std::cout << "scale = " << scalefactor<<std::endl;
	//float scalefactor = p_phoEt_data->Integral(6,8)/p_phoEt_ZG->Integral(6,8);
	float scalefactor = 1; 

	can_lepPt->cd();
	gPad->SetLogy();
	p_lepPt_data->Draw();
	p_lepPt_ZG->Add(p_lepPt_rare);
	p_lepPt_ZG->SetLineColor(kRed);
	p_lepPt_ZG->Scale(scalefactor);
	p_lepPt_ZG->Draw("same");
	p_lepPt_rare->SetLineColor(kGreen);
	p_lepPt_rare->Scale(scalefactor);
	p_lepPt_rare->Draw("same");

	TCanvas *can_ratio = new TCanvas("can_ratio","",600,600);
	can_ratio->cd();
	TH1D *ratio_pt = (TH1D*)p_lepPt_data->Clone("ratio_pt");
	ratio_pt->Divide(p_lepPt_ZG);
	ratio_pt->Draw();
	for(unsigned i(1); i <= 15; i++)std::cout << i << " " << ratio_pt->GetBinContent(i) << std::endl;

	can_lepEta->cd();
	p_lepEta_data->Draw();
	p_lepEta_ZG->Add(p_lepEta_rare);
	p_lepEta_ZG->SetLineColor(kRed);
	p_lepEta_ZG->Scale(scalefactor);
	p_lepEta_ZG->Draw("same");
	p_lepEta_rare->SetLineColor(kGreen);
	p_lepEta_rare->Scale(scalefactor);
	p_lepEta_rare->Draw("same");

	can_lepPhi->cd();
	p_lepPhi_data->Draw();
	p_lepPhi_ZG->Add(p_lepPhi_rare);
	p_lepPhi_ZG->SetLineColor(kRed);
	p_lepPhi_ZG->Scale(scalefactor);
	p_lepPhi_ZG->Draw("same");
	p_lepPhi_rare->SetLineColor(kGreen);
	p_lepPhi_rare->Scale(scalefactor);
	p_lepPhi_rare->Draw("same");

	can_sigMET->cd();
	p_sigMET_data->Draw();
	p_sigMET_ZG->Add(p_sigMET_rare);
	p_sigMET_ZG->SetLineColor(kRed);
	p_sigMET_ZG->Scale(scalefactor);
	p_sigMET_ZG->Draw("same");
	p_sigMET_rare->SetLineColor(kGreen);
	p_sigMET_rare->Scale(scalefactor);
	p_sigMET_rare->Draw("same");

	can_sigMT->cd();
	p_sigMT_data->Draw();
	p_sigMT_ZG->Add(p_sigMT_rare);
	p_sigMT_ZG->SetLineColor(kRed);
	p_sigMT_ZG->Scale(scalefactor);
	p_sigMT_ZG->Draw("same");
	p_sigMT_rare->SetLineColor(kGreen);
	p_sigMT_rare->Scale(scalefactor);
	p_sigMT_rare->Draw("same");

	can_sigMETPhi->cd();
	p_sigMETPhi_data->Draw();
	p_sigMETPhi_ZG->Add(p_sigMETPhi_rare);
	p_sigMETPhi_ZG->SetLineColor(kRed);
	p_sigMETPhi_ZG->Scale(scalefactor);
	p_sigMETPhi_ZG->Draw("same");
	p_sigMETPhi_rare->SetLineColor(kGreen);
	p_sigMETPhi_rare->Scale(scalefactor);
	p_sigMETPhi_rare->Draw("same");

	can_dPhiLepMET->cd();
	p_dPhiLepMET_data->Draw();
	p_dPhiLepMET_ZG->Add(p_dPhiLepMET_rare);
	p_dPhiLepMET_ZG->SetLineColor(kRed);
	p_dPhiLepMET_ZG->Scale(scalefactor);
	p_dPhiLepMET_ZG->Draw("same");
	p_dPhiLepMET_rare->SetLineColor(kGreen);
	p_dPhiLepMET_rare->Scale(scalefactor);
	p_dPhiLepMET_rare->Draw("same");

	can_nVertex->cd();
	p_nVertex_data->Draw();
	p_nVertex_ZG->Add(p_nVertex_rare);
	p_nVertex_ZG->SetLineColor(kRed);
	p_nVertex_ZG->Scale(scalefactor);
	p_nVertex_ZG->Draw("same");
	p_nVertex_rare->SetLineColor(kGreen);
	p_nVertex_rare->Scale(scalefactor);
	p_nVertex_rare->Draw("same");

	can_HT->cd();
	gPad->SetLogy();
	p_HT_data->Draw();
	p_HT_ZG->Add(p_HT_rare);
	p_HT_ZG->SetLineColor(kRed);
	p_HT_ZG->Scale(scalefactor);
	p_HT_ZG->Draw("same");
	p_HT_rare->SetLineColor(kGreen);
	p_HT_rare->Scale(scalefactor);
	p_HT_rare->Draw("same");

	can_nJet->cd();
	p_nJet_data->Draw();
	p_nJet_ZG->Add(p_nJet_rare);
	p_nJet_ZG->SetLineColor(kRed);
	p_nJet_ZG->Scale(scalefactor);
	p_nJet_ZG->Draw("same");
	p_nJet_rare->SetLineColor(kGreen);
	p_nJet_rare->Scale(scalefactor);
	p_nJet_rare->Draw("same");

	can_llmass->cd();
	p_llmass_data->Draw();
	p_llmass_ZG->Add(p_llmass_rare);
	p_llmass_ZG->SetLineColor(kRed);
	p_llmass_ZG->Scale(scalefactor);
	p_llmass_ZG->Draw("same");
	p_llmass_rare->SetLineColor(kGreen);
	p_llmass_rare->Scale(scalefactor);
	p_llmass_rare->Draw("same");

}


