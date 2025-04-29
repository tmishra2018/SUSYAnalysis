// g++ `root-config --cflags` ../../lib/libAnaClasses.so plot_ISRweight.C -o plot_ISRweight.exe `root-config --libs`
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TROOT.h"
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
#include "TGraphErrors.h"


#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"
#include "../../include/tdrstyle.C"

// ../../include/analysis_scalefactor.h
// change the file /uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/Background/BkgPredConfig.txt
//bool doEB=false;
bool ApplyISR = false;

void plot_ISRweight(int RunYear, bool preVFP){


	setTDRStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetErrorX(0.5);

  	gSystem->Load("../../lib/libAnaClasses.so");
	
	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	gROOT->SetBatch(kTRUE);
	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={0,50,100,150,200,250,300,800};
	TH1F *p_phoEt_data     = new TH1F("p_phoEt_data","",7,plotEtBins);
	TH1F *p_JetPt_data     = new TH1F("p_JetPt_data","",7,plotEtBins);
	TH1F *p_JetPt_total    = new TH1F("p_JetPt_total","",7,plotEtBins);
	TH1F *p_llmass_data    = new TH1F("p_llmass_data","",100,30,130);
	TH1F *p_phoEt_data_highEt   = new TH1F("p_phoEt_data_highEt","",7,plotEtBins);
	TH1F *p_JetPt_data_highEt   = new TH1F("p_JetPt_data_highEt","",7,plotEtBins);
	TH1F *p_llmass_data_highEt  = new TH1F("p_llmass_data_highEt","",100,30,130);

	TH1F *p_phoEt_ZG     = new TH1F("p_phoEt_ZG","",7,plotEtBins);
	TH1F *p_JetPt_ZG     = new TH1F("p_JetPt_ZG","",7,plotEtBins);
	TH1F *p_llmass_ZG    = new TH1F("p_llmass_ZG","",100,30,130);
	TH1F *p_phoEt_ZG_highEt   = new TH1F("p_phoEt_ZG_highEt","",7,plotEtBins);
	TH1F *p_JetPt_ZG_highEt   = new TH1F("p_JetPt_ZG_highEt","",7,plotEtBins);
	TH1F *p_llmass_ZG_highEt  = new TH1F("p_llmass_ZG_highEt","",100,30,130);

	TH1F *p_phoEt_rare     = new TH1F("p_phoEt_rare","",7,plotEtBins);
	TH1F *p_JetPt_rare     = new TH1F("p_JetPt_rare","",7,plotEtBins);
	TH1F *p_llmass_rare    = new TH1F("p_llmass_rare","",100,30,130);
	TH1F *p_phoEt_rare_highEt   = new TH1F("p_phoEt_rare_highEt","",7,plotEtBins);
	TH1F *p_JetPt_rare_highEt   = new TH1F("p_JetPt_rare_highEt","",7,plotEtBins);
	TH1F *p_llmass_rare_highEt  = new TH1F("p_llmass_rare_highEt","",100,30,130);

	TProfile *p_scalefactor  = new TProfile("p_scalefactor","p_scalefactor",50,-2.5,2.5);
	
	//************ Signal Tree **********************//
	TChain *tree = new TChain("ZTree");
	// using analysis_ISR.C
  	tree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_%d%s.root",RunYear,whichVFP.c_str()));
  	float phoEt(0);
  	float phoEta(0);
  	float phoPhi(0);
  	float lepPt(0);
  	float lepEta(0);
  	float lepPhi(0);
  	float sigMT(0);
  	float sigMET(0);
  	float dPhiLepMET(0);
  	float dRPhoLep(0);
  	float HT(0);
  	float nJet(0);
	float threeMass(0);
	float dilepMass(0);
	float JetPt=0;

  	tree->SetBranchAddress("phoEt",     &phoEt);
  	tree->SetBranchAddress("phoEta",    &phoEta);
  	tree->SetBranchAddress("phoPhi",    &phoPhi);
  	tree->SetBranchAddress("lepPt",     &lepPt);
  	tree->SetBranchAddress("lepEta",    &lepEta);
  	tree->SetBranchAddress("lepPhi",    &lepPhi);
  	tree->SetBranchAddress("sigMT",     &sigMT);
  	tree->SetBranchAddress("sigMET",    &sigMET);
  	tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  	tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  	tree->SetBranchAddress("HT",        &HT);
  	tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("threeMass", &threeMass);
	tree->SetBranchAddress("dilepMass", &dilepMass);
	tree->SetBranchAddress("ISRJetPt",  &JetPt);

  	for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
		// always set 0 as doEB=0, so both EB and EE considered
		//if(doEB && fabs(phoEta) > 1.4442)continue;
		
		//else if(!doEB && (fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4))continue;

		if(phoEt > 799)phoEt = 799;
		if(JetPt > 799)JetPt = 799;
		if(dRPhoLep < 0.8)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		double weight = 1;
		
		std::stringstream fakerate_filename;
        	fakerate_filename.str("");
        	fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-";

		if(fabs(phoEta) < 1.4442) fakerate_filename << "EB_";
		else if(fabs(phoEta) >= 1.4442 ) fakerate_filename << "EE_";

		if (RunYear==2016 and preVFP==true) fakerate_filename << "2016preVFP.txt";
		else if (RunYear==2016 and preVFP==false) fakerate_filename << "2016postVFP.txt";
		else if (RunYear==2017) fakerate_filename << "2017.txt";
		else if (RunYear==2018) fakerate_filename << "2018.txt";
	
		int i(0);
		float pt_lower(0), pt_upper(0);
	        float fakerate(0), error(0), systematic(0);
        	float truefake;
		
		if(dilepMass > 80 && dilepMass < 100){
				std::ifstream jetfake_file(fakerate_filename.str().c_str());
				if(jetfake_file.is_open()){
                			for(int i(0); i < 9; i++){
                        			jetfake_file >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truefake;
						if(phoEt > pt_lower && phoEt < pt_upper) weight = 1-fakerate;            
					}
                			jetfake_file.close();
        			}
		}
		
	 	p_JetPt_data->Fill(JetPt, weight);
		p_phoEt_data->Fill(phoEt, weight); 
		p_llmass_data->Fill(dilepMass, weight);

		//if(phoEt > 145){
		//	p_llmass_data_highEt->Fill(dilepMass, weight);
		//	p_JetPt_data_highEt->Fill(JetPt, weight);
		//	p_phoEt_data_highEt->Fill(phoEt, weight);
		//}
	}//loop on  events


	double reweightF(1);
        double Normalization=1;
	float ISRWeight(0);

	//************ Signal Tree **********************//
  	TChain *ZGtree = new TChain("ZTree");
	// using analysis_ISRMC.C
  	ZGtree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_ZGToLLG_%d%s.root",RunYear,whichVFP.c_str()));
	float ZG_MCweight(0);
	float ZG_PUweight(0);
  	float ZG_phoEt(0);
  	float ZG_phoEta(0);
  	float ZG_phoPhi(0);
  	float ZG_lepPt(0);
  	float ZG_lepEta(0);
  	float ZG_lepPhi(0);
	float ZG_bosonPt(0);
  	float ZG_sigMT(0);
  	float ZG_sigMET(0);
  	float ZG_dPhiLepMET(0);
  	float ZG_dRPhoLep(0);
  	float ZG_HT(0);
  	float ZG_nJet(0);
	float ZG_JetPt(0);
	float ZG_threeMass(0);
	float ZG_dilepMass(0);
  	std::vector<int>   *ZG_mcPID=0;
  	std::vector<float> *ZG_mcEta=0;
  	std::vector<float> *ZG_mcPhi=0;
  	std::vector<float> *ZG_mcPt=0;
  	std::vector<int>   *ZG_mcMomPID=0;
  	std::vector<int>   *ZG_mcGMomPID=0;

	// MC weight is Lumi weight
	ZGtree->SetBranchAddress("MCweight",  &ZG_MCweight); 
	ZGtree->SetBranchAddress("PUweight",  &ZG_PUweight);
  	ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
  	ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
  	ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
  	ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  	ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  	ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
	ZGtree->SetBranchAddress("bosonPt",   &ZG_bosonPt);
  	ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  	ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  	ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  	ZGtree->SetBranchAddress("dRPhoLep",  &ZG_dRPhoLep);
  	ZGtree->SetBranchAddress("HT",        &ZG_HT);
  	ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
	ZGtree->SetBranchAddress("threeMass", &ZG_threeMass);
	ZGtree->SetBranchAddress("dilepMass", &ZG_dilepMass);
	ZGtree->SetBranchAddress("JetPt",  &ZG_JetPt);
  	ZGtree->SetBranchAddress("mcPID",     &ZG_mcPID);
  	ZGtree->SetBranchAddress("mcEta",     &ZG_mcEta);
  	ZGtree->SetBranchAddress("mcPhi",     &ZG_mcPhi);
  	ZGtree->SetBranchAddress("mcPt",      &ZG_mcPt);
  	ZGtree->SetBranchAddress("mcMomPID",  &ZG_mcMomPID);
  	ZGtree->SetBranchAddress("mcGMomPID", &ZG_mcGMomPID);

  	for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		//if(doEB && fabs(ZG_phoEta) > 1.4442)continue;
		//if(!doEB && (fabs(ZG_phoEta)< 1.56 || fabs(ZG_phoEta) > 2.4))continue;
		double scalefactor = objectESF.getMuonESF(ZG_lepPt,ZG_lepEta)*objectESF.getPhotonESF(ZG_phoEt,ZG_phoEta)*objectESF.getMuonEGTRGESF(ZG_phoEt, ZG_lepPt);

		if(ZG_phoEt > 799)ZG_phoEt = 799;
		if(ZG_JetPt > 799)ZG_JetPt = 799;
		if(ZG_dRPhoLep < 0.8)continue;
		if(ZG_dilepMass < 80 || ZG_dilepMass >100)continue;
		if(RunYear==2016 && preVFP==1){
			if(ZG_JetPt < 50)reweightF = 1.07048;
	                 else if(ZG_JetPt >= 50 && ZG_JetPt < 100)reweightF  = 1.33745;
        	         else if(ZG_JetPt >= 100 && ZG_JetPt < 150)reweightF = 1.11197;
               	  	 else if(ZG_JetPt >= 150 && ZG_JetPt < 200)reweightF = 0.921799;
                	 else if(ZG_JetPt >= 200 && ZG_JetPt < 250)reweightF = 1.02533;
                 	 else if(ZG_JetPt >= 250 && ZG_JetPt < 300)reweightF = 0.986267;
                 	 else if(ZG_JetPt >= 300)reweightF = 0.857941;
                   	 Normalization = 0.853805;    }

		else if(RunYear==2016 && preVFP==0){
			if(ZG_JetPt < 50)reweightF = 1.0823;
                 	else if(ZG_JetPt >= 50 && ZG_JetPt < 100)reweightF  = 1.36705;
                 	else if(ZG_JetPt >= 100 && ZG_JetPt < 150)reweightF = 1.03422;
                 	else if(ZG_JetPt >= 150 && ZG_JetPt < 200)reweightF = 0.996796;
                 	else if(ZG_JetPt >= 200 && ZG_JetPt < 250)reweightF = 0.90385;
                 	else if(ZG_JetPt >= 250 && ZG_JetPt < 300)reweightF = 0.780053;
                 	else if(ZG_JetPt >= 300)reweightF = 0.705875;
                 	Normalization = 0.854942;    }

		else if(RunYear==2017){
			if(ZG_JetPt < 50)reweightF = 1.0457;
	                else if(ZG_JetPt >= 50 && ZG_JetPt < 100)reweightF  = 1.28507;
        	        else if(ZG_JetPt >= 100 && ZG_JetPt < 150)reweightF = 1.05021;
           	      	else if(ZG_JetPt >= 150 && ZG_JetPt < 200)reweightF = 0.895376;
                 	else if(ZG_JetPt >= 200 && ZG_JetPt < 250)reweightF = 0.986146;
                 	else if(ZG_JetPt >= 250 && ZG_JetPt < 300)reweightF = 0.792056;
                 	else if(ZG_JetPt >= 300)reweightF = 0.929084;
                 	Normalization = 0.889139;    }

		else if(RunYear==2018){
			if(ZG_JetPt < 50)reweightF = 1.17242;
                 	else if(ZG_JetPt >= 50 && ZG_JetPt < 100)reweightF  = 1.33113;
                 	else if(ZG_JetPt >= 100 && ZG_JetPt < 150)reweightF = 0.957966;
                 	else if(ZG_JetPt >= 150 && ZG_JetPt < 200)reweightF = 0.921558;
                 	else if(ZG_JetPt >= 200 && ZG_JetPt < 250)reweightF = 0.805571;
                 	else if(ZG_JetPt >= 250 && ZG_JetPt < 300)reweightF = 0.904798;
                 	else if(ZG_JetPt >= 300)reweightF = 0.828803;
                 	Normalization = 0.860178;    }
		ISRWeight = reweightF*Normalization;
		
		double weight;
		if (ApplyISR) weight =  ZG_MCweight*ZG_PUweight*scalefactor*ISRWeight;
		else weight = ZG_MCweight*ZG_PUweight*scalefactor;
 
		bool   isTruePho(false);
		double mindR(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC< ZG_mcPID->size(); iMC++){
			double dR1 = DeltaR((*ZG_mcEta)[iMC], (*ZG_mcPhi)[iMC], ZG_phoEta, ZG_phoPhi);
			double dE1 = fabs((*ZG_mcPt)[iMC] - ZG_phoEt)/ZG_phoEt;
			if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
		}
		if(mindR < 0.1){
			if((*ZG_mcPID)[phoIndex] == 22 && (fabs((*ZG_mcMomPID)[phoIndex]) == 23 || fabs((*ZG_mcMomPID)[phoIndex]) == 13 || fabs((*ZG_mcMomPID)[phoIndex])==999)){
				isTruePho=true;
			}
		}

		p_JetPt_ZG->Fill(ZG_JetPt, weight);	
		p_phoEt_ZG->Fill(ZG_phoEt, weight); 
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events


	//************ Signal Tree **********************//
  	TChain *raretree = new TChain("ZTree");
	// using analysis_ISRMC.C
  	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_TTGJets_%d%s.root",RunYear,whichVFP.c_str()));
  	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_TTJets_%d%s.root",RunYear,whichVFP.c_str()));
  	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_WWG_%d%s.root",RunYear,whichVFP.c_str()));
  	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_WZG_%d%s.root",RunYear,whichVFP.c_str()));

	float rare_MCweight(0);
	float rare_PUweight(0);
  	float rare_phoEt(0);
  	float rare_phoEta(0);
  	float rare_phoPhi(0);
  	float rare_lepPt(0);
  	float rare_lepEta(0);
  	float rare_lepPhi(0);
  	float rare_sigMT(0);
  	float rare_sigMET(0);
  	float rare_dPhiLepMET(0);
  	float rare_dRPhoLep(0);
  	float rare_HT(0);
  	float rare_nJet(0);
	float rare_threeMass(0);
	float rare_dilepMass(0);
	float rare_JetPt(0);
  	std::vector<int>   *rare_mcPID=0;
  	std::vector<float> *rare_mcEta=0;
  	std::vector<float> *rare_mcPhi=0;
  	std::vector<float> *rare_mcPt=0;
  	std::vector<int>   *rare_mcMomPID=0;
  	std::vector<int>   *rare_mcGMomPID=0;

	raretree->SetBranchAddress("MCweight",  &rare_MCweight); 
	raretree->SetBranchAddress("PUweight",  &rare_PUweight);
  	raretree->SetBranchAddress("phoEt",     &rare_phoEt);
  	raretree->SetBranchAddress("phoEta",    &rare_phoEta);
  	raretree->SetBranchAddress("phoPhi",    &rare_phoPhi);
  	raretree->SetBranchAddress("lepPt",     &rare_lepPt);
  	raretree->SetBranchAddress("lepEta",    &rare_lepEta);
  	raretree->SetBranchAddress("lepPhi",    &rare_lepPhi);
  	raretree->SetBranchAddress("sigMT",     &rare_sigMT);
  	raretree->SetBranchAddress("sigMET",    &rare_sigMET);
  	raretree->SetBranchAddress("dPhiLepMET",&rare_dPhiLepMET);
  	raretree->SetBranchAddress("dRPhoLep",  &rare_dRPhoLep);
  	raretree->SetBranchAddress("HT",        &rare_HT);
  	raretree->SetBranchAddress("nJet",      &rare_nJet);
	raretree->SetBranchAddress("JetPt",  &rare_JetPt);
	raretree->SetBranchAddress("threeMass", &rare_threeMass);
	raretree->SetBranchAddress("dilepMass", &rare_dilepMass);
  	raretree->SetBranchAddress("mcPID",    &rare_mcPID);
  	raretree->SetBranchAddress("mcEta",    &rare_mcEta);
  	raretree->SetBranchAddress("mcPhi",    &rare_mcPhi);
  	raretree->SetBranchAddress("mcPt",     &rare_mcPt);
  	raretree->SetBranchAddress("mcMomPID", &rare_mcMomPID);
  	raretree->SetBranchAddress("mcGMomPID",&rare_mcGMomPID);

  	for(unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		//if(doEB && fabs(rare_phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(rare_phoEta) < 1.56|| fabs(rare_phoEta) > 2.4))continue;
		double scalefactor = objectESF.getMuonESF(rare_lepPt,rare_lepEta)*objectESF.getPhotonESF(rare_phoEt,rare_phoEta)*objectESF.getMuonEGTRGESF(rare_phoEt, rare_lepPt);

		if(rare_phoEt > 799)rare_phoEt = 799;
		if(rare_JetPt > 799)rare_JetPt = 799;
		if(rare_dRPhoLep < 0.8)continue;
		if(rare_dilepMass < 80 || rare_dilepMass > 100)continue;
                if(RunYear==2016 && preVFP==1){
                        if(rare_JetPt < 50)reweightF = 1.07048;
                         else if(rare_JetPt >= 50 && rare_JetPt < 100)reweightF  = 1.33745;
                         else if(rare_JetPt >= 100 && rare_JetPt < 150)reweightF = 1.11197;
                         else if(rare_JetPt >= 150 && rare_JetPt < 200)reweightF = 0.921799;
                         else if(rare_JetPt >= 200 && rare_JetPt < 250)reweightF = 1.02533;
                         else if(rare_JetPt >= 250 && rare_JetPt < 300)reweightF = 0.986267;
                         else if(rare_JetPt >= 300)reweightF = 0.857941;
                         Normalization = 0.853805;    }

                else if(RunYear==2016 && preVFP==0){
                        if(rare_JetPt < 50)reweightF = 1.0823;
                        else if(rare_JetPt >= 50 && rare_JetPt < 100)reweightF  = 1.36705;
                        else if(rare_JetPt >= 100 && rare_JetPt < 150)reweightF = 1.03422;
                        else if(rare_JetPt >= 150 && rare_JetPt < 200)reweightF = 0.996796;
                        else if(rare_JetPt >= 200 && rare_JetPt < 250)reweightF = 0.90385;
                        else if(rare_JetPt >= 250 && rare_JetPt < 300)reweightF = 0.780053;
                        else if(rare_JetPt >= 300)reweightF = 0.705875;
                        Normalization = 0.854942;    }

		else if(RunYear==2017){
			if(rare_JetPt < 50)reweightF = 1.0457;
                 	else if(rare_JetPt >= 50 && rare_JetPt < 100)reweightF  = 1.28507;
                 	else if(rare_JetPt >= 100 && rare_JetPt < 150)reweightF = 1.05021;
                 	else if(rare_JetPt >= 150 && rare_JetPt < 200)reweightF = 0.895376;
                 	else if(rare_JetPt >= 200 && rare_JetPt < 250)reweightF = 0.986146;
                 	else if(rare_JetPt >= 250 && rare_JetPt < 300)reweightF = 0.792056;
                 	else if(rare_JetPt >= 300)reweightF = 0.929084;
                 	Normalization = 0.889139;    }

		else if(RunYear==2018){
			  if(rare_JetPt < 50)reweightF = 1.17242;
                        else if(rare_JetPt >= 50 && rare_JetPt < 100)reweightF  = 1.33113;
                        else if(rare_JetPt >= 100 && rare_JetPt < 150)reweightF = 0.957966;
                        else if(rare_JetPt >= 150 && rare_JetPt < 200)reweightF = 0.921558;
                        else if(rare_JetPt >= 200 && rare_JetPt < 250)reweightF = 0.805571;
                        else if(rare_JetPt >= 250 && rare_JetPt < 300)reweightF = 0.904798;
                        else if(rare_JetPt >= 300)reweightF = 0.828803;
                        Normalization = 0.860178;    }

		ISRWeight = reweightF*Normalization;
		
		double weight;
		if (ApplyISR) weight =  rare_MCweight*rare_PUweight*scalefactor*ISRWeight;
		else weight = rare_MCweight*rare_PUweight*scalefactor;

		p_phoEt_rare->Fill(rare_phoEt, weight); 
		p_JetPt_rare->Fill(rare_JetPt, weight);
		p_llmass_rare->Fill(rare_dilepMass, weight);

	}//loop on  events
	

	
	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_JetPt        = new TCanvas("can_JetPt",          "can_JetPt",600,600); 
	TCanvas *can_JetPt_alter  = new TCanvas("can_JetPt_alter",    "can_JetPt_alter",600,600); 


	p_llmass_ZG->Add(p_llmass_rare);
        float scalefactor = 1.0;
        float scalefactorhighEt = 1.0;

	
	can_phoEt->cd();
	TPad *phoEt_pad1 = new TPad("phoEt_pad1", "phoEt_pad1", 0, 0.3, 1, 1.0);
	phoEt_pad1->SetBottomMargin(0);
        phoEt_pad1->Draw();
        phoEt_pad1->cd();
        gStyle->SetOptStat(0);
        phoEt_pad1->SetLogy();
        p_phoEt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
        p_phoEt_data->SetMinimum(0.05);

	p_phoEt_data->SetLineColor(kBlack);
	p_phoEt_data->SetMarkerStyle(20);
	p_phoEt_data->SetMarkerColor(kBlack);
	p_phoEt_data->Draw("P");
	p_phoEt_ZG->Add(p_phoEt_rare);
	p_phoEt_ZG->SetLineColor(6);
	p_phoEt_ZG->SetFillColor(6);
	p_phoEt_ZG->Scale(scalefactor);
	p_phoEt_ZG->Draw("hist same");
	p_phoEt_rare->SetLineColor(8);
	p_phoEt_rare->SetFillColor(8);
	p_phoEt_rare->Scale(scalefactor);
	p_phoEt_rare->Draw("hist same");
	p_phoEt_data->Draw("EP same");

	TLegend *leg_phoEt =  new TLegend(0.6,0.7,0.9,0.9);
	leg_phoEt->SetFillStyle(0);
	leg_phoEt->AddEntry(p_phoEt_data, "Data");
	leg_phoEt->AddEntry(p_phoEt_ZG,   "Z#gamma");
	leg_phoEt->AddEntry(p_phoEt_rare, "t#bar{t},t#bar{t}#gamma,WW#gamma,WZ#gamma");   
	leg_phoEt->Draw("same");
	
	can_phoEt->cd();
	TPad *phoEt_pad2 = new TPad("phoEt_pad2", "phoEt_pad2", 0, 0.05, 1, 0.3);
	phoEt_pad2->SetTopMargin(0);
	phoEt_pad2->SetBottomMargin(0.3);
	phoEt_pad2->Draw();
	phoEt_pad2->cd();
  	TLine *flatratio_phoEt = new TLine(0,1,800,1);
	TH1F *ratio_phoEt=(TH1F*)p_phoEt_data->Clone("transfer factor");
	ratio_phoEt->SetMarkerStyle(20);
	ratio_phoEt->SetLineColor(kBlack);
	ratio_phoEt->GetXaxis()->SetRangeUser(0,800);
	ratio_phoEt->GetYaxis()->SetRangeUser(0.95,1.05);
	ratio_phoEt->GetYaxis()->SetLabelSize(12);
	ratio_phoEt->SetMinimum(0.2);
	ratio_phoEt->SetMaximum(1.7);
	ratio_phoEt->Divide(p_phoEt_ZG);
	ratio_phoEt->SetTitle("");
	ratio_phoEt->GetYaxis()->SetTitle("Data/MC");
	ratio_phoEt->Draw();
	flatratio_phoEt->Draw("same");
//	can_phoEt->SaveAs(Form("/eos/uscms/store/user/tmishra/ISRweighting/PLOT_ISRweight_phoEt_%d%s.pdf",RunYear,whichVFP.c_str()));
	
	cout<<RunYear<<"  "<<whichVFP<<endl;
	double ISRwgt_norm[11];
	double ISRwgt_alter[11];
	for(int i(1); i<=7; i++){
		ISRwgt_norm[i-1] = p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
		std::cout << "norm ratio "<<i<<" "<< ISRwgt_norm[i-1] << " stat " << p_JetPt_data->GetBinError(i)/p_JetPt_ZG->GetBinContent(i)<<std::endl;
	}
	if(ApplyISR == false){
        	cout<<"			if(ISRJetPt < 50)reweightF = "<<ISRwgt_norm[0]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF  = "<<ISRwgt_norm[1]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = "<<ISRwgt_norm[2]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = "<<ISRwgt_norm[3]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = "<<ISRwgt_norm[4]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = "<<ISRwgt_norm[5]<<";"<<endl;
        	cout<<"                 else if(ISRJetPt >= 300)reweightF = "<<ISRwgt_norm[6]<<";"<<endl;}

	gStyle->SetOptStat(0);
	can_JetPt->cd();
	TPad *JetPt_pad1 = new TPad("JetPt_pad1", "JetPt_pad1", 0, 0.3, 1, 1.0);
	JetPt_pad1->SetBottomMargin(0);
	JetPt_pad1->Draw();  
	JetPt_pad1->cd();  
	gStyle->SetOptStat(0);
	JetPt_pad1->SetLogy();
	p_JetPt_data->GetXaxis()->SetTitleOffset(0.9);
	p_JetPt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
	p_JetPt_data->SetMinimum(0.05);
	p_JetPt_data->SetLineColor(kBlack);
	p_JetPt_data->SetMarkerStyle(20);
	p_JetPt_data->SetMarkerColor(kBlack);
	p_JetPt_data->Draw("P");
	p_JetPt_ZG->Add(p_JetPt_rare);
	p_JetPt_ZG->SetLineColor(6);
	p_JetPt_ZG->SetFillStyle(1001);
	p_JetPt_ZG->SetFillColor(6);
	p_JetPt_ZG->Scale(scalefactor);
	p_JetPt_ZG->Draw("hist same");
	p_JetPt_rare->SetLineColor(8);
	p_JetPt_rare->SetFillStyle(1001);
	p_JetPt_rare->SetFillColor(8);
	p_JetPt_rare->Scale(scalefactor);
	p_JetPt_rare->Draw("hist same");
	p_JetPt_data->Draw("EP same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( JetPt_pad1,1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( JetPt_pad1,2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( JetPt_pad1,3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( JetPt_pad1,4,1, 11 );
	
	if(ApplyISR == false)  cout<<"                 Normalization = "<<p_JetPt_ZG->Integral()/p_JetPt_data->Integral()<<";    }"<<endl;
	
	cout<<"ZG integral is "<<p_JetPt_ZG->Integral()<<endl;
	cout<<"Data integral is "<<p_JetPt_data->Integral()<<endl;
	cout<<"Ratio "<<p_JetPt_ZG->Integral()/p_JetPt_data->Integral()  <<endl;
	TLegend *leg_JetPt =  new TLegend(0.6,0.7,0.9,0.9);
	leg_JetPt->SetFillStyle(0);
	leg_JetPt->AddEntry(p_JetPt_data, "Data");
	leg_JetPt->AddEntry(p_JetPt_ZG,   "Z#gamma");
	leg_JetPt->AddEntry(p_JetPt_rare, "t#bar{t},t#bar{t}#gamma,WW#gamma,WZ#gamma");   
	leg_JetPt->Draw("same");

	can_JetPt->cd();
	TPad *JetPt_pad2 = new TPad("JetPt_pad2", "JetPt_pad2", 0, 0.05, 1, 0.3);
	JetPt_pad2->SetTopMargin(0);
	JetPt_pad2->SetBottomMargin(0.4);
	JetPt_pad2->Draw();
	JetPt_pad2->cd();
  	TLine *flatratio_JetPt = new TLine(0,1,800,1);
	TH1F *ratio_JetPt=(TH1F*)p_JetPt_data->Clone("transfer factor");
	ratio_JetPt->SetMarkerStyle(20);
	ratio_JetPt->SetLineColor(kBlack);
	ratio_JetPt->GetXaxis()->SetRangeUser(0,800);
	if (ApplyISR)  ratio_JetPt->GetYaxis()->SetRangeUser(0.8,1.2);
	else ratio_JetPt->GetYaxis()->SetRangeUser(0.5,1.35);
	ratio_JetPt->GetYaxis()->SetLabelSize(12);
	ratio_JetPt->Divide(p_JetPt_ZG);
	ratio_JetPt->SetTitle("");
	ratio_JetPt->GetYaxis()->SetTitle("Data/MC");
	ratio_JetPt->Draw();
	flatratio_JetPt->Draw("same");

	TLine *ratioValue_JetPt = new TLine(0, p_JetPt_data->Integral() / p_JetPt_ZG->Integral(),
                                    800, p_JetPt_data->Integral() / p_JetPt_ZG->Integral());
	ratioValue_JetPt->SetLineColor(kRed);
	ratioValue_JetPt->Draw("same");

	TLatex *text = new TLatex();
	text->SetTextSize(0.1);
	text->DrawLatexNDC(0.7, 0.5, Form("Mean: %.2f", p_JetPt_data->Integral() / p_JetPt_ZG->Integral()));

	if (ApplyISR)	can_JetPt->SaveAs(Form("/eos/uscms/store/user/tmishra/ISRweighting/PLOT_ISRweight_%d%s-ISR-weighted.pdf",RunYear,whichVFP.c_str()));
	else	can_JetPt->SaveAs(Form("/eos/uscms/store/user/tmishra/ISRweighting/PLOT_ISRweight_%d%s.pdf",RunYear,whichVFP.c_str()));

}
int main(int argc, char** argv)
{
    bool preVFP = (atoi(argv[2]) == 1);
    plot_ISRweight(atoi(argv[1]), preVFP);
    return 0;

}
