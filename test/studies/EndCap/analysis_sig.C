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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

void analysis_sig(){

	std::ifstream configfile("BkgPredConfig.txt");
	int ichannel(1);
	int anatype(0);
	int lowMt(0);
	int highMt(-1);
	int lowMET(0);
	int highMET(-1);
	int lowPt(25);
	int highPt(-1);
	int lepIso(4);
	std::string conftype;
	double confvalue;
	if(configfile.is_open()){
  	for(int i(0); i<8; i++){ 
			configfile >> conftype >> confvalue; 
			if(conftype.find("ichannel")!=std::string::npos)ichannel = confvalue;
			if(conftype.find("anatype")!=std::string::npos)anatype = confvalue;
			if(conftype.find("lowMt")!=std::string::npos)lowMt = confvalue;
			if(conftype.find("highMt")!=std::string::npos)highMt = confvalue;
			if(conftype.find("lowMET")!=std::string::npos)lowMET = confvalue;
			if(conftype.find("highMET")!=std::string::npos)highMET = confvalue;
			if(conftype.find("lowPt")!=std::string::npos)lowPt = confvalue;
			if(conftype.find("highPt")!=std::string::npos)highPt = confvalue;
			if(conftype.find("lepIso")!=std::string::npos)lepIso = confvalue;
	  }
	}
	configfile.close();

	std::ifstream binfile("binConfig.txt");
	float METbin1(200), METbin2(300);
	if(binfile.is_open()){
		for(int i(0); i<2; i++){
			binfile >> conftype >> confvalue;
			if(conftype.find("METbin1")!=std::string::npos)METbin1= confvalue;
			if(conftype.find("METbin2")!=std::string::npos)METbin2= confvalue;
		}
	}

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;
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
	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	if(channelType==1)sigtree->Add("/uscms_data/d3/mengleis/Sep1/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==2)sigtree->Add("/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_FullEcal_HT.root");

	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float sigMETPhi(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float threeMass(0);
	float nJet(0);	
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("HTClean",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	sigtree->SetBranchAddress("threeMass", &threeMass);
	sigtree->SetBranchAddress("nJet",      &nJet);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) < 1.56 || fabs(lepEta) > 2.5)continue;
		if(fabs(phoEta) > 2.4)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		p_nJet->Fill(nJet);
		
	}        

	std::ostringstream outputname;

	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_signal";
	else if(channelType==2)outputname << "mg_signal";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";

  TCanvas *can=new TCanvas("can","",1200,800);
	can->Divide(2,3);
	can->cd(1);
	gPad->SetLogy();
	p_PhoEt->Draw();
	can->cd(2);
	gPad->SetLogy();
	p_LepPt->Draw();
	can->cd(3);
	gPad->SetLogy();
	p_MET->Draw();
	can->cd(4);
	gPad->SetLogy();
	p_Mt->Draw();
	can->cd(5);
	gPad->SetLogy();
	p_HT->Draw();


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
	outputfile->Write();
	outputfile->Close();
}


