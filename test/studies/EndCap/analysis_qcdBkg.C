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

void analysis_qcdBkg(){

	bool toDeriveScale(false);

	double factor_egQCD(0.43);
	double factorerror_egQCD(sqrt(0.011*0.011 + 0.017*0.017));
	//double factor_mgQCD(0.416174);
	double factor_mgQCD(0.4);
	double factorerror_mgQCD(sqrt(0.00044*0.00044 + 0.0229*0.0229));

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
  	for(int i(0); i<9; i++){ 
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
	double factorQCD(1);
	double factorQCDUP = 1; 
	double factorQCDDO = 1; 

	if(channelType == 1){
		factorQCD = factor_egQCD;
		factorQCDUP = factor_egQCD + factorerror_egQCD;
		factorQCDDO = factor_egQCD - factorerror_egQCD;
	}
	else if(channelType == 2){
		factorQCD = factor_mgQCD;
		factorQCDUP = factor_mgQCD + factorerror_mgQCD;
		factorQCDDO = factor_mgQCD - factorerror_mgQCD;
	}

	if(anatype == 0){
		toDeriveScale = true;
		factorQCD = 1;
		factorQCDUP = 1;
		factorQCDDO = 1;
	}

	TFile *scaleFile;
	if(channelType == 1)scaleFile = TFile::Open("qcd_eg_scale.root");
//	else if(channelType == 2)scaleFile = TFile::Open("qcd_mg_scale.root");
//	TH1D *p_scale;
//	if(channelType == 1)p_scale = (TH1D*)scaleFile->Get("p_scale");


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

	TH1D *h_qcdfakelep_norm            = new TH1D("h_qcdfakelep_norm","eventcount",9,0,9);
	TH1D *h_qcdfakelep_normup          = new TH1D("h_qcdfakelep_normup","eventcount; eventcount (GeV);",9,0,9);
	TH1D *h_qcdfakelep_normdo          = new TH1D("h_qcdfakelep_normdo","eventcount; eventcount (GeV);",9,0,9);
	TH1D *h_qcdfakelep_syserr_jes      = new TH1D("h_qcdfakelep_syserr_jes","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_jer      = new TH1D("h_qcdfakelep_syserr_jer","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_esf      = new TH1D("h_qcdfakelep_syserr_esf","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_scale    = new TH1D("h_qcdfakelep_syserr_scale","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_eleshape = new TH1D("h_qcdfakelep_syserr_eleshape","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_jetshape = new TH1D("h_qcdfakelep_syserr_jetshape","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_xs       = new TH1D("h_qcdfakelep_syserr_xs","",9,0,9);	
	TH1D *h_qcdfakelep_syserr_lumi     = new TH1D("h_qcdfakelep_syserr_lumi","",9,0,9);
	
	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *normdo_PhoEt = new TH1D("normdo_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normdo_PhoEta = new TH1D("normdo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normdo_LepPt = new TH1D("normdo_LepPt","normdo_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normdo_LepEta = new TH1D("normdo_LepEta","normdo_LepEta",60,-3,3);
	TH1D *normdo_MET = new TH1D("normdo_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normdo_Mt = new TH1D("normdo_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normdo_HT = new TH1D("normdo_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normdo_dPhiEleMET = new TH1D("normdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/Sep1/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_FullEcal_HT.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsMedium(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
	float threeMass(0);
  float HT(0);
  float nJet(0);
  
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsMedium",&fakeLepIsMedium);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	fakeEtree->SetBranchAddress("threeMass", &threeMass);
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nJet",      &nJet);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		if(ievt%1000 ==0)std::cout <<"event " << ievt << std::endl;
		//if(ichannel == 1){
		//	factorQCD = p_scale->GetBinContent(p_scale->FindBin(lepPt));
		//	factorQCDUP = p_scale->GetBinContent(p_scale->FindBin(lepPt)) + p_scale->GetBinError(p_scale->FindBin(lepPt));
		//	factorQCDDO = p_scale->GetBinContent(p_scale->FindBin(lepPt)) - p_scale->GetBinError(p_scale->FindBin(lepPt));
		//}
	  if(toDeriveScale){
		 factorQCD = 1;
		 factorQCDUP = 1;
		 factorQCDDO = 1;
	  }
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

		bool isProxy(false);
		if(channelType==1){if(fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1 && fakeLepIsMedium>0))isProxy=true;}
		if(!isProxy)continue;

		p_PhoEt->Fill(phoEt, factorQCD);
		p_PhoEta->Fill(phoEta, factorQCD);
		p_LepPt->Fill(lepPt, factorQCD);
		p_LepEta->Fill(lepEta, factorQCD);
		p_MET->Fill(sigMET, factorQCD);
		p_Mt->Fill(sigMT, factorQCD);
		p_HT->Fill(HT, factorQCD);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), factorQCD);
		p_nJet->Fill(nJet, factorQCD);	

		normup_PhoEt->Fill(phoEt, factorQCDUP);
		normup_PhoEta->Fill(phoEta,factorQCDUP);
		normup_LepPt->Fill(lepPt, factorQCDUP);
		normup_LepEta->Fill(lepEta,factorQCDUP);
		normup_MET->Fill(sigMET, factorQCDUP);
		normup_Mt->Fill(sigMT, factorQCDUP);
		normup_HT->Fill(HT, factorQCDUP);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), factorQCDUP);

		normdo_PhoEt->Fill(phoEt, factorQCDDO);
		normdo_PhoEta->Fill(phoEta,factorQCDDO);
		normdo_LepPt->Fill(lepPt, factorQCDDO);
		normdo_LepEta->Fill(lepEta,factorQCDDO);
		normdo_MET->Fill(sigMET, factorQCDDO);
		normdo_Mt->Fill(sigMT, factorQCDDO);
		normdo_HT->Fill(HT, factorQCDDO);
		normdo_dPhiEleMET->Fill(fabs(dPhiLepMET), factorQCDDO);
	}

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	

	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_qcd";
	else if(channelType==2)outputname << "mg_qcd";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << "_iso" << lepIso;
	outputname << ".root";
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
	p_nJet->Write();
	h_qcdfakelep_norm->Write();             
	h_qcdfakelep_syserr_jes->Write();       
	h_qcdfakelep_syserr_jer->Write();       
	h_qcdfakelep_syserr_esf->Write();       
	h_qcdfakelep_syserr_scale->Write();     
	h_qcdfakelep_syserr_eleshape->Write();  
	h_qcdfakelep_syserr_jetshape->Write();  
	h_qcdfakelep_syserr_xs->Write();        
	h_qcdfakelep_syserr_lumi->Write();      
	outputfile->Write();
	outputfile->Close();
}


