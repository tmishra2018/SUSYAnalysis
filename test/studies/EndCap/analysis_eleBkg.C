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

#define NTOY 1000

void analysis_eleBkg(){

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

  /**********************************/
	/*	double normfactor = par[0]; 	*/  
  /*	double slope = par[1];				*/
  /*	double constant = par[2];			*/
  /*	double index = par[3];				*/
  /*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  /**********************************/
	std::ifstream elefake_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-ByPtVtx-EE.txt");
	double scalefactor(0);
	double ptslope(0);
	double ptconstant(0);
	double ptindex(0);
	double ptcoeff(0);
	double vtxconst(0);
	double vtxslope(0);
	std::string variabletype;
	double variablevalue;
	if(elefake_file.is_open()){
  	for(int i(0); i<7; i++){ 
			elefake_file >> variabletype >> variablevalue; 
			if(variabletype.find("scalefactor")!=std::string::npos)scalefactor = variablevalue;
			else if(variabletype.find("ptslope")!=std::string::npos)ptslope = variablevalue;
			else if(variabletype.find("ptconstant")!=std::string::npos)ptconstant = variablevalue;
			else if(variabletype.find("ptindex")!=std::string::npos)ptindex = variablevalue;
			else if(variabletype.find("ptcoeff")!=std::string::npos)ptcoeff = variablevalue;
			else if(variabletype.find("vtxconst")!=std::string::npos)vtxconst = variablevalue;
			else if(variabletype.find("vtxslope")!=std::string::npos)vtxslope = variablevalue;
	  }
	}
	elefake_file.close();
	TF3 h_nominal_fakerate("h_nominal_fakerate", fakerate_func,10,1000,0,100,0,2.5,7);
	h_nominal_fakerate.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/ToyFakeRate_FullEcal.txt");
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >>  vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,2.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1, vtxconst, vtxslope); 
	  }
	}
	elefake_toyfile.close();

	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *test_PhoEt = new TH1D("test_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	TH1D *h_elefakepho_norm = new TH1D("h_elefakepho_norm","eventcount",9,0,9);
	TH1D *h_elefakepho_syserr_jes      = new TH1D("h_elefakepho_syserr_jes","",9,0,9);	
	TH1D *h_elefakepho_syserr_jer      = new TH1D("h_elefakepho_syserr_jer","",9,0,9);	
	TH1D *h_elefakepho_syserr_esf      = new TH1D("h_elefakepho_syserr_esf","",9,0,9);	
	TH1D *h_elefakepho_syserr_scale    = new TH1D("h_elefakepho_syserr_scale","",9,0,9);	
	TH1D *h_elefakepho_syserr_eleshape = new TH1D("h_elefakepho_syserr_eleshape","",9,0,9);	
	TH1D *h_elefakepho_syserr_jetshape = new TH1D("h_elefakepho_syserr_jetshape","",9,0,9);	
	TH1D *h_elefakepho_syserr_xs       = new TH1D("h_elefakepho_syserr_xs","",9,0,9);	
	TH1D *h_elefakepho_syserr_lumi     = new TH1D("h_elefakepho_syserr_lumi","",9,0,9);	
	TH1D *h_elefakepho_syserr_isr     = new TH1D("h_elefakepho_syserr_isr","",9,0,9);	

	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];
	TH1D *toy_eventcount[NTOY];
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_PhoEt_ " << ih;
		toy_PhoEt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgPtBins,bkgPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
		histname.str("");
		histname << "toy_eventcount_" << ih;
		toy_eventcount[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),9,0,9);
	}
	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
	if(channelType==1)proxytree->Add("/uscms_data/d3/mengleis/Sep1/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==2)proxytree->Add("/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_FullEcal_HT.root");

	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float sigMETPhi(0);
	float dPhiLepMET(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float HT(0);
	float nJet(0);
	
	proxytree->SetBranchAddress("phoEt",     &phoEt);
	proxytree->SetBranchAddress("phoEta",    &phoEta);
	proxytree->SetBranchAddress("phoPhi",    &phoPhi);
	proxytree->SetBranchAddress("lepPt",     &lepPt);
	proxytree->SetBranchAddress("lepEta",    &lepEta);
	proxytree->SetBranchAddress("lepPhi",    &lepPhi);
	proxytree->SetBranchAddress("sigMT",     &sigMT);
	proxytree->SetBranchAddress("sigMET",    &sigMET);
	proxytree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   &nVertex);
	proxytree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	proxytree->SetBranchAddress("HTClean",        &HT);
	proxytree->SetBranchAddress("nJet",      &nJet);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
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

		double w_ele = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
		std::cout << phoEt << " " << phoEta << " " << w_ele << std::endl;
		p_PhoEt->Fill(phoEt,w_ele);
		test_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_MET->Fill(sigMET, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_ele);
		p_nJet->Fill(nJet, w_ele);


		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_MET[it]->Fill(sigMET, toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET), toy_ele);
		}
	}

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_PhoEt->GetBinContent(ibin),1.5*p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_PhoEt[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin));
		p_PhoEt->SetBinError(ibin, totalerror);
		delete temphist;
	}
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_LepPt->GetBinContent(ibin),1.5*p_LepPt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_LepPt[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_LepPt->GetBinError(ibin)*p_LepPt->GetBinError(ibin));
		p_LepPt->SetBinError(ibin, totalerror);
		delete temphist;
	}
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_MET->GetBinContent(ibin),1.5*p_MET->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_MET[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_MET->GetBinError(ibin)*p_MET->GetBinError(ibin));
		p_MET->SetBinError(ibin, totalerror);
		delete temphist;
	}
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_Mt->GetBinContent(ibin),1.5*p_Mt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_Mt[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_Mt->GetBinError(ibin)*p_Mt->GetBinError(ibin));
		p_Mt->SetBinError(ibin, totalerror);
		delete temphist;
	}
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_dPhiEleMET->GetBinContent(ibin),1.5*p_dPhiEleMET->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_dPhiEleMET[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_dPhiEleMET->GetBinError(ibin)*p_dPhiEleMET->GetBinError(ibin));
		p_dPhiEleMET->SetBinError(ibin, totalerror);
		delete temphist;
	}
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		TH1D *temphist = new TH1D("temphist","",500,0.5*p_HT->GetBinContent(ibin),1.5*p_HT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)temphist->Fill(toy_HT[it]->GetBinContent(ibin));
		temphist->Fit("gaus");
		double syserr = temphist->GetFunction("gaus")->GetParameter(2);
		double totalerror = sqrt(syserr*syserr + p_HT->GetBinError(ibin)*p_HT->GetBinError(ibin));
		p_HT->SetBinError(ibin, totalerror);
		delete temphist;
	}
		
	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_eleBkg";
	else if(channelType==2)outputname << "mg_eleBkg";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	test_PhoEt->Sumw2();
	test_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	h_elefakepho_norm->Write();
	h_elefakepho_syserr_jes->Write();       
	h_elefakepho_syserr_jer->Write();       
	h_elefakepho_syserr_esf->Write();       
	h_elefakepho_syserr_scale->Write();     
	h_elefakepho_syserr_eleshape->Write();  
	h_elefakepho_syserr_jetshape->Write();  
	h_elefakepho_syserr_xs->Write();        
	h_elefakepho_syserr_lumi->Write();      
	h_elefakepho_syserr_isr->Write();      
	for(unsigned it(0); it < NTOY; it++){
	//	toy_PhoEt[it]->Write();
	//	toy_MET[it]->Write();
	//	toy_Mt[it]->Write();
	//	toy_LepPt[it]->Write();
	//	toy_HT[it]->Write();
		toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();
}


