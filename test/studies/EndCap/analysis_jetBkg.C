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

void analysis_jetBkg(){

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

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double randomweight_jet[1000];
	randomweight_jet[0] = 0;
	for(unsigned ir(1); ir<1000; ir++)	
		randomweight_jet[ir] = -1+ gRandom->Rndm()*2.0;

	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,1000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,1000,4);
	double jetfake_numerror[264];
	double jetfake_denerror[264];
	
	std::stringstream JetFakeRateFile;
  JetFakeRateFile.str();
	//if(channelType==1)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-DoubleEG-EE.txt";
	if(channelType==1)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-MuonEG-EE.txt";
	if(channelType==2)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-MuonEG-EE.txt";
	std::ifstream jetfakefile(JetFakeRateFile.str().c_str());
	std::string paratype;
	float paravalue;	
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_den->SetParameter(i, paravalue);
	}
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_num->SetParameter(i, paravalue);
	}
	int binnumber;
	for(int i(0); i < 264; i++){
		jetfakefile >> paratype >> binnumber >> paravalue;
		jetfake_numerror[i] = paravalue;
	//	std::cout << "error " << jetfake_numerror[i] << std::endl;
	}
	for(int i(0); i < 264; i++){
		jetfakefile >> paratype >> binnumber >>  paravalue;
		jetfake_denerror[i] = paravalue;
	//	std::cout << "error " << jetfake_denerror[i] << std::endl;
	}

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

	TH1D *h_jetfakepho_norm            = new TH1D("h_jetfakepho_norm","eventcount",9,0,9);
	TH1D *h_jetfakepho_syserr_jes      = new TH1D("h_jetfakepho_syserr_jes","",9,0,9);	
	TH1D *h_jetfakepho_syserr_jer      = new TH1D("h_jetfakepho_syserr_jer","",9,0,9);	
	TH1D *h_jetfakepho_syserr_esf      = new TH1D("h_jetfakepho_syserr_esf","",9,0,9);	
	TH1D *h_jetfakepho_syserr_scale    = new TH1D("h_jetfakepho_syserr_scale","",9,0,9);	
	TH1D *h_jetfakepho_syserr_eleshape = new TH1D("h_jetfakepho_syserr_eleshape","",9,0,9);	
	TH1D *h_jetfakepho_syserr_jetshape = new TH1D("h_jetfakepho_syserr_jetshape","",9,0,9);	
	TH1D *h_jetfakepho_syserr_xs       = new TH1D("h_jetfakepho_syserr_xs","",9,0,9);	
	TH1D *h_jetfakepho_syserr_lumi     = new TH1D("h_jetfakepho_syserr_lumi","",9,0,9);
	
	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_HT[NTOY];
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
	TH1D *p_transfer = new TH1D("p_transfer","",965,35,1000);
	for(int ibin(1); ibin <= 965; ibin++){
		p_transfer->SetBinContent(ibin, fitfunc_num->Eval(35+ibin)/fitfunc_den->Eval(35+ibin));
	}

	/************ jet tree **************************/ 
		TChain *jettree = new TChain("jetTree");
		if(channelType==1)jettree->Add("/uscms_data/d3/mengleis/Sep1/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
		if(channelType==2)jettree->Add("/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_FullEcal_HT.root");

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
		float threeMass(0);
		float HT(0);
		float nJet(0);
		
		jettree->SetBranchAddress("phoEt",     &phoEt);
		jettree->SetBranchAddress("phoEta",    &phoEta);
		jettree->SetBranchAddress("phoPhi",    &phoPhi);
		jettree->SetBranchAddress("lepPt",     &lepPt);
		jettree->SetBranchAddress("lepEta",    &lepEta);
		jettree->SetBranchAddress("lepPhi",    &lepPhi);
		jettree->SetBranchAddress("sigMT",     &sigMT);
		jettree->SetBranchAddress("sigMET",    &sigMET);
		jettree->SetBranchAddress("sigMETPhi", &sigMETPhi);
		jettree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
		jettree->SetBranchAddress("nVertex",   &nVertex);
		jettree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
		jettree->SetBranchAddress("threeMass", &threeMass);
		jettree->SetBranchAddress("HTClean",        &HT);
		jettree->SetBranchAddress("nJet",      &nJet);
	 
		for (unsigned ievt(0); ievt<jettree->GetEntries(); ++ievt){//loop on entries
			jettree->GetEntry(ievt);
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

			double w_jet(0);
			w_jet = fitfunc_num->Eval(phoEt)/fitfunc_den->Eval(phoEt);

			double jetfakeerror(0);
			for(int ipt(0); ipt < 264; ipt++){
				if(phoEt >= ipt+35 && phoEt < ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt]/fitfunc_den->Eval(phoEt)/fitfunc_den->Eval(phoEt) + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(phoEt);
			}
			if(phoEt >= 264)jetfakeerror = sqrt(jetfake_numerror[263]*jetfake_numerror[263]/fitfunc_den->Eval(300)/fitfunc_den->Eval(300) + jetfake_denerror[263]*jetfake_denerror[263]*w_jet*w_jet)/fitfunc_den->Eval(300); 
			double sysJetFakePho = jetfakeerror/w_jet;
		
			p_PhoEt->Fill(phoEt,w_jet);
			p_PhoEta->Fill(phoEta, w_jet);
			p_MET->Fill(sigMET,w_jet);
			p_Mt->Fill(sigMT, w_jet);
			p_HT->Fill(HT, w_jet);
			p_LepPt->Fill(lepPt, w_jet);
			p_LepEta->Fill(lepEta, w_jet);
			p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_jet);
			p_nJet->Fill(nJet, w_jet);

			for(unsigned ih(0); ih<NTOY; ih++){
				
				toy_PhoEt[ih]->Fill(phoEt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_MET[ih]->Fill(sigMET,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_Mt[ih]->Fill(sigMT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_HT[ih]->Fill(HT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_LepPt[ih]->Fill(lepPt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_dPhiEleMET[ih]->Fill(fabs(dPhiLepMET), w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
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
	if(channelType==1)outputname << "egamma_jetbkg";
	else if(channelType==2)outputname << "mg_jetbkg";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
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
	h_jetfakepho_syserr_jes->Write();       
	h_jetfakepho_syserr_jer->Write();       
	h_jetfakepho_syserr_esf->Write();       
	h_jetfakepho_syserr_scale->Write();     
	h_jetfakepho_syserr_eleshape->Write();  
	h_jetfakepho_syserr_jetshape->Write();  
	h_jetfakepho_syserr_xs->Write();        
	h_jetfakepho_syserr_lumi->Write();      
	for(unsigned it(0); it < NTOY; it++){
		toy_PhoEt[it]->Write();
		toy_MET[it]->Write();
		toy_Mt[it]->Write();
		toy_HT[it]->Write();
		toy_LepPt[it]->Write();
		toy_dPhiEleMET[it]->Write();
	}
	p_transfer->Write();
	outputfile->Write();
	outputfile->Close();

}


