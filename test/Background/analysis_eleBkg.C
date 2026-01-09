#include<string>
#include "../../include/analysis_commoncode.h"

#define NTOY 1000
bool useGaussFit=false;

void analysis_eleBkg(){

	SetRunConfig();
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;

  	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	float correction =1.0;
        if (RunYear == 2016 and preVFP == 1 and ichannel == 1) correction = 0.948005;
        else if (RunYear == 2016 and preVFP == 1 and ichannel == 2) correction = 0.900015;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 1) correction = 0.775942;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 2) correction = 0.730433;
        else if (RunYear == 2017 and ichannel == 1) correction = 1.39352;
        else if (RunYear == 2017 and ichannel == 2) correction = 1.03558;
        else if (RunYear == 2018 and ichannel == 1) correction = 1.35192;
        else if (RunYear == 2018 and ichannel == 2) correction = 1.20607;

  	/**********************************/
	/*	double normfactor = par[0]; 	*/  
  	/*	double slope = par[1];		*/
  	/*	double constant = par[2];	*/
  	/*	double index = par[3];		*/
  	/*	double coeff = par[4]; 		*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];	*/
  	/**********************************/
	
	std::ifstream elefake_file;
	elefake_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPtVtx-EB.txt",RunYear,whichVFP.c_str()));
	// fake rate as input
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
	TF3 h_nominal_fakerate("h_nominal_fakerate", fakerate_func,10,1000,0,100,0,1.5,7);
	h_nominal_fakerate.SetParameters(scalefactor*correction, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile;
	elefake_toyfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB.txt",RunYear,whichVFP.c_str()));
	// toy fake rate as input
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >>  vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor*correction, ptslope, ptconstant, ptindex, 1, vtxconst, vtxslope); 
	  }
	}
	elefake_toyfile.close();

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

	TH1D *p_LepPt_TT = new TH1D("p_LepPt_TT","p_LepPt",nBkgPtBins,bkgPtBins);
        TH1D *p_nJet_TT = new TH1D("p_nJet_TT","p_nJet",10,0,10);
        TH1D *p_nBJet_TT = new TH1D("p_nBJet_TT","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET_TT = new TH1D("p_dPhiEleMET_TT","dPhiEleMET",32,0,3.2);

	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];

	TH1D *toy_PhoEt_TT[NTOY];
	TH1D *toy_HT_TT[NTOY];
	TH1D *toy_MET_TT[NTOY];
	TH1D *toy_Mt_TT[NTOY];
	TH1D *toy_dPhiEleMET_TT[NTOY];
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
		histname << "toy_PhoEt_TT_ " << ih;
		toy_PhoEt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
		histname.str("");
		histname << "toy_MET_TT_" << ih;
		toy_MET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
		histname.str("");
		histname << "toy_HT_TT_" << ih;
		toy_HT_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
		histname.str("");
		histname << "toy_Mt_TT_" << ih;
		toy_Mt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_TT_" << ih;
		toy_dPhiEleMET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}
	//************ Proxy Tree **********************//
	// background estimated from data, with proxyTree
	TChain *proxytree = new TChain("proxyTree");
                if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
                if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
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
	float nJetFloat(0);
	int nJetInt(0);

	int   nBJet(0);	


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
	proxytree->SetBranchAddress("HT",        &HT);
	proxytree->SetBranchAddress("nBJet",     &nBJet);
	if (channelType == 1) proxytree->SetBranchAddress("nJet", &nJetFloat);
	else proxytree->SetBranchAddress("nJet", &nJetInt);

	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		// different MET, MT and lep pT cuts
		double w_ele = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
		// fake rate weight
		p_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_MET->Fill(sigMET, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_ele);
		if (channelType == 1) p_nJet->Fill(nJetFloat, w_ele);
                if (channelType == 2) p_nJet->Fill(nJetInt, w_ele);

		//p_nJet->Fill(nJet, w_ele);
		p_nBJet->Fill(nBJet, w_ele);
		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_MET[it]->Fill(sigMET, toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET), toy_ele);
		}

		// ttbar events when nBJets >= 1
		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  w_ele);
			p_MET_TT->Fill(sigMET,  w_ele);
			p_Mt_TT->Fill(sigMT,  w_ele);
			p_HT_TT->Fill(HT,  w_ele);
			p_LepPt_TT->Fill(lepPt,  w_ele);
                        if (channelType == 1) p_nJet_TT->Fill(nJetFloat,  w_ele);
                        if (channelType == 2) p_nJet_TT->Fill(nJetInt,  w_ele);
                        p_nBJet_TT->Fill(nBJet, w_ele);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), w_ele);
			for(unsigned it(0); it < NTOY; it++){
				double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
				toy_PhoEt_TT[it]->Fill(phoEt,toy_ele);
				toy_MET_TT[it]->Fill(sigMET, toy_ele);
				toy_HT_TT[it]->Fill(HT, toy_ele);
				toy_Mt_TT[it]->Fill(sigMT, toy_ele);
				toy_dPhiEleMET_TT[it]->Fill(fabs(dPhiLepMET), toy_ele);
			}
		}

	}


	// Error is assigned with Gauss fit
	std::vector<double> toyvec; 
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); // systematic error from Toy study
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin)); // stat + syst error
		p_PhoEt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_LepPt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_LepPt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_LepPt->GetBinError(ibin)*p_LepPt->GetBinError(ibin));
		p_LepPt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET->GetBinError(ibin)*p_MET->GetBinError(ibin));
		p_MET->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_Mt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt->GetBinError(ibin)*p_Mt->GetBinError(ibin));
		p_Mt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT->GetBinError(ibin)*p_HT->GetBinError(ibin));
		p_HT->SetBinError(ibin, totalerror);
	}
		
	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt_TT->GetBinError(ibin)*p_PhoEt_TT->GetBinError(ibin));
		p_PhoEt_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET_TT->GetBinError(ibin)*p_MET_TT->GetBinError(ibin));
		p_MET_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT_TT->GetBinError(ibin)*p_HT_TT->GetBinError(ibin));
		p_HT_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_Mt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt_TT->GetBinError(ibin)*p_Mt_TT->GetBinError(ibin));
		p_Mt_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_dPhiEleMET_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_dPhiEleMET_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_dPhiEleMET_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_dPhiEleMET_TT->GetBinError(ibin)*p_dPhiEleMET_TT->GetBinError(ibin));
		p_dPhiEleMET_TT->SetBinError(ibin, totalerror);
	}
	std::ostringstream outputname;
	outputname << "/eos/uscms/store/user/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_eleBkg";
	else if(channelType==2)outputname << "mg_eleBkg";
	if(anatype ==0 or anatype ==1)	outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname <<"_" << RunYear<<whichVFP <<".root";

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
	p_nBJet->Write();
	p_PhoEt_TT->Write();
	p_dPhiEleMET_TT->Write();
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	for(unsigned it(0); it < NTOY; it++){
		toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();
}
