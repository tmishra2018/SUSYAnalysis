#include<string>
#include "../../include/analysis_commoncode.h"
#define NTOY 1000
bool useGaussFit;

void analysis_jetBkg(){
	gROOT->SetBatch(kTRUE);
	SetRunConfig();
	setTDRStyle();
 
	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";


	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double randomweight_jet[NTOY];
	randomweight_jet[0] = 0;
	for(unsigned ir(1); ir<NTOY; ir++)	
		randomweight_jet[ir] = -1+ gRandom->Rndm()*2.0;
	// numerator and denominator function for fake rate
	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,1000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,1000,4);
	double jetfake_numerror[264];
	double jetfake_denerror[264];
	
	std::stringstream JetFakeRateFile;
  	JetFakeRateFile.str();

	if(channelType==1)JetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-DoubleEG-EB.txt";
	if(channelType==2)JetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-MuonEG-EB.txt";

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
		jetfake_numerror[i] = paravalue/2;
	}
	for(int i(0); i < 264; i++){
		jetfakefile >> paratype >> binnumber >>  paravalue;
		jetfake_denerror[i] = paravalue/2;
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
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	
	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];
	TH1D *toy_PhoEt_TT[NTOY];
	TH1D *toy_MET_TT[NTOY];
	TH1D *toy_Mt_TT[NTOY];
	TH1D *toy_HT_TT[NTOY];
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
		histname << "toy_PhoEt_TT_" << ih;
		toy_PhoEt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
		histname.str("");
		histname << "toy_MET_TT_" << ih;
		toy_MET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
		histname.str("");
		histname << "toy_Mt_TT_" << ih;
		toy_Mt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_HT_TT_" << ih;
		toy_HT_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
	}

	/************ jet tree **************************/ 
		// Background derived from data, as it is fake, jetTree
		TChain *jettree = new TChain("jetTree");

                if(channelType==1)jettree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
                if(channelType==2)jettree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	
		float phoEt(0);
		float phoEta(0);
		float phoPhi(0);
    		float phoChIso(0);
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
		int   nBJet(0);	
	
		jettree->SetBranchAddress("phoEt",     &phoEt);
		jettree->SetBranchAddress("phoEta",    &phoEta);
		jettree->SetBranchAddress("phoPhi",    &phoPhi);
    		jettree->SetBranchAddress("phoChIso",  &phoChIso);
		jettree->SetBranchAddress("lepPt",     &lepPt);
		jettree->SetBranchAddress("lepEta",    &lepEta);
		jettree->SetBranchAddress("lepPhi",    &lepPhi);
		jettree->SetBranchAddress("sigMT",     &sigMT);
		jettree->SetBranchAddress("sigMET",    &sigMET);
		jettree->SetBranchAddress("sigMETPhi", &sigMETPhi);
		jettree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
		jettree->SetBranchAddress("nVertex",   &nVertex);
		jettree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
		jettree->SetBranchAddress("HT",        &HT);
		jettree->SetBranchAddress("nJet",      &nJet);
		jettree->SetBranchAddress("nBJet",     &nBJet);
	 
		for (unsigned ievt(0); ievt<jettree->GetEntries(); ++ievt){//loop on entries
			jettree->GetEntry(ievt);
			p_PU->Fill(nVertex);
			/** cut flow *****/
			if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
			if(sigMET < lowMET)continue;
			if(highMET > 0 && sigMET > highMET)continue;
			if(sigMT < lowMt)continue;
			if(highMt > 0 && sigMT > highMt)continue;
			if(lepPt < lowPt)continue;
			if(highPt > 0 && lepPt > highPt)continue;
			// different MET, MT and lepton pT ranges
			if(phoChIso > 5)continue;

			double w_jet(0);
			// weight for fake background, central value from function
			w_jet = fitfunc_num->Eval(phoEt)/fitfunc_den->Eval(phoEt);

			double jetfakeerror(0);
			for(int ipt(0); ipt < 264; ipt++){
				// weight for fake background, errors are stored for pt value
				if(phoEt >= ipt+35 && phoEt < ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt] + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(phoEt);
			}
			if(phoEt >= 264)jetfakeerror = sqrt(jetfake_numerror[263]*jetfake_numerror[263] + jetfake_denerror[263]*jetfake_denerror[263]*w_jet*w_jet)/fitfunc_den->Eval(300); 
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
			p_nBJet->Fill(nBJet, w_jet);

			if(nBJet >= 1){
				p_PhoEt_TT->Fill(phoEt,  w_jet);
				p_MET_TT->Fill(sigMET,  w_jet);
				p_Mt_TT->Fill(sigMT,  w_jet);
				p_HT_TT->Fill(HT,  w_jet);
			}

			for(unsigned ih(0); ih<NTOY; ih++){
				toy_PhoEt[ih]->Fill(phoEt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_MET[ih]->Fill(sigMET,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_Mt[ih]->Fill(sigMT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_HT[ih]->Fill(HT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_LepPt[ih]->Fill(lepPt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_dPhiEleMET[ih]->Fill(fabs(dPhiLepMET), w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				if(nBJet >= 1){
					toy_PhoEt_TT[ih]->Fill(phoEt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
					toy_MET_TT[ih]->Fill(sigMET,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
					toy_Mt_TT[ih]->Fill(sigMT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
					toy_HT_TT[ih]->Fill(HT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				}
			}
		} 

	std::vector<double> toyvec; 
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin));
		p_PhoEt->SetBinError(ibin, totalerror);
		// total error : stat + syst error
	}
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_LepPt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_LepPt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_LepPt->GetBinError(ibin)*p_LepPt->GetBinError(ibin));
		p_LepPt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_MET[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET->GetBinError(ibin)*p_MET->GetBinError(ibin));
		p_MET->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_Mt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt->GetBinError(ibin)*p_Mt->GetBinError(ibin));
		p_Mt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_HT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT->GetBinError(ibin)*p_HT->GetBinError(ibin));
		p_HT->SetBinError(ibin, totalerror);
	}
		
	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_PhoEt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt_TT->GetBinError(ibin)*p_PhoEt_TT->GetBinError(ibin));
		p_PhoEt_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_MET_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET_TT->GetBinError(ibin)*p_MET_TT->GetBinError(ibin));
		p_MET_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_HT_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT_TT->GetBinError(ibin)*p_HT_TT->GetBinError(ibin));
		p_HT_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)
			toyvec.push_back(toy_Mt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt_TT->GetBinError(ibin)*p_Mt_TT->GetBinError(ibin));
		p_Mt_TT->SetBinError(ibin, totalerror);
	}
	std::ostringstream outputname;
	outputname << "/eos/uscms/store/user/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_jetbkg";
	else if(channelType==2)outputname << "mg_jetbkg";
	if(anatype==0 or anatype==1) outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
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
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	for(unsigned it(0); it < NTOY; it++){
	//	toy_PhoEt[it]->Write();
	//	toy_MET[it]->Write();
	//	toy_Mt[it]->Write();
	//	toy_HT[it]->Write();
	//	toy_LepPt[it]->Write();
		toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();

}


