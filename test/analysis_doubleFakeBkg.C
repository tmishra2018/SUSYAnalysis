#include "../../include/analysis_commoncode.h"
void analysis_doubleFakeBkg(){

	SetRunConfig();
	setTDRStyle();

	gSystem->Load("../../lib/libAnaClasses.so");
	int channelType = ichannel; // eg = 1; mg = 2;

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	bool toDeriveScale(false);
        if(anatype == 0)toDeriveScale = true;

	// lepton fake transfer factor (same as analysis_qcdBkg.C)
	double factorQCD(1);

        if(toDeriveScale)
                factorQCD = 1;
        else{
                if(channelType == 1){
                        if(RunYear==2016 and preVFP == 1)
                                factorQCD = factor_egQCD_2016preVFP;
                        if(RunYear==2016 and preVFP == 0)
                                factorQCD = factor_egQCD_2016postVFP;
                        if(RunYear==2017 and preVFP == 0)
                                factorQCD = factor_egQCD_2017;
                        if(RunYear==2018 and preVFP == 0)
                                factorQCD = factor_egQCD_2018;
                }
                else if(channelType == 2){
                        if(RunYear==2016 and preVFP == 1)
                                factorQCD = factor_mgQCD_2016preVFP;
                        if(RunYear==2016 and preVFP == 0)
                                factorQCD = factor_mgQCD_2016postVFP;
                        if(RunYear==2017 and preVFP == 0)
                                factorQCD = factor_mgQCD_2017;
                        if(RunYear==2018 and preVFP == 0)
                                factorQCD = factor_mgQCD_2018;
                }
        }
	// electron->photon fake transfer factor (same as analysis_eleBkg.C)
	
	float correction_ePho =1.0;
        if (RunYear == 2016 and preVFP == 1 and ichannel == 1) correction_ePho = 0.948005;
        else if (RunYear == 2016 and preVFP == 1 and ichannel == 2) correction_ePho = 0.900015;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 1) correction_ePho = 0.775942;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 2) correction_ePho = 0.730433;
        else if (RunYear == 2017 and ichannel == 1) correction_ePho = 1.39352;
        else if (RunYear == 2017 and ichannel == 2) correction_ePho = 1.03558;
        else if (RunYear == 2018 and ichannel == 1) correction_ePho = 1.35192;
        else if (RunYear == 2018 and ichannel == 2) correction_ePho = 1.20607;

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
	h_nominal_fakerate.SetParameters(scalefactor*correction_ePho, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	// jet->photon fake transfer factor (same as analysis_jetBkg.C)
	float correction_jetPho = 1.0;
        if (RunYear == 2016 and preVFP == 1 and ichannel == 1) correction_jetPho = 0.816828;
        else if (RunYear == 2016 and preVFP == 1 and ichannel == 2) correction_jetPho = 1.11658;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 1) correction_jetPho = 0.788482;
        else if (RunYear == 2016 and preVFP == 0 and ichannel == 2) correction_jetPho = 1.03158;
        else if (RunYear == 2017 and ichannel == 1) correction_jetPho = 0.866308;
        else if (RunYear == 2017 and ichannel == 2) correction_jetPho = 1.00739;
        else if (RunYear == 2018 and ichannel == 1) correction_jetPho = 0.713093;
        else if (RunYear == 2018 and ichannel == 2) correction_jetPho = 0.976136;

	TF1 *fitfunc_num = new TF1("fitfunc_num",Exp2c_Func,35,1000,6);
	TF1 *fitfunc_den = new TF1("fitfunc_den",Exp2c_Func,35,1000,6);
	
	std::stringstream JetFakeRateFile;
  	JetFakeRateFile.str();
	if(channelType==1)JetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-DoubleEG-EB.txt";
	if(channelType==2)JetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-MuonEG-EB.txt";
	std::ifstream jetfakefile(JetFakeRateFile.str().c_str());
	std::string paratype;
	float paravalue;	
	for(int i(0); i < 6; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_den->SetParameter(i, paravalue);
	}
	for(int i(0); i < 6; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_num->SetParameter(i, paravalue);
	}
	jetfakefile.close();

	//*********** histo list **********************//
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

	TH1D *p_elePhoPart = new TH1D("p_elePhoPart","ele-to-pho piece",1,0,1);
	TH1D *p_jetPhoPart = new TH1D("p_jetPhoPart","jet-to-pho piece",1,0,1);

	//************ Double Fake Tree **********************//
	TChain *doubleFaketree = new TChain("doubleFakeTree");
	if(channelType==1)doubleFaketree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	if(channelType==2)doubleFaketree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	
	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float lepMiniIso(0);
	int   lepIsStandardProxy(0);
	int   phoType(0);
	float sigMT(0);
	float sigMET(0);
	float sigMETPhi(0);
	float dPhiLepMET(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float HT(0);
	float nJet(0);
	int   nBJet(0);

	doubleFaketree->SetBranchAddress("phoEt",     &phoEt);
	doubleFaketree->SetBranchAddress("phoEta",    &phoEta);
	doubleFaketree->SetBranchAddress("phoPhi",    &phoPhi);
	doubleFaketree->SetBranchAddress("lepPt",     &lepPt);
	doubleFaketree->SetBranchAddress("lepEta",    &lepEta);
	doubleFaketree->SetBranchAddress("lepPhi",    &lepPhi);
	doubleFaketree->SetBranchAddress("lepMiniIso",&lepMiniIso);
	doubleFaketree->SetBranchAddress("lepIsStandardProxy", &lepIsStandardProxy);
	doubleFaketree->SetBranchAddress("phoType",   &phoType);
	doubleFaketree->SetBranchAddress("sigMT",     &sigMT);
	doubleFaketree->SetBranchAddress("sigMET",    &sigMET);
	doubleFaketree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	doubleFaketree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	doubleFaketree->SetBranchAddress("nVertex",   &nVertex);
	doubleFaketree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	doubleFaketree->SetBranchAddress("HT",        &HT);
	doubleFaketree->SetBranchAddress("nJet",      &nJet);
	doubleFaketree->SetBranchAddress("nBJet",     &nBJet);

	for(unsigned ievt(0); ievt < doubleFaketree->GetEntries(); ievt++){
		doubleFaketree->GetEntry(ievt);
		p_PU->Fill(nVertex);

		// same analysis-region cuts as other background macros
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
                if(sigMET < lowMET)continue;
                if(highMET > 0 && sigMET > highMET)continue;
                if(sigMT < lowMt)continue;
                if(highMt > 0 && sigMT > highMt)continue;
                if(lepPt < lowPt)continue;
                if(highPt > 0 && lepPt > highPt)continue;
                // different MET, MT, lepton pT ranges
                bool isProxy(false);
                // lepIso = 4
		if(channelType==1){if(lepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((lepMiniIso > 0.2 && lepMiniIso < lepIso*0.1))isProxy=true;}
		if(lepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;

		double w_lep(0);
		w_lep = factorQCD;

		double w_pho(0);
		if(phoType == 0){
			w_pho = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
			p_elePhoPart->Fill(0.5, w_lep*w_pho);
		}
		else if(phoType == 1){
			w_pho = correction_jetPho*fitfunc_num->Eval(phoEt)/fitfunc_den->Eval(phoEt);
			p_jetPhoPart->Fill(0.5, w_lep*w_pho);
		}
		else continue;

		double w_double = w_lep*w_pho;
		if(w_double < 0)continue;

		p_PhoEt->Fill(phoEt, w_double);
		p_PhoEta->Fill(phoEta, w_double);
		p_LepPt->Fill(lepPt, w_double);
		p_LepEta->Fill(lepEta, w_double);
		p_MET->Fill(sigMET, w_double);
		p_Mt->Fill(sigMT, w_double);
		p_HT->Fill(HT, w_double);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_double);
		p_nJet->Fill(nJet, w_double);
		p_nBJet->Fill(nBJet, w_double);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  w_double);
			p_MET_TT->Fill(sigMET,  w_double);
			p_Mt_TT->Fill(sigMT,  w_double);
			p_HT_TT->Fill(HT,  w_double);
		}
	}

	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_doubleFake";
	else if(channelType==2)outputname << "mg_doubleFake";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
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
	p_elePhoPart->Write();
	p_jetPhoPart->Write();
	outputfile->Write();
	outputfile->Close();
}
