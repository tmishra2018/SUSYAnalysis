#include "../../include/analysis_commoncode.h"
#define NTOY 1000
bool useGaussFit;

struct runinfo{
	int runN;
	int lumiN;
	Long_t eventN;
	int binN;
};

bool compareByRun(const runinfo &a, const runinfo &b)
{
		if(a.runN < b.runN)
			return true;
		else if(a.runN == b.runN)
			return a.lumiN < b.lumiN;
		else return false;
}

void pred_jetBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;
	
	std::string whichVFP;
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
	else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";
	
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double randomweight_jet[1000];
	randomweight_jet[0] = 0;
	for(unsigned ir(1); ir<1000; ir++)	
		randomweight_jet[ir] = gRandom->Gaus(0, 0.333);

	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,10000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,10000,4);
	TF1 *fitfunc_num_alt = new TF1("fitfunc_num_alt",jetfake_func,35,10000,4);
	TF1 *fitfunc_den_alt = new TF1("fitfunc_den_alt",jetfake_func,35,10000,4);
	double jetfake_numerror[265];
	double jetfake_denerror[265];
	
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
		std::cout << paratype << " " << paravalue << std::endl;
	}
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_num->SetParameter(i, paravalue);
		std::cout << paratype << " " << paravalue << std::endl;
	}
	int binnumber;
	for(int i(0); i < 265; i++){
		jetfakefile >> paratype >> binnumber >> paravalue;
		jetfake_numerror[i] = paravalue/2;
		std::cout << paratype << " " << paravalue << std::endl;
	}
	for(int i(0); i < 265; i++){
		jetfakefile >> paratype >> binnumber >>  paravalue;
		jetfake_denerror[i] = paravalue/2;
		std::cout << paratype << " " << paravalue << std::endl;
	}

	std::stringstream AltJetFakeRateFile;
  	AltJetFakeRateFile.str();
	if(channelType==1)AltJetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-DoubleEG-EB.txt";
	if(channelType==2)AltJetFakeRateFile << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-transferfactor-MuonEG-EB.txt";
	std::ifstream Altjetfakefile(AltJetFakeRateFile.str().c_str());
	for(int i(0); i < 4; i++){
		Altjetfakefile >> paratype >> paravalue;
		fitfunc_den_alt->SetParameter(i, paravalue);
	}
	for(int i(0); i < 4; i++){
		Altjetfakefile >> paratype >> paravalue;
		fitfunc_num_alt->SetParameter(i, paravalue);
	}



	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","; p_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nSigPtBins,sigPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	Double_t sysMETBin[] = {0,50,100,200,300,1000};
	Double_t sysHTBin[] = {0,100,200,350,500,1000};
	Double_t sysETBin[] = {0,50,100,200,1000};
	TH1D *p_MET_nom = new TH1D("p_MET_nom","p_MET_nom",5, sysMETBin);
	TH1D *p_MET_stat = new TH1D("p_MET_stat","",5, sysMETBin);
	TH1D *p_MET_alt = new TH1D("p_MET_alt","", 5, sysMETBin);
	TH1D *p_HT_nom = new TH1D("p_HT_nom","",5, sysHTBin);
	TH1D *p_HT_stat = new TH1D("p_HT_stat","",5, sysHTBin);
	TH1D *p_HT_alt = new TH1D("p_HT_alt","", 5, sysHTBin);
	TH1D *p_ET_nom = new TH1D("p_ET_nom","",4, sysETBin);
	TH1D *p_ET_stat = new TH1D("p_ET_stat","",4, sysETBin);
	TH1D *p_ET_alt = new TH1D("p_ET_alt","",4, sysETBin);
	
	std::vector<runinfo> sig_runV;
	sig_runV.clear();
	
	double predict_lowEt(0), predict_highEt(0), proxy_lowEt(0), proxy_highEt(0);
	double toy_predict_lowEt[NTOY];
	double toy_predict_highEt[NTOY];
	for(unsigned i(0); i<NTOY; i++){
		toy_predict_lowEt[i] = 0;
	  toy_predict_highEt[i] = 0;	
	}

	TH1D *h_jetfakepho_norm = 0;
	TH1D *h_jetfakepho_controlsample = 0;
	TH1D *h_jetfakepho_transferfactor = 0;
	TH1D *h_jetfakepho_syserr_jes = 0;
	TH1D *h_jetfakepho_syserr_jer = 0;
	TH1D *h_jetfakepho_syserr_esf = 0;
	TH1D *h_jetfakepho_syserr_scale = 0;
	TH1D *h_jetfakepho_syserr_e_to_pho = 0;
	TH1D *h_jetfakepho_syserr_j_to_pho = 0;
	TH1D *h_jetfakepho_syserr_j_to_lep = 0;
	TH1D *h_jetfakepho_syserr_xs = 0;
	TH1D *h_jetfakepho_syserr_lumi = 0;
	TH1D *h_jetfakepho_syserr_isr = 0;
	if(channelType==1){
		h_jetfakepho_norm            = new TH1D("eg_jetfakepho_norm","eventcount",NBIN,0,NBIN);
		h_jetfakepho_controlsample   = new TH1D("eg_jetfakepho_controlsample","",NBIN,0,NBIN);
		h_jetfakepho_transferfactor  = new TH1D("eg_jetfakepho_transferfactor","",NBIN,0,NBIN);
		h_jetfakepho_syserr_jes      = new TH1D("eg_jetfakepho_syserr_jes","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_jer      = new TH1D("eg_jetfakepho_syserr_jer","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_esf      = new TH1D("eg_jetfakepho_syserr_esf","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_scale    = new TH1D("eg_jetfakepho_syserr_scale","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_e_to_pho = new TH1D("eg_jetfakepho_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_j_to_pho = new TH1D("eg_jetfakepho_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_j_to_lep = new TH1D("eg_jetfakepho_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_xs       = new TH1D("eg_jetfakepho_syserr_xs","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_lumi     = new TH1D("eg_jetfakepho_syserr_lumi","",NBIN,0,NBIN);
		h_jetfakepho_syserr_isr      = new TH1D("eg_jetfakepho_syserr_isr","",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_jetfakepho_norm            = new TH1D("mg_jetfakepho_norm","eventcount",NBIN,0,NBIN);
		h_jetfakepho_controlsample   = new TH1D("mg_jetfakepho_controlsample","",NBIN,0,NBIN);
		h_jetfakepho_transferfactor  = new TH1D("mg_jetfakepho_transferfactor","",NBIN,0,NBIN);
		h_jetfakepho_syserr_jes      = new TH1D("mg_jetfakepho_syserr_jes","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_jer      = new TH1D("mg_jetfakepho_syserr_jer","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_esf      = new TH1D("mg_jetfakepho_syserr_esf","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_scale    = new TH1D("mg_jetfakepho_syserr_scale","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_e_to_pho = new TH1D("mg_jetfakepho_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_j_to_pho = new TH1D("mg_jetfakepho_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_j_to_lep = new TH1D("mg_jetfakepho_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_xs       = new TH1D("mg_jetfakepho_syserr_xs","",NBIN,0,NBIN);	
		h_jetfakepho_syserr_lumi     = new TH1D("mg_jetfakepho_syserr_lumi","",NBIN,0,NBIN);
		h_jetfakepho_syserr_isr      = new TH1D("mg_jetfakepho_syserr_isr","",NBIN,0,NBIN);
	} 
	
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
		toy_PhoEt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigEtBins,sigEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigPtBins,sigPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigMETBins, sigMETBins);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigMtBins,sigMtBins);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigHTBins, sigHTBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
		histname.str("");
		histname << "toy_eventcount_" << ih;
		toy_eventcount[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),NBIN,0,NBIN);
	}

		/************ jet tree **************************/ 
		TChain *jettree = new TChain("jetTree");
		if(channelType==1)
		jettree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
                if(channelType==2)
		jettree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s_Muon20.root",RunYear,whichVFP.c_str()));

  		int   run(0);
  		Long64_t  event(0);
  		int   lumis(0);
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
		float nJetFloat(0);
                int nJetInt(0);


		//jettree->SetBranchAddress("run",       &run);	
		//jettree->SetBranchAddress("event",     &event);
		//jettree->SetBranchAddress("lumis",     &lumis);
		jettree->SetBranchAddress("phoEt",     &phoEt);
		jettree->SetBranchAddress("phoEta",    &phoEta);
		jettree->SetBranchAddress("phoPhi",    &phoPhi);
		//jettree->SetBranchAddress("phoChIso",  &phoChIso);
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
		if (channelType == 1) jettree->SetBranchAddress("nJet", &nJetFloat);
                else jettree->SetBranchAddress("nJet", &nJetInt);

		for(unsigned i(0); i<264; i++){
			double w_jet = fitfunc_num->Eval(35+i)/fitfunc_den->Eval(35+i);
			double jetfakeerror = sqrt(jetfake_numerror[i]*jetfake_numerror[i] + jetfake_denerror[i]*jetfake_denerror[i]*w_jet*w_jet)/fitfunc_den->Eval(35+i);
			double sysJetFakePho = jetfakeerror/w_jet;	
			std::cout << "sys " << sysJetFakePho << " num " << jetfake_numerror[i]/fitfunc_den->Eval(35+i) << "  den " << jetfake_denerror[i]*w_jet/fitfunc_den->Eval(35+i) <<  std::endl;
		}

                for(unsigned i(0); i<1000; i++){
                        double w_jet = fitfunc_num->Eval(35+i)/fitfunc_den->Eval(35+i);
                        std::cout << "Et " << 35+i << "  fakerate " << w_jet << std::endl;
                }

	 
		for (unsigned ievt(0); ievt<jettree->GetEntries(); ++ievt){
			jettree->GetEntry(ievt);
			if (channelType == 1 && nJetFloat <1 ) continue; // NEW
                        if (channelType == 2 && nJetInt <1 ) continue;
			
			p_PU->Fill(nVertex);
			/** cut flow *****/
			if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
			if(sigMT < lowMt)continue;
			if(highMt > 0 && sigMT > highMt)continue;
			if(lepPt < lowPt)continue;
			if(highPt > 0 && lepPt > highPt)continue;

			//if(phoChIso > 20)continue;

			double evalphoEt(0);
			if(phoEt < 300)evalphoEt = phoEt;
			else evalphoEt = phoEt;
			double w_jet(0);
			w_jet = fitfunc_num->Eval(evalphoEt)/fitfunc_den->Eval(evalphoEt);

			double jetfakeerror(0);
			for(int ipt(0); ipt < 264; ipt++){
				if(evalphoEt >= ipt+35 && evalphoEt < ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt] + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(evalphoEt);
			}
			double sysJetFakePho = jetfakeerror/w_jet;	

			p_MET->Fill(sigMET,w_jet);

			p_MET_nom->Fill(sigMET,w_jet);
		//	if(phoChIso < 5)p_MET_stat->Fill(sigMET);
		//	if(phoChIso < 5)p_MET_alt->Fill(sigMET, fitfunc_num_alt->Eval(phoEt)/fitfunc_den_alt->Eval(phoEt));
		//	if(sigMT > 100 && sigMET > 120){
		//		p_HT_nom->Fill(HT, w_jet);
		//		p_ET_nom->Fill(phoEt,w_jet);
		//		if(phoChIso < 5)p_HT_stat->Fill(HT);
		//		if(phoChIso < 5)p_ET_stat->Fill(phoEt);
		//		//p_ET_stat->Fill(phoEt);
		//		if(phoChIso < 5)p_HT_alt->Fill(HT, fitfunc_num_alt->Eval(evalphoEt)/fitfunc_den_alt->Eval(evalphoEt));
		//		if(phoChIso < 5)p_ET_alt->Fill(phoEt, fitfunc_num_alt->Eval(evalphoEt)/fitfunc_den_alt->Eval(evalphoEt));
		//		std::cout << phoEt << " " << w_jet << " " << fitfunc_num_alt->Eval(evalphoEt)/fitfunc_den_alt->Eval(evalphoEt) << std::endl;
		//	}

			for(unsigned ih(0); ih<NTOY; ih++){
				toy_MET[ih]->Fill(sigMET,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
			}

			/** cut flow *****/
			if(sigMET < lowMET)continue;
			if(highMET > 0 && sigMET > highMET)continue;
	
			p_PhoEt->Fill(phoEt,w_jet);
			p_PhoEta->Fill(phoEta, w_jet);
			p_Mt->Fill(sigMT, w_jet);
			p_HT->Fill(HT, w_jet);
			p_LepPt->Fill(lepPt, w_jet);
			p_LepEta->Fill(lepEta, w_jet);
			p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_jet);
			if (channelType == 1) p_nJet->Fill(nJetFloat, w_jet);
                        if (channelType == 2) p_nJet->Fill(nJetInt, w_jet);

			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
			if(SigBinIndex >= 0){
				runinfo sig_runinfo;
				sig_runinfo.runN = run;
				sig_runinfo.eventN = event;
				sig_runinfo.lumiN = lumis;
				sig_runinfo.binN = SigBinIndex;
				sig_runV.push_back(sig_runinfo);
				h_jetfakepho_norm->Fill( SigBinIndex, w_jet);
				h_jetfakepho_controlsample->Fill( SigBinIndex );
		  	if(phoEt < 200){
		  		predict_lowEt += w_jet;
		  		proxy_lowEt   += 1;
		  	}
		  	else{
		  		predict_highEt += w_jet;
		  		proxy_highEt   += 1;
		  	}
			}	

			for(unsigned ih(0); ih<NTOY; ih++){
				
				toy_PhoEt[ih]->Fill(phoEt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_Mt[ih]->Fill(sigMT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_HT[ih]->Fill(HT, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_LepPt[ih]->Fill(lepPt,w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				toy_dPhiEleMET[ih]->Fill(fabs(dPhiLepMET), w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
				if(SigBinIndex >= 0){
					toy_eventcount[ih]->Fill(SigBinIndex, w_jet*(1+sysJetFakePho*randomweight_jet[ih]));
					if(phoEt < 200){
						toy_predict_lowEt[ih] += w_jet*(1+sysJetFakePho*randomweight_jet[ih]);
					}
					else{
						toy_predict_highEt[ih] += w_jet*(1+sysJetFakePho*randomweight_jet[ih]);
					}
				}
			}
		} 

	std::vector<double> toyvec; 
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin));
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
	
	for(int ibin(1); ibin < h_jetfakepho_norm->GetSize(); ibin++){
		if(h_jetfakepho_controlsample->GetBinContent(ibin) > 0){
			toyvec.clear();
			toyvec.push_back(h_jetfakepho_norm->GetBinContent(ibin));
			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_eventcount[it]->GetBinContent(ibin));
			double syserr = calcToyError( toyvec, useGaussFit, channelType);

			h_jetfakepho_transferfactor->SetBinContent(ibin,h_jetfakepho_norm->GetBinContent(ibin)/h_jetfakepho_controlsample->GetBinContent(ibin)); 
			h_jetfakepho_transferfactor->SetBinError(ibin,syserr/h_jetfakepho_norm->GetBinContent(ibin));
			h_jetfakepho_syserr_jes->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_jer->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_esf->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_scale->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_e_to_pho->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_j_to_pho->SetBinContent(ibin, syserr);
			h_jetfakepho_syserr_j_to_lep->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_xs->SetBinContent(ibin, -1);     
			h_jetfakepho_syserr_lumi->SetBinContent(ibin, -1);      
			h_jetfakepho_syserr_isr->SetBinContent(ibin, -1);     
		} 
		else{	
			double sf(1);
			if(ibin <= 9)sf = predict_lowEt/proxy_lowEt;
			else sf = predict_highEt/proxy_highEt;
			toyvec.clear();
			toyvec.push_back(sf);
			for(unsigned it(0); it < NTOY; it++)toyvec.push_back( ibin <= 9? toy_predict_lowEt[it]/proxy_lowEt: toy_predict_highEt[it]/proxy_highEt);
			double syserr = calcToyError( toyvec, useGaussFit, channelType);

			h_jetfakepho_transferfactor->SetBinContent(ibin, sf); 
			h_jetfakepho_transferfactor->SetBinError(ibin, syserr);
			h_jetfakepho_syserr_jes->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_jer->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_esf->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_scale->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_e_to_pho->SetBinContent(ibin, -1); 
			h_jetfakepho_syserr_j_to_pho->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_j_to_lep->SetBinContent(ibin, -1);
			h_jetfakepho_syserr_xs->SetBinContent(ibin, -1);     
			h_jetfakepho_syserr_lumi->SetBinContent(ibin, -1);      
			h_jetfakepho_syserr_isr->SetBinContent(ibin, -1);      
		}
	}
		
	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_jetbkg";
	else if(channelType==2)outputname << "mg_jetbkg";
	if(anatype ==0)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << "_" <<RunYear<<whichVFP<<".root";
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
	h_jetfakepho_norm->Write();             
	h_jetfakepho_controlsample->Write();
	h_jetfakepho_transferfactor->Write();
	h_jetfakepho_syserr_jes->Write();       
	h_jetfakepho_syserr_jer->Write();       
	h_jetfakepho_syserr_esf->Write();       
	h_jetfakepho_syserr_scale->Write();     
	h_jetfakepho_syserr_e_to_pho->Write();  
	h_jetfakepho_syserr_j_to_pho->Write();  
	h_jetfakepho_syserr_j_to_lep->Write();  
	h_jetfakepho_syserr_xs->Write();        
	h_jetfakepho_syserr_lumi->Write();      
	h_jetfakepho_syserr_isr->Write();      
	for(unsigned it(0); it < NTOY; it++){
	//	toy_PhoEt[it]->Write();
	//	toy_MET[it]->Write();
	//	toy_Mt[it]->Write();
	//	toy_HT[it]->Write();
	//	toy_LepPt[it]->Write();
	//	toy_dPhiEleMET[it]->Write();
	}
	p_MET_nom->Write();	
	p_MET_stat->Write();	
	p_MET_alt->Write();	
	p_HT_nom->Write();	
	p_HT_stat->Write();	
	p_HT_alt->Write();	
	p_ET_nom->Write();	
	p_ET_stat->Write();	
	p_ET_alt->Write();	

	std::sort(sig_runV.begin(), sig_runV.end(), compareByRun);
	for(unsigned i=0; i < sig_runV.size(); i++)std::cout <<  sig_runV[i].binN << " " <<  sig_runV[i].runN << " " << sig_runV[i].lumiN << " " << sig_runV[i].eventN << std::endl;	
	outputfile->Write();
	outputfile->Close();

}

