#include "../../include/analysis_commoncode.h"

#define NTOY 10000
bool useGaussFit=true;

void fill_jetBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double randomweight_jet[NTOY];
	randomweight_jet[0] = 0;

	TH1D *gau=new TH1D("gau","",100,-2,2);
	for(unsigned ir(1); ir<NTOY; ir++){
		randomweight_jet[ir] = gRandom->Gaus(0, 1);
		gau->Fill(randomweight_jet[ir]);
	}
	gau->Draw();
	

	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,1000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,1000,4);
	TF1 *fitfunc_num_alt = new TF1("fitfunc_num_alt",jetfake_func,35,1000,4);
	TF1 *fitfunc_den_alt = new TF1("fitfunc_den_alt",jetfake_func,35,1000,4);
	double jetfake_numerror[265];
	double jetfake_denerror[265];
	
	std::stringstream JetFakeRateFile;
  JetFakeRateFile.str();
	//if(channelType==1)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-DoubleEG-EB-15.txt";
	//if(channelType==2)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-MuonEG-EB-15.txt";
	if(channelType==1)JetFakeRateFile << "../script/JetFakeRate-transferfactor-DoubleEG-EB.txt";
	if(channelType==2)JetFakeRateFile << "../script/JetFakeRate-transferfactor-MuonEG-EB.txt";
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
	if(channelType==1)AltJetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-DoubleEG-EB-5.txt";
	if(channelType==2)AltJetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-MuonEG-EB-129.txt";
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
	double n_PhoEt[nSigEtBins] = {0};
	double n_LepPt[nSigPtBins] = {0};
	double n_MET[nSigMETBins]  = {0};
	double n_Mt[nSigMtBins]    = {0};
	double n_HT[nSigHTBins]    = {0};
	double n_controlPhoEt[nSigEtBins] = {0};
	double n_controlLepPt[nSigPtBins] = {0};
	double n_controlMET[nSigMETBins]  = {0};
	double n_controlMt[nSigMtBins]    = {0};
	double n_controlHT[nSigHTBins]    = {0};
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	double n_norminal[NBIN]      = {0};
	double n_controlsample[NBIN] = {0};
	double n_transferfactor[NBIN]= {0};
	double e_stat[NBIN]          = {-1};
	double e_jes[NBIN]           = {-1};
	double e_phoesf[NBIN]        = {-1};
	double e_lepesf[NBIN]        = {-1};
	double e_scale[NBIN]         = {-1};
	double e_e_to_pho[NBIN]      = {-1};
	double e_j_to_pho[NBIN]      = {-1};
	double e_j_to_lep[NBIN]      = {-1};
	double e_xs[NBIN]            = {-1};
	double e_lumi[NBIN]          = {-1};
	double e_isr[NBIN]           = {-1};
	double e_total_up[NBIN]      = {-1};
	double e_total_down[NBIN]    = {-1};

	double toy_norminal[NBIN][NTOY]    = {0}; 
	double toy_PhoEt[nSigEtBins][NTOY] = {0};
	double toy_LepPt[nSigPtBins][NTOY] = {0};
	double toy_MET[nSigMETBins][NTOY]  = {0};
	double toy_Mt[nSigMtBins][NTOY]    = {0};
	double toy_HT[nSigHTBins][NTOY]    = {0};
	TH1D *toy_dPhiEleMET[NTOY];
	std::ostringstream histname;
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}

	double n_MET_low[nSigMETBins]  = {0};
	double n_controlMET_low[nSigMETBins]  = {0};
	double toy_MET_low[nSigMETBins][NTOY]  = {0};
	double n_MET_high[nSigMETBins]  = {0};
	double n_controlMET_high[nSigMETBins]  = {0};
	double toy_MET_high[nSigMETBins][NTOY]  = {0};
	/************ jet tree **************************/ 
		TChain *jettree = new TChain("jetTree");
		if(channelType==1)jettree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal_newEta.root");
		if(channelType==2)jettree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
		//if(channelType==1)jettree->Add("/uscms_data/d3/mengleis/Combination/resTree_egsignal_DoubleEG-test.root");
		//if(channelType==2)jettree->Add("/uscms_data/d3/mengleis/Combination/resTree_mgsignal_MuonEG-test.root");

		float phoEt(0);
		float phoEta(0);
		float phoChIso(0);
		float lepPt(0);
		float lepEta(0);
		float sigMT(0);
		float sigMET(0);
		float dPhiLepMET(0);
		float HT(0);
	
		jettree->SetBranchAddress("phoEt",     &phoEt);
		jettree->SetBranchAddress("phoEta",    &phoEta);
		//jettree->SetBranchAddress("phoChIso",  &phoChIso);
		jettree->SetBranchAddress("lepPt",     &lepPt);
		jettree->SetBranchAddress("lepEta",    &lepEta);
		jettree->SetBranchAddress("sigMT",     &sigMT);
		jettree->SetBranchAddress("sigMET",    &sigMET);
		jettree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
		jettree->SetBranchAddress("HT",        &HT);

			for(int ipt(0); ipt < 264; ipt++){
				double w_jet = fitfunc_num->Eval(35+ipt)/fitfunc_den->Eval(35+ipt);
				double jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt] + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(35+ipt);
				std::cout << "err% " << 35+ipt << " " << jetfakeerror/w_jet << std::endl;
			}
	 
		for (unsigned ievt(0); ievt<jettree->GetEntries(); ++ievt){//loop on entries
			jettree->GetEntry(ievt);
			/** cut flow *****/
			if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
			if(sigMT < lowMt)continue;
			if(highMt > 0 && sigMT > highMt)continue;
			if(lepPt < lowPt)continue;
			if(highPt > 0 && lepPt > highPt)continue;

			int binPhoEt = findIndex(sigEtBins, phoEt, nSigEtBins);
			int binLepPt = findIndex(sigPtBins, lepPt, nSigPtBins);
			int binMET   = findIndex(sigMETBins,sigMET,nSigMETBins);
			int binMt    = findIndex(sigMtBins, sigMT, nSigMtBins);
			int binHT    = findIndex(sigHTBins, HT,    nSigHTBins); 
 
			if(binPhoEt < 0 || binLepPt < 0 || binMET < 0 || binMt < 0 || binHT < 0){
				std::cout << "variable below threshold" << std::endl;
				exit (EXIT_FAILURE);
			}

			double evalphoEt(0);
			if(phoEt < 299)evalphoEt = phoEt;
			else evalphoEt = 299;
			double w_jet(0);
			w_jet = fitfunc_num->Eval(evalphoEt)/fitfunc_den->Eval(evalphoEt);
			double jetfakeerror(0);
			for(int ipt(0); ipt < 264; ipt++){
				if(evalphoEt >= ipt+35 && evalphoEt <= ipt+1+35)jetfakeerror = sqrt(jetfake_numerror[ipt]*jetfake_numerror[ipt] + jetfake_denerror[ipt]*jetfake_denerror[ipt]*w_jet*w_jet)/fitfunc_den->Eval(evalphoEt);
			}
			if(jetfakeerror/w_jet < 0.001)std::cout << phoEt << " relerr " << jetfakeerror/w_jet << std::endl;
			double sysJetFakePho = jetfakeerror/w_jet;
			double w_toy[NTOY] = {0};	
			for(unsigned ih(0); ih<NTOY; ih++)w_toy[ih] = w_jet*(1+sysJetFakePho*randomweight_jet[ih]);

			n_MET[binMET] += w_jet;
			n_controlMET[binMET] += 1;
			for(unsigned it(0); it < NTOY; it++)toy_MET[binMET][it] += w_toy[it];
			if(phoEt < 200){
				n_MET_low[binMET] += w_jet;
				n_controlMET_low[binMET] += 1;
				for(unsigned it(0); it < NTOY; it++)toy_MET_low[binMET][it] += w_toy[it];
			}
			else if(phoEt >= 200){
				n_MET_high[binMET] += w_jet;
				n_controlMET_high[binMET] += 1;
				for(unsigned it(0); it < NTOY; it++)toy_MET_high[binMET][it] += w_toy[it];
			}

			/** cut flow *****/
			if(sigMET < lowMET)continue;
			if(highMET > 0 && sigMET > highMET)continue;

			n_PhoEt[binPhoEt] += w_jet;
			n_LepPt[binLepPt] += w_jet;
			n_Mt[binMt]       += w_jet;
			n_HT[binHT]       += w_jet;
			n_controlPhoEt[binPhoEt] += 1;
			n_controlLepPt[binLepPt] += 1;
			n_controlMt[binMt]       += 1;
			n_controlHT[binHT]       += 1;
			p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_jet);
			for(unsigned it(0); it < NTOY; it++){
				toy_PhoEt[binPhoEt][it] 	+= w_toy[it]; 
				toy_LepPt[binLepPt][it] 	+= w_toy[it];
				toy_Mt[binMt][it]    	+= w_toy[it];
				toy_HT[binHT][it]    	+= w_toy[it];
				toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET),w_toy[it]);
			}

			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt);
			if(SigBinIndex >=0){
				n_norminal[SigBinIndex]        += w_jet;
				n_controlsample[SigBinIndex]   += 1;
				for(unsigned it(0); it < NTOY; it++)toy_norminal[SigBinIndex][it] += w_toy[it];
			}			

		} 

	//*****  Writing Estimation ********//
	std::ostringstream outputname;
	if(anatype > 0){
		switch(anatype){
			case 1: outputname << "bkgTree_";break;	
			case 2: outputname << "validTree_"; break;
			case 3: outputname << "signalTree_"; break;
		}
		if(channelType==1)outputname << "egamma_jetBkg.dat";
		else if(channelType==2)outputname << "mg_jetBkg.dat";
		
		FILE *outputfile = fopen(outputname.str().c_str(),"w"); 

		std::vector<double> toyvec; 
	
		double ave_rate[2] = {0}, ave_relerr[2] = {0};
		for(unsigned ir(0); ir < 2; ir++){ 
			double totalN(0);
			double totalCR(0); 
			toyvec.clear();
			for(unsigned ibin(0); ibin < 9; ibin++){
				totalN += n_norminal[ibin + ir*9];
				totalCR+= n_controlsample[ibin + ir*9];
			}
			toyvec.push_back(totalN/totalCR);
			ave_rate[ir] = totalN/totalCR;
	
			for(unsigned it(0); it < NTOY; it++){
				totalN = 0;
				for(unsigned ibin(0); ibin < 9; ibin++)totalN += toy_norminal[ibin + ir*9][it];
				toyvec.push_back( totalN/totalCR ); 
			}
			double syserr = calcToyError( toyvec, useGaussFit, 100);
			
			ave_relerr[ir] = syserr/ave_rate[ir];
		}

		for(int ibin(0); ibin < NBIN; ibin++){
			if(n_norminal[ibin] > 0){
				toyvec.clear();
				toyvec.push_back(n_norminal[ibin]);
				for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_norminal[ibin][it]);
				if(ibin == 11)
					for(unsigned it(0); it < NTOY; it++)std::cout << "toy " << it << toy_norminal[ibin][it] << std::endl;
				double syserr = calcToyError( toyvec, useGaussFit, ibin);
				double staterrU = getPosErrUp(n_controlsample[ibin])*(n_norminal[ibin]/n_controlsample[ibin]);
				double staterrD = getPosErrDown(n_controlsample[ibin])*(n_norminal[ibin]/n_controlsample[ibin]);
	
				n_transferfactor[ibin] = n_norminal[ibin]/n_controlsample[ibin];
				e_stat[ibin] = sqrt(n_controlsample[ibin])/n_controlsample[ibin];
				e_j_to_pho[ibin] = syserr/n_norminal[ibin];
				e_total_up[ibin] = sqrt(syserr*syserr + staterrU*staterrU);
				e_total_down[ibin] = sqrt(syserr*syserr + staterrD*staterrD);
			}
			else{	
				n_transferfactor[ibin] = ave_rate[ibin/9];
				e_stat[ibin] = 0;
				e_j_to_pho[ibin] = ave_relerr[ibin/9];
				e_total_up[ibin] = getPosErrUp(0)*ave_rate[ibin/9];
				e_total_down[ibin] = getPosErrDown(0)*ave_rate[ibin/9];
			}

			fprintf(outputfile, "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
							n_norminal[ibin], n_controlsample[ibin], n_transferfactor[ibin], e_stat[ibin],  
							e_jes[ibin], e_phoesf[ibin], e_lepesf[ibin], e_scale[ibin], e_e_to_pho[ibin], e_j_to_pho[ibin], e_j_to_lep[ibin],
							e_xs[ibin],  e_lumi[ibin],e_isr[ibin],   e_total_up[ibin], e_total_down[ibin]); 
		}

//**		for(int ibin(0); ibin < nSigEtBins; ibin++){
//**			toyvec.clear();
//**			toyvec.push_back(n_PhoEt[ibin]);
//**			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[ibin][it]);
//**			double syserr = calcToyError( toyvec, useGaussFit);
//**			double staterrU = getPosErrUp(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
//**			double staterrD = getPosErrDown(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
//**			if(n_PhoEt[ibin] == 0){
//**				staterrU = getPosErrUp(n_controlPhoEt[ibin])*ave_rate[1];
//**				staterrD = getPosErrDown(n_controlPhoEt[ibin])*ave_rate[1];
//**			}
//**			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//**			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//**			fprintf(outputfile, "ET %10.3f %10.3f %10.3f\n", n_PhoEt[ibin], totalerrorU, totalerrorD);
//**		}
//**		for(int ibin(0); ibin < nSigMETBins; ibin++){
//**			toyvec.clear();
//**			toyvec.push_back(n_MET[ibin]);
//**			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET[ibin][it]);
//**			double syserr = calcToyError( toyvec, useGaussFit);
//**			double staterrU = getPosErrUp(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
//**			double staterrD = getPosErrDown(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
//**			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//**			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//**			fprintf(outputfile, "MET %10.3f %10.3f %10.3f\n", n_MET[ibin], totalerrorU, totalerrorD);
//**		}
//**		for(int ibin(0); ibin < nSigMETBins; ibin++){
//**			toyvec.clear();
//**			toyvec.push_back(n_MET_low[ibin]);
//**			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET_low[ibin][it]);
//**			double syserr = calcToyError( toyvec, useGaussFit);
//**			double staterrU = getPosErrUp(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
//**			double staterrD = getPosErrDown(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
//**			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//**			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//**			fprintf(outputfile, "MET_low %10.3f %10.3f %10.3f\n", n_MET_low[ibin], totalerrorU, totalerrorD);
//**		}
//**		for(int ibin(0); ibin < nSigMETBins; ibin++){
//**			toyvec.clear();
//**			toyvec.push_back(n_MET_high[ibin]);
//**			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET_high[ibin][it]);
//**			double syserr = calcToyError( toyvec, useGaussFit);
//**			double staterrU = getPosErrUp(n_controlMET_high[ibin])*(n_MET_high[ibin]/n_controlMET_high[ibin]);
//**			double staterrD = getPosErrDown(n_controlMET_high[ibin])*(n_MET_high[ibin]/n_controlMET_high[ibin]);
//**			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//**			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//**			fprintf(outputfile, "MET_high %10.3f %10.3f %10.3f\n", n_MET_high[ibin], totalerrorU, totalerrorD);
//**		}
//**		for(int ibin(0); ibin < nSigHTBins; ibin++){
//**			toyvec.clear();
//**			toyvec.push_back(n_HT[ibin]);
//**			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[ibin][it]);
//**			double syserr = calcToyError( toyvec, useGaussFit);
//**			double staterrU = getPosErrUp(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
//**			double staterrD = getPosErrDown(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
//**			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//**			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//**			fprintf(outputfile, "HT %10.3f %10.3f %10.3f\n", n_HT[ibin], totalerrorU, totalerrorD);
//**		}
//**	
	}
	if(anatype == 0){
		outputname << "controlTree_";
		if(channelType==1)outputname << "egamma_jetBkg";
		else if(channelType==2)outputname << "mg_jetBkg";
		outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << ".root";
		TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
		outputfile->cd();
		p_dPhiEleMET->Write();
		for(unsigned it(0); it < NTOY; it++){
			toy_dPhiEleMET[it]->Write();
		}
		outputfile->Write();
		outputfile->Close();
	}
	
}

