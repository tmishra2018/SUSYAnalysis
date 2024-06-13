#include "../../include/analysis_commoncode.h"

void fill_qcdBkg(){

	bool toDeriveScale(false);

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  int channelType = ichannel; // eg = 1; mg =2;
	double factorQCD(1);
	double factorQCDUP = 1; 

	if(channelType == 1){
		factorQCD = factor_egQCD;
		factorQCDUP = factor_egQCD + factorerror_egQCD;
	}
	else if(channelType == 2){
		factorQCD = factor_mgQCD;
		factorQCDUP = factor_mgQCD + factorerror_mgQCD;
	}
	if(anatype == 0){
		toDeriveScale = true;
		factorQCD = 1;
		factorQCDUP = 1;
	}

	TFile *scaleFile;
	if(channelType == 1)scaleFile = TFile::Open("qcd_eg_scale.root");
	else if(channelType == 2)scaleFile = TFile::Open("qcd_mg_scale.root");
	TH1D *p_scale;
	if(channelType == 1)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	else if(channelType == 2)p_scale = (TH1D*)scaleFile->Get("transfer_factor");

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

	double normup_norminal[NBIN]    = {0}; 
	double normup_PhoEt[nSigEtBins] = {0};
	double normup_LepPt[nSigPtBins] = {0};
	double normup_MET[nSigMETBins]  = {0};
	double normup_Mt[nSigMtBins]    = {0};
	double normup_HT[nSigHTBins]    = {0};

	double unweight_norminal[NBIN]    = {0}; 
	double unweight_PhoEt[nSigEtBins] = {0};
	double unweight_LepPt[nSigPtBins] = {0};
	double unweight_MET[nSigMETBins]  = {0};
	double unweight_Mt[nSigMtBins]    = {0};
	double unweight_HT[nSigHTBins]    = {0};
	
	double n_MET_low[nSigMETBins]  = {0};
	double normup_MET_low[nSigMETBins]  = {0};
	double unweight_MET_low[nSigMETBins]  = {0};
	double n_controlMET_low[nSigMETBins]  = {0};
	double n_MET_high[nSigMETBins]  = {0};
	double n_controlMET_high[nSigMETBins]  = {0};
	double unweight_MET_high[nSigMETBins]  = {0};
	double normup_MET_high[nSigMETBins]  = {0};
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	//if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	//if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/Combination/resTree_egsignal_DoubleEG-test.root");
	if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/Combination/resTree_mgsignal_MuonEG-test.root");
  float phoEt(0);
  float phoEta(0);
  float lepPt(0);
  float lepEta(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
  float sigMT(0);
  float sigMET(0);
  float dPhiLepMET(0);
  float HT(0);
  
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("HT",        &HT);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		if(ievt%1000 ==0)std::cout <<"event " << ievt << std::endl;

		double w_qcd = 0; 
		double w_qcd_up = 0; 
		double w_qcd_unweight = 0;

		if(channelType == 1){
			w_qcd = factorQCD*p_scale->GetBinContent(p_scale->FindBin(lepPt));
			w_qcd_up = factorQCDUP*p_scale->GetBinContent(p_scale->FindBin(lepPt));
			w_qcd_unweight = factorQCD;
		}
		else{
			w_qcd = factorQCD;
			w_qcd_up = factorQCDUP;
			w_qcd_unweight = factorQCD;
		}
	
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		
		bool isProxy(false);
		if(channelType==1){if(fakeLepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;

		int binPhoEt = findIndex(sigEtBins, phoEt, nSigEtBins);
		int binLepPt = findIndex(sigPtBins, lepPt, nSigPtBins);
		int binMET   = findIndex(sigMETBins,sigMET,nSigMETBins);
		int binMt    = findIndex(sigMtBins, sigMT, nSigMtBins);
		int binHT    = findIndex(sigHTBins, HT,    nSigHTBins); 
 
		if(binPhoEt < 0 || binLepPt < 0 || binMET < 0 || binMt < 0 || binHT < 0){
			std::cout << "variable below threshold" << std::endl;
			exit (EXIT_FAILURE);
		}

		n_MET[binMET] += w_qcd;
		n_controlMET[binMET] += 1;
		normup_MET[binMET] += w_qcd_up;
		unweight_MET[binMET] += w_qcd_unweight;
		if(phoEt < 200){
			n_MET_low[binMET] += w_qcd;
			n_controlMET_low[binMET] += 1;
			normup_MET_low[binMET] += w_qcd_up;
			unweight_MET_low[binMET] += w_qcd_unweight;
		}
		else if(phoEt >= 200){
			n_MET_high[binMET] += w_qcd;
			n_controlMET_high[binMET] += 1;
			normup_MET_high[binMET] += w_qcd_up;
			unweight_MET_high[binMET] += w_qcd_unweight;
		}

		/** cut flow *****/
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		n_PhoEt[binPhoEt] += w_qcd;
		n_LepPt[binLepPt] += w_qcd;
		n_Mt[binMt]       += w_qcd;
		n_HT[binHT]       += w_qcd;
		n_controlPhoEt[binPhoEt] += 1;
		n_controlLepPt[binLepPt] += 1;
		n_controlMt[binMt]       += 1;
		n_controlHT[binHT]       += 1;
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd);

		normup_PhoEt[binPhoEt] += w_qcd_up;
		normup_LepPt[binLepPt] += w_qcd_up;
		normup_Mt[binMt]       += w_qcd_up;
		normup_HT[binHT]       += w_qcd_up;

		unweight_PhoEt[binPhoEt] += w_qcd_unweight;
		unweight_LepPt[binLepPt] += w_qcd_unweight;
		unweight_Mt[binMt]       += w_qcd_unweight;
		unweight_HT[binHT]       += w_qcd_unweight;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt);
		if(SigBinIndex >=0){
			n_norminal[SigBinIndex]        += w_qcd;
			n_controlsample[SigBinIndex]   += 1;
			normup_norminal[SigBinIndex]   += w_qcd_up;
			unweight_norminal[SigBinIndex] += w_qcd_unweight;
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
		if(channelType==1)outputname << "egamma_qcdBkg.dat";
		else if(channelType==2)outputname << "mg_qcdBkg.dat";
		
		FILE *outputfile = fopen(outputname.str().c_str(),"w"); 
	
		for(int ibin(0); ibin < NBIN; ibin++){
			if(n_norminal[ibin] > 0){
				double syserr = fabs(unweight_norminal[ibin] - n_norminal[ibin]); 
				double staterrU = getPosErrUp(n_controlsample[ibin])*(n_norminal[ibin]/n_controlsample[ibin]);
				double staterrD = getPosErrDown(n_controlsample[ibin])*(n_norminal[ibin]/n_controlsample[ibin]);
	
				n_transferfactor[ibin] = n_norminal[ibin]/n_controlsample[ibin];
				e_stat[ibin] = sqrt(n_controlsample[ibin])/n_controlsample[ibin];
				e_j_to_lep[ibin] = syserr/n_norminal[ibin];
				e_scale[ibin] = 0.2;
				e_total_up[ibin] = sqrt(syserr*syserr + 0.2*0.2*n_norminal[ibin]*n_norminal[ibin] + staterrU*staterrU);
				e_total_down[ibin] = sqrt(syserr*syserr + 0.2*0.2*n_norminal[ibin]*n_norminal[ibin] +staterrD*staterrD);
			}
			else{	
				n_transferfactor[ibin] = factorQCD;
				e_stat[ibin] = 0;
				e_j_to_lep[ibin] = -1;
				e_scale[ibin] = 0.2;
				e_total_up[ibin] = getPosErrUp(0)*factorQCD;
				e_total_down[ibin] = getPosErrDown(0)*factorQCD;
			}

			fprintf(outputfile, "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
							n_norminal[ibin], n_controlsample[ibin], n_transferfactor[ibin], e_stat[ibin],  
							e_jes[ibin], e_phoesf[ibin], e_lepesf[ibin], e_scale[ibin], e_e_to_pho[ibin], e_j_to_pho[ibin], e_j_to_lep[ibin],
							e_xs[ibin],  e_lumi[ibin],e_isr[ibin],   e_total_up[ibin], e_total_down[ibin]); 
		}

		for(int ibin(0); ibin < nSigEtBins; ibin++){
			double syserr(0);
			syserr += pow((normup_PhoEt[ibin] - n_PhoEt[ibin]),2);
			syserr += pow((unweight_PhoEt[ibin] - n_PhoEt[ibin]),2);
			syserr = sqrt(syserr);
			double staterrU = getPosErrUp(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
			double staterrD = getPosErrDown(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
			fprintf(outputfile, "ET %10.3f %10.3f %10.3f\n", n_PhoEt[ibin], totalerrorU, totalerrorD);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normup_MET[ibin] - n_MET[ibin]),2);
			syserr += pow((unweight_MET[ibin] - n_MET[ibin]),2);
			syserr = sqrt(syserr);
			double staterrU = getPosErrUp(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
			double staterrD = getPosErrDown(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
			fprintf(outputfile, "MET %10.3f %10.3f %10.3f\n", n_MET[ibin], totalerrorU, totalerrorD);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normup_MET_low[ibin] - n_MET_low[ibin]),2);
			syserr += pow((unweight_MET_low[ibin] - n_MET_low[ibin]),2);
			syserr = sqrt(syserr);
			double staterrU = getPosErrUp(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
			double staterrD = getPosErrDown(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
			fprintf(outputfile, "MET_low %10.3f %10.3f %10.3f\n", n_MET_low[ibin], totalerrorU, totalerrorD);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normup_MET_high[ibin] - n_MET_high[ibin]),2);
			syserr += pow((unweight_MET_high[ibin] - n_MET_high[ibin]),2);
			syserr = sqrt(syserr);
			double staterrU = getPosErrUp(n_controlMET_high[ibin])*(n_MET_high[ibin]/n_controlMET_high[ibin]);
			double staterrD = getPosErrDown(n_controlMET_high[ibin])*(n_MET_high[ibin]/n_controlMET_high[ibin]);
			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
			fprintf(outputfile, "MET_high %10.3f %10.3f %10.3f\n", n_MET_high[ibin], totalerrorU, totalerrorD);
		}
//		for(int ibin(0); ibin < nSigHTBins; ibin++){
//			toyvec.clear();
//			toyvec.push_back(n_HT[ibin]);
//			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[ibin][it]);
//			double syserr = calcToyError( toyvec, useGaussFit);
//			double staterrU = getPosErrUp(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
//			double staterrD = getPosErrDown(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
//			double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
//			double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
//			fprintf(outputfile, "HT %10.3f %10.3f %10.3f\n", n_HT[ibin], totalerrorU, totalerrorD);
//		}
	
	}
	if(anatype == 0){
		outputname << "controlTree_";
		if(channelType==1)outputname << "egamma_qcdBkg";
		else if(channelType==2)outputname << "mg_qcdBkg";
		outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << ".root";
		TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
		outputfile->cd();
		p_dPhiEleMET->Write();
		outputfile->Write();
		outputfile->Close();
	}

}


