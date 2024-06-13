#include "../../include/analysis_commoncode.h"
#include "TProfile2D.h"

void fill_VGBkg(){
	
	bool toDeriveScale(false);
	SetSignalConfig();
	setTDRStyle();

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	esfScaleFactor  objectESF;
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);

	if(anatype == 0)toDeriveScale = true;

  int channelType = ichannel; // eg = 1; mg =2;
	double factorMC(1);
	double factorMCUP = factorMC*(1+0);
	if(toDeriveScale){
		factorMC = 1;
		factorMCUP = 1;
	}
	else{
		if(channelType == 1){
			factorMC = factor_egVGamma;
			factorMCUP = factor_egVGamma+factorerror_egVGamma;
		}
		else if(channelType == 2){
			factorMC = factor_mgVGamma;
			factorMCUP = factor_mgVGamma + factorerror_mgVGamma;
		}
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
	TH1D   *h_VG_control         = new TH1D("VG_control","eventcount",NBIN,0,NBIN);
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

	double jesUp_norminal[NBIN]    = {0}; 
	double jesUp_MET[nSigMETBins]  = {0};
	double jesUp_Mt[nSigMtBins]    = {0};
	double jesUp_HT[nSigHTBins]    = {0};
	TH1D *jesup_dPhiEleMET = new TH1D("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	double jesDown_norminal[NBIN]    = {0}; 
	double jesDown_MET[nSigMETBins]  = {0};
	double jesDown_Mt[nSigMtBins]    = {0};
	double jesDown_HT[nSigHTBins]    = {0};
	TH1D *jesDown_dPhiEleMET = new TH1D("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	double phoesfUp_norminal[NBIN] = {0}; 
	double lepesfUp_norminal[NBIN] = {0}; 
	double esfUp_PhoEt[nSigEtBins] = {0};
	double esfUp_LepPt[nSigPtBins] = {0};
	double esfUp_MET[nSigMETBins]  = {0};
	double esfUp_Mt[nSigMtBins]    = {0};
	double esfUp_HT[nSigHTBins]    = {0};

	double normUp_norminal[NBIN]    = {0}; 
	double normUp_PhoEt[nSigEtBins] = {0};
	double normUp_LepPt[nSigPtBins] = {0};
	double normUp_MET[nSigMETBins]  = {0};
	double normUp_Mt[nSigMtBins]    = {0};
	double normUp_HT[nSigHTBins]    = {0};

	double isrUp_norminal[NBIN]    = {0}; 
	double isrUp_PhoEt[nSigEtBins] = {0};
	double isrUp_LepPt[nSigPtBins] = {0};
	double isrUp_MET[nSigMETBins]  = {0};
	double isrUp_Mt[nSigMtBins]    = {0};
	double isrUp_HT[nSigHTBins]    = {0};

	double n_MET_low[nSigMETBins]  = {0};
	double n_controlMET_low[nSigMETBins]  = {0};
	double jesUp_MET_low[nSigMETBins]  = {0};
	double jesDown_MET_low[nSigMETBins]  = {0};
	double esfUp_MET_low[nSigMETBins]  = {0};
	double normUp_MET_low[nSigMETBins]  = {0};
	double isrUp_MET_low[nSigMETBins]  = {0};
	double n_MET_high[nSigMETBins]  = {0};
	double n_controlMET_high[nSigMETBins]  = {0};
	double jesUp_MET_high[nSigMETBins]  = {0};
	double jesDown_MET_high[nSigMETBins]  = {0};
	double esfUp_MET_high[nSigMETBins]  = {0};
	double normUp_MET_high[nSigMETBins]  = {0};
	double isrUp_MET_high[nSigMETBins]  = {0};
// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG50_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG130_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_ZG_VetoEle.root");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_DY.root");
	float crosssection(0);
	float ntotalevent(0);
	float ISRWeight(0);
	int   mcType(0);
	float PUweight(1);
  float phoEt(0);
  float phoEta(0);
  float lepPt(0);
  float lepEta(0);
  float sigMT(0);
  float sigMET(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float HT(0);
  float nJet(0);
	float llmass(0);
	float ISRPt(0); 
	float sigMETJESup(0);
	float sigMETJESdo(0);
	float sigMETJERup(0);
	float sigMETJERdo(0);
	float sigMTJESup(0);
	float sigMTJESdo(0);
	float sigMTJERup(0);
	float sigMTJERdo(0);
	float HTJESup(0);
	float HTJESdo(0);
	float dPhiLepMETJESup(0);
	float dPhiLepMETJESdo(0);
	float dPhiLepMETJERup(0);
	float dPhiLepMETJERdo(0);
  std::vector<int> *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
//  std::vector<int> *mcGMomPID=0;
//
	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("ISRWeight", &ISRWeight);
	mctree->SetBranchAddress("mcType",    &mcType);
	mctree->SetBranchAddress("PUweight",  &PUweight);
  mctree->SetBranchAddress("phoEt",     &phoEt);
  mctree->SetBranchAddress("phoEta",    &phoEta);
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("sigMT",     &sigMT);
  mctree->SetBranchAddress("sigMET",    &sigMET);
  mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  mctree->SetBranchAddress("nVertex",   &nVertex);
	mctree->SetBranchAddress("llmass",    &llmass);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("ISRJetPt",     &ISRPt);
	mctree->SetBranchAddress("sigMETJESup",     &sigMETJESup);
	mctree->SetBranchAddress("sigMETJESdo",     &sigMETJESdo);
	mctree->SetBranchAddress("sigMETJERup",     &sigMETJERup);
	mctree->SetBranchAddress("sigMETJERdo",     &sigMETJERdo);
	mctree->SetBranchAddress("sigMTJESup",      &sigMTJESup);
	mctree->SetBranchAddress("sigMTJESdo",      &sigMTJESdo);
	mctree->SetBranchAddress("sigMTJERup",      &sigMTJERup);
	mctree->SetBranchAddress("sigMTJERdo",      &sigMTJERdo);
	mctree->SetBranchAddress("HTJESup",         &HTJESup);
	mctree->SetBranchAddress("HTJESdo",         &HTJESdo);
	mctree->SetBranchAddress("dPhiLepMETJESup", &dPhiLepMETJESup);
	mctree->SetBranchAddress("dPhiLepMETJESdo", &dPhiLepMETJESdo);
	mctree->SetBranchAddress("dPhiLepMETJERup", &dPhiLepMETJERup);
	mctree->SetBranchAddress("dPhiLepMETJERdo", &dPhiLepMETJERdo);
  mctree->SetBranchAddress("mcPID",     &mcPID);
  mctree->SetBranchAddress("mcEta",     &mcEta);
  mctree->SetBranchAddress("mcPhi",     &mcPhi);
  mctree->SetBranchAddress("mcPt",      &mcPt);
  mctree->SetBranchAddress("mcMomPID",  &mcMomPID);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		if(mcType == 4 && llmass < 30)continue;
		if(mcType == 5 && llmass > 30)continue;
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		double phoESF(0),   lepESF(0);	
		double phoESFUp(0), lepESFUp(0);	
		double esfscale(0);
		if(channelType == 1){
			phoESF = objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			lepESF = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			esfscale = phoESF*lepESF;
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			phoESFUp = esfscale + sqrt(pow(s_pho_error,2) + pow(s_photrg_error, 2)); 
			lepESFUp = esfscale + sqrt(pow(s_ele_error,2) + pow(s_eletrg_error,2));
		}
		if(channelType == 2){
			phoESF = objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			lepESF = objectESF.getMuonESF(lepPt,lepEta);
			esfscale = phoESF*lepESF; 
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			phoESFUp = esfscale + sqrt(pow(s_pho_error,2) + pow(s_trg_error,2));
			lepESFUp = esfscale + s_mu_error; 
		}

		float XS_weight = 35.9*1000*crosssection/ntotalevent;

		float weight = PUweight*XS_weight*esfscale*ISRWeight*factorMC;
		float weight_phoesfUp = PUweight*XS_weight*phoESFUp*ISRWeight*factorMC;
		float weight_lepesfUp = PUweight*XS_weight*lepESFUp*ISRWeight*factorMC;
		float weight_normUp = PUweight*XS_weight*esfscale*ISRWeight*factorMCUP;
		float weight_noisr = PUweight*XS_weight*esfscale*factorMC;

		int binPhoEt = findIndex(sigEtBins, phoEt, nSigEtBins);
		int binLepPt = findIndex(sigPtBins, lepPt, nSigPtBins);
		int binMET   = findIndex(sigMETBins,sigMET,nSigMETBins);
		int binMt    = findIndex(sigMtBins, sigMT, nSigMtBins);
		int binHT    = findIndex(sigHTBins, HT,    nSigHTBins);
		int binJESUPMET = findIndex(sigMETBins,sigMETJESup,nSigMETBins);
		int binJESDOMET = findIndex(sigMETBins,sigMETJESdo,nSigMETBins);
		int binJESUPMt  = findIndex(sigMtBins, sigMTJESup, nSigMtBins);
		int binJESDOMt  = findIndex(sigMtBins, sigMTJESdo, nSigMtBins);
		int binJESUPHT  = findIndex(sigHTBins, HTJESup,    nSigHTBins);
		int binJESDOHT  = findIndex(sigHTBins, HTJESdo,    nSigHTBins);
 
		if(binPhoEt < 0 || binLepPt < 0 || binMET < 0 || binMt < 0 || binHT < 0){
			std::cout << "variable below threshold" << std::endl;
			std::cout << binPhoEt << " " << binLepPt << " " << binMET << " " << binMt << " " << binHT << std::endl;
			std::cout << lepPt << std::endl;
			exit (EXIT_FAILURE);
		}

		n_MET[binMET] += weight; 
		n_controlMET[binMET] += 1;
		if(binJESUPMET >= 0)jesUp_MET[binJESUPMET] += weight;
		if(binJESDOMET >= 0)jesDown_MET[binJESDOMET] += weight;
		normUp_MET[binMET] += weight_normUp;
		esfUp_MET[binMET]  += weight_phoesfUp;
		isrUp_MET[binMET]  += weight_noisr; 
		if(phoEt < 200){
			n_MET_low[binMET] += weight;
			n_controlMET_low[binMET] += 1;
			if(binJESUPMET >= 0)jesUp_MET_low[binMET] += weight;
			if(binJESDOMET >= 0)jesDown_MET_low[binMET] += weight;
			normUp_MET_low[binMET] += weight_normUp;
			esfUp_MET_low[binMET]  += weight_phoesfUp;
			isrUp_MET_low[binMET]  += weight_noisr; 
		}
		else if(phoEt >= 200){
			n_MET_high[binMET] += weight;
			n_controlMET_high[binMET] += 1;
			if(binJESUPMET >= 0)jesUp_MET_high[binMET] += weight;
			if(binJESDOMET >= 0)jesDown_MET_high[binMET] += weight;
			normUp_MET_high[binMET] += weight_normUp;
			esfUp_MET_high[binMET]  += weight_phoesfUp;
			isrUp_MET_high[binMET]  += weight_noisr; 
		}


		/** cut flow *****/
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		n_PhoEt[binPhoEt] += weight;
		n_LepPt[binLepPt] += weight;
		n_Mt[binMt]       += weight;
		n_HT[binHT]       += weight;
		n_controlPhoEt[binPhoEt] += 1;
		n_controlLepPt[binLepPt] += 1;
		n_controlMt[binMt]       += 1;
		n_controlHT[binHT]       += 1;
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			n_norminal[SigBinIndex]        += weight; 
			n_controlsample[SigBinIndex]   += 1;
			h_VG_control->Fill( SigBinIndex, weight);
			normUp_norminal[SigBinIndex]   += weight_normUp;
			phoesfUp_norminal[SigBinIndex] += weight_phoesfUp;
			lepesfUp_norminal[SigBinIndex] += weight_lepesfUp;
			isrUp_norminal[SigBinIndex]    += weight_noisr;
		}			
		if(Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)>=0)
			jesUp_norminal[Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)] += weight;
		if(Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt)>=0)
			jesDown_norminal[Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt)] += weight;

		normUp_PhoEt[binPhoEt] += weight_normUp;
		esfUp_PhoEt[binPhoEt]  += weight_phoesfUp;
		isrUp_PhoEt[binPhoEt]  += weight_noisr; 

		normUp_Mt[binMt] += weight_normUp;
		esfUp_Mt[binMt]  += weight_phoesfUp;
		isrUp_Mt[binMt]  += weight_noisr;
		if(binJESUPMt >= 0)jesUp_Mt[binJESUPMt] += weight;
		if(binJESDOMt >= 0)jesDown_Mt[binJESDOMt] += weight;

		normUp_HT[binHT] += weight_normUp;
		esfUp_HT[binHT]  += weight_phoesfUp;
		isrUp_HT[binHT]  += weight_noisr;
		if(binJESUPHT >= 0)jesUp_HT[binJESUPHT] += weight;
		if(binJESDOHT >= 0)jesDown_HT[binJESDOHT] += weight; 
	}

	//*****  Writing Estimation ********//
	std::ostringstream outputname;
	if(anatype > 0){
		switch(anatype){
			case 1: outputname << "bkgTree_";break;	
			case 2: outputname << "validTree_"; break;
			case 3: outputname << "signalTree_"; break;
		}
		if(channelType==1)outputname << "egamma_VGBkg.dat";
		else if(channelType==2)outputname << "mg_VGBkg.dat";

		h_VG_control->Sumw2();	
		
		FILE *outputfile = fopen(outputname.str().c_str(),"w"); 

		for(int ibin(0); ibin < NBIN; ibin++){
			double jesUperror = fabs(jesUp_norminal[ibin] - n_norminal[ibin]);
			double jesDownerror = fabs(jesDown_norminal[ibin] - n_norminal[ibin]);
			double jeserror   = max(jesUperror, jesDownerror);
			double phoesfUperror = fabs(phoesfUp_norminal[ibin] - n_norminal[ibin]);
			double lepesfUperror = fabs(lepesfUp_norminal[ibin] - n_norminal[ibin]);
			double normUperror= fabs(normUp_norminal[ibin] - n_norminal[ibin]);
			double isrUperror = fabs(isrUp_norminal[ibin] - n_norminal[ibin]);

			double syserr = sqrt(jeserror*jeserror + phoesfUperror*phoesfUperror + lepesfUperror*lepesfUperror +  normUperror*normUperror + isrUperror*isrUperror); 
			double staterrU = sqrt(n_controlsample[ibin])/n_controlsample[ibin]*n_norminal[ibin];
			double staterrD = sqrt(n_controlsample[ibin])/n_controlsample[ibin]*n_norminal[ibin]; 
				
			n_transferfactor[ibin]= -1;
			e_stat[ibin] = h_VG_control->GetBinError(ibin+1)/h_VG_control->GetBinContent(ibin+1);
			e_jes[ibin]           = jeserror/n_norminal[ibin];
			e_phoesf[ibin]        = phoesfUperror/n_norminal[ibin];
			e_lepesf[ibin]        = lepesfUperror/n_norminal[ibin];
			e_scale[ibin]         = 0.2;
			e_e_to_pho[ibin]      = -1;
			e_j_to_pho[ibin]      = -1;
			e_j_to_lep[ibin]      = -1;
			e_xs[ibin]            = -1;
			e_lumi[ibin]          = -1;
			e_isr[ibin]           = isrUperror/n_norminal[ibin];
			e_total_up[ibin] = sqrt(syserr*syserr + staterrU*staterrU);
			e_total_down[ibin] = sqrt(syserr*syserr + staterrD*staterrD);

			fprintf(outputfile, "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
							n_norminal[ibin], n_controlsample[ibin], n_transferfactor[ibin], e_stat[ibin],  
							e_jes[ibin], e_phoesf[ibin], e_lepesf[ibin], e_scale[ibin], e_e_to_pho[ibin], e_j_to_pho[ibin], e_j_to_lep[ibin],
							e_xs[ibin],  e_lumi[ibin],e_isr[ibin],   e_total_up[ibin], e_total_down[ibin]); 
		}

		for(int ibin(0); ibin < nSigEtBins; ibin++){
			double syserr(0);
			syserr += pow((normUp_PhoEt[ibin] - n_PhoEt[ibin]),2);
			syserr += pow((esfUp_PhoEt[ibin] - n_PhoEt[ibin]),2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "ET %10.3f %10.3f %10.3f\n", n_PhoEt[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normUp_MET[ibin] - n_MET[ibin]),2);
			syserr += pow((esfUp_MET[ibin]  - n_MET[ibin]),2);
			syserr += pow((isrUp_MET[ibin]  - n_MET[ibin]),2);
			double jeserror = max( fabs(jesUp_MET[ibin]  - n_MET[ibin]), fabs(jesDown_MET[ibin]  - n_MET[ibin]));
			syserr += pow( jeserror, 2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "MET %10.3f %10.3f %10.3f\n", n_MET[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normUp_MET_low[ibin] - n_MET_low[ibin]),2);
			syserr += pow((esfUp_MET_low[ibin]  - n_MET_low[ibin]),2);
			syserr += pow((isrUp_MET_low[ibin]  - n_MET_low[ibin]),2);
			double jeserror = max( fabs(jesUp_MET_low[ibin]  - n_MET_low[ibin]), fabs(jesDown_MET_low[ibin]  - n_MET_low[ibin]));
			syserr += pow( jeserror, 2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "MET_low %10.3f %10.3f %10.3f\n", n_MET_low[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((normUp_MET_high[ibin] - n_MET_high[ibin]),2);
			syserr += pow((esfUp_MET_high[ibin]  - n_MET_high[ibin]),2);
			syserr += pow((isrUp_MET_high[ibin]  - n_MET_high[ibin]),2);
			double jeserror = max( fabs(jesUp_MET_high[ibin]  - n_MET_high[ibin]), fabs(jesDown_MET_high[ibin]  - n_MET_high[ibin]));
			syserr += pow( jeserror, 2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlMET_high[ibin])*(n_MET_high[ibin]/n_controlMET_high[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "MET_high %10.3f %10.3f %10.3f\n", n_MET_high[ibin], totalerror, totalerror);
		}
	//	for(int ibin(0); ibin < nSigHTBins; ibin++){
	//		toyvec.clear();
	//		toyvec.push_back(n_HT[ibin]);
	//		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[ibin][it]);
	//		double syserr = calcToyError( toyvec, useGaussFit);
	//		double staterrU = getPosErrUp(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
	//		double staterrD = getPosErrDown(n_controlHT[ibin])*(n_HT[ibin]/n_controlHT[ibin]);
	//		double totalerrorU = sqrt(syserr*syserr + staterrU*staterrU);
	//		double totalerrorD = sqrt(syserr*syserr + staterrD*staterrD);
	//		fprintf(outputfile, "HT %10.3f %10.3f %10.3f\n", n_HT[ibin], totalerrorU, totalerrorD);
	//	}
	
	}
	if(anatype == 0){
		outputname << "controlTree_";
		if(channelType==1)outputname << "egamma_VGBkg";
		else if(channelType==2)outputname << "mg_VGBkg";
		outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << ".root";
		TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
		outputfile->cd();
		p_dPhiEleMET->Write();
		outputfile->Write();
		outputfile->Close();
	}
}


