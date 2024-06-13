#include "../../include/analysis_commoncode.h"

void fill_rareBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	esfScaleFactor  objectESF;

  int channelType = ichannel; // eg = 1; mg =2;
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
	TH1D   *h_rare_control       = new TH1D("rare_control","eventcount",NBIN,0,NBIN);
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

	double n_MET_low[nSigMETBins]  = {0};
	double n_controlMET_low[nSigMETBins]  = {0};
	double jesUp_MET_low[nSigMETBins]  = {0};
	double jesDown_MET_low[nSigMETBins]  = {0};
	double esfUp_MET_low[nSigMETBins]  = {0};
	double n_MET_high[nSigMETBins]  = {0};
	double n_controlMET_high[nSigMETBins]  = {0};
	double jesUp_MET_high[nSigMETBins]  = {0};
	double jesDown_MET_high[nSigMETBins]  = {0};
	double esfUp_MET_high[nSigMETBins]  = {0};
// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_TTG_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WWG_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WZG_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WW_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WZ_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_TT_VetoEle.root");
	float crosssection(0);
	float ntotalevent(0);
	float PUweight(1);
	int   mcType(0);
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
	int   nISRJet(0);
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
  std::vector<int> *mcStatus=0;

	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	mctree->SetBranchAddress("mcType",    &mcType);
  mctree->SetBranchAddress("phoEt",     &phoEt);
  mctree->SetBranchAddress("phoEta",    &phoEta);
  mctree->SetBranchAddress("phoPhi",    &phoPhi);
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("lepPhi",    &lepPhi);
  mctree->SetBranchAddress("sigMT",     &sigMT);
  mctree->SetBranchAddress("sigMET",    &sigMET);
  mctree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  mctree->SetBranchAddress("nVertex",   &nVertex);
  mctree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("nISRJet",   &nISRJet);
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
  mctree->SetBranchAddress("mcStatus",  &mcStatus);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(highMET > 0 && sigMET > highMET)continue;
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

		float weight = PUweight*XS_weight*esfscale;
		float weight_phoesfUp = PUweight*XS_weight*phoESFUp;
		float weight_lepesfUp = PUweight*XS_weight*lepESFUp;

		double mindRpho(0.3);
		int phoIndex(0), anyphoIndex(-1);
		int posW(0), negW(0), posZ(0);
		int decaytype(-1);
		bool phoCloseT(false);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho){mindRpho=dR; phoIndex=iMC;}
			if(dR < 0.3 && fabs((*mcPID)[iMC]) == 22)anyphoIndex = iMC;
			if(dR < 0.3 && fabs((*mcPID)[iMC]) == 5 && fabs((*mcMomPID)[iMC]) == 6)phoCloseT=true;

			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == 24 )posW = 1;
			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == -24)negW = 1;
			if( (fabs((*mcPID)[iMC]) == 11 || fabs((*mcPID)[iMC]) == 13 || fabs((*mcPID)[iMC]) == 15) && (*mcMomPID)[iMC] == 23 )posZ = 1;
		}
		if( (posW == 1 && negW ==1) || ( posW == 1 && posZ == 1) || (negW == 1 && posZ == 1))decaytype = 1; 

		bool isTruePho(false);
		if(mindRpho < 0.1){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )isTruePho = true;
		}
		else if(anyphoIndex >= 0){
			if((*mcPID)[anyphoIndex] == 22 && (fabs((*mcMomPID)[anyphoIndex]) <= 6 || fabs((*mcMomPID)[anyphoIndex]) == 21 || fabs((*mcMomPID)[anyphoIndex]) == 999 || fabs((*mcMomPID)[anyphoIndex])== 11 || fabs((*mcMomPID)[anyphoIndex])== 13 || fabs((*mcMomPID)[anyphoIndex])== 15 || fabs((*mcMomPID)[anyphoIndex])== 23 || fabs((*mcMomPID)[anyphoIndex])== 24)  )isTruePho = true;
		}

		bool isFSRPho(false);
		if(mindRpho < 0.3){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex])==11 || fabs((*mcMomPID)[phoIndex]) == 13 ||  fabs((*mcMomPID)[phoIndex]) == 15 || fabs((*mcMomPID)[phoIndex])==24 || fabs((*mcMomPID)[phoIndex])==6 ))isFSRPho = true;
			if((*mcPID)[phoIndex] == 22 &&  fabs((*mcMomPID)[phoIndex])==5 && phoCloseT )isFSRPho = true;
		}

		if(!isTruePho)continue;
		if(mcType ==11 && !isFSRPho)continue;
		if((mcType == 12 || mcType == 13) &&  decaytype < 1)continue;


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

		if(binPhoEt < 0 || binLepPt < 0 || binHT < 0){
			std::cout << "variable below threshold" << std::endl;
			std::cout << binPhoEt << " " << binLepPt << " " << binHT << std::endl;
			std::cout << lepPt << std::endl;
			exit (EXIT_FAILURE);
		}

		n_MET[binMET] += weight; 
		n_controlMET[binMET] += 1;
		if(binJESUPMET >= 0)jesUp_MET[binJESUPMET] += weight;
		if(binJESDOMET >= 0)jesDown_MET[binJESDOMET] += weight;
		esfUp_MET[binMET]  += weight_phoesfUp;
		if(phoEt < 200){
			n_MET_low[binMET] += weight;
			n_controlMET_low[binMET] += 1;
			if(binJESUPMET >= 0)jesUp_MET_low[binMET] += weight;
			if(binJESDOMET >= 0)jesDown_MET_low[binMET] += weight;
			esfUp_MET_low[binMET]  += weight_phoesfUp;
		}
		else if(phoEt >= 200){
			n_MET_high[binMET] += weight;
			n_controlMET_high[binMET] += 1;
			if(binJESUPMET >= 0)jesUp_MET_high[binMET] += weight;
			if(binJESDOMET >= 0)jesDown_MET_high[binMET] += weight;
			esfUp_MET_high[binMET]  += weight_phoesfUp;
		}

		/** cut flow *****/
		if(sigMT >= lowMt && sigMET >= lowMET){
	
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
				h_rare_control->Fill( SigBinIndex, weight);
				phoesfUp_norminal[SigBinIndex] += weight_phoesfUp;
				lepesfUp_norminal[SigBinIndex] += weight_lepesfUp;
			}		
	
			esfUp_PhoEt[binPhoEt]  += weight_phoesfUp;
			esfUp_Mt[binMt]  += weight_phoesfUp;
			esfUp_HT[binHT]  += weight_phoesfUp;
		}

		/** cut flow *****/
		if(sigMTJESup >= lowMt && sigMETJESup >= lowMET){
			if(Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)>=0)
				jesUp_norminal[Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)] += weight;
		}
		if(sigMTJESdo >= lowMt && sigMETJESdo >= lowMET){
			if(Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt)>=0)
				jesDown_norminal[Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt)] += weight;
		}

		if(sigMETJESup >= lowMET){
			if(binJESUPMt >= 0)jesUp_Mt[binJESUPMt] += weight;
			if(binJESDOMt >= 0)jesDown_Mt[binJESDOMt] += weight;
		}
		if(sigMETJESdo >= lowMET){
			if(binJESUPHT >= 0)jesUp_HT[binJESUPHT] += weight;
			if(binJESDOHT >= 0)jesDown_HT[binJESDOHT] += weight;
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
		if(channelType==1)outputname << "egamma_rareBkg.dat";
		else if(channelType==2)outputname << "mg_rareBkg.dat";

		h_rare_control->Sumw2();	
	
		FILE *outputfile = fopen(outputname.str().c_str(),"w"); 

		for(int ibin(0); ibin < NBIN; ibin++){
			double jesUperror = fabs(jesUp_norminal[ibin] - n_norminal[ibin]);
			double jesDownerror = fabs(jesDown_norminal[ibin] - n_norminal[ibin]);
			double jeserror   = max(jesUperror, jesDownerror);
			double phoesfUperror = fabs(phoesfUp_norminal[ibin] - n_norminal[ibin]);
			double lepesfUperror = fabs(lepesfUp_norminal[ibin] - n_norminal[ibin]);
			double xserror    = 0.5*n_norminal[ibin] + 0.025*n_norminal[ibin];

			double syserr = sqrt(jeserror*jeserror + phoesfUperror*phoesfUperror + lepesfUperror*lepesfUperror + xserror*xserror); 
			double staterrU = sqrt(n_controlsample[ibin])/n_controlsample[ibin]*n_norminal[ibin];
			double staterrD = sqrt(n_controlsample[ibin])/n_controlsample[ibin]*n_norminal[ibin]; 
				
			n_transferfactor[ibin]= -1;
			e_stat[ibin] = h_rare_control->GetBinError(ibin+1)/h_rare_control->GetBinContent(ibin+1);
			e_jes[ibin]           = jeserror/n_norminal[ibin];
			e_phoesf[ibin]        = phoesfUperror/n_norminal[ibin];
			e_lepesf[ibin]        = lepesfUperror/n_norminal[ibin];
			e_scale[ibin]         = -1;
			e_e_to_pho[ibin]      = -1;
			e_j_to_pho[ibin]      = -1;
			e_j_to_lep[ibin]      = -1;
			e_xs[ibin]            = 0.5;
			e_lumi[ibin]          = 0.025;
			e_isr[ibin]           = -1; 
			e_total_up[ibin] = sqrt(syserr*syserr + staterrU*staterrU);
			e_total_down[ibin] = sqrt(syserr*syserr + staterrD*staterrD);

			fprintf(outputfile, "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
							n_norminal[ibin], n_controlsample[ibin], n_transferfactor[ibin], e_stat[ibin],  
							e_jes[ibin], e_phoesf[ibin], e_lepesf[ibin], e_scale[ibin], e_e_to_pho[ibin], e_j_to_pho[ibin], e_j_to_lep[ibin],
							e_xs[ibin],  e_lumi[ibin],e_isr[ibin],   e_total_up[ibin], e_total_down[ibin]); 
		}
		for(int ibin(0); ibin < nSigEtBins; ibin++){
			double syserr(0);
			syserr += pow((0.5*n_PhoEt[ibin]),2);
			syserr += pow((esfUp_PhoEt[ibin] - n_PhoEt[ibin]),2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlPhoEt[ibin])*(n_PhoEt[ibin]/n_controlPhoEt[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "ET %10.3f %10.3f %10.3f\n", n_PhoEt[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((0.5*n_MET[ibin]),2);
			syserr += pow((esfUp_MET[ibin]  - n_MET[ibin]),2);
			double jeserror = max( fabs(jesUp_MET[ibin]  - n_MET[ibin]), fabs(jesDown_MET[ibin]  - n_MET[ibin]));
			syserr += pow( jeserror, 2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlMET[ibin])*(n_MET[ibin]/n_controlMET[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "MET %10.3f %10.3f %10.3f\n", n_MET[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((0.5*n_MET_low[ibin]),2);
			syserr += pow((esfUp_MET_low[ibin]  - n_MET_low[ibin]),2);
			double jeserror = max( fabs(jesUp_MET_low[ibin]  - n_MET_low[ibin]), fabs(jesDown_MET_low[ibin]  - n_MET_low[ibin]));
			syserr += pow( jeserror, 2);
			syserr = sqrt(syserr);
			double staterr = sqrt(n_controlMET_low[ibin])*(n_MET_low[ibin]/n_controlMET_low[ibin]);
			double totalerror = sqrt(syserr*syserr + staterr*staterr);
			fprintf(outputfile, "MET_low %10.3f %10.3f %10.3f\n", n_MET_low[ibin], totalerror, totalerror);
		}
		for(int ibin(0); ibin < nSigMETBins; ibin++){
			double syserr(0);
			syserr += pow((0.5*n_MET_high[ibin]),2);
			syserr += pow((esfUp_MET_high[ibin]  - n_MET_high[ibin]),2);
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
		if(channelType==1)outputname << "egamma_rareBkg";
		else if(channelType==2)outputname << "mg_rareBkg";
		outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << ".root";
		TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
		outputfile->cd();
		p_dPhiEleMET->Write();
		outputfile->Write();
		outputfile->Close();
	}


}


