#include "pred_weight.h"

void pred_rareBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
	bkgType bgtype = bkgType::rare;

	init_weight( ichannel, bgtype); 

  gSystem->Load("../../lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;
	//*********** histo list **********************//
	TH1F *prefit_PhoEt = (TH1F*)file_prefit->Get("p_PhoEt");
	TH1F *prefit_MET   = (TH1F*)file_prefit->Get("p_MET");
	TH1F *prefit_HT  = (TH1F*)file_prefit->Get("p_HT");

	TH1D *postfit_PhoEt = new TH1D("postfit_PhoEt","; postfit_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *postfit_MET = new TH1D("postfit_MET","; postfit_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *postfit_HT = new TH1D("postfit_HT","HT; H_{T} (GeV);",nSigHTBins, sigHTBins); 

	TProfile *profile_PhoEt = new TProfile("profile_PhoEt","; profile_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TProfile *profile_MET = new TProfile("profile_MET","; profile_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TProfile *profile_HT = new TProfile("profile_HT","HT; H_{T} (GeV);",nSigHTBins, sigHTBins); 

	TH1D *h_rare_norm;
	if(channelType==1){
		h_rare_norm            = new TH1D("eg_rare_norm","eventcount",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_rare_norm            = new TH1D("mg_rare_norm","eventcount",NBIN,0,NBIN);
	} 

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
  mctree->SetBranchAddress("mcPID",     &mcPID);
  mctree->SetBranchAddress("mcEta",     &mcEta);
  mctree->SetBranchAddress("mcPhi",     &mcPhi);
  mctree->SetBranchAddress("mcPt",      &mcPt);
  mctree->SetBranchAddress("mcMomPID",  &mcMomPID);
  mctree->SetBranchAddress("mcStatus",  &mcStatus);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);

		double scalefactor = getWeight(ichannel, bgtype, phoEt, nVertex, phoEta, lepPt, lepEta); 
		double XS_weight = getEvtWeight(RunYear,crosssection,ntotalevent);
		double weight = PUweight*XS_weight*scalefactor;

		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

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

	
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			weight = weight*fit_norm[SigBinIndex];	
		}	

		postfit_MET->Fill(sigMET, weight);
		if(sigMET <= 120)profile_MET->Fill(sigMET, 1);
		else profile_MET->Fill(sigMET, fit_error[SigBinIndex]); 
		if(sigMET > lowMET){
			postfit_PhoEt->Fill(phoEt,weight);
			postfit_HT->Fill(HT, weight);
			profile_PhoEt->Fill(phoEt, fit_error[SigBinIndex]);
			profile_HT->Fill(HT, fit_error[SigBinIndex]);

			if(SigBinIndex >=0){
				h_rare_norm->Fill( SigBinIndex, weight);
			}	
		}

	}

	for(int ibin(1); ibin < postfit_PhoEt->GetSize(); ibin++){
		double totalerror = prefit_PhoEt->GetBinError(ibin)*profile_PhoEt->GetBinContent(ibin); 
		postfit_PhoEt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < postfit_MET->GetSize(); ibin++){
		double totalerror = prefit_MET->GetBinError(ibin)*profile_MET->GetBinContent(ibin); 
		postfit_MET->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < postfit_HT->GetSize(); ibin++){
		double totalerror = prefit_HT->GetBinError(ibin)*profile_HT->GetBinContent(ibin); 
		postfit_HT->SetBinError(ibin, totalerror);
	}


	std::ostringstream outputname;
	outputname << "postTree_";
	if(channelType==1)outputname << "egamma_rareBkg";
	else if(channelType==2)outputname << "mg_rareBkg";
	outputname << ".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	postfit_PhoEt->Write();
	postfit_MET->Write();
	postfit_HT->Write();

	outputfile->Write();
	outputfile->Close();
}


