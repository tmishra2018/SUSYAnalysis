#include "../analysis_commoncode.h"

void analysis_VGBkg(){
	
	SetRunConfig();
	setTDRStyle();

  gSystem->Load("../../lib/libAnaClasses.so");
	
	esfScaleFactor  objectESF;
	bool toDeriveScale(false);

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
	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "NLOvalidTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_VGBkg";
	else if(channelType==2)outputname << "mg_VGBkg";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
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
	TH1D *p_ISRJet = new TH1D("p_ISRJet",";p_{T}^{ISR} (GeV);", 100,0,1000);

	TH1D *isrup_PhoEt = new TH1D("isrup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *isrup_PhoEta = new TH1D("isrup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *isrup_LepPt = new TH1D("isrup_LepPt","isrup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *isrup_LepEta = new TH1D("isrup_LepEta","isrup_LepEta",60,-3,3);
	TH1D *isrup_MET = new TH1D("isrup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *isrup_Mt = new TH1D("isrup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *isrup_HT = new TH1D("isrup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *isrup_dPhiEleMET = new TH1D("isrup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
//mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG35_VetoEle.root");
//mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG50_VetoEle.root");
//mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WG130_VetoEle.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WGPt35_NLO.root");
  //mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WGPt130_NLO.root");
	//mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_WGPt130_NLO.root");
	float crosssection(0);
	float ntotalevent(0);
	float ISRWeight(0);
	float pdfWeight(0);
	std::vector<float> *pdfSystWeight = 0;
	std::vector<float> *ScaleSystWeight = 0; 
	int   mcType(0);
	float PUweight(1);
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
	int   nBJet(0);
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
	mctree->SetBranchAddress("pdfWeight", &pdfWeight);
	mctree->SetBranchAddress("pdfSystWeight", &pdfSystWeight);
	mctree->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight);
	mctree->SetBranchAddress("mcType",    &mcType);
	mctree->SetBranchAddress("PUweight",  &PUweight);
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
	mctree->SetBranchAddress("llmass",    &llmass);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("nBJet",     &nBJet);
  mctree->SetBranchAddress("ISRJetPt",        &ISRPt);
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

	//   start filling ///
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);

		double scalefactor(0);
		double scalefactorup(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
			scalefactorup = scalefactor + s_error; 
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
			scalefactorup = scalefactor + s_error; 
		}

		float XS_weight = 35.87*1000*crosssection/ntotalevent;
		XS_weight = XS_weight*385.41*0.985719;
	//	if(ievt <= 62202 && phoEt > 140)continue;  // 117233=mg, 62202=eg
	//	else if(ievt >  62202 && phoEt <= 140)continue;
	//	if(ievt > 62202)XS_weight = XS_weight*0.0531572*1561498.0/27426836.0; 

		float weight = PUweight*XS_weight*scalefactor*ISRWeight*factorMC;
		float weight_noisr = PUweight*XS_weight*scalefactor*factorMC;

		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;

		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){mindRpho=dR; phoIndex=iMC;}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )istruepho=true;
		}
		if(!istruepho)continue;

		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta,weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta,weight);
		p_MET->Fill(sigMET, weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
		p_nJet->Fill(nJet, weight);
		p_nBJet->Fill(nBJet, weight);
		p_ISRJet->Fill(ISRPt, weight_noisr);
	
		isrup_PhoEt->Fill(phoEt, weight_noisr);
		isrup_PhoEta->Fill(phoEta,weight_noisr);
		isrup_LepPt->Fill(lepPt, weight_noisr);
		isrup_LepEta->Fill(lepEta,weight_noisr);
		isrup_MET->Fill(sigMET, weight_noisr);
		isrup_Mt->Fill(sigMT, weight_noisr);
		isrup_HT->Fill(HT, weight_noisr);
		isrup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_noisr);
	}

	outputfile->Write();
	outputfile->Close();

}


