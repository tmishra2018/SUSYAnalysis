#include "pred_weight.h"

void pred_qcdBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
	bkgType bgtype = bkgType::qcdfakelep;

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

	TH1D *h_qcdfakelep_norm;
	if(channelType==1){
		h_qcdfakelep_norm            = new TH1D("eg_qcdfakelep_norm","eventcount",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_qcdfakelep_norm            = new TH1D("mg_qcdfakelep_norm","eventcount",NBIN,0,NBIN);
	} 
	
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	//if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	//if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/Combination/resTree_egsignal_DoubleEG-test.root");
	if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/Combination/resTree_mgsignal_MuonEG-test.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
  float sigMT(0);
  float sigMET(0);
  int   nVertex(0);
  float HT(0);
  
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("HT",        &HT);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);

		double w_qcd = getWeight(ichannel, bgtype, phoEt, nVertex, phoEta, lepPt, lepEta); 
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
	
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			w_qcd = w_qcd*fit_norm[SigBinIndex];	
		}	

		postfit_MET->Fill(sigMET, w_qcd);
		if(sigMET <= 120)profile_MET->Fill(sigMET, 1);
		else profile_MET->Fill(sigMET, fit_error[SigBinIndex]); 
		if(sigMET > lowMET){
			postfit_PhoEt->Fill(phoEt,w_qcd);
			postfit_HT->Fill(HT, w_qcd);
			profile_PhoEt->Fill(phoEt, fit_error[SigBinIndex]);
			profile_HT->Fill(HT, fit_error[SigBinIndex]);

			if(SigBinIndex >=0){
				h_qcdfakelep_norm->Fill( SigBinIndex, w_qcd);
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
	if(channelType==1)outputname << "egamma_qcd";
	else if(channelType==2)outputname << "mg_qcd";
	outputname << ".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	postfit_PhoEt->Write();
	postfit_MET->Write();
	postfit_HT->Write();

	outputfile->Write();
	outputfile->Close();
}


