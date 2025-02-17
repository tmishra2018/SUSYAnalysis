#include "pred_weight.h"

void pred_eleBkg(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
	bkgType bgtype = bkgType::elefakepho;

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

	TH1D *h_elefakepho_norm;
	if(channelType==1){
		h_elefakepho_norm            = new TH1D("eg_elefakepho_norm","eventcount",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_elefakepho_norm            = new TH1D("mg_elefakepho_norm","eventcount",NBIN,0,NBIN);
	} 

	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
	if(channelType==1)proxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal_newEta.root");
	if(channelType==2)proxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");

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
	
	proxytree->SetBranchAddress("phoEt",     &phoEt);
	proxytree->SetBranchAddress("phoEta",    &phoEta);
	proxytree->SetBranchAddress("phoPhi",    &phoPhi);
	proxytree->SetBranchAddress("lepPt",     &lepPt);
	proxytree->SetBranchAddress("lepEta",    &lepEta);
	proxytree->SetBranchAddress("lepPhi",    &lepPhi);
	proxytree->SetBranchAddress("sigMT",     &sigMT);
	proxytree->SetBranchAddress("sigMET",    &sigMET);
	proxytree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   &nVertex);
	proxytree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	proxytree->SetBranchAddress("HT",        &HT);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		double w_ele = getWeight(ichannel, bgtype, phoEt, nVertex, phoEta, lepPt, lepEta); 
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			w_ele = w_ele*fit_norm[SigBinIndex];	
		}	

		postfit_MET->Fill(sigMET, w_ele);
		if(sigMET <= 120)profile_MET->Fill(sigMET, 1);
		else profile_MET->Fill(sigMET, fit_error[SigBinIndex]); 
		if(sigMET > lowMET){
			postfit_PhoEt->Fill(phoEt,w_ele);
			postfit_HT->Fill(HT, w_ele);
			profile_PhoEt->Fill(phoEt, fit_error[SigBinIndex]);
			profile_HT->Fill(HT, fit_error[SigBinIndex]);

			if(SigBinIndex >=0){
				h_elefakepho_norm->Fill( SigBinIndex, w_ele);
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
	if(channelType==1)outputname << "egamma_eleBkg";
	else if(channelType==2)outputname << "mg_eleBkg";
	outputname << ".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	postfit_PhoEt->Write();
	postfit_MET->Write();
	postfit_HT->Write();

	outputfile->Write();
	outputfile->Close();
}


