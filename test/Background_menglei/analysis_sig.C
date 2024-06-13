#include "../../include/analysis_commoncode.h"

void analysis_sig(){
	
	SetRunConfig();
	setTDRStyle();
        gSystem->Load("../../lib/libAnaClasses.so");

  int channelType = ichannel; // eg = 1; mg =2;
	//*********** histo list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	if(channelType==1)sigtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==2)sigtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");

	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float sigMETPhi(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float nJet(0);
	int   nBJet(0);	
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("HT",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	sigtree->SetBranchAddress("nJet",      &nJet);
	sigtree->SetBranchAddress("nBJet",     &nBJet);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		p_nJet->Fill(nJet);
		p_nBJet->Fill(nBJet);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt);
			p_MET_TT->Fill(sigMET);
			p_Mt_TT->Fill(sigMT);
			p_HT_TT->Fill(HT);
		}
	}        

	std::ostringstream outputname;

	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_signal";
	else if(channelType==2)outputname << "mg_signal";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";

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
	outputfile->Write();
	outputfile->Close();
}


