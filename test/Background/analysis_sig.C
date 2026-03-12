#include<string>
#include "../../include/analysis_commoncode.h"

void analysis_sig(){
	
	SetRunConfig();
	setTDRStyle();

  gSystem->Load("../../lib/libAnaClasses.so");
	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

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

	TH1D *p_LepPt_TT = new TH1D("p_LepPt_TT","p_LepPt",nBkgPtBins,bkgPtBins);
        TH1D *p_nJet_TT = new TH1D("p_nJet_TT","p_nJet",10,0,10);
        TH1D *p_nBJet_TT = new TH1D("p_nBJet_TT","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET_TT = new TH1D("p_dPhiEleMET_TT","dPhiEleMET",32,0,3.2);

	const int nHTSlices = 3;
        double HTsliceMin[nHTSlices] = {0.,   100., 400.};
        double HTsliceMax[nHTSlices] = {100., 400., 1e6};
        TH2D *p_MET_vs_HTslice = new TH2D("p_MET_vs_HTslice", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(1,"HT: 0-100");
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(2,"HT: 100-400");
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(3,"HT: >400");

	const int nMETSlices = 3;
        double METsliceMin[nMETSlices] = {120., 200., 400.};
        double METsliceMax[nMETSlices] = {200., 400., 1e6};
        TH2D *p_HT_vs_METslice = new TH2D("p_HT_vs_METslice", "HT in MET slices; H_{T} (GeV); MET slice", nBkgHTBins, bkgHTBins, nMETSlices, 0, nMETSlices);
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(1,"MET: 120-200");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(2,"MET: 200-400");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(3,"MET: >400");

	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	// signatree from data
                if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
                if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

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
	float nJetFloat(0);
  	int nJetInt(0);
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
	sigtree->SetBranchAddress("nBJet",     &nBJet);
	if (channelType == 1) sigtree->SetBranchAddress("nJet", &nJetFloat);
  	else sigtree->SetBranchAddress("nJet", &nJetInt);

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
		// different MET, MT, lepton pT ranges
		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		//p_nJet->Fill(nJet);
		p_nBJet->Fill(nBJet);
		if (channelType == 1) p_nJet->Fill(nJetFloat);
                if (channelType == 2) p_nJet->Fill(nJetInt);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt);
			p_MET_TT->Fill(sigMET);
			p_Mt_TT->Fill(sigMT);
			p_HT_TT->Fill(HT);
			p_LepPt_TT->Fill(lepPt);
                        if (channelType == 1) p_nJet_TT->Fill(nJetFloat);
                        if (channelType == 2) p_nJet_TT->Fill(nJetInt);
                        p_nBJet_TT->Fill(nBJet);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET));
		}
		
		if(anatype != 0){
                        int htslice = -1;
                        if(HT < 100)        htslice = 0;
                        else if(HT < 400)   htslice = 1;
                        else                htslice = 2;
                        p_MET_vs_HTslice->Fill(sigMET, htslice + 0.5);}

		if (anatype != 0) {
                        int metslice = -1;
                        if      (sigMET >= 120 && sigMET < 200) metslice = 0;
                        else if (sigMET >= 200 && sigMET < 400) metslice = 1;
                        else if (sigMET >= 400)                 metslice = 2;
                        p_HT_vs_METslice->Fill(HT, metslice + 0.5);}
	} 

        
	TH1D *h_MET_HT0 = p_MET_vs_HTslice->ProjectionX("h_MET_HT0",1,1);
        TH1D *h_MET_HT1 = p_MET_vs_HTslice->ProjectionX("h_MET_HT1",2,2);
        TH1D *h_MET_HT2 = p_MET_vs_HTslice->ProjectionX("h_MET_HT2",3,3);
        int nTotBins = nHTSlices * nBkgMETBins;
        TH1D *h_MET_HT_concat = new TH1D("h_MET_HT_concat", "MET in HT slices; MET / HT slice index;Events", nTotBins, 0, nTotBins);
        TH1D *hists[3] = {h_MET_HT0, h_MET_HT1, h_MET_HT2};
        for(int i=0;i<3;i++){
                for(int b=1;b<=nBkgMETBins;b++){
                        int newBin = i*nBkgMETBins + b;
                        h_MET_HT_concat->SetBinContent(newBin, hists[i]->GetBinContent(b));
                        h_MET_HT_concat->SetBinError  (newBin, hists[i]->GetBinError(b));
                }
        }
        delete h_MET_HT0;
        delete h_MET_HT1;
        delete h_MET_HT2;
        delete p_MET_vs_HTslice;

	TH1D *h_HT_MET0 = p_HT_vs_METslice->ProjectionX("h_HT_MET0",1,1);
        TH1D *h_HT_MET1 = p_HT_vs_METslice->ProjectionX("h_HT_MET1",2,2);
        TH1D *h_HT_MET2 = p_HT_vs_METslice->ProjectionX("h_HT_MET2",3,3);

        nTotBins = nMETSlices * nBkgHTBins;
        TH1D *h_HT_MET_concat = new TH1D("h_HT_MET_concat","HT in MET slices; HT / MET slice index; Events",nTotBins, 0, nTotBins);
        TH1D *p_HT_MET[3] = {h_HT_MET0, h_HT_MET1, h_HT_MET2};
        for (int i = 0; i < 3; i++) {
                for (int b = 1; b <= nBkgHTBins; b++) {
                        int newBin = i*nBkgHTBins + b;
                        h_HT_MET_concat->SetBinContent(newBin, p_HT_MET[i]->GetBinContent(b));
                        h_HT_MET_concat->SetBinError  (newBin, p_HT_MET[i]->GetBinError(b));
                }
        }
        delete h_HT_MET0;
        delete h_HT_MET1;
        delete h_HT_MET2;
        delete p_HT_vs_METslice;

	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_signal";
	else if(channelType==2)outputname << "mg_signal";
	if(anatype==0 or anatype==1)	outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname <<"_" << RunYear<<whichVFP <<".root";	

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	h_MET_HT_concat->Write();
	h_HT_MET_concat->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	p_nBJet->Write();
	p_PhoEt_TT->Write();
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	p_dPhiEleMET_TT->Write();
	outputfile->Write();
	outputfile->Close();
}


