// Observed events (data)
#include "../../include/analysis_commoncode.h"
void pred_sig(){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;

        std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";

	//*********** histo list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt",";p_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nSigPtBins,sigPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_MET_low = new TH1D("p_MET_low","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
        TH1D *p_MET_high = new TH1D("p_MET_high","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; H_{T} (GeV);",nSigHTBins, sigHTBins); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",NBIN,0,NBIN);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
        if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
        if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

  	int   run(0);
  	Long64_t  event(0);
	int   lumis(0);
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

 	sigtree->SetBranchAddress("run",       &run);
  	sigtree->SetBranchAddress("event",     &event);
	sigtree->SetBranchAddress("lumis",     &lumis);
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
	if (channelType == 1) sigtree->SetBranchAddress("nJet", &nJetFloat);
        else sigtree->SetBranchAddress("nJet", &nJetInt);
	sigtree->SetBranchAddress("nBJet",     &nBJet);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){
		sigtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		p_MET->Fill(sigMET);
		if(phoEt < 200)p_MET_low->Fill(sigMET);
                else if(phoEt >= 200)p_MET_high->Fill(sigMET);

		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		if (channelType == 1) p_nJet->Fill(nJetFloat);
                if (channelType == 2) p_nJet->Fill(nJetInt);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0)p_eventcount->Fill( SigBinIndex+0.5, 1);
	}        

	std::ostringstream outputname;
  	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_signal";
	else if(channelType==2)outputname << "mg_signal";
	if(anatype ==0)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << "_" <<RunYear<<whichVFP<<".root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_MET_low->Write();
        p_MET_high->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_eventcount->Write();
	p_nJet->Write();
	for(int ibin(1); ibin <= NBIN; ibin++)std::cout << ibin << " " << "bin " << p_eventcount->GetBinContent(ibin) << std::endl;
	for(int ibin(1); ibin <= NBIN; ibin++)std::cout << p_eventcount->GetBinContent(ibin) << "  ";
	outputfile->Write();
	outputfile->Close();
}


