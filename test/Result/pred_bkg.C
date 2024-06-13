// MC events
// check binConfig.txt and SigConfig.txt
// g++ `root-config --cflags` ../../lib/libAnaClasses.so pred_bkg.C -o pred_bkg.exe `root-config --libs`
#include "../../include/analysis_commoncode.h"


void pred_bkg(int RunYear, const char *Sample){

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;
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

	TH1D *p_eventcountData = new TH1D("p_eventcountData","eventcountData",NBIN,0,NBIN);
	//************ Signal Tree **********************//
	
	cout<<Sample<<endl;
		
	TChain *sigtreeData = new TChain("signalTree");
	//if(channelType==1) sigtreeData->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==1) sigtreeData->Add("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016.root");
	if(channelType==2) sigtreeData->Add("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016_NEW.root");
	
	//TChain *sigtree = new TChain("egTree");
	TChain *sigtree = new TChain("signalTree");
	if(channelType==1)
		//efakepho
		//sigtree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_%s_%d.root",Sample,RunYear));
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_%s_%d.root",Sample,RunYear));

        if(channelType==2)
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_%s_%d.root",Sample,RunYear));
	
  	int   run(0);
  	float ntotalevent(0);
	float crosssection(0);
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
	float nJet(0);
	float lumiWeight(1);

	float phoEta_(0), lepPt_(0),phoEt_(0), HT_(0), sigMET_(0), sigMT_(0);

 	sigtree->SetBranchAddress("run",       &run);
  	sigtree->SetBranchAddress("ntotalevent",&ntotalevent);
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
	sigtree->SetBranchAddress("nJet",      &nJet);
	sigtree->SetBranchAddress("crosssection",  &crosssection);

	sigtreeData->SetBranchAddress("phoEt",     &phoEt_);
	sigtreeData->SetBranchAddress("phoEta",    &phoEta_);
	sigtreeData->SetBranchAddress("lepPt",     &lepPt_);
	sigtreeData->SetBranchAddress("sigMT",     &sigMT_);
	sigtreeData->SetBranchAddress("sigMET",    &sigMET_);
  	sigtreeData->SetBranchAddress("HT",        &HT_);



	float puWeight = 1;
	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		lumiWeight = 1.0;
		//lumiWeight = 35.87*1000*crosssection/ntotalevent;
		//puWeight = getPUESF16(nVertex);
		puWeight = 1.0;
		//if(ievt==0) cout<<"XSec "<<crosssection<<"\t Lumi Wt "<<lumiWeight<<"\t pu weight "<<puWeight<<endl;
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
		p_nJet->Fill(nJet);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		//cout<<lumiWeight<<endl;
		if(SigBinIndex >=0)p_eventcount->Fill( SigBinIndex+0.5, lumiWeight*puWeight);
	}
	//cout<<"events lies in SR "<<p_PhoEt->Integral()<<endl;
	for (unsigned ievt(0); ievt<sigtreeData->GetEntries(); ++ievt){//loop on entries
		sigtreeData->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_ < 35 || fabs(phoEta_) > 1.4442)continue;
		if(sigMT_ < lowMt)continue;
		if(highMt > 0 && sigMT_ > highMt)continue;
		if(lepPt_ < lowPt)continue;
		if(highPt > 0 && lepPt_ > highPt)continue;
		if(sigMET_ < lowMET)continue;
		if(highMET > 0 && sigMET_ > highMET)continue;
		int SigBinIndex_(-1);
		SigBinIndex_ = Bin.findSignalBin(sigMET_, HT_, phoEt_); 
		if(SigBinIndex_ >=0)p_eventcountData->Fill(SigBinIndex_+0.5, 1);
	}

//	p_eventcount->Scale(p_eventcountData->Integral()/p_eventcount->Integral());


        

	std::ostringstream outputname;
  	outputname << "/uscms_data/d3/tmishra/Output/";
	anatype = 1;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_signal";
	else if(channelType==2)outputname << "mg_signal";
	if(anatype == 1)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << ".root";

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

	//cout<<"bkgs   "<<p_eventcount->Integral()<<" data evts  "<<p_eventcountData->Integral()<<" "<<p_eventcount->Integral()/p_eventcountData->Integral()<<endl;
	for(int ibin(1); ibin <= NBIN; ibin++)std::cout << ibin << " " << "bin " << p_eventcount->GetBinContent(ibin) << "\t"<<1/sqrt(p_eventcount->GetBinContent(ibin))<<std::endl;
	//for(int ibin(1); ibin <= NBIN; ibin++)std::cout << ibin << " " << "bin " << p_eventcount->GetBinContent(ibin) << "\t"<<p_eventcountData->GetBinContent(ibin)<<std::endl;
	outputfile->Write();
	outputfile->Close();
}

int main(int argc, char** argv)
{
    if(argc < 2)
      cout << "You have to provide two arguments!!\n";
    pred_bkg(atoi(argv[1]),argv[2]);
    return 0;
}
