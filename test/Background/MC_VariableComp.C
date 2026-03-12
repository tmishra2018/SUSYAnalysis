#include "../../include/analysis_commoncode.h"

void MC_VariableComp(){

   TString channel;
   Double_t plotEtBins[]={0,50,100,150,200,250,300,800};
   TString process_Trib[11] = {"WGToLNuG", "WGJet40", "WGJet130","ZGToLLG","DYJetsToLL","TTGJets", "WWG", "WZG", "WW", "WZ","TTJets"};
   TString process_Mengl[11] = {"WG35_VetoEle", "WG50_VetoEle", "WG130_VetoEle", "ZG_VetoEle", "DY", "TTG_VetoEle", "WWG_VetoEle", "WZG_VetoEle", "WW_VetoEle", "WZ_VetoEle", "TT_VetoEle"};

   Float_t         phoEt;
   Float_t         phoEta;
   Float_t         lepPt;
   Float_t         lepEta;
   Float_t         sigMT;
   Float_t         sigMET;
   Float_t         HT;
   Float_t         nJet;
   Int_t           nISRJet;
   Float_t         nISRJet_;
   Float_t         ISRJetPt;
   Float_t         dPhiLepMET;

for (int ichannel=0; ichannel<=1; ichannel++){
   if(ichannel==0) channel = "egamma";
   if(ichannel==1) channel = "mgamma";
   for (int i=0; i<=10; i++)
   {
   	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
   	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
   	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
   	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
   	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
   	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
 	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins);
   	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
   	TH1D *p_nISRJet = new TH1D("p_nISRJet","p_nISRJet",10,0,10);
   	TH1F *p_ISRJetPt     = new TH1F("p_ISRJetPt","",7,plotEtBins);
   	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2);
	
	TFile file_Trib("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_"+process_Trib[i]+"_2016postVFP.root");
	//TFile file_Trib("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_"+process_Trib[i]+".root");
   	TTree *egtree = (TTree*)file_Trib.Get("egTree");

	TFile *file;
	if (ichannel==0 and i ==0) file = new TFile("Simulation_Variables_Tri_2016postVFP.root", "RECREATE");
	else file = new TFile("Simulation_Variables_Tri_2016postVFP.root", "UPDATE");
	
	if (ichannel == 1){
		egtree->SetBranchStatus("*",1);
		egtree->SetBranchAddress("phoEt", &phoEt);
   		egtree->SetBranchAddress("phoEta", &phoEta);
   		egtree->SetBranchAddress("lepPt", &lepPt);
   		egtree->SetBranchAddress("lepEta", &lepEta);
   		egtree->SetBranchAddress("sigMT", &sigMT);
   		egtree->SetBranchAddress("sigMET", &sigMET);
   		egtree->SetBranchAddress("HT", &HT);
   		egtree->SetBranchAddress("nJet", &nJet);
   		egtree->SetBranchAddress("nISRJet", &nISRJet);
   		egtree->SetBranchAddress("ISRJetPt", &ISRJetPt);
   		egtree->SetBranchAddress("dPhiLepMET", &dPhiLepMET);
        	
		for(unsigned ievt(0); ievt<egtree->GetEntries(); ++ievt){//loop on entries
		 	egtree->GetEntry(ievt);
		 	p_PhoEt->Fill(phoEt);
		 	p_PhoEta->Fill(phoEta);
		 	p_LepPt->Fill(lepPt);
		 	p_LepEta->Fill(lepEta);
		 	p_Mt->Fill(sigMT);
		 	p_MET->Fill(sigMET);
		 	p_HT->Fill(HT);
		 	p_nJet->Fill(nJet);
		 	p_nISRJet->Fill(nISRJet);
		 	p_ISRJetPt->Fill(ISRJetPt);
		 	p_dPhiEleMET->Fill(dPhiLepMET);
		}
   		if (i ==0 ) TDirectory *chan = file->mkdir(channel);
		TDirectory *proces = file->mkdir(channel+"/"+process_Trib[i]);
		TDirectory *Dir = proces->GetDirectory(process_Trib[i]);
    		Dir->cd(); 
     		p_PhoEt->Write();
     		p_PhoEta->Write();
      		p_LepPt->Write();
      		p_LepEta->Write();
      		p_Mt->Write();
      		p_MET->Write();
      		p_HT->Write();
      		p_nJet->Write();
      		p_nISRJet->Write();
      		p_ISRJetPt->Write();
      		p_dPhiEleMET->Write();
	}
	else {
   		TTree *mgtree = (TTree*)file_Trib.Get("mgTree");

		mgtree->SetBranchStatus("*",1);
		mgtree->SetBranchAddress("phoEt", &phoEt);
   		mgtree->SetBranchAddress("phoEta", &phoEta);
   		mgtree->SetBranchAddress("lepPt", &lepPt);
   		mgtree->SetBranchAddress("lepEta", &lepEta);
   		mgtree->SetBranchAddress("sigMT", &sigMT);
   		mgtree->SetBranchAddress("sigMET", &sigMET);
   		mgtree->SetBranchAddress("HT", &HT);
   		mgtree->SetBranchAddress("nJet", &nJet);
   		mgtree->SetBranchAddress("nISRJet", &nISRJet_);
   		mgtree->SetBranchAddress("ISRJetPt", &ISRJetPt);
   		mgtree->SetBranchAddress("dPhiLepMET", &dPhiLepMET);
        	
		for(unsigned ievt(0); ievt<mgtree->GetEntries(); ++ievt){//loop on entries
		 	mgtree->GetEntry(ievt);
		 	p_PhoEt->Fill(phoEt);
		 	p_PhoEta->Fill(phoEta);
		 	p_LepPt->Fill(lepPt);
		 	p_LepEta->Fill(lepEta);
		 	p_Mt->Fill(sigMT);
		 	p_MET->Fill(sigMET);
		 	p_HT->Fill(HT);
		 	p_nJet->Fill(nJet);
		 	p_nISRJet->Fill(nISRJet_);
		 	p_ISRJetPt->Fill(ISRJetPt);
		 	p_dPhiEleMET->Fill(dPhiLepMET);
		}
   		if (i ==0) TDirectory *chan = file->mkdir(channel);
		TDirectory *proces = file->mkdir(channel+"/"+process_Trib[i]);
		TDirectory *Dir = proces->GetDirectory(process_Trib[i]);
    		Dir->cd(); 
     		p_PhoEt->Write();
     		p_PhoEta->Write();
      		p_LepPt->Write();
      		p_LepEta->Write();
      		p_Mt->Write();
      		p_MET->Write();
      		p_HT->Write();
      		p_nJet->Write();
      		p_nISRJet->Write();
      		p_ISRJetPt->Write();
      		p_dPhiEleMET->Write();
	}
	
      	file->Close();
      }}

	
}
