#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"

#include "../../../include/analysis_commoncode.h"
#include "TProfile2D.h"
#include "../analysis_PU.C"
#include <TROOT.h>

void analysis_TChiWG(){//main  
	gROOT->SetBatch(true);
	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;
	std::string whichVFP;
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";


	float lumi = 1.0, lumi_eg = 1.0, lumi_mg = 1.0;
	if(RunYear == 2016 and preVFP == 1){
    		lumi_eg = lumi_2016preVFP_DoubleEG;
    		lumi_mg = lumi_2016preVFP_MuonEG;}
        else if(RunYear == 2016 and preVFP == 0){
    		lumi_eg = lumi_2016postVFP_DoubleEG;
    		lumi_mg = lumi_2016postVFP_MuonEG;}
        else if(RunYear == 2017){
    		lumi_eg = lumi_2017_DoubleEG;
    		lumi_mg = lumi_2017_MuonEG;}
        else if(RunYear == 2018){
    		lumi_eg = lumi_2018_DoubleEG;
    		lumi_mg = lumi_2018_MuonEG;}
        else if(RunYear == 678){
    		lumi_eg = lumi_678_DoubleEG;
    		lumi_mg = lumi_678_MuonEG;}


  	gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/lib/libAnaClasses.so");

	TFile xSecFile("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/cross/susyCrossSection.root");
	TH1D *p_crosssection_tchiwg = (TH1D*)xSecFile.Get("p_charginoSec");
	TH1D *p_crosssection_t5wg   = (TH1D*)xSecFile.Get("p_gluinoxSec");
	TH1D *p_crosssection_t6wg   = (TH1D*)xSecFile.Get("p_squarkxSec");

	TChain *datachain = new TChain("signalTree");
	datachain->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s_Muon20.root",RunYear,whichVFP.c_str()));
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  	datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	std::ostringstream histname;
	//**************   T5WG  ***************************//
	TFile *file_t5wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_%d%s.root",RunYear, whichVFP.c_str()));
  	TTree *tree_t5wg = (TTree*)file_t5wg->Get("SUSYtree");
	float Mgluino_t5wg(0);
  	float Mchargino_t5wg(0);
  	float Mneutralino_t5wg(0);
	int   nVertex(0);
	tree_t5wg->SetBranchAddress("Mgluino",    &Mgluino_t5wg);  
  	tree_t5wg->SetBranchAddress("Mchargino",  &Mchargino_t5wg);
  	tree_t5wg->SetBranchAddress("Mneutralino",&Mneutralino_t5wg);
  	tree_t5wg->SetBranchAddress("nVertex",    &nVertex);

	TFile *outputfile_t5wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_%d%s.root",RunYear, whichVFP.c_str()),"RECREATE");
	outputfile_t5wg->cd();

	TH2D *p_T5WGMASS         = new TH2D("SUSYMass","",27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_T5WGselect       = new TH2D("p_T5WGselect","",27, 775.0, 2125.0, 80, 12.5, 2012.5); 
	TH2D *p_lowPU_t5wg_pass  = new TH2D("p_lowPU_t5wg_pass","", 27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_lowPU_t5wg_all   = new TH2D("p_lowPU_t5wg_all", "", 27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_highPU_t5wg_pass = new TH2D("p_highPU_t5wg_pass","",27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_highPU_t5wg_all  = new TH2D("p_highPU_t5wg_all","", 27, 775.0, 2125.0, 80, 12.5, 2012.5);

	for(int ievt(0); ievt < tree_t5wg->GetEntries(); ievt++){
		tree_t5wg->GetEntry(ievt);
    		if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
	
		double NLSPMass(0);
    		if(Mchargino_t5wg >0)NLSPMass = Mchargino_t5wg;
		else if(Mneutralino_t5wg >0)NLSPMass = Mneutralino_t5wg;
    	
		p_T5WGMASS->Fill(Mgluino_t5wg, NLSPMass);
		if(nVertex < 20)p_lowPU_t5wg_all->Fill(Mgluino_t5wg, NLSPMass, 1);
		else p_highPU_t5wg_all->Fill(Mgluino_t5wg, NLSPMass, 1);
	}

	TH2D *t5wg_h_chan_rate_nom[NBIN*2]; 
	for(int i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		t5wg_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
	}

	TH2D *t5wg_h_chan_rate_jesUp[NBIN*2]; 
	TH2D *t5wg_h_chan_rate_jesDown[NBIN*2];    
	TH2D *t5wg_h_chan_rate_jerUp[NBIN*2];    
	TH2D *t5wg_h_chan_rate_jerDown[NBIN*2];    
	TH2D *t5wg_h_chan_rate_xsUp[NBIN*2];       
	TH2D *t5wg_h_chan_rate_esfUp[NBIN*2];       
                                  
	TH2D *t5wg_h_syserr_PU = new TH2D("t5wg_h_syserr_PU","t5wg_h_syserr_PU", 27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *t5wg_h_chan_syserr_jes[NBIN*2];      
	TH2D *t5wg_h_chan_syserr_jer[NBIN*2];     
	TH2D *t5wg_h_chan_syserr_esf[NBIN*2];     
	TH2D *t5wg_h_chan_syserr_scale[NBIN*2];   
	TH2D *t5wg_h_chan_syserr_eleshape[NBIN*2]; 
	TH2D *t5wg_h_chan_syserr_jetshape[NBIN*2];
	TH2D *t5wg_h_chan_syserr_qcdshape[NBIN*2];
	TH2D *t5wg_h_chan_syserr_xs[NBIN*2];
	TH2D *t5wg_h_chan_syserr_lumi[NBIN*2];     
	TH2D *t5wg_h_chan_syserr_isr[NBIN*2];     

	for(int i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		t5wg_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		t5wg_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		t5wg_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		t5wg_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		t5wg_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		t5wg_h_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		t5wg_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		t5wg_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		t5wg_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		t5wg_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		t5wg_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		t5wg_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		t5wg_h_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		t5wg_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		t5wg_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		t5wg_h_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
	}
	
	TH1D *p_t5wg_MET_signal_mg = new TH1D("p_t5wg_MET_signal_1700_1000_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_mg  = new TH1D("p_t5wg_HT_signal_1700_1000_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_mg = new TH1D("p_t5wg_PhoEt_signal_1700_1000_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_MET_signal_eg = new TH1D("p_t5wg_MET_signal_1700_1000_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_eg  = new TH1D("p_t5wg_HT_signal_1700_1000_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_eg = new TH1D("p_t5wg_PhoEt_signal_1700_1000_eg","",nSigEtBins, sigEtBins);
		
  	TChain *mgtree_t5wg;
  	mgtree_t5wg = new TChain("mgTree","mgTree");
	mgtree_t5wg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_%d%s.root",RunYear, whichVFP.c_str()));
  	float phoEt_t5wg_mg(0);
  	float phoEta_t5wg_mg(0);
  	float lepPt_t5wg_mg(0);
  	float lepEta_t5wg_mg(0);
  	float sigMT_t5wg_mg(0);
  	float sigMET_t5wg_mg(0);
  	float HT_t5wg_mg(0);
  	float nJet_t5wg_mg(0);
	int   nVertex_t5wg_mg(0);
	float sigMETJESup_t5wg_mg(0);
	float sigMETJESdo_t5wg_mg(0);
	float sigMETJERup_t5wg_mg(0);
	float sigMETJERdo_t5wg_mg(0);
	float sigMTJESup_t5wg_mg(0);
	float sigMTJESdo_t5wg_mg(0);
	float sigMTJERup_t5wg_mg(0);
	float sigMTJERdo_t5wg_mg(0);
	float HTJESup_t5wg_mg(0);
	float HTJESdo_t5wg_mg(0);
	float gluinoMass_t5wg_mg(0);
  	float charginoMass_t5wg_mg(0);
  	float neutralinoMass_t5wg_mg(0);
  	mgtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_mg);
  	mgtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_mg);
  	mgtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("nJet",         &nJet_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("nVertex",    &nVertex_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("Mgluino",    &gluinoMass_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("Mchargino",  &charginoMass_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("Mneutralino",&neutralinoMass_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("HTJESup",    &HTJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("HTJESdo",    &HTJESdo_t5wg_mg);

	for(int ievt(0); ievt < mgtree_t5wg->GetEntries(); ievt++){
		mgtree_t5wg->GetEntry(ievt);

		double NLSPMass(0);
		if(charginoMass_t5wg_mg > 0)NLSPMass = charginoMass_t5wg_mg;
		else if(neutralinoMass_t5wg_mg > 0)NLSPMass = neutralinoMass_t5wg_mg;
		if(NLSPMass <= 0)continue;	
	
		/** cut flow *****/
		if(nJet_t5wg_mg < 1)continue; 
		if(phoEt_t5wg_mg < 35 || lepPt_t5wg_mg < 20)continue;
		if(fabs(phoEta_t5wg_mg) > 1.4442 || fabs(lepEta_t5wg_mg) > 2.5)continue;
    p_T5WGselect->Fill(gluinoMass_t5wg_mg, NLSPMass);

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_t5wg_mg,lepEta_t5wg_mg)*objectESF.getPhotonESF(phoEt_t5wg_mg, phoEta_t5wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t5wg_mg, lepPt_t5wg_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_t5wg_mg,lepEta_t5wg_mg)*objectESF.getPhotonESF(phoEt_t5wg_mg, phoEta_t5wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t5wg_mg, lepPt_t5wg_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_t5wg_mg, phoEta_t5wg_mg)*objectESF.getFastMuonESF(lepPt_t5wg_mg,lepEta_t5wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t5wg_mg, lepPt_t5wg_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_t5wg_mg, lepPt_t5wg_mg)*objectESF.getFastMuonESF(lepPt_t5wg_mg,lepEta_t5wg_mg)*objectESF.getPhotonESF(phoEt_t5wg_mg, phoEta_t5wg_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		if( fabs(gluinoMass_t5wg_mg - 1700) < 10 && fabs(NLSPMass - 1000) < 5){
			if(sigMT_t5wg_mg > 100){
				p_t5wg_MET_signal_mg->Fill(sigMET_t5wg_mg);
				if(sigMET_t5wg_mg > 120)p_t5wg_HT_signal_mg->Fill(HT_t5wg_mg);
				if(sigMET_t5wg_mg > 120)p_t5wg_PhoEt_signal_mg->Fill(phoEt_t5wg_mg);
			}
		}

		if(sigMET_t5wg_mg < 120 || sigMT_t5wg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg);
		if(SigBinIndex >=0){
			t5wg_h_chan_rate_nom[SigBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor);
			t5wg_h_chan_rate_xsUp[SigBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor); 
			t5wg_h_chan_rate_esfUp[SigBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass,scalefactorup); 

			if(nVertex_t5wg_mg < 20)p_lowPU_t5wg_pass->Fill( gluinoMass_t5wg_mg, NLSPMass, 1);
			else p_highPU_t5wg_pass->Fill( gluinoMass_t5wg_mg, NLSPMass, 1);
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_t5wg_mg, HTJESup_t5wg_mg, phoEt_t5wg_mg); 
		if(jesupBinIndex >=0){
			t5wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_t5wg_mg, HTJESdo_t5wg_mg, phoEt_t5wg_mg);
		if(jesdoBinIndex >=0){
			t5wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg); 
		if( jerupBinIndex >=0){
			t5wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg);
		if(jerdoBinIndex >= 0){
			t5wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( gluinoMass_t5wg_mg, NLSPMass, scalefactor);
		}  
	}

  TChain *egtree_t5wg;
  egtree_t5wg = new TChain("egTree","egTree");
  egtree_t5wg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_%d%s.root",RunYear, whichVFP.c_str()));
  float phoEt_t5wg_eg(0);
  float phoEta_t5wg_eg(0);
  float lepPt_t5wg_eg(0);
  float lepEta_t5wg_eg(0);
  float sigMT_t5wg_eg(0);
  float sigMET_t5wg_eg(0);
  float HT_t5wg_eg(0);
  float nJet_t5wg_eg(0);
  int   nVertex_t5wg_eg(0); 
	float sigMETJESup_t5wg_eg(0);
	float sigMETJESdo_t5wg_eg(0);
	float sigMETJERup_t5wg_eg(0);
	float sigMETJERdo_t5wg_eg(0);
	float sigMTJESup_t5wg_eg(0);
	float sigMTJESdo_t5wg_eg(0);
	float sigMTJERup_t5wg_eg(0);
	float sigMTJERdo_t5wg_eg(0);
	float HTJESup_t5wg_eg(0);
	float HTJESdo_t5wg_eg(0);
	float gluinoMass_t5wg_eg(0);
  float charginoMass_t5wg_eg(0);
  float neutralinoMass_t5wg_eg(0);
  egtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_eg);
	egtree_t5wg->SetBranchAddress("nVertex",    &nVertex_t5wg_eg);
  egtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_eg);
  egtree_t5wg->SetBranchAddress("nJet",         &nJet_t5wg_eg);
	egtree_t5wg->SetBranchAddress("Mgluino",    &gluinoMass_t5wg_eg);
  egtree_t5wg->SetBranchAddress("Mchargino",  &charginoMass_t5wg_eg);
  egtree_t5wg->SetBranchAddress("Mneutralino",&neutralinoMass_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("HTJESup",    &HTJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("HTJESdo",    &HTJESdo_t5wg_eg);

	for(int ievt(0); ievt < egtree_t5wg->GetEntries(); ievt++){
		egtree_t5wg->GetEntry(ievt);

		double NLSPMass(0);
		if(charginoMass_t5wg_eg > 0)NLSPMass = charginoMass_t5wg_eg;
		else if(neutralinoMass_t5wg_eg > 0)NLSPMass = neutralinoMass_t5wg_eg;
		if(NLSPMass <= 0)continue;	
		/** cut flow *****/
		if(nJet_t5wg_eg < 1)continue;
		if(phoEt_t5wg_eg < 35 || lepPt_t5wg_eg < 25)continue;
		if(fabs(phoEta_t5wg_eg) > 1.4442 || fabs(lepEta_t5wg_eg) > 2.5)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getPhotonESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t5wg_eg,lepEta_t5wg_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getPhotonESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t5wg_eg,lepEta_t5wg_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronESF(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t5wg_eg,lepEta_t5wg_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getFastElectronESF(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getPhotonESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t5wg_eg,phoEta_t5wg_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronESF(lepPt_t5wg_eg,lepEta_t5wg_eg)*objectESF.getPhotonESF(phoEt_t5wg_eg,phoEta_t5wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t5wg_eg,lepEta_t5wg_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		if( fabs(gluinoMass_t5wg_eg - 1700) < 10 && fabs(NLSPMass - 1000) < 5){
			if(sigMT_t5wg_eg > 100){
				p_t5wg_MET_signal_eg->Fill(sigMET_t5wg_eg);
				if(sigMET_t5wg_eg > 120)p_t5wg_HT_signal_eg->Fill(HT_t5wg_eg);
				if(sigMET_t5wg_eg > 120)p_t5wg_PhoEt_signal_eg->Fill(phoEt_t5wg_eg);
			}
		}

		if(sigMET_t5wg_eg < 120 || sigMT_t5wg_eg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg) + NBIN;
		if(SigBinIndex >=0){
			t5wg_h_chan_rate_nom[SigBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor);
			t5wg_h_chan_rate_xsUp[SigBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor); 
			t5wg_h_chan_rate_esfUp[SigBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass,scalefactorup); 

			if(nVertex_t5wg_eg < 20)p_lowPU_t5wg_pass->Fill( gluinoMass_t5wg_eg, NLSPMass, 1);
			else p_highPU_t5wg_pass->Fill( gluinoMass_t5wg_eg, NLSPMass, 1);
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_t5wg_eg, HTJESup_t5wg_eg, phoEt_t5wg_eg) + NBIN;
		if(jesupBinIndex >=0){
			t5wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_t5wg_eg, HTJESdo_t5wg_eg, phoEt_t5wg_eg) + NBIN;
		if(jesdoBinIndex >=0){
			t5wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg) + NBIN;
		if( jerupBinIndex >=0){
			t5wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg) + NBIN;	
		if(jerdoBinIndex >= 0){
			t5wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( gluinoMass_t5wg_eg, NLSPMass, scalefactor);
		}  
	}

	for(int ih(0); ih < NBIN*2; ih++){ // first 18 binning mg, next 18 for eg
		if(ih < NBIN)   lumi = lumi_mg;
                else            lumi = lumi_eg;
		for(int i(1); i < t5wg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(int j(1); j < t5wg_h_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_T5WGMASS->GetBinContent(i,j) < 1000){
					p_T5WGselect->SetBinContent(i,j,-1);
					p_T5WGMASS->SetBinContent(i,j,-1);
				}
	
				if(p_T5WGMASS->GetBinContent(i,j) <= 0){
					t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					t5wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					t5wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					t5wg_h_syserr_PU->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
				}

				else{
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(i);
					float crosssection = p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					float crosssectionUp = (crosssection+ p_crosssection_t5wg->GetBinError( p_crosssection_t5wg->FindBin(sparticleMass) ) );

					t5wg_h_chan_rate_nom[ih]->SetBinError(i,j, sqrt(t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j))*lumi*crosssection/noe);
					t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)*lumi*crosssection/noe); 
					t5wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe); 
					t5wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t5wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t5wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t5wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)*lumi*crosssectionUp/noe);
					t5wg_h_chan_rate_esfUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
									
					t5wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(t5wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t5wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(t5wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t5wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( t5wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					t5wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j));   
					t5wg_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					if(ih==0)t5wg_h_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_t5wg_pass, p_lowPU_t5wg_all, p_highPU_t5wg_pass, p_highPU_t5wg_all, p_PU_data));  
				}
			}
		} 
	} 

		for(int i(1); i < p_T5WGMASS->GetXaxis()->GetNbins() + 1; i++){
			for(int j(1); j < p_T5WGMASS->GetYaxis()->GetNbins() + 1; j++){
				float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(i);
				if( fabs(sparticleMass - 1700 ) < 10 && fabs( p_T5WGMASS->GetYaxis()->GetBinCenter(j) - 1000) < 5){
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					p_t5wg_MET_signal_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_HT_signal_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_PhoEt_signal_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_MET_signal_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_HT_signal_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_PhoEt_signal_eg->Scale(lumi_eg*crosssection/noe);
			}
		}
	}

	outputfile_t5wg->Write();
	outputfile_t5wg->Close();


	//****************   TChiWG ***************************

/*      TFile *file_tchiwg = TFile::Open("/uscms_data/d3/mengleis/test/test_TChiWg.root");
  	TTree *tree_tchiwg = (TTree*)file_tchiwg->Get("SUSYtree");
  	float Mchargino_tchiwg(0);
  	float Mneutralino_tchiwg(0);
	int   nVertex_tchiwg(0);
  	tree_tchiwg->SetBranchAddress("Mchagino",  &Mchargino_tchiwg);
  	tree_tchiwg->SetBranchAddress("Mneutralino",&Mneutralino_tchiwg);
	tree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg);

	TFile *outputfile_tchiwg = TFile::Open("signalTree_TChiWG.root","RECREATE");
	outputfile_tchiwg->cd();

	TH1D *p_TChiWGMASS = new TH1D("p_TChiWGMASS","",40,287.5,1287.5);
	TH1D *p_lowPU_tchiwg_pass  = new TH1D("p_lowPU_tchiwg_pass","", 40,287.5,1287.5); 
	TH1D *p_lowPU_tchiwg_all   = new TH1D("p_lowPU_tchiwg_all", "", 40,287.5,1287.5);
	TH1D *p_highPU_tchiwg_pass = new TH1D("p_highPU_tchiwg_pass","",40,287.5,1287.5);
	TH1D *p_highPU_tchiwg_all  = new TH1D("p_highPU_tchiwg_all","", 40,287.5,1287.5);

	TH1D *p_tchiwg_MET_signal_eg = new TH1D("p_tchiwg_MET_signal_eg","",nSigMETBins, sigMETBins);
	TH1D *p_tchiwg_HT_signal_eg  = new TH1D("p_tchiwg_HT_signal_eg", "",nSigHTBins, sigHTBins);
	TH1D *p_tchiwg_PhoEt_signal_eg=new TH1D("p_tchiwg_PhoEt_signal_eg","",nSigEtBins, sigEtBins);
	TH1D *p_tchiwg_MET_signal_mg = new TH1D("p_tchiwg_MET_signal_mg","",nSigMETBins, sigMETBins);
	TH1D *p_tchiwg_HT_signal_mg  = new TH1D("p_tchiwg_HT_signal_mg", "",nSigHTBins, sigHTBins);
	TH1D *p_tchiwg_PhoEt_signal_mg=new TH1D("p_tchiwg_PhoEt_signal_mg","",nSigEtBins, sigEtBins);

	for(int ievt(0); ievt < tree_tchiwg->GetEntries(); ievt++){
		tree_tchiwg->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(Mchargino_tchiwg >=0){
			p_TChiWGMASS->Fill(Mchargino_tchiwg);
			if(nVertex_tchiwg < 20)p_lowPU_tchiwg_all->Fill(Mchargino_tchiwg, 1);
			else p_highPU_tchiwg_all->Fill(Mchargino_tchiwg, 1);
		}
	}

  //p_T5WGMASS->Draw();
	TH1D *tchiwg_h_chan_rate_nom[NBIN*2]; 
	TH1D *tchiwg_h_chan_rate_jesUp[NBIN*2]; 
	TH1D *tchiwg_h_chan_rate_jesDown[NBIN*2];    
	TH1D *tchiwg_h_chan_rate_jerUp[NBIN*2];    
	TH1D *tchiwg_h_chan_rate_jerDown[NBIN*2];    
	TH1D *tchiwg_h_chan_rate_xsUp[NBIN*2];       
	TH1D *tchiwg_h_chan_rate_esfUp[NBIN*2];       
                                  
	TH1D *tchiwg_h_syserr_PU = new TH1D("tchiwg_h_syserr_PU","",40,287.5,1287.5);
	TH1D *tchiwg_h_chan_syserr_jes[NBIN*2];      
	TH1D *tchiwg_h_chan_syserr_jer[NBIN*2];     
	TH1D *tchiwg_h_chan_syserr_esf[NBIN*2];     
	TH1D *tchiwg_h_chan_syserr_scale[NBIN*2];   
	TH1D *tchiwg_h_chan_syserr_eleshape[NBIN*2]; 
	TH1D *tchiwg_h_chan_syserr_jetshape[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_qcdshape[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_xs[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_lumi[NBIN*2];     
	TH1D *tchiwg_h_chan_syserr_isr[NBIN*2];     

	for(int i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		tchiwg_h_chan_rate_nom[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		tchiwg_h_chan_rate_jesUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		tchiwg_h_chan_rate_jesDown[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		tchiwg_h_chan_rate_jerUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		tchiwg_h_chan_rate_jerDown[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		tchiwg_h_chan_rate_xsUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		tchiwg_h_chan_rate_esfUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		tchiwg_h_chan_syserr_jes[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		tchiwg_h_chan_syserr_jer[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		tchiwg_h_chan_syserr_esf[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		tchiwg_h_chan_syserr_scale[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		tchiwg_h_chan_syserr_eleshape[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		tchiwg_h_chan_syserr_jetshape[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		tchiwg_h_chan_syserr_qcdshape[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		tchiwg_h_chan_syserr_xs[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		tchiwg_h_chan_syserr_lumi[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		tchiwg_h_chan_syserr_isr[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
	}

  TChain *mgtree_tchiwg;
  mgtree_tchiwg = new TChain("mgTree","mgTree");
  mgtree_tchiwg->Add("/uscms_data/d3/mengleis/test/test_TChiWg.root");
  float phoEt_tchiwg_mg(0);
  float phoEta_tchiwg_mg(0);
  float lepPt_tchiwg_mg(0);
  float lepEta_tchiwg_mg(0);
  float sigMT_tchiwg_mg(0);
  float sigMET_tchiwg_mg(0);
  float HT_tchiwg_mg(0);
	int   nVertex_tchiwg_mg(0);
	float sigMETJESup_tchiwg_mg(0);
	float sigMETJESdo_tchiwg_mg(0);
	float sigMETJERup_tchiwg_mg(0);
	float sigMETJERdo_tchiwg_mg(0);
	float sigMTJESup_tchiwg_mg(0);
	float sigMTJESdo_tchiwg_mg(0);
	float sigMTJERup_tchiwg_mg(0);
	float sigMTJERdo_tchiwg_mg(0);
	float HTJESup_tchiwg_mg(0);
	float HTJESdo_tchiwg_mg(0);
	float gluinoMass_tchiwg_mg(0);
  float charginoMass_tchiwg_mg(0);
  float neutralinoMass_tchiwg_mg(0);
  mgtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("Mgluino",    &gluinoMass_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("Mchagino",  &charginoMass_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("Mneutralino",&neutralinoMass_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJESup",&sigMETJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJERup",&sigMETJERup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJESup", &sigMTJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJERup", &sigMTJERup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("HTJESup",    &HTJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("HTJESdo",    &HTJESdo_tchiwg_mg);
		
	for(int ievt(0); ievt < mgtree_tchiwg->GetEntries(); ievt++){
		mgtree_tchiwg->GetEntry(ievt);

		// cut flow *****
		if(phoEt_tchiwg_mg < 35 || lepPt_tchiwg_mg < 20)continue;
		if(fabs(phoEta_tchiwg_mg) > 1.4442 || fabs(lepEta_tchiwg_mg) > 2.5)continue;

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_tchiwg_mg,lepEta_tchiwg_mg)*objectESF.getPhotonESF(phoEt_tchiwg_mg, phoEta_tchiwg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_tchiwg_mg, lepPt_tchiwg_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_tchiwg_mg,lepEta_tchiwg_mg)*objectESF.getPhotonESF(phoEt_tchiwg_mg, phoEta_tchiwg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_tchiwg_mg, lepPt_tchiwg_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_tchiwg_mg, phoEta_tchiwg_mg)*objectESF.getFastMuonESF(lepPt_tchiwg_mg,lepEta_tchiwg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_tchiwg_mg, lepPt_tchiwg_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_tchiwg_mg, lepPt_tchiwg_mg)*objectESF.getFastMuonESF(lepPt_tchiwg_mg,lepEta_tchiwg_mg)*objectESF.getPhotonESF(phoEt_tchiwg_mg, phoEta_tchiwg_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		if( fabs(charginoMass_tchiwg_mg - 800 ) < 2){
			if(sigMT_tchiwg_mg > 100){
			  p_tchiwg_MET_signal_mg->Fill( sigMET_tchiwg_mg );
			  if(sigMET_tchiwg_mg > 120)p_tchiwg_HT_signal_mg->Fill( HT_tchiwg_mg );
			  if(sigMET_tchiwg_mg > 120)p_tchiwg_PhoEt_signal_mg->Fill( phoEt_tchiwg_mg );
      }
		}

		if(sigMET_tchiwg_mg < 120 || sigMT_tchiwg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);
		if(SigBinIndex >=0){
			tchiwg_h_chan_rate_nom[SigBinIndex]->Fill( charginoMass_tchiwg_mg, scalefactor);
			tchiwg_h_chan_rate_xsUp[SigBinIndex]->Fill( charginoMass_tchiwg_mg, scalefactor); 
			tchiwg_h_chan_rate_esfUp[SigBinIndex]->Fill( charginoMass_tchiwg_mg,scalefactorup); 

			if(nVertex_tchiwg_mg < 20)p_lowPU_tchiwg_pass->Fill( charginoMass_tchiwg_mg, 1);
			else p_highPU_tchiwg_pass->Fill(charginoMass_tchiwg_mg, 1);
		}

		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_tchiwg_mg, HTJESup_tchiwg_mg, phoEt_tchiwg_mg);
		if(jesupBinIndex >=0){
			tchiwg_h_chan_rate_jesUp[jesupBinIndex]->Fill( charginoMass_tchiwg_mg, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_tchiwg_mg, HTJESdo_tchiwg_mg, phoEt_tchiwg_mg);
		if(jesdoBinIndex >=0){
			tchiwg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( charginoMass_tchiwg_mg, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);
		if( jerupBinIndex >=0){
			tchiwg_h_chan_rate_jerUp[jerupBinIndex]->Fill(  charginoMass_tchiwg_mg, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);	
		if(jerdoBinIndex >= 0){
			tchiwg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( charginoMass_tchiwg_mg, scalefactor);
		}  
	}


  TChain *egtree_tchiwg;
  egtree_tchiwg = new TChain("egTree","egTree");
  egtree_tchiwg->Add("/uscms_data/d3/mengleis/test/test_TChiWg.root");
  float phoEt_tchiwg_eg(0);
  float phoEta_tchiwg_eg(0);
  float lepPt_tchiwg_eg(0);
  float lepEta_tchiwg_eg(0);
  float sigMT_tchiwg_eg(0);
  float sigMET_tchiwg_eg(0);
  float HT_tchiwg_eg(0);
	int   nVertex_tchiwg_eg(0);
	float sigMETJESup_tchiwg_eg(0);
	float sigMETJESdo_tchiwg_eg(0);
	float sigMETJERup_tchiwg_eg(0);
	float sigMETJERdo_tchiwg_eg(0);
	float sigMTJESup_tchiwg_eg(0);
	float sigMTJESdo_tchiwg_eg(0);
	float sigMTJERup_tchiwg_eg(0);
	float sigMTJERdo_tchiwg_eg(0);
	float HTJESup_tchiwg_eg(0);
	float HTJESdo_tchiwg_eg(0);
	float gluinoMass_tchiwg_eg(0);
  float charginoMass_tchiwg_eg(0);
  float neutralinoMass_tchiwg_eg(0);
  egtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("Mgluino",    &gluinoMass_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("Mchagino",  &charginoMass_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("Mneutralino",&neutralinoMass_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJESup",&sigMETJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJERup",&sigMETJERup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJESup", &sigMTJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJERup", &sigMTJERup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("HTJESup",    &HTJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("HTJESdo",    &HTJESdo_tchiwg_eg);
		
	for(int ievt(0); ievt < egtree_tchiwg->GetEntries(); ievt++){
		egtree_tchiwg->GetEntry(ievt);

		//  cut flow *****
		if(phoEt_tchiwg_eg < 35 || lepPt_tchiwg_eg < 25)continue;
		if(fabs(phoEta_tchiwg_eg) > 1.4442 || fabs(lepEta_tchiwg_eg) > 2.5)continue;
		if(sigMET_tchiwg_eg < 120 || sigMT_tchiwg_eg < 100)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getPhotonESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronTRGESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getPhotonESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronTRGESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronTRGESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getFastElectronESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getPhotonESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg)*objectESF.getPhotonESF(phoEt_tchiwg_eg,phoEta_tchiwg_eg)*objectESF.getFastElectronTRGESF(lepPt_tchiwg_eg,lepEta_tchiwg_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		if( fabs(charginoMass_tchiwg_eg - 800 ) < 2){
			if(sigMT_tchiwg_eg > 100){
			  p_tchiwg_MET_signal_eg->Fill( sigMET_tchiwg_eg );
			  if(sigMET_tchiwg_eg > 120)p_tchiwg_HT_signal_eg->Fill( HT_tchiwg_eg );
			  if(sigMET_tchiwg_eg > 120)p_tchiwg_PhoEt_signal_eg->Fill( phoEt_tchiwg_eg );
      }
		}

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg)+NBIN;
		if(SigBinIndex >=0){
			tchiwg_h_chan_rate_nom[SigBinIndex]->Fill( charginoMass_tchiwg_eg, scalefactor);
			tchiwg_h_chan_rate_xsUp[SigBinIndex]->Fill( charginoMass_tchiwg_eg, scalefactor); 
			tchiwg_h_chan_rate_esfUp[SigBinIndex]->Fill( charginoMass_tchiwg_eg,scalefactorup); 

			if(nVertex_tchiwg_eg < 20)p_lowPU_tchiwg_pass->Fill( charginoMass_tchiwg_eg, 1);
			else p_highPU_tchiwg_pass->Fill(charginoMass_tchiwg_eg, 1);
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_tchiwg_eg, HTJESup_tchiwg_eg, phoEt_tchiwg_eg)+NBIN;
		if(jesupBinIndex >=0){
			tchiwg_h_chan_rate_jesUp[jesupBinIndex]->Fill( charginoMass_tchiwg_eg, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_tchiwg_eg, HTJESdo_tchiwg_eg, phoEt_tchiwg_eg)+NBIN;
		if(jesdoBinIndex >=0){
			tchiwg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( charginoMass_tchiwg_eg, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg)+NBIN;
		if( jerupBinIndex >=0){
			tchiwg_h_chan_rate_jerUp[jerupBinIndex]->Fill(  charginoMass_tchiwg_eg, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg)+NBIN;	
		if(jerdoBinIndex >= 0){
			tchiwg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( charginoMass_tchiwg_eg, scalefactor);
		}  
	}

	for(int ih(0); ih < NBIN*2; ih++){
		if(ih < NBIN)   lumi = lumi_mg;
                else            lumi = lumi_eg;
		for(int i(1); i < tchiwg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
				
			if(p_TChiWGMASS->GetBinContent(i) <= 0){
				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_rate_jesUp[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_rate_jesDown[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_rate_jerUp[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_rate_jerDown[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_rate_xsUp[ih]->SetBinContent(i, -1);
						
				tchiwg_h_chan_syserr_jes[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_jer[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_esf[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_scale[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_eleshape[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_jetshape[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_qcdshape[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_xs[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_lumi[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_isr[ih]->SetBinContent(i, -1);
        if(ih ==0)tchiwg_h_syserr_PU->SetBinContent(i, -1); 
			}
			else{ 
				float noe = p_TChiWGMASS->GetBinContent(i);
				float sparticleMass = p_TChiWGMASS->GetXaxis()->GetBinCenter(i);
				float crosssection = p_crosssection_tchiwg->GetBinContent( p_crosssection_tchiwg->FindBin(sparticleMass) );
				float crosssectionUp = (crosssection+p_crosssection_tchiwg->GetBinError( p_crosssection_tchiwg->FindBin(sparticleMass) ) ); 

				tchiwg_h_chan_rate_nom[ih]->SetBinError(i, sqrt(tchiwg_h_chan_rate_nom[ih]->GetBinContent(i))*lumi*crosssection/noe);
				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)*lumi*crosssection/noe); 
				tchiwg_h_chan_rate_jesUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_jesUp[ih]->GetBinContent(i)*lumi*crosssection/noe); 
				tchiwg_h_chan_rate_jesDown[ih]->SetBinContent(i, tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i)*lumi*crosssection/noe);
				tchiwg_h_chan_rate_jerUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_jerUp[ih]->GetBinContent(i)*lumi*crosssection/noe);
				tchiwg_h_chan_rate_jerDown[ih]->SetBinContent(i, tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i)*lumi*crosssection/noe);
				tchiwg_h_chan_rate_xsUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_xsUp[ih]->GetBinContent(i)*lumi*crosssectionUp/noe);
				tchiwg_h_chan_rate_esfUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_esfUp[ih]->GetBinContent(i)*lumi*crosssection/noe);
								
				tchiwg_h_chan_syserr_jes[ih]->SetBinContent(i, max( fabs(tchiwg_h_chan_rate_jesUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)), fabs(tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) ) ); 
				tchiwg_h_chan_syserr_jer[ih]->SetBinContent(i, max( fabs(tchiwg_h_chan_rate_jerUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)), fabs(tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) ) ); 
				tchiwg_h_chan_syserr_esf[ih]->SetBinContent(i, fabs( tchiwg_h_chan_rate_esfUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)));
				tchiwg_h_chan_syserr_scale[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_eleshape[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_jetshape[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_qcdshape[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_xs[ih]->SetBinContent(i, fabs( tchiwg_h_chan_rate_xsUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) );
				tchiwg_h_chan_syserr_lumi[ih]->SetBinContent(i, 0.026*tchiwg_h_chan_rate_nom[ih]->GetBinContent(i));   
				tchiwg_h_chan_syserr_isr[ih]->SetBinContent(i, -1); 
				if(ih==0)tchiwg_h_syserr_PU->SetBinContent(i, analysis_PU(i, p_lowPU_tchiwg_pass, p_lowPU_tchiwg_all, p_highPU_tchiwg_pass, p_highPU_tchiwg_all, p_PU_data));  
			}
		} 
	} 

		for(int i(1); i < p_TChiWGMASS->GetXaxis()->GetNbins() + 1; i++){
			float sparticleMass = p_TChiWGMASS->GetXaxis()->GetBinCenter(i);
			if(p_TChiWGMASS->GetBinContent(i) > 0 && fabs(sparticleMass - 800) < 2){
				float noe = p_TChiWGMASS->GetBinContent(i);
				float crosssection = p_crosssection_tchiwg->GetBinContent( p_crosssection_tchiwg->FindBin(sparticleMass) );
				p_tchiwg_MET_signal_eg->Scale(lumi_eg*crosssection/noe);
				p_tchiwg_HT_signal_eg->Scale(lumi_eg*crosssection/noe);
				p_tchiwg_PhoEt_signal_eg->Scale(lumi_eg*crosssection/noe);
				p_tchiwg_MET_signal_mg->Scale(lumi_mg*crosssection/noe);
				p_tchiwg_HT_signal_mg->Scale(lumi_mg*crosssection/noe);
				p_tchiwg_PhoEt_signal_mg->Scale(lumi_mg*crosssection/noe);
			}
		} 

	outputfile_tchiwg->Write();
	outputfile_tchiwg->Close();

	// ************   T6WG  ***************************

  TFile *file_t6wg = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/resTree_T6WG.root");
  TTree *tree_t6wg = (TTree*)file_t6wg->Get("SUSYtree");
	float Msquark_t6wg(0);
  float Mchargino_t6wg(0);
  float Mneutralino_t6wg(0);
	tree_t6wg->SetBranchAddress("Msquark",    &Msquark_t6wg);  
  tree_t6wg->SetBranchAddress("Mchargino",  &Mchargino_t6wg);
  tree_t6wg->SetBranchAddress("Mneutralino",&Mneutralino_t6wg);

	TFile *outputfile_t6wg = TFile::Open("signalTree_T6WG.root","RECREATE");
	outputfile_t6wg->cd();

	TH2D *p_T6WGMASS = new TH2D("T6WGSUSYMass","",27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_T6WGselect = new TH2D("p_T6WGselect","",27, 775.0, 2125.0, 80, 12.5, 2012.5); 

	for(int ievt(0); ievt < tree_t6wg->GetEntries(); ievt++){
		tree_t6wg->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;

    double NLSPMass(0);
    if(Mchargino_t6wg >0)NLSPMass=Mchargino_t6wg;
		else if(Mneutralino_t6wg >0)NLSPMass=Mneutralino_t6wg;

		if(NLSPMass <= 0)continue;
		p_T6WGMASS->Fill(Msquark_t6wg, NLSPMass);
	}
  //p_T6WGMASS->Draw();
	TH2D *t6wg_h_chan_rate_nom[NBIN*2]; 
	TH2D *t6wg_h_chan_rate_jesUp[NBIN*2]; 
	TH2D *t6wg_h_chan_rate_jesDown[NBIN*2];    
	TH2D *t6wg_h_chan_rate_jerUp[NBIN*2];    
	TH2D *t6wg_h_chan_rate_jerDown[NBIN*2];    
	TH2D *t6wg_h_chan_rate_xsUp[NBIN*2];       
	TH2D *t6wg_h_chan_rate_esfUp[NBIN*2];       
                                  
	TH2D *t6wg_h_chan_syserr_jes[NBIN*2];      
	TH2D *t6wg_h_chan_syserr_jer[NBIN*2];     
	TH2D *t6wg_h_chan_syserr_esf[NBIN*2];     
	TH2D *t6wg_h_chan_syserr_scale[NBIN*2];   
	TH2D *t6wg_h_chan_syserr_eleshape[NBIN*2]; 
	TH2D *t6wg_h_chan_syserr_jetshape[NBIN*2];
	TH2D *t6wg_h_chan_syserr_qcdshape[NBIN*2];
	TH2D *t6wg_h_chan_syserr_xs[NBIN*2];
	TH2D *t6wg_h_chan_syserr_lumi[NBIN*2];     
	TH2D *t6wg_h_chan_syserr_isr[NBIN*2];     

	for(int i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		t6wg_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		t6wg_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		t6wg_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		t6wg_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		t6wg_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		t6wg_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		t6wg_h_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		t6wg_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		t6wg_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		t6wg_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		t6wg_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		t6wg_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		t6wg_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		t6wg_h_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		t6wg_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		t6wg_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		t6wg_h_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
	}
		
  TChain *mgtree_t6wg;
  mgtree_t6wg = new TChain("mgTree","mgTree");
  mgtree_t6wg->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T6WG.root");
  float phoEt_t6wg_mg(0);
  float phoEta_t6wg_mg(0);
  float lepPt_t6wg_mg(0);
  float lepEta_t6wg_mg(0);
  float sigMT_t6wg_mg(0);
  float sigMET_t6wg_mg(0);
  float HT_t6wg_mg(0);
	float sigMETJESup_t6wg_mg(0);
	float sigMETJESdo_t6wg_mg(0);
	float sigMETJERup_t6wg_mg(0);
	float sigMETJERdo_t6wg_mg(0);
	float sigMTJESup_t6wg_mg(0);
	float sigMTJESdo_t6wg_mg(0);
	float sigMTJERup_t6wg_mg(0);
	float sigMTJERdo_t6wg_mg(0);
	float HTJESup_t6wg_mg(0);
	float HTJESdo_t6wg_mg(0);
	float squarkMass_t6wg_mg(0);
  float charginoMass_t6wg_mg(0);
  float neutralinoMass_t6wg_mg(0);
  mgtree_t6wg->SetBranchAddress("phoEt",      &phoEt_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("phoEta",     &phoEta_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("lepPt",      &lepPt_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("lepEta",     &lepEta_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("sigMT",      &sigMT_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("sigMET",     &sigMET_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("HT",         &HT_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("Msquark",    &squarkMass_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("Mchargino",  &charginoMass_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("Mneutralino",&neutralinoMass_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("HTJESup",    &HTJESup_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("HTJESdo",    &HTJESdo_t6wg_mg);

	for(int ievt(0); ievt < mgtree_t6wg->GetEntries(); ievt++){
		mgtree_t6wg->GetEntry(ievt);

    double NLSPMass(0);
    if(charginoMass_t6wg_mg >0)NLSPMass=charginoMass_t6wg_mg;
		else if(neutralinoMass_t6wg_mg >0)NLSPMass=neutralinoMass_t6wg_mg;
		if(NLSPMass <= 0)continue;
		// ** cut flow *****
		if(phoEt_t6wg_mg < 35 || lepPt_t6wg_mg < 20)continue;
		if(fabs(phoEta_t6wg_mg) > 1.4442 || fabs(lepEta_t6wg_mg) > 2.5)continue;

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_t6wg_mg,lepEta_t6wg_mg)*objectESF.getPhotonESF(phoEt_t6wg_mg, phoEta_t6wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t6wg_mg, lepPt_t6wg_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_t6wg_mg,lepEta_t6wg_mg)*objectESF.getPhotonESF(phoEt_t6wg_mg, phoEta_t6wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t6wg_mg, lepPt_t6wg_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_t6wg_mg, phoEta_t6wg_mg)*objectESF.getFastMuonESF(lepPt_t6wg_mg,lepEta_t6wg_mg)*objectESF.getFastMuonEGTRGESF(phoEt_t6wg_mg, lepPt_t6wg_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_t6wg_mg, lepPt_t6wg_mg)*objectESF.getFastMuonESF(lepPt_t6wg_mg,lepEta_t6wg_mg)*objectESF.getPhotonESF(phoEt_t6wg_mg, phoEta_t6wg_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		if(sigMET_t6wg_mg < 120 || sigMT_t6wg_mg < 100)continue;
    p_T6WGselect->Fill(squarkMass_t6wg_mg, NLSPMass);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg);
		if(SigBinIndex >=0){
			t6wg_h_chan_rate_nom[SigBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor);
			t6wg_h_chan_rate_xsUp[SigBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor); 
			t6wg_h_chan_rate_esfUp[SigBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactorup); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_t6wg_mg, HTJESup_t6wg_mg, phoEt_t6wg_mg); 
		if(jesupBinIndex >=0){
			t6wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_t6wg_mg, HTJESdo_t6wg_mg, phoEt_t6wg_mg);
		if(jesdoBinIndex >=0){
			t6wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg); 
		if( jerupBinIndex >=0){
			t6wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg);
		if(jerdoBinIndex >= 0){
			t6wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( squarkMass_t6wg_mg, NLSPMass, scalefactor);
		}  
	}

  TChain *egtree_t6wg;
  egtree_t6wg = new TChain("egTree","egTree");
  egtree_t6wg->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T6WG.root");
  float phoEt_t6wg_eg(0);
  float phoEta_t6wg_eg(0);
  float lepPt_t6wg_eg(0);
  float lepEta_t6wg_eg(0);
  float sigMT_t6wg_eg(0);
  float sigMET_t6wg_eg(0);
  float HT_t6wg_eg(0);
	float sigMETJESup_t6wg_eg(0);
	float sigMETJESdo_t6wg_eg(0);
	float sigMETJERup_t6wg_eg(0);
	float sigMETJERdo_t6wg_eg(0);
	float sigMTJESup_t6wg_eg(0);
	float sigMTJESdo_t6wg_eg(0);
	float sigMTJERup_t6wg_eg(0);
	float sigMTJERdo_t6wg_eg(0);
	float HTJESup_t6wg_eg(0);
	float HTJESdo_t6wg_eg(0);
	float squarkMass_t6wg_eg(0);
  float charginoMass_t6wg_eg(0);
  float neutralinoMass_t6wg_eg(0);
  egtree_t6wg->SetBranchAddress("phoEt",      &phoEt_t6wg_eg);
  egtree_t6wg->SetBranchAddress("phoEta",     &phoEta_t6wg_eg);
  egtree_t6wg->SetBranchAddress("lepPt",      &lepPt_t6wg_eg);
  egtree_t6wg->SetBranchAddress("lepEta",     &lepEta_t6wg_eg);
  egtree_t6wg->SetBranchAddress("sigMT",      &sigMT_t6wg_eg);
  egtree_t6wg->SetBranchAddress("sigMET",     &sigMET_t6wg_eg);
  egtree_t6wg->SetBranchAddress("HT",         &HT_t6wg_eg);
	egtree_t6wg->SetBranchAddress("Msquark",    &squarkMass_t6wg_eg);
  egtree_t6wg->SetBranchAddress("Mchargino",  &charginoMass_t6wg_eg);
  egtree_t6wg->SetBranchAddress("Mneutralino",&neutralinoMass_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t6wg_eg);
	egtree_t6wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t6wg_eg);
	egtree_t6wg->SetBranchAddress("HTJESup",    &HTJESup_t6wg_eg);
	egtree_t6wg->SetBranchAddress("HTJESdo",    &HTJESdo_t6wg_eg);

	for(int ievt(0); ievt < egtree_t6wg->GetEntries(); ievt++){
		egtree_t6wg->GetEntry(ievt);

    double NLSPMass(0);
    if(charginoMass_t6wg_eg >0)NLSPMass=charginoMass_t6wg_eg;
		else if(neutralinoMass_t6wg_eg >0)NLSPMass=neutralinoMass_t6wg_eg;
		if(NLSPMass <= 0)continue;
		// ** cut flow *****
		if(phoEt_t6wg_eg < 35 || lepPt_t6wg_eg < 25)continue;
		if(fabs(phoEta_t6wg_eg) > 1.4442 || fabs(lepEta_t6wg_eg) > 2.5)continue;
		if(sigMET_t6wg_eg < 120 || sigMT_t6wg_eg < 100)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg) + NBIN;
		if(SigBinIndex >=0){
			t6wg_h_chan_rate_nom[SigBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor);
			t6wg_h_chan_rate_xsUp[SigBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor); 
			t6wg_h_chan_rate_esfUp[SigBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactorup); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_t6wg_eg, HTJESup_t6wg_eg, phoEt_t6wg_eg) + NBIN;
		if(jesupBinIndex >=0){
			t6wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_t6wg_eg, HTJESdo_t6wg_eg, phoEt_t6wg_eg) + NBIN;
		if(jesdoBinIndex >=0){
			t6wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg) + NBIN;
		if( jerupBinIndex >=0){
			t6wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg) + NBIN;	
		if(jerdoBinIndex >= 0){
			t6wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( squarkMass_t6wg_eg, NLSPMass, scalefactor);
		}  
	}


	for(int ih(0); ih < NBIN*2; ih++){
		if(ih < NBIN) 	lumi = lumi_mg;
		else 		lumi = lumi_eg;	
		for(int i(1); i < t6wg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(int j(1); j < t6wg_h_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){

				if(p_T6WGMASS->GetBinContent(i,j) < 1000){
					p_T6WGselect->SetBinContent(i,j,-1);
					p_T6WGMASS->SetBinContent(i,j,-1);
				}
		
				if(p_T6WGMASS->GetBinContent(i,j) <= 0){
					t6wg_h_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					t6wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					t6wg_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
				}
				else{
					float noe = p_T6WGMASS->GetBinContent(i,j);
					float sparticleMass = p_T6WGMASS->GetXaxis()->GetBinCenter(i);
					float crosssection = p_crosssection_t6wg->GetBinContent( p_crosssection_t6wg->FindBin(sparticleMass) );
					float crosssectionUp = (crosssection+ p_crosssection_t6wg->GetBinError( p_crosssection_t6wg->FindBin(sparticleMass) ) ); 

				  	t6wg_h_chan_rate_nom[ih]->SetBinError(i,j, sqrt(t6wg_h_chan_rate_nom[ih]->GetBinContent(i))*lumi*crosssection/noe);
					t6wg_h_chan_rate_nom[ih]->SetBinContent(i,j, t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)*lumi*crosssection/noe); 
					t6wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j, t6wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe); 
					t6wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, t6wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t6wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, t6wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t6wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j, t6wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
					t6wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, t6wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)*lumi*crosssectionUp/noe);
					t6wg_h_chan_rate_esfUp[ih]->SetBinContent(i,j, t6wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)*lumi*crosssection/noe);
									
					t6wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(t6wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t6wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t6wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(t6wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t6wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t6wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( t6wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)));
					t6wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					t6wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, fabs( t6wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)-t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					t6wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j));   
					t6wg_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 

				}
			}
		} 
	} 

	outputfile_t6wg->Write();
	outputfile_t6wg->Close(); */


  xSecFile.Close();
}
