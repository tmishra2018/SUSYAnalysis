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
#define NPDF 9

void analysis_TChiWG_noSRcuts(){//main  
	gROOT->SetBatch(true);
	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;
	std::string whichVFP;
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";
	
	float flat_pu_unc;
	if     (RunYear == 2016) flat_pu_unc = 0.022;  // 2.2
	else if(RunYear == 2017) flat_pu_unc = 0.023;  // 2.3
	else                     flat_pu_unc = 0.025;  // 2.5 (2018)


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

  	gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/lib/libAnaClasses.so");

	TFile xSecFile("/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/test/cross/ULsusyCrossSection.root");
	TH1D *p_crosssection_tchiwg = (TH1D*)xSecFile.Get("p_charginoSec");
	TH1D *p_crosssection_t5wg   = (TH1D*)xSecFile.Get("p_gluinoxSec");
	TH1D *p_crosssection_t6wg   = (TH1D*)xSecFile.Get("p_squarkxSec");

	TChain *datachain = new TChain("signalTree");
	datachain->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  	datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	std::ostringstream histname;
	//**************   T5WG  *************************** //
	TFile *file_t5wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_%d%s.root",RunYear, whichVFP.c_str()));
  	TTree *tree_t5wg = (TTree*)file_t5wg->Get("SUSYtree");
	float Mgluino_t5wg(0);
  	float Mchargino_t5wg(0);
  	float Mneutralino_t5wg(0);
	int   nVertex(0);
	std::vector<float> *ScaleSystWeight_t5wg = 0;

	tree_t5wg->SetBranchAddress("MsGsQ",    &Mgluino_t5wg);  
  	tree_t5wg->SetBranchAddress("Mchargino",  &Mchargino_t5wg);
  	tree_t5wg->SetBranchAddress("Mneutralino",&Mneutralino_t5wg);
  	tree_t5wg->SetBranchAddress("nVertex",    &nVertex);
	tree_t5wg->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight_t5wg);

	TFile *outputfile_t5wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/Tree_T5WG_%d%s.root",RunYear, whichVFP.c_str()),"RECREATE");
	outputfile_t5wg->cd();

	double xBins[] = {950, 1050, 1150, 1250, 1350, 1450, 1550, 1650, 1750, 1850, 1950, 2025, 2075, 2125 	 };
	//		, 2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525, 2575, 2625, 2675, 2725, 2775, 2825}; // cropping temperarily fot better limit plot

	int nXbins = sizeof(xBins) / sizeof(xBins[0]) - 1;

	double yBins[] = {25, 75, 125, 175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 
		1325, 1375, 1425, 1475, 1525, 1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975, 2025, 2075, 2125 };
	//	, 2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525, 2575, 2625, 2675, 2725, 2775, 2825};
	int nYbins = sizeof(yBins) / sizeof(yBins[0]) - 1;

	TH2D *p_T5WGMASS         = new TH2D("SUSYMass","",nXbins, xBins, nYbins, yBins);
	TH2D *p_T5WGselect       = new TH2D("p_T5WGselect","",nXbins, xBins, nYbins, yBins); 
	TH2D *p_T5WGMASS_pdf[9];

	for(unsigned i(0); i < 9; i++){
                histname.str("");
                histname << "p_T5WGMASS_pdf" << i;
                p_T5WGMASS_pdf[i] = new TH2D(histname.str().c_str(), histname.str().c_str(), nXbins,xBins,nYbins,yBins);
        }
	
	// NEW: PU-weighted and PU-varied 2D histograms for syserr_PU
	TH2D *t5wg_h_puNom  = new TH2D("t5wg_h_puNom", "",nXbins,xBins,nYbins,yBins);
	TH2D *t5wg_h_puUp   = new TH2D("t5wg_h_puUp",  "",nXbins,xBins,nYbins,yBins);
	TH2D *t5wg_h_puDown = new TH2D("t5wg_h_puDown","",nXbins,xBins,nYbins,yBins);

	for(int ievt(0); ievt < tree_t5wg->GetEntries(); ievt++){
		tree_t5wg->GetEntry(ievt);
    		if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		double NLSPMass(0);
    		if(Mchargino_t5wg >0)NLSPMass = Mchargino_t5wg;
		else if(Mneutralino_t5wg >0)NLSPMass = Mneutralino_t5wg;
		if(NLSPMass <= 0) continue;
		p_T5WGMASS->Fill(Mgluino_t5wg, NLSPMass);
		for(unsigned i(0); i < NPDF; i++)p_T5WGMASS_pdf[i]->Fill(Mgluino_t5wg, NLSPMass, (*ScaleSystWeight_t5wg)[i]);
	}
	
	TH2D *t5wg_h_chan_rate_nom[NBIN*2]; 
	for(int i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		t5wg_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
	}

	TH2D *t5wg_h_chan_rate_jesUp[NBIN*2]; 
	TH2D *t5wg_h_chan_rate_jesDown[NBIN*2];    
	TH2D *t5wg_h_chan_rate_jerUp[NBIN*2];    
	TH2D *t5wg_h_chan_rate_jerDown[NBIN*2];    
	TH2D *t5wg_h_chan_rate_xsUp[NBIN*2];       
	TH2D *t5wg_h_chan_rate_esfUp[NBIN*2];       
        TH2D *t5wg_h_chan_rate_puUp[NBIN*2];
        TH2D *t5wg_h_chan_rate_puDown[NBIN*2];
	TH2D *t5wg_h_chan_rate_qcdscale[NBIN*2][9];

	TH2D *t5wg_h_syserr_PU = new TH2D("t5wg_h_syserr_PU","t5wg_h_syserr_PU", nXbins, xBins, nYbins, yBins);
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
	TH2D *t5wg_h_chan_syserr_pu[NBIN*2];
	TH2D *t5wg_h_chan_syserr_qcdscale[NBIN*2];

	
	for(int i(0); i < NBIN*2; i++){
		for(unsigned is(0); is<9; is++){
			histname.str("");
			histname << "t5wg_h_chan" << i+1 << "_rate_qcdscale" << is;
	                t5wg_h_chan_rate_qcdscale[i][is] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		}

		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		t5wg_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		t5wg_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		t5wg_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		t5wg_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		t5wg_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		t5wg_h_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
                histname << "h_chan" << i+1 << "_rate_puUp";
                t5wg_h_chan_rate_puUp[i]    = new TH2D(histname.str().c_str(),histname.str().c_str(), nXbins, xBins, nYbins, yBins);
                histname.str("");
                histname << "h_chan" << i+1 << "_rate_puDown";
                t5wg_h_chan_rate_puDown[i]  = new TH2D(histname.str().c_str(),histname.str().c_str(), nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		t5wg_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		t5wg_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		t5wg_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		t5wg_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		t5wg_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		t5wg_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		t5wg_h_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		t5wg_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		t5wg_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		t5wg_h_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
                histname << "h_chan" << i+1 << "_syserr_pu";
                t5wg_h_chan_syserr_pu[i]    = new TH2D(histname.str().c_str(),histname.str().c_str(), nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdscale";
		t5wg_h_chan_syserr_qcdscale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(), nXbins, xBins, nYbins, yBins);
	}
	
	TH1D *p_t5wg_MET_signal_mg = new TH1D("p_t5wg_MET_signal_1700_1000_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_mg  = new TH1D("p_t5wg_HT_signal_1700_1000_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_mg = new TH1D("p_t5wg_PhoEt_signal_1700_1000_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_mg = new TH1D("p_t5wg_Mt_signal_1700_1000_mg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_MET_signal_eg = new TH1D("p_t5wg_MET_signal_1700_1000_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_eg  = new TH1D("p_t5wg_HT_signal_1700_1000_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_eg = new TH1D("p_t5wg_PhoEt_signal_1700_1000_eg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_eg = new TH1D("p_t5wg_Mt_signal_1700_1000_eg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_1700_1000_nom = new TH1D("p_t5wg_1700_1000_nom", "p_t5wg_1700_1000_nom", NBIN*2,0,NBIN*2);

	// ---- T5WG(1987,1900), added alongside the existing 1700,1000 point ----
	TH1D *p_t5wg_MET_signal_1987_1900_mg = new TH1D("p_t5wg_MET_signal_1987_1900_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1987_1900_mg  = new TH1D("p_t5wg_HT_signal_1987_1900_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1987_1900_mg = new TH1D("p_t5wg_PhoEt_signal_1987_1900_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1987_1900_mg = new TH1D("p_t5wg_Mt_signal_1987_1900_mg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_MET_signal_1987_1900_eg = new TH1D("p_t5wg_MET_signal_1987_1900_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1987_1900_eg  = new TH1D("p_t5wg_HT_signal_1987_1900_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1987_1900_eg = new TH1D("p_t5wg_PhoEt_signal_1987_1900_eg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1987_1900_eg = new TH1D("p_t5wg_Mt_signal_1987_1900_eg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_1987_1900_nom = new TH1D("p_t5wg_1987_1900_nom", "p_t5wg_1987_1900_nom", NBIN*2,0,NBIN*2);

	// ---- T5WG(1800,600), new mass point ----
	TH1D *p_t5wg_MET_signal_1800_600_mg = new TH1D("p_t5wg_MET_signal_1800_600_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1800_600_mg  = new TH1D("p_t5wg_HT_signal_1800_600_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1800_600_mg = new TH1D("p_t5wg_PhoEt_signal_1800_600_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1800_600_mg = new TH1D("p_t5wg_Mt_signal_1800_600_mg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_MET_signal_1800_600_eg = new TH1D("p_t5wg_MET_signal_1800_600_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1800_600_eg  = new TH1D("p_t5wg_HT_signal_1800_600_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1800_600_eg = new TH1D("p_t5wg_PhoEt_signal_1800_600_eg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1800_600_eg = new TH1D("p_t5wg_Mt_signal_1800_600_eg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_1800_600_nom = new TH1D("p_t5wg_1800_600_nom", "p_t5wg_1800_600_nom", NBIN*2,0,NBIN*2);

	// ---- T5WG(1900,100), new mass point ----
	TH1D *p_t5wg_MET_signal_1900_100_mg = new TH1D("p_t5wg_MET_signal_1900_100_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1900_100_mg  = new TH1D("p_t5wg_HT_signal_1900_100_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1900_100_mg = new TH1D("p_t5wg_PhoEt_signal_1900_100_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1900_100_mg = new TH1D("p_t5wg_Mt_signal_1900_100_mg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_MET_signal_1900_100_eg = new TH1D("p_t5wg_MET_signal_1900_100_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t5wg_HT_signal_1900_100_eg  = new TH1D("p_t5wg_HT_signal_1900_100_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t5wg_PhoEt_signal_1900_100_eg = new TH1D("p_t5wg_PhoEt_signal_1900_100_eg","",nSigEtBins, sigEtBins);
	TH1D *p_t5wg_Mt_signal_1900_100_eg = new TH1D("p_t5wg_Mt_signal_1900_100_eg","",nSigMtBins, sigMtBins);
	TH1D *p_t5wg_1900_100_nom = new TH1D("p_t5wg_1900_100_nom", "p_t5wg_1900_100_nom", NBIN*2,0,NBIN*2);
	
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
	std::vector<float> *ScaleSystWeight_t5wg_mg = 0;

	float PUweight_t5wg_mg(1), PUweightUp_t5wg_mg(1), PUweightDown_t5wg_mg(1);
  	mgtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_mg);
  	mgtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_mg);
  	mgtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("nJet",         &nJet_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("nVertex",    &nVertex_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("MsGsQ",    &gluinoMass_t5wg_mg);
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
		mgtree_t5wg->SetBranchAddress("PUweight",     &PUweight_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("PUweightUp",   &PUweightUp_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("PUweightDown", &PUweightDown_t5wg_mg);
	        mgtree_t5wg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_t5wg_mg);


	for(int ievt(0); ievt < mgtree_t5wg->GetEntries(); ievt++){
		mgtree_t5wg->GetEntry(ievt);

		double NLSPMass(0);
		if(charginoMass_t5wg_mg > 0)NLSPMass = charginoMass_t5wg_mg;
		else if(neutralinoMass_t5wg_mg > 0)NLSPMass = neutralinoMass_t5wg_mg;
		if(NLSPMass <= 0)continue;	
	
		// cut flow 
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
			p_t5wg_MET_signal_mg->Fill(sigMET_t5wg_mg);
			p_t5wg_Mt_signal_mg->Fill(sigMT_t5wg_mg);
			p_t5wg_HT_signal_mg->Fill(HT_t5wg_mg);
			p_t5wg_PhoEt_signal_mg->Fill(phoEt_t5wg_mg);
		}
		if( fabs(gluinoMass_t5wg_mg - 1987) < 10 && fabs(NLSPMass - 1900) < 5){
			p_t5wg_MET_signal_1987_1900_mg->Fill(sigMET_t5wg_mg);
			p_t5wg_Mt_signal_1987_1900_mg->Fill(sigMT_t5wg_mg);
			p_t5wg_HT_signal_1987_1900_mg->Fill(HT_t5wg_mg);
			p_t5wg_PhoEt_signal_1987_1900_mg->Fill(phoEt_t5wg_mg);
		}
		// ---- T5WG(1800,600), new mass point ----
		if( fabs(gluinoMass_t5wg_mg - 1800) < 10 && fabs(NLSPMass - 600) < 5){
			p_t5wg_MET_signal_1800_600_mg->Fill(sigMET_t5wg_mg);
			p_t5wg_Mt_signal_1800_600_mg->Fill(sigMT_t5wg_mg);
			p_t5wg_HT_signal_1800_600_mg->Fill(HT_t5wg_mg);
			p_t5wg_PhoEt_signal_1800_600_mg->Fill(phoEt_t5wg_mg);
		}
		// ---- T5WG(1900,100), new mass point ----
		if( fabs(gluinoMass_t5wg_mg - 1900) < 10 && fabs(NLSPMass - 100) < 5){
			p_t5wg_MET_signal_1900_100_mg->Fill(sigMET_t5wg_mg);
			p_t5wg_Mt_signal_1900_100_mg->Fill(sigMT_t5wg_mg);
			p_t5wg_HT_signal_1900_100_mg->Fill(HT_t5wg_mg);
			p_t5wg_PhoEt_signal_1900_100_mg->Fill(phoEt_t5wg_mg);
		}

		if(sigMET_t5wg_mg < 120 || sigMT_t5wg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg);
		if(SigBinIndex >= 0){
			// Nominal: scalefactor × PUweight (nominal)
			t5wg_h_chan_rate_nom   [SigBinIndex]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor * PUweight_t5wg_mg);
			// PU varied: only PU weight changes
			t5wg_h_chan_rate_puUp  [SigBinIndex]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor * PUweightUp_t5wg_mg);   // NEW
			t5wg_h_chan_rate_puDown[SigBinIndex]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor * PUweightDown_t5wg_mg); // NEW
			t5wg_h_chan_rate_xsUp  [SigBinIndex]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor * PUweight_t5wg_mg);
			t5wg_h_chan_rate_esfUp [SigBinIndex]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactorup * PUweight_t5wg_mg);

			for(unsigned is(0); is < 9; is++)
				t5wg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor*(*ScaleSystWeight_t5wg_mg)[is]);

			if(fabs(gluinoMass_t5wg_mg-1700)<10 && fabs(NLSPMass-1000)<5)
				p_t5wg_1700_1000_nom->Fill(SigBinIndex, scalefactor * PUweight_t5wg_mg);
			if(fabs(gluinoMass_t5wg_mg-1987)<10 && fabs(NLSPMass-1900)<5)
				p_t5wg_1987_1900_nom->Fill(SigBinIndex, scalefactor * PUweight_t5wg_mg);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_t5wg_mg, HTJESup_t5wg_mg, phoEt_t5wg_mg);
		if(rawJesup  >= 0) t5wg_h_chan_rate_jesUp  [rawJesup ]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor*PUweight_t5wg_mg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_t5wg_mg, HTJESdo_t5wg_mg, phoEt_t5wg_mg);
		if(rawJesdo  >= 0) t5wg_h_chan_rate_jesDown[rawJesdo ]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor*PUweight_t5wg_mg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg);
		if(rawJerup  >= 0) t5wg_h_chan_rate_jerUp  [rawJerup ]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor*PUweight_t5wg_mg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_t5wg_mg, HT_t5wg_mg, phoEt_t5wg_mg);
		if(rawJerdo  >= 0) t5wg_h_chan_rate_jerDown[rawJerdo ]->Fill(gluinoMass_t5wg_mg, NLSPMass, scalefactor*PUweight_t5wg_mg);
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
	float PUweight_t5wg_eg(1), PUweightUp_t5wg_eg(1), PUweightDown_t5wg_eg(1);
	        std::vector<float> *ScaleSystWeight_t5wg_eg = 0;

  egtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_eg);
	egtree_t5wg->SetBranchAddress("nVertex",    &nVertex_t5wg_eg);
  egtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_eg);
  egtree_t5wg->SetBranchAddress("nJet",         &nJet_t5wg_eg);
	egtree_t5wg->SetBranchAddress("MsGsQ",    &gluinoMass_t5wg_eg);
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
		egtree_t5wg->SetBranchAddress("PUweight",     &PUweight_t5wg_eg);     // NEW
	egtree_t5wg->SetBranchAddress("PUweightUp",   &PUweightUp_t5wg_eg);   // NEW
	egtree_t5wg->SetBranchAddress("PUweightDown", &PUweightDown_t5wg_eg); // NEW
	egtree_t5wg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_t5wg_eg);


	for(int ievt(0); ievt < egtree_t5wg->GetEntries(); ievt++){
		egtree_t5wg->GetEntry(ievt);

		double NLSPMass(0);
		if(charginoMass_t5wg_eg > 0)NLSPMass = charginoMass_t5wg_eg;
		else if(neutralinoMass_t5wg_eg > 0)NLSPMass = neutralinoMass_t5wg_eg;
		if(NLSPMass <= 0)continue;	
		// cut flow 
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
			p_t5wg_MET_signal_eg->Fill(sigMET_t5wg_eg);
			p_t5wg_Mt_signal_eg->Fill(sigMT_t5wg_eg);
			p_t5wg_HT_signal_eg->Fill(HT_t5wg_eg);
			p_t5wg_PhoEt_signal_eg->Fill(phoEt_t5wg_eg);
		}
		if( fabs(gluinoMass_t5wg_eg - 1987) < 10 && fabs(NLSPMass - 1900) < 5){
			p_t5wg_MET_signal_1987_1900_eg->Fill(sigMET_t5wg_eg);
			p_t5wg_Mt_signal_1987_1900_eg->Fill(sigMT_t5wg_eg);
			p_t5wg_HT_signal_1987_1900_eg->Fill(HT_t5wg_eg);
			p_t5wg_PhoEt_signal_1987_1900_eg->Fill(phoEt_t5wg_eg);
		}
		// ---- T5WG(1800,600), new mass point ----
		if( fabs(gluinoMass_t5wg_eg - 1800) < 10 && fabs(NLSPMass - 600) < 5){
			p_t5wg_MET_signal_1800_600_eg->Fill(sigMET_t5wg_eg);
			p_t5wg_Mt_signal_1800_600_eg->Fill(sigMT_t5wg_eg);
			p_t5wg_HT_signal_1800_600_eg->Fill(HT_t5wg_eg);
			p_t5wg_PhoEt_signal_1800_600_eg->Fill(phoEt_t5wg_eg);
		}
		// ---- T5WG(1900,100), new mass point ----
		if( fabs(gluinoMass_t5wg_eg - 1900) < 10 && fabs(NLSPMass - 100) < 5){
			p_t5wg_MET_signal_1900_100_eg->Fill(sigMET_t5wg_eg);
			p_t5wg_Mt_signal_1900_100_eg->Fill(sigMT_t5wg_eg);
			p_t5wg_HT_signal_1900_100_eg->Fill(HT_t5wg_eg);
			p_t5wg_PhoEt_signal_1900_100_eg->Fill(phoEt_t5wg_eg);
		}

		if(sigMET_t5wg_eg < 120 || sigMT_t5wg_eg < 100)continue;
				// FIX: guard on raw index BEFORE adding +NBIN offset
		int rawSigBin = Bin.findSignalBin(sigMET_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg);
		if(rawSigBin >= 0){
			int SigBinIndex = rawSigBin + NBIN;
			t5wg_h_chan_rate_nom   [SigBinIndex]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor * PUweight_t5wg_eg);
			t5wg_h_chan_rate_puUp  [SigBinIndex]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor * PUweightUp_t5wg_eg);   // NEW
			t5wg_h_chan_rate_puDown[SigBinIndex]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor * PUweightDown_t5wg_eg); // NEW
			t5wg_h_chan_rate_xsUp  [SigBinIndex]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor * PUweight_t5wg_eg);
			t5wg_h_chan_rate_esfUp [SigBinIndex]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactorup * PUweight_t5wg_eg);
			for(unsigned is(0); is < 9; is++)
				t5wg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor*(*ScaleSystWeight_t5wg_eg)[is]);
			if(fabs(gluinoMass_t5wg_eg-1700)<10 && fabs(NLSPMass-1000)<5)
				p_t5wg_1700_1000_nom->Fill(SigBinIndex, scalefactor * PUweight_t5wg_eg);
			if(fabs(gluinoMass_t5wg_eg-1987)<10 && fabs(NLSPMass-1900)<5)
				p_t5wg_1987_1900_nom->Fill(SigBinIndex, scalefactor * PUweight_t5wg_eg);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_t5wg_eg, HTJESup_t5wg_eg, phoEt_t5wg_eg);
		if(rawJesup  >= 0) t5wg_h_chan_rate_jesUp  [rawJesup+NBIN]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor*PUweight_t5wg_eg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_t5wg_eg, HTJESdo_t5wg_eg, phoEt_t5wg_eg);
		if(rawJesdo  >= 0) t5wg_h_chan_rate_jesDown[rawJesdo+NBIN]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor*PUweight_t5wg_eg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg);
		if(rawJerup  >= 0) t5wg_h_chan_rate_jerUp  [rawJerup+NBIN]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor*PUweight_t5wg_eg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_t5wg_eg, HT_t5wg_eg, phoEt_t5wg_eg);
		if(rawJerdo  >= 0) t5wg_h_chan_rate_jerDown[rawJerdo+NBIN]->Fill(gluinoMass_t5wg_eg, NLSPMass, scalefactor*PUweight_t5wg_eg);
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
					t5wg_h_chan_rate_puUp[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_rate_puDown[ih]->SetBinContent(i,j, -1);

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
					t5wg_h_chan_syserr_qcdscale[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_pu  [ih]->SetBinContent(i,j,-1); // NEW
				}

				else{
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(i);
					float crosssection = p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					float crosssectionUp = (crosssection+ p_crosssection_t5wg->GetBinError( p_crosssection_t5wg->FindBin(sparticleMass) ) );
					float sc   = lumi*crosssection/noe;
					
					t5wg_h_chan_rate_nom[ih]->SetBinError  (i,j, sqrt(t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) * sc);
					t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j) * sc);
					if(ih%3 == 0){
                                                double acc_nom = 0;
                                                acc_nom += (sc > 0 ? t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j) / sc : 0) / noe; 
                                                for(unsigned binHT(1); binHT < 3; binHT++)
                                                        acc_nom += t5wg_h_chan_rate_nom[ih+binHT]->GetBinContent(i,j) / noe;        

                                                const int scaleIdx[] = {1, 2, 3, 4, 6, 8};
                                                double pdfdiff(0);
                                                for(int is : scaleIdx){
                                                        double noepdf = p_T5WGMASS_pdf[is]->GetBinContent(i,j);
                                                        if(noepdf <= 0) continue;
                                                        double acc_var = 0;
                                                        for(unsigned binHT(0); binHT < 3; binHT++)
                                                                acc_var += t5wg_h_chan_rate_qcdscale[ih+binHT][is]->GetBinContent(i,j) / noepdf;
                                                        pdfdiff = std::max(pdfdiff, fabs(acc_var - acc_nom));
                                                }

                                                double N_raw_ih   = (sc > 0 ? t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j) / sc : 0);
                                                double N_raw_tot  = N_raw_ih;
                                                for(unsigned binHT(1); binHT < 3; binHT++)
                                                        N_raw_tot += t5wg_h_chan_rate_nom[ih+binHT]->GetBinContent(i,j);

                                                for(unsigned binHT(0); binHT < 3; binHT++){
                                                        if(t5wg_h_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0){
                                                                if(acc_nom > 0 && N_raw_tot >= 3.0)
                                                                        t5wg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i,j, pdfdiff / acc_nom);
                                                                else
                                                                        t5wg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i,j, 0.30);
                                                        }
                                                }
                                        }

					t5wg_h_chan_rate_puUp  [ih]->SetBinContent(i,j, t5wg_h_chan_rate_puUp  [ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_puDown[ih]->SetBinContent(i,j, t5wg_h_chan_rate_puDown[ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_jesUp [ih]->SetBinContent(i,j, t5wg_h_chan_rate_jesUp [ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j,t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_jerUp [ih]->SetBinContent(i,j, t5wg_h_chan_rate_jerUp [ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*sc);
					t5wg_h_chan_rate_xsUp  [ih]->SetBinContent(i,j, t5wg_h_chan_rate_xsUp  [ih]->GetBinContent(i,j)*lumi*crosssectionUp/noe);
					t5wg_h_chan_rate_esfUp [ih]->SetBinContent(i,j, t5wg_h_chan_rate_esfUp [ih]->GetBinContent(i,j)*sc);
					
					float nom     = t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j);
					float jesUp   = t5wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j);
					float jesDown = t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j);
					float jerUp   = t5wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j);
					float jerDown = t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j);

					float jesErr = std::max(fabs(jesUp-nom), fabs(jesDown-nom));
					float jerErr = std::max(fabs(jerUp-nom), fabs(jerDown-nom));

					const float LOW_STAT_YIELD = 2.0;

					if (nom > 0.0 && nom < LOW_STAT_YIELD && (jesErr  > nom || jerErr > nom)) {

    						float maxShift = std::max({nom, jesUp, jesDown, jerUp, jerDown});
    						float newNom   = 0.5 * maxShift;

    						t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, newNom);
    						float statErr = t5wg_h_chan_rate_nom[ih]->GetBinError(i,j);
    						t5wg_h_chan_rate_nom[ih]->SetBinError  (i,j, statErr * newNom / nom);

    						std::cout << "[JES/JER FIX] Year=" << RunYear << " bin=(" << i << "," << j << ")" 
							<< " nom=" << nom << " newNom=" << newNom << std::endl;
    						nom = newNom;
					}
					t5wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, std::max(fabs(jesUp-nom), fabs(jesDown-nom)));
					t5wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, std::max(fabs(jerUp-nom), fabs(jerDown-nom)));
					t5wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( t5wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					t5wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, fabs( t5wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					t5wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j));   
					t5wg_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					
					float puUp_sc   = t5wg_h_chan_rate_puUp  [ih]->GetBinContent(i,j);
					float puDown_sc = t5wg_h_chan_rate_puDown[ih]->GetBinContent(i,j);
					float puErr     = std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom));
					float flat_pu_err = flat_pu_unc * nom;
					float final_pu_err;
					if(nom >= 1.0)	final_pu_err = std::min(puErr, 3.0f * flat_pu_err);
					else final_pu_err = flat_pu_err;
					t5wg_h_chan_syserr_pu[ih]->SetBinContent(i,j, final_pu_err);
					if(ih==0) t5wg_h_syserr_PU->SetBinContent(i,j, std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom)));

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
					p_t5wg_Mt_signal_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_MET_signal_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_HT_signal_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_PhoEt_signal_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_Mt_signal_eg->Scale(lumi_eg*crosssection/noe);
					
					for(int s(1); s <= NBIN*2; s++){
                				if(s <= NBIN)   lumi = lumi_mg;
                				else            lumi = lumi_eg;
						p_t5wg_1700_1000_nom->SetBinError(s, sqrt(p_t5wg_1700_1000_nom->GetBinContent(s))*lumi*crosssection/noe);
						p_t5wg_1700_1000_nom->SetBinContent(s, p_t5wg_1700_1000_nom->GetBinContent(s)*lumi*crosssection/noe);}
				}
				// ---- T5WG(1987,1900) normalization, mirroring the 1700,1000 block above ----
				if( fabs(sparticleMass - 1987 ) < 10 && fabs( p_T5WGMASS->GetYaxis()->GetBinCenter(j) - 1900) < 5){
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					p_t5wg_MET_signal_1987_1900_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_HT_signal_1987_1900_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_PhoEt_signal_1987_1900_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_Mt_signal_1987_1900_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_MET_signal_1987_1900_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_HT_signal_1987_1900_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_PhoEt_signal_1987_1900_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_Mt_signal_1987_1900_eg->Scale(lumi_eg*crosssection/noe);

					for(int s(1); s <= NBIN*2; s++){
                				if(s <= NBIN)   lumi = lumi_mg;
                				else            lumi = lumi_eg;
						p_t5wg_1987_1900_nom->SetBinError(s, sqrt(p_t5wg_1987_1900_nom->GetBinContent(s))*lumi*crosssection/noe);
						p_t5wg_1987_1900_nom->SetBinContent(s, p_t5wg_1987_1900_nom->GetBinContent(s)*lumi*crosssection/noe);}
				}
				// ---- T5WG(1800,600) normalization, new mass point ----
				if( fabs(sparticleMass - 1800 ) < 10 && fabs( p_T5WGMASS->GetYaxis()->GetBinCenter(j) - 600) < 5){
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					p_t5wg_MET_signal_1800_600_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_HT_signal_1800_600_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_PhoEt_signal_1800_600_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_Mt_signal_1800_600_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_MET_signal_1800_600_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_HT_signal_1800_600_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_PhoEt_signal_1800_600_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_Mt_signal_1800_600_eg->Scale(lumi_eg*crosssection/noe);

					for(int s(1); s <= NBIN*2; s++){
                				if(s <= NBIN)   lumi = lumi_mg;
                				else            lumi = lumi_eg;
						p_t5wg_1800_600_nom->SetBinError(s, sqrt(p_t5wg_1800_600_nom->GetBinContent(s))*lumi*crosssection/noe);
						p_t5wg_1800_600_nom->SetBinContent(s, p_t5wg_1800_600_nom->GetBinContent(s)*lumi*crosssection/noe);}
				}
				// ---- T5WG(1900,100) normalization, new mass point ----
				if( fabs(sparticleMass - 1900 ) < 10 && fabs( p_T5WGMASS->GetYaxis()->GetBinCenter(j) - 100) < 5){
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					p_t5wg_MET_signal_1900_100_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_HT_signal_1900_100_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_PhoEt_signal_1900_100_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_Mt_signal_1900_100_mg->Scale(lumi_mg*crosssection/noe);
					p_t5wg_MET_signal_1900_100_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_HT_signal_1900_100_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_PhoEt_signal_1900_100_eg->Scale(lumi_eg*crosssection/noe);
					p_t5wg_Mt_signal_1900_100_eg->Scale(lumi_eg*crosssection/noe);

					for(int s(1); s <= NBIN*2; s++){
                				if(s <= NBIN)   lumi = lumi_mg;
                				else            lumi = lumi_eg;
						p_t5wg_1900_100_nom->SetBinError(s, sqrt(p_t5wg_1900_100_nom->GetBinContent(s))*lumi*crosssection/noe);
						p_t5wg_1900_100_nom->SetBinContent(s, p_t5wg_1900_100_nom->GetBinContent(s)*lumi*crosssection/noe);}
				}
		}
	}

        double largeunc(0),smallunc(2);
	for(int ih(0); ih < NBIN*2; ih++){
    		for(int i(1); i <= t5wg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins(); i++){
        		for(int j(1); j <= t5wg_h_chan_rate_nom[ih]->GetYaxis()->GetNbins(); j++){
                                if(p_T5WGMASS->GetBinContent(i,j) > 1000){
                                        //if(t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j) > 0.1)std::cout << " !!" << std::endl;
                                        //std::cout << "pdf ch " << ih << " mass " << t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j) << std::endl;
                                        if(t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j) > largeunc)largeunc=t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j);
                                        if(t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j) >= 0 &&t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j) < smallunc)smallunc = largeunc=t5wg_h_chan_syserr_qcdscale[ih]->GetBinContent(i,j);
                                }
                        }
                }
        }
        std::cout << "large " << largeunc << std::endl;
        std::cout << "small " << smallunc << std::endl;

	outputfile_t5wg->Write();
	outputfile_t5wg->Close();


	// ****************   TChiWG ***************************
	TFile *file_tchiwg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_%d%s.root",RunYear, whichVFP.c_str()));
  	TTree *tree_tchiwg = (TTree*)file_tchiwg->Get("SUSYtree");
  	float Mchargino_tchiwg(0);
  	float Mneutralino_tchiwg(0);
	int   nVertex_tchiwg(0);
	std::vector<float> *ScaleSystWeight_tchiwg = 0;

  	tree_tchiwg->SetBranchAddress("Mchargino",  &Mchargino_tchiwg);
  	tree_tchiwg->SetBranchAddress("Mneutralino",&Mneutralino_tchiwg);
	tree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg);
	tree_tchiwg->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight_tchiwg);

	TFile *outputfile_tchiwg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/Tree_TChiWG_%d%s.root",RunYear, whichVFP.c_str()),"RECREATE");
	outputfile_tchiwg->cd();

	TH1D *p_TChiWGMASS = new TH1D("p_TChiWGMASS","",40,287.5,1287.5);
	TH1D *p_TChiWGMASS_pdf[9];
	for(unsigned i(0); i < 9; i++){
		histname.str("");
		histname << "p_TChiWGMASS_pdf" << i;
		p_TChiWGMASS_pdf[i] = new TH1D(histname.str().c_str(), histname.str().c_str(), 40, 287.5, 1287.5);
	}

	TH1D *p_tchiwg_MET_signal_eg = new TH1D("p_tchiwg_MET_signal_eg","",nSigMETBins, sigMETBins);
	TH1D *p_tchiwg_HT_signal_eg  = new TH1D("p_tchiwg_HT_signal_eg", "",nSigHTBins, sigHTBins);
	TH1D *p_tchiwg_PhoEt_signal_eg=new TH1D("p_tchiwg_PhoEt_signal_eg","",nSigEtBins, sigEtBins);
	TH1D *p_tchiwg_Mt_signal_eg = new TH1D("p_tchiwg_Mt_signal_eg","",nSigMtBins, sigMtBins);
	TH1D *p_tchiwg_MET_signal_mg = new TH1D("p_tchiwg_MET_signal_mg","",nSigMETBins, sigMETBins);
	TH1D *p_tchiwg_HT_signal_mg  = new TH1D("p_tchiwg_HT_signal_mg", "",nSigHTBins, sigHTBins);
	TH1D *p_tchiwg_PhoEt_signal_mg=new TH1D("p_tchiwg_PhoEt_signal_mg","",nSigEtBins, sigEtBins);
	TH1D *p_tchiwg_Mt_signal_mg = new TH1D("p_tchiwg_Mt_signal_mg","",nSigMtBins, sigMtBins);
	TH1D *p_tchiwg_800_nom = new TH1D("p_tchiwg_800_nom", "p_tchiwg_800_nom", NBIN*2,0,NBIN*2);

	for(int ievt(0); ievt < tree_tchiwg->GetEntries(); ievt++){
		tree_tchiwg->GetEntry(ievt);
    		if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(Mchargino_tchiwg <= 0) continue;
		p_TChiWGMASS->Fill(Mchargino_tchiwg);
		for(unsigned i(0); i < NPDF; i++) p_TChiWGMASS_pdf[i]->Fill(Mchargino_tchiwg, (*ScaleSystWeight_tchiwg)[i]);
	}

  //p_T5WGMASS->Draw();
  	TH1D *tchiwg_h_chan_rate_nom    [NBIN*2];
	TH1D *tchiwg_h_chan_rate_puUp   [NBIN*2];
	TH1D *tchiwg_h_chan_rate_puDown [NBIN*2];
	TH1D *tchiwg_h_chan_rate_jesUp  [NBIN*2];
	TH1D *tchiwg_h_chan_rate_jesDown[NBIN*2];
	TH1D *tchiwg_h_chan_rate_jerUp  [NBIN*2];
	TH1D *tchiwg_h_chan_rate_jerDown[NBIN*2];
	TH1D *tchiwg_h_chan_rate_xsUp   [NBIN*2];
	TH1D *tchiwg_h_chan_rate_esfUp  [NBIN*2];
 	TH1D *tchiwg_h_chan_rate_qcdscale    [NBIN*2][9];

	TH1D *tchiwg_h_syserr_PU = new TH1D("tchiwg_h_syserr_PU","",40,287.5,1287.5);
	TH1D *tchiwg_h_chan_syserr_jes     [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_jer     [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_esf     [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_scale   [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_eleshape[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_jetshape[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_qcdshape[NBIN*2];
	TH1D *tchiwg_h_chan_syserr_xs      [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_lumi    [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_isr     [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_pu      [NBIN*2];
	TH1D *tchiwg_h_chan_syserr_qcdscale     [NBIN*2];

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
                histname << "h_chan" << i+1 << "_rate_puUp";
                tchiwg_h_chan_rate_puUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
                histname.str("");
                histname << "h_chan" << i+1 << "_rate_puDown";
                tchiwg_h_chan_rate_puDown[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
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
		histname.str("");
                histname << "h_chan" << i+1 << "_syserr_pu";
                tchiwg_h_chan_syserr_pu[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		for(unsigned is(0); is < 9; is++){
			histname.str("");
			histname << "tchiwg_h_chan" << i+1 << "_rate_qcdscale" << is;
			tchiwg_h_chan_rate_qcdscale[i][is] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		}
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdscale";
		tchiwg_h_chan_syserr_qcdscale[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
	}

  TChain *mgtree_tchiwg;
  mgtree_tchiwg = new TChain("mgTree","mgTree");
  mgtree_tchiwg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_%d%s.root",RunYear, whichVFP.c_str()));
  float phoEt_tchiwg_mg(0);
  float phoEta_tchiwg_mg(0);
  float lepPt_tchiwg_mg(0);
  float lepEta_tchiwg_mg(0);
  float sigMT_tchiwg_mg(0);
  float sigMET_tchiwg_mg(0);
  float HT_tchiwg_mg(0);
  float nJet_tchiwg_mg(0);
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
  float charginoMass_tchiwg_mg(0);
  float neutralinoMass_tchiwg_mg(0);
  std::vector<float> *ScaleSystWeight_tchiwg_mg = 0;
  float PUweight_tchiwg_mg(1), PUweightUp_tchiwg_mg(1), PUweightDown_tchiwg_mg(1); // NEW
  mgtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("nJet",         &nJet_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("Mchargino",  &charginoMass_tchiwg_mg);
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
	mgtree_tchiwg->SetBranchAddress("PUweight",     &PUweight_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("PUweightUp",   &PUweightUp_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("PUweightDown", &PUweightDown_tchiwg_mg);
 	mgtree_tchiwg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_tchiwg_mg);

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
			p_tchiwg_MET_signal_mg->Fill( sigMET_tchiwg_mg );
			p_tchiwg_Mt_signal_mg->Fill( sigMT_tchiwg_mg );
			p_tchiwg_HT_signal_mg->Fill( HT_tchiwg_mg );
			p_tchiwg_PhoEt_signal_mg->Fill( phoEt_tchiwg_mg );
		}

		if(sigMET_tchiwg_mg < 120 || sigMT_tchiwg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);
		if(SigBinIndex >= 0){
			tchiwg_h_chan_rate_nom   [SigBinIndex]->Fill(charginoMass_tchiwg_mg, scalefactor * PUweight_tchiwg_mg);
			tchiwg_h_chan_rate_puUp  [SigBinIndex]->Fill(charginoMass_tchiwg_mg, scalefactor * PUweightUp_tchiwg_mg);   // NEW
			tchiwg_h_chan_rate_puDown[SigBinIndex]->Fill(charginoMass_tchiwg_mg, scalefactor * PUweightDown_tchiwg_mg); // NEW
			tchiwg_h_chan_rate_xsUp  [SigBinIndex]->Fill(charginoMass_tchiwg_mg, scalefactor * PUweight_tchiwg_mg);
			tchiwg_h_chan_rate_esfUp [SigBinIndex]->Fill(charginoMass_tchiwg_mg, scalefactorup * PUweight_tchiwg_mg);
			if(fabs(charginoMass_tchiwg_mg-800)<2) p_tchiwg_800_nom->Fill(SigBinIndex, scalefactor*PUweight_tchiwg_mg);
			for(unsigned is(0); is < 9; is++)
				tchiwg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(charginoMass_tchiwg_mg, scalefactor*(*ScaleSystWeight_tchiwg_mg)[is]);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_tchiwg_mg, HTJESup_tchiwg_mg, phoEt_tchiwg_mg);
		if(rawJesup  >= 0) tchiwg_h_chan_rate_jesUp  [rawJesup ]->Fill(charginoMass_tchiwg_mg, scalefactor*PUweight_tchiwg_mg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_tchiwg_mg, HTJESdo_tchiwg_mg, phoEt_tchiwg_mg);
		if(rawJesdo  >= 0) tchiwg_h_chan_rate_jesDown[rawJesdo ]->Fill(charginoMass_tchiwg_mg, scalefactor*PUweight_tchiwg_mg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);
		if(rawJerup  >= 0) tchiwg_h_chan_rate_jerUp  [rawJerup ]->Fill(charginoMass_tchiwg_mg, scalefactor*PUweight_tchiwg_mg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_tchiwg_mg, HT_tchiwg_mg, phoEt_tchiwg_mg);
		if(rawJerdo  >= 0) tchiwg_h_chan_rate_jerDown[rawJerdo ]->Fill(charginoMass_tchiwg_mg, scalefactor*PUweight_tchiwg_mg);
	}
  TChain *egtree_tchiwg;
  egtree_tchiwg = new TChain("egTree","egTree");
  egtree_tchiwg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_%d%s.root",RunYear, whichVFP.c_str()));
  float phoEt_tchiwg_eg(0);
  float phoEta_tchiwg_eg(0);
  float lepPt_tchiwg_eg(0);
  float lepEta_tchiwg_eg(0);
  float sigMT_tchiwg_eg(0);
  float sigMET_tchiwg_eg(0);
  float HT_tchiwg_eg(0);
  float nJet_tchiwg_eg(0);
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
  	float charginoMass_tchiwg_eg(0);
  	float neutralinoMass_tchiwg_eg(0);
	std::vector<float> *ScaleSystWeight_tchiwg_eg = 0;
	float PUweight_tchiwg_eg(1), PUweightUp_tchiwg_eg(1), PUweightDown_tchiwg_eg(1);
  egtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("nJet",         &nJet_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("nVertex",    &nVertex_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("Mchargino",  &charginoMass_tchiwg_eg);
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
	egtree_tchiwg->SetBranchAddress("PUweight",     &PUweight_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("PUweightUp",   &PUweightUp_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("PUweightDown", &PUweightDown_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_tchiwg_eg);

	for(int ievt(0); ievt < egtree_tchiwg->GetEntries(); ievt++){
		egtree_tchiwg->GetEntry(ievt);

		//  cut flow *****
		if(phoEt_tchiwg_eg < 35 || lepPt_tchiwg_eg < 25)continue;
		if(fabs(phoEta_tchiwg_eg) > 1.4442 || fabs(lepEta_tchiwg_eg) > 2.5)continue;

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
			p_tchiwg_MET_signal_eg->Fill( sigMET_tchiwg_eg );
			p_tchiwg_Mt_signal_eg->Fill( sigMT_tchiwg_eg );
			p_tchiwg_HT_signal_eg->Fill( HT_tchiwg_eg );
			p_tchiwg_PhoEt_signal_eg->Fill( phoEt_tchiwg_eg );
		}
		// signal-region selection applied here (after the diagnostic fills above),
		// matching the T5WG eg / T6WG eg ordering
		if(sigMET_tchiwg_eg < 120 || sigMT_tchiwg_eg < 100)continue;
		// FIX: guard on raw index BEFORE adding +NBIN offset
		int rawSigBin = Bin.findSignalBin(sigMET_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg);
		if(rawSigBin >= 0){
			int SigBinIndex = rawSigBin + NBIN;
			tchiwg_h_chan_rate_nom   [SigBinIndex]->Fill(charginoMass_tchiwg_eg, scalefactor * PUweight_tchiwg_eg);
			tchiwg_h_chan_rate_puUp  [SigBinIndex]->Fill(charginoMass_tchiwg_eg, scalefactor * PUweightUp_tchiwg_eg);
			tchiwg_h_chan_rate_puDown[SigBinIndex]->Fill(charginoMass_tchiwg_eg, scalefactor * PUweightDown_tchiwg_eg);
			tchiwg_h_chan_rate_xsUp  [SigBinIndex]->Fill(charginoMass_tchiwg_eg, scalefactor * PUweight_tchiwg_eg);
			tchiwg_h_chan_rate_esfUp [SigBinIndex]->Fill(charginoMass_tchiwg_eg, scalefactorup * PUweight_tchiwg_eg);
			if(fabs(charginoMass_tchiwg_eg-800)<2) p_tchiwg_800_nom->Fill(SigBinIndex, scalefactor*PUweight_tchiwg_eg);
			for(unsigned is(0); is < 9; is++)
				tchiwg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(charginoMass_tchiwg_eg, scalefactor*(*ScaleSystWeight_tchiwg_eg)[is]);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_tchiwg_eg, HTJESup_tchiwg_eg, phoEt_tchiwg_eg);
		if(rawJesup  >= 0) tchiwg_h_chan_rate_jesUp  [rawJesup+NBIN]->Fill(charginoMass_tchiwg_eg, scalefactor*PUweight_tchiwg_eg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_tchiwg_eg, HTJESdo_tchiwg_eg, phoEt_tchiwg_eg);
		if(rawJesdo  >= 0) tchiwg_h_chan_rate_jesDown[rawJesdo+NBIN]->Fill(charginoMass_tchiwg_eg, scalefactor*PUweight_tchiwg_eg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg);
		if(rawJerup  >= 0) tchiwg_h_chan_rate_jerUp  [rawJerup+NBIN]->Fill(charginoMass_tchiwg_eg, scalefactor*PUweight_tchiwg_eg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_tchiwg_eg, HT_tchiwg_eg, phoEt_tchiwg_eg);
		if(rawJerdo  >= 0) tchiwg_h_chan_rate_jerDown[rawJerdo+NBIN]->Fill(charginoMass_tchiwg_eg, scalefactor*PUweight_tchiwg_eg);
	}

	for(int ih(0); ih < NBIN*2; ih++){
		if(ih < NBIN)   lumi = lumi_mg;
                else            lumi = lumi_eg;
		for(int i(1); i < tchiwg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
				
			if(p_TChiWGMASS->GetBinContent(i) <= 0){
				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_rate_puUp  [ih]->SetBinContent(i,-1);
				tchiwg_h_chan_rate_puDown[ih]->SetBinContent(i,-1);
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
				tchiwg_h_chan_syserr_pu  [ih]->SetBinContent(i,-1);
				tchiwg_h_chan_syserr_qcdscale [ih]->SetBinContent(i,-1);
				if(ih==0) tchiwg_h_syserr_PU->SetBinContent(i,-1);
			}
			else{ 
				float noe  = p_TChiWGMASS->GetBinContent(i);
				float mSp  = p_TChiWGMASS->GetXaxis()->GetBinCenter(i);
				float xs   = p_crosssection_tchiwg->GetBinContent(p_crosssection_tchiwg->FindBin(mSp));
				float xsUp = xs + p_crosssection_tchiwg->GetBinError(p_crosssection_tchiwg->FindBin(mSp));
				float sc   = lumi*xs/noe;

				tchiwg_h_chan_rate_nom[ih]->SetBinError  (i, sqrt(tchiwg_h_chan_rate_nom[ih]->GetBinContent(i))*sc);
				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)*sc);
				if(ih%3 == 0){
					// Use noepdf[0] (nominal scale, muR=muF=1) as the common normalization
					// denominator for both acc_nom and acc_var, so that both use the same
					// scale-weighted event count and the comparison is apples-to-apples.
					double noe0 = p_TChiWGMASS_pdf[0]->GetBinContent(i);
					if(noe0 <= 0) goto tchiwg_qcdscale_done;

					{
					double acc_nom = 0;
					for(unsigned binHT(0); binHT < 3; binHT++)
						acc_nom += tchiwg_h_chan_rate_qcdscale[ih+binHT][0]->GetBinContent(i) / noe0;

					// Use only the 6 envelope-defining scale indices (skip redundant variations)
					const int scaleIdx[] = {1, 2, 3, 4, 6, 8};
					double pdfdiff(0);
					for(int is : scaleIdx){
						double noepdf = p_TChiWGMASS_pdf[is]->GetBinContent(i);
						if(noepdf <= 0) continue;
						double acc_var = 0;
						for(unsigned binHT(0); binHT < 3; binHT++)
							acc_var += tchiwg_h_chan_rate_qcdscale[ih+binHT][is]->GetBinContent(i) / noepdf;
						pdfdiff = std::max(pdfdiff, fabs(acc_var - acc_nom));
					}

					bool lowStats = (noe0 < 3.0 || acc_nom <= 0);

					for(unsigned binHT(0); binHT < 3; binHT++){
						if(tchiwg_h_chan_rate_nom[ih+binHT]->GetBinContent(i) > 0){
							if(!lowStats)
								tchiwg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i, pdfdiff / acc_nom);
							else
								tchiwg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i, 0.30);
						}
					}
					}
					tchiwg_qcdscale_done:;
				}
				tchiwg_h_chan_rate_puUp  [ih]->SetBinContent(i, tchiwg_h_chan_rate_puUp  [ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_puDown[ih]->SetBinContent(i, tchiwg_h_chan_rate_puDown[ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_jesUp [ih]->SetBinContent(i, tchiwg_h_chan_rate_jesUp [ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_jesDown[ih]->SetBinContent(i,tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_jerUp [ih]->SetBinContent(i, tchiwg_h_chan_rate_jerUp [ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_jerDown[ih]->SetBinContent(i,tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i)*sc);
				tchiwg_h_chan_rate_xsUp  [ih]->SetBinContent(i, tchiwg_h_chan_rate_xsUp  [ih]->GetBinContent(i)*lumi*xsUp/noe);
				tchiwg_h_chan_rate_esfUp [ih]->SetBinContent(i, tchiwg_h_chan_rate_esfUp [ih]->GetBinContent(i)*sc);

				float nom     = tchiwg_h_chan_rate_nom[ih]->GetBinContent(i);
                                float jesUp   = tchiwg_h_chan_rate_jesUp[ih]->GetBinContent(i);
                                float jesDown = tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i);
                                float jerUp   = tchiwg_h_chan_rate_jerUp[ih]->GetBinContent(i);
                                float jerDown = tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i);

                                float jesErr = std::max(fabs(jesUp-nom), fabs(jesDown-nom));
                                float jerErr = std::max(fabs(jerUp-nom), fabs(jerDown-nom));
                                const float LOW_STAT_YIELD = 2.0;
                                if (nom > 0.0 && nom < LOW_STAT_YIELD && (jesErr  > nom || jerErr > nom)) {
                                	float maxShift = std::max({nom, jesUp, jesDown, jerUp, jerDown});
                                        float newNom   = 0.5 * maxShift;
                                        tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, newNom);
                                        float statErr = tchiwg_h_chan_rate_nom[ih]->GetBinError(i);
                                        tchiwg_h_chan_rate_nom[ih]->SetBinError  (i, statErr * newNom / nom);
                                        std::cout << "[JES/JER FIX] Year=" << RunYear << " bin=(" << i << ")"
                                                        << " nom=" << nom << " newNom=" << newNom << std::endl;
                                        nom = newNom;}
                                tchiwg_h_chan_syserr_jes[ih]->SetBinContent(i, std::max(fabs(jesUp-nom), fabs(jesDown-nom)));
                                tchiwg_h_chan_syserr_jer[ih]->SetBinContent(i, std::max(fabs(jerUp-nom), fabs(jerDown-nom)));

				tchiwg_h_chan_syserr_esf[ih]->SetBinContent(i, fabs(tchiwg_h_chan_rate_esfUp[ih]->GetBinContent(i)-nom));
				tchiwg_h_chan_syserr_scale[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_syserr_eleshape[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_syserr_jetshape[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_syserr_xs  [ih]->SetBinContent(i, fabs(tchiwg_h_chan_rate_xsUp[ih]->GetBinContent(i)-nom));
				tchiwg_h_chan_syserr_lumi[ih]->SetBinContent(i, 0.026*nom);
				tchiwg_h_chan_syserr_isr [ih]->SetBinContent(i,-1);

				float puUp_sc   = tchiwg_h_chan_rate_puUp  [ih]->GetBinContent(i);
                                float puDown_sc = tchiwg_h_chan_rate_puDown[ih]->GetBinContent(i);
                                float puErr     = std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom));
                                float flat_pu_err = flat_pu_unc * nom;
                                float final_pu_err;
                                if(nom >= 1.0)  final_pu_err = std::min(puErr, 3.0f * flat_pu_err);
                                else final_pu_err = flat_pu_err;
                                tchiwg_h_chan_syserr_pu[ih]->SetBinContent(i, final_pu_err);
				if(ih==0) tchiwg_h_syserr_PU->SetBinContent(i, std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom)));
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
				p_tchiwg_Mt_signal_eg->Scale(lumi_eg*crosssection/noe);
				p_tchiwg_MET_signal_mg->Scale(lumi_mg*crosssection/noe);
				p_tchiwg_HT_signal_mg->Scale(lumi_mg*crosssection/noe);
				p_tchiwg_PhoEt_signal_mg->Scale(lumi_mg*crosssection/noe);
				p_tchiwg_Mt_signal_mg->Scale(lumi_mg*crosssection/noe);

				for(int k(1); k <= NBIN*2; k++){
                			if(k <= NBIN)   lumi = lumi_mg;
                			else            lumi = lumi_eg;
					p_tchiwg_800_nom->SetBinError(k, sqrt(p_tchiwg_800_nom->GetBinContent(k))*lumi*crosssection/noe);
					p_tchiwg_800_nom->SetBinContent(k, p_tchiwg_800_nom->GetBinContent(k)*lumi*crosssection/noe);}

			}
		} 

	outputfile_tchiwg->Write();
	outputfile_tchiwg->Close();  

	// ************   T6WG  ***************************

  TFile *file_t6wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T6WG_%d%s.root",RunYear, whichVFP.c_str()));
  TTree *tree_t6wg = (TTree*)file_t6wg->Get("SUSYtree");
	float Msquark_t6wg(0);
  float Mchargino_t6wg(0);
  float Mneutralino_t6wg(0);
  	std::vector<float> *ScaleSystWeight_t6wg = 0;
	tree_t6wg->SetBranchAddress("MsGsQ",    &Msquark_t6wg);  
  tree_t6wg->SetBranchAddress("Mchargino",  &Mchargino_t6wg);
  tree_t6wg->SetBranchAddress("Mneutralino",&Mneutralino_t6wg);
	tree_t6wg->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight_t6wg);
	TFile *outputfile_t6wg = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/Tree_T6WG_%d%s.root",RunYear, whichVFP.c_str()),"RECREATE");
	outputfile_t6wg->cd();

	TH2D *p_T6WGMASS = new TH2D("T6WGSUSYMass","",nXbins, xBins, nYbins, yBins);
	TH2D *p_T6WGselect = new TH2D("p_T6WGselect","",nXbins, xBins, nYbins, yBins); 
	
	TH2D *p_T6WGMASS_pdf[9];
        for(unsigned i(0); i < 9; i++){
                histname.str("");
                histname << "p_T6WGMASS_pdf" << i;
                p_T6WGMASS_pdf[i] = new TH2D(histname.str().c_str(), histname.str().c_str(), nXbins,xBins,nYbins,yBins);
        }

	for(int ievt(0); ievt < tree_t6wg->GetEntries(); ievt++){
		tree_t6wg->GetEntry(ievt);
    		//if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
    		double NLSPMass(0);
    		if(Mchargino_t6wg >0)NLSPMass=Mchargino_t6wg;
		else if(Mneutralino_t6wg >0)NLSPMass=Mneutralino_t6wg;
		if(NLSPMass <= 0)continue;
		p_T6WGMASS->Fill(Msquark_t6wg, NLSPMass);
		for(unsigned i(0); i < NPDF; i++)p_T6WGMASS_pdf[i]->Fill(Msquark_t6wg, NLSPMass, (*ScaleSystWeight_t6wg)[i]);
	}
  //p_T6WGMASS->Draw();
	TH2D *t6wg_h_chan_rate_nom[NBIN*2]; 
	TH2D *t6wg_h_chan_rate_puUp   [NBIN*2];
	TH2D *t6wg_h_chan_rate_puDown [NBIN*2];
	TH2D *t6wg_h_chan_rate_jesUp[NBIN*2]; 
	TH2D *t6wg_h_chan_rate_jesDown[NBIN*2];    
	TH2D *t6wg_h_chan_rate_jerUp[NBIN*2];    
	TH2D *t6wg_h_chan_rate_jerDown[NBIN*2];    
	TH2D *t6wg_h_chan_rate_xsUp[NBIN*2];       
	TH2D *t6wg_h_chan_rate_esfUp[NBIN*2]; 
        TH2D *t6wg_h_chan_rate_qcdscale[NBIN*2][9];	
                                  
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
	TH2D *t6wg_h_chan_syserr_pu      [NBIN*2];
	TH2D *t6wg_h_chan_syserr_qcdscale[NBIN*2];
	TH2D *t6wg_h_syserr_PU = new TH2D("t6wg_h_syserr_PU","",nXbins,xBins,nYbins,yBins);

	for(int i(0); i < NBIN*2; i++){
		 for(unsigned is(0); is<9; is++){
                        histname.str("");
                        histname << "t6wg_h_chan" << i+1 << "_rate_qcdscale" << is;
                        t6wg_h_chan_rate_qcdscale[i][is] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
                }
		histname.str("");
                histname << "h_chan" << i+1 << "_syserr_qcdscale";
                t6wg_h_chan_syserr_qcdscale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(), nXbins, xBins, nYbins, yBins);

		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		t6wg_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		t6wg_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		t6wg_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		t6wg_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		t6wg_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		t6wg_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		t6wg_h_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
			
		histname.str("");
                histname << "h_chan" << i+1 << "_rate_puUp";
                t6wg_h_chan_rate_puUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
                histname << "h_chan" << i+1 << "_rate_puDown";
                t6wg_h_chan_rate_puDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		t6wg_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		t6wg_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		t6wg_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		t6wg_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		t6wg_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		t6wg_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		t6wg_h_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		t6wg_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		t6wg_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		t6wg_h_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);
		histname.str("");
                histname << "h_chan" << i+1 << "_syserr_pu";
                t6wg_h_chan_syserr_pu[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),nXbins, xBins, nYbins, yBins);

	}
	// ---- T6WG(1600,1500), added following the same pattern as T5WG's single-mass-point histograms ----
	TH1D *p_t6wg_MET_signal_1600_1500_mg = new TH1D("p_t6wg_MET_signal_1600_1500_mg","", nSigMETBins, sigMETBins);
	TH1D *p_t6wg_HT_signal_1600_1500_mg  = new TH1D("p_t6wg_HT_signal_1600_1500_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_t6wg_PhoEt_signal_1600_1500_mg = new TH1D("p_t6wg_PhoEt_signal_1600_1500_mg","",nSigEtBins, sigEtBins);
	TH1D *p_t6wg_Mt_signal_1600_1500_mg = new TH1D("p_t6wg_Mt_signal_1600_1500_mg","",nSigMtBins, sigMtBins);
	TH1D *p_t6wg_MET_signal_1600_1500_eg = new TH1D("p_t6wg_MET_signal_1600_1500_eg","", nSigMETBins, sigMETBins);
	TH1D *p_t6wg_HT_signal_1600_1500_eg  = new TH1D("p_t6wg_HT_signal_1600_1500_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_t6wg_PhoEt_signal_1600_1500_eg = new TH1D("p_t6wg_PhoEt_signal_1600_1500_eg","",nSigEtBins, sigEtBins);
	TH1D *p_t6wg_Mt_signal_1600_1500_eg = new TH1D("p_t6wg_Mt_signal_1600_1500_eg","",nSigMtBins, sigMtBins);
	TH1D *p_t6wg_1600_1500_nom = new TH1D("p_t6wg_1600_1500_nom", "p_t6wg_1600_1500_nom", NBIN*2,0,NBIN*2);

  TChain *mgtree_t6wg;
  mgtree_t6wg = new TChain("mgTree","mgTree");
  mgtree_t6wg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T6WG_%d%s.root",RunYear, whichVFP.c_str()));
  float phoEt_t6wg_mg(0);
  float phoEta_t6wg_mg(0);
  float lepPt_t6wg_mg(0);
  float lepEta_t6wg_mg(0);
  float sigMT_t6wg_mg(0);
  float sigMET_t6wg_mg(0);
  float HT_t6wg_mg(0);
  float nJet_t6wg_mg(0);
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
  	std::vector<float> *ScaleSystWeight_t6wg_mg = 0;
  	float PUweight_t6wg_mg(1), PUweightUp_t6wg_mg(1), PUweightDown_t6wg_mg(1);
  mgtree_t6wg->SetBranchAddress("phoEt",      &phoEt_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("phoEta",     &phoEta_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("lepPt",      &lepPt_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("lepEta",     &lepEta_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("sigMT",      &sigMT_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("sigMET",     &sigMET_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("HT",         &HT_t6wg_mg);
  mgtree_t6wg->SetBranchAddress("nJet",         &nJet_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("MsGsQ",    &squarkMass_t6wg_mg);
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
		mgtree_t6wg->SetBranchAddress("PUweight",     &PUweight_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("PUweightUp",   &PUweightUp_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("PUweightDown", &PUweightDown_t6wg_mg);
	mgtree_t6wg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_t6wg_mg);
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

		if( fabs(squarkMass_t6wg_mg - 1600) < 10 && fabs(NLSPMass - 1500) < 5){
			p_t6wg_MET_signal_1600_1500_mg->Fill(sigMET_t6wg_mg);
			p_t6wg_Mt_signal_1600_1500_mg->Fill(sigMT_t6wg_mg);
			p_t6wg_HT_signal_1600_1500_mg->Fill(HT_t6wg_mg);
			p_t6wg_PhoEt_signal_1600_1500_mg->Fill(phoEt_t6wg_mg);
		}

		if(sigMET_t6wg_mg < 120 || sigMT_t6wg_mg < 100)continue;
    p_T6WGselect->Fill(squarkMass_t6wg_mg, NLSPMass);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg);
		if(SigBinIndex >= 0){
			t6wg_h_chan_rate_nom   [SigBinIndex]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor * PUweight_t6wg_mg);
			t6wg_h_chan_rate_puUp  [SigBinIndex]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor * PUweightUp_t6wg_mg);
			t6wg_h_chan_rate_puDown[SigBinIndex]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor * PUweightDown_t6wg_mg);
			t6wg_h_chan_rate_xsUp  [SigBinIndex]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor * PUweight_t6wg_mg);
			t6wg_h_chan_rate_esfUp [SigBinIndex]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactorup * PUweight_t6wg_mg);
			for(unsigned is(0); is < 9; is++)
                        	t6wg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor*(*ScaleSystWeight_t6wg_mg)[is]);

			if(fabs(squarkMass_t6wg_mg-1600)<10 && fabs(NLSPMass-1500)<5)
				p_t6wg_1600_1500_nom->Fill(SigBinIndex, scalefactor * PUweight_t6wg_mg);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_t6wg_mg, HTJESup_t6wg_mg, phoEt_t6wg_mg);
		if(rawJesup >= 0) t6wg_h_chan_rate_jesUp [rawJesup]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor*PUweight_t6wg_mg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_t6wg_mg, HTJESdo_t6wg_mg, phoEt_t6wg_mg);
		if(rawJesdo >= 0) t6wg_h_chan_rate_jesDown[rawJesdo]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor*PUweight_t6wg_mg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg);
		if(rawJerup >= 0) t6wg_h_chan_rate_jerUp [rawJerup]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor*PUweight_t6wg_mg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_t6wg_mg, HT_t6wg_mg, phoEt_t6wg_mg);
		if(rawJerdo >= 0) t6wg_h_chan_rate_jerDown[rawJerdo]->Fill(squarkMass_t6wg_mg, NLSPMass, scalefactor*PUweight_t6wg_mg);
	}

  TChain *egtree_t6wg;
  egtree_t6wg = new TChain("egTree","egTree");
  egtree_t6wg->Add(Form("/uscms/home/tmishra/nobackup/signal_trees/resTree_T6WG_%d%s.root",RunYear, whichVFP.c_str()));
  float phoEt_t6wg_eg(0);
  float phoEta_t6wg_eg(0);
  float lepPt_t6wg_eg(0);
  float lepEta_t6wg_eg(0);
  float sigMT_t6wg_eg(0);
  	int   nVertex_t6wg_eg(0);
  float sigMET_t6wg_eg(0);
  float HT_t6wg_eg(0);
  float nJet_t6wg_eg(0);
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
  std::vector<float> *ScaleSystWeight_t6wg_eg = 0;
  float PUweight_t6wg_eg(1), PUweightUp_t6wg_eg(1), PUweightDown_t6wg_eg(1);

  egtree_t6wg->SetBranchAddress("phoEt",      &phoEt_t6wg_eg);
  egtree_t6wg->SetBranchAddress("phoEta",     &phoEta_t6wg_eg);
  egtree_t6wg->SetBranchAddress("lepPt",      &lepPt_t6wg_eg);
  egtree_t6wg->SetBranchAddress("lepEta",     &lepEta_t6wg_eg);
  egtree_t6wg->SetBranchAddress("sigMT",      &sigMT_t6wg_eg);
  egtree_t6wg->SetBranchAddress("sigMET",     &sigMET_t6wg_eg);
  egtree_t6wg->SetBranchAddress("HT",         &HT_t6wg_eg);
  egtree_t6wg->SetBranchAddress("nJet",         &nJet_t6wg_eg);
  	egtree_t6wg->SetBranchAddress("nVertex",    &nVertex_t6wg_eg);
	egtree_t6wg->SetBranchAddress("MsGsQ",    &squarkMass_t6wg_eg);
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
	egtree_t6wg->SetBranchAddress("PUweight",     &PUweight_t6wg_eg); 
	egtree_t6wg->SetBranchAddress("PUweightUp",   &PUweightUp_t6wg_eg); 
	egtree_t6wg->SetBranchAddress("PUweightDown", &PUweightDown_t6wg_eg);
	egtree_t6wg->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_t6wg_eg); 

	for(int ievt(0); ievt < egtree_t6wg->GetEntries(); ievt++){
		egtree_t6wg->GetEntry(ievt);

    double NLSPMass(0);
    if(charginoMass_t6wg_eg >0)NLSPMass=charginoMass_t6wg_eg;
		else if(neutralinoMass_t6wg_eg >0)NLSPMass=neutralinoMass_t6wg_eg;
		if(NLSPMass <= 0)continue;
		// ** cut flow *****
		if(phoEt_t6wg_eg < 35 || lepPt_t6wg_eg < 25)continue;
		if(fabs(phoEta_t6wg_eg) > 1.4442 || fabs(lepEta_t6wg_eg) > 2.5)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastegPhotonTRGESF(phoEt_t6wg_eg,phoEta_t6wg_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronESF(lepPt_t6wg_eg,lepEta_t6wg_eg)*objectESF.getPhotonESF(phoEt_t6wg_eg,phoEta_t6wg_eg)*objectESF.getFastElectronTRGESF(lepPt_t6wg_eg,lepEta_t6wg_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		if( fabs(squarkMass_t6wg_eg - 1600) < 10 && fabs(NLSPMass - 1500) < 5){
			p_t6wg_MET_signal_1600_1500_eg->Fill(sigMET_t6wg_eg);
			p_t6wg_Mt_signal_1600_1500_eg->Fill(sigMT_t6wg_eg);
			p_t6wg_HT_signal_1600_1500_eg->Fill(HT_t6wg_eg);
			p_t6wg_PhoEt_signal_1600_1500_eg->Fill(phoEt_t6wg_eg);
		}

		if(sigMET_t6wg_eg < 120 || sigMT_t6wg_eg < 100)continue;

		// FIX: guard on raw index BEFORE adding +NBIN offset
		int rawSigBin = Bin.findSignalBin(sigMET_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg);
		if(rawSigBin >= 0){
			int SigBinIndex = rawSigBin + NBIN;
			t6wg_h_chan_rate_nom   [SigBinIndex]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor * PUweight_t6wg_eg);
			t6wg_h_chan_rate_puUp  [SigBinIndex]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor * PUweightUp_t6wg_eg);
			t6wg_h_chan_rate_puDown[SigBinIndex]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor * PUweightDown_t6wg_eg);
			t6wg_h_chan_rate_xsUp  [SigBinIndex]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor * PUweight_t6wg_eg);
			t6wg_h_chan_rate_esfUp [SigBinIndex]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactorup * PUweight_t6wg_eg);
			for(unsigned is(0); is < 9; is++)
                                t6wg_h_chan_rate_qcdscale[SigBinIndex][is]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor*(*ScaleSystWeight_t6wg_eg)[is]);

			if(fabs(squarkMass_t6wg_eg-1600)<10 && fabs(NLSPMass-1500)<5)
				p_t6wg_1600_1500_nom->Fill(SigBinIndex, scalefactor * PUweight_t6wg_eg);
		}
		int rawJesup = Bin.findSignalBin(sigMETJESup_t6wg_eg, HTJESup_t6wg_eg, phoEt_t6wg_eg);
		if(rawJesup >= 0) t6wg_h_chan_rate_jesUp [rawJesup+NBIN]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor*PUweight_t6wg_eg);
		int rawJesdo = Bin.findSignalBin(sigMETJESdo_t6wg_eg, HTJESdo_t6wg_eg, phoEt_t6wg_eg);
		if(rawJesdo >= 0) t6wg_h_chan_rate_jesDown[rawJesdo+NBIN]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor*PUweight_t6wg_eg);
		int rawJerup = Bin.findSignalBin(sigMETJERup_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg);
		if(rawJerup >= 0) t6wg_h_chan_rate_jerUp [rawJerup+NBIN]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor*PUweight_t6wg_eg);
		int rawJerdo = Bin.findSignalBin(sigMETJERdo_t6wg_eg, HT_t6wg_eg, phoEt_t6wg_eg);
		if(rawJerdo >= 0) t6wg_h_chan_rate_jerDown[rawJerdo+NBIN]->Fill(squarkMass_t6wg_eg, NLSPMass, scalefactor*PUweight_t6wg_eg);
	}

	// ---- T6WG(1600,1500) normalization, mirroring T5WG's dedicated single-mass-point block ----
	for(int i(1); i < p_T6WGMASS->GetXaxis()->GetNbins() + 1; i++){
		for(int j(1); j < p_T6WGMASS->GetYaxis()->GetNbins() + 1; j++){
			float sparticleMass = p_T6WGMASS->GetXaxis()->GetBinCenter(i);
			if( fabs(sparticleMass - 1600 ) < 10 && fabs( p_T6WGMASS->GetYaxis()->GetBinCenter(j) - 1500) < 5){
				float noe = p_T6WGMASS->GetBinContent(i,j);
				float crosssection = p_crosssection_t6wg->GetBinContent( p_crosssection_t6wg->FindBin(sparticleMass) );
				p_t6wg_MET_signal_1600_1500_mg->Scale(lumi_mg*crosssection/noe);
				p_t6wg_HT_signal_1600_1500_mg->Scale(lumi_mg*crosssection/noe);
				p_t6wg_PhoEt_signal_1600_1500_mg->Scale(lumi_mg*crosssection/noe);
				p_t6wg_Mt_signal_1600_1500_mg->Scale(lumi_mg*crosssection/noe);
				p_t6wg_MET_signal_1600_1500_eg->Scale(lumi_eg*crosssection/noe);
				p_t6wg_HT_signal_1600_1500_eg->Scale(lumi_eg*crosssection/noe);
				p_t6wg_PhoEt_signal_1600_1500_eg->Scale(lumi_eg*crosssection/noe);
				p_t6wg_Mt_signal_1600_1500_eg->Scale(lumi_eg*crosssection/noe);

				for(int s(1); s <= NBIN*2; s++){
                				if(s <= NBIN)   lumi = lumi_mg;
                				else            lumi = lumi_eg;
					p_t6wg_1600_1500_nom->SetBinError(s, sqrt(p_t6wg_1600_1500_nom->GetBinContent(s))*lumi*crosssection/noe);
					p_t6wg_1600_1500_nom->SetBinContent(s, p_t6wg_1600_1500_nom->GetBinContent(s)*lumi*crosssection/noe);}
			}
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
										t6wg_h_chan_rate_nom   [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_puUp  [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_puDown[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_jesUp [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_jerUp [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_rate_xsUp  [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_jes [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_jer [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_esf [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_scale[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_xs  [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_isr [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_pu  [ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_qcdscale[ih]->SetBinContent(i,j, -1);
					if(ih==0) t6wg_h_syserr_PU->SetBinContent(i,j,-1);

				}
				else{
					float noe  = p_T6WGMASS->GetBinContent(i,j);
					float mSp  = p_T6WGMASS->GetXaxis()->GetBinCenter(i);
					float xs   = p_crosssection_t6wg->GetBinContent(p_crosssection_t6wg->FindBin(mSp));
					float xsUp = xs + p_crosssection_t6wg->GetBinError(p_crosssection_t6wg->FindBin(mSp));
					float sc   = lumi*xs/noe;

					t6wg_h_chan_rate_nom[ih]->SetBinError  (i,j, sqrt(t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j))*sc);
					t6wg_h_chan_rate_nom[ih]->SetBinContent(i,j, t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j)*sc);
					if(ih%3 == 0){
						// Use noepdf[0] (nominal scale, muR=muF=1) as the common normalization
						// denominator for both acc_nom and acc_var, so that both use the same
						// scale-weighted event count and the comparison is apples-to-apples.
						double noe0 = p_T6WGMASS_pdf[0]->GetBinContent(i,j);
						if(noe0 <= 0) goto t6wg_qcdscale_done;

						{
						double acc_nom = 0;
						for(unsigned binHT(0); binHT < 3; binHT++)
							acc_nom += t6wg_h_chan_rate_qcdscale[ih+binHT][0]->GetBinContent(i,j) / noe0;

						// Use only the 6 envelope-defining scale indices (skip redundant variations)
						const int scaleIdx[] = {1, 2, 3, 4, 6, 8};
						double pdfdiff(0);
						for(int is : scaleIdx){
							double noepdf = p_T6WGMASS_pdf[is]->GetBinContent(i,j);
							if(noepdf <= 0) continue;
							double acc_var = 0;
							for(unsigned binHT(0); binHT < 3; binHT++)
								acc_var += t6wg_h_chan_rate_qcdscale[ih+binHT][is]->GetBinContent(i,j) / noepdf;
							pdfdiff = std::max(pdfdiff, fabs(acc_var - acc_nom));
						}

						// Low-stats guard: count raw selected events (index 0 qcdscale ≈ nominal)
						double N_raw_tot = 0;
						for(unsigned binHT(0); binHT < 3; binHT++)
							N_raw_tot += t6wg_h_chan_rate_qcdscale[ih+binHT][0]->GetBinContent(i,j);
						// N_raw_tot here is ESF-weighted; use noe0 (sum of ScaleWeight[0]) as proxy
						// for stats: if noe0 < 3 the point is low-stats
						bool lowStats = (noe0 < 3.0 || acc_nom <= 0);

						for(unsigned binHT(0); binHT < 3; binHT++){
							if(t6wg_h_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0){
								if(!lowStats)
									t6wg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i,j, pdfdiff / acc_nom);
								else
									t6wg_h_chan_syserr_qcdscale[ih+binHT]->SetBinContent(i,j, 0.30);
							}
						}
						}
						t6wg_qcdscale_done:;
					}

					t6wg_h_chan_rate_puUp  [ih]->SetBinContent(i,j, t6wg_h_chan_rate_puUp  [ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_puDown[ih]->SetBinContent(i,j, t6wg_h_chan_rate_puDown[ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_jesUp [ih]->SetBinContent(i,j, t6wg_h_chan_rate_jesUp [ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j,t6wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_jerUp [ih]->SetBinContent(i,j, t6wg_h_chan_rate_jerUp [ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,t6wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*sc);
					t6wg_h_chan_rate_xsUp  [ih]->SetBinContent(i,j, t6wg_h_chan_rate_xsUp  [ih]->GetBinContent(i,j)*lumi*xsUp/noe);
					t6wg_h_chan_rate_esfUp [ih]->SetBinContent(i,j, t6wg_h_chan_rate_esfUp [ih]->GetBinContent(i,j)*sc);

					float nom     = t6wg_h_chan_rate_nom[ih]->GetBinContent(i,j);
                                        float jesUp   = t6wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j);
                                        float jesDown = t6wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j);
                                        float jerUp   = t6wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j);
                                        float jerDown = t6wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j);

                                        float jesErr = std::max(fabs(jesUp-nom), fabs(jesDown-nom));
                                        float jerErr = std::max(fabs(jerUp-nom), fabs(jerDown-nom));

                                        const float LOW_STAT_YIELD = 2.0;

                                        if (nom > 0.0 && nom < LOW_STAT_YIELD && (jesErr  > nom || jerErr > nom)) {

                                                float maxShift = std::max({nom, jesUp, jesDown, jerUp, jerDown});
                                                float newNom   = 0.5 * maxShift;

                                                t6wg_h_chan_rate_nom[ih]->SetBinContent(i,j, newNom);
                                                float statErr = t6wg_h_chan_rate_nom[ih]->GetBinError(i,j);
                                                t6wg_h_chan_rate_nom[ih]->SetBinError  (i,j, statErr * newNom / nom);

                                                //std::cout << "[JES/JER FIX] Year=" << RunYear << " bin=(" << i << "," << j << ")"
                                                //       << " nom=" << nom << " newNom=" << newNom << std::endl;
                                                nom = newNom;
                                        }
                                        t6wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, std::max(fabs(jesUp-nom), fabs(jesDown-nom)));
                                        t6wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, std::max(fabs(jerUp-nom), fabs(jerDown-nom)));

					t6wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, fabs(t6wg_h_chan_rate_esfUp[ih]->GetBinContent(i,j)-nom));
					t6wg_h_chan_syserr_scale[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j,-1);
					t6wg_h_chan_syserr_xs  [ih]->SetBinContent(i,j, fabs(t6wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)-nom));
					t6wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*nom);
					t6wg_h_chan_syserr_isr [ih]->SetBinContent(i,j,-1);

					float puUp_sc   = t6wg_h_chan_rate_puUp  [ih]->GetBinContent(i,j);
                                        float puDown_sc = t6wg_h_chan_rate_puDown[ih]->GetBinContent(i,j);
                                        float puErr     = std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom));
                                        float flat_pu_err = flat_pu_unc * nom;
                                        float final_pu_err;
                                        if(nom >= 1.0)  final_pu_err = std::min(puErr, 3.0f * flat_pu_err);
                                        else final_pu_err = flat_pu_err;
                                        t6wg_h_chan_syserr_pu[ih]->SetBinContent(i,j, final_pu_err);
					if(ih==0) t6wg_h_syserr_PU->SetBinContent(i,j, std::max(fabs(puUp_sc-nom), fabs(puDown_sc-nom)));

				}
			}
		} 
	} 

	outputfile_t6wg->Write();
	outputfile_t6wg->Close(); 


  xSecFile.Close();
}
