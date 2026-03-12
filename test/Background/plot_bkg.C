#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<vector>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "../../include/analysis_commoncode.h"
// plottype  1 = bkg, 2 = valid

struct StackComponent {
	TH1F *hist;
	TString label;
};

std::vector<StackComponent> BuildIntegralOrderedStack(const std::vector<StackComponent>& components, const TString& prefix){
	std::vector<StackComponent> sorted(components);
	std::sort(sorted.begin(), sorted.end(), [](const StackComponent& a, const StackComponent& b){
		return a.hist->Integral() < b.hist->Integral();
	});

	std::vector<StackComponent> cumulative;
	TH1F *running = 0;
	for(size_t i = 0; i < sorted.size(); ++i){
		TH1F *histClone = (TH1F*)sorted[i].hist->Clone(Form("%s_%zu", prefix.Data(), i));
		histClone->SetDirectory(0);
		if(running) histClone->Add(running);
		running = histClone;
		cumulative.push_back({histClone, sorted[i].label});
	}
	return cumulative;
}

bool doTT = false;
void plot_bkg(int channel,int RunYear,bool preVFP){//main  
	SetRunConfig();
	int plottype = anatype;
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
  	gSystem->Load("../../lib/libAnaClasses.so");

	//Double_t bkgEtBins[]={35,40,50,60,70,80,90,100,110,120,130,140,150,160,170,185,200,215,230,250,275,290, 305,325,345,370,400,500,800};
	Double_t bkgEtBins[]={35,50,100,150,200,250,300,500,800};

	int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
	Double_t bkgPtBins[]={25,50,75,100,125,150,200,400,800};
	int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
	//Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
	Double_t bkgMETBins[]={0,40,60,80,100,120};
	int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
	Double_t bkgMtBins[]={0,20,40,60,80,100,120,140,160,180,200,300,400,500,1000};
	int nBkgMtBins= sizeof(bkgMtBins)/sizeof(bkgMtBins[0]) -1;
	Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
	int nBkgHTBins= sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	std::string whichChannel;
	if(channel == 1) whichChannel = "egamma";
	if(channel == 2) whichChannel = "mgamma";

	std::ostringstream signame;  signame.str("");
	std::ostringstream elename;  elename.str("");
	std::ostringstream jetname;  jetname.str("");
	std::ostringstream qcdname;  qcdname.str("");
	std::ostringstream VGname;   VGname.str("");
	std::ostringstream rarename; rarename.str(""); 

	if(channel == 1){
		if(plottype == 1){
			signame << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_signal_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			elename << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_eleBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			jetname << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_jetbkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			qcdname << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_qcd_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			VGname  << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_VGBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			rarename << "/uscms_data/d3/tmishra/Background/bkgTree_egamma_rareBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
		}
		else if(plottype == 2){
			signame << "/uscms_data/d3/tmishra/Background/validTree_egamma_signal_" << RunYear<<whichVFP <<".root";
			elename << "/uscms_data/d3/tmishra/Background/validTree_egamma_eleBkg_" << RunYear<<whichVFP <<".root";
			jetname << "/uscms_data/d3/tmishra/Background/validTree_egamma_jetbkg_" << RunYear<<whichVFP <<".root";
			qcdname << "/uscms_data/d3/tmishra/Background/validTree_egamma_qcd_" << RunYear<<whichVFP <<".root";
			VGname  << "/uscms_data/d3/tmishra/Background/validTree_egamma_VGBkg_" << RunYear<<whichVFP <<".root";
			rarename << "/uscms_data/d3/tmishra/Background/validTree_egamma_rareBkg_" << RunYear<<whichVFP <<".root";
		}
	}
	else if(channel == 2){
		if(plottype == 1){
			signame << "/uscms_data/d3/tmishra/Background/bkgTree_mg_signal_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			elename << "/uscms_data/d3/tmishra/Background/bkgTree_mg_eleBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			jetname << "/uscms_data/d3/tmishra/Background/bkgTree_mg_jetbkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			qcdname << "/uscms_data/d3/tmishra/Background/bkgTree_mg_qcd_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			VGname  << "/uscms_data/d3/tmishra/Background/bkgTree_mg_VGBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			rarename << "/uscms_data/d3/tmishra/Background/bkgTree_mg_rareBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
		}
		else if(plottype == 2){
			signame << "/uscms_data/d3/tmishra/Background/validTree_mg_signal_" << RunYear<<whichVFP <<".root";
			elename << "/uscms_data/d3/tmishra/Background/validTree_mg_eleBkg_" << RunYear<<whichVFP <<".root";
			jetname << "/uscms_data/d3/tmishra/Background/validTree_mg_jetbkg_" << RunYear<<whichVFP <<".root";
			qcdname << "/uscms_data/d3/tmishra/Background/validTree_mg_qcd_" << RunYear<<whichVFP <<".root";
			VGname  << "/uscms_data/d3/tmishra/Background/validTree_mg_VGBkg_" << RunYear<<whichVFP <<".root";
			rarename << "/uscms_data/d3/tmishra/Background/validTree_mg_rareBkg_" << RunYear<<whichVFP <<".root";
		}
	}

	TFile *file_sig = TFile::Open(signame.str().c_str());
	TFile *file_ele = TFile::Open(elename.str().c_str());
	TFile *file_jet = TFile::Open(jetname.str().c_str());
	TFile *file_qcd = TFile::Open(qcdname.str().c_str());
	TFile *file_VG  = TFile::Open(VGname.str().c_str());
	TFile *file_rare= TFile::Open(rarename.str().c_str()); 

	//TFile *file_t5 = TFile::Open("signalTree_T5WG.root");
	//TFile *file_tchi=TFile::Open("signalTree_TChiWG.root");

	TGraphErrors *error_nJet = new TGraphErrors(10);
	TGraphErrors *error_dPhiEleMET = new TGraphErrors(32);
	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_MET_HT = new TGraphErrors(3 * nBkgMETBins);
	TGraphErrors *error_HT_MET = new TGraphErrors(3 * nBkgHTBins);
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *error_nBJet = new TGraphErrors(5);
	TGraphErrors *ratioerror_nJet = new TGraphErrors(10);
	TGraphErrors *ratioerror_dPhiEleMET = new TGraphErrors(32);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_MET_HT = new TGraphErrors(3 * nBkgMETBins);
	TGraphErrors *ratioerror_HT_MET = new TGraphErrors(3 * nBkgHTBins);
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_nBJet = new TGraphErrors(5);

	TH1F *p_allnJet = (TH1F*)file_sig->Get("p_nJet");
	TH1F *p_allnBJet = (TH1F*)file_sig->Get("p_nBJet");
	TH1F *p_alldPhiEleMET = (TH1F*)file_sig->Get("p_dPhiEleMET");
	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMET_HT   = (TH1F*)file_sig->Get("h_MET_HT_concat");
	TH1F *p_allHT_MET   = (TH1F*)file_sig->Get("h_HT_MET_concat");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");
	TH1F *p_allPU    = (TH1F*)file_sig->Get("p_PU");

	p_allnJet->SetMarkerColor(1);
	p_allnBJet->SetMarkerColor(1);
	p_alldPhiEleMET->SetMarkerColor(1);
	p_allPhoEt->SetMarkerColor(1);
	p_allLepPt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allMET_HT->SetMarkerColor(1);
	p_allMt->SetMarkerColor(1);  
	p_allHT->SetMarkerColor(1);
	p_allPU->SetMarkerColor(1);
	
	TH1F *p_elenJet = (TH1F*)file_ele->Get("p_nJet");
	TH1F *p_elenBJet = (TH1F*)file_ele->Get("p_nBJet");
	TH1F *p_eledPhiEleMET = (TH1F*)file_ele->Get("p_dPhiEleMET");
	TH1F *p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt");
	TH1F *p_eleLepPt = (TH1F*)file_ele->Get("p_LepPt");
	TH1F *p_eleMET   = (TH1F*)file_ele->Get("p_MET");
	TH1F *p_eleMET_HT   = (TH1F*)file_ele->Get("h_MET_HT_concat");
	TH1F *p_eleHT_MET   = (TH1F*)file_ele->Get("h_HT_MET_concat");
	TH1F *p_eleMt    = (TH1F*)file_ele->Get("p_Mt");
	TH1F *p_eleHT  = (TH1F*)file_ele->Get("p_HT");
	TH1F *p_elePU    = (TH1F*)file_ele->Get("p_PU");

	TH1F *p_jetnJet = (TH1F*)file_jet->Get("p_nJet");
	TH1F *p_jetnBJet = (TH1F*)file_jet->Get("p_nBJet");
	TH1F *p_jetdPhiEleMET = (TH1F*)file_jet->Get("p_dPhiEleMET");
	TH1F *p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt");
	TH1F *p_jetLepPt = (TH1F*)file_jet->Get("p_LepPt");
	TH1F *p_jetMET   = (TH1F*)file_jet->Get("p_MET");
	TH1F *p_jetMET_HT   = (TH1F*)file_jet->Get("h_MET_HT_concat");
	TH1F *p_jetHT_MET   = (TH1F*)file_jet->Get("h_HT_MET_concat");
	TH1F *p_jetMt    = (TH1F*)file_jet->Get("p_Mt");
	TH1F *p_jetHT  = (TH1F*)file_jet->Get("p_HT");
	TH1F *p_jetPU    = (TH1F*)file_jet->Get("p_PU");

	TH1F *p_qcdnJet = (TH1F*)file_qcd->Get("p_nJet");
	TH1F *p_qcdnBJet = (TH1F*)file_qcd->Get("p_nBJet");
	TH1F *p_qcddPhiEleMET = (TH1F*)file_qcd->Get("p_dPhiEleMET");
	TH1F *p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt");
	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	TH1F *p_qcdMET   = (TH1F*)file_qcd->Get("p_MET");
	TH1F *p_qcdMET_HT   = (TH1F*)file_qcd->Get("h_MET_HT_concat");
	TH1F *p_qcdHT_MET   = (TH1F*)file_qcd->Get("h_HT_MET_concat");
	TH1F *p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt");
	TH1F *p_qcdHT  = (TH1F*)file_qcd->Get("p_HT");
	TH1F *p_qcdPU    = (TH1F*)file_qcd->Get("p_PU");

	TH1F *p_VGnJet = (TH1F*)file_VG->Get("p_nJet");
	TH1F *p_VGnBJet = (TH1F*)file_VG->Get("p_nBJet");
	TH1F *p_VGdPhiEleMET = (TH1F*)file_VG->Get("p_dPhiEleMET");
	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMET_HT   = (TH1F*)file_VG->Get("h_MET_HT_concat");
	TH1F *p_VGHT_MET   = (TH1F*)file_VG->Get("h_HT_MET_concat");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGHT  = (TH1F*)file_VG->Get("p_HT");
	TH1F *p_VGPU    = (TH1F*)file_VG->Get("p_PU");

	TH1F *p_rarenJet = (TH1F*)file_rare->Get("p_nJet");
	TH1F *p_rarenBJet = (TH1F*)file_rare->Get("p_nBJet");
	TH1F *p_raredPhiEleMET = (TH1F*)file_rare->Get("p_dPhiEleMET");
	TH1F *p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt");
	TH1F *p_rareLepPt = (TH1F*)file_rare->Get("p_LepPt");
	TH1F *p_rareMET   = (TH1F*)file_rare->Get("p_MET");
	TH1F *p_rareMET_HT   = (TH1F*)file_rare->Get("h_MET_HT_concat");
	TH1F *p_rareHT_MET   = (TH1F*)file_rare->Get("h_HT_MET_concat");
	TH1F *p_rareMt    = (TH1F*)file_rare->Get("p_Mt");
	TH1F *p_rareHT  = (TH1F*)file_rare->Get("p_HT");
	TH1F *p_rarePU    = (TH1F*)file_rare->Get("p_PU");

	if(doTT){
		p_allnJet = (TH1F*)file_sig->Get("p_nJet_TT");
		p_allnBJet = (TH1F*)file_sig->Get("p_nBJet_TT");
		p_alldPhiEleMET = (TH1F*)file_sig->Get("p_dPhiEleMET_TT");
		p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt_TT");
		p_allMET   = (TH1F*)file_sig->Get("p_MET_TT");
		p_allMt    = (TH1F*)file_sig->Get("p_Mt_TT");
		p_allHT  = (TH1F*)file_sig->Get("p_HT_TT");
		
		p_elenJet = (TH1F*)file_ele->Get("p_nJet_TT");
		p_elenBJet = (TH1F*)file_ele->Get("p_nBJet_TT");
		p_eledPhiEleMET = (TH1F*)file_ele->Get("p_dPhiEleMET_TT");
		p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt_TT");
		p_eleMET   = (TH1F*)file_ele->Get("p_MET_TT");
		p_eleMt    = (TH1F*)file_ele->Get("p_Mt_TT");
		p_eleHT  = (TH1F*)file_ele->Get("p_HT_TT");
	
		p_jetnJet = (TH1F*)file_jet->Get("p_nJet_TT");
		p_jetnBJet = (TH1F*)file_jet->Get("p_nBJet_TT");
		p_jetdPhiEleMET = (TH1F*)file_jet->Get("p_dPhiEleMET_TT");
		p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt_TT");
		p_jetMET   = (TH1F*)file_jet->Get("p_MET_TT");
		p_jetMt    = (TH1F*)file_jet->Get("p_Mt_TT");
		p_jetHT  = (TH1F*)file_jet->Get("p_HT_TT");
	
		p_qcdnJet = (TH1F*)file_qcd->Get("p_nJet_TT");
		p_qcdnBJet = (TH1F*)file_qcd->Get("p_nBJet_TT");
		p_qcddPhiEleMET = (TH1F*)file_qcd->Get("p_dPhiEleMET_TT");
		p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt_TT");
		p_qcdMET   = (TH1F*)file_qcd->Get("p_MET_TT");
		p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt_TT");
		p_qcdHT  = (TH1F*)file_qcd->Get("p_HT_TT");
	
		p_VGnJet = (TH1F*)file_VG->Get("p_nJet_TT");
		p_VGnBJet = (TH1F*)file_VG->Get("p_nBJet_TT");
		p_VGdPhiEleMET = (TH1F*)file_VG->Get("p_dPhiEleMET_TT");
		p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt_TT");
		p_VGMET   = (TH1F*)file_VG->Get("p_MET_TT");
		p_VGMt    = (TH1F*)file_VG->Get("p_Mt_TT");
		p_VGHT  = (TH1F*)file_VG->Get("p_HT_TT");
	
		p_rarenJet = (TH1F*)file_rare->Get("p_nJet_TT");
		p_rarenBJet = (TH1F*)file_rare->Get("p_nBJet_TT");
		p_raredPhiEleMET = (TH1F*)file_rare->Get("p_dPhiEleMET_TT");
		p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt_TT");
		p_rareMET   = (TH1F*)file_rare->Get("p_MET_TT");
		p_rareMt    = (TH1F*)file_rare->Get("p_Mt_TT");
		p_rareHT  = (TH1F*)file_rare->Get("p_HT_TT");
	}

	std::ostringstream nJetplot;  nJetplot.str("");
	std::ostringstream nBJetplot;  nBJetplot.str("");
	std::ostringstream dPhiplot;  dPhiplot.str("");
	std::ostringstream etplot;    etplot.str("");
	std::ostringstream ptplot;    ptplot.str("");
	std::ostringstream metplot;   metplot.str("");
	std::ostringstream methtplot;   methtplot.str("");
	std::ostringstream htmetplot;   htmetplot.str("");
	std::ostringstream mtplot;    mtplot.str("");
	std::ostringstream htplot;    htplot.str("");
	
	if(doTT){
	   if(channel == 1){
		if(plottype == 1){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_nJet_TT_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_nBJet_TT_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_met_TT_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_metht_TT_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_htmet_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_nJet_TT_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_nBJet_TT_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_met_TT_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_metht_TT_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_htmet_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
	   }
	   else if(channel == 2){
		if(plottype == 1){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_nJet_TT_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_nBJet_TT_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_met_TT_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_metht_TT_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_htmet_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_nJet_TT_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_nBJet_TT_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_met_TT_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_metht_TT_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_htmet_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
	   }
	}
	else{
	   if(channel == 1){
		if(plottype == 1){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_nJet_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_nBJet_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_met_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_metht_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_htmet_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_nJet_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_nBJet_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_met_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_metht_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_htmet_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_ht_" << RunYear<<whichVFP <<".png";	
		}
	   }
	   else if(channel == 2){
		if(plottype == 1){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_nJet_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_nBJet_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_met_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_metht_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_htmet_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			nJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_nJet_" << RunYear<<whichVFP <<".png";
			nBJetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_nBJet_" << RunYear<<whichVFP <<".png";
			dPhiplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_leppt_" << RunYear<<whichVFP <<".png";	
			metplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_met_" << RunYear<<whichVFP <<".png";
			methtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_metht_" << RunYear<<whichVFP <<".png";
			htmetplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_htmet_" << RunYear<<whichVFP <<".png";
			mtplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/uscms_data/d3/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_ht_" << RunYear<<whichVFP <<".png";	
		}
	   }
	}
	// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	setCanvas(c_mt); 
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(mt_pad1); 
	mt_pad1->SetBottomMargin(0);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_allMt->SetMinimum(0.5);
  	p_allMt->SetMaximum(10*p_allMt->GetBinContent(p_allMt->GetMaximumBin()));
	p_allMt->GetXaxis()->SetRangeUser(0,1000);
	p_allMt->SetLineColor(1);
	p_allMt->SetMarkerStyle(20);
	p_allMt->Draw("P");
	p_VGMt->SetFillStyle(1001);
	p_VGMt->SetLineColor(kMagenta);
	p_VGMt->SetFillColor(kMagenta);
	p_rareMt->SetFillStyle(1001);
	p_rareMt->SetLineColor(kYellow-4);
	p_rareMt->SetFillColor(kYellow-4);
	p_qcdMt->SetFillStyle(1001);
	p_qcdMt->SetLineColor(kBlue);
	p_qcdMt->SetFillColor(kBlue);
	p_eleMt->SetFillStyle(1001);
	p_eleMt->SetLineColor(kRed);
	p_eleMt->SetFillColor(kRed);
	p_jetMt->SetFillStyle(1001);
	p_jetMt->SetLineColor(kGreen);
	p_jetMt->SetFillColor(kGreen);
	
	float events_eleMt, events_jetMt, events_qcdMt, events_VGMt, events_rareMt, events_allMt;
	events_eleMt =  p_eleMt->Integral();
	events_jetMt =  p_jetMt->Integral();
	events_qcdMt =  p_qcdMt->Integral();
	events_VGMt =   p_VGMt->Integral();
	events_rareMt = p_rareMt->Integral();
	events_allMt =  p_allMt->Integral();
	cout<<"MT:\t" <<events_jetMt <<"\t"<<events_eleMt<<"\t"<<events_rareMt<<"\t"<<events_VGMt<<"\t"<<events_qcdMt<<"\t"<<events_allMt<<endl;		
        std::vector<StackComponent> mtComponents = {
		{p_rareMt, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
		{p_eleMt, "e#rightarrow#gamma fakes"},
		{p_jetMt, "j#rightarrow#gamma fakes"},
		{p_qcdMt, "fake lepton"},
		{p_VGMt, "W#gamma / Z#gamma"}
	};
	std::vector<StackComponent> mtStack = BuildIntegralOrderedStack(mtComponents, "mt_stack");
	TH1F *p_totalMt = mtStack.back().hist;
	for(int ibin(1); ibin < p_totalMt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,p_totalMt->GetBinCenter(ibin), p_totalMt->GetBinContent(ibin));
		float prederror = p_totalMt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(p_totalMt->GetBinLowEdge(ibin+1)-p_totalMt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,p_totalMt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(p_totalMt->GetBinLowEdge(ibin+1)-p_totalMt->GetBinLowEdge(ibin))/2, prederror/p_totalMt->GetBinContent(ibin)); 
        }
	
	for(int idx = mtStack.size()-1; idx >= 0; --idx) mtStack[idx].hist->Draw("hist same");
	// Put histograms in a vector of pairs {integral, hist pointer}

  	error_Mt->SetFillColor(kBlack);
  	error_Mt->SetFillStyle(3345);
	error_Mt->Draw("E2 same");
	TLegend *leg_mt =  new TLegend(0.4,0.6,0.9,0.9);
	leg_mt->SetTextSize(13);
	leg_mt->SetNColumns(2);
	leg_mt->SetFillStyle(0);
	leg_mt->SetBorderSize(0);
	leg_mt->SetFillColor(0);
	leg_mt->AddEntry(p_allMt,"observed");
	for(int idx = mtStack.size()-1; idx >= 0; --idx) leg_mt->AddEntry(mtStack[idx].hist, mtStack[idx].label);
	leg_mt->AddEntry(ratioerror_Mt, "Unc");
	leg_mt->Draw("same");
	p_allMt->Draw("pe same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( mt_pad1,1,channel, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( mt_pad1,2,channel, 11 );
	else if(RunYear==2017)  		CMS_lumi( mt_pad1,3,channel, 11 );
	else if(RunYear==2018)  		CMS_lumi( mt_pad1,4,channel, 11 );

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.3);
	mt_pad2->SetTopMargin(0);
	mt_pad2->SetBottomMargin(0.4);
	mt_pad2->Draw();
	mt_pad2->cd();
  	TLine *flatratio_mt = new TLine(0,1,1000,1);
	TH1F *ratio_mt=(TH1F*)p_allMt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->Divide(p_totalMt);
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("obs./bkg.");
	ratio_mt->GetXaxis()->SetTitleOffset(0.9);
	ratio_mt->GetXaxis()->SetTitle("M_{T} (GeV)");

	ratio_mt->GetYaxis()->SetRangeUser(0.6,1.4);
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(11);
	ratio_mt->Draw("pe");
	ratioerror_Mt->SetFillColor(kBlack);
	ratioerror_Mt->SetFillStyle(3345);
	ratioerror_Mt->Draw("E2 same");
	flatratio_mt->Draw("same");

  	TLine *ratioValue_mt = new TLine(0,p_allMt->Integral()/p_VGMt->Integral(),1000,p_allMt->Integral()/p_VGMt->Integral());
	ratioValue_mt->SetLineColor(kRed);
	//ratioValue_mt->Draw("same");
	c_mt->SaveAs(mtplot.str().c_str());
 
	// ******** PhoEt ************************//
  	gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon p_{T}",600,600);
	setCanvas(c_pt); 
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(pt_pad1); 
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetMaximum(10*p_allPhoEt->GetBinContent(p_allPhoEt->GetMaximumBin()));
	p_allPhoEt->SetMinimum(0.5);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,800);
	p_allPhoEt->SetLineColor(1);
	p_allPhoEt->SetMarkerStyle(20);
	p_allPhoEt->Draw("P");
	p_VGPhoEt->SetFillStyle(1001);
	p_VGPhoEt->SetLineColor(kMagenta);
	p_VGPhoEt->SetFillColor(kMagenta);
	p_rarePhoEt->SetFillStyle(1001);
	p_rarePhoEt->SetLineColor(kYellow-4);
	p_rarePhoEt->SetFillColor(kYellow-4);
	p_qcdPhoEt->SetFillStyle(1001);
	p_qcdPhoEt->SetLineColor(kBlue);
	p_qcdPhoEt->SetFillColor(kBlue);
	p_elePhoEt->SetFillStyle(1001);
	p_elePhoEt->SetLineColor(kRed);
	p_elePhoEt->SetFillColor(kRed);
	p_jetPhoEt->SetFillStyle(1001);
	p_jetPhoEt->SetLineColor(kGreen);
	p_jetPhoEt->SetFillColor(kGreen);
	
	float events_elePhoEt, events_jetPhoEt, events_qcdPhoEt, events_VGPhoEt, events_rarePhoEt, events_allPhoEt;
	events_elePhoEt =  p_elePhoEt->Integral();
	events_jetPhoEt =  p_jetPhoEt->Integral();
	events_qcdPhoEt =  p_qcdPhoEt->Integral();
	events_VGPhoEt =   p_VGPhoEt->Integral();
	events_rarePhoEt = p_rarePhoEt->Integral();
	events_allPhoEt =  p_allPhoEt->Integral();
	cout<<"Pho pT:\t" <<events_elePhoEt <<"\t"<<events_jetPhoEt<<"\t"<<events_qcdPhoEt<<"\t"<<events_VGPhoEt<<"\t"<<events_rarePhoEt<<"\t"<<events_allPhoEt<<endl;		
	std::vector<StackComponent> phoEtComponents = {
		{p_rarePhoEt, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
		{p_elePhoEt, "e#rightarrow#gamma fakes"},
		{p_jetPhoEt, "j#rightarrow#gamma fakes"},
		{p_qcdPhoEt, "fake lepton"},
		{p_VGPhoEt, "W#gamma / Z#gamma"}
	};
	std::vector<StackComponent> phoEtStack = BuildIntegralOrderedStack(phoEtComponents, "phoet_stack");
	TH1F *p_totalPhoEt = phoEtStack.back().hist;
	p_totalPhoEt->Sumw2();
	for(int ibin(1); ibin < p_totalPhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,p_totalPhoEt->GetBinCenter(ibin), p_totalPhoEt->GetBinContent(ibin));
		float prederror = p_totalPhoEt->GetBinError(ibin);
		error_PhoEt->SetPointError(ibin-1,(p_totalPhoEt->GetBinLowEdge(ibin+1)-p_totalPhoEt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_PhoEt->SetPoint(ibin-1,p_totalPhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(p_totalPhoEt->GetBinLowEdge(ibin+1)-p_totalPhoEt->GetBinLowEdge(ibin))/2, prederror/p_totalPhoEt->GetBinContent(ibin)); 
	}	
	for(int idx = phoEtStack.size()-1; idx >= 0; --idx) phoEtStack[idx].hist->Draw("hist same");
  	error_PhoEt->SetFillColor(kBlack);
  	error_PhoEt->SetFillStyle(3345);
	error_PhoEt->Draw("E2 same");
	leg_mt->Draw("same");
	p_allPhoEt->Draw("PE same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( pt_pad1,1,channel, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( pt_pad1,2,channel, 11 );
	else if(RunYear==2017)  		CMS_lumi( pt_pad1,3,channel, 11 );
	else if(RunYear==2018)  		CMS_lumi( pt_pad1,4,channel, 11 );

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.3);
	pt_pad2->SetTopMargin(0);
	pt_pad2->SetBottomMargin(0.4);
	pt_pad2->Draw();
	pt_pad2->cd();
  	TLine *flatratio = new TLine(35,1,800,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	        ratio->GetYaxis()->SetRangeUser(0.6,1.4);

	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_totalPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("obs./bkg.");
	ratio->GetXaxis()->SetTitleOffset(0.9);
	ratio->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(11);
	ratio->Draw("ep");
	ratioerror_PhoEt->SetFillColor(kBlack);
	ratioerror_PhoEt->SetFillStyle(3345);
	ratioerror_PhoEt->Draw("E2 same");
	flatratio->Draw("same");
  	TLine *ratioValue = new TLine(35,p_allPhoEt->Integral()/p_VGPhoEt->Integral(),800,p_allPhoEt->Integral()/p_VGPhoEt->Integral());
	ratioValue->SetLineColor(kRed);
	//ratioValue->Draw("same");
	c_pt->SaveAs(etplot.str().c_str());


	// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	setCanvas(c_met); 
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->SetBottomMargin(0);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->GetYaxis()->SetRangeUser(1, 10*p_allMET->GetBinContent(p_allMET->GetMaximumBin()));
	p_allMET->SetMinimum(0.5);
	p_allMET->GetXaxis()->SetRangeUser(0,1000);
	p_allMET->SetLineColor(1);
	p_allMET->SetMarkerStyle(20);
	p_allMET->Draw("P");
	p_VGMET->SetFillStyle(1001);
	p_VGMET->SetLineColor(kMagenta);
	p_VGMET->SetFillColor(kMagenta);
	p_rareMET->SetFillStyle(1001);
	p_rareMET->SetLineColor(kYellow-4);
	p_rareMET->SetFillColor(kYellow-4);
	p_qcdMET->SetFillStyle(1001);
	p_qcdMET->SetLineColor(kBlue);
	p_qcdMET->SetFillColor(kBlue);
	p_eleMET->SetFillStyle(1001);
	p_eleMET->SetLineColor(kRed);
	p_eleMET->SetFillColor(kRed);
	p_jetMET->SetFillStyle(1001);
	p_jetMET->SetLineColor(kGreen);
	p_jetMET->SetFillColor(kGreen);
	float events_eleMET, events_jetMET, events_qcdMET, events_VGMET, events_rareMET, events_allMET;
	events_eleMET =  p_eleMET->Integral();
	events_jetMET =  p_jetMET->Integral();
	events_qcdMET =  p_qcdMET->Integral();
	events_VGMET =   p_VGMET->Integral();
	events_rareMET = p_rareMET->Integral();
	events_allMET =  p_allMET->Integral();
	cout<<"MET:\t" <<events_eleMET <<"\t"<<events_jetMET<<"\t"<<events_qcdMET<<"\t"<<events_VGMET<<"\t"<<events_rareMET<<"\t"<<events_allMET<<endl;		

	std::vector<StackComponent> metComponents = {
		{p_rareMET, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
		{p_eleMET, "e#rightarrow#gamma fakes"},
		{p_jetMET, "j#rightarrow#gamma fakes"},
		{p_qcdMET, "fake lepton"},
		{p_VGMET, "W#gamma / Z#gamma"}
	};
	std::vector<StackComponent> metStack = BuildIntegralOrderedStack(metComponents, "met_stack");
	TH1F *p_totalMET = metStack.back().hist;
	p_totalMET->Sumw2();
	for(int ibin(1); ibin < p_totalMET->GetSize(); ibin++){
		float prederror = p_totalMET->GetBinError(ibin);
		error_MET->SetPoint(ibin-1,p_totalMET->GetBinCenter(ibin), p_totalMET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(p_totalMET->GetBinLowEdge(ibin+1)-p_totalMET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,p_totalMET->GetBinCenter(ibin), 1);
		ratioerror_MET->SetPointError(ibin-1,(p_totalMET->GetBinLowEdge(ibin+1)-p_totalMET->GetBinLowEdge(ibin))/2, prederror/p_totalMET->GetBinContent(ibin));
	}
	for(int idx = metStack.size()-1; idx >= 0; --idx) metStack[idx].hist->Draw("hist same");	
		
	error_MET->SetFillColor(kBlack);
		error_MET->SetFillStyle(3345);
		error_MET->Draw("E2 same");
		leg_mt->Draw("same");
		p_allMET->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( met_pad1, 1,channel, 11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( met_pad1, 2,channel, 11 );
		else if(RunYear==2017)  		CMS_lumi( met_pad1, 3,channel, 11 );
		else if(RunYear==2018)  		CMS_lumi( met_pad1, 4,channel, 11 );

		c_met->cd();
		TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
		met_pad2->SetTopMargin(0);
		met_pad2->SetBottomMargin(0.4);
		met_pad2->Draw();
		met_pad2->cd();
		TLine *flatratio_met = new TLine(0,1,1000,1);
		TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
		ratio_met->SetLineColor(kBlack);
		ratio_met->SetMarkerStyle(20);
		ratio_met->Divide(p_totalMET);
		ratio_met->SetTitle("");
		ratio_met->GetYaxis()->SetTitle("obs./bkg.");
		ratio_met->GetXaxis()->SetTitleOffset(0.9);
		ratio_met->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
		ratio_met->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_met->GetXaxis()->SetLabelFont(63);
		ratio_met->GetXaxis()->SetLabelSize(14);
		ratio_met->GetYaxis()->SetLabelFont(63);
		ratio_met->GetYaxis()->SetLabelSize(11);
		ratio_met->Draw("pe");
		ratioerror_MET->SetFillColor(kBlack);
		ratioerror_MET->SetFillStyle(3345);
		ratioerror_MET->Draw("E2 same");
		flatratio_met->Draw("same");

		TLine *ratioValue_met = new TLine(0,p_allMET->Integral()/p_VGMET->Integral(),1000,p_allMET->Integral()/p_VGMET->Integral());
		ratioValue_met->SetLineColor(kRed);
		//ratioValue_met->Draw("same");

		c_met->SaveAs(metplot.str().c_str());

	// ******** HT ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
		setCanvas(c_HT); 
		c_HT->cd();
		TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.3, 1, 1.0);
		setTopPad(HT_pad1); 
		HT_pad1->SetBottomMargin(0);
		HT_pad1->Draw();  
		HT_pad1->cd();  
		gPad->SetLogy();
		p_allHT->GetXaxis()->SetRangeUser(0,900);
		p_allHT->SetMinimum(5);
		p_allHT->SetMaximum(10*p_allHT->GetBinContent(p_allHT->GetMaximumBin()));
		p_allHT->SetLineColor(1);
		p_allHT->SetMarkerStyle(20);
		p_allHT->Draw("P");
		p_VGHT->SetFillStyle(1001);
		p_VGHT->SetLineColor(kMagenta);
		p_VGHT->SetFillColor(kMagenta);
		p_rareHT->SetFillStyle(1001);
		p_rareHT->SetLineColor(kYellow-4);
		p_rareHT->SetFillColor(kYellow-4);
		p_qcdHT->SetFillStyle(1001);
		p_qcdHT->SetLineColor(kBlue);
		p_qcdHT->SetFillColor(kBlue);
		p_eleHT->SetFillStyle(1001);
		p_eleHT->SetLineColor(kRed);
		p_eleHT->SetFillColor(kRed);
		p_jetHT->SetFillStyle(1001);
		p_jetHT->SetLineColor(kGreen);
		p_jetHT->SetFillColor(kGreen);
		
                for(int ibin(1); ibin < p_VGHT->GetSize(); ibin++){
                        cout << ibin << "\tData\t"<< p_allHT->GetBinContent(ibin) << "\teleFake\t" <<  p_eleHT->GetBinContent(ibin) << "\t jetFake\t" << p_jetHT->GetBinContent(ibin) << "\t qcd Fake \t" << p_qcdHT->GetBinContent(ibin) << "\t VG \t" << p_VGHT->GetBinContent(ibin) << "\t rare \t" << p_rareHT->GetBinContent(ibin) << "\n";

                        cout << ibin << "\tData\t"<< p_allHT->GetBinError(ibin) << "\teleFake\t" <<  p_eleHT->GetBinError(ibin) << "\t jetFake\t" << p_jetHT->GetBinError(ibin) << "\t qcd Fake \t" << p_qcdHT->GetBinError(ibin) << "\t VG \t" << p_VGHT->GetBinError(ibin) << "\t rare \t" << p_rareHT->GetBinError(ibin) << "\n";
                        cout << ibin << "\tData\t"<< p_allHT->GetBinError(ibin)/p_allHT->GetBinContent(ibin) << "\teleFake\t" <<  p_eleHT->GetBinError(ibin)/p_eleHT->GetBinContent(ibin) << "\t jetFake\t" << p_jetHT->GetBinError(ibin)/p_jetHT->GetBinContent(ibin) << "\t qcd Fake \t" << p_qcdHT->GetBinError(ibin)/p_qcdHT->GetBinContent(ibin) << "\t VG \t" << p_VGHT->GetBinError(ibin)/p_VGHT->GetBinContent(ibin) << "\t rare \t" << p_rareHT->GetBinError(ibin)/p_rareHT->GetBinContent(ibin) << "\n";

                        cout <<"\n";
                }

		float events_eleHT, events_jetHT, events_qcdHT, events_VGHT, events_rareHT, events_allBkg, events_allHT;
		events_eleHT =  p_eleHT->Integral();
		events_jetHT =  p_jetHT->Integral();
		events_qcdHT =  p_qcdHT->Integral();
		events_VGHT =   p_VGHT->Integral();
		events_rareHT = p_rareHT->Integral();
		events_allHT =  p_allHT->Integral();
		events_allBkg = p_eleHT->Integral() + p_jetHT->Integral() + p_qcdHT->Integral() + p_VGHT->Integral() + p_rareHT->Integral();

		if(channel == 1) cout<<endl<<endl<<"eg channel :                 "<<endl;
		if(channel == 2) cout<<endl<<endl<<"mg channel :                 "<<endl;

        	cout<<"elefakePho :   "<< events_eleHT << endl;
        	cout<<"jetfakePho :   "<< events_jetHT << endl;
        	cout<<"jetfakeLep :   "<< events_qcdHT << endl;
        	cout<<"VGamma :       "<< events_VGHT << endl;
        	cout<<"Electroweak :         "<< events_rareHT <<endl;
        	cout<<"total SM backgrounds:    "<< events_allBkg  <<"\tTotal data  "<< events_allHT <<endl;
        	cout<<"Discrepancy :  "<<100*(events_allBkg-events_allHT) / events_allHT <<" % "<<endl;

		std::vector<StackComponent> htComponents = {
			{p_rareHT, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
			{p_eleHT, "e#rightarrow#gamma fakes"},
			{p_jetHT, "j#rightarrow#gamma fakes"},
			{p_qcdHT, "fake lepton"},
			{p_VGHT, "W#gamma / Z#gamma"}
		};
		std::vector<StackComponent> htStack = BuildIntegralOrderedStack(htComponents, "ht_stack");
		TH1F *p_totalHT = htStack.back().hist;
		for(int ibin(1); ibin < p_totalHT->GetSize(); ibin++){
			error_HT->SetPoint(ibin-1,p_totalHT->GetBinCenter(ibin), p_totalHT->GetBinContent(ibin));
			float prederror = p_totalHT->GetBinError(ibin);
			error_HT->SetPointError(ibin-1,(p_totalHT->GetBinLowEdge(ibin+1)-p_totalHT->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_HT->SetPoint(ibin-1,p_totalHT->GetBinCenter(ibin), 1); 
			ratioerror_HT->SetPointError(ibin-1,(p_totalHT->GetBinLowEdge(ibin+1)-p_totalHT->GetBinLowEdge(ibin))/2, prederror/p_totalHT->GetBinContent(ibin)); 
		}
		for(int idx = htStack.size()-1; idx >= 0; --idx) htStack[idx].hist->Draw("hist same");
		error_HT->SetFillColor(kBlack);
		error_HT->SetFillStyle(3345);
		error_HT->Draw("E2 same");

		leg_mt->Draw("same");
		p_allHT->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( HT_pad1,1, channel,11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( HT_pad1,2, channel,11 );
		else if(RunYear==2017)  		CMS_lumi( HT_pad1,3, channel,11 );
		else if(RunYear==2018)  		CMS_lumi( HT_pad1,4, channel,11 );

		c_HT->cd();
		TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.3);
		HT_pad2->SetTopMargin(0);
		HT_pad2->SetBottomMargin(0.4);
		HT_pad2->Draw();
		HT_pad2->cd();
		TLine *flatratio_HT = new TLine(0,1,1000,1);
		TH1F *ratio_HT=(TH1F*)p_allHT->Clone("transfer factor");
		ratio_HT->SetMarkerStyle(20);
		ratio_HT->SetLineColor(kBlack);
		ratio_HT->GetXaxis()->SetRangeUser(0,900);
		ratio_HT->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_HT->Divide(p_totalHT);
		ratio_HT->SetTitle("");
		ratio_HT->GetYaxis()->SetTitle("obs./bkg.");
		ratio_HT->GetXaxis()->SetTitleOffset(0.9);
		ratio_HT->GetXaxis()->SetTitle("H_{T} (GeV)");
		ratio_HT->GetXaxis()->SetLabelFont(63);
		ratio_HT->GetXaxis()->SetLabelSize(14);
		ratio_HT->GetYaxis()->SetLabelFont(63);
		ratio_HT->GetYaxis()->SetLabelSize(11);
		ratio_HT->Draw("ep");
		ratioerror_HT->SetFillColor(kBlack);
		ratioerror_HT->SetFillStyle(3345);
		ratioerror_HT->Draw("E2 same");
		flatratio_HT->Draw("same");

		TLine *ratioValue_ht = new TLine(0,p_allHT->Integral()/p_VGHT->Integral(),1000,p_allHT->Integral()/p_VGHT->Integral());
		ratioValue_ht->SetLineColor(kRed);
		//ratioValue_ht->Draw("same");
		c_HT->SaveAs(htplot.str().c_str());

		// ******** MET_HT ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_metht = new TCanvas("MET_HT", "MET_HT",600,600);
		setCanvas(c_metht);
		c_metht->cd();
		TPad *metht_pad1 = new TPad("metht_pad1", "metht_pad1", 0, 0.3, 1, 1.0);
		setTopPad(metht_pad1);
		metht_pad1->SetBottomMargin(0);
		metht_pad1->Draw();
		metht_pad1->cd();
		gPad->SetLogy();
		p_allMET_HT->GetYaxis()->SetRangeUser(1, 50*p_allMET_HT->GetBinContent(p_allMET_HT->GetMaximumBin()));
		p_allMET_HT->SetMinimum(0.5);
		p_allMET_HT->GetXaxis()->SetRangeUser(0,1000);
		p_allMET_HT->SetLineColor(1);
		p_allMET_HT->SetMarkerStyle(20);
		p_allMET_HT->Draw("P");
		p_VGMET_HT->SetFillStyle(1001);
		p_VGMET_HT->SetLineColor(kMagenta);
		p_VGMET_HT->SetFillColor(kMagenta);
		p_rareMET_HT->SetFillStyle(1001);
		p_rareMET_HT->SetLineColor(kYellow-4);
		p_rareMET_HT->SetFillColor(kYellow-4);
		p_qcdMET_HT->SetFillStyle(1001);
		p_qcdMET_HT->SetLineColor(kBlue);
		p_qcdMET_HT->SetFillColor(kBlue);
		p_eleMET_HT->SetFillStyle(1001);
		p_eleMET_HT->SetLineColor(kRed);
		p_eleMET_HT->SetFillColor(kRed);
		p_jetMET_HT->SetFillStyle(1001);
		p_jetMET_HT->SetLineColor(kGreen);
		p_jetMET_HT->SetFillColor(kGreen);
		float events_eleMET_HT, events_jetMET_HT, events_qcdMET_HT, events_VGMET_HT, events_rareMET_HT, events_allMET_HT;
		events_eleMET_HT =  p_eleMET_HT->Integral();
		events_jetMET_HT =  p_jetMET_HT->Integral();
		events_qcdMET_HT =  p_qcdMET_HT->Integral();
		events_VGMET_HT =   p_VGMET_HT->Integral();
		events_rareMET_HT = p_rareMET_HT->Integral();
		events_allMET_HT =  p_allMET_HT->Integral();
		cout<< "MET_HT:\t" <<events_eleMET_HT <<"\t"<<events_jetMET_HT<<"\t"<<events_qcdMET_HT<<"\t"<<events_VGMET_HT<<"\t"<<events_rareMET_HT<<"\t"<<events_allMET_HT<<endl;

		std::vector<StackComponent> methtComponents = {
			{p_rareMET_HT, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
			{p_eleMET_HT, "e#rightarrow#gamma fakes"},
			{p_jetMET_HT, "j#rightarrow#gamma fakes"},
			{p_qcdMET_HT, "fake lepton"},
			{p_VGMET_HT, "W#gamma / Z#gamma"}
		};
		std::vector<StackComponent> methtStack = BuildIntegralOrderedStack(methtComponents, "metht_stack");
		TH1F *p_totalMET_HT = methtStack.back().hist;
		p_totalMET_HT->Sumw2();
		for(int ibin(1); ibin < p_totalMET_HT->GetSize(); ibin++){
			float prederror = p_totalMET_HT->GetBinError(ibin);
			error_MET_HT->SetPoint(ibin-1,p_totalMET_HT->GetBinCenter(ibin), p_totalMET_HT->GetBinContent(ibin));
			error_MET_HT->SetPointError(ibin-1,(p_totalMET_HT->GetBinLowEdge(ibin+1)-p_totalMET_HT->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_MET_HT->SetPoint(ibin-1,p_totalMET_HT->GetBinCenter(ibin), 1); 
			ratioerror_MET_HT->SetPointError(ibin-1,(p_totalMET_HT->GetBinLowEdge(ibin+1)-p_totalMET_HT->GetBinLowEdge(ibin))/2, prederror/p_totalMET_HT->GetBinContent(ibin)); 
		}
		for(int idx = methtStack.size()-1; idx >= 0; --idx) methtStack[idx].hist->Draw("hist same");
		error_MET_HT->SetFillColor(kBlack);
		error_MET_HT->SetFillStyle(3345);
		error_MET_HT->Draw("E2 same");
		leg_mt->Draw("same");
		p_allMET_HT->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( metht_pad1, 1,channel, 11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( metht_pad1, 2,channel, 11 );
		else if(RunYear==2017)  		CMS_lumi( metht_pad1, 3,channel, 11 );
		else if(RunYear==2018)  		CMS_lumi( metht_pad1, 4,channel, 11 );

		TLatex latex;
		latex.SetNDC();
		latex.SetTextFont(42);
		latex.SetTextSize(0.035);
		latex.SetTextAlign(22);
		latex.SetTextColor(kGray+2);

		if(channel == 1) 
			latex.DrawLatex(0.30, 0.6, "e#gamma, H_{T} < 100 GeV");
		else 
			latex.DrawLatex(0.30, 0.6, "#mu#gamma, H_{T} < 100 GeV");
		latex.DrawLatex(0.55, 0.6, "100 < H_{T} < 400 GeV");
		latex.DrawLatex(0.80, 0.6, "H_{T} > 400 GeV");
		int nBins = p_allMET_HT->GetNbinsX();
		TLine line1(nBins/3, 0, nBins/3, p_allMET_HT->GetMaximum()*1.2);
		TLine line2(2*nBins/3, 0, 2*nBins/3, p_allMET_HT->GetMaximum()*1.2);
		line1.SetLineStyle(2);
		line2.SetLineStyle(2);
		line1.Draw("same");
		line2.Draw("same");
		
		c_metht->cd();
		TPad *metht_pad2 = new TPad("metht_pad2", "metht_pad2", 0, 0, 1, 0.3);
		metht_pad2->SetTopMargin(0);
		metht_pad2->SetBottomMargin(0.4);
		metht_pad2->Draw();
		metht_pad2->cd();
		TLine *flatratio_metht = new TLine(0,1,1000,1);
		TH1F *ratio_metht=(TH1F*)p_allMET_HT->Clone("transfer factor");
		//ratio_metht->GetXaxis()->SetRangeUser(120,1000);
		ratio_metht->SetLineColor(kBlack);
		ratio_metht->SetMarkerStyle(20);
		ratio_metht->Divide(p_totalMET_HT);
		//cout<< p_allMET_HT->Integral()/p_VGMET_HT->Integral()<<endl;
		ratio_metht->SetTitle("");
		ratio_metht->GetYaxis()->SetTitle("obs./bkg.");
		ratio_metht->GetXaxis()->SetTitleOffset(0.9);
		ratio_metht->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
		ratio_metht->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_metht->GetXaxis()->SetLabelFont(63);
		ratio_metht->GetXaxis()->SetLabelSize(14);
		ratio_metht->GetYaxis()->SetLabelFont(63);
		ratio_metht->GetYaxis()->SetLabelSize(11);
		ratio_metht->Draw("pe");
		ratioerror_MET_HT->SetFillColor(kBlack);
		ratioerror_MET_HT->SetFillStyle(3345);
		ratioerror_MET_HT->Draw("E2 same");
		flatratio_metht->Draw("same");

		TLine *ratioValue_metht = new TLine(0,p_allMET_HT->Integral()/p_VGMET_HT->Integral(),1000,p_allMET_HT->Integral()/p_VGMET_HT->Integral());
		ratioValue_metht->SetLineColor(kRed);
		//ratioValue_metht->Draw("same");
		c_metht->SaveAs(methtplot.str().c_str());

		// ******** HT_MET ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_htmet = new TCanvas("HT_MET", "HT_MET",600,600);
		setCanvas(c_htmet);
		c_htmet->cd();
		TPad *htmet_pad1 = new TPad("htmet_pad1", "htmet_pad1", 0, 0.3, 1, 1.0);
		setTopPad(htmet_pad1);
		htmet_pad1->SetBottomMargin(0);
		htmet_pad1->Draw();
		htmet_pad1->cd();
		gPad->SetLogy();
		p_allHT_MET->GetYaxis()->SetRangeUser(1, 50*p_allHT_MET->GetBinContent(p_allHT_MET->GetMaximumBin()));
		p_allHT_MET->SetMinimum(0.5);
		p_allHT_MET->GetXaxis()->SetRangeUser(0,1000);
		p_allHT_MET->SetLineColor(1);
		p_allHT_MET->SetMarkerStyle(20);
		p_allHT_MET->Draw("P");
		p_VGHT_MET->SetFillStyle(1001);
		p_VGHT_MET->SetLineColor(kMagenta);
		p_VGHT_MET->SetFillColor(kMagenta);
		p_rareHT_MET->SetFillStyle(1001);
		p_rareHT_MET->SetLineColor(kYellow-4);
		p_rareHT_MET->SetFillColor(kYellow-4);
		p_qcdHT_MET->SetFillStyle(1001);
		p_qcdHT_MET->SetLineColor(kBlue);
		p_qcdHT_MET->SetFillColor(kBlue);
		p_eleHT_MET->SetFillStyle(1001);
		p_eleHT_MET->SetLineColor(kRed);
		p_eleHT_MET->SetFillColor(kRed);
		p_jetHT_MET->SetFillStyle(1001);
		p_jetHT_MET->SetLineColor(kGreen);
		p_jetHT_MET->SetFillColor(kGreen);
		float events_eleHT_MET, events_jetHT_MET, events_qcdHT_MET, events_VGHT_MET, events_rareHT_MET, events_allHT_MET;
		events_eleHT_MET =  p_eleHT_MET->Integral();
		events_jetHT_MET =  p_jetHT_MET->Integral();
		events_qcdHT_MET =  p_qcdHT_MET->Integral();
		events_VGHT_MET =   p_VGHT_MET->Integral();
		events_rareHT_MET = p_rareHT_MET->Integral();
		events_allHT_MET =  p_allHT_MET->Integral();
		cout<< "HT_MET:\t" <<events_eleHT_MET <<"\t"<<events_jetHT_MET<<"\t"<<events_qcdHT_MET<<"\t"<<events_VGHT_MET<<"\t"<<events_rareHT_MET<<"\t"<<events_allHT_MET<<endl;

		std::vector<StackComponent> htmetComponents = {
			{p_rareHT_MET, "t#bar{t}#gamma / WW#gamma / WZ#gamma"},
			{p_eleHT_MET, "e#rightarrow#gamma fakes"},
			{p_jetHT_MET, "j#rightarrow#gamma fakes"},
			{p_qcdHT_MET, "fake lepton"},
			{p_VGHT_MET, "W#gamma / Z#gamma"}
		};
		std::vector<StackComponent> htmetStack = BuildIntegralOrderedStack(htmetComponents, "metht_stack");
		TH1F *p_totalHT_MET = htmetStack.back().hist;
		for(int ibin(1); ibin < p_totalHT_MET->GetSize(); ibin++){
			error_HT_MET->SetPoint(ibin-1,p_totalHT_MET->GetBinCenter(ibin), p_totalHT_MET->GetBinContent(ibin));
			float prederror = p_totalHT_MET->GetBinError(ibin);
			error_HT_MET->SetPointError(ibin-1,(p_totalHT_MET->GetBinLowEdge(ibin+1)-p_totalHT_MET->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_HT_MET->SetPoint(ibin-1,p_totalHT_MET->GetBinCenter(ibin), 1); 
			ratioerror_HT_MET->SetPointError(ibin-1,(p_totalHT_MET->GetBinLowEdge(ibin+1)-p_totalHT_MET->GetBinLowEdge(ibin))/2, prederror/p_totalHT_MET->GetBinContent(ibin)); 
		}
		for(int idx = htmetStack.size()-1; idx >= 0; --idx) htmetStack[idx].hist->Draw("hist same");
		error_HT_MET->SetFillColor(kBlack);
		error_HT_MET->SetFillStyle(3345);
		error_HT_MET->Draw("E2 same");
		leg_mt->Draw("same");
		p_allHT_MET->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( htmet_pad1, 1,channel, 11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( htmet_pad1, 2,channel, 11 );
		else if(RunYear==2017)  		CMS_lumi( htmet_pad1, 3,channel, 11 );
		else if(RunYear==2018)  		CMS_lumi( htmet_pad1, 4,channel, 11 );

		TLatex latexx;
		latexx.SetNDC();
		latexx.SetTextFont(42);
		latexx.SetTextSize(0.035);
		latexx.SetTextAlign(22);
		latexx.SetTextColor(kGray+2);

		if (channel == 1)
    			latexx.DrawLatex(0.30, 0.6, "e#gamma, 120 < p_{T}^{miss} < 200 GeV");
		else
    			latexx.DrawLatex(0.30, 0.6, "#mu#gamma, 120 < p_{T}^{miss} < 200 GeV");

		latexx.DrawLatex(0.55, 0.6, "200 < p_{T}^{miss} < 400 GeV");
		latexx.DrawLatex(0.80, 0.6, "p_{T}^{miss} > 400 GeV");
		nBins = p_allHT_MET->GetNbinsX();
		TLine line3(nBins/3, 0, nBins/3, p_allHT_MET->GetMaximum()*1.2);
		TLine line4(2*nBins/3, 0, 2*nBins/3, p_allHT_MET->GetMaximum()*1.2);
		line3.SetLineStyle(2);
		line4.SetLineStyle(2);
		line3.Draw("same");
		line4.Draw("same");
		
		c_htmet->cd();
		TPad *htmet_pad2 = new TPad("htmet_pad2", "htmet_pad2", 0, 0, 1, 0.3);
		htmet_pad2->SetTopMargin(0);
		htmet_pad2->SetBottomMargin(0.4);
		htmet_pad2->Draw();
		htmet_pad2->cd();
		TLine *flatratio_htmet = new TLine(0,1,1000,1);
		TH1F *ratio_htmet=(TH1F*)p_allHT_MET->Clone("transfer factor");
		//ratio_htmet->GetXaxis()->SetRangeUser(120,1000);
		ratio_htmet->SetLineColor(kBlack);
		ratio_htmet->SetMarkerStyle(20);
		ratio_htmet->Divide(p_totalHT_MET);
		//cout<< p_allHT_MET->Integral()/p_VGHT_MET->Integral()<<endl;
		ratio_htmet->SetTitle("");
		ratio_htmet->GetYaxis()->SetTitle("obs./bkg.");
		ratio_htmet->GetXaxis()->SetTitleOffset(0.9);
		ratio_htmet->GetXaxis()->SetTitle("H_{T} (GeV)");
		ratio_htmet->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_htmet->GetXaxis()->SetLabelFont(63);
		ratio_htmet->GetXaxis()->SetLabelSize(14);
		ratio_htmet->GetYaxis()->SetLabelFont(63);
		ratio_htmet->GetYaxis()->SetLabelSize(11);
		ratio_htmet->Draw("pe");
		ratioerror_HT_MET->SetFillColor(kBlack);
		ratioerror_HT_MET->SetFillStyle(3345);
		ratioerror_HT_MET->Draw("E2 same");
		flatratio_htmet->Draw("same");

		TLine *ratioValue_htmet = new TLine(0,p_allHT_MET->Integral()/p_VGHT_MET->Integral(),1000,p_allHT_MET->Integral()/p_VGHT_MET->Integral());
		ratioValue_htmet->SetLineColor(kRed);
		//ratioValue_htmet->Draw("same");
		c_htmet->SaveAs(htmetplot.str().c_str());

		// ******** LepPt ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
		setCanvas(c_leppt);
		c_leppt->cd();
		TPad *leppt_pad1 = new TPad("leppt_pad1", "leppt_pad1", 0, 0.3, 1, 1.0);
		setTopPad(leppt_pad1);
		leppt_pad1->SetBottomMargin(0);
		leppt_pad1->Draw();
		leppt_pad1->cd();
		gPad->SetLogy();
		p_allLepPt->SetMinimum(0.5);
		p_allLepPt->SetMaximum(10*p_allLepPt->GetBinContent(p_allLepPt->GetMaximumBin()));
	//        p_allLepPt->SetMinimum(p_rareLepPt->GetBinContent(p_rareLepPt->GetSize()));
	//        p_allLepPt->GetXaxis()->SetRangeUser(35,800);
		p_allLepPt->SetLineColor(1);
		p_allLepPt->SetMarkerStyle(20);
		p_allLepPt->Draw("P");
		p_VGLepPt->SetFillStyle(1001);
		p_VGLepPt->SetLineColor(kMagenta);
		p_VGLepPt->SetFillColor(kMagenta);
		p_rareLepPt->SetFillStyle(1001);
		p_rareLepPt->SetLineColor(kYellow-4);
		p_rareLepPt->SetFillColor(kYellow-4);
		p_qcdLepPt->SetFillStyle(1001);
		p_qcdLepPt->SetLineColor(kBlue);
		p_qcdLepPt->SetFillColor(kBlue);
		p_eleLepPt->SetFillStyle(1001);
		p_eleLepPt->SetLineColor(kRed);
		p_eleLepPt->SetFillColor(kRed);
		p_jetLepPt->SetFillStyle(1001);
		p_jetLepPt->SetLineColor(kGreen);
		p_jetLepPt->SetFillColor(kGreen);
		float events_eleLepPt, events_jetLepPt, events_qcdLepPt, events_VGLepPt, events_rareLepPt, events_allLepPt;
		events_eleLepPt =  p_eleLepPt->Integral();
		events_jetLepPt =  p_jetLepPt->Integral();
		events_qcdLepPt =  p_qcdLepPt->Integral();
		events_VGLepPt =   p_VGLepPt->Integral();
		events_rareLepPt = p_rareLepPt->Integral();
		events_allLepPt =  p_allLepPt->Integral();
		//cout<<events_eleLepPt <<"\t"<<events_jetLepPt<<"\t"<<events_qcdLepPt<<"\t"<<events_VGLepPt<<"\t"<<events_rareLepPt<<"\t"<<events_allLepPt<<endl;		
		
		std::vector<StackComponent> lepPtComponents = {
			{p_rareLepPt, "t#bar{t}#gamma/WW#gamma/WZ#gamma"},
			{p_eleLepPt, "e->#gamma fake"},
			{p_jetLepPt, "j->#gamma fake"},
			{p_qcdLepPt, "l fakes"},
			{p_VGLepPt, "W#gamma/Z#gamma"}
		};
		std::vector<StackComponent> lepPtStack = BuildIntegralOrderedStack(lepPtComponents, "leppt_stack");
		TH1F *p_totalLepPt = lepPtStack.back().hist;
		for(int ibin(1); ibin < p_totalLepPt->GetSize(); ibin++){
			error_LepPt->SetPoint(ibin-1,p_totalLepPt->GetBinCenter(ibin), p_totalLepPt->GetBinContent(ibin));
			float prederror = p_totalLepPt->GetBinError(ibin);
			error_LepPt->SetPointError(ibin-1,(p_totalLepPt->GetBinLowEdge(ibin+1)-p_totalLepPt->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_LepPt->SetPoint(ibin-1,p_totalLepPt->GetBinCenter(ibin), 1);
			ratioerror_LepPt->SetPointError(ibin-1,(p_totalLepPt->GetBinLowEdge(ibin+1)-p_totalLepPt->GetBinLowEdge(ibin))/2, prederror/p_totalLepPt->GetBinContent(ibin));
		}
		for(int idx = lepPtStack.size()-1; idx >= 0; --idx) lepPtStack[idx].hist->Draw("hist same");
		error_LepPt->SetFillColor(kBlack);
		error_LepPt->SetFillStyle(3345);
		error_LepPt->Draw("E2 same");
		TLegend *leg_leppt =  new TLegend(0.6,0.75,0.9,0.9);
		leg_leppt->SetFillStyle(0);
		gStyle->SetLegendBorderSize(1);
		gStyle->SetLegendFillColor(0);
		leg_mt->Draw("same");
		p_allLepPt->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( leppt_pad1,1, channel,11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( leppt_pad1,2, channel,11 );
		else if(RunYear==2017)  		CMS_lumi( leppt_pad1,3, channel,11 );
		else if(RunYear==2018)  		CMS_lumi( leppt_pad1,4, channel,11 );

		c_leppt->cd();
		TPad *leppt_pad2 = new TPad("leppt_pad2", "leppt_pad2", 0, 0, 1, 0.3);
		leppt_pad2->SetTopMargin(0);
		leppt_pad2->SetBottomMargin(0.4);
		leppt_pad2->Draw();
		leppt_pad2->cd();
		TLine *flatratio_leppt = new TLine(25,1,800,1);
		TH1F *ratio_leppt=(TH1F*)p_allLepPt->Clone("transfer factor");
		ratio_leppt->SetMarkerStyle(20);
		ratio_leppt->SetLineColor(kBlack);
		ratio_leppt->Divide(p_totalLepPt);
		ratio_leppt->Divide(p_VGLepPt);
		//cout<< p_allLepPt->Integral()/p_VGLepPt->Integral()<<endl;
		
		ratio_leppt->SetTitle("");
		ratio_leppt->GetYaxis()->SetTitle("obs./bkg.");
		ratio_leppt->GetXaxis()->SetTitleOffset(0.9);
		ratio_leppt->GetXaxis()->SetTitle("Lepton p_{T} (GeV)");
		ratio_leppt->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_leppt->GetXaxis()->SetLabelFont(63);
		ratio_leppt->GetXaxis()->SetLabelSize(14);
		ratio_leppt->GetYaxis()->SetLabelFont(63);
		ratio_leppt->GetYaxis()->SetLabelSize(11);
		ratio_leppt->Draw("ep");
		ratioerror_LepPt->SetFillColor(kBlack);
		ratioerror_LepPt->SetFillStyle(3345);
		ratioerror_LepPt->Draw("E2 same");
		flatratio_leppt->Draw("same");

		TLine *ratioValue_leppt = new TLine(25,p_allLepPt->Integral()/p_VGLepPt->Integral(),800,p_allLepPt->Integral()/p_VGLepPt->Integral());
		ratioValue_leppt->SetLineColor(kRed);
		//ratioValue_leppt->Draw("same");
		//c_leppt->SaveAs(ptplot.str().c_str());

		
		// ******** nJet ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_nJet = new TCanvas("nJet", "nJet",600,600);
		setCanvas(c_nJet); 
		c_nJet->cd();
		TPad *nJet_pad1 = new TPad("nJet_pad1", "nJet_pad1", 0, 0.3, 1, 1.0);
		setTopPad(nJet_pad1); 
		nJet_pad1->SetBottomMargin(0);
		nJet_pad1->Draw();  
		nJet_pad1->cd();  
		gPad->SetLogy();
		p_allnJet->SetTitle("p_{T}^{#gamma}");
		p_allnJet->SetMaximum(100*p_allnJet->GetBinContent(p_allnJet->GetMaximumBin()));
		p_allnJet->SetMinimum(0.5);
		p_allnJet->GetXaxis()->SetRangeUser(35,800);
		p_allnJet->SetLineColor(1);
		p_allnJet->SetMarkerStyle(20);
		p_allnJet->Draw("P");
		p_VGnJet->SetFillStyle(1001);
		p_VGnJet->SetLineColor(kMagenta);
		p_VGnJet->SetFillColor(kMagenta);
		p_rarenJet->SetFillStyle(1001);
		p_rarenJet->SetLineColor(kYellow-4);
		p_rarenJet->SetFillColor(kYellow-4);
		p_qcdnJet->SetFillStyle(1001);
		p_qcdnJet->SetLineColor(kBlue);
		p_qcdnJet->SetFillColor(kBlue);
		p_elenJet->SetFillStyle(1001);
		p_elenJet->SetLineColor(kRed);
		p_elenJet->SetFillColor(kRed);
		p_jetnJet->SetFillStyle(1001);
		p_jetnJet->SetLineColor(kGreen);
		p_jetnJet->SetFillColor(kGreen);
		//cout<<p_allnJet->GetBinError(3)<<"\t"<<p_elenJet->GetBinError(3)<<"\t"<<p_jetnJet->GetBinError(3)<<"\t"<<p_qcdnJet->GetBinError(3)<<"\t"<<p_VGnJet->GetBinError(3)<<endl;
		
		float events_elenJet, events_jetnJet, events_qcdnJet, events_VGnJet, events_rarenJet, events_allnJet;
		events_elenJet =  p_elenJet->Integral();
		events_jetnJet =  p_jetnJet->Integral();
		events_qcdnJet =  p_qcdnJet->Integral();
		events_VGnJet =   p_VGnJet->Integral();
		events_rarenJet = p_rarenJet->Integral();
		events_allnJet =  p_allnJet->Integral();
		cout<<events_elenJet <<"\t"<<events_jetnJet<<"\t"<<events_qcdnJet<<"\t"<<events_VGnJet<<"\t"<<events_rarenJet<<"\t"<<events_allnJet<<endl;		
		
		p_elenJet->Add(p_rarenJet); // ele 2nd
		p_jetnJet->Add(p_elenJet);  // jet 3rd
		p_qcdnJet->Add(p_jetnJet);  // qcd 4th
		p_VGnJet->Add(p_qcdnJet);   // VG  5th
		p_VGnJet->Sumw2();
		for(int ibin(1); ibin < p_VGnJet->GetSize(); ibin++){
			error_nJet->SetPoint(ibin-1,p_VGnJet->GetBinCenter(ibin), p_VGnJet->GetBinContent(ibin));
			float prederror = p_VGnJet->GetBinError(ibin);
			error_nJet->SetPointError(ibin-1,(p_VGnJet->GetBinLowEdge(ibin+1)-p_VGnJet->GetBinLowEdge(ibin))/2,prederror);
			//std::cout << p_elenJet->GetBinError(ibin) << " " << p_jetnJet->GetBinError(ibin) << " " << p_qcdnJet->GetBinError(ibin) << " " << p_rarenJet->GetBinError(ibin)  << std::endl;
			ratioerror_nJet->SetPoint(ibin-1,p_VGnJet->GetBinCenter(ibin), 1); 
			ratioerror_nJet->SetPointError(ibin-1,(p_VGnJet->GetBinLowEdge(ibin+1)-p_VGnJet->GetBinLowEdge(ibin))/2, prederror/p_VGnJet->GetBinContent(ibin)); 
			//cout<<prederror/p_VGnJet->GetBinContent(ibin)<<endl; 
		}
		p_VGnJet->Draw("hist same");
		p_qcdnJet->Draw("hist same");
		p_jetnJet->Draw("hist same");
		p_elenJet->Draw("hist same");
		p_rarenJet->Draw("hist same");
		error_nJet->SetFillColor(kBlack);
		error_nJet->SetFillStyle(3345);
		error_nJet->Draw("E2 same");
		//TLegend *leg_pt =  new TLegend(0.5,0.65,0.9,0.9);
		//leg_pt->SetFillStyle(0);
		//gStyle->SetLegendBorderSize(1);
		//gStyle->SetLegendFillColor(0);
		//leg_pt->AddEntry(p_allnJet,"observed (MT < 100 GeV)");
		//leg_pt->AddEntry(p_rarenJet,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
		//leg_pt->AddEntry(p_elenJet,"e->#gamma fake");
		//leg_pt->AddEntry(p_jetnJet,"j->#gamma fake");
		//leg_pt->AddEntry(p_qcdnJet,"j->e fake");
		//leg_pt->AddEntry(p_VGnJet, "WG/ZG");
		leg_mt->Draw("same");
		p_allnJet->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( nJet_pad1,1,channel, 11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( nJet_pad1,2,channel, 11 );
		else if(RunYear==2017)  		CMS_lumi( nJet_pad1,3,channel, 11 );
		else if(RunYear==2018)  		CMS_lumi( nJet_pad1,4,channel, 11 );

		c_nJet->cd();
		TPad *nJet_pad2 = new TPad("nJet_pad2", "nJet_pad2", 0, 0, 1, 0.3);
		nJet_pad2->SetTopMargin(0);
		nJet_pad2->SetBottomMargin(0.4);
		nJet_pad2->Draw();
		nJet_pad2->cd();
		TLine *flatratio_nJet = new TLine(0,1,10,1);
		TH1F *ratio_nJet = (TH1F*)p_allnJet->Clone("transfer factor");
				ratio_nJet->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_nJet->SetMarkerStyle(20);
		ratio_nJet->SetLineColor(kBlack);
		ratio_nJet->Divide(p_VGnJet);
		//cout<< "Discrepancy : "<<100*(p_VGnJet->Integral()-p_allnJet->Integral())/p_allnJet->Integral()<<" %" <<endl;
		ratio_nJet->SetTitle("");
		ratio_nJet->GetYaxis()->SetTitle("obs./bkg.");
		ratio_nJet->GetXaxis()->SetTitleOffset(0.9);
		ratio_nJet->GetXaxis()->SetTitle("Number of jets");
		ratio_nJet->GetXaxis()->SetLabelFont(63);
		ratio_nJet->GetXaxis()->SetLabelSize(14);
		ratio_nJet->GetYaxis()->SetLabelFont(63);
		ratio_nJet->GetYaxis()->SetLabelSize(11);
		ratio_nJet->Draw("ep");
		ratioerror_nJet->SetFillColor(kBlack);
		ratioerror_nJet->SetFillStyle(3345);
		ratioerror_nJet->Draw("E2 same");
		flatratio_nJet->Draw("same");
		TLine *ratioValue_nJet = new TLine(0,p_allnJet->Integral()/p_VGnJet->Integral(),10,p_allnJet->Integral()/p_VGnJet->Integral());
		ratioValue_nJet->SetLineColor(kRed);
		//ratioValue_nJet->Draw("same");
		//c_nJet->SaveAs(nJetplot.str().c_str());

		
		// ******** dPhi ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_dPhi = new TCanvas("dPhi", "dPhi",600,600);
		setCanvas(c_dPhi); 
		c_dPhi->cd();
		TPad *dPhi_pad1 = new TPad("dPhi_pad1", "dPhi_pad1", 0, 0.3, 1, 1.0);
		setTopPad(dPhi_pad1); 
		dPhi_pad1->SetBottomMargin(0);
		dPhi_pad1->Draw();  
		dPhi_pad1->cd();  
		gPad->SetLogy();
		p_alldPhiEleMET->SetTitle("p_{T}^{#gamma}");
		p_alldPhiEleMET->SetMaximum(100*p_alldPhiEleMET->GetBinContent(p_alldPhiEleMET->GetMaximumBin()));
		p_alldPhiEleMET->SetMinimum(0.5);
		p_alldPhiEleMET->GetXaxis()->SetRangeUser(35,800);
		p_alldPhiEleMET->SetLineColor(1);
		p_alldPhiEleMET->SetMarkerStyle(20);
		p_alldPhiEleMET->Draw("P");
		p_VGdPhiEleMET->SetFillStyle(1001);
		p_VGdPhiEleMET->SetLineColor(kMagenta);
		p_VGdPhiEleMET->SetFillColor(kMagenta);
		p_raredPhiEleMET->SetFillStyle(1001);
		p_raredPhiEleMET->SetLineColor(kYellow-4);
		p_raredPhiEleMET->SetFillColor(kYellow-4);
		p_qcddPhiEleMET->SetFillStyle(1001);
		p_qcddPhiEleMET->SetLineColor(kBlue);
		p_qcddPhiEleMET->SetFillColor(kBlue);
		p_eledPhiEleMET->SetFillStyle(1001);
		p_eledPhiEleMET->SetLineColor(kRed);
		p_eledPhiEleMET->SetFillColor(kRed);
		p_jetdPhiEleMET->SetFillStyle(1001);
		p_jetdPhiEleMET->SetLineColor(kGreen);
		p_jetdPhiEleMET->SetFillColor(kGreen);
		//cout<<p_alldPhiEleMET->GetBinError(3)<<"\t"<<p_eledPhiEleMET->GetBinError(3)<<"\t"<<p_jetdPhiEleMET->GetBinError(3)<<"\t"<<p_qcddPhiEleMET->GetBinError(3)<<"\t"<<p_VGdPhiEleMET->GetBinError(3)<<endl;
		
		float events_eledPhiEleMET, events_jetdPhiEleMET, events_qcddPhiEleMET, events_VGdPhiEleMET, events_raredPhiEleMET, events_alldPhiEleMET;
		events_eledPhiEleMET =  p_eledPhiEleMET->Integral();
		events_jetdPhiEleMET =  p_jetdPhiEleMET->Integral();
		events_qcddPhiEleMET =  p_qcddPhiEleMET->Integral();
		events_VGdPhiEleMET =   p_VGdPhiEleMET->Integral();
		events_raredPhiEleMET = p_raredPhiEleMET->Integral();
		events_alldPhiEleMET =  p_alldPhiEleMET->Integral();
		cout<<events_eledPhiEleMET <<"\t"<<events_jetdPhiEleMET<<"\t"<<events_qcddPhiEleMET<<"\t"<<events_VGdPhiEleMET<<"\t"<<events_raredPhiEleMET<<"\t"<<events_alldPhiEleMET<<endl;		
		
		p_eledPhiEleMET->Add(p_raredPhiEleMET); // ele 2nd
		p_jetdPhiEleMET->Add(p_eledPhiEleMET);  // jet 3rd
		p_qcddPhiEleMET->Add(p_jetdPhiEleMET);  // qcd 4th
		p_VGdPhiEleMET->Add(p_qcddPhiEleMET);   // VG  5th
		p_VGdPhiEleMET->Sumw2();
		for(int ibin(1); ibin < p_VGdPhiEleMET->GetSize(); ibin++){
			error_dPhiEleMET->SetPoint(ibin-1,p_VGdPhiEleMET->GetBinCenter(ibin), p_VGdPhiEleMET->GetBinContent(ibin));
			float prederror = p_VGdPhiEleMET->GetBinError(ibin);
			error_dPhiEleMET->SetPointError(ibin-1,(p_VGdPhiEleMET->GetBinLowEdge(ibin+1)-p_VGdPhiEleMET->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_dPhiEleMET->SetPoint(ibin-1,p_VGdPhiEleMET->GetBinCenter(ibin), 1); 
			ratioerror_dPhiEleMET->SetPointError(ibin-1,(p_VGdPhiEleMET->GetBinLowEdge(ibin+1)-p_VGdPhiEleMET->GetBinLowEdge(ibin))/2, prederror/p_VGdPhiEleMET->GetBinContent(ibin)); 
		}
		p_VGdPhiEleMET->Draw("hist same");
		p_qcddPhiEleMET->Draw("hist same");
		p_jetdPhiEleMET->Draw("hist same");
		p_eledPhiEleMET->Draw("hist same");
		p_raredPhiEleMET->Draw("hist same");
		error_dPhiEleMET->SetFillColor(kBlack);
		error_dPhiEleMET->SetFillStyle(3345);
		error_dPhiEleMET->Draw("E2 same");
		//TLegend *leg_pt =  new TLegend(0.5,0.65,0.9,0.9);
		//leg_pt->SetFillStyle(0);
		//gStyle->SetLegendBorderSize(1);
		//gStyle->SetLegendFillColor(0);
		//leg_pt->AddEntry(p_alldPhiEleMET,"observed (MT < 100 GeV)");
		//leg_pt->AddEntry(p_raredPhiEleMET,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
		//leg_pt->AddEntry(p_eledPhiEleMET,"e->#gamma fake");
		//leg_pt->AddEntry(p_jetdPhiEleMET,"j->#gamma fake");
		//leg_pt->AddEntry(p_qcddPhiEleMET,"j->e fake");
		//leg_pt->AddEntry(p_VGdPhiEleMET, "WG/ZG");
		leg_mt->Draw("same");
		p_alldPhiEleMET->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( dPhi_pad1,1,channel, 11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( dPhi_pad1,2,channel, 11 );
		else if(RunYear==2017)  		CMS_lumi( dPhi_pad1,3,channel, 11 );
		else if(RunYear==2018)  		CMS_lumi( dPhi_pad1,4,channel, 11 );

		c_dPhi->cd();
		TPad *dPhi_pad2 = new TPad("dPhi_pad2", "dPhi_pad2", 0, 0, 1, 0.3);
		dPhi_pad2->SetTopMargin(0);
		dPhi_pad2->SetBottomMargin(0.4);
		dPhi_pad2->Draw();
		dPhi_pad2->cd();
		TLine *flatratio_dPhi = new TLine(0,1,3.2,1);
		TH1F *ratio_dPhi = (TH1F*)p_alldPhiEleMET->Clone("transfer factor");
			ratio_dPhi->GetYaxis()->SetRangeUser(0.6,1.4);

		ratio_dPhi->SetMarkerStyle(20);
		ratio_dPhi->SetLineColor(kBlack);
		ratio_dPhi->Divide(p_VGdPhiEleMET);
		ratio_dPhi->SetTitle("");
		ratio_dPhi->GetYaxis()->SetTitle("obs./bkg.");
		ratio_dPhi->GetXaxis()->SetTitleOffset(0.9);
		ratio_dPhi->GetXaxis()->SetTitle("#Delta#phi(l, E_{T}^{miss}) (radians)");
		ratio_dPhi->GetXaxis()->SetLabelFont(63);
		ratio_dPhi->GetXaxis()->SetLabelSize(14);
		ratio_dPhi->GetYaxis()->SetLabelFont(63);
		ratio_dPhi->GetYaxis()->SetLabelSize(11);
		ratio_dPhi->Draw("ep");
		ratioerror_dPhiEleMET->SetFillColor(kBlack);
		ratioerror_dPhiEleMET->SetFillStyle(3345);
		ratioerror_dPhiEleMET->Draw("E2 same");
		flatratio_dPhi->Draw("same");
		TLine *ratioValue_dPhi = new TLine(0,p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral(),3.2,p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral());
		ratioValue_dPhi->SetLineColor(kRed);
		//ratioValue_dPhi->Draw("same");
		//c_dPhi->SaveAs(dPhiplot.str().c_str());
		
		// ******** nBJet ************************//
		gStyle->SetOptStat(0);
		TCanvas *c_nBJet = new TCanvas("nBJet", "nBJet",600,600);
		setCanvas(c_nBJet); 
		c_nBJet->cd();
		TPad *nBJet_pad1 = new TPad("nBJet_pad1", "nBJet_pad1", 0, 0.3, 1, 1.0);
		setTopPad(nBJet_pad1); 
		nBJet_pad1->SetBottomMargin(0);
		nBJet_pad1->Draw();  
		nBJet_pad1->cd();  
		gPad->SetLogy();
		p_allnBJet->GetXaxis()->SetRangeUser(0,900);
		p_allnBJet->SetMinimum(5);
		p_allnBJet->SetMaximum(10*p_allnBJet->GetBinContent(p_allnBJet->GetMaximumBin()));
		p_allnBJet->SetLineColor(1);
		p_allnBJet->SetMarkerStyle(20);
		p_allnBJet->Draw("P");
		p_VGnBJet->SetFillStyle(1001);
		p_VGnBJet->SetLineColor(kMagenta);
		p_VGnBJet->SetFillColor(kMagenta);
		p_rarenBJet->SetFillStyle(1001);
		p_rarenBJet->SetLineColor(kYellow-4);
		p_rarenBJet->SetFillColor(kYellow-4);
		p_qcdnBJet->SetFillStyle(1001);
		p_qcdnBJet->SetLineColor(kBlue);
		p_qcdnBJet->SetFillColor(kBlue);
		p_elenBJet->SetFillStyle(1001);
		p_elenBJet->SetLineColor(kRed);
		p_elenBJet->SetFillColor(kRed);
		p_jetnBJet->SetFillStyle(1001);
		p_jetnBJet->SetLineColor(kGreen);
		p_jetnBJet->SetFillColor(kGreen);
		
		p_elenBJet->Add(p_rarenBJet); // ele 2nd
		p_jetnBJet->Add(p_elenBJet);  // jet 3rd
		p_qcdnBJet->Add(p_jetnBJet);  // qcd 4th
		p_VGnBJet->Add(p_qcdnBJet);   // VG  5th
		p_VGnBJet->Sumw2();
		for(int ibin(1); ibin < p_VGnBJet->GetSize(); ibin++){
			error_nBJet->SetPoint(ibin-1,p_VGnBJet->GetBinCenter(ibin), p_VGnBJet->GetBinContent(ibin));
			float prederror = p_VGnBJet->GetBinError(ibin);
			error_nBJet->SetPointError(ibin-1,(p_VGnBJet->GetBinLowEdge(ibin+1)-p_VGnBJet->GetBinLowEdge(ibin))/2,prederror);
			ratioerror_nBJet->SetPoint(ibin-1,p_VGnBJet->GetBinCenter(ibin), 1); 
			ratioerror_nBJet->SetPointError(ibin-1,(p_VGnBJet->GetBinLowEdge(ibin+1)-p_VGnBJet->GetBinLowEdge(ibin))/2, prederror/p_VGnBJet->GetBinContent(ibin)); 
		}
		p_VGnBJet->Draw("hist same");
		p_qcdnBJet->Draw("hist same");
		p_jetnBJet->Draw("hist same");
		p_elenBJet->Draw("hist same");
		p_rarenBJet->Draw("hist same");
		error_nBJet->SetFillColor(kBlack);
		error_nBJet->SetFillStyle(3345);
		error_nBJet->Draw("E2 same");
		leg_mt->Draw("same");
		p_allnBJet->Draw("PE same");
		gPad->RedrawAxis();
		if(RunYear==2016 and preVFP == 1)  	CMS_lumi( nBJet_pad1,1, channel,11 );
		else if(RunYear==2016 and preVFP == 0)  CMS_lumi( nBJet_pad1,2, channel,11 );
		else if(RunYear==2017)  		CMS_lumi( nBJet_pad1,3, channel,11 );
		else if(RunYear==2018)  		CMS_lumi( nBJet_pad1,4, channel, 11 );

		c_nBJet->cd();
		TPad *nBJet_pad2 = new TPad("nBJet_pad2", "nBJet_pad2", 0, 0, 1, 0.3);
		nBJet_pad2->SetTopMargin(0);
		nBJet_pad2->SetBottomMargin(0.4);
		nBJet_pad2->Draw();
		nBJet_pad2->cd();
		TLine *flatratio_nBJet = new TLine(0,1,6,1);
		TH1F *ratio_nBJet=(TH1F*)p_allnBJet->Clone("transfer factor");
		ratio_nBJet->SetMarkerStyle(20);
		ratio_nBJet->SetLineColor(kBlack);
		ratio_nBJet->GetXaxis()->SetRangeUser(0,900);
		ratio_nBJet->GetYaxis()->SetRangeUser(0.6,1.4);
		ratio_nBJet->Divide(p_VGnBJet);
		ratio_nBJet->SetTitle("");
		ratio_nBJet->GetYaxis()->SetTitle("obs./bkg.");
		ratio_nBJet->GetXaxis()->SetTitleOffset(0.9);
		ratio_nBJet->GetXaxis()->SetTitle("Number of b-jets");
		ratio_nBJet->GetXaxis()->SetLabelFont(63);
		ratio_nBJet->GetXaxis()->SetLabelSize(14);
		ratio_nBJet->GetYaxis()->SetLabelFont(63);
		ratio_nBJet->GetYaxis()->SetLabelSize(11);

		ratio_nBJet->Draw("ep");
		ratioerror_nBJet->SetFillColor(kBlack);
		ratioerror_nBJet->SetFillStyle(3345);
		ratioerror_nBJet->Draw("E2 same");
		flatratio_nBJet->Draw("same");

		TLine *ratioValue_nbjet = new TLine(0,p_allnBJet->Integral()/p_VGnBJet->Integral(),5,p_allnBJet->Integral()/p_VGnBJet->Integral());
		ratioValue_nbjet->SetLineColor(kRed);
		//ratioValue_nbjet->Draw("same");
		//c_nBJet->SaveAs(nBJetplot.str().c_str());


		cout<<"\n\n\n";
		cout<<"=========================================================="<<endl;
		cout<<"\n";
		cout<< RunYear <<"  "<<whichVFP << "\tchannel : " << whichChannel <<"\t"; 
		if(plottype==1) cout<<"Control region	"<<endl;
		if(plottype==2) cout<<"Validation region	"<<endl;
		cout<<"\n\n\n";
		cout<<"MET, Data =" << p_allMET->Integral()<<"		data/MC = "<< p_allMET->Integral()/p_VGMET->Integral()<<endl;
		cout<<"MET_HT, Data =" << p_allMET_HT->Integral()<<"		data/MC = "<< p_allMET_HT->Integral()/p_VGMET_HT->Integral()<<endl;
		cout<<"MT,  Data =" << p_allMt->Integral()<<"           data/MC = "<< p_allMt->Integral()/p_VGMt->Integral()<<endl;
		cout<<"HT,  Data =" << p_allHT->Integral()<<"           data/MC = "<< p_allHT->Integral()/p_VGHT->Integral()<<endl;
		cout<<"Photon Pt, Data =" << p_allPhoEt->Integral()<<"           	data/MC = "<< p_allPhoEt->Integral()/p_VGPhoEt->Integral()<<endl;
		cout<<"Lepton Pt, Data =" << p_allLepPt->Integral()<<"           	data/MC = "<< p_allLepPt->Integral()/p_VGLepPt->Integral()<<endl;
		cout<<"dPhi, 	Data =" << p_alldPhiEleMET->Integral()<<"           	data/MC = "<< p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral()<<endl;
		cout<<"nJet, 	Data =" << p_allnJet->Integral()<<"           	data/MC = "<< p_allnJet->Integral()/p_VGnJet->Integral()<<endl;
		cout<<"\n";
		cout<<"=========================================================="<<endl;

	}
