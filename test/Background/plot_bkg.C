#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

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
#include "../../include/tdrstyle.C"

int channel = 2; // 1 = eg, 2 = mg
int plottype = 1; // 1 = bkg, 2 = valid
bool doTT = false;
int lowMET = 40;//0;
int highMET = 70;//-1;

//int RunYear=2016;
//bool preVFP=0;

void plot_bkg(int channel,int RunYear,bool preVFP){//main  
	gROOT->SetBatch(kTRUE);
  	gSystem->Load("../../lib/libAnaClasses.so");
	setTDRStyle();

	Double_t bkgEtBins[]={35,40,50,60,70,80,90,100,110,120,130,140,150,160,170,185,200,215,230,250,275,290, 305,325,345,370,400,500,800};
	int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
	Double_t bkgPtBins[]={25,50,75,100,125,150,200,400,800};
	int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
	Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
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
			signame << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_signal_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			elename << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_eleBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			jetname << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_jetbkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			qcdname << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_qcd_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			VGname  << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_VGBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			rarename << "/eos/uscms/store/user/tmishra/Background/bkgTree_egamma_rareBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
		}
		else if(plottype == 2){
			signame << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_signal_" << RunYear<<whichVFP <<".root";
			elename << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_eleBkg_" << RunYear<<whichVFP <<".root";
			jetname << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_jetbkg_" << RunYear<<whichVFP <<".root";
			qcdname << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_qcd_" << RunYear<<whichVFP <<".root";
			VGname  << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_VGBkg_" << RunYear<<whichVFP <<".root";
			rarename << "/eos/uscms/store/user/tmishra/Background/validTree_egamma_rareBkg_" << RunYear<<whichVFP <<".root";
		}
	}
	else if(channel == 2){
		if(plottype == 1){
			signame << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_signal_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			elename << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_eleBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			jetname << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_jetbkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			qcdname << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_qcd_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			VGname  << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_VGBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
			rarename << "/eos/uscms/store/user/tmishra/Background/bkgTree_mg_rareBkg_met"<< lowMET<<"_"<<highMET<<"_pt0_1000_" << RunYear<<whichVFP <<".root";
		}
		else if(plottype == 2){
			signame << "/eos/uscms/store/user/tmishra/Background/validTree_mg_signal_" << RunYear<<whichVFP <<".root";
			elename << "/eos/uscms/store/user/tmishra/Background/validTree_mg_eleBkg_" << RunYear<<whichVFP <<".root";
			jetname << "/eos/uscms/store/user/tmishra/Background/validTree_mg_jetbkg_" << RunYear<<whichVFP <<".root";
			qcdname << "/eos/uscms/store/user/tmishra/Background/validTree_mg_qcd_" << RunYear<<whichVFP <<".root";
			VGname  << "/eos/uscms/store/user/tmishra/Background/validTree_mg_VGBkg_" << RunYear<<whichVFP <<".root";
			rarename << "/eos/uscms/store/user/tmishra/Background/validTree_mg_rareBkg_" << RunYear<<whichVFP <<".root";
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

	TGraphErrors *error_dPhiEleMET = new TGraphErrors(32);
	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_dPhiEleMET = new TGraphErrors(32);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);

	TH1F *p_alldPhiEleMET = (TH1F*)file_sig->Get("p_dPhiEleMET");
	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");
	TH1F *p_allPU    = (TH1F*)file_sig->Get("p_PU");

	p_alldPhiEleMET->SetMarkerColor(1);
	p_allPhoEt->SetMarkerColor(1);
	p_allLepPt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allMt->SetMarkerColor(1);  
	p_allHT->SetMarkerColor(1);
	p_allPU->SetMarkerColor(1);
	
	TH1F *p_eledPhiEleMET = (TH1F*)file_ele->Get("p_dPhiEleMET");
	TH1F *p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt");
	TH1F *p_eleLepPt = (TH1F*)file_ele->Get("p_LepPt");
	TH1F *p_eleMET   = (TH1F*)file_ele->Get("p_MET");
	TH1F *p_eleMt    = (TH1F*)file_ele->Get("p_Mt");
	TH1F *p_eleHT  = (TH1F*)file_ele->Get("p_HT");
	TH1F *p_elePU    = (TH1F*)file_ele->Get("p_PU");

	TH1F *p_jetdPhiEleMET = (TH1F*)file_jet->Get("p_dPhiEleMET");
	TH1F *p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt");
	TH1F *p_jetLepPt = (TH1F*)file_jet->Get("p_LepPt");
	TH1F *p_jetMET   = (TH1F*)file_jet->Get("p_MET");
	TH1F *p_jetMt    = (TH1F*)file_jet->Get("p_Mt");
	TH1F *p_jetHT  = (TH1F*)file_jet->Get("p_HT");
	TH1F *p_jetPU    = (TH1F*)file_jet->Get("p_PU");

	TH1F *p_qcddPhiEleMET = (TH1F*)file_qcd->Get("p_dPhiEleMET");
	TH1F *p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt");
	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	TH1F *p_qcdMET   = (TH1F*)file_qcd->Get("p_MET");
	TH1F *p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt");
	TH1F *p_qcdHT  = (TH1F*)file_qcd->Get("p_HT");
	TH1F *p_qcdPU    = (TH1F*)file_qcd->Get("p_PU");

	TH1F *p_VGdPhiEleMET = (TH1F*)file_VG->Get("p_dPhiEleMET");
	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGHT  = (TH1F*)file_VG->Get("p_HT");
	TH1F *p_VGPU    = (TH1F*)file_VG->Get("p_PU");

	TH1F *p_raredPhiEleMET = (TH1F*)file_rare->Get("p_dPhiEleMET");
	TH1F *p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt");
	TH1F *p_rareLepPt = (TH1F*)file_rare->Get("p_LepPt");
	TH1F *p_rareMET   = (TH1F*)file_rare->Get("p_MET");
	TH1F *p_rareMt    = (TH1F*)file_rare->Get("p_Mt");
	TH1F *p_rareHT  = (TH1F*)file_rare->Get("p_HT");
	TH1F *p_rarePU    = (TH1F*)file_rare->Get("p_PU");

	if(doTT){
		p_alldPhiEleMET = (TH1F*)file_sig->Get("p_dPhiEleMET_TT");
		p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt_TT");
		p_allMET   = (TH1F*)file_sig->Get("p_MET_TT");
		p_allMt    = (TH1F*)file_sig->Get("p_Mt_TT");
		p_allHT  = (TH1F*)file_sig->Get("p_HT_TT");
		
		p_eledPhiEleMET = (TH1F*)file_ele->Get("p_dPhiEleMET_TT");
		p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt_TT");
		p_eleMET   = (TH1F*)file_ele->Get("p_MET_TT");
		p_eleMt    = (TH1F*)file_ele->Get("p_Mt_TT");
		p_eleHT  = (TH1F*)file_ele->Get("p_HT_TT");
	
		p_jetdPhiEleMET = (TH1F*)file_jet->Get("p_dPhiEleMET_TT");
		p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt_TT");
		p_jetMET   = (TH1F*)file_jet->Get("p_MET_TT");
		p_jetMt    = (TH1F*)file_jet->Get("p_Mt_TT");
		p_jetHT  = (TH1F*)file_jet->Get("p_HT_TT");
	
		p_qcddPhiEleMET = (TH1F*)file_qcd->Get("p_dPhiEleMET_TT");
		p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt_TT");
		p_qcdMET   = (TH1F*)file_qcd->Get("p_MET_TT");
		p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt_TT");
		p_qcdHT  = (TH1F*)file_qcd->Get("p_HT_TT");
	
		p_VGdPhiEleMET = (TH1F*)file_VG->Get("p_dPhiEleMET_TT");
		p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt_TT");
		p_VGMET   = (TH1F*)file_VG->Get("p_MET_TT");
		p_VGMt    = (TH1F*)file_VG->Get("p_Mt_TT");
		p_VGHT  = (TH1F*)file_VG->Get("p_HT_TT");
	
		p_raredPhiEleMET = (TH1F*)file_rare->Get("p_dPhiEleMET_TT");
		p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt_TT");
		p_rareMET   = (TH1F*)file_rare->Get("p_MET_TT");
		p_rareMt    = (TH1F*)file_rare->Get("p_Mt_TT");
		p_rareHT  = (TH1F*)file_rare->Get("p_HT_TT");
	}


//	int binnumber;
// binnumber = p_allPhoEt->GetSize()-2;  p_allPhoEt->SetBinContent(binnumber,  p_allPhoEt->GetBinContent(binnumber) +   p_allPhoEt->GetBinContent(binnumber+1) );   
// binnumber = p_allLepPt->GetSize()-2;  p_allLepPt->SetBinContent(binnumber,  p_allLepPt->GetBinContent(binnumber) +   p_allLepPt->GetBinContent(binnumber+1) ); 
// binnumber = p_allMET->GetSize()-2;    p_allMET->SetBinContent(binnumber,    p_allMET->GetBinContent(binnumber) +     p_allMET->GetBinContent(binnumber+1) );   
// binnumber = p_allMt->GetSize()-2;     p_allMt->SetBinContent(binnumber,     p_allMt->GetBinContent(binnumber) +      p_allMt->GetBinContent(binnumber+1) );    
// binnumber = p_allHT->GetSize()-2;     p_allHT->SetBinContent(binnumber,     p_allHT->GetBinContent(binnumber) +      p_allHT->GetBinContent(binnumber+1) );    
//                                                                                                                          
// binnumber = p_elePhoEt->GetSize()-2;  p_elePhoEt->SetBinContent(binnumber,  p_elePhoEt->GetBinContent(binnumber) +   p_elePhoEt->GetBinContent(binnumber+1) ); 
// binnumber = p_eleLepPt->GetSize()-2;  p_eleLepPt->SetBinContent(binnumber,  p_eleLepPt->GetBinContent(binnumber) +   p_eleLepPt->GetBinContent(binnumber+1) ); 
// binnumber = p_eleMET->GetSize()-2;    p_eleMET->SetBinContent(binnumber,    p_eleMET->GetBinContent(binnumber) +     p_eleMET->GetBinContent(binnumber+1) );   
// binnumber = p_eleMt->GetSize()-2;     p_eleMt->SetBinContent(binnumber,     p_eleMt->GetBinContent(binnumber) +      p_eleMt->GetBinContent(binnumber+1) );    
// binnumber = p_eleHT->GetSize()-2;     p_eleHT->SetBinContent(binnumber,     p_eleHT->GetBinContent(binnumber) +      p_eleHT->GetBinContent(binnumber+1) );    
//                                                                                                                                                                           
// binnumber = p_jetPhoEt->GetSize()-2;  p_jetPhoEt->SetBinContent(binnumber,  p_jetPhoEt->GetBinContent(binnumber) +   p_jetPhoEt->GetBinContent(binnumber+1) );             
// binnumber = p_jetLepPt->GetSize()-2;  p_jetLepPt->SetBinContent(binnumber,  p_jetLepPt->GetBinContent(binnumber) +   p_jetLepPt->GetBinContent(binnumber+1) );             
// binnumber = p_jetMET->GetSize()-2;    p_jetMET->SetBinContent(binnumber,    p_jetMET->GetBinContent(binnumber) +     p_jetMET->GetBinContent(binnumber+1) );               
// binnumber = p_jetMt->GetSize()-2;     p_jetMt->SetBinContent(binnumber,     p_jetMt->GetBinContent(binnumber) +      p_jetMt->GetBinContent(binnumber+1) );                
// binnumber = p_jetHT->GetSize()-2;     p_jetHT->SetBinContent(binnumber,     p_jetHT->GetBinContent(binnumber) +      p_jetHT->GetBinContent(binnumber+1) );  
//                                                                                                                                                                           
// binnumber = p_qcdPhoEt->GetSize()-2;  p_qcdPhoEt->SetBinContent(binnumber,  p_qcdPhoEt->GetBinContent(binnumber) +   p_qcdPhoEt->GetBinContent(binnumber+1) );  
// binnumber = p_qcdLepPt->GetSize()-2;  p_qcdLepPt->SetBinContent(binnumber,  p_qcdLepPt->GetBinContent(binnumber) +   p_qcdLepPt->GetBinContent(binnumber+1) );  
// binnumber = p_qcdMET->GetSize()-2;    p_qcdMET->SetBinContent(binnumber,    p_qcdMET->GetBinContent(binnumber) +     p_qcdMET->GetBinContent(binnumber+1) );    
// binnumber = p_qcdMt->GetSize()-2;     p_qcdMt->SetBinContent(binnumber,     p_qcdMt->GetBinContent(binnumber) +      p_qcdMt->GetBinContent(binnumber+1) );     
// binnumber = p_qcdHT->GetSize()-2;     p_qcdHT->SetBinContent(binnumber,     p_qcdHT->GetBinContent(binnumber) +      p_qcdHT->GetBinContent(binnumber+1) );  
//                                                                                                                                                                           
// binnumber = p_VGPhoEt->GetSize()-2;   p_VGPhoEt->SetBinContent(binnumber,   p_VGPhoEt->GetBinContent(binnumber) +    p_VGPhoEt->GetBinContent(binnumber+1) );   
// binnumber = p_VGLepPt->GetSize()-2;   p_VGLepPt->SetBinContent(binnumber,   p_VGLepPt->GetBinContent(binnumber) +    p_VGLepPt->GetBinContent(binnumber+1) );   
// binnumber = p_VGMET->GetSize()-2;     p_VGMET->SetBinContent(binnumber,     p_VGMET->GetBinContent(binnumber) +      p_VGMET->GetBinContent(binnumber+1) );     
// binnumber = p_VGMt->GetSize()-2;      p_VGMt->SetBinContent(binnumber,      p_VGMt->GetBinContent(binnumber) +       p_VGMt->GetBinContent(binnumber+1) );      
// binnumber = p_VGHT->GetSize()-2;      p_VGHT->SetBinContent(binnumber,      p_VGHT->GetBinContent(binnumber) +       p_VGHT->GetBinContent(binnumber+1) );   
//                                                                                                                                                                           
// binnumber = p_rarePhoEt->GetSize()-2; p_rarePhoEt->SetBinContent(binnumber, p_rarePhoEt->GetBinContent(binnumber) +  p_rarePhoEt->GetBinContent(binnumber+1) );            
// binnumber = p_rareLepPt->GetSize()-2; p_rareLepPt->SetBinContent(binnumber, p_rareLepPt->GetBinContent(binnumber) +  p_rareLepPt->GetBinContent(binnumber+1) );            
// binnumber = p_rareMET->GetSize()-2;   p_rareMET->SetBinContent(binnumber,   p_rareMET->GetBinContent(binnumber) +    p_rareMET->GetBinContent(binnumber+1) );              
// binnumber = p_rareMt->GetSize()-2;    p_rareMt->SetBinContent(binnumber,    p_rareMt->GetBinContent(binnumber) +     p_rareMt->GetBinContent(binnumber+1) );               
// binnumber = p_rareHT->GetSize()-2;    p_rareHT->SetBinContent(binnumber,    p_rareHT->GetBinContent(binnumber) +     p_rareHT->GetBinContent(binnumber+1) );                  
//

	std::ostringstream dPhiplot;  dPhiplot.str("");
	std::ostringstream etplot;    etplot.str("");
	std::ostringstream ptplot;    ptplot.str("");
	std::ostringstream metplot;   metplot.str("");
	std::ostringstream mtplot;    mtplot.str("");
	std::ostringstream htplot;    htplot.str("");
	
	if(doTT){
	   if(channel == 1){
		if(plottype == 1){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_met_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
	   }
	   else if(channel == 2){
		if(plottype == 1){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_dPhi_TT_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_pt_TT_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_leppt_TT_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_met_TT_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_mt_TT_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_ht_TT_" << RunYear<<whichVFP <<".png";	
		}
	   }
	}
	else{
	   if(channel == 1){
		if(plottype == 1){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_egamma_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_egamma_ht_" << RunYear<<whichVFP <<".png";	
		}
	   }
	   else if(channel == 2){
		if(plottype == 1){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"BKG_mg_ht_" << RunYear<<whichVFP <<".png";	
		}
		else if(plottype == 2){
			dPhiplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_dPhi_" << RunYear<<whichVFP <<".png";
			etplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_pt_" << RunYear<<whichVFP <<".png";
			ptplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_leppt_" << RunYear<<whichVFP <<".png";		
			metplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_met_" << RunYear<<whichVFP <<".png";
			mtplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_mt_" << RunYear<<whichVFP <<".png";
			htplot << "/eos/uscms/store/user/tmishra/Background/plots/"<<RunYear<<whichVFP<<"/"<<whichChannel<<"/"<<"VALID_mg_ht_" << RunYear<<whichVFP <<".png";	
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
	//cout<<p_eleMt->GetBinError(3)<<"\t"<<p_jetMt->GetBinError(3)<<"\t"<<p_qcdMt->GetBinError(3)<<"\t"<<p_VGMt->GetBinError(3)<<endl;
	
	float events_eleMt, events_jetMt, events_qcdMt, events_VGMt, events_rareMt, events_allMt;
	events_eleMt =  p_eleMt->Integral();
	events_jetMt =  p_jetMt->Integral();
	events_qcdMt =  p_qcdMt->Integral();
	events_VGMt =   p_VGMt->Integral();
	events_rareMt = p_rareMt->Integral();
	events_allMt =  p_allMt->Integral();
	cout<<events_eleMt <<"\t"<<events_jetMt<<"\t"<<events_qcdMt<<"\t"<<events_VGMt<<"\t"<<events_rareMt<<"\t"<<events_allMt<<endl;		

	p_eleMt->Add(p_rareMt); // ele 2nd
	p_jetMt->Add(p_eleMt);  // jet 3rd
	p_qcdMt->Add(p_jetMt);  // qcd 4th
	p_VGMt->Add(p_qcdMt);   // VG  5th
	for(int ibin(1); ibin < p_VGMt->GetSize(); ibin++){
		//p_VGMt->SetBinError(ibin-1,0);
		error_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), p_VGMt->GetBinContent(ibin));
		float prederror = p_VGMt->GetBinError(ibin);
//		prederror += p_eleMt->GetBinError(ibin);
//		prederror += p_jetMt->GetBinError(ibin);
//		prederror += p_qcdMt->GetBinError(ibin);
//		prederror += p_rareMt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2, prederror/p_VGMt->GetBinContent(ibin)); 
	}
	p_VGMt->Draw("hist same");
	p_qcdMt->Draw("hist same");
	p_jetMt->Draw("hist same");
	p_eleMt->Draw("hist same");
	p_rareMt->Draw("hist same");
  	error_Mt->SetFillColor(kBlack);
  	error_Mt->SetFillStyle(3345);
	//error_Mt->Draw("E2 same");
	TLegend *leg_mt =  new TLegend(0.4,0.6,0.9,0.9);
	leg_mt->SetTextSize(13);
	leg_mt->SetNColumns(2);
	leg_mt->SetFillStyle(0);
	leg_mt->SetBorderSize(0);
	leg_mt->SetFillColor(0);
	leg_mt->AddEntry(p_allMt,Form("observed [%.1f]",events_allMt));
	leg_mt->AddEntry(p_rareMt,Form("t#bar{t}#gamma/WW#gamma/WZ#gamma [%.1f]",events_rareMt));
	leg_mt->AddEntry(p_eleMt,Form("e->#gamma fake [%.1f]",events_eleMt));
	leg_mt->AddEntry(p_jetMt,Form("j->#gamma fake [%.1f]",events_jetMt));
	leg_mt->AddEntry(p_qcdMt,Form("j->e fake [%.1f]", events_qcdMt));
	leg_mt->AddEntry(p_VGMt, Form("WG/ZG [%.1f]",events_VGMt));
	//leg_mt->AddEntry(ratioerror_Mt, "Unc");
	leg_mt->Draw("same");
	p_allMt->Draw("E same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( mt_pad1,1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( mt_pad1,2, 11 );
	else if(RunYear==2017)  		CMS_lumi( mt_pad1,3, 11 );
	else if(RunYear==2018)  		CMS_lumi( mt_pad1,4, 11 );

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
	//ratio_mt->GetXaxis()->SetRangeUser(0,800);
	//ratio_mt->SetMinimum(0);
	//ratio_mt->SetMaximum(2);
	ratio_mt->Divide(p_VGMt);
	//cout<< p_allMt->Integral()/p_VGMt->Integral()<<endl;
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("obs./bkg.");
	ratio_mt->GetYaxis()->SetRangeUser(0.7,1.23);
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(14);
	ratio_mt->Draw("pe");
	//ratioerror_Mt->SetFillColor(kBlack);
	//ratioerror_Mt->SetFillStyle(3345);
	//ratioerror_Mt->Draw("E2 same");
	flatratio_mt->Draw("same");

  	TLine *ratioValue_mt = new TLine(0,p_allMt->Integral()/p_VGMt->Integral(),1000,p_allMt->Integral()/p_VGMt->Integral());
	ratioValue_mt->SetLineColor(kRed);
	ratioValue_mt->Draw("same");
	c_mt->SaveAs(mtplot.str().c_str());
 
	// ******** PhoEt ************************//
  	gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	setCanvas(c_pt); 
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(pt_pad1); 
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetTitle("p_{T}^{#gamma}");
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
	//cout<<p_allPhoEt->GetBinError(3)<<"\t"<<p_elePhoEt->GetBinError(3)<<"\t"<<p_jetPhoEt->GetBinError(3)<<"\t"<<p_qcdPhoEt->GetBinError(3)<<"\t"<<p_VGPhoEt->GetBinError(3)<<endl;
	
	float events_elePhoEt, events_jetPhoEt, events_qcdPhoEt, events_VGPhoEt, events_rarePhoEt, events_allPhoEt;
	events_elePhoEt =  p_elePhoEt->Integral();
	events_jetPhoEt =  p_jetPhoEt->Integral();
	events_qcdPhoEt =  p_qcdPhoEt->Integral();
	events_VGPhoEt =   p_VGPhoEt->Integral();
	events_rarePhoEt = p_rarePhoEt->Integral();
	events_allPhoEt =  p_allPhoEt->Integral();
	cout<<events_elePhoEt <<"\t"<<events_jetPhoEt<<"\t"<<events_qcdPhoEt<<"\t"<<events_VGPhoEt<<"\t"<<events_rarePhoEt<<"\t"<<events_allPhoEt<<endl;		
	
	p_elePhoEt->Add(p_rarePhoEt); // ele 2nd
	p_jetPhoEt->Add(p_elePhoEt);  // jet 3rd
	p_qcdPhoEt->Add(p_jetPhoEt);  // qcd 4th
	p_VGPhoEt->Add(p_qcdPhoEt);   // VG  5th
	//p_VGPhoEt->Sumw2();
	for(int ibin(1); ibin < p_VGPhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), p_VGPhoEt->GetBinContent(ibin));
		float prederror = p_VGPhoEt->GetBinError(ibin);
		//prederror += p_elePhoEt->GetBinError(ibin);
		//prederror += p_jetPhoEt->GetBinError(ibin);
		//prederror += p_qcdPhoEt->GetBinError(ibin);
		//prederror += p_rarePhoEt->GetBinError(ibin)*0.6;
		error_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2,prederror);
		//std::cout << p_elePhoEt->GetBinError(ibin) << " " << p_jetPhoEt->GetBinError(ibin) << " " << p_qcdPhoEt->GetBinError(ibin) << " " << p_rarePhoEt->GetBinError(ibin)  << std::endl;
		ratioerror_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2, prederror/p_VGPhoEt->GetBinContent(ibin)); 
		//cout<<prederror/p_VGPhoEt->GetBinContent(ibin)<<endl; 
	}
	p_VGPhoEt->Draw("hist same");
	p_qcdPhoEt->Draw("hist same");
	p_jetPhoEt->Draw("hist same");
	p_elePhoEt->Draw("hist same");
	p_rarePhoEt->Draw("hist same");
  	error_PhoEt->SetFillColor(kBlack);
  	error_PhoEt->SetFillStyle(3345);
	//error_PhoEt->Draw("E2 same");
	//TLegend *leg_pt =  new TLegend(0.5,0.65,0.9,0.9);
	//leg_pt->SetFillStyle(0);
	//gStyle->SetLegendBorderSize(1);
	//gStyle->SetLegendFillColor(0);
	//leg_pt->AddEntry(p_allPhoEt,"observed (MT < 100 GeV)");
	//leg_pt->AddEntry(p_rarePhoEt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	//leg_pt->AddEntry(p_elePhoEt,"e->#gamma fake");
	//leg_pt->AddEntry(p_jetPhoEt,"j->#gamma fake");
	//leg_pt->AddEntry(p_qcdPhoEt,"j->e fake");
	//leg_pt->AddEntry(p_VGPhoEt, "WG/ZG");
	leg_mt->Draw("same");
	p_allPhoEt->Draw("E same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( pt_pad1,1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( pt_pad1,2, 11 );
	else if(RunYear==2017)  		CMS_lumi( pt_pad1,3, 11 );
	else if(RunYear==2018)  		CMS_lumi( pt_pad1,4, 11 );

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.3);
	pt_pad2->SetTopMargin(0);
	pt_pad2->SetBottomMargin(0.4);
	pt_pad2->Draw();
	pt_pad2->cd();
  	TLine *flatratio = new TLine(35,1,800,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	ratio->SetMinimum(0.7);
	ratio->SetMaximum(1.23);
	//ratio->SetMinimum(0);
	//ratio->SetMaximum(2);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_VGPhoEt);
	//cout<< "Discrepancy : "<<100*(p_VGPhoEt->Integral()-p_allPhoEt->Integral())/p_allPhoEt->Integral()<<" %" <<endl;
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("obs./bkg.");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(14);
	ratio->Draw("ep");
	//ratioerror_PhoEt->SetFillColor(kBlack);
	//ratioerror_PhoEt->SetFillStyle(3345);
	//ratioerror_PhoEt->Draw("E2 same");
	flatratio->Draw("same");
  	TLine *ratioValue = new TLine(35,p_allPhoEt->Integral()/p_VGPhoEt->Integral(),800,p_allPhoEt->Integral()/p_VGPhoEt->Integral());
	ratioValue->SetLineColor(kRed);
	ratioValue->Draw("same");
	c_pt->SaveAs(etplot.str().c_str());

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
	//p_VGdPhiEleMET->Sumw2();
	for(int ibin(1); ibin < p_VGdPhiEleMET->GetSize(); ibin++){
		error_dPhiEleMET->SetPoint(ibin-1,p_VGdPhiEleMET->GetBinCenter(ibin), p_VGdPhiEleMET->GetBinContent(ibin));
		float prederror = p_VGdPhiEleMET->GetBinError(ibin);
		//prederror += p_eledPhiEleMET->GetBinError(ibin);
		//prederror += p_jetdPhiEleMET->GetBinError(ibin);
		//prederror += p_qcddPhiEleMET->GetBinError(ibin);
		//prederror += p_raredPhiEleMET->GetBinError(ibin)*0.6;
		error_dPhiEleMET->SetPointError(ibin-1,(p_VGdPhiEleMET->GetBinLowEdge(ibin+1)-p_VGdPhiEleMET->GetBinLowEdge(ibin))/2,prederror);
		//std::cout << p_eledPhiEleMET->GetBinError(ibin) << " " << p_jetdPhiEleMET->GetBinError(ibin) << " " << p_qcddPhiEleMET->GetBinError(ibin) << " " << p_raredPhiEleMET->GetBinError(ibin)  << std::endl;
		ratioerror_dPhiEleMET->SetPoint(ibin-1,p_VGdPhiEleMET->GetBinCenter(ibin), 1); 
		ratioerror_dPhiEleMET->SetPointError(ibin-1,(p_VGdPhiEleMET->GetBinLowEdge(ibin+1)-p_VGdPhiEleMET->GetBinLowEdge(ibin))/2, prederror/p_VGdPhiEleMET->GetBinContent(ibin)); 
		//cout<<prederror/p_VGdPhiEleMET->GetBinContent(ibin)<<endl; 
	}
	p_VGdPhiEleMET->Draw("hist same");
	p_qcddPhiEleMET->Draw("hist same");
	p_jetdPhiEleMET->Draw("hist same");
	p_eledPhiEleMET->Draw("hist same");
	p_raredPhiEleMET->Draw("hist same");
  	error_dPhiEleMET->SetFillColor(kBlack);
  	error_dPhiEleMET->SetFillStyle(3345);
	//error_dPhiEleMET->Draw("E2 same");
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
	p_alldPhiEleMET->Draw("E same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( dPhi_pad1,1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( dPhi_pad1,2, 11 );
	else if(RunYear==2017)  		CMS_lumi( dPhi_pad1,3, 11 );
	else if(RunYear==2018)  		CMS_lumi( dPhi_pad1,4, 11 );

	c_dPhi->cd();
	TPad *dPhi_pad2 = new TPad("dPhi_pad2", "dPhi_pad2", 0, 0, 1, 0.3);
	dPhi_pad2->SetTopMargin(0);
	dPhi_pad2->SetBottomMargin(0.4);
	dPhi_pad2->Draw();
	dPhi_pad2->cd();
  	TLine *flatratio_dPhi = new TLine(0,1,3.2,1);
	TH1F *ratio_dPhi = (TH1F*)p_alldPhiEleMET->Clone("transfer factor");
	ratio_dPhi->SetMinimum(0.7);
	ratio_dPhi->SetMaximum(1.23);
	//ratio->SetMinimum(0);
	//ratio->SetMaximum(2);
	ratio_dPhi->SetMarkerStyle(20);
	ratio_dPhi->SetLineColor(kBlack);
	ratio_dPhi->Divide(p_VGdPhiEleMET);
	//cout<< "Discrepancy : "<<100*(p_VGdPhiEleMET->Integral()-p_alldPhiEleMET->Integral())/p_alldPhiEleMET->Integral()<<" %" <<endl;
	ratio_dPhi->SetTitle("");
	ratio_dPhi->GetYaxis()->SetTitle("obs./bkg.");
	ratio_dPhi->GetXaxis()->SetTitle("#Delta#phi(l, E_{T}^{miss}) (radians)");
	ratio_dPhi->GetXaxis()->SetLabelFont(63);
	ratio_dPhi->GetXaxis()->SetLabelSize(14);
	ratio_dPhi->GetYaxis()->SetLabelFont(63);
	ratio_dPhi->GetYaxis()->SetLabelSize(14);
	ratio_dPhi->Draw("ep");
	//ratioerror_dPhiEleMET->SetFillColor(kBlack);
	//ratioerror_dPhiEleMET->SetFillStyle(3345);
	//ratioerror_dPhiEleMET->Draw("E2 same");
	flatratio_dPhi->Draw("same");
  	TLine *ratioValue_dPhi = new TLine(0,p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral(),3.2,p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral());
	ratioValue_dPhi->SetLineColor(kRed);
	ratioValue_dPhi->Draw("same");
	c_dPhi->SaveAs(dPhiplot.str().c_str());

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
	//p_allMET->GetYaxis()->SetRangeUser(0.01, 1000);
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
	//cout<<p_eleMET->GetBinError(3)<<"\t"<<p_jetMET->GetBinError(3)<<"\t"<<p_qcdMET->GetBinError(3)<<"\t"<<p_VGMET->GetBinError(3)<<endl;
	float events_eleMET, events_jetMET, events_qcdMET, events_VGMET, events_rareMET, events_allMET;
	events_eleMET =  p_eleMET->Integral();
	events_jetMET =  p_jetMET->Integral();
	events_qcdMET =  p_qcdMET->Integral();
	events_VGMET =   p_VGMET->Integral();
	events_rareMET = p_rareMET->Integral();
	events_allMET =  p_allMET->Integral();
	cout<<events_eleMET <<"\t"<<events_jetMET<<"\t"<<events_qcdMET<<"\t"<<events_VGMET<<"\t"<<events_rareMET<<"\t"<<events_allMET<<endl;		
	
	p_eleMET->Add(p_rareMET); // ele 2nd
	p_jetMET->Add(p_eleMET);  // jet 3rd
	p_qcdMET->Add(p_jetMET);  // qcd 4th
	p_VGMET->Add(p_qcdMET);   // VG  5th
	//p_VGMET->Sumw2();
	for(int ibin(1); ibin < p_VGMET->GetSize(); ibin++){
		float prederror = p_VGMET->GetBinError(ibin);
		//	prederror += p_eleMET->GetBinError(ibin);
		//	prederror += p_jetMET->GetBinError(ibin);
		//	prederror += p_qcdMET->GetBinError(ibin);
		//	prederror += p_rareMET->GetBinError(ibin);
		error_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), p_VGMET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), 1); 
		ratioerror_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2, prederror/p_VGMET->GetBinContent(ibin)); 
	}
	p_VGMET->Draw("hist same");
	p_qcdMET->Draw("hist same");
	p_jetMET->Draw("hist same");
	p_eleMET->Draw("hist same");
	p_rareMET->Draw("hist same");
  	error_MET->SetFillColor(kBlack);
  	error_MET->SetFillStyle(3345);
	//error_MET->Draw("E2 same");
	//TLegend *leg_met =  new TLegend(0.5,0.65,0.9,0.9);
	//leg_met->SetFillStyle(0);
	//gStyle->SetLegendBorderSize(1);
	//gStyle->SetLegendFillColor(0);
	//leg_met->AddEntry(p_allMET,"observed (MT < 100 GeV)");
	//leg_met->AddEntry(p_rareMET,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	//leg_met->AddEntry(p_eleMET,"e->#gamma fake");
	//leg_met->AddEntry(p_jetMET,"j->#gamma fake");
	//leg_met->AddEntry(p_qcdMET,"j->e fake");
	//leg_met->AddEntry(p_VGMET, "WG/ZG");
	leg_mt->Draw("same");
	p_allMET->Draw("E same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( met_pad1, 1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( met_pad1, 2, 11 );
	else if(RunYear==2017)  		CMS_lumi( met_pad1, 3, 11 );
	else if(RunYear==2018)  		CMS_lumi( met_pad1, 4, 11 );

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
	met_pad2->SetTopMargin(0);
	met_pad2->SetBottomMargin(0.4);
	met_pad2->Draw();
	met_pad2->cd();
  	TLine *flatratio_met = new TLine(0,1,1000,1);
	TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
	//ratio_met->GetXaxis()->SetRangeUser(120,1000);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(p_VGMET);
	//cout<< p_allMET->Integral()/p_VGMET->Integral()<<endl;
	ratio_met->SetTitle("");
	ratio_met->GetYaxis()->SetTitle("obs./bkg.");
	//ratio_met->GetYaxis()->SetRangeUser(0,2.1);
	ratio_met->GetYaxis()->SetRangeUser(0.7,1.23);
	ratio_met->GetXaxis()->SetLabelFont(63);
	ratio_met->GetXaxis()->SetLabelSize(14);
	ratio_met->GetYaxis()->SetLabelFont(63);
	ratio_met->GetYaxis()->SetLabelSize(14);
	ratio_met->Draw("pe");
	//ratioerror_MET->SetFillColor(kBlack);
	//ratioerror_MET->SetFillStyle(3345);
	//ratioerror_MET->Draw("E2 same");
	flatratio_met->Draw("same");

  	TLine *ratioValue_met = new TLine(0,p_allMET->Integral()/p_VGMET->Integral(),1000,p_allMET->Integral()/p_VGMET->Integral());
	ratioValue_met->SetLineColor(kRed);
	ratioValue_met->Draw("same");

	c_met->SaveAs(metplot.str().c_str());


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
	cout<<events_eleLepPt <<"\t"<<events_jetLepPt<<"\t"<<events_qcdLepPt<<"\t"<<events_VGLepPt<<"\t"<<events_rareLepPt<<"\t"<<events_allLepPt<<endl;		
	
	p_eleLepPt->Add(p_rareLepPt); // ele 2nd
	p_jetLepPt->Add(p_eleLepPt);  // jet 3rd
	p_qcdLepPt->Add(p_jetLepPt);  // qcd 4th
	p_VGLepPt->Add(p_qcdLepPt);   // VG  5th
	for(int ibin(1); ibin < p_VGLepPt->GetSize(); ibin++){
		error_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), p_VGLepPt->GetBinContent(ibin));
		float prederror = p_VGLepPt->GetBinError(ibin);
//		prederror += p_eleLepPt->GetBinError(ibin);
//		prederror += p_jetLepPt->GetBinError(ibin);
//		prederror += p_qcdLepPt->GetBinError(ibin);
//		prederror += p_rareLepPt->GetBinError(ibin);
		error_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), 1); 
		ratioerror_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2, prederror/p_VGLepPt->GetBinContent(ibin)); 
	}
	p_VGLepPt->Draw("hist same");
	p_qcdLepPt->Draw("hist same");
	p_jetLepPt->Draw("hist same");
	p_eleLepPt->Draw("hist same");
	p_rareLepPt->Draw("hist same");
  	error_LepPt->SetFillColor(kBlack);
  	error_LepPt->SetFillStyle(3345);
	//error_LepPt->Draw("E2 same");
	TLegend *leg_leppt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_leppt->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_leppt->AddEntry(p_allLepPt,"observed (LepPt < 100 GeV)");
	leg_leppt->AddEntry(p_rareLepPt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_leppt->AddEntry(p_eleLepPt,"e->#gamma fake");
	leg_leppt->AddEntry(p_jetLepPt,"j->#gamma fake");
	leg_leppt->AddEntry(p_qcdLepPt,"j->e fake");
	leg_leppt->AddEntry(p_VGLepPt, "WG/ZG");
	//leg_leppt->Draw("same");
	leg_mt->Draw("same");
	p_allLepPt->Draw("E same");
	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( leppt_pad1,1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( leppt_pad1,2, 11 );
	else if(RunYear==2017)  		CMS_lumi( leppt_pad1,3, 11 );
	else if(RunYear==2018)  		CMS_lumi( leppt_pad1,4, 11 );

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
	ratio_leppt->Divide(p_VGLepPt);
	//cout<< p_allLepPt->Integral()/p_VGLepPt->Integral()<<endl;
	
	ratio_leppt->SetTitle("");
	ratio_leppt->GetYaxis()->SetTitle("obs./bkg.");
	ratio_leppt->GetXaxis()->SetTitle("Lepton Pt (GeV)");
	ratio_leppt->GetYaxis()->SetRangeUser(0.7,1.23);
	//ratio_leppt->GetYaxis()->SetRangeUser(0,2.1);
	ratio_leppt->GetXaxis()->SetLabelFont(63);
	ratio_leppt->GetXaxis()->SetLabelSize(14);
	ratio_leppt->GetYaxis()->SetLabelFont(63);
	ratio_leppt->GetYaxis()->SetLabelSize(14);
	ratio_leppt->Draw("ep");
	//ratioerror_LepPt->SetFillColor(kBlack);
	//ratioerror_LepPt->SetFillStyle(3345);
	//ratioerror_LepPt->Draw("E2 same");
	flatratio_leppt->Draw("same");

  	TLine *ratioValue_leppt = new TLine(25,p_allLepPt->Integral()/p_VGLepPt->Integral(),800,p_allLepPt->Integral()/p_VGLepPt->Integral());
	ratioValue_leppt->SetLineColor(kRed);
	ratioValue_leppt->Draw("same");
	c_leppt->SaveAs(ptplot.str().c_str());

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
	
	float events_eleHT, events_jetHT, events_qcdHT, events_VGHT, events_rareHT, events_allHT;
	events_eleHT =  p_eleHT->Integral();
	events_jetHT =  p_jetHT->Integral();
	events_qcdHT =  p_qcdHT->Integral();
	events_VGHT =   p_VGHT->Integral();
	events_rareHT = p_rareHT->Integral();
	events_allHT =  p_allHT->Integral();
	cout<<events_eleHT <<"\t"<<events_jetHT<<"\t"<<events_qcdHT<<"\t"<<events_VGHT<<"\t"<<events_rareHT<<"\t"<<events_allHT<<endl;		

	p_eleHT->Add(p_rareHT); // ele 2nd
	p_jetHT->Add(p_eleHT);  // jet 3rd
	p_qcdHT->Add(p_jetHT);  // qcd 4th
	p_VGHT->Add(p_qcdHT);   // VG  5th
	for(int ibin(1); ibin < p_VGHT->GetSize(); ibin++){
		error_HT->SetPoint(ibin-1,p_VGHT->GetBinCenter(ibin), p_VGHT->GetBinContent(ibin));
		float prederror = p_VGHT->GetBinError(ibin);
		error_HT->SetPointError(ibin-1,(p_VGHT->GetBinLowEdge(ibin+1)-p_VGHT->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_HT->SetPoint(ibin-1,p_VGHT->GetBinCenter(ibin), 1); 
		ratioerror_HT->SetPointError(ibin-1,(p_VGHT->GetBinLowEdge(ibin+1)-p_VGHT->GetBinLowEdge(ibin))/2, prederror/p_VGHT->GetBinContent(ibin)); 
	}
	p_VGHT->Draw("hist same");
	p_qcdHT->Draw("hist same");
	p_jetHT->Draw("hist same");
	p_eleHT->Draw("hist same");
	p_rareHT->Draw("hist same");
	leg_mt->Draw("same");
	p_allHT->Draw("E same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)  	CMS_lumi( HT_pad1,1, 11 );
	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( HT_pad1,2, 11 );
	else if(RunYear==2017)  		CMS_lumi( HT_pad1,3, 11 );
	else if(RunYear==2018)  		CMS_lumi( HT_pad1,4, 11 );

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
	ratio_HT->GetYaxis()->SetRangeUser(0.7,1.23);
	//ratio_HT->SetMinimum(0);
	//ratio_HT->SetMaximum(2);
	ratio_HT->Divide(p_VGHT);
	//cout<< p_allHT->Integral()/p_VGHT->Integral()<<endl;
	ratio_HT->SetTitle("");
	ratio_HT->GetYaxis()->SetTitle("obs./bkg.");
	ratio_HT->GetXaxis()->SetLabelFont(63);
	ratio_HT->GetXaxis()->SetLabelSize(14);
	ratio_HT->GetYaxis()->SetLabelFont(63);
	ratio_HT->GetYaxis()->SetLabelSize(14);
	ratio_HT->Draw("ep");
	//ratioerror_HT->SetFillColor(kBlack);
	//ratioerror_HT->SetFillStyle(3345);
	//ratioerror_HT->Draw("E2 same");
	flatratio_HT->Draw("same");

  	TLine *ratioValue_ht = new TLine(0,p_allHT->Integral()/p_VGHT->Integral(),1000,p_allHT->Integral()/p_VGHT->Integral());
	ratioValue_ht->SetLineColor(kRed);
	ratioValue_ht->Draw("same");
	c_HT->SaveAs(htplot.str().c_str());

	TH1F *p_allnBJet  = (TH1F*)file_sig->Get("p_nBJet");
	TH1F *p_VGnBJet   = (TH1F*)file_VG->Get("p_nBJet");
	TH1F *p_rarenBJet = (TH1F*)file_rare->Get("p_nBJet");
	TH1F *p_elenBJet  = (TH1F*)file_ele->Get("p_nBJet");
	TH1F *p_jetnBJet  = (TH1F*)file_jet->Get("p_nBJet");
	TH1F *p_qcdnBJet  = (TH1F*)file_qcd->Get("p_nBJet");
	TCanvas *can_BJet = new TCanvas("can_BJet","",600,600);
	can_BJet->cd();
	p_allnBJet->SetLineColor(1);
	p_allnBJet->SetMarkerStyle(20);
	p_allnBJet->Draw();
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
	p_VGnBJet->Draw("hist same");
	p_qcdnBJet->Draw("hist same");
	p_jetnBJet->Draw("hist same");
	p_elenBJet->Draw("hist same");
	p_rarenBJet->Draw("hist same");
	p_allnBJet->Draw("EP same");
	cout<<"\n\n\n";
	cout<<"=========================================================="<<endl;
	cout<<"\n";
	cout<< RunYear <<"  "<<whichVFP << "\tchannel : " << whichChannel <<"\t"; 
	if(plottype==1) cout<<"Control region	"<<endl;
	if(plottype==2) cout<<"Validation region	"<<endl;
	cout<<"\n\n\n";
	cout<<"MET, 		data/MC = "<< p_allMET->Integral()/p_VGMET->Integral()<<endl;
	cout<<"MT, 		data/MC = "<< p_allMt->Integral()/p_VGMt->Integral()<<endl;
	cout<<"HT, 		data/MC = "<< p_allHT->Integral()/p_VGHT->Integral()<<endl;
	cout<<"Photon Pt, 	data/MC = "<< p_allPhoEt->Integral()/p_VGPhoEt->Integral()<<endl;
	cout<<"Lepton Pt, 	data/MC = "<< p_allLepPt->Integral()/p_VGLepPt->Integral()<<endl;
	cout<<"dPhi, 		data/MC = "<< p_alldPhiEleMET->Integral()/p_VGdPhiEleMET->Integral()<<endl;
	cout<<"\n";
	cout<<"=========================================================="<<endl;

}
