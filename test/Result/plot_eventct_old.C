// Run using root -b -q "plot_eventct.C+($NBIN)"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
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
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "../../include/tdrstyle.C"
void plot_eventct(int NBIN){//main  

	gStyle->SetOptStat(0);
	setTDRStyle();
	gStyle->SetErrorX(0);
  gSystem->Load("../../lib/libAnaClasses.so");

	gStyle->SetPadTickY(0);
	gStyle->SetPadTickX(0);
	
	TFile *mgfile_ele = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_eleBkg.root");
	TFile *mgfile_jet = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_jetbkg.root");
	TFile *mgfile_qcd = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_qcd.root");
	TFile *mgfile_VG  = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_VGBkg.root");
	TFile *mgfile_rare= TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_rareBkg.root");
	TFile *egammafile_ele = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_eleBkg.root");
	TFile *egammafile_jet = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_jetbkg.root");
	TFile *egammafile_qcd = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_qcd.root");
	TFile *egammafile_VG  = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_VGBkg.root");
	TFile *egammafile_rare= TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_rareBkg.root");

	TFile *mgfile_sig = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_mg_signal.root");
	TFile *egammafile_sig = TFile::Open("/uscms_data/d3/tmishra/Output/signalTree_egamma_signal.root");
	TH1D  *mg_sig     = (TH1D*)mgfile_sig->Get("p_eventcount");
	TH1D  *egamma_sig     = (TH1D*)egammafile_sig->Get("p_eventcount");
	TH1D  *h_sig      = new TH1D("h_sig","",NBIN*2,0,NBIN*2);
	h_sig->SetBinErrorOption(TH1::kPoisson);
  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();	
	TH1D  *h_bkg      = new TH1D("h_bkg","",NBIN*2,0,NBIN*2);
	TH1D  *h_bkg_elefakepho = new TH1D("h_bkg_elefakepho","",NBIN*2,0,NBIN*2);
	TH1D  *h_bkg_jetfakepho = new TH1D("h_bkg_jetfakepho","",NBIN*2,0,NBIN*2);
	TH1D  *h_bkg_qcdfakepho = new TH1D("h_bkg_qcdfakepho","",NBIN*2,0,NBIN*2);
	TH1D  *h_bkg_VGamma     = new TH1D("h_bkg_VGamma",  "",  NBIN*2,0,NBIN*2);
	TGraphErrors *error_bkg = new TGraphErrors(NBIN*2);	
	TGraphErrors *error_ratio = new TGraphErrors(NBIN*2);	
	for(int ibin(1); ibin <= NBIN; ibin++){
		if(egamma_sig->GetBinContent(ibin) > 0.001)h_sig->SetBinContent(ibin, egamma_sig->GetBinContent(ibin));
		else h_sig->SetBinContent(ibin, 0.001);
		if(mg_sig->GetBinContent(ibin) > 0.001)h_sig->SetBinContent(ibin+NBIN, mg_sig->GetBinContent(ibin));
		else h_sig->SetBinContent(ibin+NBIN, 0.001);
	}

	for(int ibin(1); ibin <= 2*NBIN; ibin++){
		std::cout << ibin << " " << h_sig->GetBinContent(ibin) << " low " << h_sig->GetBinErrorLow(ibin) << " up " << h_sig->GetBinErrorUp(ibin) << " sqrt " << sqrt(h_sig->GetBinContent(ibin)) <<  std::endl;
  }

	TH1D *h_elefakepho_norm           = new TH1D("h_elefakepho_norm",            "h_elefakepho_norm",           NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_controlsample  = new TH1D("h_elefakepho_controlsample",   "h_elefakepho_controlsample",  NBIN*2,0,NBIN*2);
	TH1D *h_elefakepho_transferfactor = new TH1D("h_elefakepho_transferfactor",  "h_elefakepho_transferfactor", NBIN*2,0,NBIN*2);
	TH1D *h_elefakepho_syserr_jes     = new TH1D("h_elefakepho_syserr_jes",      "h_elefakepho_syserr_jes",     NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_jer     = new TH1D("h_elefakepho_syserr_jer",      "h_elefakepho_syserr_jer",     NBIN*2,0,NBIN*2);   
	TH1D *h_elefakepho_syserr_esf     = new TH1D("h_elefakepho_syserr_esf",      "h_elefakepho_syserr_esf",     NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_scale   = new TH1D("h_elefakepho_syserr_scale",    "h_elefakepho_syserr_scale",   NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_eleshape= new TH1D("h_elefakepho_syserr_e_to_pho_syst", "h_elefakepho_syserr_eleshape",NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_jetshape= new TH1D("h_elefakepho_syserr_j_to_pho_syst", "h_elefakepho_syserr_jetshape",NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_qcdshape= new TH1D("h_elefakepho_syserr_fakelep_shape", "h_elefakepho_syserr_qcdshape",NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_xs      = new TH1D("h_elefakepho_syserr_xs",       "h_elefakepho_syserr_xs",      NBIN*2,0,NBIN*2);
	TH1D *h_elefakepho_syserr_lumi    = new TH1D("h_elefakepho_syserr_lumi",     "h_elefakepho_syserr_lumi",    NBIN*2,0,NBIN*2); 
	TH1D *h_elefakepho_syserr_isr     = new TH1D("h_elefakepho_syserr_isr",      "h_elefakepho_syserr_isr",     NBIN*2,0,NBIN*2); 
	                                                                                                          
	TH1D *h_jetfakepho_norm           = new TH1D("h_jetfakepho_norm",            "h_jetfakepho_norm",           NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_controlsample  = new TH1D("h_jetfakepho_controlsample",   "h_jetfakepho_controlsample",  NBIN*2,0,NBIN*2);
	TH1D *h_jetfakepho_transferfactor = new TH1D("h_jetfakepho_transferfactor",  "h_jetfakepho_transferfactor", NBIN*2,0,NBIN*2);
	TH1D *h_jetfakepho_syserr_jes     = new TH1D("h_jetfakepho_syserr_jes",      "h_jetfakepho_syserr_jes",     NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_jer     = new TH1D("h_jetfakepho_syserr_jer",      "h_jetfakepho_syserr_jer",     NBIN*2,0,NBIN*2);   
	TH1D *h_jetfakepho_syserr_esf     = new TH1D("h_jetfakepho_syserr_esf",      "h_jetfakepho_syserr_esf",     NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_scale   = new TH1D("h_jetfakepho_syserr_scale",    "h_jetfakepho_syserr_scale",   NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_eleshape= new TH1D("h_jetfakepho_syserr_e_to_pho_syst", "h_jetfakepho_syserr_eleshape",NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_jetshape= new TH1D("h_jetfakepho_syserr_j_to_pho_syst", "h_jetfakepho_syserr_jetshape",NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_qcdshape= new TH1D("h_jetfakepho_syserr_fakelep_shape", "h_jetfakepho_syserr_qcdshape",NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_xs      = new TH1D("h_jetfakepho_syserr_xs",       "h_jetfakepho_syserr_xs",      NBIN*2,0,NBIN*2);
	TH1D *h_jetfakepho_syserr_lumi    = new TH1D("h_jetfakepho_syserr_lumi",     "h_jetfakepho_syserr_lumi",    NBIN*2,0,NBIN*2); 
	TH1D *h_jetfakepho_syserr_isr     = new TH1D("h_jetfakepho_syserr_isr",      "h_jetfakepho_syserr_isr",     NBIN*2,0,NBIN*2); 
	                                                                                                          
	                                                                             
	TH1D *h_qcdfakelep_norm           = new TH1D("h_qcdfakelep_norm",            "h_qcdfakelep_norm",           NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_controlsample  = new TH1D("h_qcdfakelep_controlsample",   "h_qcdfakelep_controlsample",  NBIN*2,0,NBIN*2);
	TH1D *h_qcdfakelep_transferfactor = new TH1D("h_qcdfakelep_transferfactor",  "h_qcdfakelep_transferfactor", NBIN*2,0,NBIN*2);
	TH1D *h_qcdfakelep_syserr_jes     = new TH1D("h_qcdfakelep_syserr_jes",      "h_qcdfakelep_syserr_jes",     NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_jer     = new TH1D("h_qcdfakelep_syserr_jer",      "h_qcdfakelep_syserr_jer",     NBIN*2,0,NBIN*2);   
	TH1D *h_qcdfakelep_syserr_esf     = new TH1D("h_qcdfakelep_syserr_esf",      "h_qcdfakelep_syserr_esf",     NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_scale   = new TH1D("h_qcdfakelep_syserr_scale",    "h_qcdfakelep_syserr_scale",   NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_eleshape= new TH1D("h_qcdfakelep_syserr_e_to_pho_syst", "h_qcdfakelep_syserr_eleshape",NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_jetshape= new TH1D("h_qcdfakelep_syserr_j_to_pho_syst", "h_qcdfakelep_syserr_jetshape",NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_qcdshape= new TH1D("h_qcdfakelep_syserr_fakelep_shape", "h_qcdfakelep_syserr_qcdshape",NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_xs      = new TH1D("h_qcdfakelep_syserr_xs",       "h_qcdfakelep_syserr_xs",      NBIN*2,0,NBIN*2);
	TH1D *h_qcdfakelep_syserr_lumi    = new TH1D("h_qcdfakelep_syserr_lumi",     "h_qcdfakelep_syserr_lumi",    NBIN*2,0,NBIN*2); 
	TH1D *h_qcdfakelep_syserr_isr     = new TH1D("h_qcdfakelep_syserr_isr",      "h_qcdfakelep_syserr_isr",     NBIN*2,0,NBIN*2); 
	
	TH1D *h_VGamma_norm           = new TH1D("h_VGamma_norm",             "h_VGamma_norm",           NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_jes     = new TH1D("h_VGamma_syserr_jes",       "h_VGamma_syserr_jes",     NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_jer     = new TH1D("h_VGamma_syserr_jer",       "h_VGamma_syserr_jer",     NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_esf     = new TH1D("h_VGamma_syserr_esf",       "h_VGamma_syserr_esf",     NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_scale   = new TH1D("h_VGamma_syserr_scale",     "h_VGamma_syserr_scale",   NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_eleshape= new TH1D("h_VGamma_syserr_e_to_pho_syst",  "h_VGamma_syserr_eleshape",NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_jetshape= new TH1D("h_VGamma_syserr_j_to_pho_syst",  "h_VGamma_syserr_jetshape",NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_qcdshape= new TH1D("h_VGamma_syserr_fakelep_shape",  "h_VGamma_syserr_qcdshape",NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_xs      = new TH1D("h_VGamma_syserr_xs",        "h_VGamma_syserr_xs",      NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_lumi    = new TH1D("h_VGamma_syserr_lumi",      "h_VGamma_syserr_lumi",    NBIN*2,0,NBIN*2);   
	TH1D *h_VGamma_syserr_isr     = new TH1D("h_VGamma_syserr_isr",       "h_VGamma_syserr_isr",     NBIN*2,0,NBIN*2);   
	
	TH1D *h_rare_norm           = new TH1D("h_rare_norm",              "h_rare_norm",           NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_jes     = new TH1D("h_rare_syserr_jes",        "h_rare_syserr_jes",     NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_jer     = new TH1D("h_rare_syserr_jer",        "h_rare_syserr_jer",     NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_esf     = new TH1D("h_rare_syserr_esf",        "h_rare_syserr_esf",     NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_scale   = new TH1D("h_rare_syserr_scale",      "h_rare_syserr_scale",   NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_eleshape= new TH1D("h_rare_syserr_e_to_pho_syst",   "h_rare_syserr_eleshape",NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_jetshape= new TH1D("h_rare_syserr_j_to_pho_syst",   "h_rare_syserr_jetshape",NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_qcdshape= new TH1D("h_rare_syserr_fakelep_shape",   "h_rare_syserr_qcdshape",NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_xs      = new TH1D("h_rare_syserr_xs",         "h_rare_syserr_xs",      NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_lumi    = new TH1D("h_rare_syserr_lumi",       "h_rare_syserr_lumi",    NBIN*2,0,NBIN*2);   
	TH1D *h_rare_syserr_isr     = new TH1D("h_rare_syserr_isr",        "h_rare_syserr_isr",     NBIN*2,0,NBIN*2);   



	TH1D *mg_elefakepho_norm           = (TH1D*)mgfile_ele->Get("mg_elefakepho_norm");         
	TH1D *mg_elefakepho_controlsample  = (TH1D*)mgfile_ele->Get("mg_elefakepho_controlsample");
	TH1D *mg_elefakepho_transferfactor = (TH1D*)mgfile_ele->Get("mg_elefakepho_transferfactor");
	TH1D *mg_elefakepho_syserr_jes     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jes");     
	TH1D *mg_elefakepho_syserr_jer     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jer");  
	TH1D *mg_elefakepho_syserr_esf     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_esf");  
	TH1D *mg_elefakepho_syserr_scale   = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_scale");  
	TH1D *mg_elefakepho_syserr_eleshape= (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_e_to_pho");
	TH1D *mg_elefakepho_syserr_jetshape= (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_j_to_pho");
	TH1D *mg_elefakepho_syserr_qcdshape= (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_j_to_lep");
	TH1D *mg_elefakepho_syserr_xs      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_xs");
	TH1D *mg_elefakepho_syserr_lumi    = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_lumi");    
	TH1D *mg_elefakepho_syserr_isr     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_isr");    
	
	TH1D *mg_jetfakepho_norm           = (TH1D*)mgfile_jet->Get("mg_jetfakepho_norm");         
	TH1D *mg_jetfakepho_controlsample  = (TH1D*)mgfile_jet->Get("mg_jetfakepho_controlsample");
	TH1D *mg_jetfakepho_transferfactor = (TH1D*)mgfile_jet->Get("mg_jetfakepho_transferfactor");
	TH1D *mg_jetfakepho_syserr_jes     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jes");     
	TH1D *mg_jetfakepho_syserr_jer     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jer");  
	TH1D *mg_jetfakepho_syserr_esf     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_esf");  
	TH1D *mg_jetfakepho_syserr_scale   = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_scale");  
	TH1D *mg_jetfakepho_syserr_eleshape= (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_e_to_pho");
	TH1D *mg_jetfakepho_syserr_jetshape= (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_j_to_pho");
	TH1D *mg_jetfakepho_syserr_qcdshape= (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_j_to_lep");
	TH1D *mg_jetfakepho_syserr_xs      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_xs");
	TH1D *mg_jetfakepho_syserr_lumi    = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_lumi");    
	TH1D *mg_jetfakepho_syserr_isr     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_isr");    
	
	TH1D *mg_qcdfakelep_norm           = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_norm");         
	TH1D *mg_qcdfakelep_controlsample  = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_controlsample");
	TH1D *mg_qcdfakelep_transferfactor = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_transferfactor");
	TH1D *mg_qcdfakelep_syserr_jes     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jes");     
	TH1D *mg_qcdfakelep_syserr_jer     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jer");  
	TH1D *mg_qcdfakelep_syserr_esf     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_esf");  
	TH1D *mg_qcdfakelep_syserr_scale   = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_scale");  
	TH1D *mg_qcdfakelep_syserr_eleshape= (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_e_to_pho");
	TH1D *mg_qcdfakelep_syserr_jetshape= (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_j_to_pho");
	TH1D *mg_qcdfakelep_syserr_qcdshape= (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_j_to_lep");
	TH1D *mg_qcdfakelep_syserr_xs      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_xs");
	TH1D *mg_qcdfakelep_syserr_lumi    = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_lumi");    
	TH1D *mg_qcdfakelep_syserr_isr     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_isr");    
	
	TH1D *mg_VGamma_norm           = (TH1D*)mgfile_VG->Get("mg_VGamma_norm");         
	TH1D *mg_VGamma_syserr_jes     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jes");     
	TH1D *mg_VGamma_syserr_jer     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jer");  
	TH1D *mg_VGamma_syserr_esf     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_esf");  
	TH1D *mg_VGamma_syserr_scale   = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_scale");  
	TH1D *mg_VGamma_syserr_eleshape= (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_eleshape");
	TH1D *mg_VGamma_syserr_jetshape= (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jetshape");
	TH1D *mg_VGamma_syserr_qcdshape= (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_qcdshape");
	TH1D *mg_VGamma_syserr_xs      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_xs");
	TH1D *mg_VGamma_syserr_lumi    = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_lumi");    
	TH1D *mg_VGamma_syserr_isr     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_isr");    
	
	TH1D *mg_rare_norm           = (TH1D*)mgfile_rare->Get("mg_rare_norm");         
	TH1D *mg_rare_syserr_jes     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jes");     
	TH1D *mg_rare_syserr_jer     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jer");  
	TH1D *mg_rare_syserr_esf     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_esf");  
	TH1D *mg_rare_syserr_scale   = (TH1D*)mgfile_rare->Get("mg_rare_syserr_scale");  
	TH1D *mg_rare_syserr_eleshape= (TH1D*)mgfile_rare->Get("mg_rare_syserr_eleshape");
	TH1D *mg_rare_syserr_jetshape= (TH1D*)mgfile_rare->Get("mg_rare_syserr_jetshape");
	TH1D *mg_rare_syserr_qcdshape= (TH1D*)mgfile_rare->Get("mg_rare_syserr_qcdshape");
	TH1D *mg_rare_syserr_xs      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_xs");
	TH1D *mg_rare_syserr_lumi    = (TH1D*)mgfile_rare->Get("mg_rare_syserr_lumi");    
	TH1D *mg_rare_syserr_isr     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_isr");    


	TH1D *egamma_elefakepho_norm           = (TH1D*)egammafile_ele->Get("eg_elefakepho_norm");         
	TH1D *egamma_elefakepho_controlsample  = (TH1D*)egammafile_ele->Get("eg_elefakepho_controlsample");
	TH1D *egamma_elefakepho_transferfactor = (TH1D*)egammafile_ele->Get("eg_elefakepho_transferfactor");
	TH1D *egamma_elefakepho_syserr_jes     = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_jes");     
	TH1D *egamma_elefakepho_syserr_jer     = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_jer");  
	TH1D *egamma_elefakepho_syserr_esf     = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_esf");  
	TH1D *egamma_elefakepho_syserr_scale   = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_scale");  
	TH1D *egamma_elefakepho_syserr_eleshape= (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_e_to_pho");
	TH1D *egamma_elefakepho_syserr_jetshape= (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_j_to_pho");
	TH1D *egamma_elefakepho_syserr_qcdshape= (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_j_to_lep");
	TH1D *egamma_elefakepho_syserr_xs      = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_xs");
	TH1D *egamma_elefakepho_syserr_lumi    = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_lumi");    
	TH1D *egamma_elefakepho_syserr_isr     = (TH1D*)egammafile_ele->Get("eg_elefakepho_syserr_isr");    
	
	TH1D *egamma_jetfakepho_norm           = (TH1D*)egammafile_jet->Get("eg_jetfakepho_norm");         
	TH1D *egamma_jetfakepho_controlsample  = (TH1D*)egammafile_jet->Get("eg_jetfakepho_controlsample");
	TH1D *egamma_jetfakepho_transferfactor = (TH1D*)egammafile_jet->Get("eg_jetfakepho_transferfactor");
	TH1D *egamma_jetfakepho_syserr_jes     = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_jes");     
	TH1D *egamma_jetfakepho_syserr_jer     = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_jer");  
	TH1D *egamma_jetfakepho_syserr_esf     = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_esf");  
	TH1D *egamma_jetfakepho_syserr_scale   = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_scale");  
	TH1D *egamma_jetfakepho_syserr_eleshape= (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_e_to_pho");
	TH1D *egamma_jetfakepho_syserr_jetshape= (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_j_to_pho");
	TH1D *egamma_jetfakepho_syserr_qcdshape= (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_j_to_lep");
	TH1D *egamma_jetfakepho_syserr_xs      = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_xs");
	TH1D *egamma_jetfakepho_syserr_lumi    = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_lumi");    
	TH1D *egamma_jetfakepho_syserr_isr     = (TH1D*)egammafile_jet->Get("eg_jetfakepho_syserr_isr");    
	
	TH1D *egamma_qcdfakelep_norm           = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_norm");         
	TH1D *egamma_qcdfakelep_controlsample  = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_controlsample");
	TH1D *egamma_qcdfakelep_transferfactor = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_transferfactor");
	TH1D *egamma_qcdfakelep_syserr_jes     = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_jes");     
	TH1D *egamma_qcdfakelep_syserr_jer     = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_jer");  
	TH1D *egamma_qcdfakelep_syserr_esf     = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_esf");  
	TH1D *egamma_qcdfakelep_syserr_scale   = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_scale");  
	TH1D *egamma_qcdfakelep_syserr_eleshape= (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_e_to_pho");
	TH1D *egamma_qcdfakelep_syserr_jetshape= (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_j_to_pho");
	TH1D *egamma_qcdfakelep_syserr_qcdshape= (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_j_to_lep");
	TH1D *egamma_qcdfakelep_syserr_xs      = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_xs");
	TH1D *egamma_qcdfakelep_syserr_lumi    = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_lumi");    
	TH1D *egamma_qcdfakelep_syserr_isr     = (TH1D*)egammafile_qcd->Get("eg_qcdfakelep_syserr_isr");    
	
	TH1D *egamma_VGamma_norm           = (TH1D*)egammafile_VG->Get("eg_VGamma_norm");         
	TH1D *egamma_VGamma_syserr_jes     = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_jes");     
	TH1D *egamma_VGamma_syserr_jer     = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_jer");  
	TH1D *egamma_VGamma_syserr_esf     = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_esf");  
	TH1D *egamma_VGamma_syserr_scale   = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_scale");  
	TH1D *egamma_VGamma_syserr_eleshape= (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_eleshape");
	TH1D *egamma_VGamma_syserr_jetshape= (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_jetshape");
	TH1D *egamma_VGamma_syserr_qcdshape= (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_qcdshape");
	TH1D *egamma_VGamma_syserr_xs      = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_xs");
	TH1D *egamma_VGamma_syserr_lumi    = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_lumi");    
	TH1D *egamma_VGamma_syserr_isr     = (TH1D*)egammafile_VG->Get("eg_VGamma_syserr_isr");    
	
	TH1D *egamma_rare_norm           = (TH1D*)egammafile_rare->Get("eg_rare_norm");         
	TH1D *egamma_rare_syserr_jes     = (TH1D*)egammafile_rare->Get("eg_rare_syserr_jes");     
	TH1D *egamma_rare_syserr_jer     = (TH1D*)egammafile_rare->Get("eg_rare_syserr_jer");  
	TH1D *egamma_rare_syserr_esf     = (TH1D*)egammafile_rare->Get("eg_rare_syserr_esf");  
	TH1D *egamma_rare_syserr_scale   = (TH1D*)egammafile_rare->Get("eg_rare_syserr_scale");  
	TH1D *egamma_rare_syserr_eleshape= (TH1D*)egammafile_rare->Get("eg_rare_syserr_eleshape");
	TH1D *egamma_rare_syserr_jetshape= (TH1D*)egammafile_rare->Get("eg_rare_syserr_jetshape");
	TH1D *egamma_rare_syserr_qcdshape= (TH1D*)egammafile_rare->Get("eg_rare_syserr_qcdshape");
	TH1D *egamma_rare_syserr_xs      = (TH1D*)egammafile_rare->Get("eg_rare_syserr_xs");
	TH1D *egamma_rare_syserr_lumi    = (TH1D*)egammafile_rare->Get("eg_rare_syserr_lumi");    
	TH1D *egamma_rare_syserr_isr     = (TH1D*)egammafile_rare->Get("eg_rare_syserr_isr");    

	double bkgContent(0);
	double bkgError(0);
	
	for(int ibin(1); ibin <= NBIN; ibin++){
		h_elefakepho_norm->SetBinContent(ibin, egamma_elefakepho_norm->GetBinContent(ibin));
		h_elefakepho_controlsample->SetBinContent(ibin,   egamma_elefakepho_controlsample->GetBinContent(ibin));
		h_elefakepho_transferfactor->SetBinContent(ibin,  egamma_elefakepho_transferfactor->GetBinContent(ibin));
		h_elefakepho_syserr_jes->SetBinContent(ibin,      egamma_elefakepho_syserr_jes->GetBinContent(ibin));      
		h_elefakepho_syserr_jer->SetBinContent(ibin,      egamma_elefakepho_syserr_jer->GetBinContent(ibin)); 
		h_elefakepho_syserr_esf->SetBinContent(ibin,      egamma_elefakepho_syserr_esf->GetBinContent(ibin));      
		h_elefakepho_syserr_scale->SetBinContent(ibin,    egamma_elefakepho_syserr_scale->GetBinContent(ibin));    
		h_elefakepho_syserr_eleshape->SetBinContent(ibin, egamma_elefakepho_syserr_eleshape->GetBinContent(ibin));
		h_elefakepho_syserr_jetshape->SetBinContent(ibin, egamma_elefakepho_syserr_jetshape->GetBinContent(ibin));
		h_elefakepho_syserr_qcdshape->SetBinContent(ibin, egamma_elefakepho_syserr_qcdshape->GetBinContent(ibin));
		h_elefakepho_syserr_xs->SetBinContent(ibin,       egamma_elefakepho_syserr_xs->GetBinContent(ibin));       
		h_elefakepho_syserr_lumi->SetBinContent(ibin,     egamma_elefakepho_syserr_lumi->GetBinContent(ibin));     
		h_elefakepho_syserr_isr->SetBinContent(ibin,      egamma_elefakepho_syserr_isr->GetBinContent(ibin));     
																																																																																																							
		h_jetfakepho_norm->SetBinContent(ibin,            egamma_jetfakepho_norm->GetBinContent(ibin));            
		h_jetfakepho_controlsample->SetBinContent(ibin,   egamma_jetfakepho_controlsample->GetBinContent(ibin));
		h_jetfakepho_transferfactor->SetBinContent(ibin,  egamma_jetfakepho_transferfactor->GetBinContent(ibin));
		h_jetfakepho_syserr_jes->SetBinContent(ibin,      egamma_jetfakepho_syserr_jes->GetBinContent(ibin));      
		h_jetfakepho_syserr_jer->SetBinContent(ibin,      egamma_jetfakepho_syserr_jer->GetBinContent(ibin));      
		h_jetfakepho_syserr_esf->SetBinContent(ibin,      egamma_jetfakepho_syserr_esf->GetBinContent(ibin));      
		h_jetfakepho_syserr_scale->SetBinContent(ibin,    egamma_jetfakepho_syserr_scale->GetBinContent(ibin));    
		h_jetfakepho_syserr_eleshape->SetBinContent(ibin, egamma_jetfakepho_syserr_eleshape->GetBinContent(ibin));
		h_jetfakepho_syserr_jetshape->SetBinContent(ibin, egamma_jetfakepho_syserr_jetshape->GetBinContent(ibin));
		h_jetfakepho_syserr_qcdshape->SetBinContent(ibin, egamma_jetfakepho_syserr_qcdshape->GetBinContent(ibin));
		h_jetfakepho_syserr_xs->SetBinContent(ibin,       egamma_jetfakepho_syserr_xs->GetBinContent(ibin));       
		h_jetfakepho_syserr_lumi->SetBinContent(ibin,     egamma_jetfakepho_syserr_lumi->GetBinContent(ibin));     
		h_jetfakepho_syserr_isr->SetBinContent(ibin,     egamma_jetfakepho_syserr_isr->GetBinContent(ibin));     
																														
		h_qcdfakelep_norm->SetBinContent(ibin,            egamma_qcdfakelep_norm->GetBinContent(ibin));            
		h_qcdfakelep_controlsample->SetBinContent(ibin,   egamma_qcdfakelep_controlsample->GetBinContent(ibin));
		h_qcdfakelep_transferfactor->SetBinContent(ibin,  egamma_qcdfakelep_transferfactor->GetBinContent(ibin));
		h_qcdfakelep_syserr_jes->SetBinContent(ibin,      egamma_qcdfakelep_syserr_jes->GetBinContent(ibin));      
		h_qcdfakelep_syserr_jer->SetBinContent(ibin,      egamma_qcdfakelep_syserr_jer->GetBinContent(ibin));      
		h_qcdfakelep_syserr_esf->SetBinContent(ibin,      egamma_qcdfakelep_syserr_esf->GetBinContent(ibin));      
		h_qcdfakelep_syserr_scale->SetBinContent(ibin,    egamma_qcdfakelep_syserr_scale->GetBinContent(ibin));    
		h_qcdfakelep_syserr_eleshape->SetBinContent(ibin, egamma_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_jetshape->SetBinContent(ibin, egamma_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_qcdshape->SetBinContent(ibin, egamma_qcdfakelep_syserr_qcdshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_xs->SetBinContent(ibin,       egamma_qcdfakelep_syserr_xs->GetBinContent(ibin));       
		h_qcdfakelep_syserr_lumi->SetBinContent(ibin,     egamma_qcdfakelep_syserr_lumi->GetBinContent(ibin));     
		h_qcdfakelep_syserr_isr->SetBinContent(ibin,     egamma_qcdfakelep_syserr_isr->GetBinContent(ibin));     
																
		h_VGamma_norm->SetBinContent(ibin,                egamma_VGamma_norm->GetBinContent(ibin));            
		h_VGamma_syserr_jes->SetBinContent(ibin,          egamma_VGamma_syserr_jes->GetBinContent(ibin));      
		h_VGamma_syserr_jer->SetBinContent(ibin,          egamma_VGamma_syserr_jer->GetBinContent(ibin));      
		h_VGamma_syserr_esf->SetBinContent(ibin,          egamma_VGamma_syserr_esf->GetBinContent(ibin));      
		h_VGamma_syserr_scale->SetBinContent(ibin,        egamma_VGamma_syserr_scale->GetBinContent(ibin));    
		h_VGamma_syserr_eleshape->SetBinContent(ibin,     egamma_VGamma_syserr_eleshape->GetBinContent(ibin));
		h_VGamma_syserr_jetshape->SetBinContent(ibin,     egamma_VGamma_syserr_jetshape->GetBinContent(ibin));
		h_VGamma_syserr_qcdshape->SetBinContent(ibin,     egamma_VGamma_syserr_qcdshape->GetBinContent(ibin));
		h_VGamma_syserr_xs->SetBinContent(ibin,           egamma_VGamma_syserr_xs->GetBinContent(ibin));       
		h_VGamma_syserr_lumi->SetBinContent(ibin,         egamma_VGamma_syserr_lumi->GetBinContent(ibin));     
		h_VGamma_syserr_isr->SetBinContent(ibin,          egamma_VGamma_syserr_isr->GetBinContent(ibin));     
																
		h_rare_norm->SetBinContent(ibin,                  egamma_rare_norm->GetBinContent(ibin));            
		h_rare_syserr_jes->SetBinContent(ibin,            egamma_rare_syserr_jes->GetBinContent(ibin));      
		h_rare_syserr_jer->SetBinContent(ibin,            egamma_rare_syserr_jer->GetBinContent(ibin));      
		h_rare_syserr_esf->SetBinContent(ibin,            egamma_rare_syserr_esf->GetBinContent(ibin));      
		h_rare_syserr_scale->SetBinContent(ibin,          egamma_rare_syserr_scale->GetBinContent(ibin));   
		h_rare_syserr_eleshape->SetBinContent(ibin,       egamma_rare_syserr_eleshape->GetBinContent(ibin));   
		h_rare_syserr_jetshape->SetBinContent(ibin,       egamma_rare_syserr_jetshape->GetBinContent(ibin)); 
		h_rare_syserr_qcdshape->SetBinContent(ibin,       egamma_rare_syserr_qcdshape->GetBinContent(ibin)); 
		h_rare_syserr_xs->SetBinContent(ibin,             egamma_rare_syserr_xs->GetBinContent(ibin));       
		h_rare_syserr_lumi->SetBinContent(ibin,           egamma_rare_syserr_lumi->GetBinContent(ibin));     
		h_rare_syserr_isr->SetBinContent(ibin,           egamma_rare_syserr_isr->GetBinContent(ibin));     


		h_elefakepho_norm->SetBinContent(ibin+NBIN, mg_elefakepho_norm->GetBinContent(ibin));
		h_elefakepho_controlsample->SetBinContent(ibin+NBIN,   mg_elefakepho_controlsample->GetBinContent(ibin));
		h_elefakepho_transferfactor->SetBinContent(ibin+NBIN,  mg_elefakepho_transferfactor->GetBinContent(ibin));
		h_elefakepho_syserr_jes->SetBinContent(ibin+NBIN,      mg_elefakepho_syserr_jes->GetBinContent(ibin));      
		h_elefakepho_syserr_jer->SetBinContent(ibin+NBIN,      mg_elefakepho_syserr_jer->GetBinContent(ibin)); 
		h_elefakepho_syserr_esf->SetBinContent(ibin+NBIN,      mg_elefakepho_syserr_esf->GetBinContent(ibin));      
		h_elefakepho_syserr_scale->SetBinContent(ibin+NBIN,    mg_elefakepho_syserr_scale->GetBinContent(ibin));    
		h_elefakepho_syserr_eleshape->SetBinContent(ibin+NBIN, mg_elefakepho_syserr_eleshape->GetBinContent(ibin));
		h_elefakepho_syserr_jetshape->SetBinContent(ibin+NBIN, mg_elefakepho_syserr_jetshape->GetBinContent(ibin));
		h_elefakepho_syserr_qcdshape->SetBinContent(ibin+NBIN, mg_elefakepho_syserr_qcdshape->GetBinContent(ibin));
		h_elefakepho_syserr_xs->SetBinContent(ibin+NBIN,       mg_elefakepho_syserr_xs->GetBinContent(ibin));       
		h_elefakepho_syserr_lumi->SetBinContent(ibin+NBIN,     mg_elefakepho_syserr_lumi->GetBinContent(ibin));     
		h_elefakepho_syserr_isr->SetBinContent(ibin+NBIN,     mg_elefakepho_syserr_isr->GetBinContent(ibin));     
																																																																																																							
		h_jetfakepho_norm->SetBinContent(ibin+NBIN,            mg_jetfakepho_norm->GetBinContent(ibin));            
		h_jetfakepho_controlsample->SetBinContent(ibin+NBIN,   mg_jetfakepho_controlsample->GetBinContent(ibin));
		h_jetfakepho_transferfactor->SetBinContent(ibin+NBIN,  mg_jetfakepho_transferfactor->GetBinContent(ibin));
		h_jetfakepho_syserr_jes->SetBinContent(ibin+NBIN,      mg_jetfakepho_syserr_jes->GetBinContent(ibin));      
		h_jetfakepho_syserr_jer->SetBinContent(ibin+NBIN,      mg_jetfakepho_syserr_jer->GetBinContent(ibin));      
		h_jetfakepho_syserr_esf->SetBinContent(ibin+NBIN,      mg_jetfakepho_syserr_esf->GetBinContent(ibin));      
		h_jetfakepho_syserr_scale->SetBinContent(ibin+NBIN,    mg_jetfakepho_syserr_scale->GetBinContent(ibin));    
		h_jetfakepho_syserr_eleshape->SetBinContent(ibin+NBIN, mg_jetfakepho_syserr_eleshape->GetBinContent(ibin));
		h_jetfakepho_syserr_jetshape->SetBinContent(ibin+NBIN, mg_jetfakepho_syserr_jetshape->GetBinContent(ibin));
		h_jetfakepho_syserr_qcdshape->SetBinContent(ibin+NBIN, mg_jetfakepho_syserr_qcdshape->GetBinContent(ibin));
		h_jetfakepho_syserr_xs->SetBinContent(ibin+NBIN,       mg_jetfakepho_syserr_xs->GetBinContent(ibin));       
		h_jetfakepho_syserr_lumi->SetBinContent(ibin+NBIN,     mg_jetfakepho_syserr_lumi->GetBinContent(ibin));     
		h_jetfakepho_syserr_isr->SetBinContent(ibin+NBIN,     mg_jetfakepho_syserr_isr->GetBinContent(ibin));     
																														
		h_qcdfakelep_norm->SetBinContent(ibin+NBIN,            mg_qcdfakelep_norm->GetBinContent(ibin));            
		h_qcdfakelep_controlsample->SetBinContent(ibin+NBIN,   mg_qcdfakelep_controlsample->GetBinContent(ibin));
		h_qcdfakelep_transferfactor->SetBinContent(ibin+NBIN,  mg_qcdfakelep_transferfactor->GetBinContent(ibin));
		h_qcdfakelep_syserr_jes->SetBinContent(ibin+NBIN,      mg_qcdfakelep_syserr_jes->GetBinContent(ibin));      
		h_qcdfakelep_syserr_jer->SetBinContent(ibin+NBIN,      mg_qcdfakelep_syserr_jer->GetBinContent(ibin));      
		h_qcdfakelep_syserr_esf->SetBinContent(ibin+NBIN,      mg_qcdfakelep_syserr_esf->GetBinContent(ibin));      
		h_qcdfakelep_syserr_scale->SetBinContent(ibin+NBIN,    mg_qcdfakelep_syserr_scale->GetBinContent(ibin));    
		h_qcdfakelep_syserr_eleshape->SetBinContent(ibin+NBIN, mg_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_jetshape->SetBinContent(ibin+NBIN, mg_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_qcdshape->SetBinContent(ibin+NBIN, mg_qcdfakelep_syserr_qcdshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_xs->SetBinContent(ibin+NBIN,       mg_qcdfakelep_syserr_xs->GetBinContent(ibin));       
		h_qcdfakelep_syserr_lumi->SetBinContent(ibin+NBIN,     mg_qcdfakelep_syserr_lumi->GetBinContent(ibin));     
		h_qcdfakelep_syserr_isr->SetBinContent(ibin+NBIN,     mg_qcdfakelep_syserr_isr->GetBinContent(ibin));     
																
		h_VGamma_norm->SetBinContent(ibin+NBIN,                mg_VGamma_norm->GetBinContent(ibin));            
		h_VGamma_syserr_jes->SetBinContent(ibin+NBIN,          mg_VGamma_syserr_jes->GetBinContent(ibin));      
		h_VGamma_syserr_jer->SetBinContent(ibin+NBIN,          mg_VGamma_syserr_jer->GetBinContent(ibin));      
		h_VGamma_syserr_esf->SetBinContent(ibin+NBIN,          mg_VGamma_syserr_esf->GetBinContent(ibin));      
		h_VGamma_syserr_scale->SetBinContent(ibin+NBIN,        mg_VGamma_syserr_scale->GetBinContent(ibin));    
		h_VGamma_syserr_eleshape->SetBinContent(ibin+NBIN,     mg_VGamma_syserr_eleshape->GetBinContent(ibin));
		h_VGamma_syserr_jetshape->SetBinContent(ibin+NBIN,     mg_VGamma_syserr_jetshape->GetBinContent(ibin));
		h_VGamma_syserr_qcdshape->SetBinContent(ibin+NBIN,     mg_VGamma_syserr_qcdshape->GetBinContent(ibin));
		h_VGamma_syserr_xs->SetBinContent(ibin+NBIN,           mg_VGamma_syserr_xs->GetBinContent(ibin));       
		h_VGamma_syserr_lumi->SetBinContent(ibin+NBIN,         mg_VGamma_syserr_lumi->GetBinContent(ibin));     
		h_VGamma_syserr_isr->SetBinContent(ibin+NBIN,          mg_VGamma_syserr_isr->GetBinContent(ibin));     
																
		h_rare_norm->SetBinContent(ibin+NBIN,                  mg_rare_norm->GetBinContent(ibin));            
		h_rare_syserr_jes->SetBinContent(ibin+NBIN,            mg_rare_syserr_jes->GetBinContent(ibin));      
		h_rare_syserr_jer->SetBinContent(ibin+NBIN,            mg_rare_syserr_jer->GetBinContent(ibin));      
		h_rare_syserr_esf->SetBinContent(ibin+NBIN,            mg_rare_syserr_esf->GetBinContent(ibin));      
		h_rare_syserr_scale->SetBinContent(ibin+NBIN,          mg_rare_syserr_scale->GetBinContent(ibin));   
		h_rare_syserr_eleshape->SetBinContent(ibin+NBIN,       mg_rare_syserr_eleshape->GetBinContent(ibin));   
		h_rare_syserr_jetshape->SetBinContent(ibin+NBIN,       mg_rare_syserr_jetshape->GetBinContent(ibin)); 
		h_rare_syserr_qcdshape->SetBinContent(ibin+NBIN,       mg_rare_syserr_qcdshape->GetBinContent(ibin)); 
		h_rare_syserr_xs->SetBinContent(ibin+NBIN,             mg_rare_syserr_xs->GetBinContent(ibin));       
		h_rare_syserr_lumi->SetBinContent(ibin+NBIN,           mg_rare_syserr_lumi->GetBinContent(ibin));     
		h_rare_syserr_isr->SetBinContent(ibin+NBIN,           mg_rare_syserr_isr->GetBinContent(ibin));     

		h_elefakepho_norm->SetBinError(ibin,            egamma_elefakepho_norm->GetBinError(ibin));
		h_jetfakepho_norm->SetBinError(ibin,            egamma_jetfakepho_norm->GetBinError(ibin));            
		h_qcdfakelep_norm->SetBinError(ibin,            egamma_qcdfakelep_norm->GetBinError(ibin));            
		h_VGamma_norm->SetBinError(ibin,                egamma_VGamma_norm->GetBinError(ibin));            
		h_rare_norm->SetBinError(ibin,                  egamma_rare_norm->GetBinError(ibin));            
		h_elefakepho_transferfactor->SetBinError(ibin,  egamma_elefakepho_transferfactor->GetBinError(ibin));   
		h_jetfakepho_transferfactor->SetBinError(ibin,  egamma_jetfakepho_transferfactor->GetBinError(ibin));   
		h_qcdfakelep_transferfactor->SetBinError(ibin,  egamma_qcdfakelep_transferfactor->GetBinError(ibin));   
		h_elefakepho_norm->SetBinError(ibin+NBIN,          mg_elefakepho_norm->GetBinError(ibin));
		h_jetfakepho_norm->SetBinError(ibin+NBIN,          mg_jetfakepho_norm->GetBinError(ibin));            
		h_qcdfakelep_norm->SetBinError(ibin+NBIN,          mg_qcdfakelep_norm->GetBinError(ibin));            
		h_VGamma_norm->SetBinError(ibin+NBIN,              mg_VGamma_norm->GetBinError(ibin));            
		h_rare_norm->SetBinError(ibin+NBIN,                mg_rare_norm->GetBinError(ibin));            
		h_elefakepho_transferfactor->SetBinError(ibin+NBIN,  mg_elefakepho_transferfactor->GetBinError(ibin));   
		h_jetfakepho_transferfactor->SetBinError(ibin+NBIN,  mg_jetfakepho_transferfactor->GetBinError(ibin));   
		h_qcdfakelep_transferfactor->SetBinError(ibin+NBIN,  mg_qcdfakelep_transferfactor->GetBinError(ibin));   
	}

	double totalerror[36]={15.745,50.100,19.984, 2.664,11.223, 8.607, 0.269, 0.466, 0.759, 1.936, 6.924, 3.647, 1.421, 3.230, 1.795, 0.194, 0.213, 0.502,
39.987,93.704,26.113, 3.992,14.471,13.323, 0.202, 1.356, 2.037, 2.362, 7.202, 4.811, 1.805, 3.183, 1.957, 0.400, 0.234, 0.209};
	for(int ibin(1); ibin <= 2*NBIN; ibin++){
		bkgContent=0; 
		bkgError = 0;
		bkgContent = h_elefakepho_norm->GetBinContent(ibin) + h_jetfakepho_norm->GetBinContent(ibin) + h_qcdfakelep_norm->GetBinContent(ibin) +
								 h_VGamma_norm->GetBinContent(ibin)     + h_rare_norm->GetBinContent(ibin);
		double bkgqcdfakelep = h_qcdfakelep_norm->GetBinContent(ibin);
		double bkgelefakepho = h_qcdfakelep_norm->GetBinContent(ibin) + h_elefakepho_norm->GetBinContent(ibin);
		double bkgjetfakepho = h_qcdfakelep_norm->GetBinContent(ibin) + h_elefakepho_norm->GetBinContent(ibin) + h_jetfakepho_norm->GetBinContent(ibin);
		double bkgVGamma = h_elefakepho_norm->GetBinContent(ibin) + h_jetfakepho_norm->GetBinContent(ibin) + h_qcdfakelep_norm->GetBinContent(ibin) + h_VGamma_norm->GetBinContent(ibin);

		if(h_elefakepho_syserr_eleshape->GetBinContent(ibin)>0)bkgError += pow(h_elefakepho_syserr_eleshape->GetBinContent(ibin), 2);
		if(h_jetfakepho_syserr_jetshape->GetBinContent(ibin)>0)bkgError += pow(h_jetfakepho_syserr_jetshape->GetBinContent(ibin), 2);
		if(h_qcdfakelep_syserr_scale->GetBinContent(ibin)   >0)bkgError += pow(h_qcdfakelep_syserr_scale->GetBinContent(ibin), 2);
		if(h_qcdfakelep_syserr_qcdshape->GetBinContent(ibin)>0)bkgError += pow(h_qcdfakelep_syserr_qcdshape->GetBinContent(ibin), 2);
								                                          															
		if(h_VGamma_syserr_jes->GetBinContent(ibin) >0)bkgError += pow(h_VGamma_syserr_jes->GetBinContent(ibin), 2);
		if(h_VGamma_syserr_jer->GetBinContent(ibin) >0)bkgError += pow(h_VGamma_syserr_jer->GetBinContent(ibin), 2);
		if(h_VGamma_syserr_esf->GetBinContent(ibin)  >0)bkgError += pow(h_VGamma_syserr_esf->GetBinContent(ibin), 2);
		if(h_VGamma_syserr_scale->GetBinContent(ibin)>0)bkgError += pow(h_VGamma_syserr_scale->GetBinContent(ibin), 2);
		if(h_VGamma_syserr_isr->GetBinContent(ibin)  >0)bkgError += pow(h_VGamma_syserr_isr->GetBinContent(ibin), 2);
		if(h_rare_syserr_jes->GetBinContent(ibin)    >0)bkgError += pow(h_rare_syserr_jes->GetBinContent(ibin), 2);
		if(h_rare_syserr_jer->GetBinContent(ibin)    >0)bkgError += pow(h_rare_syserr_jer->GetBinContent(ibin), 2);
		if(h_rare_syserr_esf->GetBinContent(ibin)    >0)bkgError += pow(h_rare_syserr_esf->GetBinContent(ibin), 2);
		if(h_rare_syserr_xs->GetBinContent(ibin)     >0)bkgError += pow(h_rare_syserr_xs->GetBinContent(ibin), 2);
		if(h_rare_syserr_lumi->GetBinContent(ibin)   >0)bkgError += pow(h_rare_syserr_lumi->GetBinContent(ibin), 2);

		bkgError = sqrt(bkgError);
		h_bkg->SetBinContent(ibin, bkgContent);
		h_bkg_qcdfakepho->SetBinContent(ibin, bkgqcdfakelep);
		h_bkg_elefakepho->SetBinContent(ibin, bkgelefakepho);
		h_bkg_jetfakepho->SetBinContent(ibin, bkgjetfakepho);
		h_bkg_VGamma->SetBinContent(ibin, bkgVGamma);
		h_bkg->SetBinError( ibin, totalerror[ibin-1]);
		error_bkg->SetPoint(ibin-1, ibin -1 + 0.5, bkgContent);
		error_bkg->SetPointError(ibin-1, 0.5, totalerror[ibin-1]);
		error_ratio->SetPoint(ibin-1, ibin -1 + 0.5, 1);
		error_ratio->SetPointError(ibin-1, 0.5, totalerror[ibin-1]/bkgContent);
	}

	TFile *outputfile = TFile::Open("SignalSystematic.root","RECREATE");
	outputfile->cd();
	h_elefakepho_norm->Write();       
	h_elefakepho_controlsample->Write();       
	h_elefakepho_transferfactor->Write();       
	h_elefakepho_syserr_jes->Write();       
	h_elefakepho_syserr_jer->Write();       
	h_elefakepho_syserr_esf->Write();       
	h_elefakepho_syserr_scale->Write();     
	h_elefakepho_syserr_eleshape->Write();  
	h_elefakepho_syserr_jetshape->Write();  
	h_elefakepho_syserr_qcdshape->Write();  
	h_elefakepho_syserr_xs->Write();        
	h_elefakepho_syserr_lumi->Write();      
	h_elefakepho_syserr_isr->Write();      
	h_jetfakepho_norm->Write();             
	h_jetfakepho_controlsample->Write();       
	h_jetfakepho_transferfactor->Write();       
	h_jetfakepho_syserr_jes->Write();       
	h_jetfakepho_syserr_jer->Write();       
	h_jetfakepho_syserr_esf->Write();       
	h_jetfakepho_syserr_scale->Write();     
	h_jetfakepho_syserr_eleshape->Write();  
	h_jetfakepho_syserr_jetshape->Write();  
	h_jetfakepho_syserr_qcdshape->Write();  
	h_jetfakepho_syserr_xs->Write();        
	h_jetfakepho_syserr_lumi->Write();      
	h_jetfakepho_syserr_isr->Write();      
	h_qcdfakelep_norm->Write();             
	h_qcdfakelep_controlsample->Write();       
	h_qcdfakelep_transferfactor->Write();       
	h_qcdfakelep_syserr_jes->Write();       
	h_qcdfakelep_syserr_jer->Write();       
	h_qcdfakelep_syserr_esf->Write();       
	h_qcdfakelep_syserr_scale->Write();     
	h_qcdfakelep_syserr_eleshape->Write();  
	h_qcdfakelep_syserr_jetshape->Write();  
	h_qcdfakelep_syserr_qcdshape->Write();  
	h_qcdfakelep_syserr_xs->Write();        
	h_qcdfakelep_syserr_lumi->Write();      
	h_qcdfakelep_syserr_isr->Write();      
	h_VGamma_norm->Write();             
	h_VGamma_syserr_jes->Write();       
	h_VGamma_syserr_jer->Write();       
	h_VGamma_syserr_esf->Write();       
	h_VGamma_syserr_scale->Write();     
	h_VGamma_syserr_eleshape->Write();  
	h_VGamma_syserr_jetshape->Write();  
	h_VGamma_syserr_qcdshape->Write();  
	h_VGamma_syserr_xs->Write();        
	h_VGamma_syserr_lumi->Write();      
	h_VGamma_syserr_isr->Write();      
	h_rare_norm->Write();             
	h_rare_syserr_jes->Write();       
	h_rare_syserr_jer->Write();       
	h_rare_syserr_esf->Write();       
	h_rare_syserr_scale->Write();     
	h_rare_syserr_eleshape->Write();  
	h_rare_syserr_jetshape->Write();  
	h_rare_syserr_qcdshape->Write();  
	h_rare_syserr_xs->Write();        
	h_rare_syserr_lumi->Write();      
	h_rare_syserr_isr->Write();      
	h_bkg->Write();
	error_bkg->SetName("bkgerror");
	error_bkg->Write();

	//TCanvas *can=new TCanvas("can","",1200,800);
	TCanvas *can=new TCanvas("can","",1200,800);
  gStyle->SetPadLeftMargin(0.10);
  can->SetLeftMargin(0.10);
  can->SetRightMargin(0.04);
	can->SetTopMargin(0.08);
	can->SetBottomMargin(0.1);
	can->cd();
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
	setTopPad(pad1); 
	pad1->Draw();  
	pad1->cd(); 
	gPad->SetLogy();
	h_sig->SetMarkerStyle(20); 
	h_sig->GetYaxis()->SetRangeUser(0.05,1000000);
	h_sig->GetXaxis()->SetRangeUser(0,2*NBIN);
	h_sig->GetYaxis()->SetTitle("Events / bin");
	h_sig->GetYaxis()->SetTitleOffset(1.0);
	h_sig->Draw("E0P");
  h_bkg->SetFillColor(kAzure-9);
	h_bkg->SetLineColor(kAzure-9);
	h_bkg_elefakepho->SetFillColor(kMagenta-7);
	h_bkg_jetfakepho->SetFillColor(kOrange-9);
	h_bkg_qcdfakepho->SetFillColor(kYellow);
	h_bkg_VGamma->SetFillColor(kGreen-9);
	h_bkg_elefakepho->SetLineColor(kMagenta-7);
	h_bkg_jetfakepho->SetLineColor(kOrange-9);
	h_bkg_qcdfakepho->SetLineColor(kYellow);
	h_bkg_VGamma->SetLineColor(kGreen-9);
	h_bkg->Draw("hist same");
	h_bkg_VGamma->Draw("hist same");
	h_bkg_jetfakepho->Draw("hist same");
	h_bkg_elefakepho->Draw("hist same");
	h_bkg_qcdfakepho->Draw("hist same");
	h_sig->SetLineWidth(1);
	h_sig->Draw("E0P same");
	TLegend *leg =  new TLegend(0.13,0.82,0.92,0.9);
	leg->SetNColumns(7);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	h_bkg->SetMarkerSize(0);
	h_bkg_elefakepho->SetMarkerSize(0);
	h_bkg_jetfakepho->SetMarkerSize(0);
	h_bkg_qcdfakepho->SetMarkerSize(0);
	h_bkg_VGamma->SetMarkerSize(0);
	error_bkg->SetMarkerSize(0);
	error_bkg->SetLineWidth(0);
	leg->AddEntry(h_sig,"Data","ep");
	leg->AddEntry(h_bkg,"t#bar{t}#gamma / WW#gamma / WZ#gamma");
	leg->AddEntry(h_bkg_VGamma, "W#gamma / Z#gamma");
	leg->AddEntry(h_bkg_jetfakepho,"j #rightarrow #gamma misid.");
	leg->AddEntry(h_bkg_elefakepho,"e #rightarrow #gamma misid.");
	leg->AddEntry(h_bkg_qcdfakepho,"Misid. leptons");
	leg->AddEntry(error_bkg, "Unc.");
	leg->Draw("same");
  error_bkg->SetFillColor(12);
  error_bkg->SetFillStyle(3345);
	error_bkg->Draw("E2 same");

	TLine *line_egamma_pt1 = new TLine(9,0,9,150000);
	TLine *line_egamma_pt2 = new TLine(18,0,18,150000);
	TLine *line_mg_pt1 = new TLine(27,0,27,150000);
	TLine *line_mg_pt2 = new TLine(36,0,36,150000);
	TLine *line_egamma_met1= new TLine(3,0,3,10000);
	TLine *line_egamma_met2= new TLine(6,0,6,10000);
	TLine *line_egamma_met3= new TLine(12,0,12,10000);
	TLine *line_egamma_met4= new TLine(15,0,15,10000);
	TLine *line_mg_met1= new TLine(21,0,21,10000);
	TLine *line_mg_met2= new TLine(24,0,24,10000);
	TLine *line_mg_met3= new TLine(30,0,30,10000);
	TLine *line_mg_met4= new TLine(33,0,33,10000);
	line_egamma_pt1->SetLineStyle(2);
	line_egamma_pt2->SetLineStyle(2);
	line_egamma_pt2->SetLineWidth(4);
	line_mg_pt1->SetLineStyle(2);
	line_mg_pt2->SetLineStyle(2);
	line_egamma_met1->SetLineStyle(3);
	line_egamma_met2->SetLineStyle(3);
	line_egamma_met3->SetLineStyle(3);
	line_egamma_met4->SetLineStyle(3);
	line_mg_met1->SetLineStyle(3);
	line_mg_met2->SetLineStyle(3);
	line_mg_met3->SetLineStyle(3);
	line_mg_met4->SetLineStyle(3);
	line_egamma_pt1->Draw("same");	
	line_egamma_pt2->Draw("same");	
	line_mg_pt1->Draw("same");	
	line_mg_pt2->Draw("same");	
	line_egamma_met1->Draw("same");	
	line_egamma_met2->Draw("same");	
	line_egamma_met3->Draw("same");	
	line_egamma_met4->Draw("same");	
	line_mg_met1->Draw("same");	
	line_mg_met2->Draw("same");	
	line_mg_met3->Draw("same");	
	line_mg_met4->Draw("same");	
	TLatex* latex = new TLatex();
	latex->SetTextSize(0.04);
	latex->DrawLatex(1, 50000,"e#gamma, p_{T}^{#gamma} < 200 GeV");
	latex->DrawLatex(9.5,50000,"e#gamma, p_{T}^{#gamma} > 200 GeV");
	latex->DrawLatex(18.5,50000,"#mu#gamma, p_{T}^{#gamma} < 200 GeV");
	latex->DrawLatex(27.5,50000,"#mu#gamma, p_{T}^{#gamma} > 200 GeV");
	latex->SetTextSize(0.04);
	latex->DrawLatex(2.8, 15000,"p^{miss}_{T} (GeV)");
	latex->DrawLatex(11.8, 15000,"p^{miss}_{T} (GeV)");
	latex->DrawLatex(20.8, 15000,"p^{miss}_{T} (GeV)");
	latex->DrawLatex(29.8, 15000,"p^{miss}_{T} (GeV)");
	latex->SetTextSize(0.04);
	latex->DrawLatex(1, 5000, "< 200");
	latex->DrawLatex(3.2, 5000, "200--400");
	latex->DrawLatex(6.5, 5000, "> 400");
	latex->DrawLatex(9.5, 5000, "< 200");
	latex->DrawLatex(12.2, 5000, "200--400");
	latex->DrawLatex(15.5, 5000, "> 400");
	latex->DrawLatex(18.5, 5000, "< 200");
	latex->DrawLatex(21.2, 5000, "200--400");
	latex->DrawLatex(24.5, 5000, "> 400");
	latex->DrawLatex(27.5, 5000, "< 200");
	latex->DrawLatex(30.2, 5000, "200--400");
	latex->DrawLatex(33.5, 5000, "> 400");
  CMS_lumi( pad1, 0 );
 	gPad->RedrawAxis();

	can->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	TH1F *ratioframe = new TH1F("frame",";Search region;#frac{Data}{Bkg.}",2*NBIN,0,2*NBIN);
	ratioframe->GetXaxis()->SetRangeUser(0,2*NBIN);
  ratioframe->GetYaxis()->SetTitleOffset(1.0);
	ratioframe->GetYaxis()->SetRangeUser(0,4);
	std::ostringstream ratiotitle;
	for(int i(1); i <=2*NBIN; i++){
		ratiotitle.str("");
		ratiotitle << i;
		ratioframe->GetXaxis()->SetBinLabel(i,ratiotitle.str().c_str());
	}
	ratioframe->Draw();
  TLine *flatratio = new TLine(0,1,2*NBIN,1);
	double alpha = 1-0.6827;
	for(int ibin(1); ibin <= 2*NBIN; ibin++){
		int N = h_sig->GetBinContent(ibin);
		ratio->SetPoint(ibin-1, ibin-0.5, N/h_bkg->GetBinContent(ibin));
		double L = (N==0)? 0: (ROOT::Math::gamma_quantile(alpha/2, N, 1.));
		double U = (ROOT::Math::gamma_quantile(1-alpha/2, N+1, 1.));
		std::cout << ibin << " N " << N << " L " << L << " U " << U << std::endl;
    ratio->SetPointEYlow( ibin-1, (N-L)/h_bkg->GetBinContent(ibin));
		ratio->SetPointEYhigh(ibin-1, (U-N)/h_bkg->GetBinContent(ibin));
  }
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->GetYaxis()->SetRangeUser(0,4);
	ratio->Draw("P same");
  error_ratio->SetFillColor(12);
  error_ratio->SetFillStyle(3345);
	error_ratio->Draw("E2 same");
	flatratio->Draw("same");
	can->Update();
	can->SaveAs("signalCount.pdf");

	for(int ibin(1); ibin <= 2*NBIN; ibin++){  }
	h_sig->Write();
	outputfile->Close();


	TH1D *test_h = new TH1D("test_h","",NBIN,0,NBIN);
	test_h->SetBinErrorOption(TH1::kPoisson);
	for(int ibin(1); ibin <= NBIN; ibin++){test_h->SetBinContent(ibin, ibin-1);}
	for(int ibin(2); ibin <= NBIN; ibin++){
		int N = test_h->GetBinContent(ibin);
		double L = (N==0)? 0: (ROOT::Math::gamma_quantile(alpha/2, N, 1.));
		double U = (ROOT::Math::gamma_quantile(1-alpha/2, N+1, 1.));
		std::cout << ibin << " N " << N << " L " << L << " U " << U << " GetBinErrorLow " << test_h->GetBinErrorLow(ibin) << " up " << test_h->GetBinErrorUp(ibin) <<  std::endl;
  }

}
