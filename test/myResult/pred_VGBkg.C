#include "../../include/analysis_commoncode.h"
#include "invmass_binning.h"

#include "TProfile2D.h"
#include <set>
bool useVarSF = true;
// Per-mass-bin version of the "added by hand" JES/JER fix that the signal-bin loop
// applies: if a JES/JER variation exceeds 100% of the nominal yield of a bin, that bin's
// nominal yield is replaced by max(jesUp,jesDown,jerUp,jerDown)/2. Kept ON so the
// mass-bin prediction is treated like the signal-bin one; set to false to see raw yields.
bool applyLowStatJESJERFix = true;
// HT-binned V+gamma factor is replaced by the fixed (!useVarSF) factors if it is negative,
// or if its relative error (binError/binContent) exceeds this value
const double maxRelErrVarSF = 0.5;
bool applyISRwt = true;
// LO -> NLO W+gamma ISR pT shape factors (NLO/LO in ISRJetPtReco, unit area, Run 2 preselection; AN App. H.2,
// ISRReweight/compare_WG_LOvsNLO.C, /eos/.../ISRweighting/WGLOvsNLO_factors_<ch>.txt).
// Applied to W+gamma only (mcType 1-3), on top of the data/NLO-Z+gamma ISR weights (all V+gamma); they keep the
// LO preselection yield unchanged, so no normalization is needed.
// (set NO_LOTONLO=1 in the environment to switch it off, e.g. for the validation without it)
bool applyLOtoNLO = (getenv("NO_LOTONLO") == nullptr);
double LOtoNLOFactor(int channelType, float isrPt){
	static const double edges[7] = {0, 50, 100, 150, 200, 250, 300};
	static const double eg[7] = {1.06025, 0.94482, 0.89026, 0.89977, 0.89691, 0.89908, 0.87877};
	static const double mg[7] = {1.06165, 0.93595, 0.89678, 0.87963, 0.85610, 0.86764, 0.85067};
	int i = 6;
	while(i > 0 && isrPt < edges[i]) i--;
	return channelType == 1 ? eg[i] : mg[i];
}
void pred_VGBkg(){
	bool toDeriveScale(false);
	SetSignalConfig();
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
	
	esfScaleFactor  objectESF;
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);

	if(anatype == 0)toDeriveScale = true;

	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
	else whichVFP = "";

	int channelType = ichannel; // eg = 1; mg =2;
        double factorMC(1);
        double factorMCUP(1);
        double factorMC_fixed(1);   // fixed factors (fallback when HT-binned factor is negative)
        double factorMCUP_fixed(1);

        TFile *NormFile= nullptr;
        TH1D *hVGamma= nullptr;

         // make this ON when you have to derive the QCD and VGamma scales
        if(!toDeriveScale){
                if(useVarSF){
                        if(channelType == 1)NormFile = TFile::Open(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/VGammaScaleHist_eg_%d%s.root",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str() ));
                        else if(channelType == 2)NormFile = TFile::Open(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/VGammaScaleHist_mg_%d%s.root",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
                        hVGamma = (TH1D*)NormFile->Get("VGamma_ht");}
                // fixed (non-HT-binned) factors: always filled; used directly when !useVarSF,
                // and as the fallback for HT bins with negative factorMC when useVarSF
                        if(channelType == 1){
                                if(RunYear==2016 and preVFP == 1){
                                        factorMC = factor_egVGamma_2016preVFP;
                                        factorMCUP = factor_egVGamma_2016preVFP+factorerror_egVGamma_2016preVFP;}
                                if(RunYear==2016 and preVFP == 0){
                                        factorMC = factor_egVGamma_2016postVFP;
                                        factorMCUP = factor_egVGamma_2016postVFP+factorerror_egVGamma_2016postVFP;}
                                if(RunYear==2017){
                                        factorMC = factor_egVGamma_2017;
                                        factorMCUP = factor_egVGamma_2017+factorerror_egVGamma_2017;}
                                if(RunYear==2018){
                                        factorMC = factor_egVGamma_2018;
                                        factorMCUP = factor_egVGamma_2018+factorerror_egVGamma_2018;}
                        }
                        else if(channelType == 2){
                                if(RunYear==2016 and preVFP == 1){
                                        factorMC = factor_mgVGamma_2016preVFP;
                                        factorMCUP = factor_mgVGamma_2016preVFP+factorerror_mgVGamma_2016preVFP;}
                                if(RunYear==2016 and preVFP == 0){
                                        factorMC = factor_mgVGamma_2016postVFP;
                                        factorMCUP = factor_mgVGamma_2016postVFP+factorerror_mgVGamma_2016postVFP;}
                                if(RunYear==2017){
                                        factorMC = factor_mgVGamma_2017;
                                        factorMCUP = factor_mgVGamma_2017+factorerror_mgVGamma_2017;}
                                if(RunYear==2018){
                                        factorMC = factor_mgVGamma_2018;
                                        factorMCUP = factor_mgVGamma_2018+factorerror_mgVGamma_2018;}
                        }
                factorMC_fixed   = factorMC;
                factorMCUP_fixed = factorMCUP;
        }

	//*********** histo list **********************
	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_VGBkg";
	else if(channelType==2)outputname << "mg_VGBkg";
	if(anatype ==0)outputname << "met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << "_" <<RunYear<<whichVFP<<".root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	std::ostringstream histname;

	TH1D *p_PhoEt = new TH1D("p_PhoEt","; p_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nSigPtBins,sigPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	// yields per mass bin for the nominal prediction and for each variation that enters the
	// total uncertainty in plot_eventct.C (JES, JER, ESF, norm scale, QCD muR/muF, PDF),
	// filled for the same events as the signal-bin histograms below.
	const std::string chPrefix = (channelType==1) ? "eg" : "mg";
	auto newMassHist = [&](const std::string &name, const char *title = ""){
		return new TH1D((chPrefix+"_VGamma_invmass_"+name).c_str(), title, nInvMassBins, invMassMin, invMassMax);
	};
	TH1D *h_VGamma_invmass_norm     = newMassHist("norm", "m_{l#gamma};m_{l#gamma} (GeV);Events / bin");
	TH1D *h_VGamma_invmass_normUp   = newMassHist("normUp");
	TH1D *h_VGamma_invmass_esfUp    = newMassHist("esfUp");
	TH1D *h_VGamma_invmass_esfphoUp    = newMassHist("esfphoUp");
	TH1D *h_VGamma_invmass_jesUp    = newMassHist("jesUp");
	TH1D *h_VGamma_invmass_jesDown  = newMassHist("jesDown");
	TH1D *h_VGamma_invmass_jerUp    = newMassHist("jerUp");
	TH1D *h_VGamma_invmass_jerDown  = newMassHist("jerDown");
	TH1D *h_VGamma_invmass_muRup    = newMassHist("muRup");
	TH1D *h_VGamma_invmass_muRdown  = newMassHist("muRdown");
	TH1D *h_VGamma_invmass_muFup    = newMassHist("muFup");
	TH1D *h_VGamma_invmass_muFdown  = newMassHist("muFdown");
	TH1D *h_VGamma_invmass_pdfUp    = newMassHist("pdfUp");
	TH1D *h_VGamma_invmass_pdfDown  = newMassHist("pdfDown");
	// per-bin uncertainties (same conventions as h_VGamma_syserr_*: jes/jer/esf/scale are absolute
	// event yields; qcdmuR/qcdmuF/pdf are FRACTIONAL, i.e. |kappa-1|)
	TH1D *h_VGamma_invmass_syserr_jes    = newMassHist("syserr_jes");
	TH1D *h_VGamma_invmass_syserr_jer    = newMassHist("syserr_jer");
	TH1D *h_VGamma_invmass_syserr_esf    = newMassHist("syserr_esf");
	TH1D *h_VGamma_invmass_syserr_esfpho    = newMassHist("syserr_esfpho");
	TH1D *h_VGamma_invmass_syserr_scale  = newMassHist("syserr_scale");
	TH1D *h_VGamma_invmass_syserr_qcdmuR = newMassHist("syserr_qcdmuR");
	TH1D *h_VGamma_invmass_syserr_qcdmuF = newMassHist("syserr_qcdmuF");
	TH1D *h_VGamma_invmass_syserr_pdf    = newMassHist("syserr_pdf");
	h_VGamma_invmass_norm->Sumw2();

	TH1D *h_VGamma_norm = 0;
	TH1D *h_VGamma_jesUp = 0;
	TH1D *h_VGamma_jesDown = 0;
	TH1D *h_VGamma_jerUp = 0;
	TH1D *h_VGamma_jerDown = 0;
	TH1D *h_VGamma_esfUp = 0;
	TH1D *h_VGamma_esfphoUp = 0;
	TH1D *h_VGamma_normUp = 0;
	TH1D *h_VGamma_isrAlter = 0;
	TH1D *h_VGamma_nloAlter = 0;
	TH1D *h_VGamma_syserr_jes = 0;
	TH1D *h_VGamma_syserr_jer = 0;
	TH1D *h_VGamma_syserr_esf = 0;
	TH1D *h_VGamma_syserr_esfpho = 0;
	TH1D *h_VGamma_syserr_scale = 0;
	TH1D *h_VGamma_syserr_eleshape = 0;
	TH1D *h_VGamma_syserr_jetshape = 0;
	TH1D *h_VGamma_syserr_qcdshape = 0;
	TH1D *h_VGamma_syserr_xs = 0;
	TH1D *h_VGamma_syserr_lumi = 0;
	TH1D *h_VGamma_syserr_isr = 0;
	TH1D *h_VGamma_syserr_wgnlo = 0;
	TH1D *h_VGamma_puUp = 0;
	TH1D *h_VGamma_puDown = 0;
	TH1D *h_VGamma_syserr_pu = 0;
	TH1D *h_VGamma_muRup         = 0;
	TH1D *h_VGamma_muRdown       = 0;
	TH1D *h_VGamma_muFup         = 0;
	TH1D *h_VGamma_muFdown       = 0;
	TH1D *h_VGamma_pdfUp         = 0;
	TH1D *h_VGamma_pdfDown       = 0;
	TH1D *h_VGamma_syserr_qcdmuR = 0;
	TH1D *h_VGamma_syserr_qcdmuF = 0;
	TH1D *h_VGamma_syserr_pdf    = 0;

	if(channelType==1){
		h_VGamma_norm            = new TH1D("eg_VGamma_norm","eventcount",NBIN,0,NBIN);
		h_VGamma_jesUp           = new TH1D("eg_VGamma_jesUp","eventcount",NBIN,0,NBIN);
		h_VGamma_jesDown         = new TH1D("eg_VGamma_jesDown","eventcount",NBIN,0,NBIN);
		h_VGamma_jerUp           = new TH1D("eg_VGamma_jerUp","eventcount",NBIN,0,NBIN);
		h_VGamma_jerDown         = new TH1D("eg_VGamma_jerDown","eventcount",NBIN,0,NBIN);
		h_VGamma_esfUp           = new TH1D("eg_VGamma_esfUp","eventcount",NBIN,0,NBIN);
		h_VGamma_esfphoUp           = new TH1D("eg_VGamma_esfphoUp","eventcount",NBIN,0,NBIN);
		h_VGamma_normUp          = new TH1D("eg_VGamma_normUp","eventcount",NBIN,0,NBIN);
		h_VGamma_isrAlter        = new TH1D("eg_VGamma_isrAlter","eventcount",NBIN,0,NBIN);
		h_VGamma_nloAlter        = new TH1D("eg_VGamma_nloAlter","eventcount",NBIN,0,NBIN);
		h_VGamma_puUp           = new TH1D("eg_VGamma_puUp","eventcount",NBIN,0,NBIN);
		h_VGamma_puDown         = new TH1D("eg_VGamma_puDown","eventcount",NBIN,0,NBIN);
		h_VGamma_syserr_jes      = new TH1D("eg_VGamma_syserr_jes","",NBIN,0,NBIN);	
		h_VGamma_syserr_jer      = new TH1D("eg_VGamma_syserr_jer","",NBIN,0,NBIN);	
		h_VGamma_syserr_esf      = new TH1D("eg_VGamma_syserr_esf","",NBIN,0,NBIN);	
		h_VGamma_syserr_esfpho      = new TH1D("eg_VGamma_syserr_esfpho","",NBIN,0,NBIN);	
		h_VGamma_syserr_scale    = new TH1D("eg_VGamma_syserr_scale","",NBIN,0,NBIN);	
		h_VGamma_syserr_eleshape = new TH1D("eg_VGamma_syserr_eleshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_jetshape = new TH1D("eg_VGamma_syserr_jetshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_qcdshape = new TH1D("eg_VGamma_syserr_qcdshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_xs       = new TH1D("eg_VGamma_syserr_xs","",NBIN,0,NBIN);	
		h_VGamma_syserr_lumi     = new TH1D("eg_VGamma_syserr_lumi","",NBIN,0,NBIN);
		h_VGamma_syserr_isr      = new TH1D("eg_VGamma_syserr_isr","",NBIN,0,NBIN);
		h_VGamma_syserr_wgnlo    = new TH1D("eg_VGamma_syserr_wgnlo","",NBIN,0,NBIN);
		h_VGamma_syserr_pu       = new TH1D("eg_VGamma_syserr_pu","",NBIN,0,NBIN);
		h_VGamma_muRup         = new TH1D("eg_VGamma_muRup","muR up yield",NBIN,0,NBIN);
		h_VGamma_muRdown       = new TH1D("eg_VGamma_muRdown","muR down yield",NBIN,0,NBIN);
		h_VGamma_muFup         = new TH1D("eg_VGamma_muFup","muF up yield",NBIN,0,NBIN);
		h_VGamma_muFdown       = new TH1D("eg_VGamma_muFdown","muF down yield",NBIN,0,NBIN);
		h_VGamma_pdfUp         = new TH1D("eg_VGamma_pdfUp","PDF up yield",NBIN,0,NBIN);
		h_VGamma_pdfDown       = new TH1D("eg_VGamma_pdfDown","PDF down yield",NBIN,0,NBIN);
		h_VGamma_syserr_qcdmuR = new TH1D("eg_VGamma_syserr_qcdmuR","",NBIN,0,NBIN);
		h_VGamma_syserr_qcdmuF = new TH1D("eg_VGamma_syserr_qcdmuF","",NBIN,0,NBIN);
		h_VGamma_syserr_pdf    = new TH1D("eg_VGamma_syserr_pdf","",NBIN,0,NBIN);
	} 
	else if(channelType==2){
		h_VGamma_norm            = new TH1D("mg_VGamma_norm","eventcount",NBIN,0,NBIN);
		h_VGamma_jesUp           = new TH1D("mg_VGamma_jesUp","eventcount",NBIN,0,NBIN);
		h_VGamma_jesDown         = new TH1D("mg_VGamma_jesDown","eventcount",NBIN,0,NBIN);
		h_VGamma_jerUp           = new TH1D("mg_VGamma_jerUp","eventcount",NBIN,0,NBIN);
		h_VGamma_jerDown         = new TH1D("mg_VGamma_jerDown","eventcount",NBIN,0,NBIN);
		h_VGamma_esfUp           = new TH1D("mg_VGamma_esfUp","eventcount",NBIN,0,NBIN);
		h_VGamma_esfphoUp           = new TH1D("mg_VGamma_esfphoUp","eventcount",NBIN,0,NBIN);
		h_VGamma_normUp          = new TH1D("mg_VGamma_normUp","eventcount",NBIN,0,NBIN);
		h_VGamma_isrAlter        = new TH1D("mg_VGamma_isrAlter","eventcount",NBIN,0,NBIN);
		h_VGamma_nloAlter        = new TH1D("mg_VGamma_nloAlter","eventcount",NBIN,0,NBIN);
		h_VGamma_puUp            = new TH1D("mg_VGamma_puUp","eventcount",NBIN,0,NBIN);
		h_VGamma_puDown          = new TH1D("mg_VGamma_puDown","eventcount",NBIN,0,NBIN);
		h_VGamma_syserr_jes      = new TH1D("mg_VGamma_syserr_jes","",NBIN,0,NBIN);	
		h_VGamma_syserr_jer      = new TH1D("mg_VGamma_syserr_jer","",NBIN,0,NBIN);	
		h_VGamma_syserr_esf      = new TH1D("mg_VGamma_syserr_esf","",NBIN,0,NBIN);	
		h_VGamma_syserr_esfpho      = new TH1D("mg_VGamma_syserr_esfpho","",NBIN,0,NBIN);	
		h_VGamma_syserr_scale    = new TH1D("mg_VGamma_syserr_scale","",NBIN,0,NBIN);	
		h_VGamma_syserr_eleshape = new TH1D("mg_VGamma_syserr_eleshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_jetshape = new TH1D("mg_VGamma_syserr_jetshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_qcdshape = new TH1D("mg_VGamma_syserr_qcdshape","",NBIN,0,NBIN);	
		h_VGamma_syserr_xs       = new TH1D("mg_VGamma_syserr_xs","",NBIN,0,NBIN);	
		h_VGamma_syserr_lumi     = new TH1D("mg_VGamma_syserr_lumi","",NBIN,0,NBIN);
		h_VGamma_syserr_isr      = new TH1D("mg_VGamma_syserr_isr","",NBIN,0,NBIN);
		h_VGamma_syserr_wgnlo    = new TH1D("mg_VGamma_syserr_wgnlo","",NBIN,0,NBIN);
		h_VGamma_syserr_pu       = new TH1D("mg_VGamma_syserr_pu","",NBIN,0,NBIN);
		h_VGamma_muRup         = new TH1D("mg_VGamma_muRup","muR up yield",NBIN,0,NBIN);
		h_VGamma_muRdown       = new TH1D("mg_VGamma_muRdown","muR down yield",NBIN,0,NBIN);
		h_VGamma_muFup         = new TH1D("mg_VGamma_muFup","muF up yield",NBIN,0,NBIN);
		h_VGamma_muFdown       = new TH1D("mg_VGamma_muFdown","muF down yield",NBIN,0,NBIN);
		h_VGamma_pdfUp         = new TH1D("mg_VGamma_pdfUp","PDF up yield",NBIN,0,NBIN);
		h_VGamma_pdfDown       = new TH1D("mg_VGamma_pdfDown","PDF down yield",NBIN,0,NBIN);
		h_VGamma_syserr_qcdmuR = new TH1D("mg_VGamma_syserr_qcdmuR","",NBIN,0,NBIN);
		h_VGamma_syserr_qcdmuF = new TH1D("mg_VGamma_syserr_qcdmuF","",NBIN,0,NBIN);
		h_VGamma_syserr_pdf    = new TH1D("mg_VGamma_syserr_pdf","",NBIN,0,NBIN);
	} 

	TH1D *jesup_MET = new TH1D("jesup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jesup_Mt = new TH1D("jesup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *jesup_HT = new TH1D("jesup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *jesup_dPhiEleMET = new TH1D("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jesdo_MET = new TH1D("jesdo_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jesdo_Mt = new TH1D("jesdo_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jesdo_HT = new TH1D("jesdo_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *jesdo_dPhiEleMET = new TH1D("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerup_MET = new TH1D("jerup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jerup_Mt = new TH1D("jerup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jerup_dPhiEleMET = new TH1D("jerup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *jerdo_MET = new TH1D("jerdo_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *jerdo_Mt = new TH1D("jerdo_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *jerdo_dPhiEleMET = new TH1D("jerdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *scaleup_PhoEt = new TH1D("scaleup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *scaleup_PhoEta = new TH1D("scaleup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *scaleup_LepPt = new TH1D("scaleup_LepPt","scaleup_LepPt",nSigPtBins,sigPtBins);
	TH1D *scaleup_LepEta = new TH1D("scaleup_LepEta","scaleup_LepEta",60,-3,3);
	TH1D *scaleup_MET = new TH1D("scaleup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *scaleup_Mt = new TH1D("scaleup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *scaleup_HT = new TH1D("scaleup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *scaleup_dPhiEleMET = new TH1D("scaleup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nSigPtBins,sigPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *isrup_PhoEt = new TH1D("isrup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *isrup_PhoEta = new TH1D("isrup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *isrup_LepPt = new TH1D("isrup_LepPt","isrup_LepPt",nSigPtBins,sigPtBins);
	TH1D *isrup_LepEta = new TH1D("isrup_LepEta","isrup_LepEta",60,-3,3);
	TH1D *isrup_MET = new TH1D("isrup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *isrup_Mt = new TH1D("isrup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *isrup_HT = new TH1D("isrup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	TH1D *isrup_dPhiEleMET = new TH1D("isrup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *puup_PhoEt = new TH1D("puup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *puup_MET = new TH1D("puup_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *puup_Mt = new TH1D("puup_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *puup_HT = new TH1D("puup_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
    	
	TH1D *pudo_PhoEt = new TH1D("pudo_PhoEt","#gamma E_{T}; E_{T} (GeV)",nSigEtBins,sigEtBins);
	TH1D *pudo_MET = new TH1D("pudo_MET","MET; MET (GeV);",nSigMETBins, sigMETBins);
	TH1D *pudo_Mt = new TH1D("pudo_Mt","M_{T}; M_{T} (GeV);",nSigMtBins,sigMtBins); 
	TH1D *pudo_HT = new TH1D("pudo_HT","HT; HT (GeV);",nSigHTBins, sigHTBins); 
	
    	float psWeightIsrUp(1), psWeightIsrDown(1);
    	TH1D *pswIsrUp_MET   = new TH1D("pswIsrUp_MET",   "MET ISRup;   p_{T}^{miss}(GeV);", nSigMETBins, sigMETBins);
    	TH1D *pswIsrDo_MET   = new TH1D("pswIsrDo_MET",   "MET ISRdo;   p_{T}^{miss}(GeV);", nSigMETBins, sigMETBins);
	TH1D *pswIsrUp_Mt    = new TH1D("pswIsrUp_Mt",	  "M_{T}; 	M_{T} (GeV);",	     nSigMtBins,sigMtBins); 
	TH1D *pswIsrDo_Mt    = new TH1D("pswIsrDo_Mt",	  "M_{T}; 	M_{T} (GeV);",	     nSigMtBins,sigMtBins); 
    	TH1D *pswIsrUp_HT    = new TH1D("pswIsrUp_HT",    "HT ISRup;    H_{T}(GeV);",         nSigHTBins,  sigHTBins);
    	TH1D *pswIsrDo_HT    = new TH1D("pswIsrDo_HT",    "HT ISRdo;    H_{T}(GeV);",         nSigHTBins,  sigHTBins);
    	TH1D *pswIsrUp_PhoEt = new TH1D("pswIsrUp_PhoEt", "phoEt ISRup; p_{T}^{#gamma}(GeV);",nSigEtBins,  sigEtBins);
    	TH1D *pswIsrDo_PhoEt = new TH1D("pswIsrDo_PhoEt", "phoEt ISRdo; p_{T}^{#gamma}(GeV);",nSigEtBins,  sigEtBins);
    	TH1D *h_VGamma_pswIsrUp   = new TH1D((std::string(channelType==1?"eg":"mg")+"_VGamma_pswIsrUp").c_str(),   "eventcount",NBIN,0,NBIN);
    	TH1D *h_VGamma_pswIsrDown = new TH1D((std::string(channelType==1?"eg":"mg")+"_VGamma_pswIsrDown").c_str(), "eventcount",NBIN,0,NBIN);

    	float psWeightFsrUp(1), psWeightFsrDown(1);
    	TH1D *pswFsrUp_MET   = new TH1D("pswFsrUp_MET",   "MET ISRup;   p_{T}^{miss}(GeV);", nSigMETBins, sigMETBins);
    	TH1D *pswFsrDo_MET   = new TH1D("pswFsrDo_MET",   "MET ISRdo;   p_{T}^{miss}(GeV);", nSigMETBins, sigMETBins);
	TH1D *pswFsrUp_Mt    = new TH1D("pswFsrUp_Mt",	  "M_{T}; 	M_{T} (GeV);",	     nSigMtBins,sigMtBins); 
	TH1D *pswFsrDo_Mt    = new TH1D("pswFsrDo_Mt",	  "M_{T}; 	M_{T} (GeV);",	     nSigMtBins,sigMtBins); 
    	TH1D *pswFsrUp_HT    = new TH1D("pswFsrUp_HT",    "HT ISRup;    H_{T}(GeV);",         nSigHTBins,  sigHTBins);
    	TH1D *pswFsrDo_HT    = new TH1D("pswFsrDo_HT",    "HT ISRdo;    H_{T}(GeV);",         nSigHTBins,  sigHTBins);
    	TH1D *pswFsrUp_PhoEt = new TH1D("pswFsrUp_PhoEt", "phoEt ISRup; p_{T}^{#gamma}(GeV);",nSigEtBins,  sigEtBins);
    	TH1D *pswFsrDo_PhoEt = new TH1D("pswFsrDo_PhoEt", "phoEt ISRdo; p_{T}^{#gamma}(GeV);",nSigEtBins,  sigEtBins);
    	TH1D *h_VGamma_pswFsrUp   = new TH1D((std::string(channelType==1?"eg":"mg")+"_VGamma_pswFsrUp").c_str(),   "eventcount",NBIN,0,NBIN);
    	TH1D *h_VGamma_pswFsrDown = new TH1D((std::string(channelType==1?"eg":"mg")+"_VGamma_pswFsrDown").c_str(), "eventcount",NBIN,0,NBIN);
	
	// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  	TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
	
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_%d%s.root",RunYear,whichVFP.c_str()));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
  	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_%d%s.root",RunYear,whichVFP.c_str()));
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_%d%s.root",RunYear,whichVFP.c_str()));
  	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
	
	float crosssection(0);
	float ntotalevent(0);
	float ISRWeight(0);
	int   mcType(0);
	float PUweight(1);
	float PUweightUp(1);
	float PUweightDown(1);
        float phoEt(0);
        float phoEta(0);
        float phoPhi(0);
        float lepPt(0);
        float lepEta(0);
        float lepPhi(0);
        float sigMT(0);
        float sigMET(0);
        float sigMETPhi(0);
        float dPhiLepMET(0);
        int   nVertex(0);
        float dRPhoLep(0);
        float HT(0);
        float nJet(0);
	float llmass(0);
	float invmass(0); // m(lepton,photon) stored in the input tree
	float ISRJetPt(0); 
	float sigMETJESup(0);
	float sigMETJESdo(0);
	float sigMETJERup(0);
	float sigMETJERdo(0);
	float sigMTJESup(0);
	float sigMTJESdo(0);
	float sigMTJERup(0);
	float sigMTJERdo(0);
	float HTJESup(0);
	float HTJESdo(0);
	float dPhiLepMETJESup(0);
	float dPhiLepMETJESdo(0);
	float dPhiLepMETJERup(0);
	float dPhiLepMETJERdo(0);
  	std::vector<int> *mcPID=0;
  	std::vector<float> *mcEta=0;
  	std::vector<float> *mcPhi=0;
  	std::vector<float> *mcPt=0;
  	std::vector<int> *mcMomPID=0;
	//  std::vector<int> *mcGMomPID=0;

	
	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
//	mctree->SetBranchAddress("ISRWeight", &ISRWeight);
	mctree->SetBranchAddress("mcType",    &mcType);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	mctree->SetBranchAddress("PUweightUp",  &PUweightUp);
	mctree->SetBranchAddress("PUweightDown",  &PUweightDown);
	// QCD scale and PDF scalar weights (for bin-wise lnN, from analysis_VGamma_Run2.C)
	float muRweightUp(1.f), muRweightDown(1.f), muFweightUp(1.f), muFweightDown(1.f);
	float pdfWeightRMS(0.f);

	mctree->SetBranchAddress("muRweightUp",   &muRweightUp);
	mctree->SetBranchAddress("muRweightDown", &muRweightDown);
	mctree->SetBranchAddress("muFweightUp",   &muFweightUp);
	mctree->SetBranchAddress("muFweightDown", &muFweightDown);
	mctree->SetBranchAddress("pdfWeightRMS",  &pdfWeightRMS);
	mctree->SetBranchAddress("psWeightIsrUp",   &psWeightIsrUp);
    	mctree->SetBranchAddress("psWeightIsrDown", &psWeightIsrDown);
	mctree->SetBranchAddress("psWeightFsrUp",   &psWeightFsrUp);
    	mctree->SetBranchAddress("psWeightFsrDown", &psWeightFsrDown);

  	mctree->SetBranchAddress("phoEt",     &phoEt);
  	mctree->SetBranchAddress("phoEta",    &phoEta);
  	mctree->SetBranchAddress("phoPhi",    &phoPhi);
  	mctree->SetBranchAddress("lepPt",     &lepPt);
  	mctree->SetBranchAddress("lepEta",    &lepEta);
  	mctree->SetBranchAddress("lepPhi",    &lepPhi);
  	mctree->SetBranchAddress("sigMT",     &sigMT);
  	mctree->SetBranchAddress("sigMET",    &sigMET);
  	mctree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  	mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  	mctree->SetBranchAddress("nVertex",   &nVertex);
  	mctree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	mctree->SetBranchAddress("llmass",    &llmass);
	mctree->SetBranchAddress("invmass",   &invmass);
  	mctree->SetBranchAddress("HT",        &HT);
  	mctree->SetBranchAddress("nJet",      &nJet);
  	mctree->SetBranchAddress("ISRJetPt",     &ISRJetPt);
	mctree->SetBranchAddress("sigMETJESup",     &sigMETJESup);
	mctree->SetBranchAddress("sigMETJESdo",     &sigMETJESdo);
	mctree->SetBranchAddress("sigMETJERup",     &sigMETJERup);
	mctree->SetBranchAddress("sigMETJERdo",     &sigMETJERdo);
	mctree->SetBranchAddress("sigMTJESup",      &sigMTJESup);
	mctree->SetBranchAddress("sigMTJESdo",      &sigMTJESdo);
	mctree->SetBranchAddress("sigMTJERup",      &sigMTJERup);
	mctree->SetBranchAddress("sigMTJERdo",      &sigMTJERdo);
	mctree->SetBranchAddress("HTJESup",         &HTJESup);
	mctree->SetBranchAddress("HTJESdo",         &HTJESdo);
	mctree->SetBranchAddress("dPhiLepMETJESup", &dPhiLepMETJESup);
	mctree->SetBranchAddress("dPhiLepMETJESdo", &dPhiLepMETJESdo);
	mctree->SetBranchAddress("dPhiLepMETJERup", &dPhiLepMETJERup);
	mctree->SetBranchAddress("dPhiLepMETJERdo", &dPhiLepMETJERdo);
 	mctree->SetBranchAddress("mcPID",     &mcPID);
  	mctree->SetBranchAddress("mcEta",     &mcEta);
  	mctree->SetBranchAddress("mcPhi",     &mcPhi);
 	mctree->SetBranchAddress("mcPt",      &mcPt);
  	mctree->SetBranchAddress("mcMomPID",  &mcMomPID);

	// yields with the photon-pT binned V+gamma factor of bin k (0-100, 100-200, >200 GeV) at its -1 / +1 sigma and
	// the other bins nominal: one normalisation nuisance per year x channel x photon-pT bin (createDataCard.py)
	TH1D *h_VGamma_sfUp[3], *h_VGamma_sfDown[3];
	for(int k = 0; k < 3; k++){
		h_VGamma_sfUp[k]   = new TH1D(Form("%s_VGamma_sfUp_pho%d", chPrefix.c_str(), k), "", NBIN, 0, NBIN);
		h_VGamma_sfDown[k] = new TH1D(Form("%s_VGamma_sfDown_pho%d", chPrefix.c_str(), k), "", NBIN, 0, NBIN);
		h_VGamma_sfUp[k]->SetDirectory(outputfile); h_VGamma_sfDown[k]->SetDirectory(outputfile);}
	int nNoSFbin = 0;

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(0);
		double scalefactorup(0);
		double scalefactorphoup(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_error = sqrt(pow(s_ele_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
			scalefactorup = scalefactor + s_error; 
			scalefactorphoup = scalefactor + s_pho_error;
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      			double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			double s_error = sqrt(pow(s_mu_error,2) + pow(s_trg_error,2));
			scalefactorup = scalefactor + s_error; 
			scalefactorphoup = scalefactor + s_pho_error;
		}

		if(mcType == 4 && llmass < 30)continue;
		if(mcType == 5 && llmass > 30)continue;
		
		float XS_weight = 1;
		if(channelType == 1){
			if(RunYear == 2016 and preVFP == 1) 		XS_weight = lumi_2016preVFP_DoubleEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2016 and preVFP == 0) 	XS_weight = lumi_2016postVFP_DoubleEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2017) 			XS_weight = lumi_2017_DoubleEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2018) 			XS_weight = lumi_2018_DoubleEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 678) 			XS_weight = lumi_678_DoubleEG*1000*crosssection/ntotalevent;}

		else if(channelType == 2){
			if(RunYear == 2016 and preVFP == 1) 		XS_weight = lumi_2016preVFP_MuonEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2016 and preVFP == 0) 	XS_weight = lumi_2016postVFP_MuonEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2017) 			XS_weight = lumi_2017_MuonEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 2018) 			XS_weight = lumi_2018_MuonEG*1000*crosssection/ntotalevent;
        		else if(RunYear == 678) 			XS_weight = lumi_678_MuonEG*1000*crosssection/ntotalevent;}
		double reweightF(1);
		double Normalization=1;

if(channelType == 1){
	if(RunYear==2016 && preVFP==1){
		if(ISRJetPt < 50)reweightF = 1.02164;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 0.964486;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 1.10058;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.809766;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.795802;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.681508;
		else if(ISRJetPt >= 300)reweightF = 1.12425;
		Normalization = 0.989201;
	}
	else if(RunYear==2016 && preVFP==0){
		if(ISRJetPt < 50)reweightF = 0.92813;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 1.01507;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 1.03074;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.812674;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.905019;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.937368;
		else if(ISRJetPt >= 300)reweightF = 0.885853;
		Normalization = 1.06105;
	}
	else if(RunYear==2017){
		if(ISRJetPt < 50)reweightF = 0.80824;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 0.77143;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.891771;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.898236;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 1.12492;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.943156;
		else if(ISRJetPt >= 300)reweightF = 0.552922;
		Normalization = 1.23312;
	}
	else if(RunYear==2018){
		if(ISRJetPt < 50)reweightF = 0.883194;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 0.931125;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.949823;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.818615;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.854429;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 1.24516;
		else if(ISRJetPt >= 300)reweightF = 0.850772;
		Normalization = 1.12;
	}
}

else if(channelType == 2){
	if(RunYear==2016 && preVFP==1){
		if(ISRJetPt < 50)reweightF = 0.995888;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 1.08304;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 1.09747;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 1.08664;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 1.17574;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 0.806937;
		else if(ISRJetPt >= 300)reweightF = 1.12874;
		Normalization = 0.984708;
	}
	else if(RunYear==2016 && preVFP==0){
		if(ISRJetPt < 50)reweightF = 0.931823;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 1.01419;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 1.06259;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.847699;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.958803;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 1.24982;
		else if(ISRJetPt >= 300)reweightF = 0.900651;
		Normalization = 1.05421;
	}
	else if(RunYear==2017){
		if(ISRJetPt < 50)reweightF = 0.868712;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 0.874605;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.924282;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.976312;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 1.06787;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 1.22388;
		else if(ISRJetPt >= 300)reweightF = 0.676038;
		Normalization = 1.1413;
	}
	else if(RunYear==2018){
		if(ISRJetPt < 50)reweightF = 0.943713;
		else if(ISRJetPt >= 50 && ISRJetPt < 100)reweightF = 0.958948;
		else if(ISRJetPt >= 100 && ISRJetPt < 150)reweightF = 0.966778;
		else if(ISRJetPt >= 150 && ISRJetPt < 200)reweightF = 0.938759;
		else if(ISRJetPt >= 200 && ISRJetPt < 250)reweightF = 0.96244;
		else if(ISRJetPt >= 250 && ISRJetPt < 300)reweightF = 1.0375;
		else if(ISRJetPt >= 300)reweightF = 0.800688;
		Normalization = 1.05642;
	}
}


		ISRWeight = reweightF*Normalization;
		// as in the 2016 analysis: the data/Z+gamma ISR weights go to all V+gamma MC (W+gamma, Z+gamma, DY);
		// the LO -> NLO factor below is W+gamma only (mcType 1-3)
		const double nloFactor = (applyLOtoNLO && mcType >= 1 && mcType <= 3) ? LOtoNLOFactor(channelType, ISRJetPt) : 1.0;

		float weight = 1.0;
		float weight_scaleup = 1.0;
		float weight_scalephoup = 1.0;
		float weight_normup = 1.0;
		float weight_noisr = 1.0;
		float weight_nonlo = 1.0;
		float weight_puUp = 1.0;
		float weight_puDown = 1.0;

		if(!toDeriveScale and useVarSF){
			int htBin = hVGamma->FindBin(HT);
			factorMC   = hVGamma->GetBinContent(htBin);
			double factorErr = hVGamma->GetBinError(htBin);
			factorMCUP = factorMC + factorErr;
			// fall back to the fixed (!useVarSF) factors for an HT bin with a negative factor,
			// or with a relative error above maxRelErrVarSF (too few events to trust the bin)
			bool negFactor  = (factorMC < 0);
			bool largeError = (factorMC > 0 && factorErr/factorMC > maxRelErrVarSF);
			if(negFactor || largeError){
				static std::set<int> warnedHTbins;
				if(warnedHTbins.insert(htBin).second)
					cout << "WARNING: HT bin " << htBin << " factor = " << factorMC << " +- " << factorErr
					     << (negFactor ? " (negative)" : " (relative error > maxRelErrVarSF)")
					     << " -> using fixed factor " << factorMC_fixed << endl;
				factorMC   = factorMC_fixed;
				factorMCUP = factorMCUP_fixed;}}
		// photon-pT binned V+gamma factor (default; SF_BINNING=HT keeps the HT-binned factor above);
		// fixed factor for a photon-pT bin with a negative factor or a relative error above maxRelErrVarSF
		int kSF = -1;                       // photon-pT SF bin of this event (-1: no binned factor)
		double sfUpFac = 1., sfDownFac = 1.;  // +1 / -1 sigma factor relative to the nominal
		if(!toDeriveScale and usePhoPtSF){
			static TH1D *hVGpho = nullptr, *hVGphoDn = nullptr, *hVGphoUp = nullptr;
			if(!hVGpho){ TFile *fpho = TFile::Open(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/VGammaScaleHist_%s_%d%s.root",RunYear,whichVFP.c_str(),(channelType == 1 ? "eg" : "mg"),RunYear,whichVFP.c_str())); hVGpho = (TH1D*)fpho->Get("VGamma_pho");
				hVGphoDn = (TH1D*)fpho->Get("VGamma_pho_errDown"); hVGphoUp = (TH1D*)fpho->Get("VGamma_pho_errUp");
				if(!hVGphoDn || !hVGphoUp){ cout << "ERROR: VGamma_pho_errDown/Up missing, rerun plotVGammaScales.C" << endl; return; }}
			int phoBin = hVGpho->FindBin(std::min<double>(phoEt, hVGpho->GetXaxis()->GetXmax() - 1e-3));
			double sfv = hVGpho->GetBinContent(phoBin), sfe = hVGpho->GetBinError(phoBin);
			if(sfv > 0 && sfe/sfv <= maxRelErrVarSF){ factorMC = sfv; factorMCUP = sfv + sfe;
				kSF = phoBin - 1; sfUpFac = (sfv + hVGphoUp->GetBinContent(phoBin))/sfv; sfDownFac = (sfv - hVGphoDn->GetBinContent(phoBin))/sfv; }
			else { factorMC = factorMC_fixed; factorMCUP = factorMCUP_fixed; nNoSFbin++; }}

                weight = 		PUweight     * XS_weight * factorMC * scalefactor;
                weight_scaleup = 	PUweight     * XS_weight * factorMC * scalefactorup;
                weight_scalephoup = 	PUweight     * XS_weight * factorMC * scalefactorphoup;
                weight_normup = 	PUweight     * XS_weight * factorMCUP * scalefactor;
                weight_noisr = 		PUweight     * XS_weight * factorMC * scalefactor;
        	weight_puUp   = 	PUweightUp   * XS_weight * factorMC * scalefactor;
        	weight_puDown = 	PUweightDown * XS_weight * factorMC * scalefactor;
	        
		// weight_noisr: without the data/Z+gamma ISR weight (ISR uncertainty);
		// weight_nonlo: without the LO -> NLO W+gamma factor (W+gamma generator uncertainty)
		weight_nonlo = weight;
		if(applyISRwt){
			weight = weight*ISRWeight*nloFactor;
			weight_scaleup = weight_scaleup*ISRWeight*nloFactor;
			weight_scalephoup = weight_scalephoup*ISRWeight*nloFactor;
			weight_normup = weight_normup*ISRWeight*nloFactor;
			weight_puUp   = weight_puUp*ISRWeight*nloFactor;
			weight_puDown = weight_puDown*ISRWeight*nloFactor;
			weight_noisr  = weight_noisr*nloFactor;
			weight_nonlo  = weight_nonlo*ISRWeight;
		}

		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		// Mt and lepton pT cuts
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){	mindRpho=dR; phoIndex=iMC;	}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && (fabs((*mcMomPID)[phoIndex]) <= 6 || fabs((*mcMomPID)[phoIndex]) == 21 || fabs((*mcMomPID)[phoIndex]) == 999 || fabs((*mcMomPID)[phoIndex])== 11 || fabs((*mcMomPID)[phoIndex])== 13 || fabs((*mcMomPID)[phoIndex])== 15 || fabs((*mcMomPID)[phoIndex])== 23 || fabs((*mcMomPID)[phoIndex])== 24)  )istruepho=true;
		}
		if(!istruepho)continue;
		p_MET->Fill(sigMET, weight);
		jesup_MET->Fill(sigMETJESup, weight);
		jesdo_MET->Fill(sigMETJESdo, weight);
		jerup_MET->Fill(sigMETJERup, weight);
		jerdo_MET->Fill(sigMETJERdo, weight);
		scaleup_MET->Fill(sigMET, weight_scaleup);
		normup_MET->Fill(sigMET,   weight_normup);
		isrup_MET->Fill(sigMET, weight_noisr);
		
		puup_MET   ->Fill(sigMET, weight_puUp);
        	pudo_MET   ->Fill(sigMET, weight_puDown);
        	puup_HT    ->Fill(HT,    weight_puUp);
        	pudo_HT    ->Fill(HT,    weight_puDown);
        	puup_PhoEt ->Fill(phoEt, weight_puUp);
        	pudo_PhoEt ->Fill(phoEt, weight_puDown);
		puup_Mt   ->Fill(sigMT, weight_puUp);
        	pudo_Mt   ->Fill(sigMT, weight_puDown);

		float weight_pswIsrUp  = PUweight * XS_weight * scalefactor * factorMC * psWeightIsrUp;
        	float weight_pswIsrDo  = PUweight * XS_weight * scalefactor * factorMC * psWeightIsrDown;
        	pswIsrUp_MET  ->Fill(sigMET, weight_pswIsrUp);
        	pswIsrDo_MET  ->Fill(sigMET, weight_pswIsrDo);
        	pswIsrUp_Mt  ->Fill(sigMT, weight_pswIsrUp);
        	pswIsrDo_Mt  ->Fill(sigMT, weight_pswIsrDo);
        	pswIsrUp_HT   ->Fill(HT,    weight_pswIsrUp);
        	pswIsrDo_HT   ->Fill(HT,    weight_pswIsrDo);
        	pswIsrUp_PhoEt->Fill(phoEt, weight_pswIsrUp);
        	pswIsrDo_PhoEt->Fill(phoEt, weight_pswIsrDo);
		
		float weight_pswFsrUp  = PUweight * XS_weight * scalefactor * factorMC * psWeightFsrUp;
        	float weight_pswFsrDo  = PUweight * XS_weight * scalefactor * factorMC * psWeightFsrDown;
        	pswFsrUp_MET  ->Fill(sigMET, weight_pswFsrUp);
        	pswFsrDo_MET  ->Fill(sigMET, weight_pswFsrDo);
        	pswFsrUp_Mt  ->Fill(sigMT, weight_pswFsrUp);
        	pswFsrDo_Mt  ->Fill(sigMT, weight_pswFsrDo);
        	pswFsrUp_HT   ->Fill(HT,    weight_pswFsrUp);
        	pswFsrDo_HT   ->Fill(HT,    weight_pswFsrDo);
        	pswFsrUp_PhoEt->Fill(phoEt, weight_pswFsrUp);
        	pswFsrDo_PhoEt->Fill(phoEt, weight_pswFsrDo);
		/** cut flow *****/
		// MET cut
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta,weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta,weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
		p_nJet->Fill(nJet, weight);

		int SigBinIndex(-1);
		float m_lg = invmass;
		// fill events in signal bins defined in analysis_binning.h
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		//if(SigBinIndex >=0){
			h_VGamma_norm->Fill( SigBinIndex, weight);
			h_VGamma_normUp->Fill( SigBinIndex, weight_normup); 
			for(int k = 0; k < 3; k++){
				h_VGamma_sfUp[k]->Fill(SigBinIndex,   k == kSF ? weight*sfUpFac   : weight);
				h_VGamma_sfDown[k]->Fill(SigBinIndex, k == kSF ? weight*sfDownFac : weight);}
			h_VGamma_esfUp->Fill( SigBinIndex, weight_scaleup); 
			h_VGamma_esfphoUp->Fill( SigBinIndex, weight_scalephoup); 
			h_VGamma_isrAlter->Fill( SigBinIndex, weight_noisr);
			h_VGamma_nloAlter->Fill( SigBinIndex, weight_nonlo);
			h_VGamma_puUp      ->Fill(SigBinIndex, weight_puUp);
            		h_VGamma_puDown    ->Fill(SigBinIndex, weight_puDown);

			// QCD scale and PDF varied yields for bin-wise lnN
						const float wClampLo = 0.01f;
			const float wClampHi = 10.0f;
			auto clamp = [&](float v){ return std::max(wClampLo, std::min(wClampHi, v)); };
			float cMuRup   = clamp(muRweightUp);
			float cMuRdown = clamp(muRweightDown);
			float cMuFup   = clamp(muFweightUp);
			float cMuFdown = clamp(muFweightDown);
			float safePdfRMS = (std::isfinite(pdfWeightRMS) && pdfWeightRMS >= 0.f)
			                   ? std::min(pdfWeightRMS, 1.0f) : 0.f;

			h_VGamma_muRup  ->Fill(SigBinIndex, weight * cMuRup);
			h_VGamma_muRdown->Fill(SigBinIndex, weight * cMuRdown);
			h_VGamma_muFup  ->Fill(SigBinIndex, weight * cMuFup);
			h_VGamma_muFdown->Fill(SigBinIndex, weight * cMuFdown);

			// PDF varied yields
			h_VGamma_pdfUp  ->Fill(SigBinIndex, weight * (1.f + safePdfRMS));
			h_VGamma_pdfDown->Fill(SigBinIndex, weight * std::max(0.f, 1.f - safePdfRMS));

			h_VGamma_pswIsrUp  ->Fill(SigBinIndex, weight_pswIsrUp);
            		h_VGamma_pswIsrDown->Fill(SigBinIndex, weight_pswIsrDo);
			h_VGamma_pswFsrUp  ->Fill(SigBinIndex, weight_pswFsrUp);
            		h_VGamma_pswFsrDown->Fill(SigBinIndex, weight_pswFsrDo);
			                        // m(l,gamma) of the same SR events
                        h_VGamma_invmass_norm    ->Fill(m_lg, weight);
                        h_VGamma_invmass_normUp  ->Fill(m_lg, weight_normup);
                        h_VGamma_invmass_esfUp   ->Fill(m_lg, weight_scaleup);
                        h_VGamma_invmass_esfphoUp   ->Fill(m_lg, weight_scalephoup);
                        h_VGamma_invmass_muRup   ->Fill(m_lg, weight * cMuRup);
                        h_VGamma_invmass_muRdown ->Fill(m_lg, weight * cMuRdown);
                        h_VGamma_invmass_muFup   ->Fill(m_lg, weight * cMuFup);
                        h_VGamma_invmass_muFdown ->Fill(m_lg, weight * cMuFdown);
                        h_VGamma_invmass_pdfUp   ->Fill(m_lg, weight * (1.f + safePdfRMS));
                        h_VGamma_invmass_pdfDown ->Fill(m_lg, weight * std::max(0.f, 1.f - safePdfRMS));
		//}

		if(Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)>=0)h_VGamma_jesUp->Fill( Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt), weight );
		if(Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt)>=0)h_VGamma_jesDown->Fill( Bin.findSignalBin(sigMETJESdo, HTJESdo,phoEt),  weight);
		if(Bin.findSignalBin(sigMETJERup, HT, phoEt)>=0)h_VGamma_jerUp->Fill( Bin.findSignalBin(sigMETJERup, HT, phoEt),  weight);
		if(Bin.findSignalBin(sigMETJERdo, HT, phoEt)>=0)h_VGamma_jerDown->Fill( Bin.findSignalBin(sigMETJERdo, HT, phoEt), weight);

		// JES/JER only move events between signal bins (the lepton and photon, hence m(l,gamma), are unchanged)
		if(Bin.findSignalBin(sigMETJESup, HTJESup,  phoEt)>=0)h_VGamma_invmass_jesUp  ->Fill(m_lg, weight);
		if(Bin.findSignalBin(sigMETJESdo, HTJESdo,  phoEt)>=0)h_VGamma_invmass_jesDown->Fill(m_lg, weight);
		if(Bin.findSignalBin(sigMETJERup, HT,       phoEt)>=0)h_VGamma_invmass_jerUp  ->Fill(m_lg, weight);
		if(Bin.findSignalBin(sigMETJERdo, HT,       phoEt)>=0)h_VGamma_invmass_jerDown->Fill(m_lg, weight);

		jesup_Mt->Fill(sigMTJESup, weight);
		jesup_HT->Fill(HTJESup, weight);
		jesup_dPhiEleMET->Fill(fabs(dPhiLepMETJESup), weight);

		jesdo_Mt->Fill(sigMTJESdo, weight);
		jesdo_HT->Fill(HTJESdo, weight);
		jesdo_dPhiEleMET->Fill(fabs(dPhiLepMETJESdo), weight);

		jerup_Mt->Fill(sigMTJERup, weight);
		jerup_dPhiEleMET->Fill(fabs(dPhiLepMETJERup), weight);

		jerdo_Mt->Fill(sigMTJERdo, weight);
		jerdo_dPhiEleMET->Fill(fabs(dPhiLepMETJERdo), weight);

		scaleup_PhoEt->Fill(phoEt, weight_scaleup);
		scaleup_PhoEta->Fill(phoEta,weight_scaleup);
		scaleup_LepPt->Fill(lepPt, weight_scaleup);
		scaleup_LepEta->Fill(lepEta,weight_scaleup);
		scaleup_Mt->Fill(sigMT, weight_scaleup);
		scaleup_HT->Fill(HT, weight_scaleup);
		scaleup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaleup);

		normup_PhoEt->Fill(phoEt,  weight_normup);
		normup_PhoEta->Fill(phoEta,weight_normup);
		normup_LepPt->Fill(lepPt,  weight_normup);
		normup_LepEta->Fill(lepEta,weight_normup);
		normup_Mt->Fill(sigMT,     weight_normup);
		normup_HT->Fill(HT,        weight_normup);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_normup);

		isrup_PhoEt->Fill(phoEt, weight_noisr);
		isrup_PhoEta->Fill(phoEta,weight_noisr);
		isrup_LepPt->Fill(lepPt, weight_noisr);
		isrup_LepEta->Fill(lepEta,weight_noisr);
		isrup_Mt->Fill(sigMT, weight_noisr);
		isrup_HT->Fill(HT, weight_noisr);
		isrup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_noisr);
	}


	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((isrup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);    // ISR corrections
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((isrup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);    // ISR corrections
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((isrup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		double puerror = max(fabs(puup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(pudo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)));
        	syserror += pow(puerror, 2);
		double jeserror = max( fabs(jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)), fabs(jerdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((isrup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)), fabs(jerdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((isrup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)), fabs(jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_dPhiEleMET->GetBinError(ibin)* p_dPhiEleMET->GetBinError(ibin);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((normup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		syserror += pow((isrup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		double jeserror = max( fabs(jesup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jesdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		double jererror = max( fabs(jerup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)), fabs(jerdo_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)));
		syserror += pow(jeserror,2);
		syserror += pow(jererror,2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	
	float binSum=0;
	cout<< "bin nominalsig jesUp jesDown jerUp jerDown max(jesuperror, jesdoerror) max(jeruperror, jerdoerror) max_jes_jer  max_jes_jer/2.0"<<endl;
	for(int contbin(1); contbin <=NBIN; contbin++){
		binSum= binSum+h_VGamma_norm->GetBinContent(contbin);
		float nominalsig = h_VGamma_norm->GetBinContent(contbin); 

		float jesuperror = fabs(h_VGamma_jesUp->GetBinContent(contbin)- nominalsig);
		float jesdoerror = fabs(h_VGamma_jesDown->GetBinContent(contbin)- nominalsig);
		float jeruperror = fabs(h_VGamma_jerUp->GetBinContent(contbin)- nominalsig);
		float jerdoerror = fabs(h_VGamma_jerDown->GetBinContent(contbin)- nominalsig);
		float esfuperror = fabs(h_VGamma_esfUp->GetBinContent(contbin)- nominalsig); 	
		float esfphouperror = fabs(h_VGamma_esfphoUp->GetBinContent(contbin)- nominalsig); 	
		float normerror  = fabs(h_VGamma_normUp->GetBinContent(contbin)- nominalsig);		
		float isrerror   = fabs(h_VGamma_isrAlter->GetBinContent(contbin)- nominalsig);
		float wgnloerror = fabs(h_VGamma_nloAlter->GetBinContent(contbin)- nominalsig);
		float puerror     = max(fabs(h_VGamma_puUp  ->GetBinContent(contbin) - nominalsig), fabs(h_VGamma_puDown->GetBinContent(contbin) - nominalsig));

		h_VGamma_syserr_jes->SetBinContent(contbin, max(jesuperror,jesdoerror));
		h_VGamma_syserr_jer->SetBinContent(contbin, max(jeruperror,jerdoerror));
		
		// added by hand
    		if (std::max(jeruperror, jerdoerror) > nominalsig || std::max(jesuperror, jesdoerror) > nominalsig) {
        		// Correcting logic to set the bin content
        		float max_jes_jer = std::max({h_VGamma_jesUp->GetBinContent(contbin),
                                      h_VGamma_jesDown->GetBinContent(contbin),
                                      h_VGamma_jerUp->GetBinContent(contbin),
                                      h_VGamma_jerDown->GetBinContent(contbin)});
        		std::cout <<"Bin: " << contbin << "\tNominal: " << nominalsig << "\tjesUp: "
			<< h_VGamma_jesUp->GetBinContent(contbin) << "\tjesDown: "
			<< h_VGamma_jesDown->GetBinContent(contbin) << "\tjerUp: "
			<< h_VGamma_jerUp->GetBinContent(contbin) << "\tjerDown: "
			<< h_VGamma_jerDown->GetBinContent(contbin) << "\t"
                  	<< std::max(jesuperror, jesdoerror) << "\t"
                  	<< std::max(jeruperror, jerdoerror) << "\t"
                  	<< std::max(std::max(jesuperror, jesdoerror), std::max(jeruperror, jerdoerror)) << "\t" 
			<< max_jes_jer / 2.0 << std::endl;

        		h_VGamma_norm->SetBinContent(contbin, max_jes_jer / 2.0);
			nominalsig = h_VGamma_norm->GetBinContent(contbin); // update after set

        		// Recalculate errors
        		jesuperror = fabs(h_VGamma_jesUp->GetBinContent(contbin) - nominalsig);
        		jesdoerror = fabs(h_VGamma_jesDown->GetBinContent(contbin) - nominalsig);
        		jeruperror = fabs(h_VGamma_jerUp->GetBinContent(contbin) - nominalsig);
        		jerdoerror = fabs(h_VGamma_jerDown->GetBinContent(contbin) - nominalsig);
        		h_VGamma_syserr_jes->SetBinContent(contbin, std::max(jesuperror, jesdoerror));
        		h_VGamma_syserr_jer->SetBinContent(contbin, std::max(jeruperror, jerdoerror));
        		
			std::cout <<"Bin: " << contbin << "\tNominal: " << nominalsig << "\tjesUp: "
			<< h_VGamma_jesUp->GetBinContent(contbin) << "\tjesDown: "
			<< h_VGamma_jesDown->GetBinContent(contbin) << "\tjerUp: "
			<< h_VGamma_jerUp->GetBinContent(contbin) << "\tjerDown: "
			<< h_VGamma_jerDown->GetBinContent(contbin) << "\t"
                  	<< std::max(jesuperror, jesdoerror) << "\t"
                  	<< std::max(jeruperror, jerdoerror) << "\t"
                  	<< std::max(std::max(jesuperror, jesdoerror), std::max(jeruperror, jerdoerror)) << "\t" 
			<< max_jes_jer / 2.0 << std::endl;
		}

		h_VGamma_syserr_esf->SetBinContent(contbin, esfuperror);
		h_VGamma_syserr_esfpho->SetBinContent(contbin, esfphouperror);
		h_VGamma_syserr_scale->SetBinContent(contbin, normerror);
		h_VGamma_syserr_eleshape->SetBinContent(contbin, -1);
		h_VGamma_syserr_jetshape->SetBinContent(contbin, -1);
		h_VGamma_syserr_qcdshape->SetBinContent(contbin, -1);
		h_VGamma_syserr_xs->SetBinContent(contbin, -1);     
		h_VGamma_syserr_lumi->SetBinContent(contbin, -1);      
		h_VGamma_syserr_isr->SetBinContent(contbin, isrerror);
		h_VGamma_syserr_wgnlo->SetBinContent(contbin, wgnloerror);      
		h_VGamma_syserr_pu      ->SetBinContent(contbin, puerror);

		// ---- Bin-wise lnN kappas for QCD scale (muR, muF) and PDF ----
		if(nominalsig > 0){
			float kmuRup   = h_VGamma_muRup  ->GetBinContent(contbin) / nominalsig;
			float kmuRdown = h_VGamma_muRdown->GetBinContent(contbin) / nominalsig;
			float kmuFup   = h_VGamma_muFup  ->GetBinContent(contbin) / nominalsig;
			float kmuFdown = h_VGamma_muFdown->GetBinContent(contbin) / nominalsig;
			float kpdfUp   = h_VGamma_pdfUp  ->GetBinContent(contbin) / nominalsig;
			float kpdfDown = h_VGamma_pdfDown->GetBinContent(contbin) / nominalsig;
			h_VGamma_syserr_qcdmuR->SetBinContent(contbin,
				std::max(fabs(kmuRup-1.f), fabs(kmuRdown-1.f)));
			h_VGamma_syserr_qcdmuF->SetBinContent(contbin,
				std::max(fabs(kmuFup-1.f), fabs(kmuFdown-1.f)));
			h_VGamma_syserr_pdf   ->SetBinContent(contbin,
				std::max(fabs(kpdfUp-1.f), fabs(kpdfDown-1.f)));
		} else {
			h_VGamma_syserr_qcdmuR->SetBinContent(contbin, 0.f);
			h_VGamma_syserr_qcdmuF->SetBinContent(contbin, 0.f);
			h_VGamma_syserr_pdf   ->SetBinContent(contbin, 0.f);
		}
	}
	for(int contbin(1); contbin <=NBIN; contbin++){
		cout<<contbin <<"\t" << h_VGamma_norm->GetBinContent(contbin) << "\t" << h_VGamma_syserr_jes->GetBinContent(contbin) << "\t" << h_VGamma_syserr_jer->GetBinContent(contbin)<< endl;
	}

	// ---- per-mass-bin uncertainties: same recipe as the signal-bin loop above ----
	std::cout << "\n=== m(l,gamma) bins (VGamma): bin nominal jes jer esf scale muR muF pdf ===" << std::endl;
	for(int ibin(1); ibin <= nInvMassBins; ibin++){
		float nominalsig = h_VGamma_invmass_norm->GetBinContent(ibin);

		float jesuperror = fabs(h_VGamma_invmass_jesUp  ->GetBinContent(ibin) - nominalsig);
		float jesdoerror = fabs(h_VGamma_invmass_jesDown->GetBinContent(ibin) - nominalsig);
		float jeruperror = fabs(h_VGamma_invmass_jerUp  ->GetBinContent(ibin) - nominalsig);
		float jerdoerror = fabs(h_VGamma_invmass_jerDown->GetBinContent(ibin) - nominalsig);
		// ESF and norm are taken w.r.t. the raw nominal yield, as in the signal-bin loop
		float esfuperror = fabs(h_VGamma_invmass_esfUp  ->GetBinContent(ibin) - nominalsig);
		float esfphouperror = fabs(h_VGamma_invmass_esfphoUp  ->GetBinContent(ibin) - nominalsig);
		float normerror  = fabs(h_VGamma_invmass_normUp ->GetBinContent(ibin) - nominalsig);

		if(applyLowStatJESJERFix && (std::max(jeruperror, jerdoerror) > nominalsig || std::max(jesuperror, jesdoerror) > nominalsig)){
			float max_jes_jer = std::max({h_VGamma_invmass_jesUp->GetBinContent(ibin),
			                              h_VGamma_invmass_jesDown->GetBinContent(ibin),
			                              h_VGamma_invmass_jerUp->GetBinContent(ibin),
			                              h_VGamma_invmass_jerDown->GetBinContent(ibin)});
			std::cout << "mass bin " << ibin << ": JES/JER variation > nominal (" << nominalsig
			          << ") -> nominal set to " << max_jes_jer/2.0 << std::endl;
			h_VGamma_invmass_norm->SetBinContent(ibin, max_jes_jer / 2.0);
			nominalsig = h_VGamma_invmass_norm->GetBinContent(ibin);
			jesuperror = fabs(h_VGamma_invmass_jesUp  ->GetBinContent(ibin) - nominalsig);
			jesdoerror = fabs(h_VGamma_invmass_jesDown->GetBinContent(ibin) - nominalsig);
			jeruperror = fabs(h_VGamma_invmass_jerUp  ->GetBinContent(ibin) - nominalsig);
			jerdoerror = fabs(h_VGamma_invmass_jerDown->GetBinContent(ibin) - nominalsig);
		}
		h_VGamma_invmass_syserr_jes  ->SetBinContent(ibin, std::max(jesuperror, jesdoerror));
		h_VGamma_invmass_syserr_jer  ->SetBinContent(ibin, std::max(jeruperror, jerdoerror));
		h_VGamma_invmass_syserr_esf  ->SetBinContent(ibin, esfuperror);
		h_VGamma_invmass_syserr_esfpho  ->SetBinContent(ibin, esfphouperror);
		h_VGamma_invmass_syserr_scale->SetBinContent(ibin, normerror);

		float fmuR(0.f), fmuF(0.f), fpdf(0.f);
		if(nominalsig > 0){
			fmuR = std::max(fabs(h_VGamma_invmass_muRup->GetBinContent(ibin)/nominalsig - 1.f), fabs(h_VGamma_invmass_muRdown->GetBinContent(ibin)/nominalsig - 1.f));
			fmuF = std::max(fabs(h_VGamma_invmass_muFup->GetBinContent(ibin)/nominalsig - 1.f), fabs(h_VGamma_invmass_muFdown->GetBinContent(ibin)/nominalsig - 1.f));
			fpdf = std::max(fabs(h_VGamma_invmass_pdfUp->GetBinContent(ibin)/nominalsig - 1.f), fabs(h_VGamma_invmass_pdfDown->GetBinContent(ibin)/nominalsig - 1.f));
		}
		h_VGamma_invmass_syserr_qcdmuR->SetBinContent(ibin, fmuR);
		h_VGamma_invmass_syserr_qcdmuF->SetBinContent(ibin, fmuF);
		h_VGamma_invmass_syserr_pdf   ->SetBinContent(ibin, fpdf);
		std::cout << ibin << "\t" << h_VGamma_invmass_norm->GetBinContent(ibin) << "\t" << std::max(jesuperror, jesdoerror) << "\t" << std::max(jeruperror, jerdoerror)
		          << "\t" << esfuperror << "\t" << normerror << "\t" << fmuR << "\t" << fmuF << "\t" << fpdf << std::endl;
	}

	// ---- Bin-wise lnN nuisance table for Combine datacards ----
	std::cout << "\n=== Bin-wise lnN nuisance parameters (VGamma) ===" << std::endl;
	std::cout << "Bin\tNominal\tkmuRup\tkmuRdn\tkmuFup\tkmuFdn\tkpdfUp\tkpdfDn" << std::endl;
	for(int contbin(1); contbin <=NBIN; contbin++){
		float nom = h_VGamma_norm->GetBinContent(contbin);
		if(nom <= 0){ cout << contbin << "\t0\t-\t-\t-\t-\t-\t-" << endl; continue; }
		float kmuRup   = h_VGamma_muRup  ->GetBinContent(contbin) / nom;
		float kmuRdown = h_VGamma_muRdown->GetBinContent(contbin) / nom;
		float kmuFup   = h_VGamma_muFup  ->GetBinContent(contbin) / nom;
		float kmuFdown = h_VGamma_muFdown->GetBinContent(contbin) / nom;
		float kpdfUp   = h_VGamma_pdfUp  ->GetBinContent(contbin) / nom;
		float kpdfDown = h_VGamma_pdfDown->GetBinContent(contbin) / nom;
		cout << std::fixed << std::setprecision(4)
		     << contbin << "\t" << nom
		     << "\t" << kmuRup   << "\t" << kmuRdown
		     << "\t" << kmuFup   << "\t" << kmuFdown
		     << "\t" << kpdfUp   << "\t" << kpdfDown << endl;
	}
	std::cout << "\nDatacard lines (process=VGamma, one lnN per bin):" << std::endl;
	std::cout << "  QCDscale_muR_VGamma  lnN  kmuRup/kmuRdn" << std::endl;
	std::cout << "  QCDscale_muF_VGamma  lnN  kmuFup/kmuFdn" << std::endl;
	std::cout << "  pdf_VGamma           lnN  kpdfUp/kpdfDn" << std::endl;
	cout<<"total entries : "<<binSum<<endl;
	if(nNoSFbin) cout << "WARNING: " << nNoSFbin << " events with the fixed V+gamma factor: no normalisation nuisance in the cards" << endl;
	outputfile->Write();
	outputfile->Close();
}
