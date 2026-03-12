#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "RooPolynomial.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooKeysPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
//#include "../../../include/tdrstyle.C"
//#include "../include/RooCBExGaussShape.h"
#define NTOY 1000
using namespace RooFit;

int RunYear = 2016; bool ISpreVFP = true;
bool doEB = true;
const char*processName ="DY";
//const char*processName ="Data";

enum BinType{
  byPt = 0,
  byEta = 1,
  byVtx = 2,
  byHT  = 3 
};

bool isElectron(int PID, int momID){
   bool isEle;
   if(fabs(PID) == 11){ 
	   switch(momID){
	     case 1: isEle = true; break;
	     case 2: isEle = true; break;
	     case 3: isEle = true; break;
	     case 4: isEle = true; break;
	     case 5: isEle = true; break;
	     case 6: isEle = true; break;
	     case 21: isEle = true; break;
	     case 23: isEle = true; break;
	     default: isEle = false; break;
	   }
  }
  else isEle = false;

  return isEle;
}

void FitKer(int fitFunc, int inputbintype, int inputfittype, float lowercut, float uppercut, int fitrangelow, int fitrangehigh){
   gSystem->Load("../../../lib/libAnaClasses.so");
   gSystem->Load("../../../lib/libRooFitClasses.so");
	
	std::string whichVFP;
        if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
	//setTDRStyle();
        //bool useCMSShape;
	bool useExpo;
	bool useKer;
	bool useDY;
	bool useCrystalBall;
	bool usePoly;

	//bool SaveOutput;
        ifstream configfile;
        if(fitFunc==1)
		configfile.open ("configFit_Bw_ker.txt");
        if(fitFunc==2)
		configfile.open ("configFit_DY_ker.txt");
        if(fitFunc==3)
		configfile.open ("configFit_Bw_expo.txt");
        if(fitFunc==4)
		configfile.open ("configFit_cb_poly.txt");

	std::string variabletype;
	bool variablevalue;
	if(configfile.is_open()){
  	for(int i(0); i<=4; i++){ 
			configfile >> variabletype >> variablevalue; 
			if(variabletype.find("Expo")!=std::string::npos)useExpo = variablevalue;
			else if(variabletype.find("Ker")!=std::string::npos)useKer = variablevalue;
			else if(variabletype.find("DY")!=std::string::npos)useDY = variablevalue;
			else if(variabletype.find("Poly")!=std::string::npos)
						{usePoly = variablevalue; useCrystalBall = variablevalue;}
	  }
	}
	configfile.close();
 
  BinType bintype = (BinType)inputbintype;

	char lowername[5];
	if(bintype == byEta)
		sprintf(lowername, "%d.%02d", (int)lowercut, (int)(((lowercut - (int)lowercut) * 100) + 0.5));
	else sprintf(lowername, "%d", (int)lowercut);	
  char uppername[5];
	if(uppercut < 1000){
		if(bintype == byEta)sprintf(uppername, "%d.%02d", (int)uppercut, (int)((uppercut-(int)uppercut)*100));	
		else sprintf(uppername, "%d", (int)uppercut);
	}
	else sprintf(uppername, "Inf");	
  std::ostringstream histname;      
	TH1F *DY_PU  = new TH1F("DY_PU", "",100,0,100);

//************** Process Z->ee Tree ********************************************************
  // Define RooRealVars
RooRealVar mass("invmass", "Invariant Mass", fitrangelow, fitrangehigh);
RooRealVar pt("probePt", "Probe pT", 0, 500);
RooRealVar eta("probeEta", "Probe Eta", -3.0, 3.0);
RooRealVar vtx("nVertex", "nVertex", 0, 100);
RooRealVar veto("vetovalue", "Veto", 0, 2);
RooRealVar fsr("FSRveto", "FSR Veto", 0, 2);
RooRealVar ht("HT", "HT", 0, 3000);

TChain *etree = new TChain("FakeRateTree");
etree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/plot_elefakepho_%sTnP_dR05_%d%s.root",
                processName, RunYear, whichVFP.c_str()));

RooArgSet varset(mass, pt, eta, vtx, ht, veto, fsr);

std::ostringstream cut;
//cut << "(HT >= 0 && HT < 100 &&";  // check diff HT ranges
//cut << "(HT >= 100 && HT < 400 &&";  // check diff HT ranges
//cut << "(HT >= 400 &&";  // check diff HT ranges
cut << "(";

if (doEB)
    cut << "fabs(probeEta)<=1.4442";
else
    cut << "fabs(probeEta)>=1.56";

cut << " && ";
if (bintype == byPt)
    cut << "probePt>=" << lowercut << " && probePt<" << uppercut;
else if (bintype == byEta)
    cut << "fabs(probeEta)>=" << lowercut << " && fabs(probeEta)<" << uppercut;
else if (bintype == byVtx)
    cut << "nVertex>=" << lowercut << " && nVertex<" << uppercut;
else if (bintype == byHT)
    cut << "HT>=" << lowercut << " && HT<" << uppercut;

cut << " && ";
if (inputfittype == 0)
    cut << "vetovalue==0";
else
    cut << "vetovalue==1 && FSRveto==1";

cut << ")";

RooDataSet data("data", "Filtered dataset", etree, varset, cut.str().c_str());

if (data.numEntries() < 10) {
    std::cout << "WARNING: Too few entries in dataset (" << data.numEntries() << "). Skipping bin." << std::endl;
    return;
}


//************** Process mg Bkg Tree *******************************************
RooRealVar mass_bg("invmass", "Invariant Mass", 60, 120);
RooRealVar pt_bg("probePt", "probePt", 0, 300);
RooRealVar eta_bg("probeEta", "probeEta", -3.0, 3.0);
RooRealVar veto_bg("vetovalue", "vetovalue", 0, 2);
RooRealVar fsr_bg("FSRveto", "FSRveto", 0, 2);
RooRealVar nvtx_bg("nVertex", "nVertex", 0, 100);
RooRealVar ht_bg("HT", "HT", 0, 3000);

RooArgSet vars_bg(mass_bg, pt_bg, eta_bg, veto_bg, fsr_bg, nvtx_bg, ht_bg);

TChain* bgtree = new TChain("BGTree");
bgtree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/plot_bgtemplate_FullEcal_%d%s.root", RunYear, whichVFP.c_str()));

std::stringstream cut_bg;
//cut_bg << "HT >= 0 && HT < 100 &&";  // check diff HT ranges

if (doEB)
    cut_bg << "fabs(probeEta) <= 1.4442";
else
    cut_bg << "fabs(probeEta) >= 1.56";

cut_bg << " && ";

if (bintype == byPt)
    cut_bg << Form("probePt >= %.1f && probePt < %.1f", lowercut, uppercut);
else if (bintype == byEta)
    cut_bg << Form("fabs(probeEta) >= %.2f && fabs(probeEta) < %.2f", lowercut, uppercut);
else if (bintype == byVtx)
    cut_bg << Form("nVertex >= %.1f && nVertex < %.1f", lowercut, uppercut);
else if (bintype == byHT)
    cut_bg << Form("HT >= %.1f && HT < %.1f", lowercut, uppercut);

cut_bg << " && ";

if (inputfittype == 0)
    cut_bg << "vetovalue == 0";
else if (inputfittype == 1)
    cut_bg << "vetovalue == 1 && FSRveto == 1";

std::cout << "Background cut applied: " << cut_bg.str() << std::endl;
RooDataSet data_bg("data_bg", "Filtered BG dataset", bgtree, vars_bg, cut_bg.str().c_str());
if (data_bg.numEntries() < 10) {
    std::cout << "Not enough background entries in bin: " << lowercut << " - " << uppercut << std::endl;
    return;
}


//************** Process DYJet Tree *******************************************//
  TTree *newDYtree = new TTree("newDYtree","newDYtree");
  double invmass_DYsignal(0);
  newDYtree->Branch("invmass", &invmass_DYsignal);
	TH1F  *h_DYinvmass = new TH1F("h_DYinvmass","h_DYinvmass",(int)(fitrangehigh-fitrangelow),fitrangelow,fitrangehigh);

	if(useDY){
		TChain *DYtree = new TChain("FakeRateTree");
		DYtree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/plot_elefakepho_DYTnP_dR05_%d%s.root",RunYear,whichVFP.c_str()));

		float DY_invmass=0; 
		float DY_tagPt=0; 
		float DY_tagEta=0; 
		float DY_tagPhi=0; 
		float DY_probePt=0; 
		float DY_probeEta=0;
		float DY_probePhi=0;
		bool  DY_vetovalue=0;
		int   DY_nVertex=0;
		float DY_HT = 0;
		std::vector<int>   *mcPID=0;
		std::vector<float> *mcEta=0;
		std::vector<float> *mcPhi=0;
		std::vector<float> *mcPt=0;
		std::vector<int> *mcMomPID=0;
		std::vector<int> *mcGMomPID=0;
		DYtree->SetBranchAddress("invmassUncalib",   &DY_invmass); 
		DYtree->SetBranchAddress("tagPt",     &DY_tagPt);
		DYtree->SetBranchAddress("tagEta",    &DY_tagEta);
		DYtree->SetBranchAddress("tagPhi",    &DY_tagPhi);
		DYtree->SetBranchAddress("probePt",   &DY_probePt);
		DYtree->SetBranchAddress("probeEta",  &DY_probeEta);
		DYtree->SetBranchAddress("probePhi",  &DY_probePhi);
		DYtree->SetBranchAddress("vetovalue", &DY_vetovalue);
		DYtree->SetBranchAddress("nVertex",   &DY_nVertex);
		DYtree->SetBranchAddress("mcPID",			&mcPID);
		DYtree->SetBranchAddress("mcEta",			&mcEta);
		DYtree->SetBranchAddress("mcPhi",			&mcPhi);
		DYtree->SetBranchAddress("mcPt",			&mcPt);
		DYtree->SetBranchAddress("mcMomPID",	&mcMomPID);
		DYtree->SetBranchAddress("mcGMomPID",	&mcGMomPID);
		DYtree->SetBranchAddress("HT", &DY_HT);

		for(unsigned iEvt(0); iEvt < DYtree->GetEntries(); iEvt++){
			DYtree->GetEntry(iEvt);

			bool isZee(false);
			double mindRtag(0.3), mindRprobe(0.3);
			unsigned tagIndex(0), probeIndex(0);
			for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
				double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_tagEta, DY_tagPhi);
				double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_probeEta,DY_probePhi);
				double dE1 = fabs((*mcPt)[iMC] - DY_tagPt)/DY_tagPt;
				double dE2 = fabs((*mcPt)[iMC] - DY_probePt)/DY_probePt;
				if(dR1 < mindRtag && dE1 < 0.1){mindRtag=dR1; tagIndex=iMC;}
				if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
			}
			if(mindRtag < 0.1 && mindRprobe < 0.1){
				bool isZe(false),isZg(false);
				isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
				isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
				if(isZe && isZg)isZee=true; 
			}
			float keyvariable(0);
			if(bintype == byPt)keyvariable = DY_probePt;
			else if(bintype == byEta)keyvariable = fabs(DY_probeEta);
			else if(bintype == byVtx)keyvariable = DY_nVertex;
			if(keyvariable >= lowercut && keyvariable < uppercut){
				double mcPUweight;
				if(RunYear==2016 and ISpreVFP==true) mcPUweight = getPUESF16preVFP(DY_nVertex);
				else if(RunYear==2016 and ISpreVFP==false) mcPUweight = getPUESF16(DY_nVertex);
				else if(RunYear==2017) mcPUweight = getPUESF17(DY_nVertex);
				else if(RunYear==2018) mcPUweight = getPUESF18(DY_nVertex);
				DY_PU->Fill(DY_nVertex, mcPUweight);
				if(isZee){invmass_DYsignal = DY_invmass; h_DYinvmass->Fill(DY_invmass, mcPUweight); }
				if(isZee)newDYtree->Fill();		
			}
		}
		h_DYinvmass->Sumw2();
	} // Run only useDY

  //************* Construct RooFit Models *********************************************************//
  
  RooRealVar mass_axis("invmass", "M_{tag-probe}", fitrangelow, fitrangehigh);
  mass_axis.setBins(2048, "cache");

  RooDataSet BkgDataSet("BkgDataSet","BkgDataSet", bgtree, RooArgSet(mass_axis, pt_bg, eta_bg, veto_bg, fsr_bg, nvtx_bg, ht_bg), cut_bg.str().c_str());
  RooKeysPdf BkgKer("BkgKer","BkgKer",mass_axis,BkgDataSet,RooKeysPdf::MirrorBoth, 2);
	RooDataSet *DYDataSet;
	RooDataHist *datahist_DY;
	RooKeysPdf *DYKer;
	RooHistPdf *DYpdf;

	if(useDY && h_DYinvmass->GetEntries() < 50)DYDataSet = new RooDataSet("DYDataSet","DYDataSet",RooArgSet(mass_axis), RooFit::Import(*newDYtree));
	else if(useDY && h_DYinvmass->GetEntries() >= 50)datahist_DY = new RooDataHist("DYDataSet","DYDataSet", mass_axis, h_DYinvmass);

	if(useDY && h_DYinvmass->GetEntries() < 50)DYKer = new RooKeysPdf("DYKer","DYKer",mass_axis, *DYDataSet,RooKeysPdf::MirrorBoth,2);
 	else if(useDY && h_DYinvmass->GetEntries() >= 50){
		DYpdf = new RooHistPdf("DYpdf","DYpdf", mass_axis, *datahist_DY); 
		DYpdf->setInterpolationOrder(1); 
	}

  TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
	c_fitMass->SetLeftMargin(0.15);
	c_fitMass->SetBottomMargin(0.15);

  c_fitMass->cd();

    RooRealVar cb_mean("cb_mean", "mean", 91.2, 85, 97);
    RooRealVar cb_sigma("cb_sigma", "sigma", 2.0, 0.5, 10.0);
    RooRealVar cb_alpha("cb_alpha", "alpha", 1.5, 0.5, 5.0);
    RooRealVar cb_n("cb_n", "n", 5, 0.5, 20);
    RooCBShape cb_sig("cb_sig", "Crystal Ball", mass_axis, cb_mean, cb_sigma, cb_alpha, cb_n);

    // Background model
    RooRealVar c0("c0", "c0", 0.1, -10, 10);
    RooRealVar c1("c1", "c1", -0.1, -10, 10);
    RooRealVar c2("c2", "c2", 0.1, -10, 10);
    RooPolynomial poly3("poly3", "Poly3", mass_axis, RooArgList(c0, c1, c2));

    TH1F* p_invmass = new TH1F("p_invmass", "p_invmass", 60, 60, 120);
	for (int i = 0; i < data.numEntries(); ++i) {
    		const RooArgSet* row = data.get(i);
    		p_invmass->Fill(row->getRealValue("invmass"));
	}

	double slopeupper = (p_invmass->GetBinContent(31) - p_invmass->GetBinContent(1))/31.0;
	double slopelower = (p_invmass->GetBinContent((int)(fitrangehigh-fitrangelow)) - p_invmass->GetBinContent(30))/30.0;
  double inislope = (p_invmass->GetBinContent(50) - p_invmass->GetBinContent(10))/40.0;
  RooRealVar slope("slope","slope",inislope,slopelower,slopeupper);
  double inilambda = (log(p_invmass->GetBinContent(1)) - log(p_invmass->GetBinContent(40)))/(p_invmass->GetBinCenter(1) - p_invmass->GetBinCenter(40));
  RooRealVar lambda("lambda", "slope", inilambda, -10., 10.);
  RooExponential *expo = new RooExponential("expo", "exponential PDF", mass_axis, lambda);                     // Exponential pdf.

  //RooRealVar m0( "m0", "m0", 91.188);
  RooRealVar m0( "m0", "m0", 91.188, 80,100);
  RooRealVar width( "width", "width", 2.495, 0.0, 15);
  RooRealVar mean("mean", "" ,0.,-3,3);
  RooRealVar sigma("sigma", "",2.4 , 0.0, 15.0);
  RooBreitWigner bw("bw", "", mass_axis, m0, width);               // Breit Wigner
  
  RooRealVar alpha("alpha", "", 1.0, 0.0, 20.0);
  RooRealVar n("n","", 10.0, 0.0, 100.0);
  RooRealVar alpha2("2ndalpha","", 1.0, 0.0, 20.0);
  RooRealVar n2("2ndn", "", 1.0, 0.0, 20.0);
  RooDCBShape *cb;                                                 // Double-sided crystal ball
  cb = new RooDCBShape("cb","cb", mass_axis, mean, sigma, alpha, n, alpha2, n2);

  RooGaussian gauss("gs", "gs", mass_axis, mean, sigma);           // gaussian


  RooAbsPdf* signalRes = nullptr;
  if(useCrystalBall) signalRes = &cb_sig;
  else if(!useDY)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, bw, *cb);                                                  // CONVOLUTION of bw & cb
  else if(useDY && h_DYinvmass->GetEntries() < 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYKer, gauss);            // CONVOLUTION of DY and gauss
  else if(useDY && h_DYinvmass->GetEntries() >= 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYpdf, gauss);     // when events are less use Kernel estimator

  int bin80GeV = p_invmass->FindBin(80);
	int bin100GeV= p_invmass->FindBin(100);
	int bin60GeV = p_invmass->FindBin(60);
	int bin120GeV= p_invmass->FindBin(120);
  double iniSig = p_invmass->Integral(bin80GeV,bin100GeV);
  double iniBkg = p_invmass->Integral(1,bin80GeV) + p_invmass->Integral(bin100GeV,p_invmass->GetSize()); 
	if(iniBkg > p_invmass->GetEntries()/2)iniBkg = p_invmass->GetEntries()/2;
  RooRealVar nSig("nSig", "", 0.5*iniSig, 0, p_invmass->GetEntries()*1.2);
  RooRealVar nBkg("nBkg", "", iniBkg, p_invmass->Integral(1,bin60GeV)+p_invmass->Integral(bin120GeV,p_invmass->GetSize()), p_invmass->GetEntries());
  RooAddPdf *model;

  if(usePoly) model = new RooAddPdf("model", "cb+poly3", RooArgList(cb_sig, poly3), RooArgList(nSig, nBkg));
  else if(useExpo)model = new RooAddPdf("model", "", RooArgList(*expo, *signalRes),RooArgList(nBkg, nSig));                             // use expo background model
  else if(useKer)model = new RooAddPdf("model", "", RooArgList(BkgKer, *signalRes),RooArgList(nBkg, nSig));                        // use μ + prob background model

  histname.str("");
  if(inputfittype == 0) histname << "denominator ";
  else if(inputfittype == 1) histname << "numerator ";
	histname << lowername << " < ";
  if(inputbintype == 0)histname << "pt <";
  else if(inputbintype == 1)histname << "eta <";
  else if(inputbintype == 2)histname << "nvtx <";
  else if(inputbintype == 3) histname << "ht <";

  histname << uppername;
  RooPlot* mass_Frame = mass_axis.frame(RooFit::Title(histname.str().c_str()),RooFit::Bins(60));
  mass_Frame->SetStats(0);
  
  	model->fitTo(data, RooFit::Extended(), RooFit::Save());
	data.plotOn(mass_Frame);

  if(usePoly) model->plotOn(mass_Frame,
                		RooFit::Components(poly3),
                		RooFit::LineStyle(kDashed),
                		RooFit::LineColor(kBlue),
                		RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useExpo)model->plotOn(mass_Frame, RooFit::Components(*expo),                                ///  plot background
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame, RooFit::Components(BkgKer),
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));

  if (usePoly)  model->plotOn(mass_Frame,
                		RooFit::Components(RooArgSet(cb_sig, poly3)),
                		RooFit::LineStyle(kSolid),
                		RooFit::LineColor(kBlue),
                		RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));

  else if(useExpo)model->plotOn(mass_Frame,                                                           ///  plot signal + background
				 RooFit::Components(RooArgSet(*expo, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame,
				 RooFit::Components(RooArgSet(BkgKer, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  mass_Frame->SetYTitle("Events");
  mass_Frame->SetXTitle("M_{tag-probe}");
  mass_Frame->SetTitleSize(0.05,"X");
  mass_Frame->SetTitleOffset(1.1,"X");
  mass_Frame->Draw();

  double chi2 = mass_Frame->chiSquare();
  TLatex* latex = new TLatex();
  char chi2str[50];
  sprintf (chi2str, "chi2/ndf = %.02f", chi2);
  latex->DrawLatex(95,0.7*p_invmass->GetBinContent(31),chi2str);
  histname.str("");
	histname << processName << "_";
	if(useDY)histname << "DY_";                //        ------------  All four distributions
	else histname << "Bw_";
  if(useKer)histname << "ker_";
  else if(useExpo)histname << "expo_";
  if(inputbintype == 0)histname << "pt_";
  else if(inputbintype == 1)histname << "eta_";
  else if(inputbintype == 2)histname << "vtx_";
  else if(inputbintype == 3) histname << "ht_";

  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
	histname << fitrangelow << "-" << fitrangehigh << "_"; 
  histname << lowername << "_" << uppername <<".png";
  c_fitMass->SaveAs(histname.str().c_str());
  mass_axis.setRange("signal",80,101);
  RooAbsReal* igx_sig;
	igx_sig = signalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	double norminalmean = igx_sig->getVal()*(nSig.getVal());
	double norminalrms  = igx_sig->getVal()*(nSig.getError());

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  RooDataHist* dhist = new RooDataHist("dhist", "binned version", RooArgSet(mass_axis), data);
  RooHistPdf histpdf_data("histpdf_data", "histpdf_data", RooArgSet(mass_axis), *dhist);
  RooMCStudy* mcstudy = new RooMCStudy(histpdf_data,mass_axis,RooFit::FitModel(*model),RooFit::Binned(kTRUE),RooFit::Silence(),
                        RooFit::Extended(),RooFit::FitOptions(RooFit::Save(kTRUE),RooFit::PrintEvalErrors(0)));
  mcstudy->generateAndFit(NTOY,p_invmass->Integral(1, (int)(fitrangehigh-fitrangelow)));
  mcstudy->fitParDataSet().Print("v");
	std::vector<double> toyresult;
	toyresult.clear();
for(int i(0); i<mcstudy->fitParDataSet().sumEntries(); i++){
	double toynSig = mcstudy->fitParDataSet().get(i)->getRealValue("nSig",0,kFALSE);
        cout<<"toynSig : "<<toynSig<<endl;
		if(useDY){
			double v_toymean = mcstudy->fitParDataSet().get(i)->getRealValue("mean");
			double v_toysigma = mcstudy->fitParDataSet().get(i)->getRealValue("sigma");
			RooRealVar toymean("toymean", "" ,v_toymean);
			RooRealVar toysigma("toysigma", "", v_toysigma);
			RooGaussian toygauss("toygs", "", mass_axis, toymean, toysigma);
			RooFFTConvPdf *toysignalRes;
			if(h_DYinvmass->GetEntries() < 50)toysignalRes = new RooFFTConvPdf("toypdf", "",mass_axis, *DYKer, toygauss);
			else if(h_DYinvmass->GetEntries() >= 50)toysignalRes = new RooFFTConvPdf("toypdf", "",mass_axis, *DYpdf, toygauss);
			RooAbsReal* toy_sig = toysignalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
			toyresult.push_back(toynSig*toy_sig->getVal());
		}
		else{
			double v_toym0 = mcstudy->fitParDataSet().get(i)->getRealValue("m0");
			double v_toywidth = mcstudy->fitParDataSet().get(i)->getRealValue("width");
			double v_toymean = mcstudy->fitParDataSet().get(i)->getRealValue("mean");
			double v_toysigma = mcstudy->fitParDataSet().get(i)->getRealValue("sigma");
			double v_toyalpha = mcstudy->fitParDataSet().get(i)->getRealValue("alpha");
			double v_toyn = mcstudy->fitParDataSet().get(i)->getRealValue("n");
			double v_toyalpha2 = mcstudy->fitParDataSet().get(i)->getRealValue("2ndalpha");
			double v_toyn2 = mcstudy->fitParDataSet().get(i)->getRealValue("2ndn");
			RooRealVar toym0( "toym0", "", v_toym0);
			RooRealVar toywidth( "toywidth", "", v_toywidth);
			RooRealVar toymean("toymean", "" ,v_toymean );
			RooRealVar toysigma("toysigma", "", v_toysigma);
			RooRealVar toyalpha("toyalpha", "", v_toyalpha);
			RooRealVar toyn("toyn","", v_toyn);
			RooRealVar toyalpha2("toy2ndalpha","", v_toyalpha2);
			RooRealVar toyn2("toy2ndn", "", v_toyn2);
			RooBreitWigner toybw("toybw", "", mass_axis, toym0, toywidth);
			RooDCBShape *toycb;
			toycb = new RooDCBShape("toycb","", mass_axis, toymean, toysigma, toyalpha, toyn, toyalpha2, toyn2);
			RooFFTConvPdf *toysignalRes = new RooFFTConvPdf("toypdf", "",mass_axis, toybw, *toycb);
			RooAbsReal* toy_sig = toysignalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
			toyresult.push_back(toynSig*toy_sig->getVal());
		}

	}

	std::sort(toyresult.begin(), toyresult.end());
	TH1F *prob = new TH1F("nSigHist","",100, toyresult[0], toyresult[999]);
	for(unsigned ii(0); ii < toyresult.size(); ii++)prob->Fill(toyresult[ii]);
  TCanvas *toymc = new TCanvas("toymc","",600,600);
  toymc->cd();
//  p_chi2->Draw();
  prob->Draw();
  prob->Fit("gaus");
  double fitmean  =(prob->GetFunction("gaus")->GetParameter(1));
  double fitrms   =(prob->GetFunction("gaus")->GetParameter(2));

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  histname.str("");
	histname << processName << "_";
  if(useKer)histname << "ToyMC_ker_";
  else if(useExpo)histname << "ToyMC_expo_";
  if(inputbintype == 0)histname << "pt_";
  else if(inputbintype == 1)histname << "eta_";
  else if(inputbintype == 2)histname << "vtx_";
  else if(inputbintype == 3) histname << "ht_";

  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
  histname << lowername << "_" << uppername << ".png";
  toymc->SaveAs(histname.str().c_str());

  histname.str("");
  if(inputbintype == 0)histname << "pt ";
  else if(inputbintype == 1)histname << "eta ";
  else if(inputbintype == 2)histname << "vtx ";
  else if(inputbintype == 3) histname << "ht ";

  if(inputfittype == 0) histname << "den ";
  else if(inputfittype == 1) histname << "num ";
  histname << lowername << " " <<  norminalmean << " " << norminalrms << " " << fitmean << " " << fitrms << std::endl;
  //histname << lowername << " " <<  norminalmean << " " << norminalrms << " " << norminalmean << " " << norminalrms << std::endl;

        std::ostringstream textfilename;
	textfilename.str("");
	textfilename << "EleFakeRate-" << processName << "-";
	if(useDY)textfilename << "DY-";
	else textfilename << "Bw-";
	if(useKer)textfilename << "ker-";
	else textfilename << "expo-";
	if(inputbintype == 0)textfilename << "pt";
	else if(inputbintype == 1)textfilename << "eta";
	else if(inputbintype == 2)textfilename << "vtx";
	else if(inputbintype == 3) textfilename << "ht";

	if(inputfittype == 0) textfilename << "-den";
	else if(inputfittype == 1) textfilename << "-num";
	textfilename << "-" << lowercut << "-"<< uppercut<< "-" << fitrangelow << "-" << fitrangehigh <<  ".txt";
 
 // if(SaveOutput){
		ofstream myfile;
		myfile.open(textfilename.str().c_str(), std::ios_base::ate | std::ios_base::out);
		myfile << histname.str();
		myfile.close();
//  }
//  else std::cout << histname.str();
  cout<<"chi2 per ndf"<<chi2str<<endl;

}

int main(int argc, char *argv[])
{
using namespace RooFit ;
    if(argc < 8)
      std::cout << "You have to provide seven arguments!!\n";
   printf("number of arguments: %d\n", argc);
   FitKer(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atof(argv[5]),atoi(argv[6]),atoi(argv[7]));
   return 0;
}
