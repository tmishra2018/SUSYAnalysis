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
#include "../include/RooCMSShape.h"
#include "../include/RooDCBShape.h"
#include "../include/RooUserPoly.h"
#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_jet.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_tools.h"
#include "../include/analysis_mcData.h"
//#include "../include/RooCBExGaussShape.h"

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

void FitR9(int inputfittype, float lowerptcut, float upperptcut, float loweretacut, float upperetacut){//main  

  gSystem->Load("../lib/libRooFitClasses.so");
  
	bool usePoly(false);
	bool useKer(true);
	bool useDY(false);
	bool SaveOutput(true);
	float fitrangehigh = 120;
	float fitrangelow  = 60;

	char lowerptname[5];
	sprintf(lowerptname, "%d", (int)lowerptcut);	
  char upperptname[5];
	if(upperptcut < 1000){
		sprintf(upperptname, "%d", (int)upperptcut);
	}
	else sprintf(upperptname, "Inf");	

	char loweretaname[7];
	sprintf(loweretaname, "%d.%03d", loweretacut > 0? (int)loweretacut: -1*(int)loweretacut, abs((int)(((loweretacut-(int)loweretacut))*1000)));	
  char upperetaname[7];
	if(upperetacut < 1000){
		sprintf(upperetaname, "%d.%03d",upperetacut >0? (int)upperetacut: -1*(int)upperetacut, abs((int)(((upperetacut-(int)upperetacut))*1000)));	
	}
	else sprintf(upperetaname, "Inf");	

  std::ostringstream histname;      

//************** Process Z->ee Tree ********************************************************//   
  TChain *etree = new TChain("eeTree");																									
  etree->Add("/uscms_data/d3/mengleis/CMSSW_Work/CMSSW_8_0_26_patch1/src/EgammaAnalysis/TnPTreeProducer/plot_R9_data.root");

  float invmass=0; 
  float probePt=0; 
  float probeEta=0;
	float probeR9=0;
  etree->SetBranchAddress("invmass",   &invmass); 
  etree->SetBranchAddress("probeEt",   &probePt);
  etree->SetBranchAddress("probeEta",  &probeEta);
	etree->SetBranchAddress("probeR9",   &probeR9);

  TH1F*  p_invmass = new TH1F("p_invmass", "p_invmass",(int)(fitrangehigh-fitrangelow),fitrangelow,fitrangehigh);
	TTree* newetree  = new TTree("newetree","newetree"); 
	double invmass_newetree(0);
	newetree->Branch("invmass", &invmass_newetree); 

  for(unsigned iEvt(0); iEvt < etree->GetEntries(); iEvt++){
    etree->GetEntry(iEvt);
		bool vetovalue = ((fabs(probeEta) < 1.442 && probeR9 > 0.5) || (fabs(probeEta) > 1.56 && probeR9 > 0.8) );
		//bool vetovalue = probeMatchTrailing;
	
		if(probePt >= lowerptcut && probePt < upperptcut && probeEta >= loweretacut	&& probeEta < upperetacut){
			if(inputfittype == 0){p_invmass->Fill(invmass); invmass_newetree = invmass; newetree->Fill(); } 
			else if(inputfittype == 1 && vetovalue== true){p_invmass->Fill(invmass); invmass_newetree = invmass; newetree->Fill(); }
    }
  }

//************** Process mg Bkg Tree *******************************************//
  TTree *newbgtree = new TTree("newbgtree","newbgtree");
  double invmass_fordataset(0);
  newbgtree->Branch("invmass", &invmass_fordataset);

  TChain *bgtree = new TChain("BGTree");
  bgtree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_bgtemplate_FullEcal.root");
  float invmass_bg=0; 
  float probePt_bg=0; 
  float probeEta_bg=0; 
  bgtree->SetBranchAddress("invmass",   &invmass_bg); 
  bgtree->SetBranchAddress("probePt",   &probePt_bg);
  bgtree->SetBranchAddress("probeEta",  &probeEta_bg);

  TH1F* bg_pt;
  bg_pt = new TH1F("bg", "", 60, 60,120);
  for(unsigned iEvt(0); iEvt < bgtree->GetEntries(); iEvt++){
		bgtree->GetEntry(iEvt);
		if(probePt_bg >= lowerptcut && probePt_bg < upperptcut && probeEta_bg >= loweretacut	&& probeEta_bg < upperetacut){
		  bg_pt->Fill(invmass_bg); 
			invmass_fordataset=invmass_bg; 
			newbgtree->Fill();
		}
  }
	if(newbgtree->GetEntries() <= 0){
		useKer = false;
		usePoly = true;
	}


//************** Process DYJet Tree *******************************************//
  TTree *newDYtree = new TTree("newDYtree","newDYtree");
  double invmass_DYsignal(0);
  newDYtree->Branch("invmass", &invmass_DYsignal);
	TH1F  *h_DYinvmass = new TH1F("h_DYinvmass","h_DYinvmass",(int)(fitrangehigh-fitrangelow),fitrangelow,fitrangehigh);

  TChain *DYtree = new TChain("FakeRateTree");
  DYtree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_R9_DY.root");

  float DY_invmass=0; 
  float DY_tagPt=0; 
  float DY_tagEta=0; 
  float DY_tagPhi=0; 
  float DY_probePt=0; 
  float DY_probeEta=0;
  float DY_probePhi=0;
	float DY_probeR9 = 0;
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
  DYtree->SetBranchAddress("probeR9",   &DY_probeR9);
  DYtree->SetBranchAddress("mcPID",			&mcPID);
  DYtree->SetBranchAddress("mcEta",			&mcEta);
  DYtree->SetBranchAddress("mcPhi",			&mcPhi);
  DYtree->SetBranchAddress("mcPt",			&mcPt);
  DYtree->SetBranchAddress("mcMomPID",	&mcMomPID);
  DYtree->SetBranchAddress("mcGMomPID",	&mcGMomPID);

  for(unsigned iEvt(0); iEvt < DYtree->GetEntries(); iEvt++){
    DYtree->GetEntry(iEvt);

		if( (fabs(DY_probeEta) < 1.442 && DY_probeR9 < 0.5 ) || (fabs(DY_probeEta) > 1.56 && DY_probeR9 < 0.8 ))continue;
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
		if(DY_probePt >= lowerptcut && DY_probePt < upperptcut && fabs(DY_probeEta) >= loweretacut	&& fabs(DY_probeEta) < upperetacut){
      if(isZee){invmass_DYsignal = DY_invmass; h_DYinvmass->Fill(DY_invmass); }
			if(isZee)newDYtree->Fill();		
    }
  }
//************* Construct RooFit Models *********************************************************//
  RooRealVar mass_axis("invmass","invmass",fitrangelow, fitrangehigh);
	RooDataSet *ZeeDataSet;
	if(newetree->GetEntries() < 50000)ZeeDataSet = new RooDataSet("ZeeDataSet","ZeeDataSet",RooArgSet(mass_axis),RooFit::Import(*newetree));
	RooDataSet *DYDataSet;
	RooDataHist *datahist_DY;
	if(useDY && h_DYinvmass->GetEntries() < 50)DYDataSet = new RooDataSet("DYDataSet","DYDataSet",RooArgSet(mass_axis), RooFit::Import(*newDYtree));
	else if(useDY && h_DYinvmass->GetEntries() >= 50)datahist_DY = new RooDataHist("DYDataSet","DYDataSet", mass_axis, h_DYinvmass);
	RooKeysPdf *DYKer;
	RooHistPdf *DYpdf;
	if(useDY && h_DYinvmass->GetEntries() < 50)DYKer = new RooKeysPdf("DYKer","DYKer",mass_axis, *DYDataSet,RooKeysPdf::MirrorBoth,2);
 	else if(useDY && h_DYinvmass->GetEntries() >= 50){
		DYpdf = new RooHistPdf("DYpdf","DYpdf", mass_axis, *datahist_DY); 
		DYpdf->setInterpolationOrder(1); 
	}
  RooDataSet BkgDataSet("BkgDataSet","BkgDataSet",RooArgSet(mass_axis),RooFit::Import(*newbgtree));
  RooKeysPdf BkgKer("BkgKer","BkgKer",mass_axis,BkgDataSet,RooKeysPdf::MirrorBoth, 2); 

  TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
  RooDataHist datahist_data("both", "", mass_axis, p_invmass);

  c_fitMass->cd();
	double slopeupper = (p_invmass->GetBinContent(31) - p_invmass->GetBinContent(1))/31.0;
	double slopelower = (p_invmass->GetBinContent((int)(fitrangehigh-fitrangelow)) - p_invmass->GetBinContent(30))/30.0;
  double inislope = (p_invmass->GetBinContent(50) - p_invmass->GetBinContent(10))/40.0;
  RooRealVar slope("slope","slope",inislope,slopelower,slopeupper);
  RooUserPoly *pol = new RooUserPoly("px","px",mass_axis, slope);
	double inilambda = (log(p_invmass->GetBinContent(2)) - log(p_invmass->GetBinContent(48)))/(p_invmass->GetBinCenter(2) - p_invmass->GetBinCenter(48)); 
	RooRealVar lambda("lambda", "slope", inilambda, -100, 100.);
	RooExponential *expo = new RooExponential("expo", "exponential PDF", mass_axis, lambda);

  RooRealVar m0( "m0", "m0", 91.188, 80,100);
  RooRealVar width( "width", "width", 2.495, 0, 15);
  //RooRealVar m0( "m0", "m0", 91.188);
  //RooRealVar width( "width", "width", 2.495);
  RooRealVar mean("mean", "" ,0.,-1,1);
  //RooRealVar sigma("sigma", "",2.4 , 0.0, 15.0);
  RooRealVar sigma("sigma", "", 1.0, 0, 2);
  RooRealVar alpha("alpha", "", 1.0, 0.0, 20.0);
  RooRealVar n("n","", 1.0, 0.0, 20.0);
  RooRealVar alpha2("2ndalpha","", 1.0, 0.0, 20.0);
  RooRealVar n2("2ndn", "", 1.0, 0.0, 20.0);
  RooBreitWigner bw("bw", "", mass_axis, m0, width);
  RooDCBShape *cb;
  cb = new RooDCBShape("cb","cb", mass_axis, mean, sigma, alpha, n, alpha2, n2);
  RooGaussian gauss("gs", "gs", mass_axis, mean, sigma);
  RooFFTConvPdf *signalRes;
	if(!useDY)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, bw, *cb);
  else if(useDY && h_DYinvmass->GetEntries() < 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYKer, gauss);
	else if(useDY && h_DYinvmass->GetEntries() >= 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYpdf, gauss);
  double iniSig = 0.8*p_invmass->Integral(1,100);
  double iniBkg = 60*(p_invmass->GetBinContent(10)+p_invmass->GetBinContent(50));
	if(iniBkg > p_invmass->GetEntries()/2)iniBkg = p_invmass->GetEntries()/2;
  //RooRealVar nSig("nSig", "", 0.7*iniSig, 0, p_invmass->GetEntries()*1.2);
  RooRealVar nSig("nSig", "", iniSig, 0, p_invmass->GetEntries()*1.2);
  RooRealVar nBkg("nBkg", "", iniBkg, 0, p_invmass->GetEntries());
  RooAddPdf *model;
  //if(useCMSShape)model = new RooAddPdf("model", "", RooArgList(CMSShape, signalRes),RooArgList(nBkg, nSig));
  if(usePoly)model = new RooAddPdf("model", "", RooArgList(*expo, *signalRes),RooArgList(nBkg, nSig));   
  else if(useKer)model = new RooAddPdf("model", "", RooArgList(BkgKer, *signalRes),RooArgList(nBkg, nSig));   

  histname.str("");
  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
  histname << lowerptname << "_" << upperptname << loweretaname << "_" << upperetaname;
  RooPlot* mass_Frame = mass_axis.frame(RooFit::Title(histname.str().c_str()),RooFit::Bins(100));
  if(newetree->GetEntries() < 50000)model->fitTo(*ZeeDataSet);
  else model->fitTo(datahist_data);
  //if(useCMSShape)model->plotOn(mass_Frame, RooFit::Components(CMSShape),
	//			 RooFit::LineStyle(kDashed),
	//			 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  if(usePoly)model->plotOn(mass_Frame, RooFit::Components(*expo),
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame, RooFit::Components(BkgKer),
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  if(usePoly)model->plotOn(mass_Frame,
				 RooFit::Components(RooArgSet(*expo, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame,
				 RooFit::Components(RooArgSet(BkgKer, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  if(newetree->GetEntries() < 50000)ZeeDataSet->plotOn(mass_Frame);
  else datahist_data.plotOn(mass_Frame);
  mass_Frame->Draw();

  double chi2 = mass_Frame->chiSquare();
  TLatex* latex = new TLatex();
  char chi2str[50];
  sprintf (chi2str, "chi2/ndof = %f", chi2);
  latex->DrawLatex(90,0.5*p_invmass->GetMaximum(),chi2str);
  histname.str("");
	if(useDY)histname << "DY_";
	else histname << "Bw_";
  if(useKer)histname << "ker_";
  else if(usePoly)histname << "expo_";
  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
  histname << lowerptname << "_" << upperptname << "_" <<  loweretaname << "_" << upperetaname  <<"R9.pdf"; 
  c_fitMass->SaveAs(histname.str().c_str());

  mass_axis.setRange("signal",60,120);
  RooAbsReal* igx_sig;
	igx_sig = signalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	double norminalmean = igx_sig->getVal()*(nSig.getVal());
	double norminalrms  = igx_sig->getVal()*(nSig.getError());
	//if(p_invmass->Integral(1,40) < 500)norminalmean = p_invmass->Integral(10,31);

  histname.str("");
  if(inputfittype == 0) histname << "den ";
  else if(inputfittype == 1) histname << "num ";
  histname << lowerptname << " " << loweretaname << " " <<  norminalmean << " " << norminalrms << std::endl;
	
	std::ostringstream textfilename;
	textfilename.str("");
	textfilename << "eeR9";
	textfilename << fitrangelow << "-" << fitrangehigh <<  "-FullIntegral.txt";
 
  if(SaveOutput){
		ofstream myfile;
		myfile.open(textfilename.str().c_str(), std::ios_base::app | std::ios_base::out);
		myfile << histname.str();
		myfile.close();
  }
  else std::cout << histname.str();
}
