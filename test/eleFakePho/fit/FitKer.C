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
//#define NTOY 10000
using namespace RooFit ;
int RunYear = 2016;
bool ISpreVFP = false;

bool doEB = true;
//const char*processName ="Data";
const char*processName ="DY";

enum BinType{
  byPt = 0,
  byEta = 1,
  byVtx = 2,
  nonType = 3
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
	//bool SaveOutput;
        ifstream configfile;
        if(fitFunc==1)
		configfile.open ("configFit_Bw_ker.txt");
        if(fitFunc==2)
		configfile.open ("configFit_DY_ker.txt");
        if(fitFunc==3)
		configfile.open ("configFit_Bw_expo.txt");
	std::string variabletype;
	bool variablevalue;
	if(configfile.is_open()){
  	for(int i(0); i<=4; i++){ 
			configfile >> variabletype >> variablevalue; 
			if(variabletype.find("Expo")!=std::string::npos)useExpo = variablevalue;
			else if(variabletype.find("Ker")!=std::string::npos)useKer = variablevalue;
			else if(variabletype.find("DY")!=std::string::npos)useDY = variablevalue;
	  }
	}
	configfile.close();
 
  BinType bintype = (BinType)inputbintype;

	char lowername[5];
	if(bintype == byEta)sprintf(lowername, "%d.%02d", (int)lowercut, (int)((lowercut-(int)lowercut)*100));	
	else sprintf(lowername, "%d", (int)lowercut);	
  char uppername[5];
	if(uppercut < 1000){
		if(bintype == byEta)sprintf(uppername, "%d.%02d", (int)uppercut, (int)((uppercut-(int)uppercut)*100));	
		else sprintf(uppername, "%d", (int)uppercut);
	}
	else sprintf(uppername, "Inf");	

  std::ostringstream histname;      

	TH1F *Zee_PU = new TH1F("Zee_PU","",100,0,100);
	TH1F *DY_PU  = new TH1F("DY_PU", "",100,0,100);

//************** Process Z->ee Tree ********************************************************//   
  TChain *etree = new TChain("FakeRateTree");
  etree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/plot_elefakepho_%sTnP_dR05_%d%s.root",processName,RunYear,whichVFP.c_str()));
  float invmass=0; 
  float probePt=0; 
  float probeEta=0;
  bool  vetovalue=0;
	bool  FSRveto=0;
  int   nVertex=0;
  etree->SetBranchAddress("invmass",   &invmass); 
  etree->SetBranchAddress("probePt",   &probePt);
  etree->SetBranchAddress("probeEta",  &probeEta);
  etree->SetBranchAddress("vetovalue", &vetovalue);
	etree->SetBranchAddress("FSRveto",   &FSRveto);
  etree->SetBranchAddress("nVertex",   &nVertex);

  TH1F*  p_invmass = new TH1F("p_invmass", "p_invmass",(int)(fitrangehigh-fitrangelow),fitrangelow,fitrangehigh);
	TTree* newetree  = new TTree("newetree","newetree");                     // newetree is
	double invmass_newetree(0);
	newetree->Branch("invmass", &invmass_newetree); 

  TH1F *his = new TH1F("his","",100,0,100);
  etree->Draw("nVertex>>his");
  cout<< "range to consider ....   "<< his->GetMean()-his->GetRMS() << "-"<< his->GetMean()+his->GetRMS()<<endl;
  
  TH1F *hisPt = new TH1F("hisPt","",220,0,1100);
  etree->Draw("probePt>>hisPt");
  cout<< "range to consider ....   "<< hisPt->GetMean()-hisPt->GetRMS() << "-"<< hisPt->GetMean()+hisPt->GetRMS()<<endl;
  
  for(unsigned iEvt(0); iEvt < etree->GetEntries(); iEvt++){
    etree->GetEntry(iEvt);
		if(doEB && fabs(probeEta) > 1.4442)continue;
		else if(!doEB && fabs(probeEta) < 1.56)continue;	
	
    float keyvariable(0);
    if(bintype == byPt)keyvariable = probePt;
    else if(bintype == byEta)keyvariable = fabs(probeEta);
    else if(bintype == byVtx)keyvariable = nVertex;

    if(bintype == byPt and (nVertex < his->GetMean()-his->GetRMS() or nVertex > his->GetMean()+his->GetRMS() )) continue; // consider pT when 10<nVtx<25
    if(bintype == byVtx and (probePt < hisPt->GetMean()-hisPt->GetRMS() or probePt > hisPt->GetMean()+hisPt->GetRMS() )) continue; // consider nVtx when 
    
    if(keyvariable >= lowercut && keyvariable < uppercut){ 
			if(inputfittype == 0 && vetovalue== false){p_invmass->Fill(invmass); invmass_newetree = invmass; newetree->Fill(); } 
			else if(inputfittype == 1 && vetovalue== true && FSRveto == true){p_invmass->Fill(invmass); invmass_newetree = invmass; newetree->Fill(); }
			Zee_PU->Fill(nVertex);
    }
  }

//************** Process mg Bkg Tree *******************************************//
  TTree *newbgtree = new TTree("newbgtree","newbgtree");                        // newbgtree is for μ + probe template(Background Model)
  double invmass_fordataset(0);
  newbgtree->Branch("invmass", &invmass_fordataset);

  TChain *bgtree = new TChain("BGTree");
  bgtree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/plot_bgtemplate_FullEcal_%d%s.root",RunYear,whichVFP.c_str()));
  
  float invmass_bg=0; 
  float probePt_bg=0; 
  float probeEta_bg=0; 
  bool  vetovalue_bg=0;
  bool  FSRveto_bg=0;	
  int   nVertex_bg=0;
  bgtree->SetBranchAddress("invmass",   &invmass_bg); 
  bgtree->SetBranchAddress("probePt",   &probePt_bg);
  bgtree->SetBranchAddress("probeEta",   &probeEta_bg);
  bgtree->SetBranchAddress("vetovalue", &vetovalue_bg);
  bgtree->SetBranchAddress("FSRveto",   &FSRveto_bg);
  bgtree->SetBranchAddress("nVertex",   &nVertex_bg);

  TH1F* bg_pt;
  if(bintype == byPt)bg_pt = new TH1F("bg", "", 95, 55,150);
  else bg_pt = new TH1F("bg", "bg",150,50,200);
  for(unsigned iEvt(0); iEvt < bgtree->GetEntries(); iEvt++){
		bgtree->GetEntry(iEvt);
    if(doEB && fabs(probeEta_bg) > 1.4442)continue;
    else if(!doEB && fabs(probeEta_bg) < 1.56)continue;

    float keyvariable(0);
    if(bintype == byPt)keyvariable = probePt_bg;
    else if(bintype == byEta)keyvariable = fabs(probeEta_bg);
    else if(bintype == byVtx)keyvariable = nVertex_bg;

    if(bintype == byPt and (nVertex_bg < his->GetMean()-his->GetRMS()  or nVertex_bg > his->GetMean()+his->GetRMS() )) continue; // consider pT when 10<nVtx<25
    if(bintype == byVtx and (probePt_bg < hisPt->GetMean()-hisPt->GetRMS() or probePt_bg > hisPt->GetMean()+hisPt->GetRMS() )) continue; // consider nVtx when 

    if(keyvariable >= lowercut && keyvariable < uppercut){ 
			if(inputfittype == 0 && vetovalue_bg== false){ bg_pt->Fill(invmass_bg); invmass_fordataset=invmass_bg; newbgtree->Fill();}
			else if(inputfittype == 1 && vetovalue_bg== true && FSRveto_bg==true){bg_pt->Fill(invmass_bg); invmass_fordataset=invmass_bg; newbgtree->Fill();}
    }
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

    			if(bintype == byPt and (DY_nVertex < his->GetMean()-his->GetRMS()  or DY_nVertex > his->GetMean()+his->GetRMS())) continue; // consider pT when 10<nVtx<25
    			if(bintype == byVtx and (DY_probePt < hisPt->GetMean()-hisPt->GetRMS() or DY_probePt > hisPt->GetMean()+hisPt->GetRMS() )) continue; // consider nVtx when 

	//		int tagbin_x = h_elescale->GetXaxis()->FindBin(DY_tagPt);
	//		int tagbin_y = h_elescale->GetYaxis()->FindBin(fabs(DY_tagEta));
	//		int probebin_x=h_phoscale->GetXaxis()->FindBin(DY_probeEta);
	//		int probebin_y=h_phoscale->GetYaxis()->FindBin(fabs(DY_probePt));
	//		double eventweight = h_elescale->GetBinContent(h_elescale->GetBin(tagbin_x, tagbin_y))*h_eleiso->GetBinContent(h_eleiso->GetBin(tagbin_x, tagbin_y))*h_phoscale->GetBinContent(h_phoscale->GetBin(probebin_x, probebin_y));
			if(keyvariable >= lowercut && keyvariable < uppercut){
				double mcPUweight;
				if(RunYear==2016 and ISpreVFP==true) mcPUweight = getPUESF16preVFP(DY_nVertex);
				else if(RunYear==2016 and ISpreVFP==false) mcPUweight = getPUESF16(DY_nVertex);
				else if(RunYear==2017) mcPUweight = getPUESF17(DY_nVertex);
				else if(RunYear==2018) mcPUweight = getPUESF18(DY_nVertex);
				DY_PU->Fill(DY_nVertex, mcPUweight);
				if(isZee){invmass_DYsignal = DY_invmass; h_DYinvmass->Fill(DY_invmass, mcPUweight); }
			//  if(inputfittype == 0 && DY_vetovalue== false){invmass_DYsignal = DY_invmass; h_DYinvmass->Fill(DY_invmass, eventweight); }
			//	if(isZee){invmass_DYsignal = DY_invmass; h_DYinvmass->Fill(DY_invmass); }
				if(isZee)newDYtree->Fill();		
			}
		}
		h_DYinvmass->Sumw2();
	} // Run only useDY

  //************* Construct RooFit Models *********************************************************//
  RooRealVar mass_axis("invmass","M_{tag-probe}",fitrangelow, fitrangehigh);
  //	RooDataSet *TnPDataSet;
  //	if(newetree->GetEntries() < 50000) TnPDataSet = new RooDataSet("TnPDataSet","TnPDataSet",RooArgSet(mass_axis),RooFit::Import(*newetree));
  RooDataSet BkgDataSet("BkgDataSet","BkgDataSet",RooArgSet(mass_axis),RooFit::Import(*newbgtree));
  RooKeysPdf BkgKer("BkgKer","BkgKer",mass_axis,BkgDataSet,RooKeysPdf::MirrorBoth, 2);
	RooDataSet *DYDataSet;
  // RooKeysPdf : Kernel estimation, A p.d.f. that represent the shape of an external unbinned dataset as a superposition of Gaussians with equal surface
	RooDataHist *datahist_DY;
	RooKeysPdf *DYKer;
	RooHistPdf *DYpdf;

        // for less events use RooDataSet (unbinned) and for more events use RooDataHist (binned)
	if(useDY && h_DYinvmass->GetEntries() < 50)DYDataSet = new RooDataSet("DYDataSet","DYDataSet",RooArgSet(mass_axis), RooFit::Import(*newDYtree));
	else if(useDY && h_DYinvmass->GetEntries() >= 50)datahist_DY = new RooDataHist("DYDataSet","DYDataSet", mass_axis, h_DYinvmass);

        // for less events use RooDataSet->RooKeysPdf and for more events use RooDataHist->RooHistPdf
	if(useDY && h_DYinvmass->GetEntries() < 50)DYKer = new RooKeysPdf("DYKer","DYKer",mass_axis, *DYDataSet,RooKeysPdf::MirrorBoth,2);
 	else if(useDY && h_DYinvmass->GetEntries() >= 50){
		DYpdf = new RooHistPdf("DYpdf","DYpdf", mass_axis, *datahist_DY); 
		DYpdf->setInterpolationOrder(1); 
	}
/*********************** Plot Kernel Background ********************************/ 
//  TCanvas *c_ker = new TCanvas("ker","",600,600);														 
//  c_ker->cd();
//  RooPlot* kerFrame = mass_axis.frame(RooFit::Title("kernel"));
//  BkgDataSet.plotOn(kerFrame);
//  BkgKer.plotOn(kerFrame);
//  kerFrame->Draw();
//  histname.str("");
//  histname << "Ker_";
//  if(inputbintype == 0)histname << "pt_";
//  else if(inputbintype == 1)histname << "eta_";
//  else if(inputbintype == 2)histname << "vtx_";
//  if(inputfittype == 0) histname << "den_";
//  else if(inputfittype == 1) histname << "num_";
//  histname << lowername << "_" << uppername  << ".pdf"; 
//  if(useKer)c_ker->SaveAs(histname.str().c_str());
/******************************************************************************/


  TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
	c_fitMass->SetLeftMargin(0.15);
	c_fitMass->SetBottomMargin(0.15);
  RooDataHist datahist_data("both", "", mass_axis, p_invmass);

  c_fitMass->cd();
  //c_fitMass->SetLogy();



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


  RooFFTConvPdf *signalRes;                                        // can convolve any two RooAbsPdfs with fast fourier transform
	if(!useDY)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, bw, *cb);                                                  // CONVOLUTION of bw & cb
  else if(useDY && h_DYinvmass->GetEntries() < 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYKer, gauss);            // CONVOLUTION of DY and gauss
	else if(useDY && h_DYinvmass->GetEntries() >= 50)signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, *DYpdf, gauss);     // when events are less use Kernel estimator



  //double iniSig = 0.8*p_invmass->Integral(1,100);
  //double iniBkg = 60*(p_invmass->GetBinContent(10)+p_invmass->GetBinContent(50));

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

 //if(useCMSShape)model = new RooAddPdf("model", "", RooArgList(CMSShape, signalRes),RooArgList(nBkg, nSig));
  if(useExpo)model = new RooAddPdf("model", "", RooArgList(*expo, *signalRes),RooArgList(nBkg, nSig));                             // use expo background model
  else if(useKer)model = new RooAddPdf("model", "", RooArgList(BkgKer, *signalRes),RooArgList(nBkg, nSig));                        // use μ + prob background model

  histname.str("");
  if(inputfittype == 0) histname << "denominator ";
  else if(inputfittype == 1) histname << "numerator ";
	histname << lowername << " < ";
  if(inputbintype == 0)histname << "pt <";
  else if(inputbintype == 1)histname << "eta <";
  else if(inputbintype == 2)histname << "nvtx <";
  histname << uppername;
  RooPlot* mass_Frame = mass_axis.frame(RooFit::Title(histname.str().c_str()),RooFit::Bins(60));
  mass_Frame->SetStats(0);
  
  //if(newetree->GetEntries() < 50000){
  //	cout<<"fitting TnPDataSet"<<endl;
  //	model->fitTo(*TnPDataSet);
  //	}
  //else {
  	model->fitTo(datahist_data);
  //	cout<<"fitting datahist_data"<<endl;}
	
  //if(useCMSShape)model->plotOn(mass_Frame, RooFit::Components(CMSShape),
	//			 RooFit::LineStyle(kDashed),
	//			 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  if(useExpo)model->plotOn(mass_Frame, RooFit::Components(*expo),                                ///  plot background
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame, RooFit::Components(BkgKer),
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  if(useExpo)model->plotOn(mass_Frame,                                                           ///  plot signal + background
				 RooFit::Components(RooArgSet(*expo, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  else if(useKer)model->plotOn(mass_Frame,
				 RooFit::Components(RooArgSet(BkgKer, *signalRes)),
				 RooFit::LineStyle(kSolid),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  //if(newetree->GetEntries() < 50000)
  //TnPDataSet->plotOn(mass_Frame);                               ///  plot data
  //else datahist_data.plotOn(mass_Frame);
  datahist_data.plotOn(mass_Frame);

  //mass_Frame->SetMaximum(2E6);
  //mass_Frame->SetTitleSize(0.01);
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
  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
	histname << fitrangelow << "-" << fitrangehigh << "_"; 
  histname << lowername << "_" << uppername <<".png";
  c_fitMass->SaveAs(histname.str().c_str());
  // $ROOTSYS/tutorial/roofit/rf110_normintegration.C
  mass_axis.setRange("signal",80,101);
  RooAbsReal* igx_sig;
        // fraction of of p.d.f. gx_Norm[x] which is in the range named "signal"
	igx_sig = signalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	double norminalmean = igx_sig->getVal()*(nSig.getVal());
	double norminalrms  = igx_sig->getVal()*(nSig.getError());

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  RooHistPdf histpdf_data("histpdf_data","",mass_axis,datahist_data); 
  RooMCStudy* mcstudy = new RooMCStudy(histpdf_data,mass_axis,RooFit::FitModel(*model),RooFit::Binned(kTRUE),RooFit::Silence(),
                        RooFit::Extended(),RooFit::FitOptions(RooFit::Save(kTRUE),RooFit::PrintEvalErrors(0)));
//  RooChi2MCSModule chi2mod;
//  mcstudy->addModule(chi2mod);
  mcstudy->generateAndFit(NTOY,p_invmass->Integral(1, (int)(fitrangehigh-fitrangelow)));
  mcstudy->fitParDataSet().Print("v");
//  TH1* p_chi2 = mcstudy->fitParDataSet().createHistogram("chi2"); 
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
  if(inputfittype == 0) histname << "den_";
  else if(inputfittype == 1) histname << "num_";
  histname << lowername << "_" << uppername << ".png";
  toymc->SaveAs(histname.str().c_str());

  histname.str("");
  if(inputbintype == 0)histname << "pt ";
  else if(inputbintype == 1)histname << "eta ";
  else if(inputbintype == 2)histname << "vtx ";
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
