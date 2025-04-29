#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "../../../include/tdrstyle.C"
#include "TROOT.h"

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
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "RooMultiVarGaussian.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
// #include "../../../include/analysis_fakes.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
#define NTOY 100

#define MINPT 30
#define MAXPT 200
#define MAXVTX 46


int MINVTX = 0;
bool doEB = false;
//const char*processName ="Data"; bool doDrellYan = false;

const char*processName ="DY"; bool doDrellYan = true;

Double_t fakerate_ptDependence(Double_t *x, Double_t *par)
{
	double slope = par[0];
	double constant = par[1]; 
	double index = par[2];
	double coeff = 1.0; 
	double pt = TMath::Max(x[0],0.000001);
	double arg = 0;
	arg = slope*pt + constant; 
	double fitval = pow(arg, index)*coeff; 
	return fitval;
}


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

void DY_fakerates(int RunYear, bool preVFP){//main 
	gROOT->SetBatch(1);
	std::string whichVFP;
        if(RunYear==2016 and preVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
	if (!doDrellYan){
		if(RunYear==2017) MINVTX = 18; // changing MINVTX here
		else MINVTX = 12; // changing MINVTX here
	}

	std::string isEB;
	if(doEB == true) isEB = "";
	else if(doEB == false) isEB = "_EE";

	
	ofstream resultfile;
	ofstream fakeRateFile;
	ofstream etafile;
	if (doEB == true) {
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPtVtx-EB.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/result_fitFakeFunc_Data.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-ByPtVtx-EB.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/result_eta_dependence_Data.txt", RunYear, whichVFP.c_str()));
    		}} 
	else {  
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPtVtx-EE.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY_EE.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/result_fitFakeFunc_Data_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-ByPtVtx-EE.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/result_eta_dependence_Data_EE.txt", RunYear, whichVFP.c_str()));}}


	gSystem->Load("../../../lib/libAnaClasses.so");
        gSystem->Load("../../../lib/libRooFitClasses.so");
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	
	std::ifstream Pt_file, Pt_DYfile, Pt_Polfile, Eta_file, Eta_DYfile, Eta_Polfile, Vtx_file, Vtx_DYfile, Vtx_Polfile;

	/****************************   MC      *********************************/
	
	if(doDrellYan==true) {
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
     		Vtx_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-expo-vtx-60-120.txt",RunYear,whichVFP.c_str())); 
		
	}

	/****************************   Data   *********************************/
	else {
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
		if(RunYear==2016 and preVFP == true){
		       	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
			Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));}
		else{
			Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        		Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-pt-60-120.txt",RunYear,whichVFP.c_str()));}
        	


		Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-eta-60-120.txt",RunYear,whichVFP.c_str()));

		if(RunYear==2017){
        		Vtx_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-vtx-60-120_above18.txt",RunYear,whichVFP.c_str()));
        		Vtx_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-vtx-60-120_above18.txt",RunYear,whichVFP.c_str()));
        		Vtx_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-vtx-60-120_above18.txt",RunYear,whichVFP.c_str())); }
        	else {
        		Vtx_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-vtx-60-120_above12.txt",RunYear,whichVFP.c_str()));
        		Vtx_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-vtx-60-120_above12.txt",RunYear,whichVFP.c_str()));
        		Vtx_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-vtx-60-120_above12.txt",RunYear,whichVFP.c_str())); }
        }
	std::string line;
	unsigned nPtBins(0);
	unsigned nEtaBins(0);
	unsigned nVtxBins(0); 
	while(std::getline(Pt_file , line))
		nPtBins += 1;
	std::cout << "Pt_file total line : " << nPtBins << std::endl;
	Pt_file.clear(); Pt_file.seekg(0, ios::beg);
	while(std::getline(Eta_file , line))
		nEtaBins += 1;
	std::cout << "Eta_file total line : " << nEtaBins << std::endl;
	Eta_file.clear(); Eta_file.seekg(0, ios::beg);
	while(std::getline(Vtx_file , line))
		nVtxBins += 1;
	std::cout << "Vtx_file total line : " << nVtxBins << std::endl;
	Vtx_file.clear(); Vtx_file.seekg(0, ios::beg);
	
	nPtBins /= 2;
	nEtaBins /= 2;
	nVtxBins /= 2;

	double PtBins[nPtBins];
	double EtaBins[nEtaBins];
	double VtxBins[nVtxBins];
	
	double pt_den[nPtBins];
	double pt_denerror[nPtBins];
	double pt_num[nPtBins];
	double pt_numerror[nPtBins];

	double eta_den[nEtaBins];
	double eta_denerror[nEtaBins];
	double eta_num[nEtaBins];
	double eta_numerror[nEtaBins];

	double vtx_den[nVtxBins];
	double vtx_denerror[nVtxBins];
	double vtx_num[nVtxBins];
	double vtx_numerror[nVtxBins];
	std::string bintype;
	std::string numtype;
	float lowcut;
	float signal1, error1;
	float fitmean, fitrms;
	float DYmean,  DYrms;
	float Polmean, Polrms;
	TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins);
	TGraphErrors *fr_bothcount_eta= new TGraphErrors(nEtaBins);
	TGraphErrors *fr_bothcount_vtx= new TGraphErrors(nVtxBins);
	TGraphErrors *fr_pt_ratio = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_pt_ratioError = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_vtx_ratio = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_vtx_ratioError = new TGraphErrors(nVtxBins-1);

	int graphPtBins = (int)((MAXPT - MINPT)/0.5);
	int graphVtxBins = (int)((MAXVTX - MINVTX)/0.5);
	TGraphErrors *fr_pt_sigmaband = new TGraphErrors(graphPtBins);
	TGraphErrors *fr_vtx_sigmaband = new TGraphErrors(graphVtxBins);
	TGraphErrors *fr_eta_sigmaband = new TGraphErrors(nEtaBins-1);

//************************* Fill and fit the pt dependence. Calculate ratio *********************// 
	if(Pt_file.is_open() && Pt_DYfile.is_open() && Pt_Polfile.is_open()){
		for(unsigned i(0); i<nPtBins; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			PtBins[i] = lowcut;
			pt_den[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			pt_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}    
		for(unsigned i(0); i<nPtBins; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			pt_num[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			pt_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}    
	}
	for(unsigned i(0); i<nPtBins; i++){
		double fakerate = pt_num[i]/pt_den[i];
		double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
		double xvalue = (i+1 < nPtBins)? (PtBins[i]+PtBins[i+1])/2.0 : (MAXPT + PtBins[i])/2.0;
		double xerror = (i+1 < nPtBins)? (PtBins[i+1]-PtBins[i])/2.0 : (MAXPT - PtBins[i])/2.0;
		fr_bothcount_pt->SetPoint(i, xvalue, fakerate );
		fr_bothcount_pt->SetPointError(i, xerror, error);
	}


	TF1 *f1 = new TF1("f1", fakerate_ptDependence,30,1000,3);
	TF1 *fit_fakerate_pt;
	TFitResultPtr result_fitpt;
	TVirtualFitter::SetMaxIterations(1000000);
// 	for(int i(1); i<50; i++){
// 		for(int j(1); j <2500; j++){
// 			std::cout << i << " " << j << std::endl;
//			f1->SetParameter(0, i);
//			f1->SetParameter(1, 206);
//			f1->SetParameter(2, -4 + j*4.0/2500);
//			f1->SetParNames("slope","constant","index");
//			result_fitpt = fr_bothcount_pt->Fit("f1","R S");
//			if(  result_fitpt->ParError(0)/result_fitpt->Parameter(0) < 1 && result_fitpt->ParError(1)/result_fitpt->Parameter(1) <1 && result_fitpt->ParError(2)/result_fitpt->Parameter(2) < 1)std::cout << "good fit status " << fr_bothcount_pt->Fit("f1","R S") << std::endl;
// 		}
// 	}
	//*************** EB ******************//
	if(RunYear==2016){

	f1->SetParameter(0,1.0e+04);
        f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
        f1->SetParameter(1, -10000);
        f1->SetParLimits(1, -15000.0, -5000);
        f1->SetParameter(2,-0.4);
        f1->SetParLimits(2, -0.4, -0.1);

/*		f1->SetParameter(0,1.0e+04);
        	f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
        	f1->SetParameter(1, -38000);
        	f1->SetParLimits(1, -38100.0, -10000);
        	f1->SetParameter(2,-0.4);
        	f1->SetParLimits(2, -0.4, -0.1);*/


	/*f1->SetParameter(0, 846.34);
	f1->SetParLimits(0, 800, 900);
	f1->SetParameter(1, -25000);
	f1->SetParLimits(1, -26000, -24000);
	f1->SetParameter(2, -0.4);
	f1->SetParLimits(2, -0.8, -0.1);*/

	}
	if(RunYear==2017){
		f1->SetParameter(0,1.0e+04);
		f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
                f1->SetParameter(1, -20000);
                f1->SetParLimits(1, -39070.0, -10000);
                f1->SetParameter(2,-0.4);
                f1->SetParLimits(2, -0.4, -0.1);
	}
	if(RunYear==2018){
		f1->SetParameter(0,1.0e+04);
		f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
                f1->SetParameter(1, -30000);
                f1->SetParLimits(1, -38685.0, 0);
                f1->SetParameter(2,-0.4);
                f1->SetParLimits(2, -0.40, -0.1);
	}
	
	/*if(RunYear==2018){
        	f1->SetParameter(0, 9900);
        	f1->SetParLimits(0, 0, 100000);
        	f1->SetParameter(1, 206);
        	f1->SetParameter(2, -4 + 955*4.0/2500);
        	f1->SetParLimits(2, -4, 0);}*/
	//*************** MC ******************//
	if(doDrellYan==true and RunYear == 2016){
		f1->SetParameter(0, 1);
		f1->SetParameter(1, 206);
		f1->SetParameter(2, -0.86);}
	else if(doDrellYan==true and RunYear == 2017){
        	f1->SetParameter(0, 29);
        	f1->SetParLimits(0, 0, 1000);
        	f1->SetParameter(1, 206);
        	f1->SetParameter(2, -4 + 955*4.0/2500);
        	f1->SetParLimits(2, -4, 0);}

	else if(doDrellYan==true and RunYear==2018){
                f1->SetParameter(0,1.0e+04);
                f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
                f1->SetParameter(1, -30000);
                f1->SetParLimits(1, -38685.0, 0);
                f1->SetParameter(2,-0.4);
                f1->SetParLimits(2, -0.40, -0.1);}
        



	f1->SetParNames("slope","constant","index");
        result_fitpt = fr_bothcount_pt->Fit("f1","R S");

 	fit_fakerate_pt = fr_bothcount_pt->GetFunction("f1");
	resultfile <<"pt chiSquare() : "<<fit_fakerate_pt->GetChisquare()/fit_fakerate_pt->GetNDF()<<endl;
 	for(unsigned i(0); i<nPtBins; i++){
 		double fakerate = pt_num[i]/pt_den[i];
 		double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
		double xvalue = (i+1 < nPtBins)? (PtBins[i]+PtBins[i+1])/2.0 : (MAXPT + PtBins[i])/2.0;
		double xerror = (i+1 < nPtBins)? (PtBins[i+1]-PtBins[i])/2.0 : (MAXPT - PtBins[i])/2.0;
 		fr_pt_ratio->SetPoint(i, xvalue, fit_fakerate_pt->Eval(xvalue)/fakerate);
 		fr_pt_ratio->SetPointError(i, xerror, error/fit_fakerate_pt->Eval(xvalue));
 	}
 
 //********************* Fill and fit the eta dependence *****************************************************//
 
 	if(Eta_file.is_open() && Eta_DYfile.is_open() && Eta_Polfile.is_open()){
 		for(unsigned i(0); i<nEtaBins; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			eta_den[i] = fitmean;
			EtaBins[i] = lowcut;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			eta_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 		for(unsigned i(0); i<nEtaBins; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			eta_num[i] = fitmean;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			eta_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 	}
	
	
 	for(unsigned i(0); i<nEtaBins; i++){
 		double fakerate = eta_num[i]/eta_den[i];
		std::cout << "eta " << EtaBins[i] << " " << fakerate << std::endl;
		if (i == 0) {
    			if (doEB == true) {
        			etafile << "double MCetaRatesEB_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    			} else {
        			etafile << "double MCetaRatesEE_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    			}
		} else if (i == nEtaBins - 1) 
    			etafile << fakerate << "};" << endl;
		 else 
    			etafile << fakerate << ", ";
		

 		
		double error = sqrt(fakerate*fakerate*eta_denerror[i]*eta_denerror[i]/(eta_den[i]*eta_den[i])+ eta_numerror[i]*eta_numerror[i]/eta_den[i]/eta_den[i]);
		double xvalue = (i+1 < nEtaBins)? (EtaBins[i]+EtaBins[i+1])/2.0 : EtaBins[i];
		double xerror = (i+1 < nEtaBins)? (EtaBins[i+1]-EtaBins[i])/2.0 : 0.005;
 		fr_bothcount_eta->SetPoint(i,xvalue, fakerate );
 		fr_bothcount_eta->SetPointError(i, xerror, error);
 		fr_eta_sigmaband->SetPoint(i,xvalue, fakerate );
 		fr_eta_sigmaband->SetPointError(i, xerror, error);
 	}
 
}
