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

double MINHT = 0;
double MAXHT = 1000;
bool doEB = true;
const char*processName ="Data"; bool doDrellYan = false;

Double_t fakerate_ptDependence(Double_t *x, Double_t *par)
{
	double slope = par[0];
	double constant = par[1]; 
	double index = par[2];
	double coeff = 1.0; 
	double pt = TMath::Max(x[0],0.000001);
	double arg = slope * pt + constant;
	if (arg <= 0) return 0.0;
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
 
 
double fakerate_binned(double x, const double *bins, const double *num, const double *den, unsigned nbins)
{
	for(unsigned i(0); i<nbins; ++i){
		if(i+1 < nbins){
			if(x > bins[i] && x < bins[i+1]) return num[i]/den[i];
		}
		else if(x > bins[i]) return num[i]/den[i];
	}
	return (nbins > 0 && den[0] != 0) ? num[0]/den[0] : 0;
}

 
void fitFakeFunc_HT(int RunYear, bool preVFP){//main
	gROOT->SetBatch(1);
	gSystem->Load("../../../lib/libAnaClasses.so");
        gSystem->Load("../../../lib/libRooFitClasses.so");
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	std::string whichVFP;
        if(RunYear==2016 and preVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	std::string isEB;
	if(doEB == true) isEB = "";
	else if(doEB == false) isEB = "_EE";

	std::ofstream resultfile, fakeRateFile, etafile, htfile;

	if (doEB == true) {
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPt-EB.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_fitFakeFunc_Data-HT-added.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPt-EB.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_eta_dependence_Data.txt", RunYear, whichVFP.c_str()));
        		htfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_ht_dependence_Data.txt", RunYear, whichVFP.c_str()));
    		}}
	else {
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPt-EE.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY_EE.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_fitFakeFunc_Data_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateFile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPt-EE.txt", RunYear, whichVFP.c_str()));
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_eta_dependence_Data_EE.txt", RunYear, whichVFP.c_str()));
		}}

	std::ifstream Pt_file, Pt_DYfile, Pt_Polfile, Eta_file, Eta_DYfile, Eta_Polfile, HT_file, HT_DYfile, HT_Polfile;

	if(doDrellYan==true) {
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
     		HT_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-ht-60-120.txt",RunYear,whichVFP.c_str()));
       		HT_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-ht-60-120.txt",RunYear,whichVFP.c_str()));
       		HT_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-expo-ht-60-120.txt",RunYear,whichVFP.c_str()));
	}
	else{
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));

		Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-eta-60-120.txt",RunYear,whichVFP.c_str()));

		HT_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-ht-60-120.txt",RunYear,whichVFP.c_str()));
       		HT_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-ht-60-120.txt",RunYear,whichVFP.c_str()));
       		HT_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-ht-60-120.txt",RunYear,whichVFP.c_str()));
	}

 	std::string line;
 	unsigned nPtBins(0);
 	unsigned nEtaBins(0);
	unsigned nHTBins(0);
 	while(std::getline(Pt_file , line))
 		nPtBins += 1;
 	std::cout << "Pt_file total line" << nPtBins << std::endl;
 	Pt_file.clear(); Pt_file.seekg(0, ios::beg);
 	while(std::getline(Eta_file , line))
 		nEtaBins += 1;
 	std::cout << "Eta_file total line" << nEtaBins << std::endl;
 	Eta_file.clear(); Eta_file.seekg(0, ios::beg);
	while(std::getline(HT_file , line))
		nHTBins += 1;
	std::cout << "HT_file total line" << nHTBins << std::endl;
	HT_file.clear(); HT_file.seekg(0, ios::beg);

 	nPtBins /= 2;
 	nEtaBins /= 2;
	nHTBins /= 2;
 
 	double PtBins[nPtBins];
 	double EtaBins[nEtaBins];
	double HTBins[nHTBins];

 	double pt_den[nPtBins];
 	double pt_denerror[nPtBins];
 	double pt_num[nPtBins];
 	double pt_numerror[nPtBins];
 
 	double eta_den[nEtaBins];
 	double eta_denerror[nEtaBins];
 	double eta_num[nEtaBins];
 	double eta_numerror[nEtaBins];
 
	double ht_den[nHTBins];
	double ht_denerror[nHTBins];
	double ht_num[nHTBins];
	double ht_numerror[nHTBins];
 	std::string bintype;
 	std::string numtype;
 	float lowcut;
 	float signal1, error1;
 	float fitmean, fitrms;
 	float DYmean,  DYrms;
 	float Polmean, Polrms;
 	TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins);
 	TGraphErrors *fr_bothcount_eta= new TGraphErrors(nEtaBins);
	TGraphErrors *fr_bothcount_ht= new TGraphErrors(nHTBins);
	TGraphErrors *fr_pt_ratio = new TGraphErrors(nPtBins);
	TGraphErrors *fr_pt_ratioError = new TGraphErrors(nPtBins);
	TGraphErrors *fr_ht_ratio = new TGraphErrors(nHTBins);
	TGraphErrors *fr_ht_ratioError = new TGraphErrors(nHTBins);
 
 	int graphPtBins = (int)((MAXPT - MINPT)/0.5);
	int graphHTBins = nHTBins;
 	TGraphErrors *fr_pt_sigmaband = new TGraphErrors(graphPtBins);
	TGraphErrors *fr_ht_sigmaband = new TGraphErrors(graphHTBins);
 	TGraphErrors *fr_eta_sigmaband = new TGraphErrors(nEtaBins-1);
 
//************************************ Fill and fit the pt dependence. Calculate ratio ********************************************************** 
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

	TF1 *f1 = new TF1("f1", fakerate_ptDependence,MINPT,MAXPT,3);
	TF1 *fit_fakerate_pt;
	TFitResultPtr result_fitpt;
	TVirtualFitter::SetMaxIterations(1000000);

	f1->SetParameter(0,1.0e+04);
        f1->SetParLimits(0, 1.001e+02, 1.25601e+04);
        f1->SetParameter(2,-0.4);
        f1->SetParLimits(2, -0.4, -0.1);
	
	if(!doDrellYan){
		if(RunYear==2016){
        		f1->SetParameter(1, -10000);
        		f1->SetParLimits(1, -15000.0, -5000);}
		if(RunYear==2017){
                	f1->SetParameter(1, -20000);
                	f1->SetParLimits(1, -39070.0, -10000);}
		if(RunYear==2018){
                	f1->SetParameter(1, -30000);
                	f1->SetParLimits(1, -38685.0, 0);}
	}	

		if(doDrellYan){
                f1->SetParameter(1, -38000);
                f1->SetParLimits(1, -38100.0, -10000);}

	f1->SetParNames("slope","constant","index");
        result_fitpt = fr_bothcount_pt->Fit("f1","R S");

	int status = result_fitpt->Status();        // 0 means OK
	int covQual = result_fitpt->CovMatrixStatus();  // 3 is the best

	std::cout << "Fit status = " << status << std::endl;
	std::cout << "CovQual = " << covQual << std::endl;

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

  //********************************* Fill and fit the eta dependence ********************************************************************************//
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
			if(doDrellYan == false){
    				if (doEB == true) 
        				etafile << "double etaRatesEB_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    				else if (doEB == false)
        				etafile << "double etaRatesEE_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";}
			else{
    				if (doEB == true) 
        				etafile << "double MCetaRatesEB_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    				else if (doEB == false)
        				etafile << "double MCetaRatesEE_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";}
    			
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

 //******************************** Fill and fit ht dependence **********************************************************************************//
	if(HT_file.is_open() && HT_DYfile.is_open() && HT_Polfile.is_open()){
		for(unsigned i(0); i<nHTBins; i++){
			HT_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			ht_den[i] = fitmean;
			HTBins[i] = lowcut;
			HT_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			HT_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			ht_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}
		for(unsigned i(0); i<nHTBins; i++){
			HT_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			ht_num[i] = fitmean;
			HT_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			HT_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			ht_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}
	}
	for(unsigned i(0); i<nHTBins; i++){
		double fakerate = ht_num[i]/ht_den[i];
		std::cout << "ht " << HTBins[i] << " " << fakerate << std::endl;
		if (i == 0) {
			if(doDrellYan == false){
    				if (doEB == true) 
        				htfile << "double htRatesEB_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    				else if (doEB == false)
        				htfile << "double htRatesEE_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";}
			else{
    				if (doEB == true) 
        				htfile << "double MChtRatesEB_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";
    				else if (doEB == false)
        				htfile << "double MChtRatesEE_" << RunYear << whichVFP << "[] = {" << fakerate << ", ";}
    			
		} else if (i == nHTBins - 1) 
    			htfile << fakerate << "};" << endl;
		 else 
    			htfile << fakerate << ", ";
		
		double error = sqrt(fakerate*fakerate*ht_denerror[i]*ht_denerror[i]/(ht_den[i]*ht_den[i])+ ht_numerror[i]*ht_numerror[i]/ht_den[i]/ht_den[i]);
		double xvalue = (i+1 < nHTBins)? (HTBins[i]+HTBins[i+1])/2.0 : HTBins[i];
	   double xerror = (i+1 < nHTBins)? (HTBins[i+1]-HTBins[i])/2.0 : 0.5;
		fr_bothcount_ht->SetPoint(i,xvalue, fakerate );
		fr_bothcount_ht->SetPointError(i, xerror, error);
		fr_ht_sigmaband->SetPoint(i,xvalue,fakerate);
		fr_ht_sigmaband->SetPointError(i,xerror,error);
		fr_ht_ratio->SetPoint(i,xvalue,1.0);
		fr_ht_ratio->SetPointError(i,xerror,(fakerate>0)?error/fakerate:0);
		fr_ht_ratioError->SetPoint(i,xvalue,1.0);
		fr_ht_ratioError->SetPointError(i,xerror,(fakerate>0)?error/fakerate:0);
	}

 // *****************************************************************************************************************************//

 	TH1D* invmass_den = new TH1D("invmass_den", "invmass_den",100,40,140);
 	TH1D* invmass_num = new TH1D("invmass_num", "invmass_num",100,40,140);
 	TH1D* h_invmass_bg  = new TH1D("invmass_bg",  "invmass_bg", 100,40,140);
 	TH1D* invmass_prednum = new TH1D("invmass_prednum", "invmass_prednum",100,40,140);

	TChain *etree = new TChain("FakeRateTree");
	etree->Add(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_%sTnP_dR05_%d%s.root",processName,RunYear,whichVFP.c_str()));

	float invmass=0;
	float tagPt=0;
	float probePt=0;
 	float probeEta=0;
 	bool  vetovalue=0;
 	bool  FSRveto = 0;
	float probeHT=0;
	etree->SetBranchAddress("invmass",   &invmass);
 	etree->SetBranchAddress("tagPt",     &tagPt);
 	etree->SetBranchAddress("probePt",   &probePt);
 	etree->SetBranchAddress("probeEta",  &probeEta);
 	etree->SetBranchAddress("vetovalue", &vetovalue);
	etree->SetBranchAddress("HT",        &probeHT);
 	etree->SetBranchAddress("FSRveto",   &FSRveto);

 	std::vector<float> etreeEt;
	std::vector<float> etreeHT;
 	std::vector<float> etreeEta;
 	std::vector<float> etreeInvmass;
 	etreeEt.clear();
	etreeHT.clear();
 	etreeEta.clear();
 	etreeInvmass.clear();
	for(unsigned iEvt(0); iEvt < etree->GetEntries(); iEvt++){
		etree->GetEntry(iEvt);

		if(probePt < 35)continue;
		if(doEB && fabs(probeEta) > 1.4442)continue;
		else if(!doEB && (fabs(probeEta) < 1.56 || fabs(probeEta) > 2.1))continue;

		if(vetovalue== false)invmass_den->Fill(invmass);
		else if(vetovalue== true && FSRveto == true)invmass_num->Fill(invmass);

		float weight_eta(0);
 		for(unsigned i(0); i<nEtaBins; i++){
			if(i+1 < nEtaBins){
				if(fabs(probeEta) > EtaBins[i] && fabs(probeEta) < EtaBins[i+1])
					weight_eta = eta_num[i]/eta_den[i];
			}
			
			else{
				if(fabs(probeEta) > EtaBins[i])
					weight_eta = eta_num[i]/eta_den[i];
			}
			
		}
		float weight_ht = fakerate_binned(probeHT, HTBins, ht_num, ht_den, nHTBins);
		float w_ele = fit_fakerate_pt->Eval(probePt)*weight_ht*weight_eta;
 		if(vetovalue==false){
 			invmass_prednum->Fill(invmass, w_ele);
 			etreeEt.push_back(probePt);
			etreeHT.push_back(probeHT);
 			etreeEta.push_back(fabs(probeEta));
 			etreeInvmass.push_back(invmass);
		}
 	}

 	TChain *bgtree = new TChain("BGTree");
 	bgtree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_bgtemplate_FullEcal.root");
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

	for(unsigned iEvt(0); iEvt < bgtree->GetEntries(); iEvt++){
		bgtree->GetEntry(iEvt);

		if(probePt_bg < 35)continue;
		if(doEB && fabs(probeEta_bg) > 1.4442)	continue;
		else if(!doEB && (fabs(probeEta_bg) < 1.56 || fabs(probeEta_bg) > 2.1))continue;
		if(vetovalue_bg== true && FSRveto_bg==true)	{h_invmass_bg->Fill(invmass_bg);}
	}

 	TH1D  *h_DYinvmass = new TH1D("h_DYinvmass","h_DYinvmass",80,70,110);
 	TChain *DYtree = new TChain("FakeRateTree");
	DYtree->Add(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_DYTnP_dR05_%d%s.root",RunYear,whichVFP.c_str()));

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
	DYtree->SetBranchAddress("invmassUncalib", &DY_invmass);
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

		if(DY_probePt < 35)continue;
		if(doEB && fabs(DY_probeEta) > 1.4442)continue;
		else if(!doEB && (fabs(DY_probeEta) < 1.56 || fabs(DY_probeEta) > 2.1))continue;

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
		if(isZee)h_DYinvmass->Fill(DY_invmass, 1);
	}

	RooRealVar mass_axis("invmass","invmass",70,110);
	TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
	c_fitMass->cd();

	RooDataHist datahist_data("both", "", mass_axis, invmass_num);
	RooDataHist datahist_bg("bg","",mass_axis, h_invmass_bg);
	RooHistPdf  pdf_bg("pdf_bg","pdf_bg",mass_axis, datahist_bg);
	RooHistPdf  pdf_bgAlter("pdf_bgAlter","pdf_bgAlter",mass_axis, datahist_bg);
	RooDataHist *datahist_DY = new RooDataHist("DYDataSet","DYDataSet", mass_axis, h_DYinvmass);
	RooHistPdf *DYpdf = new RooHistPdf("DYpdf","DYpdf", mass_axis, *datahist_DY);
	DYpdf->setInterpolationOrder(1);

  	RooRealVar m0( "m0", "m0", 91.188,85,95);
  	RooRealVar width( "width", "width", 2.495,0, 15);
	RooRealVar mean("mean", "" ,0.);
	RooRealVar sigma("sigma", "",2.4 , 0.0, 15.0);
	RooRealVar alpha("alpha", "", 1.0, 0.0, 20.0);
	RooRealVar n("n","", 1.0, 0.0, 20.0);
	RooRealVar alpha2("2ndalpha","", 1.0, -10.0, 20.0);
	RooRealVar n2("2ndn", "", 1.0, 0.0, 20.0);
	RooBreitWigner bw("bw", "", mass_axis, m0, width);
	RooDCBShape *cb;
	cb = new RooDCBShape("cb","cb", mass_axis, mean, sigma, alpha, n, alpha2, n2);
	RooGaussian gauss("gs", "gs", mass_axis, mean, sigma);
	RooFFTConvPdf signalRes("pdf", "pdf",mass_axis, bw, *cb);
  	RooFFTConvPdf signalResDY("pdfalter", "pdfalter",mass_axis, *DYpdf, gauss);

	double iniSig = 0.2*invmass_num->Integral(1,100);
	int    lowBinNumber = invmass_num->FindBin(70.0);
	int    highBinNumber= invmass_num->FindBin(110.0);
	double iniBkg = 40*(invmass_num->GetBinContent(lowBinNumber)+invmass_num->GetBinContent(highBinNumber));

	RooRealVar nSig("nSig", "", iniSig, 0, invmass_num->GetEntries()*1.2);
	RooRealVar nSigAlter("nSigAlter", "", iniSig, 0, invmass_num->GetEntries()*1.2);
	RooRealVar nBkg("nBkg", "", iniBkg, 0, invmass_num->GetEntries());
	RooRealVar nBkgAlter("nBkgAlter", "", iniBkg, 0, invmass_num->GetEntries());
	RooAddPdf *model = new RooAddPdf("modelalter", "", RooArgList(pdf_bg, signalResDY),RooArgList(nBkg, nSig));
	RooAddPdf *modelAlter = new RooAddPdf("model", "", RooArgList(pdf_bgAlter, signalRes),RooArgList(nBkgAlter, nSigAlter));

	RooPlot* mass_Frame = mass_axis.frame(RooFit::Title("totalNum"),RooFit::Bins(40));
	model->fitTo(datahist_data);
	modelAlter->fitTo(datahist_data);
	datahist_data.plotOn(mass_Frame);
	model->plotOn(mass_Frame, RooFit::Components(pdf_bg),
	RooFit::LineStyle(kDashed),
	RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	modelAlter->plotOn(mass_Frame, RooFit::Components(pdf_bgAlter),
	RooFit::LineStyle(kDashed),
	RooFit::LineColor(kRed),
	RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	model->plotOn(mass_Frame, RooFit::Components(RooArgSet(pdf_bg, signalResDY)),
	RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	modelAlter->plotOn(mass_Frame,RooFit::Components(RooArgSet(pdf_bgAlter, signalRes)),
	RooFit::LineColor(kRed),
	RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	mass_Frame->Draw();

	if(doDrellYan)
		c_fitMass->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/fit_totalNum_DrellYan_%d%s.png",RunYear,whichVFP.c_str()));
	else
		c_fitMass->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/fit_totalNum_data_%d%s.png",RunYear,whichVFP.c_str()));
	mass_axis.setRange("signal",70,110);
	RooAbsReal* igx_sig = signalResDY.createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	std::cout << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " error = " << igx_sig->getVal()*(nSig.getError()) <<  std::endl;
 	std::cout << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	std::cout << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
	resultfile << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " error = " << igx_sig->getVal()*(nSig.getError()) <<  std::endl;
 	resultfile << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	resultfile << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	fakeRateFile << "scalefactor " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;

	RooAbsReal* igx_sig_alter = signalRes.createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	std::cout << "num alter= " <<  igx_sig_alter->getVal()*(nSigAlter.getVal()) << " error = " << igx_sig_alter->getVal()*(nSigAlter.getError()) <<  std::endl;
 	std::cout << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	std::cout << "scale factor alter= " << igx_sig_alter->getVal()*(nSigAlter.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
	resultfile << "num alter= " <<  igx_sig_alter->getVal()*(nSigAlter.getVal()) << " error = " << igx_sig_alter->getVal()*(nSigAlter.getError()) <<  std::endl;
 	resultfile << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	resultfile << "scale factor alter= " << igx_sig_alter->getVal()*(nSigAlter.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;

	TCanvas *cancompare = new TCanvas("cancompare","",600,600);
	cancompare->cd();
	RooPlot* compare_Frame = mass_axis.frame(RooFit::Title("compare"),RooFit::Bins(40));
	model->plotOn(compare_Frame,
	RooFit::Components(signalResDY),
	RooFit::LineStyle(kDotted),
	RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	invmass_prednum->Scale(igx_sig->getVal()*(nSig.getVal())*1.0/invmass_prednum->Integral(lowBinNumber,highBinNumber));
	RooDataHist datahist_prednum("both", "", mass_axis, invmass_prednum);
	datahist_prednum.plotOn(compare_Frame, RooFit::MarkerColor(kRed));
	compare_Frame->Draw();
	if(doDrellYan)cancompare->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/compare_predvsnum_DrellYan_%d%s.png",RunYear,whichVFP.c_str()));
	else cancompare->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/compare_predvsnum_data_%d%s.png",RunYear,whichVFP.c_str()));

 
 //*************************************   Toy MC **************************************************************************//

 	double totalnumError = sqrt(pow(igx_sig->getVal()*(nSig.getVal()) - igx_sig_alter->getVal()*(nSigAlter.getVal()), 2) + pow(igx_sig->getVal()*(nSig.getError()), 2));
   resultfile << "totalnum " << igx_sig->getVal()*(nSig.getVal()) <<  "  Error " << totalnumError << std::endl;
 	gRandom = new TRandom3(0);
 	gRandom->SetSeed(0);
	double random_totalnum[NTOY];
	for(unsigned ir(0); ir<NTOY; ir++)
 		random_totalnum[ir] = igx_sig->getVal()*(nSig.getVal()) + totalnumError*(-1+ gRandom->Rndm()*2.0);

	TMatrixDSym cov_pt = result_fitpt->GetCovarianceMatrix();
	result_fitpt->Print("V");
	TVectorD mean_pt(3) ;
	mean_pt(0) = result_fitpt->Parameter(0);
	mean_pt(1) = result_fitpt->Parameter(1);
	mean_pt(2) = result_fitpt->Parameter(2);
	RooRealVar central_slope("central_slope","",mean_pt(0)-result_fitpt->ParError(0), mean_pt(0)+result_fitpt->ParError(0));
	RooRealVar central_const("central_const","",mean_pt(1)-result_fitpt->ParError(1),mean_pt(1)+result_fitpt->ParError(1));
	RooRealVar central_index("central_index","",mean_pt(2)-result_fitpt->ParError(2),mean_pt(2)+result_fitpt->ParError(2));
			resultfile << "central_slope " << mean_pt(0) << " err " << result_fitpt->ParError(0) << std::endl;
			resultfile << "central_const " << mean_pt(1) << " err " << result_fitpt->ParError(1) << std::endl;
			resultfile << "central_index " << mean_pt(2) << " err " << result_fitpt->ParError(2) << std::endl;

			fakeRateFile << "ptslope " << mean_pt(0) <<  std::endl;
                        fakeRateFile << "ptconstant " << mean_pt(1)  << std::endl;
                        fakeRateFile << "ptindex " << mean_pt(2)  << std::endl;
                        fakeRateFile << "ptcoeff " << " 1.0 " << std::endl;

	RooMultiVarGaussian mvg_pt("mvg_pt","mvg_pt",RooArgList(central_slope,central_const,central_index),mean_pt,cov_pt);
	RooDataSet* toymcdata_pt = mvg_pt.generate(RooArgSet(central_slope,central_const,central_index),NTOY);
	std::ostringstream toymchistname;
	TF1 *h_toymc_pt[NTOY];

	ofstream myfile;
	if(doEB && doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_EB.txt",RunYear,whichVFP.c_str()));
        else if(!doEB && doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_FullEcal.txt",RunYear,whichVFP.c_str()));
        else if(doEB && !doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB_HT-added.txt",RunYear,whichVFP.c_str()));
        else if(!doEB && !doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_FullEcal.txt",RunYear,whichVFP.c_str()));
	TH1D *p_scalefactor = new TH1D("p_scalefactor","fake rate scale factor; scale factor;",100,1000,2000);

	float toyptvalue[graphPtBins][NTOY];
	float lowtoyptvalue[graphPtBins];
	float hightoyptvalue[graphPtBins];
	for(unsigned ii(0); ii < graphPtBins; ii++){
		lowtoyptvalue[ii] = 1;
		hightoyptvalue[ii] = 0;
	}
	float toyhtvalue[graphHTBins][NTOY];
	float lowtoyhtvalue[graphHTBins];
	float hightoyhtvalue[graphHTBins];
	for(unsigned ii(0); ii < graphHTBins; ii++){
		lowtoyhtvalue[ii] = 1;
		hightoyhtvalue[ii] = 0;
	}

	for(int i(0); i<NTOY; i++){
		double data1 = toymcdata_pt->get(i)->getRealValue("central_slope");
		double data2 = toymcdata_pt->get(i)->getRealValue("central_const");
		double data3 = toymcdata_pt->get(i)->getRealValue("central_index");
		if(data1 > 1e6 || data2 > 1e6 || data3 > 1e6)continue;
		toymchistname.str("");
		toymchistname << "h_toymc_pt_" << i;
		h_toymc_pt[i] = new TF1(toymchistname.str().c_str(), fakerate_ptDependence, 30,1000, 3);
		h_toymc_pt[i]->SetParameter(0, data1);
		h_toymc_pt[i]->SetParameter(1, data2);
		h_toymc_pt[i]->SetParameter(2, data3);


		invmass_prednum->Reset();
		for(unsigned iEvt(0); iEvt < etreeEt.size(); iEvt++){
 					float weight_eta(0);
 			 		for(unsigned i(0); i<nEtaBins; i++){
 						if(i+1 < nEtaBins){
 							if(fabs(etreeEta[iEvt]) > EtaBins[i] && fabs(etreeEta[iEvt]) < EtaBins[i+1])
								weight_eta = eta_num[i]/eta_den[i];}
						else{
 							if(fabs(etreeEta[iEvt]) > EtaBins[i])
								weight_eta = eta_num[i]/eta_den[i];
 						}
 					}
				float toy_weight_ht = fakerate_binned(etreeHT[iEvt], HTBins, ht_num, ht_den, nHTBins);
				float w_ele = h_toymc_pt[i]->Eval(etreeEt[iEvt])*toy_weight_ht*weight_eta;
				invmass_prednum->Fill(etreeInvmass[iEvt], w_ele);
			}

		if(invmass_prednum->Integral(lowBinNumber,highBinNumber) > 0 && invmass_prednum->Integral(lowBinNumber,highBinNumber) < 1e20){
				myfile << random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber) << " " << data1 << " " << data2 << " " << data3 << std::endl;
				p_scalefactor->Fill(random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber));
			for(unsigned ibin(1); ibin <= graphPtBins; ibin++){
				double estimated = h_toymc_pt[i]->Eval(MINPT+0.5*ibin);
				toyptvalue[ibin-1][i] = estimated;
				if(estimated < lowtoyptvalue[ibin-1])lowtoyptvalue[ibin-1] = estimated;
				if(estimated > hightoyptvalue[ibin-1])hightoyptvalue[ibin-1] = estimated;
			}
			for(unsigned ibin(0); ibin < nHTBins; ibin++){
					double estimated = ht_num[ibin]/ht_den[ibin];
					toyhtvalue[ibin][i] = estimated;
					if(estimated < lowtoyhtvalue[ibin])lowtoyhtvalue[ibin] = estimated;
					if(estimated > hightoyhtvalue[ibin])hightoyhtvalue[ibin] = estimated;
				}
			}
		}

      // *************************   Calculated errors  ***************************************************************************************************************//
	for(unsigned ibin(0); ibin < graphPtBins; ibin++){
		TH1D *h_toyptdis = new TH1D("h_toyptdis","",50,lowtoyptvalue[ibin],hightoyptvalue[ibin]);
		for(unsigned i(0); i < NTOY; i++)h_toyptdis->Fill(toyptvalue[ibin][i]);
		h_toyptdis->Fit("gaus");
		//float fiterror = h_toyptdis->GetFunction("gaus")->GetParameter(2);
		float fiterror = fabs(hightoyptvalue[ibin] - lowtoyptvalue[ibin])/2;
		fr_pt_sigmaband->SetPoint(ibin, MINPT+0.5*ibin,fit_fakerate_pt->Eval(MINPT+0.5*ibin));
		fr_pt_sigmaband->SetPointError(ibin, 0.25, fiterror);
		fr_pt_ratioError->SetPoint(ibin, MINPT+0.5*ibin, 1);
		fr_pt_ratioError->SetPointError(ibin, 0.25, fiterror/fit_fakerate_pt->Eval(MINPT+0.5*ibin));
		delete h_toyptdis;
	}
	    // HT is used bin-by-bin; keep direct bin errors and no HT fit uncertainty propagation.

	fr_bothcount_pt->SetMarkerStyle(20);
	fr_bothcount_pt->SetMarkerColor(kBlue);
	fr_bothcount_pt->SetLineColor(kBlue);
	fr_bothcount_pt->SetLineWidth(2);
	fr_bothcount_pt->SetFillColor(0);
	fr_bothcount_eta->SetMarkerStyle(20);
	fr_bothcount_eta->SetMarkerColor(kBlue);
	fr_bothcount_eta->SetLineColor(kBlue);
	fr_bothcount_eta->SetLineWidth(2);
	fr_bothcount_eta->SetFillColor(0);
	fr_bothcount_ht->SetMarkerStyle(20);
	fr_bothcount_ht->SetMarkerColor(kBlue);
	fr_bothcount_ht->SetLineColor(kBlue);
	fr_bothcount_ht->SetLineWidth(2);
	fr_bothcount_ht->SetFillColor(0);
	fr_pt_ratio->SetMarkerStyle(20);
	fr_pt_ratio->SetMarkerColor(kBlue);
	fr_pt_ratio->SetLineColor(kBlue);
	fr_pt_ratio->SetLineWidth(2);
	fr_pt_ratio->SetFillColor(0);
	fr_ht_ratio->SetMarkerStyle(20);
	fr_ht_ratio->SetMarkerColor(kBlue);
	fr_ht_ratio->SetLineColor(kBlue);
	fr_ht_ratio->SetLineWidth(2);
	fr_ht_ratio->SetFillColor(0);
		TLine *flatratio_pt = new TLine(30,1,200,1);
		TLine *flatratio_eta= new TLine(0,1,1.5,1);
		TLine *flatratio_ht= new TLine(MINHT,1,MAXHT,1);

	TLegend *leg = new TLegend(0.5,0.6,0.89,0.89);
	leg->AddEntry(fr_bothcount_pt, "data");
	leg->AddEntry(fit_fakerate_pt,"fit result");
	leg->AddEntry(fr_pt_sigmaband,"total uncertainty");

	TCanvas *canpt = new TCanvas("canpt","",600,600);
	canpt->cd();
	TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.3, 1, 1.0);
	canpt_pad1->SetBottomMargin(0);
	canpt_pad1->Draw();
	canpt_pad1->cd();
	TH1D *dummy_pt = new TH1D("",";Pt(GeV);fake rate",nPtBins*10,MINPT,MAXPT);
	dummy_pt->SetMaximum(0.06);
	dummy_pt->SetMinimum(0.0);
	dummy_pt->GetYaxis()->SetTitleOffset(1.4);
	dummy_pt->Draw();
	fr_bothcount_pt->Draw("EP same");
	fit_fakerate_pt->Draw("same");
	fr_pt_sigmaband->SetFillColor(15);
	fr_pt_sigmaband->SetFillStyle(3345);
	fr_pt_sigmaband->Draw("E2 same");
	fr_bothcount_pt->Draw("EP same");
	leg->Draw("same");
	if(RunYear==2016 and preVFP == 1)       CMS_lumi( canpt_pad1, 1, 1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canpt_pad1, 2, 1, 11 );
        else if(RunYear==2017)                  CMS_lumi( canpt_pad1, 3, 1, 11 );
        else if(RunYear==2018)                  CMS_lumi( canpt_pad1, 4, 1, 11 );

	canpt->cd();
	TPad *canpt_pad2 = new TPad("canpt_pad2", "pad2", 0, 0, 1, 0.3);
	canpt_pad2->SetTopMargin(0);
	canpt_pad2->SetBottomMargin(0.5);
	canpt_pad2->Draw();
	canpt_pad2->cd();
	TH1D *dummy_ptratio = new TH1D("dummy_ptratio",";p_{T}(GeV);data/fit",17,30,200);
	dummy_ptratio->SetMaximum(2);
	dummy_ptratio->SetMinimum(0);
 			dummy_ptratio->GetYaxis()->SetNdivisions(504);
	dummy_ptratio->Draw();
	fr_pt_ratio->Draw("EP same");
	fr_pt_ratioError->SetFillStyle(3345);
	fr_pt_ratioError->Draw("E2 same");
	fr_pt_ratio->Draw("EP same");
 	flatratio_pt->Draw("same");
	if(doDrellYan)canpt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_pt_systematic_DrellYan_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));
	else canpt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_pt_systematic_data_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));


	TCanvas *canht = new TCanvas("canht","",600,600);
	canht->cd();

	TPad *canht_pad1 = new TPad("canht_pad1","pad1",0,0,1,1.0);
	canht_pad1->SetBottomMargin(0.12);
	canht_pad1->SetLeftMargin(0.14);
	canht_pad1->SetRightMargin(0.04);
	canht_pad1->Draw();
	canht_pad1->cd();

	TH1D *dummy_ht = new TH1D("dummy_ht",";H_{T} (GeV);fake rate",50,MINHT,MAXHT);
	dummy_ht->SetMaximum(0.06);
	dummy_ht->SetMinimum(0.0);
	dummy_ht->GetYaxis()->SetTitleOffset(1.4);
	dummy_ht->Draw();

	fr_bothcount_ht->Draw("EP same");

	fr_ht_sigmaband->SetFillStyle(3345);
	fr_ht_sigmaband->Draw("E2 same");
	fr_bothcount_ht->Draw("EP same");
	leg->Draw("same");
	if(RunYear==2016 and preVFP == 1)       CMS_lumi( canht_pad1, 1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canht_pad1, 2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( canht_pad1, 3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( canht_pad1, 4,1, 11 );

	if(doDrellYan)canht->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_ht_systematic_DrellYan_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));
	else canht->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_ht_systematic_data_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));

	TLegend *legEta = new TLegend(0.5,0.7,0.89,0.89);
	legEta->AddEntry(fr_bothcount_eta, "data");

	TCanvas *caneta = new TCanvas("caneta","",600,600);
	caneta->cd();
	TPad *caneta_pad1 = new TPad("caneta_pad1", "pad1", 0, 0.05, 1, 1.0);
	caneta_pad1->SetBottomMargin(0.15);
	caneta_pad1->Draw();
	caneta_pad1->cd();
	TH1D *dummy_eta = new TH1D("",";|eta|;fake rate",30,0,1.5);
	dummy_eta->SetMaximum(0.08);
	dummy_eta->GetXaxis()->SetTitleOffset(0.8);
	dummy_eta->Draw();
	fr_bothcount_eta->Draw("EPL same");
	fr_eta_sigmaband->SetFillColor(15);
	fr_eta_sigmaband->SetFillStyle(3345);
	fr_eta_sigmaband->Draw("E2 same");
	fr_bothcount_eta->Draw("EPL same");
	legEta->Draw("same");
	if(RunYear==2016 and preVFP == 1)       CMS_lumi( caneta_pad1, 1, 1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( caneta_pad1, 2, 1, 11 );
        else if(RunYear==2017)                  CMS_lumi( caneta_pad1, 3, 1, 11 );
        else if(RunYear==2018)                  CMS_lumi( caneta_pad1, 4, 1, 11 );

	if(doDrellYan)caneta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_systematic_DrellYan_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));
	else caneta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_systematic_data_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));

	TCanvas *canscale = new TCanvas("scale","scale",600,600);
	p_scalefactor->Draw("hist");
	if(doDrellYan)canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_scalefactor_DrellYan_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));
	else canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_scalefactor_data_%d%s%s.png",RunYear,whichVFP.c_str(),isEB.c_str()));
 	resultfile.close();
}
