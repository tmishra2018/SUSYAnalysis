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
#define NTOY 10

#define MINETA 0
#define MAXETA 1.4
#define MINPT 30
#define MAXPT 200
#define MINVTX 0
#define MAXVTX 46
bool Sys = false;
bool doEB = true;


void FakeRates_Data_vs_DY(int RunYear, bool preVFP){//main 
        gROOT->SetBatch(1);
        std::string whichVFP;
        if(RunYear==2016 and preVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
            
	 
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	/* DY files */
	std::ifstream Pt_file(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Pt_DYfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Pt_Polfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Eta_file(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Eta_DYfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Eta_Polfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Vtx_file(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Vtx_DYfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-DY-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
        std::ifstream Vtx_Polfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DY-Bw-expo-vtx-60-120.txt",RunYear,whichVFP.c_str()));
	
	std::string line;
	unsigned nPtBins(0);
	unsigned nEtaBins(0);
	unsigned nVtxBins(0); 
	while(std::getline(Pt_file , line))
		nPtBins += 1;
	//std::cout << "Pt_file total line " << nPtBins << std::endl;
	Pt_file.clear(); Pt_file.seekg(0, ios::beg);

	while(std::getline(Eta_file , line))
		nEtaBins += 1;
	//std::cout << "Eta_file total line " << nEtaBins << std::endl;
	Eta_file.clear(); Eta_file.seekg(0, ios::beg);

	while(std::getline(Vtx_file , line))
                nVtxBins += 1;
        //std::cout << "Vtx_file total line " << nVtxBins << std::endl;
        Vtx_file.clear(); Vtx_file.seekg(0, ios::beg);
	
	nPtBins /= 2;
	nEtaBins /= 2;
	nVtxBins /= 2;
      
        //cout<<nPtBins <<" "<<nEtaBins<<" "<<nVtxBins<<" "<<endl;
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

	//************************************ Fill and fit the pt dependence. Calculate ratio **********************************************************// 

	if(Pt_file.is_open() && Pt_DYfile.is_open() && Pt_Polfile.is_open()){
		for(unsigned i(0); i<nPtBins; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			PtBins[i] = lowcut;
			pt_den[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			if(Sys==true) pt_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			if(Sys==false) pt_denerror[i] = sqrt(fitrms*fitrms);
		}    
		for(unsigned i(0); i<nPtBins; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			pt_num[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			if(Sys==true) pt_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			if(Sys==false) pt_numerror[i] = sqrt(fitrms*fitrms);
		}    
	}
	for(unsigned i(0); i<nPtBins; i++){
		double fakerate = pt_num[i]/pt_den[i];
		//std::cout << "pt " << PtBins[i] << " " << fakerate << std::endl;
		double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
		double xvalue = (i+1 < nPtBins)? (PtBins[i]+PtBins[i+1])/2.0 : (MAXPT + PtBins[i])/2.0;
		double xerror = (i+1 < nPtBins)? (PtBins[i+1]-PtBins[i])/2.0 : (MAXPT - PtBins[i])/2.0;
		fr_bothcount_pt->SetPoint(i, xvalue, fakerate );
		fr_bothcount_pt->SetPointError(i, xerror, error);
	}
	//********************************* Fill and fit the eta dependence ***********/	
	if(Eta_file.is_open() && Eta_DYfile.is_open() && Eta_Polfile.is_open()){
 		for(unsigned i(0); i<nEtaBins; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			eta_den[i] = fitmean;
			EtaBins[i] = lowcut;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			if(Sys==true) eta_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 			if(Sys==false) eta_denerror[i] = sqrt(fitrms*fitrms);
 		}    
 		for(unsigned i(0); i<nEtaBins; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			eta_num[i] = fitmean;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			if(Sys==true) eta_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 			if(Sys==false) eta_numerror[i] = sqrt(fitrms*fitrms);
 		}    
 	}
 	for(unsigned i(0); i<nEtaBins; i++){
 		double fakerate = eta_num[i]/eta_den[i];
		//std::cout << "eta " << EtaBins[i] << " " << fakerate << std::endl;
 		double error = sqrt(fakerate*fakerate*eta_denerror[i]*eta_denerror[i]/(eta_den[i]*eta_den[i])+ eta_numerror[i]*eta_numerror[i]/eta_den[i]/eta_den[i]);
		double xvalue = (i+1 < nEtaBins)? (EtaBins[i]+EtaBins[i+1])/2.0 : (MAXETA + EtaBins[i])/2.0;
                double xerror = (i+1 < nEtaBins)? (EtaBins[i+1]-EtaBins[i])/2.0 : (MAXETA - EtaBins[i])/2.0;
                fr_bothcount_eta->SetPoint(i, xvalue, fakerate );
                fr_bothcount_eta->SetPointError(i, xerror, error);

 	}    

	//******************************** Fill and fit vtx dependence **************************/
        if(Vtx_file.is_open() && Vtx_DYfile.is_open() && Vtx_Polfile.is_open()){
 		for(unsigned i(0); i<nVtxBins; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_den[i] = fitmean;
			VtxBins[i] = lowcut;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			if(Sys==true) vtx_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 			if(Sys==false) vtx_denerror[i] = sqrt(fitrms*fitrms);
 		}    
 		for(unsigned i(0); i<nVtxBins; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_num[i] = fitmean;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			if(Sys==true) vtx_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 			if(Sys==false) vtx_numerror[i] = sqrt(fitrms*fitrms);
 		}    
 	}
 	for(unsigned i(0); i<nVtxBins; i++){
 		double fakerate = vtx_num[i]/vtx_den[i];
		//std::cout << "vtx " << VtxBins[i] << " " << fakerate << std::endl;
 		double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
		double xvalue = (i+1 < nVtxBins)? (VtxBins[i]+VtxBins[i+1])/2.0 : (MAXVTX + VtxBins[i])/2.0;
    		double xerror = (i+1 < nVtxBins)? (VtxBins[i+1]-VtxBins[i])/2.0 : (MAXVTX - VtxBins[i])/2.0;
 		fr_bothcount_vtx->SetPoint(i,xvalue, fakerate );
 		fr_bothcount_vtx->SetPointError(i, xerror, error);
 	}

        		
	/* Data files */
	std::ifstream Pt_fileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Pt_DYfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Pt_PolfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-pt-60-120.txt",RunYear,whichVFP.c_str()));

	std::ifstream Eta_fileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Eta_DYfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Eta_PolfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-eta-60-120.txt",RunYear,whichVFP.c_str()));

	std::ifstream Vtx_fileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Vtx_DYfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-DY-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
	std::ifstream Vtx_PolfileData(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-Bw-expo-vtx-60-120.txt",RunYear,whichVFP.c_str()));
	

	std::string lineData;
	unsigned nPtBins_Data(0);
	unsigned nEtaBins_Data(0);
	unsigned nVtxBins_Data(0); 
	while(std::getline(Pt_fileData , lineData))
		nPtBins_Data += 1;
	//std::cout << "Pt_fileData total line " << nPtBins_Data << std::endl;
	Pt_fileData.clear(); Pt_fileData.seekg(0, ios::beg);
	while(std::getline(Eta_fileData , lineData))
		nEtaBins_Data += 1;
	//std::cout << "Eta_fileData total line " << nEtaBins_Data << std::endl;
	Eta_fileData.clear(); Eta_fileData.seekg(0, ios::beg);
	while(std::getline(Vtx_fileData , lineData))
		nVtxBins_Data += 1;
	//std::cout << "Vtx_fileData total line " << nVtxBins_Data << std::endl;
	Vtx_fileData.clear(); Vtx_fileData.seekg(0, ios::beg);
	
	nPtBins_Data /= 2;
	nEtaBins_Data /= 2;
	nVtxBins_Data /= 2;
      
        //cout<<nPtBins_Data <<" "<<nEtaBins_Data<<" "<<nVtxBins_Data<<" "<<endl;
	double PtBins_Data[nPtBins_Data];
	double EtaBins_Data[nEtaBins_Data];
	double VtxBins_Data[nVtxBins_Data];
	
	double ptData_den[nPtBins_Data];
	double ptData_denerror[nPtBins_Data];
	double ptData_num[nPtBins_Data];
	double ptData_numerror[nPtBins_Data];

	double etaData_den[nEtaBins_Data];
	double etaData_denerror[nEtaBins_Data];
	double etaData_num[nEtaBins_Data];
	double etaData_numerror[nEtaBins_Data];

	double vtxData_den[nVtxBins_Data];
	double vtxData_denerror[nVtxBins_Data];
	double vtxData_num[nVtxBins_Data];
	double vtxData_numerror[nVtxBins_Data];
	std::string bintypeData;
	std::string numtypeData;
	float lowcutData;
	float signal1_Data, error1_Data;
	float fitmeanData, fitrmsData;
	float DYmeanData,  DYrmsData;
	float PolmeanData, PolrmsData;
	TGraphErrors *fr_bothcount_ptData = new TGraphErrors(nPtBins_Data);
	TGraphErrors *fr_bothcount_etaData= new TGraphErrors(nEtaBins_Data);
	TGraphErrors *fr_bothcount_vtxData= new TGraphErrors(nVtxBins_Data);
	TGraphErrors *fr_ptData_ratio = new TGraphErrors(nVtxBins_Data-1);
	TGraphErrors *fr_ptData_ratioError = new TGraphErrors(nVtxBins_Data-1);
	TGraphErrors *fr_vtxData_ratio = new TGraphErrors(nVtxBins_Data-1);
	TGraphErrors *fr_vtxData_ratioError = new TGraphErrors(nVtxBins_Data-1);

	int graphPtBins_Data = (int)((MAXPT - MINPT)/0.5);
	int graphVtxBins_Data = (int)((MAXVTX - MINVTX)/0.5);
	TGraphErrors *fr_ptData_sigmaband = new TGraphErrors(graphPtBins_Data);
	TGraphErrors *fr_vtxData_sigmaband = new TGraphErrors(graphVtxBins_Data);
	TGraphErrors *fr_etaData_sigmaband = new TGraphErrors(nEtaBins_Data-1);
	//************************************ Fill and fit the ptData dependence. Calculate ratio **********************************************************// 

	if(Pt_fileData.is_open() && Pt_DYfileData.is_open() && Pt_PolfileData.is_open()){
		for(unsigned i(0); i<nPtBins_Data; i++){ 
			Pt_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			PtBins_Data[i] = lowcutData;
			ptData_den[i] = fitmeanData;
			Pt_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Pt_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) ptData_denerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) ptData_denerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
		for(unsigned i(0); i<nPtBins_Data; i++){ 
			Pt_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			ptData_num[i] = fitmeanData;
			Pt_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Pt_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) ptData_numerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) ptData_numerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
	}
	for(unsigned i(0); i<nPtBins_Data; i++){
		double fakerate = ptData_num[i]/ptData_den[i];
		//std::cout << "pt " << PtBins_Data[i] << " " << fakerate << std::endl;
		double error = sqrt(fakerate*fakerate*ptData_denerror[i]*ptData_denerror[i]/(ptData_den[i]*ptData_den[i])+ ptData_numerror[i]*ptData_numerror[i]/ptData_den[i]/ptData_den[i]);
		double xvalue = (i+1 < nPtBins_Data)? (PtBins_Data[i]+PtBins_Data[i+1])/2.0 : (MAXPT + PtBins_Data[i])/2.0;
		double xerror = (i+1 < nPtBins_Data)? (PtBins_Data[i+1]-PtBins_Data[i])/2.0 : (MAXPT - PtBins_Data[i])/2.0;
		fr_bothcount_ptData->SetPoint(i, xvalue, fakerate );
		fr_bothcount_ptData->SetPointError(i, xerror, error);
	}

	if(Eta_fileData.is_open() && Eta_DYfileData.is_open() && Eta_PolfileData.is_open()){
		for(unsigned i(0); i<nEtaBins_Data; i++){ 
			Eta_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			EtaBins_Data[i] = lowcutData;
			etaData_den[i] = fitmeanData;
			Eta_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Eta_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) etaData_denerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) etaData_denerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
		for(unsigned i(0); i<nEtaBins_Data; i++){ 
			Eta_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			etaData_num[i] = fitmeanData;
			Eta_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Eta_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) etaData_numerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) etaData_numerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
	}
	for(unsigned i(0); i<nEtaBins_Data; i++){
		double fakerate = etaData_num[i]/etaData_den[i];
		//std::cout << "eta " << EtaBins_Data[i] << " " << fakerate << std::endl;
		double error = sqrt(fakerate*fakerate*etaData_denerror[i]*etaData_denerror[i]/(etaData_den[i]*etaData_den[i])+ etaData_numerror[i]*etaData_numerror[i]/etaData_den[i]/etaData_den[i]);
		double xvalue = (i+1 < nEtaBins_Data)? (EtaBins_Data[i]+EtaBins_Data[i+1])/2.0 : (MAXETA + EtaBins_Data[i])/2.0;
		double xerror = (i+1 < nEtaBins_Data)? (EtaBins_Data[i+1]-EtaBins_Data[i])/2.0 : (MAXETA - EtaBins_Data[i])/2.0;
		fr_bothcount_etaData->SetPoint(i, xvalue, fakerate );
		fr_bothcount_etaData->SetPointError(i, xerror, error);
	}

	if(Vtx_fileData.is_open() && Vtx_DYfileData.is_open() && Vtx_PolfileData.is_open()){
		for(unsigned i(0); i<nVtxBins_Data; i++){ 
			Vtx_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			VtxBins_Data[i] = lowcutData;
			vtxData_den[i] = fitmeanData;
			Vtx_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Vtx_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) vtxData_denerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) vtxData_denerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
		for(unsigned i(0); i<nVtxBins_Data; i++){ 
			Vtx_fileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> fitmeanData >> fitrmsData;
			vtxData_num[i] = fitmeanData;
			Vtx_DYfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> DYmeanData >> DYrmsData;
			Vtx_PolfileData >> bintypeData >> numtypeData >> lowcutData >> signal1_Data >> error1_Data >> PolmeanData >> PolrmsData;
			double sysdiff = fabs(DYmeanData - fitmeanData) > fabs(PolmeanData - fitmeanData)? fabs(DYmeanData - fitmeanData):fabs(PolmeanData - fitmeanData);
			if(Sys==true) vtxData_numerror[i] = sqrt(fitrmsData*fitrmsData + sysdiff*sysdiff);
			if(Sys==false) vtxData_numerror[i] = sqrt(fitrmsData*fitrmsData);
		}    
	}
	for(unsigned i(0); i<nVtxBins_Data; i++){
		double fakerate = vtxData_num[i]/vtxData_den[i];
		//std::cout << "vtx " << VtxBins_Data[i] << " " << fakerate << std::endl;
		double error = sqrt(fakerate*fakerate*vtxData_denerror[i]*vtxData_denerror[i]/(vtxData_den[i]*vtxData_den[i])+ vtxData_numerror[i]*vtxData_numerror[i]/vtxData_den[i]/vtxData_den[i]);
		double xvalue = (i+1 < nVtxBins_Data)? (VtxBins_Data[i]+VtxBins_Data[i+1])/2.0 : (MAXVTX + VtxBins_Data[i])/2.0;
		double xerror = (i+1 < nVtxBins_Data)? (VtxBins_Data[i+1]-VtxBins_Data[i])/2.0 : (MAXVTX - VtxBins_Data[i])/2.0;
		fr_bothcount_vtxData->SetPoint(i, xvalue, fakerate );
		fr_bothcount_vtxData->SetPointError(i, xerror, error);
	}
			
			TCanvas *can = new TCanvas("can","",1200,1000);
			can->SetRightMargin(0.1);
			can->SetLeftMargin(0.13);
			can->cd();

                       	fr_bothcount_pt->SetMinimum(0);
                       	fr_bothcount_pt->SetMaximum(0.03);
			
                        fr_bothcount_pt->GetXaxis()->SetTitle("Pt(GeV) ");
                        fr_bothcount_pt->GetYaxis()->SetTitle("fake rate");
			fr_bothcount_pt->GetXaxis()->SetTitleSize(.05);
			fr_bothcount_pt->GetYaxis()->SetTitleSize(.05);
                       // fr_bothcount_pt->SetTitle("e #rightarrow #gamma fake rate");
                        fr_bothcount_pt->SetTitle("");
                      
                        fr_bothcount_pt->SetMarkerColor(kRed);
                        fr_bothcount_pt->SetLineColor(kRed);
                        fr_bothcount_ptData->SetMarkerColor(kBlack);
                        fr_bothcount_ptData->SetLineColor(kBlack);
			fr_bothcount_pt->SetMarkerStyle(kCircle);
			fr_bothcount_ptData->SetMarkerStyle(kStar);
			
			fr_bothcount_pt->SetMarkerSize(2);
			fr_bothcount_ptData->SetMarkerSize(3);
                        
			fr_bothcount_ptData->SetMaximum(0.05);
			fr_bothcount_ptData->SetMinimum(0.0);
			fr_bothcount_pt->Draw("Ap");
			fr_bothcount_ptData->Draw("p");
			
			TLegend *leg = new TLegend(0.71,0.69,0.88,0.84);
			leg->AddEntry(fr_bothcount_ptData,"Data");
			leg->AddEntry(fr_bothcount_pt, "Drell-Yan ");
			leg->Draw();

			if(RunYear==2016 and preVFP == 1)       CMS_lumi( can, 1, 1, 11 );
	                else if(RunYear==2016 and preVFP == 0)  CMS_lumi( can, 2, 1, 11 );
        	        else if(RunYear==2017)                  CMS_lumi( can, 3, 1, 11 );
                	else if(RunYear==2018)                  CMS_lumi( can, 4, 1, 11 );
			can->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_pt_data_vs_DY_%d%s.pdf",RunYear,whichVFP.c_str()));

			TCanvas *canEta = new TCanvas("canEta","",1200,1000);
			canEta->SetRightMargin(0.04);
			canEta->SetLeftMargin(0.14);
			canEta->cd();

                        fr_bothcount_eta->GetXaxis()->SetTitle("Eta ");
                        fr_bothcount_eta->GetYaxis()->SetTitle("fake rate");
			fr_bothcount_eta->GetXaxis()->SetTitleSize(.05);
			fr_bothcount_eta->GetYaxis()->SetTitleSize(.05);
                        fr_bothcount_eta->SetTitle("");
                      
                        fr_bothcount_eta->SetMarkerColor(kRed);
                        fr_bothcount_eta->SetLineColor(kRed);
                        fr_bothcount_etaData->SetMarkerColor(kBlack);
                        fr_bothcount_etaData->SetLineColor(kBlack);
			fr_bothcount_eta->SetMarkerStyle(kCircle);
			fr_bothcount_etaData->SetMarkerStyle(kStar);
			
			fr_bothcount_eta->SetMarkerSize(2);
			fr_bothcount_etaData->SetMarkerSize(3);
                        
			fr_bothcount_eta->SetMaximum(0.08);
			fr_bothcount_eta->SetMinimum(0);
			fr_bothcount_eta->Draw("Ap");
			fr_bothcount_etaData->Draw("p");
			TLegend *leg2 = new TLegend(0.31,0.69,0.48,0.84);
			leg2->AddEntry(fr_bothcount_etaData,"Data");
			leg2->AddEntry(fr_bothcount_eta, "Drell-Yan ");
			leg2->Draw();
			
			if(RunYear==2016 and preVFP == 1)       CMS_lumi( canEta, 1, 1, 11 );
	                else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canEta, 2, 1, 11 );
        	        else if(RunYear==2017)                  CMS_lumi( canEta, 3, 1, 11 );
                	else if(RunYear==2018)                  CMS_lumi( canEta, 4, 1, 11 );
			canEta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_data_vs_DY_%d%s.pdf",RunYear,whichVFP.c_str()));


			TCanvas *canVtx = new TCanvas("canVtx","",1200,1000);
			canVtx->SetRightMargin(0.1);
			canVtx->SetLeftMargin(0.13);
			canVtx->cd();

                       	fr_bothcount_vtx->SetMinimum(0.0);
                       	fr_bothcount_vtx->SetMaximum(0.1);
                        fr_bothcount_vtx->GetXaxis()->SetTitle("nVtx ");
                        fr_bothcount_vtx->GetYaxis()->SetTitle("fake rate");
			fr_bothcount_vtx->GetXaxis()->SetTitleSize(.05);
			fr_bothcount_vtx->GetYaxis()->SetTitleSize(.05);
                        fr_bothcount_vtx->SetTitle("");
                      
                        fr_bothcount_vtx->SetMarkerColor(kRed);
                        fr_bothcount_vtx->SetLineColor(kRed);
                        fr_bothcount_vtxData->SetMarkerColor(kBlack);
                        fr_bothcount_vtxData->SetLineColor(kBlack);
                        
			fr_bothcount_vtx->SetMarkerStyle(kCircle);
			fr_bothcount_vtxData->SetMarkerStyle(kStar);
			
			fr_bothcount_vtx->SetMarkerSize(2);
			fr_bothcount_vtxData->SetMarkerSize(3);
                        
			fr_bothcount_vtxData->SetMaximum(0.1);
			fr_bothcount_vtx->Draw("Ap");
			fr_bothcount_vtxData->Draw("p");

			TLegend *legVtx = new TLegend(0.31,0.69,0.48,0.84);
			legVtx->AddEntry(fr_bothcount_vtxData,"Data");
			legVtx->AddEntry(fr_bothcount_vtx, "Drell-Yan ");
			legVtx->Draw();
			if(RunYear==2016 and preVFP == 1)       CMS_lumi( canVtx, 1, 1, 11 );
	                else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canVtx, 2, 1, 11 );
        	        else if(RunYear==2017)                  CMS_lumi( canVtx, 3, 1, 11 );
                	else if(RunYear==2018)                  CMS_lumi( canVtx, 4, 1, 11 );

			canVtx->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_nVtx_data_vs_DY_%d%s.pdf",RunYear,whichVFP.c_str()));
}
