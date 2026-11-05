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

#define MAXPT 200
#define MAXVTX 46

int MINPT = 30;
int MINVTX = 8;
bool doEB = true;
const char*processName ="Data"; bool doDrellYan = false;
//const char*processName ="DY"; bool doDrellYan = true;

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

bool EstimatePtFitStartValues(TGraphErrors *g, double &slope0, double &const0, double &index0)
{
	int n = g->GetN();
	if(n < 2) return false;

	double bestChi2 = 1e300;
	bool found = false;
	for(double idxTry = -0.9; idxTry <= -0.1 + 1e-9; idxTry += 0.02){
		double sumX=0, sumY=0, sumXY=0, sumXX=0; int npts=0;
		for(int i=0;i<n;i++){
			double x,y; g->GetPoint(i,x,y);
			if(y <= 0) continue;
			double arg = pow(y, 1.0/idxTry);
			if(!TMath::Finite(arg)) continue;
			sumX+=x; sumY+=arg; sumXY+=x*arg; sumXX+=x*x; npts++;
		}
		if(npts < 2) continue;
		double denom = npts*sumXX - sumX*sumX;
		if(fabs(denom) < 1e-9) continue;
		double s = (npts*sumXY - sumX*sumY)/denom;
		double c = (sumY - s*sumX)/npts;
		if(s <= 0) continue; // need arg increasing with pt for a falling fake rate

		double chi2 = 0; bool ok = true;
		for(int i=0;i<n;i++){
			double x,y; g->GetPoint(i,x,y);
			double arg = s*x + c;
			if(arg <= 0){ ok=false; break; }
			double model = pow(arg, idxTry);
			chi2 += (model-y)*(model-y);
		}
		if(!ok) continue;

		if(chi2 < bestChi2){
			bestChi2 = chi2;
			slope0 = s; const0 = c; index0 = idxTry;
			found = true;
		}
	}
	return found;
}

bool FitPtDependenceFromEstimate(TGraphErrors *g, TF1 *f, double slope0, double const0, double index0, TFitResultPtr &result)
{
	double slopeLo = TMath::Max(1.0, fabs(slope0)*0.02);
	double slopeHi = fabs(slope0)*50.0 + 10.0;
	f->SetParameter(0, slope0);
	f->SetParLimits(0, slopeLo, slopeHi);

	double constSpread = 5.0*fabs(const0) + 5000.0;
	f->SetParameter(1, const0);
	f->SetParLimits(1, const0 - constSpread, const0 + constSpread);

	double indexLo = TMath::Max(-0.95, index0 - 0.5);
	double indexHi = TMath::Min(-0.02, index0 + 0.5);
	f->SetParameter(2, index0);
	f->SetParLimits(2, indexLo, indexHi);

	result = g->Fit(f->GetName(), "R S Q W");
	return (result.Get()!=nullptr) && (result->Status() == 0);
}

// Fallback for a pT slice whose nominal nVtx fit did not converge: a straight line p0 + p1*x fit on the
// correct points of g only. Points with a non-finite or non-positive value or error, or an error
// above 100%, are dropped first; if the line fit still fails, the point with the largest relative
// error is dropped and the fit repeated, as long as 2 points remain. line must be a "pol1" TF1.
// Returns a null TFitResultPtr if no fit could be made.
TFitResultPtr FitLineOnGoodPoints(TGraphErrors *g, TF1 *line, const std::string &label)
{
	std::vector<int> keep;
	for(int i=0;i<g->GetN();i++){
		double x,y; g->GetPoint(i,x,y);
		double ey = g->GetErrorY(i);
		if(TMath::Finite(y) && TMath::Finite(ey) && y > 0 && ey > 0 && ey < y) keep.push_back(i);
		else std::cout << label << ": dropping point x=" << x << " (" << y << " +- " << ey << ")" << std::endl;
	}
	while(keep.size() >= 2){
		TGraphErrors gr(keep.size());
		for(size_t j=0;j<keep.size();j++){
			double x,y; g->GetPoint(keep[j],x,y);
			gr.SetPoint(j, x, y);
			gr.SetPointError(j, g->GetErrorX(keep[j]), g->GetErrorY(keep[j]));
		}
		TFitResultPtr r = gr.Fit(line, "S Q N");
		if(r.Get() && r->Status() == 0){
			std::cout << label << ": straight-line fit on " << keep.size() << " of " << g->GetN()
			          << " points: p0=" << line->GetParameter(0) << " p1=" << line->GetParameter(1) << std::endl;
			return r;
		}
		size_t worst = 0;
		for(size_t j=1;j<keep.size();j++)
			if(g->GetErrorY(keep[j])/g->GetY()[keep[j]] > g->GetErrorY(keep[worst])/g->GetY()[keep[worst]]) worst = j;
		if(keep.size() == 2) break;
		std::cout << label << ": straight-line fit failed, dropping point x=" << g->GetX()[keep[worst]] << std::endl;
		keep.erase(keep.begin() + worst);
	}
	std::cout << label << ": no straight-line fit possible" << std::endl;
	return TFitResultPtr();
}

TGraph* ReadFakeRatePointsOnly(const std::string &path, double xAxisMax)
{
	std::ifstream fin(path);
	if(!fin.is_open()){
		std::cout<<"Cannot open "<<path<<std::endl;
		return nullptr;
	}

	std::string bintype, numtype;
	float lowcut, signal1, error1, fitmean, fitrms;

	std::string line;
	int nLines=0;
	while(std::getline(fin,line)) nLines++;
	fin.clear(); fin.seekg(0, ios::beg);
	int nBins = nLines/2;

	std::vector<double> x, den, num;
	for(int i=0;i<nBins;i++){
		fin >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		x.push_back(lowcut);
		den.push_back(fitmean);
	}
	for(int i=0;i<nBins;i++){
		fin >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		num.push_back(fitmean);
	}

	TGraph *g = new TGraph(nBins);
	for(int i=0;i<nBins;i++){
		double fr = num[i]/den[i];
		double xvalue = (i+1 < nBins) ? (x[i]+x[i+1])/2.0 : (xAxisMax + x[i])/2.0;
		g->SetPoint(i, xvalue, fr);
	}
	return g;
}

TGraphErrors* ReadFakeRateGraphWithErrors(const std::string &path, double xAxisMax)
{
	std::ifstream fin(path);
	if(!fin.is_open()){
		std::cout<<"Cannot open "<<path<<std::endl;
		return nullptr;
	}

	std::string bintype, numtype;
	float lowcut, signal1, error1, fitmean, fitrms;

	std::string line;
	int nLines=0;
	while(std::getline(fin,line)) nLines++;
	fin.clear(); fin.seekg(0, ios::beg);
	int nBins = nLines/2;

	std::vector<double> x, den, denErr, num, numErr;
	for(int i=0;i<nBins;i++){
		fin >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		x.push_back(lowcut);
		den.push_back(fitmean);
		denErr.push_back(fitrms);
	}
	for(int i=0;i<nBins;i++){
		fin >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		num.push_back(fitmean);
		numErr.push_back(fitrms);
	}

	TGraphErrors *g = new TGraphErrors(nBins);
	for(int i=0;i<nBins;i++){
		double fr = num[i]/den[i];
		double er = sqrt( fr*fr*denErr[i]*denErr[i]/(den[i]*den[i]) + numErr[i]*numErr[i]/(den[i]*den[i]) );
		double xvalue = (i+1 < nBins) ? (x[i]+x[i+1])/2.0 : (xAxisMax + x[i])/2.0;
		double xerror = (i+1 < nBins) ? (x[i+1]-x[i])/2.0 : (xAxisMax - x[i])/2.0;
		g->SetPoint(i, xvalue, fr);
		g->SetPointError(i, xerror, er);
	}
	return g;
}

// Same as ReadFakeRateGraphWithErrors, but folds in the fit-function systematic the
// way the nominal (non-sliced) vtx/pt dependence does: the central value and its
// fitrms come from the Bw-ker file, and a systematic term is added in quadrature,
// taken as the larger of |DY-ker - Bw-ker| and |Bw-expo - Bw-ker| in each bin.
TGraphErrors* ReadFakeRateGraphWithSystematics(const std::string &path_ker,
                                                const std::string &path_DYker,
                                                const std::string &path_expo,
                                                double xAxisMax)
{
	std::ifstream fin_ker(path_ker);
	std::ifstream fin_DYker(path_DYker);
	std::ifstream fin_expo(path_expo);
	if(!fin_ker.is_open()){   std::cout<<"Cannot open "<<path_ker<<std::endl;   return nullptr; }
	if(!fin_DYker.is_open()){ std::cout<<"Cannot open "<<path_DYker<<std::endl; return nullptr; }
	if(!fin_expo.is_open()){  std::cout<<"Cannot open "<<path_expo<<std::endl;  return nullptr; }

	std::string bintype, numtype;
	float lowcut, signal1, error1, fitmean, fitrms, DYmean, DYrms, Polmean, Polrms;

	std::string line;
	int nLines=0;
	while(std::getline(fin_ker,line)) nLines++;
	fin_ker.clear(); fin_ker.seekg(0, ios::beg);
	int nBins = nLines/2;

	std::vector<double> x, den, denErr, num, numErr;
	for(int i=0;i<nBins;i++){
		fin_ker   >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		fin_DYker >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean  >> DYrms;
		fin_expo  >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
		double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean) ? fabs(DYmean - fitmean) : fabs(Polmean - fitmean);
		x.push_back(lowcut);
		den.push_back(fitmean);
		denErr.push_back(sqrt(fitrms*fitrms + sysdiff*sysdiff));
	}
	for(int i=0;i<nBins;i++){
		fin_ker   >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
		fin_DYker >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean  >> DYrms;
		fin_expo  >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
		double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean) ? fabs(DYmean - fitmean) : fabs(Polmean - fitmean);
		num.push_back(fitmean);
		numErr.push_back(sqrt(fitrms*fitrms + sysdiff*sysdiff));
	}

	TGraphErrors *g = new TGraphErrors(nBins);
	for(int i=0;i<nBins;i++){
		double fr = num[i]/den[i];
		double er = sqrt( fr*fr*denErr[i]*denErr[i]/(den[i]*den[i]) + numErr[i]*numErr[i]/(den[i]*den[i]) );
		double xvalue = (i+1 < nBins) ? (x[i]+x[i+1])/2.0 : (xAxisMax + x[i])/2.0;
		double xerror = (i+1 < nBins) ? (x[i+1]-x[i])/2.0 : (xAxisMax - x[i])/2.0;
		g->SetPoint(i, xvalue, fr);
		g->SetPointError(i, xerror, er);
	}
	return g;
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

void fitFakeFunc_pT_nVtx(int RunYear, bool preVFP, bool useDY = false){//main 
	// useDY: fit the DY MC bins (DYResult*) instead of DATA, outputs go to DrellYanResult*
	doDrellYan = useDY;
	processName = useDY ? "DY" : "Data";
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
	
	if(RunYear==2017)
                MINVTX = 18; 

	std::string isEB;
	if(doEB == true) isEB = "";
	else if(doEB == false) isEB = "_EE";
	
	std::ofstream resultfile, fakeRateFile, etafile;
	std::string fakeRateBaseName;

	if (doEB == true) {
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY.txt", RunYear, whichVFP.c_str()));
        		fakeRateBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPtVtx-EB", RunYear, whichVFP.c_str());
        		fakeRateFile.open(fakeRateBaseName + ".txt");
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_fitFakeFunc_Data.txt", RunYear, whichVFP.c_str()));
        		fakeRateBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPtVtx-EB", RunYear, whichVFP.c_str());
        		fakeRateFile.open(fakeRateBaseName + ".txt");
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_eta_dependence_Data.txt", RunYear, whichVFP.c_str()));
    		}} 
	else {  
    		if (doDrellYan == true) {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_fitFakeFunc_DY_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPtVtx-EE", RunYear, whichVFP.c_str());
        		fakeRateFile.open(fakeRateBaseName + ".txt");
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/result_eta_dependence_DY_EE.txt", RunYear, whichVFP.c_str()));
    		} else {
        		resultfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_fitFakeFunc_Data_EE.txt", RunYear, whichVFP.c_str()));
        		fakeRateBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPtVtx-EE", RunYear, whichVFP.c_str());
        		fakeRateFile.open(fakeRateBaseName + ".txt");
        		etafile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/result_eta_dependence_Data_EE.txt", RunYear, whichVFP.c_str()));
		}}
	
	std::ifstream Pt_file, Pt_DYfile, Pt_Polfile, Eta_file, Eta_DYfile, Eta_Polfile, Vtx_file, Vtx_DYfile, Vtx_Polfile;
	
	if(doDrellYan==true) {
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
     		Vtx_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-DY-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-Bw-expo-vtx-60-120.txt",RunYear,whichVFP.c_str())); 
		
	}
	else{
		Pt_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));
        	Pt_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",RunYear,whichVFP.c_str()));

		Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
        	Eta_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-eta-60-120.txt",RunYear,whichVFP.c_str()));

		Vtx_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_DYfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-vtx-60-120.txt",RunYear,whichVFP.c_str()));
       		Vtx_Polfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-vtx-60-120.txt",RunYear,whichVFP.c_str())); 
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
	TGraphErrors *fr_pt_ratio = new TGraphErrors(nPtBins);
	TGraphErrors *fr_pt_ratioError = new TGraphErrors(nPtBins);
	TGraphErrors *fr_vtx_ratio = new TGraphErrors(nVtxBins);
	TGraphErrors *fr_vtx_ratioError = new TGraphErrors(nVtxBins);

	int graphPtBins = (int)((MAXPT - MINPT)/0.5);
	int graphVtxBins = (int)((MAXVTX - MINVTX)/0.5);
	TGraphErrors *fr_pt_sigmaband = new TGraphErrors(graphPtBins);
	TGraphErrors *fr_vtx_sigmaband = new TGraphErrors(graphVtxBins);
	TGraphErrors *fr_eta_sigmaband = new TGraphErrors(nEtaBins);

//************************* Fill and fit the pt dependence. Calculate ratio ********************* 
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

	// Not converged (e.g. constant stuck at its limit): refit from data-driven start values, as the nVtx slices do.
	// Fit with the point errors (no "W"), since the toys below sample from this covariance.
	if(result_fitpt->Status() != 0){
		double slope0=0, const0=0, index0=0;
		if(EstimatePtFitStartValues(fr_bothcount_pt, slope0, const0, index0)){
			std::cout << "Inclusive pT fit status=" << result_fitpt->Status() << ", refitting from start values slope="
			          << slope0 << " constant=" << const0 << " index=" << index0 << std::endl;
			f1->SetParameter(0, slope0);
			f1->SetParLimits(0, TMath::Max(1.0, fabs(slope0)*0.02), fabs(slope0)*50.0 + 10.0);
			double constSpread = 5.0*fabs(const0) + 5000.0;
			f1->SetParameter(1, const0);
			f1->SetParLimits(1, const0 - constSpread, const0 + constSpread);
			f1->SetParameter(2, index0);
			f1->SetParLimits(2, TMath::Max(-0.95, index0 - 0.5), TMath::Min(-0.02, index0 + 0.5));
			result_fitpt = fr_bothcount_pt->Fit("f1","R S");
			resultfile << "pt fit refit from data-driven start values: status=" << result_fitpt->Status() << std::endl;
		}
	}

	int status = result_fitpt->Status();        // 0 means OK
	int covQual = result_fitpt->CovMatrixStatus();  // 3 is the best

	std::cout << "Fit status = " << status << std::endl;
	std::cout << "CovQual = " << covQual << std::endl;

 	fit_fakerate_pt = fr_bothcount_pt->GetFunction("f1");
	resultfile <<"pt chiSquare() : "<<fit_fakerate_pt->GetChisquare()/fit_fakerate_pt->GetNDF()<<endl;

	double avg_fr_pt = fit_fakerate_pt->Integral(MINPT, MAXPT) / (MAXPT - MINPT);
	std::cout << "avg_fr_pt (mean fake rate over pT range) = " << avg_fr_pt << std::endl;
	resultfile << "avg_fr_pt = " << avg_fr_pt << std::endl;

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

	double avg_fr_eta = 0;
	for(unsigned i(0); i<nEtaBins; i++) avg_fr_eta += eta_num[i]/eta_den[i];
	avg_fr_eta /= (nEtaBins > 0 ? nEtaBins : 1);
	std::cout << "avg_fr_eta (mean fake rate over eta bins) = " << avg_fr_eta << std::endl;
	resultfile << "avg_fr_eta = " << avg_fr_eta << std::endl;
 
 	if(Vtx_file.is_open() && Vtx_DYfile.is_open() && Vtx_Polfile.is_open()){
 		for(unsigned i(0); i<nVtxBins; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_den[i] = fitmean;
			VtxBins[i] = lowcut;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			vtx_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 		for(unsigned i(0); i<nVtxBins; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_num[i] = fitmean;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			vtx_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 	}
 	for(unsigned i(0); i<nVtxBins; i++){
 		double fakerate = vtx_num[i]/vtx_den[i];
 		double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
		double xvalue = (i+1 < nVtxBins)? (VtxBins[i]+VtxBins[i+1])/2.0 : (MAXVTX + VtxBins[i])/2.0;
    		double xerror = (i+1 < nVtxBins)? (VtxBins[i+1]-VtxBins[i])/2.0 : (MAXVTX - VtxBins[i])/2.0;
 		fr_bothcount_vtx->SetPoint(i,xvalue, fakerate );
 		fr_bothcount_vtx->SetPointError(i, xerror, error);
 	}

 	TF1 *fit_fakerate_vtx = new TF1("vtxfake","pol1",MINVTX,46);
 	TFitResultPtr result_fitvtx;
 	//fit_fakerate_vtx->SetParameter(0,0.02);
 	//fit_fakerate_vtx->SetParameter(1,1);
 	result_fitvtx = fr_bothcount_vtx->Fit("vtxfake","S");
	resultfile <<"vtx chiSquare() : "<<fit_fakerate_vtx->GetChisquare()/fit_fakerate_vtx->GetNDF()<<endl;

	// Mean of the nVtx fit over its fit range -- turns fit_fakerate_vtx into a
	// pure shape correction (nVtx-dependent multiplicative factor centered on
	// 1) instead of a second independent normalization on top of the pT fit.
	double avg_fr_vtx = fit_fakerate_vtx->Integral(MINVTX, 46) / (46 - MINVTX);
	std::cout << "avg_fr_vtx (mean fake rate over nVtx range) = " << avg_fr_vtx << std::endl;
	resultfile << "avg_fr_vtx = " << avg_fr_vtx << std::endl;

 	for(unsigned i(0); i<nVtxBins; i++){
 		double fakerate = vtx_num[i]/vtx_den[i];
 		double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
		double xvalue = (i+1 < nVtxBins)? (VtxBins[i]+VtxBins[i+1])/2.0 : (MAXVTX + VtxBins[i])/2.0;
    		double xerror = (i+1 < nVtxBins)? (VtxBins[i+1]-VtxBins[i])/2.0 : (MAXVTX - VtxBins[i])/2.0;
 		fr_vtx_ratio->SetPoint(i,xvalue, fit_fakerate_vtx->Eval(xvalue)/fakerate);
 		fr_vtx_ratio->SetPointError(i, xerror, error/fit_fakerate_vtx->Eval(xvalue));
 	}
 
	std::vector<std::string> pTSliceDirs   = {"pT30-50","pT50-80","pT-gt80"};
	std::vector<TGraphErrors*> vtxSliceGraphs(pTSliceDirs.size(), nullptr);
	std::vector<TF1*> fit_fakerate_vtx_pTSlice(pTSliceDirs.size(), nullptr);
	std::vector<TFitResultPtr> result_vtx_pTSlice(pTSliceDirs.size()); // kept for the VtxByPtSlice toys
	const std::string sampleType = doDrellYan ? "DY" : "DATA";
	for(size_t is=0; is<pTSliceDirs.size(); is++){
		std::string sliceDir = Form("/eos/uscms/store/user/tmishra/elefakepho/%sResult%d%s_%s",sampleType.c_str(),
		                             RunYear, whichVFP.c_str(), pTSliceDirs[is].c_str());
		std::string subfile_ker   = sliceDir + Form("/EleFakeRate-%s-%s-Bw-ker-vtx-60-120.txt", sampleType.c_str(), pTSliceDirs[is].c_str());
		std::string subfile_DYker = sliceDir + Form("/EleFakeRate-%s-%s-DY-ker-vtx-60-120.txt", sampleType.c_str(), pTSliceDirs[is].c_str());
		std::string subfile_expo  = sliceDir + Form("/EleFakeRate-%s-%s-Bw-expo-vtx-60-120.txt", sampleType.c_str(), pTSliceDirs[is].c_str());
		vtxSliceGraphs[is] = ReadFakeRateGraphWithSystematics(subfile_ker, subfile_DYker, subfile_expo, MAXVTX);
		if(!vtxSliceGraphs[is]) continue;

		std::string fname = Form("vtxfake_slice_%s", pTSliceDirs[is].c_str());
		TF1 *fSlice = new TF1(fname.c_str(), "pol1", MINVTX, 46);
		TFitResultPtr resultSlice = vtxSliceGraphs[is]->Fit(fname.c_str(), "S");
		std::cout << "pT slice " << pTSliceDirs[is] << " (nVtx dependence): fit status=" << resultSlice->Status() << std::endl;
		resultfile << "pT slice " << pTSliceDirs[is] << " (nVtx dependence): fit status=" << resultSlice->Status()
		           << " chi2/ndf=" << (fSlice->GetNDF()>0 ? fSlice->GetChisquare()/fSlice->GetNDF() : -1) << std::endl;
		if(resultSlice->Status() != 0){
			resultSlice = FitLineOnGoodPoints(vtxSliceGraphs[is], fSlice, "pT slice " + pTSliceDirs[is]);
			if(resultSlice.Get()) resultfile << "pT slice " << pTSliceDirs[is] << ": refit on reduced points, p0="
			                                  << fSlice->GetParameter(0) << " p1=" << fSlice->GetParameter(1) << std::endl;
		}
		if(resultSlice.Get() && resultSlice->Status() == 0){ fit_fakerate_vtx_pTSlice[is] = fSlice; result_vtx_pTSlice[is] = resultSlice; }
	}

	std::vector<std::string> nVtxSliceDirs   = {"nVtx8-20","nVtx20-30","nVtx-gt30"};
	std::vector<double> nVtxSliceUpperEdge   = {20.0, 30.0, 1e9}; // event falls in first slice whose nVertex < edge
	std::vector<TGraphErrors*> ptSliceGraphs(nVtxSliceDirs.size(), nullptr);
	std::vector<TF1*> fit_fakerate_pt_nVtxSlice(nVtxSliceDirs.size(), nullptr);

	for(size_t is=0; is<nVtxSliceDirs.size(); is++){
		std::string sliceDir = Form("/eos/uscms/store/user/tmishra/elefakepho/%sResult%d%s_%s", sampleType.c_str(),
		                             RunYear, whichVFP.c_str(), nVtxSliceDirs[is].c_str());
		std::string subfile_ker   = sliceDir + Form("/EleFakeRate-%s-%s-Bw-ker-pt-60-120.txt", sampleType.c_str(), nVtxSliceDirs[is].c_str());
		std::string subfile_DYker = sliceDir + Form("/EleFakeRate-%s-%s-DY-ker-pt-60-120.txt", sampleType.c_str(), nVtxSliceDirs[is].c_str());
		std::string subfile_expo  = sliceDir + Form("/EleFakeRate-%s-%s-Bw-expo-pt-60-120.txt", sampleType.c_str(), nVtxSliceDirs[is].c_str());
		ptSliceGraphs[is] = ReadFakeRateGraphWithSystematics(subfile_ker, subfile_DYker, subfile_expo, MAXPT);
		if(!ptSliceGraphs[is]) continue; // stays nullptr -> falls back to nominal fit_fakerate_pt at use site

		std::string fname = Form("f1_slice_%s", nVtxSliceDirs[is].c_str());
		TF1 *fSlice = new TF1(fname.c_str(), fakerate_ptDependence, MINPT, MAXPT, 3);
		fSlice->SetParNames("slope","constant","index");

		double slope0=0, const0=0, index0=0;
		TFitResultPtr resultSlice;
		bool converged = false;
		if(EstimatePtFitStartValues(ptSliceGraphs[is], slope0, const0, index0))
			converged = FitPtDependenceFromEstimate(ptSliceGraphs[is], fSlice, slope0, const0, index0, resultSlice);

		if(!converged){
			std::cout << "nVtx slice " << nVtxSliceDirs[is] << ": data-driven start failed/did not converge, retrying with legacy fixed starting values" << std::endl;
			fSlice->SetParameter(0,1.0e+04);
			fSlice->SetParLimits(0, 1.001e+02, 1.25601e+04);
			fSlice->SetParameter(2,-0.4);
			fSlice->SetParLimits(2, -0.4, -0.1);
			if(!doDrellYan){
				if(RunYear==2016){ fSlice->SetParameter(1, -10000); fSlice->SetParLimits(1, -15000.0, -5000); }
				if(RunYear==2017){ fSlice->SetParameter(1, -20000); fSlice->SetParLimits(1, -39070.0, -10000); }
				if(RunYear==2018){ fSlice->SetParameter(1, -30000); fSlice->SetParLimits(1, -38685.0, 0); }
			} else {
				fSlice->SetParameter(1, -38000);
				fSlice->SetParLimits(1, -38100.0, -10000);
			}
			resultSlice = ptSliceGraphs[is]->Fit(fname.c_str(), "R S W");
			converged = (resultSlice->Status() == 0);
		}

		std::cout << "nVtx slice " << nVtxSliceDirs[is] << ": fit status=" << resultSlice->Status()
		          << " covQual=" << resultSlice->CovMatrixStatus() << std::endl;
		resultfile << "nVtx slice " << nVtxSliceDirs[is] << ": fit status=" << resultSlice->Status()
		           << " covQual=" << resultSlice->CovMatrixStatus()
		           << " chi2/ndf=" << (fSlice->GetNDF()>0 ? fSlice->GetChisquare()/fSlice->GetNDF() : -1) << std::endl;

		if(converged) fit_fakerate_pt_nVtxSlice[is] = fSlice;
		else std::cout << "  -> fit did not converge cleanly, falling back to nominal fit_fakerate_pt for this slice" << std::endl;
	}

	auto SelectPtFitForNVtx = [&](int nVertexVal) -> TF1* {
		for(size_t is=0; is<nVtxSliceUpperEdge.size(); is++){
			if(nVertexVal < nVtxSliceUpperEdge[is]){
				return fit_fakerate_pt_nVtxSlice[is] ? fit_fakerate_pt_nVtxSlice[is] : fit_fakerate_pt;
			}
		}
		return fit_fakerate_pt;
	};

	
	TH1D* invmass_den = new TH1D("invmass_den", "invmass_den",100,40,140);
	TH1D* invmass_num = new TH1D("invmass_num", "invmass_num",100,40,140);
	TH1D* h_invmass_bg  = new TH1D("invmass_bg",  "invmass_bg", 100,40,140);
	TH1D* invmass_prednum = new TH1D("invmass_prednum", "invmass_prednum",100,40,140);

	std::vector<double> ptSliceLo = {30.0, 50.0, 80.0};
	std::vector<double> ptSliceHi = {50.0, 80.0, 1e9};
	std::vector<TH1D*> invmass_den_slice, invmass_num_slice, h_invmass_bg_slice, invmass_prednum_slice;
	// Fix #2: a SECOND, independent predicted-yield histogram per pT slice,
	// built with the vtx-linear shape for that slice (instead of the pT-shape
	// weighting used for invmass_prednum_slice above). This is what the
	// scalefactor written alongside vtxconst/vtxslope will actually be
	// calibrated against, instead of reusing a scalefactor that was
	// calibrated with a different (pT-shape) prediction.
	std::vector<TH1D*> invmass_prednum_slice_vtxshape;
	for(size_t is=0; is<pTSliceDirs.size(); is++){
		invmass_den_slice.push_back(new TH1D(Form("invmass_den_%s",pTSliceDirs[is].c_str()), "invmass_den_slice",100,40,140));
		invmass_num_slice.push_back(new TH1D(Form("invmass_num_%s",pTSliceDirs[is].c_str()), "invmass_num_slice",100,40,140));
		h_invmass_bg_slice.push_back(new TH1D(Form("invmass_bg_%s",pTSliceDirs[is].c_str()), "invmass_bg_slice",100,40,140));
		invmass_prednum_slice.push_back(new TH1D(Form("invmass_prednum_%s",pTSliceDirs[is].c_str()), "invmass_prednum_slice",100,40,140));
		invmass_prednum_slice_vtxshape.push_back(new TH1D(Form("invmass_prednum_vtxshape_%s",pTSliceDirs[is].c_str()), "invmass_prednum_slice_vtxshape",100,40,140));
	}
	auto SelectPtSliceIndex = [&](double pt) -> int {
		for(size_t is=0; is<ptSliceLo.size(); is++)
			if(pt >= ptSliceLo[is] && pt < ptSliceHi[is]) return (int)is;
		return -1;
	};

	// Fix #1: average of the INCLUSIVE pT shape over each pT slice's own
	// range -- used to build a residual pT correction (residualPt(pt) =
	// fit_fakerate_pt(pt)/avg_fr_pt_pTSlice[is]) so the vtx-shape-by-pT-slice
	// weighting doesn't silently flatten out the real pT dependence *within*
	// a slice. Same pattern as weight_eta/avg_fr_eta already used for eta.
	std::vector<double> avg_fr_pt_pTSlice(pTSliceDirs.size(), 1.0);
	for(size_t is=0; is<pTSliceDirs.size(); is++){
		double lo = ptSliceLo[is];
		double hi = std::min(ptSliceHi[is], (double)MAXPT);
		double width = hi - lo;
		avg_fr_pt_pTSlice[is] = (width > 0) ? fit_fakerate_pt->Integral(lo, hi)/width : avg_fr_pt;
	}

	// Selects the vtx-linear fit (fit_fakerate_vtx_pTSlice, falling back to
	// the nominal fit_fakerate_vtx) for a given photon pT, mirroring
	// SelectPtFitForNVtx above.
	auto SelectVtxFitForPt = [&](double ptVal) -> TF1* {
		int idx = SelectPtSliceIndex(ptVal);
		if(idx < 0) return fit_fakerate_vtx;
		return fit_fakerate_vtx_pTSlice[idx] ? fit_fakerate_vtx_pTSlice[idx] : fit_fakerate_vtx;
	};

	std::vector<TH1D*> invmass_den_nVtxSlice, invmass_num_nVtxSlice, h_invmass_bg_nVtxSlice, invmass_prednum_nVtxSlice;
	for(size_t is=0; is<nVtxSliceDirs.size(); is++){
		invmass_den_nVtxSlice.push_back(new TH1D(Form("invmass_den_%s",nVtxSliceDirs[is].c_str()), "invmass_den_nVtxSlice",100,40,140));
		invmass_num_nVtxSlice.push_back(new TH1D(Form("invmass_num_%s",nVtxSliceDirs[is].c_str()), "invmass_num_nVtxSlice",100,40,140));
		h_invmass_bg_nVtxSlice.push_back(new TH1D(Form("invmass_bg_%s",nVtxSliceDirs[is].c_str()), "invmass_bg_nVtxSlice",100,40,140));
		invmass_prednum_nVtxSlice.push_back(new TH1D(Form("invmass_prednum_%s",nVtxSliceDirs[is].c_str()), "invmass_prednum_nVtxSlice",100,40,140));
	}
	auto SelectNVtxSliceIndex = [&](int nVertexVal) -> int {
		if(nVertexVal < MINVTX) return -1;
		for(size_t is=0; is<nVtxSliceUpperEdge.size(); is++)
			if(nVertexVal < nVtxSliceUpperEdge[is]) return (int)is;
		return -1;
	};
	
	TChain *etree = new TChain("FakeRateTree");
	etree->Add(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_%sTnP_dR05_%d%s.root",processName,RunYear,whichVFP.c_str()));

	float invmass=0; 
	float tagPt=0; 
	float probePt=0; 
	float probeEta=0;
	bool  vetovalue=0;
	bool  FSRveto = 0;
	int   nVertex=0;
	etree->SetBranchAddress("invmass",   &invmass); 
	etree->SetBranchAddress("tagPt",     &tagPt);
	etree->SetBranchAddress("probePt",   &probePt);
	etree->SetBranchAddress("probeEta",  &probeEta);
	etree->SetBranchAddress("vetovalue", &vetovalue);
	etree->SetBranchAddress("nVertex",   &nVertex);
	etree->SetBranchAddress("FSRveto",   &FSRveto);
	
	std::vector<float> etreeEt;
	std::vector<float> etreeVtx;
	std::vector<float> etreeEta;
	std::vector<float> etreeInvmass;
	etreeEt.clear();
	etreeVtx.clear(); 
	etreeEta.clear();
	etreeInvmass.clear();
	for(unsigned iEvt(0); iEvt < etree->GetEntries(); iEvt++){
		etree->GetEntry(iEvt);
		     	
		if(probePt < 35)continue;
		if(doEB && fabs(probeEta) > 1.4442)continue;
		else if(!doEB && (fabs(probeEta) < 1.56 || fabs(probeEta) > 2.1))continue;

		if(vetovalue== false)invmass_den->Fill(invmass);
		else if(vetovalue== true && FSRveto == true)invmass_num->Fill(invmass);

		int ptSliceIdx = SelectPtSliceIndex(probePt);
		if(ptSliceIdx >= 0){
			if(vetovalue== false)invmass_den_slice[ptSliceIdx]->Fill(invmass);
			else if(vetovalue== true && FSRveto == true)invmass_num_slice[ptSliceIdx]->Fill(invmass);
		}

		int nVtxSliceIdx = SelectNVtxSliceIndex(nVertex);
		if(nVtxSliceIdx >= 0){
			if(vetovalue== false)invmass_den_nVtxSlice[nVtxSliceIdx]->Fill(invmass);
			else if(vetovalue== true && FSRveto == true)invmass_num_nVtxSlice[nVtxSliceIdx]->Fill(invmass);
		}
	
		float weight_eta(0);
 		for(unsigned i(0); i<nEtaBins; i++){
			if(i+1 < nEtaBins){
				if(fabs(probeEta) > EtaBins[i] && fabs(probeEta) < EtaBins[i+1])weight_eta = eta_num[i]/eta_den[i];
			}
			else{ 
				if(fabs(probeEta) > EtaBins[i])weight_eta = eta_num[i]/eta_den[i];
			}
		}
				
		TF1* ptFitForEvent = SelectPtFitForNVtx(nVertex);
		float w_ele = ptFitForEvent->Eval(probePt)
		            * (avg_fr_eta > 0 ? weight_eta / avg_fr_eta : 1.0);

		// Fix #1/#2: independent prediction for the VtxByPtSlice scheme --
		// vtx-linear shape for this event's pT slice, times a residual pT
		// correction (normalized to average 1 across the slice, same pattern
		// as weight_eta/avg_fr_eta) so the pT dependence within the slice
		// isn't silently dropped. This is what scalefactor_slice_vtxbased
		// below is calibrated against, instead of reusing w_ele's prediction.
		float w_ele_vtxshape = 0;
		if(ptSliceIdx >= 0){
			TF1* vtxFitForEvent = SelectVtxFitForPt(probePt);
			double residualPt = (avg_fr_pt_pTSlice[ptSliceIdx] > 0)
			                   ? fit_fakerate_pt->Eval(probePt)/avg_fr_pt_pTSlice[ptSliceIdx] : 1.0;
			w_ele_vtxshape = vtxFitForEvent->Eval(nVertex) * residualPt
			                * (avg_fr_eta > 0 ? weight_eta / avg_fr_eta : 1.0);
		}

		if(vetovalue==false){
			invmass_prednum->Fill(invmass, w_ele);
			if(ptSliceIdx >= 0){
				invmass_prednum_slice[ptSliceIdx]->Fill(invmass, w_ele);
				invmass_prednum_slice_vtxshape[ptSliceIdx]->Fill(invmass, w_ele_vtxshape);
			}
			if(nVtxSliceIdx >= 0) invmass_prednum_nVtxSlice[nVtxSliceIdx]->Fill(invmass, w_ele);
			etreeEt.push_back(probePt);
			etreeVtx.push_back(nVertex);
			etreeEta.push_back(fabs(probeEta));
			etreeInvmass.push_back(invmass);
		} 
	}
	//invmass_prednum->Sumw2();
 
	TChain *bgtree = new TChain("BGTree");
	bgtree->Add(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_bgtemplate_FullEcal_%d%s.root",RunYear,whichVFP.c_str()));
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
		if(doEB && fabs(probeEta_bg) > 1.4442)continue;
		else if(!doEB && (fabs(probeEta_bg) < 1.56 || fabs(probeEta_bg) > 2.1))continue;
	     	
		if(vetovalue_bg== true && FSRveto_bg==true){
			h_invmass_bg->Fill(invmass_bg);
			int ptSliceIdx_bg = SelectPtSliceIndex(probePt_bg);
			if(ptSliceIdx_bg >= 0) h_invmass_bg_slice[ptSliceIdx_bg]->Fill(invmass_bg);
			int nVtxSliceIdx_bg = SelectNVtxSliceIndex(nVertex_bg);
			if(nVtxSliceIdx_bg >= 0) h_invmass_bg_nVtxSlice[nVtxSliceIdx_bg]->Fill(invmass_bg);
		}
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
	cout<<"integral of hist is "<< h_DYinvmass->Integral()<<endl;
	//h_DYinvmass->Sumw2();

        
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
		c_fitMass->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/fit_totalNum_DrellYan_%d%s.pdf",RunYear,whichVFP.c_str()));
	else 
		c_fitMass->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/fit_totalNum_data_%d%s.pdf",RunYear,whichVFP.c_str()));
	mass_axis.setRange("signal",70,110);
	RooAbsReal* igx_sig = signalResDY.createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	std::cout << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " error = " << igx_sig->getVal()*(nSig.getError()) <<  std::endl;
 	std::cout << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl; 
 	std::cout << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
	resultfile << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " error = " << igx_sig->getVal()*(nSig.getError()) <<  std::endl;
 	resultfile << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl; 
 	resultfile << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 	fakeRateFile << "scalefactor " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;

	mean.setConstant(true);
	sigma.setConstant(true);

	std::vector<double> scalefactor_slice(pTSliceDirs.size(), -1);
	// Fix #2: scalefactor calibrated against the vtx-shape (residual-pT-
	// corrected) prediction -- this is the one that's actually self-
	// consistent with the vtxconst/vtxslope parameterization written into
	// the same output file, and should be used wherever that
	// parameterization is applied (EvalFakeRate_VtxByPtSlice downstream).
	std::vector<double> scalefactor_slice_vtxbased(pTSliceDirs.size(), -1);
	std::vector<double> numSignalSlice_v(pTSliceDirs.size(), 0), numSignalSliceErr_v(pTSliceDirs.size(), 0); // for the VtxByPtSlice toys
	for(size_t is=0; is<pTSliceDirs.size(); is++){
		double iniSigSlice = 0.2*invmass_num_slice[is]->Integral(1,100);
		int    lowBinSlice = invmass_num_slice[is]->FindBin(70.0);
		int    highBinSlice= invmass_num_slice[is]->FindBin(110.0);
		double iniBkgSlice = 40*(invmass_num_slice[is]->GetBinContent(lowBinSlice)+invmass_num_slice[is]->GetBinContent(highBinSlice));

		RooDataHist datahist_bg_slice(Form("bg_%s",pTSliceDirs[is].c_str()), "", mass_axis, h_invmass_bg_slice[is]);
		RooHistPdf  pdf_bg_slice(Form("pdf_bg_%s",pTSliceDirs[is].c_str()), "", mass_axis, datahist_bg_slice);

		RooRealVar nSig_slice(Form("nSig_%s",pTSliceDirs[is].c_str()), "", iniSigSlice, 0, invmass_num_slice[is]->GetEntries()*1.2+1);
		RooRealVar nBkg_slice(Form("nBkg_%s",pTSliceDirs[is].c_str()), "", iniBkgSlice, 0, invmass_num_slice[is]->GetEntries()+1);
		RooAddPdf  model_slice(Form("model_%s",pTSliceDirs[is].c_str()), "", RooArgList(pdf_bg_slice, signalResDY), RooArgList(nBkg_slice, nSig_slice));

		RooDataHist datahist_data_slice(Form("both_%s",pTSliceDirs[is].c_str()), "", mass_axis, invmass_num_slice[is]);
		model_slice.fitTo(datahist_data_slice);

		// Observed signal yield for this slice -- independent of which
		// prediction (pT-shape or vtx-shape) it gets compared against below.
		double numSignalSlice = igx_sig->getVal()*nSig_slice.getVal();
		numSignalSlice_v[is]    = numSignalSlice;
		numSignalSliceErr_v[is] = igx_sig->getVal()*nSig_slice.getError();

		double predictSlice = invmass_prednum_slice[is]->Integral(lowBinNumber,highBinNumber);
		if(predictSlice > 0){
			scalefactor_slice[is] = numSignalSlice/predictSlice;
		}

		double predictSlice_vtxbased = invmass_prednum_slice_vtxshape[is]->Integral(lowBinNumber,highBinNumber);
		if(predictSlice_vtxbased > 0){
			scalefactor_slice_vtxbased[is] = numSignalSlice/predictSlice_vtxbased;
		}

		std::cout << "pT slice " << pTSliceDirs[is] << ": num=" << numSignalSlice
		          << " predict(ptshape)=" << predictSlice << " scalefactor(ptshape)=" << scalefactor_slice[is]
		          << " predict(vtxshape)=" << predictSlice_vtxbased << " scalefactor(vtxshape)=" << scalefactor_slice_vtxbased[is] << std::endl;
		resultfile << "pT slice " << pTSliceDirs[is] << ": num=" << numSignalSlice
		           << " predict(ptshape)=" << predictSlice << " scalefactor(ptshape)=" << scalefactor_slice[is]
		           << " predict(vtxshape)=" << predictSlice_vtxbased << " scalefactor(vtxshape)=" << scalefactor_slice_vtxbased[is] << std::endl;
	}

	mean.setConstant(false);
	sigma.setConstant(false);

	std::vector<std::string> ptSliceFileTag = {"_pT30-50","_pT50-80","_pT-gt80"};
	for(size_t is=0; is<pTSliceDirs.size(); is++){
		std::ofstream fakeRateFileSlice;
		fakeRateFileSlice.open(fakeRateBaseName + ptSliceFileTag[is] + ".txt");
		// Fix #2: use the scalefactor calibrated against the vtx-shape
		// (residual-pT-corrected) prediction, so it's self-consistent with
		// the vtxconst/vtxslope parameterization written below.
		fakeRateFileSlice << "scalefactor " << scalefactor_slice_vtxbased[is] << std::endl;
		if(fit_fakerate_vtx_pTSlice[is]){
			fakeRateFileSlice << "vtxconst " << fit_fakerate_vtx_pTSlice[is]->GetParameter(0) << std::endl;
			fakeRateFileSlice << "vtxslope " << fit_fakerate_vtx_pTSlice[is]->GetParameter(1) << std::endl;
		} else {
			// fell back to the nominal (full-pT-range) nVtx fit -- see the
			// per-pT-slice nVtx fitting block earlier in this function.
			fakeRateFileSlice << "vtxconst " << fit_fakerate_vtx->GetParameter(0) << std::endl;
			fakeRateFileSlice << "vtxslope " << fit_fakerate_vtx->GetParameter(1) << std::endl;
		}
		// Fix #1: this slice's average of the inclusive pT shape -- the
		// downstream consumer divides the inclusive pT shape by this number
		// to build a residual pT correction (average 1 across the slice),
		// exactly matching what calibrated scalefactor_slice_vtxbased above,
		// instead of recomputing (and potentially not matching) it.
		fakeRateFileSlice << "ptshapeavg " << avg_fr_pt_pTSlice[is] << std::endl;
		fakeRateFileSlice.close();
	}

	mean.setConstant(true);
	sigma.setConstant(true);

	std::vector<double> scalefactor_nVtxSlice(nVtxSliceDirs.size(), -1);
	for(size_t is=0; is<nVtxSliceDirs.size(); is++){
		double iniSigSlice = 0.2*invmass_num_nVtxSlice[is]->Integral(1,100);
		int    lowBinSlice = invmass_num_nVtxSlice[is]->FindBin(70.0);
		int    highBinSlice= invmass_num_nVtxSlice[is]->FindBin(110.0);
		double iniBkgSlice = 40*(invmass_num_nVtxSlice[is]->GetBinContent(lowBinSlice)+invmass_num_nVtxSlice[is]->GetBinContent(highBinSlice));

		RooDataHist datahist_bg_nVtxSlice(Form("bg_%s",nVtxSliceDirs[is].c_str()), "", mass_axis, h_invmass_bg_nVtxSlice[is]);
		RooHistPdf  pdf_bg_nVtxSlice(Form("pdf_bg_%s",nVtxSliceDirs[is].c_str()), "", mass_axis, datahist_bg_nVtxSlice);

		RooRealVar nSig_nVtxSlice(Form("nSig_%s",nVtxSliceDirs[is].c_str()), "", iniSigSlice, 0, invmass_num_nVtxSlice[is]->GetEntries()*1.2+1);
		RooRealVar nBkg_nVtxSlice(Form("nBkg_%s",nVtxSliceDirs[is].c_str()), "", iniBkgSlice, 0, invmass_num_nVtxSlice[is]->GetEntries()+1);
		RooAddPdf  model_nVtxSlice(Form("model_%s",nVtxSliceDirs[is].c_str()), "", RooArgList(pdf_bg_nVtxSlice, signalResDY), RooArgList(nBkg_nVtxSlice, nSig_nVtxSlice));

		RooDataHist datahist_data_nVtxSlice(Form("both_%s",nVtxSliceDirs[is].c_str()), "", mass_axis, invmass_num_nVtxSlice[is]);
		model_nVtxSlice.fitTo(datahist_data_nVtxSlice);

		double predictSlice = invmass_prednum_nVtxSlice[is]->Integral(lowBinNumber,highBinNumber);
		if(predictSlice > 0){
			scalefactor_nVtxSlice[is] = igx_sig->getVal()*nSig_nVtxSlice.getVal()/predictSlice;
		}
		std::cout << "nVtx slice " << nVtxSliceDirs[is] << ": num=" << igx_sig->getVal()*nSig_nVtxSlice.getVal()
		          << " predict=" << predictSlice << " scalefactor=" << scalefactor_nVtxSlice[is] << std::endl;
		resultfile << "nVtx slice " << nVtxSliceDirs[is] << ": num=" << igx_sig->getVal()*nSig_nVtxSlice.getVal()
		           << " predict=" << predictSlice << " scalefactor=" << scalefactor_nVtxSlice[is] << std::endl;
	}

	mean.setConstant(false);
	sigma.setConstant(false);

	std::vector<std::string> nVtxSliceFileTag = {"_nVtx8-20","_nVtx20-30","_nVtx-gt30"};
	for(size_t is=0; is<nVtxSliceDirs.size(); is++){
		std::ofstream fakeRateFileNVtxSlice;
		fakeRateFileNVtxSlice.open(fakeRateBaseName + nVtxSliceFileTag[is] + ".txt");
		fakeRateFileNVtxSlice << "scalefactor " << scalefactor_nVtxSlice[is] << std::endl;
		if(fit_fakerate_pt_nVtxSlice[is]){
			fakeRateFileNVtxSlice << "ptslope "    << fit_fakerate_pt_nVtxSlice[is]->GetParameter(0) << std::endl;
			fakeRateFileNVtxSlice << "ptconstant " << fit_fakerate_pt_nVtxSlice[is]->GetParameter(1) << std::endl;
			fakeRateFileNVtxSlice << "ptindex "    << fit_fakerate_pt_nVtxSlice[is]->GetParameter(2) << std::endl;
		} else {
			// fell back to the nominal (nVtx-inclusive) pT fit -- see the
			// per-nVtx-slice pT fitting block earlier in this function.
			fakeRateFileNVtxSlice << "ptslope "    << fit_fakerate_pt->GetParameter(0) << std::endl;
			fakeRateFileNVtxSlice << "ptconstant " << fit_fakerate_pt->GetParameter(1) << std::endl;
			fakeRateFileNVtxSlice << "ptindex "    << fit_fakerate_pt->GetParameter(2) << std::endl;
		}
		fakeRateFileNVtxSlice.close();
	}

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
	if(doDrellYan)cancompare->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/compare_predvsnum_DrellYan_%d%s.pdf",RunYear,whichVFP.c_str()));
	else cancompare->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/compare_predvsnum_data_%d%s.pdf",RunYear,whichVFP.c_str()));
 
//a**************************   Toy MC *********************************a//
 
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
     
     	TMatrixDSym cov_vtx = result_fitvtx->GetCovarianceMatrix(); 
     	result_fitvtx->Print("V");     
     	TVectorD mean_vtx(2) ;
     	mean_vtx(0) = result_fitvtx->Parameter(0); 
     	mean_vtx(1) = result_fitvtx->Parameter(1);
     	RooRealVar central_p0("central_p0","", mean_vtx(0)-result_fitvtx->ParError(0), mean_vtx(0)+result_fitvtx->ParError(0));
     	RooRealVar central_p1("central_p1","",mean_vtx(1)-result_fitvtx->ParError(1), mean_vtx(1)+result_fitvtx->ParError(1));
    	resultfile << "central_p0 " << mean_vtx(0) << " err " << result_fitvtx->ParError(0) << std::endl;
    	resultfile << "central_p1 " << mean_vtx(1) << " err " << result_fitvtx->ParError(1) << std::endl;

    	fakeRateFile << "vtxconst " << mean_vtx(0)  << std::endl;
		fakeRateFile << "vtxslope " << mean_vtx(1)  << std::endl;
		RooMultiVarGaussian mvg_vtx("mvg_vtx","mvg_vtx",RooArgList(central_p0,central_p1),mean_vtx,cov_vtx);
		RooDataSet* toymcdata_vtx = mvg_vtx.generate(RooArgSet(central_p0,central_p1),NTOY);
		TF1 *h_toymc_vtx[NTOY];
	     
		ofstream myfile;
		if(doEB && doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_EB.txt",RunYear,whichVFP.c_str()));
                else if(!doEB && doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_FullEcal.txt",RunYear,whichVFP.c_str()));
                else if(doEB && !doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB.txt",RunYear,whichVFP.c_str()));
                else if(!doEB && !doDrellYan) myfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_FullEcal.txt",RunYear,whichVFP.c_str()));

		TH1D *p_scalefactor = new TH1D("p_scalefactor","fake rate scale factor; scale factor;",100,1000,2000);
	     
		float toyptvalue[graphPtBins][NTOY];
		float lowtoyptvalue[graphPtBins];
		float hightoyptvalue[graphPtBins];
		for(int ii(0); ii < graphPtBins; ii++){
			lowtoyptvalue[ii] = 1;
			hightoyptvalue[ii] = 0;
		}
		float toyvtxvalue[graphVtxBins][NTOY];
		float lowtoyvtxvalue[graphVtxBins];
		float hightoyvtxvalue[graphVtxBins];
		for(int ii(0); ii < graphVtxBins; ii++){
			lowtoyvtxvalue[ii] = 1;
			hightoyvtxvalue[ii] = 0;
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
	     
			double data4 = toymcdata_vtx->get(i)->getRealValue("central_p0");
			double data5 = toymcdata_vtx->get(i)->getRealValue("central_p1");
			toymchistname.str("");
			toymchistname << "h_toymc_vtx_" << i;
			h_toymc_vtx[i] = new TF1(toymchistname.str().c_str(), "pol1", 0,50);
			h_toymc_vtx[i]->SetParameters(data4, data5);
	     
			invmass_prednum->Reset();
			for(unsigned iEvt(0); iEvt < etreeEt.size(); iEvt++){
						float weight_eta(0);
						for(unsigned i(0); i<nEtaBins; i++){
							if(i+1 < nEtaBins){
								if(fabs(etreeEta[iEvt]) > EtaBins[i] && fabs(etreeEta[iEvt]) < EtaBins[i+1])weight_eta = eta_num[i]/eta_den[i];
							}
							else{ 
								if(fabs(etreeEta[iEvt]) > EtaBins[i])weight_eta = eta_num[i]/eta_den[i];
							}
						}
				float w_ele = h_toymc_pt[i]->Eval(etreeEt[iEvt])
				            * (h_toymc_vtx[i]->Eval(etreeVtx[iEvt]) / avg_fr_vtx)
				            * (avg_fr_eta > 0 ? weight_eta / avg_fr_eta : 1.0);
				invmass_prednum->Fill(etreeInvmass[iEvt], w_ele); 
			}

//			invmass_prednum->Sumw2();
			if(invmass_prednum->Integral(lowBinNumber,highBinNumber) > 0 && invmass_prednum->Integral(lowBinNumber,highBinNumber) < 1e20){
				myfile << random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber) << " " << data1 << " " << data2 << " " << data3 <<  " " << data4 << " " << data5 << std::endl;
				p_scalefactor->Fill(random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber));
				for(int ibin(1); ibin <= graphPtBins; ibin++){
					double estimated = h_toymc_pt[i]->Eval(MINPT+0.5*ibin);
					toyptvalue[ibin-1][i] = estimated;
					if(estimated < lowtoyptvalue[ibin-1])lowtoyptvalue[ibin-1] = estimated;
					if(estimated > hightoyptvalue[ibin-1])hightoyptvalue[ibin-1] = estimated;
				}
				for(int ibin(1); ibin <= graphVtxBins; ibin++){
					double estimated = h_toymc_vtx[i]->Eval(0.5*ibin);
					toyvtxvalue[ibin-1][i] = estimated;
					if(estimated < lowtoyvtxvalue[ibin-1])lowtoyvtxvalue[ibin-1] = estimated;
					if(estimated > hightoyvtxvalue[ibin-1])hightoyvtxvalue[ibin-1] = estimated;
				}
			}
		}
	     
	     //*************************   Calculated errors  *************************************************************************************************************** //
		for(int ibin(0); ibin < graphPtBins; ibin++){
			TH1D *h_toyptdis = new TH1D("h_toyptdis","",50,lowtoyptvalue[ibin],hightoyptvalue[ibin]);
			for(unsigned i(0); i < NTOY; i++)h_toyptdis->Fill(toyptvalue[ibin][i]);
			h_toyptdis->Fit("gaus");
			//float fiterror = h_toyptdis->GetFunction("gaus")->GetParameter(2);
			float fiterror = fabs(hightoyptvalue[ibin] - lowtoyptvalue[ibin])/2;
			cout<<" fiterror  "<< fiterror<<endl;
			fr_pt_sigmaband->SetPoint(ibin, MINPT+0.5*ibin,fit_fakerate_pt->Eval(MINPT+0.5*ibin)); 
			fr_pt_sigmaband->SetPointError(ibin, 0.25, fiterror);
			fr_pt_ratioError->SetPoint(ibin, MINPT+0.5*ibin, 1);
			fr_pt_ratioError->SetPointError(ibin, 0.25, fiterror/fit_fakerate_pt->Eval(MINPT+0.5*ibin));
			cout<< ibin << " "<< fiterror/fit_fakerate_pt->Eval(MINPT+0.5*ibin) <<endl;
			delete h_toyptdis; 
		}
		for(int ibin(0); ibin < graphVtxBins; ibin++){
			TH1D *h_toyvtxdis = new TH1D("h_toyvtxdis","",50,lowtoyvtxvalue[ibin],hightoyvtxvalue[ibin]);
			for(unsigned i(0); i < NTOY; i++)h_toyvtxdis->Fill(toyvtxvalue[ibin][i]);
			h_toyvtxdis->Fit("gaus");
			//float fiterror = h_toyvtxdis->GetFunction("gaus")->GetParameter(2);
			float fiterror = (hightoyvtxvalue[ibin] - lowtoyvtxvalue[ibin])/2;
			fr_vtx_sigmaband->SetPoint(ibin, MINVTX+0.5*ibin,fit_fakerate_vtx->Eval(MINVTX+0.5*ibin)); 
			fr_vtx_sigmaband->SetPointError(ibin, 0.25, fiterror);
			fr_vtx_ratioError->SetPoint(ibin, MINVTX+0.5*ibin, 1);
			fr_vtx_ratioError->SetPointError(ibin, 0.25, fiterror/fit_fakerate_vtx->Eval(MINVTX+0.5*ibin));
			delete h_toyvtxdis;
		}

	// toy bands of the nVtx-in-pT-slices plot: per pT slice and nVtx the half spread (max-min)/2 of the toy lines
	// (the VtxByPtSlice toys below) about the nominal slice fit, as for the inclusive band above. The pT-in-nVtx-slices
	// plot has no band: its slice fits have 3 parameters for 3 points (ndf = 0), and no toys are used for them.
	std::vector<TGraphErrors*> vtxSliceBand(pTSliceDirs.size(), nullptr);
	auto MakeToyBand = [](TF1 *nominal, const std::vector<std::vector<double> > &toyPars, double xlo, double xhi, double step){
		int n = (int)((xhi - xlo)/step) + 1;
		TGraphErrors *band = new TGraphErrors(n);
		TF1 *f = (TF1*)nominal->Clone(Form("%s_toy", nominal->GetName()));
		for(int i = 0; i < n; i++){
			double x = xlo + i*step, lo = 1e9, hi = -1e9;
			for(auto &p : toyPars){ f->SetParameters(p.data()); double v = f->Eval(x); lo = std::min(lo, v); hi = std::max(hi, v); }
			band->SetPoint(i, x, nominal->Eval(x));
			band->SetPointError(i, 0, toyPars.empty() ? 0 : (hi - lo)/2);
		}
		delete f;
		return band;
	};

	//a***********************   Toy MC for the VtxByPtSlice model   ***********************a//
	// The toys above vary the old factorized model. These toys vary the model actually used in
	// the background prediction (EvalFakeRate_VtxByPtSlice): the inclusive pT shape (a, b, alpha),
	// the nVtx line (c_s, beta_s) of each pT slice within its fit covariance (same +-1 sigma
	// truncation as above), and the observed signal yield of each slice within its fit error.
	// The slice normalisation N_s is recomputed per toy as for the nominal value.
	// One line per toy: [N_s c_s beta_s <f(pT)>_s] x 3 slices, then a b alpha.
	{
		const int NTOYSLICE = 1000;   // = NTOY in pred_eleBkg_pT_nVtx.C / analysis_eleBkg_pT_nVtx.C
		const double PTLO = 35.0, PTHI = 1000.0, PTBINW = 0.5;
		const int nPtB = (int)((PTHI-PTLO)/PTBINW);
		const size_t nS = pTSliceDirs.size();
		// per slice and pT bin: sum of eta weights (S0) and of eta weight x nVtx (S1) of the
		// denominator events inside the mass window used for the normalisation
		std::vector<std::vector<double> > S0(nS, std::vector<double>(nPtB,0.)), S1(nS, std::vector<double>(nPtB,0.));
		for(unsigned iEvt(0); iEvt < etreeEt.size(); iEvt++){
			int is = SelectPtSliceIndex(etreeEt[iEvt]);
			if(is < 0) continue;
			int mbin = invmass_num->FindBin(etreeInvmass[iEvt]);
			if(mbin < lowBinNumber || mbin > highBinNumber) continue;
			float weight_eta(0);
			for(unsigned k(0); k<nEtaBins; k++){
				if(k+1 < nEtaBins){
					if(fabs(etreeEta[iEvt]) > EtaBins[k] && fabs(etreeEta[iEvt]) < EtaBins[k+1])weight_eta = eta_num[k]/eta_den[k];
				}
				else if(fabs(etreeEta[iEvt]) > EtaBins[k])weight_eta = eta_num[k]/eta_den[k];
			}
			double etaw = (avg_fr_eta > 0 ? weight_eta / avg_fr_eta : 1.0);
			int ib = std::min(nPtB-1, std::max(0, (int)((etreeEt[iEvt]-PTLO)/PTBINW)));
			S0[is][ib] += etaw;
			S1[is][ib] += etaw*etreeVtx[iEvt];
		}

		RooDataSet* toyPt = mvg_pt.generate(RooArgSet(central_slope,central_const,central_index), 2*NTOYSLICE);
		std::vector<RooRealVar*> sliceP0(nS), sliceP1(nS);
		std::vector<RooMultiVarGaussian*> sliceMvg(nS);
		std::vector<RooDataSet*> toyVtx(nS);
		for(size_t s=0; s<nS; s++){
			// slice fit, or the inclusive nVtx fit if the slice fit failed (as for the nominal file)
			TFitResultPtr r = (fit_fakerate_vtx_pTSlice[s] && result_vtx_pTSlice[s].Get()) ? result_vtx_pTSlice[s] : result_fitvtx;
			TVectorD m(2); m(0) = r->Parameter(0); m(1) = r->Parameter(1);
			sliceP0[s] = new RooRealVar(Form("slice%d_p0",(int)s), "", m(0)-r->ParError(0), m(0)+r->ParError(0));
			sliceP1[s] = new RooRealVar(Form("slice%d_p1",(int)s), "", m(1)-r->ParError(1), m(1)+r->ParError(1));
			TMatrixDSym cov = r->GetCovarianceMatrix();
			sliceMvg[s] = new RooMultiVarGaussian(Form("mvg_slice%d",(int)s), "", RooArgList(*sliceP0[s],*sliceP1[s]), m, cov);
			toyVtx[s] = sliceMvg[s]->generate(RooArgSet(*sliceP0[s],*sliceP1[s]), 2*NTOYSLICE);
		}

		std::string toyName;
		if(doEB && doDrellYan)        toyName = Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_EB_VtxByPtSlice.txt",RunYear,whichVFP.c_str());
		else if(!doEB && doDrellYan)  toyName = Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_FullEcal_VtxByPtSlice.txt",RunYear,whichVFP.c_str());
		else if(doEB && !doDrellYan)  toyName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB_VtxByPtSlice.txt",RunYear,whichVFP.c_str());
		else                          toyName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_FullEcal_VtxByPtSlice.txt",RunYear,whichVFP.c_str());
		std::ofstream toySliceFile(toyName.c_str());

		TF1 fpt_toy("fpt_toy_slice", fakerate_ptDependence, 30, 1000, 3);
		std::vector<std::vector<std::vector<double> > > vtxToyPars(nS);   // (c_s, beta_s) of the written toys
		int nWritten = 0;
		for(int it = 0; it < 2*NTOYSLICE && nWritten < NTOYSLICE; it++){
			double a0 = toyPt->get(it)->getRealValue("central_slope");
			double a1 = toyPt->get(it)->getRealValue("central_const");
			double a2 = toyPt->get(it)->getRealValue("central_index");
			fpt_toy.SetParameters(a0, a1, a2);
			std::ostringstream line;
			bool ok = true;
			for(size_t s=0; s<nS && ok; s++){
				double lo = ptSliceLo[s], hi = std::min(ptSliceHi[s], (double)MAXPT);
				double avg = (hi > lo) ? fpt_toy.Integral(lo, hi)/(hi - lo) : 0;
				double c = toyVtx[s]->get(it)->getRealValue(sliceP0[s]->GetName());
				double b = toyVtx[s]->get(it)->getRealValue(sliceP1[s]->GetName());
				double pred = 0;
				for(int ib = 0; ib < nPtB; ib++){
					if(S0[s][ib] == 0 && S1[s][ib] == 0) continue;
					pred += fpt_toy.Eval(PTLO + (ib+0.5)*PTBINW) * (c*S0[s][ib] + b*S1[s][ib]);
				}
				if(avg > 0) pred /= avg;
				double num = numSignalSlice_v[s] + gRandom->Gaus(0, numSignalSliceErr_v[s]);
				if(!(avg > 0) || !(pred > 0) || !(num > 0)){ ok = false; break; }
				line << num/pred << " " << c << " " << b << " " << avg << " ";
			}
			if(!ok) continue;
			line << a0 << " " << a1 << " " << a2;
			toySliceFile << line.str() << std::endl;
			for(size_t s=0; s<nS; s++)
				vtxToyPars[s].push_back({toyVtx[s]->get(it)->getRealValue(sliceP0[s]->GetName()), toyVtx[s]->get(it)->getRealValue(sliceP1[s]->GetName())});
			nWritten++;
		}
		// nVtx band of each pT slice (only for slices with their own fit, which are the ones drawn)
		for(size_t s=0; s<nS; s++)
			if(fit_fakerate_vtx_pTSlice[s]) vtxSliceBand[s] = MakeToyBand(fit_fakerate_vtx_pTSlice[s], vtxToyPars[s], fit_fakerate_vtx_pTSlice[s]->GetXmin(), fit_fakerate_vtx_pTSlice[s]->GetXmax(), 0.5);
		toySliceFile.close();
		std::cout << "VtxByPtSlice toys written: " << nWritten << " -> " << toyName << std::endl;
		resultfile << "VtxByPtSlice toys written: " << nWritten << std::endl;
		if(nWritten < NTOYSLICE) std::cout << "WARNING: fewer VtxByPtSlice toys than " << NTOYSLICE << std::endl;
		for(size_t s=0; s<nS; s++){ delete toyVtx[s]; delete sliceMvg[s]; delete sliceP0[s]; delete sliceP1[s]; }
		delete toyPt;
	}

	     
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
		fr_bothcount_vtx->SetMarkerStyle(20);
		fr_bothcount_vtx->SetMarkerColor(kBlue);
		fr_bothcount_vtx->SetLineColor(kBlue);
		fr_bothcount_vtx->SetLineWidth(2);
		fr_bothcount_vtx->SetFillColor(0);
		fr_pt_ratio->SetMarkerStyle(20);
		fr_pt_ratio->SetMarkerColor(kBlue);
		fr_pt_ratio->SetLineColor(kBlue);
		fr_pt_ratio->SetLineWidth(2);
		fr_pt_ratio->SetFillColor(0);
		fr_vtx_ratio->SetMarkerStyle(20);
		fr_vtx_ratio->SetMarkerColor(kBlue);
		fr_vtx_ratio->SetLineColor(kBlue);
		fr_vtx_ratio->SetLineWidth(2);
		fr_vtx_ratio->SetFillColor(0);
			TLine *flatratio_pt = new TLine(30,1,200,1);
			TLine *flatratio_eta= new TLine(0,1,1.5,1);
			TLine *flatratio_vtx= new TLine(0,1,46,1);

		std::vector<std::string> nVtxSubLabels = {"nVtx 8-20","nVtx 20-30","nVtx >30"};
		std::vector<int> nVtxSubColors         = {kBlack,kGreen+2,kOrange+7};
		std::vector<std::string> &nVtxSubDirs  = nVtxSliceDirs; // reuse the slices already fit above
		std::vector<TGraphErrors*> &ptSubGraphs = ptSliceGraphs;

		std::vector<std::string> &pTSubDirs   = pTSliceDirs; // reuse the slices already fit above (moved up)
		std::vector<std::string> pTSubLabels  = {"p_{T} 30-50","p_{T} 50-80","p_{T} >80"};
		std::vector<int> pTSubColors          = {kBlack,kGreen+2,kOrange+7};
		std::vector<TGraphErrors*> &vtxSubGraphs = vtxSliceGraphs;

		for(size_t i=0;i<ptSubGraphs.size();i++){
			if(!ptSubGraphs[i]) continue;
			ptSubGraphs[i]->SetMarkerStyle(20);
			ptSubGraphs[i]->SetMarkerColor(nVtxSubColors[i]);
			ptSubGraphs[i]->SetMarkerSize(1.1);
			ptSubGraphs[i]->SetLineColor(nVtxSubColors[i]);
			if(fit_fakerate_pt_nVtxSlice[i]) fit_fakerate_pt_nVtxSlice[i]->SetLineColor(nVtxSubColors[i]);
		}
		for(size_t i=0;i<vtxSubGraphs.size();i++){
			if(!vtxSubGraphs[i]) continue;
			vtxSubGraphs[i]->SetMarkerStyle(20);
			vtxSubGraphs[i]->SetMarkerColor(pTSubColors[i]);
			vtxSubGraphs[i]->SetMarkerSize(1.1);
			vtxSubGraphs[i]->SetLineColor(pTSubColors[i]);
			if(fit_fakerate_vtx_pTSlice[i]) fit_fakerate_vtx_pTSlice[i]->SetLineColor(pTSubColors[i]);
		}

		TLegend *legPt = new TLegend(0.5,0.55,0.89,0.89);
		for(size_t i=0;i<ptSubGraphs.size();i++)
			if(ptSubGraphs[i]) legPt->AddEntry(ptSubGraphs[i], nVtxSubLabels[i].c_str(), "p");

		TLegend *legVtx = new TLegend(0.5,0.55,0.89,0.89);
		for(size_t i=0;i<vtxSubGraphs.size();i++)
			if(vtxSubGraphs[i]) legVtx->AddEntry(vtxSubGraphs[i], pTSubLabels[i].c_str(), "p");
		legVtx->AddEntry((TObject*)0, "hatched: toy #pm1#sigma", "");
	     
		TCanvas *canpt = new TCanvas("canpt","",600,600);
		canpt->cd();
		TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.3, 1, 1.0);
		canpt_pad1->SetBottomMargin(0);
		canpt_pad1->Draw();          
		canpt_pad1->cd();          
		TH1D *dummy_pt = new TH1D("",";Pt(GeV);fake rate",nPtBins*10,MINPT,MAXPT);
		dummy_pt->SetMaximum(0.06);
		dummy_pt->SetMinimum(0.0);
		//dummy_pt->SetMaximum(fit_fakerate_pt->Eval(30)*1.2);
		dummy_pt->GetYaxis()->SetTitleOffset(1.4);
		dummy_pt->Draw();
		//fr_bothcount_pt->Draw("EP same");
		//fit_fakerate_pt->Draw("same");
		fr_pt_sigmaband->SetFillColor(15);
		fr_pt_sigmaband->SetFillStyle(3345);
		//fr_pt_sigmaband->Draw("E2 same");
		//fr_bothcount_pt->Draw("EP same");
		for(size_t i=0;i<ptSubGraphs.size();i++)
			if(ptSubGraphs[i]) ptSubGraphs[i]->Draw("P same");
		for(size_t i=0;i<fit_fakerate_pt_nVtxSlice.size();i++)
			if(fit_fakerate_pt_nVtxSlice[i]) fit_fakerate_pt_nVtxSlice[i]->Draw("same");
		legPt->Draw("same");
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
		TH1D *dummy_ptratio = new TH1D("dummy_ptratio",";p_{T}(GeV);data/fit",nPtBins*10,MINPT,MAXPT);
		dummy_ptratio->SetMaximum(1.7);
		dummy_ptratio->SetMinimum(0.5);
		//dummy_ptratio->SetMaximum(2.0);
		//dummy_ptratio->SetMinimum(0.0);
		dummy_ptratio->GetYaxis()->SetNdivisions(504);
		dummy_ptratio->Draw();
		/*fr_pt_ratio->Draw("EP same");
		fr_pt_ratioError->SetFillColor(15);
		fr_pt_ratioError->SetFillStyle(3345);
		fr_pt_ratioError->Draw("E2 same");
		fr_pt_ratio->Draw("EP same");
		flatratio_pt->Draw("same"); */
		if(doDrellYan)canpt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_pt_systematic_DrellYan_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		else canpt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_pt_systematic_data_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		  
		TCanvas *canvtx = new TCanvas("canvtx","",600,600);
		canvtx->cd();
		TPad *canvtx_pad1 = new TPad("canvtx_pad1", "pad1", 0, 0.3, 1, 1.0);
		canvtx_pad1->SetBottomMargin(0);
		canvtx_pad1->Draw();          
		canvtx_pad1->cd();          
		TH1D *dummy_vtx = new TH1D("",";nVtx;fake rate",MAXVTX,0,MAXVTX);
		dummy_vtx->SetMaximum(0.1);
	//	dummy_vtx->SetMinimum(-0.005);
		dummy_vtx->Draw();
	//	fr_bothcount_vtx->Draw("EP same");
	//	fit_fakerate_vtx->Draw("same");
		fr_vtx_sigmaband->SetFillColor(15);
		fr_vtx_sigmaband->SetFillStyle(3345);
	//	fr_vtx_sigmaband->Draw("E2 same");
	//	fr_bothcount_vtx->Draw("EP same");
		const int bandStyles[3] = {3354, 3345, 3359};
		for(size_t i=0;i<vtxSliceBand.size();i++){
			if(!vtxSliceBand[i]) continue;
			vtxSliceBand[i]->SetFillColor(pTSubColors[i]);
			vtxSliceBand[i]->SetFillStyle(bandStyles[i%3]);
			vtxSliceBand[i]->SetLineColor(pTSubColors[i]);
			vtxSliceBand[i]->Draw("E3 same");
		}
		for(size_t i=0;i<vtxSubGraphs.size();i++)
			if(vtxSubGraphs[i]) vtxSubGraphs[i]->Draw("P same");
		for(size_t i=0;i<fit_fakerate_vtx_pTSlice.size();i++)
			if(fit_fakerate_vtx_pTSlice[i]) fit_fakerate_vtx_pTSlice[i]->Draw("same");
		legVtx->Draw("same");
		if(RunYear==2016 and preVFP == 1)       CMS_lumi( canvtx_pad1, 1,1, 11 );
        	else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canvtx_pad1, 2,1, 11 );
        	else if(RunYear==2017)                  CMS_lumi( canvtx_pad1, 3,1, 11 );
        	else if(RunYear==2018)                  CMS_lumi( canvtx_pad1, 4,1, 11 );

	        
		canvtx->cd();   
		TPad *canvtx_pad2 = new TPad("canvtx_pad2", "pad2", 0, 0, 1, 0.3);
		canvtx_pad2->SetTopMargin(0);
		canvtx_pad2->SetBottomMargin(0.5);
		canvtx_pad2->Draw();
		canvtx_pad2->cd(); 	
		TH1D *dummy_vtxratio = new TH1D("dummy_vtxratio",";nVtx;data/fit",MAXVTX,0,MAXVTX);
		dummy_vtxratio->SetMaximum(1.7);
		dummy_vtxratio->SetMinimum(0.5);
				dummy_vtxratio->GetYaxis()->SetNdivisions(504);
		dummy_vtxratio->Draw();
		/*fr_vtx_ratio->Draw("EP same");
		fr_vtx_ratioError->SetFillColor(15);
		fr_vtx_ratioError->SetFillStyle(3345);
		fr_vtx_ratioError->Draw("E2 same");
		flatratio_vtx->Draw("same");
		fr_vtx_ratio->Draw("EP same"); */
		if(doDrellYan)canvtx->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_vtx_systematic_DrellYan_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		else canvtx->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_vtx_systematic_data_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));  
		
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
		if(doDrellYan)caneta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_systematic_DrellYan_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		else{
			caneta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_systematic_data_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
			caneta->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_eta_systematic_data_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		}
	     
		TCanvas *canscale = new TCanvas("scale","scale",600,600);
		p_scalefactor->Draw("hist");
		if(doDrellYan)canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_scalefactor_DrellYan_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));
		else canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_scalefactor_data_%d%s%s.pdf",RunYear,whichVFP.c_str(),isEB.c_str()));		

				resultfile.close();
				fakeRateFile.close();
}
