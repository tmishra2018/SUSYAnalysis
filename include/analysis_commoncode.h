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
#include "TVector2.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_jet.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"
#include "../include/analysis_fakes.h"
#include "../include/analysis_binning.h"
#include "../include/analysis_scalefactor.h"
#include "../include/tdrstyle.C"

//float lumi_2016preVFP_DoubleEG = 19.499256;
//float lumi_2016postVFP_DoubleEG = 16.810813;
//float lumi_2017_DoubleEG = 41.4613;
//float lumi_2018_DoubleEG = 59.816229;

//float lumi_2016preVFP_MuonEG = 19.495077;
//float lumi_2016postVFP_MuonEG = 16.810813;
//float lumi_2017_MuonEG = 27.100211; // considering only DEF eras, because of availabilty of trigger
//float lumi_2018_MuonEG = 59.810698;

int ichannel(1);
int anatype(-1);
int lowMt(0);
int highMt(-1);
int lowMET(0);
int highMET(-1);
int lowPt(25);
int highPt(-1);
int lepIso(4);

//int NBIN(32);
//float METbin1(200), METbin2(350), METbin3(500);
//float HTbin1(100),  HTbin2(300), HTbin3(500);
//float PHOETbin(200), PHOETBin2(200);

//int NBIN(24);
//float METbin1(200), METbin2(300), METbin3(400);
//float HTbin1(100),  HTbin2(250), HTbin3(400);
//float PHOETbin(100), PHOETBin2(200);

int NBIN(18);
float METbin1(200), METbin2(400);
float HTbin1(100),  HTbin2(400);
float PHOETbin(200);
int RunYear(2016);
bool preVFP(1);
//mengleis

double factor_egQCD(0.237);
double factorerror_egQCD(0.0474);
double factor_mgQCD(0.62);
double factorerror_mgQCD(0.124);

double factor_egVGamma(1.17);
double factorerror_egVGamma(0.234);
double factor_mgVGamma(1.33);
double factorerror_mgVGamma(0.266);

// for 2016 preVFP

double factor_egQCD_2016preVFP(0.398154);
double factorerror_egQCD_2016preVFP(0.032191);
double factor_mgQCD_2016preVFP(0.912074);
double factorerror_mgQCD_2016preVFP(0.0919012);

double factor_egVGamma_2016preVFP(1.54809);
double factorerror_egVGamma_2016preVFP(0.0931608);
double factor_mgVGamma_2016preVFP(1.40733);
double factorerror_mgVGamma_2016preVFP(0.0684939);

// for 2016 postVFP

double factor_egQCD_2016postVFP(0.511481);
double factorerror_egQCD_2016postVFP(0.0344402);
double factor_mgQCD_2016postVFP(0.980161);
double factorerror_mgQCD_2016postVFP(0.131539);

double factor_egVGamma_2016postVFP(1.39163);
double factorerror_egVGamma_2016postVFP(0.112449);
double factor_mgVGamma_2016postVFP(1.45883);
double factorerror_mgVGamma_2016postVFP(0.117483);

// for 2017

double factor_egQCD_2017(0.403763);
double factorerror_egQCD_2017(0.0264701);
double factor_mgQCD_2017(0.926775);
double factorerror_mgQCD_2017(0.136297);

double factor_egVGamma_2017(1.45821);
double factorerror_egVGamma_2017(0.101554);
double factor_mgVGamma_2017(1.3442);
double factorerror_mgVGamma_2017(0.11364);

// for 2018


double factor_egQCD_2018(0.357835);
double factorerror_egQCD_2018(0.0227174);
double factor_mgQCD_2018(0.902506);
double factorerror_mgQCD_2018(0.145422);

double factor_egVGamma_2018(1.73086);
double factorerror_egVGamma_2018(0.0858579);
double factor_mgVGamma_2018(1.32145);
double factorerror_mgVGamma_2018(0.0842323);

bool SetRunConfig(){
	
	std::ifstream configfile("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/test/Background/BkgPredConfig.txt");
	std::string conftype;
	double confvalue;
	if(configfile.is_open()){
  	for(int i(0); i<11; i++){ 
			configfile >> conftype >> confvalue; 
			if(conftype.find("ichannel")!=std::string::npos)ichannel = confvalue;
			if(conftype.find("anatype")!=std::string::npos)anatype = confvalue;
			if(conftype.find("lowMt")!=std::string::npos)lowMt = confvalue;
			if(conftype.find("highMt")!=std::string::npos)highMt = confvalue;
			if(conftype.find("lowMET")!=std::string::npos)lowMET = confvalue;
			if(conftype.find("highMET")!=std::string::npos)highMET = confvalue;
			if(conftype.find("lowPt")!=std::string::npos)lowPt = confvalue;
			if(conftype.find("highPt")!=std::string::npos)highPt = confvalue;
			if(conftype.find("lepIso")!=std::string::npos)lepIso = confvalue;
			if(conftype.find("RunYear")!=std::string::npos)RunYear = confvalue;
			if(conftype.find("preVFP")!=std::string::npos)preVFP = confvalue;
	  }
	}
	configfile.close();
	
	if(anatype >= 0)return true;
	else return false;
}

bool SetSignalConfig(){
	
	std::ifstream configfile("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/test/Result/SigConfig.txt");
	std::string conftype;
	double confvalue;
	if(configfile.is_open()){
  	for(int i(0); i<11; i++){ 
			configfile >> conftype >> confvalue; 
			if(conftype.find("ichannel")!=std::string::npos)ichannel = confvalue;
			if(conftype.find("anatype")!=std::string::npos)anatype = confvalue;
			if(conftype.find("lowMt")!=std::string::npos)lowMt = confvalue;
			if(conftype.find("highMt")!=std::string::npos)highMt = confvalue;
			if(conftype.find("lowMET")!=std::string::npos)lowMET = confvalue;
			if(conftype.find("highMET")!=std::string::npos)highMET = confvalue;
			if(conftype.find("lowPt")!=std::string::npos)lowPt = confvalue;
			if(conftype.find("highPt")!=std::string::npos)highPt = confvalue;
			if(conftype.find("lepIso")!=std::string::npos)lepIso = confvalue;
			if(conftype.find("RunYear")!=std::string::npos)RunYear = confvalue;
			if(conftype.find("preVFP")!=std::string::npos)preVFP = confvalue;
	  }
	}
	configfile.close();

	std::ifstream binfile("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/test/Result/binConfig.txt");
	if(binfile.is_open()){
		for(int i(0); i<6; i++){
			binfile >> conftype >> confvalue;
			if(conftype.find("NBIN")!=std::string::npos)NBIN = confvalue;
			if(conftype.find("METbin1")!=std::string::npos)METbin1= confvalue;
			if(conftype.find("METbin2")!=std::string::npos)METbin2= confvalue;
			if(conftype.find("HTbin1")!=std::string::npos)HTbin1= confvalue;
			if(conftype.find("HTbin2")!=std::string::npos)HTbin2= confvalue;
			if(conftype.find("PHOETbin")!=std::string::npos)PHOETbin= confvalue;
		}
	}
	
	if(anatype == 3)return true;
	else return false;
}

double calcToyError(vector<double> & vtoy, bool useGauss, int ic){
	
		double syserr(0);
		double normvalue = vtoy[0];
		std::sort(vtoy.begin(), vtoy.end());
		TH1D *temphist = new TH1D("temphist","",500,0.2*normvalue,1.8*normvalue);
		for(unsigned it(0); it < vtoy.size(); it++){
			temphist->Fill(vtoy[it]);
		}
		if(useGauss){
			TCanvas *can=new TCanvas("can","",600,600);
			can->cd();
			temphist->Draw();
			temphist->Fit("gaus");
			syserr = temphist->GetFunction("gaus")->GetParameter(2);
			std::ostringstream canname;
			canname.str("");
			canname << ic << ".png";
			can->SaveAs(canname.str().c_str());	
		}
		else syserr = std::max( fabs(vtoy.front() - normvalue), fabs(vtoy.back() - normvalue));
		delete temphist;

		return syserr;
}


double isrW(double isrPt){
		double reweightF = 1;
    		if(isrPt < 50)reweightF = 1.015;
    		else if(isrPt >= 50 && isrPt < 100)reweightF  = 1.110;
    		else if(isrPt >= 100 && isrPt < 150)reweightF = 0.845;
    		else if(isrPt >= 150 && isrPt < 200)reweightF = 0.715;
    		else if(isrPt >= 200 && isrPt < 250)reweightF = 0.730;
    		else if(isrPt >= 250 && isrPt < 300)reweightF = 0.732;
    		else if(isrPt >= 300)reweightF =  0.642;
		return reweightF;
}
