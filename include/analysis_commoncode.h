#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cstdlib>

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


int ichannel(1);
int anatype(-1);
int lowMt(0);
int highMt(-1);
int lowMET(0);
int highMET(-1);
int lowHT(0);
int highHT(-1);
int lowPt(25);
int highPt(-1);
int lepIso(4);
int lowPhoEt(0);    // optional photon-pT window (BkgPredConfig keys lowPhoEt/highPhoEt), used for the SF binning study
int highPhoEt(-1);
// V+gamma / fake-lepton SFs in the CR/VR/SR: photon-pT binned (VGamma_pho/QCDfake_pho) by default;
// environment variable SF_BINNING=HT selects the HT-binned factors (VGamma_ht/QCDfake_ht) instead
bool usePhoPtSF = !(getenv("SF_BINNING") && std::string(getenv("SF_BINNING")) == "HT");

int NBIN(18);
float METbin1(200), METbin2(400);
float HTbin1(100),  HTbin2(400);
float PHOETbin(200);
int RunYear(2016);
bool preVFP(1);


// inclusive V+gamma / fake-lepton factors (MET 40-70 GeV template fit, 2026-10-31 production): nominal fit and
// +1 sigma (toy 16-84% spread (+) MINOS, as VGamma/fit/plotVGammaScales.C). Used as the fixed factors / fallback
// for photon-pT bins without a stable fit (none in the 2026-10-31 production), and with SF_BINNING=HT.
double factor_mgQCD_2016preVFP(0.889385);
double factorerror_mgQCD_2016preVFP(0.0351834);
double factor_mgVGamma_2016preVFP(1.27792);
double factorerror_mgVGamma_2016preVFP(0.0333102);
double factor_mgQCD_2016postVFP(0.9233);
double factorerror_mgQCD_2016postVFP(0.0412466);
double factor_mgVGamma_2016postVFP(1.30636);
double factorerror_mgVGamma_2016postVFP(0.0436338);
double factor_mgQCD_2017(0.952944);
double factorerror_mgQCD_2017(0.0402262);
double factor_mgVGamma_2017(1.49174);
double factorerror_mgVGamma_2017(0.0520737);
double factor_mgQCD_2018(0.887432);
double factorerror_mgQCD_2018(0.0318337);
double factor_mgVGamma_2018(1.33581);
double factorerror_mgVGamma_2018(0.0374351);
double factor_egQCD_2016preVFP(0.734165);
double factorerror_egQCD_2016preVFP(0.0678965);
double factor_egQCD_2016postVFP(1.03671);
double factorerror_egQCD_2016postVFP(0.0751351);
double factor_egVGamma_2016postVFP(1.01211);
double factorerror_egVGamma_2016postVFP(0.13593);
double factor_egVGamma_2016preVFP(1.41235);
double factorerror_egVGamma_2016preVFP(0.121151);
double factor_egQCD_2017(0.765974);
double factorerror_egQCD_2017(0.057587);
double factor_egVGamma_2017(1.17028);
double factorerror_egVGamma_2017(0.13174);
double factor_egQCD_2018(0.668626);
double factorerror_egQCD_2018(0.059795);
double factor_egVGamma_2018(1.44105);
double factorerror_egVGamma_2018(0.141647);



// SUSYAnalysis area that contains the working directory: walk up until <dir>/include/analysis_commoncode.h exists.
// Works from any folder inside the area (test/..., limits/...), so a copied or moved area reads its own config files.
std::string anaBaseDir(){
	std::string dir = gSystem->WorkingDirectory();
	while(!dir.empty()){
		std::ifstream probe((dir + "/include/analysis_commoncode.h").c_str());
		if(probe.good()) return dir;
		size_t pos = dir.find_last_of('/');
		if(pos == std::string::npos) break;
		dir = dir.substr(0, pos);
	}
	std::cerr << "analysis_commoncode.h: no SUSYAnalysis area above " << gSystem->WorkingDirectory() << std::endl;
	return "";
}

bool SetRunConfig(){
	
	std::ifstream configfile((anaBaseDir() + "/test/Background/BkgPredConfig.txt").c_str());
	std::string conftype;
	double confvalue;
	if(configfile.is_open()){
  	while(configfile >> conftype >> confvalue){   // 13 standard keys, optionally lowPhoEt/highPhoEt
			if(conftype == "lowPhoEt")lowPhoEt = confvalue;
			if(conftype == "highPhoEt")highPhoEt = confvalue;
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
			if(conftype.find("lowHT")!=std::string::npos)lowHT = confvalue;
			if(conftype.find("highHT")!=std::string::npos)highHT = confvalue;
	  }
	}
	configfile.close();
	
	if(anatype >= 0)return true;
	else return false;
}

bool SetSignalConfig(){
	
	std::ifstream configfile((anaBaseDir() + "/test/myResult/SigConfig.txt").c_str());
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

	std::ifstream binfile((anaBaseDir() + "/test/myResult/binConfig.txt").c_str());
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
			temphist->Fit("gaus","Q0");
			syserr = temphist->GetFunction("gaus")->GetParameter(2);
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
