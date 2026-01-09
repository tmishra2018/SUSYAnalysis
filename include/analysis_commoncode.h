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


int ichannel(1);
int anatype(-1);
int lowMt(0);
int highMt(-1);
int lowMET(0);
int highMET(-1);
int lowPt(25);
int highPt(-1);
int lepIso(4);

int NBIN(18);
float METbin1(200), METbin2(400);
float HTbin1(100),  HTbin2(400);
float PHOETbin(200);
int RunYear(2016);
bool preVFP(1);

double factor_mgQCD_2016preVFP(0.87261);
double factorerror_mgQCD_2016preVFP(0.0820623);
double factor_mgVGamma_2016preVFP(1.35197);
double factorerror_mgVGamma_2016preVFP(0.107935);
double factor_egQCD_2016preVFP(0.491149);
double factorerror_egQCD_2016preVFP(0.12104);
double factor_egVGamma_2016preVFP(1.42088);
double factorerror_egVGamma_2016preVFP(0.311331);

double factor_mgQCD_2016postVFP(0.919834);
double factorerror_mgQCD_2016postVFP(0.173638);
double factor_mgVGamma_2016postVFP(1.47747);
double factorerror_mgVGamma_2016postVFP(0.165919);
double factor_egQCD_2016postVFP(0.674543);
double factorerror_egQCD_2016postVFP(0.105477);
double factor_egVGamma_2016postVFP(1.13285);
double factorerror_egVGamma_2016postVFP(0.334537);

double factor_mgQCD_2017(0.915692);
double factorerror_mgQCD_2017(0.254747);
double factor_mgVGamma_2017(1.49914);
double factorerror_mgVGamma_2017(0.238776);
double factor_egQCD_2017(0.457467);
double factorerror_egQCD_2017(0.112141);
double factor_egVGamma_2017(1.18835);
double factorerror_egVGamma_2017(0.432986);

double factor_mgQCD_2018(0.859031);
double factorerror_mgQCD_2018(0.128368);
double factor_mgVGamma_2018(1.36389);
double factorerror_mgVGamma_2018(0.1488);
double factor_egQCD_2018(0.411959);
double factorerror_egQCD_2018(0.0950488);
double factor_egVGamma_2018(1.53025);
double factorerror_egVGamma_2018(0.349885);

bool SetRunConfig(){
	
	std::ifstream configfile("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/Background/BkgPredConfig.txt");
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
	
	std::ifstream configfile("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/myResult/SigConfig.txt");
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

	std::ifstream binfile("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/myResult/binConfig.txt");
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
			//std::ostringstream canname;
			//canname.str("");
			//canname << ic << ".png";
			//can->SaveAs(canname.str().c_str());	
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
