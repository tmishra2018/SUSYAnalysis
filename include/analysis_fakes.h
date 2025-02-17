#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2016preVFP/result_eta_dependence_Data.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2016postVFP/result_eta_dependence_Data.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2017/result_eta_dependence_Data.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2018/result_eta_dependence_Data.txt"

#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2016preVFP/result_eta_dependence_Data_EE.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2016postVFP/result_eta_dependence_Data_EE.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2017/result_eta_dependence_Data_EE.txt"
#include	"/eos/uscms/store/user/tmishra/elefakepho/DataResult2018/result_eta_dependence_Data_EE.txt"


// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2016preVFP/result_eta_dependence_DY.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2016postVFP/result_eta_dependence_DY.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2017/result_eta_dependence_DY.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2018/result_eta_dependence_DY.txt"

// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2016preVFP/result_eta_dependence_DY_EE.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2016postVFP/result_eta_dependence_DY_EE.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2017/result_eta_dependence_DY_EE.txt"
// #include        "/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult2018/result_eta_dependence_DY_EE.txt"


std::ifstream configfile("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/myResult/SigConfig.txt");
std::string conftype;
double confvalue;

int runyear = 2016;
bool prevfp =  1;
Double_t fakerate_func(Double_t *x, Double_t *par)
{
	if(configfile.is_open()){
    		for(int i(0); i<11; i++){
                        configfile >> conftype >> confvalue;
                        if(conftype.find("RunYear")!=std::string::npos)runyear = confvalue;
                        if(conftype.find("preVFP")!=std::string::npos)prevfp = confvalue;
          	}
	}
	configfile.close();
	double weight_pt(0.02);
	double weight_nvtx(0.02);
	double weight_eta(0.02);

	double normfactor = par[0];  // scalefactor
  	double slope = par[1];
  	double constant = par[2];
  	double index = par[3];
  	double coeff = 1.0; 
	double vtx_constant = par[5];
	double vtx_slope = par[6];
	
  	double pt = TMath::Max(x[0],0.000001);
	double nvtx=x[1];
	double eta= x[2];

	double arg = 0;
	arg = slope*pt + constant; 
	double fitval = pow(arg, index)*coeff;  // (constant + slope.pt)^(-index)
	weight_pt = fitval;
	weight_nvtx = vtx_constant + vtx_slope*nvtx; // (const + slope*nvtx)

	if(eta >= 0 && eta < 1.4442){
		for(int ieta(0); ieta < 29; ieta++)
			if(eta > ieta*0.05 && eta <= (ieta+1)*0.05)
			{
				if (runyear==2016 and prevfp==1) 	weight_eta = etaRatesEB_2016preVFP[ieta];
				else if (runyear==2016 and prevfp==0) 	weight_eta = etaRatesEB_2016postVFP[ieta];
				else if (runyear==2017) 		weight_eta = etaRatesEB_2017[ieta];
				else if (runyear==2018) 		weight_eta = etaRatesEB_2018[ieta];
			}
	}
	else if(eta > 1.56 && eta <= 2.5){
		for(int ieta(0); ieta < 94; ieta++)
			if(eta > 1.56 + ieta*0.01 && eta <= 1.56+(ieta+1)*0.01){
				if (runyear==2016 and prevfp==1) 	weight_eta = etaRatesEE_2016preVFP[ieta];
				else if (runyear==2016 and prevfp==0) 	weight_eta = etaRatesEE_2016postVFP[ieta];
				else if (runyear==2017) 		weight_eta = etaRatesEE_2017[ieta];
				else if (runyear==2018) 		weight_eta = etaRatesEE_2018[ieta];
			}
	}
	else weight_eta = 0;
	double totalfakerate= normfactor*weight_pt*weight_nvtx*weight_eta;
	return totalfakerate;
}

/*Double_t mcfakerate_func(Double_t *x, Double_t *par)
{
	double weight_pt(0.02);
	double weight_nvtx(0.02);
	double weight_eta(0.02);

	double normfactor = par[0];
  	double slope = par[1];
  	double constant = par[2];
  	double index = par[3];
  	double coeff = 1.0; 
	double vtx_constant = par[5];
	double vtx_slope = par[6];
  	double pt = TMath::Max(x[0],0.000001);
	double nvtx=x[1];
	double eta= x[2];

	double arg = 0;
	arg = slope*pt + constant; 
	double fitval = pow(arg, index)*coeff; 
	weight_pt = fitval;
	weight_nvtx = vtx_constant + vtx_slope*nvtx;

	if(eta >= 0 && eta < 1.4442){
		for(int ieta(0); ieta < 29; ieta++)
			if(eta > ieta*0.05 && eta <= (ieta+1)*0.05){
				if (runyear==2016 and prevfp==1) 	weight_eta = MCetaRatesEB_2016preVFP[ieta];
				else if (runyear==2016 and prevfp==0) 	weight_eta = MCetaRatesEB_2016postVFP[ieta];
				else if (runyear==2017) 		weight_eta = MCetaRatesEB_2017[ieta];
				else if (runyear==2018) 		weight_eta = MCetaRatesEB_2018[ieta];}
	}
	else if(eta > 1.56 && eta <= 2.5){
		for(int ieta(0); ieta < 94; ieta++)
			if(eta > 1.56 + ieta*0.01 && eta <= 1.56+(ieta+1)*0.01){
				if (RunYear==2016 and preVFP==1)        	weight_eta = MCetaRatesEE_2016preVFP[ieta];
				else if (RunYear==2016 and preVFP==0)        	weight_eta = MCetaRatesEE_2016postVFP[ieta];
				else if (RunYear==2017)        			weight_eta = MCetaRatesEE_2017[ieta];
				else if (RunYear==2018)        			weight_eta = MCetaRatesEE_2018[ieta];}
	}
	else weight_eta = 0;

	double totalfakerate= normfactor*weight_pt*weight_nvtx*weight_eta;
	return totalfakerate;
}
*/

//Double_t jetfake_func(Double_t *x, Double_t *par)
//{
//	double pt_low = x[0] - 1.0/2.0;
//	double pt_high = x[0]+ 1.0/2.0;
//
//	double c1 = par[0];
//	double c2 = par[1];
//	double lamda1 = par[2];
//	double lamda2 = par[3];
//
//	double jetfakes_lowedge = c1*exp(lamda1*pt_low)/lamda1 + c2*exp(lamda2*pt_low)/lamda2;
//	double jetfakes_highedge =  c1*exp(lamda1*pt_high)/lamda1 + c2*exp(lamda2*pt_high)/lamda2;
//	//return (jetfakes_highedge + jetfakes_lowedge)/2.0*REBINSIZE;
//	return (jetfakes_highedge - jetfakes_lowedge);
//}

Double_t jetfake_func(Double_t *x, Double_t *par)
{
	double pt_low = x[0] - 1.0/2.0;
	double pt_high = x[0]+ 1.0/2.0;

	double c1 = par[0];
	double c2 = par[1];
	double lamda1 = par[2];
	double lamda2 = par[3];

	double jetfakes_lowedge = c1*exp(lamda1*pt_low)/lamda1 + c2*exp(lamda2*pt_low)/lamda2;
	double jetfakes_highedge =  c1*exp(lamda1*pt_high)/lamda1 + c2*exp(lamda2*pt_high)/lamda2;
	//return (jetfakes_highedge + jetfakes_lowedge)/2.0*REBINSIZE;
	return (jetfakes_highedge - jetfakes_lowedge);
}

