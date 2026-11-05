#include "../../../include/analysis_commoncode.h"

enum bkgType{
	elefakepho = 1,
	jetfakepho = 2,
	qcdfakelep = 3,
	VGamma = 4,
	rare = 5
};	

TFile *file_prefit;
double fit_norm[36];
double fit_error[36]; 
/**********************************/
/*	double normfactor = par[0]; 	*/  
/*	double slope = par[1];				*/
/*	double constant = par[2];			*/
/*	double index = par[3];				*/
/*	double coeff = par[4]; 				*/
/*	double vtx_constant = par[5];	*/
/*	double vtx_slope = par[6];		*/
/**********************************/
double scalefactor(-99999);
double ptslope(-99999);
double ptconstant(-99999);
double ptindex(-99999);
double ptcoeff(-99999);
double vtxconst(-99999);
double vtxslope(-99999);
TF3 h_nominal_fakerate("h_nominal_fakerate", fakerate_func,10,1000,0,100,0,1.5,7);
/**** Jet Fake Photon ****/
TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,10000,4);
TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,10000,4);
/**** Fake Lepton ****/
double factorQCD(1);
TH1D *p_scale;
/**** V+Gamma ****/
double factorMC(1);
/**** all MC ****/
esfScaleFactor  objectESF;


void init_weight(int ichannel, bkgType bkg){ 

	std::ifstream elefake_file("../../script/EleFakeRate-ByPtVtx-EB.txt");
	std::string variabletype;
	double variablevalue;
	if(elefake_file.is_open()){
  	for(int i(0); i<7; i++){ 
			elefake_file >> variabletype >> variablevalue; 
			if(variabletype.find("scalefactor")!=std::string::npos)scalefactor = variablevalue;
			else if(variabletype.find("ptslope")!=std::string::npos)ptslope = variablevalue;
			else if(variabletype.find("ptconstant")!=std::string::npos)ptconstant = variablevalue;
			else if(variabletype.find("ptindex")!=std::string::npos)ptindex = variablevalue;
			else if(variabletype.find("ptcoeff")!=std::string::npos)ptcoeff = variablevalue;
			else if(variabletype.find("vtxconst")!=std::string::npos)vtxconst = variablevalue;
			else if(variabletype.find("vtxslope")!=std::string::npos)vtxslope = variablevalue;
	  }
	}
	elefake_file.close();
	h_nominal_fakerate.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	std::stringstream JetFakeRateFile;
  JetFakeRateFile.str();
	//if(ichannel==1)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-DoubleEG-EB-15.txt";
	//if(ichannel==2)JetFakeRateFile << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-transferfactor-MuonEG-EB-15.txt";
	if(ichannel==1)JetFakeRateFile << "../../script/JetFakeRate-transferfactor-DoubleEG-EB.txt";
	if(ichannel==2)JetFakeRateFile << "../../script/JetFakeRate-transferfactor-MuonEG-EB.txt";
	std::ifstream jetfakefile(JetFakeRateFile.str().c_str());
	std::string paratype;
	float paravalue;	
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_den->SetParameter(i, paravalue);
		std::cout << paratype << " " << paravalue << std::endl;
	}
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_num->SetParameter(i, paravalue);
		std::cout << paratype << " " << paravalue << std::endl;
	}
 
	if(ichannel == 1){
		factorQCD = factor_egQCD;
	}
	else if(ichannel == 2){
		factorQCD = factor_mgQCD;
	}
	TFile *scaleFile;
	if(ichannel == 1)scaleFile = TFile::Open("../Result/qcd_eg_scale.root");
	else if(ichannel == 2)scaleFile = TFile::Open("../Result/qcd_mg_scale.root");
	if(ichannel == 1)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	else if(ichannel == 2)p_scale = (TH1D*)scaleFile->Get("transfer_factor");

	if(ichannel == 1){
		factorMC = factor_egVGamma;
	}
	else if(ichannel == 2){
		factorMC = factor_mgVGamma;
	}

	std::ostringstream postfitname;
	postfitname.str("");
	if(bkg == bkgType::elefakepho)postfitname << "NormPostFit_elefakepho.txt";
	else if(bkg == bkgType::jetfakepho)postfitname << "NormPostFit_jetfakepho.txt";
	else if(bkg == bkgType::qcdfakelep)postfitname << "NormPostFit_qcdfakelep.txt";
	else if(bkg == bkgType::VGamma)postfitname << "NormPostFit_VGamma.txt";
	else if(bkg == bkgType::rare)postfitname << "NormPostFit_rare.txt";
	
	std::ifstream postfit_file(postfitname.str().c_str());
  double postnorm, posterr;
	double prenom, preerror;
  if(postfit_file.is_open()){
    for(int i(0); i<36; i++){ 
      postfit_file >> prenom >> preerror >> postnorm >> posterr;
      fit_norm[i] = postnorm/prenom;
      fit_error[i]= posterr/preerror;
    }
  }

  std::ostringstream prefitname;  
	prefitname.str("");
	prefitname <<  "../Result/vetoDiPhoData/";
	if(ichannel == 1){
		if(bkg == bkgType::elefakepho)prefitname <<  "signalTree_egamma_eleBkg.root";
		else if(bkg == bkgType::jetfakepho)prefitname <<  "signalTree_egamma_jetbkg.root";
		else if(bkg == bkgType::qcdfakelep)prefitname <<  "signalTree_egamma_qcd.root";
		else if(bkg == bkgType::VGamma)prefitname <<  "signalTree_egamma_VGBkg.root";
		else if(bkg == bkgType::rare)prefitname <<  "signalTree_egamma_rareBkg.root";
	}
	else if(ichannel == 2){
		if(bkg == bkgType::elefakepho)prefitname <<  "signalTree_mg_eleBkg.root";
		else if(bkg == bkgType::jetfakepho)prefitname <<  "signalTree_mg_jetbkg.root";
		else if(bkg == bkgType::qcdfakelep)prefitname <<  "signalTree_mg_qcd.root";
		else if(bkg == bkgType::VGamma)prefitname <<  "signalTree_mg_VGBkg.root";
		else if(bkg == bkgType::rare)prefitname <<  "signalTree_mg_rareBkg.root";
	}
  file_prefit = TFile::Open(prefitname.str().c_str());
}

double getWeight(int ichannel, bkgType bkg, double phoEt, double nVertex, double phoEta, double lepPt, double lepEta){
	double eventweight(0);
	if(bkg == bkgType::elefakepho){
		eventweight =  h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
	}
	else if(bkg == bkgType::qcdfakelep){
		if(ichannel == 1){
			eventweight = factorQCD*p_scale->GetBinContent(p_scale->FindBin(lepPt));
		}
		else{
			eventweight = factorQCD;
		}
	}
	else if(bkg == bkgType::jetfakepho){
		if(phoEt < 300)eventweight = fitfunc_num->Eval(phoEt)/fitfunc_den->Eval(phoEt);
		else eventweight = fitfunc_num->Eval(300)/fitfunc_den->Eval(300);
	}
	else if(bkg == bkgType::VGamma){
		double scalefactor(0);
		if(ichannel == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
		}
		if(ichannel == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
		}
		eventweight = scalefactor*factorMC;
	}
	else if(bkg == bkgType::rare){
		double scalefactor(0);
		if(ichannel == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
		}
		if(ichannel == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
		}
		eventweight = scalefactor;
	}

	return eventweight;
}
