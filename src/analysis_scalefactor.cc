#include "../include/analysis_scalefactor.h"

float esfScaleFactor::getElectronESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDESF->GetXaxis();
		etaaxis= electronIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDESF->GetYaxis();
		etaaxis= electronIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleESF(1);
	if(ptOnX){
		eleESF = electronIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		//eleESF = electronIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta))*electronISOESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleESF = electronIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	//else eleESF = electronIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt))*electronISOESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));

	eleESF = eleESF*getR9ESF(pt, eta);	
	
  return eleESF;
}


float esfScaleFactor::getPhotonESF(float inputpt, float inputeta){    // new
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonIDESF->GetXaxis();
		etaaxis= photonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonIDESF->GetYaxis();
		etaaxis= photonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoESF(1);
	if(ptOnX){
		phoESF = photonIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoESF = photonIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoESF;
}


float esfScaleFactor::getMuonESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDESF->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDESF->GetXaxis();
		etaaxis= muonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDESF->GetYaxis();
		etaaxis= muonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 120)pt = 119;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.4) eta = 2.39;

	float muonESF(1);
	float muonIDESFvalue(1);
	//float muonISOESFvalue(1), muonIPESFvalue(1);
	if(ptOnX){
		muonIDESFvalue = muonIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));  // new available
		//muonISOESFvalue= muonISOESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		//muonIPESFvalue = muonIPESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		//muonESF = muonIDESFvalue*muonISOESFvalue*muonIPESFvalue;
		muonESF = muonIDESFvalue;
	}
	else{
	 muonIDESFvalue = muonIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt)); // new available
	 //muonISOESFvalue = muonISOESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 //muonIPESFvalue = muonIPESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 //muonESF = muonIDESFvalue*muonISOESFvalue*muonIPESFvalue;
	 muonESF = muonIDESFvalue;
	}
		
  return muonESF;
}

float esfScaleFactor::getegPhotonTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGESF->GetXaxis();
		etaaxis= photonTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGESF->GetYaxis();
		etaaxis= photonTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGESF(1);
	if(ptOnX){
		phoTRGESF = photonTRGESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGESF = photonTRGESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGESF;
}

float esfScaleFactor::getElectronTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGESF->GetXaxis();
		etaaxis= electronTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGESF->GetYaxis();
		etaaxis= electronTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleTRGESF(1);
	if(ptOnX){
		eleTRGESF = electronTRGESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleTRGESF = electronTRGESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleTRGESF;
}

float esfScaleFactor::getMuonEGTRGESF(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGESF->GetXaxis();
	TAxis *muptaxis  = muonTRGESF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGESF = muonTRGESF->GetBinContent(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGESF; 
}


float esfScaleFactor::getElectronESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDESF->GetXaxis();
		etaaxis= electronIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDESF->GetYaxis();
		etaaxis= electronIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleESFError(1);
	if(ptOnX){
		eleESFError = sqrt(pow(electronIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)), 2) + pow(getR9ESFError(pt, eta),2));
		//eleESFError = sqrt(pow(electronIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)), 2)+ pow(electronISOESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)), 2) + pow(getR9ESFError(pt, eta),2));
	}
	else eleESFError = sqrt(pow( electronIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)), 2) + pow(getR9ESFError(pt, eta),2));
	//else eleESFError = sqrt(pow( electronIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)), 2)+ pow(electronISOESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)),2) + pow(getR9ESFError(pt, eta),2));
	
		
  return eleESFError;
}


float esfScaleFactor::getPhotonESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonIDESF->GetXaxis();
		etaaxis= photonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonIDESF->GetYaxis();
		etaaxis= photonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoESFError(1);
	if(ptOnX){
		phoESFError = photonIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoESFError = photonIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoESFError;
}


float esfScaleFactor::getMuonESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDESF->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDESF->GetXaxis();
		etaaxis= muonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDESF->GetYaxis();
		etaaxis= muonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 120)pt = 119;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float muonESFError(1);
	float muonIDESFError(1), muonISOESFError(1), muonIPESFError(1);
	if(ptOnX){
		muonIDESFError = muonIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonISOESFError = muonISOESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonIPESFError = muonIPESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonESFError = sqrt(pow(muonIDESFError,2)+pow(muonISOESFError+muonIPESFError,2));
	}
	else{
	  muonIDESFError = muonIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonISOESFError = muonISOESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonIPESFError = muonIPESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		muonESFError = sqrt(pow(muonIDESFError,2)+pow(muonISOESFError+muonIPESFError,2));
	}
		
  return muonESFError;
}

float esfScaleFactor::getegPhotonTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGESF->GetXaxis();
		etaaxis= photonTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGESF->GetYaxis();
		etaaxis= photonTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGESFError(1);
	if(ptOnX){
		phoTRGESFError = photonTRGESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGESFError = photonTRGESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGESFError;
}

float esfScaleFactor::getElectronTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGESF->GetXaxis();
		etaaxis= electronTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGESF->GetYaxis();
		etaaxis= electronTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 500)pt = 499;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float electronTRGESFError(1);
	if(ptOnX){
		electronTRGESFError = electronTRGESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else electronTRGESFError = electronTRGESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return electronTRGESFError;
}

float esfScaleFactor::getMuonEGTRGESFError(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGESF->GetXaxis();
	TAxis *muptaxis  = muonTRGESF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGESFError = muonTRGESF->GetBinError(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGESFError; 
}


//////////////        Fast Sim     ///////////////////////////////////////


float esfScaleFactor::getFastElectronESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDFast->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDFast->GetXaxis();
		etaaxis= electronIDFast->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDFast->GetYaxis();
		etaaxis= electronIDFast->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleFast(1);
	if(ptOnX){
		eleFast = getElectronESF(pt,eta)*electronIDFast->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta))*electronISOFast->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleFast = getElectronESF(pt,eta)*electronIDFast->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt))*electronISOFast->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleFast;
}


float esfScaleFactor::getFastMuonESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDFast->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDFast->GetXaxis();
		etaaxis= muonIDFast->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDFast->GetYaxis();
		etaaxis= muonIDFast->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float muonFast(1);
	float muonIDFastvalue(1), muonISOFastvalue(1), muonIPFastvalue(1);
	if(ptOnX){
		muonIDFastvalue = muonIDFast->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonISOFastvalue= muonISOFast->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonIPFastvalue = muonIPFast->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonFast = getMuonESF(pt, eta)*muonIDFastvalue*muonISOFastvalue*muonIPFastvalue;
	}
	else{
	 muonIDFastvalue = muonIDFast->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonISOFastvalue = muonISOFast->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonIPFastvalue = muonIPFast->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonFast = getMuonESF(pt, eta)*muonIDFastvalue*muonISOFastvalue*muonIPFastvalue;
	}
		
  return muonFast;
}

float esfScaleFactor::getFastegPhotonTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGEFF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGEFF->GetXaxis();
		etaaxis= photonTRGEFF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGEFF->GetYaxis();
		etaaxis= photonTRGEFF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGEFF(1);
	if(ptOnX){
		phoTRGEFF = photonTRGEFF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGEFF = photonTRGEFF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGEFF;
}

float esfScaleFactor::getFastElectronTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGEFF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGEFF->GetXaxis();
		etaaxis= electronTRGEFF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGEFF->GetYaxis();
		etaaxis= electronTRGEFF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleTRGEFF(1);
	if(ptOnX){
		eleTRGEFF = electronTRGEFF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleTRGEFF = electronTRGEFF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleTRGEFF;
}

float esfScaleFactor::getFastMuonEGTRGESF(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGEFF->GetXaxis();
	TAxis *muptaxis  = muonTRGEFF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGEFF = muonTRGEFF->GetBinContent(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGEFF; 
}


float esfScaleFactor::getFastElectronESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDFast->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDFast->GetXaxis();
		etaaxis= electronIDFast->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDFast->GetYaxis();
		etaaxis= electronIDFast->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleFastError(1);
	if(ptOnX){
		eleFastError = pow(electronIDFast->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)),2)+ pow(electronISOFast->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)),2);
	}
	else eleFastError = pow(electronIDFast->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)),2)+pow(electronISOFast->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)),2);
	
	eleFastError += pow(getElectronESFError(pt, eta),2);
	eleFastError = sqrt(eleFastError);
  return eleFastError;
}

float esfScaleFactor::getFastMuonESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDFast->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDFast->GetXaxis();
		etaaxis= muonIDFast->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDFast->GetYaxis();
		etaaxis= muonIDFast->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float muonFastError(1);
	float muonIDFastError(1), muonISOFastError(1), muonIPFastError(1);
	if(ptOnX){
		muonIDFastError = muonIDFast->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonISOFastError = muonISOFast->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonIPFastError = muonIPFast->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonFastError = pow(muonIDFastError,2)+pow(muonISOFastError,2)+pow(muonIPFastError,2);
	}
	else{
	  muonIDFastError = muonIDFast->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonISOFastError = muonISOFast->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonIPFastError = muonIPFast->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		muonFastError = pow(muonIDFastError,2)+pow(muonISOFastError,2)+pow(muonIPFastError,2);
	}

	muonFastError += pow(getMuonESFError(pt, eta),2);
	muonFastError = sqrt(muonFastError);
		
  return muonFastError;
}

float esfScaleFactor::getFastegPhotonTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGEFF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGEFF->GetXaxis();
		etaaxis= photonTRGEFF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGEFF->GetYaxis();
		etaaxis= photonTRGEFF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGEFFError(1);
	if(ptOnX){
		phoTRGEFFError = photonTRGEFF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGEFFError = photonTRGEFF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGEFFError;
}

float esfScaleFactor::getFastElectronTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGEFF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGEFF->GetXaxis();
		etaaxis= electronTRGEFF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGEFF->GetYaxis();
		etaaxis= electronTRGEFF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float electronTRGEFFError(1);
	if(ptOnX){
		electronTRGEFFError = electronTRGEFF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else electronTRGEFFError = electronTRGEFF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return electronTRGEFFError;
}

float esfScaleFactor::getFastMuonEGTRGESFError(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGEFF->GetXaxis();
	TAxis *muptaxis  = muonTRGEFF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGEFFError = muonTRGEFF->GetBinError(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGEFFError; 
}


float esfScaleFactor::getR9ESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronR9ESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronR9ESF->GetXaxis();
		etaaxis= electronR9ESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronR9ESF->GetYaxis();
		etaaxis= electronR9ESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleESF(1);
	if(ptOnX){
		eleESF = electronR9ESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleESF = electronR9ESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleESF;
}

float esfScaleFactor::getR9ESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDESF->GetXaxis();
		etaaxis= electronIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDESF->GetYaxis();
		etaaxis= electronIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt >= 200)pt = 199;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleR9ESFError(1);
	if(ptOnX){
		eleR9ESFError = sqrt(pow(electronR9ESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta)), 2));
	}
	else eleR9ESFError = sqrt(pow(electronR9ESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt)), 2));
		
  return eleR9ESFError;
}

