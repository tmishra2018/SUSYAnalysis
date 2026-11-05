#include "../include/analysis_ele.h"
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V

bool recoEle::passSignalSelection(){
  
  bool passCut(true);
  
  if(Calibp4_.Pt() < 25.0){passCut = false; return passCut;} 
  if(!isEB() && !isEE()){passCut = false; return passCut;}
  if(!isMiniMedium()){passCut = false; return passCut;}
  else passCut = true;

  return passCut; 
}

bool recoEle::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;
    case DoubleEG2015: case MuonEG2015: case SingleMuon2015: case DoubleMuon2015: case SingleElectron2015:
      passHLT = false; break;
    case DoubleEG2016:  if(fireTrgs(21) || fireTrgs(22))passHLT = true; break; // HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    case DoubleEG2017:  if(fireTrgs(43) || fireTrgs(44))passHLT = true; break; // HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v13
    case DoubleEG2018:  if(fireTrgs(43) || fireTrgs(44))passHLT = true; break; // HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v13

    case MCDoubleEG2016: if(fireTrgs(21) || fireTrgs(22))passHLT = true; break;
    case MCDoubleEG2017: if(fireTrgs(43) || fireTrgs(44))passHLT = true; break;
    case MCDoubleEG2018: if(fireTrgs(43) || fireTrgs(44))passHLT = true; break;

    case MuonEG2016: passHLT = true; break;
    case MuonEG2017: passHLT = true; break;
    case MuonEG2018: passHLT = true; break;
   
    case MCMuonEG2016: passHLT = true; break;
    case MCMuonEG2017: passHLT = true; break;
    case MCMuonEG2018: passHLT = true; break;
    
    case SingleElectron2016: if(fireTrgs(12))passHLT = true; break; // HLT_Ele27_WPTight_Gsf
    case SingleElectron2017: if(fireTrgs(46))passHLT = true; break; // HLT_Ele35_WPTight_Gsf
    case SingleElectron2018: if(fireTrgs(13))passHLT = true; break; // HLT_Ele32_WPTight_Gsf  

    case MCSingleElectron2016: if(fireTrgs(12))passHLT = true; break; // HLT_Ele27_WPTight_Gsf
    case MCSingleElectron2017: if(fireTrgs(46))passHLT = true; break; // HLT_Ele35_WPTight_Gsf
    case MCSingleElectron2018: if(fireTrgs(13))passHLT = true; break; // HLT_Ele32_WPTight_Gsf  

    case SingleMuon2016: passHLT = true; break;
    case DoubleMuon2016: passHLT = true; break;
    case MCSingleMuon2016: passHLT = true; break;
    case MCDoubleMuon2016: passHLT = true; break;
    case MCMET2016: passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}

bool recoEle::passBasicID(){
  bool passBasic(true);
 
  if(!isEB() && !isEE()) {passBasic=false; return passBasic; }
  if(isEB()){
		if(getHoverE() > getCutValueHoverE()){passBasic=false; return passBasic;}
		if(fabs(getEoverPInv()) > 0.184){passBasic=false; return passBasic;}
		if(getMissHits() > 1){passBasic=false; return passBasic;}
		if(getConvVeto() == 0){passBasic=false; return passBasic;}
  }
  else if(isEE()){
		if(getHoverE() > getCutValueHoverE()){passBasic=false; return passBasic;}
		if(fabs(getEoverPInv()) > 0.0721){passBasic=false; return passBasic;}
		if(getMissHits() > 1){passBasic=false; return passBasic;}
		if(getConvVeto() == 0){passBasic=false; return passBasic;}
  }
  return passBasic;
}

bool recoEle::isMiniMedium(){
  bool passMiniMedium(false);
  if(!isEB() && !isEE()) return passMiniMedium;
  const double esc = Calibp4_.E();
  const double pt  = Calibp4_.Pt();
  if(esc <= 0 || pt <= 0) return passMiniMedium;
  const double absDeta = fabs(getdEtaIn());
  const double absDphi = fabs(getdPhiIn());
  const double absEopInv = fabs(getEoverPInv());
  if(isEB()){
    	const double hoeMedium = 0.046 + 1.16/esc + 0.0324*getRho()/esc;
    	passMiniMedium = (getSigma() <= 0.0106 &&
                      absDeta <= 0.0032 &&
                      absDphi <= 0.0547 &&
                      getHoverE() <= hoeMedium &&
                      getMiniIso() <= 0.1 &&
                      absEopInv <= 0.184 &&
                      getMissHits() <= 1 &&
                      getConvVeto());}
   else if(isEE()){
    	const double hoeMedium = 0.0275 + 2.52/esc + 0.183*getRho()/esc;
    	passMiniMedium = (getSigma() <= 0.0387 &&
                      absDeta <= 0.00632 &&
                      absDphi <= 0.0394 &&
                      getHoverE() <= hoeMedium &&
                      getMiniIso() <= 0.1 &&
                      absEopInv <= 0.0721 &&
                      getMissHits() <= 1 &&
                      getConvVeto());}
    return passMiniMedium;
}

bool recoEle::isMiniLoose(){
  bool passMiniLoose(true);
 
  if(!isEB() && !isEE()) {passMiniLoose=false; return passMiniLoose; }
  if(isEB()){
	if(getSigma() > 0.0112) {passMiniLoose=false; return passMiniLoose; }
	if(fabs(getdEtaIn()) > 0.00377){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getdPhiIn()) > 0.0884){passMiniLoose=false; return passMiniLoose;}
	if(getHoverE() > 0.298){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getEoverPInv()) > 0.193){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getD0()) > 0.05){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getDz()) > 0.10){passMiniLoose=false; return passMiniLoose; }
	if(getMissHits() > 1){passMiniLoose=false; return passMiniLoose;}
	if(!getConvVeto()){passMiniLoose=false; return passMiniLoose;}
	if(getMiniIso() > 0.2){passMiniLoose=false; return passMiniLoose;}
  }
  else if(isEE()){
	if(getSigma() > 0.0425) {passMiniLoose=false; return passMiniLoose; }
	if(fabs(getdEtaIn()) > 0.00674){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getdPhiIn()) > 0.169){passMiniLoose=false; return passMiniLoose;}
	if(getHoverE() > 0.101){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getEoverPInv()) > 0.111){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getD0()) > 0.10){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getDz()) > 0.20){passMiniLoose=false; return passMiniLoose; }
	if(getMissHits() > 1){passMiniLoose=false; return passMiniLoose;}
	if(!getConvVeto()){passMiniLoose=false; return passMiniLoose;}
	if(getMiniIso() > 0.2){passMiniLoose=false; return passMiniLoose;}
  }
  return passMiniLoose;
}

// from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
 bool recoEle::isFakeProxy(){
	bool passFakeProxy(false);
  	if(!isEB() && !isEE()) return passFakeProxy;

  	const double esc = Calibp4_.E();
  	const double pt  = Calibp4_.Pt();
  	if(esc <= 0 || pt <= 0) return passFakeProxy;

  	const double relIso = getCalibRelIso();
  	const double absDeta = fabs(getdEtaIn());
  	const double absDphi = fabs(getdPhiIn());
  	const double absEopInv = fabs(getEoverPInv());

  	bool passMedium(false);
  	bool passVeto(false);
 
   	if(isEB()){
    		const double hoeMedium = 0.046 + 1.16/esc + 0.0324*getRho()/esc;
    		const double hoeVeto   = 0.05  + 1.16/esc + 0.0324*getRho()/esc;

    		passMedium = (getSigma() <= 0.0106 &&
                  	absDeta <= 0.0032 &&
                  	absDphi <= 0.0547 &&
                  	getHoverE() <= hoeMedium &&
                  	relIso <= (0.0478 + 0.506/pt) &&
                  	absEopInv <= 0.184 &&
                  	getMissHits() <= 1 &&
                  	getConvVeto());

    		passVeto = (getSigma() <= 0.0126 &&
                	absDeta <= 0.00463 &&
                	absDphi <= 0.148 &&
                	getHoverE() <= hoeVeto &&
                	relIso <= (0.198 + 0.506/pt) &&
                	absEopInv <= 0.209 &&
                	getMissHits() <= 2 &&
                	getConvVeto());
   	}
   	else if(isEE()){
    			const double hoeMedium = 0.0275 + 2.52/esc + 0.183*getRho()/esc;
    			const double hoeVeto   = 0.05   + 2.54/esc + 0.183*getRho()/esc;

    			passMedium = (getSigma() <= 0.0387 &&
                  		absDeta <= 0.00632 &&
                  		absDphi <= 0.0394 &&
                  		getHoverE() <= hoeMedium &&
                  		relIso <= (0.0658 + 0.963/pt) &&
                  		absEopInv <= 0.0721 &&
                  		getMissHits() <= 1 &&
                  		getConvVeto());

    			passVeto = (getSigma() <= 0.0457 &&
                		absDeta <= 0.00814 &&
                		absDphi <= 0.19 &&
                		getHoverE() <= hoeVeto &&
                		relIso <= (0.203 + 0.963/pt) &&
                		absEopInv <= 0.132 &&
                		getMissHits() <= 3 &&
                		getConvVeto());
  	}

  	// Fake proxy: fail 94X Medium but remain within 94X Veto envelope.
  	passFakeProxy = passVeto && !passMedium;
 	return passFakeProxy;
 }

bool recoEle::isLooseFakeProxy(){
	bool passFakeProxy(true);

  if(!isEB() && !isEE()) {passFakeProxy=false; return passFakeProxy; }
  if(isEB()){
	if(getHoverE() > getCutValueHoverE()){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.184){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
  else if(isEE()){
	if(getHoverE() > getCutValueHoverE()){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.0721){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
	
  if(isEB()){
	if(getSigma() > 0.0106 || fabs(getdEtaIn()) > 0.0032 || fabs(getdPhiIn()) > 0.0547 || getMiniIso() > 0.1)passFakeProxy=true;
  	else passFakeProxy=false;
  }
  else if(isEE()){
	if(getSigma() > 0.0387 || fabs(getdEtaIn()) > 0.00632 || fabs(getdPhiIn()) > 0.0394 || getMiniIso() > 0.1)passFakeProxy=true;
  	else passFakeProxy=false;
  }

	return passFakeProxy;

}
