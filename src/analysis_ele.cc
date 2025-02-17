#include "../include/analysis_ele.h"

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
    //case SingleElectron2017: if(fireTrgs(12))passHLT = true; break; // HLT_Ele27_WPTight_Gsf

    case MCSingleElectron2016: if(fireTrgs(12))passHLT = true; break; // HLT_Ele27_WPTight_Gsf
    case MCSingleElectron2017: if(fireTrgs(46))passHLT = true; break; // HLT_Ele35_WPTight_Gsf
    case MCSingleElectron2018: if(fireTrgs(13))passHLT = true; break; // HLT_Ele32_WPTight_Gsf  
    //case MCSingleElectron2017: if(fireTrgs(12))passHLT = true; break; // HLT_Ele27_WPTight_Gsf

/*  case SingleElectron2016: if(fireTrgs(10))passHLT = true; break; // HLT_Ele27_eta2p1_WPLoose_Gsf
    case SingleElectron2017: if(fireTrgs(42))passHLT = true; break; // HLT_Ele20_eta2p1_WPLoose_Gsf
    case SingleElectron2018: if(fireTrgs(42))passHLT = true; break; // HLT_Ele20_eta2p1_WPLoose_Gsf
    case MCSingleElectron2016: if(fireTrgs(10))passHLT = true; break;
    case MCSingleElectron2017: if(fireTrgs(42))passHLT = true; break;
    case MCSingleElectron2018: if(fireTrgs(42))passHLT = true; break;
*/
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
  bool passMiniMedium(true);
 
  if(!isEB() && !isEE()) {passMiniMedium=false; return passMiniMedium; }
  if(isEB()){
	if(getSigma() > 0.0106) {passMiniMedium=false; return passMiniMedium; }
	if(fabs(getdEtaIn()) > 0.0032){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getdPhiIn()) > 0.0547){passMiniMedium=false; return passMiniMedium;}
	if(getHoverE() > getCutValueHoverE()){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getEoverPInv()) > 0.184){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getD0()) > 0.05){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getDz()) > 0.10){passMiniMedium=false; return passMiniMedium; }
	if(getMissHits() > 1){passMiniMedium=false; return passMiniMedium;}
	if(getConvVeto() == 0){passMiniMedium=false; return passMiniMedium;}
	if(getMiniIso() > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  else if(isEE()){
	if(getSigma() > 0.0387) {passMiniMedium=false; return passMiniMedium; }
	if(fabs(getdEtaIn()) > 0.00632){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getdPhiIn()) > 0.0394){passMiniMedium=false; return passMiniMedium;}
	if(getHoverE() > getCutValueHoverE()){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getEoverPInv()) > 0.0721){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getD0()) > 0.10){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getDz()) > 0.20){passMiniMedium=false; return passMiniMedium; }
	if(getMissHits() > 1){passMiniMedium=false; return passMiniMedium;}
	if(getConvVeto() == 0){passMiniMedium=false; return passMiniMedium;}
	if(getMiniIso() > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
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

// mengleis
bool recoEle::isFakeProxyOLD(){
        bool passFakeProxy(true);

  if(!isEB() && !isEE()) {passFakeProxy=false; return passFakeProxy; }
  if(isEB()){
        if(getHoverE() > 0.253){passFakeProxy=false; return passFakeProxy;}
        if(fabs(getEoverPInv()) > 0.134){passFakeProxy=false; return passFakeProxy;}
        if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
        if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
  else if(isEE()){
        if(getHoverE() > 0.0878){passFakeProxy=false; return passFakeProxy;}
        if(fabs(getEoverPInv()) > 0.13){passFakeProxy=false; return passFakeProxy;}
        if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
        if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }

  if(isEB()){
        if(getSigma() > 0.00998 || fabs(getdEtaIn()) > 0.00311 || fabs(getdPhiIn()) > 0.103 || getMiniIso() > 0.1)passFakeProxy=true;
        else passFakeProxy=false;
  }
  else if(isEE()){
        if(getSigma() > 0.0298 || fabs(getdEtaIn()) > 0.00609 || fabs(getdPhiIn()) > 0.045 || getMiniIso() > 0.1)passFakeProxy=true;
        else passFakeProxy=false;
  }
//Upper bound:  Veto
  if(isEB()){
        if(getSigma() > 0.0115 || fabs(getdEtaIn()) > 0.00749 || fabs(getdPhiIn()) > 0.228)passFakeProxy=false;
  }
  else if(isEE()){
        if(getSigma() > 0.037  || fabs(getdEtaIn()) > 0.00895 || fabs(getdPhiIn()) > 0.213)passFakeProxy=false;
  }
        return passFakeProxy;
}

// from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
bool recoEle::isFakeProxy(){
	bool passFakeProxy(true);

  if(!isEB() && !isEE()) {passFakeProxy=false; return passFakeProxy; }
  if(isEB()){
	if(getHoverE() > getCutValueHoverE()){passFakeProxy=false; return passFakeProxy;}  // this cut reduces, ele-fake-photon by 30% wrt 80X old cut 
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

//Upper bound:  Veto cuts
  if(isEB()){
	if(getSigma() > 0.0126 || fabs(getdEtaIn()) > 0.00463 || fabs(getdPhiIn()) > 0.148)passFakeProxy=false;
  }
  else if(isEE()){
	if(getSigma() > 0.0457  || fabs(getdEtaIn()) > 0.00814 || fabs(getdPhiIn()) > 0.19)passFakeProxy=false;
  }
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
