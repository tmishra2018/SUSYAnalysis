#include "../include/analysis_muon.h"

bool recoMuon::passSignalSelection(){
  bool passCut(true);
  if(p4_.Pt() < 20.0){passCut = false; return passCut;} // suggestion from convenors
  if(fabs(p4_.Eta())>2.40){passCut = false; return passCut;}
  if(!isMedium()){passCut = false; return passCut;}
  if(getMiniIso()>0.2){passCut = false; return passCut;}
  if(fabs(getD0()) > 0.02 || fabs(getDz()) > 0.1){passCut = false; return passCut;} // medium Prompt now, with dxy cut changed
  else passCut = true;
  return passCut; 
}
bool recoMuon::passSignalSelection_Not_MiniIso(){
  bool passCut(true);
  if(p4_.Pt() < 20.0){passCut = false; return passCut;} // suggestion from convenors
  if(fabs(p4_.Eta())>2.40){passCut = false; return passCut;}
  if(!isMedium()){passCut = false; return passCut;}
  if(fabs(getD0()) > 0.02 || fabs(getDz()) > 0.1){passCut = false; return passCut;}
  else passCut = true;
  return passCut; 
}

bool recoMuon::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;
    case DoubleEG2016: passHLT = true; break;
    case DoubleEG2017: passHLT = true; break;
    case DoubleEG2018: passHLT = true; break;
    case MCDoubleEG2016: passHLT = true; break;
    case MCDoubleEG2017: passHLT = true; break;
    case MCDoubleEG2018: passHLT = true; break;
 
    case MuonEG2016: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break; // HLT_Mu17_Photon || HLT_Mu17_Photon || HLT_Mu38NoFiltersNoVtx_Photon38 
    case MuonEG2017: if(fireSingleTrg(0) || fireSingleTrg(32))passHLT = true; break; // HLT_Mu17_Photon30_IsoCaloId || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL
    case MuonEG2018: if(fireSingleTrg(0) || fireSingleTrg(32))passHLT = true; break; // HLT_Mu17_Photon30_IsoCaloId || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL
 
    case MCMuonEG2016: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break;
    case MCMuonEG2017: if(fireSingleTrg(0) || fireSingleTrg(32))passHLT = true; break; // HLT_Mu17_Photon30_IsoCaloId || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL
    case MCMuonEG2018: if(fireSingleTrg(0) || fireSingleTrg(32))passHLT = true; break; // HLT_Mu17_Photon30_IsoCaloId || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL

    case SingleElectron2016: passHLT = true; break;
    case SingleElectron2017: passHLT = true; break;
    case SingleElectron2018: passHLT = true; break;

    case SingleMuon2016: if(fireSingleTrg(1) || fireSingleTrg(19))passHLT = true; break; // HLT_IsoMu24 || HLT_IsoTkMu24
    case SingleMuon2017: if(fireSingleTrg(31))passHLT = true; break; // HLT_IsoMu27
    case SingleMuon2018: if(fireSingleTrg(30))passHLT = true; break; // HLT_IsoMu24
 
    case MCSingleMuon2016: if(fireSingleTrg(1) || fireSingleTrg(19))passHLT = true; break; // HLT_IsoMu24 || HLT_IsoTkMu24
    case MCSingleMuon2017: if(fireSingleTrg(31))passHLT = true; break; // HLT_IsoMu27
    case MCSingleMuon2018: if(fireSingleTrg(30))passHLT = true; break; // HLT_IsoMu24
 
    case DoubleMuon2016: passHLT = true; break;
    case MCSingleElectron2016:passHLT = true; break;
    case MCDoubleMuon2016: passHLT = true; break;
    case MCMET2016: passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}

bool recoMuon::isFakeProxy(){
	bool passCut(true);

  if(p4_.Pt() < 20.0){passCut = false; return passCut;} // suggestion from convenors
  if(fabs(p4_.Eta())>2.40){passCut = false; return passCut;}
  //if(!isMedium()){passCut = false; return passCut;}
  if(!isLoose()){passCut = false; return passCut;}
  if(getD0() > 0.02 || getDz() > 0.1){passCut = false; return passCut;}
  if(isMedium() && getMiniIso() <= 0.2){passCut = false; return passCut;}
  else passCut = true;

  return passCut;
}
 
