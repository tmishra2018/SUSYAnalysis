#include "../include/analysis_photon.h"

bool recoPhoton::fireSingleTrg(int itrg) {
            if(((phoFiredSingleTrgs_ >> itrg)&1) == 1)return true;
            else{ return false;}
}

bool recoPhoton::fireDoubleTrg(int itrg) {
     if(((phoFiredDoubleTrgs_ >> itrg)&1) == 1)return true;
     else{ return false;}
}

bool recoPhoton::fireL1Trg(int itrg) {
     if(((phoFiredL1Trgs_ >> itrg)&1) == 1)return true;
     else{ return false;}
}

bool recoPhoton::isLoose() {
     if(((phoIDbit_ >> 0)&1)==1)return true;
     else return false;
}

bool recoPhoton::isMedium() {
     if(((phoIDbit_ >> 1)&1)==1)return true;
     else return false;
}

bool recoPhoton::isTight() {
     if(((phoIDbit_ >> 2)&1)==1)return true;
     else return false;
}

bool recoPhoton::isEB(){
     if(fabs(p4_.Eta())<1.444)return true;
     else return false;
}

bool recoPhoton::isEE(){
     if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5)return true;
     else return false;
}

bool recoPhoton::passHoverE(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
     	if(WP==1)pass = phoHoverE_ < 0.04596? true : false;
        else if(WP==2)pass = phoHoverE_ < 0.02197? true : false;
        else if(WP==3)pass = phoHoverE_ < 0.02148? true : false;
     }	
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
     	if(WP==1)pass = phoHoverE_ < 0.0590? true : false;
        else if(WP==2)pass = phoHoverE_ < 0.0326? true : false;
        else if(WP==3)pass = phoHoverE_ < 0.0321? true : false;
     }	
     return pass;
}

bool recoPhoton::passSigma(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.0106? true:false;
       else if(WP ==2)pass = phoSigmaIEtaIEta_ < 0.01015? true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ < 0.00996? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.0272? true:false;
       else if(WP==2)pass = phoSigmaIEtaIEta_ < 0.0272? true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ < 0.0271? true:false;
     } 
     return pass;
}

// WP = 1, for Loose photon
bool recoPhoton::passChIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFChIso_ < 1.694? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 1.141? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.65? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFChIso_ < 2.089? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 1.051? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.517? true:false;
     } 
     return pass;
}

// Mengleis
bool recoPhoton::passSigmaOLD(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.01031? true:false;
       else if(WP ==2)pass = phoSigmaIEtaIEta_ < 0.01022?true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ < 0.00994? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.03013? true:false;
       else if(WP==2)pass = phoSigmaIEtaIEta_ < 0.03001? true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ <0.03000? true:false;
     }
     return pass;
}


bool recoPhoton::passChIsoOLD(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFChIso_ < 1.295? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 0.441? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.202? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFChIso_ < 1.011? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 0.442? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.034? true:false;
     }
     return pass;
}

bool recoPhoton::passPhoIsoOLD(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFPhoIso_ < 3.630+0.0047*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 2.571+0.0047*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 2.362+0.0047*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFPhoIso_ < 6.641+0.0034*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 3.863+0.0034*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 2.617+0.0034*p4_.Et()? true:false;
     }
     return pass;
}


bool recoPhoton::passNeuIsoOLD(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFNeuIso_ < 10.910+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 2.725+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 0.264+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFNeuIso_ < 5.931+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 1.715+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 0.586+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
     }
     return pass;
}

bool recoPhoton::passHoverEOLD(int WP){
     bool pass(false);
     if(WP<=3)pass = phoHoverE_ < 0.05? true : false;
     return pass;
}
// mengleis finished

bool recoPhoton::passNeuIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFNeuIso_ <  24.032+0.01512*p4_.Et()+0.00002259*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 1.189+0.01512*p4_.Et()+0.00002259*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 0.317+0.1512*p4_.Et()+0.00002259*p4_.Et()*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFNeuIso_ < 19.722 + 0.0117*p4_.Et()+0.000023*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 2.718 + 0.0117*p4_.Et()+0.000023*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 2.716 + 0.0117*p4_.Et()+0.000023*p4_.Et()*p4_.Et()? true:false;
     } 
     return pass;
}

bool recoPhoton::passPhoIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFPhoIso_ < 2.876 + 0.004017*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 2.08 + 0.004017*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 2.044 + 0.004017*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFPhoIso_ < 4.162 + 0.0037*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 3.867 + 0.0037*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 3.032 + 0.0037*p4_.Et()? true:false;
     } 
     return pass;
}



bool recoPhoton::passSignalSelection(){
  
  bool passCut(true);
  
  if(Calibp4_.Et() < 35.0){passCut = false; return passCut;} 
  if(!isLoose()){passCut = false; return passCut;}
  else passCut = true;

  return passCut; 
}

bool recoPhoton::passBasicSelection(){
  bool passCut = (passHoverE(1) && passNeuIso(1) && passPhoIso(1));
	if(getChIso() > 20)passCut = false;
  return passCut;
}

bool recoPhoton::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;

    case DoubleEG2016: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break; // filters for Seeded leg
    case DoubleEG2017: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break;
    case DoubleEG2018: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break;
    case MCDoubleEG2016: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break; 
    case MCDoubleEG2017: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break;
    case MCDoubleEG2018: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break;

    // hltMu17Photon30CaloIdLL1ISOHEFilter || hltMu17Photon30CaloIdLL1ISOORHEFilter || hltMu38NoFiltersNoVtxPhoton38CaloIdLHEFilter
    case MuonEG2016: if(fireDoubleTrg(28)|| fireDoubleTrg(29) || fireDoubleTrg(30))passHLT = true; break;
    // hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter || hltMu43NoFiltersNoVtxPhoton43CaloIdLHEFilter
    case MuonEG2017:   if(fireDoubleTrg(27)|| fireDoubleTrg(36))passHLT = true; break;
    case MuonEG2018:   if(fireDoubleTrg(27)|| fireDoubleTrg(36))passHLT = true; break;
    
    case MCMuonEG2016:   if(fireDoubleTrg(28)|| fireDoubleTrg(29) || fireDoubleTrg(30))passHLT = true; break;
    case MCMuonEG2017:   if(fireDoubleTrg(27)|| fireDoubleTrg(36))passHLT = true; break;
    case MCMuonEG2018:   if(fireDoubleTrg(27)|| fireDoubleTrg(36))passHLT = true; break;
    
    case SingleElectron2016: passHLT = true; break;
    case SingleElectron2017: passHLT = true; break;
    case SingleElectron2018: passHLT = true; break;
    case SingleMuon2016: passHLT = true; break;
    case DoubleMuon2016: passHLT = true; break;
    case MET2016:  passHLT = true; break;
    case MCSingleElectron2016: passHLT = true; break;
    case MCSingleMuon2016: passHLT = true; break;
    case MCDoubleMuon2016: passHLT = true; break;
    case MCMET2016:  passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}
