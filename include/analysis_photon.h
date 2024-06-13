#ifndef ANALYSIS_PHOTON
#define ANALYSIS_PHOTON

#ifndef ANALYSIS_RAWDATA
#include "analysis_rawData.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif



class recoPhoton{
  public:
  recoPhoton(rawData& raw, int ipho){
    runtype_ = raw.runtype_;

    p4_.SetPtEtaPhiM((*raw.phoEt)[ipho],(*raw.phoEta)[ipho],(*raw.phoPhi)[ipho],0);
    Calibp4_.SetPtEtaPhiM((*raw.phoCalibEt)[ipho],(*raw.phoEta)[ipho],(*raw.phoPhi)[ipho],0);
    phohasPixelSeed_ = (*raw.phohasPixelSeed)[ipho];
    phoEleVeto_ = (*raw.phoEleVeto)[ipho];
    phoR9_ = (*raw.phoR9)[ipho];
    phoHoverE_ = (*raw.phoHoverE)[ipho];

    phoSigmaIEtaIEta_ = (*raw.phoSigmaIEtaIEtaFull5x5)[ipho];

    float eta = fabs((*raw.phoSCEta)[ipho]);
    float rho = raw.rho;

    if(eta<=1.0){
	  phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0112);
	  phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0668);
	  phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.1113);
    }
    else if(eta>1.0 && eta<= 1.479){
	  phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0108);
	  phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.1054);
	  phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.0953);
    }
    else if(eta>1.479 && eta<=2.0){
	  phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0106);
	  phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0786);
	  phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.0619);
    }
    else if(eta>2.0 && eta<=2.2){
	  phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.01002);
	  phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0233);
	  phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.0837);
    } 
    else if(eta>2.2 && eta<=2.3){
	  phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0098);
	  phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0078);
	  phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.1070);
    } 
    else if(eta>2.3 && eta<=2.4){
      phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0089);
      phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0028);
      phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.1212);
    } 
    else if(eta>2.4 && eta<=2.5){
      phoPFChIso_  = std::max(0.0, (*raw.phoPFChIso)[ipho] - rho*0.0087);
      phoPFNeuIso_ = std::max(0.0, (*raw.phoPFNeuIso)[ipho]- rho*0.0137);
      phoPFPhoIso_ = std::max(0.0, (*raw.phoPFPhoIso)[ipho]- rho*0.1466);
    } 
    else{
      phoPFChIso_  = (*raw.phoPFChIso)[ipho];
      phoPFPhoIso_ = (*raw.phoPFPhoIso)[ipho];
      phoPFNeuIso_ = (*raw.phoPFNeuIso)[ipho];
    } 
 	// Only comment these for signal files T5WG, as triggers are not stored in ntuples     
    phoFiredSingleTrgs_ = (*raw.phoFiredSingleTrgs)[ipho];
    phoFiredDoubleTrgs_ = (*raw.phoFiredDoubleTrgs)[ipho];
    phoFiredL1Trgs_ = (*raw.phoFiredL1Trgs)[ipho]; 
    phoIDbit_ = (*raw.phoIDbit)[ipho];
 }

 ~recoPhoton(){
 };
   
inline float getE() { return p4_.E();}
inline float getEt() { return p4_.Et();}
inline float getCalibEt() { return Calibp4_.Et();}
inline float getEta() { return p4_.Eta();}
inline float getPhi() {return p4_.Phi();}
inline TLorentzVector getP4() {
         TLorentzVector p = p4_;
          return p;
}
inline TLorentzVector getCalibP4() {
         TLorentzVector p = Calibp4_;
         return p;
}
inline int PixelSeed() {return phohasPixelSeed_;}
inline int EleVeto() {return phoEleVeto_;}
inline float getR9() {return phoR9_;}
inline float getHoverE() {return phoHoverE_;}
inline float getSigma() {return phoSigmaIEtaIEta_;}
inline float getChIso() {return phoPFChIso_;}
inline float getNeuIso() {return phoPFNeuIso_;}
inline float getPhoIso() {return phoPFPhoIso_;}
inline int   getDoubleTrg() { return phoFiredDoubleTrgs_;} 

bool fireSingleTrg(int itrg);
bool fireDoubleTrg(int itrg);
bool fireL1Trg(int itrg);
bool isLoose();
bool isMedium();
bool isTight();
bool isEB();
bool isEE();
bool passHoverE(int WP);
bool passSigma(int WP);
bool passChIso(int WP);
bool passNeuIso(int WP);
bool passPhoIso(int WP);
bool passSigmaOLD(int WP);
bool passChIsoOLD(int WP);
bool passHoverEOLD(int WP);
bool passNeuIsoOLD(int WP);
bool passPhoIsoOLD(int WP);

bool passSignalSelection();
bool passBasicSelection();
bool passHLTSelection();

private:
    RunType runtype_;

    TLorentzVector p4_;
    TLorentzVector Calibp4_;
    int    phohasPixelSeed_;
    int    phoEleVeto_;
    float  phoR9_;
    float  phoHoverE_;
    float  phoSigmaIEtaIEta_;
    float  phoPFChIso_;
    float  phoPFPhoIso_;
    float  phoPFNeuIso_;
    Int_t  phoFiredSingleTrgs_;
    Int_t  phoFiredDoubleTrgs_;
    Int_t  phoFiredL1Trgs_;
    UShort_t phoIDbit_;
};
   

#endif
