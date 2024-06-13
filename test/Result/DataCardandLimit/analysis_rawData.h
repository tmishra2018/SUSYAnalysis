Int_t     nVtx = 0;
float     rho =0;
float     pfMET = 0;
float     pfMETPhi = 0;
float     genMET = 0;
float     genMETPhi = 0;
float 		pfMET_T1JERUp = 0;
float 		pfMET_T1JERDo = 0;
float 		pfMET_T1JESUp = 0;
float 		pfMET_T1JESDo = 0;	
float			pfMETPhi_T1JESUp = 0;
float			pfMETPhi_T1JESDo = 0;
float			pfMETPhi_T1UESUp = 0;
float			pfMETPhi_T1UESDo = 0;
Int_t     metFilters=0;

Int_t                 nMC =0;
std::vector<int>      *mcPID =0;
std::vector<float>    *mcPt =0;
std::vector<float>    *mcMass =0;
std::vector<float>    *mcEta =0;
std::vector<float>    *mcPhi =0;
std::vector<float>    *mcE =0;
std::vector<float>    *mcEt =0;
vector<UShort_t>      *mcStatusFlag=0;
vector<int>      			*mcStatus=0;
std::vector<int>      *mcGMomPID =0;
std::vector<int>      *mcMomPID =0;
std::vector<float>    *mcMomPt =0;
std::vector<float>    *mcMomMass =0;
std::vector<float>    *mcMomEta =0;
std::vector<float>    *mcMomPhi =0;

Int_t          nPho =0;
vector<float>  *phoE = 0;
vector<float>  *phoEt = 0;
vector<float>  *phoCalibEt = 0;
vector<float>  *phoEta = 0;
vector<float>  *phoPhi = 0;
vector<float>  *phoSCEta =0;
vector<float>  *phoSCPhi =0;
vector<int>    *phohasPixelSeed =0;
vector<int>    *phoEleVeto =0;
vector<float>  *phoR9 =0;
vector<float>  *phoHoverE =0;
vector<float>  *phoSigmaIEtaIEta =0;
vector<float>  *phoSigmaIEtaIEtaFull5x5 = 0;
vector<float>  *phoPFChIso =0;
vector<float>  *phoPFPhoIso =0;
vector<float>  *phoPFNeuIso =0;
vector<Int_t>  *phoFiredSingleTrgs =0;
vector<Int_t>  *phoFiredDoubleTrgs =0;
vector<Int_t>  *phoFiredL1Trgs =0;
vector<UShort_t> *phoIDbit =0;


Int_t          nMu =0;
vector<float>  *muPt =0;
vector<float>  *muEn =0;
vector<float>  *muEta =0;
vector<float>  *muPhi =0;
vector<int>    *muCharge =0;
vector<int>    *muType =0;
vector<UShort_t> *muIDbit =0;
vector<float>  *muD0 = 0;
vector<float>  *muDz = 0;
vector<float>  *muPFMiniIso=0;
vector<float>  *muPFChIso =0;
vector<float>  *muPFPhoIso =0;
vector<float>  *muPFNeuIso =0;
vector<float>  *muPFPUIso =0;
vector<UInt_t> *muFiredTrgs =0;
vector<UInt_t> *muFiredL1Trgs = 0;

Int_t          nEle =0;
vector<int>    *eleCharge =0;
vector<float>  *eleEn =0;
vector<float>  *eleD0 =0;
vector<float>  *eleDz =0;
vector<float>  *elePt =0;
vector<float>  *eleCalibPt=0;
vector<float>  *eleCalibEn=0;
vector<float>  *eleEta =0;
vector<float>  *eleSCEta =0;
vector<float>  *elePhi =0;
vector<float>  *eleR9 =0;
vector<float>  *eleHoverE =0;
vector<float>  *eleEoverPInv =0;
vector<float>  *eleSigmaIEtaIEta =0;
vector<int>    *eleConvVeto =0;
vector<int>    *eleMissHits =0;
vector<float>  *elePFChIso =0;
vector<float>  *elePFPhoIso =0;
vector<float>  *elePFNeuIso =0;
vector<float>  *elePFPUIso =0;
vector<float>  *elePFMiniIso =0;
vector<float>  *eleTrkdxy =0;
vector<float>  *eledEtaAtVtx=0;
vector<float>  *eledPhiAtVtx=0;
vector<UInt_t>  *eleFiredSingleTrgs=0;
vector<UInt_t>  *eleFiredDoubleTrgs=0;
vector<UInt_t>  *eleFiredL1Trgs=0;
vector<UShort_t> *eleIDbit =0;

Int_t          nJet = 0;
std::vector<float>  *jetPt = 0;
std::vector<float>  *jetEn = 0;
std::vector<float>  *jetEta = 0;
std::vector<float>  *jetPhi = 0;
std::vector<float>  *jetArea = 0;
std::vector<float>  *jetJECUnc = 0;
std::vector<float>  *jetCSV2BJetTags = 0;


bool passEleSignalSelection(float eta, float sigma, float dEta, float dPhi, float HoverE, float EoverP, int MissHit, int ConvVeto, float miniIso ){
  bool passMiniMedium(true);
 
  if( fabs(eta) > 1.4442 && fabs(eta) < 1.56) {passMiniMedium=false; return passMiniMedium; }
  if(fabs(eta) < 1.4442){
		if(sigma > 0.00998) {passMiniMedium=false; return passMiniMedium; }
		if(fabs(dEta) > 0.00311){passMiniMedium=false; return passMiniMedium;}
		if(fabs(dPhi) > 0.103){passMiniMedium=false; return passMiniMedium;}
		if(HoverE > 0.253){passMiniMedium=false; return passMiniMedium;}
		if(fabs(EoverP) > 0.134){passMiniMedium=false; return passMiniMedium;}
		if(MissHit > 1){passMiniMedium=false; return passMiniMedium;}
		if(ConvVeto == 0){passMiniMedium=false; return passMiniMedium;}
		if(miniIso > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  else if(fabs(eta) > 1.56){
		if(sigma > 0.0298) {passMiniMedium=false; return passMiniMedium; }
		if(fabs(dEta) > 0.00609){passMiniMedium=false; return passMiniMedium;}
		if(fabs(dPhi) > 0.045){passMiniMedium=false; return passMiniMedium;}
		if(HoverE > 0.0878){passMiniMedium=false; return passMiniMedium;}
		if(fabs(EoverP) > 0.13){passMiniMedium=false; return passMiniMedium;}
		if(MissHit > 1){passMiniMedium=false; return passMiniMedium;}
		if(ConvVeto == 0){passMiniMedium=false; return passMiniMedium;}
		if(miniIso  > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  return passMiniMedium;
}

float DeltaPhi(float phi1, float phi2){
  float deltaPhi = phi1 - phi2;
  if(fabs(deltaPhi) > TMath::Pi()){
    if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
    else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
  }
  return deltaPhi;
}

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
