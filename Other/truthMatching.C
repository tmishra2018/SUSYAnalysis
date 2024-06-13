float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

		std::vector<int>  *mcPID=0;
		std::vector<int>  *mcMomPID=0;
		std::vector<int>  *mcGMomPID=0;
    std::vector<float> *mcPt=0;
    std::vector<float> *mcEta=0;
    std::vector<float> *mcPhi=0;
    
		tree->SetBranchAddress("mcPID",        &mcPID);
		tree->SetBranchAddress("mcPt",         &mcPt);
		tree->SetBranchAddress("mcEta",        &mcEta);
		tree->SetBranchAddress("mcPhi",        &mcPhi);
		tree->SetBranchAddress("mcGMomPID",    &mcGMomPID);
		tree->SetBranchAddress("mcMomPID",     &mcMomPID);

   // Assume the photon has eta=photonEta, phi=photonPhi, pt=photonPt

	 int matchIndex(-1);
	 for(int itMC(0); itMC < mcPID->size(); itMC++){
		 if( (*mcPt)[itMC] < 10.0)continue;
		 float dR = DeltaR( (*mcEta)[itMC], (*mcPhi)[itMC], photonEta, photonPhi);
		 float dE = fabs( (*mcPt)[itMC] - photonPt)/photonPt; 
		 if( dR < 0.1 && dE < 0.1 ){
			 matchIndex = itMC;
		 }
   }

	 bool isTruePhoton(false);
   if(matchIndex >= 0){
	   std::cout << " PID = " << (*mcPID)[matchIndex] 
               << " MomPID = " << (*mcMomPID)[matchIndex] << std::endl;

		 if((*mcPID)[matchIndex] == 22 && (fabs((*mcMomPID)[matchIndex]) == 11 || fabs((*mcMomPID)[matchIndex]) == 13))isTruePhoton=true;
   }
