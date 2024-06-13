#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "analysis_rawData.h"

void analysis_SUSY(){//main  


  TChain* es = new TChain("ggNtuplizer/EventTree");
  //es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-TChiWG_TuneCUETP8M1_RunIISummer16MiniAODv2.root");
  es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-T5Wg_TuneCUETP8M1_RunIISummer16MiniAODv2.root");

	es->SetBranchAddress("nMC",          &nMC);
	es->SetBranchAddress("mcPID",        &mcPID);
	es->SetBranchAddress("mcPt",         &mcPt);
	es->SetBranchAddress("mcMass",       &mcMass);
	es->SetBranchAddress("mcEta",        &mcEta);
	es->SetBranchAddress("mcPhi",        &mcPhi);
	es->SetBranchAddress("mcMomPID",     &mcMomPID);
	es->SetBranchAddress("mcMomPt",      &mcMomPt);
	es->SetBranchAddress("mcMomMass",    &mcMomMass);
	es->SetBranchAddress("mcMomEta",     &mcMomEta);
	es->SetBranchAddress("mcMomPhi",     &mcMomPhi);

  es->SetBranchAddress("rho",                 &rho);
  es->SetBranchAddress("pfMET",    			&pfMET);
  es->SetBranchAddress("pfMETPhi",    		&pfMETPhi);
	es->SetBranchAddress("pfMET_T1JERUp",			&pfMET_T1JERUp);
	es->SetBranchAddress("pfMET_T1JERDo",     &pfMET_T1JERDo);
	es->SetBranchAddress("pfMET_T1JESUp",     &pfMET_T1JESUp);
	es->SetBranchAddress("pfMET_T1JESDo",     &pfMET_T1JESDo);
	es->SetBranchAddress("pfMETPhi_T1JESUp",  &pfMETPhi_T1JESUp);
	es->SetBranchAddress("pfMETPhi_T1JESDo",  &pfMETPhi_T1JESDo);
	es->SetBranchAddress("pfMETPhi_T1UESUp",  &pfMETPhi_T1UESUp);
	es->SetBranchAddress("pfMETPhi_T1UESDo",  &pfMETPhi_T1UESDo);
  es->SetBranchAddress("metFilters",   		&metFilters);

  es->SetBranchAddress("nPho",                    &nPho);
  es->SetBranchAddress("phoCalibEt",              &phoCalibEt);
  es->SetBranchAddress("phoEta",                  &phoEta);
  es->SetBranchAddress("phoPhi",                  &phoPhi);
  es->SetBranchAddress("phohasPixelSeed",         &phohasPixelSeed);
  es->SetBranchAddress("phoR9Full5x5",                   &phoR9);
  es->SetBranchAddress("phoIDbit",      &phoIDbit);

  es->SetBranchAddress("nMu",           &nMu);
  es->SetBranchAddress("muPt",          &muPt);
  es->SetBranchAddress("muEta",         &muEta);
  es->SetBranchAddress("muPhi",         &muPhi);
  es->SetBranchAddress("muIDbit",       &muIDbit); 
	es->SetBranchAddress("muD0",          &muD0);
	es->SetBranchAddress("muDz",          &muDz);
  es->SetBranchAddress("muPFMiniIso",   &muPFMiniIso);

  es->SetBranchAddress("nEle",                    &nEle);
  es->SetBranchAddress("eleEn",                   &eleEn);
  es->SetBranchAddress("eleCalibEn",              &eleCalibEn);
  es->SetBranchAddress("eleD0",                   &eleD0);
  es->SetBranchAddress("eleDz",                   &eleDz);
  es->SetBranchAddress("elePt",                   &elePt);
  es->SetBranchAddress("eleCalibPt",              &eleCalibPt);
  es->SetBranchAddress("eleEta",                  &eleEta);
  es->SetBranchAddress("eleSCEta",                &eleSCEta);
  es->SetBranchAddress("elePhi",                  &elePhi);
  es->SetBranchAddress("eleR9",                   &eleR9);
  es->SetBranchAddress("eleHoverE",               &eleHoverE);
  es->SetBranchAddress("eleEoverPInv",            &eleEoverPInv);
  es->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEta);
  es->SetBranchAddress("eleConvVeto",             &eleConvVeto);
  es->SetBranchAddress("eleMissHits",             &eleMissHits);
  es->SetBranchAddress("elePFChIso",              &elePFChIso);
  es->SetBranchAddress("elePFPhoIso",             &elePFPhoIso);
  es->SetBranchAddress("elePFNeuIso",             &elePFNeuIso);
  es->SetBranchAddress("elePFPUIso",              &elePFPUIso);
  es->SetBranchAddress("elePFMiniIso",            &elePFMiniIso);
  es->SetBranchAddress("eleTrkdxy",               &eleTrkdxy);
  es->SetBranchAddress("eledEtaAtVtx",            &eledEtaAtVtx);
  es->SetBranchAddress("eledPhiAtVtx",            &eledPhiAtVtx);
  es->SetBranchAddress("eleIDbit",  &eleIDbit);

  es->SetBranchAddress("nJet",      &nJet);
  es->SetBranchAddress("jetPt",     &jetPt);
  es->SetBranchAddress("jetEta",    &jetEta);
  es->SetBranchAddress("jetPhi",    &jetPhi);
	es->SetBranchAddress("jetJECUnc", &jetJECUnc);


  std::ostringstream outputname;
  //outputname << "/uscms_data/d3/mengleis/test/resTree_TChiWG.root";
  outputname << "/uscms_data/d3/mengleis/test/resTree_T5WG_standalone.root";
  TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
  outputfile->cd();

  TTree *mctree = new TTree("SUSYtree","SUSYtree");
	float Mgluino(0);
  float Mchargino(0);
  float Mneutralino(0);

  mctree->Branch("Mgluino",        &Mgluino);
  mctree->Branch("Mchargino",      &Mchargino);
  mctree->Branch("Mneutralino",    &Mneutralino);

//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_lepPt(0);
  float eg_lepEta(0);
  float eg_lepPhi(0);
  float eg_sigMT(0);
  float eg_sigMET(0);
  float eg_sigMETPhi(0);
  float eg_dPhiLepMET(0);
  int   eg_nVertex(0);
  float eg_dRPhoLep(0);
  float eg_HT(0);
  float eg_nJet(0);
  float eg_invmass(0);
	float eg_sigMETJESup(0);
	float eg_sigMETJESdo(0);
	float eg_sigMETJERup(0);
	float eg_sigMETJERdo(0);
	float eg_sigMTJESup(0);
	float eg_sigMTJESdo(0);
	float eg_sigMTJERup(0);
	float eg_sigMTJERdo(0);
	float eg_HTJESup(0);
	float eg_HTJESdo(0);
	float eg_dPhiLepMETJESup(0);
	float eg_dPhiLepMETJESdo(0);
	float eg_dPhiLepMETJERup(0);
	float eg_dPhiLepMETJERdo(0);

  egtree->Branch("Mgluino",    &Mgluino);
  egtree->Branch("Mchargino",  &Mchargino);
  egtree->Branch("Mneutralino",&Mneutralino);
  egtree->Branch("phoEt",      &eg_phoEt);
  egtree->Branch("phoEta",     &eg_phoEta);
  egtree->Branch("phoPhi",     &eg_phoPhi);
  egtree->Branch("lepPt",      &eg_lepPt);
  egtree->Branch("lepEta",     &eg_lepEta);
  egtree->Branch("lepPhi",     &eg_lepPhi);
  egtree->Branch("sigMT",      &eg_sigMT);
  egtree->Branch("sigMET",     &eg_sigMET);
  egtree->Branch("sigMETPhi",  &eg_sigMETPhi);
  egtree->Branch("dPhiLepMET", &eg_dPhiLepMET);
  egtree->Branch("nVertex",    &eg_nVertex);
  egtree->Branch("dRPhoLep",   &eg_dRPhoLep);
  egtree->Branch("HT",         &eg_HT);
  egtree->Branch("nJet",       &eg_nJet);
  egtree->Branch("invmass",    &eg_invmass);
	egtree->Branch("sigMETJESup",&eg_sigMETJESup);
	egtree->Branch("sigMETJESdo",&eg_sigMETJESdo);
	egtree->Branch("sigMETJERup",&eg_sigMETJERup);
	egtree->Branch("sigMETJERdo",&eg_sigMETJERdo);
	egtree->Branch("sigMTJESup", &eg_sigMTJESup);
	egtree->Branch("sigMTJESdo", &eg_sigMTJESdo);
	egtree->Branch("sigMTJERup", &eg_sigMTJERup);
	egtree->Branch("sigMTJERdo", &eg_sigMTJERdo);
	egtree->Branch("dPhiLepMETJESup", &eg_dPhiLepMETJESup);
	egtree->Branch("dPhiLepMETJESdo", &eg_dPhiLepMETJESdo);
	egtree->Branch("dPhiLepMETJERup", &eg_dPhiLepMETJERup);
	egtree->Branch("dPhiLepMETJERdo", &eg_dPhiLepMETJERdo);
	egtree->Branch("HTJESup",     &eg_HTJESup);
	egtree->Branch("HTJESdo",     &eg_HTJESdo);

  TTree *mgtree = new TTree("mgTree","mgTree");
  float mg_phoEt(0);
  float mg_phoEta(0);
  float mg_phoPhi(0);
  float mg_lepPt(0);
  float mg_lepEta(0);
  float mg_lepPhi(0);
  float mg_sigMT(0);
  float mg_sigMET(0);
  float mg_sigMETPhi(0);
  float mg_dPhiLepMET(0);
	float mg_threeMass(0);
  int   mg_nVertex(0);
  float mg_dRPhoLep(0);
  float mg_HT(0);
  float mg_nJet(0);
	float mg_sigMETJESup(0);
	float mg_sigMETJESdo(0);
	float mg_sigMETJERup(0);
	float mg_sigMETJERdo(0);
	float mg_sigMTJESup(0);
	float mg_sigMTJESdo(0);
	float mg_sigMTJERup(0);
	float mg_sigMTJERdo(0);
	float mg_HTJESup(0);
	float mg_HTJESdo(0);
	float mg_dPhiLepMETJESup(0);
	float mg_dPhiLepMETJESdo(0);
	float mg_dPhiLepMETJERup(0);
	float mg_dPhiLepMETJERdo(0);

  mgtree->Branch("Mgluino",    &Mgluino);
  mgtree->Branch("Mchargino",  &Mchargino);
  mgtree->Branch("Mneutralino",&Mneutralino);
  mgtree->Branch("phoEt",      &mg_phoEt);
  mgtree->Branch("phoEta",     &mg_phoEta);
  mgtree->Branch("phoPhi",     &mg_phoPhi);
  mgtree->Branch("lepPt",      &mg_lepPt);
  mgtree->Branch("lepEta",     &mg_lepEta);
  mgtree->Branch("lepPhi",     &mg_lepPhi);
  mgtree->Branch("sigMT",      &mg_sigMT);
  mgtree->Branch("sigMET",     &mg_sigMET);
  mgtree->Branch("sigMETPhi",  &mg_sigMETPhi);
  mgtree->Branch("dPhiLepMET", &mg_dPhiLepMET);
	mgtree->Branch("threeMass",  &mg_threeMass);
  mgtree->Branch("nVertex",    &mg_nVertex);
  mgtree->Branch("dRPhoLep",   &mg_dRPhoLep);
  mgtree->Branch("HT",         &mg_HT);
  mgtree->Branch("nJet",       &mg_nJet);
	mgtree->Branch("sigMETJESup",&mg_sigMETJESup);
	mgtree->Branch("sigMETJESdo",&mg_sigMETJESdo);
	mgtree->Branch("sigMETJERup",&mg_sigMETJERup);
	mgtree->Branch("sigMETJERdo",&mg_sigMETJERdo);
	mgtree->Branch("sigMTJESup", &mg_sigMTJESup);
	mgtree->Branch("sigMTJESdo", &mg_sigMTJESdo);
	mgtree->Branch("sigMTJERup", &mg_sigMTJERup);
	mgtree->Branch("sigMTJERdo", &mg_sigMTJERdo);
	mgtree->Branch("dPhiLepMETJESup", &mg_dPhiLepMETJESup);
	mgtree->Branch("dPhiLepMETJESdo", &mg_dPhiLepMETJESdo);
	mgtree->Branch("dPhiLepMETJERup", &mg_dPhiLepMETJERup);
	mgtree->Branch("dPhiLepMETJERdo", &mg_dPhiLepMETJERdo);
	mgtree->Branch("HTJESup",     &mg_HTJESup);
	mgtree->Branch("HTJESdo",     &mg_HTJESdo);

  float MET(0);
  float METPhi(0);
	float MET_T1JERUp(0);
	float MET_T1JERDo(0);
	float MET_T1JESUp(0);
	float MET_T1JESDo(0);	
	float	METPhi_T1JESUp(0);
	float	METPhi_T1JESDo(0);
	float	METPhi_T1UESUp(0);
	float	METPhi_T1UESDo(0);
  int jetNumber(0);

  const int nEvts = es->GetEntries(); 
  std::cout << "total event : " << nEvts << std::endl;

  for (int ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		es->GetEntry(ievt);
    if(ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

		Mgluino=-1;
  	Mchargino=-1;
  	Mneutralino=-1;

    MET = pfMET;
    METPhi = pfMETPhi;
		MET_T1JERUp = pfMET_T1JERUp;
		MET_T1JERDo = pfMET_T1JERDo;
		MET_T1JESUp = pfMET_T1JESUp;
		MET_T1JESDo = pfMET_T1JESDo;
		METPhi_T1JESUp = pfMETPhi_T1JESUp;
		METPhi_T1JESDo = pfMETPhi_T1JESDo;
		METPhi_T1UESUp = pfMETPhi_T1UESUp;
		METPhi_T1UESDo = pfMETPhi_T1UESDo;


		/*********************  find gluino mass and fill the MC tree ********************************/
		Mgluino=-1;
  	Mchargino=-1;
  	Mneutralino=-1;

    for(int iMC(0); iMC < nMC;  iMC++){

        //Look for gluino
        if( (*mcMomPID)[iMC] == 1000021)Mgluino = (*mcMomMass)[iMC]; 
        //Look for neutralino
        if( (*mcPID)[iMC] == 1000022 && (*mcMomPID)[iMC] == 1000023)Mneutralino = (*mcMomMass)[iMC]; 
        //Look for chargino
        if( (*mcPID)[iMC] == 1000022 && fabs((*mcMomPID)[iMC]) == 1000024)Mchargino = (*mcMomMass)[iMC];
         
    }//loop on MC particles

    mctree->Fill();


		/******************** select signal photons ****************************************************/
		int  phoIndex(-1);
		for(int iPho(0); iPho < nPho; iPho++){
			if( (*phoR9)[iPho] < 0.5)continue;
  		if( (*phoCalibEt)[iPho] < 35.0 || fabs((*phoEta)[iPho]) > 1.4442 )continue;
      if( (((*phoIDbit)[iPho] >> 0)&1)!=1)continue;
			bool PixelVeto = ((*phohasPixelSeed)[iPho]==0)? true: false;
			bool GSFveto(true);
			bool FSRVeto(true);
			for(int ie(0); ie < nEle; ie++){
			 if(DeltaR((*phoEta)[iPho], (*phoPhi)[iPho],(*eleEta)[ie],(*elePhi)[ie]) < 0.02)GSFveto = false;
			 if(DeltaR((*phoEta)[iPho], (*phoPhi)[iPho],(*eleEta)[ie],(*elePhi)[ie]) < 0.3)FSRVeto=false;
			}
			for(int im(0); im < nMu; im++)
			 if(DeltaR((*phoEta)[iPho], (*phoPhi)[iPho],(*muEta)[im],(*muPhi)[im]) < 0.3)FSRVeto=false;
			if(GSFveto && PixelVeto && FSRVeto){
				if(phoIndex < 0){ // select the leading one
					phoIndex = iPho;
				}
			}
		}

		/******************** select signal electron ****************************************************/
		int eleIndex(-1);
		if(phoIndex >= 0){
			for(int ie(0); ie < nEle; ie++){
				if(eleIndex >= 0)continue; //already find an electron
				if( (*elePt)[ie] < 25 || fabs((*eleSCEta)[ie]) > 2.5)continue;
				if(( fabs((*eleSCEta)[ie])<1.4442 && (*eleR9)[ie] < 0.5) || (fabs((*eleSCEta)[ie])> 1.56 && (*eleR9)[ie] < 0.8))continue;
				// electron ID
				// bool passEleSignalSelection(float eta, float sigma, float dEta, float dPhi, float HoverE, float EoverP, int MissHit, int ConvVeto, float miniIso )
				// defined in analysis_raw.h
				if( passEleSignalSelection((*eleSCEta)[ie],  (*eleSigmaIEtaIEta)[ie], (*eledEtaAtVtx)[ie], (*eledPhiAtVtx)[ie], (*eleHoverE)[ie], (*eleEoverPInv)[ie], (*eleMissHits)[ie], (*eleConvVeto)[ie], (*elePFMiniIso)[ie] )){
				  eleIndex = ie;	
				}
			}
		}

		/******************** select signal muon *******************************************************/
		int muIndex(-1);
		if(phoIndex >= 0){
			for(int im(0); im < nMu; im++){
				if( muIndex >=0 )continue; //already find an muon
  			if( (*muPt)[im] < 25.0)continue;
  			if(fabs( (*muEta)[im] )>2.40)continue;
				if(fabs( (*muD0)[im] ) > 0.05 || fabs( (*muDz)[im] ) > 0.1)continue;
  			if((( (*muIDbit)[im] >> 1)&1) == 1 && (*muPFMiniIso)[im] < 0.2){
					muIndex = im;
				}
			}
		}


		/******************  Fill e+g tree ************************************************************/
		if(phoIndex >= 0 && eleIndex >=0){
			double dReg = DeltaR( (*phoEta)[phoIndex], (*phoPhi)[phoIndex], (*eleSCEta)[eleIndex], (*elePhi)[eleIndex]);
			if(dReg>0.8){
    			TLorentzVector phoP4(0,0,0,0);
					TLorentzVector eleP4(0,0,0,0);
    			phoP4.SetPtEtaPhiM((*phoCalibEt)[phoIndex],(*phoEta)[phoIndex],(*phoPhi)[phoIndex],0);
    			eleP4.SetPtEtaPhiE((*eleCalibPt)[eleIndex],(*eleSCEta)[eleIndex],(*elePhi)[eleIndex],(*eleCalibEn)[eleIndex]);
					if(fabs((phoP4 + eleP4).M() - 91.188) > 10.0){

						float deltaPhi = DeltaPhi( (*elePhi)[eleIndex], METPhi);
						float MT = sqrt(2*MET*(*eleCalibPt)[eleIndex]*(1-std::cos(deltaPhi)));
						eg_phoEt = (*phoCalibEt)[phoIndex];
						eg_phoEta= (*phoEta)[phoIndex];
						eg_phoPhi= (*phoPhi)[phoIndex];
						eg_lepPt = (*eleCalibPt)[eleIndex];
						eg_lepEta= (*eleSCEta)[eleIndex];
						eg_lepPhi= (*elePhi)[eleIndex];
						eg_sigMT = MT;
						eg_sigMET= MET;
						eg_sigMETPhi = METPhi;
						eg_dPhiLepMET = deltaPhi; 
						eg_nVertex = nVtx; 
						eg_dRPhoLep= dReg;
						eg_invmass = (phoP4 + eleP4).M();
						eg_sigMETJESup = MET_T1JESUp;
						eg_sigMETJESdo = MET_T1JESDo;
						eg_sigMETJERup = MET_T1JERUp;
						eg_sigMETJERdo = MET_T1JERDo;
						eg_dPhiLepMETJESup = DeltaPhi((*elePhi)[eleIndex], METPhi_T1JESUp);
						eg_dPhiLepMETJESdo = DeltaPhi((*elePhi)[eleIndex], METPhi_T1JESDo);
						eg_dPhiLepMETJERup = deltaPhi;
						eg_dPhiLepMETJERdo = deltaPhi;
						eg_sigMTJESup = sqrt(2*MET_T1JESUp*(*eleCalibPt)[eleIndex]*(1-std::cos(eg_dPhiLepMETJESup)));
						eg_sigMTJESdo = sqrt(2*MET_T1JESDo*(*eleCalibPt)[eleIndex]*(1-std::cos(eg_dPhiLepMETJESdo)));
						eg_sigMTJERup = sqrt(2*MET_T1JERUp*(*eleCalibPt)[eleIndex]*(1-std::cos(eg_dPhiLepMETJERup)));
						eg_sigMTJERdo = sqrt(2*MET_T1JERDo*(*eleCalibPt)[eleIndex]*(1-std::cos(eg_dPhiLepMETJERdo)));
					
						eg_nJet = 0;
						eg_HT = 0;
						eg_HTJESup = 0;
						eg_HTJESdo = 0;
						for(int iJ(0); iJ < nJet; iJ++){
							if( (*jetPt)[iJ] < 30 || fabs( (*jetEta)[iJ]) > 2.5)continue;
							if(DeltaR((*jetEta)[iJ], (*jetPhi)[iJ],(*eleEta)[eleIndex],(*elePhi)[eleIndex]) <= 0.4)continue;	
							if(DeltaR((*jetEta)[iJ], (*jetPhi)[iJ],(*phoEta)[phoIndex],(*phoPhi)[phoIndex]) <= 0.4)continue;
							eg_nJet += 1;
							eg_HT += (*jetPt)[iJ]; 
							eg_HTJESup += (*jetPt)[iJ]*(1+ (*jetJECUnc)[iJ]);
							eg_HTJESdo += (*jetPt)[iJ]*(1- (*jetJECUnc)[iJ]);
						}	

						egtree->Fill();

				}// Z mass Filter
			}//dR filter
		}// ele + pho candidate
		 

		 if(phoIndex >= 0 && muIndex >=0){
				double dRmg = DeltaR((*phoEta)[phoIndex], (*phoPhi)[phoIndex], (*muEta)[muIndex], (*muPhi)[muIndex]);
				if(dRmg>0.8){
					float deltaPhi = DeltaPhi((*muPhi)[muIndex], METPhi);
					float MT = sqrt(2*MET*(*muPt)[muIndex]*(1-std::cos(deltaPhi)));

					mg_phoEt =(*phoCalibEt)[phoIndex]; 
					mg_phoEta=(*phoEta)[phoIndex];
					mg_phoPhi=(*phoPhi)[phoIndex];
					mg_lepPt =(*muPt)[muIndex]; 
					mg_lepEta=(*muEta)[muIndex];
					mg_lepPhi=(*muPhi)[muIndex];
					mg_sigMT = MT;
					mg_sigMET= MET;
					mg_sigMETPhi = METPhi;
					mg_dPhiLepMET = deltaPhi;
					mg_nVertex = nVtx;
					mg_dRPhoLep= dRmg;
					mg_sigMETJESup = MET_T1JESUp;
					mg_sigMETJESdo = MET_T1JESDo;
					mg_sigMETJERup = MET_T1JERUp;
					mg_sigMETJERdo = MET_T1JERDo;
					mg_dPhiLepMETJESup = DeltaPhi((*muPhi)[muIndex], METPhi_T1JESUp);
					mg_dPhiLepMETJESdo = DeltaPhi((*muPhi)[muIndex], METPhi_T1JESDo);
					mg_dPhiLepMETJERup = deltaPhi;
					mg_dPhiLepMETJERdo = deltaPhi;
					mg_sigMTJESup = sqrt(2*MET_T1JESUp*(*muPt)[muIndex]*(1-std::cos(mg_dPhiLepMETJESup)));
					mg_sigMTJESdo = sqrt(2*MET_T1JESDo*(*muPt)[muIndex]*(1-std::cos(mg_dPhiLepMETJESdo)));
					mg_sigMTJERup = sqrt(2*MET_T1JERUp*(*muPt)[muIndex]*(1-std::cos(mg_dPhiLepMETJERup)));
					mg_sigMTJERdo = sqrt(2*MET_T1JERDo*(*muPt)[muIndex]*(1-std::cos(mg_dPhiLepMETJERdo)));

					mg_nJet = 0;
					mg_HT = 0;
					mg_HTJESup = 0;
					mg_HTJESdo = 0;
					for(int iJ(0); iJ < nJet; iJ++){
						if( (*jetPt)[iJ] < 30 || fabs( (*jetEta)[iJ]) > 2.5)continue;
						if(DeltaR((*jetEta)[iJ], (*jetPhi)[iJ],(*muEta)[muIndex], (*muPhi)[muIndex]) <= 0.4)continue;	
						if(DeltaR((*jetEta)[iJ], (*jetPhi)[iJ],(*phoEta)[phoIndex],(*phoPhi)[phoIndex]) <= 0.4)continue;
						mg_nJet += 1;
						mg_HT += (*jetPt)[iJ];
						mg_HTJESup += (*jetPt)[iJ]*(1+ (*jetJECUnc)[iJ]);
						mg_HTJESdo += (*jetPt)[iJ]*(1- (*jetJECUnc)[iJ]);
					}	
					mgtree->Fill();
			 }//dR Filter
		 }//Candidate Filter

	}//loop on entries


	outputfile->Write();
}
