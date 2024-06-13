// g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_SUSY.C -o analysis_SUSY.exe `root-config --libs`
#include "../include/analysis_commoncode.h"
#include "../include/analysis_cuts.h"

#define debugMCLevel 4
enum branchType{
	gg = 1,
	gZ = 2,
	gH = 3,
	gW = 4,
	ZZ = 5,
	HH = 6,
	noType = 7	
};

struct decayChain{
	std::vector<mcData>::iterator iter; 
  std::vector< std::vector<mcData>::iterator > daughter;
};

void analysis_SUSY_modified_tribeni(int RunYear, const char *channel){//main

	//gSystem->Load("../lib/libAnaClasses.so");

	TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-T5WG_TuneCUETP8M1_RunIISummer16MiniAOD.root");
	//es->Add(Form("/eos/uscms/store/group/lpcsusyhad/Tribeni/SMS-TChiWg_mChi-1000_mLSP-1/SMS-TChiWg_mChi-1000_mLSP-1_%d.root",RunYear));
	//es->Add(Form("/eos/uscms/store/group/lpcsusyhad/Tribeni/T5Wg_mG-1800_mLSP-800/T5Wg_mG-1800_mLSP-800_%d.root",RunYear));

	//es->Add("/eos/uscms/store/group/lpcsusyhad/Tribeni/Signal/TChiWG_FastSim_2018.root");
	//es->Add("/eos/uscms/store/user/mengleis/Signal/SMS-T5WG_TuneCUETP8M1_RunIISummer16MiniAOD.root");
	//es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-T5Wg_TuneCUETP8M1_RunIISummer16MiniAODv2_scan.root");
	//es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-TChiNG_BF50N50G_TuneCUETP8M1.root");
	//es->Add("/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/SMS-T6Wg_TuneCUETP8M1_RunIISummer16MiniAODv2_scan.root");

	RunType datatype(MC); 
	//outputname << "/eos/uscms/store/user/tmishra/Signal/resTree_TChiWG_2018.root";
	//outputname << "/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_test.root";
	//outputname << "/uscms_data/d3/mengleis/FullStatusOct/resTree_TChiNg_test_debug.root";
	//outputname << "/uscms_data/d3/mengleis/FullStatusOct/resTree_GMSB_test.root";
	//outputname << "/uscms_data/d3/mengleis/FullStatusOct/resTree_GMSB_M1M3.root";
	//outputname << "/uscms_data/d3/mengleis/FullStatusOct/resTree_T6WG.root";

	int SUSYtype(-1);
	std::string model = channel;
	if(model == "T5WG"){
		std::cout << "T5WG Model !" << std::endl;
		SUSYtype = 5;
	}
	else if(model == "T6WG"){
		std::cout << "T6WG Model !" << std::endl;
		SUSYtype = 6;
	}
	else if(model == "TChiWG"){
		std::cout << "TChiWG Model !" << std::endl;
		SUSYtype = 1;
	}
	else if(model == "GMSB"){
		std::cout << "GMSB Model !" << std::endl;
		SUSYtype = 2;
	}
	else{
		std::cout << "model undefined" << std::endl;
		abort();
	}

	int decaybranch;
        //TFile *output = TFile::Open(Form("/eos/uscms/store/user/tmishra/Signal/resTree_%s_%d.root",channel,RunYear),"RECREATE");
        TFile *output = TFile::Open("resTree_T5WG_SUSY_modified_.root","RECREATE");
	output->cd();
	TTree *tree = new TTree("SUSYtree","SUSYtree");
	float Mass1_(0);
	float Mass2_(0);
	float mcPhotonEt_(0);
	float mcPhotonEta_(0);
	float mcPhotonPhi_(0);
	float recoPhotonEt_(0);
	float recoPhotonEta_(0);
	float recoPhotonPhi_(0);
	float PhoR9_(0);
	float PhoHoverE_(0);
	float PhoSigma_(0);
	float PhoChIso_(0);
	float PhoNeuIso_(0);
	float PhoPhoIso_(0);
	bool  PhoPassID_;
	float mcElePt_(0);
	float mcEleEta_(0);
	float mcElePhi_(0);
	float recoEleEt_(0);
	float recoEleEta_(0);
	float recoElePhi_(0);
	float EleR9_(0);
	float EleHoverE_(0);
	float EleSigma_(0);
	float EleChIso_(0);
	float EleNeuIso_(0);
	float ElePhoIso_(0);
	float EleMiniIso_(0);
	float EledEtaIn_(0);
	float EledPhiIn_(0);
	float EleD0_(0);
	float EleDz_(0);
	float EleooEmooP_(0);
	bool  ElePassID_;
	float dRPhoEle_(0);
	float Invmass_(0);
	int   nJet_(0);
	float HT_(0);
	float sigMET_(0);
	double MT_(0), ThreeBodyMass_(0);
	int   nVertex_(0);
	
	tree->Branch("Mass1",          &Mass1_);
	tree->Branch("Mass2",          &Mass2_);
	tree->Branch("branch",         &decaybranch);
	tree->Branch("nVertex",        &nVertex_);
	tree->Branch("MT",&MT_);
	tree->Branch("ThreeBodyMass",&ThreeBodyMass_);
	tree->Branch("mcPhotonEt"     ,&mcPhotonEt_); 
	tree->Branch("mcPhotonEta"    ,&mcPhotonEta_);
	tree->Branch("mcPhotonPhi"    ,&mcPhotonPhi_);
	tree->Branch("recoPhotonEt",   &recoPhotonEt_);
	tree->Branch("recoPhotonEta",  &recoPhotonEta_);
	tree->Branch("recoPhotonPhi",  &recoPhotonPhi_);
	tree->Branch("PhoR9"          ,&PhoR9_);
	tree->Branch("PhoHoverE"      ,&PhoHoverE_);
	tree->Branch("PhoSigma"       ,&PhoSigma_);
	tree->Branch("PhoChIso"       ,&PhoChIso_);
	tree->Branch("PhoNeuIso"      ,&PhoNeuIso_);
	tree->Branch("PhoPhoIso"      ,&PhoPhoIso_);
	tree->Branch("PhoPassID"      ,&PhoPassID_);
	tree->Branch("mcElePt"        ,&mcElePt_);
	tree->Branch("mcEleEta"       ,&mcEleEta_);
	tree->Branch("mcElePhi"       ,&mcElePhi_);
	tree->Branch("recoEleEt",     &recoEleEt_);
	tree->Branch("recoEleEta",    &recoEleEta_);
	tree->Branch("recoElePhi",    &recoElePhi_);
	tree->Branch("EleR9"          ,&EleR9_);
	tree->Branch("EleHoverE"      ,&EleHoverE_);
	tree->Branch("EleSigma"       ,&EleSigma_);
	tree->Branch("EleChIso"       ,&EleChIso_);
	tree->Branch("EleNeuIso"      ,&EleNeuIso_);
	tree->Branch("ElePhoIso"      ,&ElePhoIso_);
	tree->Branch("EleMiniIso"     ,&EleMiniIso_);
	tree->Branch("EledEtaIn"      ,&EledEtaIn_);
	tree->Branch("EledPhiIn"      ,&EledPhiIn_);
	tree->Branch("EleD0"          ,&EleD0_);
	tree->Branch("EleDz"          ,&EleDz_);
	tree->Branch("EleooEmooP"     ,&EleooEmooP_);
	tree->Branch("ElePassID"      ,&ElePassID_);
	tree->Branch("dRPhoEle"       ,&dRPhoEle_);
	tree->Branch("Invmass"        ,&Invmass_);    
	tree->Branch("nJet"           ,&nJet_);
	tree->Branch("HT"             ,&HT_);
	tree->Branch("sigMET"         ,&sigMET_);
	
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
	
	egtree->Branch("Mass1",     &Mass1_);
	egtree->Branch("Mass2",     &Mass2_);
	egtree->Branch("branch",    &decaybranch);
	egtree->Branch("phoEt",     &eg_phoEt);
	egtree->Branch("phoEta",    &eg_phoEta);
	egtree->Branch("phoPhi",    &eg_phoPhi);
	egtree->Branch("lepPt",     &eg_lepPt);
	egtree->Branch("lepEta",    &eg_lepEta);
	egtree->Branch("lepPhi",    &eg_lepPhi);
	egtree->Branch("sigMT",     &eg_sigMT);
	egtree->Branch("sigMET",    &eg_sigMET);
	egtree->Branch("sigMETPhi", &eg_sigMETPhi);
	egtree->Branch("dPhiLepMET",&eg_dPhiLepMET);
	egtree->Branch("nVertex",   &eg_nVertex);
	egtree->Branch("dRPhoLep",  &eg_dRPhoLep);
	egtree->Branch("HT",        &eg_HT);
	egtree->Branch("nJet",      &eg_nJet);
	egtree->Branch("invmass",   &eg_invmass);
	egtree->Branch("sigMETJESup",     &eg_sigMETJESup);
	egtree->Branch("sigMETJESdo",     &eg_sigMETJESdo);
	egtree->Branch("sigMETJERup",     &eg_sigMETJERup);
	egtree->Branch("sigMETJERdo",     &eg_sigMETJERdo);
	egtree->Branch("sigMTJESup",      &eg_sigMTJESup);
	egtree->Branch("sigMTJESdo",      &eg_sigMTJESdo);
	egtree->Branch("sigMTJERup",      &eg_sigMTJERup);
	egtree->Branch("sigMTJERdo",      &eg_sigMTJERdo);
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
	
	mgtree->Branch("Mass1",     &Mass1_);
	mgtree->Branch("Mass2",     &Mass2_);
	mgtree->Branch("branch",    &decaybranch);
	mgtree->Branch("phoEt",     &mg_phoEt);
	mgtree->Branch("phoEta",    &mg_phoEta);
	mgtree->Branch("phoPhi",    &mg_phoPhi);
	mgtree->Branch("lepPt",     &mg_lepPt);
	mgtree->Branch("lepEta",    &mg_lepEta);
	mgtree->Branch("lepPhi",    &mg_lepPhi);
	mgtree->Branch("sigMT",     &mg_sigMT);
	mgtree->Branch("sigMET",    &mg_sigMET);
	mgtree->Branch("sigMETPhi", &mg_sigMETPhi);
	mgtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
	mgtree->Branch("threeMass", &mg_threeMass);
	mgtree->Branch("nVertex",   &mg_nVertex);
	mgtree->Branch("dRPhoLep",  &mg_dRPhoLep);
	mgtree->Branch("HT",        &mg_HT);
	mgtree->Branch("nJet",      &mg_nJet);
	mgtree->Branch("sigMETJESup",     &mg_sigMETJESup);
	mgtree->Branch("sigMETJESdo",     &mg_sigMETJESdo);
	mgtree->Branch("sigMETJERup",     &mg_sigMETJERup);
	mgtree->Branch("sigMETJERdo",     &mg_sigMETJERdo);
	mgtree->Branch("sigMTJESup",      &mg_sigMTJESup);
	mgtree->Branch("sigMTJESdo",      &mg_sigMTJESdo);
	mgtree->Branch("sigMTJERup",      &mg_sigMTJERup);
	mgtree->Branch("sigMTJERdo",      &mg_sigMTJERdo);
	mgtree->Branch("dPhiLepMETJESup", &mg_dPhiLepMETJESup);
	mgtree->Branch("dPhiLepMETJESdo", &mg_dPhiLepMETJESdo);
	mgtree->Branch("dPhiLepMETJERup", &mg_dPhiLepMETJERup);
	mgtree->Branch("dPhiLepMETJERdo", &mg_dPhiLepMETJERdo);
	mgtree->Branch("HTJESup",     &mg_HTJESup);
	mgtree->Branch("HTJESdo",     &mg_HTJESdo);
	
	rawData raw(es, datatype);
	std::vector<mcData>  MCData;
	std::vector<recoPhoton> Photon;
	std::vector<recoEle> Ele;
	std::vector<recoMuon>   Muon;
	std::vector<recoJet>   JetCollection;
	float met(0);
	float metPhi(0);
	float met_T1JERUp(0);
	float met_T1JERDo(0);
	float met_T1JESUp(0);
	float met_T1JESDo(0);	
	float	metPhi_T1JESUp(0);
	float	metPhi_T1JESDo(0);
	float	metPhi_T1UESUp(0);
	float	metPhi_T1UESDo(0);
	int nVtx(0);
	int jetNumber(0);

	const unsigned nEvts = es->GetEntries()/1000;
	//const unsigned nEvts = es->GetEntries();
	std::cout << "total event : " << nEvts << std::endl;

	for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;
		
		Mass1_ = -1;
		Mass2_ = -1;
		mcPhotonEt_ = 0;
		mcPhotonEta_ = 0;
		mcPhotonPhi_ = 0;
		recoPhotonEt_ = 0;
		recoPhotonEta_ = 0;
		recoPhotonPhi_ = 0;
		PhoR9_ = 0;
		PhoHoverE_ = 0;
		PhoSigma_ = 0;
		PhoChIso_ = 0;
		PhoNeuIso_ = 0;
		PhoPhoIso_ = 0;
		PhoPassID_ = false;
		mcElePt_ = 0;
		mcEleEta_ = 0;
		mcElePhi_ = 0;
		recoEleEt_ = 0;
		recoEleEta_ = 0;
		recoElePhi_ = 0;
		EleR9_ = 0;
		EleHoverE_ = 0;
		EleSigma_ = 0;
		EleChIso_ = 0;
		EleNeuIso_ = 0;
		ElePhoIso_ = 0;
		EleMiniIso_ = 0;
		EledEtaIn_ = 0;
		EledPhiIn_ = 0;
		EleD0_ = 0;
		EleDz_ = 0;
		EleooEmooP_ = 0;
		ElePassID_ = 0;
		dRPhoEle_ = 0;
		Invmass_ = 0;
		MT_ = 0;
		ThreeBodyMass_ = 0;
		nJet_ = 0;
		HT_ = 0;
		
		raw.GetData(es, ievt);
		MCData.clear();
		Photon.clear();
		Muon.clear();
		Ele.clear();
		JetCollection.clear();
		if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
		for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
		for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
		for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
		for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
		met = raw.pfMET;
		metPhi = raw.pfMETPhi;
		met_T1JERUp = raw.pfMET_T1JERUp;
		met_T1JERDo = raw.pfMET_T1JERDo;
		met_T1JESUp = raw.pfMET_T1JESUp;
		met_T1JESDo = raw.pfMET_T1JESDo;
		metPhi_T1JESUp = raw.pfMETPhi_T1JESUp;
		metPhi_T1JESDo = raw.pfMETPhi_T1JESDo;
		metPhi_T1UESUp = raw.pfMETPhi_T1UESUp;
		metPhi_T1UESDo = raw.pfMETPhi_T1UESDo;
		nVtx = raw.nVtx;
		// GGM
		if(SUSYtype == 2){
			int scan1 = raw.EventTag->Index("M1");
			int scan2 = raw.EventTag->Index("M2");
			if(scan2 < 1)scan2 = raw.EventTag->Index("M3");
			int scan3 = raw.EventTag->Length();
			TString m1 = (*raw.EventTag)(scan1+2, scan2-(scan1+2)-1);
			TString m2 = (*raw.EventTag)(scan2+2, scan3-(scan2+2));
			Mass1_ = m1.Atoi();
			Mass2_ = m2.Atoi();
		}
		// T5Wg, T6Wg
		else if(SUSYtype == 5 || SUSYtype == 6){
			int scan1 = raw.EventTag->First('_');
			int scan2 = raw.EventTag->Last('_');
			int scan3 = raw.EventTag->Length();
			TString m1 = (*raw.EventTag)(scan1+1, scan2-scan1-1);
			TString m2 = (*raw.EventTag)(scan2+1, scan3-scan2);
			Mass1_ = m1.Atoi();
			Mass2_ = m2.Atoi();
		}
		// TChiWg
		else{
			int scan1 = raw.EventTag->First('_');
			int scan3 = raw.EventTag->Length();
			TString m1 = (*raw.EventTag)(scan1+1, scan3-scan1);
			Mass1_ = m1.Atoi();
			Mass2_ = Mass1_;
		}
	
		nVertex_ = nVtx;
		sigMET_  = met;
		std::vector<mcData>::iterator genPho = MCData.end();
		std::vector<mcData>::iterator genEle = MCData.end();
		std::vector<mcData>::iterator genNeu = MCData.end();
		std::vector<recoPhoton>::iterator recopho = Photon.end();
		std::vector<recoEle>::iterator recoele = Ele.end();

		std::vector< pair<int, int> > finalState;
		finalState.clear();
		std::vector< std::vector<mcData>::iterator > itMCList;
		itMCList.clear();	
	//	if(debugMCLevel > 0)std::cout << std::endl;
		bool hasLepFromPho(false);	
		for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){ 

			if(debugMCLevel >= 1 && itMC->getPt() > 0.1){
				bool inList(false);
				for(unsigned imc(0); imc < itMCList.size(); imc++){
					if( itMC->getPID() == itMCList[imc]->getPID() && itMC->getMomPID() == itMCList[imc]->getMomPID() && DeltaR(itMC->getEta(), itMC->getPhi(), itMCList[imc]->getEta(), itMCList[imc]->getPhi()) < 0.3)inList=true;
				}
				if(!inList)itMCList.push_back(itMC);
			}
			// 21-25 are gluon, gamma, Z, W Higgs
			// 1000022 1000023 neutralino
			// 1000021 gluino
			if(fabs(itMC->getPID()) > 20 && fabs(itMC->getPID()) < 26 && fabs(itMC->getMomPID()) > 1000000){
				finalState.push_back( make_pair(fabs(itMC->getPID()), fabs(itMC->getMomPID())) );
			}
			// photon radiate muon
			if(fabs(itMC->getPID()) == 13 && fabs(itMC->getMomPID()) == 22)hasLepFromPho=true;
			//Look for neutralino-decayed photon
			if(itMC->isPhoton() && itMC->decayFromNeu()){
				genPho = itMC;
				// Match reco photon to signal photon
				float mindR(PI);
				for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
					float tmpdE = fabs(itpho->getEt() - genPho->getEt())/genPho->getEt();
					float tmpdR = DeltaR(itpho->getEta(), itpho->getPhi(), genPho->getEta(), genPho->getPhi());
					if(tmpdE < 0.1 && tmpdR < mindR){
						mindR = tmpdR;
						if(tmpdR < MCMatchCone)recopho = itpho;
					}
				}
			}//endif: found signal photon 
		     
			if(abs(itMC->getPID()) == NvEPID && abs(itMC->getGMomPID())== CharginoPID){
   			genNeu = itMC;
			}

			//Look for chargino-decayed electron
			if(itMC->isElectron() && abs(itMC->getMomPID())== WPID){
			
				genEle = itMC;
				// Match reco electron to signal electron
				float mindR(PI);
				for(std::vector<recoEle>::iterator itele = Ele.begin() ; itele != Ele.end(); ++itele){
					float tmpdE = fabs(itele->getPt() - genEle->getPt())/genEle->getPt();
					float tmpdR = DeltaR(itele->getEta(), itele->getPhi(), genEle->getEta(), genEle->getPhi());
					if(tmpdE < 0.1 && tmpdR < mindR){
						mindR = tmpdR;
						if(tmpdR < MCMatchCone)recoele = itele;
					}
				}
			}//endif:: found signal ELectron
		
		}//loop on MC particles
     
		bool hasBranchPho(false);
		bool hasBranchZ(false);
		bool hasBranchH(false);
		bool hasBranchW(false);
		for(int i(0); i < finalState.size(); i++){
			if(finalState[i].first == 22)hasBranchPho = true;
			if(finalState[i].first == 23)hasBranchZ = true;
			if(finalState[i].first == 24)hasBranchW = true;
			if(finalState[i].first == 25)hasBranchH = true;
		}
		if(hasBranchPho && hasBranchZ)decaybranch = branchType::gZ;
		else if(hasBranchPho && hasBranchH)decaybranch = branchType::gH;
		else if(hasBranchPho && hasBranchW)decaybranch = branchType::gW;
		else if(hasBranchPho && !hasBranchZ && !hasBranchH && !hasBranchW)decaybranch = branchType::gg;
		else decaybranch = branchType::noType;
		// Calculate HT 
		for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
			if(!itJet->passSignalSelection())continue;
			// well separated jets
			if(genPho != MCData.end())
				if(DeltaR(itJet->getEta(), itJet->getPhi(), genPho->getEta(),genPho->getPhi()) <= AK4Cone)continue;	
			if(genEle != MCData.end())
				if(DeltaR(itJet->getEta(), itJet->getPhi(), genEle->getEta(),genEle->getPhi()) <= AK4Cone)continue;
			nJet_ += 1;
			HT_ += itJet->getPt();
		}
		if(genPho != MCData.end()){
			mcPhotonEt_ = genPho->getEt();
			mcPhotonEta_= genPho->getEta();
			mcPhotonPhi_= genPho->getPhi();
		}
		if(recopho != Photon.end()){ 
			recoPhotonEt_ =recopho->getEt();
			recoPhotonEta_=recopho->getEta();
			recoPhotonPhi_=recopho->getPhi();
			PhoR9_    =recopho->getR9();
			PhoHoverE_=recopho->getHoverE();
			PhoSigma_ =recopho->getSigma();
			PhoChIso_ =recopho->getChIso();
			PhoNeuIso_=recopho->getNeuIso();
			PhoPhoIso_=recopho->getPhoIso();
			PhoPassID_=recopho->isLoose();
		}
		if(genEle != MCData.end()){
			mcElePt_ =genEle->getPt();
			mcEleEta_=genEle->getEta();
			mcElePhi_=genEle->getPhi();
		}
		if(recoele != Ele.end()){ 
			recoEleEt_ =recoele->getPt();
			recoEleEta_=recoele->getEta();
			recoElePhi_=recoele->getPhi();
			EleR9_     =recoele->getR9();
			EleHoverE_ =recoele->getHoverE();
			EleSigma_  =recoele->getSigma();
			EleChIso_  =recoele->getChIso();
			EleNeuIso_ =recoele->getNeuIso();
			ElePhoIso_ =recoele->getPhoIso();
			EleMiniIso_=recoele->getMiniIso();
			EledEtaIn_ =recoele->getdEtaIn();
			EledPhiIn_ =recoele->getdPhiIn();
			EleD0_     =recoele->getD0();
			EleDz_     =recoele->getDz();
			EleooEmooP_=recoele->getEoverPInv();
			ElePassID_ =recoele->isMiniMedium(); 
		}
		if(genPho != MCData.end() && genEle != MCData.end()){
			dRPhoEle_ = DeltaR(genPho->getEta(), genPho->getPhi(), genEle->getEta(), genEle->getPhi());
			Invmass_  = (genPho->getP4() + genEle->getP4()).M();
			MT_       = sqrt(2*met*genEle->getPt()*(1-std::cos(DeltaR(0, genEle->getPhi(), 0, metPhi))));
  		ThreeBodyMass_ = sqrt(2*met*(genPho->getP4()+genEle->getP4()).Pt()*(1-std::cos(DeltaR(0, (genPho->getP4()+genEle->getP4()).Phi(), 0, metPhi))));
		}
		// Fill the MC tree
		tree->Fill();




		/***************   Start processing the event selection ****************************************/

		// Veto di-photon events. Disabled //
		//int Nmedpho(0);
		//for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
		//	if(!itpho->isMedium())continue;
		//	if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
		//	if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
		//	if(itpho->getSigma() < 0.005)continue;
		//	if(itpho->PixelSeed() != 0)continue;
		//
		//	Nmedpho+=1;
		//}
		////if(Nmedpho >= 2)continue;

		std::vector<recoPhoton>::iterator signalPho = Photon.end(); // signal photon for e+g channel
		for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
			if(itpho->getR9() < R9EBCut)continue;
			if(!itpho->passSignalSelection())continue;
			bool PixelVeto = itpho->PixelSeed()==0? true: false;
			bool GSFveto(true);
			bool FSRVeto(true);
			for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < ElectronVetoCone)GSFveto = false; // electron veto
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < GeneralCone)FSRVeto=false;
			}
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < GeneralCone && im->getEt()>2.0)FSRVeto=false;

			if(GSFveto && PixelVeto && FSRVeto){
				if(signalPho == Photon.end())
					signalPho = itpho; // leading photon
			}
		}

		// select leading electron
		std::vector<recoEle>::iterator egsignalEle = Ele.end();
		for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
			if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;
			if(itEle->passSignalSelection()){
				if(egsignalEle == Ele.end())egsignalEle = itEle;
				else if(itEle->getPt() > egsignalEle->getPt())egsignalEle = itEle;
			}
		}

		std::vector<recoMuon>::iterator signalMu = Muon.end();
		for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
			if(itMu->passSignalSelection()){
				if(signalMu == Muon.end())signalMu = itMu;
				else if(itMu->getPt() > signalMu->getPt())signalMu = itMu;
			}
		}

		// If isolated photon and electron exist, process the e+g tree
		if(signalPho != Photon.end() && egsignalEle != Ele.end()){ 
			double dReg = DeltaR(signalPho->getEta(), signalPho->getPhi(), egsignalEle->getEta(), egsignalEle->getPhi()); 
			if(dReg > PhoLepdR){ // 0.8 dR
				if(((signalPho->getP4()+egsignalEle->getP4()).M() - ZMass) > ZMassWindow){ // 10 GeV Zmass cut

					float deltaPhi = DeltaPhi(egsignalEle->getPhi(), metPhi);
					float mt = sqrt(2*met*egsignalEle->getPt()*(1-std::cos(deltaPhi)));
					eg_phoEt = signalPho->getCalibEt();
					eg_phoEta= signalPho->getEta();
					eg_phoPhi= signalPho->getPhi();
					eg_lepPt = egsignalEle->getCalibPt();
					eg_lepEta= egsignalEle->getEta();
					eg_lepPhi= egsignalEle->getPhi();
					eg_sigMT = mt;
					eg_sigMET= met;
					eg_sigMETPhi = metPhi;
					eg_dPhiLepMET = deltaPhi; 
					eg_nVertex = nVtx; 
					eg_dRPhoLep= dReg;
					eg_invmass = (signalPho->getP4()+egsignalEle->getP4()).M();
					eg_sigMETJESup = met_T1JESUp;
					eg_sigMETJESdo = met_T1JESDo;
					eg_sigMETJERup = met_T1JERUp;
					eg_sigMETJERdo = met_T1JERDo;
					eg_dPhiLepMETJESup = DeltaPhi(egsignalEle->getPhi(), metPhi_T1JESUp);
					eg_dPhiLepMETJESdo = DeltaPhi(egsignalEle->getPhi(), metPhi_T1JESDo);
					eg_dPhiLepMETJERup = deltaPhi;
					eg_dPhiLepMETJERdo = deltaPhi;
					eg_sigMTJESup = sqrt(2*met_T1JESUp*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESup)));
					eg_sigMTJESdo = sqrt(2*met_T1JESDo*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESdo)));
					eg_sigMTJERup = sqrt(2*met_T1JERUp*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERup)));
					eg_sigMTJERdo = sqrt(2*met_T1JERDo*egsignalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERdo)));
	
					eg_nJet = 0;
					eg_HT = 0;
					eg_HTJESup = 0;
					eg_HTJESdo = 0;
					for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
						if(!itJet->passSignalSelection())continue;
						if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= AK4Cone)continue;	
						if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalEle->getEta(),egsignalEle->getPhi()) <= AK4Cone)continue;
						eg_nJet += 1;
						eg_HT += itJet->getPt();
						eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
						eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
					}	
					// Fill the e+g tree
					egtree->Fill();
				}// Z mass Filter
			}//dR filter
		}// ele + pho candidate
		 
		if( signalPho != Photon.end() && signalMu != Muon.end() && egsignalEle == Ele.end()){// mu+e+g events are counted in e+g channel
			double dRmg = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
			if(dRmg>PhoLepdR){
				float deltaPhi = DeltaPhi(signalMu->getPhi(), metPhi);
				float mt = sqrt(2*met*signalMu->getPt()*(1-std::cos(deltaPhi)));
				float threeBodyMass = sqrt(2*met*(signalPho->getP4()+ signalMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (signalPho->getP4()+signalMu->getP4()).Phi(), 0, metPhi))));
		
				mg_phoEt = signalPho->getCalibEt();
				mg_phoEta= signalPho->getEta();
				mg_phoPhi= signalPho->getPhi();
				mg_lepPt = signalMu->getPt();
				mg_lepEta= signalMu->getEta();
				mg_lepPhi= signalMu->getPhi();
				mg_sigMT = mt;
				mg_sigMET= met;
				mg_sigMETPhi = metPhi;
				mg_dPhiLepMET = deltaPhi;
				mg_threeMass = threeBodyMass;
				mg_nVertex = nVtx;
				mg_dRPhoLep= dRmg;
				mg_sigMETJESup = met_T1JESUp;
				mg_sigMETJESdo = met_T1JESDo;
				mg_sigMETJERup = met_T1JERUp;
				mg_sigMETJERdo = met_T1JERDo;
				mg_dPhiLepMETJESup = DeltaPhi(signalMu->getPhi(), metPhi_T1JESUp);
				mg_dPhiLepMETJESdo = DeltaPhi(signalMu->getPhi(), metPhi_T1JESDo);
				mg_dPhiLepMETJERup = deltaPhi;
				mg_dPhiLepMETJERdo = deltaPhi;
				mg_sigMTJESup = sqrt(2*met_T1JESUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESup)));
				mg_sigMTJESdo = sqrt(2*met_T1JESDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESdo)));
				mg_sigMTJERup = sqrt(2*met_T1JERUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERup)));
				mg_sigMTJERdo = sqrt(2*met_T1JERDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERdo)));
		
				mg_nJet = 0;
				mg_HT = 0;
				mg_HTJESup = 0;
				mg_HTJESdo = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= AK4Cone)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= AK4Cone)continue;
					mg_nJet += 1;
					mg_HT += itJet->getPt();
					mg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
					mg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
				}
				mgtree->Fill();
			}//dR Filter
		}//Candidate Filter

		std::vector< decayChain > chains;
		chains.clear();
		for(unsigned imc(0); imc<itMCList.size(); imc++){
			decayChain element;
			element.iter = itMCList[imc];
			element.daughter.clear();
			for(unsigned id(0); id < itMCList.size(); id++){
				if(itMCList[imc]->getPID() == itMCList[id]->getMomPID() && DeltaR( itMCList[imc]->getEta(), itMCList[imc]->getPhi(), itMCList[id]->getmomEta(), itMCList[id]->getmomPhi()) < 0.3)element.daughter.push_back( itMCList[id] );
			}
			//if(element.daughter.size() > 0)chains.push_back(element);
			chains.push_back(element);
		}

		std::vector< decayChain > branchOne;
		std::vector< decayChain > branchTwo;
		branchOne.clear();
		branchTwo.clear(); 

		for(unsigned ic(0); ic < chains.size(); ic++){
			if(fabs(chains[ic].iter->getMomPID()) > 100000 ){
				bool findbranch(false);
				for(unsigned ione(0); ione < branchOne.size(); ione++){
					if( fabs(chains[ic].iter->getMomPID()) == fabs(branchOne[ione].iter->getMomPID()) && DeltaR(chains[ic].iter->getmomEta(), chains[ic].iter->getmomPhi(), branchOne[ione].iter->getmomEta(), branchOne[ione].iter->getmomPhi()) < 0.3){branchOne.push_back(chains[ic]); findbranch=true;break;	}
				}
				for(unsigned itwo(0); itwo < branchTwo.size(); itwo++){
					if( fabs(chains[ic].iter->getMomPID()) == fabs(branchTwo[itwo].iter->getMomPID()) && DeltaR(chains[ic].iter->getmomEta(), chains[ic].iter->getmomPhi(), branchTwo[itwo].iter->getmomEta(), branchTwo[itwo].iter->getmomPhi()) < 0.3){branchTwo.push_back(chains[ic]); findbranch=true;break;	}
				}
				if(!findbranch && branchOne.size() == 0)branchOne.push_back(chains[ic]);
				else if(!findbranch && branchTwo.size() == 0)branchTwo.push_back(chains[ic]);
			}
		}
		
			//	if(decaybranch != 1)continue;
			//	if(!acceptMu)continue;
			//	else std::cout << "acc" << std::endl;
			//	for(unsigned imc(0); imc<itMCList.size(); imc++){
			//		printf("%10d %8.2f %8.2f %8.2f   %10d %8.2f %8.2f %8.2f \n", itMCList[imc]->getPID(), itMCList[imc]->getPt(), itMCList[imc]->getEta(), itMCList[imc]->getPhi(), itMCList[imc]->getMomPID(), itMCList[imc]->getmomPt(), itMCList[imc]->getmomEta(), itMCList[imc]->getmomPhi());
			//	}
				
				for(unsigned ione(0); ione < branchOne.size(); ione++){
					if(ione == 0)std::cout << branchOne[ione].iter->getMomPID() << "-> ";
					std::cout <<  branchOne[ione].iter->getPID() << " ";
				}
				std::cout << std::endl;
				for(unsigned ione(0); ione < branchOne.size(); ione++){
					if(branchOne[ione].daughter.size() > 0){
						std::cout << branchOne[ione].iter->getPID() << "-> ";
						for(unsigned id(0); id < branchOne[ione].daughter.size();  id++)std::cout << (branchOne[ione].daughter)[id]->getPID() << " ";
						std::cout << std::endl;
					}
				}
				for(unsigned ione(0); ione < branchTwo.size(); ione++){
					if(ione == 0)std::cout << branchTwo[ione].iter->getMomPID() << "-> ";
					std::cout <<  branchTwo[ione].iter->getPID() << " ";
				}
				std::cout << std::endl;
				for(unsigned ione(0); ione < branchTwo.size(); ione++){
					if(branchTwo[ione].daughter.size() > 0){
						std::cout << "\t" << branchTwo[ione].iter->getPID() << "-> ";
						for(unsigned id(0); id < branchTwo[ione].daughter.size();  id++)std::cout << (branchTwo[ione].daughter)[id]->getPID() << " ";
						std::cout << std::endl;
					}
				}

				bool hasPhotonConv(false);
				for(unsigned ione(0); ione < branchOne.size(); ione++){
					if(branchOne[ione].daughter.size() >= 2){
						if(branchOne[ione].iter->getPID() == 22 && fabs((branchOne[ione].daughter)[0]->getPID()) == 11 && fabs((branchOne[ione].daughter)[1]->getPID()) == 11){hasPhotonConv = true; std::cout << "conv " << branchOne[ione].iter->getPt() << " " << branchOne[ione].iter->getEta() << " " << branchOne[ione].iter->getPhi() << std::endl;}
					}
				}
				for(unsigned ione(0); ione < branchTwo.size(); ione++){
					if(branchTwo[ione].daughter.size() >= 2){
						if(branchTwo[ione].iter->getPID() == 22 && fabs((branchTwo[ione].daughter)[0]->getPID()) == 11 && fabs((branchTwo[ione].daughter)[1]->getPID()) == 11){hasPhotonConv = true; std::cout << "conv " << branchTwo[ione].iter->getPt() << " " << branchTwo[ione].iter->getEta() << " " << branchTwo[ione].iter->getPhi() << std::endl;}
					}
				}

				int NrecoPho(0);
					for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
						if(!itpho->passSignalSelection())continue;
						bool PixelVeto = itpho->PixelSeed()==0? true: false;
						if(PixelVeto ){
							NrecoPho += 1;
						}
					}
				
				
				if(decaybranch== 1){
					if(hasPhotonConv)std::cout << "conve " << NrecoPho << std::endl;
					else std::cout << "prompt " << NrecoPho << std::endl;	
				}

	//	int iNEU(0);	
	//	for(unsigned ic(0); ic < chains.size(); ic++){
	//		if(chains[ic].iter->getPID() == 1000022)iNEU+=1;
	//		if(fabs(chains[ic].iter->getPID()) > 100000 )std::cout << chains[ic].iter->getMomPID() << " "  << chains[ic].iter->getPID() << "->";
	//		else if(fabs(chains[ic].iter->getMomPID()) > 100000 && fabs(chains[ic].iter->getGMomPID()) > 100000 )std::cout << "\t\t" << chains[ic].iter->getMomPID() << " " << chains[ic].iter->getPID() << "->";
	//		else if(fabs(chains[ic].iter->getMomPID()) > 100000)std::cout <<  "\t" << chains[ic].iter->getMomPID() << " "  << chains[ic].iter->getPID() << "->";
	//		else continue;
	//		for(unsigned id(0); id < chains[ic].daughter.size();  id++)std::cout << (chains[ic].daughter)[id]->getPID() << " ";
	//		std::cout << std::endl;
	//	}
	//	if(iNEU != 2)std::cout << "NO" << std::endl;
	//	std::cout << "type: " << decaybranch << std::endl;
	
	}//loop on entries


	output->Write();
	output->Close();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_SUSY_modified_tribeni(atoi(argv[1]),argv[2]);
    return 0;
}
