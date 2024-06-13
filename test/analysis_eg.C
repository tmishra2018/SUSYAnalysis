// g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_eg.C -o analysis_eg.exe `root-config --libs`
#include "../include/analysis_commoncode.h"
#include "../include/analysis_cuts.h"
bool apply_HEMveto=false;

void analysis_eg(int RunYear, const char *Era){//main

  gSystem->Load("/uscms/home/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;

  RunType datatype;
  if(RunYear==2016) datatype = DoubleEG2016;
  if(RunYear==2017) datatype = DoubleEG2017;
  if(RunYear==2018) datatype = DoubleEG2018;
  bool  isMC(false);
  if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add(Form("/eos/uscms/store/group/lpcsusyphotons/Tribeni/DoubleEG/DoubleEG_%d%s.root",RunYear,Era));

  if(RunYear==2018) apply_HEMveto=true; 
  
  const unsigned nEvts = es->GetEntries()/10.;
  //const unsigned nEvts = es->GetEntries();
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_"<<RunYear<<Era<<".root" << std::endl;

  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,Era),"RECREATE");
  outputfile->cd();
  TH1D *p_METFilter = new TH1D("p_METFilter","",12,-2,10);	
  TH1D *p_invmass = new TH1D("p_invmass","",200,0,200);	

  int nBJet(0);
  //************ Signal Tree **********************//
  TTree *sigtree = new TTree("signalTree","signalTree");
  int   run(0);
  Long64_t  event(0);
  int   lumis(0);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
  float trailPt(0);
  float trailEta(0);
  float trailPhi(0);
  std::vector<int>   mcPID;
  std::vector<float> mcEta;
  std::vector<float> mcPhi;
  std::vector<float> mcPt;
  std::vector<int>   mcMomPID;
  std::vector<int>   mcGMomPID;

  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("lumis",     &lumis);
  sigtree->Branch("phoEt",     &phoEt);
  sigtree->Branch("phoEta",    &phoEta);
  sigtree->Branch("phoPhi",    &phoPhi);
  sigtree->Branch("lepPt",     &lepPt);
  sigtree->Branch("lepEta",    &lepEta);
  sigtree->Branch("lepPhi",    &lepPhi);
  sigtree->Branch("sigMT",     &sigMT);
  sigtree->Branch("sigMET",    &sigMET);
  sigtree->Branch("sigMETPhi", &sigMETPhi);
  sigtree->Branch("dPhiLepMET",&dPhiLepMET);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
  sigtree->Branch("nBJet",     &nBJet);
  sigtree->Branch("trailPt",   &trailPt);
  sigtree->Branch("trailEta",  &trailEta);
  sigtree->Branch("trailPhi",  &trailPhi);
  if(isMC){
  	sigtree->Branch("mcPID",     &mcPID);
  	sigtree->Branch("mcEta",     &mcEta);
  	sigtree->Branch("mcPhi",     &mcPhi);
  	sigtree->Branch("mcPt",      &mcPt);
  	sigtree->Branch("mcMomPID",  &mcMomPID);
  	sigtree->Branch("mcGMomPID", &mcGMomPID);
  }

//************ Signal Tree **********************//
  TTree *proxytree = new TTree("proxyTree","proxyTree");
  float proxyphoEt(0);
  float proxyphoEta(0);
  float proxyphoPhi(0);
  float proxylepPt(0);
  float proxylepEta(0);
  float proxylepPhi(0);
  float proxysigMT(0);
  float proxysigMET(0);
  float proxysigMETPhi(0);
  float proxydPhiLepMET(0);
  float proxythreeMass(0);
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  float proxynJet(0);
  int   proxyFSRVeto(0);
  float proxytrailPt(0);
  float proxytrailEta(0);
  float proxytrailPhi(0);
  
  proxytree->Branch("run",       &run);
  proxytree->Branch("event",     &event);
  proxytree->Branch("lumis",     &lumis);
  proxytree->Branch("phoEt",     &proxyphoEt);
  proxytree->Branch("phoEta",    &proxyphoEta);
  proxytree->Branch("phoPhi",    &proxyphoPhi);
  proxytree->Branch("lepPt",     &proxylepPt);
  proxytree->Branch("lepEta",    &proxylepEta);
  proxytree->Branch("lepPhi",    &proxylepPhi);
  proxytree->Branch("sigMT",     &proxysigMT);
  proxytree->Branch("sigMET",    &proxysigMET);
  proxytree->Branch("sigMETPhi", &proxysigMETPhi);
  proxytree->Branch("dPhiLepMET",&proxydPhiLepMET);
  proxytree->Branch("threeMass", &proxythreeMass);
  proxytree->Branch("nVertex",   &proxynVertex);
  proxytree->Branch("dRPhoLep",  &proxydRPhoLep);
  proxytree->Branch("HT",        &proxyHT);
  proxytree->Branch("nJet",      &proxynJet);
  proxytree->Branch("nBJet",     &nBJet);
  proxytree->Branch("FSRVeto",   &proxyFSRVeto);
  proxytree->Branch("trailPt",   &proxytrailPt);
  proxytree->Branch("trailEta",  &proxytrailEta);
  proxytree->Branch("trailPhi",  &proxytrailPhi);

//************ Signal Tree **********************//
  TTree *jettree = new TTree("jetTree","jetTree");
  float jetphoEt(0);
  float jetphoEta(0);
  float jetphoPhi(0);	
  float jetphoChIso(0);
  float jetphoSigma(0);
  float jetlepPt(0);
  float jetlepEta(0);
  float jetlepPhi(0);
  float jetsigMT(0);
  float jetsigMET(0);
  float jetsigMETPhi(0);
  float jetdPhiLepMET(0);
  int   jetnVertex(0);
  float jetdRPhoLep(0);
  float jetHT(0);
  float jetnJet(0);
  float jettrailPt(0);
  float jettrailEta(0);
  float jettrailPhi(0);
  
  jettree->Branch("run",       &run);
  jettree->Branch("event",     &event);
  jettree->Branch("lumis",     &lumis);
  jettree->Branch("phoEt",     &jetphoEt);
  jettree->Branch("phoEta",    &jetphoEta);
  jettree->Branch("phoPhi",    &jetphoPhi);
  jettree->Branch("phoChIso",  &jetphoChIso);
  jettree->Branch("phoSigma",  &jetphoSigma);
  jettree->Branch("lepPt",     &jetlepPt);
  jettree->Branch("lepEta",    &jetlepEta);
  jettree->Branch("lepPhi",    &jetlepPhi);
  jettree->Branch("sigMT",     &jetsigMT);
  jettree->Branch("sigMET",    &jetsigMET);
  jettree->Branch("sigMETPhi", &jetsigMETPhi);
  jettree->Branch("dPhiLepMET",&jetdPhiLepMET);
  jettree->Branch("nVertex",   &jetnVertex);
  jettree->Branch("dRPhoLep",  &jetdRPhoLep);
  jettree->Branch("HT",        &jetHT);
  jettree->Branch("nJet",      &jetnJet);
  jettree->Branch("nBJet",     &nBJet);
  jettree->Branch("trailPt",   &jettrailPt);
  jettree->Branch("trailEta",  &jettrailEta);
  jettree->Branch("trailPhi",  &jettrailPhi);
  
//*********** fake lepton *********************//
  TTree *fakeLeptree = new TTree("fakeLepTree","fakeLepTree");
  float fakeLepphoEt(0);
  float fakeLepphoEta(0);
  float fakeLepphoPhi(0);
  float fakeLepPt(0);
  float fakeLepEta(0);
  float fakeLepPhi(0);
  float fakeLepMiniIso(0);
  int   fakeLepIsStandardProxy(0);
  float fakeLepSigma(0);
  float fakeLepdEta(0);
  float fakeLepdPhi(0);
  unsigned fakeLepIndex(0);
  float fakeLepsigMT(0);
  float fakeLepsigMET(0);
  float fakeLepsigMETPhi(0);
  float fakeLepdPhiLepMET(0);
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  float fakeLepnJet(0);
  float fakeLeptrailPt(0);
  float fakeLeptrailEta(0);
  float fakeLeptrailPhi(0);
  
  
  fakeLeptree->Branch("run",       &run);
  fakeLeptree->Branch("event",     &event);
  fakeLeptree->Branch("lumis",     &lumis);
  fakeLeptree->Branch("phoEt",     &fakeLepphoEt);
  fakeLeptree->Branch("phoEta",    &fakeLepphoEta);
  fakeLeptree->Branch("phoPhi",    &fakeLepphoPhi);
  fakeLeptree->Branch("lepPt",     &fakeLepPt);
  fakeLeptree->Branch("lepEta",    &fakeLepEta);
  fakeLeptree->Branch("lepPhi",    &fakeLepPhi);
  fakeLeptree->Branch("fakeLepMiniIso",&fakeLepMiniIso);
  fakeLeptree->Branch("fakeLepIsStandardProxy", &fakeLepIsStandardProxy);
  fakeLeptree->Branch("fakeLepSigma",&fakeLepSigma);
  fakeLeptree->Branch("fakeLepdEta", &fakeLepdEta);
  fakeLeptree->Branch("fakeLepdPhi", &fakeLepdPhi);
  fakeLeptree->Branch("fakeLepIndex",&fakeLepIndex);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
  fakeLeptree->Branch("nBJet",     &nBJet);
  fakeLeptree->Branch("trailPt",   &fakeLeptrailPt);
  fakeLeptree->Branch("trailEta",  &fakeLeptrailEta);
  fakeLeptree->Branch("trailPhi",  &fakeLeptrailPhi);

//*************** for jet-photon fake rate ***********************//
	TTree *hadrontree = new TTree("hadronTree","hadronTree");
	float hadron_phoEt(0);
	float hadron_phoEta(0);
	float hadron_phoPhi(0);
	float hadron_phoSigma(0);
	float hadron_phoChIso(0);
  	std::vector<float> hadron_eleproxyEt;
	std::vector<float> hadron_eleproxyEta;
	std::vector<float> hadron_eleproxyPhi;
	std::vector<float> hadron_eleproxySigma;
	std::vector<float> hadron_eleproxyChIso;
	std::vector<int> hadron_eleproxynVertex;
	float hadron_sigMT(0);
	float hadron_sigMET(0);
	float hadron_sigMETPhi(0);
	float hadron_dPhiLepMET(0);
	int   hadron_nVertex(0);
	float hadron_HT(0);
	float hadron_nJet(0);
	std::vector<int> hadron_mcPID;
	std::vector<float> hadron_mcEta;
	std::vector<float> hadron_mcPhi;
	std::vector<float> hadron_mcPt;
	std::vector<int> hadron_mcMomPID;

	hadrontree->Branch("phoEt",     &hadron_phoEt);
	hadrontree->Branch("phoEta",    &hadron_phoEta);
	hadrontree->Branch("phoPhi",    &hadron_phoPhi);
	hadrontree->Branch("phoSigma",  &hadron_phoSigma);
	hadrontree->Branch("phoChIso",  &hadron_phoChIso);
	hadrontree->Branch("eleproxyEt",&hadron_eleproxyEt);
	hadrontree->Branch("eleproxyEta",&hadron_eleproxyEta);
	hadrontree->Branch("eleproxyPhi",&hadron_eleproxyPhi);
	hadrontree->Branch("eleproxySigma", &hadron_eleproxySigma);
	hadrontree->Branch("eleproxyChIso", &hadron_eleproxyChIso);
	hadrontree->Branch("eleproxynVertex",&hadron_eleproxynVertex);
	hadrontree->Branch("sigMT",     &hadron_sigMT);
	hadrontree->Branch("sigMET",    &hadron_sigMET);
	hadrontree->Branch("sigMETPhi", &hadron_sigMETPhi);
	hadrontree->Branch("dPhiLepMET",&hadron_dPhiLepMET);
	hadrontree->Branch("nVertex",   &hadron_nVertex);
	hadrontree->Branch("HT",        &hadron_HT);
	hadrontree->Branch("nJet",      &hadron_nJet);
	if(isMC){
		hadrontree->Branch("mcPID",     &hadron_mcPID);
		hadrontree->Branch("mcEta",     &hadron_mcEta);
		hadrontree->Branch("mcPhi",     &hadron_mcPhi);
		hadrontree->Branch("mcPt",      &hadron_mcPt);
		hadrontree->Branch("mcMomPID",  &hadron_mcMomPID);
	}
//*********** histo list **********************//
  TH1F *p_eventcount = new TH1F("p_eventcount","p_eventcount",7,0,7);

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  float MET(0);
  float METPhi(0);
  int nVtx(0);
  int jetNumber(0);
  int METFilter(0);
  logfile << "RunType: " << datatype << std::endl;

  int passHEM(0);

  std::cout << "Total evetns : " << nEvts << std::endl;
  logfile << "Total evetns : " << nEvts << std::endl;
	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

		if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if (ievt%1000000==0) logfile  << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			MCData.clear();
			Photon.clear();
			Muon.clear();
			Ele.clear();
			JetCollection.clear();
			if(isMC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			nVtx = raw.nVtx;
			run=raw.run;
			event=raw.event;
			lumis=raw.lumis;

			if(RunYear==2018 && !passHEMVeto(0,raw)) continue;
        		passHEM++;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;
			// eg combination
			if(raw.nEle < 1 || raw.nPho <1)continue;

  //    int Nmedpho(0);
  //    for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
  //      if(!itpho->isMedium())continue;
  //      if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
  //      if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
  //      if(itpho->getSigma() < 0.005)continue;
  //      if(itpho->PixelSeed() != 0)continue;

  //      Nmedpho+=1;
  //    }
  //    if(Nmedpho >= 2)continue;

			nBJet = 0;
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(itJet->getPt() < 20)continue;
				if(itJet->isBJet())nBJet+=1;
			}

			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Photon                                              **********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			std::vector< int > proxyPhoFSRVeto;
			proxyPhoCollection.clear();
			proxyPhoFSRVeto.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			// ******* for jet-photon fake rate***************//
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  hadeleproxyPhoCollection;
			hadeleproxyPhoCollection.clear();

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < R9EBCut)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				//bool passSigma = itpho->passSigmaOLD(1);
				//bool passChIso = itpho->passChIsoOLD(1);
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool photonFSRVeto(true);
				bool eleFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= ElectronVetoCone/* 0.02 */)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < GeneralCone /* 0.3 */ && ie->getEt()>2.0)photonFSRVeto=false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < GeneralCone && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > ElectronVetoCone)eleFSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < GeneralCone && im->getEt()>2.0){
						photonFSRVeto=false;
						eleFSRVeto=false;
					}
				}
				// ******** very loose, before sigma and isolation cut.  For jet-photon fake rate, and hadron proxy ************//	
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasHadronPho){
						hasHadronPho = true;
						hadronPho = itpho;
					}

					if(!passSigma || !passChIso){
						// collection for jet->photon fake
						if( (itpho->getSigma()< 0.02 && itpho->isEB()) || (itpho->getSigma()< 0.04 && itpho->isEE()) )jetPhoCollection.push_back(itpho);
					}
				}
				// collection for jets->lepton fake
				else if( (!GSFveto || !PixelVeto) )hadeleproxyPhoCollection.push_back(itpho);

        // ****************  standard ID ************************************//
				if(!itpho->passSignalSelection())continue;
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

				if((!PixelVeto || !GSFveto)){
						// collection for ele->photon fake
						proxyPhoCollection.push_back(itpho);
						proxyPhoFSRVeto.push_back( eleFSRVeto? 1:0);
				}
			}

			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Lepton                                              **********************************************/
			bool hasLep(false);
			std::vector<recoEle>::iterator signalLep = Ele.begin();
			std::vector< std::vector<recoEle>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoEle>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){

				if(itEle->getCalibPt() < 25)continue;
				//  collection for jets->Lep fake
				//if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
				//if(itEle->isFakeProxyOLD())fakeLepCollection.push_back(itEle);	
				if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;

				if(!itEle->passHLTSelection())continue;
				//if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
				if(itEle->isLooseFakeProxy())fakeLepCollection.push_back(itEle);//Loose the proxy definition	
				if(itEle->passSignalSelection()){
					proxyLepCollection.push_back(itEle);
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itEle;
					}
				}
			}

			bool hasTrail(false);
			std::vector<recoMuon>::iterator trailLep = Muon.begin();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 25)continue;
				if(itMu->passSignalSelection()){
					if(!hasTrail){
						hasTrail = true;
						trailLep = itMu;
					}
				}
			}

			// Filling ele + pho signalTree
			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi()); 
				if(dRlepphoton > 0.8){
					npassdR+=1;

					float tempdeltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float tempMT = sqrt(2*MET*signalLep->getCalibPt()*(1-std::cos(tempdeltaPhi)));
					if(fabs(signalPho->getEta()) < 1.4442 && tempMT > 100 && MET > 100)p_invmass->Fill((signalPho->getCalibP4()+signalLep->getCalibP4()).M());
					if(((signalPho->getCalibP4()+signalLep->getCalibP4()).M() - 91.188) > 10.0){

						npassZ+=1;
						p_METFilter->Fill(-2);
						p_METFilter->Fill(raw.failFilterStep(METFilter));	
						if(raw.passMETFilter(METFilter)){
							npassMETFilter +=1;

							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getCalibPt()*(1-std::cos(deltaPhi)));
							phoEt = signalPho->getCalibEt();
							phoEta= signalPho->getEta();
							phoPhi= signalPho->getPhi();
							lepPt = signalLep->getCalibPt();
							lepEta= signalLep->getEta();
							lepPhi= signalLep->getPhi();
							sigMT = MT;
							sigMET= MET;
							sigMETPhi = METPhi;
							dPhiLepMET = deltaPhi; 
							nVertex = nVtx; 
							dRPhoLep= dRlepphoton;
							if(hasTrail){
								trailPt = trailLep->getPt();
								trailEta = trailLep->getEta();
								trailPhi = trailLep->getPhi();
							}
							else{
								trailPt  = 0; 
								trailEta = 0;
								trailPhi = 0;
							}

							nJet = 0;
							HT = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
								nJet += 1;
								HT += itJet->getPt();
							}	

          	  mcPID.clear();
          	  mcEta.clear();
          	  mcPhi.clear();
          	  mcPt.clear();
          	  mcMomPID.clear();
          	  mcGMomPID.clear();
							if(isMC){
          	 	 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
          	 	   if(itMC->getEt() < 1.0)continue;
          	 	   float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
          	 	   if(mcdR < 0.3){
          	 	     mcPID.push_back(itMC->getPID());
          	 	     mcMomPID.push_back(itMC->getMomPID());
          	 	     mcGMomPID.push_back(itMC->getGMomPID());
          	 	     mcEta.push_back(itMC->getEta());
          	 	     mcPhi.push_back(itMC->getPhi());
          	 	     mcPt.push_back(itMC->getEt());
          	 	   }
          	 	 }
							}
							sigtree->Fill();
						}//MET Filter
					}// Z mass Filter
				}//dR filter
			}// ele + pho candidate
			
			// Filling ele->photon fake tree : proxyTree
			if(!hasPho){
			for(unsigned ip(0); ip < proxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = proxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyLepCollection[ie];
					double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dRlepphoton>0.8){
						if(((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(METFilter)){

								float proxy_deltaPhi = DeltaPhi(proxyEle->getPhi(), METPhi);
								float proxy_MT = sqrt(2*MET*proxyEle->getCalibPt()*(1-std::cos(proxy_deltaPhi)));
								proxyphoEt = proxyPho->getCalibEt();
								proxyphoEta= proxyPho->getEta();
								proxyphoPhi= proxyPho->getPhi();
								proxylepPt = proxyEle->getCalibPt();
								proxylepEta= proxyEle->getEta();
								proxylepPhi= proxyEle->getPhi();
								proxysigMT = proxy_MT;
								proxysigMET= MET;
								proxysigMETPhi = METPhi;
								proxydPhiLepMET = proxy_deltaPhi; 
								proxynVertex = nVtx; 
								proxydRPhoLep= dRlepphoton;
								proxyFSRVeto = proxyPhoFSRVeto[ip];
								if(hasTrail){
									proxytrailPt = trailLep->getPt();
									proxytrailEta = trailLep->getEta();
									proxytrailPhi = trailLep->getPhi();
								}
								else{
									proxytrailPt  = 0; 
									proxytrailEta = 0;
									proxytrailPhi = 0;
								}

								proxynJet = 0;
								proxyHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyEle->getEta(),proxyEle->getPhi()) <= 0.4)continue;
									proxynJet += 1;
									proxyHT += itJet->getPt();
								}
								proxytree->Fill();


                                                       }//MET Filter
						}// Z mass Filter
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection
			}

			// Filling jet->photon fake tree; jetTree
			if(!hasPho){
				for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
					for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoEle>::iterator jetEle = proxyLepCollection[ie];
					std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
					double dRlepphoton = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetEle->getEta(), jetEle->getPhi());
					if(dRlepphoton>0.8){
						if(((jetPho->getCalibP4()+jetEle->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(METFilter)){


								float jet_deltaPhi = DeltaPhi(jetEle->getPhi(), METPhi);
								float jet_MT = sqrt(2*MET*jetEle->getCalibPt()*(1-std::cos(jet_deltaPhi)));
								jetphoEt = jetPho->getCalibEt();
								jetphoEta= jetPho->getEta();
								jetphoPhi= jetPho->getPhi();
								jetphoChIso = jetPho->getChIso();
								jetphoSigma = jetPho->getSigma();
								jetlepPt = jetEle->getCalibPt();
								jetlepEta= jetEle->getEta();
								jetlepPhi= jetEle->getPhi();
								jetsigMT = jet_MT;
								jetsigMET= MET;
								jetsigMETPhi = METPhi;
								jetdPhiLepMET = jet_deltaPhi; 
								jetnVertex = nVtx; 
								jetdRPhoLep= dRlepphoton;
								if(hasTrail){
									jettrailPt = trailLep->getPt();
									jettrailEta = trailLep->getEta();
									jettrailPhi = trailLep->getPhi();
								}
								else{
									jettrailPt  = 0; 
									jettrailEta = 0;
									jettrailPhi = 0;
								}

								jetnJet = 0;
								jetHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetEle->getEta(),jetEle->getPhi()) <= 0.4)continue;
									jetnJet += 1;
									jetHT += itJet->getPt();	
								}
								jettree->Fill();
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				}} // loop on pho collection
			}

			// Filling jet->lep fake tree ; fakeLeptree
			if(hasPho && !hasLep){
				std::vector<recoPhoton>::iterator fakeLepPho = signalPho;
				for(unsigned ip(0); ip < fakeLepCollection.size(); ip++){
					std::vector<recoEle>::iterator fakeLep = fakeLepCollection[ip];
					double dRlepphoton = DeltaR(fakeLepPho->getEta(), fakeLepPho->getPhi(), fakeLep->getEta(), fakeLep->getPhi());
					if(dRlepphoton>0.8){
						if(((fakeLepPho->getCalibP4()+fakeLep->getCalibP4()).M() - 91.188) > 10.0){
							if(raw.passMETFilter(METFilter)){

								fakeLepdRPhoLep = 3;
								for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
									double dR = DeltaR(itpho->getEta(), itpho->getPhi(), fakeLep->getEta(), fakeLep->getPhi());
									if(dR < 0.05 && fabs(itpho->getEt() - fakeLep->getPt())/fakeLep->getPt() < 0.1)continue;
									if(dR < fakeLepdRPhoLep)fakeLepdRPhoLep = dR;
								}	
		
								float fakeLep_deltaPhi = DeltaPhi(fakeLep->getPhi(), METPhi);
								float fakeLep_MT = sqrt(2*MET*fakeLep->getCalibPt()*(1-std::cos(fakeLep_deltaPhi)));
								fakeLepphoEt = fakeLepPho->getCalibEt();
								fakeLepphoEta= fakeLepPho->getEta();
								fakeLepphoPhi= fakeLepPho->getPhi();
								fakeLepPt = fakeLep->getCalibPt();
								fakeLepEta= fakeLep->getEta();
								fakeLepPhi= fakeLep->getPhi();
								fakeLepMiniIso = fakeLep->getMiniIso();
								if(fakeLep->isFakeProxy())fakeLepIsStandardProxy = 1;
								else fakeLepIsStandardProxy = 0;
								fakeLepsigMT = fakeLep_MT;
								fakeLepsigMET= MET;
								fakeLepsigMETPhi = METPhi;
								fakeLepdPhiLepMET = fakeLep_deltaPhi; 
								fakeLepnVertex = nVtx; 
								fakeLepSigma = fakeLep->getSigma();
								fakeLepdEta  = fabs(fakeLep->getdEtaIn());
								fakeLepdPhi  = fabs(fakeLep->getdPhiIn());
								fakeLepIndex = ip;
								if(hasTrail){
									fakeLeptrailPt = trailLep->getPt();
									fakeLeptrailEta = trailLep->getEta();
									fakeLeptrailPhi = trailLep->getPhi();
								}
								else{
									fakeLeptrailPt  = 0; 
									fakeLeptrailEta = 0;
									fakeLeptrailPhi = 0;
								}

								fakeLepnJet = 0;
								fakeLepHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLep->getEta(),  fakeLep->getPhi()) <= 0.4)continue;
									fakeLepnJet += 1;
									fakeLepHT += itJet->getPt();
								}	
								fakeLeptree->Fill();
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}

		// Filling hadron tree
		if(hasHadronPho || hadeleproxyPhoCollection.size() > 0){
			hadron_phoEt = 0;
			hadron_phoEta = 0;
			hadron_phoPhi = 0;
			hadron_phoSigma = 0;
			hadron_phoChIso = 0;
			hadron_sigMT = 0;
			hadron_sigMET = 0;
			hadron_sigMETPhi = 0;
			hadron_dPhiLepMET = 0;
			hadron_nVertex = 0;
			hadron_HT = 0;
			hadron_nJet = 0;
  		hadron_eleproxyEt.clear();
			hadron_eleproxyEta.clear();
			hadron_eleproxyPhi.clear();
			hadron_eleproxySigma.clear();
			hadron_eleproxyChIso.clear();
			hadron_eleproxynVertex.clear();

			if(hasHadronPho && hasLep){
				double DeltaPhoLep = DeltaR(hadronPho->getEta(), hadronPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				double DoubleMass  = (hadronPho->getP4()+signalLep->getP4()).M();
				if(DeltaPhoLep > 0.8 && DoubleMass > 90 && (DoubleMass - 91.188) > 10.0){
				if(raw.passMETFilter(METFilter)){
					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
					hadron_phoEt = hadronPho->getCalibEt();
					hadron_phoEta = hadronPho->getEta();
					hadron_phoPhi = hadronPho->getPhi();
					hadron_phoSigma = hadronPho->getSigma();
					hadron_phoChIso = hadronPho->getChIso();
					hadron_sigMT = MT;
					hadron_sigMET = MET;
					hadron_sigMETPhi = METPhi;
					hadron_dPhiLepMET = deltaPhi;
					hadron_nVertex = nVtx;
					hadron_HT = 0;
					hadron_nJet = jetNumber;
				}}
			}
			for(unsigned ip(0); ip < hadeleproxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = hadeleproxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyLepCollection[ie];
					double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dRlepphoton>0.8 && ((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0 && raw.passMETFilter(METFilter)){
						hadron_eleproxyEt.push_back(proxyPho->getCalibEt());
						hadron_eleproxyEta.push_back(proxyPho->getEta());
						hadron_eleproxyPhi.push_back(proxyPho->getPhi());
						hadron_eleproxySigma.push_back(proxyPho->getSigma());
						hadron_eleproxyChIso.push_back(proxyPho->getChIso());
						hadron_eleproxynVertex.push_back(nVtx); 
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection

      hadron_mcPID.clear();
      hadron_mcEta.clear();
      hadron_mcPhi.clear();
      hadron_mcPt.clear();
      hadron_mcMomPID.clear();
			if(isMC){
       for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
         if(itMC->getPt() < 10.0)continue;
         mcPID.push_back(itMC->getPID());
         mcMomPID.push_back(itMC->getMomPID());
         mcGMomPID.push_back(itMC->getGMomPID());
         mcEta.push_back(itMC->getEta());
         mcPhi.push_back(itMC->getPhi());
         mcPt.push_back(itMC->getEt());
       }
			}

			hadrontree->Fill();
		}
	
	}//loop on  events
  cout<< sigtree->GetEntries()<<endl;
  cout<< proxytree->GetEntries()<<endl;
  cout<< jettree->GetEntries()<<endl;
  cout<< fakeLeptree->GetEntries()<<endl;
  cout<< hadrontree->GetEntries()<<endl;

  p_eventcount->GetXaxis()->SetBinLabel(1,"nTotal");
  p_eventcount->GetXaxis()->SetBinLabel(2,"npassHLT");
  p_eventcount->GetXaxis()->SetBinLabel(3,"npassPho");
  p_eventcount->GetXaxis()->SetBinLabel(4,"npassLep");
  p_eventcount->GetXaxis()->SetBinLabel(5,"npassdR");
  p_eventcount->GetXaxis()->SetBinLabel(6,"npassZ");
  p_eventcount->GetXaxis()->SetBinLabel(7,"npassMETFilter");
  p_eventcount->Fill(0.5, nTotal);
  p_eventcount->Fill(1.5, npassHLT);
  p_eventcount->Fill(2.5, npassPho);
  p_eventcount->Fill(3.5, npassLep);
  p_eventcount->Fill(4.5, npassdR);
  p_eventcount->Fill(5.5, npassZ);
  p_eventcount->Fill(6.5, npassMETFilter);
	if(RunYear==2018) logfile << "pass HEM cut:  " << passHEM*100/nEvts<<endl;
	
	outputfile->Write();
	outputfile->Close();
	logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_eg(atoi(argv[1]),argv[2]);
    return 0;
}
