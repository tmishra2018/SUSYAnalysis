// g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_eg.C -o analysis_eg.exe `root-config --libs`
#include "../include/analysis_commoncode.h"
#include "../include/analysis_cuts.h"

void analysis_eg_check_fakeLepton_tree(int RunYear, const char *Era){//main

  gSystem->Load("/uscms/home/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s_fakeLep.log",RunYear,Era),ios::trunc);

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

  const unsigned nEvts = es->GetEntries();
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_"<<RunYear<<Era<<".root" << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s_check_fakeLep3.root",RunYear,Era),"RECREATE");
  outputfile->cd();
	TH1D *p_METFilter = new TH1D("p_METFilter","",12,-2,10);	
	TH1D *p_invmass = new TH1D("p_invmass","",200,0,200);	

	int nBJet(0);
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

  std::cout << "Total evetns : " << nEvts << std::endl;
  logfile << "Total evetns : " << nEvts << std::endl;
	for (unsigned ievt(0); ievt<10000000; ++ievt){//loop on entries
	//for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

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
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= ElectronVetoCone)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < GeneralCone && ie->getEt()>2.0)photonFSRVeto=false;
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
				cout<<itEle->getCutValueHoverE()<<endl;
				if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
			//	if(itEle->isFakeProxyOLD())fakeLepCollection.push_back(itEle);	
				if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;

				if(!itEle->passHLTSelection())continue;
				//if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
				//if(itEle->isLooseFakeProxy())fakeLepCollection.push_back(itEle);//Loose the proxy definition	
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

	
	}//loop on  events
	cout<<1.*fakeLeptree->GetEntries()<<endl;
	outputfile->Write();
	outputfile->Close();
	logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_eg_check_fakeLepton_tree(atoi(argv[1]),argv[2]);
    return 0;
}
