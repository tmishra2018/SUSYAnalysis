// g++ `root-config --cflags` ../../../lib/libAnaClasses.so skimPreselection_egMC.C -o skimPreselection_egMC.exe `root-config --libs`

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

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
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"


void skimPreselection_egMC(int RunYear, const char *Sample){//main 

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/signalMC/egsignal_%s_%d.log",Sample,RunYear)); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;
  //logfile << "Loose the proxy definition: no upper bounds for photon; LooseFakeProxy for electron" << std::endl;

  RunType datatype;
  if(RunYear==2016) datatype = MCDoubleEG2016;
  if(RunYear==2017) datatype = MCDoubleEG2017;
  if(RunYear==2018) datatype = MCDoubleEG2018;
	bool  isMC(false);
	//if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
	//if(datatype == MC || datatype == MCDoubleEG2017 || datatype == MCMuonEG2017||  datatype == MCSingleElectron2017 || datatype == MCSingleMuon2017||  datatype == MCDoubleMuon2017 || datatype == MCMET2017)isMC=true;
	//if(datatype == MC || datatype == MCDoubleEG2018 || datatype == MCMuonEG2018||  datatype == MCSingleElectron2018 || datatype == MCSingleMuon2018||  datatype == MCDoubleMuon2018 || datatype == MCMET2018)isMC=true;

  TChain* es = new TChain("ggNtuplizer/EventTree");
	
	char* inputfile = new char[300];
  	sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d.root",Sample,Sample,RunYear);
  	es->Add(inputfile);

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;

  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile* outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/signalMC/egsignal_%s_%d.root",Sample,RunYear),"RECREATE");
  outputfile->cd();
  
    int mcType;
  if(strstr(inputfile, "WGToLNuG") != NULL){
                std::cout << "WGToLNuG sample !" << std::endl;
                mcType = MCType::WGJetInclusive;
  }
  else if(strstr(inputfile, "WGJet40") != NULL){
                std::cout << "WGJet40 sample !" << std::endl;
                mcType = MCType::WGJet40;
  }
  else if(strstr(inputfile, "WGJet130") != NULL){
                std::cout << "WGJet130 sample !" << std::endl;
                mcType = MCType::WGJet130;
  }
  else if(strstr(inputfile, "ZGToLLG") != NULL){
                std::cout << "ZGInclusive sample !" << std::endl;
                mcType = MCType::ZGInclusive;
  }
  else if(strstr(inputfile, "DYJetsToLL") != NULL){
                std::cout << "DYJetsToLL sample !" << std::endl;
                mcType = MCType::DYLL50;
  }
  else if(strstr(inputfile, "TTGJets") != NULL){
                std::cout << "TTGJets sample !" << std::endl;
                mcType = MCType::TTG;
  }
  else if(strstr(inputfile, "WWG") != NULL){
                std::cout << "WWG sample !" << std::endl;
                mcType = MCType::WWG;
  }
  else if(strstr(inputfile, "WZG") != NULL){
                std::cout << "WZG sample !" << std::endl;
                mcType = MCType::WZG;
  }
  else if(strstr(inputfile, "TTJets") != NULL){
                std::cout << "TTJets sample !" << std::endl;
                mcType = MCType::TT;
  }
  else if(strstr(inputfile, "WW") != NULL){
                std::cout << "WW sample !" << std::endl;
                mcType = MCType::WW;
  }
  else if(strstr(inputfile, "WZ") != NULL){
                std::cout << "WZ sample !" << std::endl;
                mcType = MCType::WZ;
  }
  else if(strstr(inputfile, "GJet") != NULL){
                std::cout << "GJet sample !" << std::endl;
                mcType = MCType::GJet;
  }
  else {
                std::cout << "not specific MC !" << std::endl;
                mcType = MCType::NOMC;
  }
  Double_t  L1ECALPrefire;
  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();
  float lumiWeight = getEvtWeight(RunYear,crosssection, ntotalevent);

  cout<<"crosssection = "<<crosssection<<endl;
  cout<<"ntotalevent = "<<ntotalevent<<endl;
  cout<<"lumiWeight = "<<lumiWeight<<endl;

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

  sigtree->Branch("crosssection",&crosssection);
  sigtree->Branch("ntotalevent", &ntotalevent);
  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("L1ECALPrefire",     &L1ECALPrefire);
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
	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

		if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if (ievt%100000==0) logfile  << " -- Processing event " << ievt << std::endl;

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
			L1ECALPrefire=raw.L1ECALPrefire;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nEle < 1 || raw.nPho <1)continue;


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
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool photonFSRVeto(true);
				bool eleFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getEt()>2.0)photonFSRVeto=false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.02)eleFSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0){
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
						if( (itpho->getSigma()< 0.02 && itpho->isEB()) || (itpho->getSigma()< 0.04 && itpho->isEE()) )jetPhoCollection.push_back(itpho);
					}
				}
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
						proxyPhoCollection.push_back(itpho);
						proxyPhoFSRVeto.push_back( eleFSRVeto? 1:0);
				}
			}





			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Lepton                                              **********************************************/
			bool hasLep(false);
			std::vector<recoEle>::iterator signalLep = Ele.begin();
			bool hasTrail(false);
			std::vector<recoEle>::iterator trailLep = Ele.begin();
			std::vector< std::vector<recoEle>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoEle>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			std::vector< std::vector<recoEle>::iterator > miniisoLep;
			miniisoLep.clear();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if(itEle->isMedium() && itEle->getPt() > 15 && itEle->getMiniIso() < 0.2)miniisoLep.push_back(itEle);
				if(itEle->getCalibPt() < 25)continue;
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;

				if(!itEle->passHLTSelection())continue;
				if(itEle->isFakeProxy())fakeLepCollection.push_back(itEle);	
				//if(itEle->isLooseFakeProxy())fakeLepCollection.push_back(itEle);//Loose the proxy definition	
				if(itEle->passSignalSelection()){
					proxyLepCollection.push_back(itEle);
					if(hasLep && !hasTrail){
						hasTrail = true;
						trailLep = itEle;
					}
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itEle;
					}
				}
			}

			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi()); 
				if(dRlepphoton > 0.8){
					npassdR+=1;
					if(fabs((signalPho->getCalibP4()+signalLep->getCalibP4()).M() - 91.188) > 10.0){

						npassZ+=1;
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
	 
	
	}//loop on  events

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

	outputfile->Write();
	outputfile->Close();
	logfile.close();
}

int main(int argc, char** argv)
{
    if(argc < 2)
      cout << "You have to provide two arguments!!\n";
    skimPreselection_egMC(atoi(argv[1]),argv[2]);
    return 0;
}
