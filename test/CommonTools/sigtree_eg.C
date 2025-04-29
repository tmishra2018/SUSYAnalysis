// g++ `root-config --cflags` ../../lib/libAnaClasses.so sigtree_eg.C -o sigtree_eg.exe `root-config --libs`
#include "../../include/analysis_commoncode.h"
#include "../../include/analysis_cuts.h"
bool apply_HEMveto=false;

void sigtree_eg(int RunYear, const char *Era /*for Data*/, bool preVFP, bool isMC, const char *Sample /*for MC*/ ){//main

  gSystem->Load("../../lib/libAnaClasses.so");
  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";


  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/logs/sigtree_eg_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;

  RunType datatype;
  TChain* es = new TChain("ggNtuplizer/EventTree");
  char* inputfile = new char[300];

  if (isMC == false){
  	if(RunYear==2016) datatype = DoubleEG2016;
  	if(RunYear==2017) datatype = DoubleEG2017;
  	if(RunYear==2018) datatype = DoubleEG2018;
  	sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_%d%s.root",RunYear,Era);
  }
  else if(isMC == true){
	if(RunYear==2016) datatype = MCDoubleEG2016;
  	if(RunYear==2017) datatype = MCDoubleEG2017;
  	if(RunYear==2018) datatype = MCDoubleEG2018;
  //	if (strstr(Sample, "DYJetsToLL") != NULL or strstr(Sample, "TTJets") != NULL or strstr(Sample, "WJetsToLNu") != NULL or  strstr(Sample, "T5Wg") != NULL or strstr(Sample, "TChiWG") != NULL)
        	sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  //	else
  //     	sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
  }

  es->Add(inputfile);
  if(RunYear==2018) apply_HEMveto=true; 
  
  const unsigned nEvts = es->GetEntries();
  logfile << "Output file: " << "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_"<<Sample<<"_"<<RunYear<<whichVFP<<".root" << std::endl;

  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_%s_%d%s.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
  outputfile->cd();
  TH1D *p_METFilter = new TH1D("p_METFilter","",12,-2,10);	
  TH1D *p_invmass = new TH1D("p_invmass","",200,0,200);	

  int nBJet(0);
  //************ Signal Tree **********************
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
  Double_t  L1ECALPrefire(1.0);

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
	sigtree->Branch("L1ECALPrefire",     &L1ECALPrefire);
  	sigtree->Branch("mcPID",     &mcPID);
  	sigtree->Branch("mcEta",     &mcEta);
  	sigtree->Branch("mcPhi",     &mcPhi);
  	sigtree->Branch("mcPt",      &mcPt);
  	sigtree->Branch("mcMomPID",  &mcMomPID);
  	sigtree->Branch("mcGMomPID", &mcGMomPID);
  }

//*********** histo list **********************
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

  std::cout << Sample <<"  "<< RunYear<<"    Total events : " << nEvts << std::endl;
  logfile << "Total events : " << nEvts << std::endl;
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
			
			if(isMC) L1ECALPrefire=raw.L1ECALPrefire;
			//std::cout<<ievt <<"  "<< L1ECALPrefire <<std::endl;
			if(RunYear==2018 && !passHEMVeto(0,raw)) continue;
        		passHEM++;

			nTotal+=1;
			//if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;
			// eg combination
			if(raw.nEle < 1 || raw.nPho <1)continue;

			nBJet = 0;
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(itJet->getPt() < 20)continue;
				if(itJet->isBJet())nBJet+=1;
			}

			/*******************************Select Photon**********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< int > proxyPhoFSRVeto;
			// ******* for jet-photon fake rate**************
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < R9EBCut)continue;
		//		if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
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
        		/****************  standard ID ************************************/
				if(!itpho->passSignalSelection())continue;
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

			}

			/************** Select Lepton *********************************/

			bool hasLep(false);
			std::vector<recoEle>::iterator signalLep = Ele.begin();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){

				if(itEle->getCalibPt() < 25)continue;
				if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;
		//		if(!itEle->passHLTSelection())continue;
				if(itEle->passSignalSelection()){
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
						p_METFilter->Fill(raw.failFilterStep(RunYear, METFilter));	
						if(raw.passMETFilter(RunYear, METFilter)){
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
  cout<<endl<<"signal tree entries: " << sigtree->GetEntries()<<endl;

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
	p_METFilter->Write();
	outputfile->Write();
	outputfile->Close();
	logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 5)
      cout << "You have to provide five arguments!!\n";
        bool preVFP = (atoi(argv[3]) == 1);
        bool isMC = (atoi(argv[4]) == 1);

    sigtree_eg(atoi(argv[1]),argv[2],preVFP,isMC,argv[5]);
    return 0;
}
