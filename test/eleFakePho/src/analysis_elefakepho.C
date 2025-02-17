// Run using g++ `root-config --cflags` analysis_elefakepho.C -o analysis_elefakepho.exe `root-config --libs`
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
#include "TRandom3.h"
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
#include "../../../src/analysis_rawData.cc"
#include "../../../src/analysis_ele.cc"
#include "../../../src/analysis_photon.cc"

void analysis_elefakepho(int RunYear, const char *Era){//main

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/logs/plot_elefakepho_DYTnP_dR05_DY_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_elefakepho()" << std::endl;

  RunType datatype(MC);                              // Run Type
  TFile *f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DYJetsToLL/DYJetsToLL_%d%s.root",RunYear,Era));
  TTree *es =(TTree*)f->Get("ggNtuplizer/EventTree");

  TFile *output = TFile::Open(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_DYTnP_dR05_%d%s.root",RunYear,Era),"RECREATE");
  output->cd();

  int   tracks(0);
  int   nVertex(0); 
  int   mcType = MCType::DYLL50;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;

  TTree *etree = new TTree("FakeRateTree","FakeRateTree");
  float tagEta_bothcount;
  float tagPhi_bothcount;
  float tagPt_bothcount;
  float tagUncalibPt_bothcount;
  float probeEta_bothcount;
  float probePhi_bothcount;
  float probePt_bothcount;
  float probeUncalibPt_bothcount;
  float invmass_bothcount;
  float invmassUncalib_bothcount;
  bool  vetovalue_bothcount;
	bool  FSRveto_bothcount;
  std::vector<int> mcPID_bothcount;
  std::vector<float> mcEta_bothcount;
  std::vector<float> mcPhi_bothcount;
  std::vector<float> mcPt_bothcount;
  std::vector<int> mcMomPID_bothcount;
  std::vector<int> mcGMomPID_bothcount;

  etree->Branch("tracks",              &tracks);
  etree->Branch("nVertex",             &nVertex);
  etree->Branch("tagEta",              &tagEta_bothcount);
  etree->Branch("tagPhi",              &tagPhi_bothcount);
  etree->Branch("tagPt",               &tagPt_bothcount);
  etree->Branch("tagUncalibPt",        &tagUncalibPt_bothcount);
  etree->Branch("probeEta",            &probeEta_bothcount);
  etree->Branch("probePhi",            &probePhi_bothcount);
  etree->Branch("probePt",             &probePt_bothcount);
  etree->Branch("probeUncalibPt",      &probeUncalibPt_bothcount);
  etree->Branch("invmass",             &invmass_bothcount);
  etree->Branch("invmassUncalib",      &invmassUncalib_bothcount);
  etree->Branch("vetovalue",           &vetovalue_bothcount);
	etree->Branch("FSRveto",       &FSRveto_bothcount);
  etree->Branch("mcPID",			   &mcPID_bothcount);
  etree->Branch("mcEta",			   &mcEta_bothcount);
  etree->Branch("mcPhi",			   &mcPhi_bothcount);
  etree->Branch("mcPt",				   &mcPt_bothcount);
  etree->Branch("mcMomPID",			   &mcMomPID_bothcount);
  etree->Branch("mcGMomPID",		   &mcGMomPID_bothcount);
  etree->Branch("mcType",			   &mcType);

  TTree *rantree = new TTree("FakeRateRandomTree","FakeRateRandomTree");
  float tagEta_random;
  float tagPhi_random;
  float tagPt_random;
  float tagUncalibPt_random;
  float probeEta_random;
  float probePhi_random;
  float probePt_random;
  float probeUncalibPt_random;
  float invmass_random;
  float invmassUncalib_random;
  bool  vetovalue_random;
	bool  FSRveto_random;
  std::vector<int> mcPID_random;
  std::vector<float> mcEta_random;
  std::vector<float> mcPhi_random;
  std::vector<float> mcPt_random;
  std::vector<int> mcMomPID_random;
  std::vector<int> mcGMomPID_random;

  rantree->Branch("tracks",              &tracks);
  rantree->Branch("nVertex",             &nVertex);
  rantree->Branch("tagEta",              &tagEta_random);
  rantree->Branch("tagPhi",              &tagPhi_random);
  rantree->Branch("tagPt",               &tagPt_random);
  rantree->Branch("tagUncalibPt",        &tagUncalibPt_random);
  rantree->Branch("probeEta",            &probeEta_random);
  rantree->Branch("probePhi",            &probePhi_random);
  rantree->Branch("probePt",             &probePt_random);
  rantree->Branch("probeUncalibPt",      &probeUncalibPt_random);
  rantree->Branch("invmass",             &invmass_random);
  rantree->Branch("invmassUncalib",      &invmassUncalib_random);
  rantree->Branch("vetovalue",           &vetovalue_random);
	rantree->Branch("FSRveto",         &FSRveto_random);
  rantree->Branch("mcPID",			     &mcPID_random);
  rantree->Branch("mcEta",			     &mcEta_random);
  rantree->Branch("mcPhi",			     &mcPhi_random);
  rantree->Branch("mcPt",				 &mcPt_random);
  rantree->Branch("mcMomPID",			 &mcMomPID_random);
  rantree->Branch("mcGMomPID",		     &mcGMomPID_random);
  rantree->Branch("mcType",			     &mcType);

  const unsigned nEvts = es->GetEntries();
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << "/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_DYTnP_dR05_"<<RunYear<<Era<<".root"<< std::endl;
  
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
  int   ntrks(0);
  int   nvtx(0);

  TRandom3 ran(0);
  int pTcut = 0;
  if (RunYear==2016) pTcut = 30;
  if (RunYear==2017) pTcut = 38;
  if (RunYear==2018) pTcut = 35;

  int passPixelIssue(0);
  std::cout << "total: " << nEvts << std::endl;
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
        if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
        raw.GetData(es, ievt);
        MCData.clear();
        Photon.clear();
        Muon.clear();
        Ele.clear();
        if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
        for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
        for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
        MET = raw.pfMET;
        METPhi = raw.pfMETPhi;
        nvtx = raw.nVtx;

        tracks = ntrks;
        nVertex = nvtx;
	
	if(RunYear==2017 && !passPixelIssue17(raw)) continue;
        if(RunYear==2018 && !passPixelIssue18(raw)) continue;
        passPixelIssue++;

        if(MET > 70.0)continue;
        //if(!raw.passHLT())continue;
				if(RunYear==2016 && ((raw.HLTEleMuX >> 4) &1) ==0)continue;  //HLT_Ele27_WPTight_Gsf_v
				if(RunYear==2017 && ((raw.HLTEleMuX >> 3) &1) ==0)continue;  //HLT_Ele35_WPTight_Gsf_v
				if(RunYear==2018 && ((raw.HLTEleMuX >> 55) &1) ==0)continue; //HLT_Ele32_WPTight_Gsf_v
        std::vector<std::vector<recoEle>::iterator> ElectronCollection;
        ElectronCollection.clear();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
           if(itEle->getCalibEt() < pTcut || fabs(itEle->getEta())>2.1)continue;                              // Tag electron selection
           //if(!itEle->passHLTSelection())continue;
					 if(RunYear==2016 && !itEle->fireTrgs(12))continue;  //HLT_Ele27_WPTight_Gsf_v
					 if(RunYear==2017 && !itEle->fireTrgs(46))continue;  //HLT_Ele35_WPTight_Gsf_v
					 if(RunYear==2018 && !itEle->fireTrgs(13))continue;  //HLT_Ele32_WPTight_Gsf_v
	   if(itEle->passSignalSelection())ElectronCollection.push_back(itEle); // Tag electron selection
		// pt > 30 GeV, medium ID, eta < 2.1 electron
        }
        
        std::vector<std::vector<recoPhoton>::iterator> signalPho;
        std::vector<bool> PhoPixelVeto;
        std::vector<bool> PhoEleVeto;
				std::vector<bool> PhoFSRVeto;
        signalPho.clear(); 
        PhoPixelVeto.clear();
        PhoEleVeto.clear();
				PhoFSRVeto.clear();
        if(ElectronCollection.size() > 0){
	      for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
	        if(itpho->getCalibEt() < 30)continue;                            // Photon pt and loose selection
	        if(itpho->isLoose()){
              bool PixelVeto = itpho->PixelSeed()==0? true: false;               // check whether pixel seed or not
              bool GSFveto(true);
		      bool FSRVeto(true);
			  for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				 //if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.05)GSFveto = false;
				 // electron match to a photon
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibEt()>2.0)FSRVeto=false;
			  }
			  for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;
					// FSR rejection; no events with lepton near photon.
				        signalPho.push_back(itpho);              // signal photons
					PhoPixelVeto.push_back(PixelVeto);
					PhoEleVeto.push_back(GSFveto);
					PhoFSRVeto.push_back(FSRVeto);
            }
	      }
        }
        
        int nTagEle = ElectronCollection.size();
	    int randomTag(0);
	    if(nTagEle>1)randomTag = ran.Integer(nTagEle);
	    for(int iTag(0); iTag < nTagEle; iTag++){                             // for # of tags (electrons)
	      for(unsigned iProbe(0); iProbe < signalPho.size(); iProbe++){       // for # of probes (signal photons)
	        double dRTagProbe = DeltaR(ElectronCollection[iTag]->getEta(), ElectronCollection[iTag]->getPhi(), signalPho[iProbe]->getEta(), signalPho[iProbe]->getPhi());
	        double dETagProbe = fabs(ElectronCollection[iTag]->getCalibEt() - signalPho[iProbe]->getCalibEt())/signalPho[iProbe]->getCalibEt();
	        double InvMass = (ElectronCollection[iTag]->getCalibP4() + signalPho[iProbe]->getCalibP4()).M();
	        if(dRTagProbe < 0.05 && dETagProbe < 0.1)continue;
	        if(InvMass < 20 || InvMass > 160)continue;
              
	        if(iTag == randomTag){
			  tagEta_random = ElectronCollection[iTag]->getEta();
			  tagPhi_random = ElectronCollection[iTag]->getPhi();
			  tagPt_random = ElectronCollection[iTag]->getCalibPt();
			  tagUncalibPt_random = ElectronCollection[iTag]->getPt();
			  probeEta_random = signalPho[iProbe]->getEta();
			  probePhi_random = signalPho[iProbe]->getPhi();
			  probePt_random = signalPho[iProbe]->getCalibEt();
			  probeUncalibPt_random = signalPho[iProbe]->getEt();
			  invmass_random = InvMass;
			  invmassUncalib_random = (ElectronCollection[iTag]->getP4() + signalPho[iProbe]->getP4()).M();
			  vetovalue_random = (PhoPixelVeto[iProbe] && PhoEleVeto[iProbe]);
				FSRveto_random = PhoFSRVeto[iProbe];
              if(datatype == MC){
				mcPID_random.clear();
				mcEta_random.clear();
				mcPhi_random.clear();
				mcPt_random.clear();
				mcMomPID_random.clear();
				mcGMomPID_random.clear();

				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				  if(itMC->getEt() < 5)continue;
				  float mcdR = DeltaR(signalPho[iProbe]->getEta(),signalPho[iProbe]->getPhi(), itMC->getEta(), itMC->getPhi());
				  float mcdRe= DeltaR(ElectronCollection[iTag]->getEta(),ElectronCollection[iTag]->getPhi(),itMC->getEta(), itMC->getPhi());
				  if(mcdR < 0.3 || mcdRe < 0.3){
					mcPID_random.push_back(itMC->getPID());
					mcMomPID_random.push_back(itMC->getMomPID());
					mcGMomPID_random.push_back(itMC->getGMomPID());
					mcEta_random.push_back(itMC->getEta());      
					mcPhi_random.push_back(itMC->getPhi());
					mcPt_random.push_back(itMC->getEt());
				  }
				}
              }
		
              rantree->Fill();
            }
			tagEta_bothcount = ElectronCollection[iTag]->getEta();
			tagPhi_bothcount = ElectronCollection[iTag]->getPhi();
			tagPt_bothcount = ElectronCollection[iTag]->getCalibPt();
			tagUncalibPt_bothcount = ElectronCollection[iTag]->getPt();
			probeEta_bothcount = signalPho[iProbe]->getEta();
			probePhi_bothcount = signalPho[iProbe]->getPhi();
			probePt_bothcount = signalPho[iProbe]->getCalibEt();
			probeUncalibPt_bothcount = signalPho[iProbe]->getEt();
			invmass_bothcount = InvMass;
			invmassUncalib_bothcount = (ElectronCollection[iTag]->getP4() + signalPho[iProbe]->getP4()).M();
			vetovalue_bothcount = (PhoPixelVeto[iProbe] && PhoEleVeto[iProbe]);
			FSRveto_bothcount = PhoFSRVeto[iProbe];

            if(vetovalue_bothcount == false){
              bool isTag(false);
              for(int iTag(0); iTag < nTagEle; iTag++){
                double testdR = DeltaR(ElectronCollection[iTag]->getEta(), ElectronCollection[iTag]->getPhi(), signalPho[iProbe]->getEta(), signalPho[iProbe]->getPhi());
                if(testdR < 0.05)isTag = true;
              }
            }

			if(datatype == MC){
			  mcPID_bothcount.clear();
			  mcEta_bothcount.clear();
			  mcPhi_bothcount.clear();
			  mcPt_bothcount.clear();
			  mcMomPID_bothcount.clear();
			  mcGMomPID_bothcount.clear();

			  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 5)continue;
				float mcdR = DeltaR(signalPho[iProbe]->getEta(),signalPho[iProbe]->getPhi(), itMC->getEta(), itMC->getPhi());
				float mcdRe= DeltaR(ElectronCollection[iTag]->getEta(),ElectronCollection[iTag]->getPhi(),itMC->getEta(), itMC->getPhi());
				if(mcdR < 0.3 || mcdRe < 0.3){
				  mcPID_bothcount.push_back(itMC->getPID());
				  mcMomPID_bothcount.push_back(itMC->getMomPID());
				  mcGMomPID_bothcount.push_back(itMC->getGMomPID());
				  mcEta_bothcount.push_back(itMC->getEta());      
				  mcPhi_bothcount.push_back(itMC->getPhi());
				  mcPt_bothcount.push_back(itMC->getEt());
				}
		  }
			}
            etree->Fill();
	   }
	 }
           

   }//loop on  events
float percent1=100*rantree->GetEntries()/nEvts;
float percent2=100*etree->GetEntries()/nEvts;
logfile << "FakeRateRandomTree events: " << rantree->GetEntries() <<"; "<<percent1<<"\%"<<std::endl;
logfile << "FakeRateTree events: " << etree->GetEntries() <<"; "<<percent2<<"\%"<<std::endl;
cout<< "FakeRateRandomTree events: " << rantree->GetEntries() <<"; "<<percent1<<"\%"<<std::endl;
cout<< "FakeRateTree events: " << etree->GetEntries() <<"; "<<percent2<<"\%"<<std::endl;
if(RunYear==2018) logfile << "pass Pixel veto " << passPixelIssue*100/nEvts<<endl;
if(RunYear==2017) logfile << "pass Pixel veto " << passPixelIssue*100/nEvts<<endl;


output->Write();
output->Close();
logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    analysis_elefakepho(atoi(argv[1]),argv[2]);
    return 0;
}
