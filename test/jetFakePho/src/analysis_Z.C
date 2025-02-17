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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"


void analysis_Z(){//main  

  gSystem->Load("../../../lib/libAnaClasses.so");

  char outputname[100] = "plotZ_DoubleMu_eg.root";
  ofstream logfile;
  logfile.open("plotZ_DoubleMu_2016Rereco.log");

  logfile << "analysis_Z()" << std::endl;

  RunType datatype(DoubleMuon2016); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016B_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016C_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016D_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016E_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016F_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016G_SepRereco.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016H_SepRereco_0.root");
  es->Add("root://cmsxrootd.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleMu_Run2016H_SepRereco_1.root");
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016B_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016C_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016D_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016E_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016F_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016G_SepRereco.root" << std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016H_SepRereco_0.root"<< std::endl;
  logfile << "root://cmsxrootd.fnal.gov//store/user/msun/2016ggNtuple/skim-DoubleMu_Run2016H_SepRereco_1.root"<< std::endl;
  

  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

//************ Signal Tree **********************//
  TTree *Zeetree = new TTree("ZeeTree","ZeeTree");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float phoSigma(0);
  float phoChIso(0);
  float invmass(0);
  float dimuonmass(0);
  bool  ismg(0);
  int   nVtx(0);
  
  Zeetree->Branch("phoEt",     &phoEt);
  Zeetree->Branch("phoEta",    &phoEta);
  Zeetree->Branch("phoPhi",    &phoPhi);
  Zeetree->Branch("phoSigma",  &phoSigma);
  Zeetree->Branch("phoChIso",  &phoChIso);
  Zeetree->Branch("invmass",   &invmass);
  Zeetree->Branch("dimuonmass",&dimuonmass);
  Zeetree->Branch("ismg",      &ismg); 
  Zeetree->Branch("nVtx",      &nVtx);
 
  TH1F *p_HLT =new TH1F("p_HLT","",60,0,60);

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "total " << nEvts << std::endl;
  logfile <<   "total " << nEvts << std::endl;
  
  rawData raw(es, datatype);
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  logfile << "RunType: " << datatype << std::endl;

  for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
    if(ievt%10000==0)std::cout << " -- Processing event " << ievt << std::endl;

        raw.GetData(es, ievt);
        Photon.clear();
        Muon.clear();
	    Ele.clear();
        for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
	    for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
        MET = raw.pfMET;
   
		if(!raw.passHLT())continue;
		if(raw.nPho <1 || raw.nMu<2)continue;
		if(MET>70)continue;

        bool foundZ(false);  
        for(std::vector<recoMuon>::iterator itLeadMu = Muon.begin(); itLeadMu!= Muon.end(); ++itLeadMu){
          if(foundZ)break;
          if(!itLeadMu->isMedium() || itLeadMu->getPt() < 10.5 || fabs(itLeadMu->getEta())>2.4)continue;
          for(std::vector<recoMuon>::iterator itMu(itLeadMu+1); itMu!= Muon.end(); itMu++){
            if(!itMu->isMedium()|| itMu->getPt() < 10.5 || fabs(itLeadMu->getEta())>2.4)continue;
            bool sameCharge = ((itLeadMu->isPosi() && itMu->isPosi()) || (!itLeadMu->isPosi() && !itMu->isPosi()));        
            if(sameCharge)continue; 

            float dimuonmass = (itLeadMu->getP4()+itMu->getP4()).M(); 
            if(dimuonmass < 35)continue;
            if(dimuonmass > 80)continue;    

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
              if(!itpho->isEB())continue;
			  bool passID = (itpho->passHoverE(1) && itpho->passNeuIso(1) && itpho->passPhoIso(1));
			  if(itpho->getEt() < 25|| !passID)continue;
			  bool PixelVeto = itpho->PixelSeed()==0? true: false;
			  bool GSFveto(true);
			  bool FSRVeto(true);
			  for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
			  }
			  for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
			  float dR1,dR2,InvMass;
			  dR1 = DeltaR(itpho->getEta(), itpho->getPhi(), itLeadMu->getEta(), itLeadMu->getPhi());
			  dR2 = DeltaR(itpho->getEta(), itpho->getPhi(), itMu->getEta(), itMu->getPhi());
			  InvMass = (itLeadMu->getP4()+itMu->getP4()+ itpho->getP4()).M(); 
			  bool farthestMu(true);
			  //if(dR1>dR2 && itLeadMu->getEt()< 21)farthestMu=false;
			  //else if(dR1<dR2 && itMu->getEt() < 21)farthestMu=false;        
			  if(itLeadMu->getEt()< 21 && itMu->getEt() < 21)farthestMu=false;
              if(!itpho->fireDoubleTrg(5) && !itpho->fireDoubleTrg(6))continue;

			  if(!GSFveto || !PixelVeto)continue; 
              if(!FSRVeto)continue;
			  //if(dimuonmass+InvMass > 180 || !farthestMu)continue;
			  if(!farthestMu)continue;
				if(itpho->getR9() < 0.5)continue;
				   
			  if(dR1 < 0.8 || dR2 < 0.8){
                if(InvMass > 60 && InvMass < 120){
                  for(unsigned iH(0); iH< 60; iH++){
                    if(((raw.HLTEleMuX >> iH) &1)!=0)p_HLT->Fill(iH);
                  }
                } 
				if(InvMass > 60 && InvMass < 120 && raw.passHLT() ){
                  ismg=(itpho->fireL1Trg(12));
				  phoEt = itpho->getCalibEt();
				  phoEta= itpho->getEta();
				  phoPhi= itpho->getPhi();
				  phoSigma=itpho->getSigma();
				  phoChIso=itpho->getChIso();
                  nVtx = raw.nVtx;
				  invmass=InvMass;
				  dimuonmass=dimuonmass;
            
                  Zeetree->Fill();

                  if(InvMass > 80 && InvMass < 100)foundZ=true;
                }
			  }
			}//loop on photons
         }//loop on second muon
       }//loop on leading muon
    }

outputfile->Write();
}


