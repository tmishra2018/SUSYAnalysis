#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>
#include "TStopwatch.h"
#include "TString.h"
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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_cuts.h"
#include "../../include/analysis_jet.h"


void access_hist_Data()
{
	gSystem->Load("../../lib/libAnaClasses.so");
// /eos/uscms/store/user/msun/Signal/SMS-TChiWG_TuneCUETP8M1_RunIISummer16MiniAODv2.root
// /eos/uscms/store/user/msun/Signal/SMS-T5WG_TuneCUETP8M1_RunIISummer16MiniAOD.root
  int   nJet_(0);
  float	HT_(0);
  TFile f("/eos/uscms/store/user/msun/Signal/SMS-TChiWG_TuneCUETP8M1_RunIISummer16MiniAODv2.root"); 
  TTree *es = (TTree*)f.Get("ggNtuplizer/EventTree");
  RunType datatype(MC);
  std::ostringstream outputname;
  outputname << "resTree_TChiWG_2016.root"; 
  int SUSYtype(-1);
        if(outputname.str().find("T5WG") != std::string::npos){
                std::cout << "T5WG Model !" << std::endl;
                SUSYtype = 5;
        }
  else if(outputname.str().find("TChi") != std::string::npos){
                std::cout << "TChiWG Model !" << std::endl;
                SUSYtype = 1;
        }

  float met(0);
  float metPhi(0);

  TH1F* hnjets=new TH1F("hnjets","",100,0,100);
  TH1F* hht=new TH1F("hht","",64,300,2000);
  TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoEle> Ele;
  std::vector<recoMuon>   Muon;
  std::vector<recoJet>   JetCollection;
  const unsigned nEvts = es->GetEntries();
  //const unsigned nEvts = es->GetEntries()/10;
  std::cout << "total event : " << nEvts << std::endl;
  	for(unsigned ievt(0); ievt<nEvts; ++ievt){
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
		std::vector<mcData>::iterator genPho = MCData.end();
		std::vector<mcData>::iterator genEle = MCData.end();
		for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
                        if(!itJet->passSignalSelection())continue;
                        //if(genPho != MCData.end())
                        //        if(DeltaR(itJet->getEta(), itJet->getPhi(), genPho->getEta(),genPho->getPhi()) <= AK4Cone)continue;
                        //if(genEle != MCData.end())
                        //        if(DeltaR(itJet->getEta(), itJet->getPhi(), genEle->getEta(),genEle->getPhi()) <= AK4Cone)continue;
                        nJet_ += 1;
                        HT_ += itJet->getPt();
                }
		hht->Fill(HT_);
		hnjets->Fill(nJet_);
 	}
	hht->Write();
 	hnjets->Write();
 	outputfile->Close();

/*
  TH1F* hnjets=new TH1F("hnjets","",40,-5,35);
  TH1F* hnelectrons=new TH1F("hnelectrons","",15,-5,10);
  TH1F* hnmuons=new TH1F("hnmuons","",15,-5,10);
  TH1F* hnphotons=new TH1F("hnphotons","",25,-5,20);

  TH1F* hphotonpt = new TH1F("hphotonpt","",100,0,1200);
  TH1F* helectronpt = new TH1F("helectronpt","",100,0,800);
  TH1F* hmuonpt = new TH1F("hmuonpt","",100,0,500);
  TH1F* hjetpt = new TH1F("hjetpt","",100,0,1400);

  TH1F* hphotonMass = new TH1F("hphotonMass","",100,-0.00008,0.00008);
  TH1F* helectronMass = new TH1F("helectronMass","",100,-0.4,0.4);
  TH1F* hjetMass = new TH1F("hjetMass","",100,0.0,140);
  TH1F* hmuonMass = new TH1F("hmuonMass","",100,0.1056,0.1058);

  Int_t nentries1 = (Int_t)t->GetEntries();

   for(Int_t i=0; i<nentries1; i++)
   {                                                                                                                                              
                  t->GetEntry(i);
		  if()
                  hnjets->Fill(nJet);
                  hnelectrons->Fill(nElectron);
                  hnmuons->Fill(nMuon);
                  hnphotons->Fill(nPhoton);
                  hphotonpt->Fill(Photon_pt[0]);
                  helectronpt->Fill(Electron_pt[0]);
                  hmuonpt->Fill(Muon_pt[0]);
                  hjetpt->Fill(Jet_pt[0]);
                  hphotonMass->Fill(Photon_mass[0]);
                  helectronMass->Fill(Electron_mass[0]);
                  hjetMass->Fill(Jet_mass[0]);
                  hmuonMass->Fill(Muon_mass[0]);
   }

   TFile *file = new TFile(Form("%s_%d_slim.root",channel,RunYear), "RECREATE");
   hnjets->Write();
   hnelectrons->Write();
   hnmuons->Write();
   hnphotons->Write();
   hphotonpt->Write();
   helectronpt->Write();
   hmuonpt->Write();
   hjetpt->Write();
   hphotonMass->Write();
   helectronMass->Write();
   hjetMass->Write();
   hmuonMass->Write();

   f.Close();
   file->Close();*/

}


