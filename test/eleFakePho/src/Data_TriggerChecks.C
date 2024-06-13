// Run using g++ `root-config --cflags` Data_TriggerChecks.C -o Data_TriggerChecks.exe `root-config --libs`
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
bool apply_HEMveto=false;
bool apply_L1=false;

void Data_TriggerChecks(int RunYear, const char *Era){//main

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/logs/plot_elefakepho_DataTnP_dR05_Data_%d%s.log",RunYear,Era),ios::trunc);

  logfile << "analysis_elefakepho()" << std::endl;

  RunType datatype;
  if(RunYear==2016) datatype = SingleElectron2016;
  if(RunYear==2017) datatype = SingleElectron2017;
  if(RunYear==2018) datatype = SingleElectron2018;
  TFile *f; 

  if(RunYear==2016 || RunYear==2017){
  	f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/Tribeni/SingleElectron/SingleElectron_%d%s.root",RunYear,Era));
	apply_L1=true;}

  if(RunYear==2018){
  	f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/Tribeni/DoubleEG/DoubleEG_%d%s.root",RunYear,Era));
	apply_HEMveto=true; }
	
   cout<<"Applying L1 prefiring prob.? "<<apply_L1<<endl;
   cout<<"Applying HEM veto? "<<apply_HEMveto<<endl;

  TTree *es =(TTree*)f->Get("ggNtuplizer/EventTree");
  const unsigned nEvts = es->GetEntries();
  std::cout << "Total event: " << nEvts << std::endl;
  
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  float MET(0);
  float METPhi(0);
  int   ntrks(0);
  int   nvtx(0);
  int passHEM(0);
  TRandom3 ran(0);
  int count = 0;
    std::cout << "total: " << nEvts << std::endl;
    for (unsigned ievt(0); ievt<2000000; ++ievt){//loop on entries
  
        if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;

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

        MET = raw.pfMET;
        METPhi = raw.pfMETPhi;
        nvtx = raw.nVtx;

        
	if(RunYear==2018 && !passHEMVeto(0,raw)) continue; 
        passHEM++;

        if(MET > 70.0)continue;
        //if(!raw.passHLT())continue;
				//if(RunYear==2016 && ((raw.HLTEleMuX >> 4) &1) ==0)continue;  //HLT_Ele27_WPTight_Gsf_v
				//if(RunYear==2017 && ((raw.HLTEleMuX >> 3) &1) ==0)continue;  //HLT_Ele35_WPTight_Gsf_v
				//if(RunYear==2018 && ((raw.HLTEleMuX >> 55) &1) ==0)continue; //HLT_Ele32_WPTight_Gsf_v
				if(RunYear==2016 && ((raw.HLTEleMuX >> 2) &1) ==0)continue;  
				if(RunYear==2017 && ((raw.HLTEleMuX >> 59) &1) ==0)continue;  
				if(RunYear==2018 && ((raw.HLTEleMuX >> 59) &1) ==0)continue; 
        std::vector<std::vector<recoEle>::iterator> ElectronCollection;
        ElectronCollection.clear();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
		// use 30 for 2016 only, use 35 for 2018 and 38 for 2017
		   if(itEle->getCalibEt() < 30 || fabs(itEle->getEta())>2.1)continue;                              // Tag electron selection
		   //if(!itEle->passHLTSelection())continue;
						 if(RunYear==2016 && !itEle->fireTrgs(10))continue;   //HLT_Ele27_WPLoose_Gsf_v
						 if(RunYear==2017 && !itEle->fireTrgs(42))continue;  
						 if(RunYear==2018 && !itEle->fireTrgs(42))continue;  
		   if(itEle->passSignalSelection())ElectronCollection.push_back(itEle); // Tag electron selection
			// pt > 30 GeV, medium ID, eta < 2.1 electron
		}
		
           
	count++;
   }//loop on  events
cout<< "FakeRateRandomTree events: " << 1.*count/2000000. <<"; "<<std::endl;
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    Data_TriggerChecks(atoi(argv[1]),argv[2]);
    return 0;
}
