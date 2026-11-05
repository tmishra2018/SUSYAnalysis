// Run using g++ `root-config --cflags` analysis_elefakepho.C -o analysis_elefakepho.exe `root-config --libs`
//
// Merged DY (MC) / Data version of the electron-fakes-photon tag-and-probe analysis.
// Usage: ./analysis_elefakepho.exe <RunYear> <Era> <Mode>
//   RunYear : 2016, 2017, 2018
//   Era     : e.g. preVFP, postVFP, "" ...
//   Mode    : "DY"   -> runs on DYJetsToLL MC
//             "Data" -> runs on SingleElectron (2016/2017) or DoubleEG (2018) data
//
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
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_cuts.h"
#include "../../../src/analysis_ele.cc"
#include "../../../src/analysis_photon.cc"

bool apply_HEMveto=false;
bool apply_L1=false;

void analysis_elefakepho(int RunYear, const char *Era, const char *Mode){//main

  TString mode(Mode);
  bool isData = mode.EqualTo("Data", TString::kIgnoreCase);
  bool isDY   = mode.EqualTo("DY",   TString::kIgnoreCase);
  if(!isData && !isDY){
    std::cout << "Unknown Mode \"" << Mode << "\": must be \"DY\" or \"Data\"" << std::endl;
    throw;
  }
  TString sampleTag = isData ? "Data" : "DY";

  ofstream logfile;
  logfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/logs/plot_elefakepho_%sTnP_dR05_%s_%d%s.log",
                     sampleTag.Data(), sampleTag.Data(), RunYear, Era), ios::trunc);

  logfile << "analysis_elefakepho()" << std::endl;

  RunType datatype;
  TFile *f;

  if(isDY){
    datatype = MC;                                     // Run Type
    f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DYJetsToLL/DYJetsToLL_%d%s.root",RunYear,Era));
  } else { // isData
    if(RunYear==2016) datatype = SingleElectron2016;
    if(RunYear==2017) datatype = SingleElectron2017;
    if(RunYear==2018) datatype = SingleElectron2018;

    if(RunYear==2016 || RunYear==2017){
      f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/SingleElectron/SingleElectron_%d%s.root",RunYear,Era));
      apply_L1=true;
    }
    if(RunYear==2018){
      f = TFile::Open(Form("/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_%d%s.root",RunYear,Era));
      apply_HEMveto=true;
    }

    cout<<"Applying L1 prefiring prob.? "<<apply_L1<<endl;
    cout<<"Applying HEM veto? "<<apply_HEMveto<<endl;
  }

  TTree *es =(TTree*)f->Get("ggNtuplizer/EventTree");

  TFile *output = TFile::Open(Form("/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_%sTnP_dR05_%d%s.root",
                                    sampleTag.Data(), RunYear, Era), "RECREATE");
  output->cd();

  int   tracks(0);
  int   nVertex(0);
  float PUweight(1);  // true-pileup weight (DY), 1 for data
  int   mcType = isDY ? MCType::DYLL50 : MCType::NOMC;
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
  etree->Branch("PUweight",            &PUweight);
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
  etree->Branch("FSRveto",             &FSRveto_bothcount);
  etree->Branch("mcPID",               &mcPID_bothcount);
  etree->Branch("mcEta",               &mcEta_bothcount);
  etree->Branch("mcPhi",               &mcPhi_bothcount);
  etree->Branch("mcPt",                &mcPt_bothcount);
  etree->Branch("mcMomPID",            &mcMomPID_bothcount);
  etree->Branch("mcGMomPID",           &mcGMomPID_bothcount);
  etree->Branch("mcType",              &mcType);
  float HT(0);
  etree->Branch("HT",        &HT);

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
  rantree->Branch("PUweight",            &PUweight);
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
  rantree->Branch("FSRveto",             &FSRveto_random);
  rantree->Branch("mcPID",               &mcPID_random);
  rantree->Branch("mcEta",               &mcEta_random);
  rantree->Branch("mcPhi",               &mcPhi_random);
  rantree->Branch("mcPt",                &mcPt_random);
  rantree->Branch("mcMomPID",            &mcMomPID_random);
  rantree->Branch("mcGMomPID",           &mcGMomPID_random);
  rantree->Branch("mcType",              &mcType);

  const unsigned nEvts = es->GetEntries();
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << "/eos/uscms/store/user/tmishra/elefakepho/files/plot_elefakepho_"
          << sampleTag << "TnP_dR05_" << RunYear << Era << ".root" << std::endl;

  // True-pileup reweighting for DY (as VGamma/src/analysis_VGamma_Run2.C): data nTrueInt profile
  // (golden JSON, 69.2 mb) over the nTrueInt (BX==0) profile of this MC sample
  TH1D *h_puRatio = 0;
  if(isDY){
    TH1D *h_mc_pu = new TH1D("h_mc_pu","MC nTrueInt",200,0,200);
    h_mc_pu->SetDirectory(0);
    es->Project("h_mc_pu", "puTrue", "puBX==0");
    TFile *f_pu = TFile::Open(Form("../../PUfiles/dataPU_%d%s.root",RunYear,Era));
    if(!f_pu || f_pu->IsZombie()){std::cout << "missing ../../PUfiles/dataPU_" << RunYear << Era << ".root" << std::endl; throw;}
    h_puRatio = (TH1D*)((TH1D*)f_pu->Get("pileup"))->Clone("h_puRatio");
    h_puRatio->SetDirectory(0);
    f_pu->Close();
    output->cd();
    h_puRatio->Scale(1.0/h_puRatio->Integral());
    h_mc_pu->Scale(1.0/h_mc_pu->Integral());
    h_puRatio->Divide(h_mc_pu);
    logfile << "true-PU weights from ../../PUfiles/dataPU_" << RunYear << Era << ".root, MC nTrueInt entries " << h_mc_pu->GetEntries() << std::endl;
  }

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
  int passHEM(0), passPixelIssue(0);
  TRandom3 ran(0);
  int pTcut = 0;
  if (RunYear==2016) pTcut = 30;
  if (RunYear==2017) pTcut = 38;
  if (RunYear==2018) pTcut = 35;

  std::cout << "total: " << nEvts << std::endl;
  const bool applyGolden = (isData);  // golden JSON lumi mask for data
  if(applyGolden && !loadGolden(RunYear)) return;
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
        if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
        raw.GetData(es, ievt);
        if(applyGolden && !passGolden(raw.run, raw.lumis)) continue;
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

        tracks = ntrks;
        nVertex = nvtx;
        PUweight = 1;
        if(h_puRatio && raw.puBX && raw.puTrue){
          for(size_t ibx = 0; ibx < raw.puBX->size(); ibx++){
            if((*raw.puBX)[ibx] != 0) continue;
            int b = std::max(1, std::min(h_puRatio->FindBin((*raw.puTrue)[ibx]), h_puRatio->GetNbinsX()));
            if(h_puRatio->GetBinContent(b) > 0) PUweight = h_puRatio->GetBinContent(b);
            break;
          }
        }

        // Detector cleaning cuts, applied to both Data and DY.
        if(RunYear==2018 && isData && !passHEMVeto(0, raw, 0.0)) continue;
        if(RunYear==2018 && isDY && !passHEMVeto(0, raw, fracHEM_eg)) continue;
        passHEM++;

        if(RunYear==2017 && !passPixelIssue17(raw)) continue;
        if(RunYear==2018 && !passPixelIssue18(raw)) continue;
        passPixelIssue++;

        if(MET > 70.0)continue;
        //if(!raw.passHLT())continue; Explicitly for each year as mentioned below.
        if(RunYear==2016 && ((raw.HLTEleMuX >> 4) &1) ==0)continue;  //HLT_Ele27_WPTight_Gsf_v
        if(RunYear==2017 && ((raw.HLTEleMuX >> 3) &1) ==0)continue;  //HLT_Ele35_WPTight_Gsf_v
        if(RunYear==2018 && ((raw.HLTEleMuX >> 55) &1) ==0)continue; //HLT_Ele32_WPTight_Gsf_v

        std::vector<std::vector<recoEle>::iterator> ElectronCollection;
        ElectronCollection.clear();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
           if(itEle->getCalibEt() < pTcut || fabs(itEle->getEta())>2.1)continue;                              // Tag electron selection
           //if(!itEle->passHLTSelection())continue;  Explicitly for each year as mentioned below.
           if(RunYear==2016 && !itEle->fireTrgs(12))continue;  //HLT_Ele27_WPTight_Gsf_v
           if(RunYear==2017 && !itEle->fireTrgs(46))continue;  //HLT_Ele35_WPTight_Gsf_v
           if(RunYear==2018 && !itEle->fireTrgs(13))continue;  //HLT_Ele32_WPTight_Gsf_v
           if(itEle->passSignalSelection())ElectronCollection.push_back(itEle); // Tag electron selection
           // pt > 30 GeV, medium ID, eta < 2.1 electron
        }

        bool hasEle(false);
        std::vector<recoEle>::iterator signalEle = Ele.begin();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
                   if(itEle->getCalibPt() < 25)continue;
                   if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;
                   if(itEle->passSignalSelection()){
                         if(!hasEle){
                                hasEle=true;
                                signalEle = itEle;
                          }
                   }
        }

        std::vector<recoPhoton>::iterator LeadingsignalPho = Photon.begin();
        bool hasLeadingSigPho(false);
        for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
                                if(itpho->getR9() < 0.5)continue;
                                if(!itpho->passBasicSelection())continue;
                                bool passSigma = itpho->passSigma(1);
                                bool passChIso = itpho->passChIso(1);
                                bool PixelVeto = itpho->PixelSeed()==0? true: false;
                                bool GSFveto(true);
                                bool FSRVeto(true);
                                for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
                                        if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
                                        if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
                                }
                                for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
                                        if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;
                                if(!itpho->passSignalSelection())continue;
                                if(GSFveto && PixelVeto && FSRVeto){
                                        if(!hasLeadingSigPho){
                                                hasLeadingSigPho=true;
                                                LeadingsignalPho = itpho;
                                        }
                                }
         }

         HT = 0;
         for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
                     if(!itJet->passSignalSelection())continue;
                     // Ele/Photon empty -> begin()==end(), dereferencing it crashed (2018A); other events unchanged
                     if(!Ele.empty() && DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(), signalEle->getPhi()) <= 0.4)continue;
                     if(!Photon.empty() && DeltaR(itJet->getEta(), itJet->getPhi(), LeadingsignalPho->getEta(), LeadingsignalPho->getPhi()) <= 0.4)continue;
                     HT += itJet->getPt();
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
if(RunYear==2018) logfile << "pass HEM cut:  " << passHEM*100/nEvts<<endl;
if(RunYear==2018) logfile << "pass Pixel veto and HEM veto " << passPixelIssue*100/nEvts<<endl;
if(RunYear==2017) logfile << "pass Pixel veto " << passPixelIssue*100/nEvts<<endl;

cout<< "FakeRateRandomTree events: " << rantree->GetEntries() <<"; "<<percent1<<"\%"<<std::endl;
cout<< "FakeRateTree events: " << etree->GetEntries() <<"; "<<percent2<<"\%"<<std::endl;
if(RunYear==2018) cout<< "pass HEM cut:  " << passHEM<<endl;

output->Write();
output->Close();
logfile.close();
}
int main(int argc, char** argv)
{
    if(argc < 4){
      cout << "You have to provide three arguments: <RunYear> <Era> <Mode: DY|Data>\n";
      return 1;
    }
    analysis_elefakepho(atoi(argv[1]),argv[2],argv[3]);
    return 0;
}
