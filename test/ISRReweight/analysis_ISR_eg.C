// analysis_ISR_eg.C
// ==================
// Derives the Z+gamma control region tree for ISR reweighting
// in the EGAMMA channel using DoubleEG data.
//
// KEY DIFFERENCE from analysis_ISR.C (mg channel):
//   - Uses DoubleEG data stream (not MuonEG)
//   - Requires 2 signal electrons + 1 signal photon
//   - M(ee) in [80,100] GeV defines the Z+gamma CR
//   - Stores ISRJetPt AND gen-level bosonPt (for MC version)
//
// Rec. 1: Separate ISR CR per channel using correct data stream
// Rec. 2: Store bosonPt for gen-level reweighting
//
// Output: resTree_ISR_data_eg_{year}.root containing ZTree

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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_jet.h"

bool passMETFilter(int filter){
    bool passfilter(true);
    for(int im(1); im <= 8; im++)
        if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}
    if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
    if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}
    return passfilter;
}

void analysis_ISR_eg(){

    gSystem->Load("../../lib/libAnaClasses.so");

    // OUTPUT — change year tag as needed
    char outputname[200] = "/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_eg_2018.root";
    ofstream logfile;
    logfile.open("resTree_ISR_data_eg_2018.log");
    logfile << "analysis_ISR_eg(): Z+gamma CR for DoubleEG channel" << std::endl;

    // ---- INPUT: DoubleEG data (change to correct ggNtuple path per year) ----
    RunType datatype(DoubleEG2018);
    bool isMC = false;
    TChain* es = new TChain("ggNtuplizer/EventTree");
    es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/2018/DoubleEG_2018A.root");
    es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/2018/DoubleEG_2018B.root");
    es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/2018/DoubleEG_2018C.root");
    es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/2018/DoubleEG_2018D.root");
    // NOTE: adjust file paths for each year as in analysis_eg.C

    const unsigned nEvts = es->GetEntries();
    logfile << "Total events: " << nEvts << std::endl;

    TFile *outputfile = TFile::Open(outputname, "RECREATE");
    outputfile->cd();

    // ---- ZTree (same structure as analysis_ISR.C for compatibility) --------
    TTree *Ztree = new TTree("ZTree","ZTree");
    float MCweight = 1.0;  // =1 for data; set to lumi*xs/N for MC version
    float ZphoEt(0), ZphoEta(0), ZphoPhi(0);
    float ZlepPt(0),  ZlepEta(0),  ZlepPhi(0);   // leading electron
    float ZtrailPt(0),ZtrailEta(0),ZtrailPhi(0); // trailing electron
    float ZsigMT(0),  ZsigMET(0),  ZsigMETPhi(0);
    float ZdPhiLepMET(0), ZthreeMass(0), ZdilepMass(0);
    int   ZnVertex(0);
    float ZdRPhoLep(0), ZHT(0), ZnJet(0);
    float ZISRJetPt(0);    // Rec.2: reco vector-sum jet pT (kept for backward compat)
    float ZbosonPt(0);     // Rec.2: NEW — gen boson pT (filled in MC version; 0 for data)

    Ztree->Branch("MCweight",   &MCweight);
    Ztree->Branch("phoEt",      &ZphoEt);
    Ztree->Branch("phoEta",     &ZphoEta);
    Ztree->Branch("phoPhi",     &ZphoPhi);
    Ztree->Branch("lepPt",      &ZlepPt);
    Ztree->Branch("lepEta",     &ZlepEta);
    Ztree->Branch("lepPhi",     &ZlepPhi);
    Ztree->Branch("trailPt",    &ZtrailPt);
    Ztree->Branch("trailEta",   &ZtrailEta);
    Ztree->Branch("trailPhi",   &ZtrailPhi);
    Ztree->Branch("sigMT",      &ZsigMT);
    Ztree->Branch("sigMET",     &ZsigMET);
    Ztree->Branch("sigMETPhi",  &ZsigMETPhi);
    Ztree->Branch("dPhiLepMET", &ZdPhiLepMET);
    Ztree->Branch("threeMass",  &ZthreeMass);
    Ztree->Branch("dilepMass",  &ZdilepMass);
    Ztree->Branch("nVertex",    &ZnVertex);
    Ztree->Branch("dRPhoLep",   &ZdRPhoLep);
    Ztree->Branch("HT",         &ZHT);
    Ztree->Branch("nJet",       &ZnJet);
    Ztree->Branch("ISRJetPt",   &ZISRJetPt);
    Ztree->Branch("bosonPt",    &ZbosonPt);   // NEW — gen boson pT

    TH1F *p_eventcount = new TH1F("p_eventcount","",7,0,7);

    rawData raw(es, datatype);
    std::vector<recoPhoton> Photon;
    std::vector<recoMuon>   Muon;
    std::vector<recoEle>    Ele;
    std::vector<recoJet>    JetCollection;
    std::vector<mcData>     MCData;
    float MET(0), METPhi(0);
    int   nVtx(0), METFilter(0);
    int   nTotal(0), npassHLT(0), npassPho(0), npassLep(0), npassZ(0);

    for(unsigned ievt(0); ievt < nEvts; ++ievt){

        if(ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;

        raw.GetData(es, ievt);
        Photon.clear(); Muon.clear(); Ele.clear(); JetCollection.clear(); MCData.clear();
        for(int i=0; i<raw.nPho; i++) Photon.push_back(recoPhoton(raw,i));
        for(int i=0; i<raw.nMu;  i++) Muon.push_back(recoMuon(raw,i));
        for(int i=0; i<raw.nEle; i++) Ele.push_back(recoEle(raw,i));
        for(int i=0; i<raw.nJet; i++) JetCollection.push_back(recoJet(raw,i));
        MET      = raw.pfMET;
        METPhi   = raw.pfMETPhi;
        METFilter= raw.metFilters;
        nVtx     = raw.nVtx;

        nTotal++;
        if(!raw.passHLT())       continue;
        if(raw.nGoodVtx < 1)     continue;
        npassHLT++;
        if(raw.nEle < 2 || raw.nPho < 1) continue;

        // ---- Select signal photon (same as analysis_eg.C) ------------------
        bool hasPho(false);
        std::vector<recoPhoton>::iterator signalPho = Photon.begin();
        for(auto itpho = Photon.begin(); itpho != Photon.end(); ++itpho){
            if(itpho->getR9() < 0.5)                    continue;
            if(!itpho->passHLTSelection())               continue;
            if(!itpho->passBasicSelection())             continue;
            bool PixelVeto = (itpho->PixelSeed() == 0);
            bool GSFveto(true), FSRVeto(true);
            for(auto ie = Ele.begin(); ie != Ele.end(); ie++){
                if(DeltaR(itpho->getEta(),itpho->getPhi(),ie->getEta(),ie->getPhi()) <= 0.02) GSFveto=false;
                if(DeltaR(itpho->getEta(),itpho->getPhi(),ie->getEta(),ie->getPhi()) <  0.3)  FSRVeto=false;
            }
            for(auto im = Muon.begin(); im != Muon.end(); im++)
                if(DeltaR(itpho->getEta(),itpho->getPhi(),im->getEta(),im->getPhi()) < 0.3 && im->getEt()>2.0) FSRVeto=false;
            if(!itpho->passSignalSelection()) continue;
            if(GSFveto && PixelVeto && FSRVeto){
                if(!hasPho){ hasPho=true; signalPho=itpho; npassPho++; }
            }
        }

        // ---- Select two signal electrons forming a Z ----------------------
        bool hasLep(false), foundZ(false);
        std::vector<recoEle>::iterator leadEle  = Ele.begin();
        std::vector<recoEle>::iterator trailEle = Ele.begin();
        for(auto it1 = Ele.begin(); it1 != Ele.end(); ++it1){
            if(it1->getCalibPt() < 25)              continue;
            if(!it1->passHLTSelection())             continue;
            if(!it1->passSignalSelection())          continue;
            if((it1->isEB() && it1->getR9()<0.5) || (it1->isEE() && it1->getR9()<0.8)) continue;
            for(auto it2 = it1+1; it2 != Ele.end(); ++it2){
                if(it2->getCalibPt() < 15)           continue;
                if(!it2->isMedium())                 continue;
                double mee = (it1->getCalibP4() + it2->getCalibP4()).M();
                if(mee > 80 && mee < 100){
                    // pick highest-pT pair
                    if(!foundZ || it1->getCalibPt() > leadEle->getCalibPt()){
                        foundZ   = true;
                        hasLep   = true;
                        leadEle  = it1;
                        trailEle = it2;
                    }
                }
            }
        }
        if(!hasLep) continue;
        npassLep++;

        // ---- Require signal photon + Z + dR + MET filter ------------------
        if(!hasPho || !foundZ) continue;
        double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(),
                                    leadEle->getEta(),   leadEle->getPhi());
        if(dRlepphoton < 0.8) continue;
        if(!passMETFilter(METFilter)) continue;
        npassZ++;

        // ---- Build jet quantities ------------------------------------------
        ZnJet    = 0; ZHT = 0;
        TLorentzVector JetVec(0,0,0,0);
        for(auto itJet = JetCollection.begin(); itJet != JetCollection.end(); ++itJet){
            if(!itJet->passSignalSelection()) continue;
            if(DeltaR(itJet->getEta(),itJet->getPhi(),leadEle->getEta(), leadEle->getPhi())  <= 0.4) continue;
            if(DeltaR(itJet->getEta(),itJet->getPhi(),trailEle->getEta(),trailEle->getPhi()) <= 0.4) continue;
            if(DeltaR(itJet->getEta(),itJet->getPhi(),signalPho->getEta(),signalPho->getPhi()) <= 0.4) continue;
            ZnJet++;
            ZHT += itJet->getPt();
            JetVec = JetVec + itJet->getP4();
        }
        ZISRJetPt = JetVec.Pt();
        ZbosonPt  = 0; // data: no gen info; MC version fills from mcData loop

        // ---- Fill ZTree ----------------------------------------------------
        float deltaPhi = DeltaPhi(leadEle->getPhi(), METPhi);
        double mee     = (leadEle->getCalibP4() + trailEle->getCalibP4()).M();

        ZphoEt      = signalPho->getCalibEt();
        ZphoEta     = signalPho->getEta();
        ZphoPhi     = signalPho->getPhi();
        ZlepPt      = leadEle->getCalibPt();
        ZlepEta     = leadEle->getEta();
        ZlepPhi     = leadEle->getPhi();
        ZtrailPt    = trailEle->getCalibPt();
        ZtrailEta   = trailEle->getEta();
        ZtrailPhi   = trailEle->getPhi();
        ZsigMT      = sqrt(2*MET*leadEle->getCalibPt()*(1-std::cos(deltaPhi)));
        ZsigMET     = MET;
        ZsigMETPhi  = METPhi;
        ZdPhiLepMET = deltaPhi;
        ZdilepMass  = mee;
        ZthreeMass  = (leadEle->getCalibP4()+trailEle->getCalibP4()+signalPho->getCalibP4()).M();
        ZnVertex    = nVtx;
        ZdRPhoLep   = dRlepphoton;

        Ztree->Fill();

    }// event loop

    p_eventcount->Fill("Total",      nTotal);
    p_eventcount->Fill("passHLT",    npassHLT);
    p_eventcount->Fill("passPho",    npassPho);
    p_eventcount->Fill("passLep",    npassLep);
    p_eventcount->Fill("passZ",      npassZ);

    outputfile->Write();
    logfile << "passZ events: " << npassZ << std::endl;
    logfile.close();
    std::cout << "Done. passZ = " << npassZ << std::endl;
}
