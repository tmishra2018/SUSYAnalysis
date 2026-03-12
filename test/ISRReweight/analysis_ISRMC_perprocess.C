// analysis_ISRMC_perprocess.C
// ============================
// Produces ISR reweighting MC trees for each relevant process
// separately: ZGToLLG, WGToLNuG, TTGJets, WWG, WZG
//
// KEY DIFFERENCES from analysis_ISRMC.C:
//   Rec.2: Stores gen-level bosonPt (Z pT or W pT from MC truth)
//          In addition to reco ISRJetPt, for gen-level reweighting
//   Rec.3: Separate output file per process — do NOT mix Wgamma
//          and Zgamma into the same reweight derivation
//   Rec.4: Rare backgrounds (WWG, WZG, TTG) get their own trees
//          so the Zgamma-derived weight is NOT applied to them
//
// Usage: change PROC_TAG and mcType before compiling, then run
//   root -b -q "analysis_ISRMC_perprocess.C+"
//
// Supported process tags:
//   "ZGToLLG"   mcType = MCType::ZGInclusive
//   "WGToLNuG"  mcType = MCType::WGJetInclusive   (or WGJet40/130)
//   "TTGJets"   mcType = MCType::TTG
//   "WWG"       mcType = MCType::WWG
//   "WZG"       mcType = MCType::WZG

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_jet.h"

// ============================================================
// USER CONFIG — change these per run
// ============================================================
const std::string PROC_TAG = "ZGToLLG";   // label for output file
const int         mcType   = MCType::ZGInclusive;
const int         RunYear  = 2018;
const bool        isWchannel = false; // true for Wgamma processes (1 mu + MET topology)
                                      // false for Zgamma (2 mu topology)
// Input file — match to your ggNtuple path
const char* INPUT_FILE = "root://cmseos.fnal.gov//store/user/msun/MCUL18/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.root";
// ============================================================

void analysis_ISRMC_perprocess(){

    gSystem->Load("../../lib/libAnaClasses.so");

    char outputname[300];
    sprintf(outputname,
        "/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_%s_%d.root",
        PROC_TAG.c_str(), RunYear);
    ofstream logfile;
    char logname[300];
    sprintf(logname, "resTree_ISR_%s_%d.log", PROC_TAG.c_str(), RunYear);
    logfile.open(logname);
    logfile << "analysis_ISRMC_perprocess(): " << PROC_TAG << " " << RunYear << std::endl;

    RunType datatype(MCMuonEG);
    TChain* es = new TChain("ggNtuplizer/EventTree");
    es->Add(INPUT_FILE);

    const unsigned nEvts = es->GetEntries();
    float PUweight(1);
    float crosssection = MC_XS[mcType];
    float MCweight     = crosssection * 59.74 * 1000.0 / nEvts; // adjust lumi per year
    logfile << "xs=" << crosssection << "  MCweight=" << MCweight
            << "  nEvts=" << nEvts << std::endl;

    TFile *outputfile = TFile::Open(outputname,"RECREATE");
    outputfile->cd();

    // ---- ZTree (extended: bosonPt + procType for per-process weighting) ----
    TTree *Ztree = new TTree("ZTree","ZTree");
    float ZphoEt(0),  ZphoEta(0),   ZphoPhi(0);
    float ZlepPt(0),  ZlepEta(0),   ZlepPhi(0);
    float ZtrailPt(0),ZtrailEta(0), ZtrailPhi(0);
    float ZsigMT(0),  ZsigMET(0),   ZsigMETPhi(0);
    float ZdPhiLepMET(0), ZthreeMass(0), ZdilepMass(0);
    int   ZnVertex(0);
    float ZdRPhoLep(0), ZHT(0), ZnJet(0);
    float ZISRJetPt(0);   // reco vector-sum jet pT (backward compat)
    float ZbosonPt(0);    // Rec.2: gen-level boson pT (Z or W pT from MC truth)
    float ZbosonPt_gen(0);// Rec.2: same as above, explicit gen label
    int   ZprocType(0);   // Rec.3/4: mcType integer for per-process selection
    std::vector<int>   Z_mcPID;
    std::vector<float> Z_mcEta, Z_mcPhi, Z_mcPt;
    std::vector<int>   Z_mcMomPID, Z_mcGMomPID, Z_mcStatus;

    Ztree->Branch("MCweight",   &MCweight);
    Ztree->Branch("PUweight",   &PUweight);
    Ztree->Branch("mcType",     &ZprocType);
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
    Ztree->Branch("ISRJetPt",   &ZISRJetPt);     // reco (for backward compat)
    Ztree->Branch("bosonPt",    &ZbosonPt);      // gen boson pT
    Ztree->Branch("bosonPt_gen",&ZbosonPt_gen);  // same, explicit
    Ztree->Branch("mcPID",      &Z_mcPID);
    Ztree->Branch("mcEta",      &Z_mcEta);
    Ztree->Branch("mcPhi",      &Z_mcPhi);
    Ztree->Branch("mcPt",       &Z_mcPt);
    Ztree->Branch("mcMomPID",   &Z_mcMomPID);
    Ztree->Branch("mcGMomPID",  &Z_mcGMomPID);
    Ztree->Branch("mcStatus",   &Z_mcStatus);

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
        MCData.clear(); Photon.clear(); Muon.clear(); Ele.clear(); JetCollection.clear();
        for(int i=0;i<raw.nMC; i++) MCData.push_back(mcData(raw,i));
        for(int i=0;i<raw.nPho;i++) Photon.push_back(recoPhoton(raw,i));
        for(int i=0;i<raw.nMu; i++) Muon.push_back(recoMuon(raw,i));
        for(int i=0;i<raw.nEle;i++) Ele.push_back(recoEle(raw,i));
        for(int i=0;i<raw.nJet;i++) JetCollection.push_back(recoJet(raw,i));
        MET      = raw.pfMET;
        METPhi   = raw.pfMETPhi;
        METFilter= raw.metFilters;
        nVtx     = raw.nVtx;
        PUweight = getPUESF(nVtx);
        ZprocType= mcType;

        nTotal++;
        if(!raw.passHLT())    continue;
        if(raw.nGoodVtx < 1)  continue;
        npassHLT++;

        // ---- Rec.2: Extract gen-level boson pT from MC truth ---------------
        // IMPORTANT: mcData wrapper only stores status=1 final-state particles,
        // so Z/W bosons (status=22 or status=62 in Pythia8) are NOT in MCData.
        // Read raw ggNtuple arrays (raw.mcPID, raw.mcStatus, raw.mcPt) directly.
        //
        // Strategy: find the LAST occurrence of PID=23 (Z) or PID=24 (W) with
        // status=22 (intermediate resonance, Pythia8 convention).
        // "Last" avoids the incoming parton-shower copies and picks the
        // on-shell boson just before it decays.
        // Fallback: if no status=22 copy found, take any status=62 or status=52.
        ZbosonPt = 0; ZbosonPt_gen = 0;
        int targetPID = isWchannel ? 24 : 23;
        float bestPt_s22 = -1, bestPt_s62 = -1, bestPt_any = -1;
        for(int imc=0; imc < raw.nMC; imc++){
            if(abs(raw.mcPID[imc]) != targetPID) continue;
            int   st = raw.mcStatus[imc];
            float pt = raw.mcPt[imc];
            if(st == 22)                 bestPt_s22 = pt;   // keep last s22
            if(st == 62 || st == 52)     bestPt_s62 = pt;
            bestPt_any = pt;                                 // keep last match
        }
        float bosonPtFound = (bestPt_s22 > 0) ? bestPt_s22
                           : (bestPt_s62 > 0) ? bestPt_s62
                           : bestPt_any;
        if(bosonPtFound > 0){
            ZbosonPt     = bosonPtFound;
            ZbosonPt_gen = bosonPtFound;
        }

        // ---- Require correct lepton multiplicity per channel ---------------
        int nGoodMu(0);
        for(auto im = Muon.begin(); im != Muon.end(); im++)
            if(im->isMedium() && im->getMiniIso()<0.2 && im->getPt()>15) nGoodMu++;

        // Zgamma topology: 2 muons
        if(!isWchannel && (nGoodMu != 2 || raw.nPho < 1)) continue;
        // Wgamma topology: exactly 1 tight muon (see plot_ISRweight_Wchannel.C)
        if(isWchannel  && (raw.nMu  < 1  || raw.nPho < 1)) continue;

        // ---- Select signal photon ------------------------------------------
        bool hasPho(false);
        std::vector<recoPhoton>::iterator signalPho = Photon.begin();
        for(auto itpho = Photon.begin(); itpho != Photon.end(); ++itpho){
            if(itpho->getR9() < 0.5)          continue;
            if(!itpho->passHLTSelection())     continue;
            if(!itpho->passBasicSelection())   continue;
            bool PixelVeto = (itpho->PixelSeed()==0);
            bool GSFveto(true), FSRVeto(true);
            for(auto ie = Ele.begin(); ie != Ele.end(); ie++){
                if(DeltaR(itpho->getEta(),itpho->getPhi(),ie->getEta(),ie->getPhi())<=0.02) GSFveto=false;
                if(DeltaR(itpho->getEta(),itpho->getPhi(),ie->getEta(),ie->getPhi())< 0.3)  FSRVeto=false;
            }
            for(auto im = Muon.begin(); im != Muon.end(); im++)
                if(DeltaR(itpho->getEta(),itpho->getPhi(),im->getEta(),im->getPhi())<0.3 && im->getEt()>2.0) FSRVeto=false;
            if(!itpho->passSignalSelection()) continue;
            if(GSFveto && PixelVeto && FSRVeto && !hasPho){
                hasPho=true; signalPho=itpho; npassPho++;
            }
        }

        // ---- Select leptons and find Z/W candidate ------------------------
        bool hasLep(false), foundZG(false);
        std::vector<recoMuon>::iterator signalLep  = Muon.begin();
        std::vector<recoMuon>::iterator trailLepIt = Muon.begin();
        for(auto itMu = Muon.begin(); itMu != Muon.end(); ++itMu){
            if(itMu->getPt() < 25)          continue;
            if(!itMu->passHLTSelection())   continue;
            if(itMu->passSignalSelection()){
                if(!hasLep){ hasLep=true; signalLep=itMu; npassLep++; }
            }
        }
        if(!hasLep || !hasPho) continue;

        // For Z topology: find second muon with M(mumu) in [80,100]
        if(!isWchannel){
            for(auto im = Muon.begin(); im != Muon.end(); ++im){
                if(im == signalLep) continue;
                if(!im->isMedium() || im->getMiniIso()>0.2) continue;
                double llmass = (im->getP4()+signalLep->getP4()).M();
                if(llmass>80 && llmass<100){ foundZG=true; trailLepIt=im; break; }
            }
            if(!foundZG) continue;
        }
        else{
            // W topology: single lepton + MET; no second muon required
            // Use a loose MT cut to select W-like events
            float deltaPhi_W = DeltaPhi(signalLep->getPhi(), METPhi);
            float MT_W = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi_W)));
            if(MT_W < 30 || MT_W > 120) continue; // W transverse mass window
            foundZG = true;
            trailLepIt = signalLep; // dummy: same as lead for W
        }

        double dRlepphoton = DeltaR(signalPho->getEta(),signalPho->getPhi(),
                                     signalLep->getEta(),signalLep->getPhi());
        if(dRlepphoton < 0.8) continue;
        if(!raw.passMETFilter(METFilter)) continue;
        npassZ++;

        // ---- Jets ----------------------------------------------------------
        ZnJet=0; ZHT=0;
        TLorentzVector JetVec(0,0,0,0);
        for(auto itJet = JetCollection.begin(); itJet != JetCollection.end(); ++itJet){
            if(!itJet->passSignalSelection()) continue;
            if(DeltaR(itJet->getEta(),itJet->getPhi(),signalLep->getEta(),signalLep->getPhi()) <=0.4) continue;
            if(!isWchannel && DeltaR(itJet->getEta(),itJet->getPhi(),trailLepIt->getEta(),trailLepIt->getPhi())<=0.4) continue;
            if(DeltaR(itJet->getEta(),itJet->getPhi(),signalPho->getEta(),signalPho->getPhi())<=0.4) continue;
            ZnJet++; ZHT += itJet->getPt();
            JetVec = JetVec + itJet->getP4();
        }
        ZISRJetPt = JetVec.Pt();

        // ---- MC matching info ----------------------------------------------
        Z_mcPID.clear(); Z_mcEta.clear(); Z_mcPhi.clear(); Z_mcPt.clear();
        Z_mcMomPID.clear(); Z_mcGMomPID.clear(); Z_mcStatus.clear();
        for(auto itMC = MCData.begin(); itMC != MCData.end(); ++itMC){
            if(itMC->getEt() < 1.0) continue;
            float dR_pho = DeltaR(signalPho->getEta(),signalPho->getPhi(),itMC->getEta(),itMC->getPhi());
            float dR_lep = DeltaR(signalLep->getEta(),signalLep->getPhi(),itMC->getEta(),itMC->getPhi());
            if(dR_pho < 0.3 || dR_lep < 0.3){
                Z_mcPID.push_back(itMC->getPID());
                Z_mcMomPID.push_back(itMC->getMomPID());
                Z_mcGMomPID.push_back(itMC->getGMomPID());
                Z_mcEta.push_back(itMC->getEta());
                Z_mcPhi.push_back(itMC->getPhi());
                Z_mcPt.push_back(itMC->getEt());
                Z_mcStatus.push_back(itMC->getStatus());
            }
        }

        // ---- Fill ----------------------------------------------------------
        float dPhi = DeltaPhi(signalLep->getPhi(), METPhi);
        double mll = isWchannel ? 0.0 : (trailLepIt->getP4()+signalLep->getP4()).M();

        ZphoEt      = signalPho->getCalibEt();
        ZphoEta     = signalPho->getEta();
        ZphoPhi     = signalPho->getPhi();
        ZlepPt      = signalLep->getPt();
        ZlepEta     = signalLep->getEta();
        ZlepPhi     = signalLep->getPhi();
        ZtrailPt    = isWchannel ? 0 : trailLepIt->getPt();
        ZtrailEta   = isWchannel ? 0 : trailLepIt->getEta();
        ZtrailPhi   = isWchannel ? 0 : trailLepIt->getPhi();
        ZsigMT      = sqrt(2*MET*signalLep->getPt()*(1-std::cos(dPhi)));
        ZsigMET     = MET;
        ZsigMETPhi  = METPhi;
        ZdPhiLepMET = dPhi;
        ZdilepMass  = mll;
        ZthreeMass  = isWchannel ? 0 : (trailLepIt->getP4()+signalLep->getP4()+signalPho->getCalibP4()).M();
        ZnVertex    = nVtx;
        ZdRPhoLep   = dRlepphoton;

        Ztree->Fill();

    }// event loop

    p_eventcount->Fill("Total",   nTotal);
    p_eventcount->Fill("passHLT", npassHLT);
    p_eventcount->Fill("passPho", npassPho);
    p_eventcount->Fill("passLep", npassLep);
    p_eventcount->Fill("passZ",   npassZ);

    outputfile->Write();
    logfile << "passZ/W events: " << npassZ << std::endl;
    logfile.close();
    std::cout << "Done. passZ/W = " << npassZ << std::endl;
}
