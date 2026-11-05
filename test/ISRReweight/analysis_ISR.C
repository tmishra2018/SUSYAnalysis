// analysis_ISR.C
// Universal data skimmer for ISR JetPt weight derivation.
// Runs both channels in a single job, writing separate output files:
//   mg     channel: MuonEG dataset,       Z -> mumu + photon
//   egamma channel: DoubleEG/EGamma data, Z -> ee   + photon
//
// Build:
//   g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISR.C \
//       -o analysis_ISR.exe `root-config --libs`
// Usage:
//   ./analysis_ISR.exe <RunYear> <Era> [channel both|mg|egamma]
//   e.g.  ./analysis_ISR.exe 2018 D
//         ./analysis_ISR.exe 2016 B-ver1
//         ./analysis_ISR.exe 2017 B egamma   (MuonEG lumi covers 2017 D-F only)

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <cmath>

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
#include "TMemFile.h"

#include "../../src/analysis_rawData.cc"
 #include "../../include/analysis_rawData.h"
 #include "../../include/analysis_photon.h"
 #include "../../include/analysis_muon.h"
 #include "../../include/analysis_jet.h"
 #include "../../include/analysis_ele.h"
 #include "../../include/analysis_mcData.h"
 #include "../../include/analysis_tools.h"

#include "../../src/analysis_ele.cc"
#include "../../src/analysis_muon.cc"
#include "../../src/analysis_photon.cc"


// ---------------------------------------------------------------------------
// Shared ZTree branch struct
// ---------------------------------------------------------------------------
struct ZBranches {
    float phoEt, phoEta, phoPhi;
    float lepPt, lepEta, lepPhi;
    float trailPt, trailEta, trailPhi;
    float sigMT, sigMET, sigMETPhi, dPhiLepMET;
    float threeMass, dilepMass;
    int   nVertex;
    float dRPhoLep, HT, nJet, ISRJetPt;
    float MCweight;

    void book(TTree* t) {
        t->Branch("MCweight",   &MCweight);
        t->Branch("phoEt",      &phoEt);
        t->Branch("phoEta",     &phoEta);
        t->Branch("phoPhi",     &phoPhi);
        t->Branch("lepPt",      &lepPt);
        t->Branch("lepEta",     &lepEta);
        t->Branch("lepPhi",     &lepPhi);
        t->Branch("trailPt",    &trailPt);
        t->Branch("trailEta",   &trailEta);
        t->Branch("trailPhi",   &trailPhi);
        t->Branch("sigMT",      &sigMT);
        t->Branch("sigMET",     &sigMET);
        t->Branch("sigMETPhi",  &sigMETPhi);
        t->Branch("dPhiLepMET", &dPhiLepMET);
        t->Branch("threeMass",  &threeMass);
        t->Branch("dilepMass",  &dilepMass);
        t->Branch("nVertex",    &nVertex);
        t->Branch("dRPhoLep",   &dRPhoLep);
        t->Branch("HT",         &HT);
        t->Branch("nJet",       &nJet);
        t->Branch("ISRJetPt",   &ISRJetPt);
    }
};

// ---------------------------------------------------------------------------
// Photon selection helper: shared between mg and egamma.
// Fills sigPho (passes signal selection) and hadPho (hadron proxy: passes
// basic but not necessarily signal selection).
// ---------------------------------------------------------------------------
void selectPhoton(std::vector<recoPhoton>& Photon,
                  std::vector<recoEle>&    Ele,
                  std::vector<recoMuon>&   Muon,
                  std::vector<recoPhoton>::iterator& sigPho,
                  std::vector<recoPhoton>::iterator& hadPho,
                  bool& hasSig, bool& hasHad)
{
    hasSig = false; hasHad = false;
    for (auto it = Photon.begin(); it != Photon.end(); ++it) {
        if (it->getR9()             < 0.5) continue;
        if (!it->passHLTSelection())        continue;
        if (!it->passBasicSelection())      continue;

        bool pixVeto = (it->PixelSeed() == 0);
        bool gsfVeto(true), fsrVeto(true);

        for (auto& ele : Ele) {
            if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <= 0.02) gsfVeto = false;
            if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <  0.3)  fsrVeto = false;
        }
        for (auto& mu : Muon)
            if (DeltaR(it->getEta(), it->getPhi(), mu.getEta(), mu.getPhi()) < 0.3 && mu.getPt() > 2.0)
                fsrVeto = false;

        // hadron proxy candidate (loose: before signal selection)
        if (gsfVeto && pixVeto && fsrVeto && !hasHad) { hasHad = true; hadPho = it; }

        if (!it->passSignalSelection()) continue;
        if (gsfVeto && pixVeto && fsrVeto && !hasSig) { hasSig = true; sigPho = it; }
    }
}

// ---------------------------------------------------------------------------
// ISR jet pT: vector sum of jets not matched to either lepton or the photon
// (dR > 0.4), i.e. the recoil of the Z+gamma system (AN Sec. 3 jet cleaning;
// the photon is also clustered into an AK4 jet). Also fills HT and nJet.
// ---------------------------------------------------------------------------
float calcISRJetPt(std::vector<recoJet>& jets,
                   float l1Eta, float l1Phi,
                   float l2Eta, float l2Phi,
                   float phoEta, float phoPhi,
                   float& HT,   float& nJet)
{
    TLorentzVector jsys(0,0,0,0);
    HT = 0; nJet = 0;
    for (auto& jet : jets) {
        if (!jet.passSignalSelection()) continue;
        if (DeltaR(jet.getEta(), jet.getPhi(), l1Eta, l1Phi) <= 0.4) continue;
        if (DeltaR(jet.getEta(), jet.getPhi(), l2Eta, l2Phi) <= 0.4) continue;
        if (DeltaR(jet.getEta(), jet.getPhi(), phoEta, phoPhi) <= 0.4) continue;
        nJet++; HT += jet.getPt(); jsys += jet.getP4();
    }
    return jsys.Pt();
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void analysis_ISR(int RunYear, const char* Era, const char* channelArg = "both") {
    gSystem->Load("../../lib/libAnaClasses.so");

    const std::string channel = channelArg;
    const bool runMg = (channel == "both" || channel == "mg");
    const bool runEG = (channel == "both" || channel == "egamma");

    RunType rtMg, rtEG;
    if (RunYear == 2016) {
        rtMg = MuonEG2016;
        rtEG = DoubleEG2016;
    } else if (RunYear == 2017) {
        rtMg = MuonEG2017;
        rtEG = DoubleEG2017;
    } else {
        rtMg = MuonEG2018;
        rtEG = DoubleEG2018;
    }

    std::ofstream log(Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/logs/resTree_ISR_data_%s_%d%s.log",
        channelArg, RunYear, Era), std::ios::trunc);
    log << "analysis_ISR  year=" << RunYear << "  era=" << Era
        << "  rtMg=" << rtMg << "  rtEG=" << rtEG << "\n";

    TChain* esMg = new TChain("ggNtuplizer/EventTree");
    if (runMg) esMg->Add(Form(
        "/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_%d%s.root",
        RunYear, Era));

    TChain* esEG = new TChain("ggNtuplizer/EventTree");
    if (runEG) esEG->Add(Form(
            "/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_%d%s.root",
            RunYear, Era));

    log << "MuonEG   entries: " << esMg->GetEntries() << "\n";
    log << "EG       entries: " << esEG->GetEntries() << "\n";
    std::cout << "MuonEG   entries: " << esMg->GetEntries() << "\n";
    std::cout << "EG       entries: " << esEG->GetEntries() << "\n";

    // ---- Output files: one per channel (skipped channel -> in-memory only) --
    TFile* fMg = runMg ? TFile::Open(Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_mg_%d%s.root",
        RunYear, Era), "RECREATE") : new TMemFile("mem_mg.root", "RECREATE");
    TFile* fEG = runEG ? TFile::Open(Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_data_egamma_%d%s.root",
        RunYear, Era), "RECREATE") : new TMemFile("mem_eg.root", "RECREATE");

    fMg->cd();
    TTree* ZtreeMg = new TTree("ZTree", "ZTree (mg)");
    ZBranches brMg; brMg.MCweight = 1.0f; brMg.book(ZtreeMg);

    fEG->cd();
    TTree* ZtreeEG = new TTree("ZTree", "ZTree (egamma)");
    ZBranches brEG; brEG.MCweight = 1.0f; brEG.book(ZtreeEG);

    // Hadron proxy tree: mg only — used by analysis_mgHadron / FitJetFake
    fMg->cd();
    TTree* hadtree = new TTree("hadronTree", "hadronTree");
    float had_phoEt(0), had_phoEta(0), had_phoPhi(0), had_phoSigma(0), had_phoChIso(0);
    hadtree->Branch("phoEt",    &had_phoEt);
    hadtree->Branch("phoEta",   &had_phoEta);
    hadtree->Branch("phoPhi",   &had_phoPhi);
    hadtree->Branch("phoSigma", &had_phoSigma);
    hadtree->Branch("phoChIso", &had_phoChIso);

    // Shared object containers
    std::vector<recoPhoton> Photon;
    std::vector<recoMuon>   Muon;
    std::vector<recoEle>    Ele;
    std::vector<recoJet>    Jets;
    float MET(0), METPhi(0);
    int   nVtx(0), METFilter(0);

    const bool applyGolden = (true);  // golden JSON lumi mask for data
    if(applyGolden && !loadGolden(RunYear)) return;
    // =========================================================================
    // LOOP A — mg channel: MuonEG data,  Z -> mu mu + gamma
    // =========================================================================
    if (runMg) {
        rawData raw(esMg, rtMg);
        const unsigned N = esMg->GetEntries();
        int nTotal(0), npassHLT(0), npassPho(0), npassLep(0), npassZ(0);

        for (unsigned iev = 0; iev < N; ++iev) {

            if (iev % 100000 == 0) {
                std::cout << "[mg] event " << iev << "\n";
                log       << "[mg] event " << iev << "\n";
            }

            raw.GetData(esMg, iev);
            if(applyGolden && !passGolden(raw.run, raw.lumis)) continue;
            Photon.clear(); Muon.clear(); Ele.clear(); Jets.clear();
            for (int i = 0; i < raw.nPho; i++) Photon.push_back(recoPhoton(raw, i));
            for (int i = 0; i < raw.nMu;  i++) Muon.push_back(recoMuon(raw, i));
            for (int i = 0; i < raw.nEle; i++) Ele.push_back(recoEle(raw, i));
            for (int i = 0; i < raw.nJet; i++) Jets.push_back(recoJet(raw, i));
            MET = raw.pfMET; METPhi = raw.pfMETPhi;
            METFilter = raw.metFilters; nVtx = raw.nVtx;

            nTotal++;
            if (!raw.passHLT())    continue;   // MuonEG HLT
            if (raw.nGoodVtx < 1) continue;
            npassHLT++;

            // Pre-filter: exactly 2 medium mini-iso muons
            int nGoodMu = 0;
            for (auto& mu : Muon)
                if (mu.isMedium() && mu.getMiniIso() < 0.2 && mu.getPt() > 15) nGoodMu++;
            if (nGoodMu != 2 || raw.nPho < 1) continue;

            // Photon selection
            std::vector<recoPhoton>::iterator sigPho = Photon.begin(), hadPho = Photon.begin();
            bool hasSigPho(false), hasHadPho(false);
            selectPhoton(Photon, Ele, Muon, sigPho, hadPho, hasSigPho, hasHadPho);

            bool hasLep = false;
            std::vector<recoMuon>::iterator sigLep = Muon.begin();
            for (auto it = Muon.begin(); it != Muon.end(); ++it) {
                if (it->getPt() < 20)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasLep) { hasLep = true; sigLep = it; npassLep++; }
            }

            // ---- Signal ZTree fill ------------------------------------------
            if (hasSigPho && hasLep) {
                float dRlp = DeltaR(sigPho->getEta(), sigPho->getPhi(),
                                    sigLep->getEta(), sigLep->getPhi());
                if (dRlp > 0.8 && raw.passMETFilter(RunYear,METFilter)) {
                    bool foundZ = false;
                    std::vector<recoMuon>::iterator trailLep = Muon.begin();
                    for (auto im = Muon.begin(); im != Muon.end(); ++im) {
                        if (im == sigLep) continue;
                        if (!im->isMedium() || im->getMiniIso() > 0.2) continue;
                        double mll = (im->getP4() + sigLep->getP4()).M();
                        if (mll > 80 && mll < 100) { foundZ = true; trailLep = im; }
                    }
                    if (foundZ) {
                        float dp        = DeltaPhi(sigLep->getPhi(), METPhi);
                        brMg.phoEt      = sigPho->getCalibEt();
                        brMg.phoEta     = sigPho->getEta();
                        brMg.phoPhi     = sigPho->getPhi();
                        brMg.lepPt      = sigLep->getPt();
                        brMg.lepEta     = sigLep->getEta();
                        brMg.lepPhi     = sigLep->getPhi();
                        brMg.trailPt    = trailLep->getPt();
                        brMg.trailEta   = trailLep->getEta();
                        brMg.trailPhi   = trailLep->getPhi();
                        brMg.sigMT      = sqrt(2*MET*sigLep->getPt()*(1-cos(dp)));
                        brMg.sigMET     = MET;
                        brMg.sigMETPhi  = METPhi;
                        brMg.dPhiLepMET = dp;
                        brMg.threeMass  = (trailLep->getP4() + sigLep->getP4()
                                           + sigPho->getCalibP4()).M();
                        brMg.dilepMass  = (trailLep->getP4() + sigLep->getP4()).M();
                        brMg.nVertex    = nVtx;
                        brMg.dRPhoLep   = dRlp;
                        brMg.ISRJetPt   = calcISRJetPt(Jets,
                                              sigLep->getEta(),   sigLep->getPhi(),
                                              trailLep->getEta(), trailLep->getPhi(),
                                              sigPho->getEta(),   sigPho->getPhi(),
                                              brMg.HT, brMg.nJet);
                        fMg->cd(); ZtreeMg->Fill(); npassPho++; npassZ++;
                    }
                }
            }

            // ---- Hadron proxy fill (mg only) --------------------------------
            if (hasHadPho && hasLep) {
                float dRlp = DeltaR(hadPho->getEta(), hadPho->getPhi(),
                                    sigLep->getEta(), sigLep->getPhi());
                if (dRlp > 0.8 && raw.passMETFilter(RunYear,METFilter)) {
                    bool foundZ = false;
                    std::vector<recoMuon>::iterator trailLep = Muon.begin();
                    for (auto im = Muon.begin(); im != Muon.end(); ++im) {
                        if (im == sigLep) continue;
                        if (!im->isMedium() || im->getMiniIso() > 0.2) continue;
                        double mll = (im->getP4() + sigLep->getP4()).M();
                        if (mll > 80 && mll < 100) { foundZ = true; trailLep = im; }
                    }
                    if (foundZ) {
                        had_phoEt    = hadPho->getCalibEt();
                        had_phoEta   = hadPho->getEta();
                        had_phoPhi   = hadPho->getPhi();
                        had_phoSigma = hadPho->getSigma();
                        had_phoChIso = hadPho->getChIso();
                        fMg->cd(); hadtree->Fill();
                    }
                }
            }
        } // end mg event loop

        log << "[mg] nTotal=" << nTotal << "  npassHLT=" << npassHLT
            << "  npassPho=" << npassPho << "  npassLep=" << npassLep
            << "  npassZ=" << npassZ << "\n";
    } // end mg scope

    // =========================================================================
    // LOOP B — egamma channel: DoubleEG/EGamma data,  Z -> ee + gamma
    // =========================================================================
    if (runEG) {
        rawData raw(esEG, rtEG);
        const unsigned N = esEG->GetEntries();
        int nTotal(0), npassHLT(0), npassPho(0), npassLep(0), npassZ(0);

        for (unsigned iev = 0; iev < N; ++iev) {

            if (iev % 100000 == 0) {
                std::cout << "[eg] event " << iev << "\n";
                log       << "[eg] event " << iev << "\n";
            }

            raw.GetData(esEG, iev);
            if(applyGolden && !passGolden(raw.run, raw.lumis)) continue;
            Photon.clear(); Muon.clear(); Ele.clear(); Jets.clear();
            for (int i = 0; i < raw.nPho; i++) Photon.push_back(recoPhoton(raw, i));
            for (int i = 0; i < raw.nMu;  i++) Muon.push_back(recoMuon(raw, i));
            for (int i = 0; i < raw.nEle; i++) Ele.push_back(recoEle(raw, i));
            for (int i = 0; i < raw.nJet; i++) Jets.push_back(recoJet(raw, i));
            MET = raw.pfMET; METPhi = raw.pfMETPhi;
            METFilter = raw.metFilters; nVtx = raw.nVtx;

            nTotal++;
            if (!raw.passHLT())    continue;   // DoubleEG / EGamma HLT
            if (raw.nGoodVtx < 1) continue;
            npassHLT++;

            // Pre-filter: exactly 2 signal electrons (as nGoodMu != 2)
            int nGoodEle = 0;
            for (auto& ele : Ele)
                if (ele.passSignalSelection() && ele.getPt() > 15) nGoodEle++;
            if (nGoodEle != 2 || raw.nPho < 1) continue;

            // Photon selection (same veto logic as mg)
            std::vector<recoPhoton>::iterator sigPho = Photon.begin(), hadPho = Photon.begin();
            bool hasSigPho(false), hasHadPho(false);
            selectPhoton(Photon, Ele, Muon, sigPho, hadPho, hasSigPho, hasHadPho);

            bool hasLep = false;
            std::vector<recoEle>::iterator sigEle = Ele.begin();
            for (auto it = Ele.begin(); it != Ele.end(); ++it) {
                if (it->getPt() < 25)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasLep) { hasLep = true; sigEle = it; npassLep++; }
            }

            // ---- Signal ZTree fill ------------------------------------------
            if (hasSigPho && hasLep) {
                float dRlp = DeltaR(sigPho->getEta(), sigPho->getPhi(),
                                    sigEle->getEta(), sigEle->getPhi());
                if (dRlp > 0.8 && raw.passMETFilter(RunYear,METFilter)) {
                    bool foundZ = false;
                    std::vector<recoEle>::iterator trailEle = Ele.begin();
                    for (auto ie = Ele.begin(); ie != Ele.end(); ++ie) {
                        if (ie == sigEle) continue;
                        if (!ie->passSignalSelection()) continue;
                        double mll = (ie->getP4() + sigEle->getP4()).M();
                        if (mll > 80 && mll < 100) { foundZ = true; trailEle = ie; }
                    }
                    if (foundZ) {
                        float dp        = DeltaPhi(sigEle->getPhi(), METPhi);
                        brEG.phoEt      = sigPho->getCalibEt();
                        brEG.phoEta     = sigPho->getEta();
                        brEG.phoPhi     = sigPho->getPhi();
                        brEG.lepPt      = sigEle->getPt();
                        brEG.lepEta     = sigEle->getEta();
                        brEG.lepPhi     = sigEle->getPhi();
                        brEG.trailPt    = trailEle->getPt();
                        brEG.trailEta   = trailEle->getEta();
                        brEG.trailPhi   = trailEle->getPhi();
                        brEG.sigMT      = sqrt(2*MET*sigEle->getPt()*(1-cos(dp)));
                        brEG.sigMET     = MET;
                        brEG.sigMETPhi  = METPhi;
                        brEG.dPhiLepMET = dp;
                        brEG.threeMass  = (trailEle->getP4() + sigEle->getP4()
                                           + sigPho->getCalibP4()).M();
                        brEG.dilepMass  = (trailEle->getP4() + sigEle->getP4()).M();
                        brEG.nVertex    = nVtx;
                        brEG.dRPhoLep   = dRlp;
                        brEG.ISRJetPt   = calcISRJetPt(Jets,
                                              sigEle->getEta(),   sigEle->getPhi(),
                                              trailEle->getEta(), trailEle->getPhi(),
                                              sigPho->getEta(),   sigPho->getPhi(),
                                              brEG.HT, brEG.nJet);
                        fEG->cd(); ZtreeEG->Fill(); npassPho++; npassZ++;
                    }
                }
            }
        } // end egamma event loop

        log << "[eg] nTotal=" << nTotal << "  npassHLT=" << npassHLT
            << "  npassPho=" << npassPho << "  npassLep=" << npassLep
            << "  npassZ=" << npassZ << "\n";
    } // end egamma scope

    fMg->cd(); fMg->Write("", TFile::kOverwrite); fMg->Close();
    fEG->cd(); fEG->Write("", TFile::kOverwrite); fEG->Close();
    delete esMg; delete esEG;
    log.close();
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <RunYear> <Era> [both|mg|egamma]\n";
        return 1;
    }
    analysis_ISR(atoi(argv[1]), argv[2], argc > 3 ? argv[3] : "both");
    return 0;
}
