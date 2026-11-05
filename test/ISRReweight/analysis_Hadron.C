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
#include "../../include/analysis_jet.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_mcData.h"

// ---------------------------------------------------------------------------
// Helper: book the three common trees (signal / eleproxy / had)
// ---------------------------------------------------------------------------
struct ProxyTreeSet {
    // signalTree branches
    float sig_phoEt, sig_phoEta, sig_phoPhi, sig_phoSigma, sig_phoChIso;
    float sig_sigMT, sig_sigMET, sig_sigMETPhi, sig_dPhiLepMET;
    int   sig_nVertex;
    float sig_HT, sig_nJet, sig_dilepMass;

    // eleproxyTree branches
    float ep_phoEt, ep_phoEta, ep_phoPhi, ep_phoSigma, ep_phoChIso;
    float ep_sigMT, ep_sigMET, ep_sigMETPhi, ep_dPhiLepMET;
    int   ep_nVertex;
    float ep_HT, ep_nJet, ep_dilepMass;

    // hadTree branches
    float had_phoEt, had_phoEta, had_phoPhi, had_phoSigma, had_phoChIso;
    float had_sigMT, had_sigMET, had_sigMETPhi, had_dPhiLepMET;
    int   had_nVertex;
    float had_HT, had_nJet, had_dilepMass;

    TTree *sigTree, *epTree, *hadTree;

    void book(TFile* f, int mcType) {
        f->cd();

        sigTree = new TTree("signalTree", "signalTree");
        sigTree->Branch("phoEt",     &sig_phoEt);
        sigTree->Branch("phoEta",    &sig_phoEta);
        sigTree->Branch("phoPhi",    &sig_phoPhi);
        sigTree->Branch("phoSigma",  &sig_phoSigma);
        sigTree->Branch("phoChIso",  &sig_phoChIso);
        sigTree->Branch("sigMT",     &sig_sigMT);
        sigTree->Branch("sigMET",    &sig_sigMET);
        sigTree->Branch("sigMETPhi", &sig_sigMETPhi);
        sigTree->Branch("dPhiLepMET",&sig_dPhiLepMET);
        sigTree->Branch("nVertex",   &sig_nVertex);
        sigTree->Branch("HT",        &sig_HT);
        sigTree->Branch("nJet",      &sig_nJet);
        sigTree->Branch("dilepMass", &sig_dilepMass);
        sigTree->Branch("mcType",    new int(mcType));   // constant branch

        epTree = new TTree("eleproxyTree", "eleproxyTree");
        epTree->Branch("phoEt",     &ep_phoEt);
        epTree->Branch("phoEta",    &ep_phoEta);
        epTree->Branch("phoPhi",    &ep_phoPhi);
        epTree->Branch("phoSigma",  &ep_phoSigma);
        epTree->Branch("phoChIso",  &ep_phoChIso);
        epTree->Branch("sigMT",     &ep_sigMT);
        epTree->Branch("sigMET",    &ep_sigMET);
        epTree->Branch("sigMETPhi", &ep_sigMETPhi);
        epTree->Branch("dPhiLepMET",&ep_dPhiLepMET);
        epTree->Branch("nVertex",   &ep_nVertex);
        epTree->Branch("HT",        &ep_HT);
        epTree->Branch("nJet",      &ep_nJet);
        epTree->Branch("dilepMass", &ep_dilepMass);

        hadTree = new TTree("hadTree", "hadTree");
        hadTree->Branch("phoEt",     &had_phoEt);
        hadTree->Branch("phoEta",    &had_phoEta);
        hadTree->Branch("phoPhi",    &had_phoPhi);
        hadTree->Branch("phoSigma",  &had_phoSigma);
        hadTree->Branch("phoChIso",  &had_phoChIso);
        hadTree->Branch("sigMT",     &had_sigMT);
        hadTree->Branch("sigMET",    &had_sigMET);
        hadTree->Branch("sigMETPhi", &had_sigMETPhi);
        hadTree->Branch("dPhiLepMET",&had_dPhiLepMET);
        hadTree->Branch("nVertex",   &had_nVertex);
        hadTree->Branch("HT",        &had_HT);
        hadTree->Branch("nJet",      &had_nJet);
        hadTree->Branch("dilepMass", &had_dilepMass);
    }
};

// ---------------------------------------------------------------------------
// Fill one of the three proxy trees from a photon iterator + kinematics.
// tree = 0 -> signal, 1 -> eleproxy, 2 -> had
// ---------------------------------------------------------------------------
static void fillProxyTree(ProxyTreeSet& ts, int tree,
                           std::vector<recoPhoton>::iterator pho,
                           float MT, float MET, float METPhi, float dp,
                           int nVtx, float HT, int nJet, float dilepMass)
{
    auto fill = [&](float& et, float& eta, float& phi, float& sigma, float& chiso,
                    float& mt_, float& met_, float& metphi_, float& dphi_,
                    int& vtx_, float& ht_, float& nj_, float& mll_, TTree* t)
    {
        et = pho->getCalibEt(); eta = pho->getEta(); phi = pho->getPhi();
        sigma = pho->getSigma(); chiso = pho->getChIso();
        mt_ = MT; met_ = MET; metphi_ = METPhi; dphi_ = dp;
        vtx_ = nVtx; ht_ = HT; nj_ = (float)nJet; mll_ = dilepMass;
        t->Fill();
    };
    if (tree == 0) fill(ts.sig_phoEt,  ts.sig_phoEta,  ts.sig_phoPhi,  ts.sig_phoSigma,  ts.sig_phoChIso,
                        ts.sig_sigMT,  ts.sig_sigMET,  ts.sig_sigMETPhi, ts.sig_dPhiLepMET,
                        ts.sig_nVertex, ts.sig_HT, ts.sig_nJet, ts.sig_dilepMass, ts.sigTree);
    if (tree == 1) fill(ts.ep_phoEt,   ts.ep_phoEta,   ts.ep_phoPhi,   ts.ep_phoSigma,   ts.ep_phoChIso,
                        ts.ep_sigMT,   ts.ep_sigMET,   ts.ep_sigMETPhi,  ts.ep_dPhiLepMET,
                        ts.ep_nVertex, ts.ep_HT, ts.ep_nJet, ts.ep_dilepMass, ts.epTree);
    if (tree == 2) fill(ts.had_phoEt,  ts.had_phoEta,  ts.had_phoPhi,  ts.had_phoSigma,  ts.had_phoChIso,
                        ts.had_sigMT,  ts.had_sigMET,  ts.had_sigMETPhi, ts.had_dPhiLepMET,
                        ts.had_nVertex, ts.had_HT, ts.had_nJet, ts.had_dilepMass, ts.hadTree);
}

// ---------------------------------------------------------------------------
// Z-pair finder for muons: returns true if a second muon makes 80<mll<100
// Same Z(ll) window and trailing-lepton cuts as analysis_ISR.C (FitJetFake.C
// keeps only 80<mll<100, so a pair found via the mllg window would be lost)
// ---------------------------------------------------------------------------
template <typename LepIter>
bool findZMuon(std::vector<recoMuon>& Muon, LepIter sigLep,
               std::vector<recoPhoton>::iterator pho,
               std::vector<recoMuon>::iterator& trailLep)
{
    for (auto im = Muon.begin(); im != Muon.end(); ++im) {
        if (im == sigLep) continue;
        if (!im->isMedium() || im->getMiniIso() > 0.2) continue;
        double mll  = (im->getP4() + sigLep->getP4()).M();
        if (mll > 80 && mll < 100) {
            trailLep = im; return true;
        }
    }
    return false;
}

// Z-pair finder for electrons
template <typename LepIter>
bool findZEle(std::vector<recoEle>& Ele, LepIter sigEle,
              std::vector<recoPhoton>::iterator pho,
              std::vector<recoEle>::iterator& trailEle)
{
    for (auto ie = Ele.begin(); ie != Ele.end(); ++ie) {
        if (ie == sigEle) continue;
        if (!ie->passSignalSelection()) continue;
        double mll  = (ie->getP4() + sigEle->getP4()).M();
        if (mll > 80 && mll < 100) {
            trailEle = ie; return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void analysis_Hadron(int RunYear, const char* Era, const char* channelArg = "mg") {

    gSystem->Load("../../lib/libAnaClasses.so");

    const std::string channel = channelArg;
    const bool isMG = (channel == "mg");

    // ---- RunType -----------------------------------------------------------
    RunType datatype;
    if (isMG) {
        if (RunYear == 2016) datatype = MuonEG2016;
        else if (RunYear == 2017) datatype = MuonEG2017;
        else                      datatype = MuonEG2018;
    } else {
        if (RunYear == 2016) datatype = DoubleEG2016;
        else if (RunYear == 2017) datatype = DoubleEG2017;
        else                      datatype = DoubleEG2018;
    }

    // ---- Log ---------------------------------------------------------------
    std::ofstream log(Form(
        "/eos/uscms/store/user/tmishra/jetfakepho/logs/plot_hadron_ISR_%s_%d%s.log",
        channelArg, RunYear, Era), std::ios::trunc);
    log << "analysis_Hadron  channel=" << channel
        << "  year=" << RunYear << "  era=" << Era
        << "  RunType=" << datatype << "\n";

    // ---- Input chain -------------------------------------------------------
    TChain* es = new TChain("ggNtuplizer/EventTree");
    if (isMG) {
        es->Add(Form(
            "/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_%d%s.root",
            RunYear, Era));
    } else {
            es->Add(Form(
                "/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_%d%s.root",
                RunYear, Era));
    }

    const unsigned nEvts = es->GetEntries();
    log   << "total entries: " << nEvts << "\n";
    std::cout << "[" << channel << "] total entries: " << nEvts << "\n";

    // ---- Output file -------------------------------------------------------
    TFile* outfile = TFile::Open(Form(
        "/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_ISR_%s_%d%s.root",
        channelArg, RunYear, Era), "RECREATE");

    ProxyTreeSet ts;
    ts.book(outfile, MCType::NOMC);

    // ---- Objects -----------------------------------------------------------
    rawData raw(es, datatype);
    std::vector<recoPhoton> Photon;
    std::vector<recoMuon>   Muon;
    std::vector<recoEle>    Ele;
    float MET(0), METPhi(0);
    int   nVtx(0), METFilter(0);
    int   nJet(0);

    // ---- Event loop --------------------------------------------------------
    const bool applyGolden = (true);  // golden JSON lumi mask for data
    if(applyGolden && !loadGolden(RunYear)) return;
    for (unsigned iev = 0; iev < nEvts; ++iev) {

        if (iev % 100000 == 0) std::cout << "[" << channel << "] event " << iev << "\n";

        raw.GetData(es, iev);
        if(applyGolden && !passGolden(raw.run, raw.lumis)) continue;
        Photon.clear(); Muon.clear(); Ele.clear();
        for (int i = 0; i < raw.nPho; i++) Photon.push_back(recoPhoton(raw, i));
        for (int i = 0; i < raw.nMu;  i++) Muon.push_back(recoMuon(raw, i));
        for (int i = 0; i < raw.nEle; i++) Ele.push_back(recoEle(raw, i));
        MET = raw.pfMET; METPhi = raw.pfMETPhi;
        METFilter = raw.metFilters; nVtx = raw.nVtx; nJet = raw.nJet;

        if (!raw.passHLT())    continue;
        if (raw.nGoodVtx < 1) continue;
        if (raw.nPho < 1)     continue;

        // ---- Classify photon candidates ------------------------------------
        // sigPho     : passes full signal selection, GSF+pixel+FSR vetoes
        // eleproxyVec: fails pixel OR GSF veto (electron-faking-photon proxy)
        // hadPhoVec  : passes basic+loose (before signal iso/sigma), GSF+pixel+FSR vetoes
        std::vector<recoPhoton>::iterator sigPho;
        bool hasSigPho = false;
        std::vector<std::vector<recoPhoton>::iterator> eleproxyVec, hadPhoVec;
        eleproxyVec.clear(); hadPhoVec.clear();

        for (auto it = Photon.begin(); it != Photon.end(); ++it) {
            if (it->getR9()             < 0.5)  continue;
            if (!it->passHLTSelection())          continue;
            if (!it->passBasicSelection())        continue;
            if (it->getCalibEt()        < 25)    continue;
            if (it->getChIso()          > 20.0)  continue;

            bool pixVeto = (it->PixelSeed() == 0);
            bool gsfVeto(true), fsrVeto(true);
            for (auto& ele : Ele) {
                if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <= 0.02) gsfVeto = false;
                if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <  0.3)  fsrVeto = false;
            }
            for (auto& mu : Muon)
                if (DeltaR(it->getEta(), it->getPhi(), mu.getEta(), mu.getPhi()) < 0.3 && mu.getPt() > 2.0)
                    fsrVeto = false;

            // Electron proxy: fails pixel veto or GSF veto
            if (!pixVeto || !gsfVeto) eleproxyVec.push_back(it);

            // Hadron proxy: passes all vetoes (loose, before signal cuts).
            // signalTree photon stays loose too: FitJetFake.C needs the
            // sigma_ietaieta sideband as its fake template.
            if (gsfVeto && pixVeto && fsrVeto) {
                hadPhoVec.push_back(it);
                if (!hasSigPho) { hasSigPho = true; sigPho = it; }
            }
        }

        // ====================================================================
        // mg channel: muon lepton
        // ====================================================================
        if (isMG) {
            // Pre-filter: exactly 2 medium mini-iso muons
            int nGoodMu = 0;
            for (auto& mu : Muon)
                if (mu.isMedium() && mu.getMiniIso() < 0.2 && mu.getPt() > 15) nGoodMu++;
            if (nGoodMu != 2) continue;

            // Leading muon: pT > 20, HLT matched (original uses 20 for this script)
            bool hasMu = false;
            std::vector<recoMuon>::iterator sigMu = Muon.begin();
            for (auto it = Muon.begin(); it != Muon.end(); ++it) {
                if (it->getPt() < 20)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasMu) { hasMu = true; sigMu = it; }
            }
            if (!hasMu) continue;

            // Signal tree
            if (hasSigPho) {
                float dR = DeltaR(sigPho->getEta(), sigPho->getPhi(), sigMu->getEta(), sigMu->getPhi());
                if (dR > 0.8) {
                    std::vector<recoMuon>::iterator trailLep = Muon.begin();
                    if (findZMuon(Muon, sigMu, sigPho, trailLep)) {
                        float dp  = DeltaPhi(sigMu->getPhi(), METPhi);
                        float MT  = sqrt(2*MET*sigMu->getPt()*(1-cos(dp)));
                        float mll = (trailLep->getP4() + sigMu->getP4()).M();
                        fillProxyTree(ts, 0, sigPho, MT, MET, METPhi, dp, nVtx, 0, nJet, mll);
                    }
                }
            }

            // Electron proxy tree
            for (auto& epPho : eleproxyVec) {
                float dR = DeltaR(epPho->getEta(), epPho->getPhi(), sigMu->getEta(), sigMu->getPhi());
                if (dR > 0.8) {
                    std::vector<recoMuon>::iterator trailLep = Muon.begin();
                    if (findZMuon(Muon, sigMu, epPho, trailLep)) {
                        float mll = (trailLep->getP4() + sigMu->getP4()).M();
                        fillProxyTree(ts, 1, epPho, 0, MET, METPhi, 0, nVtx, 0, nJet, mll);
                    }
                }
            }

            // Hadron proxy tree
            for (auto& hPho : hadPhoVec) {
                float dR = DeltaR(hPho->getEta(), hPho->getPhi(), sigMu->getEta(), sigMu->getPhi());
                if (dR > 0.8) {
                    std::vector<recoMuon>::iterator trailLep = Muon.begin();
                    if (findZMuon(Muon, sigMu, hPho, trailLep)) {
                        float mll = (trailLep->getP4() + sigMu->getP4()).M();
                        fillProxyTree(ts, 2, hPho, 0, MET, METPhi, 0, nVtx, 0, nJet, mll);
                    }
                }
            }

        // ====================================================================
        // egamma channel: electron lepton
        // ====================================================================
        } else {
            // Pre-filter: exactly 2 signal electrons (as nGoodMu != 2)
            int nGoodEle = 0;
            for (auto& ele : Ele)
                if (ele.passSignalSelection() && ele.getPt() > 15) nGoodEle++;
            if (nGoodEle != 2) continue;

            // Leading signal electron: pT > 25 (AN Sec. 3), HLT matched
            bool hasEle = false;
            std::vector<recoEle>::iterator sigEle = Ele.begin();
            for (auto it = Ele.begin(); it != Ele.end(); ++it) {
                if (it->getPt() < 25)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasEle) { hasEle = true; sigEle = it; }
            }
            if (!hasEle) continue;

            // Signal tree
            if (hasSigPho) {
                float dR = DeltaR(sigPho->getEta(), sigPho->getPhi(), sigEle->getEta(), sigEle->getPhi());
                if (dR > 0.8) {
                    std::vector<recoEle>::iterator trailEle = Ele.begin();
                    if (findZEle(Ele, sigEle, sigPho, trailEle)) {
                        float dp  = DeltaPhi(sigEle->getPhi(), METPhi);
                        float MT  = sqrt(2*MET*sigEle->getPt()*(1-cos(dp)));
                        float mll = (trailEle->getP4() + sigEle->getP4()).M();
                        fillProxyTree(ts, 0, sigPho, MT, MET, METPhi, dp, nVtx, 0, nJet, mll);
                    }
                }
            }

            // Electron proxy tree
            for (auto& epPho : eleproxyVec) {
                float dR = DeltaR(epPho->getEta(), epPho->getPhi(), sigEle->getEta(), sigEle->getPhi());
                if (dR > 0.8) {
                    std::vector<recoEle>::iterator trailEle = Ele.begin();
                    if (findZEle(Ele, sigEle, epPho, trailEle)) {
                        float mll = (trailEle->getP4() + sigEle->getP4()).M();
                        fillProxyTree(ts, 1, epPho, 0, MET, METPhi, 0, nVtx, 0, nJet, mll);
                    }
                }
            }

            // Hadron proxy tree
            for (auto& hPho : hadPhoVec) {
                float dR = DeltaR(hPho->getEta(), hPho->getPhi(), sigEle->getEta(), sigEle->getPhi());
                if (dR > 0.8) {
                    std::vector<recoEle>::iterator trailEle = Ele.begin();
                    if (findZEle(Ele, sigEle, hPho, trailEle)) {
                        float mll = (trailEle->getP4() + sigEle->getP4()).M();
                        fillProxyTree(ts, 2, hPho, 0, MET, METPhi, 0, nVtx, 0, nJet, mll);
                    }
                }
            }
        }
    } // event loop

    outfile->Write();
    delete es;
    log.close();
    std::cout << "[DONE] " << channel << " " << RunYear << Era << "\n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <RunYear> <Era> <channel mg|egamma>\n";
        return 1;
    }
    analysis_Hadron(atoi(argv[1]), argv[2], argv[3]);
    return 0;
}
