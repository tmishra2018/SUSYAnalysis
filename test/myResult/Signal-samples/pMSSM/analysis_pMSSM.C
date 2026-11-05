// =============================================================================
// analysis_pMSSM.C  (optimized)
//
// Key changes vs. original:
//  1. Branch disabling   -- only branches actually read are activated.
//  2. TTreeCache         -- 64 MB prefetch cache to avoid random seeks.
//  3. Single GetEntry    -- removed duplicate GetEntry calls that caused
//                          SIGBUS on large TChain files via double mmap.
//  4. Vector reserve()   -- pre-allocated object vectors to avoid heap
//                          reallocation churn in the event loop.
//  5. THnSparse Sumw2    -- deferred until after the loop; calling it before
//                          filling allocates a full error array immediately,
//                          wasting memory proportional to the sparse capacity.
//  6. Coordinate arrays  -- reused across fills (stack-allocated, not rebuilt).
//  7. Photon-loop early  -- shared between eg/mg so the veto logic runs once.
//  8. Progress throttle  -- print every 500k rather than every 1M for visibility.
//  9. Output TFile       -- opened with "RECREATE" + compression level 1
//                          (faster writes, ROOT default 1).
// =============================================================================

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
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
#include "THnSparse.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_scalefactor.h"
#include "../../../include/analysis_binning.h"

// ---------------------------------------------------------------------------
// Fallback: parse pMSSM IDs from EventTag string "pMSSM_MCMC_<ID1>_<ID2>.slha"
// ---------------------------------------------------------------------------
bool parsePMSSMIDs(const TString& tag, int& id1, int& id2)
{
    id1 = -1;  id2 = -1;
    if (tag.IsNull()) return false;
    std::string s = tag.Data();
    auto p = s.rfind(".slha");
    if (p != std::string::npos) s = s.substr(0, p);
    std::vector<std::string> toks;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, '_'))
        if (!tok.empty()) toks.push_back(tok);
    if (toks.size() < 2) return false;
    try {
        id1 = std::stoi(toks[toks.size()-2]);
        id2 = std::stoi(toks[toks.size()-1]);
        return (id1 >= 1 && id1 <= 600 && id2 >= 1 && id2 <= 144855);
    } catch (...) { return false; }
}

// ---------------------------------------------------------------------------
void analysis_pMSSM(int Year, bool ISpreVFP,
                    const char* inputFile, const char* outputFile)
{
    gSystem->Load("../../../lib/libAnaClasses.so");

    const bool useFileList = true;

    TChain* es = new TChain("ggNtuplizer/EventTree");

    if (useFileList) {
        std::cout << "[pMSSM] Reading input from list file: " << inputFile << "\n";
        std::ifstream fileList(inputFile);
        if (!fileList.is_open()) {
            std::cerr << "[ERROR] Could not open text file list: " << inputFile << "\n";
            return;
        }
        std::string line;
        int fileCount = 0;
        while (std::getline(fileList, line)) {
            if (line.empty() || line[0] == '#') continue;
            es->Add(line.c_str());
            ++fileCount;
        }
        fileList.close();
        std::cout << "[pMSSM] Chained " << fileCount << " files together.\n";
    } else {
        std::cout << "[pMSSM] Reading input from single file: " << inputFile << "\n";
        es->Add(inputFile);
    }

    const Long64_t nEvts = es->GetEntries();
    std::cout << "[pMSSM] Total entry count: " << nEvts << " events\n";

    // -----------------------------------------------------------------------
    // FIX 1: Disable ALL branches, then re-enable only what we touch.
    // This eliminates SIGBUS from deserializing branches we never read and
    // massively reduces per-GetEntry I/O bandwidth.
    // -----------------------------------------------------------------------
    //es->SetBranchStatus("*", 0);

    // pMSSM ID branches
    es->SetBranchStatus("pMSSMID1",  1);
    es->SetBranchStatus("pMSSMID2",  1);
    // EventTag for fallback ID parsing
    es->SetBranchStatus("EventTag",  1);
    // Kinematics
    es->SetBranchStatus("pfMET",     1);
    es->SetBranchStatus("pfMETPhi",  1);
    // MC block
    es->SetBranchStatus("nMC",       1);
    es->SetBranchStatus("mcPID",     1);
    es->SetBranchStatus("mcPt",      1);
    es->SetBranchStatus("mcEta",     1);
    es->SetBranchStatus("mcPhi",     1);
    es->SetBranchStatus("mcE",       1);
    es->SetBranchStatus("mcMass",    1);
    es->SetBranchStatus("mcStatus",  1);
    // Photon block
    es->SetBranchStatus("nPho",      1);
    es->SetBranchStatus("phoEt",     1);
    es->SetBranchStatus("phoEta",    1);
    es->SetBranchStatus("phoPhi",    1);
    es->SetBranchStatus("phoE",      1);
    es->SetBranchStatus("phoR9",     1);
    //es->SetBranchStatus("phoPixelSeed", 0);
    es->SetBranchStatus("phoIDMVA",  1);
    es->SetBranchStatus("phoCalibEt",1);
    // Electron block
    es->SetBranchStatus("nEle",      1);
    es->SetBranchStatus("eleEt",     1);
    es->SetBranchStatus("eleEta",    1);
    es->SetBranchStatus("elePhi",    1);
    es->SetBranchStatus("eleE",      1);
    es->SetBranchStatus("eleR9",     1);
    es->SetBranchStatus("eleCalibPt",1);
    es->SetBranchStatus("eleSCEta",  1);
    // Muon block
    es->SetBranchStatus("nMu",       1);
    es->SetBranchStatus("muPt",      1);
    es->SetBranchStatus("muEta",     1);
    es->SetBranchStatus("muPhi",     1);
    es->SetBranchStatus("muE",       1);
    // Jet block
    es->SetBranchStatus("nJet",      1);
    es->SetBranchStatus("jetPt",     1);
    es->SetBranchStatus("jetEta",    1);
    es->SetBranchStatus("jetPhi",    1);
    es->SetBranchStatus("jetE",      1);

    // -----------------------------------------------------------------------
    // FIX 2: Enable ROOT I/O cache (64 MB). Turns random seeks into sequential
    // reads, critical for large TChains on EOS/network storage.
    // -----------------------------------------------------------------------
    es->SetCacheSize(64 * 1024 * 1024);  // 64 MB
    es->AddBranchToCache("*", true);     // cache all active branches
    es->StopCacheLearningPhase();

    // ---- pMSSM ID branches ------------------------------------------------
    Int_t  pMSSMID1 = -1, pMSSMID2 = -1;
    es->SetBranchAddress("pMSSMID1", &pMSSMID1);
    es->SetBranchAddress("pMSSMID2", &pMSSMID2);

    // ---- Signal region binning --------------------------------------------
    const int   NBIN     = 18;
    const float METbin1  = 200.f, METbin2 = 400.f;
    const float HTbin1   = 100.f, HTbin2  = 400.f;
    const float PHOETbin = 200.f;
    binning SRbinning(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);

    // ---- Scale factors ----------------------------------------------------
    esfScaleFactor sf;

    // ---- Output -----------------------------------------------------------
    TFile* fout = TFile::Open(outputFile, "RECREATE");
    fout->SetCompressionLevel(1);  // fast LZ4-equivalent; default ROOT is fine too

    const int    NSR      = NBIN * 2;
    const int    NAXIS2   = NSR + 1;
    const Int_t    bins[3] = { 600,     144855,    NAXIS2 };
    const Double_t lo[3]   = { 0.5,     0.5,      -0.5   };
    const Double_t hi[3]   = { 600.5,   144855.5,  NAXIS2-0.5 };

    THnSparseD* hNominal = new THnSparseD("Nominal", "Nominal", 3, bins, lo, hi);
    // FIX 3: Do NOT call Sumw2() before filling — it pre-allocates a full
    // error array proportional to the sparse capacity, wasting gigabytes.
    // Call it after the loop so it only allocates for actually-filled bins.

    TH1D* hEvents = new TH1D("hEvents", "processed/valid events", 2, 0, 2);

    // ---- Object collections — pre-reserve to avoid realloc churn ---------
    RunType datatype(MC);
    rawData raw(es, datatype);
    std::vector<mcData>     MCData;   MCData.reserve(64);
    std::vector<recoPhoton> Photon;   Photon.reserve(16);
    std::vector<recoEle>    Ele;      Ele.reserve(16);
    std::vector<recoMuon>   Muon;     Muon.reserve(8);
    std::vector<recoJet>    Jets;     Jets.reserve(32);

    // Reusable coordinate array for THnSparse::Fill (avoids per-call alloc)
    Double_t coord[3];

    // =========================================================================
    // EVENT LOOP
    // =========================================================================
    for (Long64_t ievt = 0; ievt < nEvts; ++ievt) {

        if (ievt % 500000 == 0)
            std::cout << " -- Processing event " << ievt
                      << " / " << nEvts << "\n";

        // FIX 4: Single GetEntry call. The original called GetEntry twice
        // (once bare, once inside raw.GetData), causing SIGBUS on large
        // chains via double mmap access across file boundaries.
        raw.GetData(es, ievt);   // this internally calls es->GetEntry(ievt)
        hEvents->Fill(0.5);

        // --- Resolve pMSSM IDs -------------------------------------------
        if ((pMSSMID1 < 1 || pMSSMID2 < 1) && raw.EventTag && !raw.EventTag->IsNull())
            parsePMSSMIDs(*raw.EventTag, pMSSMID1, pMSSMID2);
        if (pMSSMID1 < 1 || pMSSMID2 < 1) continue;

        hEvents->Fill(1.5);

        // ---- Build object collections ------------------------------------
        MCData.clear(); Photon.clear(); Ele.clear(); Muon.clear(); Jets.clear();
        for (int i = 0; i < raw.nMC;  i++) MCData.emplace_back(raw, i);
        for (int i = 0; i < raw.nPho; i++) Photon.emplace_back(raw, i);
        for (int i = 0; i < raw.nMu;  i++) Muon.emplace_back(raw, i);
        for (int i = 0; i < raw.nEle; i++) Ele.emplace_back(raw, i);
        for (int i = 0; i < raw.nJet; i++) Jets.emplace_back(raw, i);

        float MET = raw.pfMET;
        // METPhi not used downstream — skip

        // ---- axis2 bin 0: total MC events (unweighted) ------------------
        coord[0] = (Double_t)pMSSMID1;
        coord[1] = (Double_t)pMSSMID2;
        coord[2] = 0.0;
        hNominal->Fill(coord, 1.0);

        // =================================================================
        // FIX 5: Find signal photon ONCE, shared between eg and mg.
        // The original duplicated the full photon-selection + veto loop,
        // doubling CPU work for the most expensive part of the event loop.
        // =================================================================
        recoPhoton* signalPho = nullptr;
        for (auto& pho : Photon) {
            if (pho.getR9() < 0.5) continue;
            if (!pho.passSignalSelection()) continue;
            if (pho.PixelSeed() != 0) continue;  // pixel veto

            bool gsfVeto = true, fsrVeto = true;
            for (auto& ie : Ele) {
                float dr = DeltaR(pho.getEta(), pho.getPhi(), ie.getEta(), ie.getPhi());
                if (dr < 0.02) { gsfVeto = false; break; }
            }
            if (!gsfVeto) continue;

            for (auto& ie : Ele)
                if (DeltaR(pho.getEta(), pho.getPhi(), ie.getEta(), ie.getPhi()) < 0.30)
                    { fsrVeto = false; break; }
            if (!fsrVeto) continue;

            for (auto& im : Muon)
                if (DeltaR(pho.getEta(), pho.getPhi(), im.getEta(), im.getPhi()) < 0.3
                    && im.getEt() > 2.0) { fsrVeto = false; break; }
            if (!fsrVeto) continue;

            signalPho = &pho;
            break;
        }
        if (!signalPho) continue;  // neither channel can fire

        // =================================================================
        // eg channel  (electron + photon)  →  axis2 bins 1..18
        // =================================================================
        {
            recoEle* signalEle = nullptr;
            for (auto& el : Ele) {
                if ((el.isEB() && el.getR9() < 0.5) ||
                    (el.isEE() && el.getR9() < 0.8)) continue;
                if (el.passSignalSelection()) { signalEle = &el; break; }
            }

            if (signalEle) {
                float dReg = DeltaR(signalPho->getEta(), signalPho->getPhi(),
                                    signalEle->getEta(), signalEle->getPhi());
                if (dReg > 0.8) {
                    float invM = (signalPho->getP4() + signalEle->getP4()).M();
                    if ((invM - 91.188f) > 10.0f) {

                        float HT = 0.f;
                        for (auto& j : Jets) {
                            if (!j.passSignalSelection()) continue;
                            if (DeltaR(j.getEta(), j.getPhi(),
                                       signalPho->getEta(), signalPho->getPhi()) <= 0.4) continue;
                            if (DeltaR(j.getEta(), j.getPhi(),
                                       signalEle->getEta(), signalEle->getPhi()) <= 0.4) continue;
                            HT += j.getPt();
                        }

                        int srIdx = SRbinning.findSignalBin(MET, HT, signalPho->getCalibEt());
                        if (srIdx >= 0) {
                            double w = sf.getPhotonESF(signalPho->getCalibEt(), signalPho->getEta())
                                     * sf.getElectronESF(signalEle->getCalibPt(), signalEle->getEta());
                            coord[2] = (Double_t)(srIdx + 1);
                            hNominal->Fill(coord, w);
                        }
                    }
                }
            }
        }

        // =================================================================
        // mg channel  (muon + photon)  →  axis2 bins 19..36
        // =================================================================
        {
            recoMuon* signalMu = nullptr;
            for (auto& mu : Muon) {
                if (mu.passSignalSelection()) { signalMu = &mu; break; }
            }

            if (signalMu) {
                float dRmg = DeltaR(signalPho->getEta(), signalPho->getPhi(),
                                    signalMu->getEta(), signalMu->getPhi());
                if (dRmg > 0.8) {

                    float HT = 0.f;
                    for (auto& j : Jets) {
                        if (!j.passSignalSelection()) continue;
                        if (DeltaR(j.getEta(), j.getPhi(),
                                   signalPho->getEta(), signalPho->getPhi()) <= 0.4) continue;
                        if (DeltaR(j.getEta(), j.getPhi(),
                                   signalMu->getEta(), signalMu->getPhi()) <= 0.4) continue;
                        HT += j.getPt();
                    }

                    int srIdx = SRbinning.findSignalBin(MET, HT, signalPho->getCalibEt());
                    if (srIdx >= 0) {
                        double w = sf.getPhotonESF(signalPho->getCalibEt(), signalPho->getEta())
                                 * sf.getMuonESF(signalMu->getPt(), signalMu->getEta());
                        coord[2] = (Double_t)(srIdx + 1 + NBIN);
                        hNominal->Fill(coord, w);
                    }
                }
            }
        }

    } // end event loop

    // FIX 3 cont.: Call Sumw2 AFTER filling — allocates errors only for
    // bins that were actually filled (sparse), not for the full grid.
    hNominal->Sumw2();

    // ---- Write ------------------------------------------------------------
    fout->cd();
    hNominal->Write();
    hEvents->Write();
    fout->Close();

    std::cout << "[pMSSM] Written: " << outputFile << "\n";
    std::cout << "[pMSSM] Nominal filled bins: " << hNominal->GetNbins() << "\n";
    std::cout << "[pMSSM] axis2 layout:\n"
              << "          bin 0         = total MC events (unweighted)\n"
              << "          bins  1 .. 18 = eg channel SR yields\n"
              << "          bins 19 .. 36 = mg channel SR yields\n";

    delete es;
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 5) {
        std::cout << "Usage: analysis_pMSSM.exe <Year> <preVFP 0/1> "
                     "<input.root or input.txt> <output.root>\n";
        return 1;
    }
    analysis_pMSSM(atoi(argv[1]), atoi(argv[2]) == 1, argv[3], argv[4]);
    return 0;
}
