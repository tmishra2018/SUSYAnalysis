// analysis_ISRMC.C
// Universal MC skimmer for ISR JetPt weight derivation.
// Invoked once per channel; the channel argument controls:
//   - RunType  -> HLT path emulation (MuonEG vs DoubleEG/EGamma)
//   - Lumi     -> getLumi(runtype) returns the matching integrated lumi
//   - Lepton   -> muon (mg) or electron (egamma)
//
// Output:
//   mg     -> resTree_ISR_{sample}_mg_{year}{vfp}.root
//   egamma -> resTree_ISR_{sample}_egamma_{year}{vfp}.root
//
// Build:
//   g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISRMC.C \
//       -o analysis_ISRMC.exe `root-config --libs`
// Usage:
//   ./analysis_ISRMC.exe <RunYear> <preVFP 0|1> <SampleName> <channel mg|egamma>
//   e.g.  ./analysis_ISRMC.exe 2018 0 ZGToLLG mg
//         ./analysis_ISRMC.exe 2018 0 ZGToLLG egamma
//         ./analysis_ISRMC.exe 2016 1 TTGJets  egamma

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
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"

// ---------------------------------------------------------------------------
// Map sample name -> MCType enum (defined in analysis_mcData.h)
// ---------------------------------------------------------------------------
static int getMCType(const std::string& s) {
    if (s == "ZGToLLG")  return MCType::ZGInclusive;
    if (s == "TTGJets")  return MCType::TTG;
    if (s == "TTJets")   return MCType::TT;
    if (s == "WWG")      return MCType::WWG;
    if (s == "WZG")      return MCType::WZG;
    std::cerr << "[WARN] Unknown sample: " << s << ". Defaulting to ZGInclusive.\n";
    return MCType::ZGInclusive;
}

// ---------------------------------------------------------------------------
// Channel-dependent RunType.
//   mg     -> MCMuonEG{year}   — HLT paths and lumi for MuonEG dataset
//   egamma -> MCDoubleEG{year} / MCEGamma2018 — for DoubleEG / EGamma dataset
// ---------------------------------------------------------------------------
static RunType getRunType(const std::string& channel, int RunYear) {
    if (channel == "mg") {
        if (RunYear == 2016) return MCMuonEG2016;
        if (RunYear == 2017) return MCMuonEG2017;
        return MCMuonEG2018;
    } else {
        if (RunYear == 2016) return MCDoubleEG2016;
        if (RunYear == 2017) return MCDoubleEG2017;
        return MCDoubleEG2018;   // 2018: EGamma dataset
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void analysis_ISRMC(int RunYear, int preVFP,
                    const char* sampleName,
                    const char* channelArg = "mg") {

    gSystem->Load("../../lib/libAnaClasses.so");

    const std::string sample  = sampleName;
    const std::string channel = channelArg;
    const bool        isMG    = (channel == "mg");

    const std::string vfpTag = (RunYear == 2016)
                               ? (preVFP ? "preVFP" : "postVFP") : "";

    // ---- RunType: sets HLT emulation + lumi lookup -------------------------
    const int     mcType  = getMCType(sample);
    const RunType runtype = getRunType(channel, RunYear);

    // ---- Input chain -------------------------------------------------------
    TChain* es = new TChain("ggNtuplizer/EventTree");
    char* inputfile = new char[300];
    if (strstr(sampleName, "DYJetsToLL") != NULL or strstr(sampleName, "TTJets") != NULL or strstr(sampleName, "WJetsToLNu") or strstr(sampleName, "ZGToLLG") != NULL )
        sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/%s/%s_%d%s.root",sampleName,sampleName,RunYear,vfpTag.c_str());
    else
        sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",sampleName,sampleName,RunYear,vfpTag.c_str());
    es->Add(inputfile);

    const unsigned nEvts = es->GetEntries();
    // NLO samples (ZGToLLG, TTGJets, TTJets, WWG, WZG) have negative generator weights:
    // each event counts with the sign of genWeight, normalised to N(+) - N(-)
    const long long nNeg  = es->GetEntries("genWeight < 0");
    const double    nEff  = double(nEvts) - 2.0*double(nNeg);
    const float crosssection       = MC_XS[mcType];
    // lumi of the dataset used in that channel (as Background/analysis_VGBkg.C)
    float MCweight;
    if(isMG){
        if(RunYear == 2016 and preVFP == 1)             MCweight = lumi_2016preVFP_MuonEG*1000*crosssection/nEff;
        else if(RunYear == 2016 and preVFP == 0)        MCweight = lumi_2016postVFP_MuonEG*1000*crosssection/nEff;
        else if(RunYear == 2017)                        MCweight = lumi_2017_MuonEG*1000*crosssection/nEff;
        else if(RunYear == 2018)                        MCweight = lumi_2018_MuonEG*1000*crosssection/nEff;
    } else {
        if(RunYear == 2016 and preVFP == 1)             MCweight = lumi_2016preVFP_DoubleEG*1000*crosssection/nEff;
        else if(RunYear == 2016 and preVFP == 0)        MCweight = lumi_2016postVFP_DoubleEG*1000*crosssection/nEff;
        else if(RunYear == 2017)                        MCweight = lumi_2017_DoubleEG*1000*crosssection/nEff;
        else if(RunYear == 2018)                        MCweight = lumi_2018_DoubleEG*1000*crosssection/nEff;
    }


        cout << "crosssection = " << crosssection << std::endl;
        std::cout << "Total event: " << nEvts << ", negative-weight events: " << nNeg << std::endl;

    // ---- Log ---------------------------------------------------------------
    std::ofstream log(Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/logs/resTree_ISR_%s_%s_%d%s.log",
        sampleName, channelArg, RunYear, vfpTag.c_str()), std::ios::trunc);
    log << "analysis_ISRMC  sample=" << sample
        << "  channel=" << channel
        << "  RunType=" << runtype
        << "  year="    << RunYear << vfpTag
        << "  nEvts="   << nEvts
        << "  nNeg="    << nNeg
        << "  MCweight=" << MCweight << "\n";
    std::cout << "analysis_ISRMC  sample=" << sample
              << "  channel=" << channel
              << "  RunType=" << runtype
              << "  year="    << RunYear << vfpTag
              << "  nEvts="   << nEvts
              << "  MCweight=" << MCweight << "\n";

    // ---- Output file -------------------------------------------------------
    TFile* fOut = TFile::Open(Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/resTree_ISR_%s_%s_%d%s.root",
        sampleName, channelArg, RunYear, vfpTag.c_str()), "RECREATE");
    TTree* Ztree = new TTree("ZTree", Form("ZTree %s", channelArg));

    // ---- Branch variables --------------------------------------------------
    float phoEt(0), phoEta(0), phoPhi(0);
    float lepPt(0), lepEta(0), lepPhi(0);
    float trailPt(0), trailEta(0), trailPhi(0);
    float sigMT(0), sigMET(0), sigMETPhi(0), dPhiLepMET(0);
    float threeMass(0), dilepMass(0);
    int   nVertex(0);
    float dRPhoLep(0), HT(0), nJet(0), ISRJetPt(0), bosonPt(0);
    float PUweight(1.0f);
    float _MCweight = MCweight;
    int   _mcType   = mcType;
    std::vector<int>   mcPID, mcMomPID, mcGMomPID;
    std::vector<float> mcEtaV, mcPhiV, mcPtV;

    Ztree->Branch("MCweight",   &_MCweight);
    Ztree->Branch("PUweight",   &PUweight);
    Ztree->Branch("mcType",     &_mcType);
    Ztree->Branch("phoEt",      &phoEt);
    Ztree->Branch("phoEta",     &phoEta);
    Ztree->Branch("phoPhi",     &phoPhi);
    Ztree->Branch("lepPt",      &lepPt);
    Ztree->Branch("lepEta",     &lepEta);
    Ztree->Branch("lepPhi",     &lepPhi);
    Ztree->Branch("trailPt",    &trailPt);
    Ztree->Branch("trailEta",   &trailEta);
    Ztree->Branch("trailPhi",   &trailPhi);
    Ztree->Branch("sigMT",      &sigMT);
    Ztree->Branch("sigMET",     &sigMET);
    Ztree->Branch("sigMETPhi",  &sigMETPhi);
    Ztree->Branch("dPhiLepMET", &dPhiLepMET);
    Ztree->Branch("threeMass",  &threeMass);
    Ztree->Branch("dilepMass",  &dilepMass);
    Ztree->Branch("nVertex",    &nVertex);
    Ztree->Branch("dRPhoLep",   &dRPhoLep);
    Ztree->Branch("HT",         &HT);
    Ztree->Branch("nJet",       &nJet);
    Ztree->Branch("ISRJetPt",   &ISRJetPt);
    Ztree->Branch("bosonPt",    &bosonPt);
    Ztree->Branch("mcPID",      &mcPID);
    Ztree->Branch("mcMomPID",   &mcMomPID);
    Ztree->Branch("mcGMomPID",  &mcGMomPID);
    Ztree->Branch("mcEta",      &mcEtaV);
    Ztree->Branch("mcPhi",      &mcPhiV);
    Ztree->Branch("mcPt",       &mcPtV);

    // ---- True-pileup weights (as VGamma/src/analysis_VGamma_Run2.C) --------
    // data nTrueInt profile (golden JSON, 69.2 mb) over the nTrueInt (BX==0) profile of this sample
    TH1D* h_mc_pu = new TH1D("h_mc_pu", "MC nTrueInt", 200, 0, 200);
    h_mc_pu->SetDirectory(0);
    es->Project("h_mc_pu", "puTrue", "(puBX==0)*(1 - 2*(genWeight < 0))");
    TFile* f_pu = TFile::Open(Form("../PUfiles/dataPU_%d%s.root", RunYear, vfpTag.c_str()));
    if (!f_pu || f_pu->IsZombie()) {
        std::cerr << "missing ../PUfiles/dataPU_" << RunYear << vfpTag << ".root" << std::endl;
        return;
    }
    TH1D* h_puRatio = (TH1D*)((TH1D*)f_pu->Get("pileup"))->Clone("h_puRatio");
    h_puRatio->SetDirectory(0);
    f_pu->Close();
    fOut->cd();
    h_puRatio->Scale(1.0 / h_puRatio->Integral());
    h_mc_pu->Scale(1.0 / h_mc_pu->Integral());
    h_puRatio->Divide(h_mc_pu);
    log << "true-PU weights from ../PUfiles/dataPU_" << RunYear << vfpTag
        << ".root, MC nTrueInt entries " << h_mc_pu->GetEntries() << "\n";

    // ---- rawData: RunType controls HLT + lumi ------------------------------
    rawData raw(es, runtype);
    float genWeight(1);
    es->SetBranchAddress("genWeight", &genWeight);

    std::vector<mcData>     MCData;
    std::vector<recoPhoton> Photon;
    std::vector<recoMuon>   Muon;
    std::vector<recoEle>    Ele;
    std::vector<recoJet>    Jets;
    float MET(0), METPhi(0);
    int   nVtx(0), METFilter(0);

    int nTotal(0), npassHLT(0), npassPho(0), npassLep(0), npassZ(0);

    // ---- Event loop --------------------------------------------------------
    for (unsigned iev = 0; iev < nEvts; ++iev) {

        if (iev % 100000 == 0) {
            std::cout << "[" << channel << "] event " << iev << "\n";
            log       << "[" << channel << "] event " << iev << "\n";
        }

        raw.GetData(es, iev);
        _MCweight = (genWeight < 0) ? -MCweight : MCweight;
        MCData.clear(); Photon.clear(); Muon.clear(); Ele.clear(); Jets.clear();
        for (int i = 0; i < raw.nMC;  i++) MCData.push_back(mcData(raw, i));
        for (int i = 0; i < raw.nPho; i++) Photon.push_back(recoPhoton(raw, i));
        for (int i = 0; i < raw.nMu;  i++) Muon.push_back(recoMuon(raw, i));
        for (int i = 0; i < raw.nEle; i++) Ele.push_back(recoEle(raw, i));
        for (int i = 0; i < raw.nJet; i++) Jets.push_back(recoJet(raw, i));
        MET = raw.pfMET; METPhi = raw.pfMETPhi;
        METFilter = raw.metFilters; nVtx = raw.nVtx;

        nTotal++;
        if (!raw.passHLT())    continue;   // HLT path set by runtype
        if (raw.nGoodVtx < 1) continue;
        npassHLT++;

        PUweight = 1.0f;
        if (raw.puBX && raw.puTrue) {
            for (size_t ibx = 0; ibx < raw.puBX->size(); ibx++) {
                if ((*raw.puBX)[ibx] != 0) continue;
                int b = std::max(1, std::min(h_puRatio->FindBin((*raw.puTrue)[ibx]), h_puRatio->GetNbinsX()));
                if (h_puRatio->GetBinContent(b) > 0) PUweight = h_puRatio->GetBinContent(b);
                break;
            }
        }

        // Generator-level boson pT
        bosonPt = 0;
        for (auto& mc : MCData) {
            int pid = abs(mc.getPID());
            if (pid == 23 || pid == 24) { bosonPt = mc.getPt(); break; }
        }

        // ---- Photon selection (identical for both channels) ----------------
        bool hasSigPho(false);
        std::vector<recoPhoton>::iterator sigPho = Photon.begin();
        for (auto it = Photon.begin(); it != Photon.end(); ++it) {
            if (it->getR9()             < 0.5) continue;
            if (!it->passHLTSelection())        continue;
            if (!it->passBasicSelection())      continue;
            if (!it->passSignalSelection())     continue;
            bool pixVeto = (it->PixelSeed() == 0);
            bool gsfVeto(true), fsrVeto(true);
            for (auto& ele : Ele) {
                if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <= 0.02) gsfVeto = false;
                if (DeltaR(it->getEta(), it->getPhi(), ele.getEta(), ele.getPhi()) <  0.3)  fsrVeto = false;
            }
            for (auto& mu : Muon)
                if (DeltaR(it->getEta(), it->getPhi(), mu.getEta(), mu.getPhi()) < 0.3 && mu.getPt() > 2.0)
                    fsrVeto = false;
            if (gsfVeto && pixVeto && fsrVeto && !hasSigPho) {
                hasSigPho = true; sigPho = it; npassPho++;
            }
        }
        if (!hasSigPho) continue;

        // ---- MC truth matching helper --------------------------------------
        auto fillMCVecs = [&](float lepEta_, float lepPhi_) {
            mcPID.clear(); mcMomPID.clear(); mcGMomPID.clear();
            mcEtaV.clear(); mcPhiV.clear(); mcPtV.clear();
            for (auto& mc : MCData) {
                if (mc.getEt() < 1.0) continue;
                if (DeltaR(sigPho->getEta(), sigPho->getPhi(), mc.getEta(), mc.getPhi()) < 0.3 ||
                    DeltaR(lepEta_, lepPhi_, mc.getEta(), mc.getPhi()) < 0.3) {
                    mcPID.push_back(mc.getPID());
                    mcMomPID.push_back(mc.getMomPID());
                    mcGMomPID.push_back(mc.getGMomPID());
                    mcEtaV.push_back(mc.getEta());
                    mcPhiV.push_back(mc.getPhi());
                    mcPtV.push_back(mc.getEt());
                }
            }
        };

        // ---- ISR jet pT helper --------------------------------------------
        // jets cleaned against both leptons and the photon (dR > 0.4), as in analysis_ISR.C
        auto calcISR = [&](float l1Eta, float l1Phi,
                           float l2Eta, float l2Phi) -> float {
            TLorentzVector jsys(0,0,0,0);
            HT = 0; nJet = 0;
            for (auto& jet : Jets) {
                if (!jet.passSignalSelection()) continue;
                if (DeltaR(jet.getEta(), jet.getPhi(), l1Eta, l1Phi) <= 0.4) continue;
                if (DeltaR(jet.getEta(), jet.getPhi(), l2Eta, l2Phi) <= 0.4) continue;
                if (DeltaR(jet.getEta(), jet.getPhi(), sigPho->getEta(), sigPho->getPhi()) <= 0.4) continue;
                nJet++; HT += jet.getPt(); jsys += jet.getP4();
            }
            return jsys.Pt();
        };

        // ==================================================================
        // mg channel: mu + gamma,  Z -> mu mu
        // ==================================================================
        if (isMG) {
            // Pre-filter: exactly 2 medium mini-iso muons
            int nGoodMu = 0;
            for (auto& mu : Muon)
                if (mu.isMedium() && mu.getMiniIso() < 0.2 && mu.getPt() > 15) nGoodMu++;
            if (nGoodMu != 2) continue;

            bool hasLep = false;
            std::vector<recoMuon>::iterator sigLep = Muon.begin();
            for (auto it = Muon.begin(); it != Muon.end(); ++it) {
                if (it->getPt() < 20)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasLep) { hasLep = true; sigLep = it; npassLep++; }
            }
            if (!hasLep) continue;

            float dRlp = DeltaR(sigPho->getEta(), sigPho->getPhi(),
                                sigLep->getEta(), sigLep->getPhi());
            if (dRlp <= 0.8 || !raw.passMETFilter(RunYear,METFilter)) continue;

            bool foundZ = false;
            std::vector<recoMuon>::iterator trailLep = Muon.begin();
            for (auto im = Muon.begin(); im != Muon.end(); ++im) {
                if (im == sigLep) continue;
                if (!im->isMedium() || im->getMiniIso() > 0.2) continue;
                double mll = (im->getP4() + sigLep->getP4()).M();
                if (mll > 80 && mll < 100) { foundZ = true; trailLep = im; }
            }
            if (!foundZ) continue;

            float dp    = DeltaPhi(sigLep->getPhi(), METPhi);
            phoEt       = sigPho->getCalibEt();
            phoEta      = sigPho->getEta();
            phoPhi      = sigPho->getPhi();
            lepPt       = sigLep->getPt();
            lepEta      = sigLep->getEta();
            lepPhi      = sigLep->getPhi();
            trailPt     = trailLep->getPt();
            trailEta    = trailLep->getEta();
            trailPhi    = trailLep->getPhi();
            sigMT       = sqrt(2*MET*sigLep->getPt()*(1-cos(dp)));
            sigMET      = MET;
            sigMETPhi   = METPhi;
            dPhiLepMET  = dp;
            threeMass   = (trailLep->getP4() + sigLep->getP4() + sigPho->getCalibP4()).M();
            dilepMass   = (trailLep->getP4() + sigLep->getP4()).M();
            nVertex     = nVtx;
            dRPhoLep    = dRlp;
            ISRJetPt    = calcISR(sigLep->getEta(), sigLep->getPhi(),
                                   trailLep->getEta(), trailLep->getPhi());
            fillMCVecs(sigLep->getEta(), sigLep->getPhi());
            fOut->cd(); Ztree->Fill(); npassZ++;

        // ==================================================================
        // egamma channel: e + gamma,  Z -> ee
        // ==================================================================
        } else {
            // Pre-filter: exactly 2 signal electrons (as nGoodMu != 2)
            int nGoodEle = 0;
            for (auto& ele : Ele)
                if (ele.passSignalSelection() && ele.getPt() > 15) nGoodEle++;
            if (nGoodEle != 2) continue;

            bool hasLep = false;
            std::vector<recoEle>::iterator sigEle = Ele.begin();
            for (auto it = Ele.begin(); it != Ele.end(); ++it) {
                if (it->getPt() < 25)           continue;
                if (!it->passHLTSelection())     continue;
                if (it->passSignalSelection() && !hasLep) { hasLep = true; sigEle = it; npassLep++; }
            }
            if (!hasLep) continue;

            float dRlp = DeltaR(sigPho->getEta(), sigPho->getPhi(),
                                sigEle->getEta(), sigEle->getPhi());
            if (dRlp <= 0.8 || !raw.passMETFilter(RunYear,METFilter)) continue;

            bool foundZ = false;
            std::vector<recoEle>::iterator trailEle = Ele.begin();
            for (auto ie = Ele.begin(); ie != Ele.end(); ++ie) {
                if (ie == sigEle) continue;
                if (!ie->passSignalSelection()) continue;
                double mll = (ie->getP4() + sigEle->getP4()).M();
                if (mll > 80 && mll < 100) { foundZ = true; trailEle = ie; }
            }
            if (!foundZ) continue;

            float dp    = DeltaPhi(sigEle->getPhi(), METPhi);
            phoEt       = sigPho->getCalibEt();
            phoEta      = sigPho->getEta();
            phoPhi      = sigPho->getPhi();
            lepPt       = sigEle->getPt();
            lepEta      = sigEle->getEta();
            lepPhi      = sigEle->getPhi();
            trailPt     = trailEle->getPt();
            trailEta    = trailEle->getEta();
            trailPhi    = trailEle->getPhi();
            sigMT       = sqrt(2*MET*sigEle->getPt()*(1-cos(dp)));
            sigMET      = MET;
            sigMETPhi   = METPhi;
            dPhiLepMET  = dp;
            threeMass   = (trailEle->getP4() + sigEle->getP4() + sigPho->getCalibP4()).M();
            dilepMass   = (trailEle->getP4() + sigEle->getP4()).M();
            nVertex     = nVtx;
            dRPhoLep    = dRlp;
            ISRJetPt    = calcISR(sigEle->getEta(), sigEle->getPhi(),
                                   trailEle->getEta(), trailEle->getPhi());
            fillMCVecs(sigEle->getEta(), sigEle->getPhi());
            fOut->cd(); Ztree->Fill(); npassZ++;
        }

    } // event loop

    log << "DONE  nTotal=" << nTotal
        << "  npassHLT="   << npassHLT
        << "  npassPho="   << npassPho
        << "  npassLep="   << npassLep
        << "  npassZ="     << npassZ << "\n";

    fOut->cd(); fOut->Write("", TFile::kOverwrite); fOut->Close();
    delete es;
    log.close();
    std::cout << "[DONE] " << sample << " " << channel
              << " " << RunYear << vfpTag << "\n";
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <RunYear> <preVFP 0|1> <SampleName> <channel mg|egamma>\n";
        return 1;
    }
    analysis_ISRMC(atoi(argv[1]), atoi(argv[2]), argv[3], argv[4]);
    return 0;
}
