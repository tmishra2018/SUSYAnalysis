// ============================================================================
// closure_doubleFake_v2.C
//
// FULL CLOSURE TEST FOR THE DOUBLE-FAKE BACKGROUND
// =================================================
//
// PURPOSE
// -------
// The original closure_doubleFake.C answers ONE question:
//   "What fraction of signal-region events have both objects fake?"
//
// This script answers the MORE IMPORTANT question for a CMS analysis:
//   "Does our data-driven background estimation METHOD correctly predict
//    the double-fake yield, or does it double-count or miss it?"
//
// These are different questions. The first gives you a number (f_df ~ 10%).
// The second tells you whether your existing background estimates already
// implicitly absorb that 10% or not -- i.e., whether you need to subtract
// something or assign a systematic.
//
// WHAT THIS SCRIPT IMPLEMENTS
// ----------------------------
// The full closure test has five stages, all run on MC:
//
//   STAGE 1 — Truth tagging
//     Every event is gen-matched. Each reco photon and lepton is tagged as:
//       PHO_TRUE, PHO_EFAKE, PHO_JETFAKE
//       LEP_PROMPT, LEP_FAKE
//     Combinations define: DF_ejl (e->g + fake-l), DF_jjl (jet->g + fake-l),
//     SF_pho (fake-g + prompt-l), SF_lep (prompt-g + fake-l), BOTH_PROMPT.
//
//   STAGE 2 — Transfer factor derivation in MC (mimicking the data method)
//     The analysis estimates fake backgrounds using three separate methods:
//       a) e->gamma fake rate:   f_eg  = N(phoID pass) / N(phoID fail, e-like)
//          derived in a MC proxy region (inverted pixel seed / sigma_ietaieta)
//       b) jet->gamma fake rate: f_jg  = N(phoID pass) / N(phoID fail, not e-like)
//          derived after subtracting MC e->g contribution from the CR
//       c) jet->lepton SF:       SF_jl = N(tight lep) / N(loose lep)
//          derived in a MC QCD-enriched region
//     NOTE: This script approximates these from the MC truth counts since the
//     full proxy/jetTree/fakeLepTree ntuples would be needed for the exact
//     derivation. The key structural test (double-counting check) does not
//     require the exact fake rates.
//
//   STAGE 3 — Build sequential prediction in double-fake-enriched VRs
//     Two validation regions (VR) are defined, each enriched in one DF type:
//       VR_A: jet->g + jet->l enriched
//              (photon: fail charged isolation; lepton: fail tight isolation)
//       VR_B: e->g + jet->l enriched
//              (photon: has pixel seed; lepton: fail tight isolation)
//     In each VR, the sequential prediction is:
//       N_pred = N(SF_pho in VR) * f_jl  +  N(SF_lep in VR) * f_jg
//     This is what the data-driven method would predict for double fakes.
//
//   STAGE 4 — Compare prediction vs truth in VR
//     N_truth (gen-matched double fakes in VR)  vs  N_pred (sequential method)
//     Plot as function of MET, HT, pT(gamma), pT(lep), dR(gamma,lep).
//     A ratio N_pred/N_truth flat at 1.0 = perfect closure.
//
//   STAGE 5 — Decomposition test
//     Verify each double-fake event is absorbed by exactly ONE background
//     category, not zero (missed) and not two (double-counted).
//     For each DF event: check which prediction category would pick it up.
//
// RELATION TO closure_doubleFake.C (v1)
// --------------------------------------
// v1 measures f_df = N(DF)/N(SF) in the signal region. This is the effect size.
// v2 tests whether the background method correctly handles that effect.
// Both are needed for a complete systematic treatment:
//   - v1 result: "double fakes are ~10% of single fakes in the eg channel"
//   - v2 result: "the method over/under-predicts double fakes by X%, assign
//                 X% as systematic on [affected background component]"
//
// HOW TO RUN
// ----------
// root -l -b -q 'closure_doubleFake_v2.C++(1, 2018, 0, "/path/to/ntuples")'
//   arg1: channelType (1=eg, 2=mg)
//   arg2: runYear     (2016, 2017, 2018)
//   arg3: isPreVFP    (1=preVFP, 0=postVFP/other, only relevant for 2016)
//   arg4: basePath    (directory containing resTree_*.root files)
//
// ============================================================================

#include "../../include/analysis_commoncode.h"
#include <iomanip>
#include <ctime>
#include <sstream>
#include "TROOT.h"

static const char* SCRIPT_V2_VERSION = "v2-fullclosure-2024";

// ============================================================================
// Reuse truth classifiers from closure_doubleFake.C
// (copied verbatim — identical framework logic)
// ============================================================================

enum PhoCategory { PHO_TRUE, PHO_EFAKE, PHO_JETFAKE, PHO_UNKNOWN };
enum LepCategory { LEP_PROMPT, LEP_FAKE, LEP_UNKNOWN };

// Double-fake sub-type for decomposition test
enum DFType { DF_EFAKE_JLEP,   // e->g fake + jet->l fake
              DF_JFAKE_JLEP,   // jet->g fake + jet->l fake
              DF_OTHER,        // rare (e->g + prompt-l from wrong sample, etc.)
              NOT_DF };        // not a double fake

static bool DEBUG_V2 = false;

PhoCategory classifyPhoton(float phoEta, float phoPhi, float phoEt,
                            std::vector<int>   *mcPID,
                            std::vector<float> *mcEta,
                            std::vector<float> *mcPhi,
                            std::vector<float> *mcPt,
                            std::vector<int>   *mcMomPID,
                            int /*pixelSeed*/) {

  // Step 1: TRUE PHOTON — analysis_VGamma.C istruepho
  double mindRpho = 0.3;
  unsigned phoIndex = 0;
  for(unsigned iMC = 0; iMC < mcPID->size(); iMC++){
    double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
    double dE = fabs((*mcPt)[iMC] - phoEt) / phoEt;
    if(dR < mindRpho && dE < 0.5){ mindRpho = dR; phoIndex = iMC; }
  }
  if(mindRpho < 0.2){
    int pid    = (*mcPID)[phoIndex];
    int momPID = fabs((*mcMomPID)[phoIndex]);
    if(pid == 22 && (momPID <= 6  || momPID == 21 || momPID == 999 ||
                     momPID == 11 || momPID == 13  || momPID == 15  ||
                     momPID == 23 || momPID == 24))
      return PHO_TRUE;
  }

  // Step 2: e->gamma FAKE — closure_efakepho.C isFakePho
  double mindR_eFake = 0.3, minEledR = 3.0;
  unsigned matchIndex = 0;
  for(unsigned iMC = 0; iMC < mcPID->size(); iMC++){
    double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
    if(dR < mindR_eFake){ mindR_eFake = dR; matchIndex = iMC; }
    if(fabs((*mcPID)[iMC]) == 11 && dR < minEledR) minEledR = dR;
  }
  if(mindR_eFake < 0.1 &&
     ((*mcPID)[matchIndex] == 11 || (*mcPID)[matchIndex] == -11))
    return PHO_EFAKE;
  if(minEledR < 0.02) return PHO_EFAKE;

  // Step 3: jet->gamma FAKE — closure_jetfakepho.C isFakePho
  bool isFakePho = true;
  double mindR_jFake = 0.3;
  unsigned jMatchIndex = 0;
  for(unsigned iMC = 0; iMC < mcPID->size(); iMC++){
    double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
    double dE = fabs((*mcPt)[iMC] - phoEt) / phoEt;
    if(dR < mindR_jFake && dE < 0.5){ mindR_jFake = dR; jMatchIndex = iMC; }
  }
  if(mindR_jFake < 0.1){
    int pid    = fabs((*mcPID)[jMatchIndex]);
    int momPID = fabs((*mcMomPID)[jMatchIndex]);
    if((pid == 11 || pid == 22) &&
       (momPID <= 6  || momPID == 21 || momPID == 11 || momPID == 999 ||
        momPID == 13 || momPID == 15  || momPID == 24 || momPID == 23))
      isFakePho = false;
  }
  if(isFakePho) return PHO_JETFAKE;
  return PHO_UNKNOWN;
}

LepCategory classifyLepton(float lepEta, float lepPhi, float lepPt,
                            std::vector<int>   *mcPID,
                            std::vector<float> *mcEta,
                            std::vector<float> *mcPhi,
                            std::vector<float> *mcPt,
                            std::vector<int>   *mcMomPID) {
  double mindR = 0.1;
  int    bestIdx = -1;
  for(unsigned iMC = 0; iMC < mcPID->size(); iMC++){
    double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], lepEta, lepPhi);
    double dE = fabs((*mcPt)[iMC] - lepPt) / lepPt;
    if(dR < mindR && dE < 0.3){ mindR = dR; bestIdx = (int)iMC; }
  }
  if(bestIdx >= 0){
    int pid    = abs((*mcPID)[bestIdx]);
    int momPID = fabs((*mcMomPID)[bestIdx]);
    if((pid == 11 || pid == 13) &&
       (momPID == 24 || momPID == 23 || momPID == 15)) return LEP_PROMPT;
    return LEP_FAKE;
  }
  return LEP_FAKE;
}

// ============================================================================
// Helper: build per-object filtered gen-particle views
// The eg ntuple stores all gen particles; we pre-filter to dR<0.3 of each
// reco object to match the mg ntuple behaviour and the framework convention.
// ============================================================================
void buildFilteredViews(float phoEta, float phoPhi,
                        float lepEta, float lepPhi,
                        std::vector<int>   *mcPID_in,
                        std::vector<float> *mcEta_in,
                        std::vector<float> *mcPhi_in,
                        std::vector<float> *mcPt_in,
                        std::vector<int>   *mcMomPID_in,
                        // outputs:
                        std::vector<int>   &pho_pid,
                        std::vector<float> &pho_eta,
                        std::vector<float> &pho_phi,
                        std::vector<float> &pho_pt,
                        std::vector<int>   &pho_mom,
                        std::vector<int>   &lep_pid,
                        std::vector<float> &lep_eta,
                        std::vector<float> &lep_phi,
                        std::vector<float> &lep_pt,
                        std::vector<int>   &lep_mom){
  pho_pid.clear(); pho_eta.clear(); pho_phi.clear();
  pho_pt.clear();  pho_mom.clear();
  lep_pid.clear(); lep_eta.clear(); lep_phi.clear();
  lep_pt.clear();  lep_mom.clear();
  for(unsigned i = 0; i < mcPID_in->size(); i++){
    float dRp = DeltaR(phoEta, phoPhi, (*mcEta_in)[i], (*mcPhi_in)[i]);
    float dRl = DeltaR(lepEta, lepPhi, (*mcEta_in)[i], (*mcPhi_in)[i]);
    if(dRp < 0.3f){
      pho_pid.push_back((*mcPID_in)[i]);
      pho_eta.push_back((*mcEta_in)[i]);
      pho_phi.push_back((*mcPhi_in)[i]);
      pho_pt .push_back((*mcPt_in) [i]);
      pho_mom.push_back((*mcMomPID_in)[i]);
    }
    if(dRl < 0.3f){
      lep_pid.push_back((*mcPID_in)[i]);
      lep_eta.push_back((*mcEta_in)[i]);
      lep_phi.push_back((*mcPhi_in)[i]);
      lep_pt .push_back((*mcPt_in) [i]);
      lep_mom.push_back((*mcMomPID_in)[i]);
    }
  }
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================
void closure_doubleFake_v2(){
  SetRunConfig();
  setTDRStyle();
  gSystem->Load("../../lib/libAnaClasses.so");
  gROOT->SetBatch(kTRUE);
  std::string basePath = "/eos/uscms/store/user/tmishra/VGammaNEW/";

  int channelType = ichannel; // eg=1, mg=2

  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

  std::string chainname;
  if(channelType == 1) chainname = "egTree";
  else                 chainname = "mgTree";

  TChain *qcdTree = new TChain(chainname.c_str());

  //qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT50to100_%d%s.root",RunYear,whichVFP.c_str()));
  //qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT100to200_%d%s.root",RunYear,whichVFP.c_str()));
  //qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT200to300_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT300to500_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT500to700_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT700to1000_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT1000to1500_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT1500to2000_%d%s.root",RunYear,whichVFP.c_str()));
  qcdTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_QCD_HT2000toInf_%d%s.root",RunYear,whichVFP.c_str()));
  // Category B: single-fake reference samples
  // GJets: prompt photon + fake lepton
  TChain *gjetsTree = new TChain(chainname.c_str());
  gjetsTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_GJets_DoubleEM_%d%s.root",RunYear,whichVFP.c_str()));

  // W+jets: prompt lepton + fake photon
  TChain *wjetsTree = new TChain(chainname.c_str());
  wjetsTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));
  
  TChain *ttjetsTree = new TChain(chainname.c_str());
  ttjetsTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_TTJets_%d%s.root",RunYear,whichVFP.c_str()));

  printf("  QCD=%lld  GJets=%lld  WJets=%lld entries\n\n", qcdTree->GetEntries(), gjetsTree->GetEntries(), wjetsTree->GetEntries());

  struct TreeVars {
    float crosssection, ntotalevent, PUweight;
    int   mcType;
    float phoEt, phoEta, phoPhi;
    float lepPt, lepEta, lepPhi;
    float sigMT, sigMET, dRPhoLep, HT, nJet;
    int   nBJet, nVertex;
    float sigMETPhi, dPhiLepMET;
    // photon ID discriminating variables for VR definitions
    // (used to define the VRs without separate ntuple files)
    float phoSigmaIetaIeta;   // σ_ηη : <0.0105 barrel tight
    int   phoPixelSeed;       // 0=no pixel seed (prompt), 1=has seed (e-like)
    float phoChargedIso;      // charged isolation (inverted for VR_A)
    float lepRelIso;          // lepton relative isolation
    std::vector<int>   *mcPID    = 0;
    std::vector<float> *mcEta    = 0;
    std::vector<float> *mcPhi    = 0;
    std::vector<float> *mcPt     = 0;
    std::vector<int>   *mcMomPID = 0;
  };

  auto connectBranches = [](TChain *tree, TreeVars &v){
    tree->SetBranchAddress("crosssection",    &v.crosssection);
    tree->SetBranchAddress("ntotalevent",     &v.ntotalevent);
    tree->SetBranchAddress("PUweight",        &v.PUweight);
    tree->SetBranchAddress("mcType",          &v.mcType);
    tree->SetBranchAddress("phoEt",           &v.phoEt);
    tree->SetBranchAddress("phoEta",          &v.phoEta);
    tree->SetBranchAddress("phoPhi",          &v.phoPhi);
    tree->SetBranchAddress("lepPt",           &v.lepPt);
    tree->SetBranchAddress("lepEta",          &v.lepEta);
    tree->SetBranchAddress("lepPhi",          &v.lepPhi);
    tree->SetBranchAddress("sigMT",           &v.sigMT);
    tree->SetBranchAddress("sigMET",          &v.sigMET);
    tree->SetBranchAddress("sigMETPhi",       &v.sigMETPhi);
    tree->SetBranchAddress("dPhiLepMET",      &v.dPhiLepMET);
    tree->SetBranchAddress("dRPhoLep",        &v.dRPhoLep);
    tree->SetBranchAddress("HT",              &v.HT);
    tree->SetBranchAddress("nJet",            &v.nJet);
    tree->SetBranchAddress("nBJet",           &v.nBJet);
    tree->SetBranchAddress("nVertex",         &v.nVertex);
    tree->SetBranchAddress("mcPID",           &v.mcPID);
    tree->SetBranchAddress("mcEta",           &v.mcEta);
    tree->SetBranchAddress("mcPhi",           &v.mcPhi);
    tree->SetBranchAddress("mcPt",            &v.mcPt);
    tree->SetBranchAddress("mcMomPID",        &v.mcMomPID);
    // Optional branches — present in some ntuple versions
    if(tree->GetBranch("phoSigmaIetaIeta"))
      tree->SetBranchAddress("phoSigmaIetaIeta", &v.phoSigmaIetaIeta);
    if(tree->GetBranch("phoPixelSeed"))
      tree->SetBranchAddress("phoPixelSeed",      &v.phoPixelSeed);
    if(tree->GetBranch("phoChargedIso"))
      tree->SetBranchAddress("phoChargedIso",     &v.phoChargedIso);
    if(tree->GetBranch("lepRelIso"))
      tree->SetBranchAddress("lepRelIso",         &v.lepRelIso);
  };

  TreeVars vQCD, vGJets, vWJets, vTTJets;
  connectBranches(qcdTree,   vQCD);
  connectBranches(gjetsTree, vGJets);
  connectBranches(wjetsTree, vWJets);
  connectBranches(ttjetsTree, vTTJets);

  // For each of the three regions (SR, VR_A, VR_B) and each of the four
  // truth categories (DF_ejl, DF_jjl, SF_pho, SF_lep), book key observables.
  //
  // Naming: h_{region}_{category}_{observable}
  //   region   : SR, VRA (jet->g+jet->l enriched), VRB (e->g+jet->l enriched)
  //   category : truth (DF_ejl / DF_jjl / SFpho / SFlep / BothPrompt)
  //              OR pred (sequential prediction of double fake)
  //   observable: MET / MT / phoEt / lepPt / HT / dR

  // MET bin edges — fine for shape comparison
  const int nMET = 8;
  double METbins[nMET+1] = {0,25,50,75,100,125,150,200,300};

  // HT bins
  const int nHT = 7;
  double HTbins[nHT+1] = {0,100,200,300,500,700,1000,2000};

  // pT bins (photon Et or lepton pT)
  const int nPT = 7;
  double PTbins[nPT+1] = {0,35,50,75,100,150,200,500};

  // dR(g,l) bins
  const int ndR = 6;
  double dRbins[ndR+1] = {0.8,1.0,1.5,2.0,2.5,3.0,5.0};

  // Helper macro to book a set of histograms for one category in one region
  // Returns array of 5 histograms: [MET, MT, phoEt, lepPt, HT]
  auto bookSet = [&](const char* region, const char* cat) {
    std::vector<TH1D*> hv;
    auto name = [&](const char* obs) {
      return std::string("h_") + region + "_" + cat + "_" + obs;
    };
    hv.push_back(new TH1D(name("MET").c_str(),
                  (std::string(cat)+" MET;MET (GeV);Events").c_str(),
                  nMET, METbins));
    hv.push_back(new TH1D(name("MT").c_str(),
                  (std::string(cat)+" MT;M_{T} (GeV);Events").c_str(),
                  nMET, METbins));
    hv.push_back(new TH1D(name("phoEt").c_str(),
                  (std::string(cat)+" phoEt;E_{T}^{#gamma} (GeV);Events").c_str(),
                  nPT, PTbins));
    hv.push_back(new TH1D(name("lepPt").c_str(),
                  (std::string(cat)+" lepPt;p_{T}^{l} (GeV);Events").c_str(),
                  nPT, PTbins));
    hv.push_back(new TH1D(name("HT").c_str(),
                  (std::string(cat)+" HT;HT (GeV);Events").c_str(),
                  nHT, HTbins));
    hv.push_back(new TH1D(name("dR").c_str(),
                  (std::string(cat)+" dR;#DeltaR(#gamma,l);Events").c_str(),
                  ndR, dRbins));
    for(auto h : hv){ h->Sumw2(); h->SetDirectory(0); }
    return hv;
  };

  // Truth histograms in all three regions
  auto hSR_DF_ejl  = bookSet("SR",  "DF_ejl");   // e->g + jet->l (truth)
  auto hSR_DF_jjl  = bookSet("SR",  "DF_jjl");   // jet->g + jet->l (truth)
  auto hSR_SFpho   = bookSet("SR",  "SFpho");    // fake-g + prompt-l (truth)
  auto hSR_SF_ejl  = bookSet("SR",  "SF_ejl");   // e->g
  auto hSR_SF_jjl  = bookSet("SR",  "SF_jjl");   // jet->g
  auto hSR_SFlep   = bookSet("SR",  "SFlep");    // prompt-g + fake-l (truth)
  auto hSR_Both    = bookSet("SR",  "Both");     // both prompt (sanity)

  auto hVRA_DF_ejl = bookSet("VRA", "DF_ejl");
  auto hVRA_DF_jjl = bookSet("VRA", "DF_jjl");
  auto hVRA_SFpho  = bookSet("VRA", "SFpho");
  auto hVRA_SFlep  = bookSet("VRA", "SFlep");
  auto hVRA_Both   = bookSet("VRA", "Both");

  auto hVRB_DF_ejl = bookSet("VRB", "DF_ejl");
  auto hVRB_DF_jjl = bookSet("VRB", "DF_jjl");
  auto hVRB_SFpho  = bookSet("VRB", "SFpho");
  auto hVRB_SFlep  = bookSet("VRB", "SFlep");
  auto hVRB_Both   = bookSet("VRB", "Both");

  // Sequential prediction histograms in VR_A and VR_B
  // Prediction = N(SFpho in VR) * f_jl  +  N(SFlep in VR) * f_jg
  // These are filled after the event loop using the truth histogram integrals
  // (see Stage 3 computation below)
  // For now, book prediction histograms that mirror VR_A truth
  auto hVRA_pred   = bookSet("VRA", "Pred");   // sequential prediction in VRA
  auto hVRB_pred   = bookSet("VRB", "Pred");   // sequential prediction in VRB

  // SR
  double nSR_DF_ejl(0), nSR_DF_jjl(0), nSR_SFpho(0), nSR_SFlep(0), nSR_Both(0), nSR_SF_ejl(0), nSR_SF_jjl(0);
  // VR_A (jet->g + jet->l enriched: photon fails charged iso)
  double nVRA_DF_ejl(0), nVRA_DF_jjl(0), nVRA_SFpho(0), nVRA_SFlep(0), nVRA_Both(0);
  // VR_B (e->g + jet->l enriched: photon has pixel seed)
  double nVRB_DF_ejl(0), nVRB_DF_jjl(0), nVRB_SFpho(0), nVRB_SFlep(0), nVRB_Both(0);

  // Decomposition counters:
  //   For each DF event in SR, track which background method would claim it.
  //   "Absorbed by method X" means the event satisfies the selection criteria
  //   of method X's control region definition.
  double nDF_absorbedByEfakePho(0);    // counted by analysis_eleBkg.C method
  double nDF_absorbedByJetfakePho(0);  // counted by analysis_jetBkg.C method
  double nDF_absorbedByJetfakeLep(0);  // counted by analysis_qcdBkg.C method
  double nDF_absorbedByMultiple(0);    // double-counted across methods
  double nDF_absorbedByNone(0);        // not absorbed by any method (missed)

  auto fillSet = [](std::vector<TH1D*> &hv, float MET, float MT,
                    float phoEt, float lepPt, float HT, float dR, double w){
    hv[0]->Fill(MET,   w);
    hv[1]->Fill(MT,    w);
    hv[2]->Fill(phoEt, w);
    hv[3]->Fill(lepPt, w);
    hv[4]->Fill(HT,    w);
    hv[5]->Fill(dR,    w);
  };

  // VR_A: jet->g enriched photon side (invert charged isolation)
  //       tight lepton fails tight iso → defined by lepRelIso > loose threshold
  // VR_B: e->g enriched photon side (require pixel seed on photon)
  //       same lepton sideband as VR_A
  //
  // IMPORTANT: since phoChargedIso, phoPixelSeed, lepRelIso are OPTIONAL
  // branches (may not be in all ntuple versions), we use truth-based proxies
  // when the branches are absent:
  //   VR_A proxy: truth-tag says PHO_JETFAKE (i.e., reco photon IS a jet fake)
  //               + truth-tag says LEP_FAKE
  //   VR_B proxy: truth-tag says PHO_EFAKE
  //               + truth-tag says LEP_FAKE
  // This is equivalent from the MC perspective: we are selecting exactly
  // the events we want to study in each VR, using truth instead of a reco cut.
  // When the full reco branches are available, swap to reco-based definitions.

  // Process all three samples together; weight by XS*lumi/nevents*PU

  auto processTree = [&](TChain *tree, TreeVars &v, const std::string &label){

    printf("  Processing %s (%lld entries)...\n", label.c_str(), tree->GetEntries());

    std::vector<int>   pho_pid, lep_pid;
    std::vector<float> pho_eta, lep_eta;
    std::vector<float> pho_phi, lep_phi;
    std::vector<float> pho_pt,  lep_pt;
    std::vector<int>   pho_mom, lep_mom;

    for(long long ievt = 0; ievt < tree->GetEntries(); ievt++){
      tree->GetEntry(ievt);

      // Standard kinematic selection — identical to all background scripts
      if(v.phoEt < 35 || fabs(v.phoEta) > 1.4442) continue;
      if(v.sigMET < lowMET)  continue;
      if(highMET > 0 && v.sigMET > highMET) continue;
      if(v.sigMT  < lowMt)   continue;
      if(highMt > 0 && v.sigMT > highMt)   continue;
      if(v.lepPt  < lowPt)   continue;
      if(highPt > 0 && v.lepPt > highPt)   continue;
      if(v.dRPhoLep < 0.8) continue;

      // Lumi weight only — no object SFs for fake closure study
      // SFs are derived for prompt objects (Z->ee, Z->mumu) and are not
      // meaningful for fake leptons/photons from QCD jets
      float XS_weight = 1.0;
      if(channelType == 1){
             if(RunYear == 2016 and preVFP == 1)                XS_weight = lumi_2016preVFP_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2016 and preVFP == 0)           XS_weight = lumi_2016postVFP_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2017)                           XS_weight = lumi_2017_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2018)                           XS_weight = lumi_2018_DoubleEG*1000*v.crosssection/v.ntotalevent;}

      else if(channelType == 2){
             if(RunYear == 2016 and preVFP == 1)                XS_weight = lumi_2016preVFP_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2016 and preVFP == 0)           XS_weight = lumi_2016postVFP_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2017)                           XS_weight = lumi_2017_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2018)                           XS_weight = lumi_2018_MuonEG*1000*v.crosssection/v.ntotalevent;}
      float w = v.PUweight * XS_weight;


      buildFilteredViews(v.phoEta, v.phoPhi, v.lepEta, v.lepPhi,
                         v.mcPID, v.mcEta, v.mcPhi, v.mcPt, v.mcMomPID,
                         pho_pid, pho_eta, pho_phi, pho_pt, pho_mom,
                         lep_pid, lep_eta, lep_phi, lep_pt, lep_mom);

      PhoCategory phoCat = classifyPhoton(v.phoEta, v.phoPhi, v.phoEt,
                                          &pho_pid, &pho_eta, &pho_phi,
                                          &pho_pt,  &pho_mom, 0);
      LepCategory lepCat = classifyLepton(v.lepEta, v.lepPhi, v.lepPt,
                                          &lep_pid, &lep_eta, &lep_phi,
                                          &lep_pt,  &lep_mom);

      bool phoIsEFake  = (phoCat == PHO_EFAKE);
      bool phoIsJFake  = (phoCat == PHO_JETFAKE);
      bool phoIsFake   = (phoIsEFake || phoIsJFake);
      bool phoIsTrue   = (phoCat == PHO_TRUE);
      bool lepIsFake   = (lepCat == LEP_FAKE);
      bool lepIsPrompt = (lepCat == LEP_PROMPT);

      bool isDF_ejl = phoIsEFake && lepIsFake;   // e->g + jet->l
      bool isDF_jjl = phoIsJFake && lepIsFake;   // jet->g + jet->l

      bool isSF_ejl = phoIsEFake && lepIsPrompt;   // e->g
      bool isSF_jjl = phoIsJFake && lepIsPrompt;   // jet->g

      bool isSFpho  = phoIsFake  && lepIsPrompt; // single fake photon
      bool isSFlep  = phoIsTrue  && lepIsFake;   // single fake lepton
      bool isBoth   = phoIsTrue  && lepIsPrompt; // both prompt (sanity)
      bool isDF     = isDF_ejl || isDF_jjl;

      // SR = inclusive (all ntuple events). VRs are defined by truth category only.
      // No kinematic sub-windows: the resTree ntuples already have MET>lowMET
      // applied upstream, and we cannot reopen a wider sideband from here.
      // The closure test is therefore self-consistency within the same event sample.

      float dR = v.dRPhoLep;

      if(isDF_ejl){ fillSet(hSR_DF_ejl, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_DF_ejl+=w; }
      if(isDF_jjl){ fillSet(hSR_DF_jjl, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_DF_jjl+=w; }
      if(isSFpho) { fillSet(hSR_SFpho,  v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_SFpho +=w; }
      if(isSF_ejl) { fillSet(hSR_SFpho,  v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_SF_ejl +=w; }
      if(isSF_jjl) { fillSet(hSR_SFpho,  v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_SF_jjl +=w; }
      if(isSFlep) { fillSet(hSR_SFlep,  v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_SFlep +=w; }
      if(isBoth)  { fillSet(hSR_Both,   v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nSR_Both  +=w; }

      // VRA = DF_jjl truth region; VRB = DF_ejl truth region
      // (inVRA and inVRB are exclusive by definition)
      if(isDF_jjl){ fillSet(hVRA_DF_jjl, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRA_DF_jjl+=w; }
      if(isDF_ejl){ fillSet(hVRB_DF_ejl, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRB_DF_ejl+=w; }
      // SF and Both counts are the same for both VRs (same inclusive sample)
      if(isSFpho){ fillSet(hVRA_SFpho,v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRA_SFpho+=w;
                   fillSet(hVRB_SFpho,v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRB_SFpho+=w; }
      if(isSFlep){ fillSet(hVRA_SFlep,v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRA_SFlep+=w;
                   fillSet(hVRB_SFlep,v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRB_SFlep+=w; }
      if(isBoth) { fillSet(hVRA_Both, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRA_Both +=w;
                   fillSet(hVRB_Both, v.sigMET,v.sigMT,v.phoEt,v.lepPt,v.HT,dR,w); nVRB_Both +=w; }

      // For each DF event, which background method absorbs it?
      //   e->g method (proxyTree):   absorbs DF_ejl — reco photon came from
      //                              electron, pixel seed present → in proxyTree
      //   jet->g method (jetTree):   absorbs DF_jjl — hadronic photon fake
      //                              lands in σ_ηη sideband → in jetTree
      //   jet->l method (fakeLepTree): requires hasPho=true (REAL photon).
      //                              DF events have fake photon → NOT in fakeLepTree.
      //                              → jet->l method structurally CANNOT absorb DF.
      // Expected: each DF event absorbed by exactly one method, none missed.
      if(isDF){
        bool absorbedByEfake = phoIsEFake;  // e->g method
        bool absorbedByJfake = phoIsJFake;  // jet->g method
        // absorbedByJlep = false by design (see above)
        int nAbs = (int)absorbedByEfake + (int)absorbedByJfake;
        if     (nAbs == 0) nDF_absorbedByNone      += w;
        else if(nAbs == 1){
          if(absorbedByEfake) nDF_absorbedByEfakePho  += w;
          else                nDF_absorbedByJetfakePho += w;
        }
        else                 nDF_absorbedByMultiple   += w;
      }
    } // end event loop
  }; // end processTree lambda

  printf("\n");
  processTree(qcdTree,   vQCD,   "QCD");
  processTree(gjetsTree, vGJets, "GJets");
  processTree(wjetsTree, vWJets, "WJets");
  processTree(ttjetsTree, vTTJets, "TTJets");

  // ============================================================================
  // STAGE 3: Build sequential prediction
  //
  // Since VR_A is defined as PHO_JETFAKE+LEP_FAKE and VR_B as PHO_EFAKE+LEP_FAKE,
  // these ARE the DF categories themselves. The closure test asks:
  //
  //   "If the background method applied a fake-lepton rate f_jl to the
  //    SF_pho events, would it predict the correct DF yield?"
  //
  // Transfer factors from the INCLUSIVE sample (SR = all events):
  //   f_jl = N(DF_jjl) / N(SF_pho)   [jet->g fake-photon events: how often
  //                                     does a second jet also fake the lepton?]
  //   f_el = N(DF_ejl) / N(SF_pho)   [e->g fake-photon events: same question]
  //   f_jg = N(DF_jjl) / N(SF_lep)   [prompt-g events: how often does a jet
  //
  // Sequential prediction:
  //   N_pred(VR_A) = N(SF_pho) * f_jl  [jet->g method predicts jet->l fakes]
  //   N_pred(VR_B) = N(SF_pho) * f_el  [e->g method predicts e->g+jet->l fakes]
  //
  // Note: f_jl and f_el are derived from the SAME inclusive sample —
  // this is a self-consistency test, not an extrapolation. It tests whether
  // N(DF) = N(SF) * f, i.e., whether the double-fake rate is consistent
  // with a factorised product of single-fake rates.
  // If pred = truth: fakes factorise and the sequential method is valid.
  // If pred ≠ truth: fakes are correlated and a correction is needed.
  // ============================================================================

  double nSR_DF_total = nSR_DF_ejl + nSR_DF_jjl;

  // Transfer factors
  double f_jl = (nSR_SFpho > 0) ? nSR_DF_jjl / nSR_SFpho : 0.0;  // jet->g+jet->l rate
  double f_el = (nSR_SFpho > 0) ? nSR_DF_ejl / nSR_SFpho : 0.0;  // e->g+jet->l rate
  double f_jg = (nSR_SFlep > 0) ? nSR_DF_jjl / nSR_SFlep : 0.0;  // complementary jet->g rate

  // Sequential predictions
  // VR_A truth = DF_jjl events; prediction from SF_pho * f_jl
  double nVRA_pred = nVRA_SFpho * f_jl;
  // VR_B truth = DF_ejl events; prediction from SF_pho * f_el
  double nVRB_pred = nVRB_SFpho * f_el;

  // Build prediction histograms
  for(int ih = 0; ih < 6; ih++){
    hVRA_pred[ih]->Add(hVRA_SFpho[ih], f_jl);
    hVRB_pred[ih]->Add(hVRB_SFpho[ih], f_el);
  }

  // ============================================================================
  // OUTPUT — Stage 4: Closure comparison
  // ============================================================================
  std::string sep(80,'=');
  std::string sepm(80,'-');

  printf("\n%s\n", sep.c_str());
  printf("  FULL CLOSURE TEST RESULTS  [channel: %s  year: %d%s]\n",
         chainname.c_str(), RunYear,whichVFP.c_str());
  printf("%s\n\n", sep.c_str());

  printf("  SIGNAL REGION (MET>120, MT>100) — TRUTH DECOMPOSITION\n");
  printf("  %s\n", sepm.c_str());
  printf("  %-35s %12s %12s %12s\n", "Category", "Events", "Fraction", "Type");
  printf("  %s\n", sepm.c_str());
  double nSR_total_fake = nSR_DF_total + nSR_SFpho + nSR_SFlep;
  printf("  %-35s %12.1f %11.1f%% %s\n",
         "DF_ejl  (e->g + jet->l)",
         nSR_DF_ejl, 100.*nSR_DF_ejl/std::max(nSR_total_fake,1.), "double fake");
  printf("  %-35s %12.1f %11.1f%% %s\n",
         "DF_jjl  (jet->g + jet->l)",
         nSR_DF_jjl, 100.*nSR_DF_jjl/std::max(nSR_total_fake,1.), "double fake");
  printf("  %-35s %12.1f %11.1f%% %s\n",
         "SF_pho  (fake-g + prompt-l)",
         nSR_SFpho, 100.*nSR_SFpho/std::max(nSR_total_fake,1.), "single fake");
  printf("  %-35s %12.1f %11.1f%% %s\n",
         "SF_lep  (prompt-g + fake-l)",
         nSR_SFlep, 100.*nSR_SFlep/std::max(nSR_total_fake,1.), "single fake");
  printf("  %-35s %12.1f %s\n",
         "Both prompt (sanity)", nSR_Both, "(not background)");
  printf("  %s\n", sepm.c_str());
  printf("  %-35s %12.1f\n", "Total DF", nSR_DF_total);
  printf("  %-35s %12.1f\n", "Total SF", nSR_SFpho + nSR_SFlep);
  printf("  %-35s %12.4f\n", "f_df = DF/SF", (nSR_SFpho+nSR_SFlep>0) ?
         nSR_DF_total/(nSR_SFpho+nSR_SFlep) : 0.0);

  printf("\n  DERIVED TRANSFER FACTORS (from inclusive truth counts)\n");
  printf("  %s\n", sepm.c_str());
  printf("  f_jl = N(DF_jjl) / N(SF_pho_ejl) = %.1f / %.1f = %.4f\n",
         nSR_DF_jjl, nSR_SF_ejl, nSR_DF_jjl / nSR_SF_ejl);
  printf("  f_jl = N(DF_jjl) / N(SF_pho_jjl) = %.1f / %.1f = %.4f\n",
         nSR_DF_jjl, nSR_SF_jjl, nSR_DF_jjl / nSR_SF_jjl);
  
  printf("  f_jl = N(DF_jjl) / N(SF_pho) = %.1f / %.1f = %.4f\n",
         nSR_DF_jjl, nSR_SFpho, f_jl);
  printf("  f_el = N(DF_ejl) / N(SF_pho) = %.1f / %.1f = %.4f\n",
         nSR_DF_ejl, nSR_SFpho, f_el);
  printf("  f_jg = N(DF_jjl) / N(SF_lep) = %.1f / %.1f = %.4f\n",
         nSR_DF_jjl, nSR_SFlep, f_jg);
  printf("  Factorisation check: if fakes are uncorrelated,\n");
  printf("    f_jl should equal f_jg * (lepton fake rate / photon fake rate)\n");

  printf("\n  VR_A jet->g + jet->l enriched (truth: PHO_JETFAKE + LEP_FAKE)\n");
  printf("  %s\n", sepm.c_str());
  double nVRA_DFtrue = nVRA_DF_ejl + nVRA_DF_jjl;
  printf("  N_truth (DF_ejl)     = %8.1f\n", nVRA_DF_ejl);
  printf("  N_truth (DF_jjl)     = %8.1f\n", nVRA_DF_jjl);
  printf("  N_truth (DF total)   = %8.1f\n", nVRA_DFtrue);
  printf("  N_pred  (sequential) = %8.1f  [= N(SF_pho) %.1f * f_jl %.4f]\n",
         nVRA_pred, nVRA_SFpho, f_jl);
  double ratioVRA = (nVRA_DFtrue > 0) ? nVRA_pred / nVRA_DFtrue : 0.0;
  printf("  Closure ratio pred/truth = %.3f  ", ratioVRA);
  if(fabs(ratioVRA - 1.0) < 0.15) printf("[GOOD: within 15%%]\n");
  else                             printf("[WARNING: >15%% non-closure]\n");
  double sysVRA = fabs(1.0 - ratioVRA);
  printf("  Non-closure (|1 - pred/truth|) = %.1f%%  → systematic on DF\n",
         100.*sysVRA);

  printf("\n  VR_B e->g + jet->l enriched (truth: PHO_EFAKE + LEP_FAKE)\n");
  printf("  %s\n", sepm.c_str());
  double nVRB_DFtrue = nVRB_DF_ejl + nVRB_DF_jjl;
  printf("  N_truth (DF_ejl)     = %8.1f\n", nVRB_DF_ejl);
  printf("  N_truth (DF_jjl)     = %8.1f\n", nVRB_DF_jjl);
  printf("  N_truth (DF total)   = %8.1f\n", nVRB_DFtrue);
  printf("  N_pred  (sequential) = %8.1f  [= N(SF_pho) %.1f * f_el %.4f]\n",
         nVRB_pred, nVRB_SFpho, f_el);
  double ratioVRB = (nVRB_DFtrue > 0) ? nVRB_pred / nVRB_DFtrue : 0.0;
  printf("  Closure ratio pred/truth = %.3f  ", ratioVRB);
  if(fabs(ratioVRB - 1.0) < 0.15) printf("[GOOD: within 15%%]\n");
  else                             printf("[WARNING: >15%% non-closure]\n");
  double sysVRB = fabs(1.0 - ratioVRB);
  printf("  Non-closure (|1 - pred/truth|) = %.1f%% systematic on DF\n",
         100.*sysVRB);

  double nDF_SR_total = nSR_DF_ejl + nSR_DF_jjl;
  printf("\n  DECOMPOSITION TEST (all DF events which method absorbs each?)\n");
  printf("  %s\n", sepm.c_str());
  printf("  Total DF events (inclusive)      = %8.1f\n", nDF_SR_total);
  printf("  Absorbed by e->g method ONLY     = %8.1f  (%.1f%%)\n",
         nDF_absorbedByEfakePho,
         100.*nDF_absorbedByEfakePho/std::max(nDF_SR_total,1.));
  printf("  Absorbed by jet->g method ONLY   = %8.1f  (%.1f%%)\n",
         nDF_absorbedByJetfakePho,
         100.*nDF_absorbedByJetfakePho/std::max(nDF_SR_total,1.));
  printf("  Absorbed by jet->l method ONLY   = %8.1f  (%.1f%%)\n",
         nDF_absorbedByJetfakeLep,
         100.*nDF_absorbedByJetfakeLep/std::max(nDF_SR_total,1.));
  printf("  Absorbed by MULTIPLE methods     = %8.1f  (%.1f%%) [double-count bug if >0]\n",
         nDF_absorbedByMultiple,
         100.*nDF_absorbedByMultiple/std::max(nDF_SR_total,1.));
  printf("  Absorbed by NO method (missed)   = %8.1f  (%.1f%%) [missed if >0]\n",
         nDF_absorbedByNone,
         100.*nDF_absorbedByNone/std::max(nDF_SR_total,1.));
  printf("  %s\n", sepm.c_str());
  if(nDF_absorbedByMultiple > 0.1 * nDF_SR_total)
    printf("  *** WARNING: >10%% of DF events absorbed by multiple methods "
           "double-counting in background estimate ***\n");
  if(nDF_absorbedByNone > 0.1 * nDF_SR_total)
    printf("  *** WARNING: >10%% of DF events absorbed by no method "
           "DF contribution missed in background estimate ***\n");
  if(nDF_absorbedByMultiple < 0.1 * nDF_SR_total &&
     nDF_absorbedByNone     < 0.1 * nDF_SR_total)
    printf("  DECOMPOSITION OK: each DF event absorbed by exactly one method.\n");

  double sys_final = std::max(sysVRA, sysVRB);
  printf("\n  SYSTEMATIC UNCERTAINTY RECOMMENDATION\n");
  printf("  %s\n", sepm.c_str());
  printf("  Non-closure in VR_A = %.1f%%\n", 100.*sysVRA);
  printf("  Non-closure in VR_B = %.1f%%\n", 100.*sysVRB);
  printf("  Recommended systematic on DF = max(VRA,VRB) = %.1f%%\n",
         100.*sys_final);
  printf("  Apply as normalization uncertainty on: ");
  if(channelType == 1) printf("jet->g and e->g background estimates (eg channel)\n");
  else           printf("jet->g background estimate (mg channel)\n");
  printf("%s\n\n", sep.c_str());

  // ============================================================================
  // STAGE 4 PLOTS: Prediction vs Truth in VRs as function of key observables
  // ============================================================================
  std::string outDir = "/eos/uscms/store/user/tmishra/elefakepho/Closure/";
  std::string tag = std::string(chainname) + "_" + std::to_string(RunYear);

  // Observable names for axis and file labeling
  const char* obsName[6]  = {"MET","MT","phoEt","lepPt","HT","dR"};
  const char* obsLabel[6] = {"MET (GeV)","M_{T} (GeV)","E_{T}^{#gamma} (GeV)",
                              "p_{T}^{l} (GeV)","HT (GeV)","#DeltaR(#gamma,l)"};

  // For each observable, make a 2-panel canvas: top=shapes, bottom=pred/truth
  for(int iobs = 0; iobs < 6; iobs++){

    {
      TH1D *hTruth = (TH1D*)hVRA_DF_jjl[iobs]->Clone("hVRA_truth_tmp");
      hTruth->Add(hVRA_DF_ejl[iobs]);     // total truth DF in VR_A
      TH1D *hPred  = (TH1D*)hVRA_pred[iobs]->Clone("hVRA_pred_tmp");

      TCanvas *c = new TCanvas(
        Form("c_VRA_%s_%s", obsName[iobs], tag.c_str()), "", 600, 700);
      c->Divide(1,2);

      // Top pad: overlaid shapes
      c->cd(1)->SetPad(0,0.35,1,1);
      hTruth->SetLineColor(kBlue+1);  hTruth->SetLineWidth(2);
      hPred ->SetLineColor(kRed+1);   hPred ->SetLineWidth(2);
      hPred ->SetLineStyle(2);
      double ymax = std::max(hTruth->GetMaximum(), hPred->GetMaximum()) * 1.4;
      hTruth->SetMaximum(ymax);
      hTruth->GetXaxis()->SetTitle(obsLabel[iobs]);
      hTruth->GetYaxis()->SetTitle("Events");
      hTruth->SetTitle(Form("VR_A (jet#rightarrow#gamma + jet#rightarrowl)  %s %d%s",
                            chainname.c_str(), RunYear,whichVFP.c_str()));
      hTruth->Draw("HIST E");
      hPred ->Draw("HIST SAME");
      TLegend *leg = new TLegend(0.55,0.65,0.88,0.88);
      leg->SetBorderSize(0);
      leg->AddEntry(hTruth, "MC Truth (DF)", "l");
      leg->AddEntry(hPred,  "Sequential Pred.", "l");
      leg->Draw();

      // Bottom pad: ratio
      c->cd(2)->SetPad(0,0,1,0.35);
      TH1D *hRatio = (TH1D*)hPred->Clone("hVRA_ratio_tmp");
      hRatio->Divide(hTruth);
      hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2);
      hRatio->SetMarkerStyle(20);   hRatio->SetMarkerSize(0.8);
      hRatio->GetYaxis()->SetTitle("Pred/Truth");
      hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
      hRatio->GetXaxis()->SetTitle(obsLabel[iobs]);
      hRatio->Draw("EP");
      TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                               hRatio->GetXaxis()->GetXmax(), 1.0);
      line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();

      c->SaveAs(Form("%sclosure_v2_VRA_%s_%s.png",
                     outDir.c_str(), obsName[iobs], tag.c_str()));
      delete hTruth; delete hPred; delete hRatio; delete c; delete leg; delete line;
    }

    {
      TH1D *hTruth = (TH1D*)hVRB_DF_ejl[iobs]->Clone("hVRB_truth_tmp");
      hTruth->Add(hVRB_DF_jjl[iobs]);
      TH1D *hPred  = (TH1D*)hVRB_pred[iobs]->Clone("hVRB_pred_tmp");

      TCanvas *c = new TCanvas(
        Form("c_VRB_%s_%s", obsName[iobs], tag.c_str()), "", 600, 700);
      c->Divide(1,2);
      c->cd(1)->SetPad(0,0.35,1,1);
      hTruth->SetLineColor(kBlue+1);  hTruth->SetLineWidth(2);
      hPred ->SetLineColor(kRed+1);   hPred ->SetLineWidth(2);
      hPred ->SetLineStyle(2);
      double ymax = std::max(hTruth->GetMaximum(), hPred->GetMaximum()) * 1.4;
      hTruth->SetMaximum(ymax);
      hTruth->GetXaxis()->SetTitle(obsLabel[iobs]);
      hTruth->GetYaxis()->SetTitle("Events");
      hTruth->SetTitle(Form("VR_B (e#rightarrow#gamma + jet#rightarrowl)  %s %d%s",
                            chainname.c_str(), RunYear,whichVFP.c_str()));
      hTruth->Draw("HIST E");
      hPred ->Draw("HIST SAME");
      TLegend *leg = new TLegend(0.55,0.65,0.88,0.88);
      leg->SetBorderSize(0);
      leg->AddEntry(hTruth, "MC Truth (DF)", "l");
      leg->AddEntry(hPred,  "Sequential Pred.", "l");
      leg->Draw();
      c->cd(2)->SetPad(0,0,1,0.35);
      TH1D *hRatio = (TH1D*)hPred->Clone("hVRB_ratio_tmp");
      hRatio->Divide(hTruth);
      hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2);
      hRatio->SetMarkerStyle(20);   hRatio->SetMarkerSize(0.8);
      hRatio->GetYaxis()->SetTitle("Pred/Truth");
      hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
      hRatio->GetXaxis()->SetTitle(obsLabel[iobs]);
      hRatio->Draw("EP");
      TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                               hRatio->GetXaxis()->GetXmax(), 1.0);
      line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
      c->SaveAs(Form("%sclosure_v2_VRB_%s_%s.png",
                     outDir.c_str(), obsName[iobs], tag.c_str()));
      delete hTruth; delete hPred; delete hRatio; delete c; delete leg; delete line;
    }
  }

  std::string outFile = outDir + "closure_doubleFake_v2_" + tag + ".root";
  TFile *fout = TFile::Open(outFile.c_str(), "RECREATE");
  for(int i=0;i<6;i++){
    hSR_DF_ejl[i]->Write(); hSR_DF_jjl[i]->Write();
    hSR_SFpho[i]->Write();  hSR_SFlep[i]->Write();
    hSR_Both[i]->Write();
    hVRA_DF_ejl[i]->Write(); hVRA_DF_jjl[i]->Write();
    hVRA_SFpho[i]->Write();  hVRA_SFlep[i]->Write();
    hVRA_pred[i]->Write();
    hVRB_DF_ejl[i]->Write(); hVRB_DF_jjl[i]->Write();
    hVRB_SFpho[i]->Write();  hVRB_SFlep[i]->Write();
    hVRB_pred[i]->Write();
  }
  fout->Close();
  printf("  Output written to: %s\n", outFile.c_str());
}
