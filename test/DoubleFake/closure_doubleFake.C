// ============================================================================
// closure_doubleFake.C
//
// MC STUDY: Quantification of the double-fake background
// (fake lepton AND fake photon simultaneously passing signal selection)
//
// PHYSICS MOTIVATION
// ==================
// The signal selection requires one lepton + one photon. In the analysis,
// each fake object is estimated separately by data-driven methods:
//   - e->gamma    : analysis_eleBkg.C  (proxyTree, hasPho=false)
//   - jet->gamma  : analysis_jetBkg.C  (jetTree,   hasPho=false)
//   - jet->lepton : analysis_qcdBkg.C  (fakeLepTree, hasPho=true, hasLep=false)
//
// None of these estimate the case where BOTH the lepton AND photon are fakes.
// This MC study quantifies this contribution using truth matching on simulation.
//
// TRUTH MATCHING LOGIC (consistent with existing closure scripts and VGBkg/rareBkg)
// =================================================================================
//
// PHOTON classification (mirrors closure_jetfakepho.C and closure_efakepho.C):
//   isTruePho  : mcPID==22 within dR<0.2, mother is prompt (q/g/W/Z/lepton)
//                Same logic as analysis_VGBkg.C and analysis_rareBkg.C
//   isEleFake  : best-match mcPID==11 or 211+pixel seed, or no match to prompt
//                photon — mirrors closure_efakepho.C
//   isJetFake  : no match to any prompt gen-level photon or electron within
//                dR<0.2, i.e. isHad() is true for best match — mirrors
//                closure_jetfakepho.C
//
// LEPTON classification (mirrors closure_jetfakelep.C):
//   isPromptLep : gen-level match within dR<0.1 to lepton from W/Z decay
//                 (mcMomPID==24 or 23)
//   isFakeLep   : no gen-level prompt lepton match — from jet/QCD activity
//
// SAMPLE STRATEGY
// ===============
// The double-fake process requires BOTH objects to be non-prompt.
//
//   PRIMARY (double-fake enriched): QCD HT-binned multijet
//     No prompt leptons or photons by construction. Cross section completely
//     dominates over hadronic ttbar (~100x larger in the relevant HT regime),
//     so QCD alone is sufficient. All HT bins from 300 GeV upward are included
//     to cover the full kinematic range of the selection.
//
//   REFERENCE (single fake, prompt companion):
//     - W+jets     : prompt lepton + jet->gamma or e->gamma fake
//     - Gamma+jets : prompt photon + jet->lepton fake
//     These represent the single-fake backgrounds already estimated by the
//     data-driven methods (analysis_eleBkg.C, analysis_jetBkg.C, analysis_qcdBkg.C).
//
// The double-fake FRACTION = N(both fake) / N(exactly one fake)
// is the key quantity. If this fraction is small (<5%), the existing
// single-fake estimates are sufficient. If large, a correction is needed.
//
// HOW TO USE
// ==========
// Run separately on each MC sample set. The script produces:
//   1. Event counts and fractions printed to screen / log
//   2. Histograms of key kinematic distributions split by object truth category
//   3. A ratio plot: N(doubleFake) / N(singleFake) as function of MET, MT, phoEt
//
// These ratios can then be used as a correction factor or upper limit on the
// double-fake contamination in the signal region.
//
// ============================================================================
#include "TROOT.h"
#include "../../include/analysis_commoncode.h"
static bool DEBUG_PHOCLASS = false;  // set true to print momPID diagnostics
// Photon truth categories — exclusive
enum PhoCategory { PHO_TRUE, PHO_EFAKE, PHO_JETFAKE, PHO_UNKNOWN };

// Lepton truth categories — exclusive
enum LepCategory { LEP_PROMPT, LEP_FAKE, LEP_UNKNOWN };

PhoCategory classifyPhoton(float phoEta, float phoPhi, float phoEt,
                            std::vector<int>   *mcPID,
                            std::vector<float> *mcEta,
                            std::vector<float> *mcPhi,
                            std::vector<float> *mcPt,
                            std::vector<int>   *mcMomPID,
                            int /*pixelSeed — not available in resTree*/) {

  // ── Step 1: TRUE PHOTON — analysis_VGamma.C istruepho ────────────────────
  // Scan for best match with dR<0.3 AND dE<0.5. Require winner at dR<0.2.
  // pid==22 with mother in prompt list → PHO_TRUE.
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
                     momPID == 23 || momPID == 24)){
      return PHO_TRUE;
    }
    if(DEBUG_PHOCLASS && pid == 22)
      printf("  [DEBUG] PHO_TRUE fail: phoEt=%.1f dR=%.3f dE=%.3f momPID=%d\n",
             phoEt, mindRpho, fabs((*mcPt)[phoIndex]-phoEt)/phoEt, momPID);
  }

  // ── Step 2: e->gamma FAKE — closure_efakepho.C isFakePho ─────────────────
  // Separate scan with NO dE cut (dE cut explicitly commented out in source).
  // Condition A: best match (dR<0.3, no dE cut) within dR<0.1 is gen e±.
  // Condition B: any gen e± within dR<0.02 (collinear conversion/Dalitz).
  double mindR_eFake = 0.3;
  double minEledR    = 3.0;
  unsigned matchIndex = 0;
  for(unsigned iMC = 0; iMC < mcPID->size(); iMC++){
    double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
    if(dR < mindR_eFake){ mindR_eFake = dR; matchIndex = iMC; }
    if(fabs((*mcPID)[iMC]) == 11 && dR < minEledR) minEledR = dR;
  }
  if(mindR_eFake < 0.1){
    if((*mcPID)[matchIndex] == 11 || (*mcPID)[matchIndex] == -11)
      return PHO_EFAKE;
  }
  if(minEledR < 0.02) return PHO_EFAKE;

  // ── Step 3: jet->gamma FAKE — closure_jetfakepho.C isFakePho ─────────────
  // Exclusion logic: default true (jet fake). Set false only if best match
  // (dR<0.3, dE<0.5) within dR<0.1 is a prompt e± or γ.
  // Correctly handles π⁰→γγ: pid==22 but hadronic momPID → stays jet fake.
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
  if(isFakePho){
    if(DEBUG_PHOCLASS && mindR_jFake < 0.2)
      printf("  [DEBUG] PHO_JETFAKE: phoEt=%.1f pid=%d momPID=%d dR=%.3f\n",
             phoEt, (*mcPID)[jMatchIndex], (*mcMomPID)[jMatchIndex], mindR_jFake);
    return PHO_JETFAKE;
  }

  return PHO_UNKNOWN;
}


// ============================================================================
// LEPTON TRUTH CLASSIFICATION
// ============================================================================
// Copied verbatim from closure_jetfakelep.C isFakeLep logic.
// Prompt: gen e/μ within dR<0.1, dE<0.3, mother = W(24), Z(23), or τ(15).
// Fake: everything else — no match, wrong mother, b→ℓν, QCD jet→ℓ.
//
// NOTE: dE<0.3 here matches the framework source exactly. This is tighter
// than the photon dE<0.5 because leptons are better measured than photons
// and the tighter cut reduces fake-prompt misclassification from b→ℓν.
// ============================================================================
LepCategory classifyLepton(float lepEta, float lepPhi, float lepPt,
                            std::vector<int>   *mcPID,
                            std::vector<float> *mcEta,
                            std::vector<float> *mcPhi,
                            std::vector<float> *mcPt,
                            std::vector<int>   *mcMomPID) {

  double mindR   = 0.1;
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
void closure_doubleFake(){
  SetRunConfig();
  setTDRStyle();
  gSystem->Load("../../lib/libAnaClasses.so");
  gROOT->SetBatch(kTRUE);

  int channelType = ichannel; // eg=1, mg=2
  
  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
  // ==========================================================================
  // MC sample chains
  // Two categories of samples for the study:
  //
  //   Category A — "doubleFake enriched": QCD HT-binned multijet
  //     No prompt leptons or photons. All HT bins from 300 GeV upward.
  //     QCD cross section (~100x hadronic ttbar) completely dominates the
  //     double-fake yield, so hadronic ttbar is not needed.
  //
  //   Category B — "singleFake reference": W+jets and Gamma+jets
  //     W+jets  : prompt lepton from W + fake photon (jet->gamma or e->gamma)
  //     GJets   : prompt photon + fake lepton (jet->lepton)
  //     These are the single-fake processes already estimated data-driven.
  //
  // The ratio N_A(doubleFake) / N_B(singleFake) gives the double-fake fraction.
  // ==========================================================================

  std::ostringstream chainname;
  if(channelType == 1) chainname << "egTree";
  else                 chainname << "mgTree";

  // Category A: double-fake enriched — QCD HT-binned, all bins from HT300 up
  // HT300to500 and HT500to700 dominate numerically; higher bins cover the
  // high-MET tail relevant to the signal region.
  // File naming convention follows the existing framework resTree structure.
  TChain *qcdTree = new TChain(chainname.str().c_str());

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
  TChain *gjetsTree = new TChain(chainname.str().c_str());
  gjetsTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_GJets_DoubleEM_%d%s.root",RunYear,whichVFP.c_str()));

  // W+jets: prompt lepton + fake photon
  TChain *wjetsTree = new TChain(chainname.str().c_str());
  wjetsTree->Add(Form("/eos/uscms/store/user/tmishra/VGammaNEW/resTree_VGamma_WJetsToLNu_%d%s.root",RunYear,whichVFP.c_str()));

  // ==========================================================================
  // Branch setup — identical to analysis_VGBkg.C and analysis_rareBkg.C
  // ==========================================================================
  struct TreeVars {
    float crosssection, ntotalevent, PUweight;
    int   mcType;
    float phoEt, phoEta, phoPhi;
    float lepPt, lepEta, lepPhi;
    float sigMT, sigMET, sigMETPhi, dPhiLepMET;
    int   nVertex;
    float dRPhoLep, HT, nJet;
    int   nBJet;
    std::vector<int>   *mcPID    = 0;
    std::vector<float> *mcEta   = 0;
    std::vector<float> *mcPhi   = 0;
    std::vector<float> *mcPt    = 0;
    std::vector<int>   *mcMomPID = 0;
    std::vector<int>   *mcStatus = 0;
  };

  auto connectBranches = [](TChain *tree, TreeVars &v){
    tree->SetBranchAddress("crosssection", &v.crosssection);
    tree->SetBranchAddress("ntotalevent",  &v.ntotalevent);
    tree->SetBranchAddress("PUweight",     &v.PUweight);
    tree->SetBranchAddress("mcType",       &v.mcType);
    tree->SetBranchAddress("phoEt",        &v.phoEt);
    tree->SetBranchAddress("phoEta",       &v.phoEta);
    tree->SetBranchAddress("phoPhi",       &v.phoPhi);
    tree->SetBranchAddress("lepPt",        &v.lepPt);
    tree->SetBranchAddress("lepEta",       &v.lepEta);
    tree->SetBranchAddress("lepPhi",       &v.lepPhi);
    tree->SetBranchAddress("sigMT",        &v.sigMT);
    tree->SetBranchAddress("sigMET",       &v.sigMET);
    tree->SetBranchAddress("sigMETPhi",    &v.sigMETPhi);
    tree->SetBranchAddress("dPhiLepMET",   &v.dPhiLepMET);
    tree->SetBranchAddress("nVertex",      &v.nVertex);
    tree->SetBranchAddress("dRPhoLep",     &v.dRPhoLep);
    tree->SetBranchAddress("HT",           &v.HT);
    tree->SetBranchAddress("nJet",         &v.nJet);
    tree->SetBranchAddress("nBJet",        &v.nBJet);
    tree->SetBranchAddress("mcPID",        &v.mcPID);
    tree->SetBranchAddress("mcEta",        &v.mcEta);
    tree->SetBranchAddress("mcPhi",        &v.mcPhi);
    tree->SetBranchAddress("mcPt",         &v.mcPt);
    tree->SetBranchAddress("mcMomPID",     &v.mcMomPID);
  };

  TreeVars vQCD, vGJets, vWJets;
  connectBranches(qcdTree,   vQCD);
  connectBranches(gjetsTree, vGJets);
  connectBranches(wjetsTree, vWJets);

  // ==========================================================================
  // Output histograms
  // For each observable, four histograms:
  //   _doubleFake  : both lepton and photon are fake (our target)
  //   _fakePhoOnly : prompt lepton, fake photon  (single-fake reference)
  //   _fakeLepOnly : prompt photon, fake lepton  (single-fake reference)
  //   _bothPrompt  : both prompt (sanity check / signal contamination)
  // ==========================================================================
  TH1D *h_MET_doubleFake   = new TH1D("h_MET_doubleFake",  "MET; MET (GeV);Events",  nBkgMETBins, bkgMETBins);
  TH1D *h_MET_fakePhoOnly  = new TH1D("h_MET_fakePhoOnly", "MET; MET (GeV);Events",  nBkgMETBins, bkgMETBins);
  TH1D *h_MET_fakeLepOnly  = new TH1D("h_MET_fakeLepOnly", "MET; MET (GeV);Events",  nBkgMETBins, bkgMETBins);
  TH1D *h_MET_bothPrompt   = new TH1D("h_MET_bothPrompt",  "MET; MET (GeV);Events",  nBkgMETBins, bkgMETBins);

  TH1D *h_MT_doubleFake    = new TH1D("h_MT_doubleFake",   "M_{T}; M_{T} (GeV);",    nBkgMtBins,  bkgMtBins);
  TH1D *h_MT_fakePhoOnly   = new TH1D("h_MT_fakePhoOnly",  "M_{T}; M_{T} (GeV);",    nBkgMtBins,  bkgMtBins);
  TH1D *h_MT_fakeLepOnly   = new TH1D("h_MT_fakeLepOnly",  "M_{T}; M_{T} (GeV);",    nBkgMtBins,  bkgMtBins);
  TH1D *h_MT_bothPrompt    = new TH1D("h_MT_bothPrompt",   "M_{T}; M_{T} (GeV);",    nBkgMtBins,  bkgMtBins);

  TH1D *h_phoEt_doubleFake = new TH1D("h_phoEt_doubleFake","#gamma E_{T}; E_{T} (GeV);",nBkgEtBins, bkgEtBins);
  TH1D *h_phoEt_fakePhoOnly= new TH1D("h_phoEt_fakePhoOnly","#gamma E_{T}; E_{T} (GeV);",nBkgEtBins,bkgEtBins);
  TH1D *h_phoEt_fakeLepOnly= new TH1D("h_phoEt_fakeLepOnly","#gamma E_{T}; E_{T} (GeV);",nBkgEtBins,bkgEtBins);
  TH1D *h_phoEt_bothPrompt = new TH1D("h_phoEt_bothPrompt", "#gamma E_{T}; E_{T} (GeV);",nBkgEtBins,bkgEtBins);

  TH1D *h_lepPt_doubleFake = new TH1D("h_lepPt_doubleFake","Lepton p_{T}; p_{T} (GeV);",nBkgPtBins, bkgPtBins);
  TH1D *h_lepPt_fakePhoOnly= new TH1D("h_lepPt_fakePhoOnly","Lepton p_{T}; p_{T} (GeV);",nBkgPtBins,bkgPtBins);
  TH1D *h_lepPt_fakeLepOnly= new TH1D("h_lepPt_fakeLepOnly","Lepton p_{T}; p_{T} (GeV);",nBkgPtBins,bkgPtBins);
  TH1D *h_lepPt_bothPrompt = new TH1D("h_lepPt_bothPrompt", "Lepton p_{T}; p_{T} (GeV);",nBkgPtBins,bkgPtBins);

  TH1D *h_HT_doubleFake    = new TH1D("h_HT_doubleFake",  "HT; HT (GeV);",  nBkgHTBins, bkgHTBins);
  TH1D *h_HT_fakePhoOnly   = new TH1D("h_HT_fakePhoOnly", "HT; HT (GeV);",  nBkgHTBins, bkgHTBins);
  TH1D *h_HT_fakeLepOnly   = new TH1D("h_HT_fakeLepOnly", "HT; HT (GeV);",  nBkgHTBins, bkgHTBins);
  TH1D *h_HT_bothPrompt    = new TH1D("h_HT_bothPrompt",  "HT; HT (GeV);",  nBkgHTBins, bkgHTBins);

  // Photon sub-category in the double-fake: was it e->gamma or jet->gamma?
  TH1D *h_MET_df_efake  = new TH1D("h_MET_df_efake", "doubleFake (e#rightarrow#gamma); MET (GeV);",nBkgMETBins,bkgMETBins);
  TH1D *h_MET_df_jfake  = new TH1D("h_MET_df_jfake", "doubleFake (jet#rightarrow#gamma); MET (GeV);",nBkgMETBins,bkgMETBins);

  // Global truth-category counters (summed over all samples)
  double nDoubleFake(0), nFakePhoOnly(0), nFakeLepOnly(0), nBothPrompt(0);
  double nDF_efake(0),   nDF_jfake(0);

  // Per-sample truth-category counters (all four categories per sample)
  // QCD: double-fake enriched — expected to dominate double fake
  double nQCD_df(0),  nQCD_df_efake(0), nQCD_df_jfake(0);
  double nQCD_sfPho(0), nQCD_sfLep(0), nQCD_both(0);
  // GJets: prompt photon + fake lepton reference
  double nGJ_df(0),   nGJ_df_efake(0),  nGJ_df_jfake(0);
  double nGJ_sfPho(0),  nGJ_sfLep(0),  nGJ_both(0);
  // WJets: prompt lepton + fake photon reference
  double nWJ_df(0),   nWJ_df_efake(0),  nWJ_df_jfake(0);
  double nWJ_sfPho(0),  nWJ_sfLep(0),  nWJ_both(0);

  // Legacy aliases (used by processTree lambda signature)
  double nDF_QCD(0);
  double nSF_Gjets(0), nSF_Wjets(0);

  // ==========================================================================
  // Lambda: process one tree (called for each sample chain)
  // sampleTag: 0=QCD, 1=GJets, 2=WJets
  // ==========================================================================
 auto processTree = [&](TChain *tree, TreeVars &v, int sampleTag,
                         double &nDF_local, double &nSF_local) {

    for(unsigned ievt = 0; ievt < tree->GetEntries(); ++ievt){
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
             if(RunYear == 2016 and preVFP == 1)             	XS_weight = lumi_2016preVFP_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2016 and preVFP == 0)        	XS_weight = lumi_2016postVFP_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2017)                        	XS_weight = lumi_2017_DoubleEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2018)                        	XS_weight = lumi_2018_DoubleEG*1000*v.crosssection/v.ntotalevent;}

      else if(channelType == 2){
             if(RunYear == 2016 and preVFP == 1)             	XS_weight = lumi_2016preVFP_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2016 and preVFP == 0)        	XS_weight = lumi_2016postVFP_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2017)                        	XS_weight = lumi_2017_MuonEG*1000*v.crosssection/v.ntotalevent;
             else if(RunYear == 2018)                        	XS_weight = lumi_2018_MuonEG*1000*v.crosssection/v.ntotalevent;}
      float weight    = v.PUweight * XS_weight;

      // ---- Build per-object filtered gen-particle views ----------------
      // The eg-channel ntuple stores ALL gen particles >=5 GeV (no spatial
      // pre-filter), while the mg-channel and version-A ntuples only store
      // particles within dR<0.3 of the signal objects. To guarantee identical
      // classifier behaviour across channels and ntuple versions, we build
      // filtered views here — exactly replicating the ntuple-filling logic:
      //   photon view : gen particles within dR<0.3 of reco photon
      //   lepton view : gen particles within dR<0.3 of reco lepton
      // (For mg the stored vector already satisfies this, so filtering is
      // a no-op that does not change the result.)
      std::vector<int>   pho_mcPID,    lep_mcPID;
      std::vector<float> pho_mcEta,    lep_mcEta;
      std::vector<float> pho_mcPhi,    lep_mcPhi;
      std::vector<float> pho_mcPt,     lep_mcPt;
      std::vector<int>   pho_mcMomPID, lep_mcMomPID;
      for(unsigned iMC = 0; iMC < v.mcPID->size(); iMC++){
        float dR_pho = DeltaR(v.phoEta, v.phoPhi,
                              (*v.mcEta)[iMC], (*v.mcPhi)[iMC]);
        float dR_lep = DeltaR(v.lepEta, v.lepPhi,
                              (*v.mcEta)[iMC], (*v.mcPhi)[iMC]);
        if(dR_pho < 0.3f){
          pho_mcPID.push_back((*v.mcPID)[iMC]);
          pho_mcEta.push_back((*v.mcEta)[iMC]);
          pho_mcPhi.push_back((*v.mcPhi)[iMC]);
          pho_mcPt .push_back((*v.mcPt) [iMC]);
          pho_mcMomPID.push_back((*v.mcMomPID)[iMC]);
        }
        if(dR_lep < 0.3f){
          lep_mcPID.push_back((*v.mcPID)[iMC]);
          lep_mcEta.push_back((*v.mcEta)[iMC]);
          lep_mcPhi.push_back((*v.mcPhi)[iMC]);
          lep_mcPt .push_back((*v.mcPt) [iMC]);
          lep_mcMomPID.push_back((*v.mcMomPID)[iMC]);
        }
      }

      PhoCategory phoCat = classifyPhoton(v.phoEta, v.phoPhi, v.phoEt,
                                          v.mcPID, v.mcEta, v.mcPhi, v.mcPt,
                                          v.mcMomPID, 0);

      // ---- Lepton truth classification ----
      LepCategory lepCat = classifyLepton(v.lepEta, v.lepPhi, v.lepPt,
                                          v.mcPID, v.mcEta, v.mcPhi, v.mcPt,
                                          v.mcMomPID);

      // ---- Photon truth classification ----
      //PhoCategory phoCat = classifyPhoton(v.phoEta, v.phoPhi, v.phoEt,
      //                                    &pho_mcPID, &pho_mcEta,
      //                                    &pho_mcPhi, &pho_mcPt,
      //                                    &pho_mcMomPID, 0);

      // ---- Lepton truth classification ----
      //LepCategory lepCat = classifyLepton(v.lepEta, v.lepPhi, v.lepPt,
      //                                    &lep_mcPID, &lep_mcEta,
      //                                    &lep_mcPhi, &lep_mcPt,
      //                                    &lep_mcMomPID);


      // ---- Fill histograms by truth category ----
      bool phoIsFake   = (phoCat == PHO_EFAKE || phoCat == PHO_JETFAKE);
      bool phoIsTrue   = (phoCat == PHO_TRUE);
      bool lepIsFake   = (lepCat == LEP_FAKE);
      bool lepIsPrompt = (lepCat == LEP_PROMPT);

      if(phoIsFake && lepIsFake){
        // DOUBLE FAKE — the quantity we want to measure
        h_MET_doubleFake ->Fill(v.sigMET, weight);
        h_MT_doubleFake  ->Fill(v.sigMT,  weight);
        h_phoEt_doubleFake->Fill(v.phoEt, weight);
        h_lepPt_doubleFake->Fill(v.lepPt, weight);
        h_HT_doubleFake  ->Fill(v.HT,     weight);
        nDoubleFake      += weight;
        nDF_local        += weight;
        // Sub-classify photon fake type in double-fake events
        if(phoCat == PHO_EFAKE){
          h_MET_df_efake->Fill(v.sigMET, weight);
          nDF_efake += weight;
          if(sampleTag==0){ nQCD_df += weight; nQCD_df_efake += weight; }
          if(sampleTag==1){ nGJ_df  += weight; nGJ_df_efake  += weight; }
          if(sampleTag==2){ nWJ_df  += weight; nWJ_df_efake  += weight; }
        } else {
          h_MET_df_jfake->Fill(v.sigMET, weight);
          nDF_jfake += weight;
          if(sampleTag==0){ nQCD_df += weight; nQCD_df_jfake += weight; }
          if(sampleTag==1){ nGJ_df  += weight; nGJ_df_jfake  += weight; }
          if(sampleTag==2){ nWJ_df  += weight; nWJ_df_jfake  += weight; }
        }
      }
      else if(phoIsFake && lepIsPrompt){
        // SINGLE FAKE: fake photon, prompt lepton
        h_MET_fakePhoOnly ->Fill(v.sigMET, weight);
        h_MT_fakePhoOnly  ->Fill(v.sigMT,  weight);
        h_phoEt_fakePhoOnly->Fill(v.phoEt, weight);
        h_lepPt_fakePhoOnly->Fill(v.lepPt, weight);
        h_HT_fakePhoOnly  ->Fill(v.HT,     weight);
        nFakePhoOnly      += weight;
        nSF_local         += weight;
        if(sampleTag==0) nQCD_sfPho += weight;
        if(sampleTag==1) nGJ_sfPho  += weight;
        if(sampleTag==2) nWJ_sfPho  += weight;
      }
      else if(phoIsTrue && lepIsFake){
        // SINGLE FAKE: prompt photon, fake lepton
        h_MET_fakeLepOnly ->Fill(v.sigMET, weight);
        h_MT_fakeLepOnly  ->Fill(v.sigMT,  weight);
        h_phoEt_fakeLepOnly->Fill(v.phoEt, weight);
        h_lepPt_fakeLepOnly->Fill(v.lepPt, weight);
        h_HT_fakeLepOnly  ->Fill(v.HT,     weight);
        nFakeLepOnly      += weight;
        nSF_local         += weight;
        if(sampleTag==0) nQCD_sfLep += weight;
        if(sampleTag==1) nGJ_sfLep  += weight;
        if(sampleTag==2) nWJ_sfLep  += weight;
      }
      else if(phoIsTrue && lepIsPrompt){
        // BOTH PROMPT — expected to be tiny in QCD/had-TT samples
        // (sanity check; if large → sample contamination from semi-leptonic ttbar)
        h_MET_bothPrompt ->Fill(v.sigMET, weight);
        h_MT_bothPrompt  ->Fill(v.sigMT,  weight);
        h_phoEt_bothPrompt->Fill(v.phoEt, weight);
        h_lepPt_bothPrompt->Fill(v.lepPt, weight);
        h_HT_bothPrompt  ->Fill(v.HT,     weight);
        nBothPrompt      += weight;
        if(sampleTag==0) nQCD_both += weight;
        if(sampleTag==1) nGJ_both  += weight;
        if(sampleTag==2) nWJ_both  += weight;
      }
      // PHO_UNKNOWN / LEP_UNKNOWN: not counted (tiny residual)
    }
  };
  // ==========================================================================
  // Process each sample
  // ==========================================================================
  std::cout << "\n=== Processing QCD HT-binned (double-fake enriched) ===" << std::endl;
  processTree(qcdTree,   vQCD,   0, nDF_QCD,   nSF_Gjets);  // nSF placeholder, not used for QCD

  std::cout << "\n=== Processing GJets (single-fake lepton reference) ===" << std::endl;
  double dummy1(0);
  processTree(gjetsTree, vGJets, 1, dummy1, nSF_Gjets);

  std::cout << "\n=== Processing W+jets (single-fake photon reference) ===" << std::endl;
  double dummy2(0);
  processTree(wjetsTree, vWJets, 2, dummy2, nSF_Wjets);

  // ==========================================================================
  // Compute double-fake fraction
  // Key quantity: how often does a reconstructed signal-like event (1 lep + 1 pho
  // passing all selections) arise from BOTH objects being fakes, relative to
  // the more common case where exactly one object is a fake?
  // ==========================================================================
  double nSingleFake = nFakePhoOnly + nFakeLepOnly;
  double dfFraction  = (nSingleFake > 0) ? nDoubleFake / nSingleFake : 0;
  double dfFracErr   = (nSingleFake > 0 && nDoubleFake > 0)
                       ? dfFraction * sqrt(1.0/nDoubleFake + 1.0/nSingleFake) : 0;

	
  // ==========================================================================
  // Print results — structured per-sample truth table
  // ==========================================================================
std::string ch = (channelType==1) ? "eg" : "mg";
  std::string sep80(80, '=');
  std::string sep80d(80, '-');

  std::cout << "\\n" << sep80 << std::endl;
  std::cout << "  DOUBLE FAKE MC STUDY RESULTS  [channel: " << ch << "]" << std::endl;
  std::cout << sep80 << std::endl;

  // Column widths
  std::cout << std::fixed << std::setprecision(1);
  std::cout << std::left;

  // Header
  std::cout << std::setw(32) << "  Truth Category"
            << std::setw(12) << "QCD"
            << std::setw(12) << "GJets"
            << std::setw(12) << "W+jets"
            << std::setw(12) << "TOTAL"
            << std::endl;
  std::cout << "  " << sep80d << std::endl;

  // Note on sample role
  std::cout << "  [role]"
            << std::setw(26) << ""
            << std::setw(12) << "df-enriched"
            << std::setw(12) << "sf-lep ref"
            << std::setw(12) << "sf-pho ref"
            << std::endl;
  std::cout << "  " << sep80d << std::endl;

  // Double fake row
  std::cout << std::setw(32) << "  Double fake (fake g + fake l)"
            << std::setw(12) << nQCD_df
            << std::setw(12) << nGJ_df
            << std::setw(12) << nWJ_df
            << std::setw(12) << nDoubleFake
            << std::endl;
  std::cout << std::setw(32) << "    of which e->g fake"
            << std::setw(12) << nQCD_df_efake
            << std::setw(12) << nGJ_df_efake
            << std::setw(12) << nWJ_df_efake
            << std::setw(12) << nDF_efake
            << std::endl;
  std::cout << std::setw(32) << "    of which jet->g fake"
            << std::setw(12) << nQCD_df_jfake
            << std::setw(12) << nGJ_df_jfake
            << std::setw(12) << nWJ_df_jfake
            << std::setw(12) << nDF_jfake
            << std::endl;
  std::cout << "  " << sep80d << std::endl;

  // Single fake photon row
  std::cout << std::setw(32) << "  Single fake (fake g, prompt l)"
            << std::setw(12) << nQCD_sfPho
            << std::setw(12) << nGJ_sfPho
            << std::setw(12) << nWJ_sfPho
            << std::setw(12) << nFakePhoOnly
            << std::endl;

  // Single fake lepton row
  std::cout << std::setw(32) << "  Single fake (prompt g, fake l)"
            << std::setw(12) << nQCD_sfLep
            << std::setw(12) << nGJ_sfLep
            << std::setw(12) << nWJ_sfLep
            << std::setw(12) << nFakeLepOnly
            << std::endl;
  std::cout << "  " << sep80d << std::endl;

  // Total single fake
  double nSF_QCD = nQCD_sfPho + nQCD_sfLep;
  double nSF_GJ  = nGJ_sfPho  + nGJ_sfLep;
  double nSF_WJ  = nWJ_sfPho  + nWJ_sfLep;
  std::cout << std::setw(32) << "  Total single fake"
            << std::setw(12) << nSF_QCD
            << std::setw(12) << nSF_GJ
            << std::setw(12) << nSF_WJ
            << std::setw(12) << nSingleFake
            << std::endl;
  std::cout << "  " << sep80d << std::endl;

  // Both prompt (sanity)
  std::cout << std::setw(32) << "  Both prompt (sanity check)"
            << std::setw(12) << nQCD_both
            << std::setw(12) << nGJ_both
            << std::setw(12) << nWJ_both
            << std::setw(12) << nBothPrompt
            << std::endl;
  std::cout << "  " << sep80 << std::endl;

  // Fraction
  std::cout << std::setprecision(4);
  std::cout << "  Double-fake fraction  f_df = N(df) / N(sf)" << std::endl;
  std::cout << "    N(df) = " << std::setprecision(1) << nDoubleFake
            << "  [QCD: " << nQCD_df << " | GJets: " << nGJ_df
            << " | W+jets: " << nWJ_df << "]" << std::endl;
  std::cout << "    N(sf) = " << nSingleFake
            << "  [QCD: " << nSF_QCD << " | GJets: " << nSF_GJ
            << " | W+jets: " << nSF_WJ << "]" << std::endl;
  std::cout << std::setprecision(4);
  std::cout << "    f_df  = " << dfFraction
            << " +/- " << dfFracErr
            << "  (" << dfFraction*100 << " %)" << std::endl;

  // Physics note: which samples contribute to each category
  // Data-driven contamination warning
  std::cout << std::endl;
  std::cout << "  Sample purity check:" << std::endl;
  // Double fake purity: QCD should dominate
  double dfTot = nDoubleFake > 0 ? nDoubleFake : 1;
  std::cout << std::setprecision(1);
  std::cout << "    Double fake:  QCD=" << 100*nQCD_df/dfTot << "%"
            << "  GJets=" << 100*nGJ_df/dfTot << "%"
            << "  W+jets=" << 100*nWJ_df/dfTot << "%"
            << "  (expect QCD~100%)" << std::endl;
  // Single fake (fake g) purity: W+jets should dominate
  double sfgTot = nFakePhoOnly > 0 ? nFakePhoOnly : 1;
  std::cout << "    Single fake (fake g):  QCD=" << 100*nQCD_sfPho/sfgTot << "%"
            << "  GJets=" << 100*nGJ_sfPho/sfgTot << "%"
            << "  W+jets=" << 100*nWJ_sfPho/sfgTot << "%"
            << "  (expect W+jets~100%)" << std::endl;
  // Single fake (fake l) purity: GJets should dominate
  double sflTot = nFakeLepOnly > 0 ? nFakeLepOnly : 1;
  std::cout << "    Single fake (fake l):  QCD=" << 100*nQCD_sfLep/sflTot << "%"
            << "  GJets=" << 100*nGJ_sfLep/sflTot << "%"
            << "  W+jets=" << 100*nWJ_sfLep/sflTot << "%"
            << "  (expect GJets~100%)" << std::endl;
  // Contamination warnings
  bool contamDF  = (nQCD_df/dfTot < 0.7 && nDoubleFake > 5);
  bool contamSFG = (nWJ_sfPho/sfgTot < 0.7 && nFakePhoOnly > 5);
  bool contamSFL = (nGJ_sfLep/sflTot < 0.7 && nFakeLepOnly > 5);
  if(contamDF || contamSFG || contamSFL){
    std::cout << std::endl;
    std::cout << "  *** WARNING: cross-sample contamination detected ***" << std::endl;
    if(contamDF)
      std::cout << "      Double fake not dominated by QCD -- check photon truth matching" << std::endl;
    if(contamSFG)
      std::cout << "      Single fake (fake g) not dominated by W+jets -- check photon truth matching" << std::endl;
    if(contamSFL)
      std::cout << "      Single fake (fake l) not dominated by GJets -- check lepton truth matching" << std::endl;
  } else {
    std::cout << "    Sample composition as expected." << std::endl;
  }
  std::cout << "  " << sep80 << std::endl;

  if(dfFraction < 0.05)
    std::cout << "  CONCLUSION: f_df = " << dfFraction*100
              << "% < 5% -- negligible, no correction needed." << std::endl;
  else if(dfFraction < 0.15)
    std::cout << "  CONCLUSION: f_df = " << dfFraction*100
              << "% in 5-15% range -- assign as systematic uncertainty." << std::endl;
  else
    std::cout << "  CONCLUSION: f_df = " << dfFraction*100
              << "% > 15% -- significant, explicit correction needed." << std::endl;
  std::cout << "  " << sep80 << std::endl;

 
  // ==========================================================================
  // Ratio histograms: N(doubleFake) / N(singleFake) as function of kinematics
  // These show whether the double-fake fraction is uniform or peaks in
  // specific kinematic regions (e.g. low MET, high HT) which would indicate
  // whether the signal region (high MET, high MT) is more or less contaminated
  // ==========================================================================
  TH1D *h_ratio_MET = (TH1D*)h_MET_doubleFake->Clone("h_ratio_MET_dfOverSF");
  TH1D *h_singleFake_MET = (TH1D*)h_MET_fakePhoOnly->Clone("h_singleFakeTot_MET");
  h_singleFake_MET->Add(h_MET_fakeLepOnly);
  h_ratio_MET->Divide(h_singleFake_MET);
  h_ratio_MET->GetYaxis()->SetTitle("N(double fake) / N(single fake)");

  TH1D *h_ratio_MT = (TH1D*)h_MT_doubleFake->Clone("h_ratio_MT_dfOverSF");
  TH1D *h_singleFake_MT = (TH1D*)h_MT_fakePhoOnly->Clone("h_singleFakeTot_MT");
  h_singleFake_MT->Add(h_MT_fakeLepOnly);
  h_ratio_MT->Divide(h_singleFake_MT);
  h_ratio_MT->GetYaxis()->SetTitle("N(double fake) / N(single fake)");

  TH1D *h_ratio_phoEt = (TH1D*)h_phoEt_doubleFake->Clone("h_ratio_phoEt_dfOverSF");
  TH1D *h_singleFake_phoEt = (TH1D*)h_phoEt_fakePhoOnly->Clone("h_singleFakeTot_phoEt");
  h_singleFake_phoEt->Add(h_phoEt_fakeLepOnly);
  h_ratio_phoEt->Divide(h_singleFake_phoEt);
  h_ratio_phoEt->GetYaxis()->SetTitle("N(double fake) / N(single fake)");

  // ==========================================================================
  // Canvas: kinematic distributions split by truth category
  // Allows visual inspection of whether double-fake events look like
  // the signal region or are suppressed at high MET/MT
  // ==========================================================================
  TCanvas *cMET = new TCanvas("cMET","MET by truth category",800,600);
  cMET->SetLogy();
  h_MET_fakePhoOnly ->SetLineColor(kBlue);   h_MET_fakePhoOnly ->SetLineWidth(2);
  h_MET_fakeLepOnly ->SetLineColor(kGreen+2);h_MET_fakeLepOnly ->SetLineWidth(2);
  h_MET_doubleFake  ->SetLineColor(kRed);    h_MET_doubleFake  ->SetLineWidth(2);
  h_MET_df_efake    ->SetLineColor(kOrange); h_MET_df_efake    ->SetLineWidth(2); h_MET_df_efake->SetLineStyle(2);
  h_MET_df_jfake    ->SetLineColor(kMagenta);h_MET_df_jfake    ->SetLineWidth(2); h_MET_df_jfake->SetLineStyle(2);
  h_MET_bothPrompt  ->SetLineColor(kBlack);  h_MET_bothPrompt  ->SetLineWidth(2);
  double maxMET = std::max({h_MET_fakePhoOnly->GetMaximum(), h_MET_fakeLepOnly->GetMaximum(),
                             h_MET_doubleFake->GetMaximum(), h_MET_bothPrompt->GetMaximum()});
  h_MET_fakePhoOnly->SetMaximum(maxMET*10);
  h_MET_fakePhoOnly->SetMinimum(0.01);
  h_MET_fakePhoOnly->Draw("hist");
  h_MET_fakeLepOnly->Draw("hist same");
  h_MET_doubleFake ->Draw("hist same");
  h_MET_df_efake   ->Draw("hist same");
  h_MET_df_jfake   ->Draw("hist same");
  h_MET_bothPrompt ->Draw("hist same");
  TLegend *legMET = new TLegend(0.55,0.65,0.88,0.90);
  legMET->SetBorderSize(0); legMET->SetFillStyle(0);
  legMET->AddEntry(h_MET_fakePhoOnly,  "fake #gamma + prompt l",  "l");
  legMET->AddEntry(h_MET_fakeLepOnly,  "prompt #gamma + fake l",  "l");
  legMET->AddEntry(h_MET_doubleFake,   "fake #gamma + fake l",    "l");
  legMET->AddEntry(h_MET_df_efake,     "  e#rightarrow#gamma fake","l");
  legMET->AddEntry(h_MET_df_jfake,     "  jet#rightarrow#gamma fake","l");
  legMET->AddEntry(h_MET_bothPrompt,   "both prompt (sanity)",     "l");
  legMET->Draw();
  std::ostringstream canname;
  canname.str(""); canname << "/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_doubleFake_MET_";
  if(channelType==1) canname << "eg"; else canname << "mg";
  canname << "_" << RunYear << whichVFP.c_str() << ".png";
  cMET->SaveAs(canname.str().c_str());

  TCanvas *cRatio = new TCanvas("cRatio","Double-fake / Single-fake ratio",800,600);
  cRatio->Divide(1,3);
  cRatio->cd(1); h_ratio_MET  ->Draw("EP"); h_ratio_MET  ->SetMarkerStyle(20);
  cRatio->cd(2); h_ratio_MT   ->Draw("EP"); h_ratio_MT   ->SetMarkerStyle(20);
  cRatio->cd(3); h_ratio_phoEt->Draw("EP"); h_ratio_phoEt->SetMarkerStyle(20);
  canname.str(""); canname << "/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_doubleFake_ratio_";
  if(channelType==1) canname << "eg"; else canname << "mg";
  canname << "_" << RunYear << whichVFP.c_str() << ".png";
  cRatio->SaveAs(canname.str().c_str());

  // ==========================================================================
  // Write output
  // ==========================================================================
  std::ostringstream outputname;
  outputname << "/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_doubleFake_";
  if(channelType==1) outputname << "eg"; 
  else outputname << "mg";
  outputname << "_" << RunYear <<whichVFP.c_str() <<".root";

  TFile *outputfile = TFile::Open(outputname.str().c_str(), "RECREATE");
  outputfile->cd();
  h_MET_doubleFake ->Write(); h_MET_fakePhoOnly->Write();
  h_MET_fakeLepOnly->Write(); h_MET_bothPrompt ->Write();
  h_MET_df_efake   ->Write(); h_MET_df_jfake   ->Write();
  h_MT_doubleFake  ->Write(); h_MT_fakePhoOnly ->Write();
  h_MT_fakeLepOnly ->Write(); h_MT_bothPrompt  ->Write();
  h_phoEt_doubleFake->Write();h_phoEt_fakePhoOnly->Write();
  h_phoEt_fakeLepOnly->Write();h_phoEt_bothPrompt->Write();
  h_lepPt_doubleFake->Write();h_lepPt_fakePhoOnly->Write();
  h_lepPt_fakeLepOnly->Write();h_lepPt_bothPrompt->Write();
  h_HT_doubleFake  ->Write(); h_HT_fakePhoOnly ->Write();
  h_HT_fakeLepOnly ->Write(); h_HT_bothPrompt  ->Write();
  h_ratio_MET      ->Write(); h_ratio_MT       ->Write();
  h_ratio_phoEt    ->Write();
  outputfile->Write();
  outputfile->Close();

  std::cout << "\nOutput written to: " << outputname.str() << std::endl;
}
