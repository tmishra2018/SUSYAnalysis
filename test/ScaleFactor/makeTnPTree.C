// Converts ggNtuplizer trees into flat tag-and-probe trees readable by egamma-tnp
// (workflow ElectronTagNProbeFromMiniNTuples), using the analysis electron definitions.
// g++ `root-config --cflags` ../../lib/libAnaClasses.so makeTnPTree.C -o makeTnPTree.exe `root-config --libs`
// ./makeTnPTree.exe <input.root> <output.root> <RunYear> <isMC 0/1> [maxEvents]
#include "../../include/analysis_commoncode.h"
#include "../../include/analysis_cuts.h"

// isMiniMedium() without the miniIso cut: denominator for the miniIso efficiency
bool passMiniMediumNoIso(recoEle &ele){
  if(!ele.isEB() && !ele.isEE()) return false;
  const double esc = ele.getCalibP4().E();
  if(esc <= 0) return false;
  const double absDeta = fabs(ele.getdEtaIn());
  const double absDphi = fabs(ele.getdPhiIn());
  const double absEopInv = fabs(ele.getEoverPInv());
  if(ele.isEB()){
    const double hoeMedium = 0.046 + 1.16/esc + 0.0324*ele.getRho()/esc;
    return (ele.getSigma() <= 0.0106 && absDeta <= 0.0032 && absDphi <= 0.0547 &&
            ele.getHoverE() <= hoeMedium && absEopInv <= 0.184 &&
            ele.getMissHits() <= 1 && ele.getConvVeto());
  }
  const double hoeMedium = 0.0275 + 2.52/esc + 0.183*ele.getRho()/esc;
  return (ele.getSigma() <= 0.0387 && absDeta <= 0.00632 && absDphi <= 0.0394 &&
          ele.getHoverE() <= hoeMedium && absEopInv <= 0.0721 &&
          ele.getMissHits() <= 1 && ele.getConvVeto());
}

bool passR9(recoEle &ele){
  if(ele.isEB()) return ele.getR9() > R9EBCut;
  if(ele.isEE()) return ele.getR9() > R9EECut;
  return false;
}

int main(int argc, char **argv){
  if(argc < 5){
    std::cout << "usage: " << argv[0] << " <input.root> <output.root> <RunYear> <isMC 0/1> [maxEvents]" << std::endl;
    return 1;
  }
  const char *inName = argv[1];
  const char *outName = argv[2];
  int RunYear = atoi(argv[3]);
  bool isMC = atoi(argv[4]);
  Long64_t maxEvents = argc > 5 ? atoll(argv[5]) : -1;

  // Tags must fire the single-electron trigger (matched per electron via eleFiredSingleTrgs)
  RunType datatype = MC;
  if(RunYear==2016) datatype = isMC ? MCSingleElectron2016 : SingleElectron2016;
  if(RunYear==2017) datatype = isMC ? MCSingleElectron2017 : SingleElectron2017;
  if(RunYear==2018) datatype = isMC ? MCSingleElectron2018 : SingleElectron2018;
  // Tag pT threshold above the trigger turn-on: Ele27 (2016), Ele35 (2017), Ele32 (2018)
  float tagPtCut = 35;
  if(RunYear==2016) tagPtCut = 30;
  if(RunYear==2017) tagPtCut = 38;

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add(inName);
  rawData raw(es, datatype);
  // Read only what the electron objects need: much faster on the full ggNtuples
  es->SetBranchStatus("*", 0);
  for(const char *br : {"ele*", "nEle", "rho", "run", "event", "lumis", "nGoodVtx", "nVtx"}) es->SetBranchStatus(br, 1);
  if(isMC){
    es->SetBranchStatus("puTrue", 1);
    es->SetBranchStatus("puBX", 1);
  }

  TFile *outputfile = TFile::Open(outName, "RECREATE");
  TDirectory *dir = outputfile->mkdir("tnpEleIDs");
  dir->cd();
  TTree *tree = new TTree("fitter_tree", "fitter_tree");
  float tag_Ele_pt, tag_Ele_eta, tag_Ele_phi, tag_sc_eta;
  float el_pt, el_eta, el_phi, el_sc_eta, el_r9, el_miniIso;
  int   tag_Ele_q, el_q;
  float pair_mass, truePU(0);
  int   run, lumi, nVtx;
  Long64_t event;
  int   passingMiniMediumNoIso, passingMiniMedium, passingMiniIso, passingR9, mcTrue(1);
  tree->Branch("tag_Ele_pt", &tag_Ele_pt);
  tree->Branch("tag_Ele_eta", &tag_Ele_eta);
  tree->Branch("tag_Ele_phi", &tag_Ele_phi);
  tree->Branch("tag_sc_eta", &tag_sc_eta);
  tree->Branch("tag_Ele_q", &tag_Ele_q);
  tree->Branch("el_pt", &el_pt);
  tree->Branch("el_eta", &el_eta);
  tree->Branch("el_phi", &el_phi);
  tree->Branch("el_sc_eta", &el_sc_eta);
  tree->Branch("el_q", &el_q);
  tree->Branch("el_r9", &el_r9);
  tree->Branch("el_miniIso", &el_miniIso);
  tree->Branch("pair_mass", &pair_mass);
  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("event", &event);
  tree->Branch("nVtx", &nVtx);
  tree->Branch("passingMiniMediumNoIso", &passingMiniMediumNoIso);
  tree->Branch("passingMiniMedium", &passingMiniMedium);
  tree->Branch("passingMiniIso", &passingMiniIso);
  tree->Branch("passingR9", &passingR9);
  // egamma-tnp treats a file as MC if it has an mcTrue branch
  if(isMC){
    tree->Branch("mcTrue", &mcTrue);
    tree->Branch("truePU", &truePU);
  }

  std::vector<recoEle> Ele;
  Long64_t nEvts = es->GetEntries();
  if(maxEvents > 0 && maxEvents < nEvts) nEvts = maxEvents;
  std::cout << "Processing " << nEvts << " events, RunType " << datatype << std::endl;
  Long64_t nPairs(0);

  for(Long64_t ievt(0); ievt < nEvts; ++ievt){
    if(ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
    raw.GetData(es, ievt);
    if(raw.nEle < 2 || raw.nGoodVtx < 1) continue;

    Ele.clear();
    for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}

    run = raw.run;
    lumi = raw.lumis;
    event = raw.event;
    nVtx = raw.nVtx;
    if(isMC){
      truePU = 0;
      for(unsigned ipu(0); ipu < raw.puBX->size(); ipu++){
        if((*raw.puBX)[ipu] == 0){truePU = (*raw.puTrue)[ipu]; break;}
      }
    }

    for(unsigned itag(0); itag < Ele.size(); itag++){
      recoEle &tag = Ele[itag];
      if(tag.getCalibPt() < tagPtCut) continue;
      if(!tag.isEB() && !tag.isEE()) continue;
      if(!tag.isTight()) continue;
      if(!tag.passHLTSelection()) continue;

      for(unsigned iprobe(0); iprobe < Ele.size(); iprobe++){
        if(iprobe == itag) continue;
        recoEle &probe = Ele[iprobe];
        if(!probe.isEB() && !probe.isEE()) continue;
        if(probe.getCalibPt() < 10) continue;

        pair_mass = (tag.getCalibP4() + probe.getCalibP4()).M();
        if(pair_mass < 50 || pair_mass > 130) continue;

        tag_Ele_pt = tag.getCalibPt();
        tag_Ele_eta = tag.getEta();
        tag_Ele_phi = tag.getPhi();
        tag_sc_eta = tag.getEta();
        tag_Ele_q = tag.getCharge();
        el_pt = probe.getCalibPt();
        el_eta = probe.getEta();
        el_phi = probe.getPhi();
        el_sc_eta = probe.getEta();
        el_q = probe.getCharge();
        el_r9 = probe.getR9();
        el_miniIso = probe.getMiniIso();

        passingMiniMediumNoIso = passMiniMediumNoIso(probe);
        passingMiniMedium = probe.isMiniMedium();
        passingMiniIso = probe.getMiniIso() <= 0.1;
        passingR9 = passR9(probe);

        tree->Fill();
        nPairs++;
      }
    }
  }

  std::cout << "Wrote " << nPairs << " tag-probe pairs to " << outName << std::endl;
  dir->cd();
  tree->Write();
  outputfile->Close();
  return 0;
}
