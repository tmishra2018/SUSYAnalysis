// Converts ggNtuplizer trees into flat Z->ee tag-and-probe trees with a photon probe, readable by egamma-tnp
// (workflow PhotonTagNProbeFromMiniNTuples), for the photon R9 efficiency/SF.
// Tag: tight, trigger-matched electron (as in makeTnPTree.C). Probe: the other electron reconstructed as a photon,
// passing the analysis photon selection (passSignalSelection: calib Et > 35, Fall17v2 loose), no pixel-seed veto.
// g++ `root-config --cflags` ../../lib/libAnaClasses.so makeTnPPhoTree.C -o makeTnPPhoTree.exe `root-config --libs`
// ./makeTnPPhoTree.exe <input.root> <output.root> <RunYear> <isMC 0/1> [maxEvents]
#include "../../include/analysis_commoncode.h"
#include "../../include/analysis_cuts.h"

// Analysis photon R9 cut (EB and EE), e.g. test/VGamma/src/analysis_VGamma.C
#define phoR9Cut 0.5

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
  // Read only what the electron and photon objects need
  es->SetBranchStatus("*", 0);
  for(const char *br : {"ele*", "nEle", "pho*", "nPho", "rho", "run", "event", "lumis", "nGoodVtx", "nVtx"}) es->SetBranchStatus(br, 1);
  if(isMC){
    es->SetBranchStatus("puTrue", 1);
    es->SetBranchStatus("puBX", 1);
  }

  TFile *outputfile = TFile::Open(outName, "RECREATE");
  TDirectory *dir = outputfile->mkdir("tnpPhoIDs");
  dir->cd();
  TTree *tree = new TTree("fitter_tree", "fitter_tree");
  float tag_Ele_pt, tag_Ele_eta, tag_Ele_phi, tag_sc_eta;
  float ph_et, ph_eta, ph_phi, ph_sc_eta, ph_r9;
  int   tag_Ele_q;
  float pair_mass, truePU(0);
  int   run, lumi, nVtx;
  Long64_t event;
  int   passingLoose, passingR9, ph_hasPixelSeed, mcTrue(1);
  tree->Branch("tag_Ele_pt", &tag_Ele_pt);
  tree->Branch("tag_Ele_eta", &tag_Ele_eta);
  tree->Branch("tag_Ele_phi", &tag_Ele_phi);
  tree->Branch("tag_sc_eta", &tag_sc_eta);
  tree->Branch("tag_Ele_q", &tag_Ele_q);
  tree->Branch("ph_et", &ph_et);
  tree->Branch("ph_eta", &ph_eta);
  tree->Branch("ph_phi", &ph_phi);
  tree->Branch("ph_sc_eta", &ph_sc_eta);
  tree->Branch("ph_r9", &ph_r9);
  tree->Branch("ph_hasPixelSeed", &ph_hasPixelSeed);
  tree->Branch("pair_mass", &pair_mass);
  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("event", &event);
  tree->Branch("nVtx", &nVtx);
  tree->Branch("passingLoose", &passingLoose);
  tree->Branch("passingR9", &passingR9);
  // egamma-tnp treats a file as MC if it has an mcTrue branch
  if(isMC){
    tree->Branch("mcTrue", &mcTrue);
    tree->Branch("truePU", &truePU);
  }

  std::vector<recoEle> Ele;
  std::vector<recoPhoton> Photon;
  Long64_t nEvts = es->GetEntries();
  if(maxEvents > 0 && maxEvents < nEvts) nEvts = maxEvents;
  std::cout << "Processing " << nEvts << " events, RunType " << datatype << std::endl;
  Long64_t nPairs(0);

  for(Long64_t ievt(0); ievt < nEvts; ++ievt){
    if(ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
    raw.GetData(es, ievt);
    if(raw.nEle < 1 || raw.nPho < 1 || raw.nGoodVtx < 1) continue;

    Ele.clear();
    Photon.clear();
    for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
    for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}

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

      for(unsigned iprobe(0); iprobe < Photon.size(); iprobe++){
        recoPhoton &probe = Photon[iprobe];
        if(!probe.isEB() && !probe.isEE()) continue;
        if(!probe.passSignalSelection()) continue;
        // the tag electron's own photon object
        if(DeltaR(tag.getEta(), tag.getPhi(), probe.getEta(), probe.getPhi()) < 0.3) continue;

        pair_mass = (tag.getCalibP4() + probe.getCalibP4()).M();
        if(pair_mass < 50 || pair_mass > 130) continue;

        tag_Ele_pt = tag.getCalibPt();
        tag_Ele_eta = tag.getEta();
        tag_Ele_phi = tag.getPhi();
        tag_sc_eta = tag.getEta();
        tag_Ele_q = tag.getCharge();
        ph_et = probe.getCalibEt();
        ph_eta = probe.getEta();
        ph_phi = probe.getPhi();
        ph_sc_eta = probe.getEta();
        ph_r9 = probe.getR9();
        ph_hasPixelSeed = probe.PixelSeed();

        passingLoose = probe.isLoose();
        passingR9 = probe.getR9() > phoR9Cut;

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
