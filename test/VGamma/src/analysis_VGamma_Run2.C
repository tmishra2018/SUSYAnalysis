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
 #include "TProfile2D.h"
 #include "TFileCollection.h"
 
 #include "../../../include/analysis_rawData.h"
 #include "../../../include/analysis_photon.h"
 #include "../../../include/analysis_muon.h"
 #include "../../../include/analysis_jet.h"
 #include "../../../include/analysis_ele.h"
 #include "../../../include/analysis_mcData.h"
 #include "../../../include/analysis_tools.h"

 // PU helper: get nTrueInt for BX==0
 float getNTrueInt(std::vector<float>* puTrue, std::vector<int>* puBX) {
    if (!puTrue || !puBX) return -1.f;
    for (size_t i = 0; i < puBX->size(); ++i)
        if ((*puBX)[i] == 0) return (*puTrue)[i];
    return -1.f;
 }
 
 // Safe lookup from ratio histogram
 float getPUweight(TH1D* h, float nTrueInt) {
    if (!h || nTrueInt < 0) return 1.f;
    int bin = h->FindBin(nTrueInt);
    bin = std::max(1, std::min(bin, h->GetNbinsX()));
    float w = h->GetBinContent(bin);
    return (w > 0) ? w : 1.f;
 }

 bool isHardLepton(int momID){
 		if((fabs(momID) >= 0 && fabs(momID) <= 6) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
 		else return false;
 }
 		
 bool isHardPhoton(int momID){
 		if((fabs(momID) >= 0 && fabs(momID) <= 16) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
 		else return false;
 }
 
 bool nisrMatch(float jetEta, float jetPhi, std::vector<mcData>& genParticles){

     bool matched=false;
     for(std::vector<mcData>::iterator itMC = genParticles.begin(); itMC!= genParticles.end(); itMC++){
      if(matched)break;
 			bool isPromptFinal(false);
      int momid = abs(itMC->getMomPID());
      if(abs(itMC->getPID())<=5 && (momid==6 || momid==23 || momid==24 || momid==25))isPromptFinal = true;
 			else if(abs(itMC->getPID()) == 11 && (momid<=6 || momid==23 || momid==24 || momid== 22 || momid== 15))isPromptFinal = true;
 			else if(abs(itMC->getPID()) == 13 && (momid<=6 || momid==23 || momid==24 || momid== 22 || momid== 15))isPromptFinal = true;
 			else if(abs(itMC->getPID()) == 15 && (momid<=6 || momid==23 || momid==24 || momid== 22))isPromptFinal = true;
 			else if(abs(itMC->getPID()) == 22 && (momid == 11 || momid == 13 || momid == 15 || momid==23 || momid==24))isPromptFinal = true;
 			if(!isPromptFinal)continue;
 			if(DeltaR(jetEta, jetPhi, itMC->getEta(), itMC->getPhi()) < 0.3){
 				//std::cout << "match " << itMC->getPID() << " " << itMC->getMomPID() << " status:" << itMC->getStatus() << std::endl; 
      	matched = true;
 			}
     } // Loop over MC particles
		return matched;
 } 
 
// WG NLO K-factors as a function of photon pT (GeV)
// Binning: 0-190: 1.40, 190-250: 1.37, 250-400: 1.29, 400-700: 1.26, 700-inf: 1.15
// Source: CMS-SMP-16-004, CMS-SMP-22-009
float getWGKfactor(float phoPt){
  if     (phoPt <  190.) return 1.40;
  else if(phoPt <  250.) return 1.37;
  else if(phoPt <  400.) return 1.29;
  else if(phoPt <  700.) return 1.26;
  else                   return 1.15;
}

bool apply_HEMveto=false;

 // Modes (split running on condor, see condor/submit_skim.py):
 //   ""    : whole sample in one job (local), output /uscms_data/d3/tmishra/resTree_VGamma_<sample>_<era>.root
 //   "norm": events [first, last) -> normFile with the signed MC PU histogram and h_count (events, negative genWeight)
 //   "range": events [first, last), normalisation (ntotalevent, MC PU) from normFile = hadd of all "norm" parts;
 //           output and log go to the current directory as resTree_VGamma_<sample>_<era>_<first>.{root,log}
 // VG_XROOTD=1 reads the input via root://cmseos.fnal.gov, VG_PUDIR sets the data PU directory (default ../../PUfiles).
 void analysis_VGamma_Run2(int RunYear, bool preVFP, const char *Sample, std::string mode = "",
                           long long first = 0, long long last = -1, const char *normFile = ""){//main 
   
   std::string whichVFP;
   if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
   if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
   if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
   const bool splitRun = (mode == "norm" || mode == "range");
   const std::string chunkTag = splitRun ? Form("_%lld", first) : "";
 	
   ofstream logfile;
   if(splitRun) logfile.open(Form("resTree_VGamma_%s_%d%s%s%s.log",Sample,RunYear,whichVFP.c_str(),chunkTag.c_str(),mode=="norm"?"_norm":""));
   else         logfile.open(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_%s_%d%s.log",Sample,RunYear,whichVFP.c_str()));

   logfile << "analysis_VGamma_Run2()" << std::endl;
 
   RunType datatype(MC);
   TChain* es = new TChain("ggNtuplizer/EventTree");
   char* inputfile = new char[300];
 
   if (strstr(Sample, "DYJetsToLL") or strstr(Sample, "GJets_DoubleEM") or strstr(Sample, "TTJets") or strstr(Sample, "ZGToLLG") or strstr(Sample, "QCDMuEnriched") or strstr(Sample, "QCD_DoubleEM"))
         sprintf(inputfile,"/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
   else
         sprintf(inputfile,"/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root",Sample,Sample,RunYear,whichVFP.c_str());
   if(getenv("VG_XROOTD")) es->Add(Form("root://cmseos.fnal.gov/%s", inputfile + strlen("/eos/uscms")));
   else es->Add(inputfile);
   if(RunYear==2018) apply_HEMveto=true;
 
   //if (strstr(Sample, "QCD_HT"))
   //	 sprintf(inputfile, "/eos/uscms/store/user/tmishra/QCD_HT_%d%s/%s_%d%s.root", RunYear, whichVFP.c_str(), Sample, RunYear, whichVFP.c_str());
   logfile << "VG";
   logfile << "muon MiniIso" << std::endl;
 
   const unsigned nEvts = es->GetEntries();
   logfile << "Total event: " << nEvts << std::endl;
   if(last < 0 || last > nEvts) last = nEvts;
   if(!splitRun) first = 0;
   logfile << "Event range: [" << first << ", " << last << ")" << std::endl;

    // Build MC PU histogram from nTrueInt (BX==0)
    // NLO samples (ZGToLLG, DYJetsToLL, TTGJets, WWG, WZG, ...) have negative generator weights: each event
    // counts with the sign of genWeight, normalised to N(+) - N(-) (as ISRReweight/analysis_ISRMC.C).
    // LO samples have genWeight = 1, so nothing changes for them.
    TH1D *h_mc_pu = new TH1D("h_mc_pu","MC nTrueInt",200,0,200);
    h_mc_pu->Sumw2();
    long long nNeg(0);
    if(mode == "range"){
        TFile *fn = TFile::Open(normFile);
        TH1D *hc = (TH1D*)fn->Get("h_count");
        if(!hc || (long long)hc->GetBinContent(1) != (long long)nEvts){
            std::cout << "normalisation file " << normFile << " does not cover the " << nEvts << " events of the sample" << std::endl; throw;
        }
        nNeg = (long long)hc->GetBinContent(2);
        h_mc_pu->Add((TH1D*)fn->Get("h_mc_pu"));
        fn->Close();
    }
    else{
        // Temporary local branches for the pre-scan
        std::vector<float> *pt_prescan = nullptr;
        std::vector<int>   *bx_prescan = nullptr;
        float genWeight_prescan(1);
        es->SetBranchStatus("*", 0);         // disable all for speed
        es->SetBranchStatus("puTrue", 1);
        es->SetBranchStatus("puBX",   1);
        es->SetBranchStatus("genWeight", 1);
        es->SetBranchAddress("puTrue", &pt_prescan);
        es->SetBranchAddress("puBX",   &bx_prescan);
        es->SetBranchAddress("genWeight", &genWeight_prescan);
        std::cout << "Pre-scan: building MC PU histogram..." << std::endl;
        for(long long iev=first; iev<last; ++iev){
            es->GetEntry(iev);
            const float sign = genWeight_prescan < 0 ? -1.f : 1.f;
            if(sign < 0) nNeg++;
            float nti = getNTrueInt(pt_prescan, bx_prescan);
            if(nti>=0) h_mc_pu->Fill(nti, sign);
        }
        es->SetBranchStatus("*", 1);  // re-enable all branches
        std::cout << "MC PU histogram built (" << h_mc_pu->GetEntries() << " entries)." << std::endl;
    }
    logfile << "negative genWeight events: " << nNeg << " of " << nEvts << std::endl;
    if(mode == "norm"){
        TFile *fn = new TFile(normFile, "RECREATE");
        h_mc_pu->Write();
        TH1D *hc = new TH1D("h_count", "bin 1: events, bin 2: negative genWeight", 2, 0, 2);
        hc->SetBinContent(1, last - first); hc->SetBinContent(2, nNeg); hc->Write();
        fn->Close(); logfile.close();
        return;
    }
    h_mc_pu->Scale(1.0/h_mc_pu->Integral());
 
    // Load data PU files and build ratio histograms
    const char *puDir = getenv("VG_PUDIR") ? getenv("VG_PUDIR") : "../../PUfiles";
    TFile *f_puNom  = TFile::Open(Form("%s/dataPU_%d%s.root",puDir,RunYear,whichVFP.c_str()));
    TFile *f_puUp   = TFile::Open(Form("%s/dataPU_%d%s_PUUp.root",puDir,RunYear,whichVFP.c_str()));
    TFile *f_puDown = TFile::Open(Form("%s/dataPU_%d%s_PUDown.root",puDir,RunYear,whichVFP.c_str()));

    TH1D *h_data_nom  = (TH1D*)((TH1D*)f_puNom ->Get("pileup"))->Clone("dn");
    TH1D *h_data_up   = (TH1D*)((TH1D*)f_puUp  ->Get("pileup"))->Clone("du");
    TH1D *h_data_down = (TH1D*)((TH1D*)f_puDown->Get("pileup"))->Clone("dd");
    
    h_data_nom ->Scale(1.0/h_data_nom ->Integral());
    h_data_up  ->Scale(1.0/h_data_up  ->Integral());
    h_data_down->Scale(1.0/h_data_down->Integral());
 
    TH1D *h_puRatio_nom  = (TH1D*)h_data_nom ->Clone("h_puRatio_nom");
    TH1D *h_puRatio_up   = (TH1D*)h_data_up  ->Clone("h_puRatio_up");
    TH1D *h_puRatio_down = (TH1D*)h_data_down->Clone("h_puRatio_down");
    h_puRatio_nom ->Divide(h_mc_pu);
    h_puRatio_up  ->Divide(h_mc_pu);
    h_puRatio_down->Divide(h_mc_pu);


   TFile* outputfile;
   if(splitRun) outputfile = new TFile(Form("resTree_VGamma_%s_%d%s%s.root",Sample,RunYear,whichVFP.c_str(),chunkTag.c_str()),"RECREATE");
   else         outputfile = new TFile(Form("/uscms_data/d3/tmishra/resTree_VGamma_%s_%d%s.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
   //outputfile = new TFile(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_%s_%d%s.root",Sample,RunYear,whichVFP.c_str()),"RECREATE");
   outputfile->cd();

   int mcType;
  if(strstr(inputfile, "WGToLNuG") != NULL){
                std::cout << "WGToLNuG sample !" << std::endl;
                mcType = MCType::WGJetInclusive;
  }
  else if(strstr(inputfile, "WGJet40") != NULL){
                std::cout << "WGJet40 sample !" << std::endl;
                mcType = MCType::WGJet40;
  }
  else if(strstr(inputfile, "WGJet130") != NULL){
                std::cout << "WGJet130 sample !" << std::endl;
                mcType = MCType::WGJet130;
  }
  else if(strstr(inputfile, "ZGToLLG") != NULL){
                std::cout << "ZGInclusive sample !" << std::endl;
                mcType = MCType::ZGInclusive;
  }
  else if(strstr(inputfile, "DYJetsToLL") != NULL){
                std::cout << "DYJetsToLL sample !" << std::endl;
                mcType = MCType::DYLL50;
  }
  else  if(strstr(inputfile, "TTGJets") != NULL){
                std::cout << "TTGJets sample !" << std::endl;
                mcType = MCType::TTG;
  }
  else if(strstr(inputfile, "WWG") != NULL){
                std::cout << "WWG sample !" << std::endl;
                mcType = MCType::WWG;
  }
  else if(strstr(inputfile, "WZG") != NULL){
                std::cout << "WZG sample !" << std::endl;
                mcType = MCType::WZG;
  }
  else if(strstr(inputfile, "TTJets") != NULL){
                std::cout << "TTJets sample !" << std::endl;
                mcType = MCType::TT;
  }
  else if(strstr(inputfile, "WW") != NULL){
                std::cout << "WW sample !" << std::endl;
                mcType = MCType::WW;
  }
  else if(strstr(inputfile, "WZ") != NULL){
                std::cout << "WZ sample !" << std::endl;
                mcType = MCType::WZ;
  }
  else if(strstr(inputfile, "ZZ") != NULL){
                std::cout << "ZZ sample !" << std::endl;
                mcType = MCType::ZZ;
  }
  else if(strstr(inputfile, "QCD_DoubleEM") != NULL){
                std::cout << "QCDEM40 sample !" << std::endl;
                mcType = MCType::QCDEM40;
  }
  else if(strstr(inputfile, "QCDMuEnriched") != NULL){
                std::cout << "QCDMU20 sample !" << std::endl;
                mcType = MCType::QCDMU20;
  }
  else if(strstr(inputfile, "WJetsToLNu"))  mcType = MCType::W;
  else if(strstr(inputfile, "GJet"))  mcType = MCType::GJet;
  else if(strstr(inputfile, "QCD_HT50to100"))  mcType = MCType::QCDHT50to100;
  else if(strstr(inputfile, "QCD_HT100to200"))  mcType = MCType::QCDHT100to200;
  else if(strstr(inputfile, "QCD_HT200to300"))  mcType = MCType::QCDHT200to300;
  else if(strstr(inputfile, "QCD_HT300to500"))  mcType = MCType::QCDHT300to500;
  else if(strstr(inputfile, "QCD_HT500to700"))  mcType = MCType::QCDHT500to700;
  else if(strstr(inputfile, "QCD_HT700to1000")) mcType = MCType::QCDHT700to1000;
  else if(strstr(inputfile, "QCD_HT1000to1500"))mcType = MCType::QCDHT1000to1500;
  else if(strstr(inputfile, "QCD_HT1500to2000"))mcType = MCType::QCDHT1500to2000;
  else if(strstr(inputfile, "QCD_HT2000toInf")) mcType = MCType::QCDHT2000toInf;
  
  // Determine if this is a WG sample (K-factor applies)
  bool isWGSample = (mcType == MCType::WGJetInclusive ||
                     mcType == MCType::WGJet40         ||
                     mcType == MCType::WGJet130);


   if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
   logfile << "mcType" << mcType << std::endl;
   float crosssection = MC_XS[mcType];
   cout<< "crosssection: " << crosssection << endl;
   // signed event count; crosssection carries the per-event genWeight sign, so every reader's
   // lumi*crosssection/ntotalevent is the signed NLO weight without changes to the readers
   float ntotalevent = double(nEvts) - 2.0*double(nNeg);
   float genWeight(1);
   float genWeightSign(1);
	float PUweight(1);
	float PUweightUp(1);
        float PUweightDown(1);
	float llmass(0);
	int   nBJet(0);
	float ISRWeight(1);
	float pdfWeight(0);
	std::vector<float> pdfSystWeight;
	std::vector<float> ScaleSystWeight; 
	float psWeightFsrUp(1), psWeightFsrDown(1);
    	float psWeightIsrUp(1), psWeightIsrDown(1);
	// muR/muF: ratio w.r.t. nominal (1.0 = no variation). See below for per-type logic.
	// pdfWeightRMS: Hessian quadrature sum sqrt(Σ(w_i/w_0-1)^2) for first 100 eigenvectors.
	// These replace the old pdfSystWeight and ScaleSystWeight vectors.
	float muRweightUp(1.f), muRweightDown(1.f);
	float muFweightUp(1.f), muFweightDown(1.f);
	float pdfWeightRMS(0.f);

 //************ Signal Tree **********************//
   bool passHEMVeto_eg(true); // true = passes HEM veto (or N/A for this year), false = fails
   bool passHEMVeto_mg(true);
   TTree *egtree = new TTree("egTree","egTree");
   float eg_phoEt(0);
   float eg_phoEta(0);
   float eg_phoPhi(0);
   float eg_lepPt(0);
   float eg_lepEta(0);
   float eg_lepPhi(0);
	float trailPt(0);
	float trailEta(0);
	float trailPhi(0);
   float eg_sigMT(0);
   float eg_sigMET(0);
   float eg_sigMETPhi(0);
   float eg_dPhiLepMET(0);
   int   eg_nVertex(0);
   float eg_dRPhoLep(0);
   float eg_HT(0);
   float eg_nJet(0);
	int   eg_nISRJet(0);
	float eg_ISRJetPt(0);
	float eg_ISRJetPtGen(0);
   float eg_invmass(0);
	float eg_bosonPt(0);
	float eg_sigMETJESup(0);
	float eg_sigMETJESdo(0);
	float eg_sigMETJERup(0);
	float eg_sigMETJERdo(0);
	float eg_sigMTJESup(0);
	float eg_sigMTJESdo(0);
	float eg_sigMTJERup(0);
	float eg_sigMTJERdo(0);
	float eg_HTJESup(0);
	float eg_HTJESdo(0);
	float eg_dPhiLepMETJESup(0);
	float eg_dPhiLepMETJESdo(0);
	float eg_dPhiLepMETJERup(0);
	float eg_dPhiLepMETJERdo(0);
   std::vector<int> eg_mcPID;
   std::vector<float> eg_mcEta;
   std::vector<float> eg_mcPhi;
   std::vector<float> eg_mcPt;
   std::vector<int> eg_mcMomPID;
   std::vector<int> eg_mcGMomPID;
   std::vector<int> eg_mcStatus;
 
	egtree->Branch("crosssection",&crosssection);
	egtree->Branch("ntotalevent", &ntotalevent);
	egtree->Branch("ISRWeight", &ISRWeight);
	egtree->Branch("pdfWeight", &pdfWeight);
	egtree->Branch("pdfSystWeight", &pdfSystWeight);
	egtree->Branch("ScaleSystWeight", &ScaleSystWeight);
   egtree->Branch("phoEt",     &eg_phoEt);
   egtree->Branch("phoEta",    &eg_phoEta);
   egtree->Branch("phoPhi",    &eg_phoPhi);
   egtree->Branch("lepPt",     &eg_lepPt);
   egtree->Branch("lepEta",    &eg_lepEta);
   egtree->Branch("lepPhi",    &eg_lepPhi);
	egtree->Branch("trailPt",   &trailPt);
	egtree->Branch("trailEta",  &trailEta);
	egtree->Branch("trailPhi",  &trailPhi);
	egtree->Branch("PUweight",  &PUweight);
	egtree->Branch("PUweightUp",   &PUweightUp);
        egtree->Branch("PUweightDown", &PUweightDown);
	egtree->Branch("psWeightFsrUp",  &psWeightFsrUp);
    	egtree->Branch("psWeightFsrDown",&psWeightFsrDown);
    	egtree->Branch("psWeightIsrUp",  &psWeightIsrUp);
    	egtree->Branch("psWeightIsrDown",&psWeightIsrDown);
	// QCD scale and PDF weight branches (scalar, for bin-wise lnN)
	egtree->Branch("muRweightUp",   &muRweightUp);
	egtree->Branch("muRweightDown", &muRweightDown);
	egtree->Branch("muFweightUp",   &muFweightUp);
	egtree->Branch("muFweightDown", &muFweightDown);
	egtree->Branch("pdfWeightRMS",  &pdfWeightRMS);
   egtree->Branch("sigMT",     &eg_sigMT);
   egtree->Branch("sigMET",    &eg_sigMET);
   egtree->Branch("sigMETPhi", &eg_sigMETPhi);
   egtree->Branch("dPhiLepMET",&eg_dPhiLepMET);
   egtree->Branch("nVertex",   &eg_nVertex);
   egtree->Branch("dRPhoLep",  &eg_dRPhoLep);
   egtree->Branch("HT",        &eg_HT);
   egtree->Branch("nJet",      &eg_nJet);
   egtree->Branch("nBJet",     &nBJet);
   egtree->Branch("nISRJet",   &eg_nISRJet);
	egtree->Branch("ISRJetPt",  &eg_ISRJetPt);
	egtree->Branch("ISRJetPtGen", &eg_ISRJetPtGen);
	egtree->Branch("genWeightSign", &genWeightSign);
   egtree->Branch("invmass",   &eg_invmass);
	egtree->Branch("bosonPt",   &eg_bosonPt);
	egtree->Branch("sigMETJESup",     &eg_sigMETJESup);
	egtree->Branch("sigMETJESdo",     &eg_sigMETJESdo);
	egtree->Branch("sigMETJERup",     &eg_sigMETJERup);
	egtree->Branch("sigMETJERdo",     &eg_sigMETJERdo);
	egtree->Branch("sigMTJESup",      &eg_sigMTJESup);
	egtree->Branch("sigMTJESdo",      &eg_sigMTJESdo);
	egtree->Branch("sigMTJERup",      &eg_sigMTJERup);
	egtree->Branch("sigMTJERdo",      &eg_sigMTJERdo);
	egtree->Branch("dPhiLepMETJESup", &eg_dPhiLepMETJESup);
	egtree->Branch("dPhiLepMETJESdo", &eg_dPhiLepMETJESdo);
	egtree->Branch("dPhiLepMETJERup", &eg_dPhiLepMETJERup);
	egtree->Branch("dPhiLepMETJERdo", &eg_dPhiLepMETJERdo);
	egtree->Branch("HTJESup",     &eg_HTJESup);
	egtree->Branch("HTJESdo",     &eg_HTJESdo);
   egtree->Branch("mcType",   &mcType);
   egtree->Branch("mcPID",    &eg_mcPID);
   egtree->Branch("mcEta",    &eg_mcEta);
   egtree->Branch("mcPhi",    &eg_mcPhi);
   egtree->Branch("mcPt",     &eg_mcPt);
   egtree->Branch("mcMomPID", &eg_mcMomPID);
   egtree->Branch("mcGMomPID", &eg_mcGMomPID);
   egtree->Branch("mcStatus", &eg_mcStatus);
   egtree->Branch("llmass",   &llmass);
   egtree->Branch("passHEMVeto", &passHEMVeto_eg);
 
   TTree *mgtree = new TTree("mgTree","mgTree");
   float mg_phoEt(0);
   float mg_phoEta(0);
   float mg_phoPhi(0);
   float mg_lepPt(0);
   float mg_lepEta(0);
   float mg_lepPhi(0);
   float mg_sigMT(0);
   float mg_sigMET(0);
   float mg_sigMETPhi(0);
   float mg_dPhiLepMET(0);
	float mg_threeMass(0);
   int   mg_nVertex(0);
   float mg_dRPhoLep(0);
   float mg_HT(0);
   float mg_nJet(0);
   float mg_nISRJet(0);
	float mg_ISRJetPt(0);
	float mg_ISRJetPtGen(0);
	float mg_bosonPt(0);
	float mg_sigMETJESup(0);
	float mg_sigMETJESdo(0);
	float mg_sigMETJERup(0);
	float mg_sigMETJERdo(0);
	float mg_sigMTJESup(0);
	float mg_sigMTJESdo(0);
	float mg_sigMTJERup(0);
	float mg_sigMTJERdo(0);
	float mg_HTJESup(0);
	float mg_HTJESdo(0);
	float mg_dPhiLepMETJESup(0);
	float mg_dPhiLepMETJESdo(0);
	float mg_dPhiLepMETJERup(0);
	float mg_dPhiLepMETJERdo(0);
   std::vector<int> mg_mcPID;
   std::vector<float> mg_mcEta;
   std::vector<float> mg_mcPhi;
   std::vector<float> mg_mcPt;
   std::vector<int> mg_mcMomPID;
   std::vector<int> mg_mcGMomPID;
	std::vector<int> mg_mcStatus;

	mgtree->Branch("crosssection",&crosssection);
	mgtree->Branch("ntotalevent", &ntotalevent);
	mgtree->Branch("ISRWeight", &ISRWeight);
	mgtree->Branch("pdfWeight", &pdfWeight);
	mgtree->Branch("pdfSystWeight", &pdfSystWeight);
	mgtree->Branch("ScaleSystWeight", &ScaleSystWeight);
   mgtree->Branch("phoEt",     &mg_phoEt);
   mgtree->Branch("phoEta",    &mg_phoEta);
   mgtree->Branch("phoPhi",    &mg_phoPhi);
   mgtree->Branch("lepPt",     &mg_lepPt);
   mgtree->Branch("lepEta",    &mg_lepEta);
   mgtree->Branch("lepPhi",    &mg_lepPhi);
	mgtree->Branch("PUweight",      &PUweight);
	mgtree->Branch("PUweightUp",     &PUweightUp);
        mgtree->Branch("PUweightDown",   &PUweightDown);
	mgtree->Branch("psWeightFsrUp",  &psWeightFsrUp);
    	mgtree->Branch("psWeightFsrDown",&psWeightFsrDown);
    	mgtree->Branch("psWeightIsrUp",  &psWeightIsrUp);
    	mgtree->Branch("psWeightIsrDown",&psWeightIsrDown);
	// QCD scale and PDF weight branches (scalar, for bin-wise lnN)
	mgtree->Branch("muRweightUp",   &muRweightUp);
	mgtree->Branch("muRweightDown", &muRweightDown);
	mgtree->Branch("muFweightUp",   &muFweightUp);
	mgtree->Branch("muFweightDown", &muFweightDown);
	mgtree->Branch("pdfWeightRMS",  &pdfWeightRMS);
   mgtree->Branch("sigMT",     &mg_sigMT);
   mgtree->Branch("sigMET",    &mg_sigMET);
   mgtree->Branch("sigMETPhi", &mg_sigMETPhi);
   mgtree->Branch("dPhiLepMET",&mg_dPhiLepMET);
	mgtree->Branch("threeMass", &mg_threeMass);
   mgtree->Branch("nVertex",   &mg_nVertex);
   mgtree->Branch("dRPhoLep",  &mg_dRPhoLep);
   mgtree->Branch("HT",        &mg_HT);
   mgtree->Branch("nJet",      &mg_nJet);
   mgtree->Branch("nBJet",     &nBJet);
   mgtree->Branch("nISRJet",   &mg_nISRJet);
	mgtree->Branch("ISRJetPt",  &mg_ISRJetPt);
	mgtree->Branch("ISRJetPtGen", &mg_ISRJetPtGen);
	mgtree->Branch("genWeightSign", &genWeightSign);
	mgtree->Branch("bosonPt",   &mg_bosonPt);
	mgtree->Branch("sigMETJESup",     &mg_sigMETJESup);
	mgtree->Branch("sigMETJESdo",     &mg_sigMETJESdo);
	mgtree->Branch("sigMETJERup",     &mg_sigMETJERup);
	mgtree->Branch("sigMETJERdo",     &mg_sigMETJERdo);
	mgtree->Branch("sigMTJESup",      &mg_sigMTJESup);
	mgtree->Branch("sigMTJESdo",      &mg_sigMTJESdo);
	mgtree->Branch("sigMTJERup",      &mg_sigMTJERup);
	mgtree->Branch("sigMTJERdo",      &mg_sigMTJERdo);
	mgtree->Branch("dPhiLepMETJESup", &mg_dPhiLepMETJESup);
	mgtree->Branch("dPhiLepMETJESdo", &mg_dPhiLepMETJESdo);
	mgtree->Branch("dPhiLepMETJERup", &mg_dPhiLepMETJERup);
	mgtree->Branch("dPhiLepMETJERdo", &mg_dPhiLepMETJERdo);
	mgtree->Branch("HTJESup",     &mg_HTJESup);
	mgtree->Branch("HTJESdo",     &mg_HTJESdo);
   mgtree->Branch("mcType",   &mcType);
   mgtree->Branch("mcPID",    &mg_mcPID);
   mgtree->Branch("mcEta",    &mg_mcEta);
   mgtree->Branch("mcPhi",    &mg_mcPhi);
   mgtree->Branch("mcPt",     &mg_mcPt);
   mgtree->Branch("mcMomPID", &mg_mcMomPID);
   mgtree->Branch("mcGMomPID", &mg_mcGMomPID);
   mgtree->Branch("mcStatus", &mg_mcStatus);
   mgtree->Branch("llmass",   &llmass);
   mgtree->Branch("passHEMVeto", &passHEMVeto_mg);
 
//*********** histo list **********************//
   rawData raw(es, datatype);
   es->SetBranchAddress("genWeight", &genWeight);
   std::vector<mcData>  MCData;
   std::vector<recoPhoton> Photon;
   std::vector<recoMuon>   Muon;
   std::vector<recoEle>   Ele;
   std::vector<recoJet>   JetCollection;
   float MET(0);
   float METPhi(0);
   float MET_T1JERUp(0);
   float MET_T1JERDo(0);
   float MET_T1JESUp(0);
   float MET_T1JESDo(0);	
   float	METPhi_T1JESUp(0);
   float	METPhi_T1JESDo(0);
   float	METPhi_T1UESUp(0);
   float	METPhi_T1UESDo(0);
   int nVtx(0);
   int jetNumber(0);
   int METFilter(0);
   float ISRJetPt(0);
   logfile << "RunType: " << datatype << std::endl;
   std::cout << "Total evetns : " << nEvts << std::endl;
 
   for (long long ievt(first); ievt<last; ++ievt){//loop on entries
   
       if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
 
         raw.GetData(es, ievt);
         genWeightSign = genWeight < 0 ? -1.f : 1.f;
         MCData.clear();
         Photon.clear();
         Muon.clear();
         Ele.clear();
	 JetCollection.clear();

	 passHEMVeto_eg = true;
	 passHEMVeto_mg = true;
	 if(RunYear==2018){
	 	passHEMVeto_eg = passHEMVeto(1, raw, fracHEM_eg);
	 	passHEMVeto_mg = passHEMVeto(1, raw, fracHEM_mg);
	 }

	 if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
         for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
         for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
         for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
	 for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
         
	 MET = raw.pfMET;
         METPhi = raw.pfMETPhi;
	 MET_T1JERUp = raw.pfMET_T1JERUp;
	 MET_T1JERDo = raw.pfMET_T1JERDo;
	 MET_T1JESUp = raw.pfMET_T1JESUp;
	 MET_T1JESDo = raw.pfMET_T1JESDo;
	 METPhi_T1JESUp = raw.pfMETPhi_T1JESUp;
	 METPhi_T1JESDo = raw.pfMETPhi_T1JESDo;
	 METPhi_T1UESUp = raw.pfMETPhi_T1UESUp;
         METPhi_T1UESDo = raw.pfMETPhi_T1UESDo;
         METFilter = raw.metFilters;
         nVtx = raw.nVtx;
 	
	// PU weights from nTrueInt
        // raw.puTrue and raw.puBX must be available in analysis_rawData.h
        float nTrueInt   = getNTrueInt(raw.puTrue, raw.puBX);
        PUweight         = getPUweight(h_puRatio_nom,  nTrueInt);
        PUweightUp       = getPUweight(h_puRatio_up,   nTrueInt);   // NEW
        PUweightDown     = getPUweight(h_puRatio_down, nTrueInt);   // NEW	 
 	
	psWeightFsrDown = raw.psWeightFsrDefDown;
	psWeightIsrDown	= raw.psWeightIsrDefDown;
	psWeightFsrUp   = raw.psWeightFsrDefUp;
	psWeightIsrUp	= raw.psWeightIsrDefUp;
         //if(preVFP==true and RunYear==2016) PUweight = getPUESF16preVFP(nVtx);
         //if(preVFP==false and RunYear==2016) PUweight = getPUESF16(nVtx);
         //if(RunYear==2017) PUweight = getPUESF17(nVtx);
         //if(RunYear==2018) PUweight = getPUESF18(nVtx);

 	
         if(raw.nPho <1)continue;
 
	      int Nmedpho(0);
	      for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
 	        if(!itpho->isMedium())continue;
 	        if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
 	        if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
 	        if(itpho->getSigma() < 0.005)continue;
 	        if(itpho->PixelSeed() != 0)continue;
 	
 	        Nmedpho+=1;
	      }
	      if(Nmedpho >= 2)continue;
 
				nBJet = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(itJet->getPt() < 20)continue;
					if(itJet->isBJet())nBJet+=1;
				}
				pdfWeight = raw.pdfWeight;
				pdfSystWeight.clear();     // legacy — kept empty
				ScaleSystWeight.clear();   // legacy — kept empty

				// Read scalar branches written by ggNtuplizer
				muRweightUp   = raw.muRweightUp;
				muRweightDown = raw.muRweightDown;
				muFweightUp   = raw.muFweightUp;
				muFweightDown = raw.muFweightDown;
				pdfWeightRMS  = raw.pdfWeightRMS;
				//process MC first  ****************************/
				llmass = -1;
				double WGPt = -1;
				if(datatype == MC){
				 std::vector<mcData>::iterator posMcMu; 
				 bool hasposMcMu(false);
				 std::vector<mcData>::iterator negMcMu;
				 bool hasnegMcMu(false);
				 std::vector<mcData>::iterator posMcTau; 
				 bool hasposMcTau(false);
				 std::vector<mcData>::iterator negMcTau;
				 bool hasnegMcTau(false);
				 std::vector<mcData>::iterator posMcEle; 
				 bool hasposMcEle(false);
				 std::vector<mcData>::iterator negMcEle;
				 bool hasnegMcEle(false);
				 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
 					switch(itMC->getPID()){
 						case 13: if(isHardLepton(itMC->getMomPID())){posMcMu = itMC; hasposMcMu=true;};break;
 						case -13:if(isHardLepton(itMC->getMomPID())){negMcMu = itMC; hasnegMcMu=true;};break;
 						case 15: if(isHardLepton(itMC->getMomPID())){posMcTau = itMC; hasposMcTau=true;};break;
 						case -15:if(isHardLepton(itMC->getMomPID())){negMcTau = itMC; hasnegMcTau=true;};break;
 						case 11: if(isHardLepton(itMC->getMomPID())){posMcEle = itMC; hasposMcEle=true;};break;
 						case -11:if(isHardLepton(itMC->getMomPID())){negMcEle = itMC; hasnegMcEle=true;};break;
 						case 23: WGPt = itMC->getPt(); break;
 						case 24: WGPt = itMC->getPt(); break;
 						case -24: WGPt = itMC->getPt(); break;
 						default: break;
 					}
				 }
				 if(hasposMcMu && hasnegMcMu)llmass = (posMcMu->getP4() + negMcMu->getP4()).M();
				 else if(hasposMcTau && hasnegMcTau)llmass = (posMcTau->getP4() + negMcTau->getP4()).M();
				 else if(hasposMcEle && hasnegMcEle)llmass = (posMcEle->getP4() + negMcEle->getP4()).M();
				}
 
         bool hasegPho(false);
         bool hasmgPho(false);
         std::vector<recoPhoton>::iterator egsignalPho = Photon.begin();
         std::vector<recoPhoton>::iterator mgsignalPho = Photon.begin();
       	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
					if(itpho->getR9() < 0.5)continue;
          if(!itpho->passSignalSelection())continue;
					bool PixelVeto = itpho->PixelSeed()==0? true: false;
					bool GSFveto(true);
					bool FSRVeto(true);
					for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
					 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
					}
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
						if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
					if(GSFveto && PixelVeto && FSRVeto){
					bool passEGPhoTrig = (RunYear==2016)
					                    ? (itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6))
					                    : (itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6));
					bool passMGPhoTrig = (RunYear==2016)
					                    ? (itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30))
					                    : (itpho->fireDoubleTrg(27) || itpho->fireDoubleTrg(36));

					if(passEGPhoTrig && !hasegPho){
						hasegPho=true;
						egsignalPho = itpho;
					}
					if(passMGPhoTrig && !hasmgPho){
						hasmgPho=true;
						mgsignalPho = itpho;
					}
					}
         }
      
				if(hasegPho){ 
					if(egsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasegPho=false;
					if(egsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasegPho=false;
				}
				if(hasmgPho){ 
					if(mgsignalPho->getCalibEt() > 50  && mcType==MCType::WGJetInclusive)hasmgPho=false;
					if(mgsignalPho->getCalibEt() <= 50 && mcType==MCType::WGJet40)hasmgPho=false;
				}
  
         		bool hasEle(false);
         		std::vector<recoEle>::iterator signalEle = Ele.begin();
			if(hasegPho){
         			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
          				if(hasEle)break;
					if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
					bool passEleTrig = (RunYear==2016) ? (itEle->fireTrgs(21) || itEle->fireTrgs(22))
					                                  : (itEle->fireTrgs(43) || itEle->fireTrgs(44));
					if(passEleTrig){
						if(itEle->passSignalSelection()){
							hasEle=true; 
							signalEle = itEle;
						}
					}
         			}
			}
         		bool hasMu(false);
         		std::vector<recoMuon>::iterator signalMu = Muon.begin();
        		if(hasmgPho){
				for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
					if(hasMu)break;
						bool passMuTrig = (RunYear==2016)
						               ? (itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22))
						               : (itMu->fireSingleTrg(0) || itMu->fireSingleTrg(32));
						if(passMuTrig){
							if(itMu->passSignalSelection()){
								hasMu=true; 
								signalMu = itMu;
							}
						}
				}
         		}
 
				bool hasTrail(false);
				std::vector<recoMuon>::iterator trailLep = Muon.begin();
				for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
					if(itMu->getPt() < 25)continue;
					if(itMu->passSignalSelection()){
						if(!hasTrail){
							hasTrail = true;
							trailLep = itMu;
						}
					}
				}
 

				bool hasDoubleEG(false);
 	
				ISRJetPt = 0;
				TLorentzVector JetVec(0,0,0,0);	
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))JetVec = JetVec + itJet->getP4();
				}
				ISRJetPt = JetVec.Pt();
	///*************************   eg filters *****************************//
         if(hasegPho && hasEle && (((raw.HLTPho >> 14)&1)==1)){ // Year = 2016 Year = 2017
          double dReg = DeltaR(egsignalPho->getEta(), egsignalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
	  if(dReg>0.8){
	     if (fabs((egsignalPho->getP4() + signalEle->getP4()).M() - 91.188) > 10.0) {
              if(raw.passMETFilter(RunYear, METFilter)){

								hasDoubleEG = true;
								float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
								float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
        	    						eg_phoEt = egsignalPho->getEt();
			    					eg_phoEta= egsignalPho->getEta();
                						eg_phoPhi= egsignalPho->getPhi();
			    					eg_lepPt = signalEle->getPt();
 			    					eg_lepEta= signalEle->getEta();
                						eg_lepPhi= signalEle->getPhi();
								if(hasTrail){
									trailPt = trailLep->getPt();
									trailEta = trailLep->getEta();
									trailPhi = trailLep->getPhi();
								}
								else{
									trailPt  = 0; 
									trailEta = 0;
									trailPhi = 0;
								}
 			    					eg_sigMT = MT;
			    					eg_sigMET= MET;
                						eg_sigMETPhi = METPhi;
			    					eg_dPhiLepMET = deltaPhi; 
			    					eg_nVertex = nVtx; 
			    					eg_dRPhoLep= dReg;
                						eg_invmass = (egsignalPho->getP4()+signalEle->getP4()).M();
								eg_bosonPt = WGPt;
								eg_sigMETJESup = MET_T1JESUp;
								eg_sigMETJESdo = MET_T1JESDo;
								eg_sigMETJERup = MET_T1JERUp;
								eg_sigMETJERdo = MET_T1JERDo;
								eg_dPhiLepMETJESup = DeltaPhi(signalEle->getPhi(), METPhi_T1JESUp);
								eg_dPhiLepMETJESdo = DeltaPhi(signalEle->getPhi(), METPhi_T1JESDo);
								eg_dPhiLepMETJERup = deltaPhi;
								eg_dPhiLepMETJERdo = deltaPhi;
								eg_sigMTJESup = sqrt(2*MET_T1JESUp*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESup)));
								eg_sigMTJESdo = sqrt(2*MET_T1JESDo*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJESdo)));
								eg_sigMTJERup = sqrt(2*MET_T1JERUp*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERup)));
								eg_sigMTJERdo = sqrt(2*MET_T1JERDo*signalEle->getPt()*(1-std::cos(eg_dPhiLepMETJERdo)));
 							
								eg_nJet = 0;
								eg_nISRJet = 0;
								eg_HT = 0;
								eg_HTJESup = 0;
								eg_HTJESdo = 0;
								TLorentzVector egISRVec(0,0,0,0);
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
 									
									if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))eg_nISRJet += 1;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(), signalEle->getPhi()) <= 0.4)continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), egsignalPho->getEta(),egsignalPho->getPhi()) <= 0.4)continue;	
									eg_nJet += 1;
									eg_HT += itJet->getPt();
									eg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
									eg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
									egISRVec += itJet->getP4();
								}
								// ISR pT as in the derivation (ISRReweight/analysis_ISRMC.C calcISR); old gen-matched one kept
								eg_ISRJetPt = egISRVec.Pt();
								eg_ISRJetPtGen = ISRJetPt;
								crosssection = genWeightSign*(isWGSample ? getWGKfactor(eg_phoEt)*MC_XS[mcType] : MC_XS[mcType]);
								if(datatype == MC){
									eg_mcPID.clear();
									eg_mcEta.clear();
									eg_mcPhi.clear();
									eg_mcPt.clear();
									eg_mcMomPID.clear();
									eg_mcGMomPID.clear();
									eg_mcStatus.clear();
									for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
										if(itMC->getEt() < 5.0)continue;
										eg_mcPID.push_back(itMC->getPID());
										eg_mcMomPID.push_back(itMC->getMomPID());
										eg_mcEta.push_back(itMC->getEta());      
										eg_mcPhi.push_back(itMC->getPhi());
										eg_mcPt.push_back(itMC->getEt());
										eg_mcGMomPID.push_back(itMC->getGMomPID());
										eg_mcStatus.push_back(itMC->getStatus());
									}
								}
                						egtree->Fill();
  
               }//MET Filter
             }// Z mass Filter
           }//dR filter
         }// ele + pho candidate
        
//**********************  mg filter **************************************         
  
       bool passMGHLT = (RunYear==2016) ? ((((raw.HLTEleMuX >> 8)&1)!=0) || (((raw.HLTEleMuX >> 51)&1)!=0))
                                        : ((((raw.HLTEleMuX >> 8)&1)!=0) || (((raw.HLTEleMuX >> 57)&1)!=0));
       if(hasmgPho && hasMu && passMGHLT && !hasDoubleEG){
 					double dRmg = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), signalMu->getEta(), signalMu->getPhi());
 					if(dRmg>0.8){
 						if(raw.passMETFilter(RunYear, METFilter)){ 
 
 						float deltaPhi = DeltaPhi(signalMu->getPhi(), METPhi);
 						float MT = sqrt(2*MET*signalMu->getPt()*(1-std::cos(deltaPhi)));
 						float ThreeBodyMass = sqrt(2*MET*(mgsignalPho->getP4()+ signalMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (mgsignalPho->getP4()+signalMu->getP4()).Phi(), 0, METPhi))));
 
 						mg_phoEt = mgsignalPho->getCalibEt();
 						mg_phoEta= mgsignalPho->getEta();
 						mg_phoPhi= mgsignalPho->getPhi();
 						mg_lepPt = signalMu->getPt();
 						mg_lepEta= signalMu->getEta();
 						mg_lepPhi= signalMu->getPhi();
 						mg_sigMT = MT;
 						mg_sigMET= MET;
 						mg_sigMETPhi = METPhi;
 						mg_dPhiLepMET = deltaPhi;
 						mg_threeMass = ThreeBodyMass;
 						mg_nVertex = nVtx;
 						mg_dRPhoLep= dRmg;
 						mg_bosonPt = WGPt;
 						mg_sigMETJESup = MET_T1JESUp;
 						mg_sigMETJESdo = MET_T1JESDo;
						mg_sigMETJERup = MET_T1JERUp;
                                                mg_sigMETJERdo = MET_T1JERDo;
                                                mg_dPhiLepMETJESup = DeltaPhi(signalMu->getPhi(), METPhi_T1JESUp);
                                                mg_dPhiLepMETJESdo = DeltaPhi(signalMu->getPhi(), METPhi_T1JESDo);
                                                mg_dPhiLepMETJERup = deltaPhi;
                                                mg_dPhiLepMETJERdo = deltaPhi;
                                                mg_sigMTJESup = sqrt(2*MET_T1JESUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESup)));
                                                mg_sigMTJESdo = sqrt(2*MET_T1JESDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJESdo)));
                                                mg_sigMTJERup = sqrt(2*MET_T1JERUp*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERup)));
                                                mg_sigMTJERdo = sqrt(2*MET_T1JERDo*signalMu->getPt()*(1-std::cos(mg_dPhiLepMETJERdo)));

                                                mg_nJet = 0;
                                                mg_nISRJet = 0;
                                                mg_HT = 0;
                                                mg_HTJESup = 0;
                                                mg_HTJESdo = 0;
                                                TLorentzVector mgISRVec(0,0,0,0);
                                                for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
                                                        if(!itJet->passSignalSelection())continue;
                                                        if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))mg_nISRJet += 1;
                                                        if(DeltaR(itJet->getEta(), itJet->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= 0.4)continue;
                                                        if(DeltaR(itJet->getEta(), itJet->getPhi(), mgsignalPho->getEta(),mgsignalPho->getPhi()) <= 0.4)continue;
                                                        mg_nJet += 1;
                                                        mg_HT += itJet->getPt();
                                                        mg_HTJESup += itJet->getPt()*(1+itJet->getPtUnc());
                                                        mg_HTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
                                                        mgISRVec += itJet->getP4();
                                                }
						// ISR pT as in the derivation (ISRReweight/analysis_ISRMC.C calcISR); old gen-matched one kept
 						mg_ISRJetPt = mgISRVec.Pt();
						mg_ISRJetPtGen = ISRJetPt;
						crosssection = genWeightSign*(isWGSample ? getWGKfactor(mg_phoEt)*MC_XS[mcType] : MC_XS[mcType]);
 						
						if(datatype == MC){
 						 	mg_mcPID.clear();
 						 	mg_mcEta.clear();
 						 	mg_mcPhi.clear();
 						 	mg_mcPt.clear();
 						 	mg_mcMomPID.clear();
 						 	mg_mcGMomPID.clear();
 						 	mg_mcStatus.clear();
 						 //for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
 						 //	if(itMC->getEt() < 1.0)continue;
 						 //	float mcdRmu = DeltaR(signalMu->getEta(), signalMu->getPhi(), itMC->getEta(), itMC->getPhi());
 						 //	float mcdR = DeltaR(mgsignalPho->getEta(), mgsignalPho->getPhi(), itMC->getEta(), itMC->getPhi());
 						 //	if(mcdR < 0.3 || mcdRmu < 0.3){
						 	for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
                                                                if(itMC->getEt() < 5.0)continue;
 						 		mg_mcPID.push_back(itMC->getPID());
 						 		mg_mcMomPID.push_back(itMC->getMomPID());
 						 		mg_mcGMomPID.push_back(itMC->getGMomPID());
 						 		mg_mcEta.push_back(itMC->getEta());
 						 		mg_mcPhi.push_back(itMC->getPhi());
 						 		mg_mcPt.push_back(itMC->getEt());
 								mg_mcStatus.push_back(itMC->getStatus());
 						 //	}
 						  	}
 						}
 				 		mgtree->Fill();
            }//MET Filter
          }//dR Filter
        }//Candidate Filter
  
 	}//loop on  events
 
cout<< "egtree: " << egtree->GetEntries() << "\tmgtree: "<< mgtree->GetEntries() <<endl;
outputfile->Write();
outputfile->Close();
logfile.close();
 
 }
 
 int main(int argc, char** argv)
 {
    if(argc < 4){
      cout << "Usage: ./analysis_VGamma_Run2.exe <RunYear> <preVFP:0|1> <Sample> [norm|range <first> <last> <normFile>]\n";
      return 1;
    }
    if(argc > 7) analysis_VGamma_Run2(atoi(argv[1]), atoi(argv[2]) == 1, argv[3], argv[4], atoll(argv[5]), atoll(argv[6]), argv[7]);
    else         analysis_VGamma_Run2(atoi(argv[1]), atoi(argv[2]) == 1, argv[3]);
     return 0;
 }
