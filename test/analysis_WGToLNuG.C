#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include <TROOT.h>

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
#include "TVector2.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_jet.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"


void analysis_WGToLNuG() { // main
  gSystem->Load("../lib/libAnaClasses.so");
  gROOT->SetBatch(kTRUE);

  // -------------------- Input Files --------------------
  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("/eos/uscms/store/user/tmishra/InputFilesMC/WGToLNuG/WGToLNuG_2018.root");
  es->Add("/eos/uscms/store/user/tmishra/InputFilesMC/ZGToLLG/ZGToLLG_2018.root");
  es->Add("/eos/uscms/store/user/tmishra/InputFilesMC/WGJet40/WGJet40_2018.root");
  es->Add("/eos/uscms/store/user/tmishra/InputFilesMC/WGJet130/WGJet130_2018.root");

  std::ostringstream outputname;
  outputname << "/eos/uscms/store/user/tmishra/VGamma_2018.root";
  TFile *outputfile = TFile::Open(outputname.str().c_str(), "RECREATE");
  outputfile->cd();

  const unsigned nEvts = es->GetEntries();
  std::cout << "Total events: " << nEvts << std::endl;

  RunType datatype(MC);
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoEle> Ele;
  std::vector<recoMuon> Muon;



TH1F* h_genMatchedMuonPt  = new TH1F("h_genMatchedMuonPt",  ";p_{T} [GeV];Normalised Events", 750, 50, 800);
TH1F* h_genMatchedEleCalibPt = new TH1F("h_genMatchedEleCalibPt", ";p_{T}^{e, calib} [GeV];Normalised Events", 750, 50, 800);
TH1F* h_genMatchedEleRecoPt  = new TH1F("h_genMatchedEleRecoPt",  ";p_{T}^{e, reco} [GeV];Normalised Events", 750, 50, 800);

  // -------------------- Event Loop --------------------
  for (unsigned ievt = 0; ievt < nEvts; ++ievt) {
    if (ievt % 1000000 == 0) std::cout << " -- Processing event " << ievt << std::endl;

    raw.GetData(es, ievt);
    MCData.clear();
    Muon.clear();
    Ele.clear();

    for (int iEle = 0; iEle < raw.nEle; iEle++) Ele.push_back(recoEle(raw, iEle));
    for (int iMC  = 0; iMC  < raw.nMC;  iMC++)  MCData.push_back(mcData(raw, iMC));
    for (int iMu  = 0; iMu  < raw.nMu;  iMu++)  Muon.push_back(recoMuon(raw, iMu));

    for (auto itMC = MCData.begin(); itMC != MCData.end(); ++itMC) {
      int pid = fabs(itMC->getPID());

      // -------------------- Electron match (pid == 11) --------------------
      if (pid == 11) {
        float mindR = 0.1;
        for (auto itEle = Ele.begin(); itEle != Ele.end(); ++itEle) {
          float dR = DeltaR(itMC->getEta(), itMC->getPhi(), itEle->getEta(), itEle->getPhi());
          if (dR < mindR) {
            if (itEle->getCalibPt() > 50)
              h_genMatchedEleCalibPt->Fill(itEle->getCalibPt());
            if (itEle->getPt() > 50)
              h_genMatchedEleRecoPt->Fill(itEle->getPt());
            break;
          }
        }
      }

      // -------------------- Muon match (pid == 13) --------------------
      if (pid == 13) {
        float mindR = 0.1;
        for (auto itMu = Muon.begin(); itMu != Muon.end(); ++itMu) {
          if (itMu->getPt() < 50) continue;
          float dR = DeltaR(itMC->getEta(), itMC->getPhi(), itMu->getEta(), itMu->getPhi());
          if (dR < mindR) {
            h_genMatchedMuonPt->Fill(itMu->getPt());
            break;
          }
        }
      }
    }
  } // end event loop

  // -------------------- Normalize histograms --------------------
  auto normalize = [](TH1F* h) { if (h->Integral() > 0) h->Scale(1.0 / h->Integral()); };
  normalize(h_genMatchedEleCalibPt);
  normalize(h_genMatchedEleRecoPt);
  normalize(h_genMatchedMuonPt);

  // -------------------- Save histograms --------------------
  outputfile->cd();
  h_genMatchedEleCalibPt->Write();
  h_genMatchedEleRecoPt->Write();
  h_genMatchedMuonPt->Write();

  // -------------------- Canvas 1: Reco pT (electron vs muon) --------------------
  TCanvas* c1 = new TCanvas("c1", "Reco pT: Electron vs Muon", 800, 700);
  h_genMatchedMuonPt->SetLineColor(kBlue);
  h_genMatchedEleRecoPt->SetLineColor(kRed);
  h_genMatchedMuonPt->SetStats(0);
  h_genMatchedEleRecoPt->SetStats(0);
  h_genMatchedMuonPt->Draw("HIST");
  h_genMatchedEleRecoPt->Draw("HIST SAME");
  TLegend* leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg1->AddEntry(h_genMatchedEleRecoPt, "Reco Electron", "l");
  leg1->AddEntry(h_genMatchedMuonPt, "Reco Muon", "l");
  leg1->Draw();
  c1->SetLogy();
  c1->SaveAs("RecoPt_EleMuon.png");

  // -------------------- Canvas 2: Calib pT (electron) vs Reco pT (muon) --------------------
  TCanvas* c2 = new TCanvas("c2", "Calib Electron vs Reco Muon", 800, 700);
  h_genMatchedMuonPt->SetLineColor(kBlue);
  h_genMatchedEleCalibPt->SetLineColor(kMagenta);
  h_genMatchedEleRecoPt->SetLineColor(kRed);
  h_genMatchedMuonPt->Draw("HIST");
  h_genMatchedEleCalibPt->Draw("HIST SAME");
  h_genMatchedEleRecoPt->Draw("HIST SAME");
  TLegend* leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg2->AddEntry(h_genMatchedEleCalibPt, "Calib Electron", "l");
  leg2->AddEntry(h_genMatchedEleRecoPt, "Reco Electron", "l");
  leg2->AddEntry(h_genMatchedMuonPt, "Reco Muon", "l");
  leg2->Draw();
  c2->SetLogy();
  c2->SaveAs("CalibEle_vs_RecoMuon.png");

  // -------------------- Close output file --------------------
  outputfile->Close();
}
