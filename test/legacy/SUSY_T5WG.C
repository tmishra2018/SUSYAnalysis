#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

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

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
        float deltaPhi = TMath::Abs(phi1-phi2);
        float deltaEta = eta1-eta2;
        if(deltaPhi > TMath::Pi())
        deltaPhi = TMath::TwoPi() - deltaPhi;
                return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
int RunYear = 2017;


void SUSY_T5WG(){//main  
  gROOT->SetBatch();
  TFile f(Form("T5Wg_%d_skim.root",RunYear));
  TTree * es = (TTree*)f.Get("Events/Events");
  
  UInt_t          nJet,nElectron,nPhoton,nMuon;
  Float_t Electron_eta[8], Electron_pt[8], Electron_phi[8], Electron_mass[8];
  Float_t Muon_eta[8], Muon_pt[8], Muon_phi[8], Muon_mass[8];
  Bool_t          Flag_goodVertices;
  Int_t           PV_npvsGood;  

  es->SetBranchAddress("nElectron", &nElectron);
  es->SetBranchAddress("nMuon", &nMuon);
  es->SetBranchAddress("nPhoton", &nPhoton);
  es->SetBranchAddress("nJet", &nJet);
/*  es->SetBranchAddress("Electron_pt", &Electron_pt);
  es->SetBranchAddress("Electron_eta", &Electron_eta);
  es->SetBranchAddress("Electron_phi", &Electron_phi);
  es->SetBranchAddress("Electron_mass", &Electron_mass);
  es->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
  es->SetBranchAddress("PV_npvsGood", &PV_npvsGood);

  TLorentzVector electron;
  
  TLorentzVector photon;
  Float_t Photon_eta[13], Photon_pt[13], Photon_phi[13], Photon_mass[13];
  es->SetBranchAddress("Photon_pt", &Photon_pt);
  es->SetBranchAddress("Photon_eta", &Photon_eta);
  es->SetBranchAddress("Photon_phi", &Photon_phi);
  es->SetBranchAddress("Photon_mass", &Photon_mass);

  TLorentzVector muon;
  es->SetBranchAddress("Muon_pt", &Muon_pt);
  es->SetBranchAddress("Muon_eta", &Muon_eta);
  es->SetBranchAddress("Muon_phi", &Muon_phi);
  es->SetBranchAddress("Muon_mass", &Muon_mass);

  Float_t invMass,deltaR; 
*/
  TH1F *InvMass_H_T5Wg = new TH1F ("nJet_H",";nJet;Number of events",40,-5,35);

  //TH1F * InvMass_H_T5Wg = new TH1F ("InvMass_H_T5Wg",";Invariant Mass ;Number of events",300,0,5000);
  
    std::cout << "total event : " << es->GetEntries() << std::endl;
    
    for(Int_t i=0;i<es->GetEntries();i++){
	    es->GetEntry(i);
//	    if(PV_npvsGood < 1) continue;       

//            if(nPhoton != 1) continue; 
 //           if(Photon_pt[0] <= 35) continue;

//	    if(nMuon+nElectron < 1) continue;
//	    if(Electron_pt[0] <= 25 && Muon_pt[0] <= 25) continue;

//	    electron.SetPtEtaPhiE(Electron_pt[0],Electron_eta[0],Electron_phi[0],Electron_mass[0]);
//	    photon.SetPtEtaPhiE(Photon_pt[0],Photon_eta[0],Photon_phi[0],Photon_mass[0]);
//	    invMass = (electron+photon).M(); 
//	    deltaR = DeltaR(Photon_eta[0],Photon_phi[0],Electron_eta[0],Electron_phi[0]);
	    
//	    if(deltaR <= 0.8) continue;
//	    if (!(abs(invMass)-91.188>10.0)) continue; 
	    //InvMass_H_T5Wg->Fill(abs(invMass));
	    InvMass_H_T5Wg->Fill(nJet);
    }
    



 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  TFile f1(Form("TChiWG_%d_skim.root",RunYear));
  TTree * es1 = (TTree*)f1.Get("Events/Events");
  
  UInt_t          nJet1,nElectron1,nPhoton1,nMuon1;
  Bool_t          Flag_goodVertices1;
  Int_t           PV_npvsGood1;  

  Float_t Electron_eta1[8], Electron_pt1[8], Electron_phi1[8], Electron_mass1[8];
  Float_t Muon_eta1[8], Muon_pt1[8], Muon_phi1[8], Muon_mass1[8];
  
  es1->SetBranchAddress("nJet", &nJet1);
  es1->SetBranchAddress("nElectron", &nElectron1);
  es1->SetBranchAddress("nMuon", &nMuon1);
  es1->SetBranchAddress("nPhoton", &nPhoton1);
/*
  es1->SetBranchAddress("Electron_pt", &Electron_pt1);
  es1->SetBranchAddress("Electron_eta", &Electron_eta1);
  es1->SetBranchAddress("Electron_phi", &Electron_phi1);
  es1->SetBranchAddress("Electron_mass", &Electron_mass1);
  TLorentzVector electron1;
  
  TLorentzVector photon1;
  Float_t Photon_eta1[13], Photon_pt1[13], Photon_phi1[13], Photon_mass1[13];
  es1->SetBranchAddress("Photon_pt", &Photon_pt1);
  es1->SetBranchAddress("Photon_eta", &Photon_eta1);
  es1->SetBranchAddress("Photon_phi", &Photon_phi1);
  es1->SetBranchAddress("Photon_mass", &Photon_mass1);
  es1->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices1);
  es1->SetBranchAddress("PV_npvsGood", &PV_npvsGood1);

  TLorentzVector muon1;
  es1->SetBranchAddress("Muon_pt", &Muon_pt1);
  es1->SetBranchAddress("Muon_eta", &Muon_eta1);
  es1->SetBranchAddress("Muon_phi", &Muon_phi1);
  es1->SetBranchAddress("Muon_mass", &Muon_mass1);

  Float_t invMass1, deltaR1; */
  TH1F *InvMass_H_TChiWG = new TH1F ("nJet_H",";nJet;Number of events",40,-5,35);
  //TH1F * InvMass_H_TChiWG = new TH1F ("InvMass_H_TChiWG",";Invariant Mass ;Number of events",300,0,5000);

    for(Int_t j=0;j<es1->GetEntries();j++){
	    es1->GetEntry(j);
//	    if(PV_npvsGood1 < 1) continue;       
 
//            if(nPhoton1 != 1) continue; 
//            if(Photon_pt1[0] <= 35) continue;  

//	    if(nMuon1+nElectron1 < 1) continue;
//	    if(Electron_pt1[0] <= 25 && Muon_pt1[0] <= 25) continue;

//	    electron1.SetPtEtaPhiE(Electron_pt1[0],Electron_eta1[0],Electron_phi1[0],Electron_mass1[0]);
//	    photon1.SetPtEtaPhiE(Photon_pt1[0],Photon_eta1[0],Photon_phi1[0],Photon_mass1[0]);
//	    invMass1 = (electron1+photon1).M();
//	    deltaR1 = DeltaR(Photon_eta1[0],Photon_phi1[0],Electron_eta1[0],Electron_phi1[0]);
	    
//	    if(deltaR1<=0.8) continue;
//	    if (!(abs(invMass1)-91.188>10.0)) continue; 
	    //InvMass_H_TChiWG->Fill(abs(invMass1));
	    InvMass_H_TChiWG->Fill(nJet1);
    }

    TCanvas* can = new TCanvas("can","", 1400, 1200);
    can->cd();   
    InvMass_H_T5Wg->SetStats(0);
    InvMass_H_TChiWG->SetStats(0);

    can->SetLogy();
    InvMass_H_T5Wg->Scale(1/InvMass_H_T5Wg->Integral());
    InvMass_H_TChiWG->Scale(1/InvMass_H_TChiWG->Integral());
    InvMass_H_T5Wg->SetMaximum(10);
    
    InvMass_H_T5Wg->SetLineColor(kRed);
    InvMass_H_TChiWG->SetLineColor(kBlue);
    //InvMass_H_T5Wg->SetTitle(Form("%d",RunYear));
    //InvMass_H_T5Wg->SetTitleSize(1.0);
    InvMass_H_T5Wg->Draw("hist"); 
    TPaveStats *stats = (TPaveStats*)can->GetPrimitive("stats");

    InvMass_H_TChiWG->Draw("hist same"); 

    TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    leg->SetFillColor(0);
    leg->AddEntry(InvMass_H_TChiWG,"TChiWG","l");
    leg->AddEntry(InvMass_H_T5Wg,"T5Wg","l");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();



    can->SaveAs(Form("try_%d.pdf",RunYear)); 

}
