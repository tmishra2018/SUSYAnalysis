#include"TFile.h"
#include"TH1.h"
#include"TVirtualPad.h"
void overlay(){
  // TFile *f1=new TFile("/home/vinay/work/ROOT_Files/Physics/GMSB_susy/GMSB_skims_ST_RA2b_TreesV12/SignalRegion_v1/SR_Summer16.WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  // TFile *f2=new TFile("/home/vinay/work/ROOT_Files/Physics/GMSB_susy/GMSB_skims_ST_RA2b_TreesV11/SignalRegion/SR_Spring16.WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");

  // TFile *f1=new TFile("/home/vinay/work/ROOT_Files/Physics/GMSB_susy/GMSB_skims_ST_RA2b_TreesV12/SignalRegion_v1/SR_Summer16.TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  // TFile *f2=new TFile("/home/vinay/work/ROOT_Files/Physics/GMSB_susy/GMSB_skims_ST_RA2b_TreesV11/SignalRegion/SR_Spring16.TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
// gPad->SetGridx(2);

   TCanvas *c1 = new TCanvas("c1", "c1",65,24,1680,1056);
   c1->Range(-1.981819,-20.70421,12.90736,118.2711);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetBottomMargin(0.1489776);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);

   c1->cd();
  gStyle->SetOptStat(10);  

   TFile *f1=new TFile("../rootoutput/T5bbbbZg_1800_150_FastSim_v17.root");
   TFile *f2=new TFile("../rootoutput/T5bbbbZg_1800_150_FastSim_v12.root");
   //TFile *f1=new TFile("../rootoutput/ZGZJ_NuNuG.root");
   //TFile *f2=new TFile("../rootoutput/ZGZJ_NuNuG_v12.root");
  //  TFile *f2=new TFile("/home/vinay/work/ROOT_Files/Physics/GMSB_susy/GMSB_skims_ST_RA2b_TreesV11/SignalRegion/SR_Spring16.QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  

  // TH1F *h1=(TH1F*)f1->FindObjectAny("selectBaselineYields_");
  // // TH1F *h2=(TH1F*)f1->FindObjectAny("selectBaselineYields_");

  TH1F *h1=(TH1F*)f1->FindObjectAny("cutflows");
  TH1F *h2=(TH1F*)f2->FindObjectAny("cutflows");
  // TH1F *h1=(TH1F*)f1->FindObjectAny("nJets");
  // TH1F *h2=(TH1F*)f2->FindObjectAny("nJets");
  // TH1F *h1=(TH1F*)f1->FindObjectAny("MET");
  // TH1F *h2=(TH1F*)f2->FindObjectAny("MET");
  // TH1F *h1=(TH1F*)f1->FindObjectAny("BestPhotonPt");
  // TH1F *h2=(TH1F*)f2->FindObjectAny("BestPhotonPt");
  // TH1F *h1=(TH1F*)f1->FindObjectAny("nBTags");
  // TH1F *h2=(TH1F*)f2->FindObjectAny("nBTags");
  
   // TH1F *hnew=(TH1F*)h1->Clone("h1new");
  // hnew->Divide(h2);
  // hnew->Draw("text hist");
  // h1->SetTextColor(1);
  // h2->SetTextColor(2);


   TPaveStats *p = new TPaveStats(0.6303725,0.6934461,0.8710602,0.832135,"brNDC");
   p->SetName("stats");
   p->SetBorderSize(1);
   p->SetFillColor(0);
   p->SetTextAlign(12);
   p->SetTextFont(42);
   // p->GetLineWith("Entries")->SetTextColor(1);
   p->SetTextColor(1);
   p->SetLineColor(1);
   p->Draw();
   h1->GetListOfFunctions()->Add(p);
   p->SetParent(h1);
   h1->SetLineColor(1);
   h1->GetXaxis()->SetRange(1,11);
   h1->GetXaxis()->SetNdivisions(-510);
   h1->Draw("text");


   TPaveStats *p2 = new TPaveStats(0.6289398,0.517167,0.8724928,0.6511628,"brNDC");
   p2->SetName("stats");
   p2->SetBorderSize(1);
   p2->SetFillColor(0);
   p2->SetTextAlign(12);
   p2->SetTextFont(42);
   p2->SetTextColor(2);
   p2->SetLineColor(1);
//   //   p2->GetLineWith("Name")->SetTextColor(2);
   p2->Draw();
   h2->GetListOfFunctions()->Add(p2);
   p2->SetParent(h2);
   h2->SetLineColor(2);
   h2->Draw("hist sames");





  // TPaveStats *p = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
   // h1->GetListOfFunctions()->Remove(p);
   // h1->SetStats(0);
   // p->SetName("stats");
   // p->SetBorderSize(1);
   // p->SetFillColor(0);

  // TTree *t1 = (TTree*)f1->Get("PreSelection");
  // //  TTree *t2 = (TTree*)f2->Get("PreSelection");
  // double MET1,MET2;
  // t1->SetBranchAddress("MET",&MET1);
  // t2->SetBranchAddress("MET",&MET2);
  
  // // vector<TLorentzVector> v1,v2;
  // // t1->SetBranchAddress("Photons",&(v1[0]));
  // // t2->SetBranchAddress("Photons",&(v2[0]));
  
  // TH1D *ht1   = new TH1D("ht1","ht1 distribution",500,0,5000);
  // TH1D *ht2   = new TH1D("ht2","ht2 distribution",500,0,5000);

  // Long64_t nentries1 = t1->GetEntries();
  // for (Long64_t i=0;i<nentries1;i++) {
  //   t1->GetEntry(i);
  //   ht1->Fill(MET1);
  // }
  // Long64_t nentries2 = t2->GetEntries();
  // for (Long64_t i=0;i<nentries2;i++) {
  //   t2->GetEntry(i);
  //   ht2->Fill(MET2);
  // }

  // ht1->SetLineColor(kRed);
  // ht2->SetLineColor(kBlue);
  // ht1->DrawNormalized("e0");
  // ht2->DrawNormalized("sames e0");
}
