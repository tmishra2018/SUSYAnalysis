#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<stdlib.h>

using namespace std;
void getSignalEff(TString model){
  TFile *f1= TFile::Open(model+"_MassScan.root");
  TFile *f2;
  //  //  TString pathName="hists_datacards_limitTree_"+model+"_SbinV7/FastSim_"+model+"_",name;
  TString pathName="hists_datacards_limitTree_"+model+"_SbinV7_v2/FastSim_"+model+"_",name;
  //  TString pathName="root://cmseos.fnal.gov//store/user/vhegde/GMSB_skims_ST_RA2b_TreesV12/signalSystematics/hists_"+model+"/FastSim_"+model+"_bTagSF_",name;
  TFile *fout = new TFile("SignalEff_"+model+".root","recreate");

  TH1D *h1;
  cout<<name<<endl;

  TH2D *h2_den = (TH2D*)f1->FindObjectAny("MGlMNLSP");
  TH2D *h2_num = (TH2D*)h2_den->Clone("SignalEff");

  double entries2=-100,entriesNum=0.;
  int minGlMass = 0;

  for(int i=1;i<=h2_den->GetNbinsX();i++){
    for(int j=1;j<=h2_den->GetNbinsY();j++){
      h2_num->SetBinContent(i,j,0);
      int mp=h2_den->GetXaxis()->GetBinCenter(i);
      int md=h2_den->GetYaxis()->GetBinCenter(j);
      if(mp<minGlMass) continue;
      entries2=h2_den->GetBinContent(i,j);
      if(entries2<0.1) continue;
      name = pathName+to_string(mp)+"_"+to_string(md)+".root";
      f2 = TFile::Open(name);
      //      cout<<name<<endl;
      h1 = (TH1D*)f2->Get("nHadJets");
      entriesNum = h1->GetEntries();
      h2_num->SetBinContent(i,j,entriesNum);
    }
  }
  TH2D *h2_evtLeft = (TH2D*)h2_num->Clone("EventsPassed");
  // if(h2_den->GetNbinsX() > 500){ h2_den->RebinX(50); h2_num->RebinX(50); h2_evtLeft->RebinX(50);}
  // if(h2_den->GetNbinsY() > 500){ h2_den->RebinY(50); h2_num->RebinY(50); h2_evtLeft->RebinY(50);}
  h2_num->Divide(h2_den);
  
  fout->cd();
  h2_den->Write();
  h2_evtLeft->Write();
  h2_num->Write();

  cout<<h2_den->GetNbinsX()<<" "<<h2_num->GetNbinsX()<<" "<<h2_evtLeft->GetNbinsX()<<endl;

}
