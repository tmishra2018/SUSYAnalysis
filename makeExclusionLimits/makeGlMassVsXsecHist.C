#include<iostream>
#include<iomanip>
#include<string>
#include"TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"

void makeGlMassVsXsecHist(){
  gStyle->SetOptStat(0);
  vector<double> mass,xsec,pcUnc,absUnc;
  ifstream infile("mGl_Xsecpb_pcUnc_absUnc.txt");
  string line;
  while (getline(infile, line)){
    std::istringstream iss(line);
    double mass1,xsec1,pcUnc1,absUnc1;
    iss>>mass1;
    iss>>xsec1;
    iss>>pcUnc1;
    iss>>absUnc1;
    mass.push_back(mass1);
    xsec.push_back(xsec1);
    pcUnc.push_back(pcUnc1);
    absUnc.push_back(absUnc1);
  }
  // TH1D *h1 ;//= new TH1D("mGlXsec","x:mass of gluino in GeV, bin content(y-axis) is xsec in pb and bin error is xsec uncertainty in pb",3011,-0.5,3010.5);
  // TFile *f1 = new TFile("T5bbbbZg_MassScan.root","recreate");
  TFile *f1 = TFile::Open("T5bbbbZg_MassScan.root");
  TH1D *h1 = (TH1D*)f1->FindObjectAny("mGlXsec");
  for(int i=0;i<mass.size();i++){
    h1->Fill(mass[i],xsec[i]);
    h1->SetBinError(h1->FindBin(mass[i]),absUnc[i]);
    cout<<mass[i]<<"    "<<xsec[i]<<endl;
  }

  h1->GetXaxis()->SetTitle("Gluino Mass(GeV)");
  h1->GetYaxis()->SetTitle("#sigma(pb)");
  h1->Draw();
//TFile *f1=new TFile("mGl_Xsecpb_absUnc.root","recreate");
  h1->Write();


  // if(0){
  //   int np=mass.size();
  //   if(np!=xsec.size() || np!=pcUnc.size() || np!=absUnc.size()){ cout<<"problem in getting mass and xsec"<<endl; return;}
  //   cout<<"mGl_1"<<endl<<"{";
  //   for(int i=0;i<np;i++){
  //     cout<<mass[i]<<",";
  //   }
  //   cout<<endl<<"xsec_1"<<endl<<"{";
  //   for(int i=0;i<np;i++){
  //     cout<<xsec[i]<<",";
  //   }
  //   cout<<endl<<"Uncxsec_1"<<endl<<"{";
  //   for(int i=0;i<np;i++){
  //     cout<<absUnc[i]<<",";
  //   }
  //   cout<<endl;
  //  }
}


// TFile *f1 = new TFile("T5bbbbZg_MassScan.root");
// TH1D *h1_xsechist = (TH1D*)f1->FindObjectAny("mGlXsec");
// TH2D *h2_mass = (TH2D*)f1->FindObjectAny("MGlMNLSP");
// if(!h1_xsechist){ cout<<"AHHHHH: could not find xsec hist"<<endl; return;}
// if(!h2_mass){ cout<<"AHHHHH: could not find no. of events hist"<<endl; return;}
// xsec1 = h1_xsechist->GetBinContent(h1_xsechist->FindBin(momMass));
