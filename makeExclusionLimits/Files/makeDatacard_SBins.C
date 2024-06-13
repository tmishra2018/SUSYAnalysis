#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"

using namespace std;
string getfname(const char *fname1){string fname=fname1;fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();return fname;}
void setLastBinAsOverFlow(TH1D*);

void makeDatacard_SBins(double mGl,double mNLSP,TString sigFile){
  char name[1000];
  int jmax=8;//no. of backgrounds
  int nSig=1;//1 signal
  int nFiles=nSig+jmax;
  TFile *f[nFiles];
  f[0] = new TFile(sigFile);
  f[1] = new TFile("WGJetsToLNuG_v12.root");
  f[2] = new TFile("WJetsToLNu_v12.root");
  f[3] = new TFile("TTJetsHT_v12.root");
  f[4] = new TFile("TTGJets_v12.root");
  // f[5] = new TFile("ZGZJ_NuNuG_v12.root");
  f[5] = new TFile("ZGJetsToNuNuG_v12.root");
  f[6] = new TFile("ZJetsToNuNu_v12.root");
  f[7] = new TFile("GJets_v12.root");
  f[8] = new TFile("QCD_v12.root");

  //  char histname[100]="MET_R1";
  char histname[100]="AllSBins_v7_CD";
  TH1D *hist1[nFiles];
  //  double min_binLowedge=99.9999;
  ofstream outf;
  vector<int> observation;
  for(int i=0;i<nFiles;i++){
    hist1[i] = (TH1D*)f[i]->FindObjectAny(histname);
    //setLastBinAsOverFlow(hist1[i]);
  }

  TH1D *hist=(TH1D*)f[0]->FindObjectAny(histname);
  int imax=hist->GetNbinsX();
  for(int i=1;i<=imax;i++){observation.push_back(1);}
    for(int i=1;i<=imax;i++){
      cout<<"MEt bin --> "<<i<<" bin content --> "<<hist->GetBinContent(i)<<endl;
    }
    // cout<<"-----"<<getfname(f[j]->GetName())<<"-----"<<endl;}
  cout<<"file name : "<<getfname(f[0]->GetName())<<endl;
  for(int i=1;i<=imax;i++){
   
    // string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname+"/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    cout<<name2<<endl;
    sprintf(name,"%s",name2.c_str());
    
    outf.open(name);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      "# Datacard for mGl= "<<mGl<<" mNLSP= "<<mNLSP<<endl<<
      "# - - - - - - - - - - - - - - - - - - - "<<endl<<
      "imax 1 number of channels"<<endl<<
      "jmax *  number of backgrounds('*' = automatic)"<<endl<<
      "kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl<<
      "------------"<<endl<<
      "bin "<<histname<<i<<endl<<
      "observation "<<observation[i-1]<<endl<<
      "------------"<<endl<<
      "bin ";
    for(int j=0;j<jmax+nSig;j++){outf<<histname<<i<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<" ";}
    outf<<endl<<
      "rate ";
    for(int j=0;j<jmax+nSig;j++){
      if(hist1[j]->GetBinContent(i) >= 0) outf<<hist1[j]->GetBinContent(i)<<" ";
      else outf<<"0 ";
    }
    outf<<endl<<"------------"<<endl;

    for(int j1=0;j1<jmax+nSig;j1++){
      outf<<getfname(f[j1]->GetName())<<"bin"<<i<<" lnN ";
      for(int j2=0;j2<jmax+nSig;j2++){
	if(j1==j2){
	  //	  outf<<"1.20 ";
	  if((hist1[j2]->GetBinContent(i))>0.00000001) outf<<1.2 + (hist1[j2]->GetBinError(i))/(hist1[j2]->GetBinContent(i))<<" ";
	  else outf<<"1.20 ";
	}
	else{
	  outf<<"- ";
	}
      }
      outf<<endl;
    }
    outf.close();
    //delete hist2;
  }
}


void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);

}

