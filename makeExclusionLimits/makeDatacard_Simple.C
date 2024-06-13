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

//void makeDatacard_Simple(double mParent,double mKid,TString sigFile,TString histname){
void makeDatacard_Simple(TString sigFile,TString histname){
  TString name;
  int jmax=5;//no. of backgrounds
  int nSig=1;//1 signal
  int nFiles=nSig+jmax;
  TFile *f[nFiles];
  f[0] = new TFile(sigFile+".root");
  f[1] = new TFile("ST__MC2018.root");
  f[2] = new TFile("Rare_MC2018MC2017MC2016.root");
  f[3] = new TFile("TTJets_MC2018.root");
  f[4] = new TFile("WJetsToLNu_HT_MC2018.root");
  f[5] = new TFile("ZJetsToNuNu_HT_MC2018.root");
  TString dirName="";

  //  TString histname = "METvBin_1Wt1Wm1Ht1Hm"; //"METvBin_0-1Wm1Hm";
  TH1D *hist1[nFiles];
  double min_binLowedge=249.9999;
  double scaleToLumi = 1.0;//36.0/137.0;
  ofstream outf;
  vector<int> observation;
  for(int i=0;i<nFiles;i++){
    f[i]->GetObject(dirName+histname,hist1[i]);
    setLastBinAsOverFlow(hist1[i]);
    if(scaleToLumi < 0.99999 || scaleToLumi > 1.00001){
      cout<<"!!! changing lumi (scaling) to a different value. Scale: "<<scaleToLumi<<endl;
      hist1[i]->Scale(scaleToLumi);
    }
  }

  TH1D *hist;
  f[0]->GetObject(dirName+histname,hist);

  cout<<"making datacards for "<<dirName+histname<<endl;

  int imax=hist->GetNbinsX();
  for(int i=1;i<=imax;i++){observation.push_back(1);}
  for(int i=1;i<=imax;i++){
    if(hist->GetBinLowEdge(i) < min_binLowedge){
      cout<<"skipping bins with low edge "<<hist->GetBinLowEdge(i)<<endl;
      continue;
    }
    
    TString name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    name = name2;
    
    outf.open(name);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      //      "# Datacard for mParent= "<<mParent<<" mKid= "<<mKid<<endl<<
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
    for(int j=0;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<"\t";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<"\t\t";}
    outf<<endl<<
      "rate ";
    for(int j=0;j<jmax+nSig;j++){      
      if(hist1[j]->GetBinContent(i) >= 0) outf<<hist1[j]->GetBinContent(i)<<"\t";
      else outf<<"0 ";
    }
    outf<<endl<<"------------"<<endl;

    for(int j1=0;j1<jmax+nSig;j1++){
      outf<<getfname(f[j1]->GetName())<<"bin"<<i<<" lnN ";
      for(int j2=0;j2<jmax+nSig;j2++){
	if(j1==j2){
	  //	  outf<<"1.20 ";
	  if((hist1[j2]->GetBinContent(i))>0.00000001) outf<<1.0 + (hist1[j2]->GetBinError(i))/(hist1[j2]->GetBinContent(i))<<"\t";
	  else outf<<"1.20\t";
	}
	else{
	  outf<<"-\t";
	}
      }
      outf<<endl;
    }
    outf.close();
    //delete hist2;
    TString exe_cmd = "mv " + name + " .temp.tmp; combineCards.py .temp.tmp > "+name + ";rm .temp.tmp";
    //    cout<<exe_cmd<<endl;
    //    system(exe_cmd);
  }
}


void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = sqrt( (lastBinErr*lastBinErr) + (overflErr*overflErr) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);

}

