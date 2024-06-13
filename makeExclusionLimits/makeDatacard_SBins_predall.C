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

void makeDatacard_SBins_predall(double mGl,double mNLSP,TString sigFile,char histname1[1000],char histname[1000]){
  char name[1000];
  int jmax=4;//no. of backgrounds
  int nSig=1;//1 signal
  int ndata=1;
  int nFiles=nSig+jmax+ndata;
  TFile *f[nFiles] , *flow[nFiles], *fR[jmax],*fRlow[jmax];
  f[0] = new TFile(sigFile);
  f[1] = new TFile("LL.root");
  f[2] = new TFile("FR.root");//ZGJetsToNuNuG_v17.root  
  f[3] = new TFile("Zin.root");
  f[4] = new TFile("QCD.root");
  f[jmax+1]= new TFile("Data.root");
  TString sigFile2=getfname(f[0]->GetName())+"_lowdphi.root";
  cout<<sigFile2<<endl;
  flow[0] = new TFile(sigFile2);
  flow[1] = new TFile("EWK_LL.root");                                                                                           
  flow[2] = new TFile("EWK_FR.root");                                                                                           
  flow[3] = new TFile("EWK_Zin.root");                                                                                          
  flow[4] = new TFile("QCD.root");
  flow[jmax+1] = new TFile("Data_lowdphi.root");
  fR[0]=new TFile("LL.root");
  fR[1]=new TFile("FR.root");
  fR[2]=new TFile("Zin.root");
  fR[3]=new TFile("QCD.root");

  fRlow[0]=new TFile("EWK_LL.root");
  fRlow[1]=new TFile("EWK_FR.root");
  fRlow[2]=new TFile("EWK_Zin.root");
  fRlow[3]=new TFile("QCD.root");

  TH1D *hist1[nFiles] , *hist2[nFiles] ,*hist3[nFiles], *hist4[nFiles], *Fr[jmax],*Frlow[jmax];
  //  double min_binLowedge=99.9999;
  ofstream outf;
  vector<int> observation;
  for(int i=0;i<(nSig);i++){
    //Pred in high dphi region
    hist1[i] = (TH1D*)f[i]->FindObjectAny(histname);
    //Pred in low dphi region
    hist2[i] = (TH1D*)flow[i]->FindObjectAny(histname);

  cout<<histname<<endl;
    //setLastBinAsOverFlow(hist1[i]);
  }

  //Data CR in LL & FR
  for(int i=1;i<=(jmax-1);i++){
    hist1[i] = (TH1D*)f[i]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");
    hist2[i] = (TH1D*)flow[i]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");
    hist3[i] = (TH1D*)f[i]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1");
    hist4[i] = (TH1D*)flow[i]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1");
    Fr[i-1] = (TH1D*)fR[i-1]->FindObjectAny("fr");
    Frlow[i-1] = (TH1D*)fRlow[i-1]->FindObjectAny("fr");
    cout<<" high dphi file : "<<f[i]->GetName()<<" , low dphi file : "<<flow[i]->GetName()<<endl;
  }
  hist1[jmax] = (TH1D*)f[jmax]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");
  hist2[jmax] = (TH1D*)flow[jmax]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1");
  hist3[jmax] = (TH1D*)f[jmax]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");
  hist4[jmax] = (TH1D*)flow[jmax]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1");
  Fr[jmax-1] = (TH1D*)fR[jmax-1]->FindObjectAny("fr");
  Frlow[jmax-1] = (TH1D*)fRlow[jmax-1]->FindObjectAny("fr");

  //Data in high dphi region     
  hist1[nFiles-1] = (TH1D*)f[nFiles-1]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");
  //Data in low dphi region   
  hist2[nFiles-1] = (TH1D*)flow[nFiles-1]->FindObjectAny("AllSBins_v6_CD_EW_50bin_AB");

  TH1D *hist=(TH1D*)f[0]->FindObjectAny(histname);
  int imax=hist->GetNbinsX();
  //int imax=46;
  for(int i=1;i<=imax;i++){observation.push_back(1);}
    for(int i=1;i<=imax;i++){
      cout<<"MEt bin --> "<<i<<" bin content --> "<<hist->GetBinContent(i)<<endl;
    }
    // cout<<"-----"<<getfname(f[j]->GetName())<<"-----"<<endl;}
  cout<<"file name : "<<getfname(f[0]->GetName())<<endl;

  for(int i=1;i<=imax;i++){
    int q=0;
    if(i==2 || i==9 || i==15 || i==20 || i==25 || i==30 || i==35 || i==41)
      {q=-1;
	f[jmax+1]= new TFile("Data_lowdphi.root");
	hist1[nFiles-1] = (TH1D*)f[nFiles-1]->FindObjectAny("AllSBins_v6_CD_EW_50bin_CD");
      }   
    else
      {
        f[jmax+1]= new TFile("Data.root");
        hist1[nFiles-1] = (TH1D*)f[nFiles-1]->FindObjectAny("AllSBins_v6_CD_EW_50bin_elec1_closure");

      }
    // string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname1+"/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    //    cout<<name2<<endl;
    sprintf(name,"%s",name2.c_str());
    
    outf.open(name);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      "# Datacard for mGl= "<<mGl<<" mNLSP= "<<mNLSP<<endl<<
      "# - - - - - - - - - - - - - - - - - - - "<<endl<<
      "imax * number of bins"<<endl<<
      "jmax *  number of backgrounds('*' = automatic)"<<endl<<
      "kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl<<
      "------------"<<endl<<
      "bin \t \t Bin"<<i<<" \t bin"<<i<<" \t "<<endl<<
      "observation \t"<<hist1[nFiles-1]->GetBinContent(i)<<"\t"<<hist2[nFiles-1]->GetBinContent(i)<<endl<<
      "------------"<<endl<<
      "bin \t \t \t";

    for(int j=0;j<jmax+nSig;j++){outf<<"Bin"<<i<<"\t";}
    for(int j=0;j<jmax+nSig;j++){outf<<"bin"<<i<<"\t";}

    outf<<endl<<
      "process \t \t";
    outf<<"sig \t";
    for(int j=1;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<"\t";}
    outf<<"sig \t";
    for(int j=1;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<"\t";}
    outf<<endl<<
      "process \t \t";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<"\t";}
    for(int j=0;j<jmax+nSig;j++){outf<<j<<"\t";}
    outf<<endl<<
      "rate \t \t";
    for(int j=0;j<jmax+nSig;j++){
      if(hist1[j]->GetBinContent(i) >= 0) outf<<hist1[j]->GetBinContent(i)<<"\t";
      else outf<<"0 ";
    }
    for(int j=0;j<jmax+nSig;j++){
      if(hist2[j]->GetBinContent(i) >= 0) outf<<hist2[j]->GetBinContent(i)<<"\t";
      else outf<<"0 ";

    }
    
    outf<<endl<<"------------"<<endl;
   
 
   for(int j1=1;j1<=jmax-2;j1++){
      outf<<getfname(f[j1]->GetName())<<"_Bin"<<i<<" gmN "<<hist3[j1]->GetBinContent(i)<<" - ";
      for(int j2=1;j2<=jmax-2;j2++){
        if(j1==j2){
	  outf<<Fr[j2-1]->GetBinContent(i)<<" ";
	}
	else outf<<" - ";
      }
      outf<<" - - - - - - - ";
      outf<<endl;
      outf<<getfname(f[j1]->GetName())<<"_bin"<<i<<" gmN "<<hist4[j1]->GetBinContent(i)<<" - - - - - - ";

      for(int j2=1;j2<=jmax-2;j2++){
	if(j1==j2) outf<<Frlow[j2-1]->GetBinContent(i)<<" ";
	else
	  outf<<" - ";
      }
      outf<<" - - "; 
      outf<<endl;
    }
   if(q==-1|| q==0)
     //if(q==2)
   {
     if(i>=2 && i<=19 ) outf<<getfname(f[3]->GetName())<<"_0b"<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - "<<Fr[2]->GetBinContent(i)<<" - - - - "<<Frlow[2]->GetBinContent(i)<<" - "<<endl;
     else if(i>=20 && i<=34) outf<<getfname(f[3]->GetName())<<"_1b"<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - "<<Fr[2]->GetBinContent(i)<<" - - - - "<<Frlow[2]->GetBinContent(i)<<" - "<<endl;
     else if(i>=35 || i<=40) outf<<getfname(f[3]->GetName())<<"_Wtag"<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - "<<Fr[2]->GetBinContent(i)<<" - - - - "<<Frlow[2]->GetBinContent(i)<<" - "<<endl;
     else if(i>=41 || i<=46) outf<<getfname(f[3]->GetName())<<"_Htag"<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - "<<Fr[2]->GetBinContent(i)<<" - - - - "<<Frlow[2]->GetBinContent(i)<<" - "<<endl;
   }
   //if(q==-1) outf<<getfname(f[3]->GetName())<<"_Bin"<<i<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - "<<Fr[2]->GetBinContent(i)<<" - - - - - - "<<endl;
   //else if(q==0) outf<<getfname(f[3]->GetName())<<"_Bin"<<i<<" gmN "<<hist3[3]->GetBinContent(i)<<" - - - - - - - - "<<Frlow[2]->GetBinContent(i)<<" - "<<endl;


    for(int j1=1;j1<=jmax;j1++){
      outf<<"st"<<getfname(f[j1]->GetName())<<"_Bin"<<i<<" lnN "<<" - ";
      for(int j2=1;j2<=jmax;j2++){
        if(j1==j2){
	  //          outf<<Fr[j2-1]->GetBinError(i)<<" ";
          if((Fr[j2-1]->GetBinContent(i))>0.00000001) outf<<1.0 + (Fr[j2-1]->GetBinError(i))/(Fr[j2-1]->GetBinContent(i))<<" ";
          else outf<<"1.0 ";
         }
        else outf<<" - ";
      }
      //      if(j1==4) outf<<" 1.05 - - - - - ";
      // else outf<<" - - - - - - ";
      outf<<" - - - - - ";
      outf<<endl;
      outf<<"st"<<getfname(f[j1]->GetName())<<"_bin"<<i<<" lnN "<<" - - - - - - ";
      for(int j2=1;j2<=jmax;j2++){
        if(j1==j2) {
	  if((Frlow[j2-1]->GetBinContent(i))>0.00000001) outf<<1.0 + (Frlow[j2-1]->GetBinError(i))/(Frlow[j2-1]->GetBinContent(i))<<" ";
	  else outf<<"1.0 ";}
        else
          outf<<" - ";
      }
      //      if(j1==4)outf<<" 1.05 ";
      // else outf<<" - ";
      outf<<endl;

    }
    outf<<endl<<"------------"<<endl;
    int p=1;
    if(p==1)
      {   
      if(i==2 || i==9 || i==15 || i==20 || i==25 || i==30 || i==35 || i==41)
	{
	  outf<<"rateparam"<<p<<"_CbyA   rateParam   Bin"<<i<<"   QCD 1.0"<<endl;
	  outf<<"rateparam"<<p<<"_CbyA   rateParam   bin"<<i<<"   QCD 1.0"<<endl;
	  //	p=i;
	}
      else{
	
       	outf<<"rateparam_Bin"<<i<<"   rateParam   Bin"<<i<<"   QCD (@0*@1) rateparam"<<p<<"_CbyA,rateparam_bin"<<i<<endl;
	//outf<<"rateparam_Bin"<<i<<"   rateParam   Bin"<<i<<"   QCD 1.0"<<endl;
	outf<<"rateparam_bin"<<i<<"   rateParam   bin"<<i<<"   QCD 1.0"<<endl;
      }
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

