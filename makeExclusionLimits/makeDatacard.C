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
class cls1{
  int imax;//no. of channels
  static const int jmax=6;//no. of backgrounds
  int kmax=5;//number of nuisance parameters (sources of systematical uncertainties)
  static const int nSig=1;//no. of signals
  vector<int> observation;
  char histname[100];
  TFile *f[jmax+nSig];
 
  bool skipBins=false;
  double skipBinsBelow=0,skipBinsAbove=20000;

  void writeToFile();
  string getfname(const char *fname1){string fname=fname1;fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();return fname;}
  double getbincont(int j,int i){
    TH1D *hist2=(TH1D*)f[j]->FindObjectAny(histname);
    if(i==(hist2->GetNbinsX())){return hist2->GetBinContent(i)+hist2->GetBinContent(i+1);}
    else {return hist2->GetBinContent(i);}
  }
public:
  void setInputs();
  double mGl,mNLSP;
  int cutnum=10000;
  ~cls1();
};

void makeDatacard(double mP,double mD,int cutnum){
  cls1 c1;
  c1.mGl=mP;
  c1.mNLSP=mD;
  c1.cutnum=cutnum;
  c1.setInputs();
}

void cls1::setInputs(){ 
  //f[0] = new TFile("NLSP1500.root");//f[0] must be signal. If more than one signal, be careful.Modify the whole code.
  TString inFname[jmax+nSig];char SigFName[10];
  sprintf(SigFName,"HG_NLSP%.0f_Cut",mNLSP);
  inFname[0]=SigFName         +to_string(cutnum)+".root";
  inFname[1]="GJets_Cut"      +to_string(cutnum)+".root";
  inFname[2]="QCD_Cut"        +to_string(cutnum)+".root";
  inFname[3]="TTJets_Cut"     +to_string(cutnum)+".root";
  inFname[4]="TTGJets_Cut"    +to_string(cutnum)+".root";
  inFname[5]="WJetsToLNu_Cut" +to_string(cutnum)+".root";
  inFname[6]="ZJetsToNuNu_Cut"+to_string(cutnum)+".root";
  for(int i=0;i<jmax+nSig;i++){f[i] = new TFile(inFname[i]);}
 
  cout<<"Signal is "<<f[0]->GetName()<<endl;
  sprintf(histname,"METvarBin");skipBins=true;skipBinsAbove=20000;
  //  sprintf(histname,"STvarBin");skipBins=true;skipBinsAbove=20000;
  //sprintf(histname,"nJets");skipBins=true;skipBinsBelow=4;skipBinsAbove=16;//specify the lowerEdge of the bins
  //  sprintf(histname,"nBTags");skipBins=true;skipBinsBelow=0;skipBinsAbove=5;//specify the lowerEdge of the bins
  TH1D *hist=(TH1D*)f[0]->FindObjectAny(histname);
  imax=hist->GetNbinsX();
  
  for(int j=1;j<=imax;j++){  
    double nentries=0;TH1D *hist4;
    for(int i=1;i<jmax+nSig;i++){
      hist4=(TH1D*)f[i]->FindObjectAny(histname);
      nentries+=hist4->GetBinContent(j);
      //      cout<<f[i]->GetName()<<" "<<
    }
    if(nentries<0.000001)skipBinsBelow=hist4->GetBinLowEdge(j+1);
    else break;
  }
  for(int i=0;i<imax;i++){observation.push_back(1);}
  writeToFile();
}

void cls1::writeToFile(){
  ofstream outf;
  char name[100];
  ifstream infile("specifyCuts.txt");
  string line;int lnum=0;bool usecut=false;
  char typeofCut[10]="zzzz";
  //double STmin=0,STmax=20000,METmin=0,METmax=20000,NJmin=0,NJmax=500,PtGmin=100,PtGmax=20000,dphi1cut=0,dphi2cut=0,dphi3cut=0,dphi4cut=0,dphiMET_Gcut=0;
  while (getline(infile, line)){
    std::istringstream iss(line);
    if(lnum!=0){
      if (!(iss >> usecut)) { break; } // error
      int t1;
      iss>>t1;
      if(t1==cutnum) {iss >>  typeofCut;break;}
    }
    lnum++;
  }
  infile.close();
  double statErr[jmax+nSig];
  for(int i=0;i<imax;i++){
    sprintf(name,"dataCards/NLSP%.0f_mGl%.0f_%s_Cut%i_%s_bin%i.txt",mNLSP,mGl,histname,cutnum,typeofCut,i+1);
    if(skipBins){
      TH1D *hist3=(TH1D*)f[0]->FindObjectAny(histname);
      if( ((hist3->GetBinLowEdge(i+1)) < skipBinsBelow) || ((hist3->GetBinLowEdge(i+1)) >= skipBinsAbove) ){
	cout<<"skipping bin with lowEdge "<<(hist3->GetBinLowEdge(i+1))<<" of "<<histname<<".Not creating file "<<name<<endl;
	delete hist3;
	continue;
      }
    }

    for(int p1=0;p1<jmax+nSig;p1++){
      TH1D *hist3=(TH1D*)f[p1]->FindObjectAny(histname);
      if(getbincont(p1,i+1)<=0) statErr[p1]=0.0;
      else statErr[p1]=hist3->GetBinError(i+1)/getbincont(p1,i+1);
      delete hist3;
    }
    outf.open(name);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      "# Datacard for mGl= "<<mGl<<" mNLSP= "<<mNLSP<<endl<<
      "# - - - - - - - - - - - - - - - - - - - "<<endl<<
      "imax 1 number of channels"<<endl<<
      "jmax *  number of backgrounds('*' = automatic)"<<endl<<
      "kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl<<
      "------------"<<endl<<
      "bin "<<histname<<i+1<<"Cut"<<cutnum<<endl<<
      "observation "<<observation[i]<<endl<<
      "------------"<<endl<<
      "bin ";
    for(int j=0;j<jmax+nSig;j++){outf<<histname<<i+1<<"Cut"<<cutnum<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<" ";}
    outf<<endl<<
      "rate ";
    for(int j=0;j<jmax+nSig;j++){
      if(getbincont(j,i+1) >= 0) outf<<getbincont(j,i+1)<<" ";
      else outf<<"0 ";
    }
    outf<<endl<<"------------"<<endl;

    for(int j1=0;j1<jmax+nSig;j1++){
      outf<<getfname(f[j1]->GetName())<<"bin"<<i+1<<" lnN ";
      for(int j2=0;j2<jmax+nSig;j2++){
	if(j1==j2){
	  //	  outf<<"1.20 ";
	  outf<<1.2+statErr[j2]<<" ";
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


cls1::~cls1(){
  for(int i=0;i<jmax+nSig;i++){delete f[i];}
}
