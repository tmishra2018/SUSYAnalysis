#define PlotLimits_cxx
#include "PlotLimits.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  PlotLimits ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void PlotLimits::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  int evtSurvived=0;

  TString s_data=data;

  int mGlLimit=0, mNLSPLimit=0;
  TFile *fxsec;//=TFile::Open("mGl_Xsecpb_absUnc.root");
  if(s_data.Contains("T5bbbb")) fxsec = new TFile("T5bbbbZg_MassScan.root");
  else if(s_data.Contains("T5tttt")) fxsec = new TFile("T5ttttZg_MassScan.root");
  else if(s_data.Contains("T5qqqq")) fxsec = new TFile("T5qqqqHg_MassScan.root");
  else if(s_data.Contains("T6tt")) fxsec = new TFile("T6ttZg_MassScan.root");
  else if(s_data.Contains("GGM_M1M3")){ fxsec = new TFile("GGM_M1M3_MassScan.root"); mNLSPLimit=101;}
  else if(s_data.Contains("TChiNG")){ fxsec = new TFile("TChiNG_MassScan.root");}
  else if(s_data.Contains("TChiWG")){ fxsec = new TFile("TChiWG_MassScan.root");}

  TH1D *h1_xsec;
  TH2D *h2_xsec;
  if(s_data.Contains("T6tt"))
    h1_xsec=(TH1D*)fxsec->FindObjectAny("mStopXsec");
  else if(s_data.Contains("T5bbbb") || s_data.Contains("T5tttt") || s_data.Contains("T5qqqq"))
    h1_xsec=(TH1D*)fxsec->FindObjectAny("mGlXsec");
  else if(s_data.Contains("GGM_M1M3"))
    h2_xsec=(TH2D*)fxsec->FindObjectAny("GGM_M1M3_xsec");
  else if(s_data.Contains("TChiNG"))
    h1_xsec=(TH1D*)fxsec->FindObjectAny("mEWeakino_SumN2C1_C1C1Xsec");
  else if(s_data.Contains("TChiWG"))
    h1_xsec=(TH1D*)fxsec->FindObjectAny("mEWeakino_N2C1Xsec");
  else {cout<<"AHHHHHH, could not get xsec for limits"<<endl; return;}

  vector<TString> fNames;
  ifstream filein(inputFileList);
  if(filein.is_open()){
    string line1;
    while(getline(filein,line1))
      fNames.push_back(line1);    
  }
  else cout<<"Could not open file: "<<inputFileList<<endl;
  filein.close();
  cout<<"Using Asymptotic "<<endl;
  if(mGlLimit!=0 || mNLSPLimit!=0) cout<<"!!! Careful, you are cutting on m(SUSY Mother) > "<<mGlLimit<<" m(NLSP) > "<<mNLSPLimit<<endl;
  cout<<"Taking xsec from file "<<fxsec->GetName()<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    // if (k > decade)
    //   cout << 10 * k << " %" <<endl;
    // decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    cout<<fCurrent<<" "<<fNames[fCurrent]<<endl;

    //    if(quantileExpected > 0) continue;
    int m1=mh, m2=roundf((mh-m1)*10000);
    double mGl=m1,mNLSP=round(m2),xsec=-100,xsecUnc=0;
    if(!(s_data.Contains("GGM_M1M3"))){
      if(!(s_data.Contains("TChi"))){
	xsec = h1_xsec->GetBinContent(h1_xsec->FindBin(mGl));
	xsecUnc = h1_xsec->GetBinError(h1_xsec->FindBin(mGl));
      }
      else{
	xsec = h1_xsec->GetBinContent(h1_xsec->FindBin(mNLSP));
	xsecUnc = h1_xsec->GetBinError(h1_xsec->FindBin(mNLSP));
      }
    }
    else{
      xsec = h2_xsec->GetBinContent(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP));
      xsecUnc = h2_xsec->GetBinError(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP));
      if(xsec < 1e-10){
	xsec = h2_xsec->GetBinContent(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)+1);
	xsecUnc = h2_xsec->GetBinError(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)+1);
      }
      if(xsec < 1e-10){
	xsec = h2_xsec->GetBinContent(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)-1);
	xsecUnc = h2_xsec->GetBinError(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)-1);
      }
      if(xsec < 1e-10){
	xsec = h2_xsec->GetBinContent(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)+2);
	xsecUnc = h2_xsec->GetBinError(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)+2);
      }
      if(xsec < 1e-10){
	xsec = h2_xsec->GetBinContent(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)-2);
	xsecUnc = h2_xsec->GetBinError(h2_xsec->GetXaxis()->FindBin(mGl),h2_xsec->GetYaxis()->FindBin(mNLSP)-2);
      }
      if(xsec < 1e-10){
	cout<<"I tried to get xsec. I could not. Tried 2GeV variation in mNLSP. Still did not work. Fix your code now."<<endl;
      }
    }
    //    cout<<mGl<<" "<<mNLSP<<" "<<xsec<<" +/- "<<xsecUnc<<endl;
    if(s_data.Contains("TChiWG") || s_data.Contains("TChiNG")){
      mGl=100;
      if(jentry < 5) cout<<"!! Setting gluino mass to 100 for technical reasons for this EW model. Has no physical reasoning."<<endl;
    }
    if(mGl < mGlLimit) continue;
    if(mNLSP < mNLSPLimit) continue;
    if(limit > 100. || limit < 0.0000001) cout<<"Limit:"<<limit<<" (mGl,mNLSP):"<<mGl<<" "<<mNLSP<<endl;

    if(quantileExpected < 0){
      h2_mGlmNLSP_r->Fill(mGl,mNLSP,limit);
      h2_mGlmNLSP_XsecUL->Fill(mGl,mNLSP,xsec*limit);
      h2_mGlmNLSP_XsecUL_v1->Fill(mGl,mNLSP,xsec*limit);
      h2_mGlmNLSP_r_fb->Fill(mGl,mNLSP,limit);
      h2_mGlmNLSP_XsecUL_fb->Fill(mGl,mNLSP,xsec*limit);

      h2_mGlmNLSP_rUnc->Fill(mGl,mNLSP,limitErr);
      h2_mGlmNLSP_r_XsecUp->Fill(mGl,mNLSP,( limit * (xsec/(xsec+xsecUnc)) ));
      h2_mGlmNLSP_r_XsecDn->Fill(mGl,mNLSP,( limit * (xsec/(xsec-xsecUnc)) ));
      //      cout<<"("<<mNLSP<<" , "<<xsecUnc<<")"<<endl;
      //      cout<<"("<<mNLSP<<" , "<<mNLSP<<"),"<<endl;
    }
    if(abs(quantileExpected - 0.1599999) <= 0.0001){
      h2_mGlmNLSP_16pc->Fill(mGl,mNLSP,limit);
    }
    if(abs(quantileExpected - 0.8399999) <= 0.0001){
      h2_mGlmNLSP_84pc->Fill(mGl,mNLSP,limit);
    }
    if(abs(quantileExpected - 0.5) <= 0.0001){
      h2_mGlmNLSP_median->Fill(mGl,mNLSP,limit);
      h2_mGlmNLSP_XsecULexp->Fill(mGl,mNLSP,xsec*limit);
    }
    if(abs(quantileExpected - 0.025) <= 0.0001){
      h2_mGlmNLSP_2p5pc->Fill(mGl,mNLSP,limit);
    }
    if(abs(quantileExpected - 0.975) <= 0.0001){
      h2_mGlmNLSP_97p5pc->Fill(mGl,mNLSP,limit);
    }
    //    cout<<GluinoMass<<" "<<NLSPMass<<endl;
  }
}
