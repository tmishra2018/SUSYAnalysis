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

void makeDatacard_SBinsV7_v2(int mGl,int mNLSP,TString sigFile){
  char name[100];
  const int jmax=5;//no. of backgrounds
  int nSig=1;//1 signal
  //  int nFiles=nSig+1;
  bool do_allSigSyst = 1;
  TString rateSigFile, bTagUncFile;
  string sampleName;
  if(sigFile.Contains("5bbbbZg"))     {
    rateSigFile = sigFile+"FastSim_T5bbbbZg_bTagSF_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    bTagUncFile = sigFile+"FastSim_T5bbbbZg_bTagSFup_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    sampleName = "T5bbbbZg";
  }
  else if(sigFile.Contains("5ttttZg"))     {
    rateSigFile = sigFile+"FastSim_T5ttttZg_bTagSF_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    bTagUncFile = sigFile+"FastSim_T5ttttZg_bTagSFup_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    sampleName = "T5ttttZg";
  }
  else if(sigFile.Contains("5qqqqHg"))     {
    rateSigFile = sigFile+"FastSim_T5qqqqHg_bTagSF_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    bTagUncFile = sigFile+"FastSim_T5qqqqHg_bTagSFup_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    sampleName = "T5qqqqHg";
  }
  else if(sigFile.Contains("6ttZg"))     {
    rateSigFile = sigFile+"FastSim_T6ttZg_bTagSF_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    bTagUncFile = sigFile+"FastSim_T6ttZg_bTagSFup_"+to_string(mGl)+"_"+to_string(mNLSP)+".root"; 
    sampleName = "T6ttZg";
  }
  else {cout<<"Cannot find sig rate file"<<endl;return;}
  sampleName = sampleName+"_"+to_string(mGl)+"_"+to_string(mNLSP);

  double sigXsec = 0.0, sigXsecUnc = 0.0;
  TFile *fXsec = TFile::Open("mGl_Xsecpb_absUnc.root");
  TH1D *hXsec = (TH1D*)fXsec->Get("mGlXsec");
  sigXsec = hXsec->GetBinContent(hXsec->FindBin(mGl));
  sigXsecUnc = hXsec->GetBinError(hXsec->FindBin(mGl));
  cout<<mGl<<" xsec: "<<sigXsec<<" unc: "<<sigXsecUnc<<endl;
  TFile *f[2];
  f[0] = new TFile(rateSigFile);
  f[1] = new TFile("SBinHists.root");
  
  TString rateHistName[jmax+2]={"AllSBins_v7_CD","AllSBins_v7_LElePred","AllSBins_v7_LMuPred","AllSBins_v7_FRPred","AllSBins_v7_ZGPred","AllSBins_v7_MultiJPred","AllSBins_v7_Obs"};
  //  TString processName[jmax+2]={getfname(sigFile),"LEle","LMuTau","Fake","ZG","MultiJ","dataObs"};
  TString processName[jmax+2]={sampleName,"LEle","LMuTau","Fake","ZG","MultiJ","dataObs"};
  TH1D *h_rate[jmax+2];
  for(int i=0;i<jmax+2;i++){
    if(i==0) h_rate[i] = (TH1D*)f[0]->Get(rateHistName[i]); //h_rate[i]->Draw();}
    else{
      h_rate[i] = (TH1D*)f[1]->Get(rateHistName[i]);
    }
    //h_rate[i]->Draw("same");
  }
  
  //***************** setup signal syst ****************
  TString systFile;
  TFile *f_sigSystA, *f_sigSystB;
  TH1D *h_sigSystAcentral, *h_sigJECup, *h_sigJERup, *h_sigISRup, *h_sigGenMET, *h_sigBtagUp;
  if(do_allSigSyst){
    if(sigFile.Contains("5bbbbZg")) systFile = sigFile+"FastSim_T5bbbbZg_SystA_"+to_string(mGl)+"_"+to_string(mNLSP)+".root";
    else if(sigFile.Contains("5ttttZg")) systFile = sigFile+"FastSim_T5ttttZg_SystA_"+to_string(mGl)+"_"+to_string(mNLSP)+".root";
    else if(sigFile.Contains("5qqqqHg")) systFile = sigFile+"FastSim_T5qqqqHg_SystA_"+to_string(mGl)+"_"+to_string(mNLSP)+".root";
    else if(sigFile.Contains("6ttZg")) systFile = sigFile+"FastSim_T6ttZg_SystA_"+to_string(mGl)+"_"+to_string(mNLSP)+".root";
    else {cout<<"Cannot find sig syst file"<<endl;return;}
    f_sigSystA = new TFile(systFile);
    f_sigSystB = new TFile(bTagUncFile);

    h_sigSystAcentral = (TH1D*)f_sigSystA->Get("AllSBins_v7_CD");

    h_sigJECup = (TH1D*)f_sigSystA->Get("AllSBins_v7_JECup_CD");
    h_sigJERup = (TH1D*)f_sigSystA->Get("AllSBins_v7_JERup_CD");
    h_sigGenMET = (TH1D*)f_sigSystA->Get("AllSBins_v7_genMET_CD");
    h_sigISRup = (TH1D*)f_sigSystA->Get("AllSBins_v7_ISRup_CD");

    h_sigBtagUp = (TH1D*)f_sigSystB->Get("AllSBins_v7_CD");

    h_sigJECup->Add(h_sigSystAcentral,-1); h_sigJECup->Divide(h_sigSystAcentral);
    h_sigJERup->Add(h_sigSystAcentral,-1); h_sigJERup->Divide(h_sigSystAcentral);
    h_sigISRup->Add(h_sigSystAcentral,-1); h_sigISRup->Divide(h_sigSystAcentral);
    
    h_sigGenMET->Add(h_sigSystAcentral,-1);
    h_sigGenMET->Scale(0.5);//unc is 0.5 x diff(GenMET,MET)
    h_sigGenMET->Divide(h_sigSystAcentral);

    h_sigBtagUp->Add(h_rate[0],-1);
    h_sigBtagUp->Divide(h_rate[0]);

    cout<<f_sigSystA->GetName()<<endl;
    cout<<f_sigSystB->GetName()<<endl;
    cout<<h_sigJECup->Integral()<<" "<<h_sigJERup->Integral()<<" "<<h_sigGenMET->Integral()<<" "<<h_sigISRup->Integral()<<" "<<h_sigBtagUp->Integral()<<endl;
  }
  //****************************************************
  //for lost ele
  TH1D *h_LEleCS = (TH1D*)f[1]->Get("AllSBins_v7_LEleCS");
  TH1D *h_LEleTFs = (TH1D*)f[1]->Get("LEleTFs");
  TH1D *h_LEleDr = (TH1D*)f[1]->Get("AllSBins_v7_UncTF_dRElePho");
  TH1D *h_LElePDF = (TH1D*)f[1]->Get("AllSBins_v7_UncLEleTF_pdf");
  //  TH1D *h_LElebTag = (TH1D*)f[1]->Get("AllSBins_v7_UncLEleTF_bTag");
  TH1D *h_LEleLepSF = (TH1D*)f[1]->Get("AllSBins_v7_UncLEleTF_LepSF");
  TH1D *h_LEleJEC = (TH1D*)f[1]->Get("AllSBins_v7_UncLEleTF_JEC");

  //for lost mu+had tau
  TH1D *h_LMuCS = (TH1D*)f[1]->Get("AllSBins_v7_LMuCS");
  TH1D *h_LMuTauTFs = (TH1D*)f[1]->Get("LMuTauTFs");
  TH1D *h_LMuPDF = (TH1D*)f[1]->Get("AllSBins_v7_UncLMuTF_pdf");
  //  TH1D *h_LMubTag = (TH1D*)f[1]->Get("AllSBins_v7_UncLMuTF_bTag");
  TH1D *h_LMuLepSF = (TH1D*)f[1]->Get("AllSBins_v7_UncLMuTF_LepSF");

  //for fake
  TH1D *h_FR_StatTrig = (TH1D*)f[1]->Get("AllSBins_v7_FR_StatTrig");
  TH1D *h_FR = (TH1D*)f[1]->Get("FRs");
  TH1D *h_FRSF = (TH1D*)f[1]->Get("AllSBins_v7_UncFR_SF");
  TH1D *h_FRPU = (TH1D*)f[1]->Get("AllSBins_v7_UncFR_PU");
  TH1D *h_FRISR = (TH1D*)f[1]->Get("AllSBins_v7_UncFR_ISRWt");

  //for ZG
  TH1D *h_ZGCS = (TH1D*)f[1]->Get("AllSBins_v7_ZGCS");
  TH1D *h_PureOverLLG_MC = (TH1D*)f[1]->Get("PureOverLLG_MC");
  TH1D *h_dataLLG = (TH1D*)f[1]->Get("dataLLG");
  TH1D *h_ZGTF = (TH1D*)f[1]->Get("AllSBins_v7_ZGTF");
  TH1D *h_highOrd = (TH1D*)f[1]->Get("ZG_HighOrd");

  //for MultiJ
  TH1D *h_MultiJCSraw = (TH1D*)f[1]->Get("AllSBins_v7_rawMultiJCS");
  TH1D *h_doubleR = (TH1D*)f[1]->Get("doubleRatio_MC");
  TH1D *h_HLRatio = (TH1D*)f[1]->Get("AllSBins_v7_HLRatio");
  TH1D *h_MultiJPurity = (TH1D*)f[1]->Get("MultiJLDPpurity");
  //total BG
  TH1D *h_totBG = (TH1D*)f[1]->Get("AllSBins_v7_TotalBG");

  int imax=h_rate[0]->GetNbinsX();
  ofstream outf;
  for(int i=1;i<=imax;i++){
    if( i==1 || i==7 || i==12 || i==17 || i==22 || i==27 ) continue;
    string bTagCorr, njbjCorr,nJCorr,metCorr;
    if(i<=16) bTagCorr = "A";
    else bTagCorr = "B";
    if(i<7){ njbjCorr = "R1"; nJCorr = "LJ";}
    else if(i<12){ njbjCorr = "R2"; nJCorr = "MJ";}
    else if(i<17){ njbjCorr = "R3"; nJCorr = "HJ";}
    else if(i<22){ njbjCorr = "R4"; nJCorr = "LJ";}
    else if(i<27){ njbjCorr = "R5"; nJCorr = "MJ";}
    else{          njbjCorr = "R6"; nJCorr = "HJ";}
    if     (i==2 || i==8  || i==13 || i==18 || i==23 || i==28) metCorr = "M2";
    else if(i==3 || i==9  || i==14 || i==19 || i==24 || i==29) metCorr = "M3";
    else if(i==4 || i==10 || i==15 || i==20 || i==25 || i==30) metCorr = "M4";
    else if(i==5 || i==11 || i==16 || i==21 || i==26 || i==31) metCorr = "M5";
    else if(i==6)                                              metCorr = "M7";

    //    string name2="dataCards/"+getfname(f[0]->GetName())+"_"+"bin"+to_string(i)+".txt";
    string name2="dataCards/card_"+sampleName+"_"+"bin"+to_string(i)+".txt";
    //    cout<<name2;
    sprintf(name,"%s",name2.c_str());
    cout<<name<<endl;
    outf.open(name);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      "# Datacard for mGl= "<<mGl<<" mNLSP= "<<mNLSP<<" with model "<<getfname(f[0]->GetName())<<endl<<
      "# - - - - - - - - - - - - - - - - - - - "<<endl<<
      "imax 1 number of channels"<<endl<<
      "jmax *  number of backgrounds('*' = automatic)"<<endl<<
      "kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl<<
      "------------"<<endl<<
      "bin "<<i<<endl<<
      "observation "<<h_rate[jmax+2-1]->GetBinContent(i)<<endl<<//this is for data
      //"observation "<<h_totBG->GetBinContent(i)<<endl<<//this is for exp limit calc using BG predictions only, no data
      "------------"<<endl<<
      "bin ";
    for(int j=0;j<jmax+nSig;j++){outf<<i<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<processName[j]<<"_bin"<<i<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<" ";}
    outf<<endl<<
      "rate ";
    for(int j=0;j<jmax+nSig;j++){
      if(h_rate[j]->GetBinContent(i) >= 0){
	if(j==0) outf<<0.98*h_rate[j]->GetBinContent(i)<<" ";//98% trigger effieciency. apply for signal only.
	else outf<<h_rate[j]->GetBinContent(i)<<" ";
      }
      else{cout<<"!!!! found -ve events in bin "<<i<<" of hist "<<rateHistName[j]<<". setting as 0 events"<<endl; outf<<"0 ";}
    }
    outf<<endl<<"------------"<<endl;
    
    outf<<"lumi lnN              1.025      -       -       -       -       -"<<endl;
    outf<<"sigStat_b"<<i<<" lnN       "<<1+((h_rate[0]->GetBinError(i))/(h_rate[0]->GetBinContent(i)))<<"    -       -       -       -       -"<<endl;
    outf<<"trigEff_b"<<i<<" lnN       1.02    -       -       -       -       -"<<endl;
    //    outf<<"SigXsec"<<" lnN         "<<1+sigXsecUnc/sigXsec<<"    -       -       -       -       -"<<endl;

    //--------------- Lost Ele-----------------
    outf<<"LEle_b"<<i<<" gmN "<<h_LEleCS->GetBinContent(i)<<"        -     "<<h_LEleTFs->GetBinContent(i)<<"     -       -       -       -"<<endl;
    outf<<"LEleTFUnc_b"<<i<<" lnN     -     "<<1+h_LEleTFs->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    outf<<"LEleDR_b"<<" lnN     -     "<<1+h_LEleDr->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    //    outf<<"LLeptPDF_b"<<" lnN     -     "<<1+h_LElePDF->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    // if(i<=16)
    //   outf<<"LElebTag_b"<<bTagCorr<<" lnN     -     "<<1+h_LElebTag->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    // else       
    //   outf<<"LElebTag_b"<<bTagCorr<<" lnN     -     "<<1-h_LElebTag->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    outf<<"LEleLeptSF_b"<<metCorr<<" lnN     -     "<<1+h_LEleLepSF->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;
    outf<<"LEleJEC_b"<<" lnN     -     "<<1+h_LEleJEC->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"      -       -       -       -"<<endl;

    //--------------- Lost Mu Had Tau-----------------
    outf<<"LMuTau_b"<<i<<" gmN  "<<h_LMuCS->GetBinContent(i)<<"     -       -     "<<h_LMuTauTFs->GetBinContent(i)<<"     -       -       -"<<endl;
    outf<<"LMuTauTFUnc_b"<<i<<" lnN   -       -     "<<1+h_LMuTauTFs->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"      -       -       -"<<endl;
    //    outf<<"LLeptPDF_b"<<" lnN     -     -     "<<1+h_LMuPDF->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"       -       -       -"<<endl;
    // if(i<=16)
    //   outf<<"LMubTag_b"<<bTagCorr<<" lnN     -     -    "<<1+h_LMubTag->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"      -       -       -"<<endl;
    // else       
    //   outf<<"LMubTag_b"<<bTagCorr<<" lnN     -     -    "<<1-h_LMubTag->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"      -       -       -"<<endl;
    outf<<"LMuLeptSF_b"<<metCorr<<" lnN     -     -    "<<1+h_LMuLepSF->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"     -       -       -"<<endl;
    //--------------- PDF unc for lost lepton --------------
    outf<<"LLeptPDF_b"<<" lnN     -     "<<1+h_LElePDF->GetBinError(i)/h_LEleTFs->GetBinContent(i)<<"     "<<1+h_LMuPDF->GetBinError(i)/h_LMuTauTFs->GetBinContent(i)<<"       -       -       -"<<endl;
    //--------------- Fakerate -----------------
    outf<<"FakeUncRate_b"<<i<<" lnN   -       -          -    "<<1+h_FR_StatTrig->GetBinError(i)/h_FR_StatTrig->GetBinContent(i)<<"    -       -"<<endl;
    outf<<"FakeRate_SF_b"<<" lnN   -       -          -    "<<1+h_FRSF->GetBinError(i)/h_FR->GetBinContent(i)<<"    -       -"<<endl;
    outf<<"FakeRate_PU_b"<<" lnN   -       -          -    "<<1+h_FRPU->GetBinError(i)/h_FR->GetBinContent(i)<<"    -       -"<<endl;
    //    outf<<"FakeRate_Trig_b"<<i<<" lnN   -       -          -    "<<1+h_FRTrig->GetBinError(i)/h_FR->GetBinContent(i)<<"    -       -"<<endl;
    outf<<"FakeRate_ISR_b"<<" lnN   -       -          -    "<<1+h_FRISR->GetBinError(i)/h_FR->GetBinContent(i)<<"    -       -"<<endl;
    //--------------- ZG -----------------------
    outf<<"ZG_mcStat_b"<<i<<" lnN     -       -          -      -     "<<1+(h_ZGCS->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
    outf<<"ZGTF_b"<<bTagCorr<<" gmN "<<h_dataLLG->GetBinContent(i)<<"         -       -          -      -     "<<h_ZGTF->GetBinContent(i)<<"       -"<<endl;
    outf<<"ZGPureOverLLG_b"<<bTagCorr<<" lnN  -       -          -      -     "<<1+(h_PureOverLLG_MC->GetBinError(i)/h_PureOverLLG_MC->GetBinContent(i))<<"     -"<<endl;
    if(i<=16)
      outf<<"ZG_highOrdCorr_b"<<nJCorr<<" lnN  -       -          -      -     "<<1+(h_highOrd->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
    else
      outf<<"ZG_highOrdCorr_b"<<nJCorr<<" lnN  -       -          -      -     "<<1-(h_highOrd->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
    //------------------ MultiJ -------------------
    if(h_MultiJPurity->GetBinContent(i) > 0.0)
      outf<<"MultiJ_b"<<i<<" gmN "<<h_MultiJCSraw->GetBinContent(i)<<"     -       -          -          -         -       "<<h_doubleR->GetBinContent(i)*h_HLRatio->GetBinContent(i)*h_MultiJPurity->GetBinContent(i)<<endl;
    else outf<<"MultiJ_b"<<i<<" gmN "<<h_MultiJCSraw->GetBinContent(i)<<"     -       -          -          -         -       "<<h_doubleR->GetBinContent(i)*h_HLRatio->GetBinContent(i)<<endl;
    outf<<"MultiJHL_b"<<njbjCorr<<" lnN "<<"     -       -          -          -         -       "<<1+(h_HLRatio->GetBinError(i)/h_HLRatio->GetBinContent(i))<<endl;
    outf<<"MultiJdR_b"<<njbjCorr<<" lnN "<<"     -       -          -          -         -       "<<1+(h_doubleR->GetBinError(i)/h_doubleR->GetBinContent(i))<<endl;
    if(h_MultiJPurity->GetBinContent(i) < 0.0001)
      outf<<"MultiJEW_b"<<i<<"  lnN "<<"     -       -          -          -         -        0.5/1.0"<<endl;
    else{
      double uncUp = h_MultiJPurity->GetBinContent(i)+h_MultiJPurity->GetBinError(i);
      double uncDn = h_MultiJPurity->GetBinContent(i)-h_MultiJPurity->GetBinError(i);
      if(uncUp < 1.0 && uncDn > 0.) 
	outf<<"MultiJEW_b"<<i<<"  lnN "<<"     -       -          -          -         -        "<<uncDn<<"/"<<uncUp<<endl;
      if(uncUp > 1.0) uncUp = uncUp - 1;
      if(uncDn < 0.0) uncDn = h_MultiJPurity->GetBinContent(i);
    }
    outf.close();
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

