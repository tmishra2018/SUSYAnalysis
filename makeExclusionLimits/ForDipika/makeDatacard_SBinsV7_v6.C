//Also takes into account signalcontamination
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
void clearSidebandBins(TH1D*);
void setSigSystHists(TString);
double getSigCont(double,double,TString, int);
TFile *f[2];
bool sigCont=1;
TH1D *h_sigbTag, *h_sigGenMET, *h_sigISR, *h_sigJEC, *h_sigRateAvgGenRecoMET, *h_sigUncorr;
void makeDatacard_SBinsV7_v6(double mGl,double mNLSP,TString sigFile){
  cout<<"This code needs other files: 'bTagSFupFile.root', 'genRecoMET.root', 'SignalCont_T5ttttZg_LostMu.root' , 'SignalCont_T6ttZg_LostMu.root'. If these files are absent then comment out gen MET syst and bTag SF syst. Or set useAllSigSyst = false."<<endl;
  char name[100];
  const int jmax=5;//no. of backgrounds
  int nSig=1;//1 signal
  bool useAllSigSyst = 1;
  if(!(sigFile.Contains("T5tttt") || sigFile.Contains("T6tt"))) sigCont=0;
  //  int nFiles=nSig+1;
  double sigXsec = 0.0, sigXsecUnc = 0.0;
  TFile *fXsec = TFile::Open("mGl_Xsecpb_absUnc.root");
  TH1D *hXsec = (TH1D*)fXsec->Get("mGlXsec");
  sigXsec = hXsec->GetBinContent(hXsec->FindBin(mGl));
  sigXsecUnc = hXsec->GetBinError(hXsec->FindBin(mGl));
  cout<<mGl<<" xsec: "<<sigXsec<<" unc: "<<sigXsecUnc<<endl;
  f[0] = new TFile(sigFile);
  f[1] = new TFile("SBinHists.root");

  if(useAllSigSyst) setSigSystHists(sigFile);

  TString rateHistName[jmax+2]={"AllSBins_v7_CD","AllSBins_v7_LElePred","AllSBins_v7_LMuPred","AllSBins_v7_FRPred","AllSBins_v7_ZGPred","AllSBins_v7_MultiJPred","AllSBins_v7_Obs"};
  TString processName[jmax+2]={getfname(sigFile),"LEle","LMuTau","Fake","ZG","MultiJ","dataObs"};
  TH1D *h_rate[jmax+2];

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
  TH1D *h_highOrd = (TH1D*)f[1]->Get("ZG_HighOrd");
  TH1D *h_dataLLG = (TH1D*)f[1]->Get("dataLLG");
  TH1D *h_ZGTFfinal = (TH1D*)f[1]->Get("AllSBins_v7_ZGTFfinal");
  TH1D *h_ZGTFpurity = (TH1D*)f[1]->Get("ZGTFpurity");
  TH1D *h_bTagSFZG = (TH1D*)f[1]->Get("bTagSF_ZG");

  //for MultiJ
  TH1D *h_MultiJCSraw = (TH1D*)f[1]->Get("AllSBins_v7_rawMultiJCS");
  TH1D *h_doubleR = (TH1D*)f[1]->Get("doubleRatio_MC");
  TH1D *h_HLRatio = (TH1D*)f[1]->Get("AllSBins_v7_HLRatio");
  TH1D *h_MultiJPurity = (TH1D*)f[1]->Get("MultiJLDPpurity");
  //total BG
  TH1D *h_totBG = (TH1D*)f[1]->Get("AllSBins_v7_TotalBG");
  for(int i=0;i<jmax+2;i++){
    if(i==0) h_rate[i] = (TH1D*)f[0]->Get(rateHistName[i]); //h_rate[i]->Draw();}
    else{
      h_rate[i] = (TH1D*)f[1]->Get(rateHistName[i]);
    }
    //h_rate[i]->Draw("same");
  }

  int imax=h_rate[0]->GetNbinsX();
  ofstream outf;
  for(int i=1;i<=imax;i++){
    if( i==1 || i==7 || i==12 || i==17 || i==22 || i==27 ) continue;
    string bTagCorr, njbjCorr,nJCorr,metCorr;
    if(i<=16) bTagCorr = "";
    else bTagCorr = "";
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

    string name2="dataCards/"+getfname(f[0]->GetName())+"_"+"bin"+to_string(i)+".txt";
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
	if(j==0) {
	  if(useAllSigSyst) outf<<(getSigCont(mGl,mNLSP,sigFile,i))*0.98*0.5*h_sigRateAvgGenRecoMET->GetBinContent(i)<<" ";//for singal, correct for trig eff and take avg of gen & reco MET
	  else              outf<<0.98*h_rate[j]->GetBinContent(i)<<" ";//for singal, correct for trig eff
	  //outf<<0.98*h_rate[j]->GetBinContent(i)<<" ";//for singal, correct for trig eff
	}
	else outf<<h_rate[j]->GetBinContent(i)<<" ";
	//	cout<<j<<" "<<h_rate[j]->GetBinContent(i)<<endl;
      }
      else{cout<<"!!!! found -ve events in bin "<<i<<" of hist "<<rateHistName[j]<<". setting as 0 events"<<endl; outf<<"0 ";}
    }
    outf<<endl<<"------------"<<endl;
    
    outf<<"sigLumiTrig lnN              1.032      -       -       -       -       -"<<endl;//2.5% lumi and 2% trigger
    if(h_rate[0]->GetBinContent(i) > 0.000001)    outf<<"sigStat_b"<<i<<" lnN       "<<1+((h_rate[0]->GetBinError(i))/(h_rate[0]->GetBinContent(i)))<<"    -       -       -       -       -"<<endl;
    if(useAllSigSyst){
      if(h_sigbTag->GetBinContent(i) > 0. && h_sigbTag->GetBinContent(i) < 10.)
	outf<<"sigbTag_b"  <<bTagCorr<<" lnN     "<<h_sigbTag->GetBinContent(i)  <<"    -       -       -       -       -"<<endl;
      if(h_sigGenMET->GetBinContent(i) > 1. && h_sigGenMET->GetBinContent(i) < 10.)
	outf<<"sigGenMET_b"<<metCorr <<" lnN     "<<h_sigGenMET->GetBinContent(i)<<"    -       -       -       -       -"<<endl;
      if(h_sigISR->GetBinContent(i) > 1. && h_sigISR->GetBinContent(i) < 10.)
	outf<<"sigISR_b"   <<nJCorr  <<" lnN     "<<1.0+(sqrt( pow((1.0-h_sigISR->GetBinContent(i)),2) + 0.0004 ))<<"    -       -       -       -       -"<<endl;//ISR unc & 2% scale unc
      outf<<"sigOther_b" <<i       <<" lnN     1.058    -       -       -       -       -"<<endl;//JEC(5%), JER(2%), IsoTrkVeto(2%), JetID(1%)
    }
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
    outf<<"ZGTF_b"<<" gmN "<<h_dataLLG->GetBinContent(i)<<"         -       -          -      -     "<<h_ZGTFfinal->GetBinContent(i)<<"       -"<<endl;
    outf<<"ZGTFpurity_b"<<" lnN  -       -          -      -     "<<1+(h_ZGTFpurity->GetBinError(i)/h_ZGTFpurity->GetBinContent(i))<<"     -"<<endl;
    outf<<"ZG_highOrdCorr_b"<<njbjCorr<<" lnN  -       -          -      -     "<<1+(h_highOrd->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
    if(i<=16) outf<<"ZG_bTagSF_b"<<bTagCorr<<" lnN  -       -          -      -     "<<1+(h_bTagSFZG->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
    else      outf<<"ZG_bTagSF_b"<<bTagCorr<<" lnN  -       -          -      -     "<<1-(h_bTagSFZG->GetBinError(i)/h_ZGCS->GetBinContent(i))<<"       -"<<endl;
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
      else{
	if(uncUp > 1.0) uncUp = uncUp - 1;
	if(uncDn < 0.0) uncDn = h_MultiJPurity->GetBinContent(i);
	outf<<"MultiJEW_b"<<i<<"  lnN "<<"     -       -          -          -         -        "<<uncDn<<"/"<<uncUp<<endl;
      }

    }
    outf.close();
  }
}

void setSigSystHists(TString sigFile){
  TString fileLoc = "";
  TString sfUpFile = "bTagSFupFile.root";
  double num1=0,num2=0,num3=0,den1=0,den2=0,den3=0;

  TFile *f1 = TFile::Open(sigFile);
  TFile *f2 = TFile::Open(sfUpFile);

  TH1D *h_temp = (TH1D*)f1->Get("AllSBins_v7_CD");
  TH1D *h_rate = (TH1D*)h_temp->Clone("sigRate1");
  clearSidebandBins(h_rate);

  //------------ bTag -------------------
  h_sigbTag = (TH1D*)f2->Get("AllSBins_v7_CD");
  clearSidebandBins(h_sigbTag);
  num1 = h_sigbTag->Integral(1,16);
  num2 = h_sigbTag->Integral(17,31);

  den1 = h_rate->Integral(1,16);
  den2 = h_rate->Integral(17,31);

  //  h_sigbTag->Add(h_rate,-1);
  for(int i=1;i<=h_sigbTag->GetNbinsX();i++){
    if(i<=16)
      h_sigbTag->SetBinContent(i,1.0+(abs(num1-den1)/den1));
    else
      h_sigbTag->SetBinContent(i,1.0-(abs(num2-den2)/den2));
    h_sigbTag->SetBinError(i,0);
  }
  //  for(int i=1;i<=h_sigbTag->GetNbinsX();i++) cout<<h_sigbTag->GetName()<<" "<<h_sigbTag->GetBinContent(i)<<endl;
  //------------ gen MET -------------------
  TH1D *h_met = (TH1D*)f1->Get("METvarBin_CD");
  TH1D *h_gmet = (TH1D*)f1->Get("GenMETvarBin_CD");
  h_sigGenMET = (TH1D*)h_rate->Clone("sigGenMETUnc");
  setLastBinAsOverFlow(h_met);
  setLastBinAsOverFlow(h_gmet);
  vector<float> gmetUnc;
  for(int i=1;i<=h_gmet->GetNbinsX();i++){
    if(h_gmet->GetBinLowEdge(i) < 199.99) continue;
    gmetUnc.push_back(1.0+abs(h_gmet->GetBinContent(i)-h_met->GetBinContent(i))/(h_met->GetBinContent(i)));
  }
  for(int i=1;i<=h_rate->GetNbinsX();i++){
    h_sigGenMET->SetBinContent(i,0.);
    h_sigGenMET->SetBinError(i,0.);
    if     ( i==1 || i==7  || i==12 || i==17 || i==22 || i==27 ) continue;
    else if( i==2 || i==8  || i==13 || i==18 || i==23 || i==28 ) h_sigGenMET->SetBinContent(i,gmetUnc[0]);
    else if( i==3 || i==9  || i==14 || i==19 || i==24 || i==29 ) h_sigGenMET->SetBinContent(i,gmetUnc[1]);
    else if( i==4 || i==10 || i==15 || i==20 || i==25 || i==30 ) h_sigGenMET->SetBinContent(i,gmetUnc[2]);
    else if( i==5 || i==11 || i==16 || i==21 || i==26 || i==31 ) h_sigGenMET->SetBinContent(i,gmetUnc[3]);
    else if( i==6                                              ) h_sigGenMET->SetBinContent(i,gmetUnc[4]);
  }
  TFile *f_genReco = TFile::Open("genRecoMET.root");
  h_sigRateAvgGenRecoMET = (TH1D*)f_genReco->Get("AllSBins_v7_CD");
  //  for(int i=1;i<=h_sigGenMET->GetNbinsX();i++) cout<<h_sigGenMET->GetName()<<" "<<h_sigGenMET->GetBinContent(i)<<endl;
  //------------------ ISR ------------------
  h_sigISR = (TH1D*)h_rate->Clone("sigISRUnc");
  TH1D *h_nJISRwt = (TH1D*)f1->Get("nHadJets_SBin_v7_D");
  TH1D *h_isrWts = (TH1D*)h_nJISRwt->Clone("isrWts");
  h_isrWts->Divide((TH1D*)f1->Get("nHadJets_NoISRWt_SBin_v7_D"));
  //  h_isrWts->Print("all");
  TH1D *h_isrWtUncSq = (TH1D*)f1->Get("nHadJets_ISRUncSq_SBin_v7_D");
  h_isrWtUncSq->Divide(h_nJISRwt);//sqrt(bin content) = ISR Wt unc
 
  for(int i=1;i<=h_sigISR->GetNbinsX();i++){
    h_sigISR->SetBinContent(i,0.);
    h_sigISR->SetBinError(i,0.);
    if( i==1 || i==7 || i==12 || i==17 || i==22 || i==27 ) continue;
    if(i<7)   h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(2.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(2.1)))));
    else if(i<12) h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(5.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(5.1)))));
    else if(i<17) h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(7.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(7.1)))));
    else if(i<22) h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(2.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(2.1)))));
    else if(i<27) h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(5.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(5.1)))));
    else h_sigISR->SetBinContent(i, 1.0+(sqrt(h_isrWtUncSq->GetBinContent(h_isrWtUncSq->FindBin(7.1)))/(h_isrWts->GetBinContent(h_isrWts->FindBin(7.1)))));
  }
  //  for(int i=1;i<=h_sigISR->GetNbinsX();i++) cout<<h_sigISR->GetName()<<" "<<h_sigISR->GetBinContent(i)<<endl;
}
void clearSidebandBins(TH1D* h1){
  for(int i=1;i<=h1->GetNbinsX();i++){
    if( i==1 || i==7 || i==12 || i==17 || i==22 || i==27 ){
      h1->SetBinContent(i,0);
      h1->SetBinError(i,0);
    }    
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

double getSigCont(double mGl,double mNLSP,TString sigFile, int bini){
  if(!sigCont) return 1.0;
  else{
    TFile *fLMu,*fLEle;
    if(sigFile.Contains("T6ttZg")){
      fLMu=TFile::Open("SignalCont_T6ttZg_LostMu.root");
      fLEle=TFile::Open("SignalCont_T6ttZg_LostMu.root");/////////!!!!!!!!!!!!!!!
    }
    else if(sigFile.Contains("T5ttttZg")){
      fLMu=TFile::Open("SignalCont_T5ttttZg_LostMu.root");
      fLEle=TFile::Open("SignalCont_T5ttttZg_LostMu.root");/////////!!!!!!!!!!!!!!!                                                                                     
    }
    else {cout<<"Not appying signal contamination."<<endl; return 1.0;}
    double a=0.;
    TH1D *h_LEleTFs = (TH1D*)f[1]->Get("LEleTFs");
    TH1D *h_LMuTauTFs = (TH1D*)f[1]->Get("LMuTauTFs");
    TH1D *h_sigYield = (TH1D*)f[0]->Get("AllSBins_v7_CD");

    TString name1 = "mGlmNLSP_Mu1_bin_"+to_string(bini);/////!!!!!!!!!!!!!!!!
    TH2D *h21=(TH2D*)fLEle->Get(name1);
    a = (h21->GetBinContent(h21->GetXaxis()->FindBin(mGl),h21->GetYaxis()->FindBin(mNLSP)))*(h_LEleTFs->GetBinContent(h_LEleTFs->FindBin(bini)));//lost ele
    
    name1 = "mGlmNLSP_Mu1_bin_"+to_string(bini);
    h21=(TH2D*)fLMu->Get(name1);
    a += (h21->GetBinContent(h21->GetXaxis()->FindBin(mGl),h21->GetYaxis()->FindBin(mNLSP)))*(h_LMuTauTFs->GetBinContent(h_LMuTauTFs->FindBin(bini)));//lost mu

    if(h_sigYield->GetBinContent(h_sigYield->FindBin(bini)) < 0.00001){
      return 1.0;
    }

    a = (1.0-(a/h_sigYield->GetBinContent(h_sigYield->FindBin(bini))));   
    if(a>1.0 || a< 0.0001) return 1.0;
    if(a<0){cout<<"Warning!!! signal contamination more than signal. Setting contamination to 99%"<<endl;}

    fLMu->Close();
    fLEle->Close();
    cout<<"Contamination factor:"<<a<<" ";
    return a;
  }
}
