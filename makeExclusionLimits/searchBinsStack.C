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
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=11,nBG=8;    //Specify no. of files
TFile *f[nfiles];
//int col[11]={kPink+1,kTeal+9,kOrange,kYellow,kBlue,kCyan,kGreen,kMagenta+2,kRed,kMagenta,kBlue+2};  //Specify Colors
int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kMagenta+2,kCyan,kBlue,kRed,kMagenta,kBlue+2};  //Specify Colors b's
//char name[100],name2[100];
TCanvas *c_cA=new TCanvas("MET_incl","MET for NJ>=2",1500,850);
TCanvas *c_cB=new TCanvas("MET_SBins","MET for all search bins",1500,850);

void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
TLegend *legend1=new TLegend(0.50, 0.9,  0.7, 0.70);
TLegend *legend2=new TLegend(0.7, 0.9,  0.90, 0.65);
//TLegend *legend2=new TLegend(0.6, 0.90,  0.98, 0.45);
void setLastBinAsOverFlow(TH1D*);
void searchBinsStack(){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);

  f[0] = new TFile("WGJetsToLNuG.root");
  f[1] = new TFile("TTGJets.root");
  f[2] = new TFile("TTJetsHT.root");
  f[3] = new TFile("WJetsToLNu.root");
  f[4] = new TFile("ZJetsToNuNu.root");
  f[5] = new TFile("ZGJetsToNuNuG.root");
  f[6] = new TFile("GJets.root");
  f[7] = new TFile("QCD.root");
  // f[8] = new TFile("HG_NLSP150.root");
  // f[9] = new TFile("HG_NLSP1000.root");
  //f[9] = new TFile("HG_NLSP150.root");
  f[8] = new TFile("T1bbbb_ZG_mGl1600_NLSP150.root");
  f[9] = new TFile("T1bbbb_ZG_mGl1600_NLSP1000.root");
  //  f[10] = new TFile("a.root");
  f[10] = new TFile("T1bbbb_ZG_mGl1600_NLSP1550.root");


  gStyle->SetTextSize(2);
  THStack *hs_MET=new THStack("MET_Stack","MET Stacked");
  THStack *hs_MET_SB=new THStack("MET_SB","MET for all search bins");
  //TH1D *h_R;
  TH1D *h_MET_R[nfiles];
  for(int i=0;i<nfiles;i++){
    sprintf(name,"hist_file%i",i);
    h_MET_R[i]=new TH1D(name,name,21,0.5,21.5);
  }
  vector<double> Bcnt;
  double intLumi=0.0;
  TLatex tl1;
  for(int i=0;i<nfiles;i++){
    TH1D *h_intLumi=(TH1D*)f[i]->FindObjectAny("intLumi");
    if(i==0) {
      intLumi=h_intLumi->GetMean();
      sprintf(name2, "%.2f fb^{-1}",intLumi);
    }
    else{
      if(abs(intLumi-h_intLumi->GetMean())>0.0001)
	cout<<"Integarted lumi for "<<f[i]->GetName()<<" is "<<h_intLumi->GetMean()<<" and for other files it is different"<<endl;
    }
    
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny("AllSBins");//MET_R1
    decorate(h_MET,i,f[i]->GetName());
    
    if(i<=(nBG-1))  hs_MET->Add(h_MET);
    if(i==nBG-1) {
      c_cA->cd(); 
      hs_MET->Draw("BAR HIST");
      hs_MET->Draw("HIST");
      decorate(hs_MET,i,f[i]->GetName()); 
    }
    if(i>=nBG){ 
      c_cA->cd(); 
      h_MET->Draw("HIST sames");
    }
    drawlegend(h_MET,i,f[i]->GetName());
    if(i==nfiles-1) hs_MET->SetTitle(";Bin No.;Events");

    //----------------------all search bins----------------------------
    for(int j=1;j<=3;j++){
      name3="MET_R"+to_string(j);
      TH1D *h_R=(TH1D*)f[i]->FindObjectAny(name3);
      setLastBinAsOverFlow(h_R);
      //sr_Integral+=h_R->Integral();
      for(int k=6;k<=h_R->GetNbinsX();k++){
    	h_MET_R[i]->SetBinContent( ((j-1)*7)+k-5,h_R->GetBinContent(k) );
    	h_MET_R[i]->SetBinError  ( ((j-1)*7)+k-5,h_R->GetBinError(k) );
      }
    }
    decorate(h_MET_R[i],i,f[i]->GetName());
    
    if(i<=(nBG-1))  hs_MET_SB->Add(h_MET_R[i]);
    if(i==nBG-1) {
      c_cB->cd(); 
      hs_MET_SB->Draw("BAR HIST");
      hs_MET_SB->Draw("HIST");
      decorate(hs_MET_SB,i,f[i]->GetName()); 
    }
    if(i>=nBG){ 
      c_cB->cd(); 
      h_MET_R[i]->Draw("HIST sames");
    }
    if(i==nfiles-1) hs_MET_SB->SetTitle(";Bin No.;Events");
    cout<<f[i]->GetName()<<"\t";
    for(int j=1;j<=h_MET_R[i]->GetNbinsX();j++){
      cout<<h_MET_R[i]->GetBinContent(j)<<"\t";
    }cout<<endl;
  }
  c_cA->cd(); gPad->SetLogy();legend1->Draw();legend2->Draw();
  gPad->RedrawAxis();

  textOnTop.SetTextSize(0.05);
  intLumiE.SetTextSize(0.05);
  textOnTop.DrawLatexNDC(0.1,0.91,"CMS #it{#bf{Simulation}}");
  intLumiE.SetTextSize(0.05);
  sprintf(name2,"#bf{%0.1f fb^{-1}(13TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.73,0.91,name2);

  TPaveText *text1 = new TPaveText(0.55,0.035,7.45,0.07,"NB");
  text1->AddText("0b, NJ:2to4");   text1->Draw();
  TPaveText *text2 = new TPaveText(7.55,0.035,14.45,0.07,"NB"); 
  text2->AddText("0b, NJ:>=5");  text2->Draw();
  TPaveText *text3 = new TPaveText(14.55,0.035,21.45,0.07,"NB"); 
  text3->AddText("1b, NJ:2to4");  text3->Draw();
  TPaveText *text4 = new TPaveText(21.55,0.035,28.45,0.07,"NB"); 
  text4->AddText("1b, NJ:>=5");  text4->Draw();
  TPaveText *text5 = new TPaveText(28.55,0.035,34.45,0.07,"NB"); 
  text5->AddText(">=2b, NJ:>=2");  text5->Draw();
  c_cA->SaveAs("met_inclusive.png");
  //----------------------all search bins---------------------------- 
  c_cB->cd(); gPad->SetLogy();legend1->Draw();legend2->Draw();
  textOnTop.SetTextSize(0.05);
  intLumiE.SetTextSize(0.05);
  textOnTop.DrawLatexNDC(0.1,0.91,"CMS #it{#bf{Simulation}}");
  intLumiE.SetTextSize(0.05);
  sprintf(name2,"#bf{%0.1f fb^{-1}(13TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.73,0.91,name2);
  TPaveText *text6 = new TPaveText(0.55,0.035,7.45,0.07,"NB");
  text6->AddText("NJ: 2 to 4");   text6->Draw();
  TPaveText *text7 = new TPaveText(7.55,0.035,14.45,0.07,"NB"); 
  text7->AddText("NJ: 5 or 6");  text7->Draw();
  TPaveText *text8 = new TPaveText(14.55,0.035,21.45,0.07,"NB"); 
  text8->AddText("NJ: >=7");  text8->Draw();


  c_cB->SaveAs("searchBins.png");
  cout<<"*****************************************************"<<endl;
  cout<<"Int Lumi(inv.fb) for file1:"<<setprecision(4)<<intLumi<<endl;
    
}

void decorate(THStack *hs,int i,const char* fname){
  hs->SetMinimum(0.1);
  hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.04);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  gPad->Update();
  gStyle->SetOptStat(0);
}
void decorate(TH1D* hist,int i,const char* fname){
  hist->SetLineColor(col[i]);
  if(i<nBG) {
    hist->SetFillColor(col[i]);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
  }
  else hist->SetLineWidth(2);
  hist->SetTitle(0);
  hist->GetXaxis()->SetLabelSize(.04);
  hist->GetYaxis()->SetLabelSize(.05);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  // drawlegend(hist,i,fname);
  gPad->Update();
  setLastBinAsOverFlow(hist);
  gStyle->SetOptStat(0);
  //Hlist.Add(hist);
  /*
  TPaveStats *st=(TPaveStats*)hist->FindObject("stats");
  st->SetLineColor(col[i]);
  st->SetTextColor(col[i]);
  if(i>5){
    st->SetX1NDC(0.55);
    st->SetX2NDC(0.75); 
    st->SetY1NDC(0.90-(i-6)*0.1);
    st->SetY2NDC(0.80-(i-6)*0.1);
  }
  else{
    st->SetX1NDC(0.75);
    st->SetX2NDC(0.98); 
    // st->SetY1NDC(0.90-i*0.2);
    //st->SetY2NDC(0.70-i*0.2);
    st->SetY1NDC(0.90-i*0.1);
    st->SetY2NDC(0.80-i*0.1);
  }
  */
}

void drawlegend(TH1D *hist,int i,const char* fname){
  gStyle->SetLegendBorderSize(0);
  //TLegend *legend=new TLegend(0.65,0.90-i*0.05,0.98,0.85-i*0.05);
  // TLegend *legend=new TLegend(0.65,0.90-i*0.2,0.98,0.86-i*0.2);
  //TLegend *legend;
  /*//legend with stat box
  if(i>5){
    legend=new TLegend(0.55,0.9-(i-6)*0.1,0.95,0.85-(i-6)*0.1);//for legend with stat box
  }
  else legend=new TLegend(0.55, 0.90-i*0.05,  0.98, 0.85-i*0.05);
  */
  //legend=new TLegend(0.55, 0.90-i*0.05,  0.98, 0.85-i*0.05);
  
  string lName=fname;
  lName.pop_back();lName.pop_back();lName.pop_back();lName.pop_back();lName.pop_back();
  
  if(lName=="ZJetsToNuNu"){lName="Z( #rightarrow #nu #bar{#nu})+jets";}
  else if(lName=="DYJetsToLL"){lName="DY( #rightarrow l^{+}l^{-})";}
  else if(lName=="WJetsToLNu"){lName="W( #rightarrow l  #nu)+jets";}
  else if(lName=="RareProcess"){}
  else if(lName=="TTJetsHT"){lName="t #bar{t}";}
  else if(lName=="QCD"){lName="QCD";}
  else if(lName=="GJets"){lName="#gamma +jets";}
  else if(lName=="WGJetsToLNuG"){lName="W(#rightarrow l #nu)+ #gamma";}
  else if(lName=="ZGJetsToNuNuG"){lName="Z(#rightarrow #nu #bar{#nu})+ #gamma";}
  else if(lName=="TTGJets"){lName="t #bar{t}+ #gamma";}
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP150"){lName="T5bbbb_ZG_150";}
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP1000"){lName="T5bbbb_ZG_1000";}
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP1550"){lName="T5bbbb_ZG_1550";}

  const char *l_name=lName.c_str();
  if(i<5){
    if(i<nBG)legend1->AddEntry(hist,l_name,"f");
    else legend1->AddEntry(hist,l_name,"l");
    legend1->SetTextSize(0.04);
  }
  else{
    if(i<nBG)legend2->AddEntry(hist,l_name,"f");
    else legend2->AddEntry(hist,l_name,"l");
    legend2->SetTextSize(0.04);
  }
  // legend->SetTextColor(col[i]);
  //legend->SetLineColor(col[i]);
  //if(i==nfiles-1) legend->Draw();
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
