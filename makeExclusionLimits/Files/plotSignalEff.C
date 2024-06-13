#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"TGraph2D.h"

TLatex textOnTop,intLumiE;
double intLumi = 35.9;
void plotSignalEff(TString fName){

  // const Int_t NRGBs = 5;
  // const Int_t NCont = 255;
  
  // Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  // Double_t red[NRGBs]   = { 0.50, 0.50, 1.00, 1.00, 1.00 };
  // Double_t green[NRGBs] = { 0.50, 1.00, 1.00, 0.60, 0.50 };
  // Double_t blue[NRGBs]  = { 1.00, 1.00, 0.50, 0.40, 0.50 };

  // TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  // gStyle->SetNumberContours(NCont);

  TFile *f1=new TFile(fName);
  TString modelName=f1->GetName();

  TH2D *h2_num0=(TH2D*)f1->Get("EventsPassed");
  TH2D *h2_den0=(TH2D*)f1->Get("MGlMNLSP");
  TH2D *h2_num = new TH2D("num","num",120,12.5,3012.5,120,12.5,3012.5);
  TH2D *h2_den = new TH2D("den","den",120,12.5,3012.5,120,12.5,3012.5);
  for(int i=0;i<=h2_num0->GetNbinsX();i++){
    for(int j=0;j<=h2_num0->GetNbinsY();j++){
      if(h2_num0->GetBinContent(i,j) < 0.000001) continue;
      h2_num->Fill(h2_num0->GetXaxis()->GetBinCenter(i),h2_num0->GetYaxis()->GetBinCenter(j),h2_num0->GetBinContent(i,j));
      h2_den->Fill(h2_den0->GetXaxis()->GetBinCenter(i),h2_den0->GetYaxis()->GetBinCenter(j),h2_den0->GetBinContent(i,j));
    }
  }
  
  // h2_num->RebinX(10); h2_num->RebinY(10);
  // h2_den->RebinX(10); h2_den->RebinY(10);

  h2_num->Divide(h2_den);

  gStyle->SetOptStat(0);
  TCanvas *c3=new TCanvas("c3","c3",1500,900);
  //c3->SetBottomMargin(0.12);
  // c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.15);
  h2_num->GetXaxis()->SetTitleSize(0.05);
  h2_num->GetXaxis()->SetLabelSize(0.05);
  h2_num->GetYaxis()->SetTitleSize(0.05);
  h2_num->GetYaxis()->SetLabelSize(0.05);
  //  h2_num->SetMinimum(0.10);  h2_num->SetMaximum(0.45);//T5bbbbZg,T5qqqqHg
  //  h2_num->SetMinimum(0.);  h2_num->SetMaximum(0.18);//T5tttt
  // c3->SetLogz();

  TPaveText *decayMode = new TPaveText(0.15,0.75, 0.5,0.85,"NDC");
  decayMode->SetShadowColor(0);   decayMode->SetFillColor(0);
  decayMode->SetBorderSize(0);
  if(modelName.Contains("T5qqqqHg")){
    decayMode->AddText("#bf{pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow q#bar{q}#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0}#rightarrow #gamma/H#tilde{G}}");
    modelName="T5qqqqHg";
    h2_num->GetXaxis()->SetRangeUser(1500,2600);
    h2_num->GetYaxis()->SetRangeUser(0,2499);
    h2_num->SetMinimum(0.10);  h2_num->SetMaximum(0.45);
  }
  else if(modelName.Contains("T5ttttZg")){
    decayMode->AddText("#bf{pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow t#bar{t}#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0}#rightarrow #gamma/Z#tilde{G}}");
    modelName="T5ttttZg";
    h2_num->GetXaxis()->SetRangeUser(1500,2600);
    h2_num->GetYaxis()->SetRangeUser(0,2499);
    h2_num->SetMinimum(0.);  h2_num->SetMaximum(0.18);
  }
  else if(modelName.Contains("T5bbbbZg")){
    decayMode->AddText("#bf{pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow b#bar{b}#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0}#rightarrow #gamma/Z#tilde{G}}");
    modelName="T5bbbbZg";
    h2_num->GetXaxis()->SetRangeUser(1500,2600);
    h2_num->GetYaxis()->SetRangeUser(0,2499);
    h2_num->SetMinimum(0.10);  h2_num->SetMaximum(0.45);
  }
  else if(modelName.Contains("T6ttZg")){
    decayMode->AddText("#bf{pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0}#rightarrow #gamma/Z#tilde{G}}");
    modelName="T6ttZg";
    h2_num->GetXaxis()->SetRangeUser(600,1550);
    h2_num->GetYaxis()->SetRangeUser(0,1500);
    h2_num->SetMinimum(0.);  h2_num->SetMaximum(0.25);
  }
  else decayMode->AddText(modelName);

  if(modelName!="T6ttZg")  h2_num->SetTitle(";m_{#tilde{g}}(GeV);m_{#tilde{#chi}_{1}^{0}}(GeV)");
  else                     h2_num->SetTitle(";m_{ #tilde{t}}(GeV);m_{#tilde{#chi}_{1}^{0}}(GeV)");
  h2_num->GetZaxis()->SetTitle("Acceptance x Efficiency");
  //  TGraph2D *gr2d_eff = new TGraph2D(h2_num);
  
  h2_num->Draw("colz");
  //  gr2d_eff->Draw("colz");
  decayMode->Draw();

  char name3[100];
  textOnTop.SetTextSize(0.05);
  intLumiE.SetTextSize(0.05);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation}}");
  sprintf(name3,"#bf{%0.1f fb^{-1}(13TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.6,0.91,name3);
  //  c3->SaveAs(modelName+"_AccEff.pdf");
}
