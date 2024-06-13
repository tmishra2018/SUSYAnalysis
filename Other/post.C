#include "../plotting/HttStylesNew.cc"
#include "../plotting/CMS_lumi.C"
#include "TLatex.h"
#include "TFile.h"
#include "THStack.h"


void MakePostFitPlots( TString ptRange = "5to7p5") {

  SetStyle();
  TH1::SetDefaultSumw2();

  TString xtitle = "minv (GeV)";
  TString ytitle = "Events";

  TString inputFileName = "templates_trk_isolation_" + ptRange;
  TString mlfitFileName = "mlfit_"+ptRange;
  TFile * inputs = new TFile(inputFileName+".root");
  if(!inputs){
    cout<<"File "<<inputFileName<<".root not available. Please Check."<<endl;
    exit(-1);      
  }
  TFile * mlfit  = new TFile(mlfitFileName+".root");
  if(!mlfit){
    cout<<"File "<<mlfitFileName<<".root not available. Please Check."<<endl;
    exit(-1);      
  }
  TH1F * h_Data = (TH1F*)inputs->Get("data_obs");
  TH1F * h_QCD  = (TH1F*)inputs->Get("qcd");
  TH1F * h_EWK  = (TH1F*)inputs->Get("ewk");
  TH1F * h_TT   = (TH1F*)inputs->Get("tt");
  TH1F * h_Z    = (TH1F*)inputs->Get("z");
  TH1F * h_ZTT  = (TH1F*)inputs->Get("ztt");

  TH1F * h_QCD_x   = (TH1F*)mlfit->Get("shapes_fit_s/ch1/qcd");
  TH1F * h_EWK_x   = (TH1F*)mlfit->Get("shapes_fit_s/ch1/ewk");
  TH1F * h_TT_x    = (TH1F*)mlfit->Get("shapes_fit_s/ch1/tt");
  TH1F * h_Z_x     = (TH1F*)mlfit->Get("shapes_fit_s/ch1/z");
  TH1F * h_ZTT_x   = (TH1F*)mlfit->Get("shapes_fit_s/ch1/ztt");
  TH1F * h_total_x = (TH1F*)mlfit->Get("shapes_fit_s/ch1/total");

  float WNORM  = h_ZTT_x->Integral();
  float WDATA  = h_Data->Integral() - h_QCD_x->Integral() - h_EWK_x->Integral() - h_TT_x->Integral()- h_Z_x->Integral();
  float wscale = WDATA/WNORM;
  // wscale = 1.;
  float total  =h_QCD_x->Integral()+h_EWK_x->Integral()+h_TT_x->Integral()+h_Z_x->Integral()+wscale*h_ZTT_x->Integral();

  cout << "h_QCD  : " << h_QCD_x->Integral()         << endl;
  cout << "h_EWK  : " << h_EWK_x->Integral()         << endl;
  cout << "h_Z    : " << h_Z_x->Integral()           << endl;
  cout << "h_TT   : " << h_TT_x->Integral()          << endl;
  cout << "h_ZTT  : " << wscale*h_ZTT_x->Integral()  << endl;
  cout << "SUM    : " << total                       << endl;
  cout << "DAT    : " << h_Data->Integral()          << endl;
  cout << "wscale : " << wscale                      << endl;
  
  TH1F * bkgdErr = (TH1F*)h_ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3001);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  

  int nBins = h_Data->GetNbinsX();
  for (int iB=1; iB<=nBins; ++iB) {

    h_QCD -> SetBinError(iB,0);
    h_EWK -> SetBinError(iB,0);
    h_TT  -> SetBinError(iB,0);
    h_Z   -> SetBinError(iB,0);
    h_ZTT -> SetBinError(iB,0);

    h_QCD -> SetBinContent(iB , h_QCD_x->GetBinContent(iB));
    h_EWK -> SetBinContent(iB , h_EWK_x->GetBinContent(iB));
    h_TT  -> SetBinContent(iB , h_TT_x->GetBinContent(iB));
    h_Z   -> SetBinContent(iB , h_Z_x->GetBinContent(iB));
    h_ZTT -> SetBinContent(iB , wscale*h_ZTT_x->GetBinContent(iB));

    double tot     = h_ZTT->GetBinContent(iB) + h_QCD->GetBinContent(iB) + h_EWK->GetBinContent(iB) + h_TT->GetBinContent(iB) + h_Z->GetBinContent(iB);
    double tot_err = h_total_x->GetBinError(iB);
      
    bkgdErr->SetBinContent(iB,tot);
    bkgdErr->SetBinError(iB,tot_err);

  }

  h_ZTT->SetFillColor(TColor::GetColor("#FFCC66"));
  h_QCD->SetFillColor(TColor::GetColor("#FFCCFF"));
  h_EWK->SetFillColor(TColor::GetColor("#DE5A6A"));
  h_TT ->SetFillColor(TColor::GetColor("#9999CC"));
  h_Z  ->SetFillColor(TColor::GetColor("#6F2D35"));

  THStack *stack = new THStack(ptRange,"");
  stack->Add(h_TT);
  stack->Add(h_Z);
  stack->Add(h_EWK);
  stack->Add(h_QCD);
  stack->Add(h_ZTT);

  // Set error of two first bins to zero
  //bkgdErr->SetBinError(1,0);
  //bkgdErr->SetBinError(2,0);
  //bkgdErr->SetBinError(3,0);
  //bkgdErr->SetBinError(4,0);

  h_Data->GetXaxis()->SetTitle("");
  h_Data->GetYaxis()->SetTitle(ytitle);
  h_Data->GetXaxis()->SetLabelSize(0.);
  h_Data->GetYaxis()->SetTitleOffset(1.3);
  h_Data->GetYaxis()->SetTitleSize(0.06);
  h_Data->GetYaxis()->SetLabelSize(0.06);
  h_Data->SetMarkerStyle(20);
  h_Data->SetMarkerSize(1.2);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad * upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetLogy();

  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);


  h_Data->GetXaxis()->SetRangeUser(0,200);
  h_Data->SetMinimum(100);
  h_Data  -> Draw("e1");
  stack   -> Draw("hist same");
  h_Data  -> Draw("e1 same");
  bkgdErr -> Draw("e2 same");


  TLegend * leg = new TLegend(0.6,0.45,0.84,0.78);
  SetLegendStyle(leg);
  leg->SetTextSize(0.046);
  leg->SetHeader("track p_{T} = " + ptRange + " GeV");
  leg->AddEntry(h_Data,"Data","lp");
  leg->AddEntry(h_ZTT,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(h_QCD,"QCD (data driven)","f");
  leg->AddEntry(h_EWK,"electroweak","f");
  leg->AddEntry(h_Z,"Drell-Yan","f");
  leg->AddEntry(h_TT,"t#bar{t}+single top","f");
  leg->Draw();
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(upper,4,33); 

  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1F * ratioH = (TH1F*)h_Data->Clone("ratioH");
  ratioH->GetYaxis()->SetRangeUser(0.6,1.4);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetTitle("m_{inv}^{trk,#mu} [GeV]");
  ratioH->GetYaxis()->SetTitle("Obs./Exp.");
  ratioH->GetYaxis()->CenterTitle();
  ratioH->GetXaxis()->SetTitleOffset(3.5);
  ratioH->GetXaxis()->SetLabelSize(30);

  ratioH->Divide((TH1D*)stack->GetStack()->Last());
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  for(int i=1; i<=bkgdErr->GetNbinsX(); i++){
    ratioErrH->SetBinError(i,bkgdErr->GetBinError(i)/bkgdErr->GetBinContent(i));
    if(bkgdErr->GetBinContent(i)==0) ratioErrH->SetBinError(i,0);
    ratioErrH->SetBinContent(i,1);
  }

  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.5,1.5);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.5);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  // ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.31);
  lower->SetBottomMargin(0.32);
  lower->Draw();
  lower->cd();
  lower->SetGridy();

  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);


  ratioH->Draw("e1");
  ratioErrH->Draw("e2same");
    
  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Update();
  canv1->Print("figures/"+inputFileName+"_postfit.png");
  //canv1->Print("figures/"+inputFileName+"_postfit.pdf","Portrait pdf");
  
}
