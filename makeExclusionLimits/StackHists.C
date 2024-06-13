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

const int nfiles=9,nBG=8;    //Specify no. of files
TFile *f[nfiles];
//int col[10]={kOrange,kBlue,kTeal+9,kGray+1,kCyan,kOrange-9,kYellow+2,kRed,kMagenta+2,kMagenta};  //Specify Colors
//int col[11]={kOrange,kBlue,kTeal+9,kOrange-9,kGray+1,kCyan,kMagenta+2,kPink+1,kRed,kMagenta,kBlue+2};  //Specify Colors
//int col[11]={kPink+1,kTeal+9,kOrange,kYellow,kBlue,kCyan,kGreen,kMagenta+2,kRed,kMagenta,kBlue+2};  //Specify Colors
int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kMagenta+2,kCyan,kBlue,kRed,kBlue+2,kMagenta};  //Specify Colors b's

char name[100],name2[100];
TObjArray Hlist(0);
TCanvas *c_cA=new TCanvas("basicPlots","Basic Plots",1500,850);
TCanvas *c_cB=new TCanvas("varBins","varBins",1500,850);
void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void setOverflow(TH1D*);
void printInt(TH1D*,int,const char*);
TLegend *legend=new TLegend(0.6, 0.90,  0.98, 0.45);
void StackHists(){
  
  f[0] = new TFile("WGJetsToLNuG.root");
  f[1] = new TFile("TTGJets.root");
  f[2] = new TFile("TTJetsHT.root");
  f[3] = new TFile("WJetsToLNu.root");
  f[4] = new TFile("ZJetsToNuNu.root");
  f[5] = new TFile("ZGJetsToNuNuG.root");
  f[6] = new TFile("GJets.root");
  f[7] = new TFile("QCD.root");
  //f[8] = new TFile("HG_NLSP150.root");
  //  f[9] = new TFile("HG_NLSP1000.root");
  //  f[8] = new TFile("HG_NLSP150.root");
  // f[8] = new TFile("T1bbbb_ZG_mGl1600_NLSP150.root");
  // f[9] = new TFile("T1bbbb_ZG_mGl1600_NLSP1000.root");
  // f[10] = new TFile("T1bbbb_ZG_mGl1600_NLSP1550.root");
  f[8] = new TFile("T5bbbbZg_1600_150_FastSim.root");
  //  f[10] = new TFile("T1bbbb_ZG_mGl1600_NLSP150.root");
  
  // gStyle->SetOptStat("nemri");
  gStyle->SetOptStat("ni");
  gStyle->SetTextSize(2);
  //gStyle->SetOptStat(0);
  c_cA->Divide(3,2); 
  c_cB->Divide(2,2);
  legend->SetLineWidth(0);
  //  c_cB->SetLogy();
  TFile *fout = new TFile("BGvsSig.root","RECREATE");

  THStack *hs_ST=new THStack("ST_Stack","ST Stacked");
  THStack *hs_MET=new THStack("MET_Stack","MET Stacked");
  THStack *hs_nJets=new THStack("nJets_Stack","nJets Stacked");
  THStack *hs_BTags=new THStack("BTags_Stack","BTags Stacked");
  THStack *hs_BestPhotonPt=new THStack("BestPhotonPt_Stack","Best Photon Pt Stacked");

  THStack *hs_STvBin=new THStack("ST_Stack_VarBins","ST Stacked with Variable Bin size");
  THStack *hs_METvBin=new THStack("MET_Stack_VarBins","MET Stacked with Variable Bin size");
  THStack *hs_BestPhotonPtvBin=new THStack("BestPhotonPt_Stack_VarBins","BestPhotonPt Stacked with Variable Bin size");
  
  cout<<"====================================================================================================================================================================="<<endl;
  cout.fill(' ');
  cout.width(20);
  cout<<left<<"Sample";cout.width(13);//for uncertainties replace all cout.width(13) with cout.width(25)
  cout<<left<<"[100,120)";cout.width(13);
  cout<<left<<"[120,160)";cout.width(13);
  cout<<left<<"[160,200)";cout.width(13);
  cout<<left<<"[200,270)";cout.width(13);
  cout<<left<<"[270,350)";cout.width(13);
  cout<<left<<"[350,450)";cout.width(13);
  cout<<left<<"[450,600)";cout.width(13);
  cout<<left<<"[600,750)";cout.width(13);
  cout<<left<<"[750,900)";cout.width(13);
  cout<<left<<"[900,Inf)";cout.width(13); cout<<endl;
  cout<<"====================================================================================================================================================================="<<endl;
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
    TH1D *h_ST=(TH1D*)f[i]->FindObjectAny("ST");
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny("MET");
    // TH1D *h_nJets=(TH1D*)f[i]->FindObjectAny("nJets");
    TH1D *h_nJets=(TH1D*)f[i]->FindObjectAny("nHadJets");
    TH1D *h_BTags=(TH1D*)f[i]->FindObjectAny("nBTags");
    TH1D *h_BestPhotonPt=(TH1D*)f[i]->FindObjectAny("BestPhotonPt");
    TH2D *h2_PhotonPtvsMET=(TH2D*)f[i]->FindObjectAny("BestPhotonPtvsMET");
   
    // TH1D *h_STvBin=(TH1D*)f[i]->FindObjectAny("MET_Rall");
    // TH1D *h_METvBin=(TH1D*)f[i]->FindObjectAny("MET_R2");
    // TH1D *h_BestPhotonPtvBin=(TH1D*)f[i]->FindObjectAny("MET_R3");
    TH1D *h_STvBin=(TH1D*)f[i]->FindObjectAny("STvarBin");
    TH1D *h_METvBin=(TH1D*)f[i]->FindObjectAny("METvarBin");
    TH1D *h_BestPhotonPtvBin=(TH1D*)f[i]->FindObjectAny("BestPhotonPtvarBin");

    h_ST->Rebin(5);
    h_MET->Rebin(5);
    h_BestPhotonPt->Rebin(5);
    printInt(h_MET,i,f[i]->GetName());

    decorate(h_ST,i,f[i]->GetName());
    decorate(h_MET,i,f[i]->GetName());
    decorate(h_nJets,i,f[i]->GetName());
    decorate(h_BTags,i,f[i]->GetName());
    decorate(h_BestPhotonPt,i,f[i]->GetName());

    decorate(h_STvBin,i,f[i]->GetName());
    decorate(h_METvBin,i,f[i]->GetName());
    decorate(h_BestPhotonPtvBin,i,f[i]->GetName());

    if(i<=(nBG-1)){
      hs_ST->          Add(h_ST);           
      hs_MET->         Add(h_MET);          
      hs_nJets->       Add(h_nJets);        
      hs_BTags->       Add(h_BTags);        
      hs_BestPhotonPt->Add(h_BestPhotonPt); 

      hs_STvBin ->         Add(h_STvBin);
      hs_METvBin->         Add(h_METvBin);
      hs_BestPhotonPtvBin->Add(h_BestPhotonPtvBin); 
    }
    if(i==nBG-1){
      c_cA->cd(1); hs_ST->          Draw("BAR HIST");hs_ST->          GetXaxis()->SetTitle("ST(GeV)"); decorate(hs_ST,          i,f[i]->GetName());   gPad->SetTickx();gPad->SetTicky();
      c_cA->cd(2); hs_MET->         Draw("BAR HIST");hs_MET->         GetXaxis()->SetTitle("MET(GeV)"); decorate(hs_MET,         i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      c_cA->cd(3); hs_nJets->       Draw("BAR HIST");hs_nJets->       GetXaxis()->SetTitle("no. of Jets"); decorate(hs_nJets,       i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      c_cA->cd(4); hs_BTags->       Draw("HIST BAR");hs_BTags->       GetXaxis()->SetTitle("no. of B-Tags"); decorate(hs_BTags,       i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      c_cA->cd(5); hs_BestPhotonPt->Draw("HIST BAR");hs_BestPhotonPt->GetXaxis()->SetTitle("Best Photon Pt(GeV)"); decorate(hs_BestPhotonPt,i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      //      c_cA->cd(5); hs_BestPhotonPt->Draw("HIST");hs_BestPhotonPt->GetXaxis()->SetTitle("Best Photon Pt(GeV)"); decorate(hs_BestPhotonPt,i,f[i]->GetName()); 

      c_cB->cd(1); hs_STvBin->      Draw("HIST BAR");hs_STvBin->      GetXaxis()->SetTitle("ST(GeV)");              decorate(hs_STvBin,     i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      c_cB->cd(2); hs_METvBin->     Draw("HIST BAR");hs_METvBin->     GetXaxis()->SetTitle("MET(GeV)");             decorate(hs_METvBin,    i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
      c_cB->cd(2); hs_METvBin->     Draw("HIST");//hs_METvBin->     GetXaxis()->SetTitle("MET(GeV)");             decorate(hs_METvBin,    i,f[i]->GetName());
      c_cB->cd(3); hs_BestPhotonPtvBin->Draw("HISTE BAR");hs_BestPhotonPtvBin->GetXaxis()->SetTitle("Best Photon Pt(GeV)");decorate(hs_BestPhotonPtvBin,    i,f[i]->GetName()); gPad->SetTickx();gPad->SetTicky();
    }
    if(i>=nBG){
      c_cA->cd(1); h_ST->          Draw("HISTE sames");
      c_cA->cd(2); h_MET->         Draw("HISTE sames");
      c_cA->cd(3); h_nJets->       Draw("HISTE sames");
      c_cA->cd(4); h_BTags->       Draw("HISTE sames");
      c_cA->cd(5); h_BestPhotonPt->Draw("HISTE sames");

      c_cB->cd(1); h_STvBin->      Draw("HISTE sames");
      c_cB->cd(2); h_METvBin->     Draw("HISTE sames");
      c_cB->cd(3); h_BestPhotonPtvBin->     Draw("HISTE sames");

    }

    if(i==nBG) cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout.fill(' ');
    cout.width(20);
    cout<<left<<f[i]->GetName();
    //cout.width(13);
   
    int k=0;
    for(int j=1;j<=h_METvBin->GetNbinsX();j++){
      //cout<<j<<" "<<h_METvBin->GetBinLowEdge(j)<<" "<<h_METvBin->GetBinContent(j)<<endl;
      if(h_METvBin->GetBinLowEdge(j)>=100){
	if(i<nBG){
	  if(i==0){Bcnt.push_back(h_METvBin->GetBinContent(j));}
	  else{Bcnt[k]+=h_METvBin->GetBinContent(j);k++;}
	}
	//	cout.width(8);//for uncertainties
	cout.width(13);
	//cout<<setprecision(3)<<left<<h_METvBin->GetBinContent(j);cout.width(0);cout<<" +- ";cout.width(13);cout<<left<<h_METvBin->GetBinError(j);//for uncertainties
	cout<<setprecision(3)<<left<<h_METvBin->GetBinContent(j);cout.width(13);
      }
    }//METvBinSize
    if(i==nBG-1){
      cout<<endl<<"---------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
      cout.width(20);
      cout<<left<<"Total BG";
      cout.width(13);
      for(unsigned int l=0;l<Bcnt.size();l++){cout<<setprecision(3)<<left<<Bcnt[l];cout.width(13);}
    }
    cout<<endl;
  }

  for(int i=0;i<nfiles;i++){
    TH1D *h_ST=(TH1D*)f[i]->FindObjectAny("ST");
    drawlegend(h_ST,i,f[i]->GetName());
  }
  //tl1.SetTextSize(25);
  c_cA->cd(1); gPad->SetLogy();legend->Draw(); 
  c_cA->cd(2); gPad->SetLogy();legend->Draw();
  c_cA->cd(3); gPad->SetLogy();legend->Draw();
  c_cA->cd(4); gPad->SetLogy();legend->Draw();
  c_cA->cd(5); gPad->SetLogy();legend->Draw();
  
  c_cB->cd(1); gPad->SetLogy(); //legend->Draw();
  c_cB->cd(2); gPad->SetLogy(); //legend->Draw();
  c_cB->cd(3); gPad->SetLogy(); //legend->Draw();
  c_cB->cd(4); legend->Draw();
    //}
  
  cout<<"====================================================================================================================================================================="<<endl;
  gPad->RedrawAxis();
  fout->cd();
  //Hlist.Write();
  c_cA->SaveAs("BasicPlots.png");
  //  c_cB->SaveAs("VarBinPlots.png");
  cout<<"*****************************************************"<<endl;
  cout<<"Int Lumi(inv.fb) for file1:"<<setprecision(4)<<intLumi<<endl;
}
void decorate(THStack *hs,int i,const char* fname){
  hs->SetMinimum(0.1);
  hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.04);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  gPad->Update();
  gStyle->SetOptStat(0);
}
void decorate(TH1D* hist,int i,const char* fname){
  hist->SetLineColor(col[i]);
  if(i<nBG) {
    hist->SetFillColor(col[i]);
    hist->SetLineColor(kBlack);
  }
  else hist->SetLineWidth(2);
  hist->SetTitle(0);
  hist->GetXaxis()->SetLabelSize(.04);
  hist->GetYaxis()->SetLabelSize(.05);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  // drawlegend(hist,i,fname);
  gPad->Update();
  setOverflow(hist);
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
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP150"){lName="4bZG_150";}
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP1000"){lName="4bZG_1000";}
  else if(lName=="T1bbbb_ZG_mGl1600_NLSP1550"){lName="4bZG_1550";}

  const char *l_name=lName.c_str();
  if(i<nBG)legend->AddEntry(hist,l_name,"f");
  else legend->AddEntry(hist,l_name,"l");
  legend->SetTextSize(0.05);
  // legend->SetTextColor(col[i]);
  //legend->SetLineColor(col[i]);
  //if(i==nfiles-1) legend->Draw();
}

void setOverflow(TH1D* hist){
  hist->SetBinContent(hist->GetNbinsX(), (hist->GetBinContent(hist->GetNbinsX()+1)) + (hist->GetBinContent(hist->GetNbinsX()) ) );
}

void printInt(TH1D* hist,int i,const char* fname){
  c_cA->cd(6);
  char name3[100];
  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(20);
  if(i==0)Tl.DrawText(.2, 0.9, "Sample                   Integral");
  Tl.SetTextColor(col[i]);
  float f=hist->Integral();

  string lName=fname;
  lName.pop_back();lName.pop_back();lName.pop_back();lName.pop_back();lName.pop_back();
  
  if(lName=="ZJetsToNuNu"){lName="Z( #rightarrow \nu #bar{\nu})+jets ";}
  //Tl.DrawText(.2, 0.8-0.05*i, fname  ); 
  sprintf(name, "%.2f",f);
  Tl.DrawText(.6, .8-0.05*i,   name); 
}

