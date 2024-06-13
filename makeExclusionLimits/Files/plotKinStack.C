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
const int nfiles=8,nBG=6;    //Specify no. of files
TFile *f[nfiles];
bool savePlots=1;
bool isPaper=0;
// bool v_17=false, v_12=true, v_12_vinay= false;
bool v_17=true, v_12=false, v_12_vinay= false;
// bool v_17=false, v_12=false, v_12_vinay=true;

//int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kMagenta+2,kBlue,kCyan,kRed,kBlue+2,kMagenta};  //Specify Colors b's
////int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kPink-2,kBlue,kCyan,kRed,kBlue+2,kMagenta};  //Specify Colors b's
//int col[11]={kPink-2,kTeal+9,kGreen,kYellow,kOrange,kBlue,kCyan,kRed,kBlue+2,kMagenta,kPink+1};  //Specify Colors b's
vector<int> col={kPink+1,kTeal+9,kYellow,kGray,kOrange,kCyan,kBlue,kRed,kBlue+2,kMagenta,kCyan};  //Specify Colors b's

TCanvas *c_cA=new TCanvas("kinVar","plot of a kin var",1500,900);

string getfname(const char *fname1){string fname=fname1;fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();return fname;}
void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
// TLegend *legend1=new TLegend(0.4, 0.75,  0.87, 0.88);
// TLegend *legend2=new TLegend(0.38, 0.6,  0.85, 0.72);

//TLegend *legend1=new TLegend(0.5660881,0.5702076,0.8491322,0.6996337);
TLegend *legend1=new TLegend(0.5253672,0.5299145,0.85247,0.6996337);
//TLegend *legend2=new TLegend(0.2803738,0.7350427,0.8190921,0.8669109);
TLegend *legend2=new TLegend(0.2,0.7350427,0.8190921,0.8669109);

void setLastBinAsOverFlow(TH1D*);
TH1D* setMyRange(TH1D*,double,double);

void plotKinStack(){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  double yMin=0.1,yMax = 1000;
  double xMin=0.0,xMax = 2000;
  // TString varName = "ST"; TString xLabel = "H_{T}^{#gamma} (GeV)";   int rebin=25; yMin=0.5,yMax = 1100; xMin=300.0,xMax = 3000;
  // TString varName = "METvBin2_nj2to4_nbjetnot0"; TString xLabel = "p_{T}^{miss} (GeV)";   int rebin=1; yMin=0.005,yMax = 15000; xMin=200.0,xMax = 1000;
  // TString varName = "MET"; TString xLabel = "p_{T}^{miss} (GeV)";   int rebin=10; yMin=0.1,yMax = 1000; xMin=100.0,xMax = 1400;
  // TString varName = "BestPhotonPt"; TString xLabel = "p_{T}^{#gamma} (GeV)";   int rebin=10; yMin=0.5,yMax = 100000; xMin=0.0,xMax = 1200;
  // TString varName = "PhovBin2"; TString xLabel = "p_{T}^{#gamma} (GeV)";   int rebin=1; yMin=0.4,yMax = 1000000; xMin=0.0,xMax = 1200;
   // TString varName = "nHadJets"; TString xLabel = "N_{ jets}";   int rebin=1; yMin=0.1,yMax = 1000; xMin=0.0,xMax = 14;
   // TString varName = "nJets"; TString xLabel = "N_{ jets}";   int rebin=1; yMin=0.1,yMax = 100000; xMin=1,xMax = 14;
   // TString varName = "nBTags"; TString xLabel = "N_{ b-jets}";   int rebin=1; yMin=0.1,yMax = 100000; xMin=1.0,xMax = 5;
   // TString varName = "mindPhi1dPhi2"; TString xLabel = "min(#Delta#phi_{1}, #Delta#phi_{2})";   int rebin=10; yMin=0.3,yMax = 1500; xMin=-100000.0,xMax = 100000;
  //"mindPhi1dPhi2";//"AllSBins_v7_CD";
  TString varName = "AllSBins_v7_CD"; TString xLabel = "Bin no.";   int rebin=1; yMin=0.5,yMax = 1000000; xMin=-100000,xMax = 100000;
  //TString varName = "METvarBin";
  //  TString xLabel = "p_{T}^{miss} (GeV)";//min(#Delta#Phi_{1},#Delta#Phi_{2})
  //TString varName = "dPhi_METjet1"; TString xLabel = "dPhi b/w Met & Jet1";   int rebin=1; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
  //  // "mindPhi1dPhi2";//"AllSBins_v7_CD";
  //   TString varName = "dPhi_METjet2"; TString xLabel = "dPhi b/w Met & Jet2";   int rebin=1; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
  //   TString varName = "dPhi_phojet1"; TString xLabel = "dPhi b/w pho & Jet1";   int rebin=1; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
   // TString varName = "dPhi_phojet2"; TString xLabel = "dPhi b/w pho & Jet2";   int rebin=1; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
   // TString varName = "dPhi_phoMET"; TString xLabel = "dPhi b/w pho & MET";   int rebin=1; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
   // TString varName = "mTPhoMET"; TString xLabel = "mT b/w pho & MET";   int rebin=10; yMin=0.3,yMax = 15000; xMin=-100000,xMax = 100000;
  if(v_17){
  f[0] = new TFile("TTGJets_v17.root");
  f[1] = new TFile("TTJetsHT_v17.root");//ZGJetsToNuNuG_v17.root
  f[2] = new TFile("ZGZJ_NuNuG.root");
  // f[2] = new TFile("ZJetsToNuNu_v17.root");
  // f[2] = new TFile("ZGJetsToNuNuG_v12.root");
  f[3] = new TFile("WGJetsToLNuG_v17.root");
  f[4] = new TFile("WJetsToLNu_v17.root");
  f[5] = new TFile("GJetsQCD_new_v17.root");
  col.resize(0);
  col={kGray,kTeal+9,kOrange,kRed,kCyan-1,kCyan,kBlue,kMagenta+2,kPink+1,kMagenta,kBlack};
  // f[6] = new TFile("TChiWg_0_400_FastSim_v17.root");
  // f[7] = new TFile("TChiWg_0_800_FastSim_v17.root");
  // f[8] = new TFile("TChiWg_0_1200_FastSim_v17.root");
  f[6] = new TFile("T5bbbbZg_1800_150_FastSim_v17.root");
  // f[7] = new TFile("T5bbbbZg_1800_1000_FastSim_v17.root");
  f[7] = new TFile("T5bbbbZg_1800_1750_FastSim_v17.root");
  // f[6] = new TFile("GGM_M1M3_1100_1000_FastSim.root");
  // f[7] = new TFile("GGM_M1M3_1100_200_FastSim.root");
  }
  if(v_12){
  f[0] = new TFile("TTGJets_v12.root");
  f[1] = new TFile("TTJetsHT_v12.root");//ZGJetsToNuNuG_v17.root
  // f[2] = new TFile("ZGZJ_NuNuG_v12.root");
  // f[2] = new TFile("ZJetsToNuNu_v17.root");
   f[2] = new TFile("ZGJetsToNuNuG_v12.root");
  f[3] = new TFile("WGJetsToLNuG_v12.root");
  f[4] = new TFile("WJetsToLNu_v12.root");
  f[5] = new TFile("GJetsQCD_new_v12.root");
  col.resize(0);
  col={kGray,kTeal+9,kOrange,kRed,kCyan-1,kCyan,kBlue,kMagenta+2,kPink+1,kMagenta,kBlack};
  // f[6] = new TFile("TChiWg_0_400_FastSim_v17.root");
  // f[7] = new TFile("TChiWg_0_800_FastSim_v17.root");
  // f[8] = new TFile("TChiWg_0_1200_FastSim_v17.root");
  f[6] = new TFile("T5bbbbZg_1800_150_FastSim_v12.root");
  // f[7] = new TFile("T5bbbbZg_1800_1000_FastSim_v17.root");
  f[7] = new TFile("T5bbbbZg_1800_1750_FastSim_v12.root");
  // f[6] = new TFile("GGM_M1M3_1100_1000_FastSim.root");
  // f[7] = new TFile("GGM_M1M3_1100_200_FastSim.root");
  }
  if(v_12_vinay){
  f[0] = new TFile("TTGJets_v12_vinay.root");
  f[1] = new TFile("TTJetsHT_v12_vinay.root");//ZGJetsToNuNuG_v17.root
  f[2] = new TFile("ZGZJ_NuNuG_v12_vinay.root");
  // f[2] = new TFile("ZJetsToNuNu_v17.root");
  //  f[2] = new TFile("ZGJetsToNuNuG_v12_vinay.root");
  f[3] = new TFile("WGJetsToLNuG_v12_vinay.root");
  f[4] = new TFile("WJetsToLNu_v12_vinay.root");
  f[5] = new TFile("GJetsQCD_new_v12_vinay.root");
  col.resize(0);
  col={kGray,kTeal+9,kOrange,kRed,kCyan-1,kCyan,kBlue,kMagenta+2,kPink+1,kMagenta,kBlack};
  // f[6] = new TFile("TChiWg_0_400_FastSim_v17.root");
  // f[7] = new TFile("TChiWg_0_800_FastSim_v17.root");
  // f[8] = new TFile("TChiWg_0_1200_FastSim_v17.root");
  f[6] = new TFile("T5bbbbZg_1800_150_FastSim_v12_vinay.root");
  // f[7] = new TFile("T5bbbbZg_1800_1000_FastSim_v17.root");
  f[7] = new TFile("T5bbbbZg_1800_1750_FastSim_v12_vinay.root");
  // f[6] = new TFile("GGM_M1M3_1100_1000_FastSim.root");
  // f[7] = new TFile("GGM_M1M3_1100_200_FastSim.root");
  }
  
  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","MET Stacked");
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
    
    TH1D *h_MET;
    if(i<=nBG-1) h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    if(i>=nBG) h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    
    if (varName=="AllSBins_v7_CD"){
      cout<<"filename ------->  "<<getfname(f[i]->GetName())<<endl; 
      for(int k=1; k<h_MET->GetSize()-1; k++)
	{
	  cout<<"Bin content for bin "<<k<<" ---->  "<<h_MET->GetBinContent(k)<<endl;
	}
      cout<<"==========================================="<<endl;
    }
    //    couh_MET->GetNbinsX()<<" ";
    h_MET->Rebin(rebin);
    //    h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    //    h_MET->SetMinimum(100);
    decorate(h_MET,i,f[i]->GetName());
    h_MET = setMyRange(h_MET,xMin,xMax);
    //    h_MET->GetXaxis()->SetRangeUser(xMin,xMax);
    
    if(i<=(nBG-1))  hs_var->Add(h_MET);

    if(i==nBG-1) {
      c_cA->cd();
      hs_var->Draw("BAR HIST");
      hs_var->Draw("HIST");
      hs_var->SetMinimum(yMin);
      hs_var->SetMaximum(yMax);
      decorate(hs_var,i,f[i]->GetName()); 
      if(xMin > -10000 && xMax < 10000) hs_var->GetXaxis()->SetRangeUser(xMin-0.1,0.1+xMax);
    }
    if(i>=nBG){ 
      c_cA->cd(); 
      h_MET->SetMarkerStyle(20);
      h_MET->SetMarkerColor(col[i]);
      h_MET->SetLineColor(col[i]);
      h_MET->SetLineWidth(3);
      if(i>=9)  h_MET->SetLineStyle(2);
      h_MET->Draw("hist same");
      //      h_MET->GetYaxis()->SetRangeUser(0.5,20000);
      //      h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    }
    drawlegend(h_MET,i,f[i]->GetName());
    if(i==nfiles-1){ 
      hs_var->GetXaxis()->SetTitleOffset(1.0);
      hs_var->GetXaxis()->SetTitle(xLabel); hs_var->GetYaxis()->SetTitle("Events");hs_var->SetTitle(0);
      hs_var->GetYaxis()->SetTitleOffset(.90);
      TString temp=h_MET->GetName(),temp2;
      if(temp.Contains("nHadJets") || temp.Contains("nBTags")){
	gPad->SetTickx(0);
	hs_var->GetXaxis()->SetLabelSize(0.08);
	for(int i=1;i<=h_MET->GetNbinsX();i++){
	  temp2 = to_string(i-1);
	  if(i%2==0 && temp.Contains("nHadJets")) continue;
	    hs_var->GetXaxis()->SetBinLabel(i,temp2);
	}
	//      cout<<hist->GetName()<<endl;
      }
    }
    
  }

  legend1->SetFillStyle(0); legend2->SetFillStyle(0);  
  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  legend2->SetBorderSize(0);
  legend2->SetMargin(0.12);
  c_cA->cd(); gPad->SetLogy();legend1->Draw();
  c_cA->cd(); gPad->SetLogy();legend2->Draw();
  //  gPad->RedrawAxis();
  //  hs_var->GetXaxis()->SetTitle(xLabel);
 
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  if(isPaper) textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Supplementary}}");
  else textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.7,0.91,name2);
  TLatex Tl;
  Tl.SetTextSize(0.04);
  // if(isPaper) Tl.DrawLatexNDC(0.48,0.91,"#bf{arXiv:xxxx.xxxxx}");

  if(varName == "mindPhi1dPhi2"){
    TLine *line1=new TLine( 0.3,0.11,  0.3,yMax);
    line1->Draw();
    line1->SetLineStyle(2);
    TArrow *arrow1 = new TArrow(0.3,100,1.2,100,0.01,"|>");
    arrow1->Draw();
    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatex(0.33,140,"#bf{Signal Region}");
    //  Tl.SetTextSize(0.04);
    //  Tl.DrawLatexNDC(0.48,0.91,"#bf{arXiv:xxxx.xxxxx}");
  }
  // //  c_cB->SaveAs("searchBins.png");
  //------------------------------------
  c_cA->cd(); c_cA->SetGridx(0); c_cA->SetGridy(0);
  if(varName=="AllSBins_v7_CD"){
    TLine *line1V7=new TLine( 6.5,0.5,  6.5,3000);
    TLine *line2V7=new TLine(11.5,0.5, 11.5,3000);
    TLine *line3V7=new TLine(16.5,0.5, 16.5,1000);
    TLine *line4V7=new TLine(21.5,0.5, 21.5,300);
    TLine *line5V7=new TLine(26.5,0.5, 26.5,300);
    

    c_cA->cd(); c_cA->SetGridx(0); c_cA->SetGridy(0);
    line1V7->Draw();      line2V7->Draw();  line3V7->Draw();
    line4V7->Draw();      line5V7->Draw();

    TArrow *arrow1 = new TArrow( 0.5,1000, 6.5,1000,0.01,"<|>");
    TArrow *arrow2 = new TArrow( 6.5,1000,11.5,1000,0.01,"<|>");
    TArrow *arrow3 = new TArrow(11.5,650,16.5,650,0.01,"<|>");
    TArrow *arrow4 = new TArrow(16.5,300, 21.5,300,0.01,"<|>");
    TArrow *arrow5 = new TArrow(21.5,300, 26.5,300,0.01,"<|>");
    TArrow *arrow6 = new TArrow(26.5,300, 31.5,300,0.01,"<|>");

    arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
    arrow4->Draw(); arrow5->Draw(); arrow6->Draw();

    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatex(3.,1500,"N^{ 0}_{ 2-4}");
    Tl.DrawLatex(8.5,1500,"N^{ 0}_{ 5-6}");
    Tl.DrawLatex(13.5,850,"N^{ 0}_{ #geq7}");
    Tl.DrawLatex(18.5,400,"N^{ #geq1}_{ 2-4}");
    Tl.DrawLatex(23.5,400,"N^{ #geq1}_{ 5-6}");
    Tl.DrawLatex(28.5,400,"N^{ #geq1}_{ #geq7}");
  }
  //------------------------------------
  if(savePlots){
    TString saveName = "supp_Sim_"+varName;
    TString modelName = f[6]->GetName();
    if(modelName.Contains("T5bbbb")) modelName = "T5bbbbZG";
    else if(modelName.Contains("T5qqqq")) modelName = "T5qqqqHG";
    else if(modelName.Contains("T5tttt")) modelName = "T5ttttZG";
    else if(modelName.Contains("T6tt")) modelName = "T6ttZG";
    saveName = saveName+"_"+modelName+".pdf";
    c_cA->SaveAs(saveName);
  }
}

void decorate(THStack *hs,int i,const char* fname){
  //  hs->SetMinimum(0.5);
  //hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.05);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  //  gPad->Update();
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
  hist->GetXaxis()->SetLabelSize(.06);
  hist->GetYaxis()->SetLabelSize(.06);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.06);
  // drawlegend(hist,i,fname);
  //  gPad->Update();
  setLastBinAsOverFlow(hist);
  gStyle->SetOptStat(0);
  
  //Hlist.Add(hist);
}

void drawlegend(TH1D *hist,int i,const char* fname){
  gStyle->SetLegendBorderSize(0);
 
  TString lName=fname;
  
  if(lName.Contains("ZGZJ")){lName="Z(#nu#bar{#nu}) + #gamma";}
  else if(lName.Contains("ZJets")){lName="Z(#nu#bar{#nu}) + jets";}
  else if(lName.Contains("DYJetsToLL")){lName="DY(l^{+}l^{-})";}
  else if(lName.Contains("WJetsToLNu")){lName="W(l#nu) + jets";}
  else if(lName.Contains("RareProcess")){}
  else if(lName.Contains("TTJetsHT")){lName="t #bar{t}";}
  else if(lName.Contains("WGJetsToLNuG")){lName="W(l#nu) + #gamma";}
  else if(lName.Contains("ZGJetsToNuNuG")){lName="Z(#nu#bar{#nu}) + #gamma";}
  else if(lName.Contains("TTGJets")){lName="t #bar{t} + #gamma";}
  //  else if(lName.Contains("QCD")){lName="QCD";}
  else if(lName.Contains("GJets")){lName="#gamma + jets";}
  else if(lName.Contains("Run2016")){lName="Data";}
  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="T5bbbbZg 1600, 150";}
  else if(lName.Contains("T5bbbbZg_1600_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5bbbbZG_1800_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZG_1800_1750")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}

  else if(lName.Contains("T5bbbbZg_1800_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZg_1800_1750")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}
  else if(lName.Contains("T5bbbbZg_1800_1000")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)";}
  else if(lName.Contains("T5qqqqHg_1800_150")){lName = "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5qqqqHg_1800_1750")){lName = "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_150")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_1550")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_1000")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)";}

  else if(lName.Contains("T6ttZg_1000_100")){lName = "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)";}
  else if(lName.Contains("T6ttZg_1000_900")){lName = "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 900 GeV)";}

  // else if(lName.Contains("T5bbbbZg_1600_150")){lName = "T5bbbbZG (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  // else if(lName.Contains("T5bbbbZg_1600_1550")){lName = "T5bbbbZG (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName="T5bbbbZg 1600, 1550";}
  else if(lName.Contains("T5qqqqHg_1600_1550")){lName="T5qqqqHg_1550";}
  else if(lName.Contains("T5qqqqHg_1600_150")){lName="T5qqqqHg_150";}
  else if(lName.Contains("TChiWg_0_400")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 400 GeV";}
  else if(lName.Contains("TChiWg_0_800")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 800 GeV";}
  else if(lName.Contains("TChiWg_0_1200")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 1200 GeV";}

  // const char *l_name=lName.c_str();
  if(i<nBG)legend1->AddEntry(hist,lName,"f");
  else legend2->AddEntry(hist,lName,"l");
  // legend1->SetTextSize(0.04);
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

TH1D* setMyRange(TH1D *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;
  // h1->Print("all");
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);
  //  cout<<nMax<<endl;
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
    //    cout<<":";
  }
  //  h1->Print("all");
  //  cout<<endl;
  return h1;
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);                                                                                                                      
}
