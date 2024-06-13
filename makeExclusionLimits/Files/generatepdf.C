#include"TFile.h"
#include"TH1.h"
#include"TVirtualPad.h"
void generatepdf(){
  char* plot_name = new char[200];
  char* object_name = new char[200];
  char* path = new char[200];
  char* path1 = new char[200];
  char* path2 = new char[200];
  char* rootfile = new char[200];
  char* ext = new char[200];

  // std::string c[20]={"TTJetsHT","TTGJets","WJetsToLNu","ZJetsToNuNu","QCD","GJets","WGJetsToLNuG","T5bbbbZg_1800_150_FastSim","T5bbbbZg_1800_1750_FastSim","TChiWg_0_800_FastSim"};
  std::string c[20]={"QCD"};

  for(int j=0;j<1;j++) {

  // sprintf(rootfile,"rootoutput/%s.root",plot_name);
  // TFile *f1=new TFile("rootoutput/T5bbbbZg_1800_150_FastSim.root");
  sprintf(ext,"_v17");
  //  sprintf(rootfile,"rootoutput/%s%s.root",plot_name,ext);
    sprintf(rootfile,"rootoutput/%s%s.root",c[j].c_str(),ext);
   sprintf(plot_name,"%s",c[j].c_str());

  TFile *f1=new TFile(rootfile);
  
  // std::string a[20]={"nJets","MET","BestPhotonPt","nBTags","ST","mTPhoMET"};
  std::string a[20]={"nJets_nocut","MET_nocut","BestPhotonPt_nocut","nBTags_nocut","ST_nocut","mTPhoMET_nocut"};
  std::string b[20]={"BestPhotonPtvsMET","GenPhotonvsGenMET"};

  // TH1F *h1=(TH1F*)f1->FindObjectAny("selectBaselineYields_");
  //TH1F *h2=(TH1F*)f1->FindObjectAny("selectBaselineYields_");
  // cout<<"Size of a :"<< sizeof(a)<<endl;
  // cout<<"Size of a[0] :"<< sizeof(a[0])<<endl;
  // cout<<"Size of a/a[0] :"<< sizeof(a)/sizeof(a[0])<<endl;
  for(int i=0; i<6; i++){
   TCanvas *c1 = new TCanvas("c1", "c1",65,24,1680,1056);
   c1->Range(-1.981819,-20.70421,12.90736,118.2711);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetBottomMargin(0.1489776);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);

    sprintf(object_name,"%s",a[i].c_str());
  
  //   TH1F *h1=(TH1F*)f1->FindObjectAny("cutflows");
  //   TH1F *h2=(TH1F*)f1->FindObjectAny("nJets");
  // TH1F *h3=(TH1F*)f1->FindObjectAny("MET");
  // TH1F *h4=(TH1F*)f1->FindObjectAny("BestPhotonPt");
  // TH1F *h5=(TH1F*)f1->FindObjectAny("nBTags");
  // TH1F *h6=(TH1F*)f1->FindObjectAny("ST");

  TH1F *h1=(TH1F*)f1->FindObjectAny(object_name);

  if(strcmp(a[i].c_str(),"MET")==0) h1->RebinX(5);
  if(strcmp(a[i].c_str(),"ST")==0) h1->RebinX(5);
  if(strcmp(a[i].c_str(),"BestPhotonPt")==0) h1->RebinX(5);
  c1->cd();
  h1->Draw("hist");
  // sprintf(plot_name,"%s_%s",c[j].c_str(),a[i].c_str());

  TLegend* legends = new TLegend(0.5, 0.70, 0.7, 0.8,"","brNDC"); // the numbers determine the position of the box 
  legends->SetFillColor(0); 
  legends->AddEntry(h1,a[i].c_str(),"");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
  legends->SetTextSize(0.035);
   //   legends->SetMarkerStyle(1);
  legends->Draw();

  
  sprintf(path,"plots/pdf/%s_%s_plots.pdf",plot_name,object_name);
  sprintf(path1,"plots/png/%s_%s_plots.png",plot_name,object_name);
  sprintf(path2,"plots/gif/%s_%s_plots.gif",plot_name,object_name);
  c1->SaveAs(path);
  c1->SaveAs(path1);
  c1->SaveAs(path2);
  }  

  // for(int i=0; i<2; i++){
  //  TCanvas *c2 = new TCanvas("c2", "c2",65,24,1680,1056);
  //  c2->Range(-1.981819,-20.70421,12.90736,118.2711);
  //  c2->SetFillColor(0);
  //  c2->SetBorderMode(0);
  //  c2->SetBorderSize(2);
  //  c2->SetGridx();
  //  c2->SetBottomMargin(0.1489776);
  //  c2->SetFrameBorderMode(0);
  //  c2->SetFrameBorderMode(0);

  //  sprintf(object_name,"%s",b[i].c_str());
   
  //  TH2F *h2=(TH2F*)f1->FindObjectAny(object_name);
  //  h2->Rebin2D(5,5);
  //  c2->cd();
  //  h2->Draw("colz");
  //  // sprintf(plot_name,"%s_%s",c[j].c_str(),a[i].c_str());
  //  sprintf(path,"plots/pdf/%s_%s_plots.pdf",plot_name,object_name);
  //  sprintf(path1,"plots/png/%s_%s_plots.png",plot_name,object_name);
  //  sprintf(path2,"plots/gif/%s_%s_plots.gif",plot_name,object_name);
  //  c2->SaveAs(path);
  //  c2->SaveAs(path1);
  //  c2->SaveAs(path2);
  // }  
  
  }

}
