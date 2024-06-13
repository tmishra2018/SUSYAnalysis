#include"TFile.h"
#include"TAxis.h"
#include"TH1.h"
#include"TVirtualPad.h"
void generate_overlaypdf(){
  char* plot_name = new char[200];
  char* object_name = new char[200];
  char* path = new char[200];
  char* path1 = new char[200];
  char* path2 = new char[200];
  char* rootfile = new char[200];
  char* ext = new char[200];
  // vector<int> col={kGray,kTeal+9,kOrange,kRed,kCyan-1,kCyan,kBlue,kMagenta+2,kPink+1,kMagenta,kBlack};


  // std::string c[20]={"TTJetsHT","TTGJets","WJetsToLNu","ZJetsToNuNu","QCD","GJets","WGJetsToLNuG","T5bbbbZg_1800_150_FastSim","T5bbbbZg_1800_1750_FastSim","TChiWg_0_800_FastSim"};
  // std::string c[20]={"TTGJets","TTJetsHT","ZGJetsToNuNuG","WGJetsToLNuG","WJetsToLNu","GJetsQCD_new"};
  std::string c[20]={"T5bbbbZg_1800_1750_FastSim"};//,"T5bbbbZg_1800_1750_FastSim","TChiWg_0_800_FastSim"};
  // std::string a[20]={"nJets_num1","MET_num1","BestPhotonPt_num1","nBTags_num1","ST_num1","mTPhoMET_num1"};
  // std::string c[20]={"TTGJets","TTJetsHT","ZGJetsToNuNuG","WGJetsToLNuG","WJetsToLNu","GJetsQCD_new","TChiWg_0_800_FastSim","T5bbbbZg_1800_150_FastSim","T5bbbbZg_1800_1750_FastSim"};
  // std::string a[20]={"nJets","MET","BestPhotonPt","nBTags","ST","mTPhoMET"};
  //  std::string a[20]={"nJets"};
  std::string a[20]={"MET_nj1","MET_nj2to4","MET_nj5to7","MET_n8toinf"};
  std::string b[20]={"BestPhotonPtvsMET","GenPhotonvsGenMET"};
  for(int j=0;j<4;j++) {
   TCanvas *c1 = new TCanvas("c1", "c1",65,24,1680,1056);
   c1->Range(-1.981819,-20.70421,12.90736,118.2711);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetBottomMargin(0.1489776);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
  gStyle->SetOptStat(00);  
  TLegend* legends = new TLegend(0.5, 0.5, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
  legends->SetFillColor(0); 

   for(int i=0; i<1; i++){

  // sprintf(rootfile,"rootoutput/%s.root",plot_name);
  // TFile *f1=new TFile("rootoutput/T5bbbbZg_1800_150_FastSim.root");
      sprintf(ext,"_v17");
  //  sprintf(rootfile,"rootoutput/%s%s.root",plot_name,ext);
      sprintf(rootfile,"rootoutput/%s%s.root",c[i].c_str(),ext);
      sprintf(plot_name,"%s",c[i].c_str());

  TFile *f1=new TFile(rootfile);
  

  // TH1F *h1=(TH1F*)f1->FindObjectAny("selectBaselineYields_");
  //TH1F *h2=(TH1F*)f1->FindObjectAny("selectBaselineYields_");
  // cout<<"Size of a :"<< sizeof(a)<<endl;
  // cout<<"Size of a[0] :"<< sizeof(a[0])<<endl;
  // cout<<"Size of a/a[0] :"<< sizeof(a)/sizeof(a[0])<<endl;

    sprintf(object_name,"%s",a[j].c_str());
  
  //   TH1F *h1=(TH1F*)f1->FindObjectAny("cutflows");
  //   TH1F *h2=(TH1F*)f1->FindObjectAny("nJets");
  // TH1F *h3=(TH1F*)f1->FindObjectAny("MET");
  // TH1F *h4=(TH1F*)f1->FindObjectAny("BestPhotonPt");
  // TH1F *h5=(TH1F*)f1->FindObjectAny("nBTags");
  // TH1F *h6=(TH1F*)f1->FindObjectAny("ST");

  TH1F *h1=(TH1F*)f1->FindObjectAny(object_name);

  // if(strcmp(a[j].c_str(),"MET")==0)
  //   {h1->Rebin(5);
  //    h1->GetXaxis()->SetRangeUser(100,1800);
  //    h1->GetXaxis()->SetTitle("MET (GeV)");
  //   }
  if(strcmp(a[j].c_str(),"MET")==0 ||strcmp(a[j].c_str(),"MET_nocut")==0)
    {h1->Rebin(5);
     h1->GetXaxis()->SetRangeUser(100,1800);
     //     h1->GetYaxis()->SetRangeUser(0,10000);
     h1->SetMaximum(10000);
     h1->SetMinimum(0.1);
     h1->GetXaxis()->SetTitle("MET (GeV)");
    } 
  if(strcmp(a[j].c_str(),"nJets")==0)
    {h1->Rebin(1);
     h1->GetXaxis()->SetRangeUser(0,20);
     h1->GetXaxis()->SetTitle(a[j].c_str());
    }
  if(strcmp(a[j].c_str(),"ST")==0) h1->RebinX(10);
  if(strcmp(a[j].c_str(),"BestPhotonPt")==0) h1->RebinX(10);
  if(strcmp(a[j].c_str(),"mTPhoMET")==0) h1->Rebin(10);
   // h1->Rebin(1);
   c1->cd();
   // h1->SetLineColor(col[i]);
   h1->SetLineColor(i);
   h1->SetLineWidth(2);
   // h1->GetXaxis()->SetNDivisions(2);
   if(strcmp(c[i].c_str(),"T5bbbbZg_1800_150_FastSim")==0||strcmp(c[i].c_str(),"T5bbbbZg_1800_1750_FastSim")==0||strcmp(c[i].c_str(),"TChiWg_0_800_FastSim")==0) 
     {     h1->SetLineStyle(2);   h1->SetLineWidth(3);
     }
   // else if(strcmp(c[i].c_str(),"FastSim")!=0) 
   else
     h1->SetLineStyle(1);

  c1->SetLogy();

  if(i==0)
      h1->Draw("hist");
  if(i!=0)
      h1->Draw("hist same");
// sprintf(plot_name,"%s_%s",c[j].c_str(),a[i].c_str());
  legends->AddEntry(h1,c[i].c_str(),"l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   //   legends->SetMarkerStyle(1);
  legends->Draw();

   }   


   sprintf(path,"plots/pdf/%s_overlayplots.pdf",object_name);
   sprintf(path1,"plots/png/%s_overlayplots.png",object_name);
   sprintf(path2,"plots/gif/%s_overlayplots.gif",object_name);
   c1->SaveAs(path);
   c1->SaveAs(path1);
   c1->SaveAs(path2);
  }  

}



