#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TColor.h"
#include "TStyle.h"

double xserror[]={26.3602,20.0201,15.4539,12.0956,9.61659,7.73547,6.2389,5.05005,4.16461,3.46763,2.88065,2.40183,2.04438,1.70195,1.43519,1.21937,1.04092,0.885243,0.769988,0.654544,0.564061,0.478362,0.412315,0.366257,0.314019,0.26801,0.242682,0.217618,0.175604,0.155683,0.133077,0.117638,0.107178,0.0917503,0.0828113,0.0644736,0.0579252,0.0561888,0.046491,0.0430197,0.0339561};

void plot_limit(int met1, int met2){//main 
 
   TCanvas *thecanvas = new TCanvas("thecanvas", "the canvas",4,69,600,550);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   thecanvas->Range(5.587935e-06,-3.763583,937.5,2.11013);
   thecanvas->SetFillColor(0);
   thecanvas->SetBorderMode(0);
   thecanvas->SetBorderSize(2);
   thecanvas->SetLogy();
   thecanvas->SetGridx();
   thecanvas->SetGridy();
   thecanvas->SetTickx(1);
   thecanvas->SetTicky(1);
   thecanvas->SetLeftMargin(0.16);
   thecanvas->SetRightMargin(0.04);
   thecanvas->SetTopMargin(0.07);
   thecanvas->SetBottomMargin(0.13);
   thecanvas->SetFrameFillStyle(0);
   thecanvas->SetFrameBorderMode(0);
   thecanvas->SetFrameFillStyle(0);
   thecanvas->SetFrameBorderMode(0);
   
   TH1F *frame__1 = new TH1F("frame__1","",1000,300,1300);
   frame__1->SetMinimum(0.1);
   frame__1->SetMaximum(1000);
   frame__1->SetDirectory(0);
   frame__1->SetLineStyle(0);
   frame__1->SetMarkerStyle(20);
   frame__1->GetXaxis()->SetTitle("Chargino mass [GeV]");
   frame__1->GetXaxis()->SetNdivisions(505);
   frame__1->GetXaxis()->SetLabelFont(42);
   frame__1->GetXaxis()->SetLabelOffset(0.007);
   frame__1->GetXaxis()->SetLabelSize(0.05);
   frame__1->GetXaxis()->SetTitleSize(0.06);
   frame__1->GetXaxis()->SetTitleOffset(0.9);
   frame__1->GetXaxis()->SetTitleFont(42);
   frame__1->GetYaxis()->SetTitle("95% CL upper limit on #sigma x BR [pb]");
   frame__1->GetYaxis()->SetLabelFont(42);
   frame__1->GetYaxis()->SetLabelOffset(0.007);
   frame__1->GetYaxis()->SetLabelSize(0.05);
   frame__1->GetYaxis()->SetTitleSize(0.06);
   frame__1->GetYaxis()->SetTitleOffset(1.25);
   frame__1->GetYaxis()->SetTitleFont(42);
   frame__1->GetZaxis()->SetLabelFont(42);
   frame__1->GetZaxis()->SetLabelOffset(0.007);
   frame__1->GetZaxis()->SetLabelSize(0.05);
   frame__1->GetZaxis()->SetTitleSize(0.06);
   frame__1->GetZaxis()->SetTitleFont(42);
   frame__1->Draw(" ");

   TGraph *p_theory = new TGraph(40);
   p_theory->SetName("Graph");
   p_theory->SetTitle("Graph");
   p_theory->SetFillColor(kBlue);
	 p_theory->SetLineWidth(3);

   TGraph *tchiwg_limit = new TGraph(40);
   tchiwg_limit->SetName("Graph");
   tchiwg_limit->SetTitle("Graph");
   tchiwg_limit->SetFillColor(3);

   TGraphAsymmErrors *gre = new TGraphAsymmErrors(40);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(kGreen);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(20);
   TGraphAsymmErrors *gre2 = new TGraphAsymmErrors(40);
   gre2->SetName("Graph2");
   gre2->SetTitle("Graph2");
   gre2->SetFillColor(5);
   gre2->SetLineWidth(3);
   gre2->SetMarkerStyle(20);

	 p_theory->SetLineColor(kBlue);
   p_theory->SetLineStyle(9);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000cc");
   tchiwg_limit->SetLineColor(ci);
   tchiwg_limit->SetLineStyle(9);
   tchiwg_limit->SetLineWidth(3);
   tchiwg_limit->SetMarkerStyle(20);

	 TFile xSecFile("susyxSec.root");
	 TH1F *p_crosssection = (TH1F*)xSecFile.Get("p_charginoSec");
   p_crosssection->GetXaxis()->SetTitle("M_chargino");

	 int charginomass(0);
	 float rvalue(0), sigma1up(0), sigma1lo(0), sigma2up(0), sigma2lo(0);
	 std::stringstream filename;
   filename.str("");
   filename<< "./tchiwg/limitvalue_" << met1 << "_" << met2 << ".txt";
	 std::ifstream Limit_file(filename.str().c_str());
	
	 if(Limit_file.is_open()){
		 for(unsigned i(0); i<40; i++){
			Limit_file >> charginomass >> rvalue >> sigma1up >> sigma1lo >> sigma2up >> sigma2lo;
      //std::cout << charginomass << " " << rvalue << " " << sigma1up << " " << sigma1lo << " " << sigma2up << " " << sigma2lo << std::endl;
			float theoryXS = p_crosssection->GetBinContent( p_crosssection->FindBin(charginomass));
			tchiwg_limit->SetPoint(i,charginomass,rvalue*theoryXS);
			gre->SetPoint(i, charginomass,rvalue*theoryXS);
			gre->SetPointError(i, 0,0,  (rvalue-sigma1lo)*theoryXS, (sigma1up-rvalue)*theoryXS);
			gre2->SetPoint(i, charginomass,rvalue*theoryXS);
			gre2->SetPointError(i, 0,0, (rvalue-sigma2lo)*theoryXS, (sigma2up-rvalue)*theoryXS); 
      p_theory->SetPoint(i, charginomass, theoryXS);
		}
	}
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,300,1300);
   Graph_Graph1->SetMinimum(0.1);
   Graph_Graph1->SetMaximum(1000);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineStyle(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   tchiwg_limit->SetHistogram(Graph_Graph1);
   
   tchiwg_limit->Draw("c");
  
   TH1F *Graph_GraphA = new TH1F("Graph_GraphA","Graph",100,300,1300);
   Graph_GraphA->SetMinimum(0.1);
   Graph_GraphA->SetMaximum(1000);
   Graph_GraphA->SetDirectory(0);
   Graph_GraphA->SetStats(0);
   Graph_GraphA->SetLineStyle(0);
   Graph_GraphA->SetMarkerStyle(20);
   Graph_GraphA->GetXaxis()->SetLabelFont(42);
   Graph_GraphA->GetXaxis()->SetLabelOffset(0.007);
   Graph_GraphA->GetXaxis()->SetLabelSize(0.05);
   Graph_GraphA->GetXaxis()->SetTitleSize(0.06);
   Graph_GraphA->GetXaxis()->SetTitleOffset(0.9);
   Graph_GraphA->GetXaxis()->SetTitleFont(42);
   Graph_GraphA->GetYaxis()->SetLabelFont(42);
   Graph_GraphA->GetYaxis()->SetLabelOffset(0.007);
   Graph_GraphA->GetYaxis()->SetLabelSize(0.05);
   Graph_GraphA->GetYaxis()->SetTitleSize(0.06);
   Graph_GraphA->GetYaxis()->SetTitleOffset(1.25);
   Graph_GraphA->GetYaxis()->SetTitleFont(42);
   Graph_GraphA->GetZaxis()->SetLabelFont(42);
   Graph_GraphA->GetZaxis()->SetLabelOffset(0.007);
   Graph_GraphA->GetZaxis()->SetLabelSize(0.05);
   Graph_GraphA->GetZaxis()->SetTitleSize(0.06);
   Graph_GraphA->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_GraphA);

   TH1F *Graph_GraphB = new TH1F("Graph_GraphB","Graph",100,300,1300);
   Graph_GraphB->SetMinimum(0.1);
   Graph_GraphB->SetMaximum(1000);
   Graph_GraphB->SetDirectory(0);
   Graph_GraphB->SetStats(0);
   Graph_GraphB->SetLineStyle(0);
   Graph_GraphB->SetMarkerStyle(20);
   Graph_GraphB->GetXaxis()->SetLabelFont(42);
   Graph_GraphB->GetXaxis()->SetLabelOffset(0.007);
   Graph_GraphB->GetXaxis()->SetLabelSize(0.05);
   Graph_GraphB->GetXaxis()->SetTitleSize(0.06);
   Graph_GraphB->GetXaxis()->SetTitleOffset(0.9);
   Graph_GraphB->GetXaxis()->SetTitleFont(42);
   Graph_GraphB->GetYaxis()->SetLabelFont(42);
   Graph_GraphB->GetYaxis()->SetLabelOffset(0.007);
   Graph_GraphB->GetYaxis()->SetLabelSize(0.05);
   Graph_GraphB->GetYaxis()->SetTitleSize(0.06);
   Graph_GraphB->GetYaxis()->SetTitleOffset(1.25);
   Graph_GraphB->GetYaxis()->SetTitleFont(42);
   Graph_GraphB->GetZaxis()->SetLabelFont(42);
   Graph_GraphB->GetZaxis()->SetLabelOffset(0.007);
   Graph_GraphB->GetZaxis()->SetLabelSize(0.05);
   Graph_GraphB->GetZaxis()->SetTitleSize(0.06);
   Graph_GraphB->GetZaxis()->SetTitleFont(42);
   gre2->SetHistogram(Graph_GraphB);

   gre2->Draw("3");
   gre->Draw("3");
   tchiwg_limit->Draw("c");
   p_theory->Draw("c");

	filename.str("");
	filename << "TChiWG_limit_" << met1 << "_" << met2 << ".pdf";
  thecanvas->SaveAs(filename.str().c_str());

}


