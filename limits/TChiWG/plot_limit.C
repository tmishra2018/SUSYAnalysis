// source : /uscms_data/d3/mengleis/CMSSW_Work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/DataCard/TChiWG_Unblind_Data/test_limit.C
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
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
#include "TLegendEntry.h"
#include <vector>
#include <cstdlib>
#include "../../include/tdrstyle.C"

// source: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVn2x1wino#Envelope_of_CTEQ6_6_and_MSTW2008

double xserror[]={26.3602,20.0201,15.4539,12.0956,9.61659,7.73547,6.2389,5.05005,4.16461,3.46763,2.88065,2.40183,2.04438,1.70195,1.43519,1.21937,1.04092,0.885243,0.769988,0.654544,0.564061,0.478362,0.412315,0.366257,0.314019,0.26801,0.242682,0.217618,0.175604,0.155683,0.133077,0.117638,0.107178,0.0917503,0.0828113,0.0644736,0.0579252,0.0561888,0.046491,0.0430197,0.0339561};

// ============================================================================
// Limit-file reading (merged in from makelimit.sh)
//
// EXPECTED and OBSERVED limits are both taken from the SAME combine
// AsymptoticLimits output of the UNBLINDED cards (one run per mass point),
// exactly as makeT5WGLimit.py does:
//     Expected  2.5% / 16.0% / 50.0% / 84.0% / 97.5%  -> median + +-1/2 sigma bands
//     Observed Limit                                   -> observed
// Set kExpTag = "expected" to reproduce the old behaviour (expected from a
// separate cards_Run2_splusb run) for a side-by-side comparison.
// ============================================================================
static const std::string kLimitDir = "/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2/";
static const std::string kObsTag   = "observed";   // limit_<tag>_<mass>.txt from the unblind cards
static const std::string kExpTag   = "observed";   // same file -> same run as the observed limit
static const int  kMassStart = 300;
static const int  kMassStep  = 25;
static const int  kNMass     = 40;                 // 300, 325, ..., 1275
static const bool kWriteTextFiles = true;          // also write limitvalue_Run2_{expected,observed}.txt like makelimit.sh did

struct LimitPoint { int mass; double r, s1up, s1lo, s2up, s2lo; };

static std::vector<std::string> readLines(const std::string& path, bool& ok){
   std::vector<std::string> lines;
   std::ifstream fin(path.c_str());
   ok = fin.is_open();
   std::string l;
   while(ok && std::getline(fin, l)) lines.push_back(l);
   return lines;
}

// Finds the combine line "<key>: r < X" (e.g. key = "50.0%" or "Observed Limit").
// Requires X > 0 (a failed fit gives 0/negative).
static bool findLimit(const std::vector<std::string>& lines, const std::string& key, double& val){
   const std::string pat = key + ": r <";
   for(size_t k = 0; k < lines.size(); k++){
      if(lines[k].find(pat) == std::string::npos) continue;
      val = std::atof(lines[k].c_str() + lines[k].find('<') + 1);
      return val > 0;
   }
   return false;
}

// First crossing of (r * theoryXS) through theoryXS (i.e. r = 1), linearly
// interpolated in (limit - theory) vs mass. Returns -1 if there is no crossing.
static double interpolateExclusion(const std::vector<LimitPoint>& pts, TH1F* hxs){
   bool havePrev = false;
   double m1 = 0, y1 = 0;
   for(size_t k = 0; k < pts.size(); k++){
      double xs = hxs->GetBinContent(hxs->FindBin(pts[k].mass));
      double y  = pts[k].r * xs - xs;
      if(y < 0){ m1 = pts[k].mass; y1 = y; havePrev = true; }          // excluded
      else if(havePrev) return m1 + (pts[k].mass - m1) * (-y1) / (y - y1);
   }
   return -1;
}

void plot_limit(){//main 
   gROOT->SetBatch(kTRUE);
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
   
   thecanvas->cd();
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.7);
   pad2->SetTopMargin(0);
   pad2->SetLeftMargin(0.12); 
   pad2->SetBottomMargin(0.15);
   pad2->Draw();  
   pad2->cd(); 
   gPad->SetLogy(); 
   
   TH1F *frame__1 = new TH1F("frame__1","",1000,300,1300);
   frame__1->SetMinimum(0.1);
   frame__1->SetMaximum(2000);
   frame__1->SetDirectory(0);
   frame__1->SetLineStyle(0);
   frame__1->SetMarkerStyle(20);
   frame__1->GetXaxis()->SetTitle("M_{NLSP} (GeV)");
   frame__1->GetXaxis()->SetNdivisions(505);
   frame__1->GetXaxis()->SetLabelFont(42);
   frame__1->GetXaxis()->SetLabelOffset(0.007);
   frame__1->GetXaxis()->SetLabelSize(0.05);
   frame__1->GetXaxis()->SetTitleSize(0.06);
   frame__1->GetXaxis()->SetTitleOffset(0.9);
   frame__1->GetXaxis()->SetTitleFont(42);
   frame__1->GetYaxis()->SetTitle("cross section (fb) ");
   frame__1->GetYaxis()->SetLabelFont(42);
   frame__1->GetYaxis()->SetLabelOffset(0.007);
   frame__1->GetYaxis()->SetLabelSize(0.05);
   frame__1->GetYaxis()->SetTitleSize(0.06);
   frame__1->GetYaxis()->SetTitleOffset(1.0);
   frame__1->GetYaxis()->SetTitleFont(42);
   frame__1->GetZaxis()->SetLabelFont(42);
   frame__1->GetZaxis()->SetLabelOffset(0.007);
   frame__1->GetZaxis()->SetLabelSize(0.05);
   frame__1->GetZaxis()->SetTitleSize(0.06);
   frame__1->GetZaxis()->SetTitleFont(42);
   frame__1->Draw(" ");

   TGraph *p_theory = new TGraph(40);
   p_theory->SetName("theory");
   p_theory->SetTitle("theory");
   p_theory->SetLineColor(kBlue);
   p_theory->SetFillColor(0);
	 p_theory->SetLineWidth(4);
	 p_theory->SetLineStyle(1);

   TGraph *p_theoryup = new TGraph(40);
   p_theoryup->SetName("theoryup");
   p_theoryup->SetTitle("theoryup");
   p_theoryup->SetLineColor(kBlue);
   p_theoryup->SetFillColor(0);
	 p_theoryup->SetLineWidth(2);
	 p_theoryup->SetLineStyle(2);

   TGraph *p_theorydo = new TGraph(40);
   p_theorydo->SetName("theorydo");
   p_theorydo->SetTitle("theorydo");
   p_theorydo->SetLineColor(kBlue);
   p_theorydo->SetFillColor(0);
	 p_theorydo->SetLineWidth(2);
	 p_theorydo->SetLineStyle(2);

   TGraph *tchiwg_expect = new TGraph(40);
   tchiwg_expect->SetName("Graph");
   tchiwg_expect->SetTitle("Graph");
   tchiwg_expect->SetFillColor(3);
   TGraph *tchiwg_observe = new TGraph(40);
   tchiwg_observe->SetName("Graph");
   tchiwg_observe->SetTitle("Graph");
   tchiwg_observe->SetFillColor(kBlack);

   TGraphAsymmErrors *gre = new TGraphAsymmErrors(40);
   gre->SetName("sigmaone");
   gre->SetTitle("sigmaone");
   gre->SetFillColor(kGreen);
   gre->SetLineColor(kGreen);
   gre->SetLineWidth(0);
   gre->SetMarkerStyle(20);
   TGraphAsymmErrors *gre2 = new TGraphAsymmErrors(40);
   gre2->SetName("sigma2");
   gre2->SetTitle("sigma2");
   gre2->SetFillColor(5);
   gre2->SetLineColor(5);
   gre2->SetLineWidth(0);
   gre2->SetMarkerStyle(20);

	 p_theory->SetLineColor(kBlue);
	 p_theoryup->SetLineColor(kBlue);
	 p_theorydo->SetLineColor(kBlue);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000cc");
   tchiwg_expect->SetLineColor(kBlack);
   tchiwg_expect->SetLineStyle(2);
	 tchiwg_expect->SetFillColor(0);
   tchiwg_expect->SetFillStyle(0);
   tchiwg_expect->SetLineWidth(3);
   tchiwg_expect->SetMarkerStyle(20);
   tchiwg_observe->SetLineColor(kBlack);
   tchiwg_observe->SetLineStyle(1);
	 tchiwg_observe->SetFillColor(0);
   tchiwg_observe->SetFillStyle(0);
   tchiwg_observe->SetLineWidth(3);
   tchiwg_observe->SetMarkerStyle(20);

	 TFile xSecFile("susyxSec.root");
	 TH1F *p_crosssection = (TH1F*)xSecFile.Get("p_charginoSec");
   p_crosssection->GetXaxis()->SetTitle("M_{#tilde{#chi}} (GeV)");

   // ------------------------------------------------------------------
   // Read limits (this block replaces makelimit.sh + the old
   // limitvalue_Run2_*.txt reading).
   // ------------------------------------------------------------------
   std::vector<LimitPoint> expPts, obsPts;
   int n_exp_missing = 0, n_obs_missing = 0;
   for(int i = 0; i < kNMass; i++){
      const int mass = kMassStart + i * kMassStep;
      const std::string sm = std::to_string(mass);

      // ---- expected ----
      {
         const std::string path = kLimitDir + "limit_" + kExpTag + "_" + sm + ".txt";
         bool ok = false;
         std::vector<std::string> lines = readLines(path, ok);
         LimitPoint p; p.mass = mass;
         if(!ok){
            std::cerr << "WARNING: " << path << " not found -- skipping expected for mass " << mass << std::endl;
            n_exp_missing++;
         } else if(findLimit(lines, "50.0%", p.r)  && findLimit(lines, "84.0%", p.s1up) &&
                   findLimit(lines, "16.0%", p.s1lo) && findLimit(lines, "97.5%", p.s2up) &&
                   findLimit(lines, "2.5%",  p.s2lo)){
            expPts.push_back(p);
         } else {
            std::cerr << "WARNING: incomplete expected quantiles in " << path << " -- skipping mass " << mass << std::endl;
            n_exp_missing++;
         }
      }

      // ---- observed ----
      {
         const std::string path = kLimitDir + "limit_" + kObsTag + "_" + sm + ".txt";
         bool ok = false;
         std::vector<std::string> lines = readLines(path, ok);
         LimitPoint p; p.mass = mass; p.s1up = p.s1lo = p.s2up = p.s2lo = 0;
         if(!ok){
            std::cerr << "WARNING: " << path << " not found -- skipping observed for mass " << mass << std::endl;
            n_obs_missing++;
         } else if(findLimit(lines, "Observed Limit", p.r)){
            obsPts.push_back(p);
         } else {
            std::cerr << "WARNING: no 'Observed Limit' line parsed from " << path
                      << " -- was it accidentally run with --run blind?" << std::endl;
            n_obs_missing++;
         }
      }
   }
   std::cout << "Expected: " << expPts.size() << " mass point(s) (" << n_exp_missing << " missing/skipped)." << std::endl;
   std::cout << "Observed: " << obsPts.size() << " mass point(s) (" << n_obs_missing << " missing/skipped)." << std::endl;

   if(kWriteTextFiles){   // same output as makelimit.sh
      std::ofstream fexp("limitvalue_Run2_expected.txt"), fobs("limitvalue_Run2_observed.txt");
      for(size_t k = 0; k < expPts.size(); k++)
         fexp << expPts[k].mass << " " << expPts[k].r << " " << expPts[k].s1up << " " << expPts[k].s1lo
              << " " << expPts[k].s2up << " " << expPts[k].s2lo << "\n";
      for(size_t k = 0; k < obsPts.size(); k++)
         fobs << obsPts[k].mass << " " << obsPts[k].r << "\n";
   }

   // ---- fill graphs ----
   unsigned n_exp = 0;
   for(size_t k = 0; k < expPts.size(); k++){
      const LimitPoint& p = expPts[k];
      float theoryXS = p_crosssection->GetBinContent(p_crosssection->FindBin(p.mass));
      tchiwg_expect->SetPoint(n_exp, p.mass, p.r * theoryXS);
      gre->SetPoint(n_exp, p.mass, p.r * theoryXS);
      gre->SetPointError(n_exp, 0, 0, (p.r - p.s1lo) * theoryXS, (p.s1up - p.r) * theoryXS);
      gre2->SetPoint(n_exp, p.mass, p.r * theoryXS);
      gre2->SetPointError(n_exp, 0, 0, (p.r - p.s2lo) * theoryXS, (p.s2up - p.r) * theoryXS);
      n_exp++;
   }
   tchiwg_expect->Set(n_exp); gre->Set(n_exp); gre2->Set(n_exp);   // drop unused preallocated (0,0) points

   unsigned n_obs = 0;
   for(size_t k = 0; k < obsPts.size(); k++){
      float theoryXS = p_crosssection->GetBinContent(p_crosssection->FindBin(obsPts[k].mass));
      tchiwg_observe->SetPoint(n_obs, obsPts[k].mass, obsPts[k].r * theoryXS);
      n_obs++;
   }
   tchiwg_observe->Set(n_obs);

   // Theory curve + band over the full mass grid, independent of which limit
   // files exist. xserror[] is indexed by mass (i = (mass-300)/25), NOT by the
   // position in the limit file, so a skipped mass can no longer shift the band.
   for(int i = 0; i < kNMass; i++){
      const int mass = kMassStart + i * kMassStep;
      float theoryXS = p_crosssection->GetBinContent(p_crosssection->FindBin(mass));
      p_theory->SetPoint(i, mass, theoryXS);
      p_theoryup->SetPoint(i, mass, theoryXS + xserror[i]);
      p_theorydo->SetPoint(i, mass, theoryXS - xserror[i]);
   }

   // ---- exclusion masses ----
   const double exp_excl = interpolateExclusion(expPts, p_crosssection);
   const double obs_excl = interpolateExclusion(obsPts, p_crosssection);
   if(exp_excl > 0) std::cout << "Interpolated exclusion limit (expected = theory) is at mass = " << exp_excl << " GeV" << std::endl;
   else             std::cout << "No exclusion: expected limit never crosses theory cross section." << std::endl;
   if(obs_excl > 0) std::cout << "Interpolated exclusion limit (observed = theory) is at mass = " << obs_excl << " GeV" << std::endl;
   else             std::cout << "No exclusion: observed limit never crosses theory cross section." << std::endl;

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
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.0);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   tchiwg_expect->SetHistogram(Graph_Graph1);
   
   tchiwg_expect->Draw("c");

   TH1F *Graph_GraphObs = new TH1F("Graph_GraphObs","Graph",100,300,1300);
   Graph_GraphObs->SetMinimum(0.1);
   Graph_GraphObs->SetMaximum(1000);
   Graph_GraphObs->SetDirectory(0);
   Graph_GraphObs->SetStats(0);
   Graph_GraphObs->SetLineStyle(0);
   Graph_GraphObs->SetMarkerStyle(20);
   Graph_GraphObs->GetXaxis()->SetLabelFont(42);
   Graph_GraphObs->GetXaxis()->SetLabelOffset(0.007);
   Graph_GraphObs->GetXaxis()->SetLabelSize(0.05);
   Graph_GraphObs->GetXaxis()->SetTitleSize(0.06);
   Graph_GraphObs->GetXaxis()->SetTitleOffset(0.9);
   Graph_GraphObs->GetXaxis()->SetTitleFont(42);
   Graph_GraphObs->GetYaxis()->SetLabelFont(42);
   Graph_GraphObs->GetYaxis()->SetLabelOffset(0.007);
   Graph_GraphObs->GetYaxis()->SetLabelSize(0.05);
   Graph_GraphObs->GetYaxis()->SetTitleSize(0.06);
   Graph_GraphObs->GetYaxis()->SetTitleOffset(1.0);
   Graph_GraphObs->GetYaxis()->SetTitleFont(42);
   Graph_GraphObs->GetZaxis()->SetLabelFont(42);
   Graph_GraphObs->GetZaxis()->SetLabelOffset(0.007);
   Graph_GraphObs->GetZaxis()->SetLabelSize(0.05);
   Graph_GraphObs->GetZaxis()->SetTitleSize(0.06);
   Graph_GraphObs->GetZaxis()->SetTitleFont(42);
   tchiwg_observe->SetHistogram(Graph_GraphObs);
   
   tchiwg_observe->Draw("c");
  
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
   Graph_GraphA->GetYaxis()->SetTitleOffset(1.0);
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
   Graph_GraphB->GetYaxis()->SetTitleOffset(1.0);
   Graph_GraphB->GetYaxis()->SetTitleFont(42);
   Graph_GraphB->GetZaxis()->SetLabelFont(42);
   Graph_GraphB->GetZaxis()->SetLabelOffset(0.007);
   Graph_GraphB->GetZaxis()->SetLabelSize(0.05);
   Graph_GraphB->GetZaxis()->SetTitleSize(0.06);
   Graph_GraphB->GetZaxis()->SetTitleFont(42);
   gre2->SetHistogram(Graph_GraphB);

   gre2->Draw("3");
   gre->Draw("3");
   tchiwg_expect->Draw("c");
   tchiwg_observe->Draw("c");
   p_theory->Draw("c");
   p_theoryup->Draw("c same");
   p_theorydo->Draw("c same");

   // Redraw the axes on top of the bands/curves so the frame, and ticks on
   // the top (and right) edges, are not hidden by the filled bands.
   pad2->cd();
   pad2->SetTickx(1);
   pad2->SetTicky(1);
   pad2->RedrawAxis();
   pad2->Modified();
   pad2->Update();

	 thecanvas->cd();
	 TPad *pad1 = new TPad("pad1", "pad1", 0, 0.7, 1, 1.0);
   pad1->SetTopMargin(0.3);
	 pad1->SetLeftMargin(0.12);
	 pad1->SetBottomMargin(0); 
	 pad1->Draw();  
	 pad1->cd(); 
	 TH1F *upperframe = new TH1F("upperframe","",1,0,1);
   upperframe->GetXaxis()->SetLabelSize(0);
	 upperframe->GetXaxis()->SetTickLength(0);
	 upperframe->GetXaxis()->SetTickSize(0);
	 upperframe->GetXaxis()->SetNdivisions(1);
   upperframe->GetYaxis()->SetLabelSize(0);
	 upperframe->GetYaxis()->SetTickLength(0);
   upperframe->GetYaxis()->SetTickSize(0);
	 upperframe->GetYaxis()->SetTicks("U");
	 upperframe->GetYaxis()->SetNdivisions(1);
	 upperframe->Draw();
	 CMS_lumi(pad1, 6, 1, 11);

/*	 TLatex *prelim_cnt = new TLatex(0.12,0.8, "#font[22]{CMS}");
	 TLatex *lumitex_cnt = new TLatex(0.57,0.8, "138 fb^{-1} (13 TeV)");
   prelim_cnt->SetTextSize(0.15);
	 lumitex_cnt->SetTextSize(0.15);
   prelim_cnt->SetNDC();
	 lumitex_cnt->SetNDC();
	 prelim_cnt->Draw();
	 lumitex_cnt->Draw();*/
	 TLatex *process = new TLatex(0.15, 0.5, "pp #rightarrow #tilde{#chi}^{#pm}#tilde{#chi}^{0}, #tilde{#chi}^{#pm} #rightarrow W^{#pm}#tilde{G}, #tilde{#chi}^{0} #rightarrow #gamma#tilde{G}");
   process->SetTextSize(0.12);
   process->SetNDC();
   process->Draw();

   TGraph *theo = new TGraph(2);
   theo->SetName("theo");
   theo->SetTitle("theo");
   theo->SetLineColor(kBlue);
   theo->SetLineStyle(1);
   theo->SetLineWidth(4);
   theo->SetMarkerStyle(20);
   theo->SetPoint(0, 0.5, 0.5 );
   theo->SetPoint(1, 0.6, 0.5 );

   TGraph *theoP = new TGraph(2);
   theoP->SetName("theoP");
   theoP->SetTitle("theoP");
   theoP->SetLineColor(kBlue);
   theoP->SetLineStyle(2);
   theoP->SetLineWidth(2);
   theoP->SetMarkerStyle(20);
   theoP->SetPoint(0, 0.5, 0.55 );
   theoP->SetPoint(1, 0.6, 0.55 );

   TGraph *theoM = new TGraph(2);
   theoM->SetName("theoM");
   theoM->SetTitle("theoM");
   theoM->SetLineColor(kBlue);
   theoM->SetLineStyle(2);
   theoM->SetLineWidth(2);
   theoM->SetMarkerStyle(20);
   theoM->SetPoint(0,0.5, 0.45 );
   theoM->SetPoint(1,0.6, 0.45 );

   TGraphErrors *expOne = new TGraphErrors(3);
   expOne->SetName("expOne");
   expOne->SetTitle("expOne");
   expOne->SetLineColor(1);
   expOne->SetLineStyle(2);
   expOne->SetLineWidth(4);
	 expOne->SetFillColor(kGreen);
	 expOne->SetFillStyle(1001);
   expOne->SetMarkerSize(0);
   expOne->SetPoint(0, 0.05,  0.5 );
   expOne->SetPoint(1, 0.1,  0.5 );
   expOne->SetPoint(2, 0.15,  0.5 );
   expOne->SetPointError(0, 0, 0.05 );
   expOne->SetPointError(1, 0.05, 0.05 );
   expOne->SetPointError(2, 0, 0.05 );

   TGraphErrors *expTwo = new TGraphErrors(3);
   expTwo->SetName("expTwo");
   expTwo->SetTitle("expTwo");
   expTwo->SetLineColor(1);
   expTwo->SetLineStyle(2);
   expTwo->SetLineWidth(4);
	 expTwo->SetFillColor(kYellow);
	 expTwo->SetFillStyle(1001);
   expTwo->SetMarkerSize(0);
   expTwo->SetPoint(0, 0.05,  0.25 );
   expTwo->SetPoint(1, 0.1,  0.25 );
   expTwo->SetPoint(2, 0.15,  0.25 );
   expTwo->SetPointError(0, 0, 0.05 );
   expTwo->SetPointError(1, 0.05, 0.05 );
   expTwo->SetPointError(2, 0, 0.05 );

   TGraph *obsL = new TGraph(2);
   obsL->SetName("obsL");
   obsL->SetTitle("obsL");
   obsL->SetLineColor(kBlack);
   obsL->SetLineStyle(1);
   obsL->SetLineWidth(4);
   obsL->SetMarkerStyle(20);
   obsL->SetPoint(0, 0.5, 0.25 );
   obsL->SetPoint(1, 0.6, 0.25 );

	 theo->Draw("same");
	 theoM->Draw("same");				
	 theoP->Draw("same");				
   expOne->Draw("LE2 same");
	 expTwo->Draw("LE2 same");
	 obsL->Draw("same");

	 
	 TLatex *leg_exp1 = new TLatex(0.17,0.45, "Expected #pm 1#sigma");
	 leg_exp1->SetTextSize(0.12);
   leg_exp1->Draw(); 
	 TLatex *leg_exp2 = new TLatex(0.17,0.2, "Expected #pm 2#sigma");
	 leg_exp2->SetTextSize(0.12);
   leg_exp2->Draw(); 
	 TLatex *leg_theo = new TLatex(0.62,0.45, "NLO #pm 1#sigma_{theory}");
	 leg_theo->SetTextSize(0.12);
   leg_theo->Draw(); 
	 TLatex *leg_obs = new TLatex(0.62,0.2, "Observed");
	 leg_obs->SetTextSize(0.12);
   leg_obs->Draw(); 
	 thecanvas->SaveAs("/uscms_data/d3/tmishra/Output/TChiWG_Run2_2016-17-2018.pdf");
}
