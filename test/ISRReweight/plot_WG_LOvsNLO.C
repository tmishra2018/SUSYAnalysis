// Plots the outputs of compare_WG_LOvsNLO.C (merged by run_WGLOvsNLO.sh):
// unit-area shapes of LO W+gamma, LO W+gamma with the ISR weights (data / NLO Z+gamma,
// ISRweights_<era>_<egamma|mgamma>.txt), LO W+gamma with the LO->NLO factors r, and NLO W+gamma,
// with the ratios to LO. The weighted LO curves are filled event by event in compare_WG_LOvsNLO.C
// (<var>_<ch>_<region>_LOw / _LOr), so they are valid for every variable, not only ISR pT.
// The ISR systematic is the full with/without-weights difference, so NLO/LO is covered
// if it lies within |w - 1| of unity.
//   root -l -b -q 'plot_WG_LOvsNLO.C+'
// Writes plots/WGLOvsNLO_<var>_<ch>_<region>_<era>.{png,pdf} and a bin table to stdout.

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

const char* kDir = "/eos/uscms/store/user/tmishra/ISRweighting/WGLOvsNLO";
const char* kEras[4] = {"2016preVFP", "2016postVFP", "2017", "2018"};
// lumi per era and channel (DoubleEG for eg, MuonEG for mg), as include/analysis_tools.h
const double kLumi[2][4] = {{19.499256, 16.810813, 41.4613, 59.816229},
                            {19.495077, 16.810813, 27.100211, 59.810698}};
const double kISRedges[8] = {0,50,100,150,200,250,300,1e9};

// sample histogram normalised to its own event count (bin 1 of h_count)
static TH1D* getNorm(const char* sample, const char* era, const char* hname){
	TFile* f = TFile::Open(Form("%s/WGLOvsNLO_%s_%s.root", kDir, sample, era));
	if(!f || f->IsZombie()){ std::cerr << "missing " << sample << " " << era << std::endl; return 0; }
	TH1D* h = (TH1D*)f->Get(hname)->Clone(Form("%s_%s_%s", hname, sample, era));
	h->SetDirectory(0);
	h->Scale(1.0/((TH1D*)f->Get("h_count"))->GetBinContent(1));
	f->Close();
	return h;
}

// LO (three stitched samples), LO x ISR weight w, LO x LO->NLO factor r, NLO
// for one era (iera < 4) or Run 2 (iera = 4)
static TH1D* getLO(const char* hname, int e){
	TH1D* l = 0;
	for(const char* s : {"WGToLNuG", "WGJet40", "WGJet130"}){
		TH1D* h = getNorm(s, kEras[e], hname);
		if(!l) l = h; else { l->Add(h); delete h; }
	}
	return l;
}

static void getHists(const char* var, int ch, const char* reg, int iera, TH1D*& lo, TH1D*& loISR, TH1D*& nlo, TH1D*& loR){
	lo = loISR = nlo = loR = 0;
	const char* chn = ch == 0 ? "eg" : "mg";
	TString hname = Form("%s_%s_%s", var, chn, reg);
	for(int e = 0; e < 4; e++){
		if(iera < 4 && e != iera)continue;
		TH1D* l  = getLO(hname, e);
		TH1D* li = getLO(hname + "_LOw", e);
		TH1D* lr = getLO(hname + "_LOr", e);
		TH1D* n  = getNorm("NLO", kEras[e], hname);
		// unit area per era before summing with lumi, so each era enters with its lumi share
		for(TH1D* h : {l, li, lr, n}) h->Scale(kLumi[ch][e]/h->Integral());
		if(!lo){ lo = l; loISR = li; nlo = n; loR = lr; }
		else { lo->Add(l); loISR->Add(li); nlo->Add(n); loR->Add(lr); delete l; delete li; delete n; delete lr; }
	}
	for(TH1D* h : {lo, loISR, nlo, loR}) h->Scale(1.0/h->Integral());
}

static void drawOne(const char* var, int ch, const char* reg, int iera){
	TH1D *lo, *loISR, *nlo, *loR;
	getHists(var, ch, reg, iera, lo, loISR, nlo, loR);
	const char* eraName = iera < 4 ? kEras[iera] : "Run2";

	TCanvas* c = new TCanvas("c", "", 600, 650);
	TPad* p1 = new TPad("p1", "", 0, 0.33, 1, 1);
	TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.33);
	p1->SetBottomMargin(0.02); p2->SetTopMargin(0.03); p2->SetBottomMargin(0.32);
	p1->SetLogy(); p1->Draw(); p2->Draw();

	p1->cd();
	lo->SetLineColor(kBlue+1);     lo->SetMarkerColor(kBlue+1);     lo->SetMarkerStyle(24);
	loISR->SetLineColor(kRed+1);   loISR->SetMarkerColor(kRed+1);   loISR->SetMarkerStyle(20);
	nlo->SetLineColor(kBlack);     nlo->SetMarkerColor(kBlack);     nlo->SetMarkerStyle(21);
	loR->SetLineColor(kGreen+2);   loR->SetMarkerColor(kGreen+2);   loR->SetMarkerStyle(26);
	lo->SetTitle(""); lo->GetYaxis()->SetTitle("normalised to unit area");
	lo->GetXaxis()->SetLabelSize(0);
	lo->SetMaximum(2*std::max(lo->GetMaximum(), nlo->GetMaximum()));
	lo->SetMinimum(1e-5);
	lo->Draw("E1"); loISR->Draw("E1 same"); loR->Draw("E1 same"); nlo->Draw("E1 same");
	TLegend* leg = new TLegend(0.45, 0.64, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->AddEntry(lo,    "LO W#gamma (MadGraph MLM)", "lp");
	leg->AddEntry(loISR, "LO W#gamma #times ISR weight", "lp");
	leg->AddEntry(loR,   "LO W#gamma #times LO#rightarrowNLO factor", "lp");
	leg->AddEntry(nlo,   "NLO W#gamma (aMC@NLO FxFx)", "lp");
	leg->Draw();
	TLatex t; t.SetNDC(); t.SetTextSize(0.045);
	t.DrawLatex(0.13, 0.92, Form("%s, %s, %s", ch == 0 ? "e#gamma" : "#mu#gamma", reg, eraName));

	p2->cd();
	TH1D* rNLO = (TH1D*)nlo->Clone("rNLO"); rNLO->Divide(lo);
	TH1D* rISR = (TH1D*)loISR->Clone("rISR"); rISR->Divide(lo);
	TH1D* rR = (TH1D*)loR->Clone("rR"); rR->Divide(lo);
	for(int b = 1; b <= rISR->GetNbinsX(); b++){ rISR->SetBinError(b, 0); rR->SetBinError(b, 0); }
	rNLO->SetTitle("");
	rNLO->GetYaxis()->SetTitle("ratio to LO");
	rNLO->GetYaxis()->SetRangeUser(0.4, 1.6);
	rNLO->GetYaxis()->SetNdivisions(505);
	rNLO->GetYaxis()->SetTitleSize(0.10); rNLO->GetYaxis()->SetTitleOffset(0.45); rNLO->GetYaxis()->SetLabelSize(0.09);
	rNLO->GetXaxis()->SetTitleSize(0.11); rNLO->GetXaxis()->SetLabelSize(0.09);
	rNLO->Draw("E1");
	rISR->SetLineWidth(2);
	rISR->Draw("hist same");
	rR->SetLineWidth(2);
	rR->Draw("hist same");
	TLine* one = new TLine(rNLO->GetXaxis()->GetXmin(), 1, rNLO->GetXaxis()->GetXmax(), 1);
	one->SetLineStyle(2); one->Draw();

	gSystem->mkdir(Form("%s/plots", kDir), true);
	for(const char* ext : {"png", "pdf"})
		c->SaveAs(Form("%s/plots/WGLOvsNLO_%s_%s_%s_%s.%s", kDir, var, ch == 0 ? "eg" : "mg", reg, eraName, ext));

	if(TString(var) == "ISRJetPt" || TString(var) == "ISRJetPtReco"){
		printf("%-6s %-3s %-6s | %-9s | %7s %7s %7s | %s\n", eraName, ch == 0 ? "eg" : "mg", reg, "ISR pT", "NLO/LO", "stat", "ISRw", "covered");
		for(int b = 1; b <= rNLO->GetNbinsX(); b++){
			const double r = rNLO->GetBinContent(b), e = rNLO->GetBinError(b), w = rISR->GetBinContent(b);
			const bool covered = fabs(r - 1) <= fabs(w - 1) + e;
			printf("%-6s %-3s %-6s | %3.0f-%-5.0f | %7.3f %7.3f %7.3f | %s\n", eraName, ch == 0 ? "eg" : "mg", reg,
			       rNLO->GetBinLowEdge(b), rNLO->GetBinLowEdge(b+1), r, e, w, covered ? "yes" : "NO");
		}
	}
	delete c;
	for(TH1D* h : {lo, loISR, nlo, loR, rNLO, rISR, rR}) delete h;
}

// closure of r: NLO / (LO x r) and NLO / LO in coarse bins, Run 2, for variables other than ISR pT
//   root -l -b -q 'plot_WG_LOvsNLO.C+' -e 'closureLOtoNLO()'
void closureLOtoNLO(){
	struct V { const char* name; std::vector<double> edges; };
	const std::vector<V> vars = {{"phoEt", {35, 100, 200, 300, 500, 1200}}, {"HT", {0, 100, 200, 400, 600, 1500}},
	                             {"nJet", {-0.5, 0.5, 1.5, 2.5, 3.5, 9.5}}, {"MET", {0, 70, 120, 200, 400, 900}}};
	for(int ch = 0; ch < 2; ch++)
		for(const char* reg : {"presel", "CR", "SR"})
			for(const V& v : vars){
				TH1D *lo, *loISR, *nlo, *loR;
				getHists(v.name, ch, reg, 4, lo, loISR, nlo, loR);
				printf("%-2s %-6s %-6s", ch == 0 ? "eg" : "mg", reg, v.name);
				for(size_t i = 0; i + 1 < v.edges.size(); i++){
					const int b1 = lo->FindBin(v.edges[i] + 1e-3), b2 = lo->FindBin(v.edges[i+1] - 1e-3);
					double en, er;
					const double n = nlo->IntegralAndError(b1, b2, en), l = lo->Integral(b1, b2), r = loR->IntegralAndError(b1, b2, er);
					printf(" | %g-%g: NLO/LO %.3f NLO/LOr %.3f+-%.3f", v.edges[i], v.edges[i+1], n/l, n/r, (n/r)*sqrt(pow(en/n, 2) + pow(er/r, 2)));
				}
				printf("\n");
				for(TH1D* h : {lo, loISR, nlo, loR}) delete h;
			}
}

void plot_WG_LOvsNLO(){
	gStyle->SetOptStat(0);
	for(int ch = 0; ch < 2; ch++)
		for(const char* reg : {"presel", "CR", "SR"})
			for(int iera = 0; iera < 5; iera++)
				for(const char* var : {"ISRJetPt", "ISRJetPtReco", "ISRJetPtFine", "nJet", "HT", "phoEt", "MET", "bosonPt"}){
					if(iera < 4 && !TString(var).BeginsWith("ISRJetPt"))continue;  // per-era plots only for ISR pT
					drawOne(var, ch, reg, iera);
				}
}

// LO -> NLO W+gamma ISR pT factors r = NLO/LO (unit-area shapes, Run 2 combined, preselection),
// written to /eos/uscms/store/user/tmishra/ISRweighting/WGLOvsNLO_factors_<egamma|mgamma>.txt
// in the format of the ISRweights files. The normalization line is filled when the factors are
// applied (keeps the W+gamma yield fixed), so it is written as 1 here.
// r is taken in ISRJetPtReco (jets with dR > 0.4 to the photon and lepton), the ISR pT of the derivation
// (analysis_ISRMC.C) and of the V+gamma skim since 2026-10-29; until then it was taken in the
// generator-matched ISRJetPt, which also counted the photon's own jet.
//   root -l -b -q -e '.L plot_WG_LOvsNLO.C+' -e 'writeLOtoNLOFactors()'
void writeLOtoNLOFactors(){
	for(int ch = 0; ch < 2; ch++){
		TH1D *lo, *loISR, *nlo, *loR;
		getHists("ISRJetPtReco", ch, "presel", 4, lo, loISR, nlo, loR);
		TH1D *loSR, *loISRSR, *nloSR, *loRSR;
		getHists("ISRJetPtReco", ch, "SR", 4, loSR, loISRSR, nloSR, loRSR);
		const char* chn = ch == 0 ? "egamma" : "mgamma";
		std::ofstream out(Form("/eos/uscms/store/user/tmishra/ISRweighting/WGLOvsNLO_factors_%s.txt", chn));
		out << "# LO -> NLO W+gamma ISR pT factors (NLO/LO in ISRJetPtReco, unit area, Run 2, preselection), from compare_WG_LOvsNLO.C\n";
		out << "# normalization 1\n";
		out << "# ptLow ptHigh factor [stat SR_factor SR_stat]\n";
		for(int b = 1; b <= lo->GetNbinsX(); b++){
			const double r  = nlo->GetBinContent(b)/lo->GetBinContent(b);
			const double e  = r*sqrt(pow(nlo->GetBinError(b)/nlo->GetBinContent(b), 2) + pow(lo->GetBinError(b)/lo->GetBinContent(b), 2));
			const double rs = nloSR->GetBinContent(b)/loSR->GetBinContent(b);
			const double es = rs*sqrt(pow(nloSR->GetBinError(b)/nloSR->GetBinContent(b), 2) + pow(loSR->GetBinError(b)/loSR->GetBinContent(b), 2));
			out << kISRedges[b-1] << " " << kISRedges[b] << " " << r << " " << e << " " << rs << " " << es << "\n";
			printf("%s %3.0f-%-5.0f r = %.3f +- %.3f   (SR %.3f +- %.3f, pull %.1f)\n", chn, kISRedges[b-1], b < 7 ? kISRedges[b] : 800., r, e, rs, es, (rs - r)/sqrt(e*e + es*es));
		}
	}
}
