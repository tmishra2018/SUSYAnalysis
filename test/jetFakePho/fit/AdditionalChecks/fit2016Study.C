// Fit the jet->gamma numerator and denominator with the function and method of the 2016 analysis
// (/uscms/homes/m/mengleis/work/SUSY2016/SUSYAnalysis/test/jetFakePho/fit/fitJetFunc.C) and compare
// Data/Fit with the current fits of fitJetFunc_pT_nVtx.C (exp2c, parameters from the transfer-factor txt).
//   2016 function: c1*exp(l1*pt) + c2*exp(l2*pt), integrated over a 1-GeV bin (tmpjetfake_func)
//   2016 method:   linear chi2 on new_fakesPhoEt / new_controlPhoEt over 35-300 GeV, no parameter limits;
//                  den starts from an expo fit, num from fixed per-channel values
//   "2016, 35-600": same function and method, fit range extended to 35-600 GeV
// usage: root -b -q 'fit2016Study.C+(2,2018,1)'
#include <TFile.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

namespace f16 {
double exp_safe(double x) {
	if (x > 700) return std::exp(700);
	if (x < -700) return 0.0;
	return std::exp(x);
}
double exp2c_func(double *x, double *par) {
	double pt = x[0];
	return par[0]*exp_safe(-par[2]*pt) + par[1]*exp_safe(-par[3]*(pt - par[4] - 200)) + par[5];
}
// tmpjetfake_func of the 2016 fitJetFunc.C, REBINSIZE = 1
double tmpjetfake_func(double *x, double *par) {
	double lo = x[0] - 0.5, hi = x[0] + 0.5;
	double c1 = par[0], c2 = par[1], l1 = par[2], l2 = par[3];
	return (c1*std::exp(l1*hi)/l1 + c2*std::exp(l2*hi)/l2) - (c1*std::exp(l1*lo)/l1 + c2*std::exp(l2*lo)/l2);
}

struct Metrics { double chi2[3]; int n[3]; double ratio[4]; };

Metrics evaluate(const TGraphAsymmErrors *g, TF1 *f) {
	Metrics m{};
	const double ch[4] = {35, 110, 300, 600.1};
	const double rg[5] = {35, 110, 220, 300, 600.1};
	double dat[4] = {0}, fit[4] = {0};
	for (int i = 0; i < g->GetN(); i++) {
		double x = g->GetX()[i], y = g->GetY()[i], ey = g->GetErrorYhigh(i);
		if (y <= 0 || ey <= 0 || x < 35 || x > 595) continue;   // x = 600 is the den overflow point
		double lo = x - g->GetErrorXlow(i), hi = x + g->GetErrorXhigh(i);
		for (int k = 0; k < 3; k++)
			if (x >= ch[k] && x < ch[k+1]) { m.chi2[k] += std::pow((y - f->Eval(x))/ey, 2); m.n[k]++; }
		for (int k = 0; k < 4; k++)
			if (x >= rg[k] && x < rg[k+1]) { dat[k] += y*(hi - lo); fit[k] += f->Integral(lo, hi); }
	}
	for (int k = 0; k < 4; k++) m.ratio[k] = fit[k] > 0 ? dat[k]/fit[k] : 0;
	return m;
}

// den_* or num_* exp2c parameters from the transfer-factor txt
TF1 *currentFit(const std::string &txt, const std::string &tag, const char *name) {
	std::ifstream in(txt);
	std::map<std::string, double> p;
	std::string key; double v;
	while (in >> key >> v) p[key] = v;
	TF1 *f = new TF1(name, exp2c_func, 20, 1000, 6);
	const char *k[6] = {"coeff1", "coeff2", "lambd1", "lambd2", "shift", "offset"};
	for (int i = 0; i < 6; i++) f->SetParameter(i, p[tag + "_" + k[i]]);
	return f;
}

TF1 *fit2016(TGraphAsymmErrors *g, bool isNum, int channel, double xmax, const char *name) {
	TF1 *f = new TF1(name, tmpjetfake_func, 35, xmax, 4);
	if (isNum) {
		if (channel == 1) f->SetParameters(14300, 300, -0.070, -0.04);
		else              f->SetParameters(23900, 300, -0.082, -0.01);
	} else {
		TF1 *e = new TF1(Form("%s_expo", name), "expo", 35, xmax);
		g->Fit(e, "Q R N");
		double l1 = e->GetParameter(1), c = std::exp(e->GetParameter(0))/2;
		double l2 = (std::log(g->Eval(111) - 2*c*std::exp(l1*111.0)) - std::log(g->Eval(250) - 2*c*std::exp(l1*250.0)))/(111.0 - 250.0);
		if (l2 > 0 && l2 < 1e6) f->SetParameters(c, c/10, l1, l2);
		else                    f->SetParameters(c, c/10, l1, l1);
	}
	TFitResultPtr r = g->Fit(f, "Q R S N");
	std::cout << "[FIT] " << name << " status " << int(r) << ", chi2/ndf " << f->GetChisquare() << "/" << f->GetNDF() << ", pars";
	for (int i = 0; i < 4; i++) std::cout << " " << f->GetParameter(i);
	std::cout << std::endl;
	f->SetRange(20, 1000);
	return f;
}
}

void fit2016Study(int ichannel, int RunYear, bool ISpreVFP) {
	using namespace f16;
	gROOT->SetBatch(true);
	std::string whichVFP = (RunYear == 2016) ? (ISpreVFP ? "preVFP" : "postVFP") : "";
	const char *chName = ichannel == 1 ? "DoubleEG" : "MuonEG";
	TFile *fin = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_%s_%d%s_pT_nVtx.root", chName, RunYear, whichVFP.c_str()));
	std::string txt = Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-%s-EB_pT_nVtx.txt", RunYear, whichVFP.c_str(), chName);

	TGraphAsymmErrors *gNum = (TGraphAsymmErrors*)fin->Get("new_fakesPhoEt");
	TGraphAsymmErrors *gDen = (TGraphAsymmErrors*)fin->Get("new_controlPhoEt");
	std::vector<const char*> labels = {"current", "2016, 35-300", "2016, 35-600"};
	std::vector<TF1*> num = {currentFit(txt, "num", "num_cur"), fit2016(gNum, true, ichannel, 300, "num_16"), fit2016(gNum, true, ichannel, 600, "num_16ext")};
	std::vector<TF1*> den = {currentFit(txt, "den", "den_cur"), fit2016(gDen, false, ichannel, 300, "den_16"), fit2016(gDen, false, ichannel, 600, "den_16ext")};

	printf("[TABLE] %s %d%s\n", chName, RunYear, whichVFP.c_str());
	printf("[TABLE] %-18s %13s %12s %12s | %9s %9s %9s %9s | %7s %7s %7s\n", "", "chi2/n 35-110", "110-300", "300-600",
	       "D/F35-110", "110-220", "220-300", "300-600", "TF@60", "TF@150", "TF@250");
	for (int isNum = 1; isNum >= 0; isNum--)
		for (size_t i = 0; i < labels.size(); i++) {
			TF1 *f = isNum ? num[i] : den[i];
			Metrics m = evaluate(isNum ? gNum : gDen, f);
			printf("[TABLE] %s %-14s %8.1f/%-3d %8.1f/%-3d %8.1f/%-3d | %9.3f %9.3f %9.3f %9.3f | %7.4f %7.4f %7.4f\n",
			       isNum ? "num" : "den", labels[i], m.chi2[0], m.n[0], m.chi2[1], m.n[1], m.chi2[2], m.n[2],
			       m.ratio[0], m.ratio[1], m.ratio[2], m.ratio[3],
			       num[i]->Eval(60)/den[i]->Eval(60), num[i]->Eval(150)/den[i]->Eval(150), num[i]->Eval(250)/den[i]->Eval(250));
		}

	// Data/Fit for the current and 2016 (35-300) fits, num in the top half, den in the bottom half
	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c", "", 800, 800);
	c->Divide(1, 2);
	for (int isNum = 1; isNum >= 0; isNum--) {
		c->cd(2 - isNum);
		gPad->SetGridy();
		TH1D *frame = new TH1D(Form("frame%d", isNum), Form("%s %d%s %s;p_{T} (GeV);Data / Fit", chName, RunYear, whichVFP.c_str(),
		                       isNum ? "fake photons" : "hadron proxies"), 1, 35, 600);
		frame->SetMinimum(0); frame->SetMaximum(2.5);
		frame->Draw("AXIS");
		TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.88);
		const int col[3] = {kRed, kBlue, kGreen+2};
		TGraphAsymmErrors *g = isNum ? gNum : gDen;
		for (size_t k = 0; k < labels.size(); k++) {
			TF1 *f = isNum ? num[k] : den[k];
			TGraphAsymmErrors *r = new TGraphAsymmErrors();
			int np = 0;
			for (int i = 0; i < g->GetN(); i++) {
				double x = g->GetX()[i], y = g->GetY()[i];
				if (y <= 0 || x < 35 || x > 595) continue;
				double fv = f->Eval(x);
				r->SetPoint(np, x + 0.8*k, y/fv);
				r->SetPointError(np, 0, 0, g->GetErrorYlow(i)/fv, g->GetErrorYhigh(i)/fv);
				np++;
			}
			r->SetMarkerStyle(20); r->SetMarkerSize(0.6);
			r->SetMarkerColor(col[k]); r->SetLineColor(col[k]);
			r->Draw("P same");
			leg->AddEntry(r, labels[k], "p");
		}
		TLine *one = new TLine(35, 1, 600, 1);
		one->SetLineStyle(2);
		one->Draw();
		leg->Draw();
	}
	c->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/numFitStudy/fit2016Study_%s_%d%s.png", chName, RunYear, whichVFP.c_str()));
}
