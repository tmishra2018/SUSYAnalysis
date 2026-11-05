// Compare numerator (fake-photon yield) fit variants for fitJetFunc_pT_nVtx.C.
// Reads new_fakesPhoEt from fitJetFunc_Output_<ch>_<era>_pT_nVtx.root, refits exp2c_func
// several ways and prints Data/Fit on the ratio-pad points for each.
//   A: log10 chi2 on all points (what fitJetFunc_pT_nVtx.C does now)
//   B: linear chi2 on all points, started from A
//   C: log10 chi2 on points merged to 10 GeV everywhere (1-GeV points below 110 merged)
//   D: log10 chi2 on points merged into the jet fake-rate pT bins below 220 GeV, 10 GeV above
//   E: as D below 110 GeV only (the 1-GeV points), original 10-GeV points above
// Merged errors are summed linearly: below 110 GeV they are dominated by the jet fake-rate
// error, which is the same number for every point of one fake-rate bin.
// usage: root -b -q 'numFitStudy.C+(2,2018,1)'
#include <TFile.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <cmath>

namespace nfs {
double exp_safe(double x) {
	if (x > 700) return std::exp(700);
	if (x < -700) return 0.0;
	return std::exp(x);
}
double exp2c_func(double *x, double *par) {
	double pt = x[0];
	return par[0]*exp_safe(-par[2]*pt) + par[1]*exp_safe(-par[3]*(pt - par[4] - 200)) + par[5];
}
double log_exp2c_func(double *x, double *par) {
	double v = exp2c_func(x, par);
	if (!(v > 0)) v = 1e-12;
	return std::log10(v);
}

void setLimits(TF1 *f) {
	f->SetParLimits(0, 0, 1e9);
	f->SetParLimits(1, 0, 1e5);
	f->SetParLimits(2, 0.0001, 1);
	f->SetParLimits(3, 0.0001, 1);
	f->SetParLimits(4, 0, 300);
	f->SetParLimits(5, 0, 100);
}

TGraphAsymmErrors *toLog(const TGraphAsymmErrors *g) {
	TGraphAsymmErrors *lg = new TGraphAsymmErrors();
	int np = 0;
	for (int i = 0; i < g->GetN(); i++) {
		double x = g->GetX()[i], y = g->GetY()[i];
		if (y <= 0) continue;
		lg->SetPoint(np, x, std::log10(y));
		lg->SetPointError(np, g->GetErrorXlow(i), g->GetErrorXhigh(i),
		                  g->GetErrorYlow(i)/(y*std::log(10.)), g->GetErrorYhigh(i)/(y*std::log(10.)));
		np++;
	}
	return lg;
}

// merge the points of g with x in [edges[k], edges[k+1]) into one per-GeV point
TGraphAsymmErrors *mergePoints(const TGraphAsymmErrors *g, const std::vector<double> &edges) {
	TGraphAsymmErrors *m = new TGraphAsymmErrors();
	int np = 0;
	for (size_t k = 0; k + 1 < edges.size(); k++) {
		double sumy = 0, sume = 0, w = 0;
		for (int i = 0; i < g->GetN(); i++) {
			double x = g->GetX()[i], y = g->GetY()[i];
			if (y <= 0 || x < edges[k] || x >= edges[k+1]) continue;
			double wi = g->GetErrorXlow(i) + g->GetErrorXhigh(i);
			sumy += y*wi; sume += g->GetErrorYhigh(i)*wi; w += wi;
		}
		if (w <= 0 || sumy <= 0) continue;
		double xc = 0.5*(edges[k] + edges[k+1]);
		m->SetPoint(np, xc, sumy/w);
		m->SetPointError(np, xc - edges[k], edges[k+1] - xc, sume/w, sume/w);
		np++;
	}
	return m;
}

TF1 *fitLog(TGraphAsymmErrors *g, const double *start, const char *name) {
	TGraphAsymmErrors *lg = toLog(g);
	TF1 *fl = new TF1(Form("%s_log", name), log_exp2c_func, 20, 1000, 6);
	fl->SetParameters(start);
	setLimits(fl);
	lg->Fit(fl, "Q S N");
	TFitResultPtr r = lg->Fit(fl, "Q R S N");
	TF1 *f = new TF1(name, exp2c_func, 20, 1000, 6);
	f->SetParameters(fl->GetParameters());
	std::cout << "[FIT] " << name << " status " << int(r) << ", chi2/ndf " << fl->GetChisquare() << "/" << fl->GetNDF() << ", pars";
	for (int i = 0; i < 6; i++) std::cout << " " << fl->GetParameter(i);
	std::cout << std::endl;
	return f;
}

struct Metrics { double chi2[3]; int n[3]; double ratio[4]; };

Metrics evaluate(const TGraphAsymmErrors *g, TF1 *f) {
	Metrics m{};
	const double ch[4] = {35, 110, 300, 600.1};
	const double rg[5] = {35, 110, 220, 300, 600.1};
	double dat[4] = {0}, fit[4] = {0};
	for (int i = 0; i < g->GetN(); i++) {
		double x = g->GetX()[i], y = g->GetY()[i], ey = g->GetErrorYhigh(i);
		if (y <= 0 || ey <= 0) continue;
		double lo = x - g->GetErrorXlow(i), hi = x + g->GetErrorXhigh(i);
		for (int k = 0; k < 3; k++)
			if (x >= ch[k] && x < ch[k+1]) { m.chi2[k] += std::pow((y - f->Eval(x))/ey, 2); m.n[k]++; }
		for (int k = 0; k < 4; k++)
			if (x >= rg[k] && x < rg[k+1]) { dat[k] += y*(hi - lo); fit[k] += f->Integral(lo, hi); }
	}
	for (int k = 0; k < 4; k++) m.ratio[k] = fit[k] > 0 ? dat[k]/fit[k] : 0;
	return m;
}
}

void numFitStudy(int ichannel, int RunYear, bool ISpreVFP) {
	using namespace nfs;
	gROOT->SetBatch(true);
	std::string whichVFP = (RunYear == 2016) ? (ISpreVFP ? "preVFP" : "postVFP") : "";
	const char *chName = ichannel == 1 ? "DoubleEG" : "MuonEG";
	TFile *fin = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_%s_%d%s_pT_nVtx.root", chName, RunYear, whichVFP.c_str()));
	TGraphAsymmErrors *g = (TGraphAsymmErrors*)fin->Get("new_fakesPhoEt");

	const double start[6] = {51521, 1.13974, 0.0560658, 0.0200371, 125.136, 0.0970124};
	std::vector<TF1*> fits;
	std::vector<const char*> labels = {"A log, all pts", "B linear", "C log, 10 GeV", "D log, FR bins", "E log, FR<110"};

	TF1 *fA = fitLog(g, start, "fitA");
	fits.push_back(fA);

	TF1 *fB = new TF1("fitB", exp2c_func, 35, 600, 6);
	fB->SetParameters(fA->GetParameters());
	setLimits(fB);
	TFitResultPtr rB = g->Fit(fB, "Q R S N");
	std::cout << "[FIT] fitB status " << int(rB) << ", chi2/ndf " << fB->GetChisquare() << "/" << fB->GetNDF() << std::endl;
	fits.push_back(fB);

	std::vector<double> e10;
	for (double x = 35; x < 110; x += 10) e10.push_back(x);
	for (double x = 110; x <= 600; x += 10) e10.push_back(x);
	fits.push_back(fitLog(mergePoints(g, e10), fA->GetParameters(), "fitC"));

	std::vector<double> eFR = {35, 36, 38, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150, 180, 220};
	for (double x = 230; x <= 600; x += 10) eFR.push_back(x);
	fits.push_back(fitLog(mergePoints(g, eFR), fA->GetParameters(), "fitD"));

	std::vector<double> eE = {35, 36, 38, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100};
	for (double x = 110; x <= 600; x += 10) eE.push_back(x);
	fits.push_back(fitLog(mergePoints(g, eE), fA->GetParameters(), "fitE"));

	printf("[TABLE] %s %d%s\n", chName, RunYear, whichVFP.c_str());
	printf("[TABLE] %-16s %12s %12s %12s | %9s %9s %9s %9s\n", "variant", "chi2/n 35-110", "110-300", "300-600",
	       "D/F35-110", "110-220", "220-300", "300-600");
	for (size_t i = 0; i < fits.size(); i++) {
		Metrics m = evaluate(g, fits[i]);
		printf("[TABLE] %-16s %8.1f/%-3d %8.1f/%-3d %8.1f/%-3d | %9.3f %9.3f %9.3f %9.3f\n", labels[i],
		       m.chi2[0], m.n[0], m.chi2[1], m.n[1], m.chi2[2], m.n[2], m.ratio[0], m.ratio[1], m.ratio[2], m.ratio[3]);
	}

	// ratio pads of all variants, 35-600 GeV
	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c", "", 800, 600);
	c->SetGridy();
	TH1D *frame = new TH1D("frame", Form("%s %d%s fake photons;p_{T} (GeV);Data / Fit", chName, RunYear, whichVFP.c_str()), 1, 35, 600);
	frame->SetMinimum(0); frame->SetMaximum(2.5);
	frame->Draw("AXIS");
	TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.88);
	const int col[5] = {kRed, kBlue, kGreen+2, kMagenta+1, kOrange+7};
	for (size_t k = 0; k < fits.size(); k++) {
		TGraphAsymmErrors *r = new TGraphAsymmErrors();
		int np = 0;
		for (int i = 0; i < g->GetN(); i++) {
			double x = g->GetX()[i], y = g->GetY()[i];
			if (y <= 0) continue;
			double fv = fits[k]->Eval(x);
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
	c->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/numFitStudy/numFitStudy_%s_%d%s.png", chName, RunYear, whichVFP.c_str()));
}
