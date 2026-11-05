// Estimate how merging jet fake-rate pT bins (FitJetFake_pT_nVtx.C) changes the numerator Data/Fit in
// fitJetFunc_pT_nVtx.C, without rerunning the template fits.
// The fake-photon points are (data - e->gamma) x fake rate of their pT bin. Each point is divided by its
// current fake rate and multiplied by the merged one: the yield-weighted mean of the merged bins, with
// the inverse-variance combined error (an approximation of what a merged template fit would give).
// The numerator is then refitted as in fitJetFunc_pT_nVtx.C (log10 chi2 of exp2c_func, same start values
// and limits) and Data/Fit is printed per pT region.
// usage: root -b -q 'binTuneStudy.C+(2,2018,1)'
#include <TFile.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TFitResult.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

namespace btune {
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

TF1 *fitNum(TGraphAsymmErrors *g, const char *name) {
	TGraphAsymmErrors *lg = new TGraphAsymmErrors();
	int np = 0;
	for (int i = 0; i < g->GetN(); i++) {
		double x = g->GetX()[i], y = g->GetY()[i];
		if (y <= 0) continue;
		lg->SetPoint(np, x, std::log10(y));
		lg->SetPointError(np, g->GetErrorXlow(i), g->GetErrorXhigh(i), g->GetErrorYlow(i)/(y*std::log(10.)), g->GetErrorYhigh(i)/(y*std::log(10.)));
		np++;
	}
	TF1 *fl = new TF1(Form("%s_log", name), log_exp2c_func, 20, 1000, 6);
	fl->SetParameters(51521, 1.13974, 0.0560658, 0.0200371, 125.136, 0.0970124);
	fl->SetParLimits(0, 0, 1e9);
	fl->SetParLimits(1, 0, 1e5);
	fl->SetParLimits(2, 0.0001, 1);
	fl->SetParLimits(3, 0.0001, 1);
	fl->SetParLimits(4, 0, 300);
	fl->SetParLimits(5, 0, 100);
	lg->Fit(fl, "Q S N");
	lg->Fit(fl, "Q R S N");
	TF1 *f = new TF1(name, exp2c_func, 20, 1000, 6);
	f->SetParameters(fl->GetParameters());
	return f;
}
}

void binTuneStudy(int ichannel, int RunYear, bool ISpreVFP) {
	using namespace btune;
	std::string whichVFP = (RunYear == 2016) ? (ISpreVFP ? "preVFP" : "postVFP") : "";
	const char *chName = ichannel == 1 ? "DoubleEG" : "MuonEG";
	TFile *fin = TFile::Open(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_%s_%d%s_pT_nVtx.root", chName, RunYear, whichVFP.c_str()));
	TGraphAsymmErrors *g0 = (TGraphAsymmErrors*)fin->Get("new_fakesPhoEt");

	// jet fake rates, as fitJetFunc_pT_nVtx.C reads them
	std::ifstream in(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-%s-EB_pT_nVtx.txt", RunYear, whichVFP.c_str(), chName));
	std::vector<double> lo, hi, fr, err;
	double a, b, c, d, e, f;
	while (in >> a >> b >> c >> d >> e >> f) { lo.push_back(a); hi.push_back(b); fr.push_back(c); err.push_back(d); }
	const int nfr = lo.size();
	auto frBin = [&](double x) { int k = 0; for (int i = 0; i < nfr; i++) if (x >= lo[i]) k = i; return k; };

	// (data - e->gamma) per point, and its yield per fake-rate bin
	std::vector<double> sig(g0->GetN(), 0), yield(nfr, 0);
	for (int i = 0; i < g0->GetN(); i++) {
		double x = g0->GetX()[i], y = g0->GetY()[i];
		if (y <= 0 || x < 35 || x >= 600) continue;
		int k = frBin(x);
		sig[i] = y/fr[k];
		yield[k] += sig[i]*(g0->GetErrorXlow(i) + g0->GetErrorXhigh(i));
	}

	// merging schemes: groups of [low, high) edges whose fake-rate bins are combined
	struct Scheme { const char *name; std::vector<std::pair<double,double>> merge; };
	std::vector<Scheme> schemes = {
		{"current",                {}},
		{"100-150",                {{100,150}}},
		{"150-300",                {{150,300}}},
		{"100-150, 150-300",       {{100,150},{150,300}}},
		{"100-150, 150-220",       {{100,150},{150,220}}},
		{"120-180, 180-300",       {{120,180},{180,300}}},
		{"100-150,150-220,220-300",{{100,150},{150,220}}},
		{"100-300",                {{100,300}}},
		{"120-300",                {{120,300}}},
		{"90-120, 120-300",        {{90,120},{120,300}}},
	};
	// the 7th is the same as the 5th for 220-300 (already one bin); keep it out of the table
	schemes.erase(schemes.begin() + 6);

	printf("[BINTUNE] %s %d%s\n", chName, RunYear, whichVFP.c_str());
	printf("[BINTUNE] %-24s | %8s %8s %8s %8s %8s %8s | %s\n", "merged bins", "35-100", "100-120", "120-150", "150-220", "220-300", "300-600", "merged fake rates");
	for (auto &s : schemes) {
		std::vector<double> frNew(fr), errNew(err);
		std::string frs;
		for (auto &m : s.merge) {
			double sy = 0, sfy = 0, sw = 0;
			for (int k = 0; k < nfr; k++) if (lo[k] >= m.first && hi[k] <= m.second) {
				sy += yield[k]; sfy += yield[k]*fr[k]; sw += 1.0/(err[k]*err[k]);
			}
			double frm = sfy/sy, errm = 1.0/std::sqrt(sw);
			for (int k = 0; k < nfr; k++) if (lo[k] >= m.first && hi[k] <= m.second) { frNew[k] = frm; errNew[k] = errm; }
			frs += Form("%g-%g: %.4f+-%.4f  ", m.first, m.second, frm, errm);
		}
		TGraphAsymmErrors *g = new TGraphAsymmErrors();
		int np = 0;
		for (int i = 0; i < g0->GetN(); i++) {
			double x = g0->GetX()[i], y = g0->GetY()[i];
			if (y <= 0 || x < 35 || x >= 600) continue;
			int k = frBin(x);
			double ynew = sig[i]*frNew[k];
			// keep the stat part, replace the fake-rate part of the relative error
			double rel0 = g0->GetErrorYhigh(i)/y, relfr0 = err[k]/fr[k], relfr = errNew[k]/frNew[k];
			double relstat2 = std::max(0.0, rel0*rel0 - relfr0*relfr0);
			double ey = ynew*std::sqrt(relstat2 + relfr*relfr);
			g->SetPoint(np, x, ynew);
			g->SetPointError(np, g0->GetErrorXlow(i), g0->GetErrorXhigh(i), ey, ey);
			np++;
		}
		TF1 *fnum = fitNum(g, Form("num_%d", (int)(&s - &schemes[0])));
		const double rg[7] = {35, 100, 120, 150, 220, 300, 600.1};
		double dat[6] = {0}, fit[6] = {0};
		for (int i = 0; i < g->GetN(); i++) {
			double x = g->GetX()[i], y = g->GetY()[i], l = x - g->GetErrorXlow(i), h = x + g->GetErrorXhigh(i);
			for (int r = 0; r < 6; r++) if (x >= rg[r] && x < rg[r+1]) { dat[r] += y*(h - l); fit[r] += fnum->Integral(l, h); }
		}
		printf("[BINTUNE] %-24s |", s.name);
		for (int r = 0; r < 6; r++) printf(" %8.3f", fit[r] > 0 ? dat[r]/fit[r] : 0);
		printf(" | %s\n", frs.c_str());
	}
}
