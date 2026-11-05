#!/usr/bin/env python3
# ttgamma+VVgamma (rare) shape per era: fraction of the era yield in each signal bin, from the raw pred_rareBkg.C
# histograms (<ch>_rare_norm_raw, kept by factorise_rareBkg.py), compared with the Run 2 combination used as f(bin) in
# the factorised prediction. Lower panel: era / Run 2. Prints the chi2 of each era against the other three eras.
#   python3 plot_rare_shape_eras.py [outdir]      (in the CMSSW environment, after runLimit.sh)
import math
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

OUT = "/uscms_data/d3/tmishra/Output"
outdir = sys.argv[1] if len(sys.argv) > 1 else "."
ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]
STYLE = {"2016preVFP": (ROOT.kBlue + 1, 20), "2016postVFP": (ROOT.kAzure + 7, 21), "2017": (ROOT.kGreen + 2, 22),
         "2018": (ROOT.kRed + 1, 23)}
LABEL = {"mg": "#mu#gamma", "eg": "e#gamma"}

for chfile, ch in (("mg", "mg"), ("egamma", "eg")):
    raw = {}
    for e in ERAS:
        f = ROOT.TFile.Open(f"{OUT}/signalTree_{chfile}_rareBkg_{e}.root")
        h = f.Get(f"{ch}_rare_norm_raw")
        h.SetDirectory(0)
        raw[e] = h
    nb = raw[ERAS[0]].GetNbinsX()

    def fractions(name, eras):
        h = ROOT.TH1D(name, "", nb, 0.5, nb + 0.5)
        N = sum(raw[e].Integral(1, nb) for e in eras)
        for b in range(1, nb + 1):
            h.SetBinContent(b, sum(raw[e].GetBinContent(b) for e in eras) / N)
            h.SetBinError(b, math.sqrt(sum(raw[e].GetBinError(b) ** 2 for e in eras)) / N)
        return h

    run2 = fractions(f"run2_{ch}", ERAS)
    c = ROOT.TCanvas(f"c_{ch}", "", 800, 700)
    top = ROOT.TPad("top", "", 0, 0.32, 1, 1)
    bot = ROOT.TPad("bot", "", 0, 0, 1, 0.32)
    top.SetBottomMargin(0.02); top.SetLogy(); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.32)
    for p in (top, bot):
        p.SetLeftMargin(0.12); p.SetRightMargin(0.04); p.Draw()
    top.cd()
    run2.SetLineColor(ROOT.kBlack); run2.SetLineWidth(2); run2.SetFillColor(ROOT.kGray); run2.SetFillStyle(1001)
    run2.SetMinimum(2e-4); run2.SetMaximum(40.)
    run2.GetYaxis().SetTitle("Fraction of era yield"); run2.GetXaxis().SetLabelSize(0)
    run2.GetYaxis().SetTitleSize(0.055); run2.GetYaxis().SetLabelSize(0.045); run2.GetYaxis().SetTitleOffset(0.95)
    run2.Draw("hist")
    band = run2.Clone(f"band_{ch}"); band.SetFillColor(ROOT.kGray + 2); band.SetFillStyle(3354); band.SetMarkerSize(0)
    band.Draw("E2 same")
    leg = ROOT.TLegend(0.62, 0.66, 0.94, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetNColumns(2)
    leg.AddEntry(run2, "Run 2 (f)", "f")
    ratios, keep = [], []
    print(f"== {ch}")
    for i, e in enumerate(ERAS):
        h = fractions(f"{ch}_{e}", [e])
        col, mk = STYLE[e]
        h.SetLineColor(col); h.SetMarkerColor(col); h.SetMarkerStyle(mk)
        # small horizontal offset per era, so that the points do not overlap
        g = ROOT.TGraphErrors(nb)
        for b in range(1, nb + 1):
            g.SetPoint(b - 1, b - 0.24 + 0.16 * i, h.GetBinContent(b)); g.SetPointError(b - 1, 0, h.GetBinError(b))
        g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(mk); g.SetMarkerSize(0.9)
        g.Draw("P same"); leg.AddEntry(g, e, "p"); keep += [h, g]
        r = ROOT.TGraphErrors(nb)
        for b in range(1, nb + 1):
            v = run2.GetBinContent(b)
            r.SetPoint(b - 1, b - 0.24 + 0.16 * i, h.GetBinContent(b) / v if v > 0 else 0)
            r.SetPointError(b - 1, 0, h.GetBinError(b) / v if v > 0 else 0)
        r.SetLineColor(col); r.SetMarkerColor(col); r.SetMarkerStyle(mk); r.SetMarkerSize(0.8); ratios.append(r)
        # chi2 against the other three eras (independent samples)
        rest = fractions(f"rest_{ch}_{e}", [x for x in ERAS if x != e])
        chi2 = sum((h.GetBinContent(b) - rest.GetBinContent(b)) ** 2 / (h.GetBinError(b) ** 2 + rest.GetBinError(b) ** 2)
                   for b in range(1, nb + 1) if h.GetBinError(b) ** 2 + rest.GetBinError(b) ** 2 > 0)
        print(f"  {e:12s} N = {raw[e].Integral(1, nb):7.1f}  chi2/ndf = {chi2:5.1f}/{nb - 1}  p = {ROOT.TMath.Prob(chi2, nb - 1):.2f}")
    leg.Draw()
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.055)
    lat.DrawLatex(0.15, 0.84, f"{LABEL[ch]}, t#bar{{t}}#gamma + VV#gamma (simulation)")
    lat.SetTextSize(0.042)
    lat.DrawLatex(0.25, 0.06, "p_{T}^{#gamma} < 200 GeV"); lat.DrawLatex(0.66, 0.06, "p_{T}^{#gamma} > 200 GeV")
    for x in (9.5,):
        ln = ROOT.TLine(x, 2e-4, x, 40.); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
    bot.cd()
    frame = ROOT.TH1D(f"frame_{ch}", "", nb, 0.5, nb + 0.5)
    frame.SetMinimum(0.); frame.SetMaximum(2.5)
    frame.GetXaxis().SetTitle("Search bin"); frame.GetYaxis().SetTitle("Era / Run 2")
    for ax, s in ((frame.GetXaxis(), 0.12), (frame.GetYaxis(), 0.10)):
        ax.SetTitleSize(s); ax.SetLabelSize(0.10)
    frame.GetYaxis().SetTitleOffset(0.5); frame.GetYaxis().SetNdivisions(505); frame.GetXaxis().SetTitleOffset(1.1)
    frame.Draw()
    rb = ROOT.TH1D(f"rb_{ch}", "", nb, 0.5, nb + 0.5)
    for b in range(1, nb + 1):
        v = run2.GetBinContent(b)
        rb.SetBinContent(b, 1.); rb.SetBinError(b, run2.GetBinError(b) / v if v > 0 else 0)
    rb.SetFillColor(ROOT.kGray + 2); rb.SetFillStyle(3354); rb.SetMarkerSize(0); rb.Draw("E2 same")
    one = ROOT.TLine(0.5, 1, nb + 0.5, 1); one.SetLineStyle(2); one.Draw()
    for r in ratios:
        r.Draw("P same")
    for ext in ("pdf", "png"):
        c.SaveAs(f"{outdir}/rare_shape_eras_{ch}.{ext}")
