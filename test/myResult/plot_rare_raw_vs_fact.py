#!/usr/bin/env python3
# ttgamma+VVgamma (rare) prediction per era in the signal search bins: raw pred_rareBkg.C output (<ch>_rare_norm_raw,
# points with MC stat. error) overlaid with the factorised prediction of factorise_rareBkg.py (<ch>_rare_norm,
# histogram; band = shape stat. (f) + era yield stat. (statN) in quadrature). Lower panel: raw / factorised.
# One canvas per channel with the four eras (rare_raw_vs_fact_<ch>.pdf/png) and one per era (..._<ch>_<era>.pdf/png).
#   python3 plot_rare_raw_vs_fact.py [outdir]      (in the CMSSW environment, after runLimit.sh;
#   default outdir /eos/uscms/store/user/tmishra/SUSYAnalysis_out/test/myResult/rare_raw_vs_fact)
import math
import os
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(0)

OUT = "/uscms_data/d3/tmishra/Output"
outdir = sys.argv[1] if len(sys.argv) > 1 else "/eos/uscms/store/user/tmishra/SUSYAnalysis_out/test/myResult/rare_raw_vs_fact"
os.makedirs(outdir, exist_ok=True)
ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]
LABEL = {"mg": "#mu#gamma", "eg": "e#gamma"}
YMIN, YMAX = 1e-2, 5e3
RMIN, RMAX = -0.5, 3.0
keep = []  # ROOT objects that must outlive the drawing functions


def draw_era(pad_top, pad_bot, ch, era, h_raw, h_fact, statN, text_scale=1.0, tag=""):
    nb = h_raw.GetNbinsX()
    pad_top.cd()
    # factorised: histogram + uncertainty band (shape stat. + era yield stat.)
    fact = ROOT.TH1D(f"fact_{ch}_{era}{tag}", "", nb, 0.5, nb + 0.5)
    band = ROOT.TH1D(f"band_{ch}_{era}{tag}", "", nb, 0.5, nb + 0.5)
    for b in range(1, nb + 1):
        y = h_fact.GetBinContent(b)
        fact.SetBinContent(b, y)
        band.SetBinContent(b, y)
        band.SetBinError(b, math.hypot(h_fact.GetBinError(b), statN * y))
    fact.SetLineColor(ROOT.kRed + 1); fact.SetLineWidth(2)
    band.SetFillColor(ROOT.kRed - 9); band.SetFillStyle(1001); band.SetMarkerSize(0); band.SetLineWidth(0)
    fact.SetMinimum(YMIN); fact.SetMaximum(YMAX)
    fact.GetYaxis().SetTitle("Events / bin"); fact.GetXaxis().SetLabelSize(0)
    fact.GetYaxis().SetTitleSize(0.06 * text_scale); fact.GetYaxis().SetLabelSize(0.05 * text_scale)
    fact.GetYaxis().SetTitleOffset(0.9 / text_scale)
    fact.Draw("hist")
    band.Draw("E2 same")
    fact.Draw("hist same")
    # raw: points with MC stat. error; a bin with raw <= 0 is drawn as an open marker on the x axis
    g = ROOT.TGraphAsymmErrors()
    gneg = ROOT.TGraph()
    for b in range(1, nb + 1):
        y, e = h_raw.GetBinContent(b), h_raw.GetBinError(b)
        if y > 0:
            i = g.GetN()
            g.SetPoint(i, b, y)
            g.SetPointError(i, 0, 0, min(e, y - 1.001 * YMIN), e)
        else:
            gneg.SetPoint(gneg.GetN(), b, 1.4 * YMIN)
    g.SetMarkerStyle(20); g.SetMarkerSize(0.9); g.SetLineWidth(2)
    g.Draw("P same")
    if gneg.GetN():
        gneg.SetMarkerStyle(32); gneg.SetMarkerSize(1.3); gneg.Draw("P same")
    ln = ROOT.TLine(9.5, YMIN, 9.5, YMAX); ln.SetLineStyle(2); ln.Draw()
    leg = ROOT.TLegend(0.55, 0.62, 0.94, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextSize(0.045 * text_scale)
    leg.AddEntry(g, "raw (pred_rareBkg.C)", "lep")
    leg.AddEntry(band, "factorised", "lf")
    if gneg.GetN():
        leg.AddEntry(gneg, "raw #leq 0", "p")
    leg.Draw()
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.055 * text_scale)
    lat.DrawLatex(0.15, 0.84, f"{LABEL[ch]}, {era}")
    lat.SetTextSize(0.042 * text_scale)
    lat.DrawLatex(0.15, 0.77, "t#bar{t}#gamma + VV#gamma (simulation)")
    lat.SetTextSize(0.04 * text_scale)
    lat.DrawLatex(0.25, 0.06, "p_{T}^{#gamma} < 200 GeV"); lat.DrawLatex(0.66, 0.06, "p_{T}^{#gamma} > 200 GeV")

    pad_bot.cd()
    frame = ROOT.TH1D(f"frame_{ch}_{era}{tag}", "", nb, 0.5, nb + 0.5)
    frame.SetMinimum(RMIN); frame.SetMaximum(RMAX)
    frame.GetXaxis().SetTitle("Search bin"); frame.GetYaxis().SetTitle("Raw / fact.")
    frame.GetXaxis().SetTitleSize(0.12 * text_scale); frame.GetYaxis().SetTitleSize(0.10 * text_scale)
    frame.GetXaxis().SetLabelSize(0.10 * text_scale); frame.GetYaxis().SetLabelSize(0.09 * text_scale)
    frame.GetYaxis().SetTitleOffset(0.5 / text_scale); frame.GetYaxis().SetNdivisions(505)
    frame.GetXaxis().SetTitleOffset(1.1)
    for b in range(1, nb + 1):
        frame.GetXaxis().SetBinLabel(b, str(b))
    frame.GetXaxis().SetLabelSize(0.12 * text_scale); frame.GetXaxis().LabelsOption("h")
    frame.Draw()
    rb = ROOT.TH1D(f"rb_{ch}_{era}{tag}", "", nb, 0.5, nb + 0.5)
    r = ROOT.TGraphErrors()
    for b in range(1, nb + 1):
        y = h_fact.GetBinContent(b)
        rb.SetBinContent(b, 1.)
        rb.SetBinError(b, band.GetBinError(b) / y if y > 0 else 0)
        if y > 0:
            r.SetPoint(r.GetN(), b, h_raw.GetBinContent(b) / y)
            r.SetPointError(r.GetN() - 1, 0, h_raw.GetBinError(b) / y)
    rb.SetFillColor(ROOT.kRed - 9); rb.SetFillStyle(1001); rb.SetMarkerSize(0); rb.Draw("E2 same")
    one = ROOT.TLine(0.5, 1, nb + 0.5, 1); one.SetLineStyle(2); one.Draw()
    zero = ROOT.TLine(0.5, 0, nb + 0.5, 0); zero.SetLineStyle(3); zero.Draw()
    r.SetMarkerStyle(20); r.SetMarkerSize(0.8); r.SetLineWidth(2); r.Draw("P same")
    frame.Draw("axis same")
    keep.extend([fact, band, g, gneg, ln, leg, lat, frame, rb, r, one, zero])


def split(c, x1, y1, x2, y2, name):
    top = ROOT.TPad(f"top_{name}", "", x1, y1 + 0.32 * (y2 - y1), x2, y2)
    bot = ROOT.TPad(f"bot_{name}", "", x1, y1, x2, y1 + 0.32 * (y2 - y1))
    top.SetBottomMargin(0.02); top.SetLogy(); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.32)
    for p in (top, bot):
        p.SetLeftMargin(0.12); p.SetRightMargin(0.04); c.cd(); p.Draw()
    keep.extend([top, bot])
    return top, bot


for chfile, ch in (("egamma", "eg"), ("mg", "mg")):
    data = {}
    for e in ERAS:
        f = ROOT.TFile.Open(f"{OUT}/signalTree_{chfile}_rareBkg_{e}.root")
        hr, hf, hN = f.Get(f"{ch}_rare_norm_raw"), f.Get(f"{ch}_rare_norm"), f.Get(f"{ch}_rare_statN")
        if not (hr and hf and hN):
            sys.exit(f"{f.GetName()}: missing raw/factorised histograms, run factorise_rareBkg.py first")
        for h in (hr, hf):
            h.SetDirectory(0)
        data[e] = (hr, hf, hN.GetBinContent(1))
        f.Close()

    # one canvas per era
    for e in ERAS:
        c = ROOT.TCanvas(f"c_{ch}_{e}", "", 800, 700)
        top, bot = split(c, 0, 0, 1, 1, f"{ch}_{e}")
        draw_era(top, bot, ch, e, *data[e])
        for ext in ("pdf", "png"):
            c.SaveAs(f"{outdir}/rare_raw_vs_fact_{ch}_{e}.{ext}")

    # 2x2 summary
    c = ROOT.TCanvas(f"c_{ch}_all", "", 1600, 1400)
    for i, e in enumerate(ERAS):
        x1, y1 = 0.5 * (i % 2), 0.5 * (1 - i // 2)
        top, bot = split(c, x1, y1, x1 + 0.5, y1 + 0.5, f"{ch}_{e}_all")
        draw_era(top, bot, ch, e, *data[e], tag="_all")
    for ext in ("pdf", "png"):
        c.SaveAs(f"{outdir}/rare_raw_vs_fact_{ch}.{ext}")
