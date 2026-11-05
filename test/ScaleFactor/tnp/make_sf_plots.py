#!/usr/bin/env python3
# 2D SF / efficiency maps and 1D SF-vs-pT plots for the egamma-tnp electron R9 and miniIso, and photon R9 (phoR9) results.
# Usage (cmsenv):  python3 make_sf_plots.py [outdir]     default outdir: plots/
import ROOT, sys, os
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetPaintTextFormat(".3f"); ROOT.gStyle.SetPalette(ROOT.kBird)
here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "plots")
os.makedirs(out, exist_ok=True)
years = ["2016preVFP", "2016postVFP", "2017", "2018"]
for sf in ("R9", "miniIso", "phoR9"):
    label = "photon R9" if sf == "phoR9" else sf
    ptmin = 35 if sf == "phoR9" else 25
    for y in years:
        fn = f"{here}/results/{y}/SF_{sf}_{y}_hybrid.root"
        f = ROOT.TFile.Open(fn)
        if not f: print("missing", fn); continue
        for hname, title, zr in (("EGamma_SF2D", "SF (data/MC)", (0.94, 1.04)),
                                  ("EGamma_EffData2D", "efficiency, data", (0.6, 1.0)),
                                  ("EGamma_EffMC2D", "efficiency, MC", (0.6, 1.0))):
            h = f.Get(hname)
            c = ROOT.TCanvas("c", "", 1100, 750); c.SetLogy(); c.SetRightMargin(0.14)
            h.SetTitle(f"{label} {title} {y};#eta_{{SC}};p_{{T}} [GeV]")
            h.SetMinimum(zr[0]); h.SetMaximum(zr[1]); h.SetMarkerSize(0.9)
            h.GetYaxis().SetMoreLogLabels(); h.Draw("COLZ TEXT45E" if hname == "EGamma_SF2D" else "COLZ TEXT45")
            c.SaveAs(f"{out}/{sf}_{y}_{hname.replace('EGamma_','')}.png"); c.SaveAs(f"{out}/{sf}_{y}_{hname.replace('EGamma_','')}.pdf")
        # 1D: SF vs pT, one curve per eta bin
        h = f.Get("EGamma_SF2D")
        c = ROOT.TCanvas("c1", "", 1100, 650); c.SetLogx(); c.SetGridy(); c.SetRightMargin(0.25)
        fr = c.DrawFrame(ptmin, 0.90 if sf == "R9" else 0.94, 500, 1.06); fr.SetTitle(f"{label} SF {y};p_{{T}} [GeV];SF")
        leg = ROOT.TLegend(0.76, 0.30, 0.99, 0.88); leg.SetBorderSize(0); leg.SetHeader("#eta_{SC} bin"); keep = []
        cols = [1, 2, 4, 6, 8, 9, 28, 46, 30, 38]
        for ix in range(1, h.GetNbinsX() + 1):
            lo, hi = h.GetXaxis().GetBinLowEdge(ix), h.GetXaxis().GetBinUpEdge(ix)
            if 1.4442 < abs(0.5 * (lo + hi)) < 1.566: continue  # EB-EE gap
            g = ROOT.TGraphErrors()
            for iy in range(1, h.GetNbinsY() + 1):
                n = g.GetN(); g.SetPoint(n, h.GetYaxis().GetBinCenter(iy) * (1 + 0.01 * (ix - 5)), h.GetBinContent(ix, iy)); g.SetPointError(n, 0, h.GetBinError(ix, iy))
            col = cols[(ix - 1) % len(cols)]; g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(20 + ix % 8); g.Draw("PE same")
            leg.AddEntry(g, f"[{lo:.2f},{hi:.2f}]", "lp"); keep.append(g)
        leg.Draw(); c.SaveAs(f"{out}/{sf}_{y}_SF_vs_pt.png"); c.SaveAs(f"{out}/{sf}_{y}_SF_vs_pt.pdf")
print("plots written to", out)
