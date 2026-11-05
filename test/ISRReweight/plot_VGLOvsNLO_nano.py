#!/usr/bin/env python3
# Ratios and plots from run_VGLOvsNLO_nano.sh: NLO/LO generator-level ISR pT (unit area) for W+gamma and DY,
# and signed/unsigned NLO shapes for Z+gamma and DY (the V+gamma skim ignores the genWeight sign).
#   python3 plot_VGLOvsNLO_nano.py
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
D = "/eos/uscms/store/user/tmishra/ISRweighting/VGLOvsNLO_nano"


def get(s, name):
    f = ROOT.TFile.Open(f"{D}/VGLOvsNLO_nano_{s}.root")
    h = f.Get(name).Clone(f"{s}_{name}")
    h.SetDirectory(0)
    f.Close()
    if h.Integral() > 0:
        h.Scale(1. / h.Integral())
    return h


def ratio_line(num, den):
    r = num.Clone(num.GetName() + "_r")
    r.Divide(den)
    return [(r.GetBinContent(i), r.GetBinError(i)) for i in range(1, r.GetNbinsX() + 1)]


labels = ["0-50", "50-100", "100-150", "150-200", "200-250", "250-300", ">300"]
comps = [
    ("Wgamma NLO/LO", ("WG_NLO", "isr_signed"), ("WG_LO", "isr_unsigned")),
    ("DY NLO/LO", ("DY_NLO", "isr_signed"), ("DY_LO", "isr_unsigned")),
    ("DY NLO unsigned/LO", ("DY_NLO", "isr_unsigned"), ("DY_LO", "isr_unsigned")),
    ("Zgamma NLO unsigned/signed", ("ZG_NLO", "isr_unsigned"), ("ZG_NLO", "isr_signed")),
    ("DY NLO unsigned/signed", ("DY_NLO", "isr_unsigned"), ("DY_NLO", "isr_signed")),
]
lines = ["ISR pT (GeV)            " + " ".join(f"{l:>13s}" for l in labels)]
for name, a, b in comps:
    r = ratio_line(get(*a), get(*b))
    lines.append(f"{name:24s}" + " ".join(f"{v:6.3f}+-{e:5.3f}" for v, e in r))
for s in ["WG_LO", "WG_NLO", "DY_LO", "DY_NLO", "ZG_NLO"]:
    f = ROOT.TFile.Open(f"{D}/VGLOvsNLO_nano_{s}.root")
    c = f.Get("h_count")
    hs = f.Get("isr_signed")
    hu = f.Get("isr_unsigned")
    lines.append(f"{s}: events {c.GetBinContent(1):.0f}, selected {hu.GetEntries():.0f} "
                 f"(signed sum {hs.Integral(0, hs.GetNbinsX() + 1):.0f})")
    f.Close()
txt = "\n".join(lines)
print(txt)
open(f"{D}/ratio_table.txt", "w").write(txt + "\n")

# plots: per variable, three panels (W+gamma LO vs NLO, DY LO vs NLO, Z+gamma signed vs unsigned)
for var, xt in (("isrfine", "gen ISR p_{T} (GeV)"), ("ht", "gen H_{T} (GeV)"), ("phopt", "gen photon p_{T} (GeV)")):
    c = ROOT.TCanvas("c", "", 1500, 600)
    c.Divide(3, 1)
    keep = []
    for ip, (title, h1, h2, l1, l2) in enumerate([
            ("W#gamma", ("WG_LO", f"{var}_unsigned"), ("WG_NLO", f"{var}_signed"), "LO (MLM)", "NLO (FxFx, signed)"),
            ("DY", ("DY_LO", f"{var}_unsigned"), ("DY_NLO", f"{var}_signed"), "LO (MLM)", "NLO (FxFx, signed)"),
            ("Z#gamma (NLO)", ("ZG_NLO", f"{var}_unsigned"), ("ZG_NLO", f"{var}_signed"), "unsigned", "signed")]):
        pad = c.cd(ip + 1)
        up = ROOT.TPad(f"up{ip}", "", 0, 0.3, 1, 1); up.SetBottomMargin(0.02); up.SetLogy(); up.Draw()
        lo = ROOT.TPad(f"lo{ip}", "", 0, 0, 1, 0.3); lo.SetTopMargin(0.02); lo.SetBottomMargin(0.3); lo.Draw()
        a, b = get(*h1), get(*h2)
        a.SetLineColor(ROOT.kBlue); a.SetMarkerColor(ROOT.kBlue); a.SetMarkerStyle(24)
        b.SetLineColor(ROOT.kBlack); b.SetMarkerColor(ROOT.kBlack); b.SetMarkerStyle(21)
        up.cd(); a.SetTitle(f"{title}, UL18 gen level"); a.GetXaxis().SetLabelSize(0); a.Draw("E"); b.Draw("E same")
        leg = ROOT.TLegend(0.45, 0.72, 0.88, 0.88); leg.AddEntry(a, l1, "lep"); leg.AddEntry(b, l2, "lep"); leg.Draw()
        lo.cd(); r = b.Clone(f"r{ip}{var}"); r.Divide(a); r.SetTitle(""); r.GetYaxis().SetRangeUser(0.5, 1.5)
        r.GetYaxis().SetTitle(f"{l2.split()[0]}/{l1.split()[0]}"); r.GetYaxis().SetTitleSize(0.1); r.GetYaxis().SetTitleOffset(0.5)
        r.GetYaxis().SetLabelSize(0.08); r.GetXaxis().SetLabelSize(0.1); r.GetXaxis().SetTitleSize(0.11); r.GetXaxis().SetTitle(xt)
        r.Draw("E"); one = ROOT.TLine(r.GetXaxis().GetXmin(), 1, r.GetXaxis().GetXmax(), 1); one.SetLineStyle(2); one.Draw()
        keep += [up, lo, a, b, r, leg, one]
    c.SaveAs(f"{D}/VGLOvsNLO_nano_{var}.pdf")
    c.SaveAs(f"{D}/VGLOvsNLO_nano_{var}.png")
