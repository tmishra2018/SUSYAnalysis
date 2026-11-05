#!/usr/bin/env python3
# Factorised ttgamma+VVgamma (rare) prediction in the signal bins. Run after pred_rareBkg.C for all four eras and
# before plot_eventct.C (runLimit.sh does this):
#   y(bin, era) = N(era) * f(bin),   N(era) = era yield summed over the signal bins (per channel),
#                                    f(bin) = sum_eras y(bin) / sum_eras N  (bin fraction from all Run 2 MC).
# The NLO samples (TTGJets, TTJets, WWG, WZG) carry signed generator weights; with a few events per bin and era a
# bin can come out negative (2016postVFP eg bin 9: -0.05 +- 0.44). The shape of these processes does not depend on
# the era, so f is taken from the 4x larger Run 2 sample, and N, which has hundreds of events, stays per era.
# Every variation (JES, JER, ESF, PU, muR, muF, PDF) is factorised the same way with its own N and f, and the
# syserr histograms are recomputed from them as in pred_rareBkg.C.
# Statistical uncertainty: bin error of <ch>_rare_norm = shape part (f, shared by all eras, ttGam_VVGam_stat<bin>);
# <ch>_rare_statN = relative stat. uncertainty of N per era (ttGam_VVGam_statN_<ch>_<era> in the cards).
# Raw histograms are kept as <name>_raw, so rerunning the script starts again from the pred_rareBkg.C output.
#   python3 factorise_rareBkg.py [dir]      (in the CMSSW environment; default dir /uscms_data/d3/tmishra/Output)
import sys
import ROOT

OUT = sys.argv[1] if len(sys.argv) > 1 else "/uscms_data/d3/tmishra/Output"
ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]
# total luminosity uncertainty per period from the 3-component scheme of the datacards (convert_dict.yaml:
# lumi_13TeV_1516_l, _151617_l, _15161718_l): 2016 1.23%, 2017 0.82%, 2018 0.84% (Run 2 combination: merge678.py)
LUMI_UNC = {"2016preVFP": (0.0118**2 + 0.0004**2 + 0.0035**2) ** 0.5, "2016postVFP": (0.0118**2 + 0.0004**2 + 0.0035**2) ** 0.5,
            "2017": (0.0055**2 + 0.0061**2) ** 0.5, "2018": 0.0084}
CHANNELS = {"egamma": "eg", "mg": "mg"}
VARS = ["norm", "jesUp", "jesDown", "jerUp", "jerDown", "esfUp", "esfphoUp", "puUp", "puDown",
        "muRup", "muRdown", "muFup", "muFdown", "pdfUp", "pdfDown"]

for chfile, ch in CHANNELS.items():
    files = {e: ROOT.TFile.Open(f"{OUT}/signalTree_{chfile}_rareBkg_{e}.root", "UPDATE") for e in ERAS}
    raw = {}  # raw[var][era]: histogram as written by pred_rareBkg.C
    for v in VARS:
        raw[v] = {}
        for e, f in files.items():
            h = f.Get(f"{ch}_rare_{v}_raw") or f.Get(f"{ch}_rare_{v}")
            h = h.Clone(f"{ch}_rare_{v}_raw")
            h.SetDirectory(0)
            raw[v][e] = h
    nbin = raw["norm"][ERAS[0]].GetNbinsX()

    new = {}  # new[var][era]: factorised histogram
    for v in VARS:
        N = {e: raw[v][e].Integral(1, nbin) for e in ERAS}  # signal bins only (outside events go to the underflow)
        Ntot = sum(N.values())
        new[v] = {}
        for e in ERAS:
            h = raw[v][e].Clone(f"{ch}_rare_{v}")
            h.SetDirectory(0)
            for b in range(1, nbin + 1):
                yb = sum(raw[v][x].GetBinContent(b) for x in ERAS)
                if v == "norm" and yb <= 0:
                    sys.exit(f"{ch} bin {b}: Run 2 rare yield {yb:.4f} <= 0, factorisation not possible")
                h.SetBinContent(b, N[e] * yb / Ntot)
                if v == "norm":  # shape (f) part of the stat. uncertainty, shared by all eras
                    sb = sum(raw[v][x].GetBinError(b) ** 2 for x in ERAS) ** 0.5
                    h.SetBinError(b, h.GetBinContent(b) * sb / yb)
            new[v][e] = h

    for e, f in files.items():
        f.cd()
        for v in VARS:
            raw[v][e].Write(f"{ch}_rare_{v}_raw", ROOT.TObject.kOverwrite)
            new[v][e].Write(f"{ch}_rare_{v}", ROOT.TObject.kOverwrite)
        n = {v: new[v][e] for v in VARS}
        nom = n["norm"]
        syserr = {k: f.Get(f"{ch}_rare_syserr_{k}") for k in ("jes", "jer", "esf", "esfpho", "pu", "xs", "lumi",
                                                               "qcdmuR", "qcdmuF", "pdf")}
        for b in range(1, nbin + 1):
            y = nom.GetBinContent(b)
            d = lambda v: abs(n[v].GetBinContent(b) - y)
            r = lambda v: abs(n[v].GetBinContent(b) / y - 1.)
            syserr["jes"].SetBinContent(b, max(d("jesUp"), d("jesDown")))
            syserr["jer"].SetBinContent(b, max(d("jerUp"), d("jerDown")))
            syserr["esf"].SetBinContent(b, d("esfUp"))
            syserr["esfpho"].SetBinContent(b, d("esfphoUp"))
            syserr["pu"].SetBinContent(b, max(d("puUp"), d("puDown")))
            syserr["xs"].SetBinContent(b, 0.5 * y)
            syserr["lumi"].SetBinContent(b, LUMI_UNC[e] * y)
            syserr["qcdmuR"].SetBinContent(b, max(r("muRup"), r("muRdown")))
            syserr["qcdmuF"].SetBinContent(b, max(r("muFup"), r("muFdown")))
            syserr["pdf"].SetBinContent(b, max(r("pdfUp"), r("pdfDown")))
        for h in syserr.values():
            h.Write("", ROOT.TObject.kOverwrite)
        # relative stat. uncertainty of N(era), in one bin
        hN = ROOT.TH1D(f"{ch}_rare_statN", "relative stat. uncertainty of the era yield", 1, 0, 1)
        r0 = raw["norm"][e]
        hN.SetBinContent(1, sum(r0.GetBinError(b) ** 2 for b in range(1, nbin + 1)) ** 0.5 / r0.Integral(1, nbin))
        hN.Write("", ROOT.TObject.kOverwrite)
        print(f"{ch} {e:12s} N = {r0.Integral(1, nbin):8.2f} +- {100 * hN.GetBinContent(1):.1f}%   bins: "
              + " ".join(f"{r0.GetBinContent(b):.2f}->{nom.GetBinContent(b):.2f}" for b in range(1, nbin + 1)))
        f.Close()
