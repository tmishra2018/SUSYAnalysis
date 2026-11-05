#!/usr/bin/env python3
# Compare the current V+gamma skim trees (/eos/.../VGamma/resTree_VGamma_<sample>_<era>.root, read by the background
# macros) with the new condor skim (/eos/.../VGamma/<tag>/resTree_VGamma_<sample>_<era>.root, merge_skim.sh), per
# sample, era and channel tree: entries, ntotalevent, fraction of negative crosssection (genWeight sign),
# weighted yield sum(crosssection/ntotalevent * PUweight) (the lumi factor cancels in new/old), mean ISRJetPt.
#   python3 compare_skim.py [--tag zgisr] [--samples WZ] [--eras 2018]     (in the CMSSW environment)
import argparse
import os
import ROOT
ROOT.gROOT.SetBatch(True)

SAMPLES = ["WGToLNuG", "WGJet40", "WGJet130", "TTJets", "DYJetsToLL", "TTGJets", "WWG", "WZG", "WW", "WZ", "ZZ", "ZGToLLG"]
ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]
BASE = "/eos/uscms/store/user/tmishra/VGamma"


def summary(path, tree):
    if not os.path.exists(path):
        return None
    df = ROOT.RDataFrame(tree, path)
    n = df.Count()
    w = df.Define("w_", "crosssection/ntotalevent*PUweight")
    res = dict(n=n, ntot=df.Mean("ntotalevent"), neg=df.Filter("crosssection < 0").Count(),
               yld=w.Sum("w_"), isr=df.Mean("ISRJetPt"))
    out = {k: v.GetValue() for k, v in res.items()}
    out["neg"] = out["neg"] / out["n"] if out["n"] else 0.
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="zgisr")
ap.add_argument("--samples", default=",".join(SAMPLES))
ap.add_argument("--eras", default=",".join(ERAS))
a = ap.parse_args()

hdr = (f"{'sample_era':24s} {'tree':6s} {'entries old':>11s} {'new':>9s} {'ntotalevent old':>15s} {'new':>12s} "
       f"{'neg new':>7s} {'yield new/old':>13s} {'<ISRpT> old':>11s} {'new':>6s}")
lines = [hdr]
for era in a.eras.split(","):
    for s in a.samples.split(","):
        key = f"{s}_{era}"
        for t in ("egTree", "mgTree"):
            o = summary(f"{BASE}/resTree_VGamma_{key}.root", t)
            n = summary(f"{BASE}/{a.tag}/resTree_VGamma_{key}.root", t)
            if o is None or n is None:
                lines.append(f"{key:24s} {t:6s} {'old missing' if o is None else ''} {'new missing' if n is None else ''}")
                continue
            r = n["yld"] / o["yld"] if o["yld"] else float("nan")
            lines.append(f"{key:24s} {t:6s} {o['n']:>11d} {n['n']:>9d} {o['ntot']:>15.0f} {n['ntot']:>12.0f} "
                         f"{n['neg']:>7.3f} {r:>13.3f} {o['isr']:>11.1f} {n['isr']:>6.1f}")
txt = "\n".join(lines)
print(txt)
out = f"{BASE}/{a.tag}/compare_skim.txt"
if os.path.isdir(os.path.dirname(out)):
    open(out, "w").write(txt + "\n")
    print(f"written: {out}")
