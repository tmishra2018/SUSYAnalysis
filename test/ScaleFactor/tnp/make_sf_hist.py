# Convert fit_all.py's sf_summary.json into a TH2F "EGamma_SF2D" (x = signed SC eta, y = pt),
# the layout of sf/R9_SF.root read by include/analysis_scalefactor.h.
# Run with PyROOT (cmsenv in CMSSW_14_0_7): python3 make_sf_hist.py sf_summary.json SF_R9_2017.root [--cnc]
import argparse
import json
import math
from array import array

import ROOT

ap = argparse.ArgumentParser()
ap.add_argument("summary")
ap.add_argument("output")
ap.add_argument("--cnc", action="store_true", help="use cut-and-count SFs instead of fitted ones")
ap.add_argument("--hybrid", action="store_true", help="fitted SFs, falling back to cut-and-count where the fit is unusable")
ap.add_argument("--max-err", type=float, default=0.05, help="--hybrid: fitted SF error above which cut-and-count is used")
args = ap.parse_args()


def fit_usable(r):
    sf, err = r.get("sf", float("nan")), r.get("sf_err", float("nan"))
    return math.isfinite(sf) and math.isfinite(err) and sf > 0 and err <= args.max_err

rows = json.load(open(args.summary))
eta_edges = sorted({e for r in rows for e in r["eta"]})
pt_edges = sorted({p for r in rows for p in r["pt"]})

out = ROOT.TFile(args.output, "RECREATE")
names = {"EGamma_SF2D": "sf", "EGamma_EffData2D": "eff_data", "EGamma_EffMC2D": "eff_mc"}
hists = {}
for hname in names:
    h = ROOT.TH2F(hname, ";SC #eta;p_{T} [GeV]", len(eta_edges) - 1, array("d", eta_edges), len(pt_edges) - 1, array("d", pt_edges))
    h.SetDirectory(out)
    hists[hname] = h

for r in rows:
    ix = hists["EGamma_SF2D"].GetXaxis().FindBin(0.5 * (r["eta"][0] + r["eta"][1]))
    iy = hists["EGamma_SF2D"].GetYaxis().FindBin(0.5 * (r["pt"][0] + r["pt"][1]))
    use_cnc = args.cnc or (args.hybrid and not fit_usable(r))
    src = r["cut_and_count"] if use_cnc else r.get("dcb_cms", {})
    sf = src.get("sf", float("nan"))
    err = src.get("sf_stat", float("nan")) if use_cnc else r.get("sf_err", float("nan"))
    if args.hybrid and use_cnc and math.isfinite(sf):
        print(f"cut-and-count used for eta {r['eta']} pt {r['pt']}: SF {sf:.4f} +- {err:.4f}")
    if not math.isfinite(sf) or not math.isfinite(err):
        # empty bins (EB-EE gap, excluded by the analysis electron selection): neutral SF
        sf, err = 1.0, 0.0
    hists["EGamma_SF2D"].SetBinContent(ix, iy, sf)
    hists["EGamma_SF2D"].SetBinError(ix, iy, err)
    for hname in ("EGamma_EffData2D", "EGamma_EffMC2D"):
        v = src.get(names[hname], float("nan"))
        hists[hname].SetBinContent(ix, iy, v if math.isfinite(v) else 0.0)

out.Write()
out.Close()
print("wrote", args.output)
