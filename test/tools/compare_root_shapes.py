#!/usr/bin/env python3
"""
compare_root_shapes.py
======================
Compares shapes of sigMET, HT, phoEt between:
  - *_original.root  vs  *.root
for two TTrees:  fakeLepTree  and  signalTree
across four years: 2016preVFP, 2016postVFP, 2017, 2018

Produces:
  compare_fakeLepTree.png   (12 panels)
  compare_signalTree.png    (12 panels)

Run with:
  python3 compare_root_shapes.py
  (needs PyROOT, numpy, scipy, matplotlib  ->  works in CMSSW env)
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
except ImportError:
    print("ERROR: PyROOT not available. Source CMSSW or use uproot.")
    sys.exit(1)

# ------------------------------------------------------------------ config ---
XRD_BASE = "root://cmseos.fnal.gov//eos/uscms/store/user/tmishra/eg_mg_treesData/"
YEARS    = ["2016preVFP", "2016postVFP", "2017", "2018"]
VARS     = ["sigMET", "HT", "phoEt"]
VAR_LABEL = {
    "sigMET" : "MET [GeV]",
    "HT"     : "H_{T} [GeV]",
    "phoEt"  : "Photon E_{T} [GeV]",
}
TREES = ["fakeLepTree", "signalTree"]
KS_THRESHOLD = 0.05   # p-value below this → shapes differ


# --------------------------------------------------------------- helpers ----
def file_url(year, original=False):
    tag = f"_original" if original else ""
    return XRD_BASE + f"resTree_mgsignal_MuonEG_{year}{tag}.root"


def branch_to_numpy(tfile_path, tree_name, branch_name):
    """Open a ROOT file and return a numpy array for one branch."""
    f = ROOT.TFile.Open(tfile_path)
    if not f or f.IsZombie():
        print(f"    ERROR: cannot open {tfile_path}")
        return None
    tree = f.Get(tree_name)
    if not tree:
        print(f"    ERROR: tree '{tree_name}' not found in {tfile_path}")
        f.Close()
        return None

    arr = []
    for evt in tree:
        val = getattr(evt, branch_name, None)
        if val is not None:
            arr.append(float(val))
    f.Close()
    if not arr:
        print(f"    ERROR: branch '{branch_name}' empty or not found")
        return None
    return np.array(arr, dtype=np.float32)


def plot_panel(ax, orig, new, var, year, tree_name, ks_stat, p_val):
    """Draw normalized step histograms on ax."""
    combined = np.concatenate([orig, new])
    lo, hi   = np.percentile(combined, [0.5, 99.5])
    bins     = np.linspace(lo, hi, 60)

    h_orig, edges = np.histogram(orig, bins=bins)
    h_new,  _     = np.histogram(new,  bins=bins)

    def safe_norm(h):
        s = h.sum()
        return h / s if s > 0 else h

    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.step(centers, safe_norm(h_orig), where='mid',
            color='steelblue', lw=1.8, label='original')
    ax.step(centers, safe_norm(h_new),  where='mid',
            color='tomato',    lw=1.8, label='processed', linestyle='--')

    ax.set_xlabel(VAR_LABEL.get(var, var), fontsize=9)
    ax.set_ylabel("Normalized entries", fontsize=9)
    ax.set_title(f"{year}  |  {var}", fontsize=10, fontweight='bold')
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)

    differs = (p_val < KS_THRESHOLD)
    color   = 'red' if differs else 'green'
    if differs:
        ax.set_facecolor('#fff0f0')

    ax.text(0.97, 0.95,
            f"KS = {ks_stat:.3f}\np  = {p_val:.2e}",
            ha='right', va='top', transform=ax.transAxes,
            fontsize=7.5, color=color,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))


# --------------------------------------------------------------- main loop --
for tree_name in TREES:
    print()
    print("=" * 65)
    print(f"  Tree: {tree_name}")
    print("=" * 65)

    fig, axes = plt.subplots(4, 3, figsize=(15, 18))
    fig.suptitle(
        f"Shape Comparison — {tree_name}\n"
        f"blue = *_original.root   |   red dashed = *.root\n"
        f"Red background / KS annotation = shapes differ (p < {KS_THRESHOLD})",
        fontsize=12, fontweight='bold', y=0.99
    )

    for row, year in enumerate(YEARS):
        orig_path = file_url(year, original=True)
        new_path  = file_url(year, original=False)

        for col, var in enumerate(VARS):
            ax = axes[row, col]
            tag = f"[{year}] {var}"
            print(f"\n  {tag}  ({tree_name})")

            orig = branch_to_numpy(orig_path, tree_name, var)
            new  = branch_to_numpy(new_path,  tree_name, var)

            if orig is None or new is None:
                ax.text(0.5, 0.5, "Load Error",
                        ha='center', va='center', transform=ax.transAxes,
                        fontsize=12, color='red')
                ax.set_title(f"{year} | {var}", fontsize=10, fontweight='bold')
                continue

            print(f"    original  : n={len(orig):>8,}  "
                  f"mean={orig.mean():>10.3f}  std={orig.std():>9.3f}")
            print(f"    processed : n={len(new):>8,}  "
                  f"mean={new.mean():>10.3f}  std={new.std():>9.3f}")

            ks_stat, p_val = ks_2samp(orig, new)

            plot_panel(ax, orig, new, var, year, tree_name, ks_stat, p_val)

            if p_val < KS_THRESHOLD:
                mean_diff = abs(orig.mean() - new.mean())
                std_diff  = abs(orig.std()  - new.std())
                print(f"    *** SHAPE DIFFERS ***  "
                      f"KS={ks_stat:.4f}  p={p_val:.4e}  "
                      f"Δmean={mean_diff:.3f}  Δstd={std_diff:.3f}")
            else:
                print(f"    Shape consistent.      "
                      f"KS={ks_stat:.4f}  p={p_val:.4e}")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outname = f"compare_{tree_name}.png"
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  >> Saved: {outname}")

print("\nAll done.")
