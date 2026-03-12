#!/usr/bin/env python3
"""
compare_TH1_shapes.py
=====================
Compares shapes of p_MET, p_HT, p_PhoEt (TH1D histograms stored directly
in ROOT files) between old and new versions for two samples:
  - validTree_mg_rareBkg_2018
  - validTree_mg_VGBkg_2018

Produces:
  compare_rareBkg_2018.png   (3 panels)
  compare_VGBkg_2018.png     (3 panels)

Run with:
  python3 compare_TH1_shapes.py   (inside CMSSW env or where PyROOT is available)
"""

import sys, os
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
    print("ERROR: PyROOT not available. Source CMSSW environment first.")
    sys.exit(1)

# ------------------------------------------------------------------ config ---
SAMPLES = [
    {
        "label"   : "rareBkg 2018",
        "old_file": "/uscms_data/d3/tmishra/Background/validTree_mg_rareBkg_2018_old.root",
        "new_file": "/uscms_data/d3/tmishra/Background/validTree_mg_rareBkg_2018.root",
        "outname" : "compare_rareBkg_2018.png",
    },
    {
        "label"   : "VGBkg 2018",
        "old_file": "/uscms_data/d3/tmishra/Background/validTree_mg_VGBkg_2018_old.root",
        "new_file": "/uscms_data/d3/tmishra/Background/validTree_mg_VGBkg_2018.root",
        "outname" : "compare_VGBkg_2018.png",
    },
]

HISTOGRAMS = ["p_MET", "p_HT", "p_PhoEt"]
HIST_LABEL = {
    "p_MET"   : "MET [GeV]",
    "p_HT"    : "H_{T} [GeV]",
    "p_PhoEt" : "Photon E_{T} [GeV]",
}

KS_THRESHOLD = 0.05   # p-value below this → shapes differ


# --------------------------------------------------------------- helpers ----
def th1_to_numpy(tfile_path, hist_name):
    """Return (bin_centers, bin_contents, bin_errors) for a TH1D in a file."""
    f = ROOT.TFile.Open(tfile_path)
    if not f or f.IsZombie():
        print(f"  ERROR: cannot open {tfile_path}")
        return None, None, None
    h = f.Get(hist_name)
    if not h:
        print(f"  ERROR: histogram '{hist_name}' not found in {tfile_path}")
        f.Close()
        return None, None, None
    nbins   = h.GetNbinsX()
    centers = np.array([h.GetBinCenter(i)  for i in range(1, nbins+1)])
    contents= np.array([h.GetBinContent(i) for i in range(1, nbins+1)])
    errors  = np.array([h.GetBinError(i)   for i in range(1, nbins+1)])
    f.Close()
    return centers, contents, errors


def norm(contents):
    """Normalize histogram contents to unit area."""
    s = contents.sum()
    return contents / s if s > 0 else contents


def expand_to_values(centers, contents):
    """
    Convert binned histogram back to an approximate array of values
    (integer counts rounded) so scipy KS test can be applied.
    This preserves the shape information without needing the raw data.
    """
    vals = []
    for c, n in zip(centers, contents):
        count = int(round(n))
        if count > 0:
            vals.extend([c] * count)
    return np.array(vals)


def plot_panel(ax, centers, old_norm, new_norm, old_err, new_err,
               hist_name, ks_stat, p_val):
    """Draw normalized step histograms with error bands on ax."""
    width = centers[1] - centers[0] if len(centers) > 1 else 1.0
    edges = np.append(centers - width/2, centers[-1] + width/2)

    # Step histograms
    ax.step(np.append(edges[0], edges), np.append(old_norm[0], np.append(old_norm, 0)),
            where='pre', color='steelblue', lw=2.0, label='old', zorder=3)
    ax.step(np.append(edges[0], edges), np.append(new_norm[0], np.append(new_norm, 0)),
            where='pre', color='tomato',    lw=2.0, label='new',
            linestyle='--', zorder=3)

    # Error bands
    old_norm_err = old_err / (old_err.sum() + 1e-30)  # rough propagation
    new_norm_err = new_err / (new_err.sum() + 1e-30)
    ax.fill_between(centers,
                    old_norm - old_norm_err, old_norm + old_norm_err,
                    alpha=0.20, color='steelblue', step='mid')
    ax.fill_between(centers,
                    new_norm - new_norm_err, new_norm + new_norm_err,
                    alpha=0.20, color='tomato', step='mid')

    ax.set_xlabel(HIST_LABEL.get(hist_name, hist_name), fontsize=10)
    ax.set_ylabel("Normalized entries", fontsize=10)
    ax.set_title(hist_name, fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.tick_params(labelsize=9)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=0)

    differs = (p_val < KS_THRESHOLD)
    color   = 'red' if differs else 'green'
    if differs:
        ax.set_facecolor('#fff0f0')

    ax.text(0.97, 0.95,
            f"KS = {ks_stat:.3f}\np  = {p_val:.2e}",
            ha='right', va='top', transform=ax.transAxes,
            fontsize=8.5, color=color,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=color, alpha=0.9, lw=1.2))


# --------------------------------------------------------------- main loop --
for sample in SAMPLES:
    label    = sample["label"]
    old_path = sample["old_file"]
    new_path = sample["new_file"]
    outname  = sample["outname"]

    print()
    print("=" * 65)
    print(f"  Sample: {label}")
    print(f"  OLD: {old_path}")
    print(f"  NEW: {new_path}")
    print("=" * 65)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        f"Shape Comparison — {label}\n"
        f"blue = old   |   red dashed = new   |   "
        f"Red panel = shapes differ (KS p < {KS_THRESHOLD})",
        fontsize=11, fontweight='bold'
    )

    for col, hname in enumerate(HISTOGRAMS):
        ax = axes[col]
        print(f"\n  Histogram: {hname}")

        old_centers, old_contents, old_errors = th1_to_numpy(old_path, hname)
        new_centers, new_contents, new_errors = th1_to_numpy(new_path, hname)

        if old_contents is None or new_contents is None:
            ax.text(0.5, 0.5, "Load Error",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=13, color='red')
            ax.set_title(hname, fontsize=11, fontweight='bold')
            continue

        # Basic stats
        old_integral = old_contents.sum()
        new_integral = new_contents.sum()
        old_mean = np.average(old_centers, weights=old_contents) if old_integral > 0 else 0
        new_mean = np.average(new_centers, weights=new_contents) if new_integral > 0 else 0
        old_std  = np.sqrt(np.average((old_centers - old_mean)**2, weights=old_contents)) if old_integral > 0 else 0
        new_std  = np.sqrt(np.average((new_centers - new_mean)**2, weights=new_contents)) if new_integral > 0 else 0

        print(f"    old : integral={old_integral:.2f}  mean={old_mean:.2f}  std={old_std:.2f}")
        print(f"    new : integral={new_integral:.2f}  mean={new_mean:.2f}  std={new_std:.2f}")

        # KS test on expanded arrays
        old_vals = expand_to_values(old_centers, old_contents)
        new_vals = expand_to_values(new_centers, new_contents)

        if len(old_vals) == 0 or len(new_vals) == 0:
            print(f"    WARNING: one histogram is empty, skipping KS test")
            ax.text(0.5, 0.5, "Empty histogram",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=11, color='orange')
            ax.set_title(hname, fontsize=11, fontweight='bold')
            continue

        ks_stat, p_val = ks_2samp(old_vals, new_vals)

        # Normalized contents for plotting
        old_norm_arr = norm(old_contents)
        new_norm_arr = norm(new_contents)

        plot_panel(ax, old_centers, old_norm_arr, new_norm_arr,
                   old_errors, new_errors, hname, ks_stat, p_val)

        if p_val < KS_THRESHOLD:
            print(f"    *** SHAPE DIFFERS ***")
            print(f"        KS stat        = {ks_stat:.4f}")
            print(f"        p-value        = {p_val:.4e}")
            print(f"        Delta mean     = {abs(old_mean - new_mean):.3f}")
            print(f"        Delta std      = {abs(old_std  - new_std ):.3f}")
            print(f"        Delta integral = {abs(old_integral - new_integral):.3f}  "
                  f"( {100*abs(old_integral-new_integral)/max(old_integral,1e-9):.2f}% )")
        else:
            print(f"    Shape consistent.  KS={ks_stat:.4f}  p={p_val:.4e}")

    plt.tight_layout()
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  >> Saved: {outname}")

print("\nAll done.")
