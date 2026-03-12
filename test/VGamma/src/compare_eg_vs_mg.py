#!/usr/bin/env python3
"""
compare_eg_vs_mg.py
===================
Compares shapes of p_MET, p_HT, p_PhoEt (TH1D stored directly in ROOT files)
between egamma and mg channels for all processes (2018).

Processes compared:
  VGBkg, qcd, eleBkg, rareBkg, jetbkg, signal

Produces one PNG per process (6 total), each with 3 panels (one per variable).
KS test flags shape mismatches with a print statement and red panel background.

Run with:
  python3 compare_eg_vs_mg.py   (inside CMSSW env where PyROOT is available)
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
    print("ERROR: PyROOT not available. Source CMSSW environment first.")
    sys.exit(1)

# ------------------------------------------------------------------ config ---
BASE = "/uscms_data/d3/tmishra/Background"
BASE_menglei = "/uscms_data/d3/tmishra/Background_Menglei"

PROCESSES = ["VGBkg", "qcd", "eleBkg", "rareBkg", "jetbkg", "signal"]

HISTOGRAMS = ["p_MET", "p_HT", "p_PhoEt"]
HIST_LABEL = {
    "p_MET"   : "MET [GeV]",
    "p_HT"    : "H_{T} [GeV]",
    "p_PhoEt" : "Photon E_{T} [GeV]",
}

KS_THRESHOLD = 0.05

EG_COLOR  = "#2166ac"   # blue  — egamma
MG_COLOR  = "#d6604d"   # red   — mg


# --------------------------------------------------------------- helpers ----
def th1_to_numpy(tfile_path, hist_name):
    """Return (centers, contents, errors) arrays from a TH1D in a ROOT file."""
    f = ROOT.TFile.Open(tfile_path)
    if not f or f.IsZombie():
        print(f"    ERROR: cannot open  {tfile_path}")
        return None, None, None
    h = f.Get(hist_name)
    if not h:
        print(f"    ERROR: '{hist_name}' not found in {tfile_path}")
        f.Close()
        return None, None, None
    nbins    = h.GetNbinsX()
    centers  = np.array([h.GetBinCenter(i)  for i in range(1, nbins + 1)])
    contents = np.array([h.GetBinContent(i) for i in range(1, nbins + 1)])
    errors   = np.array([h.GetBinError(i)   for i in range(1, nbins + 1)])
    f.Close()
    return centers, contents, errors


def safe_norm(contents):
    s = contents.sum()
    return contents / s if s > 0 else contents


def expand_to_values(centers, contents):
    """Expand binned histogram to pseudo-event array for KS test."""
    vals = []
    for c, n in zip(centers, contents):
        count = int(round(n))
        if count > 0:
            vals.extend([c] * count)
    return np.array(vals)


def weighted_stats(centers, contents):
    total = contents.sum()
    if total <= 0:
        return 0.0, 0.0, total
    mean = np.average(centers, weights=contents)
    std  = np.sqrt(np.average((centers - mean) ** 2, weights=contents))
    return mean, std, total


def plot_panel(ax, centers, eg_norm, mg_norm, eg_err, mg_err,
               hist_name, ks_stat, p_val):
    """Draw normalized overlaid histograms with shaded error bands."""
    if len(centers) < 2:
        ax.text(0.5, 0.5, "Too few bins", ha='center', va='center',
                transform=ax.transAxes)
        return

    width  = centers[1] - centers[0]
    edges  = np.append(centers - width / 2, centers[-1] + width / 2)

    # Step histograms — pad for matplotlib's step plotting
    def step_arrays(norm_arr):
        return np.append(edges[:-1], edges[-1]), np.append(norm_arr, 0)

    ex, ey = step_arrays(eg_norm)
    mx, my = step_arrays(mg_norm)

    ax.step(ex, ey, where='post', color=EG_COLOR, lw=2.0, label='egamma', zorder=3)
    ax.step(mx, my, where='post', color=MG_COLOR, lw=2.0, label='mg',
            linestyle='--', zorder=3)

    # Shaded error bands (rough normalized error propagation)
    eg_nerr = eg_err / (eg_err.sum() + 1e-30)
    mg_nerr = mg_err / (mg_err.sum() + 1e-30)
    ax.fill_between(centers, eg_norm - eg_nerr, eg_norm + eg_nerr,
                    alpha=0.18, color=EG_COLOR, step='mid')
    ax.fill_between(centers, mg_norm - mg_nerr, mg_norm + mg_nerr,
                    alpha=0.18, color=MG_COLOR, step='mid')

    ax.set_xlabel(HIST_LABEL.get(hist_name, hist_name), fontsize=10)
    ax.set_ylabel("Normalized entries", fontsize=10)
    ax.set_title(hist_name, fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.tick_params(labelsize=9)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(bottom=0)

    differs = (p_val < KS_THRESHOLD)
    ann_col = 'red' if differs else 'green'
    if differs:
        ax.set_facecolor('#fff0f0')

    ax.text(0.97, 0.95,
            f"KS = {ks_stat:.3f}\np  = {p_val:.2e}",
            ha='right', va='top', transform=ax.transAxes,
            fontsize=8.5, color=ann_col,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=ann_col, alpha=0.9, lw=1.2))


# --------------------------------------------------------------- main loop --
print("\n" + "=" * 70)
print("  egamma vs mg shape comparison")
print("=" * 70)

for proc in PROCESSES:
    eg_path = f"{BASE}/validTree_mg_{proc}_2016postVFP.root"
    mg_path = f"{BASE_menglei}/validTree_mg_{proc}.root"

    print(f"\n{'='*70}")
    print(f"  Process : {proc}")
    print(f"  egamma  : {eg_path}")
    print(f"  mg      : {mg_path}")
    print(f"{'='*70}")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        f"egamma  vs  mg  {proc}  (2018)\n"
        f"blue = egamma   |   red dashed = mg   |   "
        f"Red panel = shapes differ (KS p < {KS_THRESHOLD})",
        fontsize=11, fontweight='bold'
    )

    any_mismatch = False

    for col, hname in enumerate(HISTOGRAMS):
        ax = axes[col]
        print(f"\n  [{proc}]  {hname}")

        eg_centers, eg_contents, eg_errors = th1_to_numpy(eg_path, hname)
        mg_centers, mg_contents, mg_errors = th1_to_numpy(mg_path, hname)

        # Handle load failures gracefully
        if eg_contents is None or mg_contents is None:
            ax.text(0.5, 0.5, "Load Error",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=13, color='red')
            ax.set_title(hname, fontsize=11, fontweight='bold')
            continue

        eg_mean, eg_std, eg_int = weighted_stats(eg_centers, eg_contents)
        mg_mean, mg_std, mg_int = weighted_stats(mg_centers, mg_contents)

        print(f"    egamma : integral={eg_int:10.2f}  mean={eg_mean:8.2f}  std={eg_std:8.2f}")
        print(f"    mg     : integral={mg_int:10.2f}  mean={mg_mean:8.2f}  std={mg_std:8.2f}")

        eg_vals = expand_to_values(eg_centers, eg_contents)
        mg_vals = expand_to_values(mg_centers, mg_contents)

        if len(eg_vals) == 0 or len(mg_vals) == 0:
            print(f"    WARNING: one histogram is empty — skipping KS test")
            ax.text(0.5, 0.5, "Empty histogram",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=11, color='orange')
            ax.set_title(hname, fontsize=11, fontweight='bold')
            continue

        ks_stat, p_val = ks_2samp(eg_vals, mg_vals)

        plot_panel(ax,
                   eg_centers,
                   safe_norm(eg_contents),
                   safe_norm(mg_contents),
                   eg_errors, mg_errors,
                   hname, ks_stat, p_val)

        if p_val < KS_THRESHOLD:
            any_mismatch = True
            print(f"    *** SHAPE DIFFERS ***")
            print(f"        KS stat        = {ks_stat:.4f}")
            print(f"        p-value        = {p_val:.4e}")
            print(f"        Delta mean     = {abs(eg_mean - mg_mean):.3f}")
            print(f"        Delta std      = {abs(eg_std  - mg_std ):.3f}")
            print(f"        Delta integral = {abs(eg_int  - mg_int ):.2f}"
                  f"  ( {100*abs(eg_int-mg_int)/max(eg_int,1e-9):.2f}% relative )")
        else:
            print(f"    Shape consistent.   KS={ks_stat:.4f}   p={p_val:.4e}")

    if any_mismatch:
        print(f"\n  >> SUMMARY [{proc}]: At least one variable has a SHAPE MISMATCH")
    else:
        print(f"\n  >> SUMMARY [{proc}]: All three variables are shape-consistent")

    plt.tight_layout()
    outname = f"compare_eg_vs_mg_{proc}_2018.png"
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  >> Saved : {outname}")

print("\n" + "=" * 70)
print("  All processes done.")
print("=" * 70 + "\n")
