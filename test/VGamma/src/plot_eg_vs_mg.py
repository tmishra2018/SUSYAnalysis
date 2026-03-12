#!/usr/bin/env python3
"""
plot_eg_vs_mg.py
================
Produces two sets of plots:

PART 1 — eg vs mg channel comparison (48 plots total):
  For each of 4 TTrees (signalTree, proxyTree, jetTree, fakeLepTree):
    For each of 4 years (2016preVFP, 2016postVFP, 2017, 2018):
      3-panel figure comparing eg vs mg for sigMET, HT, phoEt
    → 4 trees × 4 years = 16 figures × 3 panels = 48 plots

PART 2 — Year comparison per channel (6 plots total):
  For eg channel: sigMET, HT, phoEt overlaying all 4 years  (3 plots)
  For mg channel: sigMET, HT, phoEt overlaying all 4 years  (3 plots)
  → Uses signalTree

Run with:
  python3 plot_eg_vs_mg.py   (inside CMSSW / PyROOT environment)
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import ks_2samp

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
except ImportError:
    print("ERROR: PyROOT not available. Source CMSSW environment first.")
    sys.exit(1)

# ------------------------------------------------------------------ config ---
XRD  = "root://cmseos.fnal.gov/"
BASE = "/eos/uscms/store/user/tmishra/eg_mg_treesData/"

YEARS = ["2016preVFP", "2016postVFP", "2017", "2018"]

EG_FILES = {y: f"{XRD}/{BASE}resTree_egsignal_DoubleEG_{y}.root"  for y in YEARS}
MG_FILES = {y: f"{XRD}/{BASE}resTree_mgsignal_MuonEG_{y}.root"    for y in YEARS}

TREES = ["signalTree", "proxyTree", "jetTree", "fakeLepTree"]

VARIABLES = ["sigMET", "HT", "phoEt"]
VAR_LABEL  = {"sigMET": "MET [GeV]", "HT": "H_{T} [GeV]", "phoEt": "Photon E_{T} [GeV]"}

MET_BINS = np.array([0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000], dtype=float)
HT_BINS  = np.array([0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000], dtype=float)
ET_BINS  = np.array([35,50,100,150,200,250,300,500,800], dtype=float)

VAR_EDGES = {"sigMET": MET_BINS, "HT": HT_BINS, "phoEt": ET_BINS}

KS_THRESHOLD = 0.05

# colours
EG_COL = "#2166ac"   # blue
MG_COL = "#d6604d"   # red/orange

YEAR_COLS  = {"2016preVFP":  "#1b7837",
              "2016postVFP": "#762a83",
              "2017":        "#d6604d",
              "2018":        "#2166ac"}
YEAR_LINES = {"2016preVFP": "-", "2016postVFP": "--", "2017": "-.", "2018": ":"}

def branch_to_numpy(fpath, tree_name, branch_name):
    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
        print(f"    ERROR: cannot open {fpath}")
        return None
    t = f.Get(tree_name)
    if not t:
        print(f"    ERROR: tree '{tree_name}' not found in {fpath}")
        f.Close()
        return None
    arr = []
    for ev in t:
        v = getattr(ev, branch_name, None)
        if v is not None:
            arr.append(float(v))
    f.Close()
    if not arr:
        print(f"    WARNING: '{branch_name}' empty in {tree_name} of {fpath}")
        return None
    return np.array(arr, dtype=np.float32)


def make_hist(arr, edges):
    """Return (centers, norm_density, counts, widths) using variable-width bin edges.
    Divides by bin width so the y-axis is a proper density for shape comparison."""
    counts, _ = np.histogram(arr, bins=edges)
    centers   = 0.5 * (edges[:-1] + edges[1:])
    widths    = np.diff(edges)
    total     = counts.sum()
    norm = (counts / total / widths) if total > 0 else np.zeros_like(counts, dtype=float)
    return centers, norm, counts, widths


def ks_test(arr1, arr2):
    if arr1 is None or arr2 is None or len(arr1) == 0 or len(arr2) == 0:
        return None, None
    stat, pval = ks_2samp(arr1, arr2)
    return stat, pval


def step_plot(ax, edges, norm, color, lw=1.8, ls="-", label="", alpha=1.0):
    """Draw a step histogram given bin edges and normalized density values."""
    ax.step(edges, np.append(norm, norm[-1]),
            where='post', color=color, lw=lw, ls=ls, label=label, alpha=alpha)


def annotate_ks(ax, ks_stat, p_val):
    differs = p_val is not None and p_val < KS_THRESHOLD
    col = 'red' if differs else 'green'
    if differs:
        ax.set_facecolor('#fff0f0')
    txt = (f"KS={ks_stat:.3f}\np={p_val:.1e}" if ks_stat is not None
           else "N/A")
    ax.text(0.97, 0.95, txt, ha='right', va='top',
            transform=ax.transAxes, fontsize=7.5, color=col,
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=col, alpha=0.9))


def set_auto_range(ax, edges, *norm_arrays, log=True):
    """Set x and y limits based on where data actually lives.
    x : clip to last bin with non-zero content (+1 bin margin).
    y : log scale (default) — floor at 0.3 x min positive value,
                               ceiling at 5 x peak.
        linear fallback      — 0 to 1.30 x peak."""
    last_nonzero = 0
    peak      = 0.0
    floor_val = np.inf
    for norm in norm_arrays:
        nz = np.where(norm > 0)[0]
        if len(nz):
            last_nonzero = max(last_nonzero, nz[-1])
        if norm.max() > peak:
            peak = norm.max()
        pos = norm[norm > 0]
        if len(pos) and pos.min() < floor_val:
            floor_val = pos.min()

    x_max_idx = min(last_nonzero + 2, len(edges) - 1)
    ax.set_xlim(edges[0], edges[x_max_idx])

    if log and peak > 0:
        ax.set_yscale("log")
        y_min = max(floor_val * 0.3, peak * 1e-5)
        ax.set_ylim(y_min, peak * 5.0)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.LogFormatterSciNotation(labelOnlyBase=False))
    else:
        ax.set_yscale("linear")
        ax.set_ylim(0, peak * 1.30 if peak > 0 else 1.0)


# ============================================================
# PART 1 — eg vs mg per tree per year (48 plots → 16 figures)
# ============================================================
print("\n" + "="*65)
print("  PART 1: eg vs mg comparison")
print("="*65)

for tree in TREES:
    for year in YEARS:
        eg_path = EG_FILES[year]
        mg_path = MG_FILES[year]

        print(f"\n  [{tree}]  {year}")

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        fig.suptitle(
            f"eg  vs  mg  —  {tree}  |  {year}\n"
            f"blue=eg   red dashed=mg   |   red panel = KS p < {KS_THRESHOLD}",
            fontsize=10, fontweight='bold')

        for col, var in enumerate(VARIABLES):
            ax = axes[col]
            edges = VAR_EDGES[var]

            eg_arr = branch_to_numpy(eg_path, tree, var)
            mg_arr = branch_to_numpy(mg_path, tree, var)

            if eg_arr is None or mg_arr is None:
                ax.text(0.5, 0.5, "Load Error", ha='center', va='center',
                        transform=ax.transAxes, color='red', fontsize=12)
                ax.set_title(var, fontweight='bold')
                continue

            eg_c, eg_n, _, _ = make_hist(eg_arr, edges)
            mg_c, mg_n, _, _ = make_hist(mg_arr, edges)

            step_plot(ax, edges, eg_n, EG_COL, label=f"eg  (n={len(eg_arr):,})")
            step_plot(ax, edges, mg_n, MG_COL, ls='--', label=f"mg  (n={len(mg_arr):,})")

            ks_stat, p_val = ks_test(eg_arr, mg_arr)
            annotate_ks(ax, ks_stat, p_val)

            ax.set_xlabel(VAR_LABEL[var], fontsize=9)
            ax.set_ylabel("Entries / bin width (log, normalized)", fontsize=9)
            ax.set_title(f"{var}", fontsize=10, fontweight='bold')
            ax.legend(fontsize=8)
            ax.tick_params(labelsize=8)
            set_auto_range(ax, edges, eg_n, mg_n)

            if p_val is not None and p_val < KS_THRESHOLD:
                print(f"    *** {var}: SHAPE DIFFERS  KS={ks_stat:.4f}  "
                      f"p={p_val:.2e}  "
                      f"Δmean={abs(eg_arr.mean()-mg_arr.mean()):.2f}")
            else:
                print(f"    {var}: consistent  KS={ks_stat:.4f}  p={p_val:.2e}"
                      if ks_stat is not None else f"    {var}: skipped")

        plt.tight_layout()
        outname = f"/eos/uscms/store/user/tmishra/part1_{tree}_{year}.png"
        plt.savefig(outname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    >> Saved: {outname}")


# ============================================================
# PART 2 — Year comparison per channel, all 4 TTrees
#           2 channels × 4 trees × 3 variables = 24 plots
# ============================================================
print("\n" + "="*65)
print("  PART 2: Year comparison per channel (all 4 TTrees)")
print("="*65)

for channel, files in [("eg", EG_FILES), ("mg", MG_FILES)]:
    print(f"\n  Channel: {channel}")

    for tree in TREES:
        print(f"\n    Tree: {tree}")

        for var in VARIABLES:
            edges = VAR_EDGES[var]

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.set_title(
                f"{channel.upper()} channel — {var}  |  {tree}\nAll years overlaid",
                fontsize=11, fontweight='bold')
            ax.set_xlabel(VAR_LABEL[var], fontsize=10)
            ax.set_ylabel("Entries / bin width (log, normalized)", fontsize=10)
            ax.tick_params(labelsize=9)

            arrays   = {}
            all_norms = []
            for year in YEARS:
                arr = branch_to_numpy(files[year], tree, var)
                if arr is None:
                    continue
                arrays[year] = arr
                _, norm, _, _ = make_hist(arr, edges)
                all_norms.append(norm)
                step_plot(ax, edges, norm,
                          color=YEAR_COLS[year],
                          ls=YEAR_LINES[year],
                          lw=2.0,
                          label=f"{year}  (n={len(arr):,})")

            if all_norms:
                set_auto_range(ax, edges, *all_norms)
            else:
                ax.set_xlim(edges[0], edges[-1])
                ax.set_ylim(bottom=0)
            ax.legend(fontsize=9, framealpha=0.9)

            # KS between consecutive years
            year_list = [y for y in YEARS if y in arrays]
            for i in range(len(year_list) - 1):
                ya, yb = year_list[i], year_list[i+1]
                ks, p = ks_test(arrays[ya], arrays[yb])
                flag = "*** DIFFERS ***" if p < KS_THRESHOLD else "consistent"
                print(f"      {var}  {ya} vs {yb}:  KS={ks:.4f}  p={p:.2e}  {flag}")

            plt.tight_layout()
            outname = f"/eos/uscms/store/user/tmishra/part2_{channel}_{tree}_{var}_years.png"
            plt.savefig(outname, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"      >> Saved: {outname}")

print("\n" + "="*65)
print("  All done.")
print("="*65 + "\n")
