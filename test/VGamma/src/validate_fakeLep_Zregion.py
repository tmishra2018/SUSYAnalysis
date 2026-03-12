#!/usr/bin/env python3
"""
validate_fakeLep_Zregion.py
============================
Validates the fake lepton background estimate (fakeLepTree) in the
mg (MuonEG) channel by comparing its kinematic distributions against
signalTree and proxyTree in the same data file.

Input files:
  root://cmseos.fnal.gov//eos/uscms/store/user/tmishra/eg_mg_treesData/
  resTree_mgsignal_MuonEG_{year}.root   <- mg channel (this script)
  resTree_egsignal_DoubleEG_{year}.root <- eg channel (for comparison)

All trees confirmed present: signalTree, proxyTree, fakeLepTree,
jetTree, doubleFakeTree, hadronTree.

Confirmed branches in fakeLepTree:
  phoEt, phoEta, phoPhi, lepPt, lepEta, lepPhi,
  fakeLepMiniIso, fakeLepIsStandardProxy,
  sigMT, sigMET, sigMETPhi, dPhiLepMET, threeMass,
  nVertex, dRPhoLep, HT, nJet, nBJet, run, event, lumis

PHYSICS NOTE on Z-mass veto:
  threeMass = M(gamma + mu). Z cannot decay to gamma+mu at tree level.
  There is NO Z peak in M(gamma+mu). We do NOT apply any Z-mass cut
  on this variable. The Z-veto in the signal selection (which vetoes
  M(signal_pho + signal_mu) ~ 91 GeV) removes FSR events from Z->mumu
  where the photon fakes the signal photon — that veto lives in
  analysis_mg.C on the signalTree fill, not here.

Strategy:
  Part 1 -- Shape comparison: fakeLepTree vs signalTree
            phoEt, sigMET, HT with kinematic cuts only.
            Shape differences = kinematic bias in fake muon rate.

  Part 2 -- M(gamma+mu) spectrum: signalTree vs fakeLepTree.
            No Z peak expected. Any excess in fakeLepTree at any
            mass = real Wgamma/Zgamma contamination from missing
            hasPho guard in analysis_mg.C.

  Part 3 -- fakeLepIsStandardProxy flag distribution.
            Checks what fraction of fake leptons pass the standard
            proxy definition vs the looser miniIso-based definition.

  Part 4 -- proxyTree vs signalTree phoEt.
            Spike at 35-50 GeV in proxyTree = missing HLT selection
            on proxy photon (commented out in analysis_mg.C).

  Part 5 -- eg vs mg comparison: fakeLepTree phoEt/MET/HT shapes
            between the two channels to quantify the asymmetry.
"""

import ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

# ─────────────────────────────────────────────────────────────────────────────
# USER SETTINGS — paths confirmed from .ls output
# ─────────────────────────────────────────────────────────────────────────────
EOS_BASE = "root://cmseos.fnal.gov//eos/uscms/store/user/tmishra/eg_mg_treesData"

YEARS = ["2016preVFP", "2016postVFP", "2017", "2018"]
YEAR  = "2018"   # primary year for detailed plots; change as needed

MG_FILE = f"{EOS_BASE}/resTree_mgsignal_MuonEG_{YEAR}.root"
EG_FILE = f"{EOS_BASE}/resTree_egsignal_DoubleEG_{YEAR}.root"

# Binning (same as analysis framework)
MET_BINS  = [0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000]
HT_BINS   = [0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,
             420,500,600,1000]
ET_BINS   = [35,50,100,150,200,250,300,500,800]
MASS_BINS = list(range(0, 300, 5))

# Basic kinematic cuts (mirror signal selection, no Z-mass veto on threeMass)
PHO_ET_MIN = 35.0
LEP_PT_MIN = 25.0

OUT_DIR = "fakeLep_validation_plots"
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def load_tree(fpath, treename, branches, max_events=-1):
    """Load TTree branches into numpy arrays. Supports XRootD paths."""
    print(f"  Opening {treename} from {fpath.split('/')[-1]} ...")
    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
        print(f"  [ERROR] Cannot open {fpath}")
        return None
    t = f.Get(treename)
    if not t:
        print(f"  [ERROR] TTree '{treename}' not found")
        f.Close()
        return None
    n = t.GetEntries()
    if max_events > 0:
        n = min(n, max_events)
    print(f"  Reading {n} entries ...")
    data = {b: np.zeros(n, dtype=np.float32) for b in branches}
    for i, ev in enumerate(t):
        if i >= n:
            break
        for b in branches:
            data[b][i] = getattr(ev, b, np.nan)
    f.Close()
    print(f"  Done.")
    return data


def to_density(h, bin_arr):
    bw = np.diff(bin_arr).astype(float)
    s  = float(h.sum())
    return h / s / bw if s > 0 else np.zeros_like(h, dtype=float)


def base_mask(data):
    """Kinematic selection only. No Z-mass veto on M(gamma+mu)."""
    return (
        (data["phoEt"] > PHO_ET_MIN) &
        (data["lepPt"] > LEP_PT_MIN) &
        np.isfinite(data["phoEt"]) &
        np.isfinite(data["lepPt"])
    )


def plot_two(ax, bin_arr, vals_a, vals_b, label_a, label_b,
             color_a='#2166ac', color_b='#d6604d',
             xlabel='', title='', yscale='log'):
    ha, _ = np.histogram(vals_a, bins=bin_arr)
    hb, _ = np.histogram(vals_b, bins=bin_arr)
    ax.step(bin_arr[:-1], to_density(ha, bin_arr), where='post',
            color=color_a, lw=1.8, label=label_a)
    ax.step(bin_arr[:-1], to_density(hb, bin_arr), where='post',
            color=color_b, lw=1.8, ls='--', label=label_b)
    ax.set_yscale(yscale)
    ax.set_xlim(bin_arr[0], bin_arr[-1])
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    # KS test
    if len(vals_a) > 1 and len(vals_b) > 1:
        _, p = stats.ks_2samp(vals_a, vals_b)
        col = 'darkred' if p < 0.05 else 'darkgreen'
        ax.text(0.97, 0.95, f"KS p={p:.3f}", transform=ax.transAxes,
                ha='right', va='top', fontsize=8, color=col,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        if p < 0.05:
            ax.patch.set_facecolor('#fff5f5')


# ─────────────────────────────────────────────────────────────────────────────
# Load data — mg channel
# ─────────────────────────────────────────────────────────────────────────────
branches = ["phoEt", "lepPt", "sigMET", "HT", "threeMass",
            "dRPhoLep", "nJet", "nBJet"]

print(f"\n{'='*60}")
print(f"Loading mg channel trees from {YEAR}")
print(f"{'='*60}")

sig_mg   = load_tree(MG_FILE, "signalTree",  branches)
fake_mg  = load_tree(MG_FILE, "fakeLepTree",
                     branches + ["fakeLepMiniIso", "fakeLepIsStandardProxy"])
proxy_mg = load_tree(MG_FILE, "proxyTree",   branches)

sig_mask_mg   = base_mask(sig_mg)   if sig_mg   else None
fake_mask_mg  = base_mask(fake_mg)  if fake_mg  else None
proxy_mask_mg = base_mask(proxy_mg) if proxy_mg else None

print(f"\n  signalTree   events passing cuts: {sig_mask_mg.sum() if sig_mask_mg is not None else 'N/A'}")
print(f"  fakeLepTree  events passing cuts: {fake_mask_mg.sum() if fake_mask_mg is not None else 'N/A'}")
print(f"  proxyTree    events passing cuts: {proxy_mask_mg.sum() if proxy_mask_mg is not None else 'N/A'}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 1: Shape comparison — fakeLepTree vs signalTree (mg)
#   phoEt, sigMET, HT
#   No Z-mass cut. Pure kinematic selection.
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 1: fakeLepTree vs signalTree shapes (mg channel)")
print(f"{'='*60}")

if sig_mg is not None and fake_mg is not None:
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5))
    fig1.suptitle(
        f"mg channel: fakeLepTree vs signalTree shape comparison\n"
        f"MuonEG {YEAR}  |  phoEt>{PHO_ET_MIN}, lepPt>{LEP_PT_MIN} GeV\n"
        f"No Z-mass cut applied (Z->gamma+mu is forbidden)",
        fontsize=10)

    sm  = sig_mask_mg
    fm  = fake_mask_mg

    n_sig  = sm.sum()
    n_fake = fm.sum()

    for ax, (var, bins, xlabel) in zip(axes1, [
        ("phoEt",  ET_BINS,  "phoEt (GeV)"),
        ("sigMET", MET_BINS, "sigMET (GeV)"),
        ("HT",     HT_BINS,  "HT (GeV)"),
    ]):
        plot_two(ax,
                 np.array(bins, dtype=float),
                 sig_mg[var][sm],
                 fake_mg[var][fm],
                 f"signalTree ({n_sig})",
                 f"fakeLepTree ({n_fake})",
                 xlabel=xlabel,
                 title=f"mg {YEAR} | {var}")

    plt.tight_layout()
    out1 = os.path.join(OUT_DIR, f"part1_fakeLep_vs_signal_shapes_mg_{YEAR}.png")
    plt.savefig(out1, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out1}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 2: M(gamma+mu) spectrum
#   No Z peak expected. Excess in fakeLepTree = real lepton contamination.
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 2: M(gamma+mu) spectrum — no Z peak expected")
print(f"{'='*60}")

if sig_mg is not None and fake_mg is not None:
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    bin_arr = np.array(MASS_BINS, dtype=float)

    sm = sig_mask_mg
    fm = fake_mask_mg

    plot_two(ax2, bin_arr,
             sig_mg["threeMass"][sm],
             fake_mg["threeMass"][fm],
             f"signalTree ({sm.sum()})",
             f"fakeLepTree ({fm.sum()})",
             xlabel="M(gamma+mu) [GeV]",
             title=(f"M(gamma+mu) distribution | mg MuonEG {YEAR}\n"
                    f"Z->gamma+mu FORBIDDEN: no Z peak expected.\n"
                    f"Excess in fakeLepTree = Wg/Zg contamination "
                    f"(missing hasPho guard)"))

    # Mark 91 GeV as reference only
    ax2.axvline(91.188, color='gray', ls=':', lw=1.0, alpha=0.6)
    ax2.text(92, ax2.get_ylim()[1] if ax2.get_ylim()[1]>0 else 1,
             "91 GeV\n(no peak\nexpected)", fontsize=7, color='gray', va='top')

    plt.tight_layout()
    out2 = os.path.join(OUT_DIR, f"part2_threeMass_mg_{YEAR}.png")
    plt.savefig(out2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out2}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 3: fakeLepIsStandardProxy flag
#   1 = passes standard proxy definition (inverted tight ID)
#   0 = passes only loose miniIso proxy
#   A large fraction of 0s means many fake leptons come from
#   a looser selection, which may inflate the estimate.
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 3: fakeLepIsStandardProxy flag distribution")
print(f"{'='*60}")

if fake_mg is not None:
    fm = fake_mask_mg
    flag = fake_mg["fakeLepIsStandardProxy"][fm].astype(int)
    n_standard = (flag == 1).sum()
    n_loose    = (flag == 0).sum()
    n_total    = fm.sum()

    print(f"  Standard proxy (flag=1): {n_standard}  ({100*n_standard/n_total:.1f}%)")
    print(f"  Loose miniIso  (flag=0): {n_loose}   ({100*n_loose/n_total:.1f}%)")

    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
    fig3.suptitle(
        f"fakeLepTree: standard proxy (flag=1) vs loose miniIso (flag=0)\n"
        f"mg MuonEG {YEAR}  |  standard={n_standard} ({100*n_standard/n_total:.1f}%), "
        f"loose={n_loose} ({100*n_loose/n_total:.1f}%)",
        fontsize=10)

    mask_std   = fm & (fake_mg["fakeLepIsStandardProxy"] == 1)
    mask_loose = fm & (fake_mg["fakeLepIsStandardProxy"] == 0)

    for ax, (var, bins, xlabel) in zip(axes3, [
        ("phoEt",  ET_BINS,  "phoEt (GeV)"),
        ("sigMET", MET_BINS, "sigMET (GeV)"),
        ("HT",     HT_BINS,  "HT (GeV)"),
    ]):
        plot_two(ax,
                 np.array(bins, dtype=float),
                 fake_mg[var][mask_std],
                 fake_mg[var][mask_loose],
                 f"standard proxy ({mask_std.sum()})",
                 f"loose miniIso ({mask_loose.sum()})",
                 color_a='#1a9850', color_b='#d73027',
                 xlabel=xlabel,
                 title=f"fakeLepTree proxy type | {var}")

    plt.tight_layout()
    out3 = os.path.join(OUT_DIR, f"part3_proxy_flag_mg_{YEAR}.png")
    plt.savefig(out3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out3}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 4: proxyTree vs signalTree phoEt
#   Diagnoses commented-out HLT selection in analysis_mg.C:
#     //if(!itpho->passHLTSelection())continue;
#   Expected: spike at 35-50 GeV in proxyTree only
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 4: proxyTree phoEt — HLT threshold spike check")
print(f"{'='*60}")

if sig_mg is not None and proxy_mg is not None:
    fig4, axes4 = plt.subplots(1, 3, figsize=(15, 5))
    fig4.suptitle(
        f"proxyTree vs signalTree: phoEt, sigMET, HT | mg MuonEG {YEAR}\n"
        f"Spike at 35-50 GeV in proxyTree = HLT selection commented out "
        f"in analysis_mg.C\n"
        f"(//if(!itpho->passHLTSelection())continue;)",
        fontsize=10)

    sm = sig_mask_mg
    pm = proxy_mask_mg

    for ax, (var, bins, xlabel) in zip(axes4, [
        ("phoEt",  ET_BINS,  "phoEt (GeV)"),
        ("sigMET", MET_BINS, "sigMET (GeV)"),
        ("HT",     HT_BINS,  "HT (GeV)"),
    ]):
        plot_two(ax,
                 np.array(bins, dtype=float),
                 sig_mg[var][sm],
                 proxy_mg[var][pm],
                 f"signalTree ({sm.sum()})",
                 f"proxyTree ({pm.sum()})",
                 xlabel=xlabel,
                 title=f"proxy vs signal | {var}")
        if var == "phoEt":
            ax.axvspan(35, 50, alpha=0.12, color='orange',
                       label='HLT threshold (35-50 GeV)')
            ax.legend(fontsize=7)

    plt.tight_layout()
    out4 = os.path.join(OUT_DIR, f"part4_proxyTree_phoEt_mg_{YEAR}.png")
    plt.savefig(out4, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out4}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 5: eg vs mg channel comparison — fakeLepTree shapes
#   Direct quantification of the channel asymmetry seen in validation plots.
#   mg fakeLepTree should be consistent with eg fakeLepTree if the
#   fake rate method is symmetric. Any excess in mg = overestimation.
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 5: eg vs mg fakeLepTree shape comparison")
print(f"{'='*60}")

fake_eg = load_tree(EG_FILE, "fakeLepTree", branches)

if fake_mg is not None and fake_eg is not None:
    fm_mg = fake_mask_mg
    fm_eg = base_mask(fake_eg)

    print(f"  mg fakeLepTree passing cuts: {fm_mg.sum()}")
    print(f"  eg fakeLepTree passing cuts: {fm_eg.sum()}")

    fig5, axes5 = plt.subplots(1, 3, figsize=(15, 5))
    fig5.suptitle(
        f"eg vs mg channel: fakeLepTree shape comparison\n"
        f"{YEAR}  |  phoEt>{PHO_ET_MIN}, lepPt>{LEP_PT_MIN} GeV\n"
        f"Shape difference = channel-specific bias in fake lepton estimate",
        fontsize=10)

    for ax, (var, bins, xlabel) in zip(axes5, [
        ("phoEt",  ET_BINS,  "phoEt (GeV)"),
        ("sigMET", MET_BINS, "sigMET (GeV)"),
        ("HT",     HT_BINS,  "HT (GeV)"),
    ]):
        plot_two(ax,
                 np.array(bins, dtype=float),
                 fake_eg[var][fm_eg],
                 fake_mg[var][fm_mg],
                 f"eg fakeLepTree ({fm_eg.sum()})",
                 f"mg fakeLepTree ({fm_mg.sum()})",
                 color_a='#2166ac', color_b='#d6604d',
                 xlabel=xlabel,
                 title=f"eg vs mg fakeLepTree | {var}")

    plt.tight_layout()
    out5 = os.path.join(OUT_DIR, f"part5_eg_vs_mg_fakeLep_{YEAR}.png")
    plt.savefig(out5, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out5}")


# ─────────────────────────────────────────────────────────────────────────────
# PART 6: Year-by-year yield summary
#   fakeLepTree / signalTree ratio across all years for mg channel
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PART 6: fakeLepTree / signalTree yield ratio across years")
print(f"{'='*60}")

year_labels, year_ratios = [], []

for yr in YEARS:
    fpath = f"{EOS_BASE}/resTree_mgsignal_MuonEG_{yr}.root"
    sd = load_tree(fpath, "signalTree",  ["phoEt", "lepPt"])
    fd = load_tree(fpath, "fakeLepTree", ["phoEt", "lepPt"])
    if sd is None or fd is None:
        continue
    sm = base_mask(sd)
    fm = base_mask(fd)
    ratio = fm.sum() / sm.sum() if sm.sum() > 0 else 0.0
    print(f"  {yr}: signalTree={sm.sum():7d}  fakeLepTree={fm.sum():7d}  ratio={ratio:.4f}")
    year_labels.append(yr)
    year_ratios.append(ratio)

if year_ratios:
    fig6, ax6 = plt.subplots(figsize=(8, 5))
    colors6 = ['#4dac26', '#b8e186', '#f1b6da', '#d01c8b'][:len(year_labels)]
    bars = ax6.bar(year_labels, year_ratios, color=colors6,
                   edgecolor='black', lw=0.8)
    ax6.set_ylabel("fakeLepTree / signalTree yield ratio", fontsize=11)
    ax6.set_title(
        "mg channel: fake lepton yield ratio by year\n"
        "Consistent ratio across years = stable fake rate\n"
        "Large ratio or year-dependence = systematic issue",
        fontsize=10)
    for bar, r in zip(bars, year_ratios):
        ax6.text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 0.002,
                 f"{r:.3f}", ha='center', fontsize=11, fontweight='bold')
    plt.tight_layout()
    out6 = os.path.join(OUT_DIR, "part6_fakeLep_yield_ratio_years_mg.png")
    plt.savefig(out6, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out6}")


print(f"""
====================================================================
All plots saved to: {OUT_DIR}/

Diagnostic summary:
--------------------------------------------------------------------
part1 — phoEt/MET/HT: fakeLepTree vs signalTree (mg)
         KS red background = shape bias in fake muon rate

part2 — M(gamma+mu): no Z peak expected (Z->gamma+mu forbidden)
         Excess in fakeLepTree = real lepton contamination from
         missing hasPho guard in analysis_mg.C

part3 — fakeLepIsStandardProxy flag
         Large loose (flag=0) fraction = inflated estimate from
         miniIso-based proxy definition

part4 — proxyTree phoEt: spike at 35-50 GeV = HLT selection
         commented out in analysis_mg.C for proxy photon

part5 — eg vs mg fakeLepTree direct comparison
         Shape differences = channel-specific overestimation in mg

part6 — Year-by-year yield ratio
         Inconsistency across years = year-dependent systematics
====================================================================
""")
