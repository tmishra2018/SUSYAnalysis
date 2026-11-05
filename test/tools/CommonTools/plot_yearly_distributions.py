#!/usr/bin/env python3
import os
from array import array
from collections import OrderedDict
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

YEARS = ["2016preVFP", "2016postVFP", "2017", "2018"]
YEAR_COLORS = {
    "2016preVFP": ROOT.kBlue + 1,
    "2016postVFP": ROOT.kRed + 1,
    "2017": ROOT.kGreen + 2,
    "2018": ROOT.kMagenta + 1,
}

# Requested variable-bin definitions for TTree variables
BKG_ET_BINS = array("d", [35, 50, 100, 150, 200, 250, 300, 500, 800])
BKG_MET_BINS = array("d", [0, 40, 60, 80, 100, 120, 140, 160, 180, 210, 240, 280, 320, 400, 600, 1000])
BKG_MT_BINS = array("d", [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 300, 400, 500, 1000])
BKG_HT_BINS = array("d", [0, 40, 60, 80, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300, 340, 380, 420, 500, 600, 1000])

TREE_VARIABLES = OrderedDict([
    ("phoEt", (BKG_ET_BINS, "#gamma E_{T}")),
    ("sigMET", (BKG_MET_BINS, "MET")),
    ("sigMT", (BKG_MT_BINS, "M_{T}")),
    ("HT", (BKG_HT_BINS, "HT")),
])

TREES = ["signalTree", "proxyTree", "jetTree", "fakeLepTree"]

TREE_FILES = {
    "mg": {
        "2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP.root",
        "2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP.root",
        "2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017.root",
        "2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root",
    },
    "eg": {
        "2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root",
    },
}

HIST_VARIABLES = OrderedDict([
    ("p_MET", "p_{T}^{miss}"),
    ("p_MT", "M_{T}"),
    ("p_HT", "H_{T}"),
    ("p_PhoEt", "p_{T}^{#gamma}"),
])

BACKGROUND_SAMPLES = OrderedDict([
    ("eleBkg", "e#rightarrow#gamma fake"),
    ("jetbkg", "jet#rightarrow#gamma fake"),
    ("qcd", "QCD"),
    ("VGBkg", "V#gamma"),
    ("rareBkg", "Rare"),
    ("signal", "Data"),
])

BACKGROUND_FILES = {
    "eg": "/uscms_data/d3/tmishra/Background/validTree_egamma_{sample}_{year}.root",
    "mg": "/uscms_data/d3/tmishra/Background/validTree_mg_{sample}_{year}.root",
}


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def set_positive_minimum_for_log(hist):
    # Ensure visibility on log-scale pads even for sparse distributions.
    hist.SetMinimum(1)


def draw_year_overlay(hists_by_year, title, x_title, out_path):
    canvas = ROOT.TCanvas("c", "c", 700, 600)
    canvas.SetLogy()
    legend = ROOT.TLegend(0.62, 0.72, 0.90, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    ymax = max([h.GetMaximum() for h in hists_by_year.values()] + [1.0]) * 3.0

    first = True
    for year in YEARS:
        hist = hists_by_year[year]
        hist.SetLineColor(YEAR_COLORS[year])
        hist.SetLineWidth(2)
        hist.SetMarkerColor(YEAR_COLORS[year])
        hist.SetMarkerStyle(20)
        hist.SetTitle(title)
        hist.GetXaxis().SetTitle(x_title)
        hist.GetYaxis().SetTitle("Events")
        set_positive_minimum_for_log(hist)
        hist.SetMaximum(ymax)
        draw_opt = "hist e" if first else "hist e same"
        hist.Draw(draw_opt)
        legend.AddEntry(hist, year, "l")
        first = False

    legend.Draw()
    canvas.SaveAs(out_path)


def make_tree_hist(tree, branch, bins, hist_name):
    nbins = len(bins) - 1
    hist = ROOT.TH1D(hist_name, "", nbins, bins)

    if tree:
        # Use Draw(...>>hist) so ROOT fills this explicit histogram object reliably.
        draw_expr = f"{branch}>>{hist_name}"
        tree.Draw(draw_expr, "", "goff")

    # Keep histogram alive independently of the current ROOT directory.
    hist.SetDirectory(0)
    return hist


def plot_tree_distributions(channel):
    out_dir = os.path.join("yearly_output", channel, "trees")
    ensure_dir(out_dir)

    files = {year: ROOT.TFile.Open(path) for year, path in TREE_FILES[channel].items()}

    for tree_name in TREES:
        tree_dir = os.path.join(out_dir, tree_name)
        ensure_dir(tree_dir)

        for var_name, (bins, x_title) in TREE_VARIABLES.items():
            hists = {}
            for year in YEARS:
                root_file = files[year]
                tree = root_file.Get(tree_name) if root_file else None
                hist_name = f"h_{channel}_{tree_name}_{var_name}_{year}"
                hists[year] = make_tree_hist(tree, var_name, bins, hist_name)

            draw_year_overlay(
                hists,
                f"{channel.upper()} {tree_name} : {var_name}",
                x_title,
                os.path.join(tree_dir, f"{var_name}.png"),
            )

    for root_file in files.values():
        if root_file:
            root_file.Close()


def plot_hist_distributions(channel):
    out_dir = os.path.join("yearly_output", channel, "backgrounds")
    ensure_dir(out_dir)

    for sample, sample_label in BACKGROUND_SAMPLES.items():
        sample_dir = os.path.join(out_dir, sample)
        ensure_dir(sample_dir)

        files = {
            year: ROOT.TFile.Open(BACKGROUND_FILES[channel].format(sample=sample, year=year))
            for year in YEARS
        }

        for var_name, x_title in HIST_VARIABLES.items():
            hists = {}
            for year in YEARS:
                root_file = files[year]
                source = root_file.Get(var_name) if root_file else None
                hist_name = f"h_{channel}_{sample}_{var_name}_{year}"
                if source:
                    hist = source.Clone(hist_name)
                    hist.SetDirectory(0)
                else:
                    # fallback empty histogram if source object is missing in file
                    hist = ROOT.TH1D(hist_name, "", 10, 0, 1)
                    hist.SetDirectory(0)
                    hist.Reset()
                hists[year] = hist

            draw_year_overlay(
                hists,
                f"{channel.upper()} {sample_label} : {var_name}",
                x_title,
                os.path.join(sample_dir, f"{var_name}.png"),
            )

        for root_file in files.values():
            if root_file:
                root_file.Close()


def main():
    for channel in ["eg", "mg"]:
        plot_tree_distributions(channel)
        plot_hist_distributions(channel)


if __name__ == "__main__":
    main()
