#!/usr/bin/env python3

import ROOT
import argparse
import os

ROOT.gStyle.SetOptFit(0)
ROOT.gROOT.SetBatch(True)

# Argument parsing
parser = argparse.ArgumentParser(description="Fit fake and control photon distributions.")
parser.add_argument("dataset", type=int, choices=[1, 2], help="Dataset: 1 for DoubleEG, 2 for MuonEG")
parser.add_argument("year", choices=["2016preVFP", "2016postVFP", "2017", "2018"], help="Year of data")

args = parser.parse_args()

dataset_map = {1: "DoubleEG", 2: "MuonEG"}
dataset_str = dataset_map[args.dataset]
year_str = args.year

# Construct file path
file_path = f"/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_{dataset_str}_{year_str}.root"

if not os.path.exists(file_path):
    raise FileNotFoundError(f"Input file not found: {file_path}")

f = ROOT.TFile(file_path, "READ")
g_fake = f.Get("new_fakesPhoEt")
g_ctrl = f.Get("new_controlPhoEt")

c = ROOT.TCanvas("c", "Fake vs Control Photon pT", 800, 600)
c.SetLogy()

g_fake.SetMarkerStyle(20)
g_fake.SetMarkerColor(ROOT.kRed)
g_fake.SetTitle(";p_{T} [GeV];Events/GeV")
g_ctrl.SetMarkerStyle(24)
g_ctrl.SetMarkerColor(ROOT.kBlue)

g_fake.GetXaxis().SetLimits(0, 300)
g_fake.SetMinimum(0.2)
g_fake.SetMaximum(9999)
g_fake.Draw("AP")
g_ctrl.Draw("P same")

fit_name = "exp2c"
fit = ROOT.TF1(
    fit_name,
    "[0]*exp(-[2]*x) + [1]*exp(-[3]*(x-[4]-200))  + [5]",
    20, 600
)
fit.SetParameters(51521, 1.13974, 0.0560658, 0.0200371, 125.136, 0.0970124)

g_fake.Fit(fit_name, "R")
f_fake = g_fake.GetFunction(fit_name).Clone("f_fake")
f_fake.SetLineColor(ROOT.kRed)
f_fake.SetLineWidth(2)

# Reset parameters before second fit (optional)
fit.SetParameters(1e4, 1e3, 0.05, 0.01, 1.0, 0.1)

# Fit control photon graph
g_ctrl.Fit(fit_name, "R")
f_ctrl = g_ctrl.GetFunction(fit_name).Clone("f_ctrl")
f_ctrl.SetLineColor(ROOT.kBlue)
f_ctrl.SetLineWidth(2)

# Draw fits
f_fake.Draw("same")
f_ctrl.Draw("same")

# Add legend
leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
leg.SetBorderSize(0)
leg.AddEntry(g_fake, "fake photons", "p")
leg.AddEntry(f_fake, "fit fake", "l")
leg.AddEntry(g_ctrl, "control photons", "p")
leg.AddEntry(f_ctrl, "fit control", "l")
leg.Draw()

# Save plot
output_name = f"{dataset_str}_{year_str}.png"
c.SaveAs(output_name)
c.Update()
print(f"Saved: {output_name}")
