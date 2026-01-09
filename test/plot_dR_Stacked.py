import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# --------------------------------------------------------
# Input ROOT files
# --------------------------------------------------------
sm_files = [
    "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_2017.root",
    "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_2017.root",
    "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_2017.root"
]

signal_files = {
    "T5WG": "/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_2017.root",
    "TChiWG": "/uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_2017.root"
}

# --------------------------------------------------------
# Channel definitions
# --------------------------------------------------------
channels = {
    "eg": {"hname": "h_dReg", "hname_sig": "h_dReg_t5wg", "xtitle": "#DeltaR(e, #gamma)"},
    "mg": {"hname": "h_dRmg", "hname_sig": "h_dRmg_t5wg", "xtitle": "#DeltaR(#mu, #gamma)"}
}

# --------------------------------------------------------
# Colors for SM and Signal samples
# --------------------------------------------------------
sm_colors = [ROOT.kGray + 1, ROOT.kGreen + 2, ROOT.kAzure + 1]
colors_signal = {"TChiWG": ROOT.kMagenta + 2, "T5WG": ROOT.kBlue + 1}
line_styles = {"TChiWG": 2, "T5WG": 7}

# --------------------------------------------------------
# Output directory
# --------------------------------------------------------
os.makedirs("plots_dR", exist_ok=True)

# --------------------------------------------------------
# Loop over channels
# --------------------------------------------------------
for ch_key, ch_info in channels.items():
    print(f"\n=== Processing channel: {ch_key} ({ch_info['xtitle']}) ===")

    stack = ROOT.THStack(f"stack_{ch_key}", "")
    leg = ROOT.TLegend(0.55, 0.60, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)

    hist_SM_total = None

    # -----------------------------
    # Load SM histograms
    # -----------------------------
    for i, smf in enumerate(sm_files):
        f = ROOT.TFile.Open(smf)
        if not f or f.IsZombie():
            print(f"❌ Could not open file: {smf}")
            continue

        h = f.Get(ch_info["hname"])
        if not h or not isinstance(h, ROOT.TH1):
            print(f"⚠️  Histogram {ch_info['hname']} not found in {smf}")
            f.Close()
            continue

        h = h.Clone(f"h_{i}_{ch_key}")  # Make local copy
        h.SetDirectory(0)               # Detach from file
        f.Close()

        h.SetFillColor(sm_colors[i % len(sm_colors)])
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)

        stack.Add(h)
        leg.AddEntry(h, os.path.basename(smf).replace(".root", ""), "f")

        if hist_SM_total is None:
            hist_SM_total = h.Clone(f"hist_SM_total_{ch_key}")
            hist_SM_total.Reset()
            hist_SM_total.SetDirectory(0)
        hist_SM_total.Add(h)

    # -----------------------------
    # Load Signal histograms
    # -----------------------------
    sig_hists = []
    for sig_name, sig_path in signal_files.items():
        f = ROOT.TFile.Open(sig_path)
        if not f or f.IsZombie():
            print(f"❌ Could not open signal file: {sig_path}")
            continue

        # Replace 't5wg' by actual signal key to match naming
        hname_sig = ch_info["hname_sig"].replace("t5wg", sig_name.lower())
        hsig = f.Get(hname_sig)
        if not hsig or not isinstance(hsig, ROOT.TH1):
            print(f"⚠️  Histogram {hname_sig} not found in {sig_path}")
            f.Close()
            continue

        hsig = hsig.Clone(f"h_{sig_name}_{ch_key}")
        hsig.SetDirectory(0)
        f.Close()

        hsig.SetLineColor(colors_signal[sig_name])
        hsig.SetLineStyle(line_styles[sig_name])
        hsig.SetLineWidth(3)

        # Normalize to total SM yield for visibility
        if hist_SM_total and hsig.Integral() > 0:
            scale = hist_SM_total.Integral() / hsig.Integral()
            hsig.Scale(scale)

        sig_hists.append(hsig)
        leg.AddEntry(hsig, sig_name, "l")

    # -----------------------------
    # Draw and format
    # -----------------------------
    c = ROOT.TCanvas(f"c_{ch_key}", "", 700, 600)
    c.SetLogy()

    stack.Draw("HIST")
    stack.GetXaxis().SetTitle(ch_info["xtitle"])
    stack.GetYaxis().SetTitle("Events")
    stack.SetTitle("")
    stack.SetMaximum(stack.GetMaximum() * 20)
    stack.SetMinimum(1)

    for hsig in sig_hists:
        hsig.Draw("HIST SAME")

    leg.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.15, 0.93, "#bf{CMS} Simulation (2017)  #it{Work in progress}")

    # -----------------------------
    # Save plots
    # -----------------------------
    out_png = f"plots_dR/DeltaR_{ch_key}_stacked.png"
    out_pdf = f"plots_dR/DeltaR_{ch_key}_stacked.pdf"
    c.SaveAs(out_png)
    c.SaveAs(out_pdf)
    print(f"✅ Saved plots: {out_png}, {out_pdf}")

print("\nAll done! Plots saved in 'plots_dR/' directory.\n")
