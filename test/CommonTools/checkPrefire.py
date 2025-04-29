import ROOT

# Disable batch mode for interactive plots
ROOT.gROOT.SetBatch(True)

# List of ROOT files with corresponding year tags
root_files = {
    "2016preVFP": "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2016preVFP.root",
    "2016postVFP": "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2016postVFP.root",
    "2017": "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2017.root",
    "2018": "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5Wg_2018.root"
}
#root_files = {
#    "2016preVFP": "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016preVFP.root",
#    "2016postVFP": "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2016postVFP.root",
#    "2017": "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2017.root",
#    "2018": "/uscms/home/tmishra/nobackup/signal_trees/resTree_TChiWG_2018.root"
#}


# Trees to process
trees = ["mgTree", "egTree"]

# Variables to plot
variables = ["phoEt", "sigMET", "sigMT", "HT"]
nbins, x_min, x_max = 50, 0, 2000  # Extended x-axis range

# Loop over ROOT files and create plots for each year
for year, file in root_files.items():
    for tree_name in trees:
        histograms = {}
        channel_tag = "mg" if tree_name == "mgTree" else "eg"
        
        for var in variables:
            histograms[var] = {
                "no_weight": ROOT.TH1F(f"h_{var}_no_weight_{year}_{channel_tag}", f"{var} ({year}, {channel_tag});{var};Events", nbins, x_min, x_max),
                "weighted": ROOT.TH1F(f"h_{var}_weighted_{year}_{channel_tag}", f"{var} with L1ECALPrefire ({year}, {channel_tag});{var};Events", nbins, x_min, x_max),
                "ratio": ROOT.TH1F(f"h_{var}_ratio_{year}_{channel_tag}", f";{var};Weighted / No Weight", nbins, x_min, x_max)
            }
        
        # Open ROOT file and process tree
        f = ROOT.TFile.Open(file)
        tree = f.Get(tree_name)
        
        for event in tree:
            if event.nJet >= 1 and event.sigMT > 100 and event.sigMET > 120:
                weight = event.L1ECALPrefire
                
                for var in variables:
                    value = getattr(event, var)
                    histograms[var]["no_weight"].Fill(value)
                    histograms[var]["weighted"].Fill(value, weight)
        
        f.Close()
        
        # Draw histograms
        for var in variables:
            canvas = ROOT.TCanvas(f"c_{var}_{year}_{channel_tag}", f"Canvas for {var} ({year}, {channel_tag})", 800, 800)
            canvas.Divide(1, 2)
            
            # Upper pad (distribution)
            pad1 = canvas.cd(1)
            pad1.SetLogy()
            pad1.SetPad(0, 0.3, 1, 1)
            histograms[var]["no_weight"].SetLineColor(ROOT.kRed)
            histograms[var]["weighted"].SetLineColor(ROOT.kBlue)
            histograms[var]["no_weight"].SetStats(False)
            histograms[var]["weighted"].SetStats(False)
            
            # Set Y-axis range to allow space for legend
            max_y = max(histograms[var]["no_weight"].GetMaximum(), histograms[var]["weighted"].GetMaximum()) * 1.5
            histograms[var]["no_weight"].SetMaximum(max_y)
            histograms[var]["no_weight"].Draw("HIST")
            histograms[var]["weighted"].Draw("HIST SAME")
            
            legend = ROOT.TLegend(0.5, 0.7, 0.85, 0.85)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(histograms[var]["no_weight"], "No Weight", "l")
            legend.AddEntry(histograms[var]["weighted"], "With L1ECALPrefire", "l")
            legend.Draw()
            
            # Lower pad (ratio plot)
            pad2 = canvas.cd(2)
            pad2.SetPad(0, 0, 1, 0.3)
            
            histograms[var]["ratio"].Divide(histograms[var]["weighted"], histograms[var]["no_weight"], 1, 1, "B")
            histograms[var]["ratio"].SetLineColor(ROOT.kBlack)
            histograms[var]["ratio"].SetStats(False)
            histograms[var]["ratio"].SetMinimum(0.95)
            histograms[var]["ratio"].SetMaximum(1.01)
            
            # Increase font size of axis labels and title
            histograms[var]["ratio"].GetXaxis().SetLabelSize(0.1)
            histograms[var]["ratio"].GetYaxis().SetLabelSize(0.1)
            histograms[var]["ratio"].GetYaxis().SetTitleSize(0.12)
            
            histograms[var]["ratio"].Draw("E")

            if histograms[var]["ratio"].GetEntries() > 0:
                no_weight_integral = histograms[var]["no_weight"].Integral()

                if no_weight_integral > 0:  # Avoid division by zero
                    mean_ratio = histograms[var]["weighted"].Integral() / no_weight_integral
                    line = ROOT.TLine(x_min, mean_ratio, x_max, mean_ratio)  # Define TLine object
                    line.SetLineColor(ROOT.kRed)
                    line.SetLineStyle(2)
                    line.Draw("same")
                
                    text = ROOT.TLatex()
                    text.SetTextSize(0.1)
                    text.DrawLatexNDC(0.7, 0.2, f"Mean: {mean_ratio:.2f}")
            
            canvas.Draw()
            canvas.SaveAs(f"/eos/uscms/store/user/tmishra/Signal/{var}_{channel_tag}_{year}.pdf")  # Save with channel tag and year
