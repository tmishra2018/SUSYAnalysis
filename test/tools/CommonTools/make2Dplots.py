import ROOT
import os

ROOT.gROOT.SetBatch(True)

# Output directory
out_dir = "2DPlots"
os.makedirs(out_dir, exist_ok=True)


# --- VGamma MC Samples ---
vgamma_files = {
    "2016preVFP": [
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WWG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTGJets_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WW_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTJets_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_2016preVFP.root",
    ],
    "2016postVFP": [
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WWG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTGJets_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WW_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTJets_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_2016postVFP.root",
    ],
    "2017": [
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WWG_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTGJets_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WW_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTJets_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZG_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_2017.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_2017.root",
    ],
    "2018": [
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WWG_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTGJets_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WW_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_TTJets_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZG_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet130_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGJet40_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_ZGToLLG_2018.root",
        "/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WGToLNuG_2018.root",
    ]
}

for year, files in vgamma_files.items():
    egChain = ROOT.TChain("egTree")
    mgChain = ROOT.TChain("mgTree")

    for f in files:
        egChain.Add(f)
        mgChain.Add(f)

    base_selection = "(phoEt >= 35 && abs(phoEta) <= 1.4442)"

    if year == "2017":
        pho_veto = "!(phoEta > 1.0 && phoEta < 1.5 && phoPhi > 2.7)"
        lep_veto = "!(lepEta > 1.0 && lepEta < 1.5 && lepPhi > 2.7)"
    elif year == "2018":
        pho_veto = "!(phoEta > 0.3 && phoEta < 1.2 && phoPhi > 0.4 && phoPhi < 0.8)"
        lep_veto = "!(lepEta > 0.3 && lepEta < 1.2 && lepPhi > 0.4 && lepPhi < 0.8)"
    else:
        pho_veto = lep_veto = "1"

    pho_selection = f"({base_selection}) && ({pho_veto})"
    lep_selection = f"({base_selection}) && ({lep_veto})"

    for channel, chain in [("eg", egChain), ("mg", mgChain)]:
        tag = f"MC_{'eg' if channel == 'eg' else 'mg'}_{year}"
        
        hist_pho = ROOT.TH2F(f"phoEta_phoPhi_{tag}", "", 60, -3, 3, 70, -3.5, 3.5)
        hist_lep = ROOT.TH2F(f"lepEta_lepPhi_{tag}", "", 60, -3, 3, 70, -3.5, 3.5)

        chain.Draw("phoPhi:phoEta>>" + hist_pho.GetName(), pho_selection, "goff")
        chain.Draw("lepPhi:lepEta>>" + hist_lep.GetName(), lep_selection, "goff")

        c1 = ROOT.TCanvas("c1", "", 800, 700)
        hist_pho.SetStats(0)
        hist_lep.SetStats(0)

        hist_pho.Draw("COLZ")
        c1.SaveAs(f"{out_dir}/{tag}_phoEta_phoPhi.png")

        hist_lep.Draw("COLZ")
        c1.SaveAs(f"{out_dir}/{tag}_lepEta_lepPhi.png")

        print(f"Saved VGamma plots for {tag}")

# --- Data Signal Samples ---
files_info = {
    "eg_2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
    "eg_2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
    "eg_2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
    "eg_2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root",
    "mg_2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP.root",
    "mg_2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP.root",
    "mg_2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017.root",
    "mg_2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root"
}

for tag, filepath in files_info.items():
    file = ROOT.TFile.Open(filepath)
    if not file or file.IsZombie():
        print(f"Could not open file: {filepath}")
        continue

    tree = file.Get("signalTree")
    if not tree:
        print(f"No signalTree in file: {filepath}")
        continue

    channel, year = tag.split("_")
    out_tag = f"DATA_{'eg' if channel == 'eg' else 'mg'}_{year}"
    base_selection = "(phoEt >= 35 && abs(phoEta) <= 1.4442)"

    if year == "2017":
        pho_veto = "!(phoEta > 1.0 && phoEta < 1.5 && phoPhi > 2.7)"
        lep_veto = "!(lepEta > 1.0 && lepEta < 1.5 && lepPhi > 2.7)"
    elif year == "2018":
        pho_veto = "!(phoEta > 0.3 && phoEta < 1.2 && phoPhi > 0.4 && phoPhi < 0.8)"
        lep_veto = "!(lepEta > 0.3 && lepEta < 1.2 && lepPhi > 0.4 && lepPhi < 0.8)"
    else:
        pho_veto = lep_veto = "1"

    pho_selection = f"({base_selection}) && ({pho_veto})"
    lep_selection = f"({base_selection}) && ({lep_veto})"

    hist_pho = ROOT.TH2F(f"phoEta_phoPhi_{out_tag}", "", 60, -3, 3, 70, -3.5, 3.5)
    hist_lep = ROOT.TH2F(f"lepEta_lepPhi_{out_tag}", "", 60, -3, 3, 70, -3.5, 3.5)

    hist_pho.GetXaxis().SetTitle("Photon #eta")
    hist_pho.GetYaxis().SetTitle("Photon #phi")
    hist_lep.GetXaxis().SetTitle("Lepton #eta")
    hist_lep.GetYaxis().SetTitle("Lepton #phi")


    tree.Draw("phoPhi:phoEta>>" + hist_pho.GetName(), pho_selection, "goff")
    tree.Draw("lepPhi:lepEta>>" + hist_lep.GetName(), lep_selection, "goff")

    c1 = ROOT.TCanvas("c1", "", 800, 700)
    hist_pho.SetStats(0)
    hist_lep.SetStats(0)

    hist_pho.Draw("COLZ")
    c1.SaveAs(f"{out_dir}/{out_tag}_phoEta_phoPhi.png")

    hist_lep.Draw("COLZ")
    c1.SaveAs(f"{out_dir}/{out_tag}_lepEta_lepPhi.png")

    file.Close()
    print(f"Saved Data plots for {out_tag}")

print("All plots generated successfully.")
