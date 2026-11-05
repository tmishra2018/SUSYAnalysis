import ROOT

# Enable multi-threading for faster reading over EOS
ROOT.ROOT.EnableImplicitMT(8)

datasets = [
    ('DoubleEG (e-gamma)', '/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/DoubleEG/DoubleEG_2018*.root'),
    ('MuonEG (mu-gamma)',   '/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2018*.root')
]

for channel, path in datasets:
    chain = ROOT.TChain('ggNtuplizer/EventTree')
    chain.Add(path)
    
    tot = chain.GetEntries()
    hem = chain.GetEntries('run >= 319077')
    frac = hem / tot if tot > 0 else 0
    
    print(f"\n{channel}:")
    print(f"  Total Events : {tot:,}")
    print(f"  HEM Events   : {hem:,}")
    print(f"  fracHEM      : {frac:.6f}")
