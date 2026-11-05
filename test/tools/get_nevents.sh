#!/bin/bash

# List of datasets and corresponding root files
declare -A datasets
datasets=(
["/MuonEG/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016B-ver1.root"
["/MuonEG/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016B-ver2.root"
["/MuonEG/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016C.root"
["/MuonEG/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016D.root"
["/MuonEG/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016E.root"
["/MuonEG/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016F_preVFP.root"
["/MuonEG/Run2016F-UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016F_postVFP.root"
["/MuonEG/Run2016G-UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016G.root"
["/MuonEG/Run2016H-UL2016_MiniAODv2-v2/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2016H.root"
["/MuonEG/Run2017B-UL2017_MiniAODv2-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2017B.root"
["/MuonEG/Run2017C-UL2017_MiniAODv2-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2017C.root"
["/MuonEG/Run2017D-UL2017_MiniAODv2-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2017D.root"
["/MuonEG/Run2017E-UL2017_MiniAODv2-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2017E.root"
["/MuonEG/Run2017F-UL2017_MiniAODv2-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2017F.root"
["/MuonEG/Run2018A-UL2018_MiniAODv2_GT36-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2018A.root"
["/MuonEG/Run2018B-UL2018_MiniAODv2_GT36-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2018B.root"
["/MuonEG/Run2018C-UL2018_MiniAODv2_GT36-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2018C.root"
["/MuonEG/Run2018D-UL2018_MiniAODv2_GT36-v1/MINIAOD"]="/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/MuonEG/MuonEG_2018D.root"
)

echo "Comparing DAS event count with ROOT file entries:"
echo "--------------------------------------------------------"

# Loop over each dataset
for dataset in "${!datasets[@]}"; do
    rootfile="${datasets[$dataset]}"

    # Check if the ROOT file exists
    if [[ ! -f "$rootfile" ]]; then
        echo "ERROR: ROOT file not found: $rootfile"
        continue
    fi

    # Query nevents from DAS
    das_nevents=$(dasgoclient --query="summary dataset=$dataset" | grep -o '"nevents":[0-9]*' | awk -F: '{print $2}')

    # Query nevents from ROOT file using ROOT batch mode
    root_nevents=$(root -l -b -q <<EOF | tail -1
TFile *f = TFile::Open("$rootfile");
if (!f || f->IsZombie()) { cout << "ERROR: Cannot open ROOT file" << endl; } 
else {
    TTree *t = (TTree*) f->Get("ggNtuplizer/EventTree");
    if (!t) { cout << "ERROR: TTree not found" << endl; } 
    else { cout << t->GetEntries() << endl; }
}
f->Close();
EOF
)

    # If root_nevents contains "ERROR", skip comparison
    if [[ "$root_nevents" == *"ERROR"* ]]; then
        echo "Skipping dataset due to error in ROOT file."
        continue
    fi

    # Compute difference
    diff=$((das_nevents - root_nevents))

    # Print results
    echo "Dataset: $dataset"
    echo "DAS nevents: $das_nevents"
    echo "ROOT file nevents: $root_nevents"
    echo "Difference: $diff"
    echo "--------------------------------------------------------"
done

