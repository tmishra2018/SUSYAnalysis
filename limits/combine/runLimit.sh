#!/bin/bash

export jobNumber=$1
export model=$2
export mode=$3
export model_lower=$(echo "$model" | tr '[:upper:]' '[:lower:]')
export maskBins=${4:-false}   # true -> mask the bins listed below (test run); default false
export name_list="cardsSplitted/x${jobNumber}"

#MASK_BINS=(1 2 3 4 5 6 10 11 12 13 14 15 19 20 21 22 23 24  28 29 30 31 32 33)  # MET_gt400

#MASK_BINS=(4 5 6 7 8 9 13 14 15 16 17 18 22 23 24 25 26 27  31 32 33 34 35 36)  # MET_lt200  

MASK_BINS=(1 2 3 7 8 9 10 11 12 16 17 18 19 20 21 25 26 27  28 29 30 34 35 36)  # MET_200_400

#MASK_BINS=(19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36)
#MASK_BINS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18) 
#MASK_BINS=(13 18 34 35 36) 

limitDir="limits_Run2_${mode}"

if [[ "$mode" != "splusb" && "$mode" != "unblind" ]]; then
    echo "Error: mode must be 'splusb' or 'unblind' (got '${mode}')"
    exit 1
fi

cd ${_CONDOR_SCRATCH_DIR}
ls -ltrh

# Setup CMS environment
echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "scramv1 project CMSSW CMSSW_14_1_0_pre4"
scramv1 project CMSSW CMSSW_14_1_0_pre4

echo "cd CMSSW_14_1_0_pre4/src/"
cd CMSSW_14_1_0_pre4/src/

echo "eval \`scramv1 runtime -sh\`"
eval `scramv1 runtime -sh`

# Clone HiggsAnalysis repository
echo "git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit"
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

echo "Moving cardsSplitted directory"
mv ../../cardsSplitted HiggsAnalysis/CombinedLimit

# Build CMSSW environment
scramv1 b clean
scramv1 b

echo "cd HiggsAnalysis/CombinedLimit"
cd HiggsAnalysis/CombinedLimit
echo "Processing file: $name_list"
echo "Mode: ${mode}"
echo "Mask bins: ${maskBins} (${MASK_BINS[*]})"

if [[ "$maskBins" == "true" ]]; then
    limitDir="${limitDir}_masked_MET_200_400"
fi
# Output directory is created for both masked and unmasked runs (runLimit.py also creates it at submit time)
limitPath=/store/user/tmishra/CombinedLimit/${model}/limits/${limitDir}
xrdfs root://cmseos.fnal.gov mkdir -p ${limitPath}
nFailed=0

if [[ "$maskBins" == "true" ]]; then

    cat > getMasks.py <<'PYEOF'
import sys, re, ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError
w = ROOT.TFile.Open(sys.argv[1]).Get("w")
names = [n for n in w.allVars().contentsString().split(",") if n.startswith("mask_")]
sel = []
for b in sys.argv[2:]:
    pat = re.compile(r"(^|_)bin%s(_|$)" % b)   # bin13 must not match bin130; also matches mask_bin13_2017 etc.
    sel += [n for n in names if pat.search(n[len("mask_"):])]
sys.stderr.write("mask vars in workspace: %d, matched: %d\n" % (len(names), len(sel)))
sys.stderr.write("matched: %s\n" % sel)
if not sel:
    sys.stderr.write("ERROR: no mask var matched. Available: %s\n" % names)
    sys.exit(1)
print(",".join(n + "=1" for n in sel))
PYEOF
fi

while IFS= read -r line; do
    rm -f filelist
    echo "$line" > filelist
    gluino=$(awk -F '[^0-9]*' '{print $5}' filelist)
    chargino=$(awk -F '[^0-9]*' '{print $6}' filelist)
    echo "gluino ${gluino}"
    echo "chargino ${chargino}"
    xrdcp -f root://cmseos.fnal.gov//store/user/tmishra/CombinedLimit/${model}/cards/cards_Run2_${mode}/counting_${model_lower}_${gluino}_${chargino}.txt .
    
    # Unblinded combine run
    if [[ "$maskBins" == "true" ]]; then
        text2workspace.py counting_${model_lower}_${gluino}_${chargino}.txt --channel-masks -o ws_${gluino}_${chargino}.root
        maskParams=$(python3 getMasks.py ws_${gluino}_${chargino}.root ${MASK_BINS[@]}) || { echo "ERROR: masking failed for ${gluino}_${chargino}, skipping"; rm -f ws_${gluino}_${chargino}.root counting_${model_lower}_${gluino}_${chargino}.txt; continue; }
        echo "setParameters: ${maskParams}"
        combine -M AsymptoticLimits ws_${gluino}_${chargino}.root -n ${model_lower} --setParameters ${maskParams} | tee limit_${gluino}_${chargino}.txt
        rm -f ws_${gluino}_${chargino}.root
        if grep -q ">>> 0 out of" limit_${gluino}_${chargino}.txt; then
            echo "ERROR: no channels masked for ${gluino}_${chargino}, not saving"; rm -f counting_${model_lower}_${gluino}_${chargino}.txt; continue
        fi
    else
        combine -M AsymptoticLimits counting_${model_lower}_${gluino}_${chargino}.txt -n ${model_lower} | tee limit_${gluino}_${chargino}.txt
    fi
    # Full destination file name: a missing directory is created instead of being written as a file
    if ! xrdcp -f limit_${gluino}_${chargino}.txt root://cmseos.fnal.gov/${limitPath}/limit_${gluino}_${chargino}.txt; then
        echo "ERROR: could not write ${limitPath}/limit_${gluino}_${chargino}.txt"
        nFailed=$((nFailed+1))
    fi
    rm counting_${model_lower}_${gluino}_${chargino}.txt
done < "$name_list"

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_14_1_0_pre4
if [[ $nFailed -gt 0 ]]; then
    echo "ERROR: ${nFailed} limit file(s) not written to EOS"
    exit 1
fi
echo "End.. Bye now.."
