#!/bin/bash
# Condor job: observed AsymptoticLimits for one TChiWG Run 2 card (one mass point per job).
# Args: <mass> [mode=unblind]
# Writes limit_observed_<mass>.txt to /store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2/
# (the file name and directory read by tools/TChiWG/plot_limit.C).
# First try: default combine settings (1 h timeout). If no "Observed Limit" line comes out (point 625 hung
# with the defaults on 2026-10-29), retry with --rMin 0 --rMax 5 --cminDefaultMinimizerStrategy 0.
# Submitted by runLimitTChiWG.py.

mass=$1
mode=${2:-unblind}
card=counting_tchiwg_${mass}.txt
cardPath=/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_Run2_${mode}/${card}
limitPath=/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2
out=limit_observed_${mass}.txt

cd ${_CONDOR_SCRATCH_DIR}
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 project CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src/
eval `scramv1 runtime -sh`
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
scramv1 b clean
scramv1 b -j 4
eval `scramv1 runtime -sh`
cd HiggsAnalysis/CombinedLimit

echo "Mass: ${mass}  mode: ${mode}"
if ! xrdcp -f root://cmseos.fnal.gov/${cardPath} .; then
    echo "ERROR: cannot copy ${cardPath}"; exit 1
fi

timeout 3600 combine -M AsymptoticLimits ${card} -n tchiwgobserved | tee ${out}
if ! grep -q "Observed Limit" ${out}; then
    echo "WARNING: default combine gave no observed limit for ${mass}; retrying with --rMin 0 --rMax 5"
    timeout 7200 combine -M AsymptoticLimits ${card} -n tchiwgobserved --rMin 0 --rMax 5 \
        --cminDefaultMinimizerStrategy 0 | tee ${out}
fi
if ! grep -q "Observed Limit" ${out}; then
    echo "ERROR: no observed limit for ${mass}"; exit 1
fi

xrdfs root://cmseos.fnal.gov mkdir -p ${limitPath}
if ! xrdcp -f ${out} root://cmseos.fnal.gov/${limitPath}/${out}; then
    echo "ERROR: could not write ${limitPath}/${out}"; exit 1
fi
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_14_1_0_pre4
echo "End.. Bye now.."
