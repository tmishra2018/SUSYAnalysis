#!/bin/bash
# Finishes the genWeight rerun for ZGToLLG only: waits for the ZGToLLG jobs started by
# run_ISRMC_genWeight.sh, runs any that did not finish ("[DONE]" missing from the log),
# merges, then runs steps 3-4.
# The rare MC (TTJets ~1.2 TB, TTGJets, WWG, WZG) only enters the step-3 plots, not the
# weights (data / Z+gamma, plot_ISRweight.C), so it keeps the previous trees for now.
#   nohup bash run_ISRMC_genWeight_ZGonly.sh >& /eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/ISRMC_genWeight_ZGonly.log &
set -e
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
cd SUSYAnalysis/test/ISRReweight
unset PYTHONHOME PYTHONPATH
ISRDIR=/eos/uscms/store/user/tmishra/ISRweighting
LOG=${ISRDIR}/logs_ISR/ISRMC_genWeight

while pgrep -f "analysis_ISRMC.exe .* ZGToLLG" > /dev/null; do sleep 60; done

for yv in "2018 0" "2017 0" "2016 1" "2016 0"; do
  set -- ${yv}
  for ch in mg egamma; do
    l=${LOG}/ZGToLLG_${ch}_$1_$2.log
    if ! grep -q "^\[DONE\]" ${l} 2>/dev/null; then
      echo "running ZGToLLG ${ch} $1 $2"
      ./analysis_ISRMC.exe $1 $2 ZGToLLG ${ch} > ${l} 2>&1 &
    fi
  done
done
wait
for l in ${LOG}/ZGToLLG_*.log; do
  grep -q "^\[DONE\]" ${l} || { echo "not finished: ${l}"; exit 1; }
done

for ch in mg egamma; do
  chdir=$([[ "${ch}" == "mg" ]] && echo "mgamma" || echo "egamma")
  for era in 2016preVFP 2016postVFP 2017 2018; do
    hadd -f ${ISRDIR}/${chdir}/resTree_ISR_ZGToLLG_${era}.root ${ISRDIR}/resTree_ISR_ZGToLLG_${ch}_${era}.root > /dev/null
  done
done
echo "ZGToLLG merged; running steps 3-4"
bash run_all_ISR.sh 3 4
echo "ALL DONE"
