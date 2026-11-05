#!/bin/bash
# Saturated-model GOF, T5WG 1700/1100, Run 2 + 4 eras, 500 toys; runs after phase2-impacts (one node, 8 cores).
while systemctl --user is-active --quiet phase2-impacts; do sleep 60; done
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval `scramv1 runtime -sh` 2>/dev/null; cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits
echo "[$(date '+%F %T')] GOF 500 toys start"; bash /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/gof_t5wg_1700_1100_500toys.sh > /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/log_gof500.txt 2>&1; echo "[$(date '+%F %T')] GOF done"
grep -h "p-value\|p =" /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/log_gof500.txt | tail -10
