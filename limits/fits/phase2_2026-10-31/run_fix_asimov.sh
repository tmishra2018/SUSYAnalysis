#!/bin/bash
# Rerun of the Asimov impacts that failed (no cards_<era>_bonly / cards_Run2_splusb): unblind cards, observation unused.
while systemctl --user is-active --quiet phase2-impacts || systemctl --user is-active --quiet phase2-gof; do sleep 60; done
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval `scramv1 runtime -sh` 2>/dev/null; cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits
for f in fix_t5wg_bonly fix_t5wg_Run2_splusb fix_t6wg_Run2; do echo "[$(date '+%F %T')] $f"; bash /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/$f.sh > /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/log_$f.txt 2>&1; grep -E "FAILED|failed" /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/fits/phase2_2026-10-31/log_$f.txt | head -3; done
echo "[$(date '+%F %T')] Asimov impact fixes done"
