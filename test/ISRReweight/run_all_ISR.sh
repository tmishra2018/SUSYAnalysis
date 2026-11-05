#!/usr/bin/env bash
# run_all_ISR.sh
# Runs the full ISR weight derivation in order, stopping at the first failure:
#   step 1: run_Reweight_step1.sh   skims (data, MC, hadron proxy) + hadd merge
#   step 2: runFit_step2.sh         jet->photon fake rates (~6-7 h)
#   step 3: run_ISRweight_step3.sh  derive ISR weights + weighted check plots
#   step 4: format_isrweights.py    C++ blocks for analysis_VGBkg.C / pred_VGBkg.C
#
# Usage (from any directory, any shell):
#   bash run_all_ISR.sh                 # steps 1-4
#   bash run_all_ISR.sh 2               # start at step 2
#   bash run_all_ISR.sh 2 3             # steps 2-3 only
#   MAX_JOBS=12 bash run_all_ISR.sh     # parallel jobs in step 1 (default 8)
# Long run in the background (tcsh):
#   nohup bash run_all_ISR.sh >& /eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/run_all_ISR.log &
#
# Per-step logs: /eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/step{1,2,3}.log, weight blocks: isr_blocks.txt

set -euo pipefail

FIRST="${1:-1}"
LAST="${2:-4}"

CMSSW_SRC=/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src
WORKDIR=${CMSSW_SRC}/SUSYAnalysis/test/ISRReweight
# logs go to EOS: a verbose step2.log once filled the 3 GB home quota
LOGDIR=/eos/uscms/store/user/tmishra/ISRweighting/logs_ISR

# ---- Environment (cmsenv) --------------------------------------------------
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scramv1 runtime -sh)"
set -u

cd "${WORKDIR}"
mkdir -p "${LOGDIR}"

stamp(){ date '+%F %T'; }

run_step(){
  local n="$1" label="$2"; shift 2
  if (( n < FIRST || n > LAST )); then
    echo "[$(stamp)] skip step ${n}: ${label}"
    return 0
  fi
  echo "[$(stamp)] START step ${n}: ${label}  (log: ${LOGDIR}/step${n}.log)"
  if "$@" > "${LOGDIR}/step${n}.log" 2>&1; then
    echo "[$(stamp)] DONE  step ${n}: ${label}"
  else
    echo "[$(stamp)] FAIL  step ${n}: ${label} -- see ${LOGDIR}/step${n}.log"
    tail -20 "${LOGDIR}/step${n}.log"
    exit 1
  fi
}

# ---- Guard for step 2: fake-rate fit must use the current photon ID --------
# Current loose ID (src/analysis_photon.cc, jetFakePho/fit/FitJetFake.C):
#   EB sigmaIetaIeta < 0.0106, ChIso < 1.694;  EE 0.0272, 2.089
# FitJetFake.C still carrying the 2016 values (0.0103 / 1.295) measures fakes
# for a different photon selection, so step 2 never runs with them.
if (( FIRST <= 2 && LAST >= 2 )); then
  if grep -q 'StandardCut = 0.0103' FitJetFake.C; then
    echo "[$(stamp)] STOP: FitJetFake.C still uses the 2016 photon ID cuts (StandardCut = 0.0103)."
    echo "        Update them to the current ID first."
    if (( FIRST == 1 )); then
      echo "        Running step 1 only; restart later with: bash run_all_ISR.sh 2"
      LAST=1
    else
      exit 1
    fi
  fi
fi

run_step 1 "skims + merge"                 bash run_Reweight_step1.sh
run_step 2 "jet->photon fake rates"        bash runFit_step2.sh
run_step 3 "derive + apply ISR weights"    bash run_ISRweight_step3.sh

if (( 4 >= FIRST && 4 <= LAST )); then
  echo "[$(stamp)] START step 4: format weights -> ${WORKDIR}/isr_blocks.txt"
  python3 format_isrweights.py > isr_blocks.txt
  echo "[$(stamp)] DONE  step 4"
  echo ""
  echo "Weights: /eos/uscms/store/user/tmishra/ISRweighting/ISRweights_*_{egamma,mgamma}.txt"
  echo "Next: paste isr_blocks.txt into ../Background/analysis_VGBkg.C and ../myResult/pred_VGBkg.C,"
  echo "      re-derive VGammaScaleHist (toDeriveScale=true), then pred_VGBkg.C and createDataCard.py."
fi

echo "[$(stamp)] finished steps ${FIRST}-${LAST}"
