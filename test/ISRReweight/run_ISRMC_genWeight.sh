#!/bin/bash
# Re-runs only the MC part of ISR step 1 (analysis_ISRMC, now with the genWeight sign)
# and its merge, then steps 3-4. Step 2 (jet->photon fits) uses data and GJet MC only,
# so its fake rates are unchanged and it is not rerun.
#   bash run_ISRMC_genWeight.sh   (MAX_JOBS=N to change parallelism)
set -e
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
cd SUSYAnalysis/test/ISRReweight
unset PYTHONHOME PYTHONPATH

MAX_JOBS="${MAX_JOBS:-7}"
ISRDIR=/eos/uscms/store/user/tmishra/ISRweighting
LOG=${ISRDIR}/logs_ISR/ISRMC_genWeight
mkdir -p ${LOG}

g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_ISRMC.C -o analysis_ISRMC.exe $(root-config --libs)

JOBS=${LOG}/joblist.txt
: > ${JOBS}
for sample in ZGToLLG TTJets TTGJets WWG WZG; do
  for ch in mg egamma; do
    for yv in "2018 0" "2017 0" "2016 1" "2016 0"; do
      echo "${yv} ${sample} ${ch}" >> ${JOBS}
    done
  done
done
cat ${JOBS} | xargs -P ${MAX_JOBS} -L 1 bash -c './analysis_ISRMC.exe $0 $1 $2 $3 > '${LOG}'/$2_$3_$0_$1.log 2>&1 && echo "[OK] $*" || echo "[FAIL] $*"'
if grep -q "\[FAIL\]" ${LOG}/../ISRMC_genWeight.log 2>/dev/null; then echo "failures, stopping"; exit 1; fi

for ch in mg egamma; do
  chdir=$([[ "${ch}" == "mg" ]] && echo "mgamma" || echo "egamma")
  for sample in ZGToLLG TTJets TTGJets WWG WZG; do
    for era in 2016preVFP 2016postVFP 2017 2018; do
      hadd -f ${ISRDIR}/${chdir}/resTree_ISR_${sample}_${era}.root ${ISRDIR}/resTree_ISR_${sample}_${ch}_${era}.root > /dev/null
    done
  done
done
echo "MC merged; running steps 3-4"
bash run_all_ISR.sh 3 4
echo "ALL DONE"
