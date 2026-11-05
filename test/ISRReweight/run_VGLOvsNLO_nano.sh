#!/bin/bash
# Generator-level LO vs NLO ISR pT for W+gamma, DY and Z+gamma from UL18 NanoAODv9 (compare_VG_LOvsNLO_nano.py).
# Z+gamma has no LO UL sample with m_ll > 15 GeV, so only its NLO signed vs unsigned shapes are compared.
#   nohup setsid bash run_VGLOvsNLO_nano.sh >> /eos/uscms/store/user/tmishra/ISRweighting/VGLOvsNLO_nano/run.log 2>&1 < /dev/null &
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
cd SUSYAnalysis/test/ISRReweight
MAX_JOBS="${MAX_JOBS:-8}"
OUT=/eos/uscms/store/user/tmishra/ISRweighting/VGLOvsNLO_nano
mkdir -p ${OUT}/jobs ${OUT}/logs
declare -A DS=(
  [WG_LO]=/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
  [WG_NLO]=/WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
  [DY_LO]=/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
  [DY_NLO]=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM
  [ZG_NLO]=/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM )
JOBS=${OUT}/jobs/joblist.txt; : > ${JOBS}
for s in "${!DS[@]}"; do i=0
  for f in $(dasgoclient -query="file dataset=${DS[$s]}"); do
    echo "root://cmsxrootd.fnal.gov/${f} ${OUT}/jobs/${s}_${i}.root ${OUT}/logs/${s}_${i}.log" >> ${JOBS}; i=$((i+1)); done
  echo "${s}: ${i} files"; done
echo "$(date '+%F %T') $(wc -l < ${JOBS}) jobs, ${MAX_JOBS} in parallel"
# files whose log already has [DONE] are skipped, so a stopped run can be resumed by running the script again
cat ${JOBS} | xargs -P ${MAX_JOBS} -L 1 bash -c 'grep -q DONE $2 2>/dev/null && exit 0; for t in 1 2 3; do python3 compare_VG_LOvsNLO_nano.py $0 $1 > $2 2>&1 && break; sleep 30; done; grep -q DONE $2 || echo "FAILED $0"'
for s in "${!DS[@]}"; do hadd -f ${OUT}/VGLOvsNLO_nano_${s}.root ${OUT}/jobs/${s}_*.root > /dev/null; done
echo "$(date '+%F %T') done: ${OUT}/VGLOvsNLO_nano_<sample>.root"
