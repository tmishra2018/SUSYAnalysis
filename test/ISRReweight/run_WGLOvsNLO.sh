#!/bin/bash
# LO vs NLO W+gamma ISR pT comparison (reviewer question on ISR weights).
# Splits each sample into 5M-event chunks, runs compare_WG_LOvsNLO.exe in parallel,
# then hadds per era and sample into ${OUT}/WGLOvsNLO_<sample>_<era>.root.
#   bash run_WGLOvsNLO.sh   (MAX_JOBS=N to change parallelism; SAMPLES="WGToLNuG WGJet40 WGJet130" reruns only LO)
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
cd SUSYAnalysis/test/ISRReweight

MAX_JOBS="${MAX_JOBS:-7}"
SAMPLES="${SAMPLES:-WGToLNuG WGJet40 WGJet130 NLO}"
CHUNK=5000000
OUT=/eos/uscms/store/user/tmishra/ISRweighting/WGLOvsNLO
mkdir -p ${OUT}/jobs ${OUT}/logs ${OUT}/lists

g++ $(root-config --cflags) ../../lib/libAnaClasses.so compare_WG_LOvsNLO.C -o compare_WG_LOvsNLO.exe $(root-config --libs)

declare -A NLOTAG=( [2016preVFP]=UL16MiniAODAPVv2 [2016postVFP]=UL16MiniAODv2 [2017]=UL17 [2018]=UL18 )
JOBLIST=${OUT}/jobs/joblist.txt
: > ${JOBLIST}
for era in 2016preVFP 2016postVFP 2017 2018; do
  year=${era:0:4}; vfp=0; [ "${era}" == "2016preVFP" ] && vfp=1
  list=${OUT}/lists/nlo_${era}.txt
  find /eos/uscms/store/user/tmishra/WGToLNuG_01J*RunIISummer20${NLOTAG[$era]}* -name "*.root" | sort > ${list}
  for sample in ${SAMPLES}; do
    if [ "${sample}" == "NLO" ]; then
      n=$(root -l -b -q -e "TChain c(\"ggNtuplizer/EventTree\"); std::ifstream in(\"${list}\"); std::string l; while(std::getline(in,l)) c.Add(l.c_str()); std::cout<<\"NEVT \"<<c.GetEntries()<<std::endl;" 2>/dev/null | awk '/NEVT/{print $2}')
    else
      n=$(root -l -b -q -e "TChain c(\"ggNtuplizer/EventTree\"); c.Add(\"/eos/uscms/store/user/tmishra/InputFilesMC/${sample}/${sample}_${era}.root\"); std::cout<<\"NEVT \"<<c.GetEntries()<<std::endl;" 2>/dev/null | awk '/NEVT/{print $2}')
    fi
    echo "${era} ${sample}: ${n} events"
    i=0
    for ((first=0; first<n; first+=CHUNK)); do
      last=$((first+CHUNK))
      echo "${year} ${vfp} ${sample} ${first} ${last} ${OUT}/jobs/${sample}_${era}_${i}.root ${list} ${OUT}/logs/${sample}_${era}_${i}.log" >> ${JOBLIST}
      i=$((i+1))
    done
  done
done
echo "$(wc -l < ${JOBLIST}) jobs, ${MAX_JOBS} in parallel"

cat ${JOBLIST} | xargs -P ${MAX_JOBS} -L 1 bash -c './compare_WG_LOvsNLO.exe $0 $1 $2 $3 $4 $5 $6 > $7 2>&1 || echo "FAILED $5"'

for era in 2016preVFP 2016postVFP 2017 2018; do
  for sample in ${SAMPLES}; do
    hadd -f ${OUT}/WGLOvsNLO_${sample}_${era}.root ${OUT}/jobs/${sample}_${era}_*.root > /dev/null
  done
done
echo "done: ${OUT}/WGLOvsNLO_<sample>_<era>.root"
