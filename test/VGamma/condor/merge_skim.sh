#!/bin/bash
# Merge the condor skim chunks (submit_skim.py) into one tree per sample and era:
#   /eos/uscms/store/user/tmishra/VGamma/<tag>/resTree_VGamma_<sample>_<era>.root
# Checks that every range job of the DAG left its chunk; stops otherwise.
#   bash merge_skim.sh [tag]      (default tag zgisr; run in the CMSSW environment)
set -e
TAG=${1:-zgisr}
DAG=/uscms_data/d3/tmishra/vg_condor/${TAG}/skim.dag
OUT=/eos/uscms/store/user/tmishra/VGamma/${TAG}
# chunks are listed and read over xrootd: the /eos FUSE listing can lag minutes behind newly written files
EOS=root://cmseos.fnal.gov
XOUT=/store/user/tmishra/VGamma/${TAG}
for key in $(grep "^JOB range_" $DAG | awk '{print $2}' | sed 's/^range_//; s/_[0-9]*$//' | sort -u); do
  expected=$(grep -c "^JOB range_${key}_[0-9]* " $DAG)
  chunks=$(xrdfs $EOS ls ${XOUT}/${key}/chunks | grep "/resTree_VGamma_${key}_[0-9]*\.root$" | sed "s|^|$EOS/|")
  found=$(echo "$chunks" | grep -c root)
  if [ "$found" != "$expected" ]; then echo "${key}: ${found}/${expected} chunks, not merged"; missing=1; continue; fi
  hadd -f ${OUT}/resTree_VGamma_${key}.root $chunks > /dev/null
  echo "${key}: ${found} chunks merged"
done
[ -z "$missing" ] || { echo "some samples are incomplete"; exit 1; }
