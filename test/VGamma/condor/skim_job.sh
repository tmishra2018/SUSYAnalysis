#!/bin/bash
# Condor worker for the V+gamma skim (analysis_VGamma_Run2_condor.exe), see submit_skim.py.
#   skim_job.sh norm  <year> <preVFP> <sample> <first> <last> <eosdir>   -> <eosdir>/norm/part_<first>.root
#   skim_job.sh merge <eosdir>                                            -> <eosdir>/norm.root (hadd of the parts)
#   skim_job.sh range <year> <preVFP> <sample> <first> <last> <eosdir>   -> <eosdir>/chunks/resTree_*_<first>.{root,log}
# <eosdir> is an EOS path without the /eos/uscms prefix (e.g. /store/user/tmishra/VGamma/zgisr/ZGToLLG_2018).
set -e
# outputs are copied to EOS with xrdcp and then deleted: condor would otherwise transfer every new file in the
# sandbox back to the submit directory (that filled the 3 GB home quota on 2026-10-29)
EOS=root://cmseos.fnal.gov
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scramv1 project CMSSW CMSSW_14_1_0_pre4 > /dev/null
eval `cd CMSSW_14_1_0_pre4/src && scramv1 runtime -sh`
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH VG_XROOTD=1 VG_PUDIR=$PWD/PUfiles
mode=$1
if [ "$mode" == "merge" ]; then
  d=$2
  parts=$(xrdfs $EOS ls $d/norm | grep "part_.*\.root$" | sed "s|^|$EOS/|")
  hadd -f norm.root $parts
  xrdcp -f norm.root $EOS/$d/norm.root
  rm -f norm.root
  exit 0
fi
year=$2; vfp=$3; sample=$4; first=$5; last=$6; d=$7
if [ "$mode" == "norm" ]; then
  ./analysis_VGamma_Run2_condor.exe $year $vfp $sample norm $first $last part_${first}.root
  xrdcp -f part_${first}.root $EOS/$d/norm/part_${first}.root
  xrdcp -f resTree_VGamma_*_${first}_norm.log $EOS/$d/norm/ || true
  rm -f part_${first}.root resTree_VGamma_*_${first}_norm.log
else
  xrdcp -f $EOS/$d/norm.root norm.root
  ./analysis_VGamma_Run2_condor.exe $year $vfp $sample range $first $last norm.root
  for f in resTree_VGamma_*_${first}.root resTree_VGamma_*_${first}.log; do xrdcp -f $f $EOS/$d/chunks/$f && rm -f $f; done
  rm -f norm.root
fi
