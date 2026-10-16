#!/bin/bash
# Convert the ggNtuples of one year into TnP trees with makeTnPTree.exe (4 jobs in parallel).
# usage: ./convert_year.sh <tag>    tag = 2016preVFP | 2016postVFP | 2017 | 2018
TAG=$1
HERE=$(cd $(dirname $0) && pwd)
E=/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni
O=/eos/uscms/store/user/tmishra/TnP_ggNtuple/$TAG
mkdir -p $O
case $TAG in
  2016preVFP)  YEAR=2016; DATA="SingleElectron/SingleElectron_2016B-ver1 SingleElectron/SingleElectron_2016B-ver2 SingleElectron/SingleElectron_2016C SingleElectron/SingleElectron_2016D SingleElectron/SingleElectron_2016E SingleElectron/SingleElectron_2016F_preVFP"; DY=DYJetsToLL_2016preVFP ;;
  2016postVFP) YEAR=2016; DATA="SingleElectron/SingleElectron_2016F_postVFP SingleElectron/SingleElectron_2016G SingleElectron/SingleElectron_2016H"; DY=DYJetsToLL_2016postVFP ;;
  2017)        YEAR=2017; DATA="SingleElectron/SingleElectron_2017B SingleElectron/SingleElectron_2017C SingleElectron/SingleElectron_2017D SingleElectron/SingleElectron_2017E SingleElectron/SingleElectron_2017F"; DY=DYJetsToLL_2017 ;;
  # 2018: EGamma PD (files named DoubleEG_2018*), single-electron trigger HLT_Ele32_WPTight_Gsf
  2018)        YEAR=2018; DATA="DoubleEG/DoubleEG_2018A DoubleEG/DoubleEG_2018B DoubleEG/DoubleEG_2018C DoubleEG/DoubleEG_2018D"; DY=DYJetsToLL_2018 ;;
  *) echo "unknown tag $TAG"; exit 1 ;;
esac

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $HERE/../../.. && eval `scramv1 runtime -sh` 2>/dev/null
cd $HERE/..
{ for d in $DATA; do echo "$E/$d.root $O/tnp_$(basename $d).root 0"; done
  echo "$E/DYJetsToLL/$DY.root $O/tnp_DY_${DY#DYJetsToLL_}.root 1"; } |
  xargs -P 4 -L 1 sh -c "./makeTnPTree.exe \$0 \$1 $YEAR \$2 > \$1.log 2>&1; echo \$(basename \$1): \$(tail -1 \$1.log)"
