#!/bin/bash
# Rerun the eγ, μγ and MET-PD trigger macros on data (golden JSON + new 2016 files), NJOBS at a time,
# then hadd each into its period file. DY MC outputs are unchanged and not rerun.
# usage: setsid nohup bash rerun_trigger_data.sh > /eos/uscms/store/user/tmishra/Trigger/logs_rerun/runner.log 2>&1 &
NJOBS=${NJOBS:-3}
HERE=$(cd $(dirname $0) && pwd)
F=/eos/uscms/store/user/tmishra/Trigger/files
L=/eos/uscms/store/user/tmishra/Trigger/logs_rerun
mkdir -p $L

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $HERE/../../.. && eval `scramv1 runtime -sh` 2>/dev/null
cd $HERE
echo "[$(date '+%F %T')] compiling on $(hostname)"
# analysis_egTrigger.C picks data/DY with a hardcoded "bool useData"; build the data exe from a copy with it set to true
sed 's/^bool useData = false;/bool useData = true;/' analysis_egTrigger.C > analysis_egTrigger_Data.C
grep -q '^bool useData = true;' analysis_egTrigger_Data.C || { echo "useData switch not found"; exit 1; }
g++ `root-config --cflags` analysis_egTrigger_Data.C -o analysis_egTrigger_Data.exe `root-config --libs` || exit 1
rm -f analysis_egTrigger_Data.C
g++ `root-config --cflags` analysis_mgTrigger.C     -o analysis_mgTrigger.exe     `root-config --libs` || exit 1
g++ `root-config --cflags` analysis_mgTrigger_MET.C -o analysis_mgTrigger_MET.exe `root-config --libs` || exit 1

# exe year arg output  (the macros return 0 even on failure, so success = fresh output file and no "Error:" in the log)
JOBS="
eg 2016 B-ver1 plot_egTrigger_Data_2016B-ver1
eg 2016 B-ver2 plot_egTrigger_Data_2016B-ver2
eg 2016 C plot_egTrigger_Data_2016C
eg 2016 D plot_egTrigger_Data_2016D
eg 2016 E plot_egTrigger_Data_2016E
eg 2016 F_preVFP plot_egTrigger_Data_2016F_preVFP
eg 2016 F_postVFP plot_egTrigger_Data_2016F_postVFP
eg 2016 G plot_egTrigger_Data_2016G
eg 2016 H plot_egTrigger_Data_2016H
eg 2017 B plot_egTrigger_Data_2017B
eg 2017 C plot_egTrigger_Data_2017C
eg 2017 D plot_egTrigger_Data_2017D
eg 2017 E plot_egTrigger_Data_2017E
eg 2017 F plot_egTrigger_Data_2017F
eg 2018 A plot_egTrigger_Data_2018A
eg 2018 B plot_egTrigger_Data_2018B
eg 2018 C plot_egTrigger_Data_2018C
eg 2018 D plot_egTrigger_Data_2018D
mg 2016 B-ver1 plot_MuonTrigger_Data_2016B-ver1
mg 2016 B-ver2 plot_MuonTrigger_Data_2016B-ver2
mg 2016 C plot_MuonTrigger_Data_2016C
mg 2016 D plot_MuonTrigger_Data_2016D
mg 2016 E plot_MuonTrigger_Data_2016E
mg 2016 F_preVFP plot_MuonTrigger_Data_2016F_preVFP
mg 2016 F_postVFP plot_MuonTrigger_Data_2016F_postVFP
mg 2016 G plot_MuonTrigger_Data_2016G
mg 2016 H plot_MuonTrigger_Data_2016H
mg 2017 D plot_MuonTrigger_Data_2017D
mg 2017 E plot_MuonTrigger_Data_2017E
mg 2017 F plot_MuonTrigger_Data_2017F
mg 2018 A plot_MuonTrigger_Data_2018A
mg 2018 B plot_MuonTrigger_Data_2018B
mg 2018 C plot_MuonTrigger_Data_2018C
mg 2018 D plot_MuonTrigger_Data_2018D
met 2016 1 plot_MuonTrigger_MET_2016preVFP
met 2016 0 plot_MuonTrigger_MET_2016postVFP
met 2017 1 plot_MuonTrigger_MET_2017
met 2018 1 plot_MuonTrigger_MET_2018
"
touch $L/.start
echo "$JOBS" | grep -v '^$' | xargs -P $NJOBS -L 1 bash -c '
  case $0 in eg) exe=analysis_egTrigger_Data.exe;; mg) exe=analysis_mgTrigger.exe;; met) exe=analysis_mgTrigger_MET.exe;; esac
  echo "[$(date "+%F %T")] start $3"
  ./$exe $1 $2 > '$L'/$3.log 2>&1
  if [ $? -eq 0 ] && [ '$F'/$3.root -nt '$L'/.start ] && ! grep -q "Error:" '$L'/$3.log; then
    echo "[$(date "+%F %T")] done  $3"
  else echo "[$(date "+%F %T")] FAILED $3 (see '$L'/$3.log)"; touch '$L'/$3.FAILED; fi'

if ls $L/*.FAILED >/dev/null 2>&1; then echo "[$(date '+%F %T')] some jobs failed, not merging"; exit 1; fi

echo "[$(date '+%F %T')] merging"
m() { out=$1; shift; hadd -f $F/$out.root "${@/#/$F/}" > $L/hadd_$out.log 2>&1 && echo "merged $out" || echo "FAILED hadd $out"; }
m plot_egTrigger_Data_2016preVFP  plot_egTrigger_Data_2016{B-ver1,B-ver2,C,D,E,F_preVFP}.root
m plot_egTrigger_Data_2016postVFP plot_egTrigger_Data_2016{F_postVFP,G,H}.root
m plot_egTrigger_Data_2017        plot_egTrigger_Data_2017{B,C,D,E,F}.root
m plot_egTrigger_Data_2018        plot_egTrigger_Data_2018{A,B,C,D}.root
m plot_MuonTrigger_Data_2016preVFP  plot_MuonTrigger_Data_2016{B-ver1,B-ver2,C,D,E,F_preVFP}.root
m plot_MuonTrigger_Data_2016postVFP plot_MuonTrigger_Data_2016{F_postVFP,G,H}.root
m plot_MuonTrigger_Data_2017        plot_MuonTrigger_Data_2017{D,E,F}.root
m plot_MuonTrigger_Data_2018        plot_MuonTrigger_Data_2018{A,B,C,D}.root
echo "[$(date '+%F %T')] all done"
